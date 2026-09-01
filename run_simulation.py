#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_simulation.py — Part 1: Simulation & Preliminary Checks
=============================================================
21 cm × Galaxy Cross-Correlation: HERA × Euclid (Lightcone)

This script runs the 21cmFASTv4 lightcone simulation, constructs the
galaxy density field, estimates the galaxy bias, applies Kaiser RSD,
and saves all outputs to an HDF5 file for downstream analysis.

Usage
-----
    # Direct execution (activate conda env first)
    conda activate 21cmfast
    python run_simulation.py

    # SLURM job submission
    sbatch submit_job.sh

Output
------
    outputs/lightcone_data.h5  — all simulation fields and metadata,
                                  read by plot_fields.ipynb and analysis.ipynb

References
----------
    Davies, Mesinger & Murray (2025) — arXiv:2504.17254
    Gagnon-Hartman, Davies & Mesinger (2025) — arXiv:2502.20447
    La Plante et al. (2023) — arXiv:2205.09770
    Euclid Collaboration (2022) — arXiv:2108.01201
"""

import argparse
import gc
import os
import shutil
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC (no display required)
from scipy.integrate import quad, simpson
import h5py
import warnings

warnings.filterwarnings("ignore")

# Repo root on sys.path so the shared helpers in src/ are importable when this
# script is launched from anywhere.
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Shared calibrations — the single source of truth for the UV-SFR conversion,
# the AB zero point, and 21cmFAST's star-formation timescale. Using these
# rather than local copies keeps Part 1 consistent with the analysis stage.
from src.analysis import (
    GALAXY_WEIGHTING_MODES,
    T_STAR_DEFAULT,
    effective_galaxy_bias,
    galaxy_overdensity_from_catalogue,
    select_euclid_halos,
    star_formation_timescale,
    stellar_mass_to_sfr,
)
from src.conversions import (
    Muv_to_Luv,
    cell_mass,
    sfr_to_Luv,
    sheth_tormen_bias,
    survey_area_to_box_size,
)
from src.provenance import (
    INT32_MAX,
    estimate_cache_footprint,
    RunManifest,
    estimate_catalogue_cost,
    resolve_n_threads,
)

# ── Check whether 21cmFAST is installed ──────────────────────────────────────
try:
    import py21cmfast as p21c
    from astropy.cosmology import Planck18
    import astropy.units as u
    HAS_21CMFAST = True
    print(f"21cmFAST version: {p21c.__version__}")
except ImportError:
    HAS_21CMFAST = False
    print("21cmFAST not installed — will use synthetic fields for demonstration.")
    print("Install with:  pip install 21cmFAST")

# ── HMF library ──────────────────────────────────────────────────────────────
try:
    from hmf import MassFunction
    HAS_HMF = True
except ImportError:
    HAS_HMF = False
    print("hmf not installed — galaxy bias will use the default configured value.")


# ===========================================================================
#  Command line
# ===========================================================================
# The configuration block below is the single source of truth for a real run
# and takes no arguments, by design.  The one flag here does not change any of
# it: `--smoke-test` swaps in a separate, tiny configuration *after* the block,
# so the production values stay exactly as written and exactly as documented.

def _parse_args(argv=None):
    """
    Parse the simulation script's only command-line option.

    Parameters
    ----------
    argv : sequence of str, optional
        Argument list; defaults to ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        With a single ``smoke_test`` boolean.
    """
    parser = argparse.ArgumentParser(
        description=(
            "21cmFAST lightcone + halo catalogue + galaxy field. "
            "Simulation parameters are set in the configuration block, not here."
        )
    )
    parser.add_argument(
        "--smoke-test", action="store_true",
        help="run a tiny end-to-end configuration to verify every stage "
             "executes and produces correctly-shaped output. NOT a science "
             "run: see src/smoke_test.py. Writes to outputs/smoke_test/.",
    )
    return parser.parse_args(argv)


_ARGS = _parse_args()
SMOKE_TEST = bool(_ARGS.smoke_test)


# ===========================================================================
#  ★ CONFIGURATION — ALL USER-ADJUSTABLE PARAMETERS ★
#  Edit values in this section only. All subsequent code reads from here.
#
#  These are the production values. They are documented in README.md,
#  docs/HPC.md §13 and docs/simulation_spec.md, and --smoke-test does not
#  edit them — it reassigns a handful of names in the clearly-marked block
#  immediately after the end of this section.
# ===========================================================================

# ---------------------------------------------------------------------------
#  Survey footprint  ->  simulation box geometry
# ---------------------------------------------------------------------------
# The box is sized from the survey being forecast rather than chosen by hand.
# Euclid Deep Field Fornax: 10 deg^2, centred RA 03:31:43.6, Dec -28:05:18.6.
SURVEY_AREA_DEG2 = 10.0     # Euclid Deep Field Fornax footprint  [deg^2]
SURVEY_Z_CENTRAL = 7.0      # central redshift of the analysis

# sigma_z is the *absolute* photometric redshift error, NOT sigma_z/(1+z):
# radial_smearing_length() computes sigma_r = c sigma_z / H(z) directly, and
# its docstring states the convention explicitly.  Surveys quote the
# fractional form, so it must be multiplied by (1+z) before it lands here.
#
# Source: Euclid Collaboration: Allen et al. (2026), A&A 711, A25, Sect. 3 and
# Fig. 4 -- sigma_nmad <= 0.032, the measured scatter between true and
# photometric redshift in their synthetic Euclid Deep Field catalogue test.
# This replaces the previous 0.45, which came from Euclid's *pre-launch*
# fractional requirement sigma_z/(1+z) < 0.05 -- a specification, not a
# measurement.  Converting on the standard NMAD convention:
#     sigma_z(z=7) = 0.032 * (1 + 7) = 0.256
#
# CAVEAT 1 -- field mismatch, deliberately NOT resolved.  Allen et al. give
# 0.032 as an upper bound across all three Euclid Deep Fields (EDF-N/S/F),
# and explicitly call out EDF-North as having the largest (worst) scatter.
# This pipeline's geometry is EDF-*Fornax* specifically, whose true value is
# likely smaller (better) than 0.032.  That field-specific number is not
# available as printed text in the paper -- it appears only as an annotation
# inside the Fig. 4 image, which was not extractable.  We therefore adopt the
# conservative worst-field bound, not the field-specific value.
#
# CAVEAT 2 -- normalisation ambiguity, deliberately NOT resolved.  The paper's
# caption text describes sigma_nmad as the scatter in "z_true - z_phot"
# without explicitly stating the (1+z) normalisation in the visible text, even
# though that normalisation is standard for any metric named "sigma_NMAD" in
# the photo-z literature.  The * (1+z) conversion above ASSUMES the standard
# normalised definition.  If a later full reading of Fig. 4 or Sect. 3
# contradicts this, the sigma_z = 0.256 conversion must be revisited.
#
# Alternative sources, recorded but NOT adopted -- see NUMBERS_AND_SOURCES.md
# section 5.  Neither yielded an extractable sigma_NMAD from its accessible
# text, so neither can supersede Allen et al. without a manual read:
#   - Varadaraj et al. (2026), A&A 707, A239 (arXiv:2510.00945), "Euclid:
#     Discovery of bright z ~= 7 Lyman-break galaxies in UltraVISTA and Euclid
#     COSMOS".  Selection 6.5 <= z <= 7.5 -- a much tighter match to this
#     pipeline's z_obs = 7 than Allen's 6 <= z < 8 bin, and the same LBG
#     population.  But it is the COSMOS field, not an EDF, so it trades one
#     field mismatch for another.  Mean photo-z uncertainties are in its
#     Fig. 4.
#   - Weaver et al. (2025), A&A 697, A16 (arXiv:2405.13505), "Euclid: Early
#     Release Observations -- NISP-only sources and the search for luminous
#     z = 6-8 galaxies".  Redshift range matches Allen's bin, but the fields
#     are the 1.5 deg^2 ERO 'Magnifying Lens' Abell cluster fields, and the
#     selection is NISP-only sources with no appreciable VIS flux -- a
#     specific, atypical population rather than a general EDF sample.
#
# Defined here, not in the Euclid block below, because it also sets the
# line-of-sight depth of the box via SURVEY_DELTA_Z.
photoz_uncertainty = 0.256  # sigma_z photometric redshift error (absolute)

# CHOICE, not a default: how many sigma_z of photo-z scatter the box spans.
# PHOTOZ_N_SIGMA = 1 -> delta_z = 0.512, z = 6.744-7.256, L_los = 179.3 Mpc
# PHOTOZ_N_SIGMA = 2 -> delta_z = 1.024, z = 6.488-7.512, L_los = 359.3 Mpc
#
# NOTE: because SURVEY_DELTA_Z is proportional to sigma_z and sigma_r is too,
# L_los / sigma_r == 2 * PHOTOZ_N_SIGMA *exactly*, whatever sigma_z is.  The
# photo-z kernel at the box's own fundamental mode is therefore pinned at
# W = exp(-pi^2 / 2) = 0.0072 and does NOT improve when sigma_z improves --
# only this choice moves it.  See NUMBERS_AND_SOURCES.md section 5.
# The forecast adopts +/-1 sigma; raising this widens the LOS extent but
# pushes further from the quasi-coeval regime the estimator assumes.
PHOTOZ_N_SIGMA = 1
SURVEY_DELTA_Z = 2 * PHOTOZ_N_SIGMA * photoz_uncertainty

# ---------------------------------------------------------------------------
#  Simulation grid  (derived — do not hardcode BOX_LEN/HII_DIM/DIM)
# ---------------------------------------------------------------------------
# HII_DIM follows from target_cell_size_mpc = 2.0, the resolution of the old
# 256 Mpc / 128^3 grid, so covering the footprint does not coarsen M_cell.
# 486.33 / 2.0 = 244 cells, snapped up to the next power of two for the FFTs.
SIM_BOX = survey_area_to_box_size(
    area_deg2=SURVEY_AREA_DEG2,
    z_central=SURVEY_Z_CENTRAL,
    delta_z=SURVEY_DELTA_Z,
    cosmo=None,              # -> astropy Planck18, as used for the endpoints
    target_cell_size_mpc=2.0,
)

HII_DIM = SIM_BOX.hii_dim    # 256 cells per side (low-res grid)
BOX_LEN = SIM_BOX.box_len    # 486.33 Mpc comoving, from the 10 deg^2 footprint
DIM     = SIM_BOX.dim        # 768, high-res grid for initial conditions

# ---------------------------------------------------------------------------
#  Halo sampler floor  (constrained by a 32-bit index, not by physics)
# ---------------------------------------------------------------------------
# The stochastic halo sampler populates every cell down to this mass, so the
# catalogue size scales with BOX_LEN^3 and NOT with the lightcone redshift
# range -- narrowing z_min/z_max does nothing to it.
#
# 21cmFAST's C backend indexes halo arrays with `int`.  `halo_coords` holds
# 3 * N_halos elements, so the catalogue overflows a signed 32-bit index once
# N_halos > INT_MAX/3 = 7.158e8, and the run SEGFAULTS regardless of how much
# memory the node has.  At the footprint-derived BOX_LEN = 486.33 Mpc the
# template default of 1e8 draws 9.370e8 halos -- 1.31x over the limit.
#
# 2e8 was chosen over shrinking BOX_LEN because it costs objects, not physics:
#
#   SAMPLER_MIN_MASS   halo count   star formation   int32 headroom
#         1e8            100 %          100 %            1.31  <-- segfaults
#       1.5e8             63.1 %         99.95 %          0.83
#       2e8               46.2 %         99.84 %          0.61  <-- adopted
#       3e8               28.7 %         99.44 %          0.38
#
# Star formation is retained almost entirely because M_TURN = 5e8 (section 4)
# puts exp(-M_TURN/M_h) on the stellar fraction: a 1e8 halo forms stars at
# exp(-5) = 0.7 % of its unsuppressed rate.  Raising the floor to 2e8 discards
# just over half the *objects* and 0.16 % of the *star formation*, so the
# ionizing budget, the 21 cm field and the Euclid selection (which lives at
# ~1e10-1e11 M_sun) are all effectively untouched.
#
# Shrinking BOX_LEN instead would need <= 444.6 Mpc, which breaks the
# traceability survey_area_to_box_size() exists to provide -- 486.33 Mpc
# encodes the 10 deg^2 EDF-Fornax footprint.
#
# Retained fractions are Sheth-Tormen cumulative counts from
# hmf.MassFunction(z=7, dlog10m=0.02); see NUMBERS_AND_SOURCES.md section 12
# and provenance.SAMPLER_RETAINED_FRACTION.
SAMPLER_MIN_MASS = 2e8       # halo sampler floor [M_sun]

# Multiple of the estimated cache size that must be free before the run starts.
# The estimate is an empirical extrapolation, and the run also writes its own
# HDF5 output alongside the cache, so it is not run right to the edge.
DISK_SAFETY_FACTOR = 1.25

# ---------------------------------------------------------------------------
#  Lightcone redshift range
# ---------------------------------------------------------------------------
# SMOKE-TEST SLAB — deliberately narrow. Delta z = 0.01 gives L_LOS = 3.5 Mpc
# at z ~ 7, which minimum_los_slices floors to N_z = 100 (a 0.035 Mpc LOS cell
# against a 2 Mpc transverse cell, ~57x oversampled). It is quasi-coeval, with
# negligible redshift evolution along the LOS.
#
# This is intentional and must stay this way for now: the power-spectrum
# estimator in src/analysis.py assumes statistical homogeneity along the LOS,
# which only holds for a quasi-coeval box. Widening to a true lightcone
# (z = 6.5 - 7.5, L_LOS = 350.8 Mpc, N_z = 175) requires the estimator work in
# TODO.md P0 first — otherwise the FFT is applied to unevenly sampled,
# redshift-evolving data and the resulting spectra are not trustworthy.
#
# DO NOT widen this range without doing TODO.md P0.1 and P0.2 first.
#
# The survey footprint implies a much wider range -- SIM_BOX.z_min/z_max =
# 6.744/7.256 (delta_z = 0.512, L_LOS = 179.3 Mpc) from the photo-z depth above.
# That is the range this forecast *should* run once TODO.md P0 lands; the
# transverse box size is already sized for it. Until then the slab below
# deliberately overrides it, so only BOX_LEN/HII_DIM/DIM are footprint-driven.
SURVEY_Z_MIN = SIM_BOX.z_min   # 6.744 — survey-derived, not used yet
SURVEY_Z_MAX = SIM_BOX.z_max   # 7.256 — survey-derived, not used yet

z_center = 7.0
z_width  = 1

z_min = z_center - z_width / 2   # 6.9
z_max = z_center + z_width / 2   # 7.1


# ---------------------------------------------------------------------------
#  Euclid-like survey parameters
# ---------------------------------------------------------------------------
M_UV_limit          = -22    # UV absolute magnitude cut
# Absolute UV magnitude window (more negative = brighter). Defined here
# because both the galaxy-field construction (section 3b) and the bias
# estimate (section 4) select on it.
M_UV_bright         = -26
M_UV_faint          = M_UV_limit
# photoz_uncertainty (sigma_z = 0.256, absolute) is set in the survey-footprint
# block above, because it now also sets the line-of-sight depth of the box.
# n_bar for the galaxy shot noise P_N,gal = 1/n_bar.
#
# Units: plain Mpc^-3, NOT h^3 Mpc^-3.  Every consumer -- the shot-noise term
# in src.analysis.compute_uncertainty_budget and the expected_counts draw in
# section 3b below -- uses it against volumes already in Mpc^3, and no factor
# of h^3 is applied anywhere.  The previous value carried an "[h^3 Mpc^-3]"
# label that nothing in the code acted on.
#
# Source: Euclid Collaboration: Allen et al. (2026), A&A 711, A25, "Euclid
# Quick Data Release (Q1) XXV. Hunting for luminous z>6 galaxies in the Euclid
# Deep Fields".  Table 2 (p. 9), row 6 <= z < 8, column Selected(DPL):
# N = 70,445 +/- 265 galaxies over the full 53 deg^2 survey to 26.5 AB depth.
# Divided by the comoving volume of that 53 deg^2, 6 < z < 8 shell,
# V = 3.2285e8 (h^-1 Mpc)^3:
#     n_bar = N / V = 2.18e-4 h^3 Mpc^-3 = 7.48e-5 Mpc^-3   (h = 0.7)
#
# CAVEAT -- cosmology mismatch, deliberately NOT reconciled.  That volume was
# computed with the *paper's own* cosmology (Allen et al., end of Sect. 1:
# Omega_m = 0.27, Omega_Lambda = 0.73, H0 = 70 km/s/Mpc), which is NOT the
# Planck18 cosmology this pipeline adopts elsewhere (Omega_m = 0.3111,
# H0 = 67.66).  Recomputing V under Planck18 would change n_bar, but that is a
# decision for the paper's authors, not an automatic fix -- do not silently
# reconcile the two.  See NUMBERS_AND_SOURCES.md section 2.
mean_galaxy_density = 7.48e-5   # n_bar  [Mpc^-3]  (see the block above)

# ---------------------------------------------------------------------------
#  How the galaxy overdensity field is weighted
# ---------------------------------------------------------------------------
# "lightcone_sfr" (default) — delta_gal from the lightcone `halo_sfr` field,
#     i.e. the per-cell SFR density integrated over all halos.  This is the
#     original behaviour and the only mode that evolves along the LOS.
#
# "number" — delta_gal = N / <N> - 1 from the Euclid-selected halo catalogue.
#     One unit of weight per detectable galaxy.
#
# "luminosity" — delta_gal,L = sum(L_UV) / <sum(L_UV)> - 1 from the same
#     selected catalogue, weighting each halo by its own UV luminosity
#     (L_UV = SFR / kappa_UV, Fisher et al. 2026, arXiv:2511.10741 Eq. 12;
#     the constant cancels out of an overdensity, so this mode is
#     insensitive to its value).
#
# The two catalogue modes are interchangeable: identical grid, identical
# normalisation, identical downstream handling.
#
# "number" is the DEFAULT and the right choice for a Euclid cross-correlation
# forecast: Euclid counts magnitude-limited galaxies, so the field whose power
# spectrum we measure has to be built from the same population the shot noise
# 1/n_bar describes.  Under "lightcone_sfr" it was not -- that field weights
# EVERY halo by its SFR, so P_gal described one tracer while P_N,gal described
# another, and the measured P_gal fell below its own shot noise for
# k > ~0.12 Mpc^-1, which no single population can do.  See docs/HPC.md
# section 11.15.
#
# CAVEAT, unchanged by this: the catalogue modes are built from the *coeval*
# perturbed catalogue at z_obs, so they carry no redshift evolution along the
# line of sight, while the 21 cm lightcone does.  The deposit now spans L_los
# so the two fields share a geometry exactly, but the galaxy field is still a
# z_obs snapshot.  A fully consistent treatment needs per-node halo catalogues
# selected and interpolated onto the lightcone; 21cmFAST's lightcone halobox
# cannot supply it, because every quantity it carries (count, halo_sfr,
# n_ion, ...) is already summed over all halos in a cell and so cannot have a
# per-halo magnitude cut applied to it.
GALAXY_WEIGHTING = "number"   # lightcone_sfr | number | luminosity

if GALAXY_WEIGHTING not in ("lightcone_sfr",) + GALAXY_WEIGHTING_MODES:
    raise ValueError(
        f"GALAXY_WEIGHTING must be 'lightcone_sfr', 'number' or 'luminosity', "
        f"got {GALAXY_WEIGHTING!r}"
    )

# ---------------------------------------------------------------------------
#  Galaxy bias (default; overwritten below if HMF estimation succeeds)
# ---------------------------------------------------------------------------
galaxy_bias = 8              # linear bias b (typical for high-z LBGs)

# ---------------------------------------------------------------------------
#  Cosmological parameters (Planck 2018)
# ---------------------------------------------------------------------------
OMEGA_M_0       = 0.315      # matter density parameter
HUBBLE_CONSTANT = 67.36      # H_0  [km s^-1 Mpc^-1]

# ---------------------------------------------------------------------------
#  Physical constants
# ---------------------------------------------------------------------------
SPEED_OF_LIGHT_KMS = 3e5         # c  [km s^-1]
SPEED_OF_LIGHT_MPS = 3e8         # c  [m s^-1]
F_21_MHZ           = 1420.405    # 21cm rest frequency  [MHz]
F_21_HZ            = F_21_MHZ * 1e6

# ---------------------------------------------------------------------------
#  HERA instrument
# ---------------------------------------------------------------------------
HERA_DISH_DIAMETER = 14.0        # dish diameter  [m]
integration_time   = 1000 * 3600 # total integration  [s]  (1000 h)
bandwidth          = 8e6         # per-band bandwidth  [Hz]  (8 MHz)

# ---------------------------------------------------------------------------
#  Foreground wedge
# ---------------------------------------------------------------------------
# 0.1 h Mpc^-1 at h = 0.6766 (Pober et al. 2014 "moderate"; 21cmSense default)
wedge_buffer = 0.0677            # safety margin beyond horizon line  [Mpc^-1]

# ---------------------------------------------------------------------------
#  Power spectrum binning
# ---------------------------------------------------------------------------
n_bins_perp     = 20
n_bins_parallel = 20

# ---------------------------------------------------------------------------
#  Estimator formalism  (TODO.md P0)
# ---------------------------------------------------------------------------
# "coeval" (default) — the historical formalism, valid for the quasi-coeval
#     slab this pipeline has always run: lightcone slices spaced uniformly in
#     *redshift*, and a single global mean removed from each field.  Every
#     result produced before P0 landed used this, and it reproduces them
#     exactly.
#
# "lightcone" — the P0 formalism, required once Delta z is wide enough that
#     the box evolves along the line of sight:
#       P0.1  slices spaced uniformly in *comoving distance*, so the FFT does
#             not mis-assign k_parallel  (set here, in the lightconer)
#       P0.2  per-slice mean subtraction, so the LOS ramp in <T_b> and <SFR>
#             does not alias into low k_parallel  (set here for delta_gal;
#             src/analysis.py does the same for delta_T_b)
#       P0.3  one power spectrum per frequency sub-band, each at its own
#             effective redshift        (analysis stage)
#       P0.4  the sub-band width matched to the noise bandwidth
#                                        (analysis stage)
#
# The analysis stage reads `estimator` back from the HDF5 and defaults to it,
# so the two halves cannot silently disagree; `run_pipeline.py --estimator`
# overrides it for a cached-spectra re-run.
ESTIMATOR = "lightcone"          # coeval | lightcone

if ESTIMATOR not in ("coeval", "lightcone"):
    raise ValueError(f"ESTIMATOR must be 'coeval' or 'lightcone', got {ESTIMATOR!r}")

# Both follow from ESTIMATOR; override individually to isolate one P0 item.
LIGHTCONE_SAMPLING      = "comoving" if ESTIMATOR == "lightcone" else "redshift"
GALAXY_MEAN_SUBTRACTION = "per_slice" if ESTIMATOR == "lightcone" else "global"

# ---------------------------------------------------------------------------
#  Reproducibility
# ---------------------------------------------------------------------------
# 21cmFAST's initial-conditions seed.  Recorded in the run manifest and the
# HDF5 attrs, because two runs are only comparable if this matches.
RANDOM_SEED = 42

# ---------------------------------------------------------------------------
#  Compute resources
# ---------------------------------------------------------------------------
# 21cmFAST v4 defaults N_THREADS to 1 and nothing here used to override it, so
# the 2026-08-20 run spent its entire 38 minutes on a single core (user/real =
# 0.94 in outputs/21cm_pipeline_20260820_160108.log).  Resolution order:
#
#   N_THREADS env var  ->  SLURM_CPUS_PER_TASK  ->  os.cpu_count()  ->  1
#
# The SLURM variable is preferred over cpu_count() because on a shared node
# cpu_count() reports the whole machine, not this job's allocation.

N_THREADS = resolve_n_threads()

# MINIMIZE_MEMORY trades peak RAM for extra intermediate I/O inside the C
# backend.  Off by default in 21cmFAST; on here because the footprint-derived
# box (DIM = 768) is memory-bound, not I/O-bound.  Set False to compare.
MINIMIZE_MEMORY = True

# ---------------------------------------------------------------------------
#  Output path
# ---------------------------------------------------------------------------
OUTPUT_DIR  = "outputs"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "lightcone_data.h5")
RUN_DIR     = os.path.join(OUTPUT_DIR, "runs")   # per-run parameter manifests

# ===========================================================================
#  End of configuration
# ===========================================================================


# ===========================================================================
#  ⚠ SMOKE-TEST OVERRIDES — only active with --smoke-test ⚠
# ===========================================================================
# Nothing below runs unless the flag is set.  It reassigns the grid, the LOS
# slice floor, the k-binning and the output paths to values chosen purely to
# make the pipeline finish in seconds.  They are not science values and are
# never written into the configuration block above; they live in
# src/smoke_test.py with the production value each one stands in for.
#
# The output directory changes too, so a smoke run cannot overwrite a real
# run's lightcone_data.h5.

if SMOKE_TEST:
    from src.smoke_test import (
        SMOKE_OUTPUT_DIR,
        describe_overrides,
        override as _smoke_override,
    )

    HII_DIM         = _smoke_override("HII_DIM", HII_DIM)
    BOX_LEN         = _smoke_override("BOX_LEN", BOX_LEN)
    DIM             = _smoke_override("DIM", DIM)
    n_bins_perp     = _smoke_override("n_bins_perp", n_bins_perp)
    n_bins_parallel = _smoke_override("n_bins_parallel", n_bins_parallel)
    N_THREADS       = _smoke_override("N_THREADS", N_THREADS)

    OUTPUT_DIR  = SMOKE_OUTPUT_DIR
    OUTPUT_FILE = os.path.join(OUTPUT_DIR, "lightcone_data.h5")
    RUN_DIR     = os.path.join(OUTPUT_DIR, "runs")

    print("\n" + describe_overrides() + "\n")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ===========================================================================
# 1  Derived quantities
# ===========================================================================
# Compute the lightcone geometry from the configuration above: the comoving
# line-of-sight extent L_LOS, the number of LOS slices N_z (matched
# to the transverse cell size), and the redshift array for each slice.

def hubble_parameter(z):
    """H(z) for flat ΛCDM  [km s⁻¹ Mpc⁻¹]."""
    return HUBBLE_CONSTANT * np.sqrt(OMEGA_M_0 * (1 + z)**3 + (1 - OMEGA_M_0))


# ── Comoving distances to lightcone endpoints (using Planck18 for precision) ─
if HAS_21CMFAST:
    D_min = Planck18.comoving_distance(z_min).to(u.Mpc).value
    D_max = Planck18.comoving_distance(z_max).to(u.Mpc).value
else:
    D_min, _ = quad(lambda z_: SPEED_OF_LIGHT_KMS / hubble_parameter(z_), 0, z_min)
    D_max, _ = quad(lambda z_: SPEED_OF_LIGHT_KMS / hubble_parameter(z_), 0, z_max)

L_los          = D_max - D_min   # comoving LOS extent  [Mpc]
cell_size      = BOX_LEN / HII_DIM
hires_cell_size = BOX_LEN / DIM   # initial-conditions cell size  [Mpc]

# ── Mass resolution ──────────────────────────────────────────────────────────
# The mean comoving matter mass enclosed by one cell. This is the smallest mass
# element each grid can represent: DIM sets it for the initial-conditions /
# density field, HII_DIM for the coarse ionisation and 21 cm grids.
M_cell_hires = cell_mass(hires_cell_size, OMEGA_M_0, HUBBLE_CONSTANT)
M_cell_lores = cell_mass(cell_size,       OMEGA_M_0, HUBBLE_CONSTANT)

# ── Reference redshift (centre of lightcone, used for noise / wedge) ─────────
z_obs = 0.5 * (z_min + z_max)

# ── Number of LOS slices ──────────────────────────────────────────────────────
# Normally matched to the transverse cell size; override via minimum_los_slices
# for test runs where the z-range implies fewer natural slices.
minimum_los_slices = 100
if SMOKE_TEST:                       # see the override block above
    from src.smoke_test import override as _smoke_override
    minimum_los_slices = _smoke_override("minimum_los_slices", minimum_los_slices)
N_z = max(int(round(L_los / cell_size)), minimum_los_slices)

# Provisional per-slice redshifts (low-z → high-z), used only by the
# RectilinearLightconer on the uniform-in-redshift path and by the synthetic
# fallback.  On the real 21cmFAST path this is REPLACED by
# `lightcone.lightcone_redshifts` once the run returns, because the delivered
# slice count and spacing are the lightcone's to decide, not ours to predict.
lc_redshifts = np.linspace(z_min, z_max, N_z)

# ── Node redshifts for 21cmFAST (coeval snapshots driving the physics) ────────
# Use ~10 nodes per unit redshift for good accuracy
z_pad = 0.01

n_nodes        = max(int(round(10 * ((z_max + z_pad) - (z_min - z_pad)))), 5)
node_redshifts = np.linspace((z_max + z_pad), (z_min - z_pad), n_nodes)   # high-z → low-z

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"Box         : {BOX_LEN:.1f} Mpc,  {HII_DIM}³ cells  →  cell size = {cell_size:.2f} Mpc")
print(f"Footprint   : {SURVEY_AREA_DEG2:g} deg² at z = {SURVEY_Z_CENTRAL:g}  →  "
      f"BOX_LEN = {BOX_LEN:.1f} Mpc  (Euclid Deep Field Fornax)")
# Whether the configured lightcone actually uses the footprint-derived range
# is a property of z_min/z_max, not a fixed story about them -- the banner used
# to assert "[overridden by the smoke-test slab]" unconditionally, which was
# wrong twice over: it printed on production runs, and the override comes from
# the config block's own z_min/z_max, not from --smoke-test.
_survey_range_in_use = (
    abs(z_min - SURVEY_Z_MIN) < 1e-6 and abs(z_max - SURVEY_Z_MAX) < 1e-6
)
print(f"Survey LOS  : Δz = {SURVEY_DELTA_Z:g} (±{PHOTOZ_N_SIGMA}σ_z)  →  "
      f"z = {SURVEY_Z_MIN:.3f}–{SURVEY_Z_MAX:.3f}, "
      f"L_LOS = {SIM_BOX.los_depth:.1f} Mpc  "
      + ("[in use]" if _survey_range_in_use else
         f"[NOT USED — z_min/z_max override it with "
         f"{z_min:g}–{z_max:g} below]"))
print(f"Mass res.   : {M_cell_hires:.3e} M⊙/cell (DIM={DIM}, {hires_cell_size:.3f} Mpc)  |  "
      f"{M_cell_lores:.3e} M⊙/cell (HII_DIM={HII_DIM}, {cell_size:.2f} Mpc)")
print(f"Lightcone   : z = {z_min} → {z_max}   (reference z_obs = {z_obs})")
print(f"LOS extent  : {L_los:.1f} Mpc  →  N_z = {N_z} slices")
print(f"Box shape   : ({HII_DIM}, {HII_DIM}, {N_z})")
print(f"Node redshifts: {n_nodes} nodes from z={z_max} to z={z_min}")
print(f"Euclid      : M_UV < {M_UV_limit},  σ_z = {photoz_uncertainty}")
print(f"Threads     : N_THREADS = {N_THREADS}  (of {os.cpu_count()} visible), "
      f"MINIMIZE_MEMORY = {MINIMIZE_MEMORY}")
print(f"Estimator   : {ESTIMATOR}  (LOS sampling {LIGHTCONE_SAMPLING}, "
      f"delta_gal mean {GALAXY_MEAN_SUBTRACTION})")

# ── Pre-flight halo-catalogue cost ───────────────────────────────────────────
# Extrapolated from the measured 256 Mpc run (see src/provenance.py).  The
# catalogue scales with volume, so a modest-looking change in BOX_LEN is a
# large change in memory: 256 -> 486.33 Mpc is 6.9x, which is what killed the
# 2026-08-20 run with SIGSEGV.
cost = estimate_catalogue_cost(BOX_LEN, SAMPLER_MIN_MASS)
print(f"Sampler     : SAMPLER_MIN_MASS = {SAMPLER_MIN_MASS:.2e} M⊙  →  "
      f"{cost['sampler_retained_fraction']:.1%} of halos, "
      f"{cost['int32_headroom']:.2f}x INT_MAX")
print(f"Est. halos  : {cost['n_halos_lagrangian']:.3e} drawn "
      f"({cost['n_halos_perturbed']:.3e} after perturbation) in "
      f"{cost['volume_Mpc3']:.3e} Mpc³  →  {cost['catalogue_GB']:.1f} GB on disk, "
      f"~{cost['resident_GB']:.1f} GB resident while perturbing")

if cost["int32_headroom"] > 1.0:
    # Fatal, not a warning.  This previously printed and then ran anyway,
    # segfaulting ~38 minutes into the halo stage with exit code -11 and no
    # usable output.  There is no configuration of memory or thread count that
    # makes an over-length 32-bit index work, so there is nothing to be gained
    # by continuing.
    raise SystemExit(
        f"\n  *** ABORT: halo_coords would hold "
        f"{cost['n_halos_lagrangian'] * 3:.3e} elements, "
        f"{cost['int32_headroom']:.2f}x INT_MAX ({INT32_MAX:.3e}).\n"
        f"      21cmFAST indexes halo arrays with int, so this box WILL\n"
        f"      overflow and segfault regardless of available memory.\n"
        f"      BOX_LEN = {BOX_LEN:.2f} Mpc, SAMPLER_MIN_MASS = "
        f"{SAMPLER_MIN_MASS:.2e} M_sun "
        f"(retains {cost['sampler_retained_fraction']:.1%} of halos).\n"
        f"      Fix: raise SAMPLER_MIN_MASS to >= "
        f"{SAMPLER_MIN_MASS * cost['int32_headroom']:.2e}, or reduce BOX_LEN\n"
        f"      to <= {BOX_LEN / cost['int32_headroom'] ** (1/3):.1f} Mpc.\n"
        f"      Narrowing the redshift range will NOT help: the catalogue is\n"
        f"      drawn in the full BOX_LEN^3 cube, not the lightcone slab. ***\n"
    )
elif cost["int32_headroom"] > 0.5:
    print(f"  (halo_coords at {cost['int32_headroom']:.2f}x INT_MAX — "
          f"over half the 32-bit index range)")

# ── Free-disk pre-flight ─────────────────────────────────────────────────────
# The 21cmFAST cache keeps one halo catalogue PER NODE REDSHIFT, and n_nodes is
# proportional to the lightcone's redshift span.  Widening z therefore scales
# the cache linearly -- the opposite of the int32 guard above, whose headroom
# depends only on BOX_LEN and is indifferent to the redshift range.  A z =
# 6.5-7.5 run filled a scratch filesystem overnight with errno 28 partway
# through node 7.2833, after hours of compute, so this is checked up front.
#
# `run_lightcone` is not passed a cache argument, so py21cmfast defaults to
# OutputCache(direc=Path('.')) -- the CACHE LANDS IN THE CURRENT WORKING
# DIRECTORY.  See docs/HPC.md R3.  That is the filesystem measured here.
_cache_dir = os.path.abspath(os.getcwd())
_cache = estimate_cache_footprint(BOX_LEN, n_nodes, SAMPLER_MIN_MASS)
_free_gb = shutil.disk_usage(_cache_dir).free / 1e9

print(f"Cache est.  : {_cache['total_GB']:.0f}–{_cache['total_upper_GB']:.0f} GB "
      f"({n_nodes} nodes x {_cache['per_node_GB']:.1f} GB + {_cache['ics_GB']:.1f} GB ICs)"
      f"  |  {_free_gb:.0f} GB free on {_cache_dir}")

if _free_gb < _cache["total_GB"] * DISK_SAFETY_FACTOR:
    raise SystemExit(
        f"\n  *** ABORT: not enough free disk for the 21cmFAST cache.\n"
        f"      Need >= {_cache['total_GB'] * DISK_SAFETY_FACTOR:.0f} GB "
        f"({_cache['total_GB']:.0f} GB estimated x {DISK_SAFETY_FACTOR:g} safety), "
        f"and up to {_cache['total_upper_GB']:.0f} GB if\n"
        f"      PerturbHaloField is cached per node too.  "
        f"Free: {_free_gb:.0f} GB on\n"
        f"      {_cache_dir}\n"
        f"      The cache is {n_nodes} nodes x {_cache['per_node_GB']:.1f} GB "
        f"(one halo catalogue each) + {_cache['ics_GB']:.1f} GB of ICs.\n"
        f"      Fixes, in order of cheapness:\n"
        f"        1. Delete stale cache trees from earlier runs.  Every change\n"
        f"           to BOX_LEN / SAMPLER_MIN_MASS / seed starts a NEW hashed\n"
        f"           directory, so failed runs leave full-size orphans behind.\n"
        f"        2. Narrow z_min/z_max: n_nodes = round(10 x delta_z), so the\n"
        f"           cache scales linearly with the redshift span.\n"
        f"        3. Raise SAMPLER_MIN_MASS (3e8 keeps 99.44% of the star\n"
        f"           formation and cuts the cache to "
        f"{estimate_cache_footprint(BOX_LEN, n_nodes, 3e8)['total_GB']:.0f} GB).\n"
        f"        4. Run from a larger filesystem -- the cache follows the cwd. ***\n"
    )

# ── Run manifest ─────────────────────────────────────────────────────────────
# Written now, before the expensive stages, and rewritten after every one of
# them.  A run killed by a signal cannot flush stdout, but it leaves this file
# behind with `status: running` and `stage` naming where it died.
manifest = RunManifest.create(RUN_DIR, label="sim")
manifest.record("parameters", {
    "SURVEY_AREA_DEG2": SURVEY_AREA_DEG2,
    "SURVEY_Z_CENTRAL": SURVEY_Z_CENTRAL,
    "SURVEY_DELTA_Z": SURVEY_DELTA_Z,
    "PHOTOZ_N_SIGMA": PHOTOZ_N_SIGMA,
    "HII_DIM": HII_DIM,
    "BOX_LEN": BOX_LEN,
    "DIM": DIM,
    "z_min": z_min,
    "z_max": z_max,
    "minimum_los_slices": minimum_los_slices,
    "M_UV_limit": M_UV_limit,
    "M_UV_bright": M_UV_bright,
    "M_UV_faint": M_UV_faint,
    "photoz_uncertainty": photoz_uncertainty,
    "mean_galaxy_density": mean_galaxy_density,
    "GALAXY_WEIGHTING": GALAXY_WEIGHTING,
    "OMEGA_M_0": OMEGA_M_0,
    "HUBBLE_CONSTANT": HUBBLE_CONSTANT,
    "HERA_DISH_DIAMETER": HERA_DISH_DIAMETER,
    "integration_time": integration_time,
    "bandwidth": bandwidth,
    "wedge_buffer": wedge_buffer,
    "n_bins_perp": n_bins_perp,
    "n_bins_parallel": n_bins_parallel,
    "N_THREADS": N_THREADS,
    "MINIMIZE_MEMORY": MINIMIZE_MEMORY,
    "SMOKE_TEST": SMOKE_TEST,
    "ESTIMATOR": ESTIMATOR,
    "LIGHTCONE_SAMPLING": LIGHTCONE_SAMPLING,
    "GALAXY_MEAN_SUBTRACTION": GALAXY_MEAN_SUBTRACTION,
    "random_seed": RANDOM_SEED,
    "has_21cmfast": HAS_21CMFAST,
})
manifest.record("derived", {
    "cell_size_Mpc": cell_size,
    "hires_cell_size_Mpc": hires_cell_size,
    "M_cell_hires_Msun": M_cell_hires,
    "M_cell_lores_Msun": M_cell_lores,
    "L_los_Mpc": L_los,
    "N_z": N_z,
    "z_obs": z_obs,
    "n_nodes": n_nodes,
    "box_shape": [HII_DIM, HII_DIM, N_z],
    "survey_z_min": SURVEY_Z_MIN,
    "survey_z_max": SURVEY_Z_MAX,
    "survey_los_depth_Mpc": SIM_BOX.los_depth,
})
manifest.record("cost_estimate", cost)
manifest.record("outputs", {"lightcone_data": os.path.abspath(OUTPUT_FILE)})
print(f"Run manifest: {manifest.path}")


# ===========================================================================
# 2  Run the 21cmFASTv4 lightcone simulation
# ===========================================================================
# The v4 API uses RectilinearLightconer to specify the LOS slice redshifts and
# the quantities to store, then run_lightcone to execute the simulation.
#
# Compared to a coeval run, the lightcone:
#   - Evolves the physics across the full redshift range [z_min, z_max]
#   - Stores one coeval slice per LOS cell rather than a single snapshot
#   - Produces a rectangular (HII_DIM, HII_DIM, N_z) output volume

if HAS_21CMFAST:
    print("\nRunning 21cmFAST lightcone simulation …")
    print(f"  z = {z_min} → {z_max},  N_z = {N_z} slices")
    manifest.begin_stage("lightcone")

    inputs = p21c.InputParameters.from_template(
        ["simple"],
        random_seed=RANDOM_SEED,
    )

    if SMOKE_TEST:
        # The reduced box has to stay large enough for the *unmodified* source
        # model: 21cmFAST refuses BOX_LEN < 3 x R_BUBBLE_MAX, and does so
        # inside input validation, i.e. after the run has been queued and the
        # manifest written.  Fail here instead, naming the constraint.
        from src.smoke_test import check_box_supports_the_source_model

        _ok, _why = check_box_supports_the_source_model(
            BOX_LEN, float(inputs.astro_params.R_BUBBLE_MAX)
        )
        if not _ok:
            manifest.finish("failed")
            raise ValueError(_why)
        print(f"  Smoke pre-flight: BOX_LEN {BOX_LEN:g} Mpc >= 3 x "
              f"R_BUBBLE_MAX ({3 * float(inputs.astro_params.R_BUBBLE_MAX):g} Mpc) ✓")

    inputs = inputs.clone(
        node_redshifts=node_redshifts,
        simulation_options={
            "HII_DIM":          HII_DIM,
            "BOX_LEN":          BOX_LEN,
            "DIM":              DIM,
            "N_THREADS":        N_THREADS,
            # Not a physics choice — see the config block.  Left out of this
            # clone, the "simple" template's 1e8 overflows the 32-bit halo
            # index at the footprint-derived BOX_LEN.
            "SAMPLER_MIN_MASS": SAMPLER_MIN_MASS,
        },
        matter_options={
            "USE_INTERPOLATION_TABLES": "hmf-interpolation",
            # Trades peak RAM for intermediate I/O in the C backend.
            "MINIMIZE_MEMORY": MINIMIZE_MEMORY,
        },
    )

    _QUANTITIES = ("brightness_temp", "density", "neutral_fraction", "halo_sfr")

    if LIGHTCONE_SAMPLING == "comoving":
        # TODO.md P0.1 — `between_redshifts` builds
        # `lc_distances = arange(d_min, d_max + res, res)` internally, so the
        # slices are evenly spaced in comoving distance and the FFT's single
        # assumed cell size is correct to machine precision.  The slice count
        # follows from the resolution and is not ours to floor, so
        # `minimum_los_slices` does not apply on this path.
        lightconer = p21c.RectilinearLightconer.between_redshifts(
            min_redshift=z_min,
            max_redshift=z_max,
            resolution=cell_size * u.Mpc,
            quantities=_QUANTITIES,
        )
        print(f"  LOS sampling    : uniform in comoving distance "
              f"({cell_size:.4f} Mpc resolution)")
    else:
        # Historical path: uniform in redshift, hence *not* uniform in
        # comoving distance because dD/dz = c/H(z) varies across the box.
        lightconer = p21c.RectilinearLightconer(
            lc_redshifts=lc_redshifts,
            quantities=_QUANTITIES,
        )
        print(f"  LOS sampling    : uniform in redshift ({N_z} slices)")

    lightcone = p21c.run_lightcone(
        lightconer=lightconer,
        inputs=inputs,
        include_dvdr_in_tau21=False,
        apply_rsds=False,
    )
    brightness_temp_field = lightcone.lightcones["brightness_temp"]   # (HII_DIM, HII_DIM, N_z)
    density_field         = lightcone.lightcones["density"]
    neutral_fraction      = lightcone.lightcones["neutral_fraction"]

    # Update LOS geometry from the actual simulation output.
    #
    # ALL FOUR come from the lightcone, not from the planned values above.
    # `lc_redshifts` used to be left at the line-570
    # `np.linspace(z_min, z_max, N_z)`, which was wrong twice over on the
    # comoving path:
    #   * Length.  `between_redshifts` builds
    #     `arange(d_min, d_max + res, res)`, whose inclusive endpoint can
    #     return one more slice than `round(L_los / cell_size)` predicted --
    #     185 planned vs 186 delivered for z = 6.5-7.5, which the analysis
    #     stage rejected outright.
    #   * Values.  Slices uniform in comoving DISTANCE are not uniform in
    #     REDSHIFT, because dD/dz = c/H(z) varies across the box.  Even at a
    #     matching slice count the linspace was off by up to dz = 0.023
    #     (2.3 % of the span) mid-box.  That array sets each sub-band's
    #     effective redshift (TODO.md P0.3), so the error was silent and
    #     systematic wherever the counts happened to agree.
    N_z         = lightcone.n_slices
    L_los       = lightcone.lightcone_dimensions[2]                  # actual LOS comoving size  [Mpc]
    lc_dist_Mpc = lightcone.lightcone_distances.to(u.Mpc).value      # (N_z,) comoving distances
    lc_redshifts = np.asarray(
        getattr(lightcone.lightcone_redshifts, "value",
                lightcone.lightcone_redshifts),
        dtype=float,
    )                                                                # (N_z,)

    if lc_redshifts.shape != (N_z,) or lc_dist_Mpc.shape != (N_z,):
        manifest.finish("failed")
        raise ValueError(
            f"lightcone LOS axes disagree: n_slices={N_z}, "
            f"lc_redshifts={lc_redshifts.shape}, "
            f"lc_dist_Mpc={lc_dist_Mpc.shape}"
        )

    print(f"  Lightcone shape : {lightcone.shape}")
    print(f"  LOS extent      : {L_los:.1f} Mpc")
    print(f"  <x_HI>          : {np.mean(neutral_fraction):.3f}")
    print(
        f"  dT_b            : "
        f"[{brightness_temp_field.min():.1f}, "
        f"{brightness_temp_field.max():.1f}] mK"
    )

    # ── Check matter_options (preliminary sanity check) ────────────────────
    print("\nMatter options:")
    print(inputs.matter_options)

    # Smallest halo the stochastic halo sampler populates — this sets the mass
    # resolution of the halo catalogue (and hence of the galaxy field).
    sampler_min_mass = float(inputs.simulation_options.SAMPLER_MIN_MASS)

    manifest.record("results", {
        "lightcone_shape": list(lightcone.shape),
        "L_los_actual_Mpc": float(L_los),
        "N_z_actual": int(N_z),
        "mean_neutral_fraction": float(np.mean(neutral_fraction)),
        "brightness_temp_min_mK": float(brightness_temp_field.min()),
        "brightness_temp_max_mK": float(brightness_temp_field.max()),
        "sampler_min_mass_Msun": sampler_min_mass,
    })
    print(f"  Lightcone stage done in {manifest.end_stage():.1f} s")

else:
    # ==================================================================
    #  Synthetic lightcone fallback (21cmFAST not installed)
    # ==================================================================
    print("\nGenerating synthetic lightcone fields (21cmFAST not available) …")
    np.random.seed(42)

    N      = HII_DIM
    kx     = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
    ky     = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
    kz_los = np.fft.fftfreq(N_z, d=cell_size) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz_los, indexing="ij")
    K_magnitude = np.sqrt(KX**2 + KY**2 + KZ**2)
    K_magnitude[0, 0, 0] = 1.0

    # Matter overdensity from P(k) ∝ k^−2
    ps_amplitude = K_magnitude**(-2)
    ps_amplitude[0, 0, 0] = 0.0
    random_phases = (
        np.random.standard_normal((N, N, N_z))
        + 1j * np.random.standard_normal((N, N, N_z))
    )
    density_field = np.fft.ifftn(np.sqrt(ps_amplitude / 2) * random_phases).real
    density_field = (density_field - density_field.mean()) / density_field.std() * 0.5

    # Neutral fraction via inside-out reionisation
    k_bubble = 2 * np.pi / 10.0
    smoothed = np.fft.ifftn(
        np.fft.fftn(density_field) * np.exp(-0.5 * K_magnitude**2 / k_bubble**2)
    ).real
    neutral_fraction = (smoothed < np.percentile(smoothed, 50)).astype(float)

    # Brightness temperature with redshift evolution along LOS
    T0_per_slice = 27.0 * np.sqrt((1 + lc_redshifts) / 10.0)   # (N_z,)
    brightness_temp_field = (
        T0_per_slice[np.newaxis, np.newaxis, :]
        * neutral_fraction * (1.0 + density_field)
    )

    lc_dist_Mpc = D_min + (D_max - D_min) * np.linspace(0, 1, N_z)

    # No halo sampler in the synthetic fallback — no halo mass resolution.
    sampler_min_mass = np.nan

    print(f"  Shape  : {brightness_temp_field.shape}")
    print(f"  ⟨x_HI⟩ : {np.mean(neutral_fraction):.3f}")
    print(f"  δT_b   : [{brightness_temp_field.min():.1f}, {brightness_temp_field.max():.1f}] mK")


# ===========================================================================
# 3a  Extract actual 21cmFASTv4 halo catalogue at z_obs
# ===========================================================================

if HAS_21CMFAST:
    print(f"\nExtracting halo catalogue at z = {z_obs} ...")
    manifest.begin_stage("halo_catalog")

    # 1. Compute initial conditions using the same inputs as the lightcone
    initial_conditions = p21c.compute_initial_conditions(
        inputs=inputs,
    )

    # 2. Draw the Lagrangian halo catalogue at z_obs
    halo_catalog = p21c.determine_halo_catalog(
        redshift=z_obs,
        inputs=inputs,
        initial_conditions=initial_conditions,
    )

    # 3. Perturb halos to Eulerian positions at z_obs
    perturbed_halos = p21c.perturb_halo_catalog(
        inputs=inputs,
        initial_conditions=initial_conditions,
        halo_catalog=halo_catalog,
    )

    # This is the memory high-water mark of the whole script: the Lagrangian
    # catalogue and its perturbed copy are both resident, ~28 bytes per halo
    # each, plus the DIM^3 initial conditions.  At BOX_LEN = 486.33 Mpc that
    # is ~48 GB of catalogue on top of ~7.7 GB of ICs.  Nothing below reads
    # either input again, so drop them before touching the arrays.
    del halo_catalog, initial_conditions
    gc.collect()

    # 4. Pull useful arrays (convert 21cmFAST Array objects to NumPy)
    #
    # UNIT NOTE — py21cmfast v4 computes SFR internally in M_sun s^-1
    # (scaling_relations.c: sfr_mean = stellar_mass / (t_star * t_h),
    # where t_h = t_hubble(z) is in seconds).  The `.value` attribute
    # returns this raw internal value.  We convert to M_sun yr^-1 here
    # so that the HDF5 catalogue matches its documented unit and is
    # directly comparable to observational SFR calibrations.
    _SEC_PER_YR = 365.25 * 24 * 3600          # 3.15576e7 s yr^-1
    halo_masses    = np.asarray(perturbed_halos.halo_masses.value)
    halo_coords    = np.asarray(perturbed_halos.halo_coords.value)
    stellar_masses = np.asarray(perturbed_halos.stellar_masses.value)
    sfr_cat        = np.asarray(perturbed_halos.sfr.value) * _SEC_PER_YR

    print(f"  Number of halos: {perturbed_halos.n_halos}")

    # ── Preliminary diagnostics ───────────────────────────────────────────────
    print("\nCatalogue diagnostics:")
    print(f"  halo_coords shape = {halo_coords.shape}, dtype = {halo_coords.dtype}")
    print(f"  halo_masses shape = {halo_masses.shape}, dtype = {halo_masses.dtype}")
    print(f"  N_halos           = {halo_coords.shape[0]:,}")
    print(f"  coord min/max     = {halo_coords.min():.3f}, {halo_coords.max():.3f}")
    print(f"  mass min/max      = {halo_masses.min():.3e}, {halo_masses.max():.3e} Msun")

    positive_sfr = sfr_cat[np.isfinite(sfr_cat) & (sfr_cat > 0)]
    print(f"\n  SFR diagnostics:")
    print(f"  N halos with SFR > 0 = {len(positive_sfr):,}")
    if len(positive_sfr) > 0:
        print(f"  SFR min/max = {positive_sfr.min():.3e} - {positive_sfr.max():.3e} Msun/yr")
        print(f"  SFR median  = {np.median(positive_sfr):.3e} Msun/yr")

    manifest.record("results", {
        "n_halos": int(halo_coords.shape[0]),
        "n_halos_sfr_positive": int(len(positive_sfr)),
        "halo_mass_min_Msun": float(halo_masses.min()),
        "halo_mass_max_Msun": float(halo_masses.max()),
    })
    print(f"  Halo-catalogue stage done in {manifest.end_stage():.1f} s")

else:
    # Empty arrays when 21cmFAST is unavailable
    halo_masses    = np.array([])
    halo_coords    = np.array([]).reshape(0, 3)
    stellar_masses = np.array([])
    sfr_cat        = np.array([])


# ===========================================================================
# 3b  Construct the galaxy density field
# ===========================================================================
# The lightcone stores `halo_sfr`: the per-cell SFR density integrated over
# all halos, evaluated at each LOS slice. This naturally traces UV-bright
# galaxies as in the coeval case, but now evolves continuously along the LOS.
#
# Synthetic fallback: if 21cmFAST is unavailable, a biased Poisson-sampled
# galaxy field is generated from the density field, anti-correlated with the
# 21 cm signal.

cell_volume = cell_size**3

if HAS_21CMFAST and GALAXY_WEIGHTING in GALAXY_WEIGHTING_MODES:
    # ── Catalogue-based field: number- or luminosity-weighted ─────────────
    # Both modes go through galaxy_overdensity_from_catalogue(), which
    # applies the Euclid magnitude window and deposits the survivors onto a
    # (HII_DIM, HII_DIM, N_z) grid. The only difference is the per-halo
    # weight: 1 for "number", L_UV = SFR / kappa_UV for "luminosity".
    print(f"\nConstructing {GALAXY_WEIGHTING}-weighted galaxy field "
          f"from the halo catalogue …")

    # Convert halo_coords to Mpc using the same convention as
    # figures._halo_coords_mpc(), so this deposit and the diagnostic
    # projection maps bin identical positions.
    if halo_coords.size and halo_coords.max() <= HII_DIM + 1:
        halo_coords_mpc = halo_coords * cell_size
    else:
        halo_coords_mpc = halo_coords

    galaxy_overdensity, galaxy_selection = galaxy_overdensity_from_catalogue(
        coords=halo_coords_mpc,
        sfr=sfr_cat,
        halo_masses=halo_masses,
        box_len=BOX_LEN,
        n_perp=HII_DIM,
        n_los=N_z,
        # L_los, NOT BOX_LEN.  The 21 cm lightcone's N_z slices span L_los;
        # depositing the galaxies over BOX_LEN instead would give the two
        # fields the same array shape but different physical depths, so cell
        # j of one would sit at a different comoving distance from cell j of
        # the other -- a progressive line-of-sight misalignment reaching
        # ~133 Mpc by the far edge, which would destroy the cross-correlation
        # while leaving both auto-spectra looking perfectly healthy.
        # Because the coeval cell size (BOX_LEN/HII_DIM) equals the lightcone
        # slice spacing, this takes the first L_los-deep slab of the coeval
        # box at native resolution -- exactly N_z cells, no interpolation.
        los_extent=L_los,
        weighting=GALAXY_WEIGHTING,
        M_UV_faint=M_UV_faint,
        M_UV_bright=M_UV_bright,
    )

    print(f"  Euclid window : {M_UV_bright} < M_UV < {M_UV_faint}")
    print(f"  Halos SFR > 0 : {galaxy_selection.n_valid:,}")
    print(f"  Deposited     : {galaxy_selection.n_selected:,} galaxies")
    print(f"  Grid          : {galaxy_overdensity.shape}, "
          f"LOS extent {L_los:.1f} Mpc (matches the 21 cm lightcone)")
    print(f"  Galaxy δ      : [{galaxy_overdensity.min():.2f}, {galaxy_overdensity.max():.2f}]")
    print(f"  Shot-noise n̄  : {mean_galaxy_density:.2e} Mpc⁻³  (survey parameter)")

    if galaxy_selection.n_selected == 0:
        print("  WARNING: no halo passed the magnitude cut — δ_gal is identically zero.")

elif HAS_21CMFAST:
    # ── Default: the lightcone halo_sfr field ─────────────────────────────
    print("\nConstructing galaxy density field from lightcone halobox …")

    sfr_field = lightcone.lightcones["halo_sfr"]   # (HII_DIM, HII_DIM, N_z)

    # TODO.md P0.2 — <SFR> evolves along the lightcone, so normalising by a
    # single global scalar leaves a monotonic LOS ramp in delta_gal.
    # Normalising per slice removes it by construction, and additionally makes
    # delta_gal a true overdensity *at each redshift* rather than relative to
    # the band average, which is what the bias and the Kaiser factor assume.
    # Note the ramp itself sits at k_perp = 0, which the power-spectrum binning
    # already discards (see analysis.subtract_field_mean), so this changes the
    # binned spectra very little; it is a correctness fix, not a large
    # numerical one.  "global" remains the default.
    if GALAXY_MEAN_SUBTRACTION == "per_slice":
        mean_sfr = sfr_field.mean(axis=(0, 1), keepdims=True)
        galaxy_overdensity = np.where(
            mean_sfr > 0, sfr_field / np.where(mean_sfr > 0, mean_sfr, 1.0) - 1.0, 0.0
        )
        print(f"  delta_gal norm  : per-slice <SFR> "
              f"(range {float(mean_sfr.min()):.3e} - {float(mean_sfr.max()):.3e})")
    else:
        mean_sfr = sfr_field.mean()
        if mean_sfr > 0:
            galaxy_overdensity = sfr_field / mean_sfr - 1.0
        else:
            galaxy_overdensity = np.zeros_like(sfr_field)
        print(f"  delta_gal norm  : single global <SFR> = {float(mean_sfr):.3e}")

    print(f"  SFR density  : [{sfr_field.min():.2e}, {sfr_field.max():.2e}] (internal units)")
    print(f"  Galaxy δ     : [{galaxy_overdensity.min():.2f}, {galaxy_overdensity.max():.2f}]")
    print(f"  Shot-noise n̄ : {mean_galaxy_density:.2e} Mpc⁻³  (survey parameter)")

else:
    print("\nGenerating synthetic galaxy field …")

    expected_counts = (
        mean_galaxy_density
        * cell_volume
        * np.exp(galaxy_bias * density_field)
        * (1.0 - neutral_fraction)
    )
    galaxy_counts      = np.random.poisson(np.clip(expected_counts, 0, 50)).astype(float)
    mean_count         = galaxy_counts.mean()
    galaxy_overdensity = galaxy_counts / max(mean_count, 1e-10) - 1.0

    print(f"  Mean galaxy density : {mean_count / cell_volume:.2e} h³ Mpc⁻³")
    print(f"  Total galaxies      : {int(galaxy_counts.sum())}")


# ===========================================================================
# 4  Estimate Euclid galaxy bias from luminosity-selected halo population
# ===========================================================================

# Two estimators are available, and they do not agree:
#
#   Catalogue-based (preferred) — converts each perturbed halo's own SFR to
#     M_UV, applies the Euclid window, and averages the Sheth-Tormen bias over
#     the survivors. This inherits 21cmFAST's log-normal scatter (SIGMA_STAR =
#     0.25 dex, SIGMA_SFR_LIM = 0.19 dex), so the abundant low-mass halos that
#     scatter up into the magnitude window are counted.
#
#   Analytic HMF integral (fallback) — integrates the Sheth-Tormen bias over
#     the mass function weighted by a *mean, scatter-free* scaling relation.
#     With no scatter, only rare high-mass halos make the magnitude cut, so
#     this systematically overestimates the bias.
#
# The previous version used the analytic estimator with a hardcoded 100 Myr
# star-formation timescale, inconsistent with 21cmFAST's own
# t_STAR x t_H(z) ~ 570 Myr at z = 7. That combination produced b_g = 33.4
# (see docs/project_update.md). The timescale is now taken from
# star_formation_timescale(), and the catalogue estimator is authoritative
# whenever a halo catalogue exists, so the bias stored here matches the one
# the analysis stage reports.

# ── Euclid absolute UV magnitude limits ───────────────────────────────────────
# M_UV_bright / M_UV_faint are set with the other survey parameters above.

# ── 21cmFAST star-formation timescale at the reference redshift ──────────────
t_sf_yr = star_formation_timescale(
    z_obs, t_star=T_STAR_DEFAULT,
    hubble_constant=HUBBLE_CONSTANT, omega_m=OMEGA_M_0,
)
print(f"\nStar-formation timescale: t_sf = t_STAR x t_H(z_obs)"
      f" = {T_STAR_DEFAULT} x {t_sf_yr / T_STAR_DEFAULT / 1e9:.3f} Gyr"
      f" = {t_sf_yr / 1e6:.1f} Myr")

galaxy_bias_catalog = None
galaxy_bias_hmf     = None

# ---------------------------------------------------------------------------
#  Preferred estimator: the perturbed halo catalogue itself
# ---------------------------------------------------------------------------
if HAS_HMF and len(sfr_cat) > 0:
    print("\nEstimating galaxy bias from the perturbed halo catalogue …")

    euclid_selection = select_euclid_halos(
        sfr=sfr_cat,
        halo_masses=halo_masses,
        M_UV_faint=M_UV_faint,
        M_UV_bright=M_UV_bright,
    )

    print(f"  Euclid window   : {M_UV_bright} < M_UV < {M_UV_faint}")
    print(f"  SFR window      : {euclid_selection.SFR_min:.3e}"
          f" - {euclid_selection.SFR_max:.3e} Msun/yr")
    print(f"  Halos SFR > 0   : {euclid_selection.n_valid:,}")
    print(f"  Selected halos  : {euclid_selection.n_selected:,}")

    if euclid_selection.n_selected > 0:
        bias_estimate = effective_galaxy_bias(
            selection=euclid_selection,
            z_obs=z_obs,
            hubble_constant=HUBBLE_CONSTANT,
        )
        galaxy_bias_catalog = bias_estimate.mean_bias

        print(f"  Halo mass range : "
              f"{euclid_selection.halo_masses.min():.3e}"
              f" - {euclid_selection.halo_masses.max():.3e} Msun")
        print(f"  <b_g>           : {galaxy_bias_catalog:.3f}"
              f"  (range {bias_estimate.bias_min:.3f}"
              f" - {bias_estimate.bias_max:.3f})")
    else:
        print("  No halos passed the magnitude cut — falling back to the HMF integral.")

# ---------------------------------------------------------------------------
#  Fallback / cross-check: analytic HMF integral over the mean relation
# ---------------------------------------------------------------------------
if HAS_HMF:
    L_UV_min = Muv_to_Luv(M_UV_faint)    # faint limit
    L_UV_max = Muv_to_Luv(M_UV_bright)   # bright limit

    print(f"\nAnalytic HMF cross-check")
    print(f"  L_UV range = {L_UV_min:.3e} - {L_UV_max:.3e} erg s^-1 Hz^-1")

    # ── 21cmFAST-like mean stellar-halo relation ─────────────────────────────
    F_STAR10   = 0.05
    ALPHA_STAR = 0.5
    M_TURN     = 5e8

    OMEGA_B_0 = 0.049

    def f_star(M_h):
        """Stellar fraction f_*(M_h) with exponential turnover below M_TURN."""
        return F_STAR10 * (M_h / 1e10)**ALPHA_STAR * np.exp(-M_TURN / M_h)

    def stellar_mass_model(M_h):
        """Mean stellar mass for a halo of mass M_h [Msun]."""
        return f_star(M_h) * (OMEGA_B_0 / OMEGA_M_0) * M_h

    def sfr_model(M_h):
        """Mean SFR [Msun/yr] using 21cmFAST's own t_STAR x t_H timescale."""
        return stellar_mass_to_sfr(
            stellar_mass_model(M_h), z=z_obs, t_star=T_STAR_DEFAULT,
            hubble_constant=HUBBLE_CONSTANT, omega_m=OMEGA_M_0,
        )

    def uv_luminosity(M_h):
        """Mean UV luminosity [erg s^-1 Hz^-1] (Madau & Dickinson 2014)."""
        return sfr_to_Luv(sfr_model(M_h))

    # ── Halo mass function ───────────────────────────────────────────────────
    mf = MassFunction(
        z=z_obs,
        Mmin=7,
        Mmax=13,
        dlog10m=0.02,
    )

    M_h       = mf.m
    dndlog10m = mf.dndlog10m

    # ── Halo bias: Sheth-Tormen, from the squared peak height mf.nu ──────────
    halo_bias = sheth_tormen_bias(mf.nu)

    # ── Euclid selection on the mean relation ────────────────────────────────
    L_UV     = uv_luminosity(M_h)
    selected = (L_UV >= L_UV_min) & (L_UV <= L_UV_max)

    if selected.sum() < 2:
        print("  Too few mass bins satisfy the luminosity selection — "
              "skipping the analytic cross-check.")
    else:
        logM = np.log10(M_h)

        n_gal = simpson(dndlog10m[selected], x=logM[selected])
        galaxy_bias_hmf = simpson(
            dndlog10m[selected] * halo_bias[selected],
            x=logM[selected],
        ) / n_gal

        print(f"  n_gal = {n_gal:.3e} Mpc^-3")
        print(f"  Analytic (scatter-free) galaxy bias = {galaxy_bias_hmf:.2f}")
        print(f"  Selected halo mass range: "
              f"{M_h[selected].min():.2e} - {M_h[selected].max():.2e} Msun")

# ---------------------------------------------------------------------------
#  Adopt one value
# ---------------------------------------------------------------------------
if galaxy_bias_catalog is not None:
    galaxy_bias        = galaxy_bias_catalog
    galaxy_bias_method = "halo_catalog"
elif galaxy_bias_hmf is not None:
    galaxy_bias        = galaxy_bias_hmf
    galaxy_bias_method = "hmf_analytic"
else:
    galaxy_bias_method = "configured_default"
    print(f"\nUsing configured galaxy bias = {galaxy_bias}  (hmf not available).")

print(f"\nAdopted galaxy bias b_g = {galaxy_bias:.3f}  (method: {galaxy_bias_method})")
if galaxy_bias_catalog is not None and galaxy_bias_hmf is not None:
    print(f"  Scatter-free analytic estimate would give {galaxy_bias_hmf:.3f} "
          f"({galaxy_bias_hmf / galaxy_bias_catalog:.1f}x higher) — the "
          f"difference is 21cmFAST's log-normal SFR scatter.")


# ===========================================================================
# 5  Kaiser redshift-space distortions
# ===========================================================================
# The galaxy overdensity is boosted along the LOS by coherent infall
# (Kaiser 1987):
#
#   δ_gal^(s)(k) = (1 + β μ²) δ_gal(k),   β = f/b
#
# The growth rate f and bias b are evaluated at z_obs (midpoint of the
# lightcone). The non-cubic box shape (HII_DIM, HII_DIM, N_z) is accounted
# for in the wavenumber grids. RSD must be applied in Fourier space.

# ── Growth rate and bias at z_obs ─────────────────────────────────────────
Omega_m_z   = OMEGA_M_0 * (1 + z_obs)**3 / (
    OMEGA_M_0 * (1 + z_obs)**3 + (1 - OMEGA_M_0)
)
growth_rate = Omega_m_z ** 0.55           # f(z)
beta_rsd    = growth_rate / galaxy_bias   # β = f / b

print(f"\nReference redshift z_obs = {z_obs}")
print(f"Ω_m(z_obs)  = {Omega_m_z:.4f}")
print(f"Growth rate f  = {growth_rate:.4f}")
print(f"Galaxy bias b  = {galaxy_bias}")
print(f"β = f/b        = {beta_rsd:.4f}")

# ── Wavenumber grids for non-cubic lightcone box ──────────────────────────
#   Transverse: N × N cells, spacing = cell_size [Mpc]
#   LOS       : N_z cells,   spacing = L_los / N_z [Mpc]
N       = HII_DIM
dz_cell = L_los / N_z   # LOS cell size  [Mpc]

# ── Geometry self-check ──────────────────────────────────────────────────────
# The two fields are cross-correlated cell-by-cell, so they must share a shape
# AND a physical extent.  Matching shapes alone is not enough: a galaxy field
# deposited over BOX_LEN and a 21 cm lightcone spanning L_los have identical
# arrays sitting at different comoving depths, which silently destroys the
# cross-power while leaving both auto-spectra intact.  Fail here rather than
# ship a file whose two fields disagree about where they are.
if galaxy_overdensity.shape != brightness_temp_field.shape:
    raise ValueError(
        f"galaxy_overdensity {galaxy_overdensity.shape} and "
        f"brightness_temp_field {brightness_temp_field.shape} differ; the "
        f"cross-power requires a shared grid"
    )

kx = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
ky = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
kz = np.fft.fftfreq(N_z, d=dz_cell)   * 2 * np.pi
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

k_magnitude = np.sqrt(KX**2 + KY**2 + KZ**2)
k_magnitude[0, 0, 0] = 1.0   # avoid division by zero at DC

# μ = k_∥ / |k|  (LOS along z-axis)
mu = KZ / k_magnitude
mu[0, 0, 0] = 0.0

# Apply Kaiser boost in Fourier space
kaiser_factor          = 1.0 + beta_rsd * mu**2
galaxy_overdensity_k   = np.fft.fftn(galaxy_overdensity)
galaxy_overdensity_rsd = np.fft.ifftn(galaxy_overdensity_k * kaiser_factor).real

print(f"\nBefore RSD  →  δ_gal std = {galaxy_overdensity.std():.3f}")
print(f"After  RSD  →  δ_gal std = {galaxy_overdensity_rsd.std():.3f}")
print(f"Max Kaiser boost (μ=1)   = {1 + beta_rsd:.3f}×")

galaxy_overdensity = galaxy_overdensity_rsd


# ===========================================================================
# 6  Save all outputs to HDF5
# ===========================================================================
# All arrays and metadata are saved for use in the plotting and analysis
# notebooks (Parts 2 and 3). This avoids re-running the expensive simulation.

print(f"\nSaving outputs to {OUTPUT_FILE} …")
manifest.begin_stage("write_hdf5")

with h5py.File(OUTPUT_FILE, "w") as f:

    # ── Simulation fields ─────────────────────────────────────────────────
    f.create_dataset("brightness_temp_field", data=brightness_temp_field,
                     compression="gzip", compression_opts=4)
    f.create_dataset("density_field",         data=density_field,
                     compression="gzip", compression_opts=4)
    f.create_dataset("neutral_fraction",      data=neutral_fraction,
                     compression="gzip", compression_opts=4)
    f.create_dataset("galaxy_overdensity",    data=galaxy_overdensity,
                     compression="gzip", compression_opts=4)

    # ── LOS geometry ─────────────────────────────────────────────────────
    f.create_dataset("lc_redshifts", data=lc_redshifts)
    f.create_dataset("lc_dist_Mpc",  data=lc_dist_Mpc)

    # ── Halo catalogue ────────────────────────────────────────────────────
    hc = f.create_group("halo_catalog")
    hc.create_dataset("halo_masses",    data=halo_masses)
    hc.create_dataset("halo_coords",    data=halo_coords)
    hc.create_dataset("stellar_masses", data=stellar_masses)
    hc.create_dataset("sfr",            data=sfr_cat)

    # ── Scalar metadata (read back as attrs in the notebooks) ─────────────
    f.attrs["HII_DIM"]             = HII_DIM
    f.attrs["DIM"]                 = DIM
    f.attrs["BOX_LEN"]             = BOX_LEN
    f.attrs["N_z"]                 = N_z
    f.attrs["L_los"]               = L_los
    f.attrs["cell_size"]           = cell_size
    f.attrs["hires_cell_size"]     = hires_cell_size
    f.attrs["M_cell_hires"]        = M_cell_hires
    f.attrs["M_cell_lores"]        = M_cell_lores
    f.attrs["sampler_min_mass"]    = sampler_min_mass
    f.attrs["z_min"]               = z_min
    f.attrs["z_max"]               = z_max
    f.attrs["z_obs"]               = z_obs
    f.attrs["galaxy_bias"]         = galaxy_bias
    f.attrs["galaxy_bias_method"]  = galaxy_bias_method
    if galaxy_bias_hmf is not None:
        f.attrs["galaxy_bias_hmf_analytic"] = galaxy_bias_hmf
    f.attrs["t_STAR"]              = T_STAR_DEFAULT
    f.attrs["sfr_timescale_yr"]    = t_sf_yr
    f.attrs["beta_rsd"]            = beta_rsd
    f.attrs["smoke_test"]              = SMOKE_TEST
    f.attrs["estimator"]               = ESTIMATOR
    f.attrs["lightcone_sampling"]      = LIGHTCONE_SAMPLING
    f.attrs["galaxy_mean_subtraction"] = GALAXY_MEAN_SUBTRACTION
    f.attrs["mean_galaxy_density"] = mean_galaxy_density
    f.attrs["galaxy_weighting"]    = GALAXY_WEIGHTING
    f.attrs["photoz_uncertainty"]  = photoz_uncertainty
    f.attrs["M_UV_limit"]          = M_UV_limit
    # ── Survey footprint provenance of the box geometry ──────────────────
    f.attrs["survey_area_deg2"]    = SURVEY_AREA_DEG2
    f.attrs["survey_z_central"]    = SURVEY_Z_CENTRAL
    f.attrs["survey_delta_z"]      = SURVEY_DELTA_Z
    f.attrs["photoz_n_sigma"]      = PHOTOZ_N_SIGMA
    f.attrs["survey_z_min"]        = SURVEY_Z_MIN
    f.attrs["survey_z_max"]        = SURVEY_Z_MAX
    f.attrs["survey_los_depth"]    = SIM_BOX.los_depth
    f.attrs["survey_field"]        = "Euclid Deep Field Fornax"
    f.attrs["OMEGA_M_0"]           = OMEGA_M_0
    f.attrs["HUBBLE_CONSTANT"]     = HUBBLE_CONSTANT
    f.attrs["SPEED_OF_LIGHT_KMS"]  = SPEED_OF_LIGHT_KMS
    f.attrs["SPEED_OF_LIGHT_MPS"]  = SPEED_OF_LIGHT_MPS
    f.attrs["F_21_MHZ"]            = F_21_MHZ
    f.attrs["F_21_HZ"]             = F_21_HZ
    f.attrs["HERA_DISH_DIAMETER"]  = HERA_DISH_DIAMETER
    f.attrs["integration_time"]    = integration_time
    f.attrs["bandwidth"]           = bandwidth
    f.attrs["wedge_buffer"]        = wedge_buffer
    f.attrs["n_bins_perp"]         = n_bins_perp
    f.attrs["n_bins_parallel"]     = n_bins_parallel
    # ── Run provenance: ties this file to outputs/runs/sim_<run_id>.json ──
    f.attrs["random_seed"]         = RANDOM_SEED
    f.attrs["n_threads"]           = N_THREADS
    f.attrs["minimize_memory"]     = MINIMIZE_MEMORY
    f.attrs["run_id"]              = manifest.data["run_id"]
    f.attrs["run_manifest"]        = manifest.path

print("  Datasets saved:")
print(f"    brightness_temp_field  : {brightness_temp_field.shape}")
print(f"    density_field          : {density_field.shape}")
print(f"    neutral_fraction       : {neutral_fraction.shape}")
print(f"    galaxy_overdensity     : {galaxy_overdensity.shape}")
print(f"    lc_redshifts           : {lc_redshifts.shape}")
print(f"    lc_dist_Mpc            : {lc_dist_Mpc.shape}")
print(f"    halo_catalog/halo_masses    : {halo_masses.shape}")
print(f"    halo_catalog/halo_coords    : {halo_coords.shape}")
print(f"    halo_catalog/stellar_masses : {stellar_masses.shape}")
print(f"    halo_catalog/sfr            : {sfr_cat.shape}")

manifest.end_stage()

# ── Close the run manifest ───────────────────────────────────────────────────
manifest.record("results", {
    "galaxy_bias": float(galaxy_bias),
    "galaxy_bias_method": galaxy_bias_method,
    "galaxy_bias_hmf_analytic": (
        None if galaxy_bias_hmf is None else float(galaxy_bias_hmf)
    ),
    "beta_rsd": float(beta_rsd),
    "sfr_timescale_yr": float(t_sf_yr),
})
manifest.record("outputs", {
    "lightcone_data_GB": os.path.getsize(OUTPUT_FILE) / 1e9,
})
manifest.finish("complete")

print(f"\nSimulation complete. Output written to: {OUTPUT_FILE}")
print(f"Run manifest written to: {manifest.path}")
print("Next steps:")
print("  jupyter notebook notebooks/plot_fields.ipynb   # Part 2: field plots")
print("  jupyter notebook notebooks/analysis.ipynb      # Part 3: power spectra & SNR")
