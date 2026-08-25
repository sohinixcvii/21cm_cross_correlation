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

# sigma_z is the *absolute* photometric redshift error, not sigma_z/(1+z):
# radial_smearing_length() computes sigma_r = c sigma_z / H(z) directly.
# 0.45 at z_obs = 7 corresponds to sigma_z/(1+z) = 0.056, consistent with the
# Euclid photometric requirement sigma_z/(1+z) < 0.05. The previous value of
# 0.059 was the fractional quantity used as if it were absolute, which
# understated sigma_r by a factor ~7.6.
# Defined here, not in the Euclid block below, because it now also sets the
# line-of-sight depth of the box.
photoz_uncertainty = 0.45   # sigma_z photometric redshift error (absolute)

# CHOICE, not a default: how many sigma_z of photo-z scatter the box spans.
# PHOTOZ_N_SIGMA = 1 -> delta_z = 0.90, z = 6.55-7.45, L_los = 315.6 Mpc
# PHOTOZ_N_SIGMA = 2 -> delta_z = 1.80, z = 6.10-7.90, L_los = 634.9 Mpc
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
# 6.55/7.45 (delta_z = 0.90, L_LOS = 315.6 Mpc) from the photo-z depth above.
# That is the range this forecast *should* run once TODO.md P0 lands; the
# transverse box size is already sized for it. Until then the slab below
# deliberately overrides it, so only BOX_LEN/HII_DIM/DIM are footprint-driven.
SURVEY_Z_MIN = SIM_BOX.z_min   # 6.55 — survey-derived, not used yet
SURVEY_Z_MAX = SIM_BOX.z_max   # 7.45 — survey-derived, not used yet

z_min = 6.995          # nearest redshift (low-z end of lightcone)
z_max = 7.005          # farthest redshift (high-z end)

# ---------------------------------------------------------------------------
#  Euclid-like survey parameters
# ---------------------------------------------------------------------------
M_UV_limit          = -18    # UV absolute magnitude cut
# Absolute UV magnitude window (more negative = brighter). Defined here
# because both the galaxy-field construction (section 3b) and the bias
# estimate (section 4) select on it.
M_UV_bright         = -22
M_UV_faint          = M_UV_limit
# photoz_uncertainty (sigma_z = 0.45, absolute) is set in the survey-footprint
# block above, because it now also sets the line-of-sight depth of the box.
mean_galaxy_density = 3e-3   # n_bar  [h^3 Mpc^-3]

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
#     (L_UV = SFR / kappa_UV, Madau & Dickinson 2014).
#
# The two catalogue modes are interchangeable: identical grid, identical
# normalisation, identical downstream handling.  They are built from the
# *coeval* catalogue at z_obs, which spans BOX_LEN along the LOS rather than
# L_los, so they do not carry the lightcone's redshift evolution.
GALAXY_WEIGHTING = "lightcone_sfr"   # lightcone_sfr | number | luminosity

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
ESTIMATOR = "coeval"          # coeval | lightcone

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

# Redshift and distance for each LOS slice (low-z → high-z)
lc_redshifts = np.linspace(z_min, z_max, N_z)

# ── Node redshifts for 21cmFAST (coeval snapshots driving the physics) ────────
# Use ~10 nodes per unit redshift for good accuracy
n_nodes        = max(int(round(10 * (z_max - z_min))), 5)
node_redshifts = np.linspace(z_max, z_min, n_nodes)   # high-z → low-z

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"Box         : {BOX_LEN:.1f} Mpc,  {HII_DIM}³ cells  →  cell size = {cell_size:.2f} Mpc")
print(f"Footprint   : {SURVEY_AREA_DEG2:g} deg² at z = {SURVEY_Z_CENTRAL:g}  →  "
      f"BOX_LEN = {BOX_LEN:.1f} Mpc  (Euclid Deep Field Fornax)")
print(f"Survey LOS  : Δz = {SURVEY_DELTA_Z:g} (±{PHOTOZ_N_SIGMA}σ_z)  →  "
      f"z = {SURVEY_Z_MIN:.2f}–{SURVEY_Z_MAX:.2f}, L_LOS = {SIM_BOX.los_depth:.1f} Mpc "
      f"[overridden by the smoke-test slab below]")
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
cost = estimate_catalogue_cost(BOX_LEN)
print(f"Est. halos  : {cost['n_halos_lagrangian']:.3e} drawn "
      f"({cost['n_halos_perturbed']:.3e} after perturbation) in "
      f"{cost['volume_Mpc3']:.3e} Mpc³  →  {cost['catalogue_GB']:.1f} GB on disk, "
      f"~{cost['resident_GB']:.1f} GB resident while perturbing")

if cost["int32_headroom"] > 1.0:
    print(
        f"\n  *** WARNING: halo_coords would hold "
        f"{cost['n_halos_lagrangian'] * 3:.3e} elements, "
        f"{cost['int32_headroom']:.2f}x INT_MAX ({INT32_MAX:.3e}).\n"
        f"      21cmFAST indexes halo arrays with int; this box may overflow\n"
        f"      regardless of available memory.  Reduce BOX_LEN or raise\n"
        f"      SAMPLER_MIN_MASS before committing cluster time. ***\n"
    )
elif cost["int32_headroom"] > 0.5:
    print(f"  (halo_coords at {cost['int32_headroom']:.2f}x INT_MAX — "
          f"over half the 32-bit index range)")

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
            "HII_DIM":   HII_DIM,
            "BOX_LEN":   BOX_LEN,
            "DIM":       DIM,
            "N_THREADS": N_THREADS,
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

    # Update LOS geometry from the actual simulation output
    N_z         = lightcone.n_slices
    L_los       = lightcone.lightcone_dimensions[2]                  # actual LOS comoving size  [Mpc]
    lc_dist_Mpc = lightcone.lightcone_distances.to(u.Mpc).value      # (N_z,) comoving distances

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
        los_extent=BOX_LEN,
        weighting=GALAXY_WEIGHTING,
        M_UV_faint=M_UV_faint,
        M_UV_bright=M_UV_bright,
    )

    print(f"  Euclid window : {M_UV_bright} < M_UV < {M_UV_faint}")
    print(f"  Halos SFR > 0 : {galaxy_selection.n_valid:,}")
    print(f"  Deposited     : {galaxy_selection.n_selected:,} galaxies")
    print(f"  Grid          : {galaxy_overdensity.shape}, "
          f"LOS extent {BOX_LEN:.1f} Mpc (coeval box, not L_los)")
    print(f"  Galaxy δ      : [{galaxy_overdensity.min():.2f}, {galaxy_overdensity.max():.2f}]")
    print(f"  Shot-noise n̄  : {mean_galaxy_density:.2e} h³ Mpc⁻³  (survey parameter)")

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
    print(f"  Shot-noise n̄ : {mean_galaxy_density:.2e} h³ Mpc⁻³  (survey parameter)")

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
