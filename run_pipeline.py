#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_pipeline.py — End-to-end driver for the 21 cm × galaxy cross-correlation
=============================================================================
HERA × Euclid (lightcone)

Runs the whole workflow in one command:

    Stage 1  simulation  — optionally run ``run_simulation.py``
                           (21cmFAST lightcone → ``outputs/lightcone_data.h5``)
    Stage 2  analysis    — 2D cylindrical power spectra, photo-z damping,
                           foreground-wedge excision, HERA noise, SNR, and the
                           Euclid-selected effective galaxy bias
    Stage 3  figures     — every figure from ``notebooks/plot_fields.ipynb``
                           and ``notebooks/analysis.ipynb``, written to
                           ``outputs/figures/``
    Stage 4  summary     — key numbers printed and written as JSON

Each stage can be run fresh or served from stored results, so the expensive
21cmFAST run happens only when it must:

    --sim      auto | force | skip     (auto: run only if the HDF5 is missing)
    --analysis auto | force | skip     (auto: recompute only if the cache is
                                        missing or older than the simulation)

Usage
-----
    conda activate 21cmfast

    python run_pipeline.py                     # analyse stored results, plot
    python run_pipeline.py --sim force         # re-run the simulation first
    python run_pipeline.py --plots power snr   # only the k-space figures
    python run_pipeline.py --max-halos 5000000 # cap catalogue memory

    bash submit_job.sh --sim force             # same, with timing + email

Outputs
-------
    outputs/lightcone_data.h5      simulation fields + halo catalogue
    outputs/analysis_products.h5   cached power spectra
    outputs/figures/*.png          all figures
    outputs/pipeline_summary.json  scalar results

References
----------
Davies, Mesinger & Murray (2025) — arXiv:2504.17254
Gagnon-Hartman, Davies & Mesinger (2025) — arXiv:2502.20447
La Plante et al. (2023) — arXiv:2205.09770
Euclid Collaboration (2022) — arXiv:2108.01201
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import warnings
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

# Repo root on sys.path so ``src.*`` imports work from any working directory.
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from src import analysis, figures                       # noqa: E402
from src.dataio import (                                # noqa: E402
    SimulationData,
    load_power_spectra,
    load_subband_power_spectra,
    SubbandPowerSpectra,
    load_simulation,
    products_are_stale,
    save_power_spectra,
    save_subband_power_spectra,
    save_uncertainty_budget,
)

warnings.filterwarnings("ignore")

# ── Figure groups ─────────────────────────────────────────────────────────
#   fields   lightcone field slices (Part 2)
#   halos    halo catalogue positions, masses, SFR scaling relations (Part 2)
#   scaling  UVLF, M_star–M_UV, star-forming main sequence (Part 2)
#   power    2D cylindrical power spectra (Part 3)
#   snr      per-mode SNR and photo-z damped cross-power (Part 3)
#   budget   uncertainty-budget breakdown: damping, sigma terms, wedge (Part 3)
#   bias     Euclid-selected halo masses and b_h(M, z) (Part 3)
#   euclid   post-Euclid-cut catalogue, galaxy overdensity, and its overlay
#            on the 21 cm field (Part 2/3)
PLOT_GROUPS = (
    "fields", "halos", "scaling", "euclid", "power", "snr", "budget", "bias",
)

DEFAULT_DATA = os.path.join("outputs", "lightcone_data.h5")
DEFAULT_PRODUCTS = os.path.join("outputs", "analysis_products.h5")
DEFAULT_FIGDIR = os.path.join("outputs", "figures")
DEFAULT_SUMMARY = os.path.join("outputs", "pipeline_summary.json")

SEPARATOR = "=" * 72


# ===========================================================================
#  Console helpers
# ===========================================================================

def banner(title: str, quiet: bool = False) -> None:
    """
    Print a stage banner.

    Parameters
    ----------
    title : str
        Stage title.
    quiet : bool, optional
        Suppress output.
    """
    if not quiet:
        print(f"\n{SEPARATOR}\n  {title}\n{SEPARATOR}", flush=True)


def log(message: str, quiet: bool = False) -> None:
    """
    Print a progress message.

    Parameters
    ----------
    message : str
        Text to print.
    quiet : bool, optional
        Suppress output.
    """
    if not quiet:
        print(message, flush=True)


# ===========================================================================
#  Stage 1 — simulation
# ===========================================================================

def run_simulation_stage(
    mode: str,
    data_path: str,
    script: str,
    smoke_test: bool = False,
    quiet: bool = False,
) -> bool:
    """
    Optionally run ``run_simulation.py`` as a subprocess.

    Parameters
    ----------
    mode : {'auto', 'force', 'skip'}
        ``auto`` runs the simulation only when ``data_path`` is missing;
        ``force`` always runs it; ``skip`` never does.
    data_path : str
        Expected simulation output, ``outputs/lightcone_data.h5``.
    script : str
        Path to the simulation script.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    bool
        True if the simulation was executed.

    Raises
    ------
    FileNotFoundError
        If the script is missing, or if the simulation was skipped and no
        stored output exists.
    RuntimeError
        If the simulation exits with a non-zero status.
    """
    banner("STAGE 1 — SIMULATION", quiet)

    data_exists = os.path.exists(data_path)

    if mode == "skip" or (mode == "auto" and data_exists):
        if not data_exists:
            raise FileNotFoundError(
                f"no stored simulation output at {data_path}, and --sim=skip.\n"
                "Run with --sim force to generate it."
            )
        size_gb = os.path.getsize(data_path) / 1e9
        stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(data_path)))
        log(f"  Using stored simulation output: {data_path}", quiet)
        log(f"  ({size_gb:.2f} GB, written {stamp})", quiet)
        return False

    if not os.path.exists(script):
        raise FileNotFoundError(f"simulation script not found: {script}")

    log(f"  Running {script} (this is the expensive stage) …", quiet)
    start = time.time()

    # -u on the child, always.  Its stdout is block-buffered when the
    # pipeline's own output is redirected to a file, and a child killed by a
    # signal never flushes — which is how the 2026-08-20 SIGSEGV produced a
    # log with no indication of which stage had failed.
    command = [sys.executable, "-u", script]
    if smoke_test:
        # The child owns its own reduced configuration; see src/smoke_test.py.
        command.append("--smoke-test")
    result = subprocess.run(command, cwd=REPO_ROOT)

    if result.returncode != 0:
        raise RuntimeError(
            f"{script} failed with exit code {result.returncode}"
        )

    log(f"  Simulation finished in {time.time() - start:.1f} s", quiet)
    return True


# ===========================================================================
#  Stage 2 — analysis
# ===========================================================================

def power_spectra_stage(
    mode: str,
    data: SimulationData,
    products_path: str,
    estimator: str = "coeval",
    subband_bandwidth: float = 8e6,
    quiet: bool = False,
):
    """
    Compute the 2D cylindrical power spectra, or load them from cache.

    Parameters
    ----------
    mode : {'auto', 'force', 'skip'}
        ``auto`` recomputes only when the cache is missing or stale;
        ``force`` always recomputes; ``skip`` requires a usable cache.
    data : SimulationData
        Loaded simulation.
    products_path : str
        Cache file for the spectra.
    estimator : {'coeval', 'lightcone'}, optional
        Which formalism to measure with.  ``'coeval'`` (default) takes one
        FFT over the whole box with a global mean subtracted — valid while
        the box is quasi-coeval, and what every earlier result used.
        ``'lightcone'`` applies `TODO.md` P0.2–P0.4: per-slice means, a
        Blackman-Harris taper, and one spectrum per frequency sub-band.
    subband_bandwidth : float, optional
        Target per-band bandwidth [Hz] for the lightcone estimator.  Match it
        to the noise model's bandwidth — that match is P0.4.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    PowerSpectra
        The three spectra on a shared ``(k_perp, k_parallel)`` grid.  Under
        the lightcone estimator this is the band whose effective redshift is
        closest to ``z_obs``, so downstream figures have one representative
        set; the per-band spectra are in the second return value.
    SubbandPowerSpectra or None
        The per-band spectra and geometry, or ``None`` under the coeval
        estimator.
    bool
        True if the spectra were recomputed (rather than loaded).

    Raises
    ------
    FileNotFoundError
        If ``mode='skip'`` and no cache exists.
    """
    stale = products_are_stale(products_path, data.path)

    if mode == "skip" or (mode == "auto" and not stale):
        if not os.path.exists(products_path):
            raise FileNotFoundError(
                f"no cached analysis products at {products_path}, and "
                "--analysis=skip.\nRun with --analysis force to compute them."
            )
        subbands, _ = load_subband_power_spectra(products_path)
        spectra, attrs = load_power_spectra(products_path)
        log(f"  Loaded cached power spectra: {products_path}", quiet)
        if subbands is not None:
            log(f"  Cache holds {subbands.n_bands} sub-bands "
                f"(z_eff {subbands.z_effective.min():.3f}–"
                f"{subbands.z_effective.max():.3f})", quiet)
        if stale:
            log("  WARNING: cache is older than the simulation output "
                "(use --analysis force to refresh).", quiet)
        return spectra, subbands, False

    n_bins_perp = int(data.get("n_bins_perp", 20))
    n_bins_parallel = int(data.get("n_bins_parallel", 20))
    start = time.time()

    if estimator == "lightcone":
        # TODO.md P0.3/P0.4 — one spectrum per frequency sub-band, each at its
        # own effective redshift and its own bandwidth.
        log("  Computing sub-band power spectra (lightcone estimator) …", quiet)
        bands, geometry = analysis.compute_subband_power_spectra(
            brightness_temp_field=data.brightness_temp_field,
            galaxy_overdensity=data.galaxy_overdensity,
            lc_redshifts=data.lc_redshifts,
            lc_dist_Mpc=data.lc_dist_Mpc,
            box_len_perp=data.BOX_LEN,
            bandwidth_hz=subband_bandwidth,
            n_bins_perp=n_bins_perp,
            n_bins_parallel=n_bins_parallel,
            f_21_hz=data.get("F_21_HZ", 1420.405e6),
        )
        subbands = SubbandPowerSpectra(
            bands=bands,
            z_effective=geometry.z_effective,
            z_min=geometry.z_min,
            z_max=geometry.z_max,
            frequency_min_hz=geometry.frequency_min_hz,
            frequency_max_hz=geometry.frequency_max_hz,
            bandwidth_hz=geometry.bandwidth_hz,
            los_length_mpc=geometry.los_length_mpc,
            n_slices=geometry.n_slices,
            index_ranges=geometry.index_ranges,
        )
        for index in range(subbands.n_bands):
            log(f"    band {index}: z = {geometry.z_min[index]:.3f}–"
                f"{geometry.z_max[index]:.3f} (z_eff {geometry.z_effective[index]:.3f}), "
                f"{geometry.n_slices[index]} slices, "
                f"{geometry.bandwidth_hz[index] / 1e6:.2f} MHz, "
                f"L_los {geometry.los_length_mpc[index]:.1f} Mpc", quiet)

        save_subband_power_spectra(products_path, subbands, data.path)
        # Figures need one representative set: the band closest to z_obs.
        spectra = subbands.bands[
            int(np.argmin(np.abs(geometry.z_effective - data.z_obs)))
        ]
        log(f"  Done in {time.time() - start:.1f} s → {products_path}", quiet)
        return spectra, subbands, True

    log("  Computing 2D cylindrical power spectra …", quiet)
    spectra = analysis.compute_all_power_spectra(
        brightness_temp_field=data.brightness_temp_field,
        galaxy_overdensity=data.galaxy_overdensity,
        box_len_perp=data.BOX_LEN,
        box_len_los=data.L_los,
        n_bins_perp=n_bins_perp,
        n_bins_parallel=n_bins_parallel,
    )

    save_power_spectra(products_path, spectra, data.path)
    log(f"  Done in {time.time() - start:.1f} s → {products_path}", quiet)
    return spectra, None, True


def observational_stage(
    data: SimulationData,
    spectra,
    photoz_uncertainty: Optional[float] = None,
    wedge_buffer: Optional[float] = None,
    integration_time: Optional[float] = None,
    bandwidth: Optional[float] = None,
    noise_model: str = "scaling",
    mode_weighted: bool = False,
    z_obs: Optional[float] = None,
    quiet: bool = False,
) -> analysis.UncertaintyBudget:
    """
    Compute the uncertainty budget: photo-z damping, wedge, noise, and SNR.

    Every value defaults to the corresponding root attribute of the stored
    simulation; the four overrides let a survey or instrument parameter be
    swept from the command line without re-running 21cmFAST, since none of
    them affects the simulated fields.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation (supplies survey, instrument, and cosmology
        metadata).
    spectra : PowerSpectra
        Computed power spectra.
    photoz_uncertainty : float, optional
        Absolute σ_z override.  Defaults to the ``photoz_uncertainty``
        attribute, or 0.45 — the Euclid requirement σ_z/(1+z) < 0.05 at
        z = 7.  This is **not** the fractional value; see
        :func:`src.analysis.radial_smearing_length`.
    wedge_buffer : float, optional
        Wedge margin override [Mpc^-1].  Defaults to the ``wedge_buffer``
        attribute, or 0.0677 = 0.1 h Mpc^-1 (Pober et al. 2014 "moderate").
    integration_time : float, optional
        Integration-time override [s].  Defaults to the attribute, or 1000 h.
    bandwidth : float, optional
        Bandwidth override [Hz].  Defaults to the attribute, or 8 MHz.
    noise_model : {'scaling', 'physical'}, optional
        21 cm thermal-noise model.  ``'physical'`` uses Parsons (2017) Eq. 12
        / La Plante Eq. 11 resolved through the HERA baseline distribution;
        ``'scaling'`` (default) keeps the flat historical estimate.
    mode_weighted : bool, optional
        Apply the La Plante Eq. 19 ``sqrt(N_patch dN)`` weighting when summing
        bins.  Default ``False``, preserving the historical total.
    z_obs : float, optional
        Reference redshift override.  Defaults to the simulation's own
        ``z_obs``; the sub-band estimator passes each band's effective
        redshift instead, since the wedge slope, ``T_sys`` and the photo-z
        smearing all depend on it.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    UncertaintyBudget
        Every term of the budget, from :func:`src.analysis.compute_uncertainty_budget`.
    """
    def resolve(override: Optional[float], name: str, fallback: float) -> float:
        """Command-line override, else the stored attribute, else the default."""
        return float(override) if override is not None else data.get(name, fallback)

    budget = analysis.compute_uncertainty_budget(
        spectra=spectra,
        z_obs=data.z_obs if z_obs is None else float(z_obs),
        photoz_uncertainty=resolve(
            photoz_uncertainty, "photoz_uncertainty", 0.45
        ),
        wedge_buffer=resolve(wedge_buffer, "wedge_buffer", 0.0677),
        integration_time=resolve(
            integration_time, "integration_time", 1000 * 3600
        ),
        bandwidth=resolve(bandwidth, "bandwidth", 8e6),
        mean_galaxy_density=data.get("mean_galaxy_density", 3e-3),
        dish_diameter=data.get("HERA_DISH_DIAMETER", 14.0),
        f_21_hz=data.get("F_21_HZ", 1420.405e6),
        speed_of_light_mps=data.get("SPEED_OF_LIGHT_MPS", 3e8),
        hubble_constant=data.get("HUBBLE_CONSTANT", 67.36),
        omega_m=data.get("OMEGA_M_0", 0.315),
        speed_of_light_kms=data.get("SPEED_OF_LIGHT_KMS", 3e5),
        noise_model=noise_model,
        mode_weighted=mode_weighted,
    )

    log(f"  T_sys at {budget.observed_frequency_hz / 1e6:.2f} MHz "
        f": {budget.system_temperature_mK / 1e3:.1f} K", quiet)
    log(f"  Photo-z smearing    : σ_z = {budget.photoz_uncertainty:g}  →  "
        f"σ_r = {budget.radial_smearing:.1f} Mpc", quiet)
    log(f"  Photo-z kernel      : W = {budget.photoz_kernel.ravel()[0]:.3g} "
        f"at the first k_∥ bin", quiet)
    log(f"  Horizon wedge slope : {budget.horizon_slope:.3f}  "
        f"(buffer {budget.wedge_buffer:g} Mpc⁻¹)", quiet)
    log(f"  HERA FoV slope      : {budget.fov_slope:.3f}", quiet)
    log(f"  Modes outside wedge : "
        f"{budget.outside_wedge.sum()}/{budget.outside_wedge.size} "
        f"({budget.fraction_outside_wedge:.1%})", quiet)
    # P_noise_21cm is k_perp-resolved under --noise-model physical, so report
    # its finite range rather than assuming a scalar.
    noise_21cm = np.asarray(budget.snr.P_noise_21cm, dtype=float)
    if noise_21cm.ndim == 0:
        noise_text = f"{float(noise_21cm):.4g} mK² Mpc³"
    else:
        finite = noise_21cm[np.isfinite(noise_21cm)]
        unmeasurable = int((~np.isfinite(noise_21cm)).sum())
        noise_text = (
            f"{finite.min():.4g}–{finite.max():.4g} mK² Mpc³ over k⊥"
            + (f" ({unmeasurable} bins unsampled)" if unmeasurable else "")
        )
    log(f"  P_N,21 / P_N,gal    : {noise_text} / "
        f"{budget.snr.P_noise_galaxy:.4g} Mpc³", quiet)
    log(f"  Noise model         : {budget.noise_model}"
        f"   mode-weighted: {budget.mode_weighted}", quiet)
    log(f"  σ² from cosmic var. : "
        f"{budget.cosmic_variance_fraction:.1%}", quiet)
    log(f"  Total SNR (outside wedge) : {budget.total_snr:.3g} σ", quiet)

    return budget


def subband_observational_stage(
    data: SimulationData,
    subbands,
    photoz_uncertainty: Optional[float] = None,
    wedge_buffer: Optional[float] = None,
    integration_time: Optional[float] = None,
    noise_model: str = "scaling",
    mode_weighted: bool = False,
    quiet: bool = False,
):
    """
    Uncertainty budget per frequency sub-band, combined in quadrature.

    One call to :func:`observational_stage` per band, each at that band's own
    effective redshift and its own bandwidth. The redshift matters because
    the wedge slope, ``T_sys`` and the photo-z smearing all depend on it; the
    bandwidth matters because a spectrum measured over one band and a noise
    level computed for another are not commensurate (`TODO.md` P0.4).

    Parameters
    ----------
    data : SimulationData
        Loaded simulation, for the survey and instrument metadata.
    subbands : SubbandPowerSpectra
        Per-band spectra from :func:`power_spectra_stage`.
    photoz_uncertainty, wedge_buffer, integration_time : float, optional
        Overrides passed through to :func:`observational_stage`.  ``bandwidth``
        is deliberately absent: it is set per band by the geometry.
    noise_model : {'scaling', 'physical'}, optional
        21 cm thermal-noise model.
    mode_weighted : bool, optional
        Apply the La Plante Eq. 19 weighting within each band.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    UncertaintyBudget
        The band closest to ``z_obs``, for figures and the summary's scalar
        fields.
    list of UncertaintyBudget
        One budget per band, in lightcone order.
    float
        Combined significance, ``sqrt(sum(SNR_band^2))``.
    """
    band_budgets = []
    for index, band in enumerate(subbands.bands):
        z_eff = float(subbands.z_effective[index])
        band_bandwidth = float(subbands.bandwidth_hz[index])
        log(f"\n  Band {index}: z_eff = {z_eff:.3f}, "
            f"B = {band_bandwidth / 1e6:.2f} MHz", quiet)

        band_budgets.append(
            observational_stage(
                data, band,
                z_obs=z_eff,
                photoz_uncertainty=photoz_uncertainty,
                wedge_buffer=wedge_buffer,
                integration_time=integration_time,
                bandwidth=band_bandwidth,
                noise_model=noise_model,
                mode_weighted=mode_weighted,
                quiet=quiet,
            )
        )

    band_totals = np.array([b.total_snr for b in band_budgets], dtype=float)
    combined = analysis.combine_band_snr(band_totals)

    representative = int(np.argmin(np.abs(subbands.z_effective - data.z_obs)))
    log(f"\n  Per-band SNR        : "
        + ", ".join(f"{value:.3g}" for value in band_totals), quiet)
    log(f"  Combined SNR        : {combined:.3g} σ "
        f"(quadrature over {len(band_totals)} bands)", quiet)

    return band_budgets[representative], band_budgets, combined


def bias_stage(
    data: SimulationData,
    m_uv_bright: float = -22.0,
    quiet: bool = False,
):
    """
    Euclid ``M_UV`` selection of the halo catalogue and its effective bias.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a halo catalogue.
    m_uv_bright : float, optional
        Bright-end magnitude cut.  The faint end comes from the HDF5
        ``M_UV_limit`` attribute.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    EuclidSelection or None
        The selected sample, or None if there is no catalogue.
    BiasEstimate or None
        The bias estimate, or None if ``hmf`` is unavailable or no halo
        passed the cut.
    """
    if not data.has_halo_catalog:
        log("  No halo catalogue in the simulation output — skipping.", quiet)
        return None, None

    m_uv_faint = data.get("M_UV_limit", -18.0)

    selection = analysis.select_euclid_halos(
        sfr=data.sfr,
        halo_masses=data.halo_masses,
        M_UV_faint=m_uv_faint,
        M_UV_bright=m_uv_bright,
    )

    log(f"  Euclid window       : {m_uv_bright:.1f} < M_UV < {m_uv_faint:.1f}", quiet)
    log(f"  SFR window          : {selection.SFR_min:.3e} – "
        f"{selection.SFR_max:.3e} M_sun/yr", quiet)
    log(f"  Halos with SFR > 0  : {selection.n_valid:,}", quiet)
    log(f"  Selected halos      : {selection.n_selected:,}", quiet)

    if selection.n_selected == 0:
        log("  No halos passed the magnitude cut — skipping bias estimate.", quiet)
        return selection, None

    try:
        bias = analysis.effective_galaxy_bias(
            selection=selection,
            z_obs=data.z_obs,
            hubble_constant=data.get("HUBBLE_CONSTANT", 67.36),
        )
    except ImportError:
        log("  hmf not installed — skipping the effective bias estimate.", quiet)
        return selection, None

    log(f"  Mean galaxy bias    : <b_g> = {bias.mean_bias:.3f}  "
        f"(range {bias.bias_min:.3f} – {bias.bias_max:.3f})", quiet)
    return selection, bias


# ===========================================================================
#  Stage 3 — figures
# ===========================================================================

def figure_stage(
    groups: Sequence[str],
    data: SimulationData,
    spectra,
    budget: analysis.UncertaintyBudget,
    bias,
    output_dir: str,
    fmt: str = "png",
    quiet: bool = False,
    m_uv_bright: float = -22.0,
    galaxy_weighting: str = "number",
) -> List[str]:
    """
    Render and save the requested figure groups.

    Parameters
    ----------
    groups : sequence of str
        Any of :data:`PLOT_GROUPS`.
    data : SimulationData
        Loaded simulation.
    spectra : PowerSpectra
        Computed power spectra.
    budget : UncertaintyBudget
        Output of :func:`observational_stage`.
    bias : BiasEstimate or None
        Output of :func:`bias_stage`.
    output_dir : str
        Directory for the figure files.
    fmt : str, optional
        File format (``png``, ``pdf``, ``svg``).
    quiet : bool, optional
        Suppress progress output.
    m_uv_bright : float, optional
        Bright-end Euclid magnitude cut, passed to the selection-map figure so
        it selects identically to the bias stage.
    galaxy_weighting : {'number', 'luminosity'}, optional
        Per-halo weight used to rebuild the post-cut galaxy overdensity for
        the ``euclid`` figure group.

    Returns
    -------
    list of str
        Paths of the written figures.
    """
    written: List[str] = []

    def emit(name: str, fig) -> None:
        path = figures.save_figure(fig, output_dir, name, fmt=fmt)
        written.append(path)
        log(f"  wrote {path}", quiet)

    has_catalog = data.has_halo_catalog

    if "fields" in groups:
        emit("lightcone_fields", figures.plot_lightcone_fields(data))
        emit("lightcone_slice", figures.plot_lightcone_slice(data))

    if "halos" in groups:
        if has_catalog:
            emit("halo_catalogue", figures.plot_halo_catalogue(data))
            emit("sfr_relations", figures.plot_sfr_relations(data))
        else:
            log("  (no halo catalogue — skipping halo figures)", quiet)

    if "scaling" in groups:
        if has_catalog:
            emit("uv_luminosity_function", figures.plot_uv_luminosity_function(data))
            emit("stellar_mass_muv", figures.plot_stellar_mass_muv(data))
            emit("main_sequence", figures.plot_main_sequence(data))
            emit("uv_selection_maps", figures.plot_uv_selection_maps(
                data, M_UV_bright=m_uv_bright,
            ))
        else:
            log("  (no halo catalogue — skipping scaling-relation figures)", quiet)

    if "euclid" in groups:
        if has_catalog:
            # The post-cut delta_gal is expensive to deposit for a large
            # catalogue, so build it once and hand it to both figures.
            delta_gal, selection = figures.selected_galaxy_overdensity(
                data, M_UV_bright=m_uv_bright, weighting=galaxy_weighting,
            )
            emit("euclid_selected_catalogue", figures.plot_euclid_selected_catalogue(
                data, M_UV_bright=m_uv_bright,
            ))
            emit("selected_galaxy_overdensity", figures.plot_selected_galaxy_overdensity(
                data, weighting=galaxy_weighting,
                delta_gal=delta_gal, selection=selection,
            ))
            emit("galaxy_overdensity_on_21cm", figures.plot_galaxy_overdensity_on_21cm(
                data, weighting=galaxy_weighting,
                delta_gal=delta_gal, selection=selection,
            ))
        else:
            log("  (no halo catalogue — skipping Euclid-cut figures)", quiet)

    if "power" in groups:
        emit("power_spectra_2d", figures.plot_power_spectra(
            spectra, data,
            horizon_slope=budget.horizon_slope,
            fov_slope=budget.fov_slope,
        ))
        emit("galaxy_wedge", figures.plot_galaxy_wedge(
            spectra, data,
            horizon_slope=budget.horizon_slope,
            fov_slope=budget.fov_slope,
        ))
        # Real-space counterpart.  buffer=0 draws the bare horizon boundary,
        # matching the line plot_galaxy_wedge overlays; the budget's own
        # excision adds wedge_buffer on top of it.
        emit("wedge_real_space", figures.plot_wedge_real_space(
            data, horizon_slope=budget.horizon_slope, wedge_buffer=0.0,
        ))

    if "snr" in groups:
        emit("cross_snr", figures.plot_snr(
            spectra, budget.snr, budget.P_cross_observed, data,
            horizon_slope=budget.horizon_slope,
            fov_slope=budget.fov_slope,
        ))

    if "budget" in groups:
        emit("uncertainty_budget", figures.plot_uncertainty_budget(budget, data))
        emit("photoz_suppression", figures.plot_photoz_suppression(budget, data))

    if "bias" in groups:
        if bias is not None:
            emit("galaxy_bias", figures.plot_bias_diagnostic(bias, data.z_obs))
        else:
            log("  (no bias estimate — skipping bias figure)", quiet)

    return written


# ===========================================================================
#  Stage 4 — summary
# ===========================================================================

def build_summary(
    data: SimulationData,
    spectra,
    budget: analysis.UncertaintyBudget,
    selection,
    bias,
    figure_paths: Sequence[str],
    ran_simulation: bool,
    recomputed_spectra: bool,
    estimator: str = "coeval",
    subbands=None,
    band_budgets=None,
    combined_snr: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Collect the pipeline's scalar results into a JSON-serialisable dict.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation.
    spectra : PowerSpectra
        Computed power spectra.
    budget : UncertaintyBudget
        Output of :func:`observational_stage`.
    selection : EuclidSelection or None
        Output of :func:`bias_stage`.
    bias : BiasEstimate or None
        Output of :func:`bias_stage`.
    figure_paths : sequence of str
        Figures written this run.
    ran_simulation : bool
        Whether the simulation stage executed.
    recomputed_spectra : bool
        Whether the power spectra were recomputed rather than loaded.
    estimator : {'coeval', 'lightcone'}, optional
        Which formalism produced the spectra.
    subbands : SubbandPowerSpectra, optional
        Per-band spectra, when the lightcone estimator ran.
    band_budgets : list of UncertaintyBudget, optional
        One budget per band, from :func:`subband_observational_stage`.
    combined_snr : float, optional
        Quadrature sum of the per-band totals.

    Returns
    -------
    dict
        Summary of the run.
    """
    large_scale_cross = float(np.nanmean(spectra.P_cross[:5, :5]))

    summary: Dict[str, Any] = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "data_file": os.path.abspath(data.path),
        "ran_simulation": ran_simulation,
        "recomputed_power_spectra": recomputed_spectra,
        "simulation": {
            "HII_DIM": data.HII_DIM,
            "BOX_LEN_Mpc": data.BOX_LEN,
            "N_z": data.N_z,
            "L_los_Mpc": data.L_los,
            "cell_size_Mpc": data.cell_size,
            "z_min": data.z_min,
            "z_max": data.z_max,
            "z_obs": data.z_obs,
            "mean_neutral_fraction": float(np.mean(data.neutral_fraction)),
            "mean_brightness_temp_mK": float(np.mean(data.brightness_temp_field)),
            "galaxy_bias_sim": data.get("galaxy_bias", np.nan),
            "beta_rsd": data.get("beta_rsd", np.nan),
            "n_halos_total": data.n_halos_total,
            "halo_sampling_factor": data.halo_sampling_factor,
        },
        "power_spectra": {
            "k_perp_min_Mpc-1": float(spectra.k_perp.min()),
            "k_perp_max_Mpc-1": float(spectra.k_perp.max()),
            "k_parallel_min_Mpc-1": float(spectra.k_parallel.min()),
            "k_parallel_max_Mpc-1": float(spectra.k_parallel.max()),
            "empty_bins": int(np.sum(spectra.mode_counts == 0)),
            "total_bins": int(spectra.mode_counts.size),
            "large_scale_cross_mean": large_scale_cross,
            "large_scale_anticorrelated": bool(large_scale_cross < 0),
        },
        # Every scalar term of the uncertainty budget, in the order the
        # calculation applies them: damping -> wedge -> noise -> variance.
        "uncertainty_budget": budget.as_dict(),
        "figures": [os.path.abspath(p) for p in figure_paths],
    }

    summary["estimator"] = estimator
    if subbands is not None:
        # Per-band results are the only ones with a well-defined redshift
        # under the lightcone estimator, so they are reported individually as
        # well as combined.
        summary["subbands"] = {
            "n_bands": subbands.n_bands,
            "z_effective": [float(z) for z in subbands.z_effective],
            "z_min": [float(z) for z in subbands.z_min],
            "z_max": [float(z) for z in subbands.z_max],
            "bandwidth_MHz": [float(b) / 1e6 for b in subbands.bandwidth_hz],
            "los_length_Mpc": [float(v) for v in subbands.los_length_mpc],
            "n_slices": [int(n) for n in subbands.n_slices],
            "total_snr_per_band": (
                [float(b.total_snr) for b in band_budgets]
                if band_budgets is not None else None
            ),
            "combined_total_snr": (
                float(combined_snr) if combined_snr is not None else None
            ),
            "representative_z_eff": float(
                subbands.z_effective[
                    int(np.argmin(np.abs(subbands.z_effective - data.z_obs)))
                ]
            ),
        }

    # Point back at the simulation run that produced these fields.  This
    # summary is overwritten every run; the manifest it names is not, so an
    # analysis-only run can still say which simulation its numbers came from.
    # Absent for HDF5 files written before run manifests existed.
    # h5py returns numpy scalars (np.int64, np.bool_), which json.dump cannot
    # serialise — it raises partway through and leaves a truncated, unparseable
    # file.  Coerce here rather than relying on a serialiser fallback.
    def _plain(value: Any) -> Any:
        """Convert an HDF5 attribute to a JSON-serialisable Python scalar."""
        if value is None or isinstance(value, str):
            return value
        item = getattr(value, "item", None)
        return item() if callable(item) else value

    summary["source_run"] = {
        key: _plain(data.attrs.get(key))
        for key in ("run_id", "run_manifest", "random_seed", "n_threads")
    }

    # Backwards-compatible alias: earlier summaries carried these five keys
    # under "observation", and downstream notes/scripts still read them.
    summary["observation"] = {
        "radial_smearing_Mpc": budget.radial_smearing,
        "horizon_wedge_slope": budget.horizon_slope,
        "fov_wedge_slope": budget.fov_slope,
        "modes_outside_wedge_fraction": budget.fraction_outside_wedge,
        "P_noise_21cm": budget.snr.P_noise_21cm,
        "P_noise_galaxy": budget.snr.P_noise_galaxy,
        "total_snr_sigma": budget.total_snr,
        "detection_above_5sigma": budget.detected,
    }

    if selection is not None:
        summary["euclid_selection"] = {
            "M_UV_bright": selection.M_UV_bright,
            "M_UV_faint": selection.M_UV_faint,
            "SFR_min_Msun_yr": selection.SFR_min,
            "SFR_max_Msun_yr": selection.SFR_max,
            "n_halos_sfr_positive": selection.n_valid,
            "n_halos_selected": selection.n_selected,
        }

    if bias is not None:
        summary["effective_galaxy_bias"] = {
            "mean_bias": bias.mean_bias,
            "bias_min": bias.bias_min,
            "bias_max": bias.bias_max,
            "n_selected": bias.n_selected,
        }

    return summary


def write_summary(summary: Dict[str, Any], path: str) -> None:
    """
    Write the run summary as JSON.

    Parameters
    ----------
    summary : dict
        Output of :func:`build_summary`.
    path : str
        Destination file.

    Notes
    -----
    ``default=str`` is a backstop, not the primary defence: values are coerced
    to plain Python scalars in :func:`build_summary`.  Non-finite floats are
    emitted by ``json`` as ``Infinity`` / ``NaN``, which is valid for Python's
    own ``json.load`` but not strict JSON — the ``physical`` noise model
    produces ``inf`` for unsampled ``k_perp`` bins, so consumers outside
    Python should expect them.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    # Written to a temporary file and renamed: json.dump writes incrementally,
    # so a mid-write failure would otherwise leave a truncated summary that
    # looks valid until something tries to parse it.
    temporary = f"{path}.tmp"
    with open(temporary, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    os.replace(temporary, path)


def print_report(summary: Dict[str, Any]) -> None:
    """
    Print the headline numbers of a completed run.

    Parameters
    ----------
    summary : dict
        Output of :func:`build_summary`.
    """
    sim = summary["simulation"]
    ub = summary["uncertainty_budget"]
    ps = summary["power_spectra"]

    print(f"\n{SEPARATOR}\n  PIPELINE SUMMARY\n{SEPARATOR}")
    print(f"  Box               : {sim['BOX_LEN_Mpc']:.0f} Mpc, "
          f"{sim['HII_DIM']}² × {sim['N_z']} cells")
    print(f"  Estimator         : {summary.get('estimator', 'coeval')}")
    print(f"  Lightcone         : z = {sim['z_min']} → {sim['z_max']}  "
          f"(z_obs = {sim['z_obs']})")
    print(f"  <x_HI>            : {sim['mean_neutral_fraction']:.3f}")
    print(f"  Cross-spectrum    : large-scale mean = "
          f"{ps['large_scale_cross_mean']:.3e}  "
          f"({'anti-correlated' if ps['large_scale_anticorrelated'] else 'positive'})")

    print(f"\n  Uncertainty budget")
    print(f"    Photo-z         : σ_z = {ub['photoz_uncertainty_sigma_z']:g}  →  "
          f"σ_r = {ub['radial_smearing_Mpc']:.1f} Mpc  "
          f"(W = {ub['photoz_kernel_first_bin']:.3g} at the first k_∥ bin)")
    print(f"    Wedge           : slope {ub['horizon_wedge_slope']:.3f}, "
          f"buffer {ub['wedge_buffer_Mpc-1']:g} Mpc⁻¹  →  "
          f"{ub['modes_outside_wedge']}/{ub['modes_total']} modes "
          f"({ub['modes_outside_wedge_fraction']:.1%}) usable")
    print(f"    Noise           : T_sys = {ub['system_temperature_K']:.1f} K, "
          f"P_N,21 = {ub['P_noise_21cm']:.4g} mK² Mpc³, "
          f"P_N,gal = {ub['P_noise_galaxy']:.4g} Mpc³")
    print(f"    Variance split  : "
          f"{ub['cosmic_variance_fraction']:.1%} cosmic variance, "
          f"{1.0 - ub['cosmic_variance_fraction']:.1%} noise coupling")
    print(f"    Total SNR       : {ub['total_snr_sigma']:.3g} σ  "
          f"({'detection' if ub['detection_above_5sigma'] else 'no detection'} at 5σ)")

    if "subbands" in summary:
        sub = summary["subbands"]
        print(f"\n  Sub-bands ({sub['n_bands']}, lightcone estimator)")
        print("    band   z_eff      B [MHz]   L_los [Mpc]   slices   SNR")
        for index in range(sub["n_bands"]):
            per_band = sub["total_snr_per_band"]
            snr_text = (
                f"{per_band[index]:.3g}" if per_band is not None else "--"
            )
            print(f"    {index:>4}   {sub['z_effective'][index]:.4f}   "
                  f"{sub['bandwidth_MHz'][index]:>7.2f}   "
                  f"{sub['los_length_Mpc'][index]:>11.1f}   "
                  f"{sub['n_slices'][index]:>6}   {snr_text}")
        if sub["combined_total_snr"] is not None:
            print(f"    Combined SNR  : {sub['combined_total_snr']:.3g} σ "
                  f"(quadrature over bands)")
        print(f"    Budget above is band {sub['representative_z_eff']:.4f}, "
              f"the one closest to z_obs")

    if "effective_galaxy_bias" in summary:
        bias = summary["effective_galaxy_bias"]
        print(f"  <b_g> (Euclid-selected) : {bias['mean_bias']:.3f}  "
              f"from {bias['n_selected']:,} halos")

    print(f"  Figures written   : {len(summary['figures'])}")


# ===========================================================================
#  CLI
# ===========================================================================

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse the command-line arguments.

    Parameters
    ----------
    argv : sequence of str, optional
        Argument list; defaults to ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        Parsed options.
    """
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description=(
            "Run the full 21 cm × galaxy cross-correlation pipeline: "
            "optional simulation, analysis (fresh or cached), and figures."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python run_pipeline.py                      analyse stored results\n"
            "  python run_pipeline.py --sim force          re-run the simulation first\n"
            "  python run_pipeline.py --analysis force     recompute the power spectra\n"
            "  python run_pipeline.py --plots power snr    only the k-space figures\n"
            "  python run_pipeline.py --max-halos 5000000  cap catalogue memory\n"
        ),
    )

    parser.add_argument(
        "--sim", choices=("auto", "force", "skip"), default="auto",
        help="run 21cmFAST: auto = only if the output is missing (default), "
             "force = always, skip = never",
    )
    parser.add_argument(
        "--analysis", choices=("auto", "force", "skip"), default="auto",
        help="power spectra: auto = recompute if the cache is stale (default), "
             "force = always recompute, skip = require the cache",
    )
    parser.add_argument(
        "--plots", nargs="+", default=["all"],
        choices=("all", "none") + PLOT_GROUPS,
        help="figure groups to render (default: all)",
    )
    parser.add_argument("--data", default=DEFAULT_DATA,
                        help=f"simulation HDF5 (default: {DEFAULT_DATA})")
    parser.add_argument("--products", default=DEFAULT_PRODUCTS,
                        help=f"analysis cache (default: {DEFAULT_PRODUCTS})")
    parser.add_argument("--figdir", default=DEFAULT_FIGDIR,
                        help=f"figure directory (default: {DEFAULT_FIGDIR})")
    parser.add_argument("--summary", default=DEFAULT_SUMMARY,
                        help=f"summary JSON (default: {DEFAULT_SUMMARY})")
    parser.add_argument("--sim-script", default="run_simulation.py",
                        help="simulation script (default: run_simulation.py)")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"),
                        help="figure file format (default: png)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="figure resolution (default: 200)")
    parser.add_argument(
        "--max-halos", type=int, default=0,
        help="cap on halos loaded from the catalogue, uniformly strided "
             "(default: 0 = all). Lower this if memory is tight.",
    )
    parser.add_argument("--m-uv-bright", type=float, default=-22.0,
                        help="bright-end Euclid magnitude cut (default: -22)")
    parser.add_argument(
        "--galaxy-weighting", choices=analysis.GALAXY_WEIGHTING_MODES,
        default="number",
        help="per-halo weight used to rebuild the post-cut galaxy overdensity "
             "for the 'euclid' figure group (default: number)",
    )

    parser.add_argument(
        "--noise-model", choices=("scaling", "physical"), default="scaling",
        help="21 cm thermal noise: 'scaling' (default, flat, historical) or "
             "'physical' (Parsons 2017 Eq. 12 / La Plante Eq. 11, resolved in "
             "k_perp through the HERA baseline distribution, ~10^3 larger)",
    )
    parser.add_argument(
        "--mode-weighted", action="store_true",
        help="apply the La Plante Eq. 19 sqrt(N_patch dN) mode weighting when "
             "summing bins (default off; raises the total SNR ~10x)",
    )

    budget_group = parser.add_argument_group(
        "uncertainty budget",
        "Survey and instrument overrides. None of these affects the simulated "
        "fields, so they can be swept without --sim force. Each defaults to "
        "the corresponding attribute of the stored HDF5.",
    )
    budget_group.add_argument(
        "--sigma-z", type=float, default=None, metavar="SIGMA_Z",
        help="absolute photo-z uncertainty sigma_z, NOT sigma_z/(1+z) "
             "(stored default: 0.45, the Euclid requirement at z = 7)",
    )
    budget_group.add_argument(
        "--wedge-buffer", type=float, default=None, metavar="MPC-1",
        help="foreground-wedge margin [Mpc^-1] (stored default: 0.0677 "
             "= 0.1 h/Mpc, Pober et al. 2014 'moderate')",
    )
    budget_group.add_argument(
        "--integration-time", type=float, default=None, metavar="SECONDS",
        help="HERA integration time [s] (stored default: 3.6e6 = 1000 h)",
    )
    budget_group.add_argument(
        "--bandwidth", type=float, default=None, metavar="HZ",
        help="per-band bandwidth [Hz] (stored default: 8e6)",
    )
    estimator_group = parser.add_argument_group(
        "estimator (TODO.md P0)",
        "Which power-spectrum formalism to use. 'auto' follows the stored "
        "simulation's own `estimator` attribute, so the two halves of the "
        "pipeline cannot silently disagree.",
    )
    estimator_group.add_argument(
        "--estimator", default="auto", choices=("auto", "coeval", "lightcone"),
        help="coeval: one FFT over the box, global mean, no taper (default "
             "for data written before P0). lightcone: per-slice means, "
             "Blackman-Harris taper, one spectrum per sub-band at its own "
             "effective redshift and bandwidth.",
    )
    estimator_group.add_argument(
        "--subband-bandwidth", type=float, default=None, metavar="HZ",
        help="target per-band bandwidth for --estimator lightcone [Hz]; "
             "defaults to the noise bandwidth, which is the point of P0.4",
    )

    smoke_group = parser.add_argument_group(
        "smoke test",
        "Pre-flight check: run every stage on a tiny configuration and assert "
        "the outputs have the right shapes. Not a science run.",
    )
    smoke_group.add_argument(
        "--smoke-test", action="store_true",
        help="run the whole pipeline on the reduced configuration in "
             "src/smoke_test.py, verify every stage's output shape, and print "
             "a stage-by-stage report. Redirects --data/--products/--figdir/"
             "--summary to outputs/smoke_test/ unless they are given "
             "explicitly, so a real run's products are never overwritten.",
    )

    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress output")

    return parser.parse_args(argv)


def resolve_plot_groups(requested: Sequence[str]) -> List[str]:
    """
    Expand the ``--plots`` selection into concrete figure groups.

    Parameters
    ----------
    requested : sequence of str
        Raw ``--plots`` values, possibly containing ``all`` or ``none``.

    Returns
    -------
    list of str
        Groups to render; empty when ``none`` was requested.
    """
    if "none" in requested:
        return []
    if "all" in requested:
        return list(PLOT_GROUPS)
    return [group for group in PLOT_GROUPS if group in requested]


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run the pipeline end to end.

    Parameters
    ----------
    argv : sequence of str, optional
        Argument list; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code: 0 on success, 1 on a handled failure.
    """
    args = parse_args(argv)
    quiet = args.quiet
    started = time.time()

    # ── Smoke test ────────────────────────────────────────────────────────
    # Redirect every output path into outputs/smoke_test/ unless the caller
    # named one explicitly, so a pre-flight check cannot overwrite a real
    # run's products.  Nothing here touches the production defaults.
    smoke_report = None
    if args.smoke_test:
        from src.smoke_test import (
            SMOKE_OUTPUT_DIR, SMOKE_TEST_OVERRIDES, SmokeReport,
            check_figures, check_mcmc_chain, check_power_spectra,
            check_simulation_output, check_summary, check_uncertainty_budget,
            describe_overrides,
        )

        smoke_report = SmokeReport()
        defaults = {
            "data": DEFAULT_DATA, "products": DEFAULT_PRODUCTS,
            "figdir": DEFAULT_FIGDIR, "summary": DEFAULT_SUMMARY,
        }
        for name, default in defaults.items():
            if getattr(args, name) == default:
                setattr(args, name, os.path.join(
                    SMOKE_OUTPUT_DIR, os.path.basename(default)
                ))
        if args.max_halos == 0:
            args.max_halos = SMOKE_TEST_OVERRIDES["max_halos"]["value"]

        if not quiet:
            print(describe_overrides())
            print(f"  Outputs → {SMOKE_OUTPUT_DIR}/\n")

    plot_groups = resolve_plot_groups(args.plots)
    figures.apply_plot_style(dpi=args.dpi)

    # Catalogue-dependent work: halo/scaling/bias figures and the bias stage.
    needs_catalog = bool({"halos", "scaling", "euclid", "bias"} & set(plot_groups))
    if args.smoke_test:
        # A pre-flight check must exercise the halo catalogue even with
        # --plots none, since that is one of the stages it exists to verify.
        needs_catalog = True

    try:
        # ── Stage 1 ───────────────────────────────────────────────────────
        ran_simulation = run_simulation_stage(
            mode=args.sim,
            data_path=args.data,
            script=args.sim_script,
            smoke_test=args.smoke_test,
            quiet=quiet,
        )

        # ── Load ──────────────────────────────────────────────────────────
        banner("STAGE 2 — ANALYSIS", quiet)
        log(f"  Loading {args.data} …", quiet)

        data = load_simulation(
            args.data,
            max_halos=args.max_halos,
            load_halos=needs_catalog,
        )

        log(f"  Fields   : {data.brightness_temp_field.shape}", quiet)
        log(f"  Halos    : {data.n_halos_total:,} in file, "
            f"{data.halo_masses.shape[0]:,} loaded "
            f"(1/{data.halo_sampling_factor:.0f} sampling)", quiet)

        # ── Stage 2 ───────────────────────────────────────────────────────
        # The estimator follows the stored simulation unless overridden: a
        # lightcone sampled uniformly in redshift must not be analysed as if
        # its slices were evenly spaced in comoving distance.
        estimator = args.estimator
        if estimator == "auto":
            estimator = data.get_str("estimator", "coeval")
        log(f"  Estimator: {estimator}"
            + (" (from the stored simulation)" if args.estimator == "auto"
               else " (--estimator)"), quiet)

        noise_bandwidth = (
            float(args.bandwidth) if args.bandwidth is not None
            else data.get("bandwidth", 8e6)
        )
        subband_bandwidth = (
            float(args.subband_bandwidth)
            if args.subband_bandwidth is not None else noise_bandwidth
        )

        spectra, subbands, recomputed = power_spectra_stage(
            mode=args.analysis,
            data=data,
            products_path=args.products,
            estimator=estimator,
            subband_bandwidth=subband_bandwidth,
            quiet=quiet,
        )

        band_budgets, combined_snr = None, None
        if subbands is not None:
            budget, band_budgets, combined_snr = subband_observational_stage(
                data, subbands,
                photoz_uncertainty=args.sigma_z,
                wedge_buffer=args.wedge_buffer,
                integration_time=args.integration_time,
                noise_model=args.noise_model,
                mode_weighted=args.mode_weighted,
                quiet=quiet,
            )
        else:
            budget = observational_stage(
                data, spectra,
                photoz_uncertainty=args.sigma_z,
                wedge_buffer=args.wedge_buffer,
                integration_time=args.integration_time,
                bandwidth=args.bandwidth,
                noise_model=args.noise_model,
                mode_weighted=args.mode_weighted,
                quiet=quiet,
            )
        save_uncertainty_budget(args.products, budget)
        log(f"  Uncertainty budget cached → {args.products}", quiet)

        selection, bias = (None, None)
        if needs_catalog:
            selection, bias = bias_stage(
                data, m_uv_bright=args.m_uv_bright, quiet=quiet,
            )

        # ── Stage 3 ───────────────────────────────────────────────────────
        figure_paths: List[str] = []
        if plot_groups:
            banner("STAGE 3 — FIGURES", quiet)
            figure_paths = figure_stage(
                groups=plot_groups,
                data=data,
                spectra=spectra,
                budget=budget,
                bias=bias,
                output_dir=args.figdir,
                fmt=args.format,
                quiet=quiet,
                m_uv_bright=args.m_uv_bright,
                galaxy_weighting=args.galaxy_weighting,
            )
        else:
            log("\n  Figures skipped (--plots none).", quiet)

        # ── Stage 4 ───────────────────────────────────────────────────────
        summary = build_summary(
            data=data,
            spectra=spectra,
            budget=budget,
            selection=selection,
            bias=bias,
            figure_paths=figure_paths,
            ran_simulation=ran_simulation,
            recomputed_spectra=recomputed,
            estimator=estimator,
            subbands=subbands,
            band_budgets=band_budgets,
            combined_snr=combined_snr,
        )
        write_summary(summary, args.summary)

        if not quiet:
            print_report(summary)
            print(f"  Summary written   : {args.summary}")
            print(f"  Elapsed           : {time.time() - started:.1f} s")

        # ── Smoke-test verification ───────────────────────────────────────
        # Shapes, dtypes and finiteness at every stage — an explicit check
        # rather than "the run did not crash".
        if smoke_report is not None:
            n_perp = int(data.get("n_bins_perp", 20))
            n_parallel = int(data.get("n_bins_parallel", 20))

            check_simulation_output(
                data,
                expected_hii_dim=int(data.HII_DIM),
                expected_n_z=int(data.N_z),
                report=smoke_report,
            )
            check_power_spectra(
                spectra, n_perp, n_parallel,
                report=smoke_report, subbands=subbands,
            )
            check_uncertainty_budget(
                budget, n_perp, n_parallel, report=smoke_report,
            )
            check_figures(figure_paths, report=smoke_report)
            check_summary(summary, report=smoke_report)
            check_mcmc_chain(report=smoke_report)

            print("\n" + smoke_report.render())
            if not smoke_report.passed:
                return 1

    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        if smoke_report is not None:
            smoke_report.add("pipeline", False, f"aborted: {exc}")
            print("\n" + smoke_report.render())
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
