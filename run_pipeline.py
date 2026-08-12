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
    load_simulation,
    products_are_stale,
    save_power_spectra,
)

warnings.filterwarnings("ignore")

# ── Figure groups ─────────────────────────────────────────────────────────
#   fields   lightcone field slices (Part 2)
#   halos    halo catalogue positions, masses, SFR scaling relations (Part 2)
#   scaling  UVLF, M_star–M_UV, star-forming main sequence (Part 2)
#   power    2D cylindrical power spectra (Part 3)
#   snr      per-mode SNR and photo-z damped cross-power (Part 3)
#   bias     Euclid-selected halo masses and b_h(M, z) (Part 3)
PLOT_GROUPS = ("fields", "halos", "scaling", "power", "snr", "bias")

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

    result = subprocess.run([sys.executable, script], cwd=REPO_ROOT)

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
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    PowerSpectra
        The three spectra on a shared ``(k_perp, k_parallel)`` grid.
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
        spectra, attrs = load_power_spectra(products_path)
        log(f"  Loaded cached power spectra: {products_path}", quiet)
        if stale:
            log("  WARNING: cache is older than the simulation output "
                "(use --analysis force to refresh).", quiet)
        return spectra, False

    log("  Computing 2D cylindrical power spectra …", quiet)
    start = time.time()

    spectra = analysis.compute_all_power_spectra(
        brightness_temp_field=data.brightness_temp_field,
        galaxy_overdensity=data.galaxy_overdensity,
        box_len_perp=data.BOX_LEN,
        box_len_los=data.L_los,
        n_bins_perp=int(data.get("n_bins_perp", 20)),
        n_bins_parallel=int(data.get("n_bins_parallel", 20)),
    )

    save_power_spectra(products_path, spectra, data.path)
    log(f"  Done in {time.time() - start:.1f} s → {products_path}", quiet)
    return spectra, True


def observational_stage(
    data: SimulationData,
    spectra,
    quiet: bool = False,
) -> Dict[str, Any]:
    """
    Apply photo-z damping, wedge excision, and the noise/SNR calculation.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation (supplies survey and instrument metadata).
    spectra : PowerSpectra
        Computed power spectra.
    quiet : bool, optional
        Suppress progress output.

    Returns
    -------
    dict
        Keys: ``P_cross_observed``, ``P_galaxy_observed``, ``outside_wedge``,
        ``snr`` (:class:`src.analysis.SNRResult`), ``horizon_slope``,
        ``fov_slope``, ``radial_smearing``.
    """
    z_obs = data.z_obs
    h0 = data.get("HUBBLE_CONSTANT", 67.36)
    omega_m = data.get("OMEGA_M_0", 0.315)
    c_kms = data.get("SPEED_OF_LIGHT_KMS", 3e5)

    cosmology = dict(hubble_constant=h0, omega_m=omega_m, speed_of_light_kms=c_kms)

    # ── Photo-z damping ───────────────────────────────────────────────────
    radial_smearing = analysis.radial_smearing_length(
        # Absolute sigma_z (not sigma_z/(1+z)); 0.45 at z = 7 matches the
        # Euclid requirement sigma_z/(1+z) < 0.05.
        photoz_uncertainty=data.get("photoz_uncertainty", 0.45),
        z_obs=z_obs, **cosmology,
    )
    kernel = analysis.photoz_damping_kernel(spectra.k_parallel, radial_smearing)

    p_cross_observed = spectra.P_cross * kernel
    p_galaxy_observed = spectra.P_galaxy_auto * kernel ** 2

    # ── Foreground wedge ──────────────────────────────────────────────────
    horizon_slope = analysis.horizon_wedge_slope(z_obs, **cosmology)
    fov_slope = analysis.fov_wedge_slope(
        z_obs,
        dish_diameter=data.get("HERA_DISH_DIAMETER", 14.0),
        f_21_hz=data.get("F_21_HZ", 1420.405e6),
        speed_of_light_mps=data.get("SPEED_OF_LIGHT_MPS", 3e8),
        **cosmology,
    )
    outside_wedge = analysis.foreground_wedge_mask(
        spectra.k_perp, spectra.k_parallel,
        slope=horizon_slope,
        # 0.1 h Mpc^-1 (Pober et al. 2014 "moderate"; 21cmSense default)
        buffer=data.get("wedge_buffer", 0.0677),
    )

    # ── Noise and SNR ─────────────────────────────────────────────────────
    p_noise_21cm = analysis.hera_thermal_noise_power(
        z_obs=z_obs,
        integration_time=data.get("integration_time", 1000 * 3600),
        bandwidth=data.get("bandwidth", 8e6),
        f_21_hz=data.get("F_21_HZ", 1420.405e6),
    )
    p_noise_galaxy = 1.0 / data.get("mean_galaxy_density", 3e-3)

    snr = analysis.cross_power_snr(
        P_cross_observed=p_cross_observed,
        P_21cm_auto=spectra.P_21cm_auto,
        P_galaxy_observed=p_galaxy_observed,
        P_noise_21cm=p_noise_21cm,
        P_noise_galaxy=p_noise_galaxy,
        outside_wedge=outside_wedge,
    )

    log(f"  Photo-z smearing    : σ_r = {radial_smearing:.1f} Mpc", quiet)
    log(f"  Horizon wedge slope : {horizon_slope:.3f}", quiet)
    log(f"  HERA FoV slope      : {fov_slope:.3f}", quiet)
    log(f"  Modes outside wedge : "
        f"{outside_wedge.sum()}/{outside_wedge.size} "
        f"({outside_wedge.mean():.1%})", quiet)
    log(f"  Total SNR (outside wedge) : {snr.total_snr:.1f} σ", quiet)

    return {
        "P_cross_observed": p_cross_observed,
        "P_galaxy_observed": p_galaxy_observed,
        "outside_wedge": outside_wedge,
        "snr": snr,
        "horizon_slope": horizon_slope,
        "fov_slope": fov_slope,
        "radial_smearing": radial_smearing,
    }


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
    observed: Dict[str, Any],
    bias,
    output_dir: str,
    fmt: str = "png",
    quiet: bool = False,
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
    observed : dict
        Output of :func:`observational_stage`.
    bias : BiasEstimate or None
        Output of :func:`bias_stage`.
    output_dir : str
        Directory for the figure files.
    fmt : str, optional
        File format (``png``, ``pdf``, ``svg``).
    quiet : bool, optional
        Suppress progress output.

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
        else:
            log("  (no halo catalogue — skipping scaling-relation figures)", quiet)

    if "power" in groups:
        emit("power_spectra_2d", figures.plot_power_spectra(
            spectra, data,
            horizon_slope=observed["horizon_slope"],
            fov_slope=observed["fov_slope"],
        ))

    if "snr" in groups:
        emit("cross_snr", figures.plot_snr(
            spectra, observed["snr"], observed["P_cross_observed"], data,
            horizon_slope=observed["horizon_slope"],
            fov_slope=observed["fov_slope"],
        ))

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
    observed: Dict[str, Any],
    selection,
    bias,
    figure_paths: Sequence[str],
    ran_simulation: bool,
    recomputed_spectra: bool,
) -> Dict[str, Any]:
    """
    Collect the pipeline's scalar results into a JSON-serialisable dict.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation.
    spectra : PowerSpectra
        Computed power spectra.
    observed : dict
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

    Returns
    -------
    dict
        Summary of the run.
    """
    snr = observed["snr"]
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
        "observation": {
            "radial_smearing_Mpc": observed["radial_smearing"],
            "horizon_wedge_slope": observed["horizon_slope"],
            "fov_wedge_slope": observed["fov_slope"],
            "modes_outside_wedge_fraction": float(observed["outside_wedge"].mean()),
            "P_noise_21cm": snr.P_noise_21cm,
            "P_noise_galaxy": snr.P_noise_galaxy,
            "total_snr_sigma": snr.total_snr,
            "detection_above_5sigma": bool(snr.total_snr > 5.0),
        },
        "figures": [os.path.abspath(p) for p in figure_paths],
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
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)


def print_report(summary: Dict[str, Any]) -> None:
    """
    Print the headline numbers of a completed run.

    Parameters
    ----------
    summary : dict
        Output of :func:`build_summary`.
    """
    sim = summary["simulation"]
    obs = summary["observation"]
    ps = summary["power_spectra"]

    print(f"\n{SEPARATOR}\n  PIPELINE SUMMARY\n{SEPARATOR}")
    print(f"  Box               : {sim['BOX_LEN_Mpc']:.0f} Mpc, "
          f"{sim['HII_DIM']}² × {sim['N_z']} cells")
    print(f"  Lightcone         : z = {sim['z_min']} → {sim['z_max']}  "
          f"(z_obs = {sim['z_obs']})")
    print(f"  <x_HI>            : {sim['mean_neutral_fraction']:.3f}")
    print(f"  Cross-spectrum    : large-scale mean = "
          f"{ps['large_scale_cross_mean']:.3e}  "
          f"({'anti-correlated' if ps['large_scale_anticorrelated'] else 'positive'})")
    print(f"  Modes outside wedge : {obs['modes_outside_wedge_fraction']:.1%}")
    print(f"  Total SNR         : {obs['total_snr_sigma']:.1f} σ  "
          f"({'detection' if obs['detection_above_5sigma'] else 'no detection'} at 5σ)")

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

    plot_groups = resolve_plot_groups(args.plots)
    figures.apply_plot_style(dpi=args.dpi)

    # Catalogue-dependent work: halo/scaling/bias figures and the bias stage.
    needs_catalog = bool({"halos", "scaling", "bias"} & set(plot_groups))

    try:
        # ── Stage 1 ───────────────────────────────────────────────────────
        ran_simulation = run_simulation_stage(
            mode=args.sim,
            data_path=args.data,
            script=args.sim_script,
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
        spectra, recomputed = power_spectra_stage(
            mode=args.analysis,
            data=data,
            products_path=args.products,
            quiet=quiet,
        )
        observed = observational_stage(data, spectra, quiet=quiet)

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
                observed=observed,
                bias=bias,
                output_dir=args.figdir,
                fmt=args.format,
                quiet=quiet,
            )
        else:
            log("\n  Figures skipped (--plots none).", quiet)

        # ── Stage 4 ───────────────────────────────────────────────────────
        summary = build_summary(
            data=data,
            spectra=spectra,
            observed=observed,
            selection=selection,
            bias=bias,
            figure_paths=figure_paths,
            ran_simulation=ran_simulation,
            recomputed_spectra=recomputed,
        )
        write_summary(summary, args.summary)

        if not quiet:
            print_report(summary)
            print(f"  Summary written   : {args.summary}")
            print(f"  Elapsed           : {time.time() - started:.1f} s")

    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
