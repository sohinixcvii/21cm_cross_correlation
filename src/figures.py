#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
figures.py — Figure generation for the 21 cm × galaxy cross-correlation pipeline
=================================================================================

Every figure previously drawn inline in ``notebooks/plot_fields.ipynb``
(Part 2) and ``notebooks/analysis.ipynb`` (Part 3), refactored into
single-purpose functions that return a :class:`matplotlib.figure.Figure`.

The module forces the non-interactive ``Agg`` backend on import so it is safe
on a headless compute node.  Saving is handled by :func:`save_figure`; the
plotting functions themselves never write to disk.

References
----------
Bouwens et al. (2021), AJ 162, 47 — UVLF Schechter fit at z ~ 7
Finkelstein et al. (2015), ApJ 810, 71 — UVLF Schechter fit at z ~ 7
Song et al. (2016), ApJ 825, 5 — M_star–M_UV relation
Speagle et al. (2014), ApJS 214, 15 — star-forming main sequence
Schreiber et al. (2015), A&A 575, A74 — main sequence, high-z limit
Madau & Dickinson (2014), ARA&A 52, 415 — UV–SFR calibration
"""

from __future__ import annotations

import os
from typing import Tuple

import matplotlib

matplotlib.use("Agg")   # headless-safe; must precede pyplot import

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from scipy.ndimage import generic_filter

try:  # local package import (repo root on sys.path)
    from src.analysis import BiasEstimate, SNRResult
    from src.conversions import sfr_to_Muv
    from src.dataio import PowerSpectra, SimulationData
except ImportError:  # direct import of the module (src/ on sys.path)
    from analysis import BiasEstimate, SNRResult
    from conversions import sfr_to_Muv
    from dataio import PowerSpectra, SimulationData

__all__ = [
    "apply_plot_style",
    "save_figure",
    "fill_nan_nearest",
    "eor_colormap",
    "plot_halo_catalogue",
    "plot_sfr_relations",
    "plot_lightcone_fields",
    "plot_lightcone_slice",
    "plot_uv_luminosity_function",
    "plot_stellar_mass_muv",
    "plot_main_sequence",
    "plot_power_spectra",
    "plot_snr",
    "plot_bias_diagnostic",
]


# ===========================================================================
#  Style and output helpers
# ===========================================================================

def apply_plot_style(dpi: int = 200) -> None:
    """
    Apply the project-wide matplotlib defaults.

    Mirrors the rcParams block used by the notebooks, but with a
    file-output-appropriate DPI and constrained layout enabled globally.

    Parameters
    ----------
    dpi : int, optional
        Figure resolution for saved output.
    """
    plt.rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "font.size": 12,
        "axes.labelsize": 13,
        "legend.fontsize": 10,
        "axes.grid": False,
        "figure.constrained_layout.use": True,
    })


def save_figure(
    fig: Figure,
    output_dir: str,
    name: str,
    fmt: str = "png",
    close: bool = True,
) -> str:
    """
    Save a figure to ``output_dir/name.fmt`` and (optionally) close it.

    Parameters
    ----------
    fig : Figure
        Figure to write.
    output_dir : str
        Destination directory; created if missing.
    name : str
        Base filename without extension.
    fmt : str, optional
        File format understood by matplotlib (``png``, ``pdf``, ``svg``).
    close : bool, optional
        Close the figure afterwards to release memory.

    Returns
    -------
    str
        Path of the written file.
    """
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, f"{name}.{fmt}")
    fig.savefig(path, bbox_inches="tight")
    if close:
        plt.close(fig)
    return path


def fill_nan_nearest(arr: np.ndarray) -> np.ndarray:
    """
    Replace NaN pixels with the mean of their non-NaN neighbours.

    Cosmetic only — used so empty ``(k_perp, k_parallel)`` bins do not punch
    holes in the displayed power spectra.  Never use the result for science.

    Parameters
    ----------
    arr : ndarray
        2D array, possibly containing NaNs.

    Returns
    -------
    ndarray
        Copy of ``arr`` with NaNs filled.  Returned unchanged if there are no
        NaNs, or if every value is NaN.
    """
    mask = np.isnan(arr)
    if not mask.any() or mask.all():
        return arr

    filled = arr.copy()
    for _ in range(max(arr.shape)):
        if not np.isnan(filled).any():
            break
        kernel_filled = generic_filter(filled, np.nanmean, size=3, mode="nearest")
        filled = np.where(np.isnan(filled), kernel_filled, filled)
    return filled


def eor_colormap() -> LinearSegmentedColormap:
    """
    EoR-style colormap: dark (ionised, δT_b ≈ 0) → cyan → yellow → white.

    Returns
    -------
    LinearSegmentedColormap
        Colormap named ``"EoR21"``.
    """
    nodes = [
        (0.00, (0.00, 0.00, 0.15)),
        (0.12, (0.00, 0.10, 0.55)),
        (0.30, (0.00, 0.40, 0.85)),
        (0.50, (0.20, 0.78, 0.82)),
        (0.68, (0.95, 0.90, 0.22)),
        (0.85, (0.95, 0.42, 0.00)),
        (1.00, (0.97, 0.97, 0.97)),
    ]
    return LinearSegmentedColormap.from_list("EoR21", nodes, N=256)


def _halo_coords_mpc(data: SimulationData) -> np.ndarray:
    """
    Halo coordinates in Mpc, whether the catalogue stores grid cells or Mpc.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation, including the halo catalogue.

    Returns
    -------
    ndarray
        Coordinates of shape ``(N_halos, 3)`` in Mpc.
    """
    coords = data.halo_coords
    if coords.size == 0:
        return coords
    if coords.max() <= data.HII_DIM + 1:
        return coords * data.cell_size
    return coords


# ===========================================================================
#  Part 2 — fields and halo catalogue
# ===========================================================================

def plot_halo_catalogue(
    data: SimulationData,
    max_points: int = 50_000,
    seed: int = 42,
) -> Figure:
    """
    Halo catalogue overview: projected positions, mass function, and a slice.

    Three panels — projected ``(x, y)`` positions coloured by halo mass, the
    halo mass distribution, and halos in a thin slice over the 21 cm
    brightness-temperature field.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    max_points : int, optional
        Cap on the number of halos scattered per panel.
    seed : int, optional
        RNG seed for the subsampling.

    Returns
    -------
    Figure
        The assembled 1 × 3 figure.
    """
    rng = np.random.default_rng(seed)
    coords_mpc = _halo_coords_mpc(data)
    n_halos = coords_mpc.shape[0]

    plot_idx = rng.choice(n_halos, size=min(max_points, n_halos), replace=False)
    positions = coords_mpc[plot_idx]
    masses = data.halo_masses[plot_idx]

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # ── Panel 1: projected positions ──────────────────────────────────────
    scatter = axes[0].scatter(
        positions[:, 0], positions[:, 1],
        s=np.clip((np.log10(masses) - 7.0) ** 2, 1, 30),
        c=np.log10(masses),
        alpha=0.5,
    )
    fig.colorbar(scatter, ax=axes[0], label=r"$\log_{10}(M_{\rm halo}/M_\odot)$")
    axes[0].set_xlabel("x  [Mpc]")
    axes[0].set_ylabel("y  [Mpc]")
    axes[0].set_title(f"Projected halo catalogue at $z = {data.z_obs}$")
    axes[0].set_xlim(0, data.BOX_LEN)
    axes[0].set_ylim(0, data.BOX_LEN)
    axes[0].set_aspect("equal")

    # ── Panel 2: mass distribution ────────────────────────────────────────
    positive = data.halo_masses[data.halo_masses > 0]
    axes[1].hist(np.log10(positive), bins=60)
    axes[1].set_xlabel(r"$\log_{10}(M_{\rm halo}/M_\odot)$")
    axes[1].set_ylabel("Number of halos")
    axes[1].set_title("Halo mass distribution")

    # ── Panel 3: halos over a 21 cm slice ─────────────────────────────────
    slice_index = data.HII_DIM // 2
    slice_z_mpc = slice_index * data.cell_size
    slice_width_mpc = 2 * data.cell_size

    in_slice = np.abs(coords_mpc[:, 2] - slice_z_mpc) < slice_width_mpc
    slice_indices = np.where(in_slice)[0]
    if len(slice_indices) > 30_000:
        slice_indices = rng.choice(slice_indices, size=30_000, replace=False)

    image = axes[2].imshow(
        data.brightness_temp_field[:, :, slice_index].T,
        origin="lower",
        extent=[0, data.BOX_LEN, 0, data.BOX_LEN],
        interpolation="nearest",
        alpha=0.85,
    )
    fig.colorbar(image, ax=axes[2], label=r"$\delta T_b$  [mK]")
    axes[2].scatter(
        coords_mpc[slice_indices, 0], coords_mpc[slice_indices, 1],
        s=np.clip((np.log10(data.halo_masses[slice_indices]) - 7.0) ** 2, 1, 30),
        c="white", alpha=0.45, edgecolors="none",
    )
    axes[2].set_xlabel("x  [Mpc]")
    axes[2].set_ylabel("y  [Mpc]")
    axes[2].set_title(f"Halos over 21 cm slice at $z = {data.z_obs}$")
    axes[2].set_xlim(0, data.BOX_LEN)
    axes[2].set_ylim(0, data.BOX_LEN)
    axes[2].set_aspect("equal")

    return fig


def plot_sfr_relations(
    data: SimulationData,
    max_points: int = 100_000,
    seed: int = 42,
) -> Figure:
    """
    Halo SFR distribution and its scaling with halo and stellar mass.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    max_points : int, optional
        Cap on the number of halos in the scatter panels.
    seed : int, optional
        RNG seed for the subsampling.

    Returns
    -------
    Figure
        The assembled 1 × 3 figure.
    """
    rng = np.random.default_rng(seed)

    sfr = data.sfr
    halo_masses = data.halo_masses
    stellar_masses = data.stellar_masses

    positive_sfr = sfr[np.isfinite(sfr) & (sfr > 0)]

    valid = (
        np.isfinite(halo_masses) & np.isfinite(stellar_masses) & np.isfinite(sfr)
        & (halo_masses > 0) & (stellar_masses > 0) & (sfr > 0)
    )
    valid_idx = np.where(valid)[0]
    if len(valid_idx) > max_points:
        valid_idx = rng.choice(valid_idx, size=max_points, replace=False)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5))

    axes[0].hist(np.log10(positive_sfr), bins=80)
    axes[0].set_xlabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    axes[0].set_ylabel("Number of halos")
    axes[0].set_title("Halo catalogue SFR distribution")

    axes[1].scatter(
        np.log10(halo_masses[valid_idx]), np.log10(sfr[valid_idx]),
        s=2, alpha=0.25,
    )
    axes[1].set_xlabel(r"$\log_{10}(M_{\rm halo}/M_\odot)$")
    axes[1].set_ylabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    axes[1].set_title("SFR vs. halo mass")

    axes[2].scatter(
        np.log10(stellar_masses[valid_idx]), np.log10(sfr[valid_idx]),
        s=2, alpha=0.25,
    )
    axes[2].set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    axes[2].set_ylabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    axes[2].set_title("SFR vs. stellar mass")

    return fig


def plot_lightcone_fields(data: SimulationData) -> Figure:
    """
    Transverse and line-of-sight slices of the lightcone fields.

    Left: δT_b in a transverse ``(x, y)`` slice at mid-LOS.  Centre: δT_b
    along the LOS with a redshift axis on top.  Right: the neutral fraction
    on the same LOS slice.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation.

    Returns
    -------
    Figure
        The assembled 1 × 3 figure.
    """
    mid_z = data.N_z // 2
    mid_y = data.HII_DIM // 2

    extent_xy = [0, data.BOX_LEN, 0, data.BOX_LEN]
    extent_los = [data.lc_dist_Mpc[0], data.lc_dist_Mpc[-1], 0, data.BOX_LEN]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    vmax_t = np.percentile(np.abs(data.brightness_temp_field), 99)

    im0 = axes[0].imshow(
        data.brightness_temp_field[:, :, mid_z].T,
        origin="lower", extent=extent_xy,
        cmap="RdBu_r", vmin=-vmax_t, vmax=vmax_t,
    )
    fig.colorbar(im0, ax=axes[0], label="mK")
    axes[0].set_title(
        rf"$\delta T_b$ — transverse slice "
        rf"($z \approx {data.lc_redshifts[mid_z]:.2f}$)"
    )
    axes[0].set_xlabel("$x$  [Mpc]")
    axes[0].set_ylabel("$y$  [Mpc]")

    im1 = axes[1].imshow(
        data.brightness_temp_field[:, mid_y, :],
        origin="lower", extent=extent_los, aspect="auto",
        cmap="RdBu_r", vmin=-vmax_t, vmax=vmax_t,
    )
    fig.colorbar(im1, ax=axes[1], label="mK")
    axes[1].set_title(r"$\delta T_b$ — LOS slice (redshift evolution)")
    axes[1].set_xlabel("Comoving distance  [Mpc]")
    axes[1].set_ylabel("Transverse $x$  [Mpc]")

    ax_twin = axes[1].twiny()
    ax_twin.set_xlim(axes[1].get_xlim())
    z_ticks = np.linspace(data.z_min, data.z_max, 6)
    ax_twin.set_xticks(np.interp(z_ticks, data.lc_redshifts, data.lc_dist_Mpc))
    ax_twin.set_xticklabels([f"{z:.2f}" for z in z_ticks])
    ax_twin.set_xlabel("Redshift $z$")

    im2 = axes[2].imshow(
        data.neutral_fraction[:, mid_y, :],
        origin="lower", extent=extent_los, aspect="auto",
        cmap="bone", vmin=0, vmax=1,
    )
    fig.colorbar(im2, ax=axes[2], label=r"$x_{\rm HI}$")
    axes[2].set_title("Neutral fraction — LOS slice")
    axes[2].set_xlabel("Comoving distance  [Mpc]")
    axes[2].set_ylabel("Transverse $x$  [Mpc]")

    mean_xhi = float(np.mean(data.neutral_fraction))
    fig.suptitle(
        rf"Lightcone $z = {data.z_min}$–${data.z_max}$,  "
        rf"$\langle x_{{\rm HI}} \rangle = {mean_xhi:.2f}$",
        fontsize=14,
    )
    return fig


def plot_lightcone_slice(data: SimulationData) -> Figure:
    """
    Wide-format δT_b lightcone slice in the canonical EoR style.

    The x-axis runs along the line of sight from the observer (left,
    ``z_min``) toward high redshift (right, ``z_max``); the y-axis is one
    transverse direction.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation.

    Returns
    -------
    Figure
        A single wide panel with comoving distance below and redshift above.
    """
    mid_y = data.HII_DIM // 2
    t_slice = data.brightness_temp_field[:, mid_y, :]

    # Clip at zero: the emission-only regime at z ~ 6–8.
    t_lo = max(0.0, float(np.percentile(t_slice, 1)))
    t_hi = float(np.percentile(t_slice, 99.5))

    fig, ax = plt.subplots(figsize=(16, 3.5))
    image = ax.imshow(
        t_slice,
        origin="lower",
        extent=[data.lc_dist_Mpc[0], data.lc_dist_Mpc[-1], 0, data.BOX_LEN],
        aspect="auto",
        cmap=eor_colormap(),
        vmin=t_lo, vmax=t_hi,
    )
    cbar = fig.colorbar(image, ax=ax, fraction=0.015, pad=0.01)
    cbar.set_label(r"$\delta T_b$  [mK]", fontsize=12)

    ax.set_xlabel("Comoving distance  [Mpc]", fontsize=12)
    ax.set_ylabel("Transverse  $x$  [Mpc]", fontsize=12)

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    z_ticks = np.linspace(data.z_min, data.z_max, 7)
    ax_top.set_xticks(np.interp(z_ticks, data.lc_redshifts, data.lc_dist_Mpc))
    ax_top.set_xticklabels([f"{z:.2f}" for z in z_ticks])
    ax_top.set_xlabel("Redshift  $z$", fontsize=12)

    f_21_mhz = data.get("F_21_MHZ", 1420.405)
    f_lo = f_21_mhz / (1 + data.z_max)
    f_hi = f_21_mhz / (1 + data.z_min)
    mean_xhi = float(np.mean(data.neutral_fraction))

    ax.set_title(
        rf"21 cm brightness temperature lightcone slice   "
        rf"($z = {data.z_min}$–${data.z_max}$,  "
        rf"$f_{{\rm obs}} = {f_lo:.0f}$–${f_hi:.0f}$ MHz,  "
        rf"$\langle x_{{\rm HI}} \rangle = {mean_xhi:.2f}$)",
        fontsize=12, pad=14,
    )
    return fig


# ===========================================================================
#  Part 2 — literature comparison
# ===========================================================================

def _schechter_muv(
    magnitude: np.ndarray,
    phi_star: float,
    m_star: float,
    alpha: float,
) -> np.ndarray:
    """
    Schechter function in UV-magnitude space.

    Parameters
    ----------
    magnitude : ndarray
        Absolute UV magnitudes.
    phi_star : float
        Normalisation [Mpc^-3].
    m_star : float
        Characteristic magnitude.
    alpha : float
        Faint-end slope.

    Returns
    -------
    ndarray
        Φ(M_UV) [Mpc^-3 mag^-1].

    References
    ----------
    Schechter (1976), ApJ 203, 297.
    """
    x = 10.0 ** (0.4 * (m_star - magnitude))
    return (np.log(10.0) / 2.5) * phi_star * x ** (alpha + 1.0) * np.exp(-x)


def plot_uv_luminosity_function(
    data: SimulationData,
    magnitude_bin_width: float = 0.5,
) -> Figure:
    """
    Simulated UV luminosity function against z ~ 7 Schechter fits.

    Halo SFRs are converted to ``M_UV`` with the Madau & Dickinson (2014)
    calibration and binned over the coeval box volume ``BOX_LEN^3``.  When the
    catalogue was subsampled on load, counts are rescaled by
    ``data.halo_sampling_factor``.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    magnitude_bin_width : float, optional
        Bin width in magnitudes.

    Returns
    -------
    Figure
        UVLF figure with Bouwens+21 and Finkelstein+15 overlaid.
    """
    literature = [
        dict(label=r"Bouwens+21  $z \sim 7$", phi_star=0.19e-3,
             m_star=-21.15, alpha=-2.06, color="royalblue", ls="--"),
        dict(label=r"Finkelstein+15  $z \sim 7$", phi_star=1.57e-4,
             m_star=-21.03, alpha=-2.03, color="tomato", ls="-."),
    ]

    sfr_valid = data.sfr[np.isfinite(data.sfr) & (data.sfr > 0)]
    muv_sim = sfr_to_Muv(sfr_valid)

    edges = np.arange(-25.0, -10.0 + magnitude_bin_width, magnitude_bin_width)
    centres = 0.5 * (edges[:-1] + edges[1:])

    counts, _ = np.histogram(muv_sim, bins=edges)
    counts_scaled = counts * data.halo_sampling_factor

    volume = data.BOX_LEN ** 3
    phi_sim = counts_scaled / (volume * magnitude_bin_width)
    phi_err = (
        np.sqrt(np.maximum(counts, 1)) * data.halo_sampling_factor
        / (volume * magnitude_bin_width)
    )
    keep = counts > 0

    m_uv_limit = data.get("M_UV_limit", -18.0)

    fig, ax = plt.subplots(figsize=(8, 6))

    magnitudes = np.linspace(-25.5, -11.0, 600)
    for entry in literature:
        ax.semilogy(
            magnitudes,
            _schechter_muv(magnitudes, entry["phi_star"], entry["m_star"], entry["alpha"]),
            color=entry["color"], ls=entry["ls"], lw=2.0, label=entry["label"],
        )

    ax.set_yscale("log")
    ax.errorbar(
        centres[keep], phi_sim[keep],
        yerr=[np.minimum(phi_err[keep], phi_sim[keep] * 0.9999), phi_err[keep]],
        fmt="ko", ms=5, capsize=3, capthick=1.2, elinewidth=1.0, zorder=6,
        label=rf"21cmFAST ($z = {data.z_obs:.1f}$)",
    )

    ax.axvline(m_uv_limit, color="dimgray", ls="--", lw=1.5,
               label=rf"Euclid limit ($M_{{\rm UV}} = {m_uv_limit:.0f}$)")
    ax.axvspan(m_uv_limit, -11.0, color="gray", alpha=0.08, zorder=0,
               label="Fainter than Euclid limit")

    ax.set_xlim(-25.5, -11.0)
    ax.set_ylim(1e-8, 2e-1)
    ax.set_xlabel(r"$M_{\rm UV}$  [AB mag]", fontsize=13)
    ax.set_ylabel(r"$\Phi(M_{\rm UV})$  [Mpc$^{-3}$ mag$^{-1}$]", fontsize=13)
    ax.legend(fontsize=10, loc="lower right", framealpha=0.88)
    ax.set_title(
        rf"UV luminosity function  ($z = {data.z_obs:.1f}$, "
        rf"$V_{{\rm box}} = {volume:.1e}$ Mpc$^3$)",
        fontsize=12,
    )

    n_bright = int(np.sum(muv_sim <= m_uv_limit) * data.halo_sampling_factor)
    n_total = int(len(muv_sim) * data.halo_sampling_factor)
    subsample_note = (
        "" if data.halo_sampling_factor == 1.0
        else rf"  (catalogue subsampled $\times 1/{data.halo_sampling_factor:.0f}$)"
    )
    ax.text(
        0.03, 0.04,
        rf"$N_{{\rm Euclid-bright}}$ = {n_bright:,} / $N_{{\rm total}}$ = "
        rf"{n_total:,}{subsample_note}",
        transform=ax.transAxes, fontsize=9, color="dimgray",
    )
    return fig


def _binned_median(
    x: np.ndarray,
    y: np.ndarray,
    edges: np.ndarray,
    min_count: int = 10,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Median and 16th/84th percentiles of ``y`` in bins of ``x``.

    Parameters
    ----------
    x, y : ndarray
        Paired samples of equal length.
    edges : ndarray
        Bin edges along ``x``.
    min_count : int, optional
        Minimum samples for a bin to be reported; sparser bins give NaN.

    Returns
    -------
    centres : ndarray
        Bin centres.
    median, p16, p84 : ndarray
        Per-bin statistics, NaN where the bin holds fewer than ``min_count``.
    """
    centres = 0.5 * (edges[:-1] + edges[1:])
    median = np.full(centres.shape, np.nan)
    p16 = np.full(centres.shape, np.nan)
    p84 = np.full(centres.shape, np.nan)

    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        values = y[(x >= lo) & (x < hi)]
        if values.size >= min_count:
            median[i] = np.median(values)
            p16[i] = np.percentile(values, 16)
            p84[i] = np.percentile(values, 84)

    return centres, median, p16, p84


def plot_stellar_mass_muv(
    data: SimulationData,
    max_points: int = 100_000,
    seed: int = 42,
) -> Figure:
    """
    Stellar mass versus UV magnitude, against two z ~ 7 relations.

    Overlays Song et al. (2016) and González et al. (2010); the latter lies
    ~0.2 dex higher because of its constant-star-formation-history assumption.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    max_points : int, optional
        Cap on the number of scattered halos.
    seed : int, optional
        RNG seed for the subsampling.

    Returns
    -------
    Figure
        Scatter plus binned median and 16–84 percentile band.
    """
    valid = (
        np.isfinite(data.stellar_masses) & np.isfinite(data.sfr)
        & (data.stellar_masses > 0) & (data.sfr > 0)
    )
    muv_halos = sfr_to_Muv(data.sfr[valid])
    log_mstar = np.log10(data.stellar_masses[valid])

    edges = np.arange(-24.5, -13.0 + 0.5, 0.5)
    centres, median, p16, p84 = _binned_median(muv_halos, log_mstar, edges)
    show = np.isfinite(median)

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(muv_halos), size=min(max_points, len(muv_halos)), replace=False)

    m_uv_limit = data.get("M_UV_limit", -18.0)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(muv_halos[idx], log_mstar[idx], s=1, alpha=0.04, color="0.2",
               rasterized=True, label="21cmFAST galaxies")
    ax.fill_between(centres[show], p16[show], p84[show], color="0.3", alpha=0.20,
                    lw=0, label="21cmFAST 16–84 percentile")
    ax.plot(centres[show], median[show], "o-", color="k", ms=5, lw=2,
            label=rf"21cmFAST median ($z={data.z_obs:.1f}$)")

    muv_lit = np.linspace(-24.5, -13.5, 300)
    ax.plot(muv_lit, 8.86 - 0.5 * (muv_lit + 20.0), color="royalblue", ls="--",
            lw=2, label=r"Song+16 $z\sim7$")
    # González+10 sits ~0.2 dex above Song+16: their constant-star-formation-
    # history assumption raises the inferred stellar masses.
    ax.plot(muv_lit, 9.06 - 0.5 * (muv_lit + 20.0), color="tomato", ls="-.",
            lw=2, label="González+10 " + r"$z\sim7$")

    ax.axvline(m_uv_limit, color="dimgray", ls="--", lw=1.5,
               label=rf"Euclid limit ($M_{{\rm UV}}={m_uv_limit:.0f}$)")
    ax.axvspan(m_uv_limit, -13.0, color="gray", alpha=0.08, zorder=0,
               label="Fainter than Euclid limit")

    ax.set_xlim(-24.5, -13.5)
    ax.set_ylim(5.5, 11.5)
    ax.set_xlabel(r"$M_{\rm UV}$  [AB mag]", fontsize=13)
    ax.set_ylabel(r"$\log_{10}(M_\star/M_\odot)$", fontsize=13)
    ax.set_title(rf"Stellar mass–UV magnitude relation ($z={data.z_obs:.1f}$)",
                 fontsize=12)
    ax.legend(fontsize=9, loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.25)
    return fig


def plot_main_sequence(
    data: SimulationData,
    t_star: float = 0.5,
) -> Figure:
    """
    Star-forming main sequence against Speagle+14 and Schreiber+15.

    The green line is the 21cmFAST model prediction ``SFR = M_star / t_sf``
    with ``t_sf = t_STAR × t_H(z)``; the offset from the observational
    calibrations is a model difference, not a bug (see ``docs/Low_SFR_fix.md``).

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    t_star : float, optional
        21cmFAST ``t_STAR`` parameter (fraction of the Hubble time).

    Returns
    -------
    Figure
        Binned median SFR with the 16–84 percentile band and literature lines.
    """
    from astropy.cosmology import Planck18
    import astropy.units as u

    t_hubble_yr = float((1.0 / Planck18.H(data.z_obs)).to(u.yr).value)
    t_sf_yr = t_star * t_hubble_yr
    t_age_gyr = float(Planck18.age(data.z_obs).to(u.Gyr).value)

    valid = (
        np.isfinite(data.stellar_masses) & np.isfinite(data.sfr)
        & (data.stellar_masses > 0) & (data.sfr > 0)
    )
    log_mstar = np.log10(data.stellar_masses[valid])
    log_sfr = np.log10(data.sfr[valid])

    edges = np.arange(4.0, 11.5 + 0.5, 0.5)
    centres, median, p16, p84 = _binned_median(log_mstar, log_sfr, edges, min_count=5)
    show = np.isfinite(median)

    logm_lit = np.linspace(4.5, 11.5, 300)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.fill_between(centres[show], p16[show], p84[show], color="0.3", alpha=0.18,
                    zorder=2, label="16–84 % scatter")
    ax.plot(centres[show], median[show], "o-", color="k", ms=5, lw=2, zorder=3,
            label=rf"21cmFAST median ($z = {data.z_obs:.1f}$)")

    ax.plot(logm_lit, logm_lit - np.log10(t_sf_yr), color="forestgreen", ls="-",
            lw=2.0,
            label=rf"21cmFAST model ($t_\star = {t_star}\,t_H = {t_sf_yr/1e6:.0f}$ Myr)")
    ax.plot(logm_lit,
            (0.84 - 0.026 * t_age_gyr) * logm_lit - (6.51 - 0.11 * t_age_gyr),
            color="royalblue", ls="--", lw=2,
            label=rf"Speagle+14  ($t_\mathrm{{age}} = {t_age_gyr:.2f}$ Gyr)")
    ax.plot(logm_lit, logm_lit - 8.0, color="tomato", ls="-.", lw=2,
            label=r"Schreiber+15  (sSFR $= 10$ Gyr$^{-1}$)")

    ax.text(
        0.03, 0.97,
        rf"All 21cmFAST halos: $N = {int(valid.sum()):,}$" "\n"
        rf"Model offset: $\Delta\log_{{10}}\mathrm{{SFR}} \approx "
        rf"{np.log10(1 / t_sf_yr) - np.log10(1e-8):.2f}$ dex (expected)",
        transform=ax.transAxes, fontsize=8.5, color="dimgray", va="top",
    )

    ax.set_xlim(4.5, 11.5)
    ax.set_xlabel(r"$\log_{10}(M_\star\,/\,M_\odot)$", fontsize=13)
    ax.set_ylabel(r"$\log_{10}(\mathrm{SFR}\,/\,M_\odot\,\mathrm{yr}^{-1})$",
                  fontsize=13)
    ax.legend(fontsize=9.5, loc="best", framealpha=0.88)
    ax.set_title(rf"Star-forming main sequence  ($z = {data.z_obs:.1f}$)",
                 fontsize=12)
    return fig


# ===========================================================================
#  Part 3 — power spectra and SNR
# ===========================================================================

def _add_wedge_lines(
    ax: plt.Axes,
    k_perp: np.ndarray,
    horizon_slope: float,
    fov_slope: float,
    color: str = "w",
    label: bool = False,
) -> None:
    """
    Overlay the horizon and primary-beam wedge lines on a k-space panel.

    Parameters
    ----------
    ax : Axes
        Target axes.
    k_perp : ndarray
        Transverse bin centres [Mpc^-1].
    horizon_slope, fov_slope : float
        Wedge slopes from ``src.analysis``.
    color : str, optional
        Line colour.
    label : bool, optional
        Add legend entries.
    """
    line = np.logspace(np.log10(k_perp.min()), np.log10(k_perp.max()), 200)
    ax.plot(line, line * horizon_slope, color=color, ls="-", lw=1.2, alpha=0.7,
            label="Horizon" if label else None)
    ax.plot(line, line * fov_slope, color=color, ls="--", lw=1.2, alpha=0.8,
            label="HERA FoV wedge" if label else None)


def _style_k_axes(ax: plt.Axes, k_perp: np.ndarray, k_parallel: np.ndarray) -> None:
    """
    Apply log axes, labels, and limits to a ``(k_perp, k_parallel)`` panel.

    Parameters
    ----------
    ax : Axes
        Target axes.
    k_perp, k_parallel : ndarray
        Bin centres [Mpc^-1].
    """
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$k_\perp$  [Mpc$^{-1}$]")
    ax.set_ylabel(r"$k_\parallel$  [Mpc$^{-1}$]")
    ax.set_xlim(k_perp[0], k_perp[-1])
    ax.set_ylim(k_parallel[0], k_parallel[-2])


def _signed_log(power: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Sign-preserving log of a power spectrum, plus a symmetric colour limit.

    Parameters
    ----------
    power : ndarray
        Power spectrum, already transposed for display.

    Returns
    -------
    signed_log : ndarray
        ``sign(P) × log10|P|``, with empty bins filled for display.
    clim : float
        95th-percentile absolute value, for symmetric colour limits.
    """
    amplitude = np.abs(power)
    with np.errstate(divide="ignore", invalid="ignore"):
        raw = np.where(amplitude > 0, np.log10(amplitude) * np.sign(power), 0.0)
    filled = fill_nan_nearest(raw)
    nonzero = filled[filled != 0]
    clim = float(np.nanpercentile(np.abs(nonzero), 95)) if nonzero.size else 1.0
    return filled, clim


def plot_power_spectra(
    spectra: PowerSpectra,
    data: SimulationData,
    horizon_slope: float,
    fov_slope: float,
) -> Figure:
    """
    The three 2D cylindrical power spectra with wedge lines overlaid.

    Parameters
    ----------
    spectra : PowerSpectra
        Computed 21 cm auto-, galaxy auto-, and cross-spectra.
    data : SimulationData
        Loaded simulation, for the title metadata.
    horizon_slope, fov_slope : float
        Wedge slopes from ``src.analysis``.

    Returns
    -------
    Figure
        The assembled 1 × 3 figure.
    """
    k_perp = spectra.k_perp
    k_parallel = spectra.k_parallel

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    with np.errstate(divide="ignore", invalid="ignore"):
        p21_display = fill_nan_nearest(np.log10(np.abs(spectra.P_21cm_auto.T)))
        pgal_display = fill_nan_nearest(np.log10(np.abs(spectra.P_galaxy_auto.T)))

    im0 = axes[0].pcolormesh(k_perp, k_parallel, p21_display, cmap="viridis",
                             shading="auto")
    fig.colorbar(im0, ax=axes[0], label=r"$\log_{10} |P_{21}|$  [mK² Mpc³]")
    _add_wedge_lines(axes[0], k_perp, horizon_slope, fov_slope, "w", label=True)
    axes[0].legend(loc="upper left", fontsize=8)
    axes[0].set_title(r"$P_{21}(k_\perp,\, k_\parallel)$")

    im1 = axes[1].pcolormesh(k_perp, k_parallel, pgal_display, cmap="plasma",
                             shading="auto")
    fig.colorbar(im1, ax=axes[1], label=r"$\log_{10} |P_{\rm gal}|$  [Mpc³]")
    _add_wedge_lines(axes[1], k_perp, horizon_slope, fov_slope, "w")
    axes[1].set_title(r"$P_{\rm gal}(k_\perp,\, k_\parallel)$")

    signed_log, clim = _signed_log(spectra.P_cross.T)
    im2 = axes[2].pcolormesh(k_perp, k_parallel, signed_log, cmap="RdBu_r",
                             shading="auto", vmin=-clim, vmax=clim)
    fig.colorbar(im2, ax=axes[2],
                 label=r"sign $\times$ $\log_{10} |P_{21 \times \rm gal}|$")
    _add_wedge_lines(axes[2], k_perp, horizon_slope, fov_slope, "k")
    axes[2].set_title(r"$P_{21 \times \rm gal}$  (signed)")

    for ax in axes:
        _style_k_axes(ax, k_perp, k_parallel)

    mean_xhi = float(np.mean(data.neutral_fraction))
    fig.suptitle(
        rf"2D cylindrical power spectra — lightcone $z = {data.z_min}$–${data.z_max}$   "
        rf"($\langle x_{{\rm HI}} \rangle = {mean_xhi:.2f}$,  "
        rf"Euclid $M_{{\rm UV}} < {data.get('M_UV_limit', -18.0):.0f}$)",
        fontsize=13,
    )
    return fig


def plot_snr(
    spectra: PowerSpectra,
    snr: SNRResult,
    P_cross_observed: np.ndarray,
    data: SimulationData,
    horizon_slope: float,
    fov_slope: float,
) -> Figure:
    """
    Per-mode SNR map and the photo-z damped cross-power spectrum.

    Parameters
    ----------
    spectra : PowerSpectra
        Provides the ``(k_perp, k_parallel)`` grid.
    snr : SNRResult
        Output of ``src.analysis.cross_power_snr``.
    P_cross_observed : ndarray
        Photo-z damped cross-power spectrum.
    data : SimulationData
        Loaded simulation, for the title metadata.
    horizon_slope, fov_slope : float
        Wedge slopes from ``src.analysis``.

    Returns
    -------
    Figure
        The assembled 1 × 2 figure.
    """
    k_perp = spectra.k_perp
    k_parallel = spectra.k_parallel

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    with np.errstate(divide="ignore", invalid="ignore"):
        snr_display = fill_nan_nearest(np.log10(snr.snr_per_mode.T))

    im_snr = axes[0].pcolormesh(k_perp, k_parallel, snr_display, cmap="magma",
                                shading="auto", vmin=-3, vmax=0)
    fig.colorbar(im_snr, ax=axes[0], label=r"$\log_{10}$ SNR per mode")
    _add_wedge_lines(axes[0], k_perp, horizon_slope, fov_slope, "w", label=True)
    axes[0].legend(loc="upper left", fontsize=8)
    axes[0].set_title(rf"Per-mode SNR  (total = {snr.total_snr:.1f}$\sigma$)")

    signed_log, clim = _signed_log(P_cross_observed.T)
    im_cross = axes[1].pcolormesh(k_perp, k_parallel, signed_log, cmap="RdBu_r",
                                  shading="auto", vmin=-clim, vmax=clim)
    fig.colorbar(im_cross, ax=axes[1],
                 label=r"sign $\times$ $\log_{10}|P_{21 \times \rm gal}^{\rm obs}|$")
    _add_wedge_lines(axes[1], k_perp, horizon_slope, fov_slope, "k")
    axes[1].set_title(r"Observed $P_{21 \times \rm gal}$  (photo-$z$ damped)")

    for ax in axes:
        _style_k_axes(ax, k_perp, k_parallel)

    fig.suptitle(
        rf"HERA $\times$ Euclid — lightcone $z = {data.z_min}$–${data.z_max}$   "
        rf"($\sigma_z = {data.get('photoz_uncertainty', 0.059)}$)",
        fontsize=13,
    )
    return fig


def plot_bias_diagnostic(bias: BiasEstimate, z_obs: float) -> Figure:
    """
    Mass distribution of the Euclid-selected halos with the bias curve.

    Parameters
    ----------
    bias : BiasEstimate
        Output of ``src.analysis.effective_galaxy_bias``.
    z_obs : float
        Reference redshift, for the title.

    Returns
    -------
    Figure
        Histogram of selected halo masses with ``b_h(M, z)`` on a twin axis.
    """
    fig, ax1 = plt.subplots(figsize=(7, 5))

    ax1.hist(bias.log10_selected_mass_h, bins=50, alpha=0.7, color="steelblue")
    ax1.set_xlabel(r"$\log_{10}(M_{\rm halo}\ /\ M_\odot\,h^{-1})$")
    ax1.set_ylabel("Number of selected halos")

    ax2 = ax1.twinx()
    ax2.plot(bias.log10_mass_grid, bias.bias_grid, "r--", lw=1.5,
             label=r"$b_h(M,z)$")
    ax2.set_ylabel(r"Halo bias $b_h(M,z)$", color="r")
    ax2.tick_params(axis="y", labelcolor="r")
    ax2.axhline(bias.mean_bias, color="r", lw=1, alpha=0.4, linestyle=":")

    ax1.set_title(
        rf"Euclid-selected halos and bias at $z = {z_obs}$  "
        rf"($\langle b_g \rangle = {bias.mean_bias:.2f}$)"
    )
    return fig
