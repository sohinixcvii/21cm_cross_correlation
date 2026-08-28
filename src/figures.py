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
from typing import Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")   # headless-safe; must precede pyplot import

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm, TwoSlopeNorm
from matplotlib.figure import Figure
from scipy.ndimage import gaussian_filter, generic_filter

try:  # local package import (repo root on sys.path)
    from src.analysis import (
        BiasEstimate,
        EuclidSelection,
        SNRResult,
        UncertaintyBudget,
        deposit_halo_field,
        foreground_wedge_mask,
        galaxy_overdensity_from_catalogue,
        photoz_damping_kernel,
        select_euclid_halos,
    )
    from src.analysis import (
        SphericalSpectra, NumberDensity,
        spherically_average_spectra, comoving_number_density,
    )
    from src.conversions import sfr_to_Luv, sfr_to_Muv
    from src.dataio import PowerSpectra, SimulationData
except ImportError:  # direct import of the module (src/ on sys.path)
    from analysis import (
        BiasEstimate,
        EuclidSelection,
        SNRResult,
        UncertaintyBudget,
        deposit_halo_field,
        foreground_wedge_mask,
        galaxy_overdensity_from_catalogue,
        photoz_damping_kernel,
        select_euclid_halos,
    )
    from analysis import (
        SphericalSpectra, NumberDensity,
        spherically_average_spectra, comoving_number_density,
    )
    from conversions import sfr_to_Luv, sfr_to_Muv
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
    "plot_uv_selection_maps",
    "selected_galaxy_overdensity",
    "plot_euclid_selected_catalogue",
    "plot_selected_galaxy_overdensity",
    "plot_galaxy_overdensity_on_21cm",
    "plot_power_spectra",
    "plot_galaxy_wedge",
    "plot_wedge_real_space",
    "plot_photoz_suppression",
    "plot_snr",
    "plot_uncertainty_budget",
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


#: Contour colour for overlays on a diverging (RdBu_r) map.  Gold is high
#: luminance against *both* saturated ends of that colormap, where black
#: disappears into the blue and white disappears into the pale centre.
CONTOUR_COLOR_ON_DIVERGING = "#FFC61E"

#: Contour colour for overlays on the EoR map, whose low end is near-black.
CONTOUR_COLOR_ON_EOR = "#FF3B5C"


def _outline(artist_contour, width: float = 2.0) -> None:
    """
    Give contour lines and their labels a dark outline.

    A single colour cannot contrast with every part of a diverging or
    multi-hue colormap, so the lines carry a stroke instead: the outline
    holds the line legible over the light middle of the scale, and the fill
    holds it legible over the dark ends.

    Parameters
    ----------
    artist_contour : QuadContourSet
        The contour set returned by ``Axes.contour``.
    width : float, optional
        Total stroke width in points; the line sits inside it.
    """
    stroke = [patheffects.withStroke(linewidth=width, foreground="0.15")]
    # Matplotlib >= 3.8 makes ContourSet an artist in its own right and drops
    # the `.collections` list; older versions need the per-collection loop.
    if hasattr(artist_contour, "set_path_effects"):
        artist_contour.set_path_effects(stroke)
    else:                                            # pragma: no cover - old mpl
        for collection in artist_contour.collections:
            collection.set_path_effects(stroke)


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
    # The 21 cm field is a *lightcone*: its third axis has N_z cells, which is
    # not HII_DIM.  The halo catalogue is *coeval* and spans BOX_LEN along z.
    # The two share the (x, y) plane but not the line-of-sight scale
    # (docs/HPC.md §6), so each is taken at its own midpoint rather than at a
    # shared index.  Indexing the lightcone with HII_DIM // 2 happened to stay
    # in range whenever N_z > HII_DIM / 2, and raised IndexError as soon as it
    # did not — see tests/test_figures.py.
    slice_index = data.brightness_temp_field.shape[2] // 2
    slice_z_mpc = 0.5 * data.BOX_LEN
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
    M_UV_bright: Optional[float] = None,
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
    # Literature Schechter fits at three redshifts.  The SIMULATED LF is a
    # single curve at z_obs and cannot be split into redshift bins: the stored
    # halo catalogue is the *coeval* catalogue at z_obs, with no per-halo
    # redshift, so slicing it by line-of-sight position would manufacture
    # evolution the data does not contain.  These literature curves supply the
    # redshift comparison instead, and are labelled as such.
    literature = [
        dict(label=r"Bouwens+21  $z \sim 6$", phi_star=0.29e-3,
             m_star=-20.94, alpha=-1.93, color="seagreen", ls=":"),
        dict(label=r"Bouwens+21  $z \sim 7$", phi_star=0.19e-3,
             m_star=-21.15, alpha=-2.06, color="royalblue", ls="--"),
        dict(label=r"Bouwens+21  $z \sim 8$", phi_star=0.088e-3,
             m_star=-21.03, alpha=-2.23, color="darkorange", ls=(0, (5, 1))),
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
    m_uv_bright = (
        float(M_UV_bright) if M_UV_bright is not None
        else float(data.get("M_UV_bright", -22.0))
    )

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

    # Both edges of the selection window, and the window itself.  Marking
    # only the faint cut hid the fact that the bright cut also removes
    # galaxies -- which matters a great deal once the window is narrow.
    ax.axvline(m_uv_limit, color="dimgray", ls="--", lw=1.5,
               label=rf"Faint cut ($M_{{\rm UV}} = {m_uv_limit:.1f}$)")
    ax.axvline(m_uv_bright, color="darkslateblue", ls=":", lw=1.8,
               label=rf"Bright cut ($M_{{\rm UV}} = {m_uv_bright:.1f}$)")
    ax.axvspan(m_uv_bright, m_uv_limit, color="mediumseagreen", alpha=0.13,
               zorder=0, label="Selection window")
    ax.axvspan(m_uv_limit, -11.0, color="gray", alpha=0.08, zorder=0,
               label="Fainter than the faint cut")

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

def plot_uv_selection_maps(
    data: SimulationData,
    M_UV_bright: float = -22.0,
    n_grid: int = 128,
) -> Figure:
    """
    Where the Euclid-selected galaxies actually sit, and how bright they are.

    The scaling-relation figures show *what* the selection keeps; this one
    shows *where*.  Three panels — the projected UV luminosity of the whole
    catalogue, the projected counts of the magnitude-selected sample, and the
    selected magnitude distribution against the cuts that produced it.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    M_UV_bright : float, optional
        Bright-end magnitude cut.  The faint end comes from the stored
        ``M_UV_limit`` attribute.
    n_grid : int, optional
        Transverse cells per side for the projected maps.

    Returns
    -------
    Figure
        The assembled 1 x 3 figure.
    """
    M_UV_faint = float(data.get("M_UV_limit", -18.0))

    coords = _halo_coords_mpc(data)
    sfr = np.asarray(data.sfr, dtype=float)

    # The magnitude window is applied by src.analysis, the same call the bias
    # stage makes, so this figure and the bias estimate select identically.
    selection = select_euclid_halos(
        sfr, np.asarray(data.halo_masses, dtype=float),
        M_UV_faint=M_UV_faint, M_UV_bright=M_UV_bright,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        luminosity = sfr_to_Luv(sfr)
        magnitude = sfr_to_Muv(sfr)

    finite = np.isfinite(luminosity) & (luminosity > 0)

    # The projected maps are the LOS sum of the same 3D deposit that builds
    # delta_gal (analysis.deposit_halo_field), so figure and field bin
    # identically.  n_los=1 collapses the LOS axis at bin time.
    def _project(mask: np.ndarray, weights: Optional[np.ndarray]) -> np.ndarray:
        return deposit_halo_field(
            coords[mask], box_len=data.BOX_LEN, n_perp=n_grid, n_los=1,
            weights=weights,
        )[:, :, 0]

    luminosity_map = _project(finite, luminosity[finite])

    # selection.mask indexes the valid (SFR > 0) subset; lift it back to
    # full-catalogue indices to match `coords`.
    selected_full = _lift_selection_mask(
        sfr, np.asarray(data.halo_masses, dtype=float), selection,
    )

    counts_map = _project(selected_full, None)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
    extent_xy = [0, data.BOX_LEN, 0, data.BOX_LEN]

    with np.errstate(divide="ignore", invalid="ignore"):
        luminosity_display = np.log10(np.where(luminosity_map > 0, luminosity_map, np.nan))

    im0 = axes[0].imshow(luminosity_display.T, origin="lower", extent=extent_xy,
                         cmap="inferno")
    fig.colorbar(im0, ax=axes[0], label=r"$\log_{10} \sum L_{\rm UV}$  [erg s$^{-1}$ Hz$^{-1}$]")
    axes[0].set_title("Projected UV luminosity — all halos")

    im1 = axes[1].imshow(counts_map.T, origin="lower", extent=extent_xy, cmap="viridis")
    fig.colorbar(im1, ax=axes[1], label="Selected galaxies per cell")
    axes[1].set_title(
        rf"Projected selected counts  (${M_UV_bright} \leq M_{{\rm UV}} \leq {M_UV_faint}$)"
    )

    for ax in axes[:2]:
        ax.set_xlabel("$x$  [Mpc]")
        ax.set_ylabel("$y$  [Mpc]")

    selected_mag = magnitude[selected_full]
    if selected_mag.size:
        axes[2].hist(selected_mag, bins=40, color="steelblue", alpha=0.8)
    axes[2].axvline(M_UV_faint, color="crimson", ls="--", lw=1.5, label="faint cut")
    axes[2].axvline(M_UV_bright, color="darkorange", ls="--", lw=1.5, label="bright cut")
    axes[2].set_xlabel(r"$M_{\rm UV}$")
    axes[2].set_ylabel("Number of galaxies")
    axes[2].legend(fontsize=8)
    axes[2].set_title("Selected magnitude distribution")

    fig.suptitle(
        rf"Euclid selection — {selection.n_selected:,} of {selection.n_valid:,} halos "
        rf"at $z_{{\rm obs}} = {data.z_obs}$",
        fontsize=13,
    )
    return fig


def _lift_selection_mask(
    sfr: np.ndarray,
    halo_masses: np.ndarray,
    selection: EuclidSelection,
) -> np.ndarray:
    """
    Expand a Euclid selection mask to full-catalogue indices.

    :func:`src.analysis.select_euclid_halos` masks *within* the valid
    (``SFR > 0`` and ``M > 0``) subset, so its mask is shorter than the
    catalogue and cannot index ``halo_coords`` directly.

    Parameters
    ----------
    sfr : ndarray
        Per-halo star-formation rate [M_sun yr^-1], full catalogue.
    halo_masses : ndarray
        Per-halo mass [M_sun], full catalogue.
    selection : EuclidSelection
        Selection returned for the same catalogue.

    Returns
    -------
    ndarray of bool
        Mask of shape ``(N_halos,)``, True for the selected halos.
    """
    sfr = np.asarray(sfr, dtype=float)
    halo_masses = np.asarray(halo_masses, dtype=float)

    valid = np.isfinite(sfr) & (sfr > 0) & (halo_masses > 0)
    full_mask = np.zeros(sfr.shape[0], dtype=bool)
    full_mask[np.flatnonzero(valid)[selection.mask]] = True
    return full_mask


def _display_contour_levels(
    field: np.ndarray,
    quantiles: Sequence[float] = (0.80, 0.90, 0.97),
) -> np.ndarray:
    """
    Strictly increasing contour levels at quantiles of a 2D field.

    Returns an empty array when the field is degenerate (all-NaN, or too
    few distinct values to place a contour), which the callers treat as
    "draw no contours" rather than an error.

    Parameters
    ----------
    field : ndarray
        2D array to take quantiles of.
    quantiles : sequence of float, optional
        Quantiles in ``[0, 1]``.

    Returns
    -------
    ndarray
        Sorted, unique, finite levels; possibly empty.
    """
    finite = field[np.isfinite(field)]
    if finite.size == 0 or np.ptp(finite) <= 0:
        return np.array([])
    levels = np.unique(np.quantile(finite, quantiles))
    return levels[np.isfinite(levels)]


def selected_galaxy_overdensity(
    data: SimulationData,
    M_UV_bright: float = -22.0,
    weighting: str = "number",
) -> Tuple[np.ndarray, EuclidSelection]:
    """
    Galaxy overdensity field built from the Euclid-selected halo catalogue.

    Thin wrapper over :func:`src.analysis.galaxy_overdensity_from_catalogue`
    that pins the grid to the geometry ``run_simulation.py`` uses for its
    catalogue-based ``galaxy_overdensity`` (§3b): ``HII_DIM`` transverse
    cells over ``BOX_LEN`` and ``N_z`` cells over ``BOX_LEN`` along the
    line of sight.  Computing it here — rather than reading the stored
    field — is what makes these figures "after the Euclid cut" even when
    the run used the default ``lightcone_sfr`` weighting, which applies no
    magnitude cut at all.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    M_UV_bright : float, optional
        Bright-end magnitude cut.  The faint end comes from the stored
        ``M_UV_limit`` attribute.
    weighting : {'number', 'luminosity'}, optional
        Per-halo weight, as in :func:`galaxy_overdensity_from_catalogue`.

    Returns
    -------
    delta_gal : ndarray
        Overdensity of shape ``(HII_DIM, HII_DIM, N_z)``.
    selection : EuclidSelection
        The magnitude window actually applied.

    Notes
    -----
    The line-of-sight axis spans the *coeval* box (``BOX_LEN``), not the
    lightcone extent ``L_los``: the halo catalogue is a coeval snapshot at
    ``z_obs``, and ``run_simulation.py`` makes the same choice
    (``los_extent=BOX_LEN``).  The two arrays therefore share a shape and a
    transverse plane but not an LOS scale — see ``PIPELINE.md``.

    When the catalogue was subsampled at load time (``--max-halos``), the
    field is the overdensity of that subsample: unbiased, but noisier per
    cell than the full catalogue.
    """
    M_UV_faint = float(data.get("M_UV_limit", -18.0))

    return galaxy_overdensity_from_catalogue(
        coords=_halo_coords_mpc(data),
        sfr=np.asarray(data.sfr, dtype=float),
        halo_masses=np.asarray(data.halo_masses, dtype=float),
        box_len=data.BOX_LEN,
        n_perp=data.HII_DIM,
        n_los=data.N_z,
        los_extent=data.BOX_LEN,
        weighting=weighting,
        M_UV_faint=M_UV_faint,
        M_UV_bright=M_UV_bright,
    )


def plot_euclid_selected_catalogue(
    data: SimulationData,
    M_UV_bright: float = -22.0,
    max_points: int = 50_000,
    seed: int = 42,
) -> Figure:
    """
    Galaxies, halo masses, and SFRs that survive the Euclid magnitude cut.

    The post-selection counterpart of :func:`plot_halo_catalogue` and
    :func:`plot_sfr_relations`, which show the full catalogue.  Three
    panels — the selected galaxies projected on the sky and coloured by
    ``M_UV``, the halo-mass distribution before and after the cut, and the
    SFR distribution before and after it with the equivalent SFR window
    marked.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    M_UV_bright : float, optional
        Bright-end magnitude cut.  The faint end comes from the stored
        ``M_UV_limit`` attribute.
    max_points : int, optional
        Cap on the number of galaxies scattered in the first panel.
    seed : int, optional
        RNG seed for the subsampling.

    Returns
    -------
    Figure
        The assembled 1 x 3 figure.
    """
    rng = np.random.default_rng(seed)

    M_UV_faint = float(data.get("M_UV_limit", -18.0))
    coords = _halo_coords_mpc(data)
    sfr = np.asarray(data.sfr, dtype=float)
    halo_masses = np.asarray(data.halo_masses, dtype=float)

    # The same call the bias stage and plot_uv_selection_maps make, so all
    # three figures describe one selection.
    selection = select_euclid_halos(
        sfr, halo_masses, M_UV_faint=M_UV_faint, M_UV_bright=M_UV_bright,
    )
    selected = _lift_selection_mask(sfr, halo_masses, selection)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # -- Panel 1: the selected galaxies on the sky ------------------------
    selected_idx = np.flatnonzero(selected)
    if selected_idx.size > max_points:
        selected_idx = rng.choice(selected_idx, size=max_points, replace=False)

    if selected_idx.size:
        magnitude = sfr_to_Muv(sfr[selected_idx])
        scatter = axes[0].scatter(
            coords[selected_idx, 0], coords[selected_idx, 1],
            s=np.clip((np.log10(halo_masses[selected_idx]) - 7.0) ** 2, 1, 30),
            c=magnitude,
            cmap="plasma_r",
            alpha=0.6,
            edgecolors="none",
        )
        colorbar = fig.colorbar(scatter, ax=axes[0], label=r"$M_{\rm UV}$")
        colorbar.ax.invert_yaxis()          # brighter (more negative) on top
    axes[0].set_xlabel("$x$  [Mpc]")
    axes[0].set_ylabel("$y$  [Mpc]")
    axes[0].set_xlim(0, data.BOX_LEN)
    axes[0].set_ylim(0, data.BOX_LEN)
    axes[0].set_aspect("equal")
    axes[0].set_title("Selected galaxies — projected positions")

    # -- Panel 2: halo mass, before and after the cut ---------------------
    all_masses = halo_masses[halo_masses > 0]
    mass_bins = np.linspace(
        np.log10(all_masses.min()), np.log10(all_masses.max()), 60,
    ) if all_masses.size else np.linspace(8.0, 12.0, 60)

    if all_masses.size:
        axes[1].hist(np.log10(all_masses), bins=mass_bins,
                     color="0.75", label="All halos")
    if selection.n_selected:
        axes[1].hist(np.log10(selection.halo_masses), bins=mass_bins,
                     color="crimson", alpha=0.8, label="Euclid-selected")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$\log_{10}(M_{\rm halo}/M_\odot)$")
    axes[1].set_ylabel("Number of halos")
    axes[1].legend(fontsize=9)
    axes[1].set_title("Halo mass — before and after the cut")

    # -- Panel 3: SFR, before and after the cut ---------------------------
    all_sfr = sfr[np.isfinite(sfr) & (sfr > 0)]
    sfr_bins = np.linspace(
        np.log10(all_sfr.min()), np.log10(all_sfr.max()), 80,
    ) if all_sfr.size else np.linspace(-4.0, 2.0, 80)

    if all_sfr.size:
        axes[2].hist(np.log10(all_sfr), bins=sfr_bins,
                     color="0.75", label="All halos")
    if selection.n_selected:
        axes[2].hist(np.log10(selection.sfr), bins=sfr_bins,
                     color="crimson", alpha=0.8, label="Euclid-selected")

    # The magnitude window is a pure SFR window under the Madau & Dickinson
    # calibration, so the cut lands on exactly these two SFR values.
    for bound, colour, label in (
        (selection.SFR_min, "darkorange", r"$M_{\rm UV}$ faint cut"),
        (selection.SFR_max, "navy", r"$M_{\rm UV}$ bright cut"),
    ):
        if bound > 0:
            axes[2].axvline(np.log10(bound), color=colour, ls="--", lw=1.5,
                            label=label)
    axes[2].set_yscale("log")
    axes[2].set_xlabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    axes[2].set_ylabel("Number of halos")
    axes[2].legend(fontsize=9)
    axes[2].set_title("SFR — before and after the cut")

    selected_fraction = (
        selection.n_selected / selection.n_valid if selection.n_valid else 0.0
    )
    fig.suptitle(
        rf"After the Euclid cut (${M_UV_bright} \leq M_{{\rm UV}} \leq {M_UV_faint}$) — "
        rf"{selection.n_selected:,} of {selection.n_valid:,} halos "
        rf"({selected_fraction:.2%}) at $z_{{\rm obs}} = {data.z_obs}$",
        fontsize=13,
    )
    return fig


def _overdensity_norm(field: np.ndarray) -> TwoSlopeNorm:
    """
    Diverging norm for an overdensity: floor at -1, long positive tail.

    ``delta_gal`` is bounded below by -1 (an empty cell) but unbounded
    above, so a symmetric scale saturates the whole map as soon as the
    catalogue is sparse.  Anchoring the midpoint at ``delta = 0`` and the
    top at a high percentile of the *occupied* cells keeps both the voids
    and the peaks readable.

    Parameters
    ----------
    field : ndarray
        Overdensity values.

    Returns
    -------
    TwoSlopeNorm
        Norm spanning ``[-1, vmax]`` with ``vcenter = 0``.
    """
    positive = field[field > 0]
    vmax = float(np.percentile(positive, 99)) if positive.size else 1.0
    return TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=max(vmax, 1e-3))


def _binned_mean_by(
    x: np.ndarray,
    y: np.ndarray,
    n_bins: int = 25,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Mean of ``y`` in equal-count bins of ``x``, with the standard error.

    Quantile edges rather than uniform ones, because ``delta T_b`` piles up
    at zero in the ionised cells; duplicate edges are collapsed, so a
    strongly degenerate field simply yields fewer bins.

    Parameters
    ----------
    x : ndarray
        Values to bin by.
    y : ndarray
        Values to average, same shape as ``x``.
    n_bins : int, optional
        Requested number of bins; the result may have fewer.

    Returns
    -------
    centres : ndarray
        Mean ``x`` per bin.
    means : ndarray
        Mean ``y`` per bin.
    errors : ndarray
        Standard error on ``means``.
    """
    edges = np.unique(np.quantile(x, np.linspace(0.0, 1.0, n_bins + 1)))
    if edges.size < 3:
        return (
            np.array([x.mean()]),
            np.array([y.mean()]),
            np.array([y.std() / max(np.sqrt(y.size), 1.0)]),
        )

    index = np.clip(np.digitize(x, edges[1:-1]), 0, edges.size - 2)
    counts = np.bincount(index, minlength=edges.size - 1).astype(float)
    occupied = counts > 0

    sum_x = np.bincount(index, weights=x, minlength=edges.size - 1)
    sum_y = np.bincount(index, weights=y, minlength=edges.size - 1)
    sum_yy = np.bincount(index, weights=y ** 2, minlength=edges.size - 1)

    counts = counts[occupied]
    centres = sum_x[occupied] / counts
    means = sum_y[occupied] / counts
    variance = np.maximum(sum_yy[occupied] / counts - means ** 2, 0.0)
    return centres, means, np.sqrt(variance / counts)


def _slab_mean(field: np.ndarray, centre: int, slab_cells: int) -> np.ndarray:
    """
    Mean of a 3D field over a window of line-of-sight cells.

    Parameters
    ----------
    field : ndarray
        Array of shape ``(n_x, n_y, n_los)``.
    centre : int
        Central LOS index.
    slab_cells : int
        Window width in cells; clipped to the array and to at least 1.

    Returns
    -------
    ndarray
        2D array of shape ``(n_x, n_y)``.
    """
    half = max(int(slab_cells), 1) // 2
    lo = max(centre - half, 0)
    hi = min(centre + max(int(slab_cells), 1) - half, field.shape[2])
    return field[:, :, lo:hi].mean(axis=2)


def plot_selected_galaxy_overdensity(
    data: SimulationData,
    M_UV_bright: float = -22.0,
    weighting: str = "number",
    delta_gal: Optional[np.ndarray] = None,
    selection: Optional[EuclidSelection] = None,
) -> Figure:
    """
    The galaxy overdensity field built from the Euclid-selected catalogue.

    Three panels — the line-of-sight projection of ``delta_gal``, a single
    transverse slice at mid-LOS, and the one-point distribution over every
    cell.

    The projection comes first because the selected sample is sparse: the
    Euclid window keeps a small fraction of the catalogue, so a single
    ``(HII_DIM, HII_DIM)`` slice holds far fewer galaxies than it has cells
    and is shot noise almost everywhere.  The projection is the same field
    integrated along the LOS, where the clustering is actually visible; the
    slice panel shows the raw, unprojected field the power-spectrum
    estimator is handed.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    M_UV_bright : float, optional
        Bright-end magnitude cut, used only when ``delta_gal`` is not given.
    weighting : {'number', 'luminosity'}, optional
        Per-halo weight, used only when ``delta_gal`` is not given.
    delta_gal : ndarray, optional
        Pre-computed field from :func:`selected_galaxy_overdensity`.  Pass it
        to avoid re-depositing a large catalogue for a second figure.
    selection : EuclidSelection, optional
        The selection that produced ``delta_gal``; required for the title
        when ``delta_gal`` is supplied.

    Returns
    -------
    Figure
        The assembled 1 x 3 figure.

    Notes
    -----
    The LOS axis spans the coeval box (``BOX_LEN``), not the lightcone
    ``L_los`` — see :func:`selected_galaxy_overdensity`.
    """
    if delta_gal is None or selection is None:
        delta_gal, selection = selected_galaxy_overdensity(
            data, M_UV_bright=M_UV_bright, weighting=weighting,
        )

    mid_z = delta_gal.shape[2] // 2
    projected = delta_gal.mean(axis=2)
    transverse = delta_gal[:, :, mid_z]

    extent_xy = [0, data.BOX_LEN, 0, data.BOX_LEN]

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    im0 = axes[0].imshow(
        projected.T, origin="lower", extent=extent_xy,
        cmap="RdBu_r", norm=_overdensity_norm(projected),
        interpolation="nearest",
    )
    fig.colorbar(im0, ax=axes[0], label=r"$\langle \delta_{\rm gal} \rangle_{\rm LOS}$")
    axes[0].set_title("LOS projection of the whole box")

    im1 = axes[1].imshow(
        transverse.T, origin="lower", extent=extent_xy,
        cmap="RdBu_r", norm=_overdensity_norm(delta_gal),
        interpolation="nearest",
    )
    fig.colorbar(im1, ax=axes[1], label=r"$\delta_{\rm gal}$")
    axes[1].set_title(f"Single transverse slice (LOS cell {mid_z})")

    for ax in axes[:2]:
        ax.set_xlabel("$x$  [Mpc]")
        ax.set_ylabel("$y$  [Mpc]")
        ax.set_aspect("equal")

    # -- Panel 3: the one-point distribution ------------------------------
    # A sparse count field piles up at delta = -1 and has a tail thousands
    # of times wider, so the axis has to be symlog to show both.
    flat = delta_gal.ravel()
    n_cells = flat.size
    empty_fraction = float(np.mean(flat <= -1.0 + 1e-12))
    occupancy = selection.n_selected / n_cells if n_cells else 0.0

    span = max(float(np.max(np.abs(flat))), 1.0)
    positive_bins = np.logspace(0, np.log10(span), 40)
    bins = np.concatenate(([-1.0], np.linspace(-1.0, 1.0, 21)[1:], positive_bins[1:]))

    axes[2].hist(flat, bins=np.unique(bins), color="steelblue", alpha=0.85)
    axes[2].axvline(0.0, color="k", ls="--", lw=1.2, label=r"$\delta_{\rm gal} = 0$")
    axes[2].set_xscale("symlog", linthresh=1.0)
    axes[2].set_yscale("log")
    axes[2].set_xlabel(r"$\delta_{\rm gal}$")
    axes[2].set_ylabel("Number of cells")
    axes[2].legend(fontsize=9)
    axes[2].set_title(
        f"One-point distribution\n"
        f"empty cells {empty_fraction:.1%},  "
        f"mean occupancy {occupancy:.3g} galaxies/cell"
    )

    fig.suptitle(
        rf"Galaxy overdensity after the Euclid cut "
        rf"(${selection.M_UV_bright} \leq M_{{\rm UV}} \leq {selection.M_UV_faint}$, "
        rf"{weighting}-weighted) — {selection.n_selected:,} galaxies on a "
        rf"{delta_gal.shape[0]}$\times${delta_gal.shape[1]}$\times${delta_gal.shape[2]} grid",
        fontsize=13,
    )
    return fig


def plot_galaxy_overdensity_on_21cm(
    data: SimulationData,
    M_UV_bright: float = -22.0,
    weighting: str = "number",
    slab_cells: int = 8,
    smooth_cells: float = 2.0,
    delta_gal: Optional[np.ndarray] = None,
    selection: Optional[EuclidSelection] = None,
) -> Figure:
    """
    The post-cut galaxy overdensity overlaid on the 21 cm field.

    Three panels — the 21 cm brightness temperature with ``delta_gal``
    contours on top, the same view with the two roles swapped, and the
    cell-by-cell joint distribution of the two fields with their Pearson
    correlation coefficient.

    Both maps are the mean over the *same* window of line-of-sight cells, so
    the two fields are paired exactly as ``compute_all_power_spectra`` pairs
    them; the sign of ``r`` in the third panel is therefore the real-space
    counterpart of the large-scale cross-power.  Averaging over a slab
    rather than a single cell is what makes the sparse selected sample
    legible at all.

    The overlays are transverse only: the halo catalogue is a coeval
    snapshot and the 21 cm field is a lightcone, so the two share the
    ``(x, y)`` plane and the grid, but not a line-of-sight scale — see
    :func:`selected_galaxy_overdensity` and ``PIPELINE.md``.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation with a non-empty halo catalogue.
    M_UV_bright : float, optional
        Bright-end magnitude cut, used only when ``delta_gal`` is not given.
    weighting : {'number', 'luminosity'}, optional
        Per-halo weight, used only when ``delta_gal`` is not given.
    slab_cells : int, optional
        Line-of-sight cells averaged into each map, centred on mid-LOS.
    smooth_cells : float, optional
        Gaussian smoothing width, in cells, applied to the contoured field.
        Display only — the shot noise of a sparse selected catalogue makes
        raw per-cell contours unreadable.  ``0`` disables it.
    delta_gal : ndarray, optional
        Pre-computed field from :func:`selected_galaxy_overdensity`.
    selection : EuclidSelection, optional
        The selection that produced ``delta_gal``.

    Returns
    -------
    Figure
        The assembled 1 x 3 figure.
    """
    if delta_gal is None or selection is None:
        delta_gal, selection = selected_galaxy_overdensity(
            data, M_UV_bright=M_UV_bright, weighting=weighting,
        )

    temperature = np.asarray(data.brightness_temp_field, dtype=float)
    n_los = min(delta_gal.shape[2], temperature.shape[2])
    mid_z = n_los // 2

    delta_map = _slab_mean(delta_gal, mid_z, slab_cells)
    temp_map = _slab_mean(temperature, mid_z, slab_cells)

    def _smooth(field: np.ndarray) -> np.ndarray:
        return gaussian_filter(field, smooth_cells) if smooth_cells > 0 else field

    delta_smooth = _smooth(delta_map)
    temp_smooth = _smooth(temp_map)

    extent_xy = [0, data.BOX_LEN, 0, data.BOX_LEN]
    grid_x = np.linspace(0, data.BOX_LEN, delta_map.shape[0])
    grid_y = np.linspace(0, data.BOX_LEN, delta_map.shape[1])

    vmax_t = float(np.percentile(temp_map, 99.5))

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # -- Panel 1: delta_gal contours over the 21 cm field -----------------
    im0 = axes[0].imshow(
        temp_map.T, origin="lower", extent=extent_xy,
        cmap=eor_colormap(), vmin=0.0, vmax=max(vmax_t, 1e-3),
        interpolation="nearest",
    )
    fig.colorbar(im0, ax=axes[0], label=r"$\delta T_b$  [mK]")

    gal_levels = _display_contour_levels(delta_smooth)
    if gal_levels.size:
        contour = axes[0].contour(
            grid_x, grid_y, delta_smooth.T, levels=gal_levels,
            colors=CONTOUR_COLOR_ON_EOR, linewidths=1.1, alpha=1.0,
        )
        _outline(contour)
        labels = axes[0].clabel(contour, fmt=lambda v: f"{v:.2f}", fontsize=7)
        for label in labels:
            label.set_path_effects(
                [patheffects.withStroke(linewidth=2.0, foreground="0.15")]
            )
    axes[0].set_title(r"$\delta_{\rm gal}$ contours over $\delta T_b$")

    # -- Panel 2: the same maps, roles swapped ----------------------------
    im1 = axes[1].imshow(
        delta_map.T, origin="lower", extent=extent_xy,
        cmap="RdBu_r", norm=_overdensity_norm(delta_map),
        interpolation="nearest",
    )
    fig.colorbar(im1, ax=axes[1], label=r"$\delta_{\rm gal}$")

    temp_levels = _display_contour_levels(temp_smooth)
    if temp_levels.size:
        contour = axes[1].contour(
            grid_x, grid_y, temp_smooth.T, levels=temp_levels,
            colors=CONTOUR_COLOR_ON_DIVERGING, linewidths=1.1, alpha=1.0,
        )
        _outline(contour)
        labels = axes[1].clabel(contour, fmt=lambda v: f"{v:.0f}", fontsize=7)
        for label in labels:
            label.set_path_effects(
                [patheffects.withStroke(linewidth=2.0, foreground="0.15")]
            )
    axes[1].set_title(r"$\delta T_b$ contours over $\delta_{\rm gal}$")

    for ax in axes[:2]:
        ax.set_xlabel("$x$  [Mpc]")
        ax.set_ylabel("$y$  [Mpc]")
        ax.set_aspect("equal")

    # -- Panel 3: the cell-by-cell relation over every cell ---------------
    # A 2D histogram is unreadable here: delta_gal is a sparse count field,
    # so it is quantised and ~97% of cells sit at -1.  Binning delta_gal by
    # delta T_b instead averages that shot noise away and leaves the trend
    # the cross-power measures.
    gal_flat = delta_gal[:, :, :n_los].ravel()
    temp_flat = temperature[:, :, :n_los].ravel()

    finite = np.isfinite(gal_flat) & np.isfinite(temp_flat)
    gal_flat, temp_flat = gal_flat[finite], temp_flat[finite]

    usable = (
        gal_flat.size > 0
        and np.ptp(gal_flat) > 0
        and np.ptp(temp_flat) > 0
    )
    if usable:
        centres, means, errors = _binned_mean_by(temp_flat, gal_flat)
        axes[2].errorbar(centres, means, yerr=errors, fmt="o-", ms=4, lw=1.2,
                         color="crimson", ecolor="0.5", capsize=2)
        pearson_r = float(np.corrcoef(temp_flat, gal_flat)[0, 1])
        correlation_label = rf"$r = {pearson_r:+.3f}$"
    else:
        correlation_label = "degenerate field"
    axes[2].axhline(0.0, color="k", ls="--", lw=1.0)
    axes[2].set_xlabel(r"$\delta T_b$  [mK]")
    axes[2].set_ylabel(r"$\langle \delta_{\rm gal} \rangle$")
    axes[2].set_title(
        rf"Mean $\delta_{{\rm gal}}$ per $\delta T_b$ bin  ({correlation_label})"
    )

    fig.suptitle(
        rf"Euclid-selected galaxies against the 21 cm field at "
        rf"$z_{{\rm obs}} = {data.z_obs}$  "
        rf"(${selection.M_UV_bright} \leq M_{{\rm UV}} \leq {selection.M_UV_faint}$, "
        rf"{selection.n_selected:,} galaxies; maps are the mean over "
        rf"{slab_cells} LOS cells at cell {mid_z})",
        fontsize=13,
    )
    return fig


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


def _signed_log_norm(
    power: np.ndarray,
) -> Tuple[np.ndarray, SymLogNorm]:
    """
    Masked cross-power and a symmetric log norm that preserves its sign.

    Replaces an earlier ``sign(P) * log10|P|`` transform that had three
    defects, all of them visible in the cross-power panel:

    * **Sign inversion below unity.**  ``log10|P|`` is negative for
      ``|P| < 1``, so ``sign(P) * log10|P|`` *flipped the sign* of every bin
      with ``|P| < 1``.  On a diverging colormap positive cross-power in that
      range rendered blue and negative rendered red — the colours actively
      misreported the sign wherever the signal was small.
    * **A blind spot at unity.**  ``|P| = 1`` maps to exactly 0, the centre of
      a diverging map, so those bins rendered blank white and read as missing
      data.
    * **Empty bins conflated with signal.**  ``NaN`` and exactly-zero bins were
      also mapped to 0.0, i.e. the same blank white, and the subsequent
      ``fill_nan_nearest`` call was dead code because the ``np.where`` had
      already removed every ``NaN``.

    :class:`~matplotlib.colors.SymLogNorm` handles all three: it is linear
    within ``linthresh`` of zero and logarithmic outside it, monotonic in the
    signed value throughout, so the colour is always an honest function of
    ``P``.  Empty bins are masked instead of coloured, so "no data" is
    visually distinct from "small".

    Parameters
    ----------
    power : ndarray
        Power spectrum, already transposed for display.

    Returns
    -------
    masked : ndarray
        ``power`` with non-finite bins masked.
    norm : SymLogNorm
        Symmetric-log norm covering the finite data.

    Notes
    -----
    ``linthresh`` is the 1st percentile of ``|P|``, floored at ``1e-8 x`` the
    colour limit.  It has to track the data floor rather than the peak: the
    cross-power spans roughly seven decades, and a linthresh set as a fixed
    fraction of the maximum would put most bins inside the linear region,
    washing them out to the neutral centre — the same visual failure this
    function exists to remove.
    """
    masked = np.ma.masked_invalid(np.asarray(power, dtype=float))
    finite = masked.compressed()
    nonzero = np.abs(finite[finite != 0])

    if nonzero.size == 0:
        return masked, SymLogNorm(linthresh=1.0, vmin=-1.0, vmax=1.0, base=10)

    clim = float(np.percentile(nonzero, 99))
    # The cross-power spans ~7 decades, so linthresh must sit near the FLOOR
    # of the data, not a fixed fraction of the peak -- otherwise most of the
    # range falls inside the linear region and washes out to the neutral
    # centre, reproducing the very "missing bins" look this replaced.
    linthresh = max(float(np.percentile(nonzero, 1)), clim * 1e-8)
    return masked, SymLogNorm(
        linthresh=linthresh, vmin=-clim, vmax=clim, base=10
    )


def _mathtext_float(value: float, significant: int = 3) -> str:
    """
    Render a float as mathtext, using ``a × 10^b`` for extreme exponents.

    Python's ``%g`` produces ``1.06e-111``, which mathtext typesets with a
    stray gap around the minus sign.  Splitting the exponent out avoids it.

    Parameters
    ----------
    value : float
        Number to render.
    significant : int, optional
        Significant figures in the mantissa.

    Returns
    -------
    str
        A mathtext fragment, already wrapped in ``$…$`` where needed.
    """
    formatted = f"{value:.{significant}g}"
    if "e" not in formatted:
        return formatted
    mantissa, exponent = formatted.split("e")
    return rf"${mantissa} \times 10^{{{int(exponent)}}}$"


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

    cross, cross_norm = _signed_log_norm(spectra.P_cross.T)
    cmap_cross = plt.get_cmap("RdBu_r").copy()
    cmap_cross.set_bad("0.85")          # empty bins: grey, not white
    im2 = axes[2].pcolormesh(k_perp, k_parallel, cross, cmap=cmap_cross,
                             shading="auto", norm=cross_norm)
    fig.colorbar(im2, ax=axes[2],
                 label=r"$P_{21 \times \rm gal}$  [mK Mpc³], symlog")
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


def plot_galaxy_wedge(
    spectra: PowerSpectra,
    data: SimulationData,
    horizon_slope: float,
    fov_slope: float,
) -> Figure:
    """
    The galaxy auto-power spectrum with the wedge region filled, not outlined.

    Companion to :func:`plot_power_spectra`'s middle panel.  Drawing the wedge
    as two lines leaves the eye to decide which side is excluded; shading and
    hatching the contaminated region instead makes the accessible window
    immediately legible.

    Parameters
    ----------
    spectra : PowerSpectra
        Provides ``P_galaxy_auto`` and the ``(k_perp, k_parallel)`` grid.
    data : SimulationData
        Loaded simulation, for the title metadata.
    horizon_slope, fov_slope : float
        Wedge slopes from ``src.analysis``.

    Returns
    -------
    Figure
        A single-panel figure.
    """
    k_perp = spectra.k_perp
    k_parallel = spectra.k_parallel

    fig, ax = plt.subplots(figsize=(7.5, 6))

    with np.errstate(divide="ignore", invalid="ignore"):
        pgal_display = fill_nan_nearest(np.log10(np.abs(spectra.P_galaxy_auto.T)))

    im = ax.pcolormesh(k_perp, k_parallel, pgal_display, cmap="plasma",
                       shading="auto", zorder=1)
    fig.colorbar(im, ax=ax, label=r"$\log_{10} |P_{\rm gal}|$  [Mpc³]")

    # Wedge region: everything below the horizon line.  A dark overlay rather
    # than a light one — on "plasma" a white wash brightens the already-bright
    # low-k corner and the legend swatch disappears.
    line = np.logspace(np.log10(k_perp.min()), np.log10(k_perp.max()), 200)
    horizon_line = line * horizon_slope
    ax.fill_between(
        line, k_parallel[0], horizon_line,
        where=horizon_line > k_parallel[0],
        facecolor="black", alpha=0.35, hatch="///", edgecolor="white", lw=0.0,
        zorder=2, label="Wedge (excluded)",
    )
    _add_wedge_lines(ax, k_perp, horizon_slope, fov_slope, "w", label=True)

    _style_k_axes(ax, k_perp, k_parallel)
    ax.set_title(
        rf"Galaxy power spectrum with foreground wedge, "
        rf"$z_{{\rm obs}} = {data.z_obs}$"
    )
    ax.legend(loc="upper left", fontsize=9, framealpha=0.85)
    return fig


def plot_wedge_real_space(
    data: SimulationData,
    horizon_slope: float,
    wedge_buffer: float = 0.0,
) -> Figure:
    """
    What foreground-wedge excision does to the galaxy field in real space.

    The wedge is a statement about Fourier modes, but discarding them changes
    the field itself.  This FFTs the 3D galaxy overdensity, zeroes every mode
    inside the wedge, transforms back, and shows the same line-of-sight slice
    before and after on a shared colour scale.

    The surviving modes satisfy ``k_par > slope * k_perp`` with a slope of
    order 3 at z ~ 7, i.e. they vary rapidly along the line of sight and
    slowly across it — so the filtered field appears striped across the LOS
    axis, with most of the original clumping gone.

    Parameters
    ----------
    data : SimulationData
        Supplies ``galaxy_overdensity`` and the box geometry.
    horizon_slope : float
        Wedge slope from :func:`src.analysis.horizon_wedge_slope`.
    wedge_buffer : float, optional
        Safety margin added above the wedge line [Mpc^-1].  Default 0.0, the
        bare horizon boundary ``k_par <= slope * k_perp``.  The uncertainty
        budget applies a non-zero buffer, so it discards slightly more than
        this figure shows.

    Returns
    -------
    Figure
        The assembled 1 × 2 figure.
    """
    delta_gal = np.asarray(data.galaxy_overdensity, dtype=float)
    n_x, n_y, n_z = delta_gal.shape

    # Wavenumber grids for the non-cubic (N, N, N_z) lightcone box.
    kx = np.fft.fftfreq(n_x, d=data.cell_size) * 2 * np.pi
    ky = np.fft.fftfreq(n_y, d=data.cell_size) * 2 * np.pi
    kz = np.fft.fftfreq(n_z, d=data.L_los / n_z) * 2 * np.pi

    # The wedge condition factorises: k_perp varies over the transverse plane
    # only, k_par over the LOS axis only.  So the (n_x*n_y, n_z) mask that
    # foreground_wedge_mask returns for the flattened transverse plane
    # reshapes straight back onto the 3D grid — no separate implementation of
    # the boundary is needed here.
    k_perp_plane = np.hypot(kx[:, np.newaxis], ky[np.newaxis, :])
    outside_wedge = foreground_wedge_mask(
        k_perp_plane.ravel(), np.abs(kz), horizon_slope, buffer=wedge_buffer,
    ).reshape(n_x, n_y, n_z)

    percent_excluded = 100.0 * (1.0 - outside_wedge.sum() / outside_wedge.size)

    delta_filtered = np.fft.ifftn(
        np.where(outside_wedge, np.fft.fftn(delta_gal), 0.0)
    ).real

    mid_y = n_y // 2
    extent_los = [data.lc_dist_Mpc[0], data.lc_dist_Mpc[-1], 0, data.BOX_LEN]
    vmax = float(np.percentile(np.abs(delta_gal), 99))

    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
    panels = (
        (delta_gal, r"Original  $\delta_{\rm gal}$"),
        (delta_filtered,
         rf"Wedge-filtered  ({percent_excluded:.1f}% of modes zeroed)"),
    )
    for ax, (field, title) in zip(axes, panels):
        im = ax.imshow(
            field[:, mid_y, :],
            origin="lower", extent=extent_los, aspect="auto",
            cmap="RdBu_r", vmin=-vmax, vmax=vmax,
        )
        ax.set_title(title)
        ax.set_xlabel("Comoving distance  [Mpc]")

    axes[0].set_ylabel(r"Transverse $x$  [Mpc]")
    fig.colorbar(im, ax=axes, label=r"$\delta_{\rm gal}$",
                 fraction=0.025, pad=0.02)
    fig.suptitle(
        rf"Foreground wedge in real space, $z_{{\rm obs}} = {data.z_obs}$   "
        rf"({percent_excluded:.1f}% of modes excluded)",
        fontsize=13,
    )
    return fig


def plot_photoz_suppression(
    budget: UncertaintyBudget,
    data: SimulationData,
    sigma_z_scenarios: Optional[Sequence[float]] = None,
) -> Figure:
    """
    The photo-z damping kernel swept over a range of survey scenarios.

    :func:`plot_uncertainty_budget`'s first panel shows ``W(k_par)`` for the
    single adopted ``sigma_z``.  This one sweeps it, from the spectroscopic
    limit (``sigma_z = 0``, ``sigma_r = 0``, ``W = 1`` everywhere) up to the
    adopted value, so the cost of photometric redshifts is visible as a family
    rather than a point.

    ``sigma_r`` for each scenario is obtained by scaling the budget's own
    ``radial_smearing`` — ``sigma_r = c sigma_z / H(z)`` is linear in
    ``sigma_z``, so this is exact and cannot drift from the adopted value.

    Parameters
    ----------
    budget : UncertaintyBudget
        Supplies ``k_parallel``, the adopted ``photoz_uncertainty`` and its
        ``radial_smearing``.
    data : SimulationData
        Loaded simulation, for the title metadata.
    sigma_z_scenarios : sequence of float, optional
        Absolute photo-z uncertainties to draw.  The adopted value is always
        appended if absent.  Default ``(0, 0.02, 0.05, 0.10, 0.30)``.

    Returns
    -------
    Figure
        A single-panel figure.
    """
    k_parallel = budget.k_parallel
    sigma_z_adopted = float(budget.photoz_uncertainty)

    scenarios = list(sigma_z_scenarios if sigma_z_scenarios is not None
                     else (0.0, 0.02, 0.05, 0.10, 0.30))
    if not any(np.isclose(s, sigma_z_adopted) for s in scenarios):
        scenarios.append(sigma_z_adopted)
    scenarios.sort()

    # sigma_r = c sigma_z / H(z) is linear in sigma_z, so scaling the budget's
    # own value reproduces radial_smearing_length exactly for every scenario.
    if sigma_z_adopted > 0:
        sigma_r_per_sigma_z = budget.radial_smearing / sigma_z_adopted
    else:  # degenerate configuration; nothing to scale
        sigma_r_per_sigma_z = 0.0

    colors = plt.cm.viridis_r(np.linspace(0.10, 0.95, len(scenarios)))

    # Extend below the lowest sampled bin only when 1/sigma_r falls outside it,
    # so the damping-scale marker stays on the axis.
    sigma_r_adopted = float(budget.radial_smearing)
    k_damp = 1.0 / sigma_r_adopted if sigma_r_adopted > 0 else k_parallel[0]
    k_lo = min(k_parallel[0], k_damp) * 0.7
    k_line = np.logspace(np.log10(k_lo), np.log10(k_parallel[-1]), 400)

    fig, ax = plt.subplots(figsize=(7.5, 5))

    for sigma_z, color in zip(scenarios, colors):
        sigma_r = sigma_z * sigma_r_per_sigma_z
        kernel = photoz_damping_kernel(k_line, sigma_r).ravel()

        adopted = np.isclose(sigma_z, sigma_z_adopted)
        label = rf"$\sigma_z = {sigma_z:g}$"
        label += "  (spectroscopic)" if sigma_z == 0.0 else \
                 rf",  $\sigma_r = {sigma_r:.0f}$ Mpc"
        if adopted:
            label += "  — adopted"

        ax.plot(k_line, kernel, color=color, lw=2.4 if adopted else 1.5,
                zorder=3 if adopted else 2, label=label)

    if k_damp < k_parallel[0]:
        ax.axvspan(k_lo, k_parallel[0], color="0.5", alpha=0.12, lw=0, zorder=0,
                   label=r"below the box $k_\parallel$ range")
    ax.axvline(k_damp, color="k", ls="--", lw=1.2, alpha=0.8, zorder=4)
    ax.text(
        k_damp * 1.45, 0.26,
        rf"$1/\sigma_r = {k_damp:.3f}$ Mpc$^{{-1}}$" "\n"
        rf"($\sigma_z = {sigma_z_adopted:g}$)",
        fontsize=9, va="top", ha="left", zorder=5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85),
    )

    ax.set_xscale("log")
    ax.set_xlim(k_lo, k_parallel[-1])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"$k_\parallel$  [Mpc$^{-1}$]")
    ax.set_ylabel(r"$W(k_\parallel)$")
    ax.set_title(rf"Photo-$z$ suppression, $z_{{\rm obs}} = {data.z_obs}$")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    return fig


def plot_power_spectra_1d(
    spectra: PowerSpectra,
    data: SimulationData,
    shot_noise: float,
    thermal_noise: Optional[float] = None,
    outside_wedge: Optional[np.ndarray] = None,
    n_bins: int = 12,
) -> Figure:
    """
    Spherically averaged 21 cm, galaxy and cross power, with the noise floors.

    The cylindrical spectra are the ones the SNR is built from; this collapses
    them onto ``|k|`` so the three can be read against each other and against
    the shot-noise and thermal floors on one axis.

    The galaxy shot noise ``1/n`` is drawn as a horizontal line on the galaxy
    panel — the scale at which the measured galaxy auto-power stops being
    signal and becomes counting noise.

    Parameters
    ----------
    spectra : PowerSpectra
        Cylindrical spectra.
    data : SimulationData
        Loaded simulation, for the title metadata.
    shot_noise : float
        Galaxy shot noise ``P_N,gal = 1/n`` [Mpc^3].
    thermal_noise : float, optional
        21 cm thermal noise [mK^2 Mpc^3].  A scalar under the flat model; the
        mean of a ``k_perp``-resolved array under the physical model.
    outside_wedge : ndarray, optional
        Wedge mask.  When given, a second dashed curve shows the
        foreground-clean average alongside the all-mode one.
    n_bins : int, optional
        Number of spherical bins.

    Returns
    -------
    Figure
        Three panels: 21 cm auto, galaxy auto, and cross power.

    Notes
    -----
    The 21 cm signal is anisotropic, so a spherical average is a summary for
    reading off amplitudes — not a substitute for the cylindrical spectra.
    """
    everything = spherically_average_spectra(spectra, n_bins=n_bins)
    clean = (
        None if outside_wedge is None
        else spherically_average_spectra(
            spectra, n_bins=n_bins, outside_wedge=outside_wedge
        )
    )

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    panels = (
        ("P_21cm_auto", r"$P_{21}(k)$  [mK$^2$ Mpc$^3$]",
         r"21 cm auto-power", "crimson"),
        ("P_galaxy_auto", r"$P_{\rm gal}(k)$  [Mpc$^3$]",
         r"Galaxy auto-power", "royalblue"),
        ("P_cross", r"$|P_{\times}(k)|$  [mK Mpc$^3$]",
         "21 cm $\\times$ galaxy  (negative on large scales; "
         r"$|P_\times|$ shown)", "darkviolet"),
    )

    for ax, (name, ylabel, title, colour) in zip(axes, panels):
        for source, style, tag in (
            (everything, "-", "all modes"),
            (clean, "--", "outside wedge"),
        ):
            if source is None:
                continue
            values = np.abs(getattr(source, name))
            good = np.isfinite(values) & (values > 0)
            if not good.any():
                continue
            ax.plot(source.k[good], values[good], style, color=colour,
                    lw=2.0 if style == "-" else 1.6,
                    alpha=1.0 if style == "-" else 0.75,
                    marker="o" if style == "-" else None, ms=4,
                    label=tag)

        if name == "P_galaxy_auto" and np.isfinite(shot_noise):
            ax.axhline(shot_noise, color="darkorange", ls=":", lw=2.0,
                       label=rf"$1/\bar n = {shot_noise:.3g}$ Mpc$^3$")
        if name == "P_21cm_auto" and thermal_noise is not None and (
                np.isfinite(thermal_noise)):
            ax.axhline(thermal_noise, color="dimgray", ls=":", lw=2.0,
                       label=rf"$P_{{N,21}} = {thermal_noise:.3g}$")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$k$  [Mpc$^{-1}$]", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=12)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=9, framealpha=0.9)

    fig.suptitle(
        rf"Spherically averaged power spectra  ($z = {data.z_obs:.1f}$)",
        fontsize=13,
    )
    fig.tight_layout()
    return fig


def plot_number_density(
    density: "NumberDensity",
    data: SimulationData,
) -> Figure:
    """
    Cumulative comoving number density and the shot noise it implies.

    Left: ``n(< M_UV)`` against the magnitude cut, with the adopted selection
    window and the configured ``mean_galaxy_density`` marked.  Right: the
    same as ``1/n``, the shot-noise power that enters the cross-power
    variance.

    Parameters
    ----------
    density : NumberDensity
        From :func:`src.analysis.comoving_number_density`.
    data : SimulationData
        Loaded simulation, for the title metadata.

    Returns
    -------
    Figure
        Two panels sharing the magnitude axis.
    """
    fig, (ax_n, ax_shot) = plt.subplots(1, 2, figsize=(13, 5.2))

    good = density.n_cumulative > 0
    ax_n.plot(density.M_UV[good], density.n_cumulative[good],
              color="royalblue", lw=2.2, label=r"$n(< M_{\rm UV})$, catalogue")
    ax_n.axhline(density.adopted_mean_density, color="crimson", ls="--", lw=1.8,
                 label=rf"adopted $\bar n = "
                       rf"{density.adopted_mean_density:.3g}$ Mpc$^{{-3}}$")
    if density.n_at_selection > 0:
        ax_n.axhline(density.n_at_selection, color="seagreen", ls="-.", lw=1.8,
                     label=rf"in window $= {density.n_at_selection:.3g}$ "
                           rf"Mpc$^{{-3}}$")
    ax_n.set_ylabel(r"$n(< M_{\rm UV})$  [Mpc$^{-3}$]", fontsize=12)
    ax_n.set_title("Cumulative comoving number density", fontsize=12)

    ax_shot.plot(density.M_UV[good], density.shot_noise[good],
                 color="darkorange", lw=2.2, label=r"$1/n(< M_{\rm UV})$")
    if density.adopted_mean_density > 0:
        ax_shot.axhline(1.0 / density.adopted_mean_density, color="crimson",
                        ls="--", lw=1.8,
                        label=rf"adopted $1/\bar n = "
                              rf"{1.0 / density.adopted_mean_density:.4g}$ "
                              rf"Mpc$^3$")
    ax_shot.set_ylabel(r"$1/n$  [Mpc$^3$]", fontsize=12)
    ax_shot.set_title(r"Implied shot noise $P_{N,\rm gal}$", fontsize=12)

    m_faint = float(data.get("M_UV_limit", -18.0))
    m_bright = float(data.get("M_UV_bright", -22.0))
    for ax in (ax_n, ax_shot):
        ax.axvspan(m_bright, m_faint, color="mediumseagreen", alpha=0.13,
                   zorder=0, label="Selection window")
        ax.axvline(m_faint, color="dimgray", ls="--", lw=1.3)
        ax.axvline(m_bright, color="darkslateblue", ls=":", lw=1.5)
        ax.set_yscale("log")
        ax.set_xlabel(r"$M_{\rm UV}$ cut  [AB mag]", fontsize=12)
        ax.invert_xaxis()
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=9, framealpha=0.9, loc="best")

    fig.suptitle(
        rf"Number density and shot noise  ($z = {data.z_obs:.1f}$, "
        rf"$V = {density.volume_Mpc3:.2e}$ Mpc$^3$)",
        fontsize=13,
    )
    fig.tight_layout()
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

    cross, cross_norm = _signed_log_norm(P_cross_observed.T)
    cmap_cross = plt.get_cmap("RdBu_r").copy()
    cmap_cross.set_bad("0.85")
    im_cross = axes[1].pcolormesh(k_perp, k_parallel, cross, cmap=cmap_cross,
                                  shading="auto", norm=cross_norm)
    fig.colorbar(im_cross, ax=axes[1],
                 label=r"$P_{21 \times \rm gal}^{\rm obs}$  [mK Mpc³], symlog")
    _add_wedge_lines(axes[1], k_perp, horizon_slope, fov_slope, "k")
    axes[1].set_title(r"Observed $P_{21 \times \rm gal}$  (photo-$z$ damped)")

    for ax in axes:
        _style_k_axes(ax, k_perp, k_parallel)

    fig.suptitle(
        rf"HERA $\times$ Euclid — lightcone $z = {data.z_min}$–${data.z_max}$   "
        rf"($\sigma_z = {data.get('photoz_uncertainty', 0.45)}$)",
        fontsize=13,
    )
    return fig


def plot_uncertainty_budget(
    budget: UncertaintyBudget,
    data: SimulationData,
) -> Figure:
    """
    Where the cross-power uncertainty comes from, in three panels.

    Panel 1 — the photo-z damping kernel ``W(k_par)`` against the smallest
    ``k_par`` the wedge admits, which is what makes the two cuts conflict.
    Panel 2 — ``σ_cross`` over the ``(k_perp, k_parallel)`` plane, with the
    wedge boundary overlaid.
    Panel 3 — the share of ``σ_cross²`` carried by sample variance rather than
    noise coupling, i.e. which term limits each mode.

    Parameters
    ----------
    budget : UncertaintyBudget
        Output of ``src.analysis.compute_uncertainty_budget``.
    data : SimulationData
        Loaded simulation, for the title metadata.

    Returns
    -------
    Figure
        The assembled 1 × 3 figure.
    """
    k_perp = budget.k_perp
    k_parallel = budget.k_parallel

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))

    # ── Panel 1: damping kernel vs the wedge's lowest admitted k_par ──────
    kernel = budget.photoz_kernel.ravel()
    axes[0].loglog(k_parallel, np.clip(kernel, 1e-300, None), "o-", lw=2,
                   color="C0", label=rf"$W$, $\sigma_z={budget.photoz_uncertainty:g}$")
    # Smallest k_par the wedge lets through, at the smallest k_perp.
    k_par_wedge_min = k_perp.min() * budget.horizon_slope + budget.wedge_buffer
    axes[0].axvline(k_par_wedge_min, color="crimson", ls="--", lw=1.5,
                    label=r"wedge floor at $k_\perp^{\rm min}$")
    axes[0].axhline(0.5, color="grey", ls=":", lw=1.2, label="$W = 0.5$")
    axes[0].set_xlabel(r"$k_\parallel$  [Mpc$^{-1}$]")
    axes[0].set_ylabel(r"$W(k_\parallel)$")
    axes[0].set_ylim(max(np.min(kernel[kernel > 0], initial=1e-12) * 0.1, 1e-30), 2.0)
    axes[0].legend(fontsize=8, loc="lower left")
    axes[0].set_title(rf"Photo-$z$ damping ($\sigma_r = {budget.radial_smearing:.1f}$ Mpc)")

    # ── Panel 2: sigma_cross across the plane ─────────────────────────────
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_display = fill_nan_nearest(np.log10(budget.snr.sigma_cross.T))
    im_sigma = axes[1].pcolormesh(k_perp, k_parallel, sigma_display,
                                  cmap="viridis", shading="auto")
    fig.colorbar(im_sigma, ax=axes[1], label=r"$\log_{10}\sigma_\times$")
    _add_wedge_lines(axes[1], k_perp, budget.horizon_slope, budget.fov_slope,
                     "w", label=True)
    axes[1].legend(loc="upper left", fontsize=8)
    axes[1].set_title(r"Cross-power uncertainty $\sigma_\times$")
    _style_k_axes(axes[1], k_perp, k_parallel)

    # ── Panel 3: which term dominates the variance ────────────────────────
    with np.errstate(divide="ignore", invalid="ignore"):
        denominator = (
            budget.snr.cosmic_variance_term + budget.snr.noise_coupling_term
        )
        fraction = np.where(
            denominator > 0,
            budget.snr.cosmic_variance_term / denominator,
            np.nan,
        )
    im_frac = axes[2].pcolormesh(k_perp, k_parallel, fill_nan_nearest(fraction.T),
                                 cmap="coolwarm", shading="auto", vmin=0, vmax=1)
    fig.colorbar(im_frac, ax=axes[2],
                 label=r"cosmic-variance share of $\sigma_\times^2$")
    _add_wedge_lines(axes[2], k_perp, budget.horizon_slope, budget.fov_slope, "k")
    axes[2].set_title("Sample variance (1) vs noise (0)")
    _style_k_axes(axes[2], k_perp, k_parallel)

    fig.suptitle(
        rf"Uncertainty budget — $z_{{\rm obs}} = {budget.z_obs}$, "
        rf"$\sigma_z = {budget.photoz_uncertainty:g}$, "
        rf"buffer $= {budget.wedge_buffer:g}$ Mpc$^{{-1}}$   "
        rf"({budget.fraction_outside_wedge:.1%} of modes usable, "
        rf"total SNR = {_mathtext_float(budget.total_snr)} $\sigma$)",
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
