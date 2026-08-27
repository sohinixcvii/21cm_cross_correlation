#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
analysis.py — Post-simulation analysis for the 21 cm × galaxy cross-correlation
================================================================================

Pure computation used by ``run_pipeline.py`` (and previously inlined in
``notebooks/analysis.ipynb``):

1. 2D cylindrical power spectra of the non-cubic lightcone box
2. Foreground-wedge geometry (horizon and HERA primary-beam lines)
3. Euclid photometric-redshift damping along the line of sight
4. Simplified HERA thermal noise, galaxy shot noise, and the per-mode SNR
5. The assembled uncertainty budget, :func:`compute_uncertainty_budget` —
   the single entry point that chains 2–4 in the order the notebook does
6. Euclid ``M_UV`` selection of the halo catalogue and the resulting
   effective galaxy bias

Nothing here imports ``matplotlib`` or ``py21cmfast``; plotting lives in
``src/figures.py`` and the simulation itself in ``run_simulation.py``.

References
----------
La Plante et al. (2023) — arXiv:2205.09770 (Eqs. 15–17, wedge geometry)
Davies, Mesinger & Murray (2025) — arXiv:2504.17254
Thyagarajan et al. (2015), ApJ 804, 14 (foreground wedge)
Sheth & Tormen (1999), MNRAS 308, 119 (halo bias)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.integrate import quad

try:  # local package import (repo root on sys.path)
    from src.dataio import PowerSpectra
    from src.conversions import (
        Luv_to_Muv, Luv_to_sfr, Muv_to_Luv, sfr_to_Luv, sheth_tormen_bias,
    )
except ImportError:  # direct import of the module (src/ on sys.path)
    from dataio import PowerSpectra
    from conversions import (
        Luv_to_Muv, Luv_to_sfr, Muv_to_Luv, sfr_to_Luv, sheth_tormen_bias,
    )

__all__ = [
    "hubble_parameter",
    "comoving_distance",
    "star_formation_timescale",
    "stellar_mass_to_sfr",
    "T_STAR_DEFAULT",
    "compute_cylindrical_cross_power",
    "compute_all_power_spectra",
    "ESTIMATOR_MODES",
    "MEAN_SUBTRACTION_MODES",
    "subtract_field_mean",
    "blackman_harris_taper",
    "subband_index_ranges",
    "compute_subband_power_spectra",
    "combine_band_snr",
    "SubbandGeometry",
    "horizon_wedge_slope",
    "fov_wedge_slope",
    "foreground_wedge_mask",
    "photoz_damping_kernel",
    "radial_smearing_length",
    "system_temperature",
    "hera_thermal_noise_power",
    "cross_power_snr",
    "total_snr",
    "compute_uncertainty_budget",
    "UncertaintyBudget",
    "T_RECEIVER_K",
    "T_SKY_300MHZ_K",
    "SKY_SPECTRAL_INDEX",
    "NOISE_NORMALISATION_MPC3",
    "HERA_ANTENNA_SPACING_M",
    "HERA_HEX_N_SIDE",
    "HERA_OMEGA_P_OVER_PP",
    "HERA_APERTURE_EFFICIENCY",
    "N_POLARISATIONS",
    "NOISE_EQUIVALENT_BANDWIDTH",
    "cosmological_scalar_x2y",
    "hera_beam_solid_angles",
    "hera_baseline_counts",
    "hera_thermal_noise_power_physical",
    "euclid_sfr_window",
    "select_euclid_halos",
    "effective_galaxy_bias",
    "deposit_halo_field",
    "galaxy_overdensity_from_catalogue",
    "GALAXY_WEIGHTING_MODES",
    "EuclidSelection",
    "BiasEstimate",
    "SNRResult",
]


# ===========================================================================
# 1  Background cosmology
# ===========================================================================

def hubble_parameter(
    z: float | np.ndarray,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
) -> float | np.ndarray:
    """
    Hubble parameter H(z) for a flat ΛCDM cosmology.

    Parameters
    ----------
    z : float or ndarray
        Redshift.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].  Default: Planck 2018.
    omega_m : float, optional
        Present-day matter density parameter.

    Returns
    -------
    float or ndarray
        H(z) [km s^-1 Mpc^-1].
    """
    return hubble_constant * np.sqrt(omega_m * (1.0 + z) ** 3 + (1.0 - omega_m))


def comoving_distance(
    z: float,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
) -> float:
    """
    Line-of-sight comoving distance to redshift ``z``.

    Parameters
    ----------
    z : float
        Target redshift.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].
    omega_m : float, optional
        Present-day matter density parameter.
    speed_of_light_kms : float, optional
        Speed of light [km s^-1].

    Returns
    -------
    float
        D_c(z) [Mpc].
    """
    integral, _ = quad(
        lambda z_: speed_of_light_kms / hubble_parameter(z_, hubble_constant, omega_m),
        0.0,
        z,
    )
    return float(integral)


# 21cmFAST "simple" template default: the star-formation timescale is a fixed
# fraction of the Hubble time (Park et al. 2019; scaling_relations.c).
T_STAR_DEFAULT = 0.5

# Megaparsec in kilometres, for converting 1/H(z) from Mpc s km^-1 to seconds.
_KM_PER_MPC = 3.0856775814913673e19
_SEC_PER_YR = 365.25 * 24 * 3600


def star_formation_timescale(
    z: float,
    t_star: float = T_STAR_DEFAULT,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
) -> float:
    """
    21cmFAST star-formation timescale ``t_sf = t_STAR × t_H(z)``.

    This is the timescale 21cmFAST actually uses internally to turn a stellar
    mass into a star-formation rate (``sfr = M_star / (t_STAR t_H)``), so any
    external SFR model must use it too — a hardcoded 100 Myr instead of the
    ~570 Myr this returns at z = 7 biases every derived UV luminosity by
    ~1.9 magnitudes.

    Parameters
    ----------
    z : float
        Redshift.
    t_star : float, optional
        Fraction of the Hubble time, 21cmFAST's ``t_STAR``.  Default 0.5,
        the "simple" template value.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].
    omega_m : float, optional
        Present-day matter density parameter.

    Returns
    -------
    float
        t_sf [yr].

    References
    ----------
    Park et al. (2019), MNRAS 484, 933 — Eq. 3.
    ``py21cmfast/src/scaling_relations.c``, ``get_halo_sfr()``.
    """
    h_z_kms_mpc = hubble_parameter(z, hubble_constant, omega_m)
    t_hubble_seconds = _KM_PER_MPC / h_z_kms_mpc
    return float(t_star * t_hubble_seconds / _SEC_PER_YR)


def stellar_mass_to_sfr(
    stellar_mass: float | np.ndarray,
    z: float,
    t_star: float = T_STAR_DEFAULT,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
) -> float | np.ndarray:
    """
    Convert stellar mass to star-formation rate using 21cmFAST's own timescale.

    Parameters
    ----------
    stellar_mass : float or ndarray
        Stellar mass [M_sun].
    z : float
        Redshift.
    t_star : float, optional
        21cmFAST ``t_STAR``.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].
    omega_m : float, optional
        Present-day matter density parameter.

    Returns
    -------
    float or ndarray
        SFR [M_sun yr^-1].
    """
    return stellar_mass / star_formation_timescale(
        z, t_star, hubble_constant, omega_m
    )


# ===========================================================================
# 2  Cylindrical power spectra
# ===========================================================================

def compute_cylindrical_cross_power(
    field_a: np.ndarray,
    field_b: np.ndarray,
    box_len_perp: float,
    box_len_los: float,
    n_bins_perp: int = 20,
    n_bins_parallel: int = 20,
    taper: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    2D cylindrical cross-power spectrum P(k_perp, k_parallel) of a lightcone box.

    Handles the non-cubic ``(N, N, N_z)`` geometry by using separate cell sizes
    and fundamental modes for the transverse and line-of-sight directions.

    Parameters
    ----------
    field_a, field_b : ndarray
        Real-space fields of shape ``(N, N, N_z)``.  Pass the same array twice
        for an auto-spectrum.
    box_len_perp : float
        Comoving transverse side length [Mpc].
    box_len_los : float
        Comoving line-of-sight length [Mpc].
    n_bins_perp, n_bins_parallel : int, optional
        Number of log-spaced bins along k_perp and k_parallel.
    taper : ndarray, optional
        Line-of-sight window of length ``N_z``, applied to both fields before
        the transform and divided back out of the power as ``<w^2>`` so the
        amplitude is preserved.  ``None`` (the default) is the bare FFT, which
        is what every result before the sub-band estimator used.  See
        :func:`blackman_harris_taper`.

    Returns
    -------
    k_perp_centres : ndarray
        Log bin centres along k_perp [Mpc^-1], shape ``(n_bins_perp,)``.
    k_parallel_centres : ndarray
        Log bin centres along k_parallel [Mpc^-1], shape ``(n_bins_parallel,)``.
    power_2d : ndarray
        Binned power, shape ``(n_bins_perp, n_bins_parallel)``.  Empty bins
        are ``NaN``.
    mode_counts : ndarray
        Number of Fourier modes per bin, same shape as ``power_2d``.

    Raises
    ------
    ValueError
        If the two fields have different shapes.
    """
    if field_a.shape != field_b.shape:
        raise ValueError(
            f"field shapes differ: {field_a.shape} vs {field_b.shape}"
        )

    n_perp = field_a.shape[0]   # transverse cells
    n_los = field_a.shape[2]    # LOS cells

    # ── Optional line-of-sight taper ──────────────────────────────────────
    # Applied to both fields, so the measured power carries a factor <w^2>
    # that is divided back out below.  This is the noise-equivalent-bandwidth
    # normalisation: it restores the amplitude of a statistically homogeneous
    # field exactly, and leaves the bare-FFT path bit-identical when
    # ``taper is None``.
    taper_normalisation = 1.0
    if taper is not None:
        taper = np.asarray(taper, dtype=float)
        if taper.shape != (n_los,):
            raise ValueError(
                f"taper length {taper.shape} does not match the "
                f"line-of-sight axis ({n_los},)"
            )
        taper_normalisation = float(np.mean(taper ** 2))
        if taper_normalisation <= 0.0:
            raise ValueError("taper has zero power; cannot normalise")
        window = taper[np.newaxis, np.newaxis, :]
        field_a = field_a * window
        field_b = field_b * window

    dx = box_len_perp / n_perp   # transverse cell size [Mpc]
    dz = box_len_los / n_los     # LOS cell size [Mpc]
    volume = box_len_perp ** 2 * box_len_los

    # ── Fourier transforms ────────────────────────────────────────────────
    ft_factor = dx * dx * dz
    fourier_a = np.fft.fftn(field_a) * ft_factor
    fourier_b = np.fft.fftn(field_b) * ft_factor
    power_3d = (fourier_a * np.conj(fourier_b)).real / (volume * taper_normalisation)

    # ── Wavenumber grids ──────────────────────────────────────────────────
    kx = np.fft.fftfreq(n_perp, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(n_perp, d=dx) * 2 * np.pi
    kz = np.fft.fftfreq(n_los, d=dz) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

    k_perp_3d = np.sqrt(KX ** 2 + KY ** 2)
    k_parallel_3d = np.abs(KZ)

    # ── Bin edges (log-spaced, starting below the first discrete mode) ────
    dk_perp = 2 * np.pi / box_len_perp
    dk_par = 2 * np.pi / box_len_los

    k_max_perp = np.sqrt(2) * np.abs(kx).max() * 1.05
    k_max_par = np.abs(kz).max() * 1.05

    k_perp_edges = np.logspace(
        np.log10(0.5 * dk_perp), np.log10(k_max_perp), n_bins_perp + 1
    )
    k_par_edges = np.logspace(
        np.log10(0.5 * dk_par), np.log10(k_max_par), n_bins_parallel + 1
    )

    # ── Bin the 3D power into 2D ──────────────────────────────────────────
    power_2d = np.zeros((n_bins_perp, n_bins_parallel))
    mode_counts = np.zeros_like(power_2d)

    bin_perp = np.digitize(k_perp_3d.ravel(), k_perp_edges) - 1
    bin_par = np.digitize(k_parallel_3d.ravel(), k_par_edges) - 1
    power_flat = power_3d.ravel()

    inside = (
        (bin_perp >= 0) & (bin_perp < n_bins_perp)
        & (bin_par >= 0) & (bin_par < n_bins_parallel)
    )
    np.add.at(power_2d, (bin_perp[inside], bin_par[inside]), power_flat[inside])
    np.add.at(mode_counts, (bin_perp[inside], bin_par[inside]), 1)

    power_2d = np.divide(
        power_2d,
        mode_counts,
        where=mode_counts > 0,
        out=np.full_like(power_2d, np.nan),
    )

    k_perp_centres = np.sqrt(k_perp_edges[:-1] * k_perp_edges[1:])
    k_par_centres = np.sqrt(k_par_edges[:-1] * k_par_edges[1:])

    return k_perp_centres, k_par_centres, power_2d, mode_counts


def compute_all_power_spectra(
    brightness_temp_field: np.ndarray,
    galaxy_overdensity: np.ndarray,
    box_len_perp: float,
    box_len_los: float,
    n_bins_perp: int = 20,
    n_bins_parallel: int = 20,
    mean_subtraction: str = "global",
    taper: Optional[np.ndarray] = None,
) -> PowerSpectra:
    """
    Compute the 21 cm auto-, galaxy auto-, and cross-power spectra.

    The 21 cm field is converted to fluctuations by subtracting its mean; the
    galaxy field is already an overdensity.

    With ``mean_subtraction="per_slice"`` both fields have their per-slice
    transverse mean removed instead, which is what a lightcone spanning a
    non-negligible redshift range requires — see :func:`subtract_field_mean`.

    Parameters
    ----------
    brightness_temp_field : ndarray
        δT_b lightcone [mK], shape ``(N, N, N_z)``.
    galaxy_overdensity : ndarray
        Galaxy overdensity lightcone, same shape.
    box_len_perp, box_len_los : float
        Comoving transverse and line-of-sight box lengths [Mpc].
    n_bins_perp, n_bins_parallel : int, optional
        Number of log-spaced ``(k_perp, k_parallel)`` bins.
    mean_subtraction : {'global', 'per_slice'}, optional
        How the mean is removed before transforming.  ``'global'`` (default)
        subtracts a single scalar from the 21 cm field and leaves the galaxy
        overdensity untouched, reproducing every result this pipeline produced
        before the lightcone estimator existed.  ``'per_slice'`` subtracts the
        transverse mean of each line-of-sight slice from **both** fields.
    taper : ndarray, optional
        Line-of-sight window passed through to
        :func:`compute_cylindrical_cross_power`.

    Returns
    -------
    PowerSpectra
        The three spectra on a shared k-grid, plus the mode counts.
    """
    t21_fluctuations = subtract_field_mean(brightness_temp_field, mean_subtraction)
    galaxy_field = (
        galaxy_overdensity
        if mean_subtraction == "global"
        else subtract_field_mean(galaxy_overdensity, mean_subtraction)
    )

    k_perp, k_parallel, p_21, mode_counts = compute_cylindrical_cross_power(
        t21_fluctuations, t21_fluctuations,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel, taper,
    )
    _, _, p_gal, _ = compute_cylindrical_cross_power(
        galaxy_field, galaxy_field,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel, taper,
    )
    _, _, p_cross, _ = compute_cylindrical_cross_power(
        t21_fluctuations, galaxy_field,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel, taper,
    )

    return PowerSpectra(
        k_perp=k_perp,
        k_parallel=k_parallel,
        P_21cm_auto=p_21,
        P_galaxy_auto=p_gal,
        P_cross=p_cross,
        mode_counts=mode_counts,
    )


# ===========================================================================
# 2b  Lightcone estimator — per-slice means, LOS taper, sub-bands
# ===========================================================================
#
# The estimator above was inherited from a *coeval* notebook and assumes the
# box is statistically homogeneous along the line of sight.  That holds for a
# quasi-coeval slab and fails for a true lightcone, in four ways that
# `TODO.md` P0 enumerates:
#
#   P0.1  slices uniform in redshift are not uniform in comoving distance,
#         so the FFT mis-assigns k_parallel        -> run_simulation.py
#   P0.2  a single global mean leaves a monotonic LOS ramp that aliases into
#         low k_parallel                            -> subtract_field_mean
#   P0.3  one FFT over a wide band returns a redshift-averaged spectrum with
#         an ill-defined effective redshift         -> compute_subband_power_spectra
#   P0.4  the power spectrum and the noise are computed over different
#         bandwidths                                -> compute_subband_power_spectra
#
# Everything here is opt-in.  With the defaults (`mean_subtraction="global"`,
# `taper=None`, no sub-banding) the pipeline reproduces its previous numbers
# bit for bit; `ESTIMATOR_MODES` names the two formalisms.

#: The two estimator formalisms.  ``"coeval"`` is the historical one: global
#: mean subtraction, no taper, a single FFT over the whole box.  ``"lightcone"``
#: applies P0.2-P0.4 — per-slice means, a Blackman-Harris taper, and one
#: spectrum per sub-band at its own effective redshift.
ESTIMATOR_MODES: Tuple[str, ...] = ("coeval", "lightcone")

#: How the mean is removed before transforming.  See :func:`subtract_field_mean`.
MEAN_SUBTRACTION_MODES: Tuple[str, ...] = ("global", "per_slice")


def subtract_field_mean(field: np.ndarray, mode: str = "global") -> np.ndarray:
    """
    Remove the mean of a lightcone field.

    Parameters
    ----------
    field : ndarray
        Lightcone of shape ``(N, N, N_z)``; the line of sight is the last axis.
    mode : {'global', 'per_slice'}, optional
        ``'global'`` subtracts one scalar, the mean over the whole box.
        ``'per_slice'`` subtracts the transverse mean of each line-of-sight
        slice, ``field[:, :, i] - <field[:, :, i]>``.

    Returns
    -------
    ndarray
        Fluctuations about the chosen mean, same shape as ``field``.

    Raises
    ------
    ValueError
        If ``mode`` is not one of :data:`MEAN_SUBTRACTION_MODES`, or the field
        is not 3D.

    Notes
    -----
    This is `TODO.md` P0.2.  Over a lightcone the mean brightness temperature
    evolves strongly with redshift, so a single global scalar leaves a
    monotonic ramp along the line of sight, which is not signal.  Removing the
    per-slice mean removes it by construction, at the cost of also removing the
    genuine ``k_parallel = 0`` mode — which the wedge excises anyway.

    **How much this moves the binned spectra: measured, and less than
    `TODO.md` claims.**  A ramp that is uniform across the sky has ~99 % of its
    power at ``k_perp = 0``, and
    :func:`compute_cylindrical_cross_power` bins from ``0.5 dk_perp`` upwards,
    so that column is discarded before any binning happens.  On a field whose
    only contamination is such a ramp the two modes agree to floating-point
    precision (``tests/test_analysis.py``).  The operation is still correct,
    and it matters for anything that does use ``k_perp = 0`` — real-space
    fields, figures, and any future estimator that keeps the column — but it is
    not by itself the fix for low-``k_parallel`` contamination.

    What genuinely couples the redshift evolution into non-zero ``k_perp`` is
    that ``delta_T_b = T_0(z) [1 + delta(x)]``: an evolving mean *modulates*
    the fluctuation amplitude across the band.  Removing a per-slice mean does
    not undo a per-slice gain; per-slice **normalisation** would.  That is a
    separate change and is deliberately not made here.

    The galaxy overdensity is built with its own global normalisation
    (``SFR/<SFR> - 1``), so it carries the same ramp and needs the same
    treatment; :func:`compute_all_power_spectra` applies this to both fields
    when ``mean_subtraction="per_slice"``.

    Examples
    --------
    >>> box = np.ones((2, 2, 3)) * np.arange(3)      # a pure LOS ramp
    >>> np.allclose(subtract_field_mean(box, "per_slice"), 0.0)
    True
    """
    if mode not in MEAN_SUBTRACTION_MODES:
        raise ValueError(
            f"mode must be one of {MEAN_SUBTRACTION_MODES}, got {mode!r}"
        )
    field = np.asarray(field, dtype=float)
    if field.ndim != 3:
        raise ValueError(f"expected a 3D lightcone, got shape {field.shape}")

    if mode == "global":
        return field - field.mean()
    return field - field.mean(axis=(0, 1), keepdims=True)


def blackman_harris_taper(n_slices: int) -> np.ndarray:
    """
    Four-term Blackman-Harris window for the line-of-sight axis.

    Parameters
    ----------
    n_slices : int
        Length of the line-of-sight axis.  Must be at least 2.

    Returns
    -------
    ndarray
        The window, shape ``(n_slices,)``, peaking at 1.

    Raises
    ------
    ValueError
        If ``n_slices < 2``.

    Notes
    -----
    A band of finite extent is a top-hat in frequency, and a bare FFT of it
    leaks power from the bright low-``k_parallel`` modes across the whole
    axis — the ``k_parallel^-1.5`` leakage that ``src/foregrounds.py``
    measures.  The four-term Blackman-Harris window has -92 dB sidelobes and
    is the standard choice in 21 cm power-spectrum estimation for that reason.

    The window is *not* normalised here; the estimator divides the measured
    power by ``<w^2>`` so amplitudes are preserved
    (:func:`compute_cylindrical_cross_power`).  The cost is resolution: the
    effective number of independent line-of-sight modes falls, which is why
    tapering is applied per sub-band rather than to the whole box.

    References
    ----------
    Harris, F. J. (1978), Proc. IEEE 66, 51 — window functions for harmonic
    analysis.  Applied to 21 cm delay spectra by Parsons et al. (2012a).
    """
    if n_slices < 2:
        raise ValueError(f"n_slices must be >= 2, got {n_slices}")

    a = (0.35875, 0.48829, 0.14128, 0.01168)
    n = np.arange(n_slices)
    x = 2 * np.pi * n / (n_slices - 1)
    return a[0] - a[1] * np.cos(x) + a[2] * np.cos(2 * x) - a[3] * np.cos(3 * x)


@dataclass
class SubbandGeometry:
    """
    Where each sub-band sits along the lightcone.

    Attributes
    ----------
    index_ranges : ndarray
        ``(n_bands, 2)`` integer ``[start, stop)`` slice indices.
    z_effective : ndarray
        Effective redshift of each band, from its mean observed frequency.
    z_min, z_max : ndarray
        Redshift limits of each band.
    frequency_min_hz, frequency_max_hz : ndarray
        Observed frequency limits of each band [Hz].
    bandwidth_hz : ndarray
        Frequency span actually covered by each band [Hz].  This is the value
        the noise model must use (`TODO.md` P0.4).
    los_length_mpc : ndarray
        Comoving line-of-sight extent of each band [Mpc].
    n_slices : ndarray
        Number of lightcone slices in each band.
    """

    index_ranges: np.ndarray
    z_effective: np.ndarray
    z_min: np.ndarray
    z_max: np.ndarray
    frequency_min_hz: np.ndarray
    frequency_max_hz: np.ndarray
    bandwidth_hz: np.ndarray
    los_length_mpc: np.ndarray
    n_slices: np.ndarray

    @property
    def n_bands(self) -> int:
        """Number of sub-bands."""
        return int(self.index_ranges.shape[0])


def subband_index_ranges(
    n_slices: int,
    frequency_span_hz: float,
    bandwidth_hz: float = 8e6,
    min_slices_per_band: int = 8,
) -> np.ndarray:
    """
    Split a line-of-sight axis into contiguous sub-bands of at most one bandwidth.

    Parameters
    ----------
    n_slices : int
        Length of the line-of-sight axis.
    frequency_span_hz : float
        Total observed-frequency span of the lightcone [Hz].
    bandwidth_hz : float, optional
        Target per-band bandwidth [Hz] — the one the noise model assumes.
    min_slices_per_band : int, optional
        Refuse to split so finely that a band cannot support a useful
        ``k_parallel`` axis.  The band count is reduced until every band has
        at least this many slices.

    Returns
    -------
    ndarray
        ``(n_bands, 2)`` array of ``[start, stop)`` indices covering the axis
        with no gaps and no overlap.

    Raises
    ------
    ValueError
        If any argument is non-positive.

    Notes
    -----
    The number of bands is ``ceil(span / bandwidth)``, so each band spans
    ``span / n_bands <= bandwidth``.  Slices are divided as evenly as
    possible by index rather than by frequency: with comoving-uniform sampling
    the two differ by well under a cell, and an even split keeps every band's
    ``k_parallel`` grid identical, which is what makes the per-band spectra
    stackable.

    A lightcone narrower than one bandwidth returns a single band — the
    sub-band estimator then reduces to the whole-box one, with the taper still
    applied.
    """
    if n_slices < 1:
        raise ValueError(f"n_slices must be >= 1, got {n_slices}")
    if frequency_span_hz <= 0:
        raise ValueError(
            f"frequency_span_hz must be positive, got {frequency_span_hz}"
        )
    if bandwidth_hz <= 0:
        raise ValueError(f"bandwidth_hz must be positive, got {bandwidth_hz}")
    if min_slices_per_band < 1:
        raise ValueError(
            f"min_slices_per_band must be >= 1, got {min_slices_per_band}"
        )

    n_bands = int(np.ceil(frequency_span_hz / bandwidth_hz))
    n_bands = max(1, min(n_bands, n_slices // max(min_slices_per_band, 1)))
    n_bands = max(1, n_bands)

    edges = np.linspace(0, n_slices, n_bands + 1).astype(int)
    return np.column_stack((edges[:-1], edges[1:]))


def compute_subband_power_spectra(
    brightness_temp_field: np.ndarray,
    galaxy_overdensity: np.ndarray,
    lc_redshifts: np.ndarray,
    lc_dist_Mpc: np.ndarray,
    box_len_perp: float,
    bandwidth_hz: float = 8e6,
    n_bins_perp: int = 20,
    n_bins_parallel: int = 20,
    mean_subtraction: str = "per_slice",
    taper: bool = True,
    min_slices_per_band: int = 8,
    f_21_hz: float = 1420.405e6,
) -> Tuple[list, SubbandGeometry]:
    """
    Power spectra of a lightcone, one per frequency sub-band.

    Splits the line of sight into bands of at most ``bandwidth_hz``, removes
    the per-slice mean, applies a Blackman-Harris taper within each band, and
    measures the three spectra of that band on its own comoving extent.

    Parameters
    ----------
    brightness_temp_field, galaxy_overdensity : ndarray
        Lightcone fields of shape ``(N, N, N_z)``.
    lc_redshifts : ndarray
        Redshift of each slice, shape ``(N_z,)``.
    lc_dist_Mpc : ndarray
        Comoving distance to each slice [Mpc], shape ``(N_z,)``.
    box_len_perp : float
        Transverse box length [Mpc].
    bandwidth_hz : float, optional
        Target per-band bandwidth [Hz].  Match it to the noise model's.
    n_bins_perp, n_bins_parallel : int, optional
        Binning, applied identically in every band.
    mean_subtraction : {'per_slice', 'global'}, optional
        Passed to :func:`compute_all_power_spectra`.  ``'per_slice'`` is the
        point of this estimator; ``'global'`` is available for isolating the
        effect of one P0 item at a time.
    taper : bool, optional
        Apply the Blackman-Harris window along each band's line of sight.
    min_slices_per_band : int, optional
        Floor on the slices per band; see :func:`subband_index_ranges`.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].

    Returns
    -------
    list of PowerSpectra
        One entry per band, in lightcone order.
    SubbandGeometry
        Effective redshift, frequency limits, bandwidth and comoving extent of
        each band — everything the per-band uncertainty budget needs.

    Raises
    ------
    ValueError
        If the fields, redshifts and distances disagree in length, or the
        lightcone has fewer than two slices.

    Notes
    -----
    This is `TODO.md` P0.3 and P0.4 together.  One FFT over a wide lightcone
    returns a redshift-*averaged* spectrum whose effective redshift is
    ill-defined, and mixes a band of one width with a noise model of another.
    Measuring per band fixes both: each band has its own ``z_eff``, its own
    ``k_parallel`` grid from its own comoving extent, and its own bandwidth to
    hand to :func:`compute_uncertainty_budget`.

    The effective redshift is taken from the band's mean *observed frequency*,
    ``z_eff = f_21 / <nu> - 1``, rather than the mean redshift: the estimator
    works in frequency, and the two differ by ~0.2 % over an 8 MHz band at
    z = 7.

    Bands are combined afterwards with :func:`combine_band_snr`, not here —
    the spectra themselves are reported per band, since that is the only form
    in which they have a well-defined redshift.
    """
    brightness_temp_field = np.asarray(brightness_temp_field, dtype=float)
    galaxy_overdensity = np.asarray(galaxy_overdensity, dtype=float)
    lc_redshifts = np.asarray(lc_redshifts, dtype=float).ravel()
    lc_dist_Mpc = np.asarray(lc_dist_Mpc, dtype=float).ravel()

    if brightness_temp_field.shape != galaxy_overdensity.shape:
        raise ValueError(
            f"field shapes differ: {brightness_temp_field.shape} vs "
            f"{galaxy_overdensity.shape}"
        )
    n_los = brightness_temp_field.shape[2]
    if n_los < 2:
        raise ValueError(f"need at least 2 lightcone slices, got {n_los}")
    if lc_redshifts.size != n_los or lc_dist_Mpc.size != n_los:
        raise ValueError(
            f"lc_redshifts ({lc_redshifts.size}) and lc_dist_Mpc "
            f"({lc_dist_Mpc.size}) must match the LOS axis ({n_los})"
        )

    frequencies = f_21_hz / (1.0 + lc_redshifts)
    frequency_span = float(np.abs(frequencies[-1] - frequencies[0]))

    ranges = subband_index_ranges(
        n_slices=n_los,
        frequency_span_hz=max(frequency_span, np.finfo(float).tiny),
        bandwidth_hz=bandwidth_hz,
        min_slices_per_band=min_slices_per_band,
    )

    spectra_per_band = []
    z_eff, z_lo, z_hi = [], [], []
    nu_lo, nu_hi, band_width, los_length, slice_count = [], [], [], [], []

    for start, stop in ranges:
        n_band = int(stop - start)
        band_t21 = brightness_temp_field[:, :, start:stop]
        band_gal = galaxy_overdensity[:, :, start:stop]
        band_dist = lc_dist_Mpc[start:stop]
        band_nu = frequencies[start:stop]

        # Comoving extent: the slice spacing times the number of slices, so a
        # band of n cells spans n * dz_cell rather than the distance between
        # its first and last cell *centres*.
        if n_band > 1:
            cell = float(np.mean(np.diff(band_dist)))
        else:                                       # pragma: no cover - guarded
            cell = float(np.mean(np.diff(lc_dist_Mpc)))
        length = abs(cell) * n_band

        window = blackman_harris_taper(n_band) if taper else None

        spectra_per_band.append(
            compute_all_power_spectra(
                brightness_temp_field=band_t21,
                galaxy_overdensity=band_gal,
                box_len_perp=box_len_perp,
                box_len_los=length,
                n_bins_perp=n_bins_perp,
                n_bins_parallel=n_bins_parallel,
                mean_subtraction=mean_subtraction,
                taper=window,
            )
        )

        mean_nu = float(np.mean(band_nu))
        z_eff.append(f_21_hz / mean_nu - 1.0)
        band_z = lc_redshifts[start:stop]
        z_lo.append(float(band_z.min()))
        z_hi.append(float(band_z.max()))
        nu_lo.append(float(band_nu.min()))
        nu_hi.append(float(band_nu.max()))
        band_width.append(float(np.abs(band_nu.max() - band_nu.min())))
        los_length.append(length)
        slice_count.append(n_band)

    geometry = SubbandGeometry(
        index_ranges=np.asarray(ranges, dtype=int),
        z_effective=np.asarray(z_eff, dtype=float),
        z_min=np.asarray(z_lo, dtype=float),
        z_max=np.asarray(z_hi, dtype=float),
        frequency_min_hz=np.asarray(nu_lo, dtype=float),
        frequency_max_hz=np.asarray(nu_hi, dtype=float),
        bandwidth_hz=np.asarray(band_width, dtype=float),
        los_length_mpc=np.asarray(los_length, dtype=float),
        n_slices=np.asarray(slice_count, dtype=int),
    )
    return spectra_per_band, geometry


def combine_band_snr(band_totals: np.ndarray) -> float:
    """
    Combine per-band total SNRs into one detection significance.

    Parameters
    ----------
    band_totals : array_like
        Total SNR of each sub-band.  ``NaN`` entries are ignored.

    Returns
    -------
    float
        ``sqrt(sum(SNR_band^2))``.

    Notes
    -----
    Sub-bands sample disjoint frequency ranges and therefore disjoint comoving
    volumes, so their measurements are independent and add in quadrature —
    the same rule the per-bin sum inside a band already uses
    (:func:`total_snr`).
    """
    totals = np.asarray(band_totals, dtype=float)
    return float(np.sqrt(np.nansum(totals ** 2)))

# ===========================================================================
# 3  Foreground wedge geometry
# ===========================================================================

def horizon_wedge_slope(
    z_obs: float,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
) -> float:
    """
    Slope of the foreground horizon line in ``(k_perp, k_parallel)`` space.

    The horizon wedge boundary is ``k_par = m * k_perp`` with
    ``m = D_c(z) H(z) / [c (1 + z)]``.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].
    omega_m : float, optional
        Present-day matter density parameter.
    speed_of_light_kms : float, optional
        Speed of light [km s^-1].

    Returns
    -------
    float
        Dimensionless horizon slope ``m``.

    References
    ----------
    Thyagarajan et al. (2015); La Plante et al. (2023), Eq. 10.
    """
    d_c = comoving_distance(z_obs, hubble_constant, omega_m, speed_of_light_kms)
    h_z = hubble_parameter(z_obs, hubble_constant, omega_m)
    return float(d_c * h_z / (speed_of_light_kms * (1.0 + z_obs)))


def fov_wedge_slope(
    z_obs: float,
    dish_diameter: float = 14.0,
    f_21_hz: float = 1420.405e6,
    speed_of_light_mps: float = 3e8,
    **kwargs: float,
) -> float:
    """
    Slope of the primary-beam ("FoV") wedge line for a dish of given diameter.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    dish_diameter : float, optional
        Dish diameter [m].  Default: HERA, 14 m.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].
    speed_of_light_mps : float, optional
        Speed of light [m s^-1].
    **kwargs
        Passed through to :func:`horizon_wedge_slope` (``hubble_constant``,
        ``omega_m``, ``speed_of_light_kms``).

    Returns
    -------
    float
        Dimensionless FoV wedge slope ``sin(θ_FoV) × m_horizon``.
    """
    lambda_obs = speed_of_light_mps * (1.0 + z_obs) / f_21_hz   # [m]
    theta_fov = lambda_obs / dish_diameter                       # [rad]
    return float(np.sin(theta_fov) * horizon_wedge_slope(z_obs, **kwargs))


def foreground_wedge_mask(
    k_perp: np.ndarray,
    k_parallel: np.ndarray,
    slope: float,
    buffer: float = 0.0677,
) -> np.ndarray:
    """
    Boolean mask of ``(k_perp, k_parallel)`` bins that survive wedge excision.

    Parameters
    ----------
    k_perp, k_parallel : ndarray
        Bin centres [Mpc^-1].
    slope : float
        Wedge slope, from :func:`horizon_wedge_slope`.
    buffer : float, optional
        Safety margin added above the wedge line [Mpc^-1]. Default 0.0677
        = 0.1 h Mpc^-1 at h = 0.6766 (Pober et al. 2014 "moderate" foreground
        model; the 21cmSense ``horizon_buffer`` default).

    Returns
    -------
    ndarray of bool
        Shape ``(len(k_perp), len(k_parallel))``; True outside the wedge.
    """
    k_perp_grid, k_par_grid = np.meshgrid(k_perp, k_parallel, indexing="ij")
    return k_par_grid > (k_perp_grid * slope + buffer)


# ===========================================================================
# 4  Photometric-redshift damping
# ===========================================================================

def radial_smearing_length(
    photoz_uncertainty: float,
    z_obs: float,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
) -> float:
    """
    Comoving line-of-sight smearing induced by photometric-redshift errors.

    ``σ_r = c σ_z / H(z_obs)``.

    Parameters
    ----------
    photoz_uncertainty : float
        Photometric redshift uncertainty σ_z, **absolute — not σ_z/(1+z)**.
        Surveys usually quote the fractional form, so multiply by (1 + z)
        before passing it here: Euclid's σ_z/(1+z) < 0.05 requirement means
        σ_z ≈ 0.4 at z = 7, not 0.05.
    z_obs : float
        Reference redshift.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1].
    omega_m : float, optional
        Present-day matter density parameter.
    speed_of_light_kms : float, optional
        Speed of light [km s^-1].

    Returns
    -------
    float
        σ_r [Mpc].
    """
    h_z = hubble_parameter(z_obs, hubble_constant, omega_m)
    return float(speed_of_light_kms * photoz_uncertainty / h_z)


def photoz_damping_kernel(
    k_parallel: np.ndarray,
    radial_smearing: float,
) -> np.ndarray:
    """
    Gaussian photo-z damping kernel ``W(k_par) = exp(-k_par² σ_r² / 2)``.

    The galaxy auto-spectrum is damped by ``W²`` and the cross-spectrum by
    ``W`` (only one field is smeared).

    Parameters
    ----------
    k_parallel : ndarray
        Line-of-sight bin centres [Mpc^-1], shape ``(n_par,)``.
    radial_smearing : float
        σ_r [Mpc], from :func:`radial_smearing_length`.

    Returns
    -------
    ndarray
        Kernel of shape ``(1, n_par)``, broadcastable against a
        ``(n_perp, n_par)`` power spectrum.
    """
    k_par_broadcast = np.asarray(k_parallel)[np.newaxis, :]
    return np.exp(-0.5 * k_par_broadcast ** 2 * radial_smearing ** 2)


# ===========================================================================
# 5  Noise and signal-to-noise
# ===========================================================================

@dataclass
class SNRResult:
    """
    Per-mode and cumulative cross-correlation signal-to-noise.

    Attributes
    ----------
    snr_per_mode : ndarray
        ``|P_cross| / σ_cross`` for every ``(k_perp, k_parallel)`` bin.
    snr_outside_wedge : ndarray
        Same, with wedge-contaminated bins set to ``NaN``.
    total_snr : float
        Quadrature sum over the bins outside the wedge [σ].
    sigma_cross : ndarray
        Per-mode cross-power uncertainty.
    P_noise_21cm : float or ndarray
        21 cm thermal noise power [mK^2 Mpc^3].  An array when the
        ``k_perp``-resolved instrument model was used.
    P_noise_galaxy : float
        Galaxy shot noise power [Mpc^3].
    mode_weight : ndarray or None
        ``sqrt(N_patch dN)`` applied per bin (La Plante Eq. 19), or ``None``
        when the unweighted per-bin ratio was returned.
    sigma_21cm : ndarray
        ``|P_21| + P_N,21`` — the 21 cm side of the La Plante Eq. 15 product.
    sigma_galaxy : ndarray
        ``|P_gal| + P_N,gal`` — the galaxy side of the same product.
    cosmic_variance_term : ndarray
        ``0.5 P_cross²``, the sample-variance half of ``σ_cross²``.
    noise_coupling_term : ndarray
        ``0.5 σ_21 σ_gal``, the noise half of ``σ_cross²``.
    """

    snr_per_mode: np.ndarray
    snr_outside_wedge: np.ndarray
    total_snr: float
    sigma_cross: np.ndarray
    P_noise_21cm: float | np.ndarray
    P_noise_galaxy: float
    mode_weight: Optional[np.ndarray] = None
    sigma_21cm: Optional[np.ndarray] = None
    sigma_galaxy: Optional[np.ndarray] = None
    cosmic_variance_term: Optional[np.ndarray] = None
    noise_coupling_term: Optional[np.ndarray] = None


# ── Sky + receiver temperature model ───────────────────────────────────────
# T_sys = T_rcvr + T_sky(ν), with the synchrotron sky scaled from 300 MHz.
# These are the values the notebook inlines in its noise cell; naming them
# keeps the pipeline and the notebook auditable against each other.
# This is the 21cmSense convention, not DeBoer et al. (2017) Table 2 - the two
# differ, and the difference matters. 21cmSense's calc_sense.py computes
#     Tsky = 60e3 * (3e8 / freq_Hz) ** 2.55        [mK]
# i.e. 60 K x (lambda/1 m)^2.55, identical to the form below because
# 300 MHz / nu = lambda / 1 m. DeBoer et al. (2017) Table 2 instead give
# T_sys = 100 + 120 (nu/150 MHz)^-2.55 K, whose sky term is 120 K at 150 MHz
# against this model's 60 x 2^2.55 = 352 K - a factor 2.9 colder. The
# 21cmSense values are kept: they are what the notebook inlined, and they are
# the more conservative (hotter) choice. See NUMBERS_AND_SOURCES.md section 3.
T_RECEIVER_K = 100.0        # HERA receiver temperature [K] - 21cmSense `Trx`
T_SKY_300MHZ_K = 60.0       # Galactic synchrotron sky at 300 MHz [K] - 21cmSense
SKY_SPECTRAL_INDEX = 2.55   # T_sky ∝ ν^-2.55 - Pober et al. (2013, 2014)

# The notebook's noise expression is T_sys² × 1e3 / (t_int Δν).  T_sys² / (t Δν)
# is dimensionless × mK², so this factor carries the [Mpc^3] that makes P_N
# commensurate with the measured P_21 — i.e. it stands in for the survey
# volume per mode that a full instrument model (X²YΩ'/n(k_perp)) would supply.
# It is a normalisation of the scaling estimate, not a physical constant.
# See `hera_thermal_noise_power_physical` for the instrument model it stands in
# for, and NUMBERS_AND_SOURCES.md §3 for the ~10^4 discrepancy between them.
NOISE_NORMALISATION_MPC3 = 1e3

# ── HERA instrument model (for the physical noise path) ────────────────────
# Values traceable to the references named on each line.  See
# NUMBERS_AND_SOURCES.md §3 for the full audit.

#: Antenna-to-antenna spacing of the hexagonal core [m].
#: DeBoer et al. (2017), PASP 129, 045001 — "14.6 m center-to-center spacing".
HERA_ANTENNA_SPACING_M = 14.6

#: Antennas per side of the close-packed hexagonal core.  11 per side gives
#: 3n² − 3n + 1 = 331 elements, the closest hex number to HERA's 320-element
#: split core (350 total = 320 core + 30 outriggers; DeBoer et al. 2017).
#: The outriggers contribute negligibly to the power-spectrum sensitivity and
#: are not modelled.
HERA_HEX_N_SIDE = 11

#: Ratio Ω_P / Ω_PP of the beam and squared-beam solid angles, median over the
#: CST-simulated HERA beam models across 100–200 MHz.
#: Parsons (2017), "Power Spectrum Normalizations for HERA", HERA memo —
#: prints "HERA Omega_P/OMEGA_PP 2.1752891255".  (PAPER's is 2.35.)
HERA_OMEGA_P_OVER_PP = 2.175

#: Aperture efficiency of a HERA dish.  Calibrated so that the textbook
#: relation Ω_P = λ²/A_e with A_e = η_ap π D²/4 reproduces the Ω_P ≈ 0.04 sr
#: at 150 MHz plotted in Parsons (2017) for D = 14 m: A_e = λ²/Ω_P = 100 m²
#: against a geometric area of 154 m².
HERA_APERTURE_EFFICIENCY = 0.65

#: Orthogonal polarisations combined to measure the unpolarised signal.
#: Parsons (2017) Eq. 12 — "the factor of N_pol in the denominator explicitly
#: counts the two orthogonal polarizations".
N_POLARISATIONS = 2

#: Noise-equivalent bandwidth of the line-of-sight taper.  Exactly 1 for no
#: taper, which is what `compute_cylindrical_cross_power` applies.
#: Parsons (2017) Eq. 12, ``WINDOW = 'none'``.
NOISE_EQUIVALENT_BANDWIDTH = 1.0


def cosmological_scalar_x2y(
    z_obs: float,
    f_21_hz: float = 1420.405e6,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
) -> float:
    """
    The ``X²Y`` scalar converting (sr, Hz) to comoving volume.

    ``X = D_c(z)`` maps angle to transverse comoving distance and
    ``Y = c (1+z)² / [H(z) f_21]`` maps frequency to line-of-sight comoving
    distance, so ``X²Y`` carries a visibility-space power spectrum into
    Mpc³.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology.

    Returns
    -------
    float
        ``X²Y`` [Mpc³ sr⁻¹ Hz⁻¹].  ≈ 1227 at z = 7 for Planck 2018.

    References
    ----------
    Parsons et al. (2012a), ApJ 756, 165 — Eq. 12 and the surrounding
    definitions of ``X`` and ``Y``.
    Parsons (2017), "Power Spectrum Normalizations for HERA" — Eq. 1, where
    ``X²Y`` is "a cosmological scalar with units of h⁻³ Mpc³ / (sr·Hz)".
    """
    distance = comoving_distance(
        z_obs, hubble_constant, omega_m, speed_of_light_kms
    )
    hubble_z = hubble_parameter(z_obs, hubble_constant, omega_m)
    y_factor = speed_of_light_kms * (1.0 + z_obs) ** 2 / (hubble_z * f_21_hz)
    return float(distance ** 2 * y_factor)


def hera_beam_solid_angles(
    z_obs: float,
    dish_diameter: float = 14.0,
    aperture_efficiency: float = HERA_APERTURE_EFFICIENCY,
    omega_p_over_pp: float = HERA_OMEGA_P_OVER_PP,
    f_21_hz: float = 1420.405e6,
    speed_of_light_mps: float = 3e8,
) -> Tuple[float, float]:
    """
    HERA primary-beam solid angle and the effective beam area for power spectra.

    ``Ω_P = λ²/A_e`` with ``A_e = η_ap π D²/4`` is the standard antenna
    relation.  The quantity a power-spectrum normalisation needs is not
    ``Ω_P`` but ``Ω_eff ≡ Ω_P²/Ω_PP``, where ``Ω_PP`` is the solid angle of
    the *squared* beam — using ``Ω_P`` alone is a known error in
    power-spectrum normalisation.

    Parameters
    ----------
    z_obs : float
        Reference redshift, which sets the observed wavelength.
    dish_diameter : float, optional
        Dish diameter [m].  Default: HERA, 14 m (DeBoer et al. 2017).
    aperture_efficiency : float, optional
        η_ap.  See :data:`HERA_APERTURE_EFFICIENCY`.
    omega_p_over_pp : float, optional
        Ω_P/Ω_PP.  See :data:`HERA_OMEGA_P_OVER_PP`.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].
    speed_of_light_mps : float, optional
        Speed of light [m s⁻¹].

    Returns
    -------
    omega_p : float
        Primary-beam solid angle [sr].  ≈ 0.04 at 150 MHz for a 14 m dish.
    omega_eff : float
        ``Ω_P²/Ω_PP`` [sr], the factor entering the noise power.

    References
    ----------
    Parsons et al. (2014), ApJ 788, 106 — Appendix B, the definition of
    ``Ω_eff = Ω_P²/Ω_PP`` and why the squared beam is the relevant one.
    Parsons (2017), "Power Spectrum Normalizations for HERA" — Eqs. 2, 3, 7;
    HERA ``Ω_P/Ω_PP = 2.175``, and ``Ω_P ≈ 0.04`` sr at 150 MHz.
    Thompson, Moran & Swenson (2017), *Interferometry and Synthesis in Radio
    Astronomy*, 3rd ed. — the ``Ω_P A_e = λ²`` antenna theorem.
    """
    wavelength = speed_of_light_mps * (1.0 + z_obs) / f_21_hz      # [m]
    effective_area = aperture_efficiency * np.pi * dish_diameter ** 2 / 4.0
    omega_p = wavelength ** 2 / effective_area
    return float(omega_p), float(omega_p * omega_p_over_pp)


def hera_baseline_counts(
    k_perp: np.ndarray,
    z_obs: float,
    antenna_spacing: float = HERA_ANTENNA_SPACING_M,
    hex_n_side: int = HERA_HEX_N_SIDE,
    dish_diameter: float = 14.0,
    f_21_hz: float = 1420.405e6,
    speed_of_light_mps: float = 3e8,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
) -> np.ndarray:
    """
    Number of HERA baselines sampling each ``k_perp`` bin.

    Builds the close-packed hexagonal core, forms every antenna pair, converts
    each baseline length to ``u = |b|/λ`` and then to
    ``k_perp = 2π u / D_c(z)``, and counts how many fall in each bin.  This is
    what makes the physical noise model ``k_perp``-dependent: short baselines
    are hugely redundant on a hex array, so low ``k_perp`` is sampled far more
    deeply than high ``k_perp``.

    Parameters
    ----------
    k_perp : ndarray
        Bin centres [Mpc⁻¹].  Bin edges are taken as the geometric midpoints.
    z_obs : float
        Reference redshift.
    antenna_spacing : float, optional
        Core spacing [m].  See :data:`HERA_ANTENNA_SPACING_M`.
    hex_n_side : int, optional
        Antennas per hexagon side.  See :data:`HERA_HEX_N_SIDE`.
    dish_diameter : float, optional
        Dish diameter [m]; sets the ``uv``-cell size ``D/λ``.
    f_21_hz, speed_of_light_mps : float, optional
        21 cm rest frequency [Hz] and speed of light [m s⁻¹].
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology, for ``D_c(z)``.

    Returns
    -------
    ndarray
        Mean baselines **per independent uv-cell** in each bin, shape
        ``(len(k_perp),)``.  Zero where the array has no baselines — those
        modes are unmeasurable, and the noise model returns ``inf``.

    References
    ----------
    DeBoer et al. (2017), PASP 129, 045001 — 350 elements, 14 m dishes,
    14.6 m hexagonal-core spacing.
    Parsons et al. (2012a), ApJ 756, 165 — gridding baselines into
    ``uv``-cells of the antenna footprint ``D/λ``.

    Notes
    -----
    Counting *per uv-cell* rather than per ``k_perp`` bin is the point.  Only
    baselines landing in the same ``uv``-cell — within one antenna footprint
    ``D/λ`` — sample the same Fourier mode and integrate down coherently.
    Baselines elsewhere in the bin measure *different* modes; they add to the
    mode count, not to the depth of any one mode, and are accounted for by the
    ``mode_counts`` weighting in :func:`cross_power_snr` instead.  Summing all
    baselines in a bin would count the same redundancy twice, and on the
    fiducial grid would under-predict the noise several-fold.
    """
    # ── Hexagonal close-packed layout ─────────────────────────────────────
    positions = []
    for row in range(-hex_n_side + 1, hex_n_side):
        n_in_row = 2 * hex_n_side - 1 - abs(row)
        x_offset = -(n_in_row - 1) / 2.0
        for column in range(n_in_row):
            positions.append((
                (x_offset + column) * antenna_spacing,
                row * antenna_spacing * np.sqrt(3.0) / 2.0,
            ))
    antennas = np.asarray(positions)

    # ── Every antenna pair, as a uv vector ────────────────────────────────
    wavelength = speed_of_light_mps * (1.0 + z_obs) / f_21_hz
    delta = antennas[:, None, :] - antennas[None, :, :]
    upper = np.triu_indices(antennas.shape[0], k=1)
    u_coord = delta[..., 0][upper] / wavelength
    v_coord = delta[..., 1][upper] / wavelength

    # ── Grid into uv-cells one antenna footprint (D/λ) across ─────────────
    cell = dish_diameter / wavelength
    cell_index = np.stack(
        (np.round(u_coord / cell), np.round(v_coord / cell)), axis=1
    ).astype(np.int64)
    _, inverse, per_cell = np.unique(
        cell_index, axis=0, return_inverse=True, return_counts=True
    )
    # Representative |u| of each cell: the first baseline assigned to it.
    order = np.argsort(inverse, kind="stable")
    starts = np.concatenate(([0], np.cumsum(per_cell)[:-1]))
    cell_u = np.abs(np.hypot(u_coord, v_coord))[order][starts]

    # ── |u| -> k_perp, then average the per-cell counts in each bin ───────
    distance = comoving_distance(
        z_obs, hubble_constant, omega_m, speed_of_light_kms
    )
    k_perp_of_cell = 2.0 * np.pi * cell_u / distance

    centres = np.asarray(k_perp, dtype=float)
    inner = np.sqrt(centres[:-1] * centres[1:])
    edges = np.concatenate((
        [centres[0] ** 2 / inner[0]] if inner.size else [centres[0] * 0.5],
        inner,
        [centres[-1] ** 2 / inner[-1]] if inner.size else [centres[-1] * 2.0],
    ))
    which = np.digitize(k_perp_of_cell, edges) - 1
    valid = (which >= 0) & (which < centres.size)

    totals = np.bincount(which[valid], weights=per_cell[valid],
                         minlength=centres.size)
    occupied = np.bincount(which[valid], minlength=centres.size)
    return np.divide(
        totals, occupied, out=np.zeros_like(totals), where=occupied > 0
    )


def hera_thermal_noise_power_physical(
    k_perp: np.ndarray,
    z_obs: float,
    integration_time: float,
    dish_diameter: float = 14.0,
    f_21_hz: float = 1420.405e6,
    speed_of_light_mps: float = 3e8,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
    n_polarisations: int = N_POLARISATIONS,
    noise_equivalent_bandwidth: float = NOISE_EQUIVALENT_BANDWIDTH,
    **kwargs: float,
) -> np.ndarray:
    """
    HERA thermal-noise power spectrum, resolved in ``k_perp``.

    Implements

    ``P_N(k_perp) = X²Y · Ω_eff · NEB · T_sys² / [N_pol · t_int · N_bl(k_perp)]``

    which is Parsons (2017) Eq. 12 with the per-mode integration time written
    as ``t_int × N_bl(u)`` for a redundant array, and is algebraically
    identical to La Plante et al. (2023) Eq. 11,
    ``T_sys² Ω_P² X² Y / [Ω_PP t_int N_pol N_bl(u)]``, since
    ``Ω_eff ≡ Ω_P²/Ω_PP``.

    Unlike :func:`hera_thermal_noise_power` this is **not** flat in
    ``k_perp``: the hexagonal core is highly redundant on short baselines, so
    the noise rises steeply where few baselines sample the mode.

    Parameters
    ----------
    k_perp : ndarray
        Bin centres [Mpc⁻¹].
    z_obs : float
        Reference redshift.
    integration_time : float
        Total integration time [s].
    dish_diameter : float, optional
        Dish diameter [m].
    f_21_hz, speed_of_light_mps : float, optional
        21 cm rest frequency [Hz], speed of light [m s⁻¹].
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology.
    n_polarisations : int, optional
        Orthogonal polarisations combined.
    noise_equivalent_bandwidth : float, optional
        NEB of the line-of-sight taper; 1 for no taper.
    **kwargs
        Forwarded to :func:`hera_beam_solid_angles` and
        :func:`hera_baseline_counts` (``aperture_efficiency``,
        ``omega_p_over_pp``, ``antenna_spacing``, ``hex_n_side``).

    Returns
    -------
    ndarray
        ``P_N`` per ``k_perp`` bin [mK² Mpc³], shape ``(len(k_perp), 1)`` so it
        broadcasts against a ``(n_perp, n_par)`` spectrum.  ``inf`` where no
        baseline samples the bin.

    References
    ----------
    Parsons (2017), "Power Spectrum Normalizations for HERA" — Eq. 12.
    La Plante et al. (2023), arXiv:2205.09770 — Eq. 11.
    Parsons et al. (2014), ApJ 788, 106 — Appendix B, ``Ω_eff``.

    Notes
    -----
    Still an idealisation: it assumes every baseline integrates for the full
    ``integration_time`` on the same field, ignores the ``uv``-plane rotation
    that fills in coverage, and takes a single reference redshift.  For a
    publication forecast use `21cmSense
    <https://github.com/rasg-affiliates/21cmSense>`_.
    """
    cosmology = dict(
        hubble_constant=hubble_constant,
        omega_m=omega_m,
        speed_of_light_kms=speed_of_light_kms,
    )
    beam_kwargs = {
        key: kwargs[key] for key in ("aperture_efficiency", "omega_p_over_pp")
        if key in kwargs
    }
    array_kwargs = {
        key: kwargs[key] for key in ("antenna_spacing", "hex_n_side")
        if key in kwargs
    }

    _, omega_eff = hera_beam_solid_angles(
        z_obs, dish_diameter=dish_diameter, f_21_hz=f_21_hz,
        speed_of_light_mps=speed_of_light_mps, **beam_kwargs,
    )
    x2y = cosmological_scalar_x2y(z_obs, f_21_hz=f_21_hz, **cosmology)
    t_sys_mK, _ = system_temperature(z_obs, f_21_hz)
    n_baselines = hera_baseline_counts(
        k_perp, z_obs, dish_diameter=dish_diameter, f_21_hz=f_21_hz,
        speed_of_light_mps=speed_of_light_mps, **cosmology, **array_kwargs,
    )

    numerator = (
        x2y * omega_eff * noise_equivalent_bandwidth * t_sys_mK ** 2
    )
    denominator = n_polarisations * integration_time * n_baselines
    with np.errstate(divide="ignore"):
        power = np.where(n_baselines > 0, numerator / denominator, np.inf)
    return power[:, np.newaxis]


def system_temperature(
    z_obs: float,
    f_21_hz: float = 1420.405e6,
) -> Tuple[float, float]:
    """
    HERA system temperature and observed frequency at a reference redshift.

    ``T_sys = T_rcvr + T_sky(300 MHz / ν)^2.55``, with the constants from
    :data:`T_RECEIVER_K`, :data:`T_SKY_300MHZ_K` and
    :data:`SKY_SPECTRAL_INDEX`.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].

    Returns
    -------
    t_sys_mK : float
        System temperature [mK].
    observed_frequency : float
        Redshifted 21 cm frequency [Hz].

    References
    ----------
    Pober et al. (2013), AJ 145, 65; Pober et al. (2014), ApJ 782, 66 - the
    ``T_sky = 60 K (lambda/1 m)^2.55`` sky model with ``T_rcvr = 100 K``, as
    implemented in `21cmSense <https://github.com/jpober/21cmSense>`_
    (``calc_sense.py``: ``Tsky = 60e3 * (3e8/freq)**2.55``, in mK).
    DeBoer et al. (2017), PASP 129, 045001 - Table 2 gives a different
    normalisation, ``T_sys = 100 + 120(nu/150 MHz)^-2.55`` K; see the module
    constants for why this one is kept.
    """
    observed_frequency = f_21_hz / (1.0 + z_obs)
    t_sys_kelvin = T_RECEIVER_K + T_SKY_300MHZ_K * (
        300e6 / observed_frequency
    ) ** SKY_SPECTRAL_INDEX
    return float(t_sys_kelvin * 1e3), float(observed_frequency)


def hera_thermal_noise_power(
    z_obs: float,
    integration_time: float,
    bandwidth: float,
    f_21_hz: float = 1420.405e6,
) -> float:
    """
    Simplified HERA-like 21 cm thermal noise power.

    Uses the standard sky+receiver model
    ``T_sys = 100 K + 60 K (300 MHz / ν)^2.55`` and
    ``P_N = T_sys² × 10³ / (t_int Δν)``.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    integration_time : float
        Total integration time [s].
    bandwidth : float
        Per-band bandwidth [Hz].
    f_21_hz : float, optional
        21 cm rest frequency [Hz].  Note that the source notebook hardcodes
        ``1.42e9`` in its noise cell while using ``1420.405e6`` everywhere
        else; the pipeline uses the precise value throughout, which shifts
        ``P_N`` by 0.1 %.

    Returns
    -------
    float
        Thermal noise power in the same units as the 21 cm auto-spectrum
        [mK^2 Mpc^3].  Independent of ``k`` — the estimate carries no
        baseline-density or beam information.

    Notes
    -----
    This is a scaling estimate, not a full instrument model.  For publication
    forecasts replace it with `21cmSense
    <https://github.com/rasg-affiliates/21cmSense>`_, or with La Plante et al.
    (2023) Eq. 11, which resolves ``P_N(k_perp)`` through the baseline density.
    """
    t_sys_mK, _ = system_temperature(z_obs, f_21_hz)
    return float(
        t_sys_mK ** 2 * NOISE_NORMALISATION_MPC3 / (integration_time * bandwidth)
    )


def cross_power_snr(
    P_cross_observed: np.ndarray,
    P_21cm_auto: np.ndarray,
    P_galaxy_observed: np.ndarray,
    P_noise_21cm: float | np.ndarray,
    P_noise_galaxy: float,
    outside_wedge: Optional[np.ndarray] = None,
    mode_counts: Optional[np.ndarray] = None,
    n_patch: int = 1,
) -> SNRResult:
    """
    Per-mode cross-correlation SNR and its cumulative significance.

    Follows La Plante et al. (2023), Eqs. 15–17:
    ``σ_× = sqrt{0.5 [P_×² + (|P_21| + P_N,21)(|P_gal| + P_N,gal)]}``.

    Parameters
    ----------
    P_cross_observed : ndarray
        Photo-z damped cross-power, shape ``(n_perp, n_par)``.
    P_21cm_auto : ndarray
        21 cm auto-power, same shape.
    P_galaxy_observed : ndarray
        Photo-z damped galaxy auto-power, same shape.
    P_noise_21cm : float or ndarray
        21 cm thermal noise power.  A scalar (the flat estimate from
        :func:`hera_thermal_noise_power`) or an array broadcastable against
        the spectra — e.g. the ``(n_perp, 1)`` output of
        :func:`hera_thermal_noise_power_physical`.
    P_noise_galaxy : float
        Galaxy shot noise power (``1 / n̄``).
    outside_wedge : ndarray of bool, optional
        Mask of usable modes.  When omitted, all modes are used.
    mode_counts : ndarray, optional
        Fourier modes averaged into each bin, from
        :class:`~src.dataio.PowerSpectra`.  When given, the per-bin SNR is
        weighted by ``sqrt(n_patch × mode_counts / 2)`` — La Plante Eq. 19.
        When omitted (the default) the unweighted per-bin ratio is returned,
        which is what every result produced before this option existed used.
    n_patch : int, optional
        Independent survey patches, ``N_patch`` in Eq. 19.

    Returns
    -------
    SNRResult
        Per-mode maps, the two variance terms, and the total significance.

    Notes
    -----
    La Plante's Eqs. 15–17 carry factors of the brightness-temperature
    normalisation ``T_0(z)``: ``σ_21 = (P_21 + P_N,21)/T_0²`` and the SNR is
    ``|P_× / T_0| / σ_×``.  Those factors **cancel exactly** in the ratio, so
    the form used here — which omits ``T_0`` throughout, as the source
    notebook does — gives an identical SNR.

    **Mode weighting.**  Eq. 15 gives the variance of a *single* mode; Eq. 19
    combines bins as ``ŝ = sqrt(N_patch dN) P_×/σ_×`` with ``dN`` the
    independent-mode count of Eq. 18.  ``dN`` is ``mode_counts / 2``: the FFT
    of a real field is Hermitian, so half its cells are redundant.  Checked
    numerically against Eq. 18 on the fiducial grid, the two agree to a
    median 4.7 %.  Leaving ``mode_counts`` unset omits the weighting and
    under-reports the total SNR by roughly an order of magnitude — see
    ``docs/uncertainty_budget.md`` §6.5.
    """
    sigma_21cm = np.abs(P_21cm_auto) + P_noise_21cm
    sigma_galaxy = np.abs(P_galaxy_observed) + P_noise_galaxy

    cosmic_variance_term = 0.5 * P_cross_observed ** 2
    noise_coupling_term = 0.5 * sigma_21cm * sigma_galaxy
    sigma_cross = np.sqrt(cosmic_variance_term + noise_coupling_term)

    snr_per_mode = np.abs(P_cross_observed) / sigma_cross

    # La Plante Eq. 19: s_hat = sqrt(N_patch dN) P_x / sigma_x, dN = Eq. 18.
    mode_weight = None
    if mode_counts is not None:
        mode_weight = np.sqrt(
            n_patch * np.asarray(mode_counts, dtype=float) / 2.0
        )
        snr_per_mode = snr_per_mode * mode_weight

    if outside_wedge is None:
        outside_wedge = np.ones_like(snr_per_mode, dtype=bool)

    snr_outside_wedge = np.where(outside_wedge, snr_per_mode, np.nan)

    noise_21cm = np.asarray(P_noise_21cm, dtype=float)

    return SNRResult(
        snr_per_mode=snr_per_mode,
        snr_outside_wedge=snr_outside_wedge,
        total_snr=total_snr(snr_per_mode, outside_wedge),
        sigma_cross=sigma_cross,
        P_noise_21cm=(
            float(noise_21cm) if noise_21cm.ndim == 0 else noise_21cm
        ),
        P_noise_galaxy=float(P_noise_galaxy),
        mode_weight=mode_weight,
        sigma_21cm=sigma_21cm,
        sigma_galaxy=sigma_galaxy,
        cosmic_variance_term=cosmic_variance_term,
        noise_coupling_term=noise_coupling_term,
    )


def total_snr(snr_per_mode: np.ndarray, mask: Optional[np.ndarray] = None) -> float:
    """
    Quadrature sum of the per-mode SNR over the unmasked bins.

    Parameters
    ----------
    snr_per_mode : ndarray
        Per-mode signal-to-noise.
    mask : ndarray of bool, optional
        Bins to include.  When omitted, all bins are used.

    Returns
    -------
    float
        Total detection significance [σ].
    """
    values = snr_per_mode if mask is None else snr_per_mode[mask]
    return float(np.sqrt(np.nansum(np.asarray(values) ** 2)))


# ===========================================================================
# 5b  The complete uncertainty budget
# ===========================================================================

@dataclass
class UncertaintyBudget:
    """
    Every term of the 21 cm × galaxy cross-power uncertainty budget.

    This is the pipeline's transcription of the photo-z / wedge / noise / SNR
    chain in ``21cmfast_HERAxEuclid_lightcone.ipynb`` (its "Photo-z Radial
    Smearing & Foreground Wedge Mask" and "Cross-correlation SNR map"
    sections), computed by :func:`compute_uncertainty_budget`.

    Attributes
    ----------
    k_perp, k_parallel : ndarray
        Bin centres the budget was evaluated on [Mpc^-1].
    z_obs : float
        Reference redshift.
    photoz_uncertainty : float
        Absolute σ_z used (**not** σ_z/(1+z)).
    radial_smearing : float
        σ_r = c σ_z / H(z_obs) [Mpc].
    photoz_kernel : ndarray
        W(k_par), shape ``(1, n_par)``.
    P_cross_observed : ndarray
        ``P_cross × W`` — one factor, only the galaxy side is smeared.
    P_galaxy_observed : ndarray
        ``P_gal × W²`` — both fields in the pair are smeared.
    horizon_slope, fov_slope : float
        Wedge slopes; the **horizon** slope defines the mask, the FoV slope is
        drawn on figures only.
    wedge_buffer : float
        Margin added above the wedge line [Mpc^-1].
    outside_wedge : ndarray of bool
        True where a bin survives wedge excision.
    observed_frequency_hz : float
        Redshifted 21 cm frequency [Hz].
    system_temperature_mK : float
        T_sys at that frequency [mK].
    integration_time, bandwidth : float
        Instrument configuration [s], [Hz].
    mean_galaxy_density : float
        n̄ used for the shot noise [Mpc^-3].
    noise_model : str
        Which thermal-noise model produced ``snr.P_noise_21cm`` —
        ``'scaling'`` or ``'physical'``.
    mode_weighted : bool
        Whether the La Plante Eq. 19 mode weighting was applied.
    snr : SNRResult
        Per-mode maps, variance terms, and the total significance.

    Notes
    -----
    With ``mode_weighted=False`` (the default) the per-bin ``mode_counts`` are
    **not** folded into the total: the quoted SNR is a per-bin quantity summed
    in quadrature, exactly as the notebook computes it, and is roughly an
    order of magnitude smaller than the Eq. 19 weighted total.  See
    ``docs/uncertainty_budget.md`` §6.5.
    """

    k_perp: np.ndarray
    k_parallel: np.ndarray
    z_obs: float

    photoz_uncertainty: float
    radial_smearing: float
    photoz_kernel: np.ndarray
    P_cross_observed: np.ndarray
    P_galaxy_observed: np.ndarray

    horizon_slope: float
    fov_slope: float
    wedge_buffer: float
    outside_wedge: np.ndarray

    observed_frequency_hz: float
    system_temperature_mK: float
    integration_time: float
    bandwidth: float
    mean_galaxy_density: float
    noise_model: str = "scaling"
    mode_weighted: bool = False

    snr: SNRResult = None            # type: ignore[assignment]

    # ── Derived summaries ─────────────────────────────────────────────────
    @property
    def total_snr(self) -> float:
        """Total detection significance outside the wedge [σ]."""
        return self.snr.total_snr

    @property
    def fraction_outside_wedge(self) -> float:
        """Fraction of ``(k_perp, k_parallel)`` bins surviving the wedge."""
        return float(self.outside_wedge.mean())

    @property
    def detected(self) -> bool:
        """True when the total significance exceeds the 5σ convention."""
        return bool(self.snr.total_snr > 5.0)

    @property
    def cosmic_variance_fraction(self) -> float:
        """
        Share of ``σ_cross²`` carried by sample variance, outside the wedge.

        Near 1 the budget is cosmic-variance limited (more integration time
        will not help); near 0 it is noise limited.  ``NaN`` when no usable
        bin has a finite variance.
        """
        usable = self.outside_wedge & np.isfinite(self.snr.sigma_cross)
        if not usable.any():
            return float("nan")
        cosmic = float(np.nansum(self.snr.cosmic_variance_term[usable]))
        noise = float(np.nansum(self.snr.noise_coupling_term[usable]))
        if cosmic + noise <= 0:
            return float("nan")
        return cosmic / (cosmic + noise)

    def as_dict(self) -> dict:
        """
        Scalar summary of the budget, ready for the pipeline summary JSON.

        Returns
        -------
        dict
            Every scalar term, JSON-serialisable.
        """
        return {
            "z_obs": float(self.z_obs),
            "photoz_uncertainty_sigma_z": float(self.photoz_uncertainty),
            "radial_smearing_Mpc": float(self.radial_smearing),
            "photoz_kernel_first_bin": float(self.photoz_kernel.ravel()[0]),
            "photoz_kernel_max": float(self.photoz_kernel.max()),
            "horizon_wedge_slope": float(self.horizon_slope),
            "fov_wedge_slope": float(self.fov_slope),
            "wedge_buffer_Mpc-1": float(self.wedge_buffer),
            "modes_outside_wedge": int(self.outside_wedge.sum()),
            "modes_total": int(self.outside_wedge.size),
            "modes_outside_wedge_fraction": self.fraction_outside_wedge,
            "observed_frequency_MHz": self.observed_frequency_hz / 1e6,
            "system_temperature_K": self.system_temperature_mK / 1e3,
            "integration_time_s": float(self.integration_time),
            "bandwidth_Hz": float(self.bandwidth),
            "mean_galaxy_density": float(self.mean_galaxy_density),
            "noise_model": self.noise_model,
            "mode_weighted": bool(self.mode_weighted),
            # k_perp-resolved under noise_model="physical"; reduced to its
            # finite range so the summary stays scalar and JSON-serialisable.
            "P_noise_21cm": (
                float(self.snr.P_noise_21cm)
                if np.ndim(self.snr.P_noise_21cm) == 0
                else float(np.min(self.snr.P_noise_21cm))
            ),
            "P_noise_21cm_max": (
                float(self.snr.P_noise_21cm)
                if np.ndim(self.snr.P_noise_21cm) == 0
                else float(np.max(
                    np.asarray(self.snr.P_noise_21cm)[
                        np.isfinite(self.snr.P_noise_21cm)
                    ]
                ))
            ),
            "P_noise_galaxy": self.snr.P_noise_galaxy,
            "cosmic_variance_fraction": self.cosmic_variance_fraction,
            "total_snr_sigma": self.snr.total_snr,
            "detection_above_5sigma": self.detected,
        }


def compute_uncertainty_budget(
    spectra: PowerSpectra,
    z_obs: float,
    photoz_uncertainty: float = 0.45,
    wedge_buffer: float = 0.0677,
    integration_time: float = 1000 * 3600,
    bandwidth: float = 8e6,
    mean_galaxy_density: float = 7.48e-5,
    dish_diameter: float = 14.0,
    f_21_hz: float = 1420.405e6,
    speed_of_light_mps: float = 3e8,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
    noise_model: str = "scaling",
    mode_weighted: bool = False,
    n_patch: int = 1,
) -> UncertaintyBudget:
    """
    Run the full uncertainty-budget chain on a set of power spectra.

    Applies, in order: photo-z damping of the galaxy and cross spectra,
    foreground-wedge excision, the HERA thermal-noise and galaxy shot-noise
    estimates, and the La Plante et al. (2023) per-mode variance and SNR.

    This is the single entry point for the budget — ``run_pipeline.py`` calls
    nothing else, so the notebook and the HPC run cannot drift apart.

    Parameters
    ----------
    spectra : PowerSpectra
        Undamped spectra from :func:`compute_all_power_spectra`.
    z_obs : float
        Reference redshift for the wedge, smearing, and noise.
    photoz_uncertainty : float, optional
        Absolute σ_z (**not** σ_z/(1+z)).  Default 0.45, the Euclid
        requirement σ_z/(1+z) < 0.05 evaluated at z = 7.
    wedge_buffer : float, optional
        Margin above the wedge line [Mpc^-1].  Default 0.0677 = 0.1 h Mpc^-1
        (Pober et al. 2014 "moderate"; the 21cmSense default).
    integration_time : float, optional
        Total integration [s].  Default 1000 h.
    bandwidth : float, optional
        Per-band bandwidth [Hz].
    mean_galaxy_density : float, optional
        n̄ for the shot noise ``P_N,gal = 1/n̄`` [Mpc^-3]. Euclid Collab.:
        Allen et al. (2026), A&A 711, A25, Table 2; see NUMBERS_AND_SOURCES.md
        §2 for the derivation and the cosmology caveat.
    dish_diameter : float, optional
        Dish diameter [m], for the FoV wedge line.
    f_21_hz : float, optional
        21 cm rest frequency [Hz].
    speed_of_light_mps : float, optional
        Speed of light [m s^-1], for the observed wavelength.
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology.
    noise_model : {'scaling', 'physical'}, optional
        Which 21 cm thermal-noise model to use.

        - ``'scaling'`` (default) — :func:`hera_thermal_noise_power`, flat in
          ``k`` and carrying an arbitrary Mpc³ normalisation.  Reproduces every
          number this pipeline produced before the ``'physical'`` option
          existed.
        - ``'physical'`` — :func:`hera_thermal_noise_power_physical`,
          Parsons (2017) Eq. 12 / La Plante Eq. 11, resolved in ``k_perp``
          through the HERA baseline distribution.  Roughly 10³ larger than
          the scaling estimate and, unlike it, infinite where no baseline
          samples the mode.
    mode_weighted : bool, optional
        Apply the La Plante Eq. 19 weighting ``sqrt(N_patch dN)`` when summing
        bins, using the estimator's own ``mode_counts``.  Default ``False``,
        which preserves the pre-existing unweighted total.
    n_patch : int, optional
        Independent survey patches, ``N_patch`` in Eq. 19.

    Returns
    -------
    UncertaintyBudget
        Every term of the budget, plus the SNR maps.

    Raises
    ------
    ValueError
        If ``noise_model`` is not one of the two accepted values.

    References
    ----------
    La Plante et al. (2023), arXiv:2205.09770 — Eqs. 10, 11, 15–17, 18–20.
    Parsons (2017), "Power Spectrum Normalizations for HERA" — Eq. 12.
    Pober et al. (2014), ApJ 782, 66 — the wedge buffer.

    Notes
    -----
    Both ``noise_model='physical'`` and ``mode_weighted=True`` change every
    number the budget reports, so both default to the historical behaviour.
    ``docs/uncertainty_budget.md`` §6.5 quantifies each against the papers.
    """
    if noise_model not in ("scaling", "physical"):
        raise ValueError(
            f"noise_model must be 'scaling' or 'physical', got {noise_model!r}"
        )
    cosmology = dict(
        hubble_constant=hubble_constant,
        omega_m=omega_m,
        speed_of_light_kms=speed_of_light_kms,
    )

    # ── Photo-z damping: one factor of W on the cross, two on the auto ────
    radial_smearing = radial_smearing_length(
        photoz_uncertainty=photoz_uncertainty, z_obs=z_obs, **cosmology
    )
    kernel = photoz_damping_kernel(spectra.k_parallel, radial_smearing)
    p_cross_observed = spectra.P_cross * kernel
    p_galaxy_observed = spectra.P_galaxy_auto * kernel ** 2

    # ── Foreground wedge: horizon slope masks, FoV slope is for figures ───
    horizon_slope = horizon_wedge_slope(z_obs, **cosmology)
    fov_slope = fov_wedge_slope(
        z_obs,
        dish_diameter=dish_diameter,
        f_21_hz=f_21_hz,
        speed_of_light_mps=speed_of_light_mps,
        **cosmology,
    )
    outside_wedge = foreground_wedge_mask(
        spectra.k_perp, spectra.k_parallel,
        slope=horizon_slope, buffer=wedge_buffer,
    )

    # ── Noise ─────────────────────────────────────────────────────────────
    t_sys_mK, observed_frequency = system_temperature(z_obs, f_21_hz)
    if noise_model == "physical":
        p_noise_21cm = hera_thermal_noise_power_physical(
            k_perp=spectra.k_perp,
            z_obs=z_obs,
            integration_time=integration_time,
            dish_diameter=dish_diameter,
            f_21_hz=f_21_hz,
            speed_of_light_mps=speed_of_light_mps,
            **cosmology,
        )
    else:
        p_noise_21cm = hera_thermal_noise_power(
            z_obs=z_obs,
            integration_time=integration_time,
            bandwidth=bandwidth,
            f_21_hz=f_21_hz,
        )
    p_noise_galaxy = 1.0 / mean_galaxy_density

    # ── Variance and SNR ──────────────────────────────────────────────────
    snr = cross_power_snr(
        P_cross_observed=p_cross_observed,
        P_21cm_auto=spectra.P_21cm_auto,
        P_galaxy_observed=p_galaxy_observed,
        P_noise_21cm=p_noise_21cm,
        P_noise_galaxy=p_noise_galaxy,
        outside_wedge=outside_wedge,
        mode_counts=spectra.mode_counts if mode_weighted else None,
        n_patch=n_patch,
    )

    return UncertaintyBudget(
        k_perp=spectra.k_perp,
        k_parallel=spectra.k_parallel,
        z_obs=float(z_obs),
        photoz_uncertainty=float(photoz_uncertainty),
        radial_smearing=radial_smearing,
        photoz_kernel=kernel,
        P_cross_observed=p_cross_observed,
        P_galaxy_observed=p_galaxy_observed,
        horizon_slope=horizon_slope,
        fov_slope=fov_slope,
        wedge_buffer=float(wedge_buffer),
        outside_wedge=outside_wedge,
        observed_frequency_hz=observed_frequency,
        system_temperature_mK=t_sys_mK,
        integration_time=float(integration_time),
        bandwidth=float(bandwidth),
        mean_galaxy_density=float(mean_galaxy_density),
        noise_model=noise_model,
        mode_weighted=bool(mode_weighted),
        snr=snr,
    )


# ===========================================================================
# 6  Euclid selection and effective galaxy bias
# ===========================================================================

@dataclass
class EuclidSelection:
    """
    Result of applying the Euclid ``M_UV`` window to the halo catalogue.

    Attributes
    ----------
    mask : ndarray of bool
        Selection mask over the *valid* (SFR > 0) halos.
    M_UV : ndarray
        Absolute UV magnitudes of the selected halos.
    sfr : ndarray
        Star-formation rates of the selected halos [M_sun yr^-1].
    halo_masses : ndarray
        Halo masses of the selected halos [M_sun].
    n_valid : int
        Number of halos with SFR > 0 (and mass > 0) before the cut.
    n_selected : int
        Number of halos passing the magnitude cut.
    SFR_min, SFR_max : float
        SFR bounds equivalent to the magnitude window [M_sun yr^-1].
    M_UV_bright, M_UV_faint : float
        Magnitude window used.
    """

    mask: np.ndarray
    M_UV: np.ndarray
    sfr: np.ndarray
    halo_masses: np.ndarray
    n_valid: int
    n_selected: int
    SFR_min: float
    SFR_max: float
    M_UV_bright: float
    M_UV_faint: float


@dataclass
class BiasEstimate:
    """
    Effective linear galaxy bias of the Euclid-selected halo sample.

    Attributes
    ----------
    mean_bias : float
        Catalogue-averaged Sheth-Tormen bias ⟨b_g⟩.
    bias_min, bias_max : float
        Range of per-halo bias values.
    log10_mass_grid : ndarray
        log10(M / M_sun h^-1) grid used for the bias curve.
    bias_grid : ndarray
        Interpolated ``b_h(M, z)`` on ``log10_mass_grid``.
    log10_selected_mass_h : ndarray
        log10 halo masses [M_sun/h] of the selected sample.
    n_selected : int
        Number of selected halos.
    """

    mean_bias: float
    bias_min: float
    bias_max: float
    log10_mass_grid: np.ndarray
    bias_grid: np.ndarray
    log10_selected_mass_h: np.ndarray
    n_selected: int


def euclid_sfr_window(
    M_UV_faint: float = -18.0,
    M_UV_bright: float = -22.0,
) -> Tuple[float, float]:
    """
    SFR bounds equivalent to a Euclid absolute-magnitude window.

    Parameters
    ----------
    M_UV_faint : float, optional
        Faint-end (least negative) magnitude limit.
    M_UV_bright : float, optional
        Bright-end (most negative) magnitude limit.

    Returns
    -------
    sfr_min : float
        SFR corresponding to ``M_UV_faint`` [M_sun yr^-1].
    sfr_max : float
        SFR corresponding to ``M_UV_bright`` [M_sun yr^-1].

    Notes
    -----
    Uses the adopted UV–SFR calibration κ_UV = 2.7e-29 via
    ``src.conversions`` — Fisher et al. (2026), arXiv:2511.10741, Eq. 12.

    .. warning::
       That κ_UV recovers ``SFR_100Myr`` from *rising* SFHs, while the
       ``halo_sfr`` this window is applied to is 21cmFAST's
       ``M_star / t_sf`` with ``t_sf = 570.3 Myr`` at z = 7 (Park et al.
       2019, Eq. 3).  The two are not the same quantity, so this window's
       meaning — not merely its numbers — is an open question.  See the
       ``_KAPPA_UV_MADAU14`` definition in ``src/conversions.py`` and
       ``NUMBERS_AND_SOURCES.md`` §2.
    """
    return (
        float(Luv_to_sfr(Muv_to_Luv(M_UV_faint))),
        float(Luv_to_sfr(Muv_to_Luv(M_UV_bright))),
    )


def select_euclid_halos(
    sfr: np.ndarray,
    halo_masses: np.ndarray,
    M_UV_faint: float = -18.0,
    M_UV_bright: float = -22.0,
) -> EuclidSelection:
    """
    Apply the Euclid ``M_UV`` window to a 21cmFAST halo catalogue.

    Each halo's SFR is converted to a UV luminosity and absolute magnitude
    (Madau & Dickinson 2014), then filtered by
    ``M_UV_bright <= M_UV <= M_UV_faint``.

    Parameters
    ----------
    sfr : ndarray
        Per-halo star-formation rate [M_sun yr^-1].
    halo_masses : ndarray
        Per-halo mass [M_sun], same length as ``sfr``.
    M_UV_faint, M_UV_bright : float, optional
        Magnitude window.

    Returns
    -------
    EuclidSelection
        The selected sample and its diagnostics.

    Raises
    ------
    ValueError
        If ``sfr`` and ``halo_masses`` have different lengths.
    """
    sfr = np.asarray(sfr)
    halo_masses = np.asarray(halo_masses)
    if sfr.shape[0] != halo_masses.shape[0]:
        raise ValueError(
            f"sfr and halo_masses length mismatch: "
            f"{sfr.shape[0]} vs {halo_masses.shape[0]}"
        )

    sfr_clean = np.where(np.isfinite(sfr), sfr, 0.0)
    valid = (sfr_clean > 0) & (halo_masses > 0)

    sfr_valid = sfr_clean[valid]
    mass_valid = halo_masses[valid]

    M_UV = Luv_to_Muv(sfr_to_Luv(sfr_valid))
    mask = (M_UV >= M_UV_bright) & (M_UV <= M_UV_faint)

    sfr_min, sfr_max = euclid_sfr_window(M_UV_faint, M_UV_bright)

    return EuclidSelection(
        mask=mask,
        M_UV=M_UV[mask],
        sfr=sfr_valid[mask],
        halo_masses=mass_valid[mask],
        n_valid=int(valid.sum()),
        n_selected=int(mask.sum()),
        SFR_min=sfr_min,
        SFR_max=sfr_max,
        M_UV_bright=float(M_UV_bright),
        M_UV_faint=float(M_UV_faint),
    )


def effective_galaxy_bias(
    selection: EuclidSelection,
    z_obs: float,
    hubble_constant: float = 67.36,
    dlog10m: float = 0.02,
) -> BiasEstimate:
    """
    Effective linear bias of a Euclid-selected halo sample.

    Sheth-Tormen bias is evaluated on an ``hmf`` mass-function grid at
    ``z_obs`` and interpolated onto the selected halo masses.  Halo masses
    are converted from M_sun to M_sun/h to match ``hmf``'s convention.

    Parameters
    ----------
    selection : EuclidSelection
        Output of :func:`select_euclid_halos`.
    z_obs : float
        Reference redshift.
    hubble_constant : float, optional
        H_0 [km s^-1 Mpc^-1]; ``h = H_0 / 100``.
    dlog10m : float, optional
        Mass-function grid resolution in dex.

    Returns
    -------
    BiasEstimate
        Mean bias and the bias-vs-mass curve used to derive it.

    Raises
    ------
    ImportError
        If the ``hmf`` package is not installed.
    ValueError
        If no halos passed the selection.
    """
    from hmf import MassFunction
    from scipy.interpolate import interp1d

    if selection.n_selected == 0:
        raise ValueError("no halos passed the Euclid magnitude selection")

    h_cosmo = hubble_constant / 100.0
    selected_mass_h = selection.halo_masses.astype(np.float64) * h_cosmo

    log10_m_min = float(np.floor(np.log10(selected_mass_h.min())) - 0.5)
    log10_m_max = float(np.ceil(np.log10(selected_mass_h.max())) + 0.5)

    mass_function = MassFunction(
        z=z_obs,
        Mmin=log10_m_min,
        Mmax=log10_m_max,
        dlog10m=dlog10m,
    )

    # hmf stores nu as the *squared* peak height (δ_c/σ)², which is what
    # sheth_tormen_bias expects.
    hmf_bias = sheth_tormen_bias(mass_function.nu)

    bias_interp = interp1d(
        np.log10(mass_function.m),
        hmf_bias,
        bounds_error=False,
        fill_value="extrapolate",
    )

    log10_selected_mass_h = np.log10(selected_mass_h)
    selected_biases = bias_interp(log10_selected_mass_h)

    mass_grid = np.linspace(log10_m_min, log10_m_max, 300)

    return BiasEstimate(
        mean_bias=float(np.mean(selected_biases)),
        bias_min=float(np.min(selected_biases)),
        bias_max=float(np.max(selected_biases)),
        log10_mass_grid=mass_grid,
        bias_grid=np.asarray(bias_interp(mass_grid)),
        log10_selected_mass_h=log10_selected_mass_h,
        n_selected=selection.n_selected,
    )


# ===========================================================================
# 9  Galaxy overdensity from the halo catalogue
# ===========================================================================
#
# The Euclid-selected halo catalogue can be deposited onto the simulation
# grid with two different weights, and the choice matters physically:
#
#   "number"     — every selected halo contributes 1, so delta_gal traces
#                  the *abundance* of detectable galaxies.  This is the
#                  quantity a galaxy clustering measurement counts.
#
#   "luminosity" — every selected halo contributes its own L_UV, so
#                  delta_gal traces the *UV emissivity*.  Bright halos are
#                  up-weighted, which is closer to what a flux-limited
#                  intensity-mapping style measurement responds to.
#
# L_UV comes from sfr_to_Luv() (Fisher et al. 2026, arXiv:2511.10741 Eq. 12,
# kappa_UV = 2.7e-29); this module only consumes it.  The value cancels out of
# an overdensity, so the weighting mode is insensitive to it -- see
# tests/test_galaxy_weighting.py.

GALAXY_WEIGHTING_MODES: Tuple[str, ...] = ("number", "luminosity")


def deposit_halo_field(
    coords: np.ndarray,
    box_len: float,
    n_perp: int,
    n_los: Optional[int] = None,
    los_extent: Optional[float] = None,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Deposit halos onto a regular 3D grid by nearest-cell (CIC-free) binning.

    This is the 3D generalisation of the projected map built by
    ``figures.plot_uv_selection_maps``: a plain histogram of halo positions,
    optionally weighted per halo.

    Parameters
    ----------
    coords : ndarray
        Halo positions of shape ``(N_halos, 3)`` in Mpc.
    box_len : float
        Transverse box side length [Mpc]; sets the ``x`` and ``y`` extent.
    n_perp : int
        Number of transverse cells per side.
    n_los : int, optional
        Number of line-of-sight cells.  Defaults to ``n_perp`` (cubic grid).
    los_extent : float, optional
        Line-of-sight extent [Mpc].  Defaults to ``box_len``.
    weights : ndarray, optional
        Per-halo weight of shape ``(N_halos,)``.  ``None`` means unit
        weights, i.e. a pure number count.

    Returns
    -------
    ndarray
        Grid of shape ``(n_perp, n_perp, n_los)`` holding the summed weight
        per cell.

    Raises
    ------
    ValueError
        If ``coords`` is not ``(N, 3)``, or ``weights`` has a different
        length than ``coords``.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords must have shape (N_halos, 3), got {coords.shape}")

    n_los = int(n_perp) if n_los is None else int(n_los)
    los_extent = float(box_len) if los_extent is None else float(los_extent)

    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        if weights.shape[0] != coords.shape[0]:
            raise ValueError(
                f"weights and coords length mismatch: "
                f"{weights.shape[0]} vs {coords.shape[0]}"
            )

    if coords.shape[0] == 0:
        return np.zeros((int(n_perp), int(n_perp), n_los), dtype=float)

    edges = (
        np.linspace(0.0, float(box_len), int(n_perp) + 1),
        np.linspace(0.0, float(box_len), int(n_perp) + 1),
        np.linspace(0.0, los_extent, n_los + 1),
    )
    field, _ = np.histogramdd(coords, bins=edges, weights=weights)
    return field


def galaxy_overdensity_from_catalogue(
    coords: np.ndarray,
    sfr: np.ndarray,
    halo_masses: np.ndarray,
    box_len: float,
    n_perp: int,
    n_los: Optional[int] = None,
    los_extent: Optional[float] = None,
    weighting: str = "number",
    M_UV_faint: float = -18.0,
    M_UV_bright: float = -22.0,
    apply_selection: bool = True,
) -> Tuple[np.ndarray, EuclidSelection]:
    """
    Galaxy overdensity field from the halo catalogue, number- or L_UV-weighted.

    Both modes deposit the same halos onto the same grid and normalise the
    same way; only the per-halo weight differs::

        "number"      delta_gal   = N / <N> - 1
        "luminosity"  delta_gal,L = sum(L_UV) / <sum(L_UV)> - 1

    The two results are therefore interchangeable downstream — same shape,
    same zero mean, same units (dimensionless).

    Parameters
    ----------
    coords : ndarray
        Halo positions of shape ``(N_halos, 3)`` in Mpc.
    sfr : ndarray
        Per-halo star-formation rate [M_sun yr^-1].
    halo_masses : ndarray
        Per-halo mass [M_sun].
    box_len : float
        Transverse box side length [Mpc].
    n_perp : int
        Transverse cells per side.
    n_los : int, optional
        Line-of-sight cells.  Defaults to ``n_perp``.
    los_extent : float, optional
        Line-of-sight extent [Mpc].  Defaults to ``box_len``.
    weighting : {'number', 'luminosity'}, optional
        Per-halo weight.  ``'number'`` (default) reproduces the existing
        number-count field; ``'luminosity'`` weights by ``sfr_to_Luv(sfr)``.
    M_UV_faint, M_UV_bright : float, optional
        Euclid magnitude window passed to :func:`select_euclid_halos`.
    apply_selection : bool, optional
        Apply the Euclid magnitude window before depositing.  ``False``
        deposits every halo with SFR > 0.

    Returns
    -------
    delta_gal : ndarray
        Overdensity of shape ``(n_perp, n_perp, n_los)``.  All-zero if no
        halo survives the selection or the mean weight is zero.
    selection : EuclidSelection
        The selection actually applied, for diagnostics and logging.

    Raises
    ------
    ValueError
        If ``weighting`` is not one of ``GALAXY_WEIGHTING_MODES``, or the
        catalogue arrays have inconsistent lengths.

    Notes
    -----
    ``L_UV`` is obtained from :func:`src.conversions.sfr_to_Luv`, i.e.
    ``L_UV = SFR / kappa_UV`` with ``kappa_UV = 2.7e-29`` (Fisher et al.
    2026, arXiv:2511.10741, Eq. 12).  Because that conversion is a
    constant rescaling, the
    luminosity-weighted field is identical to an SFR-weighted one; it
    differs from the number-weighted field only through the per-halo
    spread in SFR.
    """
    if weighting not in GALAXY_WEIGHTING_MODES:
        raise ValueError(
            f"weighting must be one of {GALAXY_WEIGHTING_MODES}, got {weighting!r}"
        )

    coords = np.asarray(coords, dtype=float)
    sfr = np.asarray(sfr, dtype=float)
    halo_masses = np.asarray(halo_masses, dtype=float)

    if not (coords.shape[0] == sfr.shape[0] == halo_masses.shape[0]):
        raise ValueError(
            f"catalogue length mismatch: coords={coords.shape[0]}, "
            f"sfr={sfr.shape[0]}, halo_masses={halo_masses.shape[0]}"
        )

    selection = select_euclid_halos(
        sfr, halo_masses, M_UV_faint=M_UV_faint, M_UV_bright=M_UV_bright,
    )

    # select_euclid_halos masks *within* the valid (SFR > 0, M > 0) subset,
    # so lift its mask back to full-catalogue indices before slicing coords.
    sfr_clean = np.where(np.isfinite(sfr), sfr, 0.0)
    valid = (sfr_clean > 0) & (halo_masses > 0)
    keep = np.zeros(coords.shape[0], dtype=bool)
    if apply_selection:
        keep[np.flatnonzero(valid)[selection.mask]] = True
    else:
        keep = valid

    weights = sfr_to_Luv(sfr_clean[keep]) if weighting == "luminosity" else None

    field = deposit_halo_field(
        coords[keep], box_len=box_len, n_perp=n_perp,
        n_los=n_los, los_extent=los_extent, weights=weights,
    )

    mean_weight = float(field.mean())
    if mean_weight <= 0.0:
        return np.zeros_like(field), selection

    return field / mean_weight - 1.0, selection
