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
5. Euclid ``M_UV`` selection of the halo catalogue and the resulting
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
    "horizon_wedge_slope",
    "fov_wedge_slope",
    "foreground_wedge_mask",
    "photoz_damping_kernel",
    "radial_smearing_length",
    "hera_thermal_noise_power",
    "cross_power_snr",
    "total_snr",
    "euclid_sfr_window",
    "select_euclid_halos",
    "effective_galaxy_bias",
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

    dx = box_len_perp / n_perp   # transverse cell size [Mpc]
    dz = box_len_los / n_los     # LOS cell size [Mpc]
    volume = box_len_perp ** 2 * box_len_los

    # ── Fourier transforms ────────────────────────────────────────────────
    ft_factor = dx * dx * dz
    fourier_a = np.fft.fftn(field_a) * ft_factor
    fourier_b = np.fft.fftn(field_b) * ft_factor
    power_3d = (fourier_a * np.conj(fourier_b)).real / volume

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
) -> PowerSpectra:
    """
    Compute the 21 cm auto-, galaxy auto-, and cross-power spectra.

    The 21 cm field is converted to fluctuations by subtracting its mean; the
    galaxy field is already an overdensity.

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

    Returns
    -------
    PowerSpectra
        The three spectra on a shared k-grid, plus the mode counts.
    """
    t21_fluctuations = brightness_temp_field - brightness_temp_field.mean()

    k_perp, k_parallel, p_21, mode_counts = compute_cylindrical_cross_power(
        t21_fluctuations, t21_fluctuations,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel,
    )
    _, _, p_gal, _ = compute_cylindrical_cross_power(
        galaxy_overdensity, galaxy_overdensity,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel,
    )
    _, _, p_cross, _ = compute_cylindrical_cross_power(
        t21_fluctuations, galaxy_overdensity,
        box_len_perp, box_len_los, n_bins_perp, n_bins_parallel,
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
    buffer: float = 0.02,
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
        Safety margin added above the wedge line [Mpc^-1].

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
        Photometric redshift uncertainty σ_z.
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
    P_noise_21cm : float
        21 cm thermal noise power [mK^2 Mpc^3].
    P_noise_galaxy : float
        Galaxy shot noise power [Mpc^3].
    """

    snr_per_mode: np.ndarray
    snr_outside_wedge: np.ndarray
    total_snr: float
    sigma_cross: np.ndarray
    P_noise_21cm: float
    P_noise_galaxy: float


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
    ``P_N ∝ T_sys² / (t_int Δν)``.

    Parameters
    ----------
    z_obs : float
        Reference redshift.
    integration_time : float
        Total integration time [s].
    bandwidth : float
        Per-band bandwidth [Hz].
    f_21_hz : float, optional
        21 cm rest frequency [Hz].

    Returns
    -------
    float
        Thermal noise power in the same units as the 21 cm auto-spectrum
        [mK^2 Mpc^3].

    Notes
    -----
    This is a scaling estimate, not a full instrument model.  For publication
    forecasts replace it with `21cmSense
    <https://github.com/rasg-affiliates/21cmSense>`_.
    """
    observed_frequency = f_21_hz / (1.0 + z_obs)                    # [Hz]
    system_temperature = (
        100.0 + 60.0 * (300e6 / observed_frequency) ** 2.55
    ) * 1e3                                                          # [mK]
    return float(
        system_temperature ** 2 * 1e3 / (integration_time * bandwidth)
    )


def cross_power_snr(
    P_cross_observed: np.ndarray,
    P_21cm_auto: np.ndarray,
    P_galaxy_observed: np.ndarray,
    P_noise_21cm: float,
    P_noise_galaxy: float,
    outside_wedge: Optional[np.ndarray] = None,
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
    P_noise_21cm : float
        21 cm thermal noise power.
    P_noise_galaxy : float
        Galaxy shot noise power (``1 / n̄``).
    outside_wedge : ndarray of bool, optional
        Mask of usable modes.  When omitted, all modes are used.

    Returns
    -------
    SNRResult
        Per-mode maps and the total significance.
    """
    sigma_21cm = np.abs(P_21cm_auto) + P_noise_21cm
    sigma_galaxy = np.abs(P_galaxy_observed) + P_noise_galaxy
    sigma_cross = np.sqrt(
        0.5 * (P_cross_observed ** 2 + sigma_21cm * sigma_galaxy)
    )

    snr_per_mode = np.abs(P_cross_observed) / sigma_cross

    if outside_wedge is None:
        outside_wedge = np.ones_like(snr_per_mode, dtype=bool)

    snr_outside_wedge = np.where(outside_wedge, snr_per_mode, np.nan)

    return SNRResult(
        snr_per_mode=snr_per_mode,
        snr_outside_wedge=snr_outside_wedge,
        total_snr=total_snr(snr_per_mode, outside_wedge),
        sigma_cross=sigma_cross,
        P_noise_21cm=float(P_noise_21cm),
        P_noise_galaxy=float(P_noise_galaxy),
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
    Uses the Madau & Dickinson (2014) UV–SFR calibration
    (κ_UV = 1.15e-28, Chabrier IMF) via ``src.conversions``.
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
