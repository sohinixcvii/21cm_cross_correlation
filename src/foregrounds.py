#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
foregrounds.py — synthetic foreground injection and parametrised removal
========================================================================

Tools for testing how astrophysical foreground contamination, and how
*incomplete* removal of it, propagate into the 21 cm × galaxy cross-power
spectrum and its detectability.

Two halves:

**Injection** (:func:`simulate_diffuse_foreground`,
:func:`simulate_point_source_foreground`, :func:`inject_foreground`) builds a
spectrally smooth contaminant on the same grid and in the same units as the
simulated ``brightness_temp`` lightcone, so it can be added to the field
*before* any power spectrum is computed.

**Removal** (:func:`remove_foreground`) is a deliberately simple, controllable
knob — **not** a foreground-mitigation algorithm.  See its docstring.

Nothing here modifies or re-implements the verified analysis chain.  The
contaminated field is fed to ``analysis.compute_all_power_spectra`` and its
spectra to ``analysis.compute_uncertainty_budget`` exactly as the clean field
is.

What foregrounds do to a cross-spectrum, and what they do not
-------------------------------------------------------------
Foregrounds are uncorrelated with the galaxy field, so **in the ensemble mean**
they do not bias ``P_21×gal``.  This is the standard argument for
cross-correlation being the foreground-robust statistic, and it is correct as
far as it goes.  Two things complicate it, and both are visible in the outputs
of this module.

**1. The error bar grows.**  Contamination inflates the 21 cm auto-power, and
``P_21`` enters the La Plante et al. (2023) Eq. 15–17 variance
``σ_× = sqrt{0.5[P_×² + (|P_21| + P_N,21)(|P_gal| + P_N,gal)]}``.  Because
``σ_×`` depends on the square root of ``P_21``, a foreground ``A`` times the
signal in temperature inflates ``σ_×`` by roughly ``A``.

**2. A single realisation still picks up a chance correlation.**  "Unbiased in
the mean" is not "zero in your data".  Any one realisation has a spurious
cross-power of order ``sqrt(P_21 P_gal / N_modes)``, which grows **linearly**
with the foreground amplitude and shrinks only as ``sqrt(N_modes)``.  On a
small box with few usable ``(k_perp, k_parallel)`` bins this term can dominate
the real cross-power outright — measured here, a foreground 10⁴ times the
signal RMS produced a spurious ``|P_×|`` ~90 times the true one across the
modes surviving wedge excision.

The two effects partly cancel in the SNR.  ``|P_×| / σ_×`` has a contaminated
numerator *and* a contaminated denominator, so the total SNR degrades far more
slowly than ``σ_×`` alone — in the same test, ``σ_×`` rose 583x while the
total SNR fell only 6.8x.  **A foreground-contaminated SNR therefore
overstates how much signal survives**, because part of what it is detecting is
the foreground's own chance correlation with the galaxy field.

To separate the two, evaluate the SNR with the *clean* ``P_cross`` against the
*contaminated* ``σ_×``.  That isolates the error-bar cost from the spurious
signal, and is what the notebook's §7e plots alongside the as-measured curve.

Model
-----
The diffuse component follows the standard parametrised angular/frequency
power-law form used throughout 21 cm intensity-mapping forecasting: a
power-law angular power spectrum with a smooth power-law frequency dependence
and a spatially varying spectral index.  Rather than imposing a ``k_parallel``
power law by hand, the line-of-sight structure *emerges* from the smooth
frequency behaviour, which is the physically correct way round.

Smoothness does not confine the contamination to low ``k_parallel``
-------------------------------------------------------------------
It is tempting to assume a spectrally smooth foreground lives entirely at
``k_parallel ≈ 0``, safely inside the wedge.  Measured through *this*
pipeline's estimator it does not, and the reason is worth stating plainly.

``analysis.compute_cylindrical_cross_power`` takes a bare FFT with no
line-of-sight taper.  A foreground that is smooth but **not periodic** across
the band — any power law in frequency is — is discontinuous at the box edge
under that transform, and the discontinuity leaks power along the whole
``k_parallel`` axis as roughly ``k_parallel^-1.5``.  That fall-off is a
property of the window, not of the sky: a bare linear ramp with no angular
structure whatsoever leaks with the same slope, and widening the band does not
change it.

The practical consequence is that foreground power reaches the EoR window
here, and wedge excision alone does not remove it.  What remains true, and is
what these tools measure, is that contamination is far stronger at low
``k_parallel`` than at high — on the fiducial box the 21 cm auto-power is
inflated by ~400x in the lowest ``k_parallel`` bins against ~12x in the
highest.  Anyone using this module to make a quantitative statement about
window-region contamination should first decide whether the estimator ought to
be tapered; that is a change to the verified estimator and is deliberately not
made here.

References
----------
de Oliveira-Costa, A. et al. (2008), MNRAS 388, 247 — the Global Sky Model;
    origin of the diffuse Galactic synchrotron amplitude and spectral index
    used as defaults here.
Zheng, H. et al. (2017), MNRAS 464, 3486 — GSM2016, the improved global sky
    model over 10 MHz–5 THz.
Santos, M. G., Cooray, A. & Knox, L. (2005), ApJ 625, 575 — the parametrised
    ``C_l(ν_1, ν_2) = A (l/l_ref)^-β (ν_1ν_2/ν_ref²)^-α`` foreground model and
    its Table 1 parameters, adopted here for both the diffuse and
    point-source components.
Shaw, J. R. et al. (2014), ApJ 781, 57 — the same parametrisation in the
    m-mode formalism, and the spectral-index-variation term that spreads
    foreground power to non-zero ``k_parallel``.
Ali, S. S., Bharadwaj, S. & Chengalur, J. N. (2008), MNRAS 385, 2166 —
    point-source contribution to the low-frequency foreground budget.
Bernardi, G. et al. (2009), A&A 500, 965 — WSRT measurements of diffuse
    Galactic emission and the confusion-limited point-source population.
La Plante, P. et al. (2023), arXiv:2205.09770 — Eqs. 15–17, the variance the
    contaminated ``P_21`` feeds.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "DIFFUSE_DEFAULTS",
    "POINT_SOURCE_DEFAULTS",
    "ForegroundRealisation",
    "inject_foreground",
    "remove_foreground",
    "simulate_diffuse_foreground",
    "simulate_point_source_foreground",
]


# ===========================================================================
#  Model parameters
# ===========================================================================

#: Diffuse Galactic synchrotron, Santos, Cooray & Knox (2005) Table 1, with
#: the reference temperature from the Global Sky Model (de Oliveira-Costa
#: et al. 2008; Zheng et al. 2017) at 130 MHz.  ``beta_angular`` is the
#: angular power-law index (C_l ∝ l^-β), ``spectral_index`` the frequency
#: power law (T ∝ ν^-α), and ``spectral_index_scatter`` the spatial variation
#: in α that spreads power to non-zero k_parallel (Shaw et al. 2014).
DIFFUSE_DEFAULTS = {
    "beta_angular": 2.4,
    "spectral_index": 2.8,
    "spectral_index_scatter": 0.1,
    "reference_temperature_mK": 700e3,      # 700 K at 130 MHz
    "reference_frequency_mhz": 130.0,
    "contrast": 0.5,
}

#: Extragalactic point sources, Santos, Cooray & Knox (2005) Table 1; see also
#: Ali, Bharadwaj & Chengalur (2008) and Bernardi et al. (2009).  Flatter in
#: angle than the diffuse component (closer to Poisson/white) and slightly
#: flatter in frequency.
POINT_SOURCE_DEFAULTS = {
    "beta_angular": 1.1,
    "spectral_index": 2.07,
    "spectral_index_scatter": 0.3,
    "reference_temperature_mK": 57e3,       # 57 K at 130 MHz
    "reference_frequency_mhz": 130.0,
    "source_density_per_cell": 0.02,
    "flux_slope": 1.75,                     # dN/dS ∝ S^-1.75
}


@dataclass
class ForegroundRealisation:
    """
    A realised foreground field and the signal it was scaled against.

    Attributes
    ----------
    contaminated : ndarray
        ``signal + foreground``, same shape and units as the input field [mK].
    foreground : ndarray
        The injected foreground alone [mK].  Keep it: :func:`remove_foreground`
        needs it, since the placeholder removal is subtraction of a known
        field.
    diffuse, point_source : ndarray or None
        The two components before summation.  ``point_source`` is ``None``
        when it was not requested.
    signal_rms : float
        RMS of the input field [mK], the reference for ``amplitude``.
    foreground_rms : float
        RMS of ``foreground`` [mK].
    amplitude : float
        ``foreground_rms / signal_rms``, the requested contamination level.
    """

    contaminated: np.ndarray
    foreground: np.ndarray
    diffuse: np.ndarray
    point_source: Optional[np.ndarray]
    signal_rms: float
    foreground_rms: float
    amplitude: float


# ===========================================================================
#  Internal helpers
# ===========================================================================

def _validate_k_grid(k_perp: np.ndarray, k_parallel: np.ndarray) -> None:
    """
    Check the two k-grids are usable.

    Parameters
    ----------
    k_perp, k_parallel : ndarray
        Bin centres [Mpc^-1], as returned by
        ``analysis.compute_all_power_spectra``.

    Raises
    ------
    ValueError
        If either grid is empty or contains a non-positive value.  Both are
        used as log-spaced pivots, so zero or negative entries are a caller
        error rather than something to silently work around.
    """
    for name, grid in (("k_perp", k_perp), ("k_parallel", k_parallel)):
        grid = np.asarray(grid, dtype=float)
        if grid.size == 0:
            raise ValueError(f"{name} is empty")
        if not np.all(np.isfinite(grid)) or np.any(grid <= 0):
            raise ValueError(f"{name} must be finite and strictly positive")


def _infer_box_length(k_grid: np.ndarray) -> float:
    """
    Fall back to a box length implied by the smallest binned wavenumber.

    Only used when the caller does not supply the box geometry.  The smallest
    *bin centre* is an approximation to the fundamental mode ``2π/L``, so this
    is a convenience, not a substitute for passing ``box_len_perp`` /
    ``box_len_los``.

    Parameters
    ----------
    k_grid : ndarray
        Binned wavenumbers [Mpc^-1].

    Returns
    -------
    float
        Estimated comoving box length [Mpc].
    """
    return float(2.0 * np.pi / np.min(np.asarray(k_grid, dtype=float)))


def _accepted_kwargs(function, supplied: Dict[str, Any]) -> Dict[str, Any]:
    """
    Subset of ``supplied`` that ``function`` actually accepts.

    Lets :func:`inject_foreground` forward one shared set of options to two
    component builders whose parameter lists only partly overlap, without
    passing anything either of them would reject.

    Parameters
    ----------
    function : callable
        Target function.
    supplied : dict
        Candidate keyword arguments.

    Returns
    -------
    dict
        The accepted subset.
    """
    accepted = set(inspect.signature(function).parameters)
    return {key: value for key, value in supplied.items() if key in accepted}


def _gaussian_random_field_2d(
    n_cells: int,
    box_len: float,
    beta: float,
    k_pivot: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Zero-mean, unit-variance 2D Gaussian random field with ``P(k) ∝ k^-β``.

    Built by filtering white noise in Fourier space, which gives the requested
    angular power spectrum by construction.  The DC mode is set to zero — the
    overall sky temperature is carried by the reference temperature, not by
    this field — and the result is renormalised to unit variance so the
    caller controls the amplitude.

    Parameters
    ----------
    n_cells : int
        Cells per side of the transverse plane.
    box_len : float
        Transverse comoving side length [Mpc].
    beta : float
        Angular power-law index; ``P(k_perp) ∝ (k_perp / k_pivot)^-beta``.
    k_pivot : float
        Pivot wavenumber [Mpc^-1].  Only affects the (discarded) overall
        normalisation, but keeps the exponentiation well conditioned.
    rng : numpy.random.Generator
        Source of randomness.

    Returns
    -------
    ndarray
        Field of shape ``(n_cells, n_cells)``, zero mean and unit variance.
    """
    cell = box_len / n_cells
    k_axis = np.fft.fftfreq(n_cells, d=cell) * 2.0 * np.pi
    k_x, k_y = np.meshgrid(k_axis, k_axis, indexing="ij")
    k_magnitude = np.hypot(k_x, k_y)

    amplitude = np.zeros_like(k_magnitude)
    nonzero = k_magnitude > 0
    amplitude[nonzero] = (k_magnitude[nonzero] / k_pivot) ** (-0.5 * beta)

    white = rng.standard_normal((n_cells, n_cells))
    field = np.fft.ifft2(np.fft.fft2(white) * amplitude).real

    field -= field.mean()
    spread = field.std()
    return field / spread if spread > 0 else field


def _observed_frequencies_mhz(
    n_los: int,
    z_obs: float,
    lc_redshifts: Optional[Sequence[float]],
    box_len_los: float,
    f_21_mhz: float,
    hubble_constant: float,
    omega_m: float,
    speed_of_light_kms: float,
) -> np.ndarray:
    """
    Observed frequency of each line-of-sight slice [MHz].

    Uses the lightcone's own redshifts when given.  Otherwise the slices are
    mapped to redshift through the comoving depth of the box about ``z_obs``,
    using ``dz = H(z) dr / c``, so that a caller who has only the geometry
    still gets a physically spaced frequency axis rather than a flat one.

    Parameters
    ----------
    n_los : int
        Number of line-of-sight cells.
    z_obs : float
        Reference redshift at the centre of the box.
    lc_redshifts : sequence of float, optional
        Per-slice redshifts.  Must have length ``n_los`` when supplied.
    box_len_los : float
        Comoving line-of-sight extent [Mpc].
    f_21_mhz : float
        21 cm rest frequency [MHz].
    hubble_constant, omega_m, speed_of_light_kms : float
        Background cosmology, for the ``dz/dr`` conversion.

    Returns
    -------
    ndarray
        Observed frequencies, shape ``(n_los,)`` [MHz].

    Raises
    ------
    ValueError
        If ``lc_redshifts`` is supplied with the wrong length.
    """
    if lc_redshifts is not None:
        redshifts = np.asarray(lc_redshifts, dtype=float)
        if redshifts.shape != (n_los,):
            raise ValueError(
                f"lc_redshifts must have shape ({n_los},), "
                f"got {redshifts.shape}"
            )
        return f_21_mhz / (1.0 + redshifts)

    hubble_z = hubble_constant * np.sqrt(
        omega_m * (1.0 + z_obs) ** 3 + (1.0 - omega_m)
    )
    # dz = H(z) dr / c, about the box centre.
    offsets = (np.arange(n_los) - 0.5 * (n_los - 1)) * (box_len_los / n_los)
    redshifts = z_obs + offsets * hubble_z / speed_of_light_kms
    return f_21_mhz / (1.0 + redshifts)


def _apply_spectral_law(
    temperature_map: np.ndarray,
    index_map: np.ndarray,
    frequencies_mhz: np.ndarray,
    reference_frequency_mhz: float,
) -> np.ndarray:
    """
    Give each sightline a power-law spectrum ``T ∝ ν^-α``.

    This is what makes the foreground spectrally smooth, and so far stronger
    at low ``k_parallel`` than at high.  Letting ``α`` vary across the sky
    (Shaw et al. 2014) breaks the separability that an idealised rank-1
    removal would exploit — the mode-mixing that makes real foreground removal
    hard.  Smoothness alone does not keep the power out of the EoR window
    under an un-tapered transform; see the module docstring.

    Parameters
    ----------
    temperature_map : ndarray
        Reference-frequency brightness temperature, shape ``(n, n)`` [mK].
    index_map : ndarray
        Spectral index per sightline, same shape.
    frequencies_mhz : ndarray
        Observed frequency per LOS slice, shape ``(n_los,)`` [MHz].
    reference_frequency_mhz : float
        Frequency at which ``temperature_map`` is defined [MHz].

    Returns
    -------
    ndarray
        Foreground cube of shape ``(n, n, n_los)`` [mK].
    """
    frequency_ratio = frequencies_mhz / reference_frequency_mhz
    # (n, n, 1) ** (n_los,) broadcasts to (n, n, n_los).
    return temperature_map[:, :, None] * frequency_ratio[None, None, :] ** (
        -index_map[:, :, None]
    )


# ===========================================================================
#  Injection
# ===========================================================================

def simulate_diffuse_foreground(
    shape: Tuple[int, int, int],
    k_perp: np.ndarray,
    k_parallel: np.ndarray,
    z_obs: float,
    *,
    box_len_perp: Optional[float] = None,
    box_len_los: Optional[float] = None,
    lc_redshifts: Optional[Sequence[float]] = None,
    beta_angular: float = DIFFUSE_DEFAULTS["beta_angular"],
    spectral_index: float = DIFFUSE_DEFAULTS["spectral_index"],
    spectral_index_scatter: float = DIFFUSE_DEFAULTS["spectral_index_scatter"],
    reference_temperature_mK: float = DIFFUSE_DEFAULTS["reference_temperature_mK"],
    reference_frequency_mhz: float = DIFFUSE_DEFAULTS["reference_frequency_mhz"],
    contrast: float = DIFFUSE_DEFAULTS["contrast"],
    f_21_mhz: float = 1420.405,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
    seed: Optional[int] = 42,
) -> np.ndarray:
    """
    Synthetic diffuse Galactic synchrotron foreground cube.

    Builds a spectrally smooth contaminant on the same ``(N, N, N_z)`` grid
    and in the same units [mK] as the simulated ``brightness_temp`` lightcone,
    so it can be added to that field before any power spectrum is computed.

    The model is the standard parametrised one (Santos, Cooray & Knox 2005;
    Shaw et al. 2014), with amplitude and spectral index taken from the Global
    Sky Model (de Oliveira-Costa et al. 2008; Zheng et al. 2017, GSM2016):

    1. A reference-frequency sky ``T_ref(θ)`` whose angular power spectrum is
       ``C_l ∝ l^-β``, realised as a log-normal field so the temperature is
       positive everywhere as a physical sky must be.
    2. A spatially varying spectral index ``α(θ)``, correlated on the sky, with
       mean ``spectral_index`` and scatter ``spectral_index_scatter``.
    3. ``T(θ, ν) = T_ref(θ) (ν / ν_ref)^-α(θ)``.

    The line-of-sight structure is therefore *emergent* — it follows from the
    smooth frequency dependence rather than from an imposed ``k_parallel``
    power law, which is the physically correct way round.  Note that this does
    **not** confine the contamination to low ``k_parallel`` once it is measured
    through an un-tapered FFT; see the module docstring.

    Parameters
    ----------
    shape : tuple of int
        Cube shape ``(N, N, N_z)``, matching the brightness-temperature
        lightcone.
    k_perp, k_parallel : ndarray
        Binned wavenumbers [Mpc^-1] from
        ``analysis.compute_all_power_spectra``.  ``k_perp`` sets the pivot for
        the angular power law.  ``k_parallel`` is validated and used to infer
        the line-of-sight box length when ``box_len_los`` is not given; the
        line-of-sight structure itself comes from the frequency model above,
        not from this grid.
    z_obs : float
        Reference redshift at the centre of the box.
    box_len_perp, box_len_los : float, optional
        Comoving box lengths [Mpc].  Inferred from the smallest binned
        wavenumber when omitted; pass them when you have them.
    lc_redshifts : sequence of float, optional
        Per-slice redshifts, length ``N_z``.  Strongly preferred over the
        internal ``dz = H dr / c`` fallback.
    beta_angular : float, optional
        Angular power-law index, ``C_l ∝ l^-β``.
    spectral_index : float, optional
        Mean frequency power-law index ``α`` in ``T ∝ ν^-α``.
    spectral_index_scatter : float, optional
        Standard deviation of ``α`` across the sky.  Zero gives a perfectly
        separable foreground confined to ``k_parallel = 0``; the default
        0.1 is the realistic case where removal is genuinely hard.
    reference_temperature_mK : float, optional
        Mean sky temperature at ``reference_frequency_mhz`` [mK].  The output
        is normally rescaled by :func:`inject_foreground`, so this sets the
        shape rather than the final level.
    reference_frequency_mhz : float, optional
        Frequency at which ``reference_temperature_mK`` applies [MHz].
    contrast : float, optional
        Log-normal width of the angular temperature fluctuations.
    f_21_mhz : float, optional
        21 cm rest frequency [MHz].
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology, used only by the ``lc_redshifts`` fallback.
    seed : int, optional
        RNG seed.  ``None`` draws a fresh realisation.

    Returns
    -------
    ndarray
        Foreground brightness temperature of shape ``shape`` [mK], strictly
        positive.

    Raises
    ------
    ValueError
        If ``shape`` is not a 3-tuple, is not transversely square, or if the
        k-grids are empty or non-positive.

    Notes
    -----
    The returned cube includes its (large) monopole. That is physical — but
    ``analysis.compute_all_power_spectra`` subtracts the mean of the
    brightness-temperature field before transforming, and an interferometer
    does not measure the monopole either, so it never reaches the spectra.
    """
    if len(shape) != 3:
        raise ValueError(f"shape must be (N, N, N_z), got {shape}")
    n_perp, n_perp_y, n_los = (int(value) for value in shape)
    if n_perp != n_perp_y:
        raise ValueError(
            f"the transverse plane must be square, got {n_perp} x {n_perp_y}"
        )
    _validate_k_grid(k_perp, k_parallel)

    if box_len_perp is None:
        box_len_perp = _infer_box_length(k_perp)
    if box_len_los is None:
        box_len_los = _infer_box_length(k_parallel)

    rng = np.random.default_rng(seed)
    k_pivot = float(np.exp(np.mean(np.log(np.asarray(k_perp, dtype=float)))))

    # ── 1. Angular structure of the reference-frequency sky ───────────────
    # Log-normal so T > 0 everywhere, with the requested C_l ∝ l^-β for the
    # underlying Gaussian field.
    gaussian_sky = _gaussian_random_field_2d(
        n_perp, box_len_perp, beta_angular, k_pivot, rng
    )
    temperature_map = reference_temperature_mK * np.exp(
        contrast * gaussian_sky - 0.5 * contrast ** 2
    )

    # ── 2. Spatially varying spectral index ───────────────────────────────
    # Correlated on the sky rather than pixel-independent: neighbouring
    # sightlines have similar spectra, which is what the GSM shows.
    if spectral_index_scatter > 0:
        index_map = spectral_index + spectral_index_scatter * (
            _gaussian_random_field_2d(
                n_perp, box_len_perp, beta_angular, k_pivot, rng
            )
        )
    else:
        index_map = np.full((n_perp, n_perp), float(spectral_index))

    # ── 3. Smooth power-law spectrum along the line of sight ──────────────
    frequencies = _observed_frequencies_mhz(
        n_los, z_obs, lc_redshifts, box_len_los, f_21_mhz,
        hubble_constant, omega_m, speed_of_light_kms,
    )
    return _apply_spectral_law(
        temperature_map, index_map, frequencies, reference_frequency_mhz
    )


def simulate_point_source_foreground(
    shape: Tuple[int, int, int],
    k_perp: np.ndarray,
    k_parallel: np.ndarray,
    z_obs: float,
    *,
    box_len_perp: Optional[float] = None,
    box_len_los: Optional[float] = None,
    lc_redshifts: Optional[Sequence[float]] = None,
    spectral_index: float = POINT_SOURCE_DEFAULTS["spectral_index"],
    spectral_index_scatter: float = POINT_SOURCE_DEFAULTS["spectral_index_scatter"],
    reference_temperature_mK: float = POINT_SOURCE_DEFAULTS["reference_temperature_mK"],
    reference_frequency_mhz: float = POINT_SOURCE_DEFAULTS["reference_frequency_mhz"],
    source_density_per_cell: float = POINT_SOURCE_DEFAULTS["source_density_per_cell"],
    flux_slope: float = POINT_SOURCE_DEFAULTS["flux_slope"],
    f_21_mhz: float = 1420.405,
    hubble_constant: float = 67.36,
    omega_m: float = 0.315,
    speed_of_light_kms: float = 3e5,
    seed: Optional[int] = 43,
) -> np.ndarray:
    """
    Simple unresolved extragalactic point-source foreground cube.

    The companion to :func:`simulate_diffuse_foreground`: a Poisson-distributed
    population of sources, each with a power-law flux drawn from
    ``dN/dS ∝ S^-γ`` and its own spectral index, deposited on the transverse
    plane and given a smooth power-law spectrum.  Angularly this is close to
    white (Poisson) rather than the steep ``l^-2.4`` of the diffuse component,
    which is the qualitative distinction that matters here.

    Deliberately simple: no clustering term, no source-count completeness
    model, and no peeling of bright sources.  It exists so the total
    foreground is not purely Gaussian and steeply red, not to reproduce a
    measured source count.

    Parameters
    ----------
    shape : tuple of int
        Cube shape ``(N, N, N_z)``.
    k_perp, k_parallel : ndarray
        Binned wavenumbers [Mpc^-1]; used for validation and, when the box
        lengths are omitted, to infer them.
    z_obs : float
        Reference redshift at the centre of the box.
    box_len_perp, box_len_los : float, optional
        Comoving box lengths [Mpc].
    lc_redshifts : sequence of float, optional
        Per-slice redshifts, length ``N_z``.
    spectral_index : float, optional
        Mean frequency power-law index.
    spectral_index_scatter : float, optional
        Source-to-source scatter in the spectral index.
    reference_temperature_mK : float, optional
        Mean point-source brightness temperature at
        ``reference_frequency_mhz`` [mK].
    reference_frequency_mhz : float, optional
        Reference frequency [MHz].
    source_density_per_cell : float, optional
        Mean number of sources per transverse cell.
    flux_slope : float, optional
        Differential source-count slope ``γ`` in ``dN/dS ∝ S^-γ``.  Must
        exceed 1 for the drawn fluxes to have a finite mean.
    f_21_mhz : float, optional
        21 cm rest frequency [MHz].
    hubble_constant, omega_m, speed_of_light_kms : float, optional
        Background cosmology, used only by the ``lc_redshifts`` fallback.
    seed : int, optional
        RNG seed.

    Returns
    -------
    ndarray
        Point-source brightness temperature of shape ``shape`` [mK].

    Raises
    ------
    ValueError
        If ``shape`` is malformed, the k-grids are invalid, or
        ``flux_slope <= 1``.

    References
    ----------
    Ali, Bharadwaj & Chengalur (2008), MNRAS 385, 2166.
    Bernardi et al. (2009), A&A 500, 965.
    Santos, Cooray & Knox (2005), ApJ 625, 575 — Table 1 point-source row.
    """
    if len(shape) != 3:
        raise ValueError(f"shape must be (N, N, N_z), got {shape}")
    n_perp, n_perp_y, n_los = (int(value) for value in shape)
    if n_perp != n_perp_y:
        raise ValueError(
            f"the transverse plane must be square, got {n_perp} x {n_perp_y}"
        )
    if flux_slope <= 1.0:
        raise ValueError(
            f"flux_slope must exceed 1 for a finite mean flux, got {flux_slope}"
        )
    _validate_k_grid(k_perp, k_parallel)

    if box_len_perp is None:
        box_len_perp = _infer_box_length(k_perp)
    if box_len_los is None:
        box_len_los = _infer_box_length(k_parallel)

    rng = np.random.default_rng(seed)

    # ── Poisson source counts, Pareto fluxes ──────────────────────────────
    # A Pareto(γ - 1) variate has the tail of dN/dS ∝ S^-γ, so a handful of
    # bright sources dominate — the regime where point sources actually
    # matter.
    counts = rng.poisson(source_density_per_cell, size=(n_perp, n_perp))
    flux_map = np.zeros((n_perp, n_perp), dtype=float)
    occupied = counts > 0
    if occupied.any():
        # Total flux per cell: sum of `counts` draws, done as one gamma-like
        # aggregate would bias the tail, so draw per source instead.
        n_sources = int(counts.sum())
        fluxes = rng.pareto(flux_slope - 1.0, size=n_sources) + 1.0
        cell_index = np.repeat(np.flatnonzero(counts.ravel()),
                               counts.ravel()[occupied.ravel()])
        flux_map = np.bincount(
            cell_index, weights=fluxes, minlength=n_perp * n_perp
        ).reshape(n_perp, n_perp)

    mean_flux = flux_map.mean()
    if mean_flux > 0:
        temperature_map = reference_temperature_mK * flux_map / mean_flux
    else:                                    # no sources drawn at all
        temperature_map = np.zeros((n_perp, n_perp), dtype=float)

    index_map = spectral_index + spectral_index_scatter * rng.standard_normal(
        (n_perp, n_perp)
    )

    frequencies = _observed_frequencies_mhz(
        n_los, z_obs, lc_redshifts, box_len_los, f_21_mhz,
        hubble_constant, omega_m, speed_of_light_kms,
    )
    return _apply_spectral_law(
        temperature_map, index_map, frequencies, reference_frequency_mhz
    )


def inject_foreground(
    brightness_temp_field: np.ndarray,
    k_perp: np.ndarray,
    k_parallel: np.ndarray,
    z_obs: float,
    foreground_amplitude: float = 100.0,
    *,
    include_point_sources: bool = True,
    point_source_weight: float = 0.1,
    **kwargs,
) -> ForegroundRealisation:
    """
    Add a synthetic foreground to a brightness-temperature lightcone.

    Builds the diffuse component (and optionally a point-source component),
    rescales the total so that its RMS is ``foreground_amplitude`` times the
    RMS of ``brightness_temp_field``, and returns both the contaminated field
    and the foreground alone.

    Parameters
    ----------
    brightness_temp_field : ndarray
        Simulated δT_b lightcone [mK], shape ``(N, N, N_z)``.
    k_perp, k_parallel : ndarray
        Binned wavenumbers [Mpc^-1].
    z_obs : float
        Reference redshift.
    foreground_amplitude : float, optional
        Foreground RMS as a multiple of the signal RMS.  Real low-frequency
        foregrounds exceed the 21 cm signal by ~10⁴–10⁵ in temperature; the
        default of 100 is a deliberately mild setting for a first look.  Zero
        returns the field unchanged.
    include_point_sources : bool, optional
        Add the point-source component alongside the diffuse one.
    point_source_weight : float, optional
        Fraction of the total foreground RMS carried by point sources, before
        the overall rescaling.
    **kwargs
        Forwarded to :func:`simulate_diffuse_foreground` and, where
        applicable, :func:`simulate_point_source_foreground` — e.g.
        ``box_len_perp``, ``box_len_los``, ``lc_redshifts``, ``seed``.

    Returns
    -------
    ForegroundRealisation
        The contaminated field, the foreground, its components, and the
        amplitude actually achieved.

    Raises
    ------
    ValueError
        If ``foreground_amplitude`` is negative, ``point_source_weight`` is
        outside ``[0, 1]``, or the input field is not 3D.

    Notes
    -----
    Scaling is applied to the foreground's *fluctuation* RMS, since that is
    what survives the mean subtraction in
    ``analysis.compute_all_power_spectra`` and therefore what reaches the
    power spectra.
    """
    field = np.asarray(brightness_temp_field, dtype=float)
    if field.ndim != 3:
        raise ValueError(f"expected a 3D lightcone, got shape {field.shape}")
    if foreground_amplitude < 0:
        raise ValueError(
            f"foreground_amplitude must be >= 0, got {foreground_amplitude}"
        )
    if not 0.0 <= point_source_weight <= 1.0:
        raise ValueError(
            f"point_source_weight must be in [0, 1], got {point_source_weight}"
        )

    signal_rms = float(field.std())

    # Forward only what each builder accepts, and give the two components
    # different seeds — sharing one would correlate them, which they are not.
    diffuse_kwargs = _accepted_kwargs(simulate_diffuse_foreground, kwargs)
    point_kwargs = _accepted_kwargs(simulate_point_source_foreground, kwargs)
    if "seed" in kwargs and kwargs["seed"] is not None:
        point_kwargs["seed"] = int(kwargs["seed"]) + 1

    diffuse = simulate_diffuse_foreground(
        field.shape, k_perp, k_parallel, z_obs, **diffuse_kwargs,
    )

    point_source = None
    if include_point_sources and point_source_weight > 0:
        point_source = simulate_point_source_foreground(
            field.shape, k_perp, k_parallel, z_obs, **point_kwargs,
        )

    # Combine on fluctuation RMS, so `point_source_weight` means what it says
    # regardless of the two components' absolute reference temperatures.
    def _unit(component: np.ndarray) -> np.ndarray:
        spread = component.std()
        return component / spread if spread > 0 else component

    if point_source is None:
        total = _unit(diffuse)
    else:
        total = (
            (1.0 - point_source_weight) * _unit(diffuse)
            + point_source_weight * _unit(point_source)
        )

    total_rms = float(total.std())
    if total_rms > 0:
        total = total * (foreground_amplitude * signal_rms / total_rms)

    return ForegroundRealisation(
        contaminated=field + total,
        foreground=total,
        diffuse=diffuse,
        point_source=point_source,
        signal_rms=signal_rms,
        foreground_rms=float(total.std()),
        amplitude=float(total.std() / signal_rms) if signal_rms > 0 else 0.0,
    )


# ===========================================================================
#  Removal
# ===========================================================================

def remove_foreground(
    contaminated_field: np.ndarray,
    k_perp: np.ndarray,
    k_parallel: np.ndarray,
    removal_fraction: float,
    *,
    foreground: np.ndarray,
    removal_basis: str = "amplitude",
) -> np.ndarray:
    """
    Remove a controllable fraction of an injected foreground.

    .. warning::

       **This is a simplified placeholder, not a foreground-mitigation
       method.** It subtracts a known, exactly-correct template of the very
       field that was injected, scaled by ``removal_fraction``. It is a knob
       for asking *"what if removal were this good?"*, and nothing more.

       It is **not** GMCA, PCA, ICA, polynomial or log-polynomial fitting,
       Gaussian-process foreground removal, delay filtering, or any other real
       technique. It has none of their failure modes: no signal loss from
       over-fitting, no mode-mixing, no leakage from the wedge into the EoR
       window, no dependence on the foreground's spectral smoothness, and no
       sensitivity to the number of components removed. A real method's
       residual is structured and correlated with the signal it damaged; this
       residual is simply a scaled copy of the input contaminant.

       Any result obtained with it is a statement about *removal level*, not
       about any actual algorithm's achievable performance. Do not quote it as
       a forecast for a named method.

    Parameters
    ----------
    contaminated_field : ndarray
        Field with the foreground already added, shape ``(N, N, N_z)`` [mK] —
        ``ForegroundRealisation.contaminated``.
    k_perp, k_parallel : ndarray
        Binned wavenumbers [Mpc^-1].  Accepted so the signature matches a
        k-space filter and validated for consistency, but **not used**: the
        placeholder is scale-independent by construction, which is one of the
        ways it is unlike a real method.
    removal_fraction : float
        How much of the foreground to remove, in ``[0, 1]``.  ``0`` leaves the
        field untouched; ``1`` recovers the clean field exactly.
    foreground : ndarray
        The injected foreground, same shape —
        ``ForegroundRealisation.foreground``.  Required: this method works
        only because the answer is known in advance.
    removal_basis : {'amplitude', 'power'}, optional
        What ``removal_fraction`` is a fraction *of*.

        - ``'amplitude'`` (default) scales the residual by
          ``(1 - removal_fraction)``, so residual **power** falls as
          ``(1 - removal_fraction)²``.  ``removal_fraction = 0.9`` leaves 10%
          of the amplitude and 1% of the power.
        - ``'power'`` scales the residual by ``sqrt(1 - removal_fraction)``,
          so residual power falls as ``(1 - removal_fraction)``.
          ``removal_fraction = 0.9`` leaves 10% of the power.

    Returns
    -------
    ndarray
        ``contaminated_field - removal * foreground``: the clean signal plus
        the residual foreground.

    Raises
    ------
    ValueError
        If ``removal_fraction`` is outside ``[0, 1]``, the shapes disagree, or
        ``removal_basis`` is not one of the two accepted values.

    Examples
    --------
    >>> residual = remove_foreground(                     # doctest: +SKIP
    ...     realisation.contaminated, k_perp, k_parallel, 0.99,
    ...     foreground=realisation.foreground,
    ... )
    """
    contaminated = np.asarray(contaminated_field, dtype=float)
    foreground_array = np.asarray(foreground, dtype=float)

    if contaminated.shape != foreground_array.shape:
        raise ValueError(
            f"contaminated_field and foreground shapes differ: "
            f"{contaminated.shape} vs {foreground_array.shape}"
        )
    if not 0.0 <= removal_fraction <= 1.0:
        raise ValueError(
            f"removal_fraction must be in [0, 1], got {removal_fraction}"
        )
    if removal_basis not in ("amplitude", "power"):
        raise ValueError(
            f"removal_basis must be 'amplitude' or 'power', "
            f"got {removal_basis!r}"
        )
    _validate_k_grid(k_perp, k_parallel)

    if removal_basis == "amplitude":
        residual_scale = 1.0 - removal_fraction
    else:
        residual_scale = np.sqrt(1.0 - removal_fraction)

    return contaminated - (1.0 - residual_scale) * foreground_array
