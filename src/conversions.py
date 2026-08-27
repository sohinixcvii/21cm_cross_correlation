# ===================================================================
# Useful cosmological conversions for high-redshift galaxy surveys.
# ===================================================================

# Currently includes:
# - Muv <-> Luv conversions
# - Survey volume <-> survey area conversions
# - Area conversions between square degrees and steradians
# - Volume calculation from area and redshift interval
# - Survey footprint -> 21cmFAST simulation box geometry

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
import astropy.units as u
from astropy.cosmology import Planck18
from scipy.integrate import simpson

def Muv_to_Luv(Muv):
    """
    Convert an absolute UV AB magnitude to a monochromatic UV luminosity.

    Parameters
    ----------
    Muv : float or ndarray
        Absolute UV magnitude in the AB magnitude system.
        Typical values for high-redshift galaxies are between
        approximately -15 and -25.

    Returns
    -------
    Luv : float or ndarray
        Monochromatic UV luminosity L_nu in units of

            erg s^-1 Hz^-1

    Notes
    -----
    Uses the standard AB magnitude relation

        M_AB = -2.5 log10(L_nu) + 51.60

    which assumes luminosity is measured at a distance of 10 pc.

    References
    ----------
    Oke & Gunn (1983)
    Madau & Dickinson (2014)

    Examples
    --------
    >>> Muv_to_Luv(-22.66)
    5.06e29
    """
    return 10 ** ((51.60 - Muv) / 2.5)


def Luv_to_Muv(Luv):
    """
    Convert a monochromatic UV luminosity to an absolute UV AB magnitude.

    Parameters
    ----------
    Luv : float or ndarray
        Monochromatic UV luminosity in units of

            erg s^-1 Hz^-1

    Returns
    -------
    Muv : float or ndarray
        Absolute UV magnitude in the AB magnitude system.

    Notes
    -----
    Inverts the relation

        M_AB = -2.5 log10(L_nu) + 51.60

    References
    ----------
    Oke & Gunn (1983)
    Madau & Dickinson (2014)

    Examples
    --------
    >>> Luv_to_Muv(5.06e29)
    -22.66
    """
    return 51.60 - 2.5 * np.log10(Luv)

# UV-SFR conversion factor:  SFR [M_sun yr^-1] = kappa_UV * L_UV [erg/s/Hz]
#
# ADOPTED: Fisher et al. (2026), MNRAS in press, arXiv:2511.10741, Eq. 12 --
# "REBELS-IFU: Steeply rising star formation histories and the importance of
# dust obscuration in massive z=7 galaxies revealed by multi-wavelength
# observations" (Jodrell Bank Centre for Astrophysics).
#
#     kappa_UV = (2.7 +/- 0.9) x 10^-29  M_sun/yr / (erg/s/Hz)
#
# Derived from JWST NIRSpec-fit *non-parametric* star formation histories of
# 12 massive (log(M*/M_sun) = 9-10) Lyman-break galaxies at z = 6.5-7.7 from
# the REBELS-IFU sample, whose SFHs rise steeply rather than being constant.
#
# ---------------------------------------------------------------------------
# CAVEAT 1 -- DEFINITIONAL, and the most important.  This kappa_UV is
# calibrated specifically to recover SFR_100Myr -- the SFR averaged over the
# prior 100 Myr -- from a population with RISING SFHs.  It is NOT an
# instantaneous / constant-SFH calibration.  Any downstream code that assumes
# a constant or instantaneous SFR interpretation now means something DIFFERENT
# by "SFR", not merely a different number.  This is a correctness question,
# not a rescaling.  Two places in this pipeline are affected:
#   * analysis.euclid_sfr_window() / select_euclid_halos() convert the Euclid
#     M_UV window into an SFR window and select halos on 21cmFAST's `halo_sfr`.
#     That SFR is M_star / t_sf with t_sf = t_STAR * t_H(z) = 570.3 Myr at
#     z = 7 (Park et al. 2019, Eq. 3) -- a 5.7x LONGER averaging window than
#     the 100 Myr this kappa_UV recovers, under a constant-SFH-like
#     prescription rather than a rising one.  Which halos count as
#     "Euclid-detected" therefore changes meaning, and b_g / n_bar / the
#     galaxy field follow.
#   * docs/halo_catalogue_reference.md states the pipeline's L_UV-SFR relation
#     assumes "Continuous star formation".  This value contradicts that
#     documented assumption directly.
# NOT resolved here -- flagged for human review.
#
# CAVEAT 2 -- UNCERTAINTY.  +/- 0.9 on 2.7 is a ~33% systematic uncertainty,
# substantially larger than a typical calibration constant.  Any forecast
# quoting a number derived through this factor should carry that spread, not
# the central value alone.
#
# CAVEAT 2b -- IMF, UNVERIFIED.  Fisher et al.'s Table 2 fiducial
# (constant-SFH) value uses Chabrier (2003), and the Eq. 12 rising-SFH value
# comes from the same BAGPIPES SED-fitting setup, so it is almost certainly
# also Chabrier -- but that was NOT verified.  It should be checked against
# the IMF assumed by the Park et al. (2019) SFE model used elsewhere in this
# pipeline (F_STAR10, ALPHA_STAR, t_STAR).  **This codebase does not record
# Park et al. (2019)'s IMF anywhere**, so consistency can be neither confirmed
# nor denied from the repository alone.  NOT resolved here.
#
# CAVEAT 3 -- POPULATION SPECIFICITY.  Calibrated on massive
# (log(M*/M_sun) = 9-10), rest-UV-bright REBELS galaxies.  Fisher et al.
# Sect. 6.4 note that "conversion factors applicable to the more abundant
# lower-mass galaxies at z~7... may not be quite as low."  If this pipeline's
# simulated population spans a broader or fainter mass range than the
# REBELS-IFU calibration sample, applying this value is an EXTRAPOLATION, not
# a like-for-like match.  NOT resolved here.
# ---------------------------------------------------------------------------
#
# Candidates considered and NOT adopted, recorded so the choice is visible:
#   * 7.2e-29  -- Chabrier (2003) IMF, constant SFH (Fisher et al. Table 2
#                 fiducial).
#   * 1.15e-28 -- the previous value: raw Salpeter IMF, constant SFH, Madau &
#                 Dickinson (2014), ARA&A 52, 415, rest-frame ~1500 A.
#                 Dhandha et al. (2026), MNRAS in press, arXiv:2508.13761,
#                 Eq. 18 independently confirms this 1.15e-28 -- that
#                 corroboration stands, but now corroborates a SUPERSEDED
#                 value.
#
# NOTE -- unresolved IMF labelling discrepancy.  This repository labelled
# 1.15e-28 as "Chabrier (2003)" in several places, while
# docs/halo_catalogue_reference.md calls the same relation "Salpeter-like" and
# Fisher et al. treat 1.15e-28 as the raw Salpeter value with 7.2e-29 as the
# Chabrier equivalent.  The repo labelling was therefore inconsistent and
# probably wrong.  Recorded, not silently corrected.  See
# NUMBERS_AND_SOURCES.md section 2.
#
# The name `_KAPPA_UV_MADAU14` is kept for now so no call site breaks; it is
# no longer accurate and is a rename candidate.
_KAPPA_UV_MADAU14 = 2.7e-29

# Critical density of the Universe today, in units of  M_sun Mpc^-3 h^-2.
#   rho_crit,0 = 3 H_0^2 / (8 pi G)  with  H_0 = 100 h km/s/Mpc
_RHO_CRIT_0_MSUN_MPC3_H2 = 2.77536627e11


def Luv_to_sfr(Luv, kappa_uv=_KAPPA_UV_MADAU14):
    """
    Convert a monochromatic UV luminosity to a star-formation rate.

    Parameters
    ----------
    Luv : float or ndarray
        Monochromatic UV luminosity in units of

            erg s^-1 Hz^-1

    kappa_uv : float, optional
        UV–SFR conversion factor in units of

            M_sun yr^-1 / (erg s^-1 Hz^-1)

        Default is 2.7e-29 from Fisher et al. (2026), arXiv:2511.10741,
        Eq. 12 -- a rising-SFH calibration recovering SFR_100Myr, with a
        ~33% systematic uncertainty (+/- 0.9).  See the module-level
        definition for three caveats that affect correctness, not just
        magnitude.

    Returns
    -------
    sfr : float or ndarray
        Star-formation rate in units of

            M_sun yr^-1

    Notes
    -----
    Uses the linear relation

        SFR = kappa_UV * L_UV

    References
    ----------
    Madau & Dickinson (2014), ARAA 52, 415
    Kennicutt & Evans (2012), ARAA 50, 531

    Examples
    --------
    >>> Luv_to_sfr(2.754e27)
    0.317
    """
    return kappa_uv * Luv


def sfr_to_Luv(sfr, kappa_uv=_KAPPA_UV_MADAU14):
    """
    Convert a star-formation rate to a monochromatic UV luminosity.

    Parameters
    ----------
    sfr : float or ndarray
        Star-formation rate in units of

            M_sun yr^-1

    kappa_uv : float, optional
        UV–SFR conversion factor in units of

            M_sun yr^-1 / (erg s^-1 Hz^-1)

        Default is 2.7e-29 from Fisher et al. (2026), arXiv:2511.10741,
        Eq. 12 -- a rising-SFH calibration recovering SFR_100Myr, with a
        ~33% systematic uncertainty (+/- 0.9).  See the module-level
        definition for three caveats that affect correctness, not just
        magnitude.

    Returns
    -------
    Luv : float or ndarray
        Monochromatic UV luminosity in units of

            erg s^-1 Hz^-1

    Notes
    -----
    Inverts the relation

        SFR = kappa_UV * L_UV

    References
    ----------
    Madau & Dickinson (2014), ARAA 52, 415

    Examples
    --------
    >>> sfr_to_Luv(1.0)
    8.696e27
    """
    return sfr / kappa_uv


def sheth_tormen_bias(nu_sq, delta_c=1.686, a=0.707, p=0.3):
    """
    Sheth-Tormen halo bias as a function of squared peak height.

    Parameters
    ----------
    nu_sq : float or ndarray
        Squared peak height  (δ_c / σ)²  as returned by the ``hmf`` package
        (``MassFunction.nu``).  This is **not** δ_c/σ itself.

    delta_c : float, optional
        Linear collapse threshold. Default 1.686.

    a : float, optional
        Sheth-Tormen parameter. Default 0.707.

    p : float, optional
        Sheth-Tormen parameter. Default 0.3.

    Returns
    -------
    bias : float or ndarray
        Eulerian linear halo bias  b_h(ν, z).

    Notes
    -----
    Implements the Sheth-Tormen (1999) bias formula

        b(ν̃) = 1 + (a ν̃ - 1) / δ_c + 2p / (δ_c (1 + (a ν̃)^p))

    where  ν̃ = (δ_c / σ)²  is the squared peak height.

    The ``hmf`` package (v3+) stores ``MassFunction.nu`` as this squared
    quantity, consistent with the original Sheth & Tormen (1999) notation.
    Passing ``MassFunction.nu`` directly to this function is therefore
    correct.

    References
    ----------
    Sheth & Tormen (1999), MNRAS 308, 119
    Tinker et al. (2010), ApJ 724, 878

    Examples
    --------
    >>> sheth_tormen_bias(5.0)   # nu_sq = 5 -> bias ~ 2.6
    2.63
    """
    a_nu = a * nu_sq
    return 1.0 + (a_nu - 1.0) / delta_c + (2.0 * p) / (delta_c * (1.0 + a_nu**p))


def mean_matter_density(
    omega_m: float = Planck18.Om0,
    hubble_constant: float = Planck18.H0.value,
) -> float:
    """
    Mean comoving matter density of the Universe.

    Parameters
    ----------
    omega_m : float, optional
        Matter density parameter Ω_m at z = 0. Defaults to Planck18.

    hubble_constant : float, optional
        Hubble constant H_0 in km s^-1 Mpc^-1. Defaults to Planck18.

    Returns
    -------
    rho_m : float
        Comoving mean matter density  ρ̄_m  in  M_sun Mpc^-3.

    Notes
    -----
    The comoving matter density is constant in time,

        ρ̄_m = Ω_m ρ_crit,0 = Ω_m · 3 H_0² / (8 π G),

    so a comoving cell of fixed size always encloses the same mean mass.

    Examples
    --------
    >>> f"{mean_matter_density():.3e}"
    '3.934e+10'
    """
    h_squared = (hubble_constant / 100.0) ** 2
    return omega_m * _RHO_CRIT_0_MSUN_MPC3_H2 * h_squared


def cell_mass(
    cell_size_mpc: float,
    omega_m: float = Planck18.Om0,
    hubble_constant: float = Planck18.H0.value,
) -> float:
    """
    Mean matter mass enclosed by one cubic comoving grid cell.

    This is the mass resolution of a grid-based (Eulerian) simulation: the
    smallest mass element the density field can represent.

    Parameters
    ----------
    cell_size_mpc : float
        Comoving side length of a single cubic cell, in Mpc.

    omega_m : float, optional
        Matter density parameter Ω_m at z = 0. Defaults to Planck18.

    hubble_constant : float, optional
        Hubble constant H_0 in km s^-1 Mpc^-1. Defaults to Planck18.

    Returns
    -------
    mass : float
        Mean enclosed matter mass  M_cell = ρ̄_m · L_cell³  in  M_sun.

    Examples
    --------
    >>> f"{cell_mass(2.0):.3e}"   # 256 Mpc box on a 128³ grid
    '3.147e+11'
    """
    return mean_matter_density(omega_m, hubble_constant) * cell_size_mpc**3


def sfr_to_Muv(sfr, kappa_uv=_KAPPA_UV_MADAU14):
    """
    Convert a star-formation rate directly to an absolute UV AB magnitude.

    Parameters
    ----------
    sfr : float or ndarray
        Star-formation rate in units of

            M_sun yr^-1

    kappa_uv : float, optional
        UV–SFR conversion factor in units of

            M_sun yr^-1 / (erg s^-1 Hz^-1)

        Default is 2.7e-29 from Fisher et al. (2026), arXiv:2511.10741,
        Eq. 12.  See the module-level definition for its caveats.

    Returns
    -------
    Muv : float or ndarray
        Absolute UV magnitude in the AB magnitude system.

    Notes
    -----
    Chains the two-step conversion

        SFR → L_UV = SFR / kappa_UV
        L_UV → M_UV = 51.60 - 2.5 log10(L_UV)

    References
    ----------
    Madau & Dickinson (2014), ARAA 52, 415
    Oke & Gunn (1983)

    Examples
    --------
    >>> sfr_to_Muv(1.0)
    -18.25
    """
    return Luv_to_Muv(sfr_to_Luv(sfr, kappa_uv=kappa_uv))


def mab_to_Muv(mab, z, cosmo: Optional[Any] = None):
    """
    Convert an apparent AB magnitude to an absolute UV AB magnitude.

    ``M_UV = m_AB - DM(z) + 2.5 log10(1 + z)``, where
    ``DM(z) = 5 log10(D_L / 10 pc)`` is the standard distance modulus and the
    ``+2.5 log10(1 + z)`` term is the usual bandwidth-compression correction
    that takes the observed band to the rest-frame band.

    Parameters
    ----------
    mab : float or ndarray
        Apparent magnitude in the AB system.
    z : float
        Redshift of the source.
    cosmo : astropy.cosmology.FLRW, optional
        Cosmology for the luminosity distance.  Defaults to ``Planck18``,
        matching the ``cosmo=None`` convention used elsewhere in this module.

    Returns
    -------
    Muv : float or ndarray
        Absolute UV magnitude in the AB magnitude system.

    Notes
    -----
    Matches Euclid Collaboration: Allen et al. (2026), A&A 711, A25, Eq. 15.
    Their additional ``Delta M_loss`` term is deliberately omitted: it
    corrects for Lyman-alpha forest absorption entering the H_E band, which
    only applies at z > 10 and so is irrelevant at this pipeline's z ~ 7.

    The conversion is only weakly cosmology-dependent — Planck18 and the
    Allen et al. cosmology (Omega_m = 0.27, H_0 = 70) differ by 0.033 mag at
    z = 7.

    References
    ----------
    Euclid Collaboration: Allen et al. (2026), A&A 711, A25, Eq. 15
    Oke & Gunn (1983)

    Examples
    --------
    >>> round(float(mab_to_Muv(25.0, 7.0)), 2)
    -21.98
    """
    cosmology = Planck18 if cosmo is None else cosmo
    distance_modulus = 5.0 * np.log10(
        cosmology.luminosity_distance(z).to(u.pc).value / 10.0
    )
    return mab - distance_modulus + 2.5 * np.log10(1.0 + z)


def Muv_to_mab(Muv, z, cosmo: Optional[Any] = None):
    """
    Convert an absolute UV AB magnitude to an apparent AB magnitude.

    Exact inverse of :func:`mab_to_Muv`:
    ``m_AB = M_UV + DM(z) - 2.5 log10(1 + z)``.

    Parameters
    ----------
    Muv : float or ndarray
        Absolute UV magnitude in the AB magnitude system.
    z : float
        Redshift of the source.
    cosmo : astropy.cosmology.FLRW, optional
        Cosmology for the luminosity distance.  Defaults to ``Planck18``.

    Returns
    -------
    mab : float or ndarray
        Apparent magnitude in the AB system.

    See Also
    --------
    mab_to_Muv : The forward conversion, with the full derivation and
        references.

    Examples
    --------
    >>> round(float(Muv_to_mab(-21.983, 7.0)), 2)
    25.0
    """
    cosmology = Planck18 if cosmo is None else cosmo
    distance_modulus = 5.0 * np.log10(
        cosmology.luminosity_distance(z).to(u.pc).value / 10.0
    )
    return Muv + distance_modulus - 2.5 * np.log10(1.0 + z)


def survey_area_from_volume(
    volume_mpc3,
    z_min,
    z_max,
    cosmo=None,
):
    """
    Convert a comoving survey volume into an equivalent sky area.

    Parameters
    ----------
    volume_mpc3 : float
        Total comoving survey volume in units of

            Mpc^3

    z_min : float
        Lower redshift limit of the survey.

    z_max : float
        Upper redshift limit of the survey.

    cosmo : astropy.cosmology.FLRW, optional
        Cosmology used for the conversion.
        Defaults to Planck18 cosmology.

    Returns
    -------
    area_deg2 : astropy.units.Quantity
        Survey area in square degrees.

    Notes
    -----
    A survey volume alone does not uniquely determine an area.
    The redshift interval must also be supplied.

    The conversion is performed using

        V = Ω ∫ (dV/dz/dΩ) dz

    where

        V           = comoving volume
        Ω           = solid angle
        dV/dz/dΩ    = differential comoving volume

    References
    ----------
    Hogg (1999), "Distance Measures in Cosmology"
    Astropy Cosmology Documentation

    Examples
    --------
    >>> survey_area_from_volume(
    ...     volume_mpc3=1.02e8,
    ...     z_min=6,
    ...     z_max=7
    ... )
    <Quantity ... deg2>
    """
    cosmo = Planck18 if cosmo is None else cosmo

    z_grid = np.linspace(z_min, z_max, 1000)

    dV_dz_dOmega = cosmo.differential_comoving_volume(z_grid)

    V_per_sr = (
        simpson(dV_dz_dOmega.value, x=z_grid)
        * u.Mpc**3
        / u.sr
    )

    omega = volume_mpc3 * u.Mpc**3 / V_per_sr

    return omega.to(u.deg**2)

def area_deg2_to_steradians(area_deg2):
    """
    Convert survey area from square degrees to steradians.

    Parameters
    ----------
    area_deg2 : float
        Survey area in square degrees.

    Returns
    -------
    area_sr : astropy.units.Quantity
        Survey area in steradians.
    """
    import astropy.units as u

    return (area_deg2 * u.deg**2).to(u.sr)


def volume_from_area(
    area_deg2,
    z_min,
    z_max,
    cosmo=None,
    n_z=1000,
):
    """
    Compute the comoving survey volume for a sky area and redshift interval.

    Parameters
    ----------
    area_deg2 : float
        Survey area in square degrees.

    z_min : float
        Lower redshift limit.

    z_max : float
        Upper redshift limit.

    cosmo : astropy.cosmology.FLRW, optional
        Cosmology used for the volume calculation.
        Defaults to Planck18.

    n_z : int, optional
        Number of redshift grid points used for numerical integration.
        Default is 1000.

    Returns
    -------
    volume : astropy.units.Quantity
        Comoving survey volume in Mpc^3.

    Notes
    -----
    Uses

        V = Omega * integral_zmin^zmax (dV / dz / dOmega) dz

    where Omega is the survey solid angle in steradians.
    """
    import numpy as np
    import astropy.units as u
    from astropy.cosmology import Planck18
    from scipy.integrate import simpson

    cosmo = Planck18 if cosmo is None else cosmo

    omega = area_deg2_to_steradians(area_deg2)

    z_grid = np.linspace(z_min, z_max, n_z)
    dV_dz_dOmega = cosmo.differential_comoving_volume(z_grid)

    volume_per_sr = (
        simpson(dV_dz_dOmega.value, x=z_grid)
        * u.Mpc**3
        / u.sr
    )

    volume = omega * volume_per_sr

    return volume.to(u.Mpc**3)

# ===================================================================
# Survey footprint -> 21cmFAST simulation box geometry
# ===================================================================

# Default target comoving cell size [Mpc].  The production grid before the
# footprint-driven sizing was BOX_LEN = 256 Mpc on HII_DIM = 128, i.e. exactly
# 2.0 Mpc per cell (M_cell = 3.17e11 M_sun).  Deriving HII_DIM from this
# target instead of pinning it preserves the mass resolution when the box
# grows to cover the survey footprint.
_TARGET_CELL_SIZE_MPC = 2.0

# 21cmFAST convention: the high-resolution initial-conditions grid is a fixed
# multiple of the low-resolution ionisation/21 cm grid.
_DIM_PER_HII_DIM = 3


@dataclass
class SimulationBox:
    """
    21cmFAST box geometry implied by a survey footprint.

    All lengths are comoving and in Mpc, matching 21cmFAST's ``BOX_LEN``
    convention (``simulation_options["BOX_LEN"]`` is comoving Mpc, *not*
    Mpc/h).

    Attributes
    ----------
    box_len : float
        Transverse comoving side length [Mpc] — pass straight to
        ``BOX_LEN``.  21cmFAST boxes are cubic, so this is the side of the
        coeval cube from which the lightcone is drawn.
    hii_dim : int
        Cells per side of the low-resolution ionisation / 21 cm grid.
    dim : int
        Cells per side of the high-resolution initial-conditions grid,
        ``_DIM_PER_HII_DIM * hii_dim``.
    cell_size : float
        Comoving size of one ``hii_dim`` cell [Mpc], ``box_len / hii_dim``.
    los_depth : float
        Comoving line-of-sight depth implied by ``delta_z`` [Mpc].  This is
        *not* a 21cmFAST argument: the LOS extent of a lightcone comes from
        the redshift range handed to ``RectilinearLightconer``, so use
        ``z_min`` / ``z_max`` below.  Reported so the LOS geometry can be
        checked against the box.
    z_min, z_max : float
        Redshift range spanning ``delta_z`` about the central redshift,
        ``z_central -/+ delta_z / 2``.
    transverse_area : float
        Comoving transverse area of the footprint [Mpc^2].
    solid_angle_sr : float
        Survey solid angle [sr].
    n_los_tiles : float
        ``los_depth / box_len``.  Greater than 1 means 21cmFAST must tile the
        coeval box more than once along the line of sight to fill the
        lightcone, which repeats structure on the box scale.

    Notes
    -----
    ``simulation_options`` returns the dict of 21cmFAST keyword arguments
    directly, so the caller never has to restate the mapping.
    """

    box_len: float
    hii_dim: int
    dim: int
    cell_size: float
    los_depth: float
    z_min: float
    z_max: float
    transverse_area: float
    solid_angle_sr: float
    n_los_tiles: float

    @property
    def simulation_options(self) -> dict[str, Any]:
        """
        The ``simulation_options`` mapping for ``InputParameters.clone()``.

        Returns
        -------
        dict
            ``{"HII_DIM": ..., "BOX_LEN": ..., "DIM": ...}``, ready for
            ``p21c.InputParameters.from_template(...).clone(
            simulation_options=box.simulation_options)``.
        """
        return {
            "HII_DIM": self.hii_dim,
            "BOX_LEN": self.box_len,
            "DIM": self.dim,
        }


def _next_power_of_two(n: int) -> int:
    """
    Smallest power of two greater than or equal to ``n``.

    Parameters
    ----------
    n : int
        Positive integer.

    Returns
    -------
    int
        The next power of two at or above ``n``.
    """
    if n < 1:
        raise ValueError(f"n must be a positive integer, got {n}.")
    return 1 << (n - 1).bit_length()


def survey_area_to_box_size(
    area_deg2: float,
    z_central: float,
    delta_z: float,
    cosmo: Optional[Any] = None,
    target_cell_size_mpc: float = _TARGET_CELL_SIZE_MPC,
    hii_dim: Optional[int] = None,
    snap_hii_dim_to_power_of_two: bool = True,
) -> SimulationBox:
    """
    Convert a survey footprint into a 21cmFAST simulation box geometry.

    Turns an angular survey area and a redshift interval into the comoving
    box size 21cmFAST needs, so ``BOX_LEN`` is traceable to the survey being
    forecast rather than chosen by hand.

    Parameters
    ----------
    area_deg2 : float
        Survey footprint area [deg^2].  Assumed approximately square, as in
        ``FOV_to_cMpc.transverse_comoving_size_from_area``.
    z_central : float
        Central redshift of the analysis.
    delta_z : float
        Full redshift depth of the box.  Set from the photometric redshift
        uncertainty — see Notes.
    cosmo : astropy.cosmology.FLRW, optional
        Cosmology used for the distance calculations.  Defaults to
        ``Planck18``, matching the ``cosmo=None`` convention of
        ``volume_from_area`` and ``survey_area_from_volume`` above, and the
        Planck18 comoving distances ``run_simulation.py`` §1 already uses for
        the lightcone endpoints.
    target_cell_size_mpc : float, optional
        Comoving cell size the grid should preserve [Mpc].  ``hii_dim`` is
        derived from this so that growing the box does not silently coarsen
        the mass resolution.  Default 2.0 Mpc, the resolution of the previous
        256 Mpc / 128^3 production grid (M_cell = 3.17e11 M_sun).
    hii_dim : int, optional
        Pin the grid size explicitly instead of deriving it.  Passing this
        *does* change the cell size and hence the mass resolution.
    snap_hii_dim_to_power_of_two : bool, optional
        Round a derived ``hii_dim`` up to the next power of two, which keeps
        the FFTs in ``src/analysis.py`` efficient.  Ignored when ``hii_dim``
        is given.  Default True.

    Returns
    -------
    SimulationBox
        Box geometry; ``.simulation_options`` gives the 21cmFAST kwargs.

    Raises
    ------
    ValueError
        If ``area_deg2`` or ``delta_z`` is not positive, or if ``delta_z`` is
        wide enough to reach z <= 0.

    Notes
    -----
    **Transverse extent.**  Small-angle approximation, treating the footprint
    as square:

        Omega [sr] = area_deg2 * (pi / 180)^2
        L_perp     = sqrt(Omega) * D_M(z_central)

    with ``D_M`` the comoving transverse distance.  This is the same
    construction as
    ``src/FOV_to_cMpc.py:transverse_comoving_size_from_area``.

    **Line-of-sight extent.**  By differencing the comoving distance at the
    two edges,

        L_los = D_C(z_central + delta_z/2) - D_C(z_central - delta_z/2)

    rather than ``dD_C/dz * delta_z``, matching how the codebase already
    computes the lightcone LOS extent (``run_simulation.py`` §1 and the
    notebook's derived-quantities cell both use
    ``L_los = D_C(z_max) - D_C(z_min)``).

    ``L_los`` is returned for reference but is **not** a 21cmFAST box
    argument.  21cmFAST boxes are cubic and a lightcone's LOS extent is set
    by the redshift range passed to ``RectilinearLightconer``; use ``z_min``
    and ``z_max`` for that.  When ``n_los_tiles > 1`` the coeval box is
    tiled along the line of sight, repeating structure on the ``box_len``
    scale.

    **Assumptions for the Euclid x HERA forecast** (the thesis defaults):

    - *Survey geometry.*  Euclid Deep Field Fornax, 10 deg^2, centred on
      RA 03:31:43.6, Dec -28:05:18.6.  Treated as a square footprint; the
      real field is not, so ``L_perp`` is an equivalent-square side.
    - *Central redshift.*  z = 7.
    - *Redshift depth.*  Set by the photometric redshift uncertainty,
      sigma_z = 0.256 **absolute** at z = 7.  This comes from Euclid
      Collaboration: Allen et al. (2026), A&A 711, A25, Sect. 3 / Fig. 4,
      ``sigma_nmad <= 0.032``, converted on the standard NMAD normalisation
      as ``0.032 * (1 + z)``; it is the same absolute sigma_z that
      ``src/analysis.py:radial_smearing_length`` consumes.  Two caveats on
      that value -- a field mismatch and a normalisation ambiguity -- are
      recorded in ``NUMBERS_AND_SOURCES.md`` section 5.  The multiple of
      sigma_z is a **deliberate choice, not a default of this function** —
      the caller passes ``delta_z`` explicitly.  The forecast adopts
      +/-1 sigma, ``delta_z = 2 * 0.256 = 0.512`` (z = 6.744 - 7.256,
      L_los = 179.3 Mpc); +/-2 sigma would give ``delta_z = 1.024``
      (z = 6.488 - 7.512, L_los = 359.3 Mpc).
    - *Cosmology.*  Planck18 via astropy, consistent with the comoving
      distances used elsewhere for lightcone geometry.  Note the analysis
      functions in ``src/analysis.py`` take literal H_0 = 67.36,
      Omega_m = 0.315 instead; the two differ by ~0.4% in H_0 and are not
      interchangeable at that precision.

    References
    ----------
    Hogg (1999), "Distance Measures in Cosmology" — D_M, D_C definitions.
    Euclid Collaboration, Scaramella et al. (2022), A&A 662, A112 — Euclid
    Deep Field Fornax footprint.
    Euclid Collaboration, Blanchard et al. (2020), A&A 642, A191 — photo-z
    requirement sigma_z/(1+z) < 0.05.

    Examples
    --------
    >>> box = survey_area_to_box_size(10.0, 7.0, 0.90)
    >>> f"{box.box_len:.1f}"
    '486.3'
    >>> box.hii_dim, box.dim
    (256, 768)
    """
    if area_deg2 <= 0:
        raise ValueError(f"area_deg2 must be positive, got {area_deg2}.")
    if delta_z <= 0:
        raise ValueError(f"delta_z must be positive, got {delta_z}.")
    if target_cell_size_mpc <= 0:
        raise ValueError(
            f"target_cell_size_mpc must be positive, got {target_cell_size_mpc}."
        )

    cosmo = Planck18 if cosmo is None else cosmo

    z_min = z_central - 0.5 * delta_z
    z_max = z_central + 0.5 * delta_z
    if z_min <= 0:
        raise ValueError(
            f"delta_z = {delta_z} about z_central = {z_central} reaches "
            f"z = {z_min} <= 0."
        )

    # ── Transverse extent: small-angle approximation on a square footprint ──
    solid_angle = area_deg2_to_steradians(area_deg2)
    transverse_distance = cosmo.comoving_transverse_distance(z_central)

    # Steradians are dimensionless here, so Omega * D_M^2 is an area.
    transverse_area = solid_angle.value * transverse_distance**2
    box_len = float(np.sqrt(transverse_area).to(u.Mpc).value)

    # ── Line-of-sight extent: difference of comoving distances at the edges ─
    los_depth = float(
        (
            cosmo.comoving_distance(z_max) - cosmo.comoving_distance(z_min)
        ).to(u.Mpc).value
    )

    # ── Grid: derived from the target cell size so M_cell is preserved ──────
    if hii_dim is None:
        hii_dim = int(np.ceil(box_len / target_cell_size_mpc))
        if snap_hii_dim_to_power_of_two:
            hii_dim = _next_power_of_two(hii_dim)
    elif hii_dim < 1:
        raise ValueError(f"hii_dim must be a positive integer, got {hii_dim}.")

    return SimulationBox(
        box_len=box_len,
        hii_dim=int(hii_dim),
        dim=_DIM_PER_HII_DIM * int(hii_dim),
        cell_size=box_len / hii_dim,
        los_depth=los_depth,
        z_min=float(z_min),
        z_max=float(z_max),
        transverse_area=float(transverse_area.to(u.Mpc**2).value),
        solid_angle_sr=float(solid_angle.value),
        n_los_tiles=los_depth / box_len,
    )
