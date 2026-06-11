# ===================================================================
# Useful cosmological conversions for high-redshift galaxy surveys.
# ===================================================================

# Currently includes:
# - Muv <-> Luv conversions
# - Survey volume <-> survey area conversions
# - Area conversions between square degrees and steradians
# - Volume calculation from area and redshift interval

import numpy as np
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

# Madau & Dickinson (2014) UV–SFR conversion factor
# SFR [M_sun yr^-1] = kappa_UV * L_UV [erg s^-1 Hz^-1]
# Chabrier (2003) IMF, rest-frame ~1500 Å
_KAPPA_UV_MADAU14 = 1.15e-28


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

        Default is 1.15e-28 from Madau & Dickinson (2014) for a
        Chabrier (2003) IMF.

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

        Default is 1.15e-28 from Madau & Dickinson (2014) for a
        Chabrier (2003) IMF.

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

        Default is 1.15e-28 from Madau & Dickinson (2014).

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