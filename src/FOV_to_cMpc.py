#!/usr/bin/env python3
"""
euclid_volume_converter.py

Convert a survey field-of-view and redshift range into:

1. Solid angle
2. Comoving survey volume
3. Transverse comoving area
4. Equivalent square side length

Example
-------
python euclid_volume_converter.py --area-deg2 0.53 --z-min 6.0 --z-max 7.0
"""

import argparse
import numpy as np
import astropy.units as u
from astropy.cosmology import Planck18
from scipy.integrate import simpson


def survey_volume_from_area(area_deg2, z_min, z_max, cosmology=Planck18, n_z=1000):
    area_sr = (area_deg2 * u.deg**2).to(u.sr)

    z_grid = np.linspace(z_min, z_max, n_z)
    dV_dz_dOmega = cosmology.differential_comoving_volume(z_grid)

    volume_per_sr = simpson(dV_dz_dOmega.value, x=z_grid) * u.Mpc**3 / u.sr
    volume = area_sr * volume_per_sr

    return volume.to(u.Mpc**3), area_sr, volume_per_sr


def transverse_comoving_size_from_area(area_deg2, z, cosmology=Planck18):
    """
    Convert angular survey area into an approximate transverse comoving size.

    Assumes the field is approximately square.

    Parameters
    ----------
    area_deg2 : float
        Survey area in square degrees.

    z : float
        Redshift at which the transverse comoving size is evaluated.

    cosmology : astropy.cosmology.FLRW, optional
        Cosmology used for the distance calculation.

    Returns
    -------
    side_length_mpc : astropy.units.Quantity
        Approximate transverse comoving side length in Mpc.

    transverse_area_mpc2 : astropy.units.Quantity
        Approximate transverse comoving area in Mpc^2.
    """

    area_sr = (area_deg2 * u.deg**2).to(u.sr)

    transverse_distance = cosmology.comoving_transverse_distance(z)

    # Treat steradians as dimensionless here.
    transverse_area = area_sr.value * transverse_distance**2

    side_length = np.sqrt(transverse_area)

    return side_length.to(u.Mpc), transverse_area.to(u.Mpc**2)

def main():
    parser = argparse.ArgumentParser(
        description="Cute little Euclid survey geometry converter ✨"
    )

    parser.add_argument(
        "--area-deg2",
        type=float,
        required=True,
        help="Survey field of view / area in square degrees.",
    )

    parser.add_argument(
        "--z-min",
        type=float,
        required=True,
        help="Lower redshift bound.",
    )

    parser.add_argument(
        "--z-max",
        type=float,
        required=True,
        help="Upper redshift bound.",
    )

    parser.add_argument(
        "--n-z",
        type=int,
        default=1000,
        help="Number of redshift samples for integration. Default: 1000.",
    )

    args = parser.parse_args()

    if args.z_max <= args.z_min:
        raise ValueError("z_max must be greater than z_min.")

    z_mid = 0.5 * (args.z_min + args.z_max)

    volume, area_sr, volume_per_sr = survey_volume_from_area(
        area_deg2=args.area_deg2,
        z_min=args.z_min,
        z_max=args.z_max,
        n_z=args.n_z,
    )

    side_length, transverse_area = transverse_comoving_size_from_area(
        area_deg2=args.area_deg2,
        z=z_mid,
    )

    print("\n" + "=" * 64)
    print("✨ Euclid Survey Geometry Converter ✨")
    print("=" * 64)

    print("\nInput")
    print("-----")
    print(f"Survey area              : {args.area_deg2:.6g} deg^2")
    print(f"Redshift range           : {args.z_min:.3f} – {args.z_max:.3f}")
    print(f"Midpoint redshift        : {z_mid:.3f}")
    print(f"Cosmology                : Planck18")

    print("\nAngular Geometry")
    print("----------------")
    print(f"Solid angle              : {area_sr.value:.6e} sr")

    print("\nComoving Volume")
    print("----------------")
    print(f"Volume per steradian     : {volume_per_sr.value:.6e} Mpc^3 sr^-1")
    print(f"Survey volume            : {volume.value:.6e} Mpc^3")

    print("\nTransverse Footprint")
    print("--------------------")
    print(f"Evaluated at z           : {z_mid:.3f}")
    print(f"Transverse area          : {transverse_area.value:.6e} Mpc^2")
    print(f"Equivalent square side   : {side_length.value:.3f} Mpc")

    print("\nDone ✨\n")


if __name__ == "__main__":
    main()