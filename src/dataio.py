#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dataio.py — HDF5 I/O for the 21 cm × galaxy cross-correlation pipeline
=======================================================================

Loading helpers for ``outputs/lightcone_data.h5`` (written by
``run_simulation.py``) and save/load helpers for the cached analysis
products written by ``run_pipeline.py``.

The halo catalogue produced by a 21cmFASTv4 run can hold >10^8 entries
(several GB in memory), so :func:`load_simulation` supports strided
subsampling via ``max_halos``.  Any quantity that depends on the *number*
density of halos (e.g. the UV luminosity function) must be rescaled by
:attr:`SimulationData.halo_sampling_factor`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import h5py
import numpy as np

__all__ = [
    "SimulationData",
    "PowerSpectra",
    "load_simulation",
    "save_power_spectra",
    "load_power_spectra",
    "products_are_stale",
]


# ===========================================================================
#  Containers
# ===========================================================================

@dataclass
class SimulationData:
    """
    Container for the contents of ``lightcone_data.h5``.

    Attributes
    ----------
    brightness_temp_field, density_field, neutral_fraction, galaxy_overdensity
        Lightcone fields of shape ``(HII_DIM, HII_DIM, N_z)``.
    lc_redshifts, lc_dist_Mpc : ndarray
        Per-slice redshift and comoving distance, shape ``(N_z,)``.
    halo_masses, halo_coords, stellar_masses, sfr : ndarray
        Per-halo catalogue at ``z_obs``.  Empty arrays when the simulation
        ran without 21cmFAST.
    halo_sampling_factor : float
        Number of catalogue halos represented by each loaded halo.  1.0 when
        the full catalogue was loaded, ``stride`` when subsampled.
    n_halos_total : int
        Number of halos in the file, before any subsampling.
    attrs : dict
        Raw HDF5 root attributes (grid, cosmology, survey, instrument).
    path : str
        Path the data was loaded from.
    """

    brightness_temp_field: np.ndarray
    density_field: np.ndarray
    neutral_fraction: np.ndarray
    galaxy_overdensity: np.ndarray
    lc_redshifts: np.ndarray
    lc_dist_Mpc: np.ndarray
    halo_masses: np.ndarray
    halo_coords: np.ndarray
    stellar_masses: np.ndarray
    sfr: np.ndarray
    attrs: Dict[str, Any] = field(default_factory=dict)
    halo_sampling_factor: float = 1.0
    n_halos_total: int = 0
    path: str = ""

    # ── Convenience accessors for the scalar metadata ──────────────────────
    @property
    def HII_DIM(self) -> int:
        """Transverse grid cells per side."""
        return int(self.attrs["HII_DIM"])

    @property
    def BOX_LEN(self) -> float:
        """Transverse comoving box length [Mpc]."""
        return float(self.attrs["BOX_LEN"])

    @property
    def N_z(self) -> int:
        """Number of line-of-sight slices."""
        return int(self.attrs["N_z"])

    @property
    def L_los(self) -> float:
        """Comoving line-of-sight extent [Mpc]."""
        return float(self.attrs["L_los"])

    @property
    def cell_size(self) -> float:
        """Transverse cell size [Mpc]."""
        return float(self.attrs["cell_size"])

    @property
    def z_min(self) -> float:
        """Low-redshift end of the lightcone."""
        return float(self.attrs["z_min"])

    @property
    def z_max(self) -> float:
        """High-redshift end of the lightcone."""
        return float(self.attrs["z_max"])

    @property
    def z_obs(self) -> float:
        """Reference (midpoint) redshift."""
        return float(self.attrs["z_obs"])

    @property
    def has_halo_catalog(self) -> bool:
        """True when a non-empty halo catalogue was loaded."""
        return self.halo_masses.size > 0

    def get(self, name: str, default: Optional[float] = None) -> float:
        """
        Read a scalar root attribute as ``float``.

        Parameters
        ----------
        name : str
            Attribute name, e.g. ``"photoz_uncertainty"``.
        default : float, optional
            Value returned when the attribute is absent.

        Returns
        -------
        float
            The attribute value, or ``default``.
        """
        if name not in self.attrs:
            if default is None:
                raise KeyError(f"attribute {name!r} missing from {self.path}")
            return float(default)
        return float(self.attrs[name])


@dataclass
class PowerSpectra:
    """
    2D cylindrical power spectra on a shared ``(k_perp, k_parallel)`` grid.

    Attributes
    ----------
    k_perp, k_parallel : ndarray
        Log-spaced bin centres [Mpc^-1], shapes ``(n_perp,)`` / ``(n_par,)``.
    P_21cm_auto : ndarray
        21 cm auto-power [mK^2 Mpc^3], shape ``(n_perp, n_par)``.
    P_galaxy_auto : ndarray
        Galaxy auto-power [Mpc^3].
    P_cross : ndarray
        21 cm × galaxy cross-power [mK Mpc^3]; negative on large scales.
    mode_counts : ndarray
        Number of Fourier modes averaged into each bin.
    """

    k_perp: np.ndarray
    k_parallel: np.ndarray
    P_21cm_auto: np.ndarray
    P_galaxy_auto: np.ndarray
    P_cross: np.ndarray
    mode_counts: np.ndarray


# ===========================================================================
#  Simulation output
# ===========================================================================

def load_simulation(
    path: str,
    max_halos: int = 0,
    load_fields: bool = True,
    load_halos: bool = True,
) -> SimulationData:
    """
    Load the simulation HDF5 written by ``run_simulation.py``.

    Parameters
    ----------
    path : str
        Path to ``lightcone_data.h5``.
    max_halos : int, optional
        Cap on the number of halos loaded.  ``0`` (default) loads the full
        catalogue.  When the catalogue is larger, a uniform stride is applied
        and :attr:`SimulationData.halo_sampling_factor` is set to the stride
        so number densities can be rescaled.
    load_fields : bool, optional
        When False, the ``(HII_DIM, HII_DIM, N_z)`` fields are returned as
        empty arrays.  Useful for catalogue-only work.
    load_halos : bool, optional
        When False, the halo catalogue is returned as empty arrays.  The
        catalogue can be several GB, so skip it for field-only work.

    Returns
    -------
    SimulationData
        Fields, catalogue, and metadata.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"simulation output not found: {path}\n"
            "Run the simulation stage first (python run_pipeline.py --sim force)."
        )

    empty = np.empty((0, 0, 0), dtype=np.float32)

    with h5py.File(path, "r") as f:
        if load_fields:
            brightness_temp_field = f["brightness_temp_field"][:]
            density_field = f["density_field"][:]
            neutral_fraction = f["neutral_fraction"][:]
            galaxy_overdensity = f["galaxy_overdensity"][:]
        else:
            brightness_temp_field = empty
            density_field = empty
            neutral_fraction = empty
            galaxy_overdensity = empty

        lc_redshifts = f["lc_redshifts"][:]
        lc_dist_Mpc = f["lc_dist_Mpc"][:]

        n_halos_total = int(f["halo_catalog/halo_masses"].shape[0])
        stride = 1

        if load_halos:
            if max_halos > 0 and n_halos_total > max_halos:
                stride = int(np.ceil(n_halos_total / max_halos))

            halo_masses = f["halo_catalog/halo_masses"][::stride]
            halo_coords = f["halo_catalog/halo_coords"][::stride]
            stellar_masses = f["halo_catalog/stellar_masses"][::stride]
            sfr = f["halo_catalog/sfr"][::stride]
        else:
            halo_masses = np.empty(0, dtype=np.float32)
            halo_coords = np.empty((0, 3), dtype=np.float32)
            stellar_masses = np.empty(0, dtype=np.float32)
            sfr = np.empty(0, dtype=np.float32)

        attrs = dict(f.attrs)

    return SimulationData(
        brightness_temp_field=brightness_temp_field,
        density_field=density_field,
        neutral_fraction=neutral_fraction,
        galaxy_overdensity=galaxy_overdensity,
        lc_redshifts=lc_redshifts,
        lc_dist_Mpc=lc_dist_Mpc,
        halo_masses=halo_masses,
        halo_coords=halo_coords,
        stellar_masses=stellar_masses,
        sfr=sfr,
        attrs=attrs,
        halo_sampling_factor=float(stride),
        n_halos_total=n_halos_total,
        path=path,
    )


# ===========================================================================
#  Cached analysis products
# ===========================================================================

def save_power_spectra(
    path: str,
    spectra: PowerSpectra,
    source_path: str,
) -> None:
    """
    Write computed power spectra to an HDF5 cache.

    Parameters
    ----------
    path : str
        Destination file, e.g. ``outputs/analysis_products.h5``.
    spectra : PowerSpectra
        Spectra to store.
    source_path : str
        Path of the simulation file the spectra were computed from.  Its
        modification time is stored so :func:`products_are_stale` can detect
        a re-run simulation.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    with h5py.File(path, "w") as f:
        f.create_dataset("k_perp", data=spectra.k_perp)
        f.create_dataset("k_parallel", data=spectra.k_parallel)
        f.create_dataset("P_21cm_auto", data=spectra.P_21cm_auto)
        f.create_dataset("P_galaxy_auto", data=spectra.P_galaxy_auto)
        f.create_dataset("P_cross", data=spectra.P_cross)
        f.create_dataset("mode_counts", data=spectra.mode_counts)

        f.attrs["source_path"] = os.path.abspath(source_path)
        f.attrs["source_mtime"] = (
            os.path.getmtime(source_path) if os.path.exists(source_path) else 0.0
        )


def load_power_spectra(path: str) -> Tuple[PowerSpectra, Dict[str, Any]]:
    """
    Read power spectra back from the HDF5 cache.

    Parameters
    ----------
    path : str
        Cache file written by :func:`save_power_spectra`.

    Returns
    -------
    PowerSpectra
        The stored spectra.
    dict
        The cache root attributes (``source_path``, ``source_mtime``).
    """
    with h5py.File(path, "r") as f:
        spectra = PowerSpectra(
            k_perp=f["k_perp"][:],
            k_parallel=f["k_parallel"][:],
            P_21cm_auto=f["P_21cm_auto"][:],
            P_galaxy_auto=f["P_galaxy_auto"][:],
            P_cross=f["P_cross"][:],
            mode_counts=f["mode_counts"][:],
        )
        attrs = dict(f.attrs)

    return spectra, attrs


def products_are_stale(products_path: str, source_path: str) -> bool:
    """
    Decide whether cached analysis products need recomputing.

    Parameters
    ----------
    products_path : str
        Cache file written by :func:`save_power_spectra`.
    source_path : str
        Current simulation HDF5.

    Returns
    -------
    bool
        True when the cache is missing, unreadable, or was computed from an
        older version of ``source_path``.
    """
    if not os.path.exists(products_path):
        return True
    if not os.path.exists(source_path):
        return False

    try:
        with h5py.File(products_path, "r") as f:
            cached_mtime = float(f.attrs.get("source_mtime", 0.0))
    except OSError:
        return True

    # Allow a small tolerance: filesystem mtimes can lose sub-second precision.
    return os.path.getmtime(source_path) > cached_mtime + 1.0
