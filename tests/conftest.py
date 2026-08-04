#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
conftest.py — shared pytest fixtures
=====================================

Builds a small synthetic ``lightcone_data.h5`` with the same schema as the
one written by ``run_simulation.py``, so the analysis, figure, and pipeline
tests run in well under a second without needing 21cmFAST.
"""

from __future__ import annotations

import os
import sys
from typing import Tuple

import h5py
import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from src.dataio import SimulationData, load_simulation  # noqa: E402

# Small enough to be instant, large enough that binning and percentiles work.
TINY_HII_DIM = 16
TINY_N_Z = 12
TINY_N_HALOS = 4_000


def make_tiny_fields(
    hii_dim: int = TINY_HII_DIM,
    n_z: int = TINY_N_Z,
    seed: int = 7,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate correlated synthetic lightcone fields.

    The galaxy field is built anti-correlated with the 21 cm field, as in the
    physical case, so the cross-spectrum has a well-defined sign.

    Parameters
    ----------
    hii_dim : int, optional
        Transverse cells per side.
    n_z : int, optional
        Line-of-sight slices.
    seed : int, optional
        RNG seed.

    Returns
    -------
    brightness_temp, density, neutral_fraction, galaxy_overdensity : ndarray
        Fields of shape ``(hii_dim, hii_dim, n_z)``.
    """
    rng = np.random.default_rng(seed)

    density = rng.standard_normal((hii_dim, hii_dim, n_z)) * 0.3
    neutral_fraction = (density < np.median(density)).astype(np.float64)
    brightness_temp = 27.0 * neutral_fraction * (1.0 + density)
    galaxy_overdensity = (
        -1.5 * (brightness_temp / max(brightness_temp.std(), 1e-12))
        + 0.1 * rng.standard_normal((hii_dim, hii_dim, n_z))
    )

    return brightness_temp, density, neutral_fraction, galaxy_overdensity


def write_tiny_simulation(
    path: str,
    n_halos: int = TINY_N_HALOS,
    hii_dim: int = TINY_HII_DIM,
    n_z: int = TINY_N_Z,
    seed: int = 7,
) -> str:
    """
    Write a synthetic simulation HDF5 matching the ``run_simulation.py`` schema.

    Parameters
    ----------
    path : str
        Destination file.
    n_halos : int, optional
        Number of synthetic halos.
    hii_dim : int, optional
        Transverse cells per side.
    n_z : int, optional
        Line-of-sight slices.
    seed : int, optional
        RNG seed.

    Returns
    -------
    str
        The path written.
    """
    rng = np.random.default_rng(seed)
    brightness_temp, density, neutral_fraction, galaxy = make_tiny_fields(
        hii_dim, n_z, seed
    )

    box_len = 64.0
    z_min, z_max = 6.9, 7.1
    lc_redshifts = np.linspace(z_min, z_max, n_z)
    lc_dist = np.linspace(6300.0, 6400.0, n_z)

    # Halo catalogue: power-law masses with a 21cmFAST-like SFR relation.
    # The 1e8–3e11 M_sun range puts a usable fraction of halos inside the
    # Euclid window (0.8 < SFR < 32 M_sun/yr) so the selection and bias
    # paths are genuinely exercised.
    log_masses = 8.0 + 3.5 * rng.power(0.6, n_halos)
    halo_masses = 10.0 ** log_masses
    stellar_masses = 0.05 * (halo_masses / 1e10) ** 0.5 * 0.156 * halo_masses
    sfr = stellar_masses / 5.7e8          # M_sun yr^-1, t_sf ~ 0.57 Gyr
    halo_coords = rng.uniform(0, hii_dim, size=(n_halos, 3))

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    with h5py.File(path, "w") as f:
        f.create_dataset("brightness_temp_field", data=brightness_temp.astype(np.float32))
        f.create_dataset("density_field", data=density.astype(np.float32))
        f.create_dataset("neutral_fraction", data=neutral_fraction.astype(np.float32))
        f.create_dataset("galaxy_overdensity", data=galaxy)

        f.create_dataset("lc_redshifts", data=lc_redshifts)
        f.create_dataset("lc_dist_Mpc", data=lc_dist)

        hc = f.create_group("halo_catalog")
        hc.create_dataset("halo_masses", data=halo_masses.astype(np.float32))
        hc.create_dataset("halo_coords", data=halo_coords.astype(np.float32))
        hc.create_dataset("stellar_masses", data=stellar_masses.astype(np.float32))
        hc.create_dataset("sfr", data=sfr.astype(np.float32))

        f.attrs["HII_DIM"] = hii_dim
        f.attrs["BOX_LEN"] = box_len
        f.attrs["N_z"] = n_z
        f.attrs["L_los"] = 100.0
        f.attrs["cell_size"] = box_len / hii_dim
        f.attrs["z_min"] = z_min
        f.attrs["z_max"] = z_max
        f.attrs["z_obs"] = 0.5 * (z_min + z_max)
        f.attrs["galaxy_bias"] = 8.0
        f.attrs["beta_rsd"] = 0.12
        f.attrs["mean_galaxy_density"] = 3e-3
        f.attrs["photoz_uncertainty"] = 0.059
        f.attrs["M_UV_limit"] = -18
        f.attrs["OMEGA_M_0"] = 0.315
        f.attrs["HUBBLE_CONSTANT"] = 67.36
        f.attrs["SPEED_OF_LIGHT_KMS"] = 3e5
        f.attrs["SPEED_OF_LIGHT_MPS"] = 3e8
        f.attrs["F_21_MHZ"] = 1420.405
        f.attrs["F_21_HZ"] = 1420.405e6
        f.attrs["HERA_DISH_DIAMETER"] = 14.0
        f.attrs["integration_time"] = 1000 * 3600
        f.attrs["bandwidth"] = 8e6
        f.attrs["wedge_buffer"] = 0.02
        f.attrs["n_bins_perp"] = 8
        f.attrs["n_bins_parallel"] = 8

    return path


@pytest.fixture(scope="session")
def tiny_sim_path(tmp_path_factory: pytest.TempPathFactory) -> str:
    """
    Path to a synthetic ``lightcone_data.h5`` shared by the whole test session.

    Returns
    -------
    str
        Path to the written HDF5 file.
    """
    directory = tmp_path_factory.mktemp("sim")
    return write_tiny_simulation(str(directory / "lightcone_data.h5"))


@pytest.fixture(scope="session")
def tiny_sim(tiny_sim_path: str) -> SimulationData:
    """
    Loaded :class:`src.dataio.SimulationData` for the synthetic simulation.

    Returns
    -------
    SimulationData
        The loaded synthetic simulation.
    """
    return load_simulation(tiny_sim_path)
