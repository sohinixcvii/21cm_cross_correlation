#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/dataio.py``."""

from __future__ import annotations

import os
import time

import numpy as np
import pytest

from src import dataio
from src.dataio import PowerSpectra, SimulationData

from conftest import TINY_HII_DIM, TINY_N_HALOS, TINY_N_Z, write_tiny_simulation


def test_load_simulation_returns_expected_shapes(tiny_sim: SimulationData) -> None:
    """Fields, geometry, and catalogue load with the documented shapes."""
    expected = (TINY_HII_DIM, TINY_HII_DIM, TINY_N_Z)

    assert tiny_sim.brightness_temp_field.shape == expected
    assert tiny_sim.density_field.shape == expected
    assert tiny_sim.neutral_fraction.shape == expected
    assert tiny_sim.galaxy_overdensity.shape == expected
    assert tiny_sim.lc_redshifts.shape == (TINY_N_Z,)
    assert tiny_sim.lc_dist_Mpc.shape == (TINY_N_Z,)
    assert tiny_sim.halo_masses.shape == (TINY_N_HALOS,)
    assert tiny_sim.halo_coords.shape == (TINY_N_HALOS, 3)
    assert tiny_sim.has_halo_catalog


def test_metadata_properties(tiny_sim: SimulationData) -> None:
    """The typed accessors return the HDF5 attributes."""
    assert tiny_sim.HII_DIM == TINY_HII_DIM
    assert tiny_sim.N_z == TINY_N_Z
    assert tiny_sim.BOX_LEN == pytest.approx(64.0)
    assert tiny_sim.z_obs == pytest.approx(7.0)
    assert tiny_sim.cell_size == pytest.approx(64.0 / TINY_HII_DIM)
    assert tiny_sim.halo_sampling_factor == 1.0


def test_get_falls_back_to_default(tiny_sim: SimulationData) -> None:
    """``get`` returns the default for a missing attribute, and raises without one."""
    assert tiny_sim.get("not_an_attribute", 42.0) == 42.0
    with pytest.raises(KeyError):
        tiny_sim.get("not_an_attribute")


def test_max_halos_subsamples_and_records_factor(tiny_sim_path: str) -> None:
    """Strided loading caps the halo count and reports the sampling factor."""
    data = dataio.load_simulation(tiny_sim_path, max_halos=500)

    assert data.halo_masses.shape[0] <= 500
    assert data.halo_sampling_factor > 1.0
    assert data.n_halos_total == TINY_N_HALOS
    assert data.halo_coords.shape[0] == data.halo_masses.shape[0]


def test_load_halos_false_skips_the_catalogue(tiny_sim_path: str) -> None:
    """Field-only loading returns empty catalogue arrays but keeps the count."""
    data = dataio.load_simulation(tiny_sim_path, load_halos=False)

    assert data.halo_masses.size == 0
    assert not data.has_halo_catalog
    assert data.n_halos_total == TINY_N_HALOS
    assert data.brightness_temp_field.size > 0


def test_load_fields_false_skips_the_fields(tiny_sim_path: str) -> None:
    """Catalogue-only loading returns empty field arrays."""
    data = dataio.load_simulation(tiny_sim_path, load_fields=False)

    assert data.brightness_temp_field.size == 0
    assert data.halo_masses.size == TINY_N_HALOS


def test_missing_file_raises() -> None:
    """A helpful error is raised when the simulation output is absent."""
    with pytest.raises(FileNotFoundError, match="simulation output not found"):
        dataio.load_simulation("/nonexistent/lightcone_data.h5")


def test_power_spectra_round_trip(tmp_path, tiny_sim_path: str) -> None:
    """Spectra survive a save/load cycle unchanged."""
    spectra = PowerSpectra(
        k_perp=np.linspace(0.1, 1.0, 5),
        k_parallel=np.linspace(0.1, 2.0, 4),
        P_21cm_auto=np.ones((5, 4)),
        P_galaxy_auto=np.full((5, 4), 2.0),
        P_cross=np.full((5, 4), -3.0),
        mode_counts=np.full((5, 4), 7.0),
    )
    path = str(tmp_path / "products.h5")

    dataio.save_power_spectra(path, spectra, tiny_sim_path)
    loaded, attrs = dataio.load_power_spectra(path)

    np.testing.assert_allclose(loaded.k_perp, spectra.k_perp)
    np.testing.assert_allclose(loaded.P_cross, spectra.P_cross)
    np.testing.assert_allclose(loaded.mode_counts, spectra.mode_counts)
    assert attrs["source_path"] == os.path.abspath(tiny_sim_path)


def test_products_are_stale_detects_missing_and_newer_source(
    tmp_path, tiny_sim_path: str,
) -> None:
    """A missing cache is stale; a cache older than the simulation is stale."""
    path = str(tmp_path / "products.h5")
    assert dataio.products_are_stale(path, tiny_sim_path)

    spectra = PowerSpectra(
        k_perp=np.ones(2), k_parallel=np.ones(2),
        P_21cm_auto=np.ones((2, 2)), P_galaxy_auto=np.ones((2, 2)),
        P_cross=np.ones((2, 2)), mode_counts=np.ones((2, 2)),
    )
    dataio.save_power_spectra(path, spectra, tiny_sim_path)
    assert not dataio.products_are_stale(path, tiny_sim_path)

    # Re-writing the simulation makes the cache stale.
    newer = str(tmp_path / "lightcone_data.h5")
    write_tiny_simulation(newer)
    os.utime(newer, (time.time() + 10, time.time() + 10))
    dataio.save_power_spectra(path, spectra, tiny_sim_path)
    assert dataio.products_are_stale(path, newer)
