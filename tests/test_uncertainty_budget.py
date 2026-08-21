#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for the uncertainty budget (``src/analysis.py`` §5b, ``src/dataio.py``).

The first group is a **regression lock against the source notebook**,
``21cmfast_HERAxEuclid_lightcone.ipynb``.  Its photo-z / wedge / noise / SNR
cells are transcribed here verbatim and asserted to agree with
``compute_uncertainty_budget`` to machine precision.

The notebook has since been rewritten to *import* ``compute_uncertainty_budget``
rather than carry its own copy, so the two can no longer diverge by
construction.  These tests are kept anyway: the transcription below is an
independent statement of what the formulas are, so an edit to
``src/analysis.py`` that changes the physics fails here rather than silently
propagating to both the pipeline and the notebook at once.

The constants below are therefore the notebook's **former** configuration —
sigma_z = 0.059, buffer = 0.02 — because those are the parameters its
*published* numbers were produced with, which is what makes them a usable
fixed reference point.  They are deliberately not tracked forward as the
notebook's configuration changes; see ``docs/uncertainty_budget.md`` §4.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
from scipy.integrate import quad

from src import analysis
from src.dataio import (
    PowerSpectra,
    load_uncertainty_budget,
    save_power_spectra,
    save_uncertainty_budget,
)

# ── The notebook's configuration cell as it stood when its stored outputs
#    were produced.  See the module docstring: these are a fixed historical
#    reference point, not the notebook's current values. ─────────────────────
NB_SPEED_OF_LIGHT_KMS = 3e5
NB_SPEED_OF_LIGHT_MPS = 3e8
NB_F_21_HZ = 1420.405e6
NB_DISH_DIAMETER = 14.0
NB_INTEGRATION_TIME = 1000 * 3600
NB_BANDWIDTH = 8e6
NB_MEAN_GALAXY_DENSITY = 3e-3
NB_PHOTOZ_UNCERTAINTY = 0.059
NB_WEDGE_BUFFER = 0.02
NB_Z_OBS = 7.0

# The cosmology the notebook's *stored outputs* were produced with, before the
# uncommitted switch to astropy's Planck18 in its configuration cell.
NB_HUBBLE_CONSTANT = 67.36
NB_OMEGA_M = 0.315

NB_COSMOLOGY = dict(
    hubble_constant=NB_HUBBLE_CONSTANT,
    omega_m=NB_OMEGA_M,
    speed_of_light_kms=NB_SPEED_OF_LIGHT_KMS,
)


def _nb_hubble_parameter(z: float) -> float:
    """The notebook's inline H(z), flat ΛCDM."""
    return NB_HUBBLE_CONSTANT * np.sqrt(
        NB_OMEGA_M * (1 + z) ** 3 + (1 - NB_OMEGA_M)
    )


@pytest.fixture(scope="module")
def notebook_grid() -> tuple:
    """
    The notebook's ``(k_perp, k_parallel)`` bin centres.

    Reproduces its 256 Mpc × 350.8 Mpc lightcone on a 128² × 175 grid with
    20 × 20 log bins, matching the ranges printed in the stored notebook
    output: k_perp [0.0140, 2.046], k_parallel [0.0102, 1.44] Mpc^-1.

    Returns
    -------
    tuple of ndarray
        ``(k_perp, k_parallel)``.
    """
    box_len, hii_dim = 256.0, 128
    l_los, n_z = 350.6, 175

    kx = np.fft.fftfreq(hii_dim, d=box_len / hii_dim) * 2 * np.pi
    kz = np.fft.fftfreq(n_z, d=l_los / n_z) * 2 * np.pi

    perp_edges = np.logspace(
        np.log10(0.5 * 2 * np.pi / box_len),
        np.log10(np.sqrt(2) * np.abs(kx).max() * 1.05),
        21,
    )
    par_edges = np.logspace(
        np.log10(0.5 * 2 * np.pi / l_los),
        np.log10(np.abs(kz).max() * 1.05),
        21,
    )
    return (
        np.sqrt(perp_edges[:-1] * perp_edges[1:]),
        np.sqrt(par_edges[:-1] * par_edges[1:]),
    )


@pytest.fixture(scope="module")
def notebook_spectra(notebook_grid) -> PowerSpectra:
    """
    Synthetic spectra on the notebook grid, including empty (NaN) bins.

    The values are arbitrary but fixed; what the regression tests check is
    that the *transformation chain* matches, not the input spectra.

    Returns
    -------
    PowerSpectra
        Spectra with one NaN bin, as a real run would have.
    """
    k_perp, k_parallel = notebook_grid
    shape = (k_perp.size, k_parallel.size)
    rng = np.random.default_rng(0)

    p_21 = np.abs(rng.lognormal(3.0, 1.0, shape))
    p_gal = np.abs(rng.lognormal(2.0, 1.0, shape))
    p_cross = -np.abs(rng.lognormal(2.5, 1.0, shape))
    for array in (p_21, p_gal, p_cross):
        array[3, 4] = np.nan          # an empty bin, as in a real run

    return PowerSpectra(
        k_perp=k_perp,
        k_parallel=k_parallel,
        P_21cm_auto=p_21,
        P_galaxy_auto=p_gal,
        P_cross=p_cross,
        mode_counts=np.full(shape, 100.0),
    )


@pytest.fixture(scope="module")
def notebook_budget(notebook_spectra) -> analysis.UncertaintyBudget:
    """
    The pipeline's budget, evaluated with the notebook's own parameters.

    Returns
    -------
    UncertaintyBudget
        Budget for the notebook configuration.
    """
    return analysis.compute_uncertainty_budget(
        spectra=notebook_spectra,
        z_obs=NB_Z_OBS,
        photoz_uncertainty=NB_PHOTOZ_UNCERTAINTY,
        wedge_buffer=NB_WEDGE_BUFFER,
        integration_time=NB_INTEGRATION_TIME,
        bandwidth=NB_BANDWIDTH,
        mean_galaxy_density=NB_MEAN_GALAXY_DENSITY,
        dish_diameter=NB_DISH_DIAMETER,
        f_21_hz=NB_F_21_HZ,
        speed_of_light_mps=NB_SPEED_OF_LIGHT_MPS,
        **NB_COSMOLOGY,
    )


# ===========================================================================
#  Regression against the notebook's published numbers
# ===========================================================================

def test_matches_notebook_printed_wedge_geometry() -> None:
    """
    Reproduce the geometry the notebook printed: D_c = 8821 Mpc, m = 3.151,
    theta_FoV = 6.9 deg, m_FoV = 0.379.
    """
    d_c = analysis.comoving_distance(NB_Z_OBS, **NB_COSMOLOGY)
    horizon = analysis.horizon_wedge_slope(NB_Z_OBS, **NB_COSMOLOGY)
    fov = analysis.fov_wedge_slope(
        NB_Z_OBS, dish_diameter=NB_DISH_DIAMETER, f_21_hz=NB_F_21_HZ,
        speed_of_light_mps=NB_SPEED_OF_LIGHT_MPS, **NB_COSMOLOGY,
    )

    assert d_c == pytest.approx(8821.0, abs=1.0)
    assert horizon == pytest.approx(3.151, abs=5e-4)
    assert fov == pytest.approx(0.379, abs=5e-4)


def test_matches_notebook_printed_radial_smearing() -> None:
    """The notebook printed σ_r = 20.6 Mpc for σ_z = 0.059 at z = 7."""
    smearing = analysis.radial_smearing_length(
        NB_PHOTOZ_UNCERTAINTY, NB_Z_OBS, **NB_COSMOLOGY
    )
    assert smearing == pytest.approx(20.6, abs=0.05)


def test_matches_notebook_printed_wedge_fraction(notebook_budget) -> None:
    """
    The notebook printed "Modes outside wedge : 24.2%".

    Asserted through the same ``.1%`` formatting the notebook used: the exact
    fraction is 0.2425 = 97/400, which that format renders as 24.2 %.
    """
    assert f"{notebook_budget.fraction_outside_wedge:.1%}" == "24.2%"
    assert notebook_budget.outside_wedge.sum() == 97


def test_horizon_slope_equals_the_notebooks_longhand_expression() -> None:
    """
    ``horizon_wedge_slope`` must equal the notebook's longhand form

        lambda_obs * D_c * f_21 * H(z) / (c_mps * c_kms * (1 + z)^2)

    which is the same quantity written through the observed wavelength.
    """
    d_c, _ = quad(
        lambda z_: NB_SPEED_OF_LIGHT_KMS / _nb_hubble_parameter(z_), 0, NB_Z_OBS
    )
    lambda_obs = NB_SPEED_OF_LIGHT_MPS * (1 + NB_Z_OBS) / NB_F_21_HZ
    longhand = (
        lambda_obs * d_c * NB_F_21_HZ * _nb_hubble_parameter(NB_Z_OBS)
        / (NB_SPEED_OF_LIGHT_MPS * NB_SPEED_OF_LIGHT_KMS * (1 + NB_Z_OBS) ** 2)
    )
    assert analysis.horizon_wedge_slope(NB_Z_OBS, **NB_COSMOLOGY) == pytest.approx(
        longhand, rel=1e-12
    )


def test_budget_reproduces_the_notebook_chain(
    notebook_spectra, notebook_budget
) -> None:
    """
    End-to-end lock: the notebook's damping / wedge / noise / SNR cells,
    transcribed verbatim, must agree with ``compute_uncertainty_budget``.
    """
    k_perp = notebook_spectra.k_perp
    k_parallel = notebook_spectra.k_parallel
    hz_obs = _nb_hubble_parameter(NB_Z_OBS)

    # ── Notebook: photo-z radial smearing ─────────────────────────────────
    radial_smearing = NB_SPEED_OF_LIGHT_KMS * NB_PHOTOZ_UNCERTAINTY / hz_obs
    kernel = np.exp(
        -0.5 * k_parallel[np.newaxis, :] ** 2 * radial_smearing ** 2
    )
    p_cross_observed = notebook_spectra.P_cross * kernel
    p_galaxy_observed = notebook_spectra.P_galaxy_auto * kernel ** 2

    # ── Notebook: foreground wedge mask ───────────────────────────────────
    d_c, _ = quad(
        lambda z_: NB_SPEED_OF_LIGHT_KMS / _nb_hubble_parameter(z_), 0, NB_Z_OBS
    )
    wedge_slope = d_c * hz_obs / (NB_SPEED_OF_LIGHT_KMS * (1 + NB_Z_OBS))
    perp_grid, par_grid = np.meshgrid(k_perp, k_parallel, indexing="ij")
    outside_wedge = par_grid > (perp_grid * wedge_slope + NB_WEDGE_BUFFER)

    # ── Notebook: noise and per-mode uncertainty ──────────────────────────
    observed_frequency = NB_F_21_HZ / (1 + NB_Z_OBS)
    t_sys = (100 + 60 * (300e6 / observed_frequency) ** 2.55) * 1e3
    p_noise_21cm = t_sys ** 2 * 1e3 / (NB_INTEGRATION_TIME * NB_BANDWIDTH)
    p_noise_galaxy = 1.0 / NB_MEAN_GALAXY_DENSITY

    sigma_21cm = np.abs(notebook_spectra.P_21cm_auto) + p_noise_21cm
    sigma_galaxy = np.abs(p_galaxy_observed) + p_noise_galaxy
    sigma_cross = np.sqrt(
        0.5 * (p_cross_observed ** 2 + sigma_21cm * sigma_galaxy)
    )
    snr_per_mode = np.abs(p_cross_observed) / sigma_cross
    total = np.sqrt(np.nansum(snr_per_mode[outside_wedge] ** 2))

    # ── Compare, term by term ─────────────────────────────────────────────
    assert notebook_budget.radial_smearing == pytest.approx(radial_smearing, rel=1e-12)
    assert np.array_equal(notebook_budget.photoz_kernel, kernel)
    assert np.allclose(
        notebook_budget.P_cross_observed, p_cross_observed, equal_nan=True
    )
    assert np.allclose(
        notebook_budget.P_galaxy_observed, p_galaxy_observed, equal_nan=True
    )
    assert np.array_equal(notebook_budget.outside_wedge, outside_wedge)
    assert notebook_budget.snr.P_noise_21cm == pytest.approx(p_noise_21cm, rel=1e-12)
    assert notebook_budget.snr.P_noise_galaxy == pytest.approx(p_noise_galaxy)
    assert np.allclose(
        notebook_budget.snr.sigma_cross, sigma_cross, equal_nan=True
    )
    assert np.allclose(
        notebook_budget.snr.snr_per_mode, snr_per_mode, equal_nan=True
    )
    assert notebook_budget.total_snr == pytest.approx(total, rel=1e-12)


# ===========================================================================
#  Budget structure and derived quantities
# ===========================================================================

def test_system_temperature_splits_receiver_and_sky() -> None:
    """T_sys is the receiver floor plus a synchrotron sky ∝ ν^-2.55."""
    t_sys, frequency = analysis.system_temperature(7.0, NB_F_21_HZ)

    expected_kelvin = analysis.T_RECEIVER_K + analysis.T_SKY_300MHZ_K * (
        300e6 / frequency
    ) ** analysis.SKY_SPECTRAL_INDEX

    assert frequency == pytest.approx(NB_F_21_HZ / 8.0)
    assert t_sys == pytest.approx(expected_kelvin * 1e3)
    assert t_sys > analysis.T_RECEIVER_K * 1e3     # sky dominates at 177 MHz


def test_thermal_noise_is_consistent_with_system_temperature() -> None:
    """``hera_thermal_noise_power`` is T_sys² × norm / (t Δν)."""
    t_sys, _ = analysis.system_temperature(7.0, NB_F_21_HZ)
    expected = (
        t_sys ** 2 * analysis.NOISE_NORMALISATION_MPC3
        / (NB_INTEGRATION_TIME * NB_BANDWIDTH)
    )
    noise = analysis.hera_thermal_noise_power(
        7.0, NB_INTEGRATION_TIME, NB_BANDWIDTH, NB_F_21_HZ
    )
    assert noise == pytest.approx(expected, rel=1e-12)


def test_variance_terms_sum_to_sigma_cross_squared() -> None:
    """The two halves of Eq. 15 must reconstruct σ_cross² exactly."""
    shape = (4, 5)
    result = analysis.cross_power_snr(
        np.full(shape, -3.0), np.full(shape, 1e3), np.full(shape, 1e2),
        10.0, 300.0,
    )
    total_variance = result.cosmic_variance_term + result.noise_coupling_term
    assert np.allclose(result.sigma_cross ** 2, total_variance)


def test_cosmic_variance_fraction_is_one_without_noise(notebook_spectra) -> None:
    """
    With no noise and a cross-power far above the auto-spectra, sample
    variance carries essentially the whole budget.
    """
    shape = notebook_spectra.P_cross.shape
    result = analysis.cross_power_snr(
        P_cross_observed=np.full(shape, 1e6),
        P_21cm_auto=np.zeros(shape),
        P_galaxy_observed=np.zeros(shape),
        P_noise_21cm=0.0,
        P_noise_galaxy=0.0,
    )
    assert np.all(result.noise_coupling_term == 0.0)
    assert np.allclose(result.snr_per_mode, np.sqrt(2.0))


def test_larger_sigma_z_lowers_the_total_snr(notebook_spectra) -> None:
    """Photo-z smearing can only remove signal, never add it."""
    def budget_for(sigma_z: float) -> analysis.UncertaintyBudget:
        return analysis.compute_uncertainty_budget(
            spectra=notebook_spectra, z_obs=NB_Z_OBS,
            photoz_uncertainty=sigma_z, wedge_buffer=NB_WEDGE_BUFFER,
            **NB_COSMOLOGY,
        )

    assert budget_for(0.5).total_snr < budget_for(0.05).total_snr


def test_larger_wedge_buffer_removes_modes(notebook_spectra) -> None:
    """Raising the buffer can only shrink the usable region."""
    def budget_for(buffer: float) -> analysis.UncertaintyBudget:
        return analysis.compute_uncertainty_budget(
            spectra=notebook_spectra, z_obs=NB_Z_OBS,
            photoz_uncertainty=NB_PHOTOZ_UNCERTAINTY, wedge_buffer=buffer,
            **NB_COSMOLOGY,
        )

    strict = budget_for(0.2)
    loose = budget_for(0.0)
    assert strict.outside_wedge.sum() < loose.outside_wedge.sum()
    assert np.all(loose.outside_wedge[strict.outside_wedge])


def test_as_dict_is_json_serialisable(notebook_budget) -> None:
    """Every summary value must survive ``json.dumps``."""
    import json

    payload = notebook_budget.as_dict()
    json.dumps(payload)          # raises on numpy scalars

    assert payload["modes_total"] == notebook_budget.outside_wedge.size
    assert payload["photoz_uncertainty_sigma_z"] == NB_PHOTOZ_UNCERTAINTY
    assert payload["detection_above_5sigma"] is (notebook_budget.total_snr > 5.0)


# ===========================================================================
#  Persistence
# ===========================================================================

def test_budget_round_trips_through_hdf5(
    tmp_path, notebook_spectra, notebook_budget
) -> None:
    """Saving then loading the budget preserves every map and scalar."""
    path = str(tmp_path / "analysis_products.h5")
    save_power_spectra(path, notebook_spectra, source_path=path)
    save_uncertainty_budget(path, notebook_budget)

    maps, attrs = load_uncertainty_budget(path)

    assert maps["outside_wedge"].dtype == np.bool_
    assert np.array_equal(maps["outside_wedge"], notebook_budget.outside_wedge)
    assert np.allclose(
        maps["sigma_cross"], notebook_budget.snr.sigma_cross, equal_nan=True
    )
    assert np.allclose(
        maps["P_cross_observed"], notebook_budget.P_cross_observed, equal_nan=True
    )
    assert attrs["total_snr_sigma"] == pytest.approx(notebook_budget.total_snr)
    assert attrs["modes_outside_wedge"] == notebook_budget.outside_wedge.sum()


def test_saving_the_budget_preserves_the_cached_spectra(
    tmp_path, notebook_spectra, notebook_budget
) -> None:
    """The budget is appended, not written over the spectra."""
    from src.dataio import load_power_spectra

    path = str(tmp_path / "analysis_products.h5")
    save_power_spectra(path, notebook_spectra, source_path=path)
    save_uncertainty_budget(path, notebook_budget)

    reloaded, _ = load_power_spectra(path)
    assert np.allclose(
        reloaded.P_cross, notebook_spectra.P_cross, equal_nan=True
    )


def test_resaving_the_budget_replaces_the_group(
    tmp_path, notebook_spectra, notebook_budget
) -> None:
    """A second save overwrites rather than failing on an existing group."""
    path = str(tmp_path / "analysis_products.h5")
    save_power_spectra(path, notebook_spectra, source_path=path)
    save_uncertainty_budget(path, notebook_budget)

    wider = analysis.compute_uncertainty_budget(
        spectra=notebook_spectra, z_obs=NB_Z_OBS,
        photoz_uncertainty=0.45, wedge_buffer=NB_WEDGE_BUFFER, **NB_COSMOLOGY,
    )
    save_uncertainty_budget(path, wider)

    _, attrs = load_uncertainty_budget(path)
    assert attrs["photoz_uncertainty_sigma_z"] == pytest.approx(0.45)


def test_loading_a_missing_budget_raises(tmp_path, notebook_spectra) -> None:
    """A products file without the group gives a clear error."""
    path = str(tmp_path / "analysis_products.h5")
    save_power_spectra(path, notebook_spectra, source_path=path)

    with pytest.raises(KeyError, match="uncertainty_budget"):
        load_uncertainty_budget(path)

    assert os.path.exists(path)
