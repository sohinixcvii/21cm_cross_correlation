#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/analysis.py``."""

from __future__ import annotations

import numpy as np
import pytest

from src import analysis
from src.dataio import SimulationData


# ===========================================================================
#  Cosmology
# ===========================================================================

def test_hubble_parameter_at_z0_is_h0() -> None:
    """H(0) must return H_0 exactly."""
    assert analysis.hubble_parameter(0.0, hubble_constant=67.36) == pytest.approx(67.36)


def test_hubble_parameter_increases_with_redshift() -> None:
    """H(z) is monotonically increasing for ΛCDM."""
    z = np.linspace(0, 10, 20)
    h_z = analysis.hubble_parameter(z)
    assert np.all(np.diff(h_z) > 0)


def test_comoving_distance_matches_astropy() -> None:
    """D_c(z) agrees with astropy's Planck18 to better than 1 %."""
    astropy_cosmology = pytest.importorskip("astropy.cosmology")
    expected = astropy_cosmology.Planck18.comoving_distance(7.0).value
    computed = analysis.comoving_distance(7.0)
    assert computed == pytest.approx(expected, rel=0.01)


# ===========================================================================
#  Power spectra
# ===========================================================================

def test_cross_power_shapes_and_bins() -> None:
    """Returned arrays have the requested binning and are finite where filled."""
    rng = np.random.default_rng(0)
    field = rng.standard_normal((16, 16, 12))

    k_perp, k_par, power, counts = analysis.compute_cylindrical_cross_power(
        field, field, box_len_perp=64.0, box_len_los=100.0,
        n_bins_perp=6, n_bins_parallel=5,
    )

    assert k_perp.shape == (6,)
    assert k_par.shape == (5,)
    assert power.shape == (6, 5)
    assert counts.shape == (6, 5)
    assert np.all(np.isnan(power[counts == 0]))


def test_auto_power_is_positive() -> None:
    """An auto-spectrum is non-negative in every populated bin."""
    rng = np.random.default_rng(1)
    field = rng.standard_normal((16, 16, 12))

    _, _, power, counts = analysis.compute_cylindrical_cross_power(
        field, field, 64.0, 100.0, 6, 6,
    )
    assert np.all(power[counts > 0] >= 0)


def test_cross_power_is_symmetric_in_its_fields() -> None:
    """P_ab equals P_ba for real fields."""
    rng = np.random.default_rng(2)
    field_a = rng.standard_normal((16, 16, 12))
    field_b = rng.standard_normal((16, 16, 12))

    _, _, power_ab, _ = analysis.compute_cylindrical_cross_power(
        field_a, field_b, 64.0, 100.0, 6, 6,
    )
    _, _, power_ba, _ = analysis.compute_cylindrical_cross_power(
        field_b, field_a, 64.0, 100.0, 6, 6,
    )
    np.testing.assert_allclose(power_ab, power_ba, rtol=1e-10)


def test_anticorrelated_fields_give_negative_cross_power() -> None:
    """A sign-flipped field yields a strictly negative cross-spectrum."""
    rng = np.random.default_rng(3)
    field = rng.standard_normal((16, 16, 12))

    _, _, power, counts = analysis.compute_cylindrical_cross_power(
        field, -field, 64.0, 100.0, 6, 6,
    )
    assert np.all(power[counts > 0] <= 0)


def test_white_noise_power_matches_analytic_value() -> None:
    """White noise of variance σ² has P = σ² V_cell, to within sampling noise."""
    rng = np.random.default_rng(4)
    n_perp, n_los = 32, 32
    box_perp, box_los = 64.0, 64.0
    sigma = 2.0

    field = rng.standard_normal((n_perp, n_perp, n_los)) * sigma
    _, _, power, counts = analysis.compute_cylindrical_cross_power(
        field, field, box_perp, box_los, 6, 6,
    )

    cell_volume = (box_perp / n_perp) ** 2 * (box_los / n_los)
    expected = sigma ** 2 * cell_volume
    measured = np.nanmean(power[counts > 20])

    assert measured == pytest.approx(expected, rel=0.15)


def test_field_shape_mismatch_raises() -> None:
    """Mismatched field shapes are rejected rather than broadcast."""
    with pytest.raises(ValueError, match="field shapes differ"):
        analysis.compute_cylindrical_cross_power(
            np.zeros((8, 8, 4)), np.zeros((8, 8, 5)), 64.0, 100.0,
        )


def test_compute_all_power_spectra_shares_k_grid(tiny_sim: SimulationData) -> None:
    """All three spectra land on one shared k-grid with consistent shapes."""
    spectra = analysis.compute_all_power_spectra(
        tiny_sim.brightness_temp_field,
        tiny_sim.galaxy_overdensity,
        box_len_perp=tiny_sim.BOX_LEN,
        box_len_los=tiny_sim.L_los,
        n_bins_perp=8,
        n_bins_parallel=8,
    )

    assert spectra.P_21cm_auto.shape == (8, 8)
    assert spectra.P_galaxy_auto.shape == (8, 8)
    assert spectra.P_cross.shape == (8, 8)
    assert spectra.k_perp.shape == (8,)
    assert spectra.k_parallel.shape == (8,)

    populated = spectra.mode_counts > 0
    assert np.all(spectra.P_21cm_auto[populated] >= 0)
    # The synthetic galaxy field is built anti-correlated with δT_b.
    assert np.nanmean(spectra.P_cross[populated]) < 0


# ===========================================================================
#  Wedge geometry
# ===========================================================================

def test_horizon_slope_matches_definition() -> None:
    """The horizon slope equals D_c H / [c (1+z)]."""
    z_obs = 7.0
    expected = (
        analysis.comoving_distance(z_obs)
        * analysis.hubble_parameter(z_obs)
        / (3e5 * (1 + z_obs))
    )
    assert analysis.horizon_wedge_slope(z_obs) == pytest.approx(expected)


def test_fov_slope_is_below_horizon_slope() -> None:
    """The primary-beam wedge is always inside the horizon wedge."""
    horizon = analysis.horizon_wedge_slope(7.0)
    fov = analysis.fov_wedge_slope(7.0, dish_diameter=14.0)
    assert 0 < fov < horizon


def test_wedge_mask_excises_low_k_parallel() -> None:
    """Modes below the wedge line plus buffer are masked out."""
    k_perp = np.array([0.1, 1.0])
    k_parallel = np.array([0.01, 10.0])

    mask = analysis.foreground_wedge_mask(k_perp, k_parallel, slope=1.0, buffer=0.02)

    assert mask.shape == (2, 2)
    assert not mask[0, 0]   # k_par = 0.01 < 0.1 * 1.0 + 0.02
    assert mask[0, 1]       # k_par = 10 is well outside the wedge


def test_wedge_mask_buffer_is_monotonic() -> None:
    """A larger buffer can only remove modes, never add them."""
    k_perp = np.logspace(-2, 0, 10)
    k_parallel = np.logspace(-2, 0, 10)

    small = analysis.foreground_wedge_mask(k_perp, k_parallel, 1.0, buffer=0.0)
    large = analysis.foreground_wedge_mask(k_perp, k_parallel, 1.0, buffer=0.5)

    assert large.sum() <= small.sum()
    assert np.all(large <= small)


# ===========================================================================
#  Photo-z damping
# ===========================================================================

def test_radial_smearing_scales_with_sigma_z() -> None:
    """σ_r is linear in σ_z."""
    single = analysis.radial_smearing_length(0.01, 7.0)
    double = analysis.radial_smearing_length(0.02, 7.0)
    assert double == pytest.approx(2 * single)


def test_photoz_kernel_is_unity_at_k_zero_and_decays() -> None:
    """W(0) = 1 and W decreases monotonically with k_parallel."""
    k_parallel = np.array([0.0, 0.1, 0.5, 1.0])
    kernel = analysis.photoz_damping_kernel(k_parallel, radial_smearing=10.0)

    assert kernel.shape == (1, 4)
    assert kernel[0, 0] == pytest.approx(1.0)
    assert np.all(np.diff(kernel[0]) < 0)
    assert np.all(kernel > 0)


# ===========================================================================
#  Noise and SNR
# ===========================================================================

def test_thermal_noise_falls_with_integration_time() -> None:
    """P_N ∝ 1 / t_int."""
    short = analysis.hera_thermal_noise_power(7.0, 3600.0, 8e6)
    long = analysis.hera_thermal_noise_power(7.0, 7200.0, 8e6)
    assert long == pytest.approx(short / 2)


def test_snr_increases_with_signal() -> None:
    """A larger cross-power gives a larger per-mode SNR."""
    shape = (4, 4)
    p_21 = np.full(shape, 1e3)
    p_gal = np.full(shape, 1e2)

    weak = analysis.cross_power_snr(np.full(shape, 1.0), p_21, p_gal, 10.0, 300.0)
    strong = analysis.cross_power_snr(np.full(shape, 10.0), p_21, p_gal, 10.0, 300.0)

    assert np.all(strong.snr_per_mode > weak.snr_per_mode)
    assert strong.total_snr > weak.total_snr


def test_snr_masks_wedge_modes() -> None:
    """Masked modes are NaN in the map and excluded from the total."""
    shape = (3, 3)
    mask = np.zeros(shape, dtype=bool)
    mask[0, 0] = True

    result = analysis.cross_power_snr(
        np.full(shape, 5.0), np.full(shape, 1e3), np.full(shape, 1e2),
        10.0, 300.0, outside_wedge=mask,
    )

    assert np.isnan(result.snr_outside_wedge[1, 1])
    assert not np.isnan(result.snr_outside_wedge[0, 0])
    assert result.total_snr == pytest.approx(result.snr_per_mode[0, 0])


def test_total_snr_is_quadrature_sum() -> None:
    """total_snr adds the per-mode values in quadrature and ignores NaNs."""
    snr_map = np.array([[3.0, 4.0], [np.nan, 0.0]])
    assert analysis.total_snr(snr_map) == pytest.approx(5.0)


# ===========================================================================
#  Euclid selection and bias
# ===========================================================================

def test_euclid_sfr_window_is_ordered() -> None:
    """The bright cut implies a higher SFR than the faint cut."""
    sfr_min, sfr_max = analysis.euclid_sfr_window(M_UV_faint=-18.0, M_UV_bright=-22.0)
    assert 0 < sfr_min < sfr_max


def test_select_euclid_halos_applies_magnitude_window() -> None:
    """Only halos inside the magnitude window survive, and SFR>0 is required."""
    # SFR chosen to straddle the window: 0.1, 10, 1000 M_sun/yr.
    sfr = np.array([0.0, 0.1, 10.0, 1000.0, np.nan])
    masses = np.array([1e8, 1e9, 1e10, 1e12, 1e11])

    selection = analysis.select_euclid_halos(
        sfr, masses, M_UV_faint=-18.0, M_UV_bright=-22.0,
    )

    assert selection.n_valid == 3          # zero and NaN SFR excluded
    assert selection.n_selected == selection.M_UV.size
    assert np.all(selection.M_UV >= -22.0)
    assert np.all(selection.M_UV <= -18.0)
    assert selection.halo_masses.size == selection.n_selected


def test_select_euclid_halos_rejects_length_mismatch() -> None:
    """SFR and mass arrays must be the same length."""
    with pytest.raises(ValueError, match="length mismatch"):
        analysis.select_euclid_halos(np.ones(5), np.ones(4))


def test_effective_galaxy_bias_on_synthetic_catalogue(
    tiny_sim: SimulationData,
) -> None:
    """The bias estimate runs end to end and returns a plausible b_g > 1."""
    pytest.importorskip("hmf")

    selection = analysis.select_euclid_halos(
        tiny_sim.sfr, tiny_sim.halo_masses,
        M_UV_faint=tiny_sim.get("M_UV_limit"), M_UV_bright=-22.0,
    )
    if selection.n_selected == 0:
        pytest.skip("no synthetic halos in the Euclid window")

    bias = analysis.effective_galaxy_bias(selection, z_obs=tiny_sim.z_obs)

    assert bias.mean_bias > 1.0
    assert bias.bias_min <= bias.mean_bias <= bias.bias_max
    assert bias.log10_mass_grid.shape == bias.bias_grid.shape
    assert bias.n_selected == selection.n_selected


def test_effective_galaxy_bias_rejects_empty_selection() -> None:
    """An empty selection is an error, not a silent NaN."""
    pytest.importorskip("hmf")

    empty = analysis.select_euclid_halos(
        np.array([1e-30]), np.array([1e8]),
        M_UV_faint=-30.0, M_UV_bright=-32.0,
    )
    with pytest.raises(ValueError, match="no halos passed"):
        analysis.effective_galaxy_bias(empty, z_obs=7.0)
