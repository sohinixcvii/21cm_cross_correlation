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


def test_star_formation_timescale_matches_astropy() -> None:
    """t_sf = t_STAR / H(z) agrees with astropy's Hubble time to <1 %."""
    astropy_cosmology = pytest.importorskip("astropy.cosmology")
    import astropy.units as u

    expected = 0.5 * (1.0 / astropy_cosmology.Planck18.H(7.0)).to(u.yr).value
    assert analysis.star_formation_timescale(7.0) == pytest.approx(expected, rel=0.01)


def test_star_formation_timescale_is_570_myr_at_z7() -> None:
    """The 21cmFAST 'simple' template gives ~570 Myr at z = 7, not 100 Myr.

    Regression guard: a hardcoded 100 Myr timescale was the documented cause
    of the inconsistent Part 1 galaxy bias.
    """
    t_sf = analysis.star_formation_timescale(7.0)
    assert 5.5e8 < t_sf < 5.9e8
    assert t_sf / 1e8 == pytest.approx(5.70, rel=0.02)


def test_star_formation_timescale_scales_with_t_star() -> None:
    """t_sf is linear in t_STAR and decreases with redshift."""
    assert analysis.star_formation_timescale(7.0, t_star=1.0) == pytest.approx(
        2 * analysis.star_formation_timescale(7.0, t_star=0.5)
    )
    assert analysis.star_formation_timescale(7.5) < analysis.star_formation_timescale(6.5)


def test_stellar_mass_to_sfr_inverts_the_timescale() -> None:
    """SFR = M_star / t_sf, elementwise."""
    stellar = np.array([1e8, 1e9, 1e10])
    sfr = analysis.stellar_mass_to_sfr(stellar, z=7.0)
    np.testing.assert_allclose(
        sfr, stellar / analysis.star_formation_timescale(7.0), rtol=1e-12
    )


def test_sheth_tormen_bias_expects_squared_peak_height() -> None:
    """b(nu_sq) must not re-square its argument.

    Regression guard: run_simulation.py once passed hmf's already-squared
    ``mf.nu`` into a helper that squared it again, inflating the stored
    galaxy bias from ~4.2 to 33.4.
    """
    from src.conversions import sheth_tormen_bias

    # For nu_sq = 1/a the linear term vanishes, giving an analytic value.
    a, p, delta_c = 0.707, 0.3, 1.686
    nu_sq = 1.0 / a
    expected = 1.0 + 0.0 + 2.0 * p / (delta_c * (1.0 + 1.0**p))
    assert sheth_tormen_bias(nu_sq) == pytest.approx(expected)

    # A doubly-squared argument would land far away from this.
    assert sheth_tormen_bias(nu_sq**2) != pytest.approx(expected)


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


# ===========================================================================
#  Lightcone estimator — TODO.md P0
# ===========================================================================

class TestSubtractFieldMean:
    """P0.2 — per-slice mean subtraction."""

    def test_global_mode_removes_one_scalar(self):
        field = np.arange(2 * 2 * 4, dtype=float).reshape(2, 2, 4)
        result = analysis.subtract_field_mean(field, "global")
        assert np.isclose(result.mean(), 0.0)
        assert np.isclose(np.ptp(result), np.ptp(field))

    def test_per_slice_mode_removes_a_pure_los_ramp(self):
        ramp = np.ones((4, 4, 8)) * np.arange(8)
        assert np.allclose(analysis.subtract_field_mean(ramp, "per_slice"), 0.0)
        assert not np.allclose(analysis.subtract_field_mean(ramp, "global"), 0.0)

    def test_per_slice_mode_leaves_transverse_structure(self):
        rng = np.random.default_rng(3)
        field = rng.normal(size=(8, 8, 5)) + np.arange(5) * 10.0
        result = analysis.subtract_field_mean(field, "per_slice")
        assert np.allclose(result.mean(axis=(0, 1)), 0.0, atol=1e-12)
        assert result.std() > 0.5

    def test_rejects_unknown_mode(self):
        with pytest.raises(ValueError, match="mode must be one of"):
            analysis.subtract_field_mean(np.zeros((2, 2, 2)), "nope")

    def test_rejects_non_3d_input(self):
        with pytest.raises(ValueError, match="3D lightcone"):
            analysis.subtract_field_mean(np.zeros((2, 2)), "global")


class TestBlackmanHarrisTaper:
    """P0.3 — the line-of-sight window."""

    def test_shape_and_endpoints(self):
        window = analysis.blackman_harris_taper(32)
        assert window.shape == (32,)
        assert window[0] < 1e-3 and window[-1] < 1e-3
        assert np.isclose(window.max(), 1.0, atol=0.02)

    def test_symmetric(self):
        window = analysis.blackman_harris_taper(17)
        assert np.allclose(window, window[::-1])

    def test_noise_equivalent_bandwidth(self):
        # <w^2> is the factor the estimator divides back out.  For the 4-term
        # window it tends to sum(a_i^2) with the cosine terms halved:
        # 0.35875^2 + (0.48829^2 + 0.14128^2 + 0.01168^2)/2 = 0.2580.
        expected = 0.35875 ** 2 + (
            0.48829 ** 2 + 0.14128 ** 2 + 0.01168 ** 2
        ) / 2
        assert np.isclose(np.mean(analysis.blackman_harris_taper(512) ** 2),
                          expected, atol=1e-3)

    def test_rejects_short_axis(self):
        with pytest.raises(ValueError, match="n_slices must be"):
            analysis.blackman_harris_taper(1)


class TestTaperedEstimator:
    """The taper must preserve amplitude and leave the untapered path alone."""

    def test_taper_none_is_the_historical_path(self):
        rng = np.random.default_rng(11)
        box = rng.normal(size=(8, 8, 16))
        bare = analysis.compute_cylindrical_cross_power(box, box, 32.0, 32.0)[2]
        explicit_none = analysis.compute_cylindrical_cross_power(
            box, box, 32.0, 32.0, taper=None
        )[2]
        assert np.array_equal(np.nan_to_num(bare), np.nan_to_num(explicit_none))

    def test_taper_preserves_white_noise_amplitude(self):
        rng = np.random.default_rng(12)
        box = rng.normal(size=(16, 16, 64))
        bare = analysis.compute_cylindrical_cross_power(box, box, 64.0, 64.0)[2]
        tapered = analysis.compute_cylindrical_cross_power(
            box, box, 64.0, 64.0, taper=analysis.blackman_harris_taper(64)
        )[2]
        ratio = np.nanmedian(bare / tapered)
        assert 0.7 < ratio < 1.4

    def test_rejects_taper_of_the_wrong_length(self):
        box = np.zeros((4, 4, 8))
        with pytest.raises(ValueError, match="taper length"):
            analysis.compute_cylindrical_cross_power(
                box, box, 8.0, 8.0, taper=np.ones(5)
            )


class TestSubbandIndexRanges:
    """P0.4 — splitting the line of sight to match the noise bandwidth."""

    def test_bands_tile_the_axis_without_gaps(self):
        ranges = analysis.subband_index_ranges(166, 22.28e6, 8e6)
        assert ranges[0, 0] == 0
        assert ranges[-1, 1] == 166
        assert np.all(ranges[1:, 0] == ranges[:-1, 1])

    def test_band_count_follows_the_bandwidth_ratio(self):
        assert analysis.subband_index_ranges(166, 22.28e6, 8e6).shape[0] == 3
        assert analysis.subband_index_ranges(166, 16.0e6, 8e6).shape[0] == 2

    def test_narrow_lightcone_gives_a_single_band(self):
        ranges = analysis.subband_index_ranges(100, 0.9e6, 8e6)
        assert ranges.shape == (1, 2)
        assert ranges[0].tolist() == [0, 100]

    def test_min_slices_per_band_caps_the_split(self):
        # 20 slices, 40 MHz span: 5 bands by bandwidth, but only 2 survive
        # the 8-slice floor.
        ranges = analysis.subband_index_ranges(
            20, 40e6, 8e6, min_slices_per_band=8
        )
        assert ranges.shape[0] == 2
        assert np.all(np.diff(ranges, axis=1) >= 8)

    @pytest.mark.parametrize("kwargs", [
        {"n_slices": 0},
        {"frequency_span_hz": 0.0},
        {"bandwidth_hz": -1.0},
        {"min_slices_per_band": 0},
    ])
    def test_rejects_non_physical_arguments(self, kwargs):
        base = dict(n_slices=10, frequency_span_hz=1e6, bandwidth_hz=8e6)
        base.update(kwargs)
        with pytest.raises(ValueError):
            analysis.subband_index_ranges(**base)


class TestSubbandPowerSpectra:
    """P0.3 — one spectrum per band, each with its own redshift and extent."""

    @staticmethod
    def _lightcone(n_z=48, z_min=6.5, z_max=7.5):
        rng = np.random.default_rng(5)
        t21 = rng.normal(size=(8, 8, n_z))
        gal = rng.normal(size=(8, 8, n_z))
        z = np.linspace(z_min, z_max, n_z)
        dist = 8600.0 + np.linspace(0.0, 350.0, n_z)
        return t21, gal, z, dist

    def test_band_count_and_geometry(self):
        t21, gal, z, dist = self._lightcone()
        bands, geometry = analysis.compute_subband_power_spectra(
            t21, gal, z, dist, box_len_perp=64.0, bandwidth_hz=8e6,
            min_slices_per_band=8,
        )
        assert len(bands) == geometry.n_bands == 3
        assert geometry.z_effective.shape == (3,)
        # Bands run low-z to high-z and do not overlap.
        assert np.all(np.diff(geometry.z_effective) > 0)
        assert geometry.n_slices.sum() == 48
        # Each band spans at most the requested bandwidth.
        assert np.all(geometry.bandwidth_hz <= 8e6)

    def test_each_band_has_the_requested_binning(self):
        t21, gal, z, dist = self._lightcone()
        bands, _ = analysis.compute_subband_power_spectra(
            t21, gal, z, dist, box_len_perp=64.0, bandwidth_hz=8e6,
            n_bins_perp=6, n_bins_parallel=5, min_slices_per_band=8,
        )
        for band in bands:
            assert band.P_cross.shape == (6, 5)
            assert band.k_perp.shape == (6,)
            assert band.k_parallel.shape == (5,)

    def test_effective_redshift_lies_inside_its_band(self):
        t21, gal, z, dist = self._lightcone()
        _, geometry = analysis.compute_subband_power_spectra(
            t21, gal, z, dist, box_len_perp=64.0, bandwidth_hz=8e6,
            min_slices_per_band=8,
        )
        assert np.all(geometry.z_effective >= geometry.z_min)
        assert np.all(geometry.z_effective <= geometry.z_max)

    def test_bands_sample_lower_k_parallel_than_the_whole_box(self):
        # Each band is shorter, so its fundamental k_parallel is larger.
        t21, gal, z, dist = self._lightcone()
        bands, geometry = analysis.compute_subband_power_spectra(
            t21, gal, z, dist, box_len_perp=64.0, bandwidth_hz=8e6,
            min_slices_per_band=8,
        )
        whole = analysis.compute_all_power_spectra(
            t21, gal, box_len_perp=64.0, box_len_los=350.0,
        )
        assert bands[0].k_parallel.min() > whole.k_parallel.min()
        assert np.all(geometry.los_length_mpc < 350.0)

    def test_single_band_when_the_lightcone_is_narrow(self):
        t21, gal, z, dist = self._lightcone(n_z=16, z_min=6.995, z_max=7.005)
        bands, geometry = analysis.compute_subband_power_spectra(
            t21, gal, z, dist, box_len_perp=64.0, bandwidth_hz=8e6,
            min_slices_per_band=8,
        )
        assert geometry.n_bands == 1
        assert bands[0].P_cross.shape == (20, 20)

    def test_rejects_mismatched_inputs(self):
        t21, gal, z, dist = self._lightcone(n_z=16)
        with pytest.raises(ValueError, match="must match the LOS axis"):
            analysis.compute_subband_power_spectra(
                t21, gal, z[:-1], dist, box_len_perp=64.0,
            )
        with pytest.raises(ValueError, match="field shapes differ"):
            analysis.compute_subband_power_spectra(
                t21, gal[:, :, :-1], z, dist, box_len_perp=64.0,
            )


class TestCombineBandSnr:
    """Independent bands add in quadrature."""

    def test_quadrature_sum(self):
        assert np.isclose(analysis.combine_band_snr([3.0, 4.0]), 5.0)

    def test_ignores_nan_bands(self):
        assert np.isclose(analysis.combine_band_snr([3.0, np.nan, 4.0]), 5.0)

    def test_single_band_is_itself(self):
        assert np.isclose(analysis.combine_band_snr([2.5]), 2.5)


class TestEstimatorDefaultsUnchanged:
    """The coeval formalism must survive the P0 additions untouched."""

    def test_default_mean_subtraction_is_global(self, tiny_sim):
        spectra = analysis.compute_all_power_spectra(
            tiny_sim.brightness_temp_field, tiny_sim.galaxy_overdensity,
            box_len_perp=tiny_sim.BOX_LEN, box_len_los=tiny_sim.L_los,
        )
        explicit = analysis.compute_all_power_spectra(
            tiny_sim.brightness_temp_field, tiny_sim.galaxy_overdensity,
            box_len_perp=tiny_sim.BOX_LEN, box_len_los=tiny_sim.L_los,
            mean_subtraction="global", taper=None,
        )
        assert np.array_equal(
            np.nan_to_num(spectra.P_cross), np.nan_to_num(explicit.P_cross)
        )

    def test_per_slice_zeroes_the_transverse_mean_profile(self):
        """P0.2's direct effect: no residual <T_b>(z) evolution is left."""
        rng = np.random.default_rng(21)
        field = rng.normal(size=(16, 16, 32)) + np.linspace(0.0, 40.0, 32)

        global_mode = analysis.subtract_field_mean(field, "global")
        per_slice = analysis.subtract_field_mean(field, "per_slice")

        assert np.ptp(global_mode.mean(axis=(0, 1))) > 30.0
        assert np.allclose(per_slice.mean(axis=(0, 1)), 0.0, atol=1e-12)

    def test_uniform_ramp_lives_at_k_perp_zero_and_is_already_binned_out(self):
        """
        Why P0.2 barely moves *these* spectra.

        A line-of-sight ramp that is uniform across the sky has essentially
        all its power at ``k_perp = 0``, and the log-spaced binning starts at
        ``0.5 dk_perp``, so that column never enters a bin.  Per-slice mean
        subtraction is still the correct operation — it just cannot show up
        in the binned output for this particular contaminant.  Recorded as a
        regression guard so the claim is not quietly overstated later.
        """
        rng = np.random.default_rng(21)
        signal = rng.normal(size=(16, 16, 32))
        field = signal + np.linspace(0.0, 40.0, 32)

        fourier = np.fft.fftn(field - field.mean())
        power_at_k_perp_zero = np.sum(np.abs(fourier[0, 0, :]) ** 2)
        assert power_at_k_perp_zero / np.sum(np.abs(fourier) ** 2) > 0.99

        global_mode = analysis.compute_all_power_spectra(
            field, signal, box_len_perp=32.0, box_len_los=64.0,
            mean_subtraction="global",
        )
        per_slice = analysis.compute_all_power_spectra(
            field, signal, box_len_perp=32.0, box_len_los=64.0,
            mean_subtraction="per_slice",
        )
        assert np.allclose(
            np.nan_to_num(global_mode.P_21cm_auto),
            np.nan_to_num(per_slice.P_21cm_auto),
            rtol=1e-9, atol=1e-9,
        )
