#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/foregrounds.py``.

Three things matter here and are checked directly: that the injected field is
spectrally smooth, and so contaminates low ``k_parallel`` far more than high;
that smoothness nonetheless does **not** confine it to low ``k_parallel`` once
it passes through the pipeline's un-tapered FFT; and that the removal knob is
exact at both ends — 0 % leaves the contaminated field alone, 100 %
reproduces the clean one bit for bit.
"""

from __future__ import annotations

import numpy as np
import pytest

from conftest import write_tiny_simulation
from src import analysis
from src.dataio import SimulationData, load_simulation
from src.foregrounds import (
    ForegroundRealisation,
    inject_foreground,
    remove_foreground,
    simulate_diffuse_foreground,
    simulate_point_source_foreground,
)

SHAPE = (16, 16, 12)
K_PERP = np.logspace(-2, 0, 10)
K_PAR = np.logspace(-2, 0, 10)
BOX = dict(box_len_perp=64.0, box_len_los=100.0)


@pytest.fixture(scope="module")
def fg_sim(tmp_path_factory) -> SimulationData:
    """
    A synthetic lightcone large enough for the end-to-end tests.

    The session-wide ``tiny_sim`` (16² x 12) leaves most ``(k_perp,
    k_parallel)`` bins empty, so the binned spectra are largely NaN and the
    total SNR collapses to exactly zero — which makes "contamination lowers
    the SNR" untestable.  32² x 24 populates the bins for a fraction of a
    second more.
    """
    path = write_tiny_simulation(
        str(tmp_path_factory.mktemp("fg") / "sim.h5"), hii_dim=32, n_z=24,
    )
    return load_simulation(path, load_halos=False)


def _low_over_high_k_parallel(field: np.ndarray) -> float:
    """
    Ratio of line-of-sight power in the lowest to the highest ``k_parallel``.

    A spectrally smooth field concentrates its power in the slowly varying
    LOS modes, so this ratio is large; white noise gives ~1.

    Parameters
    ----------
    field : ndarray
        Cube of shape ``(N, N, N_z)``.

    Returns
    -------
    float
        Mean power in the two lowest non-zero LOS modes divided by the mean
        in the two highest.
    """
    fluctuations = field - field.mean(axis=2, keepdims=True)
    power = np.abs(np.fft.fft(fluctuations, axis=2)) ** 2
    per_mode = power.mean(axis=(0, 1))
    n_half = per_mode.size // 2
    return float(per_mode[1:3].mean() / per_mode[n_half - 2:n_half].mean())


# ===========================================================================
#  Diffuse component
# ===========================================================================

def test_diffuse_shape_and_positivity() -> None:
    """The cube matches the requested grid and is a physical (positive) sky."""
    field = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)

    assert field.shape == SHAPE
    assert np.all(np.isfinite(field))
    assert field.min() > 0


def test_diffuse_is_spectrally_smooth() -> None:
    """
    Power sits preferentially at low ``k_parallel``.

    Compared against white noise on the same grid, which has no LOS
    preference at all.
    """
    field = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)
    white = np.random.default_rng(0).standard_normal(SHAPE)

    assert _low_over_high_k_parallel(field) > 3.0
    assert _low_over_high_k_parallel(white) < 2.0


def test_smoothness_does_not_confine_power_to_low_k_parallel() -> None:
    """
    A smooth foreground still leaks across the whole ``k_parallel`` axis.

    ``compute_cylindrical_cross_power`` applies no line-of-sight taper, so a
    smooth-but-non-periodic spectrum is discontinuous at the box edge and
    leaks as a mild power law rather than falling off a cliff.  The slope is
    set by that discontinuity, not by the sky: a bare ramp with no angular
    structure leaks identically, and widening the band does not change it.

    Documented as a test because it is easy to assume the opposite, and the
    assumption would make wedge excision look far more protective than it is.
    """
    def los_slope(field: np.ndarray) -> float:
        fluctuations = field - field.mean(axis=2, keepdims=True)
        power = (np.abs(np.fft.fft(fluctuations, axis=2)) ** 2).mean(axis=(0, 1))
        n_half = power.size // 2
        modes = np.arange(1, n_half + 1)
        return float(np.polyfit(np.log(modes), np.log(power[1:n_half + 1]), 1)[0])

    foreground = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)
    bare_ramp = np.ones((16, 16, 1)) * np.linspace(1.0, 1.04, 12)[None, None, :]

    assert los_slope(foreground) == pytest.approx(los_slope(bare_ramp), abs=0.05)
    assert -2.5 < los_slope(foreground) < -1.0


def test_zero_index_scatter_makes_every_sightline_identical() -> None:
    """
    With no spectral-index variation the foreground is exactly separable.

    ``T(θ, ν) = T_ref(θ) f(ν)`` means every sightline has the *same* shape in
    frequency, differing only by a constant.  That is the idealised case a
    rank-1 removal would clean perfectly; the default non-zero scatter is what
    makes real removal hard.
    """
    field = simulate_diffuse_foreground(
        SHAPE, K_PERP, K_PAR, 7.0, spectral_index_scatter=0.0, **BOX
    )
    profiles = field / field[:, :, :1]         # normalise each sightline

    assert np.allclose(profiles, profiles[0, 0], rtol=1e-10)


def test_index_scatter_leaks_power_up_the_line_of_sight() -> None:
    """Spectral-index variation is what breaks perfect smoothness."""
    smooth = simulate_diffuse_foreground(
        SHAPE, K_PERP, K_PAR, 7.0, spectral_index_scatter=0.0, **BOX
    )
    rough = simulate_diffuse_foreground(
        SHAPE, K_PERP, K_PAR, 7.0, spectral_index_scatter=0.5, **BOX
    )

    assert _low_over_high_k_parallel(rough) < _low_over_high_k_parallel(smooth)


def test_steeper_angular_index_gives_more_large_scale_power() -> None:
    """``C_l ∝ l^-β``: larger β puts relatively more power on large scales."""
    def large_scale_fraction(beta: float) -> float:
        field = simulate_diffuse_foreground(
            SHAPE, K_PERP, K_PAR, 7.0, beta_angular=beta, seed=3, **BOX
        )
        plane = field[:, :, 0]
        coarse = plane.reshape(4, 4, 4, 4).mean(axis=(1, 3))
        return float(coarse.var() / plane.var())

    assert large_scale_fraction(4.0) > large_scale_fraction(0.5)


def test_diffuse_uses_supplied_redshifts() -> None:
    """
    Passing the lightcone's own redshifts changes the frequency axis.

    The internal ``dz = H dr / c`` fallback is an approximation; a caller with
    ``lc_redshifts`` should get a different — and correct — answer.
    """
    inferred = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)
    supplied = simulate_diffuse_foreground(
        SHAPE, K_PERP, K_PAR, 7.0,
        lc_redshifts=np.linspace(6.5, 7.5, SHAPE[2]), **BOX,
    )

    assert not np.allclose(inferred, supplied)


def test_diffuse_is_reproducible() -> None:
    """Same seed, same field; different seed, different field."""
    first = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, seed=1, **BOX)
    same = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, seed=1, **BOX)
    other = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, seed=2, **BOX)

    assert np.array_equal(first, same)
    assert not np.allclose(first, other)


def test_box_lengths_are_inferred_when_absent() -> None:
    """The k-grids alone are enough to produce a field, if less exactly."""
    field = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0)
    assert field.shape == SHAPE and np.all(np.isfinite(field))


@pytest.mark.parametrize("bad_shape", [(16, 16), (16, 8, 12), (16, 16, 12, 2)])
def test_diffuse_rejects_bad_shapes(bad_shape) -> None:
    """Only a square-transverse 3D grid is meaningful."""
    with pytest.raises(ValueError):
        simulate_diffuse_foreground(bad_shape, K_PERP, K_PAR, 7.0, **BOX)


@pytest.mark.parametrize(
    "k_perp, k_par",
    [(np.array([]), K_PAR), (np.array([0.0, 0.1]), K_PAR),
     (K_PERP, np.array([-1.0, 0.1])), (K_PERP, np.array([np.nan, 0.1]))],
)
def test_diffuse_rejects_bad_k_grids(k_perp, k_par) -> None:
    """Non-positive or empty k-grids are a caller error, not something to fix up."""
    with pytest.raises(ValueError):
        simulate_diffuse_foreground(SHAPE, k_perp, k_par, 7.0, **BOX)


def test_diffuse_rejects_mismatched_redshifts() -> None:
    """``lc_redshifts`` must have one entry per LOS slice."""
    with pytest.raises(ValueError, match="lc_redshifts"):
        simulate_diffuse_foreground(
            SHAPE, K_PERP, K_PAR, 7.0, lc_redshifts=np.linspace(6.9, 7.1, 5), **BOX
        )


# ===========================================================================
#  Point-source component
# ===========================================================================

def test_point_sources_are_sparse_and_smooth() -> None:
    """Sparse on the sky, but still a smooth power law in frequency."""
    field = simulate_point_source_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)

    assert field.shape == SHAPE
    assert np.all(field >= 0)
    # Far fewer than half the cells host a source at the default density.
    assert (field[:, :, 0] > 0).mean() < 0.5
    assert _low_over_high_k_parallel(field) > 3.0


def test_point_sources_are_angularly_flatter_than_diffuse() -> None:
    """
    A Poisson population is close to white; the diffuse sky is steeply red.

    That contrast is the only reason to carry a separate component at all.
    """
    def large_scale_fraction(field: np.ndarray) -> float:
        plane = field[:, :, 0]
        coarse = plane.reshape(4, 4, 4, 4).mean(axis=(1, 3))
        return float(coarse.var() / plane.var())

    diffuse = simulate_diffuse_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)
    points = simulate_point_source_foreground(SHAPE, K_PERP, K_PAR, 7.0, **BOX)

    assert large_scale_fraction(points) < large_scale_fraction(diffuse)


def test_point_sources_reject_divergent_flux_slope() -> None:
    """``dN/dS ∝ S^-γ`` needs ``γ > 1`` for a finite mean flux."""
    with pytest.raises(ValueError, match="flux_slope"):
        simulate_point_source_foreground(
            SHAPE, K_PERP, K_PAR, 7.0, flux_slope=0.9, **BOX
        )


def test_point_sources_survive_an_empty_draw() -> None:
    """A density low enough to draw no sources must give zeros, not a crash."""
    field = simulate_point_source_foreground(
        SHAPE, K_PERP, K_PAR, 7.0, source_density_per_cell=0.0, **BOX
    )
    assert np.all(field == 0)


# ===========================================================================
#  Injection
# ===========================================================================

def test_injection_hits_the_requested_amplitude(fg_sim: SimulationData) -> None:
    """``foreground_amplitude`` is the foreground/signal RMS ratio, exactly."""
    result = inject_foreground(
        fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
        foreground_amplitude=250.0,
        box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
    )

    assert isinstance(result, ForegroundRealisation)
    assert result.amplitude == pytest.approx(250.0, rel=1e-9)
    assert result.foreground_rms == pytest.approx(250.0 * result.signal_rms, rel=1e-9)
    assert np.allclose(
        result.contaminated, fg_sim.brightness_temp_field + result.foreground
    )


def test_injection_of_zero_amplitude_is_a_no_op(fg_sim: SimulationData) -> None:
    """Zero contamination must leave the field untouched."""
    result = inject_foreground(
        fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
        foreground_amplitude=0.0,
        box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
    )
    assert np.allclose(result.contaminated, fg_sim.brightness_temp_field)


def test_injection_components_are_independent(fg_sim: SimulationData) -> None:
    """
    Diffuse and point-source components must not share an RNG stream.

    Passing one ``seed`` down to both would correlate two physically
    unrelated populations.
    """
    result = inject_foreground(
        fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
        include_point_sources=True, seed=7,
        box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
    )

    assert result.point_source is not None
    unit = lambda a: (a - a.mean()) / a.std()
    correlation = float(np.mean(unit(result.diffuse) * unit(result.point_source)))
    assert abs(correlation) < 0.5


def test_injection_can_skip_point_sources(fg_sim: SimulationData) -> None:
    """The diffuse component alone is a valid configuration."""
    result = inject_foreground(
        fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
        include_point_sources=False,
        box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
    )
    assert result.point_source is None


@pytest.mark.parametrize(
    "kwargs", [{"foreground_amplitude": -1.0}, {"point_source_weight": 1.5}]
)
def test_injection_validates_its_inputs(fg_sim: SimulationData, kwargs) -> None:
    """Out-of-range contamination settings are rejected."""
    with pytest.raises(ValueError):
        inject_foreground(
            fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
            box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los, **kwargs,
        )


def test_injection_rejects_a_2d_field() -> None:
    """A lightcone is 3D; anything else is a mistake worth surfacing."""
    with pytest.raises(ValueError, match="3D"):
        inject_foreground(np.zeros((16, 16)), K_PERP, K_PAR, 7.0)


# ===========================================================================
#  Removal
# ===========================================================================

@pytest.fixture(scope="module")
def realisation(fg_sim: SimulationData) -> ForegroundRealisation:
    """A 100x contaminated version of the synthetic lightcone."""
    return inject_foreground(
        fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
        foreground_amplitude=100.0, seed=42,
        box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
        lc_redshifts=fg_sim.lc_redshifts,
    )


def test_removal_endpoints_are_exact(
    fg_sim: SimulationData, realisation: ForegroundRealisation
) -> None:
    """0 % changes nothing; 100 % recovers the clean field exactly."""
    none_removed = remove_foreground(
        realisation.contaminated, K_PERP, K_PAR, 0.0,
        foreground=realisation.foreground,
    )
    all_removed = remove_foreground(
        realisation.contaminated, K_PERP, K_PAR, 1.0,
        foreground=realisation.foreground,
    )

    assert np.allclose(none_removed, realisation.contaminated)
    assert np.allclose(all_removed, fg_sim.brightness_temp_field)


@pytest.mark.parametrize("fraction", [0.0, 0.25, 0.5, 0.9, 0.99, 1.0])
def test_removal_amplitude_basis_scales_linearly(
    realisation: ForegroundRealisation, fraction: float
) -> None:
    """Default basis: residual amplitude is ``(1 - f)`` of the injected one."""
    residual = remove_foreground(
        realisation.contaminated, K_PERP, K_PAR, fraction,
        foreground=realisation.foreground,
    ) - (realisation.contaminated - realisation.foreground)

    assert residual.std() == pytest.approx(
        (1.0 - fraction) * realisation.foreground_rms, rel=1e-9, abs=1e-12
    )


@pytest.mark.parametrize("fraction", [0.0, 0.5, 0.9, 0.99])
def test_removal_power_basis_scales_as_the_square_root(
    realisation: ForegroundRealisation, fraction: float
) -> None:
    """Power basis: residual *power* is ``(1 - f)`` of the injected power."""
    residual = remove_foreground(
        realisation.contaminated, K_PERP, K_PAR, fraction,
        foreground=realisation.foreground, removal_basis="power",
    ) - (realisation.contaminated - realisation.foreground)

    assert residual.var() == pytest.approx(
        (1.0 - fraction) * realisation.foreground.var(), rel=1e-9
    )


@pytest.mark.parametrize("fraction", [-0.1, 1.1, np.nan])
def test_removal_rejects_out_of_range_fractions(
    realisation: ForegroundRealisation, fraction: float
) -> None:
    """``removal_fraction`` is a fraction; anything else is a bug upstream."""
    with pytest.raises(ValueError, match="removal_fraction"):
        remove_foreground(
            realisation.contaminated, K_PERP, K_PAR, fraction,
            foreground=realisation.foreground,
        )


def test_removal_rejects_a_mismatched_foreground(
    realisation: ForegroundRealisation,
) -> None:
    """Subtracting a differently shaped template is always a mistake."""
    with pytest.raises(ValueError, match="shapes differ"):
        remove_foreground(
            realisation.contaminated, K_PERP, K_PAR, 0.5,
            foreground=np.zeros((4, 4, 4)),
        )


def test_removal_rejects_an_unknown_basis(
    realisation: ForegroundRealisation,
) -> None:
    """A typo in ``removal_basis`` must not silently pick a default."""
    with pytest.raises(ValueError, match="removal_basis"):
        remove_foreground(
            realisation.contaminated, K_PERP, K_PAR, 0.5,
            foreground=realisation.foreground, removal_basis="pca",
        )


# ===========================================================================
#  End to end, through the unmodified analysis chain
# ===========================================================================

def _spectra(field: np.ndarray, data: SimulationData):
    """Power spectra of ``field`` against the synthetic galaxy overdensity."""
    return analysis.compute_all_power_spectra(
        field, data.galaxy_overdensity,
        box_len_perp=data.BOX_LEN, box_len_los=data.L_los,
        n_bins_perp=8, n_bins_parallel=8,
    )


def test_contamination_inflates_the_21cm_auto_power(
    fg_sim: SimulationData, realisation: ForegroundRealisation
) -> None:
    """
    The foreground shows up in ``P_21``, and preferentially at low ``k_∥``.

    This is the channel by which it reaches the cross-correlation error bar.
    """
    clean = _spectra(fg_sim.brightness_temp_field, fg_sim)
    dirty = _spectra(realisation.contaminated, fg_sim)

    low = np.nanmean(dirty.P_21cm_auto[:, :2]) / np.nanmean(clean.P_21cm_auto[:, :2])
    high = np.nanmean(dirty.P_21cm_auto[:, -2:]) / np.nanmean(clean.P_21cm_auto[:, -2:])

    assert low > 10.0
    assert low > high


def test_removal_recovers_the_clean_snr_monotonically(
    fg_sim: SimulationData, realisation: ForegroundRealisation
) -> None:
    """
    Detectability degrades with contamination and is restored as it is removed.

    Runs the real chain — ``compute_all_power_spectra`` then
    ``compute_uncertainty_budget`` — at each removal level, with no
    re-implementation of either.
    """
    def total_snr(field: np.ndarray) -> float:
        return analysis.compute_uncertainty_budget(
            spectra=_spectra(field, fg_sim),
            z_obs=fg_sim.z_obs,
            photoz_uncertainty=fg_sim.get("photoz_uncertainty", 0.059),
        ).total_snr

    clean_snr = total_snr(fg_sim.brightness_temp_field)
    snrs = [
        total_snr(remove_foreground(
            realisation.contaminated, K_PERP, K_PAR, fraction,
            foreground=realisation.foreground,
        ))
        for fraction in (0.0, 0.5, 0.9, 0.99, 1.0)
    ]

    assert snrs[0] < clean_snr                       # contamination hurts
    assert all(a <= b for a, b in zip(snrs, snrs[1:]))   # removal helps
    assert snrs[-1] == pytest.approx(clean_snr, rel=1e-9)


def test_foregrounds_move_the_cross_power_less_than_the_auto_power(
    fg_sim: SimulationData, realisation: ForegroundRealisation
) -> None:
    """
    The cross-spectrum is the more robust statistic — but not immune.

    ``P_21`` picks up the foreground directly; ``P_21×gal`` picks it up only
    through a chance correlation with the galaxy field.  The second is far
    smaller than the first, which is the usual argument for cross-correlating.
    It is not zero, though — see the scaling test below.
    """
    clean = _spectra(fg_sim.brightness_temp_field, fg_sim)
    dirty = _spectra(realisation.contaminated, fg_sim)

    cross_shift = np.nanmedian(
        np.abs(dirty.P_cross - clean.P_cross) / np.abs(clean.P_cross)
    )
    auto_shift = np.nanmedian(
        np.abs(dirty.P_21cm_auto - clean.P_21cm_auto) / np.abs(clean.P_21cm_auto)
    )

    assert cross_shift < 0.05 * auto_shift


def test_spurious_cross_correlation_grows_linearly_with_amplitude(
    fg_sim: SimulationData,
) -> None:
    """
    "Unbiased in the mean" is not "zero in your data".

    A single realisation carries a chance cross-correlation of order
    ``sqrt(P_21 P_gal / N_modes)``, so it scales **linearly** with the
    foreground amplitude rather than staying put.  At high contamination it
    can exceed the true cross-power outright, which is why the notebook plots
    a signal-only SNR beside the as-measured one: the as-measured SNR has a
    contaminated numerator and flatters the result.

    Regression guard on the physics — if this ever came out flat or quadratic,
    the injected field would not be behaving as an uncorrelated contaminant.
    """
    clean = _spectra(fg_sim.brightness_temp_field, fg_sim)

    def excess(amplitude: float):
        """Median per-bin shift in each spectrum at this contamination level."""
        realisation = inject_foreground(
            fg_sim.brightness_temp_field, K_PERP, K_PAR, fg_sim.z_obs,
            foreground_amplitude=amplitude, seed=42,
            box_len_perp=fg_sim.BOX_LEN, box_len_los=fg_sim.L_los,
            lc_redshifts=fg_sim.lc_redshifts,
        )
        dirty = _spectra(realisation.contaminated, fg_sim)
        return (
            float(np.nanmedian(np.abs(dirty.P_cross - clean.P_cross))),
            float(np.nanmedian(np.abs(dirty.P_21cm_auto - clean.P_21cm_auto))),
        )

    cross_low, auto_low = excess(1e3)
    cross_high, auto_high = excess(1e4)

    # The cross term is T_fg x delta_gal — one power of the foreground.
    assert cross_low > 0
    assert cross_high / cross_low == pytest.approx(10.0, rel=0.05)

    # The auto term is the foreground's own power — two powers of it.  The
    # contrast between the two exponents *is* the statement that foregrounds
    # reach the cross-spectrum only through a chance correlation.
    assert auto_high / auto_low == pytest.approx(100.0, rel=0.05)
