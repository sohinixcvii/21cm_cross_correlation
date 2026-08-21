"""
Tests for the galaxy-overdensity weighting modes.

Covers ``analysis.deposit_halo_field`` and
``analysis.galaxy_overdensity_from_catalogue`` — the number- versus
luminosity-weighted constructions of delta_gal from the halo catalogue.
"""

from __future__ import annotations

import numpy as np
import pytest

from src.analysis import (
    GALAXY_WEIGHTING_MODES,
    deposit_halo_field,
    galaxy_overdensity_from_catalogue,
    select_euclid_halos,
)
from src.conversions import sfr_to_Luv


BOX_LEN = 100.0
N_PERP = 8


def _catalogue(n_halos: int = 4000, seed: int = 7):
    """
    Synthetic halo catalogue spanning the Euclid magnitude window.

    Returns
    -------
    tuple of ndarray
        ``(coords, sfr, halo_masses)``.
    """
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0.0, BOX_LEN, size=(n_halos, 3))
    # Wide log-normal SFR so a decent fraction lands inside -22 < M_UV < -18.
    sfr = 10.0 ** rng.normal(1.0, 0.8, size=n_halos)
    halo_masses = 10.0 ** rng.uniform(9.0, 12.0, size=n_halos)
    return coords, sfr, halo_masses


# ---------------------------------------------------------------------------
#  deposit_halo_field
# ---------------------------------------------------------------------------

def test_deposit_unit_weights_conserves_count() -> None:
    """Unweighted deposit must place every in-box halo exactly once."""
    coords, _, _ = _catalogue()
    field = deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP)
    assert field.shape == (N_PERP, N_PERP, N_PERP)
    assert field.sum() == pytest.approx(coords.shape[0])


def test_deposit_weights_conserve_total_weight() -> None:
    """Weighted deposit must conserve the summed weight."""
    coords, sfr, _ = _catalogue()
    weights = sfr_to_Luv(sfr)
    field = deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP, weights=weights)
    assert field.sum() == pytest.approx(weights.sum(), rel=1e-10)


def test_deposit_respects_non_cubic_grid() -> None:
    """n_los controls the LOS axis independently of n_perp."""
    coords, _, _ = _catalogue()
    field = deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP, n_los=13)
    assert field.shape == (N_PERP, N_PERP, 13)


def test_deposit_empty_catalogue_returns_zeros() -> None:
    """An empty catalogue yields a correctly shaped zero grid, not an error."""
    field = deposit_halo_field(
        np.empty((0, 3)), box_len=BOX_LEN, n_perp=N_PERP, n_los=5,
    )
    assert field.shape == (N_PERP, N_PERP, 5)
    assert not field.any()


def test_deposit_rejects_bad_shapes() -> None:
    """Malformed coords or mismatched weights must raise ValueError."""
    with pytest.raises(ValueError):
        deposit_halo_field(np.zeros((10, 2)), box_len=BOX_LEN, n_perp=N_PERP)
    with pytest.raises(ValueError):
        deposit_halo_field(
            np.zeros((10, 3)), box_len=BOX_LEN, n_perp=N_PERP,
            weights=np.ones(9),
        )


# ---------------------------------------------------------------------------
#  galaxy_overdensity_from_catalogue
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("weighting", GALAXY_WEIGHTING_MODES)
def test_overdensity_has_zero_mean(weighting: str) -> None:
    """Both modes must produce a field with vanishing mean, by construction."""
    coords, sfr, halo_masses = _catalogue()
    delta, selection = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        weighting=weighting,
    )
    assert selection.n_selected > 0
    assert delta.mean() == pytest.approx(0.0, abs=1e-12)
    assert delta.min() >= -1.0


@pytest.mark.parametrize("weighting", GALAXY_WEIGHTING_MODES)
def test_overdensity_shape_is_interchangeable(weighting: str) -> None:
    """The two modes are drop-in replacements: same grid, same shape."""
    coords, sfr, halo_masses = _catalogue()
    delta, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP, n_los=11,
        weighting=weighting,
    )
    assert delta.shape == (N_PERP, N_PERP, 11)


def test_number_weighted_matches_manual_formula() -> None:
    """delta_gal = N / <N> - 1 over the Euclid-selected halos."""
    coords, sfr, halo_masses = _catalogue()
    delta, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        weighting="number",
    )

    selection = select_euclid_halos(sfr, halo_masses)
    valid = np.isfinite(sfr) & (sfr > 0) & (halo_masses > 0)
    keep = np.zeros(coords.shape[0], dtype=bool)
    keep[np.flatnonzero(valid)[selection.mask]] = True

    counts = deposit_halo_field(coords[keep], box_len=BOX_LEN, n_perp=N_PERP)
    np.testing.assert_allclose(delta, counts / counts.mean() - 1.0)


def test_luminosity_weighted_matches_manual_formula() -> None:
    """delta_gal,L = sum(L_UV) / <sum(L_UV)> - 1 over the same halos."""
    coords, sfr, halo_masses = _catalogue()
    delta, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        weighting="luminosity",
    )

    selection = select_euclid_halos(sfr, halo_masses)
    valid = np.isfinite(sfr) & (sfr > 0) & (halo_masses > 0)
    keep = np.zeros(coords.shape[0], dtype=bool)
    keep[np.flatnonzero(valid)[selection.mask]] = True

    lum = deposit_halo_field(
        coords[keep], box_len=BOX_LEN, n_perp=N_PERP, weights=sfr_to_Luv(sfr[keep]),
    )
    np.testing.assert_allclose(delta, lum / lum.mean() - 1.0)


def test_luminosity_weighting_is_scale_invariant() -> None:
    """
    A constant rescaling of L_UV cannot change delta_gal,L.

    This is what makes the kappa_UV conversion irrelevant to the field: it
    divides out of the ratio.  It also means an SFR-weighted field is the
    same field.
    """
    coords, sfr, halo_masses = _catalogue()
    delta_lum, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        weighting="luminosity",
    )

    keep = select_euclid_halos(sfr, halo_masses).mask
    valid = np.isfinite(sfr) & (sfr > 0) & (halo_masses > 0)
    full = np.zeros(coords.shape[0], dtype=bool)
    full[np.flatnonzero(valid)[keep]] = True

    sfr_field = deposit_halo_field(
        coords[full], box_len=BOX_LEN, n_perp=N_PERP, weights=sfr[full],
    )
    np.testing.assert_allclose(delta_lum, sfr_field / sfr_field.mean() - 1.0, atol=1e-12)


def test_modes_differ_when_sfr_is_not_uniform() -> None:
    """
    With a spread in SFR the two weightings must give genuinely different
    fields — otherwise the option would be pointless.
    """
    coords, sfr, halo_masses = _catalogue()
    delta_n, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP, weighting="number",
    )
    delta_l, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP, weighting="luminosity",
    )
    assert not np.allclose(delta_n, delta_l)
    # Luminosity weighting up-weights bright halos, so it is the noisier field.
    assert delta_l.std() > delta_n.std()


def test_selection_can_be_disabled() -> None:
    """apply_selection=False keeps every halo with SFR > 0."""
    coords, sfr, halo_masses = _catalogue()
    delta_all, _ = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        apply_selection=False,
    )
    delta_sel, selection = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
    )
    assert selection.n_selected < selection.n_valid
    assert not np.allclose(delta_all, delta_sel)


def test_empty_selection_returns_zero_field() -> None:
    """An impossible magnitude window gives a zero field, not a divide error."""
    coords, sfr, halo_masses = _catalogue()
    delta, selection = galaxy_overdensity_from_catalogue(
        coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        M_UV_faint=-30.0, M_UV_bright=-31.0,
    )
    assert selection.n_selected == 0
    assert not delta.any()


def test_unknown_weighting_rejected() -> None:
    """An unrecognised mode must fail loudly rather than silently default."""
    coords, sfr, halo_masses = _catalogue()
    with pytest.raises(ValueError, match="weighting must be one of"):
        galaxy_overdensity_from_catalogue(
            coords, sfr, halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
            weighting="mass",
        )


def test_catalogue_length_mismatch_rejected() -> None:
    """Inconsistent catalogue arrays must raise rather than mis-align."""
    coords, sfr, halo_masses = _catalogue()
    with pytest.raises(ValueError, match="catalogue length mismatch"):
        galaxy_overdensity_from_catalogue(
            coords, sfr[:-1], halo_masses, box_len=BOX_LEN, n_perp=N_PERP,
        )


# ---------------------------------------------------------------------------
#  Coeval cubic deposit — the notebook's UV-map path
# ---------------------------------------------------------------------------
# 21cmfast_HERAxEuclid_lightcone.ipynb §3 builds its UV luminosity and
# selected-galaxy maps from the *coeval* perturbed halo catalogue, so those
# grids are cubic (HII_DIM, HII_DIM, HII_DIM) — deliberately not the
# (HII_DIM, HII_DIM, N_z) lightcone shape used for the power spectra.  Those
# cells used to hand-roll the deposit with floor/mod + np.add.at; these tests
# pin the switch to deposit_halo_field.


def _floor_mod_deposit(coords, box_len, n_perp, weights=None):
    """
    The hand-rolled deposit the notebook used before ``deposit_halo_field``.

    Kept only as a reference implementation for the equivalence test below.
    """
    cell = box_len / n_perp
    idx = np.mod(np.floor(coords / cell).astype(int), n_perp)
    field = np.zeros((n_perp, n_perp, n_perp), dtype=float)
    np.add.at(
        field,
        (idx[:, 0], idx[:, 1], idx[:, 2]),
        1.0 if weights is None else weights,
    )
    return field


def test_cubic_deposit_matches_hand_rolled_floor_mod() -> None:
    """For in-box coords the helper reproduces the old floor/mod deposit exactly."""
    coords, sfr, _ = _catalogue()
    weights = sfr_to_Luv(sfr)

    expected = _floor_mod_deposit(coords, BOX_LEN, N_PERP, weights)
    actual = deposit_halo_field(
        coords, box_len=BOX_LEN, n_perp=N_PERP, weights=weights
    )
    assert actual.shape == (N_PERP, N_PERP, N_PERP)
    assert np.allclose(actual, expected, rtol=1e-12)


def test_cubic_deposit_counts_match_hand_rolled() -> None:
    """Same equivalence for the unweighted count grid."""
    coords, _, _ = _catalogue()
    expected = _floor_mod_deposit(coords, BOX_LEN, N_PERP)
    actual = deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP)
    assert np.array_equal(actual, expected)


def test_cubic_deposit_drops_rather_than_wraps_out_of_box_halos() -> None:
    """
    The one behavioural difference from the old path, and why the notebook
    now prints a warning: floor/mod wrapped out-of-box halos back inside
    periodically, whereas histogramdd drops them.
    """
    inside = np.array([[10.0, 10.0, 10.0]])
    outside = np.array([[BOX_LEN + 5.0, 10.0, 10.0]])
    coords = np.vstack([inside, outside])

    assert deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP).sum() == 1.0
    # The old path kept both, wrapping the stray halo around the box.
    assert _floor_mod_deposit(coords, BOX_LEN, N_PERP).sum() == 2.0


def test_cubic_deposit_is_separate_from_lightcone_shape() -> None:
    """
    The coeval map grid and the lightcone spectra grid are different shapes;
    the helper produces whichever is asked for, so they cannot be conflated.
    """
    coords, _, _ = _catalogue()
    n_z = 13
    cubic = deposit_halo_field(coords, box_len=BOX_LEN, n_perp=N_PERP)
    lightcone = deposit_halo_field(
        coords, box_len=BOX_LEN, n_perp=N_PERP, n_los=n_z, los_extent=BOX_LEN,
    )
    assert cubic.shape == (N_PERP, N_PERP, N_PERP)
    assert lightcone.shape == (N_PERP, N_PERP, n_z)
    # Same halos either way — only the LOS binning differs.
    assert cubic.sum() == pytest.approx(lightcone.sum())
