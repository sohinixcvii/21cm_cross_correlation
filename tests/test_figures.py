#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/figures.py``.

Every plotting function is exercised on the synthetic simulation: it must
build a figure and write a non-empty file without touching a display.
"""

from __future__ import annotations

import os
from dataclasses import replace

import matplotlib
import numpy as np
import pytest
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from src import analysis, figures
from src.dataio import SimulationData


@pytest.fixture(scope="module")
def spectra(tiny_sim: SimulationData):
    """Power spectra of the synthetic simulation."""
    return analysis.compute_all_power_spectra(
        tiny_sim.brightness_temp_field,
        tiny_sim.galaxy_overdensity,
        box_len_perp=tiny_sim.BOX_LEN,
        box_len_los=tiny_sim.L_los,
        n_bins_perp=8,
        n_bins_parallel=8,
    )


@pytest.fixture(scope="module")
def wedge_slopes(tiny_sim: SimulationData):
    """Horizon and FoV wedge slopes at the reference redshift."""
    return (
        analysis.horizon_wedge_slope(tiny_sim.z_obs),
        analysis.fov_wedge_slope(tiny_sim.z_obs),
    )


def assert_saves(fig: Figure, tmp_path, name: str) -> None:
    """
    Assert that a figure is a real Figure and writes a non-empty file.

    Parameters
    ----------
    fig : Figure
        Figure under test.
    tmp_path : pathlib.Path
        pytest temporary directory.
    name : str
        Base filename.
    """
    assert isinstance(fig, Figure)
    path = figures.save_figure(fig, str(tmp_path), name)
    assert os.path.exists(path)
    assert os.path.getsize(path) > 1000


# ===========================================================================
#  Helpers
# ===========================================================================

def test_backend_is_headless() -> None:
    """Importing the module must force a non-interactive backend."""
    assert matplotlib.get_backend().lower() == "agg"


def test_apply_plot_style_sets_dpi() -> None:
    """The style helper propagates the requested DPI to rcParams."""
    figures.apply_plot_style(dpi=123)
    assert matplotlib.rcParams["figure.dpi"] == 123
    assert matplotlib.rcParams["savefig.dpi"] == 123
    figures.apply_plot_style()


def test_fill_nan_nearest_fills_holes() -> None:
    """Isolated NaNs are replaced by the neighbourhood mean."""
    arr = np.arange(9, dtype=float).reshape(3, 3)
    arr[1, 1] = np.nan

    filled = figures.fill_nan_nearest(arr)

    assert not np.isnan(filled).any()
    assert filled[0, 0] == 0.0


def test_fill_nan_nearest_passes_through_clean_and_all_nan_arrays() -> None:
    """No NaNs, or all NaNs, both return the input unchanged."""
    clean = np.ones((3, 3))
    np.testing.assert_array_equal(figures.fill_nan_nearest(clean), clean)

    all_nan = np.full((3, 3), np.nan)
    assert np.isnan(figures.fill_nan_nearest(all_nan)).all()


def test_eor_colormap_endpoints() -> None:
    """The EoR colormap runs from near-black to near-white."""
    cmap = figures.eor_colormap()
    assert cmap.name == "EoR21"
    assert sum(cmap(0.0)[:3]) < 0.5
    assert sum(cmap(1.0)[:3]) > 2.5


def test_save_figure_creates_directory(tmp_path) -> None:
    """A missing output directory is created rather than raising."""
    fig = Figure()
    fig.add_subplot(111).plot([0, 1], [0, 1])

    target = tmp_path / "nested" / "figures"
    path = figures.save_figure(fig, str(target), "line", fmt="pdf")

    assert path.endswith("line.pdf")
    assert os.path.exists(path)


# ===========================================================================
#  Part 2 figures
# ===========================================================================

def test_plot_lightcone_fields(tiny_sim: SimulationData, tmp_path) -> None:
    """The three-panel field figure renders and saves."""
    assert_saves(figures.plot_lightcone_fields(tiny_sim), tmp_path, "fields")


def test_plot_lightcone_slice(tiny_sim: SimulationData, tmp_path) -> None:
    """The wide EoR lightcone slice renders and saves."""
    assert_saves(figures.plot_lightcone_slice(tiny_sim), tmp_path, "slice")


def test_plot_halo_catalogue(tiny_sim: SimulationData, tmp_path) -> None:
    """The halo catalogue overview renders and saves."""
    assert_saves(figures.plot_halo_catalogue(tiny_sim), tmp_path, "halos")


def test_plot_sfr_relations(tiny_sim: SimulationData, tmp_path) -> None:
    """The SFR scaling-relation panels render and save."""
    assert_saves(figures.plot_sfr_relations(tiny_sim), tmp_path, "sfr")


def test_plot_uv_luminosity_function(tiny_sim: SimulationData, tmp_path) -> None:
    """The UVLF figure renders and saves."""
    assert_saves(figures.plot_uv_luminosity_function(tiny_sim), tmp_path, "uvlf")


def test_plot_stellar_mass_muv(tiny_sim: SimulationData, tmp_path) -> None:
    """The stellar-mass–magnitude figure renders and saves."""
    assert_saves(figures.plot_stellar_mass_muv(tiny_sim), tmp_path, "smuv")


def test_plot_main_sequence(tiny_sim: SimulationData, tmp_path) -> None:
    """The star-forming main-sequence figure renders and saves."""
    pytest.importorskip("astropy")
    assert_saves(figures.plot_main_sequence(tiny_sim), tmp_path, "ms")


def test_uvlf_rescales_subsampled_catalogues(tiny_sim: SimulationData) -> None:
    """Subsampling must not shift the luminosity-function normalisation."""
    import copy

    subsampled = copy.copy(tiny_sim)
    subsampled.sfr = tiny_sim.sfr[::4]
    subsampled.stellar_masses = tiny_sim.stellar_masses[::4]
    subsampled.halo_masses = tiny_sim.halo_masses[::4]
    subsampled.halo_coords = tiny_sim.halo_coords[::4]
    subsampled.halo_sampling_factor = 4.0

    full_fig = figures.plot_uv_luminosity_function(tiny_sim)
    sub_fig = figures.plot_uv_luminosity_function(subsampled)

    full_points = full_fig.axes[0].collections
    sub_points = sub_fig.axes[0].collections
    assert len(full_points) == len(sub_points)

    # The plotted data points (errorbar markers) carry the normalisation.
    full_y = full_fig.axes[0].lines[2].get_ydata()
    sub_y = sub_fig.axes[0].lines[2].get_ydata()
    overlap = min(len(full_y), len(sub_y))
    assert overlap > 0
    ratio = np.nanmedian(sub_y[:overlap] / full_y[:overlap])
    assert 0.5 < ratio < 2.0


# ===========================================================================
#  Part 3 figures
# ===========================================================================

def test_plot_power_spectra(tiny_sim: SimulationData, spectra, wedge_slopes, tmp_path) -> None:
    """The 2D power-spectrum panels render and save."""
    horizon, fov = wedge_slopes
    fig = figures.plot_power_spectra(spectra, tiny_sim, horizon, fov)
    assert_saves(fig, tmp_path, "power")


def test_plot_snr(tiny_sim: SimulationData, spectra, wedge_slopes, tmp_path) -> None:
    """The SNR figure renders and saves."""
    horizon, fov = wedge_slopes

    smearing = analysis.radial_smearing_length(
        tiny_sim.get("photoz_uncertainty"), tiny_sim.z_obs
    )
    kernel = analysis.photoz_damping_kernel(spectra.k_parallel, smearing)
    p_cross_observed = spectra.P_cross * kernel

    snr = analysis.cross_power_snr(
        P_cross_observed=p_cross_observed,
        P_21cm_auto=spectra.P_21cm_auto,
        P_galaxy_observed=spectra.P_galaxy_auto * kernel ** 2,
        P_noise_21cm=analysis.hera_thermal_noise_power(tiny_sim.z_obs, 3.6e6, 8e6),
        P_noise_galaxy=1.0 / tiny_sim.get("mean_galaxy_density"),
        outside_wedge=analysis.foreground_wedge_mask(
            spectra.k_perp, spectra.k_parallel, horizon, 0.02
        ),
    )

    fig = figures.plot_snr(spectra, snr, p_cross_observed, tiny_sim, horizon, fov)
    assert_saves(fig, tmp_path, "snr")


def test_plot_uncertainty_budget(
    tiny_sim: SimulationData, spectra, tmp_path
) -> None:
    """The uncertainty-budget figure renders and saves."""
    budget = analysis.compute_uncertainty_budget(
        spectra=spectra,
        z_obs=tiny_sim.z_obs,
        photoz_uncertainty=tiny_sim.get("photoz_uncertainty"),
        wedge_buffer=tiny_sim.get("wedge_buffer"),
        mean_galaxy_density=tiny_sim.get("mean_galaxy_density"),
    )

    fig = figures.plot_uncertainty_budget(budget, tiny_sim)
    assert_saves(fig, tmp_path, "uncertainty_budget")


def test_plot_bias_diagnostic(tiny_sim: SimulationData, tmp_path) -> None:
    """The bias diagnostic renders and saves."""
    pytest.importorskip("hmf")

    selection = analysis.select_euclid_halos(
        tiny_sim.sfr, tiny_sim.halo_masses,
        M_UV_faint=tiny_sim.get("M_UV_limit"), M_UV_bright=-22.0,
    )
    if selection.n_selected == 0:
        pytest.skip("no synthetic halos in the Euclid window")

    bias = analysis.effective_galaxy_bias(selection, z_obs=tiny_sim.z_obs)
    assert_saves(figures.plot_bias_diagnostic(bias, tiny_sim.z_obs), tmp_path, "bias")


def test_plot_galaxy_wedge(
    tiny_sim: SimulationData, spectra, wedge_slopes, tmp_path
) -> None:
    """The filled-wedge galaxy panel renders, saves, and shades the wedge."""
    horizon, fov = wedge_slopes
    fig = figures.plot_galaxy_wedge(spectra, tiny_sim, horizon, fov)
    assert_saves(fig, tmp_path, "galaxy_wedge")

    # The excluded region must be drawn as a filled patch, not lines alone —
    # that is the whole point of this figure over plot_power_spectra's panel.
    ax = fig.axes[0]
    assert any(c.get_hatch() for c in ax.collections), "wedge region is not hatched"
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Wedge (excluded)" in labels
    assert "Horizon" in labels and "HERA FoV wedge" in labels


def test_plot_wedge_real_space(tiny_sim: SimulationData, wedge_slopes, tmp_path) -> None:
    """The real-space wedge figure renders and actually suppresses structure."""
    horizon, _ = wedge_slopes
    fig = figures.plot_wedge_real_space(tiny_sim, horizon)
    assert_saves(fig, tmp_path, "wedge_real_space")

    original, filtered = (fig.axes[0].images[0], fig.axes[1].images[0])

    # Removing modes can only take power out of the field.
    assert np.std(filtered.get_array()) < np.std(original.get_array())

    # Both panels must share a colour scale, or the comparison is meaningless.
    assert original.get_clim() == filtered.get_clim()

    # The excluded percentage belongs in the title.
    assert "% of modes excluded" in fig._suptitle.get_text()


def test_plot_wedge_real_space_buffer_removes_more(
    tiny_sim: SimulationData, wedge_slopes
) -> None:
    """A non-zero wedge buffer excludes at least as much as the bare line."""
    horizon, _ = wedge_slopes

    def excluded_percent(buffer: float) -> float:
        fig = figures.plot_wedge_real_space(tiny_sim, horizon, wedge_buffer=buffer)
        title = fig.axes[1].get_title()
        return float(title.split("(")[1].split("%")[0])

    assert excluded_percent(0.5) >= excluded_percent(0.0)


def _photoz_budget(tiny_sim: SimulationData, spectra):
    """Uncertainty budget for the synthetic simulation."""
    return analysis.compute_uncertainty_budget(
        spectra=spectra,
        z_obs=tiny_sim.z_obs,
        photoz_uncertainty=tiny_sim.get("photoz_uncertainty"),
        wedge_buffer=tiny_sim.get("wedge_buffer"),
        mean_galaxy_density=tiny_sim.get("mean_galaxy_density"),
    )


def test_plot_photoz_suppression(
    tiny_sim: SimulationData, spectra, tmp_path
) -> None:
    """The photo-z sweep renders, saves, and puts W = 1 at sigma_z = 0."""
    budget = _photoz_budget(tiny_sim, spectra)
    fig = figures.plot_photoz_suppression(budget, tiny_sim)
    assert_saves(fig, tmp_path, "photoz_suppression")

    ax = fig.axes[0]
    curves = [ln for ln in ax.lines if len(ln.get_xdata()) > 2]
    assert len(curves) >= 5

    # sigma_z = 0 is the spectroscopic limit: no damping anywhere.
    assert np.allclose(curves[0].get_ydata(), 1.0)

    # Larger sigma_z damps harder at every k_par, so the curves never cross.
    for lower, higher in zip(curves, curves[1:]):
        assert np.all(higher.get_ydata() <= lower.get_ydata() + 1e-12)


def test_plot_photoz_suppression_always_includes_adopted(
    tiny_sim: SimulationData, spectra
) -> None:
    """The adopted sigma_z is drawn even when the caller omits it."""
    budget = _photoz_budget(tiny_sim, spectra)
    fig = figures.plot_photoz_suppression(
        budget, tiny_sim, sigma_z_scenarios=(0.0, 0.01),
    )
    labels = " ".join(t.get_text() for t in fig.axes[0].get_legend().get_texts())
    assert "adopted" in labels
    assert f"{budget.photoz_uncertainty:g}" in labels


def test_photoz_suppression_sigma_r_matches_analysis(
    tiny_sim: SimulationData, spectra
) -> None:
    """The scaled sigma_r reproduces radial_smearing_length for every scenario."""
    budget = _photoz_budget(tiny_sim, spectra)
    scale = budget.radial_smearing / budget.photoz_uncertainty

    for sigma_z in (0.02, 0.1, 0.45):
        assert np.isclose(
            sigma_z * scale,
            analysis.radial_smearing_length(sigma_z, tiny_sim.z_obs),
        )


def test_plot_uv_selection_maps(tiny_sim: SimulationData, tmp_path) -> None:
    """The UV/selection maps render, save, and mark both magnitude cuts."""
    fig = figures.plot_uv_selection_maps(tiny_sim, M_UV_bright=-22.0)
    assert_saves(fig, tmp_path, "uv_selection_maps")

    # Two projected maps plus the magnitude histogram.
    assert len(fig.axes[0].images) == 1 and len(fig.axes[1].images) == 1

    cut_positions = sorted(
        line.get_xdata()[0] for line in fig.axes[2].lines
    )
    assert cut_positions == [-22.0, float(tiny_sim.get("M_UV_limit", -18.0))]


def test_plot_uv_selection_maps_matches_analysis_selection(
    tiny_sim: SimulationData
) -> None:
    """The map's selected count equals src.analysis's own selection."""
    selection = analysis.select_euclid_halos(
        tiny_sim.sfr, tiny_sim.halo_masses,
        M_UV_faint=tiny_sim.get("M_UV_limit", -18.0), M_UV_bright=-22.0,
    )
    fig = figures.plot_uv_selection_maps(tiny_sim, M_UV_bright=-22.0)

    # The counts map is built from the same mask, so it must sum to n_selected.
    counts = fig.axes[1].images[0].get_array()
    assert int(np.nansum(counts)) == selection.n_selected
    assert f"{selection.n_selected:,}" in fig._suptitle.get_text()


def test_plot_uv_selection_maps_handles_zero_sfr_halos(
    tiny_sim: SimulationData, tmp_path
) -> None:
    """
    Halos with SFR <= 0 must not break the selection indexing.

    ``select_euclid_halos`` returns a mask over the *valid* (SFR > 0) subset,
    not the full catalogue.  Indexing ``halo_coords`` with it directly is a
    length mismatch as soon as any halo has SFR <= 0.
    """
    sfr = np.asarray(tiny_sim.sfr, dtype=float).copy()
    sfr[::3] = 0.0                      # knock out a third of the catalogue
    zeroed = replace(tiny_sim, sfr=sfr)

    selection = analysis.select_euclid_halos(
        sfr, zeroed.halo_masses,
        M_UV_faint=zeroed.get("M_UV_limit", -18.0), M_UV_bright=-22.0,
    )
    assert selection.n_valid < sfr.size   # the guard is actually exercised

    fig = figures.plot_uv_selection_maps(zeroed, M_UV_bright=-22.0)
    counts = fig.axes[1].images[0].get_array()
    assert int(np.nansum(counts)) == selection.n_selected
    assert_saves(fig, tmp_path, "uv_selection_maps_zero_sfr")


# ===========================================================================
#  Post-Euclid-cut figures
# ===========================================================================

@pytest.fixture(scope="module")
def selected_field(tiny_sim: SimulationData):
    """Post-cut galaxy overdensity and the selection that produced it."""
    return figures.selected_galaxy_overdensity(tiny_sim, M_UV_bright=-22.0)


def test_selected_galaxy_overdensity_geometry(
    tiny_sim: SimulationData, selected_field
) -> None:
    """The rebuilt field matches run_simulation.py's catalogue-field grid."""
    delta_gal, selection = selected_field

    assert delta_gal.shape == (tiny_sim.HII_DIM, tiny_sim.HII_DIM, tiny_sim.N_z)
    assert delta_gal.min() >= -1.0 - 1e-9          # delta = N/<N> - 1
    assert np.isclose(delta_gal.mean(), 0.0, atol=1e-9)

    # It must apply the same window as the bias stage, not the stored field.
    expected = analysis.select_euclid_halos(
        tiny_sim.sfr, tiny_sim.halo_masses,
        M_UV_faint=tiny_sim.get("M_UV_limit", -18.0), M_UV_bright=-22.0,
    )
    assert selection.n_selected == expected.n_selected
    assert 0 < selection.n_selected < selection.n_valid


def test_plot_euclid_selected_catalogue(tiny_sim: SimulationData, tmp_path) -> None:
    """Galaxies, halo masses, and SFRs after the cut render and save."""
    fig = figures.plot_euclid_selected_catalogue(tiny_sim, M_UV_bright=-22.0)
    assert_saves(fig, tmp_path, "euclid_selected_catalogue")

    selection = analysis.select_euclid_halos(
        tiny_sim.sfr, tiny_sim.halo_masses,
        M_UV_faint=tiny_sim.get("M_UV_limit", -18.0), M_UV_bright=-22.0,
    )
    assert f"{selection.n_selected:,}" in fig._suptitle.get_text()

    # Panel 3 marks the SFR window equivalent to the magnitude cut.
    sfr_cuts = sorted(line.get_xdata()[0] for line in fig.axes[2].lines)
    assert np.allclose(
        sfr_cuts, sorted(np.log10([selection.SFR_min, selection.SFR_max])),
    )


def test_plot_selected_galaxy_overdensity(
    tiny_sim: SimulationData, selected_field, tmp_path
) -> None:
    """The post-cut overdensity figure renders both slices and the PDF."""
    delta_gal, selection = selected_field
    fig = figures.plot_selected_galaxy_overdensity(
        tiny_sim, delta_gal=delta_gal, selection=selection,
    )
    assert_saves(fig, tmp_path, "selected_galaxy_overdensity")

    assert len(fig.axes[0].images) == 1 and len(fig.axes[1].images) == 1
    assert f"{selection.n_selected:,}" in fig._suptitle.get_text()


def test_plot_galaxy_overdensity_on_21cm(
    tiny_sim: SimulationData, selected_field, tmp_path
) -> None:
    """The overlay renders contours over the 21 cm slice and reports r."""
    delta_gal, selection = selected_field
    fig = figures.plot_galaxy_overdensity_on_21cm(
        tiny_sim, delta_gal=delta_gal, selection=selection,
    )
    assert_saves(fig, tmp_path, "galaxy_overdensity_on_21cm")

    # Panel 1: 21 cm image with galaxy contours on top; panel 2 is the swap.
    assert len(fig.axes[0].images) == 1
    assert fig.axes[0].collections            # the contour set
    assert len(fig.axes[1].images) == 1
    assert fig.axes[1].collections
    assert "$r = " in fig.axes[2].get_title()


def test_overlay_survives_empty_selection(tiny_sim: SimulationData, tmp_path) -> None:
    """
    An empty selection gives an all-zero field, not a crash.

    ``galaxy_overdensity_from_catalogue`` returns zeros when nothing passes
    the cut, which leaves both the contour levels and the Pearson r
    degenerate; the figures must fall back rather than raise.
    """
    # Bright end fainter than the faint end: the window is empty by
    # construction, whatever the catalogue contains.
    delta_gal, selection = figures.selected_galaxy_overdensity(
        tiny_sim, M_UV_bright=-17.0,
    )
    assert selection.n_selected == 0
    assert np.all(delta_gal == 0.0)

    assert_saves(
        figures.plot_selected_galaxy_overdensity(
            tiny_sim, delta_gal=delta_gal, selection=selection,
        ),
        tmp_path, "selected_galaxy_overdensity_empty",
    )
    fig = figures.plot_galaxy_overdensity_on_21cm(
        tiny_sim, delta_gal=delta_gal, selection=selection,
    )
    assert "degenerate field" in fig.axes[2].get_title()
    assert_saves(fig, tmp_path, "galaxy_overdensity_on_21cm_empty")


# ===========================================================================
#  Lightcone vs coeval axis lengths
# ===========================================================================

def test_halo_catalogue_figure_handles_a_short_line_of_sight(tmp_path):
    """
    Regression guard for the 2026-08-25 smoke-run IndexError.

    ``plot_halo_catalogue`` indexed the lightcone's line-of-sight axis with
    ``HII_DIM // 2``.  That is an index into the *transverse* grid, and it
    stayed in range only while ``N_z > HII_DIM / 2``.  The smoke-test
    configuration (24 transverse cells, 12 slices) is the first geometry in
    this project where it does not.
    """
    from conftest import write_tiny_simulation
    from src.dataio import load_simulation

    path = write_tiny_simulation(
        str(tmp_path / "short_los.h5"), hii_dim=24, n_z=12, n_halos=500,
    )
    data = load_simulation(path)
    assert data.HII_DIM // 2 >= data.brightness_temp_field.shape[2], (
        "fixture no longer reproduces the failing geometry"
    )

    figure = figures.plot_halo_catalogue(data)
    assert figure is not None
    plt.close(figure)


def test_halo_catalogue_slice_is_taken_from_the_line_of_sight_axis(tiny_sim):
    """The plotted slice must come from the field's own third axis."""
    figure = figures.plot_halo_catalogue(tiny_sim)
    # Panel 3 holds the imshow; its data must be one transverse plane.
    image = figure.axes[2].images[0]
    n_z = tiny_sim.brightness_temp_field.shape[2]
    expected = tiny_sim.brightness_temp_field[:, :, n_z // 2].T
    assert np.allclose(image.get_array(), expected)
    plt.close(figure)


def test_every_figure_survives_the_smoke_test_geometry(tmp_path):
    """
    All 18 figures must render on the reduced grid, not just the first one.

    The 2026-08-25 smoke run failed on ``halo_catalogue``, the third figure
    written; the ones after it were never reached. This renders the whole set
    on a smoke-shaped box so a second latent shape assumption cannot hide
    behind the first.
    """
    from conftest import write_tiny_simulation
    from src.dataio import load_simulation
    import run_pipeline

    path = write_tiny_simulation(
        str(tmp_path / "smoke_shaped.h5"), hii_dim=24, n_z=12, n_halos=5_000,
    )
    data = load_simulation(path)

    spectra = analysis.compute_all_power_spectra(
        data.brightness_temp_field, data.galaxy_overdensity,
        box_len_perp=data.BOX_LEN, box_len_los=data.L_los,
        n_bins_perp=8, n_bins_parallel=8,
    )
    budget = analysis.compute_uncertainty_budget(spectra, z_obs=data.z_obs)
    _, bias = run_pipeline.bias_stage(data, m_uv_bright=-22.0, quiet=True)

    written = run_pipeline.figure_stage(
        groups=list(run_pipeline.PLOT_GROUPS),
        data=data, spectra=spectra, budget=budget, bias=bias,
        output_dir=str(tmp_path / "figures"),
        fmt="png", quiet=True, m_uv_bright=-22.0, galaxy_weighting="number",
    )

    # 18 before the 1D spectra and number-density figures were added.
    assert len(written) == 20
    assert all(os.path.exists(p) and os.path.getsize(p) > 0 for p in written)


def test_empty_bins_are_masked_not_fabricated() -> None:
    """
    Empty bins must not be filled with a neighbour's value.

    All three panels of ``power_spectra_2d`` share one k-grid and one set of
    mode counts, so their coverage is identical.  The auto-spectra used to be
    passed through ``fill_nan_nearest``, which drew empty bins as if measured
    and made the cross-power look uniquely sparse by comparison.
    """
    values = np.array([[1.0, np.nan], [100.0, 4.0]])
    masked = figures._masked_log10(values)

    assert np.ma.is_masked(masked)
    assert bool(np.ma.getmaskarray(masked)[0, 1])
    assert masked[1, 0] == pytest.approx(2.0)
    # the fabricating helper would have replaced the hole with a real number
    assert not np.isnan(figures.fill_nan_nearest(values)).any()
