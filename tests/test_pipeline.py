#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``run_pipeline.py`` — the end-to-end driver.

These never invoke 21cmFAST: the simulation stage is either skipped or, for
the ``--sim force`` path, pointed at a stub script.
"""

from __future__ import annotations

import json
import os
import sys

import h5py
import numpy as np
import pytest

import run_pipeline

from conftest import TINY_N_Z, write_tiny_simulation


@pytest.fixture()
def workspace(tmp_path):
    """
    A scratch run directory with its own copy of the synthetic simulation.

    Returns
    -------
    dict
        Paths for ``data``, ``products``, ``figdir``, and ``summary``.
    """
    data = str(tmp_path / "lightcone_data.h5")
    write_tiny_simulation(data)
    return {
        "data": data,
        "products": str(tmp_path / "analysis_products.h5"),
        "figdir": str(tmp_path / "figures"),
        "summary": str(tmp_path / "summary.json"),
    }


def base_args(workspace, *extra: str) -> list:
    """
    Build a CLI argument list pointing at the scratch workspace.

    Parameters
    ----------
    workspace : dict
        Fixture output.
    *extra
        Additional command-line flags.

    Returns
    -------
    list of str
        Argument list for :func:`run_pipeline.main`.
    """
    return [
        "--data", workspace["data"],
        "--products", workspace["products"],
        "--figdir", workspace["figdir"],
        "--summary", workspace["summary"],
        "--sim", "skip",
        *extra,
    ]


# ===========================================================================
#  CLI plumbing
# ===========================================================================

def test_default_arguments() -> None:
    """Defaults are the safe ones: never re-run the simulation blindly."""
    args = run_pipeline.parse_args([])
    assert args.sim == "auto"
    assert args.analysis == "auto"
    assert args.plots == ["all"]
    assert args.max_halos == 0


def test_resolve_plot_groups() -> None:
    """``all`` expands to every group, ``none`` to nothing, names to themselves."""
    assert run_pipeline.resolve_plot_groups(["all"]) == list(run_pipeline.PLOT_GROUPS)
    assert run_pipeline.resolve_plot_groups(["none"]) == []
    assert run_pipeline.resolve_plot_groups(["snr", "power"]) == ["power", "snr"]
    # 'none' wins over any other selection.
    assert run_pipeline.resolve_plot_groups(["none", "power"]) == []


# ===========================================================================
#  Stage behaviour
# ===========================================================================

def test_simulation_stage_skips_when_output_exists(workspace) -> None:
    """``--sim auto`` does not re-run when the HDF5 is already there."""
    ran = run_pipeline.run_simulation_stage(
        "auto", workspace["data"], "run_simulation.py", quiet=True,
    )
    assert ran is False


def test_simulation_stage_errors_when_skipped_without_data(tmp_path) -> None:
    """``--sim skip`` with no stored output is a clear error."""
    with pytest.raises(FileNotFoundError, match="no stored simulation output"):
        run_pipeline.run_simulation_stage(
            "skip", str(tmp_path / "missing.h5"), "run_simulation.py", quiet=True,
        )


def test_simulation_stage_runs_stub_script(tmp_path, workspace) -> None:
    """``--sim force`` executes the configured script."""
    stub = tmp_path / "stub_sim.py"
    marker = tmp_path / "ran.txt"
    stub.write_text(f"open({str(marker)!r}, 'w').write('ok')\n")

    ran = run_pipeline.run_simulation_stage(
        "force", workspace["data"], str(stub), quiet=True,
    )

    assert ran is True
    assert marker.exists()


def test_simulation_stage_propagates_failure(tmp_path, workspace) -> None:
    """A failing simulation script aborts the pipeline."""
    stub = tmp_path / "failing_sim.py"
    stub.write_text("import sys; sys.exit(3)\n")

    with pytest.raises(RuntimeError, match="exit code 3"):
        run_pipeline.run_simulation_stage(
            "force", workspace["data"], str(stub), quiet=True,
        )


def test_analysis_skip_without_cache_errors(workspace) -> None:
    """``--analysis skip`` with no cache is a clear error, not a silent rerun."""
    assert run_pipeline.main(base_args(workspace, "--analysis", "skip")) == 1


# ===========================================================================
#  End-to-end
# ===========================================================================

def test_full_run_writes_figures_products_and_summary(workspace) -> None:
    """A complete run produces every figure, the cache, and the summary JSON."""
    exit_code = run_pipeline.main(base_args(workspace, "--analysis", "force"))

    assert exit_code == 0
    assert os.path.exists(workspace["products"])
    assert os.path.exists(workspace["summary"])

    figure_files = sorted(os.listdir(workspace["figdir"]))
    for expected in (
        "lightcone_fields.png",
        "lightcone_slice.png",
        "halo_catalogue.png",
        "sfr_relations.png",
        "uv_luminosity_function.png",
        "stellar_mass_muv.png",
        "main_sequence.png",
        "euclid_selected_catalogue.png",
        "selected_galaxy_overdensity.png",
        "galaxy_overdensity_on_21cm.png",
        "power_spectra_2d.png",
        "cross_snr.png",
    ):
        assert expected in figure_files

    with open(workspace["summary"]) as f:
        summary = json.load(f)

    assert summary["ran_simulation"] is False
    assert summary["recomputed_power_spectra"] is True
    assert summary["simulation"]["N_z"] == 12
    assert summary["observation"]["total_snr_sigma"] >= 0.0
    assert summary["power_spectra"]["large_scale_anticorrelated"] is True
    assert len(summary["figures"]) == len(figure_files)


def test_cached_run_reuses_products(workspace) -> None:
    """A second run with ``--analysis auto`` loads the cache instead of recomputing."""
    assert run_pipeline.main(base_args(workspace, "--analysis", "force")) == 0
    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0

    with open(workspace["summary"]) as f:
        summary = json.load(f)

    assert summary["recomputed_power_spectra"] is False
    assert summary["figures"] == []


def test_plot_selection_limits_output(workspace) -> None:
    """``--plots power snr`` writes only the k-space figures and their companions."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "power", "snr")
    ) == 0

    assert sorted(os.listdir(workspace["figdir"])) == [
        "cross_snr.png", "galaxy_wedge.png", "power_spectra_1d.png",
        "power_spectra_2d.png", "wedge_real_space.png",
    ]


def test_max_halos_is_recorded_in_summary(workspace) -> None:
    """Catalogue subsampling is reported so downstream densities can be checked."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--max-halos", "500",
                  "--plots", "scaling")
    ) == 0

    with open(workspace["summary"]) as f:
        summary = json.load(f)

    assert summary["simulation"]["halo_sampling_factor"] > 1.0
    assert summary["simulation"]["n_halos_total"] == 4_000


def test_budget_overrides_default_to_the_stored_attributes(workspace) -> None:
    """With no flags, the budget uses the HDF5 attributes verbatim."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "none")
    ) == 0

    with open(workspace["summary"]) as f:
        budget = json.load(f)["uncertainty_budget"]

    # Fixture attributes: sigma_z = 0.059, wedge_buffer = 0.0677.
    assert budget["photoz_uncertainty_sigma_z"] == pytest.approx(0.059)
    assert budget["wedge_buffer_Mpc-1"] == pytest.approx(0.0677)


def test_sigma_z_override_changes_the_budget(workspace) -> None:
    """``--sigma-z`` overrides the stored attribute and damps harder."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "none",
                  "--sigma-z", "0.45")
    ) == 0

    with open(workspace["summary"]) as f:
        budget = json.load(f)["uncertainty_budget"]

    assert budget["photoz_uncertainty_sigma_z"] == pytest.approx(0.45)
    # sigma_r = c sigma_z / H(z) grows in proportion to sigma_z.
    assert budget["radial_smearing_Mpc"] == pytest.approx(157.5, rel=1e-2)


def test_wedge_buffer_override_shrinks_the_usable_region(workspace) -> None:
    """A larger ``--wedge-buffer`` leaves strictly fewer usable modes."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "none",
                  "--wedge-buffer", "0.0")
    ) == 0
    with open(workspace["summary"]) as f:
        loose = json.load(f)["uncertainty_budget"]

    assert run_pipeline.main(
        base_args(workspace, "--plots", "none", "--wedge-buffer", "0.5")
    ) == 0
    with open(workspace["summary"]) as f:
        strict = json.load(f)["uncertainty_budget"]

    assert strict["modes_outside_wedge"] < loose["modes_outside_wedge"]
    assert strict["wedge_buffer_Mpc-1"] == pytest.approx(0.5)


def test_instrument_overrides_change_the_thermal_noise(workspace) -> None:
    """``--integration-time`` feeds through to P_N,21 as 1/t."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "none",
                  "--integration-time", "3.6e6")
    ) == 0
    with open(workspace["summary"]) as f:
        short = json.load(f)["uncertainty_budget"]

    assert run_pipeline.main(
        base_args(workspace, "--plots", "none", "--integration-time", "7.2e6")
    ) == 0
    with open(workspace["summary"]) as f:
        long = json.load(f)["uncertainty_budget"]

    assert long["P_noise_21cm"] == pytest.approx(short["P_noise_21cm"] / 2)


def test_budget_is_cached_alongside_the_spectra(workspace) -> None:
    """The products file carries the budget group after a run."""
    from src.dataio import load_uncertainty_budget

    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "none")
    ) == 0

    maps, attrs = load_uncertainty_budget(workspace["products"])

    assert maps["outside_wedge"].dtype == bool
    assert "sigma_cross" in maps and "cosmic_variance_term" in maps
    assert attrs["total_snr_sigma"] >= 0.0


def test_budget_figure_is_written(workspace) -> None:
    """``--plots budget`` renders the budget figure and the photo-z sweep."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "budget")
    ) == 0

    assert sorted(os.listdir(workspace["figdir"])) == [
        "photoz_suppression.png", "uncertainty_budget.png",
    ]


def test_missing_data_file_returns_error_code(tmp_path) -> None:
    """A missing simulation file exits non-zero rather than raising."""
    exit_code = run_pipeline.main([
        "--data", str(tmp_path / "nope.h5"),
        "--sim", "skip",
        "--plots", "none",
    ])
    assert exit_code == 1


def test_pdf_output_format(workspace) -> None:
    """The figure format flag is honoured."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "power",
                  "--format", "pdf")
    ) == 0
    assert sorted(os.listdir(workspace["figdir"])) == [
        "galaxy_wedge.pdf", "power_spectra_1d.pdf", "power_spectra_2d.pdf",
        "wedge_real_space.pdf",
    ]


def test_euclid_plot_group_is_written(workspace) -> None:
    """``--plots euclid`` writes only the three post-Euclid-cut figures."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "euclid")
    ) == 0

    assert sorted(os.listdir(workspace["figdir"])) == [
        "euclid_selected_catalogue.png",
        "galaxy_overdensity_on_21cm.png",
        "selected_galaxy_overdensity.png",
    ]


def test_euclid_group_honours_galaxy_weighting(workspace) -> None:
    """``--galaxy-weighting luminosity`` reaches the overdensity figures."""
    assert run_pipeline.main(
        base_args(
            workspace, "--analysis", "force", "--plots", "euclid",
            "--galaxy-weighting", "luminosity",
        )
    ) == 0

    assert os.path.exists(
        os.path.join(workspace["figdir"], "selected_galaxy_overdensity.png")
    )


# ===========================================================================
#  Run provenance
# ===========================================================================

def test_summary_names_the_source_simulation_run(workspace) -> None:
    """
    The summary points back at the simulation run that made its data.

    ``pipeline_summary.json`` is overwritten every run; the manifest it names
    is not, so an analysis-only run can still say where its numbers came from.
    Absent from HDF5 files written before manifests existed, hence the Nones.
    """
    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0

    with open(workspace["summary"]) as f:
        summary = json.load(f)

    assert set(summary["source_run"]) == {
        "run_id", "run_manifest", "random_seed", "n_threads",
    }


def test_simulation_subprocess_runs_unbuffered(workspace, monkeypatch) -> None:
    """
    The simulation child is always launched with ``-u``.

    Its stdout is block-buffered when the pipeline's own output is redirected
    to a file, and a child killed by a signal never flushes — which is how the
    2026-08-20 SIGSEGV produced a log that could not name the failing stage.
    """
    captured = {}

    class _Result:
        returncode = 0

    def fake_run(command, **kwargs):
        captured["command"] = command
        return _Result()

    monkeypatch.setattr(run_pipeline.subprocess, "run", fake_run)
    run_pipeline.run_simulation_stage(
        mode="force", data_path=workspace["data"],
        script="run_simulation.py", quiet=True,
    )

    assert captured["command"][1] == "-u"
    assert captured["command"][2] == "run_simulation.py"


def test_summary_is_valid_json_with_provenance_attrs(workspace) -> None:
    """
    Numpy scalars from the HDF5 must not break the summary.

    ``h5py`` returns ``np.int64`` / ``np.bool_`` for the ``random_seed`` and
    ``n_threads`` attributes that ``run_simulation.py`` now writes.
    ``json.dump`` cannot serialise those and fails *partway through*, leaving
    a truncated file that only reveals itself when something parses it.  The
    original ``source_run`` change slipped past because the test fixture's
    HDF5 predates those attributes.
    """
    import h5py

    with h5py.File(workspace["data"], "a") as f:
        f.attrs["run_id"] = "20260821_120000"
        f.attrs["run_manifest"] = "/tmp/sim_20260821_120000.json"
        f.attrs["random_seed"] = np.int64(42)
        f.attrs["n_threads"] = np.int64(8)

    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0

    with open(workspace["summary"]) as f:
        summary = json.load(f)          # would raise on a truncated write

    assert summary["source_run"]["random_seed"] == 42
    assert isinstance(summary["source_run"]["n_threads"], int)
    assert summary["source_run"]["run_id"] == "20260821_120000"


def test_noise_model_and_mode_weighting_reach_the_summary(workspace) -> None:
    """Both new flags are recorded, so a stored result says how it was made."""
    assert run_pipeline.main(base_args(
        workspace, "--plots", "none", "--noise-model", "physical",
        "--mode-weighted",
    )) == 0

    with open(workspace["summary"]) as f:
        budget = json.load(f)["uncertainty_budget"]

    assert budget["noise_model"] == "physical"
    assert budget["mode_weighted"] is True
    # k_perp-resolved noise is summarised by its finite range.
    assert budget["P_noise_21cm"] < budget["P_noise_21cm_max"]


def test_defaults_are_unchanged_by_the_new_options(workspace) -> None:
    """
    The default run must still be the historical one.

    Every stored figure, summary and note in this repo predates the physical
    noise model and the mode weighting; if the defaults moved, all of them
    would be silently wrong.
    """
    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0
    with open(workspace["summary"]) as f:
        budget = json.load(f)["uncertainty_budget"]

    assert budget["noise_model"] == "scaling"
    assert budget["mode_weighted"] is False
    assert budget["P_noise_21cm"] == pytest.approx(3.749, rel=1e-3)


# ===========================================================================
#  Estimator toggle — TODO.md P0
# ===========================================================================

def test_estimator_defaults_to_coeval(workspace) -> None:
    """Data written before P0 has no `estimator` attribute; assume coeval."""
    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0

    summary = json.loads(open(workspace["summary"]).read())
    assert summary["estimator"] == "coeval"
    assert "subbands" not in summary

    with h5py.File(workspace["products"], "r") as f:
        assert "subbands" not in f


def test_lightcone_estimator_runs_end_to_end(workspace) -> None:
    """The P0 formalism produces a sub-band cache and a per-band summary."""
    assert run_pipeline.main(
        base_args(workspace, "--plots", "none", "--estimator", "lightcone")
    ) == 0

    summary = json.loads(open(workspace["summary"]).read())
    assert summary["estimator"] == "lightcone"

    sub = summary["subbands"]
    assert sub["n_bands"] >= 1
    for key in ("z_effective", "bandwidth_MHz", "los_length_Mpc",
                "n_slices", "total_snr_per_band"):
        assert len(sub[key]) == sub["n_bands"]
    assert sub["combined_total_snr"] >= 0.0
    # Bands tile the line of sight without losing a slice.
    assert sum(sub["n_slices"]) == TINY_N_Z

    with h5py.File(workspace["products"], "r") as f:
        assert f.attrs["estimator"] == "lightcone"
        assert int(f["subbands"].attrs["n_bands"]) == sub["n_bands"]


def test_lightcone_estimator_changes_nothing_when_not_requested(workspace) -> None:
    """Every P0 addition is opt-in: the default totals must not move."""
    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0
    baseline = json.loads(open(workspace["summary"]).read())

    assert run_pipeline.main(
        base_args(workspace, "--plots", "none", "--analysis", "force")
    ) == 0
    repeat = json.loads(open(workspace["summary"]).read())

    assert (baseline["uncertainty_budget"]["total_snr_sigma"]
            == repeat["uncertainty_budget"]["total_snr_sigma"])
    assert repeat["estimator"] == "coeval"


def test_estimator_follows_the_stored_simulation(workspace) -> None:
    """`--estimator auto` reads the attribute run_simulation.py wrote."""
    with h5py.File(workspace["data"], "a") as f:
        f.attrs["estimator"] = "lightcone"

    assert run_pipeline.main(base_args(workspace, "--plots", "none")) == 0
    summary = json.loads(open(workspace["summary"]).read())
    assert summary["estimator"] == "lightcone"
    assert "subbands" in summary


def test_explicit_estimator_overrides_the_stored_attribute(workspace) -> None:
    """An explicit flag beats the file, for re-analysing an existing run."""
    with h5py.File(workspace["data"], "a") as f:
        f.attrs["estimator"] = "lightcone"

    assert run_pipeline.main(
        base_args(workspace, "--plots", "none", "--estimator", "coeval")
    ) == 0
    summary = json.loads(open(workspace["summary"]).read())
    assert summary["estimator"] == "coeval"
    assert "subbands" not in summary
