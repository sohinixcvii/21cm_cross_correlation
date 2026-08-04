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

import pytest

import run_pipeline

from conftest import write_tiny_simulation


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
    """``--plots power snr`` writes only the two k-space figures."""
    assert run_pipeline.main(
        base_args(workspace, "--analysis", "force", "--plots", "power", "snr")
    ) == 0

    assert sorted(os.listdir(workspace["figdir"])) == [
        "cross_snr.png", "power_spectra_2d.png",
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
    assert os.listdir(workspace["figdir"]) == ["power_spectra_2d.pdf"]
