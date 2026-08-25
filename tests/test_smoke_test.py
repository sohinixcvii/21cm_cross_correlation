#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test_smoke_test.py — the pre-flight checker must itself be checked
===================================================================

A smoke test that passes on broken output is worse than no smoke test, so
every ``check_*`` function is exercised here on both good and deliberately
corrupted inputs.

These are ordinary unit tests: they do **not** run 21cmFAST and do not invoke
the reduced simulation configuration. They verify the checker's logic and the
CLI wiring only.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

import run_pipeline
from src import smoke_test
from src.dataio import PowerSpectra, SubbandPowerSpectra


# ===========================================================================
#  The override table
# ===========================================================================

class TestOverrideTable:
    """The reduced configuration must be explicit, auditable and inert."""

    def test_every_override_names_its_production_value(self):
        for name, entry in smoke_test.SMOKE_TEST_OVERRIDES.items():
            assert "value" in entry, f"{name} has no smoke value"
            assert "production" in entry, f"{name} does not name what it replaces"
            assert entry["why"], f"{name} has no justification"

    def test_overrides_are_smaller_than_production(self):
        """Every reduction must actually reduce — max_halos 0 means 'all'."""
        for name, entry in smoke_test.SMOKE_TEST_OVERRIDES.items():
            if name == "max_halos":
                continue
            assert entry["value"] < entry["production"], (
                f"{name} is not a reduction"
            )

    def test_dim_keeps_the_21cmfast_convention(self):
        overrides = smoke_test.SMOKE_TEST_OVERRIDES
        assert overrides["DIM"]["value"] == 3 * overrides["HII_DIM"]["value"]

    def test_override_returns_production_value_for_unlisted_names(self):
        assert smoke_test.override("integration_time", 3.6e6) == 3.6e6
        assert smoke_test.override("HII_DIM", 256) == 16

    def test_describe_overrides_names_every_parameter(self):
        text = smoke_test.describe_overrides()
        for name in smoke_test.SMOKE_TEST_OVERRIDES:
            assert name in text
        assert "NOT a science run" in text

    def test_smoke_outputs_are_isolated_from_real_ones(self):
        assert smoke_test.SMOKE_OUTPUT_DIR != "outputs"
        assert smoke_test.SMOKE_OUTPUT_DIR.startswith("outputs")


# ===========================================================================
#  Report
# ===========================================================================

class TestSmokeReport:
    def test_passes_only_when_every_check_passes(self):
        report = smoke_test.SmokeReport()
        report.add("a", True, "fine")
        assert report.passed
        report.add("b", False, "broken")
        assert not report.passed

    def test_render_marks_failures(self):
        report = smoke_test.SmokeReport()
        report.add("stage", False, "wrong shape")
        text = report.render()
        assert "FAIL" in text and "wrong shape" in text
        assert "not physical" in text


# ===========================================================================
#  Stage checks
# ===========================================================================

class TestCheckSimulationOutput:
    def test_accepts_a_well_formed_simulation(self, tiny_sim):
        report = smoke_test.check_simulation_output(
            tiny_sim,
            expected_hii_dim=int(tiny_sim.HII_DIM),
            expected_n_z=int(tiny_sim.N_z),
        )
        assert report.passed, report.render()
        assert {c.stage for c in report.checks} == {
            "21 cm lightcone", "ionisation field",
            "lightcone geometry", "halo catalogue",
        }

    def test_rejects_the_wrong_grid_size(self, tiny_sim):
        report = smoke_test.check_simulation_output(
            tiny_sim,
            expected_hii_dim=int(tiny_sim.HII_DIM) + 1,
            expected_n_z=int(tiny_sim.N_z),
        )
        assert not report.passed

    def test_rejects_a_non_finite_field(self, tiny_sim):
        import copy
        broken = copy.copy(tiny_sim)
        broken.brightness_temp_field = tiny_sim.brightness_temp_field.copy()
        broken.brightness_temp_field[0, 0, 0] = np.nan

        report = smoke_test.check_simulation_output(
            broken,
            expected_hii_dim=int(tiny_sim.HII_DIM),
            expected_n_z=int(tiny_sim.N_z),
        )
        assert not report.passed
        assert any("finite" in c.detail for c in report.checks)

    def test_rejects_a_neutral_fraction_outside_zero_to_one(self, tiny_sim):
        import copy
        broken = copy.copy(tiny_sim)
        broken.neutral_fraction = tiny_sim.neutral_fraction.copy()
        broken.neutral_fraction[0, 0, 0] = 5.0

        report = smoke_test.check_simulation_output(
            broken,
            expected_hii_dim=int(tiny_sim.HII_DIM),
            expected_n_z=int(tiny_sim.N_z),
        )
        assert not report.passed


class TestCheckPowerSpectra:
    @staticmethod
    def _spectra(n_perp=8, n_par=8):
        rng = np.random.default_rng(1)
        return PowerSpectra(
            k_perp=np.linspace(0.1, 1.0, n_perp),
            k_parallel=np.linspace(0.1, 1.0, n_par),
            P_21cm_auto=rng.normal(size=(n_perp, n_par)),
            P_galaxy_auto=rng.normal(size=(n_perp, n_par)),
            P_cross=rng.normal(size=(n_perp, n_par)),
            mode_counts=np.ones((n_perp, n_par)),
        )

    def test_accepts_the_requested_binning(self):
        report = smoke_test.check_power_spectra(self._spectra(), 8, 8)
        assert report.passed, report.render()

    def test_rejects_the_wrong_binning(self):
        report = smoke_test.check_power_spectra(self._spectra(), 20, 20)
        assert not report.passed

    def test_rejects_an_entirely_empty_grid(self):
        spectra = self._spectra()
        spectra.mode_counts = np.zeros_like(spectra.mode_counts)
        report = smoke_test.check_power_spectra(spectra, 8, 8)
        assert not report.passed
        assert any("every bin is empty" in c.detail for c in report.checks)

    def test_checks_subbands_when_present(self):
        bands = [self._spectra() for _ in range(2)]
        subbands = SubbandPowerSpectra(
            bands=bands,
            z_effective=np.array([6.7, 7.3]),
            z_min=np.array([6.55, 7.05]),
            z_max=np.array([7.05, 7.45]),
            frequency_min_hz=np.array([1.7e8, 1.8e8]),
            frequency_max_hz=np.array([1.75e8, 1.85e8]),
            bandwidth_hz=np.array([7.0e6, 7.0e6]),
            los_length_mpc=np.array([150.0, 150.0]),
            n_slices=np.array([6, 6]),
            index_ranges=np.array([[0, 6], [6, 12]]),
        )
        report = smoke_test.check_power_spectra(
            bands[0], 8, 8, subbands=subbands
        )
        assert report.passed, report.render()
        assert any(c.stage == "sub-bands" for c in report.checks)


class TestCheckUncertaintyBudget:
    def test_accepts_a_real_budget(self, tiny_sim):
        from src import analysis
        spectra = analysis.compute_all_power_spectra(
            tiny_sim.brightness_temp_field, tiny_sim.galaxy_overdensity,
            box_len_perp=tiny_sim.BOX_LEN, box_len_los=tiny_sim.L_los,
            n_bins_perp=8, n_bins_parallel=8,
        )
        budget = analysis.compute_uncertainty_budget(spectra, z_obs=tiny_sim.z_obs)
        report = smoke_test.check_uncertainty_budget(budget, 8, 8)
        assert report.passed, report.render()

    def test_rejects_the_wrong_map_shape(self, tiny_sim):
        from src import analysis
        spectra = analysis.compute_all_power_spectra(
            tiny_sim.brightness_temp_field, tiny_sim.galaxy_overdensity,
            box_len_perp=tiny_sim.BOX_LEN, box_len_los=tiny_sim.L_los,
            n_bins_perp=8, n_bins_parallel=8,
        )
        budget = analysis.compute_uncertainty_budget(spectra, z_obs=tiny_sim.z_obs)
        report = smoke_test.check_uncertainty_budget(budget, 20, 20)
        assert not report.passed


class TestCheckSummaryAndFigures:
    def test_accepts_a_complete_summary(self):
        summary = {
            "generated": "now", "data_file": "x.h5", "simulation": {},
            "power_spectra": {}, "uncertainty_budget": {},
            "estimator": "coeval", "figures": [],
        }
        assert smoke_test.check_summary(summary).passed

    def test_rejects_a_summary_missing_a_block(self):
        report = smoke_test.check_summary({"generated": "now"})
        assert not report.passed
        assert "missing keys" in report.checks[0].detail

    def test_figures_skipped_is_not_a_failure(self):
        report = smoke_test.check_figures([])
        assert report.passed
        assert "skipped" in report.checks[0].detail

    def test_rejects_a_missing_figure_file(self, tmp_path):
        report = smoke_test.check_figures([str(tmp_path / "absent.png")])
        assert not report.passed

    def test_accepts_written_figures(self, tmp_path):
        path = tmp_path / "figure.png"
        path.write_bytes(b"not really a png, but non-empty")
        assert smoke_test.check_figures([str(path)]).passed


class TestCheckMcmcChain:
    """There is no sampler in this pipeline; the hook must say so."""

    def test_absent_stage_is_reported_not_silently_passed(self):
        report = smoke_test.check_mcmc_chain()
        assert report.passed
        assert "no sampler" in report.checks[0].detail

    def test_validates_a_chain_when_one_is_supplied(self):
        chain = np.zeros((4, 10, 3))
        report = smoke_test.check_mcmc_chain(
            chain, expected_walkers=4, expected_steps=10
        )
        assert report.passed

    def test_rejects_the_wrong_chain_shape(self):
        report = smoke_test.check_mcmc_chain(
            np.zeros((4, 10, 3)), expected_walkers=8, expected_steps=10
        )
        assert not report.passed

    def test_rejects_a_non_finite_chain(self):
        chain = np.zeros((2, 5, 2))
        chain[0, 0, 0] = np.inf
        assert not smoke_test.check_mcmc_chain(chain).passed


# ===========================================================================
#  CLI wiring
# ===========================================================================

class TestSmokeTestCli:
    def test_flag_defaults_to_off(self):
        assert run_pipeline.parse_args([]).smoke_test is False

    def test_flag_parses(self):
        assert run_pipeline.parse_args(["--smoke-test"]).smoke_test is True

    def test_simulation_script_accepts_the_flag(self):
        """`run_simulation.py --smoke-test` must parse without side effects."""
        import ast
        source = open(
            os.path.join(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))), "run_simulation.py")
        ).read()
        tree = ast.parse(source)
        assert "--smoke-test" in source
        # The override block must be guarded, never unconditional.
        assert "if SMOKE_TEST:" in source
        assert isinstance(tree, ast.Module)

    def test_production_defaults_are_not_edited_by_the_override_path(self):
        """
        The documented values must appear verbatim in the config block.

        This is the regression guard for the constraint that a smoke test may
        never rewrite a production parameter in place.
        """
        source = open(
            os.path.join(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__))), "run_simulation.py")
        ).read()
        config = source.split("End of configuration")[0]
        assert "SURVEY_AREA_DEG2 = 10.0" in config
        assert "target_cell_size_mpc=2.0" in config
        assert "integration_time   = 1000 * 3600" in config
        assert "photoz_uncertainty = 0.45" in config
        assert "wedge_buffer = 0.0677" in config
        assert "n_bins_perp     = 20" in config
        assert "n_bins_parallel = 20" in config
        assert "RANDOM_SEED = 42" in config
        # And no smoke value leaked into it.
        for entry in smoke_test.SMOKE_TEST_OVERRIDES.values():
            assert f"= {entry['value']}\n" not in config


def test_smoke_test_run_verifies_every_stage(tmp_path):
    """
    End-to-end wiring check on the synthetic fixture.

    Runs the driver with ``--smoke-test`` against the tiny test simulation and
    ``--sim skip``, so 21cmFAST is never invoked: this exercises the report
    plumbing, not the reduced simulation configuration.
    """
    from conftest import write_tiny_simulation

    workspace = tmp_path / "smoke"
    workspace.mkdir()
    data_path = write_tiny_simulation(str(workspace / "lightcone_data.h5"))

    exit_code = run_pipeline.main([
        "--data", data_path,
        "--products", str(workspace / "analysis_products.h5"),
        "--figdir", str(workspace / "figures"),
        "--summary", str(workspace / "pipeline_summary.json"),
        "--sim", "skip",
        "--plots", "none",
        "--smoke-test",
        "--quiet",
    ])
    assert exit_code == 0

    summary = json.loads((workspace / "pipeline_summary.json").read_text())
    assert summary["estimator"] == "coeval"
