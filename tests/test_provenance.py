#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/provenance.py``.

The manifest's contract is that it survives a run that never gets to close
it, so most of these check what is on disk *before* ``finish()`` is called.
"""

from __future__ import annotations

import json
import os

import pytest

from src import provenance
from src.provenance import (
    INT32_MAX,
    SAMPLER_MIN_MASS_REFERENCE,
    SAMPLER_RETAINED_FRACTION,
    RunManifest,
    environment_info,
    estimate_cache_footprint,
    estimate_catalogue_cost,
    git_revision,
    package_versions,
    peak_memory_gb,
    resolve_n_threads,
    sampler_retained_fraction,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ===========================================================================
#  Thread resolution
# ===========================================================================

def test_n_threads_prefers_explicit_variable(monkeypatch) -> None:
    """``N_THREADS`` wins over the SLURM allocation."""
    monkeypatch.setenv("N_THREADS", "12")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    assert resolve_n_threads() == 12


def test_n_threads_falls_back_to_slurm(monkeypatch) -> None:
    """SLURM's allocation is used when nothing more specific is set."""
    monkeypatch.delenv("N_THREADS", raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    assert resolve_n_threads() == 8


def test_n_threads_uses_default_then_cpu_count(monkeypatch) -> None:
    """With no environment, the explicit default wins over ``cpu_count()``."""
    monkeypatch.delenv("N_THREADS", raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    assert resolve_n_threads(default=3) == 3
    assert resolve_n_threads() == max(os.cpu_count() or 1, 1)


@pytest.mark.parametrize("value", ["", "all", "8x", "-4"])
def test_n_threads_survives_malformed_values(monkeypatch, value: str) -> None:
    """
    A bad environment value must never abort a queued job.

    Empty and non-numeric values fall through; a negative one is clamped to
    the 1-thread floor rather than being passed to 21cmFAST.
    """
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setenv("N_THREADS", value)
    assert resolve_n_threads(default=2) >= 1


# ===========================================================================
#  Environment capture
# ===========================================================================

def test_git_revision_reads_this_repository() -> None:
    """The repo's own revision is discoverable."""
    revision = git_revision(REPO_ROOT)
    assert revision["commit"] is not None
    assert len(revision["commit"]) == 40
    assert isinstance(revision["dirty"], bool)


def test_git_revision_outside_a_repository(tmp_path) -> None:
    """A non-repository yields Nones, not an exception."""
    revision = git_revision(str(tmp_path))
    assert revision == {"commit": None, "branch": None, "dirty": None}


def test_package_versions_reports_numpy() -> None:
    """Importable packages get a version; missing ones get ``None``."""
    versions = package_versions()
    assert versions["numpy"] is not None
    assert "py21cmfast" in versions


def test_environment_info_is_json_serialisable() -> None:
    """Everything captured must survive a round trip through JSON."""
    info = environment_info(REPO_ROOT)
    assert info["host"] and info["python"]
    assert json.loads(json.dumps(info, default=str))["git"]["commit"] is not None


def test_peak_memory_is_positive() -> None:
    """Peak RSS is reported in GB, and this process has used some."""
    peak = peak_memory_gb()
    assert peak is not None and 0.0 < peak < 1000.0


# ===========================================================================
#  Cost estimate
# ===========================================================================

def test_cost_estimate_reproduces_the_reference_run() -> None:
    """The 256 Mpc box must return both counts of the run it was calibrated on."""
    cost = estimate_catalogue_cost(256.0)
    assert cost["n_halos_lagrangian"] == pytest.approx(136_663_818, rel=1e-9)
    assert cost["n_halos_perturbed"] == pytest.approx(114_289_081, rel=1e-9)
    assert cost["catalogue_GB"] == pytest.approx(3.83, rel=0.01)
    assert cost["int32_headroom"] < 0.2


def test_cost_estimate_matches_an_independent_run() -> None:
    """
    A 64 Mpc 21cmFAST run on 2026-08-21 produced 1,782,540 perturbed halos.

    That run is not part of the calibration, so it checks the volume scaling
    and the Lagrangian-to-perturbed ratio against something held out.
    """
    cost = estimate_catalogue_cost(64.0)
    assert cost["n_halos_perturbed"] == pytest.approx(1_782_540, rel=0.02)


def test_cost_estimate_scales_with_volume() -> None:
    """Doubling the box is 8x the halos — the sampler floor is a mass, not a cell."""
    small = estimate_catalogue_cost(100.0)
    large = estimate_catalogue_cost(200.0)
    assert large["n_halos_lagrangian"] == pytest.approx(8 * small["n_halos_lagrangian"])
    assert large["resident_GB"] == pytest.approx(8 * small["resident_GB"])


def test_cost_estimate_flags_the_box_that_crashed() -> None:
    """
    The 2026-08-20 geometry must trip the int32 guard.

    BOX_LEN = 486.33 Mpc draws ~9.4e8 halos, whose flattened ``halo_coords``
    is 1.3x INT_MAX — the run that segfaulted after 38 minutes.
    """
    cost = estimate_catalogue_cost(486.32904710223784)
    assert cost["n_halos_lagrangian"] == pytest.approx(9.37e8, rel=0.01)
    assert cost["catalogue_GB"] == pytest.approx(26.2, rel=0.02)
    assert cost["resident_GB"] == pytest.approx(48.2, rel=0.02)
    assert cost["int32_headroom"] > 1.0
    assert cost["n_halos_lagrangian"] * 3 > INT32_MAX


def test_cost_estimate_clears_a_smaller_box() -> None:
    """350 Mpc — the suggested intermediate — stays inside the index range."""
    cost = estimate_catalogue_cost(350.0)
    assert cost["int32_headroom"] < 1.0
    assert cost["catalogue_GB"] < 12.0


# ===========================================================================
#  The manifest
# ===========================================================================

def test_manifest_is_on_disk_before_anything_else_happens(tmp_path) -> None:
    """
    ``create()`` writes immediately.

    This is the whole point: a run killed by a signal cannot flush stdout or
    run an exit hook, so the parameters have to already be on disk.
    """
    manifest = RunManifest.create(str(tmp_path / "runs"), label="sim")

    assert os.path.exists(manifest.path)
    stored = json.load(open(manifest.path))
    assert stored["status"] == "running"
    assert stored["environment"]["host"]


def test_manifest_records_stage_it_died_in(tmp_path) -> None:
    """An unclosed stage is left named on disk, with the earlier ones logged."""
    manifest = RunManifest.create(str(tmp_path / "runs"))
    manifest.record("parameters", {"BOX_LEN": 486.33, "HII_DIM": 256})
    manifest.begin_stage("lightcone")
    manifest.end_stage()
    manifest.begin_stage("halo_catalog")     # never closed — the crash

    stored = json.load(open(manifest.path))
    assert stored["status"] == "running"
    assert stored["stage"] == "halo_catalog"
    assert stored["stages_completed"] == ["lightcone"]
    assert "lightcone" in stored["timings_seconds"]
    assert stored["parameters"]["BOX_LEN"] == 486.33


def test_manifest_record_merges_sections(tmp_path) -> None:
    """Repeated ``record`` calls accumulate rather than replacing."""
    manifest = RunManifest.create(str(tmp_path / "runs"))
    manifest.record("results", {"n_halos": 10})
    manifest.record("results", {"mean_neutral_fraction": 0.17})

    stored = json.load(open(manifest.path))
    assert stored["results"] == {"n_halos": 10, "mean_neutral_fraction": 0.17}


def test_manifest_finish_marks_status(tmp_path) -> None:
    """A closed manifest carries its outcome and a finish timestamp."""
    manifest = RunManifest.create(str(tmp_path / "runs"))
    manifest.finish("complete")

    stored = json.load(open(manifest.path))
    assert stored["status"] == "complete"
    assert stored["finished"] is not None
    assert stored["elapsed_seconds"] >= 0.0


def test_manifest_end_stage_without_begin_is_harmless(tmp_path) -> None:
    """Closing a stage that was never opened returns 0, not an error."""
    manifest = RunManifest.create(str(tmp_path / "runs"))
    assert manifest.end_stage() == 0.0


def test_manifest_filename_carries_the_run_id(tmp_path) -> None:
    """The run id is both the filename and a field, so files sort by time."""
    manifest = RunManifest.create(str(tmp_path / "runs"), label="sim",
                                  run_id="20260821_120000")
    assert manifest.path.endswith("sim_20260821_120000.json")
    assert json.load(open(manifest.path))["run_id"] == "20260821_120000"


def test_manifest_write_is_atomic(tmp_path) -> None:
    """No temporary file is left behind, and the JSON always parses."""
    manifest = RunManifest.create(str(tmp_path / "runs"))
    for index in range(5):
        manifest.record("results", {"step": index})
        json.load(open(manifest.path))          # must parse at every point

    assert not os.path.exists(f"{manifest.path}.tmp")


def test_manifest_survives_unserialisable_values(tmp_path) -> None:
    """
    ``default=str`` keeps provenance capture from ever failing a run.

    Numpy scalars and arrays turn up in the recorded results; none of them
    should be able to raise out of a ``record`` call.
    """
    import numpy as np

    manifest = RunManifest.create(str(tmp_path / "runs"))
    manifest.record("results", {"shape": np.array([1, 2, 3]), "x": np.float32(1.5)})

    assert json.load(open(manifest.path))["results"]["x"]


def test_calibration_constants_are_self_consistent() -> None:
    """The published calibration must match the measurement it came from."""
    assert provenance.HALOS_PER_MPC3 == pytest.approx(8.146, rel=1e-3)
    assert provenance.BYTES_PER_HALO == pytest.approx(28.0, rel=1e-3)
    assert provenance.PERTURBED_FRACTION == pytest.approx(0.836, rel=1e-3)


# ---------------------------------------------------------------------------
#  SAMPLER_MIN_MASS and the 32-bit halo index
# ---------------------------------------------------------------------------


def test_retained_fraction_is_one_at_and_below_the_reference_mass() -> None:
    """The calibration was measured at 1e8; it cannot extrapolate downward."""
    assert sampler_retained_fraction(SAMPLER_MIN_MASS_REFERENCE) == 1.0
    assert sampler_retained_fraction(1e7) == 1.0


@pytest.mark.parametrize("mass", sorted(SAMPLER_RETAINED_FRACTION))
def test_retained_fraction_reproduces_its_tabulated_points(mass) -> None:
    """Interpolation must pass through every tabulated value."""
    assert sampler_retained_fraction(mass) == pytest.approx(
        SAMPLER_RETAINED_FRACTION[mass], rel=1e-9
    )


def test_retained_fraction_decreases_with_the_floor() -> None:
    """A higher mass floor can only keep fewer halos."""
    masses = [1e8, 1.25e8, 1.5e8, 2e8, 2.5e8, 3e8]
    fractions = [sampler_retained_fraction(m) for m in masses]
    assert all(b < a for a, b in zip(fractions, fractions[1:]))


def test_retained_fraction_rejects_a_non_positive_floor() -> None:
    """A zero or negative mass floor is a configuration error, not a default."""
    for bad in (0.0, -1e8):
        with pytest.raises(ValueError, match="sampler_min_mass"):
            sampler_retained_fraction(bad)


def test_cost_estimate_default_reproduces_the_reference_calibration() -> None:
    """Omitting the floor must behave exactly as the one-argument form did."""
    assert estimate_catalogue_cost(486.33) == estimate_catalogue_cost(
        486.33, SAMPLER_MIN_MASS_REFERENCE
    )


def test_raising_the_sampler_floor_clears_the_int32_guard() -> None:
    """
    The adopted production setting must fit inside a signed 32-bit index.

    BOX_LEN = 486.33 Mpc overflows at the template's 1e8 floor (1.31x
    INT_MAX, the run that segfaulted).  SAMPLER_MIN_MASS = 2e8 is the adopted
    fix: it keeps the footprint-derived box and brings the flattened
    ``halo_coords`` to 0.61x INT_MAX.
    """
    overflowing = estimate_catalogue_cost(486.33, 1e8)
    adopted = estimate_catalogue_cost(486.33, 2e8)

    assert overflowing["int32_headroom"] > 1.0
    assert adopted["int32_headroom"] < 1.0
    assert adopted["int32_headroom"] == pytest.approx(0.605, abs=0.01)
    assert adopted["n_halos_lagrangian"] * 3 < INT32_MAX


def test_cost_estimate_reports_the_floor_it_used() -> None:
    """The floor and its retained fraction travel with the estimate."""
    cost = estimate_catalogue_cost(486.33, 2e8)
    assert cost["sampler_min_mass"] == 2e8
    assert cost["sampler_retained_fraction"] == pytest.approx(0.462, rel=1e-9)


def test_the_sampler_floor_does_not_change_the_volume() -> None:
    """The floor scales the catalogue, not the box it is drawn in."""
    for mass in (1e8, 2e8, 3e8):
        assert estimate_catalogue_cost(486.33, mass)["volume_Mpc3"] == (
            pytest.approx(486.33 ** 3)
        )


# ---------------------------------------------------------------------------
#  Cache footprint — the disk cost that DOES scale with the redshift range
# ---------------------------------------------------------------------------


def test_cache_scales_linearly_with_node_count() -> None:
    """One halo catalogue per node: twice the nodes, twice the catalogues."""
    five = estimate_cache_footprint(486.33, 5, 2e8)
    ten = estimate_cache_footprint(486.33, 10, 2e8)
    assert ten["total_GB"] - ten["ics_GB"] == pytest.approx(
        2 * (five["total_GB"] - five["ics_GB"])
    )
    assert ten["ics_GB"] == pytest.approx(five["ics_GB"])


def test_cache_per_node_is_one_catalogue() -> None:
    """The measured 3.83 GB/node at 256 Mpc is exactly one catalogue."""
    assert estimate_cache_footprint(256.0, 1, 1e8)["per_node_GB"] == (
        pytest.approx(estimate_catalogue_cost(256.0, 1e8)["catalogue_GB"])
    )
    assert estimate_cache_footprint(256.0, 1, 1e8)["per_node_GB"] == (
        pytest.approx(3.83, abs=0.01)
    )


def test_cache_reproduces_the_run_that_filled_the_disk() -> None:
    """
    z = 6.5-7.5 at BOX_LEN = 486.33 needs ~128 GB, and errno 28'd overnight.

    Ten node redshifts, each caching a 12.1 GB halo catalogue at the adopted
    SAMPLER_MIN_MASS = 2e8.
    """
    cache = estimate_cache_footprint(486.33, 10, 2e8)
    assert cache["total_GB"] == pytest.approx(128, abs=2)
    assert cache["total_upper_GB"] == pytest.approx(229, abs=3)
    assert cache["total_upper_GB"] > cache["total_GB"]


def test_raising_the_sampler_floor_shrinks_the_cache() -> None:
    """The floor is the cheapest lever on cache size after the redshift span."""
    assert (
        estimate_cache_footprint(486.33, 10, 3e8)["total_GB"]
        < estimate_cache_footprint(486.33, 10, 2e8)["total_GB"]
        < estimate_cache_footprint(486.33, 10, 1e8)["total_GB"]
    )


def test_cache_rejects_a_non_positive_node_count() -> None:
    """A lightcone with no nodes is a configuration error."""
    for bad in (0, -3):
        with pytest.raises(ValueError, match="n_nodes"):
            estimate_cache_footprint(486.33, bad)
