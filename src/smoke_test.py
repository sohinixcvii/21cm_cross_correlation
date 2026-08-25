#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
smoke_test.py — pre-flight overrides and shape checks
======================================================

A **smoke test** runs the whole pipeline end to end on a deliberately tiny
configuration, to prove that every stage executes and produces
correctly-shaped output before an HPC job is submitted. It is not a science
run and its numbers mean nothing.

.. warning::

   **Nothing in this module is a scientific parameter.** The values in
   :data:`SMOKE_TEST_OVERRIDES` are chosen to make the pipeline finish in
   seconds, and are wrong by construction: a 32 Mpc box does not sample the
   modes the forecast needs, 16 cells per side does not resolve the ionisation
   field, and a 12-slice line of sight has no usable ``k_parallel`` axis.
   Never copy a number from this file into ``run_simulation.py``'s
   configuration block, into the documentation, or into a paper.

How it stays out of the production path
---------------------------------------
This module is imported **only** when ``--smoke-test`` is passed. The
production defaults in ``run_simulation.py`` are not edited, not shadowed and
not reassigned unless that flag is set, and a smoke run writes to
``outputs/smoke_test/`` so it cannot overwrite a real run's products.

What it verifies
----------------
Not merely "it did not crash". Every stage's output is checked for shape,
dtype and finiteness by the ``check_*`` functions below, which are also
exercised by the ordinary test suite (``tests/test_smoke_test.py``) so the
checker itself cannot rot.

Stages covered, in pipeline order::

    initial conditions -> halo catalogue -> perturbed catalogue
      -> 21 cm lightcone -> galaxy overdensity -> galaxy bias -> Kaiser RSD
      -> HDF5 -> cylindrical power spectra -> photo-z / wedge / noise / SNR
      -> figures -> summary JSON

There is no MCMC stage in this project. Nothing here samples a posterior,
``emcee`` is not a dependency, and no chain is produced — so no chain shape is
asserted. :func:`check_mcmc_chain` exists as the hook to fill in when a
sampler is added; it currently reports the stage as absent rather than
passing silently.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "SMOKE_TEST_OVERRIDES",
    "SMOKE_OUTPUT_DIR",
    "describe_overrides",
    "SmokeCheck",
    "SmokeReport",
    "check_simulation_output",
    "check_power_spectra",
    "check_uncertainty_budget",
    "check_summary",
    "check_figures",
    "check_mcmc_chain",
]


#: Where a smoke run writes, so it can never overwrite a real run's outputs.
SMOKE_OUTPUT_DIR = os.path.join("outputs", "smoke_test")


#: Reduced parameters for a smoke run.  **Not science values** — see the module
#: warning.  Each entry names the production value it stands in for, so the
#: reduction is auditable and nobody has to guess what was changed.
#:
#: Chosen against the runtime drivers measured in ``docs/simulation_spec.md``:
#: the halo catalogue scales with comoving volume and the initial conditions
#: with ``DIM^3``, so shrinking the box is worth far more than shrinking
#: anything else.  ``(32 / 486.33)^3`` is a factor 3500 in catalogue size.
SMOKE_TEST_OVERRIDES: Dict[str, Dict[str, Any]] = {
    "HII_DIM": {
        "value": 16,
        "production": 256,
        "why": "ionisation/21 cm grid; cost scales as HII_DIM^3",
    },
    "BOX_LEN": {
        "value": 32.0,
        "production": 486.33,
        "why": "comoving box [Mpc]; the halo catalogue scales with its cube",
    },
    "DIM": {
        "value": 48,
        "production": 768,
        "why": "initial-conditions grid, held at the 3x HII_DIM convention",
    },
    "minimum_los_slices": {
        "value": 12,
        "production": 100,
        "why": "line-of-sight slices; sets the lightcone's third dimension",
    },
    "n_bins_perp": {
        "value": 8,
        "production": 20,
        "why": "a 16-cell transverse grid cannot fill 20 log-spaced k bins",
    },
    "n_bins_parallel": {
        "value": 8,
        "production": 20,
        "why": "as above, for the line-of-sight axis",
    },
    "max_halos": {
        "value": 200_000,
        "production": 0,
        "why": "cap on the catalogue the analysis stage loads (0 = all)",
    },
}

#: Parameters deliberately **not** overridden, with the reason.  Recorded so a
#: future reader does not assume they were forgotten.
SMOKE_TEST_UNCHANGED: Dict[str, str] = {
    "SURVEY_AREA_DEG2": (
        "the box is overridden directly instead, so the documented 10 deg^2 "
        "footprint stays exactly as written"
    ),
    "integration_time": (
        "a post-processing scalar in P_N = T_sys^2 x 1e3 / (t B); it changes "
        "no array size and costs no runtime"
    ),
    "bandwidth": "as above, and it sets the sub-band split the smoke run should exercise",
    "photoz_uncertainty": "post-processing scalar; no runtime cost",
    "wedge_buffer": "post-processing scalar; no runtime cost",
    "M_UV_limit / M_UV_bright": "selection thresholds; no runtime cost",
    "RANDOM_SEED": "kept at 42 so a smoke run is reproducible",
    "SAMPLER_MIN_MASS": (
        "left at the template default; shrinking the box already removes the "
        "catalogue cost, and raising the mass floor would change which halos "
        "exist rather than how many boxes' worth of them do"
    ),
}


def describe_overrides() -> str:
    """
    Render the override table for the run log.

    Returns
    -------
    str
        A block naming every overridden parameter, its smoke value, the
        production value it replaces, and why it was reduced.
    """
    lines = [
        "  ┌─ SMOKE TEST — reduced configuration, NOT a science run ─────────",
        "  │  parameter            smoke        production   why",
    ]
    for name, entry in SMOKE_TEST_OVERRIDES.items():
        lines.append(
            f"  │  {name:<20} {str(entry['value']):<12} "
            f"{str(entry['production']):<12} {entry['why']}"
        )
    lines.append("  │  unchanged: " + ", ".join(SMOKE_TEST_UNCHANGED))
    lines.append("  └──────────────────────────────────────────────────────────────────")
    return "\n".join(lines)


def override(name: str, production_value: Any) -> Any:
    """
    Return the smoke value for ``name``, or the production value if not reduced.

    Parameters
    ----------
    name : str
        Parameter name, as keyed in :data:`SMOKE_TEST_OVERRIDES`.
    production_value : Any
        The value the production path would use.  Returned unchanged when the
        parameter is not one of the overrides, so callers never have to
        special-case a missing key.

    Returns
    -------
    Any
        The smoke value, or ``production_value``.
    """
    entry = SMOKE_TEST_OVERRIDES.get(name)
    return production_value if entry is None else entry["value"]


# ===========================================================================
#  Report
# ===========================================================================

@dataclass
class SmokeCheck:
    """
    One stage's verdict.

    Attributes
    ----------
    stage : str
        Pipeline stage, e.g. ``"power spectra"``.
    passed : bool
        Whether every assertion for that stage held.
    detail : str
        What was checked, or what failed.
    """

    stage: str
    passed: bool
    detail: str


@dataclass
class SmokeReport:
    """
    Collected verdicts for a smoke run.

    Attributes
    ----------
    checks : list of SmokeCheck
        One entry per stage, in pipeline order.
    """

    checks: List[SmokeCheck] = field(default_factory=list)

    def add(self, stage: str, passed: bool, detail: str) -> None:
        """Record one stage's verdict."""
        self.checks.append(SmokeCheck(stage=stage, passed=passed, detail=detail))

    @property
    def passed(self) -> bool:
        """True when every recorded check passed."""
        return all(check.passed for check in self.checks)

    def render(self) -> str:
        """
        Format the report for the terminal.

        Returns
        -------
        str
            One line per stage, then a verdict line.
        """
        width = max((len(c.stage) for c in self.checks), default=10)
        lines = ["  SMOKE TEST — stage checks", "  " + "-" * (width + 40)]
        for check in self.checks:
            mark = "PASS" if check.passed else "FAIL"
            lines.append(f"  [{mark}] {check.stage:<{width}}  {check.detail}")
        lines.append("  " + "-" * (width + 40))
        lines.append(
            "  RESULT: all stages executed and produced correctly-shaped output"
            if self.passed
            else "  RESULT: FAILED — see the lines marked FAIL above"
        )
        lines.append(
            "  (shapes only; the numbers from a smoke run are not physical)"
        )
        return "\n".join(lines)


# ===========================================================================
#  Stage checks
# ===========================================================================

def _shape_of(array: Any) -> Tuple[int, ...]:
    """Shape of an array-like, or ``()`` when it has none."""
    return tuple(getattr(array, "shape", ()))


def check_simulation_output(
    data: Any,
    expected_hii_dim: int,
    expected_n_z: int,
    report: Optional[SmokeReport] = None,
) -> SmokeReport:
    """
    Verify the simulation HDF5's fields, catalogue and geometry.

    Parameters
    ----------
    data : SimulationData
        Loaded simulation, from :func:`src.dataio.load_simulation`.
    expected_hii_dim, expected_n_z : int
        Grid dimensions the smoke configuration asked for.
    report : SmokeReport, optional
        Report to append to; a new one is created when omitted.

    Returns
    -------
    SmokeReport
        The report, with one entry per checked stage.
    """
    report = report or SmokeReport()
    expected = (expected_hii_dim, expected_hii_dim, expected_n_z)

    fields = {
        "brightness_temp_field": data.brightness_temp_field,
        "density_field": data.density_field,
        "neutral_fraction": data.neutral_fraction,
        "galaxy_overdensity": data.galaxy_overdensity,
    }
    problems = [
        f"{name} has shape {_shape_of(array)}, expected {expected}"
        for name, array in fields.items()
        if _shape_of(array) != expected
    ]
    problems += [
        f"{name} is not finite everywhere"
        for name, array in fields.items()
        if _shape_of(array) == expected and not np.all(np.isfinite(array))
    ]
    report.add(
        "21 cm lightcone",
        not problems,
        "; ".join(problems) or f"4 fields at {expected}, all finite",
    )

    neutral = np.asarray(data.neutral_fraction)
    in_range = bool(neutral.min() >= -1e-6 and neutral.max() <= 1 + 1e-6)
    report.add(
        "ionisation field",
        in_range,
        f"x_HI in [{neutral.min():.3f}, {neutral.max():.3f}]"
        + ("" if in_range else " — outside [0, 1]"),
    )

    los = {
        "lc_redshifts": data.lc_redshifts,
        "lc_dist_Mpc": data.lc_dist_Mpc,
    }
    los_problems = [
        f"{name} has shape {_shape_of(array)}, expected ({expected_n_z},)"
        for name, array in los.items()
        if _shape_of(array) != (expected_n_z,)
    ]
    report.add(
        "lightcone geometry",
        not los_problems,
        "; ".join(los_problems)
        or f"{expected_n_z} slices, z = {np.min(data.lc_redshifts):.4f}"
           f"–{np.max(data.lc_redshifts):.4f}",
    )

    n_loaded = int(_shape_of(data.halo_masses)[0]) if data.halo_masses.size else 0
    catalogue_problems = []
    if n_loaded == 0 and int(data.n_halos_total) > 0:
        # The driver skips the catalogue read when no requested figure needs
        # it.  That is deliberate, and different from a catalogue that is
        # missing or malformed, so it is reported rather than failed.
        report.add(
            "halo catalogue",
            True,
            f"{int(data.n_halos_total):,} halos in the file, not loaded "
            f"(no stage requested them)",
        )
        return report
    if n_loaded == 0:
        catalogue_problems.append("no halo catalogue in the file")
    else:
        if _shape_of(data.halo_coords) != (n_loaded, 3):
            catalogue_problems.append(
                f"halo_coords is {_shape_of(data.halo_coords)}, "
                f"expected ({n_loaded}, 3)"
            )
        for name in ("stellar_masses", "sfr"):
            if _shape_of(getattr(data, name)) != (n_loaded,):
                catalogue_problems.append(
                    f"{name} is {_shape_of(getattr(data, name))}, "
                    f"expected ({n_loaded},)"
                )
        if np.any(np.asarray(data.halo_masses) <= 0):
            catalogue_problems.append("non-positive halo masses")
    report.add(
        "halo catalogue",
        not catalogue_problems,
        "; ".join(catalogue_problems)
        or f"{n_loaded:,} halos loaded of {data.n_halos_total:,}, "
           f"4 arrays consistent",
    )
    return report


def check_power_spectra(
    spectra: Any,
    expected_n_perp: int,
    expected_n_parallel: int,
    report: Optional[SmokeReport] = None,
    subbands: Any = None,
) -> SmokeReport:
    """
    Verify the cylindrical power spectra, and the sub-bands when present.

    Parameters
    ----------
    spectra : PowerSpectra
        Spectra from the analysis stage.
    expected_n_perp, expected_n_parallel : int
        Binning the smoke configuration asked for.
    report : SmokeReport, optional
        Report to append to.
    subbands : SubbandPowerSpectra, optional
        Per-band spectra, when the lightcone estimator ran.

    Returns
    -------
    SmokeReport
        The report.
    """
    report = report or SmokeReport()
    expected = (expected_n_perp, expected_n_parallel)

    problems = []
    for name in ("P_21cm_auto", "P_galaxy_auto", "P_cross", "mode_counts"):
        if _shape_of(getattr(spectra, name)) != expected:
            problems.append(
                f"{name} is {_shape_of(getattr(spectra, name))}, "
                f"expected {expected}"
            )
    if _shape_of(spectra.k_perp) != (expected_n_perp,):
        problems.append(f"k_perp is {_shape_of(spectra.k_perp)}")
    if _shape_of(spectra.k_parallel) != (expected_n_parallel,):
        problems.append(f"k_parallel is {_shape_of(spectra.k_parallel)}")

    populated = int(np.sum(np.asarray(spectra.mode_counts) > 0))
    if populated == 0:
        problems.append("every bin is empty")

    report.add(
        "power spectra",
        not problems,
        "; ".join(problems)
        or f"3 spectra at {expected}, {populated}/{expected[0] * expected[1]} "
           f"bins populated",
    )

    if subbands is not None:
        band_problems = []
        if subbands.n_bands < 1:
            band_problems.append("no sub-bands produced")
        for index, band in enumerate(subbands.bands):
            if _shape_of(band.P_cross) != expected:
                band_problems.append(
                    f"band {index} P_cross is {_shape_of(band.P_cross)}"
                )
        for name in ("z_effective", "bandwidth_hz", "los_length_mpc", "n_slices"):
            if _shape_of(getattr(subbands, name)) != (subbands.n_bands,):
                band_problems.append(f"{name} has the wrong length")
        report.add(
            "sub-bands",
            not band_problems,
            "; ".join(band_problems)
            or f"{subbands.n_bands} bands, each {expected}, "
               f"z_eff {np.min(subbands.z_effective):.3f}"
               f"–{np.max(subbands.z_effective):.3f}",
        )
    return report


def check_uncertainty_budget(
    budget: Any,
    expected_n_perp: int,
    expected_n_parallel: int,
    report: Optional[SmokeReport] = None,
) -> SmokeReport:
    """
    Verify the photo-z, wedge, noise and SNR maps.

    Parameters
    ----------
    budget : UncertaintyBudget
        Output of the observational stage.
    expected_n_perp, expected_n_parallel : int
        Binning the smoke configuration asked for.
    report : SmokeReport, optional
        Report to append to.

    Returns
    -------
    SmokeReport
        The report.
    """
    report = report or SmokeReport()
    expected = (expected_n_perp, expected_n_parallel)

    problems = []
    if _shape_of(budget.outside_wedge) != expected:
        problems.append(f"wedge mask is {_shape_of(budget.outside_wedge)}")
    if _shape_of(budget.P_cross_observed) != expected:
        problems.append(f"damped P_cross is {_shape_of(budget.P_cross_observed)}")
    if not np.isfinite(budget.radial_smearing):
        problems.append("radial smearing is not finite")
    if budget.horizon_slope <= 0:
        problems.append(f"horizon slope is {budget.horizon_slope}")
    total = float(budget.total_snr)
    if not np.isfinite(total) or total < 0:
        problems.append(f"total SNR is {total}")

    report.add(
        "uncertainty budget",
        not problems,
        "; ".join(problems)
        or f"maps at {expected}, sigma_r = {budget.radial_smearing:.1f} Mpc, "
           f"{int(np.sum(budget.outside_wedge))} bins outside the wedge, "
           f"SNR = {total:.3g}",
    )
    return report


def check_summary(
    summary: Dict[str, Any],
    report: Optional[SmokeReport] = None,
) -> SmokeReport:
    """
    Verify the summary JSON has every block the pipeline promises.

    Parameters
    ----------
    summary : dict
        Parsed ``pipeline_summary.json``.
    report : SmokeReport, optional
        Report to append to.

    Returns
    -------
    SmokeReport
        The report.
    """
    report = report or SmokeReport()
    required = ("generated", "data_file", "simulation", "power_spectra",
                "uncertainty_budget", "estimator", "figures")
    missing = [key for key in required if key not in summary]
    report.add(
        "summary JSON",
        not missing,
        f"missing keys: {missing}" if missing
        else f"{len(summary)} top-level keys, estimator "
             f"'{summary.get('estimator')}'",
    )
    return report


def check_figures(
    figure_paths: Sequence[str],
    report: Optional[SmokeReport] = None,
) -> SmokeReport:
    """
    Verify the figure files exist and are non-empty.

    Parameters
    ----------
    figure_paths : sequence of str
        Paths the figure stage reported writing.
    report : SmokeReport, optional
        Report to append to.

    Returns
    -------
    SmokeReport
        The report.  A run with ``--plots none`` records the stage as skipped
        rather than failed.
    """
    report = report or SmokeReport()
    if not figure_paths:
        report.add("figures", True, "skipped (--plots none)")
        return report

    problems = [
        f"{os.path.basename(path)} is missing or empty"
        for path in figure_paths
        if not os.path.exists(path) or os.path.getsize(path) == 0
    ]
    report.add(
        "figures",
        not problems,
        "; ".join(problems) or f"{len(figure_paths)} files written, all non-empty",
    )
    return report


def check_mcmc_chain(
    chain: Any = None,
    expected_walkers: Optional[int] = None,
    expected_steps: Optional[int] = None,
    report: Optional[SmokeReport] = None,
) -> SmokeReport:
    """
    Verify a posterior sampler's chain — **not yet applicable to this project**.

    This pipeline has no MCMC stage: nothing samples a posterior, ``emcee`` is
    not a dependency, and no chain is written. The hook exists so that adding a
    sampler later is a matter of passing its chain here rather than inventing a
    checker under time pressure, and so that the smoke report says *absent*
    rather than passing over a stage that does not exist.

    Parameters
    ----------
    chain : ndarray, optional
        Sampler output of shape ``(n_walkers, n_steps, n_parameters)``.
    expected_walkers, expected_steps : int, optional
        Dimensions the smoke configuration asked for.
    report : SmokeReport, optional
        Report to append to.

    Returns
    -------
    SmokeReport
        The report.
    """
    report = report or SmokeReport()
    if chain is None:
        report.add(
            "MCMC",
            True,
            "no sampler in this pipeline — stage absent, nothing to check",
        )
        return report

    array = np.asarray(chain)
    problems = []
    if array.ndim != 3:
        problems.append(f"chain is {array.ndim}D, expected 3D "
                        f"(walkers, steps, parameters)")
    else:
        walkers, steps, _ = array.shape
        if expected_walkers is not None and walkers != expected_walkers:
            problems.append(f"{walkers} walkers, expected {expected_walkers}")
        if expected_steps is not None and steps != expected_steps:
            problems.append(f"{steps} steps, expected {expected_steps}")
        if not np.all(np.isfinite(array)):
            problems.append("chain contains non-finite samples")

    report.add(
        "MCMC",
        not problems,
        "; ".join(problems) or f"chain {array.shape}, finite",
    )
    return report
