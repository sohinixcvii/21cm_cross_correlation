#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
provenance.py — per-run parameter manifests
============================================

Every simulation run writes one JSON manifest next to its outputs, recording
the parameters it used, the code and environment it ran under, and how far it
got.  The manifest is rewritten on **every** update rather than once at the
end, so a run that dies mid-stage — including one killed by a signal, which
Python cannot trap — still leaves a complete record of its configuration and
a ``status`` of ``"running"`` naming the stage it was in.

That last property is the point.  A crashed run's stdout is lost if it was
block-buffered; the manifest is not.

Usage
-----
::

    manifest = RunManifest.create("outputs/runs", label="sim")
    manifest.record("parameters", {"BOX_LEN": 486.33, "HII_DIM": 256})
    manifest.begin_stage("lightcone")
    ...
    manifest.end_stage()
    manifest.finish("complete")
"""

from __future__ import annotations

import json
import os
import platform
import socket
import subprocess
import sys
import time
from typing import Any, Dict, Mapping, Optional

__all__ = [
    "HALOS_PER_MPC3",
    "BYTES_PER_HALO",
    "PERTURBED_FRACTION",
    "INT32_MAX",
    "RunManifest",
    "estimate_catalogue_cost",
    "environment_info",
    "git_revision",
    "resolve_n_threads",
    "package_versions",
    "peak_memory_gb",
]

# Empirical calibration from the 12 Aug 2026 run and its 21cmFAST cache:
# 136,663,818 halos in a (256 Mpc)^3 box at z = 7 with SAMPLER_MIN_MASS = 1e8,
# stored in a 3.564 GiB HaloCatalog.h5.  Used only for the pre-flight cost
# estimate below — nothing downstream consumes it.
HALOS_PER_MPC3 = 136_663_818 / 256.0 ** 3      # 8.146 halos Mpc^-3
BYTES_PER_HALO = 3.564 * 2 ** 30 / 136_663_818  # 28.0 bytes per halo

#: Fraction of the Lagrangian catalogue that survives ``perturb_halo_catalog``.
#: 114,289,081 / 136,663,818 from the 12 Aug 256 Mpc run; independently
#: reproduced as 1,782,540 / 2,135,000 by a 64 Mpc run on 21 Aug.  The
#: Lagrangian count is the one that sets peak memory and index width, so the
#: estimates below lead with it.
PERTURBED_FRACTION = 114_289_081 / 136_663_818  # 0.836

#: Largest signed 32-bit index.  21cmFAST's C backend indexes halo arrays with
#: ``int``, so a catalogue whose flattened coordinate array exceeds this is at
#: risk of overflowing regardless of how much memory the node has.
INT32_MAX = 2 ** 31 - 1


# ===========================================================================
#  Run resources
# ===========================================================================

def resolve_n_threads(default: Optional[int] = None) -> int:
    """
    Thread count for 21cmFAST, resolved from the environment.

    Resolution order: the ``N_THREADS`` environment variable, then
    ``SLURM_CPUS_PER_TASK``, then ``default``, then ``os.cpu_count()``, then
    1.  The SLURM variable is preferred over ``cpu_count()`` because on a
    shared node the latter reports the whole machine rather than this job's
    allocation.

    Lives here, rather than in the simulation script, so that the value the
    run uses and the value the manifest records are resolved by one tested
    function.

    Parameters
    ----------
    default : int, optional
        Value to use when neither environment variable is set.  ``None``
        falls through to ``os.cpu_count()``.

    Returns
    -------
    int
        At least 1.

    Notes
    -----
    A non-numeric environment value is ignored rather than raising, so a
    malformed variable cannot abort a queued job.
    """
    for name in ("N_THREADS", "SLURM_CPUS_PER_TASK"):
        raw = os.environ.get(name)
        if not raw:
            continue
        try:
            return max(int(raw), 1)
        except ValueError:
            print(f"  (ignoring non-numeric {name}={raw!r})")

    if default is not None:
        return max(int(default), 1)
    return max(os.cpu_count() or 1, 1)


# ===========================================================================
#  Environment capture
# ===========================================================================

def git_revision(repo_root: str) -> Dict[str, Any]:
    """
    Describe the git checkout the run is using.

    Parameters
    ----------
    repo_root : str
        Directory inside the repository.

    Returns
    -------
    dict
        ``commit``, ``branch`` and ``dirty`` keys.  Values are ``None`` (and
        ``dirty`` is ``None``) when git is unavailable or the directory is not
        a repository — provenance capture must never break a run.
    """
    def _git(*args: str) -> Optional[str]:
        try:
            out = subprocess.run(
                ("git", *args), cwd=repo_root, capture_output=True,
                text=True, timeout=10,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        return out.stdout.strip() if out.returncode == 0 else None

    status = _git("status", "--porcelain")

    return {
        "commit": _git("rev-parse", "HEAD"),
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": None if status is None else bool(status),
    }


def package_versions() -> Dict[str, Optional[str]]:
    """
    Versions of the packages whose behaviour the results depend on.

    Returns
    -------
    dict
        Package name → version string, or ``None`` when not importable.
    """
    versions: Dict[str, Optional[str]] = {}
    for name in ("py21cmfast", "numpy", "scipy", "h5py", "astropy", "hmf",
                 "matplotlib"):
        try:
            module = __import__(name)
            versions[name] = str(getattr(module, "__version__", "unknown"))
        except Exception:                       # noqa: BLE001 — never fatal
            versions[name] = None
    return versions


def environment_info(repo_root: str) -> Dict[str, Any]:
    """
    Host, interpreter, code revision, and package versions.

    Parameters
    ----------
    repo_root : str
        Repository root, passed to :func:`git_revision`.

    Returns
    -------
    dict
        Everything needed to tell two runs of the same script apart.
    """
    return {
        "host": socket.gethostname(),
        "user": os.environ.get("USER") or os.environ.get("USERNAME"),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "cwd": os.getcwd(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "git": git_revision(repo_root),
        "packages": package_versions(),
    }


def peak_memory_gb() -> Optional[float]:
    """
    Peak resident set size of this process, in GB.

    Returns
    -------
    float or None
        Peak RSS, or ``None`` on platforms without ``resource``.

    Notes
    -----
    ``ru_maxrss`` is kilobytes on Linux and bytes on macOS; both are handled.
    """
    try:
        import resource
    except ImportError:                         # pragma: no cover — Windows
        return None

    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    scale = 1e9 if sys.platform == "darwin" else 1e6
    return float(max_rss) / scale


# ===========================================================================
#  Pre-flight cost estimate
# ===========================================================================

def estimate_catalogue_cost(box_len: float) -> Dict[str, float]:
    """
    Extrapolate the halo-catalogue cost of a box from a measured run.

    The halo sampler's floor is a fixed mass (``SAMPLER_MIN_MASS``), not a
    grid property, so the catalogue scales with comoving volume and is
    independent of ``HII_DIM``.  Scaling :data:`HALOS_PER_MPC3` by the volume
    therefore predicts the count for any box at the same sampler settings.

    Parameters
    ----------
    box_len : float
        Comoving box side length [Mpc].

    Returns
    -------
    dict
        ``volume_Mpc3``; ``n_halos_lagrangian``, the catalogue
        ``determine_halo_catalog`` draws, which sets peak memory and index
        width; ``n_halos_perturbed``, the smaller catalogue that survives
        ``perturb_halo_catalog`` and reaches the HDF5 — compare this one
        against a run's ``results.n_halos``; ``catalogue_GB`` for the
        Lagrangian catalogue on disk; ``resident_GB`` for both catalogues,
        as held simultaneously during perturbation; and ``int32_headroom``,
        the flattened ``halo_coords`` length as a fraction of
        :data:`INT32_MAX`.  A headroom above 1.0 means the coordinate array
        is longer than a signed 32-bit index can address.

    Notes
    -----
    An empirical extrapolation, not a model: accurate to the extent that the
    new box shares the reference run's redshift, cosmology and sampler
    settings.  Treat it as an order-of-magnitude guard, not a budget.
    """
    volume = float(box_len) ** 3
    n_lagrangian = HALOS_PER_MPC3 * volume
    catalogue_bytes = n_lagrangian * BYTES_PER_HALO

    return {
        "volume_Mpc3": volume,
        "n_halos_lagrangian": n_lagrangian,
        "n_halos_perturbed": n_lagrangian * PERTURBED_FRACTION,
        "catalogue_GB": catalogue_bytes / 1e9,
        # perturb_halo_catalog holds its input and its output simultaneously.
        "resident_GB": (1.0 + PERTURBED_FRACTION) * catalogue_bytes / 1e9,
        "int32_headroom": n_lagrangian * 3.0 / INT32_MAX,
    }


# ===========================================================================
#  The manifest
# ===========================================================================

class RunManifest:
    """
    A JSON record of one run, rewritten after every update.

    Parameters
    ----------
    path : str
        Destination file.  Its parent directory is created if missing.
    data : dict, optional
        Initial contents.  Normally built by :meth:`create`.

    Attributes
    ----------
    path : str
        The file being maintained.
    data : dict
        The current contents, mirrored to ``path`` on every mutation.
    """

    def __init__(self, path: str, data: Optional[Dict[str, Any]] = None) -> None:
        self.path = os.path.abspath(path)
        self.data: Dict[str, Any] = data if data is not None else {}
        self._started = time.time()
        self._stage_started: Optional[float] = None
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        self.flush()

    # ── construction ──────────────────────────────────────────────────────
    @classmethod
    def create(
        cls,
        output_dir: str,
        label: str = "sim",
        repo_root: Optional[str] = None,
        run_id: Optional[str] = None,
    ) -> "RunManifest":
        """
        Start a manifest named ``<label>_<run_id>.json`` in ``output_dir``.

        Parameters
        ----------
        output_dir : str
            Directory for manifests, e.g. ``outputs/runs``.
        label : str, optional
            Filename prefix identifying the producer.
        repo_root : str, optional
            Repository root for :func:`git_revision`.  Defaults to the parent
            of this file's directory.
        run_id : str, optional
            Identifier for the run.  Defaults to a local-time stamp.

        Returns
        -------
        RunManifest
            An open manifest, already written to disk with
            ``status = "running"``.
        """
        if repo_root is None:
            repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if run_id is None:
            run_id = time.strftime("%Y%m%d_%H%M%S")

        data = {
            "run_id": run_id,
            "label": label,
            "status": "running",
            "stage": None,
            "stages_completed": [],
            "started": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "finished": None,
            "elapsed_seconds": 0.0,
            "command": " ".join([os.path.basename(sys.argv[0]), *sys.argv[1:]]),
            "environment": environment_info(repo_root),
        }
        return cls(os.path.join(output_dir, f"{label}_{run_id}.json"), data)

    # ── mutation ──────────────────────────────────────────────────────────
    def record(self, section: str, values: Mapping[str, Any]) -> None:
        """
        Merge ``values`` into a named section and rewrite the file.

        Parameters
        ----------
        section : str
            Top-level key, e.g. ``"parameters"`` or ``"results"``.
        values : mapping
            Entries to merge.  Existing keys are overwritten.
        """
        self.data.setdefault(section, {}).update(dict(values))
        self.flush()

    def begin_stage(self, name: str) -> None:
        """
        Mark the start of a stage.

        A crashed run leaves this name in ``stage``, which is what identifies
        where it died.

        Parameters
        ----------
        name : str
            Stage name, e.g. ``"lightcone"``.
        """
        self.data["stage"] = name
        self._stage_started = time.time()
        self.flush()

    def end_stage(self) -> float:
        """
        Close the current stage and record how long it took.

        Returns
        -------
        float
            Stage duration in seconds; ``0.0`` if no stage was open.
        """
        name = self.data.get("stage")
        if name is None or self._stage_started is None:
            return 0.0

        elapsed = time.time() - self._stage_started
        self.data.setdefault("timings_seconds", {})[name] = round(elapsed, 2)
        self.data.setdefault("stages_completed", []).append(name)
        self.data["stage"] = None
        self._stage_started = None
        self.flush()
        return elapsed

    def finish(self, status: str = "complete") -> None:
        """
        Close the manifest.

        Parameters
        ----------
        status : str, optional
            Final status, conventionally ``"complete"`` or ``"failed"``.
        """
        self.data["status"] = status
        self.data["finished"] = time.strftime("%Y-%m-%d %H:%M:%S %Z")
        self.flush()

    def flush(self) -> None:
        """
        Write the manifest to disk, replacing any previous contents.

        Written to a temporary file and renamed, so a crash mid-write cannot
        leave a half-parsed manifest behind.
        """
        self.data["elapsed_seconds"] = round(time.time() - self._started, 2)
        peak = peak_memory_gb()
        if peak is not None:
            self.data["peak_memory_GB"] = round(peak, 3)

        temporary = f"{self.path}.tmp"
        with open(temporary, "w") as stream:
            json.dump(self.data, stream, indent=2, default=str)
            stream.write("\n")
        os.replace(temporary, self.path)
