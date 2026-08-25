# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased]

<!-- ─── Planned-run specification, 2026-08-24 ───────────────────────────── -->

### Added

- **`docs/simulation_spec.md` — the specification and cost of the run that has
  not happened yet.** `docs/HPC.md` documents the run *as the code stands*;
  this is the forward-looking companion: the production lightcone at
  *z* = 6.55–7.45 in the footprint-derived 486.33 Mpc / 256³ box, with every
  parameter, every derived number, and what it will cost to run.

  Every number carries its provenance — **[M]** measured, **[E]** extrapolated,
  **[C]** computed for the document, **[T]** a target — because the mix was
  previously invisible and the extrapolations are load-bearing.

  | Requirement | Figure |
  |---|---|
  | Cores | 16 per task **[T]** |
  | Memory | 128 GB requested; ~56 GB peak resident **[E]** |
  | Scratch | 300 GB requested; ~245 GB cache + ~19 GB output **[E]** |
  | Wall time | request 6 h; 0.5–1.5 h on 16 cores, 2–4 h on 1 core **[E]** |

  The wall-time extrapolation rests on the only end-to-end timings this project
  has: **543.5 s wall / 567.4 s CPU** for the 2026-08-12 run at 256 Mpc / 128³
  / 5 nodes, effectively single-threaded, reproduced by the 2026-08-07 run to
  within one second **[M]**. Stage 1 is 96 % of it (520.1 s); the power spectra
  took 0.8 s. Scaling is 6.86× (volume) × 1.8× (nodes) ≈ 12.3× of serial work.

  The document also carries the derived geometry of the planned range **[C]**:
  L_LOS = 315.598 Mpc, `N_z` = 166, 9 node redshifts with *z* = 7.0 landing
  exactly on node 5, 168.095–188.133 MHz (20.038 MHz), Δk⊥ = 0.01292 and
  Δk∥ = 0.01991 Mpc⁻¹; a SLURM template; and §9, which lists what is needed
  about a target cluster to turn the wall-time range into a number, with the
  commands that produce it (`lscpu`, `sinfo`, `scontrol show partition`,
  `lfs quota`, `sacct`, `seff`).

### Two findings worth recording

- **The `L_los` discrepancy does not survive the wider range.** `HPC.md` §11.1
  records the stored `L_los` (200.0 Mpc) disagreeing with the data's actual
  3.5 Mpc span by a factor 56.5. At *z* = 6.55–7.45 the recorded value —
  `N_z` × transverse cell — is **315.354 Mpc against a true 315.598, a 0.08 %
  difference [C]**. The factor 56.5 was `minimum_los_slices` flooring a
  two-slice slab to 100, not a defect in the geometry, and it disappears the
  moment the natural `N_z` binds.

- **Raising the sampler floor is a cheap way past `INT_MAX`, and now has a
  number attached.** At 486.33 Mpc the flattened `halo_coords` is 1.31×
  `INT_MAX` **[E]** — the 2026-08-20 SIGSEGV. Against a Sheth-Tormen mass
  function at *z* = 7 weighted by 21cmFAST's own
  f_★ ∝ (M/10¹⁰)^0.5 exp(−`M_TURN`/M), moving `SAMPLER_MIN_MASS` from 10⁸ to
  **2 × 10⁸ M☉ keeps 46.2 % of the halos and 99.84 % of the star formation
  [C]** — headroom 0.60, catalogue 12.1 GB. The halos it drops sit far below
  the M_TURN = 5 × 10⁸ turnover. The document still recommends the 350 Mpc
  box (headroom 0.49) as the first attempt, because that leaves the physics
  untouched and calibrates the extrapolations; the mass floor is the fallback
  if the full 10 deg² footprint is required.

### Changed

- **Costed against the target machine, 2026-08-24.** `lscpu` and `df -hT` on
  the machine the run will use replaced most of the document's [T] guesses
  with [M] facts, and one of them changes the shape of the run:

  | | |
  |---|---|
  | Host | `andromeda1.jb.man.ac.uk` — **the same machine the 2026-08-12 baseline was measured on** |
  | CPU | 2 × AMD EPYC 9374F (Zen 4), 64 physical cores / 128 logical, 3.85 GHz base |
  | RAM | 1.5 TiB, 1.4 TiB available (`free -h`) — the 56 GB peak is 1/25 of it |
  | Scratch | local NVMe/XFS: `/nvme1` 851 GB free, `/nvme4` 836 GB free |
  | Home | NFS, **144 GB free** — cannot hold the 245 GB cache |
  | Scheduler | **none** — `sinfo`, `scontrol`, `sacct`, `sacctmgr` all absent |

  Because the host is the baseline's own host, the serial scaling carries **no
  cross-machine correction** — the extrapolation is across box sizes only, and
  the sole remaining assumption in the wall-time figure is the OpenMP
  speed-up, which one smoke-test run settles.

  Consequences: the wall-time estimate drops to **15–40 min at 32 threads**
  (budget 1 h) from 0.5–1.5 h on a hypothetical 16-core node; memory stops
  being a constraint at any box size; §8 becomes a `tmux` + `numactl` launch
  recipe rather than a SLURM script (the `#SBATCH` version is kept folded away
  for a future move to CSD3); and I/O drops to ~2–4 min on local NVMe.

  **Two new hazards, both specific to an unscheduled machine.** `R7` —
  `resolve_n_threads()` falls through to `os.cpu_count()` when there is no
  `SLURM_CPUS_PER_TASK`, which is **128** here: every SMT thread of a shared
  workstation. `N_THREADS=32` must be set explicitly. `R3` — `$HOME` is NFS
  with 144 GB free, so a run launched from the home directory dies partway
  through writing a 245 GB cache.

  The 1.5 TiB does **not** unblock the 486.33 Mpc box: `halo_coords` is
  1.31 × `INT_MAX` whatever the node has. Abundant memory removes a worry the
  2026-08-20 post-mortem had to rule out; it does not remove the blocker.

- **`NUMBERS_AND_SOURCES.md`**: the two `M_cell` entries in §1 still quoted the
  retired 256 Mpc / 128³ grid (`1.175e10` / `3.173e11` M☉) and now carry the
  footprint-derived values (**`1.007e10` / `2.720e11` M☉**, `DIM` = 768 /
  `HII_DIM` = 256), with the old numbers named. New **§10** records the planned
  run's derived quantities and cost estimates with the same [M]/[E]/[C]
  provenance marks.

- **`README.md`**: `docs/simulation_spec.md` added to the documentation table
  and cross-referenced from the derived-geometry note in the configuration
  section.

<!-- ─── Mode weighting + physical HERA noise, 2026-08-21 ────────────────── -->

### Added

Both discrepancies found by the literature validation below are now
implemented, and **both are opt-in** — the defaults reproduce every number
this pipeline produced before them, so no stored result is silently
invalidated. `UncertaintyBudget` records which was used, and both reach
`pipeline_summary.json`.

- **`--mode-weighted`** — La Plante et al. (2023) Eq. 19's
  $\sqrt{N_\mathrm{patch}\,dN}$ weighting, using the estimator's own
  `mode_counts`. `dN` is `mode_counts / 2`, since the FFT of a real field is
  Hermitian; verified against Eq. 18 to a median 4.7 % on the fiducial grid.
  `cross_power_snr` gains `mode_counts` and `n_patch`, both defaulting to the
  unweighted behaviour.

- **`--noise-model physical`** — a real HERA instrument model replacing the
  flat estimate:

  $$P_N(k_\perp) = \frac{X^2 Y \,\Omega_\mathrm{eff}\, \mathrm{NEB}\, T_\mathrm{sys}^2}{N_\mathrm{pol}\, t_\mathrm{int}\, N_\mathrm{bl}(k_\perp)}$$

  This is Parsons (2017) HERA memo Eq. 12, algebraically identical to
  La Plante Eq. 11 since $\Omega_\mathrm{eff} \equiv \Omega_P^2/\Omega_{PP}$.
  Four new functions: `cosmological_scalar_x2y`, `hera_beam_solid_angles`,
  `hera_baseline_counts`, `hera_thermal_noise_power_physical`.
  `hera_thermal_noise_power` is untouched and remains the default.

  The baseline count is computed, not assumed: an 11-per-side hexagonal core
  at 14.6 m spacing (DeBoer et al. 2017) is built, every pair formed, and the
  baselines gridded into $uv$-cells one antenna footprint $D/\lambda$ across.
  **Per uv-cell, not per $k_\perp$ bin** — only baselines within one footprint
  sample the same mode and integrate down coherently; summing the whole bin
  would double-count the redundancy that the mode weighting already handles,
  and under-predicted the noise several-fold in a first pass.

  At $z = 7$ for 1000 h this gives 290 baselines per cell at
  $k_\perp = 0.010$ falling to 12 at $k_\perp = 0.123$, hence
  $P_N = 3.9\times10^3 \to 9.8\times10^4$ mK² Mpc³ — bracketing the
  $\approx 2.5\times10^4$ implied by published HERA forecasts, against the
  default estimate's 3.75. Beyond $k_\perp \approx 0.13$ Mpc⁻¹ the 292 m core
  has no baselines at all and the noise is `inf`: those modes are
  *unmeasurable*, a statement the flat estimate cannot make.

  Measured on a $48^2\times100$, 96 Mpc test run at $\sigma_z = 0.02$:

  | Configuration | Total SNR | vs default |
  |---|---|---|
  | default | 0.0500 σ | — |
  | `--mode-weighted` | 0.1420 σ | ×2.8 |
  | `--noise-model physical` | 0.0056 σ | ÷8.9 |
  | both | 0.0159 σ | ÷3.1 |

### Cited

- **`NUMBERS_AND_SOURCES.md` extended to cover every value introduced today**,
  and six previously unsourced entries traced:

  | Value | Now cited to |
  |---|---|
  | `T_RECEIVER_K`, `T_SKY_300MHZ_K`, `SKY_SPECTRAL_INDEX`, 300 MHz pivot | 21cmSense `calc_sense.py` (`Tsky = 60e3*(3e8/freq)**2.55`), from Pober et al. (2013, 2014) |
  | `dish_diameter`, array layout | DeBoer et al. (2017), PASP 129, 045001 |
  | `F_21_MHZ` | Hydrogen hyperfine frequency |

  New §7 (HERA instrument model, 6 constants + 4 formulas), §8 (foreground
  module, 11 constants + 3 formulas), §9 (run-manifest calibration, 4
  measured values, explicitly labelled internal rather than literature). Every
  constant also carries its citation inline in the code.

- **A convention difference worth knowing.** The `T_sys` model here is
  21cmSense's `T_sky = 60 K (λ/1 m)^2.55` — 352 K at 150 MHz. DeBoer et al.
  (2017) Table 2 instead give `T_sys = 100 + 120(ν/150 MHz)^-2.55` K, whose
  sky term is 120 K at 150 MHz: **a factor 2.9 colder, ~8× in $P_N$.** The
  21cmSense values are kept — they are what the notebook inlined and the more
  conservative choice — but the difference is now recorded at the constants,
  in `system_temperature`'s docstring, and in `NUMBERS_AND_SOURCES.md` §3, so
  it is not mistaken for an error.

### Fixed

- **`pipeline_summary.json` was being written truncated and unparseable** for
  any HDF5 produced by the current `run_simulation.py`. The `source_run` block
  added earlier today reads `random_seed` and `n_threads` straight from the
  HDF5, and `h5py` returns those as `np.int64`. `json.dump` cannot serialise a
  numpy scalar and raises *partway through the write*, leaving a file that
  ends mid-key. It went unnoticed because the test fixture's HDF5 predates
  those attributes, so `.get()` returned `None` and the path never ran.

  Attributes are now coerced to plain Python scalars in `build_summary`, and
  `write_summary` writes to a temporary file and renames, so a failed write
  can no longer leave a half-valid summary behind. Regression test:
  `test_summary_is_valid_json_with_provenance_attrs`, which stamps the
  attributes onto the fixture first.

<!-- ─── Literature validation of the budget, 2026-08-21 ─────────────────── -->

### Validated

- **The uncertainty budget and SNR chain checked term by term against the
  cited papers** (La Plante et al. 2023, arXiv:2205.09770; Pober et al. 2014).
  Full write-up in [`docs/uncertainty_budget.md`](docs/uncertainty_budget.md)
  §6.5. **No code changed** — the two discrepancies below would move every
  number the pipeline has produced, which is a project decision.

  | Term | Cited as | Verdict |
  |---|---|---|
  | Horizon wedge slope | La Plante Eq. 10 | ✅ exact |
  | Cross-spectrum variance | La Plante Eq. 15 | ✅ exact in form |
  | Galaxy variance | La Plante Eq. 17 | ✅ exact |
  | $T_0(z)$ cancellation | La Plante Eq. 16 | ✅ algebra confirmed |
  | $T_\mathrm{sys}$, wedge buffer | Pober et al. (2014) | ✅ correct |
  | Mode weighting | La Plante Eq. 19 | ❌ omitted |
  | Thermal noise | La Plante Eq. 11 | ❌ placeholder, ~10⁴ low |

  The horizon slope is worth recording as *identical*, not approximate: the
  paper's $m = \lambda D_c f_{21} H/[c^2(1+z)^2]$ reduces to the implemented
  $D_c H/[c(1+z)]$ on substituting $\lambda = c(1+z)/f_{21}$.

### Two discrepancies quantified

- **The Eq. 19 mode weighting is missing, and it is worth ~an order of
  magnitude.** La Plante combines bins as
  $\hat{s} = \sqrt{N_\mathrm{patch}\,dN}\,P_\times/\sigma_\times$; the pipeline
  computes $P_\times/\sigma_\times$ and stops. The omitted $dN$ is already on
  hand — `mode_counts / 2` (the FFT of a real field is Hermitian) reproduces
  Eq. 18 to a median 4.7 %. Over the 65 usable bins of the fiducial
  $128^2\times100$ grid, $\sqrt{dN}$ runs 2.0–64.1 with median 8.2, so the
  reported total SNR is low by roughly 18×. Previously logged as
  "conservative" in §7.3 without a number; the number makes it more than a
  footnote.

- **`hera_thermal_noise_power` is ~4 orders of magnitude too small, in the
  optimistic direction.** It returns 3.75 mK² Mpc³ at $z=7$ for 1000 h.
  La Plante Eq. 11 evaluated for HERA ($X^2Y = 1227$ Mpc³ sr⁻¹ Hz⁻¹,
  $\Omega_p^2/\Omega_{pp} \approx 0.19$ sr) gives $5\times10^2$–$5\times10^4$
  mK² Mpc³ depending on $N_\mathrm{bl}(u)$, and published HERA forecasts
  ($\Delta^2_N \sim 10$ mK² at $k=0.2$) independently give
  $\approx 2.5\times10^4$. At the implemented value the 21 cm side of
  $\sigma_{21} = |P_{21}| + P_{N,21}$ is sample-variance limited, which is not
  true of HERA at 1000 h.

  The docstring already called this "a scaling estimate, not a full instrument
  model" — accurate, but it does not convey a factor of $10^4$, so §7.2 now
  states the number.

  The two pull in opposite directions and do not cancel in general; applying
  both would move the fiducial total SNR up by roughly 5×.

### Fixed

- **Author list for La Plante et al. (2023)** in
  `docs/uncertainty_budget.md` §9 — recorded as "La Plante, Kaur, Battaglia
  et al."; the paper is La Plante, Mirocha, Gorce, Lidz & Parsons,
  *"Prospects for 21cm-Galaxy Cross-Correlations with HERA and the Roman
  High-Latitude Survey"*. Title and Eq. numbers now cited alongside.

<!-- ─── Foreground injection and removal, 2026-08-21 ────────────────────── -->

### Added

- **`src/foregrounds.py`** — synthetic foreground injection and a
  parametrised removal knob, so the effect of contamination and of
  *incomplete* removal on the cross-power spectrum and its detectability can
  be measured. Nothing in the verified chain is modified: contaminated fields
  go through `compute_all_power_spectra` and `compute_uncertainty_budget`
  exactly as clean ones do.

  | Function | Description |
  |---|---|
  | `simulate_diffuse_foreground` | Diffuse Galactic synchrotron cube [mK] on the lightcone grid |
  | `simulate_point_source_foreground` | Poisson point-source component, angularly flatter |
  | `inject_foreground` | Combines both, scales to `foreground_amplitude × signal RMS` |
  | `remove_foreground` | The removal knob — a placeholder, see below |

  The model is the standard parametrised angular/frequency power law
  (Santos, Cooray & Knox 2005; Shaw et al. 2014) with amplitude and spectral
  index from the Global Sky Model (de Oliveira-Costa et al. 2008; Zheng
  et al. 2017, GSM2016); point sources follow Ali, Bharadwaj & Chengalur
  (2008) and Bernardi et al. (2009). The line-of-sight structure is
  *emergent* — it follows from a smooth power-law spectrum with a spatially
  varying index, rather than from an imposed $k_\parallel$ power law.

- **Notebook §7e**, after the §7d wedge cell: injection at a chosen
  `FOREGROUND_AMPLITUDE`, a sweep over `REMOVAL_FRACTIONS`
  (0/50/90/99/99.9/100 %), and three panels — $P_{21}$ showing where the
  foreground sits and the removal working, $|P_\times|$ against
  $\sigma_\times$, and total SNR versus removal fraction.

### Two findings worth recording

- **Spectral smoothness does not confine the contamination to low
  $k_\parallel$.** `compute_cylindrical_cross_power` takes a bare FFT with no
  line-of-sight taper, so a smooth-but-non-periodic spectrum is discontinuous
  at the box edge and leaks along the whole $k_\parallel$ axis as
  $\approx k_\parallel^{-1.5}$. The slope is a property of the window, not the
  sky: a bare linear ramp with no angular structure leaks with the same slope,
  and widening the band does not change it. **Wedge excision alone therefore
  does not remove foreground power from the EoR window here.** Measured, not
  assumed — `test_smoothness_does_not_confine_power_to_low_k_parallel` pins it
  against the bare-ramp reference.

- **A foreground-contaminated SNR flatters the result.** Foregrounds are
  unbiased in the ensemble mean, but a single realisation carries a chance
  cross-correlation of order $\sqrt{P_{21}P_{\rm gal}/N_{\rm modes}}$. On the
  synthetic fixture the per-bin shift in $P_\times$ scales **linearly** with
  foreground amplitude (×10 per decade) while the shift in $P_{21}$ scales
  **quadratically** (×100) — exactly the contrast expected if the foreground
  reaches the cross-spectrum only by chance. Because $|P_\times|/\sigma_\times$
  then has a contaminated numerator *and* denominator, the total SNR degrades
  far more slowly than $\sigma_\times$ alone: at $10^4\times$ contamination
  $\sigma_\times$ rose 583× while the total SNR fell only 6.8×, and the
  as-measured SNR sat at 15 % of the clean value where the surviving *signal*
  was 0.16 %.

  §7e.3 therefore plots a **signal-only** SNR (clean $P_\times$ against
  contaminated $\sigma_\times$) beside the as-measured one; the shaded gap
  between them is a spurious detection. An earlier draft of this section
  claimed foregrounds "widen the error bar but do not move the measurement" —
  true in the mean, wrong for a single realisation at high contamination, and
  corrected throughout.

### Note on `remove_foreground`

It subtracts an exactly-correct template of the very field that was injected,
scaled by `removal_fraction`. It is **not** GMCA, PCA, ICA, polynomial or
log-polynomial fitting, Gaussian-process removal, or delay filtering, and has
none of their failure modes — no signal loss from over-fitting, no
mode-mixing, no leakage from the wedge into the window, no dependence on the
foreground's smoothness. Results are statements about *removal level*, not
about any named method's achievable performance. The docstring, the module
header, the notebook markdown and `docs/reference.md` all say so explicitly,
so it cannot be mistaken later for something more sophisticated.

`removal_basis` resolves the amplitude/power ambiguity: the default
`"amplitude"` scales the residual by $(1-f)$, so residual *power* falls as
$(1-f)^2$; `"power"` scales it by $\sqrt{1-f}$.

<!-- ─── Run manifests + crash mitigations, 2026-08-21 ───────────────────── -->

### Context — the 2026-08-20 SIGSEGV

`run_pipeline.py --sim force` died after 2,303 s with `exit code -11`
(SIGSEGV) and **no output at all** from the simulation child. Both facts are
now addressed.

The silence was buffering. `run_pipeline.py`'s `log()` uses
`print(..., flush=True)`; `run_simulation.py` had no `flush=True` anywhere, so
~8 KB of the child's progress output sat in a block buffer that the signal
discarded. The crash destroyed its own diagnostics.

The crash was scale. Commit `81e08ef` (2026-08-20 12:39, three hours earlier)
replaced the hardcoded `HII_DIM = 128 / BOX_LEN = 256 / DIM = 384` with the
footprint-derived `256 / 486.33 / 768`. Calibrated against the 2026-08-12 run
and its 21cmFAST cache (136,663,818 halos in a 3.564 GiB `HaloCatalog.h5` —
28.0 bytes/halo), that is 6.86× the volume: ~9.4 × 10⁸ halos, a 26.2 GB
catalogue, 48–52 GB resident during `perturb_halo_catalog`, and a flattened
`halo_coords` **1.31× past `INT_MAX`**. Of the 16 cached `HaloCatalog.h5`
files in the tree exactly one is unreadable, and it is the one written at
`BOX_LEN = 486.33` on a 512× smaller grid under no memory pressure — stopping
dead at 2,147,491,839 bytes, the signed 32-bit boundary. Full analysis in
[`docs/HPC.md`](docs/HPC.md) §13.5.

### Added

- **`src/provenance.py` — per-run parameter manifests.**
  `run_simulation.py` now writes `outputs/runs/sim_<run_id>.json` **before**
  the expensive stages and rewrites it after each one, so it survives a run
  that never finishes. A process killed by a signal cannot flush stdout or run
  an exit hook; the manifest is left with `"status": "running"` and `"stage"`
  naming exactly where it died.

  It records `parameters` (every configuration constant), `derived` (geometry,
  mass resolution, LOS slicing), `cost_estimate`, `environment` (host, git
  commit and dirty flag, package versions, SLURM job id), `timings_seconds`
  per stage, `peak_memory_GB`, `results`, and `outputs`. Written via a
  temporary file and `os.replace`, so a crash mid-write cannot leave
  half-parsed JSON.

- **A pre-flight halo-catalogue cost estimate** (`estimate_catalogue_cost()`),
  printed and recorded before any compute is spent. The sampler's floor is a
  fixed mass, not a grid property, so the catalogue scales with comoving
  volume — meaning a modest-looking `BOX_LEN` change is a large cost change.
  Past `INT_MAX` it warns explicitly that more memory will not help:

  ```
  Est. halos  : 9.370e+08 drawn (7.836e+08 after perturbation) in 1.150e+08 Mpc³
                →  26.2 GB on disk, ~48.2 GB resident while perturbing
    *** WARNING: halo_coords would hold 2.811e+09 elements, 1.31x INT_MAX ***
  ```

  It reports the Lagrangian and perturbed counts separately: the first sets
  peak memory and index width, the second is what reaches the HDF5 and is
  directly comparable to a run's `results.n_halos`. The 83.6 % ratio comes
  from the 256 Mpc run and was reproduced by a held-out 64 Mpc run
  (1,782,540 actual against 1,785,000 predicted).

- **`RANDOM_SEED` as a named constant** (was a literal `42` inside the
  `from_template` call), recorded in the manifest and the HDF5 attrs alongside
  `n_threads`, `minimize_memory`, `run_id` and `run_manifest`.

- **`source_run` in `pipeline_summary.json`** — `run_id`, `run_manifest`,
  `random_seed` and `n_threads` read back from the HDF5, so an analysis-only
  run names the simulation its numbers came from. `null` for files written
  before manifests existed.

- **`submit_job.sh` reports the newest manifest in its email**, including a
  `DIED IN STAGE:` line when the run did not close it.

### Changed

- **`python -u` for the simulation child**, in `submit_job.sh` and in
  `run_pipeline.py`'s `subprocess.run` call. The pipeline no longer depends on
  how it was itself launched for its child's output to reach the log.

- **`del halo_catalog, initial_conditions`** immediately after
  `perturb_halo_catalog` returns. That call is the memory high-water mark of
  the script — the Lagrangian catalogue and its perturbed copy are both
  resident — and nothing below reads either input again. Frees ~26 GB of
  catalogue and ~7.7 GB of ICs at the peak on the 486.33 Mpc box.

- **`MINIMIZE_MEMORY = True`** in `matter_options`, trading peak RAM for
  intermediate I/O in the C backend. It was `False` (21cmFAST's default) in
  every run so far.

- **`N_THREADS` is now set**, via `provenance.resolve_n_threads()`:
  `N_THREADS` env → `SLURM_CPUS_PER_TASK` → `os.cpu_count()` → 1. The SLURM
  variable is preferred over `cpu_count()` because on a shared node the latter
  reports the whole machine rather than the job's allocation. 21cmFAST's own
  default is 1, which is why the failed run showed `user/real = 0.94`. A
  non-numeric environment value is ignored rather than raised, so a malformed
  variable cannot abort a queued job.

- **`.gitignore` covers all hash-named 21cmFAST cache directories**, not just
  the one that existed when the rule was written. `a14661e5…/` — which holds
  the corrupt 2 GiB `HaloCatalog.h5` — had been showing up as untracked.

### Note

None of this makes the 486.33 Mpc box fit. It makes the next attempt
diagnosable and cheaper, and it says up front when a box cannot work. The
suggested next step is a ~350 Mpc intermediate (2.5× volume, ~340 M halos,
~9.5 GB catalogue, 0.48× `INT_MAX`) rather than the full footprint.

<!-- ─── Post-Euclid-cut figures, 2026-08-20 ─────────────────────────────── -->

### Added

- **Three post-Euclid-cut figures**, in both front ends. Everything the
  pipeline plotted about the galaxy population was either the *full* halo
  catalogue or the stored `galaxy_overdensity`; neither is what the survey
  sees. `uv_selection_maps` showed *where* the selected galaxies sit, but no
  figure showed the selected population's own distributions, its overdensity
  field, or that field against the 21 cm signal.

  | Figure | `src/figures.py` | Panels |
  |---|---|---|
  | `euclid_selected_catalogue` | `plot_euclid_selected_catalogue` | selected galaxies on the sky coloured by $M_\mathrm{UV}$; halo mass before/after the cut; SFR before/after, with the equivalent SFR window marked |
  | `selected_galaxy_overdensity` | `plot_selected_galaxy_overdensity` | LOS projection of $\delta_\mathrm{gal}$; a single transverse slice; the one-point distribution on a symlog axis |
  | `galaxy_overdensity_on_21cm` | `plot_galaxy_overdensity_on_21cm` | $\delta_\mathrm{gal}$ contours over $\delta T_b$; the same maps with the roles swapped; $\langle\delta_\mathrm{gal}\rangle$ binned by $\delta T_b$, with the cell-by-cell Pearson $r$ |

  These come with a new figure group, **`--plots euclid`** (catalogue-dependent,
  skipped with a message when the HDF5 has no catalogue), and a new
  **`--galaxy-weighting {number,luminosity}`** flag. Figure count 15 → 18.

- **`figures.selected_galaxy_overdensity()`** — the shared rebuild behind those
  last two figures. It calls `analysis.galaxy_overdensity_from_catalogue()` on
  `run_simulation.py` §3b's grid (`n_perp = HII_DIM`, `n_los = N_z`,
  `los_extent = BOX_LEN`, i.e. the *coeval* box, not `L_los`) and applies the
  Euclid window.

  **This is deliberately not the stored field.** The default
  `GALAXY_WEIGHTING = "lightcone_sfr"` builds `galaxy_overdensity` from the
  lightcone `halo_sfr` field, which applies **no magnitude cut at all** — so
  plotting the stored field would not have been "after the Euclid cut" in any
  run using the shipped default. `run_pipeline.py` deposits it once and hands
  it to both figures rather than binning a 114 M-halo catalogue twice.

- **Notebook §3c and §5c** in `21cmfast_HERAxEuclid_lightcone.ipynb`: §3c.1
  plots the post-cut galaxies/halos/SFR, §3c.2 rebuilds and plots
  $\delta_\mathrm{gal}$ (weighting switchable via
  `GALAXY_WEIGHTING_DIAGNOSTIC`), and §5c overlays it on the 21 cm field. The
  overlay lives in §5 rather than §3 because it needs `eor_cmap`, which §5b
  defines. The cells import `select_euclid_halos` and
  `galaxy_overdensity_from_catalogue` from `src.analysis` — the same calls the
  batch figures make — so the two front ends cannot drift.

### Notes on reading the new figures

- **The selected sample is shot-noise dominated on the fiducial grid.** 49,315
  of 114 M halos survive the cut: 0.03 galaxies per cell over
  $128^2 \times 100$ cells, with ~97.5 % of cells empty. A single transverse
  slice of $\delta_\mathrm{gal}$ is therefore almost pure noise, which is why
  `plot_selected_galaxy_overdensity` leads with the LOS projection and
  `plot_galaxy_overdensity_on_21cm` averages over `slab_cells = 8` LOS cells
  and Gaussian-smooths the contoured field (`smooth_cells = 2`, display only).
  The slice panel is kept because it is what the power-spectrum estimator is
  actually handed.
- **The overlays are transverse only.** The halo catalogue is a coeval snapshot
  at $z_\mathrm{obs}$ and the 21 cm field is a lightcone, so the two share the
  $(x, y)$ plane and the array shape but not an LOS scale — the same mismatch
  `run_simulation.py` already accepts in catalogue mode. The third panel's
  Pearson $r$ pairs cells exactly as `compute_all_power_spectra` does, so it is
  the real-space counterpart of the cross-power sign, mismatch included. At the
  fiducial parameters $r = -0.010$, and $\langle\delta_\mathrm{gal}\rangle$
  falls from $+0.03$ at $\delta T_b \approx 0$ to $-0.14$ in the brightest
  bin — galaxies avoid the neutral gas, as the negative large-scale
  cross-power says.

### Changed

- **`figures.plot_uv_selection_maps`** now lifts its selection mask through the
  new `_lift_selection_mask()` helper instead of inlining the
  valid-subset-to-full-catalogue index lift. Same result; the logic is no
  longer duplicated between two figures.

<!-- ─── README streamlined + quickstart, 2026-08-20 ────────────────────── -->

### Added

- **`docs/reference.md`** — the long-form companion to the README, holding the
  detail that used to sit inline there: per-notebook structure, equations and
  fiducial parameters; the figure-by-figure literature references; the
  21cmFASTv4 `HaloBox` / lightcone API notes and source-model templates; the
  `src/` function reference; the requirements table; the test-suite coverage
  table; the `%matplotlib widget` notes; and the full bibliography. Content is
  carried over unchanged apart from the corrections noted below.

- **A `## Quickstart` section in `README.md`**, covering *all* options of both
  front ends:

  | Subsection | Contents |
  |---|---|
  | Install | `env.yml`, the pinned freeze, and the HPC pointer |
  | Run the HPC pipeline | stage control, the seven `--plots` groups and the figures each writes, paths/rendering options, the four uncertainty-budget overrides, `submit_job.sh`'s environment variables, and `run_simulation.py`'s configuration constants |
  | Run the notebooks | which notebook needs a simulation, what to edit in each, and the full ★ CONFIGURATION cell of `21cmfast_HERAxEuclid_lightcone.ipynb` |
  | Outputs | the four output paths and the runtime note |

  `run_simulation.py` has **no CLI** — it is a configuration-constant script, so
  its parameters are documented as a constants table rather than as flags.

### Changed

- **`README.md` cut from 929 to ~310 lines.** The overview is now a three-row
  table (pipeline / notebooks / `src/`), and the reference material moved to
  `docs/reference.md`. The documentation index gained a row for the new doc.

### Fixed

- **`--plots` was documented without the `budget` group.** The README listed
  `{all,none,fields,halos,scaling,power,snr,bias}`; the actual choices are
  `{all,none,fields,halos,scaling,power,snr,budget,bias}`, with `budget`
  writing `uncertainty_budget` and `photoz_suppression`. The four
  uncertainty-budget overrides (`--sigma-z`, `--wedge-buffer`,
  `--integration-time`, `--bandwidth`), `--sim-script` and `--m-uv-bright` were
  undocumented entirely.

- **Notebook 3's magnitude window and redshift range.** The README recorded
  only "$M_\mathrm{UV}$ limit $< -18$"; the notebook's configuration cell sets
  `M_UV_faint = -18.0` **and** `M_UV_bright = -22.66`, which differs from
  `run_simulation.py`'s $-22$. The prose also still said the lightcone spans a
  "default $z = 6.5$–$7.5$", whereas the range has been derived from the survey
  footprint as $6.55$–$7.45$ since the footprint-driven box change.

- **The `TODO.md` index row** referred to "the Δz = 1.0 range now requires";
  the committed range is the Δz = 0.01 smoke-test slab, so the row now reads
  "a wider Δz would require".

<!-- ─── Footprint-driven simulation box size, 2026-08-19 ───────────────── -->

### Added

- **`survey_area_to_box_size()` in `src/conversions.py`** — converts a survey
  footprint (area in deg², central redshift, redshift depth) into a 21cmFAST
  box geometry, so `BOX_LEN` is traceable to the survey being forecast instead
  of chosen by hand.

  ```python
  survey_area_to_box_size(area_deg2, z_central, delta_z, cosmo=None,
                          target_cell_size_mpc=2.0, hii_dim=None,
                          snap_hii_dim_to_power_of_two=True)
  ```

  | Step | Formula | Fornax result |
  |---|---|---|
  | Transverse extent | $L_\perp = \sqrt{\Omega}\,D_M(z_c)$ (small-angle, square footprint) | 486.33 Mpc |
  | LOS depth | $L_\parallel = D_C(z_c + \Delta z/2) - D_C(z_c - \Delta z/2)$ | 315.60 Mpc |
  | Grid | $\lceil L_\perp / 2.0\,\mathrm{Mpc}\rceil$, snapped to a power of two | 244 → 256 |

  The LOS depth uses **distance differencing**, not
  $\mathrm{d}D_C/\mathrm{d}z \times \Delta z$, matching how
  `run_simulation.py` §1 and the notebook already compute `L_los`. The
  transverse step reuses the same construction as
  `src/FOV_to_cMpc.py:transverse_comoving_size_from_area`, and `cosmo=None`
  defaults to astropy `Planck18`, following the convention already used by
  `volume_from_area` and `survey_area_from_volume`.

- **`SimulationBox` dataclass**, returned by the above. Carries `box_len`,
  `hii_dim`, `dim`, `cell_size`, `los_depth`, `z_min`, `z_max`,
  `transverse_area`, `solid_angle_sr`, and `n_los_tiles`. Its
  `.simulation_options` property is the `{"HII_DIM", "BOX_LEN", "DIM"}` mapping
  `p21c.InputParameters.clone()` expects, so callers never restate it.

  `los_depth` is reported but is **not** a 21cmFAST argument — boxes are cubic
  and a lightcone's LOS extent comes from the redshift range passed to
  `RectilinearLightconer`; `z_min` / `z_max` are provided for that.
  `n_los_tiles` = $L_\parallel / L_\perp$ flags when the coeval box would be
  tiled along the LOS (0.65 for Fornax — no tiling).

- **11 tests in `tests/test_conversions.py`** covering the transverse and LOS
  formulae against independent astropy calculations, the power-of-two snap,
  mass-resolution preservation, the `simulation_options` mapping, $\sqrt{A}$
  scaling, and seven rejected non-physical inputs.

### Changed

- **`BOX_LEN` / `HII_DIM` / `DIM` are no longer hardcoded** in either
  `run_simulation.py` or `21cmfast_HERAxEuclid_lightcone.ipynb`. Both now derive
  the grid from a new survey-footprint configuration block:

  | Parameter | Before | After |
  |---|---|---|
  | `BOX_LEN` | 256.0 Mpc (hand-picked) | **486.33 Mpc** (10 deg² Fornax at $z=7$) |
  | `HII_DIM` | 128 | **256** |
  | `DIM` | 384 | **768** |
  | `HII_DIM` cell size | 2.00 Mpc | 1.90 Mpc |
  | $M_\mathrm{cell}$ (`HII_DIM`) | $3.17\times10^{11}\ M_\odot$ | $2.72\times10^{11}\ M_\odot$ |

  **Grid resolution is preserved, not changed.** `HII_DIM` is derived from
  `target_cell_size_mpc = 2.0` — the resolution of the old 256 Mpc / 128³ grid
  — rather than pinned, so covering the larger footprint does not coarsen
  $M_\mathrm{cell}$. 486.33 / 2.0 = 244 cells, snapped up to 256 for the FFTs.
  This is a **~8× increase in volume and compute** over the previous grid; see
  `TODO.md` for the storage/compute implications.

- **New configuration block** in both entry points, ahead of the grid:
  `SURVEY_AREA_DEG2 = 10.0` (Euclid Deep Field Fornax, RA 03:31:43.6,
  Dec −28:05:18.6), `SURVEY_Z_CENTRAL = 7.0`, and `PHOTOZ_N_SIGMA = 1`.
  `photoz_uncertainty` (σ_z = 0.45, absolute) **moved into this block** from the
  Euclid-parameters section below it, because it now also sets the LOS depth;
  a pointer comment marks its old location. Its value is unchanged.

  The photo-z multiple is an **explicit choice, not a silent default**:
  `PHOTOZ_N_SIGMA = 1` gives $\Delta z = 0.90$ ($z = 6.55$–$7.45$,
  $L_\parallel = 315.6$ Mpc); `= 2` would give $\Delta z = 1.80$
  ($z = 6.10$–$7.90$, $L_\parallel = 634.9$ Mpc).

- **Notebook lightcone range** is now derived: `z_min` / `z_max` come from
  `SIM_BOX`, giving 6.55–7.45 in place of the hardcoded 6.5–7.5.

- **`run_simulation.py` keeps its smoke-test slab** ($z = 6.995$–$7.005$). The
  survey-derived range is exposed as `SURVEY_Z_MIN` / `SURVEY_Z_MAX` and
  reported in the summary as overridden, but is deliberately *not* adopted —
  the existing DO-NOT-WIDEN warning tying it to `TODO.md` P0.1/P0.2 (the
  power-spectrum estimator's LOS-homogeneity assumption) still stands. Only the
  transverse box geometry is footprint-driven there for now.

- **HDF5 provenance attributes** added by `run_simulation.py`:
  `survey_area_deg2`, `survey_z_central`, `survey_delta_z`, `photoz_n_sigma`,
  `survey_z_min`, `survey_z_max`, `survey_los_depth`, and `survey_field`, so a
  stored run records the footprint its box was sized from.

- **`src/conversions.py` imports tidied** — a duplicated `import numpy as np`
  removed; `dataclasses` and `typing` imports added.

<!-- ─── Notebook UV-map deposit, 2026-08-20 ───────────────────────────── -->

### Changed

- **The notebook's UV-map cells now use `analysis.deposit_halo_field`** instead
  of hand-rolling the deposit. §3 of
  `21cmfast_HERAxEuclid_lightcone.ipynb` (cells "Convert halo SFRs to UV
  luminosities" and "Apply Euclid-like UV magnitude selection") previously
  computed `floor(coords / cell_size) % HII_DIM` and filled the grid with
  `np.add.at`, duplicating the helper `src/analysis.py` already provides. The
  binning is now shared with `run_pipeline.py`'s galaxy field, per the
  notebook's own "every formula is imported from `src/`" rule.

  Verified **bit-identical** to the previous path for in-box coordinates. One
  deliberate behavioural change: halos outside `[0, BOX_LEN]` are now *dropped*
  by `histogramdd` rather than periodically wrapped by the modulo, so the cell
  prints a warning counting them instead of silently relocating them.

- **Three write-only full-box allocations removed** from the same cells:
  `luminosity_density_grid`, `selected_luminosity_density_grid`, and
  `selected_galaxy_overdensity` were each computed and never read. At
  `HII_DIM = 256` those are ~134 MB apiece (~400 MB total). Each is replaced by
  a scalar diagnostic print plus a comment giving the one-line form, so the
  quantities remain documented without being materialised.

- **The cubic grid shape is now explained rather than left ambiguous.** These
  maps are built from the *coeval* perturbed halo catalogue
  (`perturb_halo_catalog` on the coeval ICs), so their coordinates span
  `[0, BOX_LEN)` on all three axes and a `(HII_DIM, HII_DIM, HII_DIM)` grid is
  correct — it is **not** the `(HII_DIM, HII_DIM, N_z)` lightcone shape used
  for the power spectra in §4. A comment now says so, and the note in §3
  records why `selected_galaxy_overdensity` must not be fed to the spectra:
  `delta_gal` there comes from the lightcone `halo_sfr` field, and mixing the
  two would be a shape and geometry mismatch.

### Added

- **4 tests in `tests/test_galaxy_weighting.py`** pinning the switch: equality
  with a reference floor/mod implementation for both weighted and unweighted
  deposits, the drop-versus-wrap difference for out-of-box halos, and that the
  cubic coeval grid and the `n_los` lightcone grid stay distinguishable.

<!-- ─── Luminosity-weighted galaxy overdensity, 2026-08-19 ─────────────── -->

### Added

- **`delta_gal` can now be luminosity-weighted instead of number-weighted.**
  `src/analysis.py` gains two functions and one constant:

  | New API | Purpose |
  |---|---|
  | `deposit_halo_field(coords, box_len, n_perp, n_los, los_extent, weights)` | 3D `histogramdd` deposit of a halo catalogue onto the simulation grid, with optional per-halo weights |
  | `galaxy_overdensity_from_catalogue(..., weighting=...)` | Euclid-selected `delta_gal` in either mode; returns `(field, EuclidSelection)` |
  | `GALAXY_WEIGHTING_MODES` | `("number", "luminosity")` |

  The two modes differ only in the per-halo weight, and are normalised
  identically:

  | Mode | Formula |
  |---|---|
  | `"number"` (default) | $\delta_\mathrm{gal} = N / \langle N \rangle - 1$ |
  | `"luminosity"` | $\delta_{\mathrm{gal},L} = \sum L_\mathrm{UV} / \langle \sum L_\mathrm{UV} \rangle - 1$ |

  `L_UV` comes from the existing `conversions.sfr_to_Luv()`
  ($L_\mathrm{UV} = \mathrm{SFR} / \kappa_\mathrm{UV}$, $\kappa_\mathrm{UV}
  = 1.15\times10^{-28}$, Madau & Dickinson 2014) — consumed, not
  reimplemented. Both modes return the same shape and the same zero mean, so
  they are interchangeable with no downstream change.

- **`GALAXY_WEIGHTING` flag in `run_simulation.py`** (§3b), with three values:

  | Value | Source of `delta_gal` |
  |---|---|
  | `"lightcone_sfr"` (default) | the lightcone `halo_sfr` field — **existing behaviour, unchanged** |
  | `"number"` | Euclid-selected catalogue, unit weights |
  | `"luminosity"` | Euclid-selected catalogue, `L_UV` weights |

  The selected array flows into the identical downstream path — Kaiser RSD,
  HDF5 `galaxy_overdensity`, and `compute_all_power_spectra()`. The mode is
  recorded as the root attribute `galaxy_weighting`.

- **`tests/test_galaxy_weighting.py`** — 17 tests covering weight
  conservation, non-cubic grids, the zero-mean normalisation, both formulas
  against manual recomputation, scale-invariance of the `kappa_UV` factor,
  that the two modes genuinely differ, and the error paths (bad shapes,
  unknown mode, empty selection).

### Changed

- **`figures.plot_uv_selection_maps()` now calls `deposit_halo_field()`**
  instead of its own `np.histogram2d` pair. The projected maps are the LOS
  sum of the same 3D deposit that builds `delta_gal`, so figure and field bin
  identical positions.

- **`M_UV_bright` / `M_UV_faint` moved up** in `run_simulation.py` to the
  Euclid survey-parameter block, since §3b now selects on them as well as §4.
  Values are unchanged (`-22`, `M_UV_limit`).

### Fixed

- **`plot_uv_selection_maps()` indexed the full catalogue with a mask sized to
  the valid subset.** `selection.mask` from `select_euclid_halos()` is defined
  over halos with `SFR > 0` and `M > 0`, not over the whole catalogue, so
  `coords[selection.mask]` and `magnitude[selection.mask]` raised
  `IndexError` — or silently mis-aligned — as soon as any halo had
  `SFR <= 0`. The mask is now lifted back to full-catalogue indices.

### Notes

- The pre-existing default `delta_gal` was **never** $N/\bar N - 1$: it is the
  lightcone `halo_sfr` density, `sfr_field / mean_sfr - 1`. Since
  $L_\mathrm{UV} \propto \mathrm{SFR}$, that field is *already*
  luminosity-weighted up to the constant $\kappa_\mathrm{UV}$, which divides
  out of the ratio. The number-vs-luminosity distinction therefore only
  exists on the catalogue path.
- The catalogue is **coeval** at `z_obs` and spans `BOX_LEN` along the LOS,
  while the lightcone spans `L_los`. The catalogue modes are deposited into
  an `(HII_DIM, HII_DIM, N_z)` grid for shape compatibility, but they carry
  no redshift evolution along the LOS, and their LOS cell size is
  `BOX_LEN / N_z`. See the open `L_los` discrepancy in `TODO.md`.

<!-- ─── Notebook/pipeline figure parity, 2026-08-18 ────────────────────── -->

### Added

- **The HPC pipeline gained the four figures that existed only in
  `21cmfast_HERAxEuclid_lightcone.ipynb`.** `src/figures.py` goes from 11 to
  15 figures, all wired into `run_pipeline.py`'s existing plot groups:

  | New figure | Group | Notebook origin |
  |---|---|---|
  | `photoz_suppression` | `budget` | §7b — `W(k_par)` swept over σ_z |
  | `galaxy_wedge` | `power` | §7c — `P_gal` with the wedge filled, not outlined |
  | `wedge_real_space` | `power` | §7d — wedge excision applied to the 3D field |
  | `uv_selection_maps` | `scaling` | §3 — projected UV luminosity and Euclid-selected counts |

- **Each reuses `src/analysis.py` rather than restating it.**
  `plot_galaxy_wedge` reuses `_add_wedge_lines` / `_style_k_axes`;
  `plot_wedge_real_space` reuses `foreground_wedge_mask`, reshaped onto the 3D
  FFT grid (the wedge condition factorises, so the flattened transverse plane
  gives a mask that reshapes straight back); `plot_uv_selection_maps` reuses
  `select_euclid_halos`, so the map and the bias stage select the same halos;
  `plot_photoz_suppression` reuses `photoz_damping_kernel` and obtains each
  scenario's σ_r by scaling the budget's own `radial_smearing` — exact, since
  σ_r = c σ_z / H(z) is linear in σ_z, and it cannot drift from the adopted
  value.

- `figure_stage()` now takes an explicit `m_uv_bright` argument, threaded from
  `--m-uv-bright` at the call site, so the selection-map figure cuts where the
  bias stage cuts.

### Fixed

- **Three pipeline tests asserted the exact per-group figure inventory** and so
  failed the moment a group gained a figure — which is the behaviour those
  assertions are for. `test_plot_selection_limits_output`,
  `test_budget_figure_is_written` and `test_pdf_output_format` now list the new
  expected sets (the last also needed `sorted()`, having compared an unordered
  `os.listdir` against a one-element list).

### Verified

- **Six new tests in `tests/test_figures.py`**, going past "it renders":
  the wedge region is a hatched patch and not just lines; the real-space
  filtered panel has strictly lower variance than the original and shares its
  colour limits; a non-zero wedge buffer excludes at least as much as the bare
  line; the σ_z = 0 photo-z curve is flat at W = 1 and the family never
  crosses; the adopted σ_z is drawn even when a caller omits it; the scaled
  σ_r reproduces `radial_smearing_length` exactly; and the selection map's
  counts sum to `select_euclid_halos`'s own `n_selected`.
- Suite **99 → 107 passing**.
- **Full end-to-end run against the real 2.7 GB `outputs/lightcone_data.h5`**
  (`--plots all`, into a scratch figure directory so the committed outputs were
  left alone): all **15 figures written**, no errors. The selection maps
  recover 49,315 of 114,291,212 halos with visible cosmic web in both panels.

### Parity gaps found but deliberately not closed

- **The bright-end magnitude cut differs.** The notebook uses
  `M_UV_bright = -22.66` ("obtained from collaborators"); `run_simulation.py`
  line 492 hardcodes `-22` and `--m-uv-bright` defaults to `-22.0`. The
  *capability* is at parity — the flag exists — but the default is not, and
  changing it moves the pipeline's selection, b_g and n_gal. Left as a science
  decision; `--m-uv-bright -22.66` reproduces the notebook today.
- **`sfr_mini` / `fesc_sfr` panels** (notebook §3, Plots 4–5) have no pipeline
  equivalent. Both fields are `None` under the `"simple"` template, so the
  panels never render in practice; building pipeline versions would add code
  that no current configuration exercises.
- **Not a code gap, but worth knowing:** the cached
  `outputs/lightcone_data.h5` is a *thin slab* — z = 6.995–7.005, 100 LOS
  cells spanning ~3.2 Mpc — where the notebook runs z = 6.5–7.5 over 350.8 Mpc
  with 175 cells. `wedge_real_space` is correct on it but uninformative, since
  the box holds almost no line-of-sight structure to remove. Re-running the
  simulation over the notebook's redshift range would make that figure say
  something.

<!-- ─── Real-space wedge figure §7d, 2026-08-18 ────────────────────────── -->

### Added

- **§7d of `21cmfast_HERAxEuclid_lightcone.ipynb`: the foreground wedge shown in
  real space.** §7c draws the wedge as a boundary in $(k_\perp, k_\parallel)$;
  §7d applies it — FFT the 3D galaxy overdensity, zero every mode with
  $k_\parallel \le m_{\rm horizon}k_\perp$, inverse-FFT, and show the same LOS
  slice before and after on a shared colour scale.

- **Nothing is recomputed.** The field is `galaxy_overdensity` **post-Kaiser**
  (§4 rebinds it), i.e. exactly what §6 fed to `compute_all_power_spectra`, so
  the panel really is `P_galaxy_auto`'s underlying field. The wavenumber grids
  are §4's `KX/KY/KZ`, the slope is §7's `horizon_slope`, and the slice index
  and extent are §5's `mid_y` / `extent_los` with §5's orientation convention
  (`field[:, mid_y, :]`, no `.T`, `origin="lower"` → LOS on x).

- **The mask is `src.analysis.foreground_wedge_mask`, not a hand-rolled
  comparison.** The wedge condition factorises — $k_\perp$ varies over the
  transverse plane only, $k_\parallel$ over the LOS axis only — so the
  $(N^2, N_z)$ mask the function returns for the flattened transverse plane
  reshapes straight onto the $(N, N, N_z)$ grid. Verified identical to the
  direct comparison `|KZ| > hypot(KX, KY) * horizon_slope`.

- **`buffer=0.0` is passed deliberately**, giving the bare horizon line
  $k_\parallel \le m k_\perp$ that §7c draws. §8's excision adds `wedge_buffer`
  on top, so the budget discards slightly *more* than this figure shows.

### Verified

- Executed at the **fiducial grid geometry** — $(128, 128, 175)$, 2.0 Mpc
  cells, $L_{\rm LOS} = 350.8$ Mpc, read off the run's own §1 output — in 0.2 s.
- **97.4 % of the 3D modes fall inside the bare horizon wedge**
  (2,791,938 of 2,867,200). This number is purely geometric — it depends only
  on the $k$-grid and the slope, not on the field — so it is the figure's real
  value, not an artefact of the stand-in field used to exercise the cell.
- The filtered field is real-valued to $10^{-8}$ after the round trip, and its
  variance is strictly reduced. **The retained-variance figure printed by the
  cell is field-dependent** and was only exercised against a power-law mock, so
  no value for it is quoted here; the real one appears when the notebook is run.
- A first draft of the §7d prose had the physics backwards — it claimed the
  filtered field "keeps its transverse texture but loses coherence along the
  LOS". Rendering the figure showed the opposite: the cut keeps
  $k_\parallel > 3.15\,k_\perp$, so what survives varies *rapidly along the LOS*
  and is smooth transversally, which is why the filtered panel is striped across
  the LOS axis. The markdown was corrected before the cell was finalised.

### Confirmed from the user's own run

The notebook on disk now carries execution outputs through `In[28]`, which
covers §7b and §7c added earlier today. Both **ran clean on real data**, no
errors. Notably $1/\sigma_r = 0.0064$ Mpc⁻¹ does sit below the real lowest
sampled $k_\parallel = 0.0102$ Mpc⁻¹, so §7b's axis extension does engage as
designed; $W$ at that first bin is **0.274**, less severe than the coarser
stand-in grid had suggested. `galaxy_bias = 5.394` in the live namespace,
matching the value documented in the README fix above.

<!-- ─── Optional halo-catalogue fields TypeError, 2026-08-18 ───────────── -->

### Fixed

- **`21cmfast_HERAxEuclid_lightcone.ipynb` §3's halo-catalogue plotting cell
  crashed on `np.isfinite(sfr_mini)`** with

  ```
  TypeError: ufunc 'isfinite' not supported for the input types, and the inputs
  could not be safely coerced to any supported types according to the casting
  rule ''safe''
  ```

  whenever the run had no mini-halos — which is *every* run under the `"simple"`
  template, i.e. the notebook's default. Plot 4 was unreachable.

  **Cause.** py21cmfast v4 declares optional catalogue fields as
  `_arrayfield(optional=True)` → `attrs.field(default=None)`
  (`py21cmfast/wrapper/outputs.py:888`). The attribute therefore **always
  exists** and is `None` when unpopulated, so the `hasattr(perturbed_halos,
  "sfr_mini")` guard was unconditionally True. `get_21cmfast_array` then
  returned `np.asarray(None)` — a 0-d array of dtype **object**, which is *not*
  `None` and so passed the `sfr_mini is not None` check, but which every ufunc
  rejects. Three steps, each individually reasonable, composing into a crash.

  **Fix**, both in the cell that owns the helper:

  | Change | Effect |
  |---|---|
  | `get_21cmfast_array` returns `None` for a `None` input | no more 0-d object arrays; the `is not None` guards downstream become meaningful |
  | guards use `getattr(perturbed_halos, ..., None)` and test the **value** | `hasattr` was testing the wrong thing for an attrs field |

  `fesc_sfr` (Plot 5) carried the identical latent bug and is fixed by the same
  two changes.

- **Minihalo SFR was left in the wrong unit.** Plot 4 histograms `sfr_mini` on
  the *same axes* as `sfr`, but only `sfr` got the M_sun s⁻¹ → M_sun yr⁻¹
  correction documented in `docs/Low_SFR_fix.md`. Once the TypeError was fixed
  and Plot 4 could actually run, the two distributions would have sat **7.5 dex
  apart on an axis labelled M_sun/yr**. `sfr_mini` now takes the same
  `* _SEC_PER_YR` conversion.

  `fesc_sfr` is very likely in the same internal unit, but it is plotted on its
  own axes with no unit in the label, so it was **left as returned** rather than
  converted on an assumption — flagged in a comment at the point of use.

### Verified

- The cell **executed both ways** against a stand-in catalogue: with
  `sfr_mini = None` (the reported failure — now runs clean, Plots 4–5 skipped,
  3 figures) and with `sfr_mini` populated (Plot 4 runs, 4 figures).
- Units confirmed consistent after the fix: median `sfr_mini` = 9.3×10⁻⁴ vs
  median `sfr` = 1.0×10⁻¹ M_sun yr⁻¹ — the ~2 dex minihalo offset that is
  physically expected, not the 7.5 dex unit artefact.
- The failure mode itself was reproduced directly:
  `np.asarray(None)` → `dtype=object`, `ndim=0`, `is not None` → True,
  and `np.isfinite` on it raises the reported error verbatim.
- Suite still **99 passing**; notebook valid under `nbformat` (37 cells).
- A new print reports which optional fields were populated, so the next run
  says *why* Plots 4–5 are absent instead of silently omitting them.

<!-- ─── Notebook diagnostic figures §7b/§7c, 2026-08-18 ─────────────────── -->

### Added

- **Two diagnostic figures in `21cmfast_HERAxEuclid_lightcone.ipynb`**, between
  §7 (the 2D spectra) and §8 (the uncertainty budget). Both are pure additions:
  the notebook diff is **194 inserted lines, 0 deleted**, and no pre-existing
  cell was touched or reformatted.

  | § | Figure | What it shows |
  |---|--------|---------------|
  | **7b** | Photo-$z$ suppression | $W(k_\parallel) = e^{-k_\parallel^2\sigma_r^2/2}$ for $\sigma_z \in \{0,\,0.02,\,0.05,\,0.10,\,0.30,\,0.45\}$ on log-$k_\parallel$ axes, each labelled with its $\sigma_r$, the adopted Euclid Wide case emphasised, $1/\sigma_r$ marked |
  | **7c** | $P_\mathrm{gal}$ vs. the wedge | $\log_{10}\lvert P_\mathrm{gal}\rvert$ with the wedge region **filled and hatched**, not merely outlined, plus the horizon (solid) and HERA FoV (dashed) boundaries |

- **Neither cell reimplements anything.** §7b imports
  `radial_smearing_length` and `photoz_damping_kernel` from `src.analysis` —
  the same two functions `compute_uncertainty_budget()` applies in step ① —
  and passes `**COSMOLOGY`, so the figure cannot drift from the damping
  actually adopted. §7c reuses §7's `fill_nan_nearest`, `k_perp`, `k_parallel`,
  `P_galaxy_auto`, `horizon_slope`, `fov_wedge_slope_value` and the
  `k_perp_line` / `k_par_horizon` / `k_par_fov` lines built from them. The only
  new code in either cell is the figure itself.

### Two judgement calls worth recording

- **§7b's x-axis is extended below `k_parallel[0]` when — and only when —
  $1/\sigma_r$ falls outside the sampled range.** At the adopted
  $\sigma_z = 0.45$, $\sigma_r = 157.5$ Mpc gives $1/\sigma_r = 0.0064$
  Mpc⁻¹, while the box's lowest $k_\parallel$ bin sits at $\approx 0.018$
  Mpc⁻¹. Held strictly to `k_parallel`'s range the requested marker renders
  off-canvas, so the axis extends to $0.7/\sigma_r$ and the unsampled strip is
  shaded with its own legend entry. **This is the physical content of the
  figure, not a cosmetic fix**: every mode the lightcone measures already lies
  past the damping scale, and $W = 0.018$ at the very first bin. `x_lo` is the
  single line to change to pin the axis back to `k_parallel`.
- **§7c's wedge overlay is dark (35 % black + white hatching), not light.** A
  white overlay was tried first: on `plasma` it *brightens* the already-bright
  low-$k$ corner and the legend swatch disappears against the colourbar's
  yellow end. The dark overlay dims the excluded modes so the accessible
  upper-left window reads as the signal region at a glance.

### Verified

- Both cells **executed standalone against the real `src/` functions** in the
  `21cmfast` env and their output figures inspected.
- §7b reproduces the documented smearing lengths exactly:
  $\sigma_z = 0.45 \Rightarrow \sigma_r = 157.5$ Mpc, matching
  `docs/uncertainty_budget.md` §4.2 and the notebook's own §8 prose; the
  $\sigma_z = 0$ curve is flat at $W = 1$, as it must be.
- One bug caught and fixed during that check: the per-scenario print reported
  $W$ at the *extended plotting floor* rather than at the lowest sampled bin
  (0.783 instead of 0.0184 for $\sigma_z = 0.45$). It now evaluates
  `photoz_damping_kernel(k_parallel[:1], sigma_r)` explicitly.
- §7c exercised with a NaN bin present, confirming `fill_nan_nearest` is
  reached; horizon and FoV slopes 3.151 and 0.379 at $z = 7$.
- Suite still **99 passing**; notebook valid under `nbformat` (37 cells); all
  code cells parse; all cell ids unique.

### Fixed

- **`README.md`'s notebook-3 fiducial table still listed $\sigma_z = 0.059$**,
  directly contradicting the configuration cell (0.45 since the
  absolute-vs-fractional correction) and the new §7b figure built on it. Now
  0.45, with the absolute-vs-fractional distinction spelled out inline. The
  identical row in the **archived** coeval notebook's table was deliberately
  left at 0.059 — that notebook was run with 0.059 and the table is a record of
  what was run, so it gained only a clarifying note.
- **`README.md`'s §-structure list for notebook 3** gained §7b and §7c, and its
  §8 entry — still described as "photo-$z$ damping and foreground wedge
  excision" — now matches the notebook's actual header, the single
  `compute_uncertainty_budget()` call.
- **`README.md` still advertised a fixed $b_\mathrm{gal} = 8$ for notebook 3.**
  That fallback no longer exists — `galaxy_bias = 8` is **commented out** in the
  configuration cell. The row now states that the bias is computed in-line and
  says where: **§3**, the analytic-galaxy-bias cell, at
  `galaxy_bias = galaxy_bias_hmf` — a Simpson integral of the Sheth-Tormen halo
  bias weighted by the HMF over the Euclid-selected mass bins,

  $$b_g = \frac{\int b_h(M)\,\frac{dn}{d\log M}\,d\log M}{\int \frac{dn}{d\log M}\,d\log M} \approx 5.39,$$

  consumed downstream by §4's $\beta = f/b$, and requiring `hmf`. The quoted
  5.39 is the value from the verified clean-namespace run recorded in the
  2026-08-17 entry above, not a fresh execution — reproducing it needs the full
  upstream lightcone and halo catalogue.

  The $b_\mathrm{gal} = 8$ rows for the **archived coeval notebook** and the
  **HPC pipeline** were left alone: the first records what that notebook was run
  with, and the second already documents itself correctly as a fallback
  overwritten by the halo-catalogue estimate.

### Not covered by a new test

`tests/` is unchanged. Both additions are notebook plotting cells that define
no new function — every formula they draw is already covered by
`tests/test_uncertainty_budget.py`. There is no notebook-execution test in the
suite for either of them to hook into.

<!-- ─── Notebook galaxy-bias cell debugged, 2026-08-17 ──────────────────── -->

### Fixed

- **The analytic galaxy-bias cell of `21cmfast_HERAxEuclid_lightcone.ipynb`
  carried four bugs.** It now mirrors `run_simulation.py` §4's
  "Fallback / cross-check: analytic HMF integral" and imports every
  conversion from `src/`. **`b_g` 33.66 → 5.39**, which is where high-$z$ LBGs
  at $M_\mathrm{UV} < -18$ should sit.

  | # | Bug | Effect |
  |---|---|---|
  | 1 | **`Muv_to_Luv` was used but never defined or imported anywhere in the notebook** | `NameError` on any clean kernel — the cell could only ever have run against a name left over from a deleted cell or an earlier session |
  | 2 | `sfr()` divided stellar mass by a **hardcoded 1e8 yr** | 21cmFAST's own timescale is $t_\mathrm{STAR}\,t_H(z) = 570.3$ Myr at $z = 7$, so every SFR was **5.70× too high** — 1.89 mag too bright. Exactly the failure `star_formation_timescale`'s docstring warns about |
  | 3 | `uv_luminosity()` computed `1.15e28 * sfr` | $\kappa_\mathrm{UV} = 1.15\times10^{-28}$ enters as $L = \mathrm{SFR}/\kappa$, i.e. $\times 8.696\times10^{27}$. **1.32× too luminous** |
  | 4 | local `sheth_tormen_bias_from_nu` used `a * nu**2` on `mf.nu` | `hmf` stores `nu` **already squared**, so this formed $\nu^4$. $b_h(\nu^2{=}10)$ came out **42.4 instead of 4.73** — the dominant cause of the implausible $b_g \approx 33$ |

  Bugs 2–4 all pushed in the same direction, which is why the error was large
  but the cell never looked obviously broken.

- **Three name collisions removed.** The cell bound `sfr`, `selected`,
  `stellar_mass` and `L_UV` — all of which the UV-selection cells below use for
  **per-halo catalogue arrays**. Only top-to-bottom execution order was hiding
  it; re-running this cell after them (what one does while tuning the $M_UV$
  limits) would silently replace a 114-million-element array with a
  300-element one. Renamed to `sfr_model`, `selected_mass_bins`,
  `stellar_mass_model`, `L_UV_model`, following `run_simulation.py`.

- **Caveats now printed rather than implied**: `hmf`'s $M_\odot/h$ and
  $h^3$ Mpc⁻³ conventions against a stellar-halo relation written in $M_\odot$;
  `MassFunction()` using `hmf`'s default cosmology rather than the notebook's;
  this being the scatter-free estimate that the pipeline does *not* prefer; and
  the stored `galaxy_bias = 33.389` attribute in `lightcone_data.h5`, which
  still correctly describes the stored `galaxy_overdensity` field and so was
  deliberately left unpatched.

### Verified

- The rewritten cell **executes from a clean namespace** (previously
  impossible, bug 1) and reproduces the pipeline's analytic branch exactly:
  $b_g = 5.3944$, $n_\mathrm{gal} = 2.394\times10^{-3}\ h^3$ Mpc⁻³, 62 of 300
  mass bins selected, $M_h \in [3.31\times10^{10}, 5.50\times10^{11}]\ M_\odot/h$.
- It leaves `L_UV`, `sfr`, `selected` and `stellar_mass` unbound, confirmed by
  introspection after execution.
- Suite still **99 passing**; the budget chain still bit-identical to
  `run_pipeline.py`; notebook valid under `nbformat`.

### Consequence not yet propagated

`beta_rsd = f / b_g` and the stored `galaxy_overdensity` field were built with
$b_g = 33.389$. With $b_g = 5.39$, $\beta$ rises from 0.0299 to **0.185** and
the maximum Kaiser boost from 1.03× to **1.19×**. Propagating this needs
`--sim force`; until then `lightcone_data.h5` and every number derived from it
still carry the old bias. Not yet tracked in `TODO.md`.

<!-- ─── Notebook aligned with the pipeline, 2026-08-17 ──────────────────── -->

### Changed

- **`21cmfast_HERAxEuclid_lightcone.ipynb` now imports the pipeline's code
  instead of duplicating it.** Seven cells rewritten; roughly 190 lines of
  reimplemented physics deleted. The notebook and `run_pipeline.py` now execute
  the same functions, so they cannot drift apart.

  | Notebook cell | Was | Now |
  |---|---|---|
  | imports | `sfr_to_Luv` only | + `hubble_parameter`, `comoving_distance`, `compute_all_power_spectra`, `compute_uncertainty_budget`, `horizon_wedge_slope`, `fov_wedge_slope` |
  | derived quantities | local `hubble_parameter` def | imported, called with the configured cosmology explicitly |
  | power spectra | a ~100-line local copy of `compute_cylindrical_cross_power` | `compute_all_power_spectra(...)` |
  | wedge geometry | longhand slope through $\lambda_\mathrm{obs}$ | `horizon_wedge_slope` / `fov_wedge_slope` |
  | photo-$z$ + wedge | inline kernel, mask, `quad` call | one `compute_uncertainty_budget(...)` call |
  | SNR map | inline $T_\mathrm{sys}$, $P_N$, $\sigma$, SNR | read off the returned `UncertaintyBudget` |

- **Three notebook bugs fixed**, all found by the audit recorded below:
  - `Hz_obs = Planck18.H(z).value` referenced **`z`, a stale loop variable**,
    where `z_obs` was meant. Now `hubble_parameter(z_obs, HUBBLE_CONSTANT, OMEGA_M_0)`.
    (It fed only a figure overlay, so no published number depended on it.)
  - The noise cell hardcoded `1.42e9` while the notebook's own config defines
    `F_21_HZ = 1420.405e6` — a 0.10 % shift in $P_{N,21}$. The literal is gone.
  - `wedge_buffer` 0.02 → **0.0677 Mpc⁻¹** and `photoz_uncertainty`
    0.059 → **0.45**, matching the pipeline, each with its provenance in a
    comment.

- **Cosmology reverted to literal `67.36` / `0.315`** in the configuration
  cell, matching `run_simulation.py:136-137` and hence the `HUBBLE_CONSTANT` /
  `OMEGA_M_0` attributes of `lightcone_data.h5` that `run_pipeline.py` reads.
  The uncommitted `Planck18.Om0` / `Planck18.H0.value` edit (67.66 / 0.30966)
  is undone; `Planck18` is retained for the comoving-distance **endpoints**
  only, exactly as `run_simulation.py` §1 does. **Side effect:** because this
  restores the cosmology the notebook's last real execution used, the stored
  outputs of every *untouched* cell remain valid.

- **Markdown updated** to match: the photo-$z$/wedge section is now a full
  four-step derivation of the budget including the $T_0(z)$ cancellation
  algebra and the sample-variance / noise-coupling split; the power-spectrum
  section states the estimator and flags the unused `mode_counts`; the summary
  gained a table mapping each concern to the `src/` function that implements
  it, and a table of every parameter that changed with the reason.

### Verified

- **The rewritten cells execute, and agree with the pipeline exactly.** Cells
  2, 4, 6, 24, 26, 28 and 30 were run headlessly against
  `outputs/lightcone_data.h5` (the 21cmFAST cells cannot run outside a full
  simulation, so the fields were injected) and compared to
  `run_pipeline.observational_stage` on the same data:
  $\sigma_r$, horizon slope, FoV slope, $P_{N,21}$, $P_{N,\mathrm{gal}}$, the
  surviving mode count, the per-mode SNR map, and the total
  ($1.0574485217836499\times10^{-111}$) are all **bit-identical**.
- **`src/figures.py` deliberately not imported** by the notebook: it calls
  `matplotlib.use("Agg")` at import, which would silently disable inline
  plotting. `fill_nan_nearest` therefore stays defined locally — the one piece
  of duplication left, and it is display-only.

### Note for the next run

Outputs of the seven rewritten cells were **cleared**, not regenerated —
producing them requires a full 21cmFAST lightcone run. Re-execute the notebook
top to bottom to repopulate them. The numbers it will print for the budget are
already known and recorded in `docs/uncertainty_budget.md` §8, since they do
not depend on the simulated fields.

<!-- ─── Uncertainty budget incorporated into the pipeline, 2026-08-17 ───── -->

### Verified (no change required)

- **The pipeline's transcription of the notebook's uncertainty budget is
  numerically exact.** Every term of the photo-$z$ / wedge / noise / SNR chain
  in `21cmfast_HERAxEuclid_lightcone.ipynb` (its "Photo-z Radial Smearing &
  Foreground Wedge Mask" and "Cross-correlation SNR map" cells) was
  transcribed verbatim and compared against `src/analysis.py` on the
  notebook's own grid (256 × 350.6 Mpc, 128² × 175, 20 × 20 log bins).
  **Bit-identical**: $\sigma_r$, $W(k_\parallel)$, the wedge mask (identical
  boolean array), $P_\times^\mathrm{obs}$, $P_\mathrm{gal}^\mathrm{obs}$,
  $\sigma_{21}$, $\sigma_\mathrm{gal}$, $\sigma_\times$, the per-mode SNR, and
  the total. The horizon and FoV slopes agree to $2\times10^{-16}$ relative —
  the notebook writes the slope the long way round through
  $\lambda_\mathrm{obs}$, which reduces algebraically to $D_c H/[c(1+z)]$.
  Now locked by `test_budget_reproduces_the_notebook_chain`.

- **The missing $T_0(z)$ factors are correct to omit.** La Plante Eqs. 15–17
  carry $\sigma_{21} = (P_{21}+P_{N,21})/T_0^2$ and
  $\mathrm{SNR} = |P_\times/T_0|/\sigma_\times$; neither the notebook nor the
  pipeline carries them. They **cancel exactly** in the ratio, so the SNR is
  unaffected. They do *not* cancel in $\sigma_\times$ itself, so it must not be
  quoted as a standalone error bar without reinstating $T_0$ — recorded in the
  `cross_power_snr` docstring and `docs/uncertainty_budget.md` §2.5.

- **The notebook's stored outputs predate its own configuration cell.** Its
  printed $D_c = 8821$ Mpc, $m = 3.151$, $\sigma_r = 20.6$ Mpc reproduce
  exactly under $H_0 = 67.36$, $\Omega_{m,0} = 0.315$ — but its (uncommitted)
  config cell now reads `Planck18.Om0` / `Planck18.H0.value` (0.30966 / 67.66),
  which would give 8835 Mpc, 3.1429, 20.73 Mpc. **The notebook has not been
  re-executed since that edit.** Its `Hz_obs` cell also carries a live bug —
  `Planck18.H(z).value` references `z`, a stale loop variable, where `z_obs` is
  meant; that value feeds only a figure overlay, and the wedge-mask cell
  recomputes $H(z_\mathrm{obs})$ correctly, so no science number depends on it.
  The pipeline reads cosmology from the HDF5 root attributes and is therefore
  always self-consistent with the run it describes.

### Added

- **`src/analysis.py` — `compute_uncertainty_budget` and `UncertaintyBudget`.**
  A single entry point for the whole chain (damping → wedge → noise → variance
  → SNR), returning a dataclass that holds every intermediate term plus derived
  summaries (`total_snr`, `fraction_outside_wedge`, `cosmic_variance_fraction`,
  `detected`, `as_dict()`). The notebook and the HPC run now execute the same
  code path.
- **The two halves of Eq. 15 are now separated.** `SNRResult` gained
  `sigma_21cm`, `sigma_galaxy`, `cosmic_variance_term`
  ($\tfrac12 P_\times^2$) and `noise_coupling_term`
  ($\tfrac12\sigma_{21}\sigma_\mathrm{gal}$) — the split that makes this a
  *budget* rather than one error bar, and the basis of the new
  `cosmic_variance_fraction` diagnostic. **For the stored run it is
  $2\times10^{-224}$**: the measurement is entirely noise dominated, because
  photo-$z$ damping has erased $P_\times$ wherever the wedge admits it.
- **`analysis.system_temperature(z, f_21_hz)`** returning $(T_\mathrm{sys}, \nu_\mathrm{obs})$,
  with the model's constants named — `T_RECEIVER_K` (100 K),
  `T_SKY_300MHZ_K` (60 K), `SKY_SPECTRAL_INDEX` (2.55), and
  `NOISE_NORMALISATION_MPC3` ($10^3$). The last is **not** a physical
  constant: it supplies the Mpc³ that $T_\mathrm{sys}^2/(t\Delta\nu)$ lacks,
  standing in for the per-mode survey volume La Plante Eq. 11 computes
  properly. It was previously an unnamed inline literal in both notebook and
  pipeline.
- **Four CLI overrides** on `run_pipeline.py`: `--sigma-z`, `--wedge-buffer`,
  `--integration-time`, `--bandwidth`. Each resolves CLI → HDF5 attribute →
  default. None affects the simulated fields, so all can be swept from cached
  spectra in seconds. `python run_pipeline.py --sigma-z 0.059 --wedge-buffer 0.02`
  reproduces the notebook's configuration without editing anything.
- **The budget is now persisted**, not discarded: `save_uncertainty_budget` /
  `load_uncertainty_budget` write an `uncertainty_budget` group into
  `outputs/analysis_products.h5` (10 maps + 21 scalar attrs), appended
  alongside the cached spectra and replaced in place on recomputation.
- **`figures.plot_uncertainty_budget`** and the `budget` plot group: the
  damping kernel against the lowest $k_\parallel$ the wedge admits, the
  $\sigma_\times$ map, and the sample-variance share per mode. It shows §7.1
  directly — the kernel is below $10^{-26}$ before the wedge floor is reached.
- **`docs/uncertainty_budget.md`** — the full reference: every formula with its
  provenance and evaluated number at $z = 7$, the audit table, the three
  discrepancies, the parameter/override table, the HDF5 schema, and five known
  limitations.
- **23 tests** (16 in `tests/test_uncertainty_budget.py`, plus CLI and figure
  coverage). Suite: **76 → 99 passing**.

### Changed

- **`run_pipeline.observational_stage` now contains no physics.** It resolves
  parameters and calls `compute_uncertainty_budget` once; the inline
  transcription it used to hold is gone. Verified behaviour-preserving against
  the stored run: total SNR `1.0574485217836499e-111`, $\sigma_r$ 157.478 Mpc,
  horizon slope 3.150906 — identical to the values in the previous
  `pipeline_summary.json`.
- **`pipeline_summary.json`** gained the `uncertainty_budget` block (21
  scalars). The former `observation` block is **retained as an alias** of its
  eight values, so existing notes and scripts keep working.
- **`print_report`** now prints the budget as four labelled lines — photo-$z$,
  wedge, noise, variance split — instead of two bare numbers.
- **`docs/HPC.md`** — §3.1 gained the four flags, §5.2 a pointer to the new
  document, §5.5 the variance split and the $T_0$ cancellation, §6/§7/§8 the
  new figure, JSON block, and file sizes; §11.4 now records that the notebook
  omits the mode-count factor too (so matching it is faithful, not a porting
  error), and §11.5 names the $10^3$ normalisation.
- **`README.md`, `PIPELINE.md`** — documentation table, output inventory, and
  the `analysis_products.h5` schema updated.

### Known, unchanged

- `mode_counts` is still not divided into the variance. The **notebook does
  not apply it either**, so the pipeline matching it is a faithful port; the
  quoted total SNR is a per-bin quadrature sum and therefore conservative.
  (`TODO.md` §P1.1.)
- The thermal noise is still $k$-independent. La Plante Eq. 11 with
  $X^2Y\Omega'/n(k_\perp)$ exists in `21cm_galaxy_cross_uncertainty.ipynb` and
  was deliberately **not** ported — out of scope for the HERAxEuclid notebook's
  budget. (`TODO.md` §P1.2.)
- The notebook's noise cell hardcodes `1.42e9` where its own config defines
  `F_21_HZ = 1420.405e6`; the pipeline uses the precise value throughout,
  moving $P_{N,21}$ by 0.10 % (3.752581 → 3.748786). Recorded in the
  `hera_thermal_noise_power` docstring; not propagated.
- `mean_galaxy_density` is still declared `h³ Mpc⁻³` and consumed as
  $1/\bar n$ in Mpc³ — carried over from the notebook unchanged.

<!-- ─── User-defined parameter requirements, 2026-08-13 ─────────────────── -->

### Added
- **`docs/HPC.md` §13 — user-defined parameter requirements.** A single
  checklist of every parameter the user has to (or may) set for an HPC run,
  **each with the file and line where it is set**:
  - **§13.0** the four layers parameters live in (scheduler directives →
    `submit_job.sh:26–33` → `run_pipeline.py` CLI → `run_simulation.py:82–170`)
    and the fact that a value set in one layer does not propagate backwards.
  - **§13.1** seven site-specific requirements (R1–R7) with no correct default:
    the **absent `#SBATCH` directives** (`submit_job.sh`, insert after line 24),
    `CONDA_ENV`, the **21cmFAST cache directory** (R3, see below), the
    cwd-relative `OUTPUT_DIR`, the hardcoded `EMAIL_TO`, install-time
    `PIP_CACHE_DIR`/`XDG_CACHE_HOME`, and `JOB_NAME`.
  - **§13.2** the full config block by line number, plus the parameters that
    are user-editable but sit **outside** it: `minimum_los_slices` (line 216),
    `M_UV_bright` (line 492), `OMEGA_B_0` (line 558), `random_seed` (255),
    the template name (253), `apply_rsds`/`include_dvdr_in_tau21` (278–279),
    and `N_THREADS` (never set — one core regardless of what the scheduler
    is asked for).
  - **§13.3** the four CLI flags that matter on a cluster.
  - **§13.4** which edits actually take effect and when — group A
    (14 parameters read from the **HDF5 root attrs**, so a config-block edit
    does nothing until `--sim force` or an in-place patch), group B (baked into
    the stored `galaxy_overdensity`; re-simulation only), group C (live).
  - **§13.5** eight consistency rules the code does not enforce, **§13.6** a
    pre-flight checklist.
  - Former §13 (References) renumbered to **§14**; no internal cross-reference
    pointed at either number.
### Verified (no change required)
- **Where the 56 GB 21cmFAST cache actually comes from** (py21cmfast 4.1.1,
  checked by introspection). `p21c.run_lightcone(**kwargs)` forwards to
  `generate_lightcone`, whose `cache` argument defaults to
  **`OutputCache(direc=Path('.'))`** — the *current working directory*, which
  is why the hash directory lands in the project root. `run_simulation.py:275`
  does not pass `cache`, and the user-level `p21c.config['direc']` (here
  `~/21cmFAST-cache`, holding 2.2 GB of stale **v3** flat-named files) does
  **not** govern this path. `compute_initial_conditions`,
  `determine_halo_catalog`, and `perturb_halo_catalog` (§4.4) accept no cache
  argument at all. Recorded as requirement R3 with the one-line fix
  (`cache=p21c.OutputCache("<scratch>")`); §8's "written to the project root"
  stands, now with its mechanism.

- **Noted while auditing the parameter units:** `mean_galaxy_density` is
  declared `h³ Mpc⁻³` (`run_simulation.py:126`) but consumed as
  $P_{N,\mathrm{gal}} = 1/\bar n$ and reported in Mpc³ (§5.4). If the declared
  $h³$ is meant literally the shot noise is low by $h^{-3} = 3.3\times$.
  Recorded in §13.2 as unresolved; not yet tracked in `TODO.md`.

### Changed
- **`README.md`** — the `docs/HPC.md` row now points at §13 as the
  parameter-setup checklist.
- **`docs/HPC.md`** — header block gained a "setting up a run on a new
  machine?" pointer to §13.

<!-- ─── UV coefficient audit + HDF5 attr patch, 2026-08-12 ──────────────── -->

### Verified (no change required)
- **Mass-to-light coefficient in the galaxy-bias path audited — already
  correct.** `run_simulation.py:575` `uv_luminosity()` carries no hardcoded
  coefficient; it delegates to `src.conversions.sfr_to_Luv`, which returns
  `sfr / kappa_uv` with $\kappa_\mathrm{UV} = 1.15\times10^{-28}$, i.e.
  multiplication by **8.6957 × 10²⁷** = $1/1.15\times10^{-28}$ exactly
  (Madau & Dickinson 2014). The erroneous literal `1.15e28` appears nowhere in
  the codebase — that discrepancy was already removed by the 2026-08-04
  κ_UV/AB-zero-point unification (see the entry further down this file).
- **Scope correction for the record:** `uv_luminosity()` belongs to the
  *analytic HMF cross-check* branch, which sets `galaxy_bias_hmf_analytic`
  (5.39, **not adopted**). The "Estimating galaxy bias from the perturbed halo
  catalogue" step — the one that sets the adopted `b_g` — is a separate branch
  calling `select_euclid_halos` → `Luv_to_Muv(sfr_to_Luv(...))`
  (`src/analysis.py:865`). Both route through the same shared helper. Had the
  coefficient been wrong it would have moved only the cross-check: measured
  counterfactually, `1.15e28` gives $b_g^\mathrm{analytic}$ 5.3875 → 5.1962
  (−3.6 %, not the naive 1.32×, because the coefficient shifts the selection's
  mass threshold where $b_h(M)$ varies slowly).
- **Adopted galaxy bias re-measured from the stored 114,291,212-halo
  catalogue** with current code, avoiding a multi-hour re-simulation:
  **$b_g = 4.744191$** (49,315 halos selected, range 2.826 – 9.142,
  SFR window 0.79561 – 31.674 M☉ yr⁻¹), $f = 0.997672$,
  **$\beta = f/b_g = 0.210293$**. Unchanged from the recorded 4.744 / 0.2103 —
  no downstream numbers need updating.

### Changed
- **`outputs/lightcone_data.h5` — two root attributes patched in place**
  (data untouched, file still 2.76 GB): `photoz_uncertainty` 0.059 → **0.45**
  and `wedge_buffer` 0.02 → **0.0677**, so the analysis stage picks up both
  corrections without re-simulating. A provenance string is stored in the new
  `attrs_patched` attribute. `galaxy_bias` (33.3889) and `beta_rsd` (0.029880)
  were **deliberately not patched** — the stored `galaxy_overdensity` field
  carries that Kaiser boost, so the attributes still correctly describe the
  field on disk; patching them would make the file internally inconsistent.
  Full correctness still requires `bash submit_job.sh --sim force`.

  **Post-patch pipeline output:** $\sigma_r$ = **157.478 Mpc**, modes outside
  wedge **97/400 = 24.25 %**, total SNR **1.06 × 10⁻¹¹¹ σ** (numerically zero).
  Patching bumps the file mtime, so the spectra cache is recomputed (~0.6 s);
  the spectra are unchanged (large-scale $P_\times$ mean still −5.644 × 10³).
- **`docs/HPC.md`** — §11.7 rewritten for the partial-patch state, with the
  post-patch numbers and the unchanged $b_g$/$\beta$ measurement.

<!-- ─── Photo-z uncertainty convention, 2026-08-12 ──────────────────────── -->

### Fixed
- **Photometric redshift uncertainty: $\sigma_z\ 0.059 \to 0.45$
  (science-affecting).** `src.analysis.radial_smearing_length` computes
  $\sigma_r = c\,\sigma_z / H(z)$, which requires an **absolute** $\sigma_z$.
  The configured 0.059 was the *fractional* quantity $\sigma_z/(1+z)$ that
  surveys actually quote, used as if it were absolute — understating the
  radial smearing by a factor $(1+z) = 8$ at $z_\mathrm{obs} = 7$. The adopted
  0.45 corresponds to $\sigma_z/(1+z) = 0.056$, consistent with the Euclid
  photometric requirement $\sigma_z/(1+z) < 0.05$ (which gives 0.40 at $z=7$).

  **Downstream at $z_\mathrm{obs} = 7$:**

  | Quantity | 0.059 (old) | 0.45 (current) |
  |---|---|---|
  | $\sigma_r$ | 20.65 Mpc | **157.48 Mpc** |
  | $W$ at the smallest bin $k_\parallel = 0.0176\ \mathrm{Mpc}^{-1}$ | 0.936 | **0.021** |
  | $W$ at the first mode the wedge admits (0.1118) | 0.070 | $5\times10^{-68}$ |
  | $k_\parallel$ where $W = 0.5$ | 0.0574 Mpc⁻¹ | **0.0075 Mpc⁻¹** |
  | Total SNR (cached spectra, current buffer) | 0.0048 σ | **0.0000 σ** |

  The half-power scale now sits *below* the smallest $k_\parallel$ the box can
  sample, and would need $L_\mathrm{LOS} > 840$ Mpc to reach. Combined with the
  wedge this leaves no usable modes at all: the wedge admits only
  $k_\parallel > 0.1118\ \mathrm{Mpc}^{-1}$, where the kernel is
  $\sim 10^{-68}$. This is physical, not a bug — a photometric survey with
  $\sigma_z = 0.45$ at $z = 7$ retains essentially no LOS information. Any
  forecast must use spectroscopic redshifts, a box large enough to sample
  $k_\parallel \to 0$, or the angular (2D) cross-correlation.

  **Call sites updated:** `run_simulation.py` (config, also written to the
  HDF5 root attrs), `run_pipeline.py` (`data.get` fallback), `src/figures.py`
  (figure-label fallback), and the docstring of
  `src.analysis.radial_smearing_length`, which now states the absolute-vs-
  fractional convention explicitly — the ambiguity that caused this.
  **Deliberately not updated:** `tests/conftest.py` keeps 0.059 (its 64 Mpc
  synthetic box with $L_\mathrm{los} = 100$ Mpc would give a kernel that
  underflows to zero in every bin, making the damping and SNR tests vacuous);
  a comment there records why. `21cmfast_HERAxEuclid_lightcone.ipynb` still
  hardcodes 0.059 in its own config cell and is **not** covered by this fix.

  **Not yet in effect:** `run_pipeline.py` reads `photoz_uncertainty` from the
  HDF5 root attrs, so the stored `outputs/lightcone_data.h5` (written
  2026-06-15, $\sigma_z = 0.059$) keeps the old value until
  `bash submit_job.sh --sim force`. The same holds for `wedge_buffer` and
  `galaxy_bias`.
- **`README.md`, `docs/HPC.md`, `docs/project_update.md`** — $\sigma_z$ tables,
  the $\sigma_r$ evaluation, the photo-$z$ kernel table, and the wedge/photo-z
  conflict discussion updated. Test suite: 76 passed.

<!-- ─── HPC run specification, 2026-08-12 ───────────────────────────────── -->

### Added
- **`docs/HPC.md` — complete parameter-level specification of the HPC run.**
  A single reference covering: the conda environment and pinned
  `py21cmfast 4.1.1` stack; every `submit_job.sh` setting and every
  `run_pipeline.py` flag with its default; the full `run_simulation.py`
  configuration block; the `"simple"` template's `AstroParams`,
  `MatterOptions`, and `AstroOptions` as actually instantiated in the
  environment (`SOURCE_MODEL = CHMF-SAMPLER`, `USE_TS_FLUCT = False`,
  `SAMPLER_MIN_MASS = 1e8 M☉`, `N_THREADS = 1`); the derived geometry and mass
  resolution; every analysis formula with its **evaluated value at
  $z_\mathrm{obs} = 7$** ($t_\mathrm{sf} = 570.3$ Myr, $\sigma_r = 20.647$ Mpc,
  horizon slope 3.1509, FoV slope 0.37936, $T_\mathrm{sys} = 328.6$ K,
  $P_{N,21} = 3.7488$ mK² Mpc³, $P_{N,\mathrm{gal}} = 333.33$ Mpc³,
  Euclid window ↔ SFR 0.7956–31.674 M☉ yr⁻¹); the disk footprint; and a
  one-page reference table of every number.

  **New quantitative findings recorded there:**
  - **`L_los` is recorded as 200.0 Mpc while `lc_dist_Mpc` spans 3.4999 Mpc**
    (slice spacing 0.035385 Mpc) — a factor **56.5** disagreement. The
    attribute comes from `lightcone.lightcone_dimensions[2]`, which is
    $N_z \times$ the *transverse* cell size, and it propagates into
    $\Delta k_\parallel = 2\pi/200$, the wedge mask, the photo-$z$ kernel
    argument, and the Kaiser $\mu$ grid. Previously described in the changelog
    as "200 Mpc after rounding to cell boundaries"; the slice distances show it
    is not a rounding effect. Not currently tracked in `TODO.md`.
  - **Wedge-buffer change quantified end to end.** Recomputed from the cached
    spectra: $0.02 \to 0.0677\ \mathrm{Mpc}^{-1}$ moves the outside-wedge mode
    count from 105/400 (26.2 %) to **97/400 (24.2 %)** and the total SNR from
    0.0629 σ to **0.0048 σ** — a **13×** drop from only 8 lost bins, because
    those bins are the low-$k_\parallel$ ones the photo-$z$ kernel had not yet
    damped.
  - **Wedge vs photo-$z$ conflict made explicit.** At the smallest
    $k_\perp = 0.0140\ \mathrm{Mpc}^{-1}$ the wedge admits only
    $k_\parallel > 0.1118\ \mathrm{Mpc}^{-1}$, where $W = 0.070$ — every
    surviving mode is damped to ≤ 7 % of its amplitude.
  - **21cmFAST cache footprint: ~56 GB** in the gitignored
    `d1f8b93ecb5e05f9040e32ca2a1534a2/` directory at the project root
    (920 MB `InitialConditions.h5` + **~3.6 GB per node redshift**: 18 GB for
    the 5 smoke-test nodes, 36 GB for the 10 nodes of the briefly-set
    production range). This dominates the 2.76 GB of `outputs/` and is the
    quantity to size a scratch quota against.
  - Minor: the wedge buffer is converted at $h = 0.6766$ while the run's own
    `HUBBLE_CONSTANT = 67.36` implies $h = 0.6736$ (0.4 % difference).
  - The stored `lightcone_data.h5` carries 25 root attributes — 9 fewer than
    the 34 the current script writes.

### Changed
- **`README.md`** — `docs/HPC.md` added to the documentation table.
- **`PIPELINE.md`** — closing pointer now directs to `docs/HPC.md` for the
  parameter-level specification.

<!-- ─── Foreground wedge buffer, 2026-08-07 ─────────────────────────────── -->

### Changed
- **Foreground wedge buffer: $0.02 \to 0.0677\ \mathrm{Mpc}^{-1}$
  (science-affecting).** The old value was an unsourced placeholder, carried
  unchanged since the first commit (`dc68fdf`) with only the comment "safety
  margin beyond the horizon line". It is now set to the literature standard,
  $0.1\ h\ \mathrm{Mpc}^{-1}$, converted at $h = 0.6766$ (Planck 2018).

  **Motivation.** $0.1\ h\ \mathrm{Mpc}^{-1}$ is the additive buffer of the
  "moderate" foreground model of **Pober et al. (2014)**
  ([arXiv:1310.7031](https://arxiv.org/abs/1310.7031)), and is the default
  `horizon_buffer` in [21cmSense](https://github.com/rasg-affiliates/21cmSense).
  It traces to **Parsons et al. (2012a)**
  ([arXiv:1204.4749](https://arxiv.org/abs/1204.4749)), who showed that primary
  beam chromaticity combined with the tapering function applied in delay-space
  power spectrum estimation leaks foreground power $\sim 0.15\ h\
  \mathrm{Mpc}^{-1}$ beyond the horizon line. Note that **La Plante et al.
  (2023)**, the source of our wedge slope (their Eq. 10), applies *no* buffer
  and treats the bare horizon as the maximal-contamination case; the buffer
  here is the more conservative choice.

  **Scale check at $z_\mathrm{obs} = 7$.** Converting to delay via
  $k_\parallel = 2\pi\tau f_{21} H(z) / [c(1+z)^2]$, the old buffer was only a
  $\sim 50$ ns margin, against the $\gtrsim 300$ ns at which HERA sees
  chromatic calibration wings; the new value is $\sim 170$ ns. The old buffer
  was also comparable to the lightcone's fundamental mode
  ($k_\parallel^\mathrm{min} = 0.018\ \mathrm{Mpc}^{-1}$ for
  $L_\mathrm{LOS} = 351$ Mpc at $z = 6.5$–$7.5$), i.e. it excised barely one
  bin and left the outside-wedge mode fraction and total SNR optimistic.

  **Call sites updated:** `run_simulation.py` (config, also written to the HDF5
  root attrs), `21cmfast_HERAxEuclid_lightcone.ipynb` (config cell 4),
  `src/analysis.py` (`foreground_wedge_mask` default), `run_pipeline.py`
  (`data.get` fallback), and `tests/conftest.py` (synthetic fixture attr).
  Each now carries a one-line citation comment. The wedge mask continues to
  use the *horizon* slope rather than the FoV slope, which is the configuration
  Pober et al.'s buffer was calibrated against.

  **Downstream:** the fraction of modes outside the wedge and the total
  cross-correlation SNR both decrease; any previously quoted values (e.g. the
  26.2 % / 0.1 σ figures in `docs/project_update.md`) predate this change and
  need regenerating.
- **`README.md`** — parameter table entry for the wedge buffer now states the
  value in both $\mathrm{Mpc}^{-1}$ and $h\ \mathrm{Mpc}^{-1}$ with its source;
  Pober et al. (2014) and Parsons et al. (2012a) added to the reference list.

<!-- ─── Mass resolution reporting, 2026-08-05 ───────────────────────────── -->

### Added
- **`src/conversions.py` — mass-resolution helpers:** `mean_matter_density()`
  returns the comoving mean matter density $\bar\rho_m = \Omega_m
  \rho_{\mathrm{crit},0}$, and `cell_mass()` returns the mean matter mass
  enclosed by one cubic comoving cell, $M_\mathrm{cell} = \bar\rho_m
  L_\mathrm{cell}^3$ — the grid mass resolution.
- **`notebooks/plot_fields.ipynb` — "Mass resolution" block** in the parameter
  summary cell (cell 5, directly after "Grid & box geometry"). Prints
  $\bar\rho_m$, `DIM`, the high-res cell size, the mass resolution of both
  grids, and `SAMPLER_MIN_MASS`. For the production run this is
  $1.18\times10^{10}\ M_\odot$ per `DIM` cell (0.667 Mpc),
  $3.17\times10^{11}\ M_\odot$ per `HII_DIM` cell (2.00 Mpc), and a halo
  sampler floor of $1\times10^{8}\ M_\odot$. The block falls back to
  `DIM = 3 × HII_DIM` and the halo-catalogue minimum for HDF5 files written
  before these attributes existed, so it works with the stored output.
- **`run_simulation.py` — mass resolution in the startup summary** and four new
  HDF5 attributes: `DIM`, `hires_cell_size`, `M_cell_hires`, `M_cell_lores`,
  and `sampler_min_mass` (read from
  `inputs.simulation_options.SAMPLER_MIN_MASS`; `NaN` on the synthetic
  fallback path).
- **`tests/test_conversions.py`** — 7 tests for the new helpers: $\bar\rho_m$
  against `astropy.cosmology.Planck18`, the $H_0^2$ and $L^3$ scalings, the
  production-grid values, the 27× ratio between the two grids, total-box mass
  conservation, and array input. Suite is now 76 tests, all passing.

<!-- ─── Galaxy bias + production redshift range, 2026-08-04 ─────────────── -->

### Fixed
- **`run_simulation.py` — Sheth-Tormen ν convention (science-affecting):** the
  local helper `sheth_tormen_bias_from_nu(nu)` computed `a * nu**2`, but
  `hmf`'s `MassFunction.nu` is *already* the squared peak height (δ_c/σ)².
  Squaring it again inflated ν from the range 2.37–51.7 to 5.6–2670 and the
  bias with it, producing the `galaxy_bias = 33.39` recorded in the stored
  HDF5. Reproducing the original code path returns 33.39 exactly, and the
  corrected path returns 4.23 — confirming this, **not** the Euclid bright
  limit, was the dominant cause of the anomaly that
  `docs/project_update.md` had attributed to the magnitude cut and the SFR
  timescale. This is the same convention error already fixed in
  `notebooks/analysis.ipynb` ("Fix 1") but never back-ported. The local helper
  is deleted; `src.conversions.sheth_tormen_bias` is used instead.
- **`run_simulation.py` — SFR timescale:** `sfr_model` divided the stellar mass
  by a hardcoded 100 Myr, inconsistent with 21cmFAST's internal
  `t_STAR × t_H(z)` = 570.3 Myr at z = 7. The 5.70× SFR overestimate made every
  galaxy 1.89 mag too bright at the selection step. Now uses the new
  `src.analysis.star_formation_timescale`. On the analytic path this moves
  b_g from 4.23 to 5.39.
- **`run_simulation.py` — duplicated calibrations:** the script carried its own
  `Muv_to_Luv` with a 51.63 AB zero point (vs 51.60 in `src/conversions.py`)
  and an inverse-κ_UV factor of `1.15e28` (vs the correct 8.696e27 = 1/1.15e-28,
  a 32 % discrepancy). Both are now imported from `src/conversions.py`, which
  is the single source of truth.

### Changed
- **Galaxy bias is now measured from the halo catalogue, not the mean scaling
  relation.** `run_simulation.py` calls `src.analysis.select_euclid_halos` +
  `effective_galaxy_bias` — the same estimator the analysis stage uses — so
  Part 1 and Part 3 can no longer disagree. **Adopted b_g = 4.744** (range
  2.83–9.14 over 49,315 selected halos). The scatter-free analytic HMF
  integral is retained as a printed cross-check (5.39, ~14 % high because it
  misses the low-mass halos that 21cmFAST's log-normal scatter pushes into the
  magnitude window). Consequence: **β_rsd = f/b_g goes from 0.0299 to 0.2103**,
  a ~7× stronger Kaiser boost.
- **Redshift range: widened, then reverted — deliberately left at Δz = 0.01.**
  `z_min`/`z_max` were briefly set to the production range 6.5 / 7.5
  (L_LOS = 350.8 Mpc, N_z = 175, 10 node redshifts, verified numerically
  against Planck18) and then returned to 6.995 / 7.005. Reason: the
  power-spectrum estimator in `src/analysis.py` is inherited from the coeval
  notebook and assumes statistical homogeneity along the LOS. That assumption
  holds for the quasi-coeval Δz = 0.01 slab, so configuration and formalism
  currently match; at Δz = 1.0 it fails in four measured ways (see the
  `TODO.md` entry below). Widening Δz is now **gated on `TODO.md` §P0**, and
  `run_simulation.py` carries an explicit "do not widen without P0.1/P0.2"
  comment at the config block.
- **`outputs/lightcone_data.h5` is stale — because of the bias fix alone.**
  The corrected `galaxy_bias` (33.39 → 4.744) changes β_rsd (0.030 → 0.210)
  and therefore the `galaxy_overdensity` field, which carries the Kaiser
  boost. The other three fields are unaffected. `--sim auto` will not
  regenerate while the file exists; use `bash submit_job.sh --sim force`
  — cheap, since the redshift range is unchanged. Warnings added to
  `README.md`, `PIPELINE.md`, and `docs/project_update.md`.
- **`src/figures.py` — González+10 relation re-enabled** in
  `plot_stellar_mass_muv`, alongside Song+16. It had been commented out in
  `notebooks/plot_fields.ipynb` and was carried over disabled. Its ~0.2 dex
  higher normalisation follows from the constant-star-formation-history
  assumption.
- **`docs/project_update.md` rewritten** (2026-06-15 → 2026-08-04): corrects
  the b_g diagnosis, records both estimators, adds the analysis-stage results
  the pipeline now produces, marks the superseded run's numbers as such, and
  adds §12 costing the deferred 1 Gpc box.

### Added
- **`TODO.md`** — priority-ordered outstanding work. Records, with measured
  numbers, that widening to Δz = 1.0 would **invalidate three assumptions**
  baked into the power-spectrum estimator inherited from the coeval notebook,
  none of which has been addressed — which is why the range was reverted and
  the widening (P0.5) is gated behind them: (P0.1) `lc_redshifts` is uniform in redshift,
  so the LOS comoving cell varies by **20.4 %** across the box (2.214 → 1.838
  Mpc) while the FFT assumes a single 2.005 Mpc spacing — fixable with
  `RectilinearLightconer.between_redshifts`, which spaces slices uniformly in
  comoving distance; (P0.2) `compute_all_power_spectra` subtracts a single
  global mean, leaving the ⟨T_b⟩(z) evolution as a spurious LOS ramp that
  aliases into low-k_∥; (P0.3/P0.4) one FFT over Δz = 1 returns a
  redshift-averaged spectrum with no well-defined effective redshift, and it
  spans **22.28 MHz** against the **8 MHz** the noise model assumes — a 2.8×
  mismatch, resolvable by computing the spectrum in 8 MHz sub-bands. The old
  Δz = 0.01 slab had a 0.19 % cell spread and was quasi-coeval, which is why
  none of this surfaced before.
- **`src.analysis.star_formation_timescale(z, t_star)`** — 21cmFAST's
  `t_sf = t_STAR × t_H(z)`, matching astropy's Planck18 Hubble time to 0.08 %,
  plus **`stellar_mass_to_sfr`** and the `T_STAR_DEFAULT = 0.5` constant.
- **Five new HDF5 root attributes:** `galaxy_bias_method`,
  `galaxy_bias_hmf_analytic`, `t_STAR`, `sfr_timescale_yr` (all additive;
  `dataio.get()` tolerates their absence in older files).
- **Five tests** (69 total, all passing): the timescale against astropy, its
  t_STAR/redshift scaling, a 570-Myr regression guard against the 100 Myr
  value returning, `stellar_mass_to_sfr` inversion, and a guard that
  `sheth_tormen_bias` is not double-squaring its argument.
- **`docs/project_update.md` §12 — 1 Gpc box costed and deferred.** Reaching
  the Davies et al. (2025) box at the same 2 Mpc cell means `HII_DIM` 128 →
  500, `DIM` 384 → 1500: **59.6× the volume**, ~6.8 × 10⁹ halos, a **~163 GB**
  halo catalogue (from 2.74 GB), and **13.5 GB** for a single high-resolution
  initial-conditions array. Storage, not code, is the constraint — the
  parameters are already configurable and `--max-halos` exists for catalogues
  of this size. A 512 Mpc intermediate (8× volume, ~22 GB catalogue) is
  suggested as a scaling test.

<!-- ─── Full-pipeline driver, 2026-08-04 ────────────────────────────────── -->

### Added
- **`run_pipeline.py` — end-to-end pipeline driver.** One command now runs the
  whole workflow: an optional 21cmFAST simulation, the complete analysis
  (either fresh or from stored results), all figures, and a JSON summary.
  Previously `submit_job.sh` ran only `run_simulation.py`, and everything
  downstream lived in notebook cells that had to be executed by hand.
  Each stage is independently controllable:
  - `--sim {auto,force,skip}` — `auto` (default) invokes `run_simulation.py`
    as a subprocess only when `outputs/lightcone_data.h5` is missing, so an
    expensive 21cmFAST run can never be triggered by accident.
  - `--analysis {auto,force,skip}` — `auto` recomputes the power spectra only
    when the cache is missing or older than the simulation file (mtime
    comparison against the `source_mtime` attribute).
  - `--plots` selects among the groups `fields`, `halos`, `scaling`, `power`,
    `snr`, `bias` (or `all`/`none`); plus `--format`, `--dpi`, `--data`,
    `--products`, `--figdir`, `--summary`, `--max-halos`, `--m-uv-bright`,
    `--sim-script`, `--quiet`.
- **`src/analysis.py` — the Part 3 science, extracted from
  `notebooks/analysis.ipynb`.** `compute_cylindrical_cross_power`,
  `compute_all_power_spectra`, `horizon_wedge_slope`, `fov_wedge_slope`,
  `foreground_wedge_mask`, `radial_smearing_length`, `photoz_damping_kernel`,
  `hera_thermal_noise_power`, `cross_power_snr`, `total_snr`,
  `euclid_sfr_window`, `select_euclid_halos`, `effective_galaxy_bias`, plus
  the `EuclidSelection`, `BiasEstimate`, and `SNRResult` containers. Imports
  neither matplotlib nor py21cmfast.
  - The notebook computed the horizon slope twice from two different-looking
    expressions (cells 8 and 18). They are algebraically identical — the
    λ·f₂₁ factors cancel to `D_c H / [c(1+z)]` — so this is now one function.
- **`src/figures.py` — all 10 figures from both notebooks**, as functions
  returning a `Figure`. Forces the `Agg` backend on import, so figure
  generation is safe on a headless compute node. Includes the shared helpers
  `apply_plot_style`, `save_figure`, `fill_nan_nearest`, `eor_colormap`, and
  a `_binned_median` used by both percentile-band plots.
- **`src/dataio.py` — HDF5 I/O and caching.** `load_simulation` returns a
  typed `SimulationData` with accessors for the scalar metadata, and supports
  `max_halos` (uniform strided subsampling; the resulting
  `halo_sampling_factor` rescales the UVLF normalisation so number densities
  stay correct) and `load_halos=False` / `load_fields=False` for partial
  loads. `save_power_spectra` / `load_power_spectra` / `products_are_stale`
  back the `outputs/analysis_products.h5` cache.
- **`tests/` — 64 tests, the project's first suite.** `conftest.py` writes a
  synthetic `lightcone_data.h5` with the production schema (16² × 12 cells,
  4 000 halos), so nothing in the suite needs 21cmFAST and the whole run
  takes ~20 s. Covers the analysis functions (including the analytic
  white-noise normalisation $P = \sigma^2 V_\mathrm{cell}$ and the
  cross-spectrum sign), the I/O and cache-staleness logic, every figure
  function, and the pipeline's stage control end to end (with a stub
  simulation script standing in for 21cmFAST).
- **`outputs/analysis_products.h5` and `outputs/pipeline_summary.json`** as
  new pipeline products. The summary records ⟨x_HI⟩, the large-scale
  cross-spectrum sign, σ_r, both wedge slopes, noise levels, the total SNR,
  the Euclid selection counts, and ⟨b_g⟩.

### Changed
- **`submit_job.sh` now launches `run_pipeline.py` rather than
  `run_simulation.py`**, and forwards all of its arguments verbatim
  (`bash submit_job.sh --sim force`). The timing, CPU-hour accounting, and
  `sendmail` notification are unchanged; the email now also lists the figures
  written and points at the summary JSON. `JOB_NAME`, `CONDA_ENV`, and
  `PYTHON_SCRIPT` are overridable from the environment, so
  `PYTHON_SCRIPT=run_simulation.py bash submit_job.sh` restores the old
  simulation-only behaviour.
- **`README.md` and `PIPELINE.md`** document the driver, the stage-control
  flags, the new outputs, and the expanded flowchart. The README's Testing
  section no longer says "no `tests/` directory exists yet".

### Notes
- Verified against the stored fiducial run (128² × 100 cells, 114 M halos,
  2.76 GB HDF5): 34.5 s end to end with the simulation skipped, reproducing
  the notebook results — ⟨x_HI⟩ = 0.176, anti-correlated large-scale
  cross-spectrum, 26.2 % of modes outside the wedge, total SNR = 0.1 σ,
  ⟨b_g⟩ = 4.74 from 49 315 Euclid-selected halos.
- The halo catalogue is loaded only when a requested figure or the bias stage
  needs it, so `--plots power snr` skips 2.7 GB of reads entirely.

<!-- ─── Interactive inline figures, 2026-08-04 ──────────────────────────── -->

### Changed
- **All four notebooks — interactive inline plotting and constrained layout:**
  `21cmfast_HERAxEuclid_lightcone.ipynb`, `21cm_galaxy_cross_uncertainty.ipynb`,
  `notebooks/plot_fields.ipynb` and `notebooks/analysis.ipynb` now begin their
  imports cell with `%matplotlib widget` and set
  `plt.rcParams['figure.constrained_layout.use'] = True` immediately after the
  existing `plt.rcParams.update({...})` block. Figures are therefore
  pan/zoomable inline, and layout is resolved by the constrained-layout engine
  at draw time rather than per-figure.
- **All 49 `plt.tight_layout()` calls removed** (21 lightcone, 14 uncertainty,
  11 plot_fields, 3 analysis). Constrained layout and `tight_layout()` are
  mutually exclusive — a `tight_layout()` call switches the figure back to the
  tight engine and warns — so the calls are now both redundant and harmful.
  Existing explicit `constrained_layout=True` arguments (e.g. the `GridSpec`
  figure in `21cm_galaxy_cross_uncertainty.ipynb`) were left in place; they
  agree with the new global default.

### Added
- **`ipympl` dependency** (0.10.0, pulling `ipywidgets` 8.1.8,
  `jupyterlab_widgets` 3.0.16, `widgetsnbextension` 4.0.15), which provides the
  `widget` backend. Added to `env.yml` (notebook-support block) and pinned in
  `requirements.txt`. README gains a "Figure display" subsection under Usage
  and an `ipympl` row in the Requirements table. Note that the widget backend
  also needs `ipywidgets >= 8` in the front-end environment serving JupyterLab,
  which may differ from the `21cmfast` kernel environment.

<!-- ─── Notebook consolidation, 2026-08-03 ──────────────────────────────── -->

### Fixed
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — SFR unit bug (science-affecting):**
  The notebook consumed `perturbed_halos.sfr` directly, but py21cmfast v4
  returns it in **M☉ s⁻¹**, not M☉ yr⁻¹. Every UV luminosity, absolute
  magnitude, and $M_\mathrm{UV}$-based selection derived from it was therefore
  off by a factor of $3.15576\times10^{7}$ (~7.5 dex). This is the same bug
  already diagnosed in `docs/Low_SFR_fix.md` and already fixed in
  `run_simulation.py` (`_SEC_PER_YR`, line 334) and `notebooks/plot_fields.ipynb`
  (via `src.conversions.sfr_to_Muv`) — it had simply never been back-ported to
  this notebook. Added the conversion at both points where the halo catalogue
  SFR is extracted (the halo-catalogue cell and the UV-selection cell).
  **All $M_\mathrm{UV}$-derived numbers and figures in this notebook change;
  it must be re-run.** The stored outputs are now stale.

### Changed
- **`21cmfast_HERAxEuclid.ipynb` → `_archive/21cmfast_HERAxEuclid.ipynb`:**
  Moved to `_archive/` per the project's no-deletion policy. Measured overlap:
  50% of its distinct code lines are duplicated in the lightcone notebook, and
  the remainder is largely re-implemented boilerplate (`hubble_parameter`,
  matplotlib rcParams, the cylindrical power-spectrum binner, cosmology
  constants) rather than distinct science. Retained as the only **coeval**
  reference implementation. README updated throughout: overview list,
  repository tree, section 2 heading and archive notice, section 3 cross-
  reference, the three affected Requirements rows, and the Usage block.

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — UV calibration now shared:** The
  cell-local `K_UV = 1.15e-28` and its inline commentary were replaced with
  `src.conversions.sfr_to_Luv`, imported once in the imports cell, so the
  notebook can no longer drift from `run_simulation.py` and the analysis
  notebooks. The Kennicutt (1998) alternative is retained as a comment
  documenting the `kappa_uv=1.4e-28` override.

### Removed
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — internal redundancies:**
  - Duplicate `get_21cmfast_array()` definition (defined identically in two
    cells; the second copy also had a truncated docstring). Now defined once.
  - `OMEGA_M_0 = 0.315` redefined in the galaxy-bias cell despite being set in
    the ★ CONFIGURATION cell, contradicting that cell's "edit only this cell"
    contract. Values were identical, so **no numerical change**.
  - 11 redundant re-imports of `numpy`, `matplotlib.pyplot`, and `py21cmfast`
    across six cells, all already imported in the imports cell.
  - A leftover scratch cell (`[x for x in dir(p21c) if "halo" in x.lower()]`)
    and a trailing empty cell, both after the Summary section.

  Notebook drops from 35 to 33 cells; all 12.51 MB of stored outputs, and an
  uncommitted user edit in the HMF cell, were preserved verbatim.

- **`src/__pycache__/`** — deleted. Compiled bytecode is a regenerable build
  artifact already covered by `.gitignore`; removed at the author's explicit
  request (the project's no-deletion policy is intended for source files).

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — inconsistent bright-end
  $M_\mathrm{UV}$ cut (science-affecting):** The galaxy-bias cell used
  `M_UV_bright = -22` while the UV-selection cell running after it rebound
  `M_UV_bright = -22.66`, so the Sheth-Tormen bias integral was evaluated over
  a different magnitude range than the galaxy sample it was applied to.
  **−22.66 is the correct value** (confirmed by the author); −22 was the stale
  one. Both cuts are now defined once in the ★ CONFIGURATION cell and the
  scattered redefinitions removed, so the bias integral and the selection can
  no longer diverge. **The effective galaxy bias $b_\mathrm{gal}$ changes and
  the notebook must be re-run.**

  The dead `M_UV_limit = -18` config entry — set in the configuration cell but
  only ever consumed by one `print` and one plot label, never by the selection
  logic — was folded into the new `M_UV_faint` (same value, now actually
  authoritative). Both consumers were repointed.

<!-- ─── Documentation audit, 2026-08-03 ─────────────────────────────────── -->

### Fixed
- **`README.md` — Figure group 6 literature relations mis-described:** The text
  claimed *both* the Song+16 and González+10 $M_\star$–$M_\mathrm{UV}$
  relations were "defined in the code but commented out and are not rendered".
  In `notebooks/plot_fields.ipynb` (cell `wbf9ns3xkgb`) `song2016_z7` **is**
  actively plotted with the label `Song+16 $z\sim7$`; only `gonzalez2010_z7`
  is commented out. The table now carries an explicit per-row **Status**
  column ("Plotted" / "Defined but commented out").

- **`README.md`, `PIPELINE.md` — `sbatch submit_job.sh` is not a valid
  invocation:** `submit_job.sh` contains **zero** `#SBATCH` directives (no
  `--partition`, `--time`, or `--account`), and its own usage header specifies
  `bash submit_job.sh`. Submitting it via `sbatch` would be rejected or
  silently assigned default resources. All four occurrences corrected to
  `bash submit_job.sh` (README workflow block, README Usage block, PIPELINE
  Mermaid node `A`, PIPELINE "Running it" block). A note in both documents now
  states that `#SBATCH` directives must be added before the script can be
  submitted as a true batch job. The README file table and the PIPELINE stage
  table no longer describe it as a "SLURM batch submission script".

- **`env.yml` — did not reproduce the environment:** The file listed only
  `python=3.11` and `ipykernel`, so `conda env create -f env.yml` produced an
  environment in which no notebook or script could run. Added the FFTW/GSL
  build dependencies (required *before* the pip step or the 21cmFAST C
  extension build fails with `cannot find -lfftw3f`), the scientific stack
  (`numpy`, `scipy`, `astropy`, `matplotlib`, `h5py`), `jupyter`, and a `pip:`
  section installing `21cmFAST==4.1.1` and `hmf>=3.5`. Version floors were
  chosen to be satisfied by the versions currently installed in the working
  `21cmfast` environment.

### Added
- **`README.md` — Installation section:** Previously absent entirely. Documents
  `conda env create -f env.yml`, the `py21cmfast.__version__` verification
  step, why FFTW/GSL must precede the pip install, `requirements.txt` for
  reproducing the exact pinned environment, and a pointer to
  `docs/INSTALL_21cmFASTv4.md` for the CSD3-specific quota, `CONDA_NO_PLUGINS`,
  and FFTW-linking problems that `conda env create` cannot resolve on its own.

- **`README.md` — Repository structure tree:** New top-level section giving an
  annotated file tree, marking `outputs/` and `resources/` as gitignored.

- **`README.md` — Documentation index:** New section with a table linking all
  seven companion documents. Six of them (`PIPELINE.md` beyond a single inline
  link, `CHANGELOG.md`, and four of the five `docs/*.md` files) were previously
  unreachable from the README — only `docs/Low_SFR_fix.md` was ever cited.

- **`README.md` — fiducial parameters for the HPC pipeline (§4):** Section 4
  described `run_simulation.py` only as "a refactored version of notebook 3",
  implying it inherits notebook 3's $z = 6.5$–$7.5$ range. It does not. Added a
  full parameter table plus a **thin-slab warning**: the committed config spans
  $\Delta z = 0.01$ ($z = 6.995$–$7.005$, $L_\mathrm{LOS} = 3.5$ Mpc at
  $z = 7$). The cell-size-matched slice count would be $N_z = 2$, but
  `minimum_los_slices = 100` (`run_simulation.py:170`) raises it to 100,
  yielding a 0.035 Mpc LOS cell against a 2 Mpc transverse cell — a ~57×
  line-of-sight oversampling. This is a smoke-test slab with negligible
  redshift evolution, not a science configuration; the production equivalent
  ($z = 6.5$–$7.5$) gives $L_\mathrm{LOS} = 350.8$ Mpc and $N_z = 175$. A
  cross-reference was added to `PIPELINE.md`.
  **The configuration itself was left unchanged — this is a documentation fix
  only, and the redshift range remains a science decision for the author.**

- **`README.md` — `src/FOV_to_cMpc.py` documentation:** The "Source Modules"
  section covered only `conversions.py`, leaving this module entirely
  undocumented. Added its CLI usage example, the full argument table
  (`--area-deg2`, `--z-min`, `--z-max`, `--n-z`), both importable functions
  (`survey_volume_from_area`, `transverse_comoving_size_from_area`), and a note
  that it overlaps with `conversions.volume_from_area` but additionally returns
  the intermediate solid angle and per-steradian volume.

- **`README.md` — Testing section:** Records the `tests/test_<module>.py`
  convention and the `conda run -n 21cmfast pytest tests/ -v` command, and
  states plainly that **no `tests/` directory exists yet**, so `src/` is
  currently untested. Names the round-trip identities and the
  `sheth_tormen_bias` squared-peak-height convention as first candidates.

- **`README.md` — Part 1 prerequisite note:** Clarifies that `outputs/` is
  gitignored, so a fresh clone has no `lightcone_data.h5` and Parts 2–3 fail at
  the loading cell until `run_simulation.py` has been run. Distinguishes
  "independent of the simulation run" (never imports 21cmFAST) from "ships with
  the repository" (it does not). `resources/` noted as local-only.

### Changed
- **`README.md` — notebook structure lists realigned to actual notebook
  headers:** Both lists were renumbered prose that did not match the notebooks.
  - Notebook 3: the brightness-temperature plot was listed as "7b" but is
    **5b** in the notebook, and every item from "Kaiser RSD" onward was offset
    by one (README 6 = notebook §4, README 7 = §5). Renumbered to match.
  - Notebook 1: the 10-item list did not correspond to the notebook's 11
    sections or their order. Rewritten against the real headers, including an
    explicit warning that the numbering is non-monotonic — §2c sits *after*
    §3a, and two distinct subsections are both labelled §3b.

<!-- ─── Earlier unreleased work ─────────────────────────────────────────── -->

### Added
- **`PIPELINE.md`** — new top-level document summarising the HPC pipeline
  (`submit_job.sh` → `run_simulation.py` → `outputs/lightcone_data.h5` →
  `notebooks/plot_fields.ipynb` / `notebooks/analysis.ipynb`), including a
  Mermaid flowchart, a stage table, and the HDF5 output schema.

### Fixed
- **`notebooks/plot_fields.ipynb` — Bouwens et al. (2021) journal citation:**
  All occurrences of "ApJ 908, 24" corrected to "AJ 162, 47". The paper is
  published in *The Astronomical Journal*, not *The Astrophysical Journal*.
  Affected cells: section-5 header markdown (`ev0xfh6438f`) and the
  `literature` list comment (`knrobm93rko`). Parameter values (φ\* = 0.19×10⁻³,
  M\* = −21.15, α = −2.06 from Table 5, z = 6.8) were already correct.

- **`notebooks/plot_fields.ipynb` — Finkelstein et al. (2015) Schechter
  parameters (`knrobm93rko`):** The `literature` list used values from Table 3
  (the galaxy catalogue), not Table 4 (the Schechter fits). Corrected to the
  Table 4, z = 7 values:
  - φ\*: `0.74e-3` → `1.57e-4` Mpc⁻³
  - M\*: `−20.81` → `−21.03`
  - α: `−1.87` → `−2.03`
  Comment updated from "Table 3" to "Table 4".

- **`notebooks/plot_fields.ipynb` — González citation corrected to 2010
  paper:** The function `gonzalez2011_z7` and all associated labels and
  comments cited "González et al. (2011, ApJ 736, 133)", which resolves to an
  unrelated Galactic Center paper (An et al. 2011). The intended reference is
  González et al. (2010, ApJ 713, 115), whose constant-SFH SED fitting at
  z ~ 7 motivates the ~0.2 dex higher normalisation relative to Song+16.
  Renamed to `gonzalez2010_z7`; labels updated from "González+11" to
  "González+10" throughout cells `1c28abd8`, `zcmcw9903x`, and `wbf9ns3xkgb`.

- **`README.md` — Bouwens et al. (2021) values and citation:** The UVLF table
  in Figure group 5 listed wrong journal ("ApJ 908, 24"), wrong φ\*
  (2.9×10⁻⁴), wrong M\* (−21.03), and wrong α (−2.03). Corrected to AJ 162,
  47 with Table 5 (z = 6.8) values: φ\* = 1.9×10⁻⁴, M\* = −21.15, α = −2.06.
  Corresponding References entry also corrected.

- **`README.md` — Finkelstein et al. (2015) values and table reference:**
  Figure group 5 UVLF table listed Table 3 with wrong Schechter parameters
  (φ\* = 7.4×10⁻⁴, M\* = −20.81, α = −1.87). Corrected to Table 4, z = 7:
  φ\* = 1.57×10⁻⁴, M\* = −21.03, α = −2.03.

- **`README.md` — González citation in Figure group 6 table:** "González et
  al. (2011, ApJ 736, 133)" replaced with "González et al. (2010, ApJ 713,
  115)" for the same reason as the notebook fix above.

### Changed
- **`README.md` — References section:** Added arXiv links to all entries that
  previously lacked them:

  | Paper | arXiv |
  |-------|-------|
  | Park et al. (2019) | [1809.08995](https://arxiv.org/abs/1809.08995) |
  | Bouwens et al. (2021) | [2102.07775](https://arxiv.org/abs/2102.07775) |
  | Finkelstein et al. (2015) | [1410.5439](https://arxiv.org/abs/1410.5439) |
  | Speagle et al. (2014) | [1405.2041](https://arxiv.org/abs/1405.2041) |
  | Schreiber et al. (2015) | [1409.5433](https://arxiv.org/abs/1409.5433) |
  | Song et al. (2016) | [1507.05636](https://arxiv.org/abs/1507.05636) |
  | Murray et al. (2013) | [1306.6721](https://arxiv.org/abs/1306.6721) |

  Bardeen et al. (1986), Kaiser (1987), and Oke & Gunn (1983) predate arXiv
  and have no preprint record. Two new entries added for both González papers:
  - González et al. (2010), ApJ 713, 115 — [arXiv:0909.3517](https://arxiv.org/abs/0909.3517)
  - González et al. (2011), ApJL 735, L34 — [arXiv:1008.3901](https://arxiv.org/abs/1008.3901)

### Added
- **`src/conversions.py` — UV luminosity–SFR conversions:**
  - `_KAPPA_UV_MADAU14 = 1.15e-28` — module-level constant for the Madau &
    Dickinson (2014) UV–SFR calibration factor [M☉ yr⁻¹ / (erg s⁻¹ Hz⁻¹)],
    Chabrier (2003) IMF, rest-frame ~1500 Å.
  - `Luv_to_sfr(Luv, kappa_uv=1.15e-28)` — UV luminosity [erg s⁻¹ Hz⁻¹] →
    SFR [M☉ yr⁻¹] via `SFR = κ_UV × L_UV`.
  - `sfr_to_Luv(sfr, kappa_uv=1.15e-28)` — inverse: SFR [M☉ yr⁻¹] → UV
    luminosity [erg s⁻¹ Hz⁻¹].
  - `sfr_to_Muv(sfr, kappa_uv=1.15e-28)` — convenience chain: SFR →
    `sfr_to_Luv` → `Luv_to_Muv`, giving the AB magnitude directly.

- **`src/conversions.py` — Sheth-Tormen halo bias:**
  - `sheth_tormen_bias(nu_sq, delta_c=1.686, a=0.707, p=0.3)` — Eulerian
    linear halo bias $b(\tilde\nu) = 1 + (a\tilde\nu-1)/\delta_c +
    2p/(\delta_c(1+(a\tilde\nu)^p))$ where $\tilde\nu = (\delta_c/\sigma)^2$
    is the squared peak height as returned by `hmf.MassFunction.nu`.
    Documented with an explicit warning that `hmf` ≥ 3.x stores `mf.nu` as
    the *squared* peak height (Sheth & Tormen 1999 convention), not
    $\delta_c/\sigma$.

- **`notebooks/analysis.ipynb` — halo catalogue loading (Section 0, HDF5
  cell):** Extended the HDF5 load block to read the halo catalogue stored by
  `run_simulation.py`: `sfr_cat`, `halo_masses`, `halo_coords`,
  `stellar_masses` (all under `halo_catalog/`). A count of halos with
  `SFR > 0` is printed on load. Arrays are empty when the simulation was run
  without 21cmFAST.

- **`notebooks/analysis.ipynb` — Section 4: Euclid luminosity and SFR cuts
  (three cells):**
  - **Cell 12** — Converts the HDF5 `M_UV_limit` attribute to a UV luminosity
    floor via `Muv_to_Luv` from `src/conversions.py`. Sets up `sys.path` so
    all subsequent cells can import from `src/`.
  - **Cell 13 (Version 1 — SFR bounds)** — Derives the SFR selection window
    `[SFR_min, SFR_max]` from the Euclid magnitude window
    `[M_UV_bright=-22, M_UV_limit]` using `Luv_to_sfr`. Applies the SFR
    window as a direct cut on the `sfr_cat` halo catalogue and prints the
    selected count and SFR range.
  - **Cell 14 (Version 2 — per-halo magnitude assignment)** — For each halo
    with `SFR > 0`, computes `L_UV = sfr_to_Luv(SFR)` then
    `M_UV = Luv_to_Muv(L_UV)`, and applies the Euclid magnitude window as an
    explicit $M_\mathrm{UV}$ cut. Prints selected counts and the full
    M_UV / L_UV / SFR ranges. Mathematically equivalent to Version 1;
    confirms internal consistency of the conversion chain.

- **`notebooks/analysis.ipynb` — Section 4.3: effective galaxy bias from the
  halo catalogue (cells 15–16):** Implements the `temp.py` logic with two
  physics bugs corrected (see **Fixed** below):
  - Sanitises `sfr_cat` (NaN/negative → 0), selects halos with `SFR > 0` and
    `halo_mass > 0`, converts SFR → M_UV via `sfr_to_Luv` + `Luv_to_Muv`,
    and applies the Euclid magnitude window.
  - Converts selected 21cmFAST halo masses from M☉ to M☉ h⁻¹ using
    `h = HUBBLE_CONSTANT / 100` before querying the HMF grid.
  - Builds a `MassFunction` grid spanning the selected mass range (+0.5 dex
    margin), computes `sheth_tormen_bias(mf.nu)`, and interpolates to each
    selected halo mass to obtain `selected_biases`.
  - Reports the number-weighted mean effective galaxy bias
    $\langle b_g \rangle$ and the bias range.
  - Produces a diagnostic plot: histogram of selected halo masses (left axis)
    with the Sheth-Tormen bias curve overlaid (right axis, red dashed), and a
    dotted line marking $\langle b_g \rangle$.
  - The entire block is guarded by `if len(sfr_cat) > 0` so the notebook
    runs end-to-end when 21cmFAST is unavailable.

### Fixed
- **`notebooks/analysis.ipynb` — Section 4.3: ν² double-squaring bug
  (from `temp.py`):** `hmf.MassFunction.nu` in `hmf` ≥ 3.x stores the
  *squared* peak height $\tilde\nu = (\delta_c/\sigma)^2$ (Sheth & Tormen
  1999 convention). `temp.py` treated `mf.nu` as $\delta_c/\sigma$ and
  squared it again in the bias formula, computing $a\tilde\nu^2$ instead of
  the correct $a\tilde\nu$. This caused a **4–17× overestimate** of the
  galaxy bias (e.g. $\langle b_g\rangle \approx 60$ instead of $\approx 5$
  for a typical Euclid sample at $z \sim 7$). Fixed by using
  `sheth_tormen_bias(mf.nu)` from `src/conversions.py`, which correctly
  treats its argument as $\tilde\nu$.

- **`notebooks/analysis.ipynb` — Section 4.3: halo mass unit mismatch
  (from `temp.py`):** `hmf.MassFunction.m` returns masses in M☉ h⁻¹, while
  21cmFAST `perturbed_halos.halo_masses` (and the `halo_masses` array loaded
  from the HDF5) are in M☉. `temp.py` passed M☉ masses directly to the
  log-spaced HMF grid defined in M☉ h⁻¹ units, producing a systematic
  $\log_{10}(h) \approx -0.17$ dex offset in the bias interpolation. Fixed
  by converting selected halo masses to M☉ h⁻¹ via
  `selected_mass * (HUBBLE_CONSTANT / 100)` before computing `log10_m_min`,
  `log10_m_max`, and the interpolation argument.

### Added
- **HPC lightcone pipeline** — the monolithic `21cmfast_HERAxEuclid_lightcone.ipynb`
  has been refactored into three self-contained parts for efficient cluster use:

  - **`run_simulation.py`** (Part 1 — batch script): runs the 21cmFASTv4 lightcone
    simulation, constructs the galaxy density field from `halo_sfr`, estimates the
    galaxy bias via HMF integration over the Euclid UV magnitude range, applies
    Kaiser redshift-space distortions in Fourier space, and saves all outputs to
    `outputs/lightcone_data.h5` with gzip compression. Uses `matplotlib.use("Agg")`
    for headless HPC execution. No logic changes from the notebook — same algorithms,
    same parameters, same comments.

  - **`notebooks/plot_fields.ipynb`** (Part 2 — visualisation notebook): loads
    `outputs/lightcone_data.h5` and reproduces all field plots from the original
    notebook: halo catalogue scatter plots and SFR distributions (Cells 10–11),
    the three-panel lightcone slice (Cell 18), and the wide-format EoR brightness
    temperature plot (Cell 20). Gracefully skips halo catalogue cells when no
    catalogue is available (synthetic fallback).

  - **`notebooks/analysis.ipynb`** (Part 3 — calculation notebook): loads
    `outputs/lightcone_data.h5` and performs all post-simulation calculations:
    `compute_cylindrical_cross_power` (Cell 22), foreground wedge geometry and
    power spectrum plots (Cells 23–24), photo-$z$ damping and wedge excision
    (Cell 26), per-mode SNR map and total detection significance (Cell 28), and
    the summary table (Cell 29).

  - **`submit_job.sh`**: SLURM batch submission script that activates the
    `21cmfast` conda environment and runs `run_simulation.py`. Configurable
    wall-time, memory, and partition via `#SBATCH` directives.

  All scalar parameters and metadata are stored as HDF5 attributes in
  `outputs/lightcone_data.h5` so that Parts 2 and 3 require no configuration
  beyond pointing at the output file. The kernelspec for both notebooks is
  set to `21cmfast`.



### Added
- **`src/conversions.py`** — new module of cosmological conversion utilities for
  high-redshift galaxy surveys. Functions:
  - `Muv_to_Luv(Muv)` — converts absolute UV AB magnitude to monochromatic UV
    luminosity in erg s⁻¹ Hz⁻¹ (Oke & Gunn 1983; Madau & Dickinson 2014).
  - `Luv_to_Muv(Luv)` — inverse conversion from luminosity to AB magnitude.
  - `survey_area_from_volume(volume_mpc3, z_min, z_max, cosmo=None)` — infers
    the sky area in deg² that corresponds to a given comoving survey volume
    (Mpc³) over a redshift interval, using Simpson integration of the
    differential comoving volume $\mathrm{d}V/\mathrm{d}z\,\mathrm{d}\Omega$
    (Hogg 1999; Astropy).
  - `area_deg2_to_steradians(area_deg2)` — unit conversion from deg² to sr.
  - `volume_from_area(area_deg2, z_min, z_max, cosmo=None, n_z=1000)` —
    computes the comoving survey volume (Mpc³) for a given sky area and
    redshift range; inverse of `survey_area_from_volume`.

  All functions accept scalar or array inputs and default to the Planck18
  cosmology from astropy; a custom `astropy.cosmology` object may be passed
  via the `cosmo` argument.

### Fixed
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — lightcone slice orientation:**
  Panels 2 and 3 of Section 5 previously had the LOS and transverse axes
  swapped. The `brightness_temp_field[:, mid_y, :]` and
  `neutral_fraction[:, mid_y, :]` arrays (shape `(HII_DIM, N_z)`) were
  incorrectly transposed with `.T`, placing the LOS on the y-axis and the
  transverse direction on the x-axis. Removed the transpose so that
  `imshow(origin="lower")` correctly maps rows→transverse (y) and cols→LOS (x),
  consistent with every published 21cm lightcone figure. The `extent` parameter
  was updated accordingly to `[lc_dist_Mpc[0], lc_dist_Mpc[-1], 0, BOX_LEN]`.

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — secondary redshift axis:**
  The redshift annotation used `twinx()`, adding a second y-axis on the right.
  After the orientation fix the LOS is on the x-axis, so the annotation must
  use `twiny()` with tick positions computed via
  `np.interp(z_ticks, lc_redshifts, lc_dist_Mpc)`.

### Added
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — brightness temperature evolution
  plot (Section 5b):** New wide-format (16×3.5") lightcone slice cell styled
  after the canonical Mesinger & Furlanetto (2007) figure. Features:
  - Custom EoR colourmap (`EoR21`): dark blue-black (ionised, δT_b ≈ 0) →
    blue → cyan → yellow → orange → near-white (neutral, high δT_b)
  - LOS (comoving distance) on the x-axis; transverse distance on the y-axis
  - Dual x-axes: comoving distance [Mpc] on the bottom, redshift z on the top
    (using `twiny()` with interpolated tick positions)
  - Title reports the observed frequency range
    ($f_\mathrm{obs} = F_{21}/(1+z)$) alongside z and ⟨x_HI⟩

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — configurable minimum LOS slices:**
  `minimum_los_slices` (default 100) replaces the previous hardcoded floor of
  10. For narrow redshift ranges where the natural slice count (L_LOS / cell_size)
  is small, this ensures the k_∥ grid is adequately sampled.
  Validated with z_min=6.5, z_max=6.505, 100 slices: notebook executes without
  errors, 21cmFAST returns a (128, 128, 100) box (L_LOS = 200 Mpc after rounding
  to cell boundaries), empty bins fall from 315/400 to 145/400, and the
  large-scale cross-spectrum sign is correctly negative (anti-correlated).

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — Section 3b synthetic galaxy field
  fallback:** When the lightcone `halo_sfr` field is all-zero (e.g. 21cmFAST
  not installed or CHMF-SAMPLER placed no halos), Section 3b automatically
  generates a synthetic galaxy overdensity by Poisson-sampling the matter
  density field. This allows the power-spectrum and SNR cells to run
  end-to-end without a 21cmFAST installation.

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — galaxy bias estimation (Cell 13):**
  New cell estimates the luminosity-weighted linear bias by integrating a
  Schechter UV luminosity function over the Euclid magnitude range
  $M_\mathrm{UV} \in [-24, -18]$. The result informs the default
  `galaxy_bias` value used in the Kaiser RSD correction.

### Added
- **`21cmfast_HERAxEuclid_lightcone.ipynb`** — lightcone counterpart to the
  coeval simulation notebook. Uses `RectilinearLightconer` + `run_lightcone`
  (21cmFASTv4) to produce a self-consistent lightcone over a configurable
  redshift range ($z_\mathrm{min}$–$z_\mathrm{max}$, default 6.5–7.5).
  Key differences from the coeval version:
  - Non-cubic $(N_\perp \times N_\perp \times N_z)$ box with separate transverse
    and LOS cell sizes handled throughout
  - Fields accessed via `lightcone.lightcones['field_name']`; neutral fraction
    stored as `'neutral_fraction'` (not `'xH_box'` as in coeval)
  - Visualisation includes LOS slice panels showing redshift evolution of
    $\delta T_b$ and $x_\mathrm{HI}$ alongside a transverse (x–y) slice
  - Updated `compute_cylindrical_cross_power` accepts separate `box_len_perp`
    and `box_len_los` arguments for non-cubic boxes
  - All user-adjustable parameters consolidated into a single clearly marked
    configuration cell

### Fixed
- `21cmfast_HERAxEuclid.ipynb`
  - Updated foreground and horizon wedge prescription
  - Fixed Nan appearance due to log binning in 2d power spectra plots

### Added
- `21cmfast_HERAxEuclid.ipynb` — new notebook demonstrating an end-to-end
  HERA × Euclid 21 cm–galaxy cross-correlation workflow using 21cmFASTv4
  (Davies et al. 2025, arXiv:2504.17254). Covers simulation, galaxy field
  construction, 2D cylindrical power spectra, photo-z damping, foreground
  wedge excision, and per-mode SNR estimation.

### Fixed
- **`21cmfast_HERAxEuclid.ipynb` — Cell 4 (galaxy field construction):**
  The cell previously raised an unconditional `RuntimeError` when 21cmFAST
  was installed, blocking all subsequent cells. Fixed by implementing the
  galaxy density field using `coeval.halobox.get('halo_sfr')`, the per-cell
  SFR density field provided by the 21cmFASTv4 `HaloBox` API
  (`py21cmfast >= 4.1.1`). The SFR density field is converted to an
  overdensity $\delta_{\rm gal} = \rm SFR / \langle SFR \rangle - 1$, which
  correctly traces the galaxy distribution and produces the expected
  large-scale anti-correlation with the 21 cm brightness temperature field
  (cross-spectrum sign: **NEGATIVE** on large scales ✓).

### Changed
- **`21cmfast_HERAxEuclid.ipynb` — Section 3 markdown:** Updated the list of
  coeval output fields to reflect the actual 21cmFASTv4 API: `halo_field`
  (incorrect) → `halobox` (the `HaloBox` object containing per-cell gridded
  quantities including `halo_sfr`, `n_ion`, etc.).
- **`21cmfast_HERAxEuclid.ipynb` — Section 4 markdown:** Clarified that the
  `HaloBox` API exposes cell-averaged quantities rather than individual halo
  catalogues, so a strict per-halo $M_{\rm UV}$ cut requires a lightcone
  post-processing step. Added a note explaining why the SFR-density proxy is
  a valid and physically motivated tracer of the Euclid-observable galaxy
  population.

---

## Notes on 21cmFASTv4 `HaloBox` API

In 21cmFAST v4.1+, `coeval.halobox` is a `HaloBox` object whose arrays are
accessed via `.get('<field_name>')`. Available fields include:

| Field | Description |
|-------|-------------|
| `halo_sfr` | Total SFR per cell, summed over all halos [internal units] |
| `n_ion`    | Number of ionizing photons per cell |

Individual halo positions and UV magnitudes are not exposed by this API.
For per-halo catalogues, use the raw halo output from a lightcone run.
