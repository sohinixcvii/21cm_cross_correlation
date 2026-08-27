# The Planned Simulation — Parameters, Compute and Storage

**What this document is.** The specification of the run this repository is
*about to make*: every parameter it will use, every number those parameters
imply, and what it will cost in cores, memory, scratch and wall time. It is
written to be handed to whoever grants the allocation, and to be checked
against before `sbatch` is typed.

**Verified:** 2026-08-24, against `HEAD` on branch `test/methods`
(`350a3fd`). Test suite: **245 passed**
(`conda run -n 21cmfast pytest tests/ -q`).

**Companions.** [`HPC.md`](HPC.md) is the parameter-level ground truth for the
run *as the code stands today*, including every known inconsistency;
[`uncertainty_budget.md`](uncertainty_budget.md) is the photo-*z* / wedge /
noise / SNR chain; [`../TODO.md`](../TODO.md) §P0 is the estimator work that
gates the redshift range below. This document is the forward-looking one: the
cost of the next run, not the audit of the last.

**Provenance of every number here** is marked:
**[M]** measured on a real run, **[E]** extrapolated from a measured run,
**[C]** computed from the configuration for this document, **[T]** a target or
recommendation.

---

## 0. TL;DR

One 21cmFAST v4.1.1 lightcone in a **486.33 Mpc box on a 256³ grid**
(1.90 Mpc cells), sized from the **Euclid Deep Field Fornax footprint,
10 deg² at z = 7**, spanning **z = 6.55 → 7.45** (±1 σ_z of photometric
depth, 315.6 Mpc, 166 LOS slices, 9 node redshifts), drawing a
**~9.4 × 10⁸ halo catalogue** down to a 10⁸ M☉ sampler floor, then the full
analysis chain: galaxy overdensity → Euclid selection → bias → Kaiser RSD →
cylindrical power spectra → photo-*z* damping → wedge → HERA noise → SNR →
18 figures.

Costed against the machine it will actually run on — 2 × AMD EPYC 9374F,
~1.5 TiB, local NVMe, **no scheduler** (§4.1).

| | Requirement | Against what the machine has |
|---|---|---|
| **Threads** | **32**, one socket, set explicitly **[T]** | 64 physical cores / 128 logical |
| **Memory** | ~56 GB peak resident **[E]** | ~1.5 TiB — not a constraint |
| **Disk** | ~265 GB per run **[E]** on `/nvme1` or `/nvme4` **[T]** | 851 / 836 GB free; **not** the 144 GB NFS home |
| **Wall time** | **15–40 min**, budget 1 h **[E]** | no queue, nothing billed |
| **Serial-equivalent work** | 2–4 core-hours **[E]** | — |

**Two things block submission today**, both in §7: the halo catalogue
overflows a signed 32-bit index at this box size (1.31 × `INT_MAX`), and the
power-spectrum estimator is not yet valid over a Δz = 0.9 lightcone
(`TODO.md` §P0.1–P0.4). §7 gives the numbers for each fix.

---

## 1. What is planned versus what is committed

The transverse box has already been switched to the survey footprint
(commit `81e08ef`, 2026-08-19). The redshift range has **not**: the committed
`run_simulation.py` still runs the Δz = 0.01 smoke-test slab, deliberately,
because widening it before `TODO.md` §P0 lands produces spectra that cannot be
trusted (§7.2).

| | Committed today | Planned run |
|---|---|---|
| `BOX_LEN` / `HII_DIM` / `DIM` | 486.33 Mpc / 256 / 768 | *unchanged* |
| `z_min` → `z_max` | 6.995 → 7.005 (slab) | **6.55 → 7.45** |
| `L_LOS` | 3.50 Mpc **[C]** | **315.60 Mpc** **[C]** |
| `N_z` | 100 (floor binds) | **166** (natural) |
| Node redshifts | 5 | **9** |
| Frequency span | ~0.9 MHz | **20.04 MHz** **[C]** |
| Cache | ~140 GB **[E]** | ~245 GB **[E]** |
| Where set | `run_simulation.py:173–174` | same two lines |

`SURVEY_Z_MIN` / `SURVEY_Z_MAX` (6.55 / 7.45) are already computed from the
footprint block and written to the HDF5; the slab simply overrides them.
Flipping to the planned run is a two-line edit plus `--sim force`.

---

## 2. Configuration — every parameter

Everything below the survey block is derived, not chosen. Line numbers are
`run_simulation.py` unless stated.

### 2.1 Survey footprint → box geometry

| Parameter | Line | Value | Note |
|---|---|---|---|
| `SURVEY_AREA_DEG2` | 109 | **10.0 deg²** | Euclid Deep Field Fornax, RA 03:31:43.6, Dec −28:05:18.6 |
| `SURVEY_Z_CENTRAL` | 110 | **7.0** | |
| `photoz_uncertainty` σ_z | 120 | **0.45** | **Absolute**, not σ_z/(1+z). Equals σ_z/(1+z) = 0.056 at z = 7, against Euclid's < 0.05 requirement |
| `PHOTOZ_N_SIGMA` | 127 | **1** | Explicit choice: the box spans ±1 σ_z. `= 2` gives Δz = 1.80, z = 6.10–7.90, L∥ = 634.9 Mpc **[C]** |
| `SURVEY_DELTA_Z` | 128 | **0.512** | `2 × PHOTOZ_N_SIGMA × σ_z` |
| `target_cell_size_mpc` | 141 | **2.0 Mpc** | Preserves the resolution of the old 256 Mpc / 128³ grid |

`survey_area_to_box_size()` (`src/conversions.py`) turns those into the grid:

| Step | Formula | Value **[C]** |
|---|---|---|
| Transverse extent | L⊥ = √Ω · D_M(z_c), small-angle square footprint | **486.329 Mpc** |
| LOS depth | L∥ = D_C(z_c + Δz/2) − D_C(z_c − Δz/2) | **315.598 Mpc** |
| Grid | ⌈L⊥ / 2.0⌉ = 244, snapped up to a power of two for the FFTs | `HII_DIM` = **256** |
| High-res grid | 3 × `HII_DIM`, the 21cmFAST convention | `DIM` = **768** |
| Solid angle | | 3.046 × 10⁻³ sr |
| LOS tiles | L∥ / L⊥ = 0.65 → no tiling needed | 0.649 |

### 2.2 Lightcone geometry (planned range)

| Quantity | Value | Source |
|---|---|---|
| z_min → z_max | **6.55 → 7.45** | config |
| z_obs (midpoint, used for noise/wedge/bias/RSD) | **7.0** | derived |
| D_C(6.55) / D_C(7.45) | 8647.13 / 8962.72 Mpc | **[C]** Planck18 |
| `L_los` | **315.598 Mpc** | **[C]** |
| `N_z` = round(L_los / cell), floored at 100 | **166** (floor idle) | **[C]** |
| LOS slice spacing (nominal) | **1.9012 Mpc**, against a 1.8997 Mpc transverse cell | **[C]** |
| LOS spacing, first → last cell | 2.0803 → 1.7599 Mpc, a **16.75 % spread** | **[C]** — this is `TODO.md` §P0.1, see §7.2 |
| Node redshifts | **9**, `linspace(7.45, 6.55, 9)`, step 0.1125 | **[C]** |
| z_obs falls exactly on node 5 | yes | **[C]** — §3a's `determine_halo_catalog(z=7.0)` therefore draws at a redshift the lightcone also evaluates; whether the C backend reuses that cached `HaloCatalog.h5` or recomputes it depends on the cache directory being shared. §4.5 budgets for it being recomputed |
| Observed frequency span | **168.095 – 188.133 MHz**, 20.038 MHz | **[C]** |
| Box shape written | (256, 256, 166) | **[C]** |
| `L_los` attribute recorded (`N_z` × transverse cell) | 315.354 Mpc against a true 315.598 | **[C]** — **0.08 %**, so `HPC.md` §11.1's factor-56.5 discrepancy disappears at this range. It is an artifact of `minimum_los_slices` flooring the 3.5 Mpc slab to 100 slices, not a bug in the geometry |

### 2.3 Resolution and mass

| Grid | Cell | Mass resolution **[C]** |
|---|---|---|
| `DIM` = 768 (initial conditions) | 0.6332 Mpc | **1.007 × 10¹⁰ M☉** |
| `HII_DIM` = 256 (ionisation, 21 cm) | 1.8997 Mpc | **2.720 × 10¹¹ M☉** |
| Halo sampler floor `SAMPLER_MIN_MASS` | — | **1 × 10⁸ M☉** (template default) |

Using ρ̄_m = Ω_m ρ_crit,0 = 3.967 × 10¹⁰ M☉ Mpc⁻³. The sampler resolves halos
100× below the IC cell mass via the conditional mass function.

### 2.4 Astrophysics — the `"simple"` template

`p21c.InputParameters.from_template(["simple"], random_seed=42)`, cloned with
the grid, `USE_INTERPOLATION_TABLES = "hmf-interpolation"` and
`MINIMIZE_MEMORY = True`. The values that set the cost and the physics:

| Parameter | Value | Consequence |
|---|---|---|
| `SOURCE_MODEL` | `CHMF-SAMPLER` | Discrete halo sampling — this is what makes the catalogue, and the cost, scale with volume |
| `SAMPLE_METHOD` | `MASS-LIMITED` | |
| `SAMPLER_MIN_MASS` | 1 × 10⁸ M☉ | Sets the halo count outright; see §7.1 |
| `M_TURN` | 10^8.7 = **5 × 10⁸ M☉** | Star-formation turnover; halos below it are exponentially suppressed |
| `F_STAR10` / `F_ESC10` | 10^−1.3 / 10^−1.0 | |
| `HII_EFF_FACTOR` | 30 | |
| `t_STAR` | 0.5 | → t_sf = t_★ t_H(7) = **570.3 Myr** |
| `USE_TS_FLUCT` | **False** | No spin-temperature evolution — the single biggest runtime saving in the template |
| `USE_MINI_HALOS` / `INHOMO_RECO` | False / False | |
| `PERTURB_ALGORITHM` | 2LPT | |
| `apply_rsds` | **False** | Kaiser applied analytically afterwards (`TODO.md` §P1.3) |
| `include_dvdr_in_tau21` | False | |
| `quantities` stored | `brightness_temp`, `density`, `neutral_fraction`, `halo_sfr` | |
| `random_seed` | **42** | Any other value is an independent realisation — the only route to cosmic-variance error bars, and it invalidates the whole cache |

### 2.5 Euclid survey, instrument, analysis

| Group | Parameter | Value |
|---|---|---|
| Selection | `M_UV_limit` (faint) / `M_UV_bright` | **−18** / **−22** → SFR ∈ **[0.1868, 7.4364]** M☉ yr⁻¹ at $\kappa_\mathrm{UV} = 2.7\times10^{-29}$ (Fisher et al. 2026) |
| | `mean_galaxy_density` n̄ | **7.48 × 10⁻⁵ Mpc⁻³** (plain Mpc⁻³, *not* h³ Mpc⁻³) → P_N,gal = 13 368.98 Mpc³ — Euclid Collab.: Allen et al. (2026), A&A 711, A25, Table 2; see `HPC.md` §11.10 |
| | `GALAXY_WEIGHTING` | `lightcone_sfr` (alternatives: `number`, `luminosity`) |
| | `galaxy_bias` | 8 — **fallback only**; the catalogue estimator overwrites it (4.744 at 256 Mpc **[M]**) |
| Cosmology | Ω_m, H₀, Ω_b | 0.315, 67.36, 0.049 (Planck 2018; distances use astropy `Planck18`) |
| HERA | dish, t_int, Δν | 14.0 m, 1000 h = 3.6 × 10⁶ s, 8 MHz |
| | T_sys at 177.55 MHz | 328.6 K → P_N,21 = 3.7488 mK² Mpc³ (default `scaling` noise model) |
| Wedge | `wedge_buffer` | 0.0677 Mpc⁻¹ (= 0.1 h Mpc⁻¹, Pober+2014 moderate); horizon slope 3.1509, FoV slope 0.37936 at z = 7 |
| Binning | `n_bins_perp` × `n_bins_parallel` | 20 × 20, log-spaced |
| Compute | `N_THREADS` | `N_THREADS` env → `SLURM_CPUS_PER_TASK` → `os.cpu_count()` → 1 |
| | `MINIMIZE_MEMORY` | **True** — trades peak RAM for intermediate I/O |

### 2.6 k-space the planned box delivers **[C]**

| Quantity | Planned (486.33 Mpc, L∥ = 315.6 Mpc) | Committed slab, for contrast |
|---|---|---|
| Δk⊥ = 2π/L⊥ | **0.01292 Mpc⁻¹** | 0.01292 |
| Δk∥ = 2π/L∥ | **0.01991 Mpc⁻¹** | 0.03142 (from the recorded `L_los`, §11.1 of `HPC.md`) |
| k_Nyq,⊥ / k_Nyq,∥ | 1.6537 / 1.6524 Mpc⁻¹ | 1.6537 / — |
| Fourier modes on the grid | 256 × 256 × 166 = 1.088 × 10⁷ | 1.638 × 10⁶ |

---

## 3. Outputs and their sizes **[E]**

Extrapolated from the 2026-08-12 run at 256 Mpc, which produced 114,289,081
halos in 2.74 GB.

| Path | Size | Note |
|---|---|---|
| `outputs/lightcone_data.h5` | **~19.1 GB** | 18.8 GB of it is the halo catalogue at 24 B/halo (masses, coords ×3, stellar masses, SFR — all float32) |
| — four lightcone fields | 174 MB | (256, 256, 166) float32, gzip-4 |
| — `galaxy_overdensity` | 87 MB | float64 |
| `outputs/analysis_products.h5` | ~54 KB | 3 × (20 × 20) spectra + mode counts + k-grids + the 10-map budget group |
| `outputs/figures/*.png` | ~10–40 MB | 18 figures at dpi 200; `halo_catalogue.png` scales with the number of plotted halos — use `--max-halos` |
| `outputs/pipeline_summary.json` | ~2 KB | |
| `outputs/runs/sim_<run_id>.json` | ~5 KB | The run manifest — parameters, derived geometry, cost estimate, per-stage timings, peak RAM |
| **21cmFAST cache** | **~245 GB** | See §4.5. Written to the *current working directory* unless redirected |

---

## 4. Compute, memory and storage

### 4.1 The target machine **[M]**

Reported 2026-08-24 by `lscpu` and `df -hT` on the machine the run will use.
It is a large shared workstation, **not a scheduled cluster**: `sinfo`,
`scontrol`, `sacct` and `sacctmgr` are all absent, so there is no queue, no
partition and no allocation to apply for — and no scheduler to stop two jobs
colliding.

| | Value | Consequence |
|---|---|---|
| CPU | **2 × AMD EPYC 9374F**, Zen 4, 32 cores each | 64 physical cores, 128 logical (SMT2). A frequency-optimised Genoa part — 3.85 GHz base, 4.30 GHz boost — so the measured single-core baseline needs **no speed correction** if it came from this same host |
| NUMA | 2 nodes: cores 0–31/64–95 and 32–63/96–127 | 21cmFAST's large FFTs are memory-bandwidth-bound; crossing sockets costs more than it gains unless memory is interleaved |
| RAM | **1.5 TiB**, 1.4 TiB available, 88 GiB in use by others **[M]** | The 56 GB peak is **1 in 27** of the machine. Memory is not a constraint here at any box size in §4.6. Swap is only 4 GiB, so an over-run would be an OOM kill rather than a thrash — irrelevant at this margin, but worth knowing |
| Local scratch | `/nvme1` **851 GB free**, `/nvme4` **836 GB free** (XFS on NVMe) | Where the 245 GB cache belongs. `/nvme2` (344 GB) and `/nvme3` (532 GB) are already 41–62 % used |
| Home | `/home/sdutta`, NFS, **144 GB free of 367 GB** | **Will not hold the cache** — and it is NFS, so writing 245 GB across it would dominate the runtime |
| `/dev/shm` | 756 GiB tmpfs | Tempting, and wrong: it is RAM |

Two consequences reshape the rest of this section, and both are in §7.3:
there are no `#SBATCH` directives to add (R1), and `N_THREADS` must be set by
hand, because `resolve_n_threads()` falls through to `os.cpu_count()` on a
machine with no `SLURM_CPUS_PER_TASK` — **which is 128 here**, i.e. every SMT
thread of a shared machine.

> **All confirmed 2026-08-24.** The host is `andromeda1.jb.man.ac.uk` — the
> same machine the 2026-08-12 baseline in §4.2 was measured on, so the serial
> scaling below carries **no cross-machine assumption at all**: same CPU, same
> filesystem, same environment. `free -h` gives 1.5 TiB total / 1.4 TiB
> available. `/nvme1` and `/nvme4` are writable and persistent.
>
> One assumption survives, and no survey can settle it: **how well 21cmFAST
> threads**. Everything else in §4.3 is measured.

### 4.2 The measured baseline **[M]**

Two independent end-to-end runs, both on `andromeda1.jb.man.ac.uk`, at
256 Mpc / 128³ / `DIM` 384 / 5 nodes / N_z = 100, effectively **single-threaded**
(`N_THREADS` was not set then; user/real = 1.04):

| Run | Stage 1 (simulation) | Total wall | User CPU | CPU-hours |
|---|---|---|---|---|
| 2026-08-07 | — | 544.6 s | 568.7 s | 0.1636 |
| 2026-08-12 | 520.1 s | 543.5 s | 567.4 s | 0.1632 |

Stage 2 (power spectra + Euclid selection + bias) took **0.8 s**; Stage 3
(10 figures) took ~20 s. **Stage 1 is 96 % of the cost.**

A third data point, less useful but worth recording: the 2026-08-20 attempt at
the current 486.33 Mpc box ran **2,303 s on one core** before dying with
SIGSEGV (§7.1).

### 4.3 Scaling to the planned run **[E]**

The halo sampler's floor is a fixed mass, not a grid property, so its cost
tracks comoving volume; the IC/density work tracks `DIM³`.

| Term | Factor vs baseline | Basis |
|---|---|---|
| Comoving volume (486.33/256)³ | **6.86×** | halo sampling, catalogue I/O |
| High-res grid (768/384)³ | **8×** | initial conditions, 2LPT |
| Node redshifts 9/5 | **1.8×** | per-node work repeats |
| Lightcone slices 166/100 | 1.66× | interpolation and output only |
| **Serial simulation cost** | **~12.3×** | 6.86 × 1.8, IC term folded in |

520 s × 12.3 ≈ **6,400 s ≈ 1.8 h of single-core compute**, plus initial
conditions at the larger `DIM` (~10 min) and the I/O below.

The baseline was measured on this same host, so that serial figure needs no
per-core correction — it is not an extrapolation across machines, only across
box sizes. The one open variable is how well 21cmFAST threads:

| Configuration | Expected wall **[E]** | Note |
|---|---|---|
| 1 thread (what the baseline used) | 2–3 h | the current committed behaviour if `N_THREADS` is left to chance and the machine reports 1 |
| **32 threads, one socket** **[T]** | **15–40 min** | assumes 8–14× on a memory-bandwidth-bound workload; `numactl --cpunodebind=0 --membind=0` keeps every FFT on local memory |
| 64 threads, both sockets | 10–30 min, or worse | needs `numactl --interleave=all`; measure before trusting it |
| 128 threads (`os.cpu_count()`, the accidental default) | **do not** | SMT siblings share one FP unit — bad for FFTs — and it takes the whole shared machine |

Budget **1 h**, with 2 h as the ceiling. There is no queue to satisfy and
nothing is billed, so the only cost of over-running is the machine's other
users.

The next run will settle this: `outputs/runs/sim_<run_id>.json` records
`timings_seconds` per stage and `peak_memory_GB`, whether or not the run
finishes.

I/O is now the cheap part: ~245 GB of cache to a local NVMe at ~1–2 GB/s is
**~2–4 min** **[E]**, against ~20 min on a contended shared filesystem — as
long as the cache lands on `/nvme1` or `/nvme4` and not on the NFS home
(R3). `MINIMIZE_MEMORY = True` deliberately adds more I/O, and on this machine
that is the right trade only if memory were tight; with 1.5 TiB it is not, so
**consider setting it back to `False`** and measuring.

### 4.4 Memory **[E]**

| Moment | Resident | What is held |
|---|---|---|
| Initial conditions | ~8–18 GB | 1.81 GB per `DIM³` float32 array; the written IC file is ~7.7 GB |
| `perturb_halo_catalog` — **the high-water mark** | **~48.2 GB** | Lagrangian catalogue (26.2 GB) and its perturbed copy (21.9 GB) simultaneously |
| Peak overall | **~56 GB** | the above plus the ICs, before `del halo_catalog, initial_conditions` frees them |
| After `del halo_catalog, initial_conditions` | ~22–25 GB | the perturbed catalogue, plus the one genuine NumPy copy the unit conversion makes (`sfr_cat = ... * _SEC_PER_YR`, 3.1 GB) |
| Stage 2/3 with the full catalogue | ~19 GB + float64 temporaries (~6.3 GB each) | mitigate with `--max-halos 5000000` |

**Not a constraint on this machine.** The measured 1.4 TiB available (§4.1)
swallows the 56 GB peak twenty-five times over, and the full 784 M-halo
catalogue can be loaded into the analysis stage without `--max-halos`. Note what this does *not* buy: the
2026-08-20 crash at this box size was an index overflow, not memory pressure
(§7.1), and abundant RAM does nothing for it. Keep the figure recorded anyway
— it is what a future SLURM `--mem` request would have to state.

### 4.5 Storage **[E]**

| Item | Size |
|---|---|
| `InitialConditions.h5` (`DIM` = 768) | 7.7 GB |
| Per node redshift (`HaloCatalog.h5` 26.3 GB + `PerturbedField.h5` 0.14 GB) | **26.4 GB** |
| × 9 nodes | 238 GB |
| **Cache total** | **~245 GB** |
| Output HDF5 + figures | ~19 GB |
| **Total per run** | **~265 GB** |
| **Where it goes** **[T]** | `/nvme1` (851 GB free) or `/nvme4` (836 GB free) |

Scaled from the measured cache at 256 Mpc (`InitialConditions.h5` 0.96 GB,
3.83 GB per node). The cache is the dominant footprint and the easiest thing
to under-budget; it is also **written to the working directory by default**,
so on this machine launching from `$HOME` puts 245 GB onto an NFS mount with
144 GB free — see R3 in §7.3.

Two runs' caches will not both fit on one NVMe with the outputs, so clear the
previous cache — or point the second run at the other disk — before
re-running with a different seed.

### 4.6 Cheaper variants, if the numbers above do not fit **[C]**

| Box | Area implied | Volume | Halos (Lagrangian) | Catalogue | Resident | INT_MAX headroom | Cache (9 nodes) |
|---|---|---|---|---|---|---|---|
| 256 Mpc (the old box, `DIM` = 384) | 2.8 deg² | 1.68 × 10⁷ | 1.37 × 10⁸ | 3.8 GB | 7.0 GB | 0.19 | 36 GB |
| **350 Mpc** | 5.2 deg² | 4.29 × 10⁷ | 3.49 × 10⁸ | 9.8 GB | 18.0 GB | **0.49** | 97 GB |
| 384 Mpc | 6.2 deg² | 5.66 × 10⁷ | 4.61 × 10⁸ | 12.9 GB | 23.7 GB | 0.64 | 125 GB |
| **486.33 Mpc (planned)** | 10 deg² | 1.15 × 10⁸ | 9.37 × 10⁸ | 26.2 GB | 48.2 GB | **1.31 ✗** | 245 GB |
| 1000 Mpc (`TODO.md` §P3.1) | 42 deg² | 1.00 × 10⁹ | 8.15 × 10⁹ | 228 GB | 419 GB | 11.4 ✗ | ~2.1 TB (`DIM` = 1500) |

The 350 Mpc row is the recommended first attempt at the new scale: 2.6× the
old volume, comfortably inside the 32-bit index, and it measures the scaling
before the full footprint is committed to.

**Shrinking the box saves the halo sampler, not the grid.** With the
power-of-two snap, every footprint between ~257 and 512 Mpc lands on the same
`HII_DIM` = 256 / `DIM` = 768 grid **[C]** — a 5.2 deg² field gives 350.02 Mpc
with a *finer* 1.367 Mpc cell (and a correspondingly lower 1.0 × 10¹¹ M☉
`M_cell`), not a cheaper one. The IC and FFT work is therefore unchanged
across those rows; what falls is the catalogue, which is where the memory,
the storage and the 32-bit risk all live. Pass
`snap_hii_dim_to_power_of_two=False` to `survey_area_to_box_size()` if the
2.0 Mpc target cell is wanted instead (350 Mpc → `HII_DIM` = 176,
`DIM` = 528, and the IC cost falls with it).

---

## 5. What the run will deliver

`lightcone_data.h5` with the four fields and the full catalogue; the
catalogue-measured Euclid bias b_g and β; 20 × 20 cylindrical P₂₁, P_gal and
P₂₁ₓgal; the complete uncertainty budget (photo-*z* kernel, wedge mask,
thermal and shot noise, per-bin and total SNR); 18 figures; a summary JSON;
and a run manifest recording the parameters, timings and peak memory.

Over the old box it buys: **6.9× the volume**, hence a far better-sampled
bright end of the UV luminosity function (~49 k Euclid-selected halos at
256 Mpc **[M]** → ~340 k **[E]**), a largest mode of Δk⊥ = 0.0129 Mpc⁻¹, and —
for the first time in this pipeline — genuine line-of-sight structure across
20 MHz rather than an interpolated 0.9 MHz slab.

## 6. What it will not deliver — state this up front

**It will not produce a detection, and no box size fixes that.** At
σ_z = 0.45 the radial smearing is σ_r = 157.48 Mpc **[C]**, so the photo-*z*
kernel W = exp(−k∥²σ_r²/2) evaluates to **[C]**:

| k∥ [Mpc⁻¹] | Where it comes from | W |
|---|---|---|
| 0.0199 | the largest LOS mode this box samples | 7.3 × 10⁻³ |
| 0.1084 | the first mode the wedge admits, at k⊥,min | 5.2 × 10⁻⁶⁴ |
| 0.00748 | where W = 0.5 | needs L_LOS > 840 Mpc to sample at all |

The wedge and the photo-*z* kernel have no overlap on this grid, which is why
the current pipeline reports a total SNR of ~10⁻¹¹¹ σ. That is physics, not a
bug: a photometric survey with σ_z = 0.45 at z = 7 retains essentially no
line-of-sight information. Recovering a forecast needs spectroscopic
redshifts, the angular (2D) cross-correlation, or a box long enough to reach
k∥ → 0 — not more cells. Budget this run as the simulation and estimator
infrastructure it is, and expect the detectability statement to come from a
change of observable.

---

## 7. Blockers, in the order they will bite

### 7.1 The halo catalogue overflows a 32-bit index — **[E]**, hard blocker

At 486.33 Mpc the flattened `halo_coords` array holds 2.81 × 10⁹ elements,
**1.31 × `INT_MAX`**. 21cmFAST's C backend indexes halo arrays with `int`.
This is what killed the 2026-08-20 run with SIGSEGV rather than SIGKILL, and
the supporting evidence is a cached `HaloCatalog.h5` written at this box size
on a 512× smaller grid — under no memory pressure whatsoever — that stops dead
at 2,147,491,839 bytes, the signed 32-bit boundary. **A bigger node does not
fix this.** `run_simulation.py` now prints the warning before spending any
compute.

Three ways out, with numbers:

| Fix | Halos | Catalogue | Resident | Headroom | Cost |
|---|---|---|---|---|---|
| **(a)** `SAMPLER_MIN_MASS` 1 × 10⁸ → **2 × 10⁸ M☉** | 4.33 × 10⁸ | 12.1 GB | 22.3 GB | **0.60** | Loses **0.16 % of the total star formation [C]** — halos below 2 × 10⁸ sit far under the M_TURN = 5 × 10⁸ turnover and are exponentially suppressed. Does change halo-count statistics and the SFR field's shot noise |
| **(b)** Shrink the box to **≤ 445 Mpc** (headroom 1.0), **350 Mpc** for margin | 3.49 × 10⁸ at 350 Mpc | 9.8 GB | 18.0 GB | 0.49 | Forecasts a 5.2 deg² sub-field rather than the full 10 deg² footprint |
| **(c)** Both | 1.61 × 10⁸ | 4.5 GB | 8.3 GB | 0.22 | Cheapest; furthest from the target |

**Recommendation [T]:** (b) at 350 Mpc for the first production run — it keeps
the physics untouched and its measured timings then calibrate the
extrapolations in §4 for the full footprint. Take (a) on top only if the full
10 deg² is required and the 0.16 % SFR loss is acceptable.

Note that the §4.1 machine's 1.5 TiB does **not** open route (d), "just run
the 486.33 Mpc box on a big node". The constraint is the width of a C `int`,
not the address space: the flattened `halo_coords` is 1.31 × `INT_MAX`
whatever the node has. Memory being abundant here removes a worry the
2026-08-20 post-mortem had to rule out — it does not remove the blocker.

The mass-floor fractions come from a Sheth-Tormen mass function at z = 7
weighted by 21cmFAST's own f_★ ∝ (M/10¹⁰)^0.5 exp(−M_TURN/M): 1.5 × 10⁸ keeps
63.1 % of halos and 99.95 % of the star formation; 2 × 10⁸ keeps 46.2 % and
99.84 %; 3 × 10⁸ keeps 28.7 % and 99.44 %.

### 7.2 The estimator is not lightcone-ready — `TODO.md` §P0, science blocker

Measured for the planned range **[C]**:

| # | Issue | At z = 6.55–7.45 |
|---|---|---|
| P0.1 | `lc_redshifts` is uniform in z, not in comoving distance | LOS cell runs 2.0803 → 1.7599 Mpc, a **16.75 % spread**, against the single 1.9127 Mpc the FFT assumes |
| P0.2 | A single global scalar mean is subtracted from δT_b and ⟨SFR⟩ | Leaves a monotonic LOS ramp that aliases into low-k∥ — exactly where the wedge analysis looks |
| P0.3 | One FFT over Δz = 0.9 returns a redshift-*averaged* spectrum | x_HI evolves substantially over this range; the effective redshift is ill-defined |
| P0.4 | PS bandwidth vs noise bandwidth | **20.04 MHz vs 8 MHz = 2.50× mismatch** — signal and noise are computed over different volumes |

P0.1 is a one-line change to
`RectilinearLightconer.between_redshifts(...)` and **must be done before the
production run**, because it changes what is written to disk. P0.2 is small
and localised. P0.3/P0.4 are one shared fix — sub-band the lightcone into
~2.5 chunks of 8 MHz, tapered, each reported at its own effective redshift.

### 7.3 Site requirements not yet in the repo

| # | Requirement | Status |
|---|---|---|
| **R1** | ~~`#SBATCH` directives~~ → run under `tmux`/`nohup` | **Moot on this machine: there is no scheduler** (§4.1). `submit_job.sh` running in the foreground is the correct mode here; what it needs instead is a session that survives a dropped SSH connection, and thread and NUMA binding. Template in §8 |
| **R2** | `21cmfast` conda env with py21cmfast 4.1.1 | Per [`INSTALL_21cmFASTv4.md`](INSTALL_21cmFASTv4.md). **If `py21cmfast` is missing the run does not fail — it silently falls back to synthetic fields** |
| **R3** | 21cmFAST cache on local NVMe | **Not passed.** `run_lightcone` defaults to `OutputCache(Path('.'))`, i.e. the working directory, and `~/.21cmfast/config.toml` does *not* govern it. On this machine `$HOME` is NFS with **144 GB free against a 245 GB cache** — the run will die partway. Launch from `/nvme1/<user>/21cm_run` (or `/nvme4`), or pass `cache=p21c.OutputCache("/nvme1/<user>/cache")` at `run_simulation.py:477` |
| **R4** | `outputs/` writable, ≥ 25 GB, path relative to cwd | |
| **R5** | `EMAIL_TO` in `submit_job.sh:26`, `sendmail` present | Currently hardcoded to the developer's address |
| **R6** | `PIP_CACHE_DIR`, `XDG_CACHE_HOME` on scratch *before* `pip install` | Otherwise the build dies with `Errno 122: Disk quota exceeded` |
| **R7** | `N_THREADS` | **Must be set by hand here.** `resolve_n_threads()` looks for `N_THREADS`, then `SLURM_CPUS_PER_TASK`, then falls through to `os.cpu_count()` — **128 on this machine**, every SMT thread of a shared box. Launch with `N_THREADS=32` |

---

## 8. How to launch it **[T]**

There is no scheduler, so the job is a foreground process that must outlive
the SSH session, stay on one socket, and keep its 265 GB off the NFS home.

```bash
tmux new -s 21cm                      # or: nohup ... &  — the run outlives the login

mkdir -p /nvme1/$USER/21cm_run        # R3/R4: cache and outputs are cwd-relative
cd /nvme1/$USER/21cm_run

N_THREADS=32 \
numactl --cpunodebind=0 --membind=0 \
  bash /path/to/repo/submit_job.sh --sim force
```

`N_THREADS=32` is read by `provenance.resolve_n_threads()` before anything
else, so it beats the `os.cpu_count()` fallback that would otherwise take all
128 logical CPUs. `numactl` keeps the FFT working set on one socket's memory
controllers; drop it and use `--interleave=all` with `N_THREADS=64` only after
measuring both. Everything lands under the current directory: the
`d1f8b93…/` cache, `outputs/lightcone_data.h5`, `outputs/figures/`, and the
run manifest.

Watch it from another pane with `tail -f outputs/21cm_pipeline_*.log`; the
child now runs under `python -u`, so progress reaches the log as it happens.

<details>
<summary>If the run later moves to a scheduled cluster (CSD3)</summary>

```bash
#!/bin/bash
#SBATCH --job-name=21cm_pipeline
#SBATCH --partition=<partition>
#SBATCH --account=<account>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # -> N_THREADS, via resolve_n_threads()
#SBATCH --mem=128G                  # peak ~56 GB; margin is deliberate
#SBATCH --time=06:00:00
#SBATCH --output=outputs/slurm-%j.out

cd /path/to/scratch/21cm_run
bash /path/to/repo/submit_job.sh --sim force
```

There `SLURM_CPUS_PER_TASK` sets the threads on its own and `--mem` matters;
the wall time is padded because a shared cluster's cores are usually slower
than this machine's 3.85 GHz Genoa parts.
</details>

Pre-flight, from `HPC.md` §13.7 plus this document:

```text
[ ] R3  launched from /nvme1 or /nvme4, >= 300 GB free there
[ ] R7  N_THREADS=32 exported (not left to os.cpu_count() = 128)
[ ] 7.1 pre-flight cost estimate prints NO INT_MAX warning
[ ] 7.2 P0.1 landed (np.diff(lc_dist_Mpc) constant to machine precision)
[ ] z_min/z_max edited to 6.55/7.45 and launched with --sim force
[ ] M_UV_bright (run_simulation.py:183) == --m-uv-bright
[ ] seed / template / SAMPLER_MIN_MASS change recorded in CHANGELOG.md
[ ] previous run's cache cleared, or this run pointed at the other NVMe
```

After any failure, read `outputs/runs/sim_<run_id>.json` first: `status` and
`stage` name where it died even when the log does not.

---

## 9. What is measured, and what is still open

The machine was surveyed on 2026-08-24 and is `andromeda1.jb.man.ac.uk` — the
same host every timing in §4.2 came from. Six things were needed to turn the
wall-time range into a number; five are now measured, and the sixth was
answered by its own absence.

| # | What was needed | Answer |
|---|---|---|
| 1 | CPU model and clock | **2 × AMD EPYC 9374F**, 3.85 GHz base / 4.30 boost, Zen 4 **[M]**. Same host as the baseline, so the per-core factor is exactly 1.0 |
| 2 | Cores available | **64 physical / 128 logical**; use 32 on one socket **[M]/[T]** |
| 3 | RAM | **1.5 TiB total, 1.4 TiB available** **[M]** |
| 4 | Scratch | Local NVMe, XFS: `/nvme1` 851 GB, `/nvme4` 836 GB free, both writable and persistent **[M]** |
| 5 | Partitions, wall limits, charging | **None — there is no scheduler.** `sinfo`, `scontrol`, `sacct`, `sacctmgr` are all absent **[M]**. Nothing to apply for; nothing billed; nothing stopping a collision with another user either |
| 6 | FFTW/GSL as site modules | Moot — conda-forge provides both per `env.yml`, and `module` is not available anyway |

### The one thing a survey cannot answer

**How well 21cmFAST threads on this workload.** §4.3 assumes 8–14× on 32
cores, which is a judgement about a memory-bandwidth-bound FFT and sampler,
not a measurement. It is the whole width of the 15–40 min range. The NVMe
write bandwidth (assumed ~1–2 GB/s) is a distant second, worth ±2 min.

Both fall out of one run.

### Measuring rather than predicting

There is no `sacct` here, so the run's own instrumentation is the record:

```bash
/usr/bin/time -v python run_pipeline.py --sim force     # Maximum resident set size
cat outputs/runs/sim_<run_id>.json                      # per-stage timings, peak RSS
htop -u $USER                                           # threads actually running
```

`submit_job.sh` already wraps the run in `/usr/bin/time -p` and reports
wall-clock and CPU-hours by email, and `RunManifest` records
`timings_seconds` per stage and `peak_memory_GB` even when the run dies.

**The calibration run to do first** is the committed smoke-test slab on this
machine — `N_THREADS=32 bash submit_job.sh --sim force` from `/nvme1`. It
costs a few minutes, converts the remaining **[E]** marks in §4 into **[M]**,
and exercises R3 and R7 before the production run depends on them. Compare its
`timings_seconds` against the 2026-08-12 baseline's 520 s and the speed-up
question is answered directly.

---

## 10. Reference — the planned run in one table

| Symbol | Value | Provenance |
|---|---|---|
| Survey / field | 10 deg², Euclid Deep Field Fornax, z_c = 7.0 | config |
| `BOX_LEN` / `HII_DIM` / `DIM` | 486.329 Mpc / 256 / 768 | derived from footprint |
| Cell (transverse / high-res) | 1.8997 / 0.6332 Mpc | **[C]** |
| Mass resolution (lo-res / hi-res / sampler floor) | 2.720 × 10¹¹ / 1.007 × 10¹⁰ / 1 × 10⁸ M☉ | **[C]** / template |
| z range, z_obs | 6.55 → 7.45, 7.0 | planned |
| L_LOS, N_z, LOS cell | 315.598 Mpc, 166, 1.9012 Mpc | **[C]** |
| Node redshifts | 9, 7.45 → 6.55 | **[C]** |
| Frequency span | 168.095 – 188.133 MHz (20.038 MHz) | **[C]** |
| Δk⊥ / Δk∥ / k_Nyq | 0.01292 / 0.01991 / 1.65 Mpc⁻¹ | **[C]** |
| Halos (Lagrangian / perturbed) | 9.37 × 10⁸ / 7.84 × 10⁸ | **[E]** |
| Euclid-selected halos | ~3.4 × 10⁵ | **[E]** |
| σ_z, σ_r | 0.45 (absolute), 157.48 Mpc | config / **[C]** |
| Wedge slope, buffer | 3.1509, 0.0677 Mpc⁻¹ | **[C]** |
| P_N,21 / P_N,gal | 3.7488 mK² Mpc³ / 333.33 Mpc³ | **[C]** |
| Output HDF5 | ~19.1 GB | **[E]** |
| 21cmFAST cache | ~245 GB | **[E]** |
| Peak RAM | ~56 GB | **[E]** |
| Wall time | 15–40 min at 32 threads; budget 1 h | **[E]** / **[T]** |
| INT_MAX headroom | **1.31 — blocking** | **[E]** |

---

## 11. References

Davies, J. et al. (2025), arXiv:2504.17254 — 21cmFASTv4 discrete source model ·
La Plante, P. et al. (2023), arXiv:2205.09770 — wedge geometry and
cross-spectrum variance · Pober, J. et al. (2014), arXiv:1310.7031 —
"moderate" foreground model, 0.1 h Mpc⁻¹ buffer · DeBoer, D. et al. (2017),
PASP 129, 045001 — HERA array · Park, J. et al. (2019), MNRAS 484, 933 —
t_★ t_H star-formation prescription · Sheth & Tormen (1999), MNRAS 308, 119 ·
Madau & Dickinson (2014), ARA&A 52, 415 · Euclid Collaboration (2022),
arXiv:2108.01201. Full bibliography in [`reference.md`](reference.md).
