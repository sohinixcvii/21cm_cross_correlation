# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased]

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
