# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770) and [Davies et al. (2025)](https://arxiv.org/abs/2504.17254).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Euclid/Roman Space Telescope.

The repository holds three things:

| | What | Entry point |
|---|---|---|
| **HPC pipeline** | The production path: 21cmFAST lightcone → galaxy field → power spectra → uncertainty budget → figures, driven by one command | `run_pipeline.py`, `submit_job.sh` |
| **Notebooks** | The same science interactively — an analytical-only forecast, a monolithic lightcone notebook, and the pipeline's Parts 2 and 3 | `21cm_galaxy_cross_uncertainty.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `notebooks/` |
| **`src/`** | The shared implementation both of the above import, so they cannot drift apart | `src/analysis.py`, `src/figures.py`, `src/dataio.py`, `src/conversions.py`, `src/foregrounds.py`, `src/provenance.py` |

Detail that used to sit in this file — per-notebook structure, equations and
fiducial parameters, the figure-by-figure literature references, the 21cmFAST
v4 API notes, the `src/` function reference, and the full bibliography — now
lives in [`docs/reference.md`](docs/reference.md).

---

## Quickstart

### 0. Complete, physical run:

```
conda run -n 21cmfast python run_pipeline.py \
  --sim force \
  --analysis force \
  --estimator auto \
  --noise-model physical \
  --mode-weighted \
  --plots all \
  --format pdf --dpi 300
```

### 1. Install

```bash
conda env create -f env.yml
conda activate 21cmfast
python -c "import py21cmfast; print(py21cmfast.__version__)"   # → 4.1.1
```

`21cmFAST` compiles C extensions against FFTW and GSL, which is why those come
from conda-forge *before* the pip step — `env.yml` orders this correctly. To
reproduce the exact environment of the last run, use the pinned freeze instead:
`pip install -r requirements.txt`.

**On HPC systems** follow [`docs/INSTALL_21cmFASTv4.md`](docs/INSTALL_21cmFASTv4.md)
instead — it covers the CSD3-specific failures `conda env create` will not solve
on its own (home-quota exhaustion during the pip build, `conda-libmamba-solver`
entry-point errors, FFTW linking).

Every command below assumes the `21cmfast` environment is active, or is run
through `conda run -n 21cmfast <command>`.

### 2. Run the HPC pipeline

```bash
python run_pipeline.py --smoke-test      # pre-flight: every stage, tiny box, seconds
python run_pipeline.py                   # everything; simulation runs only if the HDF5 is missing
python run_pipeline.py --sim force       # re-run the 21cmFAST simulation first
python run_pipeline.py --analysis force  # recompute the power spectra from stored fields
python run_pipeline.py --plots power snr # only the k-space figures
python run_pipeline.py --sigma-z 0.05    # sweep a survey parameter without re-simulating
python run_pipeline.py --estimator lightcone   # per-sub-band spectra (TODO.md P0)
bash submit_job.sh --sim force           # same, with timing + email notification
```

> **Part 1 must run first.** `outputs/` is gitignored, so a fresh clone has no
> `lightcone_data.h5` and every downstream stage will fail at the loading step
> until a simulation has produced it. The analysis and figure stages never
> import or re-run 21cmFAST — they read the HDF5 only — but that does not mean
> the file ships with the repository. (`resources/`, the reference PDFs, is
> likewise local-only.)

[`PIPELINE.md`](PIPELINE.md) has the stage-by-stage breakdown — which file runs
where, what it consumes, and what it writes — plus a flowchart and the
`lightcone_data.h5` schema.

**Stage control** — each stage runs fresh or from stored results, so the
expensive 21cmFAST run happens only when it must:

| Flag | `auto` (default) | `force` | `skip` |
|------|------------------|---------|--------|
| `--sim` | run only if `outputs/lightcone_data.h5` is missing | always re-run `run_simulation.py` | never run; error if there is no stored output |
| `--analysis` | recompute the power spectra only if the cache is missing or older than the simulation | always recompute | load `outputs/analysis_products.h5`; error if absent |

**Figure groups** — `--plots` takes any number of these (default `all`; `none`
skips plotting entirely):

| Group | Figures written |
|-------|-----------------|
| `fields` | `lightcone_fields`, `lightcone_slice` |
| `halos` | `halo_catalogue`, `sfr_relations` |
| `scaling` | `uv_luminosity_function`, `stellar_mass_muv`, `main_sequence`, `uv_selection_maps` |
| `euclid` | `euclid_selected_catalogue`, `selected_galaxy_overdensity`, `galaxy_overdensity_on_21cm` |
| `power` | `power_spectra_2d`, `galaxy_wedge`, `wedge_real_space` |
| `snr` | `cross_snr` |
| `budget` | `uncertainty_budget`, `photoz_suppression` |
| `bias` | `galaxy_bias` |

The `halos`, `scaling`, `euclid` and `bias` groups are skipped with a message
when the stored HDF5 carries no halo catalogue or bias estimate.

The `euclid` group is the post-selection view: the galaxies, halo masses and
SFRs that survive the $M_\mathrm{UV}$ window, the galaxy overdensity field
rebuilt from that selected catalogue alone, and that field overlaid on the
21 cm field. It rebuilds $\delta_\mathrm{gal}$ rather than reading the stored
one, because `run_simulation.py`'s default `GALAXY_WEIGHTING = "lightcone_sfr"`
builds it from the lightcone `halo_sfr` field and applies no magnitude cut at
all. Use `--galaxy-weighting luminosity` to weight the rebuild by
$L_\mathrm{UV}$ instead of galaxy counts.

**Paths and rendering:**

| Option | Default | Meaning |
|--------|---------|---------|
| `--data PATH` | `outputs/lightcone_data.h5` | simulation HDF5 |
| `--products PATH` | `outputs/analysis_products.h5` | analysis cache |
| `--figdir PATH` | `outputs/figures` | figure directory |
| `--summary PATH` | `outputs/pipeline_summary.json` | summary JSON |
| `--sim-script PATH` | `run_simulation.py` | simulation script to invoke |
| `--format {png,pdf,svg}` | `png` | figure file format |
| `--dpi N` | `200` | figure resolution |
| `--max-halos N` | `0` (all) | cap on halos loaded, uniformly strided; number densities are rescaled automatically. Lower it when memory is tight |
| `--m-uv-bright M` | `-22` | bright-end Euclid magnitude cut |
| `--galaxy-weighting {number,luminosity}` | `number` | per-halo weight used to rebuild the post-cut $\delta_\mathrm{gal}$ for the `euclid` figure group |
| `--quiet` | off | suppress progress output |

**Uncertainty-budget overrides** — none of these touches the simulated fields,
so they can be swept without `--sim force`. Each defaults to the corresponding
attribute of the stored HDF5:

| Option | Stored default | Meaning |
|--------|----------------|---------|
| `--noise-model {scaling,physical}` | `scaling` | 21 cm thermal noise. `physical` uses Parsons (2017) Eq. 12 / La Plante Eq. 11, resolved in $k_\perp$ through the HERA baseline distribution — ~10³ larger than the default flat estimate, and infinite where no baseline samples the mode |
| `--mode-weighted` | off | Apply La Plante Eq. 19's $\sqrt{N_\mathrm{patch}\,dN}$ weighting when summing bins, using the estimator's own `mode_counts`. Raises the total SNR ~10× |
| `--sigma-z σ` | `0.45` | absolute photo-$z$ uncertainty, **not** $\sigma_z/(1+z)$ — the Euclid requirement at $z = 7$ |
| `--wedge-buffer k` | `0.0677` Mpc⁻¹ | foreground-wedge margin — $0.1\ h\ \mathrm{Mpc}^{-1}$, Pober et al. (2014) "moderate" |
| `--integration-time s` | `3.6e6` (1000 h) | HERA integration time |
| `--bandwidth Hz` | `8e6` (8 MHz) | per-band bandwidth |

`python run_pipeline.py --help` prints the same list.

**`submit_job.sh`** forwards every argument verbatim to `run_pipeline.py` and
adds timing, a log file, and a `sendmail` report. It reads four environment
variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `JOB_NAME` | `21cm_pipeline` | log filename stem and email subject |
| `PYTHON_SCRIPT` | `run_pipeline.py` | set to `run_simulation.py` for simulation-only behaviour |
| `CONDA_ENV` | `21cmfast` | environment to activate |
| — | `EMAIL_TO` is set inside the script | notification recipient |

> It is a plain shell wrapper with **no `#SBATCH` directives**, so it runs in the
> foreground on whatever node invokes it. To submit it with `sbatch`, add the
> directives your cluster needs (`--partition`, `--time`, `--account`,
> `--cpus-per-task`, …) at the top first.

**Simulation parameters have no CLI.** `run_simulation.py` takes no arguments —
everything is a constant in its configuration block (lines ~95–255), and
changing any of them needs `--sim force` to take effect:

| Constant | Default | Meaning |
|----------|---------|---------|
| `SURVEY_AREA_DEG2` | `10.0` | Euclid Deep Field Fornax footprint [deg²] |
| `SURVEY_Z_CENTRAL` | `7.0` | central redshift of the analysis |
| `photoz_uncertainty` | `0.45` | absolute $\sigma_z$; also sets the LOS depth of the box |
| `PHOTOZ_N_SIGMA` | `1` | how many $\sigma_z$ the box spans — 1 → $\Delta z = 0.90$, 2 → $1.80$ |
| `z_min`, `z_max` | `6.995`, `7.005` | lightcone range — the deliberate smoke-test slab (see below) |
| `minimum_los_slices` | `100` | floor on $N_z$, overriding the cell-size-matched value |
| `M_UV_limit` / `M_UV_bright` | `-18` / `-22` | Euclid magnitude window |
| `mean_galaxy_density` | `3e-3` | $\bar n_\mathrm{gal}$ [$h^3$ Mpc⁻³] |
| `GALAXY_WEIGHTING` | `"lightcone_sfr"` | `lightcone_sfr` \| `number` \| `luminosity` — how $\delta_\mathrm{gal}$ is built |
| `ESTIMATOR` | `"coeval"` | `coeval` \| `lightcone` — which power-spectrum formalism the run is built for (see below) |
| `LIGHTCONE_SAMPLING` | derived → `"redshift"` | `redshift` \| `comoving`; `comoving` is `TODO.md` P0.1 |
| `GALAXY_MEAN_SUBTRACTION` | derived → `"global"` | `global` \| `per_slice`; `per_slice` is `TODO.md` P0.2 |
| `N_THREADS` | `N_THREADS` env → `SLURM_CPUS_PER_TASK` → `os.cpu_count()` | OpenMP threads for 21cmFAST. Its own default is **1**, which is why the 2026-08-20 run spent 38 minutes on a single core |
| `MINIMIZE_MEMORY` | `True` | Trades peak RAM for intermediate I/O in the C backend |
| `RANDOM_SEED` | `42` | 21cmFAST initial-conditions seed; recorded in the manifest and the HDF5 attrs |
| `galaxy_bias` | `8` | fallback only; overwritten by the halo-catalogue estimate (≈ 4.7) |
| `OMEGA_M_0`, `HUBBLE_CONSTANT` | `0.315`, `67.36` | Planck 2018 |
| `HERA_DISH_DIAMETER` | `14.0` | m |
| `integration_time`, `bandwidth` | `3.6e6` s, `8e6` Hz | 1000 h, 8 MHz |
| `wedge_buffer` | `0.0677` | Mpc⁻¹ |
| `n_bins_perp`, `n_bins_parallel` | `20`, `20` | $(k_\perp, k_\parallel)$ binning |

`BOX_LEN`, `HII_DIM` and `DIM` are **derived**, not hardcoded:
`survey_area_to_box_size()` turns the footprint into 486.33 Mpc / 256 / 768 at a
2.0 Mpc target cell size. [`docs/HPC.md`](docs/HPC.md) §13 is the checklist of
what must be set on a new cluster and which edits need `--sim force`;
[`docs/simulation_spec.md`](docs/simulation_spec.md) is the same configuration
projected forward to the production run — every derived number, the compute,
memory and storage it needs, and the two blockers in front of it.

> **The committed $z$ range is a deliberate thin slab**, $\Delta z = 0.01$
> ($L_\mathrm{LOS} = 3.5$ Mpc at $z = 7$). The power-spectrum estimator in
> `src/analysis.py` assumes statistical homogeneity along the LOS, which only
> holds for a quasi-coeval box, so configuration and formalism currently match.
> Widening to a true lightcone requires the estimator work in
> [`TODO.md`](TODO.md) §P0 — **now implemented, behind the `ESTIMATOR`
> toggle**, though the range itself has not been widened yet (P0.5).
> Treat the resulting SNR as a smoke-test
> number, not a forecast.

### 3. Run the notebooks

```bash
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb    # analytical forecast, no 21cmFAST
jupyter notebook 21cmfast_HERAxEuclid_lightcone.ipynb   # monolithic lightcone pipeline
jupyter notebook notebooks/plot_fields.ipynb            # Part 2 — field & catalogue figures
jupyter notebook notebooks/analysis.ipynb               # Part 3 — power spectra & SNR
```

Run all cells sequentially; each notebook is self-contained and renders its
figures inline.

| Notebook | Needs a simulation? | What to edit |
|----------|--------------------|--------------|
| `21cm_galaxy_cross_uncertainty.ipynb` | no — `numpy`/`matplotlib` only | parameters are set inline per section |
| `21cmfast_HERAxEuclid_lightcone.ipynb` | runs one itself | the **★ CONFIGURATION** cell (cell 4) — nothing else |
| `notebooks/plot_fields.ipynb` | yes | `OUTPUT_FILE` only (default `../outputs/lightcone_data.h5`) |
| `notebooks/analysis.ipynb` | yes | `OUTPUT_FILE` only; all other parameters are read from the HDF5 attributes |

**`21cmfast_HERAxEuclid_lightcone.ipynb` — the ★ CONFIGURATION cell.** These are
the notebook's equivalents of `run_simulation.py`'s constants, deliberately kept
in step with it so both describe the same experiment:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `SURVEY_AREA_DEG2`, `SURVEY_Z_CENTRAL` | `10.0`, `7.0` | survey footprint → box geometry |
| `photoz_uncertainty` | `0.45` | absolute $\sigma_z$ |
| `PHOTOZ_N_SIGMA` | `1` | box spans ±1σ → $\Delta z = 0.90$ |
| `HII_DIM`, `BOX_LEN`, `DIM` | derived → `256`, `486.33` Mpc, `768` | from `survey_area_to_box_size()` — do not hardcode |
| `z_min`, `z_max` | derived → `6.55`, `7.45` | **unlike `run_simulation.py`, the notebook uses the full footprint-implied range, not the thin slab** |
| `minimum_los_slices` | `100` | floor on $N_z$ (set in §1, not the config cell) |
| `M_UV_faint`, `M_UV_bright` | `-18.0`, `-22.66` | Euclid magnitude window |
| `mean_galaxy_density` | `3e-3` | $\bar n_\mathrm{gal}$ [$h^3$ Mpc⁻³] |
| `galaxy_bias` | commented out | computed in §3 from an HMF integral (≈ 5.39); requires `hmf` |
| `OMEGA_M_0`, `HUBBLE_CONSTANT` | `0.315`, `67.36` | literal Planck 2018 values, matching the HDF5 attributes |
| `HERA_DISH_DIAMETER` | `14.0` | m |
| `integration_time`, `bandwidth` | `3.6e6` s, `8e6` Hz | 1000 h, 8 MHz |
| `wedge_buffer` | `0.0677` | Mpc⁻¹ |
| `n_bins_perp`, `n_bins_parallel` | `20`, `20` | $(k_\perp, k_\parallel)$ binning |

**§7e — foreground contamination and removal.** `src/foregrounds.py` injects
a synthetic diffuse Galactic synchrotron foreground (plus an optional
point-source component) into `brightness_temp_field` before any power spectrum
is computed, then removes a controllable fraction of it. Each removal level is
pushed through `compute_all_power_spectra` and `compute_uncertainty_budget`
unchanged. Set `FOREGROUND_AMPLITUDE` (foreground RMS as a multiple of the
signal RMS; real foregrounds are $10^4$–$10^5$) and `REMOVAL_FRACTIONS` in the
§7e.1 cell.

Two caveats the section states in full, both worth knowing before quoting a
number from it:

- **The removal step is a placeholder, not an algorithm.** It subtracts an
  exactly-correct template of the field that was injected. It is a knob for
  *"what if removal were this good?"* — not GMCA, PCA, polynomial fitting or
  delay filtering, and it has none of their failure modes.
- **A contaminated SNR flatters the result.** Foregrounds are unbiased in the
  ensemble mean, but one realisation carries a chance cross-correlation that
  grows *linearly* with foreground amplitude, so `|P_×|/σ_×` has a
  contaminated numerator as well as denominator. §7e.3 plots a *signal-only*
  SNR beside the as-measured one; the gap is a spurious detection.

**§3c and §5c — after the Euclid cut.** The lightcone notebook mirrors the
`euclid` figure group: §3c plots the galaxies, halo masses and SFRs that
survive the magnitude window, then rebuilds $\delta_\mathrm{gal}$ from that
selected catalogue with `galaxy_overdensity_from_catalogue()`; §5c overlays
that field on the 21 cm field and reports the cell-by-cell correlation
coefficient. Set `GALAXY_WEIGHTING_DIAGNOSTIC` in the §3c.2 cell to switch
between number- and $L_\mathrm{UV}$-weighting.

Notebook figures use `%matplotlib widget`, which needs `ipywidgets >= 8` in the
JupyterLab front-end environment as well as the kernel; see
[`docs/reference.md`](docs/reference.md#figure-display-in-notebooks) if figures
come up blank.

### 4. Outputs

| Path | Contents |
|------|----------|
| `outputs/lightcone_data.h5` | Simulation fields, halo catalogue, and metadata (Part 1) |
| `outputs/analysis_products.h5` | Cached $P_{21}$, $P_\mathrm{gal}$, $P_{21\times\mathrm{gal}}$ and the $k$-grid, plus the `uncertainty_budget` group (damped spectra, wedge mask, $\sigma$ terms, per-mode SNR) |
| `outputs/figures/*.png` | The 18 figures listed in the figure-group table above |
| `outputs/pipeline_summary.json` | Scalar results: $\langle x_\mathrm{HI}\rangle$, wedge slopes, $\sigma_r$, total SNR, $\langle b_g\rangle$, selection counts. Overwritten every run |
| `outputs/runs/sim_<run_id>.json` | One manifest per simulation run: every configuration parameter, the derived geometry, the pre-flight cost estimate, code revision and package versions, per-stage timings, and peak memory. **Never overwritten** — see below |

**Run manifests.** `run_simulation.py` writes
`outputs/runs/sim_<run_id>.json` *before* it starts the expensive stages and
rewrites it after each one. A run killed by a signal cannot flush stdout or
run an exit hook, so the manifest is what survives: it is left with
`"status": "running"` and `"stage"` naming where the run died. `submit_job.sh`
summarises the newest manifest in its email, and the HDF5 records `run_id` /
`run_manifest` attributes so a stored dataset names the run that produced it.

Each simulation also prints a pre-flight halo-catalogue estimate before doing
any work, extrapolated from the measured 256 Mpc run. The catalogue scales
with comoving *volume*, so a modest-looking change in `BOX_LEN` is a large
change in cost — and if the flattened `halo_coords` would exceed `INT_MAX`,
the run says so and explains that more memory will not help:

```
Est. halos  : 9.370e+08 in 1.150e+08 Mpc³  →  26.2 GB on disk, ~52.5 GB resident
  *** WARNING: halo_coords would hold 2.811e+09 elements, 1.31x INT_MAX ***
```

A full run on the fiducial $128^2 \times 100$ lightcone (114 M halos, 2.8 GB
HDF5) takes ~35 s on a laptop when the simulation itself is skipped: ~1 s for
the three power spectra, the rest dominated by reading the halo catalogue and
rendering the catalogue-based figures. Use `--plots power snr` to skip the
catalogue entirely (it is not loaded when no figure needs it), or `--max-halos`
to cap the memory footprint.

---

## Repository structure

```
21cm_cross_correlation/
├── 21cm_galaxy_cross_uncertainty.ipynb    # Notebook 1 — analytical framework
├── 21cmfast_HERAxEuclid_lightcone.ipynb   # Notebook 3 — monolithic lightcone
├── run_pipeline.py                        # Pipeline driver — sim + analysis + figures
├── run_simulation.py                      # HPC pipeline Part 1 — simulation
├── submit_job.sh                          # Launcher + email notification
├── notebooks/
│   ├── plot_fields.ipynb                  # Part 2 — field & catalogue figures
│   └── analysis.ipynb                     # Part 3 — power spectra & SNR
├── src/
│   ├── analysis.py                        # Power spectra, wedge, noise, SNR, bias
│   ├── figures.py                         # All pipeline figures (headless)
│   ├── dataio.py                          # HDF5 loading + analysis-product cache
│   ├── conversions.py                     # Cosmological conversion utilities
│   └── FOV_to_cMpc.py                     # Survey geometry CLI
├── tests/                                 # pytest suite for src/ and the pipeline
├── docs/                                  # Reference & methodology notes
├── _archive/
│   └── 21cmfast_HERAxEuclid.ipynb         # Notebook 2 — coeval (superseded)
├── outputs/                               # Simulation products (gitignored)
├── resources/                             # Reference PDFs (gitignored)
├── env.yml                                # Conda environment specification
├── requirements.txt                       # Pinned pip freeze of the env
├── PIPELINE.md                            # HPC pipeline summary + flowchart
├── CHANGELOG.md                           # Record of all project changes
└── README.md                              # This file
```

## Smoke test — pre-flight before an HPC job

```bash
python run_pipeline.py --smoke-test          # whole pipeline, tiny box, seconds
python run_pipeline.py --smoke-test --plots none   # skip the figures
python run_simulation.py --smoke-test        # simulation stage only
```

Runs every stage — halo catalogue → 21 cm lightcone → galaxy field → bias →
Kaiser RSD → HDF5 → power spectra → photo-z/wedge/noise/SNR → figures →
summary — on a deliberately tiny configuration, and **asserts the shape of
each stage's output** rather than just checking it did not crash. It prints a
PASS/FAIL line per stage and exits non-zero if any check fails.

> **Not a science run.** Its numbers mean nothing: a 32 Mpc box does not
> sample the modes the forecast needs. Use it to prove the plumbing works
> before spending cluster time.

The reduced values live in [`src/smoke_test.py`](src/smoke_test.py) as
`SMOKE_TEST_OVERRIDES`, a module imported **only** when the flag is set — no
production default is edited, shadowed or reassigned without it, and outputs
go to `outputs/smoke_test/` so a real run's products can never be overwritten.
`SMOKE_TEST_UNCHANGED` records what was deliberately left alone and why.

## The estimator toggle (`TODO.md` P0)

The power-spectrum estimator was inherited from a *coeval* notebook and
assumes the box is statistically homogeneous along the line of sight. That
holds for the committed quasi-coeval slab and fails for a true lightcone, so
both formalisms now exist:

| | `"coeval"` (default) | `"lightcone"` |
|---|---|---|
| LOS sampling | uniform in redshift | uniform in comoving distance (P0.1) |
| Mean subtraction | one global scalar | per-slice, both fields (P0.2) |
| Transform | one FFT over the box | one per sub-band, Blackman-Harris tapered (P0.3) |
| Noise bandwidth | the configured 8 MHz | each band's own frequency span (P0.4) |
| Result | one spectrum, one SNR | one spectrum and budget per band, combined in quadrature |

`ESTIMATOR = "coeval"` reproduces every number this pipeline has produced,
bit for bit; everything P0 added is opt-in. The analysis stage reads the
`estimator` attribute back from the HDF5 and follows it, so the simulation and
the analysis cannot silently disagree — `--estimator {auto,coeval,lightcone}`
overrides that for a re-analysis from cached spectra, and
`--subband-bandwidth` sets the band width (default: the noise bandwidth,
which is the point of P0.4).

Under `lightcone`, per-band results reach `pipeline_summary.json` under
`subbands` and are printed as a table. **P0.5 — widening the redshift range
itself — is still open**; see [`TODO.md`](TODO.md).

## Documentation

| Document | Contents |
|----------|----------|
| [`docs/reference.md`](docs/reference.md) | **The long-form companion to this README** — notebook structure, equations and fiducial parameters; the literature relations overlaid on each figure; the 21cmFAST v4 `HaloBox`/lightcone API notes and source-model templates; the `src/` function reference; and the full bibliography |
| [`PIPELINE.md`](PIPELINE.md) | HPC pipeline summary, Mermaid flowchart, stage table, and the `lightcone_data.h5` schema |
| [`src/smoke_test.py`](src/smoke_test.py) | **The pre-flight smoke test** — the reduced configuration, what each override replaces and why, what was deliberately left alone, and the per-stage shape checks |
| [`docs/simulation_spec.md`](docs/simulation_spec.md) | **The planned run's specification and its cost** — every parameter and derived number for the production lightcone ($z = 6.55$–$7.45$ in the footprint-derived box), the compute / memory / scratch / wall-time requirements with their measured baselines, the two blockers that stand in front of it, a SLURM template, and what to look up about a new cluster to firm up the wall-time estimate |
| [`docs/HPC.md`](docs/HPC.md) | **Parameter-level ground truth for the HPC run** — every configuration value, derived quantity, formula with its evaluated number at $z=7$, file written, disk footprint, and known inconsistency. **§13 is the user-defined parameter checklist** |
| [`docs/uncertainty_budget.md`](docs/uncertainty_budget.md) | **The uncertainty budget, end to end** — every formula of the photo-$z$ / wedge / noise / SNR chain with its evaluated number at $z=7$, the audit against the source notebook, the CLI overrides, the HDF5 schema, and what the calculation still does not do |
| [`CHANGELOG.md`](CHANGELOG.md) | Chronological record of all changes, including corrected literature values |
| [`TODO.md`](TODO.md) | Outstanding work, priority-ordered — **including the lightcone power-spectrum corrections a wider $\Delta z$ would require** |
| [`docs/INSTALL_21cmFASTv4.md`](docs/INSTALL_21cmFASTv4.md) | Step-by-step 21cmFAST v4.1.1 install on CSD3/HPC, plus fixes for quota, conda-plugin, and FFTW linking failures |
| [`docs/Galaxy_bias_formalism.md`](docs/Galaxy_bias_formalism.md) | Methodology for the effective linear galaxy bias from a 21cmFAST halo catalogue under a Euclid $M_\mathrm{UV}$ cut |
| [`docs/halo_catalogue_reference.md`](docs/halo_catalogue_reference.md) | Field-by-field reference for the v4 halo catalogue, verified against py21cmfast v4.x |
| [`docs/Low_SFR_fix.md`](docs/Low_SFR_fix.md) | Diagnosis of the ~7.5 dex SFR offset in the main-sequence plot (the M☉ s⁻¹ unit bug) |
| [`docs/project_update.md`](docs/project_update.md) | Latest run's simulation parameters and numerical results |

## Testing

Project convention (see `CLAUDE.md`): every function in `src/` should have at
least one corresponding test, in `tests/test_<module>.py`.

```bash
conda run -n 21cmfast pytest tests/ -v
```

**Current status: 146 tests, all passing in ~22 s.** No test invokes 21cmFAST —
`tests/conftest.py` writes a synthetic `lightcone_data.h5` with the same schema
(16² × 12 cells, 4 000 halos), so the whole suite runs offline. What each file
exercises is tabulated in
[`docs/reference.md`](docs/reference.md#test-suite-coverage).

> `src/FOV_to_cMpc.py` and the magnitude/SFR half of `src/conversions.py` are
> still untested directly, though the conversion round-trips are exercised
> indirectly through the selection and UVLF tests. The remaining first
> candidates are the explicit identities `Muv_to_Luv` ↔ `Luv_to_Muv` and
> `sfr_to_Luv` ↔ `Luv_to_sfr`.

## Key references

- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — uncertainty equations, HERA noise model, and foreground wedge prescription
- **Davies, Mesinger & Murray (2025)** — [arXiv:2504.17254](https://arxiv.org/abs/2504.17254) — 21cmFASTv4 discrete source model
- **Gagnon-Hartman, Davies & Mesinger (2025)** — [arXiv:2502.20447](https://arxiv.org/abs/2502.20447) — galaxy–21 cm cross-correlation detection forecasts

The full bibliography — every relation, calibration and instrument
specification used anywhere in the project — is in
[`docs/reference.md`](docs/reference.md#bibliography).
