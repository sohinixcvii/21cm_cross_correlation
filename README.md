# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770) and [Davies et al. (2025)](https://arxiv.org/abs/2504.17254).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Euclid/Roman Space Telescope.

This project contains two active notebooks plus a refactored HPC-optimised lightcone pipeline, with one superseded notebook retained in `_archive/`:

1. **`21cm_galaxy_cross_uncertainty.ipynb`** — analytical framework using semi-analytic signal models and the La Plante et al. (2023) variance estimators.
2. **`_archive/21cmfast_HERAxEuclid.ipynb`** *(archived)* — end-to-end simulation pipeline using 21cmFASTv4 with the discrete source model (Davies et al. 2025), based on a **coeval** (single-snapshot) simulation at $z = 6.5$. Superseded by (3); retained as the coeval reference implementation.
3. **`21cmfast_HERAxEuclid_lightcone.ipynb`** — same pipeline as (2) but using a **lightcone** simulation spanning a redshift range, producing a non-cubic $(N_\perp \times N_\perp \times N_z)$ volume with continuous redshift evolution along the LOS.

### HPC lightcone pipeline (recommended for cluster use)

See [`PIPELINE.md`](PIPELINE.md) for a succinct pipeline summary with a flowchart.

The lightcone workflow has been split into self-contained parts for efficient use on HPC clusters, with `run_pipeline.py` driving all of them from one command:

| File | Purpose |
|------|---------|
| `run_pipeline.py` | **Driver** — runs the whole workflow: optional simulation, full analysis (fresh or from stored results), all figures, and a JSON summary |
| `run_simulation.py` | **Part 1** — runs the 21cmFAST lightcone, constructs the galaxy field, estimates galaxy bias, applies Kaiser RSD, and saves all outputs to `outputs/lightcone_data.h5` |
| `src/analysis.py` | **Part 3 computation** — cylindrical power spectra, wedge geometry, photo-$z$ damping, HERA noise, SNR, Euclid selection, effective galaxy bias |
| `src/figures.py` | **Parts 2–3 figures** — every plot from the two notebooks, headless (`Agg`) and written to `outputs/figures/` |
| `src/dataio.py` | HDF5 loading, halo-catalogue subsampling, and the analysis-product cache |
| `notebooks/plot_fields.ipynb` | **Part 2** — the same field visualisations, interactively (halo catalogue, SFR distributions, lightcone slices, EoR brightness temperature plot) |
| `notebooks/analysis.ipynb` | **Part 3** — the same post-simulation calculations, interactively: 2D cylindrical power spectra, photo-$z$ damping, foreground wedge excision, SNR estimation, Euclid magnitude/SFR cuts, and effective galaxy bias |
| `submit_job.sh` | Launcher for `run_pipeline.py` — activates the conda env, times the run, forwards its arguments, and emails a completion/failure report via `sendmail` |

**Workflow:**

```bash
# One command for everything — simulation runs only if the HDF5 is missing
python run_pipeline.py

# On the HPC cluster, with timing + email notification
bash submit_job.sh --sim force      # force a fresh 21cmFAST run
bash submit_job.sh                  # analyse stored results and re-plot

# Notebooks remain available for interactive exploration
jupyter notebook notebooks/plot_fields.ipynb   # Part 2: field plots
jupyter notebook notebooks/analysis.ipynb      # Part 3: power spectra & SNR
```

**Stage control.** Each stage runs fresh or from stored results, so the
expensive 21cmFAST run happens only when it must:

| Flag | `auto` (default) | `force` | `skip` |
|------|------------------|---------|--------|
| `--sim` | run only if `outputs/lightcone_data.h5` is missing | always re-run `run_simulation.py` | never run; error if there is no stored output |
| `--analysis` | recompute the power spectra only if the cache is missing or older than the simulation | always recompute | load `outputs/analysis_products.h5`; error if absent |

Other options: `--plots {all,none,fields,halos,scaling,power,snr,bias}`,
`--format {png,pdf,svg}`, `--dpi`, `--max-halos N` (uniform strided
subsampling of the halo catalogue when memory is tight — number densities are
rescaled automatically), `--data`, `--products`, `--figdir`, `--summary`.
Run `python run_pipeline.py --help` for the full list.

**Outputs:**

| Path | Contents |
|------|----------|
| `outputs/lightcone_data.h5` | Simulation fields, halo catalogue, and metadata (Part 1) |
| `outputs/analysis_products.h5` | Cached $P_{21}$, $P_\mathrm{gal}$, $P_{21\times\mathrm{gal}}$ and the $k$-grid, plus the `uncertainty_budget` group (damped spectra, wedge mask, $\sigma$ terms, per-mode SNR) |
| `outputs/figures/*.png` | 15 figures: `lightcone_fields`, `lightcone_slice`, `halo_catalogue`, `sfr_relations`, `uv_luminosity_function`, `stellar_mass_muv`, `main_sequence`, `uv_selection_maps`, `power_spectra_2d`, `galaxy_wedge`, `wedge_real_space`, `cross_snr`, `uncertainty_budget`, `photoz_suppression`, `galaxy_bias` |
| `outputs/pipeline_summary.json` | Scalar results: $\langle x_\mathrm{HI}\rangle$, wedge slopes, $\sigma_r$, total SNR, $\langle b_g\rangle$, selection counts |

> **Note:** `submit_job.sh` is a plain shell wrapper — it contains no `#SBATCH`
> directives, so it runs in the foreground on whatever node invokes it. To
> submit it through SLURM with `sbatch`, add the appropriate `#SBATCH`
> directives (`--partition`, `--time`, `--account`, `--cpus-per-task`, …) for
> your cluster at the top of the script first. All arguments passed to
> `submit_job.sh` are forwarded verbatim to `run_pipeline.py`; setting
> `PYTHON_SCRIPT=run_simulation.py` restores its old simulation-only behaviour.

The HDF5 file `outputs/lightcone_data.h5` stores all simulation fields (compressed with gzip) and scalar metadata as attributes, so Parts 2 and 3 are completely independent of the simulation run.

> **Part 1 must be run first.** `outputs/` is listed in `.gitignore`, so a fresh
> clone contains no `lightcone_data.h5` and Parts 2 and 3 will fail at the
> loading cell until `run_simulation.py` has produced it. "Independent of the
> simulation run" means they never import or re-run 21cmFAST — not that the
> HDF5 file ships with the repository. `resources/` (reference PDFs) is
> likewise local-only.

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
│   ├── conftest.py                        # Synthetic lightcone fixtures
│   ├── test_analysis.py
│   ├── test_dataio.py
│   ├── test_figures.py
│   └── test_pipeline.py
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

## Documentation

| Document | Contents |
|----------|----------|
| [`PIPELINE.md`](PIPELINE.md) | HPC pipeline summary, Mermaid flowchart, stage table, and the `lightcone_data.h5` schema |
| [`docs/HPC.md`](docs/HPC.md) | **Parameter-level ground truth for the HPC run** — every configuration value, derived quantity, formula with its evaluated number at $z=7$, file written, disk footprint, and known inconsistency. **§13 is the user-defined parameter checklist**: what you must set on a new cluster, the file and line where each parameter lives, and which edits need `--sim force` to take effect |
| [`docs/uncertainty_budget.md`](docs/uncertainty_budget.md) | **The uncertainty budget, end to end** — every formula of the photo-$z$ / wedge / noise / SNR chain with its evaluated number at $z=7$, the term-by-term audit against the source notebook (bit-identical except one rest-frequency discrepancy), the parameters and their CLI overrides, the HDF5 schema, and what the calculation still does not do |
| [`CHANGELOG.md`](CHANGELOG.md) | Chronological record of all changes, including corrected literature values |
| [`TODO.md`](TODO.md) | Outstanding work, priority-ordered — **including the lightcone power-spectrum corrections the Δz = 1.0 range now requires** |
| [`docs/INSTALL_21cmFASTv4.md`](docs/INSTALL_21cmFASTv4.md) | Step-by-step 21cmFAST v4.1.1 install on CSD3/HPC, plus fixes for quota, conda-plugin, and FFTW linking failures |
| [`docs/Galaxy_bias_formalism.md`](docs/Galaxy_bias_formalism.md) | Methodology for the effective linear galaxy bias from a 21cmFAST halo catalogue under a Euclid $M_\mathrm{UV}$ cut |
| [`docs/halo_catalogue_reference.md`](docs/halo_catalogue_reference.md) | Field-by-field reference for the v4 halo catalogue, verified against py21cmfast v4.x |
| [`docs/Low_SFR_fix.md`](docs/Low_SFR_fix.md) | Diagnosis of the ~7.5 dex SFR offset in the main-sequence plot (the M☉ s⁻¹ unit bug) |
| [`docs/project_update.md`](docs/project_update.md) | Latest run's simulation parameters and numerical results |

## Notebooks

### 1. `21cm_galaxy_cross_uncertainty.ipynb`

Implements the variance estimators (Equations 15–17) for the cross-spectrum and both auto-spectra, along with physically motivated signal models and instrumental noise, to compute per-mode and total signal-to-noise ratios.

**Structure:**

Section numbers below match the notebook's own headers. Note that the numbering
is non-monotonic in places: the Fourier-grid section is labelled **2c** but sits
*after* section 3a in the file, and there are two subsections labelled **3b**.

- **1** Imports
- **2** $T_0(z)$ — brightness-temperature scaling factor across EoR redshifts (Eq. 6)
- **2b** Photo-$z$ damping — Gaussian kernel $W(k_\parallel)$ from photometric redshift uncertainty
- **3** Proxy power spectra
  - **3a** CDM matter power spectrum (BBKS transfer function, $\sigma_8$ normalisation) and noise power spectra
  - **3b** Physical signal spectra — Lidz et al. (2009) reionization bias model
  - **2c** Fourier grid — 2D logarithmic $(k_\perp, k_\parallel)$ grid and helper functions
  - **3.1** Visualise the proxy power spectra
  - **3b** Applying photo-$z$ damping to the observed power spectra
  - **3d** HERA thermal noise — physically motivated $P_N^{21}$ from Eq. 11 with baseline density model
- **4** Variance of the galaxy power spectrum (Eq. 17)
- **5** Variance of the 21 cm power spectrum (Eq. 16)
- **6** Variance of the 21 cm × galaxy cross-power spectrum (Eq. 15)
- **7** Signal-to-noise ratio — overall detectability (7a) and the photo-$z$ uncertainty budget (7b)
- **8** Diagnostic 2-D maps
- **9** 1-D slices along the diagonal $k_\perp = k_\parallel$
- **10** Redshift evolution of the uncertainty at a fixed mode, and of SNR for different photo-$z$ errors (10b)
- **11** Summary and next steps

**Equations implemented:**

| Equation | Quantity | Source |
|----------|----------|--------|
| Eq. 6 | $T_0(z)$ — brightness-temperature scaling | La Plante et al. (2023) |
| Eq. 11 | $P_N^{21}(k_\perp, k_\parallel)$ — HERA thermal noise | La Plante et al. (2023) |
| Eq. 15 | $\sigma^2_{21,\mathrm{gal}}$ — cross-spectrum variance | La Plante et al. (2023) |
| Eq. 16 | $\sigma^2_{21}$ — 21 cm auto-spectrum variance | La Plante et al. (2023) |
| Eq. 17 | $\sigma^2_\mathrm{gal}$ — galaxy auto-spectrum variance | La Plante et al. (2023) |
| Eqs. 2–5 | Signal power spectra (bias model) | Lidz et al. (2009) |

**Fiducial parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| $z_\mathrm{obs}$ | 8.0 | Observing redshift (mid-reionization) |
| $\bar{x}_\mathrm{HI}$ | 0.5 | Mean neutral fraction |
| $b_\mathrm{gal}$ | 5.0 | Galaxy bias |
| $R_b$ | 10 $h^{-1}$ Mpc | Characteristic bubble radius |
| $N_\mathrm{ant}$ | 350 | HERA Phase II antennas |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| $\sigma_z$ | 0.05 | Fiducial photo-$z$ error |

---

### 2. `_archive/21cmfast_HERAxEuclid.ipynb` *(archived)*

> **Archived 2026-08-03.** Half of this notebook's code is duplicated in notebook 3, and its
> remaining unique content is largely re-implemented boilerplate rather than distinct science.
> It is kept in `_archive/` as the only **coeval** reference implementation and for the
> coeval-vs-lightcone comparison; it is no longer part of the active workflow.

End-to-end HERA × Euclid cross-correlation workflow using 21cmFASTv4 with the discrete source (CHMF-SAMPLER) model. Runs a **coeval** simulation at $z = 6.5$, constructs the galaxy density field from the `HaloBox` SFR proxy, and computes 2D cylindrical power spectra with foreground wedge excision and photo-$z$ damping.

**Structure:**

1. Imports and setup
2. Simulation parameters
3. Run 21cmFASTv4 coeval simulation (discrete source model, CHMF-SAMPLER)
4. Galaxy field construction from `HaloBox` SFR density; Kaiser redshift-space distortions applied in Fourier space
5. Visualise 2D slices — 21 cm brightness temperature, galaxy overdensity, and neutral fraction
6. Compute 2D cylindrical power spectra $P_{21}$, $P_\mathrm{gal}$, and $P_\mathrm{cross}$
7. Plot power spectra with foreground wedge and horizon lines
8. Photo-$z$ damping and foreground wedge excision
9. Per-mode and cumulative SNR estimation
10. Summary and references

**Equations implemented:**

| Equation | Quantity | Source |
|----------|----------|--------|
| — | Kaiser RSD: $\delta_\mathrm{gal}^{(s)}(\mathbf{k}) = (1 + \beta\mu^2)\,\delta_\mathrm{gal}(\mathbf{k})$ | Kaiser (1987) |
| Eq. 10 | Foreground wedge slope $m(z)$ | La Plante et al. (2023) |
| — | Photo-$z$ damping kernel $W(k_\parallel) = e^{-k_\parallel^2 \sigma_r^2/2}$ | — |
| Eqs. 15–17 | Cross-spectrum and auto-spectrum variances | La Plante et al. (2023) |

**Fiducial parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `HII_DIM` | 128 | Simulation grid cells per side |
| `BOX_LEN` | 256 Mpc | Comoving box size |
| $z_\mathrm{obs}$ | 6.5 | Observing redshift |
| $M_\mathrm{UV}$ limit | $< -18$ | Euclid galaxy selection |
| $\sigma_z$ | 0.059 | Euclid photo-$z$ uncertainty (as run; this archived coeval notebook predates the absolute-vs-fractional fix) |
| $\bar{n}_\mathrm{gal}$ | $3\times10^{-3}\ h^3\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density |
| $b_\mathrm{gal}$ | 8 | Galaxy bias |
| $N_\mathrm{ant}$ dish diameter | 14 m | HERA antenna diameter |
| $t_\mathrm{obs}$ | 1000 h | Integration time |

**Key outputs:**
- 2D field slices: 21 cm brightness temperature, galaxy overdensity, neutral fraction
- 2D cylindrical power spectra $P_{21}$, $P_\mathrm{gal}$, $P_\mathrm{cross}$ (signed, with wedge overlays)
- Per-mode SNR map in $(k_\perp, k_\parallel)$ space
- Cumulative detection significance outside the foreground wedge

---

### 3. `21cmfast_HERAxEuclid_lightcone.ipynb`

Lightcone counterpart to the archived coeval notebook (2). Uses `RectilinearLightconer` + `run_lightcone` to produce a self-consistent lightcone over a redshift range (default $z = 6.5$–$7.5$), then applies the same galaxy field construction, Kaiser RSD, 2D power spectrum calculation, and SNR estimation as the coeval notebook.

**Key differences from the coeval version:**

| Feature | Coeval | Lightcone |
|---|---|---|
| 21cmFAST function | `run_coeval` | `run_lightcone` |
| Box shape | $(N, N, N)$ | $(N, N, N_z)$ |
| Redshift | single snapshot | continuous range |
| LOS cell size | same as transverse | $L_\mathrm{LOS}/N_z$ |
| Field access | `coeval.brightness_temp` | `lightcone.lightcones['brightness_temp']` |

**Structure:**

Section numbers below match the notebook's own headers.

- **Imports and setup**
- **★ CONFIGURATION** — all user-adjustable parameters in one cell
- **1** Derived quantities (LOS geometry, node redshifts)
- **2** Run 21cmFASTv4 lightcone simulation
- **3** Construct galaxy density field — from lightcone `halo_sfr` (3a) with a
  synthetic Poisson-sampled fallback for zero-field cases (3b); galaxy bias
  estimated via HMF integration over the Euclid UV magnitude range
- **4** Kaiser RSD applied in Fourier space
- **5** Visualise lightcone — transverse (x–y) slice + LOS (x–z) slice with
  **LOS on the x-axis and transverse on the y-axis** (standard convention);
  secondary redshift axis on top via `twiny()`
- **5b** Brightness temperature evolution plot — wide-format (16×3.5") lightcone
  slice styled after Mesinger & Furlanetto (2007), with a custom EoR colourmap
  (dark = ionised, warm/bright = neutral), dual x-axes (comoving distance +
  redshift), and observed frequency range in the title
- **6** Compute 2D cylindrical power spectra (non-cubic box)
- **7** Plot power spectra with foreground wedge overlays
- **7b** Photo-$z$ suppression — the damping kernel
  $W(k_\parallel) = e^{-k_\parallel^2\sigma_r^2/2}$ swept over
  $\sigma_z \in \{0,\,0.02,\,0.05,\,0.10,\,0.30,\,0.45\}$, from the
  spectroscopic limit ($\sigma_r = 0$, $W \equiv 1$) to the adopted Euclid
  Wide value, with $1/\sigma_r$ marked. Both $\sigma_r$ and $W$ come from
  `src.analysis.radial_smearing_length` / `photoz_damping_kernel` — the same
  two functions `compute_uncertainty_budget()` applies in §8
- **7c** Galaxy power spectrum against the foreground wedge — $P_\mathrm{gal}$
  alone, with the wedge region **filled and hatched** rather than only
  outlined, so excluded modes are distinguishable from the accessible window at
  a glance. Reuses §7's `horizon_slope`, `fov_wedge_slope_value` and wedge lines
- **7d** The wedge in **real** space — the 3D galaxy overdensity FFT'd, every
  mode with $k_\parallel \le m_{\rm horizon} k_\perp$ zeroed, and
  inverse-transformed back; the same LOS slice shown before and after on a
  shared colour scale, titled with the percentage of modes removed. The mask is
  `src.analysis.foreground_wedge_mask` reshaped onto §4's `KX/KY/KZ` grids.
  **97.4 % of the 3D modes fall inside the bare horizon wedge** at $z = 7$ on
  the fiducial $(128, 128, 175)$ grid
- **8** The uncertainty budget — photo-$z$ damping, wedge excision, noise and
  SNR in a single `compute_uncertainty_budget()` call
- **9** Per-mode SNR map and total detection significance
- **10** Summary — coeval vs. lightcone comparison table and next-step recommendations

**Equations implemented:**

| Equation | Quantity | Source |
|----------|----------|--------|
| $H(z) = H_0\sqrt{\Omega_m(1+z)^3+(1-\Omega_m)}$ | Hubble parameter for LOS geometry | — |
| — | Kaiser RSD: $\delta_\mathrm{gal}^{(s)}(\mathbf{k}) = (1+\beta\mu^2)\,\delta_\mathrm{gal}(\mathbf{k})$ | Kaiser (1987) |
| Eq. 10 | Foreground horizon wedge slope $m(z)$ | La Plante et al. (2023) |
| — | Photo-$z$ damping kernel $W(k_\parallel) = e^{-k_\parallel^2\sigma_r^2/2}$ | — |
| Eqs. 15–17 | Cross-spectrum and auto-spectrum variances | La Plante et al. (2023) |

**Fiducial parameters:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `HII_DIM` | **256** | Simulation grid cells per side — *derived*, not hardcoded (see below) |
| `BOX_LEN` | **486.33 Mpc** | Transverse box size — *derived* from the 10 deg² Euclid Deep Field Fornax footprint at $z=7$ |
| `DIM` | `3 × HII_DIM` = **768** | High-res grid for initial conditions |
| $z_\mathrm{min}$, $z_\mathrm{max}$ | **6.55, 7.45** | Lightcone redshift range — *derived* from $\Delta z = 2\sigma_z = 0.90$ (±1σ) |
| $z_\mathrm{obs}$ | 7.0 | Reference redshift (midpoint of lightcone) |
| $M_\mathrm{UV}$ limit | $< -18$ | Euclid galaxy selection |
| $\sigma_z$ | 0.45 | Euclid photo-$z$ uncertainty, **absolute — not $\sigma_z/(1+z)$**. Euclid's $\sigma_z/(1+z) < 0.05$ requirement is $\sigma_z \approx 0.45$ at $z = 7$; the earlier 0.059 was the fractional number used as though absolute |
| $\bar{n}_\mathrm{gal}$ | $3\times10^{-3}\ h^3\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density |
| $b_\mathrm{gal}$ | **computed in-line**, $\approx 5.39$ | No fixed fallback — `galaxy_bias = 8` is commented out in the configuration cell. The value is calculated in **§3**, the analytic-galaxy-bias cell (`galaxy_bias = galaxy_bias_hmf`): a Simpson integral of the Sheth-Tormen halo bias weighted by the HMF over the Euclid-selected mass bins, $b_g = \int b_h\,\frac{dn}{d\log M}\,d\log M \big/ \int \frac{dn}{d\log M}\,d\log M$. Consumed downstream by §4's $\beta = f/b$. Requires `hmf` |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| Bandwidth | 8 MHz | HERA per-band bandwidth |

---

### 4. HPC pipeline — `run_simulation.py` + `notebooks/plot_fields.ipynb` + `notebooks/analysis.ipynb`

A refactored version of notebook 3 split into three independent parts for cluster use. See the [HPC lightcone pipeline](#hpc-lightcone-pipeline-recommended-for-cluster-use) section above.

**Fiducial parameters** (configuration block at the top of `run_simulation.py`) — these are *not* identical to notebook 3's:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `SURVEY_AREA_DEG2` | 10 deg² | Euclid Deep Field Fornax footprint (RA 03:31:43.6, Dec −28:05:18.6) |
| `PHOTOZ_N_SIGMA` | 1 | Box spans ±1σ_z of photo-z scatter → $\Delta z = 0.90$ |
| `HII_DIM` | **256** | Simulation grid cells per side — *derived*, not hardcoded |
| `BOX_LEN` | **486.33 Mpc** | Transverse box size — *derived* from the footprint at $z=7$ |
| `DIM` | `3 × HII_DIM` = **768** | High-res grid for initial conditions |
| $z_\mathrm{min}$, $z_\mathrm{max}$ | **6.995, 7.005** | Lightcone redshift range — smoke-test slab, see the note below |
| $z_\mathrm{obs}$ | 7.0 | Reference redshift (midpoint) |
| `minimum_los_slices` | 100 | Floor on $N_z$, overriding the cell-size-matched value |
| $M_\mathrm{UV}$ limit | $< -18$ | Euclid galaxy selection (bright end $-22$) |
| $\sigma_z$ | **0.45** | Euclid photo-$z$ uncertainty — **absolute**, not $\sigma_z/(1+z)$. Equals $\sigma_z/(1+z) = 0.056$ at $z=7$, consistent with the Euclid requirement $< 0.05$. Was 0.059 (the fractional value used as if absolute), which understated $\sigma_r$ by $\sim 7.6\times$ |
| $\bar{n}_\mathrm{gal}$ | $3\times10^{-3}\ h^3\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density |
| $b_\mathrm{gal}$ | 8 | Fallback only; overwritten by the halo-catalogue estimate ($\approx 4.7$) |
| $t_\star$ | 0.5 | SFR timescale as a fraction of $t_H(z)$ — 570 Myr at $z=7$ |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| Bandwidth | 8 MHz | HERA per-band bandwidth |
| Wedge buffer | 0.0677 Mpc$^{-1}$ | Safety margin beyond the horizon line — $0.1\ h\ \mathrm{Mpc}^{-1}$ at $h = 0.6766$ (Pober et al. 2014 "moderate" foreground model; the 21cmSense `horizon_buffer` default) |
| PS binning | 20 × 20 | $(k_\perp, k_\parallel)$ bins |
| Source model | `from_template(["simple"])`, `random_seed=42` | Grid-based halos, no Ts fluctuations |

> **Thin-slab test configuration — deliberate, and currently required.** The
> committed $z$ range spans only $\Delta z = 0.01$, i.e.
> $L_\mathrm{LOS} = 3.5$ Mpc at $z = 7$. The cell-size-matched slice count
> would be $N_z = 2$, but `minimum_los_slices` raises it to 100, giving a
> 0.035 Mpc LOS cell against a 2 Mpc transverse cell — a $\sim 57\times$
> oversampling along the line of sight. It produces a quasi-coeval slab with
> negligible redshift evolution.
>
> This is **not** merely a leftover. The power-spectrum estimator in
> `src/analysis.py` assumes statistical homogeneity along the LOS, which only
> holds for a quasi-coeval box, so configuration and formalism currently
> match. Widening to a true lightcone ($z = 6.5$–$7.5$,
> $L_\mathrm{LOS} = 350.8$ Mpc, $N_z = 175$) requires the estimator work in
> [`TODO.md`](TODO.md) §P0 first — otherwise the FFT is applied to unevenly
> sampled ($20.4\%$ cell-size spread), redshift-evolving data. Treat the
> resulting SNR as a smoke-test number, not a forecast.

> **The stored `outputs/lightcone_data.h5` is still stale.** The galaxy-bias
> calculation has been corrected (see
> [`docs/project_update.md`](docs/project_update.md) §4), which changes
> `galaxy_bias` from 33.39 to 4.744 and $\beta_\mathrm{rsd}$ from 0.030 to
> 0.210 — and therefore the `galaxy_overdensity` field, which carries the
> Kaiser boost. Regenerate with `bash submit_job.sh --sim force`; the default
> `--sim auto` will **not** re-run while the old file is present. With the
> redshift range back at $\Delta z = 0.01$ this is a cheap re-run.

**Galaxy bias.** Two estimators are computed at the end of Part 1:

| Estimator | $b_g$ | Basis |
|-----------|-------|-------|
| Halo catalogue (**adopted**) | 4.74 | Per-halo $M_\mathrm{UV}$ from the catalogue's own SFR, then the mean Sheth-Tormen bias over the Euclid-selected halos. Inherits 21cmFAST's log-normal SFR scatter. |
| Analytic HMF integral (cross-check) | 5.39 | Sheth-Tormen integrated over the mass function weighted by the *mean, scatter-free* scaling relation. |

Both are written to the HDF5 (`galaxy_bias`, `galaxy_bias_hmf_analytic`,
`galaxy_bias_method`). The previously recorded $b_g = 33.4$ was an artefact of
passing `hmf`'s already-squared `MassFunction.nu` into a helper that squared it
again; that path is gone.

#### `notebooks/plot_fields.ipynb` — figures and literature references

Produces seven groups of plots from `outputs/lightcone_data.h5`. The literature equations used in each are documented below.

**Figure group 1–2: Halo catalogue and SFR distributions**

No literature overlays. These are diagnostic plots (projected halo positions, halo mass histogram, halos overlaid on a $\delta T_b$ slice, SFR histogram, SFR vs $M_\mathrm{halo}$, SFR vs $M_\star$) showing the raw 21cmFAST catalogue.

---

**Figure group 3–4: Lightcone field slices and wide-format EoR lightcone**

No literature overlays. These are visualisation panels (transverse slice, LOS slice with dual distance/redshift axes, neutral fraction). The wide-format panel uses a custom EoR colourmap styled after Mesinger & Furlanetto (2007).

---

**Figure group 5: UV Luminosity Function**

Simulation data points are converted from SFR to $M_\mathrm{UV}$ and binned into a UVLF. Two Schechter-function curves from the literature are overlaid.

*SFR → $M_\mathrm{UV}$ conversion chain* (implemented in `src/conversions.sfr_to_Muv()`):

| Step | Equation | Reference |
|------|----------|-----------|
| SFR → $L_\mathrm{UV}$ | $L_\mathrm{UV} = \mathrm{SFR} / \kappa_\mathrm{UV}$, $\kappa_\mathrm{UV} = 1.15\times10^{-28}\ M_\odot\ \mathrm{yr}^{-1}\ /\ (\mathrm{erg\ s}^{-1}\ \mathrm{Hz}^{-1})$ | Madau & Dickinson (2014), Chabrier IMF, ~1500 Å |
| $L_\mathrm{UV}$ → $M_\mathrm{UV}$ | $M_\mathrm{AB} = 51.60 - 2.5\log_{10}(L_\nu)$ | Oke & Gunn (1983) |

*Schechter function form* (implemented inline as `schechter_muv()`):

$$\Phi(M) = \frac{\ln 10}{2.5}\,\phi_\star\,10^{0.4(\alpha+1)(M_\star - M)}\,\exp\!\left[-10^{0.4(M_\star - M)}\right]$$

Reference: Schechter (1976, ApJ 203, 297)

*Literature Schechter parameters at $z \sim 7$* (hardcoded in the `literature` list):

| Reference | $\phi_\star$ [Mpc$^{-3}$] | $M_\star$ | $\alpha$ | Source table |
|-----------|--------------------------|-----------|----------|-------------|
| Bouwens et al. (2021, AJ 162, 47) | $1.9\times10^{-4}$ | $-21.15$ | $-2.06$ | Table 5, single-Schechter, $z=6.8$ |
| Finkelstein et al. (2015, ApJ 810, 71) | $1.57\times10^{-4}$ | $-21.03$ | $-2.03$ | Table 4, $z=7$ |

---

**Figure group 6: Stellar Mass – UV Magnitude Relation**

The simulation median and 16–84th percentile scatter band are plotted. Two literature relations are defined in the code; only the first is currently rendered:

| Reference | Relation (as coded) | Status | Note |
|-----------|---------------------|--------|------|
| Song et al. (2016, ApJ 825, 5) | $\log_{10}(M_\star/M_\odot) = 8.86 - 0.5\,(M_\mathrm{UV} + 20)$ | **Plotted** | Anchored at $M_\mathrm{UV}=-21 \to \log_{10} M_\star = 9.36$; slope from their Figure 5 |
| González et al. (2010, ApJ 713, 115) | $\log_{10}(M_\star/M_\odot) = 9.06 - 0.5\,(M_\mathrm{UV} + 20)$ | Defined but **commented out** | $\sim 0.2$ dex higher normalisation from constant-SFH SED assumption |

The $M_\mathrm{UV}$ axis uses the same `sfr_to_Muv()` chain as Figure group 5.

---

**Figure group 7: Star-forming Main Sequence**

Simulation median and scatter band (all halos, no magnitude cut) are plotted alongside three curves:

| Curve | Equation (as coded) | Reference |
|-------|---------------------|-----------|
| Speagle+14 | $\log_{10}\mathrm{SFR} = (0.84 - 0.026\,t)\,\log_{10}M_\star - (6.51 - 0.11\,t)$, $t$ = cosmic age [Gyr] via `Planck18.age(z_obs)` | Speagle et al. (2014, ApJS 214, 15), Eq. 28 |
| Schreiber+15 | $\log_{10}\mathrm{SFR} = \log_{10}M_\star - 8$ (i.e., $\mathrm{sSFR} = 10\ \mathrm{Gyr}^{-1}$) | Schreiber et al. (2015, A&A 575, A74), high-$z$ limit |
| 21cmFAST model | $\log_{10}\mathrm{SFR} = \log_{10}M_\star - \log_{10}t_\mathrm{sf}$, $t_\mathrm{sf} = t_\star\,t_H(z)$, $t_\star = 0.5$ | Park et al. (2019), Eq. 3; `scaling_relations.c` |

The 21cmFAST model line has no free parameters — it is a pure prediction from the simulation's SFR prescription. At $z = 7$, $t_H \approx 1.14$ Gyr gives $t_\mathrm{sf} \approx 570$ Myr and $\mathrm{sSFR} \approx 1.75\ \mathrm{Gyr}^{-1}$, which is $\sim 0.76$ dex below Schreiber+15. This offset is a physical model difference (the 21cmFAST `simple` template assumes a longer star-formation timescale), not a numerical error.

---

## 21cmFASTv4 `HaloBox` API Notes

In 21cmFAST v4.1+, `coeval.halobox` is a `HaloBox` object whose arrays are accessed via `.get('<field_name>')`. In lightcone runs, the same fields are stored in `lightcone.lightcones['<field_name>']`. Available fields include:

| Field | Description |
|-------|-------------|
| `halo_sfr` | Total SFR per cell, summed over all halos [internal units — absolute value cancels in $\delta_\mathrm{gal}$] |
| `n_ion` | Number of ionising photons per cell |

Individual halo positions, stellar masses, and SFRs are not exposed via the lightcone `lightcones` dict. The per-halo catalogue can be retrieved after the simulation using `determine_halo_catalog` + `perturb_halo_catalog` (see `run_simulation.py`, Section 3a).

> **SFR unit warning (py21cmfast v4):** `perturbed_halos.sfr.value` returns the SFR in **M☉ s⁻¹** (the internal C unit from `scaling_relations.c`), not M☉ yr⁻¹ as stated in the documentation. Multiply by `365.25 × 24 × 3600 = 3.15576 × 10⁷` to convert to M☉ yr⁻¹. `halo_masses` and `stellar_masses` are unaffected. See `docs/Low_SFR_fix.md` for full diagnosis and verification.

The galaxy overdensity is constructed as $\delta_\mathrm{gal} = \mathrm{SFR}/\langle\mathrm{SFR}\rangle - 1$, which correctly traces the galaxy distribution and produces the expected large-scale anti-correlation with the 21 cm field (negative cross-spectrum on large scales).

### Lightcone API

```python
import numpy as np
import py21cmfast as p21c
from astropy.cosmology import Planck18
import astropy.units as u

# Compute LOS geometry
N_z          = int(round(L_los / cell_size))
lc_redshifts = np.linspace(z_min, z_max, N_z)
node_redshifts = np.linspace(z_max, z_min, n_nodes)   # high-z → low-z

inputs = p21c.InputParameters.from_template(["simple"], random_seed=42)
inputs = inputs.clone(node_redshifts=node_redshifts, ...)

lightconer = p21c.RectilinearLightconer(
    lc_redshifts=lc_redshifts,
    quantities=("brightness_temp", "density", "neutral_fraction", "halo_sfr"),
)
lightcone = p21c.run_lightcone(lightconer=lightconer, inputs=inputs)

# Access fields — shape (HII_DIM, HII_DIM, N_z)
brightness_temp = lightcone.lightcones["brightness_temp"]
neutral_frac    = lightcone.lightcones["neutral_fraction"]   # note: not "xH_box"
halo_sfr        = lightcone.lightcones["halo_sfr"]

L_los = lightcone.lightcone_dimensions[2]   # actual LOS comoving size [Mpc]
```

### 21cmFAST source model templates

`InputParameters.from_template()` accepts a template name or list of names (later entries override earlier ones). The template sets all physics toggle flags; size templates can be stacked on top of a physics template.

**Available physics templates:**

| Template | Key physics enabled | Typical use |
|---|---|---|
| `simple` | Grid-based halos only | Fast tests, reionisation morphology without Ts |
| `latest` | Grid-based halos + Ts fluctuations + inhomogeneous recombinations | Production runs |
| `latest-discrete` | Same as `latest` + discrete (sampled) halo field | Highest fidelity |
| `minihalos` | `latest` + molecularly cooled / PopIII minihalos | X-ray / Lyman-Werner background studies |
| `Park19` | Exact Park et al. (2019) fiducial | Comparison runs |
| `Qin20` | Exact Qin et al. (2020) fiducial | Comparison runs |

`simple` specifically disables: `USE_TS_FLUCT`, `INHOMO_RECO`, `CELL_RECOMB`, `USE_MINI_HALOS`; sets `HII_FILTER = 'sharp-k'` and `SOURCE_MODEL = 'E-INTEGRAL'`.

**Available size templates** (stack on top of a physics template): `size-tiny` (32 cells, 48 Mpc), `size-small` (64 cells, 92 Mpc), `size-medium` (256 cells, 384 Mpc), `size-gpc` (~1 Gpc).

```python
# Current project configuration — fast but no Ts fluctuations:
inputs = p21c.InputParameters.from_template(["simple"], random_seed=42)

# Full physics (Ts fluctuations + recombinations), same grid:
inputs = p21c.InputParameters.from_template(["latest"], random_seed=42)

# Full physics + override grid size:
inputs = p21c.InputParameters.from_template(["latest", "size-medium"], random_seed=42)

# List all available templates:
from py21cmfast._templates import list_templates
print(list_templates())
```

---

## Source Modules

### `src/conversions.py`

Cosmological conversion utilities for high-redshift galaxy surveys. Import
individual functions as needed:

```python
from src.conversions import (
    Muv_to_Luv, Luv_to_Muv,
    Luv_to_sfr, sfr_to_Luv, sfr_to_Muv,
    sheth_tormen_bias,
    mean_matter_density, cell_mass,
    survey_area_from_volume,
    area_deg2_to_steradians,
    volume_from_area,
)
```

**Magnitude–luminosity conversions**

| Function | Description |
|----------|-------------|
| `Muv_to_Luv(Muv)` | Absolute UV AB magnitude → monochromatic luminosity [erg s⁻¹ Hz⁻¹] |
| `Luv_to_Muv(Luv)` | Monochromatic luminosity [erg s⁻¹ Hz⁻¹] → absolute UV AB magnitude |

**UV luminosity – SFR conversions** (Madau & Dickinson 2014, Chabrier IMF, $\kappa_\mathrm{UV} = 1.15\times10^{-28}\ M_\odot\ \mathrm{yr}^{-1}\ /\ (\mathrm{erg\ s}^{-1}\ \mathrm{Hz}^{-1})$)

| Function | Description |
|----------|-------------|
| `Luv_to_sfr(Luv, kappa_uv=1.15e-28)` | UV luminosity [erg s⁻¹ Hz⁻¹] → SFR [M☉ yr⁻¹] |
| `sfr_to_Luv(sfr, kappa_uv=1.15e-28)` | SFR [M☉ yr⁻¹] → UV luminosity [erg s⁻¹ Hz⁻¹] |
| `sfr_to_Muv(sfr, kappa_uv=1.15e-28)` | SFR [M☉ yr⁻¹] → absolute UV AB magnitude (chains `sfr_to_Luv` → `Luv_to_Muv`) |

**Halo bias**

| Function | Description |
|----------|-------------|
| `sheth_tormen_bias(nu_sq, delta_c=1.686, a=0.707, p=0.3)` | Sheth-Tormen (1999) Eulerian halo bias. `nu_sq` is $(δ_c/σ)^2$ as returned by `hmf.MassFunction.nu` — **not** $δ_c/σ$ |

> **`hmf` convention note:** `MassFunction.nu` in `hmf` ≥ 3.x stores the
> *squared* peak height $(δ_c/σ)^2$, consistent with the original Sheth &
> Tormen (1999) notation. Pass `mf.nu` directly to `sheth_tormen_bias` — do
> not square it again.

**Mass resolution**

| Function | Description |
|----------|-------------|
| `mean_matter_density(omega_m, hubble_constant)` | Comoving mean matter density $\bar\rho_m = \Omega_m \rho_{\mathrm{crit},0}$ [M☉ Mpc⁻³]. Defaults to Planck18 |
| `cell_mass(cell_size_mpc, omega_m, hubble_constant)` | Mean matter mass enclosed by one cubic comoving cell, $M_\mathrm{cell} = \bar\rho_m L_\mathrm{cell}^3$ [M☉] — the grid mass resolution |

For the production grid (`BOX_LEN = 486.33` Mpc, `HII_DIM = 256`, `DIM = 768`,
$\Omega_m = 0.315$, $H_0 = 67.36$) — sized from the survey footprint by
`survey_area_to_box_size`, which derives `HII_DIM` from a 2.0 Mpc target cell
so growing the box does not coarsen the mass resolution:

| Grid | Cell size | Mass resolution |
|------|-----------|-----------------|
| `DIM` (initial conditions / density field) | 0.633 Mpc | $1.01\times10^{10}\ M_\odot$ per cell |
| `HII_DIM` (ionisation, 21 cm brightness) | 1.90 Mpc | $2.72\times10^{11}\ M_\odot$ per cell |
| Halo catalogue (`SAMPLER_MIN_MASS`) | — | $1\times10^{8}\ M_\odot$ (smallest sampled halo) |

> The halo catalogue is *not* limited by the grid mass resolution: 21cmFAST's
> stochastic halo sampler populates each cell down to `SAMPLER_MIN_MASS`, so
> the galaxy field resolves halos ~100× lighter than a single `DIM` cell. The
> grid masses set the resolution of the *density* field, not the halo field.
> All three values are printed by the parameter-summary cell of
> `notebooks/plot_fields.ipynb` and stored as HDF5 attributes
> (`M_cell_hires`, `M_cell_lores`, `sampler_min_mass`).

**Survey geometry conversions**

| Function | Description |
|----------|-------------|
| `survey_area_from_volume(volume_mpc3, z_min, z_max, cosmo=None)` | Comoving volume [Mpc³] → survey area [deg²] |
| `area_deg2_to_steradians(area_deg2)` | Survey area [deg²] → [sr] |
| `volume_from_area(area_deg2, z_min, z_max, cosmo=None, n_z=1000)` | Survey area [deg²] → comoving volume [Mpc³] |
| `survey_area_to_box_size(area_deg2, z_central, delta_z, cosmo=None, target_cell_size_mpc=2.0, hii_dim=None, snap_hii_dim_to_power_of_two=True)` | Survey footprint → `SimulationBox` carrying 21cmFAST's `BOX_LEN` / `HII_DIM` / `DIM` |

All functions accept scalar or array inputs. Volume–area conversions use
Simpson integration of the differential comoving volume
$\mathrm{d}V/\mathrm{d}z\,\mathrm{d}\Omega$ (Hogg 1999) and default to the
Planck18 cosmology; pass a custom `astropy.cosmology` object via `cosmo` to
override.

**Sizing the simulation box from the survey footprint**

`survey_area_to_box_size` replaces the hardcoded `BOX_LEN = 256` Mpc /
`HII_DIM = 128` grid with a box traceable to the survey being forecast. It
returns a `SimulationBox` dataclass whose `.simulation_options` property is the
`{"HII_DIM", "BOX_LEN", "DIM"}` mapping 21cmFAST's
`InputParameters.clone(simulation_options=...)` expects.

```python
from src.conversions import survey_area_to_box_size

box = survey_area_to_box_size(area_deg2=10.0, z_central=7.0, delta_z=0.90)
inputs = inputs.clone(simulation_options=box.simulation_options)
```

| Step | Formula | Result (Fornax) |
|------|---------|-----------------|
| Transverse extent | $L_\perp = \sqrt{\Omega}\;D_M(z_c)$, small-angle, square footprint | **486.33 Mpc** → `BOX_LEN` |
| Line-of-sight depth | $L_\parallel = D_C(z_c + \Delta z/2) - D_C(z_c - \Delta z/2)$ | **315.60 Mpc** |
| Grid | $N = \lceil L_\perp / 2.0\,\mathrm{Mpc} \rceil$, snapped to a power of two | 244 → **256** = `HII_DIM`, `DIM` = 768 |

The LOS depth is deliberately computed by **differencing comoving distances**
rather than $\mathrm{d}D_C/\mathrm{d}z \times \Delta z$, matching how
`run_simulation.py` §1 and the notebook already compute `L_los`.

$L_\parallel$ is **not** a 21cmFAST argument — boxes are cubic, and a
lightcone's LOS extent comes from the redshift range handed to
`RectilinearLightconer`. Use the returned `z_min` / `z_max` for that. The
returned `n_los_tiles` = $L_\parallel / L_\perp$ flags when the coeval box
would have to be tiled along the LOS (0.65 for Fornax, so no tiling).

**Assumptions, for the thesis writeup** (all recorded in the docstring):

- *Survey geometry* — Euclid Deep Field Fornax, 10 deg², centred RA 03:31:43.6,
  Dec −28:05:18.6. Treated as square, so $L_\perp$ is an equivalent-square side.
- *Redshift depth* — set by $\sigma_z = 0.45$ **absolute** at $z = 7$, from
  Euclid's fractional requirement $\sigma_z/(1+z) < 0.05$. The multiple of
  $\sigma_z$ is a **deliberate choice, passed explicitly**, not a default: the
  forecast adopts ±1σ ($\Delta z = 0.90$); ±2σ would give $\Delta z = 1.80$,
  $L_\parallel = 634.9$ Mpc.
- *Cosmology* — astropy `Planck18`, following the `cosmo=None` convention of the
  other survey-geometry functions here and the Planck18 distances
  `run_simulation.py` already uses for lightcone endpoints. Note `src/analysis.py`
  takes literal $H_0 = 67.36$, $\Omega_m = 0.315$ instead; the two differ by
  ~0.4 % in $H_0$ and are not interchangeable at that precision.
- *Resolution* — `HII_DIM` is derived from `target_cell_size_mpc = 2.0`, the
  resolution of the old 256 Mpc / 128³ grid, so covering the larger footprint
  does not silently coarsen $M_\mathrm{cell}$. Passing `hii_dim` explicitly
  overrides this and *does* change the mass resolution.

**Dependencies:** `numpy`, `astropy`, `scipy`

### `src/analysis.py` — galaxy overdensity weighting

$\delta_\mathrm{gal}$ can be built from the Euclid-selected halo catalogue with
either of two per-halo weights. Both deposit the same halos onto the same grid
and normalise the same way, so they are drop-in replacements for one another:

| Mode | Formula | Traces |
|------|---------|--------|
| `"number"` (default) | $\delta_\mathrm{gal} = N / \langle N \rangle - 1$ | the *abundance* of detectable galaxies |
| `"luminosity"` | $\delta_{\mathrm{gal},L} = \sum L_\mathrm{UV} / \langle \sum L_\mathrm{UV} \rangle - 1$ | the *UV emissivity*, up-weighting bright halos |

| Function | Description |
|----------|-------------|
| `deposit_halo_field(coords, box_len, n_perp, n_los=None, los_extent=None, weights=None)` | Deposit halos onto an `(n_perp, n_perp, n_los)` grid, summing `weights` (unit weights if `None`) |
| `galaxy_overdensity_from_catalogue(coords, sfr, halo_masses, box_len, n_perp, ..., weighting="number")` | Apply the Euclid $M_\mathrm{UV}$ window, deposit, and normalise → `(delta_gal, EuclidSelection)` |
| `GALAXY_WEIGHTING_MODES` | `("number", "luminosity")` |

$L_\mathrm{UV}$ comes from `conversions.sfr_to_Luv()`
($L_\mathrm{UV} = \mathrm{SFR}/\kappa_\mathrm{UV}$, Madau & Dickinson 2014).
Because that is a constant rescaling, it divides out of the ratio — the
luminosity-weighted field is identical to an SFR-weighted one, and differs
from the number-weighted field only through the per-halo spread in SFR.

> **Two grid shapes, deliberately.** `deposit_halo_field` is used at both
> shapes and they must not be conflated:
>
> | Field | Shape | Built from |
> |---|---|---|
> | UV luminosity / selected-galaxy maps (notebook §3) | `(HII_DIM, HII_DIM, HII_DIM)` — **cubic** | the *coeval* perturbed halo catalogue, whose coords span `[0, BOX_LEN)` on all three axes |
> | $\delta_\mathrm{gal}$ for the power spectra (notebook §4) | `(HII_DIM, HII_DIM, N_z)` | the *lightcone* `halo_sfr` field |
>
> The cubic maps are diagnostics only. Feeding one to the power-spectrum
> estimator would be both a shape and a geometry mismatch, since the LOS axis
> of the lightcone spans $L_\mathrm{los}$ with $N_z$ cells, not `BOX_LEN` with
> `HII_DIM`. Halos outside `[0, BOX_LEN]` are dropped by the underlying
> `histogramdd` (and counted in a printed warning), not wrapped periodically.

**Selecting a mode in `run_simulation.py`** — set the `GALAXY_WEIGHTING`
constant in the Euclid survey-parameter block:

| Value | Source of `delta_gal` |
|-------|-----------------------|
| `"lightcone_sfr"` (default) | the lightcone `halo_sfr` field, `sfr_field / mean_sfr - 1` — the original behaviour, and the only mode that evolves along the line of sight |
| `"number"` | Euclid-selected catalogue, unit weights |
| `"luminosity"` | Euclid-selected catalogue, $L_\mathrm{UV}$ weights |

The chosen array flows through the identical downstream path (Kaiser RSD →
HDF5 `galaxy_overdensity` → `compute_all_power_spectra()`), and the mode is
recorded as the root attribute `galaxy_weighting`.

> **Line-of-sight caveat:** the halo catalogue is *coeval* at `z_obs` and spans
> `BOX_LEN` along the LOS, whereas the lightcone spans `L_los`. The catalogue
> modes are binned into an `(HII_DIM, HII_DIM, N_z)` grid so the shape matches
> downstream, but they carry no redshift evolution along the LOS and their LOS
> cell size is `BOX_LEN / N_z`.

### `src/FOV_to_cMpc.py`

Standalone command-line utility that converts a survey field of view and
redshift range into comoving survey geometry — solid angle, comoving volume,
transverse comoving area, and the equivalent cubic box side length. Useful for
checking that a simulation `BOX_LEN` is representative of the intended Euclid
survey footprint.

```bash
conda run -n 21cmfast python src/FOV_to_cMpc.py --area-deg2 0.53 --z-min 6.0 --z-max 7.0
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--area-deg2` | yes | Survey field of view / area [deg²] |
| `--z-min` | yes | Lower redshift bound |
| `--z-max` | yes | Upper redshift bound |
| `--n-z` | no (default 1000) | Redshift samples for Simpson integration |

It also exposes two importable functions:

| Function | Description |
|----------|-------------|
| `survey_volume_from_area(area_deg2, z_min, z_max, cosmology=Planck18, n_z=1000)` | → `(volume [Mpc³], solid angle [sr], volume per steradian)` |
| `transverse_comoving_size_from_area(area_deg2, z, cosmology=Planck18)` | → `(side length [Mpc], transverse area [Mpc²])` at a single redshift |

**Dependencies:** `numpy`, `astropy`, `scipy`

> **Overlap note:** `survey_volume_from_area` here and `volume_from_area` in
> `conversions.py` compute the same quantity by the same method;
> `FOV_to_cMpc.py` additionally returns the intermediate solid angle and
> per-steradian volume, and wraps everything in a CLI.

---

## Installation

All commands must run inside the `21cmfast` conda environment.

```bash
# Create the environment (installs fftw/gsl, then pip-installs 21cmFAST)
conda env create -f env.yml
conda activate 21cmfast

# Verify the 21cmFAST C extensions built correctly
python -c "import py21cmfast; print(py21cmfast.__version__)"   # → 4.1.1
```

`21cmFAST` compiles C extensions at install time and links against FFTW and
GSL, which is why those must come from conda-forge *before* the pip step —
`env.yml` orders this correctly.

**On HPC systems**, follow [`docs/INSTALL_21cmFASTv4.md`](docs/INSTALL_21cmFASTv4.md)
instead. It covers the CSD3-specific problems that `conda env create` will not
solve on its own: home-directory quota exhaustion during the pip build
(redirect `PIP_CACHE_DIR`/`XDG_CACHE_HOME` to scratch), `conda-libmamba-solver`
entry-point failures (`CONDA_NO_PLUGINS=true`), and FFTW linking errors.

To reproduce the exact environment this project was last run in, use the pinned
freeze instead of the version ranges in `env.yml`:

```bash
pip install -r requirements.txt
```

## Requirements

| Package | Used in |
|---------|---------|
| `numpy` | All notebooks, `src/conversions.py` |
| `matplotlib` | All notebooks |
| `scipy` | `21cmfast_HERAxEuclid_lightcone.ipynb`, `notebooks/analysis.ipynb`, `src/conversions.py` |
| `py21cmfast >= 4.1.1` | `21cmfast_HERAxEuclid_lightcone.ipynb`, `run_simulation.py` |
| `astropy` | `21cmfast_HERAxEuclid_lightcone.ipynb`, `src/conversions.py` |
| `hmf` | `run_simulation.py`, `notebooks/analysis.ipynb` |
| `h5py` | `run_simulation.py`, `notebooks/plot_fields.ipynb`, `notebooks/analysis.ipynb` |
| `ipympl` | All notebooks — backs the `%matplotlib widget` interactive inline backend |

The analytical notebook (`21cm_galaxy_cross_uncertainty.ipynb`) requires only `numpy` and `matplotlib`; all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).

## Usage

```bash
# Activate the 21cmfast conda environment first
conda activate 21cmfast

# Analytical framework (no external dependencies beyond numpy/matplotlib)
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb

# 21cmFAST lightcone simulation pipeline — monolithic version
jupyter notebook 21cmfast_HERAxEuclid_lightcone.ipynb

# HPC-optimised lightcone pipeline — recommended for cluster use
python run_pipeline.py                                 # everything, one command
python run_pipeline.py --sim force                     # re-run the simulation first
python run_pipeline.py --plots power snr               # only the k-space figures
python run_pipeline.py --max-halos 5000000             # cap catalogue memory
bash submit_job.sh --sim force                         # same, with timing + email

# Interactive exploration of the same results
jupyter notebook notebooks/plot_fields.ipynb           # Part 2: field plots
jupyter notebook notebooks/analysis.ipynb              # Part 3: power spectra & SNR
```

A full run on the fiducial $128^2 \times 100$ lightcone (114 M halos, 2.8 GB
HDF5) takes ~35 s on a laptop when the simulation itself is skipped: ~1 s for
the three power spectra, the rest dominated by reading the halo catalogue and
rendering the catalogue-based figures. Use `--plots power snr` to skip the
catalogue entirely (it is not loaded when no figure needs it), or
`--max-halos` to cap the memory footprint.

Run all cells sequentially. All notebooks are self-contained and generate all figures inline. The simulation notebooks cache 21cmFAST outputs to disk on first run. The HPC pipeline saves simulation outputs to `outputs/lightcone_data.h5` for independent loading by the analysis notebooks.

### Figure display

Every notebook opens its imports cell with

```python
%matplotlib widget
plt.rcParams['figure.constrained_layout.use'] = True
```

`%matplotlib widget` (provided by `ipympl`) renders figures as interactive
inline canvases — pan, zoom, and cursor readout without leaving the notebook —
and requires `ipywidgets >= 8` in whichever environment serves the front end,
not just in the kernel. If figures come up blank in JupyterLab, that front-end
environment is the thing to check; falling back to `%matplotlib inline` gives
static figures and is otherwise harmless.

Constrained layout is enabled globally, so spacing is resolved at draw time and
individual figures no longer call `plt.tight_layout()`. The two are mutually
exclusive: adding a `tight_layout()` call back to a cell will disable
constrained layout for that figure and emit a warning.

## Testing

Project convention (see `CLAUDE.md`): every function in `src/` should have at
least one corresponding test, in `tests/test_<module>.py`. Run the suite with:

```bash
conda run -n 21cmfast pytest tests/ -v
```

**Current status: 142 tests, all passing in ~21 s.** No test invokes 21cmFAST —
`tests/conftest.py` writes a synthetic `lightcone_data.h5` with the same schema
(16² × 12 cells, 4 000 halos), so the whole suite runs offline.

| File | Covers |
|------|--------|
| `test_analysis.py` | Cosmology helpers, the cylindrical power-spectrum estimator (shapes, symmetry, sign, and the analytic white-noise normalisation $P = \sigma^2 V_\mathrm{cell}$), wedge geometry, photo-$z$ kernel, thermal noise, SNR, Euclid selection, effective bias |
| `test_dataio.py` | HDF5 loading, metadata accessors, halo subsampling, field/catalogue skipping, product-cache round trip and staleness detection |
| `test_figures.py` | Headless backend, NaN filling, colormap, and that every one of the 10 figure functions renders and writes a non-empty file |
| `test_pipeline.py` | CLI parsing, each stage's `auto`/`force`/`skip` behaviour (with a stub simulation script), and end-to-end runs checking figures, cache reuse, and the summary JSON |
| `test_galaxy_weighting.py` | The number- and luminosity-weighted `delta_gal` constructions: weight conservation, non-cubic grids, zero-mean normalisation, both formulas against manual recomputation, $\kappa_\mathrm{UV}$ scale-invariance, and the error paths |
| `test_conversions.py` | Mass-resolution helpers: $\bar\rho_m$ against astropy, $H_0^2$ and $L^3$ scalings, the production-grid values, and total-box mass conservation |

> `src/FOV_to_cMpc.py` and the magnitude/SFR half of `src/conversions.py` are
> still untested directly, though the conversion round-trips are exercised
> indirectly through the selection and UVLF tests. The remaining first
> candidates are the explicit identities `Muv_to_Luv` ↔ `Luv_to_Muv` and
> `sfr_to_Luv` ↔ `Luv_to_sfr`.

## References

- **Davies, Mesinger & Murray (2025)** — [arXiv:2504.17254](https://arxiv.org/abs/2504.17254) — 21cmFASTv4 discrete source model
- **Gagnon-Hartman, Davies & Mesinger (2025)** — [arXiv:2502.20447](https://arxiv.org/abs/2502.20447) — galaxy–21 cm cross-correlation detection forecasts
- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — uncertainty equations, HERA noise model, and foreground wedge prescription
- **Park et al. (2019)**, MNRAS, 484, 933 — [arXiv:1809.08995](https://arxiv.org/abs/1809.08995) — 21cmFAST source model parameterisation (stellar-halo relation, SFR timescale `t_STAR × t_H`)
- **Euclid Collaboration (2022)** — [arXiv:2108.01201](https://arxiv.org/abs/2108.01201) — Euclid survey specifications
- **Bouwens et al. (2021)**, AJ, 162, 47 — [arXiv:2102.07775](https://arxiv.org/abs/2102.07775) — UV luminosity function at $z \sim 2$–$9$
- **Finkelstein et al. (2015)**, ApJ, 810, 71 — [arXiv:1410.5439](https://arxiv.org/abs/1410.5439) — UV luminosity function at $z \sim 4$–$8$
- **Speagle et al. (2014)**, ApJS, 214, 15 — [arXiv:1405.2041](https://arxiv.org/abs/1405.2041) — star-forming main sequence calibration (Eq. 28)
- **Schreiber et al. (2015)**, A&A, 575, A74 — [arXiv:1409.5433](https://arxiv.org/abs/1409.5433) — star-forming main sequence at high redshift
- **Song et al. (2016)**, ApJ, 825, 5 — [arXiv:1507.05636](https://arxiv.org/abs/1507.05636) — stellar mass–UV magnitude relation at $z \sim 4$–$8$
- **González et al. (2010)**, ApJ, 713, 115 — [arXiv:0909.3517](https://arxiv.org/abs/0909.3517) — stellar mass density and sSFR at $z \sim 7$ (constant-SFH SED fitting)
- **González et al. (2011)**, ApJL, 735, L34 — [arXiv:1008.3901](https://arxiv.org/abs/1008.3901) — galaxy stellar mass functions and $M_\star$–$L_\mathrm{UV}$ relation at $z \sim 4$–$7$
- **Lidz et al. (2009)**, ApJ, 690, 252 — [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — physical signal power spectra
- **Bardeen, Bond, Kaiser & Szalay (1986)**, ApJ, 304, 15 — BBKS transfer function
- **Kaiser (1987)**, MNRAS, 227, 1 — redshift-space distortions
- **DeBoer et al. (2017)**, PASP — [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA instrument specifications
- **Pober et al. (2014)**, ApJ, 782, 66 — [arXiv:1310.7031](https://arxiv.org/abs/1310.7031) — "moderate" foreground model; source of the $0.1\ h\ \mathrm{Mpc}^{-1}$ wedge buffer and the 21cmSense `horizon_buffer` default
- **Parsons et al. (2012a)**, ApJ, 756, 165 — [arXiv:1204.4749](https://arxiv.org/abs/1204.4749) — beam chromaticity and delay-taper leakage extending foreground power $\sim 0.15\ h\ \mathrm{Mpc}^{-1}$ beyond the horizon
- **Planck Collaboration (2020)**, A&A, 641, A6 — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — cosmological parameters
- **Hogg (1999)** — [arXiv:astro-ph/9905116](https://arxiv.org/abs/astro-ph/9905116) — comoving distance and volume formulae
- **Oke & Gunn (1983)**, ApJ, 266, 713 — AB magnitude system
- **Madau & Dickinson (2014)**, ARA&A, 52, 415 — [arXiv:1403.0007](https://arxiv.org/abs/1403.0007) — UV luminosity–SFR calibration ($\kappa_\mathrm{UV}$, Chabrier IMF)
- **Sheth & Tormen (1999)**, MNRAS, 308, 119 — [arXiv:astro-ph/9901122](https://arxiv.org/abs/astro-ph/9901122) — halo mass function and bias formula
- **Murray, Robotham & Power (2013)**, Astron. Comput., 3, 23 — [arXiv:1306.6721](https://arxiv.org/abs/1306.6721) — `hmf` halo mass function code
