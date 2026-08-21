# The HPC Run — Complete Specification

**What this document is.** A single reference for everything the HPC run
actually does: every configuration value, every derived quantity, every
formula with its evaluated number at $z_\mathrm{obs} = 7$, every file written,
and every known inconsistency between what the code computes and what the
stored output contains.

**Verified:** 2026-08-12, against commit `2c55f98` (branch `test/methods`).
Test suite: **76 passed** (`conda run -n 21cmfast pytest tests/ -q`).

**Companion documents.**
[`uncertainty_budget.md`](uncertainty_budget.md) is the dedicated reference for
the photo-z / wedge / noise / SNR chain of §§5.2–5.5, including its audit
against the source notebook; [`PIPELINE.md`](../PIPELINE.md) is the short
version with the flowchart; [`docs/project_update.md`](project_update.md) is the
science narrative and result history; [`TODO.md`](../TODO.md) is the
outstanding work. This document is the parameter-level ground truth.

**Setting up a run on a new machine?** Start at
[§13](#13-user-defined-parameter-requirements) — every parameter that is yours
to choose, with the file and line where it is set, what the code will not
choose for you, and which edits silently do nothing until the simulation is
re-run.

---

## 0. TL;DR — what one HPC run is

One command, `bash submit_job.sh --sim force`, runs a 21cmFAST v4.1.1
lightcone in a **256 Mpc box on a 128³ grid** over a **deliberately narrow
redshift slab, $z = 6.995 \to 7.005$**, draws a **114-million-halo catalogue**
at $z_\mathrm{obs} = 7.0$, builds a galaxy overdensity field from the halo SFR,
measures the **Euclid-selected galaxy bias ($b_g = 4.744$)**, applies a
**Kaiser boost ($\beta = 0.2103$)**, then computes **2D cylindrical power
spectra on a 20 × 20 $(k_\perp, k_\parallel)$ grid**, applies **photo-$z$
damping ($\sigma_z = 0.45 \Rightarrow \sigma_r = 157.5$ Mpc)**,
**foreground-wedge excision (slope 3.151,
buffer 0.0677 Mpc⁻¹)**, and a **HERA thermal-noise / shot-noise SNR**, writes
**10 figures** and a **JSON summary**, and emails a report.

It produces **~2.8 GB** of HDF5 output and leaves a **~56 GB** 21cmFAST cache
on disk. The resulting total cross-correlation SNR is **≪ 1 σ** — this is a
smoke-test configuration, not a forecast (§11).

---

## 1. Environment

| Item | Value |
|---|---|
| Conda environment | `21cmfast` (mandatory — see `CLAUDE.md`) |
| Python | 3.11 |
| `py21cmfast` | **4.1.1** (pip, compiles C extensions at install) |
| C dependencies | `fftw`, `gsl` from conda-forge, installed **before** 21cmFAST |
| `hmf` | ≥ 3.5 (Murray, Robotham & Power 2013) — halo mass function / bias |
| Scientific stack | numpy ≥ 2.0, scipy ≥ 1.14, astropy ≥ 7.0, matplotlib ≥ 3.9, h5py ≥ 3.11 |
| Environment spec | `env.yml` (conda), `requirements.txt` (pinned pip freeze) |
| Matplotlib backend | `Agg`, forced in both `run_simulation.py` and `src/figures.py` — no display required |
| Threading | `SimulationOptions.N_THREADS = 1` (template default, **not** overridden) |
| Install notes | [`docs/INSTALL_21cmFASTv4.md`](INSTALL_21cmFASTv4.md) — CSD3 quota, conda-plugin, and FFTW-linking fixes |

On CSD3 the pip and XDG caches must be redirected to scratch before install
(`PIP_CACHE_DIR`, `XDG_CACHE_HOME`), otherwise the build fails with
`OSError: [Errno 122] Disk quota exceeded`.

---

## 2. Entry point — `submit_job.sh`

A plain shell wrapper. **It contains no `#SBATCH` directives**, so it runs in
the foreground on whatever node invokes it; to submit through SLURM, add the
cluster's directives (`--partition`, `--time`, `--account`, `--cpus-per-task`)
at the top first.

| Setting | Value | Override |
|---|---|---|
| `EMAIL_TO` | `sohinidutta97@gmail.com` | edit script |
| `JOB_NAME` | `21cm_pipeline` | `JOB_NAME=... bash submit_job.sh` |
| `PYTHON_SCRIPT` | `run_pipeline.py` | `PYTHON_SCRIPT=run_simulation.py` → simulation only |
| `CONDA_ENV` | `21cmfast` | `CONDA_ENV=...` |
| Log file | `outputs/${JOB_NAME}_$(date +%Y%m%d_%H%M%S).log` | — |
| Summary read for the email | `outputs/pipeline_summary.json` | — |
| Figure inventory read for the email | `outputs/figures` | — |

What it does, in order:

1. `mkdir -p outputs`, opens a timestamped log.
2. `source "$(conda info --base)/etc/profile.d/conda.sh"` then
   `conda activate 21cmfast`.
3. Runs `/usr/bin/time -p python run_pipeline.py "$@"`, appending stdout and
   stderr to the log. **All arguments are forwarded verbatim.**
4. Computes wall-clock runtime (`HH:MM:SS`) and CPU-hours by parsing the
   `user` and `sys` lines of `/usr/bin/time -p` and dividing by 3600.
5. Sends a `sendmail -t` report: job name, status (`COMPLETE`/`FAILED`), host,
   start/finish times, runtime, CPU-hours, command, log path, exit code, the
   figure count and listing, and the **last 40 lines** of the log.
6. Exits with the Python exit code.

> No `.log` files are currently retained in `outputs/`, so no measured
> wall-clock or CPU-hour figures from a real cluster run exist in the repo.
> Timings quoted below are per-stage estimates from the code paths and file
> sizes, not from a recorded job.

---

## 3. Driver — `run_pipeline.py`

### 3.1 Command-line interface (all defaults)

| Flag | Default | Meaning |
|---|---|---|
| `--sim` | `auto` | `auto` = run `run_simulation.py` only if the HDF5 is missing; `force` = always; `skip` = never (errors if absent) |
| `--analysis` | `auto` | `auto` = recompute spectra only if the cache is missing or older than the simulation; `force` = always; `skip` = require the cache |
| `--plots` | `all` | any of `all none fields halos scaling euclid power snr budget bias` |
| `--data` | `outputs/lightcone_data.h5` | simulation HDF5 |
| `--products` | `outputs/analysis_products.h5` | power-spectrum cache |
| `--figdir` | `outputs/figures` | figure directory |
| `--summary` | `outputs/pipeline_summary.json` | summary JSON |
| `--sim-script` | `run_simulation.py` | Stage-1 script |
| `--format` | `png` | `png`, `pdf`, `svg` |
| `--dpi` | `200` | applied via `figures.apply_plot_style(dpi)` |
| `--max-halos` | `0` (= all) | uniform stride subsampling of the catalogue; number densities rescaled by `halo_sampling_factor` |
| `--m-uv-bright` | `-22.0` | bright-end Euclid magnitude cut |
| `--sigma-z` | HDF5 attr (0.45) | absolute photo-$z$ uncertainty $\sigma_z$ — **not** $\sigma_z/(1+z)$ |
| `--wedge-buffer` | HDF5 attr (0.0677) | foreground-wedge margin [Mpc⁻¹] |
| `--integration-time` | HDF5 attr (3.6 × 10⁶) | HERA integration time [s] |
| `--bandwidth` | HDF5 attr (8 × 10⁶) | per-band bandwidth [Hz] |
| `--quiet` | off | suppress progress output |

The last four are **uncertainty-budget overrides**. Each resolves
CLI flag → HDF5 root attribute → hardcoded default, and none of them affects
the simulated fields, so all four can be swept from cached spectra in seconds
without `--sim force`. See [`uncertainty_budget.md`](uncertainty_budget.md) §5.

Cache staleness is decided by `src.dataio.products_are_stale`, comparing the
`source_mtime` attribute of `analysis_products.h5` against the mtime of
`lightcone_data.h5`.

The halo catalogue is loaded **only** if `halos`, `scaling`, or `bias` is in
the requested plot groups (`needs_catalog`) — `--plots power snr` therefore
skips a 2.7 GB read entirely.

### 3.2 Failure handling

`FileNotFoundError`, `RuntimeError`, and `ValueError` are caught in `main()`,
printed to stderr as `ERROR: ...`, and returned as exit code 1 — which
`submit_job.sh` turns into a `FAILED` email. Anything else propagates as a
traceback.

---

## 4. Stage 1 — `run_simulation.py`

Invoked as a subprocess (`sys.executable run_simulation.py`, `cwd=REPO_ROOT`).
This is the only expensive stage.

### 4.1 Configuration block (verbatim, lines 89–164)

**Grid**

| Parameter | Value |
|---|---|
| `HII_DIM` | 128 |
| `BOX_LEN` | 256.0 Mpc |
| `DIM` | `3 × HII_DIM` = **384** |
| Transverse cell size | 256/128 = **2.0 Mpc** |
| High-res IC cell size | 256/384 = **0.6667 Mpc** |

**Redshift range**

| Parameter | Value |
|---|---|
| `z_min` | **6.995** |
| `z_max` | **7.005** |
| `z_obs` | 7.0 (midpoint; used for noise, wedge, bias, RSD) |
| `minimum_los_slices` | 100 |

**Euclid survey**

| Parameter | Value |
|---|---|
| `M_UV_limit` (faint) | **−18** |
| `M_UV_bright` (§4 of the script) | **−22** |
| `photoz_uncertainty` $\sigma_z$ | **0.45** (absolute, not $\sigma_z/(1+z)$ — see §11.8) |
| `mean_galaxy_density` $\bar n$ | **3 × 10⁻³ h³ Mpc⁻³** |
| `galaxy_bias` (fallback only) | 8 |

**Cosmology (Planck 2018)** — note these are the *script's* constants; the
comoving distances to the lightcone endpoints use `astropy.cosmology.Planck18`
directly when 21cmFAST is available.

| Parameter | Value |
|---|---|
| `OMEGA_M_0` | 0.315 |
| `HUBBLE_CONSTANT` | 67.36 km s⁻¹ Mpc⁻¹ ($h = 0.6736$) |
| `OMEGA_B_0` (§4 analytic path only) | 0.049 |
| `SPEED_OF_LIGHT_KMS` / `_MPS` | 3 × 10⁵ / 3 × 10⁸ |
| `F_21_MHZ` | 1420.405 MHz |

**HERA instrument**

| Parameter | Value |
|---|---|
| `HERA_DISH_DIAMETER` | 14.0 m |
| `integration_time` | 1000 h = **3.6 × 10⁶ s** |
| `bandwidth` | **8 MHz** |

**Wedge and binning**

| Parameter | Value |
|---|---|
| `wedge_buffer` | **0.0677 Mpc⁻¹** = 0.1 $h$ Mpc⁻¹ at $h = 0.6766$ (Pober et al. 2014 "moderate"; 21cmSense `horizon_buffer` default) |
| `n_bins_perp` × `n_bins_parallel` | 20 × 20 |
| Output | `outputs/lightcone_data.h5` |

### 4.2 Derived geometry (§1 of the script)

Computed with `Planck18.comoving_distance`:

| Quantity | Value |
|---|---|
| $D_c(z_\mathrm{min} = 6.995)$ | 8809.810 Mpc |
| $D_c(z_\mathrm{max} = 7.005)$ | 8813.310 Mpc |
| $L_\mathrm{LOS}$ (requested) | **3.4998 Mpc** |
| Natural $N_z = \mathrm{round}(L/2.0)$ | 2 → **floored to 100** by `minimum_los_slices` |
| Actual LOS slice spacing (stored `lc_dist_Mpc`) | **0.03539 Mpc** (~57× oversampled vs the 2 Mpc transverse cell) |
| `lc_redshifts` | `np.linspace(6.995, 7.005, 100)` — uniform in $z$, **not** in comoving distance |
| Node redshifts | $n = \max(\mathrm{round}(10\Delta z), 5) =$ **5**, `np.linspace(7.005, 6.995, 5)` (high-$z$ → low-$z$) |

**Mass resolution** (`src.conversions.cell_mass`, $\bar\rho_m = \Omega_m
\rho_{\mathrm{crit},0} = 3.967 \times 10^{10}\ M_\odot\ \mathrm{Mpc}^{-3}$):

| Grid | Cell | Mass resolution |
|---|---|---|
| `DIM` = 384 (initial conditions) | 0.6667 Mpc | **1.175 × 10¹⁰ M☉** |
| `HII_DIM` = 128 (ionisation, 21 cm) | 2.0 Mpc | **3.173 × 10¹¹ M☉** |
| Halo sampler floor `SAMPLER_MIN_MASS` | — | **1 × 10⁸ M☉** |

### 4.3 21cmFAST invocation (§2)

```python
inputs = p21c.InputParameters.from_template(["simple"], random_seed=42)
inputs = inputs.clone(
    node_redshifts=node_redshifts,                     # 5 nodes, 7.005 → 6.995
    simulation_options={"HII_DIM": 128, "BOX_LEN": 256.0, "DIM": 384},
    matter_options={"USE_INTERPOLATION_TABLES": "hmf-interpolation"},
)
lightconer = p21c.RectilinearLightconer(
    lc_redshifts=lc_redshifts,                          # 100 slices
    quantities=("brightness_temp", "density", "neutral_fraction", "halo_sfr"),
)
lightcone = p21c.run_lightcone(
    lightconer=lightconer, inputs=inputs,
    include_dvdr_in_tau21=False,                        # no dv/dr in tau_21
    apply_rsds=False,                                   # RSDs applied analytically in §5
)
```

**Random seed: 42** — the run is fully reproducible.

**"simple" template parameters actually in force** (verified by instantiating
the template in the environment):

*Source model (`AstroParams`)*

| Parameter | Value | Meaning |
|---|---|---|
| `F_STAR10` | −1.3 (i.e. 0.05) | stellar fraction at 10¹⁰ M☉ |
| `ALPHA_STAR` | 0.5 | SHMR power-law slope |
| `M_TURN` | 8.7 (5.01 × 10⁸ M☉) | turnover mass |
| `F_ESC10` / `ALPHA_ESC` | −1.0 (0.1) / −0.5 | escape fraction |
| `t_STAR` | **0.5** | SFR timescale as a fraction of $t_H(z)$ |
| `SIGMA_STAR` | 0.25 dex | log-normal stellar-mass scatter |
| `SIGMA_SFR_LIM` / `SIGMA_SFR_INDEX` | 0.19 dex / −0.12 | SFR scatter |
| `HII_EFF_FACTOR` | 30.0 | ionising efficiency |
| `R_BUBBLE_MAX` / `R_BUBBLE_MIN` | 15.0 / 0.6204 Mpc | excursion-set radii |
| `CLUMPING_FACTOR` | 2.0 | |
| `L_X` / `NU_X_THRESH` | 40.5 / 500 eV | X-ray heating |

*Matter options*

| Option | Value |
|---|---|
| `SOURCE_MODEL` | **CHMF-SAMPLER** (discrete halos, Davies et al. 2025) |
| `SAMPLE_METHOD` | MASS-LIMITED |
| `HMF` | ST (Sheth-Tormen) |
| `POWER_SPECTRUM` | EH (Eisenstein & Hu) |
| `PERTURB_ALGORITHM` | 2LPT |
| `FILTER` / `HALO_FILTER` | spherical-tophat |
| `USE_INTERPOLATION_TABLES` | **hmf-interpolation** (explicitly set) |
| `PERTURB_ON_HIGH_RES` | False |

*Astro options*

| Option | Value | Consequence |
|---|---|---|
| `USE_TS_FLUCT` | **False** | no spin-temperature fluctuations — the main speed-up |
| `USE_MINI_HALOS` | False | |
| `INHOMO_RECO` | False | |
| `HII_FILTER` | sharp-k | |
| `PHOTON_CONS_TYPE` | no-photoncons | |
| `USE_CMB_HEATING` / `USE_LYA_HEATING` / `USE_X_RAY_HEATING` | True | |

**Fields retrieved:** `brightness_temp` [mK], `density`, `neutral_fraction`,
`halo_sfr`, each shape **(128, 128, 100)**.

After the run the script **overwrites** its own geometry from the lightcone
object: `N_z = lightcone.n_slices`, `L_los = lightcone.lightcone_dimensions[2]`,
`lc_dist_Mpc = lightcone.lightcone_distances`. See §11.1 — this is where the
200 Mpc / 3.5 Mpc discrepancy enters.

**Synthetic fallback.** If `py21cmfast` is not importable, the script generates
a synthetic lightcone instead (seed 42): a $P(k) \propto k^{-2}$ Gaussian
density field normalised to $\sigma = 0.5$, a 50th-percentile inside-out
ionisation field smoothed at $k_\mathrm{bubble} = 2\pi/10\ \mathrm{Mpc}^{-1}$,
$\delta T_b = 27\sqrt{(1+z)/10}\, x_\mathrm{HI}(1+\delta)$ mK, and a
biased Poisson galaxy field. `sampler_min_mass` is then `NaN`. This path exists
for CI and laptops; **it is not what runs on the cluster.**

### 4.4 Halo catalogue (§3a)

```python
initial_conditions = p21c.compute_initial_conditions(inputs=inputs)
halo_catalog       = p21c.determine_halo_catalog(redshift=7.0, inputs=inputs, ...)
perturbed_halos    = p21c.perturb_halo_catalog(inputs=inputs, ...)  # Lagrangian → Eulerian
```

**Unit correction.** py21cmfast v4 computes SFR internally in $M_\odot\
\mathrm{s}^{-1}$ (`scaling_relations.c`: `sfr = M_star / (t_STAR · t_H)` with
$t_H$ in seconds). The script multiplies by
`_SEC_PER_YR = 365.25 × 24 × 3600 = 3.15576 × 10⁷` so the stored catalogue is
in $M_\odot\ \mathrm{yr}^{-1}$. (This is the fix documented in
[`docs/Low_SFR_fix.md`](Low_SFR_fix.md).)

Recorded catalogue, at $z_\mathrm{obs} = 7.0$, seed 42, 256 Mpc box:

| Quantity | Value |
|---|---|
| Total halos | **114,291,212** |
| Halo mass range | 1.0 × 10⁸ – 1.77 × 10¹² M☉ |
| Stellar mass range | 29 – 2.79 × 10¹¹ M☉ |
| SFR range | 2.2 × 10⁻¹⁰ – 488 M☉ yr⁻¹ |
| SFR median | 1.17 × 10⁻⁵ M☉ yr⁻¹ |
| On-disk size (4 arrays, float32) | **2.74 GB** |

### 4.5 Galaxy overdensity field (§3b)

$$\delta_\mathrm{gal}(\mathbf{x}) = \frac{\mathrm{SFR}(\mathbf{x})}{\langle \mathrm{SFR}\rangle} - 1$$

built from the lightcone's `halo_sfr` field, using a **single global mean**
over the whole box (see §11.3). Shape (128, 128, 100).

### 4.6 Galaxy bias (§4) — two estimators

Both are computed; the catalogue estimator is **adopted** whenever a catalogue
exists.

**Common calibration chain** (`src/conversions.py`, Madau & Dickinson 2014 +
Oke & Gunn 1983):

$$L_\mathrm{UV} = \mathrm{SFR}/\kappa_\mathrm{UV},\quad
\kappa_\mathrm{UV} = 1.15\times10^{-28},\qquad
M_\mathrm{UV} = 51.60 - 2.5\log_{10} L_\mathrm{UV}$$

The Euclid window $-22 \le M_\mathrm{UV} \le -18$ is therefore equivalent to
**SFR ∈ [0.7956, 31.674] M☉ yr⁻¹**.

**Star-formation timescale** (`src.analysis.star_formation_timescale`,
Park et al. 2019 Eq. 3):

$$t_\mathrm{sf} = t_\star\, t_H(z) = 0.5 \times 1.1406\ \mathrm{Gyr} =
\mathbf{570.3\ Myr}\ \ (z = 7)$$

| Estimator | $b_g$ | Basis |
|---|---|---|
| **Halo catalogue (adopted)** | **4.744** (range 2.83 – 9.14) | Per-halo $M_\mathrm{UV}$ from the catalogue's own SFR → Euclid window → mean Sheth-Tormen bias over survivors. Inherits 21cmFAST's log-normal scatter. Halos selected: **49,315** (~0.04 % of the catalogue) |
| Analytic HMF integral (cross-check) | 5.39 | `hmf.MassFunction(z=7, Mmin=7, Mmax=13, dlog10m=0.02)`, Sheth-Tormen bias weighted by `dndlog10m` over the *mean, scatter-free* relation $f_\star = 0.05 (M_h/10^{10})^{0.5} e^{-5\times10^8/M_h}$, $M_\star = f_\star (\Omega_b/\Omega_m) M_h$ |

Bias formula (`src.conversions.sheth_tormen_bias`, $\delta_c = 1.686$,
$a = 0.707$, $p = 0.3$), taking `hmf`'s **already-squared** peak height:

$$b(\tilde\nu) = 1 + \frac{a\tilde\nu - 1}{\delta_c} +
\frac{2p}{\delta_c\left(1 + (a\tilde\nu)^p\right)},\qquad
\tilde\nu = (\delta_c/\sigma)^2$$

Squaring `mf.nu` a second time is what produced the historical
$b_g = 33.39$; that path is gone (see [`project_update.md`](project_update.md)
§4).

### 4.7 Kaiser RSD (§5)

| Quantity | Value at $z_\mathrm{obs} = 7$ |
|---|---|
| $\Omega_m(z)$ | 0.99577 |
| $f = \Omega_m(z)^{0.55}$ | **0.99767** |
| $\beta = f / b_g$ | **0.21030** |
| Max boost $(1+\beta)$ at $\mu = 1$ | **1.210×** |

Applied in Fourier space, $\delta^s_\mathrm{gal}(\mathbf{k}) =
(1 + \beta\mu^2)\,\delta_\mathrm{gal}(\mathbf{k})$, on a **non-cubic** $k$-grid:
transverse spacing 2.0 Mpc, LOS spacing `dz_cell = L_los / N_z` (§11.1).

### 4.8 HDF5 written (§6)

`outputs/lightcone_data.h5`, gzip level 4 on the four fields.

| Dataset | Shape | dtype | Unit |
|---|---|---|---|
| `brightness_temp_field` | (128, 128, 100) | float32 | mK |
| `density_field` | (128, 128, 100) | float32 | overdensity |
| `neutral_fraction` | (128, 128, 100) | float32 | [0, 1] |
| `galaxy_overdensity` | (128, 128, 100) | float64 | dimensionless, Kaiser-boosted |
| `lc_redshifts` | (100,) | float64 | — |
| `lc_dist_Mpc` | (100,) | float64 | Mpc |
| `halo_catalog/halo_masses` | (114291212,) | float32 | M☉ |
| `halo_catalog/halo_coords` | (114291212, 3) | float32 | Mpc |
| `halo_catalog/stellar_masses` | (114291212,) | float32 | M☉ |
| `halo_catalog/sfr` | (114291212,) | float32 | M☉ yr⁻¹ |

**Root attributes written by the current code** (34; `galaxy_bias_hmf_analytic`
only when the analytic cross-check succeeds): `HII_DIM`, `DIM`,
`BOX_LEN`, `N_z`, `L_los`, `cell_size`, `hires_cell_size`, `M_cell_hires`,
`M_cell_lores`, `sampler_min_mass`, `z_min`, `z_max`, `z_obs`, `galaxy_bias`,
`galaxy_bias_method`, `galaxy_bias_hmf_analytic`, `t_STAR`, `sfr_timescale_yr`,
`beta_rsd`, `mean_galaxy_density`, `photoz_uncertainty`, `M_UV_limit`,
`OMEGA_M_0`, `HUBBLE_CONSTANT`, `SPEED_OF_LIGHT_KMS`, `SPEED_OF_LIGHT_MPS`,
`F_21_MHZ`, `F_21_HZ`, `HERA_DISH_DIAMETER`, `integration_time`, `bandwidth`,
`wedge_buffer`, `n_bins_perp`, `n_bins_parallel`.

---

## 5. Stage 2 — analysis (`src/analysis.py`)

### 5.1 Cylindrical power spectra

`compute_all_power_spectra` subtracts a **single global scalar mean** from
$\delta T_b$ (the galaxy field is already an overdensity), then calls
`compute_cylindrical_cross_power` three times: $P_{21}$, $P_\mathrm{gal}$,
$P_{21\times\mathrm{gal}}$.

Estimator, for a $(N, N, N_z)$ box with $dx = L_\perp/N$, $dz = L_\parallel/N_z$,
$V = L_\perp^2 L_\parallel$:

$$P(\mathbf{k}) = \frac{\big[\tilde{A}(\mathbf{k})\,\tilde{B}^*(\mathbf{k})\big]_\mathrm{Re}}{V},
\qquad \tilde{A} = \mathrm{FFT}(A)\cdot dx^2 dz$$

Binning: **log-spaced**, from $0.5\,\Delta k$ to $1.05\,k_\mathrm{Nyq}$
($\times\sqrt2$ transverse), 20 × 20 bins, geometric bin centres. Empty bins
are `NaN`.

Grid actually produced (with $L_\perp = 256$ Mpc, $L_\parallel = 200$ Mpc — see
§11.1):

| Quantity | Value |
|---|---|
| $\Delta k_\perp = 2\pi/256$ | 0.02454 Mpc⁻¹ |
| $\Delta k_\parallel = 2\pi/200$ | 0.03142 Mpc⁻¹ |
| $k_\mathrm{Nyq}$ (both axes) | $\pi/2$ = 1.5708 Mpc⁻¹ |
| $k_\perp$ bin centres | 0.0140 → 2.0457 Mpc⁻¹ |
| $k_\parallel$ bin centres | 0.0176 → 1.4682 Mpc⁻¹ |
| Fourier modes binned | 1,621,917 of 1,638,400 |
| Empty bins | **145 / 400** |

Cached to `outputs/analysis_products.h5` (21 KB) with `source_path` and
`source_mtime` attributes for staleness detection.

### 5.2 Photo-$z$ damping

> §§5.2–5.5 are one calculation, assembled by
> `analysis.compute_uncertainty_budget` and returned as an `UncertaintyBudget`.
> `run_pipeline.observational_stage` holds no physics — it resolves parameters
> and calls that function once. Full treatment, including the term-by-term
> audit against `21cmfast_HERAxEuclid_lightcone.ipynb`, in
> [`uncertainty_budget.md`](uncertainty_budget.md).

$$\sigma_r = \frac{c\,\sigma_z}{H(z_\mathrm{obs})} =
\frac{3\times10^5 \times 0.45}{857.26} = \mathbf{157.48\ Mpc},
\qquad W(k_\parallel) = e^{-k_\parallel^2\sigma_r^2/2}$$

Applied as $P_\times \to P_\times W$ (one field smeared) and
$P_\mathrm{gal} \to P_\mathrm{gal} W^2$.

Evaluated on the actual $k_\parallel$ grid — with $\sigma_z = 0.45$ the kernel
is **already collapsed at the first bin**:

| $k_\parallel$ [Mpc⁻¹] | 0.0176 | 0.0223 | 0.0281 | 0.0355 | 0.0448 | 0.0565 |
|---|---|---|---|---|---|---|
| $W$ ($\sigma_z = 0.45$, current) | 0.021 | 2.1 × 10⁻³ | 5.6 × 10⁻⁵ | 1.7 × 10⁻⁷ | 1.6 × 10⁻¹¹ | 6.6 × 10⁻¹⁸ |
| $W$ ($\sigma_z = 0.059$, old) | 0.936 | 0.900 | 0.845 | 0.765 | 0.652 | 0.507 |

$W = 0.5$ now falls at $k_\parallel = 0.0075$ Mpc⁻¹, **below the smallest bin
the box can sample** (0.0176) — reaching it would need
$L_\mathrm{LOS} > 840$ Mpc. See §11.8.

### 5.3 Foreground wedge

$$m_\mathrm{horizon} = \frac{D_c(z) H(z)}{c(1+z)},\qquad
m_\mathrm{FoV} = \sin\!\left(\frac{\lambda_\mathrm{obs}}{D_\mathrm{dish}}\right) m_\mathrm{horizon}$$

| Quantity | Value at $z = 7$ |
|---|---|
| $D_c(7)$ (script cosmology) | 8821.33 Mpc |
| $H(7)$ | 857.26 km s⁻¹ Mpc⁻¹ |
| $\nu_\mathrm{obs} = 1420.405/(1+z)$ | 177.55 MHz |
| **Horizon slope** | **3.1509** |
| **HERA FoV slope** (14 m dish) | **0.37936** |
| Buffer | **0.0677 Mpc⁻¹** |

Mask: `k_par > k_perp × 3.1509 + 0.0677`. The **horizon** slope is used for
the mask (the configuration Pober et al.'s buffer was calibrated against); the
FoV slope is drawn on the figures only.

> Minor unit inconsistency: the buffer is 0.1 $h$ Mpc⁻¹ converted at
> $h = 0.6766$ (astropy's `Planck18`), while the run's own
> `HUBBLE_CONSTANT = 67.36` implies $h = 0.6736$ and hence 0.06736 Mpc⁻¹.
> A 0.4 % difference — noted for completeness, not worth acting on.

**The wedge and the photo-$z$ kernel are in direct conflict.** At the smallest
$k_\perp = 0.0140$, the wedge admits only $k_\parallel > 0.1118$ Mpc⁻¹, where
$W = 0.070$ at the old $\sigma_z = 0.059$ — every surviving mode already damped
to ≤ 7 % of its amplitude, and the structural reason the SNR is tiny (§10). At
the corrected $\sigma_z = 0.45$ the same mode has $W = 5\times10^{-68}$: the
wedge and the photo-$z$ kernel have **no overlap at all** on this $k$-grid.

### 5.4 Noise

**21 cm thermal noise** (`hera_thermal_noise_power`) — a scaling estimate, not
an instrument model:

$$T_\mathrm{sys} = 100\ \mathrm{K} + 60\ \mathrm{K}\left(\frac{300\ \mathrm{MHz}}{\nu}\right)^{2.55},
\qquad P_N = \frac{T_\mathrm{sys}^2 \times 10^3}{t_\mathrm{int}\,\Delta\nu}$$

| Quantity | Value |
|---|---|
| $T_\mathrm{sys}$ at 177.55 MHz | **328.6 K** |
| $t_\mathrm{int}$ | 3.6 × 10⁶ s (1000 h) |
| $\Delta\nu$ | 8 × 10⁶ Hz |
| $P_{N,21}$ | **3.7488 mK² Mpc³** |

**Galaxy shot noise:** $P_{N,\mathrm{gal}} = 1/\bar n = 1/(3\times10^{-3}) =
\mathbf{333.33\ Mpc^3}$.

### 5.5 SNR (La Plante et al. 2023, Eqs. 15–17)

$$\sigma_\times = \sqrt{\tfrac12\left[P_\times^2 +
\left(|P_{21}| + P_{N,21}\right)\left(|P_\mathrm{gal}| + P_{N,\mathrm{gal}}\right)\right]},
\qquad \mathrm{SNR} = \frac{|P_\times|}{\sigma_\times}$$

Total = quadrature sum over bins outside the wedge (`np.nansum`).
**`mode_counts` is *not* divided into the variance** — see §11.4.

$\sigma_\times^2$ is stored as its two halves, `cosmic_variance_term`
($\tfrac12 P_\times^2$) and `noise_coupling_term`
($\tfrac12\sigma_{21}\sigma_\mathrm{gal}$); their ratio
(`cosmic_variance_fraction`) says whether a mode is sample-variance or noise
limited. For the stored run it is $2\times10^{-224}$ — entirely noise
dominated, because photo-$z$ damping has erased $P_\times$ wherever the wedge
admits it.

La Plante's equations carry factors of $T_0(z)$ that are absent here. They
**cancel exactly** in the SNR ratio, so the omission is correct as long as
$\sigma_\times$ is not quoted as a standalone error bar —
[`uncertainty_budget.md`](uncertainty_budget.md) §2.5 gives the algebra.

### 5.6 Euclid selection and effective bias (`bias_stage`)

Re-runs `select_euclid_halos` + `effective_galaxy_bias` on the loaded
catalogue with `M_UV_bright = -22` (CLI) and `M_UV_faint` from the HDF5
`M_UV_limit` attribute. `hmf` grid: `dlog10m = 0.02`, mass range
$\lfloor\log_{10}M_\mathrm{min}\rfloor - 0.5$ to
$\lceil\log_{10}M_\mathrm{max}\rceil + 0.5$ in $M_\odot/h$
($h = H_0/100 = 0.6736$), bias interpolated onto the selected halos with
`scipy.interpolate.interp1d(..., fill_value="extrapolate")`.

---

## 6. Stage 3 — figures (`src/figures.py`)

`Agg` backend, `dpi = 200` (both `figure.dpi` and `savefig.dpi`), format `png`.

| Group | Figures | Needs catalogue |
|---|---|---|
| `fields` | `lightcone_fields`, `lightcone_slice` | no |
| `halos` | `halo_catalogue`, `sfr_relations` | yes |
| `scaling` | `uv_luminosity_function`, `stellar_mass_muv`, `main_sequence`, `uv_selection_maps` | yes |
| `euclid` | `euclid_selected_catalogue`, `selected_galaxy_overdensity`, `galaxy_overdensity_on_21cm` | yes |
| `power` | `power_spectra_2d`, `galaxy_wedge`, `wedge_real_space` | no |
| `snr` | `cross_snr` | no |
| `budget` | `uncertainty_budget`, `photoz_suppression` | no |
| `bias` | `galaxy_bias` | yes |

**18 figures total.**

The `euclid` group does not read the stored `galaxy_overdensity`. It rebuilds
the field with `analysis.galaxy_overdensity_from_catalogue()` on
`run_simulation.py` §3b's grid (`n_perp = HII_DIM`, `n_los = N_z`,
`los_extent = BOX_LEN` — the *coeval* box, not `L_los`), applying the Euclid
window, because the stored field is written under the default
`GALAXY_WEIGHTING = "lightcone_sfr"` and carries no magnitude cut at all. The
per-halo weight is `--galaxy-weighting` (`number` by default, or `luminosity`).
The deposit is done once per run and shared by both overdensity figures.

Two consequences worth knowing before reading these plots:

- **The selected sample is shot-noise dominated on this grid.** At the
  fiducial parameters 49,315 of 114 M halos survive the cut, i.e. 0.03
  galaxies per cell over $128^2\times100$ cells; ~97.5 % of cells are empty
  and a single transverse slice is almost pure noise. Hence the LOS
  projection panel, the slab averaging in the overlay
  (`slab_cells = 8`), and the display-only Gaussian smoothing of the
  contoured field (`smooth_cells = 2`).
- **The overlays are transverse only.** The catalogue is coeval and the 21 cm
  field is a lightcone, so the two share the $(x, y)$ plane and the array
  shape but not an LOS scale — the same mismatch `run_simulation.py` already
  accepts in catalogue mode (§11 and `PIPELINE.md`). The third panel's Pearson
  $r$ pairs cells exactly as `compute_all_power_spectra` does, so it is the
  real-space counterpart of the cross-power sign, mismatch included.

Literature overlays hardcoded in the figure code:

| Figure | Overlay | Parameters |
|---|---|---|
| UVLF | Bouwens et al. (2021), Table 5, $z = 6.8$ | $\phi_\star = 1.9\times10^{-4}$, $M_\star = -21.15$, $\alpha = -2.06$ |
| UVLF | Finkelstein et al. (2015), Table 4, $z = 7$ | $\phi_\star = 1.57\times10^{-4}$, $M_\star = -21.03$, $\alpha = -2.03$ |
| $M_\star$–$M_\mathrm{UV}$ | Song et al. (2016) | $\log_{10}M_\star = 8.86 - 0.5(M_\mathrm{UV} + 20)$ |
| $M_\star$–$M_\mathrm{UV}$ | González et al. (2010) | $\log_{10}M_\star = 9.06 - 0.5(M_\mathrm{UV} + 20)$ |
| Main sequence | Speagle et al. (2014) Eq. 28 | $(0.84 - 0.026t)\log_{10}M_\star - (6.51 - 0.11t)$, $t$ from `Planck18.age(7)` |
| Main sequence | Schreiber et al. (2015) high-$z$ | $\log_{10}\mathrm{SFR} = \log_{10}M_\star - 8$ (sSFR = 10 Gyr⁻¹) |
| Main sequence | 21cmFAST model (no free parameters) | $\log_{10}\mathrm{SFR} = \log_{10}M_\star - \log_{10}(t_\star t_H)$ → sSFR = 1.75 Gyr⁻¹, 0.76 dex below Schreiber+15 |

When `--max-halos` is used, UVLF counts are rescaled by
`halo_sampling_factor` and the panel is annotated with the stride.

---

## 7. Stage 4 — summary JSON

`outputs/pipeline_summary.json`, written by `build_summary` / `write_summary`
and echoed by `print_report`. Keys: `generated`, `data_file`,
`ran_simulation`, `recomputed_power_spectra`, and the blocks `simulation`
(15 fields), `power_spectra` (7), **`uncertainty_budget` (21)**,
`euclid_selection` (6), `effective_galaxy_bias` (4), `figures` (absolute
paths).

`uncertainty_budget` is `UncertaintyBudget.as_dict()` — every scalar of the
budget, listed in [`uncertainty_budget.md`](uncertainty_budget.md) §6.2. The
older `observation` block (8 fields) is **retained as an alias** of the same
values so existing notes and scripts keep working.

---

## 8. Disk footprint

| Path | Size | Note |
|---|---|---|
| `outputs/lightcone_data.h5` | **2.76 GB** | 2.74 GB of it is the halo catalogue; the four gzip-compressed fields are ~46 MB |
| `outputs/analysis_products.h5` | 54 KB | 3 × (20 × 20) spectra + mode counts + $k$-grids, plus the 10-map `uncertainty_budget` group |
| `outputs/figures/*.png` | ~6 MB | 15 files; `halo_catalogue.png` alone is 3.5 MB |
| `outputs/pipeline_summary.json` | 1.4 KB | |
| `d1f8b93ecb5e05f9040e32ca2a1534a2/42/` | **56 GB** | 21cmFAST cache, written to the **project root** (gitignored) |

The cache dominates and is easy to miss when sizing a scratch quota. Current
contents:

| Cache entry | Size |
|---|---|
| `InitialConditions.h5` | 920 MB |
| `039ee.../` — 5 nodes at $z$ = 6.9950–7.0050 (the smoke-test run) | 18 GB |
| `ad4fd.../` — 10 nodes at $z$ = 6.5000–7.5000 (the briefly-set production range) | 36 GB |
| Per node redshift, e.g. `7.0000/` | 3.6 GB (`HaloCatalog.h5` 3.83 GB uncompressed + `PerturbedField.h5` 17 MB) |

**Per-node cost is ~3.6 GB**, so the production range (10 nodes) costs ~36 GB
of cache on top of the ICs. Budget scratch accordingly.

---

## 9. Reproducing a run

```bash
conda activate 21cmfast          # mandatory

# Full run on the cluster, with timing + email
bash submit_job.sh --sim force

# Analysis + figures only, from the stored HDF5
bash submit_job.sh

# Simulation only (old behaviour)
PYTHON_SCRIPT=run_simulation.py bash submit_job.sh

# Selected variants
python run_pipeline.py --analysis force        # recompute the spectra
python run_pipeline.py --plots power snr       # k-space figures only (skips the 2.7 GB catalogue read)
python run_pipeline.py --plots none            # numbers only
python run_pipeline.py --max-halos 5000000     # cap catalogue memory
python run_pipeline.py --format pdf --dpi 300  # publication output

# Tests
conda run -n 21cmfast pytest tests/ -v         # 76 tests
```

---

## 10. Last recorded results

From `outputs/pipeline_summary.json`, generated **2026-08-04 15:54:56** against
the HDF5 written **2026-06-15** (`ran_simulation: false`,
`recomputed_power_spectra: false`).

> **These numbers are superseded on three counts.** The stored HDF5 predates
> the galaxy-bias fix ($b_g = 33.39$, $\beta = 0.0299$ instead of 4.744 /
> 0.2103, so `galaxy_overdensity` carries the wrong Kaiser boost), and the
> summary predates both the wedge-buffer change (0.02 → 0.0677 Mpc⁻¹, §11.7)
> and the $\sigma_z$ correction (0.059 → 0.45, §11.8 — which alone takes
> $\sigma_r$ from 20.6 to 157.5 Mpc and the SNR to zero). They are recorded
> here because they are the only end-to-end numbers that exist.

| Quantity | Value |
|---|---|
| Box | 256 Mpc, 128² × 100 cells, `L_los` 200.0 Mpc, cell 2.0 Mpc |
| Redshifts | 6.995 → 7.005, $z_\mathrm{obs} = 7.0$ |
| $\langle x_\mathrm{HI}\rangle$ | **0.1764** |
| $\langle \delta T_b\rangle$ | 3.152 mK |
| `galaxy_bias` (stale) | 33.389 |
| `beta_rsd` (stale) | 0.029880 |
| Halos in file / loaded | 114,291,212 / all (`halo_sampling_factor` 1.0) |
| $k_\perp$ range | 0.013992 – 2.0457 Mpc⁻¹ |
| $k_\parallel$ range | 0.017646 – 1.4682 Mpc⁻¹ |
| Empty bins | 145 / 400 |
| Large-scale $P_\times$ mean (first 5 × 5 bins) | **−5644.3** → anti-correlated ✓ |
| $\sigma_r$ | 20.647 Mpc |
| Horizon / FoV slope | 3.1509 / 0.37936 |
| Modes outside wedge (buffer 0.02) | 105 / 400 = 26.25 % |
| $P_{N,21}$ / $P_{N,\mathrm{gal}}$ | 3.7488 mK² Mpc³ / 333.33 Mpc³ |
| **Total SNR** | **0.0629 σ** — no detection |

**Effect of the wedge-buffer change alone**, recomputed here from the cached
(stale-bias) spectra:

| Buffer | Modes outside wedge | Total SNR |
|---|---|---|
| 0.02 Mpc⁻¹ (old) | 105 / 400 = 26.2 % | 0.0629 σ |
| **0.0677 Mpc⁻¹ (current)** | **97 / 400 = 24.2 %** | **0.0048 σ** |

Only 8 bins are lost, but the SNR falls **13×** — the excised bins are the low
-$k_\parallel$ ones where the photo-$z$ kernel had not yet collapsed (§5.2).
The surviving modes contribute almost nothing.

---

## 11. Known inconsistencies and caveats

Ordered by how much they affect the numbers above.

### 11.1 `L_los` is recorded as 200 Mpc; the data spans 3.5 Mpc

The single most important thing to know about the current output.

| Source | LOS extent | Implied slice spacing |
|---|---|---|
| `lc_dist_Mpc` (the actual slice distances) | 8810.002 → 8813.502 = **3.4999 Mpc** | **0.035385 Mpc** |
| `L_los` attribute, from `lightcone.lightcone_dimensions[2]` | **200.0 Mpc** | 2.0 Mpc |

The attribute equals $N_z \times$ *transverse* cell size (100 × 2.0), not the
distance actually spanned by the slices — the two disagree by a factor
**56.5**. Everything downstream uses the 200 Mpc value:
`compute_cylindrical_cross_power(box_len_los=data.L_los)` sets
$dz = 2.0$ Mpc and $\Delta k_\parallel = 2\pi/200$, and §5 of
`run_simulation.py` builds its Kaiser $k_z$ grid the same way. So the whole
$k_\parallel$ axis — and with it the wedge mask, the photo-$z$ kernel argument,
and the RSD $\mu$ — is scaled by that factor relative to the sampling of the
stored data. `CHANGELOG.md` records the 200 Mpc value as "after rounding to
cell boundaries"; the 0.0354 Mpc spacing in `lc_dist_Mpc` shows it is not a
rounding effect. **Resolve this before any $k_\parallel$ number is quoted as
physical.** It is not currently tracked in `TODO.md`.

### 11.2 The redshift slab is a smoke test, on purpose

$\Delta z = 0.01$ gives a quasi-coeval slab, which is what the estimator
assumes. The 3.5 Mpc extent means essentially no independent $k_\parallel$
modes exist, and the 100 slices are interpolated from only **5** node
redshifts — LOS structure below the 2 Mpc cell is interpolation, not signal.
Widening to $z = 6.5$–7.5 ($L_\mathrm{LOS} = 350.83$ Mpc, $N_z = 175$ natural,
10 nodes, 22.28 MHz) is gated on `TODO.md` §P0.1–P0.4. Do **not** widen it
first.

### 11.3 Estimator issues that P0 will fix

| # | Issue | Measured impact |
|---|---|---|
| P0.1 | `lc_redshifts` is uniform in $z$, not in comoving distance | 0.19 % cell-size spread now; **20.4 %** at $\Delta z = 1$ |
| P0.2 | A single global scalar mean is subtracted from $\delta T_b$ (and $\langle\mathrm{SFR}\rangle$ for $\delta_\mathrm{gal}$) | leaves a monotonic LOS ramp that aliases into low-$k_\parallel$ — exactly where the wedge analysis looks |
| P0.3 | One FFT over the full range returns a redshift-*averaged* spectrum | ill-defined effective redshift at $\Delta z = 1$ |
| P0.4 | PS bandwidth vs noise bandwidth | 22.28 MHz vs 8 MHz = **2.8× mismatch** at $\Delta z = 1$ |

### 11.4 `mode_counts` is computed but never used in the SNR

Binned, cached, and consumed only as an "empty bins" diagnostic. La Plante
Eqs. 15–17 divide the per-bin variance by the number of independent modes;
without it $\sigma_\times$ is overestimated and the total SNR is biased
**low**. (`TODO.md` §P1.1;
[`uncertainty_budget.md`](uncertainty_budget.md) §7.3.) The source notebook
does not apply the factor either, so the pipeline matching it here is
faithful, not a porting error.

### 11.5 The thermal noise is a scaling estimate

$T_\mathrm{sys}^2/(t_\mathrm{int}\Delta\nu)$ with no baseline density
$n(k_\perp)$ and no $X^2 Y \Omega'$ cosmological factors. The proper
La Plante Eq. 11 form already exists in
`21cm_galaxy_cross_uncertainty.ipynb` and needs moving into `src/`; 21cmSense
is the publication-grade alternative. (`TODO.md` §P1.2;
[`uncertainty_budget.md`](uncertainty_budget.md) §7.2.) The $10^3$ in
$P_N = T_\mathrm{sys}^2 \times 10^3/(t\,\Delta\nu)$ is now the named constant
`NOISE_NORMALISATION_MPC3`: it carries the Mpc³ the expression otherwise
lacks, standing in for the per-mode survey volume a full model would compute.

### 11.6 `apply_rsds=False`

21cmFAST's self-consistent RSDs are disabled and a Kaiser boost is applied
analytically afterwards. Enabling them would remove the $\beta$ approximation
— but the analytic step must then be *removed*, not skipped, or RSDs are
applied twice. (`TODO.md` §P1.3.)

### 11.7 The stored HDF5 predates the current code

`outputs/lightcone_data.h5` was written 2026-06-15 and carries 25 root
attributes — missing 9 of the 34 the current script writes.

> **Partially patched 2026-08-12.** `photoz_uncertainty` and `wedge_buffer`
> were rewritten **in place** so the analysis stage picks up the corrected
> values without a multi-hour re-simulation. A provenance string is recorded
> in the file's own `attrs_patched` attribute. `galaxy_bias` and `beta_rsd`
> were **deliberately left alone**: the stored `galaxy_overdensity` field was
> built with the $b_g = 33.39$ Kaiser boost, so those attributes still
> correctly describe the field on disk, and patching them would make the file
> internally inconsistent. Only `--sim force` fixes that.

| Attribute | Stored file | Current code |
|---|---|---|
| `galaxy_bias` | **33.3889** (still stale) | 4.744 |
| `beta_rsd` | **0.029880** (still stale) | 0.21030 |
| `wedge_buffer` | 0.02 → **patched to 0.0677** | 0.0677 |
| `photoz_uncertainty` | 0.059 → **patched to 0.45** | 0.45 |
| `DIM`, `hires_cell_size`, `M_cell_hires`, `M_cell_lores`, `sampler_min_mass` | absent | written |
| `galaxy_bias_method`, `galaxy_bias_hmf_analytic`, `t_STAR`, `sfr_timescale_yr` | absent | written |

`--sim auto` (the default) will **not** regenerate it while the file exists.
Use `bash submit_job.sh --sim force`. With $\Delta z = 0.01$ this is a cheap
re-run, and the halo catalogue is unchanged (same box, same $z_\mathrm{obs}$,
same seed 42) — only the bias-derived quantities move.

Patching the attributes updates the file's mtime, which marks
`analysis_products.h5` stale, so the next run recomputes the power spectra
(~0.6 s). The spectra do not depend on either patched attribute and come out
bit-identical (large-scale $P_\times$ mean unchanged at −5.644 × 10³).

**Post-patch pipeline output** (`--sim skip --plots none`, 2026-08-12):

| Quantity | Value |
|---|---|
| $\sigma_r$ | **157.478 Mpc** |
| Modes outside wedge | **97 / 400 = 24.25 %** |
| **Total SNR** | **1.06 × 10⁻¹¹¹ σ** — numerically zero |

The adopted $b_g$ is **unaffected** by all of this: recomputed from the stored
114,291,212-halo catalogue with current code it is **4.744191** (49,315 halos
selected, range 2.826 – 9.142), giving $\beta = f/b_g = $ **0.210293**. The
catalogue holds raw halo masses and SFRs, which no attribute patch or
post-processing parameter touches.

### 11.8 $\sigma_z$ was the fractional value used as an absolute one

**Corrected 2026-08-12: `photoz_uncertainty` 0.059 → 0.45.**

`radial_smearing_length` computes $\sigma_r = c\sigma_z/H(z)$, which requires
an **absolute** $\sigma_z$. The configured 0.059 is the *fractional* quantity
$\sigma_z/(1+z)$ that surveys actually quote — used as if it were absolute, it
understated the smearing by a factor $(1+z) = 8$ (0.059 → 0.472; the adopted
0.45 corresponds to $\sigma_z/(1+z) = 0.056$, consistent with Euclid's
requirement $\sigma_z/(1+z) < 0.05$, i.e. 0.40 at $z = 7$).

| Quantity | $\sigma_z = 0.059$ (old) | $\sigma_z = 0.45$ (current) |
|---|---|---|
| $\sigma_r$ | 20.65 Mpc | **157.48 Mpc** |
| $W$ at the smallest bin, $k_\parallel = 0.0176$ | 0.936 | **0.021** |
| $W$ at the first mode the wedge admits, 0.1118 | 0.070 | 5 × 10⁻⁶⁸ |
| $k_\parallel$ where $W = 0.5$ | 0.0574 Mpc⁻¹ | **0.0075 Mpc⁻¹** |
| Total SNR (cached spectra, current buffer) | 0.0048 σ | **0.0000 σ** |

**Consequence:** on the current geometry the cross-correlation is
identically unrecoverable — $\sigma_r = 157$ Mpc exceeds the box's own
$L_\mathrm{LOS}$ under either interpretation (§11.1), and the half-power scale
$k_\parallel = 0.0075$ Mpc⁻¹ needs $L_\mathrm{LOS} > 840$ Mpc to be sampled at
all. This is physical, not a bug: a photometric survey with $\sigma_z = 0.45$
at $z = 7$ retains essentially no line-of-sight information. Any future
forecast has to either use spectroscopic redshifts, work at
$k_\parallel \to 0$ with a box large enough to sample it, or drop to the
angular (2D) cross-correlation.

**Call sites updated:** `run_simulation.py` (config, also written to the HDF5
root attrs), `run_pipeline.py` (`data.get` fallback), `src/figures.py` (label
fallback), and the docstring of `src.analysis.radial_smearing_length`, which
now states the convention explicitly. **Not** updated: `tests/conftest.py`
keeps 0.059 deliberately (its 64 Mpc synthetic box would have a kernel that
underflows to zero in every bin, making the damping tests vacuous), and
`21cmfast_HERAxEuclid_lightcone.ipynb` still hardcodes 0.059 in its own
config cell.

**Stored HDF5 unaffected until re-run** — see §11.7: `run_pipeline.py` reads
`photoz_uncertainty` from the file's root attrs, so the new value is only a
fallback until `bash submit_job.sh --sim force`.

### 11.9 Box size

256 Mpc is why the UVLF bright end is noisy and why only ~49 k halos pass the
Euclid cut. Scaling to the Davies et al. (2025) 1 Gpc box is a storage
problem, not a code change (`BOX_LEN` and `HII_DIM` are already config
variables): volume ×59.6, halo catalogue 2.74 GB → ~163 GB, one high-res IC
array 0.23 GB → 13.5 GB of RAM. An intermediate 512 Mpc box (`HII_DIM = 256`,
8× volume, ~22 GB catalogue) would measure the scaling first. (`TODO.md` §P3.1,
[`project_update.md`](project_update.md) §12.)

---

## 12. Reference table — every number in one place

| Symbol | Value | Where set |
|---|---|---|
| `HII_DIM` / `DIM` | 128 / 384 | `run_simulation.py` §config |
| `BOX_LEN` | 256.0 Mpc | " |
| Cell size (transverse / high-res) | 2.0 / 0.6667 Mpc | derived |
| Mass resolution (lo-res / hi-res / sampler) | 3.173 × 10¹¹ / 1.175 × 10¹⁰ / 1 × 10⁸ M☉ | derived / 21cmFAST |
| $z_\mathrm{min}$, $z_\mathrm{max}$, $z_\mathrm{obs}$ | 6.995, 7.005, 7.0 | config |
| $N_z$ / `minimum_los_slices` | 100 / 100 | config |
| $L_\mathrm{LOS}$ requested / recorded | 3.4998 / 200.0 Mpc | derived / §11.1 |
| Node redshifts | 5, 7.005 → 6.995 | derived |
| Random seed | 42 | config |
| $\Omega_m$, $H_0$, $\Omega_b$ | 0.315, 67.36, 0.049 | config |
| $M_\mathrm{UV}$ window | −22 ≤ $M_\mathrm{UV}$ ≤ −18 | config / CLI |
| Equivalent SFR window | 0.7956 – 31.674 M☉ yr⁻¹ | derived |
| $\sigma_z$ (absolute) | 0.45 | config |
| $\bar n_\mathrm{gal}$ | 3 × 10⁻³ h³ Mpc⁻³ | config |
| $\kappa_\mathrm{UV}$ / AB zero point | 1.15 × 10⁻²⁸ / 51.60 | `src/conversions.py` |
| $t_\star$ / $t_H(7)$ / $t_\mathrm{sf}$ | 0.5 / 1.1406 Gyr / 570.3 Myr | template / derived |
| $b_g$ adopted / analytic / fallback | 4.744 / 5.39 / 8 | derived / config |
| Halos total / Euclid-selected | 114,291,212 / 49,315 | run output |
| $f$, $\beta$, max Kaiser boost | 0.99767, 0.21030, 1.210× | derived |
| $\sigma_r$ | 157.48 Mpc (was 20.647 at $\sigma_z = 0.059$) | derived |
| Horizon / FoV wedge slope | 3.1509 / 0.37936 | derived |
| Wedge buffer | 0.0677 Mpc⁻¹ (= 0.1 $h$ Mpc⁻¹) | config |
| $D_\mathrm{dish}$, $t_\mathrm{int}$, $\Delta\nu$ | 14 m, 1000 h, 8 MHz | config |
| $\nu_\mathrm{obs}$, $T_\mathrm{sys}$ | 177.55 MHz, 328.6 K | derived |
| $P_{N,21}$ / $P_{N,\mathrm{gal}}$ | 3.7488 mK² Mpc³ / 333.33 Mpc³ | derived |
| PS bins | 20 × 20, log-spaced | config |
| $\langle x_\mathrm{HI}\rangle$ | 0.1764 | run output |
| Total SNR | 0.0629 σ (as recorded) → 0.0048 σ (current buffer) → **0.0000 σ** (current buffer + $\sigma_z$) | §10, §11.8 |

---

## 13. User-defined parameter requirements

Everything a user has to (or may) set, **and the exact file and line where each
one is set**. §12 lists the values; this section lists the *inputs* — what is
yours to choose, what the code will not choose for you, and what silently has
no effect until the simulation is re-run.

### 13.0 The four layers

Parameters live in four places. They are **not** interchangeable, and a value
set in one layer does not propagate backwards into an earlier one.

| Layer | File | Controls | How to change it |
|---|---|---|---|
| 1. Scheduler | `submit_job.sh`, above line 25 | wall time, partition, cores, memory | edit the file (currently **empty** — see R1) |
| 2. Job wrapper | `submit_job.sh:26–33` | conda env, driver script, email, log/output paths | env var at launch, or edit |
| 3. Driver CLI | `run_pipeline.py:719–756` | staging, figures, catalogue subsampling, bright-end cut | command-line flag — no edit needed |
| 4. Simulation config | `run_simulation.py:82–170` | grid, redshifts, cosmology, survey, instrument, wedge, binning | edit the file, then re-run with `--sim force` (§13.4) |

### 13.1 Required before the first run on a new cluster

These have no value that is correct at an arbitrary site — the committed
defaults describe the developer's own setup.

| # | Requirement | Where to set it | Current value | What you must do |
|---|---|---|---|---|
| **R1** | SLURM directives | `submit_job.sh`, insert after line 24 (before `EMAIL_TO`) | **none present** | Add `#SBATCH --partition/--time/--account/--cpus-per-task/--mem`. Without them the script runs in the **foreground on whatever node invokes it** — on a login node that is a ~56 GB, multi-hour job in the wrong place. |
| **R2** | Conda environment name | `submit_job.sh:29` (`CONDA_ENV`), or `CONDA_ENV=... bash submit_job.sh` | `21cmfast` | The env must already exist and contain `py21cmfast` 4.1.1 — build it per [`INSTALL_21cmFASTv4.md`](INSTALL_21cmFASTv4.md). If `py21cmfast` is missing the run does **not** fail; it silently falls back to the synthetic lightcone (§4.3). |
| **R3** | 21cmFAST cache directory | `run_simulation.py:275`, `run_lightcone(...)` — add `cache=p21c.OutputCache("<scratch>")`. **Not currently passed**, and `~/.21cmfast/config.toml`'s `direc` does **not** control it. | `generate_lightcone`'s default is `OutputCache(direc=Path('.'))` — i.e. the **current working directory**, hence the project root (`d1f8b93.../`) | Point it at scratch, or launch from a scratch cwd. It is **~56 GB** for the current smoke test (§8) and dominates the footprint; home-quota clusters (CSD3) abort mid-run without this. Verified against py21cmfast 4.1.1: `p21c.config['direc']` (here `~/21cmFAST-cache`) governs other entry points, **not** the lightcone driver. `compute_initial_conditions` / `determine_halo_catalog` / `perturb_halo_catalog` (§4.4) take no cache argument at all. |
| **R4** | Output directory | `run_simulation.py:169` (`OUTPUT_DIR`), `run_pipeline.py` `--data/--products/--figdir/--summary` | `outputs/`, **relative to the working directory** | Must be writable and hold ~2.8 GB per run. Because it is relative, launching from a different cwd writes somewhere else. |
| **R5** | Notification address | `submit_job.sh:26` (`EMAIL_TO`) | hardcoded `sohinidutta97@gmail.com` | Change it, or the report goes to someone else. Requires `sendmail` on the compute node; if absent only the mail step fails — `EXIT_CODE` is captured before it, so the job still reports its true status. |
| **R6** | Install-time cache redirection | environment, **before** `pip install`, not in this repo | — | `PIP_CACHE_DIR`, `XDG_CACHE_HOME` → scratch, or the build dies with `Errno 122: Disk quota exceeded` (§1). |
| **R7** | Job name | `submit_job.sh:27` (`JOB_NAME`), or `JOB_NAME=... bash submit_job.sh` | `21cm_pipeline` | Only affects the log filename and the email subject. Safe to leave. |

### 13.2 Simulation parameters — `run_simulation.py`, config block

The block is delimited by the `★ CONFIGURATION — ALL USER-ADJUSTABLE
PARAMETERS ★` banner (lines 82–170). Everything below it derives from these.

**Grid**

| Parameter | Line | Value | Units | Constraint / consequence |
|---|---|---|---|---|
| `HII_DIM` | 92 | 128 | cells | Sets the transverse cell `BOX_LEN/HII_DIM`; keep it at 2 Mpc when changing `BOX_LEN`. Cost scales as `HII_DIM³`. |
| `BOX_LEN` | 93 | 256.0 | Mpc | Sets the largest mode, $\Delta k_\perp = 2\pi/\mathrm{BOX\_LEN}$. See §11.9 for the 512 Mpc / 1 Gpc scaling. |
| `DIM` | 94 | `3 × HII_DIM` = 384 | cells | The 3× ratio is the 21cmFAST convention; RAM per high-res array is `DIM³ × 4 B` (0.23 GB now, **13.5 GB** at `HII_DIM = 500`). |

**Redshift range**

| Parameter | Line | Value | Constraint / consequence |
|---|---|---|---|
| `z_min` | 112 | 6.995 | **Gated.** Widening $\Delta z$ requires `TODO.md` §P0.1–P0.4 first; the estimator assumes LOS homogeneity (§11.2). Production target is 6.5 / 7.5. |
| `z_max` | 113 | 7.005 | " |
| `minimum_los_slices` | **216** — *outside* the config block, in §1 | 100 | Floor on `N_z`. It **binds** at the current $\Delta z$ (natural $N_z = 2$), which is what produces the 57× LOS oversampling. At $\Delta z = 1$ the natural 175 wins and this becomes idle. |

`z_obs`, `L_los`, `N_z`, `lc_redshifts`, `node_redshifts` are all **derived**
(lines 200–225) — do not set them by hand.

**Euclid survey**

| Parameter | Line | Value | Units | Constraint / consequence |
|---|---|---|---|---|
| `M_UV_limit` | 118 | −18 | AB mag | Faint-end cut. Written to HDF5 and re-read by `run_pipeline.py:399` as `M_UV_faint`. |
| `M_UV_bright` | **492** — §4 of the script, *not* in the config block | −22 | AB mag | Bright-end cut for the **simulation-stage** bias. The analysis stage takes its own value from `--m-uv-bright` instead, so the two can silently disagree. Must be more negative than `M_UV_limit`. |
| `photoz_uncertainty` | 125 | 0.45 | — | **Absolute $\sigma_z$, not $\sigma_z/(1+z)$** (§11.8). Entering the survey-quoted fractional value here understates $\sigma_r$ by $(1+z)$. |
| `mean_galaxy_density` | 126 | 3 × 10⁻³ | declared `h³ Mpc⁻³` | Consumed as $P_{N,\mathrm{gal}} = 1/\bar n$ and reported in Mpc³ (§5.4). If the declared $h³$ is meant literally, the shot noise is low by $h^{-3} = 3.3\times$. Unresolved; not tracked in `TODO.md`. |
| `galaxy_bias` | 131 | 8 | — | **Fallback only.** Overwritten by the catalogue estimator whenever a halo catalogue exists (§4.6); it is used only if both the catalogue and `hmf` paths fail. |

**Cosmology and constants**

| Parameter | Line | Value | Note |
|---|---|---|---|
| `OMEGA_M_0` | 136 | 0.315 | Used for `hubble_parameter()`, `cell_mass`, the wedge, and $f = \Omega_m(z)^{0.55}$. Comoving distances use `astropy` `Planck18` instead when 21cmFAST is importable — **change both** if you move off Planck 2018. |
| `HUBBLE_CONSTANT` | 137 | 67.36 | $h = 0.6736$; the wedge buffer was converted at astropy's $h = 0.6766$ (0.4 % inconsistency, §5.3). |
| `OMEGA_B_0` | **558** — inside the analytic-bias branch | 0.049 | Only affects the `galaxy_bias_hmf_analytic` cross-check, not the adopted $b_g$. |
| `SPEED_OF_LIGHT_KMS` / `_MPS` | 142–143 | 3e5 / 3e8 | Rounded, not CODATA. |
| `F_21_MHZ` / `F_21_HZ` | 144–145 | 1420.405 / derived | Do not edit. |

**HERA instrument**

| Parameter | Line | Value | Constraint / consequence |
|---|---|---|---|
| `HERA_DISH_DIAMETER` | 150 | 14.0 m | Sets the FoV wedge slope only (drawn, not masked — §5.3). |
| `integration_time` | 151 | 1000 × 3600 s | $P_N \propto 1/t_\mathrm{int}$. |
| `bandwidth` | 152 | 8 × 10⁶ Hz | **Must match the band the power spectrum is measured over.** It does not today at $\Delta z = 1$ (2.8× mismatch, `TODO.md` §P0.4). |

**Wedge and binning**

| Parameter | Line | Value | Constraint / consequence |
|---|---|---|---|
| `wedge_buffer` | 158 | 0.0677 Mpc⁻¹ | **In Mpc⁻¹, not $h$ Mpc⁻¹.** Literature quotes 0.1 $h$ Mpc⁻¹ (Pober+2014); convert before entering. The SNR is extremely sensitive to it — 0.02 → 0.0677 costs 13× (§10). |
| `n_bins_perp` / `n_bins_parallel` | 163–164 | 20 / 20 | 145 of the 400 bins are already empty; raising these makes that worse. |
| `OUTPUT_DIR` / `OUTPUT_FILE` | 169–170 | `outputs/lightcone_data.h5` | Must agree with `run_pipeline.py --data`. |

**Hardcoded, but legitimately user-editable** (not in the config block — edit
in place, and record it, because nothing else reports the change):

| Item | Line | Value | Why you might change it |
|---|---|---|---|
| `random_seed` | 255 | 42 | Any other value gives an independent realisation — the *only* way to get cosmic-variance error bars. Changing it invalidates the 56 GB cache. |
| Template name | 253–254 | `["simple"]` | Selects the whole astrophysical source model (§4.3). A different template changes `USE_TS_FLUCT` and the runtime by a large factor. |
| `quantities` | 272 | 4 fields | Adding fields increases the HDF5 size; all four are consumed downstream. |
| `apply_rsds` | 279 | `False` | If set `True`, the analytic Kaiser step in §5 must be **deleted**, not skipped, or RSDs apply twice (`TODO.md` §P1.3). |
| `include_dvdr_in_tau21` | 278 | `False` | " |
| `N_THREADS` | **never set** | template default = 1 | Single-threaded today. Set it via `simulation_options` in the `clone()` at line 258 if you request more than one core in R1. |

### 13.3 Driver parameters — `run_pipeline.py` CLI

Full table in §3.1. The four that matter for an HPC run:

| Flag | Default | Why it matters on a cluster |
|---|---|---|
| `--sim` | `auto` | `auto` will **not** regenerate a stale HDF5 while the file exists (§11.7). Use `force` after any §13.2 edit. |
| `--max-halos` | `0` (all) | The catalogue is 2.74 GB on disk and loads in full. Set this if the node's memory is tight; densities are rescaled by `halo_sampling_factor`. |
| `--plots` | `all` | `--plots power snr` skips the catalogue read entirely; `--plots none` gives numbers only. |
| `--m-uv-bright` | −22.0 | The analysis-stage bright cut. Keep it equal to `M_UV_bright` at `run_simulation.py:492`. |
| `--galaxy-weighting` | `number` | Only affects the `euclid` figure group's rebuilt δ_gal. It does **not** change the stored field or any power spectrum. |
| `--noise-model` | `scaling` | `physical` swaps the flat thermal-noise estimate for the HERA baseline-density model (Parsons 2017 Eq. 12). ~10³ larger, and `inf` where the array has no baselines. Changes every SNR number. |
| `--mode-weighted` | off | Applies La Plante Eq. 19's `sqrt(N_patch dN)` weighting. Raises the total SNR ~10×. Changes every SNR number. |

### 13.4 Which edits actually take effect, and when

This is the most common way to get a wrong result. The analysis stage reads its
parameters **from the HDF5 root attributes written by Stage 1**, not from
`run_simulation.py`.

| Group | Parameters | Editing the config block… |
|---|---|---|
| **A — read from HDF5 attrs** (fallback used only if the attribute is absent) | `n_bins_perp`, `n_bins_parallel`, `HUBBLE_CONSTANT`, `OMEGA_M_0`, `SPEED_OF_LIGHT_KMS`, `SPEED_OF_LIGHT_MPS`, `photoz_uncertainty`, `HERA_DISH_DIAMETER`, `F_21_HZ`, `wedge_buffer`, `integration_time`, `bandwidth`, `mean_galaxy_density`, `M_UV_limit` | …has **no effect** until the file is rewritten (`--sim force`) or the attribute is patched in place. See §11.7 for the patch route and its provenance record. Fallbacks are at `run_pipeline.py:264–339`. |
| **B — baked into the stored field** | `galaxy_bias`, `beta_rsd` (the Kaiser boost is applied to `galaxy_overdensity` before writing), and all geometry (`HII_DIM`, `BOX_LEN`, `z_min`, `z_max`, `minimum_los_slices`, seed, template) | …requires a **full re-simulation**. Patching these attributes makes the file internally inconsistent — do not. |
| **C — live** | `--m-uv-bright`, `--galaxy-weighting`, `--plots`, `--format`, `--dpi`, `--max-halos`, `--data/--products/--figdir/--summary` | …takes effect on the next `run_pipeline.py` invocation, ~1.6 s from cached spectra. |

### 13.5 Run manifests and the 2026-08-20 SIGSEGV

`run_simulation.py` writes `outputs/runs/sim_<run_id>.json` before it starts
the lightcone and rewrites it after every stage (`src/provenance.py`). It is
the only record that survives a run killed by a signal, because such a run
cannot flush stdout or run an exit hook.

**What the 2026-08-20 failure looked like, and why it was undiagnosable.**
`run_pipeline.py --sim force` reported `exit code -11` (SIGSEGV) after 2,303 s
with no output at all from the child. The cause of the *silence* was
buffering: `run_pipeline.py`'s own `log()` uses `print(..., flush=True)`, while
`run_simulation.py` had no `flush=True` anywhere, so ~8 KB of the child's
progress output sat in a block buffer that the signal discarded. Both the
`python -u` in `submit_job.sh` and the `-u` in the pipeline's `subprocess.run`
call exist for this.

**What the failure was.** Commit `81e08ef` (2026-08-20 12:39) replaced the
hardcoded `HII_DIM = 128 / BOX_LEN = 256 / DIM = 384` with the
footprint-derived `256 / 486.33 / 768`. Measured against the 2026-08-12 run
and its 21cmFAST cache (136,663,818 halos in a 3.564 GiB `HaloCatalog.h5`,
i.e. 28.0 bytes/halo):

| | 256 Mpc (2026-08-12, OK) | 486.33 Mpc (2026-08-20, SIGSEGV) |
|---|---|---|
| Comoving volume | 1.68 × 10⁷ Mpc³ | 1.15 × 10⁸ Mpc³ (6.86×) |
| Lagrangian halos | 136,663,818 | ~9.37 × 10⁸ |
| Catalogue on disk | 3.83 GB | 26.2 GB |
| Both catalogues resident | 7.0 GB | 48–52 GB |
| IC file (`DIM` 384 → 768) | 0.96 GB | 7.7 GB |
| `halo_coords` vs `INT_MAX` | 0.19× | **1.31×** |

Two mechanisms are consistent with SIGSEGV rather than SIGKILL: an unchecked
`malloc` returning NULL in the C backend, or a 32-bit index overflow. The
evidence favours the latter. Of the 16 cached `HaloCatalog.h5` files in the
working tree, exactly one is unreadable — and it is the one written at
`BOX_LEN = 486.33` with `HII_DIM = 32, DIM = 96`, a grid 512× smaller and
under no memory pressure at all. Its size is 2,147,491,839 bytes: the signed
32-bit boundary. **If it is the overflow, a larger node will not fix it** —
the catalogue has to shrink.

Hence `estimate_catalogue_cost()` and the pre-flight warning: the estimate is
volume-scaled from the measured run and is printed, and recorded in the
manifest, before any compute is spent.

**Mitigations now in place**

| Change | Where | Effect |
|---|---|---|
| `python -u` on the child | `submit_job.sh`, `run_pipeline.py:199` | The next failure names its stage |
| Run manifest | `src/provenance.py`, `run_simulation.py` | Parameters and stage-reached survive a crash |
| `del halo_catalog, initial_conditions` | `run_simulation.py`, after `perturb_halo_catalog` | Drops ~26 GB of catalogue and ~7.7 GB of ICs at the high-water mark |
| `MINIMIZE_MEMORY = True` | `matter_options` | Trades peak RAM for intermediate I/O |
| `N_THREADS` resolution | `provenance.resolve_n_threads()` | 21cmFAST's default is 1; the failed run used one core for 38 minutes |
| Pre-flight cost estimate | `estimate_catalogue_cost()` | Warns past `INT_MAX` before any compute is spent |

None of these makes the 486.33 Mpc box fit. They make the next attempt
diagnosable and cheaper, and they say up front when a box cannot work.

### 13.6 Consistency rules the code does not enforce

Nothing checks these; violating one produces a plausible-looking wrong number.

1. `BOX_LEN / HII_DIM` should stay at the 2 Mpc cell — the mass resolution
   tables in §4.2 and the `minimum_los_slices` logic both assume it.
2. `DIM = 3 × HII_DIM` — the 21cmFAST convention; break it deliberately or not
   at all.
3. `M_UV_bright < M_UV_limit` (more negative), and `run_simulation.py:492` must
   match `--m-uv-bright`.
4. `photoz_uncertainty` is **absolute**; `wedge_buffer` is in **Mpc⁻¹**.
5. `bandwidth` must equal the frequency span the power spectrum is measured
   over (`TODO.md` §P0.4).
6. Widening `z_min`/`z_max` is gated on `TODO.md` §P0.1–P0.4 (§11.2).
7. `OUTPUT_FILE` must equal `run_pipeline.py --data`.
8. Requesting more than one core in R1 does nothing unless `N_THREADS` is also
   set (§13.2).

### 13.7 Pre-flight checklist

```text
[ ] R1  #SBATCH directives added to submit_job.sh
[ ] R2  conda env `21cmfast` exists; `python -c "import py21cmfast"` succeeds
[ ] R3  21cmFAST cache pointed at scratch, >= 60 GB free (>= 200 GB for 1 Gpc)
[ ] R4  outputs/ writable, >= 3 GB free; launching from the repo root
[ ] R5  EMAIL_TO updated; sendmail present (or accept a failed mail step)
[ ] --- if any run_simulation.py value changed ---
[ ] Group A/B edit  =>  launched with `--sim force`, not the default `auto`
[ ] M_UV_bright (line 492) == --m-uv-bright
[ ] seed / template change recorded in CHANGELOG.md
[ ] --- if BOX_LEN changed (§13.5) ---
[ ] Pre-flight estimate printed no INT_MAX warning
[ ] Node RAM >= the manifest's cost_estimate.resident_GB, plus the IC file
[ ] N_THREADS set (its own default is 1)
```

After any failed run, read `outputs/runs/sim_<run_id>.json` first: `status`
and `stage` name where it died even when the log does not.

---

## 14. References

Davies, J. et al. (2025), arXiv:2504.17254 — 21cmFASTv4 discrete source model ·
Gagnon-Hartman, Davies & Mesinger (2025), arXiv:2502.20447 ·
La Plante, P. et al. (2023), arXiv:2205.09770 — wedge geometry, cross-spectrum
variance (Eqs. 10, 11, 15–17) ·
Pober, J. et al. (2014), arXiv:1310.7031 — "moderate" foreground model,
0.1 $h$ Mpc⁻¹ buffer · Parsons, A. et al. (2012a), arXiv:1204.4749 ·
Thyagarajan, N. et al. (2015), ApJ 804, 14 ·
Park, J. et al. (2019), MNRAS 484, 933 — $t_\star t_H$ SFR prescription ·
Sheth & Tormen (1999), MNRAS 308, 119 — halo bias, $\nu$ convention ·
Madau & Dickinson (2014), ARA&A 52, 415 — UV–SFR calibration ·
Oke & Gunn (1983) — AB magnitudes ·
Bouwens et al. (2021), AJ 162, 47 · Finkelstein et al. (2015), ApJ 810, 71 ·
Song et al. (2016), ApJ 825, 5 · González et al. (2010), ApJ 713, 115 ·
Speagle et al. (2014), ApJS 214, 15 · Schreiber et al. (2015), A&A 575, A74 ·
Euclid Collaboration (2022), arXiv:2108.01201.
