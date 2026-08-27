# Reference

Detailed reference material for the 21 cm × galaxy cross-correlation project.
This is the long-form companion to [`README.md`](../README.md), which carries
the overview and the quickstart; everything here is the detail that used to
live inline there.

**Contents**

- [Notebooks](#notebooks) — structure, equations, and fiducial parameters of each notebook and of the HPC pipeline
- [Figure literature references](#notebooksplot_fieldsipynb--figures-and-literature-references) — the published relations overlaid on each figure
- [21cmFASTv4 `HaloBox` API notes](#21cmfastv4-halobox-api-notes) — field access, the SFR unit warning, the lightcone API, and the source-model templates
- [Source modules](#source-modules) — the public surface of `src/conversions.py`, `src/analysis.py`, and `src/FOV_to_cMpc.py`
- [Requirements](#requirements) — which package is needed by which entry point
- [Test suite coverage](#test-suite-coverage) — what each `tests/test_*.py` file exercises
- [Figure display in notebooks](#figure-display-in-notebooks) — `%matplotlib widget` and constrained layout
- [Bibliography](#bibliography)

For the parameter-level ground truth of an HPC run see [`HPC.md`](HPC.md); for
the photo-$z$ / wedge / noise / SNR chain see
[`uncertainty_budget.md`](uncertainty_budget.md).

---

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

Lightcone counterpart to the archived coeval notebook (2). Uses `RectilinearLightconer` + `run_lightcone` to produce a self-consistent lightcone over a redshift range (derived from the survey footprint: $z = 6.55$–$7.45$), then applies the same galaxy field construction, Kaiser RSD, 2D power spectrum calculation, and SNR estimation as the coeval notebook.

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
- **3c** After the Euclid cut — **3c.1** the galaxies, halo masses and SFRs that
  survive the $M_\mathrm{UV}$ window (`select_euclid_halos`), with the halo-mass
  and SFR distributions shown before and after the cut and the equivalent SFR
  window marked; **3c.2** $\delta_\mathrm{gal}$ rebuilt from the selected
  catalogue alone with `galaxy_overdensity_from_catalogue()` on
  `run_simulation.py` §3b's grid (`n_perp = HII_DIM`, `n_los = N_z`,
  `los_extent = BOX_LEN`), shown as an LOS projection, a single transverse
  slice, and a symlog one-point distribution. Weighting switchable via
  `GALAXY_WEIGHTING_DIAGNOSTIC`. **This is not §3a's field**, which comes from
  the lightcone `halo_sfr` and applies no magnitude cut
- **4** Kaiser RSD applied in Fourier space
- **5** Visualise lightcone — transverse (x–y) slice + LOS (x–z) slice with
  **LOS on the x-axis and transverse on the y-axis** (standard convention);
  secondary redshift axis on top via `twiny()`
- **5b** Brightness temperature evolution plot — wide-format (16×3.5") lightcone
  slice styled after Mesinger & Furlanetto (2007), with a custom EoR colourmap
  (dark = ionised, warm/bright = neutral), dual x-axes (comoving distance +
  redshift), and observed frequency range in the title
- **5c** §3c.2's post-cut $\delta_\mathrm{gal}$ overlaid on the 21 cm field —
  $\delta_\mathrm{gal}$ contours over $\delta T_b$, the same maps with the roles
  swapped, and $\langle\delta_\mathrm{gal}\rangle$ binned by $\delta T_b$ with
  the cell-by-cell Pearson $r$. Both maps are the mean over the *same* window of
  LOS cells, so the pairing matches `compute_all_power_spectra`'s. Transverse
  only: the catalogue is coeval and the field is a lightcone (see the
  line-of-sight caveat under `src/analysis.py` below). Placed in §5 rather than
  §3 because it uses §5b's `eor_cmap`
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
- **7e** Foreground contamination and its removal — a synthetic diffuse
  Galactic synchrotron foreground (plus point sources) injected into
  `brightness_temp_field` before any spectrum is computed, then removed at
  0/50/90/99/99.9/100 %. Each level goes through `compute_all_power_spectra`
  and `compute_uncertainty_budget` unchanged. Three panels: $P_{21}$ showing
  where the foreground sits and the removal working, $|P_\times|$ against
  $\sigma_\times$, and total SNR versus removal fraction — plotted both as
  measured and signal-only, since the as-measured curve has a contaminated
  numerator. See `src/foregrounds.py` above; the removal step is a
  placeholder, not an algorithm
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
| $M_\mathrm{UV}$ window | $-22.66 < M_\mathrm{UV} < -18.0$ | Euclid galaxy selection — `M_UV_bright` / `M_UV_faint`, both from collaborators. `run_simulation.py` uses a bright end of $-22$ instead |
| $\sigma_z$ | 0.45 | Euclid photo-$z$ uncertainty, **absolute — not $\sigma_z/(1+z)$**. Euclid's $\sigma_z/(1+z) < 0.05$ requirement is $\sigma_z \approx 0.45$ at $z = 7$; the earlier 0.059 was the fractional number used as though absolute |
| $\bar{n}_\mathrm{gal}$ | $3\times10^{-3}\ h^3\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density |
| $b_\mathrm{gal}$ | **computed in-line**, $\approx 5.39$ | No fixed fallback — `galaxy_bias = 8` is commented out in the configuration cell. The value is calculated in **§3**, the analytic-galaxy-bias cell (`galaxy_bias = galaxy_bias_hmf`): a Simpson integral of the Sheth-Tormen halo bias weighted by the HMF over the Euclid-selected mass bins, $b_g = \int b_h\,\frac{dn}{d\log M}\,d\log M \big/ \int \frac{dn}{d\log M}\,d\log M$. Consumed downstream by §4's $\beta = f/b$. Requires `hmf` |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| Bandwidth | 8 MHz | HERA per-band bandwidth |

---

### 4. HPC pipeline — `run_simulation.py` + `notebooks/plot_fields.ipynb` + `notebooks/analysis.ipynb`

A refactored version of notebook 3 split into three independent parts for cluster use. See the [Quickstart](../README.md#quickstart) in the README for how to run it.

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
| $\sigma_z$ | **0.256** | Euclid photo-$z$ uncertainty — **absolute**, not $\sigma_z/(1+z)$. Euclid Collab.: Allen et al. (2026), A&A 711, A25, Sect. 3 / Fig. 4: measured $\sigma_{\rm nmad} \le 0.032$, converted as $0.032\times(1+z)$. Supersedes 0.45, which was the *pre-launch requirement* $\sigma_z/(1+z) < 0.05$ rather than a measurement. Two caveats — worst-field bound, and the NMAD normalisation assumption — in `NUMBERS_AND_SOURCES.md` §5 |
| $\bar{n}_\mathrm{gal}$ | $7.48\times10^{-5}\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density — Euclid Collab.: Allen et al. (2026), A&A 711, A25, Table 2, $N/V$. **Plain Mpc⁻³, not $h^3$ Mpc⁻³.** Derived under the paper's own cosmology, not Planck18; see `NUMBERS_AND_SOURCES.md` §2 |
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
> [`TODO.md`](../TODO.md) §P0 first — otherwise the FFT is applied to unevenly
> sampled ($20.4\%$ cell-size spread), redshift-evolving data. Treat the
> resulting SNR as a smoke-test number, not a forecast.

> **The stored `outputs/lightcone_data.h5` is still stale.** The galaxy-bias
> calculation has been corrected (see
> [`project_update.md`](project_update.md) §4), which changes
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
| SFR → $L_\mathrm{UV}$ | $L_\mathrm{UV} = \mathrm{SFR} / \kappa_\mathrm{UV}$, $\kappa_\mathrm{UV} = (2.7 \pm 0.9)\times10^{-29}\ M_\odot\ \mathrm{yr}^{-1}\ /\ (\mathrm{erg\ s}^{-1}\ \mathrm{Hz}^{-1})$ | Fisher et al. (2026), arXiv:2511.10741, Eq. 12 — rising-SFH calibration recovering $\mathrm{SFR}_{100\,\mathrm{Myr}}$, **not** a constant/instantaneous SFR. ~33 % systematic. Three caveats in `NUMBERS_AND_SOURCES.md` §2 |
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

**UV luminosity – SFR conversions** (Fisher et al. 2026, arXiv:2511.10741, Eq. 12, $\kappa_\mathrm{UV} = (2.7 \pm 0.9)\times10^{-29}\ M_\odot\ \mathrm{yr}^{-1}\ /\ (\mathrm{erg\ s}^{-1}\ \mathrm{Hz}^{-1})$)

> **These recover $\mathrm{SFR}_{100\,\mathrm{Myr}}$ from rising SFHs**, not an
> instantaneous or constant SFR. Substituting this value changes what "SFR"
> *means* downstream, not just its magnitude — see `NUMBERS_AND_SOURCES.md`
> §2, Caveat 1.

| Function | Description |
|----------|-------------|
| `Luv_to_sfr(Luv, kappa_uv=2.7e-29)` | UV luminosity [erg s⁻¹ Hz⁻¹] → SFR [M☉ yr⁻¹] |
| `sfr_to_Luv(sfr, kappa_uv=2.7e-29)` | SFR [M☉ yr⁻¹] → UV luminosity [erg s⁻¹ Hz⁻¹] |
| `sfr_to_Muv(sfr, kappa_uv=2.7e-29)` | SFR [M☉ yr⁻¹] → absolute UV AB magnitude (chains `sfr_to_Luv` → `Luv_to_Muv`) |

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
| Halo catalogue (`SAMPLER_MIN_MASS`) | — | $2\times10^{8}\ M_\odot$ (smallest sampled halo; raised from $1\times10^{8}$ — see `HPC.md` §11.13) |

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

**The `euclid` figure group rebuilds the field regardless of that setting.**
`figures.selected_galaxy_overdensity()` calls
`galaxy_overdensity_from_catalogue()` itself, at `n_perp = HII_DIM`,
`n_los = N_z`, `los_extent = BOX_LEN`, with the per-halo weight taken from
`--galaxy-weighting`. It has to: under the default `"lightcone_sfr"` the stored
field carries no magnitude cut at all, so plotting it would not show anything
"after the Euclid cut". The rebuild is a figure-stage diagnostic — it never
reaches the power spectra or the stored HDF5.

> **Line-of-sight caveat:** the halo catalogue is *coeval* at `z_obs` and spans
> `BOX_LEN` along the LOS, whereas the lightcone spans `L_los`. The catalogue
> modes are binned into an `(HII_DIM, HII_DIM, N_z)` grid so the shape matches
> downstream, but they carry no redshift evolution along the LOS and their LOS
> cell size is `BOX_LEN / N_z`.

### `src/foregrounds.py` — foreground injection and parametrised removal

Adds a synthetic foreground to the simulated brightness-temperature lightcone
before any power spectrum is computed, and removes a controllable fraction of
it, so the notebook's §7e can measure the cost of contamination and of
incomplete removal. Consumes the verified analysis chain; modifies none of it.

| Function / class | Description |
|---|---|
| `simulate_diffuse_foreground(shape, k_perp, k_parallel, z_obs, **kw)` | Diffuse Galactic synchrotron cube [mK] on the lightcone grid |
| `simulate_point_source_foreground(...)` | Poisson point-source component, angularly flatter than the diffuse one |
| `inject_foreground(field, k_perp, k_parallel, z_obs, foreground_amplitude, ...)` | Combines both, scales to `foreground_amplitude × signal RMS`, returns a `ForegroundRealisation` |
| `remove_foreground(contaminated, k_perp, k_parallel, removal_fraction, *, foreground, removal_basis)` | Placeholder removal knob — see the warning below |
| `DIFFUSE_DEFAULTS`, `POINT_SOURCE_DEFAULTS` | Model parameters with their sources |

**The model.** A power-law angular power spectrum $C_\ell \propto \ell^{-\beta}$
with a smooth power-law spectrum
$T(\theta,\nu) = T_{\rm ref}(\theta)(\nu/\nu_{\rm ref})^{-\alpha(\theta)}$,
$\alpha$ varying across the sky. The reference sky is log-normal so the
temperature is positive everywhere. The line-of-sight structure is *emergent*
— it follows from the smooth frequency dependence rather than from an imposed
$k_\parallel$ power law.

| Parameter | Diffuse | Point source | Source |
|---|---|---|---|
| $\beta$ (angular) | 2.4 | 1.1 | Santos, Cooray & Knox (2005) Table 1 |
| $\alpha$ (frequency) | 2.8 | 2.07 | Santos+2005; GSM |
| $\sigma_\alpha$ | 0.1 | 0.3 | Shaw et al. (2014) |
| $T_{\rm ref}$ at 130 MHz | 700 K | 57 K | de Oliveira-Costa+2008; Zheng+2017 (GSM2016) |

> **`remove_foreground` is a placeholder, not a method.** It subtracts an
> exactly-correct template of the very field that was injected, scaled by
> `removal_fraction`. It is **not** GMCA, PCA, ICA, polynomial or
> log-polynomial fitting, Gaussian-process removal, or delay filtering, and it
> has none of their failure modes: no signal loss from over-fitting, no
> mode-mixing, no leakage from the wedge into the window. Results are
> statements about *removal level*, not about any named method's achievable
> performance. `removal_basis` selects whether `removal_fraction` is a
> fraction of the foreground amplitude (default; residual power falls as
> $(1-f)^2$) or of its power.

**Two findings the module documents and tests.**

*Smoothness does not confine the contamination to low $k_\parallel$.*
`compute_cylindrical_cross_power` takes a bare FFT with no line-of-sight
taper, so a smooth-but-non-periodic spectrum is discontinuous at the box edge
and leaks along the whole $k_\parallel$ axis as $\approx k_\parallel^{-1.5}$.
That slope comes from the window, not the sky — a bare ramp with no angular
structure leaks identically, and widening the band does not change it. Wedge
excision alone therefore does not remove foreground power from the EoR window
here.

*A contaminated SNR flatters the result.* Foregrounds are unbiased in the
ensemble mean, but a single realisation carries a chance cross-correlation of
order $\sqrt{P_{21}P_{\rm gal}/N_{\rm modes}}$. Measured on the synthetic
fixture, the per-bin shift in $P_\times$ scales **linearly** with foreground
amplitude (×10 per decade) while the shift in $P_{21}$ scales
**quadratically** (×100) — exactly the contrast expected if the foreground
reaches the cross-spectrum only through a chance correlation. Because
$|P_\times|/\sigma_\times$ then has a contaminated numerator *and*
denominator, the total SNR degrades far more slowly than $\sigma_\times$
alone: in one test $\sigma_\times$ rose 583× while the total SNR fell only
6.8×. Evaluate the SNR with the clean $P_\times$ against the contaminated
$\sigma_\times$ to separate the two; §7e.3 plots both.

**Dependencies:** `numpy` only.

### `src/provenance.py` — run manifests and pre-flight costing

Every simulation run writes one JSON manifest, `outputs/runs/sim_<run_id>.json`,
recording what it was configured to do and how far it got.

The manifest is rewritten on **every** update, not once at the end. That is the
whole design constraint: a process killed by a signal cannot flush stdout or
run an exit hook, so anything written only at exit is lost. A crashed run
leaves the manifest with `"status": "running"` and `"stage"` naming where it
died. Writes go through a temporary file and `os.replace`, so a crash
mid-write cannot leave half-parsed JSON.

| Function / class | Description |
|---|---|
| `RunManifest.create(output_dir, label="sim", repo_root=None, run_id=None)` | Start a manifest; written to disk immediately with `status = "running"` |
| `RunManifest.record(section, values)` | Merge into a top-level section and rewrite |
| `RunManifest.begin_stage(name)` / `.end_stage()` | Mark a stage; `end_stage` logs its duration and returns it |
| `RunManifest.finish(status)` | Close with `"complete"` or `"failed"` |
| `resolve_n_threads(default=None)` | `N_THREADS` env → `SLURM_CPUS_PER_TASK` → `default` → `os.cpu_count()` → 1. Non-numeric values are ignored, never raised |
| `estimate_catalogue_cost(box_len)` | Pre-flight halo-catalogue estimate (see below) |
| `environment_info(repo_root)` | Host, platform, interpreter, SLURM ids, git revision, package versions |
| `git_revision(repo_root)` | `commit` / `branch` / `dirty`; all `None` outside a repository |
| `peak_memory_gb()` | Peak RSS, handling `ru_maxrss` being KB on Linux and bytes on macOS |

**Manifest sections:** `parameters`, `derived`, `cost_estimate`, `environment`,
`timings_seconds`, `peak_memory_GB`, `results`, `outputs`, plus `status`,
`stage` and `stages_completed`.

**Pre-flight costing.** `estimate_catalogue_cost` extrapolates from a measured
run rather than a model: 136,663,818 halos in a (256 Mpc)³ box at $z = 7$ with
`SAMPLER_MIN_MASS = 1e8`, stored in a 3.564 GiB `HaloCatalog.h5` — 8.146 halos
Mpc⁻³ at 28.0 bytes each. The sampler's floor is a fixed *mass*, not a grid
property, so the catalogue scales with comoving volume and is independent of
`HII_DIM`.

It returns `n_halos_lagrangian` (what `determine_halo_catalog` draws — the
count that sets peak memory and index width), `n_halos_perturbed` (the 83.6 %
that survives `perturb_halo_catalog` and reaches the HDF5, so this is the one
to compare against a run's `results.n_halos`), `catalogue_GB`, `resident_GB`
(both catalogues, as held simultaneously during perturbation), and
`int32_headroom` — the flattened `halo_coords` length as a fraction of
`INT_MAX`. The 83.6 % ratio comes from the 256 Mpc run and was reproduced
independently by a 64 Mpc run (1,782,540 perturbed halos against a predicted
1,785,000).

A headroom above 1.0 means the coordinate array is longer than a signed 32-bit
index can address, and `run_simulation.py` says so before spending any
compute. See [`HPC.md`](HPC.md) §13.5 for why that matters and what it is
calibrated against.

**Dependencies:** standard library only.

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
| *(stdlib only)* | `src/provenance.py` |
| `numpy` only | `src/foregrounds.py` |
| `ipympl` | All notebooks — backs the `%matplotlib widget` interactive inline backend |

The analytical notebook (`21cm_galaxy_cross_uncertainty.ipynb`) requires only `numpy` and `matplotlib`; all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).


---

## Test suite coverage

`conda run -n 21cmfast pytest tests/ -v` — 146 tests, all passing in ~22 s. No
test invokes 21cmFAST: `tests/conftest.py` writes a synthetic
`lightcone_data.h5` with the same schema (16² × 12 cells, 4 000 halos), so the
whole suite runs offline.

| File | Covers |
|------|--------|
| `test_analysis.py` | Cosmology helpers, the cylindrical power-spectrum estimator (shapes, symmetry, sign, and the analytic white-noise normalisation $P = \sigma^2 V_\mathrm{cell}$), wedge geometry, photo-$z$ kernel, thermal noise, SNR, Euclid selection, effective bias |
| `test_dataio.py` | HDF5 loading, metadata accessors, halo subsampling, field/catalogue skipping, product-cache round trip and staleness detection |
| `test_figures.py` | Headless backend, NaN filling, colormap, and that every one of the figure functions renders and writes a non-empty file |
| `test_foregrounds.py` | Foreground shape/positivity, spectral smoothness and the leakage that survives it, the angular and frequency power laws, removal exactness at 0 % and 100 %, both removal bases, and end-to-end SNR degradation and recovery through the unmodified analysis chain |
| `test_provenance.py` | Thread resolution from the environment (including malformed values), git/package capture, the volume-scaled cost estimate against both the reference run and the box that crashed, and that a manifest left mid-stage records the stage it died in |
| `test_pipeline.py` | CLI parsing, each stage's `auto`/`force`/`skip` behaviour (with a stub simulation script), and end-to-end runs checking figures, cache reuse, and the summary JSON |
| `test_galaxy_weighting.py` | The number- and luminosity-weighted `delta_gal` constructions: weight conservation, non-cubic grids, zero-mean normalisation, both formulas against manual recomputation, $\kappa_\mathrm{UV}$ scale-invariance, and the error paths |
| `test_conversions.py` | Mass-resolution helpers: $\bar\rho_m$ against astropy, $H_0^2$ and $L^3$ scalings, the production-grid values, and total-box mass conservation |

> `src/FOV_to_cMpc.py` and the magnitude/SFR half of `src/conversions.py` are
> still untested directly, though the conversion round-trips are exercised
> indirectly through the selection and UVLF tests. The remaining first
> candidates are the explicit identities `Muv_to_Luv` ↔ `Luv_to_Muv` and
> `sfr_to_Luv` ↔ `Luv_to_sfr`.

---

## Figure display in notebooks


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


---

## Bibliography


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
- **Fisher et al. (2026)**, MNRAS in press — [arXiv:2511.10741](https://arxiv.org/abs/2511.10741) — *REBELS-IFU*: rising SFHs at z = 7; **the adopted $\kappa_\mathrm{UV}$** (Eq. 12)
- **Madau & Dickinson (2014)**, ARA&A, 52, 415 — [arXiv:1403.0007](https://arxiv.org/abs/1403.0007) — UV luminosity–SFR calibration; **superseded** source of the previous $\kappa_\mathrm{UV} = 1.15\times10^{-28}$
- **Dhandha et al. (2026)**, MNRAS in press — [arXiv:2508.13761](https://arxiv.org/abs/2508.13761) — Eq. 18 independently confirms $1.15\times10^{-28}$; corroborates the **superseded** value
- **Sheth & Tormen (1999)**, MNRAS, 308, 119 — [arXiv:astro-ph/9901122](https://arxiv.org/abs/astro-ph/9901122) — halo mass function and bias formula
- **Murray, Robotham & Power (2013)**, Astron. Comput., 3, 23 — [arXiv:1306.6721](https://arxiv.org/abs/1306.6721) — `hmf` halo mass function code
