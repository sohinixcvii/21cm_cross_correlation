# 21cm × Galaxy Cross-Correlation — Project Status Update

**Date:** 2026-08-04
**Supersedes:** the 2026-06-15 update (see §11 for what changed and why)

---

> **Read this first.** The galaxy-bias calculation in `run_simulation.py` has
> been corrected, so `outputs/lightcone_data.h5` from 2026-06-15 is **stale**:
> its `galaxy_bias` attribute and its `galaxy_overdensity` field (which
> carries the Kaiser boost) are both wrong. The numbers in §2–§8 below come
> from that superseded run and are retained for comparison. A fresh run is
> required:
>
> ```bash
> bash submit_job.sh --sim force
> ```
>
> `--sim auto` (the default) will **not** re-run it — the file still exists.
> The redshift range is unchanged, so this is a cheap re-run.

---

## 1. Simulation Parameters

| Parameter | Value | Note |
|---|---|---|
| Code | py21cmfast 4.1.1 ("simple" template) | |
| Grid size | HII_DIM = 128, DIM = 384 | |
| Box side length | BOX_LEN = 256.0 comoving Mpc | see §12 for the 1 Gpc target |
| Cell size | 2.0 Mpc | transverse |
| Redshift range | z = 6.995 → 7.005 | Δz = 0.01 — **deliberate smoke test** |
| Reference redshift | z_obs = 7.0 | lightcone midpoint |
| LOS slices | N_z = 100 | floored by `minimum_los_slices` |
| LOS extent | ~3.5 Mpc physical (200 Mpc reported) | |
| LOS cell size | 0.035 Mpc | ~57× oversampled vs transverse |
| Node redshifts | 5 | |
| Random seed | 42 | |

### Why the range is still narrow

A production range (z = 6.5 → 7.5, L_LOS = 350.8 Mpc, N_z = 175) was set on
2026-08-04 and then **deliberately reverted**. The power-spectrum estimator in
`src/analysis.py` was inherited from the coeval notebook and assumes
statistical homogeneity along the line of sight. That holds for a
quasi-coeval slab, so at Δz = 0.01 the configuration and the formalism match.
At Δz = 1.0 it fails in four measurable ways — non-uniform comoving sampling
(20.4 % cell-size spread), a spurious LOS ramp from global mean subtraction, a
redshift-averaged spectrum with no well-defined effective redshift, and a
2.8× mismatch between the spectrum bandwidth (22.28 MHz) and the noise model
(8 MHz).

Widening Δz is therefore **gated on the estimator work in `TODO.md` §P0**,
where all four are documented with numbers and fixes. Until then, treat the
SNR as a smoke-test diagnostic rather than a forecast.

### Source model parameters ("simple" template defaults)

| Parameter | Value | Description |
|---|---|---|
| F_STAR10 | 0.05 (log₁₀ = −1.3) | Star formation fraction at 10¹⁰ M☉ |
| ALPHA_STAR | 0.5 | Power-law slope of SHMR |
| t_STAR | 0.5 | Fraction of Hubble time for SFR timescale |
| M_TURN | 5.01 × 10⁸ M☉ (log₁₀ = 8.7) | Turnover mass (exponential suppression) |
| F_ESC10 | 0.1 (log₁₀ = −1.0) | Escape fraction at 10¹⁰ M☉ |
| ALPHA_ESC | −0.5 | Escape fraction slope |
| SIGMA_STAR | 0.25 dex | Log-normal scatter in stellar mass |
| SIGMA_SFR_LIM | 0.19 dex | Floor scatter in SFR |

### Cosmology (Planck 2018)

| Parameter | Value |
|---|---|
| Ω_m | 0.315 |
| Ω_Λ | 0.685 |
| H₀ | 67.36 km/s/Mpc |
| Ω_b | 0.049 |

---

## 2. Halo Catalogue *(superseded run)*

Extracted at z_obs = 7.0 using `determine_halo_catalog` + `perturb_halo_catalog`.

| Quantity | Value |
|---|---|
| Total halos | 114,291,212 |
| Halo mass range | 1.0 × 10⁸ – 1.77 × 10¹² M☉ |
| Stellar mass range | 29 – 2.79 × 10¹¹ M☉ |
| SFR range | 2.2 × 10⁻¹⁰ – 488 M☉ yr⁻¹ |
| SFR median | 1.17 × 10⁻⁵ M☉ yr⁻¹ |

These numbers carry over unchanged to the re-run: the catalogue is drawn from
the same 256 Mpc box, at the same z_obs, with the same random seed. Only the
bias-derived quantities change (§4).

---

## 3. Euclid Survey Cuts

| Parameter | Value |
|---|---|
| UV magnitude window | −22 ≤ M_UV ≤ −18 (AB mag) |
| Photo-z uncertainty | **σ_z = 0.45** (absolute, not σ_z/(1+z); was 0.059 — corrected 2026-08-12, see `docs/HPC.md` §11.8) |
| Target galaxy number density | n̄ = 3 × 10⁻³ h³ Mpc⁻³ |
| Halos in the window *(superseded run)* | **49,315** |
| Halos brighter than −18, no bright cut | 49,621 |
| Fraction Euclid-detected | ~0.04 % of all halos |

The 306-halo difference between the two counts is the population brighter than
M_UV = −22, which the UVLF figure includes but the bias selection excludes.

The selection uses the Madau & Dickinson (2014) calibration throughout, via
`src/conversions.py`:

```
L_UV [erg/s/Hz] = SFR [M☉/yr] / κ_UV,   κ_UV = 1.15 × 10⁻²⁸
M_UV = −2.5 log₁₀(L_UV) + 51.60
```

---

## 4. Galaxy Bias — **corrected**

### 4.1 What the previous update got wrong

The 2026-06-15 update recorded `b_g = 33.39` and attributed it to "the Euclid
bright limit being too restrictive" plus "a manual SFR timescale (100 Myr)".
**That diagnosis was incorrect.** Both effects are real but minor. The
dominant cause was a ν-convention error in `run_simulation.py`'s local bias
helper:

```python
# The original helper — note that nu is squared AGAIN
def sheth_tormen_bias_from_nu(nu, delta_c=1.686, a=0.707, p=0.3):
    return 1.0 + (a * nu**2 - 1.0) / delta_c + ...
```

`hmf`'s `MassFunction.nu` is **already** the squared peak height (δ_c/σ)².
Squaring it a second time inflates ν from the range 2.37–51.7 to 5.6–2670,
and the bias with it. This is the same convention error that was diagnosed
and fixed in `notebooks/analysis.ipynb` ("Fix 1") but never back-ported to
`run_simulation.py`.

Re-running the original code path reproduces the stored value exactly:

| Bias helper | Result |
|---|---|
| `mf.nu` squared again (original bug) | **33.39** ← matches the stored HDF5 attribute |
| `mf.nu` used correctly | 4.23 |

### 4.2 Current implementation

`run_simulation.py` now uses `src.conversions.sheth_tormen_bias`, which takes
the squared peak height directly, and computes two estimates:

| Estimator | b_g | Notes |
|---|---|---|
| **Halo catalogue (adopted)** | **4.744** | Per-halo M_UV from the catalogue's own SFR, then the mean Sheth-Tormen bias over the survivors. Inherits 21cmFAST's log-normal scatter (σ_star = 0.25 dex, σ_SFR = 0.19 dex). Range 2.83–9.14. |
| Analytic HMF integral (cross-check) | 5.39 | Sheth-Tormen integrated over the HMF weighted by the *mean, scatter-free* scaling relation. |

Both are stored: `galaxy_bias` (adopted), `galaxy_bias_hmf_analytic`, and
`galaxy_bias_method` as HDF5 attributes.

The catalogue estimator is preferred because the scatter matters: abundant
low-mass halos scatter up into the magnitude window, and they carry lower
bias. The scatter-free analytic version sees only rare high-mass halos, so it
runs ~14 % high.

### 4.3 SFR timescale

The 100 Myr timescale is also fixed. `src.analysis.star_formation_timescale`
returns 21cmFAST's own value:

```
t_sf = t_STAR × t_H(z) = 0.5 × 1.141 Gyr = 570.3 Myr   (at z = 7)
```

Using 100 Myr overestimates SFR by 5.70×, i.e. makes every galaxy 1.89 mag too
bright. On the analytic path this shifts b_g from 4.23 to 5.39.

### 4.4 Downstream effect on RSD

| Quantity | Superseded | Current |
|---|---|---|
| Galaxy bias b_g | 33.39 | 4.744 |
| Growth rate f = Ω_m(z)^0.55 | 0.9977 | 0.9977 |
| β = f / b_g | 0.0299 | **0.2103** |
| Max Kaiser boost (μ = 1) | 1.030× | **1.210×** |

The Kaiser boost is therefore ~7× stronger than the stored run applied — one
more reason the existing HDF5 must be regenerated.

---

## 5. Predicted Observable Plots

All figures are now produced programmatically by `src/figures.py` and written
to `outputs/figures/` by `run_pipeline.py`; the notebooks render the same
content inline.

### 5.1 UV Luminosity Function

- **Data**: all halos with SFR > 0 (no magnitude pre-filtering)
- **Literature**: Schechter fits from Bouwens+21 (φ* = 0.19 × 10⁻³,
  M* = −21.15, α = −2.06) and Finkelstein+15 (φ* = 1.57 × 10⁻⁴,
  M* = −21.03, α = −2.03), both at z ~ 7
- **Reference line**: Euclid limit at M_UV = −18

The 21cmFAST UVLF tracks the Schechter fits at the bright end but turns over
faintward of M_UV ≈ −14 where the M_TURN suppression acts. Bright-end
statistics are limited by the 256 Mpc box volume (see §12).

### 5.2 Stellar Mass – UV Magnitude Relation

- **Data**: all halos with M★ > 0 and SFR > 0, binned median + 16–84 % band
- **Literature**: Song+16 (`log M★ = 8.86 − 0.5(M_UV + 20)`) and González+10
  (`log M★ = 9.06 − 0.5(M_UV + 20)`) — **both now plotted**; the González
  line was previously commented out. Its ~0.2 dex higher normalisation
  follows from its constant-star-formation-history assumption.

### 5.3 Star-forming Main Sequence

- **Data**: all halos with M★ > 0 and SFR > 0 (no Euclid magnitude cut)
- **Three reference lines**:

| Line | Colour | Formula |
|---|---|---|
| 21cmFAST model | Green solid | log₁₀(SFR) = log₁₀(M★) − log₁₀(t_STAR × t_H) |
| Speagle+14 | Blue dashed | (0.84 − 0.026 t_age) log₁₀(M★) − (6.51 − 0.11 t_age) |
| Schreiber+15 | Red dash-dot | log₁₀(SFR) = log₁₀(M★) − 8  (sSFR = 10 Gyr⁻¹) |

Simulation data lie along the green model line, ~0.76 dex below Schreiber+15.
This offset is physical (different SFR timescales), not a modelling error.

---

## 6. SFR Model Comparison at z = 7

| Reference | sSFR | Timescale | Basis |
|---|---|---|---|
| Schreiber+15 high-z | 10 Gyr⁻¹ | 100 Myr | UV/Hα observations |
| Speagle+14 | ~8 Gyr⁻¹ | ~125 Myr | UV+IR compilation |
| 21cmFAST "simple" | **1.75 Gyr⁻¹** | **570 Myr** | t_STAR × t_H(z), Park+2019 |

The 21cmFAST model is calibrated to reproduce the reionisation history, not
the normalisation of the observed star-forming main sequence. The ~0.76 dex
offset is expected and documented in `docs/Low_SFR_fix.md`.

One coincidence worth flagging: the Schreiber+15 timescale is 100 Myr, the
same number that was hardcoded into the old bias model. They are not the same
quantity — 21cmFAST's internal SFR uses 570 Myr, and a bias estimate must use
the simulation's own value, not an observational one.

---

## 7. Analysis Stage Results *(superseded run)*

From `outputs/pipeline_summary.json`, produced by `run_pipeline.py` against the
2026-06-15 HDF5. The 21 cm quantities are unaffected by the bias fix; the
cross-spectrum and SNR will shift once the stronger Kaiser boost (§4.4) is
applied to a regenerated `galaxy_overdensity`.

| Quantity | Value |
|---|---|
| Mean neutral fraction ⟨x_HI⟩ | 0.176 |
| Large-scale cross-spectrum mean | −5.644 × 10³ (anti-correlated ✓) |
| Photo-z smearing σ_r | 20.6 Mpc (**157.5 Mpc** at the corrected σ_z = 0.45) |
| Horizon wedge slope | 3.151 |
| HERA FoV wedge slope | 0.379 |
| Modes outside the wedge | 105 / 400 (26.2 %) |
| **Total cross-correlation SNR** | **0.1 σ** (no detection) |

The 0.1 σ reflects the smoke-test geometry, not a physical forecast: a 3.5 Mpc
physical LOS extent gives essentially no independent k_∥ modes. Two known
issues also suppress it — see §10.

---

## 8. HDF5 Output (`outputs/lightcone_data.h5`)

| Dataset | Shape | Unit |
|---|---|---|
| `brightness_temp_field` | (128, 128, 100) | mK |
| `density_field` | (128, 128, 100) | dimensionless overdensity |
| `neutral_fraction` | (128, 128, 100) | [0, 1] |
| `galaxy_overdensity` | (128, 128, 100) | dimensionless (Kaiser RSD applied) |
| `lc_redshifts` | (100,) | — |
| `lc_dist_Mpc` | (100,) | Mpc |
| `halo_catalog/halo_masses` | (N_halos,) | M☉ |
| `halo_catalog/halo_coords` | (N_halos, 3) | Mpc |
| `halo_catalog/stellar_masses` | (N_halos,) | M☉ |
| `halo_catalog/sfr` | (N_halos,) | M☉ yr⁻¹ (unit-corrected) |

New root attributes this update: `galaxy_bias_method`,
`galaxy_bias_hmf_analytic`, `t_STAR`, `sfr_timescale_yr`.

---

## 9. Pipeline Status

| Step | File | Status |
|---|---|---|
| Driver — sim + analysis + figures + summary | `run_pipeline.py` | Complete |
| Part 1: Lightcone simulation + halo catalogue | `run_simulation.py` | Complete; **needs re-running** after the config and bias fixes |
| Part 2/3 computation | `src/analysis.py`, `src/figures.py`, `src/dataio.py` | Complete, 69 tests passing |
| Part 2: Field visualisation | `notebooks/plot_fields.ipynb` | Interactive front-end, current |
| Part 3: Power spectra & SNR | `notebooks/analysis.ipynb` | Interactive front-end; now also run headlessly by the pipeline |

---

## 10. Known Issues and Open Questions

1. **Mode counts are computed but unused in the SNR.** `mode_counts` is
   binned, cached in `analysis_products.h5`, and only ever consumed as an
   "empty bins" diagnostic. La Plante Eqs. 15–17 divide the per-bin variance
   by the number of modes; without that factor σ is overestimated and the
   total SNR is biased low. Needs checking against the paper.

2. **Thermal noise is a scaling estimate.**
   `src.analysis.hera_thermal_noise_power` uses T_sys²/(t_int Δν) with
   T_sys = 100 K + 60 K (300 MHz/ν)^2.55. The proper La Plante Eq. 11 form —
   with the baseline density n(k_⊥) and the X²Y Ω′ cosmological factors — is
   already implemented in `21cm_galaxy_cross_uncertainty.ipynb` and should be
   moved into `src/`. 21cmSense is the publication-grade alternative.

3. **`apply_rsds=False`** in `run_lightcone`; Kaiser RSD is applied
   analytically afterwards. Enabling 21cmFAST's self-consistent RSDs would
   remove the need for the β approximation entirely.

4. **Box size — deferred, see §12.**

5. **Bias model assumptions** (`docs/Galaxy_bias_formalism.md`): no dust
   attenuation, no scatter in the SFR–UV relation, deterministic halo
   occupation, bias determined by halo mass alone.

---

## 11. Changes Since the 2026-06-15 Update

| Change | Effect |
|---|---|
| Sheth-Tormen ν convention fixed in `run_simulation.py` | b_g 33.39 → 4.23 (analytic) |
| SFR timescale 100 Myr → t_STAR × t_H = 570 Myr | b_g 4.23 → 5.39 (analytic) |
| Bias now measured from the halo catalogue, not the mean relation | **adopted b_g = 4.744** |
| κ_UV and AB zero point unified via `src/conversions.py` | removes a 51.63 vs 51.60 zero-point and a 1.15e28 vs 8.696e27 κ discrepancy |
| Redshift range widened, then reverted | Stays at 6.995–7.005; widening is gated on `TODO.md` §P0 |
| β_rsd | 0.0299 → 0.2103 |
| González+10 M★–M_UV line re-enabled | figure now shows both literature relations |
| `run_pipeline.py` + `src/` modules + `tests/` added | analysis runs headlessly and is tested |

---

## 12. Future Work — 1 Gpc Box (deferred, needs HPC allocation)

Davies et al. (2025) use a 1 Gpc box at the same 2 Mpc cell size. The current
256 Mpc box is why the UVLF bright end is noisy (§5.1) and why only ~49 k
halos survive the Euclid cut. Scaling up is **deferred**: it is a compute and
storage problem, not a code change — `HII_DIM` and `BOX_LEN` are already
configuration variables at the top of `run_simulation.py`, and the pipeline's
`--max-halos` flag exists precisely for catalogues of this size.

| Quantity | Current (256 Mpc) | Target (1 Gpc) | Factor |
|---|---|---|---|
| `BOX_LEN` | 256 Mpc | 1000 Mpc | 3.9× |
| `HII_DIM` | 128 | 500 | — |
| `DIM` (high-res IC grid) | 384 | 1500 | — |
| Comoving volume | 1.7 × 10⁷ Mpc³ | 1.0 × 10⁹ Mpc³ | **59.6×** |
| Low-res cells | 2.10 × 10⁶ | 1.25 × 10⁸ | 59.6× |
| High-res IC cells | 5.66 × 10⁷ | 3.38 × 10⁹ | 59.6× |
| Halos (volume-scaled) | 1.14 × 10⁸ | **~6.8 × 10⁹** | 59.6× |
| Halo catalogue on disk | 2.74 GB | **~163 GB** | 59.6× |
| One lightcone field (N_z = 100) | 6.6 MB | 100 MB | 15.3× |
| Four lightcone fields | 46 MB | 0.70 GB | 15.3× |
| One high-res IC array | 0.23 GB | **13.5 GB** | 59.6× |

Storage is the binding constraint: the halo catalogue alone goes from 2.7 GB
to ~163 GB, and a single high-resolution initial-conditions array needs
13.5 GB of RAM. Both exceed a workstation. Before attempting this:

- Confirm the scratch quota can hold ~200 GB per run (see
  `docs/INSTALL_21cmFASTv4.md` for the CSD3 quota notes).
- Consider storing only the Euclid-selected halos rather than the full
  catalogue — at ~0.04 % selected that is ~65 MB instead of ~163 GB.
- Budget the run against the 128³ wall-clock timings before requesting an
  allocation.

An intermediate 512 Mpc box (`HII_DIM = 256`, 8× the volume, ~22 GB catalogue)
would test the scaling behaviour at a fraction of the cost.

---

## 13. Key Files

| File | Purpose |
|---|---|
| `run_pipeline.py` | Driver: simulation (optional) + analysis + figures + summary |
| `run_simulation.py` | Part 1: 21cmFAST lightcone + halo catalogue + Kaiser RSD |
| `src/analysis.py` | Power spectra, wedge, photo-z, noise, SNR, Euclid selection, bias |
| `src/figures.py` | All pipeline figures (headless) |
| `src/dataio.py` | HDF5 loading + analysis-product cache |
| `src/conversions.py` | SFR↔L_UV↔M_UV conversions, Sheth-Tormen bias, survey geometry |
| `src/FOV_to_cMpc.py` | Angular survey area to comoving Mpc conversion |
| `notebooks/plot_fields.ipynb` | Part 2 interactive front-end |
| `notebooks/analysis.ipynb` | Part 3 interactive front-end |
| `docs/Low_SFR_fix.md` | Diagnosis and fix for the SFR unit bug |
| `docs/Galaxy_bias_formalism.md` | Galaxy bias theory and implementation |
| `docs/halo_catalogue_reference.md` | 21cmFAST v4 halo catalogue API reference |

---

## 14. References

- Park, J. et al. (2019, MNRAS 484, 933) — 21cmFAST source model parameterisation
- Sheth, R. K. & Tormen, G. (1999, MNRAS 308, 119) — halo bias, ν convention
- Davies, J. et al. (2025, arXiv:2504.17254) — 21cmFASTv4 discrete source model
- La Plante, P. et al. (2023, arXiv:2205.09770) — cross-spectrum variance, wedge
- Bouwens, R. J. et al. (2021, AJ 162, 47) — UVLF at z~6–8
- Finkelstein, S. L. et al. (2015, ApJ 810, 71) — UVLF at z~7–8
- Speagle, J. S. et al. (2014, ApJS 214, 15) — star-forming main sequence
- Schreiber, C. et al. (2015, A&A 575, A74) — main sequence at high-z
- Madau, P. & Dickinson, M. (2014, ARA&A 52, 415) — UV–SFR calibration
- Song, M. et al. (2016, ApJ 825, 5) — M★–M_UV relation at z~7
- González, V. et al. (2010, ApJ 713, 115) — M★–M_UV at z~7
