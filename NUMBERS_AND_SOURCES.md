# Numbers and Sources

Audit of every hardcoded physical constant, calibration coefficient, and formula
in the pipeline, with its source citation.

**Scope.** One pass over the physics-relevant modules: `run_simulation.py`,
`run_pipeline.py`, `src/analysis.py`, `src/conversions.py`, `src/FOV_to_cMpc.py`.
Test files, `src/figures.py` (plotting only), `_archive/`, and `__pycache__/`
were skipped.

**Reading the Source column.** A citation means the value is traceable to that
reference. *Internal to pipeline config* means the value is this project's own
choice, not drawn from an external paper — legitimate, but not citable. **Source
not yet confirmed** means the value appears in the code and was not on the
supplied reference list; it is flagged rather than assigned a guessed citation.

---

## 1. Halo catalogue / source model

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `ALPHA_STAR` | `0.5` | `run_simulation.py` §4, `f_star()` | Park et al. (2019) SFE model, via Mesinger & Furlanetto (2007) and Murray et al. (2020) |
| `a` (Sheth-Tormen) | `0.707` | `conversions.sheth_tormen_bias()` | **Source not yet confirmed** — Sheth-Tormen (1999) fitting parameter, not on the supplied list |
| `delta_c` | `1.686` | `conversions.sheth_tormen_bias()` | **Source not yet confirmed** — linear collapse threshold, not on the supplied list |
| `dlog10m` | `0.02` | `analysis.effective_galaxy_bias()`, `run_simulation.py` §4 `MassFunction` | Internal to pipeline config — mass-function grid resolution [dex] |
| `F_STAR10` | `0.05` | `run_simulation.py` §4, `f_star()` | Park et al. (2019) SFE model, via Mesinger & Furlanetto (2007) and Murray et al. (2020) |
| `M_cell` (density grid, `DIM`) | `1.175e10` M☉ | `run_simulation.py` §1, `cell_mass(hires_cell_size, …)` | Internal to pipeline config — grid resolution mass |
| `M_cell` (21 cm / ionisation grid, `HII_DIM`) | `3.173e11` M☉ | `run_simulation.py` §1, `cell_mass(cell_size, …)` | Internal to pipeline config — grid resolution mass |
| `M_TURN` | `5e8` M☉ | `run_simulation.py` §4, `f_star()` | Star-formation turnover mass threshold (Park et al. 2019 SFE model) |
| `M_UV_bright` | `-22` | `run_simulation.py` config; `run_pipeline.py` `--m-uv-bright` | **Source not yet confirmed** — bright-end selection cut |
| `M_UV_limit` / `M_UV_faint` | `-18` | `run_simulation.py` config; `analysis.select_euclid_halos()` | **Source not yet confirmed** — Euclid faint magnitude cut |
| `OMEGA_B_0` | `0.049` | `run_simulation.py` §4, `stellar_mass_model()` | **Source not yet confirmed** — baryon density parameter |
| `p` (Sheth-Tormen) | `0.3` | `conversions.sheth_tormen_bias()` | **Source not yet confirmed** — Sheth-Tormen (1999) fitting parameter, not on the supplied list |
| `SAMPLER_MIN_MASS` | `1e8` M☉ | `run_simulation.py` §3a, read from `inputs.simulation_options` | Internal to pipeline config — halo sampler minimum resolved mass |
| `SIGMA_SFR_LIM` | `0.19` dex | `run_simulation.py` §4 (comment only; set inside 21cmFAST) | **Source not yet confirmed** — 21cmFAST log-normal SFR scatter |
| `SIGMA_STAR` | `0.25` dex | `run_simulation.py` §4 (comment only; set inside 21cmFAST) | **Source not yet confirmed** — 21cmFAST log-normal stellar-mass scatter |
| `T_STAR_DEFAULT` (`t_STAR`) | `0.5` | `analysis.star_formation_timescale()`; `run_simulation.py` §4 | Park et al. (2019), MNRAS 484, 933 — Eq. 3; SFR timescale as a fraction of the Hubble time |
| `USE_MINI_HALOS` | `False` | 21cmFAST `"simple"` template default (not set explicitly) | Internal to pipeline config — Population II sources only |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Sheth-Tormen bias | `b(ν̃) = 1 + (aν̃ − 1)/δ_c + 2p / (δ_c(1 + (aν̃)^p))` | `conversions.sheth_tormen_bias()` | **Source not yet confirmed** — Sheth & Tormen (1999), not on the supplied list |
| SFR timescale | `t_sf = t_STAR × t_H(z)` | `analysis.star_formation_timescale()` | Park et al. (2019), MNRAS 484, 933 — Eq. 3 |
| Stellar fraction | `f_*(M_h) = F_STAR10 (M_h/1e10)^ALPHA_STAR exp(−M_TURN/M_h)` | `run_simulation.py` §4, `f_star()` | Park et al. (2019) SFE model, via Mesinger & Furlanetto (2007) and Murray et al. (2020) |

---

## 2. Galaxy field

| Parameter | Value | Where used | Source |
|---|---|---|---|
| AB zero point | `51.60` | `conversions.Muv_to_Luv()`, `conversions.Luv_to_Muv()` | **Source not yet confirmed** — AB magnitude system zero point, not on the supplied list |
| `galaxy_bias` (fallback `b`) | `8` | `run_simulation.py` config; overwritten when the HMF estimate succeeds | **Source not yet confirmed** — "typical for high-z LBGs" per the inline comment |
| growth-rate exponent | `0.55` | `run_simulation.py` §5, `f(z) = Ω_m(z)^0.55` | **Source not yet confirmed** — standard ΛCDM growth-index approximation |
| `κ_UV` (`_KAPPA_UV_MADAU14`) | `1.15e-28` (M☉/yr)/(erg/s/Hz) | `conversions.sfr_to_Luv()`, `Luv_to_sfr()`, `sfr_to_Muv()` | Madau & Dickinson (2014), ApJ — UV-to-SFR calibration, used inverted as `L_UV = SFR/κ_UV` |
| `mean_galaxy_density` (`n̄`) | `3e-3` h³ Mpc⁻³ | `run_simulation.py` config; `analysis.compute_uncertainty_budget()` shot noise `P_N,gal = 1/n̄` | **Source not yet confirmed** — survey number density |
| `_RHO_CRIT_0_MSUN_MPC3_H2` | `2.77536627e11` M☉ Mpc⁻³ h⁻² | `conversions.mean_matter_density()` | **Source not yet confirmed** — critical density `ρ_crit,0 = 3H_0²/(8πG)` at `H_0 = 100h` |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| AB magnitude ↔ luminosity | `M_AB = −2.5 log10(L_ν) + 51.60` | `conversions.Muv_to_Luv()`, `Luv_to_Muv()` | **Source not yet confirmed** — not on the supplied list |
| Kaiser RSD boost | `δ_gal^(s)(k) = (1 + βμ²) δ_gal(k)`, `β = f/b` | `run_simulation.py` §5 | **Source not yet confirmed** — attributed to Kaiser (1987) in the inline comment, not on the supplied list |
| Luminosity-weighted overdensity | `δ_gal,L = Σ L_UV / ⟨Σ L_UV⟩ − 1` | `analysis.galaxy_overdensity_from_catalogue(weighting="luminosity")` | Internal to pipeline — consumes the Madau & Dickinson (2014) `κ_UV` conversion |
| Number-weighted overdensity | `δ_gal = N / ⟨N⟩ − 1` | `analysis.galaxy_overdensity_from_catalogue(weighting="number")` | Internal to pipeline |
| UV-to-SFR | `L_UV = SFR / κ_UV` | `conversions.sfr_to_Luv()` | Madau & Dickinson (2014), ApJ |

---

## 3. Noise model

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `bandwidth` (`B`, `Δν`) | `8e6` Hz (8 MHz) | `run_simulation.py` config; `analysis.hera_thermal_noise_power()` | **Source not yet confirmed** — per-band bandwidth |
| `dish_diameter` (`HERA_DISH_DIAMETER`) | `14.0` m | `run_simulation.py` config; `analysis.fov_wedge_slope()` | **Source not yet confirmed** — HERA dish diameter |
| `F_21_MHZ` / `f_21_hz` | `1420.405` MHz / `1420.405e6` Hz | `run_simulation.py` config; `analysis.system_temperature()`, `fov_wedge_slope()` | **Source not yet confirmed** — 21 cm rest frequency |
| `integration_time` (`t_obs`) | `1000 × 3600` s (1000 h) | `run_simulation.py` config; `analysis.hera_thermal_noise_power()` | **Source not yet confirmed** — total integration time |
| `NOISE_NORMALISATION_MPC3` | `1e3` | `analysis.hera_thermal_noise_power()` | Internal simplification — carries the `[Mpc³]` that makes `P_N` commensurate with `P_21`; stands in for the survey volume per mode. Explicitly *not* a physical constant |
| `SKY_SPECTRAL_INDEX` | `2.55` | `analysis.system_temperature()` | **Source not yet confirmed** — Galactic synchrotron spectral index |
| `T_RECEIVER_K` | `100.0` K | `analysis.system_temperature()` | **Source not yet confirmed** — HERA receiver temperature |
| `T_SKY_300MHZ_K` | `60.0` K | `analysis.system_temperature()` | **Source not yet confirmed** — Galactic synchrotron sky at 300 MHz |
| sky reference frequency | `300e6` Hz | `analysis.system_temperature()` | **Source not yet confirmed** — pivot frequency for the sky model |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Per-mode variance and SNR | Eqs. 15–17 | `analysis.cross_power_snr()`, `compute_uncertainty_budget()` | La Plante et al. (2023), ApJ, arXiv:2205.09770 |
| System temperature | `T_sys = T_rcvr + T_sky (300 MHz / ν)^2.55` | `analysis.system_temperature()` | **Source not yet confirmed** — component values not on the supplied list |
| Thermal noise power | `P_N = T_sys² × 1e3 / (t_obs × B)` | `analysis.hera_thermal_noise_power()` | Internal simplification — flat, `k`-independent. **Not** the full baseline-density interferometric model |

> **Flag — formalism referenced but not used.** The full baseline-density
> interferometric noise model (involving `N_ant` / `η_ap`, resolving
> `P_N(k_perp)` through the baseline density) is **not implemented anywhere in
> this pipeline**. It survives only as a forward-reference in prose:
> `src/analysis.py:705` (docstring: "or with La Plante et al. (2023) Eq. 11,
> which resolves `P_N(k_perp)` through the baseline density"),
> `docs/HPC.md:781`, and `docs/project_update.md:323`. The implemented noise is
> the flat `T_sys²×1e3/(t_obs B)` scaling estimate above. Any statement that
> this pipeline uses the `N_ant`/`η_ap` formalism is incorrect.

---

## 4. Foreground treatment

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `n_bins_parallel` | `20` | `run_simulation.py` config; `analysis.compute_cylindrical_cross_power()` | Internal to pipeline config — power-spectrum binning |
| `n_bins_perp` | `20` | `run_simulation.py` config; `analysis.compute_cylindrical_cross_power()` | Internal to pipeline config — power-spectrum binning |
| `wedge_buffer` | `0.0677` Mpc⁻¹ (= 0.1 h Mpc⁻¹ at h = 0.6766) | `run_simulation.py` config; `analysis.foreground_wedge_mask()`; `run_pipeline.py --wedge-buffer` | Pober et al. (2014), ApJ 782, 66 — "moderate" foreground model; matches 21cmSense's default `horizon_buffer` |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| FoV wedge slope | primary-beam-limited variant of the horizon slope | `analysis.fov_wedge_slope()` | Thyagarajan et al. (2015), ApJ 804, 14; La Plante et al. (2023), Eq. 10 |
| Horizon wedge slope | `m(z) = D_C(z) H(z) / [c(1+z)]` | `analysis.horizon_wedge_slope()` | Thyagarajan et al. (2015), ApJ 804, 14; same formula in La Plante et al. (2023), Eq. 10 |
| Wedge mask | `k_∥ > m(z) k_⊥ + buffer` | `analysis.foreground_wedge_mask()` | Thyagarajan et al. (2015) for the slope; Pober et al. (2014) for the buffer |

---

## 5. Photo-z damping

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `PHOTOZ_N_SIGMA` | `1` | `run_simulation.py` config — sets `SURVEY_DELTA_Z = 2 × N_σ × σ_z` | Internal to pipeline config — explicitly labelled a choice, not a default |
| `photoz_uncertainty` (`σ_z`) | `0.45` (absolute, at `z_obs = 7`) | `run_simulation.py` config; `analysis.radial_smearing_length()`; `run_pipeline.py --sigma-z` | Euclid fractional photo-z requirement `σ_z/(1+z) < 0.05`, evaluated at z = 7 (gives `σ_z/(1+z) = 0.056`) |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Damping kernel | `W(k_∥) = exp(−½ k_∥² σ_r²)` | `analysis.photoz_damping_kernel()` | This project's paper, Eq. 27 |
| Radial smearing length | `σ_r = c σ_z / H(z)` | `analysis.radial_smearing_length()` | **Source not yet confirmed** — standard photo-z-to-comoving conversion, not on the supplied list |

---

## 6. Survey geometry and cosmology

Cross-cutting values that do not belong to a single pipeline stage.

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `_DIM_PER_HII_DIM` | `3` | `conversions.survey_area_to_box_size()` | Internal to pipeline config — 21cmFAST convention for the high-res grid ratio |
| `HUBBLE_CONSTANT` (`H_0`) | `67.36` km s⁻¹ Mpc⁻¹ | `run_simulation.py` config; default across `src/analysis.py` | **Source not yet confirmed** — labelled "Planck 2018" in the inline comment |
| `minimum_los_slices` | `100` | `run_simulation.py` §1 | Internal to pipeline config — LOS slice floor for the smoke-test slab |
| `n_nodes` scaling | `10` nodes per unit redshift (min 5) | `run_simulation.py` §1 | Internal to pipeline config — 21cmFAST node redshift density |
| `OMEGA_M_0` (`Ω_m`) | `0.315` | `run_simulation.py` config; default across `src/analysis.py` | **Source not yet confirmed** — labelled "Planck 2018" in the inline comment |
| `random_seed` | `42` | `run_simulation.py` §2, `InputParameters.from_template()` | Internal to pipeline config |
| `SPEED_OF_LIGHT_KMS` / `_MPS` | `3e5` km s⁻¹ / `3e8` m s⁻¹ | `run_simulation.py` config; defaults across `src/analysis.py` | Physical constant (rounded to 2 s.f.; the exact value is 2.998e5 km s⁻¹) |
| `SURVEY_AREA_DEG2` | `10.0` deg² | `run_simulation.py` config → `conversions.survey_area_to_box_size()` | EDF-Fornax field definition (ESA / Euclid Consortium, 2019 field announcement) |
| `SURVEY_Z_CENTRAL` | `7.0` | `run_simulation.py` config | Internal to pipeline config — central analysis redshift |
| `_TARGET_CELL_SIZE_MPC` | `2.0` Mpc | `conversions.survey_area_to_box_size()` | Internal to pipeline config — preserves the mass resolution of the former 256 Mpc / 128³ grid |
| `z_min`, `z_max` | `6.995`, `7.005` | `run_simulation.py` config | Internal to pipeline config — deliberate quasi-coeval smoke-test slab (see `TODO.md` P0) |
| `_KM_PER_MPC` | `3.0856775814913673e19` | `analysis.star_formation_timescale()` | Physical constant — Mpc in km |
| `_SEC_PER_YR` | `365.25 × 24 × 3600` | `analysis.star_formation_timescale()`; `run_simulation.py` §3a | Physical constant — Julian year in seconds |
| EDF-Fornax centre | RA 03:31:43.6, Dec −28:05:18.6 | `run_simulation.py:101`, `conversions.py:718` — **comments only**, never computed on | EDF-Fornax field definition (ESA / Euclid Consortium, 2019 field announcement) |

---

## Values on the reference list that do not appear in the code

| Value | Source | Status |
|---|---|---|
| HERA declination stripe, δ ≈ −30° ± 5° | La Plante et al. (2023), ApJ | **Not present** anywhere in the codebase — no `.py` or `.md` file references a declination. The EDF-Fornax centre (Dec −28:05:18.6) sits inside this stripe, but the pipeline never checks or uses the overlap |

---

## Summary of flags

**Sources not yet confirmed (24).** Concentrated in three places: the
Sheth-Tormen bias parameters (`delta_c`, `a`, `p` and the bias formula itself),
the HERA instrument and sky-temperature values (`T_RECEIVER_K`,
`T_SKY_300MHZ_K`, `SKY_SPECTRAL_INDEX`, `dish_diameter`, `bandwidth`,
`integration_time`), and the cosmology (`OMEGA_M_0`, `HUBBLE_CONSTANT`, labelled
"Planck 2018" in comments but not cited). Also the AB zero point `51.60`, the
magnitude cuts `M_UV_limit = -18` / `M_UV_bright = -22`, `OMEGA_B_0 = 0.049`,
`mean_galaxy_density = 3e-3`, the fallback `galaxy_bias = 8`, the growth-index
exponent `0.55`, the Kaiser (1987) RSD formula, `σ_r = cσ_z/H(z)`, and the
21cmFAST scatter values `SIGMA_STAR` / `SIGMA_SFR_LIM`.

**Internal, non-citable by design (16).** `SAMPLER_MIN_MASS`, both `M_cell`
values, `NOISE_NORMALISATION_MPC3`, the binning counts, grid-geometry constants,
the smoke-test redshift slab, and the two `δ_gal` normalisation formulas.

**One formalism referenced but not implemented.** The baseline-density
interferometric noise model (`N_ant` / `η_ap`) appears only in docstrings and
docs prose — see the flag in §3.
