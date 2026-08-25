# Numbers and Sources

Audit of every hardcoded physical constant, calibration coefficient, and formula
in the pipeline, with its source citation.

**Scope.** One pass over the physics-relevant modules: `run_simulation.py`,
`run_pipeline.py`, `src/analysis.py`, `src/conversions.py`, `src/FOV_to_cMpc.py`,
`src/foregrounds.py`, `src/provenance.py`.
Test files, `src/figures.py` (plotting only), `_archive/`, and `__pycache__/`
were skipped. The post-Euclid-cut figures added on 2026-08-21 introduce no new
physical constants — their `slab_cells`, `smooth_cells` and contour quantiles
are display-only choices, documented as such in the function docstrings and
deliberately not listed here.

**Updated 2026-08-24.** The two `M_cell` entries in §1 were still quoting the
retired 256 Mpc / 128³ grid and now carry the footprint-derived values; §12 is
the derived numbers of the planned production run
([`docs/simulation_spec.md`](docs/simulation_spec.md)); §§10–11, added
2026-08-25, cover the lightcone estimator and the smoke-test overrides.

**Updated 2026-08-21.** Several entries previously marked *Source not yet
confirmed* have been traced and are now cited: the `T_sys` model (21cmSense),
the dish diameter and array layout (DeBoer et al. 2017), and the 21 cm rest
frequency. Sections 7–9 are new, covering the physical noise model, the
foreground module, and the run-manifest calibration.

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
| `M_cell` (density grid, `DIM`) | `1.007e10` M☉ | `run_simulation.py` §1, `cell_mass(hires_cell_size, …)` | Internal to pipeline config — grid resolution mass. **Updated 2026-08-24**: `DIM` = 768 at `BOX_LEN` = 486.33 Mpc (0.6332 Mpc cell); was `1.175e10` on the retired 256 Mpc / `DIM` = 384 grid |
| `M_cell` (21 cm / ionisation grid, `HII_DIM`) | `2.720e11` M☉ | `run_simulation.py` §1, `cell_mass(cell_size, …)` | Internal to pipeline config — grid resolution mass. **Updated 2026-08-24**: `HII_DIM` = 256 at `BOX_LEN` = 486.33 Mpc (1.8997 Mpc cell); was `3.173e11` on the retired 256 Mpc / 128³ grid |
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
| `dish_diameter` (`HERA_DISH_DIAMETER`) | `14.0` m | `run_simulation.py` config; `analysis.fov_wedge_slope()`, `hera_beam_solid_angles()`, `hera_baseline_counts()` | DeBoer et al. (2017), PASP 129, 045001 — "14-m parabolic dishes" |
| `F_21_MHZ` / `f_21_hz` | `1420.405` MHz / `1420.405e6` Hz | `run_simulation.py` config; `analysis.system_temperature()`, `fov_wedge_slope()`, `cosmological_scalar_x2y()` | Hyperfine transition of neutral hydrogen, 1420.405751768 MHz (CODATA / Hellwig et al. 1970); the pipeline rounds to 1420.405 |
| `integration_time` (`t_obs`) | `1000 × 3600` s (1000 h) | `run_simulation.py` config; `analysis.hera_thermal_noise_power()` | **Source not yet confirmed** — total integration time |
| `NOISE_NORMALISATION_MPC3` | `1e3` | `analysis.hera_thermal_noise_power()` | Internal simplification — carries the `[Mpc³]` that makes `P_N` commensurate with `P_21`; stands in for the survey volume per mode. Explicitly *not* a physical constant |
| `SKY_SPECTRAL_INDEX` | `2.55` | `analysis.system_temperature()` | Pober et al. (2013), AJ 145, 65; Pober et al. (2014), ApJ 782, 66 — as implemented in 21cmSense `calc_sense.py`: `Tsky = 60e3 * (3e8/freq)**2.55` |
| `T_RECEIVER_K` | `100.0` K | `analysis.system_temperature()` | 21cmSense `calc_sense.py` (`Trx`, from the array file); also the `T_rcvr` term of DeBoer et al. (2017) Table 2, `T_sys = 100 + 120(ν/150 MHz)^-2.55` K |
| `T_SKY_300MHZ_K` | `60.0` K | `analysis.system_temperature()` | 21cmSense `calc_sense.py`: `Tsky = 60e3 * (3e8/freq)**2.55` mK, i.e. 60 K × (λ/1 m)^2.55. **Differs from DeBoer et al. (2017) Table 2** — see the flag below |
| sky reference frequency | `300e6` Hz | `analysis.system_temperature()` | 21cmSense convention — 300 MHz is where λ = 1 m, so `(300 MHz/ν)^2.55 = (λ/1 m)^2.55` |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Per-mode variance and SNR | Eqs. 15–17 | `analysis.cross_power_snr()`, `compute_uncertainty_budget()` | La Plante et al. (2023), ApJ, arXiv:2205.09770 |
| System temperature | `T_sys = T_rcvr + T_sky (300 MHz / ν)^2.55` | `analysis.system_temperature()` | Pober et al. (2013, 2014), as implemented in 21cmSense. Gives 328.6 K at z = 7 (ν = 177.6 MHz) |
| Per-bin mode weighting | `ŝ = √(N_patch dN) P_× / σ_×`, `dN = k_⊥² k_∥ V (2π)⁻² dln k_⊥ dln k_∥` | `analysis.cross_power_snr(mode_counts=…)` | La Plante et al. (2023), Eqs. 18–19. `dN` is `mode_counts / 2`; verified against Eq. 18 to a median 4.7 % on the fiducial grid |
| Thermal noise power (default) | `P_N = T_sys² × 1e3 / (t_obs × B)` | `analysis.hera_thermal_noise_power()` | Internal simplification — flat, `k`-independent, and ~10⁴ below the instrument model. **Not** the baseline-density model |
| Thermal noise power (physical) | `P_N(k_⊥) = X²Y Ω_eff · NEB · T_sys² / (N_pol t_int N_bl(k_⊥))` | `analysis.hera_thermal_noise_power_physical()` | Parsons (2017), "Power Spectrum Normalizations for HERA", Eq. 12 — algebraically identical to La Plante et al. (2023) Eq. 11 since `Ω_eff ≡ Ω_P²/Ω_PP` |

> **Flag resolved 2026-08-21 — the baseline-density model now exists.** It was
> previously referenced in prose but implemented nowhere. It is now
> `analysis.hera_thermal_noise_power_physical()`, selected with
> `--noise-model physical` (§7). **The default remains the flat
> `T_sys²×1e3/(t_obs B)` scaling estimate**, so any result not explicitly
> produced with `--noise-model physical` still uses the simplification, which
> is ~10⁴ times smaller than the instrument model.

> **Flag — the `T_sys` normalisation is convention-dependent.** This pipeline
> uses 21cmSense's `T_sky = 60 K (λ/1 m)^2.55`, giving 352 K at 150 MHz.
> DeBoer et al. (2017) Table 2 give `T_sys = 100 + 120(ν/150 MHz)^-2.55` K,
> whose sky term is 120 K at 150 MHz — **a factor 2.9 colder**. Noise power
> scales as `T_sys²`, so the choice is worth ~8× in `P_N`. The 21cmSense
> values are kept because they are what the notebook inlined and are the more
> conservative choice; recorded here so the difference is not mistaken for an
> error.

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

## 7. HERA instrument model (physical noise path)

Added 2026-08-21 with `analysis.hera_thermal_noise_power_physical()`. Used only
when `--noise-model physical` is passed; the default path does not touch these.

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `HERA_ANTENNA_SPACING_M` | `14.6` m | `analysis.hera_baseline_counts()` | DeBoer et al. (2017), PASP 129, 045001 — "14.6 m center-to-center spacing" of the hexagonal core |
| `HERA_APERTURE_EFFICIENCY` (`η_ap`) | `0.65` | `analysis.hera_beam_solid_angles()` | Calibrated so `Ω_P = λ²/A_e` with `A_e = η_ap πD²/4` reproduces the `Ω_P ≈ 0.04` sr at 150 MHz plotted in Parsons (2017); gives `A_e = 100` m² against a geometric 154 m² |
| `HERA_HEX_N_SIDE` | `11` | `analysis.hera_baseline_counts()` | Approximates HERA's 320-element core: `3n²−3n+1 = 331` at `n = 11` is the nearest hex number. DeBoer et al. (2017) — 350 total = 320 core + 30 outriggers; outriggers not modelled |
| `HERA_OMEGA_P_OVER_PP` (`Ω_P/Ω_PP`) | `2.175` | `analysis.hera_beam_solid_angles()` | Parsons (2017), "Power Spectrum Normalizations for HERA" — prints `HERA Omega_P/OMEGA_PP 2.1752891255`, the median over CST beam models. (PAPER's is 2.35) |
| `N_POLARISATIONS` (`N_pol`) | `2` | `analysis.hera_thermal_noise_power_physical()` | Parsons (2017), Eq. 12 — "the factor of `N_pol` … explicitly counts the two orthogonal polarizations" |
| `NOISE_EQUIVALENT_BANDWIDTH` (NEB) | `1.0` | `analysis.hera_thermal_noise_power_physical()` | Parsons (2017), Eq. 12 with `WINDOW = 'none'`. Exactly 1 because `compute_cylindrical_cross_power()` applies no line-of-sight taper |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Antenna theorem | `Ω_P A_e = λ²` | `analysis.hera_beam_solid_angles()` | Thompson, Moran & Swenson (2017), *Interferometry and Synthesis in Radio Astronomy*, 3rd ed. |
| Cosmological scalar | `X²Y = D_c(z)² · c(1+z)²/[H(z) f_21]` | `analysis.cosmological_scalar_x2y()` | Parsons et al. (2012a), ApJ 756, 165; Parsons (2017), Eq. 1. **1227 Mpc³ sr⁻¹ Hz⁻¹ at z = 7** |
| Effective beam area | `Ω_eff ≡ Ω_P²/Ω_PP` | `analysis.hera_beam_solid_angles()` | Parsons et al. (2014), ApJ 788, 106, Appendix B — using `Ω_P` alone is a known power-spectrum normalisation error |
| Baseline gridding | baselines binned into `uv`-cells of side `D/λ` | `analysis.hera_baseline_counts()` | Parsons et al. (2012a) — only baselines within one antenna footprint sample the same mode and integrate down coherently |

**Values this produces at z = 7 (ν = 177.6 MHz), 1000 h**

| Quantity | Value |
|---|---|
| `X²Y` | 1227 Mpc³ sr⁻¹ Hz⁻¹ |
| `Ω_P`, `Ω_eff` | 0.0285 sr, 0.0621 sr |
| Baselines per `uv`-cell | 290 at `k_⊥ = 0.010` → 12 at `k_⊥ = 0.123` |
| `P_N,21` | 3.9 × 10³ → 9.8 × 10⁴ mK² Mpc³ |
| `P_N,21` beyond `k_⊥ ≈ 0.13` | `inf` — the 292 m core has no baselines there |

Cross-check: a published 1000 h HERA forecast of `Δ²_N ~ 10` mK² at `k = 0.2`
is `P_N = 2π²Δ²/k³ ≈ 2.5 × 10⁴` mK² Mpc³ — inside the range above. The default
scaling estimate gives 3.75 mK² Mpc³, roughly 10⁴ lower.

> **Idealisations that remain.** Every baseline is assumed to integrate for the
> full `t_int` on one field; `uv`-plane rotation is ignored; the outriggers are
> not modelled, so `k_⊥ > 0.13` Mpc⁻¹ is reported as unmeasurable when in
> reality the outriggers sample it sparsely. For a publication forecast use
> [21cmSense](https://github.com/rasg-affiliates/21cmSense).

---

## 8. Foreground injection and removal

`src/foregrounds.py`, added 2026-08-21. Notebook §7e only; not part of the
batch pipeline.

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `beta_angular` (diffuse) | `2.4` | `DIFFUSE_DEFAULTS`, `simulate_diffuse_foreground()` | Santos, Cooray & Knox (2005), ApJ 625, 575 — Table 1, Galactic synchrotron angular index `C_ℓ ∝ ℓ^-β` |
| `beta_angular` (point source) | `1.1` | `POINT_SOURCE_DEFAULTS` | Santos, Cooray & Knox (2005), Table 1, point-source row |
| `contrast` | `0.5` | `DIFFUSE_DEFAULTS` | Internal to module config — log-normal width of the angular temperature fluctuations |
| `flux_slope` (`γ`) | `1.75` | `POINT_SOURCE_DEFAULTS` | Differential source count `dN/dS ∝ S^-γ`; Ali, Bharadwaj & Chengalur (2008), MNRAS 385, 2166 |
| `reference_frequency_mhz` | `130.0` MHz | both defaults | Santos, Cooray & Knox (2005) — the `ν_ref` of their Table 1 parametrisation |
| `reference_temperature_mK` (diffuse) | `700e3` mK (700 K at 130 MHz) | `DIFFUSE_DEFAULTS` | Global Sky Model: de Oliveira-Costa et al. (2008), MNRAS 388, 247; Zheng et al. (2017), MNRAS 464, 3486 (GSM2016) |
| `reference_temperature_mK` (point source) | `57e3` mK (57 K at 130 MHz) | `POINT_SOURCE_DEFAULTS` | Santos, Cooray & Knox (2005), Table 1 point-source amplitude |
| `source_density_per_cell` | `0.02` | `POINT_SOURCE_DEFAULTS` | Internal to module config — sets the sparsity of the drawn population |
| `spectral_index` (diffuse, `α`) | `2.8` | `DIFFUSE_DEFAULTS` | Santos, Cooray & Knox (2005), Table 1; consistent with the GSM synchrotron index |
| `spectral_index` (point source) | `2.07` | `POINT_SOURCE_DEFAULTS` | Santos, Cooray & Knox (2005), Table 1 |
| `spectral_index_scatter` (`σ_α`) | `0.1` (diffuse), `0.3` (point) | both defaults | Shaw et al. (2014), ApJ 781, 57 — spatial variation of the spectral index, the term that breaks separability |

**Formulas**

| Formula | Expression | Where | Source |
|---|---|---|---|
| Foreground spectrum | `T(θ,ν) = T_ref(θ) (ν/ν_ref)^-α(θ)` | `_apply_spectral_law()` | Santos, Cooray & Knox (2005); Shaw et al. (2014) |
| Parametrised covariance | `C_ℓ(ν₁,ν₂) = A (ℓ/ℓ_ref)^-β (ν₁ν₂/ν_ref²)^-α` | model basis for the above | Santos, Cooray & Knox (2005), their §3 |
| Removal | `field − f × foreground` (amplitude basis) | `remove_foreground()` | **Not a literature method.** A deliberately simple placeholder; see the docstring |

> **Flag — `remove_foreground` is not citable as a method.** It subtracts an
> exactly-correct template of the injected field. It is not GMCA, PCA, ICA,
> polynomial fitting, Gaussian-process removal or delay filtering, and results
> from it describe a *removal level*, not any named algorithm's performance.

> **Flag — smoothness does not confine the foreground to low `k_∥`.** Measured
> through this pipeline's un-tapered FFT, a smooth-but-non-periodic spectrum
> leaks along the whole `k_∥` axis as ≈`k_∥^-1.5`. The slope comes from the
> window, not the sky: a bare ramp with no angular structure leaks identically.

---

## 9. Run-manifest calibration

`src/provenance.py`, added 2026-08-21. Used only for the pre-flight
halo-catalogue cost estimate; nothing downstream consumes it.

| Parameter | Value | Where used | Source |
|---|---|---|---|
| `BYTES_PER_HALO` | `28.0` bytes | `estimate_catalogue_cost()` | Measured from this project's own 2026-08-12 run: a 3.564 GiB `HaloCatalog.h5` holding 136,663,818 halos |
| `HALOS_PER_MPC3` | `8.146` Mpc⁻³ | `estimate_catalogue_cost()` | Measured: 136,663,818 halos in a (256 Mpc)³ box at z = 7 with `SAMPLER_MIN_MASS = 1e8` |
| `INT32_MAX` | `2³¹ − 1` | `estimate_catalogue_cost()` | Largest signed 32-bit index; 21cmFAST's C backend indexes halo arrays with `int` |
| `PERTURBED_FRACTION` | `0.836` | `estimate_catalogue_cost()` | Measured: 114,289,081 / 136,663,818 from the 2026-08-12 run, reproduced independently as 1,782,540 / 2,135,000 by a 64 Mpc run on 2026-08-21 |

These are **internal empirical calibrations, not literature values**. They
extrapolate this project's own measurements and are accurate only for runs
sharing the reference run's redshift, cosmology and sampler settings.

---

## 10. Lightcone estimator (TODO.md P0)

`src/analysis.py`, added 2026-08-25.  Active only under
`ESTIMATOR = "lightcone"`; the coeval path uses none of them.

| Parameter | Value | Where used | Source |
|---|---|---|---|
| Blackman-Harris coefficients | `0.35875, 0.48829, 0.14128, 0.01168` | `blackman_harris_taper()` | Harris, F. J. (1978), Proc. IEEE 66, 51 — Table 1, the 4-term −92 dB window; standard in 21 cm delay-spectrum work via Parsons et al. (2012), ApJ 756, 165 |
| Taper normalisation ⟨w²⟩ | `0.2580` (asymptotic; `0.25746` at N = 512) | `compute_cylindrical_cross_power()` | Derived — `a₀² + (a₁² + a₂² + a₃²)/2`, the noise-equivalent bandwidth that restores the amplitude of a homogeneous field |
| `min_slices_per_band` | `8` | `subband_index_ranges()` | Internal to pipeline config — floor below which a band has no usable k_∥ axis |
| Band count | `ceil(frequency span / bandwidth)` | `subband_index_ranges()` | Derived — guarantees each band spans ≤ the noise bandwidth (P0.4) |
| Effective redshift | `z_eff = f_21 / ⟨ν_band⟩ − 1` | `compute_subband_power_spectra()` | Derived — the mean *observed frequency*, not the mean redshift; the two differ by ~0.2 % over 8 MHz at z = 7 |
| Band combination | `sqrt(Σ SNR_band²)` | `combine_band_snr()` | Derived — sub-bands sample disjoint comoving volumes, so their measurements are independent |

---

## 11. Smoke-test overrides

`src/smoke_test.py`, added 2026-08-25.  **Not scientific values**, and never
active unless `--smoke-test` is passed.  Listed here so the audit is complete
and so nobody mistakes one for a production parameter.

| Parameter | Smoke value | Production value | Source |
|---|---|---|---|
| `HII_DIM` / `BOX_LEN` / `DIM` | `16` / `32.0` Mpc / `48` | `256` / `486.33` Mpc / `768` | Internal to smoke-test config — chosen for runtime, not physics |
| `minimum_los_slices` | `12` | `100` | " |
| `n_bins_perp` / `n_bins_parallel` | `8` / `8` | `20` / `20` | " |
| `max_halos` | `200000` | `0` (all) | " |

Everything else — the survey footprint, integration time, bandwidth, σ_z,
the wedge buffer, the magnitude cuts and the random seed — is **unchanged**
under `--smoke-test`; `SMOKE_TEST_UNCHANGED` records why in each case.

---

## 12. Planned production run — derived numbers

Computed 2026-08-24 for [`docs/simulation_spec.md`](docs/simulation_spec.md),
which specifies the run at *z* = 6.55–7.45 in the footprint-derived box. These
are **derived**, not configured: nothing in the code stores them, and they
change the moment `z_min`/`z_max`, `BOX_LEN` or the footprint does. Marked
**[C]** computed here, **[E]** extrapolated from a measured run,
**[M]** measured.

| Quantity | Value | How obtained |
|---|---|---|
| `BOX_LEN` / `HII_DIM` / `DIM` | 486.329 Mpc / 256 / 768 | `conversions.survey_area_to_box_size(10.0, 7.0, 0.90)` **[C]** |
| Transverse / high-res cell | 1.8997 / 0.6332 Mpc | `BOX_LEN` / grid **[C]** |
| Survey LOS depth L∥ | 315.598 Mpc | `Planck18.comoving_distance` differenced over 6.55–7.45 **[C]** |
| `N_z` / LOS slice spacing | 166 / 1.9012 Mpc | `round(L_los / cell)`, floor idle **[C]** |
| LOS spacing spread, first → last cell | 2.0803 → 1.7599 Mpc = 16.75 % | `np.diff` of `Planck18.comoving_distance(linspace(6.55, 7.45, 166))` **[C]** — quantifies `TODO.md` §P0.1 at the planned range |
| Node redshifts | 9, step 0.1125, *z* = 7.0 exactly on node 5 | `max(round(10 Δz), 5)` **[C]** |
| Observed frequency span | 168.095 – 188.133 MHz = 20.038 MHz | `F_21_MHZ / (1 + z)` **[C]** — a 2.50× mismatch against the 8 MHz noise `bandwidth` (`TODO.md` §P0.4) |
| Δk⊥ / Δk∥ | 0.01292 / 0.01991 Mpc⁻¹ | 2π/L **[C]** |
| k_Nyq,⊥ / k_Nyq,∥ | 1.6537 / 1.6524 Mpc⁻¹ | π/cell **[C]** |
| Recorded `L_los` vs true LOS span | 315.354 vs 315.598 Mpc = 0.08 % | `N_z` × transverse cell **[C]** — `HPC.md` §11.1's 56.5× discrepancy is an artifact of the floored slab and does not recur here |
| Photo-*z* kernel W at k∥ = 0.0199 / 0.1084 Mpc⁻¹ | 7.3 × 10⁻³ / 5.2 × 10⁻⁶⁴ | exp(−k∥²σ_r²/2) at σ_r = 157.48 Mpc **[C]** |
| Halos, Lagrangian / perturbed | 9.37 × 10⁸ / 7.84 × 10⁸ | `provenance.estimate_catalogue_cost(486.33)` **[E]** |
| Catalogue on disk / resident while perturbing | 26.2 / 48.2 GB | same **[E]** |
| `int32_headroom` | **1.31 — over `INT_MAX`** | same **[E]** |
| Box length at headroom 1.0 / 0.5 | 444.6 / 352.9 Mpc | inverting `HALOS_PER_MPC3` **[C]** |
| Halo count retained by `SAMPLER_MIN_MASS` 1.5 / 2 / 3 × 10⁸ M☉ | 63.1 % / 46.2 % / 28.7 % | Sheth-Tormen `hmf.MassFunction(z=7, dlog10m=0.02)` cumulative counts **[C]** |
| Star formation retained by the same floors | 99.95 % / 99.84 % / 99.44 % | same, weighted by M · f_★ with f_★ ∝ (M/10¹⁰)^0.5 exp(−`M_TURN`/M) **[C]** |
| Stored HDF5 | ~19.1 GB (18.8 GB catalogue at 24 B/halo) | **[E]** |
| 21cmFAST cache | ~245 GB (7.7 GB ICs + 26.4 GB × 9 nodes) | scaled from the measured 0.96 GB / 3.83 GB-per-node at 256 Mpc **[E]** |
| Peak resident memory | ~56 GB | **[E]** |
| Serial simulation cost vs the 256 Mpc baseline | ~12.3× → 2–4 h on one core | 6.86 (volume) × 1.8 (nodes) applied to a measured 520 s Stage 1 **[E]** |
| Measured baseline wall / CPU (256 Mpc, 128³, 5 nodes, 1 thread) | 543.5 s / 567.4 s = 0.163 CPU-h | 2026-08-12 run, reproduced by 2026-08-07 to within 1 s **[M]** |
| Target machine CPU | 2 × AMD EPYC 9374F, 64 physical / 128 logical cores, 3.85 GHz base | `lscpu`, 2026-08-24 **[M]** |
| Target machine RAM | 1.5 TiB total, 1.4 TiB available, 4 GiB swap | `free -h`, 2026-08-24 **[M]** |
| Target machine host | `andromeda1.jb.man.ac.uk` — the same host as the 2026-08-12 baseline | confirmed 2026-08-24 **[M]**; the serial runtime scaling therefore carries no cross-machine correction |
| Target machine scratch | `/nvme1` 851 GB free, `/nvme4` 836 GB free (XFS on NVMe); `$HOME` NFS, 144 GB free | `df -hT`, 2026-08-24 **[M]** |
| Target machine scheduler | none — `sinfo`/`scontrol`/`sacct`/`sacctmgr` absent | 2026-08-24 **[M]** |
| Planned-run wall time | 15–40 min at 32 threads; budget 1 h | 6,400 s serial ÷ an assumed 8–14× OpenMP speed-up **[E]** |

**[M]** rows are citable as this project's own measurements. **[E]** rows
inherit the caveat on §9: they extrapolate one run's sampler settings,
redshift and cosmology, and are guards, not budgets.

---

## Values on the reference list that do not appear in the code

| Value | Source | Status |
|---|---|---|
| HERA declination stripe, δ ≈ −30° ± 5° | La Plante et al. (2023), ApJ | **Not present** anywhere in the codebase — no `.py` or `.md` file references a declination. The EDF-Fornax centre (Dec −28:05:18.6) sits inside this stripe, but the pipeline never checks or uses the overlap |

---

## Summary of flags

**Newly cited on 2026-08-21 (6).** `T_RECEIVER_K`, `T_SKY_300MHZ_K`,
`SKY_SPECTRAL_INDEX` and the 300 MHz pivot are traced to 21cmSense
(`calc_sense.py`, `Tsky = 60e3 * (3e8/freq)**2.55`), themselves from Pober
et al. (2013, 2014); `dish_diameter` and the array layout to DeBoer et al.
(2017); `F_21_MHZ` to the hydrogen hyperfine frequency. Sections 7–9 add 21
further cited values for the HERA instrument model and the foreground module.

**Sources not yet confirmed (18).** Concentrated in two places: the
Sheth-Tormen bias parameters (`delta_c`, `a`, `p` and the bias formula itself)
and the cosmology (`OMEGA_M_0`, `HUBBLE_CONSTANT`, labelled "Planck 2018" in
comments but not cited). Also the AB zero point `51.60`, the magnitude cuts
`M_UV_limit = -18` / `M_UV_bright = -22`, `OMEGA_B_0 = 0.049`,
`mean_galaxy_density = 3e-3`, `bandwidth`, `integration_time`, the fallback
`galaxy_bias = 8`, the growth-index exponent `0.55`, the Kaiser (1987) RSD
formula, `σ_r = cσ_z/H(z)`, and the 21cmFAST scatter values `SIGMA_STAR` /
`SIGMA_SFR_LIM`.

**Internal, non-citable by design (22).** `SAMPLER_MIN_MASS`, both `M_cell`
values, `NOISE_NORMALISATION_MPC3`, the binning counts, grid-geometry
constants, the smoke-test redshift slab, the two `δ_gal` normalisation
formulas, the four `provenance.py` calibrations (§9), and the foreground
module's `contrast` and `source_density_per_cell` (§8).

**Formalism now implemented.** The baseline-density interferometric noise
model, previously referenced only in prose, is
`analysis.hera_thermal_noise_power_physical()` (§7) — but it is **opt-in**
(`--noise-model physical`), and the default remains the flat scaling estimate.

**Two conventions that differ from a cited source, deliberately.** The
`T_sys` normalisation follows 21cmSense rather than DeBoer et al. (2017)
Table 2 (a factor 2.9 in the sky term, ~8× in `P_N`); and `mode_weighted`
defaults to `False`, omitting La Plante Eqs. 18–19. Both are flagged where
they are used, and both are opt-in to change.
