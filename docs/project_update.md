# 21cm × Galaxy Cross-Correlation — Project Status Update

**Date:** 2026-06-15

---

## 1. Simulation Parameters

| Parameter | Value |
|---|---|
| Code | py21cmfast 4.1.1 ("simple" template) |
| Grid size | HII_DIM = 128, DIM = 384 |
| Box side length | BOX_LEN = 256.0 comoving Mpc |
| Cell size | 2.0 Mpc |
| Redshift range (lightcone) | z = 6.995 → 7.005 |
| Reference redshift | z_obs = 7.0 |
| LOS slices | N_z = 100 (minimum enforced) |
| LOS extent | ~200 Mpc |
| Random seed | 42 |
| Mean neutral fraction | ⟨x_HI⟩ = 0.176 |
| Brightness temp range | 0.0 – 22.2 mK |

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

## 2. Halo Catalogue

Extracted at z_obs = 7.0 using `determine_halo_catalog` + `perturb_halo_catalog`.

| Quantity | Value |
|---|---|
| Total halos | 114,291,212 |
| Halo mass range | 1.0 × 10⁸ – 1.77 × 10¹² M☉ |
| Stellar mass range | 29 – 2.79 × 10¹¹ M☉ |
| SFR range | 2.2 × 10⁻¹⁰ – 488 M☉ yr⁻¹ |
| SFR median | 1.17 × 10⁻⁵ M☉ yr⁻¹ |

---

## 3. Euclid Survey Cuts

| Parameter | Value |
|---|---|
| UV magnitude limit | -22.66 <= M_UV ≤ −18 (AB mag) |
| Photo-z uncertainty | σ_z = 0.059 |
| Target galaxy number density | n̄ = 3 × 10⁻³ h³ Mpc⁻³ |
| Euclid-bright galaxies (simulation) | **49,621** |
| Fraction Euclid-detected | ~0.04% of all halos |

The selection is applied as M_UV ≤ −18, where M_UV is computed from
SFR via the Madau & Dickinson (2014) calibration:

```
L_UV [erg/s/Hz] = SFR [M☉/yr] / κ_UV,   κ_UV = 1.15 × 10⁻²⁸
M_UV = −2.5 log₁₀(L_UV) + 51.60
```

---

## 4. Galaxy Bias

Estimated via Sheth-Tormen bias integrated over the luminosity-selected
halo mass function (HMF from `hmf.MassFunction`):

```
b_g = ∫ (dn/dlogM) b_h(M) φ(M) dlogM  /  ∫ (dn/dlogM) φ(M) dlogM
```

where φ(M) = 1 if L_UV(M) ∈ [L_UV(−22), L_UV(−18)], else 0.

| Quantity | Value |
|---|---|
| Galaxy bias b_g | 33.39 |
| Selected halo mass range | 9.1 × 10⁹ – 10¹¹ M☉ |
| n_gal | 2.37 × 10⁻² Mpc⁻³ |
| Growth rate f = Ω_m(z)^0.55 | 0.996 |
| β = f/b | 0.030 |

> **Note:** b_g = 33.39 is anomalously high. This likely results from the
> Euclid bright limit (M_UV = −22) being too restrictive relative to the
> 21cmFAST UVLF at this box size, selecting only very rare, high-bias halos.
> The bias calculation also uses a manual SFR timescale (100 Myr) that differs
> from the 21cmFAST internal model (t_STAR × t_H ≈ 570 Myr). This should be
> revisited with a consistent SFR model.

---

## 5. Predicted Observable Plots (from `notebooks/plot_fields.ipynb`)

### 5.1 UV Luminosity Function (Section 5)

- **Data**: All 21cmFAST halos with SFR > 0 (no magnitude pre-filtering)
- **x-axis**: M_UV [AB mag], range −25.5 to −11.0
- **y-axis**: Φ(M_UV) [Mpc⁻³ mag⁻¹], log scale 10⁻⁸ – 0.2
- **Literature**: Schechter fits from Bouwens+21 (z~7, φ* = 2.9×10⁻⁴,
  M* = −21.03, α = −2.03) and Finkelstein+15 (z~7)
- **Reference line**: Euclid limit at M_UV = −18 (vertical dashed)
- **Annotation**: N_Euclid-bright = 49,621 / N_total = 114,291,212

The 21cmFAST UVLF tracks the Schechter function at the bright end but
turns over at the faint end where M_TURN exponential suppression acts.
Bright-end statistics are limited by the 256 Mpc box volume.

### 5.2 Stellar Mass – UV Magnitude Relation (Section 6)

- **Data**: All halos with M★ > 0 and SFR > 0
- **x-axis**: M_UV [AB mag]
- **y-axis**: log₁₀(M★/M☉), binned median + 16–84% scatter
- **Literature**: Song+16 (`log M★ = 8.86 − 0.5(M_UV + 20)`) and
  González+11 (`log M★ = 9.06 − 0.5(M_UV + 20)`) — currently commented
  out in the notebook, can be re-enabled

### 5.3 Star-forming Main Sequence (Section 7)

- **Data**: All halos with M★ > 0 and SFR > 0 (no Euclid magnitude cut)
- **x-axis**: log₁₀(M★/M☉), range 4.5 – 11.5
- **y-axis**: log₁₀(SFR / M☉ yr⁻¹), binned median + 16–84% scatter
- **Three reference lines**:

| Line | Colour | Formula |
|---|---|---|
| 21cmFAST model | Green solid | log₁₀(SFR) = log₁₀(M★) − log₁₀(t_STAR × t_H) |
| Speagle+14 | Blue dashed | (0.84 − 0.026 t_age) log₁₀(M★) − (6.51 − 0.11 t_age) |
| Schreiber+15 | Red dash-dot | log₁₀(SFR) = log₁₀(M★) − 8  (sSFR = 10 Gyr⁻¹) |

Expected result: simulation data lie along the green 21cmFAST model line,
~0.76 dex below Schreiber+15. This offset is physical (different SFR
timescales), not a modelling error. Top-1000 most massive halos verify
sSFR ≈ 1.58 Gyr⁻¹ vs model prediction 1.75 Gyr⁻¹ (~10% consistent with
log-normal SFR scatter of σ = 0.19 dex).

---

## 6. SFR Model Comparison at z = 7

| Reference | sSFR | Timescale | Basis |
|---|---|---|---|
| Schreiber+15 high-z | 10 Gyr⁻¹ | 100 Myr | UV/Hα observations |
| Speagle+14 | ~8 Gyr⁻¹ | ~125 Myr | UV+IR compilation |
| 21cmFAST "simple" | **1.75 Gyr⁻¹** | **570 Myr** | t_STAR × t_H(z), Park+2018 |

The 21cmFAST model is calibrated to reproduce the reionisation history
(UV emissivity and ionising photon rate), not to match the normalisation of
the observed star-forming main sequence. The ~0.76 dex offset is expected
and is documented in `docs/Low_SFR_fix.md`.

---

## 7. Kaiser RSD

Applied in Fourier space to the galaxy overdensity field derived from the
`halo_sfr` lightcone:

```
δ_gal^(s)(k) = (1 + β μ²) δ_gal(k),   β = f / b_g = 0.030
```

| Quantity | Value |
|---|---|
| Growth rate f = Ω_m(z)^0.55 | 0.996 |
| β = f/b | 0.030 |
| Max Kaiser boost (μ = 1) | 1.030× |
| δ_gal std before RSD | 4.535 |
| δ_gal std after RSD | 4.537 |

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
| `halo_catalog/halo_masses` | (114,291,212,) | M☉ |
| `halo_catalog/halo_coords` | (114,291,212, 3) | Mpc |
| `halo_catalog/stellar_masses` | (114,291,212,) | M☉ |
| `halo_catalog/sfr` | (114,291,212,) | M☉ yr⁻¹ (unit-corrected) |

---

## 9. Pipeline Status

| Step | File | Status |
|---|---|---|
| Part 1: Lightcone simulation + halo catalogue | `run_simulation.py` | Complete (SFR unit fix applied) |
| Part 2: Field visualisation | `notebooks/plot_fields.ipynb` | Updated (Sections 5–7 corrected) |
| Part 3: Power spectra & SNR | `notebooks/analysis.ipynb` | Not yet run |

---

## 10. Known Issues and Open Questions

1. **Galaxy bias b_g = 33.39**: Anomalously high. The Euclid bright cut
   (M_UV = −22 to −18) selects very rare, high-bias halos in a 256 Mpc box.
   The bias model uses a manual SFR timescale (100 Myr) inconsistent with the
   21cmFAST internal model (570 Myr). Recommend revisiting with consistent
   F_STAR10 and t_STAR × t_H parameterisation.

2. **Literature lines in Section 6 commented out**: Song+16 and González+11
   reference lines for M★–M_UV are disabled. Re-enable if comparison needed.

3. **Narrow lightcone redshift range**: z = 6.995–7.005 (Δz = 0.01) is very
   narrow for most science cases. The 200 Mpc LOS output is driven by the
   `minimum_los_slices = 100` override, not by physical redshift extent.
   Consider widening z range for a more representative lightcone.

---

## 11. Key Files

| File | Purpose |
|---|---|
| `run_simulation.py` | Part 1: 21cmFAST lightcone + halo catalogue + Kaiser RSD |
| `notebooks/plot_fields.ipynb` | Part 2: Visualisation (UVLF, M★–MUV, SFR main sequence, fields) |
| `notebooks/analysis.ipynb` | Part 3: Power spectra and SNR estimation |
| `src/conversions.py` | SFR↔L_UV↔M_UV conversions, Sheth-Tormen bias, survey geometry |
| `src/FOV_to_cMpc.py` | Angular survey area to comoving Mpc conversion |
| `docs/Low_SFR_fix.md` | Diagnosis and fix for SFR unit bug |
| `docs/Galaxy_bias_formalism.md` | Galaxy bias theory and implementation |
| `docs/halo_catalogue_reference.md` | 21cmFAST v4 halo catalogue API reference |

---

## 12. References

- Park, J. et al. (2018, MNRAS 484, 933) — 21cmFAST source model parameterisation
- Bouwens, R. J. et al. (2021, ApJ 908, 24) — UVLF at z~6–8
- Finkelstein, S. L. et al. (2015, ApJ 810, 71) — UVLF at z~7–8
- Speagle, J. S. et al. (2014, ApJS 214, 15) — Star-forming main sequence
- Schreiber, C. et al. (2015, A&A 575, A74) — Main sequence at high-z
- Madau, P. & Dickinson, M. (2014, ARA&A 52, 415) — UV–SFR calibration
- Song, M. et al. (2016, ApJ 825, 5) — M★–M_UV relation at z~7
- González, V. et al. (2011, ApJ 735, L34) — M★–M_UV at z~7
