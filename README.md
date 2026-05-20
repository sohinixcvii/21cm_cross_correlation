# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770) and [Davies et al. (2025)](https://arxiv.org/abs/2504.17254).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Euclid/Roman Space Telescope.

This project contains two complementary notebooks:

1. **`21cm_galaxy_cross_uncertainty.ipynb`** — analytical framework using semi-analytic signal models and the La Plante et al. (2023) variance estimators.
2. **`21cmfast_HERAxEuclid.ipynb`** — end-to-end simulation pipeline using 21cmFASTv4 with the discrete source model (Davies et al. 2025), targeting the HERA × Euclid science case.

## Notebooks

### 1. `21cm_galaxy_cross_uncertainty.ipynb`

Implements the variance estimators (Equations 15–17) for the cross-spectrum and both auto-spectra, along with physically motivated signal models and instrumental noise, to compute per-mode and total signal-to-noise ratios.

**Structure:**

1. Imports and setup
2. $T_0(z)$ — brightness-temperature scaling factor across EoR redshifts
3. Photo-$z$ damping — Gaussian kernel $W(k_\parallel)$ from photometric redshift uncertainty
4. Fourier grid — 2D logarithmic $(k_\perp, k_\parallel)$ grid and helper functions
5. Signal spectra — CDM matter power spectrum (BBKS transfer function, $\sigma_8$ normalisation) combined with the Lidz et al. (2009) reionization bias model
6. HERA thermal noise — physically motivated $P_N^{21}$ from Eq. 11 with baseline density model
7. Variance estimators — Eqs. 15–17 for all three spectra
8. SNR and detectability — per-mode SNR maps, cumulative SNR vs $k_\mathrm{max}$, and overall detection significance
9. Photo-$z$ impact — SNR degradation as a function of $\sigma_z$
10. Redshift evolution — uncertainty and SNR trends across $6 \leq z \leq 12$

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

### 2. `21cmfast_HERAxEuclid.ipynb`

End-to-end HERA × Euclid cross-correlation workflow using 21cmFASTv4 with the discrete source (CHMF-SAMPLER) model. Runs a coeval simulation at $z = 6.5$, constructs the galaxy density field from the `HaloBox` SFR proxy, and computes 2D cylindrical power spectra with foreground wedge excision and photo-$z$ damping.

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
| $\sigma_z$ | 0.059 | Euclid photo-$z$ uncertainty |
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

## 21cmFASTv4 `HaloBox` API Notes

In 21cmFAST v4.1+, `coeval.halobox` is a `HaloBox` object whose arrays are accessed via `.get('<field_name>')`. Available fields include:

| Field | Description |
|-------|-------------|
| `halo_sfr` | Total SFR per cell, summed over all halos [internal units] |
| `n_ion` | Number of ionising photons per cell |

Individual halo positions and UV magnitudes are not exposed by this API. For per-halo catalogues, use the raw halo output from a lightcone run.

The galaxy overdensity is constructed as $\delta_\mathrm{gal} = \mathrm{SFR}/\langle\mathrm{SFR}\rangle - 1$, which correctly traces the galaxy distribution and produces the expected large-scale anti-correlation with the 21 cm field (negative cross-spectrum on large scales).

---

## Requirements

| Package | Used in |
|---------|---------|
| `numpy` | Both notebooks |
| `matplotlib` | Both notebooks |
| `scipy` | `21cmfast_HERAxEuclid.ipynb` |
| `py21cmfast >= 4.1.1` | `21cmfast_HERAxEuclid.ipynb` |
| `astropy` | `21cmfast_HERAxEuclid.ipynb` |

The analytical notebook (`21cm_galaxy_cross_uncertainty.ipynb`) requires only `numpy` and `matplotlib`; all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).

## Usage

```bash
# Analytical framework (no external dependencies)
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb

# 21cmFAST simulation pipeline (requires py21cmfast >= 4.1.1)
jupyter notebook 21cmfast_HERAxEuclid.ipynb
```

Run all cells sequentially. Both notebooks are self-contained and generate all figures inline. The simulation notebook caches 21cmFAST outputs to disk on first run.

## References

- **Davies, Mesinger & Murray (2025)** — [arXiv:2504.17254](https://arxiv.org/abs/2504.17254) — 21cmFASTv4 discrete source model
- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — uncertainty equations, HERA noise model, and foreground wedge prescription
- **Lidz et al. (2009)**, ApJ, 690, 252 — [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — physical signal power spectra
- **Bardeen, Bond, Kaiser & Szalay (1986)**, ApJ, 304, 15 — BBKS transfer function
- **Kaiser (1987)**, MNRAS, 227, 1 — redshift-space distortions
- **DeBoer et al. (2017)**, PASP — [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA instrument specifications
- **Planck Collaboration (2020)**, A&A, 641, A6 — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — cosmological parameters
