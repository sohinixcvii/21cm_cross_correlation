# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Roman Space Telescope.

This notebook implements the variance estimators (Equations 15–17) for the cross-spectrum and both auto-spectra, along with physically motivated signal models and instrumental noise, to compute per-mode and total signal-to-noise ratios.

## Equations Implemented

| Equation | Quantity | Source |
|----------|----------|--------|
| Eq. 6 | $T_0(z)$ — brightness-temperature scaling | La Plante et al. (2023) |
| Eq. 11 | $P_N^{21}(k_\perp, k_\parallel)$ — HERA thermal noise | La Plante et al. (2023) |
| Eq. 15 | $\sigma^2_{21,\mathrm{gal}}$ — cross-spectrum variance | La Plante et al. (2023) |
| Eq. 16 | $\sigma^2_{21}$ — 21 cm auto-spectrum variance | La Plante et al. (2023) |
| Eq. 17 | $\sigma^2_\mathrm{gal}$ — galaxy auto-spectrum variance | La Plante et al. (2023) |
| Eqs. 2–5 | Signal power spectra (bias model) | Lidz et al. (2009) |

## Notebook Structure

1. **Imports and setup**
2. **$T_0(z)$** — brightness-temperature scaling factor across EoR redshifts
3. **Photo-$z$ damping** — Gaussian kernel $W(k_\parallel)$ from photometric redshift uncertainty
4. **Fourier grid** — 2D logarithmic $(k_\perp, k_\parallel)$ grid and helper functions
5. **Signal spectra** — CDM matter power spectrum (BBKS transfer function, $\sigma_8$ normalisation) combined with the Lidz et al. (2009) reionization bias model
6. **HERA thermal noise** — physically motivated $P_N^{21}$ from Eq. 11 with baseline density model
7. **Variance estimators** — Eqs. 15–17 for all three spectra
8. **SNR and detectability** — per-mode SNR maps, cumulative SNR vs $k_\mathrm{max}$, and overall detection significance
9. **Photo-$z$ impact** — SNR degradation as a function of $\sigma_z$
10. **Redshift evolution** — uncertainty and SNR trends across $6 \leq z \leq 12$

## Fiducial Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| $z_\mathrm{obs}$ | 8.0 | Observing redshift (mid-reionization) |
| $\bar{x}_\mathrm{HI}$ | 0.5 | Mean neutral fraction |
| $b_\mathrm{gal}$ | 5.0 | Galaxy bias |
| $R_b$ | 10 $h^{-1}$ Mpc | Characteristic bubble radius |
| $N_\mathrm{ant}$ | 350 | HERA Phase II antennas |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| $\sigma_z$ | 0.05 | Fiducial photo-$z$ error |

## Requirements

The notebook uses only standard scientific Python:

```
numpy
matplotlib
```

No external cosmology libraries (e.g. CAMB, astropy) are required — all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).

## Usage

```bash
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb
```

Run all cells sequentially. The notebook is self-contained and generates all figures inline.

## Future Extensions

- Replace BBKS with CAMB/CLASS transfer functions for exact BAO features
- Use [21cmSense](https://github.com/jpober/21cmSense) for the true HERA baseline density
- Apply foreground wedge masking
- Replace the bias model with 21cmFAST simulation outputs
- Extend to a Fisher matrix forecast for $\bar{x}_\mathrm{HI}(z)$ and $R_b(z)$

## References

- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — primary reference; uncertainty equations and HERA noise model
- **Lidz et al. (2009)**, ApJ, 690, 252 — [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — physical signal power spectra
- **Bardeen, Bond, Kaiser & Szalay (1986)**, ApJ, 304, 15 — BBKS transfer function
- **DeBoer et al. (2017)**, PASP — [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA instrument specifications
- **Planck Collaboration (2020)**, A&A, 641, A6 — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — cosmological parameters
