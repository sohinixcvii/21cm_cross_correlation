# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770) and [Davies et al. (2025)](https://arxiv.org/abs/2504.17254).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Euclid/Roman Space Telescope.

This project contains three complementary notebooks:

1. **`21cm_galaxy_cross_uncertainty.ipynb`** — analytical framework using semi-analytic signal models and the La Plante et al. (2023) variance estimators.
2. **`21cmfast_HERAxEuclid.ipynb`** — end-to-end simulation pipeline using 21cmFASTv4 with the discrete source model (Davies et al. 2025), based on a **coeval** (single-snapshot) simulation at $z = 6.5$.
3. **`21cmfast_HERAxEuclid_lightcone.ipynb`** — same pipeline as (2) but using a **lightcone** simulation spanning a redshift range, producing a non-cubic $(N_\perp \times N_\perp \times N_z)$ volume with continuous redshift evolution along the LOS.

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

### 3. `21cmfast_HERAxEuclid_lightcone.ipynb`

Lightcone counterpart to notebook 2. Uses `RectilinearLightconer` + `run_lightcone` to produce a self-consistent lightcone over a redshift range (default $z = 6.5$–$7.5$), then applies the same galaxy field construction, Kaiser RSD, 2D power spectrum calculation, and SNR estimation as the coeval notebook.

**Key differences from the coeval version:**

| Feature | Coeval | Lightcone |
|---|---|---|
| 21cmFAST function | `run_coeval` | `run_lightcone` |
| Box shape | $(N, N, N)$ | $(N, N, N_z)$ |
| Redshift | single snapshot | continuous range |
| LOS cell size | same as transverse | $L_\mathrm{LOS}/N_z$ |
| Field access | `coeval.brightness_temp` | `lightcone.lightcones['brightness_temp']` |

**Structure:**

1. Imports and setup
2. **Configuration cell** — all user-adjustable parameters in one place
3. Derived quantities (LOS geometry, node redshifts)
4. Run 21cmFASTv4 lightcone simulation
5. Construct galaxy density field — from lightcone `halo_sfr` (3a) with a
   synthetic Poisson-sampled fallback for zero-field cases (3b); galaxy bias
   estimated via HMF integration over the Euclid UV magnitude range
6. Kaiser RSD applied in Fourier space
7. Visualise lightcone — transverse (x–y) slice + LOS (x–z) slice with
   **LOS on the x-axis and transverse on the y-axis** (standard convention);
   secondary redshift axis on top via `twiny()`
7b. **Brightness temperature evolution plot** — wide-format (16×3.5") lightcone
   slice styled after Mesinger & Furlanetto (2007), with a custom EoR colourmap
   (dark = ionised, warm/bright = neutral), dual x-axes (comoving distance +
   redshift), and observed frequency range in the title
8. Compute 2D cylindrical power spectra (non-cubic box)
9. Plot power spectra with foreground wedge overlays
10. Photo-$z$ damping and foreground wedge excision
11. Per-mode SNR map and total detection significance
12. Summary — coeval vs. lightcone comparison table and next-step recommendations

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
| `HII_DIM` | 128 | Simulation grid cells per side |
| `BOX_LEN` | 256 Mpc | Transverse box size |
| $z_\mathrm{min}$, $z_\mathrm{max}$ | 6.5, 7.5 | Lightcone redshift range |
| $z_\mathrm{obs}$ | 7.0 | Reference redshift (midpoint of lightcone) |
| $M_\mathrm{UV}$ limit | $< -18$ | Euclid galaxy selection |
| $\sigma_z$ | 0.059 | Euclid photo-$z$ uncertainty |
| $\bar{n}_\mathrm{gal}$ | $3\times10^{-3}\ h^3\ \mathrm{Mpc}^{-3}$ | Mean galaxy number density |
| $b_\mathrm{gal}$ | 8 | Galaxy bias (estimated via HMF integration) |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| Bandwidth | 8 MHz | HERA per-band bandwidth |

---

## 21cmFASTv4 `HaloBox` API Notes

In 21cmFAST v4.1+, `coeval.halobox` is a `HaloBox` object whose arrays are accessed via `.get('<field_name>')`. In lightcone runs, the same fields are stored in `lightcone.lightcones['<field_name>']`. Available fields include:

| Field | Description |
|-------|-------------|
| `halo_sfr` | Total SFR per cell, summed over all halos [internal units] |
| `n_ion` | Number of ionising photons per cell |

Individual halo positions and UV magnitudes are not exposed by this API. A strict per-halo $M_\mathrm{UV}$ cut requires post-processing the raw halo catalogue from a lightcone run.

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

---

## Source Modules

### `src/conversions.py`

Cosmological conversion utilities for high-redshift galaxy surveys. Import
individual functions as needed:

```python
from src.conversions import (
    Muv_to_Luv,
    Luv_to_Muv,
    survey_area_from_volume,
    area_deg2_to_steradians,
    volume_from_area,
)
```

| Function | Description |
|----------|-------------|
| `Muv_to_Luv(Muv)` | Absolute UV AB magnitude → monochromatic luminosity [erg s⁻¹ Hz⁻¹] |
| `Luv_to_Muv(Luv)` | Monochromatic luminosity [erg s⁻¹ Hz⁻¹] → AB magnitude |
| `survey_area_from_volume(volume_mpc3, z_min, z_max, cosmo=None)` | Comoving volume [Mpc³] → survey area [deg²] |
| `area_deg2_to_steradians(area_deg2)` | Survey area [deg²] → [sr] |
| `volume_from_area(area_deg2, z_min, z_max, cosmo=None, n_z=1000)` | Survey area [deg²] → comoving volume [Mpc³] |

All functions accept scalar or array inputs. Volume–area conversions use
Simpson integration of the differential comoving volume
$\mathrm{d}V/\mathrm{d}z\,\mathrm{d}\Omega$ (Hogg 1999) and default to the
Planck18 cosmology; pass a custom `astropy.cosmology` object via `cosmo` to
override.

**Dependencies:** `numpy`, `astropy`, `scipy`

---

## Requirements

| Package | Used in |
|---------|---------|
| `numpy` | All notebooks, `src/conversions.py` |
| `matplotlib` | All notebooks |
| `scipy` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `src/conversions.py` |
| `py21cmfast >= 4.1.1` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb` |
| `astropy` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `src/conversions.py` |

The analytical notebook (`21cm_galaxy_cross_uncertainty.ipynb`) requires only `numpy` and `matplotlib`; all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).

## Usage

```bash
# Analytical framework (no external dependencies)
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb

# 21cmFAST coeval simulation pipeline (requires py21cmfast >= 4.1.1)
jupyter notebook 21cmfast_HERAxEuclid.ipynb

# 21cmFAST lightcone simulation pipeline (requires py21cmfast >= 4.1.1)
jupyter notebook 21cmfast_HERAxEuclid_lightcone.ipynb
```

Run all cells sequentially. All notebooks are self-contained and generate all figures inline. The simulation notebooks cache 21cmFAST outputs to disk on first run.

## References

- **Davies, Mesinger & Murray (2025)** — [arXiv:2504.17254](https://arxiv.org/abs/2504.17254) — 21cmFASTv4 discrete source model
- **Gagnon-Hartman, Davies & Mesinger (2025)** — [arXiv:2502.20447](https://arxiv.org/abs/2502.20447) — galaxy–21 cm cross-correlation detection forecasts
- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — uncertainty equations, HERA noise model, and foreground wedge prescription
- **Euclid Collaboration (2022)** — [arXiv:2108.01201](https://arxiv.org/abs/2108.01201) — Euclid survey specifications
- **Lidz et al. (2009)**, ApJ, 690, 252 — [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — physical signal power spectra
- **Bardeen, Bond, Kaiser & Szalay (1986)**, ApJ, 304, 15 — BBKS transfer function
- **Kaiser (1987)**, MNRAS, 227, 1 — redshift-space distortions
- **DeBoer et al. (2017)**, PASP — [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA instrument specifications
- **Planck Collaboration (2020)**, A&A, 641, A6 — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — cosmological parameters
- **Hogg (1999)** — [arXiv:astro-ph/9905116](https://arxiv.org/abs/astro-ph/9905116) — comoving distance and volume formulae
- **Oke & Gunn (1983)**, ApJ, 266, 713 — AB magnitude system
- **Madau & Dickinson (2014)**, ARA&A, 52, 415 — [arXiv:1403.0007](https://arxiv.org/abs/1403.0007) — UV luminosity and star formation rate density
