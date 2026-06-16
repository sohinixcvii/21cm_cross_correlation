# 21 cm – Galaxy Cross-Correlation: Uncertainty Budget

Forecasting the detectability of the 21 cm × galaxy cross-power spectrum during the Epoch of Reionization, following the framework of [La Plante et al. (2023)](https://arxiv.org/abs/2205.09770) and [Davies et al. (2025)](https://arxiv.org/abs/2504.17254).

## Overview

During reionization, neutral hydrogen (H I) emits 21 cm radiation that is **anti-correlated** with the galaxy density field: overdense regions host ionising galaxies and become 21 cm-dark, while underdense regions remain neutral and 21 cm-bright. The cross-power spectrum $P_{21\times\mathrm{gal}}(k_\perp, k_\parallel)$ quantifies this anti-correlation and is a key science target for HERA + Euclid/Roman Space Telescope.

This project contains three complementary notebooks plus a refactored HPC-optimised lightcone pipeline:

1. **`21cm_galaxy_cross_uncertainty.ipynb`** — analytical framework using semi-analytic signal models and the La Plante et al. (2023) variance estimators.
2. **`21cmfast_HERAxEuclid.ipynb`** — end-to-end simulation pipeline using 21cmFASTv4 with the discrete source model (Davies et al. 2025), based on a **coeval** (single-snapshot) simulation at $z = 6.5$.
3. **`21cmfast_HERAxEuclid_lightcone.ipynb`** — same pipeline as (2) but using a **lightcone** simulation spanning a redshift range, producing a non-cubic $(N_\perp \times N_\perp \times N_z)$ volume with continuous redshift evolution along the LOS.

### HPC lightcone pipeline (recommended for cluster use)

The lightcone workflow has been split into three self-contained parts for efficient use on HPC clusters:

| File | Purpose |
|------|---------|
| `run_simulation.py` | **Part 1** — runs the 21cmFAST lightcone, constructs the galaxy field, estimates galaxy bias, applies Kaiser RSD, and saves all outputs to `outputs/lightcone_data.h5` |
| `notebooks/plot_fields.ipynb` | **Part 2** — loads the HDF5 output and visualises the simulated fields (halo catalogue, SFR distributions, lightcone slices, EoR brightness temperature plot) |
| `notebooks/analysis.ipynb` | **Part 3** — loads the HDF5 output and performs all post-simulation calculations: 2D cylindrical power spectra, photo-$z$ damping, foreground wedge excision, SNR estimation, Euclid magnitude/SFR cuts, and effective galaxy bias from the halo catalogue |
| `submit_job.sh` | SLURM batch submission script for Part 1 |

**Workflow:**

```bash
# On the HPC cluster — submit the simulation as a batch job
sbatch submit_job.sh

# Locally (or in a Jupyter session on the cluster) — visualise and analyse
jupyter notebook notebooks/plot_fields.ipynb   # Part 2: field plots
jupyter notebook notebooks/analysis.ipynb      # Part 3: power spectra & SNR
```

The HDF5 file `outputs/lightcone_data.h5` stores all simulation fields (compressed with gzip) and scalar metadata as attributes, so Parts 2 and 3 are completely independent of the simulation run.

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
| $b_\mathrm{gal}$ | 8 | Galaxy bias fallback; Sheth-Tormen HMF integral over the Euclid magnitude range used if `hmf` is installed |
| $t_\mathrm{obs}$ | 1000 h | Integration time |
| Bandwidth | 8 MHz | HERA per-band bandwidth |

---

### 4. HPC pipeline — `run_simulation.py` + `notebooks/plot_fields.ipynb` + `notebooks/analysis.ipynb`

A refactored version of notebook 3 split into three independent parts for cluster use. See the [HPC lightcone pipeline](#hpc-lightcone-pipeline-recommended-for-cluster-use) section above.

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
| Bouwens et al. (2021, ApJ 908, 24) | $2.9\times10^{-4}$ | $-21.03$ | $-2.03$ | Table 5, single-Schechter |
| Finkelstein et al. (2015, ApJ 810, 71) | $7.4\times10^{-4}$ | $-20.81$ | $-1.87$ | Table 3 |

---

**Figure group 6: Stellar Mass – UV Magnitude Relation**

The simulation median and 16–84th percentile scatter band are plotted. Two literature relations are **defined in the code but commented out** and are not rendered in the current notebook:

| Reference | Relation (as coded) | Note |
|-----------|---------------------|------|
| Song et al. (2016, ApJ 825, 5) | $\log_{10}(M_\star/M_\odot) = 8.86 - 0.5\,(M_\mathrm{UV} + 20)$ | Anchored at $M_\mathrm{UV}=-21 \to \log_{10} M_\star = 9.36$; slope from their Figure 5 |
| González et al. (2011, ApJ 736, 133) | $\log_{10}(M_\star/M_\odot) = 9.06 - 0.5\,(M_\mathrm{UV} + 20)$ | $\sim 0.2$ dex higher normalisation from constant-SFH SED assumption |

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

**Survey geometry conversions**

| Function | Description |
|----------|-------------|
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
| `scipy` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `notebooks/analysis.ipynb`, `src/conversions.py` |
| `py21cmfast >= 4.1.1` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `run_simulation.py` |
| `astropy` | `21cmfast_HERAxEuclid.ipynb`, `21cmfast_HERAxEuclid_lightcone.ipynb`, `src/conversions.py` |
| `hmf` | `run_simulation.py`, `notebooks/analysis.ipynb` |
| `h5py` | `run_simulation.py`, `notebooks/plot_fields.ipynb`, `notebooks/analysis.ipynb` |

The analytical notebook (`21cm_galaxy_cross_uncertainty.ipynb`) requires only `numpy` and `matplotlib`; all cosmological calculations use analytic fitting formulae (BBKS transfer function, Carroll et al. growth factor).

## Usage

```bash
# Activate the 21cmfast conda environment first
conda activate 21cmfast

# Analytical framework (no external dependencies beyond numpy/matplotlib)
jupyter notebook 21cm_galaxy_cross_uncertainty.ipynb

# 21cmFAST coeval simulation pipeline (requires py21cmfast >= 4.1.1)
jupyter notebook 21cmfast_HERAxEuclid.ipynb

# 21cmFAST lightcone simulation pipeline — monolithic version
jupyter notebook 21cmfast_HERAxEuclid_lightcone.ipynb

# HPC-optimised lightcone pipeline — recommended for cluster use
sbatch submit_job.sh                                   # Part 1: run simulation
jupyter notebook notebooks/plot_fields.ipynb           # Part 2: field plots
jupyter notebook notebooks/analysis.ipynb              # Part 3: power spectra & SNR
```

Run all cells sequentially. All notebooks are self-contained and generate all figures inline. The simulation notebooks cache 21cmFAST outputs to disk on first run. The HPC pipeline saves simulation outputs to `outputs/lightcone_data.h5` for independent loading by the analysis notebooks.

## References

- **Davies, Mesinger & Murray (2025)** — [arXiv:2504.17254](https://arxiv.org/abs/2504.17254) — 21cmFASTv4 discrete source model
- **Gagnon-Hartman, Davies & Mesinger (2025)** — [arXiv:2502.20447](https://arxiv.org/abs/2502.20447) — galaxy–21 cm cross-correlation detection forecasts
- **La Plante et al. (2023)** — [arXiv:2205.09770](https://arxiv.org/abs/2205.09770) — uncertainty equations, HERA noise model, and foreground wedge prescription
- **Park et al. (2019)**, MNRAS, 484, 933 — 21cmFAST source model parameterisation (stellar-halo relation, SFR timescale `t_STAR × t_H`)
- **Euclid Collaboration (2022)** — [arXiv:2108.01201](https://arxiv.org/abs/2108.01201) — Euclid survey specifications
- **Bouwens et al. (2021)**, ApJ, 908, 24 — UV luminosity function at $z \sim 6$–$8$
- **Finkelstein et al. (2015)**, ApJ, 810, 71 — UV luminosity function at $z \sim 7$–$8$
- **Speagle et al. (2014)**, ApJS, 214, 15 — star-forming main sequence calibration (Eq. 28)
- **Schreiber et al. (2015)**, A&A, 575, A74 — star-forming main sequence at high redshift
- **Song et al. (2016)**, ApJ, 825, 5 — stellar mass–UV magnitude relation at $z \sim 7$
- **Lidz et al. (2009)**, ApJ, 690, 252 — [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — physical signal power spectra
- **Bardeen, Bond, Kaiser & Szalay (1986)**, ApJ, 304, 15 — BBKS transfer function
- **Kaiser (1987)**, MNRAS, 227, 1 — redshift-space distortions
- **DeBoer et al. (2017)**, PASP — [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA instrument specifications
- **Planck Collaboration (2020)**, A&A, 641, A6 — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209) — cosmological parameters
- **Hogg (1999)** — [arXiv:astro-ph/9905116](https://arxiv.org/abs/astro-ph/9905116) — comoving distance and volume formulae
- **Oke & Gunn (1983)**, ApJ, 266, 713 — AB magnitude system
- **Madau & Dickinson (2014)**, ARA&A, 52, 415 — [arXiv:1403.0007](https://arxiv.org/abs/1403.0007) — UV luminosity–SFR calibration ($\kappa_\mathrm{UV}$, Chabrier IMF)
- **Sheth & Tormen (1999)**, MNRAS, 308, 119 — [arXiv:astro-ph/9901122](https://arxiv.org/abs/astro-ph/9901122) — halo mass function and bias formula
- **Murray, Robotham & Power (2013)**, Astron. Comput., 3, 23 — `hmf` halo mass function code
