# Low SFR in 21cmFAST Halo Catalogue — Diagnosis and Fix

## Summary

The star-forming main sequence plot in `notebooks/plot_fields.ipynb` (Section 7)
showed SFR values ~7.5 dex below the Speagle+14 and Schreiber+15 literature
calibrations.  Investigation revealed two separate issues:

1. **Unit bug (code fix required):** `perturbed_halos.sfr.value` in py21cmfast v4
   returns the SFR in **M☉ s⁻¹** (the internal C-code unit), not M☉ yr⁻¹ as
   stated in the documentation.  No conversion was applied before saving to HDF5,
   so all stored SFR values were a factor of `3.156 × 10⁷` (seconds per year) too
   small.

2. **Model physics (expected, not a bug):** Even with correct units, the 21cmFAST
   "simple" template uses a star-formation timescale `t_sf = t_STAR × t_H(z)` with
   `t_STAR = 0.5` and `t_H(z=7) ≈ 1.14 Gyr`, giving `t_sf ≈ 570 Myr`.  This is
   ~5.7× longer than the constant 100 Myr timescale assumed by Schreiber+15, so the
   21cmFAST SFR is intrinsically ~0.76 dex below the observed main sequence.  This
   is a deliberate calibration choice for reionisation modelling (Park+2018), not an
   error.

---

## 1. Unit Bug

### Root Cause

In `py21cmfast/src/scaling_relations.c`, the per-halo SFR is computed as:

```c
sfr_mean = stellar_mass / (consts->t_star * consts->t_h);
```

where `t_h = t_hubble(redshift)` is evaluated in **seconds** (the internal
time unit of the C code).  `stellar_mass` is in M☉.  Therefore `sfr_mean` has
units of **M☉ s⁻¹**.

The Python wrapper stores this value in the `sfr` attribute of the perturbed halo
catalogue.  The `.value` property returns the raw numerical value without any unit
conversion, so it is still in M☉ s⁻¹.

The documentation (and `docs/halo_catalogue_reference.md`) incorrectly states the
unit as M☉ yr⁻¹.

### Diagnosis

Numerical verification using the simulation data (`outputs/lightcone_data.h5`,
z_obs = 7.0, 114 M halos):

| Quantity | Raw `.value` (M☉ yr⁻¹ assumed) | Corrected (× 3.156 × 10⁷) |
|---|---|---|
| SFR max | 1.546 × 10⁻⁵ | 487.8 M☉ yr⁻¹ |
| M★ max  | 2.794 × 10¹¹ M☉ | — |
| t_eff = M★/SFR | **1.81 × 10¹⁶ yr** (impossible) | **5.73 × 10⁸ yr ≈ 573 Myr** ✓ |
| Offset from model | −7.77 dex | ≈ 0 dex ✓ |

The corrected `t_eff = 573 Myr` matches `t_STAR × t_H(z=7) = 0.5 × 1.141 Gyr =
570 Myr` to within 0.6%, confirming the unit conversion is the only issue.

### Fix Applied

**`run_simulation.py`, section 3a (halo catalogue extraction):**

```python
# Before (wrong — returns M_sun s^-1):
sfr_cat = np.asarray(perturbed_halos.sfr.value)

# After (correct — converts to M_sun yr^-1):
_SEC_PER_YR = 365.25 * 24 * 3600          # 3.15576e7 s yr^-1
sfr_cat = np.asarray(perturbed_halos.sfr.value) * _SEC_PER_YR
```

`halo_masses` and `stellar_masses` do **not** require conversion: their `.value`
attributes are in M☉ as expected (verified by checking M★/M_h ratios against
the stellar-halo model).

---

## 2. Model Physics Offset (Expected)

### Why 21cmFAST SFR ≠ Observed SFR

Even after correcting the unit bug, the 21cmFAST SFR lies **~0.76 dex below**
the Schreiber+15 and Speagle+14 calibrations.  This is expected.

#### 21cmFAST source model (Park+2018 Eq. 3)

```
SFR = M★ / (t_STAR × t_H(z))
```

Default parameter: `t_STAR = 0.5`

At z = 7.0 (Planck 2018 cosmology):

| Quantity | Value |
|---|---|
| Hubble time t_H = 1/H(z) | 1.141 Gyr |
| Effective timescale t_sf = 0.5 × t_H | **570 Myr** |
| Expected sSFR = 1/t_sf | **1.75 Gyr⁻¹** |

#### Observational calibrations at z ~ 7

| Reference | sSFR assumed | t_sf implied |
|---|---|---|
| Schreiber+15 high-z | 10 Gyr⁻¹ | 100 Myr |
| Speagle+14 (Eq. 28) | ~8 Gyr⁻¹ at z=7 | ~125 Myr |
| **21cmFAST "simple"** | **1.75 Gyr⁻¹** | **570 Myr** |

The expected offset: Δ log₁₀(SFR) = log₁₀(1.75) − log₁₀(10) ≈ **−0.76 dex**.

#### Physical reason

Speagle+14 and Schreiber+15 calibrate the **current (instantaneous) SFR** from
UV or Hα luminosity, which traces the last ~10–100 Myr of star formation.  The
21cmFAST source model uses a parametric sSFR tied to the **Hubble time** at each
redshift, scaled by `t_STAR`.  This is calibrated to reproduce the reionisation
history (UV emissivity, ionising photon rate), not to match the normalisation of
the observed star-forming main sequence.

The two sSFR definitions are physically distinct; the offset is *not* an error in
the 21cmFAST model.

### How the Notebook Plot Handles This

Section 7 of `notebooks/plot_fields.ipynb` now plots three reference lines:

1. **Green solid — "21cmFAST model"**: `log₁₀(SFR) = log₁₀(M★) − log₁₀(t_STAR × t_H)`
   computed analytically from `t_STAR = 0.5` and `t_H` from Planck18.  The
   simulation median should lie on this line (after the unit fix).

2. **Blue dashed — Speagle+14**: observational UV+IR main sequence, using the
   cosmic age `t_age` at z_obs.

3. **Red dash-dot — Schreiber+15**: constant sSFR = 10 Gyr⁻¹.

An annotation on the plot shows the expected offset in dex.

---

## 3. Changes Made

| File | Change |
|---|---|
| `run_simulation.py` | Multiply `sfr.value` by `_SEC_PER_YR = 3.156 × 10⁷` when extracting halo SFR |
| `notebooks/plot_fields.ipynb` cell `1p261xgkd4b` | Section 7 replaced duplicate UVLF code with star-forming main sequence; added 21cmFAST model line and offset annotation; added sSFR verification print |

---

## 4. Verification

Diagnostic run after applying the fix (`run_simulation.py`, py21cmfast 4.1.1, z_obs = 7.0,
256 Mpc box, 114,291,212 halos):

| Quantity | Before fix | After fix | Expected |
|---|---|---|---|
| SFR min | 6.8 × 10⁻¹⁸ M☉ yr⁻¹ | 2.2 × 10⁻¹⁰ M☉ yr⁻¹ | — |
| SFR max | 1.5 × 10⁻⁵ M☉ yr⁻¹ | **488 M☉ yr⁻¹** | physically reasonable |
| sSFR (all halos, median) | ~0 Gyr⁻¹ | 0.93 Gyr⁻¹ | 1.75 Gyr⁻¹ * |
| sSFR (top-1000 M★, median) | ~0 Gyr⁻¹ | **1.58 Gyr⁻¹** | **1.75 Gyr⁻¹** |
| M★/SFR (top-1000, median) | — | **634 Myr** | **571 Myr** |
| M_UV range | +30 to +24 (inverted, nonsensical) | **−25 to +6** | sensible |
| N(M_UV ≤ −18) Euclid-bright | ~0 | **49,621** | realistic |
| Offset from model (log space) | −7.77 dex | ~−0.07 dex | 0 dex |

\* The all-halo median sSFR (0.93 Gyr⁻¹) is lower than the model prediction (1.75 Gyr⁻¹)
because the bulk of the 114 M halos are near the minimum mass (10⁸ M☉) where the
M_TURN = 5 × 10⁸ M☉ exponential suppression `exp(−M_TURN/M_h)` actively reduces SFR
below the unsuppressed relation.  The top-1000 massive halos (where suppression is
negligible) agree with the model to within ~10%, consistent with the log-normal
scatter (`σ_SFR = 0.19 dex`, `σ_STAR = 0.25 dex`) applied by 21cmFAST.

The UV magnitude range (−25 to +6) and the Euclid-bright count (49,621 galaxies)
are physically sensible after the fix.

---

## 5. References

- Park, J. et al. (2018, MNRAS 484, 933) — 21cmFAST source model parameterisation
- Schreiber, C. et al. (2015, A&A 575, A74) — star-forming main sequence at high z
- Speagle, J. S. et al. (2014, ApJS 214, 15) — Eq. 28, main sequence calibration
- Madau, P. & Dickinson, M. (2014, ARA&A 52, 415) — UV–SFR calibration
