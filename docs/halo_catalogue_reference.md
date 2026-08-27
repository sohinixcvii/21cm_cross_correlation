# 21cmFAST v4 Halo Catalogue Reference

**Project:** HERA × Euclid Cross-Correlation Pipeline  
**Last Updated:** 2026-06-11  
**Verified Against:** py21cmfast v4.x (Davies et al. 2025)

---

# Purpose

This document records the physical meaning, units, provenance, and usage of the halo catalogue quantities produced by 21cmFAST v4.

The goal is to provide a single authoritative reference for future analyses involving:

- Galaxy selection
- UV luminosity functions
- Euclid mock surveys
- Galaxy bias calculations
- 21 cm–galaxy cross-correlations

---

# Which Halo Catalogue?

21cmFAST v4 can generate several halo-related data products.

For cross-correlation work we primarily use:

```python
PerturbedHaloCatalog
```

This catalogue contains halo properties after the density field has been evolved and halos have been displaced into their final large-scale structure positions.

---

# Accessing the Catalogue

Typical workflow:

```python
halo_catalog = p21c.determine_halo_catalog(...)
```

or

```python
halo_catalog = p21c.perturb_halo_catalog(...)
```

Available fields:

```python
[
    "fesc_sfr",
    "halo_coords",
    "halo_masses",
    "n_halos",
    "sfr",
    "sfr_mini",
    "stellar_masses",
]
```

---

# Catalogue Fields

---

## halo_coords

### Description

Comoving positions of halos within the simulation volume.

### Units

```text
Mpc (comoving)
```

### Shape

```python
(N_halos, 3)
```

### Components

```python
x, y, z
```

Each coordinate gives the comoving position of the halo inside the simulation box.

### Notes

The values are already physical coordinates within the box.

They are **not grid indices** and should **not** be multiplied by:

```python
BOX_LEN / HII_DIM
```

again.

Example:

```text
254.6
114.0
191.6
```

in a

```text
BOX_LEN = 256 Mpc
```

simulation corresponds to a halo near the box boundary.

### Used For

- Galaxy maps
- Spatial clustering
- Cross-power spectra
- Halo visualisation

---

## halo_masses

### Description

Dark matter halo masses.

### Units

```text
Msun
```

(solar masses)

### Shape

```python
(N_halos,)
```

### Typical Values

```text
10^8 – 10^12 Msun
```

depending on redshift and source model.

### Used For

- Halo mass functions
- Halo occupation models
- Stellar-mass assignment
- Galaxy bias calculations

---

## stellar_masses

### Description

Stellar mass associated with each halo.

### Units

```text
Msun
```

### Shape

```python
(N_halos,)
```

### Physical Meaning

21cmFAST computes the stellar mass using the galaxy formation model described in:

- Park et al. (2019)
- Davies et al. (2025)

This quantity represents the total stellar mass contained within the halo.

### Used For

- Stellar mass functions
- Galaxy selection
- UV luminosity calculations

---

## sfr

### Description

Population II star formation rate.

### Units

```text
Msun yr^-1
```

### Shape

```python
(N_halos,)
```

### Physical Meaning

The rate at which stars are currently forming inside each halo.

This is the primary quantity used for UV luminosity calculations.

### Used For

- UV luminosity functions
- Galaxy catalogues
- Euclid mock observations
- Ionising emissivity calculations

---

## sfr_mini

### Description

Minihalo (Population III) star formation rate.

### Units

```text
Msun yr^-1
```

### Shape

```python
(N_halos,)
```

### Physical Meaning

Star formation occurring in molecular-cooling minihalos.

These halos are below the atomic cooling threshold and are relevant primarily during Cosmic Dawn.

### Notes

Many simulations may contain:

```text
sfr_mini = 0
```

for most halos.

### Used For

- Population III studies
- Cosmic Dawn modelling
- Early X-ray source calculations

---

## fesc_sfr

### Description

Escape-fraction weighted star formation rate.

### Definition

```text
fesc_sfr = fesc × SFR
```

where:

```text
fesc
```

is the ionising-photon escape fraction.

### Units

```text
Msun yr^-1
```

### Shape

```python
(N_halos,)
```

### Physical Meaning

Ionising emissivity scales approximately with:

```text
fesc × SFR
```

so 21cmFAST stores this quantity directly.

### Used For

- Reionisation calculations
- Ionising photon budgets
- Source emissivity calculations

---

## n_halos

### Description

Total number of halos in the catalogue.

### Units

```text
dimensionless
```

### Type

```python
int
```

### Example

```python
print(catalog.n_halos)
```

---

# Recommended Quantities for Euclid Mock Catalogues

For Euclid cross-correlation studies the most important quantities are:

```python
halo_coords
halo_masses
stellar_masses
sfr
```

The typical workflow is:

```text
halo mass
      ↓
stellar mass
      ↓
star formation rate
      ↓
UV luminosity
      ↓
absolute UV magnitude
      ↓
Euclid magnitude cut
      ↓
galaxy catalogue
```

---

# UV Luminosity Conversion

A commonly used conversion is:

L_UV ≈ 1.15 × 10^28 × SFR

> **This page describes the historical relation, which this pipeline no longer
> uses.** The adopted value is now `kappa_UV = 2.7e-29` (Fisher et al. 2026,
> arXiv:2511.10741, Eq. 12), i.e. `L_UV ≈ 3.7 × 10^28 × SFR`. Critically, the
> assumptions listed below — **"Continuous star formation"** and a
> **"Salpeter-like IMF"** — are exactly what the new value does *not* assume:
> it is calibrated on **rising** SFHs to recover `SFR_100Myr`, on a Chabrier
> setup. This conflict is recorded as Caveat 1 in `NUMBERS_AND_SOURCES.md` §2
> and is **unresolved**, pending human review.

where:

```text
L_UV
```

has units:

```text
erg s^-1 Hz^-1
```

and

```text
SFR
```

has units:

```text
Msun yr^-1
```

This relation assumes:

- Continuous star formation
- Salpeter-like IMF
- Solar/sub-solar metallicity

The exact normalisation depends on SPS assumptions.

---

# Units of all parameters (tabled)

| Quantity | Units | Description |
|-----------|--------|-------------|
| `halo_coords` | Mpc (comoving) | Three-dimensional halo positions `(x, y, z)` within the simulation box. |
| `halo_masses` | Msun | Dark matter halo masses. |
| `stellar_masses` | Msun | Total stellar mass assigned to each halo by the source model. |
| `sfr` | Msun yr^-1 | Population II star formation rate. Used for UV luminosity calculations. |
| `sfr_mini` | Msun yr^-1 | Population III (minihalo) star formation rate. |
| `fesc_sfr` | Msun yr^-1 | Escape-fraction-weighted star formation rate, `f_esc × SFR`. |
| `n_halos` | Dimensionless | Total number of halos in the catalogue. |

---

# Caveats

## Not Every Halo Hosts an Observable Galaxy

The halo catalogue contains all simulated halos.

Only a subset will satisfy:

```text
M_UV < M_lim
```

for a given survey.

---

## Halo Mass ≠ Stellar Mass

Do not use:

```text
Mhalo = Mstar
```

These are distinct quantities.

---

## Halo Positions Are Already Physical

Do not rescale:

```python
halo_coords
```

by:

```python
BOX_LEN / HII_DIM
```

unless you have explicitly confirmed they are stored as cell indices.

For 21cmFAST v4 perturbed halo catalogues they are typically stored as comoving coordinates.

---

# References

## Primary 21cmFAST v4 Paper

Davies, J. E. et al. (2025)

**21cmFAST v4: A Python-C Framework for Simulating the Cosmic 21-cm Signal**

Primary reference for:

- Halo catalogue generation
- Source models
- Galaxy properties
- SFR prescriptions

---

## Source Model

Park, J. et al. (2019)

**Inferring the astrophysics of reionization and cosmic dawn from galaxy luminosity functions and the 21-cm signal**

Primary reference for:

- Stellar mass prescriptions
- Star formation rates
- UV luminosity modelling
- Escape fractions

---

## Original 21cmFAST Framework

Mesinger, A., Furlanetto, S., & Cen, R. (2011)

**21cmFAST: A Fast, Semi-Numerical Simulation of the High-Redshift 21-cm Signal**

Original description of:

- Halo sampling
- Density evolution
- Reionisation framework

---

# Verification Checklist

Before using a new halo catalogue:

```python
print(catalog.n_halos)

print(catalog.halo_masses.min(),
      catalog.halo_masses.max())

print(catalog.sfr.min(),
      catalog.sfr.max())

print(catalog.stellar_masses.min(),
      catalog.stellar_masses.max())
```

Verify that:

- Halo masses are in a physically plausible range.
- SFR values are positive.
- Stellar masses are smaller than halo masses.
- Coordinates lie within the simulation box.

Only then proceed to galaxy selection and bias calculations.