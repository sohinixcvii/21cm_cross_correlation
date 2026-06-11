# Galaxy Bias Calculation from 21cmFAST Halo Catalogues

## Overview

This document describes the methodology used to estimate the effective linear galaxy bias for a Euclid-like galaxy sample generated from a 21cmFASTv4 simulation.

The calculation begins with the halo catalogue produced by 21cmFAST, applies an observational selection based on UV magnitude limits, and computes the mean bias of the selected galaxy population.

---

# Scientific Motivation

The galaxy bias parameter describes how strongly galaxies trace the underlying matter density field.

In linear theory,

b_g = δ_g / δ_m

where:

* δ_g is the galaxy overdensity field
* δ_m is the matter overdensity field
* b_g is the linear galaxy bias

Massive galaxies tend to form in massive dark matter halos, which are themselves more strongly clustered than the matter field. Consequently, galaxy samples are generally biased tracers of the underlying matter distribution.

For a flux-limited survey such as Euclid, the observed galaxy sample contains only galaxies brighter than a limiting magnitude. The resulting galaxy bias therefore depends on the luminosity selection function.

---

# Input Data

The calculation uses the following quantities from the 21cmFAST perturbed halo catalogue:

| Quantity            | Attribute        | Units          |
| ------------------- | ---------------- | -------------- |
| Halo mass           | `halo_masses`    | Msun           |
| Star formation rate | `sfr`            | Msun yr^-1     |
| Halo positions      | `halo_coords`    | Mpc (comoving) |
| Stellar mass        | `stellar_masses` | Msun           |

Only `halo_masses` and `sfr` are required for the bias calculation.

---

# Step 1: Convert SFR to UV Luminosity

The UV luminosity is estimated from the star formation rate using the standard UV-SFR calibration.

Following Kennicutt (1998) and Madau & Dickinson (2014),

L_UV = SFR / K_UV

where

K_UV = 1.15 × 10^-28

with units

K_UV = Msun yr^-1 / (erg s^-1 Hz^-1)

and therefore

L_UV has units erg s^-1 Hz^-1.

This relation assumes:

* continuous star formation
* a standard stellar IMF
* negligible dust attenuation

at rest-frame UV wavelengths around 1500 Å.

---

# Step 2: Convert UV Luminosity to Absolute Magnitude

UV luminosity is converted to an AB absolute magnitude using

M_UV = -2.5 log10(L_UV) + 51.60

where

* L_UV is in erg s^-1 Hz^-1
* M_UV is the absolute AB magnitude

This is simply the standard AB magnitude definition evaluated at a distance of 10 pc.

---

# Step 3: Apply the Euclid Selection Function

The Euclid survey is approximated as a simple magnitude-limited sample.

Only galaxies satisfying

-22.66 ≤ M_UV ≤ -18

are retained.

The lower limit corresponds to the brightest galaxies in the survey volume, while the upper limit represents the faint-end detection threshold.

This step defines the galaxy selection function

φ(M_h, z)

which specifies the probability that a halo of mass M_h hosts an observable galaxy.

In the current implementation,

φ(M_h, z) = 1

for galaxies passing the magnitude cut, and

φ(M_h, z) = 0

otherwise.

---

# Step 4: Compute Halo Bias

For each halo mass, a linear halo bias is assigned using the Sheth-Tormen formalism.

The halo bias is

b_h(M,z)

and is calculated from the peak-height parameter

ν = δ_c / σ(M,z)

where

* δ_c = 1.686
* σ(M,z) is the linear density variance

The Sheth-Tormen bias relation is

b_h(ν) =
1 +
(aν² − 1)/δ_c +
2p/[δ_c(1 + (aν²)^p)]

with

* a = 0.707
* p = 0.3

This bias model is evaluated using the same cosmology and redshift as the simulation.

The implementation uses the `hmf` package to compute ν(M,z) and interpolate the corresponding bias values.

---

# Step 5: Compute the Effective Galaxy Bias

The effective galaxy bias is formally

b_g =
∫ (dn/dM) b_h(M,z) φ(M,z) dM
/
∫ (dn/dM) φ(M,z) dM

where

* dn/dM is the halo mass function
* b_h(M,z) is the halo bias
* φ(M,z) is the galaxy selection function

This expression represents a number-weighted average of halo bias over all observable galaxies.

---

# Catalogue-Based Estimator

Since the full halo catalogue is available, the integral can be evaluated directly using the selected halos.

For N selected galaxies,

b_g =
(1/N) Σ b_h(M_i,z)

where

M_i

is the host halo mass of galaxy i.

This approach avoids any assumptions about the analytic form of the selection function and uses the actual simulated galaxy population.

---

# Output

The final quantity reported is

<b_g>

the mean effective linear galaxy bias of the Euclid-selected galaxy sample at the simulation redshift.

Typical values are expected to lie between

b_g ≈ 3–10

for bright galaxies at

z ≈ 6–8

depending on the luminosity threshold.

---

# Assumptions and Limitations

The current implementation assumes:

1. UV luminosity depends only on instantaneous SFR.

2. Dust attenuation is neglected.

3. The Kennicutt UV-SFR calibration remains valid at high redshift.

4. Halo occupation is deterministic once the magnitude cut is applied.

5. Galaxy bias is entirely determined by halo mass.

Future improvements could include:

* dust corrections
* scatter in the SFR–UV relation
* abundance matching
* stellar-mass–halo-mass relations
* halo occupation distributions (HODs)
* direct measurement of bias from power spectra

---

# References

## Galaxy Bias

Mo, H. J. & White, S. D. M. (1996)
"A high-redshift bias model"
MNRAS, 282, 347

Sheth, R. K., Mo, H. J., & Tormen, G. (2001)
"Ellipsoidal collapse and an improved model for the number and spatial distribution of dark matter haloes"
MNRAS, 323, 1

---

## Halo Mass Functions

Murray, S. G., Power, C., & Robotham, A. S. G. (2013)
"HMF: Halo Mass Function calculations"
Astronomy and Computing, 3–4, 23

---

## UV Luminosity Calibration

Kennicutt, R. C. (1998)
"Star Formation in Galaxies Along the Hubble Sequence"
ARA&A, 36, 189

Madau, P. & Dickinson, M. (2014)
"Cosmic Star Formation History"
ARA&A, 52, 415

---

## 21cmFAST

Park, J. et al. (2019)
"21cmFAST: A Fast, Semi-Numerical Simulation of the High-Redshift 21-cm Signal"
MNRAS, 484, 933

Davies, J. E. et al. (2025)
21cmFAST v4 documentation and source model implementation.
