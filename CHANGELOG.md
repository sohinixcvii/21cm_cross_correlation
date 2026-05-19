# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased]

### Fixed
- `21cmfast_HERAxEuclid.ipynb`
  - Updated foreground and horizon wedge prescription
  - Fixed Nan appearance due to log binning in 2d power spectra plots 

### Added
- `21cmfast_HERAxEuclid.ipynb` — new notebook demonstrating an end-to-end
  HERA × Euclid 21 cm–galaxy cross-correlation workflow using 21cmFASTv4
  (Davies et al. 2025, arXiv:2504.17254). Covers simulation, galaxy field
  construction, 2D cylindrical power spectra, photo-z damping, foreground
  wedge excision, and per-mode SNR estimation.

### Fixed
- **`21cmfast_HERAxEuclid.ipynb` — Cell 4 (galaxy field construction):**
  The cell previously raised an unconditional `RuntimeError` when 21cmFAST
  was installed, blocking all subsequent cells. Fixed by implementing the
  galaxy density field using `coeval.halobox.get('halo_sfr')`, the per-cell
  SFR density field provided by the 21cmFASTv4 `HaloBox` API
  (`py21cmfast >= 4.1.1`). The SFR density field is converted to an
  overdensity $\delta_{\rm gal} = \rm SFR / \langle SFR \rangle - 1$, which
  correctly traces the galaxy distribution and produces the expected
  large-scale anti-correlation with the 21 cm brightness temperature field
  (cross-spectrum sign: **NEGATIVE** on large scales ✓).

### Changed
- **`21cmfast_HERAxEuclid.ipynb` — Section 3 markdown:** Updated the list of
  coeval output fields to reflect the actual 21cmFASTv4 API: `halo_field`
  (incorrect) → `halobox` (the `HaloBox` object containing per-cell gridded
  quantities including `halo_sfr`, `n_ion`, etc.).
- **`21cmfast_HERAxEuclid.ipynb` — Section 4 markdown:** Clarified that the
  `HaloBox` API exposes cell-averaged quantities rather than individual halo
  catalogues, so a strict per-halo $M_{\rm UV}$ cut requires a lightcone
  post-processing step. Added a note explaining why the SFR-density proxy is
  a valid and physically motivated tracer of the Euclid-observable galaxy
  population.

---

## Notes on 21cmFASTv4 `HaloBox` API

In 21cmFAST v4.1+, `coeval.halobox` is a `HaloBox` object whose arrays are
accessed via `.get('<field_name>')`. Available fields include:

| Field | Description |
|-------|-------------|
| `halo_sfr` | Total SFR per cell, summed over all halos [internal units] |
| `n_ion`    | Number of ionizing photons per cell |

Individual halo positions and UV magnitudes are not exposed by this API.
For per-halo catalogues, use the raw halo output from a lightcone run.
