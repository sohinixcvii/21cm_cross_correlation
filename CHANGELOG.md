# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

---

## [Unreleased]

### Added
- **HPC lightcone pipeline** — the monolithic `21cmfast_HERAxEuclid_lightcone.ipynb`
  has been refactored into three self-contained parts for efficient cluster use:

  - **`run_simulation.py`** (Part 1 — batch script): runs the 21cmFASTv4 lightcone
    simulation, constructs the galaxy density field from `halo_sfr`, estimates the
    galaxy bias via HMF integration over the Euclid UV magnitude range, applies
    Kaiser redshift-space distortions in Fourier space, and saves all outputs to
    `outputs/lightcone_data.h5` with gzip compression. Uses `matplotlib.use("Agg")`
    for headless HPC execution. No logic changes from the notebook — same algorithms,
    same parameters, same comments.

  - **`notebooks/plot_fields.ipynb`** (Part 2 — visualisation notebook): loads
    `outputs/lightcone_data.h5` and reproduces all field plots from the original
    notebook: halo catalogue scatter plots and SFR distributions (Cells 10–11),
    the three-panel lightcone slice (Cell 18), and the wide-format EoR brightness
    temperature plot (Cell 20). Gracefully skips halo catalogue cells when no
    catalogue is available (synthetic fallback).

  - **`notebooks/analysis.ipynb`** (Part 3 — calculation notebook): loads
    `outputs/lightcone_data.h5` and performs all post-simulation calculations:
    `compute_cylindrical_cross_power` (Cell 22), foreground wedge geometry and
    power spectrum plots (Cells 23–24), photo-$z$ damping and wedge excision
    (Cell 26), per-mode SNR map and total detection significance (Cell 28), and
    the summary table (Cell 29).

  - **`submit_job.sh`**: SLURM batch submission script that activates the
    `21cmfast` conda environment and runs `run_simulation.py`. Configurable
    wall-time, memory, and partition via `#SBATCH` directives.

  All scalar parameters and metadata are stored as HDF5 attributes in
  `outputs/lightcone_data.h5` so that Parts 2 and 3 require no configuration
  beyond pointing at the output file. The kernelspec for both notebooks is
  set to `21cmfast`.



### Added
- **`src/conversions.py`** — new module of cosmological conversion utilities for
  high-redshift galaxy surveys. Functions:
  - `Muv_to_Luv(Muv)` — converts absolute UV AB magnitude to monochromatic UV
    luminosity in erg s⁻¹ Hz⁻¹ (Oke & Gunn 1983; Madau & Dickinson 2014).
  - `Luv_to_Muv(Luv)` — inverse conversion from luminosity to AB magnitude.
  - `survey_area_from_volume(volume_mpc3, z_min, z_max, cosmo=None)` — infers
    the sky area in deg² that corresponds to a given comoving survey volume
    (Mpc³) over a redshift interval, using Simpson integration of the
    differential comoving volume $\mathrm{d}V/\mathrm{d}z\,\mathrm{d}\Omega$
    (Hogg 1999; Astropy).
  - `area_deg2_to_steradians(area_deg2)` — unit conversion from deg² to sr.
  - `volume_from_area(area_deg2, z_min, z_max, cosmo=None, n_z=1000)` —
    computes the comoving survey volume (Mpc³) for a given sky area and
    redshift range; inverse of `survey_area_from_volume`.

  All functions accept scalar or array inputs and default to the Planck18
  cosmology from astropy; a custom `astropy.cosmology` object may be passed
  via the `cosmo` argument.

### Fixed
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — lightcone slice orientation:**
  Panels 2 and 3 of Section 5 previously had the LOS and transverse axes
  swapped. The `brightness_temp_field[:, mid_y, :]` and
  `neutral_fraction[:, mid_y, :]` arrays (shape `(HII_DIM, N_z)`) were
  incorrectly transposed with `.T`, placing the LOS on the y-axis and the
  transverse direction on the x-axis. Removed the transpose so that
  `imshow(origin="lower")` correctly maps rows→transverse (y) and cols→LOS (x),
  consistent with every published 21cm lightcone figure. The `extent` parameter
  was updated accordingly to `[lc_dist_Mpc[0], lc_dist_Mpc[-1], 0, BOX_LEN]`.

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — secondary redshift axis:**
  The redshift annotation used `twinx()`, adding a second y-axis on the right.
  After the orientation fix the LOS is on the x-axis, so the annotation must
  use `twiny()` with tick positions computed via
  `np.interp(z_ticks, lc_redshifts, lc_dist_Mpc)`.

### Added
- **`21cmfast_HERAxEuclid_lightcone.ipynb` — brightness temperature evolution
  plot (Section 5b):** New wide-format (16×3.5") lightcone slice cell styled
  after the canonical Mesinger & Furlanetto (2007) figure. Features:
  - Custom EoR colourmap (`EoR21`): dark blue-black (ionised, δT_b ≈ 0) →
    blue → cyan → yellow → orange → near-white (neutral, high δT_b)
  - LOS (comoving distance) on the x-axis; transverse distance on the y-axis
  - Dual x-axes: comoving distance [Mpc] on the bottom, redshift z on the top
    (using `twiny()` with interpolated tick positions)
  - Title reports the observed frequency range
    ($f_\mathrm{obs} = F_{21}/(1+z)$) alongside z and ⟨x_HI⟩

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — configurable minimum LOS slices:**
  `minimum_los_slices` (default 100) replaces the previous hardcoded floor of
  10. For narrow redshift ranges where the natural slice count (L_LOS / cell_size)
  is small, this ensures the k_∥ grid is adequately sampled.
  Validated with z_min=6.5, z_max=6.505, 100 slices: notebook executes without
  errors, 21cmFAST returns a (128, 128, 100) box (L_LOS = 200 Mpc after rounding
  to cell boundaries), empty bins fall from 315/400 to 145/400, and the
  large-scale cross-spectrum sign is correctly negative (anti-correlated).

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — Section 3b synthetic galaxy field
  fallback:** When the lightcone `halo_sfr` field is all-zero (e.g. 21cmFAST
  not installed or CHMF-SAMPLER placed no halos), Section 3b automatically
  generates a synthetic galaxy overdensity by Poisson-sampling the matter
  density field. This allows the power-spectrum and SNR cells to run
  end-to-end without a 21cmFAST installation.

- **`21cmfast_HERAxEuclid_lightcone.ipynb` — galaxy bias estimation (Cell 13):**
  New cell estimates the luminosity-weighted linear bias by integrating a
  Schechter UV luminosity function over the Euclid magnitude range
  $M_\mathrm{UV} \in [-24, -18]$. The result informs the default
  `galaxy_bias` value used in the Kaiser RSD correction.

### Added
- **`21cmfast_HERAxEuclid_lightcone.ipynb`** — lightcone counterpart to the
  coeval simulation notebook. Uses `RectilinearLightconer` + `run_lightcone`
  (21cmFASTv4) to produce a self-consistent lightcone over a configurable
  redshift range ($z_\mathrm{min}$–$z_\mathrm{max}$, default 6.5–7.5).
  Key differences from the coeval version:
  - Non-cubic $(N_\perp \times N_\perp \times N_z)$ box with separate transverse
    and LOS cell sizes handled throughout
  - Fields accessed via `lightcone.lightcones['field_name']`; neutral fraction
    stored as `'neutral_fraction'` (not `'xH_box'` as in coeval)
  - Visualisation includes LOS slice panels showing redshift evolution of
    $\delta T_b$ and $x_\mathrm{HI}$ alongside a transverse (x–y) slice
  - Updated `compute_cylindrical_cross_power` accepts separate `box_len_perp`
    and `box_len_los` arguments for non-cubic boxes
  - All user-adjustable parameters consolidated into a single clearly marked
    configuration cell

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
