# TODO

Outstanding work for the 21 cm × Euclid galaxy cross-correlation pipeline,
priority-ordered. Each item states the problem, where it lives, and how to
verify the fix.

**Last updated:** 2026-08-25

> **Context:** the pipeline plumbing (`run_pipeline.py`, `src/analysis.py`,
> `src/figures.py`, `src/dataio.py`, `tests/`) is complete and tested. What
> remains is science. See `CHANGELOG.md` for what has already been fixed and
> `docs/project_update.md` for the current numbers.

---

## P0 — Make the estimator lightcone-ready, then widen Δz

> **P0.1–P0.4 are implemented (2026-08-25), behind a toggle.**
> `ESTIMATOR = "coeval"` (the default) is the historical formalism and
> reproduces every earlier number exactly; `ESTIMATOR = "lightcone"` applies
> all four fixes. Set it in `run_simulation.py`'s configuration block, or per
> analysis run with `run_pipeline.py --estimator lightcone`. **P0.5 — actually
> widening Δz — is still open**, and is now the only thing standing between
> this pipeline and a production lightcone.

**Current state: `run_simulation.py` runs a deliberately narrow smoke-test
slab, z = 6.995 → 7.005 (Δz = 0.01).** A production range was briefly set to
z = 6.5 → 7.5 and then **reverted on purpose**, because widening Δz without
the work below produces power spectra that cannot be trusted.

The power-spectrum estimator in `src/analysis.py` was inherited from the
**coeval** notebook and assumes the box is statistically homogeneous along the
line of sight. At Δz = 0.01 the box is quasi-coeval, so that assumption holds
and the estimator is appropriate — configuration and formalism currently
match. At Δz = 1.0 it breaks in the four ways below, each measured against
Planck18 for z = 6.5 → 7.5 with N_z = 175.

**P0.1–P0.4 are done** (2026-08-25); each subsection below records what was
implemented and how to verify it. They were the prerequisites for P0.5, which
remains.

> **What follows.** Each P0 subsection keeps its original problem statement —
> that is the record of *why* the change was made — followed by an
> "Implemented" note at the top saying what was actually done. The remaining
> work is P0.5: flip `ESTIMATOR` and the redshift range, and re-run with
> `--sim force` (P0.1 changes what is written to disk).

### These issues are already live in notebook 3

`21cmfast_HERAxEuclid_lightcone.ipynb` (the monolithic lightcone notebook)
runs **z = 6.5 → 7.5 today**, builds `lc_redshifts = np.linspace(z_min, z_max,
N_z)`, and carries its own copy of `compute_cylindrical_cross_power`. So it
has all four problems below, right now — unlike the pipeline, which is
quasi-coeval and therefore self-consistent.

Its power spectra should still be treated with that caution: the fixes landed
in `src/analysis.py`, and the notebook carries its own copy. It should be
switched to the shared implementation rather than fixed in parallel — the
shared one now takes `mean_subtraction`, `taper` and the sub-band split as
arguments, so the notebook only needs to call it.

### Caveat on the current slab

Even at Δz = 0.01 the smoke-test geometry limits what the spectra mean: N_z =
100 slices are interpolated from only 5 node redshifts across 3.5 Mpc, so LOS
structure below the 2 Mpc cell size is interpolation, not signal. The k_∥ axis
therefore extends well past the scale where independent information exists.
This is why the total SNR of 0.1 σ is not a physical forecast. It is a
limitation of the configuration, not a bug, and P0.5 is what fixes it.

### P0.1 Non-uniform LOS sampling breaks the FFT — ✅ done

**Implemented:** `LIGHTCONE_SAMPLING = "comoving"` in `run_simulation.py`
builds the lightconer with `RectilinearLightconer.between_redshifts(...)` at
`resolution = cell_size`, so slices are evenly spaced in comoving distance.
`"redshift"` keeps the historical `linspace(z_min, z_max, N_z)`. Both are
recorded in the HDF5 as `lightcone_sampling`. Note that `minimum_los_slices`
does not apply on the comoving path — the slice count follows from the
resolution.

`run_simulation.py` builds `lc_redshifts = np.linspace(z_min, z_max, N_z)` —
uniform in *redshift*, therefore **not** uniform in comoving distance, because
`dD/dz = c/H(z)` varies across the box.

Measured over z = 6.5 → 7.5 with N_z = 175 (Planck18) — i.e. what would happen
if Δz were widened today:

| Quantity | Value |
|---|---|
| First cell (low z) | 2.2137 Mpc |
| Last cell (high z) | 1.8380 Mpc |
| **Spread** | **20.4 %** |
| What `compute_cylindrical_cross_power` assumes | a single uniform 2.0047 Mpc |

For contrast, the old Δz = 0.01 slab had a 0.19 % spread — invisible.

An FFT on unevenly sampled data mis-assigns k_∥ and leaks power between modes.

**Fix:** build the lightconer with equally spaced comoving slices. The v4 API
has a constructor for exactly this:

```python
# instead of RectilinearLightconer(lc_redshifts=..., ...)
lightconer = p21c.RectilinearLightconer.between_redshifts(
    min_redshift=z_min,
    max_redshift=z_max,
    resolution=cell_size * units.Mpc,   # matches the transverse cell
    quantities=("brightness_temp", "density", "neutral_fraction", "halo_sfr"),
)
```

It computes `lc_distances = np.arange(d_min, d_max + res, res)` internally.

- **Where:** `run_simulation.py` §2 (lightconer construction)
- **Verify:** `np.diff(lc_dist_Mpc)` should be constant to machine precision.
  Add an assertion or a test on the written HDF5.

### P0.2 Global mean subtraction leaves a spurious LOS gradient — ✅ done

**Implemented:** `analysis.subtract_field_mean(field, mode)` with
`mode="per_slice"`, applied to *both* fields by
`compute_all_power_spectra(..., mean_subtraction="per_slice")`, and to
δ_gal at construction time by `GALAXY_MEAN_SUBTRACTION = "per_slice"` in
`run_simulation.py`.

**Measured, and smaller than this section claimed.** A line-of-sight ramp that
is uniform across the sky puts ~99 % of its power at k_⊥ = 0, and the
log-spaced binning starts at 0.5 Δk_⊥ — so that column never enters a bin, and
the two modes agree to floating-point precision on such a field
(`tests/test_analysis.py::test_uniform_ramp_lives_at_k_perp_zero_and_is_already_binned_out`).
The operation is still correct and matters for anything that uses k_⊥ = 0, but
it is **not** by itself the fix for low-k_∥ contamination. What genuinely
couples redshift evolution into non-zero k_⊥ is that δT_b = T_0(z)[1 + δ]: an
evolving mean *modulates* the fluctuation amplitude. Removing a per-slice mean
does not undo a per-slice gain — per-slice **normalisation** would. Left as a
follow-up rather than folded in silently; see P1.4.

`src/analysis.compute_all_power_spectra` does:

```python
t21_fluctuations = brightness_temp_field - brightness_temp_field.mean()
```

A single global scalar. But ⟨T_b⟩ evolves strongly with redshift across the
lightcone, so the residual carries a large monotonic LOS ramp. That ramp is
not signal — it aliases directly into low-k_∥ power, which is precisely where
the foreground-wedge analysis looks.

**Fix:** subtract the per-slice mean,
`T_b(x, y, z) − ⟨T_b⟩(z_slice)`, i.e. subtract along axes (0, 1) keeping the
LOS axis. Decide explicitly whether the galaxy field needs the same treatment
(it is already an overdensity, but `δ_gal` is built from a global
`sfr_field.mean()` in `run_simulation.py` §3b and has the same issue).

- **Where:** `src/analysis.py:compute_all_power_spectra`;
  `run_simulation.py` §3b for the galaxy field
- **Verify:** the k_∥ → 0 power should drop substantially; compare the
  P_21(k_∥) profile before and after.

### P0.3 The signal evolves across the box — a single FFT mixes epochs — ✅ done

**Implemented:** `analysis.compute_subband_power_spectra(...)` splits the
line of sight into contiguous sub-bands, applies a four-term Blackman-Harris
taper within each (`analysis.blackman_harris_taper`, amplitude restored by
dividing the power by ⟨w²⟩), and returns one `PowerSpectra` per band together
with a `SubbandGeometry` carrying each band's effective redshift, comoving
extent and bandwidth. `run_pipeline.py` then computes one uncertainty budget
per band at its own z_eff and combines the totals with
`analysis.combine_band_snr` — a quadrature sum, since the bands sample
disjoint volumes. The effective redshift comes from the band's mean *observed
frequency*, not its mean redshift.

Even with P0.1 and P0.2 fixed, one FFT over Δz = 1.0 returns a
redshift-*averaged* power spectrum with an ill-defined effective redshift. In
the EoR, x_HI and hence P_21 change by a large factor over Δz = 1, so the
result is not the power spectrum "at z = 7".

**Fix (standard approach):** subdivide the lightcone into redshift chunks and
compute the power spectrum per chunk, then report each at its own effective
redshift. Apply a taper (e.g. Blackman-Harris) along the LOS within each chunk
to control spectral leakage from the finite band.

The natural chunk width is already in the config — see P0.4.

- **Where:** `src/analysis.py`, new function alongside
  `compute_all_power_spectra`; `run_pipeline.py` would loop over chunks
- **References:** Datta et al. (2012), MNRAS 424, 1877; Mondal et al. (2018),
  MNRAS 474, 1390; La Plante et al. (2023), arXiv:2205.09770

### P0.4 The PS bandwidth and the noise bandwidth disagree — ✅ done

**Implemented:** each sub-band's own frequency span is passed to
`compute_uncertainty_budget(bandwidth=...)`, so signal and noise are computed
over the same volume by construction. The band count is
`ceil(span / bandwidth)`, defaulting to the noise model's own bandwidth;
`--subband-bandwidth` overrides it. A lightcone narrower than one bandwidth
returns a single band, so the smoke-test slab still works.

The noise model uses `bandwidth = 8e6` Hz (8 MHz). The power spectrum is now
measured over the whole lightcone, which spans:

| Quantity | Value |
|---|---|
| Frequency range (z = 7.5 → 6.5) | 167.11 – 189.39 MHz |
| **Total bandwidth** | **22.28 MHz** |
| Bandwidth assumed by the noise | 8.00 MHz |
| **Mismatch** | **2.8×** |

So the signal and the noise are currently computed over different volumes,
which makes the SNR internally inconsistent.

**Fix:** compute the power spectrum in ~8 MHz sub-bands (≈ 2.8 chunks across
the box), matching the noise assumption. This resolves P0.3 and P0.4 together
and mirrors what real HERA analyses do.

- **Where:** `src/analysis.py`, `run_pipeline.py`
- **Verify:** total SNR should become a quadrature sum over bands, each with a
  well-defined effective redshift.

### P0.5 Widen Δz to a true lightcone — **now unblocked**

P0.1–P0.4 are in place, so this is the remaining step. Switch
`run_simulation.py` to the production range **and** to the lightcone
estimator:

```python
ESTIMATOR = "lightcone"     # currently "coeval"
```

```python
z_min = 6.5            # currently 6.995
z_max = 7.5            # currently 7.005
```

Verified geometry for that range (computed against Planck18, no 21cmFAST run
required — see the derivation in `CHANGELOG.md`):

| Quantity | Smoke test (current) | Production (target) |
|---|---|---|
| Δz | 0.01 | 1.0 |
| L_LOS | 3.50 Mpc | **350.83 Mpc** |
| N_z | 2 natural → 100 (floor binds) | **175 natural** (floor idle) |
| LOS cell | 0.035 Mpc (57× oversampled) | 2.005 Mpc (matches transverse) |
| Node redshifts | 5 | 10 |
| Frequency span | ~0.9 MHz | 22.28 MHz |

Cost: ten full coeval evaluations instead of five, and 175 LOS slices instead
of 100. Run it as `bash submit_job.sh --sim force` on the cluster, not on a
workstation.

- P0.1 changes what gets written to the HDF5, so this needs `--sim force`,
  not the default `--sim auto`.
- The footprint-derived box now implies z = 6.55–7.45 rather than 6.5–7.5
  (Δz = 0.90, L_LOS = 315.6 Mpc, N_z = 166, 9 nodes, 20.04 MHz — see
  `docs/simulation_spec.md`), so the table above is superseded by that
  document's §2.2.
- **Blocked on the `INT_MAX` halo-catalogue overflow** at the current box
  size; `docs/simulation_spec.md` §7.1 has the three ways out and recommends
  a 350 Mpc first run.

---

## P1 — Known analysis gaps

### P1.1 Mode counts are computed but never used in the SNR

`compute_cylindrical_cross_power` bins `mode_counts`, `dataio` caches it in
`analysis_products.h5` — and it is only ever consumed as an "empty bins"
diagnostic (`run_pipeline.py:592`). La Plante Eqs. 15–17 divide the per-bin
variance by the number of independent modes. Without that factor σ is
overestimated and the total SNR is biased **low**.

- **Where:** `src/analysis.py:cross_power_snr`
- **Action:** check the exact form against La Plante et al. (2023) Eqs. 15–17,
  then apply the 1/N_modes weighting.
- **Note:** inherited faithfully from `notebooks/analysis.ipynb` — not a
  regression, but now fixable in one place.

### P1.2 Thermal noise is a scaling estimate

`src.analysis.hera_thermal_noise_power` uses
`T_sys² / (t_int Δν)` with `T_sys = 100 K + 60 K (300 MHz/ν)^2.55`.

The proper La Plante Eq. 11 form — with the baseline density `n(k_⊥)` and the
`X² Y Ω′` cosmological factors — **is already implemented** in
`21cm_galaxy_cross_uncertainty.ipynb` (see its `n_baselines()` function). It
needs moving into `src/` rather than writing from scratch.

- **Where:** `src/analysis.py`; source material in notebook 1
- **Publication-grade alternative:**
  [21cmSense](https://github.com/rasg-affiliates/21cmSense)

### P1.4 Per-slice normalisation, not just per-slice mean subtraction

P0.2 removes the evolving *mean*; it does not remove the evolving *gain*.
Since δT_b = T_0(z)[1 + δ(x)], the redshift evolution multiplies the
fluctuations, and that modulation does reach non-zero k_⊥ where per-slice mean
subtraction cannot touch it. The fix is to normalise each slice by its own
mean rather than subtract it, which changes what the amplitude of P_21 means
and so is a deliberate, separate decision.

- **Where:** `src/analysis.py:subtract_field_mean`, plus a matching option in
  `compute_all_power_spectra`
- **Verify:** on a synthetic box with an imposed T_0(z) ramp, the recovered
  P_21 should match the un-ramped one to within sample variance

### P1.3 `apply_rsds=False`

`run_simulation.py:230` disables 21cmFAST's own redshift-space distortions;
Kaiser RSD is applied analytically afterwards in §5 using β = f/b_g. Enabling
`apply_rsds=True` would make the RSDs self-consistent and remove the need for
the β approximation entirely.

- **Caveat:** the analytic Kaiser step would then have to be removed, not
  merely skipped, or RSDs would be applied twice.

---

## P2 — Coverage and infrastructure

### P2.1 `src/conversions.py` and `src/FOV_to_cMpc.py` are untested

The project convention (`CLAUDE.md`) is a test per function. These two modules
predate `tests/` and are only exercised indirectly. The round-trip identities
are the obvious start:

- `Muv_to_Luv` ↔ `Luv_to_Muv`
- `sfr_to_Luv` ↔ `Luv_to_sfr`
- `survey_area_from_volume` ↔ `volume_from_area`

`sheth_tormen_bias` already has a regression guard in `tests/test_analysis.py`
(the ν-convention test).

### P2.2 Notebook 1 lives outside the pipeline

`21cm_galaxy_cross_uncertainty.ipynb` implements the full La Plante variance
framework (Eqs. 15–17), the BBKS transfer function, and the Lidz+09 bias
model. None of it is importable or tested, and it now duplicates `T0`, `H_z`,
and `W_photoz` against `src/analysis.py`.

Extracting it into `src/` would fix P1.2 at the same time.

### P2.3 No parameter-sweep driver

`notebooks/analysis.ipynb` calls for sweeps over σ_z and the magnitude limits
to optimise survey design. The hooks exist (`--m-uv-bright`, and cached
spectra make re-analysis ~1.6 s) but nothing loops over them.

- **Suggested:** a `--sweep` mode writing one summary JSON row per parameter
  combination.

---

## P3 — Deferred (compute-bound, needs an HPC allocation)

### P3.1 1 Gpc box

Davies et al. (2025) use a 1 Gpc box at the same 2 Mpc cell. The current
256 Mpc box is why the UVLF bright end is noisy and why only ~49 k halos
survive the Euclid cut.

**This is a storage and compute problem, not a code change** — `BOX_LEN` and
`HII_DIM` are already configuration variables, and `--max-halos` exists for
catalogues of this size.

| Quantity | Current (256 Mpc) | Target (1 Gpc) | Factor |
|---|---|---|---|
| `HII_DIM` / `DIM` | 128 / 384 | 500 / 1500 | — |
| Comoving volume | 1.7 × 10⁷ Mpc³ | 1.0 × 10⁹ Mpc³ | **59.6×** |
| Halos (volume-scaled) | 1.14 × 10⁸ | ~6.8 × 10⁹ | 59.6× |
| Halo catalogue on disk | 2.74 GB | **~163 GB** | 59.6× |
| One high-res IC array | 0.23 GB | **13.5 GB** | 59.6× |

Storage binds first. Before attempting:

- Confirm ~200 GB scratch quota per run (`docs/INSTALL_21cmFASTv4.md`).
- Consider storing only Euclid-selected halos (~0.04 % → ~65 MB, not 163 GB).
- Try a 512 Mpc intermediate first (`HII_DIM = 256`, 8× volume, ~22 GB
  catalogue) to measure the scaling.

Full costing in `docs/project_update.md` §12.

---

## P4 — Modelling extensions

From `docs/Galaxy_bias_formalism.md` ("Assumptions and Limitations"). The
current bias model assumes UV luminosity depends only on instantaneous SFR,
neglects dust, treats halo occupation as deterministic, and derives bias
purely from halo mass.

Possible extensions, roughly in order of impact:

- Dust attenuation (matters most at the bright end, where Euclid selects)
- Scatter in the SFR–UV relation
- Measure bias directly from `P_gal / P_matter` rather than from the HMF —
  the pipeline already computes both auto-spectra
- Abundance matching / stellar-mass–halo-mass relations
- Halo occupation distributions (HODs)

---

## Completed

Kept for context; see `CHANGELOG.md` for detail.

- ✅ End-to-end driver (`run_pipeline.py`) with `auto`/`force`/`skip` staging
- ✅ Notebook science extracted into tested modules (`src/`, 69 tests)
- ✅ Sheth-Tormen ν convention fixed — `b_g` 33.39 → 4.23 analytic
- ✅ SFR timescale 100 Myr → `t_STAR × t_H` = 570 Myr
- ✅ Bias measured from the halo catalogue (adopted `b_g` = 4.744);
  β_rsd 0.0299 → 0.2103
- ✅ κ_UV and AB zero point unified via `src/conversions.py`
- ✅ González+10 M★–M_UV relation re-enabled
- ✅ `docs/project_update.md` rewritten with corrected diagnosis
