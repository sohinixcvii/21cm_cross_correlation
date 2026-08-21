# The Uncertainty Budget — Complete Specification

**What this document is.** The parameter-level reference for the 21 cm ×
galaxy cross-power uncertainty budget as the HPC pipeline computes it: every
formula with its provenance and its evaluated number at $z_\mathrm{obs} = 7$,
the term-by-term audit against the notebook the calculation came from, the
three discrepancies that audit found, what is written to disk, and what the
calculation still does not do.

**Verified:** 2026-08-17, branch `test/methods`. Test suite: **99 passed**
(`conda run -n 21cmfast pytest tests/ -q`), of which 23 cover the budget — 16
in `tests/test_uncertainty_budget.py`, six CLI tests in `tests/test_pipeline.py`,
and one figure test.

**Source of the calculation.** `21cmfast_HERAxEuclid_lightcone.ipynb` —
its *"Photo-z Radial Smearing & Foreground Wedge Mask"* and
*"Cross-correlation SNR map"* sections. The formalism is La Plante et al.
(2023), Eqs. 15–17.

**Companion documents.** [`HPC.md`](HPC.md) is the parameter-level ground
truth for the whole run (§5.2–5.5 cover this budget in the context of the
stored simulation); [`../PIPELINE.md`](../PIPELINE.md) is the short version;
[`../TODO.md`](../TODO.md) is the outstanding work.

---

## 0. TL;DR — what the budget is

Four transformations applied in order to the three power spectra, ending in a
single detection significance:

```
P_21, P_gal, P_×  ──①─►  photo-z damping   P_×W ,  P_gal W²
                  ──②─►  wedge excision    97 of 400 bins survive
                  ──③─►  noise             P_N,21 = 3.75 ,  P_N,gal = 333.3
                  ──④─►  variance + SNR    σ_×  ,  total SNR
```

One function does all four: **`src.analysis.compute_uncertainty_budget`**.
It returns an `UncertaintyBudget` holding every intermediate term, so nothing
in the chain is recomputed or reinterpreted downstream.

At the stored run's configuration ($z_\mathrm{obs} = 7$, $\sigma_z = 0.45$,
buffer $= 0.0677$ Mpc⁻¹) the answer is **$1.06 \times 10^{-111}\sigma$** — not
a detection, and the reason is structural, not statistical: see §7.1.

---

## 1. Where it lives

| Concern | Location |
|---|---|
| Individual formulas | `src/analysis.py` §§3–5 — `radial_smearing_length`, `photoz_damping_kernel`, `horizon_wedge_slope`, `fov_wedge_slope`, `foreground_wedge_mask`, `system_temperature`, `hera_thermal_noise_power`, `cross_power_snr`, `total_snr` |
| The assembled chain | `src/analysis.py` §5b — `compute_uncertainty_budget`, `UncertaintyBudget` |
| Driver stage | `run_pipeline.py` — `observational_stage` |
| Persistence | `src/dataio.py` — `save_uncertainty_budget`, `load_uncertainty_budget` |
| Figure | `src/figures.py` — `plot_uncertainty_budget` |
| Interactive use | `21cmfast_HERAxEuclid_lightcone.ipynb` §§6–9 — imports the functions above; defines none of them |
| Tests | `tests/test_uncertainty_budget.py` (16), `tests/test_pipeline.py` (6), `tests/test_figures.py` (1) |

`observational_stage` contains **no physics**. It resolves parameters
(CLI → HDF5 attribute → default) and calls `compute_uncertainty_budget` once.
That is deliberate: the notebook and the HPC run now execute the same code
path, so they cannot drift apart the way they had before this was
consolidated.

---

## 2. The chain, step by step

All numbers below are evaluated at $z_\mathrm{obs} = 7.0$ with the stored
run's cosmology, $H_0 = 67.36$ km s⁻¹ Mpc⁻¹ and $\Omega_{m,0} = 0.315$, giving

$$H(7) = 857.26\ \mathrm{km\,s^{-1}\,Mpc^{-1}}, \qquad D_c(7) = 8821.33\ \mathrm{Mpc}.$$

### 2.1 Photo-$z$ radial smearing (step ①)

A photometric redshift error $\sigma_z$ becomes a comoving line-of-sight
position error:

$$\sigma_r = \frac{c\,\sigma_z}{H(z_\mathrm{obs})}$$

| $\sigma_z$ | $\sigma_r$ | Provenance |
|---|---|---|
| 0.45 | **157.48 Mpc** | current default — Euclid's $\sigma_z/(1+z) < 0.05$ at $z = 7$ |
| 0.059 | 20.65 Mpc | the notebook's **former** value, since corrected (see §4.2) |

> **$\sigma_z$ here is absolute, not $\sigma_z/(1+z)$.** Surveys quote the
> fractional form; multiply by $(1+z)$ before passing it in. This was a real
> error in the original transcription — see `HPC.md` §11.8.

In Fourier space the smearing is a Gaussian window on $k_\parallel$ only
(transverse modes are untouched, and the 21 cm radial coordinate comes from
frequency, which is essentially exact):

$$W(k_\parallel) = \exp\!\left(-\tfrac12 k_\parallel^2 \sigma_r^2\right)$$

$$P_\times^\mathrm{obs} = P_\times\,W, \qquad
P_\mathrm{gal}^\mathrm{obs} = P_\mathrm{gal}\,W^2$$

**One factor of $W$ on the cross-spectrum, two on the galaxy auto-spectrum**,
because only one of the two fields in the cross pair carries photo-z error.
This is the concrete sense in which cross-correlation degrades more gracefully
than the galaxy auto-spectrum.

### 2.2 Foreground wedge (step ②)

Spectrally smooth foregrounds contaminate a wedge below a straight line
through the origin:

$$k_\parallel \le m(z)\,k_\perp + k_\mathrm{buffer}, \qquad
m(z) = \frac{D_c(z)\,H(z)}{c\,(1+z)}$$

| Quantity | Value at $z = 7$ |
|---|---|
| $\lambda_\mathrm{obs} = c(1+z)/\nu_{21}$ | 1.6897 m |
| **Horizon slope $m$** | **3.150906** |
| $\theta_\mathrm{FoV} = \lambda_\mathrm{obs}/D_\mathrm{dish}$ | 6.92° |
| **HERA FoV slope** $\sin\theta_\mathrm{FoV}\,m$ | **0.379360** |
| Buffer | 0.0677 Mpc⁻¹ |
| **Modes surviving** | **97 / 400 (24.2 %)** |

The **horizon** slope defines the mask — that is the configuration Pober et
al.'s buffer was calibrated against. The FoV slope is drawn on figures only
and never masks anything.

The notebook writes the same slope the long way round, through the observed
wavelength:

$$m = \frac{\lambda_\mathrm{obs}\, D_c\, \nu_{21}\, H(z)}{c_\mathrm{m/s}\, c_\mathrm{km/s}\,(1+z)^2}$$

Substituting $\lambda_\mathrm{obs}\nu_{21} = c_\mathrm{m/s}(1+z)$ recovers
$D_c H / [c(1+z)]$ exactly. `test_horizon_slope_equals_the_notebooks_longhand_expression`
asserts the two agree to $10^{-12}$ relative.

### 2.3 Noise (step ③)

**21 cm thermal noise** — a scaling estimate, not an instrument model:

$$T_\mathrm{sys} = T_\mathrm{rcvr} + T_\mathrm{sky}\!\left(\frac{300\ \mathrm{MHz}}{\nu_\mathrm{obs}}\right)^{2.55},
\qquad P_{N,21} = \frac{T_\mathrm{sys}^2 \times 10^3}{t_\mathrm{int}\,\Delta\nu}$$

| Quantity | Value | Constant |
|---|---|---|
| $\nu_\mathrm{obs} = \nu_{21}/(1+z)$ | 177.550625 MHz | |
| $T_\mathrm{rcvr}$ | 100 K | `T_RECEIVER_K` |
| $T_\mathrm{sky}$ at 300 MHz | 60 K | `T_SKY_300MHZ_K` |
| Spectral index | 2.55 | `SKY_SPECTRAL_INDEX` |
| $T_\mathrm{sys}$ | **328.58 K** | |
| $t_\mathrm{int}$ | 3.6 × 10⁶ s (1000 h) | |
| $\Delta\nu$ | 8 × 10⁶ Hz | |
| $P_{N,21}$ | **3.7488 mK² Mpc³** | |

The $10^3$ is `NOISE_NORMALISATION_MPC3`. It is **not a physical constant**:
$T_\mathrm{sys}^2/(t\,\Delta\nu)$ has units of mK², so this factor supplies the
Mpc³ that makes $P_{N,21}$ commensurate with the measured $P_{21}$ — it stands
in for the per-mode survey volume that a full model
($X^2 Y \Omega' / n(k_\perp)$, La Plante Eq. 11) would compute properly. It was
inlined and unnamed in the notebook; naming it makes the approximation visible
rather than accidental. See §7.2.

**Galaxy shot noise:** $P_{N,\mathrm{gal}} = 1/\bar n = 1/(3\times10^{-3}) =
\mathbf{333.33\ Mpc^3}$.

### 2.4 Variance — La Plante et al. (2023), Eqs. 15–17 (step ④)

$$\sigma_{21} = |P_{21}| + P_{N,21}, \qquad
\sigma_\mathrm{gal} = |P_\mathrm{gal}^\mathrm{obs}| + P_{N,\mathrm{gal}}$$

$$\boxed{\;\sigma_\times^2 = \underbrace{\tfrac12 \left(P_\times^\mathrm{obs}\right)^2}_{\text{sample variance}}
\;+\; \underbrace{\tfrac12\,\sigma_{21}\,\sigma_\mathrm{gal}}_{\text{noise coupling}}\;}$$

The two halves are stored separately as `cosmic_variance_term` and
`noise_coupling_term`, and that split is what makes this a *budget* rather
than a single error bar. Their origins differ:

1. **Sample variance** — the estimator scatters because the survey contains a
   finite number of independent modes. More integration time cannot reduce it.
2. **Noise coupling** — thermal and shot noise from each field leak into the
   cross measurement. This *does* fall with integration time and survey density.

`UncertaintyBudget.cosmic_variance_fraction` reports the first term's share of
the total, summed over usable modes: near 1 the measurement is cosmic-variance
limited, near 0 it is noise limited. **For the stored run it is
$2 \times 10^{-224}$** — utterly noise dominated, because photo-z damping has
crushed $P_\times^\mathrm{obs}$ to nothing (§7.1).

### 2.5 Signal-to-noise

$$\mathrm{SNR}(k_\perp, k_\parallel) = \frac{\left|P_\times^\mathrm{obs}\right|}{\sigma_\times},
\qquad
\mathrm{SNR}_\mathrm{total} = \sqrt{\sum_{\text{outside wedge}} \mathrm{SNR}^2}$$

Bins inside the wedge are set to `NaN` and excluded by `np.nansum`; empty
$(k_\perp, k_\parallel)$ bins are already `NaN` from the estimator and drop out
the same way.

#### Why $T_0(z)$ is absent, and why that is correct

La Plante's equations carry the brightness-temperature normalisation $T_0(z)$:
$\sigma_{21} = (P_{21} + P_{N,21})/T_0^2$, the cross-spectrum enters as
$P_\times/T_0$, and $\mathrm{SNR} = |P_\times/T_0| / \sigma_\times$. Neither
the notebook nor this pipeline carries those factors. **They cancel exactly:**

$$\sigma_\times^2 = \tfrac12\left[\left(\frac{P_\times}{T_0}\right)^2 + \frac{(P_{21}+P_{N,21})}{T_0^2}\left(P_\mathrm{gal}+P_{N,\mathrm{gal}}\right)\right]
= \frac{1}{T_0^2}\cdot\tfrac12\left[P_\times^2 + \sigma_{21}\sigma_\mathrm{gal}\right]$$

$$\Rightarrow\quad \mathrm{SNR} = \frac{|P_\times|/T_0}{\frac{1}{T_0}\sqrt{\tfrac12[\cdots]}} = \frac{|P_\times|}{\sqrt{\tfrac12[\cdots]}}$$

so the $T_0$-free form gives an identical SNR. The individual $\sigma$ values
are *not* identical — they differ by $T_0^{-1}$ or $T_0^{-2}$ — so if
$\sigma_\times$ is ever quoted as an absolute error bar rather than used in the
ratio, $T_0(z)$ must be reinstated. This is recorded in the docstring of
`cross_power_snr`.

---

## 3. Provenance audit — notebook vs pipeline

> **The notebook now imports this code rather than duplicating it.** As of
> 2026-08-17, `21cmfast_HERAxEuclid_lightcone.ipynb` calls
> `compute_all_power_spectra` and `compute_uncertainty_budget` from
> `src/analysis.py` directly, so the two implementations are one implementation
> and the audit below is history rather than an ongoing risk. It is retained
> because it is the evidence that the consolidation preserved the physics, and
> because the regression tests still encode it.

Every term of the notebook's chain was transcribed verbatim and compared
against the pipeline's implementation on the notebook's own grid
(256 Mpc × 350.6 Mpc, 128² × 175, 20 × 20 log bins). Result:

| Quantity | Agreement |
|---|---|
| Horizon wedge slope | exact (2 × 10⁻¹⁶ relative) |
| HERA FoV wedge slope | exact (2 × 10⁻¹⁶ relative) |
| $\sigma_r$ | **bit-identical** |
| $W(k_\parallel)$ | **bit-identical** |
| Wedge mask | **identical boolean array** |
| $P_\times^\mathrm{obs}$, $P_\mathrm{gal}^\mathrm{obs}$ | **bit-identical** |
| $\sigma_{21}$, $\sigma_\mathrm{gal}$, $\sigma_\times$ | **bit-identical** |
| Per-mode SNR | **bit-identical** |
| Total SNR | **bit-identical** |
| $P_{N,21}$ | differs by 0.10 % — see §4.1 |

`test_budget_reproduces_the_notebook_chain` locks this end to end: it holds an
independent transcription of the notebook's cells and asserts term-by-term
agreement, so a future edit to `src/analysis.py` that changes the physics will
fail the suite rather than pass silently.

The pipeline additionally reproduces the notebook's **printed** values —
$D_c = 8821$ Mpc, $m = 3.151$, $\theta = 6.9°$, $m_\mathrm{FoV} = 0.379$,
$\sigma_r = 20.6$ Mpc, 24.2 % of modes outside the wedge — asserted in four
separate tests.

Verified once more after the consolidation, by executing the notebook's
rewritten cells against `outputs/lightcone_data.h5` and comparing to
`run_pipeline.observational_stage` on the same fields: $\sigma_r$, both wedge
slopes, both noise powers, the surviving mode count, the per-mode SNR map and
the total ($1.0574485217836499\times10^{-111}$) are **bit-identical**.

---

## 4. The three discrepancies the audit found — and how each was closed

All three have now been fixed **in the notebook**, which no longer carries its
own copy of the calculation. Each subsection records the finding and its
resolution.

### 4.1 The notebook's noise cell used a different 21 cm rest frequency

Its configuration cell defines `F_21_HZ = 1420.405e6` and its wedge cell used
it, but its noise cell hardcoded `1.42e9`:

```python
observed_frequency = 1.42e9 / (1 + z_obs)   # notebook, noise cell
```

Since $T_\mathrm{sky} \propto \nu^{-2.55}$, the 0.03 % frequency shift moves
$P_{N,21}$ by 0.10 %: 3.752581 (notebook) versus **3.748786** (pipeline).

**Resolution:** the pipeline uses `F_21_HZ` consistently everywhere, which is
the correct value, and is recorded in the docstring of
`hera_thermal_noise_power`. **The notebook's hardcoded literal is gone** — its
noise now comes from `compute_uncertainty_budget`, which is passed the
configured `F_21_HZ`.

### 4.2 Corrected $\sigma_z$ and wedge buffer

| Parameter | Formerly | Now (both) | Why |
|---|---|---|---|
| $\sigma_z$ | 0.059 | **0.45** | 0.059 was $\sigma_z/(1+z)$ used as absolute $\sigma_z$; the Euclid requirement at $z = 7$ is $\approx 0.45$ (`HPC.md` §11.8) |
| $k_\mathrm{buffer}$ | 0.02 Mpc⁻¹ | **0.0677 Mpc⁻¹** | 0.1 $h$ Mpc⁻¹ at $h = 0.6766$ — Pober et al. (2014) "moderate", the 21cmSense default |

Both corrections are recorded in `CHANGELOG.md` and patched into the stored
HDF5 attributes. **The notebook's configuration cell now carries the corrected
values too**, with the provenance in comments. Both remain overridable from the
command line, so the former configuration can be reproduced without editing
anything:

```bash
python run_pipeline.py --sigma-z 0.059 --wedge-buffer 0.02
```

which yields $\sigma_r = 20.6$ Mpc and 26.2 % of modes usable, against
157.5 Mpc and 24.2 % at the corrected values.

> The two runs are **not** directly comparable in total SNR: the notebook
> integrates $z = 6.5 \to 7.5$ ($L_\mathrm{LOS} = 350.6$ Mpc, $N_z = 175$)
> while the stored simulation is a $\Delta z = 0.01$ smoke-test slab. The
> $k_\parallel$ grids differ, so the mode sets differ. Only the per-formula
> quantities transfer.

### 4.3 The notebook's stored outputs predated its own configuration

The notebook's printed values —
$D_c = 8821$ Mpc, $m = 3.151$, $\sigma_r = 20.6$ Mpc — were produced with
$H_0 = 67.36$, $\Omega_{m,0} = 0.315$. Its configuration cell had been changed
(uncommitted) to

```python
OMEGA_M_0       = Planck18.Om0        # 0.30966
HUBBLE_CONSTANT = Planck18.H0.value   # 67.66
```

which would give $D_c = 8835$ Mpc, $m = 3.1429$, $\sigma_r = 20.73$ Mpc — but
the notebook had **not been re-executed since that edit**, so its displayed
numbers and its configuration disagreed.

Its `Hz_obs` cell also carried a live bug — `Planck18.H(z).value` referenced
`z`, a stale loop variable, where `z_obs` was meant. That value fed only the
figure overlay; the wedge-mask cell recomputed $H(z_\mathrm{obs})$ correctly,
so no science number depended on it.

**Resolution:** the notebook's configuration cell was reverted to the literal
`67.36` / `0.315` that `run_simulation.py:136-137` uses and writes to the HDF5
root attributes — the values `run_pipeline.py` reads. `Planck18` is retained
for the comoving-distance **endpoints** only, exactly as `run_simulation.py` §1
does. The `Hz_obs` bug is fixed: the cell now calls
`hubble_parameter(z_obs, HUBBLE_CONSTANT, OMEGA_M_0)`.

A useful side effect: because the cosmology went back to the values the
notebook's last real execution used, the stored outputs of the **untouched**
cells (the 21cmFAST run, the halo catalogue, the galaxy field, the Kaiser RSD)
remain valid. Only the seven cells that were rewritten had their outputs
cleared.

> **The regression tests deliberately assert against the former configuration**
> ($\sigma_z = 0.059$, buffer $= 0.02$, $H_0 = 67.36$, $\Omega_{m,0} = 0.315$),
> because those are the parameters the notebook's *published* numbers were
> produced with. They are a fixed historical reference point for the formulas,
> not a claim about the notebook's current configuration — see the module
> docstring of `tests/test_uncertainty_budget.py`.

---

## 5. Parameters

Each value is resolved **CLI flag → HDF5 root attribute → hardcoded default**.
None of them affects the simulated fields, so all can be swept without
`--sim force`.

| Parameter | CLI flag | HDF5 attribute | Default | Units |
|---|---|---|---|---|
| Photo-$z$ uncertainty | `--sigma-z` | `photoz_uncertainty` | 0.45 | absolute $\sigma_z$ |
| Wedge buffer | `--wedge-buffer` | `wedge_buffer` | 0.0677 | Mpc⁻¹ |
| Integration time | `--integration-time` | `integration_time` | 3.6 × 10⁶ | s |
| Bandwidth | `--bandwidth` | `bandwidth` | 8 × 10⁶ | Hz |
| Mean galaxy density | — | `mean_galaxy_density` | 3 × 10⁻³ | see §7.4 |
| Dish diameter | — | `HERA_DISH_DIAMETER` | 14.0 | m |
| 21 cm rest frequency | — | `F_21_HZ` | 1420.405 × 10⁶ | Hz |
| $H_0$, $\Omega_{m,0}$ | — | `HUBBLE_CONSTANT`, `OMEGA_M_0` | 67.36, 0.315 | |

Examples:

```bash
# Reproduce the notebook's configuration
python run_pipeline.py --sigma-z 0.059 --wedge-buffer 0.02

# How much would a deeper HERA campaign help?
python run_pipeline.py --integration-time 1.08e7 --plots none

# Optimistic foregrounds, no buffer at all
python run_pipeline.py --wedge-buffer 0.0 --plots budget
```

---

## 6. Outputs

### 6.1 `outputs/analysis_products.h5` — group `uncertainty_budget`

Appended to the same file that caches the power spectra, so one HDF5 holds
both the raw spectra and everything derived from them. Re-saving replaces the
group, so the budget can be recomputed at a new $\sigma_z$ without touching
the cached spectra.

| Dataset | Shape | Contents |
|---|---|---|
| `photoz_kernel` | (1, n_par) | $W(k_\parallel)$ |
| `P_cross_observed` | (n_perp, n_par) | $P_\times W$ |
| `P_galaxy_observed` | (n_perp, n_par) | $P_\mathrm{gal} W^2$ |
| `outside_wedge` | (n_perp, n_par) | wedge mask, stored `uint8`, restored `bool` |
| `sigma_21cm` | (n_perp, n_par) | $|P_{21}| + P_{N,21}$ |
| `sigma_galaxy` | (n_perp, n_par) | $|P_\mathrm{gal}^\mathrm{obs}| + P_{N,\mathrm{gal}}$ |
| `cosmic_variance_term` | (n_perp, n_par) | $\tfrac12 (P_\times^\mathrm{obs})^2$ |
| `noise_coupling_term` | (n_perp, n_par) | $\tfrac12\sigma_{21}\sigma_\mathrm{gal}$ |
| `sigma_cross` | (n_perp, n_par) | $\sigma_\times$ |
| `snr_per_mode` | (n_perp, n_par) | $|P_\times^\mathrm{obs}|/\sigma_\times$ |

Group attributes carry the 21 scalars of `UncertaintyBudget.as_dict()`.

```python
from src.dataio import load_uncertainty_budget
maps, attrs = load_uncertainty_budget("outputs/analysis_products.h5")
print(attrs["total_snr_sigma"], attrs["cosmic_variance_fraction"])
```

### 6.2 `outputs/pipeline_summary.json` — key `uncertainty_budget`

The same 21 scalars. The former `observation` key is **retained as an alias**
of its eight original fields, so notes and scripts that read
`summary["observation"]["total_snr_sigma"]` keep working.

### 6.3 `outputs/figures/uncertainty_budget.png` — plot group `budget`

Three panels: the damping kernel against the lowest $k_\parallel$ the wedge
admits; $\sigma_\times$ across the $(k_\perp, k_\parallel)$ plane with the
wedge lines overlaid; and the sample-variance share of $\sigma_\times^2$,
showing which term limits each mode.

### 6.4 Console

`observational_stage` prints $T_\mathrm{sys}$, $\sigma_z \to \sigma_r$, $W$ at
the first bin, both wedge slopes, the surviving mode count, both noise powers,
the variance split, and the total SNR.

---

## 6.5 Validation against the cited literature

Checked term by term against the source papers on 2026-08-21. Three formulae
reproduce their references exactly; two documented simplifications are larger
than "simplification" suggests, and are quantified here.

| Term | Cited as | Verdict |
|---|---|---|
| Horizon wedge slope | La Plante Eq. 10 | ✅ **exact** |
| Cross-spectrum variance | La Plante Eq. 15 | ✅ **exact in form** |
| Galaxy variance $\sigma_\mathrm{gal}$ | La Plante Eq. 17 | ✅ **exact** |
| $T_0(z)$ cancellation | La Plante Eq. 16 | ✅ **algebra confirmed** |
| $T_\mathrm{sys}$ model | Pober et al. (2014) | ✅ **standard**; 328.6 K at $z=7$ |
| Wedge buffer 0.0677 Mpc⁻¹ | Pober et al. (2014) "moderate" | ✅ **correct** (0.1 $h$ Mpc⁻¹) |
| Mode weighting | La Plante Eq. 19 | ⚙️ **implemented, opt-in** (`--mode-weighted`) |
| Thermal noise $P_{N,21}$ | La Plante Eq. 11 | ⚙️ **implemented, opt-in** (`--noise-model physical`) |

**Horizon slope (Eq. 10).** The paper writes
$m(z) = \lambda(z) D_c(z) f_{21} H(z) / [c^2 (1+z)^2]$. Substituting
$\lambda(z) = c(1+z)/f_{21}$ gives $m = D_c H / [c(1+z)]$ — precisely
`horizon_wedge_slope`. Not an approximation; the same expression.

**Variance (Eqs. 15–17).** Eq. 15 is
$\sigma^2_{21\times\mathrm{gal}} = \tfrac12[P^2_{21\times\mathrm{gal}} + \sigma_{21}\sigma_\mathrm{gal}]$,
matching `cross_power_snr` exactly. Eq. 16 carries $1/T_0(z)^2$ on the 21 cm
side and Eq. 17 carries none on the galaxy side, exactly as §2.5 assumes; the
cancellation argument there is confirmed.

**Mode weighting (Eqs. 18–20) — the significant one.** Eq. 15 is written *per
mode*; the paper combines bins through
$\hat{s} = \sqrt{N_\mathrm{patch}\,dN}\;P_\times/\sigma_\times$ (Eq. 19) with
$dN = k_\perp^2 k_\parallel V_\mathrm{survey} (2\pi)^{-2}\,d\ln k_\perp\, d\ln k_\parallel$
(Eq. 18). This pipeline computes $P_\times/\sigma_\times$ and stops.

The omitted $dN$ is **already available**: `PowerSpectra.mode_counts` divided
by 2 (the FFT of a real field is Hermitian, so half its cells are redundant)
reproduces Eq. 18 to a median 4.7 % on the fiducial grid — they are the same
quantity. On the stored $128^2\times100$, 256 Mpc geometry at $z=7$, over the
65 usable bins outside the wedge:

| $\sqrt{dN}$ | min 2.0 | median 8.2 | max 64.1 |
|---|---|---|---|

so the reported total SNR is low by roughly **an order of magnitude** (≈18×
if the per-bin SNR were uniform). Directionally safe — the pipeline
under-claims — but too large to leave as a footnote.

**Thermal noise (Eq. 11).** `hera_thermal_noise_power` returns
$T_\mathrm{sys}^2 \times 10^3 / (t_\mathrm{int}\Delta\nu) = 3.75$ mK² Mpc³ at
$z = 7$ for 1000 h. Eq. 11 is
$P^\mathrm{noise}_{21} = T^2_\mathrm{sys}\Omega_p^2 X^2 Y / [\Omega_{pp} t_\mathrm{int} N_\mathrm{pol} N_\mathrm{bl}(u)]$.
Evaluated for HERA at $z = 7$ ($X^2Y = 1227$ Mpc³ sr⁻¹ Hz⁻¹,
$\Omega_p^2/\Omega_{pp} \approx 0.19$ sr, 1000 h, 2 polarisations):

| $N_\mathrm{bl}(u)$ | 50 | 200 | 1000 | 5000 |
|---|---|---|---|---|
| $P_{N,21}$ [mK² Mpc³] | 4.9 × 10⁴ | 1.2 × 10⁴ | 2.4 × 10³ | 4.9 × 10² |

Cross-checked independently against published HERA forecasts
($\Delta^2_N \sim 10$ mK² at $k = 0.2$ for 1000 h $\Rightarrow$
$P_N = 2\pi^2\Delta^2/k^3 \approx 2.5\times10^4$ mK² Mpc³). Both routes land
near $10^4$ mK² Mpc³ against the implemented 3.75 — **the placeholder is
roughly four orders of magnitude too small, and in the optimistic direction.**
It also carries no $k_\perp$ dependence, so it cannot reproduce the steep rise
in noise at high $k_\perp$ where few baselines sample the mode.

Consequence for $\sigma_{21} = |P_{21}| + P_{N,21}$: at 3.75 the term is
negligible beside $P_{21}$, so the 21 cm side of the budget looks
sample-variance limited. With a literature-scale $P_{N,21}$ it dominates by
$\sim$10×, raising $\sigma_\times$ by $\sim$3–4× and lowering the SNR by the
same factor.

**Net effect of the two.** They pull in opposite directions and do *not*
cancel in general — the mode weighting is a property of the survey volume and
binning, the noise of the instrument and integration time.

### Both are now implemented, and both are opt-in

| Flag | Default | Effect |
|---|---|---|
| `--mode-weighted` | off | Applies Eq. 19's $\sqrt{N_\mathrm{patch}\,dN}$ using the estimator's own `mode_counts` |
| `--noise-model physical` | `scaling` | Swaps `hera_thermal_noise_power` for `hera_thermal_noise_power_physical` |

The defaults reproduce every number this pipeline produced before they
existed, so no stored result is silently invalidated. Measured on a
$48^2 \times 100$, 96 Mpc test run at $\sigma_z = 0.02$:

| Configuration | Total SNR | vs default |
|---|---|---|
| default | 0.0500 σ | — |
| `--mode-weighted` | 0.1420 σ | ×2.8 |
| `--noise-model physical` | 0.0056 σ | ÷8.9 |
| both | 0.0159 σ | ÷3.1 |

The physical model also reports which bins the array cannot measure at all
(`inf` noise), which the flat estimate cannot express.

`UncertaintyBudget` records `noise_model` and `mode_weighted`, and both reach
`pipeline_summary.json`, so a stored result always says which of the four
combinations produced it.

---

## 7. Known limitations

### 7.1 The wedge and the photo-$z$ kernel do not overlap

At the smallest $k_\perp = 0.0140$ Mpc⁻¹ the wedge admits only
$k_\parallel > 0.112$ Mpc⁻¹. At $\sigma_z = 0.45$, $W$ there is $5\times10^{-68}$.
**Every mode that survives foreground excision has already been erased by
photo-z smearing**, which is why the total SNR is $10^{-111}$ rather than
merely small. This is a physical statement about the configuration, not a
numerical artefact, and it is the finding the budget exists to surface.
Panel 1 of the budget figure shows it directly. `HPC.md` §5.3 and §11.8 give
the fuller treatment.

### 7.2 The thermal noise is $k$-independent

`hera_thermal_noise_power` returns a scalar, and — as quantified in §6.5 —
one about **10⁴ times smaller** than La Plante Eq. 11 or published HERA
forecasts give. A real interferometer's noise also rises steeply at high
$k_\perp$ where fewer baselines sample the mode. La Plante Eq. 11 resolves
both via $X^2 Y \Omega' / n(k_\perp)$; a worked implementation exists in
`21cm_galaxy_cross_uncertainty.ipynb` but is deliberately **not** ported here — this document's scope is the HERAxEuclid
notebook's budget. For publication forecasts, replace with
[21cmSense](https://github.com/rasg-affiliates/21cmSense). See also
`HPC.md` §11.5.

### 7.3 `mode_counts` is computed but not used

`PowerSpectra.mode_counts` records how many Fourier modes were averaged into
each bin. La Plante's Eqs. 15–17 are written per *single* mode; Eqs. 18–20
combine them with a $\sqrt{N_\mathrm{patch}\,dN}$ weighting. Neither the
notebook nor the pipeline applies that factor, so **the quoted total SNR is a
per-bin quadrature sum and is conservative** — a mode-weighted total would be
larger.

§6.5 quantifies it: `mode_counts / 2` *is* La Plante's $dN$ (agreeing to 4.7 %
on the fiducial grid), and the missing $\sqrt{dN}$ runs from 2.0 to 64.1 with
a median of 8.2 over the usable bins — roughly an order of magnitude in the
total. Tracked in `HPC.md` §11.4.

### 7.4 `mean_galaxy_density` has an unresolved unit

Declared `h³ Mpc⁻³` at `run_simulation.py:126` but consumed as
$P_{N,\mathrm{gal}} = 1/\bar n$ and reported in Mpc³. If the declared $h^3$ is
meant literally, the shot noise is low by $h^{-3} = 3.3\times$. Carried over
from the notebook unchanged; recorded in `HPC.md` §13.2.

### 7.5 A single reference redshift

The whole budget is evaluated at $z_\mathrm{obs}$, the lightcone midpoint.
Harmless for the stored $\Delta z = 0.01$ slab; for the $\Delta z = 1.0$
production range $H(z)$ varies by ~16 % across the box, so $\sigma_r$ and the
wedge slope should evolve along the line of sight. This is part of the P0 work
in `TODO.md`.

---

## 8. Reproducing and verifying

```bash
conda run -n 21cmfast pytest tests/test_uncertainty_budget.py -v   # 16 tests
conda run -n 21cmfast pytest tests/ -q                             # 99 tests

# Recompute the budget from cached spectra — seconds, no simulation
conda run -n 21cmfast python run_pipeline.py --plots budget
```

Last recorded result (stored simulation, default parameters, 2026-08-17):

| Quantity | Value |
|---|---|
| $\sigma_r$ | 157.478 Mpc |
| $W$ at the first $k_\parallel$ bin | 0.02104 |
| Horizon slope / FoV slope | 3.150906 / 0.379360 |
| Modes outside wedge | 97 / 400 (24.25 %) |
| $T_\mathrm{sys}$ | 328.580 K |
| $P_{N,21}$ / $P_{N,\mathrm{gal}}$ | 3.74879 mK² Mpc³ / 333.333 Mpc³ |
| Cosmic-variance fraction | 2.02 × 10⁻²²⁴ |
| **Total SNR** | **1.057 × 10⁻¹¹¹ σ** |

---

## 9. References

- **La Plante, Mirocha, Gorce, Lidz & Parsons (2023)**, ApJ 944, 59 —
  [arXiv:2205.09770](https://arxiv.org/abs/2205.09770), *"Prospects for
  21cm-Galaxy Cross-Correlations with HERA and the Roman High-Latitude
  Survey"* — Eqs. 15–17, the variance formalism; Eq. 10, the wedge slope;
  Eq. 11, the thermal-noise model not ported here; Eqs. 18–20, the mode
  weighting not ported here.
- **Lidz et al. (2009)**, ApJ 690, 252 —
  [arXiv:0806.1055](https://arxiv.org/abs/0806.1055) — original derivation of
  the cross-spectrum variance.
- **Pober et al. (2014)**, ApJ 782, 66 —
  [arXiv:1310.7031](https://arxiv.org/abs/1310.7031) — foreground models and
  the 0.1 $h$ Mpc⁻¹ horizon buffer.
- **Thyagarajan et al. (2015)**, ApJ 804, 14 — foreground wedge geometry.
- **Parsons et al. (2012)**, ApJ 756, 165 —
  [arXiv:1204.4749](https://arxiv.org/abs/1204.4749) — the $X$ and $Y$
  comoving conversion factors.
- **DeBoer et al. (2017)**, PASP 129, 045001 —
  [arXiv:1606.07473](https://arxiv.org/abs/1606.07473) — HERA specifications.
- **Euclid Collaboration (2022)** —
  [arXiv:2108.01201](https://arxiv.org/abs/2108.01201) — photometric redshift
  requirements.
