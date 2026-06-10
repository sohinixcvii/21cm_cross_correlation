#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_simulation.py — Part 1: Simulation & Preliminary Checks
=============================================================
21 cm × Galaxy Cross-Correlation: HERA × Euclid (Lightcone)

This script runs the 21cmFASTv4 lightcone simulation, constructs the
galaxy density field, estimates the galaxy bias, applies Kaiser RSD,
and saves all outputs to an HDF5 file for downstream analysis.

Usage
-----
    # Direct execution (activate conda env first)
    conda activate 21cmfast
    python run_simulation.py

    # SLURM job submission
    sbatch submit_job.sh

Output
------
    outputs/lightcone_data.h5  — all simulation fields and metadata,
                                  read by plot_fields.ipynb and analysis.ipynb

References
----------
    Davies, Mesinger & Murray (2025) — arXiv:2504.17254
    Gagnon-Hartman, Davies & Mesinger (2025) — arXiv:2502.20447
    La Plante et al. (2023) — arXiv:2205.09770
    Euclid Collaboration (2022) — arXiv:2108.01201
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for HPC (no display required)
from scipy.integrate import quad, simpson
import h5py
import warnings

warnings.filterwarnings("ignore")

# ── Check whether 21cmFAST is installed ──────────────────────────────────────
try:
    import py21cmfast as p21c
    from astropy.cosmology import Planck18
    import astropy.units as u
    HAS_21CMFAST = True
    print(f"21cmFAST version: {p21c.__version__}")
except ImportError:
    HAS_21CMFAST = False
    print("21cmFAST not installed — will use synthetic fields for demonstration.")
    print("Install with:  pip install 21cmFAST")

# ── HMF library ──────────────────────────────────────────────────────────────
try:
    from hmf import MassFunction
    HAS_HMF = True
except ImportError:
    HAS_HMF = False
    print("hmf not installed — galaxy bias will use the default configured value.")


# ===========================================================================
#  ★ CONFIGURATION — ALL USER-ADJUSTABLE PARAMETERS ★
#  Edit values in this section only. All subsequent code reads from here.
# ===========================================================================

# ---------------------------------------------------------------------------
#  Simulation grid
# ---------------------------------------------------------------------------
HII_DIM = 128        # cells per side (low-res grid)
BOX_LEN = 256.0      # comoving box side length [Mpc]
DIM     = 3 * HII_DIM  # high-res grid for initial conditions

# ---------------------------------------------------------------------------
#  Lightcone redshift range
# ---------------------------------------------------------------------------
z_min = 6.995          # nearest redshift (low-z end of lightcone)
z_max = 7.005          # farthest redshift (high-z end)

# ---------------------------------------------------------------------------
#  Euclid-like survey parameters
# ---------------------------------------------------------------------------
M_UV_limit          = -18    # UV absolute magnitude cut
photoz_uncertainty  = 0.059  # sigma_z photometric redshift error
mean_galaxy_density = 3e-3   # n_bar  [h^3 Mpc^-3]

# ---------------------------------------------------------------------------
#  Galaxy bias (default; overwritten below if HMF estimation succeeds)
# ---------------------------------------------------------------------------
galaxy_bias = 8              # linear bias b (typical for high-z LBGs)

# ---------------------------------------------------------------------------
#  Cosmological parameters (Planck 2018)
# ---------------------------------------------------------------------------
OMEGA_M_0       = 0.315      # matter density parameter
HUBBLE_CONSTANT = 67.36      # H_0  [km s^-1 Mpc^-1]

# ---------------------------------------------------------------------------
#  Physical constants
# ---------------------------------------------------------------------------
SPEED_OF_LIGHT_KMS = 3e5         # c  [km s^-1]
SPEED_OF_LIGHT_MPS = 3e8         # c  [m s^-1]
F_21_MHZ           = 1420.405    # 21cm rest frequency  [MHz]
F_21_HZ            = F_21_MHZ * 1e6

# ---------------------------------------------------------------------------
#  HERA instrument
# ---------------------------------------------------------------------------
HERA_DISH_DIAMETER = 14.0        # dish diameter  [m]
integration_time   = 1000 * 3600 # total integration  [s]  (1000 h)
bandwidth          = 8e6         # per-band bandwidth  [Hz]  (8 MHz)

# ---------------------------------------------------------------------------
#  Foreground wedge
# ---------------------------------------------------------------------------
wedge_buffer = 0.02              # safety margin beyond horizon line  [Mpc^-1]

# ---------------------------------------------------------------------------
#  Power spectrum binning
# ---------------------------------------------------------------------------
n_bins_perp     = 20
n_bins_parallel = 20

# ---------------------------------------------------------------------------
#  Output path
# ---------------------------------------------------------------------------
OUTPUT_DIR  = "outputs"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "lightcone_data.h5")

# ===========================================================================
#  End of configuration
# ===========================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ===========================================================================
# 1  Derived quantities
# ===========================================================================
# Compute the lightcone geometry from the configuration above: the comoving
# line-of-sight extent L_LOS, the number of LOS slices N_z (matched
# to the transverse cell size), and the redshift array for each slice.

def hubble_parameter(z):
    """H(z) for flat ΛCDM  [km s⁻¹ Mpc⁻¹]."""
    return HUBBLE_CONSTANT * np.sqrt(OMEGA_M_0 * (1 + z)**3 + (1 - OMEGA_M_0))


# ── Comoving distances to lightcone endpoints (using Planck18 for precision) ─
if HAS_21CMFAST:
    D_min = Planck18.comoving_distance(z_min).to(u.Mpc).value
    D_max = Planck18.comoving_distance(z_max).to(u.Mpc).value
else:
    D_min, _ = quad(lambda z_: SPEED_OF_LIGHT_KMS / hubble_parameter(z_), 0, z_min)
    D_max, _ = quad(lambda z_: SPEED_OF_LIGHT_KMS / hubble_parameter(z_), 0, z_max)

L_los     = D_max - D_min   # comoving LOS extent  [Mpc]
cell_size = BOX_LEN / HII_DIM

# ── Reference redshift (centre of lightcone, used for noise / wedge) ─────────
z_obs = 0.5 * (z_min + z_max)

# ── Number of LOS slices ──────────────────────────────────────────────────────
# Normally matched to the transverse cell size; override via minimum_los_slices
# for test runs where the z-range implies fewer natural slices.
minimum_los_slices = 100
N_z = max(int(round(L_los / cell_size)), minimum_los_slices)

# Redshift and distance for each LOS slice (low-z → high-z)
lc_redshifts = np.linspace(z_min, z_max, N_z)

# ── Node redshifts for 21cmFAST (coeval snapshots driving the physics) ────────
# Use ~10 nodes per unit redshift for good accuracy
n_nodes        = max(int(round(10 * (z_max - z_min))), 5)
node_redshifts = np.linspace(z_max, z_min, n_nodes)   # high-z → low-z

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"Box         : {BOX_LEN:.0f} Mpc,  {HII_DIM}³ cells  →  cell size = {cell_size:.1f} Mpc")
print(f"Lightcone   : z = {z_min} → {z_max}   (reference z_obs = {z_obs})")
print(f"LOS extent  : {L_los:.1f} Mpc  →  N_z = {N_z} slices")
print(f"Box shape   : ({HII_DIM}, {HII_DIM}, {N_z})")
print(f"Node redshifts: {n_nodes} nodes from z={z_max} to z={z_min}")
print(f"Euclid      : M_UV < {M_UV_limit},  σ_z = {photoz_uncertainty}")


# ===========================================================================
# 2  Run the 21cmFASTv4 lightcone simulation
# ===========================================================================
# The v4 API uses RectilinearLightconer to specify the LOS slice redshifts and
# the quantities to store, then run_lightcone to execute the simulation.
#
# Compared to a coeval run, the lightcone:
#   - Evolves the physics across the full redshift range [z_min, z_max]
#   - Stores one coeval slice per LOS cell rather than a single snapshot
#   - Produces a rectangular (HII_DIM, HII_DIM, N_z) output volume

if HAS_21CMFAST:
    print("\nRunning 21cmFAST lightcone simulation …")
    print(f"  z = {z_min} → {z_max},  N_z = {N_z} slices")

    inputs = p21c.InputParameters.from_template(
        ["simple"],
        random_seed=42,
    )

    inputs = inputs.clone(
        node_redshifts=node_redshifts,
        simulation_options={
            "HII_DIM": HII_DIM,
            "BOX_LEN": BOX_LEN,
            "DIM":     DIM,
        },
        matter_options={
            "USE_INTERPOLATION_TABLES": "hmf-interpolation",
        },
    )

    lightconer = p21c.RectilinearLightconer(
        lc_redshifts=lc_redshifts,
        quantities=("brightness_temp", "density", "neutral_fraction", "halo_sfr"),
    )

    lightcone = p21c.run_lightcone(
        lightconer=lightconer,
        inputs=inputs,
        include_dvdr_in_tau21=False,
        apply_rsds=False,
    )
    brightness_temp_field = lightcone.lightcones["brightness_temp"]   # (HII_DIM, HII_DIM, N_z)
    density_field         = lightcone.lightcones["density"]
    neutral_fraction      = lightcone.lightcones["neutral_fraction"]

    # Update LOS geometry from the actual simulation output
    N_z         = lightcone.n_slices
    L_los       = lightcone.lightcone_dimensions[2]                  # actual LOS comoving size  [Mpc]
    lc_dist_Mpc = lightcone.lightcone_distances.to(u.Mpc).value      # (N_z,) comoving distances

    print(f"  Lightcone shape : {lightcone.shape}")
    print(f"  LOS extent      : {L_los:.1f} Mpc")
    print(f"  <x_HI>          : {np.mean(neutral_fraction):.3f}")
    print(
        f"  dT_b            : "
        f"[{brightness_temp_field.min():.1f}, "
        f"{brightness_temp_field.max():.1f}] mK"
    )

    # ── Check matter_options (preliminary sanity check) ────────────────────
    print("\nMatter options:")
    print(inputs.matter_options)

else:
    # ==================================================================
    #  Synthetic lightcone fallback (21cmFAST not installed)
    # ==================================================================
    print("\nGenerating synthetic lightcone fields (21cmFAST not available) …")
    np.random.seed(42)

    N      = HII_DIM
    kx     = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
    ky     = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
    kz_los = np.fft.fftfreq(N_z, d=cell_size) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz_los, indexing="ij")
    K_magnitude = np.sqrt(KX**2 + KY**2 + KZ**2)
    K_magnitude[0, 0, 0] = 1.0

    # Matter overdensity from P(k) ∝ k^−2
    ps_amplitude = K_magnitude**(-2)
    ps_amplitude[0, 0, 0] = 0.0
    random_phases = (
        np.random.standard_normal((N, N, N_z))
        + 1j * np.random.standard_normal((N, N, N_z))
    )
    density_field = np.fft.ifftn(np.sqrt(ps_amplitude / 2) * random_phases).real
    density_field = (density_field - density_field.mean()) / density_field.std() * 0.5

    # Neutral fraction via inside-out reionisation
    k_bubble = 2 * np.pi / 10.0
    smoothed = np.fft.ifftn(
        np.fft.fftn(density_field) * np.exp(-0.5 * K_magnitude**2 / k_bubble**2)
    ).real
    neutral_fraction = (smoothed < np.percentile(smoothed, 50)).astype(float)

    # Brightness temperature with redshift evolution along LOS
    T0_per_slice = 27.0 * np.sqrt((1 + lc_redshifts) / 10.0)   # (N_z,)
    brightness_temp_field = (
        T0_per_slice[np.newaxis, np.newaxis, :]
        * neutral_fraction * (1.0 + density_field)
    )

    lc_dist_Mpc = D_min + (D_max - D_min) * np.linspace(0, 1, N_z)

    print(f"  Shape  : {brightness_temp_field.shape}")
    print(f"  ⟨x_HI⟩ : {np.mean(neutral_fraction):.3f}")
    print(f"  δT_b   : [{brightness_temp_field.min():.1f}, {brightness_temp_field.max():.1f}] mK")


# ===========================================================================
# 3a  Extract actual 21cmFASTv4 halo catalogue at z_obs
# ===========================================================================

if HAS_21CMFAST:
    print(f"\nExtracting halo catalogue at z = {z_obs} ...")

    # 1. Compute initial conditions using the same inputs as the lightcone
    initial_conditions = p21c.compute_initial_conditions(
        inputs=inputs,
    )

    # 2. Draw the Lagrangian halo catalogue at z_obs
    halo_catalog = p21c.determine_halo_catalog(
        redshift=z_obs,
        inputs=inputs,
        initial_conditions=initial_conditions,
    )

    # 3. Perturb halos to Eulerian positions at z_obs
    perturbed_halos = p21c.perturb_halo_catalog(
        inputs=inputs,
        initial_conditions=initial_conditions,
        halo_catalog=halo_catalog,
    )

    # 4. Pull useful arrays (convert 21cmFAST Array objects to NumPy)
    halo_masses    = np.asarray(perturbed_halos.halo_masses.value)
    halo_coords    = np.asarray(perturbed_halos.halo_coords.value)
    stellar_masses = np.asarray(perturbed_halos.stellar_masses.value)
    sfr_cat        = np.asarray(perturbed_halos.sfr.value)

    print(f"  Number of halos: {perturbed_halos.n_halos}")

    # ── Preliminary diagnostics ───────────────────────────────────────────────
    print("\nCatalogue diagnostics:")
    print(f"  halo_coords shape = {halo_coords.shape}, dtype = {halo_coords.dtype}")
    print(f"  halo_masses shape = {halo_masses.shape}, dtype = {halo_masses.dtype}")
    print(f"  N_halos           = {halo_coords.shape[0]:,}")
    print(f"  coord min/max     = {halo_coords.min():.3f}, {halo_coords.max():.3f}")
    print(f"  mass min/max      = {halo_masses.min():.3e}, {halo_masses.max():.3e} Msun")

    positive_sfr = sfr_cat[np.isfinite(sfr_cat) & (sfr_cat > 0)]
    print(f"\n  SFR diagnostics:")
    print(f"  N halos with SFR > 0 = {len(positive_sfr):,}")
    if len(positive_sfr) > 0:
        print(f"  SFR min/max = {positive_sfr.min():.3e} - {positive_sfr.max():.3e} Msun/yr")
        print(f"  SFR median  = {np.median(positive_sfr):.3e} Msun/yr")

else:
    # Empty arrays when 21cmFAST is unavailable
    halo_masses    = np.array([])
    halo_coords    = np.array([]).reshape(0, 3)
    stellar_masses = np.array([])
    sfr_cat        = np.array([])


# ===========================================================================
# 3b  Construct the galaxy density field
# ===========================================================================
# The lightcone stores `halo_sfr`: the per-cell SFR density integrated over
# all halos, evaluated at each LOS slice. This naturally traces UV-bright
# galaxies as in the coeval case, but now evolves continuously along the LOS.
#
# Synthetic fallback: if 21cmFAST is unavailable, a biased Poisson-sampled
# galaxy field is generated from the density field, anti-correlated with the
# 21 cm signal.

cell_volume = cell_size**3

if HAS_21CMFAST:
    print("\nConstructing galaxy density field from lightcone halobox …")

    sfr_field = lightcone.lightcones["halo_sfr"]   # (HII_DIM, HII_DIM, N_z)

    mean_sfr = sfr_field.mean()
    if mean_sfr > 0:
        galaxy_overdensity = sfr_field / mean_sfr - 1.0
    else:
        galaxy_overdensity = np.zeros_like(sfr_field)

    print(f"  SFR density  : [{sfr_field.min():.2e}, {sfr_field.max():.2e}] (internal units)")
    print(f"  Galaxy δ     : [{galaxy_overdensity.min():.2f}, {galaxy_overdensity.max():.2f}]")
    print(f"  Shot-noise n̄ : {mean_galaxy_density:.2e} h³ Mpc⁻³  (survey parameter)")

else:
    print("\nGenerating synthetic galaxy field …")

    expected_counts = (
        mean_galaxy_density
        * cell_volume
        * np.exp(galaxy_bias * density_field)
        * (1.0 - neutral_fraction)
    )
    galaxy_counts      = np.random.poisson(np.clip(expected_counts, 0, 50)).astype(float)
    mean_count         = galaxy_counts.mean()
    galaxy_overdensity = galaxy_counts / max(mean_count, 1e-10) - 1.0

    print(f"  Mean galaxy density : {mean_count / cell_volume:.2e} h³ Mpc⁻³")
    print(f"  Total galaxies      : {int(galaxy_counts.sum())}")


# ===========================================================================
# 4  Estimate Euclid galaxy bias from luminosity-selected halo population
# ===========================================================================

if HAS_HMF:
    # --------------------------------------------------------------------------
    # Euclid absolute UV magnitude limits
    # More negative = brighter
    # --------------------------------------------------------------------------
    M_UV_bright = -22
    M_UV_faint  = -18

    def Muv_to_Luv(Muv):
        # L_UV in erg s^-1 Hz^-1
        return 10**((51.63 - Muv) / 2.5)

    L_UV_min = Muv_to_Luv(M_UV_faint)    # faint limit
    L_UV_max = Muv_to_Luv(M_UV_bright)   # bright limit

    print(f"\nL_UV range = {L_UV_min:.3e} - {L_UV_max:.3e} erg s^-1 Hz^-1")

    # --------------------------------------------------------------------------
    # 21cmFAST-like stellar-halo relation
    # --------------------------------------------------------------------------
    F_STAR10   = 0.05
    ALPHA_STAR = 0.5
    M_TURN     = 5e8

    OMEGA_B_0 = 0.049

    def f_star(M_h):
        return F_STAR10 * (M_h / 1e10)**ALPHA_STAR * np.exp(-M_TURN / M_h)

    def stellar_mass_model(M_h):
        return f_star(M_h) * (OMEGA_B_0 / OMEGA_M_0) * M_h

    def sfr_model(M_h):
        # approximate star-formation timescale = 100 Myr
        return stellar_mass_model(M_h) / 1e8

    def uv_luminosity(M_h):
        return 1.15e28 * sfr_model(M_h)

    # --------------------------------------------------------------------------
    # Halo mass function
    # --------------------------------------------------------------------------
    mf = MassFunction(
        z=z_obs,
        Mmin=7,
        Mmax=13,
        dlog10m=0.02,
    )

    M_h       = mf.m
    dndlog10m = mf.dndlog10m

    # --------------------------------------------------------------------------
    # Halo bias: Sheth-Tormen style analytic bias using mf.nu
    # --------------------------------------------------------------------------
    def sheth_tormen_bias_from_nu(nu, delta_c=1.686, a=0.707, p=0.3):
        return 1.0 + (a * nu**2 - 1.0) / delta_c + (
            2.0 * p / (delta_c * (1.0 + (a * nu**2)**p))
        )

    halo_bias = sheth_tormen_bias_from_nu(mf.nu)

    # --------------------------------------------------------------------------
    # Euclid selection
    # --------------------------------------------------------------------------
    L_UV     = uv_luminosity(M_h)
    selected = (L_UV >= L_UV_min) & (L_UV <= L_UV_max)

    if selected.sum() < 2:
        raise RuntimeError(
            "Too few halos satisfy the luminosity selection. "
            "Check M_UV limits or the stellar-halo/SFR parameters."
        )

    # --------------------------------------------------------------------------
    # Effective galaxy number density and bias
    # --------------------------------------------------------------------------
    logM = np.log10(M_h)

    n_gal = simpson(
        dndlog10m[selected],
        x=logM[selected],
    )

    galaxy_bias = simpson(
        dndlog10m[selected] * halo_bias[selected],
        x=logM[selected],
    ) / n_gal

    print(f"  n_gal = {n_gal:.3e} Mpc^-3")
    print(f"  Estimated Euclid galaxy bias = {galaxy_bias:.2f}")
    print(
        f"  Selected halo mass range: "
        f"{M_h[selected].min():.2e} - {M_h[selected].max():.2e} Msun"
    )
else:
    print(f"\nUsing configured galaxy bias = {galaxy_bias}  (hmf not available).")


# ===========================================================================
# 5  Kaiser redshift-space distortions
# ===========================================================================
# The galaxy overdensity is boosted along the LOS by coherent infall
# (Kaiser 1987):
#
#   δ_gal^(s)(k) = (1 + β μ²) δ_gal(k),   β = f/b
#
# The growth rate f and bias b are evaluated at z_obs (midpoint of the
# lightcone). The non-cubic box shape (HII_DIM, HII_DIM, N_z) is accounted
# for in the wavenumber grids. RSD must be applied in Fourier space.

# ── Growth rate and bias at z_obs ─────────────────────────────────────────
Omega_m_z   = OMEGA_M_0 * (1 + z_obs)**3 / (
    OMEGA_M_0 * (1 + z_obs)**3 + (1 - OMEGA_M_0)
)
growth_rate = Omega_m_z ** 0.55           # f(z)
beta_rsd    = growth_rate / galaxy_bias   # β = f / b

print(f"\nReference redshift z_obs = {z_obs}")
print(f"Ω_m(z_obs)  = {Omega_m_z:.4f}")
print(f"Growth rate f  = {growth_rate:.4f}")
print(f"Galaxy bias b  = {galaxy_bias}")
print(f"β = f/b        = {beta_rsd:.4f}")

# ── Wavenumber grids for non-cubic lightcone box ──────────────────────────
#   Transverse: N × N cells, spacing = cell_size [Mpc]
#   LOS       : N_z cells,   spacing = L_los / N_z [Mpc]
N       = HII_DIM
dz_cell = L_los / N_z   # LOS cell size  [Mpc]

kx = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
ky = np.fft.fftfreq(N,   d=cell_size) * 2 * np.pi
kz = np.fft.fftfreq(N_z, d=dz_cell)   * 2 * np.pi
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

k_magnitude = np.sqrt(KX**2 + KY**2 + KZ**2)
k_magnitude[0, 0, 0] = 1.0   # avoid division by zero at DC

# μ = k_∥ / |k|  (LOS along z-axis)
mu = KZ / k_magnitude
mu[0, 0, 0] = 0.0

# Apply Kaiser boost in Fourier space
kaiser_factor          = 1.0 + beta_rsd * mu**2
galaxy_overdensity_k   = np.fft.fftn(galaxy_overdensity)
galaxy_overdensity_rsd = np.fft.ifftn(galaxy_overdensity_k * kaiser_factor).real

print(f"\nBefore RSD  →  δ_gal std = {galaxy_overdensity.std():.3f}")
print(f"After  RSD  →  δ_gal std = {galaxy_overdensity_rsd.std():.3f}")
print(f"Max Kaiser boost (μ=1)   = {1 + beta_rsd:.3f}×")

galaxy_overdensity = galaxy_overdensity_rsd


# ===========================================================================
# 6  Save all outputs to HDF5
# ===========================================================================
# All arrays and metadata are saved for use in the plotting and analysis
# notebooks (Parts 2 and 3). This avoids re-running the expensive simulation.

print(f"\nSaving outputs to {OUTPUT_FILE} …")

with h5py.File(OUTPUT_FILE, "w") as f:

    # ── Simulation fields ─────────────────────────────────────────────────
    f.create_dataset("brightness_temp_field", data=brightness_temp_field,
                     compression="gzip", compression_opts=4)
    f.create_dataset("density_field",         data=density_field,
                     compression="gzip", compression_opts=4)
    f.create_dataset("neutral_fraction",      data=neutral_fraction,
                     compression="gzip", compression_opts=4)
    f.create_dataset("galaxy_overdensity",    data=galaxy_overdensity,
                     compression="gzip", compression_opts=4)

    # ── LOS geometry ─────────────────────────────────────────────────────
    f.create_dataset("lc_redshifts", data=lc_redshifts)
    f.create_dataset("lc_dist_Mpc",  data=lc_dist_Mpc)

    # ── Halo catalogue ────────────────────────────────────────────────────
    hc = f.create_group("halo_catalog")
    hc.create_dataset("halo_masses",    data=halo_masses)
    hc.create_dataset("halo_coords",    data=halo_coords)
    hc.create_dataset("stellar_masses", data=stellar_masses)
    hc.create_dataset("sfr",            data=sfr_cat)

    # ── Scalar metadata (read back as attrs in the notebooks) ─────────────
    f.attrs["HII_DIM"]             = HII_DIM
    f.attrs["BOX_LEN"]             = BOX_LEN
    f.attrs["N_z"]                 = N_z
    f.attrs["L_los"]               = L_los
    f.attrs["cell_size"]           = cell_size
    f.attrs["z_min"]               = z_min
    f.attrs["z_max"]               = z_max
    f.attrs["z_obs"]               = z_obs
    f.attrs["galaxy_bias"]         = galaxy_bias
    f.attrs["beta_rsd"]            = beta_rsd
    f.attrs["mean_galaxy_density"] = mean_galaxy_density
    f.attrs["photoz_uncertainty"]  = photoz_uncertainty
    f.attrs["M_UV_limit"]          = M_UV_limit
    f.attrs["OMEGA_M_0"]           = OMEGA_M_0
    f.attrs["HUBBLE_CONSTANT"]     = HUBBLE_CONSTANT
    f.attrs["SPEED_OF_LIGHT_KMS"]  = SPEED_OF_LIGHT_KMS
    f.attrs["SPEED_OF_LIGHT_MPS"]  = SPEED_OF_LIGHT_MPS
    f.attrs["F_21_MHZ"]            = F_21_MHZ
    f.attrs["F_21_HZ"]             = F_21_HZ
    f.attrs["HERA_DISH_DIAMETER"]  = HERA_DISH_DIAMETER
    f.attrs["integration_time"]    = integration_time
    f.attrs["bandwidth"]           = bandwidth
    f.attrs["wedge_buffer"]        = wedge_buffer
    f.attrs["n_bins_perp"]         = n_bins_perp
    f.attrs["n_bins_parallel"]     = n_bins_parallel

print("  Datasets saved:")
print(f"    brightness_temp_field  : {brightness_temp_field.shape}")
print(f"    density_field          : {density_field.shape}")
print(f"    neutral_fraction       : {neutral_fraction.shape}")
print(f"    galaxy_overdensity     : {galaxy_overdensity.shape}")
print(f"    lc_redshifts           : {lc_redshifts.shape}")
print(f"    lc_dist_Mpc            : {lc_dist_Mpc.shape}")
print(f"    halo_catalog/halo_masses    : {halo_masses.shape}")
print(f"    halo_catalog/halo_coords    : {halo_coords.shape}")
print(f"    halo_catalog/stellar_masses : {stellar_masses.shape}")
print(f"    halo_catalog/sfr            : {sfr_cat.shape}")

print(f"\nSimulation complete. Output written to: {OUTPUT_FILE}")
print("Next steps:")
print("  jupyter notebook notebooks/plot_fields.ipynb   # Part 2: field plots")
print("  jupyter notebook notebooks/analysis.ipynb      # Part 3: power spectra & SNR")
