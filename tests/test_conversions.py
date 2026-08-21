#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``src/conversions.py``."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.cosmology import Planck18

from src.conversions import (
    SimulationBox,
    cell_mass,
    mean_matter_density,
    survey_area_to_box_size,
)


def test_mean_matter_density_matches_astropy() -> None:
    """rho_m = Omega_m * rho_crit,0 agrees with astropy to 0.1%."""
    import astropy.units as u

    expected = (
        Planck18.Om0 * Planck18.critical_density0.to(u.Msun / u.Mpc**3).value
    )
    assert mean_matter_density() == pytest.approx(expected, rel=1e-3)


def test_mean_matter_density_scales_with_h_squared() -> None:
    """rho_crit,0 ∝ H_0², so doubling H_0 quadruples the density."""
    rho_1 = mean_matter_density(omega_m=0.3, hubble_constant=70.0)
    rho_2 = mean_matter_density(omega_m=0.3, hubble_constant=140.0)
    assert rho_2 == pytest.approx(4.0 * rho_1, rel=1e-12)


def test_cell_mass_scales_with_volume() -> None:
    """M_cell ∝ L_cell³."""
    assert cell_mass(2.0) == pytest.approx(8.0 * cell_mass(1.0), rel=1e-12)


def test_cell_mass_production_grid() -> None:
    """The 256 Mpc / 128³ production grid gives ~3.1e11 M_sun per cell."""
    mass = cell_mass(256.0 / 128, omega_m=0.315, hubble_constant=67.36)
    assert mass == pytest.approx(3.17e11, rel=1e-2)


def test_cell_mass_hires_grid_is_27x_finer() -> None:
    """DIM = 3 * HII_DIM makes the IC cell mass 27x smaller."""
    lores = cell_mass(256.0 / 128, omega_m=0.315, hubble_constant=67.36)
    hires = cell_mass(256.0 / 384, omega_m=0.315, hubble_constant=67.36)
    assert lores / hires == pytest.approx(27.0, rel=1e-12)


def test_box_mass_equals_ncells_times_cell_mass() -> None:
    """Summing all cells recovers the total mass in the box."""
    box_len, n_cells = 256.0, 128
    total = cell_mass(box_len)
    assert n_cells**3 * cell_mass(box_len / n_cells) == pytest.approx(
        total, rel=1e-12
    )


def test_cell_mass_accepts_array_input() -> None:
    """Vectorised over cell sizes."""
    sizes = np.array([1.0, 2.0, 4.0])
    masses = cell_mass(sizes)
    assert masses.shape == sizes.shape
    assert np.all(np.diff(masses) > 0)


# ---------------------------------------------------------------------------
#  survey_area_to_box_size
# ---------------------------------------------------------------------------

# Euclid Deep Field Fornax at z = 7, delta_z = 2 * sigma_z = 0.90 (+/-1 sigma).
FORNAX_AREA_DEG2 = 10.0
FORNAX_Z_CENTRAL = 7.0
FORNAX_DELTA_Z = 0.90


def test_survey_area_to_box_size_fornax_transverse() -> None:
    """10 deg^2 at z = 7 gives a ~486 Mpc transverse box under Planck18."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert box.box_len == pytest.approx(486.33, rel=1e-3)


def test_survey_area_to_box_size_matches_small_angle_formula() -> None:
    """L_perp = sqrt(Omega) * D_M(z), computed independently from astropy."""
    import astropy.units as u

    omega_sr = (FORNAX_AREA_DEG2 * u.deg**2).to(u.sr).value
    d_m = Planck18.comoving_transverse_distance(FORNAX_Z_CENTRAL).to(u.Mpc).value
    expected = np.sqrt(omega_sr) * d_m

    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert box.box_len == pytest.approx(expected, rel=1e-12)


def test_survey_area_to_box_size_los_uses_distance_differencing() -> None:
    """L_los = D_C(z + dz/2) - D_C(z - dz/2), the codebase's convention."""
    import astropy.units as u

    z_lo = FORNAX_Z_CENTRAL - FORNAX_DELTA_Z / 2
    z_hi = FORNAX_Z_CENTRAL + FORNAX_DELTA_Z / 2
    expected = (
        Planck18.comoving_distance(z_hi) - Planck18.comoving_distance(z_lo)
    ).to(u.Mpc).value

    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert box.los_depth == pytest.approx(expected, rel=1e-12)
    assert box.los_depth == pytest.approx(315.60, rel=1e-3)
    assert (box.z_min, box.z_max) == pytest.approx((6.55, 7.45))


def test_survey_area_to_box_size_preserves_mass_resolution() -> None:
    """The derived grid keeps M_cell within a factor 1.2 of the old 2 Mpc grid."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert box.hii_dim == 256
    assert box.dim == 768
    assert box.cell_size == pytest.approx(1.90, rel=1e-2)

    old_grid = cell_mass(256.0 / 128, omega_m=0.315, hubble_constant=67.36)
    new_grid = cell_mass(box.cell_size, omega_m=0.315, hubble_constant=67.36)
    assert new_grid / old_grid == pytest.approx(1.0, abs=0.2)


def test_survey_area_to_box_size_hii_dim_snaps_to_power_of_two() -> None:
    """244 cells at 2 Mpc snap up to 256; unsnapped keeps the exact ceil."""
    kwargs = dict(
        area_deg2=FORNAX_AREA_DEG2,
        z_central=FORNAX_Z_CENTRAL,
        delta_z=FORNAX_DELTA_Z,
        target_cell_size_mpc=2.0,
    )
    assert survey_area_to_box_size(**kwargs).hii_dim == 256
    assert (
        survey_area_to_box_size(
            **kwargs, snap_hii_dim_to_power_of_two=False
        ).hii_dim
        == 244
    )


def test_survey_area_to_box_size_explicit_hii_dim_is_respected() -> None:
    """Pinning hii_dim overrides the target cell size."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z, hii_dim=128
    )
    assert box.hii_dim == 128
    assert box.dim == 384
    assert box.cell_size == pytest.approx(box.box_len / 128, rel=1e-12)


def test_survey_area_to_box_size_simulation_options_mapping() -> None:
    """The 21cmFAST kwargs come out under the names clone() expects."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    options = box.simulation_options
    assert set(options) == {"HII_DIM", "BOX_LEN", "DIM"}
    assert options["HII_DIM"] == box.hii_dim
    assert options["BOX_LEN"] == box.box_len
    assert options["DIM"] == box.dim
    assert isinstance(options["BOX_LEN"], float)
    assert isinstance(options["HII_DIM"], int)


def test_survey_area_to_box_size_scales_as_sqrt_area() -> None:
    """L_perp is proportional to sqrt(area), so 4x area doubles the box."""
    small = survey_area_to_box_size(2.5, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z)
    large = survey_area_to_box_size(10.0, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z)
    assert large.box_len == pytest.approx(2.0 * small.box_len, rel=1e-12)


def test_survey_area_to_box_size_fornax_needs_no_los_tiling() -> None:
    """L_los < BOX_LEN, so the coeval box is not repeated along the LOS."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert box.n_los_tiles == pytest.approx(box.los_depth / box.box_len)
    assert box.n_los_tiles < 1.0


def test_survey_area_to_box_size_returns_dataclass() -> None:
    """The return type is the documented SimulationBox."""
    box = survey_area_to_box_size(
        FORNAX_AREA_DEG2, FORNAX_Z_CENTRAL, FORNAX_DELTA_Z
    )
    assert isinstance(box, SimulationBox)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(area_deg2=0.0, z_central=7.0, delta_z=0.9),
        dict(area_deg2=-1.0, z_central=7.0, delta_z=0.9),
        dict(area_deg2=10.0, z_central=7.0, delta_z=0.0),
        dict(area_deg2=10.0, z_central=7.0, delta_z=-0.9),
        dict(area_deg2=10.0, z_central=0.2, delta_z=0.9),   # reaches z <= 0
        dict(area_deg2=10.0, z_central=7.0, delta_z=0.9, target_cell_size_mpc=0),
        dict(area_deg2=10.0, z_central=7.0, delta_z=0.9, hii_dim=0),
    ],
)
def test_survey_area_to_box_size_rejects_bad_input(kwargs) -> None:
    """Non-physical inputs raise rather than returning a silent nonsense box."""
    with pytest.raises(ValueError):
        survey_area_to_box_size(**kwargs)
