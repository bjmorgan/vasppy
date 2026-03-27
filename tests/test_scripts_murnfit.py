"""Tests for vasppy.scripts.murnfit module."""

import numpy as np
import numpy.testing as npt
import pytest

from vasppy.scripts.murnfit import murnaghan, objective, lstsq_fit


class TestMurnaghan:
    """Tests for the Murnaghan equation of state function."""

    def test_at_equilibrium_volume(self) -> None:
        """Energy at V0 should equal E0."""
        e0, b0, bp, v0 = -10.0, 100.0, 4.0, 50.0
        result = murnaghan(v0, e0, b0, bp, v0)
        assert result == pytest.approx(e0)

    def test_known_values(self) -> None:
        """Murnaghan energy for a known set of parameters."""
        e0, b0, bp, v0 = -8.0, 50.0, 3.5, 40.0
        vol = 42.0
        result = murnaghan(vol, e0, b0, bp, v0)
        # Energy should be higher than E0 away from equilibrium
        assert result > e0

    def test_symmetric_deviation(self) -> None:
        """Energy increases for both expansion and compression."""
        e0, b0, bp, v0 = -10.0, 100.0, 4.0, 50.0
        e_compressed = murnaghan(45.0, e0, b0, bp, v0)
        e_expanded = murnaghan(55.0, e0, b0, bp, v0)
        assert e_compressed > e0
        assert e_expanded > e0

    def test_vectorised(self) -> None:
        """Function should work with numpy arrays."""
        e0, b0, bp, v0 = -10.0, 100.0, 4.0, 50.0
        volumes = np.array([45.0, 50.0, 55.0])
        energies = murnaghan(volumes, e0, b0, bp, v0)
        assert energies.shape == (3,)
        assert energies[1] == pytest.approx(e0)


class TestObjective:
    """Tests for the least-squares objective function."""

    def test_zero_residual_at_exact_fit(self) -> None:
        """Residuals should be zero when data matches the model exactly."""
        pars = (-10.0, 100.0, 4.0, 50.0)
        volumes = np.array([45.0, 48.0, 50.0, 52.0, 55.0])
        energies = murnaghan(volumes, *pars)
        residuals = objective(pars, volumes, energies)
        npt.assert_array_almost_equal(residuals, np.zeros(5))


class TestLstsqFit:
    """Tests for the least-squares fitting wrapper."""

    def test_recovers_known_parameters(self) -> None:
        """Fit should recover parameters from synthetic Murnaghan data."""
        true_pars = (-10.0, 100.0, 4.0, 50.0)
        volumes = np.linspace(44.0, 56.0, 20)
        energies = murnaghan(volumes, *true_pars)
        fitted, _ = lstsq_fit(volumes, energies)
        assert fitted[0] == pytest.approx(true_pars[0], abs=1e-4)  # e0
        assert fitted[3] == pytest.approx(true_pars[3], abs=1e-4)  # v0
