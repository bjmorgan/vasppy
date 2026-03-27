"""Characterisation tests for proc_poscar script."""

import subprocess
import sys
from pathlib import Path

import pytest

POSCAR_NACL = str(Path(__file__).resolve().parents[1] / "test_data" / "POSCAR_NaCl")


def _run_proc_poscar(*args: str) -> str:
    """Run proc_poscar as a subprocess and return its stdout."""
    result = subprocess.run(
        [sys.executable, "-m", "vasppy.scripts.proc_poscar", POSCAR_NACL, *args],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


class TestProcPoscarDefault:
    """Tests for default (direct coordinates) output."""

    def test_default_output_contains_title(self) -> None:
        output = _run_proc_poscar()
        assert output.startswith("NaCl rock salt structure\n")

    def test_default_output_lattice_vectors(self) -> None:
        lines = _run_proc_poscar().splitlines()
        # Scaling factor
        assert lines[1].strip() == "1.0"
        # Lattice vectors (diagonal 2.82 Angstrom cubic cell)
        for i in range(2, 5):
            values = [float(v) for v in lines[i].split()]
            assert len(values) == 3

    def test_default_output_species_and_counts(self) -> None:
        lines = _run_proc_poscar().splitlines()
        assert lines[5].split() == ["Na", "Cl"]
        assert lines[6].split() == ["1", "1"]

    def test_default_output_coordinate_type(self) -> None:
        lines = _run_proc_poscar().splitlines()
        assert lines[7].strip() == "Direct"

    def test_default_output_num_coordinate_lines(self) -> None:
        lines = _run_proc_poscar().splitlines()
        # 2 atoms => 2 coordinate lines after the header (8 header lines)
        coord_lines = lines[8:]
        assert len(coord_lines) == 2


class TestProcPoscarCartesian:
    """Tests for Cartesian coordinate output."""

    def test_cartesian_coordinate_type(self) -> None:
        lines = _run_proc_poscar("-t", "cartesian").splitlines()
        assert lines[7].strip() == "Cartesian"

    def test_cartesian_coordinates_values(self) -> None:
        lines = _run_proc_poscar("-t", "cartesian").splitlines()
        # First atom at origin
        coords_0 = [float(v) for v in lines[8].split()]
        assert coords_0 == pytest.approx([0.0, 0.0, 0.0])
        # Second atom at (1.41, 1.41, 1.41)
        coords_1 = [float(v) for v in lines[9].split()]
        assert coords_1 == pytest.approx([1.41, 1.41, 1.41])


class TestProcPoscarSupercell:
    """Tests for supercell generation."""

    def test_supercell_lattice_doubled(self) -> None:
        lines = _run_proc_poscar("-s", "2", "2", "2").splitlines()
        a_vector = [float(v) for v in lines[2].split()]
        assert a_vector[0] == pytest.approx(5.64)

    def test_supercell_atom_counts(self) -> None:
        lines = _run_proc_poscar("-s", "2", "2", "2").splitlines()
        counts = [int(v) for v in lines[6].split()]
        assert counts == [8, 8]

    def test_supercell_num_coordinate_lines(self) -> None:
        lines = _run_proc_poscar("-s", "2", "2", "2").splitlines()
        coord_lines = lines[8:]
        assert len(coord_lines) == 16


class TestProcPoscarBohr:
    """Tests for Bohr conversion."""

    def test_bohr_lattice_scaled(self) -> None:
        lines = _run_proc_poscar("-b").splitlines()
        a_vector = [float(v) for v in lines[2].split()]
        # 2.82 / 0.529177211 = 5.329...
        assert a_vector[0] == pytest.approx(5.3290276705, abs=1e-6)

    def test_bohr_fractional_coords_unchanged(self) -> None:
        """Fractional coordinates should not change with Bohr conversion."""
        lines = _run_proc_poscar("-b").splitlines()
        coords_0 = [float(v) for v in lines[8].split()]
        assert coords_0 == pytest.approx([0.0, 0.0, 0.0])
        coords_1 = [float(v) for v in lines[9].split()]
        assert coords_1 == pytest.approx([0.5, 0.5, 0.5])


class TestProcPoscarOutputOptions:
    """Tests for various output formatting options."""

    def test_coordinates_only(self) -> None:
        output = _run_proc_poscar("-c")
        lines = output.splitlines()
        # Should only have coordinate lines, no header
        assert len(lines) == 2

    def test_label_position_4(self) -> None:
        lines = _run_proc_poscar("-l", "4").splitlines()
        # Labels should appear as suffix on coordinate lines
        assert lines[8].rstrip().endswith("Na")
        assert lines[9].rstrip().endswith("Cl")

    def test_label_position_1(self) -> None:
        lines = _run_proc_poscar("-l", "1").splitlines()
        # Labels should appear as prefix on coordinate lines
        assert lines[8].lstrip().startswith("Na")
        assert lines[9].lstrip().startswith("Cl")

    def test_numbered_atoms(self) -> None:
        lines = _run_proc_poscar("-n").splitlines()
        assert lines[8].rstrip().endswith("1")
        assert lines[9].rstrip().endswith("2")

    def test_selective_dynamics_T(self) -> None:
        lines = _run_proc_poscar("--selective", "T").splitlines()
        assert "Selective Dynamics" in lines[7]
        assert "T T T" in lines[9]

    def test_selective_dynamics_F(self) -> None:
        lines = _run_proc_poscar("--selective", "F").splitlines()
        assert "Selective Dynamics" in lines[7]
        assert "F F F" in lines[9]
