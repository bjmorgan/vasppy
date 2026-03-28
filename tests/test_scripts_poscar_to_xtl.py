"""Integration tests for the poscar_to_xtl script."""

from pathlib import Path

from vasppy.scripts.poscar_to_xtl import poscar_to_xtl_output


POSCAR_NACL = Path(__file__).parent / "test_data" / "POSCAR_NaCl"


def test_poscar_to_xtl_output_contains_cell_keyword():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    assert "CELL" in output


def test_poscar_to_xtl_output_contains_lattice_parameters():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    # NaCl POSCAR has a cubic cell with a = b = c = 2.82 Å, all angles 90°
    assert "2.82" in output
    assert "90." in output


def test_poscar_to_xtl_output_contains_symmetry_label():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    assert "Symmetry label P1" in output


def test_poscar_to_xtl_output_contains_atoms_header():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    assert "ATOMS" in output
    assert "NAME" in output


def test_poscar_to_xtl_output_contains_species_labels():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    assert "Na" in output
    assert "Cl" in output


def test_poscar_to_xtl_output_contains_fractional_coordinates():
    output = poscar_to_xtl_output(str(POSCAR_NACL))
    # Na is at (0, 0, 0) and Cl is at (0.5, 0.5, 0.5)
    assert "0.00000" in output
    assert "0.50000" in output
