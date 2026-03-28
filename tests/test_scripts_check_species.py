"""Integration tests for the check_species script."""

from pathlib import Path

from vasppy.scripts.check_species import unique_species_from_structure


POSCAR_NACL = Path(__file__).parent / "test_data" / "POSCAR_NaCl"


def test_unique_species_from_nacl_poscar():
    species = unique_species_from_structure(str(POSCAR_NACL))
    assert species == ["Na", "Cl"]


def test_unique_species_returns_list():
    species = unique_species_from_structure(str(POSCAR_NACL))
    assert isinstance(species, list)


def test_unique_species_preserves_order():
    """Species should appear in the order they first appear in the structure."""
    species = unique_species_from_structure(str(POSCAR_NACL))
    assert species[0] == "Na"
    assert species[1] == "Cl"
