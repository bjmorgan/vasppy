import unittest
from vasppy.vasprun import Vasprun, parse_varray, parse_structure
from vasppy.vasprun import structure_from_structure_data
from unittest.mock import Mock, patch, call, PropertyMock
from lxml import etree
from io import BytesIO
from pymatgen.core import Lattice, Structure
import numpy as np


def vasprun_from_xml_string(xml_string):
    dummy_xml_file = BytesIO(xml_string.encode("ascii"))
    return Vasprun(dummy_xml_file)


class TestVasprun(unittest.TestCase):
    def test_vasprun_is_initialised(self):
        dummy_xml = "<root>data</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        self.assertEqual(etree.tostring(vasprun.doc).decode("ascii"), dummy_xml)
        for element in vasprun.doc.iter():
            self.assertEqual(element.tag, "root")
            self.assertEqual(element.text, "data")
        self.assertIsNone(vasprun._atom_names)
        self.assertIsNone(vasprun._structures)

    def test_parse_atom_names(self):
        dummy_xml = (
            "<root>\n"
            "  <atominfo>\n"
            "    <array name='atoms'>\n"
            "      <set>\n"
            "        <rc><c>Nb</c><c>   1</c></rc>\n"
            "        <rc><c>O</c><c>    2</c></rc>\n"
            "        <rc><c>F</c><c>    3</c></rc>\n"
            "        <rc><c>F</c><c>    3</c></rc>\n"
            "      </set>\n"
            "    </array>\n"
            "  </atominfo>\n"
            "</root>"
        )
        vasprun = vasprun_from_xml_string(dummy_xml)
        self.assertEqual(vasprun.parse_atom_names(), ["Nb", "O", "F", "F"])

    def test_parse_atom_names_raises_valueerror_if_no_atominfo_found(self):
        dummy_xml = "<root>\n" "  <foo>bar</foo>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        with self.assertRaises(ValueError):
            vasprun.parse_atom_names()

    def test_parse_atom_names_raises_valueerror_if_no_atomname(self):
        dummy_xml = "<root>\n" "  <atominfo>\n" "  </atominfo>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        with self.assertRaises(ValueError):
            vasprun.parse_atom_names()

    def test_parse_varray_parses_type_equals_float(self):
        elem = etree.Element("varray")
        v1 = etree.SubElement(elem, "v")
        v2 = etree.SubElement(elem, "v")
        v1.text = "0.1   0.2   0.3"
        v2.text = "0.4   0.5   0.6"
        expected_data = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]
        self.assertEqual(parse_varray(elem), expected_data)

    def test_parse_varray_parses_type_equals_int(self):
        elem = etree.Element("varray")
        elem.set("type", "int")
        v1 = etree.SubElement(elem, "v")
        v2 = etree.SubElement(elem, "v")
        v1.text = "1   2   3"
        v2.text = "4   5   6"
        expected_data = [[1, 2, 3], [4, 5, 6]]
        self.assertEqual(parse_varray(elem), expected_data)

    def test_parse_varray_parses_type_equals_logical(self):
        elem = etree.Element("varray")
        elem.set("type", "logical")
        v1 = etree.SubElement(elem, "v")
        v2 = etree.SubElement(elem, "v")
        v1.text = "T   T   F"
        v2.text = "F   T   F"
        expected_data = [[True, True, False], [False, True, False]]
        self.assertEqual(parse_varray(elem), expected_data)

    def test_atom_names_calls_parse_atom_names_if_atom_names_is_unset(self):
        dummy_xml = "<root>\n" "  <atominfo>\n" "  </atominfo>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        example_atom_names = ["A", "B", "C"]
        vasprun.parse_atom_names = Mock(return_value=example_atom_names)
        atom_names = vasprun.atom_names
        self.assertEqual(atom_names, example_atom_names)
        vasprun.parse_atom_names.assert_called_once()

    def test_atom_names_returns_cached_value_if_atom_names_is_set(self):
        dummy_xml = "<root>\n" "  <atominfo>\n" "  </atominfo>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        example_atom_names = ["A", "B", "C"]
        vasprun._atom_names = example_atom_names
        vasprun.parse_atom_names = Mock(return_value=example_atom_names)
        atom_names = vasprun.atom_names
        self.assertEqual(atom_names, example_atom_names)
        vasprun.parse_atom_names.assert_not_called()

    def test_structures_calls_parse_structures_if_structures_is_unset(self):
        dummy_xml = "<root>\n" "  <structure>\n" "  </structure>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        example_structures = [Mock(Structure), Mock(Structure)]
        vasprun.parse_structures = Mock(return_value=example_structures)
        structures = vasprun.structures
        self.assertEqual(structures, example_structures)
        vasprun.parse_structures.assert_called_once()

    def test_structure_returns_cached_value_if_structures_is_set(self):
        dummy_xml = "<root>\n" "  <structure>\n" "  </structure>\n" "</root>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        example_structures = [Mock(Structure), Mock(Structure)]
        vasprun.parse_structures = Mock(return_value=example_structures)
        vasprun._structures = example_structures
        structures = vasprun.structures
        self.assertEqual(structures, example_structures)
        vasprun.parse_structures.assert_not_called()

    def test_parse_structure(self):
        # construct a <structure> XML tree
        elem = etree.Element("structure")
        c_elem = etree.SubElement(elem, "crystal")
        v1 = etree.SubElement(c_elem, "varray")
        v11 = etree.SubElement(v1, "v")
        v12 = etree.SubElement(v1, "v")
        v13 = etree.SubElement(v1, "v")
        v11.text = "10.0 0.0 0.0"
        v12.text = "0.0 10.0 0.0"
        v13.text = "0.0 0.0 10.0"
        v2 = etree.SubElement(elem, "varray")
        v2.set("name", "positions")
        v21 = etree.SubElement(v2, "v")
        v22 = etree.SubElement(v2, "v")
        v21.text = "0.1 0.2 0.3"
        v22.text = "0.4 0.5 0.6"
        expected_lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        expected_frac_coords = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]
        with patch("vasppy.vasprun.parse_varray") as mock_parse_varray:
            mock_parse_varray.side_effect = [expected_lattice, expected_frac_coords]
            structure_dict = parse_structure(elem)
        self.assertEqual(structure_dict["lattice"], expected_lattice)
        self.assertEqual(structure_dict["frac_coords"], expected_frac_coords)
        self.assertIsNone(structure_dict["selective_dynamics"])
        mock_parse_varray.assert_has_calls([call(v1), call(v2)])

    # TODO test_parse_structure_with_selective_dynamics
    # Need to generate an example input vasprun.xml with selective dynamics
    # to construct a corresponding etree for running the test.

    def test_structure_from_structure_data(self):
        lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        atom_names = ["Li", "Fe", "O"]
        frac_coords = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
        structure = structure_from_structure_data(
            lattice=lattice, atom_names=atom_names, frac_coords=frac_coords
        )
        self.assertEqual(structure.lattice, Lattice(lattice))
        self.assertEqual([s.species_string for s in structure], atom_names)
        np.testing.assert_array_equal(structure.frac_coords, frac_coords)

    def test_parse_structures(self):
        dummy_xml = "<modeling>\n" "  data" "</modeling>"
        vasprun = vasprun_from_xml_string(dummy_xml)
        c1 = etree.SubElement(vasprun.doc, "calculation")
        c2 = etree.SubElement(vasprun.doc, "calculation")
        s1 = etree.SubElement(c1, "structure")
        s2 = etree.SubElement(c2, "structure")
        expected_lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
        expected_frac_coords = [
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            [[0.11, 0.21, 0.31], [0.41, 0.51, 0.61]],
        ]
        atom_names = ["Li", "Cl"]
        structure_data_1 = {
            "lattice": expected_lattice,
            "frac_coords": expected_frac_coords[0],
            "atom_names": atom_names,
        }
        structure_data_2 = {
            "lattice": expected_lattice,
            "frac_coords": expected_frac_coords[1],
            "atom_names": atom_names,
        }
        with patch(
            "vasppy.vasprun.Vasprun.atom_names", new_callable=PropertyMock
        ) as mock_atom_names:
            mock_atom_names.return_value = atom_names
            with patch("vasppy.vasprun.parse_structure") as mock_parse_structure:
                mock_parse_structure.side_effect = [
                    {
                        "lattice": expected_lattice,
                        "frac_coords": expected_frac_coords[0],
                    },
                    {
                        "lattice": expected_lattice,
                        "frac_coords": expected_frac_coords[1],
                    },
                ]
                with patch(
                    "vasppy.vasprun.structure_from_structure_data"
                ) as mock_structure_from_structure_data:
                    mock_structure_from_structure_data.side_effect = ["foo", "bar"]
                    structures = vasprun.parse_structures()
        self.assertEqual(structures, ["foo", "bar"])
        mock_parse_structure.assert_has_calls([call(s1), call(s2)])
        mock_structure_from_structure_data.assert_has_calls(
            [call(**structure_data_1), call(**structure_data_2)]
        )

    def test_frac_coords(self):
        with patch(
            "vasppy.vasprun.Vasprun.structures", new_callable=PropertyMock
        ) as mock_structures:
            structures = [Mock(spec=Structure), Mock(spec=Structure)]
            frac_coords = [
                [[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                [[0.6, 0.7, 0.8], [0.9, 0.0, 0.1]],
            ]
            structures[0].frac_coords = frac_coords[0]
            structures[1].frac_coords = frac_coords[1]
            mock_structures.return_value = structures
            dummy_xml = "<modeling>\n" "  data" "</modeling>"
            vasprun = vasprun_from_xml_string(dummy_xml)
            np.testing.assert_array_equal(vasprun.frac_coords, np.array(frac_coords))

    def test_cart_coords(self):
        with patch(
            "vasppy.vasprun.Vasprun.structures", new_callable=PropertyMock
        ) as mock_structures:
            structures = [Mock(spec=Structure), Mock(spec=Structure)]
            cart_coords = [
                [[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                [[0.6, 0.7, 0.8], [0.9, 0.0, 0.1]],
            ]
            structures[0].cart_coords = cart_coords[0]
            structures[1].cart_coords = cart_coords[1]
            mock_structures.return_value = structures
            dummy_xml = "<modeling>\n" "  data" "</modeling>"
            vasprun = vasprun_from_xml_string(dummy_xml)
            np.testing.assert_array_equal(vasprun.cart_coords, np.array(cart_coords))

    def test_forces(self):
        dummy_xml = (
            "<modeling>\n"
            "  <calculation>\n"
            "    <varray name='forces'>\n"
            "    </varray>\n"
            "  </calculation>\n"
            "  <calculation>\n"
            "    <varray name='forces'>\n"
            "    </varray>\n"
            "  </calculation>\n"
            "</modeling>"
        )
        expected_forces = [
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            [[0.7, 0.8, 0.9], [1.0, 0.0, 0.1]],
        ]
        vasprun = vasprun_from_xml_string(dummy_xml)
        with patch("vasppy.vasprun.parse_varray") as mock_parse_varray:
            mock_parse_varray.side_effect = expected_forces
            forces = vasprun.forces
        self.assertEqual(mock_parse_varray.call_count, 2)
        np.testing.assert_array_equal(forces, np.array(expected_forces))

    def test_forces_is_none_if_no_forces_in_vasprun(self):
        dummy_xml = (
            "<modeling>\n"
            "  <calculation>\n"
            "  </calculation>\n"
            "  <calculation>\n"
            "  </calculation>\n"
            "</modeling>"
        )
        vasprun = vasprun_from_xml_string(dummy_xml)
        with patch("vasppy.vasprun.parse_varray"):
            forces = vasprun.forces
        self.assertIsNone(forces)


if __name__ == "__main__":
    unittest.main()
