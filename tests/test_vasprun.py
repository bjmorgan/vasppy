import unittest
from vasppy.vasprun import Vasprun, parse_varray
from unittest.mock import Mock, patch, call, mock_open
from lxml import etree
from io import BytesIO

def vasprun_from_xml_string(xml_string):
    dummy_xml_file = BytesIO(xml_string.encode('ascii'))
    return Vasprun(dummy_xml_file)

class TestVasprun(unittest.TestCase):

    def test_vasprun_is_initialised(self):
        dummy_xml = "<root>data</root>" 
        vasprun = vasprun_from_xml_string(dummy_xml)
        self.assertEqual(etree.tostring(vasprun.doc).decode('ascii'), dummy_xml)
        for element in vasprun.doc.iter():
            self.assertEqual(element.tag, 'root')
            self.assertEqual(element.text, 'data')

    def test_atom_names(self):
        dummy_xml = ("<root>\n"
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
                     "</root>")
        vasprun = vasprun_from_xml_string(dummy_xml)
        self.assertEqual(vasprun.atom_names(), ['Nb', 'O', 'F', 'F'])

    def test_atom_names_raises_valueerror_if_no_atominfo_found(self):
        dummy_xml = ("<root>\n"
                     "  <foo>bar</foo>\n"
                     "</root>")
        vasprun = vasprun_from_xml_string(dummy_xml)
        with self.assertRaises(ValueError):
            vasprun.atom_names()
   
    def test_atom_names_raises_valueerror_if_no_atomname(self):
        dummy_xml = ("<root>\n"
                     "  <atominfo>\n"
                     "  </atominfo>\n"
                     "</root>")
        vasprun = vasprun_from_xml_string(dummy_xml)
        with self.assertRaises(ValueError):
            vasprun.atom_names()
  
    def test_parse_varray_parses_type_equals_float(self):
        elem = etree.Element("varray")
        v1 = etree.SubElement(elem, "v")
        v2 = etree.SubElement(elem, "v")
        v1.text = "0.1   0.2   0.3"
        v2.text = "0.4   0.5   0.6"
        expected_data = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]
        self.assertEqual(parse_varray(elem), expected_data)

    def test_parse_varray_parses_type_equals_float(self):
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


if __name__ == '__main__':
    unittest.main()
