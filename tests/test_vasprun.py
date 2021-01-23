import unittest
from vasppy.vasprun import Vasprun
from unittest.mock import Mock, patch, call, mock_open
from lxml import etree
from io import BytesIO

class TestVasprun(unittest.TestCase):

    def test_vasprun_is_initialised(self):
        dummy_xml = "<root>data</root>" 
        dummy_xml_file = BytesIO(dummy_xml.encode('ascii'))
        vasprun = Vasprun(dummy_xml_file)
        self.assertEqual(etree.tostring(vasprun.doc).decode('ascii'), dummy_xml)
        for element in vasprun.doc.iter():
            self.assertEqual(element.tag, 'root')
            self.assertEqual(element.text, 'data')

    def test_species(self):
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
        dummy_xml_file = BytesIO(dummy_xml.encode('ascii'))
        vasprun = Vasprun(dummy_xml_file)
        self.assertEqual(vasprun.species(), ['Nb', 'O', 'F', 'F'])

    def test_species_raises_valueerror_if_no_atominfo_found(self):
        dummy_xml = ("<root>\n"
                     "  <foo>bar</foo>\n"
                     "</root>")
        dummy_xml_file = BytesIO(dummy_xml.encode('ascii'))
        vasprun = Vasprun(dummy_xml_file)
        with self.assertRaises(ValueError):
            vasprun.species()
   
    def test_species_raises_valueerror_if_no_atomname(self):
        dummy_xml = ("<root>\n"
                     "  <atominfo>\n"
                     "  </atominfo>\n"
                     "</root>")
        dummy_xml_file = BytesIO(dummy_xml.encode('ascii'))
        vasprun = Vasprun(dummy_xml_file)
        with self.assertRaises(ValueError):
            vasprun.species()
   
if __name__ == '__main__':
    unittest.main()
