import unittest
import numpy as np

from vasppy.atom import Atom


class AtomTestCase(unittest.TestCase):
    def test_init_atom(self):
        label = "A"
        r = np.array([0.1, 0.2, 0.3])
        atom = Atom(label=label, r=r)
        self.assertEqual(atom.label, label)
        np.testing.assert_array_equal(atom.r, r)


if __name__ == "__main__":
    unittest.main()
