import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from vasppy.doscar import Doscar, pdos_column_names


def _make_doscar_string(
    n_atoms: int = 2,
    n_points: int = 3,
    efermi: float = 5.0,
    include_pdos: bool = True,
) -> str:
    """Build a minimal DOSCAR file as a string.

    The total DOS block has 5 columns: energy up down int_up int_down.
    Each atomic pDOS block has 1 + 9*2 = 19 columns (energy + 9 channels * 2 spins)
    for lmax=2, ispin=2, lorbit=11.
    """
    header_lines = [
        f"  {n_atoms}  {n_atoms}  0  0\n",
        " 0.00000000E+00 0.00000000E+00 0.00000000E+00\n",
        " 1.00000000\n",
        "  CAR\n",
        " unknown system\n",
        f"   -10.000    10.000   {n_points}  {efermi:.4f}  1.0000\n",
    ]
    # Total DOS block
    tdos_lines = []
    for i in range(n_points):
        e = -10.0 + i * 10.0  # -10, 0, 10
        tdos_lines.append(f"  {e:.4f}  {float(i + 1):.4f}  {float(i + 1) * 0.5:.4f}  0.0000  0.0000\n")

    if not include_pdos:
        return "".join(header_lines + tdos_lines)

    # pDOS blocks for each atom
    pdos_blocks = []
    for atom_i in range(n_atoms):
        # Separator line (same format as total DOS header)
        pdos_blocks.append(f"   -10.000    10.000   {n_points}  {efermi:.4f}  1.0000\n")
        for i in range(n_points):
            e = -10.0 + i * 10.0
            # 9 channels * 2 spins = 18 values, each set to (atom_i+1)*0.1*(channel+1)
            vals = []
            for ch in range(9):
                up_val = (atom_i + 1) * 0.01 * (ch + 1)
                down_val = up_val * 0.5
                vals.extend([up_val, down_val])
            vals_str = "  ".join(f"{v:.4f}" for v in vals)
            pdos_blocks.append(f"  {e:.4f}  {vals_str}\n")

    return "".join(header_lines + tdos_lines + pdos_blocks)


def _write_doscar_file(
    n_atoms: int = 2,
    n_points: int = 3,
    include_pdos: bool = True,
) -> str:
    """Write synthetic DOSCAR data to a temporary file and return its path.

    The caller is responsible for deleting the file when finished.
    """
    import tempfile

    content = _make_doscar_string(
        n_atoms=n_atoms,
        n_points=n_points,
        include_pdos=include_pdos,
    )
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".DOSCAR", delete=False)
    tmp.write(content)
    tmp.close()
    return tmp.name


def _create_doscar(
    n_atoms: int = 2,
    n_points: int = 3,
    include_pdos: bool = True,
    read_pdos: bool = True,
    ispin: int = 2,
    lmax: int = 2,
    lorbit: int = 11,
) -> Doscar:
    """Create a Doscar object from synthetic data, using a temporary file.

    Note: the temporary file is deleted after initialisation. Methods that
    re-read the file (e.g. read_atomic_dos_as_df) will not work on the
    returned object. Use _write_doscar_file directly for those tests.
    """
    import os

    path = _write_doscar_file(n_atoms=n_atoms, n_points=n_points, include_pdos=include_pdos)
    try:
        doscar = Doscar(
            path,
            ispin=ispin,
            lmax=lmax,
            lorbit=lorbit,
            read_pdos=read_pdos,
        )
    finally:
        os.unlink(path)
    return doscar


class TestPdosColumnNames(unittest.TestCase):

    def test_lmax2_ispin1(self):
        names = pdos_column_names(lmax=2, ispin=1)
        self.assertEqual(names[0], "energy")
        self.assertEqual(len(names), 10)  # energy + 9

    def test_lmax2_ispin2(self):
        names = pdos_column_names(lmax=2, ispin=2)
        self.assertEqual(names[0], "energy")
        self.assertEqual(len(names), 19)  # energy + 9*2
        self.assertIn("s_up", names)
        self.assertIn("s_down", names)

    def test_lmax3_ispin1(self):
        names = pdos_column_names(lmax=3, ispin=1)
        self.assertEqual(len(names), 17)  # energy + 16

    def test_unsupported_lmax_raises(self):
        with self.assertRaises(ValueError):
            pdos_column_names(lmax=1, ispin=1)


class TestDoscarInit(unittest.TestCase):

    def test_reads_header(self):
        doscar = _create_doscar(read_pdos=False)
        self.assertEqual(doscar.number_of_atoms, 2)
        self.assertEqual(doscar.number_of_data_points, 3)
        self.assertAlmostEqual(doscar.efermi, 5.0)

    def test_spin_orbit_coupling_raises(self):
        import tempfile
        import os
        content = _make_doscar_string()
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".DOSCAR", delete=False)
        tmp.write(content)
        tmp.close()
        try:
            with self.assertRaises(NotImplementedError):
                Doscar(tmp.name, spin_orbit_coupling=True)
        finally:
            os.unlink(tmp.name)


class TestReadTotalDos(unittest.TestCase):

    def test_energy_array_extracted(self):
        doscar = _create_doscar(read_pdos=False)
        self.assertEqual(len(doscar.energy), 3)
        np.testing.assert_allclose(doscar.energy, [-10.0, 0.0, 10.0])

    def test_tdos_does_not_contain_energy_column(self):
        """Bug fix: df.drop('energy') result was not assigned, so
        the energy column was incorrectly retained in self.tdos."""
        doscar = _create_doscar(read_pdos=False)
        self.assertNotIn("energy", doscar.tdos.columns)


class TestReadAtomicDosValidation(unittest.TestCase):

    def test_atom_number_zero_raises_valueerror(self):
        """Bug fix: bitwise & instead of logical and meant atom_number=0
        was not rejected."""
        doscar = _create_doscar()
        with self.assertRaises(ValueError):
            doscar.read_atomic_dos_as_df(0)

    def test_atom_number_too_large_raises_valueerror(self):
        doscar = _create_doscar(n_atoms=2)
        with self.assertRaises(ValueError):
            doscar.read_atomic_dos_as_df(3)

    def test_atom_number_negative_raises_valueerror(self):
        doscar = _create_doscar()
        with self.assertRaises(ValueError):
            doscar.read_atomic_dos_as_df(-1)

    def test_valid_atom_number_succeeds(self):
        import os
        path = _write_doscar_file(n_atoms=2)
        try:
            doscar = Doscar(path, read_pdos=True)
            df = doscar.read_atomic_dos_as_df(1)
            self.assertIsInstance(df, pd.DataFrame)
            self.assertNotIn("energy", df.columns)
        finally:
            os.unlink(path)


class TestReadProjectedDos(unittest.TestCase):

    def test_pdos_shape(self):
        doscar = _create_doscar(n_atoms=2, n_points=3)
        self.assertIsNotNone(doscar.pdos)
        # shape: (n_atoms, n_points, n_channels, ispin)
        self.assertEqual(doscar.pdos.shape, (2, 3, 9, 2))


class TestNoPdosData(unittest.TestCase):

    def test_doscar_without_pdos_read_pdos_false(self):
        """When read_pdos=False, pdos should be None."""
        doscar = _create_doscar(include_pdos=False, read_pdos=False)
        self.assertIsNone(doscar.pdos)

    def test_doscar_without_pdos_read_pdos_true_raises(self):
        """Requesting pDOS when none exists should raise."""
        with self.assertRaises(ValueError):
            _create_doscar(include_pdos=False, read_pdos=True)


class TestPdosSelect(unittest.TestCase):

    def setUp(self):
        self.doscar = _create_doscar(n_atoms=2, n_points=3)

    def test_select_all(self):
        result = self.doscar.pdos_select()
        self.assertEqual(result.shape, (2, 3, 9, 2))

    def test_select_spin_up(self):
        result = self.doscar.pdos_select(spin="up")
        self.assertEqual(result.shape[3], 1)

    def test_select_spin_down(self):
        result = self.doscar.pdos_select(spin="down")
        self.assertEqual(result.shape[3], 1)

    def test_select_spin_both(self):
        result = self.doscar.pdos_select(spin="both")
        self.assertEqual(result.shape[3], 2)

    def test_select_invalid_spin_raises(self):
        with self.assertRaises(ValueError):
            self.doscar.pdos_select(spin="invalid")

    def test_select_spin_raises_for_non_spin_polarised(self):
        doscar = _create_doscar(n_atoms=2, n_points=3, ispin=1)
        with self.assertRaises(ValueError):
            doscar.pdos_select(spin="up")

    def test_select_s_orbital(self):
        result = self.doscar.pdos_select(l="s")
        self.assertEqual(result.shape[2], 1)

    def test_select_p_orbital(self):
        result = self.doscar.pdos_select(l="p")
        self.assertEqual(result.shape[2], 3)

    def test_select_d_orbital(self):
        result = self.doscar.pdos_select(l="d")
        self.assertEqual(result.shape[2], 5)

    def test_select_invalid_l_raises(self):
        with self.assertRaises(ValueError):
            self.doscar.pdos_select(l="g")

    def test_select_p_with_specific_m(self):
        """Bug fix: valid_m_values['p'] ordering must match column order
        so that m-value selection maps to the correct channel indices."""
        # Our test data has different values per channel, so we can
        # verify the correct channel is selected.
        result_x = self.doscar.pdos_select(l="p", m=["x"])
        result_y = self.doscar.pdos_select(l="p", m=["y"])
        result_z = self.doscar.pdos_select(l="p", m=["z"])
        # Each should select exactly one channel
        self.assertEqual(result_x.shape[2], 1)
        self.assertEqual(result_y.shape[2], 1)
        self.assertEqual(result_z.shape[2], 1)
        # The channels should correspond to different data
        # p_y is channel 1, p_z is channel 2, p_x is channel 3
        # In our fixture, channel values scale with (ch+1)
        # So p_y(ch1) has factor 2, p_z(ch2) has factor 3, p_x(ch3) has factor 4
        assert self.doscar.pdos is not None
        expected_y = self.doscar.pdos[:, :, 1:2, :]  # channel 1
        expected_z = self.doscar.pdos[:, :, 2:3, :]  # channel 2
        expected_x = self.doscar.pdos[:, :, 3:4, :]  # channel 3
        np.testing.assert_array_equal(result_y, expected_y)
        np.testing.assert_array_equal(result_z, expected_z)
        np.testing.assert_array_equal(result_x, expected_x)

    def test_select_d_with_specific_m(self):
        result_xy = self.doscar.pdos_select(l="d", m=["xy"])
        assert self.doscar.pdos is not None
        expected = self.doscar.pdos[:, :, 4:5, :]
        np.testing.assert_array_equal(result_xy, expected)

    def test_select_atoms(self):
        result = self.doscar.pdos_select(atoms=[0])
        self.assertEqual(result.shape[0], 1)


class TestDeprecatedPandasApi(unittest.TestCase):

    def test_no_delim_whitespace_warning(self):
        """Ensure no FutureWarning from deprecated delim_whitespace usage."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            _create_doscar(n_atoms=1, n_points=3)


class TestDoscarWithSpecies(unittest.TestCase):
    """Tests for Doscar initialised with a species list."""

    def test_species_is_stored(self):
        import os
        path = _write_doscar_file(n_atoms=2)
        try:
            doscar = Doscar(path, species=["Fe", "O"])
            self.assertEqual(doscar.species, ["Fe", "O"])
        finally:
            os.unlink(path)

    def test_ispin1_stored(self):
        doscar = _create_doscar(ispin=1, read_pdos=False)
        self.assertEqual(doscar.ispin, 1)

    def test_lmax_stored(self):
        doscar = _create_doscar(lmax=2, read_pdos=False)
        self.assertEqual(doscar.lmax, 2)

    def test_lorbit_stored(self):
        doscar = _create_doscar(lorbit=11, read_pdos=False)
        self.assertEqual(doscar.lorbit, 11)


class TestNumberOfChannels(unittest.TestCase):

    def test_lmax2_lorbit11(self):
        doscar = _create_doscar(lmax=2, read_pdos=False)
        self.assertEqual(doscar.number_of_channels, 9)

    def test_lmax3_lorbit11(self):
        import os
        # Build a minimal f-states DOSCAR: 17 columns per energy point (energy + 16 channels * 2 spins)
        # Use the same helper but override lmax when creating the Doscar
        path = _write_doscar_file(n_atoms=1, n_points=3, include_pdos=False)
        try:
            doscar = Doscar(path, lmax=3, read_pdos=False)
            self.assertEqual(doscar.number_of_channels, 16)
        finally:
            os.unlink(path)

    def test_unsupported_lorbit_raises(self):
        doscar = _create_doscar(lorbit=11, read_pdos=False)
        doscar.lorbit = 10
        with self.assertRaises(NotImplementedError):
            _ = doscar.number_of_channels


class TestProcessHeader(unittest.TestCase):

    def test_process_header_sets_efermi(self):
        doscar = _create_doscar(read_pdos=False)
        self.assertAlmostEqual(doscar.efermi, 5.0)

    def test_process_header_sets_number_of_data_points(self):
        doscar = _create_doscar(n_points=5, read_pdos=False)
        self.assertEqual(doscar.number_of_data_points, 5)


class TestPdosSelectNoPdos(unittest.TestCase):
    """Test that pdos_select raises TypeError when pdos is not available."""

    def test_raises_typeerror_when_pdos_is_none(self):
        doscar = _create_doscar(read_pdos=False)
        with self.assertRaises(TypeError):
            doscar.pdos_select()


class TestPdosSelectAtomsList(unittest.TestCase):

    def test_atom_as_int_selects_single_atom(self):
        doscar = _create_doscar(n_atoms=2)
        result = doscar.pdos_select(atoms=1)
        self.assertEqual(result.shape[0], 1)


class TestPdosSelectFOrbitals(unittest.TestCase):
    """Tests for f-orbital selection — requires lmax=3 pDOS data."""

    def _make_f_doscar_string(self) -> str:
        """Build a minimal DOSCAR with lmax=3, ispin=2 pDOS data."""
        n_atoms, n_points, efermi = 1, 3, 0.0
        header_lines = [
            f"  {n_atoms}  {n_atoms}  0  0\n",
            " 0.00000000E+00 0.00000000E+00 0.00000000E+00\n",
            " 1.00000000\n",
            "  CAR\n",
            " unknown system\n",
            f"   -10.000    10.000   {n_points}  {efermi:.4f}  1.0000\n",
        ]
        tdos_lines = []
        for i in range(n_points):
            e = -10.0 + i * 10.0
            tdos_lines.append(f"  {e:.4f}  {float(i+1):.4f}  {float(i+1)*0.5:.4f}  0.0  0.0\n")
        pdos_blocks = [f"   -10.000    10.000   {n_points}  {efermi:.4f}  1.0000\n"]
        for i in range(n_points):
            e = -10.0 + i * 10.0
            # 16 channels * 2 spins = 32 values
            vals = [f"{(ch + 1) * 0.01:.4f}" for ch in range(32)]
            pdos_blocks.append(f"  {e:.4f}  " + "  ".join(vals) + "\n")
        return "".join(header_lines + tdos_lines + pdos_blocks)

    def _create_f_doscar(self) -> Doscar:
        import os
        import tempfile
        content = self._make_f_doscar_string()
        tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".DOSCAR", delete=False)
        tmp.write(content)
        tmp.close()
        try:
            doscar = Doscar(tmp.name, lmax=3, ispin=2, read_pdos=True)
        finally:
            os.unlink(tmp.name)
        return doscar

    def test_f_orbital_all_channels(self):
        doscar = self._create_f_doscar()
        result = doscar.pdos_select(l="f")
        self.assertEqual(result.shape[2], 7)

    def test_f_orbital_specific_m(self):
        doscar = self._create_f_doscar()
        result = doscar.pdos_select(l="f", m=["xyz"])
        self.assertEqual(result.shape[2], 1)


class TestPdosSum(unittest.TestCase):

    def setUp(self):
        self.doscar = _create_doscar(n_atoms=2, n_points=3)

    def test_pdos_sum_returns_1d_array(self):
        result = self.doscar.pdos_sum()
        self.assertEqual(result.ndim, 1)
        self.assertEqual(result.shape[0], 3)

    def test_pdos_sum_spin_up(self):
        result = self.doscar.pdos_sum(spin="up")
        self.assertEqual(result.ndim, 1)

    def test_pdos_sum_single_atom(self):
        result_all = self.doscar.pdos_sum()
        result_atom0 = self.doscar.pdos_sum(atoms=[0])
        # Sum over all atoms should be at least as large as one atom
        self.assertTrue(np.all(result_all >= result_atom0 - 1e-10))


class TestPlotPdos(unittest.TestCase):

    def setUp(self):
        path = _write_doscar_file(n_atoms=2)
        self.doscar = Doscar(path, species=["Fe", "O"])
        self._path = path
        self._remove_after = True

    def tearDown(self):
        import os
        if self._remove_after:
            try:
                os.unlink(self._path)
            except FileNotFoundError:
                pass

    def test_plot_pdos_returns_figure(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos()
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_with_axes_returns_none(self):
        import matplotlib
        matplotlib.use("Agg")
        fig_ext, ax = plt.subplots()
        result = self.doscar.plot_pdos(ax=ax)
        self.assertIsNone(result)
        plt.close("all")

    def test_plot_pdos_with_xrange(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(xrange=(-5.0, 5.0))
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_with_to_plot(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(to_plot={"Fe": ["s"], "O": ["s"]})
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_with_scaling(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(
            to_plot={"Fe": ["s"], "O": ["s"]},
            scaling={"Fe": {"s": 2.0}},
        )
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_legend_inside(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(legend_pos="upper right")
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_no_total_dos(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(plot_total_dos=False)
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_with_title(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(title="Test Title")
        self.assertIsNotNone(fig)
        plt.close("all")

    def test_plot_pdos_with_ymax(self):
        import matplotlib
        matplotlib.use("Agg")
        fig = self.doscar.plot_pdos(ymax=10.0)
        self.assertIsNotNone(fig)
        plt.close("all")


if __name__ == "__main__":
    unittest.main()
