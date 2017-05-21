import unittest
from unittest.mock import Mock, patch, mock_open
from collections import Counter

from vasppy.calculation import Calculation, delta_E, delta_stoichiometry, energy_string_to_float, import_calculations_from_file

import numpy as np

class CalculationTestCase( unittest.TestCase ):

    def test_calculation_is_initialised( self ):
        this_title = 'A'
        this_energy = -100.0
        this_stoichiometry = { 'B': 1, 'C': 2 }
        calculation = Calculation( title=this_title, energy=this_energy, stoichiometry=this_stoichiometry )
        self.assertEqual( this_title, calculation.title )
        self.assertEqual( this_energy, calculation.energy )
        self.assertEqual( this_stoichiometry, calculation.stoichiometry )

    def test_mul( self ):
        title = 'A'
        energy = -100.0
        stoichiometry = { 'B': 1, 'C': 2 }
        calculation = Calculation( title=title, energy=energy, stoichiometry=stoichiometry )
        calculation.scale_stoichiometry = Mock( return_value={ 'B': 2, 'C': 4 } )
        new_calculation = calculation * 2
        self.assertEqual( title, new_calculation.title )
        self.assertEqual( energy * 2.0, new_calculation.energy )
        self.assertEqual( { 'B': 2, 'C': 4 }, new_calculation.stoichiometry )
        calculation.scale_stoichiometry.assert_called_with( 2 )

    def test_truediv( self ):
        title = 'A'
        energy = -100.0
        stoichiometry = { 'B': 4, 'C': 2 }
        calculation = Calculation( title=title, energy=energy, stoichiometry=stoichiometry )        
        with patch( 'vasppy.calculation.Calculation.__mul__' ) as mock_mul:
            new_calculation = calculation / 2
            mock_mul.assert_called_with( 0.5 )

    def test_scale_stoichiometry( self ):
        title = 'A'
        energy = -100.0
        stoichiometry = { 'B': 1, 'C': 2 }
        calculation = Calculation( title=title, energy=energy, stoichiometry=stoichiometry )
        self.assertEqual( calculation.scale_stoichiometry( 2 ), { 'B': 2, 'C': 4 } )
     
class CalculationSupportFunctionsTestCase( unittest.TestCase ):

    def test_delta_E( self ):
        titles = [ 'A', 'B', 'C' ]
        energies = [ -50.5, -23.2, -10.1 ]
        stoichiometries = [ { 'B': 1, 'C': 2 }, { 'B': 1, 'C': 1 }, { 'C': 1 } ]
        calculations = [ Calculation( title=t, energy=e, stoichiometry=s ) for t, e, s in zip( titles, energies, stoichiometries ) ]
        self.assertAlmostEqual( delta_E( reactants=[ calculations[0] ], products=calculations[1:3] ), +17.2 ) 

    @patch( 'vasppy.calculation.delta_stoichiometry' )
    def test_delta_E_raises_value_error_if_not_balanced( self, mock_delta_stoichiometry ):
        titles = [ 'A', 'B', 'C' ]
        energies = [ -50.5, -23.2, -10.1 ]
        stoichiometries = [ { 'B': 1, 'C': 2 }, { 'B': 1, 'C': 1 }, { 'C': 2 } ]
        calculations = [ Calculation( title=t, energy=e, stoichiometry=s ) for t, e, s in zip( titles, energies, stoichiometries ) ]
        with self.assertRaises( ValueError ):
            delta_E( reactants=[ calculations[0] ], products=calculations[1:3] )
        mock_delta_stoichiometry.assert_called_with( [ calculations[0] ], calculations[1:3] )

    def test_delta_stoichiometry( self ):
        titles = [ 'A', 'B', 'C' ]
        energies = [ -50.5, -23.2, -10.1 ]
        stoichiometries = [ { 'B': 1, 'C': 2 }, { 'D': 1, 'C': 1 }, { 'C': 1 } ]
        calculations = [ Calculation( title=t, energy=e, stoichiometry=s ) for t, e, s in zip( titles, energies, stoichiometries ) ]
        self.assertEqual( delta_stoichiometry( reactants=[ calculations[0] ], products=calculations[1:3] ), { 'B': -1, 'D': 1 } )

    def test_energy_string_to_float( self ):
        test_strings = { '-1.2345 eV': -1.2345,
                          '0.2341 eV':  0.2341 }
        for k, v in test_strings.items():
            self.assertEqual( energy_string_to_float( k ), v )

    @patch( 'vasppy.calculation.energy_string_to_float' )
    @patch( 'vasppy.calculation.Calculation' )
    def test_import_calculation_from_file( self, mock_Calculation, mock_energy_converter ):
        mock_energy_converter.side_effect = [ -0.2414 ]
        example_yaml = """
        title: this_calculation
        stoichiometry:
            - A: 2
            - B: 4
        energy: -0.2414 eV
        """
        with patch( 'builtins.open', mock_open( read_data=example_yaml ), create=True ) as m:
            import_calculations_from_file( 'example_file' )
            mock_Calculation.assert_called_with( energy=-0.2414, stoichiometry=Counter({'B': 4, 'A': 2}), title='this_calculation' )
            mock_energy_converter.assert_called_with( '-0.2414 eV' )
 
if __name__ == '__main__':
    unittest.main()
