import unittest
from unittest.mock import Mock, patch, mock_open, call

from vasppy.outcar import final_energy_from_outcar, potcar_eatom_list_from_outcar

import numpy as np

class OutcarTestCase( unittest.TestCase ):

    def test_final_energy_from_outcar( self ):
        example_file = """energy without entropy =    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy  without entropy=    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy without entropy =    -2997.63289805  energy(sigma->0) =    -2997.63289805\n"""    
        with patch( 'builtins.open', mock_open( read_data=example_file ), create=True ) as m:
            self.assertEqual( final_energy_from_outcar(), -2997.63289805 )

    def test_final_energy_from_outcar_with_filename(self):
        example_file = """energy without entropy =    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy  without entropy=    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy without entropy =    -2997.63289805  energy(sigma->0) =    -2997.63289805\n"""    
        with patch('builtins.open', mock_open(read_data=example_file), create=True) as m:
            final_energy_from_outcar(filename='foo')
            self.assertEqual( m.mock_calls[0], call('foo') )

    def test_potcar_eatom_list_from_potcar( self ):
        example_file = """energy of atom  1       EATOM=-1042.3781\n
                           kinetic energy error for atom=    0.0024 (will be added to EATOM!!)\n
                           energy of atom  2       EATOM= -432.3788\n
                           kinetic energy error for atom=    0.0224 (will be added to EATOM!!)\n
                           energy of atom  3       EATOM= -659.6475\n
                           kinetic energy error for atom=    0.0354 (will be added to EATOM!!)\n"""
        with patch( 'builtins.open', mock_open( read_data=example_file ), create=True ) as m:
            self.assertEqual( potcar_eatom_list_from_outcar(), [ -1042.3781, -432.3788, -659.6475 ] )

if __name__ == '__main__':
    unittest.main()

