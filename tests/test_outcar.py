import unittest
from unittest.mock import Mock, patch, mock_open

import vasppy.outcar 
import numpy as np

class OutcarTestCase( unittest.TestCase ):

    def test_final_energy_from_outcar( self ):
        example_file = """energy without entropy =    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy  without entropy=    -2997.63294724  energy(sigma->0) =    -2997.63294724\n
                       energy without entropy =    -2997.63289805  energy(sigma->0) =    -2997.63289805\n"""    
        with patch( 'builtins.open', mock_open( read_data=example_file), create=True ) as m:
            self.assertEqual( vasppy.outcar.final_energy_from_outcar(), -2997.63289805 )

if __name__ == '__main__':
    unittest.main()
