import numpy as np
import sys
import copy

class Doscar:

    number_of_header_lines = 6

    def __init__( self, filename ):
        self.filename = filename
        self.read_header()

    def read_header( self ):
        self.header = []
        with open( self.filename, 'r' ) as file_in:
            for i in range( Doscar.number_of_header_lines ):
                self.header.append( file_in.readline() )
        self.process_header()
 
    def process_header( self ):
        self.number_of_data_points = int( self.header[5].split()[2] )
        self.fermi_energy = float( self.header[5].split()[3] )

    def read_total_dos( self ): # assumes spin_polarised
        start_to_read = Doscar.number_of_header_lines
        stop_reading  = start_to_read + self.number_of_data_points
        data_array = self.read_lines_to_numpy_array( start_to_read, stop_reading )
        return Total_DOS( data = data_array, spin_polarised = True ) 

    def read_atomic_dos( self, atom_number, align_fermi_energy = True ): # currently assume spin-polarised, no-SO-coupling, no f-states
        assert atom_number > 0 
        start_to_read = Doscar.number_of_header_lines + atom_number * ( self.number_of_data_points + 1 )
        stop_reading  = start_to_read + self.number_of_data_points 
        data_array = self.read_lines_to_numpy_array( start_to_read, stop_reading )
        new_dos = Atomic_DOS( data = data_array, spin_polarised = True, maximum_l_quantum_number = 2 )
        if align_fermi_energy:
            new_dos.shift_energies( -self.fermi_energy )
        return new_dos

    def read_lines_to_numpy_array( self, start, end ):
        data = []
        line_numbers = range( start, end )
        with open( self.filename, 'r' ) as file_in:
            for i, line in enumerate( file_in ):
                if i in line_numbers:
                    data.append( [ float( s ) for s in line.split() ] )
        return np.array( data )

class DOS:

    def __init__( self, data, spin_polarised ):
        assert type( data ) is np.ndarray
        self.energies = data[:,0] 
        self.densities = data[:,1:]
        self.spin_polarised = spin_polarised

    def data( self ):
        return np.concatenate( ( self.energies_as_2D_array(), self.densities ), axis = 1 )

    def data_using_columns( self, columns ):
        return np.concatenate( ( self.energies_as_2D_array(), self.densities[ :, columns ] ), axis = 1 )

    def write( self, filename = None, fmt = '%.4e', invert_down_spin = True ):
        if invert_down_spin:
            output_data = self.data()
            output_data[ :, 2::2 ] *= -1.0
        else:
            output_data = self.data()
        if filename == None:
            np.savetxt( sys.stdout.buffer, output_data, fmt = fmt )
        else:
            np.savetxt( filename, output_data, fmt = fmt )

    def energies_as_2D_array( self ):
        return self.energies.reshape( -1, 1 )

    def __add__( self, other ):
        """Add two densities of states.
           The energy ranges for each DOS must be equal,
           and the DOS data arrays must be the same size."""
        assert ( self.energies == other.energies ).all(), "DOS energies are not equal"
        new_dos = copy.deepcopy( self )
        new_dos.densities = np.add( self.densities, other.densities )
        return new_dos

    def __radd__( self, other ):
        if other == 0:
            return self
        else:
            return self + other

    def __sub__( self, other ):
        """Subtract one densities of states from another.
           The energy ranges for each DOS must be equal,
           and the DOS data arrays must be the same size."""
        assert ( self.energies == other.energies ).all(), "DOS energies are not equal"
        new_dos = copy.deepcopy( self )
        new_dos.densities = np.subtract( self.densities, other.densities )
        return new_dos

    def __str__( self ):
        return str( self.data() )

    def shift_energies( self, shift ):
        self.energies += shift
        return self

    def shift_densities( self, shift ):
        self.densities += shift
        return self

    def shift_densities_at_set( self, shift, set ):
        self.densities[ :, set ] += shift
        return self

    def scale_densities( self, scale ):
        self.densities *= scale
        return self

    @property
    def up( self ):
        assert self.spin_polarised == True
        new_dos = copy.deepcopy( self )
        new_dos.spin_polarised = False
        new_dos.densities = self.densities[ :, 0::2 ]
        return new_dos

    @property
    def down( self ):
        assert self.spin_polarised == True
        new_dos = copy.deepcopy( self )
        new_dos.spin_polarised = False
        new_dos.densities = self.densities[ :, 1::2 ]
        return new_dos

    @property
    def sum( self, columns = None ):
        if columns == None:
            columns = list( range( 0, self.densities.shape[1] ) )
        return np.sum( self.densities[ :, columns ], axis = 1 ).reshape( -1, 1 )

class Total_DOS( DOS ):

    def __init__( self, data, spin_polarised ):
        super( Total_DOS, self ).__init__( data, spin_polarised )

class Atomic_DOS( DOS ):

    def __init__( self, data, spin_polarised, maximum_l_quantum_number ):
        self.maximum_l_quantum_number = maximum_l_quantum_number
        super( Atomic_DOS, self ).__init__( data, spin_polarised )

    @property
    def s( self ):
        return self.specific_angular_momentum( 0 )

    @property
    def p( self ):
        return self.specific_angular_momentum( 1 )

    @property
    def d( self ):
        return self.specific_angular_momentum( 2 )

    @property
    def f( self ):
        return( self.specific_angular_momentum( 3 ) )

    def specific_angular_momentum( self, l ):
        if self.spin_polarised:
            columns = { 0 : list( range(0, 2) ),
                        1 : list( range(2, 8) ),
                        2 : list( range(8, 18) ),
                        3 : list( range(18, 32) ) }
        else:
            columns = { 0 : list( range(0, 1) ),
                        1 : list( range(1, 4) ),
                        2 : list( range(4, 9) ),
                        3 : list( range(9, 16) ) }
        new_dos = copy.deepcopy( self )
        new_dos.densities = self.densities[ :, columns[l] ]
        return( new_dos )

    def project_am( self ):
        new_dos = Summed_Atomic_DOS( self.data(), 
                                     spin_polarised = self.spin_polarised, 
                                     maximum_l_quantum_number = self.maximum_l_quantum_number )
        if self.spin_polarised:
            if self.maximum_l_quantum_number == 3:
                new_densities = ( self.s, 
                                  self.p.up.sum, 
                                  self.p.down.sum, 
                                  self.d.up.sum,
                                  self.d.down.sum, 
                                  self.f.up.sum,
                                  self.f.down.sum )
            else:
                new_densities = ( self.s.up.sum,
                                  self.s.down.sum,
                                  self.p.up.sum, 
                                  self.p.down.sum, 
                                  self.d.up.sum,
                                  self.d.down.sum )
        else:
            if self.maximum_l_quantum_number == 3:
                new_densities = ( self.s. 
                                  self.p.sum, 
                                  self.d.sum, 
                                  self.f.sum )
            else:
                new_densities = ( self.s, 
                                  self.p.sum, 
                                  self.d.sum )
        new_dos.densities = np.concatenate( new_densities, axis = 1 )
        return new_dos

class Summed_Atomic_DOS( Atomic_DOS ):

    def __init__( self, data, spin_polarised, maximum_l_quantum_number ):
        # self.maximum_l_quantum_number = maximum_l_quantum_number
        super( Summed_Atomic_DOS, self ).__init__( data, spin_polarised, maximum_l_quantum_number )

    def specific_angular_momentum( self, l ):
        if self.spin_polarised:
            columns = { 0 : list( range(0, 2) ),
                        1 : list( range(2, 4) ),
                        2 : list( range(4, 6) ),
                        3 : list( range(6, 8) ) }
        else:
            columns = { 0 : list( range(0, 1) ),
                        1 : list( range(1, 2) ),
                        2 : list( range(2, 3) ),
                        3 : list( range(3, 4) ) }
        new_dos = copy.deepcopy( self )
        new_dos.densities = self.densities[ :, columns[l] ]
        return( new_dos )
