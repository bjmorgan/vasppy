import numpy as np
import pandas as pd

class Doscar:
    '''
    Contains all the data in a VASP DOSCAR file, and methods for manipulating this.
    '''

    number_of_header_lines = 6

    def __init__( self, filename, ispin=2, lmax=2, lorbit=11, spin_orbit_coupling=False, read_pdos=True ):
        self.filename = filename
        self.ispin = ispin
        self.lmax = lmax
        self.spin_orbit_coupling = spin_orbit_coupling
        self.lorbit = lorbit
        self.pdos = None
        self.read_header()
        self.read_total_dos()
        if read_pdos:
            try:
                self.read_projected_dos()
            except:
                raise
        
    @property
    def number_of_channels( self ):
        if self.lorbit == 11:
            if self.lmax == 2:
                return 9
        raise notImplementedError

    def read_header( self ):
        self.header = []
        with open( self.filename, 'r' ) as file_in:
            for i in range( Doscar.number_of_header_lines ):
                self.header.append( file_in.readline() )
        self.process_header()

    def process_header( self ):
        self.number_of_atoms = int( self.header[0].split()[0] )
        self.number_of_data_points = int( self.header[5].split()[2] )
        self.efermi = float( self.header[5].split()[3] )
        
    def read_total_dos( self ): # assumes spin_polarised
        start_to_read = Doscar.number_of_header_lines
        df = pd.read_csv( self.filename, 
                          skiprows=start_to_read, 
                          nrows=self.number_of_data_points,
                          delim_whitespace=True, 
                          names=[ 'energy', 'up', 'down', 'int_up', 'int_down' ],
                          index_col=False )
        self.energy = df.energy.values
        df.drop( 'energy', axis=1 )
        self.tdos = df
        
    def read_atomic_dos_as_df( self, atom_number ): # currently assume spin-polarised, no-SO-coupling, no f-states
        assert atom_number > 0 & atom_number <= self.number_of_atoms
        start_to_read = Doscar.number_of_header_lines + atom_number * ( self.number_of_data_points + 1 )
        df = pd.read_csv( self.filename,
                          skiprows=start_to_read,
                          nrows=self.number_of_data_points,
                          delim_whitespace=True,
                          names=[ 'energy', 's_up', 's_down', 
                                  'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                                  'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z2-r2_up', 'd_z2-r2_down',
                                  'd_xz_up', 'd_xz_down', 'd_x2-y2_up', 'd_x2-y2_down' ],
                          index_col=False )
        return df.drop('energy', axis=1)
    
    def read_projected_dos( self ):
        """
        Read the projected density of states data into """
        pdos_list = []
        for i in range( self.number_of_atoms ):
            df = self.read_atomic_dos_as_df( i+1 )
            pdos_list.append( df )
        self.pdos = np.vstack( [ np.array( df ) for df in pdos_list ] ).reshape( 
            self.number_of_atoms, self.number_of_data_points, self.number_of_channels, self.ispin )
        
    def pdos_select( self, atoms=None, spin=None, l=None, m=None ):
        """
        Returns a subset of the projected density of states array.
        """
        valid_m_values = { 's': [],
                           'p': [ 'x', 'y', 'z' ],
                           'd': [ 'xy', 'yz', 'z2-r2', 'xz', 'x2-y2' ] }
        if not atoms:
            atom_idx = list(range( self.number_of_atoms ))
        else:
            atom_idx = atoms
        to_return = self.pdos[ atom_idx, :, :, : ]
        if not spin:
            spin_idx = list(range( self.ispin ))
        elif spin is 'up':
            spin_idx = [0]
        elif spin is 'down':
            spin_idx = [1]
        elif spin is 'both':
            spin_idx = [0,1]
        else:
            raise ArgumentError
        to_return = to_return[ :, :, :, spin_idx ]
        if not l:
            channel_idx = list(range( self.number_of_channels ))
        elif l == 's':
            channel_idx = [ 0 ]
        elif l == 'p':
            if not m:
                channel_idx = [ 1, 2, 3 ]
            else:
                channel_idx = [ i for i, v in enumerate( valid_m_values['p'] ) if v in m ]
        elif l == 'd':
            if not m:
                channel_idx = [ 4, 5, 6, 7, 8 ]
            else:
                channel_idx = [ i for i, v in enumerate( valid_m_values['d'] ) if v in m ]
        return to_return[ :, :, channel_idx, : ]
    
    def pdos_sum( self, atoms=None, spin=None, l=None, m=None ):
        return np.sum( self.pdos_select( atoms=atoms, spin=spin, l=l, m=m ), axis=(0,2,3) )

