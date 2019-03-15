import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd

tableau_grey = '#bab0ac'

def pdos_column_names( lmax, ispin ):
    if lmax == 2:
        names = [ 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2' ]
    elif lmax == 3:
        names = [ 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2',
                  'f_y(3x2-y2)', 'f_xyz', 'f_yz2', 'f_z3', 'f_xz2', 'f_z(x2-y2)', 'f_x(x2-3y2)' ]
    else:
        raise ValueError( 'lmax value not supported' )
    if ispin == 2:
        all_names = []
        for n in names:
            all_names.extend( [ '{}_up'.format(n), '{}_down'.format(n) ] )
    else:
        all_names = names
    all_names.insert( 0, 'energy' )
    return all_names

class Doscar:
    '''
    Contains all the data in a VASP DOSCAR file, and methods for manipulating this.
    '''

    number_of_header_lines = 6

    def __init__( self, filename, ispin=2, lmax=2, lorbit=11, spin_orbit_coupling=False, read_pdos=True, species=None ):
        '''
        Create a Doscar object from a VASP DOSCAR file.

        Args:
            filename (str): Filename of the VASP DOSCAR file to read.
            ispin (optional:int): ISPIN flag. 
                Set to 1 for non-spin-polarised or 2 for spin-polarised calculations.
                Default = 2.
            lmax (optional:int): Maximum l angular momentum. (d=2, f=3). Default = 2.
            lorbit (optional:int): The VASP LORBIT flag. (Default=11).
            spin_orbit_coupling (optional:bool): Spin-orbit coupling (Default=False).
            read_pdos (optional:bool): Set to True to read the atom-projected density of states (Default=True).
            species (optional:list(str)): List of atomic species strings, e.g. [ 'Fe', 'Fe', 'O', 'O', 'O' ].
                Default=None.
        '''
        self.filename = filename
        self.ispin = ispin
        self.lmax = lmax
        self.spin_orbit_coupling = spin_orbit_coupling
        if self.spin_orbit_coupling:
            raise NotImplementedError( 'Spin-orbit coupling is not yet implemented' )
        self.lorbit = lorbit
        self.pdos = None
        self.species = species
        self.read_header()
        self.read_total_dos()
        if read_pdos:
            try:
                self.read_projected_dos()
            except:
                raise
        # if species is set, should check that this is consistent with the number of entries in the
        # projected_dos dataset
        
    @property
    def number_of_channels( self ):
        if self.lorbit == 11:
            return { 2: 9, 3: 16 }[ self.lmax ]
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
                          names=pdos_column_names( lmax=self.lmax, ispin=self.ispin ),
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

        Args:
            atoms (int or list(int)): Atom numbers to include in the selection. Atom numbers count from 1. 
                                   Default is to select all atoms.
            spin (str): Select up or down, or both spin channels to include in the selection.
                        Accepted options are 'up', 'down', and 'both'. Default is to select both spins.
            l (str): Select one angular momentum to include in the selectrion.
                     Accepted options are 's', 'p', 'd', and 'f'. Default is to include all l-values.
                     Setting `l` and not setting `m` will return all projections for that angular momentum value.
            m (list(str)): Select one or more m-values. Requires `l` to be set. 
                           The accepted values depend on the value of `l`:
                           `l='s'`: Only one projection. Not set.
                           `l='p'`: One or more of [ 'x', 'y', 'z' ]
                           `l='d'`: One or more of [ 'xy', 'yz', 'z2-r2', 'xz', 'x2-y2' ]
                           `l='f'`: One or more of [ 'y(3x2-y2)', 'xyz', 'yz2', 'z3', 'xz2', 'z(x2-y2)', 'x(x2-3y2)' ]

        Returns:
            np.array: A 4-dimensional numpy array containing the selected pdos values. 
            The array dimensions are [ atom_no, energy_value, lm-projection, spin ]

        """
        valid_m_values = { 's': [],
                           'p': [ 'x', 'y', 'z' ],
                           'd': [ 'xy', 'yz', 'z2-r2', 'xz', 'x2-y2' ],
                           'f': [ 'y(3x2-y2)', 'xyz', 'yz2', 'z3', 'xz2', 'z(x2-y2)', 'x(x2-3y2)' ] }
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
            raise ValueError( "valid spin values are 'up', 'down', and 'both'. The default is 'both'" )
        to_return = to_return[ :, :, :, spin_idx ]
        if not l:
            channel_idx = list(range( self.number_of_channels ))
        elif l == 's':
            channel_idx = [ 0 ]
        elif l == 'p':
            if not m:
                channel_idx = [ 1, 2, 3 ]
            else: # TODO this looks like it should be i+1
                channel_idx = [ i+1 for i, v in enumerate( valid_m_values['p'] ) if v in m ]
        elif l == 'd':
            if not m:
                channel_idx = [ 4, 5, 6, 7, 8 ]
            else: # TODO this looks like it should be i+4
                channel_idx = [ i+4 for i, v in enumerate( valid_m_values['d'] ) if v in m ]
        elif l == 'f':
            if not m:
                channel_idx = [ 9, 10, 11, 12, 13, 14, 15 ]
            else: # TODO this looks like it should be i+9
                channel_idx = [ i+9 for i, v in enumerate( valid_m_values['f'] ) if v in m ]
        else:
            raise ValueError
        return to_return[ :, :, channel_idx, : ]
    
    def pdos_sum( self, atoms=None, spin=None, l=None, m=None ):
        return np.sum( self.pdos_select( atoms=atoms, spin=spin, l=l, m=m ), axis=(0,2,3) )

    def plot_pdos(self, ax=None, to_plot=None, colors=None, 
                  plot_total_dos=True, xrange=None, ymax=None, 
                  scaling=None, split=False, title=None, title_loc='center',
                  labels=True, title_fontsize=16, legend_pos='outside'):
        if not ax:
            fig, ax = plt.subplots(1, 1, figsize=(8.0,3.0))
        else:
            fig = None
        if not colors:
            colors = mcd.TABLEAU_COLORS        
        color_iterator = (c for c in colors)
        
        if not scaling:
            scaling = []
            
        if xrange:
            e_range = (self.energy >= xrange[0]) & (self.energy <= xrange[1])
        else:
            e_range = np.ma.make_mask( self.energy )
            
        auto_ymax = 0.0
            
        if not to_plot:
            to_plot = {}
            for s in set( self.species ):
                to_plot[s] = ['s', 'p', 'd']
                if self.lmax == 3:
                    to_plot[s].append('f')
                    
        for species in to_plot.keys():
            index = [i for i, s in enumerate(self.species) if s is species]
            for state in to_plot[species]:
                assert state in ['s', 'p', 'd', 'f']
                color = next( color_iterator )
                label = '{} {}'.format(species, state)
                up_dos = self.pdos_sum(atoms=index, l=state, spin='up')[e_range]
                down_dos = self.pdos_sum(atoms=index, l=state, spin='down')[e_range]
                if species in scaling:
                    if state in scaling[species]:
                        up_dos *= scaling[species][state]
                        down_dos *= scaling[species][state]
                        label = r'{} {} $\times${}'.format( species, state, scaling[species][state] )
                auto_ymax = max( [ auto_ymax, up_dos.max(), down_dos.max() ] )
                ax.plot(self.energy[e_range], up_dos, label=label, c=color)
                ax.plot(self.energy[e_range], down_dos * -1.0,  c=color)
        if plot_total_dos:
            ax.fill_between(self.energy[e_range], self.tdos.up.values[e_range], 
                            self.tdos.down.values[e_range] * -1.0, facecolor=tableau_grey, alpha=0.2)
            auto_ymax = max( [ auto_ymax, self.tdos.up.values[e_range].max(), self.tdos.down.values[e_range].max() ] )
    
        if xrange:
            ax.set_xlim( xrange[0], xrange[1] )
            
        if not ymax:
            ymax = 1.1 * auto_ymax
        ax.set_ylim(-ymax*1.1,ymax*1.1)
        if legend_pos == 'outside':
            ax.legend(bbox_to_anchor=(1.01, 1.04), loc='upper left')
        else:
            ax.legend( loc=legend_pos )
        if labels:
            ax.set_xlabel( 'Energy [eV]')
        ax.axhline(y=0, c='lightgrey')
        ax.axes.grid( False, axis='y' )
    
        ax.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left=False,      # ticks along the left edge are off
            right=False,         # ticks along the right edge are off
            labelleft=False) # labels along the left edge are off
        
        if title:
            ax.set_title( title, loc=title_loc, fontdict={'fontsize': title_fontsize} )
    
        return fig
