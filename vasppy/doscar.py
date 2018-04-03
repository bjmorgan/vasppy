import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd

tableau_grey = '#bab0ac'

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
        if self.lmax == 3:
            raise NotImplementedError( 'f-orbital projects DOSCARs are not yet supported' )
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

    def plot_pdos(self, ax=None, to_plot=None, colors=None, 
                  plot_total_dos=True, xrange=None, ymax=None, 
                  scaling=None, split=False, title=None, title_loc='center',
                  labels=True, title_fontsize=16):
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
        ax.legend(bbox_to_anchor=(1.01, 1.04), loc='upper left')
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
