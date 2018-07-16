#! /usr/bin/env python3

# Adapted from http://kitchingroup.cheme.cmu.edu/blog/2013/02/18/Nonlinear-curve-fitting/

import glob
import numpy as np
import pandas as pd
from scipy.optimize import leastsq
import argparse
import warnings

warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.outputs import UnconvergedVASPWarning

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from vasppy.poscar import Poscar
from vasppy.summary import find_vasp_calculations
from vasppy.utils import match_filename

def parse_args():
    parser = argparse.ArgumentParser(description='Perform a Murnaghan equation of state fit across VASP subdirectories')
    parser.add_argument( '-p', '--plot', action='store_true', help='generate murn.pdf plot of fit' )
    parser.add_argument( '-v', '--verbose', action='store_true', help='verbose output' )
    args = parser.parse_args()
    return args

def read_vasprun( filename ):
    return Vasprun( filename, parse_potcar_file=False, parse_dos=False, parse_eigen=False )

def read_data( verbose=True ):
    dir_list = find_vasp_calculations()
    if not dir_list:
        raise ValueError( 'Did not find any subdirectories containing vasprun.xml or vasprun.xml.gz files' )
    data = []
    for d in dir_list:
        converged = True
        try:
            with warnings.catch_warnings(record=True) as w:
                vasprun = read_vasprun( match_filename( d + 'vasprun.xml' ) )
                for warning in w:
                    if isinstance( warning.message, UnconvergedVASPWarning ):
                        converged = False
                    else:
                        print( warning.message )
        except:
            continue 
        poscar = Poscar.from_file( d + 'POSCAR' )
        data.append( [ poscar.scaling, 
                       vasprun.final_structure.volume, 
                       vasprun.final_energy,
                       converged ] )
    column_titles = [ 'scaling', 'volume', 'energy', 'converged' ]
    df = pd.DataFrame( data, columns=column_titles ).sort_values( by='scaling' )
    df = df.reset_index( drop=True )
    df['scaling_factor'] = df.volume / df.scaling**3
    scaling_factor_round = 4
    if verbose:
        print( df.to_string(index=False) )
    if len( set( df.scaling_factor.round( scaling_factor_round ) ) ) != 1:
        raise ValueError( "POSCAR scaling factors and volumes are inconsistent" )
    return df

def murnaghan( vol, e0, b0, bp, v0 ):
    """
    Calculate the energy as a function of volume, using the Murnaghan equation of state
    [Murnaghan, Proc. Nat. Acad. Sci. 30, 244 (1944)]
    https://en.wikipedia.org/wiki/Murnaghan_equation_of_state
    cf. Fu and Ho, Phys. Rev. B 28, 5480 (1983).

    Args:
        vol (float): this volume.
        e0 (float):  energy at the minimum-energy volume, E0.
        b0 (float):  bulk modulus at the minimum-energy volume, B0.
        bp (float):  pressure-derivative of the bulk modulus at the minimum-energy volume, B0'.
        v0 (float):  volume at the minimum-energy volume, V0.
        
    Returns:
        (float): The energy at this volume. 
    """
    energy = e0 + b0 * vol / bp * (((v0 / vol)**bp) / (bp - 1) + 1) - v0 * b0 / (bp - 1.0)
    return energy

def objective( pars, x, y ):
    err =  y - murnaghan( x, *pars )
    return err

def lstsq_fit( volumes, energies ):
    e_min = energies.min()
    v_min = volumes[ np.argwhere( energies == e_min )[0][0] ]
    x0 = [ e_min, 2.0, 10.0, v_min ] #initial guess of parameters
    plsq = leastsq( objective, x0, args=( volumes, energies ) )
    return plsq

def make_plot( df, fit_params ):
    v_min = df.volume.min()*0.99
    v_max = df.volume.max()*1.01
    v_fitting = np.linspace( v_min, v_max, num=50 )
    e_fitting = murnaghan( v_fitting, *fit_params )
    plt.figure( figsize=(8.0,6.0) )
    # plot converged data points
    loc = df.converged
    plt.plot( df[loc].volume, df[loc].energy, 'o' )
    # plot unconverged data points
    loc = [ not b for b in df.converged ]
    plt.plot( df[loc].volume, df[loc].energy, 'o', c='grey' )
    # plot fitted equation of state curve
    plt.plot( v_fitting, e_fitting, '--' )
    plt.xlabel( r'volume [$\mathrm{\AA}^3$]' )
    plt.ylabel( r'energy [eV]' )
    plt.tight_layout()
    plt.savefig( 'murn.pdf' )

def fit( verbose=False, plot=False ):
    df = read_data( verbose=verbose )
    e0, b0, bp, v0 = lstsq_fit( np.array( df.volume ), np.array( df.energy  ) )[0]
    if plot:
        make_plot( df, ( e0, b0, bp, v0 ) )
    print( "E0: {:.4f}".format( e0 ) )
    print( "V0: {:.4f}".format( v0 ) )
    print( "opt. scaling: {:.5f}".format( ( v0 / df.scaling_factor.mean() )**(1/3) ) )
    
def main():
    args = parse_args()
    fit( verbose=args.verbose, plot=args.plot )

if __name__ == '__main__':
    main()
