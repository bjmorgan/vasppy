#! /usr/bin/env python3

# Adapted from http://kitchingroup.cheme.cmu.edu/blog/2013/02/18/Nonlinear-curve-fitting/

import glob
import numpy as np
import pandas as pd
from scipy.optimize import leastsq
from pymatgen.io.vasp import Vasprun
from vasppy import Poscar
from vasppy.summary import find_vasp_calculations
import argparse

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Perform a Murnaghan equation of state fit across VASP subdirectories')
    parser.add_argument( '-p', '--plot', action='store_true', help='generate murn.pdf plot of fit' )
    args = parser.parse_args()
    return args

def read_data( verbose=True ):
    dir_list = find_vasp_calculations()
    data = []
    for d in dir_list:
        try:
            vasprun = Vasprun( d + 'vasprun.xml', parse_potcar_file=False )
        except:
            continue
        poscar = Poscar.from_file( d + 'POSCAR' )
        data.append( np.array( [ poscar.scaling, vasprun.final_structure.volume, vasprun.final_energy ] ) )
    df = pd.DataFrame( data, columns=[ 'scaling', 'volume', 'energy' ] ).sort_values( by='scaling' )
    df = df.reset_index( drop=True )
    df['scaling_factor'] = df.volume / df.scaling**3
    scaling_factor_round = 5
    if verbose:
        print( df )
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

def fit( volumes, energies ):
    e_min = energies.min()
    v_min = volumes[ np.argwhere( energies == e_min )[0][0] ]
    x0 = [ e_min, 2.0, 10.0, v_min ] #initial guess of parameters
    plsq = leastsq( objective, x0, args=( volumes, energies ) )
    return plsq

def make_plot( volumes, energies, fit_params ):
    v_min = volumes.min()*0.99
    v_max = volumes.max()*1.01
    v_fitting = np.linspace( v_min, v_max, num=50 )
    e_fitting = murnaghan( v_fitting, *fit_params )
    plt.figure( figsize=(8.0,6.0) )
    plt.plot( volumes, energies, 'o' )
    plt.plot( v_fitting, e_fitting, '--' )
    plt.xlabel( 'volume [A^3]' )
    plt.ylabel( 'energy [eV]' )
    plt.tight_layout()
    plt.savefig( 'murn.pdf' )

if __name__ == '__main__':
    args = parse_args()
    df = read_data()
    e0, b0, bp, v0 = fit( np.array( df.volume ), np.array( df.energy  ) )[0]
    if args.plot:
        make_plot( df.volume, df.energy, ( e0, b0, bp, v0 ) )
    print( "E0: {:.4f}".format( e0 ) )
    print( "V0: {:.4f}".format( v0 ) )
    print( "opt. scaling: {:.5f}".format( ( v0 / df.scaling_factor.mean() )**(1/3) ) )
