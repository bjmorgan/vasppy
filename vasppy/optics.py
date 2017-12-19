"""
functions for working with optical properties from vasprun.xml
"""

from math import pi, sqrt
import numpy as np
from scipy.constants import physical_constants, speed_of_light

eV_to_recip_cm = 1.0/(physical_constants['Planck constant in eV s'][0]*speed_of_light*1e2)


def diagonalise_matrix( m ):
    eigvals, eigvecs = np.linalg.eig(m)
    print( eigvals )
    print( eigvecs )
    print( np.dot(eigvecs.T, m) )
    print( np.dot( np.dot( eigvecs.T, m ), eigvals ) )
    
def matrix_eigvals( m ):
    eigvals, eigvecs = np.linalg.eig(m)
    return eigvals
       
def to_matrix( xx, yy, zz, xy, yz, xz ):
    matrix = np.array( [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]] )
    return matrix

def plot_dielectric_functions( dielectric, ax=None ):
    real_dielectric = np.array( [ matrix_eigvals( to_matrix( *e ) ) for e in dielectric[1] ] )
    imag_dielectric = np.array( [ matrix_eigvals( to_matrix( *e ) ) for e in dielectric[2] ] )
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6.0,3.0))
    else:
        fig = None
    ax.plot( dielectric[0], np.mean( real_dielectric, axis=1 ), '-', zorder=2 ) # better to pass in v.dielectric
    ax.plot( dielectric[0], np.mean( imag_dielectric, axis=1 ), '-', zorder=2 )
    ax.set_xlim([0,8])
    ax.set_ylim([0,5])
    return fig

def absorption_coefficient( dielectric ):
    energies_in_eV = np.array( dielectric[0] )
    real_dielectric = np.array( [ matrix_eigvals( to_matrix( *e ) ) for e in dielectric[1] ] )
    imag_dielectric = np.array( [ matrix_eigvals( to_matrix( *e ) ) for e in dielectric[2] ] )
    epsilon_1 = np.mean( real_dielectric, axis=1 )
    epsilon_2 = np.mean( imag_dielectric, axis=1 )
    return ( 2.0 * np.sqrt(2.0)*pi*eV_to_recip_cm*energies_in_eV
                 * np.sqrt( -epsilon_1 + np.sqrt( epsilon_1**2 + epsilon_2**2 ) ) )
