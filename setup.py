"""
vasppy: Python utilities for working with VASP inputs and outputs
"""

from setuptools import setup, find_packages
from vasppy import __version__ as VERSION

readme = 'README.md'
long_description = open( readme ).read()

scripts = [ 'check_species',
            'murnfit', 
            'vasp_summary', 
            'poscar_to_cif', 
            'potcar_spec',
            'effective_mass',
            'fat_bands',
            'pimaim_to_poscar',
            'pimaim_to_xtl',
            'poscar_sort',
            'poscar_to_pimaim',
            'poscar_to_xtl',
            'proc_poscar',
            'rotate_poscar',
            'spacegroup',
            'vasp_grid',
            'xdatcar_to_disp',
            'xdatcar_to_poscart',
            'xdatcar_to_rdf' ]

setup(
    name='vasppy',
    version=VERSION,
    description='Python utilities for working with VASP inputs and outputs',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Benjamin J. Morgan',
    author_email='bjm42@bath.ac.uk',
    url='https://github.com/bjmorgan/vasppy', 
    download_url='https://github.com/bjmorgan/vasppy/archive/{}.tar.gz'.format( VERSION ),
    keywords=['vasp'], # keywords
    packages=find_packages( exclude=['docs', 'tests*'] ),
    package_data={ 'vasppy': ['data/*.yaml'] },
    entry_points={ 'console_scripts': [
                       '{} = vasppy.scripts.{}:main'.format( s, s ) for s in scripts ] },
    license='MIT',
    install_requires=[ 'monty',
                       'numpy',
                       'pandas',
                       'pymatgen',
                       'PyYAML', 
                       'coverage==4.3.4',
                       'codeclimate-test-reporter',
                       'fortranformat' ],
    python_requires='>=3.5'
    )
