from setuptools import setup, find_packages

readme = 'README.md'
try:
    import pypandoc
    long_description = pypandoc.convert( readme, 'rst')
except ImportError:
    long_description = open( readme ).read()

from vasppy import __version__ as VERSION

scripts = [ 'murnfit', 
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
            'xdatcar_to_rdf' ]

config = {
    'name': 'vasppy',
    'description': 'A library for manipulating VASP input / output files',
    'long_description': long_description,
    'author': 'Benjamin J. Morgan',
    'author_email':'bjm42@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/vasppy', 
    'download_url': "https://github.com/bjmorgan/vasppy/archive/%s.tar.gz" % (VERSION),
    'version': VERSION,
    'keywords': ['vasp'], # keywords
    'packages': find_packages( exclude=['docs', 'tests*']),
    'package_data': {'vasppy': ['data/*.yaml']},
    'entry_points': { 'console_scripts': [
                          '{} = vasppy.scripts.{}:main'.format( s, s ) for s in scripts ] },
    'license': 'MIT',
    'install_requires': [ 'numpy',
                          'pandas',
                          'pymatgen',
                          'PyYAML', 
                          'coverage==4.3.4',
                          'codeclimate-test-reporter' ]
}

setup(**config)
