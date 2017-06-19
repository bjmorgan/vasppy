import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = 'README.md'
try:
    import pypandoc
    long_description = pypandoc.convert( readme, 'rst')
except ImportError:
    long_description = open( readme ).read()

from vasppy import __version__ as VERSION

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
    'packages': ['vasppy'], 
    'license': 'MIT',
    'install_requires': [ 'numpy',
                          'pymatgen',
                          'PyYAML', 
                          'coverage==4.3.4',
                          'codeclimate-test-reporter' ]
}

setup(**config)
