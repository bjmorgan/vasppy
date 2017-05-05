from setuptools import setup

setup(
    name = 'vasppy',
    packages = ['vasppy'], # this must be the same as the name above
    version = '0.1.0',
    description = 'A library for manipulating VASP input / output files',
    author = 'Benjamin J. Morgan',
    author_email = 'bjm42@bath.ac.uk',
    url = 'https://github.com/bjmorgan/vasppy',   # use the URL to the github repo
    download_url = 'https://github.com/bjmorgan/vasppy/tarball/0.1.0',
    keywords = ['vasp'], # keywords
    install_requires = [
        'numpy',
        'pymatgen',
        'PyYAML', 
        'coverage',
        'codeclimate-test-reporter'
    ],
    classifiers = []
)
