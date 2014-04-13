from distutils.core import setup
    setup(
        name = 'vasppy',
        packages = ['vasppy'], # this must be the same as the name above
        version = '0.0.1',
        description = 'A library for manipulating VASP input / output files',
        author = 'Benjamin J. Morgan',
        author_email = 'bmorgan@liv.ac.uk',
        url = 'https://github.com/bjmorgan/vasppy',   # use the URL to the github repo
        download_url = 'https://github.com/bjmorgan/vasppy/tarball/0.1', # I'll explain this in a second
        keywords = ['vasp'], # arbitrary keywords
        classifiers = [],
    )
