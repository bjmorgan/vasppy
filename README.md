# vasppy - a Python suite for manipulating VASP files

[![DOI](https://zenodo.org/badge/17946870.svg)](https://zenodo.org/badge/latestdoi/17946870)
[![Build Status](https://travis-ci.org/bjmorgan/vasppy.svg?branch=master)](https://travis-ci.org/bjmorgan/vasppy)
[![Test Coverage](https://codeclimate.com/github/bjmorgan/vasppy/badges/coverage.svg)](https://codeclimate.com/github/bjmorgan/vasppy/coverage)

`vasppy` is a suite of Python tools and scripts written in Python for manipulating and processing [VASP](https://www.vasp.at/) input and output files.

## Installation of executable scripts

The suite of executable scripts can be installed by running `install.py`. This creates symbolic links with a default destination of `$HOME/bin/`. For more options, such as installing selected scripts, or uninstalling (removing symbolic links) ) use `./install.py --help`

## Tests

Automated testing of the latest build happens [here](https://travis-ci.org/bjmorgan/vasppy).

Manual tests can be run using
```
python3 -m unittest discover
```

## Contributors

Benjamin J. Morgan
Lucy Whalley
