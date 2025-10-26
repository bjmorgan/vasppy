checkforce
==========

The ``checkforce`` script checks the convergence of geometry optimisations in VASP.

Synopsis
--------

.. code-block:: bash

    checkforce [-h] [-o OUTCAR] [-c CONVERGENCE] [-v | -a]

Description
-----------

Analyses forces from a VASP OUTCAR file and reports convergence statistics. 
By default, reports statistics for the last ionic step only.

Options
-------

.. option:: -h, --help

   Show help message and exit.

.. option:: -o OUTCAR, --outcar OUTCAR

   The filepath of the OUTCAR file to be parsed. Default is "OUTCAR".

.. option:: -c CONVERGENCE, --convergence CONVERGENCE

   Set force convergence threshold. Default is to read the convergence from the OUTCAR file (EDIFFG).

.. option:: -v, --verbose

   Verbose output. Show convergence status for all atoms in the last ionic step.

.. option:: -a, --all

   Print summary data for every ionic step.

.. note::
   The ``-v`` and ``-a`` options are mutually exclusive.

Examples
--------

Check convergence using default OUTCAR::

    checkforce

Check with a specific OUTCAR file::

    checkforce -o path/to/OUTCAR

Set a custom convergence threshold::

    checkforce -c 0.01

View summary for all ionic steps::

    checkforce -a

Output shows:
    - **Max Force**: Maximum force magnitude across all atoms
    - **Non-opt**: Number of atoms not yet converged
    - **Mean Excess**: Mean excess force above convergence threshold

Verbose output for last step::

    checkforce -v

Output Format
-------------

Default output (last step only)::

    remainder:  0.001234
    maximum:    0.023456
    non-opt:    3 / 24

All steps output (``-a`` flag)::

    Max Force      Non-opt        Mean Excess    
    0.156789       12             0.012345       
    0.089012       8              0.006789       
    0.023456       3              0.001234       
    0.009876       0              0.000000
