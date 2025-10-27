convergence_testing
===================

The ``convergence_testing`` script generates a series of VASP input files for systematic convergence testing of plane-wave cutoff energy (ENCUT) and k-point density (KSPACING).

Synopsis
--------

.. code-block:: bash

	convergence_testing [-h] -i INCAR -p POSCAR 
					   (--pseudopotentials PSEUDO [PSEUDO ...] | --potcar-file POTCAR_FILE)
					   [-e EMIN EMAX ESTEP] [-k KMIN KMAX KSTEP]
					   [-d DIRECTORY] [--base-encut BASE_ENCUT]
					   [--base-kspacing BASE_KSPACING]
					   [--job-script JOB_SCRIPT] [--dry-run]

Description
-----------

Generates a complete set of VASP input files for convergence testing. The script creates two series of calculations:

1. **ENCUT convergence**: Tests varying plane-wave cutoff energies with fixed k-point density
2. **KSPACING convergence**: Tests varying k-point densities with fixed plane-wave cutoff

For each test point, the script creates a directory containing:
   - POSCAR (crystal structure)
   - INCAR (calculation parameters)
   - POTCAR (pseudopotentials)
   - Optionally, a job submission script

This systematic approach ensures that both computational parameters are properly converged for accurate DFT calculations.

Required Arguments
------------------

.. option:: -i INCAR, --incar INCAR

   Path to template INCAR file. This file should contain all VASP parameters except ENCUT and KSPACING, which will be set automatically for each test.

.. option:: -p POSCAR, --poscar POSCAR

   Path to POSCAR file containing the crystal structure.

.. option:: --pseudopotentials PSEUDO [PSEUDO ...]

   Specify pseudopotentials to use, e.g., ``Li_sv Mn_pv O``. 
   
   Mutually exclusive with ``--potcar-file``.
   
   .. note::
	  This option requires Pymatgen to be configured for VASP POTCARs. 
	  See the `Pymatgen documentation <https://pymatgen.org/installation.html#potcar-setup>`_ for setup instructions.

.. option:: --potcar-file POTCAR_FILE

   Path to an existing POTCAR file to use for all calculations.
   
   Mutually exclusive with ``--pseudopotentials``.

Optional Arguments
------------------

.. option:: -e EMIN EMAX ESTEP, --encut EMIN EMAX ESTEP

   Set the minimum, maximum, and step size for ENCUT (in eV).
   
   Default: ``100 700 50``

.. option:: -k KMIN KMAX KSTEP, --kspacing KMIN KMAX KSTEP

   Set the minimum, maximum, and step size for KSPACING (in Å⁻¹).
   
   Default: ``0.1 0.8 0.02``

.. option:: -d DIRECTORY, --directory DIRECTORY

   Specify the output directory for generated files.
   
   Default: ``./convergence_testing``

.. option:: --base-encut BASE_ENCUT

   Fixed ENCUT value (in eV) to use for KSPACING convergence tests.
   
   Default: ``400``

.. option:: --base-kspacing BASE_KSPACING

   Fixed KSPACING value (in Å⁻¹) to use for ENCUT convergence tests.
   
   Default: ``0.3``

.. option:: --job-script JOB_SCRIPT

   Path to a job submission script to copy into each test directory. Useful for automatically setting up calculations on HPC clusters.

.. option:: --dry-run

   Show what would be created without actually creating any files or directories. Useful for checking the setup before running.

.. option:: -h, --help

   Show help message and exit.

Examples
--------

Basic usage with pseudopotential specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate convergence tests for Li₂O using automatic pseudopotential selection::

	convergence_testing -i INCAR -p POSCAR --pseudopotentials Li_sv O

Basic usage with existing POTCAR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generate convergence tests using an existing POTCAR file::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR

Custom convergence ranges
^^^^^^^^^^^^^^^^^^^^^^^^^^

Test ENCUT from 200 to 600 eV in steps of 50 eV, and KSPACING from 0.15 to 0.5 Å⁻¹ in steps of 0.05 Å⁻¹::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR \
					   -e 200 600 50 \
					   -k 0.15 0.5 0.05

Custom output directory
^^^^^^^^^^^^^^^^^^^^^^^

Generate files in a specific directory::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR \
					   -d ./my_convergence_tests

Adjust base parameters
^^^^^^^^^^^^^^^^^^^^^^^

Use ENCUT=500 eV for KSPACING tests and KSPACING=0.25 Å⁻¹ for ENCUT tests::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR \
					   --base-encut 500 \
					   --base-kspacing 0.25

Include job submission script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy a job script into each test directory for easy submission::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR \
					   --job-script submit.sh

Preview without creating files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check what will be created before running::

	convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR --dry-run

Output Structure
----------------

The script creates the following directory structure::

	convergence_testing/
	├── ENCUT/
	│   ├── 100/
	│   │   ├── INCAR    (with ENCUT=100, KSPACING=0.3)
	│   │   ├── POSCAR
	│   │   ├── POTCAR
	│   │   └── submit.sh (if --job-script provided)
	│   ├── 150/
	│   │   ├── INCAR    (with ENCUT=150, KSPACING=0.3)
	│   │   ├── POSCAR
	│   │   ├── POTCAR
	│   │   └── submit.sh
	│   └── ...
	└── KSPACING/
		├── 0.1/
		│   ├── INCAR    (with ENCUT=400, KSPACING=0.1)
		│   ├── POSCAR
		│   ├── POTCAR
		│   └── submit.sh
		├── 0.12/
		│   ├── INCAR    (with ENCUT=400, KSPACING=0.12)
		│   ├── POSCAR
		│   ├── POTCAR
		│   └── submit.sh
		└── ...

Each subdirectory contains a complete set of VASP input files ready for submission.

Workflow Integration
--------------------

Typical workflow for convergence testing:

1. **Prepare template files**::

	   # Create a template INCAR with all parameters except ENCUT and KSPACING
	   # Ensure POSCAR contains your structure
	   # Have POTCAR ready or know which pseudopotentials to use

2. **Generate test directories**::

	   convergence_testing -i INCAR -p POSCAR --potcar-file POTCAR \
						  --job-script submit.sh

3. **Submit calculations**::

	   # For ENCUT tests
	   for dir in convergence_testing/ENCUT/*/; do
		   cd "$dir"
		   sbatch submit.sh  # or your job submission command
		   cd -
	   done
	   
	   # For KSPACING tests
	   for dir in convergence_testing/KSPACING/*/; do
		   cd "$dir"
		   sbatch submit.sh
		   cd -
	   done

4. **Analyse results**:

   After calculations complete, extract energies from each directory and plot convergence curves to determine optimal parameters.

.. tip::
   Start with the ``--dry-run`` flag to verify the setup before creating files.

.. tip::
   The ENCUT and KSPACING ranges can be refined after initial tests. Run a coarse initial sweep, then focus on the convergence region with finer steps.

.. warning::
   The output directory must not already exist. The script will raise an error if the directory exists to prevent accidental overwriting.

Common Use Cases
----------------

Standard convergence test
^^^^^^^^^^^^^^^^^^^^^^^^^

For most systems, the default parameters provide a good starting point::

	convergence_testing -i INCAR -p POSCAR --pseudopotentials Li_sv O

High-accuracy requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For high-accuracy calculations (e.g., formation energies), test a wider ENCUT range::

	convergence_testing -i INCAR -p POSCAR --pseudopotentials Li_sv O \
					   -e 200 800 50 \
					   --base-encut 600

Large unit cells
^^^^^^^^^^^^^^^^

For large cells, test finer k-point grids (lower KSPACING values)::

	convergence_testing -i INCAR -p POSCAR --pseudopotentials Li_sv O \
					   -k 0.2 0.8 0.05 \
					   --base-kspacing 0.4

Notes
-----

- **ENCUT**: Represents the plane-wave cutoff energy in eV. Higher values give more accurate results but increase computational cost.

- **KSPACING**: Represents the minimum allowed spacing between k-points in Å⁻¹. Lower values mean denser k-point grids and higher accuracy but greater computational cost.

- The script uses reciprocal lattice vectors to determine appropriate KSPACING values for your specific crystal structure.

- Each convergence series tests one parameter whilst holding the other fixed at its base value, allowing independent assessment of convergence for both parameters.

See Also
--------

- :doc:`checkforce` - Check geometry optimisation convergence
- `VASP Manual: ENCUT <https://www.vasp.at/wiki/index.php/ENCUT>`_
- `VASP Manual: KSPACING <https://www.vasp.at/wiki/index.php/KSPACING>`_