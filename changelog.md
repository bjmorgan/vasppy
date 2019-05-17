2019-17-05:
- Substantial rewrite of the `Procar` class:
    - The preferred ways to read `PROCAR` files are now `Procar.from_file()` and `Procar.from_files()`. `Procar.read_from_file()` is deprecated as a public method and will be removed in a future version. 
    - `Procar.from_files()` takes a list of file path strings, and collates each PROCAR in sequence. This is useful where e.g. a band structure calculation has been split over multiple VASP calculations.
    - `Procar.from_file()` now has an optional `select_zero_weighted_kpoints` argument. Setting this to `True` will return a `Procar` object that only has those k-points with zero weights. This is useful for processing data from hybrid DFT calculations, where some k-points might have zero weightings.
    - `Procar.from_files()` inherits all keyword arguments  from `Procar.from_file()`.
    - `Procar` objects can be concatenated, using normal addition operators.
    - Specific k-points can be selected, as a new `Procar` object, using the `select_k_points()` method.

2018-07-13:
- `murnfit` script now reads from unzipped `vasprun.xml` and gzipped `vasprun.xml.gz` files.
- `murnfit` script now takes `-v`, `--verbose` argument to produce a verbose output.
- `murnfit` script plots (using `-p` or `--plot` arguments) now use correct typography for axis labels (using matplotlib LaTeX support)
