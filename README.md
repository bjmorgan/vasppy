# vasppy - a python package for manipulating VASP files

Modules:
- install
- grid
- poscar
- poscar_to_xtl
- super
- vasp_grid

To install exectuable scripts (copied as symbolic links) use `./install.py`  
For installation options use `./install.py --help`

The following modules can be run as command line applications:

## install.py

    usage: install.py [-h] [--prefix PREFIX] [--list] [-f] [-r] [name [name ...]]
    
    Installs exectuable python scripts Default location is $HOME/bin
    
    positional arguments:
      name             name(s) of scripts to install
    
    optional arguments:
      -h, --help       show this help message and exit
      --prefix PREFIX  specify directory to install files to. (default $HOME/bin/)
      --list           list installable scripts
      -f, --force      force files to be overwritten on installation
      -r, --remove     remove files

## poscar_to_xtl.py

    usage: poscar_to_xtl.py [-h] poscar
    
    Converts a VASP POSCAR file to the .xtl file format
    
    positional arguments:
      poscar      filename of the VASP POSCAR to be processed
    
    optional arguments:
      -h, --help  show this help message and exit

## super.py

    usage: super.py [-h] [-l {1,4}] [-c] [-t {c,cartesian,d,direct}] [-g] [-s h k l]
                 poscar
    
    Manipulates VASP POSCAR files
    
    positional arguments:
      poscar                filename of the VASP POSCAR to be processed
    
    optional arguments:
      -h, --help            show this help message and exit
      -l {1,4}, --label {1,4}
                            label coordinates with atom name at position {1,4}
      -c, --coordinates-only
                            only output coordinates
      -t {c,cartesian,d,direct}, --coordinate-type {c,cartesian,d,direct}
                            specify coordinate type for output
                            {(c)artesian|(d)irect} [default = (d)irect]
      -g, --group           group atoms within supercell
      -s h k l, --s

## vasp_grid

    usage: vasp_grid [-h] gridfile
    
    z-projection of a VASP (grid format) file
    
    positional arguments:
      gridfile              filename of the VASP (grid format) file to be
                            processed
 
    optional arguments:
      -h, --help            show this help message and exit
      -p {x,y,z}, --projection {x,y,z}
                            output averaged projection perpendicular to [x,y,z]
