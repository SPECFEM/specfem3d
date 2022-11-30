GEOCUBIT
========

GeoCubit is a Python library wrapping around the Cubit Python Interface.

Main author: Emanuele Casarotti, INGV, Roma, Italy.

Copyright (c) 2011 and 2017, Istituto Nazionale di Geofisica e Vulcanologia.

Improvements by: Elliott Sales de Andrade, Department of Physics, University of Toronto, Canada.


## Description

It aims at facilitating the meshing process in some common problems encountered in seismic or more generally acoustic wave propagation.
In particular, it is focused on the meshing requests of Specfem3D and it is helpful for such tedious tasks as:

• Creation of geophysical surfaces and volumes (ex. topography).

• Mesh of layered volumes with hexahedral.

• Creation of an anisotropic mesh suitable for cases where some alluvial basins (or slow velocity zones) are present.

• It can be used as a serial or parallel process. The parallel meshing capabilities are fun- damental for large geophysical problems (ex. mesh of Southern California using SRTM topography).

GeoCubit can be used inside the graphical interface of Trelis (i.e. as a Python object in the script tab) or as Unix command.

## Requirements

GeoCubit requires:
- Coreform Cubit 15+ (https://coreform.com/products/coreform-cubit/, see also cubit.sandia.gov)
- python 2.7
- numpy 1.0+

It is possible that GEOCUBIT continues to work for CUBIT 15+ but it is not supported anymore.

## Installation & Configuration

On Linux/Unix system: in order to have the possibility to import the cubit libray included in Cubit in a **python script** you should set up the path for Cubit

```
export CUBITDIR=$CUBITHOME
export CUBITLIB=$CUBITDIR/bin:$CUBITDIR/structure
export PYTHONPATH=$PYTHONPATH:$CUBITDIR/bin:$CUBITDIR/structure
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUBITDIR/bin
export PATH=$PATH:$CUBITDIR/bin
```

You could install GEOCUBIT

`python setup.py install`

or

```
export PYTHONPATH=$PYTHONPATH:[YourGeocubitDir]
export PATH=$PATH:[YourGeocubitDir]
```



On MacOSX 10.10, Trelis points the system version of python2.7 (/Library/Python/2.7/site-packages)
so the version of python required is `/usr/bin/python2.7`

You could install GEOCUBIT with administrator privilege (for running script in the python tab of the Cubit GUI)

`/usr/bin/python2.7 setup.py install --install-lib=/Library/Python/2.7/site-packages/`

and in addition if you want to run meshing scripts in python (outside of the GUI)   

```
export CUBITDIR=$CUBITHOME
export PYTHONPATH=$CUBITDIR/contents/MacOS:$PYTHONPATH
export CUBITLIB=$CUBITDIR/contents/MacOS
export DYLD_LIBRARY_PATH=$CUBITDIR/contents/MacOS:${DYLD_LIBRARY_PATH}
export PATH=$PATH:$CUBITDIR:$CUBITDIR/Contents/MacOS
```

The meshing scripts works for `/usr/bin/python2.7`.

A simple "hello world" meshing script is:

```
import cubit
cubit.init([""])
cubit.cmd('brick x 10')
cubit.cmd('mesh vol all')
```

For defining the absorbing boundary surfaces in a 4 (near parallel) side meshed volume:

```
from geocubitlib import exportlib
exportlib.collect(outdir='.',outfilename='mymesh')
```

The mesh will be exported in exodus format with name defined by the optional parameter `outfilename`:

For saving the mesh in Specfem3D format:

```
exportlib.e2SEM(outdir='.')
```

#### VERY IMPORTANT:
**Note from Emanuele Casarotti about how to set up the CUBIT (now called TRELIS) and GEOCUBIT meshing packages:**

you should set the pythonpath as in the example "CUBIT_GEOCUBIT/examples/homogeneous_halfspace" (basically installing geocubit....)

CUBIT 14 for mac
```
python2.5 setup.py install --install-lib=/Library/Python/2.5/site-packages/
```

TRELIS 14 for mac
```
python2.6 setup.py install --install-lib=/Library/Python/2.6/site-packages/
```
on a Linux machine

check that the version of cubit 32/64 bit matches the python version...

in order to have the possibility to import cubit in a python script.

Set CUBITHOME to the directory where you installed CUBIT
```
export CUBITDIR=$CUBITHOME
export CUBITLIB=$CUBITDIR/bin:$CUBITDIR/structure
export PYTHONPATH=$PYTHONPATH:$CUBITDIR/bin:$CUBITDIR/structure
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUBITDIR/bin
export PATH=$PATH:$CUBITDIR/bin
```
type "source setpaths.sh" in CUBIT_GEOCUBIT/ directory



## News / Updates


**Update by Emanuele Casarotti, INGV, Roma, Italy, January 2017**:

I have finally concluded the review of geocubit for compatibility
issue with the last version of Trelis 16.
I have submitted to the repository the changes and I have tested the examples.
I'm waitinfg for travis.

As reported there some examples that fails: the CPML and the
tomography (we need to modify the starting toy model).

This is the list of updates made on Jan 5, 2017:

- compatibility with Coreform Cubit 15+ and 16+

- new merging chunk function (faster and more stable) for Coreform Cubit 15+ (before it was only for cubit 12.2)

- pep8 beautification

- bug fixes (including hex27 support in addition to hex8)

- alpha version of netcdf mesh for specfem3d



## License

```
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
```
