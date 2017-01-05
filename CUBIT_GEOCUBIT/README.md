GEOCUBIT
======================

GeoCubit is a Python library wrapping around the Cubit Python Interface.

Main author: Emanuele Casarotti, INGV, Roma, Italy.

Copyright (c) 2011 Istituto Nazionale di Geofisica e Vulcanologia.

Improvements by: Elliott Sales de Andrade, Department of Physics, University of Toronto, Canada.

```
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
```

It aims at facilitating the meshing process in some common problems encountered in seismic or more generally acoustic wave propagation.
In particular, it is focused on the meshing requests of Specfem3D and it is helpful for such tedious tasks as:

• Creation of geophysical surfaces and volumes (ex. topography).

• Mesh of layered volumes with hexahedral.

• Creation of an anisotropic mesh suitable for cases where some alluvial basins (or slow velocity zones) are present.

• It can be used as a serial or parallel process. The parallel meshing capabilities are fun- damental for large geophysical problems (ex. mesh of Southern California using SRTM topography).

GeoCubit can be used inside the graphical interface of Trelis (i.e. as a Python object in the script tab) or as Unix command.

## Requirements:

- Trelis 15+ (www.csimsoft.com, see also cubit.sandia.gov)
- python 2.7
- numpy 1.0+

It is possible that GEOCUBIT continues to work for CUBIT 15+ but it is not supported anymore.

## Configuration:

On Linux/Unix system: in order to have the possibility to import the cubit libray included in Trelis in a **python script** you should set up the path for CUBIT/Trelis

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

You could install GEOCUBIT with administrator privilege (for running script in the python tab of the Trelis GUI)

`/usr/bin/python2.7 setup.py install --install-lib=/Library/Python/2.7/site-packages/`

and in addition if you want to run meshing scripts in python (outside of the GUI)   

```
export CUBITDIR=$CUBITHOME #for example /Applications/Trelis-16.0.app
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

