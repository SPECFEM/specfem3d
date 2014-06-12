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

GeoCubit can be used inside the graphical interface of Cubit (i.e. as a Python object in the script tab) or as Unix command.

See also Cubit (now also sold as "Trelis") here [cubit.sandia.gov]

