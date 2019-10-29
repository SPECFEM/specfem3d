#!/usr/bin/env python
"""
GEOCUBIT.py
this file is part of GEOCUBIT

Created by Emanuele Casarotti
Copyright (c) 2011 Istituto Nazionale di Geofisica e Vulcanologia

---------------------------------------------------------------------------
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
---------------------------------------------------------------------------

GEOCUBIT requires:

- CUBIT 15+ - www.cubit.sandia.gov
- python 2.7

- numpy 1.0+  - http://downloads.sourceforge.net/numpy

"""
from __future__ import print_function

import sys

print("GEOCUBIT")

# version info
python_major_version = sys.version_info[0]
python_minor_version = sys.version_info[1]
print("Python version: ","{}.{}".format(python_major_version,python_minor_version))

# in case importing menu fails due to import utilities errors to find,
# this will add the geocubitlib/ folder to the sys.path:
import geocubitlib
sys.path.append(geocubitlib.__path__[0])
#print(sys.path)
print('')

import geocubitlib.menu as menu
import geocubitlib.start as start

mpiflag, iproc, numproc, mpi = start.start_mpi()
if menu.build_surface or menu.build_volume or menu.meshing:
    cubit = start.start_cubit(init=True)
else:
    cubit = start.start_cubit()

if __name__ == '__main__':

    print('VERSION 4.1')

    # GEOMETRY
    if menu.build_surface:
        from geocubitlib import surfaces
        surfaces.surfaces()
    if menu.build_volume:
        from geocubitlib import volumes
        volumes.volumes()

    # MESHING
    if menu.meshing:
        from geocubitlib import mesh_volume
        mesh_volume.mesh()

    # MERGING and EXPORTING
    if menu.collect:
        from geocubitlib.exportlib import collect_new
        try:
            output = menu.output
        except:
            output = 'totalmesh_merged'
        output = output.upper()
        #
        collect_new(menu.cpuxmin, menu.cpuxmax, menu.cpuymin, menu.cpuymax,
                    menu.cpux, menu.cpuy, menu.cubfiles, menu.ckbound_method1,
                    menu.ckbound_method2, menu.merge_tolerance,
                    curverefining=menu.curverefining, outfilename=output,
                    qlog=menu.qlog, export2SPECFEM3D=menu.export2SPECFEM3D,
                    listblock=menu.listblock, listflag=menu.listflag,
                    outdir=menu.SPECFEM3D_output_dir, add_sea=menu.add_sea,
                    decimate=menu.decimate,
                    cpml=menu.cpml, cpml_size=menu.cpml_size,
                    top_absorbing=menu.top_absorbing, hex27=menu.hex27,
                    check_merging=menu.check_merging,
                    save_cubfile=menu.save_cubfile,
                    starting_tolerance=menu.starting_tolerance)

    if menu.export2SPECFEM3D and not menu.collect:
        from geocubitlib.exportlib import e2SEM
        print(menu.cubfiles)
        print('hex27 ', menu.hex27)
        e2SEM(files=menu.cubfiles,
              listblock=menu.listblock,
              listflag=menu.listflag,
              outdir=menu.SPECFEM3D_output_dir,
              hex27=menu.hex27)
