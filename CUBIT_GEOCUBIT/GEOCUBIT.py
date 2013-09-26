#!/usr/bin/env python
"""
GEOCUBIT.py
this file is part of GEOCUBIT

Created by Emanuele Casarotti 
Copyright (c) 2011 Istituto Nazionale di Geofisica e Vulcanologia 

---------------------------------------------------------------------------
 GEOCUBIT is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GEOCUBIT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GEOCUBIT.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------

GEOCUBIT requires:

- CUBIT 12.2 - www.cubit.sandia.gov
- python 2.5 
- numpy 1.0+  - http://downloads.sourceforge.net/numpy

  -- optional for parallel mesh
- pympi - http://downloads.sourceforge.net/pympi

"""
import geocubitlib.menu as menu
import geocubitlib.start as start
mpiflag,iproc,numproc,mpi   = start.start_mpi()
if menu.build_surface or menu.build_volume or menu.meshing:
    cubit                   = start.start_cubit(init=True)
else:
    cubit                   = start.start_cubit()

if __name__ == '__main__':
    
    #GEOMETRY
    if menu.build_surface:              
        from geocubitlib import surfaces
        surfaces.surfaces()                 
    if menu.build_volume:                   
        from geocubitlib import volumes     
        volumes.volumes()                   
               
    #MESHING
    if menu.meshing:
        from geocubitlib import mesh_volume
        mesh_volume.mesh()             
    
    #MERGING and EXPORTING
    if menu.collect:
        from geocubitlib.exportlib import collect
        try:
            output=menu.output
        except:
            output='totalmesh_merged'
        output=output.upper()
        #
        collect(menu.cpuxmin,menu.cpuxmax,menu.cpuymin,menu.cpuymax,menu.cpux,menu.cpuy,menu.cubfiles,menu.ckbound_method1,menu.ckbound_method2,menu.merge_tolerance,curverefining=menu.curverefining,outfilename=output,qlog=menu.qlog,export2SPECFEM3D=menu.export2SPECFEM3D,listblock=menu.listblock,listflag=menu.listflag,outdir=menu.SPECFEM3D_output_dir,add_sea=menu.add_sea,decimate=menu.decimate,cpml=menu.cpml,cpml_size=menu.cpml_size,top_absorbing=menu.top_absorbing,hex27=menu.hex27)
        
    if menu.export2SPECFEM3D and not menu.collect:
        from geocubitlib.exportlib import e2SEM
        print menu.cubfiles
        print 'hex27 ',menu.hex27
        e2SEM(files=menu.cubfiles,listblock=menu.listblock,listflag=menu.listflag,outdir=menu.SPECFEM3D_output_dir,hex27=menu.hex27)
