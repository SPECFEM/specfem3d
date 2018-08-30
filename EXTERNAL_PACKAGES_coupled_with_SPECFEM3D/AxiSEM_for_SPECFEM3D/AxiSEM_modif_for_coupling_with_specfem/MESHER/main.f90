!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
program gllmesh

  use data_grid
  use data_bkgrdmodel, only: have_fluid,have_solid

  use meshgen, only: generate_skeleton    ! creates mesh skeleton
  use gllmeshgen                                    ! creates complete GLL mesh
  use input
  use numbering
  use pdb, only: create_pdb
  use test_bkgrdmodel
  use discont_meshing
  use mesh_info
  use parallelization
  use data_mesh

  use data_time
  use clocks_mod, only: tick
  use clocks_wrapper_mesher, only: start_clock, end_clock

  implicit none

  call start_clock !clocks

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Need to load/predefine:
  ! 1) points per wavelength
  ! 2) dominant period
  ! 3) bkgrdmodel & discontinuity radii
  call read_params ! input
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  write(*,*)'MAIN: creating subregions/discontinuity model..........'; call flush(6)
  call create_subregions ! discont_meshing

  southern = .true.

  write(*,*)'MAIN: generating skeleton..............................'; call flush(6)
  iclock02 = tick()
  call generate_skeleton ! meshgen
  iclock02 = tick(id=idold02, since=iclock02)

  write(*,*)'MAIN: creating gllmesh.................................'; call flush(6)
  iclock03 = tick()
  call create_gllmesh ! gllmeshgen
  iclock03 = tick(id=idold03, since=iclock03)
  !call test_mapping   ! gllmeshgen

  write(*,*)'MAIN: glob-glob numbering..............................'; call flush(6)
  iclock04 = tick()
  call define_global_global_numbering ! numbering
  iclock04 = tick(id=idold04, since=iclock04)

  write(*,*)'MAIN: defining regions.................................'; call flush(6)
  call define_regions ! mesh_info

  write(*,*)'MAIN: test model.......................................'; call flush(6)
  iclock06 = tick()
  call bkgrdmodel_testing ! test_bkgrdmodel
  iclock06 = tick(id=idold06, since=iclock06)

  ! Here starts the distinction between solid and fluid regions
  write(*,*)'MAIN: define subregions................................'; call flush(6)
  call def_fluid_regions ! mesh_info
  call def_solid_regions ! mesh_info
  call extract_fluid_solid_submeshes ! gllmeshgen

  write(*,*)'MAIN: glob-slob/flob numbering.........................'; call flush(6)
  iclock08 = tick()
  if (have_fluid) &
       call define_global_flobal_numbering  ! numbering
  if (have_solid) &
       call define_global_slobal_numbering  ! numbering
  iclock08 = tick(id=idold08, since=iclock08)

  ! Boundary matrices: find corresponding element neighbors
  write(*,*)'MAIN: boundaries.......................................'; call flush(6)
  call define_boundaries  ! mesh_info

  ! Parallelization
  write(*,*)'MAIN: domain decomposition.............................'; call flush(6)
  call create_domain_decomposition !parallelization

  write(*,*)'MAIN: creating parallel database.......................'; call flush(6)
  iclock11 = tick()
  call create_pdb ! pdb
  iclock11 = tick(id=idold11, since=iclock11)

  ! clean up
  call empty_data_mesh

  call end_clock ! clocks

  print *
  write(*,*)'....DONE WITH MESHER !'

end program gllmesh
!=========================================================================================
