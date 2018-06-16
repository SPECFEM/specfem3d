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
program axisem

  use data_proc, only: nproc, mynum, appnproc, appmynum, lpr, procstrg
  use data_io, only: dump_xdmf, use_netcdf, verbose
  use nc_routines, only: nc_end_output, nc_finish_prepare
  use nc_snapshots, only: nc_close_snapfile
  use data_source, only: isim
  use data_mesh, only: do_mesh_tests
  use parameters, only: open_local_output_file, readin_parameters, &
                                    read_inparam_basic_verbosity
  use get_mesh, only: read_db
  use def_grid, only: init_grid, mesh_tests, deallocate_preloop_arrays
  use time_evol_wave, only: prepare_waves, time_loop
  use commun, only: pinit, pend, barrier
  use meshes_io, only: finish_xdmf_xml
  use data_io, only: verbose, define_io_appendix
  use clocks_wrapper_solver, only: start_clock, end_clock

  !!!! Modifs SB
  use coupling_mod
  !!!! End modifs

  implicit none

  call set_ftz() ! ftz.c, set flush to zero
  call pinit ! commun
  call read_inparam_basic_verbosity ! parameters
  if (lpr .and. verbose >= 1) write(*,'(/,a,/)') ' MAIN: Welcome to AxiSEM!'

  call define_io_appendix(appmynum, mynum)
  call define_io_appendix(appnproc, nproc)

  call open_local_output_file ! parameters, open file for processor-specific screen output
  call start_clock !clocks_wrapper_solver

  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Reading parameters..................................'
  call readin_parameters ! parameters

  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Reading mesh database...............................'
  call read_db  ! get_mesh

  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Initializing grid...................................'
  call init_grid ! def_grid

  if (do_mesh_tests) then
     if (lpr .and. verbose >= 1) &
         write(*,*) 'MAIN: Testing the mesh....................................'
     call mesh_tests ! def_grid
  endif

  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Starting wave preparation...........................'
  call prepare_waves ! time_evol_wave


!!!!  BEGIN TEST MODIFS COUPLING (SB)
  if (coupling) then
     call initialize_coupling
     call barrier
  endif
!!!!! END MODFIS

  ! Deallocate all the large arrays that are not needed in the time loop,
  ! specifically those from data_mesh_preloop and data_pointwise
  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Deallocating arrays not needed in the time loop.....'
  call deallocate_preloop_arrays ! def_grid

  if (use_netcdf) then
     if (lpr .and. verbose >= 1) &
        write(*,*) 'MAIN: Finish preparation of NetCDF file...................'
     call nc_finish_prepare
  endif

  call barrier ! Just making sure we're all ready to rupture...

  if (lpr .and. verbose >= 1) &
     write(*,*) 'MAIN: Starting wave propagation...........................'
  call time_loop ! time_evol_wave

  if (use_netcdf) then
     if (lpr .and. verbose >= 1) &
        write(*,*) 'MAIN: Flush and close all netcdf files ...................'
     call nc_end_output ! Dump receiver seismograms to finalize netcdf output
  endif

  if (dump_xdmf) then
     if (lpr .and. verbose >= 1) write(*,*)'MAIN: Finishing xdmf xml file...'
     call finish_xdmf_xml ! meshes_io
     call nc_close_snapfile ! nc_snapshots
  endif

  call end_clock ! clocks_wrapper_solver

  call pend ! commun

!!!! SB coupling
  if (coupling) then
     if (lpr) call finalize_coupling
  endif
!!!! SB


  if (lpr)         write(*,*)  ' ==  ==  ==  == =PROGRAM axisem FINISHED ==  ==  ==  ==  ==  == ='
  if (verbose > 1) write(69,*) ' ==  ==  ==  == =PROGRAM axisem FINISHED ==  ==  ==  ==  ==  == ='

end program axisem
!=========================================================================================
