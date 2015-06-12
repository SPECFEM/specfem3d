!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
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
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!===================
program axisem 
!===================

  use data_proc,      only : nproc, mynum, appnproc, appmynum, lpr, procstrg
  use data_io,        only : dump_xdmf, use_netcdf, verbose
  use nc_routines,    only : nc_end_output, nc_finish_prepare
  use data_source,    only : isim,num_simul
  use data_mesh,      only : do_mesh_tests
  use parameters,     only : open_local_output_file, readin_parameters, &
                             read_inparam_basic_verbosity
  use get_mesh,       only : read_db 
  use def_grid,       only : init_grid, mesh_tests, deallocate_preloop_arrays
  use time_evol_wave, only : prepare_waves, time_loop
  use commun,         only : pinit, pend, barrier
  use meshes_io,      only : finish_xdmf_xml
  use data_io,        only : verbose
  
  !!!!! MODIFS TEST COUPLING (SB)
  use coupling_mod
  !!!!! END MODIFS

  implicit none

  call set_ftz() ! ftz.c, set flush to zero
  call pinit ! commun
  call read_inparam_basic_verbosity ! parameters
  if (lpr .and. verbose >= 1) write(6,'(/,a,/)') ' MAIN: Welcome to AxiSEM!'

  call define_io_appendix(appmynum,mynum)
  call define_io_appendix(appnproc,nproc)
  
  !call get_mesh_params !Get very basic mesh params, including pol. order, needed later
  call open_local_output_file ! parameters, open file for processor-specific screen output 
  call start_clock !clocks

  if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Reading parameters..................................'
  call readin_parameters ! parameters
  
  if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Reading mesh database...............................'
  call read_db  ! get_mesh

  if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Initializing grid...................................'
  call init_grid ! def_grid

  if (do_mesh_tests) then
     if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Testing the mesh....................................'
     call mesh_tests ! def_grid
  endif 

  if (num_simul .ne. 1) then
     write(6,*) 'ERROR: implementation of multiple simulations within one run'
     write(6,*) '       not finished yet.'
     write(6,*) '       For now set number of simulations in inparam to 1, splitting '
     write(6,*) '       to the different sources is then done by the submit script.'
     stop 2
  endif


  do isim=1, num_simul

     if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Starting wave preparation...........................'
     call prepare_waves ! time_evol_wave

!!!!! BEGIN TEST MODIFS COUPLING (SB)
     call read_boundary_coordinates
     call barrier
     !stop
!!!!! END MODFIS



     ! Deallocate all the large arrays that are not needed in the time loop,
     ! specifically those from data_mesh_preloop and data_pointwise
     if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Deallocating arrays not needed in the time loop.....'
     call deallocate_preloop_arrays ! def_grid
 
     if (use_netcdf) then
        if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Finish preparation of NetCDF file...................'
        call nc_finish_prepare
     endif
    
     call barrier ! Just making sure we're all ready to rupture...
     
     if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Starting wave propagation...........................'
     call time_loop ! time_evol_wave
  enddo

  if (use_netcdf) then
     if (lpr .and. verbose >= 1) write(6,*) 'MAIN: Flush and close all netcdf files ...................'
     call nc_end_output ! Dump receiver seismograms to finalize netcdf output 
  end if
  
  if (dump_xdmf) then
     if (lpr .and. verbose >= 1) write(6,*)'MAIN: Finishing xdmf xml file...'
     call finish_xdmf_xml() ! meshes_io
  endif

  call end_clock ! clocks

  call pend ! commun

  !! VM VM coupling
  if (lpr) call finalize_coupling 
  !! VM VM

  if(lpr) write(6,*) '=========PROGRAM axisem FINISHED============='
  if (verbose > 1) write(69,*) '=========PROGRAM axisem FINISHED============='

!=======================
end program axisem
!=======================


!-----------------------------------------------------------------------------
subroutine start_clock
  !
  ! Driver routine to start the timing, using the clocks_mod module.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use data_time,  only : idcomm, iddump, idmpi, idmpiws, idmpiwf, idnbio, idold, &
                         idstiff, idanelts, idanelst
  use data_proc,  only : lpr, mynum
  use clocks_mod, only : clock_id, clocks_init
  use data_io,    only : verbose
  
  implicit none
  
  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime) 
  if (lpr) write(6,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Simulation started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)

  if (verbose > 1) write(69,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

  if (verbose > 1) then
      call clocks_init(mynum)
  else
      call clocks_init()
  endif

  idold    = clock_id('Time loop routine')
  idcomm   = clock_id('Assembly/MPI routines')
  idmpi    = clock_id(' > Only MPI routine')
  idmpiws  = clock_id(' > Only solid MPI_WAIT')
  idmpiwf  = clock_id(' > Only fluid MPI_WAIT')
  idstiff  = clock_id('Stiffness routine')
  idanelst = clock_id(' > Anelastic stiffness routine')
  idanelts = clock_id('Anelastic time step routine')
  iddump   = clock_id('Dump routine')
  idnbio   = clock_id('Non Blocking IO red light')

end subroutine start_clock
!=============================================================================

!-----------------------------------------------------------------------------
subroutine end_clock 
  !
  ! Wapper routine to end timing and display clock informations.
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  use clocks_mod, only : clocks_exit
  use data_proc,  only : mynum

  implicit none

  if (mynum==0) then
     write(6,*)
     write(6,"(10x,'Summary of timing measurements:')")
     write(6,*)
  endif

  call clocks_exit(mynum)

  if (mynum==0) write(6,*)

end subroutine end_clock
!=============================================================================

!-----------------------------------------------------------------------------
subroutine define_io_appendix(app,iproc)
  !
  ! Defines the 4 digit character string appended to any 
  ! data or io file related to process myid. 
  !
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  implicit none
  integer, intent(in)           :: iproc
  character(len=4), intent(out) :: app
  
  write(app,"(I4.4)") iproc

end subroutine define_io_appendix
!=============================================================================
