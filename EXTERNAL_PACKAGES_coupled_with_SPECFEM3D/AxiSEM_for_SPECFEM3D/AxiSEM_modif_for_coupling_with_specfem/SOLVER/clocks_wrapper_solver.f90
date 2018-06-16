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
!> Module to do the SOLVER-specific initialization of the clocks module
!! i.e. defining the SOLVER-specific clocks
module clocks_wrapper_solver

  implicit none
  private

  public                   :: start_clock, end_clock

contains

!-----------------------------------------------------------------------------------------
!> Driver routine to start the timing, using the clocks_mod module.
subroutine start_clock

  use clocks_mod, only: clocks_init, clock_id
  use data_time, only: idcomm, iddump, idmpi, idmpiws, idmpiwf, idnbio, idold, &
                         idstiff, idanelts, idanelst
  use data_proc, only: lpr, mynum
  use data_io, only: verbose

  implicit none

  character(len=8)  :: mydate
  character(len=10) :: mytime

  call date_and_time(mydate,mytime)
  if (lpr) write(*,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Simulation started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)

  if (verbose > 1) write(69,11) mydate(5:6), mydate(7:8), mydate(1:4), &
                                mytime(1:2), mytime(3:4)

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
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Wapper routine to end timing and display clock informations.
subroutine end_clock

  use clocks_mod, only: clocks_exit
  use data_proc, only: mynum

  implicit none

  if (mynum == 0) then
     write(*,*)
     write(*,"(10x,'Summary of timing measurements:')")
     write(*,*)
  endif

  call clocks_exit(mynum)

  if (mynum == 0) write(*,*)

end subroutine end_clock
!-----------------------------------------------------------------------------------------

end module clocks_wrapper_solver
!=========================================================================================
