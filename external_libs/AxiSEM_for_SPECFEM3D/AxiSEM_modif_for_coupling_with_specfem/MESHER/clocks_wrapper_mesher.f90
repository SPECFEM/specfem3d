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
!> Module to do the MESHER-specific initialization of the clocks module
!! i.e. defining the MESHER-specific clocks
module clocks_wrapper_mesher

  implicit none
  private

  public                   :: start_clock, end_clock

contains

!-----------------------------------------------------------------------------------------
!> Driver routine to start the timing, using the clocks_mod module.
subroutine start_clock

  use data_time
  use clocks_mod, only: clock_id, clocks_init

  implicit none

  character(len=8)    :: mydate
  character(len=10)   :: mytime

  call date_and_time(mydate,mytime)
  write(*,11) mydate(5:6), mydate(7:8), mydate(1:4), mytime(1:2), mytime(3:4)

11 format('     Meshing started on ', A2,'/',A2,'/',A4,' at ', A2,'h ',A2,'min',/)


  call clocks_init(0)

  idold01 = clock_id('mergesort')
  idold05 = clock_id('mergesort - ind')
  idold02 = clock_id('generate_skeleton')
  idold03 = clock_id('create_gllmesh')
  idold04 = clock_id('define_global_global_numbering')
  idold06 = clock_id('bkgrdmodel_testing')
  idold08 = clock_id('glob-slob/flob numbering')
  idold07 = clock_id('get_global no loop')
  idold09 = clock_id('get_global in loop')
  idold11 = clock_id('create_pdb')

  idold12 = clock_id('define_glocal_numbering')
  idold13 = clock_id('define_sflocal_numbering')
  idold14 = clock_id('generate_serendipity_per_proc')

end subroutine start_clock
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Wapper routine to end timing and display clock informations.
subroutine end_clock

  use clocks_mod, only: clocks_exit

  implicit none

  write(*,*)
  write(*,"(10x,'Summary of timing measurements:')")
  write(*,*)

  call clocks_exit(0)


end subroutine end_clock
!-----------------------------------------------------------------------------------------

end module clocks_wrapper_mesher
!=========================================================================================
