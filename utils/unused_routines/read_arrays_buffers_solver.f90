!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine read_arrays_buffers_solver(myrank, &
     iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
     npoin2D_xi,npoin2D_eta, &
     NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer myrank

  integer npoin2D_xi,npoin2D_eta
  integer NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

  character(len=256) LOCAL_PATH

  integer, dimension(NPOIN2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NPOIN2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  integer npoin2D_xi_mesher,npoin2D_eta_mesher

  double precision xdummy,ydummy,zdummy

! processor identification
  character(len=256) prname

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolleft_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 350  continue
  read(IIN,*) iboolleft_xi(npoin2D_xi),xdummy,ydummy,zdummy
  if(iboolleft_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 350
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi = npoin2D_xi - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_xi read')
  close(IIN)

! read iboolright_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 360  continue
  read(IIN,*) iboolright_xi(npoin2D_xi),xdummy,ydummy,zdummy
  if(iboolright_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 360
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi = npoin2D_xi - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_xi read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '# of points in MPI buffers along xi npoin2D_xi = ', &
                                npoin2D_xi
    write(IMAIN,*) '# of array elements transferred npoin2D_xi*NDIM = ', &
                                npoin2D_xi*NDIM
    write(IMAIN,*)
  endif

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read 2-D addressing for summation between slices along eta with MPI

! read iboolleft_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 370  continue
  read(IIN,*) iboolleft_eta(npoin2D_eta),xdummy,ydummy,zdummy
  if(iboolleft_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 370
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta = npoin2D_eta - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_eta read')
  close(IIN)

! read iboolright_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 380  continue
  read(IIN,*) iboolright_eta(npoin2D_eta),xdummy,ydummy,zdummy
  if(iboolright_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 380
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta = npoin2D_eta - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_eta read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '# of points in MPI buffers along eta npoin2D_eta = ', &
                                npoin2D_eta
    write(IMAIN,*) '# of array elements transferred npoin2D_eta*NDIM = ', &
                                npoin2D_eta*NDIM
    write(IMAIN,*)
  endif

  end subroutine read_arrays_buffers_solver

