!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

subroutine read_sea99_s_model(SEA99M_V)

  implicit none

  include "constants.h"

! sea99_s_model_variables
  type sea99_s_model_variables
    sequence
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
 end type sea99_s_model_variables

 type (sea99_s_model_variables) SEA99M_V
! sea99_s_model_variables

  integer :: i,ia,io,j

!----------------------- choose input file:  ------------------
! relative anomaly


  open(1,file='DATA/Lebedev_sea99/sea99_dvsvs')

!----------------------- read input file:  ------------------

  do i = 1, 6
     read(1,*)
  enddo
  read(1,*) SEA99M_V%sea99_ndep
  read(1,*) (SEA99M_V%sea99_depth(i), i = 1, SEA99M_V%sea99_ndep)
  read(1,*)
  read(1,*) SEA99M_V%alatmin, SEA99M_V%alatmax
  read(1,*) SEA99M_V%alonmin, SEA99M_V%alonmax
  read(1,*) SEA99M_V%sea99_ddeg,SEA99M_V%sea99_nlat,SEA99M_V%sea99_nlon
  if (SEA99M_V%sea99_nlat /= nint((SEA99M_V%alatmax-SEA99M_V%alatmin)/SEA99M_V%sea99_ddeg)+1) then
     stop 'alatmin,alatmax,sea99_nlat'
  endif
  if (SEA99M_V%sea99_nlon /= nint((SEA99M_V%alonmax-SEA99M_V%alonmin)/SEA99M_V%sea99_ddeg)+1) then
     stop 'alonmin,alonmax,sea99_nlon'
  endif
  read(1,*)
  do j = 1, SEA99M_V%sea99_ndep
     do ia = 1, SEA99M_V%sea99_nlat
        read (1,*) (SEA99M_V%sea99_vs(ia,io,j), io = 1, SEA99M_V%sea99_nlon)
     enddo
  enddo

end subroutine read_sea99_s_model

subroutine sea99_s_model(radius,theta,phi,dvs,SEA99M_V)

  implicit none

  include "constants.h"

! sea99_s_model_variables
  type sea99_s_model_variables
     sequence
     integer :: sea99_ndep
     integer :: sea99_nlat
     integer :: sea99_nlon
     double precision :: sea99_ddeg
     double precision :: alatmin
     double precision :: alatmax
     double precision :: alonmin
     double precision :: alonmax
     double precision :: sea99_vs(100,100,100)
     double precision :: sea99_depth(100)
  end type sea99_s_model_variables

  type (sea99_s_model_variables) SEA99M_V
! sea99_s_model_variables

  integer :: id1,i,ilat,ilon
  double precision :: alat1,alon1,radius,theta,phi,dvs
  double precision :: xxx,yyy,dep,pla,plo,xd1,dd1,dd2,ddd(2)
 !----------------------- depth in the model ------------------
  dep=R_EARTH_KM*(R_UNIT_SPHERE - radius)
  pla=90.0d0 - theta/DEGREES_TO_RADIANS
  plo=phi/DEGREES_TO_RADIANS
  if (dep .le. SEA99M_V%sea99_depth(1)) then
     id1 = 1
     xd1 = 0
  else if (dep .ge. SEA99M_V%sea99_depth(SEA99M_V%sea99_ndep)) then
     id1 = SEA99M_V%sea99_ndep
     xd1 = 0
  else
     do i = 2, SEA99M_V%sea99_ndep
        if (dep .le. SEA99M_V%sea99_depth(i)) then
           id1 = i-1
           xd1 = (dep-SEA99M_V%sea99_depth(i-1)) / (SEA99M_V%sea99_depth(i) - SEA99M_V%sea99_depth(i-1))
           go to 1
        endif
     enddo
  endif
1 continue

!----------------------- value at a point ---------------------
!----- approximate interpolation, OK for the (dense) 1-degree sampling ------

  ilat = int((pla - SEA99M_V%alatmin)/SEA99M_V%sea99_ddeg) + 1
  ilon = int((plo - SEA99M_V%alonmin)/SEA99M_V%sea99_ddeg) + 1
  alat1 = SEA99M_V%alatmin + (ilat-1)*SEA99M_V%sea99_ddeg
  alon1 = SEA99M_V%alonmin + (ilon-1)*SEA99M_V%sea99_ddeg

  do i = 1, 2
     xxx = (pla-alat1)/SEA99M_V%sea99_ddeg
     yyy = SEA99M_V%sea99_vs(ilat+1,ilon,id1+i-1)-SEA99M_V%sea99_vs(ilat,ilon,id1+i-1)
     dd1 = SEA99M_V%sea99_vs(ilat,ilon,id1+i-1) + yyy*xxx
     yyy = SEA99M_V%sea99_vs(ilat+1,ilon+1,id1+i-1)-SEA99M_V%sea99_vs(ilat,ilon+1,id1+i-1)
     dd2 = SEA99M_V%sea99_vs(ilat,ilon+1,id1+i-1) + yyy*xxx
     xxx = (plo-alon1)/SEA99M_V%sea99_ddeg
     yyy = dd2 - dd1
     ddd(i) = dd1 + yyy*xxx
  enddo
  dvs = ddd(1) + (ddd(2)-ddd(1)) * xd1
  if(dvs>1.d0) dvs=0.0d0

end subroutine sea99_s_model


