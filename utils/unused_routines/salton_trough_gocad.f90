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
!=======================================================================

subroutine read_salton_sea_model(vp_array)

  implicit none

  include 'constants.h'
  include 'constants_gocad.h'

  real :: vp_array(GOCAD_ST_NU,GOCAD_ST_NV,GOCAD_ST_NW)
  integer :: ios, reclen

  character(len=256) SALTON_SEA_MODEL_FILE

  reclen=(GOCAD_ST_NU * GOCAD_ST_NV * GOCAD_ST_NW) * 4
  call get_value_string(SALTON_SEA_MODEL_FILE, &
                        'model.SALTON_SEA_MODEL_FILE', &
                        'DATA/st_3D_block_harvard/regrid3_vel_p.bin')
  open(11,file=SALTON_SEA_MODEL_FILE,status='old',action='read',form='unformatted',access='direct',recl=reclen,iostat=ios)
  if (ios /= 0) then
    print *, 'iostat = ', ios
    stop 'Error opening file'
  endif
  read(11,rec=1,iostat=ios) vp_array
  if (ios /= 0) stop 'Error reading vp_array'
  close(11)

end subroutine read_salton_sea_model


subroutine vx_xyz2uvw(xmesh, ymesh, zmesh, uc, vc, wc)


  implicit none
  include 'constants.h'

  double precision :: xmesh, ymesh, zmesh, uc, vc, wc

  uc = (GOCAD_ST_NU-1) * ((xmesh -  GOCAD_ST_O_X) * GOCAD_ST_V_Y - (ymesh - GOCAD_ST_O_Y) * GOCAD_ST_V_X)  &
             / (GOCAD_ST_U_X * GOCAD_ST_V_Y - GOCAD_ST_U_Y * GOCAD_ST_V_X)
  vc = (GOCAD_ST_NV-1) * ((ymesh - GOCAD_ST_O_Y) - uc * GOCAD_ST_U_Y/(GOCAD_ST_NU-1) ) / GOCAD_ST_V_Y
  wc = (GOCAD_ST_NW-1) * (zmesh - GOCAD_ST_O_Z) / GOCAD_ST_W_Z

end subroutine vx_xyz2uvw


subroutine vx_xyz_interp(uc,vc,wc, vp, vs, rho, vp_array)

  implicit none
  include 'constants.h'

  double precision :: uc,vc,wc, vp, vs, rho
  real :: vp_array(GOCAD_ST_NU,GOCAD_ST_NV,GOCAD_ST_NW)

  integer :: i,j,k,ixi,ieta,iga
  real :: v1, v2, v3, v4, v5, v6, v7, v8, xi, eta, ga, vi
  double precision :: zmesh
  real,parameter :: eps = 1.0e-3


  i = uc + 1
  j = vc + 1
  k = wc + 1

  xi = uc + 1 - i
  eta = vc + 1- j
  ga = wc + 1 -k

  ixi = nint(xi)
  ieta = nint(eta)
  iga = nint(ga)

!  print *, 'gc = ', i, j, k
!  print *, 'xi, eta, ga = ', xi, eta, ga
!  print *, 'ixi, ieta, iga = ', ixi, ieta, iga


  if (i > 0 .or. i < GOCAD_ST_NU  .or. j > 0 .or. j < GOCAD_ST_NV .or. k > 0 .or. k < GOCAD_ST_NW) then
    v1 = vp_array(i,j,k)
    v2 = vp_array(i+1,j,k)
    v3 = vp_array(i+1,j+1,k)
    v4 = vp_array(i,j+1,k)
    v5 = vp_array(i,j,k+1)
    v6 = vp_array(i+1,j,k+1)
    v7 = vp_array(i+1,j+1,k+1)
    v8 = vp_array(i,j+1,k+1)
    vi = vp_array(i+ixi,j+ieta,k+iga)
!    print *, v1, v2, v3, v4, v5, v6, v7, v8
    if ((v1 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v2 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v3 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v4 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v5 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v6 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v7 - GOCAD_ST_NO_DATA_VALUE) > eps .and. &
               (v8 - GOCAD_ST_NO_DATA_VALUE) > eps )  then
      vp = dble(&
                 v1 * (1-xi) * (1-eta) * (1-ga) +&
                 v2 * xi * (1-eta) * (1-ga) +&
                 v3 * xi * eta * (1-ga) +&
                 v4 * (1-xi) * eta * (1-ga) +&
                 v5 * (1-xi) * (1-eta) * ga +&
                 v6 * xi * (1-eta) * ga +&
                 v7 * xi * eta * ga +&
                 v8 * (1-xi) * eta * ga)
    else if ((vi - GOCAD_ST_NO_DATA_VALUE) > eps) then
      vp = dble(vi)
!    else if ((v1 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v1)
!    else if ((v2 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v2)
!    else if ((v3 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v3)
!    else if ((v4 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v4)
!    else if ((v5 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v5)
!    else if ((v6 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v6)
!    else if ((v7 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v7)
!    else if ((v7 - GOCAD_ST_NO_DATA_VALUE) > eps) then
!      vp = dble(v8)
    else
      vp = GOCAD_ST_NO_DATA_VALUE
    endif
    zmesh = wc / (GOCAD_ST_NW - 1) * GOCAD_ST_W_Z + GOCAD_ST_O_Z
    if (zmesh > -8500.)  then
      vs = vp / (2 - (0.27*zmesh/(-8500)))
    else
      vs = vp/1.73
    endif
    if (vp > 2160.) then
      rho = vp/3 + 1280.
    else
      rho = 2000.
    endif
  else
    rho = GOCAD_ST_NO_DATA_VALUE
    vp = GOCAD_ST_NO_DATA_VALUE
    vs = GOCAD_ST_NO_DATA_VALUE
  endif

end subroutine vx_xyz_interp

