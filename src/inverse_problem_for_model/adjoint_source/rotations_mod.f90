!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

module rotations_mod

  use interpolation_mod, only: si, sp, di, dp, cp, hp, deg2rad, rad2deg
  implicit none


  real(kind=dp), dimension(3,3), private :: rotmat, rotmat_t
  real(kind=dp), dimension(3,3), private :: rotc, rotc_t
  real(kind=dp), parameter,      private :: rt=6371000._dp

contains

!================================================================================
! Define rotation matrix for mesh
  subroutine define_mesh_rotation_matrix(lat0,lon0,azi0)

    real(kind=dp),    intent(in) :: lat0, lon0, azi0
    real(kind=dp)                :: lat,   lon,  azi

    !* Pass in radian
    lat = deg2rad * lat0
    lon = deg2rad * lon0
    azi = deg2rad * azi0

    !* Define mesh rotation matrix pass Cartesian geocentric coordinates
    !                              from real mesh position    ( lat0, lon0, azi0)
    !                              to reference mesh position (lat=0,lon=0,azi=0)
    rotmat(1,1) =  cos(lat) * cos(lon)
    rotmat(1,2) =  cos(lat) * sin(lon)
    rotmat(1,3) =  sin(lat)
    rotmat(2,1) = -sin(lon) * cos(azi) - sin(azi) * sin(lat) * cos(lon)
    rotmat(2,2) =  cos(lon) * cos(azi) - sin(azi) * sin(lat) * sin(lon)
    rotmat(2,3) =  sin(azi) * cos(lat)
    rotmat(3,1) =  sin(lon) * sin(azi) - cos(azi) * sin(lat) * cos(lon)
    rotmat(3,2) = -cos(lon) * sin(azi) - cos(azi) * sin(lat) * sin(lon)
    rotmat(3,3) =  cos(azi) * cos(lat)

    !* Define reverse rotation
    rotmat_t    = transpose(rotmat)

  end subroutine define_mesh_rotation_matrix
!--------------------------------------------------------------------------------

!================================================================================
! Rotate data from mesh coordinates (vx,vy,vz) to earth coordinates (vn, ve, vz)
  subroutine rotate_comp_mesh2glob(vx, vy, vz, stalat, stalon, nt, nsta, vz2, vn, ve)

    integer(kind=si), intent(in) :: nt, nsta
    real(kind=cp),  dimension(nsta),   intent(in) :: stalat, stalon

    real(kind=cp), dimension(nsta,nt),  intent(in) :: vx, vy, vz
    real(kind=cp), dimension(nsta,nt), intent(out) :: vn, ve, vz2

    real(kind=dp), dimension(nsta,nt) :: X1, Y1, Z1
    real(kind=dp), dimension(nsta,nt) :: X2, Y2, Z2

    real(kind=dp), dimension(nsta) :: lat, lon

    integer(kind=si) :: ista

    !* Pass in radian
    lat = deg2rad * stalat
    lon = deg2rad * stalon

    do ista=1,nsta

       !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
       X2(ista,:) = vz(ista,:)
       Y2(ista,:) = vx(ista,:)
       Z2(ista,:) = vy(ista,:)

       !* From local mesh to global eath
       X1(ista,:) = rotmat_t(1,1)*X2(ista,:) + rotmat_t(1,2)*Y2(ista,:) &
            + rotmat_t(1,3)*Z2(ista,:)
       Y1(ista,:) = rotmat_t(2,1)*X2(ista,:) + rotmat_t(2,2)*Y2(ista,:) &
            + rotmat_t(2,3)*Z2(ista,:)
       Z1(ista,:) = rotmat_t(3,1)*X2(ista,:) + rotmat_t(3,2)*Y2(ista,:) &
            + rotmat_t(3,3)*Z2(ista,:)

       !* Define rotation matrix
       rotc(1,1) =  cos(lat(ista)) * cos(lon(ista))
       rotc(1,2) =  cos(lat(ista)) * sin(lon(ista))
       rotc(1,3) =  sin(lat(ista))
       rotc(2,1) = -sin(lat(ista)) * cos(lon(ista))
       rotc(2,2) = -sin(lat(ista)) * sin(lon(ista))
       rotc(2,3) =  cos(lat(ista))
       rotc(3,1) = -sin(lon(ista))
       rotc(3,2) =  cos(lon(ista))
       rotc(3,3) =  0._dp

       !* Data in geographic coordinate at real position
       vz2(ista,:) = rotc(1,1)*X1(ista,:) + rotc(1,2)*Y1(ista,:) &
            + rotc(1,3)*Z1(ista,:)
       vn(ista,:)  = rotc(2,1)*X1(ista,:) + rotc(2,2)*Y1(ista,:) &
            + rotc(2,3)*Z1(ista,:)
       ve(ista,:)  = rotc(3,1)*X1(ista,:) + rotc(3,2)*Y1(ista,:) &
            + rotc(3,3)*Z1(ista,:)

    enddo

  end subroutine rotate_comp_mesh2glob
!--------------------------------------------------------------------------------


!================================================================================
! Rotate data from mesh coordinates (vx,vy,vz) to earth coordinates (vn, ve, vz)
  subroutine rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)

    integer(kind=si), intent(in) :: nt, nsta
    real(kind=cp),  dimension(nsta),   intent(in) :: stalat, stalon

    real(kind=cp), dimension(nsta,nt), intent(out) :: vx, vy, vz
    real(kind=cp), dimension(nsta,nt),  intent(in) :: vn, ve, vz2

    real(kind=dp), dimension(nsta,nt) :: X1, Y1, Z1
    real(kind=dp), dimension(nsta,nt) :: X2, Y2, Z2

    real(kind=dp), dimension(nsta) :: lat, lon

    integer(kind=si) :: ista

    !* Pass in radian
    lat = deg2rad * stalat
    lon = deg2rad * stalon

    do ista=1,nsta

       !* Define rotation matrix
       rotc(1,1) =  cos(lat(ista)) * cos(lon(ista))
       rotc(1,2) =  cos(lat(ista)) * sin(lon(ista))
       rotc(1,3) =  sin(lat(ista))
       rotc(2,1) = -sin(lat(ista)) * cos(lon(ista))
       rotc(2,2) = -sin(lat(ista)) * sin(lon(ista))
       rotc(2,3) =  cos(lat(ista))
       rotc(3,1) = -sin(lon(ista))
       rotc(3,2) =  cos(lon(ista))
       rotc(3,3) =  0._dp
       rotc_t = transpose(rotc)

       !* Data in geographic coordinate at real position
       X1(ista,:) = rotc_t(1,1)*vz2(ista,:) + rotc_t(1,2)*vn(ista,:) &
            + rotc_t(1,3)*ve(ista,:)
       Y1(ista,:) = rotc_t(2,1)*vz2(ista,:) + rotc_t(2,2)*vn(ista,:) &
            + rotc_t(2,3)*ve(ista,:)
       Z1(ista,:) = rotc_t(3,1)*vz2(ista,:) + rotc_t(3,2)*vn(ista,:) &
            + rotc_t(3,3)*ve(ista,:)

       !* From local mesh to global eath
       X2(ista,:) = rotmat(1,1)*X1(ista,:) + rotmat(1,2)*Y1(ista,:) &
            + rotmat(1,3)*Z1(ista,:)
       Y2(ista,:) = rotmat(2,1)*X1(ista,:) + rotmat(2,2)*Y1(ista,:) &
            + rotmat(2,3)*Z1(ista,:)
       Z2(ista,:) = rotmat(3,1)*X1(ista,:) + rotmat(3,2)*Y1(ista,:) &
            + rotmat(3,3)*Z1(ista,:)

       !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
       vz(ista,:) = X2(ista,:)
       vx(ista,:) = Y2(ista,:)
       vy(ista,:) = Z2(ista,:)

    enddo

  end subroutine rotate_comp_glob2mesh
!--------------------------------------------------------------------------------


!================================================================================
! Rotate data from mesh coordinates (vx,vy,vz) to earth coordinates (vn, ve, vz)
  subroutine local_mesh_coordinate(lat, lon, ele, x, y, z)

    real(kind=dp),  intent(in) :: lat, lon, ele
    real(kind=dp), intent(out) :: x, y, z

    real(kind=dp) :: X1, Y1, Z1, X2, Y2, Z2

    !* Get Cartesian coordinates
    call sph2cart(rt+ele, lat, lon, X1, Y1, Z1)

    !* From local mesh to global eath
    X2 = rotmat(1,1)*X1 + rotmat(1,2)*Y1 + rotmat(1,3)*Z1
    Y2 = rotmat(2,1)*X1 + rotmat(2,2)*Y1 + rotmat(2,3)*Z1
    Z2 = rotmat(3,1)*X1 + rotmat(3,2)*Y1 + rotmat(3,3)*Z1

    !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
    z = X2 - rt
    x = Y2
    y = Z2

  end subroutine local_mesh_coordinate
!--------------------------------------------------------------------------------


!================================================================================
! Rotate data from mesh coordinates (vx,vy,vz) to earth coordinates (vn, ve, vz)
  subroutine global_earth_coordinate(x, y, z, lat, lon, ele)

    real(kind=dp),  intent(in) :: x, y, z
    real(kind=dp), intent(out) :: lat, lon, ele

    real(kind=dp) :: X1, Y1, Z1, X2, Y2, Z2

    !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
    X2 = z + rt
    Y2 = x
    Z2 = y

    !* From local mesh to global eath
    X1 = rotmat_t(1,1)*X2 + rotmat_t(1,2)*Y2 + rotmat_t(1,3)*Z2
    Y1 = rotmat_t(2,1)*X2 + rotmat_t(2,2)*Y2 + rotmat_t(2,3)*Z2
    Z1 = rotmat_t(3,1)*X2 + rotmat_t(3,2)*Y2 + rotmat_t(3,3)*Z2

    !* Get spherical coordinates
    call cart2sph(X1, Y1, Z1, ele, lat, lon)
    ele = ele - rt

  end subroutine global_earth_coordinate
!--------------------------------------------------------------------------------

!================================================================================
! Spherical to Cartesian coordinates
  subroutine sph2cart(ele,lat,lon,x,y,z)

    real(kind=dp), intent(out) :: x, y, z
    real(kind=dp),  intent(in) :: lat, lon, ele

    real(kind=dp) :: la, lo

    la = deg2rad * lat
    lo = deg2rad * lon

    x = ele * cos(la) * cos(lo)
    Y = ele * cos(la) * sin(lo)
    Z = ele * sin(la)

  end subroutine sph2cart
!--------------------------------------------------------------------------------

!================================================================================
! Spherical to Cartesian coordinates
  subroutine cart2sph(x, y, z, ele, lat, lon)

    real(kind=dp),  intent(in) :: x, y, z
    real(kind=dp), intent(out) :: lat, lon, ele

    ele = sqrt(x*x + y*y + z*z)
    lat = rad2deg * asin(z/ele)
    lon = rad2deg * atan2(y,x)

  end subroutine cart2sph
!--------------------------------------------------------------------------------

!================================================================================
! Earth curvature (spherical)
  subroutine z_earth_curv(x,y,z)

    real(kind=dp),  intent(in) :: x, y
    real(kind=dp), intent(out) :: z

    z = sqrt(rt*rt - x*x - y*y)

  end subroutine z_earth_curv
!--------------------------------------------------------------------------------

!================================================================================
! Rotation of components
  subroutine rotate_ZNE_to_ZRT(vz,vn,ve,vz2,vr,vt,nrec,nt,bazi)

    integer(kind=si), intent(in)                :: nt, nrec
    real(kind=cp), dimension(nrec),   intent(in) :: bazi

    real(kind=cp), dimension(nrec,nt),  intent(in) :: vz,  vn, ve
    real(kind=cp), dimension(nrec,nt), intent(out) :: vz2, vr, vt

    real(kind=dp), dimension(nrec) :: baz

    integer(kind=si) :: it

    baz = deg2rad * bazi

    do it = 1, nt
       vr(:,it) = -ve(:,it) * sin(baz(:)) - vn(:,it) * cos(baz(:))
       vt(:,it) = -ve(:,it) * cos(baz(:)) + vn(:,it) * sin(baz(:))
       vz2(:,it) = vz(:,it)
    enddo

  end subroutine rotate_ZNE_to_ZRT

  subroutine rotate_ZRT_to_ZNE(vz2,vr,vt,vz,vn,ve,nrec,nt,bazi)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=cp), dimension(nrec),   intent(in) :: bazi

    real(kind=cp), dimension(nrec,nt),  intent(in) :: vz2, vr, vt
    real(kind=cp), dimension(nrec,nt), intent(out) :: vz,  vn, ve

    real(kind=dp), dimension(nrec) :: baz
    integer(kind=si) :: it

    baz = deg2rad * bazi

    do it = 1, nt
       ve(:,it) = -vr(:,it) * sin(baz(:)) - vt(:,it) * cos(baz(:))
       vn(:,it) = -vr(:,it) * cos(baz(:)) + vt(:,it) * sin(baz(:))
       vz(:,it) = vz2(:,it)
    enddo

  end subroutine rotate_ZRT_to_ZNE

  subroutine rotate_ZNE_to_LQT(vz,vn,ve,vl,vq,vt,nrec,nt,bazi,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=cp), dimension(nrec),   intent(in) :: bazi, inci

    real(kind=cp), dimension(nrec,nt),  intent(in) :: vz, vn, ve
    real(kind=cp), dimension(nrec,nt), intent(out) :: vl, vq, vt

    real(kind=dp), dimension(nrec) :: baz, inc
    integer(kind=si) :: it

    baz = deg2rad * bazi
    inc = deg2rad * inci

    do it = 1, nt
       vl(:,it) =  vz(:,it)*cos(inc(:)) - ve(:,it)*sin(baz(:))*sin(inc(:)) - vn(:,it)*cos(baz(:))*sin(inc(:))
       vq(:,it) =  vz(:,it)*sin(inc(:)) + ve(:,it)*sin(baz(:))*cos(inc(:)) + vn(:,it)*cos(baz(:))*cos(inc(:))
       vt(:,it) =              - ve(:,it)*cos(baz(:))          + vn(:,it)*sin(baz(:))
    enddo

  end subroutine rotate_ZNE_to_LQT

  subroutine rotate_LQT_to_ZNE(vl,vq,vt,vz,vn,ve,nrec,nt,bazi,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=cp),  dimension(nrec), intent(in)    :: bazi, inci

    real(kind=cp), dimension(nrec,nt), intent(out) :: vz, vn, ve
    real(kind=cp), dimension(nrec,nt),  intent(in) :: vl, vq, vt

    real(kind=dp),  dimension(nrec) :: baz, inc
    integer(kind=si) :: it

    baz = deg2rad * bazi
    inc = deg2rad * inci

    do it = 1, nt
       vz(:,it) =  vl(:,it)*cos(inc(:))          + vq(:,it)*sin(inc(:))
       ve(:,it) = -vl(:,it)*sin(baz(:))*sin(inc(:)) + vq(:,it)*sin(baz(:))*cos(inc(:)) - vt(:,it)*cos(baz(:))
       vn(:,it) = -vl(:,it)*cos(baz(:))*sin(inc(:)) + vq(:,it)*cos(baz(:))*cos(inc(:)) + vt(:,it)*sin(baz(:))
    enddo

  end subroutine rotate_LQT_to_ZNE

  subroutine rotate_ZRT_to_LQT(vz,vr,vt,vl,vq,vt2,nrec,nt,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=cp),  dimension(nrec), intent(in)    :: inci

    real(kind=cp), dimension(nrec,nt),  intent(in) :: vz, vr, vt
    real(kind=cp), dimension(nrec,nt), intent(out) :: vl, vq, vt2

    real(kind=dp), dimension(nrec) :: inc
    integer(kind=si) :: it

    inc = deg2rad * inci

    do it = 1, nt
       vl(:,it)  = vz(:,it) * cos(inc(:)) + vr(:,it) * sin(inc(:))
       vq(:,it) = vz(:,it) * sin(inc(:)) - vr(:,it) * cos(inc(:))
       vt2(:,it) = vt(:,it)
    enddo

  end subroutine rotate_ZRT_to_LQT

  subroutine rotate_LQT_to_ZRT(vl,vq,vt2,vz,vr,vt,nrec,nt,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=cp), dimension(nrec), intent(in)    :: inci

    real(kind=cp), dimension(nrec,nt), intent(out) :: vz, vr, vt
    real(kind=cp), dimension(nrec,nt),  intent(in) :: vl, vq, vt2

    real(kind=dp), dimension(nrec) :: inc
    integer(kind=si) :: it

    inc = deg2rad * inci

    do it = 1, nt
       vz(:,it) = vl(:,it) * cos(inc(:)) + vq(:,it) * sin(inc(:))
       vr(:,it) = vl(:,it) * sin(inc(:)) - vq(:,it) * cos(inc(:))
       vt(:,it) = vt2(:,it)
    enddo

  end subroutine rotate_LQT_to_ZRT
!--------------------------------------------------------------------------------

end module
