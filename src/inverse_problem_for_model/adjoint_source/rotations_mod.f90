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
    real(kind=dp),    intent(in) :: stalat, stalon

    real(kind=dp), dimension(nsta,nt),  intent(in) :: vx, vy, vz
    real(kind=dp), dimension(nsta,nt), intent(out) :: vn, ve, vz2

    real(kind=dp), dimension(nsta,nt) :: X1, Y1, Z1
    real(kind=dp), dimension(nsta,nt) :: X2, Y2, Z2

    real(kind=dp) :: lat, lon

    !* Pass in radian
    lat = deg2rad * stalat
    lon = deg2rad * stalon

    !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
    X2 = vz
    Y2 = vx
    Z2 = vy

    !* From local mesh to global eath
    X1 = rotmat_t(1,1)*X2 + rotmat_t(1,2)*Y2 + rotmat_t(1,3)*Z2
    Y1 = rotmat_t(2,1)*X2 + rotmat_t(2,2)*Y2 + rotmat_t(2,3)*Z2
    Z1 = rotmat_t(3,1)*X2 + rotmat_t(3,2)*Y2 + rotmat_t(3,3)*Z2

    !* Define rotation matrix
    rotc(1,1) =  cos(lat) * cos(lon)
    rotc(1,2) =  cos(lat) * sin(lon)
    rotc(1,3) =  sin(lat)
    rotc(2,1) = -sin(lat) * cos(lon)
    rotc(2,2) = -sin(lat) * sin(lon)
    rotc(2,3) =  cos(lat)
    rotc(3,1) = -sin(lon)
    rotc(3,2) =  cos(lon)
    rotc(3,3) =  0._dp

    !* Data in geographic coordinate at real position
    vz2 = rotc(1,1)*X1 + rotc(1,2)*Y1 + rotc(1,3)*Z1
    vn  = rotc(2,1)*X1 + rotc(2,2)*Y1 + rotc(2,3)*Z1
    ve  = rotc(3,1)*X1 + rotc(3,2)*Y1 + rotc(3,3)*Z1

  end subroutine rotate_comp_mesh2glob
!--------------------------------------------------------------------------------


!================================================================================
! Rotate data from mesh coordinates (vx,vy,vz) to earth coordinates (vn, ve, vz)
  subroutine rotate_comp_glob2mesh(vz2, vn, ve, stalat, stalon, nt, nsta, vx, vy, vz)

    integer(kind=si), intent(in) :: nt, nsta
    real(kind=dp),    intent(in) :: stalat, stalon

    real(kind=dp), dimension(nsta,nt), intent(out) :: vx, vy, vz
    real(kind=dp), dimension(nsta,nt),  intent(in) :: vn, ve, vz2

    real(kind=dp), dimension(nsta,nt) :: X1, Y1, Z1
    real(kind=dp), dimension(nsta,nt) :: X2, Y2, Z2

    real(kind=dp) :: lat, lon

    !* Pass in radian
    lat = deg2rad * stalat
    lon = deg2rad * stalon

    !* Define rotation matrix
    rotc(1,1) =  cos(lat) * cos(lon)
    rotc(1,2) =  cos(lat) * sin(lon)
    rotc(1,3) =  sin(lat)
    rotc(2,1) = -sin(lat) * cos(lon)
    rotc(2,2) = -sin(lat) * sin(lon)
    rotc(2,3) =  cos(lat)
    rotc(3,1) = -sin(lon)
    rotc(3,2) =  cos(lon)
    rotc(3,3) =  0._dp
    rotc_t = transpose(rotc)

    !* Data in geographic coordinate at real position
    X1 = rotc_t(1,1)*vz2 + rotc_t(1,2)*vn + rotc_t(1,3)*ve
    Y1 = rotc_t(2,1)*vz2 + rotc_t(2,2)*vn + rotc_t(2,3)*ve
    Z1 = rotc_t(3,1)*vz2 + rotc_t(3,2)*vn + rotc_t(3,3)*ve

    !* From local mesh to global eath
    X2 = rotmat(1,1)*X1 + rotmat(1,2)*Y1 + rotmat(1,3)*Z1
    Y2 = rotmat(2,1)*X1 + rotmat(2,2)*Y1 + rotmat(2,3)*Z1
    Z2 = rotmat(3,1)*X1 + rotmat(3,2)*Y1 + rotmat(3,3)*Z1

    !* Equivalence of local and global Cartesian coordinates at (lat=0, lon=0)
    vz = X2
    vx = Y2
    vy = Z2

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

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp),    intent(in) :: bazi
        
    real(kind=dp), dimension(nrec,nt),  intent(in) :: vz,  vn, ve
    real(kind=dp), dimension(nrec,nt), intent(out) :: vz2, vr, vt

    real(kind=dp) :: baz

    baz = deg2rad * bazi
    vr = -ve * sin(baz) - vn * cos(baz)
    vt = -ve * cos(baz) + vn * sin(baz)
    vz2 = vz

  end subroutine rotate_ZNE_to_ZRT

  subroutine rotate_ZRT_to_ZNE(vz2,vr,vt,vz,vn,ve,nrec,nt,bazi)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp),    intent(in) :: bazi
    
    real(kind=dp), dimension(nrec,nt),  intent(in) :: vz2, vr, vt
    real(kind=dp), dimension(nrec,nt), intent(out) :: vz,  vn, ve

    real(kind=dp) :: baz

    baz = deg2rad * bazi
    ve = -vr * sin(baz) - vt * cos(baz)
    vn = -vr * cos(baz) + vt * sin(baz)
    vz = vz2

  end subroutine rotate_ZRT_to_ZNE

  subroutine rotate_ZNE_to_LQT(vz,vn,ve,vl,vq,vt,nrec,nt,bazi,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp),    intent(in) :: bazi, inci

    real(kind=dp), dimension(nrec,nt),  intent(in) :: vz, vn, ve
    real(kind=dp), dimension(nrec,nt), intent(out) :: vl, vq, vt

    real(kind=dp) :: baz, inc

    baz = deg2rad * bazi
    inc = deg2rad * inci

    vl =  vz*cos(inc) - ve*sin(baz)*sin(inc) - vn*cos(baz)*sin(inc)
    vq =  vz*sin(inc) + ve*sin(baz)*cos(inc) + vn*cos(baz)*cos(inc)
    vt =              - ve*cos(baz)          + vn*sin(baz)

  end subroutine rotate_ZNE_to_LQT

  subroutine rotate_LQT_to_ZNE(vl,vq,vt,vz,vn,ve,nrec,nt,bazi,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp), intent(in)    :: bazi, inci

    real(kind=dp), dimension(nrec,nt), intent(out) :: vz, vn, ve
    real(kind=dp), dimension(nrec,nt),  intent(in) :: vl, vq, vt

    real(kind=dp) :: baz, inc

    baz = deg2rad * bazi
    inc = deg2rad * inci

    vz =  vl*cos(inc)          + vq*sin(inc)
    ve = -vl*sin(baz)*sin(inc) + vq*sin(baz)*cos(inc) - vt*cos(baz)
    vn = -vl*cos(baz)*sin(inc) + vq*cos(baz)*cos(inc) + vt*sin(baz)

  end subroutine rotate_LQT_to_ZNE

  subroutine rotate_ZRT_to_LQT(vz,vr,vt,vl,vq,vt2,nrec,nt,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp), intent(in)    :: inci
    
    real(kind=dp), dimension(nrec,nt),  intent(in) :: vz, vr, vt
    real(kind=dp), dimension(nrec,nt), intent(out) :: vl, vq, vt2
    
    real(kind=dp) :: inc

    inc = deg2rad * inci

    vl  = vz * cos(inc) + vr * sin(inc)
    vq  = vz * sin(inc) - vr * cos(inc)
    vt2 = vt

  end subroutine rotate_ZRT_to_LQT

  subroutine rotate_LQT_to_ZRT(vl,vq,vt2,vz,vr,vt,nrec,nt,inci)

    integer(kind=si), intent(in) :: nt, nrec
    real(kind=dp), intent(in)    :: inci
        
    real(kind=dp), dimension(nrec,nt), intent(out) :: vz, vr, vt
    real(kind=dp), dimension(nrec,nt),  intent(in) :: vl, vq, vt2

    real(kind=dp) :: inc

    inc = deg2rad * inci

    vz = vl * cos(inc) + vq * sin(inc)
    vr = vl * sin(inc) - vq * cos(inc)
    vt = vt2
    
  end subroutine rotate_LQT_to_ZRT
!--------------------------------------------------------------------------------

end module
