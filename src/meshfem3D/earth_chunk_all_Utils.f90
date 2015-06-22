!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine earth_chunk_ReadIasp91(vp,vs,rho,rb,n)

  implicit none

  integer i,j,n,iunit,nlay,nco(n),ifanis
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n)
  real fref,vph,vsh,qm,qk,eta

  character(len=80) text
  character(len=2) cnlay
  character(len=11) format_to_use

  do i=1,n
     !qm(i)=0.d0
        !qk(i)=0.d0
     rb(i)=0.d0
     !iflso(i)=0
     nco(i)=0
     do j=1,4
        rho(i,j)=0.d0
        vp(i,j)=0.d0
        !vph(i,j)=0.d0
        vs(i,j)=0.d0
        !vsh(i,j)=0.d0
        !eta(i,j)=0.d0
     enddo
  enddo
  iunit=26
  open(unit=iunit,file='iasp91',status='old')


1 read(iunit,'(a72)') text
  if (text(1:1)=='#') then
     goto 1
  endif
  backspace iunit


  read(iunit,'(i2)') nlay                ! Number of layers

  write(cnlay,'(i2)') nlay
  format_to_use='('//cnlay//'i2)'                 ! Number of polynomial
  read(iunit,format_to_use) (nco(i),i=1,nlay)     ! coefficients for each layer

  read(iunit,*) fref               ! reference frequency of Qs in Hertz
  read(iunit,*) ifanis             ! Transversal isotropic? 1=y, else=n
  read(iunit,'(1x/1x/)')


  do i = 1, nlay

     !read(iunit,*) rb(i-1),rho(i,1),vpv(i,1),vph(i,1),vsv(i,1),vsh(i,1),qm(i),qk(i),eta(i,1)
     read(iunit,*) rb(i),rho(i,1),vp(i,1),vph,vs(i,1),vsh,qm,qk,eta
     !write(*,*) i,rb(i),rho(i,1)
     do j = 2, nco(i)
        read(iunit,*) rho(i,j),vp(i,j),vph,vs(i,j),vsh,eta
        !write(*,*) i,j,rho(i,j)
     enddo
     read(iunit,'(1x)')
  enddo
  i = nlay+1
  read(iunit,*) rb(i)
  j = 1
  rho(i,j) =  rho(i-1,j)
  vp(i,j) = vp(i-1,j)
  vs(i,j) = vs(i-1,j)

  return

  end subroutine earth_chunk_ReadIasp91

!
!===========================================================================
!

  subroutine Read_dsm_model(model_file,vp,vs,rho,rb,n)

  implicit none

  integer i,n,iunit,nco(n)
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n),eta(4),vrmin,vrmax
  real vph(4),vsh(4),qm,qk
  integer nzone

  character(len=250) model_file

  rb    = 0.d0
  rho   = 0.d0
  vp    = 0.d0
  vs    = 0.d0
  nco   = 0
  iunit = 26

  open(unit=iunit,file=trim(model_file),status='old',action='read')

  read(iunit,*) nzone

  do i=1, nzone
     read(iunit,*) vrmin, vrmax, &
          rho(i,1), rho(i,2), rho(i,3), rho(i,4), &
          vp(i,1), vp(i,2), vp(i,3), vp(i,4), &
          vph(1), vph(2), vph(3), vph(4), &
          vs(i,1), vs(i,2), vs(i,3), vs(i,4), &
          vsh(1), vsh(2), vsh(3), vsh(4), &
          eta(1), eta(2), eta(3), eta(4),&
          qm, qk
          rb(i)=vrmin
  enddo

  i        = nzone+1
  rb(i)    = vrmax
  vp(i,:)  = vp(i-1,:)
  vs(i,:)  = vs(i-1,:)
  rho(i,:) = rho(i-1,:)

  close(iunit)

  end subroutine Read_dsm_model

!
!=======================================================================================================
!
! compute the Euler angles and the associated rotation matrix

  subroutine euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  implicit none

  !include "constants.h"

  double precision rotation_matrix(3,3)
  double precision CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  double precision alpha,beta,gamma
  double precision sina,cosa,sinb,cosb,sing,cosg

  double precision DEGREES_TO_RADIANS

  DEGREES_TO_RADIANS = 3.141592653589793d0/180.d0


! compute colatitude and longitude and convert to radians
  alpha = CENTER_LONGITUDE_IN_DEGREES * DEGREES_TO_RADIANS
  beta = (90.0d0 - CENTER_LATITUDE_IN_DEGREES) * DEGREES_TO_RADIANS
  gamma = GAMMA_ROTATION_AZIMUTH * DEGREES_TO_RADIANS

  sina = dsin(alpha)
  cosa = dcos(alpha)
  sinb = dsin(beta)
  cosb = dcos(beta)
  sing = dsin(gamma)
  cosg = dcos(gamma)

! define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
  rotation_matrix(1,3) = sinb*cosa
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
  rotation_matrix(2,3) = sinb*sina
  rotation_matrix(3,1) = -cosg*sinb
  rotation_matrix(3,2) = sing*sinb
  rotation_matrix(3,3) = cosb

  end subroutine euler_angles

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth,updown)

  implicit none

  integer NGLLX,NGLLY,NGLLZ,nel_depth,iz,Ndepth
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision profondeur
  integer current_layer(0:nel_depth-1),ilayer,k
  integer updown(NGLLZ) !! will be also used for VM coupling with AxiSEM

  updown(:) = 0
  if (ilayer ==  current_layer(iz)) then

    do k=2,NGLLZ
      profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
      write(27,*) profondeur/1000., ilayer-1,1
      Ndepth = Ndepth + 1
      updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM
    enddo

  else ! new layer

     k=1
     profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
     if (ilayer==0) then
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1

        updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM

     else
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,-1
        Ndepth=Ndepth+1

        updown(k) = -1 !! for new output mesh files and VM coupling with AxiSEM

     endif
     do k=2,NGLLZ ! on duplique le dernier point
        profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1

        updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM

     enddo


  endif

  end subroutine write_gllz_points


!=======================================================================================================
!
!=======================================================================================================
!
! Useless for the moment ==> maybe we will need it later
!
!!$! To have one and only file, who give the Spherical coordinate on ALL the GLL points
!!$! on the surface of the 3D chunk, for the new DSM coupling (light version using 2D chunk)
!!$!
!!$
!!$  subroutine cartesian_product_to_r_theta_phi_on_chunk_surface_GLL(MESH,deg2rad)
!!$
!!$  use constants, only: R_EARTH_KM
!!$
!!$  implicit none
!!$
!!$  character(len=10)  :: MESH
!!$  double precision   :: deg2rad
!!$  integer            :: np_r, np_xmin, np_xmax, np_ymin, np_ymax, np_zmin, recflag1, recflag2, i, j, np_surf, ios
!!$  double precision   :: rec_val, xmin_val1, xmin_val2, xmax_val1, xmax_val2, ymin_val1, ymin_val2
!!$  double precision   :: ymax_val1, ymax_val2, zmin_val1, zmin_val2, zmin_fix, x, y ,z, R, R_m, latrad, lgrad
!!$
!!$  open(unit=10,file=trim(MESH)//'recdepth',action='read',status='unknown',iostat=ios)
!!$  open(unit=11,file=trim(MESH)//'stxmin',action='read',status='unknown',iostat=ios)
!!$  open(unit=12,file=trim(MESH)//'stxmax',action='read',status='unknown',iostat=ios)
!!$  open(unit=13,file=trim(MESH)//'stymin',action='read',status='unknown',iostat=ios)
!!$  open(unit=14,file=trim(MESH)//'stymax',action='read',status='unknown',iostat=ios)
!!$  open(unit=15,file=trim(MESH)//'stzmin',action='read',status='unknown',iostat=ios)
!!$
!!$  open(unit=20,file=trim(MESH)//'chunk_surface_GLL_r_theta_phi.out',status='unknown',iostat=ios)
!!$
!!$  read(10,*) np_r
!!$
!!$  read(11,*) np_xmin
!!$  read(12,*) np_xmax
!!$  read(13,*) np_ymin
!!$  read(14,*) np_ymax
!!$  read(15,*) np_zmin
!!$
!!$  np_surf = np_r*(np_xmin + np_xmax + np_ymin + np_ymax) + np_zmin
!!$
!!$  write(20,*) np_surf
!!$
!!$  do i=1,np_r
!!$
!!$    rewind(11)
!!$    read(11,*)
!!$    rewind(12)
!!$    read(12,*)
!!$    rewind(13)
!!$    read(13,*)
!!$    rewind(14)
!!$    read(14,*)
!!$
!!$    read(10,*) rec_val, recflag1, recflag2
!!$
!!$    R = dabs(R_EARTH_KM - rec_val) !! kM

!!$    do j=1,np_xmin
!!$
!!$      read(11,*) xmin_val1, xmin_val2
!!$      write(20,*) R, xmin_val1, xmin_val2
!!$
!!$ enddo
!!$
!!$    do j=1,np_xmax
!!$
!!$      read(12,*) xmax_val1, xmax_val2
!!$      write(20,*) R, xmax_val1, xmax_val2
!!$
!!$    enddo
!!$
!!$    do j=1,np_ymin
!!$
!!$      read(13,*) ymin_val1, ymin_val2
!!$      write(20,*) R, ymin_val1, ymin_val2
!!$
!!$    enddo
!!$
!!$    do j=1,np_ymax
!!$
!!$      read(14,*) ymax_val1, ymax_val2
!!$      write(20,*) R, ymax_val1, ymax_val2
!!$
!!$ enddo
!!$
!!$ if (i == np_r) zmin_fix = rec_val !! maximal depth
!!$
!!$  enddo
!!$
!!$  rewind(15)v
!!$  read(15,*)
!!$
!!$  R = dabs(R_EARTH_KM - zmin_fix) !! kM
!!$
!!$  do j=1,np_zmin
!!$
!!$    read(15,*) zmin_val1, zmin_val2
!!$    write(20,*) R, zmin_val1, zmin_val2
!!$
!!$  enddo
!!$
!!$  close(10)
!!$  close(11)
!!$  close(12)
!!$  close(13)
!!$  close(14)
!!$  close(15)
!!$  close(20)
!!$
!!$  end subroutine cartesian_product_to_r_theta_phi_on_chunk_surface_GLL
!
!=======================================================================================================
!
!=======================================================================================================
!
!! Used (among other) for VM coupling with AxiSEM

  subroutine write_all_chunk_surface_GLL_in_spherical_and_cartesian_coords(xstore,ystore,zstore, &
                                                      deg2rad,ilayer,iboun,ispec,nspec,longitud, &
                                                      latitud,radius,rotation_matrix,updown)

  use constants, only: NGLLX, NGLLY, NGLLZ, EXTERNAL_CODE_IS_DSM, EXTERNAL_CODE_IS_AXISEM

  use shared_parameters, only: EXTERNAL_CODE_TYPE

  implicit none

  integer ispec, nspec, ilayer
  integer i, j, k, imin, imax, jmin, jmax, kmin, kmax
  integer updown(NGLLZ)
  logical :: iboun(6,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: longitud, latitud, radius
  double precision rotation_matrix(3,3)
  double precision deg2rad

!
!---- CF 'earth_chunk_HEX8_Mesher' and 'earth_chunk_HEX27_Mesher' to see the name of files whose units are 91 and 92
!

1000 format(3f30.10)

!-- all gll in geographical coordinates

  call cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)

!
!-- xmin ----
!

  if (iboun(1,ispec)) then

    imin = 1
    imax = 1
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          ! CF 'earth_chunk_HEX8_Mesher' and 'earth_chunk_HEX27_Mesher' to see files whose units are 91 and 92

          if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,1,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- xmax ----
!

  if (iboun(2,ispec)) then

    imin = NGLLX
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,2,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- ymin ----
!

  if (iboun(3,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = 1
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,3,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- ymax ----
!

  if (iboun(4,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = NGLLY
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,4,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- zmin ----
!

  if (iboun(5,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = 1

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,5,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

  end subroutine write_all_chunk_surface_GLL_in_spherical_and_cartesian_coords

!
!=======================================================================================================
!

  subroutine cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM

  implicit none

  integer i, j, igll, jgll, kgll

  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: longitud, latitud, radius
  double precision rotation_matrix(3,3)
  double precision vector_ori(3), vector_rotated(3)
  double precision rayon, x, y, z, long, lati, deg2rad

!
!----
!

  do kgll=1,NGLLZ
    do jgll=1,NGLLY
      do igll=1,NGLLX

        vector_ori(1) = xstore(igll,jgll,kgll)
        vector_ori(2) = ystore(igll,jgll,kgll)
        vector_ori(3) = zstore(igll,jgll,kgll)
        !write(*,*)  vector_ori
        !write(*,*) rotation_matrix

        do i = 1,NDIM

          vector_rotated(i) = 0.d0

          do j = 1,NDIM

            !write(*,*)  vector_rotated(i),rotation_matrix(i,j),vector_ori(j)
            vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)

          enddo
        enddo

        x     = vector_rotated(1)
        y     = vector_rotated(2)
        z     = vector_rotated(3)
        rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

        long  = datan2(y,x)
        lati  = dasin(z/rayon)

        longitud(igll,jgll,kgll) = long/deg2rad
        latitud(igll,jgll,kgll)  = lati/deg2rad
        radius(igll,jgll,kgll)   = rayon/1000.d0

      enddo
    enddo
  enddo

  end subroutine cartesian2spheric


!=======================================================================================================
!
!=======================================================================================================

  subroutine write_recdepth_dsm(Ndepth,R_EARTH,MESH)

  implicit none

  integer Ndepth,i
  double precision R_EARTH,prof
  double precision, allocatable :: z(:)
  integer, allocatable :: zindex(:),ziflag(:)
  integer ilayer,flag
  character(len=10) MESH

  open(27,file=trim(MESH)//'.recdepth')
  allocate(zindex(Ndepth),ziflag(Ndepth))
  allocate(z(Ndepth))

  do i=1,Ndepth
     read(27,*) prof,ilayer,flag
     z(Ndepth-i+1)=R_EARTH/1000.d0-prof
     zindex(Ndepth-i+1)=ilayer
     ziflag(Ndepth-i+1)=flag
  enddo
  close(27)

  open(27,file=trim(MESH)//'recdepth')
  write(27,*) Ndepth
  i=1
  write(27,*) z(i),zindex(i),ziflag(i)
  do i=2,Ndepth-1
     if (ziflag(i-1) == -1 ) then
        write(27,*) z(i),zindex(i),-1
     else
         write(27,*) z(i),zindex(i),1
     endif
  enddo
  i=Ndepth
  write(27,*) z(i),zindex(i),ziflag(i)
  close(27)

end subroutine write_recdepth_dsm

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(1,jgll,NGLLZ)
     vector_ori(2)=ystore(1,jgll,NGLLZ)
     vector_ori(3)=zstore(1,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique a geographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(28,*) long/deg2rad,lati/deg2rad !,rayon/1000
      !write(38,'()') 1,(NGLLY-1)*jy_elm+jgll
       write(49,*)
       write(49,*)     vector_ori(:)
       write(49,*)     vector_rotated(:)

  enddo

  end subroutine write_stxmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(NGLLX,jgll,NGLLZ)
     vector_ori(2)=ystore(NGLLX,jgll,NGLLZ)
     vector_ori(3) =zstore(NGLLX,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique a geographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(29,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stxmax

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

   if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,1,NGLLZ)
     vector_ori(2)=ystore(jgll,1,NGLLZ)
     vector_ori(3) =zstore(jgll,1,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique a geographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(30,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,NGLLY,NGLLZ)
     vector_ori(2)=ystore(jgll,NGLLY,NGLLZ)
     vector_ori(3) =zstore(jgll,NGLLY,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      ! passage de geocentique a geographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      write(31,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymax

!=======================================================================================================
!
!=======================================================================================================

  subroutine store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,&
             lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,igll,jgll,i,j
  integer ilon,ilat,iglob,jglob,nlat_dsm,nlon_dsm
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  double precision lon_zmin(nlon_dsm,nlat_dsm),lat_zmin(nlon_dsm,nlat_dsm)


  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY
     do igll=1,NGLLX
        vector_ori(1)=xstore(igll,jgll,1)
        vector_ori(2)=ystore(igll,jgll,1)
        vector_ori(3) =zstore(igll,jgll,1)

        do i = 1,NDIM
           vector_rotated(i) = 0.d0
           do j = 1,NDIM
              vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
           enddo
        enddo
        x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
        rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

        long=atan2(y,x)
        lati=asin(z/rayon)

      ! passage de geocentique a geographique
      !!theta = PI/2.D0 - lati
      ! convert the geocentric colatitude to a geographic colatitude
      !!colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
      !!lati =   PI/2.0d0 - colat

      !write(28,*) xstore(1,jgll,NGLLZ), ystore(1,jgll,NGLLZ), zstore(1,jgll,NGLLZ)!x,y !long/deg2rad,lati/deg2rad
      !write(31,*) long/deg2rad,lati/deg2rad !,rayon/1000
        iglob=(ilon)*(NGLLX-1)+igll
        jglob=(ilat)*(NGLLY-1)+jgll
        lon_zmin(iglob,jglob)= long/deg2rad
        lat_zmin(iglob,jglob)= lati/deg2rad
        !write(32,'(3f20.10)') xstore(igll,jgll,1)/1000.d0, ystore(igll,jgll,1)/1000.d0,zstore(igll,jgll,1)/1000.d0
        !write(32,*) xstore(igll,jgll,NGLLZ), ystore(igll,igll,NGLLZ),zstore(igll,jgll,NGLLZ)
     enddo
  enddo

  end subroutine store_zmin_points

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stzmin(x,y,nx,ny,MESH)

  implicit none

  integer i,j,nx,ny
  double precision x(nx,ny),y(nx,ny)
  character(len=10) MESH

  open(27,file=trim(MESH)//'stzmin')
  write(27,*) nx*ny
  do j=1,ny
     do i=1,nx
        write(27,*) x(i,j),y(i,j)
     enddo
  enddo
  close(27)

  end subroutine write_stzmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_Igm_file(iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il)

  implicit none

  integer iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il
  integer i,j
  do j=1,NGLL2
     do i=1,NGLL1
        write(iunit,*) i,j,ispec2D,(NGLL1-1)*ie+i,(NGLL2-1)*je+j+js,il
     enddo
  enddo

  end subroutine write_Igm_file

!=======================================================================================================
!
!=======================================================================================================

  subroutine  compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  implicit none

  double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
  double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

  ! je met le chunk en 0,0
  axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
  ! rotation de l'azimuth du chunk
  axe_rotation(1)=1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)
  ! on met le chunk a la bonne latitude
  axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)
  ! on met le chunk a la bonne longitude
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
  call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)
  ! rotation resultante
  call compose4matrix(rotation_matrix,R00,R0,R1,R2)

  end subroutine compute_rotation_matrix

!=======================================================================================================
!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
! Vadim Monteiller Mars 2013
!
!-------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (en degres)
! cette matrice est en complexe
!
!=======================================================================================================
!
  subroutine rotation_matrix_axe(R,axe,theta)

  implicit none

  double precision axe(3),theta,pi,deg2rad
  double precision R(3,3)
  double precision c,s,ux,uy,uz,norme_axe

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! on normalise l'axe
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! composantes de l'axe
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! on calcule le cos et sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! matrice de rotation complexe
  R(1,1)=(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=(ux*uy*(1.d0-c)+uy*s)

  R(2,1)=(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=(uy**2+(1.d0-uy**2)*c)
  R(2,3)=(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=(uz**2+(1.d0-uz**2)*c)

  write(49,*) ' MATRICE ROTATION '
  write(49,*) R(1,:)
  write(49,*) R(2,:)
  write(49,*) R(3,:)
  write(49,*)

  end subroutine rotation_matrix_axe

!=======================================================================================================
!
! R=R2*R1*R0
!
!=======================================================================================================

  subroutine compose4matrix(R,R00,R0,R1,R2)

  implicit none

  double precision R(3,3),R0(3,3),R1(3,3),R2(3,3),R00(3,3),Rtmp(3,3)
  integer i,j,k


  R(:,:)=0.d0
  ! multiplication R=R0*R00
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R0(i,k)*R00(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R1*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R1(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R2*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R2(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  write(49,*) ' MATRICE ROTATION COMPLETE '
  write(49,*) R(1,:)
  write(49,*) R(2,:)
  write(49,*) R(3,:)
  write(49,*)

  end subroutine compose4matrix

!------------------------------------------------------------------------------
! rotation pour passer d'un repere local a un autre
!===========================================================================!
!

subroutine Lyfnd(r,rb,n,i)
  implicit none
  integer i,n
  double precision r,rb(n)
  i=1
  do while (r > rb(i) )
     i = i + 1
  enddo
  i = i - 1
  return
end subroutine Lyfnd

function IsNewLayer(x,r,n)
  implicit none
  integer IsNewLayer,n,i
  double precision x,r(n)
  IsNewLayer = 0
  ! ce test fonctionne que si les mailles sont suffisament petites !! ATTENTION
  do i = 1, n-1
     !write(*,*) x,r(i),x-r(i)
     if (abs(x-r(i))<1.d-10) then
        !write(*,*) 'its layer'
        IsNewLayer = 1
        return
     endif
  enddo
end function IsNewLayer


subroutine StorePoint(z,k,zc)
  implicit none
  integer k
  double precision z(*),zc
  if (k==0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k)==zc) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePoint

subroutine StorePointZ(z,k,zc,NoInter)
  implicit none
  integer k
  double precision z(*),zc
  logical NoInter
  if (k==0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k)==zc.and.NoInter) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePointZ

!!$subroutine FindNum(i,X,X0,n)
!!$  implicit none
!!$  integer i,n
!!$  real X(n),X0
!!$  do i=1,n
!!$     if (X0==X(n)) return
!!$  enddo
!!$  write(*,*) ' warrnig FindnNum ',X0
!!$  stop
!!$end subroutine FindNum
 subroutine CalGridProf(ProfForGemini,Niveau_elm,zlayer,nlayer,NEX_GAMMA,Z_DEPTH_BLOCK)

  implicit none
  integer NEX_GAMMA,nlayer,nbbloc(100000),Niveau_elm(0:NEX_GAMMA-1)
  double precision ProfForGemini(0:NEX_GAMMA-1,3),zlayer(nlayer)
  double precision Z_DEPTH_BLOCK,zpoint(100000),zz(100000)
  double precision epsillon
  integer nb, n, i,j,k,ilayer,ilay,nd,niveau
  double precision p, pas, longeur
  logical test

  epsillon=1d-3
   nbbloc(:)=0
   ! point de depart
   zpoint(1)=zlayer(nlayer) - Z_DEPTH_BLOCK
   write(*,*) zlayer(nlayer) ,  Z_DEPTH_BLOCK
   !! niveau de depart
   call FindLayer_for_earth_chunk_mesh(ilayer,zlayer,zpoint(1),nlayer)
   write(*,*) '              INITIALISATION calcul du niveau de depart : '
   write(*,*)
   write(*,*) 'zlayer : ', zlayer
   write(*,*) 'premier point : '   , zpoint(1),ilayer
    write(*,*)

  !!


  !! on compte le nombre d'elements par niveau
  i = 1
  k = ilayer - 1
  nb = 0
  do while (zpoint(i)<zlayer(nlayer))
    i = i + 1
    k = k + 1
    zpoint(i) = zlayer(k)
  enddo

  nb = i
  nd = i-1
  longeur = zlayer(nlayer) - zpoint(1)


  do i=1,nb-1

     pas = zpoint(i+1) - zpoint(i)
     p = NEX_GAMMA * pas / longeur

    if (p < 0.8d0) then
        n = 1
    else
        n = max(int(p),2)
    endif
    !write(*,*) 'n :',n

    nbbloc(i)=n
  enddo






  do j=1,nb-1
    write(*,*) j,nbbloc(j)
  enddo


  !! on elimine les blocs en trop
   write(*,*) 'SUM ',sum(nbbloc)

   nb = sum(nbbloc)


   do while (nb > NEX_GAMMA)

      k  =  1
      test = .true.

    do  while (test)

         j =  maxval(nbbloc)
         ! on cherche l'indice du max

         if (j == nbbloc(k)) then
            nbbloc(k ) = nbbloc(k) -1
            test = .false.
         endif

         k = k + 1

      enddo

      nb = sum(nbbloc)
      write(*,*) 'nb, ',nb,NEX_GAMMA
   enddo

  !!
  longeur = zlayer(nlayer) - zpoint(1)
  k=1
  zz(k)=zpoint(1)
  !zpoint(nb+1)=zlayer(nlayer)
  do i=1,nd
     !write(*,*) i,nbbloc(i)
     pas = (zpoint(i+1) - zpoint(i)) / nbbloc(i)
     write(*,*) i,nbbloc(i),pas
     do while (zz(k) < zpoint(i+1) - epsillon)
        k = k + 1
        zz(k) = zz(k-1) + pas
        write(*,*) zz(k), zpoint(i+1)
     enddo
  enddo




   do ilay=1,NEX_GAMMA

      ProfForGemini(ilay-1,1)  =  zz(ilay)
      ProfForGemini(ilay-1,2)  =  zz(ilay+1)
      ProfForGemini(ilay-1,3)  = 0.5d0 * (zz(ilay) + zz(ilay+1))

      call FindLayer_for_earth_chunk_mesh(niveau,zlayer, ProfForGemini(ilay-1,3),nlayer)
      Niveau_elm(ilay-1)=niveau
      write(*,'(i5,2f15.3,i10)') ilay,zz(ilay),zz(ilay+1),niveau
   enddo

   do ilay=0,NEX_GAMMA-1
     !write(*,*) ProfForGemini(ilay,1), ProfForGemini(ilay,2), ProfForGemini(ilay,3)
   enddo


 end subroutine CalGridProf

 subroutine  FindLayer_for_earth_chunk_mesh(i,z,r,n)
   implicit none
   integer i,n
   double precision z(n),r

   if (r>z(n) .or. r<z(1)) then
    write(*,*) 'STOP :: point ouside grid'
    stop
   endif
   i = 1
   do while (r > z(i))
     i = i + 1
   enddo


 end subroutine FindLayer_for_earth_chunk_mesh





!=====================================================================

  subroutine get_global1(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD,UTM_X_MIN,UTM_X_MAX)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  integer NGNOD

  integer npointot
  integer nspec,nglob
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  double precision UTM_X_MIN,UTM_X_MAX

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision SMALLVALTOL

! define geometrical tolerance based upon typical size of the model
  SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) ' SMALLVALTOL  ',SMALLVALTOL
! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers (!! VM changed NGLLCUBE (as in Specfem3D Basin Version 1.1) to NGNOD !!)
  do ispec=1,nspec
    ieoff = NGNOD * (ispec - 1)
    do ilocnum = 1,NGNOD
      loc(ilocnum + ieoff) = ilocnum + ieoff
    enddo
  enddo

  ifseg(:)  = .false.
  nseg      = 1
  ifseg(1)  = .true.
  ninseg(1) = npointot

  do j=1,3 !,NDIM

! sort within each segment
  ioff = 1

  do iseg=1,nseg

    if (j == 1) then
      call rank(xp(ioff),ind,ninseg(iseg))
    else if (j == 2) then
      call rank(yp(ioff),ind,ninseg(iseg))
    else
      call rank(zp(ioff),ind,ninseg(iseg))
    endif

    call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))

    ioff = ioff + ninseg(iseg)

  enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
  if (j == 1) then

    do i=2,npointot
      if (dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  else if (j == 2) then

    do i=2,npointot
      if (dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  else

    do i=2,npointot
      if (dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  endif

! count up number of different segments
  nseg = 0

  do i=1,npointot
    if (ifseg(i)) then
      nseg = nseg + 1
      ninseg(nseg) = 1
    else
      ninseg(nseg) = ninseg(nseg) + 1
    endif
  enddo
  enddo

! assign global node numbers (now sorted lexicographically)
  ig = 0

  do i=1,npointot
    if (ifseg(i)) ig = ig + 1

    iglob(loc(i)) = ig
  enddo

  nglob = ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global1


! -----------------------------------

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   endif
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF (A(IND(j))<A(IND(j+1))) j=j+1
      endif
      IF (q<A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   goto 200
   endif
   IND(I)=INDX
  goto 100
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all


!!$!=====================================================================
!!$
!!$! 3D shape functions for 8-node element
!!$
!!$  subroutine get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ)
!!$
!!$  implicit none
!!$
!!$  !include "constants.h"
!!$
!!$  integer myrank,NGNOD,NGLLX,NGLLY,NGLLZ
!!$  integer, parameter :: NDIM=3
!!$
!!$! Gauss-Lobatto-Legendre points of integration
!!$  double precision xigll(NGLLX)
!!$  double precision yigll(NGLLY)
!!$  double precision zigll(NGLLZ)
!!$
!!$! 3D shape functions and their derivatives
!!$  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
!!$  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
!!$
!!$  integer i,j,k,ia
!!$
!!$! location of the nodes of the 3D quadrilateral elements
!!$  double precision xi,eta,gamma
!!$  double precision ra1,ra2,rb1,rb2,rc1,rc2
!!$
!!$! for checking the 3D shape functions
!!$  double precision sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma
!!$
!!$  double precision, parameter :: ONE_EIGHTH = 0.125d0, ZERO = 0.d0, one=1.d0,TINYVAL = 1.d-9
!!$
!!$
!!$! check that the parameter file is correct
!!$  !myrank=0
!!$  if (NGNOD /= 8) call exit_MPI(myrank,'elements should have 8 control nodes')
!!$
!!$! ***
!!$! *** create 3D shape functions and jacobian
!!$! ***
!!$
!!$!--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
!!$
!!$  do i=1,NGLLX
!!$  do j=1,NGLLY
!!$  do k=1,NGLLZ
!!$
!!$  xi = xigll(i)
!!$  eta = yigll(j)
!!$  gamma = zigll(k)
!!$
!!$  ra1 = one + xi
!!$  ra2 = one - xi
!!$
!!$  rb1 = one + eta
!!$  rb2 = one - eta
!!$
!!$  rc1 = one + gamma
!!$  rc2 = one - gamma
!!$
!!$  shape3D(1,i,j,k) = ONE_EIGHTH*ra2*rb2*rc2
!!$  shape3D(2,i,j,k) = ONE_EIGHTH*ra1*rb2*rc2
!!$  shape3D(3,i,j,k) = ONE_EIGHTH*ra1*rb1*rc2
!!$  shape3D(4,i,j,k) = ONE_EIGHTH*ra2*rb1*rc2
!!$  shape3D(5,i,j,k) = ONE_EIGHTH*ra2*rb2*rc1
!!$  shape3D(6,i,j,k) = ONE_EIGHTH*ra1*rb2*rc1
!!$  shape3D(7,i,j,k) = ONE_EIGHTH*ra1*rb1*rc1
!!$  shape3D(8,i,j,k) = ONE_EIGHTH*ra2*rb1*rc1
!!$
!!$  dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
!!$  dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
!!$  dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
!!$  dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
!!$  dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
!!$  dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
!!$  dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
!!$  dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1
!!$
!!$  dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
!!$  dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
!!$  dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
!!$  dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
!!$  dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
!!$  dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
!!$  dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
!!$  dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1
!!$
!!$  dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
!!$  dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
!!$  dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
!!$  dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
!!$  dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
!!$  dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
!!$  dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
!!$  dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1
!!$
!!$  enddo
!!$  enddo
!!$  enddo
!!$
!!$!--- check the shape functions and their derivatives
!!$
!!$  do i=1,NGLLX
!!$  do j=1,NGLLY
!!$  do k=1,NGLLZ
!!$
!!$  sumshape = ZERO
!!$  sumdershapexi = ZERO
!!$  sumdershapeeta = ZERO
!!$  sumdershapegamma = ZERO
!!$
!!$  do ia=1,NGNOD
!!$    sumshape = sumshape + shape3D(ia,i,j,k)
!!$    sumdershapexi = sumdershapexi + dershape3D(1,ia,i,j,k)
!!$    sumdershapeeta = sumdershapeeta + dershape3D(2,ia,i,j,k)
!!$    sumdershapegamma = sumdershapegamma + dershape3D(3,ia,i,j,k)
!!$  enddo
!!$
!!$! sum of shape functions should be one
!!$! sum of derivative of shape functions should be zero
!!$  if (abs(sumshape-one) >  TINYVAL) call exit_MPI(myrank,'error shape functions')
!!$  if (abs(sumdershapexi) >  TINYVAL) call exit_MPI(myrank,'error derivative xi shape functions')
!!$  if (abs(sumdershapeeta) >  TINYVAL) call exit_MPI(myrank,'error derivative eta shape functions')
!!$  if (abs(sumdershapegamma) >  TINYVAL) call exit_MPI(myrank,'error derivative gamma shape functions')
!!$
!!$  enddo
!!$  enddo
!!$  enddo
!!$
!!$  end subroutine get_shape3D !

!
!!! subroutine exit_MPI(myrank,error_msg)

!!!  implicit none
!!!  integer myrank
!!!  character (len=*) error_msg
!!! write(*,*) error_msg
!!!  stop
!!!end subroutine exit_MPI

!
subroutine calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
  implicit none
  integer NGNOD,NGLLX,NGLLY,NGLLZ
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),xmesh,ymesh,zmesh
  double precision, parameter :: ZERO = 0.d0
  integer ia,i,j,k

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           ! compute mesh coordinates
           xmesh = ZERO
           ymesh = ZERO
           zmesh = ZERO
           do ia=1,NGNOD
              xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
              ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
              zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
           enddo
           xstore(i,j,k) = xmesh
           ystore(i,j,k) = ymesh
           zstore(i,j,k) = zmesh
        enddo
     enddo
  enddo
end subroutine calc_gll_points
