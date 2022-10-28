!
!            ROUTINES FOR 3D ROTATIONS
!
!       Vadim Monteiller Mars 2013
!
!-
module cart2geogr_mod

  ! chnuk parameters
  double precision ZREF
  double precision LON_CENTER,LAT_CENTER,AZI_CHUNK
  double precision rotation_matrix(3,3)

  integer NDIM
  parameter(NDIM=3)

  double precision pi
  parameter(pi=3.14159265359d0)

  double precision REARTH
  parameter(REARTH=6371000.d0)

contains

  subroutine read_chunk_parameters()

    integer i,n,m
    character(len=5) Cdummy
    open(10,file='../MESH/model_1D.in')
    read(10,*) n,m
    do i=1,n*m
       read(10,'(a)') Cdummy
    enddo
    read(10,*)  ZREF
    read(10,*) LON_CENTER,LAT_CENTER,AZI_CHUNK

  end subroutine read_chunk_parameters

  subroutine  cart2geogr(x,y,z,xlon,ylat,depth,ksi,eta)
    implicit none
    !include "constants.h"

    double precision x,y,z,xlon,ylat,depth
    double precision ksi,eta
    double precision vector_ori(3),vector_rotated(3),prof
    integer i,j


1000 format(3(f12.4,1x),5x,2(f15.5,1x),5x,4(f15.5,1x))
    vector_ori(1)=x
    vector_ori(2)=y
    vector_ori(3)=z+ZREF

    prof=dsqrt(vector_ori(1)**2 + vector_ori(2)**2 + vector_ori(3)**2)
    ksi = atan2( vector_ori(1),vector_ori(3))*180.d0/pi
    eta = atan2(vector_ori(2),vector_ori(3))*180.d0/pi

    do i = 1,NDIM
       vector_rotated(i) = 0.d0
       do j = 1,NDIM
          vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
       enddo
    enddo

    depth =  dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)
    xlon = atan2(vector_rotated(2),vector_rotated(1))
    ylat = asin(vector_rotated(3)/depth)
    xlon = xlon*180.d0/pi
    ylat = ylat*180.d0/pi

    !write(*,1000) x/1000.,y/1000.,z/1000.,ksi,eta,xlon,ylat,depth
!
!
    depth = REARTH - depth

  end subroutine cart2geogr

  subroutine  compute_rotation_matrix(rotation_matrix,lon_center_chunk,lat_center_chunk, chunk_azi)

    implicit none
    double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk,chunk_azi
    double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

    ! put chunk in (0,0) from (90,0)
    axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! chunk in (0,0)
    ! azimuth chunk rotation
    axe_rotation(1)=1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)
    ! put chunk in the right latitude
    axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
    call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)
    ! put chunk in the rigth longitude
    axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
    call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)
    ! compose 4 rotation matrises
    call compose4matrix(rotation_matrix,R00,R0,R1,R2)


  end subroutine compute_rotation_matrix

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! 3D rotation matrix with axis  "axe" and angle theta (degrees)
!
subroutine rotation_matrix_axe(R,axe,theta)
  implicit none
  double precision axe(3),theta,pi,deg2rad
  double precision R(3,3)
  double precision c,s,ux,uy,uz,norme_axe
  integer i,j

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! axis normalization
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! rotation axis
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! computing cos and sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! rotation matrix
  R(1,1)=(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=(ux*uz*(1.d0-c)+uy*s)

  R(2,1)=(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=(uy**2+(1.d0-uy**2)*c)
  R(2,3)=(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=(uz**2+(1.d0-uz**2)*c)


end subroutine rotation_matrix_axe

!-------------------------------------------------------------------------------------------------------------------------------------------------------

! R=R2*R1*R0
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


end subroutine compose4matrix

end module cart2geogr_mod
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
!

