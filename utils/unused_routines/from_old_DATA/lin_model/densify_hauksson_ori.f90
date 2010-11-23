
program jfdkfd

  !! DK DK regrid Hauksson on regular grid in lat-long for So-Cal

  implicit none

  include "../../constants.h"

  ! UTM projection zone
  integer, parameter :: UTM_PROJECTION_ZONE = 11
  logical, parameter :: SUPPRESS_UTM_PROJECTION = .false.

  !! DK DK densification factor
  integer, parameter :: N_DENSIF = 10

  integer, parameter :: NX_ORI_HAUKSSON = 29, NY_ORI_HAUKSSON = 45

  double precision, dimension(NX_ORI_HAUKSSON,NY_ORI_HAUKSSON) :: rlat_ori,rlon_ori
  double precision, dimension(16,NX_ORI_HAUKSSON,NY_ORI_HAUKSSON) :: value_ori

  integer i,j,iloop,jloop,ipoin,k
  double precision xi,eta
  double precision rx1,ry1,v1(16),rx2,ry2,v2(16),rx3,ry3,v3(16),rx4,ry4,v4(16)
  double precision rx,ry,val(16),utm_x,utm_y

  ! read original titled non-evenly spaced Hauksson grid
  do j=1,NY_ORI_HAUKSSON
     do i=1,NX_ORI_HAUKSSON
        read(*,*) rlon_ori(i,j),rlat_ori(i,j),(value_ori(k,i,j),k=1,16)
     enddo
  enddo

  ! write denser grid using interpolation
  ipoin = 0
  do j=1,NY_ORI_HAUKSSON-1
     do i=1,NX_ORI_HAUKSSON-1
        rx1 = rlon_ori(i,j)
        ry1 = rlat_ori(i,j)
        v1(:) = value_ori(:,i,j)

        rx2 = rlon_ori(i+1,j)
        ry2 = rlat_ori(i+1,j)
        v2(:) = value_ori(:,i+1,j)

        rx3 = rlon_ori(i+1,j+1)
        ry3 = rlat_ori(i+1,j+1)
        v3(:) = value_ori(:,i+1,j+1)

        rx4 = rlon_ori(i,j+1)
        ry4 = rlat_ori(i,j+1)
        v4(:) = value_ori(:,i,j+1)

        do jloop = 1,N_DENSIF
           do iloop = 1,N_DENSIF
              ipoin = ipoin + 1
              xi = dble(iloop-1)/dble(N_DENSIF-1)
              eta = dble(jloop-1)/dble(N_DENSIF-1)
              ! use bilinear interpolation to compute value
              ! this is approximate because cell is not exactly a rectangle, but that's okay
              ! because it's almost a rectangle, and Hauksson model is approximate
              rx = rx1*(1.d0-xi)*(1.d0-eta) + rx2*xi*(1.d0-eta) + rx3*xi*eta + rx4*(1.d0-xi)*eta
              ry = ry1*(1.d0-xi)*(1.d0-eta) + ry2*xi*(1.d0-eta) + ry3*xi*eta + ry4*(1.d0-xi)*eta
              val(:) = v1(:)*(1.d0-xi)*(1.d0-eta) + v2(:)*xi*(1.d0-eta) + v3(:)*xi*eta + v4(:)*(1.d0-xi)*eta

              call utm_geo(rx,ry,utm_x,utm_y,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

              print *,utm_x,utm_y,val(:)
           enddo
        enddo

     enddo
  enddo

  open(unit=13,file='header.h',status='unknown')
  write(13,*) 'integer, parameter :: NLINES_HAUKSSON_DENSER = ',ipoin
  close(13)

end program jfdkfd

!! DK DK add UTM projection routine
include "../../utm_geo.f90"

