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

  subroutine compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               r_top,r_bottom,ner,ilayer,ichunk,rotation_matrix,NCHUNKS,&
               INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)

  implicit none

  include "constants.h"

  double precision, dimension(NGNOD) :: xelm,yelm,zelm,offset_x,offset_y,offset_z

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  integer, intent(in) :: iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA, &
                   NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ner,ilayer,ichunk,NCHUNKS

  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,r_top,r_bottom

  logical :: INCLUDE_CENTRAL_CUBE
  integer :: NUMBER_OF_MESH_LAYERS

! local variables
  integer :: i,j,ignod

  double precision :: xi,eta,gamma,x,y,x_,y_,z,rgb,rgt,rn
  double precision :: x_bot,y_bot,z_bot
  double precision :: x_top,y_top,z_top

  double precision, dimension(NDIM) :: vector_ori,vector_rotated

  double precision :: ratio_xi, ratio_eta, fact_xi, fact_eta, &
                      fact_xi_,fact_eta_

  double precision, parameter :: PI_OVER_TWO = PI / 2.d0


! this to avoid compilation warnings
  x_=0
  y_=0

! loop on all the nodes in this element
  do ignod = 1,NGNOD

    if(ilayer == NUMBER_OF_MESH_LAYERS .and. INCLUDE_CENTRAL_CUBE) then
! case of the inner core
      ratio_xi = ((iproc_xi + offset_x(ignod)/dble(NEX_PER_PROC_XI))/dble(NPROC_XI))
      fact_xi = 2.d0*ratio_xi-1.d0

      ratio_eta = ((iproc_eta + offset_y(ignod)/dble(NEX_PER_PROC_ETA))/dble(NPROC_ETA))
      fact_eta = 2.d0*ratio_eta-1.d0

      fact_xi_ = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * fact_xi)
      fact_eta_ = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * fact_eta)

! uncomment the following lines to have more regular surface mesh (better aspect ratio for each element)
! uncomment the corresponding lines in the else condition of this if statement too.
! note that the ratio bigger_edge_size/smaller_edge_size for the surface mesh is a bit higher (1.43 vs 1.41)

!       fact_xi_= (3.d0*fact_xi+4.d0*fact_xi_)/7.d0
!       fact_eta_= (3.d0*fact_eta+4.d0*fact_eta_)/7.d0

      xi = PI_OVER_TWO*fact_xi
      eta = PI_OVER_TWO*fact_eta

      gamma = ONE / sqrt(ONE + fact_xi_**2 + fact_eta_**2)
      rgt = (r_top / R_EARTH)*gamma

! coordinates of the edge extremity on the central cube surface
      x_bot = ((r_bottom / R_EARTH) / sqrt(3.d0))* fact_xi * (1 + cos(eta)*CENTRAL_CUBE_INFLATE_FACTOR / PI)
      y_bot = ((r_bottom / R_EARTH) / sqrt(3.d0)) * fact_eta * (1 + cos(xi)*CENTRAL_CUBE_INFLATE_FACTOR / PI)
      z_bot = ((r_bottom / R_EARTH) / sqrt(3.d0)) * (1 + (cos(xi) + cos(eta))*CENTRAL_CUBE_INFLATE_FACTOR / PI)

! coordinates of the edge extremity on the ICB
      x_top = fact_xi_*rgt
      y_top = fact_eta_*rgt
      z_top = rgt

      rn = offset_z(ignod) / dble(ner)
      x = x_top*rn + x_bot*(ONE-rn)
      y = y_top*rn + y_bot*(ONE-rn)
      z = z_top*rn + z_bot*(ONE-rn)

      select case (ichunk)
        case(CHUNK_AB)
          xelm(ignod) = -y
          yelm(ignod) = x
          zelm(ignod) = z
        case(CHUNK_AB_ANTIPODE)
          xelm(ignod) = -y
          yelm(ignod) = -x
          zelm(ignod) = -z
        case(CHUNK_AC)
          xelm(ignod) = -y
          yelm(ignod) = -z
          zelm(ignod) = x
        case(CHUNK_AC_ANTIPODE)
          xelm(ignod) = -y
          yelm(ignod) = z
          zelm(ignod) = -x
        case(CHUNK_BC)
          xelm(ignod) = -z
          yelm(ignod) = y
          zelm(ignod) = x
        case(CHUNK_BC_ANTIPODE)
          xelm(ignod) = z
          yelm(ignod) = -y
          zelm(ignod) = x
        case default
          stop 'incorrect chunk number in compute_coord_main_mesh'
      end select
!       write(IMAIN,*) x,' ',y,' ',z
    else

! uncomment the following lines to have more regular surface mesh (better aspect ratio for each element)
! note that the ratio bigger_edge_size/smaller_edge_size for the surface mesh is a bit higher (1.43 vs 1.41)
!       ratio_xi = ((iproc_xi + offset_x(ignod)/dble(NEX_PER_PROC_XI))/dble(NPROC_XI))*tan(ANGULAR_WIDTH_XI_RAD/2.d0)
!       x_ = 2.d0*ratio_xi-tan(ANGULAR_WIDTH_XI_RAD/2.d0)
!       ratio_eta = ((iproc_eta + offset_y(ignod)/dble(NEX_PER_PROC_ETA))/dble(NPROC_ETA))*tan(ANGULAR_WIDTH_ETA_RAD/2.d0)
!       y_ = 2.d0*ratio_eta-tan(ANGULAR_WIDTH_ETA_RAD/2.d0)

      ratio_xi = ((iproc_xi + offset_x(ignod)/dble(NEX_PER_PROC_XI))/dble(NPROC_XI))
      x = 2.d0*ratio_xi-1

      ratio_eta = ((iproc_eta + offset_y(ignod)/dble(NEX_PER_PROC_ETA))/dble(NPROC_ETA))
      y = 2.d0*ratio_eta-1

      x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)
      y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

! uncomment the following lines to have more regular surface mesh (better aspect ratio for each element)
! note that the ratio bigger_edge_size/smaller_edge_size for the surface mesh is a bit higher (1.43 vs 1.41)
!       x= (3.d0*x_+4.d0*x)/7.d0
!       y= (3.d0*y_+4.d0*y)/7.d0

      gamma = ONE / sqrt(ONE + x*x + y*y)

      rgt = (r_top / R_EARTH)*gamma
      rgb = (r_bottom / R_EARTH)*gamma

    ! define the mesh points on the top and the bottom in the six regions of the cubed shpere
      select case (ichunk)

        case(CHUNK_AB)

          x_top = -y*rgt
          y_top = x*rgt
          z_top = rgt

          x_bot = -y*rgb
          y_bot = x*rgb
          z_bot = rgb

        case(CHUNK_AB_ANTIPODE)

          x_top = -y*rgt
          y_top = -x*rgt
          z_top = -rgt

          x_bot = -y*rgb
          y_bot = -x*rgb
          z_bot = -rgb

        case(CHUNK_AC)

          x_top = -y*rgt
          y_top = -rgt
          z_top = x*rgt

          x_bot = -y*rgb
          y_bot = -rgb
          z_bot = x*rgb

        case(CHUNK_AC_ANTIPODE)

          x_top = -y*rgt
          y_top = rgt
          z_top = -x*rgt

          x_bot = -y*rgb
          y_bot = rgb
          z_bot = -x*rgb

        case(CHUNK_BC)

          x_top = -rgt
          y_top = y*rgt
          z_top = x*rgt

          x_bot = -rgb
          y_bot = y*rgb
          z_bot = x*rgb

        case(CHUNK_BC_ANTIPODE)

          x_top = rgt
          y_top = -y*rgt
          z_top = x*rgt

          x_bot = rgb
          y_bot = -y*rgb
          z_bot = x*rgb

        case default
          stop 'incorrect chunk number in compute_coord_main_mesh'

      end select

    ! rotate the chunk to the right location if we do not mesh the full Earth
      if(NCHUNKS /= 6) then

    ! rotate bottom
        vector_ori(1) = x_bot
        vector_ori(2) = y_bot
        vector_ori(3) = z_bot
        do i = 1,NDIM
          vector_rotated(i) = ZERO
          do j = 1,NDIM
            vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
          enddo
        enddo
        x_bot = vector_rotated(1)
        y_bot = vector_rotated(2)
        z_bot = vector_rotated(3)

    ! rotate top
        vector_ori(1) = x_top
        vector_ori(2) = y_top
        vector_ori(3) = z_top
        do i = 1,NDIM
          vector_rotated(i) = ZERO
          do j = 1,NDIM
            vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
          enddo
        enddo
        x_top = vector_rotated(1)
        y_top = vector_rotated(2)
        z_top = vector_rotated(3)

      endif

    ! compute the position of the point
      rn = offset_z(ignod) / dble(ner)
      xelm(ignod) = x_top*rn + x_bot*(ONE-rn)
      yelm(ignod) = y_top*rn + y_bot*(ONE-rn)
      zelm(ignod) = z_top*rn + z_bot*(ONE-rn)

    endif
  enddo
!   if(ilayer == NUMBER_OF_MESH_LAYERS .and. INCLUDE_CENTRAL_CUBE) write(IMAIN,*)
  end subroutine compute_coord_main_mesh

!---------------------------------------------------------------------------

!! DK DK create value of arrays xgrid ygrid and zgrid in the central cube without storing them

  subroutine compute_coord_central_cube(ix,iy,iz, &
                  xgrid_central_cube,ygrid_central_cube,zgrid_central_cube, &
                  iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube,ny_central_cube,nz_central_cube,radius_cube)

  implicit none

  include "constants.h"

  integer :: ix,iy,iz,iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube,ny_central_cube,nz_central_cube

  double precision :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube,radius_cube

! local variables
  double precision :: ratio_x,ratio_y,ratio_z
  double precision :: fact_x,fact_y,fact_z,xi,eta,gamma
  double precision, parameter :: PI_OVER_TWO = PI / 2.d0

! the slice extends to the entire cube along Z
! but only to current block along X and Y
  ratio_x = (dble(iproc_xi) + dble(ix)/dble(2*nx_central_cube)) / dble(NPROC_XI)
  ratio_y = (dble(iproc_eta) + dble(iy)/dble(2*ny_central_cube)) / dble(NPROC_ETA)
  ratio_z = dble(iz)/dble(2*nz_central_cube)

  if(abs(ratio_x) > 1.001d0 .or. abs(ratio_y) > 1.001d0 .or. abs(ratio_z) > 1.001d0) stop 'wrong ratio in central cube'

! use a "flat" cubed sphere to create the central cube

! map ratio to [-1,1] and then map to real radius
! then add deformation
  fact_x = 2.d0*ratio_x-1.d0
  fact_y = 2.d0*ratio_y-1.d0
  fact_z = 2.d0*ratio_z-1.d0

  xi = PI_OVER_TWO*fact_x;
  eta = PI_OVER_TWO*fact_y;
  gamma = PI_OVER_TWO*fact_z;

  xgrid_central_cube = radius_cube * fact_x * (1 + (cos(eta)+cos(gamma))*CENTRAL_CUBE_INFLATE_FACTOR / PI);
  ygrid_central_cube = radius_cube * fact_y * (1 + (cos(xi)+cos(gamma))*CENTRAL_CUBE_INFLATE_FACTOR / PI);
  zgrid_central_cube = radius_cube * fact_z * (1 + (cos(xi)+cos(eta))*CENTRAL_CUBE_INFLATE_FACTOR / PI);

  end subroutine compute_coord_central_cube

