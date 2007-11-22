!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! create mesh quality data for the slice, to be recombined in postprocessing

  subroutine write_AVS_DX_mesh_quality_data(prname,nspec, &
                 xstore,ystore,zstore,kappastore,mustore,rhostore)

  implicit none

  include "constants.h"

  integer nspec
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: kappastore,mustore,rhostore

  integer ispec

  double precision, dimension(8) :: xelm,yelm,zelm,vp,vs

  integer iface,icorner,jcorner

  double precision vectorA_x,vectorA_y,vectorA_z
  double precision vectorB_x,vectorB_y,vectorB_z
  double precision norm_A,norm_B,angle_vectors
  double precision distmin,distmax,dist,dist1,dist2,dist3,dist4
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision vs_min,vp_max

! for stability and number of points per wavelength
  double precision highest_frequency_source,stability,points_per_wavelength

! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision percent_GLL(NGLL_MAX_STABILITY)

! topology of faces of cube for skewness
  integer faces_topo(6,6)

! processor identification
  character(len=150) prname

! data and element files identical to AVS_DXpoints.txt and AVS_DXelements.txt
! created in regular AVS or DX routine, therefore not created again here

! define percentage of smallest distance between GLL points for NGLL points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

! check that the degree is not above the threshold for list of percentages
  if(NGLLX > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

! define topology of faces of cube for skewness

! face 1
  faces_topo(1,1) = 1
  faces_topo(1,2) = 2
  faces_topo(1,3) = 6
  faces_topo(1,4) = 5

! face 2
  faces_topo(2,1) = 2
  faces_topo(2,2) = 3
  faces_topo(2,3) = 7
  faces_topo(2,4) = 6

! face 3
  faces_topo(3,1) = 4
  faces_topo(3,2) = 3
  faces_topo(3,3) = 7
  faces_topo(3,4) = 8

! face 4
  faces_topo(4,1) = 1
  faces_topo(4,2) = 5
  faces_topo(4,3) = 8
  faces_topo(4,4) = 4

! face 5
  faces_topo(5,1) = 1
  faces_topo(5,2) = 2
  faces_topo(5,3) = 3
  faces_topo(5,4) = 4

! face 6
  faces_topo(6,1) = 5
  faces_topo(6,2) = 6
  faces_topo(6,3) = 7
  faces_topo(6,4) = 8

! define wraparound for angles for skewness calculation
  faces_topo(:,5) = faces_topo(:,1)
  faces_topo(:,6) = faces_topo(:,2)

! writing mesh quality data for each element
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='unknown')

! number of elements in AVS or DX file
  write(10,*) nspec

! output global AVS or DX elements
  do ispec=1,nspec

! define the coordinates of the 8 corners of the element
     xelm(1) = xstore(1,1,1,ispec)
     yelm(1) = ystore(1,1,1,ispec)
     zelm(1) = zstore(1,1,1,ispec)
     vs(1) = sqrt(mustore(1,1,1,ispec) / rhostore(1,1,1,ispec))
     vp(1) = sqrt(kappastore(1,1,1,ispec)/rhostore(1,1,1,ispec) + 4.d0*vs(1)*vs(1)/3.d0)

     xelm(2) = xstore(NGLLX,1,1,ispec)
     yelm(2) = ystore(NGLLX,1,1,ispec)
     zelm(2) = zstore(NGLLX,1,1,ispec)
     vs(2) = sqrt(mustore(NGLLX,1,1,ispec) / rhostore(NGLLX,1,1,ispec))
     vp(2) = sqrt(kappastore(NGLLX,1,1,ispec)/rhostore(NGLLX,1,1,ispec) + 4.d0*vs(2)*vs(2)/3.d0)

     xelm(3) = xstore(NGLLX,NGLLY,1,ispec)
     yelm(3) = ystore(NGLLX,NGLLY,1,ispec)
     zelm(3) = zstore(NGLLX,NGLLY,1,ispec)
     vs(3) = sqrt(mustore(NGLLX,NGLLY,1,ispec) / rhostore(NGLLX,NGLLY,1,ispec))
     vp(3) = sqrt(kappastore(NGLLX,NGLLY,1,ispec)/rhostore(NGLLX,NGLLY,1,ispec) + 4.d0*vs(3)*vs(3)/3.d0)

     xelm(4) = xstore(1,NGLLY,1,ispec)
     yelm(4) = ystore(1,NGLLY,1,ispec)
     zelm(4) = zstore(1,NGLLY,1,ispec)
     vs(4) = sqrt(mustore(1,NGLLY,1,ispec) / rhostore(1,NGLLY,1,ispec))
     vp(4) = sqrt(kappastore(1,NGLLY,1,ispec)/rhostore(1,NGLLY,1,ispec) + 4.d0*vs(4)*vs(4)/3.d0)

     xelm(5) = xstore(1,1,NGLLZ,ispec)
     yelm(5) = ystore(1,1,NGLLZ,ispec)
     zelm(5) = zstore(1,1,NGLLZ,ispec)
     vs(5) = sqrt(mustore(1,1,NGLLZ,ispec) / rhostore(1,1,NGLLZ,ispec))
     vp(5) = sqrt(kappastore(1,1,NGLLZ,ispec)/rhostore(1,1,NGLLZ,ispec) + 4.d0*vs(5)*vs(5)/3.d0)

     xelm(6) = xstore(NGLLX,1,NGLLZ,ispec)
     yelm(6) = ystore(NGLLX,1,NGLLZ,ispec)
     zelm(6) = zstore(NGLLX,1,NGLLZ,ispec)
     vs(6) = sqrt(mustore(NGLLX,1,NGLLZ,ispec) / rhostore(NGLLX,1,NGLLZ,ispec))
     vp(6) = sqrt(kappastore(NGLLX,1,NGLLZ,ispec)/rhostore(NGLLX,1,NGLLZ,ispec) + 4.d0*vs(6)*vs(6)/3.d0)

     xelm(7) = xstore(NGLLX,NGLLY,NGLLZ,ispec)
     yelm(7) = ystore(NGLLX,NGLLY,NGLLZ,ispec)
     zelm(7) = zstore(NGLLX,NGLLY,NGLLZ,ispec)
     vs(7) = sqrt(mustore(NGLLX,NGLLY,NGLLZ,ispec) / rhostore(NGLLX,NGLLY,NGLLZ,ispec))
     vp(7) = sqrt(kappastore(NGLLX,NGLLY,NGLLZ,ispec)/rhostore(NGLLX,NGLLY,NGLLZ,ispec) + 4.d0*vs(7)*vs(7)/3.d0)

     xelm(8) = xstore(1,NGLLY,NGLLZ,ispec)
     yelm(8) = ystore(1,NGLLY,NGLLZ,ispec)
     zelm(8) = zstore(1,NGLLY,NGLLZ,ispec)
     vs(8) = sqrt(mustore(1,NGLLY,NGLLZ,ispec) / rhostore(1,NGLLY,NGLLZ,ispec))
     vp(8) = sqrt(kappastore(1,NGLLY,NGLLZ,ispec)/rhostore(1,NGLLY,NGLLZ,ispec) + 4.d0*vs(8)*vs(8)/3.d0)

! compute minimum Vs and maximum Vp
     vs_min = minval(vs(:))
     vp_max = maxval(vp(:))

! compute equiangle skewness (as defined in Fluent/Gambit manual)
     equiangle_skewness = - HUGEVAL
     do iface = 1,6
       do icorner = 1,4

! first vector of angle
         vectorA_x = xelm(faces_topo(iface,icorner)) - xelm(faces_topo(iface,icorner+1))
         vectorA_y = yelm(faces_topo(iface,icorner)) - yelm(faces_topo(iface,icorner+1))
         vectorA_z = zelm(faces_topo(iface,icorner)) - zelm(faces_topo(iface,icorner+1))

! second vector of angle
         vectorB_x = xelm(faces_topo(iface,icorner+2)) - xelm(faces_topo(iface,icorner+1))
         vectorB_y = yelm(faces_topo(iface,icorner+2)) - yelm(faces_topo(iface,icorner+1))
         vectorB_z = zelm(faces_topo(iface,icorner+2)) - zelm(faces_topo(iface,icorner+1))

! norm of vectors A and B
         norm_A = dsqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
         norm_B = dsqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

! angle formed by the two vectors
         angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

! compute equiangle skewness
         equiangle_skewness = dmax1(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

       enddo
     enddo

! compute edge aspect ratio using the 8 corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     do icorner = 1,8
       do jcorner = icorner + 1,8
         dist = dsqrt((xelm(jcorner) - xelm(icorner))**2 + (yelm(jcorner) - yelm(icorner))**2 + (zelm(jcorner) - zelm(icorner))**2)
         distmin = dmin1(distmin,dist)
         distmax = dmax1(distmax,dist)
       enddo
     enddo
     edge_aspect_ratio = distmax / distmin

! compute stability value and number of points per shortest S wavelength
! stability is multiplied by time step in check_mesh_quality_AVS_DX.f90
!! DK DK incomplete formula for now, frequency of source hard wired
!! DK DK should instead take from cmt_file, read and multiply
!! DK DK directly later in check_mesh_quality.f90
   highest_frequency_source = 1.
   stability = vp_max / (distmin * percent_GLL(NGLLX))
   points_per_wavelength = NGLLX * vs_min / (distmax * highest_frequency_source)

! compute diagonal aspect ratio
     dist1 = dsqrt((xelm(1) - xelm(7))**2 + (yelm(1) - yelm(7))**2 + (zelm(1) - zelm(7))**2)
     dist2 = dsqrt((xelm(2) - xelm(8))**2 + (yelm(2) - yelm(8))**2 + (zelm(2) - zelm(8))**2)
     dist3 = dsqrt((xelm(3) - xelm(5))**2 + (yelm(3) - yelm(5))**2 + (zelm(3) - zelm(5))**2)
     dist4 = dsqrt((xelm(4) - xelm(6))**2 + (yelm(4) - yelm(6))**2 + (zelm(4) - zelm(6))**2)
     distmin = dmin1(distmin,dist1,dist2,dist3,dist4)
     distmax = dmax1(distmax,dist1,dist2,dist3,dist4)
     diagonal_aspect_ratio = distmax / distmin

! write mesh quality information for each element
   write(10,*) ispec,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength

  enddo

  close(10)

  end subroutine write_AVS_DX_mesh_quality_data

