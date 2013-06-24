!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

! read a 3D CUBIT mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element


! Dimitri Komatitsch, University of Pau, France, March 2009.
! Modified by Pieyre Le Loher

!! DK DK
!! DK DK this routine could be improved by computing the mean in addition to min and max of ratios
!! DK DK

  subroutine check_mesh_quality(myrank,VP_MAX,NGLOB,NSPEC,x,y,z,ibool, &
                                CREATE_VTK_FILES,prname)

  implicit none

  include "constants.h"

  integer :: myrank

  double precision :: VP_MAX           ! maximum vp in volume block id 3

  integer :: NGLOB                    ! number of nodes
  integer :: NSPEC

  double precision, dimension(NGLOB) :: x,y,z
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC) :: ibool

  logical :: CREATE_VTK_FILES
  character(len=256) prname

  ! local parameters
  integer :: ispec,ispec_min_edge_length,ispec_max_edge_length,ispec_max_skewness, &
       ispec_max_skewness_MPI,skewness_max_rank,NSPEC_ALL_SLICES

  ! for quality of mesh
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision :: equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision :: equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision :: distance_min,distance_max
  double precision :: distmin,distmax

  double precision :: equiangle_skewness_max_MPI,edge_aspect_ratio_max_MPI,diagonal_aspect_ratio_max_MPI
  double precision :: distance_min_MPI,distance_max_MPI
  ! for stability
  double precision :: dt_suggested,dt_suggested_max,dt_suggested_max_MPI
  double precision :: stability,stability_min,stability_max,max_CFL_stability_limit
  double precision :: stability_max_MPI

  ! For MPI maxloc reduction
  double precision, dimension(2) :: buf_maxloc_send,buf_maxloc_recv

  ! for histogram
  integer, parameter :: NCLASS = 20
  integer classes_skewness(0:NCLASS-1)
  integer classes_skewnessMPI(0:NCLASS-1)
  integer :: iclass
  double precision :: current_percent,total_percent

  ! to export elements that have a certain skewness range to OpenDX
  !integer :: ntotspecAVS_DX
  !logical :: USE_OPENDX

  !character(len=256):: line
  integer,dimension(1) :: tmp_ispec_max_skewness,tmp_ispec_max_skewness_MPI

  ! debug: for vtk output
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp1
  integer:: ier,ipoin

  ! parameters read from parameter file
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,NGNOD,NGNOD2D
  double precision DT
  double precision HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML
  logical ATTENUATION,USE_OLSEN_ATTENUATION, &
          APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE
  logical STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE
  logical ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  logical PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID
  character(len=256) LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH
  integer NPROC
  integer MOVIE_TYPE,IMODEL


  if (myrank == 0) then
     write(IMAIN,*) '**************************'
     write(IMAIN,*) 'Checking mesh quality'
     write(IMAIN,*) '**************************'
     write(IMAIN,*)
     write(IMAIN,*) 'start computing the minimum and maximum edge size'
  endif

  ! read the parameter file
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                        PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  if(NGNOD /= 8) stop 'error: check_mesh_quality only supports NGNOD == 8 for now'

  ! ************* compute min and max of skewness and ratios ******************

  ! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL
  stability_min = + HUGEVAL
  distance_min = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL
  stability_max = - HUGEVAL
  distance_max = - HUGEVAL
  dt_suggested_max = HUGEVAL

  ispec_min_edge_length = -1
  ispec_max_edge_length = -1
  ispec_max_skewness = -1

  ! debug: for vtk output
  if( CREATE_VTK_FILES ) then
    allocate(tmp1(NSPEC),stat=ier)
    if( ier /= 0 ) stop 'error allocating array tmp'
    tmp1(:) = 0.0
  endif


  ! loop on all the elements
  do ispec = 1,NSPEC

     call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,dt_suggested, &
                                    equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio, &
                                    stability,distmin,distmax)

     ! store element number in which the edge of minimum or maximum length is located
     if(distmin < distance_min) ispec_min_edge_length = ispec
     if(distmax > distance_max) ispec_max_edge_length = ispec
     if(equiangle_skewness > equiangle_skewness_max) ispec_max_skewness = ispec

     if( CREATE_VTK_FILES ) then
        if(CUSTOM_REAL == SIZE_REAL) then
          tmp1(ispec) = sngl(equiangle_skewness)
        else
          tmp1(ispec) = equiangle_skewness
        endif
     endif

     ! compute minimum and maximum of quality numbers
     equiangle_skewness_min = min(equiangle_skewness_min,equiangle_skewness)
     edge_aspect_ratio_min = min(edge_aspect_ratio_min,edge_aspect_ratio)
     diagonal_aspect_ratio_min = min(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
     stability_min = min(stability_min,stability)
     distance_min = min(distance_min,distmin)

     equiangle_skewness_max = max(equiangle_skewness_max,equiangle_skewness)
     edge_aspect_ratio_max = max(edge_aspect_ratio_max,edge_aspect_ratio)
     diagonal_aspect_ratio_max = max(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
     stability_max = max(stability_max,stability)
     dt_suggested_max = min(dt_suggested_max,dt_suggested)
     distance_max = max(distance_max,distmax)

  enddo

  call sync_all()

  call min_all_dp(distance_min,distance_min_MPI)
  call max_all_dp(distance_max,distance_max_MPI)

  buf_maxloc_send(1) = equiangle_skewness_max
  buf_maxloc_send(2) = myrank
  call maxloc_all_dp(buf_maxloc_send,buf_maxloc_recv)
  equiangle_skewness_max_MPI = buf_maxloc_recv(1)
  skewness_max_rank = int(buf_maxloc_recv(2))


  call max_all_dp(stability,stability_max_MPI)
  call min_all_dp(dt_suggested_max,dt_suggested_max_MPI)

  call max_all_dp(edge_aspect_ratio_max,edge_aspect_ratio_max_MPI)
  call max_all_dp(diagonal_aspect_ratio_max,diagonal_aspect_ratio_max_MPI)

  if((myrank == skewness_max_rank) .and. (myrank /= 0)) then
     tmp_ispec_max_skewness(1) = ispec_max_skewness
     call send_i_t(tmp_ispec_max_skewness,1,0)
  endif

  if(myrank == 0) then

     if(skewness_max_rank /= myrank) then
        call recv_i_t(tmp_ispec_max_skewness_MPI,1,skewness_max_rank)
        ispec_max_skewness_MPI = tmp_ispec_max_skewness_MPI(1)
     else
        ispec_max_skewness_MPI = ispec_max_skewness
     endif

  write(IMAIN,*) 'done processing '

  write(IMAIN,*)
  write(IMAIN,*) '------------'
  write(IMAIN,*) 'mesh quality parameter definitions:'
  write(IMAIN,*)
  write(IMAIN,*) 'equiangle skewness: 0. perfect,  1. bad'
  write(IMAIN,*) 'skewness max deviation angle: 0. perfect,  90. bad'
  write(IMAIN,*) 'edge aspect ratio: 1. perfect,  above 1. gives stretching factor'
  write(IMAIN,*) 'diagonal aspect ratio: 1. perfect,  above 1. gives stretching factor'
  write(IMAIN,*) '------------'

  write(IMAIN,*)
  write(IMAIN,*) 'minimum length of an edge in the whole mesh (m) = ',distance_min_MPI!,' in element ',ispec_min_edge_length
  write(IMAIN,*)
  write(IMAIN,*) 'maximum length of an edge in the whole mesh (m) = ',distance_max_MPI!,' in element ',ispec_max_edge_length
  write(IMAIN,*)
  write(IMAIN,*) '***'
  write(IMAIN,*) '*** max equiangle skewness = ',equiangle_skewness_max_MPI,' in element ',ispec_max_skewness_MPI, &
       ' of slice ',skewness_max_rank
  write(IMAIN,*) '***'
  write(IMAIN,*)
  write(IMAIN,*) 'max deviation angle from a right angle (90 degrees) is therefore = ',90.*equiangle_skewness_max_MPI
  write(IMAIN,*)
  write(IMAIN,*) 'worst angle in the mesh is therefore ',90.*(1. - equiangle_skewness_max_MPI)
  write(IMAIN,*) 'or ',180. - 90.*(1. - equiangle_skewness_max_MPI),' degrees'
  write(IMAIN,*)
  write(IMAIN,*) 'max edge aspect ratio = ',edge_aspect_ratio_max_MPI
  write(IMAIN,*)
  write(IMAIN,*) 'max diagonal aspect ratio = ',diagonal_aspect_ratio_max_MPI
  write(IMAIN,*)
  write(IMAIN,*) '***'
  write(IMAIN,'(a50,f13.8)') ' *** Maximum suggested time step for simulation = ',dt_suggested_max_MPI
  write(IMAIN,*) '***'
  write(IMAIN,*) '*** max stability = ',stability_max_MPI
  write(IMAIN,*) '*** computed using VP_MAX = ',VP_MAX
  write(IMAIN,*) '***'

  ! max stability CFL value
!! DK DK we could probably increase this, now that the Stacey conditions have been fixed
  max_CFL_stability_limit = 0.55d0 !! DK DK increased this    0.48d0

  if(stability_max_MPI >= max_CFL_stability_limit) then
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) ' WARNING, that value is above the upper CFL limit of ',max_CFL_stability_limit
     write(IMAIN,*) 'therefore the run should be unstable'
     write(IMAIN,*) 'You can try to reduce the time step'
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) '*********************************************'
  else
     write(IMAIN,*) 'that value is below the upper CFL limit of ',max_CFL_stability_limit
     write(IMAIN,*) 'therefore the run should be stable'
  endif
  write(IMAIN,*)

  ! create statistics about mesh quality
  write(IMAIN,*) 'creating histogram and statistics of mesh quality'
  endif

  ! erase histogram of skewness
  classes_skewness(:) = 0

  ! loop on all the elements
  do ispec = 1,NSPEC

     call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,dt_suggested, &
                                      equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio, &
                                      stability,distmin,distmax)

     ! store skewness in histogram
     iclass = int(equiangle_skewness * dble(NCLASS))
     if(iclass < 0) iclass = 0
     if(iclass > NCLASS-1) iclass = NCLASS-1
     classes_skewness(iclass) = classes_skewness(iclass) + 1

  enddo

  ! sum skewness results in all processes
  do iclass = 0,NCLASS-1
     call sum_all_i(classes_skewness(iclass),classes_skewnessMPI(iclass))
  enddo

  call sum_all_i(NSPEC,NSPEC_ALL_SLICES)

  if(myrank == 0) then
  ! create histogram of skewness and save in Gnuplot file
  write(IMAIN,*)
  write(IMAIN,*) 'histogram of skewness (0. good - 1. bad):'
  write(IMAIN,*)
  total_percent = 0.
  open(unit=14,file=OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))//'mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASS-1
     current_percent = 100.*dble(classes_skewnessMPI(iclass))/dble(NSPEC_ALL_SLICES)
     total_percent = total_percent + current_percent
     write(IMAIN,*) real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewnessMPI(iclass),&
          ' ',sngl(current_percent),' %'
     write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

  ! create script for Gnuplot histogram file
  open(unit=14,file=OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))// &
       'plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term x11'
  write(14,*) '#set term gif'
  write(14,*) '#set output "mesh_quality_histogram.gif"'
  write(14,*)
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASS)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  ! display warning if maximum skewness is too high
  if(equiangle_skewness_max >= 0.75d0) then
     write(IMAIN,*)
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) ' WARNING, mesh is bad (max skewness >= 0.75)'
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*) '*********************************************'
     write(IMAIN,*)
  endif

  if(total_percent < 99.9d0 .or. total_percent > 100.1d0) then
     write(IMAIN,*) 'total percentage = ',total_percent,' %'
     stop 'total percentage should be 100%'
  endif

  endif

  ! debug: for vtk output
  if( CREATE_VTK_FILES ) then
    ! vtk file output
    open(66,file=prname(1:len_trim(prname))//'skewness.vtk',status='unknown')
    write(66,'(a)') '# vtk DataFile Version 3.1'
    write(66,'(a)') 'material model VTK file'
    write(66,'(a)') 'ASCII'
    write(66,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(66, '(a,i12,a)') 'POINTS ', nglob, ' float'
    do ipoin = 1,nglob
      write(66,*) sngl(x(ipoin)),sngl(y(ipoin)),sngl(z(ipoin))
    enddo
    write(66,*) ""

    ! note: indices for vtk start at 0
    write(66,'(a,i12,i12)') "CELLS ",nspec,nspec*9
    do ispec=1,nspec
      write(66,'(9i12)') 8, &
            ibool(1,ispec)-1,ibool(2,ispec)-1,ibool(4,ispec)-1,ibool(3,ispec)-1,&
            ibool(5,ispec)-1,ibool(6,ispec)-1,ibool(8,ispec)-1,ibool(7,ispec)-1
    enddo
    write(66,*) ""

    ! type: hexahedra
    write(66,'(a,i12)') "CELL_TYPES ",nspec
    write(66,'(6i12)') (12,ispec=1,nspec)
    write(66,*) ""

    write(66,'(a,i12)') "CELL_DATA ",nspec
    write(66,'(a)') "SCALARS skewness float"
    write(66,'(a)') "LOOKUP_TABLE default"
    do ispec = 1,nspec
      write(66,*) tmp1(ispec)
    enddo
    write(66,*) ""
    close(66)

    deallocate(tmp1)
  endif


  end subroutine check_mesh_quality

!
!=====================================================================
!

! create mesh quality data for a given 3D spectral element

  subroutine create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,dt_suggested, &
                                        equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio, &
                                        stability,distmin,distmax)

  implicit none

  include "constants.h"
  include "constants_meshfem3D.h"

  integer :: NSPEC,NGLOB

  double precision, dimension(NGLOB) :: x,y,z

  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC) :: ibool
  integer :: ispec

  ! for stability
  double precision :: VP_MAX,dt_suggested
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision :: stability,distmin,distmax

  ! local parameters
  double precision, dimension(NGNOD_EIGHT_CORNERS) :: xelm,yelm,zelm
  double precision :: vectorA_x,vectorA_y,vectorA_z
  double precision :: vectorB_x,vectorB_y,vectorB_z
  double precision :: norm_A,norm_B,angle_vectors
  double precision :: dist,dist1,dist2,dist3,dist4

  ! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision, dimension(NGLL_MAX_STABILITY) :: percent_GLL

  ! topology of faces of cube for skewness
  integer faces_topo(6,6)
  integer :: iface,icorner,i

  integer,parameter :: true_NGLLX = 5


  ! store the corners of this element for the skewness routine
  do i = 1,NGNOD_EIGHT_CORNERS
     xelm(i) = x(ibool(i,ispec))
     yelm(i) = y(ibool(i,ispec))
     zelm(i) = z(ibool(i,ispec))
  enddo

  ! define percentage of smallest distance between GLL points for NGLL points
  ! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(1) = 100.d0
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
  if(NGLLX_M > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

  ! define topology of faces of cube for skewness

  ! face 1
  faces_topo(1,1) = 5
  faces_topo(1,2) = 1
  faces_topo(1,3) = 2
  faces_topo(1,4) = 6

  ! face 2
  faces_topo(2,1) = 1
  faces_topo(2,2) = 3
  faces_topo(2,3) = 4
  faces_topo(2,4) = 2

  ! face 3
  faces_topo(3,1) = 7
  faces_topo(3,2) = 3
  faces_topo(3,3) = 4
  faces_topo(3,4) = 8

  ! face 4
  faces_topo(4,1) = 5
  faces_topo(4,2) = 6
  faces_topo(4,3) = 8
  faces_topo(4,4) = 7

  ! face 5
  faces_topo(5,1) = 5
  faces_topo(5,2) = 1
  faces_topo(5,3) = 3
  faces_topo(5,4) = 7

  ! face 6
  faces_topo(6,1) = 6
  faces_topo(6,2) = 2
  faces_topo(6,3) = 4
  faces_topo(6,4) = 8

  ! define wraparound for angles for skewness calculation
  faces_topo(:,5) = faces_topo(:,1)
  faces_topo(:,6) = faces_topo(:,2)

  ! compute equiangle skewness (as defined in Fluent/Gambit manual)
  ! and compute edge aspect ratio using the corners of the element
  distmin = + HUGEVAL
  distmax = - HUGEVAL
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
        norm_A = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
        norm_B = sqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

        ! angle formed by the two vectors
        angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

        ! compute equiangle skewness
        equiangle_skewness = max(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

        ! compute min and max size of an edge
        dist = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)

        distmin = min(distmin,dist)
        distmax = max(distmax,dist)

     enddo
  enddo

  ! compute edge aspect ratio
  edge_aspect_ratio = distmax / distmin

  !stability = delta_t * VP_MAX / (distmin * percent_GLL(true_NGLLX))

  dt_suggested = ((1.d0 - 0.02d0)*0.48d0) * (distmin * percent_GLL(true_NGLLX)) / VP_MAX
  stability = dt_suggested * VP_MAX / (distmin * percent_GLL(true_NGLLX))

  ! compute diagonal aspect ratio
  dist1 = sqrt((xelm(5) - xelm(4))**2 + (yelm(5) - yelm(4))**2 + (zelm(5) - zelm(4))**2)
  dist2 = sqrt((xelm(1) - xelm(8))**2 + (yelm(1) - yelm(8))**2 + (zelm(1) - zelm(8))**2)
  dist3 = sqrt((xelm(3) - xelm(6))**2 + (yelm(3) - yelm(6))**2 + (zelm(3) - zelm(6))**2)
  dist4 = sqrt((xelm(7) - xelm(2))**2 + (yelm(7) - yelm(2))**2 + (zelm(7) - zelm(2))**2)
  diagonal_aspect_ratio = max(dist1,dist2,dist3,dist4) / min(dist1,dist2,dist3,dist4)

  end subroutine create_mesh_quality_data_3D

