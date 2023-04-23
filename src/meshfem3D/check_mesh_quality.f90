!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

! read a 3D CUBIT mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Pau, France, March 2009.
! Modified by Pieyre Le Loher

!! DK DK
!! DK DK this routine could be improved by computing the mean in addition to min and max of ratios
!! DK DK

  subroutine check_mesh_quality(VP_MAX,NGLOB,NSPEC,x,y,z,ibool, &
                                CREATE_VTK_FILES,prname)

  use constants, only: IMAIN,myrank,MAX_STRING_LEN,CUSTOM_REAL,HUGEVAL,OUTPUT_FILES
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  implicit none

  double precision :: VP_MAX           ! maximum vp in volume block id 3

  integer,intent(in) :: NGLOB                    ! number of nodes
  integer,intent(in) :: NSPEC

  double precision, dimension(NGLOB),intent(in) :: x,y,z
  integer,dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),intent(in) :: ibool

  logical,intent(in) :: CREATE_VTK_FILES
  character(len=MAX_STRING_LEN) :: prname

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
  integer :: classes_skewness(0:NCLASS-1)
  integer :: classes_skewnessMPI(0:NCLASS-1)
  integer :: iclass
  double precision :: current_percent,total_percent

  integer,dimension(1) :: tmp_ispec_max_skewness,tmp_ispec_max_skewness_MPI

  ! debug: for vtk output
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp1
  integer:: ier

  ! user output
  if (myrank == 0) then
     write(IMAIN,*) '**************************'
     write(IMAIN,*) 'Checking mesh quality'
     write(IMAIN,*) '**************************'
     write(IMAIN,*)
     write(IMAIN,*) 'start computing the minimum and maximum edge size'
     call flush_IMAIN()
  endif

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
  if (CREATE_VTK_FILES) then
    allocate(tmp1(NSPEC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1345')
    if (ier /= 0) stop 'error allocating array tmp'
    tmp1(:) = 0.0
  endif


  ! loop on all the elements
  do ispec = 1,NSPEC
     call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,dt_suggested, &
                                      equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio, &
                                      stability,distmin,distmax)

     ! store element number in which the edge of minimum or maximum length is located
     if (distmin < distance_min) ispec_min_edge_length = ispec
     if (distmax > distance_max) ispec_max_edge_length = ispec
     if (equiangle_skewness > equiangle_skewness_max) ispec_max_skewness = ispec

     if (CREATE_VTK_FILES) then
       tmp1(ispec) = real(equiangle_skewness,kind=CUSTOM_REAL)
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

  call synchronize_all()

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

  if ((myrank == skewness_max_rank) .and. (myrank /= 0)) then
     tmp_ispec_max_skewness(1) = ispec_max_skewness
     call send_i_t(tmp_ispec_max_skewness,1,0)
  endif

  if (myrank == 0) then

    if (skewness_max_rank /= myrank) then
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
    write(IMAIN,*) '*** Max CFL stability condition of the time scheme (must be below about 0.55 or so) = ',stability_max_MPI
    write(IMAIN,*) '*** computed using the maximum P wave velocity = ',VP_MAX
    write(IMAIN,*) '***'
    call flush_IMAIN()

    ! max stability CFL value
    !! DK DK we could probably increase this, now that the Stacey conditions have been fixed
    max_CFL_stability_limit = 0.55d0 !! DK DK increased this    0.48d0

    if (stability_max_MPI >= max_CFL_stability_limit) then
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
    write(IMAIN,*) 'creating histogram of mesh quality'
    call flush_IMAIN()
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
     if (iclass < 0) iclass = 0
     if (iclass > NCLASS-1) iclass = NCLASS-1
     classes_skewness(iclass) = classes_skewness(iclass) + 1

  enddo

  ! sum skewness results in all processes
  do iclass = 0,NCLASS-1
     call sum_all_i(classes_skewness(iclass),classes_skewnessMPI(iclass))
  enddo

  call sum_all_i(NSPEC,NSPEC_ALL_SLICES)

  if (myrank == 0) then
    ! create histogram of skewness and save in Gnuplot file
    write(IMAIN,*)
    write(IMAIN,*) 'histogram of skewness (0. good - 1. bad):'
    write(IMAIN,*)
    call flush_IMAIN()

    total_percent = 0.

    ! histogram data file
    filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'mesh_quality_histogram.txt'

    open(unit=14,file=trim(filename),status='unknown')
    do iclass = 0,NCLASS-1
      current_percent = 100.*dble(classes_skewnessMPI(iclass))/dble(NSPEC_ALL_SLICES)
      total_percent = total_percent + current_percent
      write(IMAIN,*) real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewnessMPI(iclass), &
                     ' ',sngl(current_percent),' %'
      write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
    enddo
    write(IMAIN,*)
    close(14)

    ! create script for Gnuplot histogram file
    filename = OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//'plot_mesh_quality_histogram.gnu'

    open(unit=14,file=trim(filename),status='unknown')
    write(14,*) 'set terminal x11'
    write(14,*) '#set terminal wxt'
    write(14,*) '#set terminal gif'
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
    if (equiangle_skewness_max >= 0.75d0) then
      write(IMAIN,*)
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*) ' WARNING, mesh is bad (max skewness >= 0.75)'
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*) '*********************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    if (total_percent < 99.9d0 .or. total_percent > 100.1d0) then
       write(IMAIN,*) 'total percentage = ',total_percent,' %'
       stop 'total percentage should be 100%'
    endif
    call flush_IMAIN()
  endif

  ! debug: for vtk output
  if (CREATE_VTK_FILES) then
    filename = prname(1:len_trim(prname))//'skewness.vtk'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'plotting skewness to VTK-file: ',trim(filename)
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! vtk file output
    call write_VTK_data_elem_cr_meshfem(nspec,nglob,x,y,z,ibool,tmp1,filename)

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

  use constants, only: NGNOD_EIGHT_CORNERS,PI,HUGEVAL
  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  implicit none

  integer,intent(in) :: NSPEC,NGLOB

  double precision, dimension(NGLOB),intent(in) :: x,y,z

  integer, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,NSPEC),intent(in) :: ibool
  integer,intent(in) :: ispec

  ! for stability
  double precision :: VP_MAX,dt_suggested
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision :: stability,distmin,distmax

  ! local parameters
  double precision :: vectorA_x,vectorA_y,vectorA_z
  double precision :: vectorB_x,vectorB_y,vectorB_z
  double precision :: norm_A,norm_B,angle_vectors,argument_of_arccos
  double precision :: dist,dist1,dist2,dist3,dist4

  ! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision, dimension(NGLL_MAX_STABILITY) :: percent_GLL

  ! topology of faces of cube for skewness
  integer :: faces_topo(6,6)
  integer :: iface,icorner

  ! topology of the elements
  double precision, dimension(NGNOD_EIGHT_CORNERS) :: xelm,yelm,zelm
  integer, dimension(NGNOD_EIGHT_CORNERS) :: anchor_iax,anchor_iay,anchor_iaz
  integer :: ia,iglob

  integer,parameter :: true_NGLLX = 5

  ! sets up node addressing
  call hex_nodes_anchor_ijk_NGLL(NGNOD_EIGHT_CORNERS,anchor_iax,anchor_iay,anchor_iaz,NGLLX_M,NGLLY_M,NGLLZ_M)

  ! store the corners of this element for the skewness routine
  do ia = 1,NGNOD_EIGHT_CORNERS
    ! assumes NGLLX_M == NGLLY_M == NGLLZ_M == 2 == NGNOD_EIGHT_CORNERS
    !xelm(ia) = x(ibool(i,ispec))
    !yelm(ia) = y(ibool(i,ispec))
    !zelm(ia) = z(ibool(i,ispec))
    ! general
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
    xelm(ia) = x(iglob)
    yelm(ia) = y(iglob)
    zelm(ia) = z(iglob)
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
  if (NGLLX_M > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

  ! define topology of faces of cube for skewness
  ! (see hex_nodes.f90)
  !
  !           5 * -------- * 8
  !            /|         /|
  !         6 * -------- *7|
  !           |1* - - -  | * 4                 z
  !           |/         |/                   |
  !           * -------- *                    o--> y
  !          2           3                 x /
  !
  ! face 1 - left
  faces_topo(1,1) = 5
  faces_topo(1,2) = 1
  faces_topo(1,3) = 2
  faces_topo(1,4) = 6

  ! face 2 - bottom
  faces_topo(2,1) = 1
  faces_topo(2,2) = 4
  faces_topo(2,3) = 3
  faces_topo(2,4) = 2

  ! face 3 - right
  faces_topo(3,1) = 7
  faces_topo(3,2) = 3
  faces_topo(3,3) = 4
  faces_topo(3,4) = 8

  ! face 4 - top
  faces_topo(4,1) = 5
  faces_topo(4,2) = 6
  faces_topo(4,3) = 7
  faces_topo(4,4) = 8

  ! face 5 - back
  faces_topo(5,1) = 5
  faces_topo(5,2) = 1
  faces_topo(5,3) = 4
  faces_topo(5,4) = 8

  ! face 6 - front
  faces_topo(6,1) = 6
  faces_topo(6,2) = 2
  faces_topo(6,3) = 3
  faces_topo(6,4) = 7

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
        argument_of_arccos = (vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B)

        ! debug
        !if (iface == 2) &
        !  print *,'debug:',argument_of_arccos,dacos(argument_of_arccos),(2*dacos(argument_of_arccos)-PI)/PI,equiangle_skewness, &
        !        'icorner',icorner,iface,norm_A,norm_B

        ! compute equiangle skewness
        if (abs(argument_of_arccos) <= 0.9999999d0) then
           angle_vectors = dacos(argument_of_arccos)
           equiangle_skewness = max(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)
        else
           angle_vectors = 0.d0
           equiangle_skewness = 1.d0
        endif

        ! compute min and max size of an edge
        dist = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)

        distmin = min(distmin,dist)
        distmax = max(distmax,dist)

     enddo
  enddo

  ! compute edge aspect ratio
  edge_aspect_ratio = distmax / distmin

  dt_suggested = ((1.d0 - 0.02d0)*0.48d0) * (distmin * percent_GLL(true_NGLLX)) / VP_MAX
  stability = dt_suggested * VP_MAX / (distmin * percent_GLL(true_NGLLX))

  ! compute diagonal aspect ratio
  dist1 = sqrt((xelm(5) - xelm(4))**2 + (yelm(5) - yelm(4))**2 + (zelm(5) - zelm(4))**2)
  dist2 = sqrt((xelm(1) - xelm(8))**2 + (yelm(1) - yelm(8))**2 + (zelm(1) - zelm(8))**2)
  dist3 = sqrt((xelm(3) - xelm(6))**2 + (yelm(3) - yelm(6))**2 + (zelm(3) - zelm(6))**2)
  dist4 = sqrt((xelm(7) - xelm(2))**2 + (yelm(7) - yelm(2))**2 + (zelm(7) - zelm(2))**2)
  diagonal_aspect_ratio = max(dist1,dist2,dist3,dist4) / min(dist1,dist2,dist3,dist4)

  end subroutine create_mesh_quality_data_3D

