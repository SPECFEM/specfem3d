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

! read an external mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Pau, France, March 2009 and CNRS, Marseille, France, June 2015 and February 2016.

  program check_mesh_quality

  use constants

  implicit none

!------------------------------------------------------------------------------------------------

  integer, parameter :: NGNOD = 8                        ! number of control nodes for hexahedral elements (can only be 8 or 27)

  character(len=*), parameter :: nodes_coords_file = 'MESH/nodes_coords_file'
  character(len=*), parameter :: mesh_file = 'MESH/mesh_file'

  double precision, parameter :: delta_t = 1.d0          ! fictitious value used for compatibility with the subroutine call only
  double precision, parameter :: VP_MAX  = 1.d0          ! fictitious value used for compatibility with the subroutine call only

!------------------------------------------------------------------------------------------------

  integer :: NGLOB                    ! number of nodes
  integer :: NSPEC                    ! number of elements

  double precision, dimension(:), allocatable :: x,y,z

  integer, dimension(:,:), allocatable :: ibool

  integer :: i,ispec,iread,iformat,iabove_or_below,ispec_min_edge_length,ispec_max_edge_length, &
             ispec_begin,ispec_end,ier

  double precision :: xtmp,ytmp,ztmp
  integer :: n1,n2,n3,n4,n5,n6,n7,n8

! for quality of mesh
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,value_to_use
  double precision :: equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision :: equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision :: threshold_AVS_DX_min,threshold_AVS_DX_max,threshold_AVS_DX_size_to_use
  double precision :: distance_min,distance_max,distance_mean
  double precision :: distmin,distmax,distmean,min_of_distmean,max_of_distmean

! for stability
  double precision :: stability

! for histogram
  integer, parameter :: NCLASS = 40
  integer, dimension(0:NCLASS-1) :: classes_of_histogram_skewness,classes_of_histogram_meansize
  integer :: iclass
  double precision :: current_percent,total_percent

! to export elements that have a certain range to OpenDX
  integer :: ntotspecAVS_DX
  logical :: DISPLAY_HISTOGRAM_DISTMEAN

  if (NGNOD /= 8) then
    print *,'error: check_mesh_quality only supports NGNOD == 8 for now'
    stop 'thus if NGNOD == 27, just run the solver without checking the mesh with this program'
  endif

  print *
  print *,'This program will produce histograms of mesh quality.'
  print *
  print *,'1 = also output elements above a certain skewness threshold in OpenDX format in addition to histograms'
  print *,'2 = also output elements above or below a certain element size in OpenDX format in addition to histograms'
  print *,'3 = do not output any OpenDX file, only create histograms'
  print *
  print *,'enter value:'
  read(5,*) iformat

  if (iformat < 1 .or. iformat > 3) stop 'input error, exiting...'

  if (iformat /= 3) then

    if (iformat == 1) then

      ! ask the user to imput the range of skewness to use to select the elements
      print *,'enter skewness threshold (between 0. and 0.99) above which all elements will be displayed:'
      print *,'   (beware, entering 0. will display the whole mesh, since all elements will be above that)'
      read(5,*) threshold_AVS_DX_min
      if (threshold_AVS_DX_min < 0.d0) threshold_AVS_DX_min = 0.d0
      if (threshold_AVS_DX_min > 0.99999d0) threshold_AVS_DX_min = 0.99999d0

      !!!!!!!!  print *,'enter maximum skewness threshold (between 0. and 1.):'
      !!!!!!!!!!!!!  read(5,*) threshold_AVS_DX_max
      threshold_AVS_DX_max = 0.99999d0 ! we impose to display all elements above the threshold instead
      if (threshold_AVS_DX_max < 0.d0) threshold_AVS_DX_max = 0.d0
      if (threshold_AVS_DX_max > 0.99999d0) threshold_AVS_DX_max = 0.99999d0

    else if (iformat == 2) then

      print *
      print *,'1 = output elements ABOVE a certain element size'
      print *,'2 = output elements BELOW a certain element size'
      print *
      print *,'enter value:'
      read(5,*) iabove_or_below

      if (iabove_or_below < 1 .or. iabove_or_below > 2) stop 'input error, exiting...'

      ! ask the user to imput the range of size to use to select the elements
      print *,'enter the threshold element size to use:'
      read(5,*) threshold_AVS_DX_size_to_use
      if (threshold_AVS_DX_size_to_use < ZERO) stop 'input error, exiting...'

    else
      stop 'error: incorrect value to use was entered'
    endif

  endif

! read the mesh
  print *
  print *,'start reading the node coordinate file: ',nodes_coords_file(1:len_trim(nodes_coords_file))

  open(unit=10,file=nodes_coords_file,status='old',action='read')

  read(10,*) NGLOB
  print *,'  number of points: ',NGLOB

  allocate(x(NGLOB))
  allocate(y(NGLOB))
  allocate(z(NGLOB))

  x(:) = 0.d0
  y(:) = 0.d0
  z(:) = 0.d0

  do i = 1,NGLOB

    ! gets node ID and position
    read(10,*,iostat=ier) iread,xtmp,ytmp,ztmp

    ! check
    if (ier /= 0) then
      print *,'error point read:',i,iread,xtmp,ytmp,ztmp
      stop 'error while reading points'
    endif

    ! checks if out-of-range
    if (iread < 1 .or. iread > NGLOB) then
      print *,'error at i,iread = ',i,iread
      stop 'wrong ID input for a point'
    endif

    ! stores locations
    x(iread) = xtmp
    y(iread) = ytmp
    z(iread) = ztmp

  enddo

  close(10)

  print *,'xmin, xmax of mesh read = ',minval(x),maxval(x)
  print *,'ymin, ymax of mesh read = ',minval(y),maxval(y)
  print *,'zmin, zmax of mesh read = ',minval(z),maxval(z)
  print *

  print *,'start reading the mesh topology file: ',mesh_file(1:len_trim(mesh_file))

  open(unit=10,file=mesh_file,status='old',action='read')

  read(10,*) NSPEC
  print *,'  number of elements: ',NSPEC

  allocate(ibool(NGNOD,NSPEC))
  ibool(:,:) = 0

  do i = 1,NSPEC

      ! gets element connection nodes
      read(10,*,iostat=ier) iread,n1,n2,n3,n4,n5,n6,n7,n8

      ! check
      if (ier /= 0) then
        print *,'error element read:',i,iread,n1,n2,n3,n4,n5,n6,n7,n8
        stop 'error while reading element connectivity'
      endif

      if (iread < 1 .or. iread > NSPEC .or. min(n1,n2,n3,n4,n5,n6,n7,n8) < 1 .or. max(n1,n2,n3,n4,n5,n6,n7,n8) > NGLOB) then
        print *,'error at i,iread = ',i,iread
        stop 'wrong input ID for an element'
      endif

      ! stores element nodes
      ibool(1,iread) = n1
      ibool(2,iread) = n2
      ibool(3,iread) = n3
      ibool(4,iread) = n4
      ibool(5,iread) = n5
      ibool(6,iread) = n6
      ibool(7,iread) = n7
      ibool(8,iread) = n8

  enddo

  close(10)

  print *
  print *,'start computing the minimum and maximum edge size'

! ************* compute min and max of skewness and ratios ******************

! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL
  distance_min = + HUGEVAL
  min_of_distmean = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL
  distance_max = - HUGEVAL
  max_of_distmean = - HUGEVAL

  distance_mean = ZERO

  ispec_min_edge_length = -1
  ispec_max_edge_length = -1

! loop on all the elements
  do ispec = 1,NSPEC

    if (mod(ispec,100000) == 0) print *,'processed ',ispec,' elements out of ',NSPEC

    call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

! store element number in which the edge of minimum or maximum length is located
    if (distmin < distance_min) ispec_min_edge_length = ispec
    if (distmax > distance_max) ispec_max_edge_length = ispec

! compute minimum and maximum of quality numbers
    equiangle_skewness_min = min(equiangle_skewness_min,equiangle_skewness)
    edge_aspect_ratio_min = min(edge_aspect_ratio_min,edge_aspect_ratio)
    diagonal_aspect_ratio_min = min(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
    distance_min = min(distance_min,distmin)

    equiangle_skewness_max = max(equiangle_skewness_max,equiangle_skewness)
    edge_aspect_ratio_max = max(edge_aspect_ratio_max,edge_aspect_ratio)
    diagonal_aspect_ratio_max = max(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
    distance_max = max(distance_max,distmax)

    distance_mean = distance_mean + distmean
    min_of_distmean = min(min_of_distmean,distmean)
    max_of_distmean = max(max_of_distmean,distmean)

  enddo

! compute the mean distance
  distance_mean = distance_mean / dble(NSPEC)

  print *,'done processing ',NSPEC,' elements out of ',NSPEC

  DISPLAY_HISTOGRAM_DISTMEAN = .true.
  if (abs(max_of_distmean - min_of_distmean) < 1.d-3*distmean) then
    print *
    print *,'Your input mesh seems to be perfect, i.e. all mean distances are equal;'
    print *,'Will thus not display any histogram of mean distance in the mesh, since it would lead to division by zero.'
    print *
    DISPLAY_HISTOGRAM_DISTMEAN = .false.
  endif

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions:'
  print *
  print *,'equiangle skewness: 0. perfect,  1. bad'
  print *,'skewness max deviation angle: 0. perfect,  90. bad'
  print *,'edge aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'------------'

  print *
  print *,'minimum length of an edge in the whole mesh (m) = ',distance_min,' in element ',ispec_min_edge_length
  print *,'mean length of an edge in the whole mesh (m) = ',distance_mean
  print *,'maximum length of an edge in the whole mesh (m) = ',distance_max,' in element ',ispec_max_edge_length
  print *
  print *,'max equiangle skewness = ',equiangle_skewness_max
! print *,'min equiangle skewness = ',equiangle_skewness_min
  print *
  print *,'max deviation angle from a right angle (90 degrees) is therefore = ',90.*equiangle_skewness_max
  print *
  print *,'worst angle in the mesh is therefore ',90.*(1. - equiangle_skewness_max)
  print *,'or ',180. - 90.*(1. - equiangle_skewness_max),' degrees'
  print *
  print *,'max edge aspect ratio = ',edge_aspect_ratio_max
! print *,'min edge aspect ratio = ',edge_aspect_ratio_min
  print *
  print *,'max diagonal aspect ratio = ',diagonal_aspect_ratio_max
! print *,'min diagonal aspect ratio = ',diagonal_aspect_ratio_min
  print *

  if (iformat == 2) then
    if (iabove_or_below == 1) then
! output elements ABOVE a certain element size
      threshold_AVS_DX_min = threshold_AVS_DX_size_to_use
      threshold_AVS_DX_max = distance_max
    else
! output elements BELOW a certain element size
      threshold_AVS_DX_min = distance_min
      threshold_AVS_DX_max = threshold_AVS_DX_size_to_use
    endif
  endif

!---------------------------------------------------------------

! create statistics about mesh quality

  print *
  print *,'creating histogram of mesh quality'
  print *

! erase histograms
  classes_of_histogram_skewness(:) = 0
  classes_of_histogram_meansize(:) = 0

! loop on all the elements
  do ispec = 1,NSPEC
    call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

! store skewness in histogram
    iclass = int(equiangle_skewness * dble(NCLASS))
    if (iclass < 0) iclass = 0
    if (iclass > NCLASS-1) iclass = NCLASS-1
    classes_of_histogram_skewness(iclass) = classes_of_histogram_skewness(iclass) + 1

! store mean size in histogram
    if (DISPLAY_HISTOGRAM_DISTMEAN) then
      iclass = int(((distmean - min_of_distmean) / (max_of_distmean - min_of_distmean)) * dble(NCLASS))
      if (iclass < 0) iclass = 0
      if (iclass > NCLASS-1) iclass = NCLASS-1
      classes_of_histogram_meansize(iclass) = classes_of_histogram_meansize(iclass) + 1
    endif

  enddo

!---------------------------------------------------------------

! create histogram of skewness and save it in a Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file='OUTPUT_FILES/mesh_quality_histogram_skewness.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_of_histogram_skewness(iclass))/dble(NSPEC)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_of_histogram_skewness(iclass),' ', &
               sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

! display warning if maximum skewness is too high
  if (equiangle_skewness_max > 0.85d0) then
    print *
    print *,'********************************************'
    print *,'********************************************'
    print *,' WARNING, mesh is bad (max skewness > 0.85)'
    print *,'********************************************'
    print *,'********************************************'
    print *
  else if (equiangle_skewness_max > 0.75d0) then
    print *
    print *,'******************************************************************'
    print *,'******************************************************************'
    print *,' WARNING, mesh is maybe not fully optimized (max skewness > 0.75)'
    print *,'******************************************************************'
    print *,'******************************************************************'
    print *
  endif

  if (total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    print *,'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  endif

!---------------------------------------------------------------

  if (DISPLAY_HISTOGRAM_DISTMEAN) then

! create histogram of mean distance and save it in a Gnuplot file
  print *
  print *,'histogram of mean element size:'
  print *
  total_percent = 0.
  open(unit=14,file='OUTPUT_FILES/mesh_quality_histogram_meansize.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_of_histogram_meansize(iclass))/dble(NSPEC)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS) * (max_of_distmean - min_of_distmean) + min_of_distmean),' - ', &
            real((iclass+1)/dble(NCLASS) * (max_of_distmean - min_of_distmean) + min_of_distmean), &
            classes_of_histogram_meansize(iclass),' ',sngl(current_percent),' %'
    write(14,*) sngl(0.5d0*(dble(iclass/dble(NCLASS)) + dble((iclass+1)/dble(NCLASS))) &
                 * (max_of_distmean - min_of_distmean) + min_of_distmean),' ',sngl(current_percent)
  enddo
  close(14)

  if (total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    print *,'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  endif

  endif

!---------------------------------------------------------------

! create script for Gnuplot histogram files
  open(unit=14,file='OUTPUT_FILES/plot_mesh_quality_histograms.gnu',status='unknown')
  write(14,*) 'set term wxt'
  write(14,*) '#set term gif'

  if (DISPLAY_HISTOGRAM_DISTMEAN) then
    write(14,*) '#set output "mesh_quality_histogram_meansize.gif"'
    write(14,*) 'set xrange [',sngl(min_of_distmean),':',sngl(max_of_distmean),']'
    write(14,*) 'set boxwidth ',sngl((max_of_distmean - min_of_distmean)/dble(NCLASS))
    write(14,*) 'set xlabel "Range of mean size of element"'
    write(14,*) 'set ylabel "Percentage of elements (%)"'
    write(14,*) 'plot "mesh_quality_histogram_meansize.txt" with boxes'
    write(14,*) 'pause -1 "hit any key..."'
  endif

  write(14,*) '#set output "mesh_quality_histogram_skewness.gif"'
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',sngl(1.d0/dble(NCLASS))
  write(14,*) 'set xlabel "Range of skewness"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram_skewness.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'

  close(14)

!---------------------------------------------------------------

! ************* create OpenDX file with elements in a certain selected range

  if (iformat /= 3) then

  print *
  print *,'creating OpenDX file with subset of elements in selected range'
  print *,'between ',threshold_AVS_DX_min,' and ',threshold_AVS_DX_max
  print *

  if (threshold_AVS_DX_min > threshold_AVS_DX_max) stop 'error: incorrect display range selected'

! ************* count number of elements in selected range *************

! erase number of elements belonging to selected range for AVS_DX
  ntotspecAVS_DX = 0

! loop on all the elements
  do ispec = 1,NSPEC

    call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

    if (iformat == 1) then
      value_to_use = equiangle_skewness
    else if (iformat == 2) then
      value_to_use = distmean
    else
      stop 'error: incorrect value to use was entered'
    endif

! check if element belongs to selected range
    if (value_to_use >= threshold_AVS_DX_min .and. value_to_use <= threshold_AVS_DX_max) ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  if (ntotspecAVS_DX == 0) then
    stop 'no elements in selected range, no file created'
  else
    print *
    print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX selected range ',threshold_AVS_DX_min,threshold_AVS_DX_max
    print *
  endif

  open(unit=11,file='OUTPUT_FILES/DX_mesh_quality.dx',status='unknown')

! ************* generate points ******************

! write OpenDX header
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NGLOB,' data follows'

! write all the points
  do i = 1,NGLOB
    write(11,*) sngl(x(i)),sngl(y(i)),sngl(z(i))
  enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape ',NGNOD,' items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  ispec_begin = 1
  ispec_end = NSPEC

  do ispec = ispec_begin,ispec_end

    call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

    if (iformat == 1) then
      value_to_use = equiangle_skewness
    else if (iformat == 2) then
      value_to_use = distmean
    else
      stop 'error: incorrect value to use was entered'
    endif

! check if element needs to be output
    if (value_to_use >= threshold_AVS_DX_min .and. value_to_use <= threshold_AVS_DX_max) then
! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
      write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool(4,ispec)-1, ibool(1,ispec)-1, ibool(8,ispec)-1, ibool(5,ispec)-1, &
            ibool(3,ispec)-1, ibool(2,ispec)-1, ibool(7,ispec)-1, ibool(6,ispec)-1
!     print *,'element ',ispec,' belongs to the range'
    endif

  enddo

! ************* generate element data values ******************

! output OpenDX header for data
  write(11,*) 'attribute "element type" string "cubes"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  do ispec = ispec_begin,ispec_end

    call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

    if (iformat == 1) then
      value_to_use = equiangle_skewness
    else if (iformat == 2) then
      value_to_use = distmean
    else
      stop 'error: incorrect value to use was entered'
    endif

! check if element needs to be output
    if (value_to_use >= threshold_AVS_DX_min .and. value_to_use <= threshold_AVS_DX_max) write(11,*) sngl(value_to_use)

  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

! close OpenDX file
  close(11)

  endif

  print *

  end program check_mesh_quality

!
!=====================================================================
!

! create mesh quality data for a given 3D spectral element

  subroutine create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NGLOB,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax,distmean)

  use constants

  implicit none

  integer, parameter :: NGNOD = 8                        ! hexahedral elements

  integer :: iface,icorner,ispec,NSPEC,NGLOB,i

  double precision, dimension(NGLOB) :: x,y,z

  integer, dimension(NGNOD,NSPEC) :: ibool

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  double precision :: argument_of_arccos
  double precision :: vectorA_x,vectorA_y,vectorA_z
  double precision :: vectorB_x,vectorB_y,vectorB_z
  double precision :: norm_A,norm_B,angle_vectors
  double precision :: distmin,distmax,distmean,dist,dist1,dist2,dist3,dist4
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio

  integer :: count_contributions

! for stability
  double precision :: stability,VP_MAX,delta_t

! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision, dimension(NGLL_MAX_STABILITY) :: percent_GLL

! topology of faces of cube for skewness
  integer faces_topo(6,6)

! store the corners of this element for the skewness routine
  do i = 1,NGNOD
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
  if (NGLLX > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

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

! compute equiangle skewness (as defined in Fluent/Gambit manual)
! and compute edge aspect ratio using the corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     distmean = ZERO
     count_contributions = 0
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

! sanity check
         if (norm_A <= ZERO .or. norm_B <= ZERO) then
           print *,'error detected in element ',ispec,' out of ',NSPEC
           print *,'error: negative of null norm found, norm_A, norm_B = ',norm_A, norm_B
           stop 'error in the norm found'
         endif

! angle formed by the two vectors
         argument_of_arccos = (vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B)

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

         count_contributions = count_contributions + 1
         distmean = distmean + dist

       enddo
     enddo

! compute the mean distance
   distmean = distmean / count_contributions

! compute edge aspect ratio
   edge_aspect_ratio = distmax / distmin

   stability = delta_t * VP_MAX / (distmin * percent_GLL(NGLLX))

! compute diagonal aspect ratio
   dist1 = sqrt((xelm(1) - xelm(7))**2 + (yelm(1) - yelm(7))**2 + (zelm(1) - zelm(7))**2)
   dist2 = sqrt((xelm(2) - xelm(8))**2 + (yelm(2) - yelm(8))**2 + (zelm(2) - zelm(8))**2)
   dist3 = sqrt((xelm(3) - xelm(5))**2 + (yelm(3) - yelm(5))**2 + (zelm(3) - zelm(5))**2)
   dist4 = sqrt((xelm(4) - xelm(6))**2 + (yelm(4) - yelm(6))**2 + (zelm(4) - zelm(6))**2)
   diagonal_aspect_ratio = max(dist1,dist2,dist3,dist4) / min(dist1,dist2,dist3,dist4)

  end subroutine create_mesh_quality_data_3D

