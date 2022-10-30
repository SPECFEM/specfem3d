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

! read an external mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Pau, France, March 2009 and CNRS, Marseille, France, 2015, 2016, 2017.

  program check_mesh_quality

  use constants, only: NDIM,HUGEVAL,ZERO

  implicit none

!------------------------------------------------------------------------------------------------

  character(len=*), parameter :: nodes_coords_file = 'MESH/nodes_coords_file'
  character(len=*), parameter :: mesh_file = 'MESH/mesh_file'

!------------------------------------------------------------------------------------------------

  integer :: NGLOB                    ! number of nodes
  integer :: NSPEC                    ! number of elements
  integer :: NGNOD                    ! number of control nodes for hexahedral elements (can only be 8 or 27)

  double precision, dimension(:), allocatable :: x,y,z

  integer, dimension(:,:), allocatable :: ibool

  integer :: i,ia,ispec,iread,iformat,iabove_or_below,ispec_min_edge_length,ispec_max_edge_length, &
             ispec_begin,ispec_end,ier,itype_of_hex

  double precision :: xtmp,ytmp,ztmp

! for quality of mesh
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,value_to_use
  double precision :: equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision :: equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision :: threshold_AVS_DX_min,threshold_AVS_DX_max,threshold_AVS_DX_size_to_use
  double precision :: distance_min,distance_max,distance_mean
  double precision :: distmin,distmax,distmean,min_of_distmean,max_of_distmean

! for histogram
  integer, parameter :: NCLASS = 40
  integer, dimension(0:NCLASS-1) :: classes_of_histogram_skewness,classes_of_histogram_meansize
  integer :: iclass
  double precision :: current_percent,total_percent

! to export elements that have a certain range to OpenDX
  integer :: ntotspecAVS_DX
  logical :: DISPLAY_HISTOGRAM_DISTMEAN

! Gauss-Lobatto-Legendre points of integration, to check for negative Jacobians
! we are lazy here and hardwire a fixed value of NGLL = 5 even if someone changes NGLL in setup/constants.h.in
! in order to be able to hardwire the position of the GLL points for NGLL = 5 below in this routine
! instead of having to link with a bunch of libraries to compute them; this works fine because
! we will only use the element corners in this program, i.e. the GLL points that are always -1 and +1,
! not the others, and thus we do not care. However that is a bit ugly and we should probably
! update the Makefile one day to just link with the GLL libraries of SPECFEM...
! this would not change the behavior of this program though...
  integer, parameter :: local_NGLLX_always_5 = 5,local_NGLLY_always_5 = 5,local_NGLLZ_always_5 = 5
  double precision xigll(local_NGLLX_always_5)
  double precision yigll(local_NGLLY_always_5)
  double precision zigll(local_NGLLZ_always_5)

! 3D shape function derivatives, to check for negative Jacobians
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision, dimension(:), allocatable :: xelm,yelm,zelm

  double precision :: jacobian
  logical :: found_a_negative_jacobian

  print *
  print *,'This program will produce histograms of mesh quality.'
  print *

  print *,'1 = the mesh contains HEX8 elements'
  print *,'2 = the mesh contains HEX27 elements'
  print *,'3 = exit'
  read(*,*) itype_of_hex
  if (itype_of_hex /= 1 .and. itype_of_hex /= 2) stop 'exiting...'
  if (itype_of_hex == 1) then
    NGNOD = 8
  else
    NGNOD = 27
  endif
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

! hardwire GLL point location values to avoid having to link with a long library to compute them
  xigll(:) = (/ -1.d0 , -0.654653670707977d0 , 0.d0 , 0.654653670707977d0 , 1.d0 /)
  yigll(:) = xigll(:)
  zigll(:) = xigll(:)

  allocate(dershape3D(NDIM,NGNOD,local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1094')

! compute the derivatives of the 3D shape functions for a 8-node or 27-node element
  call local_version_of_get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD, &
                                 local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

  allocate(xelm(NGNOD),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1095')
  allocate(yelm(NGNOD),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1096')
  allocate(zelm(NGNOD),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1097')

  if (NGNOD == 8) then
    print *
    print *,'reading HEX8 input mesh file...'
  else
    print *
    print *,'reading HEX27 input mesh file...'
  endif

! read the mesh
  print *
  print *,'start reading the node coordinate file: ',nodes_coords_file(1:len_trim(nodes_coords_file))

  open(unit=10,file=nodes_coords_file,status='old',action='read')

  read(10,*) NGLOB
  print *,'  number of points: ',NGLOB

  allocate(x(NGLOB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1098')
  allocate(y(NGLOB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1099')
  allocate(z(NGLOB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1100')

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

  allocate(ibool(NGNOD,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1101')
  ibool(:,:) = 0

  do i = 1,NSPEC

      ! gets element connection nodes
      read(10,*,iostat=ier) iread,(ibool(ia,iread),ia=1,NGNOD)

      ! check
      if (ier /= 0) then
        print *,'error element read:',i,iread
        stop 'error while reading element connectivity'
      endif

      if (iread < 1 .or. iread > NSPEC .or. minval(ibool(:,iread)) < 1 .or. maxval(ibool(:,iread)) > NGLOB) then
        print *,'error at i,iread = ',i,iread
        stop 'wrong input ID for an element'
      endif

  enddo

  close(10)

! check the mesh read to make sure it contains no negative Jacobians
  do ispec = 1,nspec
! check the element for a negative Jacobian
      do ia = 1,NGNOD
        xelm(ia) = x(ibool(ia,ispec))
        yelm(ia) = y(ibool(ia,ispec))
        zelm(ia) = z(ibool(ia,ispec))
      enddo

      call local_version_of_calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian, &
                                          NDIM,NGNOD, &
                                          local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5,jacobian)

      if (found_a_negative_jacobian) then
        print *,'detected an element with negative Jacobian in the input mesh: element ',ispec, &
                     ' in which the Jacobian is ',jacobian
        stop 'error: the mesh read contains a negative Jacobian!'
      endif
  enddo
  print *
  print *,'input mesh successfully checked for negative Jacobians, it contains none'

! ************* compute min and max of skewness and ratios ******************

  print *
  print *,'start computing the minimum and maximum edge size'

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

    call local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

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
    call local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

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
  write(14,*) 'set terminal x11'
  write(14,*) '#set terminal wxt'
  write(14,*) '#set terminal gif'

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

    call local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

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

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  ispec_begin = 1
  ispec_end = NSPEC

  do ispec = ispec_begin,ispec_end

    call local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

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

    call local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

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

  subroutine local_version_of_create_mesh_quality_data_3D(x,y,z,ibool,ispec,NGNOD,NSPEC,NGLOB, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax,distmean)

  use constants

  implicit none

  integer :: NGNOD,iface,icorner,ispec,NSPEC,NGLOB,i

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

! topology of faces of cube for skewness
  integer faces_topo(6,6)

! store the corners of this element for the skewness routine
  do i = 1,NGNOD
    xelm(i) = x(ibool(i,ispec))
    yelm(i) = y(ibool(i,ispec))
    zelm(i) = z(ibool(i,ispec))
  enddo

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

! compute diagonal aspect ratio
   dist1 = sqrt((xelm(1) - xelm(7))**2 + (yelm(1) - yelm(7))**2 + (zelm(1) - zelm(7))**2)
   dist2 = sqrt((xelm(2) - xelm(8))**2 + (yelm(2) - yelm(8))**2 + (zelm(2) - zelm(8))**2)
   dist3 = sqrt((xelm(3) - xelm(5))**2 + (yelm(3) - yelm(5))**2 + (zelm(3) - zelm(5))**2)
   dist4 = sqrt((xelm(4) - xelm(6))**2 + (yelm(4) - yelm(6))**2 + (zelm(4) - zelm(6))**2)
   diagonal_aspect_ratio = max(dist1,dist2,dist3,dist4) / min(dist1,dist2,dist3,dist4)

  end subroutine local_version_of_create_mesh_quality_data_3D

!
!=====================================================================
!

! 3D shape functions for 8-node or 27-node element

  subroutine local_version_of_get_shape3D(dershape3D,xigll,yigll,zigll,NGNOD, &
                                 local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

  use constants

  implicit none

  integer NGNOD

  integer :: local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(local_NGLLX_always_5)
  double precision yigll(local_NGLLY_always_5)
  double precision zigll(local_NGLLZ_always_5)

! 3D shape function derivatives
  double precision dershape3D(NDIM,NGNOD,local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

  integer i,j,k,ia

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

! for checking the 3D shape functions
  double precision sumdershapexi,sumdershapeeta,sumdershapegamma

  double precision, parameter :: ONE_EIGHTH = 0.125d0

! check that the parameter file is correct
  if (NGNOD /= 8 .and. NGNOD /= 27) stop 'volume elements should have 8 or 27 control nodes'

! ***
! *** create 3D shape functions and jacobian
! ***

  do i=1,local_NGLLX_always_5
    do j=1,local_NGLLY_always_5
      do k=1,local_NGLLZ_always_5

        xi = xigll(i)
        eta = yigll(j)
        gamma = zigll(k)

        !--- case of a 3D 8-node element (Dhatt-Touzot p. 115)
        if (NGNOD == 8) then

          ra1 = one + xi
          ra2 = one - xi

          rb1 = one + eta
          rb2 = one - eta

          rc1 = one + gamma
          rc2 = one - gamma

          dershape3D(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
          dershape3D(1,2,i,j,k) = ONE_EIGHTH*rb2*rc2
          dershape3D(1,3,i,j,k) = ONE_EIGHTH*rb1*rc2
          dershape3D(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
          dershape3D(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
          dershape3D(1,6,i,j,k) = ONE_EIGHTH*rb2*rc1
          dershape3D(1,7,i,j,k) = ONE_EIGHTH*rb1*rc1
          dershape3D(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1

          dershape3D(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
          dershape3D(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
          dershape3D(2,3,i,j,k) = ONE_EIGHTH*ra1*rc2
          dershape3D(2,4,i,j,k) = ONE_EIGHTH*ra2*rc2
          dershape3D(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
          dershape3D(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
          dershape3D(2,7,i,j,k) = ONE_EIGHTH*ra1*rc1
          dershape3D(2,8,i,j,k) = ONE_EIGHTH*ra2*rc1

          dershape3D(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
          dershape3D(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
          dershape3D(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
          dershape3D(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
          dershape3D(3,5,i,j,k) = ONE_EIGHTH*ra2*rb2
          dershape3D(3,6,i,j,k) = ONE_EIGHTH*ra1*rb2
          dershape3D(3,7,i,j,k) = ONE_EIGHTH*ra1*rb1
          dershape3D(3,8,i,j,k) = ONE_EIGHTH*ra2*rb1

        else

          ! note: put further initialization for NGNOD == 27 into subroutine
          !       to avoid compilation errors in case NGNOD == 8
          call local_version_of_get_shape3D_27(NGNOD,dershape3D,xi,eta,gamma,i,j,k, &
                                 local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

        endif

      enddo
    enddo
  enddo

!--- check the shape functions and their derivatives

  do i=1,local_NGLLX_always_5
    do j=1,local_NGLLY_always_5
      do k=1,local_NGLLZ_always_5

        sumdershapexi = ZERO
        sumdershapeeta = ZERO
        sumdershapegamma = ZERO

        do ia=1,NGNOD
          sumdershapexi = sumdershapexi + dershape3D(1,ia,i,j,k)
          sumdershapeeta = sumdershapeeta + dershape3D(2,ia,i,j,k)
          sumdershapegamma = sumdershapegamma + dershape3D(3,ia,i,j,k)
        enddo

        ! sum of shape functions should be one
        ! sum of derivative of shape functions should be zero
        if (abs(sumdershapexi) > TINYVAL) stop 'error in xi derivative of 3D shape functions'
        if (abs(sumdershapeeta) > TINYVAL) stop 'error in eta derivative of 3D shape functions'
        if (abs(sumdershapegamma) > TINYVAL) stop 'error in gamma derivative of 3D shape functions'

      enddo
    enddo
  enddo

  end subroutine local_version_of_get_shape3D

!
!-------------------------------------------------------------------------------------------------
!

!--- case of a 3D 27-node element

  subroutine local_version_of_get_shape3D_27(NGNOD,dershape3D,xi,eta,gamma,i,j,k, &
                                 local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

  use constants

  implicit none

  integer :: NGNOD,i,j,k

  integer :: local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5

! 3D shape function derivatives
  double precision dershape3D(NDIM,NGNOD,local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

! location of the nodes of the 3D hexahedra elements
  double precision xi,eta,gamma
  double precision l1xi,l2xi,l3xi,l1eta,l2eta,l3eta,l1gamma,l2gamma,l3gamma
  double precision l1pxi,l2pxi,l3pxi,l1peta,l2peta,l3peta,l1pgamma,l2pgamma,l3pgamma

  l1xi=HALF*xi*(xi-ONE)
  l2xi=ONE-xi**2
  l3xi=HALF*xi*(xi+ONE)

  l1pxi=xi-HALF
  l2pxi=-TWO*xi
  l3pxi=xi+HALF

  l1eta=HALF*eta*(eta-ONE)
  l2eta=ONE-eta**2
  l3eta=HALF*eta*(eta+ONE)

  l1peta=eta-HALF
  l2peta=-TWO*eta
  l3peta=eta+HALF

  l1gamma=HALF*gamma*(gamma-ONE)
  l2gamma=ONE-gamma**2
  l3gamma=HALF*gamma*(gamma+ONE)

  l1pgamma=gamma-HALF
  l2pgamma=-TWO*gamma
  l3pgamma=gamma+HALF

  ! corner nodes

  dershape3D(1,1,i,j,k)=l1pxi*l1eta*l1gamma
  dershape3D(1,2,i,j,k)=l3pxi*l1eta*l1gamma
  dershape3D(1,3,i,j,k)=l3pxi*l3eta*l1gamma
  dershape3D(1,4,i,j,k)=l1pxi*l3eta*l1gamma
  dershape3D(1,5,i,j,k)=l1pxi*l1eta*l3gamma
  dershape3D(1,6,i,j,k)=l3pxi*l1eta*l3gamma
  dershape3D(1,7,i,j,k)=l3pxi*l3eta*l3gamma
  dershape3D(1,8,i,j,k)=l1pxi*l3eta*l3gamma

  dershape3D(2,1,i,j,k)=l1xi*l1peta*l1gamma
  dershape3D(2,2,i,j,k)=l3xi*l1peta*l1gamma
  dershape3D(2,3,i,j,k)=l3xi*l3peta*l1gamma
  dershape3D(2,4,i,j,k)=l1xi*l3peta*l1gamma
  dershape3D(2,5,i,j,k)=l1xi*l1peta*l3gamma
  dershape3D(2,6,i,j,k)=l3xi*l1peta*l3gamma
  dershape3D(2,7,i,j,k)=l3xi*l3peta*l3gamma
  dershape3D(2,8,i,j,k)=l1xi*l3peta*l3gamma

  dershape3D(3,1,i,j,k)=l1xi*l1eta*l1pgamma
  dershape3D(3,2,i,j,k)=l3xi*l1eta*l1pgamma
  dershape3D(3,3,i,j,k)=l3xi*l3eta*l1pgamma
  dershape3D(3,4,i,j,k)=l1xi*l3eta*l1pgamma
  dershape3D(3,5,i,j,k)=l1xi*l1eta*l3pgamma
  dershape3D(3,6,i,j,k)=l3xi*l1eta*l3pgamma
  dershape3D(3,7,i,j,k)=l3xi*l3eta*l3pgamma
  dershape3D(3,8,i,j,k)=l1xi*l3eta*l3pgamma

  ! midside nodes

  dershape3D(1,9,i,j,k)=l2pxi*l1eta*l1gamma
  dershape3D(1,10,i,j,k)=l3pxi*l2eta*l1gamma
  dershape3D(1,11,i,j,k)=l2pxi*l3eta*l1gamma
  dershape3D(1,12,i,j,k)=l1pxi*l2eta*l1gamma
  dershape3D(1,13,i,j,k)=l1pxi*l1eta*l2gamma
  dershape3D(1,14,i,j,k)=l3pxi*l1eta*l2gamma
  dershape3D(1,15,i,j,k)=l3pxi*l3eta*l2gamma
  dershape3D(1,16,i,j,k)=l1pxi*l3eta*l2gamma
  dershape3D(1,17,i,j,k)=l2pxi*l1eta*l3gamma
  dershape3D(1,18,i,j,k)=l3pxi*l2eta*l3gamma
  dershape3D(1,19,i,j,k)=l2pxi*l3eta*l3gamma
  dershape3D(1,20,i,j,k)=l1pxi*l2eta*l3gamma

  dershape3D(2,9,i,j,k)=l2xi*l1peta*l1gamma
  dershape3D(2,10,i,j,k)=l3xi*l2peta*l1gamma
  dershape3D(2,11,i,j,k)=l2xi*l3peta*l1gamma
  dershape3D(2,12,i,j,k)=l1xi*l2peta*l1gamma
  dershape3D(2,13,i,j,k)=l1xi*l1peta*l2gamma
  dershape3D(2,14,i,j,k)=l3xi*l1peta*l2gamma
  dershape3D(2,15,i,j,k)=l3xi*l3peta*l2gamma
  dershape3D(2,16,i,j,k)=l1xi*l3peta*l2gamma
  dershape3D(2,17,i,j,k)=l2xi*l1peta*l3gamma
  dershape3D(2,18,i,j,k)=l3xi*l2peta*l3gamma
  dershape3D(2,19,i,j,k)=l2xi*l3peta*l3gamma
  dershape3D(2,20,i,j,k)=l1xi*l2peta*l3gamma

  dershape3D(3,9,i,j,k)=l2xi*l1eta*l1pgamma
  dershape3D(3,10,i,j,k)=l3xi*l2eta*l1pgamma
  dershape3D(3,11,i,j,k)=l2xi*l3eta*l1pgamma
  dershape3D(3,12,i,j,k)=l1xi*l2eta*l1pgamma
  dershape3D(3,13,i,j,k)=l1xi*l1eta*l2pgamma
  dershape3D(3,14,i,j,k)=l3xi*l1eta*l2pgamma
  dershape3D(3,15,i,j,k)=l3xi*l3eta*l2pgamma
  dershape3D(3,16,i,j,k)=l1xi*l3eta*l2pgamma
  dershape3D(3,17,i,j,k)=l2xi*l1eta*l3pgamma
  dershape3D(3,18,i,j,k)=l3xi*l2eta*l3pgamma
  dershape3D(3,19,i,j,k)=l2xi*l3eta*l3pgamma
  dershape3D(3,20,i,j,k)=l1xi*l2eta*l3pgamma

  ! side center nodes

  dershape3D(1,21,i,j,k)=l2pxi*l2eta*l1gamma
  dershape3D(1,22,i,j,k)=l2pxi*l1eta*l2gamma
  dershape3D(1,23,i,j,k)=l3pxi*l2eta*l2gamma
  dershape3D(1,24,i,j,k)=l2pxi*l3eta*l2gamma
  dershape3D(1,25,i,j,k)=l1pxi*l2eta*l2gamma
  dershape3D(1,26,i,j,k)=l2pxi*l2eta*l3gamma

  dershape3D(2,21,i,j,k)=l2xi*l2peta*l1gamma
  dershape3D(2,22,i,j,k)=l2xi*l1peta*l2gamma
  dershape3D(2,23,i,j,k)=l3xi*l2peta*l2gamma
  dershape3D(2,24,i,j,k)=l2xi*l3peta*l2gamma
  dershape3D(2,25,i,j,k)=l1xi*l2peta*l2gamma
  dershape3D(2,26,i,j,k)=l2xi*l2peta*l3gamma

  dershape3D(3,21,i,j,k)=l2xi*l2eta*l1pgamma
  dershape3D(3,22,i,j,k)=l2xi*l1eta*l2pgamma
  dershape3D(3,23,i,j,k)=l3xi*l2eta*l2pgamma
  dershape3D(3,24,i,j,k)=l2xi*l3eta*l2pgamma
  dershape3D(3,25,i,j,k)=l1xi*l2eta*l2pgamma
  dershape3D(3,26,i,j,k)=l2xi*l2eta*l3pgamma

  ! center node

  dershape3D(1,27,i,j,k)=l2pxi*l2eta*l2gamma
  dershape3D(2,27,i,j,k)=l2xi*l2peta*l2gamma
  dershape3D(3,27,i,j,k)=l2xi*l2eta*l2pgamma

  end subroutine local_version_of_get_shape3D_27

!
!=====================================================================
!

  subroutine local_version_of_calc_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian, &
                                            NDIM,NGNOD, &
                                            local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5,jacobian)

  implicit none

  integer :: NDIM,NGNOD,local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5

  logical :: found_a_negative_jacobian

  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision dershape3D(NDIM,NGNOD,local_NGLLX_always_5,local_NGLLY_always_5,local_NGLLZ_always_5)

  integer :: i,j,k,ia
  double precision :: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision :: jacobian

  double precision, parameter :: ZERO = 0.d0

  found_a_negative_jacobian = .false.

! do k=1,local_NGLLZ_always_5
!   do j=1,local_NGLLY_always_5
!     do i=1,local_NGLLX_always_5
! for this CPML mesh extrusion routine it is sufficient to test the 8 corners of each element to reduce the cost
! because we just want to detect if the element is flipped or not, and if so flip it back
  do k=1,local_NGLLZ_always_5,local_NGLLZ_always_5-1
    do j=1,local_NGLLY_always_5,local_NGLLY_always_5-1
      do i=1,local_NGLLX_always_5,local_NGLLX_always_5-1

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)

        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)

        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
      enddo

      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
      if (jacobian <= ZERO) then
        found_a_negative_jacobian = .true.
        return
      endif

      enddo
    enddo
  enddo

  end subroutine local_version_of_calc_jacobian

