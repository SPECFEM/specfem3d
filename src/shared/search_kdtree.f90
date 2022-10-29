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

!------------------------------------------------------------------------------
!
! kdtree - nearest neighbor search
!
!------------------------------------------------------------------------------

module kdtree_search

  ! kd-tree for searching nearest neighbors

  use constants, only: IMAIN

  private

  ! single tree node
  type :: kdtree_node
    ! id (used for debugging tree)
    !integer :: id
    ! index of associated data point
    integer :: ipoint
    ! index of lower and upper bounds of point range
    integer :: ibound_lower,ibound_upper

    ! cut dimension index ( 1/2/3 for 3D point set )
    integer :: idim
    ! cut value
    double precision :: cut_value
    ! bounds value on cut dimension
    !double precision :: cut_min,cut_max

    ! child nodes in sub level
    type (kdtree_node), pointer :: left, right
  end type kdtree_node

  ! data associated kd-tree structure (root node)
  ! note: in general, a save attribute should be added to this root pointer
  !       to have the tree variable stored globally
  !       however, putting it here into the module declaration will have the same effect
  type (kdtree_node), pointer :: kdtree

  ! kdtree arrays
  ! total number of tree nodes
  integer :: kdtree_num_nodes = 0
  ! tree node locations
  double precision, dimension(:,:),allocatable,target :: kdtree_nodes_location
  ! associated node index
  integer, dimension(:),allocatable :: kdtree_nodes_index

  ! n-search result arrays
  integer :: kdtree_search_num_nodes = 0
  ! n-search result node index
  integer, dimension(:),allocatable :: kdtree_search_index
  ! sorts result node index array (increasing order)
  logical, parameter :: DO_SORT_RESULTS = .true.

  ! info output
  logical :: be_verbose = .false.

  !---------------------------------------------------------------
  ! public routines
  public :: kdtree_setup
  public :: kdtree_find_nearest_neighbor
  public :: kdtree_find_nearest_n_neighbors
  public :: kdtree_count_nearest_n_neighbors
  public :: kdtree_get_nearest_n_neighbors
  public :: kdtree_count_nearest_n_neighbors_ellip
  public :: kdtree_get_nearest_n_neighbors_ellip
  public :: kdtree_delete
  public :: kdtree_set_verbose

  ! public parameters/arrays
  public :: kdtree_num_nodes
  public :: kdtree_nodes_location
  public :: kdtree_nodes_index
  public :: kdtree_search_num_nodes
  public :: kdtree_search_index
  !---------------------------------------------------------------

contains

! example:
!
! creates kd-tree for searching
!  .. prepare point array kdtree_nodes_location,kdtree_num_nodes
!  call kdtree_setup()
!
! finds closest point
!  .. do-loop points at xyz
!     call kdtree_find_nearest_neighbor(xyz,iglob_min,dist_min)
!  .. enddo
!
! deletes search tree
!  call kdtree_delete()


  subroutine kdtree_setup()

  ! sets up the kd-tree structure
  !
  ! needs:
  !   kdtree_num_nodes  - total number of points
  !   kdtree_nodes_location   - 3D array of points
  !
  ! returns:
  !   creates internal tree representation
  !
  use constants, only: myrank

  implicit none

  ! local parameters
  integer :: npoints
  integer, dimension(:), allocatable :: points_index
  double precision,dimension(:,:), pointer :: points_data

  ! tree statistics
  integer :: depth
  integer :: numnodes,maxdepth
  integer :: i,ier

  ! timing
  real :: ct_start,ct_end

  ! test search
  double precision, dimension(3) :: xyz_target
  integer :: ipoint_min
  double precision :: dist_min

  !------------------------------------------------------

  ! debugging: performs a test search
  logical,parameter :: DEBUG = .false.

  !------------------------------------------------------

  ! checks
  if (kdtree_num_nodes <= 0 ) stop 'Error creating kdtree with zero nodes is invalid'
  if (.not. allocated(kdtree_nodes_location) ) stop 'Error array kdtree_nodes_location not allocated yet'

  ! timing
  call cpu_time(ct_start)

  ! number of data points
  npoints = kdtree_num_nodes

  ! 3D point coordinates
  points_data => kdtree_nodes_location(:,:)

  if (be_verbose) then
    write(IMAIN,*)
    write(IMAIN,*) 'kd-tree:'
    write(IMAIN,*) '  total data points: ',npoints
    !write(IMAIN,*) '  box boundaries   : x min/max = ',minval(points_data(1,:)),maxval(points_data(1,:))
    !write(IMAIN,*) '                     y min/max = ',minval(points_data(2,:)),maxval(points_data(2,:))
    !write(IMAIN,*) '                     z min/max = ',minval(points_data(3,:)),maxval(points_data(3,:))
    call flush_IMAIN()
  endif

  ! theoretical number of node for totally balanced tree
  numnodes = npoints
  i = npoints
  do while ( i >= 1 )
    i = i / 2
    numnodes = numnodes + i
    ! integer 32-bit limitation (integer overflow)
    if (numnodes > 2147483646 - i ) stop 'Error number of nodes might exceed integer limit'
  enddo
  if (be_verbose) then
    write(IMAIN,*) '  theoretical   number of nodes: ',numnodes
    write(IMAIN,*) '               tree memory size: ',( numnodes * 32 )/1024./1024.,'MB'
    call flush_IMAIN()
  endif

  ! local ordering
  allocate(points_index(kdtree_num_nodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1226')
  if (ier /= 0) stop 'Error allocating array points_index'
  points_index(:) = 0

  ! initial point ordering
  do i = 1,npoints
    points_index(i) = i
  enddo

  ! builds tree structure
  nullify(kdtree)
  depth = 0
  numnodes = 0
  maxdepth = -1

  call create_kdtree(npoints,points_data,points_index,kdtree, &
                     depth,1,npoints,numnodes,maxdepth)

  ! checks root tree node
  if (.not. associated(kdtree) ) stop 'Error creation of kd-tree failed'

  if (be_verbose) then
    write(IMAIN,*) '  actual        number of nodes: ',numnodes
    ! tree node size: 4 (idim) + 8 (cut_value) + 4 (ipoint) + 2*4 (ibound_**) + 2*4 (left,right) = 32 bytes
    write(IMAIN,*) '               tree memory size: ',( numnodes * 32 )/1024./1024.,'MB'
    write(IMAIN,*) '  maximum depth   : ',maxdepth

    ! timing
    call cpu_time(ct_end)
    write(IMAIN,*) '  creation timing : ',ct_end - ct_start, '(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! debugging
  if (DEBUG) then
    ! outputs tree
    if (be_verbose) then
      numnodes = 0
      call print_kdtree(npoints,points_data,points_index,kdtree,numnodes)
    endif

    ! test search
    print *,'search tree: rank ',myrank
    xyz_target(1) = 0.13261298835277557
    xyz_target(2) = -8.4083788096904755E-002
    xyz_target(3) = 0.97641450166702271

    print *,'search : ',xyz_target(:)

    ipoint_min = -1
    dist_min = 1.d30

    call find_nearest_node(npoints,points_data,kdtree,xyz_target,ipoint_min,dist_min)

    dist_min = sqrt(dist_min)
    print *,'found : ',ipoint_min,'distance:',dist_min

    if (ipoint_min < 1 ) stop 'Error search kd-tree found no point'

    print *,'target  : ',xyz_target(:)
    print *,'nearest : ',points_data(:,ipoint_min)
    print *
    ! safety stop
    stop 'kdtree_setup safety stop'
  endif

  ! frees temporary arrays
  deallocate(points_index)
  nullify(points_data)

  end subroutine kdtree_setup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

  ! kd-tree nearest neighbor search
  !
  ! input:
  !
  ! returns: global index iglob_min and distance dist_min to nearest point

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target

  double precision, intent(out) :: dist_min
  integer, intent(out) :: iglob_min

  ! local parameters
  integer :: ipoint_min

  ! initializes
  ipoint_min = -1
  iglob_min = -1
  dist_min = 1.d30

  ! searches closest node in kd-tree
  call find_nearest_node(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target,ipoint_min,dist_min)

  if (ipoint_min < 1 .or. ipoint_min > kdtree_num_nodes ) stop 'Error search kd-tree found no point'

  ! gets global index
  iglob_min = kdtree_nodes_index(ipoint_min)

  ! checks global index
  if (iglob_min < 1 ) stop 'Error minimum location has wrong global index in kdtree_find_nearest_neighbor'

  ! returns distance (non-squared)
  dist_min = sqrt( dist_min )

  ! debug
  !if (be_verbose) then
  !  print *,'target  : ',xyz_target(:)
  !  print *,'nearest : ',kdtree_nodes_location(:,ipoint_min),'distance:',dist_min,'(m)',ipoint_min,iglob_min
  !endif

  end subroutine kdtree_find_nearest_neighbor

!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_find_nearest_n_neighbors(xyz_target,search_radius,num_nodes)

  ! kd-tree n-search
  !
  ! input: target location (xyz_target) and search radius (search_radius)
  !
  ! returns: number of nodes within search radius (num_nodes),
  !          and
  !          fills array kdtree_search_index with search result

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target
  double precision, intent(in) :: search_radius

  integer, intent(out) :: num_nodes

  ! local parameters
  double precision :: r_squared

  ! initializes
  num_nodes = 0

  ! checks tree
  if (kdtree_num_nodes <= 0) return

  ! checks radius
  if (search_radius < 0.d0) &
    stop 'Error kd-tree search radius is invalid (negative value)'

  ! checks if result array length non-zero
  if (kdtree_search_num_nodes <= 0) &
    stop 'Please set kdtree_search_num_nodes and allocate kdtree_search_index first, before finding closest n neighbors'

  ! checks result array
  if (.not. allocated(kdtree_search_index)) &
    stop 'Please allocate kdtree_search_index first, before finding closest n neighbors'

  ! initializes search results
  kdtree_search_index(:) = 0

  ! distances are squared
  r_squared = search_radius**2

  ! searches closest nodes in kd-tree and returns results in array kdtree_search_index(:)
  call find_nearest_n_nodes(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target,r_squared,num_nodes,.true.)

  ! checks result
  !if (num_nodes /= kdtree_search_num_nodes) &
  !  stop 'Error kd-tree search paramater kdtree_search_num_nodes is different to number of nodes found'

  ! checks result bounds
  if (minval(kdtree_search_index(1:num_nodes)) <= 0) &
    stop 'Error kd-tree search found invalid search result'

  ! sorts values in increasing order
  if (DO_SORT_RESULTS) then
    call heap_sort( num_nodes, kdtree_search_index(1:num_nodes) )
  endif

  end subroutine kdtree_find_nearest_n_neighbors


!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_count_nearest_n_neighbors(xyz_target,search_radius,num_nodes)

  ! kd-tree n-search
  !
  ! input: target location (xyz_target) and search radius (search_radius)
  !
  ! returns: number of nodes within search radius (num_nodes),

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target
  double precision, intent(in) :: search_radius

  integer, intent(out) :: num_nodes

  ! local parameters
  double precision :: r_squared

  ! initializes
  num_nodes = 0

  ! checks tree
  if (kdtree_num_nodes <= 0) return

  ! checks radius
  if (search_radius < 0.d0) stop 'Error kd-tree search radius is invalid (negative value)'

  ! distances are squared always
  r_squared = search_radius**2

  ! searches closest nodes in kd-tree (only counts nodes and returns result in num_nodes)
  call find_nearest_n_nodes(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target,r_squared,num_nodes,.false.)

  end subroutine kdtree_count_nearest_n_neighbors


!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_get_nearest_n_neighbors(xyz_target,search_radius,num_nodes_get)

  ! kd-tree n-search
  !
  ! input: target location (xyz_target), search radius (search_radius), number of nodes to fill (num_nodes_get)
  !
  ! returns: fills array kdtree_search_index with search result

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target
  double precision, intent(in) :: search_radius

  integer, intent(inout) :: num_nodes_get

  ! local parameters
  double precision :: r_squared
  integer :: num_nodes

  ! initializes
  num_nodes = 0

  ! checks if anything to do
  if (num_nodes_get <= 0) return

  ! checks radius
  if (search_radius < 0.d0) &
    stop 'Error kd-tree search radius is invalid (negative value)'

  ! checks if result array length non-zero
  if (kdtree_search_num_nodes <= 0) &
    stop 'Please set kdtree_search_num_nodes and allocate kdtree_search_index first, before getting closest n neighbors'

  ! checks result array
  if (.not. allocated(kdtree_search_index)) &
    stop 'Please allocate kdtree_search_index first, before getting closest n neighbors'

  ! checks if num_nodes_get limited by search array size
  if (kdtree_search_num_nodes < num_nodes_get) then
    print *,'Warning: Requested number of n-nodes bigger than actual number of search result kdtree_search_num_nodes'
  endif

  ! initializes search results
  kdtree_search_index(:) = 0

  ! distances are squared
  r_squared = search_radius**2

  ! searches closest nodes in kd-tree and returns results in array kdtree_search_index(:)
  call find_nearest_n_nodes(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target,r_squared,num_nodes,.true.)

  ! updates return results
  if (num_nodes < num_nodes_get) num_nodes_get = num_nodes

  ! returns only requested number of nodes
  if (num_nodes > num_nodes_get) then
    kdtree_search_index(num_nodes_get+1:num_nodes) = 0
  endif

  ! checks result bounds
  if (minval(kdtree_search_index(1:num_nodes_get)) <= 0) &
    stop 'Error kd-tree search found invalid search result'

  ! sorts values in increasing order
  if (DO_SORT_RESULTS) then
    call heap_sort( num_nodes_get, kdtree_search_index(1:num_nodes_get) )
  endif

  end subroutine kdtree_get_nearest_n_neighbors


!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_count_nearest_n_neighbors_ellip(xyz_target,dist_v,dist_h,num_nodes)

  ! kd-tree n-search
  !
  ! input: target location (xyz_target) and search ellipsoid (dist_v/dist_h)
  !
  ! returns: number of nodes within search radius (num_nodes),

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target
  double precision, intent(in) :: dist_v,dist_h

  integer, intent(out) :: num_nodes

  ! local parameters
  double precision :: r_squared_v,r_squared_h

  ! initializes
  num_nodes = 0

  ! checks tree
  if (kdtree_num_nodes <= 0) return

  ! checks radius
  if (dist_v < 0.d0 .or. dist_h < 0.d0) &
    stop 'Error kd-tree search ellipsoid is invalid (negative value)'

  ! distances are squared
  r_squared_v = dist_v**2
  r_squared_h = dist_h**2

  ! searches closest nodes in kd-tree (only counts nodes and returns result in num_nodes)
  call find_nearest_n_nodes_ellip(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target, &
                                  r_squared_v,r_squared_h,num_nodes,.false.)

  end subroutine kdtree_count_nearest_n_neighbors_ellip


!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_get_nearest_n_neighbors_ellip(xyz_target,dist_v,dist_h,num_nodes_get)

  ! kd-tree n-search
  !
  ! input: target location (xyz_target), search ellipsoid (dist_v/disth),
  !        number of nodes to fill (num_nodes_get)
  !
  ! returns: fills array kdtree_search_index with search result

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target
  double precision, intent(in) :: dist_v,dist_h

  integer, intent(inout) :: num_nodes_get

  ! local parameters
  double precision :: r_squared_v,r_squared_h
  integer :: num_nodes

  ! initializes
  num_nodes = 0

  ! checks if anything to do
  if (num_nodes_get <= 0) return

  ! checks radius
  if (dist_v < 0.d0 .or. dist_h < 0.d0) &
    stop 'Error kd-tree search ellipsoid is invalid (negative value)'

  ! checks if result array length non-zero
  if (kdtree_search_num_nodes <= 0) &
    stop 'Please set kdtree_search_num_nodes and allocate kdtree_search_index first, before getting closest n neighbors'

  ! checks result array
  if (.not. allocated(kdtree_search_index)) &
    stop 'Please allocate kdtree_search_index first, before getting closest n neighbors'

  ! checks if num_nodes_get limited by search array size
  if (kdtree_search_num_nodes < num_nodes_get) then
    print *,'Warning: Requested number of n-nodes bigger than actual number of search result kdtree_search_num_nodes'
  endif

  ! initializes search results
  kdtree_search_index(:) = 0

  ! distances are squared
  r_squared_v = dist_v**2
  r_squared_h = dist_h**2

  ! searches closest nodes in kd-tree and returns results in array kdtree_search_index(:)
  call find_nearest_n_nodes_ellip(kdtree_num_nodes,kdtree_nodes_location,kdtree,xyz_target, &
                                  r_squared_v,r_squared_h,num_nodes,.true.)

  ! updates return results
  if (num_nodes < num_nodes_get) num_nodes_get = num_nodes

  ! returns only requested number of nodes
  if (num_nodes > num_nodes_get) then
    kdtree_search_index(num_nodes_get+1:num_nodes) = 0
  endif

  ! checks result bounds
  if (minval(kdtree_search_index(1:num_nodes_get)) <= 0) &
    stop 'Error kd-tree search found invalid search result'

  ! sorts values in increasing order
  if (DO_SORT_RESULTS) then
    call heap_sort( num_nodes_get, kdtree_search_index(1:num_nodes_get) )
  endif

  end subroutine kdtree_get_nearest_n_neighbors_ellip


!
!-------------------------------------------------------------------------------------------------
!

  recursive subroutine create_kdtree(npoints,points_data,points_index,node, &
                                     depth,ibound_lower,ibound_upper,numnodes,maxdepth)

  ! creates node in kd-tree structure

  implicit none

  integer, intent(in) :: npoints
  double precision, dimension(3,npoints), intent(in) :: points_data
  integer,dimension(npoints), intent(inout) :: points_index

  type (kdtree_node), pointer :: node    ! pointers in standard Fortran90 cannot have intent(..) attribute

  integer,intent(in) :: depth
  integer,intent(in) :: ibound_lower,ibound_upper

  integer,intent(inout) :: numnodes,maxdepth

  ! local parameters
  double precision :: cut_value
  double precision :: range,range_max,min,max
  integer :: i,ier,idim
  integer :: iloc,ilower,iupper
  integer :: l,u

  ! note: compiling with intel ifort version 18.0.1/19.1.0 and optimizations like -xHost -O2 or -xHost -O3 flags
  !       can lead to issues with the deallocate(workindex) statement below:
  !         *** Error in `./bin/xspecfem3D': double free or corruption (!prev): 0x00000000024f1610 ***
  !
  !       this might be due to a more aggressive optimization which leads to a change of the instruction set
  !       and the memory being free twice.
  !       a way to avoid this is by removing -xHost from FLAGS_CHECK = .. in Makefile
  !       or to use a pointer array instead of an allocatable array
  !
  ! integer,dimension(:),allocatable :: workindex
  integer,dimension(:),pointer :: workindex

  ! checks if anything to sort
  if (ibound_lower > ibound_upper) then
    nullify(node)
    return
  endif

  ! creates new node
  allocate(node,stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1227')
  if (ier /= 0) then
    print *,'Error creating node: ',numnodes
    stop 'Error allocating kd-tree node'
  endif

  ! initializes new node
  node%idim = -1
  node%ipoint = -1
  node%cut_value = 0.d0

  nullify(node%left)
  nullify(node%right)

  ! tree statistics
  numnodes = numnodes + 1
  if (maxdepth < depth) maxdepth = depth

  !node%id = numnodes

  ! checks if final node
  if (ibound_lower == ibound_upper) then
    node%idim = 0
    node%ipoint = points_index(ibound_lower)
    ! done with this node
    return
  endif

  ! sets cut dimension index

  ! version 1: varies between 1 and 3 depending on depth for 3D data set
  ! (leads to some unneccessary unbalanced nodes)
  !idim = mod(depth,3) + 1
  ! determines cut value
  ! range in this dimension
  !min = HUGEVAL
  !max = - HUGEVAL
  !do i = ibound_lower,ibound_upper
  !  iloc = points_index(i)
  !  val = points_data(idim,iloc)
  !  if (val < min ) min = val
  !  if (val > max ) max = val
  !enddo
  !min = minval(points_data(idim,points_index(ibound_lower:ibound_upper)))
  !max = maxval(points_data(idim,points_index(ibound_lower:ibound_upper)))
  !cut_value = 0.5d0 * ( min + max )

  ! version 2: selects cut dimension where biggest range is occurring
  ! (better balances tree)
  cut_value = 0.d0
  idim = -1
  range_max = 0.d0
  do i = 1,3
    min = minval(points_data(i,points_index(ibound_lower:ibound_upper)))
    max = maxval(points_data(i,points_index(ibound_lower:ibound_upper)))
    range = max - min
    ! sets cut dimension where data has maximum range
    if (range > range_max) then
      range_max = range
      idim = i
      cut_value = 0.5d0 * ( min + max )
      ! stores bounds of cut dimension
      !node%cut_min = min
      !node%cut_max = max
    endif
  enddo
  ! default dimension
  if (idim < 1) then
    ! in case we have two identical points:
    !   ibound_lower < ibound_upper but min == max value,
    ! thus zero range and idim,cut_value not set yet

    ! debug
    !print *,'create_kdtree: ',ibound_lower,ibound_upper
    !print *,'create_kdtree: data 1 min/max ', &
    !minval(points_data(1,points_index(ibound_lower:ibound_upper))),maxval(points_data(1,points_index(ibound_lower:ibound_upper)))
    !print *,'create_kdtree: data 2 min/max ', &
    !minval(points_data(2,points_index(ibound_lower:ibound_upper))),maxval(points_data(2,points_index(ibound_lower:ibound_upper)))
    !print *,'create_kdtree: data 3 min/max ', &
    !minval(points_data(3,points_index(ibound_lower:ibound_upper))),maxval(points_data(3,points_index(ibound_lower:ibound_upper)))

    ! default
    idim = 1
    cut_value = 0.d0
  endif
  ! sets node values
  node%idim = idim
  node%cut_value = cut_value

  !debug
  !print *,'index ',numnodes,'dim:',idim,'range:',ibound_lower,ibound_upper
  !print *,'  data:',points_data(idim,points_index(ibound_lower)),points_data(idim,points_index(ibound_upper))
  !print *,'  min/max:',min,max,'cut value:',cut_value

  ! temporary index array for sorting
  allocate(workindex(ibound_upper - ibound_lower + 1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1228')
  if (ier /= 0) stop 'Error allocating workindex array'
  workindex(:) = 0

  ! sorts point indices
  ! to have all points with value < cut_value on left side, all others to the right
  ilower = 0
  iupper = 0
  do i = ibound_lower,ibound_upper
    iloc = points_index(i)
    ! checks index
    if (iloc < 1) stop 'Error invalid iloc index in create_kdtree() routine'
    ! sorts tree
    if (points_data(idim,iloc) < cut_value) then
      ilower = ilower + 1
      workindex(ilower) = iloc
    else
      iupper = iupper + 1
      workindex(ibound_upper - ibound_lower + 2 - iupper) = iloc
    endif
  enddo
  !debug
  !print *,'  ilower/iupper:',ilower,iupper

  ! identical points: split first left, second right to balance tree
   if (ilower == 0) then
     ilower = 1
     iupper = (ibound_upper - ibound_lower)
   else if (iupper == 0) then
     ilower = (ibound_upper - ibound_lower)
     iupper = 1
   endif

  ! checks if we catched all
  if (ilower + iupper /= ibound_upper - ibound_lower + 1 ) stop 'Error sorting data points invalid'
  if (ilower == 0 .and. iupper == 0 .and. npoints > 1 ) stop 'Error confusing node counts, please check kdtree...'

  ! replaces index range with new sorting order
  points_index(ibound_lower:ibound_upper) = workindex(:)

  ! note: compiling with intel ifort version 18.0.1/19.1.0 and optimizations like -xHost -O2 or -xHost -O3 flags
  !       can lead to issues with the deallocate(workindex) statement below:
  !         *** Error in `./bin/xspecfem3D': double free or corruption (!prev): 0x00000000024f1610 ***
  !
  !       this might be due to a more aggressive optimization which leads to a change of the instruction set
  !       and the memory being free twice.
  !       a way to avoid this is by removing -xHost from FLAGS_CHECK = .. in Makefile
  !       or to use a pointer array instead of an allocatable array

  ! frees temporary array
  deallocate(workindex)

  ! lower hemisphere
  if (ilower > 0) then
    ! lower/upper bounds
    l = ibound_lower
    u = ibound_lower + (ilower - 1)
    if (l < 1 .or. u > npoints ) stop 'Error lower hemisphere tree bounds'
    ! adds new node
    call create_kdtree(npoints,points_data,points_index,node%left,depth+1,l,u,numnodes,maxdepth)
  endif

  ! upper hemisphere points
  if (iupper > 0) then
    ! lower/upper bounds
    l = ibound_upper - (iupper - 1)
    u = ibound_upper
    if (l < 1 .or. u > npoints ) stop 'Error upper hemisphere tree bounds'
    ! adds new node
    call create_kdtree(npoints,points_data,points_index,node%right,depth+1,l,u,numnodes,maxdepth)
  endif

  end subroutine create_kdtree

!
!-------------------------------------------------------------------------------------------------
!

  recursive subroutine print_kdtree(npoints,points_data,points_index,node,numnodes)

  ! prints out all final nodes in kd-tree structure

  implicit none

  integer, intent(in) :: npoints
  double precision, dimension(3,npoints),intent(in) :: points_data
  integer,dimension(npoints), intent(in) :: points_index

  type (kdtree_node), pointer :: node  ! pointers in standard Fortran90 cannot have intent(..) attribute

  integer, intent(inout) :: numnodes

  ! local parameters
  integer, parameter :: OUTPUT_LENGTH = 50

  ! checks if valid pointer (must have been nullified initially to be able to check with associated())
  if (.not. associated(node) ) return

  ! statistics
  numnodes = numnodes + 1
  if (numnodes == 1) then
    print *,'printing kd-tree: total number of points      = ',npoints
    !print *,'         index array = ',points_index(:)
  endif

  ! outputs infos for a final node
  if (.not. associated(node%left) .and. .not. associated(node%right)) then
    ! checks info
    if (node%idim /= 0) then
      print *,'problem kd-tree node:',node%idim,node%ipoint,numnodes
      print *,'point x/y/z: ',points_data(:,node%ipoint)
      stop 'Error kd-tree node not correct'
    endif

    ! outputs infos
    if (numnodes < OUTPUT_LENGTH) &
      print *,'node:',numnodes,'index:',node%ipoint,' x/y/z = ',points_data(:,node%ipoint)
  else
    ! outputs infos
    if (numnodes < OUTPUT_LENGTH) &
      print *,'node:',numnodes,'dim:',node%idim,'cut = ',node%cut_value
  endif

  ! checks child nodes
  if (associated(node%left)) then
    call print_kdtree(npoints,points_data,points_index,node%left,numnodes)
  endif
  if (associated(node%right)) then
    call print_kdtree(npoints,points_data,points_index,node%right,numnodes)
  endif

  end subroutine print_kdtree

!
!-------------------------------------------------------------------------------------------------
!

  subroutine kdtree_delete()

  ! deletes all (child) nodes in this given tree

  implicit none

  ! deletes tree starting with top node
  call delete_node(kdtree)

  end subroutine kdtree_delete

!
!-------------------------------------------------------------------------------------------------
!

  recursive subroutine delete_node(node)

  ! deletes all (child) nodes in this given tree

  implicit none

  type (kdtree_node), pointer :: node

  ! delete left hemisphere
  if (associated(node%left)) then
     call delete_node(node%left)
     nullify (node%left)
  endif

  ! deletes right hemisphere
  if (associated(node%right)) then
     call delete_node(node%right)
     nullify (node%right)
  endif

  deallocate(node)

  end subroutine delete_node

!
!-------------------------------------------------------------------------------------------------
!

  recursive subroutine find_nearest_node(npoints,points_data,node,xyz_target,ipoint_min,dist_min)

  ! searches for node point closest to given location
  implicit none
  integer, intent(in) :: npoints
  double precision, dimension(3,npoints), intent(in) :: points_data

  type (kdtree_node), pointer :: node  ! pointers in standard Fortran90 cannot have intent(..) attribute

  double precision, dimension(3), intent(in) :: xyz_target

  integer, intent(inout) :: ipoint_min
  double precision, intent(inout) :: dist_min

  ! local parameters
  double precision :: dist

  ! debug
  !if (node%idim == 0) then
  !  print *,'node',node%id,points_data(:,node%ipoint)
  !else
  !  print *,'node',node%id,node%idim,node%cut_value
  !endif
  !if (ipoint_min > 0) &
  !  print *,'node distance',node%id,ipoint_min,dist_min

  ! in case this is a final node
  if (.not. associated(node%left) .and. .not. associated(node%right)) then
    ! checks node
    if (node%idim /= 0 ) stop 'Error searched node is not final node'
    if (node%ipoint < 1 ) stop 'Error searched node has wrong point index'

    ! squared distance to associated data point
    dist = get_distance_squared(xyz_target,points_data(1,node%ipoint))

    ! note: using <= instead of < for comparison. both would be fine, but the first leads to identical location result
    !       as with a brute force search, if the target location is exactly on a shared GLL point.
    !       the latter would choose a different element and lead to slightly different seismograms - not sure though why...
    !       it obviously matters if the source point is shared between different elements and the source contribution added by
    !       only a single element. for such cases, we might need to spread the source contribution to all shared elements.
    if (dist <= dist_min) then
      ! debug
      !if (ipoint_min < 1) then
      !  print *,'new node distance',node%id,node%ipoint,dist
      !else
      !  print *,'     new distance',node%id,node%ipoint,dist
      !endif
      ! stores minimum point
      dist_min = dist
      ipoint_min = node%ipoint
    endif

    ! done
    return
  endif

  ! checks cut dimension
  if (node%idim < 1 .or. node%idim > 3 ) stop 'Error searched node has invalid cut dimension'

  ! compares cut value
  if (xyz_target(node%idim) < node%cut_value) then
    ! finds closer node in lower hemisphere
    if (associated(node%left)) then
      call find_nearest_node(npoints,points_data,node%left,xyz_target,ipoint_min,dist_min)
    endif
  else
    ! finds closer node in upper hemisphere
    if (associated(node%right)) then
      call find_nearest_node(npoints,points_data,node%right,xyz_target,ipoint_min,dist_min)
    endif
  endif

  ! at this point, dist_min is the distance to the closest point in the initial hemisphere search
  ! we might need to search in other hemisphere as well if distances are closer

  ! squared distance to cut plane
  dist = ( xyz_target(node%idim) - node%cut_value )**2

  if (xyz_target(node%idim) < node%cut_value) then
    if (associated(node%right)) then
      ! checks right node as a final node
      if (node%right%idim == 0) then
        dist = get_distance_squared(xyz_target,points_data(1,node%right%ipoint))
        if (dist <= dist_min) then
          ! stores minimum point
          dist_min = dist
          ipoint_min = node%right%ipoint
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist < dist_min) then
        call find_nearest_node(npoints,points_data,node%right,xyz_target,ipoint_min,dist_min)
      endif
    endif
  else
    if (associated(node%left)) then
      ! checks left node as a final node
      if (node%left%idim == 0) then
        dist = get_distance_squared(xyz_target,points_data(1,node%left%ipoint))
        if (dist <= dist_min) then
          ! stores minimum point
          dist_min = dist
          ipoint_min = node%left%ipoint
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist < dist_min) then
        call find_nearest_node(npoints,points_data,node%left,xyz_target,ipoint_min,dist_min)
      endif
    endif
  endif

  end subroutine find_nearest_node


!
!-------------------------------------------------------------------------------------------------
!


  recursive subroutine find_nearest_n_nodes(npoints,points_data,node,xyz_target,r_squared,num_nodes,fill_index)

  ! searches for all node points within a given radius to given location

  implicit none
  integer, intent(in) :: npoints
  double precision, dimension(3,npoints), intent(in) :: points_data

  type (kdtree_node), pointer :: node   ! pointers in standard Fortran90 cannot have intent(..) attribute

  double precision,dimension(3), intent(in) :: xyz_target

  double precision, intent(in) :: r_squared
  integer, intent(inout) :: num_nodes

  logical, intent(in) :: fill_index

  ! local parameters
  double precision :: dist
  double precision,dimension(3) :: xyz

  ! checks a final node
  if (.not. associated(node%left) .and. .not. associated(node%right)) then
    ! checks node
    if (node%idim /= 0 ) stop 'Error searched node is not final node'
    if (node%ipoint < 1 ) stop 'Error searched node has wrong point index'

    ! node location
    xyz(:) = points_data(:,node%ipoint)

    ! squared distance to associated data point
    dist = get_distance_squared(xyz_target,xyz)
    if (dist <= r_squared) then
      ! debug
      !print *,'     new node: ',node%ipoint,'distance = ',dist,'radius = ',r_squared
      ! counts point
      num_nodes = num_nodes + 1

      ! adds point
      if (fill_index) then
        if (num_nodes <= kdtree_search_num_nodes) then
          kdtree_search_index(num_nodes) = node%ipoint
        endif
      endif
    endif

    ! done
    return
  endif

  ! checks cut dimension
  if (node%idim < 1 .or. node%idim > 3 ) stop 'Error searched node has invalid cut dimension'

  ! compares cut value
  if (xyz_target(node%idim) < node%cut_value) then
    ! finds closer node in lower hemisphere
    if (associated(node%left)) then
      call find_nearest_n_nodes(npoints,points_data,node%left,xyz_target,r_squared,num_nodes,fill_index)
    endif
  else
    ! finds closer node in upper hemisphere
    if (associated(node%right)) then
      call find_nearest_n_nodes(npoints,points_data,node%right,xyz_target,r_squared,num_nodes,fill_index)
    endif
  endif

  ! at this point, dist_min is the distance to the closest point in the initial hemisphere search
  ! we might need to search in other hemisphere as well if distances are closer

  ! squared distance to cut plane
  dist = ( xyz_target(node%idim) - node%cut_value )**2

  if (xyz_target(node%idim) < node%cut_value) then
    if (associated(node%right)) then
      ! checks right node as a final node
      if (node%right%idim == 0) then
        xyz(:) = points_data(:,node%right%ipoint)
        dist = get_distance_squared(xyz_target,xyz)
        if (dist <= r_squared) then
          ! counts point
          num_nodes = num_nodes + 1
          ! adds point
          if (fill_index) then
            if (num_nodes <= kdtree_search_num_nodes) then
              kdtree_search_index(num_nodes) = node%right%ipoint
            endif
          endif
          ! done
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist <= r_squared) then
        call find_nearest_n_nodes(npoints,points_data,node%right,xyz_target,r_squared,num_nodes,fill_index)
      endif
    endif
  else
    if (associated(node%left)) then
      ! checks left node as a final node
      if (node%left%idim == 0) then
        xyz(:) = points_data(:,node%left%ipoint)
        dist = get_distance_squared(xyz_target,xyz)
        if (dist <= r_squared) then
          ! counts point
          num_nodes = num_nodes + 1
          ! adds point
          if (fill_index) then
            if (num_nodes <= kdtree_search_num_nodes) then
              kdtree_search_index(num_nodes) = node%left%ipoint
            endif
          endif
          ! done
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist <= r_squared) then
        call find_nearest_n_nodes(npoints,points_data,node%left,xyz_target,r_squared,num_nodes,fill_index)
      endif
    endif
  endif

  end subroutine find_nearest_n_nodes

!
!-------------------------------------------------------------------------------------------------
!


  recursive subroutine find_nearest_n_nodes_ellip(npoints,points_data,node,xyz_target, &
                                                  r_squared_v,r_squared_h,num_nodes,fill_index)

  ! searches for all node points within a given search ellipsoid to given location

  implicit none
  integer, intent(in) :: npoints
  double precision, dimension(3,npoints), intent(in) :: points_data

  type (kdtree_node), pointer :: node    ! pointers in standard Fortran90 cannot have intent(..) attribute

  double precision, dimension(3), intent(in) :: xyz_target

  double precision, intent(in) :: r_squared_v,r_squared_h
  integer, intent(inout) :: num_nodes

  logical, intent(in) :: fill_index

  ! local parameters
  double precision :: dist,dist_v,dist_h
  double precision,dimension(3) :: xyz

  ! checks a final node
  if (.not. associated(node%left) .and. .not. associated(node%right)) then
    ! checks node
    if (node%idim /= 0 ) stop 'Error searched node is not final node'
    if (node%ipoint < 1 ) stop 'Error searched node has wrong point index'

    ! node location
    xyz(:) = points_data(:,node%ipoint)

    ! squared distance to associated data point
    call get_distance_ellip(xyz_target,xyz,dist_v,dist_h)
    if (dist_v <= r_squared_v .and. dist_h <= r_squared_h) then
      ! debug
      !print *,'     new node: ',node%ipoint,'distance = ',dist,'radius = ',r_squared
      ! counts point
      num_nodes = num_nodes + 1

      ! adds point
      if (fill_index) then
        if (num_nodes <= kdtree_search_num_nodes) then
          kdtree_search_index(num_nodes) = node%ipoint
        endif
      endif
    endif

    ! done
    return
  endif

  ! checks cut dimension
  if (node%idim < 1 .or. node%idim > 3 ) stop 'Error searched node has invalid cut dimension'

  ! compares cut value
  if (xyz_target(node%idim) < node%cut_value) then
    ! finds closer node in lower hemisphere
    if (associated(node%left)) then
      call find_nearest_n_nodes_ellip(npoints,points_data,node%left,xyz_target, &
                                      r_squared_v,r_squared_h,num_nodes,fill_index)
    endif
  else
    ! finds closer node in upper hemisphere
    if (associated(node%right)) then
      call find_nearest_n_nodes_ellip(npoints,points_data,node%right,xyz_target, &
                                      r_squared_v,r_squared_h,num_nodes,fill_index)
    endif
  endif

  ! at this point, dist_min is the distance to the closest point in the initial hemisphere search
  ! we might need to search in other hemisphere as well if distances are closer

  ! squared distance to cut plane
  dist = ( xyz_target(node%idim) - node%cut_value )**2

  if (xyz_target(node%idim) < node%cut_value) then
    if (associated(node%right)) then
      ! checks right node as a final node
      if (node%right%idim == 0) then
        xyz(:) = points_data(:,node%right%ipoint)
        call get_distance_ellip(xyz_target,xyz,dist_v,dist_h)
        if (dist_v <= r_squared_v .and. dist_h <= r_squared_h) then
          ! counts point
          num_nodes = num_nodes + 1
          ! adds point
          if (fill_index) then
            if (num_nodes <= kdtree_search_num_nodes) then
              kdtree_search_index(num_nodes) = node%right%ipoint
            endif
          endif
          ! done
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist <= r_squared_v .or. dist <= r_squared_h) then
        call find_nearest_n_nodes_ellip(npoints,points_data,node%right,xyz_target, &
                                        r_squared_v,r_squared_h,num_nodes,fill_index)
      endif
    endif
  else
    if (associated(node%left)) then
      ! checks left node as a final node
      if (node%left%idim == 0) then
        xyz(:) = points_data(:,node%left%ipoint)
        call get_distance_ellip(xyz_target,xyz,dist_v,dist_h)
        if (dist_v <= r_squared_v .and. dist_h <= r_squared_h) then
          ! counts point
          num_nodes = num_nodes + 1
          ! adds point
          if (fill_index) then
            if (num_nodes <= kdtree_search_num_nodes) then
              kdtree_search_index(num_nodes) = node%left%ipoint
            endif
          endif
          ! done
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist <= r_squared_v .or. dist <= r_squared_h) then
        call find_nearest_n_nodes_ellip(npoints,points_data,node%left,xyz_target, &
                                        r_squared_v,r_squared_h,num_nodes,fill_index)
      endif
    endif
  endif

  end subroutine find_nearest_n_nodes_ellip

!
!-------------------------------------------------------------------------------------------------
!


  subroutine kdtree_set_verbose()

  implicit none

  ! sets verbosity on
  be_verbose = .true.

  end subroutine kdtree_set_verbose


!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_distance_squared(xyz0,xyz1)

  implicit none
  double precision,dimension(3),intent(in) :: xyz0,xyz1

  ! local parameters
  double precision :: dist
  !double precision,dimension(3) :: xyz

  ! calculates distance (squared) between 2 points
  dist = sum((xyz0(:) - xyz1(:))**2)

  ! explicit alternative calculation
  !xyz(:) = xyz0(:) - xyz1(:)
  !dist = xyz(1) * xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)

  ! returns distance squared
  get_distance_squared = dist

  end function get_distance_squared

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_distance_ellip(xyz0,xyz1,dist_v,dist_h)

  implicit none
  double precision,dimension(3),intent(in) :: xyz0,xyz1
  double precision,intent(out) :: dist_v,dist_h

  ! local parameters
  double precision :: r0,r1
  double precision :: theta,ratio

  ! vertical distance (squared)
  r0 = sqrt( sum( xyz0(:)**2 ) ) ! length of first position vector
  r1 = sqrt( sum( xyz1(:)**2 ) )

  dist_v = (r1 - r0)*(r1 - r0)

  ! only for flat earth with z in depth: dist_v = sqrt( (cz(ispec2)-cz0(ispec))** 2)

  ! epicentral distance
  ! (accounting for spherical curvature)
  ! calculates distance of circular segment
  ! angle between r0 and r1 in radian
  ! given by dot-product of two vectors
  ratio = sum(xyz0(:)*xyz1(:)) / (r0 * r1)

  ! checks boundaries of ratio (due to numerical inaccuracies)
  if (ratio > 1.d0) ratio = 1.d0
  if (ratio < -1.d0) ratio = -1.d0

  theta = acos( ratio )

  ! segment length at heigth of r1 (squared)
  dist_h = (r1 * theta)*(r1 * theta)

  end subroutine get_distance_ellip

end module
