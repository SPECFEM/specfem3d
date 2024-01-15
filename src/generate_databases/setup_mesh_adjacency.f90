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


  subroutine setup_mesh_adjacency()

! setups mesh adjacency array to search element neighbors for point searches

  use constants, only: myrank, &
    NDIM,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,IMAIN,CUSTOM_REAL,MAX_STRING_LEN, &
    DO_BRUTE_FORCE_POINT_SEARCH

  use generate_databases_par, only: NSPEC_AB,NGLOB_AB,ibool,NPROC,prname

  use create_regions_mesh_ext_par, only: xstore => xstore_unique, ystore => ystore_unique, zstore => zstore_unique

  ! mesh adjacency for point searches
  use generate_databases_par, only: neighbors_xadj,neighbors_adjncy,num_neighbors_all

  use kdtree_search, only: kdtree_setup,kdtree_delete, &
    kdtree_nodes_location,kdtree_nodes_index,kdtree_num_nodes, &
    kdtree_count_nearest_n_neighbors,kdtree_get_nearest_n_neighbors, &
    kdtree_search_index,kdtree_search_num_nodes

  use fault_generate_databases, only: ANY_FAULT_IN_THIS_PROC

  implicit none

  ! local parameters
  ! maximum number of neighbors
  !integer,parameter :: MAX_NEIGHBORS = 50         ! maximum number of neighbors
  integer,parameter :: MAX_NEIGHBORS = 300         ! maximum number of neighbors (with neighbor of neighbors)

  ! temporary
  integer,dimension(:),allocatable :: tmp_adjncy  ! temporary adjacency
  integer :: inum_neighbor

  ! coordinates of element midpoints
  double precision, dimension(:,:), allocatable :: xyz_midpoints

  ! timer MPI
  double precision :: time1,tCPU
  double precision, external :: wtime

  integer :: ielem_counter,num_elements_actual_max
  integer :: num_neighbors,num_neighbors_max
  integer :: ispec_ref,ispec,iglob,ier !icorner,jj

  double precision :: xyz_target(NDIM)
  double precision :: dist_squared,dist_squared_max
  double precision :: maximal_elem_size_squared

  real(kind=CUSTOM_REAL) :: distance_min_glob,distance_max_glob
  real(kind=CUSTOM_REAL) :: elemsize_min_glob,elemsize_max_glob
  real(kind=CUSTOM_REAL) :: x_min_glob,x_max_glob
  real(kind=CUSTOM_REAL) :: y_min_glob,y_max_glob
  real(kind=CUSTOM_REAL) :: z_min_glob,z_max_glob

  ! for all the elements in contact with the reference element
  integer, dimension(:,:), allocatable :: ibool_corner
  logical, dimension(:), allocatable :: mask_ispec
  logical, dimension(:), allocatable :: flag_topological  ! making array allocatable, otherwise will crash for large meshes
  logical, dimension(8) :: flag_iglob_corner

  ! kd-tree search
  integer :: nsearch_points,inodes
  integer :: num_elements,num_elements_max
  ! alternative: to avoid allocating/deallocating search index arrays, though there is hardly a speedup
  !integer, parameter :: max_search_points = 20000
  double precision :: r_search
  logical :: use_kdtree_search

  logical :: is_neighbor

  ! neighbors of neighbors
  integer :: ielem,ii,jj
  integer :: num_neighbor_neighbors,num_neighbor_neighbors_max

  ! note: we add direct neighbors plus neighbors of neighbors.
  !       for very coarse meshes, the initial location guesses especially around doubling layers can be poor such that we need
  !       to enlarge the search of neighboring elements.
  logical, parameter :: ADD_NEIGHBOR_OF_NEIGHBORS = .true.

  ! debugging
  character(len=MAX_STRING_LEN) :: filename
  integer,dimension(:),allocatable :: tmp_num_neighbors
  logical,parameter :: DEBUG_VTK_OUTPUT = .false.

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '     mesh adjacency:'
    write(IMAIN,*) '     total number of elements in this slice  = ',NSPEC_AB
    write(IMAIN,*)
    write(IMAIN,*) '     maximum number of neighbors allowed     = ',MAX_NEIGHBORS
    write(IMAIN,*) '     minimum array memory required per slice = ',((MAX_NEIGHBORS + 1) * NSPEC_AB * 4)/1024./1024.,"(MB)"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time1 = wtime()

  ! prepares midpoints coordinates
  allocate(xyz_midpoints(NDIM,NSPEC_AB),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array xyz_midpoints')
  xyz_midpoints(:,:) = 0.d0

  ! store x/y/z coordinates of center point
  do ispec = 1,NSPEC_AB
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)
    xyz_midpoints(1,ispec) =  dble(xstore(iglob))
    xyz_midpoints(2,ispec) =  dble(ystore(iglob))
    xyz_midpoints(3,ispec) =  dble(zstore(iglob))
  enddo

  ! compute typical size of elements
  ! gets mesh dimensions
  call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                            x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
                            elemsize_min_glob,elemsize_max_glob, &
                            distance_min_glob,distance_max_glob)

  ! adjacency arrays
  !
  ! how to use:
  !  num_neighbors = neighbors_xadj(ispec+1)-neighbors_xadj(ispec)
  !  do i = 1,num_neighbors
  !    ! get neighbor
  !    ispec_neighbor = neighbors_adjncy(neighbors_xadj(ispec) + i)
  !    ..
  !  enddo
  allocate(neighbors_xadj(NSPEC_AB + 1),stat=ier)
  if (ier /= 0) stop 'Error allocating xadj'
  neighbors_xadj(:) = 0

  ! temporary helper array
  allocate(tmp_adjncy(MAX_NEIGHBORS * NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_adjncy'
  tmp_adjncy(:) = 0

  ! element mask
  allocate(mask_ispec(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating mask_ispec array'
  mask_ispec(:) = .false.

  ! since we only need to check corner points for the adjacency,
  ! we build an extra ibool array with corner points only for faster accessing
  allocate(ibool_corner(8,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool_corner array'
  ibool_corner(:,:) = 0
  do ispec = 1,NSPEC_AB
    ibool_corner(1,ispec) = ibool(1,1,1,ispec)
    ibool_corner(2,ispec) = ibool(NGLLX,1,1,ispec)
    ibool_corner(3,ispec) = ibool(NGLLX,NGLLY,1,ispec)
    ibool_corner(4,ispec) = ibool(1,NGLLY,1,ispec)

    ibool_corner(5,ispec) = ibool(1,1,NGLLZ,ispec)
    ibool_corner(6,ispec) = ibool(NGLLX,1,NGLLZ,ispec)
    ibool_corner(7,ispec) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    ibool_corner(8,ispec) = ibool(1,NGLLY,NGLLZ,ispec)
  enddo

  ! topological flags to find neighboring elements
  allocate(flag_topological(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating flag_topological array'
  flag_topological(:) = .false.

  ! use 10 times the distance as a criterion for point detection
  maximal_elem_size_squared = (10. * elemsize_max_glob)**2

  ! select search type
  if (DO_BRUTE_FORCE_POINT_SEARCH) then
    use_kdtree_search = .false.
  else
    ! for kdtree search, the number of elements per slice should be large enough,
    ! otherwise a simple brute force search is faster.
    if (NSPEC_AB > 20000) then
      use_kdtree_search = .true.
    else
      use_kdtree_search = .false.
    endif
  endif

  ! kd-tree search
  if (use_kdtree_search) then
    ! kd-tree search
    ! search radius around element midpoints
    !
    ! note: we take here 3.5 times the maximum element size to include also neighbors of neighbor elements
    !       - since at low resolutions NEX and large element sizes, this search radius needs to be large enough, and,
    !       - due to doubling layers (elements at depth will become bigger)
    !
    !       the search radius r_search given as routine argument must be non-squared
    r_search = 3.5 * elemsize_max_glob

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '     using kd-tree search radius             = ',sngl(r_search)
      call flush_IMAIN()
    endif

    ! kd-tree setup for adjacency search
    !
    ! uses only element midpoint location
    kdtree_num_nodes = NSPEC_AB

    ! allocates tree arrays
    allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
    allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

    ! prepares search arrays, each element takes its internal GLL points for tree search
    kdtree_nodes_index(:) = 0
    kdtree_nodes_location(:,:) = 0.0
    ! adds tree nodes
    inodes = 0
    do ispec = 1,NSPEC_AB
      ! sets up tree nodes
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)

      ! counts nodes
      inodes = inodes + 1
      if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

      ! adds node index (index points to same ispec for all internal GLL points)
      kdtree_nodes_index(inodes) = ispec

      ! adds node location
      kdtree_nodes_location(1,inodes) = xstore(iglob)
      kdtree_nodes_location(2,inodes) = ystore(iglob)
      kdtree_nodes_location(3,inodes) = zstore(iglob)
    enddo
    if (inodes /= kdtree_num_nodes ) stop 'Error index inodes does not match nnodes_local'

    ! alternative: to avoid allocating/deallocating search index arrays, though there is hardly a speedup
    !allocate(kdtree_search_index(max_search_points),stat=ier)
    !if (ier /= 0) stop 'Error allocating array kdtree_search_index'

    ! creates kd-tree for searching
    call kdtree_setup()
  endif

  ! gets maximum number of neighbors
  inum_neighbor = 0
  num_neighbors_max = 0
  num_neighbors_all = 0
  num_neighbor_neighbors_max = 0

  num_elements_max = 0
  num_elements_actual_max = 0
  dist_squared_max = 0.d0

  do ispec_ref = 1,NSPEC_AB
    ! flagging corners
    flag_topological(:) = .false.
    mask_ispec(:) = .false.

    ! mark the eight corners of the initial guess element
    flag_topological(ibool_corner(:,ispec_ref)) = .true.
    ! or w/ explicit ibool
    !flag_topological(ibool(1,1,1,ispec_ref)) = .true.
    !flag_topological(ibool(NGLLX,1,1,ispec_ref)) = .true.
    !flag_topological(ibool(NGLLX,NGLLY,1,ispec_ref)) = .true.
    !flag_topological(ibool(1,NGLLY,1,ispec_ref)) = .true.
    !flag_topological(ibool(1,1,NGLLZ,ispec_ref)) = .true.
    !flag_topological(ibool(NGLLX,1,NGLLZ,ispec_ref)) = .true.
    !flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec_ref)) = .true.
    !flag_topological(ibool(1,NGLLY,NGLLZ,ispec_ref)) = .true.

    ! marks element for checking elements
    mask_ispec(ispec_ref) = .true.

    ! midpoint for search radius
    iglob = ibool(MIDX,MIDY,MIDZ,ispec_ref)
    xyz_target(1) = xstore(iglob)
    xyz_target(2) = ystore(iglob)
    xyz_target(3) = zstore(iglob)

    if (use_kdtree_search) then
      ! looks only at elements in kd-tree search radius

      ! gets number of tree points within search radius
      ! (within search sphere)
      call kdtree_count_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

      ! debug
      !print *,'  total number of search elements: ',nsearch_points,'ispec',ispec_ref

      ! alternative: limits search results
      !if (nsearch_points > max_search_points) nsearch_points = max_search_points
      ! or stop if exceeded
      !if (nsearch_points > max_search_points) stop 'Error increase max_search_points'

      ! sets number of search nodes to get
      kdtree_search_num_nodes = nsearch_points

      ! allocates search index - dynamic instead of setting max_search_points
      allocate(kdtree_search_index(kdtree_search_num_nodes),stat=ier)
      if (ier /= 0) stop 'Error allocating array kdtree_search_index'

      ! gets closest n points around target (within search sphere)
      call kdtree_get_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

      ! loops over search radius
      num_elements = nsearch_points
    else
      ! loops over all other elements to find closest neighbors
      num_elements = NSPEC_AB
    endif

    ! statistics
    if (num_elements > num_elements_max) num_elements_max = num_elements

    ! counts number of neighbors
    num_neighbors = 0
    ielem_counter = 0

    ! loops over all other elements to find closest neighbors
    do ielem = 1,num_elements
      ! gets element index
      if (use_kdtree_search) then
        ! kd-tree search radius
        ! gets search point/element index
        ispec = kdtree_search_index(ielem)
        ! checks index
        if (ispec < 1 .or. ispec > NSPEC_AB) stop 'Error element index is invalid'
      else
        ! loops over all elements
        ispec = ielem
      endif

      ! skip reference element
      if (ispec == ispec_ref) cycle

      ! exclude elements that are too far from target
      ! distance to reference element
      dist_squared = (xyz_target(1) - xyz_midpoints(1,ispec))*(xyz_target(1) - xyz_midpoints(1,ispec)) &
                   + (xyz_target(2) - xyz_midpoints(2,ispec))*(xyz_target(2) - xyz_midpoints(2,ispec)) &
                   + (xyz_target(3) - xyz_midpoints(3,ispec))*(xyz_target(3) - xyz_midpoints(3,ispec))

      !  we compare squared distances instead of distances themselves to significantly speed up calculations
      if (dist_squared > maximal_elem_size_squared) cycle

      ielem_counter = ielem_counter + 1

      ! checks if element has a corner iglob from reference element
      is_neighbor = .false.

      ! corner flags
      flag_iglob_corner(:) = flag_topological(ibool_corner(:,ispec))
      ! or w/ explicit ibool
      !flag_iglob_corner(1) = flag_topological(ibool(1,1,1,ispec))
      !flag_iglob_corner(2) = flag_topological(ibool(NGLLX,1,1,ispec))
      !flag_iglob_corner(3) = flag_topological(ibool(NGLLX,NGLLY,1,ispec))
      !flag_iglob_corner(4) = flag_topological(ibool(1,NGLLY,1,ispec))
      !flag_iglob_corner(5) = flag_topological(ibool(1,1,NGLLZ,ispec))
      !flag_iglob_corner(6) = flag_topological(ibool(NGLLX,1,NGLLZ,ispec))
      !flag_iglob_corner(7) = flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      !flag_iglob_corner(8) = flag_topological(ibool(1,NGLLY,NGLLZ,ispec))

      ! checks if corner also has reference element
      if (any(flag_iglob_corner(:) .eqv. .true.)) then
        is_neighbor = .true.
      endif

      ! counts neighbors to reference element
      if (is_neighbor) then
        ! marks element for checking later on
        mask_ispec(ispec) = .true.

        ! adds to adjacency
        inum_neighbor = inum_neighbor + 1
        ! checks
        if (inum_neighbor > MAX_NEIGHBORS * NSPEC_AB) stop 'Error maximum neighbors exceeded'

        ! adds element
        tmp_adjncy(inum_neighbor) = ispec

        ! for statistics
        num_neighbors = num_neighbors + 1

        ! maximum distance to reference element
        if (dist_squared > dist_squared_max) dist_squared_max = dist_squared
      endif
    enddo

    ! again loop to get neighbors of neighbors
    if (ADD_NEIGHBOR_OF_NEIGHBORS) then
      ! counter for statistics
      num_neighbor_neighbors = 0

      ! flag first neighbor elements corner
      do ii = 1,num_neighbors
        ! get neighbor
        ispec = tmp_adjncy(inum_neighbor - num_neighbors + ii)

        ! mark the eight corners of the initial guess element
        flag_topological(ibool_corner(:,ispec)) = .true.
        ! or w/ explicit ibool
        !flag_topological(ibool(1,1,1,ispec)) = .true.
        !flag_topological(ibool(NGLLX,1,1,ispec)) = .true.
        !flag_topological(ibool(NGLLX,NGLLY,1,ispec)) = .true.
        !flag_topological(ibool(1,NGLLY,1,ispec)) = .true.
        !flag_topological(ibool(1,1,NGLLZ,ispec)) = .true.
        !flag_topological(ibool(NGLLX,1,NGLLZ,ispec)) = .true.
        !flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec)) = .true.
        !flag_topological(ibool(1,NGLLY,NGLLZ,ispec)) = .true.

        ! check element
        if (.not. mask_ispec(ispec)) stop 'error element flag topological'
      enddo

      ! full loop to get neighbors of neighbors
      do ielem = 1,num_elements
        ! gets element index
        if (use_kdtree_search) then
          ! kd-tree search radius
          ! gets search point/element index
          ispec = kdtree_search_index(ielem)
        else
          ! loops over all elements
          ispec = ielem
        endif

        ! skip reference element
        if (ispec == ispec_ref) cycle

        ! check if already added
        if (mask_ispec(ispec)) cycle

        ! exclude elements that are too far from target
        ! distance to reference element
        dist_squared = (xyz_target(1) - xyz_midpoints(1,ispec))*(xyz_target(1) - xyz_midpoints(1,ispec)) &
                     + (xyz_target(2) - xyz_midpoints(2,ispec))*(xyz_target(2) - xyz_midpoints(2,ispec)) &
                     + (xyz_target(3) - xyz_midpoints(3,ispec))*(xyz_target(3) - xyz_midpoints(3,ispec))

        ! exclude elements that are too far from target
        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        if (dist_squared > maximal_elem_size_squared) cycle

        ! checks if element has a corner iglob from reference element
        is_neighbor = .false.

        ! corner flags
        flag_iglob_corner(:) = flag_topological(ibool_corner(:,ispec))
        ! or w/ explicit ibool
        !flag_iglob_corner(1) = flag_topological(ibool(1,1,1,ispec))
        !flag_iglob_corner(2) = flag_topological(ibool(NGLLX,1,1,ispec))
        !flag_iglob_corner(3) = flag_topological(ibool(NGLLX,NGLLY,1,ispec))
        !flag_iglob_corner(4) = flag_topological(ibool(1,NGLLY,1,ispec))
        !flag_iglob_corner(5) = flag_topological(ibool(1,1,NGLLZ,ispec))
        !flag_iglob_corner(6) = flag_topological(ibool(NGLLX,1,NGLLZ,ispec))
        !flag_iglob_corner(7) = flag_topological(ibool(NGLLX,NGLLY,NGLLZ,ispec))
        !flag_iglob_corner(8) = flag_topological(ibool(1,NGLLY,NGLLZ,ispec))

        ! checks if corner also has reference element
        if (any(flag_iglob_corner(:) .eqv. .true.)) then
          is_neighbor = .true.
        endif

        ! counts neighbors to reference element
        if (is_neighbor) then
          ! marks element
          mask_ispec(ispec) = .true.

          ! adds to adjacency
          inum_neighbor = inum_neighbor + 1
          ! checks
          if (inum_neighbor > MAX_NEIGHBORS * NSPEC_AB) stop 'Error maximum neighbors with neighbors of neighbors exceeded'

          ! adds element
          tmp_adjncy(inum_neighbor) = ispec

          ! for statistics
          num_neighbors = num_neighbors + 1
          num_neighbor_neighbors = num_neighbor_neighbors + 1

          ! maximum distance to reference element
          if (dist_squared > dist_squared_max) dist_squared_max = dist_squared
        endif
      enddo

      ! statistics
      if (num_neighbor_neighbors > num_neighbor_neighbors_max) num_neighbor_neighbors_max = num_neighbor_neighbors
    endif

    ! statistics
    if (num_neighbors > num_neighbors_max) num_neighbors_max = num_neighbors
    if (ielem_counter > num_elements_actual_max) num_elements_actual_max = ielem_counter

    ! adjacency indexing
    neighbors_xadj(ispec_ref + 1) = inum_neighbor
    ! how to use:
    !num_neighbors = neighbors_xadj(ispec+1)-neighbors_xadj(ispec)
    !do i = 1,num_neighbors
    !  ! get neighbor
    !  ispec_neighbor = neighbors_adjncy(neighbors_xadj(ispec) + i)
    !enddo

    ! frees kdtree search array
    if (use_kdtree_search) then
      ! dynamic instead of setting fixed max_search_points
      deallocate(kdtree_search_index)
    endif

    ! user output progress
    if (myrank == 0) then
      if (mod(ispec_ref,max(NSPEC_AB/10,1)) == 0) then
        tCPU = wtime() - time1
        ! elapsed
        write(IMAIN,*) "    ",int(ispec_ref/(max(NSPEC_AB/10,1)) * 10)," %", &
                       " - elapsed time:",sngl(tCPU),"s"
        ! flushes file buffer for main output file (IMAIN)
        call flush_IMAIN()
      endif
    endif
  enddo ! ispec_ref

  ! frees temporary array
  deallocate(flag_topological)
  deallocate(mask_ispec)
  deallocate(ibool_corner)

  ! debug: for vtk output
  if (DEBUG_VTK_OUTPUT) then
    ! number of neighbors
    allocate(tmp_num_neighbors(NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating tmp_num_neighbors array'
    ! fills temporary array
    do ispec_ref = 1,NSPEC_AB
      ! gets number of neighbors
      num_neighbors = neighbors_xadj(ispec_ref+1) - neighbors_xadj(ispec_ref)
      tmp_num_neighbors(ispec_ref) = num_neighbors
    enddo
    filename = trim(prname) // 'mesh_neighbors'
    call write_VTK_data_elem_i(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,tmp_num_neighbors,filename)
    if (myrank == 0) then
      write(IMAIN,*) '     written file: ',trim(filename)//'.vtk'
      call flush_IMAIN()
    endif
    deallocate(tmp_num_neighbors)
  endif

  ! check if element has neighbors
  ! note: in case of a fault in this slice (splitting nodes) and/or scotch paritioning
  !       it can happen that an element has no neighbors
  if (NPROC == 1 .and. (.not. ANY_FAULT_IN_THIS_PROC)) then
    ! checks if neighbors were found
    do ispec_ref = 1,NSPEC_AB
      ! gets number of neighbors
      num_neighbors = neighbors_xadj(ispec_ref+1) - neighbors_xadj(ispec_ref)

      ! element should have neighbors, otherwise mesh is probably invalid
      if (num_neighbors == 0 .and. NSPEC_AB > 1) then
        ! midpoint
        iglob = ibool(MIDX,MIDY,MIDZ,ispec_ref)
        xyz_target(1) = xstore(iglob)
        xyz_target(2) = ystore(iglob)
        xyz_target(3) = zstore(iglob)
        ! error info
        print *,'Error: rank ',myrank,' - element ',ispec_ref,'has no neighbors!'
        print *,'  element midpoint location: ',xyz_target(:)
        print *,'  maximum element size     : ',elemsize_max_glob
        print *,'  typical search distance  : ',maximal_elem_size_squared
        print *,'  kdtree search            : ',use_kdtree_search
        if (use_kdtree_search) &
          print *,'  kd-tree r_search         : ',r_search
        print *,'  maximum search elements  : ',num_elements_max
        call exit_MPI(myrank,'Error adjacency invalid')
      endif
    enddo
  endif

  ! total number of neighbors
  num_neighbors_all = inum_neighbor

  ! allocates compacted array
  allocate(neighbors_adjncy(num_neighbors_all),stat=ier)
  if (ier /= 0) stop 'Error allocating neighbors_adjncy'

  neighbors_adjncy(1:num_neighbors_all) = tmp_adjncy(1:num_neighbors_all)

  ! checks
  if (minval(neighbors_adjncy(:)) < 1 .or. maxval(neighbors_adjncy(:)) > NSPEC_AB) stop 'Invalid adjncy array'

  ! frees temporary array
  deallocate(tmp_adjncy)
  deallocate(xyz_midpoints)

  if (use_kdtree_search) then
    ! frees current tree memory
    ! deletes tree arrays
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
    ! alternative: to avoid allocating/deallocating search index arrays by using max_search_points
    !deallocate(kdtree_search_index)
    ! deletes search tree nodes
    call kdtree_delete()
  endif

  ! checks neighbors for out-of-bounds indexing and duplicates
  do ispec_ref = 1,NSPEC_AB
    ! loops over neighbors
    num_neighbors = neighbors_xadj(ispec_ref+1) - neighbors_xadj(ispec_ref)
    do ii = 1,num_neighbors
      ! get neighbor entry
      ielem = neighbors_xadj(ispec_ref) + ii

      ! checks entry index
      if (ielem < 1 .or. ielem > num_neighbors_all) &
        stop 'Invalid element index in neighbors_xadj array'

      ! get neighbor element
      ispec = neighbors_adjncy(ielem)

      ! checks element index
      if (ispec < 1 .or. ispec > NSPEC_AB) &
        stop 'Invalid ispec index in neighbors_adjncy array'

      ! loops over all other neighbors
      do jj = ii+1,num_neighbors
        ! checks for duplicate
        if (neighbors_adjncy(neighbors_xadj(ispec_ref) + jj) == ispec) &
          stop 'Invalid ispec duplicate found in neighbors_adjncy array'
      enddo
    enddo
  enddo

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of neighbor detection
    tCPU = wtime() - time1
    write(IMAIN,*)
    write(IMAIN,*) '     maximum search elements                                      = ',num_elements_max
    write(IMAIN,*) '     maximum of actual search elements (after distance criterion) = ',num_elements_actual_max
    write(IMAIN,*)
    write(IMAIN,*) '     estimated maximum element size            = ',elemsize_max_glob
    write(IMAIN,*) '     maximum distance between neighbor centers = ',real(sqrt(dist_squared_max),kind=CUSTOM_REAL)
    if (use_kdtree_search) then
      if (sqrt(dist_squared_max) > r_search - 0.5*elemsize_max_glob) then
        write(IMAIN,*) '***'
        write(IMAIN,*) '*** Warning: consider increasing the kd-tree search radius to improve this neighbor setup ***'
        write(IMAIN,*) '***'
      endif
    endif
    write(IMAIN,*)
    write(IMAIN,*) '     maximum neighbors found per element = ',num_neighbors_max
    if (ADD_NEIGHBOR_OF_NEIGHBORS) then
      write(IMAIN,*) '         (maximum neighbor of neighbors) = ',num_neighbor_neighbors_max
    endif
    write(IMAIN,*) '     total number of neighbors           = ',num_neighbors_all
    write(IMAIN,*)
    write(IMAIN,*) '     Elapsed time for detection of neighbors in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine setup_mesh_adjacency
