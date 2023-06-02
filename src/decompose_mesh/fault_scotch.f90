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

module fault_scotch

  use constants, only: MAX_STRING_LEN, TINYVAL, HUGEVAL, IN_DATA_FILES, NDIM
  use shared_parameters, only: NGNOD2D

  ! LTS mode
  use decompose_mesh_par, only: LTS_MODE, num_p_level, p_level, ispec_p_refine

  ! kd-tree for searching nearest neighbors
  use kdtree_search, only: &
    kdtree_num_nodes, kdtree_nodes_location, kdtree_nodes_index, &
    kdtree_setup, kdtree_find_nearest_neighbor, kdtree_delete

  implicit none

  private

  type fault_type
    private
    integer :: nspec
    integer, dimension(:), pointer  :: ispec1, ispec2 !, iface1, iface2
    integer, dimension(:,:), pointer  :: inodes1, inodes2
  end type fault_type

  type(fault_type), allocatable, save :: faults(:)
  double precision, dimension(:,:), allocatable, save :: nodes_coords_open

  ! fault flag (set to .true. in case Par_file_faults defines faults)
  logical, save :: ANY_FAULT = .false.

  !unused so far...
  !integer, parameter :: long = SELECTED_INT_KIND(18)

  ! tolerance of split node distance:
  !   must be larger than the fault offset in the mesh,
  !   but smaller than the smallest element size
  ! relative to minium element size on fault surface
  ! (for example: 0.01 -> 1% of minimum fault surface element size)
  double precision, parameter :: FAULT_GAP_TOLERANCE = 0.01d0

  !------------------------------------------------------------------------
  ! parallel faults

  ! for parallel faults, chooses a greedy scheme (.true.) or a heuristic balancing (.false.)
  logical,parameter :: FAULT_BALANCE_GREEDY = .false.

  ! adds elements only touching fault by a node/edge to fault partitions
  logical,parameter :: ADD_FAULT_TOUCH_ELEMENTS = .true.

  ! node search algorithm (on opposite fault surface)
  logical,parameter :: FAULT_NODE_SEARCH_BY_KDTREE = .true.

  !------------------------------------------------------------------------

  public :: read_fault_files, fault_repartition, close_faults, write_fault_database, &
            write_fault_database_mpi, save_nodes_coords, nodes_coords_open, ANY_FAULT, bcast_faults_mpi

CONTAINS

!==========================================================================================

  subroutine read_fault_files(localpath_name)

  use constants, only: IIN_PAR

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: localpath_name

  ! local parameters
  integer :: nbfaults, iflt, ier
  character(len=MAX_STRING_LEN) :: filename

  filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES)) // 'Par_file_faults'

  open(unit=IIN_PAR,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier == 0) then
    read(IIN_PAR,*) nbfaults
  else
    nbfaults = 0
    print *
    print *, 'Par_file_faults not found: assuming that there are no faults'
    print *
  endif
  close(IIN_PAR)

  ANY_FAULT = (nbfaults > 0)
  if (.not. ANY_FAULT)  return

  ! user output
  print *,'faults : '
  print *,'  number of faults = ',nbfaults

  allocate(faults(nbfaults),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
  ! NOTE: asumes that the fault ids follow a contiguous numbering, starting at 1, with unit increment
  !       The user must assign that numbering during mesh generation
  do iflt = 1 , nbfaults
    call read_single_fault_file(faults(iflt),iflt,localpath_name)
  enddo

  end subroutine read_fault_files


!---------------------------------------------------------------------------------------------------

  subroutine read_single_fault_file(f,ifault,localpath_name)

  use constants, only: IIN_FLT

  implicit none
  type(fault_type), intent(inout) :: f
  character(len=MAX_STRING_LEN), intent(in) :: localpath_name

  character(len=MAX_STRING_LEN) :: filename
  integer,intent(in) :: ifault
  character(len=5) :: NTchar
  integer :: e,ier,nspec_side1,nspec_side2

  write(NTchar,'(I5)') ifault
  NTchar = adjustl(NTchar)

  filename = localpath_name(1:len_trim(localpath_name))//'/fault_file_'//&
             NTchar(1:len_trim(NTchar))//'.dat'

  ! trim left whitespace
  filename = adjustl(filename)

  ! reads fault elements and nodes
 ! File format:
 ! Line 1:
 !   number_of_elements_in_side_1   number_of_elements_in_side_2
 ! Then for all elements that have a face on side 1:
 !   #id_element #id_global_node1 .. #id_global_node4
 ! Then the same for side 2.
 ! Note: element ids start at 1, not 0 (see cubit2specfem3d.py)
  open(unit=IIN_FLT, file=trim(filename), status='old', form='formatted', iostat = ier)
  if (ier /= 0) then
    write(*,*) 'Fatal error: file '//trim(filename)//' not found'
    write(*,*) 'Abort'
    stop
  endif

  read(IIN_FLT,*) nspec_side1,nspec_side2
  if (nspec_side1 /= nspec_side2) stop 'Number of fault nodes at do not match.'
  f%nspec = nspec_side1
  allocate(f%ispec1(f%nspec), &
           f%ispec2(f%nspec), &
           f%inodes1(NGNOD2D,f%nspec), &
           f%inodes2(NGNOD2D,f%nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 82')

  ! reads in element ids and element surface node ids (indices starting at 1)
  do e = 1,f%nspec
    read(IIN_FLT,*) f%ispec1(e),f%inodes1(:,e)
  enddo
  do e = 1,f%nspec
    read(IIN_FLT,*) f%ispec2(e),f%inodes2(:,e)
  enddo

  ! If we ever figure out how to export "ifaces" from CUBIT:
  !allocate(f%iface1(f%nspec))
  !allocate(f%iface2(f%nspec))
  !do e=1,f%nspec
  !  read(IIN_FLT,*) f%ispec1(e),f%ispec2(e),f%iface1(e),f%iface2(e)
  !enddo

  close(IIN_FLT)

  print *,'  reading fault: elements =',f%nspec

  end subroutine read_single_fault_file


! ---------------------------------------------------------------------------------------------------

! Saving nodes_coords to be used in the solver for ibool_fault_side1 and side2

  subroutine save_nodes_coords(nodes_coords,nnodes)

  implicit none
  integer, intent(in) :: nnodes
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords
  integer :: ier

  allocate(nodes_coords_open(NDIM,nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 83')
  nodes_coords_open(:,:) = nodes_coords(:,:)

  end subroutine save_nodes_coords

! ---------------------------------------------------------------------------------------------------

  subroutine close_faults(nodes_coords,nnodes)

  implicit none
  integer, intent(in) :: nnodes
  double precision, dimension(NDIM,nnodes), intent(inout) :: nodes_coords

  integer  :: i

  do i = 1,size(faults)
    call close_fault_single(faults(i),nodes_coords,nnodes)
  enddo

  end subroutine close_faults

! ---------------------------------------------------------------------------------------------------

!jpa: to do this much faster:
!     1. create a list of unique nodes from inodes1 and inodes2,
!          inodes1_u = unique(isort1)
!          inodes2_u = unique(isort2)
!     2. sort the nodes by coordinates. Now both faces correspond.
!          [coord,k1] = sort(coords(inodes1_u))
!          k1 = inodes1_u(k1)
!          [coord,k2] = sort(coords(inodes2_u))
!          k2 = inodes2_u(k2)
!     3. set the coordinates on both sides equal to their average
!          coords(k1) = 0.5*( coords(k1)+coords(k2) )
!          coords(k2) = coords(k1)
!
! daniel: not sure if above suggestion would work:
!         a) if we are sure that all fault nodes are correct
!            and have their corresponding node on the opposite fault, this might work.
!            it would imply that meshing both surfaces was done correctly, i.e. by using
!            the same meshing size and meshing algorithm type, and finally hoping that both surfaces are
!            meshed in identical ways since their geometries are only slightly different;
!         b) the lexicographically sorting sorts in increasing x-, y- then z-direction; since we are "opening"
!            the two fault surfaces by a small amount in an (unspecified) fault direction, this sorting
!            order might work sometimes, sometimes it might shuffle around points (when crossing zero-origin)
!         -> maybe better check the minimum distances between matching nodes again to be sure

  subroutine close_fault_single(f,nodes_coords,nnodes)

  implicit none
  type(fault_type), intent(in) :: f
  integer, intent(in)  :: nnodes
  double precision, dimension(NDIM,nnodes), intent(inout) :: nodes_coords

  double precision, dimension(NDIM) :: xyz_1, xyz_2, xyz

  double precision :: dist
  double precision :: dist_min,dist_min_glob,elem_size_min
  integer :: iglob1, iglob2, i, j, k1, k2
  integer :: iglob_min
  logical :: found_it
  logical,dimension(:),allocatable :: mask_iglob

  ! kd-tree search
  integer :: iloc,ier

  ! timing
  real :: ct_start,ct_end

  print *,'  closing fault: elements =',f%nspec

  ! gets minimum fault surface element size
  dist_min = HUGEVAL
  do i = 1,f%nspec
    ! loops over corners of surface element on fault surface 1
    do k2 = 1,NGNOD2D
      iglob2 = f%inodes1(k2,i)
      xyz_2 = nodes_coords(:,iglob2)
      ! distances of this reference corner to all other corners
      do k1 = k2+1,NGNOD2D
        iglob1 = f%inodes1(k1,i)
        xyz_1 = nodes_coords(:,iglob1)
        xyz = xyz_2-xyz_1
        dist = xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)  ! distance squared - to be faster
        if (dist < dist_min) then
          dist_min = dist
        endif
      enddo
    enddo
  enddo
  dist_min = sqrt(dist_min)
  print *,'  fault-up   surface: minimum element size =',dist_min

  ! on opposite fault surface
  dist_min = HUGEVAL
  do i = 1,f%nspec
    ! loops over corners of surface element on fault surface 2
    do k2 = 1,NGNOD2D
      iglob2 = f%inodes2(k2,i)
      xyz_2 = nodes_coords(:,iglob2)
      ! distances of this reference corner to all other corners
      do k1 = k2+1,NGNOD2D
        iglob1 = f%inodes2(k1,i)
        xyz_1 = nodes_coords(:,iglob1)
        xyz = xyz_2-xyz_1
        dist = xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)  ! distance squared - to be faster
        if (dist < dist_min) then
          dist_min = dist
        endif
      enddo
    enddo
  enddo
  dist_min = sqrt(dist_min)
  print *,'  fault-down surface: minimum element size =',dist_min

  ! sets minimum size of fault surface elements
  elem_size_min = dist_min

  ! user output
  if (FAULT_NODE_SEARCH_BY_KDTREE) then
    print *,'  using kd-tree search for matching fault nodes'
  else
    print *,'  using brute-force search for matching fault nodes'
  endif

  ! timing
  call cpu_time(ct_start)

  ! allocates mask for iglob values
  allocate(mask_iglob(nnodes),stat=ier)
  if (ier /= 0) stop 'Error allocating array mask_iglob'
  mask_iglob(:) = .false.

  ! prepares search arrays for opposite fault
  if (FAULT_NODE_SEARCH_BY_KDTREE) then
    ! kd-tree search
    ! counts number of global points on fault
    mask_iglob(:) = .false.
    iloc = 0
    do i = 1,f%nspec
      do k1 = 1,NGNOD2D
        iglob1 = f%inodes1(k1,i)
        ! checks if global point already added
        if (mask_iglob(iglob1) .eqv. .true.) cycle
        ! new point
        iloc = iloc + 1
        mask_iglob(iglob1) = .true.
      enddo
    enddo
    if (iloc > f%nspec * NGNOD2D) stop 'Error index iloc does not match total number'
    if (iloc < 1) stop 'Error index iloc is zero'

    ! set number of tree nodes
    kdtree_num_nodes = iloc

    ! allocates tree arrays
    allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2096')
    if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
    kdtree_nodes_location(:,:) = 0.0

    allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2097')
    if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'
    kdtree_nodes_index(:) = 0

    ! prepares search arrays
    ! fills kd-tree arrays
    mask_iglob(:) = .false.
    iloc = 0
    do i = 1,f%nspec
      do k1 = 1,NGNOD2D
        iglob1 = f%inodes1(k1,i)
        ! checks if global point already added
        if (mask_iglob(iglob1) .eqv. .true.) cycle

        ! new point
        iloc = iloc + 1
        if (iloc > kdtree_num_nodes) stop 'Error index iloc bigger than kdtree_num_nodes'
        mask_iglob(iglob1) = .true.

        ! adds point
        kdtree_nodes_index(iloc) = iglob1
        kdtree_nodes_location(:,iloc) = nodes_coords(:,iglob1)
      enddo
    enddo
    if (iloc /= kdtree_num_nodes) stop 'Error index iloc does not match kdtree_num_nodes'

    ! creates kd-tree for searching
    call kdtree_setup()

  endif

  ! searches closest matching node on opposite fault surface and merges coordinates
  dist_min_glob = -1.0d0
  mask_iglob(:) = .false.

  do i = 1,f%nspec
    do k2 = 1,NGNOD2D
      iglob2 = f%inodes2(k2,i)

      ! checks if this node has been considered already; if so, it goes to the next corner
      ! note: we loop over corners of all surface elements, i.e. some elements could share corners
      if (mask_iglob(iglob2) .eqv. .true.) cycle

      xyz_2(:) = nodes_coords(:,iglob2)

      found_it = .false.

      ! search matching node on opposite fault surface
      if (FAULT_NODE_SEARCH_BY_KDTREE) then
        ! kd-tree search
        ! finds closest point in opposite fault (inodes1)
        call kdtree_find_nearest_neighbor(xyz_2,iglob_min,dist_min)

        ! merges closest points coordinates
        if (dist_min <= FAULT_GAP_TOLERANCE * elem_size_min) then
          iglob1 = iglob_min
          xyz_1 = nodes_coords(:,iglob1)
          xyz(:) = (xyz_1(:) + xyz_2(:)) * 0.5d0
          nodes_coords(:,iglob2) = xyz
          nodes_coords(:,iglob1) = xyz
          found_it = .true.
          ! sets flag that node has been merged and will be treated only once
          mask_iglob(iglob2) = .true.
        else
          ! user output
          print *,'warning: point ',iglob2,' - minimum distance to closest node exceeds tolerance:', &
                 dist_min,FAULT_GAP_TOLERANCE * elem_size_min
          print *
        endif

      else
        ! brute-force search
        do j = 1,f%nspec
          do k1 = 1,NGNOD2D
            iglob1 = f%inodes1(k1,j)
            xyz_1(:) = nodes_coords(:,iglob1)
            xyz(:) = xyz_2(:) - xyz_1(:)
            dist_min = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))

            !jpa: Closing nodes that are already closed is not a problem
            !jpa: I process them again to leave the loop as early as possible
            !jpa: and to test if a facing node was found (see below).

            if (dist_min <= FAULT_GAP_TOLERANCE * elem_size_min) then
              ! nodes match and are close enough
              ! closing split node gaps by assigning the mid-point location between the two nodes
              xyz(:) = (xyz_1(:) + xyz_2(:)) * 0.5d0
              nodes_coords(:,iglob2) = xyz(:)
              nodes_coords(:,iglob1) = xyz(:)
              found_it = .true.
              ! sets flag that node has been merged and will be treated only once
              mask_iglob(iglob2) = .true.
              ! done
              exit
            endif
          enddo
          if (found_it) exit
        enddo

      endif

      ! largest minimum
      if (dist_min > dist_min_glob) dist_min_glob = dist_min

      ! jpa: If the two fault sides have been meshed independently they might not match. Test it here:
      if (.not. found_it) then
        print *,'Error: point not found in element',i,'of',f%nspec,' with iglob',iglob2
        print *,'  point x/y/z:',nodes_coords(:,iglob2)
        print *,'  found x/y/z:',nodes_coords(:,iglob_min)
        print *,'  minimum distance = ',sngl(dist_min),' is bigger than tolerance ',FAULT_GAP_TOLERANCE * elem_size_min
        print *,'       with tolerance ',FAULT_GAP_TOLERANCE,'of minimum element size ',elem_size_min
        stop 'Inconsistent fault mesh: corresponding node in the other fault face was not found'
      endif

    enddo
  enddo

  ! timing
  call cpu_time(ct_end)
  print *,'  search timing: ',ct_end - ct_start, '(s)'
  print *,'  opposite node maximum gap distance =',dist_min_glob

  ! deletes search tree
  if (FAULT_NODE_SEARCH_BY_KDTREE) then
    ! deletes tree arrays
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
    ! deletes search tree nodes
    call kdtree_delete()
  endif

  end subroutine close_fault_single

!===================================================================================================
! Lexicographic reordering of fault elements (based on their centroid)
! to make sure both sides are ordered in the same way
! and hence elements facing each other have the same index

  subroutine reorder_fault_elements(nodes_coords,nnodes)

  implicit none
  integer, intent(in)  :: nnodes
  double precision,dimension(NDIM,nnodes), intent(in) :: nodes_coords

  integer :: i

  do i = 1,size(faults)
    call reorder_fault_elements_single(faults(i),nodes_coords,nnodes)
  enddo

  end subroutine reorder_fault_elements

! ---------------------------------------------------------------------------------------------------

  subroutine reorder_fault_elements_single(f,nodes_coords,nnodes)

  implicit none
  type(fault_type), intent(inout) :: f
  integer, intent(in) :: nnodes
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords

  double precision, dimension(NDIM,f%nspec) :: xyz_c
  integer :: k,e,iglob
  integer, dimension(f%nspec) :: new_index_list

  ! compute element-face centers for fault side 1
  xyz_c(:,:) = 0.d0
  do e = 1,f%nspec
    do k = 1,4
      iglob = f%inodes1(k,e)
      xyz_c(:,e) = xyz_c(:,e) + nodes_coords(:,iglob)
    enddo
  enddo
  xyz_c(:,:) = 0.25d0 * xyz_c(:,:)

  ! reorder
  call lex_order(xyz_c,new_index_list,f%nspec)

  f%ispec1(:) = f%ispec1(new_index_list(:))
  f%inodes1(:,:) = f%inodes1(:,new_index_list(:))

  ! repeat for fault side 2
  xyz_c(:,:) = 0.d0
  do e = 1,f%nspec
    do k = 1,4
      iglob = f%inodes2(k,e)
      xyz_c(:,e) = xyz_c(:,e) + nodes_coords(:,iglob)
    enddo
  enddo
  xyz_c(:,:) = 0.25d0 * xyz_c(:,:)

  ! reorder
  call lex_order(xyz_c,new_index_list,f%nspec)

  f%ispec2(:) = f%ispec2(new_index_list(:))
  f%inodes2(:,:) = f%inodes2(:,new_index_list(:))

  end subroutine reorder_fault_elements_single

! ---------------------------------------------------------------------------------------------------

  subroutine lex_order(xyz_c,locval,nspec)

  implicit none
  integer, intent(in) :: nspec
  integer, intent(out) :: locval(nspec)
  double precision, intent(in) :: xyz_c(3,nspec)

  double precision, dimension(nspec) :: xp,yp,zp
  integer, dimension(nspec) :: ninseg,ibool,iglob
  integer :: nglob
  logical :: ifseg(nspec)
  double precision :: xtol
  double precision :: maxdist,dist

  xp(:) = xyz_c(1,:)
  yp(:) = xyz_c(2,:)
  zp(:) = xyz_c(3,:)

  ! define geometrical tolerance based upon typical size of the model
  !xtol = 1.d-10 * maxval( maxval(xyz_c,2) - minval(xyz_c,2) )

  maxdist = TINYVAL
  ! x-direction
  dist = maxval(xp(:)) - minval(xp(:))
  if (dist > maxdist) maxdist = dist
  ! y-direction
  dist = maxval(yp(:)) - minval(yp(:))
  if (dist > maxdist) maxdist = dist
  ! z-direction
  dist = maxval(zp(:)) - minval(zp(:))
  if (dist > maxdist) maxdist = dist

  xtol = 1.d-10 * maxdist


  call sort_array_coordinates(nspec,xp,yp,zp,ibool,iglob,locval,ifseg, &
                              nglob,ninseg,xtol)

  end subroutine lex_order

!===================================================================================================
  !--------------------------------------------------
  ! Repartitioning : two coupled elements on fault side1/side2 are transfered to the same partition
  !--------------------------------------------------

  subroutine fault_repartition(nelmnts, nnodes, elmnts, nsize, nproc, part, NGNOD, nodes_coords, elmnts_load)

  use constants, only: PARALLEL_FAULT
  use shared_parameters, only: LOCAL_PATH

  implicit none
  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, NGNOD
  integer, dimension(0:NGNOD*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords
  integer, dimension(0:nelmnts-1), intent(in)       :: elmnts_load

  ! debug vtk-output
  integer,dimension(:,:),allocatable :: tmp_elmnts
  integer, dimension(:),allocatable :: tmp_part
  double precision,dimension(:),allocatable :: xstore_dummy,ystore_dummy,zstore_dummy
  character(len=MAX_STRING_LEN) :: filename
  integer :: iflt,i,e,e1,e2,iglob

  ! debugging: outputs vtk-file with fault elements
  logical,parameter :: DEBUG_VTK = .false.

  if (PARALLEL_FAULT) then
    ! fault can be parallelized
    call fault_repartition_parallel(nelmnts,part, nodes_coords, nnodes, nproc, elmnts_load, elmnts, NGNOD, nsize)
  else
    ! moves all fault elements to the same partition (proc=0)
    call fault_repartition_not_parallel(nelmnts, nnodes, elmnts, nsize, nproc, part, NGNOD)
  endif

  ! debug: for vtk output
  if (DEBUG_VTK) then
    ! temporary arrays
    allocate(tmp_elmnts(NGNOD,nelmnts))
    do e = 1,nelmnts
      do i = 1,NGNOD
        iglob = (i-1) + (e-1)*NGNOD
        tmp_elmnts(i,e) = elmnts(iglob)+1
      enddo
    enddo

    allocate(xstore_dummy(nnodes),ystore_dummy(nnodes),zstore_dummy(nnodes))
    xstore_dummy(:) = nodes_coords(1,:)
    ystore_dummy(:) = nodes_coords(2,:)
    zstore_dummy(:) = nodes_coords(3,:)

    allocate(tmp_part(0:nelmnts-1))
    tmp_part(:) = -1

    ! fills in actual partition number for fault elements
    do iflt = 1,size(faults)
      do e = 1,faults(iflt)%nspec
        e1 = faults(iflt)%ispec1(e) - 1
        e2 = faults(iflt)%ispec2(e) - 1
        tmp_part(e1) = part(e1)
        tmp_part(e2) = part(e2)
      enddo
    enddo

    ! outputs part array valid only on fault elements
    filename = trim(LOCAL_PATH) // '/fault_partitioning_part_array'
    call write_VTK_data_ngnod_elem_i(nelmnts,nnodes,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
                                     tmp_elmnts,tmp_part,filename)
    print *,'  written file: ',trim(filename)//'.vtk'
    print *

    ! outputs final part array - same as part_array.vtk output in decompose_mesh.F90
    !filename = trim(LOCAL_PATH) // '/partitioning_part_array'
    !call write_VTK_data_ngnod_elem_i(nelmnts,nnodes,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
    !                             tmp_elmnts,part,filename)
    !print *,'  written file: ',trim(filename)//'.vtk'
    !print *

    ! frees temporary arrays
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
    deallocate(tmp_elmnts)
    deallocate(tmp_part)
  endif

  end subroutine fault_repartition

!---------------------------------------------------------------------------------------------------

!     part(e) = index of the processor to which element #e is assigned.
!               Fault elements and neighbors are assigned to the same processor.
!               Part, once modified, will be input for write_partition_database.

  subroutine fault_repartition_not_parallel(nelmnts, nnodes, elmnts, nsize, nproc, part, NGNOD)

  implicit none
  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, NGNOD
  integer, dimension(0:NGNOD*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part

  ! local parameters
  integer, dimension(0:nnodes-1) :: nnodes_elmnts
  integer, dimension(0:nsize*nnodes-1) :: nodes_elmnts
  integer  :: i,j,ipart,nproc_null,nproc_null_final
  integer  :: k1, k2, k,e,iflt,inode,ier
  integer, dimension(:), allocatable :: elem_proc_null
  integer, dimension(:),allocatable :: part_elem_fault
  integer :: nelem_fault,nelem_fault_added
  integer :: icounter,ielem,ielem_added

  ! user output
  print *,'fault repartitioning: (non-parallel version)'

  ! downloading processor 0
  nproc_null = count( part(:) == 0 )
  print *,'  elements in slice proc 0 redistributed in [{nproc}- nproc0] :',nproc_null

  ! Fault zone layer = the set of elements that contain at least one fault node
  print *,"  fault zone layer :"
  print *, "     ",size(faults),"number of faults"

  ! for each element touching fault, lists its partition
  allocate(part_elem_fault(0:nelmnts-1))

  ! List of elements per node
  !  nnodes_elmnts(i) = number of elements containing node #i (i=0:nnodes-1)
  !  nodes_elmnts(nsize*i:nsize*i+nnodes_elmnts(i)-1) = index of elements (starting at 0) containing node #i
  !  nsize = maximun number of elements in a node.
  !  NGNOD = nodes per element.
  nnodes_elmnts(:) = 0
  nodes_elmnts(:)  = 0
  do i = 0, NGNOD*nelmnts-1
    nodes_elmnts(elmnts(i)*nsize + nnodes_elmnts(elmnts(i))) = i/NGNOD
    nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
  enddo

  ! counts all elements touching faults
  part_elem_fault(:) = -1
  do iflt = 1,size(faults)
    do e = 1,faults(iflt)%nspec
      do k = 1,4
        ! fault nodes 1
        inode = faults(iflt)%inodes1(k,e)-1  ! node index, starting at 0
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part_elem_fault( nodes_elmnts(k1:k2) ) = part( nodes_elmnts(k1:k2) )
        ! fault nodes 2
        inode = faults(iflt)%inodes2(k,e)-1
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part_elem_fault( nodes_elmnts(k1:k2) ) = part( nodes_elmnts(k1:k2) )
      enddo
    enddo
  enddo

  ! user output
  print *,"     initial fault element distribution:"
  do i = 0,nproc-1
    print *,"     proc ",i,": fault elements = ",count(part_elem_fault(:) == i)," out of ",count(part(:) == i)
  enddo
  nelem_fault = count( part_elem_fault(:) >= 0 )
  print *,"     total number of fault elements: ",nelem_fault

  ! number of fault elements which will be added to process 0
  nelem_fault_added = nelem_fault - count( part_elem_fault(:) == 0 )

  ! checks if anything to do
  if (nproc == 1) return

  ! distributs elements in slice 0 to higher procs
  if (nproc_null /= 0) then
    allocate(elem_proc_null(nproc_null),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 84')
    ! Filling up proc = 0 elements
    nproc_null = 0
    do i = 0,nelmnts-1
      if (part(i) == 0) then
        nproc_null = nproc_null + 1
        elem_proc_null(nproc_null) = i
      endif
    enddo

    if (.true.) then
      ! Redistributing proc-0 elements on the rest of processors
      ipart = 1
      do i = 1, nproc_null
        part(elem_proc_null(i)) = ipart
        if (ipart == nproc-1) ipart = 0
        ipart = ipart +1
      enddo
    else
      ! Redistributing same number of proc-0 elements as total fault elements on the rest of processors
      ! select first next partition to put elements into
      ipart = 1
      do i = 1,nproc-1
        icounter = count( part_elem_fault(:) == i )
        if (icounter > 0) then
          ipart = i
          exit
        endif
      enddo

      ! loops over all elements in process 0
      icounter = 0
      ielem_added = 0
      do i = 1, nproc_null
        ! starts with last element in process 0
        ielem = elem_proc_null(nproc_null - i + 1)

        ! skips if this is already a fault element
        if (part_elem_fault(ielem) == 0) cycle

        ! re-assigns element to next partition
        part(ielem) = ipart

        ! updates counters
        icounter = icounter + 1
        ielem_added = ielem_added + 1

        ! stops loop after counter reached number of fault elements which will be put new into slice of process 0
        if (icounter >= nelem_fault_added) exit

        ! sets next partition number (if bin is full)
        !if (icounter > 1 .and. mod(icounter-1 , nelem_fault/(nproc-1) ) == 0) then
        !  if (ipart == nproc-1) ipart = 0
        !  ipart = ipart + 1
        !  !print *,"counter",icounter,"with bin:",nelem_fault/(nproc-1),"  next partition:",ipart
        !endif

        ! sets next partition if all fault elements in ipart-partition have been replaced
        if (ielem_added >= count(part_elem_fault(:) == ipart)) then
          ! select next partition to put elements into
          do j = ipart+1,nproc-1
            if (count( part_elem_fault(:) == j ) > 0) then
              ipart = j
              exit
            endif
          enddo
          !print *,"counter",icounter,"with new bin:",count(part_elem_fault(:) == ipart ),"  next partition:",ipart
          ! resets elem counter
          ielem_added = 0
        endif
      enddo
      print *,"     re-assigned elements:",icounter
    endif

    deallocate(elem_proc_null)
  endif

  ! repartitions elements
  ! puts all fault elements into process 0 slice

  ! counts all elements touching faults
  part_elem_fault(:) = -1

  do iflt = 1,size(faults)
    do e = 1,faults(iflt)%nspec
      do k = 1,4
        ! 1. fault side
        inode = faults(iflt)%inodes1(k,e)-1  ! node index, starting at 0
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        ! moves elements to proc0
        part( nodes_elmnts(k1:k2) ) = 0
        part_elem_fault( nodes_elmnts(k1:k2) ) = part( nodes_elmnts(k1:k2) )
        ! 2. fault side
        inode = faults(iflt)%inodes2(k,e)-1
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        ! moves elements to proc0
        part( nodes_elmnts(k1:k2) ) = 0
        part_elem_fault( nodes_elmnts(k1:k2) ) = part( nodes_elmnts(k1:k2) )
      enddo
    enddo
  enddo

  nproc_null_final = count( part(:) == 0 )

  ! user output
  print *,"     final fault element distribution:"
  do i = 0,nproc-1
    print *,"     proc ",i,": fault elements = ",count(part_elem_fault(:) == i),"out of",count(part(:) == i)
  enddo
  print *,"     ",nproc_null_final,"final elements in slice of process 0"
  print *

  deallocate(part_elem_fault)

  end subroutine fault_repartition_not_parallel

! ---------------------------------------------------------------------------------------------------

  subroutine fault_repartition_parallel(nelmnts, part, nodes_coords, nnodes, nproc, elmnts_load, elmnts, NGNOD, nsize)

  implicit none
  integer, intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(inout) :: part
  integer, intent(in) :: nnodes,nproc
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords

  integer, dimension(0:nelmnts-1), intent(in) :: elmnts_load ! note: the array bounds have changed from 0 to nelmnts-1
  integer, intent(in) :: NGNOD
  integer, dimension(NGNOD,0:nelmnts-1), intent(in) :: elmnts
  integer, intent(in) :: nsize

  ! local parameters
  integer :: i,iflt,e,e1,e2,e3
  integer :: proc1,proc2,proc,proc_new
  integer :: counter,num_fault_elem

  integer, dimension(:), allocatable :: part_bin,part_bin0
  integer, dimension(:), allocatable :: part_load,part_load0
  integer :: maxproc

  ! elements connected on fault surfaces
  integer, dimension(:,:), allocatable :: fault_elements_connected

  logical :: found

  ! adds elements only touching fault by a node/edge to fault partitions
  integer :: iglob,iglob1,ispec,inode,ier
  logical,dimension(:),allocatable :: mask_ispec
  integer :: nloop
  ! element list
  integer, dimension(:),allocatable :: nnodes_elmnts
  integer, dimension(:),allocatable :: nodes_elmnts
  integer :: k,k1,ndiff
  logical,dimension(:),allocatable :: mask_iglob

  ! LTS
  integer, dimension(:,:),allocatable :: num_pelem_part
  integer, dimension(:),allocatable :: p_lookup
  integer :: p,ilevel

  ! user output
  print *,'fault repartitioning: (parallel version)'

  ! checks if anything to do
  if (size(faults) < 1) return

  ! Reorder both fault sides so that elements facing each other have the same index
  call reorder_fault_elements(nodes_coords,nnodes)

  num_fault_elem = 0
  do iflt = 1,size(faults)
    num_fault_elem = num_fault_elem + faults(iflt)%nspec
  enddo
  print *,'  total number of fault elements: ',num_fault_elem
  print *

  ! checks if anything to do
  if (nproc == 1) return

  ! bins holding number of elements per partition
  maxproc = maxval( part(:) )
  allocate( part_bin(0:maxproc),part_bin0(0:maxproc),stat=ier )
  if (ier /= 0) stop 'Error allocating arrays part_bin,...'
  allocate( part_load(0:maxproc),part_load0(0:maxproc),stat=ier )
  if (ier /= 0) stop 'Error allocating arrays part_bin,...'

  part_bin(:) = 0
  part_bin0(:) = 0
  part_load(:) = 0
  part_load0(:) = 0
  do e = 0,nelmnts-1
    proc1 = part(e)
    if (proc1 < 0 .or. proc1 > maxproc) stop 'Error fault repartition in parallel, process id on element invalid'
    ! counts number of elements for this partition
    part_bin(proc1) = part_bin(proc1) + 1
    ! counts total load for this partition
    part_load(proc1) = part_load(proc1) + elmnts_load(e)
  enddo

  ! stores initial distribution (elements/loads)
  part_bin0(:) = part_bin(:)
  part_load0(:) = part_load(:)

  ! checks
  if (sum(part_bin0(:)) /= nelmnts) then
    print *,'Error: parititioned elements not equal: ','partitioned = ',sum(part_bin0(:)),'total = ',sum(part_bin(:))
    stop 'Error fault repartition in parallel, number of elements in partitions invalid'
  endif

  ! LTS balancing of fault element for each p-level
  if (LTS_MODE) then
    ! array for number of p-elements per p-level
    allocate(num_pelem_part(0:maxproc,num_p_level), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating num_pelem_part in fault repartitioning'
    num_pelem_part(:,:) = 0

    ! p_lookup maps p -> ilevel
    allocate(p_lookup(maxval(p_level(:))),stat=ier)
    if (ier /= 0) stop 'Error allocating p_lookup in fault repartitioning'
    p_lookup(:) = 0
    do ilevel = 1,num_p_level
      p = p_level(ilevel)
      p_lookup(p) = ilevel
    enddo

    ! counts number of p-elements per partition
    do e = 0,nelmnts-1
      ! assigned partition
      proc1 = part(e)
      ! gets p from ispec_p_refine (starts at 1)
      p = ispec_p_refine(e+1)
      ! gets p-level index
      ilevel = p_lookup(p)
      ! counts number of p-elements
      num_pelem_part(proc1,ilevel) = num_pelem_part(proc1,ilevel) + 1
    enddo
  endif

  ! mask
  allocate(mask_ispec(0:nelmnts-1),stat=ier)
  if (ier /= 0) stop 'Error allocating ispec_is_done'

  mask_ispec(:) = .false.
  part_bin(:) = 0

  ! array to hold list of connected elements
  allocate(fault_elements_connected(size(faults),0:nelmnts-1),stat=ier)
  if (ier /= 0) stop 'Error allocating array fault_elements_connected'
  fault_elements_connected(:,:) = -1

  ! setup connected elements by faults, to handle triple junctions (intersections between two faults)
  call setup_connected_elements(nelmnts,mask_ispec,fault_elements_connected)

  ! repartitions elements

  ! way 1: put elements to minimum partition number
  if (FAULT_BALANCE_GREEDY) then
    !JPA loop over all faults
    !JPA loop over all fault element pairs
    !JPA assign both elements to the processor with lowest rank among the pair

    ! user output
    print *,'  greedy scheme'

    !JPA do it twice, to handle triple junctions (intersections between two faults)
    if (size(faults) > 1) then
      nloop = size(faults)
    else
      nloop = 1
    endif
    do i = 1,nloop
      counter = 0
      do iflt = 1,size(faults)
        do e = 1,faults(iflt)%nspec
          e1 = faults(iflt)%ispec1(e) - 1  ! part array between [0,nelmnts-1]
          e2 = faults(iflt)%ispec2(e) - 1
          proc1 = part(e1)
          proc2 = part(e2)

          ! checks if we need to re-assign element
          if (proc1 == proc2) cycle

          ! moves elements to lower partition of the two
          proc = min(proc1,proc2)
          part(e1) = proc
          part(e2) = proc

          ! statistics
          part_bin(proc) = part_bin(proc) + 1
          ! load change
          if (proc == proc1) then
            part_load(proc1) = part_load(proc1) + elmnts_load(e2)
            part_load(proc2) = part_load(proc2) - elmnts_load(e2)
          else
            part_load(proc2) = part_load(proc2) + elmnts_load(e1)
            part_load(proc1) = part_load(proc1) - elmnts_load(e1)
          endif
          counter = counter + 1
        enddo
      enddo
      print *,'  repartition loop',i,' changed ',counter,' coupled fault elements out of ',num_fault_elem
    enddo
    print *

  else
    ! way 2:
    ! tries to balance out fault elements

    ! user output
    print *,'  heuristic scheme'

    ! only need to loop once
    do e1 = 0,nelmnts-1
      ! only looks at reference elements
      if (.not. mask_ispec(e1)) cycle

      ! proc of reference element
      proc1 = part(e1)
      proc_new = proc1

      ! flag used to decide if we have to re-assign partitions
      found = .true.
      do iflt=1,size(faults)
        e2 = fault_elements_connected(iflt,e1)
        if (e2 >= 0) then
          proc2 = part(e2)
          ! checks if we need to re-assign element
          if (proc1 /= proc2) then
            found = .false.

            ! chooses partition
            ! default: chooses minimum partition
            proc_new = min(proc_new,proc2)

            ! selects partition with less elements to "balance" number of elements within 1%
            ! note: this "balancing" could be improved...
            if (LTS_MODE) then
              ! LTS: balances partitions within each p-level
              !
              ! gets p-level
              ilevel = p_lookup(ispec_p_refine(e2+1))
              ! sets new reference partition based on element count for this p-level
              ndiff = num_pelem_part(proc1,ilevel) - num_pelem_part(proc2,ilevel)
              if (ndiff > int( 0.01*num_pelem_part(proc2,ilevel) )) then
                proc_new = proc2
              else
                proc_new = proc1
              endif
            else
              ! sets new reference partition based on element count
              ndiff = part_bin(proc1) - part_bin(proc2)
              if (ndiff > int( 0.01*part_bin(proc2) )) then
                proc_new = proc2
              else
                proc_new = proc1
              endif
            endif

            ! sets new reference partition based on load
            !ndiff = part_load(proc_new) - part_load(proc2)
            !if (ndiff > int( 0.01*part_load(proc2) )) then
            !  proc_new = proc2
            !endif
          endif
        endif
      enddo

      ! checks if we need to re-assign
      if (found) cycle

      ! needs to re-assign partitions, adds elements to reference proc_new
      ! adds reference element
      part(e1) = proc_new
      if (proc1 /= proc_new) then
        ! element count change
        part_bin(proc1) = part_bin(proc1) - 1
        part_bin(proc_new) = part_bin(proc_new) + 1
        ! load change
        part_load(proc1) = part_load(proc1) - elmnts_load(e1)
        part_load(proc_new) = part_load(proc_new) + elmnts_load(e1)
        ! LTS
        if (LTS_MODE) then
          ilevel = p_lookup(ispec_p_refine(e1+1))
          num_pelem_part(proc1,ilevel) = num_pelem_part(proc1,ilevel) - 1
          num_pelem_part(proc_new,ilevel) = num_pelem_part(proc_new,ilevel) + 1
        endif
      endif
      ! add all other connected elements
      do iflt = 1,size(faults)
        e2 = fault_elements_connected(iflt,e1)
        if (e2 >= 0) then
          proc2 = part(e2)
          part(e2) = proc_new
          if (proc2 /= proc_new) then
            ! element count change
            part_bin(proc2) = part_bin(proc2) - 1
            part_bin(proc_new) = part_bin(proc_new) + 1
            ! load change
            part_load(proc2) = part_load(proc2) - elmnts_load(e2)
            part_load(proc_new) = part_load(proc_new) + elmnts_load(e2)
            ! LTS
            if (LTS_MODE) then
              ilevel = p_lookup(ispec_p_refine(e2+1))
              num_pelem_part(proc2,ilevel) = num_pelem_part(proc2,ilevel) - 1
              num_pelem_part(proc_new,ilevel) = num_pelem_part(proc_new,ilevel) + 1
            endif
          endif
        endif
      enddo
    enddo

  endif ! FAULT_BALANCE_GREEDY

  ! frees temporary arrays
  deallocate(fault_elements_connected)
  if (LTS_MODE) then
    deallocate(p_lookup)
    deallocate(num_pelem_part)
  endif

  ! makes sure, that elements which are only touching fault nodes (edge/corner) are also in the same partition
  if (ADD_FAULT_TOUCH_ELEMENTS) then
    ! user output
    print *,'  adding elements touching faults'

    ! allocates temporary arrays
    allocate(nnodes_elmnts(nnodes), &
             nodes_elmnts(nsize*nnodes),stat=ier)
    if (ier /= 0) stop 'Error allocating temporary mask_iglob'
    nnodes_elmnts(:) = 0
    nodes_elmnts(:) = 0

    ! builds list of elements per node
    do ispec = 0, nelmnts-1
      do inode = 1, NGNOD
        ! note: elmnts is defined from 0 to nspec-1 ( = nelmnts-1) and
        !       at this stage it stores shifted global indices starting from 0 to nnodes - 1
        ! node index (between 1 and nnodes)
        iglob = elmnts(inode,ispec) + 1
        if (iglob < 1 .or. iglob > nnodes) then
          print *,'Error: node index',iglob,'should be between 1 and',nnodes
          stop 'Error node index exceeds bounds in fault_repartition_parallel'
        endif
        ! sets element into node list
        nodes_elmnts( ((iglob-1)*nsize + 1) + nnodes_elmnts(iglob)) = ispec
        ! increases element counter
        nnodes_elmnts(iglob) = nnodes_elmnts(iglob) + 1
      enddo
    enddo

    ! masks out elements which have been assigned to faults already
    mask_ispec(:) = .false.
    do iflt = 1,size(faults)
      do e = 1, faults(iflt)%nspec
        ! note: index in part(:) starts from 0 to nelmnts (= nspec) by definition above, mask will do the same
        e1 = faults(iflt)%ispec1(e) - 1
        mask_ispec(e1) = .true.
        e2 = faults(iflt)%ispec2(e) - 1
        mask_ispec(e2) = .true.
      enddo
    enddo

    allocate(mask_iglob(nnodes),stat=ier)
    if (ier /= 0) stop 'Error allocating temporary mask_iglob'
    mask_iglob(:) = .false.

    ! masks nodes belonging to faults
    do iflt = 1,size(faults)
      do e = 1, faults(iflt)%nspec
        do k = 1,NGNOD2D
          ! nodes from fault surface 'up'
          iglob = faults(iflt)%inodes1(k,e)
          if (iglob < 1 .or. iglob > nnodes) stop 'Error iglob node index on fault surface up exceeds bounds'
          mask_iglob(iglob) = .true.
          ! nodes from fault surface 'down'
          iglob = faults(iflt)%inodes2(k,e)
          if (iglob < 1 .or. iglob > nnodes) stop 'Error iglob node index on fault surface down exceeds bounds'
          mask_iglob(iglob) = .true.
        enddo
      enddo
    enddo

    ! searches elements touching faults
    do iglob = 1,nnodes
      ! checks if node belongs to a fault
      if (mask_iglob(iglob) .eqv. .true.) then

        ! gets partition from (first encountered) fault element connected to this node
        found = .false.
        proc_new = -1
        do e = 1, nnodes_elmnts(iglob)
          ! note: element ids start from 0
          e1 = nodes_elmnts( (iglob-1)*nsize + e )
          if (mask_ispec(e1) .eqv. .true.) then
            proc_new = part(e1)
            found = .true.
            exit
          endif
        enddo

        ! checks if a fault element was found
        if (.not. found) stop 'Error fault node belongs to no fault elements'

        ! sets all other elements to same partition
        do e = 1, nnodes_elmnts(iglob)
          ! note: element ids start from 0
          e1 = nodes_elmnts( (iglob-1)*nsize + e )
          if (mask_ispec(e1) .eqv. .false.) then
            ! original partition number assigned to this element
            proc1 = part(e1)

            ! search if element shares an edge on the fault
            do inode = 1,NGNOD
              ! node index (between 1 and nnodes)
              iglob1 = elmnts(inode,e1) + 1

              ! skips current reference node
              if (iglob1 == iglob) cycle

              ! checks index
              if (iglob1 < 1 .or. iglob1 > nnodes) then
                print *,'Error: node index',iglob,'should be between 1 and',nnodes
                stop 'Error node index of touching element exceeds bounds in fault_repartition_parallel'
              endif

              ! checks if this node belongs also to a fault
              if (mask_iglob(iglob1) .eqv. .true.) then
                ! search for common element in node iglob and iglob1
                found = .false.
                do k1 = 1, nnodes_elmnts(iglob1)
                  e2 = nodes_elmnts( (iglob1 - 1)*nsize + k1 )
                  if (mask_ispec(e2)) then
                    ! element e2 is a fault element
                    ! looks for e2 also in reference node list
                    do k = 1,nnodes_elmnts(iglob)
                      e3 = nodes_elmnts( (iglob - 1)*nsize + k )
                      if (e2 == e3) then
                        found = .true.
                        proc_new = part(e2)
                        exit
                      endif
                    enddo
                  endif
                  ! checks if we need to continue search
                  if (found) exit
                enddo
              endif
            enddo

            ! sets new partition
            part(e1) = proc_new

            ! updates accounts for load/element change in partition
            if (proc1 /= proc_new) then
              ! element count change
              part_bin(proc1) = part_bin(proc1) - 1
              part_bin(proc_new) = part_bin(proc_new) + 1
              ! load change
              part_load(proc1) = part_load(proc1) - elmnts_load(e1)
              part_load(proc_new) = part_load(proc_new) + elmnts_load(e1)
            endif
          endif
        enddo
      endif
    enddo

    ! frees temporary arrays
    deallocate(nnodes_elmnts,nodes_elmnts)
    deallocate(mask_iglob)
  endif ! ADD_FAULT_TOUCH_ELEMENTS

  ! user output
  do i = 0,maxval(part(:))
    print *,'  partition :',i,'initial elements =',part_bin0(i),' - fault elements added =',part_bin(i)
    !print *,'                            initial load     =',part_load0(i),' - load change =',part_load(i) - part_load0(i)
  enddo
  print *

  ! frees memory
  deallocate(part_bin,part_bin0)
  deallocate(mask_ispec)

  end subroutine fault_repartition_parallel

! ---------------------------------------------------------------------------------------------------

  subroutine setup_connected_elements(nelmnts,mask_ispec,fault_elements_connected)

  implicit none

  integer, intent(in) :: nelmnts
  logical,dimension(0:nelmnts-1),intent(inout) :: mask_ispec
  integer,dimension(size(faults),0:nelmnts-1) :: fault_elements_connected

  ! local parameters
  integer :: i,j,iflt,e,e1,e2
  logical :: found

  ! initializes array
  fault_elements_connected(:,:) = -1

  ! checks if anything to do
  if (size(faults) < 1) return

  ! setup connected elements by faults, to handle triple junctions (intersections between two faults)
  mask_ispec(:) = .false.
  do e = 1,faults(1)%nspec
    e1 = faults(1)%ispec1(e) - 1
    ! puts associated e2 element to connected list
    e2 = faults(1)%ispec2(e) - 1
    fault_elements_connected(1,e1) = e2
    mask_ispec(e1) = .true.
  enddo

  ! adds connected elements from other faults
  if (size(faults) > 1) then
    do iflt = 2,size(faults)
      do e = 1,faults(iflt)%nspec
        e1 = faults(iflt)%ispec1(e) - 1
        e2 = faults(iflt)%ispec2(e) - 1
        ! searches if element e1 has been added as reference element to connected list
        if (mask_ispec(e1)) then
          ! e1 is already a reference element, adds e2 as connected
          fault_elements_connected(iflt,e1) = e2
        else if (mask_ispec(e2)) then
          ! e2 is already a reference element, adds e1 as connected
          fault_elements_connected(iflt,e2) = e1
        else
          ! search elements in previous faults
          found = .false.
          do i = 0,nelmnts-1
            do j = 1,iflt-1
              if (fault_elements_connected(j,i) == e1) then
                ! adds e2 as newly associated element
                fault_elements_connected(iflt,i) = e2
                found = .true.
                exit
              else if (fault_elements_connected(j,i) == e2) then
                ! adds e1 as newly associated element
                fault_elements_connected(iflt,i) = e2
                found = .true.
                exit
              endif
            enddo
            if (found) exit
          enddo
          if (.not. found) then
            ! fault elements not connected to other faults, adds new pair
            fault_elements_connected(iflt,e1) = e2
            mask_ispec(e1) = .true.
          endif
        endif
      enddo
    enddo
  endif

  end subroutine setup_connected_elements

! ---------------------------------------------------------------------------------------------------
! See subroutine write_boundaries_database in part_decompose_mesh_SCOTCH.f90
!
! File format:
! one block for each fault
! first line of each block = number of fault elements in this processor
! next lines: #id_(element containing the face) #id_node1_face .. #id_node4_face
! first for all faces on side 1, then side 2

  subroutine write_fault_database(IIN_database, iproc, nelmnts, glob2loc_elmnts, &
                                  glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, part)

  implicit none
  integer, intent(in) :: IIN_database
  integer, intent(in) :: iproc
  integer, intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(in)  :: part
  integer, dimension(0:nelmnts-1), intent(in)  :: glob2loc_elmnts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_nparts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_parts
  integer, dimension(0:), intent(in) :: glob2loc_nodes

  integer :: i,j,k,iflt,e
  integer :: nspec_fault_1,nspec_fault_2
  integer :: loc_nodes(NGNOD2D),inodes(NGNOD2D)

  do iflt = 1,size(faults)
    ! get number of fault elements in this partition
    nspec_fault_1 = count( part(faults(iflt)%ispec1(:)-1) == iproc )
    nspec_fault_2 = count( part(faults(iflt)%ispec2(:)-1) == iproc )

    if (nspec_fault_1 /= nspec_fault_2) then
      print *, 'Fault # ',iflt,', proc # ',iproc
      print *, '  ispec1 : ', nspec_fault_1
      print *, '  ispec2 : ', nspec_fault_2
      print *, 'Fatal error: Number of fault elements do not coincide. Abort.'
      stop
    endif
    write(IIN_database) nspec_fault_1

    ! if no fault element in this partition, move to next fault
    if (nspec_fault_1 == 0) cycle

    ! export fault element data, side 1
    do i = 1,faults(iflt)%nspec
      e = faults(iflt)%ispec1(i)
      if (part(e-1) == iproc) then
        inodes = faults(iflt)%inodes1(:,i)
        do k = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(inodes(k)-1), glob2loc_nodes_nparts(inodes(k))-1
            if (glob2loc_nodes_parts(j) == iproc) then
              loc_nodes(k) = glob2loc_nodes(j) + 1
            endif
          enddo
        enddo
        write(IIN_database) glob2loc_elmnts(e-1)+1, loc_nodes
      endif
    enddo

    ! export fault element data, side 2
    do i = 1,faults(iflt)%nspec
      e = faults(iflt)%ispec2(i)
      if (part(e-1) == iproc) then
        inodes = faults(iflt)%inodes2(:,i)
        do k = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(inodes(k)-1), glob2loc_nodes_nparts(inodes(k))-1
            if (glob2loc_nodes_parts(j) == iproc) then
              loc_nodes(k) = glob2loc_nodes(j)+1
            endif
          enddo
        enddo
        write(IIN_database) glob2loc_elmnts(e-1)+1, loc_nodes
      endif
    enddo
  enddo

  end subroutine write_fault_database


! ---------------------------------------------------------------------------------------------------
! Below functions are only used for the MPI version
! ---------------------------------------------------------------------------------------------------

  subroutine write_fault_database_mpi(IIN_database, myrank, nE, glob2loc_elmnt, ipart, &
                                      nnodes, glob2loc_nodes, NGNOD2D)

  implicit none
  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: myrank
  integer, intent(in)  :: nE, nnodes, NGNOD2D
  integer, dimension(nE), intent(in)  :: ipart
  integer, dimension(nE), intent(in)  :: glob2loc_elmnt
  integer, dimension(nnodes), intent(in) :: glob2loc_nodes

  !integer, intent(in)  :: nnodes_loc
  !double precision, dimension(NDIM,nnodes), intent(in)  :: nodes_coords

  integer :: iflt, i, inode
  integer :: iE, iE_loc
  integer :: nspec_fault_1,nspec_fault_2
  integer :: inodes(NGNOD2D), node_loc(NGNOD2D)

  do iflt=1,size(faults)
    ! get number of fault elements in this partition
    nspec_fault_1 = count( ipart(faults(iflt)%ispec1) == myrank + 1)
    nspec_fault_2 = count( ipart(faults(iflt)%ispec2) == myrank + 1)

    if (nspec_fault_1 /= nspec_fault_2) then
      print *, 'Fault # ',iflt,', proc # ',myrank
      print *, '  ispec1 : ', nspec_fault_1
      print *, '  ispec2 : ', nspec_fault_2
      print *, 'Fatal error: Number of fault elements do not coincide. Abort.'
      stop
    endif
    write(IIN_database) nspec_fault_1

    ! if no fault element in this partition, move to next fault
    if (nspec_fault_1 == 0) cycle

    do i=1, faults(iflt)%nspec
      iE = faults(iflt)%ispec1(i)
      if (ipart(iE) /= myrank +1) cycle
      iE_loc=glob2loc_elmnt(iE)

      inodes = faults(iflt)%inodes1(:,i)
      do inode = 1, NGNOD2D
        node_loc(inode) = glob2loc_nodes(inodes(inode))
      enddo
      write(IIN_database) iE_loc, node_loc(1:NGNOD2D)
    enddo

    do i=1, faults(iflt)%nspec
      iE = faults(iflt)%ispec2(i)
      if (ipart(iE) /= myrank +1) cycle
      iE_loc=glob2loc_elmnt(iE)

      inodes = faults(iflt)%inodes2(:,i)
      do inode = 1, NGNOD2D
        node_loc(inode) = glob2loc_nodes(inodes(inode))
      enddo
      write(IIN_database) iE_loc, node_loc(1:NGNOD2D)
    enddo
  enddo

  end subroutine write_fault_database_mpi

! ---------------------------------------------------------------------------------------------------

  subroutine bcast_faults_mpi(myrank)

  implicit none
  integer, intent(in)  :: myrank
  integer :: nbfaults  , iflt, ier, nspec

  if (myrank == 0) then
    nbfaults = size(faults)
  endif
  call bcast_all_singlei(nbfaults)

  if (myrank /= 0) then
    allocate(faults(nbfaults),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
  endif

  do iflt = 1 , nbfaults
    call bcast_all_singlei(faults(iflt)%nspec)

    nspec = faults(iflt)%nspec
    if (myrank /= 0) then
      allocate(faults(iflt)%ispec1(nspec),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
      allocate(faults(iflt)%ispec2(nspec),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
      allocate(faults(iflt)%inodes1(NGNOD2D,nspec),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
      allocate(faults(iflt)%inodes2(NGNOD2D,nspec),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('Error allocating array 78')
    endif
    call bcast_all_i(faults(iflt)%ispec1,  nspec)
    call bcast_all_i(faults(iflt)%ispec2,  nspec)
    call bcast_all_i(faults(iflt)%inodes1, nspec*NGNOD2D)
    call bcast_all_i(faults(iflt)%inodes2, nspec*NGNOD2D)
  enddo

  end subroutine bcast_faults_mpi

end module fault_scotch

