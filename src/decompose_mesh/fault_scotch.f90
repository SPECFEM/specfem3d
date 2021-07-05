!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  use constants, only: MAX_STRING_LEN, IN_DATA_FILES, NDIM
  use shared_parameters, only: NGNOD2D

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

  logical, save :: ANY_FAULT = .false.

  integer, parameter :: long = SELECTED_INT_KIND(18)

  double precision, parameter :: FAULT_GAP_TOLERANCE = 1.0d0
                              ! must be larger than the fault offset in the mesh,
                              ! but smaller than the smallest element size

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

  open(unit=IIN_PAR,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file_faults',status='old',action='read',iostat=ier)
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

  allocate(faults(nbfaults),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
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

  ! reads fault elements and nodes
 ! File format:
 ! Line 1:
 !   number_of_elements_in_side_1   number_of_elements_in_side_2
 ! Then for all elements that have a face on side 1:
 !   #id_element #id_global_node1 .. #id_global_node4
 ! Then the same for side 2.
 ! Note: element ids start at 1, not 0 (see cubit2specfem3d.py)
  open(unit=IIN_FLT, file=filename, status='old', form='formatted', iostat = ier)
  if (ier /= 0) then
    write(*,*) 'Fatal error: file '//filename//' not found'
    write(*,*) 'Abort'
    stop
  endif

  read(IIN_FLT,*) nspec_side1,nspec_side2
  if (nspec_side1 /= nspec_side2) stop 'Number of fault nodes at do not match.'
  f%nspec = nspec_side1
  allocate(f%ispec1(f%nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 79')
  allocate(f%ispec2(f%nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 80')
  allocate(f%inodes1(NGNOD2D,f%nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 81')
  allocate(f%inodes2(NGNOD2D,f%nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 82')

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

  end subroutine read_single_fault_file


! ---------------------------------------------------------------------------------------------------

! Saving nodes_coords to be used in the solver for ibool_fault_side1 and side2

  subroutine save_nodes_coords(nodes_coords,nnodes)

  implicit none
  integer, intent(in) :: nnodes
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords
  integer :: ier

  allocate(nodes_coords_open(NDIM,nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 83')
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

  subroutine close_fault_single(f,nodes_coords,nnodes)

  implicit none
  type(fault_type), intent(in) :: f
  integer, intent(in)  :: nnodes
  double precision, dimension(NDIM,nnodes), intent(inout) :: nodes_coords

  double precision, dimension(NDIM) :: xyz_1, xyz_2, xyz

  double precision :: dist
  integer :: iglob1, iglob2, i, j, k1, k2
  logical :: found_it

  do i = 1,f%nspec
    do k2 = 1,NGNOD2D
      iglob2 = f%inodes2(k2,i)
      found_it = .false.
      xyz_2(:) = nodes_coords(:,iglob2)

      do j = 1,f%nspec
        do k1=1,NGNOD2D
          iglob1 = f%inodes1(k1,j)
          xyz_1(:) = nodes_coords(:,iglob1)
          xyz(:) = xyz_2(:) - xyz_1(:)
          dist = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))

          !jpa: Closing nodes that are already closed is not a problem
          !jpa: I process them again to leave the loop as early as possible
          !jpa: and to test if a facing node was found (see below).

          if (dist <= FAULT_GAP_TOLERANCE) then
            xyz(:) =  (xyz_1(:) + xyz_2(:))*0.5d0
            nodes_coords(:,iglob2) = xyz(:)
            nodes_coords(:,iglob1) = xyz(:)
            found_it = .true.
            exit
          endif

        enddo
        if (found_it) exit
      enddo

      ! jpa: If the two fault sides have been meshed independently they might not match. Test it here:
      if (.not. found_it) stop 'Inconsistent fault mesh: corresponding node in the other fault face was not found'

    enddo
  enddo

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
  xyz_c(:,:) = xyz_c(:,:) / 4.d0

  ! reorder
  call lex_order(xyz_c,new_index_list,f%nspec)
  f%ispec1 = f%ispec1(new_index_list)
  f%inodes1 = f%inodes1(:,new_index_list)

  ! repeat for fault side 2
  xyz_c = 0d0
  do e = 1,f%nspec
    do k = 1,4
      iglob = f%inodes2(k,e)
      xyz_c(:,e) = xyz_c(:,e) + nodes_coords(:,iglob)
    enddo
  enddo
  xyz_c(:,:) = xyz_c(:,:) / 4.d0

  ! reorder
  call lex_order(xyz_c,new_index_list,f%nspec)
  f%ispec2 = f%ispec2(new_index_list)
  f%inodes2 = f%inodes2(:,new_index_list)

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

  xp = xyz_c(1,:)
  yp = xyz_c(2,:)
  zp = xyz_c(3,:)

  ! define geometrical tolerance based upon typical size of the model
  xtol = 1.d-10 * maxval( maxval(xyz_c,2) - minval(xyz_c,2) )

  call sort_array_coordinates(nspec,xp,yp,zp,ibool,iglob,locval,ifseg, &
                              nglob,ninseg,xtol)

  end subroutine lex_order

!===================================================================================================
  !--------------------------------------------------
  ! Repartitioning : two coupled elements on fault side1/side2 are transfered to the same partition
  !--------------------------------------------------

  subroutine fault_repartition(nelmnts, nnodes, elmnts, nsize, nproc, part, esize, nodes_coords)

  use constants, only: PARALLEL_FAULT

  implicit none
  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, esize
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords

  if (PARALLEL_FAULT) then
    ! fault can be parallelized
    call fault_repartition_parallel(nelmnts,part, nodes_coords, nnodes, nproc)
  else
    ! move all fault elements to the same partition (proc=0)
    call fault_repartition_not_parallel(nelmnts, nnodes, elmnts, nsize, nproc, part, esize)
  endif

  end subroutine fault_repartition

!---------------------------------------------------------------------------------------------------

!     part(e) = index of the processor to which element #e is assigned.
!               Fault elements and neighbors are assigned to the same processor.
!               Part, once modified, will be input for write_partition_database.

  subroutine fault_repartition_not_parallel(nelmnts, nnodes, elmnts, nsize, nproc, part, esize)

  implicit none
  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, esize
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part

  integer, dimension(0:nnodes-1) :: nnodes_elmnts
  integer, dimension(0:nsize*nnodes-1) :: nodes_elmnts
  integer  :: i,ipart,nproc_null,nproc_null_final
  integer  :: k1, k2, k,e,iflt,inode,ier
  integer, dimension(:), allocatable :: elem_proc_null

  ! user output
  print *,'fault repartitioning: (non-parallel version)'

 ! downloading processor 0
  nproc_null = count( part(:) == 0)
  print *,'  elements in slice proc 0 redistributed in [{nproc}- nproc0] :',nproc_null

  ! distributs elements in slice 0 to higher procs
  if (nproc > 1) then
    if (nproc_null /= 0) then
      allocate(elem_proc_null(nproc_null),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 84')
      ! Filling up proc = 0 elements
      nproc_null = 0
      do i = 0,nelmnts-1
        if (part(i) == 0) then
          nproc_null = nproc_null + 1
          elem_proc_null(nproc_null) = i
        endif
      enddo
      ! Redistributing proc-0 elements on the rest of processors
      ipart = 1
      if (nproc > 1) then
        do i = 1, nproc_null
          part(elem_proc_null(i)) = ipart
          if (ipart == nproc-1) ipart = 0
          ipart = ipart +1
        enddo
      endif
      deallocate(elem_proc_null)
    endif
  endif

! Fault zone layer = the set of elements that contain at least one fault node
  print *,"  fault zone layer :"

! List of elements per node
!  nnodes_elmnts(i) = number of elements containing node #i (i=0:nnodes-1)
!  nodes_elmnts(nsize*i:nsize*i+nnodes_elmnts(i)-1) = index of elements (starting at 0) containing node #i
!  nsize = maximun number of elements in a node.
!  esize = nodes per element.

  ! repartitions elements
  if (nproc > 1) then
    nnodes_elmnts(:) = 0
    nodes_elmnts(:)  = 0

    do i = 0, esize*nelmnts-1
      nodes_elmnts(elmnts(i)*nsize + nnodes_elmnts(elmnts(i))) = i/esize
      nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
    enddo

    do iflt = 1,size(faults)
      do e = 1,faults(iflt)%nspec
        do k = 1,4
          ! 1. fault side
          inode = faults(iflt)%inodes1(k,e)-1  ! node index, starting at 0
          k1 = nsize*inode
          k2 = k1 + nnodes_elmnts(inode) -1
          ! moves elements to proc0
          part( nodes_elmnts(k1:k2) ) = 0
          ! 2. fault side
          inode = faults(iflt)%inodes2(k,e)-1
          k1 = nsize*inode
          k2 = k1 + nnodes_elmnts(inode) -1
          ! moves elements to proc0
          part( nodes_elmnts(k1:k2) ) = 0
        enddo
      enddo
    enddo
  endif

  nproc_null_final = count( part(:) == 0)
  print *,'  final number of elements in slice proc 0 (containing fault surfaces):',nproc_null
  print *

  end subroutine fault_repartition_not_parallel

! ---------------------------------------------------------------------------------------------------

  subroutine fault_repartition_parallel(nelmnts, part, nodes_coords, nnodes, nproc)

  implicit none
  integer, intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(inout) :: part
  integer, intent(in) :: nnodes,nproc
  double precision, dimension(NDIM,nnodes), intent(in) :: nodes_coords

  integer  :: i,iflt,e,e1,e2,proc1,proc2,counter,num_fault_elem

  ! user output
  print *,'fault repartitioning: (parallel version)'

 ! Reorder both fault sides so that elements facing each other have the same index
  call reorder_fault_elements(nodes_coords,nnodes)

  num_fault_elem = 0
  do iflt = 1,size(faults)
    num_fault_elem = num_fault_elem + faults(iflt)%nspec
  enddo
  print *,'  total number of fault elements: ',num_fault_elem
  print *

  ! repartitions elements
  if (nproc > 1) then
!JPA loop over all faults
!JPA loop over all fault element pairs
!JPA assign both elements to the processor with lowest rank among the pair

    do i = 1,2 !JPA do it twice, to handle triple junctions (intersections between two faults)
      counter = 0
      do iflt = 1,size(faults)
        do e = 1,faults(iflt)%nspec
          e1 = faults(iflt)%ispec1(e) - 1  ! part array between [0,nelmnts-1]
          e2 = faults(iflt)%ispec2(e) - 1
          proc1 = part(e1)
          proc2 = part(e2)
          ! moves elements to lower partition of the two
          part(e1) = min(proc1,proc2)
          part(e2) = part(e1)
          counter = counter + 1
        enddo
      enddo
      print *,'  repartition loop',i,' changed ',counter,' coupled fault elements out of ',num_fault_elem
    enddo
    print *
  endif

  end subroutine fault_repartition_parallel

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

  integer  :: i,j,k,iflt,e
  integer  :: nspec_fault_1,nspec_fault_2
  integer :: loc_nodes(NGNOD2D),inodes(NGNOD2D)

  do iflt = 1,size(faults)

   ! get number of fault elements in this partition
    nspec_fault_1 = count( part(faults(iflt)%ispec1-1) == iproc )
    nspec_fault_2 = count( part(faults(iflt)%ispec2-1) == iproc )

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
        do k=1,NGNOD2D
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
        do k=1,NGNOD2D
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
             nnodes, glob2loc_nodes, NGNOD2D, nnodes_loc, nodes_coords)
  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: myrank
  integer, intent(in)  :: nE, nnodes, NGNOD2D
  integer, dimension(nE), intent(in)  :: ipart
  integer, dimension(nE), intent(in)  :: glob2loc_elmnt
  integer, dimension(nnodes), intent(in) :: glob2loc_nodes

  integer, intent(in)  :: nnodes_loc
  double precision, dimension(NDIM,nnodes), intent(in)  :: nodes_coords

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


  subroutine bcast_faults_mpi(myrank)

  integer, intent(in)  :: myrank
  integer :: nbfaults  , iflt, ier, nspec

  if (myrank == 0) then
      nbfaults = size(faults)
  endif
  call bcast_all_singlei(nbfaults)
  if (myrank /= 0 ) then
      allocate(faults(nbfaults),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
  endif

  do iflt = 1 , nbfaults
     call bcast_all_singlei(faults(iflt)%nspec)
     nspec = faults(iflt)%nspec
     if (myrank /= 0 ) then
        allocate(faults(iflt)%ispec1(nspec),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
        allocate(faults(iflt)%ispec2(nspec),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
        allocate(faults(iflt)%inodes1(NGNOD2D,nspec),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
        allocate(faults(iflt)%inodes2(NGNOD2D,nspec),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating array 78')
     endif
     call bcast_all_i(faults(iflt)%ispec1,  nspec)
     call bcast_all_i(faults(iflt)%ispec2,  nspec)
     call bcast_all_i(faults(iflt)%inodes1, nspec*NGNOD2D)
     call bcast_all_i(faults(iflt)%inodes2, nspec*NGNOD2D)
  enddo

  end subroutine bcast_faults_mpi

end module fault_scotch

