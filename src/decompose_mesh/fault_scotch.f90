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

module fault_scotch

  use constants, only: MAX_STRING_LEN, IN_DATA_FILES
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

  logical, parameter :: PARALLEL_FAULT = .true.
 ! NOTE: PARALLEL_FAULT has to be the same
 !       in fault_solver_common.f90, fault_generate_databases.f90 and fault_scotch.f90

  integer, parameter :: long = SELECTED_INT_KIND(18)

  double precision, parameter :: FAULT_GAP_TOLERANCE = 1.0d0
                              ! must be larger than the fault offset in the mesh,
                              ! but smaller than the smallest element size

  public :: read_fault_files, fault_repartition, close_faults, write_fault_database, &
            save_nodes_coords, nodes_coords_open, ANY_FAULT

CONTAINS
!==========================================================================================

  subroutine read_fault_files(localpath_name)

  character(len=MAX_STRING_LEN), intent(in) :: localpath_name
  integer :: nbfaults, iflt, ier

  open(101,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file_faults',status='old',action='read',iostat=ier)
  if (ier == 0) then
    read(101,*) nbfaults
  else
    nbfaults = 0
    print *, 'Par_file_faults not found: assume no faults'
  endif
  close(101)

  ANY_FAULT = (nbfaults > 0)
  if (.not. ANY_FAULT)  return

  allocate(faults(nbfaults))
 ! NOTE: asumes that the fault ids follow a contiguous numbering, starting at 1, with unit increment
 !       The user must assign that numbering during mesh generation
  do iflt = 1 , nbfaults
   call read_single_fault_file(faults(iflt),iflt,localpath_name)
  enddo

  end subroutine read_fault_files


!---------------------------------------------------------------------------------------------------

  subroutine read_single_fault_file(f,ifault,localpath_name)

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
  open(unit=101, file=filename, status='old', form='formatted', iostat = ier)
  if (ier /= 0) then
    write(6,*) 'Fatal error: file '//filename//' not found'
    write(6,*) 'Abort'
    stop
  endif

  read(101,*) nspec_side1,nspec_side2
  if (nspec_side1 /= nspec_side2) stop 'Number of fault nodes at do not match.'
  f%nspec = nspec_side1
  allocate(f%ispec1(f%nspec))
  allocate(f%ispec2(f%nspec))
  allocate(f%inodes1(NGNOD2D,f%nspec))
  allocate(f%inodes2(NGNOD2D,f%nspec))
  do e=1,f%nspec
    read(101,*) f%ispec1(e),f%inodes1(:,e)
  enddo
  do e=1,f%nspec
    read(101,*) f%ispec2(e),f%inodes2(:,e)
  enddo
 ! If we ever figure out how to export "ifaces" from CUBIT:
  !allocate(f%iface1(f%nspec))
  !allocate(f%iface2(f%nspec))
  !do e=1,f%nspec
  !  read(101,*) f%ispec1(e),f%ispec2(e),f%iface1(e),f%iface2(e)
  !enddo

  close(101)

  end subroutine read_single_fault_file


! ---------------------------------------------------------------------------------------------------
! Saving nodes_coords to be used in the solver for ibool_fault_side1 and side2
   subroutine save_nodes_coords(nodes_coords,nnodes)

   integer, intent(in) :: nnodes
   double precision, dimension(3,nnodes), intent(in) :: nodes_coords

   allocate(nodes_coords_open(3,nnodes))
   nodes_coords_open = nodes_coords

   end subroutine save_nodes_coords

! ---------------------------------------------------------------------------------------------------

  subroutine close_faults(nodes_coords,nnodes)

  integer, intent(in) :: nnodes
  double precision, dimension(3,nnodes), intent(inout) :: nodes_coords

  integer  :: i

  do i=1,size(faults)
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

  type(fault_type), intent(in) :: f
  integer, intent(in)  :: nnodes
  double precision, dimension(3,nnodes), intent(inout) :: nodes_coords

  double precision, dimension(3) :: xyz_1, xyz_2, xyz

  double precision :: dist
  integer :: iglob1, iglob2, i, j, k1, k2
  logical :: found_it

  do i = 1,f%nspec
    do k2=1,NGNOD2D
      iglob2 = f%inodes2(k2,i)
      found_it = .false.
      xyz_2 = nodes_coords(:,iglob2)

      do j = 1,f%nspec
        do k1=1,NGNOD2D
          iglob1 = f%inodes1(k1,j)
          xyz_1 = nodes_coords(:,iglob1)
          xyz = xyz_2-xyz_1
          dist = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))

         !jpa: Closing nodes that are already closed is not a problem
         !jpa: I process them again to leave the loop as early as possible
         !jpa: and to test if a facing node was found (see below).

          if (dist <= FAULT_GAP_TOLERANCE) then
            xyz =  (xyz_1 + xyz_2)*0.5d0
            nodes_coords(:,iglob2) = xyz
            nodes_coords(:,iglob1) = xyz
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

  integer, intent(in)  :: nnodes
  double precision,dimension(3,nnodes), intent(in) :: nodes_coords

  integer :: i

  do i=1,size(faults)
    call reorder_fault_elements_single(faults(i),nodes_coords,nnodes)
  enddo

  end subroutine reorder_fault_elements

! ---------------------------------------------------------------------------------------------------
  subroutine reorder_fault_elements_single(f,nodes_coords,nnodes)

  type(fault_type), intent(inout) :: f
  integer, intent(in) :: nnodes
  double precision, dimension(3,nnodes), intent(in) :: nodes_coords

  double precision, dimension(3,f%nspec) :: xyz_c
  integer :: k,e,iglob
  integer, dimension(f%nspec) :: new_index_list

  ! compute element-face centers for fault side 1
  xyz_c = 0d0
  do e=1,f%nspec
    do k=1,4
      iglob = f%inodes1(k,e)
      xyz_c(:,e) = xyz_c(:,e) + nodes_coords(:,iglob)
    enddo
  enddo
  xyz_c = xyz_c / 4d0
  ! reorder
  call lex_order(xyz_c,new_index_list,f%nspec)
  f%ispec1 = f%ispec1(new_index_list)
  f%inodes1 = f%inodes1(:,new_index_list)

  ! repeat for fault side 2
  xyz_c = 0d0
  do e=1,f%nspec
    do k=1,4
      iglob = f%inodes2(k,e)
      xyz_c(:,e) = xyz_c(:,e) + nodes_coords(:,iglob)
    enddo
  enddo
  xyz_c = xyz_c / 4d0
  call lex_order(xyz_c,new_index_list,f%nspec)
  f%ispec2 = f%ispec2(new_index_list)
  f%inodes2 = f%inodes2(:,new_index_list)

  end subroutine reorder_fault_elements_single

! ---------------------------------------------------------------------------------------------------
  subroutine lex_order(xyz_c,locval,nspec)

  integer, intent(in) :: nspec
  integer, intent(out) :: locval(nspec)
  double precision, intent(in) :: xyz_c(3,nspec)

  double precision, dimension(nspec) :: xp,yp,zp
  integer, dimension(nspec) :: ninseg,ibool,iglob
  integer :: nglob
  logical :: ifseg(nspec)
  double precision :: xtol

  xp=xyz_c(1,:)
  yp=xyz_c(2,:)
  zp=xyz_c(3,:)

  ! define geometrical tolerance based upon typical size of the model
  xtol = 1.d-10 * maxval( maxval(xyz_c,2) - minval(xyz_c,2) )

  call sort_array_coordinates(nspec,xp,yp,zp,ibool,iglob,locval,ifseg, &
                              nglob,ninseg,xtol)

  end subroutine lex_order

!===================================================================================================
  !--------------------------------------------------
  ! Repartitioning : two coupled elements on fault side1/side2 are transfered to the same partition
  !--------------------------------------------------

  subroutine fault_repartition (nelmnts, nnodes, elmnts, nsize, nproc, part, esize, nodes_coords)

  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, esize
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part
  double precision, dimension(3,nnodes), intent(in) :: nodes_coords

  if (PARALLEL_FAULT) then
    call fault_repartition_parallel (nelmnts,part,nodes_coords,nnodes)
  else
    ! move all fault elements to the same partition (proc=0)
    call fault_repartition_not_parallel (nelmnts, nnodes, elmnts, nsize, nproc, part, esize)
  endif

  end subroutine fault_repartition

!---------------------------------------------------------------------------------------------------

!     part(e) = index of the processor to which element #e is assigned.
!               Fault elements and neighbors are assigned to the same processor.
!               Part, once modified, will be input for write_partition_database.

  subroutine fault_repartition_not_parallel (nelmnts, nnodes, elmnts, nsize, nproc, part, esize)

  integer, intent(in) :: nelmnts,nsize
  integer, intent(in) :: nnodes, nproc, esize
  integer, dimension(0:esize*nelmnts-1), intent(in) :: elmnts
  integer, dimension(0:nelmnts-1), intent(inout)    :: part

  integer, dimension(0:nnodes-1) :: nnodes_elmnts
  integer, dimension(0:nsize*nnodes-1) :: nodes_elmnts
  integer  :: i,ipart,nproc_null,nproc_null_final
  integer  :: k1, k2, k,e,iflt,inode
  integer, dimension(:), allocatable :: elem_proc_null

 ! downloading processor 0
  nproc_null = count( part == 0)

  print *, 'Elements proc = 0 redistributed in [{nproc}- nproc0] :'
  print *, nproc_null

  if (nproc_null /= 0) then

    allocate(elem_proc_null(nproc_null))
   ! Filling up proc = 0 elements
    nproc_null = 0
    do i = 0,nelmnts-1
      if (part(i) == 0) then
        nproc_null = nproc_null + 1
        elem_proc_null(nproc_null) = i
      endif
    enddo
   ! Redistributing proc-0 elements on the rest of processors
    ipart=1
    if (nproc > 1) then
      do i = 1, nproc_null
        part(elem_proc_null(i)) = ipart
        if (ipart == nproc-1) ipart = 0
        ipart = ipart +1
      enddo
    endif
    deallocate(elem_proc_null)
  endif

! Fault zone layer = the set of elements that contain at least one fault node
  print *, "Fault zone layer :"

! List of elements per node
!  nnodes_elmnts(i) = number of elements containing node #i (i=0:nnodes-1)
!  nodes_elmnts(nsize*i:nsize*i+nnodes_elmnts(i)-1) = index of elements (starting at 0) containing node #i
!  nsize = maximun number of elements in a node.
!  esize = nodes per element.

  nnodes_elmnts(:) = 0
  nodes_elmnts(:)  = 0

  do i = 0, esize*nelmnts-1
    nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
    nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
  enddo

  do iflt=1,size(faults)
    do e=1,faults(iflt)%nspec
      do k=1,4

        inode = faults(iflt)%inodes1(k,e)-1  ! node index, starting at 0
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part( nodes_elmnts(k1:k2) ) = 0
        inode = faults(iflt)%inodes2(k,e)-1
        k1 = nsize*inode
        k2 = k1 + nnodes_elmnts(inode) -1
        part( nodes_elmnts(k1:k2) ) = 0

      enddo
    enddo
  enddo

  nproc_null_final = count( part == 0)
  print *, nproc_null_final

  end subroutine fault_repartition_not_parallel

! ---------------------------------------------------------------------------------------------------
  subroutine fault_repartition_parallel (nelmnts, part, nodes_coords,nnodes)

  integer, intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(inout) :: part
  integer, intent(in) :: nnodes
  double precision, dimension(3,nnodes), intent(in) :: nodes_coords

  integer  :: i,iflt,e,e1,e2,proc1,proc2

 ! Reorder both fault sides so that elements facing each other have the same index
  call reorder_fault_elements(nodes_coords,nnodes)

!JPA loop over all faults
!JPA loop over all fault element pairs
!JPA assign both elements to the processor with lowest rank among the pair

  do i=1,2 !JPA do it twice, to handle triple junctions (intersections between two faults)
    do iflt=1,size(faults)
      do e=1,faults(iflt)%nspec
        e1 = faults(iflt)%ispec1(e) - 1
        e2 = faults(iflt)%ispec2(e) - 1
        proc1 = part(e1)
        proc2 = part(e2)
        part(e1) = min(proc1,proc2)
        part(e2) = part(e1)
      enddo
    enddo
  enddo

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

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer, intent(in) :: nelmnts
  integer, dimension(0:nelmnts-1), intent(in)  :: part
  integer, dimension(0:nelmnts-1), intent(in)  :: glob2loc_elmnts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_nparts
  integer, dimension(0:), intent(in) :: glob2loc_nodes_parts
  integer, dimension(0:), intent(in) :: glob2loc_nodes

  integer  :: i,j,k,iflt,e
  integer  :: nspec_fault_1,nspec_fault_2
  integer :: loc_nodes(NGNOD2D),inodes(NGNOD2D)

  do iflt=1,size(faults)

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
    do i=1,faults(iflt)%nspec
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
    do i=1,faults(iflt)%nspec
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

end module fault_scotch

