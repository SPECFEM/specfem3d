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

! Parallel heuristic mesh decomposer
!
! The decomposition is performed in the three Cartesian directions, given a
! desired number of partitions in each direction. The domain should be a box.
!
! Vadim Monteiller, CNRS Marseille, France, February 2016

program xdecompose_mesh_mpi

  use module_mesh
  use module_database
  use module_partition

  implicit none

  ! proc numbers for MPI
  integer             :: myrank,sizeprocs
  logical, parameter  :: BROADCAST_AFTER_READ=.true.

  ! number of proc in each direction
  integer             :: npartX, npartY, npartZ

  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! 1/ read Par_file
  call read_parameter_file(myrank, BROADCAST_AFTER_READ)

  ! 2/ read mesh
  if (myrank == 0) then

     print *,'         Parallel heuristic mesh decomposer'
     print *
     print *,'The decomposition is performed in the three Cartesian directions, given a'
     print *,'desired number of partitions in each direction. The domain should be a box.'
     print *
     print *,'Works for a topologically-rectangular domain in which the user has to'
     print *,'specify the number of slices to decompose the mesh into in the X, Y and Z directions.'
     print *,'Load imbalance is handled by taking into account the cost of the calculations inside each type'
     print *,'of spectral elements: acoustic, elastic, viscoelastic, ...'
     print *
     print *
     print *,'Usage: mpirun -np nparts xdecompose_mesh_mpi npartX npartY npartZ'
     print *
     print *,'  where'
     print *,'      npartX = number of mesh partitions to be created in the X direction'
     print *,'      npartY = number of mesh partitions to be created in the Y direction'
     print *,'      npartZ = number of mesh partitions to be created in the Z direction'
     print *
     print *,'  and nparts **MUST BE** equal to npartX * npartY * npartZ'
     print *
     print *
     print *,'      input_directory = directory containing mesh files mesh_file,nodes_coords_file,.. = MESH'
     print *,'      output_directory = directory for output files proc***_Databases = defined in Par_file'
     print *
     print *

     open(27,file=trim(LOCAL_PATH)//'/../output_parallel_mesh_decomposer.txt')
     write(27,*)
     write(27,*) '*********************************************************************'
     write(27,*) 'Parallel mesh decomposer, by Vadim Monteiller, CNRS Marseille, France'
     write(27,*) '*********************************************************************'
     write(27,*)
     write(27,*) ' READING MESH FILES '
     call read_mesh_files()
     write(27,*) ' COMPUTE COMPUTATIONAL LOAD of EACH ELEMENT'
     call compute_load_elemnts()
     write(27,*) ' READING Npart in each Cartesian direction'
     call read_tmp_in_file(npartX, npartY, npartZ)
  endif

  call bcast_all_singlei(npartX)
  call bcast_all_singlei(npartY)
  call bcast_all_singlei(npartZ)

  if (sizeprocs /= npartX * npartY * npartZ) stop 'error: nparts **MUST BE** equal to npartX * npartY * npartZ'

  ! 3/ Heuristic mesh partition
  if (myrank == 0) then
     write(27,*) ' DECOMPOSING MESH '
     call decompose_mesh(elmnts_glob, nodes_coords_glob, load_elmnts, nspec_glob, nnodes_glob, npartX, npartY, npartZ)
  endif

  call bcast_all_singlei(nspec_glob)
  if (myrank > 0) allocate(ipart(nspec_glob))
  call bcast_all_i(ipart, nspec_glob)
  call send_partition_mesh_to_all(myrank, ipart, npartX*npartY*npartZ)
  call send_mesh_to_all(myrank)

  ! 4/ write partitionned mesh in one file per proc
  call prepare_database(myrank, elmnts, nspec)
  call write_database(myrank, ipart, elmnts, nodes_coords, elmnts_glob,  mat, mat_prop, undef_mat_prop, &
       count_def_mat, count_undef_mat, ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, &
       ibelm_bottom, ibelm_top, nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, nodes_ibelm_ymax, &
       nodes_ibelm_bottom, nodes_ibelm_top, cpml_to_spec, cpml_regions, is_cpml, ibelm_moho, nodes_ibelm_moho, &
       nspec_glob, nnodes, nspec2D_xmin, nspec2D_xmax,nspec2D_ymin, &
       nspec2D_ymax, nspec2D_bottom, nspec2D_top, nspec_cpml, nspec2D_moho)

  if (myrank == 0) then
     write(27,*) ' END OF MESH DECOMPOSER'
     close(27)
  endif

  ! MPI finish
  call finalize_mpi()

end program xdecompose_mesh_mpi

!-------------------------------------------
!
!
! send partitioned arrays form the master to all other process
!
!
!--------------------------------------------

subroutine send_partition_mesh_to_all(myrank, ipart, npart)

!
!  this should take a long time because we
!  use send from master to each other
!  in sequential way
!

  use module_mesh
  use module_database

  implicit none
  integer,                               intent(in) :: myrank, npart
  integer,   dimension(nspec_glob),      intent(in) :: ipart

  integer                                           :: iE, kE, iE_loc
  integer                                           :: nE
  integer                                           :: irank
  integer                                           :: ivertex, inode, inode_loc
  logical                                           :: not_stored
  integer,           dimension(:),     allocatable  :: nE_irank,  nnodes_in_partition
  integer,           dimension(:,:),   allocatable  :: buffer_to_send
  integer,           dimension(:),     allocatable  :: buffer_to_send1
  integer,           dimension(:),     allocatable  :: buffer_to_send2
  integer,           dimension(:),     allocatable  :: buffer_to_send3
  double precision,  dimension(:,:),   allocatable  :: dp_buffer_to_send
  integer,           dimension(:),     allocatable  :: nelmnts_by_node_glob
  integer,           dimension(:,:),   allocatable  :: elmnts_by_node_glob
  logical,           dimension(:),     allocatable  :: stored_node
   integer,           dimension(1)                  :: nE_loc_for_shut_up_compiler

  ! glob dimension
  call bcast_all_singlei(nspec_glob)
  call bcast_all_singlei(nnodes_glob)
  call bcast_all_singlei(count_def_mat)
  call bcast_all_singlei(count_undef_mat)
  call bcast_all_singlei(nspec2D_xmin)
  call bcast_all_singlei(nspec2D_xmax)
  call bcast_all_singlei(nspec2D_ymin)
  call bcast_all_singlei(nspec2D_ymax)
  call bcast_all_singlei(nspec2D_bottom)
  call bcast_all_singlei(nspec2D_top)
  call bcast_all_singlei(nspec_cpml)
  call bcast_all_singlei(nspec2D_moho)

  nE = nspec_glob

  !! ----------------------------- COMPUTE ELEMENT BY NODE CONNECTIVITY ---------------
  if (myrank == 0) then
     ! evaluate max element by node
     allocate(nelmnts_by_node_glob(nnodes_glob))
     nelmnts_by_node_glob(:)=0
     do iE=1,nE
        do inode=1,NGNOD
           ivertex = elmnts_glob(inode,iE)
           nelmnts_by_node_glob(ivertex) = nelmnts_by_node_glob(ivertex) + 1
        enddo
     enddo

     ! compute elments by node table
     max_elmnts_by_node = maxval(nelmnts_by_node_glob)
     nelmnts_by_node_glob(:)=0
     allocate(elmnts_by_node_glob(max_elmnts_by_node, nnodes_glob))  !!
     elmnts_by_node_glob(:,:)=-1
     do iE=1,nE
        do inode=1,NGNOD
           ivertex = elmnts_glob(inode,iE)
           nelmnts_by_node_glob(ivertex) = nelmnts_by_node_glob(ivertex) + 1
           elmnts_by_node_glob(nelmnts_by_node_glob(ivertex),ivertex) = iE
        enddo
     enddo
  endif
  call bcast_all_singlei(max_elmnts_by_node)

  ! count the number of nodes in each partition
  allocate(nnodes_in_partition(npart))
  allocate(stored_node(npart))
  nnodes_in_partition(:)=0
  if (myrank == 0 ) then
     do inode = 1, nnodes_glob
        stored_node(:)=.false.
        do kE = 1, nelmnts_by_node_glob(inode)
           iE = elmnts_by_node_glob(kE, inode)
           do irank = 0, npart -1
              if (ipart(iE) == irank +1 .and. .not. stored_node(irank+1)) then
                 nnodes_in_partition(irank+1) =  nnodes_in_partition(irank+1) + 1
                 stored_node(irank+1) = .true.
              endif
           enddo
        enddo
     enddo
  endif
  call bcast_all_i(nnodes_in_partition, npart)
  nnodes = nnodes_in_partition(myrank+1)

  ! splitted dual mesh, not really distribued because of use glob numbering
  allocate(elmnts_by_node(max_elmnts_by_node,nnodes), nelmnts_by_node(nnodes))
  elmnts_by_node(:,:)=-1

  ! global to local numbering
  allocate(loc2glob_nodes(nnodes), glob2loc_nodes(nnodes_glob))
  allocate(nodes_coords(3,nnodes))
  nnodes_loc = nnodes
  glob2loc_nodes(:)=-1
  loc2glob_nodes(:)=-1

  if (myrank == 0) then
     irank = 0
     inode_loc = 0
     do inode = 1, nnodes_glob
        not_stored=.true.
        do kE = 1, nelmnts_by_node_glob(inode)
           iE = elmnts_by_node_glob(kE, inode)
           if (ipart(iE) == irank +1 .and. not_stored) then
              inode_loc=inode_loc+1
              elmnts_by_node(:,inode_loc) = elmnts_by_node_glob(:, inode)
              nelmnts_by_node(inode_loc) = nelmnts_by_node_glob(inode)
              glob2loc_nodes(inode) = inode_loc
              loc2glob_nodes(inode_loc) = inode
              nodes_coords(:,inode_loc) = nodes_coords_glob(:,inode)
              not_stored=.false.
           endif
        enddo
     enddo

     do irank = 1, npart - 1

        allocate(buffer_to_send(max_elmnts_by_node,nnodes_in_partition(irank+1)))
        allocate(buffer_to_send1(nnodes_in_partition(irank+1)))
        allocate(buffer_to_send2(nnodes_in_partition(irank+1)))
        allocate(buffer_to_send3(nnodes_glob))
        allocate(dp_buffer_to_send(3,nnodes_in_partition(irank+1)))
        inode_loc = 0
        do inode = 1, nnodes_glob
           not_stored=.true.
           do kE = 1, nelmnts_by_node_glob(inode)
              iE = elmnts_by_node_glob(kE, inode)
              if (ipart(iE) == irank +1 .and. not_stored) then
                 inode_loc=inode_loc+1
                 buffer_to_send(:,inode_loc) = elmnts_by_node_glob(:, inode)
                 buffer_to_send1(inode_loc) = nelmnts_by_node_glob(inode)
                 buffer_to_send2(inode_loc) = inode
                 buffer_to_send3(inode) = inode_loc
                 dp_buffer_to_send(:,inode_loc) = nodes_coords_glob(:,inode)
                 not_stored=.false.
              endif
           enddo
        enddo


        call send_i(buffer_to_send,  max_elmnts_by_node*nnodes_in_partition(irank+1), irank, irank + 9997)
        call send_i(buffer_to_send1, nnodes_in_partition(irank+1), irank, irank + 9987)
        call send_i(buffer_to_send2, nnodes_in_partition(irank+1), irank, irank + 9977)
        call send_i(buffer_to_send3, nnodes_glob, irank, irank + 9967)
        call send_dp(dp_buffer_to_send, 3*nnodes_in_partition(irank+1), irank, irank + 9998)

        deallocate(buffer_to_send)
        deallocate(buffer_to_send1)
        deallocate(buffer_to_send2)
        deallocate(buffer_to_send3)
        deallocate(dp_buffer_to_send)

     enddo

  else

     call recv_i(elmnts_by_node, max_elmnts_by_node*nnodes, 0, myrank  + 9997)
     call recv_i(nelmnts_by_node, nnodes, 0, myrank  + 9987)
     call recv_i(loc2glob_nodes, nnodes, 0, myrank  + 9977)
     call recv_i(glob2loc_nodes, nnodes_glob, 0, myrank  + 9967)
     call recv_dp(nodes_coords, 3*nnodes, 0, myrank + 9998)

  endif




 !!----------------------------- COMPUTE LOC2GLOB_ELMNT and GLOB2LOC_ELMNT -------------------------
  ! count element in my partition
  nE_loc=0
  do iE=1,nE !! loop on all element in partition
     if (ipart(iE) == myrank +1 ) then
        nE_loc = nE_loc + 1
     endif
  enddo

  allocate(loc2glob_elmnt(nE_loc), glob2loc_elmnt(nE))
  glob2loc_elmnt(:)=-1

  ! list of element in my partition
  iE_loc=0
  do iE=1,nE !! loop on all element in partition
     if (ipart(iE) == myrank +1 ) then
        iE_loc = iE_loc + 1
        loc2glob_elmnt(iE_loc) = iE
        glob2loc_elmnt(iE) = iE_loc
     endif
  enddo
  allocate(nE_irank(npart))
  nE_irank(:)=0
  nE_loc_for_shut_up_compiler(1) = nE_loc
  call gather_all_all_i(nE_loc_for_shut_up_compiler, 1, nE_irank, 1, npart)

  nspec=nE_loc !! global varailble to be saved
  allocate(elmnts(NGNOD,nE_loc))
  if ( myrank == 0 ) then

     irank = 0
     iE_loc=0

     do iE = 1, nE
        if (ipart(iE) == irank +1 ) then
           iE_loc = iE_loc + 1
           elmnts(1:NGNOD,iE_loc) = elmnts_glob(1:NGNOD, iE)
        endif
     enddo

     do irank = 1, npart - 1
        allocate(buffer_to_send(NGNOD, nE_irank(irank+1)))
        iE_loc=0
        do iE = 1, nE
           if (ipart(iE) == irank +1 ) then
              iE_loc = iE_loc + 1
              buffer_to_send(1:NGNOD,iE_loc) = elmnts_glob(1:NGNOD, iE)
           endif
        enddo

        call send_i(buffer_to_send, NGNOD*nE_irank(irank+1), irank, irank + 9999)
        deallocate(buffer_to_send)
     enddo

  else
     call recv_i(elmnts, NGNOD*nE_loc, 0, myrank  + 9999)
  endif

  if (myrank == 0 ) then
     deallocate(nodes_coords_glob)
     deallocate(nelmnts_by_node_glob)
  endif

end subroutine send_partition_mesh_to_all

!------------------------------------------------------
!
!
! send whole arrays from master to other
! (it's small arrays that not need to be distributed
! execpt for elmnts_glob)
!
!
!-------------------------------------------------------

subroutine send_mesh_to_all(myrank)

  use module_mesh

  integer,               intent(in) :: myrank
  integer                           :: ier

  if (myrank > 0) then
     allocate(elmnts_glob(NGNOD,nspec_glob),stat=ier)
     if (ier /= 0) then
        write(*,*) 'Error ', myrank, NGNOD,nspec_glob
        stop 'Error allocating array elmnts'
     endif

     allocate(mat(2,nspec_glob),stat=ier)
     if (ier /= 0) stop 'Error allocating array mat'

     allocate(num_material(1:nspec_glob),stat=ier)
     if (ier /= 0) stop 'Error allocating array num_material'

     allocate(mat_prop(16,count_def_mat),stat=ier)
     if (ier /= 0) stop 'Error allocating array mat_prop'

     allocate(undef_mat_prop(6,count_undef_mat),stat=ier)
     if (ier /= 0) stop 'Error allocating array undef_mat_prop'

     allocate(ibelm_xmin(nspec2D_xmin),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_xmin'
     allocate(nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmin'

     allocate(ibelm_xmax(nspec2D_xmax),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_xmax'
     allocate(nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_xmax'

     allocate(ibelm_ymin(nspec2D_ymin),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_ymin'
     allocate(nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymin'

     allocate(ibelm_ymax(nspec2D_ymax),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_ymax'
     allocate(nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_ymax'

     allocate(ibelm_bottom(nspec2D_bottom),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_bottom'
     allocate(nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_bottom'

     allocate(ibelm_top(nspec2D_top),stat=ier)
     if (ier /= 0) stop 'Error allocating array ibelm_top'
     allocate(nodes_ibelm_top(NGNOD2D,nspec2D_top),stat=ier)
     if (ier /= 0) stop 'Error allocating array nodes_ibelm_top'

     allocate(cpml_to_spec(nspec_cpml),stat=ier)
     if (ier /= 0) stop 'Error allocating array CPML_to_spec'
     ! C-PML regions (see below)
     allocate(cpml_regions(nspec_cpml),stat=ier)
     if (ier /= 0) stop 'Error allocating array CPML_regions'

      allocate(is_cpml(nspec_glob),stat=ier)
      if (ier /= 0) stop 'Error allocating array is_CPML'

      allocate(ibelm_moho(nspec2D_moho),stat=ier)
      if (ier /= 0) stop 'Error allocating array ibelm_moho'
      allocate(nodes_ibelm_moho(NGNOD2D,nspec2D_moho),stat=ier)
      if (ier /= 0) stop 'Error allocating array nodes_ibelm_moho'

  endif

  call bcast_all_dp(mat_prop,16*count_def_mat)
  call bcast_all_i(elmnts_glob, NGNOD*nspec_glob)
  call bcast_all_i(mat, 2*nspec_glob)
  call bcast_all_i(num_material, nspec_glob)
  call bcast_all_i(ibelm_xmin,nspec2D_xmin)
  call bcast_all_i(ibelm_xmax,nspec2D_xmax)
  call bcast_all_i(ibelm_ymin,nspec2D_ymin)
  call bcast_all_i(ibelm_ymax,nspec2D_ymax)
  call bcast_all_i(ibelm_bottom,nspec2D_bottom)
  call bcast_all_i(ibelm_top,nspec2D_top)
  call bcast_all_i(ibelm_moho,nspec2D_moho)
  call bcast_all_i(nodes_ibelm_xmin,NGNOD2D*nspec2D_xmin)
  call bcast_all_i(nodes_ibelm_xmax,NGNOD2D*nspec2D_xmax)
  call bcast_all_i(nodes_ibelm_ymin,NGNOD2D*nspec2D_ymin)
  call bcast_all_i(nodes_ibelm_ymax,NGNOD2D*nspec2D_ymax)
  call bcast_all_i(nodes_ibelm_bottom,NGNOD2D*nspec2D_bottom)
  call bcast_all_i(nodes_ibelm_top,NGNOD2D*nspec2D_top)
  call bcast_all_i(nodes_ibelm_moho,NGNOD2D*nspec2D_moho)
  call bcast_all_i(cpml_to_spec, nspec_cpml)
  call bcast_all_i(cpml_regions, nspec_cpml)
  call bcast_all_ch_array(undef_mat_prop,6*count_undef_mat,MAX_STRING_LEN)
  call bcast_all_l_array(is_cpml, nspec_glob)


end subroutine send_mesh_to_all

!--------------------------------
!
! read the desired decomposition in each direction
!
!--------------------------------

subroutine read_tmp_in_file(npartX,npartY,npartZ)

  integer,               intent(inout) :: npartX,npartY,npartZ
  character(len=512)                   :: arg(3)
  integer                              :: i

  do i=1,3
     call get_command_argument(i,arg(i))
  enddo

  read(arg(1),*) npartX
  read(arg(2),*) npartY
  read(arg(3),*) npartZ

end subroutine read_tmp_in_file
