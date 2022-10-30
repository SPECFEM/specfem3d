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


  subroutine check_valence()

! checks valence of nodes

  use decompose_mesh_par

  implicit none

  ! local paramters
  logical, dimension(:), allocatable :: mask_nodes_elmnts
  integer, dimension(:), allocatable :: used_nodes_elmnts
  integer :: max_neighbor,inode,ispec,id,ier

  ! number of faces per element.
  integer, parameter  :: nfaces = 6

  ! allocate temporary array
  allocate(used_nodes_elmnts(nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 107')
  if (ier /= 0) stop 'Error allocating array used_nodes_elmnts'

  allocate(mask_nodes_elmnts(nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array mask_nodes_elmnts')
  if (ier /= 0) stop 'error allocating array mask_nodes_elmnts'

  used_nodes_elmnts(:) = 0
  mask_nodes_elmnts(:) = .false.

  ! loops over all element nodes
  do ispec = 1, nspec
    do inode = 1, NGNOD
      ! node index (between 1 and nnodes)
      id = elmnts(inode,ispec)
      if ( id < 1 .or. id > nnodes ) stop 'Error node index exceeds bounds in check_valence'

      ! masks common nodes
      mask_nodes_elmnts(id) = .true.

      ! counts how many elements share this node
      used_nodes_elmnts(id) = used_nodes_elmnts(id) + 1
    enddo
  enddo

  ! user output
  print *, 'node valence:'
  print *, '  min = ',minval(used_nodes_elmnts(:)),' max = ', maxval(used_nodes_elmnts(:))

  ! checks nodes
  do inode = 1, nnodes
    if (.not. mask_nodes_elmnts(inode)) then
      print *,'Error: node ',inode,' has valence ',used_nodes_elmnts(inode),' mask ',mask_nodes_elmnts(inode)
      stop 'ERROR: found some unused nodes (weird, but not necessarily fatal; your mesher may have created extra nodes&
              & or your mesh contains HEX27 elements while NGNOD in Par_file is set to 8).'
    endif
  enddo

  ! max number of elements that contain the same node
  nsize = maxval(used_nodes_elmnts(:))

  ! frees temporary arrays
  deallocate(mask_nodes_elmnts,used_nodes_elmnts)

! majoration (overestimate) of the maximum number of neighbors per element
!! DK DK nfaces is a constant equal to 6 (number of faces of a cube).
!! DK DK I have no idea how this formula works; it was designed by Nicolas Le Goff
  sup_neighbor = NGNOD_EIGHT_CORNERS * nsize - (NGNOD_EIGHT_CORNERS + (NGNOD_EIGHT_CORNERS/2 - 1)*nfaces)

  ! checks that no underestimation
  if (sup_neighbor < nsize) sup_neighbor = nsize

  ! counts maximum number of possible neighbor elements
  call get_max_neighbor(max_neighbor)

  ! user output
  print *, 'neighbors:'
  print *, '  nsize = ',nsize
  print *, '  valence: sup_neighbor = ', sup_neighbor, 'max_neighbor = ',max_neighbor

  if (max_neighbor > sup_neighbor) stop 'found max_neighbor > sup_neighbor in check_valence()'

  ! sets supneighbor to actual number to avoid allocating too much memory for larger meshes
  sup_neighbor = max_neighbor

  end subroutine check_valence

!
!---------------------------------------------------------------------------------------------
!

  subroutine get_max_neighbor(max_neighbor)

! checks maximum number of neighbors

  use decompose_mesh_par, only: NGNOD,NGNOD_EIGHT_CORNERS,nspec,nnodes,nsize,elmnts

  implicit none

  integer, intent(out) :: max_neighbor

  ! local parameters
  integer, dimension(nnodes) :: nnodes_elmnts
  integer, dimension(nsize*nnodes) :: nodes_elmnts

  integer :: ispec,inode,id
  integer :: k,m
  integer :: elem_target

  ! note: a regular hex-mesh would have a maximum node valence of 8, and a maximum of 26 neighbors
  !       for a squeezed mesh, this might be higher
  !
  ! nsize = max number of elements that contain the same node
  integer :: elem_neighbors(NGNOD_EIGHT_CORNERS * nsize)
  integer :: num_neighbors

  logical :: is_neighbor

  ! initializes
  max_neighbor = 0

  ! builds list of elements per node
  nnodes_elmnts(:) = 0
  nodes_elmnts(:) = 0
  do ispec = 1, nspec
    do inode = 1, NGNOD
      ! node index (between 1 and nnodes)
      id = elmnts(inode,ispec)
      if ( id < 1 .or. id > nnodes ) stop 'error node index exceeds bounds in get_MAX_NEIGHBORS'

      ! sets element into node list
      nodes_elmnts( ((id-1)*nsize + 1) + nnodes_elmnts(id)) = ispec
      ! increases element counter
      nnodes_elmnts(id) = nnodes_elmnts(id) + 1
    enddo
  enddo

  ! counts all neighbors which can be connected through a single node
  ! loops over all element nodes
  do ispec = 1, nspec
    ! list of all neighbors
    elem_neighbors(:) = 0
    num_neighbors = 0

    ! loops over all nodes of element
    do inode = 1, NGNOD
      ! node index ( from 0 to nnodes-1)
      id = elmnts(inode,ispec)

      ! finds all neighbors
      ! loops over all elements assigned to this node
      do k = 0,nnodes_elmnts(id) - 1
        ! element associated with this node
        elem_target = nodes_elmnts( ((id-1)*nsize + 1) + k )
        if ( elem_target < 1 .or. elem_target > nspec ) stop 'error elem_target exceeds bounds'

        ! skips reference to own element
        if ( elem_target == ispec ) cycle

        ! searches element in neighbors list
        is_neighbor = .false.
        do m = 1,num_neighbors
          if ( elem_neighbors(m) == elem_target) then
            is_neighbor = .true.
            exit
          endif
        enddo
        ! adds as a new neighbor
        if (.not. is_neighbor ) then
          ! increases count
          num_neighbors = num_neighbors + 1
          if ( num_neighbors > NGNOD_EIGHT_CORNERS * nsize) stop 'error num_neigbors exceeds bounds'
          ! inserts into temporary list
          elem_neighbors(num_neighbors) = elem_target
        endif
      enddo
    enddo

    ! sets maximum number of neighbors
    if ( num_neighbors > max_neighbor) max_neighbor = num_neighbors

  enddo

  end subroutine get_max_neighbor

