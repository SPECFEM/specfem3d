!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

module part_decompose_mesh_SCOTCH

  implicit none

! Useful kind types
  integer ,parameter :: short = SELECTED_INT_KIND(4), long = SELECTED_INT_KIND(18)

! Number of nodes per elements.
  integer, parameter  :: ESIZE = 8

! Number of faces per element.
  integer, parameter  :: nfaces = 6

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! acoustic-elastic load balancing:
! assumes that elastic at least ~4 times more expensive than acoustic
  integer, parameter :: ACOUSTIC_LOAD = 1
  integer, parameter :: ELASTIC_LOAD = 4

!  include './constants_decompose_mesh_SCOTCH.h'

contains

  !-----------------------------------------------
  ! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
  !-----------------------------------------------
  subroutine mesh2dual_ncommonnodes(nelmnts, nnodes, nsize, sup_neighbour, elmnts,&
                        xadj, adjncy, &
                        nnodes_elmnts, nodes_elmnts, &
                        max_neighbour, ncommonnodes)

    integer(long), intent(in)  :: nelmnts
    integer, intent(in)  :: nnodes
    integer(long), intent(in)  :: nsize
    integer(long), intent(in)  :: sup_neighbour
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts

    integer, dimension(0:nelmnts)  :: xadj
    integer, dimension(0:sup_neighbour*nelmnts-1)  :: adjncy
    integer, dimension(0:nnodes-1)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1)  :: nodes_elmnts
    integer, intent(out) :: max_neighbour
    integer, intent(in)  :: ncommonnodes

    ! local parameters
    integer  :: i, j, k, l, m, nb_edges
    logical  ::  is_neighbour
    integer  :: num_node, n
    integer  :: elem_base, elem_target
    integer  :: connectivity


    ! initializes
    xadj(:) = 0
    adjncy(:) = 0
    nnodes_elmnts(:) = 0
    nodes_elmnts(:) = 0
    nb_edges = 0

    ! list of elements per node
    do i = 0, esize*nelmnts-1
       nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
       nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
    end do

    ! checking which elements are neighbours ('ncommonnodes' criteria)
    do j = 0, nnodes-1
       do k = 0, nnodes_elmnts(j)-1
          do l = k+1, nnodes_elmnts(j)-1

             connectivity = 0
             elem_base = nodes_elmnts(k+j*nsize)
             elem_target = nodes_elmnts(l+j*nsize)
             do n = 1, esize
                num_node = elmnts(esize*elem_base+n-1)
                do m = 0, nnodes_elmnts(num_node)-1
                   if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                      connectivity = connectivity + 1
                   end if
                end do
             end do

             if ( connectivity >=  ncommonnodes) then

                is_neighbour = .false.

                do m = 0, xadj(nodes_elmnts(k+j*nsize))
                   if ( .not.is_neighbour ) then
                      if ( adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour+m) == nodes_elmnts(l+j*nsize) ) then
                         is_neighbour = .true.
                      end if
                   end if
                end do
                if ( .not.is_neighbour ) then
                   adjncy(nodes_elmnts(k+j*nsize)*sup_neighbour &
                          + xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)

                   xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                   if (xadj(nodes_elmnts(k+j*nsize))>sup_neighbour) &
                    stop 'ERROR : too much neighbours per element, modify the mesh.'

                   adjncy(nodes_elmnts(l+j*nsize)*sup_neighbour &
                          + xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)

                   xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                   if (xadj(nodes_elmnts(l+j*nsize))>sup_neighbour) &
                    stop 'ERROR : too much neighbours per element, modify the mesh.'
                end if
             end if
          end do
       end do
    end do

    max_neighbour = maxval(xadj)

    ! making adjacency arrays compact (to be used for partitioning)
    do i = 0, nelmnts-1
       k = xadj(i)
       xadj(i) = nb_edges
       do j = 0, k-1
          adjncy(nb_edges) = adjncy(i*sup_neighbour+j)
          nb_edges = nb_edges + 1
       end do
    end do

    xadj(nelmnts) = nb_edges


  end subroutine mesh2dual_ncommonnodes



  !--------------------------------------------------
  ! construct local numbering for the elements in each partition
  !--------------------------------------------------
  subroutine build_glob2loc_elmnts(nelmnts, part, glob2loc_elmnts,nparts)

    integer(long), intent(in)  :: nelmnts
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(:), pointer  :: glob2loc_elmnts

    integer  :: num_glob, num_part, nparts
    integer, dimension(0:nparts-1)  :: num_loc
    integer :: ier

    ! allocates local numbering array
    allocate(glob2loc_elmnts(0:nelmnts-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array glob2loc_elmnts'

    ! initializes number of local elements per partition
    do num_part = 0, nparts-1
       num_loc(num_part) = 0
    end do

    ! local numbering
    do num_glob = 0, nelmnts-1
       ! gets partition
       num_part = part(num_glob)
       ! increments local numbering of elements (starting with 0,1,2,...)
       glob2loc_elmnts(num_glob) = num_loc(num_part)
       num_loc(num_part) = num_loc(num_part) + 1
    end do


  end subroutine build_glob2loc_elmnts



  !--------------------------------------------------
  ! construct local numbering for the nodes in each partition
  !--------------------------------------------------
  subroutine build_glob2loc_nodes(nelmnts, nnodes, nsize, nnodes_elmnts, nodes_elmnts, part, &
       glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes,nparts)

!    include './constants_decompose_mesh_SCOTCH.h'

    integer(long), intent(in)  :: nelmnts, nsize
    integer, intent(in)  :: nnodes
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(0:nnodes-1), intent(in)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1), intent(in)  :: nodes_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes

    integer  :: num_node
    integer  :: el
    integer  ::  num_part
    integer  ::  size_glob2loc_nodes,nparts
    integer, dimension(0:nparts-1)  :: parts_node
    integer, dimension(0:nparts-1)  :: num_parts
    integer :: ier

    allocate(glob2loc_nodes_nparts(0:nnodes),stat=ier)
    if( ier /= 0 ) stop 'error allocating array glob2loc_nodes_nparts'

    size_glob2loc_nodes = 0
    parts_node(:) = 0

    do num_node = 0, nnodes-1
       glob2loc_nodes_nparts(num_node) = size_glob2loc_nodes
       do el = 0, nnodes_elmnts(num_node)-1
          parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1
       end do

       do num_part = 0, nparts-1
          if ( parts_node(num_part) == 1 ) then
             size_glob2loc_nodes = size_glob2loc_nodes + 1
             parts_node(num_part) = 0
          end if
       end do

    end do

    glob2loc_nodes_nparts(nnodes) = size_glob2loc_nodes

    allocate(glob2loc_nodes_parts(0:glob2loc_nodes_nparts(nnodes)-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array glob2loc_nodes_parts'
    allocate(glob2loc_nodes(0:glob2loc_nodes_nparts(nnodes)-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array glob2loc_nodes'

    glob2loc_nodes(0) = 0

    parts_node(:) = 0
    num_parts(:) = 0
    size_glob2loc_nodes = 0


    do num_node = 0, nnodes-1
       do el = 0, nnodes_elmnts(num_node)-1
          parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1
       end do
       do num_part = 0, nparts-1

          if ( parts_node(num_part) == 1 ) then
             glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
             glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
             size_glob2loc_nodes = size_glob2loc_nodes + 1
             num_parts(num_part) = num_parts(num_part) + 1
             parts_node(num_part) = 0
          end if

       end do
    end do


  end subroutine build_glob2loc_nodes



  !--------------------------------------------------
  ! build interfaces between partitions.
  ! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
  ! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
  ! 5/ second node, if relevant.

  ! interface ignores acoustic and elastic elements

  ! Elements with undefined material are considered as elastic elements.
  !--------------------------------------------------
   subroutine build_interfaces(nelmnts, sup_neighbour, part, elmnts, xadj, adjncy, &
                              tab_interfaces, tab_size_interfaces, ninterfaces, &
                              nparts)

    integer(long), intent(in)  :: nelmnts, sup_neighbour
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts
    integer, dimension(0:nelmnts), intent(in)  :: xadj
    integer, dimension(0:sup_neighbour*nelmnts-1), intent(in)  :: adjncy
    integer, dimension(:),pointer  :: tab_size_interfaces, tab_interfaces
    integer, intent(out)  :: ninterfaces

    integer, intent(in)  :: nparts

    ! local parameters
    integer :: num_part, num_part_bis, el, el_adj, num_interface, num_edge, ncommon_nodes, &
         num_node, num_node_bis
    integer :: i, j
    integer :: ier

    ! counts number of interfaces between partitions
    ninterfaces = 0
    do  i = 0, nparts-1
       do j = i+1, nparts-1
          ninterfaces = ninterfaces + 1
       end do
    end do

    allocate(tab_size_interfaces(0:ninterfaces),stat=ier)
    if( ier /= 0 ) stop 'error allocating array tab_size_interfaces'
    tab_size_interfaces(:) = 0

    num_interface = 0
    num_edge = 0

! determines acoustic/elastic elements based upon given vs velocities
! and counts same elements for each interface
    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                ! looks at all neighbor elements
                do el_adj = xadj(el), xadj(el+1)-1
                   ! adds element if neighbor element lies in next partition
                   if ( part(adjncy(el_adj)) == num_part_bis ) then
                      num_edge = num_edge + 1
                   end if

                end do
             end if
          end do
          ! stores number of elements at interface
          tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge
          num_edge = 0
          num_interface = num_interface + 1

       end do
    end do


! stores element indices for elements from above search at each interface
    num_interface = 0
    num_edge = 0

    allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*7-1)),stat=ier)
    if( ier /= 0 ) stop 'error allocating array tab_interfaces'
    tab_interfaces(:) = 0

    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                do el_adj = xadj(el), xadj(el+1)-1
                   ! adds element if in adjacent partition
                   if ( part(adjncy(el_adj)) == num_part_bis ) then
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+0) = el
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+1) = adjncy(el_adj)
                      ncommon_nodes = 0
                      do num_node = 0, esize-1
                         do num_node_bis = 0, esize-1
                            if ( elmnts(el*esize+num_node) == elmnts(adjncy(el_adj)*esize+num_node_bis) ) then
                               tab_interfaces(tab_size_interfaces(num_interface)*7 &
                                              +num_edge*7+3+ncommon_nodes) &
                                    = elmnts(el*esize+num_node)
                               ncommon_nodes = ncommon_nodes + 1
                            end if
                         end do
                      end do
                      if ( ncommon_nodes > 0 ) then
                         tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+2) = ncommon_nodes
                      else
                         print *, "Error while building interfaces!", ncommon_nodes
                      end if
                      num_edge = num_edge + 1
                   end if
                end do
             end if

          end do
          num_edge = 0
          num_interface = num_interface + 1
       end do
    end do

  end subroutine build_interfaces


  !--------------------------------------------------
  ! build interfaces between partitions.
  ! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
  ! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
  ! 5/ second node, if relevant.

  ! No interface between acoustic and elastic elements.

  ! Elements with undefined material are considered as elastic elements.
  !--------------------------------------------------
   subroutine build_interfaces_no_ac_el_sep(nelmnts, &
                              sup_neighbour, part, elmnts, xadj, adjncy, &
                              tab_interfaces, tab_size_interfaces, ninterfaces, &
                              nb_materials, cs_material, num_material,nparts)

    integer, intent(in)  :: nb_materials,nparts
    integer(long), intent(in)  :: nelmnts, sup_neighbour
    integer, dimension(0:nelmnts-1), intent(in)  :: part
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts
    integer, dimension(0:nelmnts), intent(in)  :: xadj
    integer, dimension(0:sup_neighbour*nelmnts-1), intent(in)  :: adjncy
    integer, dimension(:),pointer  :: tab_size_interfaces, tab_interfaces
    integer, intent(out)  :: ninterfaces
    integer, dimension(1:nelmnts), intent(in)  :: num_material
    ! vs velocities
    double precision, dimension(1:nb_materials), intent(in)  :: cs_material

    ! local parameters
    integer  :: num_part, num_part_bis, el, el_adj, num_interface, num_edge, ncommon_nodes, &
         num_node, num_node_bis
    integer  :: i, j
    logical  :: is_acoustic_el, is_acoustic_el_adj
    integer :: ier

    ! counts number of interfaces between partitions
    ninterfaces = 0
    do  i = 0, nparts-1
       do j = i+1, nparts-1
          ninterfaces = ninterfaces + 1
       end do
    end do

    allocate(tab_size_interfaces(0:ninterfaces),stat=ier)
    if( ier /= 0 ) stop 'error allocating array tab_size_interfaces'
    tab_size_interfaces(:) = 0

    num_interface = 0
    num_edge = 0

! determines acoustic/elastic elements based upon given vs velocities
! and counts same elements for each interface
    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                ! determines whether element is acoustic or not
                if(num_material(el+1) > 0) then
                   if ( cs_material(num_material(el+1)) < TINYVAL) then
                      is_acoustic_el = .true.
                   else
                      is_acoustic_el = .false.
                   end if
                else
                   is_acoustic_el = .false.
                end if
                ! looks at all neighbor elements
                do el_adj = xadj(el), xadj(el+1)-1
                   ! determines whether neighbor element is acoustic or not
                   if(num_material(adjncy(el_adj)+1) > 0) then
                      if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                         is_acoustic_el_adj = .true.
                      else
                         is_acoustic_el_adj = .false.
                      end if
                   else
                      is_acoustic_el_adj = .false.
                   end if
                   ! adds element if neighbor element has same material acoustic/not-acoustic
                   ! and lies in next partition
                   if ( (part(adjncy(el_adj)) == num_part_bis) .and. &
                       (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                      num_edge = num_edge + 1
                   end if
                end do
             end if
          end do
          ! stores number of elements at interface
          tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge
          num_edge = 0
          num_interface = num_interface + 1

       end do
    end do


! stores element indices for elements from above search at each interface
    num_interface = 0
    num_edge = 0

    allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*7-1)),stat=ier)
    if( ier /= 0 ) stop 'error allocating array tab_interfaces'
    tab_interfaces(:) = 0

    do num_part = 0, nparts-1
       do num_part_bis = num_part+1, nparts-1
          do el = 0, nelmnts-1
             if ( part(el) == num_part ) then
                if(num_material(el+1) > 0) then
                   if ( cs_material(num_material(el+1)) < TINYVAL) then
                      is_acoustic_el = .true.
                   else
                      is_acoustic_el = .false.
                   end if
                else
                   is_acoustic_el = .false.
                end if
                do el_adj = xadj(el), xadj(el+1)-1
                   if(num_material(adjncy(el_adj)+1) > 0) then
                      if ( cs_material(num_material(adjncy(el_adj)+1)) < TINYVAL) then
                         is_acoustic_el_adj = .true.
                      else
                         is_acoustic_el_adj = .false.
                      end if
                   else
                      is_acoustic_el_adj = .false.
                   end if
                   if ( (part(adjncy(el_adj)) == num_part_bis) .and. &
                       (is_acoustic_el .eqv. is_acoustic_el_adj) ) then
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+0) = el
                      tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+1) = adjncy(el_adj)
                      ncommon_nodes = 0
                      do num_node = 0, esize-1
                         do num_node_bis = 0, esize-1
                            if ( elmnts(el*esize+num_node) == elmnts(adjncy(el_adj)*esize+num_node_bis) ) then
                               tab_interfaces(tab_size_interfaces(num_interface)*7 &
                                             +num_edge*7+3+ncommon_nodes) &
                                    = elmnts(el*esize+num_node)
                               ncommon_nodes = ncommon_nodes + 1
                            end if
                         end do
                      end do
                      if ( ncommon_nodes > 0 ) then
                         tab_interfaces(tab_size_interfaces(num_interface)*7+num_edge*7+2) = ncommon_nodes
                      else
                         print *, "Error while building interfaces!", ncommon_nodes
                      end if
                      num_edge = num_edge + 1
                   end if
                end do
             end if

          end do
          num_edge = 0
          num_interface = num_interface + 1
       end do
    end do

  end subroutine build_interfaces_no_ac_el_sep



  !--------------------------------------------------
  ! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_glob2loc_nodes_database(IIN_database, iproc, npgeo, &
                  nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                  glob2loc_nodes, nnodes, num_phase)

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: nnodes, iproc, num_phase
    integer, intent(inout)  :: npgeo

    double precision, dimension(3,nnodes)  :: nodes_coords
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes

    integer  :: i, j

    if ( num_phase == 1 ) then
    ! counts number of points in partition
       npgeo = 0
       do i = 0, nnodes-1
          do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
             if ( glob2loc_nodes_parts(j) == iproc ) then
                npgeo = npgeo + 1

             end if

          end do
       end do
    else
    ! writes out point coordinates
       do i = 0, nnodes-1
          do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
             if ( glob2loc_nodes_parts(j) == iproc ) then
                write(IIN_database,*) glob2loc_nodes(j)+1, nodes_coords(1,i+1), &
                                      nodes_coords(2,i+1), nodes_coords(3,i+1)
             end if
          end do
       end do
    end if

  end subroutine Write_glob2loc_nodes_database


  !--------------------------------------------------
  ! Write material properties in the Database
  !--------------------------------------------------
  subroutine write_material_props_database(IIN_database,count_def_mat, &
                                count_undef_mat, mat_prop, undef_mat_prop)

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: count_def_mat,count_undef_mat
    double precision, dimension(6,count_def_mat)  :: mat_prop
    character (len=30), dimension(6,count_undef_mat) :: undef_mat_prop
    integer  :: i

    write(IIN_database,*)  count_def_mat,count_undef_mat
    do i = 1, count_def_mat
      ! database material definition
      !
      ! format:  #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
      !
      ! (note that this order of the properties is different than the input in nummaterial_velocity_file)
      !
       write(IIN_database,*) mat_prop(1,i), mat_prop(2,i), mat_prop(3,i), &
                            mat_prop(4,i), mat_prop(5,i), mat_prop(6,i)
    end do
    do i = 1, count_undef_mat
       write(IIN_database,*) trim(undef_mat_prop(1,i)),' ',trim(undef_mat_prop(2,i)),' ', &
                            trim(undef_mat_prop(3,i)),' ',trim(undef_mat_prop(4,i)),' ', &
                            trim(undef_mat_prop(5,i)),' ',trim(undef_mat_prop(6,i))
    end do

  end subroutine  write_material_props_database


  !--------------------------------------------------
  ! Write elements on boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_boundaries_database(IIN_database, iproc, nelmnts, nspec2D_xmin, nspec2D_xmax, &
                        nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                        ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                        ibelm_ymax, ibelm_bottom, ibelm_top, &
                        nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                        nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                        glob2loc_elmnts, glob2loc_nodes_nparts, &
                        glob2loc_nodes_parts, glob2loc_nodes, part )

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: iproc
    integer(long), intent(in)  :: nelmnts
    integer, intent(in)  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
      nspec2D_ymax, nspec2D_bottom, nspec2D_top
    integer, dimension(nspec2D_xmin), intent(in) :: ibelm_xmin
    integer, dimension(nspec2D_xmax), intent(in) :: ibelm_xmax
    integer, dimension(nspec2D_ymin), intent(in) :: ibelm_ymin
    integer, dimension(nspec2D_ymax), intent(in) :: ibelm_ymax
    integer, dimension(nspec2D_bottom), intent(in) :: ibelm_bottom
    integer, dimension(nspec2D_top), intent(in) :: ibelm_top

    integer, dimension(4,nspec2D_xmin), intent(in) :: nodes_ibelm_xmin
    integer, dimension(4,nspec2D_xmax), intent(in) :: nodes_ibelm_xmax
    integer, dimension(4,nspec2D_ymin), intent(in) :: nodes_ibelm_ymin
    integer, dimension(4,nspec2D_ymax), intent(in) :: nodes_ibelm_ymax
    integer, dimension(4,nspec2D_bottom), intent(in) :: nodes_ibelm_bottom
    integer, dimension(4,nspec2D_top), intent(in) :: nodes_ibelm_top
    integer, dimension(:), pointer :: glob2loc_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    integer, dimension(1:nelmnts)  :: part

    ! local parameters
    integer  :: i,j
    integer  :: loc_node1, loc_node2, loc_node3, loc_node4
    integer  :: loc_nspec2D_xmin,loc_nspec2D_xmax,loc_nspec2D_ymin, &
               loc_nspec2D_ymax,loc_nspec2D_bottom,loc_nspec2D_top

    ! counts number of elements for boundary at xmin, xmax, ymin, ymax, bottom, top in this partition
    loc_nspec2D_xmin = 0
    do i=1,nspec2D_xmin
       if(part(ibelm_xmin(i)) == iproc) then
          loc_nspec2D_xmin = loc_nspec2D_xmin + 1
       end if
    end do
    write(IIN_database,*) 1, loc_nspec2D_xmin
    loc_nspec2D_xmax = 0
    do i=1,nspec2D_xmax
       if(part(ibelm_xmax(i)) == iproc) then
          loc_nspec2D_xmax = loc_nspec2D_xmax + 1
       end if
    end do
    write(IIN_database,*) 2, loc_nspec2D_xmax
    loc_nspec2D_ymin = 0
    do i=1,nspec2D_ymin
       if(part(ibelm_ymin(i)) == iproc) then
          loc_nspec2D_ymin = loc_nspec2D_ymin + 1
       end if
    end do
    write(IIN_database,*) 3, loc_nspec2D_ymin
    loc_nspec2D_ymax = 0
    do i=1,nspec2D_ymax
       if(part(ibelm_ymax(i)) == iproc) then
          loc_nspec2D_ymax = loc_nspec2D_ymax + 1
       end if
    end do
    write(IIN_database,*) 4, loc_nspec2D_ymax
    loc_nspec2D_bottom = 0
    do i=1,nspec2D_bottom
       if(part(ibelm_bottom(i)) == iproc) then
          loc_nspec2D_bottom = loc_nspec2D_bottom + 1
       end if
    end do
    write(IIN_database,*) 5, loc_nspec2D_bottom
    loc_nspec2D_top = 0
    do i=1,nspec2D_top
       if(part(ibelm_top(i)) == iproc) then
          loc_nspec2D_top = loc_nspec2D_top + 1
       end if
    end do
    write(IIN_database,*) 6, loc_nspec2D_top

    ! outputs element index and element node indices
    ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
    !          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
    !          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus
    !          we need to have the arg of glob2loc_elmnts start at 0 ==> glob2loc_nodes(ibelm_** -1 )
    do i=1,nspec2D_xmin
       if(part(ibelm_xmin(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmin(1,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmin(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmin(2,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmin(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmin(3,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmin(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmin(4,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmin(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_xmin(i)-1)+1, &
                                loc_node1, loc_node2, loc_node3, loc_node4
       end if
    end do

    do i=1,nspec2D_xmax
       if(part(ibelm_xmax(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmax(1,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmax(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmax(2,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmax(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmax(3,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmax(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmax(4,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmax(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_xmax(i)-1)+1, &
                                loc_node1, loc_node2, loc_node3, loc_node4
       end if
    end do

    do i=1,nspec2D_ymin
       if(part(ibelm_ymin(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymin(1,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymin(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymin(2,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymin(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymin(3,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymin(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymin(4,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymin(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_ymin(i)-1)+1, &
                                loc_node1, loc_node2, loc_node3, loc_node4
       end if
    end do

    do i=1,nspec2D_ymax
       if(part(ibelm_ymax(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymax(1,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymax(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymax(2,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymax(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymax(3,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymax(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymax(4,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymax(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_ymax(i)-1)+1, &
                               loc_node1, loc_node2, loc_node3, loc_node4
       end if
    end do

    do i=1,nspec2D_bottom
       if(part(ibelm_bottom(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_bottom(1,i)-1), &
                 glob2loc_nodes_nparts(nodes_ibelm_bottom(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_bottom(2,i)-1), &
                 glob2loc_nodes_nparts(nodes_ibelm_bottom(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_bottom(3,i)-1), &
                 glob2loc_nodes_nparts(nodes_ibelm_bottom(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_bottom(4,i)-1), &
                 glob2loc_nodes_nparts(nodes_ibelm_bottom(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_bottom(i)-1)+1, loc_node1, loc_node2, loc_node3, loc_node4
       end if
    end do

    do i=1,nspec2D_top
       if(part(ibelm_top(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_top(1,i)-1), glob2loc_nodes_nparts(nodes_ibelm_top(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_top(2,i)-1), glob2loc_nodes_nparts(nodes_ibelm_top(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_top(3,i)-1), glob2loc_nodes_nparts(nodes_ibelm_top(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_top(4,i)-1), glob2loc_nodes_nparts(nodes_ibelm_top(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_top(i)-1)+1, loc_node1, loc_node2, loc_node3, loc_node4
       end if

    end do

  end subroutine write_boundaries_database


  !--------------------------------------------------
  ! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_partition_database(IIN_database, iproc, nspec, nelmnts, elmnts, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, &
                                      part, num_modele, ngnod, num_phase)

!    include './constants_decompose_mesh_SCOTCH.h'

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: num_phase, iproc
    integer(long), intent(in)  :: nelmnts
    integer, intent(inout)  :: nspec
    integer, dimension(0:nelmnts-1)  :: part
    integer, dimension(0:esize*nelmnts-1)  :: elmnts
    integer, dimension(:), pointer :: glob2loc_elmnts
    integer, dimension(2,nelmnts)  :: num_modele
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    integer, intent(in)  :: ngnod

    integer  :: i,j,k
    integer, dimension(0:ngnod-1)  :: loc_nodes

    if ( num_phase == 1 ) then
    ! counts number of spectral elements in this partition
       nspec = 0
       do i = 0, nelmnts-1
          if ( part(i) == iproc ) then
             nspec = nspec + 1
          end if
       end do

    else
    ! writes out element corner indices
       do i = 0, nelmnts-1
          if ( part(i) == iproc ) then

             do j = 0, ngnod-1
                do k = glob2loc_nodes_nparts(elmnts(i*ngnod+j)), glob2loc_nodes_nparts(elmnts(i*ngnod+j)+1)-1

                   if ( glob2loc_nodes_parts(k) == iproc ) then
                      loc_nodes(j) = glob2loc_nodes(k)
                   end if
                end do

             end do

             ! format:
             ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
             write(IIN_database,*) glob2loc_elmnts(i)+1, num_modele(1,i+1), &
                                  num_modele(2,i+1),(loc_nodes(k)+1, k=0,ngnod-1)
          end if
       end do
    end if


  end subroutine write_partition_database



  !--------------------------------------------------
  ! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_interfaces_database(IIN_database, tab_interfaces, tab_size_interfaces, &
                              iproc, ninterfaces, &
                              my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, &
                              glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                              glob2loc_nodes, num_phase, nparts)

!    include './constants_decompose_mesh_SCOTCH.h'

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: iproc
    integer, intent(in)  :: ninterfaces, nparts
    integer, intent(inout)  :: my_ninterface
    integer, dimension(:), pointer  :: tab_size_interfaces
    integer, dimension(:), pointer  :: tab_interfaces
    integer, dimension(0:ninterfaces-1), intent(inout)  :: my_interfaces
    integer, dimension(0:ninterfaces-1), intent(inout)  :: my_nb_interfaces
    integer, dimension(:), pointer  :: glob2loc_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes

    integer, dimension(4)  :: local_nodes
    integer  :: local_elmnt
    integer  :: num_phase

    integer  :: i, j, k, l
    integer  :: num_interface

    integer  :: count_faces

    num_interface = 0

    if ( num_phase == 1 ) then
    ! counts number of interfaces to neighbouring partitions
       my_interfaces(:) = 0
       my_nb_interfaces(:) = 0

       ! double loops over all partitions
       do i = 0, nparts-1
          do j = i+1, nparts-1
             ! only counts if specified partition (iproc) appears and interface elements increment
             if ( (tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
                  (i == iproc .or. j == iproc) ) then
                ! sets flag
                my_interfaces(num_interface) = 1
                ! sets number of elements on interface
                my_nb_interfaces(num_interface) = tab_size_interfaces(num_interface+1) &
                                            - tab_size_interfaces(num_interface)
             end if
             num_interface = num_interface + 1
          end do
       end do
       my_ninterface = sum(my_interfaces(:))

    else
    ! writes out MPI interface elements
      do i = 0, nparts-1
         do j = i+1, nparts-1
            if ( my_interfaces(num_interface) == 1 ) then
               if ( i == iproc ) then
                  write(IIN_database,*) j, my_nb_interfaces(num_interface)
               else
                  write(IIN_database,*) i, my_nb_interfaces(num_interface)
               end if

               count_faces = 0
               do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
                  if ( i == iproc ) then
                     local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+0))+1
                  else
                     local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+1))+1
                  end if

!!$                  if ( tab_interfaces(k*7+2) == 1 ) then
!!$                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
!!$                          glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
!!$                        if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                           local_nodes(1) = glob2loc_nodes(l)+1
!!$                        end if
!!$                     end do
!!$
!!$                     write(IIN_database,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), -1
!!$                  else
!!$                     if ( tab_interfaces(k*7+2) == 2 ) then
!!$                        do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
!!$                             glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
!!$                           if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                              local_nodes(1) = glob2loc_nodes(l)+1
!!$                           end if
!!$                        end do
!!$                        do l = glob2loc_nodes_nparts(tab_interfaces(k*7+4)), &
!!$                           glob2loc_nodes_nparts(tab_interfaces(k*7+4)+1)-1
!!$                           if ( glob2loc_nodes_parts(l) == iproc ) then
!!$                              local_nodes(2) = glob2loc_nodes(l)+1
!!$                           end if
!!$                        end do
!!$                        write(IIN_database,*) local_elmnt, tab_interfaces(k*7+2), local_nodes(1), local_nodes(2)
!!$                     else
!!$                        write(IIN_database,*) "erreur_write_interface_", tab_interfaces(k*7+2)
!!$                     end if
!!$                  end if
                  select case (tab_interfaces(k*7+2))
                  case (1)
                     ! single point element
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(1) = glob2loc_nodes(l)+1
                        end if
                     end do
                     write(IIN_database,*) local_elmnt, tab_interfaces(k*7+2), &
                                        local_nodes(1), -1, -1, -1
                  case (2)
                     ! edge element
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(1) = glob2loc_nodes(l)+1
                        end if
                     end do
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+4)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+4)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(2) = glob2loc_nodes(l)+1
                        end if
                     end do
                     write(IIN_database,*) local_elmnt, tab_interfaces(k*7+2), &
                                        local_nodes(1), local_nodes(2), -1, -1
                  case (4)
                     ! face element
                     count_faces = count_faces + 1
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(1) = glob2loc_nodes(l)+1
                        end if
                     end do
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+4)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+4)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(2) = glob2loc_nodes(l)+1
                        end if
                     end do
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+5)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+5)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(3) = glob2loc_nodes(l)+1
                        end if
                     end do
                     do l = glob2loc_nodes_nparts(tab_interfaces(k*7+6)), &
                          glob2loc_nodes_nparts(tab_interfaces(k*7+6)+1)-1
                        if ( glob2loc_nodes_parts(l) == iproc ) then
                           local_nodes(4) = glob2loc_nodes(l)+1
                        end if
                     end do
                     write(IIN_database,*) local_elmnt, tab_interfaces(k*7+2), &
                          local_nodes(1), local_nodes(2),local_nodes(3), local_nodes(4)
                  case default
                     print *, "error in write_interfaces_database!", tab_interfaces(k*7+2), iproc
                  end select
               end do

               ! outputs infos
               !print*,'  partition MPI interface:',iproc,num_interface
               !print*,'    element faces: ',count_faces

            end if

            num_interface = num_interface + 1
         end do
      end do

   end if

  end subroutine write_interfaces_database

  !--------------------------------------------------
  ! Write elements on surface boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_moho_surface_database(IIN_database, iproc, nelmnts, &
                        glob2loc_elmnts, glob2loc_nodes_nparts, &
                        glob2loc_nodes_parts, glob2loc_nodes, part, &
                        nspec2D_moho,ibelm_moho,nodes_ibelm_moho)

    integer, intent(in)  :: IIN_database
    integer, intent(in)  :: iproc
    integer(long), intent(in)  :: nelmnts

    integer, dimension(:), pointer :: glob2loc_elmnts
    integer, dimension(:), pointer  :: glob2loc_nodes_nparts
    integer, dimension(:), pointer  :: glob2loc_nodes_parts
    integer, dimension(:), pointer  :: glob2loc_nodes
    integer, dimension(1:nelmnts)  :: part

    integer ,intent(in) :: nspec2D_moho
    integer ,dimension(nspec2D_moho), intent(in) :: ibelm_moho
    integer, dimension(4,nspec2D_moho), intent(in) :: nodes_ibelm_moho

    integer  :: i,j
    integer  :: loc_node1, loc_node2, loc_node3, loc_node4
    integer  :: loc_nspec2D_moho

    ! counts number of elements for moho surface in this partition
    ! optional moho
    loc_nspec2D_moho = 0
    do i=1,nspec2D_moho
       if(part(ibelm_moho(i)) == iproc) then
          loc_nspec2D_moho = loc_nspec2D_moho + 1
       end if
    end do
    ! checks if anything to do
    if( loc_nspec2D_moho == 0 ) return

    ! format: #surface_id, #number of elements
    write(IIN_database,*) 7, loc_nspec2D_moho

    ! outputs element index and element node indices
    ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
    !          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
    !          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus
    !          we need to have the arg of glob2loc_elmnts start at 0 ==> glob2loc_nodes(ibelm_** -1 )

    ! optional moho
    do i=1,nspec2D_moho
       if(part(ibelm_moho(i)) == iproc) then
          do j = glob2loc_nodes_nparts(nodes_ibelm_moho(1,i)-1), glob2loc_nodes_nparts(nodes_ibelm_moho(1,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node1 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_moho(2,i)-1), glob2loc_nodes_nparts(nodes_ibelm_moho(2,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node2 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_moho(3,i)-1), glob2loc_nodes_nparts(nodes_ibelm_moho(3,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node3 = glob2loc_nodes(j)+1
             end if
          end do
          do j = glob2loc_nodes_nparts(nodes_ibelm_moho(4,i)-1), glob2loc_nodes_nparts(nodes_ibelm_moho(4,i))-1
             if (glob2loc_nodes_parts(j) == iproc ) then
                loc_node4 = glob2loc_nodes(j)+1
             end if
          end do
          write(IIN_database,*) glob2loc_elmnts(ibelm_moho(i)-1)+1, &
                              loc_node1, loc_node2, loc_node3, loc_node4
       end if

    end do

  end subroutine write_moho_surface_database



  !--------------------------------------------------
  ! loading : sets weights for acoustic/elastic elements to account for different
  !               expensive calculations in specfem simulations
  !--------------------------------------------------

  subroutine acoustic_elastic_load (elmnts_load,nelmnts,count_def_mat,count_undef_mat, &
                                    num_material,mat_prop,undef_mat_prop)
  !
  ! note:
  !   acoustic material = domainID 1  (stored in mat_prop(6,..) )
  !   elastic material    = domainID 2
  !
    implicit none

    integer(long),intent(in) :: nelmnts
    integer, intent(in)  :: count_def_mat,count_undef_mat

    ! load weights
    integer,dimension(1:nelmnts),intent(out) :: elmnts_load

    ! materials
    integer, dimension(1:nelmnts), intent(in)  :: num_material
    double precision, dimension(6,count_def_mat),intent(in)  :: mat_prop
    character (len=30), dimension(6,count_undef_mat),intent(in) :: undef_mat_prop

    ! local parameters
    logical, dimension(-count_undef_mat:count_def_mat)  :: is_acoustic, is_elastic
    integer  :: i,el,idomain_id

    ! initializes flags
    is_acoustic(:) = .false.
    is_elastic(:) = .false.

    ! sets acoustic/elastic flags for defined materials
    do i = 1, count_def_mat
       idomain_id = mat_prop(6,i)
       ! acoustic material has idomain_id 1
       if (idomain_id == 1 ) then
          is_acoustic(i) = .true.
       endif
       ! elastic material has idomain_id 2
       if (idomain_id == 2 ) then
          is_elastic(i) = .true.
       endif
    enddo

    ! sets acoustic/elastic flags for undefined materials
    do i = 1, count_undef_mat
       read(undef_mat_prop(6,i),*) idomain_id
       ! acoustic material has idomain_id 1
       if (idomain_id == 1 ) then
          is_acoustic(-i) = .true.
       endif
       ! elastic material has idomain_id 2
       if (idomain_id == 2 ) then
          is_elastic(-i) = .true.
       endif
    enddo


    ! sets weights for elements
    do el = 0, nelmnts-1
      ! acoustic element (cheap)
      if ( is_acoustic(num_material(el+1)) ) then
        elmnts_load(el+1) = elmnts_load(el+1)*ACOUSTIC_LOAD
      endif
      ! elastic element (expensive)
      if ( is_elastic(num_material(el+1)) ) then
        elmnts_load(el+1) = elmnts_load(el+1)*ELASTIC_LOAD
      endif
    enddo

  end subroutine acoustic_elastic_load


  !--------------------------------------------------
  ! Repartitioning : two coupled acoustic/elastic elements are transfered to the same partition
  !--------------------------------------------------

  subroutine acoustic_elastic_repartitioning (nelmnts, nnodes, elmnts, &
                        nb_materials, num_material, mat_prop, &
                        sup_neighbour, nsize, &
                        nproc, part, nfaces_coupled, faces_coupled)

    implicit none

    integer(long),intent(in) :: nelmnts
    integer, intent(in)  :: nnodes, nproc, nb_materials
    integer(long), intent(in) :: sup_neighbour,nsize

    integer, dimension(1:nelmnts), intent(in)  :: num_material

    double precision, dimension(6,nb_materials),intent(in)  :: mat_prop

    integer, dimension(0:nelmnts-1)  :: part
    integer, dimension(0:esize*nelmnts-1)  :: elmnts

    integer, intent(out)  :: nfaces_coupled
    integer, dimension(:,:), pointer  :: faces_coupled


    logical, dimension(nb_materials)  :: is_acoustic, is_elastic

    ! neighbors
    integer, dimension(:), allocatable  :: xadj
    integer, dimension(:), allocatable  :: adjncy
    integer, dimension(:), allocatable  :: nnodes_elmnts
    integer, dimension(:), allocatable  :: nodes_elmnts
    integer  :: max_neighbour

    integer  :: i, iface, ier
    integer  :: el, el_adj
    logical  :: is_repartitioned

    ! sets acoustic/elastic flags for materials
    is_acoustic(:) = .false.
    is_elastic(:) = .false.
    do i = 1, nb_materials
       if (mat_prop(6,i) == 1 ) then
          is_acoustic(i) = .true.
       endif
       if (mat_prop(6,i) == 2 ) then
          is_elastic(i) = .true.
       endif
    enddo

    ! gets neighbors by 4 common nodes (face)
    allocate(xadj(0:nelmnts),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xadj'
    allocate(adjncy(0:sup_neighbour*nelmnts-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array adjncy'
    allocate(nnodes_elmnts(0:nnodes-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nnodes_elmnts'
    allocate(nodes_elmnts(0:nsize*nnodes-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_elmnts'
    !call mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts,4)
    call mesh2dual_ncommonnodes(nelmnts, nnodes, nsize, sup_neighbour, &
                                elmnts, xadj, adjncy, nnodes_elmnts, &
                                nodes_elmnts, max_neighbour, 4)

    ! counts coupled elements
    nfaces_coupled = 0
    do el = 0, nelmnts-1
       if ( is_acoustic(num_material(el+1)) ) then
          do el_adj = xadj(el), xadj(el+1) - 1
             if ( is_elastic(num_material(adjncy(el_adj)+1)) ) then
                nfaces_coupled = nfaces_coupled + 1
             endif
          enddo
       endif
    enddo

    ! coupled elements
    allocate(faces_coupled(2,nfaces_coupled),stat=ier)
    if( ier /= 0 ) stop 'error allocating array faces_coupled'

    ! stores elements indices
    nfaces_coupled = 0
    do el = 0, nelmnts-1
       if ( is_acoustic(num_material(el+1)) ) then
          do el_adj = xadj(el), xadj(el+1) - 1
             if ( is_elastic(num_material(adjncy(el_adj)+1)) ) then
                nfaces_coupled = nfaces_coupled + 1
                faces_coupled(1,nfaces_coupled) = el
                faces_coupled(2,nfaces_coupled) = adjncy(el_adj)
             endif
          enddo
       endif
    enddo

    ! puts coupled elements into same partition
    do i = 1, nfaces_coupled*nproc
       is_repartitioned = .false.
       do iface = 1, nfaces_coupled
          if ( part(faces_coupled(1,iface)) /= part(faces_coupled(2,iface)) ) then
             if ( part(faces_coupled(1,iface)) < part(faces_coupled(2,iface)) ) then
                part(faces_coupled(2,iface)) = part(faces_coupled(1,iface))
             else
                part(faces_coupled(1,iface)) = part(faces_coupled(2,iface))
             endif
             is_repartitioned = .true.
          endif
       enddo
       if ( .not. is_repartitioned ) then
          exit
       endif
    enddo

 end subroutine acoustic_elastic_repartitioning

  !--------------------------------------------------
  ! Repartitioning : two coupled moho surface elements are transfered to the same partition
  !--------------------------------------------------

  subroutine moho_surface_repartitioning (nelmnts, nnodes, elmnts, &
                        sup_neighbour, nsize, nproc, part, &
                        nspec2D_moho,ibelm_moho,nodes_ibelm_moho)

    implicit none

    ! number of (spectral) elements  ( <-> nspec )
    integer(long),intent(in) :: nelmnts

    ! number of (global) nodes, number or processes
    integer, intent(in)  :: nnodes, nproc

    ! maximum number of neighours and max number of elements-that-contain-the-same-node
    integer(long), intent(in) :: sup_neighbour,nsize

    ! partition index on each element
    integer, dimension(0:nelmnts-1)  :: part

    ! mesh element indexing
    ! ( elmnts(esize,nspec) )
    integer, dimension(0:esize*nelmnts-1)  :: elmnts

    ! moho surface
    integer ,intent(in) :: nspec2D_moho
    integer ,dimension(nspec2D_moho), intent(in) :: ibelm_moho
    integer, dimension(4,nspec2D_moho), intent(in) :: nodes_ibelm_moho

    ! local parameters
    integer :: nfaces_coupled
    integer, dimension(:,:), pointer  :: faces_coupled

    logical, dimension(:),allocatable  :: is_moho,node_is_moho

    ! for neighbors
    integer, dimension(:), allocatable  :: xadj
    integer, dimension(:), allocatable  :: adjncy
    integer, dimension(:), allocatable  :: nnodes_elmnts
    integer, dimension(:), allocatable  :: nodes_elmnts
    integer  :: max_neighbour

    integer  :: i, j, iface, inode, ispec2D, counter
    integer  :: el, el_adj
    logical  :: is_repartitioned
    integer :: ier

    ! temporary flag arrays
    ! element ids start from 0
    allocate( is_moho(0:nelmnts-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array is_moho'
    ! node ids start from 0
    allocate( node_is_moho(0:nnodes-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array node_is_moho'

    is_moho(:) = .false.
    node_is_moho(:) = .false.

    ! sets moho flags for known elements
    do ispec2D = 1, nspec2D_moho
      ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
      el = ibelm_moho(ispec2D) - 1
      is_moho(el) = .true.

      ! sets node flags
      do j=1,4
        ! note: assumes that node indices in nodes_ibelm_* arrays are in the range from 1 to nodes
        inode = nodes_ibelm_moho(j,ispec2D) - 1
        node_is_moho(inode) = .true.
      enddo
    enddo

    ! checks if element has moho surface
    do el = 0, nelmnts-1
      if( is_moho(el) ) cycle

      ! loops over all element corners
      counter = 0
      do i=0,esize-1
        ! note: assumes that node indices in elmnts array are in the range from 0 to nodes-1
        inode = elmnts(el*esize+i)
        if( node_is_moho(inode) ) counter = counter + 1
      enddo

      ! sets flag if it has a surface
      if( counter == 4 ) is_moho(el) = .true.
    enddo

    ! statistics output
    counter = 0
    do el=0, nelmnts-1
     if ( is_moho(el) ) counter = counter + 1
    enddo
    print*,'  moho elements = ',counter

    ! gets neighbors by 4 common nodes (face)
    ! contains number of adjacent elements (neighbours)
    allocate(xadj(0:nelmnts),stat=ier)
    if( ier /= 0 ) stop 'error allocating array xadj'
    ! contains all element id indices of adjacent elements
    allocate(adjncy(0:sup_neighbour*nelmnts-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array adjncy'
    allocate(nnodes_elmnts(0:nnodes-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nnodes_elmnts'
    allocate(nodes_elmnts(0:nsize*nnodes-1),stat=ier)
    if( ier /= 0 ) stop 'error allocating array nodes_elmnts'

    call mesh2dual_ncommonnodes(nelmnts, nnodes, nsize, sup_neighbour, &
                        elmnts, xadj, adjncy, nnodes_elmnts, &
                        nodes_elmnts, max_neighbour, 4)

    ! counts coupled elements
    nfaces_coupled = 0
    do el = 0, nelmnts-1
       if ( is_moho(el) ) then
          do el_adj = xadj(el), xadj(el+1) - 1
            ! increments counter if it contains face
            if( is_moho(adjncy(el_adj)) ) nfaces_coupled = nfaces_coupled + 1
          enddo
       endif
    enddo

    ! coupled elements
    allocate(faces_coupled(2,nfaces_coupled),stat=ier)
    if( ier /= 0 ) stop 'error allocating array faces_coupled'

    ! stores elements indices
    nfaces_coupled = 0
    do el = 0, nelmnts-1
       if ( is_moho(el) ) then
          do el_adj = xadj(el), xadj(el+1) - 1
             if ( is_moho(adjncy(el_adj)) ) then
                nfaces_coupled = nfaces_coupled + 1
                faces_coupled(1,nfaces_coupled) = el
                faces_coupled(2,nfaces_coupled) = adjncy(el_adj)
             endif
          enddo
       endif
    enddo

    ! puts coupled elements into same partition
    do i = 1, nfaces_coupled*nproc
       is_repartitioned = .false.
       do iface = 1, nfaces_coupled
          if ( part(faces_coupled(1,iface)) /= part(faces_coupled(2,iface)) ) then
             ! coupled moho elements are in different partitions
             if ( part(faces_coupled(1,iface)) < part(faces_coupled(2,iface)) ) then
                part(faces_coupled(2,iface)) = part(faces_coupled(1,iface))
             else
                part(faces_coupled(1,iface)) = part(faces_coupled(2,iface))
             endif
             is_repartitioned = .true.
          endif
       enddo
       if ( .not. is_repartitioned ) then
          exit
       endif
    enddo

 end subroutine moho_surface_repartitioning


end module part_decompose_mesh_SCOTCH

