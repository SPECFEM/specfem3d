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

module module_database

  use shared_parameters, only: NGNOD, NGNOD2D, LOCAL_PATH, MSL => MAX_STRING_LEN

  integer                                     :: nE_loc
  integer, dimension(:),  allocatable         :: loc2glob_elmnt
  integer, dimension(:),  allocatable         :: glob2loc_elmnt

  integer                                     :: nnodes_loc
  integer, dimension(:),  allocatable         :: loc2glob_nodes
  integer, dimension(:),  allocatable         :: glob2loc_nodes

  integer                                     :: size_adjacency
  integer, dimension(:),  allocatable         :: adjcy, id_adjcy

  integer                                     :: max_elmnts_by_node
  integer, dimension(:),  allocatable         :: nelmnts_by_node
  integer, dimension(:,:),allocatable         :: elmnts_by_node

  integer, private, parameter                 :: IIN_database=49
  integer, private                            :: ier
  integer, private, dimension(:), allocatable :: node_loc

contains

!----------------------------------
!
!  set up the arrays for database
!
!----------------------------------

  subroutine prepare_database(myrank,  elmnts, nE)

    implicit none
    integer,                         intent(in) :: myrank , nE
    integer,   dimension(NGNOD,nE),  intent(in) :: elmnts

    call compute_adjcy_table(myrank,  elmnts, nE)
    if (myrank == 0) write(27,*) ' COMPUTE ADJACENCY '
    allocate(node_loc(NGNOD2D))

  end subroutine prepare_database

!----------------------------------
!
! write database for xgenerate_database
!
!----------------------------------

  subroutine write_database(myrank, ipart, elmnts, nodes_coords, elmnts_glob,  num_modele,  mat_prop, &
       undef_mat_prop, count_def_mat, count_undef_mat, ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, &
       ibelm_bottom, ibelm_top, nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, nodes_ibelm_ymax, &
       nodes_ibelm_bottom, nodes_ibelm_top, cpml_to_spec, cpml_regions, is_cpml, ibelm_moho, nodes_ibelm_moho, &
       nE, nnodes, nspec2D_xmin, nspec2D_xmax,nspec2D_ymin, &
       nspec2D_ymax, nspec2D_bottom, nspec2D_top, nspec_cpml, nspec2D_moho)

    implicit none

    integer,                                              intent(in)  :: myrank
    integer,                                              intent(in)  :: nnodes, nE
    integer,                                              intent(in)  :: count_def_mat
    integer,                                              intent(in)  :: count_undef_mat
    integer,                                              intent(in)  :: nspec2D_xmin, nspec2D_xmax
    integer,                                              intent(in)  :: nspec2D_ymin, nspec2D_ymax
    integer,                                              intent(in)  :: nspec2D_bottom, nspec2D_top
    integer,                                              intent(in)  :: nspec_cpml, nspec2D_moho
    integer,            dimension(nE),                    intent(in)  :: ipart
    integer,            dimension(NGNOD,nE),              intent(in)  :: elmnts_glob
    integer,            dimension(NGNOD,nE_loc),          intent(in)  :: elmnts
    integer,            dimension(2,nE),                  intent(in)  :: num_modele
    integer,            dimension(nspec2D_xmin),          intent(in)  :: ibelm_xmin
    integer,            dimension(nspec2D_xmax),          intent(in)  :: ibelm_xmax
    integer,            dimension(nspec2D_ymin),          intent(in)  :: ibelm_ymin
    integer,            dimension(nspec2D_ymax),          intent(in)  :: ibelm_ymax
    integer,            dimension(nspec2D_bottom),        intent(in)  :: ibelm_bottom
    integer,            dimension(nspec2D_top),           intent(in)  :: ibelm_top
    integer,            dimension(nspec2D_moho),          intent(in)  :: ibelm_moho
    integer,            dimension(NGNOD2D,nspec2D_xmin),  intent(in)  :: nodes_ibelm_xmin
    integer,            dimension(NGNOD2D,nspec2D_xmax),  intent(in)  :: nodes_ibelm_xmax
    integer,            dimension(NGNOD2D,nspec2D_ymin),  intent(in)  :: nodes_ibelm_ymin
    integer,            dimension(NGNOD2D,nspec2D_ymax),  intent(in)  :: nodes_ibelm_ymax
    integer,            dimension(NGNOD2D,nspec2D_bottom),intent(in)  :: nodes_ibelm_bottom
    integer,            dimension(NGNOD2D,nspec2D_top),   intent(in)  :: nodes_ibelm_top
    integer,            dimension(NGNOD2D,nspec2D_moho),  intent(in)  :: nodes_ibelm_moho
    integer,            dimension(nspec_cpml),            intent(in)  :: cpml_to_spec, cpml_regions
    double precision,   dimension(16,count_def_mat),      intent(in)  :: mat_prop
    character(len=MSL), dimension(6,count_undef_mat),     intent(in)  :: undef_mat_prop
    double precision,   dimension(3,nnodes),              intent(in)  :: nodes_coords
    logical,            dimension(nE),                    intent(in)  :: is_cpml

    character(len=20)                                                 :: prname
    integer                                                           :: i1, i2, in1, in2
    integer                                                           :: inode, i, k
    integer                                                           :: islice, kE
    integer                                                           :: iE, iE_loc, iE_n
    integer                                                           :: nspec_cpml_local
    integer                                                           :: loc_nspec2D_moho
    integer                                                           :: npart, nb_stored_slice
    integer                                                           :: MAX_ELEMENT_B_PARTITION
    integer                                                           :: INUM_NEIGH_PART
    integer,            dimension(:),       allocatable               :: islice_neigh, istored, liste_comm_nodes
    integer,            dimension(:),       allocatable               :: num_element_in_boundary_partition
    integer,            dimension(:),       allocatable               :: ie_bnd_stored, loc_elmnt
    integer,            dimension(:,:,:),   allocatable               :: my_interfaces_ext_mesh

    integer :: ncp

    ! opens output file
    write(prname, "(i6.6,'_Database')") myrank
    open(unit=IIN_database,file=LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/proc'//prname, &
         status='unknown', action='write', form='unformatted', iostat = ier)
    if (ier /= 0) then
       print *,'Error file open:',LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/proc'//prname
       print *
       print *,'check if path exists:',LOCAL_PATH(1:len_trim(LOCAL_PATH))
       stop 'Error file open Database'
    endif

    ! write nodes coordinates located in my partition -----
    write(IIN_database) nnodes_loc
    do inode = 1, nnodes_loc
       k = inode !loc2glob_nodes(inode)
       write(IIN_database) inode,nodes_coords(1,k),nodes_coords(2,k),nodes_coords(3,k)
    enddo

    ! write material properties in my partition -----
    write(IIN_database)  count_def_mat,count_undef_mat
    do i=1, count_def_mat
       write(IIN_database) mat_prop(1:16,i)
    enddo
    do i = 1, count_undef_mat
       write(IIN_database) undef_mat_prop(1:6,i)
    enddo

    ! write element connectivity in my partition -----
    write(IIN_database) nE_loc
    allocate(loc_elmnt(NGNOD))
    do iE_loc = 1, nE_loc
       iE = iE_loc ! loc2glob_elmnt(iE_loc)
       do inode =1, NGNOD
          loc_elmnt(inode) = glob2loc_nodes(elmnts(inode, iE))
       enddo
       iE = loc2glob_elmnt(iE_loc)
       write(IIN_database) iE_loc, num_modele(1:2,iE), loc_elmnt(1:NGNOD)
    enddo

    ! counnt element boundary in my partition -----
    call count_my_boundary(myrank, ipart, ibelm_xmin,   nspec2D_xmin,   nE, 1)
    call count_my_boundary(myrank, ipart, ibelm_xmax,   nspec2D_xmax,   nE, 2)
    call count_my_boundary(myrank, ipart, ibelm_ymin,   nspec2D_ymin,   nE, 3)
    call count_my_boundary(myrank, ipart, ibelm_ymax,   nspec2D_ymax,   nE, 4)
    call count_my_boundary(myrank, ipart, ibelm_bottom, nspec2D_bottom, nE, 5)
    call count_my_boundary(myrank, ipart, ibelm_top,    nspec2D_top,    nE, 6)

    ! write boundary element in my partition -----
    call write_my_boundary(myrank, ipart, ibelm_xmin,   nodes_ibelm_xmin,   nspec2D_xmin,   nE)
    call write_my_boundary(myrank, ipart, ibelm_xmax,   nodes_ibelm_xmax,   nspec2D_xmax,   nE)
    call write_my_boundary(myrank, ipart, ibelm_ymin,   nodes_ibelm_ymin,   nspec2D_ymin,   nE)
    call write_my_boundary(myrank, ipart, ibelm_ymax,   nodes_ibelm_ymax,   nspec2D_ymax,   nE)
    call write_my_boundary(myrank, ipart, ibelm_bottom, nodes_ibelm_bottom, nspec2D_bottom, nE)
    call write_my_boundary(myrank, ipart, ibelm_top,    nodes_ibelm_top,    nspec2D_top,    nE)

    ! write CPML elements  -----
    ! writes number of C-PML elements in the global mesh
    write(IIN_database) nspec_cpml

    if (nspec_cpml > 0) then

       ! writes number of C-PML elements in this partition
       nspec_cpml_local = 0
       do i=1,nspec_cpml
          if (ipart(CPML_to_spec(i)) == myrank +1) then
             nspec_cpml_local = nspec_cpml_local + 1
          endif
       enddo

       write(IIN_database) nspec_cpml_local

       ! writes C-PML regions and C-PML spectral elements global indexing
       do i=1,nspec_cpml
          ! #id_cpml_regions = 1 : X_surface C-PML
          ! #id_cpml_regions = 2 : Y_surface C-PML
          ! #id_cpml_regions = 3 : Z_surface C-PML
          ! #id_cpml_regions = 4 : XY_edge C-PML
          ! #id_cpml_regions = 5 : XZ_edge C-PML
          ! #id_cpml_regions = 6 : YZ_edge C-PML
          ! #id_cpml_regions = 7 : XYZ_corner C-PML
          !
          ! format: #id_cpml_element #id_cpml_regions
          if (ipart(CPML_to_spec(i)) == myrank+1) then
             write(IIN_database) glob2loc_elmnt(CPML_to_spec(i)), CPML_regions(i)
          endif
       enddo

       ! writes mask of C-PML elements for all elements in this partition
       ncp=0
       do i=1,nE
          if (ipart(i) == myrank+1) then
             write(IIN_database) is_CPML(i)
             if (is_CPML(i)) ncp=ncp+1
          endif
       enddo
       !write(*,*) 'CPML ', myrank, ncp, nspec_cpml_local, nspec_cpml
    endif


    ! write MPI interfaces  -----
    npart=maxval(ipart(:))
    allocate(istored(npart),islice_neigh(npart))
    allocate(num_element_in_boundary_partition(npart))
    istored(:)=0
    islice_neigh(:)=0
    num_element_in_boundary_partition(:)=0
    nb_stored_slice=0
    do iE = 1, nE  !! loop over all elements
       if (ipart(iE) == myrank+1) then  !! in my partition
          iE_loc = glob2loc_elmnt(iE)
          !write(*,*) myrank, iE_loc
          do i = id_adjcy(iE_loc-1)+1,  id_adjcy(iE_loc)    !! loop on all my neigbhours
             iE_n = adjcy(i)
             k = 0
             if (ipart(iE_n) /= myrank+1) then      !! it's not my element
                if (istored(ipart(iE_n)) == 0 ) then
                   nb_stored_slice = nb_stored_slice + 1
                   istored(ipart(iE_n)) =  nb_stored_slice
                   islice_neigh(nb_stored_slice) = ipart(iE_n) !!istored(ipart(iE_n))
                endif
                num_element_in_boundary_partition(ipart(iE_n)) =  num_element_in_boundary_partition(ipart(iE_n)) + 1
             endif
          enddo
       endif
    enddo

    max_element_b_partition = maxval(num_element_in_boundary_partition)
    !write(*,*) myrank, max_element_b_partition,nb_stored_slice
    allocate( my_interfaces_ext_mesh(6,max_element_b_partition,nb_stored_slice))
    allocate(liste_comm_nodes(NGNOD2D), ie_bnd_stored(npart))
    ie_bnd_stored(:)=0
    my_interfaces_ext_mesh(:,:,:)=-99
    do iE = 1, nE  !! loop over all elements
       if (ipart(iE) == myrank+1) then !! in my partition
          iE_loc = glob2loc_elmnt(iE)
          do i = id_adjcy(iE_loc-1)+1,  id_adjcy(iE_loc)    !! loop on all my neigbhours
             iE_n = adjcy(i)
             k = 0
             inum_neigh_part = istored(ipart(iE_n))
             liste_comm_nodes(:) =  -1
             if (ipart(iE_n) /= myrank+1) then    !! it's not my element
                do i1 = 1, NGNOD                    !! loop over edges on my element
                   in1 = elmnts(i1,iE_loc)
                   do i2 = 1, NGNOD                 !! loop over edges on my neighbour element
                      in2 = elmnts_glob(i2,iE_n)    !! pb on utilise elmnts_glob
                      if (in1 == in2 ) then         !! it's common edge
                         k = k + 1
                         liste_comm_nodes(k) = glob2loc_nodes(in1)  !! store it in local numbering
                      endif
                   enddo
                enddo
             endif  !! on doit inclure ce if dessous ?
             if (k > 0) then
                ie_bnd_stored(inum_neigh_part) = ie_bnd_stored(inum_neigh_part) + 1
                kE = ie_bnd_stored(inum_neigh_part)
                my_interfaces_ext_mesh(1,kE ,inum_neigh_part) = glob2loc_elmnt(iE)
                my_interfaces_ext_mesh(2,kE ,inum_neigh_part) = k
                my_interfaces_ext_mesh(3:6,kE ,inum_neigh_part) = liste_comm_nodes(1:NGNOD2D)
             endif
          enddo
       endif
    enddo

    write(IIN_database) nb_stored_slice, max_element_b_partition
    do islice = 1, nb_stored_slice
       !write(27,*) 'ii ' , myrank, islice_neigh(islice)-1, num_element_in_boundary_partition(islice_neigh(islice))
       write(IIN_database) islice_neigh(islice)-1, num_element_in_boundary_partition(islice_neigh(islice))
       do iE = 1, num_element_in_boundary_partition(islice_neigh(islice))
          write(IIN_database)  my_interfaces_ext_mesh(1:6, iE, islice)
       enddo
    enddo

    ! write MOHO
    ! optional moho
    loc_nspec2D_moho = 0
    do i=1,nspec2D_moho
       if (ipart(ibelm_moho(i)) == myrank+1) then
          loc_nspec2D_moho = loc_nspec2D_moho + 1
       endif
    enddo
    ! checks if anything to do
    if (loc_nspec2D_moho > 0) then

       ! format: #surface_id, #number of elements
       write(IIN_database) 7, loc_nspec2D_moho

       ! outputs element index and element node indices
       ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
       !          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
       !          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus
       !          we need to have the arg of glob2loc_elmnts start at 0, and thus we use glob2loc_nodes(ibelm_** -1)

       ! optional moho
       do i=1,nspec2D_moho
          if (ipart(ibelm_moho(i)) == myrank+1) then
             do inode = 1,NGNOD2D
                node_loc(inode) = glob2loc_nodes(nodes_ibelm_moho(inode,i))
             enddo
             write(IIN_database) glob2loc_elmnt(ibelm_moho(i)), node_loc(1:NGNOD2D)
          endif
       enddo
    endif
    close(IIN_database)

  end subroutine write_database

!----------------------------------
!
!  count elements in boundary
!
!----------------------------------

  subroutine count_my_boundary(myrank, ipart, ibelm, nspec2D, nE, iflag)

    implicit none
    integer,                                   intent(in) :: myrank, iflag
    integer,                                   intent(in) :: nspec2D, nE
    integer,   dimension(nE),                  intent(in) :: ipart
    integer,   dimension(nspec2D),             intent(in) :: ibelm

    integer ::  i, nspec2D_loc, iE

    nspec2D_loc = 0
    do i=1, nspec2D
       iE=ibelm(i)
       if (ipart(iE) == myrank +1) nspec2D_loc = nspec2D_loc + 1
    enddo
    write(IIN_database) iflag, nspec2D_loc

  end subroutine count_my_boundary

!----------------------------------
!
!  write elements in boundary
!
!----------------------------------

  subroutine  write_my_boundary(myrank, ipart, ibelm, nodes_ibelm, nspec2D, nE)
    implicit none

    integer,                                   intent(in) :: myrank, nspec2D, nE
    integer,   dimension(nE),                  intent(in) :: ipart
    integer,   dimension(nspec2D),             intent(in) :: ibelm
    integer,   dimension(NGNOD2D,nspec2D),     intent(in) :: nodes_ibelm

    integer :: i,  iE, iE_loc, inode

    do i=1, nspec2D
       iE=ibelm(i)
       iE_loc=glob2loc_elmnt(iE)
       do inode = 1, NGNOD2D
          node_loc(inode) = glob2loc_nodes(nodes_ibelm(inode, i))
       enddo
       if (ipart(iE) == myrank +1) write(IIN_database) iE_loc, node_loc(1:NGNOD2D)
    enddo

  end subroutine write_my_boundary

!----------------------------------
!
!  distributed adjacency table
!
!----------------------------------

  subroutine compute_adjcy_table(myrank,  elmnts, nE)

    implicit none

    integer,                         intent(in) :: myrank , nE
    integer,   dimension(NGNOD,nE),  intent(in) :: elmnts
    integer,   parameter                        :: NGNOD_EIGHT_CORNERS=8
    integer                                     :: iE, iE_loc
    integer                                     :: kE, New_element
    integer                                     :: inode, ivertex
    integer                                     :: ivertex_local
    integer                                     :: nb_element_stored
    integer                                     :: max_element_to_store
    integer,   dimension(:),        allocatable :: stored_elements
    integer,   dimension(:),        allocatable :: nb_neigh


    !! --------------- COMPUTE ADJACENCY TABLE ------------------------------------------
    max_element_to_store = 8 * max_elmnts_by_node  ! over estimate the number of
                                                   ! neighbors element
    ! evalute the size of the adjacency table
    allocate(stored_elements(max_element_to_store))
    allocate(nb_neigh(nE_loc))
    do iE_loc=1,nE_loc !! loop on all element in partition
       iE = loc2glob_elmnt(iE_loc)
       nb_element_stored = 0
       stored_elements(:)=0
       do inode = 1, NGNOD_EIGHT_CORNERS  !! loop only on the corner of the element
          ivertex = elmnts(inode,iE_loc)
          ivertex_local=glob2loc_nodes(ivertex)
          do kE = 1, nelmnts_by_node(ivertex_local)  !! loop on all ivertex connected elements
             New_element = elmnts_by_node(kE,ivertex_local)
             call store_new_element(New_element, nb_element_stored, stored_elements)
          enddo
       enddo
       nb_neigh(iE_loc) = nb_element_stored  !! number of neighbour of iE_loc element
    enddo
    size_adjacency = sum(nb_neigh(:))

    allocate(adjcy(size_adjacency),id_adjcy(0:nE_loc))
    id_adjcy(0) = 0
    do iE_loc=1, nE_loc !! loop on all element in partition
       iE = loc2glob_elmnt(iE_loc)
       nb_element_stored = 0
       stored_elements(:)=0
       do inode = 1, NGNOD_EIGHT_CORNERS  !! loop only on the corner of the element
          ivertex = elmnts(inode,iE_loc)  !!
          ivertex_local = glob2loc_nodes(ivertex)
          do kE = 1, nelmnts_by_node(ivertex_local)  !! loop on all ivertex connected elements
             New_element = elmnts_by_node(kE,ivertex_local)
             call store_new_element(New_element, nb_element_stored, stored_elements)
          enddo
       enddo
       id_adjcy(iE_loc) = id_adjcy(iE_loc-1) + nb_element_stored  !! number of neighbour of iE_loc element
       adjcy(id_adjcy(iE_loc-1)+1:id_adjcy(iE_loc))=stored_elements(1:nb_element_stored)
    enddo

    ! deallocate temporary arrays to save memmory
    deallocate(nelmnts_by_node, stored_elements, nb_neigh, elmnts_by_node)

    if (myrank == 0) write(27,*) ' END OF ADJACY TABLE COMPUTATION'

  end subroutine compute_adjcy_table

!----------------------------------
!
! store new item in array
!
!----------------------------------

  subroutine store_new_element(new_indx, nb, indx)
    implicit none
    integer,               intent(in)    :: new_indx
    integer,               intent(inout) :: nb
    integer, dimension(:), intent(inout) :: indx
    integer                              :: i
    do i=1,nb
       if (new_indx == indx(i)) return
    enddo
    nb = nb + 1
    indx(nb) = new_indx

  end subroutine store_new_element


end module module_database
