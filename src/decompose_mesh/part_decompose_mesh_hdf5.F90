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

module part_decompose_mesh_hdf5

  use phdf5_utils

  use constants, only: MAX_STRING_LEN,NGNOD2D_FOUR_CORNERS,NGNOD_EIGHT_CORNERS

  implicit none


contains
  !--------------------------------------------------
  ! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_glob2loc_nodes_database_h5(gname_proc, h5, iproc, npgeo, &
                                           nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                           glob2loc_nodes, nnodes, num_phase)

  use constants, only: NDIM

  implicit none

  ! phf5 io object
  type(h5io) :: h5

  character(len=*), intent(in)     :: gname_proc
  integer, intent(in)     :: nnodes, iproc, num_phase
  integer, intent(inout)  :: npgeo

  double precision, dimension(NDIM,nnodes)  :: nodes_coords
  integer, dimension(:), pointer            :: glob2loc_nodes_nparts
  integer, dimension(:), pointer            :: glob2loc_nodes_parts
  integer, dimension(:), pointer            :: glob2loc_nodes

  integer  :: i, j, count

  ! for attribute (storing nnodes npgeo)
  character(len=10)              :: aname = "nnodes_loc"

  ! for glob2loc_nodes_this_proc
  integer, dimension(npgeo)            :: glob2loc_nodes_this_proc
  character(len=40)                    :: gdsetname
  ! for nodes_coords_this_proc
  double precision, dimension(3,npgeo) :: nodes_coords_this_proc
  character(len=40)                    :: ndsetname


  if (num_phase == 1) then
  ! counts number of points in partition
    npgeo = 0
    do i = 0, nnodes-1
      do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
        if (glob2loc_nodes_parts(j) == iproc) then
          npgeo = npgeo + 1
        endif
      enddo
    enddo

   else
    ! prepare the temporal arrays for hdf5 write.
    count = 1
    do i = 0, nnodes-1
      do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
        if (glob2loc_nodes_parts(j) == iproc) then
          ! store values in temporal arrays
          glob2loc_nodes_this_proc(count) = glob2loc_nodes(j)+1
          nodes_coords_this_proc(1,count) = nodes_coords(1,i+1)
          nodes_coords_this_proc(2,count) = nodes_coords(2,i+1)
          nodes_coords_this_proc(3,count) = nodes_coords(3,i+1)

          count = count +1
        endif
      enddo
    enddo

    ! group names
    gdsetname = "glob2loc_nodes"
    ndsetname = "nodes_coords"
 
    ! open the target group
    call h5_open_group(h5, gname_proc)

    ! create a dataset for partial glob2loc_nodes array
    call h5_write_dataset_1d_i(h5, gdsetname, glob2loc_nodes_this_proc)
    ! add attributes before closing this dataset
    ! create an attribute for npgeo
    call h5_add_attribute_i(h5,aname, (/npgeo/))
    call h5_close_dataset(h5)

    ! create a dataset for partial nnodes_coords
    call h5_write_dataset_2d_d(h5, ndsetname, nodes_coords_this_proc)
    call h5_close_dataset(h5)

  endif

  end subroutine write_glob2loc_nodes_database_h5


  !--------------------------------------------------
  ! Write material properties in the Database
  !--------------------------------------------------
  subroutine write_material_props_database_h5(gname_material,h5,count_def_mat, &
                                           count_undef_mat, mat_prop, undef_mat_prop)

  implicit none
  ! phf5 io object
  type(h5io) :: h5
  character(len=*), intent(in) :: gname_material
  integer, intent(in)  :: count_def_mat,count_undef_mat ! stored as attributions
  double precision, dimension(17,count_def_mat)  :: mat_prop ! stored as a dataset
  character(len=MAX_STRING_LEN), dimension(6,count_undef_mat) :: undef_mat_prop ! stored as a

  ! for attribute count_def_mat and count_undef_mat
  character(len=13)              :: m_aname = "count_def_mat"
  character(len=15)              :: u_aname = "count_undef_mat"
 
  ! for dataset mat_prop, undef_mat_prop
  character(len=40)              :: mdsetname = "mat_prop"
  character(len=40)              :: udsetname = "undef_mat_prop"

  ! prepare hdf5 write
  ! open the target group
  call h5_open_group(h5, gname_material)

  ! create a dataset for mat_prop
  ! database material definition
  !
  ! format:  #rho  #vp  #vs  #q_kappa  #q_mu  #anisotropy_flag  #domain_id  #q_kappa   for (visco)elastic and acoustic
  !
  ! format:  #rhos,#rhof,#phi,#tort,#eta,0,#domain_id,#kxx,#kxy,#kxz,#kyy,#kyz,#kzz,
  !          #kappas,#kappaf,#kappafr,#eta,#mufr for poroelastic
  !
  ! (note that this order of the properties is different than the input in nummaterial_velocity_file)
  !
  call h5_write_dataset_2d_d(h5, mdsetname, mat_prop)
  ! create an attribute for count_def_mat
  call h5_add_attribute_i(h5, m_aname, (/count_def_mat/))
  call h5_close_dataset(h5)

  ! create a dataset for undef_mat_prop
  call h5_write_dataset_2d_c(h5, udsetname, undef_mat_prop)
  call h5_add_attribute_i(h5, u_aname, (/count_undef_mat/))
  call h5_close_dataset(h5)

  call h5_close_group(h5)
  end subroutine  write_material_props_database_h5


  !--------------------------------------------------
  ! Write elements on boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_boundaries_database_h5(gname_proc, h5, iproc, nspec, nspec2D_xmin, nspec2D_xmax, &
                                       nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                       ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                       ibelm_ymax, ibelm_bottom, ibelm_top, &
                                       nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                       nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                       glob2loc_elmnts, glob2loc_nodes_nparts, &
                                       glob2loc_nodes_parts, glob2loc_nodes, part, NGNOD2D)

  implicit none
  ! phf5 io object
  type(h5io) :: h5

  character(len=*), intent(in)  :: gname_proc
  integer, intent(in)  :: iproc
  integer, intent(in)  :: nspec
  integer, intent(in)  :: NGNOD2D
  integer, intent(in)  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, &
    nspec2D_ymax, nspec2D_bottom, nspec2D_top

  integer, dimension(nspec2D_xmin),   intent(in) :: ibelm_xmin
  integer, dimension(nspec2D_xmax),   intent(in) :: ibelm_xmax
  integer, dimension(nspec2D_ymin),   intent(in) :: ibelm_ymin
  integer, dimension(nspec2D_ymax),   intent(in) :: ibelm_ymax
  integer, dimension(nspec2D_bottom), intent(in) :: ibelm_bottom
  integer, dimension(nspec2D_top),    intent(in) :: ibelm_top

  integer, dimension(NGNOD2D,nspec2D_xmin),   intent(in) :: nodes_ibelm_xmin
  integer, dimension(NGNOD2D,nspec2D_xmax),   intent(in) :: nodes_ibelm_xmax
  integer, dimension(NGNOD2D,nspec2D_ymin),   intent(in) :: nodes_ibelm_ymin
  integer, dimension(NGNOD2D,nspec2D_ymax),   intent(in) :: nodes_ibelm_ymax
  integer, dimension(NGNOD2D,nspec2D_bottom), intent(in) :: nodes_ibelm_bottom
  integer, dimension(NGNOD2D,nspec2D_top),    intent(in) :: nodes_ibelm_top
  integer, dimension(:), pointer :: glob2loc_elmnts
  integer, dimension(:), pointer :: glob2loc_nodes_nparts
  integer, dimension(:), pointer :: glob2loc_nodes_parts
  integer, dimension(:), pointer :: glob2loc_nodes
  integer, dimension(1:nspec)    :: part

  ! local parameters
  integer :: i,j,inode,ier,count
  integer, dimension(NGNOD2D) :: loc_node
  integer :: loc_nspec2D_xmin,loc_nspec2D_xmax,loc_nspec2D_ymin, &
             loc_nspec2D_ymax,loc_nspec2D_bottom,loc_nspec2D_top

  ! for n_elms_on_bound
  integer, dimension(6)          :: n_elms_on_bound
  character(len=40)              :: ndsetname
  ! for glob2loc_elms
  integer, dimension(:,:), allocatable :: glob2loc_elms_this_proc
  character(len=40)                    :: gdsetname
  integer, dimension(2)                :: gdims

  ! counts number of elements for boundary at xmin, xmax, ymin, ymax, bottom, top in this partition
  loc_nspec2D_xmin = 0
  do i=1,nspec2D_xmin
     if (part(ibelm_xmin(i)) == iproc) then
        loc_nspec2D_xmin = loc_nspec2D_xmin + 1
     endif
  enddo
  !write(IIN_database) 1, loc_nspec2D_xmin
  n_elms_on_bound(1) = loc_nspec2D_xmin

  loc_nspec2D_xmax = 0
  do i=1,nspec2D_xmax
     if (part(ibelm_xmax(i)) == iproc) then
        loc_nspec2D_xmax = loc_nspec2D_xmax + 1
     endif
  enddo
  !write(IIN_database) 2, loc_nspec2D_xmax
  n_elms_on_bound(2) = loc_nspec2D_xmax

  loc_nspec2D_ymin = 0
  do i=1,nspec2D_ymin
     if (part(ibelm_ymin(i)) == iproc) then
        loc_nspec2D_ymin = loc_nspec2D_ymin + 1
     endif
  enddo
  !write(IIN_database) 3, loc_nspec2D_ymin
  n_elms_on_bound(3) = loc_nspec2D_ymin

  loc_nspec2D_ymax = 0
  do i=1,nspec2D_ymax
     if (part(ibelm_ymax(i)) == iproc) then
        loc_nspec2D_ymax = loc_nspec2D_ymax + 1
     endif
  enddo
  !write(IIN_database) 4, loc_nspec2D_ymax
  n_elms_on_bound(4) = loc_nspec2D_ymax

  loc_nspec2D_bottom = 0
  do i=1,nspec2D_bottom
     if (part(ibelm_bottom(i)) == iproc) then
        loc_nspec2D_bottom = loc_nspec2D_bottom + 1
     endif
  enddo
  !write(IIN_database) 5, loc_nspec2D_bottom
  n_elms_on_bound(5) = loc_nspec2D_bottom

  loc_nspec2D_top = 0
  do i=1,nspec2D_top
     if (part(ibelm_top(i)) == iproc) then
        loc_nspec2D_top = loc_nspec2D_top + 1
     endif
  enddo
  !write(IIN_database) 6, loc_nspec2D_top
  n_elms_on_bound(6) = loc_nspec2D_top


  ! allocate the array size of glob2loc_elms_this_proc here using the element counts
  gdims = (/1+NGNOD2D,loc_nspec2d_xmin  +loc_nspec2d_xmax &
                     +loc_nspec2d_ymin  +loc_nspec2d_ymax &
                     +loc_nspec2D_bottom+loc_nspec2d_top/)
  allocate(glob2loc_elms_this_proc(gdims(1),gdims(2)),stat=ier)
  if (ier /= 0) stop 'Error allocating array glob2loc_elms_this_proc'

  count = 1 ! initialilze the counter for glob2loc_elms_this_proc

  ! outputs element index and element node indices
  ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
  !          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
  !          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus
  !          we need to have the arg of glob2loc_elmnts start at 0, and thus we use glob2loc_nodes(ibelm_** -1)
  do i=1,nspec2D_xmin
     if (part(ibelm_xmin(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmin(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmin(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_xmin(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_xmin(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
     endif
  enddo

  do i=1,nspec2D_xmax
     if (part(ibelm_xmax(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_xmax(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_xmax(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_xmax(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_xmax(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
      endif
  enddo

  do i=1,nspec2D_ymin
     if (part(ibelm_ymin(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymin(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymin(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_ymin(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_ymin(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
      endif
  enddo

  do i=1,nspec2D_ymax
     if (part(ibelm_ymax(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_ymax(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_ymax(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_ymax(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_ymax(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
     endif
  enddo

  do i=1,nspec2D_bottom
     if (part(ibelm_bottom(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_bottom(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_bottom(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_bottom(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_bottom(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
     endif
  enddo

  do i=1,nspec2D_top
     if (part(ibelm_top(i)) == iproc) then
        do inode = 1,NGNOD2D
          do j = glob2loc_nodes_nparts(nodes_ibelm_top(inode,i)-1), &
                  glob2loc_nodes_nparts(nodes_ibelm_top(inode,i))-1
             if (glob2loc_nodes_parts(j) == iproc) then
                loc_node(inode) = glob2loc_nodes(j)+1
             endif
          enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_top(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        glob2loc_elms_this_proc(1,count) = glob2loc_elmnts(ibelm_top(i)-1)+1
        do inode = 1, NGNOD2D
          glob2loc_elms_this_proc(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
     endif
  enddo

  ! group names
  ndsetname = "n_elms_on_bound"
  gdsetname = "glob2loc_elms"

  ! prepare hdf5 write
  call h5_open_group(h5, gname_proc)

  ! write n_elms_on_bound
  call h5_write_dataset_1d_i(h5, ndsetname, n_elms_on_bound)
  call h5_close_dataset(h5)

  ! write glob2loc_elms_this_proc
  call h5_write_dataset_2d_i(h5, gdsetname, glob2loc_elms_this_proc)
  call h5_close_dataset(h5)

  ! close hdf5
  call h5_close_group(h5)

  ! deallocate arrays
  deallocate(glob2loc_elms_this_proc,stat=ier); if (ier /= 0) stop 'Error deallocating array glob2loc_elms_this_proc'

  end subroutine write_boundaries_database_h5

  !--------------------------------------------------
  ! Write C-PML elements indices, CPML-regions and thickness of C-PML layer
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_cpml_database_h5(gname_proc, h5, iproc, nspec, nspec_cpml, CPML_to_spec, &
                                 CPML_regions, is_CPML, glob2loc_elmnts, part)

  implicit none
  ! phf5 io object
  type(h5io) :: h5

  character(len=*), intent(in)  :: gname_proc
  integer, intent(in)  :: iproc
  integer, intent(in)  :: nspec
  integer, intent(in)  :: nspec_cpml
  integer              :: nspec_local

  integer, dimension(nspec_cpml), intent(in) :: CPML_to_spec
  integer, dimension(nspec_cpml), intent(in) :: CPML_regions

  logical, dimension(nspec), intent(in) :: is_CPML

  integer, dimension(:), pointer :: glob2loc_elmnts

  integer, dimension(1:nspec), intent(in) :: part

  ! local parameters
  integer :: i,nspec_cpml_local,count,count2,ier

  ! for attribute nspec_cpml, nspec_cpml_local
  character(len=18)              :: ndsetname = "nspec_cpml_globloc"
  integer, dimension(2)          :: ncpmls

  ! for dataset elements_cpml
  integer, dimension(:,:), allocatable  :: elements_cpml
  character(len=40)                     :: edsetname

  ! for dataset if_cpml
  integer, dimension(:), allocatable  :: if_cpml
  character(len=40)              :: idsetname

  ! group names
  edsetname = "elements_cpml"
  idsetname = "if_cpml"

  ! prepare hdf5 write
  call h5_open_group(h5, gname_proc)

  if (nspec_cpml > 0) then
    ! dump number of C-PML elements in this partition
    nspec_cpml_local = 0
    do i=1,nspec_cpml
       if (part(CPML_to_spec(i)) == iproc) then
          nspec_cpml_local = nspec_cpml_local + 1
       endif
    enddo

    ! dump attributevalues
    ncpmls(1) = nspec_cpml
    ncpmls(2) = nspec_cpml_local

    count  = 1 ! initialize counter for elements_cpml array
    count2 = 1 ! initialize conuter for if_cpml

    allocate(elements_cpml(2,nspec_cpml_local),stat=ier)
    if (ier /= 0) stop 'Error allocating array elements_cpml'

    ! dump C-PML regions and C-PML spectral elements global indexing
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
       if (part(CPML_to_spec(i)) == iproc) then
         elements_cpml(1,count) = glob2loc_elmnts(CPML_to_spec(i)-1)+1
         elements_cpml(2,count) = CPML_regions(i)
         count = count + 1
       endif
    enddo

    ! dump mask of C-PML elements for all elements in this partition
    ! count number of element in this iproc
    do i=1,nspec
       if (part(i) == iproc) then
          !write(IIN_database) is_CPML(i)
          count2 = count2 + 1
       endif
    enddo

    nspec_local = count2
    allocate(if_cpml(nspec_local),stat=ier)
    if (ier /= 0) stop 'Error allocating array if_cpml'

    count2 = 1 ! reinitialize counter2 to reuse it below
    do i=1,nspec
       if (part(i) == iproc) then
         if (is_CPML(i)) then
            if_cpml(count2) = 1
         else
            if_cpml(count2) = 0
         endif
         count2 = count2 + 1
       endif
    enddo
    
    ! create a dataset for elements_cpml
    call h5_write_dataset_2d_i(h5, edsetname, elements_cpml)
    call h5_close_dataset(h5)

    ! create a dataset for if_cpml
    ! strangely here H5T_NATIVE_HBOOL may not be used.
    call h5_write_dataset_1d_i(h5, idsetname, if_cpml)
    call h5_close_dataset(h5)

    ! deallocate local arrays
    deallocate(elements_cpml,stat=ier); if (ier /= 0) stop 'Error deallocating array elements_cpml'
    deallocate(if_cpml,stat=ier);       if (ier /= 0) stop 'Error deallocating array if_cpml'

  else ! dummy dataset for no cpml case
    ncpmls(1) = 0
    ncpmls(2) = 0
   
  endif

  ! create a dataset for nspec_cpml and nspec_cpml_local
  call h5_write_dataset_1d_i(h5, ndsetname, ncpmls)
  call h5_close_dataset(h5)
 
  ! close hdf5
  call h5_close_group(h5)

  end subroutine write_cpml_database_h5


  !--------------------------------------------------
  ! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_partition_database_h5(gname_proc, h5, iproc, nspec_local, nspec, elmnts, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, &
                                      part, num_modele, NGNOD, num_phase)

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none
  ! phf5 io object
  type(h5io) :: h5

  character(len=*), intent(in)  :: gname_proc
  integer, intent(in)           :: iproc
  integer, intent(inout)        :: nspec_local

  integer, intent(in)                  :: nspec
  integer, intent(in)                  :: NGNOD
  integer, dimension(0:NGNOD*nspec-1)  :: elmnts

  integer, dimension(:), pointer :: glob2loc_elmnts
  integer, dimension(:), pointer :: glob2loc_nodes_nparts
  integer, dimension(:), pointer :: glob2loc_nodes_parts
  integer, dimension(:), pointer :: glob2loc_nodes

  integer, dimension(0:nspec-1)  :: part
  integer, dimension(2,nspec)    :: num_modele

  integer, intent(in)  :: num_phase

  ! local parameters
  integer  :: i,j,k,count
  integer, dimension(0:NGNOD-1)  :: loc_nodes

  ! for attribute (storing nnodes npgeo)
  character(len=11)              :: aname = "nspec_local"
  ! for elm_conn
  !integer, dimension(NGNOD, nspec_local) :: elm_conn
  integer, dimension(:, :), allocatable :: elm_conn
  character(len=40)              :: dsetname = "elm_conn"
  ! for elem_conn_xdmf (mesh visualization)
  integer, dimension(:, :), allocatable :: elm_conn_xdmf
  character(len=40)              :: dxsetname = "elm_conn_xdmf"
  ! for mat_mesh
  integer, dimension(:,:), allocatable :: mat_mesh
  character(len=40)              :: mdsetname = "mat_mesh"
  ! for ispec_local
  integer, dimension(:), allocatable :: ispec_local
  character(len=40)              :: idsetname = "ispec_local"

  if (num_phase == 1) then
     ! counts number of spectral elements in this partition
     nspec_local = 0
     do i = 0, nspec-1
        if (part(i) == iproc) then
           nspec_local = nspec_local + 1
        endif
     enddo

  else
    allocate(elm_conn(NGNOD, nspec_local))
    allocate(elm_conn_xdmf(1+NGNOD, nspec_local)) ! 1 additional col for cell type id
    allocate(mat_mesh(2, nspec_local))
    allocate(ispec_local(nspec_local))
  
    ! prepare a temporal array to be recorded in hdf5
    ! element corner indices
    count = 1
    do i = 0, nspec-1
      if (part(i) == iproc) then

        do j = 0, NGNOD-1
          do k = glob2loc_nodes_nparts(elmnts(i*NGNOD+j)), glob2loc_nodes_nparts(elmnts(i*NGNOD+j)+1)-1

            if (glob2loc_nodes_parts(k) == iproc) then
              loc_nodes(j) = glob2loc_nodes(k)
            endif
          enddo
        enddo

         ispec_local(count) = glob2loc_elmnts(i)+1
         mat_mesh(1,count) = num_modele(1,i+1)
         mat_mesh(2,count) = num_modele(2,i+1)
         do k = 0, NGNOD-1
           elm_conn(k+1,count) = loc_nodes(k)+1
         enddo
         ! elm_conn_xdmf
         if (NGNOD == 8) then
            elm_conn_xdmf(1,count) = 9 ! cell type id xdmf
            elm_conn_xdmf(2,count) = loc_nodes(0)!+1 elm id starts 0
            elm_conn_xdmf(3,count) = loc_nodes(1)
            elm_conn_xdmf(4,count) = loc_nodes(2)
            elm_conn_xdmf(5,count) = loc_nodes(3)
            elm_conn_xdmf(6,count) = loc_nodes(4)
            elm_conn_xdmf(7,count) = loc_nodes(5)
            elm_conn_xdmf(8,count) = loc_nodes(6)
            elm_conn_xdmf(9,count) = loc_nodes(7)
         else ! NGNOD = 27
            elm_conn_xdmf(1,count)  = 50 ! cell type id xdmf
            elm_conn_xdmf(2,count)  = loc_nodes(0)
            elm_conn_xdmf(3,count)  = loc_nodes(1)
            elm_conn_xdmf(4,count)  = loc_nodes(2)
            elm_conn_xdmf(5,count)  = loc_nodes(3)
            elm_conn_xdmf(6,count)  = loc_nodes(4)
            elm_conn_xdmf(7,count)  = loc_nodes(5)
            elm_conn_xdmf(8,count)  = loc_nodes(6)
            elm_conn_xdmf(9,count)  = loc_nodes(7)
            elm_conn_xdmf(10,count) = loc_nodes(8)
            elm_conn_xdmf(11,count) = loc_nodes(9)
            elm_conn_xdmf(12,count) = loc_nodes(10)
            elm_conn_xdmf(13,count) = loc_nodes(11)
            elm_conn_xdmf(14,count) = loc_nodes(12)
            elm_conn_xdmf(15,count) = loc_nodes(13)
            elm_conn_xdmf(16,count) = loc_nodes(14)
            elm_conn_xdmf(17,count) = loc_nodes(15)
            elm_conn_xdmf(18,count) = loc_nodes(16)
            elm_conn_xdmf(19,count) = loc_nodes(17)
            elm_conn_xdmf(20,count) = loc_nodes(18)
            elm_conn_xdmf(21,count) = loc_nodes(19)
            elm_conn_xdmf(22,count) = loc_nodes(26)
            elm_conn_xdmf(23,count) = loc_nodes(20)
            elm_conn_xdmf(24,count) = loc_nodes(25)
            elm_conn_xdmf(25,count) = loc_nodes(24)
            elm_conn_xdmf(26,count) = loc_nodes(22)
            elm_conn_xdmf(27,count) = loc_nodes(21)
            elm_conn_xdmf(28,count) = loc_nodes(23)
         endif

         count = count + 1
         ! writes out to file Numglob2loc_elmn.txt
         if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) write(124,*) i+1,glob2loc_elmnts(i)+1,iproc
      endif
    enddo

    ! prepare hdf5 write
    call h5_open_group(h5, gname_proc)
 
    ! create a dataset for elm_conn
    call h5_write_dataset_2d_i(h5, dsetname, elm_conn)
 
    ! create an attribute for npgeo
    call h5_add_attribute_i(h5, aname, (/nspec_local/))
    call h5_close_dataset(h5)

    ! create a dataset for elm_conn_xdmf
    call h5_write_dataset_2d_i(h5, dxsetname, elm_conn_xdmf)
    call h5_close_dataset(h5)

    ! write mat_mesh
    call h5_write_dataset_2d_i(h5, mdsetname, mat_mesh)
    call h5_close_dataset(h5)

    ! write ispec_local
    call h5_write_dataset_1d_i(h5, idsetname, ispec_local)
    call h5_close_dataset(h5)

    ! close group
    call h5_close_group(h5)

    deallocate(elm_conn,elm_conn_xdmf,mat_mesh,ispec_local)

  endif

  end subroutine write_partition_database_h5



  !--------------------------------------------------
  ! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_interfaces_database_h5(gname_proc, h5, tab_interfaces, tab_size_interfaces, &
                              iproc, ninterfaces, &
                              my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, &
                              glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                              glob2loc_nodes, num_phase, nparts)

  implicit none
  type(h5io) :: h5

  character(len=*), intent(in)  :: gname_proc
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

  integer  :: i, j, k, l, count_nb, count_elm, ier
  integer  :: num_interface

  integer  :: count_faces

  ! variables for my_ninterfaces_and_maxval
  integer, dimension(2)          :: num_interface_and_max
  character(len=40)              :: ndsetname 
 
	! my_nb_interfaces in subgroup num_interface (2,my_ninterfaces)
  integer, dimension(:,:), allocatable :: num_neighbors_elmnts
  character(len=40)                  :: mdsetname
 
  ! my_interfaces
  integer, dimension(:,:), allocatable :: neighbors_elmnts
  character(len=40)                    :: idsetname 
 
  num_interface = 0

  if (num_phase == 1) then
    ! counts number of interfaces to neighboring partitions
     my_interfaces(:) = 0
     my_nb_interfaces(:) = 0

     ! double loops over all partitions
     do i = 0, nparts-1
        do j = i+1, nparts-1
           ! only counts if specified partition (iproc) appears and interface elements increment
           if ((tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
                (i == iproc .or. j == iproc)) then
              ! sets flag
              my_interfaces(num_interface) = 1
              ! sets number of elements on interface
              my_nb_interfaces(num_interface) = tab_size_interfaces(num_interface+1) &
                                          - tab_size_interfaces(num_interface)
           endif
           num_interface = num_interface + 1
        enddo
     enddo
     my_ninterface = sum(my_interfaces(:))

  else
    count_nb  = 1 ! count for num_neighbors_elmnts
    count_elm = 1 ! count for neighbors_elmnts

    ! allocate tempral array size
    allocate(num_neighbors_elmnts(2,my_ninterface),stat=ier)
    if (ier /= 0) stop 'Error allocating array num_neighbors_elmnts'
    allocate(neighbors_elmnts(6, sum(my_nb_interfaces(:))))
    if (ier /= 0) stop 'Error allocating array neighbors_elmnts'

    ! writes out MPI interface elements
    do i = 0, nparts-1
      do j = i+1, nparts-1
        if (my_interfaces(num_interface) == 1) then
          if (i == iproc) then
            !write(IIN_database) j, my_nb_interfaces(num_interface)
            num_neighbors_elmnts(1,count_nb) = j
          else
            !write(IIN_database) i, my_nb_interfaces(num_interface)
            num_neighbors_elmnts(1,count_nb) = i
          endif
          num_neighbors_elmnts(2,count_nb) = my_nb_interfaces(num_interface)
          count_nb = count_nb + 1

          count_faces = 0
          do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
             if (i == iproc) then
                local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+0))+1
             else
                local_elmnt = glob2loc_elmnts(tab_interfaces(k*7+1))+1
             endif

             select case (tab_interfaces(k*7+2))
             case (1)
                ! single point element
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(1) = glob2loc_nodes(l)+1
                   endif
                enddo
                !write(IIN_database) local_elmnt, tab_interfaces(k*7+2), &
                !                   local_nodes(1), -1, -1, -1
                neighbors_elmnts(1,count_elm) = local_elmnt
                neighbors_elmnts(2,count_elm) = tab_interfaces(k*7+2)
                neighbors_elmnts(3,count_elm) = local_nodes(1) 
                neighbors_elmnts(4,count_elm) = -1
                neighbors_elmnts(5,count_elm) = -1
                neighbors_elmnts(6,count_elm) = -1

             case (2)
                ! edge element
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(1) = glob2loc_nodes(l)+1
                   endif
                enddo
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+4)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+4)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(2) = glob2loc_nodes(l)+1
                   endif
                enddo
                !write(IIN_database) local_elmnt, tab_interfaces(k*7+2), &
                !                   local_nodes(1), local_nodes(2), -1, -1
                neighbors_elmnts(1,count_elm) = local_elmnt
                neighbors_elmnts(2,count_elm) = tab_interfaces(k*7+2)
                neighbors_elmnts(3,count_elm) = local_nodes(1) 
                neighbors_elmnts(4,count_elm) = local_nodes(2) 
                neighbors_elmnts(5,count_elm) = -1
                neighbors_elmnts(6,count_elm) = -1

             case (4)
                ! face element
                count_faces = count_faces + 1
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+3)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+3)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(1) = glob2loc_nodes(l)+1
                   endif
                enddo
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+4)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+4)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(2) = glob2loc_nodes(l)+1
                   endif
                enddo
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+5)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+5)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(3) = glob2loc_nodes(l)+1
                   endif
                enddo
                do l = glob2loc_nodes_nparts(tab_interfaces(k*7+6)), &
                     glob2loc_nodes_nparts(tab_interfaces(k*7+6)+1)-1
                   if (glob2loc_nodes_parts(l) == iproc) then
                      local_nodes(4) = glob2loc_nodes(l)+1
                   endif
                enddo
                !write(IIN_database) local_elmnt, tab_interfaces(k*7+2), &
                !     local_nodes(1), local_nodes(2),local_nodes(3), local_nodes(4)
                neighbors_elmnts(1,count_elm) = local_elmnt
                neighbors_elmnts(2,count_elm) = tab_interfaces(k*7+2)
                neighbors_elmnts(3,count_elm) = local_nodes(1) 
                neighbors_elmnts(4,count_elm) = local_nodes(2) 
                neighbors_elmnts(5,count_elm) = local_nodes(3) 
                neighbors_elmnts(6,count_elm) = local_nodes(4) 

             case default
                print *, "fatal error in write_interfaces_database:", tab_interfaces(k*7+2), iproc
                stop "fatal error in write_interfaces_database"
             end select
          
          count_elm = count_elm + 1
          enddo

        endif

         num_interface = num_interface + 1
      enddo
    enddo

    ! hdf5 write
    call h5_open_group(h5, gname_proc)
 
    ! group names
    ndsetname = "my_ninterface_and_max"
    mdsetname = "my_nb_interfaces"
    idsetname = "my_interfaces"

    ! write my_ninterface and maxval
    if (my_ninterface == 0) then
      num_interface_and_max = (/my_ninterface, 0/)
    else
      num_interface_and_max = (/my_ninterface, maxval(my_nb_interfaces)/)
    endif
    call h5_write_dataset_1d_i(h5, ndsetname, num_interface_and_max)
    call h5_close_dataset(h5)

    ! write my_nb_interfaces
    call h5_write_dataset_2d_i(h5, mdsetname, num_neighbors_elmnts)
    call h5_close_dataset(h5)

    ! write my_interfaces
    call h5_write_dataset_2d_i(h5, idsetname, neighbors_elmnts)
    call h5_close_dataset(h5)
 
    ! close hdf5
    call h5_close_group(h5)

    ! allocate temporal array size
    deallocate(num_neighbors_elmnts,stat=ier); if (ier /= 0) stop 'Error deallocating array num_neighbors_elmnts'
    deallocate(neighbors_elmnts,stat=ier);     if (ier /= 0) stop 'Error deallocating array neighbors_elmnts'
 
  endif

  end subroutine write_interfaces_database_h5

  !--------------------------------------------------
  ! Write elements on surface boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_moho_surface_database_h5(gname_proc, h5, iproc, nspec, &
                                         glob2loc_elmnts, glob2loc_nodes_nparts, &
                                         glob2loc_nodes_parts, glob2loc_nodes, part, &
                                         nspec2D_moho,ibelm_moho,nodes_ibelm_moho, NGNOD2D)

  implicit none
  type(h5io) :: h5

  character(len=*), intent(in)  :: gname_proc
  integer, intent(in)           :: iproc
  integer, intent(in)           :: nspec
  integer, intent(in)           :: NGNOD2D

  integer, dimension(:), pointer :: glob2loc_elmnts
  integer, dimension(:), pointer :: glob2loc_nodes_nparts
  integer, dimension(:), pointer :: glob2loc_nodes_parts
  integer, dimension(:), pointer :: glob2loc_nodes
  integer, dimension(1:nspec)    :: part

  integer, intent(in)                                  :: nspec2D_moho
  integer, dimension(nspec2D_moho), intent(in)         :: ibelm_moho
  integer, dimension(NGNOD2D,nspec2D_moho), intent(in) :: nodes_ibelm_moho

  integer                     :: i,j,inode,count,ier
  integer, dimension(NGNOD2D) :: loc_node
  integer                     :: loc_nspec2D_moho

  ! for loc_nspec_2d_homo attribute
  character(len=8)               :: aname = "loc_moho"
  integer         , dimension(2) :: moho_attr

  ! homo_elements dataset
  integer, dimension(:,:), allocatable :: loc_moho_temp
  character(len=40)                    :: dsetname ! may be in "/group/dataset_name" or in subgroup "/group/subgroup/.../dataset_name"

  ! group names
  dsetname = "moho_elms"

  ! counts number of elements for moho surface in this partition
  ! optional moho
  loc_nspec2D_moho = 0
  do i=1,nspec2D_moho
     if (part(ibelm_moho(i)) == iproc) then
        loc_nspec2D_moho = loc_nspec2D_moho + 1
     endif
  enddo
  ! checks if anything to do
  if (loc_nspec2D_moho == 0) return

  ! format: #surface_id, #number of elements
  !write(IIN_database) 7, loc_nspec2D_moho
  moho_attr = (/7, loc_nspec2D_moho/)

  ! outputs element index and element node indices
  ! note: assumes that element indices in ibelm_* arrays are in the range from 1 to nspec
  !          (this is assigned by CUBIT, if this changes the following indexing must be changed as well)
  !          while glob2loc_elmnts(.) is shifted from 0 to nspec-1  thus
  !          we need to have the arg of glob2loc_elmnts start at 0, and thus we use glob2loc_nodes(ibelm_** -1)

  ! optional moho
  allocate(loc_moho_temp(1+NGNOD2D,loc_nspec2D_moho),stat=ier)
  if (ier /= 0) stop 'Error allocating array loc_moho_temp'
  count = 1

  do i=1,nspec2D_moho
     if (part(ibelm_moho(i)) == iproc) then
        do inode = 1,NGNOD2D
        do j = glob2loc_nodes_nparts(nodes_ibelm_moho(inode,i)-1), &
                glob2loc_nodes_nparts(nodes_ibelm_moho(inode,i))-1
           if (glob2loc_nodes_parts(j) == iproc) then
              loc_node(inode) = glob2loc_nodes(j)+1
           endif
        enddo
        enddo
        !write(IIN_database) glob2loc_elmnts(ibelm_moho(i)-1)+1, (loc_node(inode), inode = 1,NGNOD2D)
        loc_moho_temp(1,count) = glob2loc_elmnts(ibelm_moho(i)-1)+1
        do inode = 1, NGNOD2D
          loc_moho_temp(inode+1,count) = loc_node(inode)
        enddo
        count = count + 1
      endif
  enddo

  ! hdf5 file output
  call h5_open_group(h5, gname_proc)

  ! create datasets and write
  call h5_write_dataset_2d_i(h5, dsetname, loc_moho_temp)
  ! attribute
  call h5_add_attribute_i(h5, aname, moho_attr)
  call h5_close_dataset(h5)

  ! close hdf5 objects
  call h5_close_group(h5)

  deallocate(loc_moho_temp,stat=ier); if (ier /= 0) stop 'Error deallocating array loc_moho_temp'
 
  end subroutine write_moho_surface_database_h5


end module part_decompose_mesh_hdf5

