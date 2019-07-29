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

  use hdf5

  use constants, only: MAX_STRING_LEN,NGNOD2D_FOUR_CORNERS,NGNOD_EIGHT_CORNERS

  implicit none

! group names in hdf5
character(len=14), parameter :: gname_material   = "material_props" ! Group name for material properties

contains
  !--------------------------------------------------
  ! hd5 file initialization (create the file, groups(nodes etc) and sub-groups(proc_n)
  !--------------------------------------------------
  subroutine hdf5_initialize(IIN_database, nparts)
    implicit none
    character(len=*), intent(in)  :: IIN_database
    integer, intent(in) :: nparts

    integer(HID_T)  :: file_id, group_id
    integer         :: ipart, error
    character(len=5),  parameter :: gname_proc_head  = "proc_"          ! Group name for each processor
    character(len=11) :: gname_proc
    character(len=10) :: tempstr


    ! create file
    call h5fcreate_f(IIN_database, H5F_ACC_TRUNC_F, file_id, error)
    if (error /= 0) then
     print *,'Error file open:',IIN_database
     print *
     print *,'check if path exists:',IIN_database
     stop 'Error file open Database'
    endif

    ! create a group for material
    call h5gcreate_f(file_id, gname_material, group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 2'
    call h5gclose_f(group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 3'
 
    ! processor groups
    do ipart = 0, nparts-1
      write(tempstr, "(i6.6)") ipart
      gname_proc = gname_proc_head//trim(tempstr)
      call h5gcreate_f(file_id, gname_proc, group_id, error)
      if (error == 0) write(*,*) 'hdf5 passed at 4'
      call h5gclose_f(group_id, error)
      if (error == 0) write(*,*) 'hdf5 passed at 5'
   enddo

    ! close file
    call h5fclose_f(file_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 6'
 
  end subroutine hdf5_initialize


  !--------------------------------------------------
  ! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_glob2loc_nodes_database_h5(IIN_database, gname_proc, file_id, iproc, npgeo, &
                                           nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                                           glob2loc_nodes, nnodes, num_phase)

  use constants, only: NDIM

  implicit none
  character(len=*), intent(in)     :: IIN_database
  character(len=*), intent(in)     :: gname_proc
  integer, intent(in)     :: nnodes, iproc, num_phase
  integer, intent(inout)  :: npgeo

  double precision, dimension(NDIM,nnodes)  :: nodes_coords
  integer, dimension(:), pointer            :: glob2loc_nodes_nparts
  integer, dimension(:), pointer            :: glob2loc_nodes_parts
  integer, dimension(:), pointer            :: glob2loc_nodes

  integer  :: i, j, count, ier

  ! variables for hdf5 write
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error
  integer                    :: dim = 3 ! dimension

  ! for attribute (storing nnodes npgeo)
  integer(HID_T)                 :: aspace_id, atype_id, attr_id ! attribute id
  integer                        :: arank = 1 ! attribute rank
  integer(HSIZE_T), dimension(1) :: adims = (/2/) ! attribute dimension
  integer(HSIZE_T)               :: taglen = 10 ! loc_nnodes
  character(len=10)              :: aname = "nnodes_loc"
  integer(HSIZE_T), dimension(1) :: data_dims = 2
  ! for glob2loc_nodes_this_proc
  !integer, dimension(:), allocatable            :: glob2loc_nodes_this_proc
  integer, dimension(npgeo)           :: glob2loc_nodes_this_proc
  integer(HID_T)                 :: gdspace_id, gdset_id
  integer                        :: grank = 1
  integer(HSIZE_T), dimension(1) :: gdims
  character(len=40)              :: gdsetname
  ! for nodes_coords_this_proc
  double precision, dimension(3,npgeo) :: nodes_coords_this_proc
  integer(HID_T)                 :: ndspace_id, ndset_id
  integer                        :: nrank = 2
  integer(HSIZE_T), dimension(2) :: ndims
  character(len=40)              :: ndsetname
 
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
 
    ! prepare hdf5 write
    call h5gopen_f(file_id, gname_proc, group_id, error) ! group open
    if (error == 0) write(*,*) 'hdf5 passed at 10'
 
    ! create a dataset for partial glob2loc_nodes array
    gdims = (/npgeo/)
    call h5screate_simple_f(grank, gdims, gdspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 11'
    call h5dcreate_f(group_id, gdsetname, H5T_NATIVE_INTEGER, gdspace_id, gdset_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 12'
    call h5dwrite_f(gdset_id, H5T_NATIVE_INTEGER, glob2loc_nodes_this_proc, gdims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 13'
    ! create a dataset for partial nnodes_coords
    ndims = (/dim,npgeo/)
    call h5screate_simple_f(nrank, ndims, ndspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 14'
    call h5dcreate_f(group_id, ndsetname, H5T_NATIVE_DOUBLE, ndspace_id, ndset_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 15'
    call h5dwrite_f(ndset_id, H5T_NATIVE_DOUBLE, nodes_coords_this_proc, ndims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 16'
 
    ! create an attribute for npgeo
    call h5screate_simple_f(arank, adims, aspace_id, error) ! for a integernvalue
    if (error == 0) write(*,*) 'hdf5 passed at 17'
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! for a string tag of attribute value
    if (error == 0) write(*,*) 'hdf5 passed at 18'
    call h5tset_size_f(atype_id, taglen, error)
    if (error == 0) write(*,*) 'hdf5 passed at 19'
    call h5acreate_f(gdset_id, aname, atype_id, aspace_id, attr_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 20'
    call h5awrite_f(attr_id, atype_id, npgeo, data_dims, error) ! write
    if (error == 0) write(*,*) 'hdf5 passed at 21'
 
    ! close hdf5
    call h5sclose_f(gdspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 22'
    call h5sclose_f(ndspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 23'
    call h5sclose_f(aspace_id,  error)
    if (error == 0) write(*,*) 'hdf5 passed at 24'
    call h5dclose_f(gdset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 25'
    call h5dclose_f(ndset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 26'
    call h5aclose_f(attr_id,    error)
    if (error == 0) write(*,*) 'hdf5 passed at 27'
    call h5tclose_f(atype_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 28'
 
    call h5gclose_f(group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 29'

  endif

  end subroutine write_glob2loc_nodes_database_h5


  !--------------------------------------------------
  ! Write material properties in the Database
  !--------------------------------------------------
  subroutine write_material_props_database_h5(IIN_database,file_id,count_def_mat, &
                                           count_undef_mat, mat_prop, undef_mat_prop)

  implicit none
  character(len=*), intent(in)  :: IIN_database
  integer, intent(in)  :: count_def_mat,count_undef_mat ! stored as attributions
  double precision, dimension(17,count_def_mat)  :: mat_prop ! stored as a dataset
  character(len=MAX_STRING_LEN), dimension(6,count_undef_mat) :: undef_mat_prop ! stored as a
  integer  :: i

  ! variables for hdf5 write
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error

  ! for attribute count_def_mat and count_undef_mat
  integer(HID_T)                 :: m_aspace_id, m_atype_id, m_attr_id, &
                                    u_aspace_id, u_atype_id, u_attr_id  ! attribute id
  integer                        :: arank = 1 ! attribute rank
  integer(HSIZE_T), dimension(1) :: adims = (/2/) ! attribute dimension
  integer(HSIZE_T)               :: m_taglen = 13 ! mat_
  integer(HSIZE_T)               :: u_taglen = 15 ! mat_
  character(len=13)              :: m_aname = "count_def_mat"
  character(len=15)              :: u_aname = "count_undef_mat"
  integer(HSIZE_T), dimension(1) :: data_dims = 2
 
  ! for dataset mat_prop, undef_mat_prop
  integer(HID_T)                 :: m_dspace_id, m_dset_id, &
                                    u_dspace_id, u_dset_id
  integer                        :: mrank = 2
  integer(HSIZE_T), dimension(2) :: mdims
  character(len=40)              :: mdsetname = "mat_prop"
  integer                        :: urank = 2
  integer(HSIZE_T), dimension(2) :: udims
  character(len=40)              :: udsetname = "undef_mat_prop"

  ! prepare hdf5 write
  call h5gopen_f(file_id, gname_material, group_id, error) ! group open
  if (error == 0) write(*,*) 'hdf5 passed at 34'
 
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
  mdims = (/17,count_def_mat/)
  call h5screate_simple_f(mrank, mdims, m_dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 35'
  call h5dcreate_f(group_id, mdsetname, h5t_native_double, m_dspace_id, m_dset_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 36'
  call h5dwrite_f(m_dset_id, h5t_native_double, mat_prop, mdims, error)
  if (error == 0) write(*,*) 'hdf5 passed at 37'
  ! create a dataset for undef_mat_prop
  mdims = (/6,count_undef_mat/)
  call h5screate_simple_f(urank, udims, u_dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 38'
  call h5dcreate_f(group_id, udsetname, h5t_native_double, u_dspace_id, u_dset_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 39'
  call h5dwrite_f(u_dset_id, h5t_native_double, undef_mat_prop, udims, error)
  if (error == 0) write(*,*) 'hdf5 passed at 40'
 
  ! create an attribute for count_def_mat
  call h5screate_simple_f(arank, adims, m_aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 41'
  call h5tcopy_f(h5t_native_character, m_atype_id, error) ! for a string tag of attribute value
  if (error == 0) write(*,*) 'hdf5 passed at 42'
  call h5tset_size_f(m_atype_id, m_taglen, error)
  if (error == 0) write(*,*) 'hdf5 passed at 43'
  call h5acreate_f(m_dset_id, m_aname, m_atype_id, m_aspace_id, m_attr_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 44'
  call h5awrite_f(m_attr_id, m_atype_id, count_def_mat, data_dims, error) ! write data
  if (error == 0) write(*,*) 'hdf5 passed at 45'
  ! create an attribute for count_undef_mat
  call h5screate_simple_f(arank, adims, u_aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 46'
  call h5tcopy_f(h5t_native_character, u_atype_id, error) ! for a string tag of attribute value
  if (error == 0) write(*,*) 'hdf5 passed at 47'
  call h5tset_size_f(u_atype_id, u_taglen, error)
  if (error == 0) write(*,*) 'hdf5 passed at 48'
  call h5acreate_f(u_dset_id, u_aname, u_atype_id, u_aspace_id, u_attr_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 49'
  call h5awrite_f(u_attr_id, u_atype_id, count_undef_mat, data_dims, error) ! write data
  if (error == 0) write(*,*) 'hdf5 passed at 50'
 
  ! close hdf5
  call h5sclose_f(m_dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 51'
  call h5sclose_f(m_aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 52'
  call h5dclose_f(m_dset_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 53'
  call h5aclose_f(m_attr_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 54'
  call h5tclose_f(m_atype_id,  error)
  if (error == 0) write(*,*) 'hdf5 passed at 55'
  call h5sclose_f(u_dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 56'
  call h5sclose_f(u_aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 57'
  call h5dclose_f(u_dset_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 58'
  call h5aclose_f(u_attr_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 59'
  call h5tclose_f(u_atype_id,  error)
  if (error == 0) write(*,*) 'hdf5 passed at 60'
 
  call h5gclose_f(group_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 61'
 
  end subroutine  write_material_props_database_h5


  !--------------------------------------------------
  ! Write elements on boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_boundaries_database_h5(IIN_database, gname_proc, file_id, iproc, nspec, nspec2D_xmin, nspec2D_xmax, &
                                       nspec2D_ymin, nspec2D_ymax, nspec2D_bottom, nspec2D_top, &
                                       ibelm_xmin, ibelm_xmax, ibelm_ymin, &
                                       ibelm_ymax, ibelm_bottom, ibelm_top, &
                                       nodes_ibelm_xmin, nodes_ibelm_xmax, nodes_ibelm_ymin, &
                                       nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top, &
                                       glob2loc_elmnts, glob2loc_nodes_nparts, &
                                       glob2loc_nodes_parts, glob2loc_nodes, part, NGNOD2D)

  implicit none
  character(len=*), intent(in)  :: IIN_database
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

  ! variables for hdf5 write
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error

  ! for n_elms_on_bound
  integer, dimension(6)          :: n_elms_on_bound
  integer(HID_T)                 :: ndspace_id, ndset_id
  integer                        :: nrank = 1
  integer(HSIZE_T), dimension(1) :: ndims = (/6/)
  character(len=40)              :: ndsetname
  ! for glob2loc_elms
  integer, dimension(:,:), allocatable :: glob2loc_elms_this_proc
  integer(HID_T)                       :: gdspace_id, gdset_id
  integer                              :: grank = 2
  integer(HSIZE_T), dimension(2)       :: gdims
  character(len=40)                    :: gdsetname

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
  call h5gopen_f(file_id, gname_proc, group_id, error) ! group open
  if (error == 0) write(*,*) 'hdf5 passed at 66'
 
  ! write n_elms_on_bound
  call h5screate_simple_f(nrank, ndims, ndspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 67'
  call h5dcreate_f(group_id, ndsetname, H5T_NATIVE_INTEGER, ndspace_id, ndset_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 68'
  call h5dwrite_f(ndset_id, H5T_NATIVE_INTEGER, n_elms_on_bound, ndims, error)
  if (error == 0) write(*,*) 'hdf5 passed at 69'
 
  ! write glob2loc_elms_this_proc
  call h5screate_simple_f(grank, gdims, gdspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 70'
  call h5dcreate_f(group_id, gdsetname, H5T_NATIVE_INTEGER, gdspace_id, gdset_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 71'
  call h5dwrite_f(gdset_id, H5T_NATIVE_INTEGER, glob2loc_elms_this_proc, gdims, error)
  if (error == 0) write(*,*) 'hdf5 passed at 72'
 
  ! close hdf5
  call h5sclose_f(gdspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 73'
  call h5sclose_f(ndspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 74'
  call h5dclose_f(gdset_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 75'
  call h5dclose_f(ndset_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 76'
  
  call h5gclose_f(group_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 77'
 
  ! deallocate arrays
  deallocate(glob2loc_elms_this_proc,stat=ier); if (ier /= 0) stop 'Error deallocating array glob2loc_elms_this_proc'

  end subroutine write_boundaries_database_h5

  !--------------------------------------------------
  ! Write C-PML elements indices, CPML-regions and thickness of C-PML layer
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_cpml_database_h5(IIN_database, gname_proc, file_id, iproc, nspec, nspec_cpml, CPML_to_spec, &
                                 CPML_regions, is_CPML, glob2loc_elmnts, part)

  implicit none
  character(len=*), intent(in)  :: IIN_database
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

  ! hdf5 variables
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)  :: group_id
  integer         ::  error

  ! for attribute nspec_cpml, nspec_cpml_local
  integer(HID_T)                 :: aspace_id, atype_id, attr_id ! attribute id
  integer                        :: arank = 1 ! attribute rank
  integer(HSIZE_T), dimension(1) :: adims = (/2/) ! attribute dimension
  integer(HSIZE_T)               :: taglen = 20 ! nspec_cpml
  character(len=18)              :: aname = "nspec_cpml_globloc"
  integer(HSIZE_T), dimension(1) :: data_dims = 2
  integer, dimension(2)          :: attr_data
 
  ! for dataset elements_cpml
  integer, dimension(:,:), allocatable  :: elements_cpml
  integer(HID_T)                 :: edspace_id, edset_id
  integer                        :: edrank = 2
  integer(HSIZE_T), dimension(2) :: eddims
  character(len=40)              :: edsetname

  ! for dataset if_cpml
  !logical, dimension(:), allocatable  :: if_cpml
  integer, dimension(:), allocatable  :: if_cpml
  integer(HID_T)                 :: idspace_id, idset_id
  integer                        :: idrank = 1
  integer(HSIZE_T), dimension(1) :: iddims
  character(len=40)              :: idsetname

  ! group names
  edsetname = "elements_cpml"
  idsetname = "if_cpml"

  if (nspec_cpml > 0) then
    ! dump number of C-PML elements in this partition
    nspec_cpml_local = 0
    do i=1,nspec_cpml
       if (part(CPML_to_spec(i)) == iproc) then
          nspec_cpml_local = nspec_cpml_local + 1
       endif
    enddo

    ! dump attributevalues
    attr_data(1) = nspec_cpml
    attr_data(2) = nspec_cpml_local

    count  = 1 ! initialize counter for elements_cpml array
    count2 = 1 ! initialize conuter for if_cpml

    allocate(elements_cpml(2,nspec_cpml_local),stat=ier)
    if (ier /= 0) stop 'Error allocating array elements_cpml'
    eddims = (/2,nspec_cpml_local/)

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
    iddims = (/nspec_local/)
    count2 = 1 ! reinitialize counter2 to reuse it below
    do i=1,nspec
       if (part(i) == iproc) then
         if_cpml(count2) = is_CPML(i) 
         count2 = count2 + 1
       endif
    enddo
 
    ! prepare hdf5 write
    call h5gopen_f(file_id, gname_proc, group_id, error) ! group open
    if (error == 0) write(*,*) 'hdf5 passed at 82'
 
    ! create a dataset for elements_cpml
    call h5screate_simple_f(edrank, eddims, edspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 83'
    call h5dcreate_f(group_id, edsetname, H5T_NATIVE_INTEGER, edspace_id, edset_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 84'
    call h5dwrite_f(edset_id, H5T_NATIVE_INTEGER, elements_cpml, eddims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 85'
    
    ! create a dataset for if_cpml
    call h5screate_simple_f(idrank, iddims, idspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 86'
    ! strangely here H5T_NATIVE_HBOOL may not be used. Need to test with other HDF5 versions
    !call h5dcreate_f(group_id, idsetname, H5T_NATIVE_HBOOL, idspace_id, idset_id, error)
    !call h5dwrite_f(idset_id, H5T_NATIVE_HBOOL, if_cpml, iddims, error)
    call h5dcreate_f(group_id, idsetname, H5T_NATIVE_INTEGER, idspace_id, idset_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 87'
    call h5dwrite_f(idset_id, H5T_NATIVE_INTEGER, if_cpml, iddims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 88'
 
    ! create an attribute for nspec_cpml and nspec_cpml_local
    call h5screate_simple_f(arank, adims, aspace_id, error) ! for a integernvalue
    if (error == 0) write(*,*) 'hdf5 passed at 89'
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! for a string tag of attribute value
    if (error == 0) write(*,*) 'hdf5 passed at 90'
    call h5tset_size_f(atype_id, taglen, error)
    if (error == 0) write(*,*) 'hdf5 passed at 91'
    call h5acreate_f(edset_id, aname, atype_id, aspace_id, attr_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 92'
    call h5awrite_f(attr_id, atype_id, attr_data, data_dims, error) ! write
    if (error == 0) write(*,*) 'hdf5 passed at 93'
 
    ! close hdf5
    call h5sclose_f(edspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 94'
    call h5sclose_f(idspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 95'
    call h5sclose_f(aspace_id,  error)
    if (error == 0) write(*,*) 'hdf5 passed at 96'
    call h5dclose_f(edset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 97'
    call h5dclose_f(idset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 98'
    call h5aclose_f(attr_id,    error)
    if (error == 0) write(*,*) 'hdf5 passed at 99'
    call h5tclose_f(atype_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 100'
 
    call h5gclose_f(group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 101'
 
    ! deallocate local arrays
    deallocate(elements_cpml,stat=ier); if (ier /= 0) stop 'Error deallocating array elements_cpml'
    deallocate(if_cpml,stat=ier);       if (ier /= 0) stop 'Error deallocating array if_cpml'

  endif

  end subroutine write_cpml_database_h5


  !--------------------------------------------------
  ! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_partition_database_h5(IIN_database, gname_proc, file_id, iproc, nspec_local, nspec, elmnts, &
                                      glob2loc_elmnts, glob2loc_nodes_nparts, &
                                      glob2loc_nodes_parts, glob2loc_nodes, &
                                      part, num_modele, NGNOD, num_phase)

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,MESH_A_CHUNK_OF_THE_EARTH

  implicit none

  character(len=*), intent(in)  :: IIN_database
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
  integer  :: i,j,k,count,ier
  integer, dimension(0:NGNOD-1)  :: loc_nodes

  ! variables for hdf5 write
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error
  integer                    :: dim = 3 ! dimension

  ! for attribute (storing nnodes npgeo)
  integer(HID_T)                 :: aspace_id, atype_id, attr_id ! attribute id
  integer                        :: arank = 1 ! attribute rank
  integer(HSIZE_T), dimension(1) :: adims = (/2/) ! attribute dimension
  integer(HSIZE_T)               :: taglen = 10 ! loc_nnodes
  character(len=10)              :: aname = "nspec_local"
  integer(HSIZE_T), dimension(1) :: data_dims = 2
  ! for elm_conn
  integer, dimension(3+NGNOD, nspec_local) :: elm_conn
  integer(HID_T)                 :: dspace_id, dset_id
  integer                        :: drank = 2
  integer(HSIZE_T), dimension(2) :: ddims
  character(len=40)              :: dsetname

  
  if (num_phase == 1) then
     ! counts number of spectral elements in this partition
     nspec_local = 0
     do i = 0, nspec-1
        if (part(i) == iproc) then
           nspec_local = nspec_local + 1
        endif
     enddo

  else
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

         ! store all data to the temporal array
         ! format:
         ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
         ! or
         ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id27
         !write(IIN_database) glob2loc_elmnts(i)+1,num_modele(1,i+1),num_modele(2,i+1),(loc_nodes(k)+1, k=0,NGNOD-1)
         elm_conn(1,count) = glob2loc_elmnts(i)+1
         elm_conn(2,count) = num_modele(1,i+1)
         elm_conn(3,count) = num_modele(2,i+1)
         do k = 0, NGNOD-1
           elm_conn(4+k,count) = loc_nodes(k)+1
         enddo
         count = count + 1
         ! writes out to file Numglob2loc_elmn.txt
         if (COUPLE_WITH_INJECTION_TECHNIQUE .or. MESH_A_CHUNK_OF_THE_EARTH) write(124,*) i+1,glob2loc_elmnts(i)+1,iproc
      endif
    enddo

    ! group names
    dsetname = "elm_conn"
    ! prepare hdf5 write
    call h5gopen_f(file_id, gname_proc, group_id, error) ! group open
    if (error == 0) write(*,*) 'hdf5 passed at 106'
 
    ! create a dataset for elm_conn
    ddims = (/3+NGNOD, nspec_local/)
    call h5screate_simple_f(drank, ddims, dspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 107'
    call h5dcreate_f(group_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 108'
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, elm_conn, ddims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 109'
 
    ! create an attribute for npgeo
    call h5screate_simple_f(arank, adims, aspace_id, error) ! for a integernvalue
    if (error == 0) write(*,*) 'hdf5 passed at 110'
    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error) ! for a string tag of attribute value
    if (error == 0) write(*,*) 'hdf5 passed at 111'
    call h5tset_size_f(atype_id, taglen, error)
    if (error == 0) write(*,*) 'hdf5 passed at 112'
    call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 113'
    call h5awrite_f(attr_id, atype_id, nspec_local, data_dims, error) ! write
    if (error == 0) write(*,*) 'hdf5 passed at 114'
 
    ! close hdf5
    call h5sclose_f(dspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 115'
    call h5sclose_f(aspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 116'
    call h5dclose_f(dset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 117'
    call h5aclose_f(attr_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 118'
    call h5tclose_f(atype_id,  error)
    if (error == 0) write(*,*) 'hdf5 passed at 119'
 
    call h5gclose_f(group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 120'
 
  endif

  end subroutine write_partition_database_h5



  !--------------------------------------------------
  ! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_interfaces_database_h5(IIN_database, gname_proc, file_id, tab_interfaces, tab_size_interfaces, &
                              iproc, ninterfaces, &
                              my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, &
                              glob2loc_nodes_nparts, glob2loc_nodes_parts, &
                              glob2loc_nodes, num_phase, nparts)

  character(len=*), intent(in)  :: IIN_database
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

  ! variables for hdf5 write
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error

  ! variables for my_ninterfaces_and_maxval
  integer, dimension(2)          :: num_interface_and_max
  integer(HID_T)                 :: ndspace_id, ndset_id
  integer                        :: nrank = 1
  integer(HSIZE_T), dimension(1) :: ndims = (/2/)
  character(len=40)              :: ndsetname 
 
	! my_nb_interfaces in subgroup num_interface (2,my_ninterfaces)
  integer, dimension(:,:), allocatable :: num_neighbors_elmnts
  integer(HID_T)                     :: mdspace_id, mdset_id
  integer                            :: mrank = 2
  integer(HSIZE_T), dimension(2)     :: mdims
  character(len=40)                  :: mdsetname
 
  ! my_interfaces
  integer, dimension(:,:), allocatable :: neighbors_elmnts
  integer(HID_T)                       :: idspace_id, idset_id
  integer                              :: irank = 2
  integer(HSIZE_T), dimension(2)       :: idims
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
    ! initialize Fortran interface
    call h5gopen_f(file_id, gname_proc, group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 125'
 
    ! group names
    ndsetname = "my_ninterface_and_max"
    mdims = (/2,my_interfaces/)
    mdsetname = "my_nb_interfaces"
    idims = (/i,sum(my_nb_interfaces(:))/)
    idsetname = "my_interfaces"


    ! write my_ninterface andmaxval
    if (my_ninterface == 0) then
      num_interface_and_max = (/my_ninterface, 0/)
    else
      num_interface_and_max = (/my_ninterface, maxval(my_nb_interfaces)/)
    endif
    call h5screate_simple_f(nrank, ndims, ndspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 126'
    call h5dcreate_f(group_id, ndsetname, H5T_NATIVE_INTEGER, ndspace_id, ndset_id, error) ! a group with group_id need to be opened before.
    if (error == 0) write(*,*) 'hdf5 passed at 127'
    call h5dwrite_f(ndset_id, H5T_NATIVE_INTEGER, num_interface_and_max, ndims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 128'
 
    ! write my_nb_interfaces
    call h5screate_simple_f(mrank, mdims, mdspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 129'
    call h5dcreate_f(group_id, mdsetname, H5T_NATIVE_INTEGER, mdspace_id, mdset_id, error) ! a group with group_id need to be opened before.
    if (error == 0) write(*,*) 'hdf5 passed at 130'
    call h5dwrite_f(mdset_id, H5T_NATIVE_INTEGER, num_neighbors_elmnts, mdims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 131'
 
    ! write my_nb_interfaces
    call h5screate_simple_f(irank, idims, idspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 132'
    call h5dcreate_f(group_id, idsetname, H5T_NATIVE_INTEGER, idspace_id, idset_id, error) ! a group with group_id need to be opened before.
    if (error == 0) write(*,*) 'hdf5 passed at 133'
    call h5dwrite_f(idset_id, H5T_NATIVE_INTEGER, neighbors_elmnts, idims, error)
    if (error == 0) write(*,*) 'hdf5 passed at 134'
 
    ! close hdf5
    call h5sclose_f(ndspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 135'
    call h5dclose_f(ndset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 136'
    call h5sclose_f(mdspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 137'
    call h5dclose_f(mdset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 138'
    call h5sclose_f(idspace_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 139'
    call h5dclose_f(idset_id,   error)
    if (error == 0) write(*,*) 'hdf5 passed at 140'
 
    call h5gclose_f(group_id, error)
    if (error == 0) write(*,*) 'hdf5 passed at 141'
 
    ! allocate tempral array size
    deallocate(num_neighbors_elmnts,stat=ier); if (ier /= 0) stop 'Error deallocating array num_neighbors_elmnts'
    deallocate(neighbors_elmnts,stat=ier);     if (ier /= 0) stop 'Error deallocating array neighbors_elmnts'
 
  endif

  end subroutine write_interfaces_database_h5

  !--------------------------------------------------
  ! Write elements on surface boundaries (and their four nodes on boundaries)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_moho_surface_database_h5(IIN_database, gname_proc, file_id, iproc, nspec, &
                                         glob2loc_elmnts, glob2loc_nodes_nparts, &
                                         glob2loc_nodes_parts, glob2loc_nodes, part, &
                                         nspec2D_moho,ibelm_moho,nodes_ibelm_moho, NGNOD2D)

  character(len=*), intent(in)  :: IIN_database
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

  ! hdf5 variables
  integer(HID_T), intent(in) :: file_id
  integer(HID_T)             :: group_id
  integer                    :: error

  ! for loc_nspec_2d_homo attribute
  integer(HID_T)                 :: aspace_id, atype_id, attr_id
  integer                        :: arank = 1 ! attribute rank
  integer(HSIZE_T), dimension(1) :: adims = (/2/) ! attribute dimension
  integer(HSIZE_T)               :: taglen = 8
  character(len=8)               :: aname = "loc_moho"
  integer(HSIZE_T), dimension(1) :: attr_data_dims = 2
  integer         , dimension(2) :: moho_attr

  ! homo_elements dataset
  integer, dimension(:,:), allocatable :: loc_moho_temp
  integer(HID_T)                       :: dspace_id, dset_id
  integer                              :: rank = 2
  integer(HSIZE_T), dimension(2)       :: dims
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
  call h5gopen_f(file_id, gname_proc, group_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 144'
 
  ! create datasets and write
  call h5screate_simple_f(rank, dims, dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 145'
  call h5dcreate_f(group_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, error) ! a group with group_id need to be opened before.
  if (error == 0) write(*,*) 'hdf5 passed at 146'
  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, loc_moho_temp, dims, error)
  if (error == 0) write(*,*) 'hdf5 passed at 147'
  ! attribute
  call h5screate_simple_f(arank, adims, aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 148'
  call h5tcopy_f(h5t_native_character, atype_id, error) ! for a string tag of attribute value
  if (error == 0) write(*,*) 'hdf5 passed at 149'
  call h5tset_size_f(atype_id, taglen, error)
  if (error == 0) write(*,*) 'hdf5 passed at 150'
  call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)   ! dset_id need to opened before
  if (error == 0) write(*,*) 'hdf5 passed at 151'
  call h5awrite_f(attr_id, atype_id, moho_attr, attr_data_dims, error) ! attr_data need to be prepared before
  if (error == 0) write(*,*) 'hdf5 passed at 152'
 
  ! close hdf5 objects
  call h5sclose_f(dspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 153'
  call h5dclose_f(dset_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 154'
  call h5aclose_f(attr_id,   error)
  if (error == 0) write(*,*) 'hdf5 passed at 155'
  call h5tclose_f(atype_id,  error)
  if (error == 0) write(*,*) 'hdf5 passed at 156'
  call h5sclose_f(aspace_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 157'
 
  call h5gclose_f(group_id, error)
  if (error == 0) write(*,*) 'hdf5 passed at 158'
 
  deallocate(loc_moho_temp,stat=ier); if (ier /= 0) stop 'Error deallocating array loc_moho_temp'
 
  end subroutine write_moho_surface_database_h5


end module part_decompose_mesh_hdf5

