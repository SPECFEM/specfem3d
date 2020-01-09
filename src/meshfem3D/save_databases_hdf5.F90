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

  subroutine save_databases_hdf5(prname,nspec,nglob,iproc_xi,iproc_eta, &
                            NPROC_XI,NPROC_ETA,addressing,iMPIcut_xi,iMPIcut_eta, &
                            ibool,nodes_coords,ispec_material_id, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            NMATERIALS,material_properties, &
                            nspec_CPML,CPML_to_spec,CPML_regions,is_CPML, &
                            xstore, ystore, zstore)

  use constants, only: MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,SAVE_MESH_AS_CUBIT,NDIM,myrank
  use constants_meshfem3D, only: NGLLX_M,NGLLY_M,NGLLZ_M,NUMBER_OF_MATERIAL_PROPERTIES

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,NGNOD,NGNOD2D,LOCAL_PATH

  use meshfem3d_par, only: NPROC

  use phdf5_utils
  use my_mpi


  implicit none

  ! number of spectral elements in each block
  integer nspec

  ! number of vertices in each block
  integer nglob

  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8
  integer iproc, iproc_xi,iproc_eta
  integer NPROC_XI,NPROC_ETA
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)
  integer addressing(0:NPROC_XI-1,0:NPROC_ETA-1)

  ! arrays with the mesh
  integer ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec)
  double precision :: nodes_coords(nglob,NDIM)

  !! VM VM add all GLL points for Axisem coupling
  double precision, dimension(NGLLX_M,NGLLY_M,NGLLZ_M,nspec) :: xstore, ystore, zstore

  integer ispec_material_id(nspec)

  ! boundary parameters locator
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM)
  integer ibelm_top(NSPEC2D_TOP)

  ! material properties
  integer :: NMATERIALS
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
  double precision , dimension(NMATERIALS,NUMBER_OF_MATERIAL_PROPERTIES) :: material_properties
  double precision , dimension(17,NMATERIALS) :: mat_prop


  ! CPML
  integer, intent(in) :: nspec_CPML
  integer, dimension(nspec_CPML), intent(in) :: CPML_to_spec,CPML_regions
  logical, dimension(nspec), intent(in) :: is_CPML
  integer :: nspec_CPML_total,ispec_CPML

  double precision , dimension(17) :: matpropl
  integer :: i,ispec,iglob,ier

  ! name of the database files
  character(len=MAX_STRING_LEN) :: prname

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  character(len=MAX_STRING_LEN) :: IIN_database_hdf5

  integer :: ndef,nundef
  integer :: mat_id,domain_id
  ! there was a potential bug here if nspec is big
  !integer,dimension(2,nspec) :: material_index
  integer,dimension(:,:),allocatable            :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop


  ! mpi variables
  integer :: info, comm

  ! phdf5 io object
  type(h5io) :: h5

  ! name for material data group
  character(len=14), parameter :: gname_material   = "material_props" ! Group name for material properties
  ! processor dependent group names
  character(len=5)  :: gname_proc_head = "proc_"
  character(len=11) :: gname_proc
  character(len=10) :: tempstr
  character(len=64) :: group_name
  character(len=64) :: dset_name, attrname
 
  ! material
  ! for attribute count_def_mat and count_undef_mat
  character(len=13)              :: m_aname = "count_def_mat"
  character(len=15)              :: u_aname = "count_undef_mat"
 
  ! for dataset mat_prop, undef_mat_prop
  character(len=40)              :: mdsetname = "mat_prop"
  character(len=40)              :: udsetname = "undef_mat_prop"

  ! for node coords
  integer, dimension(nglob)               :: glob2loc_nodes_this_proc ! dummy
  double precision, dimension(NDIM,nglob) :: nodes_coords_this_proc

  ! element
  integer, dimension(NGNOD,nspec)    :: elm_conn
  integer, dimension(1+NGNOD,nspec)  :: elm_conn_xdmf
  integer, dimension(2,nspec)        :: mat_mesh
  integer, dimension(nspec)          :: ispec_local

  ! boundary
  integer, dimension(8)          :: n_elms_on_bound
  integer, dimension(1+NGNOD2D,nspec2D_xmin   &
                              +nspec2D_xmax   &
                              +nspec2D_ymin   &
                              +nspec2D_ymax   &
                              +NSPEC2D_BOTTOM &
                              +NSPEC2D_TOP) :: glob2loc_elms_this_proc

  ! cpml
  integer, dimension(2)            :: ncpmls
  integer, dimension(2,nspec_CPML) :: elements_cpml
  integer, dimension(nspec)        :: if_cpml

  ! interfaces
  integer, dimension(2) :: num_interface_and_max
  integer, dimension(:,:), allocatable :: num_neighbors_elmnts
  integer, dimension(:,:), allocatable :: neighbors_elmnts
  integer :: count1, count2


  ! opens output file
  IIN_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/'//'Database.h5'

  ! initialize hdf5 io
  h5 = h5io()
  call h5_init(h5, IIN_database_hdf5)
  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)
  ! set mpi info
  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)

  ! assignes material index
  ! format: (1,ispec) = #material_id , (2,ispec) = #material_definition
  allocate(material_index(2,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1346')
  if (ier /= 0) stop 'Error allocating array material_index'
  material_index (:,:) = 0
  do ispec = 1, nspec
    ! material id
    material_index(1,ispec) = ispec_material_id(ispec)
    ! material definition: 1 = interface type / 2 = tomography type
    if (ispec_material_id(ispec) > 0) then
      ! dummy value, not used any further
      material_index(2,ispec) = 1
    else
      ! negative material ids
      ! by default, assumes tomography model material
      ! (note: interface type not implemented yet...)
      material_index(2,ispec) = 2
    endif
  enddo

  ! Materials properties
  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do i = 1,NMATERIALS
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    mat_id = material_properties(i,8)
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1

    mat_prop(1,i) = material_properties(i,1)
    mat_prop(2,i) = material_properties(i,2)
    mat_prop(3,i) = material_properties(i,3)
    mat_prop(4,i) = material_properties(i,4)
    mat_prop(5,i) = material_properties(i,5)
    mat_prop(6,i) = material_properties(i,6)
    mat_prop(7,i) = material_properties(i,7)
  enddo
  !debug
  !print *,'materials def/undef: ',ndef,nundef

  ! writes out undefined materials
  do i = 1,NMATERIALS
    domain_id = material_properties(i,7)
    mat_id = material_properties(i,8)
    if (mat_id < 0) then
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      undef_mat_prop(:,:) = ''
      ! material id
      write(undef_mat_prop(1,1),*) mat_id
      ! name
      undef_mat_prop(2,1) = 'tomography'
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        undef_mat_prop(3,1) = 'acoustic'
      case (IDOMAIN_ELASTIC)
        undef_mat_prop(3,1) = 'elastic'
      end select
      ! default name
      undef_mat_prop(4,1) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      write(undef_mat_prop(5,1),*) 0
      ! domain-id
      write(undef_mat_prop(6,1),*) domain_id
      ! debug
      !print *,'undef mat: ',undef_mat_prop
    endif
  enddo

  if (myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(h5)
    ! prepare groups in hdf5 file
    ! material data group
    call h5_create_group(h5, gname_material)

    call h5_open_group(h5, gname_material)

    call h5_write_dataset_2d_d(h5, mdsetname, mat_prop)
    ! create an attribute for count_def_mat
    call h5_add_attribute_i(h5, m_aname, (/ndef/))
    call h5_close_dataset(h5)

    ! create a dataset for undef_mat_prop
    call h5_write_dataset_2d_c(h5, udsetname, undef_mat_prop)
    call h5_add_attribute_i(h5, u_aname, (/nundef/))
    call h5_close_dataset(h5)

    call h5_close_group(h5)
    call h5_close_file(h5)

  endif ! end write material info

  ! create hdf5 file and write datasets independently
  ! create file, groups and datasets of all mpi rank
  call h5_open_file_p(h5)
  ! create group
  do iproc = 0, NPROC-1
    write(tempstr, "(i6.6)") iproc
    group_name = gname_proc_head // trim(tempstr)

    call h5_create_group(h5, group_name)
  enddo
  call synchronize_all()
  call h5_close_file(h5)
  ! open file
!  call h5_open_file_p_collect(h5) ! collective writing is faster but mpi error occurs
! when using more than 128 cpus
  call h5_open_file_p(h5)

  ! global nodes
  do iglob=1,nglob
    nodes_coords_this_proc(1,iglob) = nodes_coords(iglob,1)
    nodes_coords_this_proc(2,iglob) = nodes_coords(iglob,2)
    nodes_coords_this_proc(3,iglob) = nodes_coords(iglob,3)
    glob2loc_nodes_this_proc(iglob) = 0 ! dummy
  enddo

  dset_name = "nodes_coords"
  call h5_write_dataset_p_2d_d(h5, dset_name, nodes_coords_this_proc)
  dset_name = "glob2loc_nodes"
  call h5_write_dataset_p_1d_i(h5, dset_name, glob2loc_nodes_this_proc)
  dset_name = "nnodes_loc"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nglob/))


  ! spectral elements
  do ispec=1,nspec
    ispec_local(ispec) = ispec ! <- ispec_local is dummy

    mat_mesh(1,ispec) = material_index(1,ispec)
    mat_mesh(2,ispec) = material_index(2,ispec) ! <- mat_mesh

    elm_conn(1,ispec) = ibool(1,1,1,ispec)
    elm_conn(2,ispec) = ibool(2,1,1,ispec)
    elm_conn(3,ispec) = ibool(2,2,1,ispec)
    elm_conn(4,ispec) = ibool(1,2,1,ispec)
    elm_conn(5,ispec) = ibool(1,1,2,ispec)
    elm_conn(6,ispec) = ibool(2,1,2,ispec)
    elm_conn(7,ispec) = ibool(2,2,2,ispec)
    elm_conn(8,ispec) = ibool(1,2,2,ispec)

    elm_conn_xdmf(1,ispec) = 9
    elm_conn_xdmf(2,ispec) = ibool(1,1,1,ispec)
    elm_conn_xdmf(3,ispec) = ibool(2,1,1,ispec)
    elm_conn_xdmf(4,ispec) = ibool(2,2,1,ispec)
    elm_conn_xdmf(5,ispec) = ibool(1,2,1,ispec)
    elm_conn_xdmf(6,ispec) = ibool(1,1,2,ispec)
    elm_conn_xdmf(7,ispec) = ibool(2,1,2,ispec)
    elm_conn_xdmf(8,ispec) = ibool(2,2,2,ispec)
    elm_conn_xdmf(9,ispec) = ibool(1,2,2,ispec) ! <- elm_conn, elm_conn_xdmf
  enddo

  ! create a dataset for elm_conn
  dset_name = "elm_conn"
  call h5_write_dataset_p_2d_i(h5, dset_name, elm_conn)

  dset_name = "nspec_local"
  call h5_write_dataset_p_1d_i(h5, dset_name, (/nspec/))

  ! create a dataset for elm_conn_xdmf
  dset_name = "elm_conn_xdmf"
  call h5_write_dataset_p_2d_i(h5, dset_name, elm_conn_xdmf)

  ! write mat_mesh
  dset_name = "mat_mesh"
  call h5_write_dataset_p_2d_i(h5, dset_name, mat_mesh)

  ! write ispec_local
  dset_name = "ispec_local"
  call h5_write_dataset_p_1d_i(h5, dset_name, ispec_local)

  ! Boundaries
  n_elms_on_bound(1) = nspec2D_xmin
  n_elms_on_bound(2) = nspec2D_xmax
  n_elms_on_bound(3) = nspec2D_ymin
  n_elms_on_bound(4) = nspec2D_ymax
  n_elms_on_bound(5) = NSPEC2D_BOTTOM
  n_elms_on_bound(6) = NSPEC2D_TOP    ! <- n_elms_on_bound

  count1 = 1

  do i=1,nspec2D_xmin
    glob2loc_elms_this_proc(1,count1) = ibelm_xmin(i)
    glob2loc_elms_this_proc(2,count1) = ibool(1,1,1,ibelm_xmin(i))
    glob2loc_elms_this_proc(3,count1) = ibool(1,NGLLY_M,1,ibelm_xmin(i))
    glob2loc_elms_this_proc(4,count1) = ibool(1,1,NGLLZ_M,ibelm_xmin(i))
    glob2loc_elms_this_proc(5,count1) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_xmin(i))
    count1 = count1 + 1
  enddo
  do i=1,nspec2D_xmax
    glob2loc_elms_this_proc(1,count1) = ibelm_xmax(i)
    glob2loc_elms_this_proc(2,count1) = ibool(NGLLX_M,1,1,ibelm_xmax(i))
    glob2loc_elms_this_proc(3,count1) = ibool(NGLLX_M,NGLLY_M,1,ibelm_xmax(i))
    glob2loc_elms_this_proc(4,count1) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_xmax(i))
    glob2loc_elms_this_proc(5,count1) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_xmax(i))
    count1 = count1 + 1
  enddo
  do i=1,nspec2D_ymin
    glob2loc_elms_this_proc(1,count1) = ibelm_ymin(i)
    glob2loc_elms_this_proc(2,count1) = ibool(1,1,1,ibelm_ymin(i))
    glob2loc_elms_this_proc(3,count1) = ibool(NGLLX_M,1,1,ibelm_ymin(i))
    glob2loc_elms_this_proc(4,count1) = ibool(1,1,NGLLZ_M,ibelm_ymin(i))
    glob2loc_elms_this_proc(5,count1) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_ymin(i))
    count1 = count1 + 1
  enddo
  do i=1,nspec2D_ymax
    glob2loc_elms_this_proc(1,count1) = ibelm_ymax(i)
    glob2loc_elms_this_proc(2,count1) = ibool(NGLLX_M,NGLLY_M,1,ibelm_ymax(i))
    glob2loc_elms_this_proc(3,count1) = ibool(1,NGLLY_M,1,ibelm_ymax(i))
    glob2loc_elms_this_proc(4,count1) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
    glob2loc_elms_this_proc(5,count1) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_ymax(i))
    count1 = count1 + 1
  enddo
  do i=1,NSPEC2D_BOTTOM
    glob2loc_elms_this_proc(1,count1) = ibelm_bottom(i)
    glob2loc_elms_this_proc(2,count1) = ibool(1,1,1,ibelm_bottom(i))
    glob2loc_elms_this_proc(3,count1) = ibool(NGLLX_M,1,1,ibelm_bottom(i))
    glob2loc_elms_this_proc(4,count1) = ibool(NGLLX_M,NGLLY_M,1,ibelm_bottom(i))
    glob2loc_elms_this_proc(5,count1) = ibool(1,NGLLY_M,1,ibelm_bottom(i))
    count1 = count1 + 1
  enddo
  do i=1,NSPEC2D_TOP
    glob2loc_elms_this_proc(1,count1) = ibelm_top(i)
    glob2loc_elms_this_proc(2,count1) = ibool(1,1,NGLLZ_M,ibelm_top(i))
    glob2loc_elms_this_proc(3,count1) = ibool(NGLLX_M,1,NGLLZ_M,ibelm_top(i))
    glob2loc_elms_this_proc(4,count1) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ibelm_top(i))
    glob2loc_elms_this_proc(5,count1) = ibool(1,NGLLY_M,NGLLZ_M,ibelm_top(i))
    count1 = count1 + 1
  enddo
  ! <- glob2loc_elms

  ! write n_elms_on_bound
  dset_name = "n_elms_on_bound"
  call h5_write_dataset_p_1d_i(h5, dset_name, n_elms_on_bound)

  ! write glob2loc_elms_this_proc
  dset_name = "glob2loc_elms"
  call h5_write_dataset_p_2d_i(h5, dset_name, glob2loc_elms_this_proc)


  ! CPML
  call sum_all_i(nspec_CPML,nspec_CPML_total)
  call synchronize_all()
  call bcast_all_singlei(nspec_CPML_total)
  call synchronize_all()

  ncpmls(1) = nspec_CPML_total
  ncpmls(2) = nspec_CPML
  if (nspec_CPML_total > 0) then
     do ispec_CPML=1,nspec_CPML
        elements_cpml(1,ispec_CPML) = CPML_to_spec(ispec_CPML)
        elements_cpml(2,ispec_CPML) = CPML_regions(ispec_CPML)
     enddo
     do ispec=1,nspec
        if(is_CPML(ispec)) then
          if_cpml(ispec) = 1
        else
          if_cpml(ispec) = 0
        endif
     enddo

    ! create a dataset for elements_cpml
    dset_name = "elements_cpml"
    call h5_write_dataset_p_2d_i(h5, dset_name, elements_cpml)
    ! create a dataset for if_cpml
      ! strangely here H5T_NATIVE_HBOOL may not be used.
    dset_name = "if_cpml"
    call h5_write_dataset_p_1d_i(h5, dset_name, if_cpml)

  endif
  ! create a dataset for nspec_cpml and nspec_cpml_local
  dset_name = "nspec_cpml_globloc"
  call h5_write_dataset_p_1d_i(h5, dset_name, ncpmls)
 




  ! MPI Interfaces
  if (NPROC_XI >= 2 .or. NPROC_ETA >= 2) then
    ! determines number of MPI interfaces for each slice
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.

    ! slices at model boundaries
    if (iproc_xi == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(W) = .false.
    endif
    if (iproc_xi == NPROC_XI-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(E) = .false.
    endif
    if (iproc_eta == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(S) = .false.
    endif
    if (iproc_eta == NPROC_ETA-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(N) = .false.
    endif

    ! slices in middle of model
    if ((interfaces(W) .eqv. .true.) .and. (interfaces(N) .eqv. .true.)) then
      interfaces(NW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(N) .eqv. .true.) .and. (interfaces(E) .eqv. .true.)) then
      interfaces(NE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(E) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
      interfaces(SE) = .true.
      nb_interfaces =  nb_interfaces +1
    endif
    if ((interfaces(W) .eqv. .true.) .and. (interfaces(S) .eqv. .true.)) then
      interfaces(SW) = .true.
      nb_interfaces =  nb_interfaces +1
    endif

    nspec_interface(:) = 0
    if (interfaces(W))   nspec_interface(W)  = count(iMPIcut_xi(1,:) .eqv. .true.)
    if (interfaces(E))   nspec_interface(E)  = count(iMPIcut_xi(2,:) .eqv. .true.)
    if (interfaces(S))   nspec_interface(S)  = count(iMPIcut_eta(1,:) .eqv. .true.)
    if (interfaces(N))   nspec_interface(N)  = count(iMPIcut_eta(2,:) .eqv. .true.)
    if (interfaces(NW))  nspec_interface(NW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(NE))  nspec_interface(NE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(2,:) .eqv. .true.))
    if (interfaces(SE))  nspec_interface(SE) = count((iMPIcut_xi(2,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))
    if (interfaces(SW))  nspec_interface(SW) = count((iMPIcut_xi(1,:) .eqv. .true.) .and. (iMPIcut_eta(1,:) .eqv. .true.))

    nspec_interfaces_max = maxval(nspec_interface)

    ! write my_ninterface and maxval
    num_interface_and_max = (/nb_interfaces, nspec_interfaces_max/)

    ! allocate tempral array size
    allocate(num_neighbors_elmnts(2,nb_interfaces),stat=ier)
    if (ier /= 0) stop 'Error allocating array num_neighbors_elmnts'
    allocate(neighbors_elmnts(6, sum(nspec_interface(:))))
    if (ier /= 0) stop 'Error allocating array neighbors_elmnts'

    count1 = 1
    count2 = 1

    if (interfaces(W)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi-1,iproc_eta)
      num_neighbors_elmnts(2,count1) = nspec_interface(W)
      count1 = count1+1
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,2,1,ispec)
          neighbors_elmnts(5,count2) = ibool(1,1,2,ispec)
          neighbors_elmnts(6,count2) = ibool(1,2,2,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(E)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi+1,iproc_eta)
      num_neighbors_elmnts(2,count1) = nspec_interface(E)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(2,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(2,2,1,ispec)
          neighbors_elmnts(5,count2) = ibool(2,1,2,ispec)
          neighbors_elmnts(6,count2) = ibool(2,2,2,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(S)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi,iproc_eta-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(S)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(2,1,1,ispec)
          neighbors_elmnts(5,count2) = ibool(1,1,2,ispec)
          neighbors_elmnts(6,count2) = ibool(2,1,2,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(N)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi,iproc_eta+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(N)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(2,2,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,2,1,ispec)
          neighbors_elmnts(5,count2) = ibool(2,2,2,ispec)
          neighbors_elmnts(6,count2) = ibool(1,2,2,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(NW)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi-1,iproc_eta+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(NW)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(1,2,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,2,2,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi+1,iproc_eta+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(NE)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(2,2,1,ispec)
          neighbors_elmnts(4,count2) = ibool(2,2,2,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi+1,iproc_eta-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(SE)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(2,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(2,1,2,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi-1,iproc_eta-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(SW)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,1,2,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

  else

    ! only one slice, no MPI interfaces
    num_interface_and_max = (/0, 0/)
 
  endif

  dset_name = "my_ninterface_and_max"
  call h5_write_dataset_p_1d_i(h5, dset_name, num_interface_and_max)

  ! write my_nb_interfaces
  dset_name = "my_nb_interfaces"
  call h5_write_dataset_p_2d_i(h5, dset_name, num_neighbors_elmnts)

  ! write my_interfaces
  dset_name = "my_interfaces"
  call h5_write_dataset_p_2d_i(h5, dset_name, neighbors_elmnts)
 
  deallocate(num_neighbors_elmnts, neighbors_elmnts)

  call synchronize_all()
  call h5_close_file(h5)

  ! CUBIT output
  if (SAVE_MESH_AS_CUBIT) then
    ! only for single process at the moment
    if (NPROC_XI == 1 .and. NPROC_ETA == 1) then
      !! VM VM add outputs as CUBIT
      call save_output_mesh_files_as_cubit(nspec,nglob, ibool,nodes_coords, ispec_material_id, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NMATERIALS,material_properties, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)
      ! output for AxiSEM coupling
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        call save_output_mesh_files_for_coupled_model(nspec, &
                                           nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                           ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                           xstore,ystore,zstore)
      endif
    endif
  endif

  deallocate(material_index)

  end subroutine save_databases_hdf5

