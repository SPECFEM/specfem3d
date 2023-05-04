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

  subroutine save_databases_hdf5(nspec,nglob, &
                                 iMPIcut_xi,iMPIcut_eta, &
                                 nodes_coords,ispec_material_id, &
                                 nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                 ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top)


  use constants, only: NDIM
  use meshfem_par, only: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

#ifdef USE_HDF5
  use constants, only: IMAIN,MAX_STRING_LEN,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    myrank

  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,NGNOD,NGNOD2D,LOCAL_PATH

  use meshfem_par, only: ibool,xstore,ystore,zstore, &
    NPROC, &
    NMATERIALS,material_properties,material_properties_undef, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    addressing, &
    iproc_xi_current,iproc_eta_current, &
    NPROC_XI,NPROC_ETA, &
    SAVE_MESH_AS_CUBIT

  use manager_hdf5
#endif

  implicit none

  ! number of spectral elements in each block
  integer, intent(in) :: nspec

  ! number of vertices in each block
  integer, intent(in) :: nglob

  logical, intent(in) :: iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  ! arrays with the mesh
  double precision, intent(in) :: nodes_coords(nglob,NDIM)
  integer, intent(in) :: ispec_material_id(nspec)

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax
  integer, intent(in) :: ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer, intent(in) :: ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer, intent(in) :: ibelm_bottom(NSPEC2D_BOTTOM)
  integer, intent(in) :: ibelm_top(NSPEC2D_TOP)

#ifdef USE_HDF5
  ! local parameters
  ! MPI Cartesian topology
  ! E for East (= XI_MIN), W for West (= XI_MAX), S for South (= ETA_MIN), N for North (= ETA_MAX)
  integer, parameter :: W=1,E=2,S=3,N=4,NW=5,NE=6,SE=7,SW=8

  ! CPML
  integer :: nspec_CPML_total,ispec_CPML

  ! material properties
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
  double precision , dimension(17,NMATERIALS) :: mat_prop

  ! for MPI interfaces
  integer ::  nb_interfaces,nspec_interfaces_max
  logical, dimension(8) ::  interfaces
  integer, dimension(8) ::  nspec_interface

  integer :: ndef,nundef
  integer :: mat_id,domain_id
  ! there was a potential bug here if nspec is big
  integer,dimension(:,:),allocatable            :: material_index
  character(len=MAX_STRING_LEN), dimension(6,1) :: undef_mat_prop

  ! MPI variables
  integer :: info, comm

  ! name for material data group
  character(len=14), parameter :: gname_material   = "material_props" ! Group name for material properties
  ! processor dependent group names
  character(len=64) :: dset_name

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
  integer, dimension(6)          :: n_elms_on_bound
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

  ! temporary array for local nodes (either NGNOD2D or NGNOD)
  integer, dimension(NGNOD) :: loc_node
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz
  integer :: ia,inode

  ! offset arrays
  integer, dimension(0:NPROC-1) :: offset_nnodes
  integer, dimension(0:NPROC-1) :: offset_nelems
  integer, dimension(0:NPROC-1) :: offset_nelems_cpml
  integer, dimension(0:NPROC-1) :: offset_n_elms_bounds
  integer, dimension(0:NPROC-1) :: offset_n_elms_interface
  integer, dimension(0:NPROC-1) :: offset_nb_interfaces

  integer :: i,ispec,iglob,ier

  ! opens output file
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/'//'Database.h5'

  ! get MPI parameters
  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initialize hdf5 io
  call h5_initialize()
  ! set MPI info
  call h5_set_mpi_info(comm, info, myrank, NPROC)

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
    ! material_properties(:,:) array:
    !   first dimension  : material_id
    !   second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    !
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'mesh files:'
    write(IMAIN,*) '  saving file : ',trim(name_database_hdf5)
    write(IMAIN,*) '  using HDF5 file format'
    call flush_IMAIN()
  endif

  ! writes out defined materials
  do i = 1,NMATERIALS
    ! material properties format: #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id  #material_id
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id > 0) then
      ! pad dummy zeros to fill up 17 entries
      mat_prop(:,i) = 0.d0
      select case(domain_id)
      case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC)
        ! material properties format:
        !#(1)rho  #(2)vp #(3)vs #(4)Q_Kappa #(5)Q_mu #(6)anisotropy_flag #(7)domain_id #(8)mat_id
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag,material_domain_id
        !
        ! skipping mat_id, not needed
        mat_prop(1:7,i) = material_properties(i,1:7)
      case (IDOMAIN_POROELASTIC)
        ! material properties format:
        !#(1)rho_s #(2)rho_f #(3)phi #(4)tort #(5)eta #(6)0 #(7)domain_id #(8)mat_id
        !            .. #(9)kxx #(10)kxy #(11)kxz #(12)kyy #(13)kyz #(14)kzz #(15)kappa_s #(16)kappa_f #(17)kappa_fr #(18)mu_fr
        !
        ! output format for xgenerate_database (same as for cubit/trelis inputs):
        !   rhos,rhof,phi,tort,eta,0,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
        mat_prop(1:7,i) = material_properties(i,1:7)
        ! skipping mat_id, not needed
        mat_prop(8:17,i) = material_properties(i,9:18)
      end select
    endif
  enddo

  ! writes out undefined materials
  do i = 1,NMATERIALS
    domain_id = int(material_properties(i,7))
    mat_id = int(material_properties(i,8))
    if (mat_id < 0) then
      ! format:
      ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
      ! format example tomography: -1 tomography elastic tomography_model.xyz 0 2
      undef_mat_prop(:,:) = ''
      ! material id
      write(undef_mat_prop(1,1),*) mat_id
      ! keyword/domain/filename
      undef_mat_prop(2,1) = material_properties_undef(i,1)  ! type-keyword : tomography/interface
      undef_mat_prop(3,1) = material_properties_undef(i,2)  ! domain-name  : acoustic/elastic/poroelastic
      undef_mat_prop(4,1) = material_properties_undef(i,3)  ! tomo-filename: tomography_model**.xyz
      ! checks consistency between domain-name and domain_id
      select case (domain_id)
      case (IDOMAIN_ACOUSTIC)
        if (trim(undef_mat_prop(3,1)) /= 'acoustic') stop 'Error in undef_mat_prop acoustic domain'
      case (IDOMAIN_ELASTIC)
        if (trim(undef_mat_prop(3,1)) /= 'elastic')  stop 'Error in undef_mat_prop elastic domain'
      case (IDOMAIN_POROELASTIC)
        if (trim(undef_mat_prop(3,1)) /= 'poroelastic')  stop 'Error in undef_mat_prop poroelastic domain'
      end select
      ! default name if none given
      if (trim(undef_mat_prop(4,1)) == "") undef_mat_prop(4,1) = 'tomography_model.xyz'
      ! default tomo-id (unused)
      write(undef_mat_prop(5,1),*) 0
      ! domain-id
      write(undef_mat_prop(6,1),*) domain_id
      ! debug
      !print *,'undef mat: ',undef_mat_prop
    endif
  enddo

  ! share the offset info
  !
  ! offset_nnodes
  ! offset_nelems
  ! offset_n_elms_bounds
  ! offset_n_elms_interface (later)
  ! offset_nb_interfaces (late)
  call gather_all_all_singlei(nglob, offset_nnodes, NPROC)
  call gather_all_all_singlei(nspec, offset_nelems, NPROC)
  call gather_all_all_singlei(nspec_CPML, offset_nelems_cpml, NPROC)
  call gather_all_all_singlei( nspec2D_xmin   &
                              +nspec2D_xmax   &
                              +nspec2D_ymin   &
                              +nspec2D_ymax   &
                              +NSPEC2D_BOTTOM &
                              +NSPEC2D_TOP, offset_n_elms_bounds, NPROC)

  if (myrank == 0) then
    ! create a hdf5 file
    call h5_create_file(name_database_hdf5)
    ! prepare groups in hdf5 file
    ! material data group
    call h5_create_group(gname_material)

    call h5_open_group(gname_material)

    call h5_write_dataset(mdsetname, mat_prop)
    ! create an attribute for count_def_mat
    call h5_add_attribute_i(m_aname, (/ndef/))
    call h5_close_dataset()

    ! create a dataset for undef_mat_prop
    call h5_write_dataset(udsetname, undef_mat_prop)
    call h5_add_attribute_i(u_aname, (/nundef/))
    call h5_close_dataset()

    call h5_close_group()

    ! create datasets

    ! offset arrays
    dset_name = "offset_nnodes"
    call h5_write_dataset_no_group(dset_name, offset_nnodes)
    dset_name = "offset_nelems"
    call h5_write_dataset_no_group(dset_name, offset_nelems)
    dset_name = "offset_nelems_cpml"
    call h5_write_dataset_no_group(dset_name, offset_nelems_cpml)
    dset_name = "offset_n_elms_bounds"
    call h5_write_dataset_no_group(dset_name, offset_n_elms_bounds)

    ! database arrays
    dset_name = "nodes_coords"
    call h5_create_dataset_gen(dset_name,(/3,sum(offset_nnodes(:))/), 2, 8) ! double
    dset_name = "glob2loc_nodes"
    call h5_create_dataset_gen(dset_name,(/sum(offset_nnodes(:))/), 1, 1)
    dset_name = "nnodes_loc"
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "elm_conn"
    call h5_create_dataset_gen(dset_name,(/NGNOD,sum(offset_nelems(:))/), 2, 1)
    dset_name = "elm_conn_xdmf"
    call h5_create_dataset_gen(dset_name,(/NGNOD+1,sum(offset_nelems(:))/), 2, 1)
    dset_name = "nspec_local"
    call h5_create_dataset_gen(dset_name,(/NPROC/), 1, 1)
    dset_name = "mat_mesh"
    call h5_create_dataset_gen(dset_name,(/2,sum(offset_nelems(:))/), 2, 1)
    dset_name = "ispec_local"
    call h5_create_dataset_gen(dset_name,(/sum(offset_nelems(:))/), 1, 1)
    dset_name = "n_elms_on_bound"
    call h5_create_dataset_gen(dset_name,(/6*NPROC/), 1, 1)
    dset_name = "glob2loc_elms"
    call h5_create_dataset_gen(dset_name, (/1+NGNOD2D,sum(offset_n_elms_bounds(:))/),2,1)
    dset_name = "elements_cpml"
    call h5_create_dataset_gen(dset_name, (/2,sum(offset_nelems_cpml(:))/),2,1)
    dset_name = "if_cpml"
    call h5_create_dataset_gen(dset_name, (/sum(offset_nelems(:))/),1,1)
    dset_name = "nspec_cpml_globloc"
    call h5_create_dataset_gen(dset_name, (/2*NPROC/),1,1)
    dset_name = "my_ninterface_and_max"
    call h5_create_dataset_gen(dset_name,(/NPROC*2/), 1, 1)

    call h5_close_file()

  endif ! end write material info

  call synchronize_all()

  ! open file
  call h5_open_file_p_collect(name_database_hdf5) ! collective writing is faster but MPI error occurs
  ! when using more than 128 cpus
  !call h5_open_file_p(name_database_hdf5)

  ! global nodes
  do iglob=1,nglob
    nodes_coords_this_proc(1,iglob) = nodes_coords(iglob,1)
    nodes_coords_this_proc(2,iglob) = nodes_coords(iglob,2)
    nodes_coords_this_proc(3,iglob) = nodes_coords(iglob,3)
    glob2loc_nodes_this_proc(iglob) = 0 ! dummy
  enddo

  dset_name = "nodes_coords"
  call h5_write_dataset_collect_hyperslab(dset_name, nodes_coords_this_proc, (/0,sum(offset_nnodes(0:myrank-1))/), .true.)
  dset_name = "glob2loc_nodes"
  call h5_write_dataset_collect_hyperslab(dset_name, glob2loc_nodes_this_proc, (/sum(offset_nnodes(0:myrank-1))/), .true.)
  dset_name = "nnodes_loc"
  call h5_write_dataset_collect_hyperslab(dset_name, (/nglob/), (/myrank/), .true.)

  ! sets up node addressing
  call hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX_M,NGLLY_M,NGLLZ_M)

  ! spectral elements
  do ispec = 1,nspec
    ispec_local(ispec) = ispec ! < - ispec_local is dummy

    ! gets anchor nodes
    do ia = 1,NGNOD
      iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
      loc_node(ia) = iglob
    enddo

    mat_mesh(1,ispec) = material_index(1,ispec)
    mat_mesh(2,ispec) = material_index(2,ispec) ! < - mat_mesh

    elm_conn_xdmf(1,ispec) = 9
    do ia = 1, NGNOD
      elm_conn(ia,ispec) = loc_node(ia) ! for generate_databases
      elm_conn_xdmf(ia+1,ispec) = loc_node(ia)-1 ! for hdf5 visualization starting 0 in vtk format
    enddo

  enddo

  ! create a dataset for elm_conn
  dset_name = "elm_conn"
  call h5_write_dataset_collect_hyperslab(dset_name, elm_conn, (/0,sum(offset_nelems(0:myrank-1))/),.true.)
  dset_name = "nspec_local"
  call h5_write_dataset_collect_hyperslab(dset_name, (/nspec/), (/myrank/), .true.)
  ! create a dataset for elm_conn_xdmf
  dset_name = "elm_conn_xdmf"
  call h5_write_dataset_collect_hyperslab(dset_name, elm_conn_xdmf,(/0,sum(offset_nelems(0:myrank-1))/),.true.)

  ! write mat_mesh
  dset_name = "mat_mesh"
  call h5_write_dataset_collect_hyperslab(dset_name, mat_mesh, (/0,sum(offset_nelems(0:myrank-1))/),.true.)

  ! write ispec_local
  dset_name = "ispec_local"
  call h5_write_dataset_collect_hyperslab(dset_name, ispec_local, (/sum(offset_nelems(0:myrank-1))/),.true.)

  ! Boundaries
  n_elms_on_bound(1) = nspec2D_xmin
  n_elms_on_bound(2) = nspec2D_xmax
  n_elms_on_bound(3) = nspec2D_ymin
  n_elms_on_bound(4) = nspec2D_ymax
  n_elms_on_bound(5) = NSPEC2D_BOTTOM
  n_elms_on_bound(6) = NSPEC2D_TOP    ! < - n_elms_on_bound

  count1 = 1

  do i = 1,nspec2D_xmin
    ispec = ibelm_xmin(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on xmin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmin'

    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo

  do i = 1,nspec2D_xmax
    ispec = ibelm_xmax(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on xmax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iax(ia) == NGLLX_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for xmax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for xmax'

    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo

  do i = 1,nspec2D_ymin
    ispec = ibelm_ymin(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on ymin
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymin'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymin'
    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo

  do i = 1,nspec2D_ymax
    ispec = ibelm_ymax(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on ymax
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iay(ia) == NGLLY_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for ymax'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for ymax'
    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo

  do i = 1,NSPEC2D_BOTTOM
    ispec = ibelm_bottom(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on bottom
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == 1) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for bottom'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for bottom'
    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo

  do i = 1,NSPEC2D_TOP
    ispec = ibelm_top(i)
    glob2loc_elms_this_proc(1,count1) = ispec
    ! gets anchor nodes on top
    inode = 0
    do ia = 1,NGNOD
      if (anchor_iaz(ia) == NGLLZ_M) then
        iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
        inode = inode + 1
        if (inode > NGNOD2D) stop 'inode index exceeds NGNOD2D for top'
        loc_node(inode) = iglob
      endif
    enddo
    if (inode /= NGNOD2D) stop 'Invalid number of inodes found for top'
    do inode=1,NGNOD2D
      glob2loc_elms_this_proc(inode+1,count1) = loc_node(inode)
    enddo
    count1 = count1 + 1
  enddo
  ! <- glob2loc_elms

  ! write n_elms_on_bound
  dset_name = "n_elms_on_bound"
  call h5_write_dataset_collect_hyperslab(dset_name, n_elms_on_bound, (/6*myrank/), .true.)
  ! write glob2loc_elms_this_proc
  dset_name = "glob2loc_elms"
  call h5_write_dataset_collect_hyperslab(dset_name, glob2loc_elms_this_proc, &
            (/0,sum(offset_n_elms_bounds(0:myrank-1))/),.true.)

  ! CPML
  !
  ! note: check with routine write_cpml_database() to produce identical output
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
        if (is_CPML(ispec)) then
          if_cpml(ispec) = 1
        else
          if_cpml(ispec) = 0
        endif
     enddo

    ! create a dataset for elements_cpml
    dset_name = "elements_cpml"
    call h5_write_dataset_collect_hyperslab(dset_name, elements_cpml, (/0,sum(offset_nelems_cpml(0:myrank-1))/), .true.)

    ! create a dataset for if_cpml
      ! strangely here H5T_NATIVE_HBOOL may not be used.
    dset_name = "if_cpml"
    call h5_write_dataset_collect_hyperslab(dset_name, if_cpml, (/sum(offset_nelems(0:myrank-1))/), .true.)

  endif

  ! create a dataset for nspec_cpml and nspec_cpml_local
  dset_name = "nspec_cpml_globloc"
  call h5_write_dataset_collect_hyperslab(dset_name, ncpmls, (/2*myrank/), .true.)

  ! MPI Interfaces
  !
  ! note: check with routine write_interfaces_database() to produce identical output
  if (NPROC_XI > 1 .or. NPROC_ETA > 1) then
    ! determines number of MPI interfaces for each slice
    nb_interfaces = 4
    interfaces(W:N) = .true.
    interfaces(NW:SW) = .false.

    ! slices at model boundaries
    if (iproc_xi_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(W) = .false.
    endif
    if (iproc_xi_current == NPROC_XI-1) then
      nb_interfaces =  nb_interfaces -1
      interfaces(E) = .false.
    endif
    if (iproc_eta_current == 0) then
      nb_interfaces =  nb_interfaces -1
      interfaces(S) = .false.
    endif
    if (iproc_eta_current == NPROC_ETA-1) then
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


    ! prepare two more offset arrays for
    ! offset_nb_interfaces
    ! offset_n_elms_interface
    call gather_all_all_singlei(nb_interfaces,offset_nb_interfaces,NPROC)
    call gather_all_all_singlei(sum(nspec_interface(:)),offset_n_elms_interface,NPROC)

    dset_name = "offset_n_elms_interface"
    call h5_write_dataset_no_group(dset_name, offset_n_elms_interface)
    dset_name = "offset_nb_interfaces"
    call h5_write_dataset_no_group(dset_name, offset_nb_interfaces)

    ! create dataset by main process
    call h5_close_file_p()

    if (myrank == 0) then
      call h5_open_file(name_database_hdf5)
      dset_name = "my_nb_interfaces"
      call h5_create_dataset_gen(dset_name,(/2, sum(offset_nb_interfaces(:))/), 2, 1)
      dset_name = "my_interfaces"
      call h5_create_dataset_gen(dset_name,(/6, sum(offset_n_elms_interface(:))/), 2, 1)
      call h5_close_file()
    endif

    call synchronize_all()

    call h5_open_file_p_collect(name_database_hdf5)

    count1 = 1
    count2 = 1

    ! face elements
    if (interfaces(W)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current-1,iproc_eta_current)
      num_neighbors_elmnts(2,count1) = nspec_interface(W)
      count1 = count1+1
      do ispec = 1,nspec
        if (iMPIcut_xi(1,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,NGLLY_M,1,ispec)
          neighbors_elmnts(5,count2) = ibool(1,1,NGLLZ_M,ispec)
          neighbors_elmnts(6,count2) = ibool(1,NGLLY_M,NGLLZ_M,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(E)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current+1,iproc_eta_current)
      num_neighbors_elmnts(2,count1) = nspec_interface(E)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_xi(2,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(NGLLX_M,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(NGLLX_M,NGLLY_M,1,ispec)
          neighbors_elmnts(5,count2) = ibool(NGLLX_M,1,NGLLZ_M,ispec)
          neighbors_elmnts(6,count2) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(S)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current,iproc_eta_current-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(S)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_eta(1,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(NGLLX_M,1,1,ispec)
          neighbors_elmnts(5,count2) = ibool(1,1,NGLLZ_M,ispec)
          neighbors_elmnts(6,count2) = ibool(NGLLX_M,1,NGLLZ_M,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(N)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current,iproc_eta_current+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(N)
      count1 = count1 + 1
      do ispec = 1,nspec
        if (iMPIcut_eta(2,ispec)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 4
          neighbors_elmnts(3,count2) = ibool(NGLLX_M,NGLLY_M,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,NGLLY_M,1,ispec)
          neighbors_elmnts(5,count2) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
          neighbors_elmnts(6,count2) = ibool(1,NGLLY_M,NGLLZ_M,ispec)
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(NW)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current-1,iproc_eta_current+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(NW)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(1,NGLLY_M,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,NGLLY_M,NGLLZ_M,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(NE)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current+1,iproc_eta_current+1)
      num_neighbors_elmnts(2,count1) = nspec_interface(NE)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(2,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(NGLLX_M,NGLLY_M,1,ispec)
          neighbors_elmnts(4,count2) = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(SE)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current+1,iproc_eta_current-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(SE)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(2,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(NGLLX_M,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(NGLLX_M,1,NGLLZ_M,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    if (interfaces(SW)) then
      num_neighbors_elmnts(1,count1) = addressing(iproc_xi_current-1,iproc_eta_current-1)
      num_neighbors_elmnts(2,count1) = nspec_interface(SW)
      count1 = count1 + 1
      do ispec = 1,nspec
        if ((iMPIcut_xi(1,ispec) .eqv. .true.) .and. (iMPIcut_eta(1,ispec) .eqv. .true.)) then
          neighbors_elmnts(1,count2) = ispec
          neighbors_elmnts(2,count2) = 2
          neighbors_elmnts(3,count2) = ibool(1,1,1,ispec)
          neighbors_elmnts(4,count2) = ibool(1,1,NGLLZ_M,ispec)
          neighbors_elmnts(5,count2) = -1
          neighbors_elmnts(6,count2) = -1
          count2 = count2 + 1
        endif
      enddo
    endif

    ! write my_nb_interfaces
    dset_name = "my_nb_interfaces"
    call h5_write_dataset_collect_hyperslab(dset_name, num_neighbors_elmnts, &
                                                 (/0,sum(offset_nb_interfaces(0:myrank-1))/), .true.)

    ! write my_interfaces
    dset_name = "my_interfaces"
    call h5_write_dataset_collect_hyperslab(dset_name, neighbors_elmnts, &
                                                 (/0, sum(offset_n_elms_interface(0:myrank-1))/), .true.)
    deallocate(num_neighbors_elmnts, neighbors_elmnts)

  else
    ! single process execution, no MPI boundaries
    num_interface_and_max = (/0, 0/)

    ! prepare dummy offset arrays for
    ! offset_nb_interfaces
    ! offset_n_elms_interface
    call gather_all_all_singlei(0,offset_nb_interfaces,NPROC)
    call gather_all_all_singlei(0,offset_n_elms_interface,NPROC)

    call h5_close_file_p()

    if (myrank == 0) then
      call h5_open_file(name_database_hdf5)
      dset_name = "offset_n_elms_interface"
      call h5_write_dataset_no_group(dset_name, offset_n_elms_interface)
      dset_name = "offset_nb_interfaces"
      call h5_write_dataset_no_group(dset_name, offset_nb_interfaces)
      call h5_close_file()
    endif

    call h5_open_file_p_collect(name_database_hdf5)

  endif

  dset_name = "my_ninterface_and_max"
  call h5_write_dataset_collect_hyperslab(dset_name, num_interface_and_max, (/myrank*2/), .true.)

  call h5_close_file_p()
  call h5_finalize()

  call synchronize_all()

  ! CUBIT output
  if (SAVE_MESH_AS_CUBIT) then
    ! only for single process at the moment
    if (NPROC_XI /= 1 .and. NPROC_ETA /= 1) then
      print *,'Error: SAVE_MESH_AS_CUBIT output requires NPROC_XI == NPROC_ETA == 1'
      print *,'       using NPROC_XI = ',NPROC_XI,' and NPROC_ETA = ',NPROC_ETA
      print *,'Please update your Mesh_Par_file and re-run the mesher...'
      stop 'Invalid NPROC_XI and/or NPROC_ETA for SAVE_MESH_AS_CUBIT output'
    endif

    !! VM VM add outputs as CUBIT
    call save_output_mesh_files_as_cubit(nspec,nglob, &
                                         nodes_coords, ispec_material_id, &
                                         nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                         ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                         nspec_CPML_total)

    ! output for AxiSEM coupling
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      call save_output_mesh_files_for_coupled_model(nspec, &
                                                    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                                    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                                                    xstore,ystore,zstore)
    endif
  endif

  deallocate(material_index)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done mesh files'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  ! no HDF5 compilation support
  ! to avoid compiler warnings
  logical :: l_dummy
  double precision :: d_dummy
  integer :: i_dummy

  d_dummy = nodes_coords(1,1)
  i_dummy = ispec_material_id(1)
  l_dummy = (iMPIcut_xi(1,1) .eqv. iMPIcut_eta(1,1))
  l_dummy = (ibelm_xmin(1) == ibelm_xmax(1))
  l_dummy = (ibelm_ymin(1) == ibelm_ymax(1))
  l_dummy = (ibelm_bottom(1) == ibelm_top(1))
  l_dummy = (nspec2D_xmin == nspec2D_xmax)
  l_dummy = (nspec2D_ymin == nspec2D_ymax)

  ! user output
  print *
  print *, "Error: HDF5 routine save_databases_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 save_databases_hdf5(): called without compilation support'

#endif

  end subroutine save_databases_hdf5

