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
!


  subroutine read_partition_files_hdf5

! reads in Databases.h5 file from each process
  use phdf5_utils
  use generate_databases_par
  use constants, only: IIN_database_hdf5

  implicit none
  integer :: num_xmin,num_xmax,num_ymin,num_ymax,num_top,num_bottom,num
  integer :: num_cpml
  integer :: num_moho
  integer :: i,j,count,error


  ! group names
  character(len=14) :: material_gname = "material_props"
  character(len=5)  :: gname_proc_head = "proc_"
  character(len=11) :: gname_proc
  character(len=10) :: tempstr
  character(len=30) :: dsetname, attrname
  integer, dimension(1) :: attr_data
  integer, dimension(2) :: attr_data2

  ! for boundary read
  integer, dimension(6) :: dset_n_bound
  integer, dimension(:,:), allocatable :: dset_alloc

  ! for interface
  integer, dimension(:), allocatable :: dset_alloc_1d

  ! offset arrays
  integer, dimension(0:NPROC-1) :: offset_nnodes
  integer, dimension(0:NPROC-1) :: offset_nelems
  integer, dimension(0:NPROC-1) :: offset_nelems_cpml
  integer, dimension(0:NPROC-1) :: offset_n_elms_bounds
  integer, dimension(0:NPROC-1) :: offset_n_elms_interface
  integer, dimension(0:NPROC-1) :: offset_nb_interfaces

  ! mpi variables
  integer :: info, comm

  ! hdf5 utility
  type(h5io) :: h5
  h5 = h5io()

! read databases about external mesh simulation
  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)
!
  ! opens Database file
  prname = "/Database.h5"
  IIN_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//prname

  call h5_init(h5, IIN_database_hdf5)

  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)
  !call h5_open_file_p(h5)
  call h5_open_file_p_collect(h5)

! read physical properties of the materials at first
! as this group is only out side of prop_n groups
! added poroelastic properties and filled with 0 the last 10 entries for elastic/acoustic
!
  ! open material parameters' data group
  !call h5_open_group_p(h5, material_gname)
  call h5_open_group(h5, material_gname)

  ! open dataset mat_prop
  ! read attribute nmat_ext_mesh
  dsetname = "mat_prop"
  attrname = "count_def_mat"
  call h5_read_attribute_p(h5,attrname,dsetname,attr_data)
  nmat_ext_mesh = attr_data(1)

  ! read data mat_prop
  ! allocate array for nodes coords
  allocate(mat_prop(17,nmat_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 585')
  if (error /= 0) stop 'Error allocating array mat_prop'
  mat_prop(:,:) = 0.d0
  call h5_read_dataset_p_2d_d(h5, dsetname, mat_prop)

  if (myrank == 0) then
    write(IMAIN,*) 'defined materials    : ',nmat_ext_mesh
  endif

  ! open dataset undef_mat_prop
  ! read attribute  nundefMat_ext_mesh
  dsetname = "undef_mat_prop"
  attrname = "count_undef_mat"
  call h5_read_attribute_p(h5,attrname,dsetname,attr_data)
  nundefMat_ext_mesh = attr_data(1)
  !print *, "num undef mat ext", nundefMat_ext_mesh

  if (nundefMat_ext_mesh /= 0) then ! error when trying to read 0 row Dataset in H5
    ! read undef_mat_prop
    ! allocate udef_mat_prop space
    allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array 585')
    if (error /= 0) stop 'Error allocating array mat_prop'
    call h5_read_dataset_p_2d_c(h5, dsetname, undef_mat_prop)
  endif

  call h5_close_group(h5)

  if (myrank == 0) then
    write(IMAIN,*) 'undefined materials  : ',nundefMat_ext_mesh
  endif


!
! read offset information
  dsetname = "offset_nnodes"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_nnodes, (/0/), .true.)
  dsetname = "offset_nelems"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_nelems, (/0/), .true.)
  dsetname = "offset_nelems_cpml"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_nelems_cpml, (/0/), .true.)
  dsetname = "offset_n_elms_bounds"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_n_elms_bounds, (/0/), .true.)
  dsetname = "offset_n_elms_interface"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_n_elms_interface, (/0/), .true.)
  dsetname = "offset_nb_interfaces"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, offset_nb_interfaces, (/0/), .true.)

!
!
! start reading processor dependent data
!

  !
  ! read nnodes_ext_mesh and nodes_coords_ext_mesh
  !

  dsetname = "nnodes_loc"
  call h5_read_dataset_scalar_i_collect_hyperslab(h5, dsetname, nnodes_ext_mesh,(/myrank/),.true.)
  n_control_node = nnodes_ext_mesh

  ! open and read dataset nodes_coords
  dsetname = "nodes_coords"
  ! allocate array for nodes coords
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 584')
  if (error /= 0) stop 'Error allocating array nodes_coords_ext_mesh'
  call h5_read_dataset_2d_d_collect_hyperslab(h5, dsetname, nodes_coords_ext_mesh, (/0,sum(offset_nnodes(0:myrank-1))/),.true.)

  call sum_all_i(nnodes_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'external mesh points : ',num
  endif

  !
  ! read element indexing, nelmnts_ext_mesh and mat_ext_mesh
  !

  ! open dataset elm_conn
  ! read attribute nspec_local == nelmnts_ext_mesh
  dsetname = "nspec_local"
  call h5_read_dataset_scalar_i_collect_hyperslab(h5, dsetname, nelmnts_ext_mesh, (/myrank/), .true.)

  !nelmnts_ext_mesh = attr_data(1)

  ! memory allocation
  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 587')
  if (error /= 0) stop 'Error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 588')
  if (error /= 0) stop 'Error allocating array mat_ext_mesh'

  ! read elmnts_ext_mesh == elm_conn
  dsetname = "elm_conn"
  call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, elmnts_ext_mesh, (/0,sum(offset_nelems(0:myrank-1))/),.true.)

  ! open and read dataset mat_mesh == mat_ext_mesh
  dsetname = "mat_mesh"
  call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, mat_ext_mesh, (/0,sum(offset_nelems(0:myrank-1))/),.true.)

  NSPEC_AB = nelmnts_ext_mesh

  call sum_all_i(nspec_ab,num)
  if (myrank == 0) then
    write(IMAIN,*) 'total number of spectral elements: ',num
  endif

!
! reads absorbing/free-surface boundaries
!

! open and read n_elms_on_bouund (nspec2D_xmin, xmax, ymin, ymax, button_ext, top_ext)
  dsetname = "n_elms_on_bound"
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, dset_n_bound, (/6*myrank/), .true.)

  nspec2D_xmin       = dset_n_bound(1)
  nspec2D_xmax       = dset_n_bound(2)
  nspec2D_ymin       = dset_n_bound(3)
  nspec2D_ymax       = dset_n_bound(4)
  nspec2D_bottom_ext = dset_n_bound(5)
  nspec2D_top_ext    = dset_n_bound(6)
  NSPEC2D_BOTTOM     = nspec2D_bottom_ext
  NSPEC2D_TOP        = nspec2D_top_ext

! memory arrocation
  allocate(ibelm_xmin(nspec2D_xmin),nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 589')
  if (error /= 0) stop 'Error allocating array ibelm_xmin etc.'
  allocate(ibelm_xmax(nspec2D_xmax),nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 590')
  if (error /= 0) stop 'Error allocating array ibelm_xmax etc.'
  allocate(ibelm_ymin(nspec2D_ymin),nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 591')
  if (error /= 0) stop 'Error allocating array ibelm_ymin'
  allocate(ibelm_ymax(nspec2D_ymax),nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 592')
  if (error /= 0) stop 'Error allocating array ibelm_ymax etc.'
  allocate(ibelm_bottom(nspec2D_bottom_ext),nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom_ext),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 593')
  if (error /= 0) stop 'Error allocating array ibelm_bottom etc.'
  allocate(ibelm_top(nspec2D_top_ext),nodes_ibelm_top(NGNOD2D,nspec2D_top_ext),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 594')
  if (error /= 0) stop 'Error allocating array ibelm_top etc.'

! open and read database glob2loc_elms
  dsetname = "glob2loc_elms"
  allocate(dset_alloc(1+NGNOD2D,sum(dset_n_bound(:))),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array glob2loc_elms')
  if (error /= 0) stop 'Error allocating array glob2loc_elms'
  call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, dset_alloc, &
            (/0,sum(offset_n_elms_bounds(0:myrank-1))/),.true.)

  ibelm_xmin         = dset_alloc(1          ,                           1:dset_n_bound(1))
  ibelm_xmax         = dset_alloc(1          ,dset_n_bound(1)       +1:sum(dset_n_bound(1:2)))
  ibelm_ymin         = dset_alloc(1          ,sum(dset_n_bound(1:2))+1:sum(dset_n_bound(1:3)))
  ibelm_ymax         = dset_alloc(1          ,sum(dset_n_bound(1:3))+1:sum(dset_n_bound(1:4)))
  ibelm_bottom       = dset_alloc(1          ,sum(dset_n_bound(1:4))+1:sum(dset_n_bound(1:5)))
  ibelm_top          = dset_alloc(1          ,sum(dset_n_bound(1:5))+1:sum(dset_n_bound(1:6)))
  nodes_ibelm_xmin   = dset_alloc(2:NGNOD2D+1,                           1:dset_n_bound(1))
  nodes_ibelm_xmax   = dset_alloc(2:NGNOD2D+1,dset_n_bound(1)       +1:sum(dset_n_bound(1:2)))
  nodes_ibelm_ymin   = dset_alloc(2:NGNOD2D+1,sum(dset_n_bound(1:2))+1:sum(dset_n_bound(1:3)))
  nodes_ibelm_ymax   = dset_alloc(2:NGNOD2D+1,sum(dset_n_bound(1:3))+1:sum(dset_n_bound(1:4)))
  nodes_ibelm_bottom = dset_alloc(2:NGNOD2D+1,sum(dset_n_bound(1:4))+1:sum(dset_n_bound(1:5)))
  nodes_ibelm_top    = dset_alloc(2:NGNOD2D+1,sum(dset_n_bound(1:5))+1:sum(dset_n_bound(1:6)))
  deallocate(dset_alloc)
  call sum_all_i(nspec2D_xmin,num_xmin)
  call sum_all_i(nspec2D_xmax,num_xmax)
  call sum_all_i(nspec2D_ymin,num_ymin)
  call sum_all_i(nspec2D_ymax,num_ymax)
  call sum_all_i(nspec2D_top_ext,num_top)
  call sum_all_i(nspec2D_bottom_ext,num_bottom)

  if (myrank == 0) then
    write(IMAIN,*) 'absorbing boundaries: '
    write(IMAIN,*) '  xmin,xmax : ',num_xmin,num_xmax
    write(IMAIN,*) '  ymin,ymax : ',num_ymin,num_ymax
    write(IMAIN,*) '  bottom,top: ',num_bottom,num_top
    write(IMAIN,*)
    call flush_IMAIN()
  endif

!
! read cpml conditions
!

  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml     = 0

  ! open and read nspec_cpml_tot = nspec_cpml_globloc(1), npsec_cpml = nspec_cpml_globloc(2)
  dsetname = "nspec_cpml_globloc"
  allocate(dset_alloc_1d(2),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array nspec_cpml_globloc')
  if (error /= 0) stop 'Error allocating array nspec_cpml_globloc'
  call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, dset_alloc_1d, (/2*myrank/), .true.)

  nspec_cpml_tot = dset_alloc_1d(1)
  nspec_cpml     = dset_alloc_1d(2)

  deallocate(dset_alloc_1d)

  ! array is_CPML gets always allocated to be usable in "if" checks
  allocate(is_CPML(NSPEC_AB),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 597')
  if (error /= 0) stop 'Error allocating array is_CPML'
  is_CPML(:) = .false.

  if (nspec_cpml_tot > 0) then
    if (myrank == 0) then
      write(IMAIN,*) '  number of C-PML spectral elements in this partition: ',nspec_cpml
      call flush_IMAIN()
    endif

    ! PML
    if (myrank == 0) then
      write(IMAIN,*) '  total number of C-PML elements in the global mesh: ',nspec_cpml_tot
      call flush_IMAIN()
    endif

    call sum_all_i(nspec_cpml,num_cpml)
    ! checks that the sum of C-PML elements over all partitions is correct
    if (myrank == 0 .and. nspec_cpml_tot /= num_cpml) stop 'Error while summing C-PML elements over all partitions'

    ! reads C-PML regions and C-PML spectral elements global indexing
    allocate(CPML_to_spec(nspec_cpml),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array 595')
    if (error /= 0) stop 'Error allocating array CPML_to_spec'
    allocate(CPML_regions(nspec_cpml),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array 596')
    if (error /= 0) stop 'Error allocating array CPML_regions'

    ! read elemnts_cpml and copy array to CPML_to_spec (1) and CPML_region (2)
    dsetname = "elements_cpml"
    allocate(dset_alloc(2,nspec_cpml),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array glob2loc_elms')
    if (error /= 0) stop 'Error allocating array glob2loc_elms'
    call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, dset_alloc, (/0,sum(offset_nelems_cpml(0:myrank-1))/), .true.)

    CPML_to_spec = dset_alloc(1,:)
    CPML_regions = dset_alloc(2,:)

    deallocate(dset_alloc)

    ! open and read if_cpml
    dsetname = "if_cpml"
    allocate(dset_alloc_1d(NSPEC_AB),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array if_cpml')
    if (error /= 0) stop 'Error allocating array if_cpml'
    call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, dset_alloc_1d, (/sum(offset_nelems(0:myrank-1))/), .true.)

    ! convert integer array to logical array
    do i = 1, NSPEC_AB
      if (dset_alloc_1d(i) /= 0) then
        is_CPML(i) = .true.
      endif
    enddo

    deallocate(dset_alloc_1d)

  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

!
! read interfaces
!
  if (NPROC > 1) then
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    ! from dsetname my_ninterface_and_max in hdf5
    dsetname = "my_ninterface_and_max"
    allocate(dset_alloc_1d(2),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array dset_dim_alloc_1d')
    if (error /= 0) stop 'Error allocating array dset_dim_alloc_1d'
    call h5_read_dataset_1d_i_collect_hyperslab(h5, dsetname, dset_alloc_1d, (/myrank*2/), .true.)


    ! read num_interfaces_ext_mesh
    num_interfaces_ext_mesh     = dset_alloc_1d(1)
    ! read max_interface_size_ext_mesh
    max_interface_size_ext_mesh = dset_alloc_1d(2)

    deallocate(dset_alloc_1d)

  else
    num_interfaces_ext_mesh = 0
    max_interface_size_ext_mesh = 0
  endif

  ! allocates interfaces
  allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 598')
  if (error /= 0) stop 'Error allocating array my_neighbors_ext_mesh'
  allocate(my_nelmnts_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 599')
  if (error /= 0) stop 'Error allocating array my_nelmnts_neighbors_ext_mesh'
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 600')
  if (error /= 0) stop 'Error allocating array my_interfaces_ext_mesh'
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 601')
  if (error /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh'
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=error)
  if (error /= 0) call exit_MPI_without_rank('error allocating array 602')
  if (error /= 0) stop 'Error allocating array nibool_interfaces_ext_mesh'

  ! read my_neighbors_ext_mesh and my_nelmnts_neighbors_ext_mesh
  if (NPROC > 1) then
    dsetname = "my_nb_interfaces"
    allocate(dset_alloc(2, num_interfaces_ext_mesh),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array my_nb_interfaces')
    if (error /= 0) stop 'Error allocating array my_nb_interfaces'

    call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, dset_alloc, &
                                                 (/0,sum(offset_nb_interfaces(0:myrank-1))/), .true.)


    my_neighbors_ext_mesh         = dset_alloc(1,:)
    my_nelmnts_neighbors_ext_mesh = dset_alloc(2,:)

    deallocate(dset_alloc)

    ! read my_interfaces_ext_mesh
    dsetname = "my_interfaces"
    allocate(dset_alloc(6, sum(my_nelmnts_neighbors_ext_mesh(:))),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array my_nb_interfaces')
    if (error /= 0) stop 'Error allocating array my_nb_interfaces'
    call h5_read_dataset_2d_i_collect_hyperslab(h5, dsetname, dset_alloc, &
                                                 (/0, sum(offset_n_elms_interface(0:myrank-1))/), .true.)
  endif

  ! loops over MPI interfaces with other partitions
  count = 1
  do num_interface = 1, num_interfaces_ext_mesh
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    ! loops over interface elements
    do i = 1, my_nelmnts_neighbors_ext_mesh(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5)...
      !
      ! interface types:
      !     1  -  corner point only
      !     2  -  element edge
      !     4  -  element face
      my_interfaces_ext_mesh(1,i,num_interface) = dset_alloc(1,count)
      my_interfaces_ext_mesh(2,i,num_interface) = dset_alloc(2,count)
      my_interfaces_ext_mesh(3,i,num_interface) = dset_alloc(3,count)
      my_interfaces_ext_mesh(4,i,num_interface) = dset_alloc(4,count)
      my_interfaces_ext_mesh(5,i,num_interface) = dset_alloc(5,count)
      my_interfaces_ext_mesh(6,i,num_interface) = dset_alloc(6,count)

      count = count + 1
    enddo
  enddo

  if (NPROC > 1) deallocate(dset_alloc)

  call sum_all_i(num_interfaces_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'number of MPI partition interfaces: ',num
  endif

  ! optional moho
  if (SAVE_MOHO_MESH) then
    ! checks if additional line exists
    ! open dataset moho_elms
    dsetname = "moho_elms"
    ! read attribute loc_moho == boundary_number, nspec2D_moho_ext
    attrname = "loc_moho"
    call h5_read_attribute_p(h5,attrname,dsetname,attr_data2)

    boundary_number  = attr_data2(1)
    nspec2D_moho_ext = attr_data2(2)

    ! checks total number of elements
    call sum_all_i(nspec2D_moho_ext,num_moho)
    if (num_moho == 0) call exit_mpi(myrank,'Error no moho mesh in database')

    ! read moho_elms(1,:) == ibelm_moho and moho_elms(2:,:) == nodes_ibelm_moho()
    dsetname = "moho_elms"
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array 603')
    if (error /= 0) stop 'Error allocating array ibelm_moho etc.'
    allocate(dset_alloc(1+NGNOD2D, nspec2D_moho_ext),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array moho_elms')
    if (error /= 0) stop 'Error allocating array moho_elms'
    call h5_read_dataset_p_2d_i(h5, dsetname, dset_alloc)

    ibelm_moho       = dset_alloc(1,:)
    nodes_ibelm_moho = dset_alloc(2:,:)

    deallocate(dset_alloc)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'moho surfaces: ',num_moho
    endif

  else
    ! allocate dummy array
    nspec2D_moho_ext = 0
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=error)
    if (error /= 0) call exit_MPI_without_rank('error allocating array 604')
    if (error /= 0) stop 'Error allocating dumy array ibelm_moho etc.'
  endif

  ! close hdf5
  call h5_close_file_p(h5)
  call h5_destructor(h5)

end subroutine read_partition_files_hdf5
