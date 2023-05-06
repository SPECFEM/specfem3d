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

! for HDF5 file I/O

  subroutine read_partition_files_hdf5()

! reads in Databases.h5 file from each process

#ifdef USE_HDF5

  use generate_databases_par
  use manager_hdf5

  implicit none
  ! local parameters
  integer :: num_xmin,num_xmax,num_ymin,num_ymax,num_top,num_bottom,num
  integer :: num_cpml
  integer :: num_moho
  integer :: i,count,ier

  ! group names
  character(len=14) :: material_gname = "material_props"
  character(len=11) :: gname_proc
  character(len=64) :: group_name
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

  ! MPI variables
  integer :: info, comm

  ! read databases about external mesh simulation
  ! opens Database file
  prname = "/Database.h5"
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//prname

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading database file: ',trim(name_database_hdf5)
    write(IMAIN,*) '  using HDF5 file format'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI parameters
  call world_get_comm(comm)
  call world_get_info_null(info)

  ! initializes HDF5
  call h5_initialize()
  call h5_set_mpi_info(comm, info, myrank, NPROC)

  ! opens file
  !call h5_open_file_p(name_database_hdf5)
  call h5_open_file_p_collect(name_database_hdf5)

! read physical properties of the materials at first
! as this group is only out side of prop_n groups
! added poroelastic properties and filled with 0 the last 10 entries for elastic/acoustic

  ! open material parameters' data group
  !call h5_open_group_p(material_gname)
  call h5_open_group(material_gname)

  ! open dataset mat_prop
  ! read attribute nmat_ext_mesh
  dsetname = "mat_prop"
  attrname = "count_def_mat"
  call h5_read_attribute_p(attrname,dsetname,attr_data)
  nmat_ext_mesh = attr_data(1)

  ! read data mat_prop
  ! allocate array for nodes coords
  allocate(mat_prop(17,nmat_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 585')
  if (ier /= 0) stop 'Error allocating array mat_prop'
  mat_prop(:,:) = 0.d0

  call h5_read_dataset_p(dsetname, mat_prop)

  if (myrank == 0) then
    write(IMAIN,*) 'defined materials    : ',nmat_ext_mesh
  endif
  call synchronize_all()

  ! open dataset undef_mat_prop
  ! read attribute  nundefMat_ext_mesh
  dsetname = "undef_mat_prop"
  attrname = "count_undef_mat"
  call h5_read_attribute_p(attrname,dsetname,attr_data)

  nundefMat_ext_mesh = attr_data(1)
  !print *, "num undef mat ext", nundefMat_ext_mesh
  ! error when trying to read 0 row Dataset in H5
  if (nundefMat_ext_mesh /= 0) then
    ! read undef_mat_prop
    ! allocate udef_mat_prop space
    allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 585')
    if (ier /= 0) stop 'Error allocating array mat_prop'
    call h5_read_dataset_p(dsetname, undef_mat_prop)
  endif

  call h5_close_group()

  if (myrank == 0) then
    write(IMAIN,*) 'undefined materials  : ',nundefMat_ext_mesh
  endif
  call synchronize_all()

  ! read offset information
  dsetname = "offset_nnodes"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_nnodes, (/0/), .true.)
  dsetname = "offset_nelems"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_nelems, (/0/), .true.)
  dsetname = "offset_nelems_cpml"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_nelems_cpml, (/0/), .true.)
  dsetname = "offset_n_elms_bounds"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_n_elms_bounds, (/0/), .true.)
  dsetname = "offset_n_elms_interface"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_n_elms_interface, (/0/), .true.)
  dsetname = "offset_nb_interfaces"
  call h5_read_dataset_collect_hyperslab(dsetname, offset_nb_interfaces, (/0/), .true.)

  !
  ! start reading processor dependent data
  !

  ! read nnodes_ext_mesh and nodes_coords_ext_mesh
  dsetname = "nnodes_loc"
  call h5_read_dataset_scalar_collect_hyperslab(dsetname, nnodes_ext_mesh,(/myrank/),.true.)

  ! stores number of mesh points in this slice for output in checkmesh.xmf
  xdmf_mesh_nnodes = nnodes_ext_mesh

  ! open and read dataset nodes_coords
  dsetname = "nodes_coords"
  ! allocate array for nodes coords
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 584')
  if (ier /= 0) stop 'Error allocating array nodes_coords_ext_mesh'
  nodes_coords_ext_mesh(:,:) = 0.0

  call h5_read_dataset_collect_hyperslab(dsetname, nodes_coords_ext_mesh, (/0,sum(offset_nnodes(0:myrank-1))/),.true.)

  call sum_all_i(nnodes_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'external mesh points : ',num
  endif

  ! read element indexing, nelmnts_ext_mesh and mat_ext_mesh

  ! open dataset elm_conn
  ! read attribute nspec_local == nelmnts_ext_mesh
  dsetname = "nspec_local"
  call h5_read_dataset_scalar_collect_hyperslab(dsetname, nelmnts_ext_mesh, (/myrank/), .true.)

  !nelmnts_ext_mesh = attr_data(1)

  ! memory allocation
  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 587')
  if (ier /= 0) stop 'Error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 588')
  if (ier /= 0) stop 'Error allocating array mat_ext_mesh'
  elmnts_ext_mesh(:,:) = 0; mat_ext_mesh(:,:) = 0

  ! read elmnts_ext_mesh == elm_conn
  dsetname = "elm_conn"
  call h5_read_dataset_collect_hyperslab(dsetname, elmnts_ext_mesh, (/0,sum(offset_nelems(0:myrank-1))/),.true.)

  ! open and read dataset mat_mesh == mat_ext_mesh
  dsetname = "mat_mesh"
  call h5_read_dataset_collect_hyperslab(dsetname, mat_ext_mesh, (/0,sum(offset_nelems(0:myrank-1))/),.true.)

  NSPEC_AB = nelmnts_ext_mesh

  call sum_all_i(nspec_ab,num)
  if (myrank == 0) then
    write(IMAIN,*) 'total number of spectral elements: ',num
  endif
  call synchronize_all()

  ! Boundaries
  ! reads absorbing/free-surface boundaries
  ! open and read n_elms_on_bouund (nspec2D_xmin, xmax, ymin, ymax, button_ext, top_ext)
  dsetname = "n_elms_on_bound"
  call h5_read_dataset_collect_hyperslab(dsetname, dset_n_bound, (/6*myrank/), .true.)

  nspec2D_xmin       = dset_n_bound(1)
  nspec2D_xmax       = dset_n_bound(2)
  nspec2D_ymin       = dset_n_bound(3)
  nspec2D_ymax       = dset_n_bound(4)
  nspec2D_bottom_ext = dset_n_bound(5)
  nspec2D_top_ext    = dset_n_bound(6)
  NSPEC2D_BOTTOM     = nspec2D_bottom_ext
  NSPEC2D_TOP        = nspec2D_top_ext

  ! memory arrocation
  if (nspec2D_xmin > 0) then
    allocate(ibelm_xmin(nspec2D_xmin),nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 589')
    if (ier /= 0) stop 'Error allocating array ibelm_xmin etc.'
  else
    allocate(ibelm_xmin(1),nodes_ibelm_xmin(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_xmin(:) = 0; nodes_ibelm_xmin(:,:) = 0

  if (nspec2D_xmax > 0) then
    allocate(ibelm_xmax(nspec2D_xmax),nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 590')
    if (ier /= 0) stop 'Error allocating array ibelm_xmax etc.'
  else
    allocate(ibelm_xmax(1),nodes_ibelm_xmax(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_xmax(:) = 0; nodes_ibelm_xmax(:,:) = 0

  if (nspec2D_ymin > 0) then
    allocate(ibelm_ymin(nspec2D_ymin),nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 591')
    if (ier /= 0) stop 'Error allocating array ibelm_ymin'
  else
    allocate(ibelm_ymin(1),nodes_ibelm_ymin(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_ymin(:) = 0; nodes_ibelm_ymin(:,:) = 0

  if (nspec2D_ymax > 0) then
    allocate(ibelm_ymax(nspec2D_ymax),nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 592')
    if (ier /= 0) stop 'Error allocating array ibelm_ymax etc.'
  else
    allocate(ibelm_ymax(1),nodes_ibelm_ymax(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_ymax(:) = 0; nodes_ibelm_ymax(:,:) = 0

  if (nspec2D_bottom_ext > 0) then
    allocate(ibelm_bottom(nspec2D_bottom_ext),nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom_ext),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 593')
    if (ier /= 0) stop 'Error allocating array ibelm_bottom etc.'
  else
    allocate(ibelm_bottom(1),nodes_ibelm_bottom(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_bottom(:) = 0; nodes_ibelm_bottom(:,:) = 0

  if (nspec2D_top_ext > 0) then
    allocate(ibelm_top(nspec2D_top_ext),nodes_ibelm_top(NGNOD2D,nspec2D_top_ext),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 594')
    if (ier /= 0) stop 'Error allocating array ibelm_top etc.'
  else
    allocate(ibelm_top(1),nodes_ibelm_top(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
  endif
  ibelm_top(:) = 0; nodes_ibelm_top(:,:) = 0

  ! open and read database glob2loc_elms
  dsetname = "glob2loc_elms"
  allocate(dset_alloc(1+NGNOD2D,sum(dset_n_bound(:))),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array glob2loc_elms')
  if (ier /= 0) stop 'Error allocating array glob2loc_elms'

  call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc, &
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
  call synchronize_all()

!
! read cpml conditions
!

  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml     = 0

  ! open and read nspec_cpml_tot = nspec_cpml_globloc(1), npsec_cpml = nspec_cpml_globloc(2)
  dsetname = "nspec_cpml_globloc"
  allocate(dset_alloc_1d(2),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array nspec_cpml_globloc')
  if (ier /= 0) stop 'Error allocating array nspec_cpml_globloc'
  call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc_1d, (/2*myrank/), .true.)

  nspec_cpml_tot = dset_alloc_1d(1)
  nspec_cpml     = dset_alloc_1d(2)

  deallocate(dset_alloc_1d)

  if (myrank == 0) then
    write(IMAIN,*) '  total number of C-PML elements in the global mesh: ',nspec_cpml_tot
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! array is_CPML gets always allocated to be usable in "if" checks
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 597')
  if (ier /= 0) stop 'Error allocating array is_CPML'
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
    allocate(CPML_to_spec(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 595')
    if (ier /= 0) stop 'Error allocating array CPML_to_spec'
    CPML_to_spec(:) = 0
    allocate(CPML_regions(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 596')
    if (ier /= 0) stop 'Error allocating array CPML_regions'
    CPML_regions(:) = 0

    ! read elemnts_cpml and copy array to CPML_to_spec (1) and CPML_region (2)
    dsetname = "elements_cpml"
    allocate(dset_alloc(2,nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array glob2loc_elms')
    if (ier /= 0) stop 'Error allocating array glob2loc_elms'
    call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc, (/0,sum(offset_nelems_cpml(0:myrank-1))/), .true.)

    CPML_to_spec = dset_alloc(1,:)
    CPML_regions = dset_alloc(2,:)

    deallocate(dset_alloc)

    ! open and read if_cpml
    dsetname = "if_cpml"
    allocate(dset_alloc_1d(NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array if_cpml')
    if (ier /= 0) stop 'Error allocating array if_cpml'
    call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc_1d, (/sum(offset_nelems(0:myrank-1))/), .true.)

    ! convert integer array to logical array
    do i = 1,NSPEC_AB
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
  call synchronize_all()

  ! MPI interfaces (between different partitions)
  ! read interfaces
  if (NPROC > 1) then
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    ! from dsetname my_ninterface_and_max in hdf5
    dsetname = "my_ninterface_and_max"
    allocate(dset_alloc_1d(2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array dset_dim_alloc_1d')
    if (ier /= 0) stop 'Error allocating array dset_dim_alloc_1d'
    call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc_1d, (/myrank*2/), .true.)

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
  if (num_interfaces_ext_mesh > 0) then
    allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 598')
    if (ier /= 0) stop 'Error allocating array my_neighbors_ext_mesh'
    allocate(my_nelmnts_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 599')
    if (ier /= 0) stop 'Error allocating array my_nelmnts_neighbors_ext_mesh'
    allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 600')
    if (ier /= 0) stop 'Error allocating array my_interfaces_ext_mesh'
    allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 601')
    if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh'
    allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 602')
    if (ier /= 0) stop 'Error allocating array nibool_interfaces_ext_mesh'
  else
    ! dummy allocations
    allocate(my_neighbors_ext_mesh(1), &
             my_nelmnts_neighbors_ext_mesh(1), &
             my_interfaces_ext_mesh(1,1,1), &
             ibool_interfaces_ext_mesh(1,1), &
             nibool_interfaces_ext_mesh(1),stat=ier)
    if (ier /= 0) stop 'Error allocating neighbors arrays'
  endif
  my_neighbors_ext_mesh(:) = -1; my_nelmnts_neighbors_ext_mesh(:) = 0
  my_interfaces_ext_mesh(:,:,:) = -1; ibool_interfaces_ext_mesh(:,:) = 0; nibool_interfaces_ext_mesh(:) = 0

  ! read my_neighbors_ext_mesh and my_nelmnts_neighbors_ext_mesh
  if (NPROC > 1) then
    dsetname = "my_nb_interfaces"
    allocate(dset_alloc(2, num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array my_nb_interfaces')
    if (ier /= 0) stop 'Error allocating array my_nb_interfaces'

    call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc, &
                                                 (/0,sum(offset_nb_interfaces(0:myrank-1))/), .true.)

    my_neighbors_ext_mesh         = dset_alloc(1,:)
    my_nelmnts_neighbors_ext_mesh = dset_alloc(2,:)

    deallocate(dset_alloc)

    ! read my_interfaces_ext_mesh
    dsetname = "my_interfaces"
    allocate(dset_alloc(6, sum(my_nelmnts_neighbors_ext_mesh(:))),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array my_nb_interfaces')
    if (ier /= 0) stop 'Error allocating array my_nb_interfaces'
    call h5_read_dataset_collect_hyperslab(dsetname, dset_alloc, &
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
  call synchronize_all()

  ! from CUBIT/trelis decompose database
  ! decomposer adds additional mesh infos for Moho surface in case choosen
  !
  ! optional moho
  if (SAVE_MOHO_MESH) then
    ! checks if additional line exists
    ! open dataset moho_elms
    dsetname = "moho_elms"

    ! create group name for the database of the current proc
    write(gname_proc,"('proc',i6.6,'_')") myrank
    group_name = trim(gname_proc) // trim(dsetname)     ! proc000000_moho_elms, ..

    call h5_open_group(group_name)

    ! read attribute loc_moho == boundary_number, nspec2D_moho_ext
    attrname = "loc_moho"
    call h5_read_attribute_p(attrname,dsetname,attr_data2)

    boundary_number  = attr_data2(1)
    nspec2D_moho_ext = attr_data2(2)

    ! checks total number of elements
    call sum_all_i(nspec2D_moho_ext,num_moho)
    if (num_moho == 0) call exit_mpi(myrank,'Error no moho mesh in database')

    ! read moho_elms(1,:) == ibelm_moho and moho_elms(2:,:) == nodes_ibelm_moho()
    dsetname = "moho_elms"
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 603')
    if (ier /= 0) stop 'Error allocating array ibelm_moho etc.'
    allocate(dset_alloc(1+NGNOD2D, nspec2D_moho_ext),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array moho_elms')
    if (ier /= 0) stop 'Error allocating array moho_elms'
    call h5_read_dataset_p(dsetname, dset_alloc)

    ibelm_moho       = dset_alloc(1,:)
    nodes_ibelm_moho = dset_alloc(2:,:)

    call h5_close_group()

    deallocate(dset_alloc)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'moho surfaces: ',num_moho
    endif
    call synchronize_all()
  else
    ! allocate dummy array
    nspec2D_moho_ext = 0
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 604')
    if (ier /= 0) stop 'Error allocating dumy array ibelm_moho etc.'
  endif

  ! close hdf5
  call h5_close_file_p()
  call h5_finalize()

#else
  ! no HDF5 compilation support
  ! user output
  print *
  print *, "Error: HDF5 routine read_partition_files_hdf5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 read_partition_files_hdf5(): called without compilation support'

#endif

end subroutine read_partition_files_hdf5
