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

!==============================================================================
!> Reads in Database.bp file
  subroutine read_partition_files_adios()

  use generate_databases_par, only: MAX_STRING_LEN,LOCAL_PATH,IMAIN,myrank, &
    NDIM,NGLLX,NGNOD,NGNOD2D,NPROC,NSPEC_AB,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    nspec2D_top_ext,nspec2D_bottom_ext,nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    nodes_coords_ext_mesh,mat_prop,undef_mat_prop,nundefMat_ext_mesh, &
    elmnts_ext_mesh,nelmnts_ext_mesh,mat_ext_mesh, &
    nmat_ext_mesh,nnodes_ext_mesh, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax,nodes_ibelm_bottom,nodes_ibelm_top

  ! MPI interfaces
  use generate_databases_par, only: num_interfaces_ext_mesh,my_neighbors_ext_mesh,my_interfaces_ext_mesh, &
    ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
    my_nelmnts_neighbors_ext_mesh,max_interface_size_ext_mesh

  ! CPML
  use generate_databases_par, only: is_CPML,CPML_to_spec,CPML_regions,nspec_cpml,nspec_cpml_tot

  use adios_read_mod
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  implicit none

  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8) :: local_dim_nodes_coords,    local_dim_matpropl, &
    local_dim_material_index,  local_dim_elmnts_mesh, &
    local_dim_ibelm_xmin,      local_dim_nodes_ibelm_xmin, &
    local_dim_ibelm_xmax,      local_dim_nodes_ibelm_xmax, &
    local_dim_ibelm_ymin,      local_dim_nodes_ibelm_ymin, &
    local_dim_ibelm_ymax,      local_dim_nodes_ibelm_ymax, &
    local_dim_ibelm_bottom,    local_dim_nodes_ibelm_bottom, &
    local_dim_ibelm_top,       local_dim_nodes_ibelm_top, &
    local_dim_neighbors_mesh,  local_dim_num_elmnts_mesh, &
    local_dim_interfaces_mesh
  integer(kind=8) :: local_dim_cpml_to_spec,local_dim_cpml_regions,local_dim_is_cpml
  integer(kind=8) :: local_dim
  integer :: comm

  !statistics
  integer :: num
  !integer :: num_xmin, num_xmax, num_ymin, num_ymax, num_top, num_bottom, num_spec, num_int
  integer,parameter :: NUNDEFMAT_MAX = 10
  character(len=6*NUNDEFMAT_MAX*MAX_STRING_LEN) :: undef_matpropl
  integer :: icount,j,is,ie
  integer :: ier

  undef_matpropl(:) = ''
  sel_num = 0

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/Database.bp"

  ! gets MPI communicator
  comm = comm_adios

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)
  call adios_read_open_file (handle, database_name, 0, comm, ier)
  if (ier /= 0) call abort_mpi()

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  !------------------------------------------------------------------'
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "nglob", 0, 1, nnodes_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "nspec", 0, 1, nelmnts_ext_mesh, ier)

  ! read physical properties of the materials
  ! added poroelastic properties and filled with 0 the last 10 entries
  ! for elastic/acoustic
  call adios_schedule_read(handle, sel, "nmaterials", 0, 1, nmat_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "nundef_materials", 0, 1, nundefMat_ext_mesh, ier)

  call adios_schedule_read(handle, sel, "nspec2d_xmin", 0, 1, nspec2D_xmin, ier)
  call adios_schedule_read(handle, sel, "nspec2d_xmax", 0, 1, nspec2D_xmax, ier)
  call adios_schedule_read(handle, sel, "nspec2d_ymin", 0, 1, nspec2D_ymin, ier)
  call adios_schedule_read(handle, sel, "nspec2d_ymax", 0, 1, nspec2D_ymax, ier)
  call adios_schedule_read(handle, sel, "nspec2d_bottom", 0, 1, nspec2D_bottom_ext, ier)
  call adios_schedule_read(handle, sel, "nspec2d_top", 0, 1, nspec2D_top_ext, ier)

  ! MPI interfaces between different partitions
  num_interfaces_ext_mesh = 0
  max_interface_size_ext_mesh = 0
  if (NPROC > 1) then
    ! format: #number_of_MPI_interfaces
    !         #maximum_number_of_elements_on_each_interface
    call adios_schedule_read(handle, sel, "nb_interfaces", 0, 1, num_interfaces_ext_mesh, ier)
    call adios_schedule_read(handle, sel, "nspec_interfaces_max", 0, 1, max_interface_size_ext_mesh, ier)
  endif

  ! CPML
  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml = 0
  call adios_schedule_read(handle, sel, "nspec_cpml_total", 0, 1, nspec_cpml_tot, ier)
  call adios_schedule_read(handle, sel, "nspec_cpml", 0, 1, nspec_cpml, ier)

  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  ! sets nspec parameters
  NSPEC_AB = nelmnts_ext_mesh
  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext

  ! user output
  call sum_all_i(nnodes_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'external mesh points : ',num
    write(IMAIN,*) 'defined materials    : ',nmat_ext_mesh
    write(IMAIN,*) 'undefined materials  : ',nundefMat_ext_mesh
  endif
  call synchronize_all()

  call sum_all_i(nspec_ab,num)
  if (myrank == 0) then
    write(IMAIN,*) 'total number of spectral elements: ',num
  endif
  call synchronize_all()

  call sum_all_i(num_interfaces_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'number of MPI partition interfaces: ',num
  endif
  call synchronize_all()

  !---------------------------------------------.
  ! Allocate arrays with previously read values |
  !---------------------------------------------'
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 668')
  if (ier /= 0) stop 'error allocating array nodes_coords_ext_mesh'
  nodes_coords_ext_mesh(:,:) = 0.d0

  allocate(mat_prop(17,nmat_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 669')
  if (ier /= 0) stop 'error allocating array mat_prop'
  mat_prop(:,:) = 0.d0

  allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 670')
  if (ier /= 0) stop 'error allocating array undef_mat_prop'

  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 671')
  if (ier /= 0) stop 'error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 672')
  if (ier /= 0) stop 'error allocating array mat_ext_mesh'
  elmnts_ext_mesh(:,:) = 0; mat_ext_mesh(:,:) = 0

  ! boundaries
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

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 685')
    if (ier /= 0) stop 'error allocating array my_neighbors_ext_mesh'
    allocate(my_nelmnts_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 686')
    if (ier /= 0) stop 'error allocating array my_nelmnts_neighbors_ext_mesh'
    allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 687')
    if (ier /= 0) stop 'error allocating array my_interfaces_ext_mesh'
    allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 688')
    if (ier /= 0) stop 'error allocating array ibool_interfaces_ext_mesh'
    allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 689')
    if (ier /= 0) stop 'error allocating array nibool_interfaces_ext_mesh'
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
  my_interfaces_ext_mesh(:,:,:) = 0; ibool_interfaces_ext_mesh(:,:) = 0; nibool_interfaces_ext_mesh(:) = 0

  ! CPML
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 597')
  if (ier /= 0) stop 'Error allocating array is_CPML'
  is_CPML(:) = .false.

  if (nspec_cpml > 0) then
    ! reads C-PML regions and C-PML spectral elements global indexing
    allocate(CPML_to_spec(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 595')
    if (ier /= 0) stop 'Error allocating array CPML_to_spec'
    CPML_to_spec(:) = 0
    allocate(CPML_regions(nspec_cpml),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 596')
    if (ier /= 0) stop 'Error allocating array CPML_regions'
    CPML_regions(:) = 0
  else
    ! dummy allocations
    allocate(CPML_to_spec(1),CPML_regions(1))
  endif

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  call adios_get_scalar(handle, "nodes_coords/local_dim", local_dim_nodes_coords, ier)

  if (nmat_ext_mesh /= 0) then
    call adios_get_scalar(handle, "matpropl/local_dim", local_dim_matpropl, ier)
  endif

  call adios_get_scalar(handle, "material_index/local_dim", local_dim_material_index, ier)
  call adios_get_scalar(handle, "elmnts_mesh/local_dim", local_dim_elmnts_mesh, ier)

  call adios_get_scalar(handle, "ibelm_xmin/local_dim", local_dim_ibelm_xmin, ier)
  call adios_get_scalar(handle, "ibelm_xmax/local_dim", local_dim_ibelm_xmax, ier)
  call adios_get_scalar(handle, "ibelm_ymin/local_dim", local_dim_ibelm_ymin, ier)
  call adios_get_scalar(handle, "ibelm_ymax/local_dim", local_dim_ibelm_ymax, ier)
  call adios_get_scalar(handle, "ibelm_bottom/local_dim", local_dim_ibelm_bottom, ier)
  call adios_get_scalar(handle, "ibelm_top/local_dim", local_dim_ibelm_top, ier)

  call adios_get_scalar(handle, "nodes_ibelm_xmin/local_dim", local_dim_nodes_ibelm_xmin, ier)
  call adios_get_scalar(handle, "nodes_ibelm_xmax/local_dim", local_dim_nodes_ibelm_xmax, ier)
  call adios_get_scalar(handle, "nodes_ibelm_ymin/local_dim", local_dim_nodes_ibelm_ymin, ier)
  call adios_get_scalar(handle, "nodes_ibelm_ymax/local_dim", local_dim_nodes_ibelm_ymax, ier)
  call adios_get_scalar(handle, "nodes_ibelm_bottom/local_dim", local_dim_nodes_ibelm_bottom, ier)
  call adios_get_scalar(handle, "nodes_ibelm_top/local_dim", local_dim_nodes_ibelm_top, ier)

  if (num_interfaces_ext_mesh > 0) then
    call adios_get_scalar(handle, "neighbors_mesh/local_dim", local_dim_neighbors_mesh, ier)
    call adios_get_scalar(handle, "num_elmnts_mesh/local_dim", local_dim_num_elmnts_mesh, ier)
    call adios_get_scalar(handle, "interfaces_mesh/local_dim", local_dim_interfaces_mesh, ier)
  endif

  if (nspec_cpml > 0) then
    call adios_get_scalar(handle, "CPML_to_spec/local_dim", local_dim_cpml_to_spec, ier)
    call adios_get_scalar(handle, "CPML_regions/local_dim", local_dim_cpml_regions, ier)
    call adios_get_scalar(handle, "is_CPML/local_dim", local_dim_is_cpml, ier)
  endif

  !---------------------------.
  ! Read arrays from Database |
  !---------------------------'
  start(1) = local_dim_nodes_coords * myrank
  count(1) = NDIM * nnodes_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  ! Remember: "nodes_coords" has been transposed in meshfem3D
  call adios_schedule_read(handle, sel, "nodes_coords/array", 0, 1, nodes_coords_ext_mesh, ier)

  if (nmat_ext_mesh /= 0) then
    start(1) = local_dim_matpropl * myrank
    count(1) = 17 * nmat_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    ! Remember: "matpropl" has been transposed in meshfem3D
    call adios_schedule_read(handle, sel, "matpropl/array", 0, 1, mat_prop, ier)
  endif

  if (nundefMat_ext_mesh /= 0) then
    ! checks for sufficient buffer length
    if (nundefMat_ext_mesh > NUNDEFMAT_MAX) stop 'Error number of undefined materials exceeds maximum NUNDEFMAT_MAX'
    ! note: arrays not working yet for strings, uses "local" variable
    local_dim = 6 * nundefMat_ext_mesh * MAX_STRING_LEN
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(handle, sel, "undef_mat_prop", 0, 1, undef_matpropl(1:local_dim), ier)
  endif

  start(1) = local_dim_material_index * myrank
  count(1) = 2 * nelmnts_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "material_index/array", 0, 1, mat_ext_mesh, ier)

  start(1) = local_dim_elmnts_mesh * myrank
  count(1) = NGNOD * nelmnts_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "elmnts_mesh/array", 0, 1, elmnts_ext_mesh, ier)

  start(1) = local_dim_ibelm_xmin * myrank
  count(1) = nspec2D_xmin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_xmin/array", 0, 1, ibelm_xmin, ier)

  start(1) = local_dim_ibelm_xmax * myrank
  count(1) = nspec2D_xmax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_xmax/array", 0, 1, ibelm_xmax, ier)

  start(1) = local_dim_ibelm_ymin * myrank
  count(1) = nspec2D_ymin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_ymin/array", 0, 1, ibelm_ymin, ier)

  start(1) = local_dim_ibelm_ymax * myrank
  count(1) = nspec2D_ymax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_ymax/array", 0, 1, ibelm_ymax, ier)

  start(1) = local_dim_ibelm_bottom * myrank
  count(1) = nspec2D_bottom_ext
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_bottom/array", 0, 1, ibelm_bottom, ier)

  start(1) = local_dim_ibelm_top * myrank
  count(1) = nspec2D_top_ext
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_top/array", 0, 1, ibelm_top, ier)

  ! needs surface boundary points
  if (nspec2d_xmin > 0) then
    start(1) = local_dim_nodes_ibelm_xmin * myrank
    count(1) = nspec2D_xmin * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_xmin/array", 0, 1, nodes_ibelm_xmin, ier)
  endif

  if (nspec2d_xmax > 0) then
    start(1) = local_dim_nodes_ibelm_xmax * myrank
    count(1) = nspec2D_xmax * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_xmax/array", 0, 1, nodes_ibelm_xmax, ier)
  endif

  if (nspec2d_ymin > 0) then
    start(1) = local_dim_nodes_ibelm_ymin * myrank
    count(1) = nspec2D_ymin * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_ymin/array", 0, 1, nodes_ibelm_ymin, ier)
  endif

  if (nspec2d_ymax > 0) then
    start(1) = local_dim_nodes_ibelm_ymax * myrank
    count(1) = nspec2D_ymax * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_ymax/array", 0, 1, nodes_ibelm_ymax, ier)
  endif

  if (nspec2d_bottom > 0) then
    start(1) = local_dim_nodes_ibelm_bottom * myrank
    count(1) = nspec2D_bottom_ext * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_bottom/array", 0, 1, nodes_ibelm_bottom, ier)
  endif

  if (nspec2d_top > 0) then
    start(1) = local_dim_nodes_ibelm_top * myrank
    count(1) = nspec2D_top_ext * NGNOD2D
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "nodes_ibelm_top/array", 0, 1, nodes_ibelm_top, ier)
  endif

  ! CPML
  if (nspec_cpml > 0) then
    start(1) = local_dim_cpml_to_spec * myrank
    count(1) = nspec_cpml
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "CPML_to_spec/array", 0, 1, CPML_to_spec, ier)

    start(1) = local_dim_cpml_regions * myrank
    count(1) = nspec_cpml
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "CPML_regions/array", 0, 1, CPML_regions, ier)

    start(1) = local_dim_is_cpml * myrank
    count(1) = NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "is_CPML/array", 0, 1, is_CPML, ier)
  endif

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    start(1) = local_dim_neighbors_mesh * myrank
    count(1) = num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "neighbors_mesh/array", 0, 1, my_neighbors_ext_mesh, ier)

    start(1) = local_dim_num_elmnts_mesh * myrank
    count(1) = num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "num_elmnts_mesh/array", 0, 1, my_nelmnts_neighbors_ext_mesh, ier)

    start(1) = local_dim_interfaces_mesh * myrank
    count(1) = 6 * num_interfaces_ext_mesh * max_interface_size_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(handle, sel, "interfaces_mesh/array", 0, 1, my_interfaces_ext_mesh, ier)
  endif

  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  ! fills undef_mat_prop 2D array
  undef_mat_prop(:,:) = ''
  if (nundefMat_ext_mesh > 0) then
    do icount = 1,nundefMat_ext_mesh
      do j = 1,6
        is = (icount-1)*6*MAX_STRING_LEN + (j-1)*MAX_STRING_LEN + 1
        ie = is + MAX_STRING_LEN
        undef_mat_prop(j,icount) = undef_matpropl(is:ie)
      enddo
    enddo
  endif

  ! debug
  !if (myrank == 1) print *,myrank,'undef_matpropl: ',trim(undef_matpropl)
  !do icount = 0,NPROC-1
  !  print *
  !  if (myrank == icount) then
  !    print *,myrank,'undef_mat_prop: '
  !    do j = 1,6
  !      print *,'  index ',j,': ',trim(undef_mat_prop(j,1))
  !    enddo
  !    call flush(6)
  !  endif
  !  call synchronize_all()
  !enddo

  ! TODO SAVE_MOHO_MESH w/ ADIOS (not in meshfem)

  !-----------------------------.
  ! Fetch and print come stats. |
  !-----------------------------'
  !call sum_all_i(nspec_ab,num_spec)
  !call sum_all_i(nnodes_ext_mesh,num)
  !call sum_all_i(nspec2D_xmin,num_xmin)
  !call sum_all_i(nspec2D_xmax,num_xmax)
  !call sum_all_i(nspec2D_ymin,num_ymin)
  !call sum_all_i(nspec2D_ymax,num_ymax)
  !call sum_all_i(nspec2D_top_ext,num_top)
  !call sum_all_i(nspec2D_bottom_ext,num_bottom)
  !call sum_all_i(num_interfaces_ext_mesh,num_int)
  !if (myrank == 0) then
    !write(IMAIN,*) ' total number of spectral elements: ',num_spec
    !write(IMAIN,*) '  external mesh points: ',num
    !write(IMAIN,*) '  defined materials: ',nmat_ext_mesh
    !write(IMAIN,*) '  undefined materials: ',nundefMat_ext_mesh
    !write(IMAIN,*) '  absorbing boundaries: '
    !write(IMAIN,*) '    xmin,xmax: ',num_xmin,num_xmax
    !write(IMAIN,*) '    ymin,ymax: ',num_ymin,num_ymax
    !write(IMAIN,*) '    bottom,top: ',num_bottom,num_top
    !write(IMAIN,*) '  number of MPI partition interfaces: ',num_int
  !endif
  !call synchronize_all()

end subroutine read_partition_files_adios
