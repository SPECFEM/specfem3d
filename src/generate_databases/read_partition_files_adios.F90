!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
!

!==============================================================================
!> Reads in Database.bp file
subroutine read_partition_files_adios()

  use mpi
  use adios_read_mod
  use generate_databases_par

  implicit none

  character(len=256) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count

  integer :: local_dim_nodes_coords,    local_dim_matpropl,           &
             local_dim_material_index,  local_dim_elmnts_mesh,        &
             local_dim_ibelm_xmin,      local_dim_nodes_ibelm_xmin,   &
             local_dim_ibelm_xmax,      local_dim_nodes_ibelm_xmax,   &
             local_dim_ibelm_ymin,      local_dim_nodes_ibelm_ymin,   &
             local_dim_ibelm_ymax,      local_dim_nodes_ibelm_ymax,   &
             local_dim_ibelm_bottom,    local_dim_nodes_ibelm_bottom, &
             local_dim_ibelm_top,       local_dim_nodes_ibelm_top,    &
             local_dim_neighbours_mesh, local_dim_num_elmnts_mesh,    &
             local_dim_interfaces_mesh

  integer :: num_xmin, num_xmax, num_ymin, num_ymax, num_top, num_bottom,num, &
             num_spec, num_int

  sel_num = 0

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = adjustl(LOCAL_PATH)
  database_name = database_name(1:len_trim(database_name)) // "/Database.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                               "verbose=1", ier)
  call adios_read_open_file (handle, database_name, 0, MPI_COMM_WORLD, ier)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  call adios_get_scalar(handle, "nodes_coords/local_dim", &
                        local_dim_nodes_coords, ier)
  call adios_get_scalar(handle, "matpropl/local_dim", &
                        local_dim_matpropl, ier)
  call adios_get_scalar(handle, "material_index/local_dim", &
                        local_dim_material_index, ier)
  call adios_get_scalar(handle, "elmnts_mesh/local_dim", &
                        local_dim_elmnts_mesh, ier)

  call adios_get_scalar(handle, "ibelm_xmin/local_dim", &
                        local_dim_ibelm_xmin, ier)
  call adios_get_scalar(handle, "nodes_ibelm_xmin/local_dim", &
                        local_dim_nodes_ibelm_xmin, ier)
  call adios_get_scalar(handle, "ibelm_xmax/local_dim", &
                        local_dim_ibelm_xmax, ier)
  call adios_get_scalar(handle, "nodes_ibelm_xmax/local_dim", &
                        local_dim_nodes_ibelm_xmax, ier)
  call adios_get_scalar(handle, "ibelm_ymin/local_dim", &
                        local_dim_ibelm_ymin, ier)
  call adios_get_scalar(handle, "nodes_ibelm_ymin/local_dim", &
                        local_dim_nodes_ibelm_ymin, ier)
  call adios_get_scalar(handle, "ibelm_ymax/local_dim", &
                        local_dim_ibelm_ymax, ier)
  call adios_get_scalar(handle, "nodes_ibelm_ymax/local_dim", &
                        local_dim_nodes_ibelm_ymax, ier)
  call adios_get_scalar(handle, "ibelm_bottom/local_dim", &
                        local_dim_ibelm_bottom, ier)
  call adios_get_scalar(handle, "nodes_ibelm_bottom/local_dim", &
                        local_dim_nodes_ibelm_bottom, ier)
  call adios_get_scalar(handle, "ibelm_top/local_dim", &
                        local_dim_ibelm_top, ier)
  call adios_get_scalar(handle, "nodes_ibelm_top/local_dim", &
                        local_dim_nodes_ibelm_top, ier)

  call adios_get_scalar(handle, "neighbours_mesh/local_dim", &
                        local_dim_neighbours_mesh, ier)
  call adios_get_scalar(handle, "num_elmnts_mesh/local_dim", &
                        local_dim_num_elmnts_mesh, ier)
  call adios_get_scalar(handle, "interfaces_mesh/local_dim", &
                        local_dim_interfaces_mesh, ier)

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  !------------------------------------------------------------------'
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "/nglob", 0, 1, nnodes_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "/nspec", 0, 1, nelmnts_ext_mesh, ier)
  ! read physical properties of the materials
  ! added poroelastic properties and filled with 0 the last 10 entries
  ! for elastic/acoustic
  call adios_schedule_read(handle, sel, "/nmaterials", 0, 1, nmat_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "/nundef_materials", 0, 1, &
                           nundefMat_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_xmin", 0, 1, nspec2D_xmin, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_xmax", 0, 1, nspec2D_xmax, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_ymin", 0, 1, nspec2D_ymin, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_ymax", 0, 1, nspec2D_ymax, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_bottom", 0, 1, &
                           nspec2D_bottom_ext, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_top", 0, 1, &
                           nspec2D_top_ext, ier)
  ! MPI interfaces between different partitions
  num_interfaces_ext_mesh = 0
  max_interface_size_ext_mesh = 0
  if( NPROC > 1 ) then
    ! format: #number_of_MPI_interfaces
    !         #maximum_number_of_elements_on_each_interface
    call adios_schedule_read(handle, sel, "/nb_interfaces", 0, 1, &
                             num_interfaces_ext_mesh, ier)
    call adios_schedule_read(handle, sel, "/nspec_interfaces_max", 0, 1, &
                             max_interface_size_ext_mesh, ier)
  endif
  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)

  NSPEC_AB = nelmnts_ext_mesh
  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext

  !---------------------------------------------.
  ! Allocate arrays with previously read values |
  !---------------------------------------------'
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array nodes_coords_ext_mesh'
  allocate(materials_ext_mesh(16,nmat_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array materials_ext_mesh'
  allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array undef_mat_prop'
  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array mat_ext_mesh'

  allocate(ibelm_xmin(nspec2D_xmin),&
           nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_xmin etc.'
  allocate(ibelm_xmax(nspec2D_xmax),&
           nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_xmax etc.'
  allocate(ibelm_ymin(nspec2D_ymin), &
           nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_ymin'
  allocate(ibelm_ymax(nspec2D_ymax),&
           nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_ymax etc.'
  allocate(ibelm_bottom(nspec2D_bottom_ext),&
           nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom_ext),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_bottom etc.'
  allocate(ibelm_top(nspec2D_top_ext), &
           nodes_ibelm_top(NGNOD2D,nspec2D_top_ext),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_top etc.'
  ! allocates interfaces
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_neighbours_ext_mesh'
  allocate(my_nelmnts_neighbours_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_nelmnts_neighbours_ext_mesh'
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh, &
                                  num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_interfaces_ext_mesh'
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh, &
                                    num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool_interfaces_ext_mesh'
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array nibool_interfaces_ext_mesh'

  !---------------------------.
  ! Read arrays from Database |
  !---------------------------'
  start(1) = local_dim_nodes_coords * myrank
  count(1) = 3 * nnodes_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  ! Remember: "nodes_coords" has been transposed in meshfem3D
  call adios_schedule_read(handle, sel, "nodes_coords/array", 0, 1, &
                           nodes_coords_ext_mesh, ier)

  start(1) = local_dim_matpropl* myrank
  count(1) = 16 * nmat_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  ! Remember: "matpropl" has been transposed in meshfem3D
  call adios_schedule_read(handle, sel, "matpropl/array", 0, 1, &
                           materials_ext_mesh, ier)

  start(1) = local_dim_material_index * myrank
  count(1) = 2 * nelmnts_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "material_index/array", 0, 1, &
                           mat_ext_mesh, ier)

  start(1) = local_dim_elmnts_mesh * myrank
  count(1) = NGNOD * nelmnts_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "elmnts_mesh/array", 0, 1, &
                           elmnts_ext_mesh, ier)

  start(1) = local_dim_ibelm_xmin * myrank
  count(1) = nspec2D_xmin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_xmin/array", 0, 1, &
                           ibelm_xmin, ier)

  start(1) = local_dim_nodes_ibelm_xmin * myrank
  count(1) = nspec2D_xmin * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_xmin/array", 0, 1, &
                           nodes_ibelm_xmin, ier)
  start(1) = local_dim_ibelm_xmax * myrank
  count(1) = nspec2D_xmax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_xmax/array", 0, 1, &
                           ibelm_xmax, ier)
  start(1) = local_dim_nodes_ibelm_xmax * myrank
  count(1) = nspec2D_xmax * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_xmax/array", 0, 1, &
                           nodes_ibelm_xmax, ier)
  start(1) = local_dim_ibelm_ymin * myrank
  count(1) = nspec2D_ymin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_ymin/array", 0, 1, &
                           ibelm_ymin, ier)
  start(1) = local_dim_nodes_ibelm_ymin * myrank
  count(1) = nspec2D_ymin * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_ymin/array", 0, 1, &
                           nodes_ibelm_ymin, ier)
  start(1) = local_dim_ibelm_ymax * myrank
  count(1) = nspec2D_ymax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_ymax/array", 0, 1, &
                           ibelm_ymax, ier)
  start(1) = local_dim_nodes_ibelm_ymax * myrank
  count(1) = nspec2D_ymax * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_ymax/array", 0, 1, &
                           nodes_ibelm_ymax, ier)
  start(1) = local_dim_ibelm_bottom * myrank
  count(1) = nspec2D_bottom_ext
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_bottom/array", 0, 1, &
                           ibelm_bottom, ier)
  start(1) = local_dim_nodes_ibelm_bottom * myrank
  count(1) = nspec2D_bottom_ext * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_bottom/array", 0, 1, &
                           nodes_ibelm_bottom, ier)
  start(1) = local_dim_ibelm_top * myrank
  count(1) = nspec2D_top_ext
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "ibelm_top/array", 0, 1, &
                           ibelm_top, ier)
  start(1) = local_dim_nodes_ibelm_top * myrank
  count(1) = nspec2D_top_ext * NGNOD2D
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "nodes_ibelm_top/array", 0, 1, &
                           nodes_ibelm_top, ier)

  start(1) = local_dim_neighbours_mesh * myrank
  count(1) = num_interfaces_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "neighbours_mesh/array", 0, 1, &
                           my_neighbours_ext_mesh, ier)
  start(1) = local_dim_num_elmnts_mesh * myrank
  count(1) = num_interfaces_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "num_elmnts_mesh/array", 0, 1, &
                           my_nelmnts_neighbours_ext_mesh, ier)
  start(1) = local_dim_interfaces_mesh * myrank
  count(1) = 6 * num_interfaces_ext_mesh * max_interface_size_ext_mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(handle, sel, "interfaces_mesh/array", 0, 1, &
                           my_interfaces_ext_mesh, ier)

  call adios_perform_reads(handle, ier)
  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)

  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml = 0

  ! TODO make 'undef_mat_prop' ADIOS compliant.
  ! TODO CPML in meshfem w/ and w/o ADIOS
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
  !if(myrank == 0) then
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
  !call sync_all()

end subroutine read_partition_files_adios
