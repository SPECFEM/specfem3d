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


  subroutine read_mesh_for_init_ADIOS()

  use constants, only: MAX_STRING_LEN,myrank

  use specfem_par, only: LOCAL_PATH, NSPEC_AB,NGLOB_AB,NSPEC_IRREGULAR

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  implicit none
  ! Local variables
  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle, sel
  integer :: comm
  integer :: ier

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/external_mesh.bp"

  ! gets MPI communicator
  comm = comm_adios

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)
  call check_adios_err(myrank,ier)

  call adios_read_open_file (handle, database_name, 0, comm, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0) call abort_mpi()

  !------------------------------------.
  ! Read variables from the adios file |
  !------------------------------------'
  call adios_selection_writeblock(sel, myrank)

  call adios_schedule_read(handle, sel, "nspec", 0, 1, NSPEC_AB, ier)
  call check_adios_err(myrank,ier)
  call adios_schedule_read(handle, sel, "nglob", 0, 1, NGLOB_AB, ier)
  call check_adios_err(myrank,ier)
  call adios_schedule_read(handle, sel, "nspec_irregular", 0, 1, NSPEC_IRREGULAR, ier)
  call check_adios_err(myrank,ier)

  !--------------------------------------------.
  ! Perform the reads and close the adios file |
  !--------------------------------------------'
  call adios_perform_reads(handle, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0) call abort_mpi()

  ! frees selection
  call adios_selection_delete(sel)

  call adios_read_close(handle,ier)
  call check_adios_err(myrank,ier)

  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  end subroutine read_mesh_for_init_ADIOS

!==============================================================================

  subroutine read_mesh_databases_adios()

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie

  implicit none
  integer :: ier !,inum

  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num,isel
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer(kind=8) :: local_dim_ibool, local_dim_irregular_element_number, &
    local_dim_x_global, local_dim_y_global, local_dim_z_global, &
    local_dim_xixstore, local_dim_xiystore, local_dim_xizstore, &
    local_dim_etaxstore, local_dim_etaystore, local_dim_etazstore, &
    local_dim_gammaxstore, local_dim_gammaystore, local_dim_gammazstore, &
    local_dim_jacobianstore, &
    local_dim_kappastore, local_dim_mustore, local_dim_rhostore, &
    local_dim_ispec_is_acoustic, local_dim_ispec_is_elastic, &
    local_dim_ispec_is_poroelastic, local_dim_rmass, &
    local_dim_rmass_ocean_load, local_dim_rmass_acoustic, &
    local_dim_rho_vp, local_dim_rho_vs, local_dim_abs_boundary_ispec, &
    local_dim_abs_boundary_ijk, local_dim_abs_boundary_jacobian2Dw, &
    local_dim_abs_boundary_normal, local_dim_ibelm_xmin, &
    local_dim_ibelm_ymin, local_dim_ibelm_bottom, &
    local_dim_ibelm_top, local_dim_free_surface_ispec, &
    local_dim_free_surface_ijk, local_dim_free_surface_jacobian2Dw, &
    local_dim_free_surface_normal, local_dim_coupling_ac_el_ispec, &
    local_dim_coupling_ac_el_ijk, &
    local_dim_coupling_ac_el_jacobian2Dw, &
    local_dim_coupling_ac_el_normal, local_dim_my_neighbors_ext_mesh, &
    local_dim_nibool_interfaces_ext_mesh, &
    local_dim_ibool_interfaces_ext_mesh, &
    local_dim_ispec_is_inner, local_dim_phase_ispec_inner_acoustic, &
    local_dim_phase_ispec_inner_elastic, local_dim_ibelm_xmax, &
    local_dim_ibelm_ymax, local_dim_rmass_solid_poroelastic, &
    local_dim_rmass_fluid_poroelastic, local_dim_rhoarraystore, &
    local_dim_kappaarraystore, local_dim_permstore, &
    local_dim_etastore, local_dim_tortstore, local_dim_phistore, &
    local_dim_rho_vpI, local_dim_rho_vpII, local_dim_rho_vsI, &
    local_dim_CPML_regions, local_dim_CPML_to_spec, local_dim_is_CPML, &
    local_dim_d_store_x, local_dim_d_store_y, local_dim_d_store_z, &
    local_dim_k_store_x, local_dim_k_store_y, local_dim_k_store_z, &
    local_dim_alpha_store_x, local_dim_alpha_store_y, &
    local_dim_alpha_store_z, &
    local_dim_points_interface_PML_acoustic, &
    local_dim_points_interface_PML_elastic, &
    local_dim_rmassx, local_dim_rmassy, local_dim_rmassz, &
    local_dim_rmassz_acoustic, local_dim_coupling_el_po_ispec, &
    local_dim_coupling_po_el_ispec, local_dim_coupling_el_po_ijk, &
    local_dim_coupling_po_el_ijk, &
    local_dim_coupling_el_po_jacobian2Dw, &
    local_dim_coupling_el_po_normal, local_dim_c11store, &
    local_dim_phase_ispec_inner_poroelastic, &
    local_dim_num_elem_colors_acoustic, &
    local_dim_num_elem_colors_elastic, local_dim_coupling_ac_po_ispec, &
    local_dim_coupling_ac_po_ijk, &
    local_dim_coupling_ac_po_jacobian2Dw, &
    local_dim_coupling_ac_po_normal, &
    local_dim_ispec_is_surface_external_mesh, &
    local_dim_iglob_is_surface_external_mesh

  integer :: comm
  integer :: nspec_ext,nglob_ext,nspec_irregular_ext

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "Reading mesh databases..."
    write(IMAIN,*) "  reads ADIOS mesh file: external_mesh.bp"
    write(IMAIN,*) "  from directory       : ",trim(LOCAL_PATH)
    call flush_IMAIN()
  endif

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0
  database_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/external_mesh.bp"

  ! gets MPI communicator
  comm = comm_adios

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)
  call check_adios_err(myrank,ier)

  call adios_read_open_file (handle, database_name, 0, comm, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0) call abort_mpi()

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  ! ONLY NSPEC_AB and NGLOB_AB
  !------------------------------------------------------------------'
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "nspec", 0, 1, nspec_ext, ier)
  call adios_schedule_read(handle, sel, "nglob", 0, 1, nglob_ext, ier)
  call adios_schedule_read(handle, sel, "nspec_irregular", 0, 1, nspec_irregular_ext, ier)

  call adios_perform_reads(handle, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0) call abort_mpi()

  ! checks if NSPEC,.. values match with initialization
  if (NSPEC_AB /= nspec_ext) &
    call exit_MPI(myrank,"Error: NGLOB_AB differs from initial read_mesh_for_init_ADIOS()")
  if (NGLOB_AB /= nglob_ext) &
    call exit_MPI(myrank,"Error: NGLOB_AB differs from initial read_mesh_for_init_ADIOS()")
  if (NSPEC_IRREGULAR /= nspec_irregular_ext) &
    call exit_MPI(myrank,"Error: NGLOB_AB differs from initial read_mesh_for_init_ADIOS()")

  !----------------------------------------------.
  ! Fetch values to compute the simulation type. |
  !----------------------------------------------'
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  call adios_get_scalar(handle, "ispec_is_acoustic/local_dim", local_dim_ispec_is_acoustic,ier)
  call adios_get_scalar(handle, "ispec_is_elastic/local_dim", local_dim_ispec_is_elastic,ier)
  call adios_get_scalar(handle, "ispec_is_poroelastic/local_dim", local_dim_ispec_is_poroelastic,ier)
  call check_adios_err(myrank,ier)

  start(1) = local_dim_ispec_is_acoustic * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ispec_is_acoustic/array", 0, 1, ispec_is_acoustic, ier)
  call adios_schedule_read(handle, sel, "ispec_is_elastic/array", 0, 1, ispec_is_elastic, ier)
  call adios_schedule_read(handle, sel, "ispec_is_poroelastic/array", 0, 1, ispec_is_poroelastic, ier)
  call check_adios_err(myrank,ier)

  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)
  call check_adios_err(myrank,ier)
  if (ier /= 0) call abort_mpi()

  ! all processes will have acoustic_simulation set if any flag is .true.
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  ! elastic simulation
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  ! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )

  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))

  ! checks simulation types are valid
  if ((.not. ACOUSTIC_SIMULATION) .and. &
      (.not. ELASTIC_SIMULATION) .and. &
      (.not. POROELASTIC_SIMULATION )) then
    call exit_mpi(myrank,'error no simulation type defined')
  endif

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  !------------------------------------------------------------------'
  sel_num = 1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)

  NSPEC_CPML = 0
  if (PML_CONDITIONS) then
    call adios_schedule_read(handle, sel, "nspec_cpml", 0, 1, nspec_cpml, ier)
    call adios_schedule_read(handle, sel, "CPML_width_x", 0, 1, CPML_width_x, ier)
    call adios_schedule_read(handle, sel, "CPML_width_y", 0, 1, CPML_width_y, ier)
    call adios_schedule_read(handle, sel, "CPML_width_x", 0, 1, CPML_width_z, ier)
    if (nspec_cpml > 0) then
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        call adios_schedule_read(handle, sel, "nglob_interface_PML_acoustic", 0, 1, nglob_interface_PML_acoustic, ier)
        call adios_schedule_read(handle, sel, "nglob_interface_PML_elastic", 0, 1, nglob_interface_PML_elastic, ier)
      endif
    endif
  endif

  call adios_schedule_read(handle, sel, "num_abs_boundary_faces", 0, 1, num_abs_boundary_faces, ier)

  call adios_schedule_read(handle, sel, "nspec2d_xmin", 0, 1, nspec2d_xmin, ier)
  call adios_schedule_read(handle, sel, "nspec2d_xmax", 0, 1, nspec2d_xmax, ier)
  call adios_schedule_read(handle, sel, "nspec2d_ymin", 0, 1, nspec2d_ymin, ier)
  call adios_schedule_read(handle, sel, "nspec2d_ymax", 0, 1, nspec2d_ymax, ier)
  call adios_schedule_read(handle, sel, "nspec2d_bottom", 0, 1, nspec2d_bottom, ier)
  call adios_schedule_read(handle, sel, "nspec2d_top", 0, 1, nspec2d_top, ier)

  call adios_schedule_read(handle, sel, "num_free_surface_faces", 0, 1, num_free_surface_faces, ier)
  call adios_schedule_read(handle, sel, "num_coupling_ac_el_faces", 0, 1, num_coupling_ac_el_faces, ier)
  call adios_schedule_read(handle, sel, "num_coupling_ac_po_faces", 0, 1, num_coupling_ac_po_faces, ier)
  call adios_schedule_read(handle, sel, "num_coupling_el_po_faces", 0, 1, num_coupling_el_po_faces, ier)

  call adios_schedule_read(handle, sel, "num_interfaces_ext_mesh", 0, 1, num_interfaces_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "max_nibool_interfaces_ext_mesh", 0, 1, max_nibool_interfaces_ext_mesh, ier)

  if (ACOUSTIC_SIMULATION) then
    call adios_schedule_read(handle, sel, "nspec_inner_acoustic", 0, 1, nspec_inner_acoustic, ier)
    call adios_schedule_read(handle, sel, "nspec_outer_acoustic", 0, 1, nspec_outer_acoustic, ier)
    call adios_schedule_read(handle, sel, "num_phase_ispec_acoustic", 0, 1, num_phase_ispec_acoustic, ier)
  endif

  if (ELASTIC_SIMULATION) then
    call adios_schedule_read(handle, sel, "nspec_inner_elastic", 0, 1, nspec_inner_elastic, ier)
    call adios_schedule_read(handle, sel, "nspec_outer_elastic", 0, 1, nspec_outer_elastic, ier)
    call adios_schedule_read(handle, sel, "num_phase_ispec_elastic", 0, 1, num_phase_ispec_elastic, ier)
  endif

  if (POROELASTIC_SIMULATION) then
    call adios_schedule_read(handle, sel, "nspec_inner_poroelastic", 0, 1, nspec_inner_poroelastic, ier)
    call adios_schedule_read(handle, sel, "nspec_outer_poroelastic", 0, 1, nspec_outer_poroelastic, ier)
    call adios_schedule_read(handle, sel, "num_phase_ispec_poroelastic", 0, 1, num_phase_ispec_poroelastic, ier)
  endif

  num_colors_outer_acoustic = 0
  num_colors_inner_acoustic = 0
  num_colors_outer_elastic = 0
  num_colors_inner_elastic = 0
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      call adios_schedule_read(handle, sel, "num_colors_outer_acoustic", 0, 1, num_colors_outer_acoustic, ier)
      call adios_schedule_read(handle, sel, "num_colors_outer_acoustic", 0, 1, num_colors_inner_acoustic, ier)
    endif
    if (ELASTIC_SIMULATION) then
      call adios_schedule_read(handle, sel, "num_colors_outer_elastic", 0, 1, num_colors_outer_elastic, ier)
      call adios_schedule_read(handle, sel, "num_colors_outer_elastic", 0, 1, num_colors_inner_elastic, ier)
    endif
  endif

  call adios_schedule_read(handle, sel, "nfaces_surface", 0, 1, nfaces_surface, ier)

  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  call adios_get_scalar(handle, "ibool/local_dim", local_dim_ibool,ier)

  call adios_get_scalar(handle, "x_global/local_dim", local_dim_x_global,ier)
  call adios_get_scalar(handle, "y_global/local_dim", local_dim_y_global,ier)
  call adios_get_scalar(handle, "z_global/local_dim", local_dim_z_global,ier)

  call adios_get_scalar(handle, "xixstore/local_dim", local_dim_xixstore,ier)
  call adios_get_scalar(handle, "xiystore/local_dim", local_dim_xiystore,ier)
  call adios_get_scalar(handle, "xizstore/local_dim", local_dim_xizstore,ier)
  call adios_get_scalar(handle, "etaxstore/local_dim", local_dim_etaxstore,ier)
  call adios_get_scalar(handle, "etaystore/local_dim", local_dim_etaystore,ier)
  call adios_get_scalar(handle, "etazstore/local_dim", local_dim_etazstore,ier)
  call adios_get_scalar(handle, "gammaxstore/local_dim", local_dim_gammaxstore,ier)
  call adios_get_scalar(handle, "gammaystore/local_dim", local_dim_gammaystore,ier)
  call adios_get_scalar(handle, "gammazstore/local_dim", local_dim_gammazstore,ier)
  call adios_get_scalar(handle, "jacobianstore/local_dim", local_dim_jacobianstore,ier)

  call adios_get_scalar(handle, "kappastore/local_dim", local_dim_kappastore,ier)
  call adios_get_scalar(handle, "mustore/local_dim", local_dim_mustore,ier)
  call adios_get_scalar(handle, "rhostore/local_dim", local_dim_rhostore,ier)

  call adios_get_scalar(handle, "irregular_element_number/local_dim", local_dim_irregular_element_number,ier)

  ! single values
  call adios_schedule_read(handle, sel, "jacobian_regular", 0, 1, jacobian_regular, ier)
  call adios_schedule_read(handle, sel, "xix_regular", 0, 1, xix_regular, ier)

  if (ACOUSTIC_SIMULATION) then
    call adios_get_scalar(handle, "rmass_acoustic/local_dim", local_dim_rmass_acoustic,ier)
  endif
  if (ELASTIC_SIMULATION) then
    call adios_get_scalar(handle, "rmass/local_dim", local_dim_rmass,ier)

    if (APPROXIMATE_OCEAN_LOAD) then
      call adios_get_scalar(handle, "rmass_ocean_load/local_dim", local_dim_rmass_ocean_load,ier)
    endif
    call adios_get_scalar(handle, "rho_vp/local_dim", local_dim_rho_vp,ier)
    call adios_get_scalar(handle, "rho_vs/local_dim", local_dim_rho_vs,ier)
  endif
  if (POROELASTIC_SIMULATION) then
    call adios_get_scalar(handle, "rmass_solid_poroelastic/local_dim", local_dim_rmass_solid_poroelastic,ier)
    call adios_get_scalar(handle, "rmass_fluid_poroelastic/local_dim", local_dim_rmass_fluid_poroelastic,ier)
    call adios_get_scalar(handle, "rhoarraystore/local_dim", local_dim_rhoarraystore,ier)
    call adios_get_scalar(handle, "kappaarraystore/local_dim", local_dim_kappaarraystore,ier)
    call adios_get_scalar(handle, "permstore/local_dim", local_dim_permstore,ier)
    call adios_get_scalar(handle, "etastore/local_dim", local_dim_etastore,ier)
    call adios_get_scalar(handle, "tortstore/local_dim", local_dim_tortstore,ier)
    call adios_get_scalar(handle, "phistore/local_dim", local_dim_phistore,ier)
    call adios_get_scalar(handle, "rho_vpI/local_dim", local_dim_rho_vpI,ier)
    call adios_get_scalar(handle, "rho_vpII/local_dim", local_dim_rho_vpII,ier)
    call adios_get_scalar(handle, "rho_vsI/local_dim", local_dim_rho_vsI,ier)
  endif
  if (PML_CONDITIONS) then
    if (nspec_cpml > 0) then
      call adios_get_scalar(handle, "CPML_regions/local_dim", local_dim_CPML_regions, ier)
      call adios_get_scalar(handle, "CPML_to_spec/local_dim", local_dim_CPML_to_spec, ier)
      call adios_get_scalar(handle, "is_CPML/local_dim", local_dim_is_CPML, ier)
      call adios_get_scalar(handle, "d_store_x/local_dim", local_dim_d_store_x, ier)
      call adios_get_scalar(handle, "d_store_y/local_dim", local_dim_d_store_y, ier)
      call adios_get_scalar(handle, "d_store_z/local_dim", local_dim_d_store_z, ier)
      call adios_get_scalar(handle, "k_store_x/local_dim", local_dim_k_store_x, ier)
      call adios_get_scalar(handle, "k_store_y/local_dim", local_dim_k_store_y, ier)
      call adios_get_scalar(handle, "k_store_z/local_dim", local_dim_k_store_z, ier)
      call adios_get_scalar(handle, "alpha_store_x/local_dim", local_dim_alpha_store_x, ier)
      call adios_get_scalar(handle, "alpha_store_y/local_dim", local_dim_alpha_store_y, ier)
      call adios_get_scalar(handle, "alpha_store_z/local_dim", local_dim_alpha_store_z, ier)
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          call adios_get_scalar(handle, "points_interface_PML_acoustic/local_dim", &
                                local_dim_points_interface_PML_acoustic, ier)
        endif
        if (nglob_interface_PML_elastic > 0) then
          call adios_get_scalar(handle, "points_interface_PML_elastic/local_dim", &
                                local_dim_points_interface_PML_elastic, ier)
        endif
      endif
    endif
  endif

  ! absorbing boundaries
  if (num_abs_boundary_faces > 0) then
    call adios_get_scalar(handle, "abs_boundary_ispec/local_dim", local_dim_abs_boundary_ispec,ier)
    call adios_get_scalar(handle, "abs_boundary_ijk/local_dim", local_dim_abs_boundary_ijk,ier)
    call adios_get_scalar(handle, "abs_boundary_jacobian2Dw/local_dim", local_dim_abs_boundary_jacobian2Dw,ier)
    call adios_get_scalar(handle, "abs_boundary_normal/local_dim", local_dim_abs_boundary_normal,ier)
    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! mass matrix contributions
      if (ELASTIC_SIMULATION) then
        call adios_get_scalar(handle, "rmassx/local_dim", local_dim_rmassx, ier)
        call adios_get_scalar(handle, "rmassy/local_dim", local_dim_rmassy, ier)
        call adios_get_scalar(handle, "rmassz/local_dim", local_dim_rmassz, ier)
      endif
      if (ACOUSTIC_SIMULATION) then
        call adios_get_scalar(handle, "rmassz_acoustic/local_dim", local_dim_rmassz_acoustic, ier)
      endif
    endif
  endif

  if (nspec2d_xmin > 0) then
    call adios_get_scalar(handle, "ibelm_xmin/local_dim", local_dim_ibelm_xmin,ier)
  endif
  if (nspec2d_xmax > 0) then
    call adios_get_scalar(handle, "ibelm_xmax/local_dim", local_dim_ibelm_xmax,ier)
  endif
  if (nspec2d_ymin > 0) then
    call adios_get_scalar(handle, "ibelm_ymin/local_dim", local_dim_ibelm_ymin,ier)
  endif
  if (nspec2d_ymax > 0) then
    call adios_get_scalar(handle, "ibelm_ymax/local_dim", local_dim_ibelm_ymax,ier)
  endif
  if (nspec2d_bottom > 0) then
    call adios_get_scalar(handle, "ibelm_bottom/local_dim", local_dim_ibelm_bottom,ier)
  endif
  if (nspec2d_top > 0) then
    call adios_get_scalar(handle, "ibelm_top/local_dim", local_dim_ibelm_top,ier)
  endif

  if (num_free_surface_faces > 0) then
    call adios_get_scalar(handle, "free_surface_ispec/local_dim", local_dim_free_surface_ispec,ier)
    call adios_get_scalar(handle, "free_surface_ijk/local_dim", local_dim_free_surface_ijk,ier)
    call adios_get_scalar(handle, "free_surface_jacobian2Dw/local_dim", local_dim_free_surface_jacobian2Dw,ier)
    call adios_get_scalar(handle, "free_surface_normal/local_dim", local_dim_free_surface_normal,ier)
  endif
  if (num_coupling_ac_el_faces > 0) then
    call adios_get_scalar(handle, "coupling_ac_el_ispec/local_dim", local_dim_coupling_ac_el_ispec,ier)
    call adios_get_scalar(handle, "coupling_ac_el_ijk/local_dim", local_dim_coupling_ac_el_ijk,ier)
    call adios_get_scalar(handle, "coupling_ac_el_jacobian2Dw/local_dim", local_dim_coupling_ac_el_jacobian2Dw,ier)
    call adios_get_scalar(handle, "coupling_ac_el_normal/local_dim", local_dim_coupling_ac_el_normal,ier)
  endif
  if (num_coupling_ac_po_faces > 0) then
    call adios_get_scalar(handle, "coupling_ac_po_ispec/local_dim", local_dim_coupling_ac_po_ispec, ier)
    call adios_get_scalar(handle, "coupling_ac_po_ijk/local_dim", local_dim_coupling_ac_po_ijk, ier)
    call adios_get_scalar(handle, "coupling_ac_po_jacobian2Dw/local_dim", local_dim_coupling_ac_po_jacobian2Dw, ier)
    call adios_get_scalar(handle, "coupling_ac_po_normal/local_dim", local_dim_coupling_ac_po_normal, ier)
  endif
  if (num_coupling_el_po_faces > 0) then
    call adios_get_scalar(handle, "coupling_el_po_ispec/local_dim", local_dim_coupling_el_po_ispec, ier)
    call adios_get_scalar(handle, "coupling_po_el_ispec/local_dim", local_dim_coupling_po_el_ispec, ier)
    call adios_get_scalar(handle, "coupling_el_po_ijk/local_dim", local_dim_coupling_el_po_ijk, ier)
    call adios_get_scalar(handle, "coupling_po_el_ijk/local_dim", local_dim_coupling_po_el_ijk, ier)
    call adios_get_scalar(handle, "coupling_el_po_jacobian2Dw/local_dim", local_dim_coupling_el_po_jacobian2Dw, ier)
    call adios_get_scalar(handle, "coupling_el_po_normal/local_dim", local_dim_coupling_el_po_normal, ier)
  endif
  if (num_interfaces_ext_mesh > 0) then
    call adios_get_scalar(handle, "my_neighbors_ext_mesh/local_dim", local_dim_my_neighbors_ext_mesh,ier)
    call adios_get_scalar(handle, "nibool_interfaces_ext_mesh/local_dim", local_dim_nibool_interfaces_ext_mesh,ier)
    call adios_get_scalar(handle, "ibool_interfaces_ext_mesh_dummy/local_dim", local_dim_ibool_interfaces_ext_mesh, ier)
  endif
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    call adios_get_scalar(handle, "c11store/local_dim", local_dim_c11store, ier)
  endif
  call adios_get_scalar(handle, "ispec_is_inner/local_dim", local_dim_ispec_is_inner,ier)
  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic > 0) then
      call adios_get_scalar(handle, "phase_ispec_inner_acoustic/local_dim", local_dim_phase_ispec_inner_acoustic,ier)
    endif
  endif
  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic > 0) then
      call adios_get_scalar(handle, "phase_ispec_inner_elastic/local_dim", local_dim_phase_ispec_inner_elastic,ier)
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic > 0) then
      call adios_get_scalar(handle, "phase_ispec_inner_poroelastic/local_dim", local_dim_phase_ispec_inner_poroelastic,ier)
    endif
  endif
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      call adios_get_scalar(handle, "num_elem_colors_acoustic/local_dim", local_dim_num_elem_colors_acoustic, ier)
    endif
    if (ELASTIC_SIMULATION) then
      call adios_get_scalar(handle, "num_elem_colors_elastic/local_dim", local_dim_num_elem_colors_elastic, ier)
    endif
  endif
  ! for mesh surface
  call adios_get_scalar(handle, "ispec_is_surface_external_mesh/local_dim", local_dim_ispec_is_surface_external_mesh, ier)
  call adios_get_scalar(handle, "iglob_is_surface_external_mesh/local_dim", local_dim_iglob_is_surface_external_mesh, ier)

!TODO
#if 1
  !---------------------------------------------.
  ! Allocate arrays with previously read values |
  !---------------------------------------------'
  if (ACOUSTIC_SIMULATION) then
    ! potentials
    allocate(potential_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1791')
    if (ier /= 0) stop 'error allocating array potential_acoustic'
    allocate(potential_dot_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1792')
    if (ier /= 0) stop 'error allocating array potential_dot_acoustic'
    allocate(potential_dot_dot_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1793')
    if (ier /= 0) stop 'error allocating array potential_dot_dot_acoustic'
    !if (SIMULATION_TYPE /= 1) then
    !  allocate(potential_acoustic_adj_coupling(NGLOB_AB),stat=ier) ! not used yet
    !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1794')
    !  if (ier /= 0) stop 'error allocating array potential_acoustic_adj_coupling'
    !endif
    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1795')
    if (ier /= 0) stop 'error allocating array rmass_acoustic'

    ! initializes mass matrix contribution
    allocate(rmassz_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1796')
    if (ier /= 0) stop 'error allocating array rmassz_acoustic'
    rmassz_acoustic(:) = 0._CUSTOM_REAL
  endif

!TODO
#endif
  ! elastic simulation
  if (ELASTIC_SIMULATION) then
!TODO
#if 1
    ! displacement,velocity,acceleration
    allocate(displ(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1798')
    if (ier /= 0) stop 'error allocating array displ'
    allocate(veloc(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1799')
    if (ier /= 0) stop 'error allocating array veloc'
    allocate(accel(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1800')
    if (ier /= 0) stop 'error allocating array accel'
    if (SIMULATION_TYPE /= 1) then
      allocate(accel_adj_coupling(NDIM,NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1801')
      if (ier /= 0) stop 'error allocating array accel_adj_coupling'
    endif

    ! allocates mass matrix
    allocate(rmass(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1802')

    if (ier /= 0) stop 'error allocating array rmass'
    ! initializes mass matrix contributions
    allocate(rmassx(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1803')
    allocate(rmassy(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1804')
    allocate(rmassz(NGLOB_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1805')
    if (ier /= 0) stop 'error allocating array rmassx,rmassy,rmassz'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL

    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1806')
    if (ier /= 0) stop 'error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1807')
    if (ier /= 0) stop 'error allocating array rho_vs'
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1808')
    allocate(c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1809')
    allocate(c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1810')
    allocate(c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1811')
    allocate(c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1812')
    allocate(c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1813')
    allocate(c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1814')
    allocate(c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1815')
    allocate(c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1816')
    allocate(c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1817')
    allocate(c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1818')
    allocate(c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1819')
    allocate(c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1820')
    allocate(c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1821')
    allocate(c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1822')
    allocate(c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1823')
    allocate(c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1824')
    allocate(c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1825')
    allocate(c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1826')
    allocate(c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1827')
    allocate(c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1828')
    if (ier /= 0) stop 'error allocating array c11store etc.'

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    allocate(R_xx(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1829')
    allocate(R_yy(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1830')
    allocate(R_xy(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1831')
    allocate(R_xz(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1832')
    allocate(R_yz(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1833')
    if (ier /= 0) stop 'error allocating array R_xx etc.'

    ! needed for attenuation and/or kernel computations
    allocate(epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1834')
    allocate(epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1835')
    allocate(epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1836')
    allocate(epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1837')
    allocate(epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1838')
    allocate(epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1839')
    if (ier /= 0) stop 'error allocating array epsilondev_xx etc.'

    allocate(R_trace(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1840')
    if (ier /= 0) stop 'error allocating array R_trace etc.'

    ! note: needed for some subroutine arguments
    allocate(epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1841')
    if (ier /= 0) stop 'error allocating array epsilon_trace_over_3'

    ! needed for attenuation
    allocate(factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1842')
    if (ier /= 0) stop 'error allocating array factor_common'

    allocate(factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1843')
    if (ier /= 0) stop 'error allocating array factor_common_kappa'

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1844')
      if (ier /= 0) stop 'error allocating array rmass_ocean_load'
    else
      ! dummy allocation
      allocate(rmass_ocean_load(1),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1845')
      if (ier /= 0) stop 'error allocating dummy array rmass_ocean_load'
    endif
! TODO
#endif
  else
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif

  if (POROELASTIC_SIMULATION) then

    if (GPU_MODE) &
      call exit_mpi(myrank,'POROELASTICITY not supported by GPU mode yet...')

    ! displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(displs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1846')
    if (ier /= 0) stop 'error allocating array displs_poroelastic'
    allocate(velocs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1847')
    if (ier /= 0) stop 'error allocating array velocs_poroelastic'
    allocate(accels_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1848')
    if (ier /= 0) stop 'error allocating array accels_poroelastic'
    allocate(displw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1849')
    if (ier /= 0) stop 'error allocating array displw_poroelastic'
    allocate(velocw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1850')
    if (ier /= 0) stop 'error allocating array velocw_poroelastic'
    allocate(accelw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1851')
    if (ier /= 0) stop 'error allocating array accelw_poroelastic'

    allocate(rmass_solid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1852')
    if (ier /= 0) stop 'error allocating array rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1853')
    if (ier /= 0) stop 'error allocating array rmass_fluid_poroelastic'

    allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1854')
    allocate(kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1855')
    allocate(etastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1856')
    allocate(tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1857')
    allocate(phistore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1858')
    allocate(permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1859')
    allocate(rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1860')
    allocate(rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1861')
    allocate(rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1862')
    if (ier /= 0) stop 'error allocating array poroelastic properties'

    ! needed for kernel computations
    allocate(epsilonsdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1863')
    allocate(epsilonsdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1864')
    allocate(epsilonsdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1865')
    allocate(epsilonsdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1866')
    allocate(epsilonsdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1867')
    allocate(epsilonwdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1868')
    allocate(epsilonwdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1869')
    allocate(epsilonwdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1870')
    allocate(epsilonwdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1871')
    allocate(epsilonwdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1872')
    if (ier /= 0) stop 'error allocating array epsilonsdev_xx etc.'

    allocate(epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1873')
    allocate(epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1874')
    if (ier /= 0) stop 'error allocating array epsilons_trace_over_3 etc.'
  else
    ! dummy allocations (needed for subroutine arguments)
    allocate(rhoarraystore(2,1,1,1,1), &
             kappaarraystore(3,1,1,1,1), &
             etastore(1,1,1,1), &
             tortstore(1,1,1,1), &
             phistore(1,1,1,1), &
             permstore(6,1,1,1,1), &
             rho_vpI(1,1,1,1), &
             rho_vpII(1,1,1,1), &
             rho_vsI(1,1,1,1))
  endif

  ! C-PML absorbing boundary conditions
  ! we allocate this array even when PMLs are absent because we need it in logical tests in "if" statements
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1875')
  if (ier /= 0) stop 'error allocating array is_CPML'
  ! make sure there are no PMLs by default,
  ! and then below if NSPEC_CPML > 0 we will read the real flags for this mesh from the disk
  is_CPML(:) = .false.

  if (PML_CONDITIONS) then
    if (NSPEC_CPML > 0) then
      allocate(CPML_regions(NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1876')
      if (ier /= 0) stop 'error allocating array CPML_regions'
      allocate(CPML_to_spec(NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1877')
      if (ier /= 0) stop 'error allocating array CPML_to_spec'
      allocate(d_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1878')
      if (ier /= 0) stop 'error allocating array d_store_x'
      allocate(d_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1879')
      if (ier /= 0) stop 'error allocating array d_store_y'
      allocate(d_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1880')
      if (ier /= 0) stop 'error allocating array d_store_z'
      allocate(K_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1881')
      if (ier /= 0) stop 'error allocating array K_store_x'
      allocate(K_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1882')
      if (ier /= 0) stop 'error allocating array K_store_y'
      allocate(K_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1883')
      if (ier /= 0) stop 'error allocating array K_store_z'
      allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1884')
      if (ier /= 0) stop 'error allocating array alpha_store_x'
      allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1885')
      if (ier /= 0) stop 'error allocating array alpha_store_y'
      allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1886')
      if (ier /= 0) stop 'error allocating array alpha_store_z'

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1887')
          if (ier /= 0) stop 'error allocating array points_interface_PML_acoustic'
        endif
        if (nglob_interface_PML_elastic > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1888')
          if (ier /= 0) stop 'error allocating array points_interface_PML_elastic'
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  allocate(abs_boundary_ispec(num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1889')
  allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1890')
  allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1891')
  allocate(abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1892')
  if (ier /= 0) stop 'error allocating array abs_boundary_ispec etc.'

  allocate(ibelm_xmin(nspec2D_xmin),ibelm_xmax(nspec2D_xmax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1893')
  allocate(ibelm_ymin(nspec2D_ymin),ibelm_ymax(nspec2D_ymax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1894')
  allocate(ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1895')
  if (ier /= 0) stop 'error allocating arrays ibelm_xmin,ibelm_xmax etc.'

  ! free surface
  allocate(free_surface_ispec(num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1896')
  allocate(free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1897')
  allocate(free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1898')
  allocate(free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1899')
  if (ier /= 0) stop 'error allocating arrays free_surface_ispec etc.'

  ! acoustic-elastic coupling surface
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1900')
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1901')
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1902')
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1903')
  if (ier /= 0) stop 'error allocating array coupling_ac_el_normal etc.'

  ! acoustic-poroelastic coupling surface
  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1904')
  allocate(coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1905')
  allocate(coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1906')
  allocate(coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1907')
  if (ier /= 0) stop 'error allocating array coupling_ac_po_normal etc.'

  ! elastic-poroelastic coupling surface
  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1908')
  allocate(coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1909')
  allocate(coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1910')
  allocate(coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1911')
  allocate(coupling_el_po_ispec(num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1912')
  allocate(coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1913')
  if (ier /= 0) stop 'error allocating array coupling_el_po_normal etc.'

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1914')
    allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1915')
    if (ier /= 0) stop 'error allocating array my_neighbors_ext_mesh etc.'
    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1916')
    if (ier /= 0) stop 'error allocating array ibool_interfaces_ext_mesh'
  else
    ! no interfaces
    max_nibool_interfaces_ext_mesh = 0
    allocate(my_neighbors_ext_mesh(1),nibool_interfaces_ext_mesh(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1915')
    allocate(ibool_interfaces_ext_mesh(1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1917')
  endif

  ! inner / outer elements
  allocate(ispec_is_inner(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1918')
  if (ier /= 0) stop 'error allocating array ispec_is_inner'

  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic < 0) stop 'error acoustic simulation:' // &
                                    'num_phase_ispec_acoustic is < zero'
    allocate(phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1919')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_acoustic'
  endif

  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic < 0) stop 'error elastic simulation:' // &
                                   'num_phase_ispec_elastic is < zero'
    allocate(phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1920')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_elastic'
  endif

  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic < 0) &
      stop 'error poroelastic simulation:num_phase_ispec_poroelastic is < zero'
    allocate(phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1921')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_poroelastic'
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! acoustic domain colors
    if (ACOUSTIC_SIMULATION) then
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1922')
      if (ier /= 0) stop 'error allocating num_elem_colors_acoustic array'
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1923')
      if (ier /= 0) stop 'error allocating num_elem_colors_elastic array'
    endif
  else
    ! allocates dummy arrays
    if (ACOUSTIC_SIMULATION) then
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1924')
      if (ier /= 0) stop 'error allocating num_elem_colors_acoustic array'
    endif
    if (ELASTIC_SIMULATION) then
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1925')
      if (ier /= 0) stop 'error allocating num_elem_colors_elastic array'
    endif
  endif

  ! for mesh surface
  allocate(ispec_is_surface_external_mesh(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1926')
  allocate(iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1927')
  if (ier /= 0) stop 'error allocating array for mesh surface'


  !-----------------------------------.
  ! Read arrays from external_mesh.bp |
  !-----------------------------------'
  start(1) = local_dim_ibool * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibool/array", 0, 1, ibool, ier)

  start(1) = local_dim_x_global * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "x_global/array", 0, 1, xstore, ier)
  call adios_schedule_read(handle, sel, "y_global/array", 0, 1, ystore, ier)
  call adios_schedule_read(handle, sel, "z_global/array", 0, 1, zstore, ier)

  start(1) = local_dim_irregular_element_number * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "irregular_element_number/array", 0, 1, irregular_element_number, ier)
  call check_adios_err(myrank,ier)

  start(1) = local_dim_xixstore * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_IRREGULAR
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "xixstore/array", 0, 1, xixstore, ier)
  call adios_schedule_read(handle, sel, "xiystore/array", 0, 1, xiystore, ier)
  call adios_schedule_read(handle, sel, "xizstore/array", 0, 1, xizstore, ier)
  call adios_schedule_read(handle, sel, "etaxstore/array", 0, 1, etaxstore, ier)
  call adios_schedule_read(handle, sel, "etaystore/array", 0, 1, etaystore, ier)
  call adios_schedule_read(handle, sel, "etazstore/array", 0, 1, etazstore, ier)
  call adios_schedule_read(handle, sel, "gammaxstore/array", 0, 1, gammaxstore, ier)
  call adios_schedule_read(handle, sel, "gammaystore/array", 0, 1, gammaystore, ier)
  call adios_schedule_read(handle, sel, "gammazstore/array", 0, 1, gammazstore, ier)
  call adios_schedule_read(handle, sel, "jacobianstore/array", 0, 1, jacobianstore, ier)

  start(1) = local_dim_kappastore * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "kappastore/array", 0, 1, kappastore, ier)
  call adios_schedule_read(handle, sel, "mustore/array", 0, 1, mustore, ier)
  call adios_schedule_read(handle, sel, "rhostore/array", 0, 1, rhostore, ier)

  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  if (ACOUSTIC_SIMULATION) then
    start(1) = local_dim_rmass_acoustic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass_acoustic/array", 0, 1, rmass_acoustic, ier)
  endif

  if (ELASTIC_SIMULATION) then
    start(1) = local_dim_rmass * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass/array", 0, 1, rmass, ier)

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      start(1) = local_dim_rmass_ocean_load * myrank
      count_ad(1) = NGLOB_AB !nglob_ocean
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "rmass_ocean_load/array", 0, 1, rmass_ocean_load, ier)
    endif

    !pll material parameters for stacey conditions
    start(1) = local_dim_rho_vp * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rho_vp/array", 0, 1, rho_vp, ier)
    call adios_schedule_read(handle, sel, "rho_vs/array", 0, 1, rho_vs, ier)
  endif

  if (POROELASTIC_SIMULATION) then
    start(1) = local_dim_rmass_solid_poroelastic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass_solid_poroelastic/array", 0, 1, rmass_solid_poroelastic, ier)
    call adios_schedule_read(handle, sel, "rmass_fluid_poroelastic/array", 0, 1, rmass_fluid_poroelastic, ier)

    start(1) = local_dim_rhoarraystore * myrank
    count_ad(1) = 2 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rhoarraystore/array", 0, 1, rhoarraystore, ier)

    start(1) = local_dim_kappaarraystore* myrank
    count_ad(1) = 3 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "kappaarraystore/array", 0, 1, kappaarraystore, ier)

    start(1) = local_dim_permstore * myrank
    count_ad(1) =  6 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "permstore/array", 0, 1, permstore, ier)

    start(1) = local_dim_etastore * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "etastore/array", 0, 1, etastore, ier)
    call adios_schedule_read(handle, sel, "tortstore/array", 0, 1, tortstore, ier)
    call adios_schedule_read(handle, sel, "phistore/array", 0, 1, phistore, ier)
    call adios_schedule_read(handle, sel, "rho_vpI/array", 0, 1, rho_vpI, ier)
    call adios_schedule_read(handle, sel, "rho_vpII/array", 0, 1, rho_vpII, ier)
    call adios_schedule_read(handle, sel, "rho_vsI/array", 0, 1, rho_vsI, ier)
  endif

  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    if (NSPEC_CPML > 0) then
      start(1) = local_dim_CPML_regions * myrank
      count_ad(1) = nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "CPML_regions/array", 0, 1, CPML_regions, ier)
      call adios_schedule_read(handle, sel, "CPML_to_spec/array", 0, 1, CPML_to_spec, ier)

      start(1) = local_dim_is_cpml * myrank
      count_ad(1) = NSPEC_AB
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "is_CPML/array", 0, 1, is_CPML, ier)

      start(1) = local_dim_d_store_x * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "d_store_x/array", 0, 1, d_store_x, ier)
      call adios_schedule_read(handle, sel, "d_store_y/array", 0, 1, d_store_y, ier)
      call adios_schedule_read(handle, sel, "d_store_z/array", 0, 1, d_store_z, ier)
      call adios_schedule_read(handle, sel, "k_store_x/array", 0, 1, k_store_x, ier)
      call adios_schedule_read(handle, sel, "k_store_y/array", 0, 1, k_store_y, ier)
      call adios_schedule_read(handle, sel, "k_store_z/array", 0, 1, k_store_z, ier)
      call adios_schedule_read(handle, sel, "alpha_store_x/array", 0, 1, alpha_store_x, ier)
      call adios_schedule_read(handle, sel, "alpha_store_y/array", 0, 1, alpha_store_y, ier)
      call adios_schedule_read(handle, sel, "alpha_store_z/array", 0, 1, alpha_store_z, ier)

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          start(1) = local_dim_points_interface_PML_acoustic* myrank
          count_ad(1) = nglob_interface_PML_acoustic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call adios_selection_boundingbox (sel , 1, start, count_ad)
          call adios_schedule_read(handle, sel, "points_interface_PML_acoustic/array", &
                                   0, 1, points_interface_PML_acoustic , ier)
        endif
        if (nglob_interface_PML_elastic > 0) then
          start(1) = local_dim_points_interface_PML_elastic* myrank
          count_ad(1) = nglob_interface_PML_elastic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call adios_selection_boundingbox (sel , 1, start, count_ad)
          call adios_schedule_read(handle, sel, "points_interface_PML_elastic/array", &
                                   0, 1, points_interface_PML_elastic , ier)
        endif
      endif
    endif
  endif

  if (num_abs_boundary_faces > 0) then
    start(1) = local_dim_abs_boundary_ispec * myrank
    count_ad(1) = num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "abs_boundary_ispec/array", 0, 1, abs_boundary_ispec, ier)

    start(1) = local_dim_abs_boundary_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "abs_boundary_ijk/array", 0, 1, abs_boundary_ijk, ier)

    start(1) = local_dim_abs_boundary_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "abs_boundary_jacobian2Dw/array", 0, 1, abs_boundary_jacobian2Dw, ier)

    start(1) = local_dim_abs_boundary_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "abs_boundary_normal/array", 0, 1, abs_boundary_normal, ier)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        start(1) = local_dim_rmassx * myrank
        count_ad(1) = NGLOB_AB  ! == nglob_xy in generate_databse
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "rmassx/array", 0, 1, rmassx, ier)
        call adios_schedule_read(handle, sel, "rmassy/array", 0, 1, rmassy, ier)
        call adios_schedule_read(handle, sel, "rmassz/array", 0, 1, rmassz, ier)
      endif
      if (ACOUSTIC_SIMULATION) then
        start(1) = local_dim_rmassz_acoustic * myrank
        count_ad(1) = NGLOB_AB ! == nglob_xy in generate_databse
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "rmassz_acoustic/array", 0, 1, rmassz_acoustic, ier)
      endif
    endif
  endif

  if (nspec2d_xmin > 0) then
    start(1) = local_dim_ibelm_xmin * myrank
    count_ad(1) = nspec2D_xmin
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_xmin/array", 0, 1, ibelm_xmin, ier)
  endif
  if (nspec2d_xmax > 0) then
    start(1) = local_dim_ibelm_xmax * myrank
    count_ad(1) = nspec2D_xmax
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_xmax/array", 0, 1, ibelm_xmax, ier)
  endif
  if (nspec2d_ymin > 0) then
    start(1) = local_dim_ibelm_ymin * myrank
    count_ad(1) = nspec2D_ymin
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_ymin/array", 0, 1, ibelm_ymin, ier)
  endif
  if (nspec2d_ymax > 0) then
    start(1) = local_dim_ibelm_ymax * myrank
    count_ad(1) = nspec2D_ymax
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_ymax/array", 0, 1, ibelm_ymax, ier)
  endif
  if (nspec2d_bottom > 0) then
    start(1) = local_dim_ibelm_bottom * myrank
    count_ad(1) = nspec2D_bottom
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_bottom/array", 0, 1, ibelm_bottom, ier)
  endif
  if (nspec2d_top > 0) then
    start(1) = local_dim_ibelm_top * myrank
    count_ad(1) = nspec2D_top
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_top/array", 0, 1, ibelm_top, ier)
  endif

  ! free surface
  if (num_free_surface_faces > 0) then
    start(1) = local_dim_free_surface_ispec * myrank
    count_ad(1) = num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ispec/array", 0, 1, free_surface_ispec, ier)

    start(1) = local_dim_free_surface_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ijk/array", 0, 1, free_surface_ijk, ier)

    start(1) = local_dim_free_surface_ijk* myrank
    count_ad(1) = NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ijk/array", 0, 1, free_surface_ijk, ier)

    start(1) = local_dim_free_surface_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_normal/array", 0, 1, free_surface_normal, ier)
  endif

  ! acoustic-elastic coupling surface
  if (num_coupling_ac_el_faces > 0) then
    start(1) = local_dim_coupling_ac_el_ispec * myrank
    count_ad(1) = num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_ispec/array", 0, 1, coupling_ac_el_ispec, ier)

    start(1) = local_dim_coupling_ac_el_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_ijk/array", 0, 1, coupling_ac_el_ijk, ier)

    start(1) = local_dim_coupling_ac_el_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_jacobian2Dw/array", 0, 1, coupling_ac_el_jacobian2Dw, ier)
    start(1) = local_dim_coupling_ac_el_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_normal/array", 0, 1, coupling_ac_el_normal, ier)
  endif

  ! acoustic-poroelastic coupling surface
  if (num_coupling_ac_po_faces > 0) then
    start(1) = local_dim_coupling_ac_po_ispec * myrank
    count_ad(1) = num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_ispec/array", 0, 1, coupling_ac_po_ispec, ier)

    start(1) = local_dim_coupling_ac_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_ijk/array", 0, 1, coupling_ac_po_ijk, ier)

    start(1) = local_dim_coupling_ac_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_jacobian2Dw/array", 0, 1, coupling_ac_po_jacobian2Dw, ier)

    start(1) = local_dim_coupling_ac_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_normal/array", 0, 1, coupling_ac_po_normal, ier)
  endif

  ! elastic-poroelastic coupling surface
  if (num_coupling_el_po_faces > 0) then
    start(1) = local_dim_coupling_el_po_ispec * myrank
    count_ad(1) = num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_ispec/array", 0, 1, coupling_el_po_ispec, ier)
    call adios_schedule_read(handle, sel, "coupling_po_el_ispec/array", 0, 1, coupling_po_el_ispec, ier)

    start(1) = local_dim_coupling_el_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_ijk/array", 0, 1, coupling_el_po_ijk, ier)
    call adios_schedule_read(handle, sel, "coupling_po_el_ijk/array", 0, 1, coupling_po_el_ijk, ier)

    start(1) = local_dim_coupling_el_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_jacobian2Dw/array", 0, 1, coupling_el_po_jacobian2Dw, ier)

    start(1) = local_dim_coupling_el_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_normal/array", 0, 1, coupling_el_po_normal, ier)
  endif

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    start(1) = local_dim_my_neighbors_ext_mesh * myrank
    count_ad(1) = num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "my_neighbors_ext_mesh/array", 0, 1, my_neighbors_ext_mesh, ier)
    call adios_schedule_read(handle, sel, "nibool_interfaces_ext_mesh/array", 0, 1, nibool_interfaces_ext_mesh, ier)

    start(1) = local_dim_ibool_interfaces_ext_mesh * myrank
    count_ad(1) = max_nibool_interfaces_ext_mesh * num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibool_interfaces_ext_mesh_dummy/array", 0, 1, ibool_interfaces_ext_mesh, ier)
  endif

  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    start(1) = local_dim_c11store * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_aniso
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "c11store/array", 0, 1, c11store, ier)
    call adios_schedule_read(handle, sel, "c12store/array", 0, 1, c12store, ier)
    call adios_schedule_read(handle, sel, "c13store/array", 0, 1, c13store, ier)
    call adios_schedule_read(handle, sel, "c14store/array", 0, 1, c14store, ier)
    call adios_schedule_read(handle, sel, "c15store/array", 0, 1, c15store, ier)
    call adios_schedule_read(handle, sel, "c16store/array", 0, 1, c16store, ier)
    call adios_schedule_read(handle, sel, "c22store/array", 0, 1, c22store, ier)
    call adios_schedule_read(handle, sel, "c23store/array", 0, 1, c23store, ier)
    call adios_schedule_read(handle, sel, "c24store/array", 0, 1, c24store, ier)
    call adios_schedule_read(handle, sel, "c25store/array", 0, 1, c25store, ier)
    call adios_schedule_read(handle, sel, "c26store/array", 0, 1, c26store, ier)
    call adios_schedule_read(handle, sel, "c33store/array", 0, 1, c33store, ier)
    call adios_schedule_read(handle, sel, "c34store/array", 0, 1, c34store, ier)
    call adios_schedule_read(handle, sel, "c35store/array", 0, 1, c35store, ier)
    call adios_schedule_read(handle, sel, "c36store/array", 0, 1, c36store, ier)
    call adios_schedule_read(handle, sel, "c44store/array", 0, 1, c44store, ier)
    call adios_schedule_read(handle, sel, "c45store/array", 0, 1, c45store, ier)
    call adios_schedule_read(handle, sel, "c46store/array", 0, 1, c46store, ier)
    call adios_schedule_read(handle, sel, "c55store/array", 0, 1, c55store, ier)
    call adios_schedule_read(handle, sel, "c56store/array", 0, 1, c56store, ier)
    call adios_schedule_read(handle, sel, "c66store/array", 0, 1, c66store, ier)
  endif

  ! inner / outer elements
  start(1) = local_dim_ispec_is_inner * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ispec_is_inner/array", 0, 1, ispec_is_inner, ier)

  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic > 0) then
      start(1) = local_dim_phase_ispec_inner_acoustic * myrank
      count_ad(1) = num_phase_ispec_acoustic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_acoustic/array", 0, 1, phase_ispec_inner_acoustic, ier)
    endif
  endif

  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic > 0) then
      start(1) = local_dim_phase_ispec_inner_elastic * myrank
      count_ad(1) = num_phase_ispec_elastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_elastic/array", 0, 1, phase_ispec_inner_elastic, ier)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic > 0) then
      start(1) = local_dim_phase_ispec_inner_poroelastic * myrank
      count_ad(1) = num_phase_ispec_poroelastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_poroelastic/array", 0, 1, phase_ispec_inner_poroelastic, ier)
    endif
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! acoustic domain colors
    if (ACOUSTIC_SIMULATION) then
      start(1) = local_dim_num_elem_colors_acoustic * myrank
      count_ad(1) = num_colors_outer_acoustic + num_colors_inner_acoustic
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "num_elem_colors_acoustic/array", 0, 1, num_elem_colors_acoustic, ier)
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      start(1) = local_dim_num_elem_colors_elastic * myrank
      count_ad(1) = num_colors_outer_elastic + num_colors_inner_elastic
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "num_elem_colors_elastic/array", 0, 1, num_elem_colors_elastic, ier)
    endif
  endif

  ! mesh surface arrays
  start(1) = local_dim_ispec_is_surface_external_mesh * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ispec_is_surface_external_mesh/array", 0, 1, ispec_is_surface_external_mesh, ier)

  start(1) = local_dim_iglob_is_surface_external_mesh * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "iglob_is_surface_external_mesh/array", 0, 1, iglob_is_surface_external_mesh, ier)

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  ! checks
  if (sel_num > 256) then
    print *,'Error: sel_num ',sel_num,'too big'
    stop 'Error number of selection exceeds array bounds'
  endif

  ! frees selection
  do isel = 1,sel_num
    sel => selections(isel)
    call adios_selection_delete(sel)
  enddo

  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  !call read_mesh_databases2()
  !call check_mesh_database()

  ! debug
  !call sum_all_i(num_interfaces_ext_mesh,inum)
  !if (myrank == 0) then
  !  write(IMAIN,*) 'number of MPI partition interfaces: ',inum
  !  write(IMAIN,*)
  !endif

  ! MPI communications
  ! acoustic wavefield buffers
  if (ACOUSTIC_SIMULATION) then
    allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1930')
    allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh*NB_RUNS_ACOUSTIC_GPU,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1931')

    allocate(request_send_scalar_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1934')
    allocate(request_recv_scalar_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1935')

    buffer_send_scalar_ext_mesh(:,:) = 0.0_CUSTOM_REAL; buffer_recv_scalar_ext_mesh(:,:) = 0.0_CUSTOM_REAL
    request_send_scalar_ext_mesh(:) = 0; request_recv_scalar_ext_mesh(:) = 0
  endif

  ! elastic wavefield buffers
  if (ELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1928')
    allocate(buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1929')

    allocate(request_send_vector_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1932')
    allocate(request_recv_vector_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1933')

    buffer_send_vector_ext_mesh(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh(:,:,:) = 0.0_CUSTOM_REAL
    request_send_vector_ext_mesh(:) = 0; request_recv_vector_ext_mesh(:) = 0
  endif

  ! poroelastic wavefield buffers
  if (POROELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1936')
    allocate(buffer_recv_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1937')
    allocate(buffer_send_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1938')
    allocate(buffer_recv_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1939')

    allocate(request_send_vector_ext_mesh_s(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1940')
    allocate(request_recv_vector_ext_mesh_s(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1941')
    allocate(request_send_vector_ext_mesh_w(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1942')
    allocate(request_recv_vector_ext_mesh_w(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1943')
    if (ier /= 0) stop 'error allocating array buffer_send_vector_ext_mesh etc.'

    buffer_send_vector_ext_mesh_s(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh_s(:,:,:) = 0.0_CUSTOM_REAL
    buffer_send_vector_ext_mesh_w(:,:,:) = 0.0_CUSTOM_REAL; buffer_recv_vector_ext_mesh_w(:,:,:) = 0.0_CUSTOM_REAL
    request_send_vector_ext_mesh_s(:) = 0; request_recv_vector_ext_mesh_s(:) = 0
    request_send_vector_ext_mesh_w(:) = 0; request_recv_vector_ext_mesh_w(:) = 0
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  done"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases_adios


!-------------------------------------------------------------------------------

!> Reads in moho meshes

  subroutine read_mesh_databases_moho_adios()

  use adios_read_mod
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num,isel
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer(kind=8) :: local_dim_ibelm_moho_bot,  local_dim_ibelm_moho_top, &
    local_dim_ijk_moho_bot,    local_dim_ijk_moho_top, &
    local_dim_normal_moho_bot, local_dim_normal_moho_top, &
    local_dim_is_moho_bot,     local_dim_is_moho_top

  integer :: ier
  integer :: comm

  ! always needed to be allocated for routine arguments
  allocate( is_moho_top(NSPEC_BOUN),is_moho_bot(NSPEC_BOUN),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1944')
  if (ier /= 0) stop 'Error allocating array is_moho_top etc.'

  ! checks if anything to do
  if (ELASTIC_SIMULATION .and. SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then

    !-------------------------------------.
    ! Open ADIOS Database file, read mode |
    !-------------------------------------'
    sel_num = 0

    database_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/moho.bp"

    ! gets MPI communicator
    comm = comm_adios

    call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)
    call adios_read_open_file (handle, database_name, 0, comm, ier)
    if (ier /= 0) call abort_mpi()

    !------------------------------------------------------------------.
    ! Get scalar values. Might be differents for different processors. |
    ! Hence the selection writeblock.                                  |
    ! ONLY NSPEC_AB and NGLOB_AB
    !------------------------------------------------------------------'
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(handle, sel, "nspec2d_moho", 0, 1, NSPEC2D_MOHO, ier)

    call adios_perform_reads(handle, ier)
    if (ier /= 0) call abort_mpi()

    !----------------------------------------------.
    ! Fetch values to compute the simulation type. |
    !----------------------------------------------'
    call adios_get_scalar(handle, "ibelm_moho_bot/local_dim", local_dim_ibelm_moho_bot ,ier)
    call adios_get_scalar(handle, "ibelm_moho_top/local_dim", local_dim_ibelm_moho_top ,ier)

    call adios_get_scalar(handle, "ijk_moho_bot/local_dim", local_dim_ijk_moho_bot ,ier)
    call adios_get_scalar(handle, "ijk_moho_top/local_dim", local_dim_ijk_moho_top ,ier)

    call adios_get_scalar(handle,"normal_moho_bot /local_dim", local_dim_normal_moho_bot ,ier)
    call adios_get_scalar(handle, "normal_moho_top/local_dim", local_dim_normal_moho_top ,ier)

    call adios_get_scalar(handle, "is_moho_bot/local_dim", local_dim_is_moho_bot ,ier)
    call adios_get_scalar(handle, "is_moho_top/local_dim", local_dim_is_moho_top ,ier)

    !---------------------------------------------.
    ! Allocate arrays with previously read values |
    !---------------------------------------------'
    allocate(ibelm_moho_bot(NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1945')
    allocate(ibelm_moho_top(NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1946')
    allocate(normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1947')
    allocate(normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1948')
    allocate(ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1949')
    allocate(ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1950')
    if (ier /= 0) stop 'error allocating array ibelm_moho_bot etc.'

    !-----------------------------------.
    ! Read arrays from external_mesh.bp |
    !-----------------------------------'
    start(1) = local_dim_ibelm_moho_bot * myrank
    count_ad(1) = NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_moho_bot/array", 0, 1, ibelm_moho_bot, ier)

    start(1) = local_dim_ibelm_moho_top * myrank
    count_ad(1) = NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibelm_moho_top/array", 0, 1, ibelm_moho_top, ier)

    start(1) = local_dim_ijk_moho_bot * myrank
    count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ijk_moho_bot/array", 0, 1, ijk_moho_bot, ier)

    start(1) = local_dim_ijk_moho_top * myrank
    count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ijk_moho_top/array", 0, 1, ijk_moho_top, ier)

    start(1) = local_dim_normal_moho_bot * myrank
    count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "normal_moho_bot/array", 0, 1, normal_moho_bot, ier)

    start(1) = local_dim_normal_moho_top * myrank
    count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "normal_moho_top/array", 0, 1, normal_moho_top, ier)

    start(1) = local_dim_is_moho_bot * myrank
    count_ad(1) = NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "is_moho_bot/array", 0, 1, is_moho_bot, ier)

    start(1) = local_dim_is_moho_top * myrank
    count_ad(1) = NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "is_moho_top/array", 0, 1, is_moho_top, ier)

    !---------------------------------------------------------------.
    ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
    !---------------------------------------------------------------'
    call adios_perform_reads(handle, ier)
    if (ier /= 0) call abort_mpi()

    ! frees selection
    do isel = 1,sel_num
      sel => selections(sel_num)
      call adios_selection_delete(sel)
    enddo

    call adios_read_close(handle,ier)
    call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
    if (ier /= 0 ) stop 'Error adios read finalize'

  else
    ! dummy
    NSPEC2D_MOHO = 1
  endif

  ! moho boundary
  if (ELASTIC_SIMULATION) then
    ! always needed to be allocated for routine arguments
    allocate(dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1951')
    allocate(dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1952')
    allocate(b_dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1953')
    allocate(b_dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1954')
    if (ier /= 0) stop 'Error allocating array dsdx_top etc.'
  endif

  end subroutine read_mesh_databases_moho_adios
