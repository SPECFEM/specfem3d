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


  subroutine read_mesh_for_init_ADIOS()

  use constants, only: MAX_STRING_LEN,myrank

  use specfem_par, only: LOCAL_PATH, NSPEC_AB,NGLOB_AB,NSPEC_IRREGULAR

  use adios_helpers_mod
  use manager_adios

  implicit none
  ! Local variables
  character(len=MAX_STRING_LEN) :: database_name

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = get_adios_filename(trim(LOCAL_PATH) // "/external_mesh")

  ! initiate new group
  call init_adios_group(myadios_group,"SolverReaderInit")

  ! setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

  !------------------------------------.
  ! Read variables from the adios file |
  !------------------------------------'
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec", NSPEC_AB)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nglob", NGLOB_AB)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_irregular", NSPEC_IRREGULAR)

  !--------------------------------------------.
  ! Perform the reads and close the adios file |
  !--------------------------------------------'
  call read_adios_perform(myadios_file)

  ! closes adios file & cleans/removes group object
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"SolverReaderInit")

  end subroutine read_mesh_for_init_ADIOS

!==============================================================================

  subroutine read_mesh_databases_adios()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_movie

  use pml_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  character(len=MAX_STRING_LEN) :: database_name

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num,isel
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad
  integer :: ier !,inum

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

  integer :: nspec_ext,nglob_ext,nspec_irregular_ext

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0

  database_name = get_adios_filename(trim(LOCAL_PATH) // "/external_mesh")

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Reading mesh databases...'
    write(IMAIN,*) '  from directory: ',trim(LOCAL_PATH)
    write(IMAIN,*) '  database file : ',trim(database_name)
#if defined(USE_ADIOS)
    write(IMAIN,*) '  using ADIOS1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  using ADIOS2 file format'
#endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initiate new group
  call init_adios_group(myadios_group,"SolverReader")

  ! setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  ! ONLY NSPEC_AB and NGLOB_AB
  !------------------------------------------------------------------'
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec", nspec_ext)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nglob", nglob_ext)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_irregular", nspec_irregular_ext)

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
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ispec_is_acoustic", local_dim_ispec_is_acoustic)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ispec_is_elastic", local_dim_ispec_is_elastic)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ispec_is_poroelastic", local_dim_ispec_is_poroelastic)

  start(1) = local_dim_ispec_is_acoustic * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "ispec_is_acoustic/array", ispec_is_acoustic)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "ispec_is_elastic/array", ispec_is_elastic)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "ispec_is_poroelastic/array", ispec_is_poroelastic)

  ! Perform the read, so we can use the values.
  call read_adios_perform(myadios_file)

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
  ! number of elastic elements in this partition
  nspec_poroelastic = count(ispec_is_poroelastic(:))

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  simulation w/ acoustic    domain: ',ACOUSTIC_SIMULATION
    write(IMAIN,*) '  simulation w/ elastic     domain: ',ELASTIC_SIMULATION
    write(IMAIN,*) '  simulation w/ poroelastic domain: ',POROELASTIC_SIMULATION
    write(IMAIN,*)
    write(IMAIN,*) '  slice 0 has:'
    write(IMAIN,*) '  number of elements acoustic   :',nspec_acoustic
    write(IMAIN,*) '  number of elements elastic    :',nspec_elastic
    write(IMAIN,*) '  number of elements poroelastic:',nspec_poroelastic
    call flush_IMAIN()
  endif

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
  NSPEC_CPML = 0
  if (PML_CONDITIONS) then
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_cpml", nspec_cpml)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "CPML_width_x", CPML_width_x)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "CPML_width_y", CPML_width_y)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "CPML_width_z", CPML_width_z)

    if (nspec_cpml > 0) then
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        call read_adios_scalar(myadios_file, myadios_group, myrank, "nglob_interface_PML_acoustic", nglob_interface_PML_acoustic)
        call read_adios_scalar(myadios_file, myadios_group, myrank, "nglob_interface_PML_elastic", nglob_interface_PML_elastic)
      endif
    endif
  endif

  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_abs_boundary_faces", num_abs_boundary_faces)

  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_xmin", nspec2d_xmin)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_xmax", nspec2d_xmax)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_ymin", nspec2d_ymin)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_ymax", nspec2d_ymax)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_bottom", nspec2d_bottom)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_top", nspec2d_top)


  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_free_surface_faces", num_free_surface_faces)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_coupling_ac_el_faces", num_coupling_ac_el_faces)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_coupling_ac_po_faces", num_coupling_ac_po_faces)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_coupling_el_po_faces", num_coupling_el_po_faces)

  call read_adios_scalar(myadios_file, myadios_group, myrank, "num_interfaces_ext_mesh", num_interfaces_ext_mesh)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "max_nibool_interfaces_ext_mesh", max_nibool_interfaces_ext_mesh)

  if (ACOUSTIC_SIMULATION) then
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_inner_acoustic", nspec_inner_acoustic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_outer_acoustic", nspec_outer_acoustic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "num_phase_ispec_acoustic", num_phase_ispec_acoustic)
  endif

  if (ELASTIC_SIMULATION) then
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_inner_elastic", nspec_inner_elastic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_outer_elastic", nspec_outer_elastic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "num_phase_ispec_elastic", num_phase_ispec_elastic)
  endif

  if (POROELASTIC_SIMULATION) then
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_inner_poroelastic", nspec_inner_poroelastic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec_outer_poroelastic", nspec_outer_poroelastic)
    call read_adios_scalar(myadios_file, myadios_group, myrank, "num_phase_ispec_poroelastic", num_phase_ispec_poroelastic)
  endif

  num_colors_outer_acoustic = 0
  num_colors_inner_acoustic = 0
  num_colors_outer_elastic = 0
  num_colors_inner_elastic = 0
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      call read_adios_scalar(myadios_file, myadios_group, myrank, "num_colors_outer_acoustic", num_colors_outer_acoustic)
      call read_adios_scalar(myadios_file, myadios_group, myrank, "num_colors_outer_acoustic", num_colors_inner_acoustic)
    endif
    if (ELASTIC_SIMULATION) then
      call read_adios_scalar(myadios_file, myadios_group, myrank, "num_colors_outer_elastic", num_colors_outer_elastic)
      call read_adios_scalar(myadios_file, myadios_group, myrank, "num_colors_outer_elastic", num_colors_inner_elastic)
    endif
  endif

  call read_adios_scalar(myadios_file, myadios_group, myrank, "nfaces_surface", nfaces_surface)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  ! note: adios_get_scalar here retrieves the same local_dim for everyone (from writer rank 0)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibool", local_dim_ibool)

  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "x_global", local_dim_x_global)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "y_global", local_dim_y_global)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "z_global", local_dim_z_global)

  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "xixstore", local_dim_xixstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "xiystore", local_dim_xiystore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "xizstore", local_dim_xizstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "etaxstore", local_dim_etaxstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "etaystore", local_dim_etaystore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "etazstore", local_dim_etazstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "gammaxstore", local_dim_gammaxstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "gammaystore", local_dim_gammaystore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "gammazstore", local_dim_gammazstore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "jacobianstore", local_dim_jacobianstore)

  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "kappastore", local_dim_kappastore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "mustore", local_dim_mustore)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rhostore", local_dim_rhostore)

  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                   "irregular_element_number", local_dim_irregular_element_number)

  ! single values
  call read_adios_scalar(myadios_file, myadios_group, myrank, "jacobian_regular", jacobian_regular)
  call read_adios_scalar(myadios_file, myadios_group, myrank, "xix_regular", xix_regular)

  if (ACOUSTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmass_acoustic", local_dim_rmass_acoustic)
  endif
  if (ELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmass", local_dim_rmass)

    if (APPROXIMATE_OCEAN_LOAD) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmass_ocean_load", local_dim_rmass_ocean_load)
    endif
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho_vp", local_dim_rho_vp)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho_vs", local_dim_rho_vs)
  endif
  if (POROELASTIC_SIMULATION) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "rmass_solid_poroelastic", local_dim_rmass_solid_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "rmass_fluid_poroelastic", local_dim_rmass_fluid_poroelastic)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rhoarraystore", local_dim_rhoarraystore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "kappaarraystore", local_dim_kappaarraystore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "permstore", local_dim_permstore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "etastore", local_dim_etastore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "tortstore", local_dim_tortstore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "phistore", local_dim_phistore)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho_vpI", local_dim_rho_vpI)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho_vpII", local_dim_rho_vpII)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho_vsI", local_dim_rho_vsI)
  endif
  if (PML_CONDITIONS) then
    if (nspec_cpml > 0) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "CPML_regions", local_dim_CPML_regions)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "CPML_to_spec", local_dim_CPML_to_spec)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "is_CPML", local_dim_is_CPML)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "d_store_x", local_dim_d_store_x)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "d_store_y", local_dim_d_store_y)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "d_store_z", local_dim_d_store_z)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "k_store_x", local_dim_k_store_x)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "k_store_y", local_dim_k_store_y)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "k_store_z", local_dim_k_store_z)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "alpha_store_x", local_dim_alpha_store_x)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "alpha_store_y", local_dim_alpha_store_y)
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "alpha_store_z", local_dim_alpha_store_z)
      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                           "points_interface_PML_acoustic", local_dim_points_interface_PML_acoustic)
        endif
        if (nglob_interface_PML_elastic > 0) then
          call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                           "points_interface_PML_elastic", local_dim_points_interface_PML_elastic)
        endif
      endif
    endif
  endif

  ! absorbing boundaries
  if (num_abs_boundary_faces > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "abs_boundary_ispec", local_dim_abs_boundary_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "abs_boundary_ijk", local_dim_abs_boundary_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "abs_boundary_jacobian2Dw", local_dim_abs_boundary_jacobian2Dw)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "abs_boundary_normal", local_dim_abs_boundary_normal)
    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! mass matrix contributions
      if (ELASTIC_SIMULATION) then
        call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmassx", local_dim_rmassx)
        call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmassy", local_dim_rmassy)
        call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmassz", local_dim_rmassz)
      endif
      if (ACOUSTIC_SIMULATION) then
        call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rmassz_acoustic", local_dim_rmassz_acoustic)
      endif
    endif
  endif

  if (nspec2d_xmin > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_xmin", local_dim_ibelm_xmin)
  endif
  if (nspec2d_xmax > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_xmax", local_dim_ibelm_xmax)
  endif
  if (nspec2d_ymin > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_ymin", local_dim_ibelm_ymin)
  endif
  if (nspec2d_ymax > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_ymax", local_dim_ibelm_ymax)
  endif
  if (nspec2d_bottom > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_bottom", local_dim_ibelm_bottom)
  endif
  if (nspec2d_top > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_top", local_dim_ibelm_top)
  endif

  if (num_free_surface_faces > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "free_surface_ispec", local_dim_free_surface_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "free_surface_ijk", local_dim_free_surface_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "free_surface_jacobian2Dw", local_dim_free_surface_jacobian2Dw)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "free_surface_normal", local_dim_free_surface_normal)
  endif
  if (num_coupling_ac_el_faces > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_el_ispec", local_dim_coupling_ac_el_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_el_ijk", local_dim_coupling_ac_el_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_el_jacobian2Dw", local_dim_coupling_ac_el_jacobian2Dw)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_el_normal", local_dim_coupling_ac_el_normal)
  endif
  if (num_coupling_ac_po_faces > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_po_ispec", local_dim_coupling_ac_po_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_po_ijk", local_dim_coupling_ac_po_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_po_jacobian2Dw", local_dim_coupling_ac_po_jacobian2Dw)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_ac_po_normal", local_dim_coupling_ac_po_normal)
  endif
  if (num_coupling_el_po_faces > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_el_po_ispec", local_dim_coupling_el_po_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_po_el_ispec", local_dim_coupling_po_el_ispec)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_el_po_ijk", local_dim_coupling_el_po_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_po_el_ijk", local_dim_coupling_po_el_ijk)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_el_po_jacobian2Dw", local_dim_coupling_el_po_jacobian2Dw)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "coupling_el_po_normal", local_dim_coupling_el_po_normal)
  endif
  if (num_interfaces_ext_mesh > 0) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "my_neighbors_ext_mesh", local_dim_my_neighbors_ext_mesh)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "nibool_interfaces_ext_mesh", local_dim_nibool_interfaces_ext_mesh)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                     "ibool_interfaces_ext_mesh_dummy", local_dim_ibool_interfaces_ext_mesh)
  endif
  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "c11store", local_dim_c11store)
  endif
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ispec_is_inner", local_dim_ispec_is_inner)
  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic > 0) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                       "phase_ispec_inner_acoustic", local_dim_phase_ispec_inner_acoustic)
    endif
  endif
  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic > 0) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                       "phase_ispec_inner_elastic", local_dim_phase_ispec_inner_elastic)
    endif
  endif
  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic > 0) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                       "phase_ispec_inner_poroelastic", local_dim_phase_ispec_inner_poroelastic)
    endif
  endif
  if (USE_MESH_COLORING_GPU) then
    if (ACOUSTIC_SIMULATION) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                       "num_elem_colors_acoustic", local_dim_num_elem_colors_acoustic)
    endif
    if (ELASTIC_SIMULATION) then
      call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                       "num_elem_colors_elastic", local_dim_num_elem_colors_elastic)
    endif
  endif
  ! for mesh surface
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                   "ispec_is_surface_external_mesh", local_dim_ispec_is_surface_external_mesh)
  call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, &
                                   "iglob_is_surface_external_mesh", local_dim_iglob_is_surface_external_mesh)

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
    potential_acoustic(:) = 0.0_CUSTOM_REAL; potential_dot_acoustic(:) = 0.0_CUSTOM_REAL
    potential_dot_dot_acoustic(:) = 0.0_CUSTOM_REAL

    !if (SIMULATION_TYPE /= 1) then
    !  allocate(potential_acoustic_adj_coupling(NGLOB_AB),stat=ier) ! not used yet
    !  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1794')
    !  if (ier /= 0) stop 'error allocating array potential_acoustic_adj_coupling'
    !endif
    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1795')
    if (ier /= 0) stop 'error allocating array rmass_acoustic'
    rmass_acoustic(:) = 0.0_CUSTOM_REAL

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
    displ(:,:) = 0.0_CUSTOM_REAL; veloc(:,:) = 0.0_CUSTOM_REAL; accel(:,:) = 0.0_CUSTOM_REAL

    if (SIMULATION_TYPE /= 1) then
      allocate(accel_adj_coupling(NDIM,NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1801')
      if (ier /= 0) stop 'error allocating array accel_adj_coupling'
      accel_adj_coupling(:,:) = 0.0_CUSTOM_REAL
    endif

    ! allocates mass matrix
    allocate(rmass(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1802')
    rmass(:) = 0.0_CUSTOM_REAL

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
    rho_vp(:,:,:,:) = 0.0_CUSTOM_REAL; rho_vs(:,:,:,:) = 0.0_CUSTOM_REAL

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
    c11store(:,:,:,:) = 0.0_CUSTOM_REAL; c12store(:,:,:,:) = 0.0_CUSTOM_REAL; c13store(:,:,:,:) = 0.0_CUSTOM_REAL
    c14store(:,:,:,:) = 0.0_CUSTOM_REAL; c15store(:,:,:,:) = 0.0_CUSTOM_REAL; c16store(:,:,:,:) = 0.0_CUSTOM_REAL
    c22store(:,:,:,:) = 0.0_CUSTOM_REAL; c23store(:,:,:,:) = 0.0_CUSTOM_REAL; c24store(:,:,:,:) = 0.0_CUSTOM_REAL
    c25store(:,:,:,:) = 0.0_CUSTOM_REAL; c26store(:,:,:,:) = 0.0_CUSTOM_REAL; c33store(:,:,:,:) = 0.0_CUSTOM_REAL
    c34store(:,:,:,:) = 0.0_CUSTOM_REAL; c35store(:,:,:,:) = 0.0_CUSTOM_REAL; c36store(:,:,:,:) = 0.0_CUSTOM_REAL
    c44store(:,:,:,:) = 0.0_CUSTOM_REAL; c45store(:,:,:,:) = 0.0_CUSTOM_REAL; c46store(:,:,:,:) = 0.0_CUSTOM_REAL
    c55store(:,:,:,:) = 0.0_CUSTOM_REAL; c56store(:,:,:,:) = 0.0_CUSTOM_REAL; c66store(:,:,:,:) = 0.0_CUSTOM_REAL

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
    R_xx(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_yy(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_xy(:,:,:,:,:) = 0.0_CUSTOM_REAL
    R_xz(:,:,:,:,:) = 0.0_CUSTOM_REAL; R_yz(:,:,:,:,:) = 0.0_CUSTOM_REAL

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
    epsilondev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilondev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilondev_yz(:,:,:,:) = 0.0_CUSTOM_REAL; epsilondev_trace(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(R_trace(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1840')
    if (ier /= 0) stop 'error allocating array R_trace etc.'
    R_trace(:,:,:,:,:) = 0.0_CUSTOM_REAL

    ! note: needed for some subroutine arguments
    allocate(epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1841')
    if (ier /= 0) stop 'error allocating array epsilon_trace_over_3'
    epsilon_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL

    ! needed for attenuation
    allocate(factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1842')
    if (ier /= 0) stop 'error allocating array factor_common'
    factor_common(:,:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1843')
    if (ier /= 0) stop 'error allocating array factor_common_kappa'
    factor_common_kappa(:,:,:,:,:) = 0.0_CUSTOM_REAL

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1844')
      if (ier /= 0) stop 'error allocating array rmass_ocean_load'
      rmass_ocean_load(:) = 0.0_CUSTOM_REAL
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
    displs_poroelastic(:,:) = 0.0_CUSTOM_REAL; velocs_poroelastic(:,:) = 0.0_CUSTOM_REAL
    accels_poroelastic(:,:) = 0.0_CUSTOM_REAL
    displw_poroelastic(:,:) = 0.0_CUSTOM_REAL; velocw_poroelastic(:,:) = 0.0_CUSTOM_REAL
    accelw_poroelastic(:,:) = 0.0_CUSTOM_REAL

    allocate(rmass_solid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1852')
    if (ier /= 0) stop 'error allocating array rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1853')
    if (ier /= 0) stop 'error allocating array rmass_fluid_poroelastic'
    rmass_solid_poroelastic(:) = 0.0_CUSTOM_REAL; rmass_fluid_poroelastic(:) = 0.0_CUSTOM_REAL

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
    rhoarraystore(:,:,:,:,:) = 0.0_CUSTOM_REAL; kappaarraystore(:,:,:,:,:) = 0.0_CUSTOM_REAL
    etastore(:,:,:,:) = 0.0_CUSTOM_REAL; tortstore(:,:,:,:) = 0.0_CUSTOM_REAL
    phistore(:,:,:,:) = 0.0_CUSTOM_REAL; permstore(:,:,:,:,:) = 0.0_CUSTOM_REAL
    rho_vpI(:,:,:,:) = 0.0_CUSTOM_REAL; rho_vpII(:,:,:,:) = 0.0_CUSTOM_REAL
    rho_vsI(:,:,:,:) = 0.0_CUSTOM_REAL

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
    epsilonsdev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonsdev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonsdev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonsdev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonsdev_yz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_xx(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonwdev_yy(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_xy(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonwdev_xz(:,:,:,:) = 0.0_CUSTOM_REAL
    epsilonwdev_yz(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate(epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1873')
    allocate(epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1874')
    if (ier /= 0) stop 'error allocating array epsilons_trace_over_3 etc.'
    epsilons_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL; epsilonw_trace_over_3(:,:,:,:) = 0.0_CUSTOM_REAL

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
      CPML_regions(:) = 0; CPML_to_spec(:) = 0

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
      d_store_x(:,:,:,:) = 0.0_CUSTOM_REAL; d_store_y(:,:,:,:) = 0.0_CUSTOM_REAL; d_store_z(:,:,:,:) = 0.0_CUSTOM_REAL
      K_store_x(:,:,:,:) = 0.0_CUSTOM_REAL; K_store_y(:,:,:,:) = 0.0_CUSTOM_REAL; K_store_z(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_x(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_y(:,:,:,:) = 0.0_CUSTOM_REAL
      alpha_store_z(:,:,:,:) = 0.0_CUSTOM_REAL

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1887')
          if (ier /= 0) stop 'error allocating array points_interface_PML_acoustic'
          points_interface_PML_acoustic(:) = 0
        endif
        if (nglob_interface_PML_elastic > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 1888')
          if (ier /= 0) stop 'error allocating array points_interface_PML_elastic'
          points_interface_PML_elastic(:) = 0
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
  abs_boundary_ispec(:) = 0; abs_boundary_ijk(:,:,:) = 0
  abs_boundary_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL; abs_boundary_normal(:,:,:) = 0.0_CUSTOM_REAL

  allocate(ibelm_xmin(nspec2D_xmin),ibelm_xmax(nspec2D_xmax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1893')
  allocate(ibelm_ymin(nspec2D_ymin),ibelm_ymax(nspec2D_ymax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1894')
  allocate(ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1895')
  if (ier /= 0) stop 'error allocating arrays ibelm_xmin,ibelm_xmax etc.'
  ibelm_xmin(:) = 0; ibelm_xmax(:) = 0
  ibelm_ymin(:) = 0; ibelm_ymax(:) = 0
  ibelm_bottom(:) = 0; ibelm_top(:) = 0

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
  free_surface_ispec(:) = 0; free_surface_ijk(:,:,:) = 0
  free_surface_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL; free_surface_normal(:,:,:) = 0.0_CUSTOM_REAL

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
  coupling_ac_el_ispec(:) = 0; coupling_ac_el_ijk(:,:,:) = 0
  coupling_ac_el_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_ac_el_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

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
  coupling_ac_po_ispec(:) = 0; coupling_ac_po_ijk(:,:,:) = 0
  coupling_ac_po_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_ac_po_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

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
  coupling_el_po_ispec(:) = 0; coupling_el_po_ijk(:,:,:) = 0
  coupling_po_el_ispec(:) = 0; coupling_po_el_ijk(:,:,:) = 0
  coupling_el_po_normal(:,:,:) = 0.0_CUSTOM_REAL; coupling_el_po_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1914')
    allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1915')
    if (ier /= 0) stop 'error allocating array my_neighbors_ext_mesh etc.'
    my_neighbors_ext_mesh(:) = -1; nibool_interfaces_ext_mesh(:) = 0

    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1916')
    if (ier /= 0) stop 'error allocating array ibool_interfaces_ext_mesh'
    ibool_interfaces_ext_mesh(:,:) = 0
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
  ispec_is_inner(:) = .false.

  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic < 0) stop 'error acoustic simulation:' // &
                                    'num_phase_ispec_acoustic is < zero'
    allocate(phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1919')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0
  endif

  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic < 0) stop 'error elastic simulation:' // &
                                   'num_phase_ispec_elastic is < zero'
    allocate(phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1920')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0
  endif

  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic < 0) &
      stop 'error poroelastic simulation:num_phase_ispec_poroelastic is < zero'
    allocate(phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1921')
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_poroelastic'
    phase_ispec_inner_poroelastic(:,:) = 0
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! acoustic domain colors
    if (ACOUSTIC_SIMULATION) then
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1922')
      if (ier /= 0) stop 'error allocating num_elem_colors_acoustic array'
      num_elem_colors_acoustic(:) = 0
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1923')
      if (ier /= 0) stop 'error allocating num_elem_colors_elastic array'
      num_elem_colors_elastic(:) = 0
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
  ispec_is_surface_external_mesh(:) = .false.; iglob_is_surface_external_mesh(:) = .false.

  !-----------------------------------.
  ! Read arrays from external_mesh.bp |
  !-----------------------------------'
  start(1) = local_dim_ibool * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibool/array", ibool)

  start(1) = local_dim_x_global * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "x_global/array", xstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "y_global/array", ystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "z_global/array", zstore)

  start(1) = local_dim_irregular_element_number * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "irregular_element_number/array", irregular_element_number)

  start(1) = local_dim_xixstore * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_IRREGULAR
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "xixstore/array", xixstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "xiystore/array", xiystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "xizstore/array", xizstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "etaxstore/array", etaxstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "etaystore/array", etaystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "etazstore/array", etazstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "gammaxstore/array", gammaxstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "gammaystore/array", gammaystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "gammazstore/array", gammazstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "jacobianstore/array", jacobianstore)

  start(1) = local_dim_kappastore * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "kappastore/array", kappastore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "mustore/array", mustore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rhostore/array", rhostore)

  ! perform reading
  call read_adios_perform(myadios_file)

  if (ACOUSTIC_SIMULATION) then
    start(1) = local_dim_rmass_acoustic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmass_acoustic/array", rmass_acoustic)
  endif

  if (ELASTIC_SIMULATION) then
    start(1) = local_dim_rmass * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmass/array", rmass)

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      start(1) = local_dim_rmass_ocean_load * myrank
      count_ad(1) = NGLOB_AB !nglob_ocean
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmass_ocean_load/array", rmass_ocean_load)
    endif

    !pll material parameters for stacey conditions
    start(1) = local_dim_rho_vp * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rho_vp/array", rho_vp)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rho_vs/array", rho_vs)
  endif

  if (POROELASTIC_SIMULATION) then
    start(1) = local_dim_rmass_solid_poroelastic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "rmass_solid_poroelastic/array", rmass_solid_poroelastic)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "rmass_fluid_poroelastic/array", rmass_fluid_poroelastic)

    start(1) = local_dim_rhoarraystore * myrank
    count_ad(1) = 2 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rhoarraystore/array", rhoarraystore)

    start(1) = local_dim_kappaarraystore* myrank
    count_ad(1) = 3 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "kappaarraystore/array", kappaarraystore)

    start(1) = local_dim_permstore * myrank
    count_ad(1) =  6 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "permstore/array", permstore)

    start(1) = local_dim_etastore * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "etastore/array", etastore)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "tortstore/array", tortstore)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "phistore/array", phistore)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rho_vpI/array", rho_vpI)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rho_vpII/array", rho_vpII)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rho_vsI/array", rho_vsI)
  endif

  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    if (NSPEC_CPML > 0) then
      start(1) = local_dim_CPML_regions * myrank
      count_ad(1) = nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "CPML_regions/array", CPML_regions)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "CPML_to_spec/array", CPML_to_spec)

      start(1) = local_dim_is_cpml * myrank
      count_ad(1) = NSPEC_AB
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "is_CPML/array", is_CPML)

      start(1) = local_dim_d_store_x * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "d_store_x/array", d_store_x)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "d_store_y/array", d_store_y)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "d_store_z/array", d_store_z)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "k_store_x/array", k_store_x)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "k_store_y/array", k_store_y)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "k_store_z/array", k_store_z)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "alpha_store_x/array", alpha_store_x)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "alpha_store_y/array", alpha_store_y)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "alpha_store_z/array", alpha_store_z)

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (nglob_interface_PML_acoustic > 0) then
          start(1) = local_dim_points_interface_PML_acoustic* myrank
          count_ad(1) = nglob_interface_PML_acoustic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call set_selection_boundingbox(sel, start, count_ad)
          call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                         "points_interface_PML_acoustic/array", points_interface_PML_acoustic)
        endif
        if (nglob_interface_PML_elastic > 0) then
          start(1) = local_dim_points_interface_PML_elastic* myrank
          count_ad(1) = nglob_interface_PML_elastic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call set_selection_boundingbox(sel, start, count_ad)
          call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                         "points_interface_PML_elastic/array", points_interface_PML_elastic)
        endif
      endif
    endif
  endif

  if (num_abs_boundary_faces > 0) then
    start(1) = local_dim_abs_boundary_ispec * myrank
    count_ad(1) = num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "abs_boundary_ispec/array", abs_boundary_ispec)

    start(1) = local_dim_abs_boundary_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "abs_boundary_ijk/array", abs_boundary_ijk)

    start(1) = local_dim_abs_boundary_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "abs_boundary_jacobian2Dw/array", abs_boundary_jacobian2Dw)

    start(1) = local_dim_abs_boundary_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_abs_boundary_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "abs_boundary_normal/array", abs_boundary_normal)

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        start(1) = local_dim_rmassx * myrank
        count_ad(1) = NGLOB_AB  ! == nglob_xy in generate_databse
        sel_num = sel_num+1
        sel => selections(sel_num)
        call set_selection_boundingbox(sel, start, count_ad)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmassx/array", rmassx)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmassy/array", rmassy)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "rmassz/array", rmassz)
      endif
      if (ACOUSTIC_SIMULATION) then
        start(1) = local_dim_rmassz_acoustic * myrank
        count_ad(1) = NGLOB_AB ! == nglob_xy in generate_databse
        sel_num = sel_num+1
        sel => selections(sel_num)
        call set_selection_boundingbox(sel, start, count_ad)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                       "rmassz_acoustic/array", rmassz_acoustic)
      endif
    endif
  endif

  if (nspec2d_xmin > 0) then
    start(1) = local_dim_ibelm_xmin * myrank
    count_ad(1) = nspec2D_xmin
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_xmin/array", ibelm_xmin)
  endif
  if (nspec2d_xmax > 0) then
    start(1) = local_dim_ibelm_xmax * myrank
    count_ad(1) = nspec2D_xmax
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_xmax/array", ibelm_xmax)
  endif
  if (nspec2d_ymin > 0) then
    start(1) = local_dim_ibelm_ymin * myrank
    count_ad(1) = nspec2D_ymin
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_ymin/array", ibelm_ymin)
  endif
  if (nspec2d_ymax > 0) then
    start(1) = local_dim_ibelm_ymax * myrank
    count_ad(1) = nspec2D_ymax
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_ymax/array", ibelm_ymax)
  endif
  if (nspec2d_bottom > 0) then
    start(1) = local_dim_ibelm_bottom * myrank
    count_ad(1) = nspec2D_bottom
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_bottom/array", ibelm_bottom)
  endif
  if (nspec2d_top > 0) then
    start(1) = local_dim_ibelm_top * myrank
    count_ad(1) = nspec2D_top
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ibelm_top/array", ibelm_top)
  endif

  ! free surface
  if (num_free_surface_faces > 0) then
    start(1) = local_dim_free_surface_ispec * myrank
    count_ad(1) = num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "free_surface_ispec/array", free_surface_ispec)

    start(1) = local_dim_free_surface_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "free_surface_ijk/array", free_surface_ijk)

    start(1) = local_dim_free_surface_ijk* myrank
    count_ad(1) = NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "free_surface_ijk/array", free_surface_ijk)

    start(1) = local_dim_free_surface_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "free_surface_normal/array", free_surface_normal)
  endif

  ! acoustic-elastic coupling surface
  if (num_coupling_ac_el_faces > 0) then
    start(1) = local_dim_coupling_ac_el_ispec * myrank
    count_ad(1) = num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_el_ispec/array", coupling_ac_el_ispec)

    start(1) = local_dim_coupling_ac_el_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_el_ijk/array", coupling_ac_el_ijk)

    start(1) = local_dim_coupling_ac_el_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_el_jacobian2Dw/array", coupling_ac_el_jacobian2Dw)
    start(1) = local_dim_coupling_ac_el_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_el_normal/array", coupling_ac_el_normal)
  endif

  ! acoustic-poroelastic coupling surface
  if (num_coupling_ac_po_faces > 0) then
    start(1) = local_dim_coupling_ac_po_ispec * myrank
    count_ad(1) = num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_po_ispec/array", coupling_ac_po_ispec)

    start(1) = local_dim_coupling_ac_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_po_ijk/array", coupling_ac_po_ijk)

    start(1) = local_dim_coupling_ac_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_po_jacobian2Dw/array", coupling_ac_po_jacobian2Dw)

    start(1) = local_dim_coupling_ac_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_ac_po_normal/array", coupling_ac_po_normal)
  endif

  ! elastic-poroelastic coupling surface
  if (num_coupling_el_po_faces > 0) then
    start(1) = local_dim_coupling_el_po_ispec * myrank
    count_ad(1) = num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_el_po_ispec/array", coupling_el_po_ispec)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_po_el_ispec/array", coupling_po_el_ispec)

    start(1) = local_dim_coupling_el_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_el_po_ijk/array", coupling_el_po_ijk)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_po_el_ijk/array", coupling_po_el_ijk)

    start(1) = local_dim_coupling_el_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_el_po_jacobian2Dw/array", coupling_el_po_jacobian2Dw)

    start(1) = local_dim_coupling_el_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "coupling_el_po_normal/array", coupling_el_po_normal)
  endif

  ! MPI interfaces
  if (num_interfaces_ext_mesh > 0) then
    start(1) = local_dim_my_neighbors_ext_mesh * myrank
    count_ad(1) = num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "my_neighbors_ext_mesh/array", my_neighbors_ext_mesh)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "nibool_interfaces_ext_mesh/array", nibool_interfaces_ext_mesh)

    start(1) = local_dim_ibool_interfaces_ext_mesh * myrank
    count_ad(1) = max_nibool_interfaces_ext_mesh * num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "ibool_interfaces_ext_mesh_dummy/array", ibool_interfaces_ext_mesh)
  endif

  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    start(1) = local_dim_c11store * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_aniso
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c11store/array", c11store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c12store/array", c12store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c13store/array", c13store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c14store/array", c14store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c15store/array", c15store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c16store/array", c16store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c22store/array", c22store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c23store/array", c23store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c24store/array", c24store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c25store/array", c25store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c26store/array", c26store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c33store/array", c33store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c34store/array", c34store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c35store/array", c35store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c36store/array", c36store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c44store/array", c44store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c45store/array", c45store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c46store/array", c46store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c55store/array", c55store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c56store/array", c56store)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "c66store/array", c66store)
  endif

  ! inner / outer elements
  start(1) = local_dim_ispec_is_inner * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, "ispec_is_inner/array", ispec_is_inner)

  if (ACOUSTIC_SIMULATION) then
    if (num_phase_ispec_acoustic > 0) then
      start(1) = local_dim_phase_ispec_inner_acoustic * myrank
      count_ad(1) = num_phase_ispec_acoustic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                     "phase_ispec_inner_acoustic/array", phase_ispec_inner_acoustic)
    endif
  endif

  if (ELASTIC_SIMULATION) then
    if (num_phase_ispec_elastic > 0) then
      start(1) = local_dim_phase_ispec_inner_elastic * myrank
      count_ad(1) = num_phase_ispec_elastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                     "phase_ispec_inner_elastic/array", phase_ispec_inner_elastic)
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    if (num_phase_ispec_poroelastic > 0) then
      start(1) = local_dim_phase_ispec_inner_poroelastic * myrank
      count_ad(1) = num_phase_ispec_poroelastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                     "phase_ispec_inner_poroelastic/array", phase_ispec_inner_poroelastic)
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
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                     "num_elem_colors_acoustic/array", num_elem_colors_acoustic)
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      start(1) = local_dim_num_elem_colors_elastic * myrank
      count_ad(1) = num_colors_outer_elastic + num_colors_inner_elastic
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count_ad)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                     "num_elem_colors_elastic/array", num_elem_colors_elastic)
    endif
  endif

  ! mesh surface arrays
  start(1) = local_dim_ispec_is_surface_external_mesh * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "ispec_is_surface_external_mesh/array", ispec_is_surface_external_mesh)

  start(1) = local_dim_iglob_is_surface_external_mesh * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count_ad)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                 "iglob_is_surface_external_mesh/array", iglob_is_surface_external_mesh)

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call read_adios_perform(myadios_file)

  ! checks
  if (sel_num > 256) then
    print *,'Error: sel_num ',sel_num,'too big'
    stop 'Error number of selection exceeds array bounds'
  endif

  ! frees selection
  do isel = 1,sel_num
    sel => selections(isel)
    call delete_adios_selection(sel)
  enddo

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"SolverReader")

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

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  use adios_helpers_mod
  use manager_adios

  implicit none

  character(len=MAX_STRING_LEN) :: database_name

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num,isel
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer(kind=8) :: local_dim_ibelm_moho_bot,  local_dim_ibelm_moho_top, &
    local_dim_ijk_moho_bot,    local_dim_ijk_moho_top, &
    local_dim_normal_moho_bot, local_dim_normal_moho_top, &
    local_dim_is_moho_bot,     local_dim_is_moho_top

  integer :: ier

  ! always needed to be allocated for routine arguments
  allocate( is_moho_top(NSPEC_BOUN),is_moho_bot(NSPEC_BOUN),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1944')
  if (ier /= 0) stop 'Error allocating array is_moho_top etc.'
  is_moho_top(:) = .false.; is_moho_bot(:) = .false.

  ! checks if anything to do
  if (ELASTIC_SIMULATION .and. SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then

    !-------------------------------------.
    ! Open ADIOS Database file, read mode |
    !-------------------------------------'
    sel_num = 0

    database_name = get_adios_filename(trim(LOCAL_PATH) // "/moho")

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  reading moho database file: ',trim(database_name)
#if defined(USE_ADIOS)
      write(IMAIN,*) '  using ADIOS1 file format'
#elif defined(USE_ADIOS2)
      write(IMAIN,*) '  using ADIOS2 file format'
#endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! initiate new group
    call init_adios_group(myadios_group,"SolverReaderMoho")

    ! setup the ADIOS library to read the file
    call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

    !------------------------------------------------------------------.
    ! Get scalar values. Might be differents for different processors. |
    ! Hence the selection writeblock.                                  |
    ! ONLY NSPEC_AB and NGLOB_AB
    !------------------------------------------------------------------'
    call read_adios_scalar(myadios_file, myadios_group, myrank, "nspec2d_moho", NSPEC2D_MOHO)

    !----------------------------------------------.
    ! Fetch values to compute the simulation type. |
    !----------------------------------------------'
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_moho_bot", local_dim_ibelm_moho_bot)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ibelm_moho_top", local_dim_ibelm_moho_top)

    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ijk_moho_bot", local_dim_ijk_moho_bot)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "ijk_moho_top", local_dim_ijk_moho_top)

    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "normal_moho_bot", local_dim_normal_moho_bot)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "normal_moho_top", local_dim_normal_moho_top)

    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "is_moho_bot", local_dim_is_moho_bot)
    call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "is_moho_top", local_dim_is_moho_top)

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
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "ibelm_moho_bot/array", ibelm_moho_bot)

    start(1) = local_dim_ibelm_moho_top * myrank
    count_ad(1) = NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "ibelm_moho_top/array", ibelm_moho_top)

    start(1) = local_dim_ijk_moho_bot * myrank
    count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "ijk_moho_bot/array", ijk_moho_bot)

    start(1) = local_dim_ijk_moho_top * myrank
    count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "ijk_moho_top/array", ijk_moho_top)

    start(1) = local_dim_normal_moho_bot * myrank
    count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "normal_moho_bot/array", normal_moho_bot)

    start(1) = local_dim_normal_moho_top * myrank
    count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "normal_moho_top/array", normal_moho_top)

    start(1) = local_dim_is_moho_bot * myrank
    count_ad(1) = NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "is_moho_bot/array", is_moho_bot)

    start(1) = local_dim_is_moho_top * myrank
    count_ad(1) = NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count_ad)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count_ad, &
                                   "is_moho_top/array", is_moho_top)

    !---------------------------------------------------------------.
    ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
    !---------------------------------------------------------------'
    call read_adios_perform(myadios_file)

    ! frees selection
    do isel = 1,sel_num
      sel => selections(sel_num)
      call delete_adios_selection(sel)
    enddo

    ! closes default file and finalizes read method
    call close_file_adios_read_and_finalize_method(myadios_file)
    call delete_adios_group(myadios_group,"SolverReaderMoho")
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
