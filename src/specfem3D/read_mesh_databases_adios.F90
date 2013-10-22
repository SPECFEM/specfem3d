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
! United States and French Government Sponsorship Acknowledged.

!==============================================================================
subroutine read_mesh_for_init(nspec, nglob)

  use mpi
  use adios_read_mod
  use specfem_par, only : myrank, LOCAL_PATH

  implicit none
  ! Paramters
  integer, intent(inout) :: nspec, nglob
  ! Local variables
  character(len=256) :: database_name
  integer(kind=8) :: handle, sel
  integer         :: ier

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = adjustl(LOCAL_PATH)
  database_name = database_name(1:len_trim(database_name)) // "/external_mesh.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                               "verbose=1", ier)
  call adios_read_open_file (handle, database_name, 0, MPI_COMM_WORLD, ier)

  !------------------------------------.
  ! Read variables from the adios file |
  !------------------------------------'
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "/nspec", 0, 1, nspec, ier)
  call adios_schedule_read(handle, sel, "/nglob", 0, 1, nglob, ier)

  !--------------------------------------------.
  ! Perform the reads and close the adios file |
  !--------------------------------------------'
  call adios_perform_reads(handle, ier)
  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)

end subroutine read_mesh_for_init

!==============================================================================
subroutine read_mesh_databases_adios()

  use mpi
  use adios_read_mod

  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  integer :: ier,inum

  character(len=256) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer :: local_dim_ibool, local_dim_x_global, local_dim_y_global,          &
             local_dim_z_global, local_dim_xixstore, local_dim_xiystore,       &
             local_dim_xizstore, local_dim_etaxstore, local_dim_etaystore,     &
             local_dim_etazstore, local_dim_gammaxstore,                       &
             local_dim_gammaystore, local_dim_gammazstore,                     &
             local_dim_jacobianstore, local_dim_kappastore,                    &
             local_dim_mustore, local_dim_rhostore,                            &
             local_dim_ispec_is_acoustic, local_dim_ispec_is_elastic,          &
             local_dim_ispec_is_poroelastic, local_dim_rmass,                  &
             local_dim_rmass_ocean_load, local_dim_rmass_acoustic,             &
             local_dim_rmass_elastic,local_dim_rho_vp,                         &
             local_dim_rho_vs, local_dim_abs_boundary_ispec,                   &
             local_dim_abs_boundary_ijk, local_dim_abs_boundary_jacobian2Dw,   &
             local_dim_abs_boundary_normal, local_dim_ibelm_xmin,              &
             local_dim_ibelm_ymin, local_dim_ibelm_bottom,                     &
             local_dim_ibelm_top, local_dim_free_surface_ispec,                &
             local_dim_free_surface_ijk, local_dim_free_surface_jacobian2Dw,   &
             local_dim_free_surface_normal, local_dim_coupling_ac_el_ispec,    &
             local_dim_coupling_ac_el_ijk,                                     &
             local_dim_coupling_ac_el_jacobian2Dw,                             &
             local_dim_coupling_ac_el_normal, local_dim_my_neighbours_ext_mesh,&
             local_dim_nibool_interfaces_ext_mesh,                             &
             local_dim_ibool_interfaces_ext_mesh,                              &
             local_dim_ispec_is_inner, local_dim_phase_ispec_inner_acoustic,   &
             local_dim_phase_ispec_inner_elastic, local_dim_ibelm_xmax,        &
             local_dim_ibelm_ymax, local_dim_rmass_solid_poroelastic,          &
             local_dim_rmass_fluid_poroelastic, local_dim_rhoarraystore,       &
             local_dim_kappaarraystore, local_dim_permstore,                   &
             local_dim_etastore, local_dim_tortstore, local_dim_phistore,      &
             local_dim_rho_vpI, local_dim_rho_vpII, local_dim_rho_vsI,         &
             local_dim_CPML_regions, local_dim_CPML_to_spec, local_dim_is_CPML,&
             local_dim_d_store_x, local_dim_d_store_y, local_dim_d_store_z,    &
             local_dim_k_store_x, local_dim_k_store_y, local_dim_k_store_z,    &
             local_dim_alpha_store, local_dim_points_interface_PML_acoustic,   &
             local_dim_points_interface_PML_elastic,                           &
             local_dim_rmassx, local_dim_rmassy, local_dim_rmassz,             &
             local_dim_rmassz_acoustic, local_dim_coupling_el_po_ispec,        &
             local_dim_coupling_po_el_ispec, local_dim_coupling_el_po_ijk,     &
             local_dim_coupling_po_el_ijk,                                     &
             local_dim_coupling_el_po_jacobian2Dw,                             &
             local_dim_coupling_el_po_normal, local_dim_c11store,              &
             local_dim_phase_ispec_inner_poroelastic,                          &
             local_dim_num_elem_colors_acoustic,                               &
             local_dim_num_elem_colors_elastic, local_dim_coupling_ac_po_ispec,&
             local_dim_coupling_ac_po_ijk,                                     &
             local_dim_coupling_ac_po_jacobian2Dw,                             &
             local_dim_coupling_ac_po_normal

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0

  database_name = adjustl(LOCAL_PATH)
  database_name = database_name(1:len_trim(database_name)) // "/external_mesh.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                               "verbose=1", ier)
  call adios_read_open_file (handle, database_name, 0, MPI_COMM_WORLD, ier)

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  ! ONLY NSPEC_AB and NGLOB_AB
  !------------------------------------------------------------------'
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "/nspec", 0, 1, NSPEC_AB, ier)
  call adios_schedule_read(handle, sel, "/nglob", 0, 1, NGLOB_AB, ier)
  call adios_perform_reads(handle, ier)

  !----------------------------------------------.
  ! Fetch values to compute the simulation type. |
  !----------------------------------------------'
  sel_num = 0
  call adios_get_scalar(handle, "ispec_is_acoustic/local_dim",&
                        local_dim_ispec_is_acoustic,ier)
  call adios_get_scalar(handle, "ispec_is_elastic/local_dim",&
                        local_dim_ispec_is_elastic,ier)
  call adios_get_scalar(handle, "ispec_is_poroelastic/local_dim",&
                        local_dim_ispec_is_poroelastic,ier)

  start(1) = local_dim_ispec_is_acoustic * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ispec_is_acoustic/array", 0, 1, &
                           ispec_is_acoustic, ier)
  call adios_schedule_read(handle, sel, "ispec_is_elastic/array", 0, 1, &
                           ispec_is_elastic, ier)
  call adios_schedule_read(handle, sel, "ispec_is_poroelastic/array", 0, 1, &
                           ispec_is_poroelastic, ier)
  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)
  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  ! all processes will have acoustic_simulation set if any flag is .true.
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))
  ! elastic simulation
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  ! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )

  ! checks simulation types are valid
  if( (.not. ACOUSTIC_SIMULATION ) .and. &
      (.not. ELASTIC_SIMULATION ) .and. &
      (.not. POROELASTIC_SIMULATION ) ) then
     call exit_mpi(myrank,'error no simulation type defined')
  endif

  ! outputs total element numbers
  call sum_all_i(count(ispec_is_acoustic(:)),inum)
  if( myrank == 0 ) then
    write(IMAIN,*) 'total acoustic elements    :',inum
  endif
  call sum_all_i(count(ispec_is_elastic(:)),inum)
  if( myrank == 0 ) then
    write(IMAIN,*) 'total elastic elements     :',inum
  endif
  call sum_all_i(count(ispec_is_poroelastic(:)),inum)
  if( myrank == 0 ) then
    write(IMAIN,*) 'total poroelastic elements :',inum
    call flush_IMAIN()
  endif

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  !------------------------------------------------------------------'
  sel_num = 1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)

  NSPEC_CPML = 0
  if( PML_CONDITIONS ) then
    call adios_schedule_read(handle, sel, "/nspec_cpml", 0, 1, nspec_cpml, ier)
    call adios_schedule_read(handle, sel, "/CPML_width_x", 0, 1, &
                             CPML_width_x, ier)
    call adios_schedule_read(handle, sel, "/CPML_width_y", 0, 1, &
                             CPML_width_y, ier)
    call adios_schedule_read(handle, sel, "/CPML_width_x", 0, 1, &
                             CPML_width_z, ier)
    if( nspec_cpml > 0 ) then
      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) &
          .or. SIMULATION_TYPE == 3) then
        call adios_schedule_read(handle, sel, "/nglob_interface_PML_acoustic", &
                                 0, 1, nglob_interface_PML_acoustic, ier)
        call adios_schedule_read(handle, sel, "/nglob_interface_PML_elastic", &
                                 0, 1, nglob_interface_PML_elastic, ier)
      endif
    endif
  endif

  call adios_schedule_read(handle, sel, "/num_abs_boundary_faces", 0, 1, &
                           num_abs_boundary_faces, ier)

  call adios_schedule_read(handle, sel, "/nspec2d_xmin", 0, 1, &
                           nspec2d_xmin, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_xmax", 0, 1, &
                           nspec2d_xmax, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_ymin", 0, 1, &
                           nspec2d_ymin, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_ymax", 0, 1, &
                           nspec2d_ymax, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_bottom", 0, 1, &
                           nspec2d_bottom, ier)
  call adios_schedule_read(handle, sel, "/nspec2d_top", 0, 1, &
                           nspec2d_top, ier)

  call adios_schedule_read(handle, sel, "/num_free_surface_faces", 0, 1, &
                           num_free_surface_faces, ier)
  call adios_schedule_read(handle, sel, "/num_coupling_ac_el_faces", 0, 1, &
                           num_coupling_ac_el_faces, ier)
  call adios_schedule_read(handle, sel, "/num_coupling_ac_po_faces", 0, 1, &
                           num_coupling_ac_po_faces, ier)
  call adios_schedule_read(handle, sel, "/num_coupling_el_po_faces", 0, 1, &
                           num_coupling_el_po_faces, ier)
  call adios_schedule_read(handle, sel, "/num_interfaces_ext_mesh", 0, 1, &
                           num_interfaces_ext_mesh, ier)
  call adios_schedule_read(handle, sel, "/max_nibool_interfaces_ext_mesh", 0, 1, &
                           max_nibool_interfaces_ext_mesh, ier)

  if( ACOUSTIC_SIMULATION ) then
    call adios_schedule_read(handle, sel, "/nspec_inner_acoustic", 0, 1, &
                             nspec_inner_acoustic, ier)
    call adios_schedule_read(handle, sel, "/nspec_outer_acoustic", 0, 1, &
                             nspec_outer_acoustic, ier)
    call adios_schedule_read(handle, sel, "/num_phase_ispec_acoustic", 0, 1, &
                             num_phase_ispec_acoustic, ier)
  endif

  if( ELASTIC_SIMULATION ) then
    call adios_schedule_read(handle, sel, "/nspec_inner_elastic", 0, 1, &
                             nspec_inner_elastic, ier)
    call adios_schedule_read(handle, sel, "/nspec_outer_elastic", 0, 1, &
                             nspec_outer_elastic, ier)
    call adios_schedule_read(handle, sel, "/num_phase_ispec_elastic", 0, 1, &
                             num_phase_ispec_elastic, ier)
  endif

  if( POROELASTIC_SIMULATION) then
    call adios_schedule_read(handle, sel, "/nspec_inner_poroelastic", 0, 1, &
                             nspec_inner_poroelastic, ier)
    call adios_schedule_read(handle, sel, "/nspec_outer_poroelastic", 0, 1, &
                             nspec_outer_poroelastic, ier)
    call adios_schedule_read(handle, sel, "/num_phase_ispec_poroelastic", 0, 1,&
                             num_phase_ispec_poroelastic, ier)
  endif

  num_colors_outer_acoustic = 0
  num_colors_inner_acoustic = 0
  num_colors_outer_elastic = 0
  num_colors_inner_elastic = 0
  if( USE_MESH_COLORING_GPU ) then
    if( ACOUSTIC_SIMULATION ) then
      call adios_schedule_read(handle, sel, "/num_colors_outer_acoustic", &
                               0, 1, num_colors_outer_acoustic, ier)
      call adios_schedule_read(handle, sel, "/num_colors_outer_acoustic", &
                               0, 1, num_colors_inner_acoustic, ier)
    endif
    if( ELASTIC_SIMULATION ) then
      call adios_schedule_read(handle, sel, "/num_colors_outer_elastic", &
                               0, 1, num_colors_outer_elastic, ier)
      call adios_schedule_read(handle, sel, "/num_colors_outer_elastic", &
                               0, 1, num_colors_inner_elastic, ier)
    endif
  endif
  ! Perform the read, so we can use the values.
  call adios_perform_reads(handle, ier)


  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  call adios_get_scalar(handle, "ibool/local_dim",&
                        local_dim_ibool,ier)
  call adios_get_scalar(handle, "x_global/local_dim",&
                        local_dim_x_global,ier)
  call adios_get_scalar(handle, "y_global/local_dim",&
                        local_dim_y_global,ier)
  call adios_get_scalar(handle, "z_global/local_dim",&
                        local_dim_z_global,ier)
  call adios_get_scalar(handle, "xixstore/local_dim",&
                        local_dim_xixstore,ier)
  call adios_get_scalar(handle, "xiystore/local_dim",&
                        local_dim_xiystore,ier)
  call adios_get_scalar(handle, "xizstore/local_dim",&
                        local_dim_xizstore,ier)
  call adios_get_scalar(handle, "etaxstore/local_dim",&
                        local_dim_etaxstore,ier)
  call adios_get_scalar(handle, "etaystore/local_dim",&
                        local_dim_etaystore,ier)
  call adios_get_scalar(handle, "etazstore/local_dim",&
                        local_dim_etazstore,ier)
  call adios_get_scalar(handle, "gammaxstore/local_dim",&
                        local_dim_gammaxstore,ier)
  call adios_get_scalar(handle, "gammaystore/local_dim",&
                        local_dim_gammaystore,ier)
  call adios_get_scalar(handle, "gammazstore/local_dim",&
                        local_dim_gammazstore,ier)
  call adios_get_scalar(handle, "jacobianstore/local_dim",&
                        local_dim_jacobianstore,ier)
  call adios_get_scalar(handle, "kappastore/local_dim",&
                        local_dim_kappastore,ier)
  call adios_get_scalar(handle, "mustore/local_dim",&
                        local_dim_mustore,ier)
  call adios_get_scalar(handle, "rhostore/local_dim",&
                        local_dim_rhostore,ier)
  if( ACOUSTIC_SIMULATION ) then
    call adios_get_scalar(handle, "rmass_acoustic/local_dim",&
                          local_dim_rmass_acoustic,ier)
  endif
  if( ELASTIC_SIMULATION ) then
    call adios_get_scalar(handle, "rmass/local_dim",&
                          local_dim_rmass,ier)

    if( APPROXIMATE_OCEAN_LOAD) then
      call adios_get_scalar(handle, "rmass_ocean_load/local_dim",&
                            local_dim_rmass_ocean_load,ier)
    endif
    call adios_get_scalar(handle, "rho_vp/local_dim",&
                          local_dim_rho_vp,ier)
    call adios_get_scalar(handle, "rho_vs/local_dim",&
                          local_dim_rho_vs,ier)
  endif
  if( POROELASTIC_SIMULATION ) then
    call adios_get_scalar(handle, "rmass_solid_poroelastic/local_dim",&
                          local_dim_rmass_solid_poroelastic,ier)
    call adios_get_scalar(handle, "rmass_fluid_poroelastic/local_dim",&
                          local_dim_rmass_fluid_poroelastic,ier)
    call adios_get_scalar(handle, "rhoarraystore/local_dim",&
                          local_dim_rhoarraystore,ier)
    call adios_get_scalar(handle, "kappaarraystore/local_dim",&
                          local_dim_kappaarraystore,ier)
    call adios_get_scalar(handle, "permstore/local_dim",&
                          local_dim_permstore,ier)
    call adios_get_scalar(handle, "etastore/local_dim",&
                          local_dim_etastore,ier)
    call adios_get_scalar(handle, "tortstore/local_dim",&
                          local_dim_tortstore,ier)
    call adios_get_scalar(handle, "phistore/local_dim",&
                          local_dim_phistore,ier)
    call adios_get_scalar(handle, "rho_vpI/local_dim",&
                          local_dim_rho_vpI,ier)
    call adios_get_scalar(handle, "rho_vpII/local_dim",&
                          local_dim_rho_vpII,ier)
    call adios_get_scalar(handle, "rho_vsI/local_dim",&
                          local_dim_rho_vsI,ier)
  endif
  if( PML_CONDITIONS ) then
    if( nspec_cpml > 0 ) then
      call adios_get_scalar(handle, "CPML_regions/local_dim",&
                            local_dim_CPML_regions, ier)
      call adios_get_scalar(handle, "CPML_to_spec/local_dim",&
                            local_dim_CPML_to_spec, ier)
      call adios_get_scalar(handle, "is_CPML/local_dim",&
                            local_dim_is_CPML, ier)
      call adios_get_scalar(handle, "d_store_x/local_dim",&
                            local_dim_d_store_x, ier)
      call adios_get_scalar(handle, "d_store_y/local_dim",&
                            local_dim_d_store_y, ier)
      call adios_get_scalar(handle, "d_store_z/local_dim",&
                            local_dim_d_store_z, ier)
      call adios_get_scalar(handle, "k_store_x/local_dim",&
                            local_dim_k_store_x, ier)
      call adios_get_scalar(handle, "k_store_y/local_dim",&
                            local_dim_k_store_y, ier)
      call adios_get_scalar(handle, "k_store_z/local_dim",&
                            local_dim_k_store_z, ier)
      call adios_get_scalar(handle, "alpha_store/local_dim",&
                            local_dim_alpha_store, ier)
      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) &
          .or. SIMULATION_TYPE == 3) then
        if(nglob_interface_PML_acoustic > 0) then
          call adios_get_scalar(handle, "points_interface_PML_acoustic/local_dim",&
                                local_dim_points_interface_PML_acoustic, ier)
        endif
        if(nglob_interface_PML_elastic > 0) then
          call adios_get_scalar(handle, "points_interface_PML_elastic/local_dim",&
                                local_dim_points_interface_PML_elastic, ier)
        endif
      endif
    endif
  endif

  if(PML_CONDITIONS)then
    if( num_abs_boundary_faces > 0 ) then
      call adios_get_scalar(handle, "abs_boundary_ispec/local_dim",&
                            local_dim_abs_boundary_ispec,ier)
      call adios_get_scalar(handle, "abs_boundary_ijk/local_dim",&
                            local_dim_abs_boundary_ijk,ier)
      call adios_get_scalar(handle, "abs_boundary_jacobian2Dw/local_dim",&
                            local_dim_abs_boundary_jacobian2Dw,ier)
      call adios_get_scalar(handle, "abs_boundary_normal/local_dim",&
                            local_dim_abs_boundary_normal,ier)
    endif
  else
    if( num_abs_boundary_faces > 0 ) then
      call adios_get_scalar(handle, "abs_boundary_ispec/local_dim",&
                            local_dim_abs_boundary_ispec,ier)
      call adios_get_scalar(handle, "abs_boundary_ijk/local_dim",&
                            local_dim_abs_boundary_ijk,ier)
      call adios_get_scalar(handle, "abs_boundary_jacobian2Dw/local_dim",&
                            local_dim_abs_boundary_jacobian2Dw,ier)
      call adios_get_scalar(handle, "abs_boundary_normal/local_dim",&
                            local_dim_abs_boundary_normal,ier)
      if( STACEY_ABSORBING_CONDITIONS ) then
        ! store mass matrix contributions
        if(ELASTIC_SIMULATION ) then
          call adios_get_scalar(handle, "rmassx/local_dim",&
                                local_dim_rmassx, ier)
          call adios_get_scalar(handle, "rmassy/local_dim",&
                                local_dim_rmassy, ier)
          call adios_get_scalar(handle, "rmassz/local_dim",&
                                local_dim_rmassz, ier)
        endif
        if(ACOUSTIC_SIMULATION) then
          call adios_get_scalar(handle, "rmassz_acoustic/local_dim",&
                                local_dim_rmassz_acoustic, ier)
        endif
      endif
    endif
  endif

  call adios_get_scalar(handle, "ibelm_xmin/local_dim",&
                        local_dim_ibelm_xmin,ier)
  call adios_get_scalar(handle, "ibelm_xmax/local_dim",&
                        local_dim_ibelm_xmax,ier)
  call adios_get_scalar(handle, "ibelm_ymin/local_dim",&
                        local_dim_ibelm_ymin,ier)
  call adios_get_scalar(handle, "ibelm_ymax/local_dim",&
                        local_dim_ibelm_ymax,ier)
  call adios_get_scalar(handle, "ibelm_bottom/local_dim",&
                        local_dim_ibelm_bottom,ier)
  call adios_get_scalar(handle, "ibelm_top/local_dim",&
                        local_dim_ibelm_top,ier)

  if( num_free_surface_faces > 0 ) then
    call adios_get_scalar(handle, "free_surface_ispec/local_dim",&
                          local_dim_free_surface_ispec,ier)
    call adios_get_scalar(handle, "free_surface_ijk/local_dim",&
                          local_dim_free_surface_ijk,ier)
    call adios_get_scalar(handle, "free_surface_jacobian2Dw/local_dim",&
                          local_dim_free_surface_jacobian2Dw,ier)
    call adios_get_scalar(handle, "free_surface_normal/local_dim",&
                          local_dim_free_surface_normal,ier)
  endif
  if( num_coupling_ac_el_faces > 0 ) then
    call adios_get_scalar(handle, "coupling_ac_el_ispec/local_dim",&
                          local_dim_coupling_ac_el_ispec,ier)
    call adios_get_scalar(handle, "coupling_ac_el_ijk/local_dim",&
                          local_dim_coupling_ac_el_ijk,ier)
    call adios_get_scalar(handle, "coupling_ac_el_jacobian2Dw/local_dim",&
                          local_dim_coupling_ac_el_jacobian2Dw,ier)
    call adios_get_scalar(handle, "coupling_ac_el_normal/local_dim",&
                          local_dim_coupling_ac_el_normal,ier)
  endif
  if( num_coupling_ac_po_faces > 0 ) then
    call adios_get_scalar(handle, "coupling_ac_po_ispec/local_dim",&
                          local_dim_coupling_ac_po_ispec, ier)
    call adios_get_scalar(handle, "coupling_ac_po_ijk/local_dim",&
                          local_dim_coupling_ac_po_ijk, ier)
    call adios_get_scalar(handle, "coupling_ac_po_jacobian2Dw/local_dim",&
                          local_dim_coupling_ac_po_jacobian2Dw, ier)
    call adios_get_scalar(handle, "coupling_ac_po_normal/local_dim",&
                          local_dim_coupling_ac_po_normal, ier)
  endif
  if( num_coupling_el_po_faces > 0 ) then
    call adios_get_scalar(handle, "coupling_el_po_ispec/local_dim",&
                          local_dim_coupling_el_po_ispec, ier)
    call adios_get_scalar(handle, "coupling_po_el_ispec/local_dim",&
                          local_dim_coupling_po_el_ispec, ier)
    call adios_get_scalar(handle, "coupling_el_po_ijk/local_dim",&
                          local_dim_coupling_el_po_ijk, ier)
    call adios_get_scalar(handle, "coupling_po_el_ijk/local_dim",&
                          local_dim_coupling_po_el_ijk, ier)
    call adios_get_scalar(handle, "coupling_el_po_jacobian2Dw/local_dim",&
                          local_dim_coupling_el_po_jacobian2Dw, ier)
    call adios_get_scalar(handle, "coupling_el_po_normal/local_dim",&
                          local_dim_coupling_el_po_normal, ier)
  endif
  if( num_interfaces_ext_mesh > 0 ) then
    call adios_get_scalar(handle, "my_neighbours_ext_mesh/local_dim",&
                          local_dim_my_neighbours_ext_mesh,ier)
    call adios_get_scalar(handle, "nibool_interfaces_ext_mesh/local_dim",&
                          local_dim_nibool_interfaces_ext_mesh,ier)
    call adios_get_scalar(handle, "ibool_interfaces_ext_mesh_dummy/local_dim",&
                          local_dim_ibool_interfaces_ext_mesh, ier)
  endif
  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
    call adios_get_scalar(handle, "c11store/local_dim",&
                          local_dim_c11store, ier)
  endif
  call adios_get_scalar(handle, "ispec_is_inner/local_dim",&
                        local_dim_ispec_is_inner,ier)
  if( ACOUSTIC_SIMULATION ) then
    if(num_phase_ispec_acoustic > 0 ) then
      call adios_get_scalar(handle, "phase_ispec_inner_acoustic/local_dim",&
                            local_dim_phase_ispec_inner_acoustic,ier)
    endif
  endif
  if( ELASTIC_SIMULATION ) then
    if(num_phase_ispec_elastic > 0 ) then
      call adios_get_scalar(handle, "phase_ispec_inner_elastic/local_dim",&
                            local_dim_phase_ispec_inner_elastic,ier)
    endif
  endif
  if( POROELASTIC_SIMULATION ) then
    if(num_phase_ispec_poroelastic > 0 ) then
      call adios_get_scalar(handle, "phase_ispec_inner_poroelastic/local_dim",&
                            local_dim_phase_ispec_inner_poroelastic,ier)
    endif
  endif
  if( USE_MESH_COLORING_GPU ) then
    if( ACOUSTIC_SIMULATION ) then
      call adios_get_scalar(handle, "num_elem_colors_acoustic/local_dim",&
                            local_dim_num_elem_colors_acoustic, ier)
    endif
    if( ELASTIC_SIMULATION ) then
      call adios_get_scalar(handle, "num_elem_colors_elastic/local_dim",&
                            local_dim_num_elem_colors_elastic, ier)
    endif
  endif

!TODO
#if 1
  !---------------------------------------------.
  ! Allocate arrays with previously read values |
  !---------------------------------------------'
  if( ACOUSTIC_SIMULATION ) then
    ! potentials
    allocate(potential_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array potential_acoustic'
    allocate(potential_dot_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array potential_dot_acoustic'
    allocate(potential_dot_dot_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array potential_dot_dot_acoustic'
    if( SIMULATION_TYPE /= 1 ) then
      allocate(potential_acoustic_adj_coupling(NGLOB_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating array potential_acoustic_adj_coupling'
    endif
    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmass_acoustic'

    ! initializes mass matrix contribution
    allocate(rmassz_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmassz_acoustic'
    rmassz_acoustic(:) = 0._CUSTOM_REAL
  endif

  ! this array is needed for acoustic simulations but also for elastic
  ! simulations with CPML, thus we now allocate it and read it in all
  ! cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array rhostore'

!TODO
#endif
  ! elastic simulation
  if( ELASTIC_SIMULATION ) then
!TODO
#if 1
    ! displacement,velocity,acceleration
    allocate(displ(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array displ'
    allocate(veloc(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array veloc'
    allocate(accel(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array accel'
    if( SIMULATION_TYPE /= 1 ) then
      allocate(accel_adj_coupling(NDIM,NGLOB_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating array accel_adj_coupling'
    endif

    ! allocates mass matrix
    allocate(rmass(NGLOB_AB),stat=ier)

    if( ier /= 0 ) stop 'error allocating array rmass'
    ! initializes mass matrix contributions
    allocate(rmassx(NGLOB_AB), &
             rmassy(NGLOB_AB), &
             rmassz(NGLOB_AB), &
             stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmassx,rmassy,rmassz'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL

    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_vs'
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
            c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
    if( ier /= 0 ) stop 'error allocating array c11store etc.'

    ! note: currently, they need to be defined, as they are used in the
    !       routine arguments for compute_forces_viscoelastic_Deville()
    allocate(R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
            R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
            R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
            R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
            R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS),stat=ier)
    if( ier /= 0 ) stop 'error allocating array R_xx etc.'

    ! needed for attenuation and/or kernel computations
    allocate(epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if( ier /= 0 ) stop 'error allocating array epsilondev_xx etc.'

    allocate(R_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa,N_SLS),&
             epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa),stat=ier)
    if( ier /= 0 ) stop 'error allocating array R_trace etc.'

    ! note: needed for argument of deville routine
    allocate(epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array epsilon_trace_over_3'

    ! needed for attenuation
    allocate(one_minus_sum_beta(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB), &
            factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array one_minus_sum_beta etc.'

    allocate(one_minus_sum_beta_kappa(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa), &
             factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa),stat=ier)
    if( ier /= 0 ) stop 'error allocating array one_minus_sum_beta_kappa etc.'

    if( APPROXIMATE_OCEAN_LOAD ) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating array rmass_ocean_load'
    else
      ! dummy allocation
      allocate(rmass_ocean_load(1),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy array rmass_ocean_load'
    endif
! TODO
#endif
  else
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif

  if( POROELASTIC_SIMULATION ) then

    if( GPU_MODE ) &
        call exit_mpi(myrank,'POROELASTICITY not supported by GPU mode yet...')

    ! displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(displs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array displs_poroelastic'
    allocate(velocs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array velocs_poroelastic'
    allocate(accels_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array accels_poroelastic'
    allocate(displw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array displw_poroelastic'
    allocate(velocw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array velocw_poroelastic'
    allocate(accelw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array accelw_poroelastic'

    allocate(rmass_solid_poroelastic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmass_fluid_poroelastic'

    allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             etastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             phistore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array poroelastic properties'

    ! needed for kernel computations
    allocate(epsilonsdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonsdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonsdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonsdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonsdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonwdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonwdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonwdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonwdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonwdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array epsilonsdev_xx etc.'

    allocate(epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array epsilons_trace_over_3 etc.'
  endif

  ! C-PML absorbing boundary conditions
  if( PML_CONDITIONS ) then
    allocate(is_CPML(NSPEC_AB),stat=ier)
    if(ier /= 0) stop 'error allocating array is_CPML'

    ! make sure there are no PMLs by default,
    ! and then below if NSPEC_CPML > 0 we will need the real flags
    ! for this mesh from the disk
    is_CPML(:) = .false.

    if( NSPEC_CPML > 0 ) then
      allocate(CPML_regions(NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array CPML_regions'
      allocate(CPML_to_spec(NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array CPML_to_spec'
      allocate(d_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array d_store_x'
      allocate(d_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array d_store_y'
      allocate(d_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array d_store_z'
      allocate(K_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array K_store_x'
      allocate(K_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array K_store_y'
      allocate(K_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array K_store_z'
      allocate(alpha_store(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array alpha_store'

      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if(nglob_interface_PML_acoustic > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if(ier /= 0) stop 'error allocating array points_interface_PML_acoustic'
        endif
        if(nglob_interface_PML_elastic > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if(ier /= 0) stop 'error allocating array points_interface_PML_elastic'
        endif
      endif
    endif
  else
    ! allocate with a dummy size of zero just to be able to use this array
    ! as argument in subroutine calls
    allocate(is_CPML(0),stat=ier)
  endif

  ! absorbing boundary surface
  allocate(abs_boundary_ispec(num_abs_boundary_faces), &
          abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
          abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
          abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array abs_boundary_ispec etc.'

  if (OLD_TEST_TO_FIX_ONE_DAY) then
    ! VM for new method
    !! DK DK for VM VM: these two arrays are undeclared, thus I comment them out i
    ! for now otherwise the code does not compile
    !! VM VM : I already declared these two array in the specfem_par module
    allocate(Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces))
    allocate(Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces))
    open(unit=IIN_veloc_dsm,file=dsmname(1:len_trim(dsmname))//'vel.bin',status='old', &
         action='read',form='unformatted',iostat=ier)
    open(unit=IIN_tract_dsm,file=dsmname(1:len_trim(dsmname))//'tract.bin',status='old', &
         action='read',form='unformatted',iostat=ier)
  else
     allocate(Veloc_dsm_boundary(1,1,1,1))
     allocate(Tract_dsm_boundary(1,1,1,1))
  endif

  allocate(ibelm_xmin(nspec2D_xmin),ibelm_xmax(nspec2D_xmax), &
       ibelm_ymin(nspec2D_ymin),ibelm_ymax(nspec2D_ymax), &
       ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error allocating arrays ibelm_xmin,ibelm_xmax etc.'

  ! free surface
  allocate(free_surface_ispec(num_free_surface_faces), &
          free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
          free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
          free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if(ier /= 0) stop 'error allocating arrays free_surface_ispec etc.'

  ! acoustic-elastic coupling surface
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_normal etc.'

  ! acoustic-poroelastic coupling surface
  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_normal etc.'

  ! elastic-poroelastic coupling surface
  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_ispec(num_coupling_el_po_faces), &
          coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_normal etc.'

  ! MPI interfaces
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh), &
          nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_neighbours_ext_mesh etc.'
  if( num_interfaces_ext_mesh > 0 ) then
    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh, &
                                       num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibool_interfaces_ext_mesh'
  else
    max_nibool_interfaces_ext_mesh = 0
    allocate(ibool_interfaces_ext_mesh(0,0),stat=ier)
  endif

  ! inner / outer elements
  allocate(ispec_is_inner(NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ispec_is_inner'

  if( ACOUSTIC_SIMULATION ) then
    if( num_phase_ispec_acoustic < 0 ) stop 'error acoustic simulation:' // &
                                    'num_phase_ispec_acoustic is < zero'
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_acoustic'
  endif

  if( ELASTIC_SIMULATION ) then
    if( num_phase_ispec_elastic < 0 ) stop 'error elastic simulation:' // &
                                   'num_phase_ispec_elastic is < zero'
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_elastic'
  endif

  if( POROELASTIC_SIMULATION ) then
    if( num_phase_ispec_poroelastic < 0 ) stop 'error poroelastic simulation:'&
                                       'num_phase_ispec_poroelastic is < zero'
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2), &
              stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_poroelastic'
  endif

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! acoustic domain colors
    if( ACOUSTIC_SIMULATION ) then
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic &
                                      + num_colors_inner_acoustic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'
    endif
    ! elastic domain colors
    if( ELASTIC_SIMULATION ) then
      read(27) num_colors_outer_elastic,num_colors_inner_elastic

      allocate(num_elem_colors_elastic(num_colors_outer_elastic &
                                     + num_colors_inner_elastic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'
    endif
  else
    ! allocates dummy arrays
    if( ACOUSTIC_SIMULATION ) then
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic &
                                      + num_colors_inner_acoustic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'
    endif
    if( ELASTIC_SIMULATION ) then
      allocate(num_elem_colors_elastic(num_colors_outer_elastic &
                                     + num_colors_inner_elastic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'
    endif
  endif

  !-----------------------------------.
  ! Read arrays from external_mesh.bp |
  !-----------------------------------'
  start(1) = local_dim_ibool * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibool/array", 0, 1, &
                           ibool, ier)
  start(1) = local_dim_x_global * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "x_global/array", 0, 1, &
                           xstore, ier)
  start(1) = local_dim_y_global * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "y_global/array", 0, 1, &
                           ystore, ier)
  start(1) = local_dim_z_global * myrank
  count_ad(1) = NGLOB_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "z_global/array", 0, 1, &
                           zstore, ier)

  start(1) = local_dim_xixstore * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "xixstore/array", 0, 1, &
                           xix, ier)
  call adios_schedule_read(handle, sel, "xiystore/array", 0, 1, &
                           xiy, ier)
  call adios_schedule_read(handle, sel, "xizstore/array", 0, 1, &
                           xiz, ier)
  call adios_schedule_read(handle, sel, "etaxstore/array", 0, 1, &
                           etax, ier)
  call adios_schedule_read(handle, sel, "etaystore/array", 0, 1, &
                           etay, ier)
  call adios_schedule_read(handle, sel, "etazstore/array", 0, 1, &
                           etaz, ier)
  call adios_schedule_read(handle, sel, "gammaxstore/array", 0, 1, &
                           gammax, ier)
  call adios_schedule_read(handle, sel, "gammaystore/array", 0, 1, &
                           gammay, ier)
  call adios_schedule_read(handle, sel, "gammazstore/array", 0, 1, &
                           gammaz, ier)
  call adios_schedule_read(handle, sel, "jacobianstore/array", 0, 1, &
                           jacobian, ier)
  call adios_schedule_read(handle, sel, "kappastore/array", 0, 1, &
                           kappastore, ier)
  call adios_schedule_read(handle, sel, "mustore/array", 0, 1, &
                           mustore, ier)
  call adios_schedule_read(handle, sel, "rhostore/array", 0, 1, &
                           rhostore, ier)

  if( ACOUSTIC_SIMULATION ) then
    start(1) = local_dim_rmass_acoustic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass_acoustic/array", 0, 1, &
                             rmass_acoustic, ier)
  endif

  if( ELASTIC_SIMULATION ) then
    start(1) = local_dim_rmass * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass/array", 0, 1, &
                             rmass, ier)

    if( APPROXIMATE_OCEAN_LOAD ) then
      ! ocean mass matrix
      start(1) = local_dim_rmass_ocean_load * myrank
      count_ad(1) = NGLOB_AB !nglob_ocean
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "rmass_ocean_load/array", 0, 1, &
                               rmass_ocean_load, ier)
    endif

    !pll material parameters for stacey conditions
    start(1) = local_dim_rho_vp * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rho_vp/array", 0, 1, &
                             rho_vp, ier)
    call adios_schedule_read(handle, sel, "rho_vs/array", 0, 1, &
                             rho_vs, ier)
  endif

  if( POROELASTIC_SIMULATION ) then
    start(1) = local_dim_rmass_solid_poroelastic * myrank
    count_ad(1) = NGLOB_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rmass_solid_poroelastic/array", &
                             0, 1, rmass_solid_poroelastic, ier)
    call adios_schedule_read(handle, sel, "rmass_fluid_poroelastic/array", &
                             0, 1, rmass_fluid_poroelastic, ier)

    start(1) = local_dim_rhoarraystore * myrank
    count_ad(1) = 2 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "rhoarraystore/array", 0, 1, &
                             rhoarraystore, ier)

    start(1) = local_dim_kappaarraystore* myrank
    count_ad(1) = 3 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "kappaarraystore/array", 0, 1, &
                             kappaarraystore, ier)

    start(1) = local_dim_permstore * myrank
    count_ad(1) =  6 * NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "permstore/array", 0, 1, &
                             permstore, ier)

    start(1) = local_dim_etastore * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC_AB
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "etastore/array", 0, 1, &
                             etastore, ier)
    call adios_schedule_read(handle, sel, "tortstore/array", 0, 1, &
                             tortstore, ier)
    call adios_schedule_read(handle, sel, "phistore/array", 0, 1, &
                             phistore, ier)
    call adios_schedule_read(handle, sel, "rho_vpI/array", 0, 1, &
                             rho_vpI, ier)
    call adios_schedule_read(handle, sel, "rho_vpII/array", 0, 1, &
                             rho_vpII, ier)
    call adios_schedule_read(handle, sel, "rho_vsI/array", 0, 1, &
                             rho_vsI, ier)
  endif

  ! C-PML absorbing boundary conditions
  if( PML_CONDITIONS ) then
    if( NSPEC_CPML > 0 ) then

      start(1) = local_dim_CPML_regions * myrank
      count_ad(1) = nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "CPML_regions/array", 0, 1, &
                               CPML_regions, ier)
      call adios_schedule_read(handle, sel, "CPML_to_spec/array", 0, 1, &
                               CPML_to_spec, ier)

      start(1) = local_dim_is_cpml * myrank
      count_ad(1) = NSPEC_AB
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "is_CPML/array", 0, 1, &
                               is_CPML, ier)

      start(1) = local_dim_d_store_x * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_cpml
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "d_store_x/array", 0, 1, &
                               d_store_x, ier)
      call adios_schedule_read(handle, sel, "d_store_y/array", 0, 1, &
                               d_store_y, ier)
      call adios_schedule_read(handle, sel, "d_store_z/array", 0, 1, &
                               d_store_z, ier)
      call adios_schedule_read(handle, sel, "k_store_x/array", 0, 1, &
                               k_store_x, ier)
      call adios_schedule_read(handle, sel, "k_store_y/array", 0, 1, &
                               k_store_y, ier)
      call adios_schedule_read(handle, sel, "k_store_z/array", 0, 1, &
                               k_store_z, ier)
      call adios_schedule_read(handle, sel, "alpha_store/array", 0, 1, &
                               alpha_store, ier)

      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if(nglob_interface_PML_acoustic > 0) then
          start(1) = local_dim_points_interface_PML_acoustic* myrank
          count_ad(1) = nglob_interface_PML_acoustic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call adios_selection_boundingbox (sel , 1, start, count_ad)
          call adios_schedule_read(handle, sel, &
                                   "points_interface_PML_acoustic/array", &
                                   0, 1, points_interface_PML_acoustic , ier)
        endif
        if(nglob_interface_PML_elastic > 0) then
          start(1) = local_dim_points_interface_PML_elastic* myrank
          count_ad(1) = nglob_interface_PML_elastic
          sel_num = sel_num+1
          sel => selections(sel_num)
          call adios_selection_boundingbox (sel , 1, start, count_ad)
          call adios_schedule_read(handle, sel, &
                                   "points_interface_PML_elastic/array", &
                                   0, 1, points_interface_PML_elastic , ier)
        endif
      endif
    endif
  endif

  if(PML_CONDITIONS)then
     if( num_abs_boundary_faces > 0 ) then
        start(1) = local_dim_abs_boundary_ispec * myrank
        count_ad(1) = num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_ispec/array", &
                                 0, 1, abs_boundary_ispec, ier)

        start(1) = local_dim_abs_boundary_ijk * myrank
        count_ad(1) = 3 * NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_ijk/array", &
                                 0, 1, abs_boundary_ijk, ier)

        start(1) = local_dim_abs_boundary_jacobian2Dw * myrank
        count_ad(1) = NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_jacobian2Dw/array", &
                                 0, 1, abs_boundary_jacobian2Dw, ier)

        start(1) = local_dim_abs_boundary_normal * myrank
        count_ad(1) = NDIM * NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_normal/array", &
                                 0, 1, abs_boundary_normal, ier)
     endif
  else
     if( num_abs_boundary_faces > 0 ) then
        start(1) = local_dim_abs_boundary_ispec * myrank
        count_ad(1) = num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_ispec/array", &
                                 0, 1, abs_boundary_ispec, ier)

        start(1) = local_dim_abs_boundary_ijk * myrank
        count_ad(1) = 3 * NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_ijk/array", &
                                 0, 1, abs_boundary_ijk, ier)

        start(1) = local_dim_abs_boundary_jacobian2Dw * myrank
        count_ad(1) = NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_jacobian2Dw/array", &
                                 0, 1, abs_boundary_jacobian2Dw, ier)

        start(1) = local_dim_abs_boundary_normal * myrank
        count_ad(1) = NDIM * NGLLSQUARE * num_abs_boundary_faces
        sel_num = sel_num+1
        sel => selections(sel_num)
        call adios_selection_boundingbox (sel , 1, start, count_ad)
        call adios_schedule_read(handle, sel, "abs_boundary_normal/array", &
                                 0, 1, abs_boundary_normal, ier)

       if( STACEY_ABSORBING_CONDITIONS ) then
          ! store mass matrix contributions
          if(ELASTIC_SIMULATION) then
            start(1) = local_dim_rmassx * myrank
            count_ad(1) = NGLOB_AB  ! == nglob_xy in generate_databse
            sel_num = sel_num+1
            sel => selections(sel_num)
            call adios_selection_boundingbox (sel , 1, start, count_ad)
            call adios_schedule_read(handle, sel, "rmassx/array", 0, 1, &
                                     rmassx, ier)
            call adios_schedule_read(handle, sel, "rmassy/array", 0, 1, &
                                     rmassy, ier)
            call adios_schedule_read(handle, sel, "rmassz/array", 0, 1, &
                                     rmassz, ier)
          endif
          if(ACOUSTIC_SIMULATION) then
            start(1) = local_dim_rmassz_acoustic * myrank
            count_ad(1) = NGLOB_AB ! == nglob_xy in generate_databse
            sel_num = sel_num+1
            sel => selections(sel_num)
            call adios_selection_boundingbox (sel , 1, start, count_ad)
            call adios_schedule_read(handle, sel, "rmassx/array", 0, 1, &
                                     rmassx, ier)
            call adios_schedule_read(handle, sel, "rmassy/array", 0, 1, &
                                     rmassy, ier)
            call adios_schedule_read(handle, sel, "rmassz_acoustic/array", &
                                     0, 1, rmassz_acoustic, ier)
          endif
       endif
     endif
  endif

  start(1) = local_dim_ibelm_xmin * myrank
  count_ad(1) = nspec2D_xmin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_xmin/array", 0, 1, &
                           ibelm_xmin, ier)

  start(1) = local_dim_ibelm_xmax * myrank
  count_ad(1) = nspec2D_xmax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_xmax/array", 0, 1, &
                           ibelm_xmax, ier)

  start(1) = local_dim_ibelm_ymin * myrank
  count_ad(1) = nspec2D_ymin
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_ymin/array", 0, 1, &
                           ibelm_ymin, ier)

  start(1) = local_dim_ibelm_ymax * myrank
  count_ad(1) = nspec2D_ymax
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_ymax/array", 0, 1, &
                           ibelm_ymax, ier)

  start(1) = local_dim_ibelm_bottom * myrank
  count_ad(1) = nspec2D_bottom
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_bottom/array", 0, 1, &
                           ibelm_bottom, ier)

  start(1) = local_dim_ibelm_top * myrank
  count_ad(1) = nspec2D_top
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_top/array", 0, 1, &
                           ibelm_top, ier)

  ! free surface
  if( num_free_surface_faces > 0 ) then

    start(1) = local_dim_free_surface_ispec * myrank
    count_ad(1) = num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ispec/array", 0, 1, &
                             free_surface_ispec, ier)

    start(1) = local_dim_free_surface_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ijk/array", 0, 1, &
                             free_surface_ijk, ier)

    start(1) = local_dim_free_surface_ijk* myrank
    count_ad(1) = NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_ijk/array", 0, 1, &
                             free_surface_ijk, ier)

    start(1) = local_dim_free_surface_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_free_surface_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "free_surface_normal/array", 0, 1, &
                             free_surface_normal, ier)
  endif

  ! acoustic-elastic coupling surface
  if( num_coupling_ac_el_faces > 0 ) then

    start(1) = local_dim_coupling_ac_el_ispec * myrank
    count_ad(1) = num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_ispec/array", 0, 1, &
                             coupling_ac_el_ispec, ier)

    start(1) = local_dim_coupling_ac_el_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_ijk/array", 0, 1, &
                             coupling_ac_el_ijk, ier)

    start(1) = local_dim_coupling_ac_el_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_jacobian2Dw/array", 0, 1, &
                             coupling_ac_el_jacobian2Dw, ier)
    start(1) = local_dim_coupling_ac_el_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_el_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_el_normal/array", 0, 1, &
                             coupling_ac_el_normal, ier)
  endif

  ! acoustic-poroelastic coupling surface
  if( num_coupling_ac_po_faces > 0 ) then

    start(1) = local_dim_coupling_ac_po_ispec * myrank
    count_ad(1) = num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_ispec/array", 0, 1, &
                             coupling_ac_po_ispec, ier)

    start(1) = local_dim_coupling_ac_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_ijk/array", 0, 1, &
                             coupling_ac_po_ijk, ier)

    start(1) = local_dim_coupling_ac_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_jacobian2Dw/array", 0, 1, &
                             coupling_ac_po_jacobian2Dw, ier)

    start(1) = local_dim_coupling_ac_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_ac_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_ac_po_normal/array", 0, 1, &
                             coupling_ac_po_normal, ier)
  endif

  ! elastic-poroelastic coupling surface
  if( num_coupling_el_po_faces > 0 ) then

    start(1) = local_dim_coupling_el_po_ispec * myrank
    count_ad(1) = num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_ispec/array", 0, 1, &
                             coupling_el_po_ispec, ier)
    call adios_schedule_read(handle, sel, "coupling_po_el_ispec/array", 0, 1, &
                             coupling_po_el_ispec, ier)

    start(1) = local_dim_coupling_el_po_ijk * myrank
    count_ad(1) = 3 * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_ijk/array", 0, 1, &
                             coupling_el_po_ijk, ier)
    call adios_schedule_read(handle, sel, "coupling_po_el_ijk/array", 0, 1, &
                             coupling_po_el_ijk, ier)

    start(1) = local_dim_coupling_el_po_jacobian2Dw * myrank
    count_ad(1) = NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_jacobian2Dw/array", 0, 1, &
                             coupling_el_po_jacobian2Dw, ier)

    start(1) = local_dim_coupling_el_po_normal * myrank
    count_ad(1) = NDIM * NGLLSQUARE * num_coupling_el_po_faces
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "coupling_el_po_normal/array", 0, 1, &
                             coupling_el_po_normal, ier)
  endif

  ! MPI interfaces
  if( num_interfaces_ext_mesh > 0 ) then
    start(1) = local_dim_my_neighbours_ext_mesh * myrank
    count_ad(1) = num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "my_neighbours_ext_mesh/array", 0, 1, &
                             my_neighbours_ext_mesh, ier)
    call adios_schedule_read(handle, sel, "nibool_interfaces_ext_mesh/array", 0, 1, &
                             nibool_interfaces_ext_mesh, ier)

    start(1) = local_dim_ibool_interfaces_ext_mesh * myrank
    count_ad(1) = max_nibool_interfaces_ext_mesh * num_interfaces_ext_mesh
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "ibool_interfaces_ext_mesh_dummy/array", 0, 1, &
                             ibool_interfaces_ext_mesh, ier)
  endif

  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
    start(1) = local_dim_c11store * myrank
    count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec_aniso
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count_ad)
    call adios_schedule_read(handle, sel, "c11store/array", 0, 1, &
                             c11store, ier)
    call adios_schedule_read(handle, sel, "c12store/array", 0, 1, &
                             c12store, ier)
    call adios_schedule_read(handle, sel, "c13store/array", 0, 1, &
                             c13store, ier)
    call adios_schedule_read(handle, sel, "c14store/array", 0, 1, &
                             c14store, ier)
    call adios_schedule_read(handle, sel, "c15store/array", 0, 1, &
                             c15store, ier)
    call adios_schedule_read(handle, sel, "c16store/array", 0, 1, &
                             c16store, ier)
    call adios_schedule_read(handle, sel, "c22store/array", 0, 1, &
                             c22store, ier)
    call adios_schedule_read(handle, sel, "c23store/array", 0, 1, &
                             c23store, ier)
    call adios_schedule_read(handle, sel, "c24store/array", 0, 1, &
                             c24store, ier)
    call adios_schedule_read(handle, sel, "c25store/array", 0, 1, &
                             c25store, ier)
    call adios_schedule_read(handle, sel, "c26store/array", 0, 1, &
                             c26store, ier)
    call adios_schedule_read(handle, sel, "c33store/array", 0, 1, &
                             c33store, ier)
    call adios_schedule_read(handle, sel, "c34store/array", 0, 1, &
                             c34store, ier)
    call adios_schedule_read(handle, sel, "c35store/array", 0, 1, &
                             c35store, ier)
    call adios_schedule_read(handle, sel, "c36store/array", 0, 1, &
                             c36store, ier)
    call adios_schedule_read(handle, sel, "c44store/array", 0, 1, &
                             c44store, ier)
    call adios_schedule_read(handle, sel, "c45store/array", 0, 1, &
                             c45store, ier)
    call adios_schedule_read(handle, sel, "c46store/array", 0, 1, &
                             c46store, ier)
    call adios_schedule_read(handle, sel, "c55store/array", 0, 1, &
                             c55store, ier)
    call adios_schedule_read(handle, sel, "c56store/array", 0, 1, &
                             c56store, ier)
    call adios_schedule_read(handle, sel, "c66store/array", 0, 1, &
                             c66store, ier)
  endif

  ! inner / outer elements
  start(1) = local_dim_ispec_is_inner * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ispec_is_inner/array", 0, 1, &
                           ispec_is_inner, ier)

  if( ACOUSTIC_SIMULATION ) then
    if(num_phase_ispec_acoustic > 0 ) then
      start(1) = local_dim_phase_ispec_inner_acoustic * myrank
      count_ad(1) = num_phase_ispec_acoustic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_acoustic/array", 0, 1, &
                               phase_ispec_inner_acoustic, ier)
    endif
  endif

  if( ELASTIC_SIMULATION ) then
    if(num_phase_ispec_elastic > 0 ) then
      start(1) = local_dim_phase_ispec_inner_elastic * myrank
      count_ad(1) = num_phase_ispec_elastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_elastic/array", 0, 1, &
                               phase_ispec_inner_elastic, ier)
    endif
  endif

  if( POROELASTIC_SIMULATION ) then
    if(num_phase_ispec_poroelastic > 0 ) then
      start(1) = local_dim_phase_ispec_inner_poroelastic * myrank
      count_ad(1) = num_phase_ispec_poroelastic * 2
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "phase_ispec_inner_poroelastic/array", 0, 1, &
                               phase_ispec_inner_poroelastic, ier)
    endif
  endif

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! acoustic domain colors
    if( ACOUSTIC_SIMULATION ) then
      start(1) = local_dim_num_elem_colors_acoustic * myrank
      count_ad(1) = num_colors_outer_acoustic + num_colors_inner_acoustic
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "num_elem_colors_acoustic/array", 0, 1, &
                               num_elem_colors_acoustic, ier)
    endif
    ! elastic domain colors
    if( ELASTIC_SIMULATION ) then
      start(1) = local_dim_num_elem_colors_elastic * myrank
      count_ad(1) = num_colors_outer_elastic + num_colors_inner_elastic
      sel_num = sel_num+1
      sel => selections(sel_num)
      call adios_selection_boundingbox (sel , 1, start, count_ad)
      call adios_schedule_read(handle, sel, "num_elem_colors_elastic/array", 0, 1, &
                               num_elem_colors_elastic, ier)
    endif
  endif

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call adios_perform_reads(handle, ier)
  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)

  !call read_mesh_databases2()
  !call check_mesh_database()

  ! debug
  !call sum_all_i(num_interfaces_ext_mesh,inum)
  !if(myrank == 0) then
  !  write(IMAIN,*) 'number of MPI partition interfaces: ',inum
  !  write(IMAIN,*)
  !endif

  ! MPI communications
  allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    request_send_vector_ext_mesh(num_interfaces_ext_mesh), &
    request_recv_vector_ext_mesh(num_interfaces_ext_mesh), &
    request_send_scalar_ext_mesh(num_interfaces_ext_mesh), &
    request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), &
    buffer_send_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_recv_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_send_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    buffer_recv_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
    request_send_vector_ext_mesh_s(num_interfaces_ext_mesh), &
    request_recv_vector_ext_mesh_s(num_interfaces_ext_mesh), &
    request_send_vector_ext_mesh_w(num_interfaces_ext_mesh), &
    request_recv_vector_ext_mesh_w(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array buffer_send_vector_ext_mesh etc.'

  ! gets model dimensions
  minl = minval( xstore )
  maxl = maxval( xstore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LONGITUDE_MIN = min_all
  LONGITUDE_MAX = max_all

  minl = minval( ystore )
  maxl = maxval( ystore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LATITUDE_MIN = min_all
  LATITUDE_MAX = max_all

  ! checks courant criteria on mesh
  if( ELASTIC_SIMULATION ) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                              ibool,xstore,ystore,zstore, &
                              kappastore,mustore,rho_vp,rho_vs, &
                              DT,model_speed_max,min_resolved_period, &
                              LOCAL_PATH,SAVE_MESH_FILES)

  else if( POROELASTIC_SIMULATION ) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    DT,model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  else if( ACOUSTIC_SIMULATION ) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_vs'
    rho_vp = sqrt( kappastore / rhostore ) * rhostore
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                              ibool,xstore,ystore,zstore, &
                              kappastore,mustore,rho_vp,rho_vs, &
                              DT,model_speed_max,min_resolved_period, &
                              LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  endif

  ! reads adjoint parameters
  call read_mesh_databases_adjoint()

end subroutine read_mesh_databases_adios


!-------------------------------------------------------------------------------
!> Reads in moho meshes
subroutine read_moho_mesh_adjoint_adios()
  use mpi
  use adios_read_mod

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none

  character(len=256) :: database_name
  integer(kind=8) :: handle

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num
  integer(kind=8), pointer :: sel => null()
  integer(kind=8), dimension(1) :: start, count_ad

  integer :: local_dim_ibelm_moho_bot,  local_dim_ibelm_moho_top,  &
             local_dim_ijk_moho_bot,    local_dim_ijk_moho_top,    &
             local_dim_normal_moho_bot, local_dim_normal_moho_top, &
             local_dim_is_moho_bot,     local_dim_is_moho_top

  integer :: ier

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  sel_num = 0

  database_name = adjustl(LOCAL_PATH)
  database_name = database_name(1:len_trim(database_name)) // "/moho.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                               "verbose=1", ier)
  call adios_read_open_file (handle, database_name, 0, MPI_COMM_WORLD, ier)

  !------------------------------------------------------------------.
  ! Get scalar values. Might be differents for different processors. |
  ! Hence the selection writeblock.                                  |
  ! ONLY NSPEC_AB and NGLOB_AB
  !------------------------------------------------------------------'
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(handle, sel, "/nspec2d_moho", 0, 1, &
                           NSPEC2D_MOHO, ier)
  call adios_perform_reads(handle, ier)

  !----------------------------------------------.
  ! Fetch values to compute the simulation type. |
  !----------------------------------------------'
  sel_num = 0
  call adios_get_scalar(handle, "ibelm_moho_bot/local_dim",&
                        local_dim_ibelm_moho_bot ,ier)
  call adios_get_scalar(handle, "ibelm_moho_top/local_dim",&
                        local_dim_ibelm_moho_top ,ier)

  call adios_get_scalar(handle, "ijk_moho_bot/local_dim",&
                        local_dim_ijk_moho_bot ,ier)
  call adios_get_scalar(handle, "ijk_moho_top/local_dim",&
                        local_dim_ijk_moho_top ,ier)

  call adios_get_scalar(handle,"normal_moho_bot /local_dim",&
                        local_dim_normal_moho_bot ,ier)
  call adios_get_scalar(handle, "normal_moho_top/local_dim",&
                        local_dim_normal_moho_top ,ier)

  call adios_get_scalar(handle, "is_moho_bot/local_dim",&
                        local_dim_is_moho_bot ,ier)
  call adios_get_scalar(handle, "is_moho_top/local_dim",&
                        local_dim_is_moho_top ,ier)

  !---------------------------------------------.
  ! Allocate arrays with previously read values |
  !---------------------------------------------'
  allocate(ibelm_moho_bot(NSPEC2D_MOHO), &
          ibelm_moho_top(NSPEC2D_MOHO), &
          normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
          normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
          ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO), &
          ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibelm_moho_bot etc.'

  !-----------------------------------.
  ! Read arrays from external_mesh.bp |
  !-----------------------------------'
  start(1) = local_dim_ibelm_moho_bot * myrank
  count_ad(1) = NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_moho_bot/array", 0, 1, &
                           ibelm_moho_bot, ier)
  start(1) = local_dim_ibelm_moho_top * myrank
  count_ad(1) = NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ibelm_moho_top/array", 0, 1, &
                           ibelm_moho_top, ier)

  start(1) = local_dim_ijk_moho_bot * myrank
  count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ijk_moho_bot/array", 0, 1, &
                           ijk_moho_bot, ier)
  start(1) = local_dim_ijk_moho_top * myrank
  count_ad(1) = 3 * NGLLSQUARE * NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "ijk_moho_top/array", 0, 1, &
                           ijk_moho_top, ier)

  start(1) = local_dim_normal_moho_bot * myrank
  count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "normal_moho_bot/array", 0, 1, &
                           normal_moho_bot, ier)
  start(1) = local_dim_normal_moho_top * myrank
  count_ad(1) = NDIM * NGLLSQUARE * NSPEC2D_MOHO
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "normal_moho_top/array", 0, 1, &
                           normal_moho_top, ier)

  start(1) = local_dim_is_moho_bot * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "is_moho_bot/array", 0, 1, &
                           is_moho_bot, ier)
  start(1) = local_dim_is_moho_top * myrank
  count_ad(1) = NSPEC_AB
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(handle, sel, "is_moho_top/array", 0, 1, &
                           is_moho_top, ier)

  !---------------------------------------------------------------.
  ! Perform the reads and close the ADIOS 'external_mesh.bp' file |
  !---------------------------------------------------------------'
  call adios_perform_reads(handle, ier)
  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)

end subroutine read_moho_mesh_adjoint_adios
