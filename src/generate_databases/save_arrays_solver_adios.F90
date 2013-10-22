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

!==============================================================================
!> \file save_arrays_solver_adios.F90
!!
!! \author MPBL
!==============================================================================

!==============================================================================
!> \def STRINGIFY_VAR(a)
!! Macro taking a variable and returning the stringified variable and
!! the variable itself.
!! STRINGIFY_VAR(x) expand as:
!!   "x", x
!! x being the variable name inside the code.
#ifdef __INTEL_COMPILER
#define STRINGIFY_VAR(a) #a, a
#else
#define STRINGIFY_VAR(a) "a", a
#endif

! for external mesh

subroutine save_arrays_solver_ext_mesh_adios(nspec, nglob,                   &
                                        APPROXIMATE_OCEAN_LOAD,              &
                                       ibool, num_interfaces_ext_mesh,       &
                                       my_neighbours_ext_mesh,               &
                                       nibool_interfaces_ext_mesh,           &
                                       max_interface_size_ext_mesh,          &
                                       ibool_interfaces_ext_mesh,            &
                                       SAVE_MESH_FILES,ANISOTROPY)

  use generate_databases_par, only: nspec_cpml, CPML_width_x, CPML_width_y, &
                                    CPML_width_z, CPML_to_spec,             &
                                    CPML_regions, is_CPML, nspec_cpml_tot,  &
                                    d_store_x, d_store_y, d_store_z,        &
                                    k_store_x, k_store_y, k_store_z,        &
                                    alpha_store,                            &
                                    nspec2D_xmin, nspec2D_xmax,             &
                                    nspec2D_ymin, nspec2D_ymax,             &
                                    NSPEC2D_BOTTOM, NSPEC2D_TOP,            &
                                    ibelm_xmin, ibelm_xmax,ibelm_ymin,      &
                                    ibelm_ymax, ibelm_bottom, ibelm_top,    &
                                    PML_CONDITIONS,                         &
                                    !for adjoint tomography
                                    SIMULATION_TYPE, SAVE_FORWARD,          &
                                    mask_ibool_interior_domain,             &
                                    nglob_interface_PML_acoustic,           &
                                    points_interface_PML_acoustic,          &
                                    nglob_interface_PML_elastic,            &
                                    points_interface_PML_elastic,           &
                                    STACEY_ABSORBING_CONDITIONS,            &
                                    LOCAL_PATH, myrank, sizeprocs,          &
                                    nspec_ab
  use mpi
  use adios_helpers_mod
  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec,nglob
  ! ocean load
  logical :: APPROXIMATE_OCEAN_LOAD
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  ! MPI interfaces
  integer :: num_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer :: max_interface_size_ext_mesh
  integer, dimension(NGLLX * NGLLX * max_interface_size_ext_mesh, &
                     num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  logical :: SAVE_MESH_FILES
  logical :: ANISOTROPY

  ! local parameters
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  integer :: ier,i

  !--- Local parameters for ADIOS ---
  character(len=256) :: output_name
  character(len=64), parameter :: group_name  = "SPECFEM3D_EXTERNAL_MESH"
  integer(kind=8) :: group, handle
  integer(kind=8) :: groupsize, totalsize
  integer :: local_dim

  !--- Variables to allreduce - wmax stands for world_max
  integer :: nglob_wmax, nspec_wmax, nspec_cpml_wmax,              &
             CPML_width_x_wmax, CPML_width_y_wmax, CPML_width_z_wmax,          &
             nglob_interface_PML_acoustic_wmax,                                &
             nglob_interface_PML_elastic_wmax, num_abs_boundary_faces_wmax,    &
             nspec2d_xmin_wmax, nspec2d_xmax_wmax,                             &
             nspec2d_ymin_wmax, nspec2d_ymax_wmax,                             &
             nspec2d_bottom_wmax, nspec2d_top_wmax,                            &
             num_free_surface_faces_wmax, num_coupling_ac_el_faces_wmax,       &
             num_coupling_ac_po_faces_wmax, num_coupling_el_po_faces_wmax,     &
             num_interfaces_ext_mesh_wmax, max_interface_size_ext_mesh_wmax,   &
             nspec_inner_acoustic_wmax, nspec_outer_acoustic_wmax,             &
             num_phase_ispec_acoustic_wmax, nspec_inner_elastic_wmax,          &
             nspec_outer_elastic_wmax, num_phase_ispec_elastic_wmax,           &
             nspec_inner_poroelastic_wmax, nspec_outer_poroelastic_wmax,       &
             num_phase_ispec_poroelastic_wmax, num_colors_outer_acoustic_wmax, &
             num_colors_inner_acoustic_wmax, num_colors_outer_elastic_wmax,    &
             num_colors_inner_elastic_wmax, nglob_dummy_wmax,                  &
             nglob_ocean_wmax, nglob_xy_wmax, nspec_ab_wmax, nspec_aniso_wmax, &
             max_nibool_interfaces_ext_mesh_wmax

  integer, parameter :: num_vars = 40
  integer, dimension(num_vars) :: max_global_values

  !---------------------------.
  ! Setup the values to write |
  !---------------------------'
  !MPI interfaces
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))

  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh, &
                                           num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = &
         ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo

  !-----------------------------------------------------------------.
  ! Get maximum value for each variable used to define a local_dim. |
  ! ADIOS write equally sized chunks for each processor.            |
  !-----------------------------------------------------------------'
  ! Filling a temporary array to avoid doing allreduces for each var.
  max_global_values(1)  = nglob
  max_global_values(2)  = nspec
  max_global_values(3)  = max_nibool_interfaces_ext_mesh
  max_global_values(4)  = nspec_cpml
  max_global_values(5)  = CPML_width_x
  max_global_values(6)  = CPML_width_y
  max_global_values(7)  = CPML_width_z
  max_global_values(8)  = nglob_interface_PML_acoustic
  max_global_values(9)  = nglob_interface_PML_elastic
  max_global_values(10) = num_abs_boundary_faces
  max_global_values(11) = nspec2d_xmin
  max_global_values(12) = nspec2d_xmax
  max_global_values(13) = nspec2d_ymin
  max_global_values(14) = nspec2d_ymax
  max_global_values(15) = nspec2d_bottom
  max_global_values(16) = nspec2d_top
  max_global_values(17) = num_free_surface_faces
  max_global_values(18) = num_coupling_ac_el_faces
  max_global_values(19) = num_coupling_ac_po_faces
  max_global_values(20) = num_coupling_el_po_faces
  max_global_values(21) = num_interfaces_ext_mesh
  max_global_values(22) = max_interface_size_ext_mesh
  max_global_values(23) = nspec_inner_acoustic
  max_global_values(24) = nspec_outer_acoustic
  max_global_values(25) = num_phase_ispec_acoustic
  max_global_values(26) = nspec_inner_elastic
  max_global_values(27) = nspec_outer_elastic
  max_global_values(28) = num_phase_ispec_elastic
  max_global_values(29) = nspec_inner_poroelastic
  max_global_values(30) = nspec_outer_poroelastic
  max_global_values(31) = num_phase_ispec_poroelastic
  max_global_values(32) = num_colors_outer_acoustic
  max_global_values(33) = num_colors_inner_acoustic
  max_global_values(34) = num_colors_outer_elastic
  max_global_values(35) = num_colors_inner_elastic
  max_global_values(36) = nglob_dummy
  max_global_values(37) = nglob_ocean
  max_global_values(38) = nglob_xy
  max_global_values(39) = nspec_ab
  max_global_values(40) = nspec_aniso

  call MPI_Allreduce(MPI_IN_PLACE, max_global_values, num_vars, &
                     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  if( ier /= 0 ) call exit_MPI(myrank,'Allreduce to get max values failed.')

  nglob_wmax                          = max_global_values(1)
  nspec_wmax                          = max_global_values(2)
  max_nibool_interfaces_ext_mesh_wmax = max_global_values(3)
  nspec_cpml_wmax                     = max_global_values(4)
  CPML_width_x_wmax                   = max_global_values(5)
  CPML_width_y_wmax                   = max_global_values(6)
  CPML_width_z_wmax                   = max_global_values(7)
  nglob_interface_PML_acoustic_wmax   = max_global_values(8)
  nglob_interface_PML_elastic_wmax    = max_global_values(9)
  num_abs_boundary_faces_wmax         = max_global_values(10)
  nspec2d_xmin_wmax                   = max_global_values(11)
  nspec2d_xmax_wmax                   = max_global_values(12)
  nspec2d_ymin_wmax                   = max_global_values(13)
  nspec2d_ymax_wmax                   = max_global_values(14)
  nspec2d_bottom_wmax                 = max_global_values(15)
  nspec2d_top_wmax                    = max_global_values(16)
  num_free_surface_faces_wmax         = max_global_values(17)
  num_coupling_ac_el_faces_wmax       = max_global_values(18)
  num_coupling_ac_po_faces_wmax       = max_global_values(19)
  num_coupling_el_po_faces_wmax       = max_global_values(20)
  num_interfaces_ext_mesh_wmax        = max_global_values(21)
  max_interface_size_ext_mesh_wmax    = max_global_values(22)
  nspec_inner_acoustic_wmax           = max_global_values(23)
  nspec_outer_acoustic_wmax           = max_global_values(24)
  num_phase_ispec_acoustic_wmax       = max_global_values(25)
  nspec_inner_elastic_wmax            = max_global_values(26)
  nspec_outer_elastic_wmax            = max_global_values(27)
  num_phase_ispec_elastic_wmax        = max_global_values(28)
  nspec_inner_poroelastic_wmax        = max_global_values(29)
  nspec_outer_poroelastic_wmax        = max_global_values(30)
  num_phase_ispec_poroelastic_wmax    = max_global_values(31)
  num_colors_outer_acoustic_wmax      = max_global_values(32)
  num_colors_inner_acoustic_wmax      = max_global_values(33)
  num_colors_outer_elastic_wmax       = max_global_values(34)
  num_colors_inner_elastic_wmax       = max_global_values(35)
  nglob_dummy_wmax                    = max_global_values(36)
  nglob_ocean_wmax                    = max_global_values(37)
  nglob_xy_wmax                       = max_global_values(38)
  nspec_ab_wmax                       = max_global_values(39)
  nspec_aniso_wmax                    = max_global_values(40)

  ! save arrays for the solver to run.
  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/external_mesh.bp"
  call adios_declare_group(group, group_name, "", 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nglob))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(ibool))

  local_dim = nglob_dummy_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "x_global", xstore_dummy)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "y_global", ystore_dummy)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "z_global", zstore_dummy)

  ! this array is needed for acoustic simulations but also for elastic
  ! simulations with CPML, thus we allocate it and read it in all cases
  ! (whether the simulation is acoustic, elastic, or acoustic/elastic)
  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(xixstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(xiystore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(xizstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(etaxstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(etaystore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(etazstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(gammaxstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(gammaystore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(gammazstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(jacobianstore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(kappastore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(mustore))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "", STRINGIFY_VAR(rhostore))

  local_dim = nspec_wmax
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",               &
                                   STRINGIFY_VAR(ispec_is_acoustic))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",               &
                                   STRINGIFY_VAR(ispec_is_elastic))
  call define_adios_global_array1D(group, groupsize, &
                                   local_dim, "",               &
                                   STRINGIFY_VAR(ispec_is_poroelastic))

  ! acoustic
  if( ACOUSTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rmass_acoustic))
  endif

  ! elastic
  if( ELASTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "", STRINGIFY_VAR(rmass))

    if( APPROXIMATE_OCEAN_LOAD) then
      local_dim = nglob_ocean_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(rmass_ocean_load))
    endif
    !pll Stacey
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "", STRINGIFY_VAR(rho_vp))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "", STRINGIFY_VAR(rho_vs))
  endif

  ! poroelastic
  if( POROELASTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rmass_solid_poroelastic))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rmass_fluid_poroelastic))
    local_dim = 2 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rhoarraystore))
    local_dim = 3 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(kappaarraystore))
    local_dim = 6 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(permstore))
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(etastore))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(tortstore))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(phistore))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rho_vpI))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rho_vpII))
    call define_adios_global_array1D(group, groupsize, &
                                     local_dim, "",               &
                                     STRINGIFY_VAR(rho_vsI))
  endif

  ! C-PML absorbing boundary conditions
  if( PML_CONDITIONS ) then
    call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec_cpml))
    call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(CPML_width_x))
    call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(CPML_width_y))
    call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(CPML_width_z))
    if( nspec_cpml > 0 ) then
      local_dim = nspec_cpml_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(CPML_regions))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(CPML_to_spec))
      local_dim = nspec_ab_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(is_CPML))
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_cpml_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(d_store_x))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(d_store_y))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(d_store_z))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(k_store_x))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(k_store_y))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "", STRINGIFY_VAR(k_store_z))
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(alpha_store))
      ! -----------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior
      ! computational domain
      ! -----------------------------------------------------------------------
      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) &
          .or. SIMULATION_TYPE == 3) then
        call define_adios_scalar(group, groupsize, "", &
                                 STRINGIFY_VAR(nglob_interface_PML_acoustic))
        call define_adios_scalar(group, groupsize, "", &
                                 STRINGIFY_VAR(nglob_interface_PML_elastic))
        if(nglob_interface_PML_acoustic > 0) then
          local_dim = nglob_interface_PML_acoustic_wmax
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                  STRINGIFY_VAR(points_interface_PML_acoustic))
        endif
        if(nglob_interface_PML_elastic > 0) then
          local_dim = nglob_interface_PML_elastic_wmax
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                  STRINGIFY_VAR(points_interface_PML_elastic))
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_abs_boundary_faces))
  if(PML_CONDITIONS)then
    if( num_abs_boundary_faces > 0 ) then
      local_dim = num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_ispec))
      local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_ijk))
      local_dim = NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_jacobian2Dw))
      local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_normal))
    endif
  else
    if( num_abs_boundary_faces > 0 ) then
      local_dim = num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_ispec))
      local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_ijk))
      local_dim = NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_jacobian2Dw))
      local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_wmax
      call define_adios_global_array1D(group, groupsize, &
                                       local_dim, "",               &
                                       STRINGIFY_VAR(abs_boundary_normal))
      if( STACEY_ABSORBING_CONDITIONS ) then
        ! store mass matrix contributions
        if(ELASTIC_SIMULATION ) then
          local_dim = nglob_xy_wmax
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                           STRINGIFY_VAR(rmassx))
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                           STRINGIFY_VAR(rmassy))
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                           STRINGIFY_VAR(rmassz))
        endif
        if(ACOUSTIC_SIMULATION) then
          local_dim = nglob_xy_wmax
          call define_adios_global_array1D(group, groupsize, &
                                           local_dim, "",               &
                                           STRINGIFY_VAR(rmassz_acoustic))
        endif
      endif
    endif
  endif

  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_xmin))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_xmax))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_ymin))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_ymax))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_bottom))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec2d_top))

  if (nspec2d_xmin .ne. 0) then
    local_dim = nspec2d_xmin_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_xmin))
  endif
  if (nspec2d_xmax .ne. 0) then
    local_dim = nspec2d_xmax_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_xmax))
  endif
  if (nspec2d_ymin .ne. 0) then
    local_dim = nspec2d_ymin_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_ymin))
  endif
  if (nspec2d_ymax .ne. 0) then
    local_dim = nspec2d_ymax_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_ymax))
  endif
  if (nspec2d_bottom .ne. 0) then
    local_dim = nspec2d_bottom_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_bottom))
  endif
  if (nspec2d_top .ne. 0) then
    local_dim = nspec2d_top_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(ibelm_top))
  endif

  ! free surface
  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_free_surface_faces))
  if( num_free_surface_faces > 0 ) then
    local_dim = num_free_surface_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(free_surface_ispec))
    local_dim = 3 * NGLLSQUARE * num_free_surface_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(free_surface_ijk))
    local_dim = NGLLSQUARE * num_free_surface_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(free_surface_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_free_surface_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(free_surface_normal))
  endif

  ! acoustic-elastic coupling surface
  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_coupling_ac_el_faces))
  if( num_coupling_ac_el_faces > 0 ) then
    local_dim = num_coupling_ac_el_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_el_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_el_ijk))
    local_dim = NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, "", &
                                     STRINGIFY_VAR(coupling_ac_el_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_el_normal))
  endif

  ! acoustic-poroelastic coupling surface
  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_coupling_ac_po_faces))
  if( num_coupling_ac_po_faces > 0 ) then
    local_dim = num_coupling_ac_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_po_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_po_ijk))
    local_dim = NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                  "", STRINGIFY_VAR(coupling_ac_po_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_ac_po_normal))
  endif

  ! elastic-poroelastic coupling surface
  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_coupling_el_po_faces))
  if( num_coupling_el_po_faces > 0 ) then
    local_dim = num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_el_po_ispec))
    local_dim = num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_po_el_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_el_po_ijk))
    local_dim = 3 * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_po_el_ijk))
    local_dim = NGLLSQUARE * num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                  "", STRINGIFY_VAR(coupling_el_po_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(coupling_el_po_normal))
  endif

  call define_adios_scalar(group, groupsize, "", &
                           STRINGIFY_VAR(num_interfaces_ext_mesh))
  if( num_interfaces_ext_mesh > 0 ) then
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(max_nibool_interfaces_ext_mesh))
    local_dim = num_interfaces_ext_mesh_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(my_neighbours_ext_mesh))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                 "", STRINGIFY_VAR(nibool_interfaces_ext_mesh))
    local_dim = max_nibool_interfaces_ext_mesh_wmax &
              * num_interfaces_ext_mesh_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                             "", STRINGIFY_VAR(ibool_interfaces_ext_mesh_dummy))
  endif

  ! anisotropy
  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_aniso_wmax
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c11store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c12store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c13store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c14store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c15store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c16store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c22store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c23store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c24store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c25store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c26store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c33store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c34store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c35store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c36store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c44store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c45store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c46store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c55store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c56store))
    call define_adios_global_array1D(group, groupsize, local_dim, &
                                     "", STRINGIFY_VAR(c66store))
  endif

  ! inner/outer elements
  local_dim = nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(ispec_is_inner))

  if( ACOUSTIC_SIMULATION ) then
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_inner_acoustic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_outer_acoustic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(num_phase_ispec_acoustic))
    if(num_phase_ispec_acoustic > 0 ) then
      local_dim = num_phase_ispec_acoustic_wmax * 2
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                 "", STRINGIFY_VAR(phase_ispec_inner_acoustic))
    endif
  endif

  if( ELASTIC_SIMULATION ) then
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_inner_elastic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_outer_elastic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(num_phase_ispec_elastic))
    if(num_phase_ispec_elastic > 0 ) then
      local_dim = num_phase_ispec_elastic_wmax * 2
      call define_adios_global_array1D(group, groupsize, local_dim, &
                                 "", STRINGIFY_VAR(phase_ispec_inner_elastic))
    endif
  endif

  if( POROELASTIC_SIMULATION ) then
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_inner_poroelastic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(nspec_outer_poroelastic))
    call define_adios_scalar(group, groupsize, "", &
                             STRINGIFY_VAR(num_phase_ispec_poroelastic))
    if(num_phase_ispec_poroelastic > 0 ) then
      local_dim = num_phase_ispec_poroelastic_wmax * 2
      call define_adios_global_array1D(group, groupsize, local_dim, &
                              "", STRINGIFY_VAR(phase_ispec_inner_poroelastic))
    endif
  endif

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    if( ACOUSTIC_SIMULATION ) then
      call define_adios_scalar(group, groupsize, "", &
                               STRINGIFY_VAR(num_colors_outer_acoustic))
      call define_adios_scalar(group, groupsize, "", &
                               STRINGIFY_VAR(num_colors_inner_acoustic))
      local_dim = num_colors_outer_acoustic_wmax &
                + num_colors_inner_acoustic_wmax
      call define_adios_global_array1D(group, groupsize, local_dim, &
                              "", STRINGIFY_VAR(num_elem_colors_acoustic))
    endif
    if( ELASTIC_SIMULATION ) then
      call define_adios_scalar(group, groupsize, "", &
                               STRINGIFY_VAR(num_colors_outer_elastic))
      call define_adios_scalar(group, groupsize, "", &
                               STRINGIFY_VAR(num_colors_inner_elastic))
      local_dim = num_colors_outer_elastic_wmax &
                + num_colors_inner_elastic_wmax
      call define_adios_global_array1D(group, groupsize, local_dim, &
                              "", STRINGIFY_VAR(num_elem_colors_elastic))
    endif
  endif

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call adios_open(handle, group_name, output_name, "w", &
                  MPI_COMM_WORLD, ier);
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, STRINGIFY_VAR(nspec), ier)
  call adios_write(handle, STRINGIFY_VAR(nglob), ier)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibool))

  local_dim = nglob_dummy_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "x_global", xstore_dummy)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "y_global", ystore_dummy)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "z_global", zstore_dummy)

  ! this array is needed for acoustic simulations but also for elastic
  ! simulations with CPML, thus we allocate it and read it in all cases
  ! (whether the simulation is acoustic, elastic, or acoustic/elastic)
  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(xixstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(xiystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(xizstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(etaxstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(etaystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(etazstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(gammaxstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(gammaystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(gammazstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(jacobianstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(kappastore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(mustore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                   local_dim, STRINGIFY_VAR(rhostore))

  local_dim = nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ispec_is_acoustic))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ispec_is_elastic))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ispec_is_poroelastic))

  ! acoustic
  if( ACOUSTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rmass_acoustic))
  endif

  ! elastic
  if( ELASTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     local_dim, STRINGIFY_VAR(rmass))
    if( APPROXIMATE_OCEAN_LOAD) then
      local_dim = nglob_ocean_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(rmass_ocean_load))
    endif
    !pll Stacey
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rho_vp))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rho_vs))
  endif

! poroelastic
  if( POROELASTIC_SIMULATION ) then
    local_dim = nglob_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rmass_solid_poroelastic))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rmass_fluid_poroelastic))

    local_dim = 2 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rhoarraystore))

    local_dim = 3 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(kappaarraystore))

    local_dim = 6 * NGLLX * NGLLY * NGLLZ * nspec_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(permstore))

    local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(etastore))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(tortstore))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(phistore))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rho_vpI))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rho_vpII))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rho_vsI))
  endif

! C-PML absorbing boundary conditions
  if( PML_CONDITIONS ) then
    call adios_write(handle, STRINGIFY_VAR(nspec_cpml), ier)
    call adios_write(handle, STRINGIFY_VAR(CPML_width_x), ier)
    call adios_write(handle, STRINGIFY_VAR(CPML_width_y), ier)
    call adios_write(handle, STRINGIFY_VAR(CPML_width_z), ier)
    if( nspec_cpml > 0 ) then
      local_dim = nspec_cpml_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(CPML_regions))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(CPML_to_spec))

      local_dim = nspec_ab_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(is_CPML))

      local_dim = NGLLX * NGLLY * NGLLZ * nspec_cpml_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(d_store_x))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(d_store_y))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(d_store_z))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(k_store_x))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(k_store_y))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(k_store_z))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(alpha_store))
      ! -----------------------------------------------------------------------
      ! for adjoint tomography
      ! save the array stored the points on interface between PML and interior
      ! computational domain
      ! -----------------------------------------------------------------------
      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) &
          .or. SIMULATION_TYPE == 3) then
        call adios_write(handle, &
                         STRINGIFY_VAR(nglob_interface_PML_acoustic), ier)
        call adios_write(handle, &
                         STRINGIFY_VAR(nglob_interface_PML_elastic), ier)
        if(nglob_interface_PML_acoustic > 0) then
          local_dim = nglob_interface_PML_acoustic_wmax
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                          local_dim,                  &
                                  STRINGIFY_VAR(points_interface_PML_acoustic))
        endif
        if(nglob_interface_PML_elastic > 0) then
          local_dim = nglob_interface_PML_elastic_wmax
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                           local_dim,                 &
                                  STRINGIFY_VAR(points_interface_PML_elastic))
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  call adios_write(handle, STRINGIFY_VAR(num_abs_boundary_faces), ier)
  if(PML_CONDITIONS)then
    if( num_abs_boundary_faces > 0 ) then
      local_dim = num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_ispec))
      local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_ijk))
      local_dim = NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_jacobian2Dw))
      local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_normal))
    endif
  else
    if( num_abs_boundary_faces > 0 ) then
      local_dim = num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_ispec))
      local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_ijk))
      local_dim = NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_jacobian2Dw))
      local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                       local_dim,                 &
                                       STRINGIFY_VAR(abs_boundary_normal))
      if( STACEY_ABSORBING_CONDITIONS ) then
        ! store mass matrix contributions
        if(ELASTIC_SIMULATION ) then
          local_dim = nglob_xy_wmax
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                           local_dim,                 &
                                           STRINGIFY_VAR(rmassx))
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                           local_dim,                 &
                                           STRINGIFY_VAR(rmassy))
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                           local_dim,                 &
                                           STRINGIFY_VAR(rmassz))
        endif
        if(ACOUSTIC_SIMULATION) then
          local_dim = nglob_xy_wmax
          call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                           local_dim,                 &
                                           STRINGIFY_VAR(rmassz_acoustic))
        endif
      endif
    endif
  endif

  call adios_write(handle, STRINGIFY_VAR(nspec2d_xmin), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_xmax), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_ymin), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_ymax), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_bottom), ier)
  call adios_write(handle, STRINGIFY_VAR(nspec2d_top), ier)

  if (nspec2d_xmin .ne. 0) then
      local_dim = nspec2d_xmin_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(ibelm_xmin))
  endif
  if (nspec2d_xmax .ne. 0) then
    local_dim = nspec2d_xmax_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(ibelm_xmax))
  endif
  if (nspec2d_ymin .ne. 0) then
    local_dim = nspec2d_ymin_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(ibelm_ymin))
  endif
  if (nspec2d_ymax .ne. 0) then
    local_dim = nspec2d_ymax_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(ibelm_ymax))
  endif
  if (nspec2d_bottom .ne. 0) then
    local_dim = nspec2d_bottom_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(ibelm_bottom))
  endif
  if (nspec2d_top .ne. 0) then
    local_dim = nspec2d_top_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(ibelm_top))
  endif

  ! free surface
  call adios_write(handle, STRINGIFY_VAR(num_free_surface_faces), ier)
  if( num_free_surface_faces > 0 ) then
    local_dim = num_free_surface_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(free_surface_ispec))
    local_dim = 3 * NGLLSQUARE * num_free_surface_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(free_surface_ijk))
    local_dim = NGLLSQUARE * num_free_surface_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(free_surface_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_free_surface_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(free_surface_normal))
  endif

  ! acoustic-elastic coupling surface
  call adios_write(handle, STRINGIFY_VAR(num_coupling_ac_el_faces), ier)
  if( num_coupling_ac_el_faces > 0 ) then
    local_dim = num_coupling_ac_el_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_el_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_el_ijk))
    local_dim = NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                STRINGIFY_VAR(coupling_ac_el_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_ac_el_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_el_normal))
  endif

  ! acoustic-poroelastic coupling surface
  call adios_write(handle, STRINGIFY_VAR(num_coupling_ac_po_faces), ier)
  if( num_coupling_ac_po_faces > 0 ) then
    local_dim = num_coupling_ac_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_po_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_po_ijk))
    local_dim = NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  STRINGIFY_VAR(coupling_ac_po_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_ac_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_ac_po_normal))
  endif

  ! elastic-poroelastic coupling surface
  call adios_write(handle, STRINGIFY_VAR(num_coupling_el_po_faces), ier)
  if( num_coupling_el_po_faces > 0 ) then
    local_dim = num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_el_po_ispec))
    local_dim = num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_po_el_ispec))
    local_dim = 3 * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_el_po_ijk))
    local_dim = 3 * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_po_el_ijk))
    local_dim = NGLLSQUARE * num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  STRINGIFY_VAR(coupling_el_po_jacobian2Dw))
    local_dim = NDIM * NGLLSQUARE * num_coupling_el_po_faces_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(coupling_el_po_normal))
  endif

  call adios_write(handle, STRINGIFY_VAR(num_interfaces_ext_mesh), ier)
  if( num_interfaces_ext_mesh > 0 ) then
    call adios_write(handle, STRINGIFY_VAR(max_nibool_interfaces_ext_mesh), ier)
    local_dim = num_interfaces_ext_mesh_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(my_neighbours_ext_mesh))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 STRINGIFY_VAR(nibool_interfaces_ext_mesh))
    local_dim = max_nibool_interfaces_ext_mesh_wmax &
              * num_interfaces_ext_mesh_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                             STRINGIFY_VAR(ibool_interfaces_ext_mesh_dummy))
  endif

! anisotropy
  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_aniso_wmax
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c11store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c12store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c13store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c14store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c15store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c16store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c22store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c23store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c24store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c25store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c26store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c33store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c34store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c35store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c36store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c44store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c45store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c46store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c55store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c56store))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(c66store))
  endif

  ! inner/outer elements
  local_dim = nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ispec_is_inner))

  if( ACOUSTIC_SIMULATION ) then
    call adios_write(handle, STRINGIFY_VAR(nspec_inner_acoustic), ier)
    call adios_write(handle, STRINGIFY_VAR(nspec_outer_acoustic), ier)
    call adios_write(handle, STRINGIFY_VAR(num_phase_ispec_acoustic), ier)
    if(num_phase_ispec_acoustic > 0 ) then
      local_dim = num_phase_ispec_acoustic_wmax * 2
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 STRINGIFY_VAR(phase_ispec_inner_acoustic))
    endif
  endif

  if( ELASTIC_SIMULATION ) then
    call adios_write(handle, STRINGIFY_VAR(nspec_inner_elastic), ier)
    call adios_write(handle, STRINGIFY_VAR(nspec_outer_elastic), ier)
    call adios_write(handle, STRINGIFY_VAR(num_phase_ispec_elastic), ier)
    if(num_phase_ispec_elastic > 0 ) then
      local_dim = num_phase_ispec_elastic_wmax * 2
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 STRINGIFY_VAR(phase_ispec_inner_elastic))
    endif
  endif

  if( POROELASTIC_SIMULATION ) then
    call adios_write(handle, STRINGIFY_VAR(nspec_inner_poroelastic), ier)
    call adios_write(handle, STRINGIFY_VAR(nspec_outer_poroelastic), ier)
    call adios_write(handle, STRINGIFY_VAR(num_phase_ispec_poroelastic), ier)
    if(num_phase_ispec_poroelastic > 0 ) then
      local_dim = num_phase_ispec_poroelastic_wmax * 2
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              STRINGIFY_VAR(phase_ispec_inner_poroelastic))
    endif
  endif

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    if( ACOUSTIC_SIMULATION ) then
      call adios_write(handle, STRINGIFY_VAR(num_colors_outer_acoustic), ier)
      call adios_write(handle, STRINGIFY_VAR(num_colors_inner_acoustic), ier)
      local_dim = num_colors_outer_acoustic_wmax &
                + num_colors_inner_acoustic_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              STRINGIFY_VAR(num_elem_colors_acoustic))
    endif
    if( ELASTIC_SIMULATION ) then
      call adios_write(handle, STRINGIFY_VAR(num_colors_outer_elastic), ier)
      call adios_write(handle, STRINGIFY_VAR(num_colors_inner_elastic), ier)
      local_dim = num_colors_outer_elastic_wmax &
                + num_colors_inner_elastic_wmax
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              STRINGIFY_VAR(num_elem_colors_elastic))
    endif
  endif

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, "", ier)
  call adios_close(handle, ier)

  ! stores arrays in binary files
  if( SAVE_MESH_FILES ) then
    call save_arrays_solver_files_adios(nspec,nglob,ibool, nspec_wmax, &
                                        nglob_wmax)

    ! debug: saves 1. MPI interface
    !if( num_interfaces_ext_mesh >= 1 ) then
    !  filename = prname(1:len_trim(prname))//'MPI_1_points'
    !  call write_VTK_data_points(nglob, &
    !                    xstore_dummy,ystore_dummy,zstore_dummy, &
    !                    ibool_interfaces_ext_mesh_dummy(1:nibool_interfaces_ext_mesh(1),1), &
    !                    nibool_interfaces_ext_mesh(1), &
    !                    filename)
    !endif
  endif

  ! cleanup
  deallocate(ibool_interfaces_ext_mesh_dummy,stat=ier)
  if( ier /= 0 ) stop 'error deallocating array ibool_interfaces_ext_mesh_dummy'

  if( nspec_cpml_tot > 0 ) then
     deallocate(CPML_to_spec,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_to_spec'
     deallocate(CPML_regions,stat=ier); if( ier /= 0 ) stop 'error deallocating array CPML_regions'
     deallocate(is_CPML,stat=ier); if( ier /= 0 ) stop 'error deallocating array is_CPML'
  endif

  if( PML_CONDITIONS ) then
     deallocate(d_store_x,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_x'
     deallocate(d_store_y,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_y'
     deallocate(d_store_z,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_z'
     deallocate(k_store_x,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_x'
     deallocate(k_store_y,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_y'
     deallocate(k_store_z,stat=ier); if( ier /= 0 ) stop 'error deallocating array d_store_z'
     deallocate(alpha_store,stat=ier); if( ier /= 0 ) stop 'error deallocating array alpha_store'
     if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
       deallocate(mask_ibool_interior_domain,stat=ier)
       if(ier /= 0) stop 'error deallocating array mask_ibool_interior_domain'

       if(nglob_interface_PML_acoustic > 0) then
         deallocate(points_interface_PML_acoustic,stat=ier)
         if( ier /= 0 ) stop 'error deallocating array points_interface_PML_acoustic'
       endif

       if(nglob_interface_PML_elastic > 0) then
         deallocate(points_interface_PML_elastic,stat=ier)
         if( ier /= 0 ) stop 'error deallocating array points_interface_PML_elastic'
       endif
     endif
  endif

end subroutine save_arrays_solver_ext_mesh_adios


!------------------------------------------------------------------------------

subroutine save_arrays_solver_files_adios(nspec,nglob,ibool, nspec_wmax, &
                                          nglob_wmax)

  use mpi
  use generate_databases_par, only: myrank, LOCAL_PATH,     &
                                    xstore, ystore, zstore, &
                                    sizeprocs
  use create_regions_mesh_ext_par
  use adios_helpers_mod

  implicit none

  integer :: nspec, nglob
  integer :: nspec_wmax, nglob_wmax
  ! mesh coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vp_tmp, &
                                                             vs_tmp, rho_tmp
  integer :: ier
  integer, dimension(:), allocatable :: iglob_tmp

  !--- Local parameters for ADIOS ---
  character(len=256) :: output_name
  integer(kind=8) :: group, handle
  integer(kind=8) :: groupsize, totalsize
  integer :: local_dim
  character(len=64), parameter :: group_name_coords  = "SPECFEM3D_MESH_COORDS"
  character(len=64), parameter :: group_name_values = "SPECFEM3D_MODEL_VALUES"

  if( myrank == 0) then
    write(IMAIN,*) '     saving mesh files for VisIt. ADIOS format'
    call flush_IMAIN()
  endif

#if 1
  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/mesh_coordinates.bp"
  call adios_declare_group(group, group_name_coords, "", 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(ngllx))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nglly))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(ngllz))

  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nglob))

  local_dim = nglob_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "x_global", xstore_dummy)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "y_global", ystore_dummy)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "z_global", zstore_dummy)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(xstore))
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(ystore))
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", STRINGIFY_VAR(ibool))

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call adios_open(handle, group_name_coords, output_name, "w", &
                  MPI_COMM_WORLD, ier);
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, STRINGIFY_VAR(ngllx), ier)
  call adios_write(handle, STRINGIFY_VAR(nglly), ier)
  call adios_write(handle, STRINGIFY_VAR(ngllz), ier)

  call adios_write(handle, STRINGIFY_VAR(nspec), ier)
  call adios_write(handle, STRINGIFY_VAR(nglob), ier)

  local_dim = nglob_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "x_global", xstore_dummy)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "y_global", ystore_dummy)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "z_global", zstore_dummy)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(xstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(ibool))

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, "", ier)
  call adios_close(handle, ier)
#endif
#if 1
  !----------------------------------.
  ! Set up the model values to write |
  !----------------------------------'
  allocate( vp_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if( ier /= 0 ) stop 'error allocating array '
  allocate( vs_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if( ier /= 0 ) stop 'error allocating array '
  allocate( rho_tmp(NGLLX,NGLLY,NGLLZ,nspec), stat=ier); if( ier /= 0 ) stop 'error allocating array '
  ! vp (for checking the mesh and model)
  !minimum = minval( abs(rho_vp) )
  !if( minimum(1) /= 0.0 ) then
  !  vp_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
  !else
  !  vp_tmp = 0.0
  !endif
  vp_tmp = 0.0
  where( rho_vp /= 0._CUSTOM_REAL ) vp_tmp = (FOUR_THIRDS * mustore + kappastore) / rho_vp
  ! vs (for checking the mesh and model)
  !minimum = minval( abs(rho_vs) )
  !if( minimum(1) /= 0.0 ) then
  !  vs_tmp = mustore / rho_vs
  !else
  !  vs_tmp = 0.0
  !endif
  vs_tmp = 0.0
  where( rho_vs /= 0._CUSTOM_REAL )  vs_tmp = mustore / rho_vs
  ! outputs density model for check
  rho_tmp = 0.0
  where( rho_vp /= 0._CUSTOM_REAL ) rho_tmp = rho_vp**2 / (FOUR_THIRDS * mustore + kappastore)

  !-----------------------------------.
  ! Setup ADIOS for the current group |
  !-----------------------------------'
  groupsize = 0
  output_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) // "/model_values.bp"
  call adios_declare_group(group, group_name_values, "", 1, ier)
  call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)

  !------------------------.
  ! Define ADIOS Variables |
  !------------------------'
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(ngllx))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nglly))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(ngllz))

  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nspec))
  call define_adios_scalar(group, groupsize, "", STRINGIFY_VAR(nglob))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "vp", vp_tmp)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "vs", vs_tmp)
  call define_adios_global_array1D(group, groupsize, local_dim, &
                                   "", "rho", rho_tmp)

  !------------------------------------------------------------.
  ! Open an handler to the ADIOS file and setup the group size |
  !------------------------------------------------------------'
  call adios_open(handle, group_name_values, output_name, "w", &
                  MPI_COMM_WORLD, ier);
  call adios_group_size (handle, groupsize, totalsize, ier)

  !------------------------------------------.
  ! Write previously defined ADIOS variables |
  !------------------------------------------'
  call adios_write(handle, STRINGIFY_VAR(ngllx), ier)
  call adios_write(handle, STRINGIFY_VAR(nglly), ier)
  call adios_write(handle, STRINGIFY_VAR(ngllz), ier)

  call adios_write(handle, STRINGIFY_VAR(nspec), ier)
  call adios_write(handle, STRINGIFY_VAR(nglob), ier)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_wmax
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "vp", vp_tmp)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "vs", vs_tmp)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   "rho", rho_tmp)

  !----------------------------------.
  ! Perform the actual write to disk |
  !----------------------------------'
  call adios_set_path(handle, "", ier)
  call adios_close(handle, ier)

  deallocate(vp_tmp)
  deallocate(vs_tmp)
  deallocate(rho_tmp)
#endif

end subroutine save_arrays_solver_files_adios
