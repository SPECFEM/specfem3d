!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine read_mesh_databases()

  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: database_name

  ! sets file name
  call create_name_database(prname,myrank,LOCAL_PATH)
  database_name = prname(1:len_trim(prname))//'external_mesh.bin'

! start reading the databases

! info about external mesh simulation
  if (I_should_read_the_database) then
    open(unit=27,file=trim(database_name),status='old', &
       action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open database file: ',trim(database_name)
      call exit_mpi(myrank,'Error opening database file')
    endif
  endif

  if (I_should_read_the_database) then
    read(27) NSPEC_AB
    read(27) NGLOB_AB

    read(27) ibool

    read(27) xstore
    read(27) ystore
    read(27) zstore

    read(27) xix
    read(27) xiy
    read(27) xiz
    read(27) etax
    read(27) etay
    read(27) etaz
    read(27) gammax
    read(27) gammay
    read(27) gammaz
    read(27) jacobian

    read(27) kappastore
    read(27) mustore

    read(27) ispec_is_acoustic
    read(27) ispec_is_elastic
    read(27) ispec_is_poroelastic
  endif

  call bcast_all_i_for_database(NSPEC_AB, 1)
  call bcast_all_i_for_database(NGLOB_AB, 1)
  call bcast_all_i_for_database(ibool(1,1,1,1), size(ibool))
  if (size(xstore) > 0) call bcast_all_cr_for_database(xstore(1), size(xstore))
  if (size(ystore) > 0) call bcast_all_cr_for_database(ystore(1), size(ystore))
  if (size(zstore) > 0) call bcast_all_cr_for_database(zstore(1), size(zstore))
  call bcast_all_cr_for_database(xix(1,1,1,1), size(xix))
  call bcast_all_cr_for_database(xiy(1,1,1,1), size(xiy))
  call bcast_all_cr_for_database(xiz(1,1,1,1), size(xiz))
  call bcast_all_cr_for_database(etax(1,1,1,1), size(etax))
  call bcast_all_cr_for_database(etay(1,1,1,1), size(etay))
  call bcast_all_cr_for_database(etaz(1,1,1,1), size(etaz))
  call bcast_all_cr_for_database(gammax(1,1,1,1), size(gammax))
  call bcast_all_cr_for_database(gammay(1,1,1,1), size(gammay))
  call bcast_all_cr_for_database(gammaz(1,1,1,1), size(gammaz))
  call bcast_all_cr_for_database(jacobian(1,1,1,1), size(jacobian))
  call bcast_all_cr_for_database(kappastore(1,1,1,1), size(kappastore))
  call bcast_all_cr_for_database(mustore(1,1,1,1), size(mustore))
  if (size(ispec_is_acoustic) > 0) &
    call bcast_all_l_for_database(ispec_is_acoustic(1), size(ispec_is_acoustic))
  if (size(ispec_is_elastic) > 0) &
    call bcast_all_l_for_database(ispec_is_elastic(1), size(ispec_is_elastic))
  if (size(ispec_is_poroelastic) > 0) &
    call bcast_all_l_for_database(ispec_is_poroelastic(1), size(ispec_is_poroelastic))

  ! acoustic
  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  ! all processes will have acoustic_simulation set if any flag is .true.
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
  if (ACOUSTIC_SIMULATION) then
    ! potentials
    allocate(potential_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array potential_acoustic'
    allocate(potential_dot_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array potential_dot_acoustic'
    allocate(potential_dot_dot_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array potential_dot_dot_acoustic'
    if (SIMULATION_TYPE /= 1) then
      allocate(potential_acoustic_adj_coupling(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array potential_acoustic_adj_coupling'
    endif
    ! mass matrix, density
    allocate(rmass_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rmass_acoustic'
    if (I_should_read_the_database) read(27) rmass_acoustic
    if (size(rmass_acoustic) > 0) call bcast_all_cr_for_database(rmass_acoustic(1), size(rmass_acoustic))

    ! initializes mass matrix contribution
    allocate(rmassz_acoustic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rmassz_acoustic'
    rmassz_acoustic(:) = 0._CUSTOM_REAL
  endif

! this array is needed for acoustic simulations but also for elastic simulations with CPML,
! thus we now allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array rhostore'
  if (I_should_read_the_database) read(27) rhostore
  call bcast_all_cr_for_database(rhostore(1,1,1,1), size(rhostore))

  ! elastic
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))

  ! elastic simulation
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  if (ELASTIC_SIMULATION) then
    ! displacement,velocity,acceleration
    allocate(displ(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array displ'
    allocate(veloc(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array veloc'
    allocate(accel(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array accel'
    if (SIMULATION_TYPE /= 1) then
      allocate(accel_adj_coupling(NDIM,NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array accel_adj_coupling'
    endif

    ! allocates mass matrix
    allocate(rmass(NGLOB_AB),stat=ier)

    if (ier /= 0) stop 'Error allocating array rmass'
    ! initializes mass matrix contributions
    allocate(rmassx(NGLOB_AB), &
             rmassy(NGLOB_AB), &
             rmassz(NGLOB_AB), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating array rmassx,rmassy,rmassz'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL

    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_vs'
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
    if (ier /= 0) stop 'Error allocating array c11store etc.'

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    allocate(R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS),stat=ier)
    if (ier /= 0) stop 'Error allocating array R_xx etc.'

    ! needed for attenuation and/or kernel computations
    allocate(epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) stop 'Error allocating array epsilondev_xx etc.'

    allocate(R_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array R_trace etc.'

    ! note: needed for some subroutine arguments
    allocate(epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array epsilon_trace_over_3'

    ! needed for attenuation
    allocate(one_minus_sum_beta(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB), &
             factor_common(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array one_minus_sum_beta etc.'

    allocate(one_minus_sum_beta_kappa(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB), &
             factor_common_kappa(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array one_minus_sum_beta_kappa etc.'

    ! reads mass matrices
    if (I_should_read_the_database) read(27,iostat=ier) rmass
    call bcast_all_cr_for_database(rmass(1), size(rmass))
    if (ier /= 0) stop 'Error reading in array rmass'

    if (APPROXIMATE_OCEAN_LOAD) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array rmass_ocean_load'
      if (I_should_read_the_database) read(27) rmass_ocean_load
      if (size(rmass_ocean_load) > 0) call bcast_all_cr_for_database(rmass_ocean_load(1), size(rmass_ocean_load))
    else
      ! dummy allocation
      allocate(rmass_ocean_load(1),stat=ier)
      if (ier /= 0) stop 'Error allocating dummy array rmass_ocean_load'
    endif

    !pll material parameters for stacey conditions
    if (I_should_read_the_database) read(27,iostat=ier) rho_vp
    if (size(rho_vp) > 0) call bcast_all_cr_for_database(rho_vp(1,1,1,1), size(rho_vp))
    if (ier /= 0) stop 'Error reading in array rho_vp'
    if (I_should_read_the_database) read(27,iostat=ier) rho_vs
    if (size(rho_vs) > 0) call bcast_all_cr_for_database(rho_vs(1,1,1,1), size(rho_vs))
    if (ier /= 0) stop 'Error reading in array rho_vs'

  else
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif

  ! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )
  if (POROELASTIC_SIMULATION) then

    if (GPU_MODE) call exit_mpi(myrank,'POROELASTICITY not supported by GPU mode yet...')

    ! displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(displs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array displs_poroelastic'
    allocate(velocs_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array velocs_poroelastic'
    allocate(accels_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array accels_poroelastic'
    allocate(displw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array displw_poroelastic'
    allocate(velocw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array velocw_poroelastic'
    allocate(accelw_poroelastic(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array accelw_poroelastic'

    allocate(rmass_solid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rmass_solid_poroelastic'
    allocate(rmass_fluid_poroelastic(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array rmass_fluid_poroelastic'

    allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             etastore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             phistore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
             rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array poroelastic properties'

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
    if (ier /= 0) stop 'Error allocating array epsilonsdev_xx etc.'

    allocate(epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array epsilons_trace_over_3 etc.'

    if (I_should_read_the_database) then
      read(27) rmass_solid_poroelastic
      read(27) rmass_fluid_poroelastic
      read(27) rhoarraystore
      read(27) kappaarraystore
      read(27) etastore
      read(27) tortstore
      read(27) permstore
      read(27) phistore
      read(27) rho_vpI
      read(27) rho_vpII
      read(27) rho_vsI
    endif
    if (size(rmass_solid_poroelastic) > 0) &
      call bcast_all_cr_for_database(rmass_solid_poroelastic(1), size(rmass_solid_poroelastic))
    if (size(rmass_fluid_poroelastic) > 0) &
      call bcast_all_cr_for_database(rmass_fluid_poroelastic(1), size(rmass_fluid_poroelastic))
    if (size(rhoarraystore) > 0) &
      call bcast_all_cr_for_database(rhoarraystore(1,1,1,1,1), size(rhoarraystore))
    if (size(kappaarraystore) > 0) &
      call bcast_all_cr_for_database(kappaarraystore(1,1,1,1,1), size(kappaarraystore))
    if (size(etastore) > 0) &
      call bcast_all_cr_for_database(etastore(1,1,1,1), size(etastore))
    if (size(tortstore) > 0) &
      call bcast_all_cr_for_database(tortstore(1,1,1,1), size(tortstore))
    if (size(permstore) > 0) &
      call bcast_all_cr_for_database(permstore(1,1,1,1,1), size(permstore))
    if (size(phistore) > 0) &
      call bcast_all_cr_for_database(phistore(1,1,1,1), size(phistore))
    if (size(rho_vpI) > 0) &
      call bcast_all_cr_for_database(rho_vpI(1,1,1,1), size(rho_vpI))
    if (size(rho_vpII) > 0) &
      call bcast_all_cr_for_database(rho_vpII(1,1,1,1), size(rho_vpII))
    if (size(rho_vsI) > 0) &
      call bcast_all_cr_for_database(rho_vsI(1,1,1,1), size(rho_vsI))
  endif

  ! checks simulation types are valid
  if ((.not. ACOUSTIC_SIMULATION) .and. &
      (.not. ELASTIC_SIMULATION) .and. &
      (.not. POROELASTIC_SIMULATION)) then
    if (I_should_read_the_database) close(27)
    call exit_mpi(myrank,'Error no simulation type defined')
  endif

  ! C-PML absorbing boundary conditions
  ! we allocate this array even when PMLs are absent because we need it in logical tests in "if" statements
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array is_CPML'

! make sure there are no PMLs by default,
! and then below if NSPEC_CPML > 0 we will read the real flags for this mesh from the disk
  is_CPML(:) = .false.
  NSPEC_CPML = 0

  if (PML_CONDITIONS) then
    if (I_should_read_the_database) then
      read(27) NSPEC_CPML
      read(27) CPML_width_x
      read(27) CPML_width_y
      read(27) CPML_width_z
      read(27) min_distance_between_CPML_parameter
    endif
    call bcast_all_i_for_database(NSPEC_CPML, 1)
    call bcast_all_cr_for_database(CPML_width_x, 1)
    call bcast_all_cr_for_database(CPML_width_y, 1)
    call bcast_all_cr_for_database(CPML_width_z, 1)
    call bcast_all_cr_for_database(min_distance_between_CPML_parameter, 1)

    if (NSPEC_CPML > 0) then
      allocate(CPML_regions(NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array CPML_regions'
      allocate(CPML_to_spec(NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array CPML_to_spec'

      allocate(d_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array d_store_x'
      allocate(d_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array d_store_y'
      allocate(d_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array d_store_z'
      allocate(K_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array K_store_x'
      allocate(K_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array K_store_y'
      allocate(K_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array K_store_z'
      allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array alpha_store'
      allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array alpha_store'
      allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if (ier /= 0) stop 'Error allocating array alpha_store'

      if (I_should_read_the_database) then
        read(27) CPML_regions
        read(27) CPML_to_spec
        read(27) is_CPML
        read(27) d_store_x
        read(27) d_store_y
        read(27) d_store_z
        read(27) k_store_x
        read(27) k_store_y
        read(27) k_store_z
        read(27) alpha_store_x
        read(27) alpha_store_y
        read(27) alpha_store_z
      endif
      if (size(CPML_regions) > 0) call bcast_all_i_for_database(CPML_regions(1), size(CPML_regions))
      if (size(CPML_to_spec) > 0) call bcast_all_i_for_database(CPML_to_spec(1), size(CPML_to_spec))
      if (size(is_CPML) > 0) call bcast_all_l_for_database(is_CPML(1), size(is_CPML))
      if (size(d_store_x) > 0) call bcast_all_cr_for_database(d_store_x(1,1,1,1), size(d_store_x))
      if (size(d_store_y) > 0) call bcast_all_cr_for_database(d_store_y(1,1,1,1), size(d_store_y))
      if (size(d_store_z) > 0) call bcast_all_cr_for_database(d_store_z(1,1,1,1), size(d_store_z))
      if (size(k_store_x) > 0) call bcast_all_cr_for_database(k_store_x(1,1,1,1), size(k_store_x))
      if (size(k_store_y) > 0) call bcast_all_cr_for_database(k_store_y(1,1,1,1), size(k_store_y))
      if (size(k_store_z) > 0) call bcast_all_cr_for_database(k_store_z(1,1,1,1), size(k_store_z))
      if (size(alpha_store_x) > 0) call bcast_all_cr_for_database(alpha_store_x(1,1,1,1), size(alpha_store_x))
      if (size(alpha_store_y) > 0) call bcast_all_cr_for_database(alpha_store_y(1,1,1,1), size(alpha_store_y))
      if (size(alpha_store_z) > 0) call bcast_all_cr_for_database(alpha_store_z(1,1,1,1), size(alpha_store_z))

      if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        if (I_should_read_the_database) read(27) nglob_interface_PML_acoustic
        call bcast_all_i_for_database(nglob_interface_PML_acoustic, 1)
        if (I_should_read_the_database) read(27) nglob_interface_PML_elastic
        call bcast_all_i_for_database(nglob_interface_PML_elastic, 1)
        if (nglob_interface_PML_acoustic > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if (ier /= 0) stop 'Error allocating array points_interface_PML_acoustic'
          if (I_should_read_the_database) read(27) points_interface_PML_acoustic
          if (size(points_interface_PML_acoustic) > 0) &
            call bcast_all_i_for_database(points_interface_PML_acoustic(1), size(points_interface_PML_acoustic))
        endif
        if (nglob_interface_PML_elastic > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if (ier /= 0) stop 'Error allocating array points_interface_PML_elastic'
          if (I_should_read_the_database) read(27) points_interface_PML_elastic
          if (size(points_interface_PML_elastic) > 0) &
            call bcast_all_i_for_database(points_interface_PML_elastic(1), size(points_interface_PML_elastic))
        endif
      endif
    endif
  endif

  ! absorbing boundary surface
  if (I_should_read_the_database) read(27) num_abs_boundary_faces
  call bcast_all_i_for_database(num_abs_boundary_faces, 1)

  ! checks
  if (num_abs_boundary_faces < 0) then
    print *,'read_mesh_databases: reading in negative num_abs_boundary_faces ',num_abs_boundary_faces,'...resetting to zero'
    num_abs_boundary_faces = 0
  endif
  allocate(abs_boundary_ispec(num_abs_boundary_faces), &
           abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

#ifdef DEBUG_COUPLED
    include "../../../add_to_read_mesh_databases.F90"
#endif

  if (num_abs_boundary_faces > 0) then
    if (I_should_read_the_database) then
      read(27) abs_boundary_ispec
      read(27) abs_boundary_ijk
      read(27) abs_boundary_jacobian2Dw
      read(27) abs_boundary_normal
    endif
    if (size(abs_boundary_ispec) > 0) &
      call bcast_all_i_for_database(abs_boundary_ispec(1), size(abs_boundary_ispec))
    if (size(abs_boundary_ijk) > 0) &
      call bcast_all_i_for_database(abs_boundary_ijk(1,1,1), size(abs_boundary_ijk))
    if (size(abs_boundary_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(abs_boundary_jacobian2Dw(1,1), size(abs_boundary_jacobian2Dw))
    if (size(abs_boundary_normal) > 0) &
      call bcast_all_cr_for_database(abs_boundary_normal(1,1,1), size(abs_boundary_normal))

    if (STACEY_ABSORBING_CONDITIONS .and. (.not. PML_CONDITIONS)) then
      ! store mass matrix contributions
      if (ELASTIC_SIMULATION) then
        if (I_should_read_the_database) then
          read(27) rmassx
          read(27) rmassy
          read(27) rmassz
        endif
        if (size(rmassx) > 0) call bcast_all_cr_for_database(rmassx(1), size(rmassx))
        if (size(rmassy) > 0) call bcast_all_cr_for_database(rmassy(1), size(rmassy))
        if (size(rmassz) > 0) call bcast_all_cr_for_database(rmassz(1), size(rmassz))
      endif
      if (ACOUSTIC_SIMULATION) then
        if (I_should_read_the_database) read(27) rmassz_acoustic
        if (size(rmassz_acoustic) > 0) call bcast_all_cr_for_database(rmassz_acoustic(1), size(rmassz_acoustic))
      endif
    endif
  endif

  if (I_should_read_the_database) then
    read(27) nspec2D_xmin
    read(27) nspec2D_xmax
    read(27) nspec2D_ymin
    read(27) nspec2D_ymax
    read(27) NSPEC2D_BOTTOM
    read(27) NSPEC2D_TOP
  endif
  call bcast_all_i_for_database(nspec2D_xmin, 1)
  call bcast_all_i_for_database(nspec2D_xmax, 1)
  call bcast_all_i_for_database(nspec2D_ymin, 1)
  call bcast_all_i_for_database(nspec2D_ymax, 1)
  call bcast_all_i_for_database(NSPEC2D_BOTTOM, 1)
  call bcast_all_i_for_database(NSPEC2D_TOP, 1)

  allocate(ibelm_xmin(nspec2D_xmin),ibelm_xmax(nspec2D_xmax), &
           ibelm_ymin(nspec2D_ymin),ibelm_ymax(nspec2D_ymax), &
           ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibelm_xmin,ibelm_xmax etc.'
  if (I_should_read_the_database) then
    read(27) ibelm_xmin
    read(27) ibelm_xmax
    read(27) ibelm_ymin
    read(27) ibelm_ymax
    read(27) ibelm_bottom
    read(27) ibelm_top
  endif
  if (size(ibelm_xmin) > 0) call bcast_all_i_for_database(ibelm_xmin(1), size(ibelm_xmin))
  if (size(ibelm_xmax) > 0) call bcast_all_i_for_database(ibelm_xmax(1), size(ibelm_xmax))
  if (size(ibelm_ymin) > 0) call bcast_all_i_for_database(ibelm_ymin(1), size(ibelm_ymin))
  if (size(ibelm_ymax) > 0) call bcast_all_i_for_database(ibelm_ymax(1), size(ibelm_ymax))
  if (size(ibelm_bottom) > 0) call bcast_all_i_for_database(ibelm_bottom(1), size(ibelm_bottom))
  if (size(ibelm_top) > 0) call bcast_all_i_for_database(ibelm_top(1), size(ibelm_top))

  ! free surface
  if (I_should_read_the_database) read(27) num_free_surface_faces
  call bcast_all_i_for_database(num_free_surface_faces, 1)
  allocate(free_surface_ispec(num_free_surface_faces), &
           free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
           free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
           free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays free_surface_ispec etc.'
  if (num_free_surface_faces > 0) then
    if (I_should_read_the_database) then
      read(27) free_surface_ispec
      read(27) free_surface_ijk
      read(27) free_surface_jacobian2Dw
      read(27) free_surface_normal
    endif
    if (size(free_surface_ispec) > 0) &
      call bcast_all_i_for_database(free_surface_ispec(1), size(free_surface_ispec))
    if (size(free_surface_ijk) > 0) &
      call bcast_all_i_for_database(free_surface_ijk(1,1,1), size(free_surface_ijk))
    if (size(free_surface_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(free_surface_jacobian2Dw(1,1), size(free_surface_jacobian2Dw))
    if (size(free_surface_normal) > 0) &
      call bcast_all_cr_for_database(free_surface_normal(1,1,1), size(free_surface_normal))
  endif

  ! acoustic-elastic coupling surface
  if (I_should_read_the_database) read(27) num_coupling_ac_el_faces
  call bcast_all_i_for_database(num_coupling_ac_el_faces, 1)
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces), &
           coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces), &
           coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces), &
           coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if (ier /= 0) stop 'Error allocating array coupling_ac_el_normal etc.'
  if (num_coupling_ac_el_faces > 0) then
    if (I_should_read_the_database) then
      read(27) coupling_ac_el_ispec
      read(27) coupling_ac_el_ijk
      read(27) coupling_ac_el_jacobian2Dw
      read(27) coupling_ac_el_normal
    endif
    if (size(coupling_ac_el_ispec) > 0) &
      call bcast_all_i_for_database(coupling_ac_el_ispec(1), size(coupling_ac_el_ispec))
    if (size(coupling_ac_el_ijk) > 0) &
      call bcast_all_i_for_database(coupling_ac_el_ijk(1,1,1), size(coupling_ac_el_ijk))
    if (size(coupling_ac_el_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_ac_el_jacobian2Dw(1,1), size(coupling_ac_el_jacobian2Dw))
    if (size(coupling_ac_el_normal) > 0) &
      call bcast_all_cr_for_database(coupling_ac_el_normal(1,1,1), size(coupling_ac_el_normal))
  endif

  ! acoustic-poroelastic coupling surface
  if (I_should_read_the_database) read(27) num_coupling_ac_po_faces
  call bcast_all_i_for_database(num_coupling_ac_po_faces, 1)
  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces), &
           coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces), &
           coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces), &
           coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if (ier /= 0) stop 'Error allocating array coupling_ac_po_normal etc.'
  if (num_coupling_ac_po_faces > 0) then
    if (I_should_read_the_database) then
      read(27) coupling_ac_po_ispec
      read(27) coupling_ac_po_ijk
      read(27) coupling_ac_po_jacobian2Dw
      read(27) coupling_ac_po_normal
    endif
    if (size(coupling_ac_po_ispec) > 0) &
      call bcast_all_i_for_database(coupling_ac_po_ispec(1), size(coupling_ac_po_ispec))
    if (size(coupling_ac_po_ijk) > 0) &
      call bcast_all_i_for_database(coupling_ac_po_ijk(1,1,1), size(coupling_ac_po_ijk))
    if (size(coupling_ac_po_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_ac_po_jacobian2Dw(1,1), size(coupling_ac_po_jacobian2Dw))
    if (size(coupling_ac_po_normal) > 0) &
      call bcast_all_cr_for_database(coupling_ac_po_normal(1,1,1), size(coupling_ac_po_normal))
  endif

  ! elastic-poroelastic coupling surface
  if (I_should_read_the_database) read(27) num_coupling_el_po_faces
  call bcast_all_i_for_database(num_coupling_el_po_faces, 1)
  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces), &
           coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces), &
           coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
           coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
           coupling_el_po_ispec(num_coupling_el_po_faces), &
           coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if (ier /= 0) stop 'Error allocating array coupling_el_po_normal etc.'
  if (num_coupling_el_po_faces > 0) then
    if (I_should_read_the_database) then
      read(27) coupling_el_po_ispec
      read(27) coupling_po_el_ispec
      read(27) coupling_el_po_ijk
      read(27) coupling_po_el_ijk
      read(27) coupling_el_po_jacobian2Dw
      read(27) coupling_el_po_normal
    endif
    if (size(coupling_el_po_ispec) > 0) &
      call bcast_all_i_for_database(coupling_el_po_ispec(1), size(coupling_el_po_ispec))
    if (size(coupling_po_el_ispec) > 0) &
      call bcast_all_i_for_database(coupling_po_el_ispec(1), size(coupling_po_el_ispec))
    if (size(coupling_el_po_ijk) > 0) &
      call bcast_all_i_for_database(coupling_el_po_ijk(1,1,1), size(coupling_el_po_ijk))
    if (size(coupling_po_el_ijk) > 0) &
      call bcast_all_i_for_database(coupling_po_el_ijk(1,1,1), size(coupling_po_el_ijk))
    if (size(coupling_el_po_jacobian2Dw) > 0) &
      call bcast_all_cr_for_database(coupling_el_po_jacobian2Dw(1,1), size(coupling_el_po_jacobian2Dw))
    if (size(coupling_el_po_normal) > 0) &
      call bcast_all_cr_for_database(coupling_el_po_normal(1,1,1), size(coupling_el_po_normal))
  endif

  ! MPI interfaces
  if (I_should_read_the_database) read(27) num_interfaces_ext_mesh
  call bcast_all_i_for_database(num_interfaces_ext_mesh, 1)
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh), &
           nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'Error allocating array my_neighbours_ext_mesh etc.'
  if (num_interfaces_ext_mesh > 0) then
    if (I_should_read_the_database) read(27) max_nibool_interfaces_ext_mesh
    call bcast_all_i_for_database(max_nibool_interfaces_ext_mesh, 1)
    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh'
    if (I_should_read_the_database) then
      read(27) my_neighbours_ext_mesh
      read(27) nibool_interfaces_ext_mesh
      read(27) ibool_interfaces_ext_mesh
    endif
    if (size(my_neighbours_ext_mesh) > 0) &
      call bcast_all_i_for_database(my_neighbours_ext_mesh(1), size(my_neighbours_ext_mesh))
    if (size(nibool_interfaces_ext_mesh) > 0) &
      call bcast_all_i_for_database(nibool_interfaces_ext_mesh(1), size(nibool_interfaces_ext_mesh))
    if (size(ibool_interfaces_ext_mesh) > 0) &
      call bcast_all_i_for_database(ibool_interfaces_ext_mesh(1,1), size(ibool_interfaces_ext_mesh))
  else
    max_nibool_interfaces_ext_mesh = 0
    allocate(ibool_interfaces_ext_mesh(0,0),stat=ier)
  endif

  if (ELASTIC_SIMULATION .and. ANISOTROPY) then
    if (I_should_read_the_database) then
      read(27) c11store
      read(27) c12store
      read(27) c13store
      read(27) c14store
      read(27) c15store
      read(27) c16store
      read(27) c22store
      read(27) c23store
      read(27) c24store
      read(27) c25store
      read(27) c26store
      read(27) c33store
      read(27) c34store
      read(27) c35store
      read(27) c36store
      read(27) c44store
      read(27) c45store
      read(27) c46store
      read(27) c55store
      read(27) c56store
      read(27) c66store
    endif
    if (size(c11store) > 0) call bcast_all_cr_for_database(c11store(1,1,1,1), size(c11store))
    if (size(c12store) > 0) call bcast_all_cr_for_database(c12store(1,1,1,1), size(c12store))
    if (size(c13store) > 0) call bcast_all_cr_for_database(c13store(1,1,1,1), size(c13store))
    if (size(c14store) > 0) call bcast_all_cr_for_database(c14store(1,1,1,1), size(c14store))
    if (size(c15store) > 0) call bcast_all_cr_for_database(c15store(1,1,1,1), size(c15store))
    if (size(c16store) > 0) call bcast_all_cr_for_database(c16store(1,1,1,1), size(c16store))
    if (size(c22store) > 0) call bcast_all_cr_for_database(c22store(1,1,1,1), size(c22store))
    if (size(c23store) > 0) call bcast_all_cr_for_database(c23store(1,1,1,1), size(c23store))
    if (size(c24store) > 0) call bcast_all_cr_for_database(c24store(1,1,1,1), size(c24store))
    if (size(c25store) > 0) call bcast_all_cr_for_database(c25store(1,1,1,1), size(c25store))
    if (size(c26store) > 0) call bcast_all_cr_for_database(c26store(1,1,1,1), size(c26store))
    if (size(c33store) > 0) call bcast_all_cr_for_database(c33store(1,1,1,1), size(c33store))
    if (size(c34store) > 0) call bcast_all_cr_for_database(c34store(1,1,1,1), size(c34store))
    if (size(c35store) > 0) call bcast_all_cr_for_database(c35store(1,1,1,1), size(c35store))
    if (size(c36store) > 0) call bcast_all_cr_for_database(c36store(1,1,1,1), size(c36store))
    if (size(c44store) > 0) call bcast_all_cr_for_database(c44store(1,1,1,1), size(c44store))
    if (size(c45store) > 0) call bcast_all_cr_for_database(c45store(1,1,1,1), size(c45store))
    if (size(c46store) > 0) call bcast_all_cr_for_database(c46store(1,1,1,1), size(c46store))
    if (size(c55store) > 0) call bcast_all_cr_for_database(c55store(1,1,1,1), size(c55store))
    if (size(c56store) > 0) call bcast_all_cr_for_database(c56store(1,1,1,1), size(c56store))
    if (size(c66store) > 0) call bcast_all_cr_for_database(c66store(1,1,1,1), size(c66store))
  endif

  ! inner / outer elements
  allocate(ispec_is_inner(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array ispec_is_inner'
  if (I_should_read_the_database) read(27) ispec_is_inner
  if (size(ispec_is_inner) > 0) call bcast_all_l_for_database(ispec_is_inner(1), size(ispec_is_inner))

  if (ACOUSTIC_SIMULATION) then
    if (I_should_read_the_database) then
      read(27) nspec_inner_acoustic,nspec_outer_acoustic
      read(27) num_phase_ispec_acoustic
    endif
    call bcast_all_i_for_database(nspec_inner_acoustic, 1)
    call bcast_all_i_for_database(nspec_outer_acoustic, 1)
    call bcast_all_i_for_database(num_phase_ispec_acoustic, 1)
    if (num_phase_ispec_acoustic < 0) stop 'Error acoustic simulation: num_phase_ispec_acoustic is < zero'
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_acoustic'
    if (num_phase_ispec_acoustic > 0) then
      if (I_should_read_the_database) read(27) phase_ispec_inner_acoustic
      if (size(phase_ispec_inner_acoustic) > 0) &
            call bcast_all_i_for_database(phase_ispec_inner_acoustic(1,1), size(phase_ispec_inner_acoustic))
    endif
  endif

  if (ELASTIC_SIMULATION) then
    if (I_should_read_the_database) then
      read(27) nspec_inner_elastic,nspec_outer_elastic
      read(27) num_phase_ispec_elastic
    endif
    call bcast_all_i_for_database(nspec_inner_elastic, 1)
    call bcast_all_i_for_database(nspec_outer_elastic, 1)
    call bcast_all_i_for_database(num_phase_ispec_elastic, 1)
    if (num_phase_ispec_elastic < 0) stop 'Error elastic simulation: num_phase_ispec_elastic is < zero'
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_elastic'
    if (num_phase_ispec_elastic > 0) then
      if (I_should_read_the_database) read(27) phase_ispec_inner_elastic
      if (size(phase_ispec_inner_elastic) > 0) &
        call bcast_all_i_for_database(phase_ispec_inner_elastic(1,1), size(phase_ispec_inner_elastic))
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    if (I_should_read_the_database) then
      read(27) nspec_inner_poroelastic,nspec_outer_poroelastic
      read(27) num_phase_ispec_poroelastic
    endif
    call bcast_all_i_for_database(nspec_inner_poroelastic, 1)
    call bcast_all_i_for_database(nspec_outer_poroelastic, 1)
    call bcast_all_i_for_database(num_phase_ispec_poroelastic, 1)
    if (num_phase_ispec_poroelastic < 0) stop 'Error poroelastic simulation: num_phase_ispec_poroelastic is < zero'
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0) stop 'Error allocating array phase_ispec_inner_poroelastic'
    if (num_phase_ispec_poroelastic > 0) then
      if (I_should_read_the_database) read(27) phase_ispec_inner_poroelastic
      if (size(phase_ispec_inner_poroelastic) > 0) &
        call bcast_all_i_for_database(phase_ispec_inner_poroelastic(1,1), size(phase_ispec_inner_poroelastic))
    endif
  endif

! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! acoustic domain colors
    if (ACOUSTIC_SIMULATION) then
      if (I_should_read_the_database) read(27) num_colors_outer_acoustic,num_colors_inner_acoustic
      call bcast_all_i_for_database(num_colors_outer_acoustic, 1)
      call bcast_all_i_for_database(num_colors_inner_acoustic, 1)

      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) stop 'Error allocating num_elem_colors_acoustic array'

      if (I_should_read_the_database) read(27) num_elem_colors_acoustic
      if (size(num_elem_colors_acoustic) > 0) &
        call bcast_all_i_for_database(num_elem_colors_acoustic(1), size(num_elem_colors_acoustic))
    endif
    ! elastic domain colors
    if (ELASTIC_SIMULATION) then
      if (I_should_read_the_database) read(27) num_colors_outer_elastic,num_colors_inner_elastic
      call bcast_all_i_for_database(num_colors_outer_elastic, 1)
      call bcast_all_i_for_database(num_colors_inner_elastic, 1)

      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) stop 'Error allocating num_elem_colors_elastic array'

      if (I_should_read_the_database) read(27) num_elem_colors_elastic
      if (size(num_elem_colors_elastic) > 0) &
        call bcast_all_i_for_database(num_elem_colors_elastic(1), size(num_elem_colors_elastic))
    endif
  else
    ! allocates dummy arrays
    if (ACOUSTIC_SIMULATION) then
      num_colors_outer_acoustic = 0
      num_colors_inner_acoustic = 0
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if (ier /= 0) stop 'Error allocating num_elem_colors_acoustic array'
    endif
    if (ELASTIC_SIMULATION) then
      num_colors_outer_elastic = 0
      num_colors_inner_elastic = 0
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) stop 'Error allocating num_elem_colors_elastic array'
    endif
  endif
  if (I_should_read_the_database) close(27)

  ! MPI communications
  if (ACOUSTIC_SIMULATION) then
    allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_scalar_ext_mesh(num_interfaces_ext_mesh), &
             request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating array buffer_send_scalar_ext_mesh,.. for acoustic simulations'
  endif
  if (ELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh(num_interfaces_ext_mesh), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating array buffer_send_vector_ext_mesh,.. for elastic simulations'
  endif
  if (POROELASTIC_SIMULATION) then
    allocate(buffer_send_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_send_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh_s(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh_s(num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh_w(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh_w(num_interfaces_ext_mesh), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating array buffer_send_vector_ext_mesh_s,.. for poroelastic simulations'
  endif

  end subroutine read_mesh_databases

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_moho()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  integer :: ier

  ! always needed to be allocated for routine arguments
  allocate( is_moho_top(NSPEC_BOUN),is_moho_bot(NSPEC_BOUN),stat=ier)
  if (ier /= 0) stop 'Error allocating array is_moho_top etc.'

  ! checks if anything to do
  if (ELASTIC_SIMULATION .and. SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ! boundary elements
    if (I_should_read_the_database) then
      open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='old', &
         action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open ibelm_moho file: ',prname(1:len_trim(prname))//'ibelm_moho.bin'
        call exit_mpi(myrank,'Error opening ibelm_moho file')
      endif
    endif

    if (I_should_read_the_database) read(27) NSPEC2D_MOHO
    call bcast_all_i_for_database(NSPEC2D_MOHO, 1)

    ! allocates arrays for moho mesh
    allocate(ibelm_moho_bot(NSPEC2D_MOHO), &
             ibelm_moho_top(NSPEC2D_MOHO), &
             normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
             normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
             ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO), &
             ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) stop 'Error allocating array ibelm_moho_bot etc.'

    if (I_should_read_the_database) then
      read(27) ibelm_moho_top
      read(27) ibelm_moho_bot
      read(27) ijk_moho_top
      read(27) ijk_moho_bot
    endif
    if (size(ibelm_moho_top) > 0) call bcast_all_i_for_database(ibelm_moho_top(1), size(ibelm_moho_top))
    if (size(ibelm_moho_bot) > 0) call bcast_all_i_for_database(ibelm_moho_bot(1), size(ibelm_moho_bot))
    if (size(ijk_moho_top) > 0) call bcast_all_i_for_database(ijk_moho_top(1,1,1), size(ijk_moho_top))
    if (size(ijk_moho_bot) > 0) call bcast_all_i_for_database(ijk_moho_bot(1,1,1), size(ijk_moho_bot))

    if (I_should_read_the_database) close(27)

    ! normals
    if (I_should_read_the_database) then
      open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='old', &
         action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open normal_moho file: ',prname(1:len_trim(prname))//'normal_moho.bin'
        call exit_mpi(myrank,'Error opening normal_moho file')
      endif
    endif

    if (I_should_read_the_database) then
      read(27) normal_moho_top
      read(27) normal_moho_bot
    endif
    if (size(normal_moho_top) > 0) call bcast_all_cr_for_database(normal_moho_top(1,1,1), size(normal_moho_top))
    if (size(normal_moho_bot) > 0) call bcast_all_cr_for_database(normal_moho_bot(1,1,1), size(normal_moho_bot))
    if (I_should_read_the_database) close(27)

    ! flags
    if (I_should_read_the_database) then
      open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='old', &
         action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open is_moho file: ',prname(1:len_trim(prname))//'is_moho.bin'
        call exit_mpi(myrank,'Error opening is_moho file')
      endif
    endif

    if (I_should_read_the_database) then
      read(27) is_moho_top
      read(27) is_moho_bot
    endif
    if (size(is_moho_top) > 0) call bcast_all_l_for_database(is_moho_top(1), size(is_moho_top))
    if (size(is_moho_bot) > 0) call bcast_all_l_for_database(is_moho_bot(1), size(is_moho_bot))

    if (I_should_read_the_database) close(27)

  else
    ! dummy
    NSPEC2D_MOHO = 1
  endif

  ! moho boundary
  if (ELASTIC_SIMULATION) then
    ! always needed to be allocated for routine arguments
    allocate(dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             b_dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             b_dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if (ier /= 0) stop 'Error allocating array dsdx_top etc.'
  endif

  end subroutine read_mesh_databases_moho


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_adjoint()

! reads in Moho meshes

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  integer :: ier

  ! allocates adjoint arrays for elastic simulations
  if (ELASTIC_SIMULATION .and. SIMULATION_TYPE == 3) then
    ! backward displacement,velocity,acceleration fields
    allocate(b_displ(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_displ'
    allocate(b_veloc(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_veloc'
    allocate(b_accel(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_accel'

    ! adjoint kernels

    ! primary, isotropic kernels
    ! density kernel
    allocate(rho_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_kl'

    if (ANISOTROPIC_KL) then
      ! anisotropic kernels
      allocate(cijkl_kl(21,NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array cijkl_kl'
      !dummy
      allocate(mu_kl(1,1,1,1))
      allocate(kappa_kl(1,1,1,1))
    else
      ! shear modulus kernel
      allocate(mu_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array mu_kl'
      ! compressional modulus kernel
      allocate(kappa_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array kappa_kl'
      !dummy
      allocate(cijkl_kl(1,1,1,1,1))
    endif

    ! noise source strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      allocate(sigma_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array sigma_kl'
    endif

    ! preconditioner
    if (APPROXIMATE_HESS_KL) then
      allocate(hess_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array hess_kl'
    else
      ! dummy allocation
      allocate(hess_kl(0,0,0,0),stat=ier)
      if (ier /= 0) stop 'Error allocating dummy array hess_kl'
    endif

    ! MPI handling
    allocate(b_request_send_vector_ext_mesh(num_interfaces_ext_mesh), &
             b_request_recv_vector_ext_mesh(num_interfaces_ext_mesh), &
             b_buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             b_buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_request_send_vector_ext_mesh etc.'

    ! allocates attenuation solids
    allocate(b_R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_R_xx etc.'

    ! note: these arrays are needed for attenuation and/or kernel computations
    allocate(b_epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             b_epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             b_epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             b_epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
             b_epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_epsilondev_xx etc.'
    ! needed for kernel computations
    allocate(b_epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_epsilon_trace_over_3'

    ! allocates attenuation solids for considering kappa
    allocate(b_R_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_R_trace etc.'

    ! Moho kernel
    if (SAVE_MOHO_MESH) then
      allocate( moho_kl(NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
      if (ier /= 0) stop 'Error allocating array moho_kl'
    endif
  else
    ! dummy allocation
    allocate(b_displ(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_displ'
    allocate(b_veloc(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_veloc'
    allocate(b_accel(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_accel'
  endif

  ! allocates adjoint arrays for acoustic simulations
  if (ACOUSTIC_SIMULATION .and. SIMULATION_TYPE == 3) then

    ! backward potentials
    allocate(b_potential_acoustic(NGLOB_ADJOINT), &
             b_potential_dot_acoustic(NGLOB_ADJOINT), &
             b_potential_dot_dot_acoustic(NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_potential_acoustic etc.'

    ! kernels
    allocate(rho_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             rhop_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             kappa_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             alpha_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array rho_ac_kl etc.'

    ! preconditioner
    if (APPROXIMATE_HESS_KL) then
      allocate(hess_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating array hess_ac_kl'
    else
      ! dummy allocation
      allocate(hess_ac_kl(0,0,0,0),stat=ier)
      if (ier /= 0) stop 'Error allocating dummy array hess_ac_kl'
    endif

    ! MPI handling
    allocate(b_request_send_scalar_ext_mesh(num_interfaces_ext_mesh), &
             b_request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), &
             b_buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_request_send_scalar_ext_mesh'

  else

    ! backward potentials
    allocate(b_potential_acoustic(1), &
             b_potential_dot_acoustic(1), &
             b_potential_dot_dot_acoustic(1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_potential_acoustic etc.'

    ! kernels
    allocate(rho_ac_kl(1,1,1,1), &
             rhop_ac_kl(1,1,1,1), &
             kappa_ac_kl(1,1,1,1), &
             alpha_ac_kl(1,1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array rho_ac_kl etc.'

    ! MPI handling
    allocate(b_request_send_scalar_ext_mesh(1), &
             b_request_recv_scalar_ext_mesh(1), &
             b_buffer_send_scalar_ext_mesh(1,1), &
             b_buffer_recv_scalar_ext_mesh(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_request_send_scalar_ext_mesh etc.'

  endif

  ! allocates adjoint arrays for poroelastic simulations
  if (POROELASTIC_SIMULATION .and. SIMULATION_TYPE == 3) then
    ! backward displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(b_displs_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_displs_poroelastic'
    allocate(b_velocs_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_velocs_poroelastic'
    allocate(b_accels_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_accels_poroelastic'
    allocate(b_displw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_displw_poroelastic'
    allocate(b_velocw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_velocw_poroelastic'
    allocate(b_accelw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_accelw_poroelastic'

    ! adjoint kernels

    ! primary, isotropic kernels
    allocate(rhot_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             rhof_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             sm_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             eta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array rhot_kl etc.'
    allocate(mufr_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array mufr_kl'
    allocate(B_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             C_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             M_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array B_kl etc.'

    ! density, isotropic kernels
    allocate(rhob_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             rhofb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             phi_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array rhob_kl etc.'
    allocate(mufrb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array mufrb_kl'
    allocate(Bb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             Cb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             Mb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array Bb_kl etc.'

    ! wavespeed, isotropic kernels
    allocate(rhobb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             rhofbb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             phib_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             ratio_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array rhobb_kl etc.'
    allocate(cs_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array cs_kl'
    allocate(cpI_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             cpII_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if (ier /= 0) stop 'Error allocating array cpI_kl etc.'

    ! MPI handling
    allocate(b_request_send_vector_ext_meshs(num_interfaces_ext_mesh), &
             b_request_recv_vector_ext_meshs(num_interfaces_ext_mesh), &
             b_buffer_send_vector_ext_meshs(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             b_buffer_recv_vector_ext_meshs(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_request_send_vector_ext_meshs etc.'

    allocate(b_request_send_vector_ext_meshw(num_interfaces_ext_mesh), &
             b_request_recv_vector_ext_meshw(num_interfaces_ext_mesh), &
             b_buffer_send_vector_ext_meshw(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             b_buffer_recv_vector_ext_meshw(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_request_send_vector_ext_meshw etc.'

    ! arrays needed for kernel computations
    allocate(b_epsilonsdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonsdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonsdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonsdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonsdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonwdev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonwdev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonwdev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonwdev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonwdev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_epsilonsdev_xx etc.'

    allocate(b_epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
             b_epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating array b_epsilons_trace_over_3 etc.'

  else ! dummy arrays

    ! backward displacement,velocity,acceleration for the solid (s) & fluid (w)
    ! phases
    allocate(b_displs_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_displs_poroelastic'
    allocate(b_velocs_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_velocs_poroelastic'
    allocate(b_accels_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_accels_poroelastic'
    allocate(b_displw_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_displw_poroelastic'
    allocate(b_velocw_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_velocw_poroelastic'
    allocate(b_accelw_poroelastic(1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array b_accelw_poroelastic'

    ! adjoint kernels

    ! primary, isotropic kernels
    allocate(rhot_kl(1,1,1,1), &
             rhof_kl(1,1,1,1), &
             sm_kl(1,1,1,1), &
             eta_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array rhot_kl etc.'
    allocate(mufr_kl(1,1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array mufr_kl'
    allocate(B_kl(1,1,1,1), &
             C_kl(1,1,1,1), &
             M_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array B_kl etc.'

    ! density, isotropic kernels
    allocate(rhob_kl(1,1,1,1), &
             rhofb_kl(1,1,1,1), &
             phi_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array rhob_kl etc.'
    allocate(mufrb_kl(1,1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array mufrb_kl'
    allocate(Bb_kl(1,1,1,1), &
             Cb_kl(1,1,1,1), &
             Mb_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array Bb_kl etc.'

    ! wavespeed, isotropic kernels
    allocate(rhobb_kl(1,1,1,1), &
             rhofbb_kl(1,1,1,1), &
             phib_kl(1,1,1,1), &
             ratio_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array rhobb_kl etc.'
    allocate(cs_kl(1,1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array cs_kl'
    allocate(cpI_kl(1,1,1,1), &
             cpII_kl(1,1,1,1), stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array cpI_kl etc.'

  endif

  end subroutine read_mesh_databases_adjoint

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_mesh_for_init()

! reads in the value of NSPEC_AB and NGLOB_AB

  use specfem_par

  implicit none
  ! Local variables
  integer :: ier
  character(len=MAX_STRING_LEN) :: database_name

  ! sets file name
  call create_name_database(prname,myrank,LOCAL_PATH)
  database_name = prname(1:len_trim(prname))//'external_mesh.bin'

  if (I_should_read_the_database) then
    open(unit=IIN,file=trim(database_name),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open database file: ',trim(database_name)
      call exit_mpi(myrank,'Error opening database file')
    endif
  endif

  if (I_should_read_the_database) then
    read(IIN) NSPEC_AB
    read(IIN) NGLOB_AB
  endif
  call bcast_all_i_for_database(NSPEC_AB, 1)
  call bcast_all_i_for_database(NGLOB_AB, 1)

  if (I_should_read_the_database) close(IIN)

  end subroutine read_mesh_for_init

