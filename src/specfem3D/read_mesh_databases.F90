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

  subroutine read_mesh_databases()

  use pml_par

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none

  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  integer :: ier,inum

! start reading the databases

! info about external mesh simulation
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=27,file=prname(1:len_trim(prname))//'external_mesh.bin',status='old',&
      action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error: could not open database '
    print*,'path: ',prname(1:len_trim(prname))//'external_mesh.bin'
    call exit_mpi(myrank,'error opening database')
  endif

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

  ! acoustic
  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  ! all processes will have acoustic_simulation set if any flag is .true.
  call any_all_l( ANY(ispec_is_acoustic), ACOUSTIC_SIMULATION )
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
    read(27) rmass_acoustic

    ! initializes mass matrix contribution
    allocate(rmassz_acoustic(NGLOB_AB),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rmassz_acoustic'
    rmassz_acoustic(:) = 0._CUSTOM_REAL
  endif

! this array is needed for acoustic simulations but also for elastic simulations with CPML,
! thus we now allocate it and read it in all cases (whether the simulation is acoustic, elastic, or acoustic/elastic)
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array rhostore'
  read(27) rhostore

  ! elastic
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))

  ! elastic simulation
  call any_all_l( ANY(ispec_is_elastic), ELASTIC_SIMULATION )
  if( ELASTIC_SIMULATION ) then
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

    ! note: currently, they need to be defined, as they are used in the routine arguments
    !          for compute_forces_viscoelastic_Deville()
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

    ! reads mass matrices
    read(27,iostat=ier) rmass
    if( ier /= 0 ) stop 'error reading in array rmass'

    if( APPROXIMATE_OCEAN_LOAD ) then
      ! ocean mass matrix
      allocate(rmass_ocean_load(NGLOB_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating array rmass_ocean_load'
      read(27) rmass_ocean_load
    else
      ! dummy allocation
      allocate(rmass_ocean_load(1),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy array rmass_ocean_load'
    endif

    !pll material parameters for stacey conditions
    read(27,iostat=ier) rho_vp
    if( ier /= 0 ) stop 'error reading in array rho_vp'
    read(27,iostat=ier) rho_vs
    if( ier /= 0 ) stop 'error reading in array rho_vs'

!! DK DK rhostore is now allocated and read in all cases (see above)
!   ! checks if rhostore is available for gravity
!   if( GRAVITY ) then
!
!     if( .not. ACOUSTIC_SIMULATION ) then
!       ! rho array needed for gravity
!       allocate(rhostore(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
!       if( ier /= 0 ) stop 'error allocating array rhostore'
!
!       ! extract rho information from mu = rho * vs * vs and rho_vs = rho * vs
!       rhostore = 0.0_CUSTOM_REAL
!       where( mustore > TINYVAL )
!         rhostore = (rho_vs*rho_vs) / mustore
!       endwhere
!
!       ! note: the construct below leads to a segmentation fault (ifort v11.1). not sure why...
!       !          (where statement - standard fortran 95)
!       !where( mustore > TINYVAL )
!       !  rhostore = (rho_vs*rho_vs) / mustore
!       !elsewhere
!       !  rhostore = 0.0_CUSTOM_REAL
!       !endwhere
!     endif
!   endif
  else
    ! no elastic attenuation & anisotropy
    ATTENUATION = .false.
    ANISOTROPY = .false.
  endif

  ! poroelastic
  call any_all_l( ANY(ispec_is_poroelastic), POROELASTIC_SIMULATION )
  if( POROELASTIC_SIMULATION ) then

    if( GPU_MODE ) call exit_mpi(myrank,'POROELASTICITY not supported by GPU mode yet...')

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

  ! checks simulation types are valid
  if( (.not. ACOUSTIC_SIMULATION ) .and. &
     (.not. ELASTIC_SIMULATION ) .and. &
     (.not. POROELASTIC_SIMULATION ) ) then
     close(27)
     call exit_mpi(myrank,'error no simulation type defined')
  endif

  ! C-PML absorbing boundary conditions
  NSPEC_CPML = 0
  if( PML_CONDITIONS ) then
    read(27) NSPEC_CPML
    read(27) CPML_width_x
    read(27) CPML_width_y
    read(27) CPML_width_z

    allocate(is_CPML(NSPEC_AB),stat=ier)
    if(ier /= 0) stop 'error allocating array is_CPML'

! make sure there are no PMLs by default,
! and then below if NSPEC_CPML > 0 we will need the real flags for this mesh from the disk
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
      allocate(alpha_store_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array alpha_store'
      allocate(alpha_store_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array alpha_store'
      allocate(alpha_store_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
      if(ier /= 0) stop 'error allocating array alpha_store'

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

      if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
        read(27) nglob_interface_PML_acoustic
        read(27) nglob_interface_PML_elastic
        if(nglob_interface_PML_acoustic > 0) then
          allocate(points_interface_PML_acoustic(nglob_interface_PML_acoustic),stat=ier)
          if(ier /= 0) stop 'error allocating array points_interface_PML_acoustic'
          read(27) points_interface_PML_acoustic
        endif
        if(nglob_interface_PML_elastic > 0) then
          allocate(points_interface_PML_elastic(nglob_interface_PML_elastic),stat=ier)
          if(ier /= 0) stop 'error allocating array points_interface_PML_elastic'
          read(27) points_interface_PML_elastic
        endif
      endif
    endif
  else
    ! allocate with a dummy size of zero just to be able to use this array as argument in subroutine calls
    allocate(is_CPML(0),stat=ier)
  endif

  ! absorbing boundary surface
  read(27) num_abs_boundary_faces

  ! checks
  if( num_abs_boundary_faces < 0 ) then
    print*,'read_mesh_databases: reading in negative num_abs_boundary_faces ',num_abs_boundary_faces,'...resetting to zero'
    num_abs_boundary_faces = 0
  endif
  allocate(abs_boundary_ispec(num_abs_boundary_faces), &
          abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
          abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
          abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array abs_boundary_ispec etc.'

  if (OLD_TEST_TO_FIX_ONE_DAY) then
     ! VM for new method
!! DK DK for VM VM: these two arrays are undeclared, thus I comment them out for now otherwise the code does not compile
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

  if(PML_CONDITIONS)then
     if( num_abs_boundary_faces > 0 ) then
       read(27) abs_boundary_ispec
       read(27) abs_boundary_ijk
       read(27) abs_boundary_jacobian2Dw
       read(27) abs_boundary_normal
     endif
  else
     if( num_abs_boundary_faces > 0 ) then
       read(27) abs_boundary_ispec
       read(27) abs_boundary_ijk
       read(27) abs_boundary_jacobian2Dw
       read(27) abs_boundary_normal
       if( STACEY_ABSORBING_CONDITIONS ) then
          ! store mass matrix contributions
          if(ELASTIC_SIMULATION) then
           read(27) rmassx
           read(27) rmassy
           read(27) rmassz
          endif
          if(ACOUSTIC_SIMULATION) then
           read(27) rmassz_acoustic
          endif
       endif
     endif
  endif

  read(27) nspec2D_xmin
  read(27) nspec2D_xmax
  read(27) nspec2D_ymin
  read(27) nspec2D_ymax
  read(27) NSPEC2D_BOTTOM
  read(27) NSPEC2D_TOP

  allocate(ibelm_xmin(nspec2D_xmin),ibelm_xmax(nspec2D_xmax), &
       ibelm_ymin(nspec2D_ymin),ibelm_ymax(nspec2D_ymax), &
       ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error allocating arrays ibelm_xmin,ibelm_xmax etc.'
  read(27) ibelm_xmin
  read(27) ibelm_xmax
  read(27) ibelm_ymin
  read(27) ibelm_ymax
  read(27) ibelm_bottom
  read(27) ibelm_top

  ! free surface
  read(27) num_free_surface_faces
  allocate(free_surface_ispec(num_free_surface_faces), &
          free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
          free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
          free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if(ier /= 0) stop 'error allocating arrays free_surface_ispec etc.'
  if( num_free_surface_faces > 0 ) then
    read(27) free_surface_ispec
    read(27) free_surface_ijk
    read(27) free_surface_jacobian2Dw
    read(27) free_surface_normal
  endif

  ! acoustic-elastic coupling surface
  read(27) num_coupling_ac_el_faces
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces), &
          coupling_ac_el_ispec(num_coupling_ac_el_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_el_normal etc.'
  if( num_coupling_ac_el_faces > 0 ) then
    read(27) coupling_ac_el_ispec
    read(27) coupling_ac_el_ijk
    read(27) coupling_ac_el_jacobian2Dw
    read(27) coupling_ac_el_normal
  endif

  ! acoustic-poroelastic coupling surface
  read(27) num_coupling_ac_po_faces
  allocate(coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces), &
          coupling_ac_po_ispec(num_coupling_ac_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_ac_po_normal etc.'
  if( num_coupling_ac_po_faces > 0 ) then
    read(27) coupling_ac_po_ispec
    read(27) coupling_ac_po_ijk
    read(27) coupling_ac_po_jacobian2Dw
    read(27) coupling_ac_po_normal
  endif

  ! elastic-poroelastic coupling surface
  read(27) num_coupling_el_po_faces
  allocate(coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces), &
          coupling_el_po_ispec(num_coupling_el_po_faces), &
          coupling_po_el_ispec(num_coupling_el_po_faces),stat=ier)
  if( ier /= 0 ) stop 'error allocating array coupling_el_po_normal etc.'
  if( num_coupling_el_po_faces > 0 ) then
    read(27) coupling_el_po_ispec
    read(27) coupling_po_el_ispec
    read(27) coupling_el_po_ijk
    read(27) coupling_po_el_ijk
    read(27) coupling_el_po_jacobian2Dw
    read(27) coupling_el_po_normal
  endif

  ! MPI interfaces
  read(27) num_interfaces_ext_mesh
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh), &
          nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array my_neighbours_ext_mesh etc.'
  if( num_interfaces_ext_mesh > 0 ) then
    read(27) max_nibool_interfaces_ext_mesh
    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array ibool_interfaces_ext_mesh'
    read(27) my_neighbours_ext_mesh
    read(27) nibool_interfaces_ext_mesh
    read(27) ibool_interfaces_ext_mesh
  else
    max_nibool_interfaces_ext_mesh = 0
    allocate(ibool_interfaces_ext_mesh(0,0),stat=ier)
  endif

  if( ELASTIC_SIMULATION .and. ANISOTROPY ) then
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

  ! inner / outer elements
  allocate(ispec_is_inner(NSPEC_AB),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ispec_is_inner'
  read(27) ispec_is_inner

  if( ACOUSTIC_SIMULATION ) then
    read(27) nspec_inner_acoustic,nspec_outer_acoustic
    read(27) num_phase_ispec_acoustic
    if( num_phase_ispec_acoustic < 0 ) stop 'error acoustic simulation: num_phase_ispec_acoustic is < zero'
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_acoustic'
    if(num_phase_ispec_acoustic > 0 ) read(27) phase_ispec_inner_acoustic
  endif

  if( ELASTIC_SIMULATION ) then
    read(27) nspec_inner_elastic,nspec_outer_elastic
    read(27) num_phase_ispec_elastic
    if( num_phase_ispec_elastic < 0 ) stop 'error elastic simulation: num_phase_ispec_elastic is < zero'
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_elastic'
    if(num_phase_ispec_elastic > 0 ) read(27) phase_ispec_inner_elastic
  endif

  if( POROELASTIC_SIMULATION ) then
    read(27) nspec_inner_poroelastic,nspec_outer_poroelastic
    read(27) num_phase_ispec_poroelastic
    if( num_phase_ispec_poroelastic < 0 ) stop 'error poroelastic simulation: num_phase_ispec_poroelastic is < zero'
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_poroelastic'
    if(num_phase_ispec_poroelastic > 0 ) read(27) phase_ispec_inner_poroelastic
  endif

! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! acoustic domain colors
    if( ACOUSTIC_SIMULATION ) then
      read(27) num_colors_outer_acoustic,num_colors_inner_acoustic

      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'

      read(27) num_elem_colors_acoustic
    endif
    ! elastic domain colors
    if( ELASTIC_SIMULATION ) then
      read(27) num_colors_outer_elastic,num_colors_inner_elastic

      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'

      read(27) num_elem_colors_elastic
    endif
  else
    ! allocates dummy arrays
    if( ACOUSTIC_SIMULATION ) then
      num_colors_outer_acoustic = 0
      num_colors_inner_acoustic = 0
      allocate(num_elem_colors_acoustic(num_colors_outer_acoustic + num_colors_inner_acoustic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_acoustic array'
    endif
    if( ELASTIC_SIMULATION ) then
      num_colors_outer_elastic = 0
      num_colors_inner_elastic = 0
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if( ier /= 0 ) stop 'error allocating num_elem_colors_elastic array'
    endif
  endif
  close(27)

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

  ! debug
  !call sum_all_i(num_interfaces_ext_mesh,inum)
  !if(myrank == 0) then
  !  write(IMAIN,*) 'number of MPI partition interfaces: ',inum
  !  write(IMAIN,*)
  !endif

  ! MPI communications
  if( ACOUSTIC_SIMULATION ) then
    allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_scalar_ext_mesh(num_interfaces_ext_mesh), &
             request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), &
             stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_scalar_ext_mesh,.. for acoustic simulations'
  endif
  if( ELASTIC_SIMULATION ) then
    allocate(buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh(num_interfaces_ext_mesh), &
             stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_vector_ext_mesh,.. for elastic simulations'
  endif
  if( POROELASTIC_SIMULATION ) then
    allocate(buffer_send_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh_s(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_send_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             buffer_recv_vector_ext_mesh_w(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh_s(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh_s(num_interfaces_ext_mesh), &
             request_send_vector_ext_mesh_w(num_interfaces_ext_mesh), &
             request_recv_vector_ext_mesh_w(num_interfaces_ext_mesh), &
             stat=ier)
    if( ier /= 0 ) stop 'error allocating array buffer_send_vector_ext_mesh_s,.. for poroelastic simulations'
  endif

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

  end subroutine read_mesh_databases


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_adjoint()

! reads in moho meshes

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none

  integer :: ier

  ! allocates adjoint arrays for elastic simulations
  if( ELASTIC_SIMULATION .and. SIMULATION_TYPE == 3 ) then
    ! backward displacement,velocity,acceleration fields
    allocate(b_displ(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_displ'
    allocate(b_veloc(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_veloc'
    allocate(b_accel(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_accel'

    ! adjoint kernels

    ! primary, isotropic kernels
    ! density kernel
    allocate(rho_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_kl'

    if (ANISOTROPIC_KL) then
      ! anisotropic kernels
      allocate(cijkl_kl(21,NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array cijkl_kl'
      !dummy
      allocate(mu_kl(1,1,1,1))
      allocate(kappa_kl(1,1,1,1))
    else
      ! shear modulus kernel
      allocate(mu_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array mu_kl'
      ! compressional modulus kernel
      allocate(kappa_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array kappa_kl'
      !dummy
      allocate(cijkl_kl(1,1,1,1,1))
    endif

    ! noise source strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      allocate(sigma_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array sigma_kl'
    endif

    ! preconditioner
    if ( APPROXIMATE_HESS_KL ) then
      allocate(hess_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array hess_kl'
    else
      ! dummy allocation
      allocate(hess_kl(0,0,0,0),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy array hess_kl'
    endif

    ! MPI handling
    allocate(b_request_send_vector_ext_mesh(num_interfaces_ext_mesh), &
             b_request_recv_vector_ext_mesh(num_interfaces_ext_mesh), &
             b_buffer_send_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
             b_buffer_recv_vector_ext_mesh(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_request_send_vector_ext_mesh etc.'

    ! allocates attenuation solids
    allocate(b_R_xx(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_yy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_xy(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_xz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS), &
             b_R_yz(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_R_xx etc.'

    ! note: these arrays are needed for attenuation and/or kernel computations
    allocate(b_epsilondev_xx(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            b_epsilondev_yy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            b_epsilondev_xy(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            b_epsilondev_xz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY), &
            b_epsilondev_yz(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_epsilon_dev_xx etc.'
    ! needed for kernel computations
    allocate(b_epsilon_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_epsilon_trace_over_3'

    ! allocates attenuation solids for considering kappa
    allocate(b_R_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa,N_SLS),&
             b_epsilondev_trace(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_R_trace etc.'

  else
    ! modification: Camille Mazoyer
    ! dummy allocation
    allocate(b_displ(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_displ'
    allocate(b_veloc(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_veloc'
    allocate(b_accel(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_accel'

  endif

  ! allocates adjoint arrays for acoustic simulations
  if( ACOUSTIC_SIMULATION .and. SIMULATION_TYPE == 3 ) then

    ! backward potentials
    allocate(b_potential_acoustic(NGLOB_ADJOINT), &
            b_potential_dot_acoustic(NGLOB_ADJOINT), &
            b_potential_dot_dot_acoustic(NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_potential_acoustic etc.'

    ! kernels
    allocate(rho_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            rhop_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            kappa_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            alpha_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array rho_ac_kl etc.'

    ! preconditioner
    if ( APPROXIMATE_HESS_KL ) then
      allocate(hess_ac_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
      if( ier /= 0 ) stop 'error allocating array hess_ac_kl'
    else
      ! dummy allocation
      allocate(hess_ac_kl(0,0,0,0),stat=ier)
      if( ier /= 0 ) stop 'error allocating dummy array hess_ac_kl'
    endif

    ! MPI handling
    allocate(b_request_send_scalar_ext_mesh(num_interfaces_ext_mesh), &
      b_request_recv_scalar_ext_mesh(num_interfaces_ext_mesh), &
      b_buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
      b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_request_send_scalar_ext_mesh'

  else

    ! backward potentials
    allocate(b_potential_acoustic(1), &
            b_potential_dot_acoustic(1), &
            b_potential_dot_dot_acoustic(1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_potential_acoustic etc.'

    ! kernels
    allocate(rho_ac_kl(1,1,1,1), &
            rhop_ac_kl(1,1,1,1), &
            kappa_ac_kl(1,1,1,1), &
            alpha_ac_kl(1,1,1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array rho_ac_kl etc.'

    ! MPI handling
    allocate(b_request_send_scalar_ext_mesh(1), &
            b_request_recv_scalar_ext_mesh(1), &
            b_buffer_send_scalar_ext_mesh(1,1), &
            b_buffer_recv_scalar_ext_mesh(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_request_send_scalar_ext_mesh etc.'

  endif

  ! ADJOINT moho
  ! moho boundary
  if( ELASTIC_SIMULATION ) then
    ! always needed to be allocated for routine arguments
    allocate( is_moho_top(NSPEC_BOUN),is_moho_bot(NSPEC_BOUN),stat=ier)
    if( ier /= 0 ) stop 'error allocating array is_moho_top etc.'

    if( SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3 ) then

      if (ADIOS_FOR_MESH) then
        call read_moho_mesh_adjoint_adios()
      else
        ! boundary elements
        !open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='unknown',form='unformatted')
        open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='old',&
              form='unformatted',iostat=ier)
        if( ier /= 0 ) then
          print*,'error: could not open ibelm_moho '
          print*,'path: ',prname(1:len_trim(prname))//'ibelm_moho.bin'
          call exit_mpi(myrank,'error opening ibelm_moho')
        endif

        read(27) NSPEC2D_MOHO

        ! allocates arrays for moho mesh
        allocate(ibelm_moho_bot(NSPEC2D_MOHO), &
                ibelm_moho_top(NSPEC2D_MOHO), &
                normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
                normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
                ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO), &
                ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
        if( ier /= 0 ) stop 'error allocating array ibelm_moho_bot etc.'

        read(27) ibelm_moho_top
        read(27) ibelm_moho_bot
        read(27) ijk_moho_top
        read(27) ijk_moho_bot

        close(27)

        ! normals
        open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='old',&
              form='unformatted',iostat=ier)
        if( ier /= 0 ) then
          print*,'error: could not open normal_moho '
          print*,'path: ',prname(1:len_trim(prname))//'normal_moho.bin'
          call exit_mpi(myrank,'error opening normal_moho')
        endif

        read(27) normal_moho_top
        read(27) normal_moho_bot
        close(27)

        ! flags
        open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='old',&
              form='unformatted',iostat=ier)
        if( ier /= 0 ) then
          print*,'error: could not open is_moho '
          print*,'path: ',prname(1:len_trim(prname))//'is_moho.bin'
          call exit_mpi(myrank,'error opening is_moho')
        endif

        read(27) is_moho_top
        read(27) is_moho_bot

        close(27)
      endif

      ! moho kernel
      allocate( moho_kl(NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
      if( ier /= 0 ) stop 'error allocating array moho_kl'
      moho_kl = 0._CUSTOM_REAL

    else
      NSPEC2D_MOHO = 1
    endif

    ! always needed to be allocated for routine arguments
    allocate( dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             b_dsdx_top(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO), &
             b_dsdx_bot(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO),stat=ier)
    if( ier /= 0 ) stop 'error allocating array dsdx_top etc.'
  endif

  ! allocates adjoint arrays for poroelastic simulations
  if( POROELASTIC_SIMULATION .and. SIMULATION_TYPE == 3 ) then
    ! backward displacement,velocity,acceleration for the solid (s) & fluid (w) phases
    allocate(b_displs_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_displs_poroelastic'
    allocate(b_velocs_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_velocs_poroelastic'
    allocate(b_accels_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_accels_poroelastic'
    allocate(b_displw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_displw_poroelastic'
    allocate(b_velocw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_velocw_poroelastic'
    allocate(b_accelw_poroelastic(NDIM,NGLOB_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_accelw_poroelastic'

    ! adjoint kernels

    ! primary, isotropic kernels
    allocate(rhot_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            rhof_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            sm_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            eta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array rhot_kl etc.'
    allocate(mufr_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array mufr_kl'
    allocate(B_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            C_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            M_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array B_kl etc.'

    ! density, isotropic kernels
    allocate(rhob_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            rhofb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            phi_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array rhob_kl etc.'
    allocate(mufrb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array mufrb_kl'
    allocate(Bb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            Cb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            Mb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array Bb_kl etc.'

    ! wavespeed, isotropic kernels
    allocate(rhobb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            rhofbb_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            phib_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            ratio_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array rhobb_kl etc.'
    allocate(cs_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array cs_kl'
    allocate(cpI_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            cpII_kl(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), stat=ier)
    if( ier /= 0 ) stop 'error allocating array cpI_kl etc.'

    ! MPI handling
    allocate(b_request_send_vector_ext_meshs(num_interfaces_ext_mesh), &
      b_request_recv_vector_ext_meshs(num_interfaces_ext_mesh), &
      b_buffer_send_vector_ext_meshs(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
      b_buffer_recv_vector_ext_meshs(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_request_send_vector_ext_meshs etc.'

    allocate(b_request_send_vector_ext_meshw(num_interfaces_ext_mesh), &
      b_request_recv_vector_ext_meshw(num_interfaces_ext_mesh), &
      b_buffer_send_vector_ext_meshw(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh), &
      b_buffer_recv_vector_ext_meshw(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_request_send_vector_ext_meshw etc.'

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
    if( ier /= 0 ) stop 'error allocating array b_epsilonsdev_xx etc.'

    allocate(b_epsilons_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT), &
            b_epsilonw_trace_over_3(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),stat=ier)
    if( ier /= 0 ) stop 'error allocating array b_epsilons_trace_over_3 etc.'

  else ! dummy arrays

    ! backward displacement,velocity,acceleration for the solid (s) & fluid (w)
    ! phases
    allocate(b_displs_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_displs_poroelastic'
    allocate(b_velocs_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_velocs_poroelastic'
    allocate(b_accels_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_accels_poroelastic'
    allocate(b_displw_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_displw_poroelastic'
    allocate(b_velocw_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_velocw_poroelastic'
    allocate(b_accelw_poroelastic(1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array b_accelw_poroelastic'

    ! adjoint kernels

    ! primary, isotropic kernels
    allocate(rhot_kl(1,1,1,1), &
            rhof_kl(1,1,1,1), &
            sm_kl(1,1,1,1), &
            eta_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array rhot_kl etc.'
    allocate(mufr_kl(1,1,1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array mufr_kl'
    allocate(B_kl(1,1,1,1), &
            C_kl(1,1,1,1), &
            M_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array B_kl etc.'

    ! density, isotropic kernels
    allocate(rhob_kl(1,1,1,1), &
            rhofb_kl(1,1,1,1), &
            phi_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array rhob_kl etc.'
    allocate(mufrb_kl(1,1,1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array mufrb_kl'
    allocate(Bb_kl(1,1,1,1), &
            Cb_kl(1,1,1,1), &
            Mb_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array Bb_kl etc.'

    ! wavespeed, isotropic kernels
    allocate(rhobb_kl(1,1,1,1), &
            rhofbb_kl(1,1,1,1), &
            phib_kl(1,1,1,1), &
            ratio_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array rhobb_kl etc.'
    allocate(cs_kl(1,1,1,1),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array cs_kl'
    allocate(cpI_kl(1,1,1,1), &
            cpII_kl(1,1,1,1), stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array cpI_kl etc.'

  endif

  end subroutine read_mesh_databases_adjoint
