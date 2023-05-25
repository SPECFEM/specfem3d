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

  subroutine create_mass_matrices(nglob)

! returns precomputed mass matrix in rmass array

  use constants, only: CUSTOM_REAL

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION, &
    PML_CONDITIONS, STACEY_ABSORBING_CONDITIONS, DT

  ! global indices
  use generate_databases_par, only: nspec => NSPEC_AB, ibool

  use create_regions_mesh_ext_par

  ! PML
  use generate_databases_par, only: is_CPML,nspec_cpml,CPML_regions,CPML_to_spec, &
    d_store_x,d_store_y,d_store_z,K_store_x,K_store_y,K_store_z

  implicit none

  integer,intent(in) :: nglob

  ! local parameters
  integer :: ier

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! allocates memory
    allocate(rmass(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 660')
    if (ier /= 0) call exit_MPI_without_rank('error allocating array rmass')
    rmass(:) = 0._CUSTOM_REAL

    ! returns elastic mass matrix
    if (PML_CONDITIONS) then
      call define_mass_matrices_pml_elastic(nglob,nspec,nspec_irregular,DT,ibool,rhostore, &
                                            jacobianstore,irregular_element_number,jacobian_regular, &
                                            wxgll,wygll,wzgll,ispec_is_elastic, &
                                            nspec_cpml,is_CPML,CPML_regions,CPML_to_spec, &
                                            d_store_x,d_store_y,d_store_z, &
                                            K_store_x,K_store_y,K_store_z, &
                                            rmass)
    else
      call define_mass_matrices_elastic(nglob,nspec,nspec_irregular,ibool,rhostore, &
                                        jacobianstore,irregular_element_number,jacobian_regular, &
                                        wxgll,wygll,wzgll,ispec_is_elastic, &
                                        rmass)
    endif
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    ! allocates memory
    allocate(rmass_acoustic(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 661')
    if (ier /= 0) call exit_MPI_without_rank('error allocating array rmass_acoustic')
    rmass_acoustic(:) = 0._CUSTOM_REAL

    ! returns acoustic mass matrix
    if (PML_CONDITIONS) then
      call define_mass_matrices_pml_acoustic(nglob,nspec,nspec_irregular,DT,ibool,kappastore, &
                                             jacobianstore,irregular_element_number,jacobian_regular, &
                                             wxgll,wygll,wzgll,ispec_is_acoustic, &
                                             nspec_cpml,is_CPML,CPML_regions,CPML_to_spec, &
                                             d_store_x,d_store_y,d_store_z, &
                                             K_store_x,K_store_y,K_store_z, &
                                             rmass_acoustic)
    else
      call define_mass_matrices_acoustic(nglob,nspec,nspec_irregular,ibool,kappastore, &
                                         jacobianstore,irregular_element_number,jacobian_regular, &
                                         wxgll,wygll,wzgll,ispec_is_acoustic, &
                                         rmass_acoustic)
    endif
  endif

  ! poroelastic domains
  if (POROELASTIC_SIMULATION) then
    ! allocates memory
    allocate(rmass_solid_poroelastic(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 662')
    if (ier /= 0) call exit_MPI_without_rank('error in allocate rmass_solid_poroelastic')
    allocate(rmass_fluid_poroelastic(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 663')
    if (ier /= 0) call exit_MPI_without_rank('error in allocate rmass_fluid_poroelastic')
    rmass_solid_poroelastic(:) = 0._CUSTOM_REAL
    rmass_fluid_poroelastic(:) = 0._CUSTOM_REAL

    ! poroelastic mass matrices
    call define_mass_matrices_poroelastic(nglob,nspec,nspec_irregular,ibool,rhoarraystore,phistore,tortstore, &
                                          jacobianstore,irregular_element_number,jacobian_regular, &
                                          wxgll,wygll,wzgll,ispec_is_poroelastic, &
                                          rmass_solid_poroelastic,rmass_fluid_poroelastic)
  endif

  ! Stacey absorbing conditions (adds C*deltat/2 contribution to the mass matrices on Stacey edges)
  if (STACEY_ABSORBING_CONDITIONS) call create_mass_matrices_Stacey(nglob)

  ! ocean load mass matrix
  call create_mass_matrices_ocean_load(nglob)

  end subroutine create_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(nglob)

! returns precomputed mass matrix in rmass array

  use constants, only: myrank,CUSTOM_REAL,IMAIN

  use shared_parameters, only: APPROXIMATE_OCEAN_LOAD

  use generate_databases_par, only: NX_TOPO,NY_TOPO,itopo_bathy

  ! global indices
  use generate_databases_par, only: nspec => NSPEC_AB, ibool

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nglob

  ! local parameters
  integer :: ier

  ! creates ocean load mass matrix
  if (APPROXIMATE_OCEAN_LOAD) then

    if (myrank == 0) then
      write(IMAIN,*) '  ...creating ocean load mass matrix '
    endif

    ! adding ocean load mass matrix at ocean bottom
    NGLOB_OCEAN = nglob
    allocate(rmass_ocean_load(NGLOB_OCEAN),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 664')
    if (ier /= 0) stop 'error allocating array rmass_ocean_load'

    ! create ocean load mass matrix for degrees of freedom at ocean bottom
    rmass_ocean_load(:) = 0._CUSTOM_REAL

    call define_mass_matrices_ocean_load(nglob,nspec,ibool,xstore_unique,ystore_unique,zstore_unique, &
                                         num_free_surface_faces,free_surface_ispec,free_surface_ijk,free_surface_jacobian2Dw, &
                                         NX_TOPO,NY_TOPO,itopo_bathy, &
                                         ispec_is_elastic,rmass_ocean_load)

    ! adds regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

    ! allocate dummy array if no oceans
    NGLOB_OCEAN = 1
    allocate(rmass_ocean_load(NGLOB_OCEAN),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 665')
    if (ier /= 0) stop 'error allocating dummy array rmass_ocean_load'

  endif

  end subroutine create_mass_matrices_ocean_load

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_Stacey(nglob)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use constants, only: CUSTOM_REAL

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, DT

  ! global indices
  use generate_databases_par, only: nspec => NSPEC_AB, ibool

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nglob

  ! local parameters
  integer :: ier

  ! checks if anything to do
  if (num_abs_boundary_faces > 0) then
    nglob_xy = nglob
  else
    nglob_xy = 1
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    allocate( rmassx(nglob_xy), rmassy(nglob_xy), rmassz(nglob_xy), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 666')
    if (ier /= 0) stop 'error in allocate 21'
    rmassx(:) = 0._CUSTOM_REAL
    rmassy(:) = 0._CUSTOM_REAL
    rmassz(:) = 0._CUSTOM_REAL

    if (num_abs_boundary_faces > 0) then
      call define_mass_matrices_Stacey_elastic(nglob,nspec,DT,ibool,rho_vp,rho_vs, &
                                               num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                               abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                               ispec_is_elastic, &
                                               rmassx, rmassy, rmassz)
    endif
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    allocate( rmassz_acoustic(nglob_xy), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 667')
    if (ier /= 0) stop 'error in allocate 22'
    rmassz_acoustic(:) = 0._CUSTOM_REAL

    if (num_abs_boundary_faces > 0) then
      call define_mass_matrices_Stacey_acoustic(nglob,nspec,DT,ibool,rho_vp, &
                                                num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                                abs_boundary_jacobian2Dw, &
                                                ispec_is_acoustic, &
                                                rmassz_acoustic)
    endif
  endif

  end subroutine create_mass_matrices_Stacey

