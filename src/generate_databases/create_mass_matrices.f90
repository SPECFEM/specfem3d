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

! returns precomputed mass matrix in rmass arrays

  use constants, only: CUSTOM_REAL,IMAIN,myrank

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION

  ! global indices
  use generate_databases_par, only: nspec => NSPEC_AB, ibool

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nglob

  ! local parameters
  integer :: ier

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! user info
    if (myrank == 0) then
      write(IMAIN,*) '     elastic mass matrix'
    endif

    ! allocates memory
    allocate(rmass_elastic(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 660')
    if (ier /= 0) call exit_MPI_without_rank('error allocating array rmass_elastic')
    rmass_elastic(:) = 0._CUSTOM_REAL

    ! defines mass matrix on all elastic elements
    call define_mass_matrices_elastic(nglob,nspec,nspec_irregular,ibool,rhostore, &
                                      jacobianstore,irregular_element_number,jacobian_regular, &
                                      wxgll,wygll,wzgll,ispec_is_elastic, &
                                      rmass_elastic)
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    ! user info
    if (myrank == 0) then
      write(IMAIN,*) '     acoustic mass matrix'
    endif

    ! allocates memory
    allocate(rmass_acoustic(nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 661')
    if (ier /= 0) call exit_MPI_without_rank('error allocating array rmass_acoustic')
    rmass_acoustic(:) = 0._CUSTOM_REAL

    ! defines mass matrix on all acoustic elements
    call define_mass_matrices_acoustic(nglob,nspec,nspec_irregular,ibool,kappastore, &
                                       jacobianstore,irregular_element_number,jacobian_regular, &
                                       wxgll,wygll,wzgll,ispec_is_acoustic, &
                                       rmass_acoustic)
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

    ! user info
    if (myrank == 0) then
      write(IMAIN,*) '     poroelastic mass matrix'
    endif

    ! poroelastic mass matrices
    call define_mass_matrices_poroelastic(nglob,nspec,nspec_irregular,ibool,rhoarraystore,phistore,tortstore, &
                                          jacobianstore,irregular_element_number,jacobian_regular, &
                                          wxgll,wygll,wzgll,ispec_is_poroelastic, &
                                          rmass_solid_poroelastic,rmass_fluid_poroelastic)
  endif

  ! ocean load mass matrix
  call create_mass_matrices_ocean_load(nglob)

  end subroutine create_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(nglob)

! compute mass matrix contribution in rmass_ocean_load array

  use constants, only: myrank,CUSTOM_REAL,IMAIN

  use shared_parameters, only: ELASTIC_SIMULATION,APPROXIMATE_OCEAN_LOAD

  use generate_databases_par, only: NX_TOPO,NY_TOPO,itopo_bathy

  ! global indices
  use generate_databases_par, only: nspec => NSPEC_AB, ibool

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nglob

  ! local parameters
  integer :: ier

  ! creates ocean load mass matrix (only for elastic domains)
  if (APPROXIMATE_OCEAN_LOAD .and. ELASTIC_SIMULATION) then
    ! user info
    if (myrank == 0) then
      write(IMAIN,*) '     ocean load mass matrix '
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
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass_elastic(:)
  else
    ! allocate dummy array if no oceans
    NGLOB_OCEAN = 1
    allocate(rmass_ocean_load(NGLOB_OCEAN),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 665')
    if (ier /= 0) stop 'error allocating dummy array rmass_ocean_load'
  endif

  end subroutine create_mass_matrices_ocean_load
