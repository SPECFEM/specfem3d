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


  subroutine prepare_mass_matrices()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing mass matrices"
    call flush_IMAIN()
  endif

  ! synchronize all the processes before assembling the mass matrix
  ! to make sure all the nodes have finished to read their databases
  call synchronize_all()

  ! sets up mass matrices
  if (ELASTIC_SIMULATION) then
    ! switches to three-component mass matrix
    ! (rmassz was read in from database file)
    rmassx(:) = rmassz(:)
    rmassy(:) = rmassz(:)
  endif

  ! PML absorbing conditions (adds C*deltat/2 to the mass matrices in PML elements)
  if (PML_CONDITIONS) then
    call prepare_mass_matrices_PML()
  endif

  ! Stacey absorbing conditions (adds C*deltat/2 contribution to the mass matrices on Stacey edges)
  if (STACEY_ABSORBING_CONDITIONS) then
    call prepare_mass_matrices_Stacey()
  endif

  ! for exact rotation (adds rotational C*deltat/2 to the elastic mass matrix)
  if (ROTATION) then
    call prepare_mass_matrices_rotation()
  endif

  ! LTS mass matrices
  if (LTS_MODE) call lts_prepare_mass_matrices()

  ! the mass matrices need to be assembled with MPI here once and for all
  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_acoustic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_acoustic <= 0._CUSTOM_REAL) rmass_acoustic = 1._CUSTOM_REAL

    ! checks mass matrix
    if (minval(rmass_acoustic) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmass_acoustic')

    ! mass matrix inversion
    ! for efficiency, invert final mass matrix once and for all on each slice
    rmass_acoustic(:) = 1._CUSTOM_REAL / rmass_acoustic(:)
  endif

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! assemble mass matrix
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassx, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassy, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmassz, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmassx <= 0._CUSTOM_REAL) rmassx = 1._CUSTOM_REAL
    where(rmassy <= 0._CUSTOM_REAL) rmassy = 1._CUSTOM_REAL
    where(rmassz <= 0._CUSTOM_REAL) rmassz = 1._CUSTOM_REAL

    ! checks mass matrix
    if (minval(rmassx) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmassx')
    if (minval(rmassy) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmassy')
    if (minval(rmassz) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmassz')

    ! mass matrix inversion
    ! for efficiency, invert final mass matrix once and for all on each slice
    rmassx(:) = 1._CUSTOM_REAL / rmassx(:)
    rmassy(:) = 1._CUSTOM_REAL / rmassy(:)
    rmassz(:) = 1._CUSTOM_REAL / rmassz(:)

    ! ocean load
    if (APPROXIMATE_OCEAN_LOAD) then
      call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_ocean_load, &
                                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                        my_neighbors_ext_mesh)
      where(rmass_ocean_load <= 0._CUSTOM_REAL) rmass_ocean_load = 1._CUSTOM_REAL
      ! checks mass matrix
      if (minval(rmass_ocean_load) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmass_ocean_load')
      ! mass matrix inversion
      ! for efficiency, invert final mass matrix once and for all on each slice
      rmass_ocean_load(:) = 1._CUSTOM_REAL / rmass_ocean_load(:)
    endif
  endif

  ! poroelastic domains
  if (POROELASTIC_SIMULATION) then
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_solid_poroelastic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,rmass_fluid_poroelastic, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                      my_neighbors_ext_mesh)

    ! fills mass matrix with fictitious non-zero values to make sure it can be inverted globally
    where(rmass_solid_poroelastic <= 0._CUSTOM_REAL) rmass_solid_poroelastic = 1._CUSTOM_REAL
    where(rmass_fluid_poroelastic <= 0._CUSTOM_REAL) rmass_fluid_poroelastic = 1._CUSTOM_REAL

    ! checks mass matrix
    if (minval(rmass_solid_poroelastic) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmass_solid_poroelastic')
    if (minval(rmass_fluid_poroelastic) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for rmass_fluid_poroelastic')

    ! mass matrix inversion
    ! for efficiency, invert final mass matrix once and for all on each slice
    rmass_solid_poroelastic(:) = 1._CUSTOM_REAL / rmass_solid_poroelastic(:)
    rmass_fluid_poroelastic(:) = 1._CUSTOM_REAL / rmass_fluid_poroelastic(:)
  endif

  ! LTS mass matrices
  if (LTS_MODE) call lts_prepare_mass_matrices_invert()

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_mass_matrices_PML()

! modifies mass matrix in PML elements

  !use constants, only: IMAIN,myrank
  !use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, DT
  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic,rmassx,rmassy,rmassz
  use specfem_par_acoustic, only: ispec_is_acoustic,rmass_acoustic
  use pml_par

  implicit none

  ! checks if anything to do
  if (.not. PML_CONDITIONS) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  adding PML contributions'
  endif

  ! checks if anything to do in this slice
  if (NSPEC_CPML == 0) return

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    call add_mass_matrices_pml_elastic(NGLOB_AB,NSPEC_AB,nspec_irregular,DT,ibool,rhostore, &
                                          jacobianstore,irregular_element_number,jacobian_regular, &
                                          wxgll,wygll,wzgll,ispec_is_elastic, &
                                          NSPEC_CPML,is_CPML,CPML_regions,CPML_to_spec, &
                                          d_store_x,d_store_y,d_store_z, &
                                          K_store_x,K_store_y,K_store_z, &
                                          rmassx,rmassy,rmassz)
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call add_mass_matrices_pml_acoustic(NGLOB_AB,NSPEC_AB,nspec_irregular,DT,ibool,kappastore, &
                                           jacobianstore,irregular_element_number,jacobian_regular, &
                                           wxgll,wygll,wzgll,ispec_is_acoustic, &
                                           NSPEC_CPML,is_CPML,CPML_regions,CPML_to_spec, &
                                           d_store_x,d_store_y,d_store_z, &
                                           K_store_x,K_store_y,K_store_z, &
                                           rmass_acoustic)
  endif

  end subroutine prepare_mass_matrices_PML

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_mass_matrices_Stacey()

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix on Stacey edges;
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic,rmassx,rmassy,rmassz,rho_vp,rho_vs
  use specfem_par_acoustic, only: ispec_is_acoustic,rmass_acoustic

  implicit none

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return

  ! only for Newmark time schemes
  if (USE_LDDRK) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  adding Stacey contributions'
  endif

  ! checks if anything to do in this slice
  if (num_abs_boundary_faces == 0) return

  ! adds Stacey contributions to mass matrices
  ! elastic domains
  if (ELASTIC_SIMULATION) then
    call add_mass_matrices_Stacey_elastic(NGLOB_AB,NSPEC_AB,DT,ibool,rho_vp,rho_vs, &
                                             num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                             abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                             ispec_is_elastic, &
                                             rmassx,rmassy,rmassz)
  endif

  ! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    call add_mass_matrices_Stacey_acoustic(NGLOB_AB,NSPEC_AB,DT,ibool,rho_vp, &
                                              num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                              abs_boundary_jacobian2Dw, &
                                              ispec_is_acoustic, &
                                              rmass_acoustic)
  endif

  end subroutine prepare_mass_matrices_Stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_mass_matrices_rotation()

! in the case of rotation, add C*deltat/2 contribution to the mass matrix;
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use specfem_par
  use specfem_par_elastic, only: ispec_is_elastic,rmassx,rmassy,rmassz

  implicit none

  ! local parameters
  double precision :: weight,jacobianl
  double precision :: facx,facy,facz
  integer :: i,j,k,ispec,ispec_irreg,iglob

  ! checks if anything to do
  if (.not. ROTATION) return

  ! only applies to elastic domain
  if (.not. ELASTIC_SIMULATION) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  adding rotation contributions'
  endif

  ! C * dt/2 == [ 2 (omega x veloc) ] * dt/2 == (omega x veloc) * dt
  facx = (ROTATION_OMEGA(2) - ROTATION_OMEGA(3)) * DT
  facy = (ROTATION_OMEGA(3) - ROTATION_OMEGA(1)) * DT
  facz = (ROTATION_OMEGA(1) - ROTATION_OMEGA(2)) * DT

  ! adds rotation contributions
  do ispec = 1,NSPEC_AB
    ! elastic domain
    if (ispec_is_elastic(ispec)) then

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) jacobianl = jacobian_regular

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)

            weight = wxgll(i)*wygll(j)*wzgll(k)
            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            ! contribution: C*delta/2 == 2 * omega * 1/2 * DT * jacobian * weight
            rmassx(iglob) = rmassx(iglob) + real(facx * jacobianl * weight, kind=CUSTOM_REAL)
            rmassy(iglob) = rmassy(iglob) + real(facy * jacobianl * weight, kind=CUSTOM_REAL)
            rmassz(iglob) = rmassz(iglob) + real(facz * jacobianl * weight, kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    endif
  enddo

  end subroutine prepare_mass_matrices_rotation

