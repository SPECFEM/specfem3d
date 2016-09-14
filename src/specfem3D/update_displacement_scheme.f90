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


  subroutine update_displacement_scheme()

! explicit Newmark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par

  implicit none

  ! time marching

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displacement_acoustic()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displacement_elastic()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displacement_poroelastic()

  ! adjoint simulations: moho kernel
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    ispec2D_moho_top = 0
    ispec2D_moho_bot = 0
  endif

  end subroutine update_displacement_scheme

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_acoustic()

! updates acoustic potentials

  use specfem_par
  use specfem_par_acoustic
  use pml_par

  implicit none

  !local prameters
  integer :: ispec_cpml,ispec,i,j,k,iglob

  ! Newmark time marching

  if (.not. GPU_MODE) then
    ! wavefields on CPU
    ! updates (forward) acoustic potentials
    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
       do ispec_cpml=1,NSPEC_CPML
          ispec = CPML_to_spec(ispec_cpml)
          do i = 1, NGLLX
             do j = 1, NGLLY
                do k = 1, NGLLZ
                   iglob = ibool(i,j,k,ispec)
                   PML_potential_acoustic_old(i,j,k,ispec_cpml) = potential_acoustic(iglob) &
                       + deltatover2 * (1._CUSTOM_REAL - 2.0_CUSTOM_REAL*theta) * potential_dot_acoustic(iglob) &
                       + deltatsqover2 * (1._CUSTOM_REAL - theta) * potential_dot_dot_acoustic(iglob)
                enddo
             enddo
          enddo
       enddo
    endif
    potential_acoustic(:) = potential_acoustic(:) + &
                            deltat * potential_dot_acoustic(:) + &
                            deltatsqover2 * potential_dot_dot_acoustic(:)
    potential_dot_acoustic(:) = potential_dot_acoustic(:) + &
                                deltatover2 * potential_dot_dot_acoustic(:)
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      do ispec_cpml=1,NSPEC_CPML
         ispec = CPML_to_spec(ispec_cpml)
         do i = 1, NGLLX
            do j = 1, NGLLY
               do k = 1, NGLLZ
                  iglob = ibool(i,j,k,ispec)
                  PML_potential_acoustic_new(i,j,k,ispec_cpml) = potential_acoustic(iglob) &
                       + deltatover2 * (1._CUSTOM_REAL - 2.0_CUSTOM_REAL*theta) * potential_dot_acoustic(iglob)
               enddo
            enddo
         enddo
      enddo
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      ! updates acoustic backward/reconstructed fields
      if (PML_CONDITIONS) then
        if (nglob_interface_PML_acoustic > 0) then
          call read_potential_on_pml_interface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic, &
                                               nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)
        endif
      endif
      b_potential_acoustic(:) = b_potential_acoustic(:) + &
                                b_deltat * b_potential_dot_acoustic(:) + &
                                b_deltatsqover2 * b_potential_dot_dot_acoustic(:)
      b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + &
                                    b_deltatover2 * b_potential_dot_dot_acoustic(:)
      b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    endif

  else
    ! wavefields on GPU
    ! check
    if (SIMULATION_TYPE == 3) then
      if (PML_CONDITIONS) then
        call exit_MPI(myrank,'acoustic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
      endif
    endif

    ! updates acoustic potentials
    call it_update_displacement_ac_cuda(Mesh_pointer,deltat,deltatsqover2,deltatover2,b_deltat,b_deltatsqover2,b_deltatover2)
  endif ! GPU_MODE

  end subroutine update_displacement_acoustic


!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_elastic()

  ! updates elastic wavefields

  use specfem_par
  use specfem_par_elastic
  use pml_par

  implicit none

  ! local parameters
  integer :: ispec_cpml,ispec,i,j,k,iglob

  ! Newmark time marching

  if (.not. GPU_MODE) then
    ! wavefields on CPU

    ! updates elastic displacement and velocity
    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
        do ispec_cpml=1,NSPEC_CPML
           ispec = CPML_to_spec(ispec_cpml)
           do k = 1, NGLLZ
              do j = 1, NGLLY
                 do i = 1, NGLLX
                    iglob = ibool(i,j,k,ispec)
                    PML_displ_old(:,i,j,k,ispec_cpml)=displ(:,iglob) &
                        + deltatover2 * (1._CUSTOM_REAL - 2._CUSTOM_REAL*theta) * veloc(:,iglob) &
                        + deltatsqover2 * (1._CUSTOM_REAL - theta) * accel(:,iglob)

                 enddo
              enddo
           enddo
        enddo
    endif
    displ(:,:) = displ(:,:) + deltat * veloc(:,:) + deltatsqover2 * accel(:,:)
    veloc(:,:) = veloc(:,:) + deltatover2 * accel(:,:)
    if (SIMULATION_TYPE /= 1) accel_adj_coupling(:,:) = accel(:,:)
    accel(:,:) = 0._CUSTOM_REAL
    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
        do ispec_cpml=1,NSPEC_CPML
           ispec = CPML_to_spec(ispec_cpml)
           do k = 1, NGLLZ
              do j = 1, NGLLY
                 do i = 1, NGLLX
                    iglob = ibool(i,j,k,ispec)
                    PML_displ_new(:,i,j,k,ispec_cpml) = displ(:,iglob)&
                        + deltatover2 * (1._CUSTOM_REAL - 2.0_CUSTOM_REAL*theta) * veloc(:,iglob)

                 enddo
              enddo
           enddo
        enddo
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      ! elastic backward fields
      if (PML_CONDITIONS) then
        if (nglob_interface_PML_elastic > 0) then
          call read_field_on_pml_interface(b_accel,b_veloc,b_displ,nglob_interface_PML_elastic, &
                                           b_PML_field,b_reclen_PML_field)
        endif
      endif
      b_displ(:,:) = b_displ(:,:) + b_deltat * b_veloc(:,:) + b_deltatsqover2 * b_accel(:,:)
      b_veloc(:,:) = b_veloc(:,:) + b_deltatover2 * b_accel(:,:)
      b_accel(:,:) = 0._CUSTOM_REAL
    endif

  else
    ! wavefields on GPU

    ! check
    if (SIMULATION_TYPE == 3) then
      if (PML_CONDITIONS) then
        if (nglob_interface_PML_elastic > 0) then
          call exit_MPI(myrank,'elastic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
        endif
      endif
    endif

    ! updates elastic displacement and velocity
    ! Includes SIM_TYPE 1 & 3 (for noise tomography)
    call update_displacement_cuda(Mesh_pointer,deltat,deltatsqover2,deltatover2,b_deltat,b_deltatsqover2,b_deltatover2)
  endif ! GPU_MODE

  end subroutine update_displacement_elastic

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_poroelastic()

! updates poroelastic wavefields

  use specfem_par
  use specfem_par_poroelastic

  implicit none

  ! Newmark time marching

  if (.not. GPU_MODE) then
    ! wavefields on CPU

    ! updates poroelastic displacements and velocities
    ! solid phase
    displs_poroelastic(:,:) = displs_poroelastic(:,:) + deltat * velocs_poroelastic(:,:) + &
                              deltatsqover2 * accels_poroelastic(:,:)
    velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + deltatover2 * accels_poroelastic(:,:)
    accels_poroelastic(:,:) = 0._CUSTOM_REAL

    ! fluid phase
    displw_poroelastic(:,:) = displw_poroelastic(:,:) + deltat * velocw_poroelastic(:,:) + &
                              deltatsqover2 * accelw_poroelastic(:,:)
    velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + deltatover2 * accelw_poroelastic(:,:)
    accelw_poroelastic(:,:) = 0._CUSTOM_REAL

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) then
      ! poroelastic backward fields
      ! solid phase
      b_displs_poroelastic(:,:) = b_displs_poroelastic(:,:) + b_deltat * b_velocs_poroelastic(:,:) + &
                                  b_deltatsqover2 * b_accels_poroelastic(:,:)
      b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + b_deltatover2 * b_accels_poroelastic(:,:)
      b_accels_poroelastic(:,:) = 0._CUSTOM_REAL

      ! fluid phase
      b_displw_poroelastic(:,:) = b_displw_poroelastic(:,:) + b_deltat * b_velocw_poroelastic(:,:) + &
                                  b_deltatsqover2 * b_accelw_poroelastic(:,:)
      b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + b_deltatover2 * b_accelw_poroelastic(:,:)
      b_accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    endif

  else
    ! wavefields on GPU
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif ! GPU_MODE

  end subroutine update_displacement_poroelastic


