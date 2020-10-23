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


  subroutine update_displ_Newmark()

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
  !
  ! only updates forward fields
  !
  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displacement_acoustic()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displacement_elastic()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displacement_poroelastic()

  end subroutine update_displ_Newmark

!
!--------------------------------------------------------------------------------------------------------------
!


  subroutine update_displ_Newmark_backward()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! time marching
  !
  ! updates only reconstructed/backward fields

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displacement_acoustic_backward()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displacement_elastic_backward()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displacement_poroelastic_backward()

  end subroutine update_displ_Newmark_backward


!--------------------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!--------------------------------------------------------------------------------------------------------------

  subroutine update_displacement_acoustic()

! updates potentials (forward fields only)

  use specfem_par
  use specfem_par_acoustic
  use pml_par

  implicit none

  if (.not. GPU_MODE) then
    ! wavefields on CPU
    ! PML store old field
    if (PML_CONDITIONS) call update_displ_acoustic_PML(PML_potential_acoustic_old, &
                                                       potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                                       deltatover2,deltatsqover2)

    ! Newmark predictor step for forward potentials
    call update_displ_acoustic(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                               deltat,deltatover2,deltatsqover2)

    ! PML store new updated field
    if (PML_CONDITIONS) call update_displ_acoustic_PML(PML_potential_acoustic_new, &
                                                       potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                                       deltatover2,deltatsqover2)
  else
    ! wavefields on GPU
    ! updates acoustic potentials
    call update_displacement_ac_cuda(Mesh_pointer,deltat,deltatover2,deltatsqover2,1) ! 1 == forward fields
  endif

  end subroutine update_displacement_acoustic

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_acoustic_backward()

! updates acoustic potentials (backward fields only)

  use specfem_par
  use specfem_par_acoustic
  use pml_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! Newmark time marching
  if (.not. GPU_MODE) then
    ! wavefields on CPU
    ! adjoint simulations
    ! PML updates acoustic backward/reconstructed fields on boundary
    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      ! reconstructs backpropagated forward wavefield on boundary for kernel simulations with PML
      if (nglob_interface_PML_acoustic > 0) then
        call read_potential_on_pml_interface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic, &
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)
      endif
    endif

    ! Newmark predictor step for backward potentials
    call update_displ_acoustic(b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                               b_deltat,b_deltatover2,b_deltatsqover2)
  else
    ! wavefields on GPU
    ! check
    if (PML_CONDITIONS) then
      call exit_MPI(myrank,'acoustic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
    endif

    ! updates acoustic backward potentials
    call update_displacement_ac_cuda(Mesh_pointer,b_deltat,b_deltatover2,b_deltatsqover2,3) ! 3 == backward fields
  endif ! GPU_MODE

  end subroutine update_displacement_acoustic_backward


!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displ_acoustic(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                   deltat,deltatover2,deltatsqover2)

! updates acoustic potentials

  use specfem_par, only: CUSTOM_REAL,NGLOB_AB

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_acoustic, &
                                                               potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), intent(in) :: deltat,deltatover2,deltatsqover2

  ! Newmark time marching
  potential_acoustic(:) = potential_acoustic(:) + &
                          deltat * potential_dot_acoustic(:) + &
                          deltatsqover2 * potential_dot_dot_acoustic(:)
  potential_dot_acoustic(:) = potential_dot_acoustic(:) + &
                              deltatover2 * potential_dot_dot_acoustic(:)
  potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL

  end subroutine update_displ_acoustic


!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displ_acoustic_PML(PML_potential_acoustic, &
                                       potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                       deltatover2,deltatsqover2)

! updates acoustic potentials in PML

  use specfem_par, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,ibool,NGLOB_AB,PML_CONDITIONS
  !use specfem_par_acoustic
  use pml_par, only: NSPEC_CPML,CPML_to_spec,theta

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),intent(out) :: PML_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), intent(in) :: deltatover2,deltatsqover2

  !local prameters
  integer :: ispec_cpml,ispec,i,j,k,iglob

  ! checks if anything to do
  if (.not. PML_CONDITIONS) return
  if (NSPEC_CPML == 0) return

  ! updates (forward) acoustic potentials
  do ispec_cpml = 1,NSPEC_CPML
    ispec = CPML_to_spec(ispec_cpml)
    do i = 1, NGLLX
       do j = 1, NGLLY
          do k = 1, NGLLZ
             iglob = ibool(i,j,k,ispec)
             PML_potential_acoustic(i,j,k,ispec_cpml) = potential_acoustic(iglob) &
                                  + deltatover2 * (1._CUSTOM_REAL - 2.0_CUSTOM_REAL*theta) * potential_dot_acoustic(iglob) &
                                  + deltatsqover2 * (1._CUSTOM_REAL - theta) * potential_dot_dot_acoustic(iglob)
          enddo
       enddo
    enddo
  enddo

  end subroutine update_displ_acoustic_PML


!--------------------------------------------------------------------------------------------------------------
!
! elastic domains
!
!--------------------------------------------------------------------------------------------------------------

  subroutine update_displacement_elastic()

! updates elastic wavefields (forward fields only)

  use specfem_par
  use specfem_par_elastic
  use pml_par

  implicit none

  ! Newmark time marching
  if (.not. GPU_MODE) then
    ! wavefields on CPU

    ! PML stores old wavefield
    if (PML_CONDITIONS) call update_displ_elastic_PML(PML_displ_old,displ,veloc,accel,deltatover2,deltatsqover2)

    ! coupling with adjoint wavefields
    if (SIMULATION_TYPE /= 1) accel_adj_coupling(:,:) = accel(:,:)

    call update_displ_elastic(displ,veloc,accel,deltat,deltatover2,deltatsqover2)

    ! PML new wavefield
    if (PML_CONDITIONS) call update_displ_elastic_PML(PML_displ_new,displ,veloc,accel,deltatover2,deltatsqover2)
  else
    ! wavefields on GPU
    ! updates elastic displacement and velocity
    call update_displacement_cuda(Mesh_pointer,deltat,deltatover2,deltatsqover2,1) ! 1 == forward
  endif ! GPU_MODE

  end subroutine update_displacement_elastic

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_elastic_backward()

! updates elastic wavefields (backward fields only)

  use specfem_par
  use specfem_par_elastic
  use pml_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! Newmark time marching
  if (.not. GPU_MODE) then
    ! wavefields on CPU

    ! updates elastic displacement and velocity
    ! elastic backward fields
    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      if (nglob_interface_PML_elastic > 0) then
        call read_field_on_pml_interface(b_accel,b_veloc,b_displ,nglob_interface_PML_elastic, &
                                         b_PML_field,b_reclen_PML_field)
      endif
    endif

    ! adjoint simulations
    call update_displ_elastic(b_displ,b_veloc,b_accel,b_deltat,b_deltatover2,b_deltatsqover2)
  else
    ! wavefields on GPU
    if (PML_CONDITIONS) then
      if (nglob_interface_PML_elastic > 0) then
        call exit_MPI(myrank,'elastic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
      endif
    endif

    ! updates elastic backward displacement and velocity
    call update_displacement_cuda(Mesh_pointer,b_deltat,b_deltatover2,b_deltatsqover2,3) ! 3 == backward
  endif ! GPU_MODE

  end subroutine update_displacement_elastic_backward



!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displ_elastic(displ,veloc,accel,deltat,deltatover2,deltatsqover2)

! updates elastic wavefields in PML

  use specfem_par, only: NDIM,CUSTOM_REAL,NGLOB_AB,FORCE_VECTORIZATION_VAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), intent(in) :: deltat,deltatover2,deltatsqover2

  ! local parameters
  integer :: i

  ! Newmark time marching
  ! updates elastic displacement and velocity
  if (FORCE_VECTORIZATION_VAL) then
! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( NGLOB_AB, displ, veloc, accel, &
!$OMP deltat, deltatsqover2, deltatover2 ) &
!$OMP PRIVATE(i)
!$OMP DO
    do i = 1,NDIM * NGLOB_AB
      displ(i,1) = displ(i,1) + deltat * veloc(i,1) + deltatsqover2 * accel(i,1)
      veloc(i,1) = veloc(i,1) + deltatover2 * accel(i,1)
      accel(i,1) = 0._CUSTOM_REAL
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  else
    do i = 1,NGLOB_AB
      displ(:,i) = displ(:,i) + deltat * veloc(:,i) + deltatsqover2 * accel(:,i)
      veloc(:,i) = veloc(:,i) + deltatover2 * accel(:,i)
      accel(:,i) = 0._CUSTOM_REAL
    enddo
  endif

  end subroutine update_displ_elastic


!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displ_elastic_PML(PML_displ,displ,veloc,accel,deltatover2,deltatsqover2)

! updates elastic wavefields in PML

  use specfem_par, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,ibool,NGLOB_AB,PML_CONDITIONS
  !use specfem_par_acoustic
  use pml_par, only: NSPEC_CPML,CPML_to_spec,theta

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),intent(out) :: PML_displ
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), intent(in) :: deltatover2,deltatsqover2

  !local prameters
  integer :: ispec_cpml,ispec,i,j,k,iglob

  ! checks if anything to do
  if (.not. PML_CONDITIONS) return
  if (NSPEC_CPML == 0) return

  ! Newmark time marching
  ! updates elastic displacement and velocity
  do ispec_cpml = 1,NSPEC_CPML
    ispec = CPML_to_spec(ispec_cpml)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          PML_displ(:,i,j,k,ispec_cpml) = displ(:,iglob) &
                                        + deltatover2 * (1._CUSTOM_REAL - 2._CUSTOM_REAL*theta) * veloc(:,iglob) &
                                        + deltatsqover2 * (1._CUSTOM_REAL - theta) * accel(:,iglob)
        enddo
      enddo
    enddo
  enddo

  end subroutine update_displ_elastic_PML

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_poroelastic()

! updates poroelastic wavefields (forward fields only)

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

  else
    ! wavefields on GPU
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif ! GPU_MODE

  end subroutine update_displacement_poroelastic

!
!--------------------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_poroelastic_backward()

! updates poroelastic wavefields (backward fields only)

  use specfem_par
  use specfem_par_poroelastic

  implicit none

  ! Newmark time marching
  if (.not. GPU_MODE) then
    ! wavefields on CPU

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

  else
    ! wavefields on GPU
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif ! GPU_MODE

  end subroutine update_displacement_poroelastic_backward

