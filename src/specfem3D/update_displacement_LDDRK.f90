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


  subroutine update_displ_lddrk()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic

  implicit none

  if (ACOUSTIC_SIMULATION) then
    potential_dot_dot_acoustic = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) potential_dot_dot_acoustic = VERYSMALLVAL
  endif

  if (ELASTIC_SIMULATION) then
    accel = 0._CUSTOM_REAL
    if (FIX_UNDERFLOW_PROBLEM) accel = VERYSMALLVAL
  endif

  end subroutine update_displ_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_displ_lddrk_backward()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic

  end subroutine update_displ_lddrk_backward

!
!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_potential_dot_acoustic_lddrk()

! updates acceleration, velocity and displacement in acoustic region (outer core)

  use specfem_par
  use specfem_par_acoustic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! forward wavefields
  call update_acoustic_lddrk(NGLOB_AB,NGLOB_AB_LDDRK, &
                             potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                             potential_acoustic_lddrk,potential_dot_acoustic_lddrk, &
                             deltat,alpha,beta)

  end subroutine update_potential_dot_acoustic_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_acoustic_lddrk(NGLOB,NGLOB_LDDRK, &
                                   potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                   potential_acoustic_lddrk,potential_dot_acoustic_lddrk, &
                                   deltat,alpha,beta)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: NGLOB,NGLOB_LDDRK

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(inout) :: &
            potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_LDDRK),intent(inout) :: &
            potential_acoustic_lddrk,potential_dot_acoustic_lddrk

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

  if (NGLOB /= NGLOB_LDDRK) stop 'error in definition of NGLOB_LDDRK'

  ! low-memory Runge-Kutta scheme

  do i = 1,NGLOB
    ! low-memory Runge-Kutta: intermediate storage wavefields
    potential_dot_acoustic_lddrk(i) = alpha * potential_dot_acoustic_lddrk(i) + deltat * potential_dot_dot_acoustic(i)
    potential_acoustic_lddrk(i) = alpha * potential_acoustic_lddrk(i) + deltat * potential_dot_acoustic(i)
    ! updates wavefields
    potential_dot_acoustic(i) = potential_dot_acoustic(i) + beta * potential_dot_acoustic_lddrk(i)
    potential_acoustic(i) = potential_acoustic(i) + beta * potential_acoustic_lddrk(i)
  enddo

  end subroutine update_acoustic_lddrk

!
!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_lddrk()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use specfem_par
  use specfem_par_elastic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! forward wavefields
  call update_elastic_lddrk(NGLOB_AB,NGLOB_AB_LDDRK,displ,veloc,accel, &
                            displ_lddrk,veloc_lddrk,deltat,alpha,beta)

  end subroutine update_veloc_elastic_lddrk

!
!-------------------------------------------------------------------------------------------------
!

! not used yet...
!  subroutine update_veloc_elastic_lddrk_backward()
!
!! updates acceleration,velocity and displacement in elastic regions
!
!  implicit none
!
!  end subroutine update_veloc_elastic_lddrk_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_elastic_lddrk(NGLOB,NGLOB_LDDRK,displ,veloc,accel, &
                                  displ_lddrk,veloc_lddrk,deltat,alpha,beta)

  use constants, only: CUSTOM_REAL,NDIM

  implicit none

  integer,intent(in) :: NGLOB,NGLOB_LDDRK

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_LDDRK),intent(inout) :: displ_lddrk,veloc_lddrk

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

  if (NGLOB /= NGLOB_LDDRK) stop 'error in definition of NGLOB_LDDRK'

  ! low-memory Runge-Kutta scheme

  do i = 1,NGLOB
    ! low-memory Runge-Kutta: intermediate storage wavefields
    veloc_lddrk(:,i) = alpha * veloc_lddrk(:,i) + deltat * accel(:,i)
    displ_lddrk(:,i) = alpha * displ_lddrk(:,i) + deltat * veloc(:,i)
    ! updates wavefields
    veloc(:,i) = veloc(:,i) + beta * veloc_lddrk(:,i)
    displ(:,i) = displ(:,i) + beta * displ_lddrk(:,i)
  enddo

  end subroutine update_elastic_lddrk


