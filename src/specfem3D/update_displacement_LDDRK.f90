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

  subroutine update_displ_lddrk()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic

  implicit none

  if (ACOUSTIC_SIMULATION) then
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
  call update_elastic_lddrk(NGLOB_AB,displ,veloc,accel,deltat,&
                            displ_lddrk,veloc_lddrk,alpha,beta)

  end subroutine update_veloc_elastic_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_lddrk_backward()

! updates acceleration,velocity and displacement in elastic regions

  end subroutine update_veloc_elastic_lddrk_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_elastic_lddrk(NGLOB,displ,veloc,accel,deltat, &
                                  displ_lddrk,veloc_lddrk,alpha,beta)

  use constants,only: CUSTOM_REAL,NDIM

  implicit none

  integer,intent(in) :: NGLOB

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ_lddrk,veloc_lddrk

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

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


