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

  subroutine prepare_wavefields()

! initializes arrays

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! initialize acoustic arrays to zero
  if (ACOUSTIC_SIMULATION) then
    potential_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_acoustic(:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) potential_acoustic(:) = VERYSMALLVAL
  endif

  ! initialize elastic arrays to zero/verysmallvall
  if (ELASTIC_SIMULATION) then
    displ(:,:) = 0._CUSTOM_REAL
    veloc(:,:) = 0._CUSTOM_REAL
    accel(:,:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) displ(:,:) = VERYSMALLVAL
  endif

  ! initialize poroelastic arrays to zero/verysmallvall
  if (POROELASTIC_SIMULATION) then
    displs_poroelastic(:,:) = 0._CUSTOM_REAL
    velocs_poroelastic(:,:) = 0._CUSTOM_REAL
    accels_poroelastic(:,:) = 0._CUSTOM_REAL
    displw_poroelastic(:,:) = 0._CUSTOM_REAL
    velocw_poroelastic(:,:) = 0._CUSTOM_REAL
    accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    ! put negligible initial value to avoid very slow underflow trapping
    if (FIX_UNDERFLOW_PROBLEM) displs_poroelastic(:,:) = VERYSMALLVAL
    if (FIX_UNDERFLOW_PROBLEM) displw_poroelastic(:,:) = VERYSMALLVAL
  endif

  end subroutine prepare_wavefields

