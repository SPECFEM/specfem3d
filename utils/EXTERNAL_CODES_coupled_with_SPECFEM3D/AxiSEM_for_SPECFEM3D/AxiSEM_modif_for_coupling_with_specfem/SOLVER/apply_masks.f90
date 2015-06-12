!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!=====================
  module apply_masks 
!=====================

  use global_parameters, only: zero, realkind

  implicit none

  public :: apply_axis_mask_onecomp, apply_axis_mask_twocomp
  public :: apply_axis_mask_threecomp, apply_axis_mask_scal
  private

contains

  ! These routine applies a mask by retaining in the array those components which do not 
  ! belong to the axis of rotation. It sets to zero the axial components of the array, 
  ! for the non-axisymmetric components of the variables have to vanish on the axis of 
  ! rotation
!-----------------------------------------------------------------------------------------
pure subroutine apply_axis_mask_scal(u, nel, ax_array, nax_array)
  ! for a scalar array

  integer, intent(in)               :: nel, nax_array
  real(kind=realkind),intent(inout) :: u(0:,0:,:)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,:,ax_array(ielem)) = zero
  end do

end subroutine apply_axis_mask_scal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine apply_axis_mask_onecomp(u, nel, ax_array, nax_array)
  ! for the first component of the array

  integer, intent(in)               :: nel,nax_array
  real(kind=realkind),intent(inout) :: u(0:,0:,:,:)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,:,ax_array(ielem),1) = zero
  end do

end subroutine apply_axis_mask_onecomp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine apply_axis_mask_twocomp(u, nel, ax_array, nax_array)
  ! for the 2nd and 3rd component of the array
  
  integer, intent(in)               :: nel, nax_array
  real(kind=realkind),intent(inout) :: u(0:,0:,:,:)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,:,ax_array(ielem),2:3) = zero
  end do

end subroutine apply_axis_mask_twocomp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine apply_axis_mask_threecomp(u, nel, ax_array, nax_array)
  ! for the all components of the array

  integer, intent(in)               :: nel, nax_array
  real(kind=realkind),intent(inout) :: u(0:,0:,:,:)
  integer, intent(in)               :: ax_array(nax_array)
  integer                           :: ielem

  do ielem = 1, nax_array
     u(0,:,ax_array(ielem),1:3) = zero
  end do

end subroutine apply_axis_mask_threecomp

!========================
  end module apply_masks 
!========================
