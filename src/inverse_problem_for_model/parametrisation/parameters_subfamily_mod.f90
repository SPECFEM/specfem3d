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

module parameters_subfamily_mod

  implicit none

contains

  !================================================================================
  ! Chain rules declarations for log(par)
  function parameter_to_logarithm_of_parameter(param_in) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,:) = log(param_in(:,:,:,:,:))
  end function parameter_to_logarithm_of_parameter

  function gradient_for_logarithm_of_parameter(grad_in,param_in) result(grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,:) = grad_in(:,:,:,:,:) * param_in(:,:,:,:,:)
  end function gradient_of_logarithm_of_parameter

  function logarithm_of_parameter_to_parameter(param_in) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,:) = exp(param_in(:,:,:,:,:))
  end function logarithm_of_parameter_to_parameter
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Chain rules declarations for par/par0
  function parameter_to_adimensional_parameter(param_in, param_ref) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    integer :: npar, ipar
    npar = size(param_in,1)
    do ipar = 1, npar
       param_out(ipar,:,:,:,:) = param_in(ipar,:,:,:,:) / param_ref(ipar)
    enddo
  end function parameter_to_adimensional_parameter

  function gradient_for_adimentional_parameter(grad_in, param_ref) result(grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    integer :: npar, ipar
    npar = size(param_in,1)
    do ipar = 1, npar
       grad_out(ipar,:,:,:,:) = grad_in(ipar,:,:,:,:) * param_ref(idim)
    enddo
  end function gradient_for_adimentional_parameter

  function adimentional_parameter_to_parameter(param_in, param_ref) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    npar = size(param_in,1)
    do ipar = 1, npar
       param_out(:,:,:,:,:) = param_in(:,:,:,:,:) / param_ref(idim)
    enddo
  end function adimentional_parameter_to_parameter
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Chain rule declaration for log(par/par0)r
  function parameter_to_logarithm_of_adimensional_parameter(param_in, param_ref) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    integer :: npar, ipar
    npar = size(param_in,1)
    do ipar = 1, npar
       param_out(ipar,:,:,:,:) = log(param_in(ipar,:,:,:,:) / param_ref(ipar))
    enddo
  end function parameter_to_logarithm_of_adimensional_parameter

  function gradient_for_logarithm_of_adimensional_parameter(param_in, grad_in, param_ref) result(grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in, grad_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    integer :: npar, ipar
    npar = size(param_in,1)
    do ipar = 1, npar
       param_out(ipar,:,:,:,:) = param_in(ipar,:,:,:,:) * param_ref(ipar) * grad_in(ipar,:,:,:,:)
    enddo
  end function gradient_for_logarithm_of_adimensional_parameter

  function logarithm_of_adimentional_parameter_to_parameter(param_in, param_ref) result(param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:),         allocatable, intent(in)  :: param_ref
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    integer :: npar, ipar
    npar = size(param_in,1)
    do ipar = 1, npar
       param_out(ipar,:,:,:,:) = exp(param_in(ipar,:,:,:,:)) * param_ref(ipar)
    enddo
  end function logarithm_of_adimentional_parameter_to_parameter
  !--------------------------------------------------------------------------------

end module parameters_subfamily_mod
