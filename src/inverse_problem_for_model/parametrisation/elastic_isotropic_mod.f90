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

module elastic_isotropic_mod

  implicit none

contains

  !================================================================================
  ! Parametrisation rho vp vs
  function rho_lambda_mu_to_rho_vp_vs(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)) / param_in(:,:,:,:,1))
    param_out(:,:,:,:,3) = sqrt(param_in(:,:,:,:,3) / param_in(:,:,:,:,1))
  end function rho_lambda_mu_to_rho_vp_vs

  function gradient_for_rho_vp_vs(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1) &
                        + grad_in(:,:,:,:,2) * (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2) &
                        + grad_in(:,:,:,:,3) * (param_in(:,:,:,:,3)**2)
    grad_out(:,:,:,:,2) = grad_in(:,:,:,:,2) * 2.* param_in(:,:,:,:,1) * param_in(:,:,:,:,2)
    grad_out(:,:,:,:,3) = grad_in(:,:,:,:,3) * 2.* param_in(:,:,:,:,1) * param_in(:,:,:,:,3) &
                        - grad_in(:,:,:,:,2) * 4.* param_in(:,:,:,:,1) * param_in(:,:,:,:,3)
  end function gradient_for_rho_vp_vs

  function rho_vp_vs_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,1) * (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2)
    param_out(:,:,:,:,3) = param_in(:,:,:,:,1) * param_in(:,:,:,:,3)**2
  end function rho_vp_vs_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Parametrisation rho ip is
  function rho_lambda_mu_to_rho_ip_is(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)   ! rho
    param_out(:,:,:,:,2) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)) * param_in(:,:,:,:,1))
    param_out(:,:,:,:,2) = sqrt(param_in(:,:,:,:,3) * param_in(:,:,:,:,1))
  end function rho_lambda_mu_to_rho_ip_is

  function gradient_for_rho_ip_is(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1) &
         - grad_in(:,:,:,:,2) * (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2) / param_in(:,:,:,:,1)**2&
         - grad_in(:,:,:,:,3) *  param_in(:,:,:,:,3)**2 / param_in(:,:,:,:,1)**2
    grad_out(:,:,:,:,2) = grad_in(:,:,:,:,2) * 2.* param_in(:,:,:,:,2) / param_in(:,:,:,:,1)
    grad_out(:,:,:,:,3) = grad_in(:,:,:,:,3) * 2.* param_in(:,:,:,:,3) / param_in(:,:,:,:,1) &
                        - grad_in(:,:,:,:,2) * 4.* param_in(:,:,:,:,3) / param_in(:,:,:,:,1)
  end function gradient_for_rho_ip_is

  function rho_ip_is_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2) / param_in(:,:,:,:,1)
    param_out(:,:,:,:,3) = param_in(:,:,:,:,3)**2 / param_in(:,:,:,:,1)
  end function rho_ip_is_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rho, Vb, Vs
  function rho_lambda_mu_to_rho_vb_vs(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)/3.) / param_in(:,:,:,:,1))
    param_out(:,:,:,:,3) = sqrt(param_in(:,:,:,:,3) / param_in(:,:,:,:,1))
  end function rho_lambda_mu_to_rho_vb_vs

  function gradient_for_rho_vb_vs(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1) &
                        + grad_in(:,:,:,:,2) * (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2 /3.) &
                        + grad_in(:,:,:,:,3) * (param_in(:,:,:,:,3)**2)
    grad_out(:,:,:,:,2) = grad_in(:,:,:,:,2) * 2.* param_in(:,:,:,:,1) * param_in(:,:,:,:,2)
    grad_out(:,:,:,:,3) = grad_in(:,:,:,:,3) * 2.* param_in(:,:,:,:,1) * param_in(:,:,:,:,3) &
                        - grad_in(:,:,:,:,2) * 4.* param_in(:,:,:,:,1) * param_in(:,:,:,:,3) / 3.
  end function gradient_for_rho_vb_vs

  function rho_vb_vs_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,1) * (param_in(:,:,:,:,2)**2 - 2.*param_in(:,:,:,:,3)**2 /3.)
    param_out(:,:,:,:,3) = param_in(:,:,:,:,1) * param_in(:,:,:,:,3)**2
  end function rho_vb_vs_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rho, kappa, lambda
  function rho_lambda_mu_to_rho_kappa_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3))/ 3.
    param_out(:,:,:,:,3) = param_in(:,:,:,:,3)
  end function rho_lambda_mu_to_rho_kappa_mu

  function gradient_for_rho_kappa_mu(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1)
    grad_out(:,:,:,:,2) = grad_in(:,:,:,:,2)  +  3. * grad_in(:,:,:,:,3) / 2.
    grad_out(:,:,:,:,3) = grad_in(:,:,:,:,3)  -  2. * grad_in(:,:,:,:,2) / 3.
  end function gradient_for_rho_vp_vs

  function rho_kappa_mu_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,2) - 2.*param_in(:,:,:,:,3))/ 3.
    param_out(:,:,:,:,3) = param_in(:,:,:,:,3)
  end function rho_kappa_mu_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rho, vp, vp/vs
  function rho_lambda_mu_to_rho_vp_vpovervs(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)) / param_in(:,:,:,:,1))
    param_out(:,:,:,:,3) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)) / param_in(:,:,:,:,3))
  end function rho_lambda_mu_to_rho_vp_vpovervs

  function gradient_for_rho_vp_vpovervs(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1) &
         + grad_in(:,:,:,:,2)*(param_in(:,:,:,:,2)**2 * (param_in(:,:,:,:,3)**2 -2.)) / param_in(:,:,:,:,3)**2 &
         + grad_in(:,:,:,:,3) * param_in(:,:,:,:,2)**2  / param_in(:,:,:,:,3)**2
    grad_out(:,:,:,:,2) = 2. * grad_in(:,:,:,:,3) * param_in(:,:,:,:,1) * param_in(:,:,:,:,2) / param_in(:,:,:,:,3)**2 &
         + 2.*param_in(:,:,:,:,1) * param_in(:,:,:,:,2)*(param_in(:,:,:,:,3)**2 -2.) * grad_in(:,:,:,:,2) / param_in(:,:,:,:,3)**2
    grad_out(:,:,:,:,3) = -2. * grad_in(:,:,:,:,3) * param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 / param_in(:,:,:,:,3)**3 &
         + 4.*param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 * grad_in(:,:,:,:,2) / param_in(:,:,:,:,3)**3
    !!! WARNING THIS ONE DESERVES TO BE RECHECKED !!!
  end function gradient_for_rho_vp_vpovervs

  function rho_vp_vpovervs_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 * (param_in(:,:,:,:,3)**2 -2.) / param_in(:,:,:,:,3)**2
    param_out(:,:,:,:,3) = param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 / param_in(:,:,:,:,3)**2
  end function rho_vp_vpovervs_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Rho, vs, vp/vs
   function rho_lambda_mu_to_rho_vpovervs_vs(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = sqrt((param_in(:,:,:,:,2) + 2.*param_in(:,:,:,:,3)) / param_in(:,:,:,:,3))
    param_out(:,:,:,:,3) = sqrt(param_in(:,:,:,:,3) / param_in(:,:,:,:,1))
  end function rho_lambda_mu_to_rho_vpovervs_vs

  function gradient_for_rho_vpovervs_vs(param_in,grad_in,grad_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: grad_in, param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: grad_out
    grad_out(:,:,:,:,1) = grad_in(:,:,:,:,1) &
         + grad_in(:,:,:,:,2) * (param_in(:,:,:,:,3)**2 * (param_in(:,:,:,:,2)**2 -2.)) &
         + grad_in(:,:,:,:,3) * param_in(:,:,:,:,3)**2
    grad_out(:,:,:,:,2) = 2. * grad_in(:,:,:,:,3) * param_in(:,:,:,:,1) * param_in(:,:,:,:,3) &
         + 2.*param_in(:,:,:,:,1) * param_in(:,:,:,:,3)*(param_in(:,:,:,:,3)**2 -2.) * grad_in(:,:,:,:,2)
    grad_out(:,:,:,:,3) = 2.*param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 * grad_in(:,:,:,:,2) * param_in(:,:,:,:,3)
    !!! WARNING THIS ONE DESERVES TO BE RECHECKED !!!
  end function gradient_for_rho_vpovervs_vs

  function rho_vpovervs_vs_to_rho_lambda_mu(param_in,param_out)
    real, dimension(:,:,:,:,:), allocatable, intent(in)  :: param_in
    real, dimension(:,:,:,:,:), allocatable, intent(out) :: param_out
    param_out(:,:,:,:,1) = param_in(:,:,:,:,1)
    param_out(:,:,:,:,2) = param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2 * (param_in(:,:,:,:,3)**2 -2.)
    param_out(:,:,:,:,3) = param_in(:,:,:,:,1) * param_in(:,:,:,:,2)**2
  end function rho_vpovervs_vs_to_rho_lambda_mu
  !--------------------------------------------------------------------------------

end module elastic_isotropic_mod
