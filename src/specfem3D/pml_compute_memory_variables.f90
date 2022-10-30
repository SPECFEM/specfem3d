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


subroutine pml_compute_memory_variables_elastic(ispec,ispec_CPML, &
                                                tempx1,tempy1,tempz1, &
                                                tempx2,tempy2,tempz2, &
                                                tempx3,tempy3,tempz3, &
                                                PML_dux_dxl, PML_dux_dyl, PML_dux_dzl, &
                                                PML_duy_dxl, PML_duy_dyl, PML_duy_dzl, &
                                                PML_duz_dxl, PML_duz_dyl, PML_duz_dzl, &
                                                PML_dux_dxl_old, PML_dux_dyl_old, PML_dux_dzl_old, &
                                                PML_duy_dxl_old, PML_duy_dyl_old, PML_duy_dzl_old, &
                                                PML_duz_dxl_old, PML_duz_dyl_old, PML_duz_dzl_old, &
                                                PML_dux_dxl_new, PML_dux_dyl_new, PML_dux_dzl_new, &
                                                PML_duy_dxl_new, PML_duy_dyl_new, PML_duy_dzl_new, &
                                                PML_duz_dxl_new, PML_duz_dyl_new, PML_duz_dzl_new, &
                                                rmemory_dux_dxl_x, rmemory_duy_dyl_x, rmemory_duz_dzl_x, &
                                                rmemory_dux_dyl_x, rmemory_dux_dzl_x, rmemory_duz_dxl_x, rmemory_duy_dxl_x, &
                                                rmemory_dux_dxl_y, rmemory_duz_dzl_y, rmemory_duy_dyl_y, &
                                                rmemory_duy_dxl_y, rmemory_duy_dzl_y, rmemory_duz_dyl_y, rmemory_dux_dyl_y, &
                                                rmemory_dux_dxl_z, rmemory_duy_dyl_z, rmemory_duz_dzl_z, &
                                                rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duy_dzl_z, rmemory_dux_dzl_z)

! calculates C-PML elastic memory variables and computes stress sigma

! second-order accurate convolution term calculation from equation (21) of
! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
! Anisotropic-medium PML for vector FETD with modified basis functions,
! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS

  use specfem_par, only: xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,jacobianstore, &
                         kappastore,mustore,irregular_element_number, &
                         jacobian_regular,xix_regular

  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z, &
                     d_store_x,d_store_y,d_store_z,alpha_store_x,alpha_store_y,alpha_store_z, &
                     convolution_coef_acoustic_alpha,convolution_coef_acoustic_beta

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempz1,tempz2,tempz3

  ! derivatives of ux, uy and uz with respect to x, y and z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_dux_dxl,PML_dux_dyl,PML_dux_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duy_dxl,PML_duy_dyl,PML_duy_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duz_dxl,PML_duz_dyl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_dux_dxl_new,PML_dux_dyl_new,PML_dux_dzl_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duy_dxl_new,PML_duy_dyl_new,PML_duy_dzl_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: PML_duz_dxl_new,PML_duz_dyl_new,PML_duz_dzl_new

  ! memory variable
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),intent(inout) ::  &
                          rmemory_dux_dxl_x, rmemory_dux_dyl_x, rmemory_dux_dzl_x, &
                          rmemory_duy_dxl_y, rmemory_duy_dyl_y, rmemory_duy_dzl_y, &
                          rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duz_dzl_z

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),intent(inout) ::  &
                          rmemory_duy_dyl_x, rmemory_duz_dzl_x, rmemory_duz_dxl_x, rmemory_duy_dxl_x, &
                          rmemory_dux_dxl_y, rmemory_duz_dzl_y, rmemory_duz_dyl_y, rmemory_dux_dyl_y, &
                          rmemory_dux_dxl_z, rmemory_duy_dyl_z, rmemory_duy_dzl_z, rmemory_dux_dzl_z
  ! local parameters
  integer :: i,j,k,ispec_irreg
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul,kappal
  real(kind=CUSTOM_REAL) :: duxdxl_x,duxdyl_x,duxdzl_x,duydxl_x,duydyl_x,duzdxl_x,duzdzl_x
  real(kind=CUSTOM_REAL) :: duxdxl_y,duxdyl_y,duydxl_y,duydyl_y,duydzl_y,duzdyl_y,duzdzl_y
  real(kind=CUSTOM_REAL) :: duxdxl_z,duxdzl_z,duydyl_z,duydzl_z,duzdxl_z,duzdyl_z,duzdzl_z
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9      ! L231
  real(kind=CUSTOM_REAL) :: A10,A11,A12,A13  ! L132
  real(kind=CUSTOM_REAL) :: A14,A15,A16,A17  ! L123
  real(kind=CUSTOM_REAL) :: A18,A19 ! L1
  real(kind=CUSTOM_REAL) :: A20,A21 ! L2
  real(kind=CUSTOM_REAL) :: A22,A23 ! L3
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  integer :: CPML_region_local
  real(kind=CUSTOM_REAL) :: kappa_x,kappa_y,kappa_z,d_x,d_y,d_z,alpha_x,alpha_y,alpha_z

  CPML_region_local = CPML_regions(ispec_CPML)

  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg == 0) jacobianl = jacobian_regular

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)

        !---------------------- A6, A7, A8, A9 --------------------------
        call lijk_parameter_computation(kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y,kappa_x,d_x,alpha_x, &
                                        CPML_region_local,231, &
                                        A6,A7,A8,A9)

        ! coefficients
        ! alpha_z
        coef0_1 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

        ! beta_x = alpha_x + d_x / kappa_x
        coef0_3 = convolution_coef_acoustic_beta(1,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(2,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(3,i,j,k,ispec_CPML)

        rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
               + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl_old(i,j,k) * coef2_1
        rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
               + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl_old(i,j,k) * coef2_1
        rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
               + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl_old(i,j,k) * coef2_1

        rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
               + PML_dux_dxl_new(i,j,k) * coef1_2 + PML_dux_dxl_old(i,j,k) * coef2_2
        rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
               + PML_duy_dxl_new(i,j,k) * coef1_2 + PML_duy_dxl_old(i,j,k) * coef2_2
        rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
               + PML_duz_dxl_new(i,j,k) * coef1_2 + PML_duz_dxl_old(i,j,k) * coef2_2

        rmemory_dux_dxl_x(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,3) &
               + PML_dux_dxl_new(i,j,k) * coef1_3 + PML_dux_dxl_old(i,j,k) * coef2_3
        rmemory_duy_dxl_y(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,3) &
               + PML_duy_dxl_new(i,j,k) * coef1_3 + PML_duy_dxl_old(i,j,k) * coef2_3
        rmemory_duz_dxl_z(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,3) &
               + PML_duz_dxl_new(i,j,k) * coef1_3 + PML_duz_dxl_old(i,j,k) * coef2_3


        !---------------------- A10,A11,A12,A13 --------------------------
        call lijk_parameter_computation(kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y, &
                                        CPML_region_local,132, &
                                        A10,A11,A12,A13)

        ! coefficients
        ! alpha_x
        coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

        ! alpha_z
        coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

        ! beta_y = alpha_y + d_y / kappa_y
        coef0_3 = convolution_coef_acoustic_beta(4,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(5,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(6,i,j,k,ispec_CPML)

        rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
               + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl_old(i,j,k) * coef2_1
        rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
               + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl_old(i,j,k) * coef2_1
        rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
               + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl_old(i,j,k) * coef2_1

        rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
               + PML_dux_dyl_new(i,j,k) * coef1_2 + PML_dux_dyl_old(i,j,k) * coef2_2
        rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
               + PML_duy_dyl_new(i,j,k) * coef1_2 + PML_duy_dyl_old(i,j,k) * coef2_2
        rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
               + PML_duz_dyl_new(i,j,k) * coef1_2 + PML_duz_dyl_old(i,j,k) * coef2_2

        rmemory_dux_dyl_x(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,3) &
               + PML_dux_dyl_new(i,j,k) * coef1_3 + PML_dux_dyl_old(i,j,k) * coef2_3
        rmemory_duy_dyl_y(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,3) &
               + PML_duy_dyl_new(i,j,k) * coef1_3 + PML_duy_dyl_old(i,j,k) * coef2_3
        rmemory_duz_dyl_z(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,3) &
               + PML_duz_dyl_new(i,j,k) * coef1_3 + PML_duz_dyl_old(i,j,k) * coef2_3

        !---------------------- A14,A15,A16,A17 --------------------------
        call lijk_parameter_computation(kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                        CPML_region_local,123, &
                                        A14,A15,A16,A17)

        ! coefficients
        ! alpha_x
        coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

        ! beta_z = alpha_z + d_z / kappa_z
        coef0_3 = convolution_coef_acoustic_beta(7,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(8,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(9,i,j,k,ispec_CPML)

        rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
               + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl_old(i,j,k) * coef2_1
        rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
               + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl_old(i,j,k) * coef2_1
        rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
               + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl_old(i,j,k) * coef2_1

        rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
               + PML_dux_dzl_new(i,j,k) * coef1_2 + PML_dux_dzl_old(i,j,k) * coef2_2
        rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
               + PML_duy_dzl_new(i,j,k) * coef1_2 + PML_duy_dzl_old(i,j,k) * coef2_2
        rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
               + PML_duz_dzl_new(i,j,k) * coef1_2 + PML_duz_dzl_old(i,j,k) * coef2_2

        rmemory_dux_dzl_x(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,3) &
               + PML_dux_dzl_new(i,j,k) * coef1_3 + PML_dux_dzl_old(i,j,k) * coef2_3
        rmemory_duy_dzl_y(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,3) &
               + PML_duy_dzl_new(i,j,k) * coef1_3 + PML_duy_dzl_old(i,j,k) * coef2_3
        rmemory_duz_dzl_z(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,3) &
               + PML_duz_dzl_new(i,j,k) * coef1_3 + PML_duz_dzl_old(i,j,k) * coef2_3

        !---------------------- A18 and A19 --------------------------
        call lx_parameter_computation(kappa_x,d_x,alpha_x, &
                                      CPML_region_local,A18,A19)

        ! coefficients
        ! alpha_x
        coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

        rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
               + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl_old(i,j,k) * coef2_1

        rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
               + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl_old(i,j,k) * coef2_1

        rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
               + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl_old(i,j,k) * coef2_1

        rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
               + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl_old(i,j,k) * coef2_1

        !---------------------- A20 and A21 --------------------------
        call ly_parameter_computation(kappa_y,d_y,alpha_y, &
                                      CPML_region_local,A20,A21)

        ! coefficients
        ! alpha_y
        coef0_1 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

        rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
               + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl_old(i,j,k) * coef2_1

        rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
               + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl_old(i,j,k) * coef2_1

        rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
               + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl_old(i,j,k) * coef2_1

        rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
               + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl_old(i,j,k) * coef2_1

        !---------------------- A22 and A23 --------------------------
        call lz_parameter_computation(kappa_z,d_z,alpha_z, &
                                      CPML_region_local,A22,A23)

        ! coefficients
        ! alpha_z
        coef0_1 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

        rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
               + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl_old(i,j,k) * coef2_1

        rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
               + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl_old(i,j,k) * coef2_1

        rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
               + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl_old(i,j,k) * coef2_1

        rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
               + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl_old(i,j,k) * coef2_1

        ! derivatives
        duxdxl_x = A6 * PML_dux_dxl(i,j,k) + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
                 + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) + A9 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,3)
        duxdyl_x = A10 * PML_dux_dyl(i,j,k) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
                   + A12 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) + A13 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,3)
        duxdzl_x = A14 * PML_dux_dzl(i,j,k) + A15 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
                 + A16 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) + A17 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,3)

        duzdzl_x = A20 * PML_duz_dzl(i,j,k) + A21 * rmemory_duz_dzl_x(i,j,k,ispec_CPML)
        duzdxl_x = A20 * PML_duz_dxl(i,j,k) + A21 * rmemory_duz_dxl_x(i,j,k,ispec_CPML)
        duydyl_x = A22 * PML_duy_dyl(i,j,k) + A23 * rmemory_duy_dyl_x(i,j,k,ispec_CPML)
        duydxl_x = A22 * PML_duy_dxl(i,j,k) + A23 * rmemory_duy_dxl_x(i,j,k,ispec_CPML)

        duydxl_y = A6 * PML_duy_dxl(i,j,k) + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
                 + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) + A9 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,3)
        duydyl_y = A10 * PML_duy_dyl(i,j,k) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
                 + A12 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) + A13 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,3)
        duydzl_y = A14 * PML_duy_dzl(i,j,k) + A15 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
                 + A16 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) + A17 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,3)

        duzdzl_y = A18 * PML_duz_dzl(i,j,k) + A19 * rmemory_duz_dzl_y(i,j,k,ispec_CPML)
        duzdyl_y = A18 * PML_duz_dyl(i,j,k) + A19 * rmemory_duz_dyl_y(i,j,k,ispec_CPML)
        duxdyl_y = A22 * PML_dux_dyl(i,j,k) + A23 * rmemory_dux_dyl_y(i,j,k,ispec_CPML)
        duxdxl_y = A22 * PML_dux_dxl(i,j,k) + A23 * rmemory_dux_dxl_y(i,j,k,ispec_CPML)

        duzdxl_z = A6 * PML_duz_dxl(i,j,k) + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
                 + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) + A9 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,3)
        duzdyl_z = A10 * PML_duz_dyl(i,j,k) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
                 + A12 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) + A13 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,3)
        duzdzl_z = A14 * PML_duz_dzl(i,j,k) + A15 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
                 + A16 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) + A17 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,3)

        duydzl_z = A18 * PML_duy_dzl(i,j,k) + A19 * rmemory_duy_dzl_z(i,j,k,ispec_CPML)
        duydyl_z = A18 * PML_duy_dyl(i,j,k) + A19 * rmemory_duy_dyl_z(i,j,k,ispec_CPML)
        duxdzl_z = A20 * PML_dux_dzl(i,j,k) + A21 * rmemory_dux_dzl_z(i,j,k,ispec_CPML)
        duxdxl_z = A20 * PML_dux_dxl(i,j,k) + A21 * rmemory_dux_dxl_z(i,j,k,ispec_CPML)

        ! elastic parameters
        kappal = kappastore(i,j,k,ispec)
        mul = mustore(i,j,k,ispec)

        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2.0_CUSTOM_REAL*mul

        ! compute stress sigma (non-symmetric)
        sigma_xx = lambdalplus2mul*duxdxl_x + lambdal*duydyl_x + lambdal*duzdzl_x
        sigma_yx = mul*duxdyl_x + mul*duydxl_x
        sigma_zx = mul*duzdxl_x + mul*duxdzl_x

        sigma_xy = mul*duxdyl_y + mul*duydxl_y
        sigma_yy = lambdal*duxdxl_y + lambdalplus2mul*duydyl_y + lambdal*duzdzl_y
        sigma_zy = mul*duzdyl_y + mul*duydzl_y

        sigma_xz = mul*duzdxl_z + mul*duxdzl_z
        sigma_yz = mul*duzdyl_z + mul*duydzl_z
        sigma_zz = lambdal*duxdxl_z + lambdal*duydyl_z + lambdalplus2mul*duzdzl_z

        if (ispec_irreg /= 0) then
          ! irregular element
          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)
          jacobianl = jacobianstore(i,j,k,ispec_irreg)

          ! form dot product with test vector, non-symmetric form (which
          ! is useful in the case of PML)
          tempx1(i,j,k) = jacobianl * (sigma_xx * xixl + sigma_yx * xiyl + sigma_zx * xizl) ! this goes to accel_x
          tempy1(i,j,k) = jacobianl * (sigma_xy * xixl + sigma_yy * xiyl + sigma_zy * xizl) ! this goes to accel_y
          tempz1(i,j,k) = jacobianl * (sigma_xz * xixl + sigma_yz * xiyl + sigma_zz * xizl) ! this goes to accel_z

          tempx2(i,j,k) = jacobianl * (sigma_xx * etaxl + sigma_yx * etayl + sigma_zx * etazl) ! this goes to accel_x
          tempy2(i,j,k) = jacobianl * (sigma_xy * etaxl + sigma_yy * etayl + sigma_zy * etazl) ! this goes to accel_y
          tempz2(i,j,k) = jacobianl * (sigma_xz * etaxl + sigma_yz * etayl + sigma_zz * etazl) ! this goes to accel_z

          tempx3(i,j,k) = jacobianl * (sigma_xx * gammaxl + sigma_yx * gammayl + sigma_zx * gammazl) ! this goes to accel_x
          tempy3(i,j,k) = jacobianl * (sigma_xy * gammaxl + sigma_yy * gammayl + sigma_zy * gammazl) ! this goes to accel_y
          tempz3(i,j,k) = jacobianl * (sigma_xz * gammaxl + sigma_yz * gammayl + sigma_zz * gammazl) ! this goes to accel_z
        else
          !regular element
          tempx1(i,j,k) = jacobianl * sigma_xx * xix_regular ! this goes to accel_x
          tempy1(i,j,k) = jacobianl * sigma_xy * xix_regular ! this goes to accel_y
          tempz1(i,j,k) = jacobianl * sigma_xz * xix_regular ! this goes to accel_z

          tempx2(i,j,k) = jacobianl * sigma_yx * xix_regular ! this goes to accel_x
          tempy2(i,j,k) = jacobianl * sigma_yy * xix_regular ! this goes to accel_y
          tempz2(i,j,k) = jacobianl * sigma_yz * xix_regular ! this goes to accel_z

          tempx3(i,j,k) = jacobianl * sigma_zx * xix_regular ! this goes to accel_x
          tempy3(i,j,k) = jacobianl * sigma_zy * xix_regular ! this goes to accel_y
          tempz3(i,j,k) = jacobianl * sigma_zz * xix_regular ! this goes to accel_z
        endif

      enddo
    enddo
  enddo

end subroutine pml_compute_memory_variables_elastic

!
!=====================================================================
!

subroutine pml_compute_memory_variables_acoustic(ispec,ispec_CPML, &
                                                 temp1,temp2,temp3, &
                                                 rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                                                 PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                                                 PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                                                 PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new)

! calculates C-PML elastic memory variables and computes stress sigma

! second-order accurate convolution term calculation from equation (21) of
! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
! Anisotropic-medium PML for vector FETD with modified basis functions,
! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,jacobianstore, &
                         rhostore,irregular_element_number, &
                         jacobian_regular,xix_regular

  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z, &
                     d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z, &
                     convolution_coef_acoustic_alpha,convolution_coef_acoustic_beta

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: temp1,temp2,temp3

  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),intent(inout) :: &
    rmemory_dpotential_dxl, rmemory_dpotential_dyl, rmemory_dpotential_dzl

  ! derivatives of potential with respect to x, y and z
  ! in computation potential_acoustic at "n" time step is used
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: &
    PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl

  ! in computation of PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old
  ! we replace potential_acoustic with potential_acoustic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: &
    PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old

  ! we replace potential_acoustic at "n" time step with
  ! we replace potential_acoustic with potential_acoustic_old with potential_acoustic_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: &
    PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: rho_invl_jacob
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9      ! L231
  real(kind=CUSTOM_REAL) :: A10,A11,A12,A13  ! L132
  real(kind=CUSTOM_REAL) :: A14,A15,A16,A17  ! L123
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  integer :: CPML_region_local
  integer :: i,j,k,ispec_irreg

  real(kind=CUSTOM_REAL) :: kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: d_x,d_y,d_z,alpha_x,alpha_y,alpha_z

  CPML_region_local = CPML_regions(ispec_CPML)

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)

        !---------------------- A6, A7, A8, A9 --------------------------
        call lijk_parameter_computation(kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y,kappa_x,d_x,alpha_x, &
                                        CPML_region_local,231, &
                                        A6,A7,A8,A9)

        ! coefficients
        ! alpha_z
        coef0_1 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

        ! beta_x = alpha_x + d_x / kappa_x
        coef0_3 = convolution_coef_acoustic_beta(1,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(2,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(3,i,j,k,ispec_CPML)

        rmemory_dpotential_dxl(1,i,j,k,ispec_CPML) = coef0_1 * rmemory_dpotential_dxl(1,i,j,k,ispec_CPML) + &
                coef1_1 * PML_dpotential_dxl_new(i,j,k) + coef2_1 * PML_dpotential_dxl_old(i,j,k)

        rmemory_dpotential_dxl(2,i,j,k,ispec_CPML) = coef0_2 * rmemory_dpotential_dxl(2,i,j,k,ispec_CPML) + &
                coef1_2 * PML_dpotential_dxl_new(i,j,k) + coef2_2 * PML_dpotential_dxl_old(i,j,k)

        rmemory_dpotential_dxl(3,i,j,k,ispec_CPML) = coef0_3 * rmemory_dpotential_dxl(3,i,j,k,ispec_CPML) + &
                coef1_3 * PML_dpotential_dxl_new(i,j,k) + coef2_3 * PML_dpotential_dxl_old(i,j,k)

        !---------------------- A10,A11,A12,A13 --------------------------
        call lijk_parameter_computation(kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y, &
                                        CPML_region_local,132, &
                                        A10,A11,A12,A13)

        ! coefficients
        ! alpha_x
        coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

        ! alpha_z
        coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

        ! beta_y = alpha_y + d_y / kappa_y
        coef0_3 = convolution_coef_acoustic_beta(4,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(5,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(6,i,j,k,ispec_CPML)

        rmemory_dpotential_dyl(1,i,j,k,ispec_CPML) = coef0_1 * rmemory_dpotential_dyl(1,i,j,k,ispec_CPML) + &
                coef1_1 * PML_dpotential_dyl_new(i,j,k) + coef2_1 * PML_dpotential_dyl_old(i,j,k)

        rmemory_dpotential_dyl(2,i,j,k,ispec_CPML) = coef0_2 * rmemory_dpotential_dyl(2,i,j,k,ispec_CPML) + &
                coef1_2 * PML_dpotential_dyl_new(i,j,k) + coef2_2 * PML_dpotential_dyl_old(i,j,k)

        rmemory_dpotential_dyl(3,i,j,k,ispec_CPML) = coef0_3 * rmemory_dpotential_dyl(3,i,j,k,ispec_CPML) + &
                coef1_3 * PML_dpotential_dyl_new(i,j,k) + coef2_3 * PML_dpotential_dyl_old(i,j,k)

        !---------------------- A14,A15,A16,A17 --------------------------
        call lijk_parameter_computation(kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                        CPML_region_local,123, &
                                        A14,A15,A16,A17)

        ! coefficients
        ! alpha_x
        coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
        coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
        coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

        ! alpha_y
        coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
        coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
        coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

        ! beta_z = alpha_z + d_z / kappa_z
        coef0_3 = convolution_coef_acoustic_beta(7,i,j,k,ispec_CPML)
        coef1_3 = convolution_coef_acoustic_beta(8,i,j,k,ispec_CPML)
        coef2_3 = convolution_coef_acoustic_beta(9,i,j,k,ispec_CPML)

        rmemory_dpotential_dzl(1,i,j,k,ispec_CPML) = coef0_1 * rmemory_dpotential_dzl(1,i,j,k,ispec_CPML) + &
                coef1_1 * PML_dpotential_dzl_new(i,j,k) + coef2_1 * PML_dpotential_dzl_old(i,j,k)

        rmemory_dpotential_dzl(2,i,j,k,ispec_CPML) = coef0_2 * rmemory_dpotential_dzl(2,i,j,k,ispec_CPML) + &
                coef1_2 * PML_dpotential_dzl_new(i,j,k) + coef2_2 * PML_dpotential_dzl_old(i,j,k)

        rmemory_dpotential_dzl(3,i,j,k,ispec_CPML) = coef0_3 * rmemory_dpotential_dzl(3,i,j,k,ispec_CPML) + &
                coef1_3 * PML_dpotential_dzl_new(i,j,k) + coef2_3 * PML_dpotential_dzl_old(i,j,k)

        ! derivatives
        dpotentialdxl(i,j,k) = A6 * PML_dpotential_dxl(i,j,k)  + &
                        A7 * rmemory_dpotential_dxl(1,i,j,k,ispec_CPML) + &
                        A8 * rmemory_dpotential_dxl(2,i,j,k,ispec_CPML) + &
                        A9 * rmemory_dpotential_dxl(3,i,j,k,ispec_CPML)
        dpotentialdyl(i,j,k) = A10 * PML_dpotential_dyl(i,j,k) + &
                        A11 * rmemory_dpotential_dyl(1,i,j,k,ispec_CPML) + &
                        A12 * rmemory_dpotential_dyl(2,i,j,k,ispec_CPML) + &
                        A13 * rmemory_dpotential_dyl(3,i,j,k,ispec_CPML)
        dpotentialdzl(i,j,k) = A14 * PML_dpotential_dzl(i,j,k) + &
                        A15 * rmemory_dpotential_dzl(1,i,j,k,ispec_CPML) + &
                        A16 * rmemory_dpotential_dzl(2,i,j,k,ispec_CPML) + &
                        A17 * rmemory_dpotential_dzl(3,i,j,k,ispec_CPML)
      enddo
    enddo
  enddo

  ispec_irreg = irregular_element_number(ispec)
  if (ispec_irreg /= 0) then
    ! irregular element
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)
          jacobianl = jacobianstore(i,j,k,ispec_irreg)

          rho_invl_jacob = jacobianl / rhostore(i,j,k,ispec)

          temp1(i,j,k) = rho_invl_jacob &
            * (xixl*dpotentialdxl(i,j,k) + xiyl*dpotentialdyl(i,j,k) + xizl*dpotentialdzl(i,j,k))
          temp2(i,j,k) = rho_invl_jacob &
            * (etaxl*dpotentialdxl(i,j,k) + etayl*dpotentialdyl(i,j,k) + etazl*dpotentialdzl(i,j,k))
          temp3(i,j,k) = rho_invl_jacob &
            * (gammaxl*dpotentialdxl(i,j,k) + gammayl*dpotentialdyl(i,j,k) + gammazl*dpotentialdzl(i,j,k))
        enddo
      enddo
    enddo
  else
    ! regular element
    jacobianl = jacobian_regular
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          rho_invl_jacob = jacobianl / rhostore(i,j,k,ispec)

          temp1(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdxl(i,j,k)
          temp2(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdyl(i,j,k)
          temp3(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdzl(i,j,k)
        enddo
      enddo
    enddo
  endif

end subroutine pml_compute_memory_variables_acoustic

!
!=====================================================================
!

subroutine pml_compute_memory_variables_acoustic_elastic(ispec_CPML,iface,iglob,i,j,k, &
                                                         displ_x,displ_y,displ_z,displ, &
                                                         num_coupling_ac_el_faces,rmemory_coupling_ac_el_displ)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: NGLOB_AB

  use pml_par, only: CPML_regions,k_store_x,k_store_y,k_store_z,d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z, &
                     convolution_coef_acoustic_alpha, &
                     PML_displ_old,PML_displ_new

  implicit none

  integer, intent(in) :: ispec_CPML,iface,iglob,num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL),intent(out) :: displ_x,displ_y,displ_z
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2),intent(inout) :: &
    rmemory_coupling_ac_el_displ

  ! local parameters
  integer :: i,j,k,CPML_region_local
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A_12,A_13,A_14
  real(kind=CUSTOM_REAL) :: kappa_x,kappa_y,kappa_z,d_x,d_y,d_z,alpha_x,alpha_y,alpha_z

  CPML_region_local = CPML_regions(ispec_CPML)
  kappa_x = k_store_x(i,j,k,ispec_CPML)
  kappa_y = k_store_y(i,j,k,ispec_CPML)
  kappa_z = k_store_z(i,j,k,ispec_CPML)
  d_x = d_store_x(i,j,k,ispec_CPML)
  d_y = d_store_y(i,j,k,ispec_CPML)
  d_z = d_store_z(i,j,k,ispec_CPML)
  alpha_x = alpha_store_x(i,j,k,ispec_CPML)
  alpha_y = alpha_store_y(i,j,k,ispec_CPML)
  alpha_z = alpha_store_z(i,j,k,ispec_CPML)


  call lxy_interface_parameter_computation(kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                           CPML_region_local,23, &
                                           A_12,A_13,A_14)

  ! coefficients
  ! alpha_y
  coef0_1 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

  ! alpha_z
  coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

  ! displ_x
  rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) + &
                                                  coef1_1 * PML_displ_new(1,i,j,k,ispec_CPML) + &
                                                  coef2_1 * PML_displ_old(1,i,j,k,ispec_CPML)

  rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) + &
                                                  coef1_2 * PML_displ_new(1,i,j,k,ispec_CPML) + &
                                                  coef2_2 * PML_displ_old(1,i,j,k,ispec_CPML)

  displ_x = A_12 * displ(1,iglob) + A_13 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) + &
                                    A_14 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2)

  call lxy_interface_parameter_computation(kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z, &
                                           CPML_region_local,13, &
                                           A_12,A_13,A_14)

  ! coefficients
  ! alpha_x
  coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

  ! alpha_z
  coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

  ! displ_y
  rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) + &
                                                  coef1_1 * PML_displ_new(2,i,j,k,ispec_CPML) + &
                                                  coef2_1 * PML_displ_old(2,i,j,k,ispec_CPML)

  rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) + &
                                                  coef1_2 * PML_displ_new(2,i,j,k,ispec_CPML) + &
                                                  coef2_2 * PML_displ_old(2,i,j,k,ispec_CPML)

  displ_y = A_12 * displ(2,iglob) + A_13 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) + &
                                    A_14 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2)

  call lxy_interface_parameter_computation(kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y, &
                                           CPML_region_local,12, &
                                           A_12,A_13,A_14)
  ! coefficients
  ! alpha_x
  coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

  ! alpha_y
  coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

  ! displ_z
  rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) + &
                                                  coef1_1 * PML_displ_new(3,i,j,k,ispec_CPML) + &
                                                  coef2_1 * PML_displ_old(3,i,j,k,ispec_CPML)

  rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) + &
                                                  coef1_2 * PML_displ_new(3,i,j,k,ispec_CPML) + &
                                                  coef2_2 * PML_displ_old(3,i,j,k,ispec_CPML)


  displ_z = A_12 * displ(3,iglob) + A_13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) + &
                                    A_14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)

end subroutine pml_compute_memory_variables_acoustic_elastic

!
!=====================================================================
!

subroutine pml_compute_memory_variables_elastic_acoustic(ispec_CPML,iface,iglob,i,j,k, &
                                                         pressure_x,pressure_y,pressure_z, &
                                                         potential_acoustic,potential_dot_acoustic, &
                                                         potential_dot_dot_acoustic,num_coupling_ac_el_faces, &
                                                         rmemory_coupling_el_ac_potential)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: NGLOB_AB

  use pml_par, only: CPML_regions,k_store_x,k_store_y,k_store_z,d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z, &
                     convolution_coef_acoustic_alpha, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  integer, intent(in) :: ispec_CPML,iface,iglob,num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL),intent(out) :: pressure_x,pressure_y,pressure_z
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic, &
                                                             potential_dot_acoustic,potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2) :: &
                                    rmemory_coupling_el_ac_potential

  ! local parameters
  integer :: i,j,k,CPML_region_local
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A_0, A_1, A_2, A_3, A_4
  real(kind=CUSTOM_REAL) :: kappa_x,kappa_y,kappa_z,d_x,d_y,d_z,alpha_x,alpha_y,alpha_z

  CPML_region_local = CPML_regions(ispec_CPML)
  kappa_x = k_store_x(i,j,k,ispec_CPML)
  kappa_y = k_store_y(i,j,k,ispec_CPML)
  kappa_z = k_store_z(i,j,k,ispec_CPML)
  d_x = d_store_x(i,j,k,ispec_CPML)
  d_y = d_store_y(i,j,k,ispec_CPML)
  d_z = d_store_z(i,j,k,ispec_CPML)
  alpha_x = alpha_store_x(i,j,k,ispec_CPML)
  alpha_y = alpha_store_y(i,j,k,ispec_CPML)
  alpha_z = alpha_store_z(i,j,k,ispec_CPML)

  call l_interface_parameter_computation(kappa_y, d_y, alpha_y, kappa_z, d_z, alpha_z, &
                                         CPML_region_local, 23, &
                                         A_0, A_1, A_2, A_3, A_4)

  ! coefficients
  ! alpha_y
  coef0_1 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

  ! alpha_z
  coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) + &
                                                    coef1_1 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_1 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,2) + &
                                                    coef1_2 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_2 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  pressure_x = A_0 * potential_dot_dot_acoustic(iglob) + A_1 * potential_dot_acoustic(iglob)  + A_2 * potential_acoustic(iglob) &
             + A_3 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) + A_4 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,2)

  call l_interface_parameter_computation(kappa_x, d_x, alpha_x, kappa_z, d_z, alpha_z, &
                                         CPML_region_local, 13, &
                                         A_0, A_1, A_2, A_3, A_4)

  ! coefficients
  ! alpha_x
  coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

  ! alpha_z
  coef0_2 = convolution_coef_acoustic_alpha(7,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(8,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(9,i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) + &
                                                    coef1_1 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_1 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,2) + &
                                                    coef1_2 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_2 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  pressure_y = A_0 * potential_dot_dot_acoustic(iglob) + A_1 * potential_dot_acoustic(iglob)  + A_2 * potential_acoustic(iglob) &
             + A_3 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) + A_4 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,2)

  call l_interface_parameter_computation(kappa_x, d_x, alpha_x, kappa_y, d_y, alpha_y, &
                                         CPML_region_local, 12, &
                                         A_0, A_1, A_2, A_3, A_4)

  ! coefficients
  ! alpha_x
  coef0_1 = convolution_coef_acoustic_alpha(1,i,j,k,ispec_CPML)
  coef1_1 = convolution_coef_acoustic_alpha(2,i,j,k,ispec_CPML)
  coef2_1 = convolution_coef_acoustic_alpha(3,i,j,k,ispec_CPML)

  ! alpha_y
  coef0_2 = convolution_coef_acoustic_alpha(4,i,j,k,ispec_CPML)
  coef1_2 = convolution_coef_acoustic_alpha(5,i,j,k,ispec_CPML)
  coef2_2 = convolution_coef_acoustic_alpha(6,i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(3,i,j,k,iface,1) + &
                                                    coef1_1 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_1 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(3,i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(3,i,j,k,iface,2) + &
                                                    coef1_2 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_2 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  pressure_z = A_0 * potential_dot_dot_acoustic(iglob) + A_1 * potential_dot_acoustic(iglob)  + A_2 * potential_acoustic(iglob) &
             + A_3 * rmemory_coupling_el_ac_potential(3,i,j,k,iface,1) + A_4 * rmemory_coupling_el_ac_potential(3,i,j,k,iface,2)


end subroutine pml_compute_memory_variables_elastic_acoustic

!
!=====================================================================
!

subroutine lijk_parameter_computation(kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                      CPML_region_local,index_ijk, &
                                      A_0,A_1,A_2,A_3)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: kappa_x,d_x,alpha_x, &
                                        kappa_y,d_y,alpha_y, &
                                        kappa_z,d_z,alpha_z

  integer, intent(in) :: CPML_region_local,index_ijk

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1,A_2,A_3

  ! local variables
  real(kind=CUSTOM_REAL) :: beta_x,beta_y,beta_z
  real(kind=CUSTOM_REAL) :: bar_A_0,bar_A_1,bar_A_2,bar_A_3

  integer :: CPML_X_ONLY_TEMP,CPML_Y_ONLY_TEMP,CPML_Z_ONLY_TEMP, &
             CPML_XY_ONLY_TEMP,CPML_XZ_ONLY_TEMP,CPML_YZ_ONLY_TEMP

  select case (index_ijk)
  case (123)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XY_ONLY_TEMP = CPML_XY_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (132)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_Y_ONLY
    CPML_XY_ONLY_TEMP = CPML_XZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (231)
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XY_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_XY_ONLY
  case default
    stop 'Invalid index in lijk_parameter_computation, index_ijk must be equal to 123 or 132 or 231'
  end select

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y
  beta_z = alpha_z + d_z / kappa_z

  !unused, stored now in pre-computed arrays
  !call compute_convolution_coef(alpha_x, coef0_1, coef1_1, coef2_1)
  !call compute_convolution_coef(alpha_y, coef0_2, coef1_2, coef2_2)
  !call compute_convolution_coef(beta_z,  coef0_3, coef1_3, coef2_3)

  if (CPML_region_local == CPML_XYZ) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter .and. &
        abs(alpha_x-beta_z) >= min_distance_between_CPML_parameter .and. &
        abs(alpha_y-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y / kappa_z

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x-alpha_z) * (alpha_x-beta_x) * (alpha_x-beta_y) / ((alpha_x-alpha_y) * (alpha_x-beta_z))
      bar_A_2 = - bar_A_0 * (alpha_y-alpha_z) * (alpha_y-beta_x) * (alpha_y-beta_y) / ((alpha_y-alpha_x) * (alpha_y-beta_z))
      bar_A_3 = - bar_A_0 * (beta_z-alpha_z)  * (beta_z-beta_x)  * (beta_z-beta_y)  / ((beta_z-alpha_x)  * (beta_z-alpha_y))

    else
      stop 'Error in lijk_parameter_computation XYZ'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY_TEMP) then

    if (abs(alpha_y-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_y / kappa_z

      !----------------A1,2,3-------------------------
      bar_A_1 = 0._CUSTOM_REAL
      bar_A_2 = - bar_A_0 * (alpha_y - alpha_z) * (alpha_y - beta_y) / (alpha_y-beta_z)
      bar_A_3 = - bar_A_0 * (beta_z - alpha_z)  * (beta_z - beta_y)  / (beta_z-alpha_y)

    else
      stop 'Error in lijk_parameter_computation YZ'
    endif

  else if (CPML_region_local == CPML_XZ_ONLY_TEMP) then

    if (abs(alpha_x-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x / kappa_z

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - alpha_z) * (alpha_x - beta_x) / (alpha_x-beta_z)
      bar_A_2 = 0._CUSTOM_REAL
      bar_A_3 = - bar_A_0 * (beta_z - alpha_z)  * (beta_z - beta_x)  / (beta_z-alpha_x)

    else
      stop 'Error in lijk_parameter_computation XZ'
    endif

  else if (CPML_region_local == CPML_XY_ONLY_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / (alpha_y-alpha_x)
      bar_A_3 = 0._CUSTOM_REAL

    else
      stop 'Error in lijk_parameter_computation XY'
    endif
  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL
    bar_A_3 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)
    bar_A_3 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = 1._CUSTOM_REAL / kappa_z

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = 0._CUSTOM_REAL
    bar_A_3 = - bar_A_0 * (beta_z - alpha_z)

  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2
  A_3 = bar_A_3

end subroutine lijk_parameter_computation

!
!=====================================================================
!

subroutine lx_parameter_computation(kappa_x,d_x,alpha_x, &
                                    CPML_region_local,A_0,A_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: kappa_x,d_x,alpha_x
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1

  !local variable
  real(kind=CUSTOM_REAL) :: beta_x

  beta_x = alpha_x + d_x / kappa_x

  if (CPML_region_local == CPML_XYZ) then
  !----------------A0-------------------------
    A_0 = kappa_x
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_x - beta_x)

  else if (CPML_region_local == CPML_YZ_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_XZ_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_x
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_x - beta_x)

  else if (CPML_region_local == CPML_XY_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_x
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_x - beta_x)

  else if (CPML_region_local == CPML_X_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_x
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_x - beta_x)

  else if (CPML_region_local == CPML_Y_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Z_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  endif

  ! unused, now uses pre-computed arrays
  !call compute_convolution_coef(alpha_x, coef0_1, coef1_1, coef2_1)

end subroutine lx_parameter_computation

!
!=====================================================================
!

subroutine ly_parameter_computation(kappa_y,d_y,alpha_y, &
                                    CPML_region_local,A_0,A_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: kappa_y,d_y,alpha_y
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1

  !local variable
  real(kind=CUSTOM_REAL) :: beta_y

  beta_y = alpha_y + d_y / kappa_y

  if (CPML_region_local == CPML_XYZ) then
  !----------------A0-------------------------
    A_0 = kappa_y
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_YZ_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_y
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_XZ_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_XY_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_y
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_X_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_y
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_Z_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  endif

  !unused, now uses pre-computed arrays
  !call compute_convolution_coef(alpha_y, coef0_1, coef1_1, coef2_1)

end subroutine ly_parameter_computation

!
!=====================================================================
!

subroutine lz_parameter_computation(kappa_z,d_z,alpha_z, &
                                    CPML_region_local,A_0,A_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: kappa_z,d_z,alpha_z
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1

  !local variable
  real(kind=CUSTOM_REAL) :: beta_z

  beta_z = alpha_z + d_z / kappa_z

  if (CPML_region_local == CPML_XYZ) then
  !----------------A0-------------------------
    A_0 = kappa_z
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_z - beta_z)

  else if (CPML_region_local == CPML_YZ_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_z
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_z - beta_z)

  else if (CPML_region_local == CPML_XZ_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_z
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_z - beta_z)

  else if (CPML_region_local == CPML_XY_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_X_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY) then
  !----------------A0-------------------------
    A_0 = 1._CUSTOM_REAL
  !----------------A1-------------------------
    A_1 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Z_ONLY) then
  !----------------A0-------------------------
    A_0 = kappa_z
  !----------------A1-------------------------
    A_1 = - A_0 * (alpha_z - beta_z)

  endif

  ! unused, now uses pre-computed arrays
  !call compute_convolution_coef(alpha_z, coef0_1, coef1_1, coef2_1)

end subroutine lz_parameter_computation

!
!=====================================================================
!

subroutine lxy_interface_parameter_computation(kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y, &
                                               CPML_region_local,index_ijk, &
                                               A_0,A_1,A_2)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x, &
                            kappa_y,d_y,alpha_y
  integer, intent(in) :: CPML_region_local,index_ijk

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1,A_2

  ! local variables
  real(kind=CUSTOM_REAL) :: bar_A_0,bar_A_1,bar_A_2,beta_x,beta_y

  integer :: CPML_X_ONLY_TEMP,CPML_Y_ONLY_TEMP,CPML_Z_ONLY_TEMP, &
             CPML_XY_ONLY_TEMP,CPML_XZ_ONLY_TEMP,CPML_YZ_ONLY_TEMP

  select case(index_ijk)
  case (12)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XY_ONLY_TEMP = CPML_XY_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (13)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_Y_ONLY
    CPML_XY_ONLY_TEMP = CPML_XZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (23)
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XY_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_XY_ONLY
  case default
    stop 'In lxy_interface_parameter_computation index_ijk must be equal to 12 or 13 or 23'
  end select

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y

  !unused, now pre-computed arrays
  !call compute_convolution_coef(alpha_x, coef0_1, coef1_1, coef2_1)
  !call compute_convolution_coef(alpha_y, coef0_2, coef1_2, coef2_2)

  if (CPML_region_local == CPML_XYZ) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / (alpha_y-alpha_x)
    else
      stop 'Error in lxy_interface_parameter_computation XYZ'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_XZ_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_XY_ONLY_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / (alpha_y-alpha_x)
    else
      stop 'Error in lxy_interface_parameter_computation XY'
    endif

  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = 1._CUSTOM_REAL

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = 0._CUSTOM_REAL

  else
    stop 'Error in lxy_interface_parameter_computation'
  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2

end subroutine lxy_interface_parameter_computation

!
!=====================================================================
!

subroutine l_interface_parameter_computation(kappa_x, d_x, alpha_x, &
                                             kappa_y, d_y, alpha_y, &
                                             CPML_region_local,index_ijk, &
                                             A_0, A_1, A_2, A_3, A_4)

  use constants, only: CUSTOM_REAL, CPML_XYZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x,beta_x, &
                            kappa_y,d_y,alpha_y,beta_y
  integer, intent(in) :: CPML_region_local,index_ijk

  real(kind=CUSTOM_REAL), intent(out) :: A_0, A_1, A_2, A_3, A_4

  !local variable
  real(kind=CUSTOM_REAL) :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4

  integer :: CPML_X_ONLY_TEMP,CPML_Y_ONLY_TEMP,CPML_Z_ONLY_TEMP, &
             CPML_XY_ONLY_TEMP,CPML_XZ_ONLY_TEMP,CPML_YZ_ONLY_TEMP

  select case(index_ijk)
  case (12)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XY_ONLY_TEMP = CPML_XY_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (13)
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_Y_ONLY
    CPML_XY_ONLY_TEMP = CPML_XZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
  case (23)
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XY_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_XY_ONLY
  case default
    stop 'Invalid index in l_interface_parameter_computation, index_ijk must be equal to 12 or 13 or 23'
  end select

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y

  !unused, now uses pre-computed arrays
  !call compute_convolution_coef(alpha_x, coef0_x, coef1_x, coef2_x)
  !call compute_convolution_coef(alpha_y, coef0_y, coef1_y, coef2_y)

  if (CPML_region_local == CPML_XYZ) then

    if ( abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter) then
      bar_A_0 = kappa_x * kappa_y
      bar_A_1 = bar_A_0 * (beta_x + beta_y - alpha_x - alpha_y)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_y - alpha_y - alpha_x) &
              + bar_A_0 * (beta_y - alpha_y) * (- alpha_y)
      bar_A_3 = bar_A_0 * alpha_x**2 &
              * (beta_x - alpha_x) * (beta_y - alpha_x) / (alpha_y - alpha_x)
      bar_A_4 = bar_A_0 * alpha_y**2 &
              * (beta_x - alpha_y) * (beta_y - alpha_y) / (alpha_x - alpha_y)
    else
      stop 'Error occured in l_parameter_computation in CPML_XYZ region'
    endif

  else if (CPML_region_local == CPML_XY_ONLY_TEMP) then

    if (abs( alpha_x - alpha_y ) >= min_distance_between_CPML_parameter) then
      bar_A_0 = kappa_x * kappa_y
      bar_A_1 = bar_A_0 * (beta_x + beta_y - alpha_x - alpha_y)
      bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_y - alpha_y - alpha_x) &
              - bar_A_0 * (beta_y - alpha_y) * alpha_y
      bar_A_3 = bar_A_0 * alpha_x**2 &
              * (beta_x - alpha_x) * (beta_y - alpha_x) / (alpha_y - alpha_x)
      bar_A_4 = bar_A_0 * alpha_y**2 &
              * (beta_x - alpha_y) * (beta_y - alpha_y) / (alpha_x - alpha_y)
    else
      stop 'Error occured in l_parameter_computation in CPML_XY_ONLY region'
    endif

  else if (CPML_region_local == CPML_XZ_ONLY_TEMP) then

    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (- alpha_x)
    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_YZ_ONLY_TEMP) then

    bar_A_0 = kappa_y
    bar_A_1 = bar_A_0 * (beta_y- alpha_y)
    bar_A_2 = bar_A_0 * (beta_y - alpha_y) * (- alpha_y)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = bar_A_0 * alpha_y**2 * (beta_y - alpha_y)

  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = - bar_A_0 * alpha_x * (beta_x - alpha_x)
    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0._CUSTOM_REAL

  else if (CPML_region_local == CPML_Y_ONLY_TEMP) then

    bar_A_0 = kappa_y
    bar_A_1 = bar_A_0 * (beta_y - alpha_y)
    bar_A_2 = - bar_A_0 * alpha_y * (beta_y - alpha_y)
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = bar_A_0 * alpha_y**2 * (beta_y - alpha_y)

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    bar_A_0 = 1._CUSTOM_REAL
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = 0._CUSTOM_REAL
    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = 0._CUSTOM_REAL

  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2
  A_3 = bar_A_3
  A_4 = bar_A_4

end subroutine l_interface_parameter_computation

