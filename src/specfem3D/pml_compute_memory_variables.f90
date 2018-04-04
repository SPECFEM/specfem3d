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
!
! United States and French Government Sponsorship Acknowledged.

subroutine pml_compute_memory_variables_elastic(ispec,ispec_CPML,tempx1,tempy1,tempz1,tempx2,tempy2,tempz2, &
                                    tempx3,tempy3,tempz3, &
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

  use specfem_par, only: deltat,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
                         kappastore,mustore,irregular_element_number, &
                         jacobian_regular,xix_regular

  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z, &
                     d_store_x,d_store_y,d_store_z,alpha_store_x,alpha_store_y,alpha_store_z, &
                     PML_dux_dxl, PML_dux_dyl, PML_dux_dzl, PML_duy_dxl, PML_duy_dyl, PML_duy_dzl, &
                     PML_duz_dxl, PML_duz_dyl, PML_duz_dzl, &
                     PML_dux_dxl_old, PML_dux_dyl_old, PML_dux_dzl_old, &
                     PML_duy_dxl_old, PML_duy_dyl_old, PML_duy_dzl_old, &
                     PML_duz_dxl_old, PML_duz_dyl_old, PML_duz_dzl_old, &
                     PML_dux_dxl_new, PML_dux_dyl_new, PML_dux_dzl_new, &
                     PML_duy_dxl_new, PML_duy_dyl_new, PML_duy_dzl_new, &
                     PML_duz_dxl_new, PML_duz_dyl_new, PML_duz_dzl_new

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3) ::  &
                          rmemory_dux_dxl_x, rmemory_dux_dyl_x, rmemory_dux_dzl_x, &
                          rmemory_duy_dxl_y, rmemory_duy_dyl_y, rmemory_duy_dzl_y, &
                          rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duz_dzl_z

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML) ::  &
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
        kappal = kappastore(i,j,k,ispec)
        mul = mustore(i,j,k,ispec)
        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2.0_CUSTOM_REAL*mul

        if (ispec_irreg /= 0) then !irregular element

          xixl = xix(i,j,k,ispec_irreg)
          xiyl = xiy(i,j,k,ispec_irreg)
          xizl = xiz(i,j,k,ispec_irreg)
          etaxl = etax(i,j,k,ispec_irreg)
          etayl = etay(i,j,k,ispec_irreg)
          etazl = etaz(i,j,k,ispec_irreg)
          gammaxl = gammax(i,j,k,ispec_irreg)
          gammayl = gammay(i,j,k,ispec_irreg)
          gammazl = gammaz(i,j,k,ispec_irreg)
          jacobianl = jacobian(i,j,k,ispec_irreg)

        endif

        !---------------------- A6, A7, A8, A9 --------------------------
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)
        call lijk_parameter_computation(deltat, &
                                        kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y,kappa_x,d_x,alpha_x, &
                                        CPML_region_local,231,A6,A7,A8,A9, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)
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
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)
        call lijk_parameter_computation(deltat, &
                                        kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y, &
                                        CPML_region_local,132,A10,A11,A12,A13, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)
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
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)
        call lijk_parameter_computation(deltat, &
                                        kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                        CPML_region_local,123,A14,A15,A16,A17, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)

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
        kappa_x = k_store_x(i,j,k,ispec_CPML)
        d_x = d_store_x(i,j,k,ispec_CPML)
        alpha_x = alpha_store_x(i,j,k,ispec_CPML)
        call lx_parameter_computation(deltat,kappa_x,d_x,alpha_x, &
                                      CPML_region_local,A18,A19, &
                                      coef0_1,coef1_1,coef2_1)

        rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
               + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl_old(i,j,k) * coef2_1

        rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
               + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl_old(i,j,k) * coef2_1

        rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
               + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl_old(i,j,k) * coef2_1

        rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
               + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl_old(i,j,k) * coef2_1

        !---------------------- A20 and A21 --------------------------
        kappa_y = k_store_y(i,j,k,ispec_CPML)
        d_y = d_store_y(i,j,k,ispec_CPML)
        alpha_y = alpha_store_y(i,j,k,ispec_CPML)
        call ly_parameter_computation(deltat,kappa_y,d_y,alpha_y, &
                                      CPML_region_local,A20,A21, &
                                      coef0_1,coef1_1,coef2_1)

        rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
               + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl_old(i,j,k) * coef2_1

        rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
               + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl_old(i,j,k) * coef2_1

        rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
               + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl_old(i,j,k) * coef2_1

        rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
               + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl_old(i,j,k) * coef2_1

        !---------------------- A22 and A23 --------------------------
        kappa_z = k_store_z(i,j,k,ispec_CPML)
        d_z = d_store_z(i,j,k,ispec_CPML)
        alpha_z = alpha_store_z(i,j,k,ispec_CPML)
        call lz_parameter_computation(deltat,kappa_z,d_z,alpha_z, &
                                      CPML_region_local,A22,A23, &
                                      coef0_1,coef1_1,coef2_1)

        rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
               + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl_old(i,j,k) * coef2_1

        rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
               + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl_old(i,j,k) * coef2_1

        rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
               + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl_old(i,j,k) * coef2_1

        rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
               + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl_old(i,j,k) * coef2_1

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

        ! compute stress sigma
        sigma_xx = lambdalplus2mul*duxdxl_x + lambdal*duydyl_x + lambdal*duzdzl_x
        sigma_yx = mul*duxdyl_x + mul*duydxl_x
        sigma_zx = mul*duzdxl_x + mul*duxdzl_x

        sigma_xy = mul*duxdyl_y + mul*duydxl_y
        sigma_yy = lambdal*duxdxl_y + lambdalplus2mul*duydyl_y + lambdal*duzdzl_y
        sigma_zy = mul*duzdyl_y + mul*duydzl_y

        sigma_xz = mul*duzdxl_z + mul*duxdzl_z
        sigma_yz = mul*duzdyl_z + mul*duydzl_z
        sigma_zz = lambdal*duxdxl_z + lambdal*duydyl_z + lambdalplus2mul*duzdzl_z

        if (ispec_irreg /= 0) then ! irregular element

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
        else !regular element
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
subroutine pml_compute_memory_variables_acoustic(ispec,ispec_CPML,temp1,temp2,temp3, &
                                                 rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                                                 PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                                                 PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                                                 PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
                         deltat,rhostore,irregular_element_number, &
                         jacobian_regular,xix_regular

  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z, &
                     d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: temp1,temp2,temp3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),intent(inout) :: &
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
  integer :: i,j,k,ispec_irreg
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: rho_invl_jacob
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9      ! L231
  real(kind=CUSTOM_REAL) :: A10,A11,A12,A13  ! L132
  real(kind=CUSTOM_REAL) :: A14,A15,A16,A17  ! L123
  real(kind=CUSTOM_REAL) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2,coef0_3,coef1_3,coef2_3
  integer :: CPML_region_local

  real(kind=CUSTOM_REAL) :: kappa_x,kappa_y,kappa_z
  real(kind=CUSTOM_REAL) :: d_x,d_y,d_z,alpha_x,alpha_y,alpha_z


  CPML_region_local = CPML_regions(ispec_CPML)
  ispec_irreg = irregular_element_number(ispec)

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
        call lijk_parameter_computation(deltat, &
                                        kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y,kappa_x,d_x,alpha_x, &
                                        CPML_region_local,231,A6,A7,A8,A9, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)

        rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + &
                coef1_1 * PML_dpotential_dxl_new(i,j,k) + coef2_1 * PML_dpotential_dxl_old(i,j,k)

        rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) + &
                coef1_2 * PML_dpotential_dxl_new(i,j,k) + coef2_2 * PML_dpotential_dxl_old(i,j,k)

        rmemory_dpotential_dxl(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,3) + &
                coef1_3 * PML_dpotential_dxl_new(i,j,k) + coef2_3 * PML_dpotential_dxl_old(i,j,k)

        !---------------------- A10,A11,A12,A13 --------------------------
        call lijk_parameter_computation(deltat, &
                                        kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z,kappa_y,d_y,alpha_y, &
                                        CPML_region_local,132,A10,A11,A12,A13, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)

        rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + &
                coef1_1 * PML_dpotential_dyl_new(i,j,k) + coef2_1 * PML_dpotential_dyl_old(i,j,k)

        rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) + &
                coef1_2 * PML_dpotential_dyl_new(i,j,k) + coef2_2 * PML_dpotential_dyl_old(i,j,k)

        rmemory_dpotential_dyl(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,3) + &
                coef1_3 * PML_dpotential_dyl_new(i,j,k) + coef2_3 * PML_dpotential_dyl_old(i,j,k)

        !---------------------- A14,A15,A16,A17 --------------------------
        call lijk_parameter_computation(deltat, &
                                        kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                        CPML_region_local,123,A14,A15,A16,A17, &
                                        coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                        coef0_3,coef1_3,coef2_3)

        rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + &
                coef1_1 * PML_dpotential_dzl_new(i,j,k) + coef2_1 * PML_dpotential_dzl_old(i,j,k)

        rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) + &
                coef1_2 * PML_dpotential_dzl_new(i,j,k) + coef2_2 * PML_dpotential_dzl_old(i,j,k)

        rmemory_dpotential_dzl(i,j,k,ispec_CPML,3) = coef0_3 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,3) + &
                coef1_3 * PML_dpotential_dzl_new(i,j,k) + coef2_3 * PML_dpotential_dzl_old(i,j,k)

        ! derivatives
        dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k)  + &
                        A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + &
                        A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) + &
                        A9 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,3)
        dpotentialdyl = A10 * PML_dpotential_dyl(i,j,k) + &
                        A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + &
                        A12 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) + &
                        A13 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,3)
        dpotentialdzl = A14 * PML_dpotential_dzl(i,j,k) + &
                        A15 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + &
                        A16 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) + &
                        A17 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,3)

        if (ispec_irreg /= 0) then !irregular element

          xixl = xix(i,j,k,ispec_irreg)
          xiyl = xiy(i,j,k,ispec_irreg)
          xizl = xiz(i,j,k,ispec_irreg)
          etaxl = etax(i,j,k,ispec_irreg)
          etayl = etay(i,j,k,ispec_irreg)
          etazl = etaz(i,j,k,ispec_irreg)
          gammaxl = gammax(i,j,k,ispec_irreg)
          gammayl = gammay(i,j,k,ispec_irreg)
          gammazl = gammaz(i,j,k,ispec_irreg)
          jacobianl = jacobian(i,j,k,ispec_irreg)
          rho_invl_jacob = jacobianl / rhostore(i,j,k,ispec)

          temp1(i,j,k) = rho_invl_jacob * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
          temp2(i,j,k) = rho_invl_jacob * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
          temp3(i,j,k) = rho_invl_jacob * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

        else !regular element

          rho_invl_jacob = jacobian_regular / rhostore(i,j,k,ispec)
          temp1(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdxl
          temp2(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdyl
          temp3(i,j,k) = rho_invl_jacob * xix_regular * dpotentialdzl

        endif

      enddo
    enddo
  enddo

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

  use specfem_par, only: NGLOB_AB,deltat
  use pml_par, only: CPML_regions,k_store_x,k_store_y,k_store_z,d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z,PML_displ_old,PML_displ_new
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

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


  call lxy_interface_parameter_computation(deltat,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                           CPML_region_local,23,A_12,A_13,A_14, &
                                           coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)
  ! displ_x
  rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) + &
                                                  coef1_1 * PML_displ_new(1,i,j,k,ispec_CPML) + &
                                                  coef2_1 * PML_displ_old(1,i,j,k,ispec_CPML)

  rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) + &
                                                  coef1_2 * PML_displ_new(1,i,j,k,ispec_CPML) + &
                                                  coef2_2 * PML_displ_old(1,i,j,k,ispec_CPML)

  displ_x = A_12 * displ(1,iglob) + A_13 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) + &
                                    A_14 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2)

  call lxy_interface_parameter_computation(deltat,kappa_x,d_x,alpha_x,kappa_z,d_z,alpha_z, &
                                           CPML_region_local,13,A_12,A_13,A_14, &
                                           coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)

  ! displ_y
  rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) + &
                                                  coef1_1 * PML_displ_new(2,i,j,k,ispec_CPML) + &
                                                  coef2_1 * PML_displ_old(2,i,j,k,ispec_CPML)

  rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) + &
                                                  coef1_2 * PML_displ_new(2,i,j,k,ispec_CPML) + &
                                                  coef2_2 * PML_displ_old(2,i,j,k,ispec_CPML)

  displ_y = A_12 * displ(2,iglob) + A_13 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) + &
                                    A_14 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2)

  call lxy_interface_parameter_computation(deltat,kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y, &
                                           CPML_region_local,12,A_12,A_13,A_14, &
                                           coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)

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

  use specfem_par, only: NGLOB_AB,deltat
  use pml_par, only: CPML_regions,k_store_x,k_store_y,k_store_z,d_store_x,d_store_y,d_store_z, &
                     alpha_store_x,alpha_store_y,alpha_store_z, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

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

  logical :: FIRST_ORDER_CONVOLUTION

  FIRST_ORDER_CONVOLUTION = .false.

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

  call l_interface_parameter_computation(deltat, kappa_y, d_y, alpha_y, kappa_z, d_z, alpha_z, &
                                         CPML_region_local, 23, A_0, A_1, A_2, A_3, A_4, &
                                         coef0_1,coef1_1,coef2_1, coef0_2,coef1_2,coef2_2)

  rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) + &
                                                    coef1_1 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_1 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,2) + &
                                                    coef1_2 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_2 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  pressure_x = A_0 * potential_dot_dot_acoustic(iglob) + A_1 * potential_dot_acoustic(iglob)  + A_2 * potential_acoustic(iglob) &
             + A_3 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,1) + A_4 * rmemory_coupling_el_ac_potential(1,i,j,k,iface,2)

  call l_interface_parameter_computation(deltat, kappa_x, d_x, alpha_x, kappa_z, d_z, alpha_z, &
                                         CPML_region_local, 13, A_0, A_1, A_2, A_3, A_4, &
                                         coef0_1,coef1_1,coef2_1, coef0_2,coef1_2,coef2_2)

  rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) + &
                                                    coef1_1 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_1 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  rmemory_coupling_el_ac_potential(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,2) + &
                                                    coef1_2 * PML_potential_acoustic_new(i,j,k,ispec_CPML) + &
                                                    coef2_2 * PML_potential_acoustic_old(i,j,k,ispec_CPML)

  pressure_y = A_0 * potential_dot_dot_acoustic(iglob) + A_1 * potential_dot_acoustic(iglob)  + A_2 * potential_acoustic(iglob) &
             + A_3 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,1) + A_4 * rmemory_coupling_el_ac_potential(2,i,j,k,iface,2)

  call l_interface_parameter_computation(deltat, kappa_x, d_x, alpha_x, kappa_y, d_y, alpha_y, &
                                         CPML_region_local, 12, A_0, A_1, A_2, A_3, A_4, &
                                         coef0_1,coef1_1,coef2_1, coef0_2,coef1_2,coef2_2)

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
subroutine lijk_parameter_computation(deltat,kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y,kappa_z,d_z,alpha_z, &
                                      CPML_region_local,index_ijk,A_0,A_1,A_2,A_3, &
                                      coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                      coef0_3,coef1_3,coef2_3)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL), intent(in) :: kappa_x,d_x,alpha_x, &
                                        kappa_y,d_y,alpha_y, &
                                        kappa_z,d_z,alpha_z

  integer, intent(in) :: CPML_region_local,index_ijk

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1,A_2,A_3
  real(kind=CUSTOM_REAL), intent(out) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2, &
                                         coef0_3,coef1_3,coef2_3

  ! local variables
  real(kind=CUSTOM_REAL) :: beta_x,beta_y,beta_z
  real(kind=CUSTOM_REAL) :: bar_A_0,bar_A_1,bar_A_2,bar_A_3

  integer :: CPML_X_ONLY_TEMP,CPML_Y_ONLY_TEMP,CPML_Z_ONLY_TEMP, &
             CPML_XY_ONLY_TEMP,CPML_XZ_ONLY_TEMP,CPML_YZ_ONLY_TEMP,CPML_XYZ_TEMP

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.

  if (index_ijk == 123) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XY_ONLY_TEMP = CPML_XY_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else if (index_ijk == 132) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_Y_ONLY
    CPML_XY_ONLY_TEMP = CPML_XZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else if (index_ijk == 231) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XY_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else
    stop 'In lijk_parameter_computation index_ijk must be equal to 123 or 132 or 231'
  endif

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y
  beta_z = alpha_z + d_z / kappa_z

  if (CPML_region_local == CPML_XYZ_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter .and. &
        abs(alpha_x-beta_z) >= min_distance_between_CPML_parameter .and. &
        abs(alpha_y-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y / kappa_z

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - alpha_z) * (alpha_x - beta_x) * (alpha_x - beta_y) / &
                ((alpha_x-alpha_y) * (alpha_x-beta_z))
      bar_A_2 = - bar_A_0 * (alpha_y - alpha_z) * (alpha_y - beta_x) * (alpha_y - beta_y) / &
                ((alpha_y-alpha_x) * (alpha_y-beta_z))
      bar_A_3 = - bar_A_0 * (beta_z - alpha_z) * (beta_z - beta_x) * (beta_z - beta_y) / &
                ((beta_z-alpha_x) * (beta_z-alpha_y))

      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
      call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)

    else
      stop 'error in lijk_parameter_computation'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY_TEMP) then

    if (abs(alpha_y-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_y / kappa_z

      !----------------A1,2,3-------------------------
      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      bar_A_1 = 0._CUSTOM_REAL
      bar_A_2 = - bar_A_0 * (alpha_y - alpha_z) * (alpha_y - beta_y) / &
                (alpha_y-beta_z)
      bar_A_3 = - bar_A_0 * (beta_z - alpha_z) * (beta_z - beta_y) / &
                (beta_z-alpha_y)

      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
      call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)
    else
      stop 'error in lijk_parameter_computation'
    endif

  else if (CPML_region_local == CPML_XZ_ONLY_TEMP) then

    if (abs(alpha_x-beta_z) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x / kappa_z

      !----------------A1,2,3-------------------------
      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      bar_A_1 = - bar_A_0 * (alpha_x - alpha_z) * (alpha_x - beta_x) / &
                (alpha_x-beta_z)
      bar_A_2 = 0._CUSTOM_REAL
      bar_A_3 = - bar_A_0 * (beta_z - alpha_z) * (beta_z - beta_x) / &
                (beta_z-alpha_x)

      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
      call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)
    else
      stop 'error in lijk_parameter_computation'
    endif

  else if (CPML_region_local == CPML_XY_ONLY_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2,3-------------------------
      beta_x = alpha_x + d_x / kappa_x
      beta_y = alpha_y + d_y / kappa_y
      beta_z = alpha_z + d_z / kappa_z

      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / &
                (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / &
                (alpha_y-alpha_x)
      bar_A_3 = 0._CUSTOM_REAL

      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
      call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)
    else
      stop 'error in lijk_parameter_computation'
    endif
  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL
    bar_A_3 = 0._CUSTOM_REAL

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
    call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)

  else if (CPML_region_local == CPML_Y_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)
    bar_A_3 = 0._CUSTOM_REAL

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
    call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = 1._CUSTOM_REAL / kappa_z

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = 0._CUSTOM_REAL
    bar_A_3 = - bar_A_0 * (beta_z - alpha_z)

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
    call compute_convolution_coef(beta_z, deltat, coef0_3, coef1_3, coef2_3)

  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2
  A_3 = bar_A_3

end subroutine lijk_parameter_computation
!
!=====================================================================
!
!
subroutine lx_parameter_computation(deltat,kappa_x,d_x,alpha_x, &
                                    CPML_region_local,A_0,A_1, &
                                    coef0_1,coef1_1,coef2_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL), intent(in) :: kappa_x,d_x,alpha_x
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1
  real(kind=CUSTOM_REAL), intent(out) :: coef0_1,coef1_1,coef2_1

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

  call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)

end subroutine lx_parameter_computation
!
!=====================================================================
!
subroutine ly_parameter_computation(deltat,kappa_y,d_y,alpha_y, &
                                    CPML_region_local,A_0,A_1, &
                                    coef0_1,coef1_1,coef2_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL), intent(in) :: kappa_y,d_y,alpha_y
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1
  real(kind=CUSTOM_REAL), intent(out) :: coef0_1,coef1_1,coef2_1

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

  call compute_convolution_coef(alpha_y, deltat, coef0_1, coef1_1, coef2_1)

end subroutine ly_parameter_computation
!
!=====================================================================
!
!
!=====================================================================
!
!
subroutine lz_parameter_computation(deltat,kappa_z,d_z,alpha_z, &
                                    CPML_region_local,A_0,A_1, &
                                    coef0_1,coef1_1,coef2_1)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL), intent(in) :: kappa_z,d_z,alpha_z
  integer, intent(in) :: CPML_region_local

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1
  real(kind=CUSTOM_REAL), intent(out) :: coef0_1,coef1_1,coef2_1

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

  call compute_convolution_coef(alpha_z, deltat, coef0_1, coef1_1, coef2_1)

end subroutine lz_parameter_computation
!
!=====================================================================
!
subroutine lxy_interface_parameter_computation(deltat,kappa_x,d_x,alpha_x,kappa_y,d_y,alpha_y, &
                                      CPML_region_local,index_ijk,A_0,A_1,A_2, &
                                      coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY, &
                       CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ
  use pml_par, only: min_distance_between_CPML_parameter

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: deltat
  real(kind=CUSTOM_REAL) :: kappa_x,d_x,alpha_x, &
                            kappa_y,d_y,alpha_y
  integer, intent(in) :: CPML_region_local,index_ijk

  real(kind=CUSTOM_REAL), intent(out) :: A_0,A_1,A_2
  real(kind=CUSTOM_REAL), intent(out) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2

  ! local variables
  real(kind=CUSTOM_REAL) :: bar_A_0,bar_A_1,bar_A_2,beta_x,beta_y

  integer :: CPML_X_ONLY_TEMP,CPML_Y_ONLY_TEMP,CPML_Z_ONLY_TEMP, &
             CPML_XY_ONLY_TEMP,CPML_XZ_ONLY_TEMP,CPML_YZ_ONLY_TEMP,CPML_XYZ_TEMP

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.

  if (index_ijk == 12) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XY_ONLY_TEMP = CPML_XY_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else if (index_ijk == 13) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Y_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_Y_ONLY
    CPML_XY_ONLY_TEMP = CPML_XZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_YZ_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else if (index_ijk == 23) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Y_ONLY_TEMP = CPML_Y_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XY_ONLY_TEMP = CPML_YZ_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
    CPML_YZ_ONLY_TEMP = CPML_XY_ONLY
    CPML_XYZ_TEMP = CPML_XYZ
  else
    stop 'In lxy_interface_parameter_computation index_ijk must be equal to 12 or 13 or 23'
  endif

  beta_x = alpha_x + d_x / kappa_x
  beta_y = alpha_y + d_y / kappa_y

  if (CPML_region_local == CPML_XYZ_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / &
                (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / &
                (alpha_y-alpha_x)


      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
    else
      stop 'error in lxy_interface_parameter_computation'
    endif

  else if (CPML_region_local == CPML_YZ_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)

  else if (CPML_region_local == CPML_XZ_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)

  else if (CPML_region_local == CPML_XY_ONLY_TEMP) then

    if (abs(alpha_x-alpha_y) >= min_distance_between_CPML_parameter) then
      !----------------A0-------------------------
      bar_A_0 = kappa_x * kappa_y

      !----------------A1,2,3-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - beta_x) * (alpha_x - beta_y) / &
                (alpha_x-alpha_y)
      bar_A_2 = - bar_A_0 * (alpha_y - beta_x) * (alpha_y - beta_y) / &
                (alpha_y-alpha_x)

      call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
      call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)
    else
      stop 'error in lxy_interface_parameter_computation'
    endif

  else if (CPML_region_local == CPML_X_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_x

    !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)

  else if (CPML_region_local == CPML_Y_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = kappa_y

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (alpha_y - beta_y)

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then

    !----------------A0-------------------------
    bar_A_0 = 1._CUSTOM_REAL

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = 0._CUSTOM_REAL

    call compute_convolution_coef(alpha_x, deltat, coef0_1, coef1_1, coef2_1)
    call compute_convolution_coef(alpha_y, deltat, coef0_2, coef1_2, coef2_2)

  else
    stop 'error in lxy_interface_parameter_computation'
  endif

  A_0 = bar_A_0
  A_1 = bar_A_1
  A_2 = bar_A_2

end subroutine lxy_interface_parameter_computation

