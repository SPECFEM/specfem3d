!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: it,deltat, &
                         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
                         kappastore,mustore
  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z,&
                     d_store_x,d_store_y,d_store_z,alpha_store, &
                     PML_dux_dxl, PML_dux_dyl, PML_dux_dzl, PML_duy_dxl, PML_duy_dyl, PML_duy_dzl, &
                     PML_duz_dxl, PML_duz_dyl, PML_duz_dzl, &
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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,2) ::  &
                          rmemory_dux_dxl_x, rmemory_dux_dyl_x, rmemory_dux_dzl_x, &
                          rmemory_duy_dxl_y, rmemory_duy_dyl_y, rmemory_duy_dzl_y, &
                          rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duz_dzl_z  

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML) ::  &
                          rmemory_duy_dyl_x, rmemory_duz_dzl_x, rmemory_duz_dxl_x, rmemory_duy_dxl_x, &
                          rmemory_dux_dxl_y, rmemory_duz_dzl_y, rmemory_duz_dyl_y, rmemory_dux_dyl_y, &
                          rmemory_dux_dxl_z, rmemory_duy_dyl_z, rmemory_duy_dzl_z, rmemory_dux_dzl_z
  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul,kappal
  real(kind=CUSTOM_REAL) :: duxdxl_x,duxdyl_x,duxdzl_x,duydxl_x,duydyl_x,duzdxl_x,duzdzl_x
  real(kind=CUSTOM_REAL) :: duxdxl_y,duxdyl_y,duydxl_y,duydyl_y,duydzl_y,duzdyl_y,duzdzl_y
  real(kind=CUSTOM_REAL) :: duxdxl_z,duxdzl_z,duydyl_z,duydzl_z,duzdxl_z,duzdyl_z,duzdzl_z
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17 ! for convolution of strain(complex)
  real(kind=CUSTOM_REAL) :: A18,A19,A20 ! for convolution of strain(simple)

  do k=1,NGLLZ
     do j=1,NGLLY
         do i=1,NGLLX
            kappal = kappastore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2.0d0*mul
            xixl = xix(i,j,k,ispec)
            xiyl = xiy(i,j,k,ispec)
            xizl = xiz(i,j,k,ispec)
            etaxl = etax(i,j,k,ispec)
            etayl = etay(i,j,k,ispec)
            etazl = etaz(i,j,k,ispec)
            gammaxl = gammax(i,j,k,ispec)
            gammayl = gammay(i,j,k,ispec)
            gammazl = gammaz(i,j,k,ispec)
            jacobianl = jacobian(i,j,k,ispec)

            if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- X-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = 1.d0 / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = - d_store_x(i,j,k,ispec_CPML) / (k_store_x(i,j,k,ispec_CPML)**2)

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dxl_new(i,j,k) * coef1_2 + PML_dux_dxl(i,j,k) * coef2_2

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dxl_new(i,j,k) * coef1_2 + PML_duy_dxl(i,j,k) * coef2_2

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dxl_new(i,j,k) * coef1_2 + PML_duz_dxl(i,j,k) * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9  = k_store_x(i,j,k,ispec_CPML)
              A10 = d_store_x(i,j,k,ispec_CPML)
              A11 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML)
              A13 = d_store_x(i,j,k,ispec_CPML)
              A14 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A15 and A16 --------------------------
              A15 = k_store_x(i,j,k,ispec_CPML)
              A16 = d_store_x(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1

              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1

              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              !---------------------- A17 and A18 --------------------------
              A17 = 1.d0
              A18 = 0.0

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A19 and A20 --------------------------
              A19 = 1.d0
              A20 = 0.0

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = 0.d0

            else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Y-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML)
              A7 = d_store_y(i,j,k,ispec_CPML)
              A8 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A9, A10 and A11 --------------------------
              A9 = 1.d0/k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = - d_store_y(i,j,k,ispec_CPML) / (k_store_y(i,j,k,ispec_CPML) ** 2)

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dyl_new(i,j,k) * coef1_2 + PML_dux_dyl(i,j,k) * coef2_2

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dyl_new(i,j,k) * coef1_2 + PML_duy_dyl(i,j,k) * coef2_2

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dyl_new(i,j,k) * coef1_2 + PML_duz_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_y(i,j,k,ispec_CPML)
              A13 = d_store_y(i,j,k,ispec_CPML)
              A14 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A15 and A16 --------------------------
              A15 = 1.d0
              A16 = 0.d0

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A17 and A18 --------------------------
              A17 = k_store_y(i,j,k,ispec_CPML)
              A18 = d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1

              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1

              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

              !---------------------- A19 and A20--------------------------
              A19 = 1.d0
              A20 = 0.0

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = 0.d0

            else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- Z-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML)
              A7 = d_store_z(i,j,k,ispec_CPML)
              A8 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_z(i,j,k,ispec_CPML)
              A10 = d_store_z(i,j,k,ispec_CPML)
              A11 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A12, A13 and A14 --------------------------
              A12 = 1.0 / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = - d_store_z(i,j,k,ispec_CPML) / (k_store_z(i,j,k,ispec_CPML) ** 2)

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dzl_new(i,j,k) * coef1_2 + PML_dux_dzl(i,j,k) * coef2_2

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dzl_new(i,j,k) * coef1_2 + PML_duy_dzl(i,j,k) * coef2_2

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dzl_new(i,j,k) * coef1_2 + PML_duz_dzl(i,j,k) * coef2_2

              !---------------------- A15 and A16 --------------------------
              A15 = 1.d0
              A16 = 0.d0

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A17 and A18 --------------------------
              A17 = 1.d0
              A18 = 0.d0

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A19 and A20 --------------------------
              A19 = k_store_z(i,j,k,ispec_CPML)
              A20 = d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1

              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1

              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

            else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- XY-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / k_store_x(i,j,k,ispec_CPML)**2

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dxl_new(i,j,k) * coef1_2 + PML_dux_dxl(i,j,k) * coef2_2

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dxl_new(i,j,k) * coef1_2 + PML_duy_dxl(i,j,k) * coef2_2

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dxl_new(i,j,k) * coef1_2 + PML_duz_dxl(i,j,k) * coef2_2


              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) - &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / k_store_y(i,j,k,ispec_CPML)**2

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dyl_new(i,j,k) * coef1_2 + PML_dux_dyl(i,j,k) * coef2_2

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dyl_new(i,j,k) * coef1_2 + PML_duy_dyl(i,j,k) * coef2_2

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dyl_new(i,j,k) * coef1_2 + PML_duz_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                    + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                    + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dzl(i,j,k) * it*deltat * coef2_2

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dzl(i,j,k) * it*deltat * coef2_2

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dzl(i,j,k) * it*deltat * coef2_2

              !---------------------- A15 and A16 --------------------------
              A15 = k_store_x(i,j,k,ispec_CPML)
              A16 = d_store_x(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1

              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1

              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              !---------------------- A17 and A18 --------------------------
              A17 = k_store_y(i,j,k,ispec_CPML)
              A18 = d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1

              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1

              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

              !---------------------- A19 and A20--------------------------
              A19 = 1.d0
              A20 = 0.0

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = 0.d0

            else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- XZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = ( d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / k_store_x(i,j,k,ispec_CPML)**2

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dxl_new(i,j,k) * coef1_2 + PML_dux_dxl(i,j,k) * coef2_2

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dxl_new(i,j,k) * coef1_2 + PML_duy_dxl(i,j,k) * coef2_2

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dxl_new(i,j,k) * coef1_2 + PML_duz_dxl(i,j,k) * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A10 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                   + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dyl(i,j,k) * it*deltat * coef2_2

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dyl(i,j,k) * it*deltat * coef2_2

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dyl(i,j,k) * it*deltat * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   - d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / k_store_z(i,j,k,ispec_CPML)**2

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dzl_new(i,j,k) * coef1_2 + PML_dux_dzl(i,j,k) * coef2_2

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dzl_new(i,j,k) * coef1_2 + PML_duy_dzl(i,j,k) * coef2_2

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dzl_new(i,j,k) * coef1_2 + PML_duz_dzl(i,j,k) * coef2_2

              !---------------------- A15 and A16 --------------------------
              A15 = k_store_x(i,j,k,ispec_CPML)
              A16 = d_store_x(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1

              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1

              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              !---------------------- A17 and A18 --------------------------
              A17 = 1.0d0
              A18 = 0.d0

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = 0.d0
              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A19 and A20 --------------------------
              A19 = k_store_z(i,j,k,ispec_CPML)
              A20 = d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1

              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1

              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

            else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- YZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A7 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   + it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1
              rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dxl(i,j,k) * it*deltat * coef2_2

              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dxl(i,j,k) * it*deltat * coef2_2

              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dxl(i,j,k) * it*deltat * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) -&
                   d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / k_store_y(i,j,k,ispec_CPML)**2

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dyl_new(i,j,k) * coef1_2 + PML_dux_dyl(i,j,k) * coef2_2

              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dyl_new(i,j,k) * coef1_2 + PML_duy_dyl(i,j,k) * coef2_2

              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dyl_new(i,j,k) * coef1_2 + PML_duz_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) -&
                   d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / k_store_z(i,j,k,ispec_CPML)**2

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                   + PML_dux_dzl_new(i,j,k) * coef1_2 + PML_dux_dzl(i,j,k) * coef2_2

              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                   + PML_duy_dzl_new(i,j,k) * coef1_2 + PML_duy_dzl(i,j,k) * coef2_2

              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                   + PML_duz_dzl_new(i,j,k) * coef1_2 + PML_duz_dzl(i,j,k) * coef2_2

              !---------------------- A15 and A16 --------------------------
              A15 = 1.0d0
              A16 = 0.0d0

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = 0.d0

              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = 0.d0

              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = 0.d0

              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = 0.d0

              !---------------------- A17 and A18 --------------------------
              A17 = k_store_y(i,j,k,ispec_CPML)
              A18 = d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1

              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1

              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

              !---------------------- A19 and A20--------------------------
              A19 = k_store_z(i,j,k,ispec_CPML)
              A20 = d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1

              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1

              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

            else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then

              !------------------------------------------------------------------------------
              !---------------------------- XYZ-corner C-PML --------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              if( abs(d_store_x(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A7 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_x(i,j,k,ispec_CPML)
                 A8 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                      d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) * &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) - &
                      d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)**2)
              else
                 A7 = ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
                      d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                      k_store_x(i,j,k,ispec_CPML) + &
                      it*deltat * d_store_y(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_x(i,j,k,ispec_CPML)
                 A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1
              rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1
              rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1

              if( abs(d_store_x(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dxl_new(i,j,k) * coef1_2 + PML_dux_dxl(i,j,k) * coef2_2
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dxl_new(i,j,k) * coef1_2 + PML_duy_dxl(i,j,k) * coef2_2
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dxl_new(i,j,k) * coef1_2 + PML_duz_dxl(i,j,k) * coef2_2
              else
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dxl(i,j,k) * it*deltat * coef2_2
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dxl(i,j,k) * it*deltat * coef2_2
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dxl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dxl(i,j,k) * it*deltat * coef2_2
              endif

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              if( abs(d_store_y(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A10 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_y(i,j,k,ispec_CPML)
                 A11 = ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)**2)
              else
                 A10 = ( d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                       + d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                       k_store_y(i,j,k,ispec_CPML) + &
                       it*deltat * d_store_x(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_y(i,j,k,ispec_CPML)
                 A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1
              rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1
              rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1

              if( abs(d_store_y(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dyl_new(i,j,k) * coef1_2 + PML_dux_dyl(i,j,k) * coef2_2
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dyl_new(i,j,k) * coef1_2 + PML_duy_dyl(i,j,k) * coef2_2
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dyl_new(i,j,k) * coef1_2 + PML_duz_dyl(i,j,k) * coef2_2
              else
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dyl(i,j,k) * it*deltat * coef2_2
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dyl(i,j,k) * it*deltat * coef2_2
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dyl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dyl(i,j,k) * it*deltat * coef2_2
              endif

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              if( abs(d_store_z(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A13 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)/d_store_z(i,j,k,ispec_CPML)
                 A14 = ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)**2)
              else
                 A13 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                       + d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) ) / &
                       k_store_z(i,j,k,ispec_CPML) + &
                       it*deltat * d_store_x(i,j,k,ispec_CPML)*d_store_y(i,j,k,ispec_CPML)/k_store_z(i,j,k,ispec_CPML)
                 A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)

              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1
              rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1
              rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              if( abs(d_store_z(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dzl_new(i,j,k) * coef1_2 + PML_dux_dzl(i,j,k) * coef2_2
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dzl_new(i,j,k) * coef1_2 + PML_duy_dzl(i,j,k) * coef2_2
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dzl_new(i,j,k) * coef1_2 + PML_duz_dzl(i,j,k) * coef2_2
              else
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &
                      + PML_dux_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_dux_dzl(i,j,k) * it*deltat * coef2_2
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &
                      + PML_duy_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_duy_dzl(i,j,k) * it*deltat * coef2_2
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &
                      + PML_duz_dzl_new(i,j,k) * it*deltat * coef1_2 + PML_duz_dzl(i,j,k) * it*deltat * coef2_2
              endif

              !---------------------- A15 and A16 --------------------------
              A15 = k_store_x(i,j,k,ispec_CPML)
              A16 = d_store_x(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML) &
                   + PML_duz_dyl_new(i,j,k) * coef1_1 + PML_duz_dyl(i,j,k) * coef2_1

              rmemory_duy_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dzl_new(i,j,k) * coef1_1 + PML_duy_dzl(i,j,k) * coef2_1

              rmemory_duy_dyl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              !---------------------- A17 and A18 --------------------------
              A17 = k_store_y(i,j,k,ispec_CPML)
              A18 = d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duz_dzl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dzl_new(i,j,k) * coef1_1 + PML_duz_dzl(i,j,k) * coef2_1

              rmemory_duz_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duz_dxl_new(i,j,k) * coef1_1 + PML_duz_dxl(i,j,k) * coef2_1

              rmemory_dux_dzl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dzl_new(i,j,k) * coef1_1 + PML_dux_dzl(i,j,k) * coef2_1

              rmemory_dux_dxl_z(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

              !---------------------- A19 and A20 --------------------------
              A19 = k_store_z(i,j,k,ispec_CPML)
              A20 = d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0))/ bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_duy_dyl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dyl_new(i,j,k) * coef1_1 + PML_duy_dyl(i,j,k) * coef2_1

              rmemory_duy_dxl_x(i,j,k,ispec_CPML) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML) &
                   + PML_duy_dxl_new(i,j,k) * coef1_1 + PML_duy_dxl(i,j,k) * coef2_1

              rmemory_dux_dyl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dyl_new(i,j,k) * coef1_1 + PML_dux_dyl(i,j,k) * coef2_1

              rmemory_dux_dxl_y(i,j,k,ispec_CPML) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML) &
                   + PML_dux_dxl_new(i,j,k) * coef1_1 + PML_dux_dxl(i,j,k) * coef2_1

            else
              stop 'wrong PML flag in PML memory variable calculation routine'
            endif

            duxdxl_x = A6 * PML_dux_dxl(i,j,k)  &
                 + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
            duxdyl_x = A9 * PML_dux_dyl(i,j,k)  &
                 + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
            duxdzl_x = A12 * PML_dux_dzl(i,j,k)  &
                 + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
            duzdzl_x = A17 * PML_duz_dzl(i,j,k) + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML)
            duzdxl_x = A17 * PML_duz_dxl(i,j,k) + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML)
            duydyl_x = A19 * PML_duy_dyl(i,j,k) + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML)
            duydxl_x = A19 * PML_duy_dxl(i,j,k) + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML)

            duydxl_y = A6 * PML_duy_dxl(i,j,k)  &
                 + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2)
            duydyl_y = A9 * PML_duy_dyl(i,j,k)  &
                 + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
            duydzl_y = A12 * PML_duy_dzl(i,j,k)  &
                 + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
            duzdzl_y = A15 * PML_duz_dzl(i,j,k) + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML)
            duzdyl_y = A15 * PML_duz_dyl(i,j,k) + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML)
            duxdyl_y = A19 * PML_dux_dyl(i,j,k) + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML)
            duxdxl_y = A19 * PML_dux_dxl(i,j,k) + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML)

            duzdxl_z = A6 * PML_duz_dxl(i,j,k)  &
                 + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2)
            duzdyl_z = A9 * PML_duz_dyl(i,j,k)  &
                 + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2)
            duzdzl_z = A12 * PML_duz_dzl(i,j,k)  &
                 + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)
            duydzl_z = A15 * PML_duy_dzl(i,j,k) + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML)
            duydyl_z = A15 * PML_duy_dyl(i,j,k) + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML)
            duxdzl_z = A17 * PML_dux_dzl(i,j,k) + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML)
            duxdxl_z = A17 * PML_dux_dxl(i,j,k) + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML)

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

            ! form dot product with test vector, non-symmetric form
            tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
            tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
            tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
            tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
            tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

            tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
            tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
            tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

          enddo
      enddo
  enddo

end subroutine pml_compute_memory_variables_elastic

!
!=====================================================================
!
!

subroutine pml_compute_memory_variables_acoustic(ispec,ispec_CPML,temp1,temp2,temp3,&
                                                 rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&
                         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian,&
                         it,deltat,rhostore
  use pml_par, only: NSPEC_CPML,CPML_regions,k_store_x,k_store_y,k_store_z,&
                     d_store_x,d_store_y,d_store_z,alpha_store,&
                     PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl,&
                     PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,FOUR_THIRDS, &
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,2) :: &
                               rmemory_dpotential_dxl, rmemory_dpotential_dyl, rmemory_dpotential_dzl

  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: rho_invl_jacob,rhoin_jacob_jk,rhoin_jacob_ik,rhoin_jacob_ij
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9,A10,A11,A12,A13,A14

  do k=1,NGLLZ
     do j=1,NGLLY
         do i=1,NGLLX
            xixl = xix(i,j,k,ispec)
            xiyl = xiy(i,j,k,ispec)
            xizl = xiz(i,j,k,ispec)
            etaxl = etax(i,j,k,ispec)
            etayl = etay(i,j,k,ispec)
            etazl = etaz(i,j,k,ispec)
            gammaxl = gammax(i,j,k,ispec)
            gammayl = gammay(i,j,k,ispec)
            gammazl = gammaz(i,j,k,ispec)
            jacobianl = jacobian(i,j,k,ispec)
            rho_invl_jacob = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec) * jacobianl
            rhoin_jacob_jk = rho_invl_jacob * wgllwgll_yz(j,k)
            rhoin_jacob_ik = rho_invl_jacob * wgllwgll_xz(i,k)
            rhoin_jacob_ij = rho_invl_jacob * wgllwgll_xy(i,j)

            if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- X-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = 1.d0 / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = - d_store_x(i,j,k,ispec_CPML) / (k_store_x(i,j,k,ispec_CPML)**2)

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_2 + PML_dpotential_dxl(i,j,k) * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9  = k_store_x(i,j,k,ispec_CPML)
              A10 = d_store_x(i,j,k,ispec_CPML)
              A11 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_1 + PML_dpotential_dyl(i,j,k) * coef2_1
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML)
              A13 = d_store_x(i,j,k,ispec_CPML)
              A14 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_1 + PML_dpotential_dzl(i,j,k) * coef2_1
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = 0.d0

            else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
              !------------------------------------------------------------------------------
              !---------------------------- Y-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML)
              A7 = d_store_y(i,j,k,ispec_CPML)
              A8 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_1 + PML_dpotential_dxl(i,j,k) * coef2_1
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A9, A10 and A11 --------------------------
              A9 = 1.d0/k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = - d_store_y(i,j,k,ispec_CPML) / (k_store_y(i,j,k,ispec_CPML) ** 2)

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_2 + PML_dpotential_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_y(i,j,k,ispec_CPML)
              A13 = d_store_y(i,j,k,ispec_CPML)
              A14 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_1 + PML_dpotential_dzl(i,j,k) * coef2_1
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = 0.d0

            else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- Z-surface C-PML ---------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML)
              A7 = d_store_z(i,j,k,ispec_CPML)
              A8 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_1 + PML_dpotential_dxl(i,j,k) * coef2_1
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_z(i,j,k,ispec_CPML)
              A10 = d_store_z(i,j,k,ispec_CPML)
              A11 = 0.d0

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_1 + PML_dpotential_dyl(i,j,k) * coef2_1
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = 0.d0

              !---------------------- A12, A13 and A14 --------------------------
              A12 = 1.0 / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = - d_store_z(i,j,k,ispec_CPML) / (k_store_z(i,j,k,ispec_CPML) ** 2)

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_2 + PML_dpotential_dzl(i,j,k) * coef2_2

            else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- XY-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / k_store_x(i,j,k,ispec_CPML)**2

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_2 + PML_dpotential_dxl(i,j,k) * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) - &
                   d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / k_store_y(i,j,k,ispec_CPML)**2

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_2 + PML_dpotential_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                    + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                    + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_1 + PML_dpotential_dzl(i,j,k) * coef2_1
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dzl_new(i,j,k) * it*deltat * coef1_2 &
                   + PML_dpotential_dzl(i,j,k) * it*deltat * coef2_2

            else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- XZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              A7 = 0.d0
              A8 = ( d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                   d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / k_store_x(i,j,k,ispec_CPML)**2

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_2 + PML_dpotential_dxl(i,j,k) * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              A10 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                    + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                    + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_1 + PML_dpotential_dyl(i,j,k) * coef2_1
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dyl_new(i,j,k) * it*deltat * coef1_2 &
                   + PML_dpotential_dyl(i,j,k) * it*deltat * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   - d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / k_store_z(i,j,k,ispec_CPML)**2

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_2 + PML_dpotential_dzl(i,j,k) * coef2_2

            else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then

              !------------------------------------------------------------------------------
              !---------------------------- YZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              A7 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                   + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                   + it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_1 + PML_dpotential_dxl(i,j,k) * coef2_1
              rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dxl_new(i,j,k) * it*deltat * coef1_2 &
                   + PML_dpotential_dxl(i,j,k) * it*deltat * coef2_2

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              A10 = 0.d0
              A11 = ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) -&
                   d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / k_store_y(i,j,k,ispec_CPML)**2

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_2 + PML_dpotential_dyl(i,j,k) * coef2_2

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              A13 = 0.d0
              A14 = ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) -&
                   d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / k_store_z(i,j,k,ispec_CPML)**2

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) / bb
                 coef2_2 = ( 1.d0 - exp(-bb * deltat/2.d0) ) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
              rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_2 + PML_dpotential_dzl(i,j,k) * coef2_2

            else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then

              !------------------------------------------------------------------------------
              !---------------------------- XYZ-corner C-PML --------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              if( abs(d_store_x(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A7 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_x(i,j,k,ispec_CPML)
                 A8 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                      d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) * &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) - &
                      d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)**2)
              else
                 A7 = (d_store_z(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML)+ &
                      d_store_y(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML)) / &
                      k_store_x(i,j,k,ispec_CPML) + &
                      it*deltat * d_store_y(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_x(i,j,k,ispec_CPML)
                 A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dxl_new(i,j,k) * coef1_1 + PML_dpotential_dxl(i,j,k) * coef2_1

              if(abs(d_store_x(i,j,k,ispec_CPML))> 1.d-5)then
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k) * coef1_2 + PML_dpotential_dxl(i,j,k) * coef2_2
              else
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k) * it*deltat * coef1_2 &
                      + PML_dpotential_dxl(i,j,k) * it*deltat * coef2_2
              endif

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              if( abs(d_store_y(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A10 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_y(i,j,k,ispec_CPML)
                 A11 = ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)**2)
              else
                 A10 = (d_store_z(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML) &
                       +d_store_x(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML)) / &
                       k_store_y(i,j,k,ispec_CPML) + &
                       it*deltat * d_store_x(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_y(i,j,k,ispec_CPML)
                 A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dyl_new(i,j,k) * coef1_1 + PML_dpotential_dyl(i,j,k) * coef2_1

              if(abs(d_store_y(i,j,k,ispec_CPML))> 1.d-5)then
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k) * coef1_2 + PML_dpotential_dyl(i,j,k) * coef2_2
              else
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k) * it*deltat * coef1_2 &
                      + PML_dpotential_dyl(i,j,k) * it*deltat * coef2_2
              endif

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              if( abs(d_store_z(i,j,k,ispec_CPML)) > 1.d-5 ) then
                 A13 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)/d_store_z(i,j,k,ispec_CPML)
                 A14 = ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)**2)
              else
                 A13 = (d_store_y(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML)&
                       +d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML)) / &
                       k_store_z(i,j,k,ispec_CPML) + &
                       it*deltat * d_store_x(i,j,k,ispec_CPML)*d_store_y(i,j,k,ispec_CPML)/k_store_z(i,j,k,ispec_CPML)
                 A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)
              endif

              bb = alpha_store(i,j,k,ispec_CPML)
              coef0_1 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0
              endif

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)
              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              endif

              rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                   + PML_dpotential_dzl_new(i,j,k) * coef1_1 + PML_dpotential_dzl(i,j,k) * coef2_1

              if(abs(d_store_z(i,j,k,ispec_CPML))> 1.d-5)then
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k) * coef1_2 + PML_dpotential_dzl(i,j,k) * coef2_2
              else
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k) * it*deltat * coef1_2 &
                      + PML_dpotential_dzl(i,j,k) * it*deltat * coef2_2
              endif

            else
              stop 'wrong PML flag in PML memory variable calculation routine'
            endif

            dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k)  &
                 + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
            dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k)  &
                 + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
            dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k)  &
                 + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
            temp1(i,j,k) = rhoin_jacob_jk * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
            temp2(i,j,k) = rhoin_jacob_ik * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
            temp3(i,j,k) = rhoin_jacob_ij * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

          enddo
      enddo
  enddo

end subroutine pml_compute_memory_variables_acoustic

!
!=====================================================================
!
!

subroutine pml_compute_memory_variables_acoustic_elastic(ispec_CPML,iface,iglob,i,j,k,&
                                                         displ_x,displ_y,displ_z,displ,veloc,&
                                                         num_coupling_ac_el_faces,rmemory_coupling_ac_el_displ)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat
  use pml_par,only : CPML_regions,k_store_x,k_store_y,k_store_z,d_store_x,d_store_y,d_store_z,alpha_store
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,&
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec_CPML,iface,iglob,num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB), intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2) :: & 
                                    rmemory_coupling_ac_el_displ

  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9,A10,A11,A12,A13,A14


  if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- X-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------

    ! displ_x
    A6 = 1.d0
    A7 = 0.d0
    A8 = 0.d0

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = 0.d0
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = 0.d0

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2)  

    ! displ_y
    A9 = k_store_x(i,j,k,ispec_CPML)
    A10 = d_store_x(i,j,k,ispec_CPML)
    A11 = 0.d0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = 0.d0

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML)
    A14 = 0.d0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = 0.d0

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
    !------------------------------------------------------------------------------
    !---------------------------- Y-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------

    ! displ_x
    A6 = k_store_y(i,j,k,ispec_CPML)
    A7 = d_store_y(i,j,k,ispec_CPML)
    A8 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = 0.d0

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2)  

    ! displ_y
    A9 = 1.d0
    A10 = 0.d0
    A11 = 0.d0

    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = 0.d0
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = 0.d0

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_y(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                      + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = 0.d0

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- Z-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------

    ! displ_x
    A6 = k_store_z(i,j,k,ispec_CPML)
    A7 = d_store_z(i,j,k,ispec_CPML)
    A8 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = 0.d0

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2)  

    ! displ_y
    A9 = k_store_z(i,j,k,ispec_CPML)
    A10 = d_store_z(i,j,k,ispec_CPML)
    A11 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = 0.d0

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = 1.d0
    A13 = 0.d0
    A14 = 0.d0

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = 0.d0
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = 0.d0

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)


  else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- XY-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------
    ! displ_x
    A6 = k_store_y(i,j,k,ispec_CPML)
    A7 = d_store_y(i,j,k,ispec_CPML)
    A8 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) +  deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = 0.d0

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) 

    ! displ_y
    A9 = k_store_x(i,j,k,ispec_CPML)
    A10 = d_store_x(i,j,k,ispec_CPML)
    A11 = 0.d0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = 0.d0

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
          + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
          + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
    A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
         + (displ(3,iglob)) * it*deltat * coef2_2

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)


  else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- XZ-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------

    ! displ_x
    A6 = k_store_z(i,j,k,ispec_CPML)
    A7 = d_store_z(i,j,k,ispec_CPML)
    A8 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = 0.d0

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) 

    ! displ_y
    A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
    A10 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
          + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
          + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
    A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1


    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                      + (displ(2,iglob)) * it*deltat * coef2_2

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = 0.d0

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)


  else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- YZ-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------

    ! displ_x
    A6 = k_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A7 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
         + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
         + it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
    A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_2 + (displ(1,iglob)) * coef2_2

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) 

    ! displ_y
    A9 = k_store_z(i,j,k,ispec_CPML)
    A10 = d_store_z(i,j,k,ispec_CPML)
    A11 = 0.d0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = 0.d0

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_y(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = 0.d0

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)


  else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then

    !------------------------------------------------------------------------------
    !---------------------------- XYZ-corner C-PML --------------------------------
    !------------------------------------------------------------------------------
    ! displ_x
    A6 = k_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A7 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
         + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
         + it*deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
    A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_1 + (displ(1,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) &
                      + (displ(1,iglob) + deltat * veloc(1,iglob)) * coef1_2 + (displ(1,iglob)) * coef2_2

    displ_x = A6 * displ(1,iglob) + A7 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,1) &
                                  + A8 * rmemory_coupling_ac_el_displ(1,i,j,k,iface,2) 

    ! displ_y
    A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
    A10 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
         + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
         + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
    A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1


    rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * coef1_1 + (displ(2,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) &
                      + (displ(2,iglob) + deltat * veloc(2,iglob)) * it*deltat * coef1_2 &
                      + (displ(2,iglob)) * it*deltat * coef2_2

    displ_y = A9 * displ(2,iglob) + A10 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,1) &
                                  + A11 * rmemory_coupling_ac_el_displ(2,i,j,k,iface,2) 

    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
          + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
          + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
    A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) = coef0_1 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * coef1_1 + (displ(3,iglob)) * coef2_1
    rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) = coef0_2 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2) &
         + (displ(3,iglob) + deltat * veloc(3,iglob)) * it*deltat * coef1_2 &
         + (displ(3,iglob)) * it*deltat * coef2_2

    displ_z = A12 * displ(3,iglob) + A13 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,1) &
                                   + A14 * rmemory_coupling_ac_el_displ(3,i,j,k,iface,2)
  else
    stop 'wrong PML flag in PML memory variable calculation routine'
  endif

end subroutine pml_compute_memory_variables_acoustic_elastic

!
!=====================================================================
!
subroutine pml_compute_memory_variables_elastic_acoustic(ispec_CPML,iface,iglob,i,j,k,&
                                                         pressure,potential_acoustic,potential_dot_acoustic,&
                                                         num_coupling_ac_el_faces,rmemory_coupling_el_ac_potential)
  ! calculates C-PML elastic memory variables and computes stress sigma

  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use specfem_par, only: NGLOB_AB,it,deltat
  use pml_par,only : CPML_regions,k_store_x,k_store_y,d_store_x,d_store_y,alpha_store
  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,&
                       CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ

  implicit none

  integer, intent(in) :: ispec_CPML,iface,iglob,num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: pressure
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB), intent(in) :: potential_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2) :: & 
                                    rmemory_coupling_el_ac_potential

  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A12,A13,A14


  if( CPML_regions(ispec_CPML) == CPML_X_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- X-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------

    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML)
    A14 = 0.d0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                      + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
                      + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = 0.d0

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_Y_ONLY ) then
    !------------------------------------------------------------------------------
    !---------------------------- Y-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------
    ! displ_z
    A12 = k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_y(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                      + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
                      + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = 0.d0

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2)  


  else if( CPML_regions(ispec_CPML) == CPML_Z_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- Z-surface C-PML ---------------------------------
    !------------------------------------------------------------------------------
    ! displ_z
    A12 = 1.d0
    A13 = 0.d0
    A14 = 0.d0

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = 0.d0
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = 0.d0

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_XY_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- XY-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------
    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
          + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
          + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
    A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
         + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
         + (potential_acoustic(iglob)) * it*deltat * coef2_2

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_XZ_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- XZ-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------
    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)

    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
         + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = 0.d0

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_YZ_ONLY ) then

    !------------------------------------------------------------------------------
    !---------------------------- YZ-edge C-PML -----------------------------------
    !------------------------------------------------------------------------------

    ! displ_z
    A12 = k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_y(i,j,k,ispec_CPML)
    A14 = 0.0

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
         + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = 0.d0

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 


  else if( CPML_regions(ispec_CPML) == CPML_XYZ ) then

    !------------------------------------------------------------------------------
    !---------------------------- XYZ-corner C-PML --------------------------------
    !------------------------------------------------------------------------------
    ! displ_z
    A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
    A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
          + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
          + it*deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
    A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)

    bb = alpha_store(i,j,k,ispec_CPML)
    coef0_1 = exp(-bb * deltat)

    if( abs(bb) > 1.d-5 ) then
       coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
       coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
    else
       coef1_1 = deltat/2.0d0
       coef2_1 = deltat/2.0d0
    endif

    coef0_2 = coef0_1
    coef1_2 = coef1_1
    coef2_2 = coef2_1

    rmemory_coupling_el_ac_potential(i,j,k,iface,1) = coef0_1 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * coef1_1 &
         + (potential_acoustic(iglob)) * coef2_1
    rmemory_coupling_el_ac_potential(i,j,k,iface,2) = coef0_2 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) &
         + (potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob)) * it*deltat * coef1_2 &
         + (potential_acoustic(iglob)) * it*deltat * coef2_2

    pressure = A12 * potential_acoustic(iglob) + A13 * rmemory_coupling_el_ac_potential(i,j,k,iface,1) &
                                               + A14 * rmemory_coupling_el_ac_potential(i,j,k,iface,2) 
  else
    stop 'wrong PML flag in PML memory variable calculation routine'
  endif

end subroutine pml_compute_memory_variables_elastic_acoustic

