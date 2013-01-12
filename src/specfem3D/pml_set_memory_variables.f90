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

subroutine pml_set_memory_variables(ispec,ispec_CPML,deltat,jacobianl,tempx1,tempy1,tempz1,tempx2,tempy2,tempz2, &
                                    tempx3,tempy3,tempz3,sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz, &
                                    sigma_yx,sigma_zx,sigma_zy,lambdal,mul,lambdalplus2mul,xixl,xiyl,xizl, &
                                    etaxl,etayl,etazl,gammaxl,gammayl,gammazl)

  ! calculates C-PML elastic memory variables and computes stress sigma

  use specfem_par, only: it
  use specfem_par_elastic, only: ispec_is_elastic
  use specfem_par_acoustic, only: ispec_is_acoustic
  use pml_par
  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer, intent(in) :: ispec,ispec_CPML

  real(kind=CUSTOM_REAL), intent(in) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL), intent(in) :: deltat,xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: tempz1,tempz2,tempz3

  ! local parameters
  integer :: i,j,k

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: duxdxl_x,duxdyl_x,duxdzl_x,duydxl_x,duydyl_x,duzdxl_x,duzdzl_x
  real(kind=CUSTOM_REAL) :: duxdxl_y,duxdyl_y,duydxl_y,duydyl_y,duydzl_y,duzdyl_y,duzdzl_y
  real(kind=CUSTOM_REAL) :: duxdxl_z,duxdzl_z,duydyl_z,duydzl_z,duzdxl_z,duzdyl_z,duzdzl_z
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: bb,coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  real(kind=CUSTOM_REAL) :: A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17 ! for convolution of strain(complex)
  real(kind=CUSTOM_REAL) :: A18,A19,A20 ! for convolution of strain(simple)

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           if( CPML_regions(ispec_CPML) == 1 ) then
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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_2

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2)
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_2
                 
                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) 
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = 0.d0

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) 
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) 
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = 0.d0

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
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
                 end if

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2)

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2)
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2)


                 !---------------------- A17 and A18 --------------------------
                 A17 = 1.d0
                 A18 = 0.0

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2)
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2)  

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) 
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) 


                 !---------------------- A19 and A20 --------------------------
                 A19 = 1.d0
                 A20 = 0.0

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) 
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) 
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2)
              endif

           elseif( CPML_regions(ispec_CPML) == 2 ) then
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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) 
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) 
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                   rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) & 
                   + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_1
                   rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = 0.d0

                   dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                        + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_2

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) & 
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = 0.d0

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
                 !---------------------- A15 and A16 --------------------------
                 A15 = 1.d0
                 A16 = 0.d0

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) 

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) 
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2)
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2)
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) 


                 !---------------------- A19 and A20--------------------------
                 A19 = 1.d0
                 A20 = 0.0

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) 
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2)
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) 
              endif

           elseif( CPML_regions(ispec_CPML) == 3 ) then
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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) 
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2)
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = 0.d0

                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) 
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = 0.d0

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_2

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) 
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
                 !---------------------- A15 and A16 --------------------------
                 A15 = 1.d0
                 A16 = 0.d0

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2)
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2)

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) 
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2)

                 !---------------------- A17 and A18 --------------------------
                 A17 = 1.d0
                 A18 = 0.d0

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2)
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2)
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) 

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
                 end if

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2)
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2)

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2)
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2)
              endif

           elseif( CPML_regions(ispec_CPML) == 4 ) then
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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_2

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) 
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
              endif

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
              end if


              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_2

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
              endif

              !---------------------- A12, A13 and A14 --------------------------
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A13 = d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      + d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                      + (it+0.0) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)
              elseif( ispec_is_acoustic(ispec) ) then
                 A13 = d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) &
                      + d_store_y(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML) &
                      + (it+0.5)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_y(i,j,k,ispec_CPML)
              endif
              A14 = - d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML) 

              bb = alpha_store(i,j,k,ispec_CPML)

              coef0_1 = exp(-bb * deltat) 

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb 
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0 
                 coef2_1 = deltat/2.0d0 
              end if

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_dux_dzl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 & 
                      + PML_duy_dzl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_duz_dzl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) *(it+0.5)*deltat * coef1_2 &
                      + PML_dpotential_dzl(i,j,k,ispec_CPML) *(it-0.5)*deltat * coef2_2

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
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
                 end if

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2)

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2)
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) 
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2)  

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) 
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2)

                 !---------------------- A19 and A20--------------------------
                 A19 = 1.d0
                 A20 = 0.0

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2)
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2)

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2)
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2)
              endif

           elseif(CPML_regions(ispec_CPML)==5) then
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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_2

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) 
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) 
              endif

              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A10 = d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                      + (it+0.0) * deltat * d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              elseif( ispec_is_acoustic(ispec) ) then
                 A10 = d_store_x(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML) &
                      + (it+0.5)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)
              endif
              A11 = - d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)

              coef0_1 = exp(-bb * deltat) 

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb 
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0
                 coef2_1 = deltat/2.0d0 
              end if

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1


              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_dux_dyl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 & 
                      + PML_duy_dyl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_duz_dyl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) *(it+0.5)*deltat* coef1_2 &
                      + PML_dpotential_dyl(i,j,k,ispec_CPML) *(it-0.5)*deltat* coef2_2

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_2

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
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
                 end if

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) 

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2)
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) 

                 !---------------------- A17 and A18 --------------------------
                 A17 = 1.0d0
                 A18 = 0.d0

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) 
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) 
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) 
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2)

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) 
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2)
              endif

           elseif( CPML_regions(ispec_CPML) == 6 ) then
              !------------------------------------------------------------------------------
              !---------------------------- YZ-edge C-PML -----------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)
              if( ispec_is_elastic(ispec) ) then
                 A7 = d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      + (it+0.0) * deltat * d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)
              elseif( ispec_is_acoustic(ispec) ) then
                 A7 = d_store_y(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML) &
                      + d_store_z(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) &
                      + (it+0.5)*deltat*d_store_y(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)
              endif
              A8 = - d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)

              bb = alpha_store(i,j,k,ispec_CPML)

              coef0_1 = exp(-bb * deltat) 

              if( abs(bb) > 1.d-5 ) then
                 coef1_1 = (1.d0 - exp(-bb * deltat/2.d0)) / bb 
                 coef2_1 = (1.d0 - exp(-bb * deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_1 = deltat/2.0d0 
                 coef2_1 = deltat/2.0d0 
              end if

              coef0_2 = coef0_1
              coef1_2 = coef1_1
              coef2_2 = coef2_1

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_dux_dxl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_duy_dxl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                      + PML_duz_dxl(i,j,k,ispec_CPML) * (it-0.0) * deltat * coef2_2

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2)
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * (it+0.5)*deltat* coef1_2 &
                      + PML_dpotential_dxl(i,j,k,ispec_CPML) *(it-0.5)*deltat* coef2_2

                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
              endif

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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_2

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) 
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) 
              endif                 

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
              end if

              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_2

                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_2

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_2

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
                 !---------------------- A15 and A16 --------------------------
                 A15 = 1.0d0
                 A16 = 0.0d0

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) 
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) 

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = 0.d0
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) 
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) 

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
                 end if

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2)
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2)  

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2)
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2)
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2)

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2)
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) 
              endif

           elseif( CPML_regions(ispec_CPML) == 7 ) then
              !------------------------------------------------------------------------------
              !---------------------------- XYZ-corner C-PML --------------------------------
              !------------------------------------------------------------------------------

              !---------------------- A6, A7 and A8 --------------------------
              A6 = k_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) 
              if( abs(d_store_x(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                 A7 = d_store_y(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_x(i,j,k,ispec_CPML)
                 A8 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) - &
                      d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) * &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) - &
                      d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) / &
                      ( d_store_x(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML)**2)
              else
                 if( ispec_is_elastic(ispec) ) then
                    A7 = ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) + &
                         d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                         k_store_x(i,j,k,ispec_CPML) + &
                         (it+0.0)*deltat*d_store_y(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_x(i,j,k,ispec_CPML)
                 elseif( ispec_is_acoustic(ispec) ) then
                    A7 = (d_store_z(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML)+ &
                         d_store_y(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML)) / &
                         k_store_x(i,j,k,ispec_CPML) + &
                         (it+0.5)*deltat*d_store_y(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_x(i,j,k,ispec_CPML)
                 endif
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
              end if

              bb = d_store_x(i,j,k,ispec_CPML) / k_store_x(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)

              coef0_2 = exp(-bb * deltat)

              if( abs(bb) > 1.d-5 ) then 
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              end if


              if( ispec_is_elastic(ispec) ) then
                 
                 rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1

                 if( abs(d_store_x(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                    rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_dux_dxl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duy_dxl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dxl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duz_dxl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                 endif

                 duxdxl_x = A6 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,1) + A8 * rmemory_dux_dxl_x(i,j,k,ispec_CPML,2)
                 duydxl_y = A6 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,1) + A8 * rmemory_duy_dxl_y(i,j,k,ispec_CPML,2)
                 duzdxl_z = A6 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,1) + A8 * rmemory_duz_dxl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then
                 
                 rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_1

                 if(abs(d_store_x(i,j,k,ispec_CPML)).gt. 1.d-5)then
                    rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dxl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dxl_new(i,j,k,ispec_CPML) * (it+0.5)*deltat* coef1_2 &
                         + PML_dpotential_dxl(i,j,k,ispec_CPML) * (it-0.5)*deltat* coef2_2
                 endif

                 dpotentialdxl = A6 * PML_dpotential_dxl(i,j,k,ispec_CPML)  &
                      + A7 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,1) + A8 * rmemory_dpotential_dxl(i,j,k,ispec_CPML,2)
              endif


              !---------------------- A9, A10 and A11 --------------------------
              A9 = k_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML)  
              if( abs(d_store_y(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                 A10 = d_store_x(i,j,k,ispec_CPML) * d_store_z(i,j,k,ispec_CPML)/d_store_y(i,j,k,ispec_CPML)
                 A11 = ( d_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) ) / &
                      ( d_store_y(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML)**2)
              else
                 if( ispec_is_elastic(ispec) ) then
                    A10 = ( d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                         + d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                         k_store_y(i,j,k,ispec_CPML) + &
                         (it+0.0)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_y(i,j,k,ispec_CPML)
                 elseif( ispec_is_acoustic(ispec) ) then
                    A10 = (d_store_z(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML) &
                         +d_store_x(i,j,k,ispec_CPML)*k_store_z(i,j,k,ispec_CPML)) / &
                         k_store_y(i,j,k,ispec_CPML) + &
                         (it+0.5)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_z(i,j,k,ispec_CPML)/k_store_y(i,j,k,ispec_CPML)
                 endif
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
              end if

              bb = d_store_y(i,j,k,ispec_CPML) / k_store_y(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)

              coef0_2 = exp(-bb * deltat) 

              if( abs(bb) > 1.d-5 ) then 
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0 
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1

                 if( abs(d_store_y(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                    rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_dux_dyl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duy_dyl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dyl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duz_dyl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                 endif

                 duxdyl_x = A9 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,1) + A11 * rmemory_dux_dyl_x(i,j,k,ispec_CPML,2)
                 duydyl_y = A9 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,1) + A11 * rmemory_duy_dyl_y(i,j,k,ispec_CPML,2)
                 duzdyl_z = A9 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,1) + A11 * rmemory_duz_dyl_z(i,j,k,ispec_CPML,2) 

              elseif( ispec_is_acoustic(ispec) ) then
                 rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_1

                 if(abs(d_store_y(i,j,k,ispec_CPML)).gt. 1.d-5)then
                    rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dyl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dyl_new(i,j,k,ispec_CPML) * (it+0.5)*deltat* coef1_2 &
                         + PML_dpotential_dyl(i,j,k,ispec_CPML) * (it-0.5)*deltat* coef2_2
                 endif

                 dpotentialdyl = A9 * PML_dpotential_dyl(i,j,k,ispec_CPML)  &
                      + A10 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,1) + A11 * rmemory_dpotential_dyl(i,j,k,ispec_CPML,2)
              endif

              !---------------------- A12, A13 and A14 -------------------------- 
              A12 = k_store_x(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML)  
              if( abs(d_store_z(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                 A13 = d_store_x(i,j,k,ispec_CPML) * d_store_y(i,j,k,ispec_CPML)/d_store_z(i,j,k,ispec_CPML)
                 A14 = ( d_store_x(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) &
                      - d_store_z(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) ) * &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_y(i,j,k,ispec_CPML) &
                      - d_store_y(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML) ) / &
                      ( d_store_z(i,j,k,ispec_CPML) * k_store_z(i,j,k,ispec_CPML)**2)
              else
                 if( ispec_is_elastic(ispec) ) then
                    A13 = ( d_store_y(i,j,k,ispec_CPML) * k_store_x(i,j,k,ispec_CPML) &
                         + d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML) ) / &
                         k_store_z(i,j,k,ispec_CPML) + &
                         (it+0.0)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_y(i,j,k,ispec_CPML)/k_store_z(i,j,k,ispec_CPML)
                 elseif( ispec_is_acoustic(ispec) ) then
                    A13 = (d_store_y(i,j,k,ispec_CPML)*k_store_x(i,j,k,ispec_CPML)&
                         +d_store_x(i,j,k,ispec_CPML)*k_store_y(i,j,k,ispec_CPML)) / &
                         k_store_z(i,j,k,ispec_CPML) + &
                         (it+0.5)*deltat*d_store_x(i,j,k,ispec_CPML)*d_store_y(i,j,k,ispec_CPML)/k_store_z(i,j,k,ispec_CPML)
                 endif
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
              end if

              bb = d_store_z(i,j,k,ispec_CPML) / k_store_z(i,j,k,ispec_CPML) + alpha_store(i,j,k,ispec_CPML)

              coef0_2 = exp(-bb * deltat) 

              if( abs(bb) > 1.d-5 ) then 
                 coef1_2 = (1.d0 - exp(-bb * deltat/2.d0)) / bb
                 coef2_2 = (1.d0 - exp(-bb* deltat/2.d0)) * exp(-bb * deltat/2.d0) / bb
              else
                 coef1_2 = deltat/2.0d0
                 coef2_2 = deltat/2.0d0
              end if

              if( ispec_is_elastic(ispec) ) then

                 rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1

                 if( abs(d_store_z(i,j,k,ispec_CPML)) .gt. 1.d-5 ) then
                    rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_2
                    rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2) &  
                         + PML_dux_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_dux_dzl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2) &  
                         + PML_duy_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duy_dzl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                    rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2) &  
                         + PML_duz_dzl_new(i,j,k,ispec_CPML) * (it+0.0) * deltat * coef1_2 &
                         + PML_duz_dzl(i,j,k,ispec_CPML) * (it-0.0)*deltat * coef2_2
                 endif

                 duxdzl_x = A12 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,1) + A14 * rmemory_dux_dzl_x(i,j,k,ispec_CPML,2)
                 duydzl_y = A12 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,1) + A14 * rmemory_duy_dzl_y(i,j,k,ispec_CPML,2)
                 duzdzl_z = A12 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,1) + A14 * rmemory_duz_dzl_z(i,j,k,ispec_CPML,2)

              elseif( ispec_is_acoustic(ispec) ) then

                 rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) &
                      + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_1

                 if(abs(d_store_z(i,j,k,ispec_CPML)).gt. 1.d-5)then
                    rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * coef1_2 + PML_dpotential_dzl(i,j,k,ispec_CPML) * coef2_2
                 else
                    rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) = coef0_2 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2) &
                         + PML_dpotential_dzl_new(i,j,k,ispec_CPML) * (it+0.5)*deltat* coef1_2 &
                         + PML_dpotential_dzl(i,j,k,ispec_CPML) * (it-0.5)*deltat* coef2_2
                 endif

                 dpotentialdzl = A12 * PML_dpotential_dzl(i,j,k,ispec_CPML)  &
                      + A13 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,1) + A14 * rmemory_dpotential_dzl(i,j,k,ispec_CPML,2)
              endif

              if( ispec_is_elastic(ispec) ) then
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
                 end if

                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_y = A15 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dzl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_y(i,j,k,ispec_CPML,2)
                 duzdyl_y = A15 * PML_duz_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duz_dyl_y(i,j,k,ispec_CPML,1) + rmemory_duz_dyl_y(i,j,k,ispec_CPML,2)

                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duydzl_z = A15 * PML_duy_dzl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dzl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dzl_z(i,j,k,ispec_CPML,2)
                 duydyl_z = A15 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A16 * rmemory_duy_dyl_z(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duz_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duz_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duzdzl_x = A17 * PML_duz_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dzl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dzl_x(i,j,k,ispec_CPML,2) 
                 duzdxl_x = A17 * PML_duz_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_duz_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duz_dxl_x(i,j,k,ispec_CPML,2) 

                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dzl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dzl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dzl_z(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_z(i,j,k,ispec_CPML,2) = 0.d0

                 duxdzl_z = A17 * PML_dux_dzl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dzl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dzl_z(i,j,k,ispec_CPML,2)
                 duxdxl_z = A17 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A18 * rmemory_dux_dxl_z(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_z(i,j,k,ispec_CPML,2)

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
                 end if

                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) &  
                      + PML_duy_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_duy_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_duy_dxl_x(i,j,k,ispec_CPML,2) = 0.d0

                 duydyl_x = A19 * PML_duy_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dyl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dyl_x(i,j,k,ispec_CPML,2) 
                 duydxl_x = A19 * PML_duy_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_duy_dxl_x(i,j,k,ispec_CPML,1) + rmemory_duy_dxl_x(i,j,k,ispec_CPML,2)

                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dyl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dyl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) = 0.d0

                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) = coef0_1 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) &  
                      + PML_dux_dxl_new(i,j,k,ispec_CPML) * coef1_1 + PML_dux_dxl(i,j,k,ispec_CPML) * coef2_1
                 rmemory_dux_dxl_y(i,j,k,ispec_CPML,2) = 0.d0

                 duxdyl_y = A19 * PML_dux_dyl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dyl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dyl_y(i,j,k,ispec_CPML,2) 
                 duxdxl_y = A19 * PML_dux_dxl(i,j,k,ispec_CPML)  &
                      + A20 * rmemory_dux_dxl_y(i,j,k,ispec_CPML,1) + rmemory_dux_dxl_y(i,j,k,ispec_CPML,2)
              endif

           endif ! CPML_regions

           if( ispec_is_elastic(ispec) ) then
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
           endif

        enddo
     enddo
  enddo

end subroutine pml_set_memory_variables
