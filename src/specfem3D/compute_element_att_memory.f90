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

subroutine compute_element_att_memory_second_order_rk(ispec,alphaval,betaval,gammaval,NSPEC_AB,kappastore,mustore, &
                          NSPEC_ATTENUATION_AB,factor_common_kappa, &
                          R_trace,epsilondev_trace,epsilondev_trace_loc, &
                          factor_common,R_xx,R_yy,R_xy,R_xz,R_yz, &
                          NSPEC_STRAIN_ONLY,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                          epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)

  use constants, only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ

  implicit none

  integer :: ispec,NSPEC_AB,NSPEC_ATTENUATION_AB,NSPEC_STRAIN_ONLY
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: R_trace
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: epsilondev_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: &
            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc

! local parameters
  integer :: i_sls,i,j,k
  real(kind=CUSTOM_REAL) :: factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1

  do i_sls = 1,N_SLS
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          alphaval_loc = alphaval(i_sls)
          betaval_loc = betaval(i_sls)
          gammaval_loc = gammaval(i_sls)

          ! term in trace
          factor_loc = kappastore(i,j,k,ispec) * factor_common_kappa(i_sls,i,j,k,ispec)

          Sn   = factor_loc * epsilondev_trace(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_trace_loc(i,j,k)
          R_trace(i,j,k,ispec,i_sls) = alphaval_loc * R_trace(i,j,k,ispec,i_sls) + &
                                    betaval_loc * Sn + gammaval_loc * Snp1

          ! term in xx yy zz xy xz yz
          factor_loc = mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

          ! term in xx
          Sn   = factor_loc * epsilondev_xx(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_xx_loc(i,j,k)
          R_xx(i,j,k,ispec,i_sls) = alphaval_loc * R_xx(i,j,k,ispec,i_sls) + &
                     betaval_loc * Sn + gammaval_loc * Snp1
          ! term in yy
          Sn   = factor_loc * epsilondev_yy(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_yy_loc(i,j,k)
          R_yy(i,j,k,ispec,i_sls) = alphaval_loc * R_yy(i,j,k,ispec,i_sls) + &
                     betaval_loc * Sn + gammaval_loc * Snp1

          ! term in zz not computed since zero trace

          ! term in xy
          Sn   = factor_loc * epsilondev_xy(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_xy_loc(i,j,k)
          R_xy(i,j,k,ispec,i_sls) = alphaval_loc * R_xy(i,j,k,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
          ! term in xz
          Sn   = factor_loc * epsilondev_xz(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_xz_loc(i,j,k)
          R_xz(i,j,k,ispec,i_sls) = alphaval_loc * R_xz(i,j,k,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
          ! term in yz
          Sn   = factor_loc * epsilondev_yz(i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_yz_loc(i,j,k)
          R_yz(i,j,k,ispec,i_sls) = alphaval_loc * R_yz(i,j,k,ispec,i_sls) + &
                     betaval_loc * Sn + gammaval_loc * Snp1
        enddo
      enddo
    enddo
  enddo   ! end of loop on memory variables

end subroutine compute_element_att_memory_second_order_rk

!
!--------------------------------------------------------------------------------------------
!

subroutine compute_element_att_memory_lddrk(ispec,deltat,NSPEC_AB,kappastore,mustore, &
                          NSPEC_ATTENUATION_AB,factor_common_kappa, &
                          R_trace,epsilondev_trace_loc, &
                          NSPEC_ATTENUATION_AB_LDDRK,R_trace_lddrk, &
                          factor_common,R_xx,R_yy,R_xy,R_xz,R_yz, &
                          R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                          epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc, &
                          epsilondev_xz_loc,epsilondev_yz_loc)

  use constants, only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ,ALPHA_LDDRK,BETA_LDDRK
  use specfem_par, only: istage
  use specfem_par_elastic, only: tau_sigma

  implicit none

  integer :: ispec,NSPEC_AB,NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_LDDRK
  real(kind=CUSTOM_REAL) :: deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: R_trace
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS) :: R_trace_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS) :: &
            R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc

! local parameters
  integer :: i_sls,i,j,k
  real(kind=CUSTOM_REAL) :: factor_loc,Snp1

  do i_sls = 1,N_SLS
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! term in trace
          factor_loc = kappastore(i,j,k,ispec) * factor_common_kappa(i_sls,i,j,k,ispec)
          Snp1   = factor_loc * epsilondev_trace_loc(i,j,k)

          R_trace_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_trace_lddrk(i,j,k,ispec,i_sls) + &
            deltat * (Snp1 - R_trace(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))
          R_trace(i,j,k,i_SLS,ispec) = R_trace(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_trace_lddrk(i,j,k,i_SLS,ispec)

          ! term in xx yy zz xy xz yz
          factor_loc = mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

          ! term in xx
          Snp1   = factor_loc * epsilondev_xx_loc(i,j,k)
          R_xx_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_xx_lddrk(i,j,k,ispec,i_sls) + &
              deltat * (Snp1 - R_xx(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))

          ! term in yy
          Snp1   = factor_loc * epsilondev_yy_loc(i,j,k)
          R_yy_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_yy_lddrk(i,j,k,ispec,i_sls) + &
              deltat * (Snp1 - R_yy(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))

          ! term in zz not computed since zero trace

          ! term in xy
          Snp1   = factor_loc * epsilondev_xy_loc(i,j,k)
          R_xy_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_xy_lddrk(i,j,k,ispec,i_sls) + &
              deltat * (Snp1 - R_xy(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))

          ! term in xz
          Snp1   = factor_loc * epsilondev_xz_loc(i,j,k)
          R_xz_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_xz_lddrk(i,j,k,ispec,i_sls) + &
              deltat * (Snp1 - R_xz(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))

          ! term in yz
          Snp1   = factor_loc * epsilondev_yz_loc(i,j,k)
          R_yz_lddrk(i,j,k,ispec,i_sls) = ALPHA_LDDRK(istage) * R_yz_lddrk(i,j,k,ispec,i_sls) + &
              deltat * (Snp1 - R_yz(i,j,k,ispec,i_sls)*(1.0_CUSTOM_REAL/tau_sigma(i_SLS)))

          R_xx(i,j,k,i_SLS,ispec) = R_xx(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(i,j,k,i_SLS,ispec)
          R_yy(i,j,k,i_SLS,ispec) = R_yy(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_yy_lddrk(i,j,k,i_SLS,ispec)
          R_xy(i,j,k,i_SLS,ispec) = R_xy(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_xy_lddrk(i,j,k,i_SLS,ispec)
          R_xz(i,j,k,i_SLS,ispec) = R_xz(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_xz_lddrk(i,j,k,i_SLS,ispec)
          R_yz(i,j,k,i_SLS,ispec) = R_yz(i,j,k,i_SLS,ispec) + BETA_LDDRK(istage) * R_yz_lddrk(i,j,k,i_SLS,ispec)
        enddo
      enddo
    enddo
  enddo   ! end of loop on memory variables

end subroutine compute_element_att_memory_lddrk
