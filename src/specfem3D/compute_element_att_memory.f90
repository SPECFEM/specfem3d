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

subroutine compute_element_att_memory_second_order_rk(ispec,alphaval,betaval,gammaval,NSPEC_AB,kappastore,mustore, &
                          NSPEC_ATTENUATION_AB, &
                          factor_common_kappa, &
                          R_trace,epsilondev_trace,epsilondev_trace_loc, &
                          factor_common,R_xx,R_yy,R_xy,R_xz,R_yz, &
                          NSPEC_STRAIN_ONLY,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                          epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)

  use constants, only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: ispec,NSPEC_AB,NSPEC_ATTENUATION_AB,NSPEC_STRAIN_ONLY
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(in) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(in) :: factor_common

  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(inout) :: &
            R_trace,R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(in) :: &
            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(in) :: epsilondev_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc

! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: Sn,Snp1
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: factor_loc

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! bulk attenuation
        ! term in trace
        factor_loc(:) = kappastore(i,j,k,ispec)* factor_common_kappa(:,i,j,k,ispec)

        Sn   = epsilondev_trace(i,j,k,ispec)
        Snp1   = epsilondev_trace_loc(i,j,k)
        R_trace(:,i,j,k,ispec) = alphaval(:) * R_trace(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

        ! shear attenuation
        ! term in xx yy zz xy xz yz
        factor_loc(:) = mustore(i,j,k,ispec) * factor_common(:,i,j,k,ispec)

        ! term in xx
        Sn   = epsilondev_xx(i,j,k,ispec)
        Snp1   = epsilondev_xx_loc(i,j,k)
        R_xx(:,i,j,k,ispec) = alphaval(:) * R_xx(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

        ! term in yy
        Sn   = epsilondev_yy(i,j,k,ispec)
        Snp1   = epsilondev_yy_loc(i,j,k)
        R_yy(:,i,j,k,ispec) = alphaval(:) * R_yy(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

        ! term in zz not computed since zero trace

        ! term in xy
        Sn   = epsilondev_xy(i,j,k,ispec)
        Snp1   = epsilondev_xy_loc(i,j,k)
        R_xy(:,i,j,k,ispec) = alphaval(:) * R_xy(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

        ! term in xz
        Sn   = epsilondev_xz(i,j,k,ispec)
        Snp1   = epsilondev_xz_loc(i,j,k)
        R_xz(:,i,j,k,ispec) = alphaval(:) * R_xz(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

        ! term in yz
        Sn   = epsilondev_yz(i,j,k,ispec)
        Snp1   = epsilondev_yz_loc(i,j,k)
        R_yz(:,i,j,k,ispec) = alphaval(:) * R_yz(:,i,j,k,ispec) + factor_loc(:) * (betaval(:) * Sn + gammaval(:) * Snp1)

      enddo
    enddo
  enddo

end subroutine compute_element_att_memory_second_order_rk

!
!--------------------------------------------------------------------------------------------
!

subroutine compute_element_att_memory_lddrk(ispec,deltat,NSPEC_AB,kappastore,mustore, &
                          NSPEC_ATTENUATION_AB, &
                          factor_common_kappa, &
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

  integer,intent(in) :: ispec,NSPEC_AB,NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_LDDRK
  real(kind=CUSTOM_REAL),intent(in) :: deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(in) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(in) :: factor_common

  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(inout) :: R_trace
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),intent(inout) :: R_trace_lddrk

  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(inout) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),intent(inout) :: &
            R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc

! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: Snp1
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: factor_loc

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! term in trace
        factor_loc(:) = kappastore(i,j,k,ispec) * factor_common_kappa(:,i,j,k,ispec)
        Snp1   = epsilondev_trace_loc(i,j,k)

        R_trace_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_trace_lddrk(:,i,j,k,ispec) + &
          deltat * (factor_loc(:)*Snp1 - R_trace(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))
        R_trace(:,i,j,k,ispec) = R_trace(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_trace_lddrk(:,i,j,k,ispec)

        ! term in xx yy zz xy xz yz
        factor_loc(:) = mustore(i,j,k,ispec) * factor_common(:,i,j,k,ispec)

        ! term in xx
        Snp1   = epsilondev_xx_loc(i,j,k)
        R_xx_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_xx_lddrk(:,i,j,k,ispec) + &
            deltat * (factor_loc(:)*Snp1 - R_xx(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))

        ! term in yy
        Snp1   = epsilondev_yy_loc(i,j,k)
        R_yy_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_yy_lddrk(:,i,j,k,ispec) + &
            deltat * (factor_loc(:)*Snp1 - R_yy(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))

        ! term in zz not computed since zero trace

        ! term in xy
        Snp1   = epsilondev_xy_loc(i,j,k)
        R_xy_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_xy_lddrk(:,i,j,k,ispec) + &
            deltat * (factor_loc(:)*Snp1 - R_xy(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))

        ! term in xz
        Snp1   = epsilondev_xz_loc(i,j,k)
        R_xz_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_xz_lddrk(:,i,j,k,ispec) + &
            deltat * (factor_loc(:)*Snp1 - R_xz(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))

        ! term in yz
        Snp1   = epsilondev_yz_loc(i,j,k)
        R_yz_lddrk(:,i,j,k,ispec) = ALPHA_LDDRK(istage) * R_yz_lddrk(:,i,j,k,ispec) + &
            deltat * (factor_loc(:)*Snp1 - R_yz(:,i,j,k,ispec)*(1.0_CUSTOM_REAL/tau_sigma(:)))

        R_xx(:,i,j,k,ispec) = R_xx(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(:,i,j,k,ispec)
        R_yy(:,i,j,k,ispec) = R_yy(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_yy_lddrk(:,i,j,k,ispec)
        R_xy(:,i,j,k,ispec) = R_xy(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_xy_lddrk(:,i,j,k,ispec)
        R_xz(:,i,j,k,ispec) = R_xz(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_xz_lddrk(:,i,j,k,ispec)
        R_yz(:,i,j,k,ispec) = R_yz(:,i,j,k,ispec) + BETA_LDDRK(istage) * R_yz_lddrk(:,i,j,k,ispec)
      enddo
    enddo
  enddo

end subroutine compute_element_att_memory_lddrk
