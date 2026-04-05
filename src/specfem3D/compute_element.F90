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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


!--------------------------------------------------------------------------------------------
!
! isotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_iso(ispec,ispec_irreg, &
                                 minus_g,minus_deriv_gravity,rho_s_H, &
                                 xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                                 gammaxstore,gammaystore,gammazstore,jacobianstore, &
                                 duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
                                 wgll_cube, &
                                 kappastore,mustore, &
                                 ibool, &
                                 R_xx,R_yy,R_xy,R_xz,R_yz,R_trace, &
                                 tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                 dummyx_loc,dummyy_loc,dummyz_loc)

! isotropic element in viscoelastic domain

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS

  use shared_parameters, only: ATTENUATION, GRAVITY, MOVIE_VOLUME_STRESS

  use specfem_par_coupling, only: do_save_coupling_wavefield

  use specfem_par, only: &
    NSPEC => NSPEC_AB, &
    NGLOB => NGLOB_AB, &
    NSPEC_ATTENUATION => NSPEC_ATTENUATION_AB

  use specfem_par, only: xix_regular, jacobian_regular

  ! PML
  use pml_par, only: is_CPML

  ! movie
  use specfem_par_movie, only: stress_xx,stress_yy,stress_zz,stress_xy,stress_xz,stress_yz

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec,ispec_irreg

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,jacobianstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(in) :: &
    duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! isotropic properties
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: kappastore,mustore

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION),intent(in) :: R_trace

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: minus_g,minus_deriv_gravity
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: R_trace_kappa_sum,R_xx_sum,R_yy_sum

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

  real(kind=CUSTOM_REAL), parameter :: FOUR_THIRDS = 4.0_CUSTOM_REAL / 3.0_CUSTOM_REAL

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  !
  ! compute stress for isotropic element
  !
  DO_LOOP_IJK
    ! elastic parameters
    kappal = kappastore(INDEX_IJK,ispec)
    mul = mustore(INDEX_IJK,ispec)

    lambdalplus2mul = kappal + FOUR_THIRDS * mul
    lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

    ! compute stress sigma
    sigma_xx(INDEX_IJK) = lambdalplus2mul * duxdxl(INDEX_IJK) + lambdal * (duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
    sigma_yy(INDEX_IJK) = lambdalplus2mul * duydyl(INDEX_IJK) + lambdal * (duxdxl(INDEX_IJK) + duzdzl(INDEX_IJK))
    sigma_zz(INDEX_IJK) = lambdalplus2mul * duzdzl(INDEX_IJK) + lambdal * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK))

    sigma_xy(INDEX_IJK) = mul * (duxdyl(INDEX_IJK) + duydxl(INDEX_IJK))
    sigma_xz(INDEX_IJK) = mul * (duzdxl(INDEX_IJK) + duxdzl(INDEX_IJK))
    sigma_yz(INDEX_IJK) = mul * (duzdyl(INDEX_IJK) + duydzl(INDEX_IJK))

    ! subtract memory variables if attenuation
    if (ATTENUATION .and. .not. is_CPML(ispec)) then
      R_xx_sum = sum(R_xx(:,INDEX_IJK,ispec))
      R_yy_sum = sum(R_yy(:,INDEX_IJK,ispec))
      R_trace_kappa_sum = sum(R_trace(:,INDEX_IJK,ispec))

      ! in case no bulk attenuation is desired:
      !R_trace_kappa_sum = 0.0

      sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_sum - R_trace_kappa_sum
      sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_sum - R_trace_kappa_sum
      sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_sum + R_yy_sum - R_trace_kappa_sum
      sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - sum(R_xy(:,INDEX_IJK,ispec))
      sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - sum(R_xz(:,INDEX_IJK,ispec))
      sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - sum(R_yz(:,INDEX_IJK,ispec))
    endif

    ! define symmetric components of sigma
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! stores stress for movie output
  ! and SPECFEM coupling injection technique to compute traction on boundary point
  if (MOVIE_VOLUME_STRESS .or. do_save_coupling_wavefield) then
    ! store stress tensor
    stress_xx(:,:,:,ispec) = sigma_xx(:,:,:)
    stress_yy(:,:,:,ispec) = sigma_yy(:,:,:)
    stress_zz(:,:,:,ispec) = sigma_zz(:,:,:)
    stress_xy(:,:,:,ispec) = sigma_xy(:,:,:)
    stress_xz(:,:,:,ispec) = sigma_xz(:,:,:)
    stress_yz(:,:,:,ispec) = sigma_yz(:,:,:)
  endif

  ! compute non-symmetric terms for gravity
  if (GRAVITY) then
    call compute_element_gravity(ispec,ispec_irreg,NSPEC,NGLOB,ibool, &
                                 jacobianstore,wgll_cube, &
                                 minus_g,minus_deriv_gravity, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product with test vector
  if (.not. is_CPML(ispec)) then
    if (ispec_irreg /= 0) then
      ! irregular element
      DO_LOOP_IJK
        xixl = xixstore(INDEX_IJK,ispec_irreg)
        xiyl = xiystore(INDEX_IJK,ispec_irreg)
        xizl = xizstore(INDEX_IJK,ispec_irreg)
        etaxl = etaxstore(INDEX_IJK,ispec_irreg)
        etayl = etaystore(INDEX_IJK,ispec_irreg)
        etazl = etazstore(INDEX_IJK,ispec_irreg)
        gammaxl = gammaxstore(INDEX_IJK,ispec_irreg)
        gammayl = gammaystore(INDEX_IJK,ispec_irreg)
        gammazl = gammazstore(INDEX_IJK,ispec_irreg)
        jacobianl = jacobianstore(INDEX_IJK,ispec_irreg)

        ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
        tempx1(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * xixl + sigma_yx(INDEX_IJK) * xiyl + sigma_zx(INDEX_IJK) * xizl)
        tempy1(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * xixl + sigma_yy(INDEX_IJK) * xiyl + sigma_zy(INDEX_IJK) * xizl)
        tempz1(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * xixl + sigma_yz(INDEX_IJK) * xiyl + sigma_zz(INDEX_IJK) * xizl)

        tempx2(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * etaxl + sigma_yx(INDEX_IJK) * etayl + sigma_zx(INDEX_IJK) * etazl)
        tempy2(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * etaxl + sigma_yy(INDEX_IJK) * etayl + sigma_zy(INDEX_IJK) * etazl)
        tempz2(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * etaxl + sigma_yz(INDEX_IJK) * etayl + sigma_zz(INDEX_IJK) * etazl)

        tempx3(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * gammaxl + sigma_yx(INDEX_IJK) * gammayl + sigma_zx(INDEX_IJK) * gammazl)
        tempy3(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * gammaxl + sigma_yy(INDEX_IJK) * gammayl + sigma_zy(INDEX_IJK) * gammazl)
        tempz3(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * gammaxl + sigma_yz(INDEX_IJK) * gammayl + sigma_zz(INDEX_IJK) * gammazl)
      ENDDO_LOOP_IJK
    else
      ! regular element
      DO_LOOP_IJK
        ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
        tempx1(INDEX_IJK) = jacobian_regular * sigma_xx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy1(INDEX_IJK) = jacobian_regular * sigma_xy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz1(INDEX_IJK) = jacobian_regular * sigma_xz(INDEX_IJK) * xix_regular ! this goes to accel_z

        tempx2(INDEX_IJK) = jacobian_regular * sigma_yx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy2(INDEX_IJK) = jacobian_regular * sigma_yy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz2(INDEX_IJK) = jacobian_regular * sigma_yz(INDEX_IJK) * xix_regular ! this goes to accel_z

        tempx3(INDEX_IJK) = jacobian_regular * sigma_zx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy3(INDEX_IJK) = jacobian_regular * sigma_zy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz3(INDEX_IJK) = jacobian_regular * sigma_zz(INDEX_IJK) * xix_regular ! this goes to accel_z
      ENDDO_LOOP_IJK
    endif
  endif


  end subroutine compute_element_iso


!--------------------------------------------------------------------------------------------
!
! anisotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_aniso(ispec,ispec_irreg, &
                                   minus_g,minus_deriv_gravity,rho_s_H, &
                                   xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                                   gammaxstore,gammaystore,gammazstore,jacobianstore, &
                                   duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
                                   wgll_cube, &
                                   c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                   c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                   c36store,c44store,c45store,c46store,c55store,c56store,c66store, &                                   
                                   ibool, &
                                   R_xx,R_yy,R_xy,R_xz,R_yz,R_trace, &
                                   tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                   dummyx_loc,dummyy_loc,dummyz_loc)

! fully anisotropic element in viscoelastic domain

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS

  use shared_parameters, only: ATTENUATION, GRAVITY, MOVIE_VOLUME_STRESS

  use specfem_par_coupling, only: do_save_coupling_wavefield

  use specfem_par, only: &
    NSPEC => NSPEC_AB, &
    NGLOB => NGLOB_AB, &
    NSPEC_ATTENUATION => NSPEC_ATTENUATION_AB

  use specfem_par_elastic, only: NSPEC_ANISO

  use specfem_par, only: xix_regular, jacobian_regular

  ! PML
  use pml_par, only: is_CPML

  ! movie
  use specfem_par_movie, only: stress_xx,stress_yy,stress_zz,stress_xy,stress_xz,stress_yz

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec,ispec_irreg

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,jacobianstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ), intent(in) :: &
    duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),intent(in) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION),intent(in) :: R_trace

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: minus_g,minus_deriv_gravity
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: R_trace_kappa_sum,R_xx_sum,R_yy_sum

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  !
  ! compute stress for isotropic element
  !
  DO_LOOP_IJK
    ! full anisotropic case, stress calculations
    c11 = c11store(INDEX_IJK,ispec)
    c12 = c12store(INDEX_IJK,ispec)
    c13 = c13store(INDEX_IJK,ispec)
    c14 = c14store(INDEX_IJK,ispec)
    c15 = c15store(INDEX_IJK,ispec)
    c16 = c16store(INDEX_IJK,ispec)
    c22 = c22store(INDEX_IJK,ispec)
    c23 = c23store(INDEX_IJK,ispec)
    c24 = c24store(INDEX_IJK,ispec)
    c25 = c25store(INDEX_IJK,ispec)
    c26 = c26store(INDEX_IJK,ispec)
    c33 = c33store(INDEX_IJK,ispec)
    c34 = c34store(INDEX_IJK,ispec)
    c35 = c35store(INDEX_IJK,ispec)
    c36 = c36store(INDEX_IJK,ispec)
    c44 = c44store(INDEX_IJK,ispec)
    c45 = c45store(INDEX_IJK,ispec)
    c46 = c46store(INDEX_IJK,ispec)
    c55 = c55store(INDEX_IJK,ispec)
    c56 = c56store(INDEX_IJK,ispec)
    c66 = c66store(INDEX_IJK,ispec)

    ! precompute some sums to save CPU time
    duxdyl_plus_duydxl = duxdyl(INDEX_IJK) + duydxl(INDEX_IJK)
    duzdxl_plus_duxdzl = duzdxl(INDEX_IJK) + duxdzl(INDEX_IJK)
    duzdyl_plus_duydzl = duzdyl(INDEX_IJK) + duydzl(INDEX_IJK)

    ! compute stress sigma
    sigma_xx(INDEX_IJK) = c11 * duxdxl(INDEX_IJK) + c16 * duxdyl_plus_duydxl + c12 * duydyl(INDEX_IJK) + &
                          c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * duzdzl(INDEX_IJK)
    sigma_yy(INDEX_IJK) = c12 * duxdxl(INDEX_IJK) + c26 * duxdyl_plus_duydxl + c22 * duydyl(INDEX_IJK) + &
                          c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * duzdzl(INDEX_IJK)
    sigma_zz(INDEX_IJK) = c13 * duxdxl(INDEX_IJK) + c36 * duxdyl_plus_duydxl + c23 * duydyl(INDEX_IJK) + &
                          c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * duzdzl(INDEX_IJK)
    sigma_xy(INDEX_IJK) = c16 * duxdxl(INDEX_IJK) + c66 * duxdyl_plus_duydxl + c26 * duydyl(INDEX_IJK) + &
                          c56 * duzdxl_plus_duxdzl + c46 * duzdyl_plus_duydzl + c36 * duzdzl(INDEX_IJK)
    sigma_xz(INDEX_IJK) = c15 * duxdxl(INDEX_IJK) + c56 * duxdyl_plus_duydxl + c25 * duydyl(INDEX_IJK) + &
                          c55 * duzdxl_plus_duxdzl + c45 * duzdyl_plus_duydzl + c35 * duzdzl(INDEX_IJK)
    sigma_yz(INDEX_IJK) = c14 * duxdxl(INDEX_IJK) + c46 * duxdyl_plus_duydxl + c24 * duydyl(INDEX_IJK) + &
                          c45 * duzdxl_plus_duxdzl + c44 * duzdyl_plus_duydzl + c34 * duzdzl(INDEX_IJK)

    ! subtract memory variables if attenuation
    if (ATTENUATION .and. .not. is_CPML(ispec)) then
      R_xx_sum = sum(R_xx(:,INDEX_IJK,ispec))
      R_yy_sum = sum(R_yy(:,INDEX_IJK,ispec))
      R_trace_kappa_sum = sum(R_trace(:,INDEX_IJK,ispec))

      ! in case no bulk attenuation is desired:
      !R_trace_kappa_sum = 0.0

      sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_sum - R_trace_kappa_sum
      sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_sum - R_trace_kappa_sum
      sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_sum + R_yy_sum - R_trace_kappa_sum
      sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - sum(R_xy(:,INDEX_IJK,ispec))
      sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - sum(R_xz(:,INDEX_IJK,ispec))
      sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - sum(R_yz(:,INDEX_IJK,ispec))
    endif

    ! define symmetric components of sigma
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! stores stress for movie output
  ! and SPECFEM coupling injection technique to compute traction on boundary point
  if (MOVIE_VOLUME_STRESS .or. do_save_coupling_wavefield) then
    ! store stress tensor
    stress_xx(:,:,:,ispec) = sigma_xx(:,:,:)
    stress_yy(:,:,:,ispec) = sigma_yy(:,:,:)
    stress_zz(:,:,:,ispec) = sigma_zz(:,:,:)
    stress_xy(:,:,:,ispec) = sigma_xy(:,:,:)
    stress_xz(:,:,:,ispec) = sigma_xz(:,:,:)
    stress_yz(:,:,:,ispec) = sigma_yz(:,:,:)
  endif

  ! compute non-symmetric terms for gravity
  if (GRAVITY) then
    call compute_element_gravity(ispec,ispec_irreg,NSPEC,NGLOB,ibool, &
                                 jacobianstore,wgll_cube, &
                                 minus_g,minus_deriv_gravity, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product with test vector
  if (.not. is_CPML(ispec)) then
    if (ispec_irreg /= 0) then
      ! irregular element
      DO_LOOP_IJK
        xixl = xixstore(INDEX_IJK,ispec_irreg)
        xiyl = xiystore(INDEX_IJK,ispec_irreg)
        xizl = xizstore(INDEX_IJK,ispec_irreg)
        etaxl = etaxstore(INDEX_IJK,ispec_irreg)
        etayl = etaystore(INDEX_IJK,ispec_irreg)
        etazl = etazstore(INDEX_IJK,ispec_irreg)
        gammaxl = gammaxstore(INDEX_IJK,ispec_irreg)
        gammayl = gammaystore(INDEX_IJK,ispec_irreg)
        gammazl = gammazstore(INDEX_IJK,ispec_irreg)
        jacobianl = jacobianstore(INDEX_IJK,ispec_irreg)

        ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
        tempx1(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * xixl + sigma_yx(INDEX_IJK) * xiyl + sigma_zx(INDEX_IJK) * xizl)
        tempy1(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * xixl + sigma_yy(INDEX_IJK) * xiyl + sigma_zy(INDEX_IJK) * xizl)
        tempz1(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * xixl + sigma_yz(INDEX_IJK) * xiyl + sigma_zz(INDEX_IJK) * xizl)

        tempx2(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * etaxl + sigma_yx(INDEX_IJK) * etayl + sigma_zx(INDEX_IJK) * etazl)
        tempy2(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * etaxl + sigma_yy(INDEX_IJK) * etayl + sigma_zy(INDEX_IJK) * etazl)
        tempz2(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * etaxl + sigma_yz(INDEX_IJK) * etayl + sigma_zz(INDEX_IJK) * etazl)

        tempx3(INDEX_IJK) = jacobianl * &
                            (sigma_xx(INDEX_IJK) * gammaxl + sigma_yx(INDEX_IJK) * gammayl + sigma_zx(INDEX_IJK) * gammazl)
        tempy3(INDEX_IJK) = jacobianl * &
                            (sigma_xy(INDEX_IJK) * gammaxl + sigma_yy(INDEX_IJK) * gammayl + sigma_zy(INDEX_IJK) * gammazl)
        tempz3(INDEX_IJK) = jacobianl * &
                            (sigma_xz(INDEX_IJK) * gammaxl + sigma_yz(INDEX_IJK) * gammayl + sigma_zz(INDEX_IJK) * gammazl)
      ENDDO_LOOP_IJK
    else
      ! regular element
      DO_LOOP_IJK
        ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
        tempx1(INDEX_IJK) = jacobian_regular * sigma_xx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy1(INDEX_IJK) = jacobian_regular * sigma_xy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz1(INDEX_IJK) = jacobian_regular * sigma_xz(INDEX_IJK) * xix_regular ! this goes to accel_z

        tempx2(INDEX_IJK) = jacobian_regular * sigma_yx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy2(INDEX_IJK) = jacobian_regular * sigma_yy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz2(INDEX_IJK) = jacobian_regular * sigma_yz(INDEX_IJK) * xix_regular ! this goes to accel_z

        tempx3(INDEX_IJK) = jacobian_regular * sigma_zx(INDEX_IJK) * xix_regular ! this goes to accel_x
        tempy3(INDEX_IJK) = jacobian_regular * sigma_zy(INDEX_IJK) * xix_regular ! this goes to accel_y
        tempz3(INDEX_IJK) = jacobian_regular * sigma_zz(INDEX_IJK) * xix_regular ! this goes to accel_z
      ENDDO_LOOP_IJK
    endif
  endif


  end subroutine compute_element_aniso

!--------------------------------------------------------------------------------------------
!
! helper functions
!
!--------------------------------------------------------------------------------------------

!
! please leave this routine in this file, to help compilers inlining this function...
!

  subroutine compute_element_gravity(ispec,ispec_irreg,NSPEC,NGLOB,ibool, &
                                          jacobianstore, wgll_cube, &
                                          minus_g,minus_deriv_gravity, &
                                          dummyx_loc,dummyy_loc,dummyz_loc, &
                                          sigma_xx,sigma_yy,sigma_zz, &
                                          sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                          rho_s_H)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES INLINE :: compute_element_gravity
#else
! cray
!DIR$ INLINEALWAYS compute_element_gravity
#endif

! computes non-symmetric stress terms for gravity

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM
  use specfem_par, only: rhostore, jacobian_regular

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  integer,intent(in) :: ispec,ispec_irreg
  integer,intent(in) :: NSPEC,NGLOB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
!  real(kind=CUSTOM_REAL), dimension(3,NGLOB),intent(in) :: rstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: jacobianstore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: minus_g,minus_deriv_gravity

!  double precision, dimension(NRAD_GRAVITY),intent(in) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H

  ! local parameters
  ! for gravity
  real(kind=CUSTOM_REAL) :: rhol
  
!  double precision :: dphi,dtheta
!  double precision :: radius,rho,minus_g,minus_dg
!  double precision :: minus_g_over_radius,minus_dg_plus_g_over_radius
!  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
!  double precision :: cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  real(kind=CUSTOM_REAL) :: factor,sx_l,sy_l,sz_l,gxl,gyl,gzl,jacobianl
  real(kind=CUSTOM_REAL) :: Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl

!  integer :: int_radius
  integer :: iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! minimum radius in inner core (to avoid zero radius)
  !double precision, parameter :: MINIMUM_RADIUS_INNER_CORE = 100.d0 / R_PLANET

  ! computes non-symmetric terms for gravity
  DO_LOOP_IJK
    ! use mesh coordinates to get theta and phi
    ! x y and z contain r theta and phi
    iglob = ibool(INDEX_IJK,ispec)

    ! density
    rhol = rhostore(INDEX_IJK,ispec)

    ! Cartesian components of the gravitational acceleration
    ! assumes g acts in vertical direction only
    gxl = 0._CUSTOM_REAL            ! (minus_g * sin_theta * cos_phi) * rho
    gyl = 0._CUSTOM_REAL            ! (minus_g * sin_theta * sin_phi) * rho
    gzl = minus_g(iglob) * rhol     ! (minus_g * cos_theta) * rho

    ! Cartesian components of gradient of gravitational acceleration
    ! get displacement and multiply by density to compute G tensor
    sx_l = dummyx_loc(INDEX_IJK)
    sy_l = dummyy_loc(INDEX_IJK)
    sz_l = dummyz_loc(INDEX_IJK)

    ! compute G tensor from s . g and add to sigma (not symmetric)
    !
    ! note: Komatitsch & Tromp 2002, Spectral-element simulations of global seismic wave propagation â II. Three-dimensional
    !       models, oceans, rotation and self-gravitation, GJI, 150, 303-318
    !       https://doi.org/10.1046/j.1365-246X.2002.01716.x
    !
    !       G is defined as G = \rho [ s g - (s \cdot g) I ]. G is non-symmetric.
    !
    !       Here, the contribution added to the elastic stress tensor (sigma, or T in the paper) is:
    !          (s \cdot (\rho g)) I - s (\rho g) = \rho [ (s \cdot g) I - s g ]
    !
    !       In index notation the contribution here is:
    !          \rho [ (s_k g_k) \delta_ij - s_i g_j ]
    !       That is, the contribution added is - G.
    !
    !       This will lead to a formulation of the weak form for stress as:
    !          - int_\Omega \nabla w : (T - G) d\Omega
    !       Note that the sign of G is different to the expression in the paper, where the derivation of the weak form seems
    !       to contain a sign mistake in the (T - G) term.
    !
    !       contribution \rho [ (s_k g_k) \delta_ij - s_i g_j ]:
    !       for example: xx (i=1, j=1): (sx * gx + sy * gy + sz * gz) * 1 - sx gx == sy * gy + sz * gz
    !                    yy (i=2, j=2): (sx * gx + sy * gy + sz * gz) * 1 - sy gy == sx * gx + sz * gz
    !                    zz (i=3, j=3): (sx * gx + sy * gy + sz * gz) * 1 - sz gz == sx * gx + sy * gy
    !                    xy (i=1, j=2): (sx * gx + sy * gy + sz * gz) * 0 - sx gy == - sx * gy
    !                    ..
    !
    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) + sy_l * gyl + sz_l * gzl
    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) + sx_l * gxl + sz_l * gzl
    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + sx_l * gxl + sy_l * gyl

    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - sx_l * gyl
    sigma_yx(INDEX_IJK) = sigma_yx(INDEX_IJK) - sy_l * gxl

    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - sx_l * gzl
    sigma_zx(INDEX_IJK) = sigma_zx(INDEX_IJK) - sz_l * gxl

    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - sy_l * gzl
    sigma_zy(INDEX_IJK) = sigma_zy(INDEX_IJK) - sz_l * gyl

    ! H term contribution
    ! note: this computes term \rho s \cdot H
    !       Since H is defined as H = \nabla g and g = - \nabla \Psi with the gravitational potential \Psi,
    !       the resulting tensor H = \nabla \nabla \Psi' must be symmetric (using \Psi' == -\Psi).
    !
    !       (And for a symmetric tensor H, the product s \cdot H == H^T \cdot s == H \cdot s,
    !        with H^T being the transpose of H)
    !
    !       Note that the H term in the Komatitsch & Tromp 2002 paper seems to have a sign mistake as well, and should be
    !          + \int_\Omega \rho s \cdot H \cdot w d\Omega
    !
    !       contribution v = s \cdot H' with H' = \rho H, in index notation v_i = s_j H'_ji:
    !       for example: vx: sx * Hxx + sy * Hyx + sz * Hzx == sx * Hxx + sy * Hxy + sz * Hxz (since Hyx == Hxy, Hzx == Hxz)
    !                    vy: sx * Hxy + sy * Hyy + sz * Hzy == sx * Hxy + sy * Hyy + sz * Hyz (since Hyx == Hxy, Hzy == Hyz)
    !                    vz: sx * Hxz + sy * Hyz + sz * Hzz
    !
    !
    ! H-matrix
    ! Cartesian components of gradient of gravitational acceleration
    ! obtained from spherical components
    !minus_g_over_radius = minus_g / radius
    !minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius
    !Hxxl = (minus_g_over_radius * (cos_phi_sq * cos_theta_sq + sin_phi_sq) + cos_phi_sq * minus_dg * sin_theta_sq)
    !Hyyl = (minus_g_over_radius * (cos_phi_sq + cos_theta_sq * sin_phi_sq) + minus_dg * sin_phi_sq * sin_theta_sq)
    !Hzzl = (cos_theta_sq * minus_dg + minus_g_over_radius * sin_theta_sq)
    !Hxyl = (cos_phi * minus_dg_plus_g_over_radius * sin_phi * sin_theta_sq)
    !Hxzl = (cos_phi * cos_theta * minus_dg_plus_g_over_radius * sin_theta)
    !Hyzl = (cos_theta * minus_dg_plus_g_over_radius * sin_phi * sin_theta)

    ! vertical: colatitude theta == 0 -> cos_theta = 1 & sin_theta = 0
    ! the general expression
    !Hxxl = (minus_g_over_radius * (cos_phi_sq * cos_theta_sq + sin_phi_sq) + cos_phi_sq * minus_dg * sin_theta_sq)
    !Hyyl = (minus_g_over_radius * (cos_phi_sq + cos_theta_sq * sin_phi_sq) + minus_dg * sin_phi_sq * sin_theta_sq)
    !Hzzl = (cos_theta_sq * minus_dg + minus_g_over_radius * sin_theta_sq)
    !Hxyl = (cos_phi * minus_dg_plus_g_over_radius * sin_phi * sin_theta_sq)
    !Hxzl = (cos_phi * cos_theta * minus_dg_plus_g_over_radius * sin_theta)
    !Hyzl = (cos_theta * minus_dg_plus_g_over_radius * sin_phi * sin_theta)
    ! becomes
    !Hxxl = minus_g_over_radius * (cos_phi_sq + sin_phi_sq) = minus_g_over_radius
    !Hyyl = minus_g_over_radius * (cos_phi_sq + sin_phi_sq) = minus_g_over_radius
    !Hzzl = minus_dg
    !Hxyl = 0
    !Hxzl = 0
    !Hyzl = 0
    ! note that the components Hxx == Hyy == - g/r account for the curvature of the Earth
    ! that is if we move horizontally by dx, then the "down" direction tilts slightly which create a change in the x-component
    ! and y-component correspondingly.
    !
    ! flat Earth: assumes g only acts in negative z-direction
    Hxxl = 0._CUSTOM_REAL
    Hyyl = 0._CUSTOM_REAL
    Hzzl = minus_deriv_gravity(iglob)
    Hxyl = 0._CUSTOM_REAL
    Hxzl = 0._CUSTOM_REAL
    Hyzl = 0._CUSTOM_REAL

    ! precompute vector
    if (ispec_irreg /= 0) then
      jacobianl = jacobianstore(INDEX_IJK,ispec)
    else
      jacobianl = jacobian_regular
    endif

    factor = jacobianl * wgll_cube(INDEX_IJK) * rhol

    rho_s_H(1,INDEX_IJK) = 0._CUSTOM_REAL     ! factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
    rho_s_H(2,INDEX_IJK) = 0._CUSTOM_REAL     ! factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
    rho_s_H(3,INDEX_IJK) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)
  ENDDO_LOOP_IJK

  end subroutine compute_element_gravity
