!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

subroutine compute_forces_elastic_noDev( iphase, &
                        NSPEC_AB,NGLOB_AB,displ,accel,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool,&
                        ATTENUATION,&
                        one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
                        NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy,&
                        epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT,&
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic  )

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM, &
                      N_SLS,SAVE_MOHO_MESH, &
                      ONE_THIRD,FOUR_THIRDS

  implicit none

  !include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,accel

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        kappastore,mustore,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! communication overlap
!  logical, dimension(NSPEC_AB) :: ispec_is_inner
!  logical :: phase_is_inner

! memory variables and standard linear solids for attenuation
  logical :: ATTENUATION
  logical :: COMPUTE_AND_STORE_STRAIN
  integer :: NSPEC_STRAIN_ONLY, NSPEC_ADJOINT
  integer :: NSPEC_ATTENUATION_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
       R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: &
       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: epsilon_trace_over_3

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store

!  logical,dimension(NSPEC_AB) :: ispec_is_elastic
  integer :: iphase
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic
  integer, dimension(num_phase_ispec_elastic,2) :: phase_ispec_inner_elastic


! adjoint simulations
  integer :: SIMULATION_TYPE
  integer :: NSPEC_BOUN,NSPEC2D_MOHO

  ! moho kernel
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO):: &
    dsdx_top,dsdx_bot
  logical,dimension(NSPEC_BOUN) :: is_moho_top,is_moho_bot
  integer :: ispec2D_moho_top, ispec2D_moho_bot

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xx_loc, &
       epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL) factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
  real(kind=CUSTOM_REAL) templ

  integer i_SLS
  integer ispec,iglob,ispec_p,num_elements
  integer i,j,k,l


  if( iphase == 1 ) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! adjoint simulations: moho kernel
    ! note: call this only once
    if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
      if (is_moho_top(ispec)) then
        ispec2D_moho_top = ispec2D_moho_top + 1
      else if (is_moho_bot(ispec)) then
        ispec2D_moho_bot = ispec2D_moho_bot + 1
      endif
    endif

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0.
          tempx2l = 0.
          tempx3l = 0.

          tempy1l = 0.
          tempy2l = 0.
          tempy3l = 0.

          tempz1l = 0.
          tempz2l = 0.
          tempz3l = 0.

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1l = tempx1l + displ(1,iglob)*hp1
            tempy1l = tempy1l + displ(2,iglob)*hp1
            tempz1l = tempz1l + displ(3,iglob)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2l = tempx2l + displ(1,iglob)*hp2
            tempy2l = tempy2l + displ(2,iglob)*hp2
            tempz2l = tempz2l + displ(3,iglob)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3l = tempx3l + displ(1,iglob)*hp3
            tempy3l = tempy3l + displ(2,iglob)*hp3
            tempz3l = tempz3l + displ(3,iglob)*hp3
          enddo

          ! get derivatives of ux, uy and uz with respect to x, y and z
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

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          ! adjoint simulations: save strain on the Moho boundary
          if (SAVE_MOHO_MESH ) then
            if (is_moho_top(ispec)) then
              dsdx_top(1,1,i,j,k,ispec2D_moho_top) = duxdxl
              dsdx_top(1,2,i,j,k,ispec2D_moho_top) = duxdyl
              dsdx_top(1,3,i,j,k,ispec2D_moho_top) = duxdzl
              dsdx_top(2,1,i,j,k,ispec2D_moho_top) = duydxl
              dsdx_top(2,2,i,j,k,ispec2D_moho_top) = duydyl
              dsdx_top(2,3,i,j,k,ispec2D_moho_top) = duydzl
              dsdx_top(3,1,i,j,k,ispec2D_moho_top) = duzdxl
              dsdx_top(3,2,i,j,k,ispec2D_moho_top) = duzdyl
              dsdx_top(3,3,i,j,k,ispec2D_moho_top) = duzdzl
            else if (is_moho_bot(ispec)) then
              dsdx_bot(1,1,i,j,k,ispec2D_moho_bot) = duxdxl
              dsdx_bot(1,2,i,j,k,ispec2D_moho_bot) = duxdyl
              dsdx_bot(1,3,i,j,k,ispec2D_moho_bot) = duxdzl
              dsdx_bot(2,1,i,j,k,ispec2D_moho_bot) = duydxl
              dsdx_bot(2,2,i,j,k,ispec2D_moho_bot) = duydyl
              dsdx_bot(2,3,i,j,k,ispec2D_moho_bot) = duydzl
              dsdx_bot(3,1,i,j,k,ispec2D_moho_bot) = duzdxl
              dsdx_bot(3,2,i,j,k,ispec2D_moho_bot) = duzdyl
              dsdx_bot(3,3,i,j,k,ispec2D_moho_bot) = duzdzl
            endif
          endif

          ! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

          ! computes deviatoric strain attenuation and/or for kernel calculations
          if (COMPUTE_AND_STORE_STRAIN) then
            templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
            if( SIMULATION_TYPE == 3 ) epsilon_trace_over_3(i,j,k,ispec) = templ
            epsilondev_xx_loc(i,j,k) = duxdxl - templ
            epsilondev_yy_loc(i,j,k) = duydyl - templ
            epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl
            epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl
            epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl
          endif

          kappal = kappastore(i,j,k,ispec)
          mul = mustore(i,j,k,ispec)

          if(ATTENUATION) then
            ! use unrelaxed parameters if attenuation
            mul  = mul * one_minus_sum_beta(i,j,k,ispec)
          endif

          ! full anisotropic case, stress calculations
          if(ANISOTROPY) then
            c11 = c11store(i,j,k,ispec)
            c12 = c12store(i,j,k,ispec)
            c13 = c13store(i,j,k,ispec)
            c14 = c14store(i,j,k,ispec)
            c15 = c15store(i,j,k,ispec)
            c16 = c16store(i,j,k,ispec)
            c22 = c22store(i,j,k,ispec)
            c23 = c23store(i,j,k,ispec)
            c24 = c24store(i,j,k,ispec)
            c25 = c25store(i,j,k,ispec)
            c26 = c26store(i,j,k,ispec)
            c33 = c33store(i,j,k,ispec)
            c34 = c34store(i,j,k,ispec)
            c35 = c35store(i,j,k,ispec)
            c36 = c36store(i,j,k,ispec)
            c44 = c44store(i,j,k,ispec)
            c45 = c45store(i,j,k,ispec)
            c46 = c46store(i,j,k,ispec)
            c55 = c55store(i,j,k,ispec)
            c56 = c56store(i,j,k,ispec)
            c66 = c66store(i,j,k,ispec)

            sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                      c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl
            sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                      c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl
            sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                      c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl
            sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                      c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl
            sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                      c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl
            sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                      c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

          else

            ! isotropic case
            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2.*mul

            ! compute stress sigma
            sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
            sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
            sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

            sigma_xy = mul*duxdyl_plus_duydxl
            sigma_xz = mul*duzdxl_plus_duxdzl
            sigma_yz = mul*duzdyl_plus_duydzl

          endif ! ANISOTROPY

          ! subtract memory variables if attenuation
          if(ATTENUATION) then
             do i_sls = 1,N_SLS
                R_xx_val = R_xx(i,j,k,ispec,i_sls)
                R_yy_val = R_yy(i,j,k,ispec,i_sls)
                sigma_xx = sigma_xx - R_xx_val
                sigma_yy = sigma_yy - R_yy_val
                sigma_zz = sigma_zz + R_xx_val + R_yy_val
                sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
                sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
                sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
             enddo
          endif

          ! form dot product with test vector, symmetric form
          tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl)
          tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl)
          tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

          tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl)
          tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl)
          tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

          tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl)
          tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl)
          tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

        enddo
      enddo
    enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0.
          tempy1l = 0.
          tempz1l = 0.

          tempx2l = 0.
          tempy2l = 0.
          tempz2l = 0.

          tempx3l = 0.
          tempy3l = 0.
          tempz3l = 0.

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac3 = hprimewgll_zz(l,k)
            tempx3l = tempx3l + tempx3(i,j,l)*fac3
            tempy3l = tempy3l + tempy3(i,j,l)*fac3
            tempz3l = tempz3l + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh
          iglob = ibool(i,j,k,ispec)

          accel(1,iglob) = accel(1,iglob) - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
          accel(2,iglob) = accel(2,iglob) - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
          accel(3,iglob) = accel(3,iglob) - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

          !  update memory variables based upon the Runge-Kutta scheme
          if(ATTENUATION) then

             ! use Runge-Kutta scheme to march in time
             do i_sls = 1,N_SLS

                factor_loc = mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)

                alphaval_loc = alphaval(i_sls)
                betaval_loc = betaval(i_sls)
                gammaval_loc = gammaval(i_sls)

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

             enddo   ! end of loop on memory variables

          endif  !  end attenuation


        enddo
      enddo
    enddo

    ! save deviatoric strain for Runge-Kutta scheme
    if ( COMPUTE_AND_STORE_STRAIN ) then
      epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
    endif

  enddo  ! spectral element loop


end subroutine compute_forces_elastic_noDev

