!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
                        ATTENUATION,USE_OLSEN_ATTENUATION,&
                        one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
                        NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_xx,epsilondev_yy,epsilondev_xy,&
                        epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
                        rho_vs,&
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                        b_displ,b_accel,kappa_kl,mu_kl,deltat, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ATT_AND_KERNEL,&
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot,b_dsdx_top,b_dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                        b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                        b_epsilondev_xz,b_epsilondev_yz, &
                        b_alphaval,b_betaval,b_gammaval,&
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic,&
                        phase_ispec_inner_elastic  )

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM, &
                      NUM_REGIONS_ATTENUATION,N_SLS,SAVE_MOHO_MESH, &
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
  logical :: ATTENUATION,USE_OLSEN_ATTENUATION
  integer :: NSPEC_ATTENUATION_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: iflag_attenuation_store
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: factor_common, alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
       R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: &
       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vs

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
  integer :: NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ATT_AND_KERNEL
  integer :: NGLOB_ADJOINT,NSPEC_ADJOINT

  ! moho kernel
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO):: &
    dsdx_top,dsdx_bot,b_dsdx_top,b_dsdx_bot
  logical,dimension(NSPEC_BOUN) :: is_moho_top,is_moho_bot
  integer :: ispec2D_moho_top, ispec2D_moho_bot
    
  ! adjoint memory variables
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS) :: &
       b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL) :: &
       b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: b_alphaval,b_betaval,b_gammaval

  ! adjoint wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_ADJOINT):: b_displ,b_accel
  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: &
    mu_kl, kappa_kl
  real(kind=CUSTOM_REAL) :: deltat
  
!adjoint

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer ispec,iglob,ispec_p,num_elements
  integer i,j,k,l

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
  real(kind=CUSTOM_REAL) epsilon_trace_over_3
  real(kind=CUSTOM_REAL) vs_val

  integer i_SLS,iselected

  ! adjoint backward arrays
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    b_tempx1,b_tempx2,b_tempx3,b_tempy1,b_tempy2,b_tempy3,b_tempz1,b_tempz2,b_tempz3
  real(kind=CUSTOM_REAL):: dsxx,dsxy,dsxz,dsyy,dsyz,dszz
  real(kind=CUSTOM_REAL):: b_duxdxl,b_duxdyl,b_duxdzl,b_duydxl,b_duydyl,b_duydzl,b_duzdxl,b_duzdyl,b_duzdzl
  real(kind=CUSTOM_REAL):: b_duxdxl_plus_duydyl,b_duxdxl_plus_duzdzl,b_duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL):: b_duxdyl_plus_duydxl,b_duzdxl_plus_duxdzl,b_duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL):: b_dsxx,b_dsxy,b_dsxz,b_dsyy,b_dsyz,b_dszz
  real(kind=CUSTOM_REAL):: b_sigma_xx,b_sigma_yy,b_sigma_zz,b_sigma_xy,b_sigma_xz,b_sigma_yz
  real(kind=CUSTOM_REAL):: kappa_k, mu_k
  real(kind=CUSTOM_REAL) b_tempx1l,b_tempx2l,b_tempx3l
  real(kind=CUSTOM_REAL) b_tempy1l,b_tempy2l,b_tempy3l
  real(kind=CUSTOM_REAL) b_tempz1l,b_tempz2l,b_tempz3l
  ! local adjoint attenuation
  real(kind=CUSTOM_REAL) b_alphaval_loc,b_betaval_loc,b_gammaval_loc,b_Sn,b_Snp1
  real(kind=CUSTOM_REAL) b_epsilon_trace_over_3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: b_epsilondev_xx_loc, &
       b_epsilondev_yy_loc, b_epsilondev_xy_loc, b_epsilondev_xz_loc, b_epsilondev_yz_loc
  real(kind=CUSTOM_REAL) b_R_xx_val,b_R_yy_val  
  ! adjoint

  if( iphase == 1 ) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  do ispec_p = 1,num_elements

!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

!      if( ispec_is_elastic(ispec) ) then

        ispec = phase_ispec_inner_elastic(ispec_p,iphase)


        ! adjoint simulations: moho kernel
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

              if (SIMULATION_TYPE == 3) then
                b_tempx1l = 0.
                b_tempx2l = 0.
                b_tempx3l = 0.

                b_tempy1l = 0.
                b_tempy2l = 0.
                b_tempy3l = 0.

                b_tempz1l = 0.
                b_tempz2l = 0.
                b_tempz3l = 0.
              endif

              do l=1,NGLLX
                hp1 = hprime_xx(i,l)
                iglob = ibool(l,j,k,ispec)
                tempx1l = tempx1l + displ(1,iglob)*hp1
                tempy1l = tempy1l + displ(2,iglob)*hp1
                tempz1l = tempz1l + displ(3,iglob)*hp1
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx1l = b_tempx1l + b_displ(1,iglob)*hp1
                  b_tempy1l = b_tempy1l + b_displ(2,iglob)*hp1
                  b_tempz1l = b_tempz1l + b_displ(3,iglob)*hp1
                endif ! adjoint
    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                hp2 = hprime_yy(j,l)
                iglob = ibool(i,l,k,ispec)
                tempx2l = tempx2l + displ(1,iglob)*hp2
                tempy2l = tempy2l + displ(2,iglob)*hp2
                tempz2l = tempz2l + displ(3,iglob)*hp2
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx2l = b_tempx2l + b_displ(1,iglob)*hp2
                  b_tempy2l = b_tempy2l + b_displ(2,iglob)*hp2
                  b_tempz2l = b_tempz2l + b_displ(3,iglob)*hp2
                endif ! adjoint
    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                hp3 = hprime_zz(k,l)
                iglob = ibool(i,j,l,ispec)
                tempx3l = tempx3l + displ(1,iglob)*hp3
                tempy3l = tempy3l + displ(2,iglob)*hp3
                tempz3l = tempz3l + displ(3,iglob)*hp3
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx3l = b_tempx3l + b_displ(1,iglob)*hp3
                  b_tempy3l = b_tempy3l + b_displ(2,iglob)*hp3
                  b_tempz3l = b_tempz3l + b_displ(3,iglob)*hp3
                endif ! adjoint
                
                
              enddo

    !         get derivatives of ux, uy and uz with respect to x, y and z
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
              if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
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

              kappal = kappastore(i,j,k,ispec)
              mul = mustore(i,j,k,ispec)

              ! adjoint simulations
              if (SIMULATION_TYPE == 3) then
                dsxx = duxdxl
                dsxy = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
                dsxz = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
                dsyy = duydyl
                dsyz = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
                dszz = duzdzl

                b_duxdxl = xixl*b_tempx1l + etaxl*b_tempx2l + gammaxl*b_tempx3l
                b_duxdyl = xiyl*b_tempx1l + etayl*b_tempx2l + gammayl*b_tempx3l
                b_duxdzl = xizl*b_tempx1l + etazl*b_tempx2l + gammazl*b_tempx3l

                b_duydxl = xixl*b_tempy1l + etaxl*b_tempy2l + gammaxl*b_tempy3l
                b_duydyl = xiyl*b_tempy1l + etayl*b_tempy2l + gammayl*b_tempy3l
                b_duydzl = xizl*b_tempy1l + etazl*b_tempy2l + gammazl*b_tempy3l

                b_duzdxl = xixl*b_tempz1l + etaxl*b_tempz2l + gammaxl*b_tempz3l
                b_duzdyl = xiyl*b_tempz1l + etayl*b_tempz2l + gammayl*b_tempz3l
                b_duzdzl = xizl*b_tempz1l + etazl*b_tempz2l + gammazl*b_tempz3l

                b_duxdxl_plus_duydyl = b_duxdxl + b_duydyl
                b_duxdxl_plus_duzdzl = b_duxdxl + b_duzdzl
                b_duydyl_plus_duzdzl = b_duydyl + b_duzdzl
                b_duxdyl_plus_duydxl = b_duxdyl + b_duydxl
                b_duzdxl_plus_duxdzl = b_duzdxl + b_duxdzl
                b_duzdyl_plus_duydzl = b_duzdyl + b_duydzl

                b_dsxx =  b_duxdxl
                b_dsxy = 0.5_CUSTOM_REAL * b_duxdyl_plus_duydxl
                b_dsxz = 0.5_CUSTOM_REAL * b_duzdxl_plus_duxdzl
                b_dsyy =  b_duydyl
                b_dsyz = 0.5_CUSTOM_REAL * b_duzdyl_plus_duydzl
                b_dszz =  b_duzdzl

                ! isotropic adjoint kernels: bulk (kappa) and shear (mu) kernels
                kappa_k = (duxdxl + duydyl + duzdzl) *  (b_duxdxl + b_duydyl + b_duzdzl)
                mu_k = dsxx * b_dsxx + dsyy * b_dsyy + dszz * b_dszz + &
                      2._CUSTOM_REAL * (dsxy * b_dsxy + dsxz * b_dsxz + dsyz * b_dsyz) &
                      - ONE_THIRD * kappa_k
                      
                kappa_kl(i,j,k,ispec) = kappa_kl(i,j,k,ispec) + deltat * kappa_k
                mu_kl(i,j,k,ispec) = mu_kl(i,j,k,ispec) + 2._CUSTOM_REAL * deltat * mu_k

                if (SAVE_MOHO_MESH) then
                  if (is_moho_top(ispec)) then
                    b_dsdx_top(1,1,i,j,k,ispec2D_moho_top) = b_duxdxl
                    b_dsdx_top(1,2,i,j,k,ispec2D_moho_top) = b_duxdyl
                    b_dsdx_top(1,3,i,j,k,ispec2D_moho_top) = b_duxdzl
                    b_dsdx_top(2,1,i,j,k,ispec2D_moho_top) = b_duydxl
                    b_dsdx_top(2,2,i,j,k,ispec2D_moho_top) = b_duydyl
                    b_dsdx_top(2,3,i,j,k,ispec2D_moho_top) = b_duydzl
                    b_dsdx_top(3,1,i,j,k,ispec2D_moho_top) = b_duzdxl
                    b_dsdx_top(3,2,i,j,k,ispec2D_moho_top) = b_duzdyl
                    b_dsdx_top(3,3,i,j,k,ispec2D_moho_top) = b_duzdzl
                  else if (is_moho_bot(ispec)) then
                    b_dsdx_bot(1,1,i,j,k,ispec2D_moho_bot) = b_duxdxl
                    b_dsdx_bot(1,2,i,j,k,ispec2D_moho_bot) = b_duxdyl
                    b_dsdx_bot(1,3,i,j,k,ispec2D_moho_bot) = b_duxdzl
                    b_dsdx_bot(2,1,i,j,k,ispec2D_moho_bot) = b_duydxl
                    b_dsdx_bot(2,2,i,j,k,ispec2D_moho_bot) = b_duydyl
                    b_dsdx_bot(2,3,i,j,k,ispec2D_moho_bot) = b_duydzl
                    b_dsdx_bot(3,1,i,j,k,ispec2D_moho_bot) = b_duzdxl
                    b_dsdx_bot(3,2,i,j,k,ispec2D_moho_bot) = b_duzdyl
                    b_dsdx_bot(3,3,i,j,k,ispec2D_moho_bot) = b_duzdzl
                  endif
                endif
              endif ! adjoint

              if(ATTENUATION) then
                ! compute deviatoric strain
                epsilon_trace_over_3 = ONE_THIRD * (duxdxl + duydyl + duzdzl)
                epsilondev_xx_loc(i,j,k) = duxdxl - epsilon_trace_over_3
                epsilondev_yy_loc(i,j,k) = duydyl - epsilon_trace_over_3
                epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl
                epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl
                epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl
                
                ! adjoint simulations                                
                if (SIMULATION_TYPE == 3) then
                  b_epsilon_trace_over_3 = ONE_THIRD * (b_duxdxl + b_duydyl + b_duzdzl)
                  b_epsilondev_xx_loc(i,j,k) = b_duxdxl - b_epsilon_trace_over_3
                  b_epsilondev_yy_loc(i,j,k) = b_duydyl - b_epsilon_trace_over_3
                  b_epsilondev_xy_loc(i,j,k) = 0.5 * b_duxdyl_plus_duydxl
                  b_epsilondev_xz_loc(i,j,k) = 0.5 * b_duzdxl_plus_duxdzl
                  b_epsilondev_yz_loc(i,j,k) = 0.5 * b_duzdyl_plus_duydzl
                endif ! adjoint

                ! uses scaling rule similar to Olsen et al. (2003) or mesh flag
                if(USE_OLSEN_ATTENUATION) then
                  vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
                  call get_attenuation_model_olsen( vs_val, iselected )
                else
                  ! iflag from (CUBIT) mesh      
                  iselected = iflag_attenuation_store(i,j,k,ispec)                
                endif

                ! use unrelaxed parameters if attenuation
                mul = mul * one_minus_sum_beta(iselected)
                 
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
                !if(ATTENUATION .and. not_fully_in_bedrock(ispec)) then
                !   mul = c44
                !   c11 = c11 + FOUR_THIRDS * minus_sum_beta * mul
                !   c12 = c12 - TWO_THIRDS * minus_sum_beta * mul
                !   c13 = c13 - TWO_THIRDS * minus_sum_beta * mul
                !   c22 = c22 + FOUR_THIRDS * minus_sum_beta * mul
                !   c23 = c23 - TWO_THIRDS * minus_sum_beta * mul
                !   c33 = c33 + FOUR_THIRDS * minus_sum_beta * mul
                !   c44 = c44 + minus_sum_beta * mul
                !   c55 = c55 + minus_sum_beta * mul
                !   c66 = c66 + minus_sum_beta * mul
                !endif

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

                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_sigma_xx = c11*b_duxdxl + c16*b_duxdyl_plus_duydxl + c12*b_duydyl + &
                       c15*b_duzdxl_plus_duxdzl + c14*b_duzdyl_plus_duydzl + c13*b_duzdzl
                  b_sigma_yy = c12*b_duxdxl + c26*b_duxdyl_plus_duydxl + c22*b_duydyl + &
                       c25*b_duzdxl_plus_duxdzl + c24*b_duzdyl_plus_duydzl + c23*b_duzdzl
                  b_sigma_zz = c13*b_duxdxl + c36*b_duxdyl_plus_duydxl + c23*b_duydyl + &
                       c35*b_duzdxl_plus_duxdzl + c34*b_duzdyl_plus_duydzl + c33*b_duzdzl
                  b_sigma_xy = c16*b_duxdxl + c66*b_duxdyl_plus_duydxl + c26*b_duydyl + &
                       c56*b_duzdxl_plus_duxdzl + c46*b_duzdyl_plus_duydzl + c36*b_duzdzl
                  b_sigma_xz = c15*b_duxdxl + c56*b_duxdyl_plus_duydxl + c25*b_duydyl + &
                       c55*b_duzdxl_plus_duxdzl + c45*b_duzdyl_plus_duydzl + c35*b_duzdzl
                  b_sigma_yz = c14*b_duxdxl + c46*b_duxdyl_plus_duydxl + c24*b_duydyl + &
                       c45*b_duzdxl_plus_duxdzl + c44*b_duzdyl_plus_duydzl + c34*b_duzdzl
                endif ! adjoint
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

                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_sigma_xx = lambdalplus2mul*b_duxdxl + lambdal*b_duydyl_plus_duzdzl
                  b_sigma_yy = lambdalplus2mul*b_duydyl + lambdal*b_duxdxl_plus_duzdzl
                  b_sigma_zz = lambdalplus2mul*b_duzdzl + lambdal*b_duxdxl_plus_duydyl                
                  b_sigma_xy = mul*b_duxdyl_plus_duydxl
                  b_sigma_xz = mul*b_duzdxl_plus_duxdzl
                  b_sigma_yz = mul*b_duzdyl_plus_duydzl
                endif !adjoint

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

                    ! adjoint simulations
                    if (SIMULATION_TYPE == 3) then
                      b_R_xx_val = b_R_xx(i,j,k,ispec,i_sls)
                      b_R_yy_val = b_R_yy(i,j,k,ispec,i_sls)
                      b_sigma_xx = b_sigma_xx - b_R_xx_val
                      b_sigma_yy = b_sigma_yy - b_R_yy_val
                      b_sigma_zz = b_sigma_zz + b_R_xx_val + b_R_yy_val
                      b_sigma_xy = b_sigma_xy - b_R_xy(i,j,k,ispec,i_sls)
                      b_sigma_xz = b_sigma_xz - b_R_xz(i,j,k,ispec,i_sls)
                      b_sigma_yz = b_sigma_yz - b_R_yz(i,j,k,ispec,i_sls)
                    endif !adjoint                    
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

              ! adjoint simulations
              if (SIMULATION_TYPE == 3) then
                b_tempx1(i,j,k) = jacobianl * (b_sigma_xx*xixl + b_sigma_xy*xiyl + b_sigma_xz*xizl)
                b_tempy1(i,j,k) = jacobianl * (b_sigma_xy*xixl + b_sigma_yy*xiyl + b_sigma_yz*xizl)
                b_tempz1(i,j,k) = jacobianl * (b_sigma_xz*xixl + b_sigma_yz*xiyl + b_sigma_zz*xizl)
                b_tempx2(i,j,k) = jacobianl * (b_sigma_xx*etaxl + b_sigma_xy*etayl + b_sigma_xz*etazl)
                b_tempy2(i,j,k) = jacobianl * (b_sigma_xy*etaxl + b_sigma_yy*etayl + b_sigma_yz*etazl)
                b_tempz2(i,j,k) = jacobianl * (b_sigma_xz*etaxl + b_sigma_yz*etayl + b_sigma_zz*etazl)
                b_tempx3(i,j,k) = jacobianl * (b_sigma_xx*gammaxl + b_sigma_xy*gammayl + b_sigma_xz*gammazl)
                b_tempy3(i,j,k) = jacobianl * (b_sigma_xy*gammaxl + b_sigma_yy*gammayl + b_sigma_yz*gammazl)
                b_tempz3(i,j,k) = jacobianl * (b_sigma_xz*gammaxl + b_sigma_yz*gammayl + b_sigma_zz*gammazl)
              endif !adjoint

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

              ! adjoint simulations
              if (SIMULATION_TYPE == 3) then
                b_tempx1l = 0.
                b_tempy1l = 0.
                b_tempz1l = 0.
                b_tempx2l = 0.
                b_tempy2l = 0.
                b_tempz2l = 0.
                b_tempx3l = 0.
                b_tempy3l = 0.
                b_tempz3l = 0.
              endif !adjoint

              do l=1,NGLLX
                fac1 = hprimewgll_xx(l,i)
                tempx1l = tempx1l + tempx1(l,j,k)*fac1
                tempy1l = tempy1l + tempy1(l,j,k)*fac1
                tempz1l = tempz1l + tempz1(l,j,k)*fac1
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx1l = b_tempx1l + b_tempx1(l,j,k)*fac1
                  b_tempy1l = b_tempy1l + b_tempy1(l,j,k)*fac1
                  b_tempz1l = b_tempz1l + b_tempz1(l,j,k)*fac1
                endif                
                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                fac2 = hprimewgll_yy(l,j)
                tempx2l = tempx2l + tempx2(i,l,k)*fac2
                tempy2l = tempy2l + tempy2(i,l,k)*fac2
                tempz2l = tempz2l + tempz2(i,l,k)*fac2
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx2l = b_tempx2l + b_tempx2(i,l,k)*fac2
                  b_tempy2l = b_tempy2l + b_tempy2(i,l,k)*fac2
                  b_tempz2l = b_tempz2l + b_tempz2(i,l,k)*fac2
                endif                
                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                fac3 = hprimewgll_zz(l,k)
                tempx3l = tempx3l + tempx3(i,j,l)*fac3
                tempy3l = tempy3l + tempy3(i,j,l)*fac3
                tempz3l = tempz3l + tempz3(i,j,l)*fac3
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  b_tempx3l = b_tempx3l + b_tempx3(i,j,l)*fac3
                  b_tempy3l = b_tempy3l + b_tempy3(i,j,l)*fac3
                  b_tempz3l = b_tempz3l + b_tempz3(i,j,l)*fac3
                endif                
              enddo

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)

    ! sum contributions from each element to the global mesh

              iglob = ibool(i,j,k,ispec)

              accel(1,iglob) = accel(1,iglob) - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
              accel(2,iglob) = accel(2,iglob) - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
              accel(3,iglob) = accel(3,iglob) - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

              ! adjoint simulations
              if (SIMULATION_TYPE == 3) then
                b_accel(1,iglob) = b_accel(1,iglob) - (fac1*b_tempx1l + fac2*b_tempx2l + fac3*b_tempx3l)
                b_accel(2,iglob) = b_accel(2,iglob) - (fac1*b_tempy1l + fac2*b_tempy2l + fac3*b_tempy3l)
                b_accel(3,iglob) = b_accel(3,iglob) - (fac1*b_tempz1l + fac2*b_tempz2l + fac3*b_tempz3l)
              endif !adjoint

              !  update memory variables based upon the Runge-Kutta scheme
              if(ATTENUATION) then
                 
                 ! use Runge-Kutta scheme to march in time
                 do i_sls = 1,N_SLS

                    ! get coefficients for that standard linear solid
                    if( USE_OLSEN_ATTENUATION ) then
                      vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
                      call get_attenuation_model_olsen( vs_val, iselected )
                    else
                      iselected = iflag_attenuation_store(i,j,k,ispec)
                    endif
                    
                    factor_loc = mustore(i,j,k,ispec) * factor_common(iselected,i_sls)
                    
                    alphaval_loc = alphaval(iselected,i_sls)
                    betaval_loc = betaval(iselected,i_sls)
                    gammaval_loc = gammaval(iselected,i_sls)
                    
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
                    
                    !adjoint simulations
                    if (SIMULATION_TYPE == 3) then
                      b_alphaval_loc = b_alphaval(iselected,i_sls)
                      b_betaval_loc = b_betaval(iselected,i_sls)
                      b_gammaval_loc = b_gammaval(iselected,i_sls)
                      ! term in xx
                      b_Sn   = factor_loc * b_epsilondev_xx(i,j,k,ispec)
                      b_Snp1   = factor_loc * b_epsilondev_xx_loc(i,j,k)
                      b_R_xx(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xx(i,j,k,ispec,i_sls) + &
                                            b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                      ! term in yy
                      b_Sn   = factor_loc * b_epsilondev_yy(i,j,k,ispec)
                      b_Snp1   = factor_loc * b_epsilondev_yy_loc(i,j,k)
                      b_R_yy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yy(i,j,k,ispec,i_sls) + &
                                            b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                      ! term in zz not computed since zero trace
                      ! term in xy
                      b_Sn   = factor_loc * b_epsilondev_xy(i,j,k,ispec)
                      b_Snp1   = factor_loc * b_epsilondev_xy_loc(i,j,k)
                      b_R_xy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xy(i,j,k,ispec,i_sls) + &
                                            b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                      ! term in xz
                      b_Sn   = factor_loc * b_epsilondev_xz(i,j,k,ispec)
                      b_Snp1   = factor_loc * b_epsilondev_xz_loc(i,j,k)
                      b_R_xz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xz(i,j,k,ispec,i_sls) + &
                                            b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                      ! term in yz
                      b_Sn   = factor_loc * b_epsilondev_yz(i,j,k,ispec)
                      b_Snp1   = factor_loc * b_epsilondev_yz_loc(i,j,k)
                      b_R_yz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yz(i,j,k,ispec,i_sls) + &
                                            b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    endif !adjoint

                 enddo   ! end of loop on memory variables

              endif  !  end attenuation


            enddo
          enddo
        enddo

        ! save deviatoric strain for Runge-Kutta scheme
        if(ATTENUATION) then
          epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
          epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
          epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
          epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
          epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
          ! adjoint simulations
          if (SIMULATION_TYPE == 3) then
            b_epsilondev_xx(:,:,:,ispec) = b_epsilondev_xx_loc(:,:,:)
            b_epsilondev_yy(:,:,:,ispec) = b_epsilondev_yy_loc(:,:,:)
            b_epsilondev_xy(:,:,:,ispec) = b_epsilondev_xy_loc(:,:,:)
            b_epsilondev_xz(:,:,:,ispec) = b_epsilondev_xz_loc(:,:,:)
            b_epsilondev_yz(:,:,:,ispec) = b_epsilondev_yz_loc(:,:,:)
          endif !adjoint
        endif
!      endif ! ispec_is_elastic
!    endif ! if (ispec_is_inner(ispec) .eqv. phase_is_inner)

  enddo  ! spectral element loop

! forces in elastic media calculated in compute_forces_elastic...
!! adding source
!  do isource = 1,NSOURCES
!
!  if (ispec_is_inner(ispec_selected_source(isource)) .eqv. phase_is_inner) then
!
!  if(USE_FORCE_POINT_SOURCE) then
!
!!   add the source (only if this proc carries the source)
!    if(myrank == islice_selected_source(isource)) then
!
!      iglob = ibool(nint(xi_source(isource)), &
!           nint(eta_source(isource)), &
!           nint(gamma_source(isource)), &
!           ispec_selected_source(isource))
!      f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
!      t0 = 1.2d0/f0
!
!  if (it == 1 .and. myrank == 0) then
!    print *,'using a source of dominant frequency ',f0
!    print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
!    print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
!  endif
!
!      ! we use nu_source(:,3) here because we want a source normal to the surface.
!      ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
!      !accel(:,iglob) = accel(:,iglob) + &
!      !     sngl(nu_source(:,3,isource) * 10000000.d0 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
!      !     exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
!    accel(:,iglob) = accel(:,iglob) + &
!           sngl(nu_source(:,3,isource) * 1.d10 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
!           exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
!
!    endif
!  endif
!
!  endif
!
!  enddo

end subroutine compute_forces_elastic_noDev

