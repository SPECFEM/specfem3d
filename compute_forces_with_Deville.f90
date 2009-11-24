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


subroutine compute_forces_with_Deville( phase_is_inner ,NSPEC_AB,NGLOB_AB, &
                                    displ,accel, &
                                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                    hprime_xx,hprime_xxT, &
                                    hprimewgll_xx,hprimewgll_xxT, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    kappastore,mustore,jacobian,ibool, &
                                    ispec_is_inner, &
                                    ATTENUATION,USE_OLSEN_ATTENUATION, &
                                    one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
                                    NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                    epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                    epsilondev_xz,epsilondev_yz,iflag_attenuation_store, &
                                    rho_vs, &
                                    ANISOTROPY,NSPEC_ANISO, &
                                    c11store,c12store,c13store,c14store,c15store,c16store,&
                                    c22store,c23store,c24store,c25store,c26store,c33store,&
                                    c34store,c35store,c36store,c44store,c45store,c46store,&
                                    c55store,c56store,c66store, &
                                    ispec_is_elastic )

! computes elastic tensor term

  implicit none

  include "constants.h"

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

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

  logical,dimension(NSPEC_AB) :: ispec_is_elastic

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points,E2_m1_m2_5points,E3_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)
  equivalence(newtempy1,E2_m1_m2_5points)
  equivalence(newtempz1,E3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)
  equivalence(newtempy3,E2_mxm_m2_m1_5points)
  equivalence(newtempz3,E3_mxm_m2_m1_5points)

! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xx_loc, &
       epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL) factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
  real(kind=CUSTOM_REAL) epsilon_trace_over_3
  real(kind=CUSTOM_REAL) vs_val

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal

! local anisotropy parameters
  real(kind=CUSTOM_REAL) c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  
  integer i_SLS,iselected

  integer ispec,iglob
  integer i,j,k
  
! loops over all elements
  do ispec = 1,NSPEC_AB

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      if( ispec_is_elastic(ispec) ) then
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                dummyx_loc(i,j,k) = displ(1,iglob)
                dummyy_loc(i,j,k) = displ(2,iglob)
                dummyz_loc(i,j,k) = displ(3,iglob)
            enddo
          enddo
        enddo

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    ! call mxm_m1_m2_5points(hprime_xx,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1)
        do j=1,m2
          do i=1,m1
            C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                                  hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                                  hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                                  hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                                  hprime_xx(i,5)*B1_m1_m2_5points(5,j)

            C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                                  hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                                  hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                                  hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                                  hprime_xx(i,5)*B2_m1_m2_5points(5,j)

            C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                                  hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                                  hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                                  hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                                  hprime_xx(i,5)*B3_m1_m2_5points(5,j)
          enddo
        enddo

    !   call mxm_m1_m1_5points(dummyx_loc(1,1,k),dummyy_loc(1,1,k),dummyz_loc(1,1,k), &
    !          hprime_xxT,tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k))
        do j=1,m1
          do i=1,m1
    ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
            do k = 1,NGLLX
              tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                            dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                            dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                            dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                            dummyx_loc(i,5,k)*hprime_xxT(5,j)

              tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                            dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                            dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                            dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                            dummyy_loc(i,5,k)*hprime_xxT(5,j)

              tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                            dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                            dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                            dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                            dummyz_loc(i,5,k)*hprime_xxT(5,j)
            enddo
          enddo
        enddo

    ! call mxm_m2_m1_5points(dummyx_loc,dummyy_loc,dummyz_loc,tempx3,tempy3,tempz3)
        do j=1,m1
          do i=1,m2
            C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                      A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                      A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                      A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                      A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

            C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                      A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                      A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                      A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                      A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

            C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                      A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                      A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                      A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                      A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
          enddo
        enddo

        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX

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

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

  ! precompute some sums to save CPU time
              duxdxl_plus_duydyl = duxdxl + duydyl
              duxdxl_plus_duzdzl = duxdxl + duzdzl
              duydyl_plus_duzdzl = duydyl + duzdzl
              duxdyl_plus_duydxl = duxdyl + duydxl
              duzdxl_plus_duxdzl = duzdxl + duxdzl
              duzdyl_plus_duydzl = duzdyl + duydzl

              kappal = kappastore(i,j,k,ispec)
              mul = mustore(i,j,k,ispec)

  ! attenuation           
              if(ATTENUATION) then
                ! compute deviatoric strain
                epsilon_trace_over_3 = ONE_THIRD * (duxdxl + duydyl + duzdzl)
                epsilondev_xx_loc(i,j,k) = duxdxl - epsilon_trace_over_3
                epsilondev_yy_loc(i,j,k) = duydyl - epsilon_trace_over_3
                epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl
                epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl
                epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl
                                
                !if (SIMULATION_TYPE == 3) then
                ! b_epsilon_trace_over_3 = ONE_THIRD * (b_duxdxl + b_duydyl + b_duzdzl)
                ! b_epsilondev_xx_loc(i,j,k) = b_duxdxl - b_epsilon_trace_over_3
                ! b_epsilondev_yy_loc(i,j,k) = b_duydyl - b_epsilon_trace_over_3
                ! b_epsilondev_xy_loc(i,j,k) = 0.5 * b_duxdyl_plus_duydxl
                ! b_epsilondev_xz_loc(i,j,k) = 0.5 * b_duzdxl_plus_duxdzl
                ! b_epsilondev_yz_loc(i,j,k) = 0.5 * b_duzdyl_plus_duydzl
                !endif

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

                !if (SIMULATION_TYPE == 3) then
                ! b_sigma_xx = c11*b_duxdxl + c16*b_duxdyl_plus_duydxl + c12*b_duydyl + &
                !       c15*b_duzdxl_plus_duxdzl + c14*b_duzdyl_plus_duydzl + c13*b_duzdzl
                ! b_sigma_yy = c12*b_duxdxl + c26*b_duxdyl_plus_duydxl + c22*b_duydyl + &
                !       c25*b_duzdxl_plus_duxdzl + c24*b_duzdyl_plus_duydzl + c23*b_duzdzl
                ! b_sigma_zz = c13*b_duxdxl + c36*b_duxdyl_plus_duydxl + c23*b_duydyl + &
                !       c35*b_duzdxl_plus_duxdzl + c34*b_duzdyl_plus_duydzl + c33*b_duzdzl
                ! b_sigma_xy = c16*b_duxdxl + c66*b_duxdyl_plus_duydxl + c26*b_duydyl + &
                !       c56*b_duzdxl_plus_duxdzl + c46*b_duzdyl_plus_duydzl + c36*b_duzdzl
                ! b_sigma_xz = c15*b_duxdxl + c56*b_duxdyl_plus_duydxl + c25*b_duydyl + &
                !       c55*b_duzdxl_plus_duxdzl + c45*b_duzdyl_plus_duydzl + c35*b_duzdzl
                ! b_sigma_yz = c14*b_duxdxl + c46*b_duxdyl_plus_duydxl + c24*b_duydyl + &
                !       c45*b_duzdxl_plus_duxdzl + c44*b_duzdyl_plus_duydzl + c34*b_duzdzl
                !endif
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

                !if (SIMULATION_TYPE == 3) then
                ! b_sigma_xx = lambdalplus2mul*b_duxdxl + lambdal*b_duydyl_plus_duzdzl
                ! b_sigma_yy = lambdalplus2mul*b_duydyl + lambdal*b_duxdxl_plus_duzdzl
                ! b_sigma_zz = lambdalplus2mul*b_duzdzl + lambdal*b_duxdxl_plus_duydyl
                !
                ! b_sigma_xy = mul*b_duxdyl_plus_duydxl
                ! b_sigma_xz = mul*b_duzdxl_plus_duxdzl
                ! b_sigma_yz = mul*b_duzdyl_plus_duydzl
                !endif

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
               
                !if (SIMULATION_TYPE == 3) then
                !  b_R_xx_val = b_R_xx(i,j,k,ispec,i_sls)
                !  b_R_yy_val = b_R_yy(i,j,k,ispec,i_sls)
                !  b_sigma_xx = b_sigma_xx - b_R_xx_val
                !  b_sigma_yy = b_sigma_yy - b_R_yy_val
                !  b_sigma_zz = b_sigma_zz + b_R_xx_val + b_R_yy_val
                !  b_sigma_xy = b_sigma_xy - b_R_xy(i,j,k,ispec,i_sls)
                !  b_sigma_xz = b_sigma_xz - b_R_xz(i,j,k,ispec,i_sls)
                !  b_sigma_yz = b_sigma_yz - b_R_yz(i,j,k,ispec,i_sls)
                !endif
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

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    ! call mxm_m1_m2_5points(hprimewgll_xxT,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1)
        do j=1,m2
          do i=1,m1
            E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
                                  hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
                                  hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
                                  hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
                                  hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)

            E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
                                  hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
                                  hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
                                  hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
                                  hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)

            E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
                                  hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
                                  hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
                                  hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
                                  hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
          enddo
        enddo

    !   call mxm_m1_m1_5points(tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k), &
    !         hprimewgll_xx,newtempx2(1,1,k),newtempy2(1,1,k),newtempz2(1,1,k))
        do i=1,m1
          do j=1,m1
    ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
            do k = 1,NGLLX
              newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                               tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                               tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                               tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                               tempx2(i,5,k)*hprimewgll_xx(5,j)

              newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
                               tempy2(i,2,k)*hprimewgll_xx(2,j) + &
                               tempy2(i,3,k)*hprimewgll_xx(3,j) + &
                               tempy2(i,4,k)*hprimewgll_xx(4,j) + &
                               tempy2(i,5,k)*hprimewgll_xx(5,j)

              newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
                               tempz2(i,2,k)*hprimewgll_xx(2,j) + &
                               tempz2(i,3,k)*hprimewgll_xx(3,j) + &
                               tempz2(i,4,k)*hprimewgll_xx(4,j) + &
                               tempz2(i,5,k)*hprimewgll_xx(5,j)
            enddo
          enddo
        enddo

    ! call mxm_m2_m1_5points(tempx3,tempy3,tempz3,hprimewgll_xx,newtempx3,newtempy3,newtempz3)
        do j=1,m1
          do i=1,m2
            E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                      C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                      C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                      C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                      C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

            E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                      C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                      C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                      C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                      C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

            E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                      C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                      C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                      C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                      C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
          enddo
        enddo

        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)

    ! sum contributions from each element to the global mesh using indirect addressing
              iglob = ibool(i,j,k,ispec)
              accel(1,iglob) = accel(1,iglob) - fac1*newtempx1(i,j,k) - &
                                fac2*newtempx2(i,j,k) - fac3*newtempx3(i,j,k)
              accel(2,iglob) = accel(2,iglob) - fac1*newtempy1(i,j,k) - &
                                fac2*newtempy2(i,j,k) - fac3*newtempy3(i,j,k)
              accel(3,iglob) = accel(3,iglob) - fac1*newtempz1(i,j,k) - &
                                fac2*newtempz2(i,j,k) - fac3*newtempz3(i,j,k)

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
                    
                    !if (SIMULATION_TYPE == 3) then
                    !  b_alphaval_loc = b_alphaval(iselected,i_sls)
                    !  b_betaval_loc = b_betaval(iselected,i_sls)
                    !  b_gammaval_loc = b_gammaval(iselected,i_sls)
                    !  ! term in xx
                    !  b_Sn   = factor_loc * b_epsilondev_xx(i,j,k,ispec)
                    !  b_Snp1   = factor_loc * b_epsilondev_xx_loc(i,j,k)
                    !  b_R_xx(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xx(i,j,k,ispec,i_sls) + &
                    !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    !  ! term in yy
                    !  b_Sn   = factor_loc * b_epsilondev_yy(i,j,k,ispec)
                    !  b_Snp1   = factor_loc * b_epsilondev_yy_loc(i,j,k)
                    !  b_R_yy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yy(i,j,k,ispec,i_sls) + &
                    !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    !  ! term in zz not computed since zero trace
                    !  ! term in xy
                    !  b_Sn   = factor_loc * b_epsilondev_xy(i,j,k,ispec)
                    !  b_Snp1   = factor_loc * b_epsilondev_xy_loc(i,j,k)
                    !  b_R_xy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xy(i,j,k,ispec,i_sls) + &
                    !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    !  ! term in xz
                    !  b_Sn   = factor_loc * b_epsilondev_xz(i,j,k,ispec)
                    !  b_Snp1   = factor_loc * b_epsilondev_xz_loc(i,j,k)
                    !  b_R_xz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xz(i,j,k,ispec,i_sls) + &
                    !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    !  ! term in yz
                    !  b_Sn   = factor_loc * b_epsilondev_yz(i,j,k,ispec)
                    !  b_Snp1   = factor_loc * b_epsilondev_yz_loc(i,j,k)
                    !  b_R_yz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yz(i,j,k,ispec,i_sls) + &
                    !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
                    !endif

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
          !if (SIMULATION_TYPE == 3) then
          !  b_epsilondev_xx(:,:,:,ispec) = b_epsilondev_xx_loc(:,:,:)
          !  b_epsilondev_yy(:,:,:,ispec) = b_epsilondev_yy_loc(:,:,:)
          !  b_epsilondev_xy(:,:,:,ispec) = b_epsilondev_xy_loc(:,:,:)
          !  b_epsilondev_xz(:,:,:,ispec) = b_epsilondev_xz_loc(:,:,:)
          !  b_epsilondev_yz(:,:,:,ispec) = b_epsilondev_yz_loc(:,:,:)
          !endif         
        endif

      endif ! ispec_is_elastic
      
    endif ! if (ispec_is_inner(ispec) .eqv. phase_is_inner)

  enddo  ! spectral element loop

end subroutine compute_forces_with_Deville

!
!-------------------------------------------------------------------------------------------------
!

!subroutine compute_forces_with_Deville_noanisotropy( phase_is_inner ,NSPEC_AB,NGLOB_AB, &
!                                    displ,accel, &
!                                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!                                    hprime_xx,hprime_xxT, &
!                                    hprimewgll_xx,hprimewgll_xxT, &
!                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
!                                    kappastore,mustore,jacobian,ibool, &
!                                    ispec_is_inner, &
!                                    ATTENUATION,USE_OLSEN_ATTENUATION, &
!                                    one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
!                                    NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
!                                    epsilondev_xx,epsilondev_yy,epsilondev_xy, &
!                                    epsilondev_xz,epsilondev_yz,iflag_attenuation_store, &
!                                    rho_vs )
!
!! computes elastic tensor term
!
!  implicit none
!
!  include "constants.h"
!!  include values created by the mesher
!!  include "OUTPUT_FILES/values_from_mesher.h"
!
!  integer :: NSPEC_AB,NGLOB_AB
!
!! displacement and acceleration
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,accel
!
!! arrays with mesh parameters per slice
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
!        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
!        kappastore,mustore,jacobian
!
!! array with derivatives of Lagrange polynomials and precalculated products
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz
!
!! communication overlap
!  logical, dimension(NSPEC_AB) :: ispec_is_inner
!  logical :: phase_is_inner
!
!! memory variables and standard linear solids for attenuation    
!  logical :: ATTENUATION,USE_OLSEN_ATTENUATION
!  integer :: NSPEC_ATTENUATION_AB
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: iflag_attenuation_store
!  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: one_minus_sum_beta
!  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: factor_common, alphaval,betaval,gammaval
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
!       R_xx,R_yy,R_xy,R_xz,R_yz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: &
!       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vs
!
!! local parameters
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
!    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
!    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
!
!! manually inline the calls to the Deville et al. (2002) routines
!  real(kind=4), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
!  real(kind=4), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
!  real(kind=4), dimension(m1,m2) :: E1_m1_m2_5points,E2_m1_m2_5points,E3_m1_m2_5points
!
!  equivalence(dummyx_loc,B1_m1_m2_5points)
!  equivalence(dummyy_loc,B2_m1_m2_5points)
!  equivalence(dummyz_loc,B3_m1_m2_5points)
!  equivalence(tempx1,C1_m1_m2_5points)
!  equivalence(tempy1,C2_m1_m2_5points)
!  equivalence(tempz1,C3_m1_m2_5points)
!  equivalence(newtempx1,E1_m1_m2_5points)
!  equivalence(newtempy1,E2_m1_m2_5points)
!  equivalence(newtempz1,E3_m1_m2_5points)
!
!  real(kind=4), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
!  real(kind=4), dimension(m2,m1) :: C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
!  real(kind=4), dimension(m2,m1) :: E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points
!
!  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
!  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
!  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
!  equivalence(tempx3,C1_mxm_m2_m1_5points)
!  equivalence(tempy3,C2_mxm_m2_m1_5points)
!  equivalence(tempz3,C3_mxm_m2_m1_5points)
!  equivalence(newtempx3,E1_mxm_m2_m1_5points)
!  equivalence(newtempy3,E2_mxm_m2_m1_5points)
!  equivalence(newtempz3,E3_mxm_m2_m1_5points)
!
!! local attenuation parameters
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xx_loc, &
!       epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
!  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
!  real(kind=CUSTOM_REAL) factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
!  real(kind=CUSTOM_REAL) epsilon_trace_over_3
!  real(kind=CUSTOM_REAL) vs_val
!
!  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
!  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
!
!  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
!  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
!
!  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz
!
!  real(kind=CUSTOM_REAL) fac1,fac2,fac3
!
!  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
!  real(kind=CUSTOM_REAL) kappal
!  
!  integer i_SLS,iselected
!
!  integer ispec,iglob
!  integer i,j,k
!  
!! loops over all elements
!  do ispec = 1,NSPEC_AB
!
!    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
!
!      do k=1,NGLLZ
!        do j=1,NGLLY
!          do i=1,NGLLX
!              iglob = ibool(i,j,k,ispec)
!              dummyx_loc(i,j,k) = displ(1,iglob)
!              dummyy_loc(i,j,k) = displ(2,iglob)
!              dummyz_loc(i,j,k) = displ(3,iglob)
!          enddo
!        enddo
!      enddo
!
!  ! subroutines adapted from Deville, Fischer and Mund, High-order methods
!  ! for incompressible fluid flow, Cambridge University Press (2002),
!  ! pages 386 and 389 and Figure 8.3.1
!  ! call mxm_m1_m2_5points(hprime_xx,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1)
!      do j=1,m2
!        do i=1,m1
!          C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
!                                hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
!                                hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
!                                hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
!                                hprime_xx(i,5)*B1_m1_m2_5points(5,j)
!
!          C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
!                                hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
!                                hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
!                                hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
!                                hprime_xx(i,5)*B2_m1_m2_5points(5,j)
!
!          C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
!                                hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
!                                hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
!                                hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
!                                hprime_xx(i,5)*B3_m1_m2_5points(5,j)
!        enddo
!      enddo
!
!  !   call mxm_m1_m1_5points(dummyx_loc(1,1,k),dummyy_loc(1,1,k),dummyz_loc(1,1,k), &
!  !          hprime_xxT,tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k))
!      do j=1,m1
!        do i=1,m1
!  ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
!          do k = 1,NGLLX
!            tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
!                          dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
!                          dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
!                          dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
!                          dummyx_loc(i,5,k)*hprime_xxT(5,j)
!
!            tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
!                          dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
!                          dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
!                          dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
!                          dummyy_loc(i,5,k)*hprime_xxT(5,j)
!
!            tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
!                          dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
!                          dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
!                          dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
!                          dummyz_loc(i,5,k)*hprime_xxT(5,j)
!          enddo
!        enddo
!      enddo
!
!  ! call mxm_m2_m1_5points(dummyx_loc,dummyy_loc,dummyz_loc,tempx3,tempy3,tempz3)
!      do j=1,m1
!        do i=1,m2
!          C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                                    A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                                    A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                                    A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                                    A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!
!          C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                                    A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                                    A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                                    A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                                    A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!
!          C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                                    A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                                    A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                                    A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                                    A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!        enddo
!      enddo
!
!      do k=1,NGLLZ
!        do j=1,NGLLY
!          do i=1,NGLLX
!
!  !         get derivatives of ux, uy and uz with respect to x, y and z
!            xixl = xix(i,j,k,ispec)
!            xiyl = xiy(i,j,k,ispec)
!            xizl = xiz(i,j,k,ispec)
!            etaxl = etax(i,j,k,ispec)
!            etayl = etay(i,j,k,ispec)
!            etazl = etaz(i,j,k,ispec)
!            gammaxl = gammax(i,j,k,ispec)
!            gammayl = gammay(i,j,k,ispec)
!            gammazl = gammaz(i,j,k,ispec)
!            jacobianl = jacobian(i,j,k,ispec)
!
!            duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
!            duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
!            duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)
!
!            duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
!            duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
!            duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)
!
!            duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
!            duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
!            duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)
!
!! precompute some sums to save CPU time
!            duxdxl_plus_duydyl = duxdxl + duydyl
!            duxdxl_plus_duzdzl = duxdxl + duzdzl
!            duydyl_plus_duzdzl = duydyl + duzdzl
!            duxdyl_plus_duydxl = duxdyl + duydxl
!            duzdxl_plus_duxdzl = duzdxl + duxdzl
!            duzdyl_plus_duydzl = duzdyl + duydzl
!
!            kappal = kappastore(i,j,k,ispec)
!            mul = mustore(i,j,k,ispec)
!
!! attenuation           
!            if(ATTENUATION) then
!              ! compute deviatoric strain
!              epsilon_trace_over_3 = ONE_THIRD * (duxdxl + duydyl + duzdzl)
!              epsilondev_xx_loc(i,j,k) = duxdxl - epsilon_trace_over_3
!              epsilondev_yy_loc(i,j,k) = duydyl - epsilon_trace_over_3
!              epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl
!              epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl
!              epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl
!                              
!              !if (SIMULATION_TYPE == 3) then
!              ! b_epsilon_trace_over_3 = ONE_THIRD * (b_duxdxl + b_duydyl + b_duzdzl)
!              ! b_epsilondev_xx_loc(i,j,k) = b_duxdxl - b_epsilon_trace_over_3
!              ! b_epsilondev_yy_loc(i,j,k) = b_duydyl - b_epsilon_trace_over_3
!              ! b_epsilondev_xy_loc(i,j,k) = 0.5 * b_duxdyl_plus_duydxl
!              ! b_epsilondev_xz_loc(i,j,k) = 0.5 * b_duzdxl_plus_duxdzl
!              ! b_epsilondev_yz_loc(i,j,k) = 0.5 * b_duzdyl_plus_duydzl
!              !endif
!
!              ! uses scaling rule similar to Olsen et al. (2003) or mesh flag
!              if(USE_OLSEN_ATTENUATION) then
!                vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
!                call get_attenuation_model_olsen( vs_val, iselected )
!              else
!                ! iflag from (CUBIT) mesh      
!                iselected = iflag_attenuation_store(i,j,k,ispec)                
!              endif
!
!              ! use unrelaxed parameters if attenuation
!              mul = mul * one_minus_sum_beta(iselected)
!               
!            endif
!
!! isotropic case
!            lambdalplus2mul = kappal + FOUR_THIRDS * mul
!            lambdal = lambdalplus2mul - 2.*mul
!
!            ! compute stress sigma
!            sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
!            sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
!            sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl
!
!            sigma_xy = mul*duxdyl_plus_duydxl
!            sigma_xz = mul*duzdxl_plus_duxdzl
!            sigma_yz = mul*duzdyl_plus_duydzl
!
!            !if (SIMULATION_TYPE == 3) then
!            ! b_sigma_xx = lambdalplus2mul*b_duxdxl + lambdal*b_duydyl_plus_duzdzl
!            ! b_sigma_yy = lambdalplus2mul*b_duydyl + lambdal*b_duxdxl_plus_duzdzl
!            ! b_sigma_zz = lambdalplus2mul*b_duzdzl + lambdal*b_duxdxl_plus_duydyl
!            !
!            ! b_sigma_xy = mul*b_duxdyl_plus_duydxl
!            ! b_sigma_xz = mul*b_duzdxl_plus_duxdzl
!            ! b_sigma_yz = mul*b_duzdyl_plus_duydzl
!            !endif
!            
!            ! subtract memory variables if attenuation
!            if(ATTENUATION) then
!              do i_sls = 1,N_SLS
!                R_xx_val = R_xx(i,j,k,ispec,i_sls)
!                R_yy_val = R_yy(i,j,k,ispec,i_sls)
!                sigma_xx = sigma_xx - R_xx_val
!                sigma_yy = sigma_yy - R_yy_val
!                sigma_zz = sigma_zz + R_xx_val + R_yy_val
!                sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
!                sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
!                sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
!              enddo
!             
!              !if (SIMULATION_TYPE == 3) then
!              !  b_R_xx_val = b_R_xx(i,j,k,ispec,i_sls)
!              !  b_R_yy_val = b_R_yy(i,j,k,ispec,i_sls)
!              !  b_sigma_xx = b_sigma_xx - b_R_xx_val
!              !  b_sigma_yy = b_sigma_yy - b_R_yy_val
!              !  b_sigma_zz = b_sigma_zz + b_R_xx_val + b_R_yy_val
!              !  b_sigma_xy = b_sigma_xy - b_R_xy(i,j,k,ispec,i_sls)
!              !  b_sigma_xz = b_sigma_xz - b_R_xz(i,j,k,ispec,i_sls)
!              !  b_sigma_yz = b_sigma_yz - b_R_yz(i,j,k,ispec,i_sls)
!              !endif
!            endif
!      
!  ! form dot product with test vector, symmetric form
!            tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl)
!            tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl)
!            tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)
!
!            tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl)
!            tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl)
!            tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)
!
!            tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl)
!            tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl)
!            tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)
!
!          enddo
!        enddo
!      enddo
!
!  ! subroutines adapted from Deville, Fischer and Mund, High-order methods
!  ! for incompressible fluid flow, Cambridge University Press (2002),
!  ! pages 386 and 389 and Figure 8.3.1
!  ! call mxm_m1_m2_5points(hprimewgll_xxT,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1)
!      do j=1,m2
!        do i=1,m1
!          E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)
!
!          E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)
!
!          E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
!        enddo
!      enddo
!
!  !   call mxm_m1_m1_5points(tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k), &
!  !         hprimewgll_xx,newtempx2(1,1,k),newtempy2(1,1,k),newtempz2(1,1,k))
!      do i=1,m1
!        do j=1,m1
!  ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
!          do k = 1,NGLLX
!            newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempx2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempx2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempx2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempx2(i,5,k)*hprimewgll_xx(5,j)
!
!            newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempy2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempy2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempy2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempy2(i,5,k)*hprimewgll_xx(5,j)
!
!            newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempz2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempz2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempz2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempz2(i,5,k)*hprimewgll_xx(5,j)
!          enddo
!        enddo
!      enddo
!
!  ! call mxm_m2_m1_5points(tempx3,tempy3,tempz3,hprimewgll_xx,newtempx3,newtempy3,newtempz3)
!      do j=1,m1
!        do i=1,m2
!          E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!
!          E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!
!          E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!        enddo
!      enddo
!
!      do k=1,NGLLZ
!        do j=1,NGLLY
!          do i=1,NGLLX
!
!            fac1 = wgllwgll_yz(j,k)
!            fac2 = wgllwgll_xz(i,k)
!            fac3 = wgllwgll_xy(i,j)
!
!  ! sum contributions from each element to the global mesh using indirect addressing
!            iglob = ibool(i,j,k,ispec)
!            accel(1,iglob) = accel(1,iglob) - fac1*newtempx1(i,j,k) - &
!                              fac2*newtempx2(i,j,k) - fac3*newtempx3(i,j,k)
!            accel(2,iglob) = accel(2,iglob) - fac1*newtempy1(i,j,k) - &
!                              fac2*newtempy2(i,j,k) - fac3*newtempy3(i,j,k)
!            accel(3,iglob) = accel(3,iglob) - fac1*newtempz1(i,j,k) - &
!                              fac2*newtempz2(i,j,k) - fac3*newtempz3(i,j,k)
!
!            !  update memory variables based upon the Runge-Kutta scheme
!            if(ATTENUATION) then
!               
!               ! use Runge-Kutta scheme to march in time
!               do i_sls = 1,N_SLS
!
!                  ! get coefficients for that standard linear solid
!                  if( USE_OLSEN_ATTENUATION ) then
!                    vs_val = mustore(i,j,k,ispec) / rho_vs(i,j,k,ispec)
!                    call get_attenuation_model_olsen( vs_val, iselected )
!                  else
!                    iselected = iflag_attenuation_store(i,j,k,ispec)
!                  endif
!                  
!                  factor_loc = mustore(i,j,k,ispec) * factor_common(iselected,i_sls)
!                  alphaval_loc = alphaval(iselected,i_sls)
!                  betaval_loc = betaval(iselected,i_sls)
!                  gammaval_loc = gammaval(iselected,i_sls)
!                  
!                  ! term in xx
!                  Sn   = factor_loc * epsilondev_xx(i,j,k,ispec)
!                  Snp1   = factor_loc * epsilondev_xx_loc(i,j,k)
!                  R_xx(i,j,k,ispec,i_sls) = alphaval_loc * R_xx(i,j,k,ispec,i_sls) + &
!                                    betaval_loc * Sn + gammaval_loc * Snp1
!    
!                  ! term in yy
!                  Sn   = factor_loc * epsilondev_yy(i,j,k,ispec)
!                  Snp1   = factor_loc * epsilondev_yy_loc(i,j,k)
!                  R_yy(i,j,k,ispec,i_sls) = alphaval_loc * R_yy(i,j,k,ispec,i_sls) + &
!                                    betaval_loc * Sn + gammaval_loc * Snp1
!
!                  ! term in zz not computed since zero trace
!                  
!                  ! term in xy
!                  Sn   = factor_loc * epsilondev_xy(i,j,k,ispec)
!                  Snp1   = factor_loc * epsilondev_xy_loc(i,j,k)
!                  R_xy(i,j,k,ispec,i_sls) = alphaval_loc * R_xy(i,j,k,ispec,i_sls) + &
!                                    betaval_loc * Sn + gammaval_loc * Snp1
!                
!                  ! term in xz
!                  Sn   = factor_loc * epsilondev_xz(i,j,k,ispec)
!                  Snp1   = factor_loc * epsilondev_xz_loc(i,j,k)
!                  R_xz(i,j,k,ispec,i_sls) = alphaval_loc * R_xz(i,j,k,ispec,i_sls) + &
!                                    betaval_loc * Sn + gammaval_loc * Snp1
!
!                  ! term in yz
!                  Sn   = factor_loc * epsilondev_yz(i,j,k,ispec)
!                  Snp1   = factor_loc * epsilondev_yz_loc(i,j,k)
!                  R_yz(i,j,k,ispec,i_sls) = alphaval_loc * R_yz(i,j,k,ispec,i_sls) + &
!                                    betaval_loc * Sn + gammaval_loc * Snp1
!                  
!                  !if (SIMULATION_TYPE == 3) then
!                  !  b_alphaval_loc = b_alphaval(iselected,i_sls)
!                  !  b_betaval_loc = b_betaval(iselected,i_sls)
!                  !  b_gammaval_loc = b_gammaval(iselected,i_sls)
!                  !  ! term in xx
!                  !  b_Sn   = factor_loc * b_epsilondev_xx(i,j,k,ispec)
!                  !  b_Snp1   = factor_loc * b_epsilondev_xx_loc(i,j,k)
!                  !  b_R_xx(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xx(i,j,k,ispec,i_sls) + &
!                  !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
!                  !  ! term in yy
!                  !  b_Sn   = factor_loc * b_epsilondev_yy(i,j,k,ispec)
!                  !  b_Snp1   = factor_loc * b_epsilondev_yy_loc(i,j,k)
!                  !  b_R_yy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yy(i,j,k,ispec,i_sls) + &
!                  !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
!                  !  ! term in zz not computed since zero trace
!                  !  ! term in xy
!                  !  b_Sn   = factor_loc * b_epsilondev_xy(i,j,k,ispec)
!                  !  b_Snp1   = factor_loc * b_epsilondev_xy_loc(i,j,k)
!                  !  b_R_xy(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xy(i,j,k,ispec,i_sls) + &
!                  !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
!                  !  ! term in xz
!                  !  b_Sn   = factor_loc * b_epsilondev_xz(i,j,k,ispec)
!                  !  b_Snp1   = factor_loc * b_epsilondev_xz_loc(i,j,k)
!                  !  b_R_xz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_xz(i,j,k,ispec,i_sls) + &
!                  !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
!                  !  ! term in yz
!                  !  b_Sn   = factor_loc * b_epsilondev_yz(i,j,k,ispec)
!                  !  b_Snp1   = factor_loc * b_epsilondev_yz_loc(i,j,k)
!                  !  b_R_yz(i,j,k,ispec,i_sls) = b_alphaval_loc * b_R_yz(i,j,k,ispec,i_sls) + &
!                  !                        b_betaval_loc * b_Sn + b_gammaval_loc * b_Snp1
!                  !endif
!
!               enddo   ! end of loop on memory variables
!
!            endif  !  end attenuation
!
!          enddo
!        enddo
!      enddo
!
!      ! save deviatoric strain for Runge-Kutta scheme
!      if(ATTENUATION) then
!        epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
!        epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
!        epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
!        epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
!        epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
!        !if (SIMULATION_TYPE == 3) then
!        !  b_epsilondev_xx(:,:,:,ispec) = b_epsilondev_xx_loc(:,:,:)
!        !  b_epsilondev_yy(:,:,:,ispec) = b_epsilondev_yy_loc(:,:,:)
!        !  b_epsilondev_xy(:,:,:,ispec) = b_epsilondev_xy_loc(:,:,:)
!        !  b_epsilondev_xz(:,:,:,ispec) = b_epsilondev_xz_loc(:,:,:)
!        !  b_epsilondev_yz(:,:,:,ispec) = b_epsilondev_yz_loc(:,:,:)
!        !endif         
!      endif
!
!    endif ! if (ispec_is_inner(ispec) .eqv. phase_is_inner)
!
!  enddo  ! spectral element loop
!
!end subroutine compute_forces_with_Deville_noanisotropy




!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! subroutines adapted from Deville, Fischer and Mund, High-order methods
!! for incompressible fluid flow, Cambridge University Press (2002),
!! pages 386 and 389 and Figure 8.3.1
!
!  subroutine old_mxm_m1_m2_5points(A,B1,B2,B3,C1,C2,C3)
!
!  implicit none
!
!  include "constants.h"
!
!  real(kind=4), dimension(m1,NGLLX) :: A
!  real(kind=4), dimension(NGLLX,m2) :: B1,B2,B3
!  real(kind=4), dimension(m1,m2) :: C1,C2,C3
!
!  integer :: i,j
!
!  do j=1,m2
!    do i=1,m1
!
!      C1(i,j) = A(i,1)*B1(1,j) + &
!                A(i,2)*B1(2,j) + &
!                A(i,3)*B1(3,j) + &
!                A(i,4)*B1(4,j) + &
!                A(i,5)*B1(5,j)
!
!      C2(i,j) = A(i,1)*B2(1,j) + &
!                A(i,2)*B2(2,j) + &
!                A(i,3)*B2(3,j) + &
!                A(i,4)*B2(4,j) + &
!                A(i,5)*B2(5,j)
!
!      C3(i,j) = A(i,1)*B3(1,j) + &
!                A(i,2)*B3(2,j) + &
!                A(i,3)*B3(3,j) + &
!                A(i,4)*B3(4,j) + &
!                A(i,5)*B3(5,j)
!
!    enddo
!  enddo
!
!  end subroutine old_mxm_m1_m2_5points
!
!!---------
!
!  subroutine old_mxm_m1_m1_5points(A1,A2,A3,B,C1,C2,C3)
!
!  implicit none
!
!  include "constants.h"
!
!  real(kind=4), dimension(m1,NGLLX) :: A1,A2,A3
!  real(kind=4), dimension(NGLLX,m1) :: B
!  real(kind=4), dimension(m1,m1) :: C1,C2,C3
!
!  integer :: i,j
!
!  do j=1,m1
!    do i=1,m1
!
!      C1(i,j) = A1(i,1)*B(1,j) + &
!                A1(i,2)*B(2,j) + &
!                A1(i,3)*B(3,j) + &
!                A1(i,4)*B(4,j) + &
!                A1(i,5)*B(5,j)
!
!      C2(i,j) = A2(i,1)*B(1,j) + &
!                A2(i,2)*B(2,j) + &
!                A2(i,3)*B(3,j) + &
!                A2(i,4)*B(4,j) + &
!                A2(i,5)*B(5,j)
!
!      C3(i,j) = A3(i,1)*B(1,j) + &
!                A3(i,2)*B(2,j) + &
!                A3(i,3)*B(3,j) + &
!                A3(i,4)*B(4,j) + &
!                A3(i,5)*B(5,j)
!
!    enddo
!  enddo
!
!  end subroutine old_mxm_m1_m1_5points
!
!!---------
!
!  subroutine old_mxm_m2_m1_5points(A1,A2,A3,B,C1,C2,C3)
!
!  implicit none
!
!  include "constants.h"
!
!  real(kind=4), dimension(m2,NGLLX) :: A1,A2,A3
!  real(kind=4), dimension(NGLLX,m1) :: B
!  real(kind=4), dimension(m2,m1) :: C1,C2,C3
!
!  integer :: i,j
!
!  do j=1,m1
!    do i=1,m2
!
!      C1(i,j) = A1(i,1)*B(1,j) + &
!                A1(i,2)*B(2,j) + &
!                A1(i,3)*B(3,j) + &
!                A1(i,4)*B(4,j) + &
!                A1(i,5)*B(5,j)
!
!      C2(i,j) = A2(i,1)*B(1,j) + &
!                A2(i,2)*B(2,j) + &
!                A2(i,3)*B(3,j) + &
!                A2(i,4)*B(4,j) + &
!                A2(i,5)*B(5,j)
!
!      C3(i,j) = A3(i,1)*B(1,j) + &
!                A3(i,2)*B(2,j) + &
!                A3(i,3)*B(3,j) + &
!                A3(i,4)*B(4,j) + &
!                A3(i,5)*B(5,j)
!
!    enddo
!  enddo
!
!  end subroutine old_mxm_m2_m1_5points
!


!subroutine compute_forces_with_Deville(phase_is_inner,NSPEC_AB,NGLOB_AB,&
!                      ATTENUATION,USE_OLSEN_ATTENUATION,displ,accel, &
!                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!                      hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
!                      kappastore,mustore,jacobian,ibool,ispec_is_inner, &
!                      NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
!                      xi_source,eta_source,gamma_source,nu_source, &
!                      hdur,hdur_gaussian,t_cmt,dt,stf,t0,sourcearrays, &
!                      one_minus_sum_beta,factor_common,alphaval,betaval,gammaval, &
!                      NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
!                      epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,iflag_attenuation_store, &
!                      ABSORBING_CONDITIONS, &
!                      abs_boundary_normal,abs_boundary_jacobian2Dw, &
!                      abs_boundary_ijk,abs_boundary_ispec, &
!                      num_abs_boundary_faces, &
!                      veloc,rho_vp,rho_vs)
!
!  implicit none
!
!  include "constants.h"
!!  include values created by the mesher
!!  include "OUTPUT_FILES/values_from_mesher.h"
!
!  integer :: NSPEC_AB,NGLOB_AB
!
!! displacement and acceleration
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,accel
!
!! arrays with mesh parameters per slice
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
!        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
!        kappastore,mustore,jacobian
!
!! array with derivatives of Lagrange polynomials and precalculated products
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz
!
!! communication overlap
!  logical, dimension(NSPEC_AB) :: ispec_is_inner
!  logical :: phase_is_inner
!
!! source
!  integer :: NSOURCES,myrank,it
!  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
!  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
!  double precision, dimension(3,3,NSOURCES) :: nu_source
!  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,t_cmt 
!  double precision :: dt
!  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays 
!
!!  integer :: isource
!  double precision :: t0 
!  double precision :: stf 
!
!! memory variables and standard linear solids for attenuation  
!  logical :: ATTENUATION,USE_OLSEN_ATTENUATION
!  integer :: NSPEC_ATTENUATION_AB
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: iflag_attenuation_store
!  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: one_minus_sum_beta
!  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: factor_common, alphaval,betaval,gammaval
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
!       R_xx,R_yy,R_xy,R_xz,R_yz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: &
!       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
!  
!! Stacey conditions
!  logical  :: ABSORBING_CONDITIONS
!!  integer  :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,nspec2D_top
!!  integer  :: NSPEC2DMAX_XMIN_XMAX_ext,NSPEC2DMAX_YMIN_YMAX_ext
!!  integer, dimension(nspec2D_xmin) :: ibelm_xmin
!!  integer, dimension(nspec2D_xmax) :: ibelm_xmax
!!  integer, dimension(nspec2D_ymin) :: ibelm_ymin
!!  integer, dimension(nspec2D_ymax) :: ibelm_ymax
!!  integer, dimension(nspec2D_bottom) :: ibelm_bottom
!!  integer, dimension(nspec2D_top) :: ibelm_top
!!  integer :: ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin),ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax), &
!!            ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin),ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax), &
!!            ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom),ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top)  
!!  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_ext) :: nimin,nimax,nkmin_eta
!!  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_ext) :: njmin,njmax,nkmin_xi
!
!
!
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs
!  
!!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmin) :: jacobian2D_xmin
!!  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nspec2D_xmax) :: jacobian2D_xmax
!!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymin) :: jacobian2D_ymin
!!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec2D_ymax) :: jacobian2D_ymax
!!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM) :: jacobian2D_bottom
!!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_top) :: jacobian2D_top
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmin) :: normal_xmin
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_xmax) :: normal_xmax
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_ymin) :: normal_ymin
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nspec2D_ymax) :: normal_ymax
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM) :: normal_bottom
!!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_top) :: normal_top
!
!  integer :: num_abs_boundary_faces
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_abs_boundary_faces) :: abs_boundary_normal
!  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_abs_boundary_faces) :: abs_boundary_jacobian2Dw
!  integer, dimension(3,NGLLSQUARE,num_abs_boundary_faces) :: abs_boundary_ijk
!  integer, dimension(num_abs_boundary_faces) :: abs_boundary_ispec
!
!
!! computes elastic stiffness term
!  call compute_forces_add_elastic_term(NSPEC_AB,NGLOB_AB,&
!                                ATTENUATION,USE_OLSEN_ATTENUATION,displ,accel, &
!                                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!                                hprime_xx,hprime_xxT,&
!                                hprimewgll_xx,hprimewgll_xxT,&
!                                wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
!                                kappastore,mustore,jacobian,ibool,ispec_is_inner,phase_is_inner, &
!                                one_minus_sum_beta,factor_common,alphaval,betaval,gammaval,&
!                                NSPEC_ATTENUATION_AB,R_xx,R_yy,R_xy,R_xz,R_yz, &
!                                epsilondev_xx,epsilondev_yy,epsilondev_xy,&
!                                epsilondev_xz,epsilondev_yz,iflag_attenuation_store,&
!                                rho_vs )
!
!
!end subroutine compute_forces_with_Deville

!
!-------------------------------------------------------------------------------------------------
!

