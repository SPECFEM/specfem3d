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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


  subroutine compute_forces_viscoelastic(iphase,deltat, &
                                         displ,veloc,accel, &
                                         alphaval,betaval,gammaval, &
                                         R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                         R_trace_lddrk, &
                                         R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                         epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                         epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                         backward_simulation)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,ONE_THIRD,FOUR_THIRDS, &
    m1,m2

  use fault_solver_common, only: Kelvin_Voigt_eta,USE_KELVIN_VOIGT_DAMPING

  use specfem_par, only: SAVE_MOHO_MESH,USE_LDDRK, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,jacobianstore, &
                         NSPEC_AB,NGLOB_AB, &
                         hprime_xx,hprime_xxT, &
                         hprime_yy,hprime_yyT, &
                         hprime_zz,hprime_zzT, &
                         hprimewgll_xx,hprimewgll_xxT, &
                         hprimewgll_yy,hprimewgll_zz, &
                         kappastore,mustore,ibool, &
                         ATTENUATION, &
                         NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_LDDRK, &
                         ANISOTROPY,SIMULATION_TYPE, &
                         NSPEC_ADJOINT,is_moho_top,is_moho_bot, &
                         irregular_element_number,xix_regular,jacobian_regular

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D
  !or: use specfem_par, only: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  use specfem_par_elastic, only: c11store,c12store,c13store,c14store,c15store,c16store, &
                                 c22store,c23store,c24store,c25store,c26store,c33store, &
                                 c34store,c35store,c36store,c44store,c45store,c46store, &
                                 c55store,c56store,c66store,factor_common, &
                                 factor_common_kappa,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                                 dsdx_top,dsdx_bot, &
                                 ispec2D_moho_top,ispec2D_moho_bot, &
                                 nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  use pml_par, only: is_CPML,NSPEC_CPML

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  integer,intent(in) :: iphase
  real(kind=CUSTOM_REAL), intent(in) :: deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(inout) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB),intent(inout) :: R_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilondev_trace

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),intent(out) :: epsilon_trace_over_3

  ! lddrk for update the memory variables
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),intent(inout) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK),intent(inout) :: R_trace_lddrk

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

! local parameters

  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...

  ! arrays for elemental computations inside a given spectral element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempy1,newtempy2,newtempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempz1,newtempz2,newtempz3

  ! attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_att,tempx2_att,tempx3_att, &
            tempy1_att,tempy2_att,tempy3_att, &
            tempz1_att,tempz2_att,tempz3_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att,dummyy_loc_att,dummyz_loc_att
  real(kind=CUSTOM_REAL) :: duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att
  real(kind=CUSTOM_REAL) :: duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

  ! faults
  real(kind=CUSTOM_REAL) :: eta

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) :: R_trace_kappa_sum,R_xx_sum,R_yy_sum
  real(kind=CUSTOM_REAL) :: templ

  integer :: ispec,iglob,ispec_p,num_elements
  integer :: ispec_irreg,ispec2D
  integer :: i,j,k,l
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

! openmp solver
!$OMP PARALLEL if (num_elements > 100) &
!$OMP DEFAULT(NONE) &
!$OMP SHARED( &
!$OMP num_elements,ibool, &
!$OMP iphase,phase_ispec_inner_elastic, &
!$OMP irregular_element_number,jacobian_regular,xix_regular, &
!$OMP displ,veloc,accel, &
!$OMP is_CPML,backward_simulation, &
!$OMP xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore, &
!$OMP kappastore,mustore, &
!$OMP Kelvin_Voigt_eta,USE_KELVIN_VOIGT_DAMPING, &
!$OMP deltat, &
!$OMP c11store,c12store,c13store,c14store,c15store,c16store, &
!$OMP c22store,c23store,c24store,c25store,c26store,c33store, &
!$OMP c34store,c35store,c36store,c44store,c45store,c46store, &
!$OMP c55store,c56store,c66store, &
!$OMP factor_common,factor_common_kappa, &
!$OMP COMPUTE_AND_STORE_STRAIN,ATTENUATION,ANISOTROPY,SIMULATION_TYPE, &
!$OMP R_xx,R_yy,R_xy,R_xz,R_yz,R_trace, &
!$OMP epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,epsilondev_trace,epsilon_trace_over_3, &
!$OMP USE_LDDRK,R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk,R_trace_lddrk, &
!$OMP NSPEC_AB,NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_LDDRK,NSPEC_STRAIN_ONLY, &
!$OMP SAVE_MOHO_MESH,dsdx_top,dsdx_bot,ispec2D_moho_top,ispec2D_moho_bot,is_moho_top,is_moho_bot &
!$OMP ) &
!$OMP PRIVATE( &
!$OMP ispec_p,ispec,ispec_irreg,i,j,k,l,iglob,ispec2D, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#endif
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl,eta, &
!$OMP duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
!$OMP duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
!$OMP sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
!$OMP c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
!$OMP lambdal,mul,lambdalplus2mul,kappal, &
!$OMP hp1,hp2,hp3,fac1,fac2,fac3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc, &
!$OMP tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l, &
!$OMP dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
!$OMP duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att, &
!$OMP duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att, &
!$OMP duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att, &
!$OMP tempx1_att,tempx2_att,tempx3_att,tempy1_att,tempy2_att,tempy3_att,tempz1_att,tempz2_att,tempz3_att, &
!$OMP epsilondev_trace_loc, epsilondev_xx_loc,epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc, &
!$OMP R_trace_kappa_sum,R_xx_sum,R_yy_sum,templ &
!$OMP ) &
!$OMP FIRSTPRIVATE( &
!$OMP hprime_xx,hprime_xxT,hprimewgll_xxT,hprimewgll_xx, &
!$OMP hprime_yy,hprime_yyT,hprimewgll_yy, &
!$OMP hprime_zz,hprime_zzT,hprimewgll_zz, &
!$OMP wgllwgll_yz_3D,wgllwgll_xz_3D,wgllwgll_xy_3D, &
!$OMP alphaval,betaval,gammaval &
!$OMP )

  ! loop over spectral elements
!$OMP DO
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! selects element contribution
    ! PML elements will be computed later
    if (is_CPML(ispec) .and. .not. backward_simulation) cycle

    ! no PML elements from here on

    ! stores displacment values in local array
    if (USE_KELVIN_VOIGT_DAMPING) then
      ! Kelvin Voigt damping: artificial viscosity around dynamic faults
      eta = Kelvin_Voigt_eta(ispec)
      if (is_CPML(ispec) .and. eta /= 0._CUSTOM_REAL) stop 'you cannot put a fault inside a PML layer'
      ! note: this loop will not fully vectorize because it contains a dependency
      !       (through indirect addressing with array ibool())
      !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
      !       which helps the compiler to unroll the innermost loop
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob) + eta * veloc(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob) + eta * veloc(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob) + eta * veloc(3,iglob)
          enddo
        enddo
      enddo
    else
      ! displacement only (without damping)
      ! note: this loop will not fully vectorize because it contains a dependency
      !       (through indirect addressing with array ibool())
      !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
      !       which helps the compiler to unroll the innermost loop
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
    endif

    !------------------------------------------------------------------------------
    !---------------------computation of strain in element-------------------------
    !------------------------------------------------------------------------------

    ! derivative along x, y, z
    ! first double loop over GLL points to compute and store gradients

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for temp1
    ! computes 2. matrix multiplication for temp2
    ! computes 3. matrix multiplication for temp3
    select case (NGLLX)
    case (5)
      call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm5_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
    case (6)
      call mxm6_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm6_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm6_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
    case (7)
      call mxm7_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm7_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm7_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
    case (8)
      call mxm8_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm8_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm8_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
    case default
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL

            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL

            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL
            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l = 1,NGLLX
              hp1 = hprime_xx(i,l)
              tempx1l = tempx1l + dummyx_loc(l,j,k) * hp1
              tempy1l = tempy1l + dummyy_loc(l,j,k) * hp1
              tempz1l = tempz1l + dummyz_loc(l,j,k) * hp1

              hp2 = hprime_yy(j,l)
              tempx2l = tempx2l + dummyx_loc(i,l,k) * hp2
              tempy2l = tempy2l + dummyy_loc(i,l,k) * hp2
              tempz2l = tempz2l + dummyz_loc(i,l,k) * hp2

              hp3 = hprime_zz(k,l)
              tempx3l = tempx3l + dummyx_loc(i,j,l) * hp3
              tempy3l = tempy3l + dummyy_loc(i,j,l) * hp3
              tempz3l = tempz3l + dummyz_loc(i,j,l) * hp3
            enddo
            tempx1(i,j,k) = tempx1l
            tempx2(i,j,k) = tempx2l
            tempx3(i,j,k) = tempx3l

            tempy1(i,j,k) = tempy1l
            tempy2(i,j,k) = tempy2l
            tempy3(i,j,k) = tempy3l

            tempz1(i,j,k) = tempz1l
            tempz2(i,j,k) = tempz2l
            tempz3(i,j,k) = tempz3l
          enddo
        enddo
      enddo
    end select

    if (COMPUTE_AND_STORE_STRAIN) then
      if (ATTENUATION .and. .not. is_CPML(ispec)) then
        ! use first order Taylor expansion of displacement for local storage of stresses
        ! at this current time step, to fix attenuation in a consistent way
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              iglob = ibool(i,j,k,ispec)
              dummyx_loc_att(i,j,k) = deltat * veloc(1,iglob)
              dummyy_loc_att(i,j,k) = deltat * veloc(2,iglob)
              dummyz_loc_att(i,j,k) = deltat * veloc(3,iglob)
            enddo
          enddo
        enddo

        call compute_strain_in_element_att(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                                           tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                                           tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                                           dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                                           hprime_xxT,hprime_yyT,hprime_zzT)
      endif
    endif

    ! grad(u)
    ispec_irreg = irregular_element_number(ispec)
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

        duxdxl(INDEX_IJK) = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
        duxdyl(INDEX_IJK) = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
        duxdzl(INDEX_IJK) = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

        duydxl(INDEX_IJK) = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
        duydyl(INDEX_IJK) = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
        duydzl(INDEX_IJK) = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

        duzdxl(INDEX_IJK) = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
        duzdyl(INDEX_IJK) = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
        duzdzl(INDEX_IJK) = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)
      ENDDO_LOOP_IJK
    else
      ! regular element
      DO_LOOP_IJK
        duxdxl(INDEX_IJK) = xix_regular*tempx1(INDEX_IJK)
        duxdyl(INDEX_IJK) = xix_regular*tempx2(INDEX_IJK)
        duxdzl(INDEX_IJK) = xix_regular*tempx3(INDEX_IJK)

        duydxl(INDEX_IJK) = xix_regular*tempy1(INDEX_IJK)
        duydyl(INDEX_IJK) = xix_regular*tempy2(INDEX_IJK)
        duydzl(INDEX_IJK) = xix_regular*tempy3(INDEX_IJK)

        duzdxl(INDEX_IJK) = xix_regular*tempz1(INDEX_IJK)
        duzdyl(INDEX_IJK) = xix_regular*tempz2(INDEX_IJK)
        duzdzl(INDEX_IJK) = xix_regular*tempz3(INDEX_IJK)
      ENDDO_LOOP_IJK
    endif

    ! adjoint simulations: moho kernel
    if (SAVE_MOHO_MESH) then
      ! saves strains on moho boundary elements
      if (SIMULATION_TYPE == 3) then
        if (is_moho_top(ispec)) then
          ! gets boundary element index
          ispec2D = ispec2D_moho_top(ispec)
          dsdx_top(1,1,INDEX_IJK,ispec2D) = duxdxl(INDEX_IJK)
          dsdx_top(1,2,INDEX_IJK,ispec2D) = duxdyl(INDEX_IJK)
          dsdx_top(1,3,INDEX_IJK,ispec2D) = duxdzl(INDEX_IJK)
          dsdx_top(2,1,INDEX_IJK,ispec2D) = duydxl(INDEX_IJK)
          dsdx_top(2,2,INDEX_IJK,ispec2D) = duydyl(INDEX_IJK)
          dsdx_top(2,3,INDEX_IJK,ispec2D) = duydzl(INDEX_IJK)
          dsdx_top(3,1,INDEX_IJK,ispec2D) = duzdxl(INDEX_IJK)
          dsdx_top(3,2,INDEX_IJK,ispec2D) = duzdyl(INDEX_IJK)
          dsdx_top(3,3,INDEX_IJK,ispec2D) = duzdzl(INDEX_IJK)
        else if (is_moho_bot(ispec)) then
          ! gets boundary element index
          ispec2D = ispec2D_moho_bot(ispec)
          dsdx_bot(1,1,INDEX_IJK,ispec2D) = duxdxl(INDEX_IJK)
          dsdx_bot(1,2,INDEX_IJK,ispec2D) = duxdyl(INDEX_IJK)
          dsdx_bot(1,3,INDEX_IJK,ispec2D) = duxdzl(INDEX_IJK)
          dsdx_bot(2,1,INDEX_IJK,ispec2D) = duydxl(INDEX_IJK)
          dsdx_bot(2,2,INDEX_IJK,ispec2D) = duydyl(INDEX_IJK)
          dsdx_bot(2,3,INDEX_IJK,ispec2D) = duydzl(INDEX_IJK)
          dsdx_bot(3,1,INDEX_IJK,ispec2D) = duzdxl(INDEX_IJK)
          dsdx_bot(3,2,INDEX_IJK,ispec2D) = duzdyl(INDEX_IJK)
          dsdx_bot(3,3,INDEX_IJK,ispec2D) = duzdzl(INDEX_IJK)
        endif
      endif
    endif ! SAVE_MOHO_MESH

    if (COMPUTE_AND_STORE_STRAIN) then
      if (ATTENUATION) then
        ! temporary variables used for fixing attenuation in a consistent way
        if (.not. is_CPML(ispec)) then
          DO_LOOP_IJK
            if (ispec_irreg /= 0) then
              ! irregular element
              xixl = xixstore(INDEX_IJK,ispec_irreg)
              xiyl = xiystore(INDEX_IJK,ispec_irreg)
              xizl = xizstore(INDEX_IJK,ispec_irreg)
              etaxl = etaxstore(INDEX_IJK,ispec_irreg)
              etayl = etaystore(INDEX_IJK,ispec_irreg)
              etazl = etazstore(INDEX_IJK,ispec_irreg)
              gammaxl = gammaxstore(INDEX_IJK,ispec_irreg)
              gammayl = gammaystore(INDEX_IJK,ispec_irreg)
              gammazl = gammazstore(INDEX_IJK,ispec_irreg)

              duxdxl_att = xixl * tempx1_att(INDEX_IJK) + etaxl * tempx2_att(INDEX_IJK) + gammaxl * tempx3_att(INDEX_IJK)
              duxdyl_att = xiyl * tempx1_att(INDEX_IJK) + etayl * tempx2_att(INDEX_IJK) + gammayl * tempx3_att(INDEX_IJK)
              duxdzl_att = xizl * tempx1_att(INDEX_IJK) + etazl * tempx2_att(INDEX_IJK) + gammazl * tempx3_att(INDEX_IJK)

              duydxl_att = xixl * tempy1_att(INDEX_IJK) + etaxl * tempy2_att(INDEX_IJK) + gammaxl * tempy3_att(INDEX_IJK)
              duydyl_att = xiyl * tempy1_att(INDEX_IJK) + etayl * tempy2_att(INDEX_IJK) + gammayl * tempy3_att(INDEX_IJK)
              duydzl_att = xizl * tempy1_att(INDEX_IJK) + etazl * tempy2_att(INDEX_IJK) + gammazl * tempy3_att(INDEX_IJK)

              duzdxl_att = xixl * tempz1_att(INDEX_IJK) + etaxl * tempz2_att(INDEX_IJK) + gammaxl * tempz3_att(INDEX_IJK)
              duzdyl_att = xiyl * tempz1_att(INDEX_IJK) + etayl * tempz2_att(INDEX_IJK) + gammayl * tempz3_att(INDEX_IJK)
              duzdzl_att = xizl * tempz1_att(INDEX_IJK) + etazl * tempz2_att(INDEX_IJK) + gammazl * tempz3_att(INDEX_IJK)
            else
              ! regular element
              duxdxl_att = xix_regular * tempx1_att(INDEX_IJK)
              duxdyl_att = xix_regular * tempx2_att(INDEX_IJK)
              duxdzl_att = xix_regular * tempx3_att(INDEX_IJK)

              duydxl_att = xix_regular * tempy1_att(INDEX_IJK)
              duydyl_att = xix_regular * tempy2_att(INDEX_IJK)
              duydzl_att = xix_regular * tempy3_att(INDEX_IJK)

              duzdxl_att = xix_regular * tempz1_att(INDEX_IJK)
              duzdyl_att = xix_regular * tempz2_att(INDEX_IJK)
              duzdzl_att = xix_regular * tempz3_att(INDEX_IJK)
            endif

            ! precompute some sums to save CPU time
            duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
            duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
            duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

            ! compute deviatoric strain
            templ = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
            if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
            epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
            epsilondev_xx_loc(INDEX_IJK) = duxdxl_att - templ
            epsilondev_yy_loc(INDEX_IJK) = duydyl_att - templ
            epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl_att
            epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl_att
            epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl_att
          ENDDO_LOOP_IJK
        endif
      else
        ! non-attenuation case
        DO_LOOP_IJK
          ! computes deviatoric strain attenuation and/or for kernel calculations
          templ = ONE_THIRD * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
          if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
          epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
          epsilondev_xx_loc(INDEX_IJK) = duxdxl(INDEX_IJK) - templ
          epsilondev_yy_loc(INDEX_IJK) = duydyl(INDEX_IJK) - templ
          epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duxdyl(INDEX_IJK) + duydxl(INDEX_IJK))
          epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdxl(INDEX_IJK) + duxdzl(INDEX_IJK))
          epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdyl(INDEX_IJK) + duydzl(INDEX_IJK))
        ENDDO_LOOP_IJK
      endif
    endif ! COMPUTE_AND_STORE_STRAIN

    ! stresses
    DO_LOOP_IJK
      ! precompute some sums to save CPU time
      duxdyl_plus_duydxl = duxdyl(INDEX_IJK) + duydxl(INDEX_IJK)
      duzdxl_plus_duxdzl = duzdxl(INDEX_IJK) + duxdzl(INDEX_IJK)
      duzdyl_plus_duydzl = duzdyl(INDEX_IJK) + duydzl(INDEX_IJK)

      ! computes either isotropic or anisotropic element stresses
      if (ANISOTROPY) then
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

        sigma_xx = c11 * duxdxl(INDEX_IJK) + c16 * duxdyl_plus_duydxl + c12 * duydyl(INDEX_IJK) + &
                   c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * duzdzl(INDEX_IJK)
        sigma_yy = c12 * duxdxl(INDEX_IJK) + c26 * duxdyl_plus_duydxl + c22 * duydyl(INDEX_IJK) + &
                   c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * duzdzl(INDEX_IJK)
        sigma_zz = c13 * duxdxl(INDEX_IJK) + c36 * duxdyl_plus_duydxl + c23 * duydyl(INDEX_IJK) + &
                   c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * duzdzl(INDEX_IJK)
        sigma_xy = c16 * duxdxl(INDEX_IJK) + c66 * duxdyl_plus_duydxl + c26 * duydyl(INDEX_IJK) + &
                   c56 * duzdxl_plus_duxdzl + c46 * duzdyl_plus_duydzl + c36 * duzdzl(INDEX_IJK)
        sigma_xz = c15 * duxdxl(INDEX_IJK) + c56 * duxdyl_plus_duydxl + c25 * duydyl(INDEX_IJK) + &
                   c55 * duzdxl_plus_duxdzl + c45 * duzdyl_plus_duydzl + c35 * duzdzl(INDEX_IJK)
        sigma_yz = c14 * duxdxl(INDEX_IJK) + c46 * duxdyl_plus_duydxl + c24 * duydyl(INDEX_IJK) + &
                   c45 * duzdxl_plus_duxdzl + c44 * duzdyl_plus_duydzl + c34 * duzdzl(INDEX_IJK)

      else
        ! isotropic case
        kappal = kappastore(INDEX_IJK,ispec)
        mul = mustore(INDEX_IJK,ispec)

        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

        ! compute stress sigma
        sigma_xx = lambdalplus2mul * duxdxl(INDEX_IJK) + lambdal * (duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
        sigma_yy = lambdalplus2mul * duydyl(INDEX_IJK) + lambdal * (duxdxl(INDEX_IJK) + duzdzl(INDEX_IJK))
        sigma_zz = lambdalplus2mul * duzdzl(INDEX_IJK) + lambdal * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK))

        sigma_xy = mul * duxdyl_plus_duydxl
        sigma_xz = mul * duzdxl_plus_duxdzl
        sigma_yz = mul * duzdyl_plus_duydzl
      endif ! ANISOTROPY

      ! subtract memory variables if attenuation
      if (ATTENUATION .and. .not. is_CPML(ispec)) then
        R_xx_sum = sum(R_xx(:,INDEX_IJK,ispec))
        R_yy_sum = sum(R_yy(:,INDEX_IJK,ispec))
        R_trace_kappa_sum = sum(R_trace(:,INDEX_IJK,ispec))

        ! in case no bulk attenuation is desired:
        !R_trace_kappa_sum = 0.0

        sigma_xx = sigma_xx - R_xx_sum - R_trace_kappa_sum
        sigma_yy = sigma_yy - R_yy_sum - R_trace_kappa_sum
        sigma_zz = sigma_zz + R_xx_sum + R_yy_sum - R_trace_kappa_sum
        sigma_xy = sigma_xy - sum(R_xy(:,INDEX_IJK,ispec))
        sigma_xz = sigma_xz - sum(R_xz(:,INDEX_IJK,ispec))
        sigma_yz = sigma_yz - sum(R_yz(:,INDEX_IJK,ispec))
      endif

      if (.not. is_CPML(ispec)) then
        ! define symmetric components of sigma
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        ! dot product with test vector
        if (ispec_irreg /= 0) then
          ! irregular element
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
          tempx1(INDEX_IJK) = jacobianl * (sigma_xx * xixl + sigma_yx * xiyl + sigma_zx * xizl) ! this goes to accel_x
          tempy1(INDEX_IJK) = jacobianl * (sigma_xy * xixl + sigma_yy * xiyl + sigma_zy * xizl) ! this goes to accel_y
          tempz1(INDEX_IJK) = jacobianl * (sigma_xz * xixl + sigma_yz * xiyl + sigma_zz * xizl) ! this goes to accel_z

          tempx2(INDEX_IJK) = jacobianl * (sigma_xx * etaxl + sigma_yx * etayl + sigma_zx * etazl) ! this goes to accel_x
          tempy2(INDEX_IJK) = jacobianl * (sigma_xy * etaxl + sigma_yy * etayl + sigma_zy * etazl) ! this goes to accel_y
          tempz2(INDEX_IJK) = jacobianl * (sigma_xz * etaxl + sigma_yz * etayl + sigma_zz * etazl) ! this goes to accel_z

          tempx3(INDEX_IJK) = jacobianl * (sigma_xx * gammaxl + sigma_yx * gammayl + sigma_zx * gammazl) ! this goes to accel_x
          tempy3(INDEX_IJK) = jacobianl * (sigma_xy * gammaxl + sigma_yy * gammayl + sigma_zy * gammazl) ! this goes to accel_y
          tempz3(INDEX_IJK) = jacobianl * (sigma_xz * gammaxl + sigma_yz * gammayl + sigma_zz * gammazl) ! this goes to accel_z
        else
          !regular element
          jacobianl = jacobian_regular

          ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
          tempx1(INDEX_IJK) = jacobianl * sigma_xx * xix_regular ! this goes to accel_x
          tempy1(INDEX_IJK) = jacobianl * sigma_xy * xix_regular ! this goes to accel_y
          tempz1(INDEX_IJK) = jacobianl * sigma_xz * xix_regular ! this goes to accel_z

          tempx2(INDEX_IJK) = jacobianl * sigma_yx * xix_regular ! this goes to accel_x
          tempy2(INDEX_IJK) = jacobianl * sigma_yy * xix_regular ! this goes to accel_y
          tempz2(INDEX_IJK) = jacobianl * sigma_yz * xix_regular ! this goes to accel_z

          tempx3(INDEX_IJK) = jacobianl * sigma_zx * xix_regular ! this goes to accel_x
          tempy3(INDEX_IJK) = jacobianl * sigma_zy * xix_regular ! this goes to accel_y
          tempz3(INDEX_IJK) = jacobianl * sigma_zz * xix_regular ! this goes to accel_z
        endif
      endif
    ENDDO_LOOP_IJK

    ! second double-loop over GLL to compute all the terms

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    ! computes 2. matrix multiplication for tempx2,..
    ! computes 3. matrix multiplication for newtempx3,..
    select case (NGLLX)
    case (5)
      call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm5_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (6)
      call mxm6_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm6_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm6_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (7)
      call mxm7_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm7_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm7_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (8)
      call mxm8_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm8_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm8_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case default
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            tempx1l = 0._CUSTOM_REAL
            tempy1l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL

            tempx2l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL

            tempx3l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL
            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              fac1 = hprimewgll_xx(l,i)
              tempx1l = tempx1l + tempx1(l,j,k) * fac1
              tempy1l = tempy1l + tempy1(l,j,k) * fac1
              tempz1l = tempz1l + tempz1(l,j,k) * fac1

              fac2 = hprimewgll_yy(l,j)
              tempx2l = tempx2l + tempx2(i,l,k) * fac2
              tempy2l = tempy2l + tempy2(i,l,k) * fac2
              tempz2l = tempz2l + tempz2(i,l,k) * fac2

              fac3 = hprimewgll_zz(l,k)
              tempx3l = tempx3l + tempx3(i,j,l) * fac3
              tempy3l = tempy3l + tempy3(i,j,l) * fac3
              tempz3l = tempz3l + tempz3(i,j,l) * fac3
            enddo
            newtempx1(i,j,k) = tempx1l
            newtempy1(i,j,k) = tempy1l
            newtempz1(i,j,k) = tempz1l

            newtempx2(i,j,k) = tempx2l
            newtempy2(i,j,k) = tempy2l
            newtempz2(i,j,k) = tempz2l

            newtempx3(i,j,k) = tempx3l
            newtempy3(i,j,k) = tempy3l
            newtempz3(i,j,k) = tempz3l
          enddo
        enddo
      enddo
    end select

    ! adds final contributions
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          fac1 = wgllwgll_yz_3D(i,j,k) !or wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz_3D(i,j,k) !or wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy_3D(i,j,k) !or wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh using indirect addressing
!$OMP ATOMIC
          accel(1,iglob) = accel(1,iglob) - (fac1 * newtempx1(i,j,k) + fac2 * newtempx2(i,j,k) + fac3 * newtempx3(i,j,k))
!$OMP ATOMIC
          accel(2,iglob) = accel(2,iglob) - (fac1 * newtempy1(i,j,k) + fac2 * newtempy2(i,j,k) + fac3 * newtempy3(i,j,k))
!$OMP ATOMIC
          accel(3,iglob) = accel(3,iglob) - (fac1 * newtempz1(i,j,k) + fac2 * newtempz2(i,j,k) + fac3 * newtempz3(i,j,k))
        enddo
      enddo
    enddo

    !  update memory variables based upon the Runge-Kutta scheme
    if (ATTENUATION .and. .not. is_CPML(ispec)) then
      ! use Runge-Kutta scheme to march in time
      if (USE_LDDRK) then
        call compute_element_att_memory_lddrk(ispec,deltat,NSPEC_AB,kappastore,mustore, &
               NSPEC_ATTENUATION_AB,factor_common_kappa, &
               R_trace,epsilondev_trace_loc, &
               NSPEC_ATTENUATION_AB_LDDRK,R_trace_lddrk, &
               factor_common,R_xx,R_yy,R_xy,R_xz,R_yz, &
               R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc, &
               epsilondev_xz_loc,epsilondev_yz_loc)
      else
        ! use Runge-Kutta scheme to march in time
        call compute_element_att_memory_second_order_rk(ispec,alphaval,betaval,gammaval, &
               NSPEC_AB,kappastore,mustore,NSPEC_ATTENUATION_AB,factor_common_kappa, &
               R_trace,epsilondev_trace,epsilondev_trace_loc, &
               factor_common,R_xx,R_yy,R_xy,R_xz,R_yz, &
               NSPEC_STRAIN_ONLY,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
               epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
      endif
    endif

    ! save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN) then
      if (ATTENUATION .and. .not. is_CPML(ispec)) epsilondev_trace(:,:,:,ispec) = epsilondev_trace_loc(:,:,:)
      epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
    endif

  enddo  ! spectral element loop
!$OMP ENDDO
!$OMP END PARALLEL

  ! adds contributions for PML elements
  if (NSPEC_CPML > 0 .and. .not. backward_simulation) then
    call compute_forces_viscoelastic_PML(iphase, &
                                         displ,veloc,accel, &
                                         epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                         epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                         backward_simulation)
  endif


  contains

  !
  !---------------
  !

  ! put the code used for computation of strain in element in a subroutine

    subroutine compute_strain_in_element_att(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                                             tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                                             tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                                             dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,hprime_yyT,hprime_zzT)

    use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: &
      tempx1_att,tempx2_att,tempx3_att, &
      tempy1_att,tempy2_att,tempy3_att, &
      tempz1_att,tempz2_att,tempz3_att

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: &
      tempx1,tempx2,tempx3, &
      tempy1,tempy2,tempy3, &
      tempz1,tempz2,tempz3

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xxT
    real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY),intent(in) :: hprime_yyT
    real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ),intent(in) :: hprime_zzT

    ! local variables
    integer :: i,j,k,l
    real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
    real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
    real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
    real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l

    ! use first order Taylor expansion of displacement for local storage of stresses
    ! at this current time step, to fix attenuation in a consistent way
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          tempx1l = tempx1(i,j,k)
          tempy1l = tempy1(i,j,k)
          tempz1l = tempz1(i,j,k)

          tempx2l = tempx2(i,j,k)
          tempy2l = tempy2(i,j,k)
          tempz2l = tempz2(i,j,k)

          tempx3l = tempx3(i,j,k)
          tempy3l = tempy3(i,j,k)
          tempz3l = tempz3(i,j,k)
          ! we can merge these loops because NGLLX = NGLLY = NGLLZ
          do l=1,NGLLX
            hp1 = hprime_xxT(l,i)
            tempx1l = tempx1l + dummyx_loc(l,j,k) * hp1
            tempy1l = tempy1l + dummyy_loc(l,j,k) * hp1
            tempz1l = tempz1l + dummyz_loc(l,j,k) * hp1

            hp2 = hprime_yyT(l,j)
            tempx2l = tempx2l + dummyx_loc(i,l,k) * hp2
            tempy2l = tempy2l + dummyy_loc(i,l,k) * hp2
            tempz2l = tempz2l + dummyz_loc(i,l,k) * hp2

            hp3 = hprime_zzT(l,k)
            tempx3l = tempx3l + dummyx_loc(i,j,l) * hp3
            tempy3l = tempy3l + dummyy_loc(i,j,l) * hp3
            tempz3l = tempz3l + dummyz_loc(i,j,l) * hp3
          enddo
          tempx1_att(i,j,k) = tempx1l
          tempx2_att(i,j,k) = tempx2l
          tempx3_att(i,j,k) = tempx3l

          tempy1_att(i,j,k) = tempy1l
          tempy2_att(i,j,k) = tempy2l
          tempy3_att(i,j,k) = tempy3l

          tempz1_att(i,j,k) = tempz1l
          tempz2_att(i,j,k) = tempz2l
          tempz3_att(i,j,k) = tempz3l
        enddo
      enddo
    enddo

    end subroutine compute_strain_in_element_att


  end subroutine compute_forces_viscoelastic

!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_forces_viscoelastic_PML(iphase, &
                                             displ,veloc,accel, &
                                             epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                             epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                             backward_simulation)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD,m1,m2

  use fault_solver_common, only: Kelvin_Voigt_eta,USE_KELVIN_VOIGT_DAMPING

  use specfem_par, only: xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                         NGLOB_AB, &
                         hprime_xx,hprime_xxT, &
                         hprime_yy,hprime_zz, &
                         hprimewgll_xx,hprimewgll_xxT, &
                         hprimewgll_yy,hprimewgll_zz, &
                         ibool, &
                         SIMULATION_TYPE,NSPEC_ADJOINT, &
                         irregular_element_number,xix_regular

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D
  !or: use specfem_par, only: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  use specfem_par_elastic, only: COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                                 nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  use pml_par, only: is_CPML,spec_to_CPML, &
                     rmemory_dux_dxl_x,rmemory_duy_dyl_x,rmemory_duz_dzl_x, &
                     rmemory_dux_dyl_x,rmemory_dux_dzl_x,rmemory_duz_dxl_x,rmemory_duy_dxl_x, &
                     rmemory_dux_dxl_y,rmemory_duz_dzl_y,rmemory_duy_dyl_y, &
                     rmemory_duy_dxl_y,rmemory_duy_dzl_y,rmemory_duz_dyl_y,rmemory_dux_dyl_y, &
                     rmemory_dux_dxl_z,rmemory_duy_dyl_z,rmemory_duz_dzl_z, &
                     rmemory_duz_dxl_z,rmemory_duz_dyl_z,rmemory_duy_dzl_z,rmemory_dux_dzl_z, &
                     rmemory_displ_elastic,PML_displ_old,PML_displ_new

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: epsilondev_trace
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT),intent(out) :: epsilon_trace_over_3

  integer,intent(in) :: iphase

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

! local parameters

  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...

  ! arrays for elemental computations inside a given spectral element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempy1,newtempy2,newtempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempz1,newtempz2,newtempz3

  ! faults
  real(kind=CUSTOM_REAL) :: eta

  ! local C-PML absorbing boundary conditions parameters
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_new,dummyy_loc_new,dummyz_loc_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_old,dummyy_loc_old,dummyz_loc_old

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_new,tempx2_new,tempx3_new, &
            tempy1_new,tempy2_new,tempy3_new, &
            tempz1_new,tempz2_new,tempz3_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_old,tempx2_old,tempx3_old, &
            tempy1_old,tempy2_old,tempy3_old, &
            tempz1_old,tempz2_old,tempz3_old

  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...

  ! derivatives of ux, uy and uz with respect to x, y and z
  ! in PML_du* computation displ at "n" time step is used
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dux_dxl,PML_dux_dyl,PML_dux_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duy_dxl,PML_duy_dyl,PML_duy_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duz_dxl,PML_duz_dyl,PML_duz_dzl
  ! in PML_du*_old computation we replace displ with displ_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old
  ! in PML_du*_new computation we replace displ with displ_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dux_dxl_new,PML_dux_dyl_new,PML_dux_dzl_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duy_dxl_new,PML_duy_dyl_new,PML_duy_dzl_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_duz_dxl_new,PML_duz_dyl_new,PML_duz_dzl_new

  ! stores C-PML contribution to update acceleration to the global mesh
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: accel_elastic_CPML

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) :: tempx1l_old,tempx2l_old,tempx3l_old
  real(kind=CUSTOM_REAL) :: tempy1l_old,tempy2l_old,tempy3l_old
  real(kind=CUSTOM_REAL) :: tempz1l_old,tempz2l_old,tempz3l_old

  real(kind=CUSTOM_REAL) :: tempx1l_new,tempx2l_new,tempx3l_new
  real(kind=CUSTOM_REAL) :: tempy1l_new,tempy2l_new,tempy3l_new
  real(kind=CUSTOM_REAL) :: tempz1l_new,tempz2l_new,tempz3l_new

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  ! local strain/attenuation parameters
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    epsilondev_xx_loc, epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) :: templ

  integer :: ispec,iglob,ispec_p,ispec_irreg,num_elements
  integer :: i,j,k,l
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

! openmp solver
!$OMP PARALLEL if (num_elements > 100) &
!$OMP DEFAULT(NONE) &
!$OMP SHARED( &
!$OMP num_elements,ibool, &
!$OMP iphase,phase_ispec_inner_elastic, &
!$OMP irregular_element_number,xix_regular, &
!$OMP displ,veloc,accel, &
!$OMP is_CPML,backward_simulation, &
!$OMP xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
!$OMP Kelvin_Voigt_eta,USE_KELVIN_VOIGT_DAMPING, &
!$OMP COMPUTE_AND_STORE_STRAIN,SIMULATION_TYPE, &
!$OMP epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
!$OMP spec_to_CPML,PML_displ_old,PML_displ_new, &
!$OMP rmemory_dux_dxl_x,rmemory_duy_dyl_x,rmemory_duz_dzl_x, &
!$OMP rmemory_dux_dyl_x,rmemory_dux_dzl_x,rmemory_duz_dxl_x,rmemory_duy_dxl_x, &
!$OMP rmemory_dux_dxl_y,rmemory_duz_dzl_y,rmemory_duy_dyl_y, &
!$OMP rmemory_duy_dxl_y,rmemory_duy_dzl_y,rmemory_duz_dyl_y,rmemory_dux_dyl_y, &
!$OMP rmemory_dux_dxl_z,rmemory_duy_dyl_z,rmemory_duz_dzl_z, &
!$OMP rmemory_duz_dxl_z,rmemory_duz_dyl_z,rmemory_duy_dzl_z,rmemory_dux_dzl_z, &
!$OMP rmemory_displ_elastic &
!$OMP ) &
!$OMP PRIVATE( &
!$OMP ispec_p,ispec,ispec_irreg,i,j,k,l,iglob,ispec_CPML, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#endif
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,eta, &
!$OMP duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
!$OMP hp1,hp2,hp3,fac1,fac2,fac3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc, &
!$OMP tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l, &
!$OMP tempx1l_old,tempx2l_old,tempx3l_old,tempy1l_old,tempy2l_old,tempy3l_old,tempz1l_old,tempz2l_old,tempz3l_old, &
!$OMP tempx1l_new,tempx2l_new,tempx3l_new,tempy1l_new,tempy2l_new,tempy3l_new,tempz1l_new,tempz2l_new,tempz3l_new, &
!$OMP dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,dummyx_loc_old,dummyy_loc_old,dummyz_loc_old, &
!$OMP tempx1_new,tempx2_new,tempx3_new,tempy1_new,tempy2_new,tempy3_new,tempz1_new,tempz2_new,tempz3_new, &
!$OMP tempx1_old,tempx2_old,tempx3_old,tempy1_old,tempy2_old,tempy3_old,tempz1_old,tempz2_old,tempz3_old, &
!$OMP PML_dux_dxl,PML_dux_dyl,PML_dux_dzl,PML_duy_dxl,PML_duy_dyl,PML_duy_dzl,PML_duz_dxl,PML_duz_dyl,PML_duz_dzl, &
!$OMP PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old, &
!$OMP PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old, &
!$OMP PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old, &
!$OMP PML_dux_dxl_new,PML_dux_dyl_new,PML_dux_dzl_new, &
!$OMP PML_duy_dxl_new,PML_duy_dyl_new,PML_duy_dzl_new, &
!$OMP PML_duz_dxl_new,PML_duz_dyl_new,PML_duz_dzl_new, &
!$OMP accel_elastic_CPML, &
!$OMP epsilondev_xx_loc,epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc,templ &
!$OMP ) &
!$OMP FIRSTPRIVATE( &
!$OMP hprime_xx,hprime_xxT,hprimewgll_xxT,hprimewgll_xx, &
!$OMP hprime_yy,hprimewgll_yy, &
!$OMP hprime_zz,hprimewgll_zz, &
!$OMP wgllwgll_yz_3D,wgllwgll_xz_3D,wgllwgll_xy_3D &
!$OMP )

  ! loop over spectral elements
!$OMP DO
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! In backward_simulation involved in SIMULATION_TYPE == 3,
    ! we only use the stored value on edge of PML interface.
    ! Thus no computation needs to be done in the PML region in this case.

    ! only PML elements (in forward runs)
    if (.not. (is_CPML(ispec) .and. .not. backward_simulation)) cycle

    ! Keeping in mind, currently we implement a PML derived based on elastic wave equation
    ! for visco-elastic wave simulation. Thus it is not a PML anymore,
    ! because the solution of PML equation derived based on elastic wave equation is not perfectly matched
    ! with solution of visco-elastic wave equation along the PML interface.
    ! you can see PML in this case as a sponge layer.
    !
    ! Based on limited numerical experiments, that PML implementation will work more or less OK anyway if Q > 70.

    ! stores displacement values in local array

    ! checks
    if (USE_KELVIN_VOIGT_DAMPING) then
      ! Kelvin Voigt damping: artificial viscosity around dynamic faults
      eta = Kelvin_Voigt_eta(ispec)
      if (is_CPML(ispec) .and. eta /= 0._CUSTOM_REAL) stop 'Error: you cannot put a fault inside a PML layer'
    endif

    ! displacement only (without damping)
    ! note: this loop will not fully vectorize because it contains a dependency
    !       (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ(1,iglob)
          dummyy_loc(i,j,k) = displ(2,iglob)
          dummyz_loc(i,j,k) = displ(3,iglob)
        enddo
      enddo
    enddo

    ! additional old and new arrays
    ispec_CPML = spec_to_CPML(ispec)

    DO_LOOP_IJK
      dummyx_loc_old(INDEX_IJK) = PML_displ_old(1,INDEX_IJK,ispec_CPML)
      dummyy_loc_old(INDEX_IJK) = PML_displ_old(2,INDEX_IJK,ispec_CPML)
      dummyz_loc_old(INDEX_IJK) = PML_displ_old(3,INDEX_IJK,ispec_CPML)

      dummyx_loc_new(INDEX_IJK) = PML_displ_new(1,INDEX_IJK,ispec_CPML)
      dummyy_loc_new(INDEX_IJK) = PML_displ_new(2,INDEX_IJK,ispec_CPML)
      dummyz_loc_new(INDEX_IJK) = PML_displ_new(3,INDEX_IJK,ispec_CPML)
    ENDDO_LOOP_IJK

    !------------------------------------------------------------------------------
    !---------------------computation of strain in element-------------------------
    !------------------------------------------------------------------------------

    ! derivative along x, y, z
    ! first double loop over GLL points to compute and store gradients

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for temp1
    ! computes 2. matrix multiplication for temp2
    ! computes 3. matrix multiplication for temp3
    select case (NGLLX)
    case (5)
      call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm5_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
      ! old
      call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,tempx1_old,tempy1_old,tempz1_old,m2)
      call mxm5_3comp_3dmat_single(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m1,hprime_xxT,m1, &
                                   tempx2_old,tempy2_old,tempz2_old,m1)
      call mxm5_3comp_singleB(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m2,hprime_xxT,tempx3_old,tempy3_old,tempz3_old,m1)
      ! new
      call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,tempx1_new,tempy1_new,tempz1_new,m2)
      call mxm5_3comp_3dmat_single(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m1,hprime_xxT,m1, &
                                   tempx2_new,tempy2_new,tempz2_new,m1)
      call mxm5_3comp_singleB(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m2,hprime_xxT,tempx3_new,tempy3_new,tempz3_new,m1)
    case (6)
      call mxm6_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm6_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm6_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
      ! old
      call mxm6_3comp_singleA(hprime_xx,m1,dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,tempx1_old,tempy1_old,tempz1_old,m2)
      call mxm6_3comp_3dmat_single(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m1,hprime_xxT,m1, &
                                   tempx2_old,tempy2_old,tempz2_old,m1)
      call mxm6_3comp_singleB(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m2,hprime_xxT,tempx3_old,tempy3_old,tempz3_old,m1)
      ! new
      call mxm6_3comp_singleA(hprime_xx,m1,dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,tempx1_new,tempy1_new,tempz1_new,m2)
      call mxm6_3comp_3dmat_single(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m1,hprime_xxT,m1, &
                                   tempx2_new,tempy2_new,tempz2_new,m1)
      call mxm6_3comp_singleB(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m2,hprime_xxT,tempx3_new,tempy3_new,tempz3_new,m1)

    case (7)
      call mxm7_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm7_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm7_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
      ! old
      call mxm7_3comp_singleA(hprime_xx,m1,dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,tempx1_old,tempy1_old,tempz1_old,m2)
      call mxm7_3comp_3dmat_single(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m1,hprime_xxT,m1, &
                                   tempx2_old,tempy2_old,tempz2_old,m1)
      call mxm7_3comp_singleB(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m2,hprime_xxT,tempx3_old,tempy3_old,tempz3_old,m1)
      ! new
      call mxm7_3comp_singleA(hprime_xx,m1,dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,tempx1_new,tempy1_new,tempz1_new,m2)
      call mxm7_3comp_3dmat_single(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m1,hprime_xxT,m1, &
                                   tempx2_new,tempy2_new,tempz2_new,m1)
      call mxm7_3comp_singleB(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m2,hprime_xxT,tempx3_new,tempy3_new,tempz3_new,m1)

    case (8)
      call mxm8_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
      call mxm8_3comp_3dmat_single(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,m1)
      call mxm8_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
      ! old
      call mxm8_3comp_singleA(hprime_xx,m1,dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,tempx1_old,tempy1_old,tempz1_old,m2)
      call mxm8_3comp_3dmat_single(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m1,hprime_xxT,m1, &
                                   tempx2_old,tempy2_old,tempz2_old,m1)
      call mxm8_3comp_singleB(dummyx_loc_old,dummyy_loc_old,dummyz_loc_old,m2,hprime_xxT,tempx3_old,tempy3_old,tempz3_old,m1)
      ! new
      call mxm8_3comp_singleA(hprime_xx,m1,dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,tempx1_new,tempy1_new,tempz1_new,m2)
      call mxm8_3comp_3dmat_single(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m1,hprime_xxT,m1, &
                                   tempx2_new,tempy2_new,tempz2_new,m1)
      call mxm8_3comp_singleB(dummyx_loc_new,dummyy_loc_new,dummyz_loc_new,m2,hprime_xxT,tempx3_new,tempy3_new,tempz3_new,m1)

    case default
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL

            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL

            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            ! old
            tempx1l_old = 0._CUSTOM_REAL
            tempx2l_old = 0._CUSTOM_REAL
            tempx3l_old = 0._CUSTOM_REAL

            tempy1l_old = 0._CUSTOM_REAL
            tempy2l_old = 0._CUSTOM_REAL
            tempy3l_old = 0._CUSTOM_REAL

            tempz1l_old = 0._CUSTOM_REAL
            tempz2l_old = 0._CUSTOM_REAL
            tempz3l_old = 0._CUSTOM_REAL

            ! new
            tempx1l_new = 0._CUSTOM_REAL
            tempx2l_new = 0._CUSTOM_REAL
            tempx3l_new = 0._CUSTOM_REAL

            tempy1l_new = 0._CUSTOM_REAL
            tempy2l_new = 0._CUSTOM_REAL
            tempy3l_new = 0._CUSTOM_REAL

            tempz1l_new = 0._CUSTOM_REAL
            tempz2l_new = 0._CUSTOM_REAL
            tempz3l_new = 0._CUSTOM_REAL
            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l = 1,NGLLX
              hp1 = hprime_xx(i,l)
              tempx1l = tempx1l + dummyx_loc(l,j,k) * hp1
              tempy1l = tempy1l + dummyy_loc(l,j,k) * hp1
              tempz1l = tempz1l + dummyz_loc(l,j,k) * hp1
              ! old
              tempx1l_old = tempx1l_old + dummyx_loc_old(l,j,k) * hp1
              tempy1l_old = tempy1l_old + dummyy_loc_old(l,j,k) * hp1
              tempz1l_old = tempz1l_old + dummyz_loc_old(l,j,k) * hp1
              ! new
              tempx1l_new = tempx1l_new + dummyx_loc_new(l,j,k) * hp1
              tempy1l_new = tempy1l_new + dummyy_loc_new(l,j,k) * hp1
              tempz1l_new = tempz1l_new + dummyz_loc_new(l,j,k) * hp1

              hp2 = hprime_yy(j,l)
              tempx2l = tempx2l + dummyx_loc(i,l,k) * hp2
              tempy2l = tempy2l + dummyy_loc(i,l,k) * hp2
              tempz2l = tempz2l + dummyz_loc(i,l,k) * hp2
              ! old
              tempx2l_old = tempx2l_old + dummyx_loc_old(i,l,k) * hp2
              tempy2l_old = tempy2l_old + dummyy_loc_old(i,l,k) * hp2
              tempz2l_old = tempz2l_old + dummyz_loc_old(i,l,k) * hp2
              ! new
              tempx2l_new = tempx2l_new + dummyx_loc_new(i,l,k) * hp2
              tempy2l_new = tempy2l_new + dummyy_loc_new(i,l,k) * hp2
              tempz2l_new = tempz2l_new + dummyz_loc_new(i,l,k) * hp2

              hp3 = hprime_zz(k,l)
              tempx3l = tempx3l + dummyx_loc(i,j,l) * hp3
              tempy3l = tempy3l + dummyy_loc(i,j,l) * hp3
              tempz3l = tempz3l + dummyz_loc(i,j,l) * hp3
              ! old
              tempx3l_old = tempx3l_old + dummyx_loc_old(i,j,l) * hp3
              tempy3l_old = tempy3l_old + dummyy_loc_old(i,j,l) * hp3
              tempz3l_old = tempz3l_old + dummyz_loc_old(i,j,l) * hp3
              ! new
              tempx3l_new = tempx3l_new + dummyx_loc_new(i,j,l) * hp3
              tempy3l_new = tempy3l_new + dummyy_loc_new(i,j,l) * hp3
              tempz3l_new = tempz3l_new + dummyz_loc_new(i,j,l) * hp3
            enddo
            tempx1(i,j,k) = tempx1l
            tempx2(i,j,k) = tempx2l
            tempx3(i,j,k) = tempx3l

            tempy1(i,j,k) = tempy1l
            tempy2(i,j,k) = tempy2l
            tempy3(i,j,k) = tempy3l

            tempz1(i,j,k) = tempz1l
            tempz2(i,j,k) = tempz2l
            tempz3(i,j,k) = tempz3l

            ! old
            tempx1_old(i,j,k) = tempx1l_old
            tempx2_old(i,j,k) = tempx2l_old
            tempx3_old(i,j,k) = tempx3l_old

            tempy1_old(i,j,k) = tempy1l_old
            tempy2_old(i,j,k) = tempy2l_old
            tempy3_old(i,j,k) = tempy3l_old

            tempz1_old(i,j,k) = tempz1l_old
            tempz2_old(i,j,k) = tempz2l_old
            tempz3_old(i,j,k) = tempz3l_old

            ! new
            tempx1_new(i,j,k) = tempx1l_new
            tempx2_new(i,j,k) = tempx2l_new
            tempx3_new(i,j,k) = tempx3l_new

            tempy1_new(i,j,k) = tempy1l_new
            tempy2_new(i,j,k) = tempy2l_new
            tempy3_new(i,j,k) = tempy3l_new

            tempz1_new(i,j,k) = tempz1l_new
            tempz2_new(i,j,k) = tempz2l_new
            tempz3_new(i,j,k) = tempz3l_new
          enddo
        enddo
      enddo
    end select

    ! grad(u)
    ispec_irreg = irregular_element_number(ispec)
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

        duxdxl(INDEX_IJK) = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
        duxdyl(INDEX_IJK) = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
        duxdzl(INDEX_IJK) = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

        duydxl(INDEX_IJK) = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
        duydyl(INDEX_IJK) = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
        duydzl(INDEX_IJK) = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

        duzdxl(INDEX_IJK) = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
        duzdyl(INDEX_IJK) = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
        duzdzl(INDEX_IJK) = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)
      ENDDO_LOOP_IJK
    else
      ! regular element
      DO_LOOP_IJK
        duxdxl(INDEX_IJK) = xix_regular*tempx1(INDEX_IJK)
        duxdyl(INDEX_IJK) = xix_regular*tempx2(INDEX_IJK)
        duxdzl(INDEX_IJK) = xix_regular*tempx3(INDEX_IJK)

        duydxl(INDEX_IJK) = xix_regular*tempy1(INDEX_IJK)
        duydyl(INDEX_IJK) = xix_regular*tempy2(INDEX_IJK)
        duydzl(INDEX_IJK) = xix_regular*tempy3(INDEX_IJK)

        duzdxl(INDEX_IJK) = xix_regular*tempz1(INDEX_IJK)
        duzdyl(INDEX_IJK) = xix_regular*tempz2(INDEX_IJK)
        duzdzl(INDEX_IJK) = xix_regular*tempz3(INDEX_IJK)
      ENDDO_LOOP_IJK
    endif

    ! no moho boundary strain will be stored in PML
    ! (hardly valid inside the PML anyway)
    ! adjoint simulations: moho kernel
    !if (SAVE_MOHO_MESH) then
    !endif

    ! computes deviatoric strain for kernel calculations
    ! (maybe not really needed, but will keep for now based on a "pure" acoustic element contribution)
    if (COMPUTE_AND_STORE_STRAIN) then
      ! non-attenuation case
      DO_LOOP_IJK
        templ = ONE_THIRD * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
        if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
        !epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ ! not needed for PML elements
        epsilondev_xx_loc(INDEX_IJK) = duxdxl(INDEX_IJK) - templ
        epsilondev_yy_loc(INDEX_IJK) = duydyl(INDEX_IJK) - templ
        epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duxdyl(INDEX_IJK) + duydxl(INDEX_IJK))
        epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdxl(INDEX_IJK) + duxdzl(INDEX_IJK))
        epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdyl(INDEX_IJK) + duydzl(INDEX_IJK))
      ENDDO_LOOP_IJK
    endif

    ! stores derivatives of ux, uy and uz with respect to x, y and z
    DO_LOOP_IJK
      PML_dux_dxl(INDEX_IJK) = duxdxl(INDEX_IJK)
      PML_dux_dyl(INDEX_IJK) = duxdyl(INDEX_IJK)
      PML_dux_dzl(INDEX_IJK) = duxdzl(INDEX_IJK)

      PML_duy_dxl(INDEX_IJK) = duydxl(INDEX_IJK)
      PML_duy_dyl(INDEX_IJK) = duydyl(INDEX_IJK)
      PML_duy_dzl(INDEX_IJK) = duydzl(INDEX_IJK)

      PML_duz_dxl(INDEX_IJK) = duzdxl(INDEX_IJK)
      PML_duz_dyl(INDEX_IJK) = duzdyl(INDEX_IJK)
      PML_duz_dzl(INDEX_IJK) = duzdzl(INDEX_IJK)
    ENDDO_LOOP_IJK

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

        ! old
        PML_dux_dxl_old(INDEX_IJK) = &
            xixl * tempx1_old(INDEX_IJK) + etaxl * tempx2_old(INDEX_IJK) + gammaxl * tempx3_old(INDEX_IJK)
        PML_dux_dyl_old(INDEX_IJK) = &
            xiyl * tempx1_old(INDEX_IJK) + etayl * tempx2_old(INDEX_IJK) + gammayl * tempx3_old(INDEX_IJK)
        PML_dux_dzl_old(INDEX_IJK) = &
            xizl * tempx1_old(INDEX_IJK) + etazl * tempx2_old(INDEX_IJK) + gammazl * tempx3_old(INDEX_IJK)

        PML_duy_dxl_old(INDEX_IJK) = &
            xixl * tempy1_old(INDEX_IJK) + etaxl * tempy2_old(INDEX_IJK) + gammaxl * tempy3_old(INDEX_IJK)
        PML_duy_dyl_old(INDEX_IJK) = &
            xiyl * tempy1_old(INDEX_IJK) + etayl * tempy2_old(INDEX_IJK) + gammayl * tempy3_old(INDEX_IJK)
        PML_duy_dzl_old(INDEX_IJK) = &
            xizl * tempy1_old(INDEX_IJK) + etazl * tempy2_old(INDEX_IJK) + gammazl * tempy3_old(INDEX_IJK)

        PML_duz_dxl_old(INDEX_IJK) = &
            xixl * tempz1_old(INDEX_IJK) + etaxl * tempz2_old(INDEX_IJK) + gammaxl * tempz3_old(INDEX_IJK)
        PML_duz_dyl_old(INDEX_IJK) = &
            xiyl * tempz1_old(INDEX_IJK) + etayl * tempz2_old(INDEX_IJK) + gammayl * tempz3_old(INDEX_IJK)
        PML_duz_dzl_old(INDEX_IJK) = &
             xizl * tempz1_old(INDEX_IJK) + etazl * tempz2_old(INDEX_IJK) + gammazl * tempz3_old(INDEX_IJK)

        ! new
        PML_dux_dxl_new(INDEX_IJK) = &
            xixl * tempx1_new(INDEX_IJK) + etaxl * tempx2_new(INDEX_IJK) + gammaxl * tempx3_new(INDEX_IJK)
        PML_dux_dyl_new(INDEX_IJK) = &
            xiyl * tempx1_new(INDEX_IJK) + etayl * tempx2_new(INDEX_IJK) + gammayl * tempx3_new(INDEX_IJK)
        PML_dux_dzl_new(INDEX_IJK) = &
            xizl * tempx1_new(INDEX_IJK) + etazl * tempx2_new(INDEX_IJK) + gammazl * tempx3_new(INDEX_IJK)

        PML_duy_dxl_new(INDEX_IJK) = &
            xixl * tempy1_new(INDEX_IJK) + etaxl * tempy2_new(INDEX_IJK) + gammaxl * tempy3_new(INDEX_IJK)
        PML_duy_dyl_new(INDEX_IJK) = &
            xiyl * tempy1_new(INDEX_IJK) + etayl * tempy2_new(INDEX_IJK) + gammayl * tempy3_new(INDEX_IJK)
        PML_duy_dzl_new(INDEX_IJK) = &
            xizl * tempy1_new(INDEX_IJK) + etazl * tempy2_new(INDEX_IJK) + gammazl * tempy3_new(INDEX_IJK)

        PML_duz_dxl_new(INDEX_IJK) = &
            xixl * tempz1_new(INDEX_IJK) + etaxl * tempz2_new(INDEX_IJK) + gammaxl * tempz3_new(INDEX_IJK)
        PML_duz_dyl_new(INDEX_IJK) = &
            xiyl * tempz1_new(INDEX_IJK) + etayl * tempz2_new(INDEX_IJK) + gammayl * tempz3_new(INDEX_IJK)
        PML_duz_dzl_new(INDEX_IJK) = &
            xizl * tempz1_new(INDEX_IJK) + etazl * tempz2_new(INDEX_IJK) + gammazl * tempz3_new(INDEX_IJK)
      ENDDO_LOOP_IJK
    else
      ! regular element
      DO_LOOP_IJK
        ! old
        PML_dux_dxl_old(INDEX_IJK) = xix_regular * tempx1_old(INDEX_IJK)
        PML_dux_dyl_old(INDEX_IJK) = xix_regular * tempx2_old(INDEX_IJK)
        PML_dux_dzl_old(INDEX_IJK) = xix_regular * tempx3_old(INDEX_IJK)

        PML_duy_dxl_old(INDEX_IJK) = xix_regular * tempy1_old(INDEX_IJK)
        PML_duy_dyl_old(INDEX_IJK) = xix_regular * tempy2_old(INDEX_IJK)
        PML_duy_dzl_old(INDEX_IJK) = xix_regular * tempy3_old(INDEX_IJK)

        PML_duz_dxl_old(INDEX_IJK) = xix_regular * tempz1_old(INDEX_IJK)
        PML_duz_dyl_old(INDEX_IJK) = xix_regular * tempz2_old(INDEX_IJK)
        PML_duz_dzl_old(INDEX_IJK) = xix_regular * tempz3_old(INDEX_IJK)

        ! new
        PML_dux_dxl_new(INDEX_IJK) = xix_regular * tempx1_new(INDEX_IJK)
        PML_dux_dyl_new(INDEX_IJK) = xix_regular * tempx2_new(INDEX_IJK)
        PML_dux_dzl_new(INDEX_IJK) = xix_regular * tempx3_new(INDEX_IJK)

        PML_duy_dxl_new(INDEX_IJK) = xix_regular * tempy1_new(INDEX_IJK)
        PML_duy_dyl_new(INDEX_IJK) = xix_regular * tempy2_new(INDEX_IJK)
        PML_duy_dzl_new(INDEX_IJK) = xix_regular * tempy3_new(INDEX_IJK)

        PML_duz_dxl_new(INDEX_IJK) =  xix_regular * tempz1_new(INDEX_IJK)
        PML_duz_dyl_new(INDEX_IJK) =  xix_regular * tempz2_new(INDEX_IJK)
        PML_duz_dzl_new(INDEX_IJK) =  xix_regular * tempz3_new(INDEX_IJK)
      ENDDO_LOOP_IJK
    endif

    ! In backward_simulation involved in SIMULATION_TYPE == 3,
    ! we only use the stored value on edge of PML interface.
    ! Thus no computation needs to be done in the PML region in this case.

    ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
    call pml_compute_memory_variables_elastic(ispec,ispec_CPML, &
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

    ! calculates contribution from each C-PML element to update acceleration
    call pml_compute_accel_contribution_elastic(ispec,ispec_CPML,displ,veloc, &
                                                accel_elastic_CPML,rmemory_displ_elastic)

    ! add "special" contributions from each element to the global values
    ! will be added at the end with together with "normal" contribution
    !do k = 1,NGLLZ
    !  do j = 1,NGLLY
    !    do i = 1,NGLLX
    !      iglob = ibool(i,j,k,ispec)
    !      accel(1,iglob) = accel(1,iglob) - accel_elastic_CPML(1,i,j,k)
    !      accel(2,iglob) = accel(2,iglob) - accel_elastic_CPML(2,i,j,k)
    !      accel(3,iglob) = accel(3,iglob) - accel_elastic_CPML(3,i,j,k)
    !    enddo
    !  enddo
    !enddo

    ! second double-loop over GLL to compute all the terms

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    ! computes 2. matrix multiplication for tempx2,..
    ! computes 3. matrix multiplication for newtempx3,..
    select case (NGLLX)
    case (5)
      call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm5_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (6)
      call mxm6_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm6_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm6_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (7)
      call mxm7_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm7_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm7_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case (8)
      call mxm8_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      call mxm8_3comp_3dmat_single(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,m1)
      call mxm8_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
    case default
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            tempx1l = 0._CUSTOM_REAL
            tempy1l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL

            tempx2l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL

            tempx3l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              fac1 = hprimewgll_xx(l,i)
              tempx1l = tempx1l + tempx1(l,j,k) * fac1
              tempy1l = tempy1l + tempy1(l,j,k) * fac1
              tempz1l = tempz1l + tempz1(l,j,k) * fac1

              fac2 = hprimewgll_yy(l,j)
              tempx2l = tempx2l + tempx2(i,l,k) * fac2
              tempy2l = tempy2l + tempy2(i,l,k) * fac2
              tempz2l = tempz2l + tempz2(i,l,k) * fac2

              fac3 = hprimewgll_zz(l,k)
              tempx3l = tempx3l + tempx3(i,j,l) * fac3
              tempy3l = tempy3l + tempy3(i,j,l) * fac3
              tempz3l = tempz3l + tempz3(i,j,l) * fac3
            enddo

            newtempx1(i,j,k) = tempx1l
            newtempy1(i,j,k) = tempy1l
            newtempz1(i,j,k) = tempz1l

            newtempx2(i,j,k) = tempx2l
            newtempy2(i,j,k) = tempy2l
            newtempz2(i,j,k) = tempz2l

            newtempx3(i,j,k) = tempx3l
            newtempy3(i,j,k) = tempy3l
            newtempz3(i,j,k) = tempz3l
          enddo
        enddo
      enddo
    end select

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          fac1 = wgllwgll_yz_3D(i,j,k) !or wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz_3D(i,j,k) !or wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy_3D(i,j,k) !or wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh using indirect addressing
!$OMP ATOMIC
          accel(1,iglob) = accel(1,iglob) &
            - (fac1 * newtempx1(i,j,k) + fac2 * newtempx2(i,j,k) + fac3 * newtempx3(i,j,k) &
                + accel_elastic_CPML(1,i,j,k))
!$OMP ATOMIC
          accel(2,iglob) = accel(2,iglob) &
            - (fac1 * newtempy1(i,j,k) + fac2 * newtempy2(i,j,k) + fac3 * newtempy3(i,j,k) &
                + accel_elastic_CPML(2,i,j,k))
!$OMP ATOMIC
          accel(3,iglob) = accel(3,iglob) &
            - (fac1 * newtempz1(i,j,k) + fac2 * newtempz2(i,j,k) + fac3 * newtempz3(i,j,k) &
                + accel_elastic_CPML(3,i,j,k))
        enddo
      enddo
    enddo

    ! save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN) then
      ! if (ATTENUATION .and. .not. is_CPML(ispec)) epsilondev_trace(:,:,:,ispec) = epsilondev_trace_loc(:,:,:)
      epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
    endif

  enddo  ! spectral element loop
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_forces_viscoelastic_PML




!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices (5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications is in general slower
!
! please leave the routines here to help compilers inline the code

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_singleA
! cray
!DIR$ INLINEALWAYS mxm5_3comp_singleA

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same A matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA

  !-------------

  subroutine mxm6_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm6_3comp_singleA
! cray
!DIR$ INLINEALWAYS mxm6_3comp_singleA

! 3 different arrays for x/y/z-components, 2-dimensional arrays (36,6)/(6,36), same A matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j) &
               + A(i,6) * B1(6,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j) &
               + A(i,6) * B2(6,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j) &
               + A(i,6) * B3(6,j)
    enddo
  enddo

  end subroutine mxm6_3comp_singleA

  !-------------

  subroutine mxm7_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm7_3comp_singleA
! cray
!DIR$ INLINEALWAYS mxm7_3comp_singleA

! 3 different arrays for x/y/z-components, 2-dimensional arrays (49,7)/(7,49), same A matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j) &
               + A(i,6) * B1(6,j) &
               + A(i,7) * B1(7,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j) &
               + A(i,6) * B2(6,j) &
               + A(i,7) * B2(7,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j) &
               + A(i,6) * B3(6,j) &
               + A(i,7) * B3(7,j)
    enddo
  enddo

  end subroutine mxm7_3comp_singleA

  !-------------

  subroutine mxm8_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm8_3comp_singleA
! cray
!DIR$ INLINEALWAYS mxm8_3comp_singleA

! 3 different arrays for x/y/z-components, 2-dimensional arrays (64,8)/(8,64), same A matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,8),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(8,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j) &
               + A(i,6) * B1(6,j) &
               + A(i,7) * B1(7,j) &
               + A(i,8) * B1(8,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j) &
               + A(i,6) * B2(6,j) &
               + A(i,7) * B2(7,j) &
               + A(i,8) * B2(8,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j) &
               + A(i,6) * B3(6,j) &
               + A(i,7) * B3(7,j) &
               + A(i,8) * B3(8,j)
    enddo
  enddo

  end subroutine mxm8_3comp_singleA


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_singleB
! cray
!DIR$ INLINEALWAYS mxm5_3comp_singleB

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleB

  !-------------

  subroutine mxm6_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm6_3comp_singleB
! cray
!DIR$ INLINEALWAYS mxm6_3comp_singleB

! 3 different arrays for x/y/z-components, 2-dimensional arrays (36,6)/(6,36), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j) &
               + A1(i,6) * B(6,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j) &
               + A2(i,6) * B(6,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j) &
               + A3(i,6) * B(6,j)
    enddo
  enddo

  end subroutine mxm6_3comp_singleB

  !-------------

  subroutine mxm7_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm6_3comp_singleB
! cray
!DIR$ INLINEALWAYS mxm6_3comp_singleB

! 3 different arrays for x/y/z-components, 2-dimensional arrays (49,7)/(7,49), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j) &
               + A1(i,6) * B(6,j) &
               + A1(i,7) * B(7,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j) &
               + A2(i,6) * B(6,j) &
               + A2(i,7) * B(7,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j) &
               + A3(i,6) * B(6,j) &
               + A3(i,7) * B(7,j)
    enddo
  enddo

  end subroutine mxm7_3comp_singleB

  !-------------

  subroutine mxm8_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm6_3comp_singleB
! cray
!DIR$ INLINEALWAYS mxm6_3comp_singleB

! 3 different arrays for x/y/z-components, 2-dimensional arrays (64,8)/(8,64), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,8),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(8,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j) &
               + A1(i,6) * B(6,j) &
               + A1(i,7) * B(7,j) &
               + A1(i,8) * B(8,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j) &
               + A2(i,6) * B(6,j) &
               + A2(i,7) * B(7,j) &
               + A2(i,8) * B(8,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j) &
               + A3(i,6) * B(6,j) &
               + A3(i,7) * B(7,j) &
               + A3(i,8) * B(8,j)
    enddo
  enddo

  end subroutine mxm8_3comp_singleB


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_3dmat_single(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_3dmat_single
! cray
!DIR$ INLINEALWAYS mxm5_3comp_3dmat_single

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3comp_3dmat_single

  !-------------

  subroutine mxm6_3comp_3dmat_single(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm6_3comp_3dmat_single
! cray
!DIR$ INLINEALWAYS mxm6_3comp_3dmat_single

! 3 different arrays for x/y/z-components, 3-dimensional arrays (6,6,6), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,6,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(6,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j) &
                   + A1(i,6,k) * B(6,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j) &
                   + A2(i,6,k) * B(6,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j) &
                   + A3(i,6,k) * B(6,j)
      enddo
    enddo
  enddo

  end subroutine mxm6_3comp_3dmat_single

  !-------------

  subroutine mxm7_3comp_3dmat_single(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm7_3comp_3dmat_single
! cray
!DIR$ INLINEALWAYS mxm7_3comp_3dmat_single

! 3 different arrays for x/y/z-components, 3-dimensional arrays (7,7,7), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,7,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(7,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j) &
                   + A1(i,6,k) * B(6,j) &
                   + A1(i,7,k) * B(7,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j) &
                   + A2(i,6,k) * B(6,j) &
                   + A2(i,7,k) * B(7,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j) &
                   + A3(i,6,k) * B(6,j) &
                   + A3(i,7,k) * B(7,j)
      enddo
    enddo
  enddo

  end subroutine mxm7_3comp_3dmat_single

  !-------------

  subroutine mxm8_3comp_3dmat_single(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
!DIR$ ATTRIBUTES FORCEINLINE :: mxm8_3comp_3dmat_single
! cray
!DIR$ INLINEALWAYS mxm8_3comp_3dmat_single

! 3 different arrays for x/y/z-components, 3-dimensional arrays (8,8,8), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,8,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(8,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j) &
                   + A1(i,6,k) * B(6,j) &
                   + A1(i,7,k) * B(7,j) &
                   + A1(i,8,k) * B(8,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j) &
                   + A2(i,6,k) * B(6,j) &
                   + A2(i,7,k) * B(7,j) &
                   + A2(i,8,k) * B(8,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j) &
                   + A3(i,6,k) * B(6,j) &
                   + A3(i,7,k) * B(7,j) &
                   + A3(i,8,k) * B(8,j)
      enddo
    enddo
  enddo

  end subroutine mxm8_3comp_3dmat_single

