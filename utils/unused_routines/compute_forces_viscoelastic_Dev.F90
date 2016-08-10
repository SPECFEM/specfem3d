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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


! Deville routine for NGLL == 5, 6 or 7 (default in constants.h is 5)

  subroutine compute_forces_viscoelastic_Dev(iphase,NSPEC_AB,NGLOB_AB, &
                                    displ,veloc,accel, &
                                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                    hprime_xx,hprime_xxT, &
                                    hprimewgll_xx,hprimewgll_xxT, &
                                    wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                                    kappastore,mustore,jacobian,ibool, &
                                    ATTENUATION,deltat,PML_CONDITIONS, &
                                    one_minus_sum_beta,factor_common, &
                                    one_minus_sum_beta_kappa,factor_common_kappa, &
                                    alphaval,betaval,gammaval, &
                                    NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_Kappa, &
                                    R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                    epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                    epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                    ANISOTROPY,NSPEC_ANISO, &
                                    c11store,c12store,c13store,c14store,c15store,c16store, &
                                    c22store,c23store,c24store,c25store,c26store,c33store, &
                                    c34store,c35store,c36store,c44store,c45store,c46store, &
                                    c55store,c56store,c66store, &
                                    SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                                    NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT, &
                                    is_moho_top,is_moho_bot, &
                                    dsdx_top,dsdx_bot, &
                                    ispec2D_moho_top,ispec2D_moho_bot, &
                                    num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                                    phase_ispec_inner_elastic,backward_simulation)


! computes elastic tensor term

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,NGLLCUBE, &
                       N_SLS,ONE_THIRD,FOUR_THIRDS,m1,m2, &
                       MAKE_HOOKE_LAW_WEAKLY_NONLINEAR,A,B,C,A_over_4,B_over_2
  use fault_solver_dynamic, only: Kelvin_Voigt_eta

  use specfem_par, only: FULL_ATTENUATION_SOLID,SAVE_MOHO_MESH

  use pml_par, only: is_CPML,spec_to_CPML,accel_elastic_CPML,NSPEC_CPML, &
                     PML_dux_dxl,PML_dux_dyl,PML_dux_dzl,PML_duy_dxl,PML_duy_dyl,PML_duy_dzl, &
                     PML_duz_dxl,PML_duz_dyl,PML_duz_dzl, &
                     PML_dux_dxl_old,PML_dux_dyl_old,PML_dux_dzl_old, &
                     PML_duy_dxl_old,PML_duy_dyl_old,PML_duy_dzl_old, &
                     PML_duz_dxl_old,PML_duz_dyl_old,PML_duz_dzl_old, &
                     PML_dux_dxl_new,PML_dux_dyl_new,PML_dux_dzl_new, &
                     PML_duy_dxl_new,PML_duy_dyl_new,PML_duy_dzl_new, &
                     PML_duz_dxl_new,PML_duz_dyl_new,PML_duz_dzl_new, &
                     rmemory_dux_dxl_x,rmemory_duy_dyl_x,rmemory_duz_dzl_x, &
                     rmemory_dux_dyl_x,rmemory_dux_dzl_x,rmemory_duz_dxl_x,rmemory_duy_dxl_x, &
                     rmemory_dux_dxl_y,rmemory_duz_dzl_y,rmemory_duy_dyl_y, &
                     rmemory_duy_dxl_y,rmemory_duy_dzl_y,rmemory_duz_dyl_y,rmemory_dux_dyl_y, &
                     rmemory_dux_dxl_z,rmemory_duy_dyl_z,rmemory_duz_dzl_z, &
                     rmemory_duz_dxl_z,rmemory_duz_dyl_z,rmemory_duy_dzl_z,rmemory_dux_dzl_z, &
                     rmemory_displ_elastic,displ_old,displ_new

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ,veloc,accel

! time step
  real(kind=CUSTOM_REAL) :: deltat

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
            kappastore,mustore,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NGLLZ) :: wgllwgll_xz_3D
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NGLLZ) :: wgllwgll_yz_3D

! memory variables and standard linear solids for attenuation
  logical :: ATTENUATION
  logical :: COMPUTE_AND_STORE_STRAIN
  integer :: NSPEC_STRAIN_ONLY, NSPEC_ADJOINT
  integer :: NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_Kappa
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_Kappa) :: one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_Kappa) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
            R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_Kappa,N_SLS) :: R_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: &
            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_Kappa) :: epsilondev_trace
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: epsilon_trace_over_3

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store

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

  ! C-PML absorbing boundary conditions
  logical :: PML_CONDITIONS

  ! CPML adjoint
  logical :: backward_simulation

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  ! attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_att,tempx2_att,tempx3_att,tempy1_att,tempy2_att,tempy3_att,tempz1_att,tempz2_att,tempz3_att
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att,dummyy_loc_att,dummyz_loc_att
  real(kind=CUSTOM_REAL) :: duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att
  real(kind=CUSTOM_REAL) :: duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att

  ! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc,epsilondev_xx_loc, &
    epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc

  real(kind=CUSTOM_REAL) :: R_xx_val1,R_yy_val1,R_xx_val2,R_yy_val2,R_xx_val3,R_yy_val3, &
                            R_trace_val1,R_trace_val2,R_trace_val3
  real(kind=CUSTOM_REAL) :: factor_loc,alphaval_loc,betaval_loc,gammaval_loc
  real(kind=CUSTOM_REAL) :: Sn,Snp1
  real(kind=CUSTOM_REAL) :: templ

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  integer :: i_SLS,imodulo_N_SLS
  integer :: ispec,iglob,ispec_p,num_elements

  real(kind=CUSTOM_REAL) :: eta

  real(kind=CUSTOM_REAL) :: epsilon_trace,epsilon_trace_squared
  real(kind=CUSTOM_REAL) :: mul_plus_A_over_4,lambdal_plus_B,lambdal_over_two_plus_B_over_2

  ! C-PML absorbing boundary conditions
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_att_new,tempx2_att_new,tempx3_att_new, &
            tempy1_att_new,tempy2_att_new,tempy3_att_new, &
            tempz1_att_new,tempz2_att_new,tempz3_att_new

#ifdef FORCE_VECTORIZATION
! this will (purposely) give out-of-bound array accesses if run through range checking,
! thus use only for production runs with no bound checking
  integer :: ijk
#else
  integer :: i,j,k
#endif

  imodulo_N_SLS = mod(N_SLS,3)

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  do ispec_p = 1,num_elements


! arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
!
! hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
!                                           no counts for operations on indices in do-loops (?)
!
!                                           counts accesses to global memory, but no shared/cache memory or register loads/stores
!                                           float/real has 4 bytes

! hand-counts: floating-point operations FLOP, DRAM accesses in BYTES
!              for "simplest kernel" (isotropic without attenuation, dynamic fault, etc.)
!              and for single element, assuming NGLLX == NGLLY == NGLLZ == 5

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! adjoint simulations: moho kernel
    if (SIMULATION_TYPE == 3 .and. SAVE_MOHO_MESH) then
      if (is_moho_top(ispec)) then
        ispec2D_moho_top = ispec2D_moho_top + 1
      else if (is_moho_bot(ispec)) then
        ispec2D_moho_bot = ispec2D_moho_bot + 1
      endif
    endif ! adjoint

! counts:
! 0 FLOP
!
! 1 float = 4 BYTE

    ! Kelvin Voigt damping: artificial viscosity around dynamic faults

    ! stores displacment values in local array
    if (allocated(Kelvin_Voigt_eta)) then
      eta = Kelvin_Voigt_eta(ispec)
      if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
        if (is_CPML(ispec) .and. eta /= 0._CUSTOM_REAL) then
          stop 'you cannot put fault in PML region'
        endif
      endif
      DO_LOOP_IJK
        iglob = ibool(INDEX_IJK,ispec)
        dummyx_loc(INDEX_IJK) = displ(1,iglob) + eta * veloc(1,iglob)
        dummyy_loc(INDEX_IJK) = displ(2,iglob) + eta * veloc(2,iglob)
        dummyz_loc(INDEX_IJK) = displ(3,iglob) + eta * veloc(3,iglob)
      ENDDO_LOOP_IJK
    else
      DO_LOOP_IJK
        iglob = ibool(INDEX_IJK,ispec)
        dummyx_loc(INDEX_IJK) = displ(1,iglob)
        dummyy_loc(INDEX_IJK) = displ(2,iglob)
        dummyz_loc(INDEX_IJK) = displ(3,iglob)
      ENDDO_LOOP_IJK
    endif

! counts:
! + 0 FLOP
!
! + NGLLX * NGLLY * NGLLZ * ( 1 + 3) float = 2000 BYTE

    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        ! In backward_simulation involved in SIMULATION_TYPE == 3,
        ! we only use the stored value on edge of PML interface.
        ! Thus no compuation need to do in PML region in this case.
        if (.not. backward_simulation) then
          DO_LOOP_IJK
            iglob = ibool(INDEX_IJK,ispec)
            dummyx_loc_att(INDEX_IJK) = displ_old(1,iglob)
            dummyy_loc_att(INDEX_IJK) = displ_old(2,iglob)
            dummyz_loc_att(INDEX_IJK) = displ_old(3,iglob)
            dummyx_loc_att_new(INDEX_IJK) = displ_new(1,iglob)
            dummyy_loc_att_new(INDEX_IJK) = displ_new(2,iglob)
            dummyz_loc_att_new(INDEX_IJK) = displ_new(3,iglob)
          ENDDO_LOOP_IJK
        endif
      endif
    endif

    ! use first order Taylor expansion of displacement for local storage of stresses
    ! at this current time step, to fix attenuation in a consistent way
    if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
      if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
        ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
        ! because array is_CPML() is unallocated when PML_CONDITIONS is false
        if (.not. is_CPML(ispec)) then
          ! Keeping in mind, currently we implement a PML derived based on elastic wave equation
          ! for visco-elastic wave simulation. Thus it is not a PML anymore,
          ! because the solution of PML equation derived based on elastic wave equation is not perfectly mathced
          ! with solution of visco-elastic wave equation along the PML interface.
          ! you can seen PML in this case as a sponge layer.
          ! Due to limited numerical experiment, that PML implementation will work in case Q > 70.
          DO_LOOP_IJK
            iglob = ibool(INDEX_IJK,ispec)
            dummyx_loc_att(INDEX_IJK) = deltat * veloc(1,iglob)
            dummyy_loc_att(INDEX_IJK) = deltat * veloc(2,iglob)
            dummyz_loc_att(INDEX_IJK) = deltat * veloc(3,iglob)
          ENDDO_LOOP_IJK
        endif
      else
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          dummyx_loc_att(INDEX_IJK) = deltat * veloc(1,iglob)
          dummyy_loc_att(INDEX_IJK) = deltat * veloc(2,iglob)
          dummyz_loc_att(INDEX_IJK) = deltat * veloc(3,iglob)
        ENDDO_LOOP_IJK
      endif
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for tempx1,..
    call mxm_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)

! counts:
! + m1 * m2 * 3 * 9 = 5 * 25 * 3 * 9 = 3375 FLOP
!
! + m1 * 5 float = 100 BYTE  (hprime_xx once, assuming B3_** in cache)

    ! computes 2. matrix multiplication for tempx2,..
    call mxm_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,NGLLX)

! counts:
! + m1 * m1 * NGLLX * 3 * 9 = 5 * 5 * 5 * 3 * 9 = 3375 FLOP
!
! + m1 * 5 float = 100 BYTE  (hprime_xxT once, assuming dummy*_** in cache)

    ! computes 3. matrix multiplication for tempx1,..
    call mxm_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)

! counts:
! + m1 * m2 * 3 * 9 = 5 * 25 * 3 * 9 = 3375 FLOP
!
! + 0 BYTE  (assuming A3_**, hprime_xxT in cache)

    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        if (.not. backward_simulation) then
          ! computes 1. matrix multiplication for tempx1,..
          call mxm_3comp_singleA(hprime_xx,m1,dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                                 tempx1_att,tempy1_att,tempz1_att,m2)
          ! computes 2. matrix multiplication for tempx2,..
          call mxm_3comp_3dmat_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m1,hprime_xxT,m1, &
                                       tempx2_att,tempy2_att,tempz2_att,NGLLX)
          ! computes 3. matrix multiplication for tempx1,..
          call mxm_3comp_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m2,hprime_xxT, &
                                 tempx3_att,tempy3_att,tempz3_att,m1)

          ! computes 1. matrix multiplication for tempx1,..
          call mxm_3comp_singleA(hprime_xx,m1,dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new, &
                                 tempx1_att_new,tempy1_att_new,tempz1_att_new,m2)
          ! computes 2. matrix multiplication for tempx2,..
          call mxm_3comp_3dmat_singleB(dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new,m1,hprime_xxT,m1, &
                                       tempx2_att_new,tempy2_att_new,tempz2_att_new,NGLLX)
          ! computes 3. matrix multiplication for tempx1,..
          call mxm_3comp_singleB(dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new,m2,hprime_xxT, &
                                 tempx3_att_new,tempy3_att_new,tempz3_att_new,m1)

        endif
      endif
    endif

    if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
      ! it is noteworthy here that if ATTENUATION == .true., COMPUTE_AND_STORE_STRAIN == .true.
      if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
        ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
        ! because array is_CPML() is unallocated when PML_CONDITIONS is false
        if (.not. is_CPML(ispec)) then
          ! temporary variables used for fixing attenuation in a consistent way

          ! computes 1. matrix multiplication for tempx1,..
          call mxm_3comp_singleA(hprime_xx,m1,dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                                 tempx1_att,tempy1_att,tempz1_att,m2)

          tempx1_att(:,:,:) = tempx1_att(:,:,:) + tempx1(:,:,:)
          tempy1_att(:,:,:) = tempy1_att(:,:,:) + tempy1(:,:,:)
          tempz1_att(:,:,:) = tempz1_att(:,:,:) + tempz1(:,:,:)

          ! computes 2. matrix multiplication for tempx2,..
          call mxm_3comp_3dmat_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m1,hprime_xxT,m1, &
                                       tempx2_att,tempy2_att,tempz2_att,NGLLX)

          tempx2_att(:,:,:) = tempx2_att(:,:,:) + tempx2(:,:,:)
          tempy2_att(:,:,:) = tempy2_att(:,:,:) + tempy2(:,:,:)
          tempz2_att(:,:,:) = tempz2_att(:,:,:) + tempz2(:,:,:)

          ! computes 3. matrix multiplication for tempx1,..
          call mxm_3comp_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m2,hprime_xxT, &
                                 tempx3_att,tempy3_att,tempz3_att,m1)

          tempx3_att(:,:,:) = tempx3_att(:,:,:) + tempx3(:,:,:)
          tempy3_att(:,:,:) = tempy3_att(:,:,:) + tempy3(:,:,:)
          tempz3_att(:,:,:) = tempz3_att(:,:,:) + tempz3(:,:,:)
        endif
      else
        ! temporary variables used for fixing attenuation in a consistent way

        ! computes 1. matrix multiplication for tempx1,..
        call mxm_3comp_singleA(hprime_xx,m1,dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                               tempx1_att,tempy1_att,tempz1_att,m2)

        tempx1_att(:,:,:) = tempx1_att(:,:,:) + tempx1(:,:,:)
        tempy1_att(:,:,:) = tempy1_att(:,:,:) + tempy1(:,:,:)
        tempz1_att(:,:,:) = tempz1_att(:,:,:) + tempz1(:,:,:)

        ! computes 2. matrix multiplication for tempx2,..
        call mxm_3comp_3dmat_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m1,hprime_xxT,m1, &
                                     tempx2_att,tempy2_att,tempz2_att,NGLLX)

        tempx2_att(:,:,:) = tempx2_att(:,:,:) + tempx2(:,:,:)
        tempy2_att(:,:,:) = tempy2_att(:,:,:) + tempy2(:,:,:)
        tempz2_att(:,:,:) = tempz2_att(:,:,:) + tempz2(:,:,:)

        ! computes 3. matrix multiplication for tempx1,..
        call mxm_3comp_singleB(dummyx_loc_att,dummyy_loc_att,dummyz_loc_att,m2,hprime_xxT, &
                               tempx3_att,tempy3_att,tempz3_att,m1)

        tempx3_att(:,:,:) = tempx3_att(:,:,:) + tempx3(:,:,:)
        tempy3_att(:,:,:) = tempy3_att(:,:,:) + tempy3(:,:,:)
        tempz3_att(:,:,:) = tempz3_att(:,:,:) + tempz3(:,:,:)
      endif
    endif

    !
    ! computes either isotropic or fully anisotropic (visco-)elastic elements
    !
    DO_LOOP_IJK

      ! get derivatives of ux, uy and uz with respect to x, y and z
      xixl = xix(INDEX_IJK,ispec)
      xiyl = xiy(INDEX_IJK,ispec)
      xizl = xiz(INDEX_IJK,ispec)
      etaxl = etax(INDEX_IJK,ispec)
      etayl = etay(INDEX_IJK,ispec)
      etazl = etaz(INDEX_IJK,ispec)
      gammaxl = gammax(INDEX_IJK,ispec)
      gammayl = gammay(INDEX_IJK,ispec)
      gammazl = gammaz(INDEX_IJK,ispec)
      jacobianl = jacobian(INDEX_IJK,ispec)

! counts:
! + 0 FLOP
!
! + NGLLX * NGLLY * NGLLZ * 10 float = 5000 BYTE  (assuming A3_**, hprime_xxT in cache)

      duxdxl = xixl * tempx1(INDEX_IJK) + etaxl * tempx2(INDEX_IJK) + gammaxl * tempx3(INDEX_IJK)
      duxdyl = xiyl * tempx1(INDEX_IJK) + etayl * tempx2(INDEX_IJK) + gammayl * tempx3(INDEX_IJK)
      duxdzl = xizl * tempx1(INDEX_IJK) + etazl * tempx2(INDEX_IJK) + gammazl * tempx3(INDEX_IJK)

      duydxl = xixl * tempy1(INDEX_IJK) + etaxl * tempy2(INDEX_IJK) + gammaxl * tempy3(INDEX_IJK)
      duydyl = xiyl * tempy1(INDEX_IJK) + etayl * tempy2(INDEX_IJK) + gammayl * tempy3(INDEX_IJK)
      duydzl = xizl * tempy1(INDEX_IJK) + etazl * tempy2(INDEX_IJK) + gammazl * tempy3(INDEX_IJK)

      duzdxl = xixl * tempz1(INDEX_IJK) + etaxl * tempz2(INDEX_IJK) + gammaxl * tempz3(INDEX_IJK)
      duzdyl = xiyl * tempz1(INDEX_IJK) + etayl * tempz2(INDEX_IJK) + gammayl * tempz3(INDEX_IJK)
      duzdzl = xizl * tempz1(INDEX_IJK) + etazl * tempz2(INDEX_IJK) + gammazl * tempz3(INDEX_IJK)

      ! stores derivatives of ux, uy and uz with respect to x, y and z
      if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
        if (is_CPML(ispec)) then
          ! In backward_simulation involved in SIMULATION_TYPE == 3,
          ! we only use the stored value on edge of PML interface.
          ! Thus no compuation need to do in PML region in this case.
          if (.not. backward_simulation) then
            PML_dux_dxl(INDEX_IJK) = duxdxl
            PML_dux_dyl(INDEX_IJK) = duxdyl
            PML_dux_dzl(INDEX_IJK) = duxdzl

            PML_duy_dxl(INDEX_IJK) = duydxl
            PML_duy_dyl(INDEX_IJK) = duydyl
            PML_duy_dzl(INDEX_IJK) = duydzl

            PML_duz_dxl(INDEX_IJK) = duzdxl
            PML_duz_dyl(INDEX_IJK) = duzdyl
            PML_duz_dzl(INDEX_IJK) = duzdzl

            PML_dux_dxl_old(INDEX_IJK) = &
                xixl * tempx1_att(INDEX_IJK) + etaxl * tempx2_att(INDEX_IJK) + gammaxl * tempx3_att(INDEX_IJK)
            PML_dux_dyl_old(INDEX_IJK) = &
                xiyl * tempx1_att(INDEX_IJK) + etayl * tempx2_att(INDEX_IJK) + gammayl * tempx3_att(INDEX_IJK)
            PML_dux_dzl_old(INDEX_IJK) = &
                xizl * tempx1_att(INDEX_IJK) + etazl * tempx2_att(INDEX_IJK) + gammazl * tempx3_att(INDEX_IJK)

            PML_duy_dxl_old(INDEX_IJK) = &
                xixl * tempy1_att(INDEX_IJK) + etaxl * tempy2_att(INDEX_IJK) + gammaxl * tempy3_att(INDEX_IJK)
            PML_duy_dyl_old(INDEX_IJK) = &
                xiyl * tempy1_att(INDEX_IJK) + etayl * tempy2_att(INDEX_IJK) + gammayl * tempy3_att(INDEX_IJK)
            PML_duy_dzl_old(INDEX_IJK) = &
                xizl * tempy1_att(INDEX_IJK) + etazl * tempy2_att(INDEX_IJK) + gammazl * tempy3_att(INDEX_IJK)

            PML_duz_dxl_old(INDEX_IJK) = &
                xixl * tempz1_att(INDEX_IJK) + etaxl * tempz2_att(INDEX_IJK) + gammaxl * tempz3_att(INDEX_IJK)
            PML_duz_dyl_old(INDEX_IJK) = &
                xiyl * tempz1_att(INDEX_IJK) + etayl * tempz2_att(INDEX_IJK) + gammayl * tempz3_att(INDEX_IJK)
            PML_duz_dzl_old(INDEX_IJK) = &
                 xizl * tempz1_att(INDEX_IJK) + etazl * tempz2_att(INDEX_IJK) + gammazl * tempz3_att(INDEX_IJK)

            PML_dux_dxl_new(INDEX_IJK) = &
                xixl * tempx1_att_new(INDEX_IJK) + etaxl * tempx2_att_new(INDEX_IJK) + gammaxl * tempx3_att_new(INDEX_IJK)
            PML_dux_dyl_new(INDEX_IJK) = &
                xiyl * tempx1_att_new(INDEX_IJK) + etayl * tempx2_att_new(INDEX_IJK) + gammayl * tempx3_att_new(INDEX_IJK)
            PML_dux_dzl_new(INDEX_IJK) = &
                xizl * tempx1_att_new(INDEX_IJK) + etazl * tempx2_att_new(INDEX_IJK) + gammazl * tempx3_att_new(INDEX_IJK)

            PML_duy_dxl_new(INDEX_IJK) = &
                xixl * tempy1_att_new(INDEX_IJK) + etaxl * tempy2_att_new(INDEX_IJK) + gammaxl * tempy3_att_new(INDEX_IJK)
            PML_duy_dyl_new(INDEX_IJK) = &
                xiyl * tempy1_att_new(INDEX_IJK) + etayl * tempy2_att_new(INDEX_IJK) + gammayl * tempy3_att_new(INDEX_IJK)
            PML_duy_dzl_new(INDEX_IJK) = &
                xizl * tempy1_att_new(INDEX_IJK) + etazl * tempy2_att_new(INDEX_IJK) + gammazl * tempy3_att_new(INDEX_IJK)

            PML_duz_dxl_new(INDEX_IJK) = &
                xixl * tempz1_att_new(INDEX_IJK) + etaxl * tempz2_att_new(INDEX_IJK) + gammaxl * tempz3_att_new(INDEX_IJK)
            PML_duz_dyl_new(INDEX_IJK) = &
                xiyl * tempz1_att_new(INDEX_IJK) + etayl * tempz2_att_new(INDEX_IJK) + gammayl * tempz3_att_new(INDEX_IJK)
            PML_duz_dzl_new(INDEX_IJK) = &
                xizl * tempz1_att_new(INDEX_IJK) + etazl * tempz2_att_new(INDEX_IJK) + gammazl * tempz3_att_new(INDEX_IJK)
          endif

          if (COMPUTE_AND_STORE_STRAIN) then
            templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
            if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
            if (FULL_ATTENUATION_SOLID) epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
            epsilondev_xx_loc(INDEX_IJK) = duxdxl - templ
            epsilondev_yy_loc(INDEX_IJK) = duydyl - templ
            epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duxdyl + duydxl)
            epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdxl + duxdzl)
            epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * (duzdyl + duydzl)
          endif

        endif
      endif

! counts:
! + NGLLX * NGLLY * NGLLZ * 9 * 5 = 5625 FLOP
!
! + 0 BYTE  (assuming temp*_** in cache)

      ! save strain on the Moho boundary
      if (SIMULATION_TYPE == 3 .and. SAVE_MOHO_MESH) then
        if (is_moho_top(ispec)) then
          dsdx_top(1,1,INDEX_IJK,ispec2D_moho_top) = duxdxl
          dsdx_top(1,2,INDEX_IJK,ispec2D_moho_top) = duxdyl
          dsdx_top(1,3,INDEX_IJK,ispec2D_moho_top) = duxdzl
          dsdx_top(2,1,INDEX_IJK,ispec2D_moho_top) = duydxl
          dsdx_top(2,2,INDEX_IJK,ispec2D_moho_top) = duydyl
          dsdx_top(2,3,INDEX_IJK,ispec2D_moho_top) = duydzl
          dsdx_top(3,1,INDEX_IJK,ispec2D_moho_top) = duzdxl
          dsdx_top(3,2,INDEX_IJK,ispec2D_moho_top) = duzdyl
          dsdx_top(3,3,INDEX_IJK,ispec2D_moho_top) = duzdzl
        else if (is_moho_bot(ispec)) then
          dsdx_bot(1,1,INDEX_IJK,ispec2D_moho_bot) = duxdxl
          dsdx_bot(1,2,INDEX_IJK,ispec2D_moho_bot) = duxdyl
          dsdx_bot(1,3,INDEX_IJK,ispec2D_moho_bot) = duxdzl
          dsdx_bot(2,1,INDEX_IJK,ispec2D_moho_bot) = duydxl
          dsdx_bot(2,2,INDEX_IJK,ispec2D_moho_bot) = duydyl
          dsdx_bot(2,3,INDEX_IJK,ispec2D_moho_bot) = duydzl
          dsdx_bot(3,1,INDEX_IJK,ispec2D_moho_bot) = duzdxl
          dsdx_bot(3,2,INDEX_IJK,ispec2D_moho_bot) = duzdyl
          dsdx_bot(3,3,INDEX_IJK,ispec2D_moho_bot) = duzdzl
        endif
      endif

      ! precompute some sums to save CPU time
      duxdxl_plus_duydyl = duxdxl + duydyl
      duxdxl_plus_duzdzl = duxdxl + duzdzl
      duydyl_plus_duzdzl = duydyl + duzdzl
      duxdyl_plus_duydxl = duxdyl + duydxl
      duzdxl_plus_duxdzl = duzdxl + duxdzl
      duzdyl_plus_duydzl = duzdyl + duydzl

! counts:
! + NGLLX * NGLLY * NGLLZ * 6 * 1 = 750 FLOP
!
! + 0 BYTE  (assuming registers)

      if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
            ! temporary variables used for fixing attenuation in a consistent way
        if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
          ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
          ! because array is_CPML() is unallocated when PML_CONDITIONS is false
          if (.not. is_CPML(ispec)) then
            ! temporary variables used for fixing attenuation in a consistent way
            duxdxl_att = xixl * tempx1_att(INDEX_IJK) + etaxl * tempx2_att(INDEX_IJK) + gammaxl * tempx3_att(INDEX_IJK)
            duxdyl_att = xiyl * tempx1_att(INDEX_IJK) + etayl * tempx2_att(INDEX_IJK) + gammayl * tempx3_att(INDEX_IJK)
            duxdzl_att = xizl * tempx1_att(INDEX_IJK) + etazl * tempx2_att(INDEX_IJK) + gammazl * tempx3_att(INDEX_IJK)

            duydxl_att = xixl * tempy1_att(INDEX_IJK) + etaxl * tempy2_att(INDEX_IJK) + gammaxl * tempy3_att(INDEX_IJK)
            duydyl_att = xiyl * tempy1_att(INDEX_IJK) + etayl * tempy2_att(INDEX_IJK) + gammayl * tempy3_att(INDEX_IJK)
            duydzl_att = xizl * tempy1_att(INDEX_IJK) + etazl * tempy2_att(INDEX_IJK) + gammazl * tempy3_att(INDEX_IJK)

            duzdxl_att = xixl * tempz1_att(INDEX_IJK) + etaxl * tempz2_att(INDEX_IJK) + gammaxl * tempz3_att(INDEX_IJK)
            duzdyl_att = xiyl * tempz1_att(INDEX_IJK) + etayl * tempz2_att(INDEX_IJK) + gammayl * tempz3_att(INDEX_IJK)
            duzdzl_att = xizl * tempz1_att(INDEX_IJK) + etazl * tempz2_att(INDEX_IJK) + gammazl * tempz3_att(INDEX_IJK)

            ! precompute some sums to save CPU time
            duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
            duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
            duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

            ! compute deviatoric strain
            templ = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
            if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
            if (FULL_ATTENUATION_SOLID) epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
            epsilondev_xx_loc(INDEX_IJK) = duxdxl_att - templ
            epsilondev_yy_loc(INDEX_IJK) = duydyl_att - templ
            epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl_att
            epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl_att
            epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl_att
          endif
        else
          ! temporary variables used for fixing attenuation in a consistent way
          duxdxl_att = xixl * tempx1_att(INDEX_IJK) + etaxl * tempx2_att(INDEX_IJK) + gammaxl * tempx3_att(INDEX_IJK)
          duxdyl_att = xiyl * tempx1_att(INDEX_IJK) + etayl * tempx2_att(INDEX_IJK) + gammayl * tempx3_att(INDEX_IJK)
          duxdzl_att = xizl * tempx1_att(INDEX_IJK) + etazl * tempx2_att(INDEX_IJK) + gammazl * tempx3_att(INDEX_IJK)

          duydxl_att = xixl * tempy1_att(INDEX_IJK) + etaxl * tempy2_att(INDEX_IJK) + gammaxl * tempy3_att(INDEX_IJK)
          duydyl_att = xiyl * tempy1_att(INDEX_IJK) + etayl * tempy2_att(INDEX_IJK) + gammayl * tempy3_att(INDEX_IJK)
          duydzl_att = xizl * tempy1_att(INDEX_IJK) + etazl * tempy2_att(INDEX_IJK) + gammazl * tempy3_att(INDEX_IJK)

          duzdxl_att = xixl * tempz1_att(INDEX_IJK) + etaxl * tempz2_att(INDEX_IJK) + gammaxl * tempz3_att(INDEX_IJK)
          duzdyl_att = xiyl * tempz1_att(INDEX_IJK) + etayl * tempz2_att(INDEX_IJK) + gammayl * tempz3_att(INDEX_IJK)
          duzdzl_att = xizl * tempz1_att(INDEX_IJK) + etazl * tempz2_att(INDEX_IJK) + gammazl * tempz3_att(INDEX_IJK)

          ! precompute some sums to save CPU time
          duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
          duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
          duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

          ! compute deviatoric strain
          templ = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
          if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
          if (FULL_ATTENUATION_SOLID) epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
          epsilondev_xx_loc(INDEX_IJK) = duxdxl_att - templ
          epsilondev_yy_loc(INDEX_IJK) = duydyl_att - templ
          epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl_att
          epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl_att
          epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl_att
        endif
      else
        ! computes deviatoric strain attenuation and/or for kernel calculations
        if (COMPUTE_AND_STORE_STRAIN) then
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          if (SIMULATION_TYPE == 3) epsilon_trace_over_3(INDEX_IJK,ispec) = templ
          if (FULL_ATTENUATION_SOLID) epsilondev_trace_loc(INDEX_IJK) = 3._CUSTOM_REAL * templ
          epsilondev_xx_loc(INDEX_IJK) = duxdxl - templ
          epsilondev_yy_loc(INDEX_IJK) = duydyl - templ
          epsilondev_xy_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          epsilondev_xz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          epsilondev_yz_loc(INDEX_IJK) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
        endif
      endif

      kappal = kappastore(INDEX_IJK,ispec)
      mul = mustore(INDEX_IJK,ispec)

! counts:
! + 0 FLOP
!
! + NGLLX * NGLLY * NGLLZ * 2 float = 1000 BYTE

      ! attenuation
      if (ATTENUATION) then
        if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
          ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
          ! because array is_CPML() is unallocated when PML_CONDITIONS is false
          if (.not. is_CPML(ispec)) then
            ! use unrelaxed parameters if attenuation
            mul  = mul * one_minus_sum_beta(INDEX_IJK,ispec)
            if (FULL_ATTENUATION_SOLID) kappal  = kappal * one_minus_sum_beta_kappa(INDEX_IJK,ispec)
          endif
        else
          ! use unrelaxed parameters if attenuation
          mul  = mul * one_minus_sum_beta(INDEX_IJK,ispec)
          if (FULL_ATTENUATION_SOLID) kappal  = kappal * one_minus_sum_beta_kappa(INDEX_IJK,ispec)
        endif
      endif

      ! full anisotropic case, stress calculations
      if (ANISOTROPY) then
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
        lambdal = lambdalplus2mul - 2._CUSTOM_REAL*mul

        ! compute stress sigma
        if (.not. MAKE_HOOKE_LAW_WEAKLY_NONLINEAR) then

          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

        else

        ! in the case of weakly nonlinear materials the stress tensor becomes non-symmetric, see for instance equation (A9) in
        ! D. L. Johnson, S. Kostek and A. N. Norris, Nonlinear tube waves, J. Acoust. Soc. Am. vol. 96, p. 1829-1843 (1994).

        ! TODO: the current (explicit Newmark) time scheme will likely *NOT* work for such a nonlinear material,
        ! the CFL in such a material will vary with time and thus the time step should be adapted while the simulation is running.
        ! To validate this with a better time scheme the easiest or best thing to do would be to perform some tests in SPECFEM1D
        ! and compare the results against a known solution for instance for the 1D Riemann problem
        ! and/or the 1D shock-wave tube problem.

          epsilon_trace = duxdxl + duydyl + duzdzl
          epsilon_trace_squared = epsilon_trace * epsilon_trace

          mul_plus_A_over_4 = mul + A_over_4
          lambdal_plus_B = lambdal + B
          lambdal_over_two_plus_B_over_2 = 0.5_CUSTOM_REAL * lambdal + B_over_2

          sigma_xx = lambdal * epsilon_trace + mul * (duxdxl + duxdxl) + B * epsilon_trace * duxdxl + &
                     lambdal_plus_B * epsilon_trace * duxdxl + C * epsilon_trace_squared + A_over_4 * &
                     (duxdxl * duxdxl + duxdyl * duydxl + duxdzl * duzdxl) + mul_plus_A_over_4 * &
                     (duxdxl * duxdxl + duxdxl * duxdxl + duxdxl * duxdxl + duxdyl * duxdyl + &
                     duydxl * duydxl + duxdyl * duydxl + duxdzl * duxdzl + duzdxl * duzdxl + &
                     duxdzl * duzdxl) + lambdal_over_two_plus_B_over_2 * (duxdxl * duxdxl + &
                     duydxl * duydxl + duzdxl * duzdxl + duxdyl * duxdyl + duydyl * duydyl + &
                     duzdyl * duzdyl + duxdzl * duxdzl + duydzl * duydzl + duzdzl * duzdzl) + &
                     B_over_2 * (duxdxl * duxdxl + duydxl * duxdyl + duzdxl * duxdzl + &
                     duxdyl * duydxl + duydyl * duydyl + duzdyl * duydzl + &
                     duxdzl * duzdxl + duydzl * duzdyl + duzdzl * duzdzl)

          sigma_xy = mul * (duxdyl + duydxl) + B * epsilon_trace * duydxl + lambdal_plus_B * &
                     epsilon_trace * duxdyl + A_over_4 * (duydxl * duxdxl + duydyl * duydxl + &
                     duydzl * duzdxl) + mul_plus_A_over_4 * (duxdxl * duydxl + duxdxl * duxdyl + &
                     duxdxl * duxdyl + duxdyl * duydyl + duydxl * duydyl + duxdyl * duydyl + &
                     duxdzl * duydzl + duzdxl * duzdyl + duxdzl * duzdyl)

          sigma_xz = mul * (duxdzl + duzdxl) + B * epsilon_trace * duzdxl + lambdal_plus_B * &
                     epsilon_trace * duxdzl + A_over_4 * (duzdxl * duxdxl + duzdyl * duydxl + &
                     duzdzl * duzdxl) + mul_plus_A_over_4 * (duxdxl * duzdxl + duxdxl * duxdzl + &
                     duxdxl * duxdzl + duxdyl * duzdyl + duydxl * duydzl + duxdyl * duydzl + &
                     duxdzl * duzdzl + duzdxl * duzdzl + duxdzl * duzdzl)

          sigma_yx = mul * (duydxl + duxdyl) + B * epsilon_trace * duxdyl + lambdal_plus_B * &
                     epsilon_trace * duydxl + A_over_4 * (duxdxl * duxdyl + duxdyl * duydyl + &
                     duxdzl * duzdyl) + mul_plus_A_over_4 * (duydxl * duxdxl + duxdyl * duxdxl + &
                     duydxl * duxdxl + duydyl * duxdyl + duydyl * duydxl + duydyl * duydxl + &
                     duydzl * duxdzl + duzdyl * duzdxl + duydzl * duzdxl)

          sigma_yy = lambdal * epsilon_trace + mul * (duydyl + duydyl) + B * epsilon_trace * duydyl + &
                     lambdal_plus_B * epsilon_trace * duydyl + C * epsilon_trace_squared + A_over_4 * &
                     (duydxl * duxdyl + duydyl * duydyl + duydzl * duzdyl) + mul_plus_A_over_4 * &
                     (duydxl * duydxl + duxdyl * duxdyl + duydxl * duxdyl + duydyl * duydyl + &
                     duydyl * duydyl + duydyl * duydyl + duydzl * duydzl + duzdyl * duzdyl + &
                     duydzl * duzdyl) + lambdal_over_two_plus_B_over_2 * (duxdxl * duxdxl + &
                     duydxl * duydxl + duzdxl * duzdxl + duxdyl * duxdyl + duydyl * duydyl + &
                     duzdyl * duzdyl + duxdzl * duxdzl + duydzl * duydzl + duzdzl * duzdzl) + &
                     B_over_2 * (duxdxl * duxdxl + duydxl * duxdyl + duzdxl * duxdzl + &
                     duxdyl * duydxl + duydyl * duydyl + duzdyl * duydzl + &
                     duxdzl * duzdxl + duydzl * duzdyl + duzdzl * duzdzl)

          sigma_yz = mul * (duydzl + duzdyl) + B * epsilon_trace * duzdyl + lambdal_plus_B * &
                     epsilon_trace * duydzl + A_over_4 * (duzdxl * duxdyl + duzdyl * duydyl + &
                     duzdzl * duzdyl) + mul_plus_A_over_4 * (duydxl * duzdxl + duxdyl * duxdzl + &
                     duydxl * duxdzl + duydyl * duzdyl + duydyl * duydzl + duydyl * duydzl + &
                     duydzl * duzdzl + duzdyl * duzdzl + duydzl * duzdzl)

          sigma_zx = mul * (duzdxl + duxdzl) + B * epsilon_trace * duxdzl + lambdal_plus_B * &
                     epsilon_trace * duzdxl + A_over_4 * (duxdxl * duxdzl + duxdyl * duydzl + &
                     duxdzl * duzdzl) + mul_plus_A_over_4 * (duzdxl * duxdxl + duxdzl * duxdxl + &
                     duzdxl * duxdxl + duzdyl * duxdyl + duydzl * duydxl + duzdyl * duydxl + &
                     duzdzl * duxdzl + duzdzl * duzdxl + duzdzl * duzdxl)

          sigma_zy = mul * (duzdyl + duydzl) + B * epsilon_trace * duydzl + lambdal_plus_B * &
                     epsilon_trace * duzdyl + A_over_4 * (duydxl * duxdzl + duydyl * duydzl + &
                     duydzl * duzdzl) + mul_plus_A_over_4 * (duzdxl * duydxl + duxdzl * duxdyl + &
                     duzdxl * duxdyl + duzdyl * duydyl + duydzl * duydyl + duzdyl * duydyl + &
                     duzdzl * duydzl + duzdzl * duzdyl + duzdzl * duzdyl)

          sigma_zz = lambdal * epsilon_trace + mul * (duzdzl + duzdzl) + B * epsilon_trace * &
                     duzdzl + lambdal_plus_B * epsilon_trace * duzdzl + C * epsilon_trace_squared + &
                     A_over_4 * (duzdxl * duxdzl + duzdyl * duydzl + duzdzl * duzdzl) + &
                     mul_plus_A_over_4 * (duzdxl * duzdxl + duxdzl * duxdzl + duzdxl * duxdzl + &
                     duzdyl * duzdyl + duydzl * duydzl + duzdyl * duydzl + duzdzl * duzdzl + &
                     duzdzl * duzdzl + duzdzl * duzdzl) + lambdal_over_two_plus_B_over_2 * (duxdxl * duxdxl + &
                     duydxl * duydxl + duzdxl * duzdxl + duxdyl * duxdyl + duydyl * duydyl + duzdyl * duzdyl + &
                     duxdzl * duxdzl + duydzl * duydzl + duzdzl * duzdzl) + B_over_2 * (duxdxl * duxdxl + &
                     duydxl * duxdyl + duzdxl * duxdzl + duxdyl * duydxl + duydyl * duydyl + duzdyl * duydzl + &
                     duxdzl * duzdxl + duydzl * duzdyl + duzdzl * duzdzl)

        endif

      endif ! of if (ANISOTROPY)

! counts:
! + NGLLX * NGLLY * NGLLZ * 16 =  2000 FLOP
!
! + 0 BYTE

      ! subtract memory variables if attenuation
      if (ATTENUATION) then
! way 1
!            do i_sls = 1,N_SLS
!              R_xx_val = R_xx(INDEX_IJK,ispec,i_sls)
!              R_yy_val = R_yy(INDEX_IJK,ispec,i_sls)
!              sigma_xx = sigma_xx - R_xx_val
!              sigma_yy = sigma_yy - R_yy_val
!              sigma_zz = sigma_zz + R_xx_val + R_yy_val
!              sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls)
!              sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls)
!              sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls)
!            enddo

! way 2
! note: this should help compilers to pipeline the code and make better use of the cache;
!          depending on compilers, it can further decrease the computation time by ~ 30%.
!          by default, N_SLS = 3, therefore we take steps of 3
        if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
          ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
          ! because array is_CPML() is unallocated when PML_CONDITIONS is false
          if (.not. is_CPML(ispec)) then
            if (imodulo_N_SLS >= 1) then
              do i_sls = 1,imodulo_N_SLS
                if (FULL_ATTENUATION_SOLID) then
                  R_trace_val1 = R_trace(INDEX_IJK,ispec,i_sls)
                else
                  R_trace_val1 = 0.
                endif
                R_xx_val1 = R_xx(INDEX_IJK,ispec,i_sls)
                R_yy_val1 = R_yy(INDEX_IJK,ispec,i_sls)
                sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
                sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
                sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
                sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls)
                sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls)
                sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls)
              enddo
            endif

            if (N_SLS >= imodulo_N_SLS+1) then
              do i_sls = imodulo_N_SLS+1,N_SLS,3
                if (FULL_ATTENUATION_SOLID) then
                  R_trace_val1 = R_trace(INDEX_IJK,ispec,i_sls)
                else
                  R_trace_val1 = 0.
                endif
                R_xx_val1 = R_xx(INDEX_IJK,ispec,i_sls)
                R_yy_val1 = R_yy(INDEX_IJK,ispec,i_sls)
                sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
                sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
                sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
                sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls)
                sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls)
                sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls)
                if (FULL_ATTENUATION_SOLID) then
                  R_trace_val2 = R_trace(INDEX_IJK,ispec,i_sls+1)
                else
                  R_trace_val2 = 0.
                endif
                R_xx_val2 = R_xx(INDEX_IJK,ispec,i_sls+1)
                R_yy_val2 = R_yy(INDEX_IJK,ispec,i_sls+1)
                sigma_xx = sigma_xx - R_xx_val2 - R_trace_val2
                sigma_yy = sigma_yy - R_yy_val2 - R_trace_val2
                sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2 - R_trace_val2
                sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls+1)
                sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls+1)
                sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls+1)

                if (FULL_ATTENUATION_SOLID) then
                  R_trace_val3 = R_trace(INDEX_IJK,ispec,i_sls+2)
                else
                  R_trace_val3 = 0.
                endif
                R_xx_val3 = R_xx(INDEX_IJK,ispec,i_sls+2)
                R_yy_val3 = R_yy(INDEX_IJK,ispec,i_sls+2)
                sigma_xx = sigma_xx - R_xx_val3 - R_trace_val3
                sigma_yy = sigma_yy - R_yy_val3 - R_trace_val3
                sigma_zz = sigma_zz + R_xx_val3 + R_yy_val3 - R_trace_val3
                sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls+2)
                sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls+2)
                sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls+2)
              enddo
            endif
          endif
        else
          if (imodulo_N_SLS >= 1) then
            do i_sls = 1,imodulo_N_SLS
              if (FULL_ATTENUATION_SOLID) then
                R_trace_val1 = R_trace(INDEX_IJK,ispec,i_sls)
              else
                R_trace_val1 = 0.
              endif
              R_xx_val1 = R_xx(INDEX_IJK,ispec,i_sls)
              R_yy_val1 = R_yy(INDEX_IJK,ispec,i_sls)
              sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
              sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
              sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
              sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls)
              sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls)
              sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls)
            enddo
          endif

          if (N_SLS >= imodulo_N_SLS+1) then
            do i_sls = imodulo_N_SLS+1,N_SLS,3
              if (FULL_ATTENUATION_SOLID) then
                R_trace_val1 = R_trace(INDEX_IJK,ispec,i_sls)
              else
                R_trace_val1 = 0.
              endif
              R_xx_val1 = R_xx(INDEX_IJK,ispec,i_sls)
              R_yy_val1 = R_yy(INDEX_IJK,ispec,i_sls)
              sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
              sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
              sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
              sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls)
              sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls)
              sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls)
              if (FULL_ATTENUATION_SOLID) then
                R_trace_val2 = R_trace(INDEX_IJK,ispec,i_sls+1)
              else
                R_trace_val2 = 0.
              endif
              R_xx_val2 = R_xx(INDEX_IJK,ispec,i_sls+1)
              R_yy_val2 = R_yy(INDEX_IJK,ispec,i_sls+1)
              sigma_xx = sigma_xx - R_xx_val2 - R_trace_val2
              sigma_yy = sigma_yy - R_yy_val2 - R_trace_val2
              sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2 - R_trace_val2
              sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls+1)
              sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls+1)
              sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls+1)

              if (FULL_ATTENUATION_SOLID) then
                R_trace_val3 = R_trace(INDEX_IJK,ispec,i_sls+2)
              else
                R_trace_val3 = 0.
              endif
              R_xx_val3 = R_xx(INDEX_IJK,ispec,i_sls+2)
              R_yy_val3 = R_yy(INDEX_IJK,ispec,i_sls+2)
              sigma_xx = sigma_xx - R_xx_val3 - R_trace_val3
              sigma_yy = sigma_yy - R_yy_val3 - R_trace_val3
              sigma_zz = sigma_zz + R_xx_val3 + R_yy_val3 - R_trace_val3
              sigma_xy = sigma_xy - R_xy(INDEX_IJK,ispec,i_sls+2)
              sigma_xz = sigma_xz - R_xz(INDEX_IJK,ispec,i_sls+2)
              sigma_yz = sigma_yz - R_yz(INDEX_IJK,ispec,i_sls+2)
            enddo
          endif
        endif

      endif

      ! define symmetric components of sigma
      if (.not. MAKE_HOOKE_LAW_WEAKLY_NONLINEAR) then
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz
      endif

      ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
      tempx1(INDEX_IJK) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
      tempy1(INDEX_IJK) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
      tempz1(INDEX_IJK) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

      tempx2(INDEX_IJK) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
      tempy2(INDEX_IJK) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
      tempz2(INDEX_IJK) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

      tempx3(INDEX_IJK) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
      tempy3(INDEX_IJK) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
      tempz3(INDEX_IJK) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

! counts:
! + NGLLX * NGLLY * NGLLZ * 9 * 6 = 6750 FLOP
!
! + NGLLX * NGLLY * NGLLZ * 9 float = 4500 BYTE (temp* stores)

    ENDDO_LOOP_IJK

    if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
      ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        ! In backward_simulation involved in SIMULATION_TYPE == 3,
        ! we only use the stored value on edge of PML interface.
        ! Thus no compuation need to do in PML region in this case.
        if (.not. backward_simulation) then
          ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
          ispec_CPML = spec_to_CPML(ispec)
          call pml_compute_memory_variables_elastic(ispec,ispec_CPML, &
                   tempx1,tempy1,tempz1, &
                   tempx2,tempy2,tempz2, &
                   tempx3,tempy3,tempz3, &
                   rmemory_dux_dxl_x, rmemory_duy_dyl_x, rmemory_duz_dzl_x, &
                   rmemory_dux_dyl_x, rmemory_dux_dzl_x, rmemory_duz_dxl_x, rmemory_duy_dxl_x, &
                   rmemory_dux_dxl_y, rmemory_duz_dzl_y, rmemory_duy_dyl_y, &
                   rmemory_duy_dxl_y, rmemory_duy_dzl_y, rmemory_duz_dyl_y, rmemory_dux_dyl_y, &
                   rmemory_dux_dxl_z, rmemory_duy_dyl_z, rmemory_duz_dzl_z, &
                   rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duy_dzl_z, rmemory_dux_dzl_z)

          ! calculates contribution from each C-PML element to update acceleration
          call pml_compute_accel_contribution_elastic(ispec,ispec_CPML,displ,veloc,rmemory_displ_elastic)
        endif
      endif
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    call mxm_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)

! counts:
! + m1 * m2 * 3 * 9 = 3375 FLOP
!
! + m1 * 5 float = 100 BYTE (hprimewgll_xxT once, assumes E3*, C1* in cache)

    ! computes 2. matrix multiplication for tempx2,..
    call mxm_3comp_3dmat_singleB(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,NGLLX)

! counts:
! + m1 * m1 * NGLLX * 3 * 9 = 3375 FLOP
!
! + m1 * 5 float = 100 BYTE (hprimewgll_xx once, assumes E3*, C1* in cache)

    ! computes 3. matrix multiplication for newtempx3,..
    call mxm_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)

! counts:
! + m1 * m2 * 3 * 9 = 3375 FLOP
!
! + 0 BYTE (assumes E1*, C1*, hprime* in cache)

    ! sums contributions
    DO_LOOP_IJK
      fac1 = wgllwgll_yz_3D(INDEX_IJK)
      fac2 = wgllwgll_xz_3D(INDEX_IJK)
      fac3 = wgllwgll_xy_3D(INDEX_IJK)

      ! sum contributions from each element to the global mesh using indirect addressing
      iglob = ibool(INDEX_IJK,ispec)
      accel(1,iglob) = accel(1,iglob) &
        - fac1*newtempx1(INDEX_IJK) - fac2*newtempx2(INDEX_IJK) - fac3*newtempx3(INDEX_IJK)
      accel(2,iglob) = accel(2,iglob) &
        - fac1*newtempy1(INDEX_IJK) - fac2*newtempy2(INDEX_IJK) - fac3*newtempy3(INDEX_IJK)
      accel(3,iglob) = accel(3,iglob) &
        - fac1*newtempz1(INDEX_IJK) - fac2*newtempz2(INDEX_IJK) - fac3*newtempz3(INDEX_IJK)

! counts:
! + NGLLX * NGLLY * NGLLZ * 3 * 6 = 2250 FLOP
!
! + NGLLX * NGLLY * 3 float = 300 BYTE (wgllwgll once)
! + NGLLX * NGLLY * NGLLZ * (1 + 3) float = 2000 BYTE (ibool & accel, assumes newtemp* in cache)

      !  update memory variables based upon the Runge-Kutta scheme
      if (ATTENUATION) then
        if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
          ! do not merge the line "if (is_CPML(ispec)) then" with the above if statement using an " .and. " statement
          ! because array is_CPML() is unallocated when PML_CONDITIONS is false
          if (.not. is_CPML(ispec)) then
            ! use Runge-Kutta scheme to march in time
            do i_sls = 1,N_SLS

              alphaval_loc = alphaval(i_sls)
              betaval_loc = betaval(i_sls)
              gammaval_loc = gammaval(i_sls)

              if (FULL_ATTENUATION_SOLID) then
                ! term in trace
                factor_loc = kappastore(INDEX_IJK,ispec) * factor_common_kappa(i_sls,INDEX_IJK,ispec)

                Sn   = factor_loc * epsilondev_trace(INDEX_IJK,ispec)
                Snp1   = factor_loc * epsilondev_trace_loc(INDEX_IJK)
                R_trace(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_trace(INDEX_IJK,ispec,i_sls) + &
                                  betaval_loc * Sn + gammaval_loc * Snp1
              endif

              ! term in xx yy zz xy xz yz
              factor_loc = mustore(INDEX_IJK,ispec) * factor_common(i_sls,INDEX_IJK,ispec)

              ! term in xx
              Sn   = factor_loc * epsilondev_xx(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_xx_loc(INDEX_IJK)
              R_xx(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xx(INDEX_IJK,ispec,i_sls) + &
                             betaval_loc * Sn + gammaval_loc * Snp1
              ! term in yy
              Sn   = factor_loc * epsilondev_yy(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_yy_loc(INDEX_IJK)
              R_yy(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_yy(INDEX_IJK,ispec,i_sls) + &
                             betaval_loc * Sn + gammaval_loc * Snp1
              ! term in zz not computed since zero trace
              ! term in xy
              Sn   = factor_loc * epsilondev_xy(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_xy_loc(INDEX_IJK)
              R_xy(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xy(INDEX_IJK,ispec,i_sls) + &
                             betaval_loc * Sn + gammaval_loc * Snp1
              ! term in xz
              Sn   = factor_loc * epsilondev_xz(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_xz_loc(INDEX_IJK)
              R_xz(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xz(INDEX_IJK,ispec,i_sls) + &
                             betaval_loc * Sn + gammaval_loc * Snp1
              ! term in yz
              Sn   = factor_loc * epsilondev_yz(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_yz_loc(INDEX_IJK)
              R_yz(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_yz(INDEX_IJK,ispec,i_sls) + &
                             betaval_loc * Sn + gammaval_loc * Snp1
            enddo   ! end of loop on memory variables
          endif
        else
          ! use Runge-Kutta scheme to march in time
          do i_sls = 1,N_SLS

            alphaval_loc = alphaval(i_sls)
            betaval_loc = betaval(i_sls)
            gammaval_loc = gammaval(i_sls)

            if (FULL_ATTENUATION_SOLID) then
              ! term in trace
              factor_loc = kappastore(INDEX_IJK,ispec) * factor_common_kappa(i_sls,INDEX_IJK,ispec)

              Sn   = factor_loc * epsilondev_trace(INDEX_IJK,ispec)
              Snp1   = factor_loc * epsilondev_trace_loc(INDEX_IJK)
              R_trace(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_trace(INDEX_IJK,ispec,i_sls) + &
                                betaval_loc * Sn + gammaval_loc * Snp1
            endif

            ! term in xx yy zz xy xz yz
            factor_loc = mustore(INDEX_IJK,ispec) * factor_common(i_sls,INDEX_IJK,ispec)

            ! term in xx
            Sn   = factor_loc * epsilondev_xx(INDEX_IJK,ispec)
            Snp1   = factor_loc * epsilondev_xx_loc(INDEX_IJK)
            R_xx(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xx(INDEX_IJK,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
            ! term in yy
            Sn   = factor_loc * epsilondev_yy(INDEX_IJK,ispec)
            Snp1   = factor_loc * epsilondev_yy_loc(INDEX_IJK)
            R_yy(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_yy(INDEX_IJK,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
            ! term in zz not computed since zero trace
            ! term in xy
            Sn   = factor_loc * epsilondev_xy(INDEX_IJK,ispec)
            Snp1   = factor_loc * epsilondev_xy_loc(INDEX_IJK)
            R_xy(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xy(INDEX_IJK,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
            ! term in xz
            Sn   = factor_loc * epsilondev_xz(INDEX_IJK,ispec)
            Snp1   = factor_loc * epsilondev_xz_loc(INDEX_IJK)
            R_xz(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_xz(INDEX_IJK,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
            ! term in yz
            Sn   = factor_loc * epsilondev_yz(INDEX_IJK,ispec)
            Snp1   = factor_loc * epsilondev_yz_loc(INDEX_IJK)
            R_yz(INDEX_IJK,ispec,i_sls) = alphaval_loc * R_yz(INDEX_IJK,ispec,i_sls) + &
                           betaval_loc * Sn + gammaval_loc * Snp1
          enddo   ! end of loop on memory variables
        endif

      endif  !  end of if attenuation

    ENDDO_LOOP_IJK

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then

        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          accel(1,iglob) = accel(1,iglob) - accel_elastic_CPML(1,INDEX_IJK)
          accel(2,iglob) = accel(2,iglob) - accel_elastic_CPML(2,INDEX_IJK)
          accel(3,iglob) = accel(3,iglob) - accel_elastic_CPML(3,INDEX_IJK)
        ENDDO_LOOP_IJK
      endif
    endif

    ! save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN) then
      if (FULL_ATTENUATION_SOLID) epsilondev_trace(:,:,:,ispec) = epsilondev_trace_loc(:,:,:)
      epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
    endif

! counts:
! + 0 FLOP
!
! + 0 BYTE

! counts:
! -----------------
! total of: 37625 FLOP per element
!
!           15204 BYTE DRAM accesses per block
!
! arithmetic intensity: 37625 FLOP / 15204 BYTES ~ 2.5 FLOP/BYTE
! -----------------


  enddo  ! spectral element loop

  contains

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
! note: the matrix-matrix multiplications are used for very small matrices ( default 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code
  subroutine mxm_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays e.g. (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,NGLLX),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(NGLLX,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! matrix-matrix multiplication wrapper

  if (NGLLX == 5) then
    call mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)
  else if (NGLLX == 6) then
    call mxm6_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)
  else if (NGLLX == 7) then
    call mxm7_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)
  endif

  end subroutine mxm_3comp_singleA

  !-------------

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

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

! 3 different arrays for x/y/z-components, 2-dimensional arrays (36,6)/(6,36), same B matrix for all 3 component arrays

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

! 3 different arrays for x/y/z-components, 2-dimensional arrays (49,7)/(7,49), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
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


!--------------------------------------------------------------------------------------------

  subroutine mxm_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays e.g. (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,NGLLX),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(NGLLX,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! matrix-matrix multiplication wrapper
  if (NGLLX == 5) then
    call mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)
  else if (NGLLX == 6) then
    call mxm6_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)
  else if (NGLLX == 7) then
    call mxm7_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)
  endif

  end subroutine mxm_3comp_singleB

  !-------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

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


!--------------------------------------------------------------------------------------------

  subroutine mxm_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 3-dimensional arrays e.g. (5,5,5), same B matrix for all 3 component arrays

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,NGLLX,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(NGLLX,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! matrix-matrix multiplication wrapper
  if (NGLLX == 5) then
    call mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)
  else if (NGLLX == 6) then
    call mxm6_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)
  else if (NGLLX == 7) then
    call mxm7_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)
  endif

  end subroutine mxm_3comp_3dmat_singleB

  !-------------

  subroutine mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

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

  end subroutine mxm5_3comp_3dmat_singleB

  !-------------

  subroutine mxm6_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

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

  end subroutine mxm6_3comp_3dmat_singleB

  !-------------

  subroutine mxm7_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

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

  end subroutine mxm7_3comp_3dmat_singleB

  end subroutine compute_forces_viscoelastic_Dev

