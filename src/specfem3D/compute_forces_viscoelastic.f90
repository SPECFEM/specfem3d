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

subroutine compute_forces_viscoelastic(iphase, &
                        NSPEC_AB,NGLOB_AB,displ,veloc,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,deltat,PML_CONDITIONS, &
                        one_minus_sum_beta,factor_common, &
                        one_minus_sum_beta_kappa,factor_common_kappa, &
                        alphaval,betaval,gammaval, &
                        NSPEC_ATTENUATION_AB, &
                        R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        NSPEC_ATTENUATION_AB_LDDRK,R_trace_lddrk, &
                        R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
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

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,ONE_THIRD,FOUR_THIRDS

  use fault_solver_dynamic, only: Kelvin_Voigt_eta

  use specfem_par, only: SAVE_MOHO_MESH,USE_LDDRK

  use pml_par, only: is_CPML,spec_to_CPML,accel_elastic_CPML, &
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
                     rmemory_displ_elastic,PML_displ_old,PML_displ_new

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! memory variables and standard linear solids for attenuation
  logical :: ATTENUATION
  logical :: COMPUTE_AND_STORE_STRAIN
  integer :: NSPEC_STRAIN_ONLY, NSPEC_ADJOINT
  integer :: NSPEC_ATTENUATION_AB
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
            R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: R_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: &
            epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: epsilondev_trace
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: epsilon_trace_over_3

! lddrk for update the memory variables
  integer :: NSPEC_ATTENUATION_AB_LDDRK
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS) :: &
            R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_LDDRK,N_SLS) :: R_trace_lddrk

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
  integer :: i_SLS,imodulo_N_SLS
  integer :: ispec,iglob,ispec_p,num_elements
  integer :: i,j,k,l

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

  ! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc, epsilondev_xx_loc, &
            epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) :: R_trace_val1,R_xx_val1,R_yy_val1
  real(kind=CUSTOM_REAL) :: R_trace_val2,R_xx_val2,R_yy_val2
  real(kind=CUSTOM_REAL) :: R_trace_val3,R_xx_val3,R_yy_val3
  real(kind=CUSTOM_REAL) :: templ

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
            newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL) duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att
  real(kind=CUSTOM_REAL) duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_att,tempx2_att,tempx3_att, &
            tempy1_att,tempy2_att,tempy3_att, &
            tempz1_att,tempz2_att,tempz3_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att,dummyy_loc_att,dummyz_loc_att

  real(kind=CUSTOM_REAL) :: eta

  ! local C-PML absorbing boundary conditions parameters
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
            tempx1_att_new,tempx2_att_new,tempx3_att_new, &
            tempy1_att_new,tempy2_att_new,tempy3_att_new, &
            tempz1_att_new,tempz2_att_new,tempz3_att_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: zero_array

  imodulo_N_SLS = mod(N_SLS,3)

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

! generate an array equal to zero
  zero_array(:,:,:) = 0._CUSTOM_REAL

  do ispec_p = 1,num_elements

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

    ! Kelvin Voigt damping: artificial viscosity around dynamic faults

    ! stores displacment values in local array
    !ZN ZN Do we need a stop statement somewhere in case of
    !ZN ZN SIMULATION_TYPE == 3 and allocated(Kelvin_Voigt_eta) = .true. ??
    if (allocated(Kelvin_Voigt_eta)) then

      eta = Kelvin_Voigt_eta(ispec)

      if (is_CPML(ispec) .and. eta /= 0._CUSTOM_REAL) stop 'you cannot put a fault inside a PML layer'

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

    if (is_CPML(ispec)) then
        ! In backward_simulation involved in SIMULATION_TYPE == 3,
        ! we only use the stored value on edge of PML interface.
        ! Thus no computation needs to be done in the PML region in this case.
        if (.not. backward_simulation) then
          ispec_CPML = spec_to_CPML(ispec)
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                 dummyx_loc_att(i,j,k) = PML_displ_old(1,i,j,k,ispec_CPML)
                 dummyy_loc_att(i,j,k) = PML_displ_old(2,i,j,k,ispec_CPML)
                 dummyz_loc_att(i,j,k) = PML_displ_old(3,i,j,k,ispec_CPML)
                 dummyx_loc_att_new(i,j,k) = PML_displ_new(1,i,j,k,ispec_CPML)
                 dummyy_loc_att_new(i,j,k) = PML_displ_new(2,i,j,k,ispec_CPML)
                 dummyz_loc_att_new(i,j,k) = PML_displ_new(3,i,j,k,ispec_CPML)
              enddo
            enddo
          enddo
        endif
    endif

    ! use first order Taylor expansion of displacement for local storage of stresses
    ! at this current time step, to fix attenuation in a consistent way
    if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN .and. .not. is_CPML(ispec)) then
          ! Keeping in mind, currently we implement a PML derived based on elastic wave equation
          ! for visco-elastic wave simulation. Thus it is not a PML anymore,
          ! because the solution of PML equation derived based on elastic wave equation is not perfectly mathced
          ! with solution of visco-elastic wave equation along the PML interface.
          ! you can seen PML in this case as a sponge layer.
          ! Based on limited numerical experiments, that PML implementation will work more or less OK anyway if Q > 70.
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
    endif

    !------------------------------------------------------------------------------
    !---------------------computation of strain in element-------------------------
    !------------------------------------------------------------------------------

    call compute_strain_in_element( &
                 tempx1,tempx2,tempx3,zero_array,zero_array,zero_array, &
                 tempy1,tempy2,tempy3,zero_array,zero_array,zero_array, &
                 tempz1,tempz2,tempz3,zero_array,zero_array,zero_array, &
                 dummyx_loc,dummyy_loc,dummyz_loc, &
                 hprime_xx,hprime_yy,hprime_zz)

    if (is_CPML(ispec)) then
        if (.not. backward_simulation) then
          call compute_strain_in_element( &
                       tempx1_att,tempx2_att,tempx3_att,zero_array,zero_array,zero_array, &
                       tempy1_att,tempy2_att,tempy3_att,zero_array,zero_array,zero_array, &
                       tempz1_att,tempz2_att,tempz3_att,zero_array,zero_array,zero_array, &
                       dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                       hprime_xx,hprime_yy,hprime_zz)

          call compute_strain_in_element( &
                       tempx1_att_new,tempx2_att_new,tempx3_att_new,zero_array,zero_array,zero_array, &
                       tempy1_att_new,tempy2_att_new,tempy3_att_new,zero_array,zero_array,zero_array, &
                       tempz1_att_new,tempz2_att_new,tempz3_att_new,zero_array,zero_array,zero_array, &
                       dummyx_loc_att_new,dummyy_loc_att_new,dummyz_loc_att_new, &
                       hprime_xx,hprime_yy,hprime_zz)
        endif
    endif

    if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN .and. .not. is_CPML(ispec)) then
        call compute_strain_in_element( &
                     tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                     tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                     tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                     dummyx_loc_att,dummyy_loc_att,dummyz_loc_att, &
                     hprime_xx,hprime_yy,hprime_zz)
    endif

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
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

          duxdxl = xixl * tempx1(i,j,k) + etaxl * tempx2(i,j,k) + gammaxl * tempx3(i,j,k)
          duxdyl = xiyl * tempx1(i,j,k) + etayl * tempx2(i,j,k) + gammayl * tempx3(i,j,k)
          duxdzl = xizl * tempx1(i,j,k) + etazl * tempx2(i,j,k) + gammazl * tempx3(i,j,k)

          duydxl = xixl * tempy1(i,j,k) + etaxl * tempy2(i,j,k) + gammaxl * tempy3(i,j,k)
          duydyl = xiyl * tempy1(i,j,k) + etayl * tempy2(i,j,k) + gammayl * tempy3(i,j,k)
          duydzl = xizl * tempy1(i,j,k) + etazl * tempy2(i,j,k) + gammazl * tempy3(i,j,k)

          duzdxl = xixl * tempz1(i,j,k) + etaxl * tempz2(i,j,k) + gammaxl * tempz3(i,j,k)
          duzdyl = xiyl * tempz1(i,j,k) + etayl * tempz2(i,j,k) + gammayl * tempz3(i,j,k)
          duzdzl = xizl * tempz1(i,j,k) + etazl * tempz2(i,j,k) + gammazl * tempz3(i,j,k)

          ! stores derivatives of ux, uy and uz with respect to x, y and z
          if (is_CPML(ispec)) then
              ! In backward_simulation involved in SIMULATION_TYPE == 3,
              ! we only use the stored value on edge of PML interface.
              ! Thus no computation needs to be done in the PML region in this case.
              if (.not. backward_simulation) then
                PML_dux_dxl(i,j,k) = duxdxl
                PML_dux_dyl(i,j,k) = duxdyl
                PML_dux_dzl(i,j,k) = duxdzl

                PML_duy_dxl(i,j,k) = duydxl
                PML_duy_dyl(i,j,k) = duydyl
                PML_duy_dzl(i,j,k) = duydzl

                PML_duz_dxl(i,j,k) = duzdxl
                PML_duz_dyl(i,j,k) = duzdyl
                PML_duz_dzl(i,j,k) = duzdzl

                PML_dux_dxl_old(i,j,k) = &
                    xixl * tempx1_att(i,j,k) + etaxl * tempx2_att(i,j,k) + gammaxl * tempx3_att(i,j,k)
                PML_dux_dyl_old(i,j,k) = &
                    xiyl * tempx1_att(i,j,k) + etayl * tempx2_att(i,j,k) + gammayl * tempx3_att(i,j,k)
                PML_dux_dzl_old(i,j,k) = &
                    xizl * tempx1_att(i,j,k) + etazl * tempx2_att(i,j,k) + gammazl * tempx3_att(i,j,k)

                PML_duy_dxl_old(i,j,k) = &
                    xixl * tempy1_att(i,j,k) + etaxl * tempy2_att(i,j,k) + gammaxl * tempy3_att(i,j,k)
                PML_duy_dyl_old(i,j,k) = &
                    xiyl * tempy1_att(i,j,k) + etayl * tempy2_att(i,j,k) + gammayl * tempy3_att(i,j,k)
                PML_duy_dzl_old(i,j,k) = &
                    xizl * tempy1_att(i,j,k) + etazl * tempy2_att(i,j,k) + gammazl * tempy3_att(i,j,k)

                PML_duz_dxl_old(i,j,k) = &
                    xixl * tempz1_att(i,j,k) + etaxl * tempz2_att(i,j,k) + gammaxl * tempz3_att(i,j,k)
                PML_duz_dyl_old(i,j,k) = &
                    xiyl * tempz1_att(i,j,k) + etayl * tempz2_att(i,j,k) + gammayl * tempz3_att(i,j,k)
                PML_duz_dzl_old(i,j,k) = &
                     xizl * tempz1_att(i,j,k) + etazl * tempz2_att(i,j,k) + gammazl * tempz3_att(i,j,k)

                PML_dux_dxl_new(i,j,k) = &
                    xixl * tempx1_att_new(i,j,k) + etaxl * tempx2_att_new(i,j,k) + gammaxl * tempx3_att_new(i,j,k)
                PML_dux_dyl_new(i,j,k) = &
                    xiyl * tempx1_att_new(i,j,k) + etayl * tempx2_att_new(i,j,k) + gammayl * tempx3_att_new(i,j,k)
                PML_dux_dzl_new(i,j,k) = &
                    xizl * tempx1_att_new(i,j,k) + etazl * tempx2_att_new(i,j,k) + gammazl * tempx3_att_new(i,j,k)

                PML_duy_dxl_new(i,j,k) = &
                    xixl * tempy1_att_new(i,j,k) + etaxl * tempy2_att_new(i,j,k) + gammaxl * tempy3_att_new(i,j,k)
                PML_duy_dyl_new(i,j,k) = &
                    xiyl * tempy1_att_new(i,j,k) + etayl * tempy2_att_new(i,j,k) + gammayl * tempy3_att_new(i,j,k)
                PML_duy_dzl_new(i,j,k) = &
                    xizl * tempy1_att_new(i,j,k) + etazl * tempy2_att_new(i,j,k) + gammazl * tempy3_att_new(i,j,k)

                PML_duz_dxl_new(i,j,k) = &
                    xixl * tempz1_att_new(i,j,k) + etaxl * tempz2_att_new(i,j,k) + gammaxl * tempz3_att_new(i,j,k)
                PML_duz_dyl_new(i,j,k) = &
                    xiyl * tempz1_att_new(i,j,k) + etayl * tempz2_att_new(i,j,k) + gammayl * tempz3_att_new(i,j,k)
                PML_duz_dzl_new(i,j,k) = &
                    xizl * tempz1_att_new(i,j,k) + etazl * tempz2_att_new(i,j,k) + gammazl * tempz3_att_new(i,j,k)
              endif

              if (COMPUTE_AND_STORE_STRAIN) then
                templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
                if (SIMULATION_TYPE == 3) epsilon_trace_over_3(i,j,k,ispec) = templ
                epsilondev_trace_loc(i,j,k) = 3._CUSTOM_REAL * templ
                epsilondev_xx_loc(i,j,k) = duxdxl - templ
                epsilondev_yy_loc(i,j,k) = duydyl - templ
                epsilondev_xy_loc(i,j,k) = 0.5_CUSTOM_REAL * (duxdyl + duydxl)
                epsilondev_xz_loc(i,j,k) = 0.5_CUSTOM_REAL * (duzdxl + duxdzl)
                epsilondev_yz_loc(i,j,k) = 0.5_CUSTOM_REAL * (duzdyl + duydzl)
              endif

          endif

          ! save strain on the Moho boundary
          if (SIMULATION_TYPE == 3 .and. SAVE_MOHO_MESH) then
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

          if (ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
            ! temporary variables used for fixing attenuation in a consistent way

              if (.not. is_CPML(ispec)) then ! 333333333333333
                duxdxl_att = xixl * tempx1_att(i,j,k) + etaxl * tempx2_att(i,j,k) + gammaxl * tempx3_att(i,j,k)
                duxdyl_att = xiyl * tempx1_att(i,j,k) + etayl * tempx2_att(i,j,k) + gammayl * tempx3_att(i,j,k)
                duxdzl_att = xizl * tempx1_att(i,j,k) + etazl * tempx2_att(i,j,k) + gammazl * tempx3_att(i,j,k)

                duydxl_att = xixl * tempy1_att(i,j,k) + etaxl * tempy2_att(i,j,k) + gammaxl * tempy3_att(i,j,k)
                duydyl_att = xiyl * tempy1_att(i,j,k) + etayl * tempy2_att(i,j,k) + gammayl * tempy3_att(i,j,k)
                duydzl_att = xizl * tempy1_att(i,j,k) + etazl * tempy2_att(i,j,k) + gammazl * tempy3_att(i,j,k)

                duzdxl_att = xixl * tempz1_att(i,j,k) + etaxl * tempz2_att(i,j,k) + gammaxl * tempz3_att(i,j,k)
                duzdyl_att = xiyl * tempz1_att(i,j,k) + etayl * tempz2_att(i,j,k) + gammayl * tempz3_att(i,j,k)
                duzdzl_att = xizl * tempz1_att(i,j,k) + etazl * tempz2_att(i,j,k) + gammazl * tempz3_att(i,j,k)

                ! precompute some sums to save CPU time
                duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
                duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
                duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

                ! compute deviatoric strain
                templ = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
                if (SIMULATION_TYPE == 3) epsilon_trace_over_3(i,j,k,ispec) = templ
                epsilondev_trace_loc(i,j,k) = 3._CUSTOM_REAL * templ
                epsilondev_xx_loc(i,j,k) = duxdxl_att - templ
                epsilondev_yy_loc(i,j,k) = duydyl_att - templ
                epsilondev_xy_loc(i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl_att
                epsilondev_xz_loc(i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl_att
                epsilondev_yz_loc(i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl_att
              endif ! 333333333333333

          else
            ! computes deviatoric strain attenuation and/or for kernel calculations
            if (COMPUTE_AND_STORE_STRAIN) then
              templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
              if (SIMULATION_TYPE == 3) epsilon_trace_over_3(i,j,k,ispec) = templ
              epsilondev_trace_loc(i,j,k) = 3._CUSTOM_REAL * templ
              epsilondev_xx_loc(i,j,k) = duxdxl - templ
              epsilondev_yy_loc(i,j,k) = duydyl - templ
              epsilondev_xy_loc(i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
              epsilondev_xz_loc(i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
              epsilondev_yz_loc(i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
            endif
          endif ! 111111111111111111

          kappal = kappastore(i,j,k,ispec)
          mul = mustore(i,j,k,ispec)

          ! use unrelaxed parameters if attenuation
          if (ATTENUATION .and. .not. is_CPML(ispec)) then
            mul  = mul * one_minus_sum_beta(i,j,k,ispec)
            kappal = kappal * one_minus_sum_beta_kappa(i,j,k,ispec)
          endif

          ! full anisotropic case, stress calculations
          if (ANISOTROPY) then
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

            sigma_xx = c11 * duxdxl + c16 * duxdyl_plus_duydxl + c12 * duydyl + &
                       c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * duzdzl
            sigma_yy = c12 * duxdxl + c26 * duxdyl_plus_duydxl + c22 * duydyl + &
                       c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * duzdzl
            sigma_zz = c13 * duxdxl + c36 * duxdyl_plus_duydxl + c23 * duydyl + &
                       c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * duzdzl
            sigma_xy = c16 * duxdxl + c66 * duxdyl_plus_duydxl + c26 * duydyl + &
                       c56 * duzdxl_plus_duxdzl + c46 * duzdyl_plus_duydzl + c36 * duzdzl
            sigma_xz = c15 * duxdxl + c56 * duxdyl_plus_duydxl + c25 * duydyl + &
                       c55 * duzdxl_plus_duxdzl + c45 * duzdyl_plus_duydzl + c35 * duzdzl
            sigma_yz = c14 * duxdxl + c46 * duxdyl_plus_duydxl + c24 * duydyl + &
                       c45 * duzdxl_plus_duxdzl + c44 * duzdyl_plus_duydzl + c34 * duzdzl

          else

            ! isotropic case
            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

            ! compute stress sigma
            sigma_xx = lambdalplus2mul * duxdxl + lambdal * duydyl_plus_duzdzl
            sigma_yy = lambdalplus2mul * duydyl + lambdal * duxdxl_plus_duzdzl
            sigma_zz = lambdalplus2mul * duzdzl + lambdal * duxdxl_plus_duydyl

            sigma_xy = mul * duxdyl_plus_duydxl
            sigma_xz = mul * duzdxl_plus_duxdzl
            sigma_yz = mul * duzdyl_plus_duydzl

          endif ! ANISOTROPY

          ! subtract memory variables if attenuation
          if (ATTENUATION .and. .not. is_CPML(ispec)) then
! old way 1
!             do i_sls = 1,N_SLS
!               R_xx_val = R_xx(i,j,k,ispec,i_sls)
!               R_yy_val = R_yy(i,j,k,ispec,i_sls)
!               sigma_xx = sigma_xx - R_xx_val
!               sigma_yy = sigma_yy - R_yy_val
!               sigma_zz = sigma_zz + R_xx_val + R_yy_val
!               sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
!               sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
!               sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
!             enddo

! new way 2
! note: this helps compilers to pipeline the code and make better use of the cache;
!       depending on compilers, it can further decrease the computation time by ~ 30%.
!       by default, N_SLS = 3, therefore we take steps of 3
                if (imodulo_N_SLS >= 1) then
                  do i_sls = 1,imodulo_N_SLS
                    R_trace_val1 = R_trace(i,j,k,ispec,i_sls)
                    R_xx_val1 = R_xx(i,j,k,ispec,i_sls)
                    R_yy_val1 = R_yy(i,j,k,ispec,i_sls)
                    sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
                    sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
                    sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
                    sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
                    sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
                    sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
                  enddo
                endif

                if (N_SLS >= imodulo_N_SLS+1) then
                  do i_sls = imodulo_N_SLS+1,N_SLS,3
                    R_trace_val1 = R_trace(i,j,k,ispec,i_sls)
                    R_xx_val1 = R_xx(i,j,k,ispec,i_sls)
                    R_yy_val1 = R_yy(i,j,k,ispec,i_sls)
                    sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
                    sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
                    sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
                    sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
                    sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
                    sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
                    R_trace_val2 = R_trace(i,j,k,ispec,i_sls+1)
                    R_xx_val2 = R_xx(i,j,k,ispec,i_sls+1)
                    R_yy_val2 = R_yy(i,j,k,ispec,i_sls+1)
                    sigma_xx = sigma_xx - R_xx_val2 - R_trace_val2
                    sigma_yy = sigma_yy - R_yy_val2 - R_trace_val2
                    sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2 - R_trace_val2
                    sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls+1)
                    sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls+1)
                    sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls+1)

                    R_trace_val3 = R_trace(i,j,k,ispec,i_sls+2)
                    R_xx_val3 = R_xx(i,j,k,ispec,i_sls+2)
                    R_yy_val3 = R_yy(i,j,k,ispec,i_sls+2)
                    sigma_xx = sigma_xx - R_xx_val3 - R_trace_val3
                    sigma_yy = sigma_yy - R_yy_val3 - R_trace_val3
                    sigma_zz = sigma_zz + R_xx_val3 + R_yy_val3 - R_trace_val3
                    sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls+2)
                    sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls+2)
                    sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls+2)
                  enddo
                endif

          endif

            if (.not. is_CPML(ispec)) then

              ! define symmetric components of sigma
              sigma_yx = sigma_xy
              sigma_zx = sigma_xz
              sigma_zy = sigma_yz

              ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
              tempx1(i,j,k) = jacobianl * (sigma_xx * xixl + sigma_yx * xiyl + sigma_zx * xizl) ! this goes to accel_x
              tempy1(i,j,k) = jacobianl * (sigma_xy * xixl + sigma_yy * xiyl + sigma_zy * xizl) ! this goes to accel_y
              tempz1(i,j,k) = jacobianl * (sigma_xz * xixl + sigma_yz * xiyl + sigma_zz * xizl) ! this goes to accel_z

              tempx2(i,j,k) = jacobianl * (sigma_xx * etaxl + sigma_yx * etayl + sigma_zx * etazl) ! this goes to accel_x
              tempy2(i,j,k) = jacobianl * (sigma_xy * etaxl + sigma_yy * etayl + sigma_zy * etazl) ! this goes to accel_y
              tempz2(i,j,k) = jacobianl * (sigma_xz * etaxl + sigma_yz * etayl + sigma_zz * etazl) ! this goes to accel_z

              tempx3(i,j,k) = jacobianl * (sigma_xx * gammaxl + sigma_yx * gammayl + sigma_zx * gammazl) ! this goes to accel_x
              tempy3(i,j,k) = jacobianl * (sigma_xy * gammaxl + sigma_yy * gammayl + sigma_zy * gammazl) ! this goes to accel_y
              tempz3(i,j,k) = jacobianl * (sigma_xz * gammaxl + sigma_yz * gammayl + sigma_zz * gammazl) ! this goes to accel_z
            endif

        enddo ! of the triple loop on i,j,k
      enddo
    enddo

    if (is_CPML(ispec) .and. .not. backward_simulation) then
        ! In backward_simulation involved in SIMULATION_TYPE == 3,
        ! we only use the stored value on edge of PML interface.
        ! Thus no computation needs to be done in the PML region in this case.
          ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
          ispec_CPML = spec_to_CPML(ispec)
          call pml_compute_memory_variables_elastic(ispec,ispec_CPML,tempx1,tempy1,tempz1, &
                   tempx2,tempy2,tempz2,tempx3,tempy3,tempz3, &
                   rmemory_dux_dxl_x, rmemory_duy_dyl_x, rmemory_duz_dzl_x, &
                   rmemory_dux_dyl_x, rmemory_dux_dzl_x, rmemory_duz_dxl_x, rmemory_duy_dxl_x, &
                   rmemory_dux_dxl_y, rmemory_duz_dzl_y, rmemory_duy_dyl_y, &
                   rmemory_duy_dxl_y, rmemory_duy_dzl_y, rmemory_duz_dyl_y, rmemory_dux_dyl_y, &
                   rmemory_dux_dxl_z, rmemory_duy_dyl_z, rmemory_duz_dzl_z, &
                   rmemory_duz_dxl_z, rmemory_duz_dyl_z, rmemory_duy_dzl_z, rmemory_dux_dzl_z)

          ! calculates contribution from each C-PML element to update acceleration
          call pml_compute_accel_contribution_elastic(ispec,ispec_CPML,displ,veloc,rmemory_displ_elastic)
    endif

    ! second double-loop over GLL to compute all the terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          newtempx1(i,j,k) = 0._CUSTOM_REAL
          newtempy1(i,j,k) = 0._CUSTOM_REAL
          newtempz1(i,j,k) = 0._CUSTOM_REAL

          newtempx2(i,j,k) = 0._CUSTOM_REAL
          newtempy2(i,j,k) = 0._CUSTOM_REAL
          newtempz2(i,j,k) = 0._CUSTOM_REAL

          newtempx3(i,j,k) = 0._CUSTOM_REAL
          newtempy3(i,j,k) = 0._CUSTOM_REAL
          newtempz3(i,j,k) = 0._CUSTOM_REAL

          ! we can merge these loops because NGLLX = NGLLY = NGLLZ
          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            newtempx1(i,j,k) = newtempx1(i,j,k) + tempx1(l,j,k) * fac1
            newtempy1(i,j,k) = newtempy1(i,j,k) + tempy1(l,j,k) * fac1
            newtempz1(i,j,k) = newtempz1(i,j,k) + tempz1(l,j,k) * fac1

            fac2 = hprimewgll_yy(l,j)
            newtempx2(i,j,k) = newtempx2(i,j,k) + tempx2(i,l,k) * fac2
            newtempy2(i,j,k) = newtempy2(i,j,k) + tempy2(i,l,k) * fac2
            newtempz2(i,j,k) = newtempz2(i,j,k) + tempz2(i,l,k) * fac2

            fac3 = hprimewgll_zz(l,k)
            newtempx3(i,j,k) = newtempx3(i,j,k) + tempx3(i,j,l) * fac3
            newtempy3(i,j,k) = newtempy3(i,j,k) + tempy3(i,j,l) * fac3
            newtempz3(i,j,k) = newtempz3(i,j,k) + tempz3(i,j,l) * fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions from each element to the global mesh using indirect addressing
          iglob = ibool(i,j,k,ispec)
          accel(1,iglob) = accel(1,iglob) - (fac1 * newtempx1(i,j,k) + fac2 * newtempx2(i,j,k) + fac3 * newtempx3(i,j,k))
          accel(2,iglob) = accel(2,iglob) - (fac1 * newtempy1(i,j,k) + fac2 * newtempy2(i,j,k) + fac3 * newtempy3(i,j,k))
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

      if (is_CPML(ispec) .and. .not. backward_simulation) then
        ! In backward_simulation involved in SIMULATION_TYPE == 3,
        ! we only use the stored value on edge of PML interface.
        ! Thus no computation needs to be done in the PML region in this case.
          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)
                accel(1,iglob) = accel(1,iglob) - accel_elastic_CPML(1,i,j,k)
                accel(2,iglob) = accel(2,iglob) - accel_elastic_CPML(2,i,j,k)
                accel(3,iglob) = accel(3,iglob) - accel_elastic_CPML(3,i,j,k)
              enddo
            enddo
          enddo
      endif

    ! save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN) then
      epsilondev_trace(:,:,:,ispec) = epsilondev_trace_loc(:,:,:)
      epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
    endif

  enddo  ! spectral element loop

end subroutine compute_forces_viscoelastic

!
!---------------
!

! put the code used for computation of strain in element in a subroutine

subroutine compute_strain_in_element(tempx1_att,tempx2_att,tempx3_att,tempx1,tempx2,tempx3, &
                                            tempy1_att,tempy2_att,tempy3_att,tempy1,tempy2,tempy3, &
                                            tempz1_att,tempz2_att,tempz3_att,tempz1,tempz2,tempz3, &
                                            dummyx_loc,dummyy_loc,dummyz_loc,hprime_xx,hprime_yy,hprime_zz)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_att,tempx2_att,tempx3_att, &
                                                          tempy1_att,tempy2_att,tempy3_att, &
                                                          tempz1_att,tempz2_att,tempz3_att

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3, &
                                                          tempy1,tempy2,tempy3, &
                                                          tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  ! local variables
  integer :: i,j,k,l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  tempx1_att(:,:,:) = tempx1(:,:,:)
  tempx2_att(:,:,:) = tempx2(:,:,:)
  tempx3_att(:,:,:) = tempx3(:,:,:)

  tempy1_att(:,:,:) = tempy1(:,:,:)
  tempy2_att(:,:,:) = tempy2(:,:,:)
  tempy3_att(:,:,:) = tempy3(:,:,:)

  tempz1_att(:,:,:) = tempz1(:,:,:)
  tempz2_att(:,:,:) = tempz2(:,:,:)
  tempz3_att(:,:,:) = tempz3(:,:,:)

  ! use first order Taylor expansion of displacement for local storage of stresses
  ! at this current time step, to fix attenuation in a consistent way
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

        ! we can merge these loops because NGLLX = NGLLY = NGLLZ
        do l=1,NGLLX
          hp1 = hprime_xx(i,l)
          tempx1_att(i,j,k) = tempx1_att(i,j,k) + dummyx_loc(l,j,k) * hp1
          tempy1_att(i,j,k) = tempy1_att(i,j,k) + dummyy_loc(l,j,k) * hp1
          tempz1_att(i,j,k) = tempz1_att(i,j,k) + dummyz_loc(l,j,k) * hp1

          hp2 = hprime_yy(j,l)
          tempx2_att(i,j,k) = tempx2_att(i,j,k) + dummyx_loc(i,l,k) * hp2
          tempy2_att(i,j,k) = tempy2_att(i,j,k) + dummyy_loc(i,l,k) * hp2
          tempz2_att(i,j,k) = tempz2_att(i,j,k) + dummyz_loc(i,l,k) * hp2

          hp3 = hprime_zz(k,l)
          tempx3_att(i,j,k) = tempx3_att(i,j,k) + dummyx_loc(i,j,l) * hp3
          tempy3_att(i,j,k) = tempy3_att(i,j,k) + dummyy_loc(i,j,l) * hp3
          tempz3_att(i,j,k) = tempz3_att(i,j,k) + dummyz_loc(i,j,l) * hp3
        enddo

      enddo
    enddo
  enddo

end subroutine compute_strain_in_element

