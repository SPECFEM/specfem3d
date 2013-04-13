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

subroutine compute_forces_viscoelastic_noDev(iphase, &
                        NSPEC_AB,NGLOB_AB,displ,veloc,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,deltat,PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                        one_minus_sum_beta,factor_common, &
                        one_minus_sum_beta_kappa,factor_common_kappa, &
                        alphaval,betaval,gammaval, &
                        NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa, &
                        R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,&
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
                        phase_ispec_inner_elastic,ispec_is_elastic)

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,SAVE_MOHO_MESH,ONE_THIRD,FOUR_THIRDS,IOUT
  use pml_par
  use fault_solver_dynamic, only : Kelvin_Voigt_eta
  use specfem_par, only : FULL_ATTENUATION_SOLID

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
  integer :: NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB) :: factor_common
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa) :: one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa) :: factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB,N_SLS) :: &
       R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa,N_SLS) :: R_trace

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY) :: &
       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATTENUATION_AB_kappa) :: epsilondev_trace
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
  logical :: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(nspec2D_xmin) :: ibelm_xmin
  integer, dimension(nspec2D_xmax) :: ibelm_xmax
  integer, dimension(nspec2D_ymin) :: ibelm_ymin
  integer, dimension(nspec2D_ymax) :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top

   logical, dimension(NSPEC_AB) :: ispec_is_elastic

! local parameters
  integer :: i_SLS,imodulo_N_SLS
  integer :: ispec,ispec2D,iglob,ispec_p,num_elements
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local attenuation parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_loc, epsilondev_xx_loc, &
       epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
  real(kind=CUSTOM_REAL) :: R_trace_val1,R_xx_val1,R_yy_val1
  real(kind=CUSTOM_REAL) :: R_trace_val2,R_xx_val2,R_yy_val2
  real(kind=CUSTOM_REAL) :: R_trace_val3,R_xx_val3,R_yy_val3
  real(kind=CUSTOM_REAL) :: factor_loc,alphaval_loc,betaval_loc,gammaval_loc,Sn,Snp1
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
    tempx1_att,tempx2_att,tempx3_att,tempy1_att,tempy2_att,tempy3_att,tempz1_att,tempz2_att,tempz3_att
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc_att,dummyy_loc_att,dummyz_loc_att

  real(kind=CUSTOM_REAL) :: eta

  ! local C-PML absorbing boundary conditions parameters
  integer :: ispec_CPML

  imodulo_N_SLS = mod(N_SLS,3)

  ! choses inner/outer elements
  if( iphase == 1 ) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  do ispec_p = 1,num_elements

        ! returns element id from stored element list
        ispec = phase_ispec_inner_elastic(ispec_p,iphase)

        ! adjoint simulations: moho kernel
        if( SIMULATION_TYPE == 3 .and. SAVE_MOHO_MESH ) then
          if (is_moho_top(ispec)) then
            ispec2D_moho_top = ispec2D_moho_top + 1
          else if (is_moho_bot(ispec)) then
            ispec2D_moho_bot = ispec2D_moho_bot + 1
          endif
        endif ! adjoint

       ! Kelvin Voigt damping: artificial viscosity around dynamic faults

        ! stores displacment values in local array
        if (allocated(Kelvin_Voigt_eta)) then
          eta = Kelvin_Voigt_eta(ispec)
          do k=1,NGLLZ
            do j=1,NGLLY
              do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                dummyx_loc(i,j,k) = displ(1,iglob) + eta*veloc(1,iglob)
                dummyy_loc(i,j,k) = displ(2,iglob) + eta*veloc(2,iglob)
                dummyz_loc(i,j,k) = displ(3,iglob) + eta*veloc(3,iglob)
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

        ! use first order Taylor expansion of displacement for local storage of stresses
        ! at this current time step, to fix attenuation in a consistent way
        if(ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
           do k=1,NGLLZ
              do j=1,NGLLY
                 do i=1,NGLLX
                    iglob = ibool(i,j,k,ispec)
                    dummyx_loc_att(i,j,k) = deltat*veloc(1,iglob)
                    dummyy_loc_att(i,j,k) = deltat*veloc(2,iglob)
                    dummyz_loc_att(i,j,k) = deltat*veloc(3,iglob)
                 enddo
              enddo
           enddo
        else if(PML_CONDITIONS) then
           ! do not merge this second line with the first using an ".and." statement
           ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
           if(CPML_mask_ibool(ispec)) then
              do k=1,NGLLZ
                 do j=1,NGLLY
                    do i=1,NGLLX
                       iglob = ibool(i,j,k,ispec)
                       dummyx_loc_att(i,j,k) = deltat*veloc(1,iglob)
                       dummyy_loc_att(i,j,k) = deltat*veloc(2,iglob)
                       dummyz_loc_att(i,j,k) = deltat*veloc(3,iglob)
                    enddo
                 enddo
              enddo
           endif
        endif

     do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

          tempx1(i,j,k) = 0._CUSTOM_REAL
          tempx2(i,j,k) = 0._CUSTOM_REAL
          tempx3(i,j,k) = 0._CUSTOM_REAL

          tempy1(i,j,k) = 0._CUSTOM_REAL
          tempy2(i,j,k) = 0._CUSTOM_REAL
          tempy3(i,j,k) = 0._CUSTOM_REAL

          tempz1(i,j,k) = 0._CUSTOM_REAL
          tempz2(i,j,k) = 0._CUSTOM_REAL
          tempz3(i,j,k) = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1
            tempy1(i,j,k) = tempy1(i,j,k) + dummyy_loc(l,j,k)*hp1
            tempz1(i,j,k) = tempz1(i,j,k) + dummyz_loc(l,j,k)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2
            tempy2(i,j,k) = tempy2(i,j,k) + dummyy_loc(i,l,k)*hp2
            tempz2(i,j,k) = tempz2(i,j,k) + dummyz_loc(i,l,k)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
            tempy3(i,j,k) = tempy3(i,j,k) + dummyy_loc(i,j,l)*hp3
            tempz3(i,j,k) = tempz3(i,j,k) + dummyz_loc(i,j,l)*hp3
          enddo

          if(ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
             tempx1_att(i,j,k) = tempx1(i,j,k)
             tempx2_att(i,j,k) = tempx2(i,j,k)
             tempx3_att(i,j,k) = tempx3(i,j,k)

             tempy1_att(i,j,k) = tempy1(i,j,k)
             tempy2_att(i,j,k) = tempy2(i,j,k)
             tempy3_att(i,j,k) = tempy3(i,j,k)

             tempz1_att(i,j,k) = tempz1(i,j,k)
             tempz2_att(i,j,k) = tempz2(i,j,k)
             tempz3_att(i,j,k) = tempz3(i,j,k)

             ! use first order Taylor expansion of displacement for local storage of stresses
             ! at this current time step, to fix attenuation in a consistent way
             do l=1,NGLLX
                hp1 = hprime_xx(i,l)
                tempx1_att(i,j,k) = tempx1_att(i,j,k) + dummyx_loc_att(l,j,k)*hp1
                tempy1_att(i,j,k) = tempy1_att(i,j,k) + dummyy_loc_att(l,j,k)*hp1
                tempz1_att(i,j,k) = tempz1_att(i,j,k) + dummyz_loc_att(l,j,k)*hp1

!!! can merge these loops because NGLLX = NGLLY = NGLLZ
                hp2 = hprime_yy(j,l)
                tempx2_att(i,j,k) = tempx2_att(i,j,k) + dummyx_loc_att(i,l,k)*hp2
                tempy2_att(i,j,k) = tempy2_att(i,j,k) + dummyy_loc_att(i,l,k)*hp2
                tempz2_att(i,j,k) = tempz2_att(i,j,k) + dummyz_loc_att(i,l,k)*hp2

!!! can merge these loops because NGLLX = NGLLY = NGLLZ
                hp3 = hprime_zz(k,l)
                tempx3_att(i,j,k) = tempx3_att(i,j,k) + dummyx_loc_att(i,j,l)*hp3
                tempy3_att(i,j,k) = tempy3_att(i,j,k) + dummyy_loc_att(i,j,l)*hp3
                tempz3_att(i,j,k) = tempz3_att(i,j,k) + dummyz_loc_att(i,j,l)*hp3
             enddo

          else if(PML_CONDITIONS) then
             ! do not merge this second line with the first using an ".and." statement
             ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
             if(CPML_mask_ibool(ispec)) then
                tempx1_att(i,j,k) = tempx1(i,j,k)
                tempx2_att(i,j,k) = tempx2(i,j,k)
                tempx3_att(i,j,k) = tempx3(i,j,k)

                tempy1_att(i,j,k) = tempy1(i,j,k)
                tempy2_att(i,j,k) = tempy2(i,j,k)
                tempy3_att(i,j,k) = tempy3(i,j,k)

                tempz1_att(i,j,k) = tempz1(i,j,k)
                tempz2_att(i,j,k) = tempz2(i,j,k)
                tempz3_att(i,j,k) = tempz3(i,j,k)

                ! use first order Taylor expansion of displacement for local storage of stresses
                ! at this current time step, to fix attenuation in a consistent way
                do l=1,NGLLX
                   hp1 = hprime_xx(i,l)
                   tempx1_att(i,j,k) = tempx1_att(i,j,k) + dummyx_loc_att(l,j,k)*hp1
                   tempy1_att(i,j,k) = tempy1_att(i,j,k) + dummyy_loc_att(l,j,k)*hp1
                   tempz1_att(i,j,k) = tempz1_att(i,j,k) + dummyz_loc_att(l,j,k)*hp1

!!! can merge these loops because NGLLX = NGLLY = NGLLZ
                   hp2 = hprime_yy(j,l)
                   tempx2_att(i,j,k) = tempx2_att(i,j,k) + dummyx_loc_att(i,l,k)*hp2
                   tempy2_att(i,j,k) = tempy2_att(i,j,k) + dummyy_loc_att(i,l,k)*hp2
                   tempz2_att(i,j,k) = tempz2_att(i,j,k) + dummyz_loc_att(i,l,k)*hp2

!!! can merge these loops because NGLLX = NGLLY = NGLLZ
                   hp3 = hprime_zz(k,l)
                   tempx3_att(i,j,k) = tempx3_att(i,j,k) + dummyx_loc_att(i,j,l)*hp3
                   tempy3_att(i,j,k) = tempy3_att(i,j,k) + dummyy_loc_att(i,j,l)*hp3
                   tempz3_att(i,j,k) = tempz3_att(i,j,k) + dummyz_loc_att(i,j,l)*hp3
                enddo
             endif
          endif

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

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

              ! stores derivatives of ux, uy and uz with respect to x, y and z
              if (PML_CONDITIONS) then
                 ! do not merge this second line with the first using an ".and." statement
                 ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
                 if(CPML_mask_ibool(ispec)) then
                    ispec_CPML = spec_to_CPML(ispec)

                    PML_dux_dxl(i,j,k,ispec_CPML) = duxdxl
                    PML_dux_dyl(i,j,k,ispec_CPML) = duxdyl
                    PML_dux_dzl(i,j,k,ispec_CPML) = duxdzl

                    PML_duy_dxl(i,j,k,ispec_CPML) = duydxl
                    PML_duy_dyl(i,j,k,ispec_CPML) = duydyl
                    PML_duy_dzl(i,j,k,ispec_CPML) = duydzl

                    PML_duz_dxl(i,j,k,ispec_CPML) = duzdxl
                    PML_duz_dyl(i,j,k,ispec_CPML) = duzdyl
                    PML_duz_dzl(i,j,k,ispec_CPML) = duzdzl
                 endif
              endif

              ! save strain on the Moho boundary
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

              if(ATTENUATION .and. COMPUTE_AND_STORE_STRAIN) then
                 ! temporary variables used for fixing attenuation in a consistent way
                 duxdxl_att = xixl*tempx1_att(i,j,k) + etaxl*tempx2_att(i,j,k) + gammaxl*tempx3_att(i,j,k)
                 duxdyl_att = xiyl*tempx1_att(i,j,k) + etayl*tempx2_att(i,j,k) + gammayl*tempx3_att(i,j,k)
                 duxdzl_att = xizl*tempx1_att(i,j,k) + etazl*tempx2_att(i,j,k) + gammazl*tempx3_att(i,j,k)

                 duydxl_att = xixl*tempy1_att(i,j,k) + etaxl*tempy2_att(i,j,k) + gammaxl*tempy3_att(i,j,k)
                 duydyl_att = xiyl*tempy1_att(i,j,k) + etayl*tempy2_att(i,j,k) + gammayl*tempy3_att(i,j,k)
                 duydzl_att = xizl*tempy1_att(i,j,k) + etazl*tempy2_att(i,j,k) + gammazl*tempy3_att(i,j,k)

                 duzdxl_att = xixl*tempz1_att(i,j,k) + etaxl*tempz2_att(i,j,k) + gammaxl*tempz3_att(i,j,k)
                 duzdyl_att = xiyl*tempz1_att(i,j,k) + etayl*tempz2_att(i,j,k) + gammayl*tempz3_att(i,j,k)
                 duzdzl_att = xizl*tempz1_att(i,j,k) + etazl*tempz2_att(i,j,k) + gammazl*tempz3_att(i,j,k)

                 ! precompute some sums to save CPU time
                 duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
                 duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
                 duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

                 ! compute deviatoric strain
                 templ = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
                 if( SIMULATION_TYPE == 3 ) epsilon_trace_over_3(i,j,k,ispec) = templ
                 if(FULL_ATTENUATION_SOLID) epsilondev_trace_loc(i,j,k) = 3.0 * templ
                 epsilondev_xx_loc(i,j,k) = duxdxl_att - templ
                 epsilondev_yy_loc(i,j,k) = duydyl_att - templ
                 epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl_att
                 epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl_att
                 epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl_att

              else if(PML_CONDITIONS) then
                    ! do not merge this second line with the first using an ".and." statement
                    ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
                    if(CPML_mask_ibool(ispec)) then
                       PML_dux_dxl_new(i,j,k,ispec_CPML) = &
                            xixl*tempx1_att(i,j,k) + etaxl*tempx2_att(i,j,k) + gammaxl*tempx3_att(i,j,k)
                       PML_dux_dyl_new(i,j,k,ispec_CPML) = &
                            xiyl*tempx1_att(i,j,k) + etayl*tempx2_att(i,j,k) + gammayl*tempx3_att(i,j,k)
                       PML_dux_dzl_new(i,j,k,ispec_CPML) = &
                            xizl*tempx1_att(i,j,k) + etazl*tempx2_att(i,j,k) + gammazl*tempx3_att(i,j,k)

                       PML_duy_dxl_new(i,j,k,ispec_CPML) = &
                            xixl*tempy1_att(i,j,k) + etaxl*tempy2_att(i,j,k) + gammaxl*tempy3_att(i,j,k)
                       PML_duy_dyl_new(i,j,k,ispec_CPML) = &
                            xiyl*tempy1_att(i,j,k) + etayl*tempy2_att(i,j,k) + gammayl*tempy3_att(i,j,k)
                       PML_duy_dzl_new(i,j,k,ispec_CPML) = &
                            xizl*tempy1_att(i,j,k) + etazl*tempy2_att(i,j,k) + gammazl*tempy3_att(i,j,k)

                       PML_duz_dxl_new(i,j,k,ispec_CPML) = &
                            xixl*tempz1_att(i,j,k) + etaxl*tempz2_att(i,j,k) + gammaxl*tempz3_att(i,j,k)
                       PML_duz_dyl_new(i,j,k,ispec_CPML) = &
                            xiyl*tempz1_att(i,j,k) + etayl*tempz2_att(i,j,k) + gammayl*tempz3_att(i,j,k)
                       PML_duz_dzl_new(i,j,k,ispec_CPML) = &
                            xizl*tempz1_att(i,j,k) + etazl*tempz2_att(i,j,k) + gammazl*tempz3_att(i,j,k)
                    endif

              else
                 ! computes deviatoric strain attenuation and/or for kernel calculations
                 if (COMPUTE_AND_STORE_STRAIN) then
                    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
                    if( SIMULATION_TYPE == 3 ) epsilon_trace_over_3(i,j,k,ispec) = templ
                    if(FULL_ATTENUATION_SOLID) epsilondev_trace_loc(i,j,k) = 3.0 * templ
                    epsilondev_xx_loc(i,j,k) = duxdxl - templ
                    epsilondev_yy_loc(i,j,k) = duydyl - templ
                    epsilondev_xy_loc(i,j,k) = 0.5 * duxdyl_plus_duydxl
                    epsilondev_xz_loc(i,j,k) = 0.5 * duzdxl_plus_duxdzl
                    epsilondev_yz_loc(i,j,k) = 0.5 * duzdyl_plus_duydzl
                 endif
              endif

              kappal = kappastore(i,j,k,ispec)
              mul = mustore(i,j,k,ispec)

              ! attenuation
              if(ATTENUATION) then
                ! use unrelaxed parameters if attenuation
                mul  = mul * one_minus_sum_beta(i,j,k,ispec)
                if(FULL_ATTENUATION_SOLID) kappal  = kappal * one_minus_sum_beta_kappa(i,j,k,ispec)
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
! way 1
!                do i_sls = 1,N_SLS
!                  R_xx_val = R_xx(i,j,k,ispec,i_sls)
!                  R_yy_val = R_yy(i,j,k,ispec,i_sls)
!                  sigma_xx = sigma_xx - R_xx_val
!                  sigma_yy = sigma_yy - R_yy_val
!                  sigma_zz = sigma_zz + R_xx_val + R_yy_val
!                  sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
!                  sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
!                  sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
!                enddo

! way 2
! note: this should help compilers to pipeline the code and make better use of the cache;
!          depending on compilers, it can further decrease the computation time by ~ 30%.
!          by default, N_SLS = 3, therefore we take steps of 3
              if(imodulo_N_SLS >= 1) then
                do i_sls = 1,imodulo_N_SLS
                  if(FULL_ATTENUATION_SOLID) then
                    R_trace_val1 = R_trace(i,j,k,ispec,i_sls)
                  else
                    R_trace_val1 = 0.
                  endif
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

              if(N_SLS >= imodulo_N_SLS+1) then
                do i_sls = imodulo_N_SLS+1,N_SLS,3
                  if(FULL_ATTENUATION_SOLID) then
                    R_trace_val1 = R_trace(i,j,k,ispec,i_sls)
                  else
                    R_trace_val1 = 0.
                  endif
                  R_xx_val1 = R_xx(i,j,k,ispec,i_sls)
                  R_yy_val1 = R_yy(i,j,k,ispec,i_sls)
                  sigma_xx = sigma_xx - R_xx_val1 - R_trace_val1
                  sigma_yy = sigma_yy - R_yy_val1 - R_trace_val1
                  sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1 - R_trace_val1
                  sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls)
                  sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls)
                  sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls)
                  if(FULL_ATTENUATION_SOLID) then
                    R_trace_val2 = R_trace(i,j,k,ispec,i_sls+1)
                  else
                    R_trace_val2 = 0.
                  endif
                  R_xx_val2 = R_xx(i,j,k,ispec,i_sls+1)
                  R_yy_val2 = R_yy(i,j,k,ispec,i_sls+1)
                  sigma_xx = sigma_xx - R_xx_val2 - R_trace_val2
                  sigma_yy = sigma_yy - R_yy_val2 - R_trace_val2
                  sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2 - R_trace_val2
                  sigma_xy = sigma_xy - R_xy(i,j,k,ispec,i_sls+1)
                  sigma_xz = sigma_xz - R_xz(i,j,k,ispec,i_sls+1)
                  sigma_yz = sigma_yz - R_yz(i,j,k,ispec,i_sls+1)

                  if(FULL_ATTENUATION_SOLID) then
                    R_trace_val3 = R_trace(i,j,k,ispec,i_sls+2)
                  else
                    R_trace_val3 = 0.
                  endif
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

!! DK DK to Jo, to debug CPML, 22 March 2013:

!! DK DK comment from DK DK, 22 March 2013, for Jo and Zhinan to debug CPML:
!! DK DK are you sure about this "if" statement below? because I am surprised to see
!! DK DK that when PML_CONDITIONS is on then you do not compute the tempx, tempy, tempz arrays
!! DK DK (even in non-PML elements!!), even though such arrays are needed below;
!! DK DK shouldn't there be at least a "if (is_CPML(ispec))" test as well here, or something like that?
              if (PML_CONDITIONS) then
                 ! do not merge this second line with the first using an ".and." statement
                 ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
                 if(.not.CPML_mask_ibool(ispec)) then
                    ! define symmetric components of sigma
                    sigma_yx = sigma_xy
                    sigma_zx = sigma_xz
                    sigma_zy = sigma_yz

                    ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
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
              else
                 ! define symmetric components of sigma
                 sigma_yx = sigma_xy
                 sigma_zx = sigma_xz
                 sigma_zy = sigma_yz

                 ! form dot product with test vector, non-symmetric form (which is useful in the case of PML)
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

    if (PML_CONDITIONS) then
       ! do not merge this second line with the first using an ".and." statement
       ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
       if(CPML_mask_ibool(ispec)) then
          ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
          call pml_compute_memory_variables(ispec,ispec_CPML,deltat,tempx1,tempy1,tempz1,tempx2,tempy2,tempz2, &
               tempx3,tempy3,tempz3,NSPEC_AB,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian)

          ! calculates contribution from each C-PML element to update acceleration
          call pml_compute_accel_contribution(ispec,ispec_CPML,deltat,nspec_AB,jacobian,accel_elastic_CPML)
       endif
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

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            newtempx1(i,j,k) = newtempx1(i,j,k) + tempx1(l,j,k)*fac1
            newtempy1(i,j,k) = newtempy1(i,j,k) + tempy1(l,j,k)*fac1
            newtempz1(i,j,k) = newtempz1(i,j,k) + tempz1(l,j,k)*fac1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac2 = hprimewgll_yy(l,j)
            newtempx2(i,j,k) = newtempx2(i,j,k) + tempx2(i,l,k)*fac2
            newtempy2(i,j,k) = newtempy2(i,j,k) + tempy2(i,l,k)*fac2
            newtempz2(i,j,k) = newtempz2(i,j,k) + tempz2(i,l,k)*fac2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            fac3 = hprimewgll_zz(l,k)
            newtempx3(i,j,k) = newtempx3(i,j,k) + tempx3(i,j,l)*fac3
            newtempy3(i,j,k) = newtempy3(i,j,k) + tempy3(i,j,l)*fac3
            newtempz3(i,j,k) = newtempz3(i,j,k) + tempz3(i,j,l)*fac3
          enddo

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

          ! updates acceleration with contribution from each C-PML element
          if (PML_CONDITIONS) then
             ! do not merge this second line with the first using an ".and." statement
             ! because array CPML_mask_ibool() is unallocated when PML_CONDITIONS is false
             if(CPML_mask_ibool(ispec)) then
                accel(1,iglob) = accel(1,iglob) - accel_elastic_CPML(1,i,j,k,ispec_CPML)
                accel(2,iglob) = accel(2,iglob) - accel_elastic_CPML(2,i,j,k,ispec_CPML)
                accel(3,iglob) = accel(3,iglob) - accel_elastic_CPML(3,i,j,k,ispec_CPML)
             endif
          endif

          !  update memory variables based upon the Runge-Kutta scheme
          if(ATTENUATION) then

             ! use Runge-Kutta scheme to march in time
             do i_sls = 1,N_SLS

                alphaval_loc = alphaval(i_sls)
                betaval_loc = betaval(i_sls)
                gammaval_loc = gammaval(i_sls)

                if(FULL_ATTENUATION_SOLID) then
                   ! term in trace
                   factor_loc = kappastore(i,j,k,ispec) * factor_common_kappa(i_sls,i,j,k,ispec)

                   Sn   = factor_loc * epsilondev_trace(i,j,k,ispec)
                   Snp1   = factor_loc * epsilondev_trace_loc(i,j,k)
                   R_trace(i,j,k,ispec,i_sls) = alphaval_loc * R_trace(i,j,k,ispec,i_sls) + &
                        betaval_loc * Sn + gammaval_loc * Snp1
                endif

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
             enddo   ! end of loop on memory variables

          endif  !  end of if attenuation

        enddo
      enddo
    enddo

        ! save deviatoric strain for Runge-Kutta scheme
        if ( COMPUTE_AND_STORE_STRAIN ) then
          if(FULL_ATTENUATION_SOLID) epsilondev_trace(:,:,:,ispec) = epsilondev_trace_loc(:,:,:)
          epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
          epsilondev_yy(:,:,:,ispec) = epsilondev_yy_loc(:,:,:)
          epsilondev_xy(:,:,:,ispec) = epsilondev_xy_loc(:,:,:)
          epsilondev_xz(:,:,:,ispec) = epsilondev_xz_loc(:,:,:)
          epsilondev_yz(:,:,:,ispec) = epsilondev_yz_loc(:,:,:)
        endif

  enddo  ! spectral element loop

!! DK DK to Jo, to debug CPML, 22 March 2013:

!! DK DK I think that there is an error in the loops below, you should also check if ispec is a CPML element,
!! DK DK and also if ispec is an elastic or viscoelastic element (and NOT for instance an acoustic element)
!! DK DK
!! DK DK thus test something like:   if (is_CPML(ispec) .and. elastic(ispec)) then
!! DK DK or something like that

  ! C-PML boundary
  if( PML_CONDITIONS ) then
     ! xmin
     do ispec2D=1,nspec2D_xmin
        ispec = ibelm_xmin(ispec2D)

        if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
           i = 1

           do k=1,NGLLZ
              do j=1,NGLLY
                 iglob = ibool(i,j,k,ispec)

                 accel(1,iglob) = 0._CUSTOM_REAL
                 accel(2,iglob) = 0._CUSTOM_REAL
                 accel(3,iglob) = 0._CUSTOM_REAL

                 veloc(1,iglob) = 0._CUSTOM_REAL
                 veloc(2,iglob) = 0._CUSTOM_REAL
                 veloc(3,iglob) = 0._CUSTOM_REAL

                 displ(1,iglob) = 0._CUSTOM_REAL
                 displ(2,iglob) = 0._CUSTOM_REAL
                 displ(3,iglob) = 0._CUSTOM_REAL
              enddo
           enddo
        endif
     enddo

     ! xmax
     do ispec2D=1,nspec2D_xmax
        ispec = ibelm_xmax(ispec2D)

        if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
           i = NGLLX

           do k=1,NGLLZ
              do j=1,NGLLY
                 iglob = ibool(i,j,k,ispec)

                 accel(1,iglob) = 0._CUSTOM_REAL
                 accel(2,iglob) = 0._CUSTOM_REAL
                 accel(3,iglob) = 0._CUSTOM_REAL

                 veloc(1,iglob) = 0._CUSTOM_REAL
                 veloc(2,iglob) = 0._CUSTOM_REAL
                 veloc(3,iglob) = 0._CUSTOM_REAL

                 displ(1,iglob) = 0._CUSTOM_REAL
                 displ(2,iglob) = 0._CUSTOM_REAL
                 displ(3,iglob) = 0._CUSTOM_REAL
              enddo
           enddo
        endif
     enddo

     ! ymin
     do ispec2D=1,nspec2D_ymin
        ispec = ibelm_ymin(ispec2D)

        if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
           j = 1

           do k=1,NGLLZ
              do i=1,NGLLX
                 iglob = ibool(i,j,k,ispec)

                 accel(1,iglob) = 0._CUSTOM_REAL
                 accel(2,iglob) = 0._CUSTOM_REAL
                 accel(3,iglob) = 0._CUSTOM_REAL

                 veloc(1,iglob) = 0._CUSTOM_REAL
                 veloc(2,iglob) = 0._CUSTOM_REAL
                 veloc(3,iglob) = 0._CUSTOM_REAL

                 displ(1,iglob) = 0._CUSTOM_REAL
                 displ(2,iglob) = 0._CUSTOM_REAL
                 displ(3,iglob) = 0._CUSTOM_REAL
              enddo
           enddo
        endif
     enddo

     ! ymax
     do ispec2D=1,nspec2D_ymax
        ispec = ibelm_ymax(ispec2D)

        if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
           j = NGLLY

           do k=1,NGLLZ
              do i=1,NGLLX
                 iglob = ibool(i,j,k,ispec)

                 accel(1,iglob) = 0._CUSTOM_REAL
                 accel(2,iglob) = 0._CUSTOM_REAL
                 accel(3,iglob) = 0._CUSTOM_REAL

                 veloc(1,iglob) = 0._CUSTOM_REAL
                 veloc(2,iglob) = 0._CUSTOM_REAL
                 veloc(3,iglob) = 0._CUSTOM_REAL

                 displ(1,iglob) = 0._CUSTOM_REAL
                 displ(2,iglob) = 0._CUSTOM_REAL
                 displ(3,iglob) = 0._CUSTOM_REAL
              enddo
           enddo
        endif
     enddo

     ! bottom (zmin)
     do ispec2D=1,NSPEC2D_BOTTOM
        ispec = ibelm_bottom(ispec2D)

        if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
           k = 1

           do j=1,NGLLY
              do i=1,NGLLX
                 iglob = ibool(i,j,k,ispec)

                 accel(1,iglob) = 0._CUSTOM_REAL
                 accel(2,iglob) = 0._CUSTOM_REAL
                 accel(3,iglob) = 0._CUSTOM_REAL

                 veloc(1,iglob) = 0._CUSTOM_REAL
                 veloc(2,iglob) = 0._CUSTOM_REAL
                 veloc(3,iglob) = 0._CUSTOM_REAL

                 displ(1,iglob) = 0._CUSTOM_REAL
                 displ(2,iglob) = 0._CUSTOM_REAL
                 displ(3,iglob) = 0._CUSTOM_REAL
              enddo
           enddo
        endif
     enddo

     ! top (zmax)
     if(PML_INSTEAD_OF_FREE_SURFACE)then
       do ispec2D=1,NSPEC2D_BOTTOM
          ispec = ibelm_top(ispec2D)

          if(CPML_mask_ibool(ispec) .and. ispec_is_elastic(ispec)) then
             k = NGLLZ

             do j=1,NGLLY
                do i=1,NGLLX
                   iglob = ibool(i,j,k,ispec)

                   accel(1,iglob) = 0._CUSTOM_REAL
                   accel(2,iglob) = 0._CUSTOM_REAL
                   accel(3,iglob) = 0._CUSTOM_REAL

                   veloc(1,iglob) = 0._CUSTOM_REAL
                   veloc(2,iglob) = 0._CUSTOM_REAL
                   veloc(3,iglob) = 0._CUSTOM_REAL

                   displ(1,iglob) = 0._CUSTOM_REAL
                   displ(2,iglob) = 0._CUSTOM_REAL
                   displ(3,iglob) = 0._CUSTOM_REAL
                enddo
             enddo
          endif
       enddo
     endif

  endif ! if( PML_CONDITIONS )

end subroutine compute_forces_viscoelastic_noDev

