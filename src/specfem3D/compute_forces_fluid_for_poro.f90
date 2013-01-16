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


  subroutine compute_forces_poro_fluid_part( iphase, &
                        NSPEC_AB,NGLOB_AB,displw_poroelastic,accelw_poroelastic,&
                        velocw_poroelastic,displs_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,&
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll,  &
                        kappaarraystore,rhoarraystore,mustore,etastore,permstore, &
                        phistore,tortstore,jacobian,ibool,&
                        epsilonwdev_xx,epsilonwdev_yy,epsilonwdev_xy,&
                        epsilonwdev_xz,epsilonwdev_yz,epsilonw_trace_over_3, &
                        SIMULATION_TYPE,NSPEC_ADJOINT, &
                        num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic,&
                        phase_ispec_inner_poroelastic )

! compute forces for the fluid poroelastic part

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM, &
                      N_SLS, &
                      ONE_THIRD,FOUR_THIRDS

  implicit none

  integer :: iphase
  integer :: NSPEC_AB,NGLOB_AB

! adjoint simulations
  integer :: SIMULATION_TYPE
  !integer :: NSPEC_BOUN
  integer :: NSPEC_ADJOINT
! adjoint wavefields
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: &
!    mufr_kl, B_kl

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displw_poroelastic,accelw_poroelastic,&
                                                      velocw_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displs_poroelastic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: &
       epsilonwdev_xx,epsilonwdev_yy,epsilonwdev_xy,epsilonwdev_xz,epsilonwdev_yz
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: epsilonw_trace_over_3

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        mustore,etastore,phistore,tortstore,jacobian
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappaarraystore
  real(kind=CUSTOM_REAL), dimension(6,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: permstore

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLY) :: wygll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

!  logical,dimension(NSPEC_AB) :: ispec_is_elastic
  integer :: num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic
  integer, dimension(num_phase_ispec_poroelastic,2) :: phase_ispec_inner_poroelastic

!---
!--- local variables
!---

  integer :: ispec,i,j,k,l,iglob,num_elements,ispec_p

! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
    tempx1p,tempx2p,tempx3p,tempy1p,tempy2p,tempy3p,tempz1p,tempz2p,tempz3p
!    b_tempx1,b_tempx2,b_tempx3,b_tempy1,b_tempy2,b_tempy3,b_tempz1,b_tempz2,b_tempz3, &
!    b_tempx1p,b_tempx2p,b_tempx3p,b_tempy1p,b_tempy2p,b_tempy3p,b_tempz1p,b_tempz2p,b_tempz3p

  real(kind=CUSTOM_REAL) :: duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl
  real(kind=CUSTOM_REAL) :: dwxdxl,dwydxl,dwzdxl,dwxdyl,dwydyl,dwzdyl,dwxdzl,dwydzl,dwzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl_plus_duzdzl,dwxdxl_plus_dwydyl_plus_dwzdzl

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  real(kind=CUSTOM_REAL) tempx1ls,tempx2ls,tempx3ls,tempx1lw,tempx2lw,tempx3lw
  real(kind=CUSTOM_REAL) tempy1ls,tempy2ls,tempy3ls,tempy1lw,tempy2lw,tempy3lw
  real(kind=CUSTOM_REAL) tempz1ls,tempz2ls,tempz3ls,tempz1lw,tempz2lw,tempz3lw

!  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
!  real(kind=CUSTOM_REAL) :: dsxx,dsxy,dsxz,dsyy,dsyz,dszz
!  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxy,b_dsxz,b_dsyy,b_dsyz,b_dszz
!  real(kind=CUSTOM_REAL) :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl
!  real(kind=CUSTOM_REAL) :: dwxx,dwxy,dwxz,dwyy,dwyz,dwzz
!  real(kind=CUSTOM_REAL) :: b_dwxx,b_dwxy,b_dwxz,b_dwyy,b_dwyz,b_dwzz
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: sigmap
!  real(kind=CUSTOM_REAL) :: b_sigma_xx,b_sigma_yy,b_sigma_zz,b_sigma_xy,b_sigma_xz,b_sigma_yz
!  real(kind=CUSTOM_REAL) :: b_sigmap
!  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,vxf,vzf,vnf,rho_vpI,rho_vpII,rho_vs,tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D

! viscous attenuation (poroelastic media)
  real(kind=CUSTOM_REAL), dimension(6) :: bl_relaxed


! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) ::  xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

! material properties of the poroelastic medium
!  real(kind=CUSTOM_REAL) :: mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed
  real(kind=CUSTOM_REAL) :: kappal_s,rhol_s
  real(kind=CUSTOM_REAL) :: etal_f,kappal_f,rhol_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr,phil,tortl,viscodampx,viscodampy,viscodampz
  real(kind=CUSTOM_REAL) :: permlxx,permlxy,permlxz,permlyz,permlyy,permlzz,&
                            invpermlxx,invpermlxy,invpermlxz,invpermlyz,invpermlyy,invpermlzz,detk
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,rhol_bar

  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G
!  real(kind=CUSTOM_REAL) :: cpIsquare,cpIIsquare,cssquare,cpIl,cpIIl,csl

! for attenuation
!  real(kind=CUSTOM_REAL) :: Un,Unp1,tauinv,Sn,Snp1,theta_n,theta_np1,tauinvsquare,tauinvcube,tauinvUn

! compute Grad(displs_poroelastic) at time step n for attenuation
!  if(TURN_ATTENUATION_ON) call compute_gradient_attenuation(displs_poroelastic,dux_dxl_n,duz_dxl_n, &
!      dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,poroelastic,hprime_xx,hprime_zz,nspec,npoin)


  if( iphase == 1 ) then
    num_elements = nspec_outer_poroelastic
  else
    num_elements = nspec_inner_poroelastic
  endif

! loop over spectral elements
  do ispec_p = 1,num_elements

        ispec = phase_ispec_inner_poroelastic(ispec_p,iphase)

! first double loop over GLL points to compute and store gradients
    do k=1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

! get poroelastic parameters of current local GLL
    phil = phistore(i,j,k,ispec)
    tortl = tortstore(i,j,k,ispec)
!solid properties
    kappal_s = kappaarraystore(1,i,j,k,ispec)
    rhol_s = rhoarraystore(1,i,j,k,ispec)
!fluid properties
    kappal_f = kappaarraystore(2,i,j,k,ispec)
    rhol_f = rhoarraystore(2,i,j,k,ispec)
!frame properties
    mul_fr = mustore(i,j,k,ispec)
    kappal_fr = kappaarraystore(3,i,j,k,ispec)
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
!The RHS has the form : div T -phi/c div T_f + phi/ceta_fk^-1.partial t w
!where T = G:grad u_s + C_biot div w I
!and T_f = C_biot div u_s I + M_biot div w I
      mul_G = mul_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
      lambdalplus2mul_G = lambdal_G + 2._CUSTOM_REAL*mul_G

! derivative along x,y,z for u_s and w
              tempx1ls = 0.
              tempx2ls = 0.
              tempx3ls = 0.

              tempy1ls = 0.
              tempy2ls = 0.
              tempy3ls = 0.

              tempz1ls = 0.
              tempz2ls = 0.
              tempz3ls = 0.

              tempx1lw = 0.
              tempx2lw = 0.
              tempx3lw = 0.

              tempy1lw = 0.
              tempy2lw = 0.
              tempy3lw = 0.

              tempz1lw = 0.
              tempz2lw = 0.
              tempz3lw = 0.

! first double loop over GLL points to compute and store gradients
          do l = 1,NGLLX
                hp1 = hprime_xx(i,l)
                iglob = ibool(l,j,k,ispec)
                tempx1ls = tempx1ls + displs_poroelastic(1,iglob)*hp1
                tempy1ls = tempy1ls + displs_poroelastic(2,iglob)*hp1
                tempz1ls = tempz1ls + displs_poroelastic(3,iglob)*hp1
                tempx1lw = tempx1lw + displw_poroelastic(1,iglob)*hp1
                tempy1lw = tempy1lw + displw_poroelastic(2,iglob)*hp1
                tempz1lw = tempz1lw + displw_poroelastic(3,iglob)*hp1
    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                hp2 = hprime_yy(j,l)
                iglob = ibool(i,l,k,ispec)
                tempx2ls = tempx2ls + displs_poroelastic(1,iglob)*hp2
                tempy2ls = tempy2ls + displs_poroelastic(2,iglob)*hp2
                tempz2ls = tempz2ls + displs_poroelastic(3,iglob)*hp2
                tempx2lw = tempx2lw + displw_poroelastic(1,iglob)*hp2
                tempy2lw = tempy2lw + displw_poroelastic(2,iglob)*hp2
                tempz2lw = tempz2lw + displw_poroelastic(3,iglob)*hp2
    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                hp3 = hprime_zz(k,l)
                iglob = ibool(i,j,l,ispec)
                tempx3ls = tempx3ls + displs_poroelastic(1,iglob)*hp3
                tempy3ls = tempy3ls + displs_poroelastic(2,iglob)*hp3
                tempz3ls = tempz3ls + displs_poroelastic(3,iglob)*hp3
                tempx3lw = tempx3lw + displw_poroelastic(1,iglob)*hp3
                tempy3lw = tempy3lw + displw_poroelastic(2,iglob)*hp3
                tempz3lw = tempz3lw + displw_poroelastic(3,iglob)*hp3
          enddo

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

! derivatives of displacement
              duxdxl = xixl*tempx1ls + etaxl*tempx2ls + gammaxl*tempx3ls
              duxdyl = xiyl*tempx1ls + etayl*tempx2ls + gammayl*tempx3ls
              duxdzl = xizl*tempx1ls + etazl*tempx2ls + gammazl*tempx3ls

              duydxl = xixl*tempy1ls + etaxl*tempy2ls + gammaxl*tempy3ls
              duydyl = xiyl*tempy1ls + etayl*tempy2ls + gammayl*tempy3ls
              duydzl = xizl*tempy1ls + etazl*tempy2ls + gammazl*tempy3ls

              duzdxl = xixl*tempz1ls + etaxl*tempz2ls + gammaxl*tempz3ls
              duzdyl = xiyl*tempz1ls + etayl*tempz2ls + gammayl*tempz3ls
              duzdzl = xizl*tempz1ls + etazl*tempz2ls + gammazl*tempz3ls

              dwxdxl = xixl*tempx1lw + etaxl*tempx2lw + gammaxl*tempx3lw
              dwxdyl = xiyl*tempx1lw + etayl*tempx2lw + gammayl*tempx3lw
              dwxdzl = xizl*tempx1lw + etazl*tempx2lw + gammazl*tempx3lw

              dwydxl = xixl*tempy1lw + etaxl*tempy2lw + gammaxl*tempy3lw
              dwydyl = xiyl*tempy1lw + etayl*tempy2lw + gammayl*tempy3lw
              dwydzl = xizl*tempy1lw + etazl*tempy2lw + gammazl*tempy3lw

              dwzdxl = xixl*tempz1lw + etaxl*tempz2lw + gammaxl*tempz3lw
              dwzdyl = xiyl*tempz1lw + etayl*tempz2lw + gammayl*tempz3lw
              dwzdzl = xizl*tempz1lw + etazl*tempz2lw + gammazl*tempz3lw

    ! precompute some sums to save CPU time
              duxdxl_plus_duydyl_plus_duzdzl = duxdxl + duydyl + duzdzl
              dwxdxl_plus_dwydyl_plus_dwzdzl = dwxdxl + dwydyl + dwzdzl
              duxdxl_plus_duydyl = duxdxl + duydyl
              duxdxl_plus_duzdzl = duxdxl + duzdzl
              duydyl_plus_duzdzl = duydyl + duzdzl
              duxdyl_plus_duydxl = duxdyl + duydxl
              duzdxl_plus_duxdzl = duzdxl + duxdzl
              duzdyl_plus_duydzl = duzdyl + duydzl

! compute stress tensor (include attenuation or anisotropy if needed)

!  if(VISCOATTENUATION) then
!chris:check

! Dissipation only controlled by frame share attenuation in poroelastic (see Morency & Tromp, GJI 2008).
! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

!  else

! no attenuation
    sigma_xx = lambdalplus2mul_G*duxdxl + lambdal_G*duydyl_plus_duzdzl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl
    sigma_yy = lambdalplus2mul_G*duydyl + lambdal_G*duxdxl_plus_duzdzl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl
    sigma_zz = lambdalplus2mul_G*duzdzl + lambdal_G*duxdxl_plus_duydyl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl

    sigma_xy = mul_G*duxdyl_plus_duydxl
    sigma_xz = mul_G*duzdxl_plus_duxdzl
    sigma_yz = mul_G*duzdyl_plus_duydzl

    sigmap = C_biot*duxdxl_plus_duydyl_plus_duzdzl + M_biot*dwxdxl_plus_dwydyl_plus_dwzdzl

          if(SIMULATION_TYPE == 3) then ! kernels calculation
    epsilonw_trace_over_3(i,j,k,ispec) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilonwdev_xx(i,j,k,ispec) = duxdxl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilonwdev_yy(i,j,k,ispec) = duydyl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilonwdev_xy(i,j,k,ispec) = 0.5 * duxdyl_plus_duydxl
    epsilonwdev_xz(i,j,k,ispec) = 0.5 * duzdxl_plus_duxdzl
    epsilonwdev_yz(i,j,k,ispec) = 0.5 * duzdyl_plus_duydzl
          endif
!  endif !if(VISCOATTENUATION)

! weak formulation term based on stress tensor (non-symmetric form)
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

          tempx1p(i,j,k) = jacobianl * sigmap*xixl
          tempy1p(i,j,k) = jacobianl * sigmap*xiyl
          tempz1p(i,j,k) = jacobianl * sigmap*xizl

          tempx2p(i,j,k) = jacobianl * sigmap*etaxl
          tempy2p(i,j,k) = jacobianl * sigmap*etayl
          tempz2p(i,j,k) = jacobianl * sigmap*etazl

          tempx3p(i,j,k) = jacobianl * sigmap*gammaxl
          tempy3p(i,j,k) = jacobianl * sigmap*gammayl
          tempz3p(i,j,k) = jacobianl * sigmap*gammazl

        enddo
      enddo
    enddo

!
! second double-loop over GLL to compute all the terms
!
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

              tempx1ls = 0.
              tempy1ls = 0.
              tempz1ls = 0.

              tempx2ls = 0.
              tempy2ls = 0.
              tempz2ls = 0.

              tempx3ls = 0.
              tempy3ls = 0.
              tempz3ls = 0.

              tempx1lw = 0.
              tempy1lw = 0.
              tempz1lw = 0.

              tempx2lw = 0.
              tempy2lw = 0.
              tempz2lw = 0.

              tempx3lw = 0.
              tempy3lw = 0.
              tempz3lw = 0.

              do l=1,NGLLX
                fac1 = hprimewgll_xx(l,i)
                tempx1ls = tempx1ls + tempx1(l,j,k)*fac1
                tempy1ls = tempy1ls + tempy1(l,j,k)*fac1
                tempz1ls = tempz1ls + tempz1(l,j,k)*fac1
                tempx1lw = tempx1lw + tempx1p(l,j,k)*fac1
                tempy1lw = tempy1lw + tempy1p(l,j,k)*fac1
                tempz1lw = tempz1lw + tempz1p(l,j,k)*fac1
                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                fac2 = hprimewgll_yy(l,j)
                tempx2ls = tempx2ls + tempx2(i,l,k)*fac2
                tempy2ls = tempy2ls + tempy2(i,l,k)*fac2
                tempz2ls = tempz2ls + tempz2(i,l,k)*fac2
                tempx2lw = tempx2lw + tempx2p(i,l,k)*fac2
                tempy2lw = tempy2lw + tempy2p(i,l,k)*fac2
                tempz2lw = tempz2lw + tempz2p(i,l,k)*fac2
                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

                !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                fac3 = hprimewgll_zz(l,k)
                tempx3ls = tempx3ls + tempx3(i,j,l)*fac3
                tempy3ls = tempy3ls + tempy3(i,j,l)*fac3
                tempz3ls = tempz3ls + tempz3(i,j,l)*fac3
                tempx3lw = tempx3lw + tempx3p(i,j,l)*fac3
                tempy3lw = tempy3lw + tempy3p(i,j,l)*fac3
                tempz3lw = tempz3lw + tempz3p(i,j,l)*fac3
              enddo

              fac1 = wgllwgll_yz(j,k)
              fac2 = wgllwgll_xz(i,k)
              fac3 = wgllwgll_xy(i,j)

! get poroelastic parameters of current local GLL
    phil = phistore(i,j,k,ispec)
!solid properties
    rhol_s = rhoarraystore(1,i,j,k,ispec)
!fluid properties
    rhol_f = rhoarraystore(2,i,j,k,ispec)
!frame properties
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

    ! sum contributions from each element to the global mesh

              iglob = ibool(i,j,k,ispec)


    accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + ( fac1*(rhol_f/rhol_bar*tempx1ls - tempx1lw) &
           + fac2*(rhol_f/rhol_bar*tempx2ls - tempx2lw) + fac3*(rhol_f/rhol_bar*tempx3ls - tempx3lw) )

    accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + ( fac1*(rhol_f/rhol_bar*tempy1ls - tempy1lw) &
           + fac2*(rhol_f/rhol_bar*tempy2ls - tempy2lw) + fac3*(rhol_f/rhol_bar*tempy3ls - tempy3lw) )

    accelw_poroelastic(3,iglob) = accelw_poroelastic(3,iglob) + ( fac1*(rhol_f/rhol_bar*tempz1ls - tempz1lw) &
           + fac2*(rhol_f/rhol_bar*tempz2ls - tempz2lw) + fac3*(rhol_f/rhol_bar*tempz3ls - tempz3lw) )


!
!---- viscous damping
!
! add + phi/tort eta_f k^-1 dot(w)

    etal_f = etastore(i,j,k,ispec)

      if(etal_f >0.d0) then

    permlxx = permstore(1,i,j,k,ispec)
    permlxy = permstore(2,i,j,k,ispec)
    permlxz = permstore(3,i,j,k,ispec)
    permlyy = permstore(4,i,j,k,ispec)
    permlyz = permstore(5,i,j,k,ispec)
    permlzz = permstore(6,i,j,k,ispec)

! calcul of the inverse of k
    detk = permlxz*(permlxy*permlyz-permlxz*permlyy) &
         - permlxy*(permlxy*permlzz-permlyz*permlxz) &
         + permlxx*(permlyy*permlzz-permlyz*permlyz)

    if(detk /= 0.d0) then
     invpermlxx = (permlyy*permlzz-permlyz*permlyz)/detk
     invpermlxy = (permlxz*permlyz-permlxy*permlzz)/detk
     invpermlxz = (permlxy*permlyz-permlxz*permlyy)/detk
     invpermlyy = (permlxx*permlzz-permlxz*permlxz)/detk
     invpermlyz = (permlxy*permlxz-permlxx*permlyz)/detk
     invpermlzz = (permlxx*permlyy-permlxy*permlxy)/detk
    else
      stop 'Permeability matrix is not inversible'
    endif

! relaxed viscous coef
          bl_relaxed(1) = etal_f*invpermlxx
          bl_relaxed(2) = etal_f*invpermlxy
          bl_relaxed(3) = etal_f*invpermlxz
          bl_relaxed(4) = etal_f*invpermlyy
          bl_relaxed(5) = etal_f*invpermlyz
          bl_relaxed(6) = etal_f*invpermlzz

!    if(VISCOATTENUATION) then
!          bl_unrelaxed(1) = etal_f*invpermlxx*theta_e/theta_s
!          bl_unrelaxed(2) = etal_f*invpermlxz*theta_e/theta_s
!          bl_unrelaxed(3) = etal_f*invpermlzz*theta_e/theta_s
!    endif

!     do k = 1,NGLLZ
!      do j = 1,NGLLY
!        do i = 1,NGLLX

!              iglob = ibool(i,j,k,ispec)

!     if(VISCOATTENUATION) then
! compute the viscous damping term with the unrelaxed viscous coef and add memory variable
!      viscodampx = velocw_poroelastic(1,iglob)*bl_unrelaxed(1) + velocw_poroelastic(2,iglob)*bl_unrelaxed(2)&
!                  - rx_viscous(i,j,ispec)
!      viscodampz = velocw_poroelastic(1,iglob)*bl_unrelaxed(2) + velocw_poroelastic(2,iglob)*bl_unrelaxed(3)&
!                  - rz_viscous(i,j,ispec)
!     else

! no viscous attenuation
      viscodampx = velocw_poroelastic(1,iglob)*bl_relaxed(1) + velocw_poroelastic(2,iglob)*bl_relaxed(2) + &
                   velocw_poroelastic(3,iglob)*bl_relaxed(3)
      viscodampy = velocw_poroelastic(1,iglob)*bl_relaxed(2) + velocw_poroelastic(2,iglob)*bl_relaxed(4) + &
                   velocw_poroelastic(3,iglob)*bl_relaxed(5)
      viscodampz = velocw_poroelastic(1,iglob)*bl_relaxed(3) + velocw_poroelastic(2,iglob)*bl_relaxed(5) + &
                   velocw_poroelastic(3,iglob)*bl_relaxed(6)
!     endif

     accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - wxgll(i)*wygll(j)*wzgll(k)*jacobian(i,j,k,ispec)*&
              viscodampx
     accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - wxgll(i)*wygll(j)*wzgll(k)*jacobian(i,j,k,ispec)*&
              viscodampy
     accelw_poroelastic(3,iglob) = accelw_poroelastic(3,iglob) - wxgll(i)*wygll(j)*wzgll(k)*jacobian(i,j,k,ispec)*&
              viscodampz

! if isolver == 1 .and. save_forward then b_viscodamp is saved in compute_forces_poro_fluid_part.f90
!          if(isolver == 2) then ! kernels calculation
!        b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + phil/tortl*b_viscodampx(iglob)
!        b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + phil/tortl*b_viscodampz(iglob)
!          endif

!        enddo
!      enddo
!     enddo

         endif ! if(etal_f >0.d0) then

        enddo ! second loop over the GLL points
      enddo
    enddo

    enddo ! end of loop over all spectral elements


  end subroutine compute_forces_poro_fluid_part

