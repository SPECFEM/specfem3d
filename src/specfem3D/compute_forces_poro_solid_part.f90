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

  subroutine compute_forces_poro_solid_part(iphase, &
                                            displs_poroelastic,accels_poroelastic, &
                                            displw_poroelastic,velocw_poroelastic, &
                                            epsilonsdev_xx,epsilonsdev_yy,epsilonsdev_xy, &
                                            epsilonsdev_xz,epsilonsdev_yz,epsilons_trace_over_3)

! compute forces for the solid poroelastic part

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ONE_THIRD

  use specfem_par, only: NGLOB_AB, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,jacobianstore, &
                         hprime_xx,hprime_yy,hprime_zz, &
                         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wxgll,wygll,wzgll, &
                         SIMULATION_TYPE,NSPEC_ADJOINT,ibool,mustore, &
                         irregular_element_number,xix_regular,jacobian_regular

  use specfem_par_poroelastic, only: kappaarraystore,rhoarraystore,etastore,permstore, &
                                     phistore,tortstore,nspec_inner_poroelastic, &
                                     nspec_outer_poroelastic,phase_ispec_inner_poroelastic
  implicit none

  integer :: iphase

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displs_poroelastic,accels_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displw_poroelastic,velocw_poroelastic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: &
       epsilonsdev_xx,epsilonsdev_yy,epsilonsdev_xy,epsilonsdev_xz,epsilonsdev_yz
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: epsilons_trace_over_3

  ! local parameters
  integer :: ispec,i,j,k,l,iglob,num_elements,ispec_p,ispec_irreg

  ! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
    tempx1p,tempx2p,tempx3p,tempy1p,tempy2p,tempy3p,tempz1p,tempz2p,tempz3p

  real(kind=CUSTOM_REAL) :: duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl
  real(kind=CUSTOM_REAL) :: dwxdxl,dwydxl,dwzdxl,dwxdyl,dwydyl,dwzdyl,dwxdzl,dwydzl,dwzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl_plus_duzdzl,dwxdxl_plus_dwydyl_plus_dwzdzl

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL) :: tempx1ls,tempx2ls,tempx3ls,tempx1lw,tempx2lw,tempx3lw
  real(kind=CUSTOM_REAL) :: tempy1ls,tempy2ls,tempy3ls,tempy1lw,tempy2lw,tempy3lw
  real(kind=CUSTOM_REAL) :: tempz1ls,tempz2ls,tempz3ls,tempz1lw,tempz2lw,tempz3lw

  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: sigmap

! viscous attenuation (poroelastic media)
  real(kind=CUSTOM_REAL), dimension(6) :: bl_relaxed

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) ::  xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

! material properties of the poroelastic medium
  real(kind=CUSTOM_REAL) :: kappal_s,rhol_s
  real(kind=CUSTOM_REAL) :: etal_f,kappal_f,rhol_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr,phil,tortl,viscodampx,viscodampy,viscodampz
  real(kind=CUSTOM_REAL) :: permlxx,permlxy,permlxz,permlyz,permlyy,permlzz, &
                            invpermlxx,invpermlxy,invpermlxz,invpermlyz,invpermlyy,invpermlzz,detk
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,rhol_bar
  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G

  if (iphase == 1) then
    num_elements = nspec_outer_poroelastic
  else
    num_elements = nspec_inner_poroelastic
  endif

! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_poroelastic(ispec_p,iphase)

    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg == 0) jacobianl = jacobian_regular
 !
    ! first double loop over GLL points to compute and store gradients
    !
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

          !The RHS has the form : div T -phi/c div T_f + phi/ceta_f_k^-1.partial t w
          !where T = G:grad u_s + C_biot div w I
          !and T_f = C_biot div u_s I + M_biot div w I
          mul_G = mul_fr
          lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
          lambdalplus2mul_G = lambdal_G + 2._CUSTOM_REAL*mul_G

          ! derivative along x,y,z for u_s and w
          tempx1ls = 0.0_CUSTOM_REAL
          tempx2ls = 0.0_CUSTOM_REAL
          tempx3ls = 0.0_CUSTOM_REAL

          tempy1ls = 0.0_CUSTOM_REAL
          tempy2ls = 0.0_CUSTOM_REAL
          tempy3ls = 0.0_CUSTOM_REAL

          tempz1ls = 0.0_CUSTOM_REAL
          tempz2ls = 0.0_CUSTOM_REAL
          tempz3ls = 0.0_CUSTOM_REAL

          tempx1lw = 0.0_CUSTOM_REAL
          tempx2lw = 0.0_CUSTOM_REAL
          tempx3lw = 0.0_CUSTOM_REAL

          tempy1lw = 0.0_CUSTOM_REAL
          tempy2lw = 0.0_CUSTOM_REAL
          tempy3lw = 0.0_CUSTOM_REAL

          tempz1lw = 0.0_CUSTOM_REAL
          tempz2lw = 0.0_CUSTOM_REAL
          tempz3lw = 0.0_CUSTOM_REAL

!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here
!! DK DK Oct 2018: we could (and should) use the Deville matrix products instead here

          ! first double loop over GLL points to compute and store gradients
          ! we can merge these loops because NGLLX = NGLLY = NGLLZ
          do l = 1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1ls = tempx1ls + displs_poroelastic(1,iglob)*hp1
            tempy1ls = tempy1ls + displs_poroelastic(2,iglob)*hp1
            tempz1ls = tempz1ls + displs_poroelastic(3,iglob)*hp1
            tempx1lw = tempx1lw + displw_poroelastic(1,iglob)*hp1
            tempy1lw = tempy1lw + displw_poroelastic(2,iglob)*hp1
            tempz1lw = tempz1lw + displw_poroelastic(3,iglob)*hp1

            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2ls = tempx2ls + displs_poroelastic(1,iglob)*hp2
            tempy2ls = tempy2ls + displs_poroelastic(2,iglob)*hp2
            tempz2ls = tempz2ls + displs_poroelastic(3,iglob)*hp2
            tempx2lw = tempx2lw + displw_poroelastic(1,iglob)*hp2
            tempy2lw = tempy2lw + displw_poroelastic(2,iglob)*hp2
            tempz2lw = tempz2lw + displw_poroelastic(3,iglob)*hp2

            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3ls = tempx3ls + displs_poroelastic(1,iglob)*hp3
            tempy3ls = tempy3ls + displs_poroelastic(2,iglob)*hp3
            tempz3ls = tempz3ls + displs_poroelastic(3,iglob)*hp3
            tempx3lw = tempx3lw + displw_poroelastic(1,iglob)*hp3
            tempy3lw = tempy3lw + displw_poroelastic(2,iglob)*hp3
            tempz3lw = tempz3lw + displw_poroelastic(3,iglob)*hp3
          enddo

          ! get derivatives of ux, uy and uz with respect to x, y and z
          if (ispec_irreg /= 0) then
            !irregular element
            xixl = xixstore(i,j,k,ispec_irreg)
            xiyl = xiystore(i,j,k,ispec_irreg)
            xizl = xizstore(i,j,k,ispec_irreg)
            etaxl = etaxstore(i,j,k,ispec_irreg)
            etayl = etaystore(i,j,k,ispec_irreg)
            etazl = etazstore(i,j,k,ispec_irreg)
            gammaxl = gammaxstore(i,j,k,ispec_irreg)
            gammayl = gammaystore(i,j,k,ispec_irreg)
            gammazl = gammazstore(i,j,k,ispec_irreg)
            jacobianl = jacobianstore(i,j,k,ispec_irreg)

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

          else
            !regular element
            ! derivatives of displacement
            duxdxl = xix_regular*tempx1ls
            duxdyl = xix_regular*tempx2ls
            duxdzl = xix_regular*tempx3ls

            duydxl = xix_regular*tempy1ls
            duydyl = xix_regular*tempy2ls
            duydzl = xix_regular*tempy3ls

            duzdxl = xix_regular*tempz1ls
            duzdyl = xix_regular*tempz2ls
            duzdzl = xix_regular*tempz3ls

            dwxdxl = xix_regular*tempx1lw
            dwxdyl = xix_regular*tempx2lw
            dwxdzl = xix_regular*tempx3lw

            dwydxl = xix_regular*tempy1lw
            dwydyl = xix_regular*tempy2lw
            dwydzl = xix_regular*tempy3lw

            dwzdxl = xix_regular*tempz1lw
            dwzdyl = xix_regular*tempz2lw
            dwzdzl = xix_regular*tempz3lw
          endif

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

          !  if (VISCOATTENUATION) then
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

          !  endif !if (ATTENUATION)

          if (SIMULATION_TYPE == 3) then ! kernels calculation
            epsilons_trace_over_3(i,j,k,ispec) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
            epsilonsdev_xx(i,j,k,ispec) = duxdxl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
            epsilonsdev_yy(i,j,k,ispec) = duydyl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
            epsilonsdev_xy(i,j,k,ispec) = 0.5 * duxdyl_plus_duydxl
            epsilonsdev_xz(i,j,k,ispec) = 0.5 * duzdxl_plus_duxdzl
            epsilonsdev_yz(i,j,k,ispec) = 0.5 * duzdyl_plus_duydzl
          endif

          ! weak formulation term based on stress tensor (non-symmetric form)
          ! define symmetric components of sigma
          sigma_yx = sigma_xy
          sigma_zx = sigma_xz
          sigma_zy = sigma_yz

          ! form dot product with test vector, non-symmetric form (which is
          ! useful in the case of PML)
          if (ispec_irreg /= 0) then
            ! irregular element
            tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
            tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
            tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
            tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
            tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

            tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
            tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
            tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

            tempx1p(i,j,k) = jacobianl * sigmap * xixl
            tempy1p(i,j,k) = jacobianl * sigmap * xiyl
            tempz1p(i,j,k) = jacobianl * sigmap * xizl

            tempx2p(i,j,k) = jacobianl * sigmap * etaxl
            tempy2p(i,j,k) = jacobianl * sigmap * etayl
            tempz2p(i,j,k) = jacobianl * sigmap * etazl

            tempx3p(i,j,k) = jacobianl * sigmap * gammaxl
            tempy3p(i,j,k) = jacobianl * sigmap * gammayl
            tempz3p(i,j,k) = jacobianl * sigmap * gammazl

          else
            ! regular element
            tempx1(i,j,k) = jacobianl * sigma_xx * xix_regular
            tempy1(i,j,k) = jacobianl * sigma_xy * xix_regular
            tempz1(i,j,k) = jacobianl * sigma_xz * xix_regular

            tempx2(i,j,k) = jacobianl * sigma_yx * xix_regular
            tempy2(i,j,k) = jacobianl * sigma_yy * xix_regular
            tempz2(i,j,k) = jacobianl * sigma_yz * xix_regular

            tempx3(i,j,k) = jacobianl * sigma_zx * xix_regular
            tempy3(i,j,k) = jacobianl * sigma_zy * xix_regular
            tempz3(i,j,k) = jacobianl * sigma_zz * xix_regular

            tempx1p(i,j,k) = jacobianl * sigmap * xix_regular
            tempy1p(i,j,k) = 0.0_CUSTOM_REAL
            tempz1p(i,j,k) = 0.0_CUSTOM_REAL

            tempx2p(i,j,k) = 0.0_CUSTOM_REAL
            tempy2p(i,j,k) = jacobianl * sigmap * xix_regular
            tempz2p(i,j,k) = 0.0_CUSTOM_REAL

            tempx3p(i,j,k) = 0.0_CUSTOM_REAL
            tempy3p(i,j,k) = 0.0_CUSTOM_REAL
            tempz3p(i,j,k) = jacobianl * sigmap * xix_regular
          endif

        enddo
      enddo
    enddo

    !
    ! second double-loop over GLL to compute all the terms
    !
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1ls = 0.0_CUSTOM_REAL
          tempy1ls = 0.0_CUSTOM_REAL
          tempz1ls = 0.0_CUSTOM_REAL

          tempx2ls = 0.0_CUSTOM_REAL
          tempy2ls = 0.0_CUSTOM_REAL
          tempz2ls = 0.0_CUSTOM_REAL

          tempx3ls = 0.0_CUSTOM_REAL
          tempy3ls = 0.0_CUSTOM_REAL
          tempz3ls = 0.0_CUSTOM_REAL

          tempx1lw = 0.0_CUSTOM_REAL
          tempy1lw = 0.0_CUSTOM_REAL
          tempz1lw = 0.0_CUSTOM_REAL

          tempx2lw = 0.0_CUSTOM_REAL
          tempy2lw = 0.0_CUSTOM_REAL
          tempz2lw = 0.0_CUSTOM_REAL

          tempx3lw = 0.0_CUSTOM_REAL
          tempy3lw = 0.0_CUSTOM_REAL
          tempz3lw = 0.0_CUSTOM_REAL

          ! we can merge these loops because NGLLX = NGLLY = NGLLZ
          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1ls = tempx1ls + tempx1(l,j,k)*fac1
            tempy1ls = tempy1ls + tempy1(l,j,k)*fac1
            tempz1ls = tempz1ls + tempz1(l,j,k)*fac1
            tempx1lw = tempx1lw + tempx1p(l,j,k)*fac1
            tempy1lw = tempy1lw + tempy1p(l,j,k)*fac1
            tempz1lw = tempz1lw + tempz1p(l,j,k)*fac1

            fac2 = hprimewgll_yy(l,j)
            tempx2ls = tempx2ls + tempx2(i,l,k)*fac2
            tempy2ls = tempy2ls + tempy2(i,l,k)*fac2
            tempz2ls = tempz2ls + tempz2(i,l,k)*fac2
            tempx2lw = tempx2lw + tempx2p(i,l,k)*fac2
            tempy2lw = tempy2lw + tempy2p(i,l,k)*fac2
            tempz2lw = tempz2lw + tempz2p(i,l,k)*fac2

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

          phil = phistore(i,j,k,ispec)
          tortl = tortstore(i,j,k,ispec)

          ! sum contributions from each element to the global mesh

          iglob = ibool(i,j,k,ispec)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) &
                                      - ( fac1*(tempx1ls - phil/tortl*tempx1lw) &
                                        + fac2*(tempx2ls - phil/tortl*tempx2lw) &
                                        + fac3*(tempx3ls - phil/tortl*tempx3lw) )

          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) &
                                      - ( fac1*(tempy1ls - phil/tortl*tempy1lw) &
                                        + fac2*(tempy2ls - phil/tortl*tempy2lw) &
                                        + fac3*(tempy3ls - phil/tortl*tempy3lw) )

          accels_poroelastic(3,iglob) = accels_poroelastic(3,iglob) &
                                      - ( fac1*(tempz1ls - phil/tortl*tempz1lw) &
                                        + fac2*(tempz2ls - phil/tortl*tempz2lw) &
                                        + fac3*(tempz3ls - phil/tortl*tempz3lw) )

          !
          !---- viscous damping
          !
          ! add + phi/tort eta_f k^-1 dot(w)

          etal_f = etastore(i,j,k,ispec)

          if (etal_f > 0.0_CUSTOM_REAL) then

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

            if (detk /= 0.d0) then
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

            ! if (VISCOATTENUATION) then
            !   bl_unrelaxed(1) = etal_f*invpermlxx*theta_e/theta_s
            !   bl_unrelaxed(2) = etal_f*invpermlxz*theta_e/theta_s
            !   bl_unrelaxed(3) = etal_f*invpermlzz*theta_e/theta_s
            ! endif

            ! do k = 1,NGLLZ
            !  do j = 1,NGLLY
            !    do i = 1,NGLLX

            !      iglob = ibool(i,j,k,ispec)

            !      if (VISCOATTENUATION) then
            !      ! compute the viscous damping term with the unrelaxed viscous coef and add memory variable
            !        viscodampx = velocw_poroelastic(1,iglob)*bl_unrelaxed(1) &
            !                   + velocw_poroelastic(2,iglob)*bl_unrelaxed(2) &
            !                   - rx_viscous(i,j,ispec)
            !        viscodampz = velocw_poroelastic(1,iglob)*bl_unrelaxed(2) &
            !                   + velocw_poroelastic(2,iglob)*bl_unrelaxed(3) &
            !                   - rz_viscous(i,j,ispec)
            !      else

            ! no viscous attenuation
            viscodampx = velocw_poroelastic(1,iglob)*bl_relaxed(1) &
                       + velocw_poroelastic(2,iglob)*bl_relaxed(2) &
                       + velocw_poroelastic(3,iglob)*bl_relaxed(3)
            viscodampy = velocw_poroelastic(1,iglob)*bl_relaxed(2) &
                       + velocw_poroelastic(2,iglob)*bl_relaxed(4) &
                       + velocw_poroelastic(3,iglob)*bl_relaxed(5)
            viscodampz = velocw_poroelastic(1,iglob)*bl_relaxed(3) &
                       + velocw_poroelastic(2,iglob)*bl_relaxed(5) &
                       + velocw_poroelastic(3,iglob)*bl_relaxed(6)
            !      endif

            if (ispec_irreg /= 0) jacobianl = jacobianstore(i,j,k,ispec_irreg)

            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) &
                                        + phil/tortl*wxgll(i)*wygll(j)*wzgll(k)*jacobianl*viscodampx
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) &
                                        + phil/tortl*wxgll(i)*wygll(j)*wzgll(k)*jacobianl*viscodampy
            accels_poroelastic(3,iglob) = accels_poroelastic(3,iglob) &
                                        + phil/tortl*wxgll(i)*wygll(j)*wzgll(k)*jacobianl*viscodampz

            !      ! if isolver == 1 .and. save_forward then b_viscodamp is saved in compute_forces_poro_fluid_part.f90
            !      if (isolver == 2) then ! kernels calculation
            !        b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + phil/tortl*b_viscodampx(iglob)
            !        b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + phil/tortl*b_viscodampz(iglob)
            !      endif

            !      enddo
            !    enddo
            !  enddo

          endif ! if (etal_f > 0.0) then

        enddo ! second loop over the GLL points
      enddo
    enddo

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_poro_solid_part

