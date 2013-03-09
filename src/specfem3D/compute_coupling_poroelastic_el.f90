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

! for poroelastic solver

  subroutine compute_coupling_poroelastic_el(NSPEC_AB,NGLOB_AB,ibool,&
                        displs_poroelastic,accels_poroelastic,displw_poroelastic,&
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz,&
                        kappaarraystore,rhoarraystore,mustore, &
                        phistore,tortstore,jacobian,&
                        displ,kappastore, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store,&
                        c22store,c23store,c24store,c25store,c26store,c33store,&
                        c34store,c35store,c36store,c44store,c45store,c46store,&
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                        num_coupling_el_po_faces, &
                        coupling_el_po_ispec,coupling_po_el_ispec, &
                        coupling_el_po_ijk,coupling_po_el_ijk, &
                        coupling_el_po_normal, &
                        coupling_el_po_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! returns the updated accelerations array: accels_poroelatsic (& accelw_poroelastic )

  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

! displacements, etc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displs_poroelastic,accels_poroelastic,&
                                                      displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  integer :: SIMULATION_TYPE
  integer :: NGLOB_ADJOINT,NSPEC_ADJOINT

! elastic-poroelastic coupling surface
  integer :: num_coupling_el_po_faces
  real(kind=CUSTOM_REAL) :: coupling_el_po_normal(NDIM,NGLLSQUARE,num_coupling_el_po_faces)
  real(kind=CUSTOM_REAL) :: coupling_el_po_jacobian2Dw(NGLLSQUARE,num_coupling_el_po_faces)
  integer :: coupling_el_po_ijk(3,NGLLSQUARE,num_coupling_el_po_faces)
  integer :: coupling_po_el_ijk(3,NGLLSQUARE,num_coupling_el_po_faces)
  integer :: coupling_el_po_ispec(num_coupling_el_po_faces)
  integer :: coupling_po_el_ispec(num_coupling_el_po_faces)

! properties
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        phistore,tortstore,jacobian
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappaarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        kappastore,mustore

! anisotropy
  logical :: ANISOTROPY
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO) :: &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
            c55store,c56store,c66store

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! local parameters
  real(kind=CUSTOM_REAL) :: sigmap,sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,phil,tortl,rhol_bar
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL) :: kappal_s
  real(kind=CUSTOM_REAL) :: kappal_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot
  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G

! local anisotropy parameters
  real(kind=CUSTOM_REAL) c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  integer :: iface,igll,ispec_po,ispec_el,iglob,iglob_po,iglob_el
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL) tempx1ls,tempx2ls,tempx3ls,tempx1lw,tempx2lw,tempx3lw
  real(kind=CUSTOM_REAL) tempy1ls,tempy2ls,tempy3ls,tempy1lw,tempy2lw,tempy3lw
  real(kind=CUSTOM_REAL) tempz1ls,tempz2ls,tempz3ls,tempz1lw,tempz2lw,tempz3lw

  real(kind=CUSTOM_REAL) :: duxdxl,duydxl,duzdxl,duxdyl,duydyl,duzdyl,duxdzl,duydzl,duzdzl
  real(kind=CUSTOM_REAL) :: dwxdxl,dwydxl,dwzdxl,dwxdyl,dwydyl,dwzdyl,dwxdzl,dwydzl,dwzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl_plus_duzdzl,dwxdxl_plus_dwydyl_plus_dwzdzl

  real(kind=CUSTOM_REAL) hp1,hp2,hp3

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

! loops on all coupling faces
  do iface = 1,num_coupling_el_po_faces

    ! gets corresponding poro/elastic spectral element
    ispec_po = coupling_el_po_ispec(iface)
    ispec_el = coupling_po_el_ispec(iface)

    if( ispec_is_inner(ispec_po) .eqv. phase_is_inner ) then

      ! loops over common GLL points
      do igll = 1, NGLLSQUARE

    !-----------------------
    ! from the elastic side
    !-----------------------
        i = coupling_po_el_ijk(1,igll,iface)
        j = coupling_po_el_ijk(2,igll,iface)
        k = coupling_po_el_ijk(3,igll,iface)

        ! gets global index of this common GLL point
        ! (note: should be the same as for corresponding i',j',k',ispec_poroelastic or ispec_elastic )
        iglob_el = ibool(i,j,k,ispec_el)

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
            iglob = ibool(l,j,k,ispec_el)
            tempx1l = tempx1l + displ(1,iglob)*hp1
            tempy1l = tempy1l + displ(2,iglob)*hp1
            tempz1l = tempz1l + displ(3,iglob)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec_el)
            tempx2l = tempx2l + displ(1,iglob)*hp2
            tempy2l = tempy2l + displ(2,iglob)*hp2
            tempz2l = tempz2l + displ(3,iglob)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec_el)
            tempx3l = tempx3l + displ(1,iglob)*hp3
            tempy3l = tempy3l + displ(2,iglob)*hp3
            tempz3l = tempz3l + displ(3,iglob)*hp3
          enddo

          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec_el)
          xiyl = xiy(i,j,k,ispec_el)
          xizl = xiz(i,j,k,ispec_el)
          etaxl = etax(i,j,k,ispec_el)
          etayl = etay(i,j,k,ispec_el)
          etazl = etaz(i,j,k,ispec_el)
          gammaxl = gammax(i,j,k,ispec_el)
          gammayl = gammay(i,j,k,ispec_el)
          gammazl = gammaz(i,j,k,ispec_el)
          jacobianl = jacobian(i,j,k,ispec_el)

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          ! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

          kappal = kappastore(i,j,k,ispec_el)
          mul = mustore(i,j,k,ispec_el)

          ! full anisotropic case, stress calculations
          if(ANISOTROPY) then
            c11 = c11store(i,j,k,ispec_el)
            c12 = c12store(i,j,k,ispec_el)
            c13 = c13store(i,j,k,ispec_el)
            c14 = c14store(i,j,k,ispec_el)
            c15 = c15store(i,j,k,ispec_el)
            c16 = c16store(i,j,k,ispec_el)
            c22 = c22store(i,j,k,ispec_el)
            c23 = c23store(i,j,k,ispec_el)
            c24 = c24store(i,j,k,ispec_el)
            c25 = c25store(i,j,k,ispec_el)
            c26 = c26store(i,j,k,ispec_el)
            c33 = c33store(i,j,k,ispec_el)
            c34 = c34store(i,j,k,ispec_el)
            c35 = c35store(i,j,k,ispec_el)
            c36 = c36store(i,j,k,ispec_el)
            c44 = c44store(i,j,k,ispec_el)
            c45 = c45store(i,j,k,ispec_el)
            c46 = c46store(i,j,k,ispec_el)
            c55 = c55store(i,j,k,ispec_el)
            c56 = c56store(i,j,k,ispec_el)
            c66 = c66store(i,j,k,ispec_el)

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


    !-----------------------
    ! from the poroelastic side
    !-----------------------
        i = coupling_el_po_ijk(1,igll,iface)
        j = coupling_el_po_ijk(2,igll,iface)
        k = coupling_el_po_ijk(3,igll,iface)

        iglob_po = ibool(i,j,k,ispec_po)
        if (iglob_el /= iglob_po) stop 'poroelastic-elastic coupling error'

! get poroelastic parameters of current local GLL
    phil = phistore(i,j,k,ispec_po)
    tortl = tortstore(i,j,k,ispec_po)
!solid properties
    kappal_s = kappaarraystore(1,i,j,k,ispec_po)
    rhol_s = rhoarraystore(1,i,j,k,ispec_po)
!fluid properties
    kappal_f = kappaarraystore(2,i,j,k,ispec_po)
    rhol_f = rhoarraystore(2,i,j,k,ispec_po)
!frame properties
    mul_fr = mustore(i,j,k,ispec_po)
    kappal_fr = kappaarraystore(3,i,j,k,ispec_po)
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)

!T = G:grad u_s + C_biot div w I
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
                iglob = ibool(l,j,k,ispec_po)
                tempx1ls = tempx1ls + displs_poroelastic(1,iglob)*hp1
                tempy1ls = tempy1ls + displs_poroelastic(2,iglob)*hp1
                tempz1ls = tempz1ls + displs_poroelastic(3,iglob)*hp1
                tempx1lw = tempx1lw + displw_poroelastic(1,iglob)*hp1
                tempy1lw = tempy1lw + displw_poroelastic(2,iglob)*hp1
                tempz1lw = tempz1lw + displw_poroelastic(3,iglob)*hp1
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                  ! to do
                  stop 'compute_coupling_poroelastic_el() : adjoint run not implemented yet'

                  ! dummy to avoid compiler warnings
                  iglob = NGLOB_ADJOINT
                  iglob = NSPEC_ADJOINT
                endif ! adjoint

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do
    !l=1,NGLLY
                hp2 = hprime_yy(j,l)
                iglob = ibool(i,l,k,ispec_po)
                tempx2ls = tempx2ls + displs_poroelastic(1,iglob)*hp2
                tempy2ls = tempy2ls + displs_poroelastic(2,iglob)*hp2
                tempz2ls = tempz2ls + displs_poroelastic(3,iglob)*hp2
                tempx2lw = tempx2lw + displw_poroelastic(1,iglob)*hp2
                tempy2lw = tempy2lw + displw_poroelastic(2,iglob)*hp2
                tempz2lw = tempz2lw + displw_poroelastic(3,iglob)*hp2
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                endif ! adjoint
    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

    !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do
    !l=1,NGLLZ
                hp3 = hprime_zz(k,l)
                iglob = ibool(i,j,l,ispec_po)
                tempx3ls = tempx3ls + displs_poroelastic(1,iglob)*hp3
                tempy3ls = tempy3ls + displs_poroelastic(2,iglob)*hp3
                tempz3ls = tempz3ls + displs_poroelastic(3,iglob)*hp3
                tempx3lw = tempx3lw + displw_poroelastic(1,iglob)*hp3
                tempy3lw = tempy3lw + displw_poroelastic(2,iglob)*hp3
                tempz3lw = tempz3lw + displw_poroelastic(3,iglob)*hp3
                ! adjoint simulations
                if (SIMULATION_TYPE == 3) then
                endif ! adjoint
          enddo

              xixl = xix(i,j,k,ispec_po)
              xiyl = xiy(i,j,k,ispec_po)
              xizl = xiz(i,j,k,ispec_po)
              etaxl = etax(i,j,k,ispec_po)
              etayl = etay(i,j,k,ispec_po)
              etazl = etaz(i,j,k,ispec_po)
              gammaxl = gammax(i,j,k,ispec_po)
              gammayl = gammay(i,j,k,ispec_po)
              gammazl = gammaz(i,j,k,ispec_po)
              jacobianl = jacobian(i,j,k,ispec_po)

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

! Dissipation only controlled by frame share attenuation in poroelastic (see
! Morency & Tromp, GJI 2008).
! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a
! linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611
! (1988).

!  else

! no attenuation
    sigma_xx = sigma_xx + &
            lambdalplus2mul_G*duxdxl + lambdal_G*duydyl_plus_duzdzl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl
    sigma_yy = sigma_yy + &
            lambdalplus2mul_G*duydyl + lambdal_G*duxdxl_plus_duzdzl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl
    sigma_zz = sigma_zz + &
            lambdalplus2mul_G*duzdzl + lambdal_G*duxdxl_plus_duydyl + C_biot*dwxdxl_plus_dwydyl_plus_dwzdzl

    sigma_xy = sigma_xy + mul_G*duxdyl_plus_duydxl
    sigma_xz = sigma_xz + mul_G*duzdxl_plus_duxdzl
    sigma_yz = sigma_yz + mul_G*duzdyl_plus_duydzl

    sigmap = C_biot*duxdxl_plus_duydyl_plus_duzdzl + M_biot*dwxdxl_plus_dwydyl_plus_dwzdzl


        ! gets associated normal on GLL point
        ! (note convention: pointing outwards of poroelastic element)
        nx = coupling_el_po_normal(1,igll,iface)
        ny = coupling_el_po_normal(2,igll,iface)
        nz = coupling_el_po_normal(3,igll,iface)

        ! gets associated, weighted 2D jacobian
        ! (note: should be the same for poroelastic and elastic element)
        jacobianw = coupling_el_po_jacobian2Dw(igll,iface)

        ! continuity of displacement and traction on global point
        !
        ! note: continuity of displacement is enforced after the velocity update
! contribution to the solid phase
        accels_poroelastic(1,iglob_po) = accels_poroelastic(1,iglob_po) + jacobianw*&
            ( (sigma_xx*nx + sigma_xy*ny + sigma_xz*nz)/2.d0 - phil/tortl*sigmap*nx )

        accels_poroelastic(2,iglob_po) = accels_poroelastic(2,iglob_po) + jacobianw*&
            ( (sigma_xy*nx + sigma_yy*ny + sigma_yz*nz)/2.d0 - phil/tortl*sigmap*ny )

        accels_poroelastic(3,iglob_po) = accels_poroelastic(3,iglob_po) + jacobianw*&
            ( (sigma_xz*nx + sigma_yz*ny + sigma_zz*nz)/2.d0 - phil/tortl*sigmap*nz )
! contribution to the fluid phase
! w = 0

      enddo ! igll

    endif

  enddo ! iface

end subroutine compute_coupling_poroelastic_el

