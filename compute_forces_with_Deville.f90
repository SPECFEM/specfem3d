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

subroutine compute_forces_with_Deville(NSPEC_AB,NGLOB_AB,displ,accel,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
     hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
     kappastore,mustore,jacobian,ibool,ispec_is_inner,phase_is_inner, &
     NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,xi_source,eta_source,gamma_source,nu_source,hdur,dt)

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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(3,3,NSOURCES) :: nu_source
  double precision, dimension(NSOURCES) :: hdur
  double precision :: dt

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer ispec,iglob
  integer i,j,k

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal

  integer :: isource
  double precision :: t0,f0

  do ispec = 1,NSPEC_AB

  if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

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

!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1
  call mxm_m1_m2_5points(hprime_xx,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1)

  do k = 1,NGLLX
    call mxm_m1_m1_5points(dummyx_loc(1,1,k),dummyy_loc(1,1,k),dummyz_loc(1,1,k), &
           hprime_xxT,tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k))
  enddo

  call mxm_m2_m1_5points(dummyx_loc,dummyy_loc,dummyz_loc,hprime_xxT,tempx3,tempy3,tempz3)

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

          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2.*mul

! compute stress sigma
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

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

!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1
  call mxm_m1_m2_5points(hprimewgll_xxT,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1)

  do k = 1,NGLLX
    call mxm_m1_m1_5points(tempx2(1,1,k),tempy2(1,1,k),tempz2(1,1,k), &
          hprimewgll_xx,newtempx2(1,1,k),newtempy2(1,1,k),newtempz2(1,1,k))
  enddo

  call mxm_m2_m1_5points(tempx3,tempy3,tempz3,hprimewgll_xx,newtempx3,newtempy3,newtempz3)

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

! sum contributions from each element to the global mesh using indirect addressing
          iglob = ibool(i,j,k,ispec)
          accel(1,iglob) = accel(1,iglob) - fac1*newtempx1(i,j,k) - fac2*newtempx2(i,j,k) - fac3*newtempx3(i,j,k)
          accel(2,iglob) = accel(2,iglob) - fac1*newtempy1(i,j,k) - fac2*newtempy2(i,j,k) - fac3*newtempy3(i,j,k)
          accel(3,iglob) = accel(3,iglob) - fac1*newtempz1(i,j,k) - fac2*newtempz2(i,j,k) - fac3*newtempz3(i,j,k)

        enddo
      enddo
    enddo

  endif ! if (ispec_is_inner(ispec) .eqv. phase_is_inner)

  enddo  ! spectral element loop

! adding source
  do isource = 1,NSOURCES

  if (ispec_is_inner(ispec_selected_source(isource)) .eqv. phase_is_inner) then

  if(USE_FORCE_POINT_SOURCE) then

!   add the source (only if this proc carries the source)
    if(myrank == islice_selected_source(isource)) then

      iglob = ibool(nint(xi_source(isource)), &
           nint(eta_source(isource)), &
           nint(gamma_source(isource)), &
           ispec_selected_source(isource))
      f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
      t0 = 1.2d0/f0

  if (it == 1 .and. myrank == 0) then
    print *,'using a source of dominant frequency ',f0
    print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
    print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
  endif

      ! we use nu_source(:,3) here because we want a source normal to the surface.
      ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
      !accel(:,iglob) = accel(:,iglob) + &
      !     sngl(nu_source(:,3,isource) * 10000000.d0 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
      !     exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))
    accel(:,iglob) = accel(:,iglob) + &
           sngl(nu_source(:,3,isource) * 1.d10 * (1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)) * &
           exp(-PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0)))

    endif
  endif

  endif

  enddo

end subroutine compute_forces_with_Deville


!! DK DK subroutines adapted from Deville, Fischer and Mund, High-order methods
!! DK DK for incompressible fluid flow, Cambridge University Press (2002),
!! DK DK pages 386 and 389 and Figure 8.3.1

  subroutine mxm_m1_m2_5points(A,B1,B2,B3,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(m1,NGLLX) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1,B2,B3
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1,C2,C3

  integer :: i,j

  do j=1,m2
    do i=1,m1

      C1(i,j) = A(i,1)*B1(1,j) + &
                A(i,2)*B1(2,j) + &
                A(i,3)*B1(3,j) + &
                A(i,4)*B1(4,j) + &
                A(i,5)*B1(5,j)

      C2(i,j) = A(i,1)*B2(1,j) + &
                A(i,2)*B2(2,j) + &
                A(i,3)*B2(3,j) + &
                A(i,4)*B2(4,j) + &
                A(i,5)*B2(5,j)

      C3(i,j) = A(i,1)*B3(1,j) + &
                A(i,2)*B3(2,j) + &
                A(i,3)*B3(3,j) + &
                A(i,4)*B3(4,j) + &
                A(i,5)*B3(5,j)

    enddo
  enddo

  end subroutine mxm_m1_m2_5points

!---------

  subroutine mxm_m1_m1_5points(A1,A2,A3,B,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(m1,NGLLX) :: A1,A2,A3
  real(kind=CUSTOM_REAL), dimension(NGLLX,m1) :: B
  real(kind=CUSTOM_REAL), dimension(m1,m1) :: C1,C2,C3

  integer :: i,j

  do j=1,m1
    do i=1,m1

      C1(i,j) = A1(i,1)*B(1,j) + &
                A1(i,2)*B(2,j) + &
                A1(i,3)*B(3,j) + &
                A1(i,4)*B(4,j) + &
                A1(i,5)*B(5,j)

      C2(i,j) = A2(i,1)*B(1,j) + &
                A2(i,2)*B(2,j) + &
                A2(i,3)*B(3,j) + &
                A2(i,4)*B(4,j) + &
                A2(i,5)*B(5,j)

      C3(i,j) = A3(i,1)*B(1,j) + &
                A3(i,2)*B(2,j) + &
                A3(i,3)*B(3,j) + &
                A3(i,4)*B(4,j) + &
                A3(i,5)*B(5,j)

    enddo
  enddo

  end subroutine mxm_m1_m1_5points

!---------

  subroutine mxm_m2_m1_5points(A1,A2,A3,B,C1,C2,C3)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1,A2,A3
  real(kind=CUSTOM_REAL), dimension(NGLLX,m1) :: B
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1,C2,C3

  integer :: i,j

  do j=1,m1
    do i=1,m2

      C1(i,j) = A1(i,1)*B(1,j) + &
                A1(i,2)*B(2,j) + &
                A1(i,3)*B(3,j) + &
                A1(i,4)*B(4,j) + &
                A1(i,5)*B(5,j)

      C2(i,j) = A2(i,1)*B(1,j) + &
                A2(i,2)*B(2,j) + &
                A2(i,3)*B(3,j) + &
                A2(i,4)*B(4,j) + &
                A2(i,5)*B(5,j)

      C3(i,j) = A3(i,1)*B(1,j) + &
                A3(i,2)*B(2,j) + &
                A3(i,3)*B(3,j) + &
                A3(i,4)*B(4,j) + &
                A3(i,5)*B(5,j)

    enddo
  enddo

  end subroutine mxm_m2_m1_5points

