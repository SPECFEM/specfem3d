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

  subroutine compute_element_strain(ispec,nglob,displ,epsilondev_loc,eps_trace_over_3_loc)

! computes strain for single element

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,ONE_THIRD

  use specfem_par, only: ibool, hprime_xx, hprime_yy, hprime_zz, &
    xixstore, xiystore, xizstore, etaxstore, etaystore, etazstore, gammaxstore, gammaystore, gammazstore, &
    irregular_element_number, xix_regular

  implicit none

  ! element id
  integer,intent(in) :: ispec

  integer,intent(in) :: NGLOB

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ),intent(out) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: eps_trace_over_3_loc

  ! local parameters
  integer :: iglob
  integer :: i,j,k,l,ispec_irreg

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  ispec_irreg = irregular_element_number(ispec)

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l = 1,NGLLX
          hp1 = hprime_xx(i,l)
          iglob = ibool(l,j,k,ispec)
          tempx1l = tempx1l + displ(1,iglob)*hp1
          tempy1l = tempy1l + displ(2,iglob)*hp1
          tempz1l = tempz1l + displ(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
          hp2 = hprime_yy(j,l)
          iglob = ibool(i,l,k,ispec)
          tempx2l = tempx2l + displ(1,iglob)*hp2
          tempy2l = tempy2l + displ(2,iglob)*hp2
          tempz2l = tempz2l + displ(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
          hp3 = hprime_zz(k,l)
          iglob = ibool(i,j,l,ispec)
          tempx3l = tempx3l + displ(1,iglob)*hp3
          tempy3l = tempy3l + displ(2,iglob)*hp3
          tempz3l = tempz3l + displ(3,iglob)*hp3
        enddo

        ! get derivatives of ux, uy and uz with respect to x, y and z
        if (ispec_irreg /= 0) then
          ! irregular element
          xixl = xixstore(i,j,k,ispec_irreg)
          xiyl = xiystore(i,j,k,ispec_irreg)
          xizl = xizstore(i,j,k,ispec_irreg)
          etaxl = etaxstore(i,j,k,ispec_irreg)
          etayl = etaystore(i,j,k,ispec_irreg)
          etazl = etazstore(i,j,k,ispec_irreg)
          gammaxl = gammaxstore(i,j,k,ispec_irreg)
          gammayl = gammaystore(i,j,k,ispec_irreg)
          gammazl = gammazstore(i,j,k,ispec_irreg)

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l
        else
          ! regular element
          duxdxl = xix_regular*tempx1l
          duxdyl = xix_regular*tempx2l
          duxdzl = xix_regular*tempx3l

          duydxl = xix_regular*tempy1l
          duydyl = xix_regular*tempy2l
          duydzl = xix_regular*tempy3l

          duzdxl = xix_regular*tempz1l
          duzdyl = xix_regular*tempz2l
          duzdzl = xix_regular*tempz3l
        endif

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        ! compute deviatoric strain
        eps_trace_over_3_loc(i,j,k) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        epsilondev_loc(1,i,j,k) = duxdxl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(2,i,j,k) = duydyl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_strain
