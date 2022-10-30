!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  program serial_specfem3D

  implicit none

!!!!!!!!!!
!!!!!!!!!! NGLLX, NGLLY and NGLLZ are set equal to 5,
!!!!!!!!!! therefore each element contains NGLLX * NGLLY * NGLLZ = 125 points.
!!!!!!!!!!

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "DATABASES_FOR_SOLVER/values_from_mesher_f90.h"

! precision of the calculations (change to 8 for double precision)
  integer, parameter :: CUSTOM_REAL = 4

  real(kind=CUSTOM_REAL), parameter :: pi = 3.141592653589793_CUSTOM_REAL

  integer, parameter :: IIN = 40

! number of GLL integration points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! 3-D simulation
  integer, parameter :: NDIM = 3

! global displacement, velocity and acceleration vectors
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ

! global diagonal mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: rmass_fictitious_inverse,assembled_gradient_to_average

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,kappav,muv

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  integer, parameter :: myrank = 0

  integer :: ispec,iglob,i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl,hp1,hp2,hp3
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

! estimate of total memory size used
  print *
  print *,'NSPEC = ',NSPEC
  print *,'NGLOB = ',NGLOB
  print *

! read the mesh from external file
  call read_arrays_solver(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        kappav,muv,ibool,rmass_fictitious_inverse,myrank,xstore,ystore,zstore)

  open(unit=IIN,file='DATABASES_FOR_SOLVER/matrices.dat',status='old')
  do j=1,NGLLY
    do i=1,NGLLX
      read(IIN,*) hprime_xx(i,j)
      read(IIN,*) hprimewgll_xx(i,j)
      read(IIN,*) wgllwgll_yz(i,j)
      read(IIN,*) wgllwgll_xz(i,j)
      read(IIN,*) wgllwgll_xy(i,j)
    enddo
  enddo
  close(IIN)

! clear the fictitious mass matrix
  rmass_fictitious_inverse(:) = 0._CUSTOM_REAL

! put in that fictitious mass matrix the valence of each mesh point i.e. the number of elements between which this point is shared;
! this is computed automatically by adding 1 from each point, and thus at the end the sum of all the 1s is equal to the
! number of elements that have this point in common
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            rmass_fictitious_inverse(iglob) = rmass_fictitious_inverse(iglob) + 1._CUSTOM_REAL
        enddo
      enddo
    enddo
  enddo

  print *
  print *,'minimum and maximum valence in the mesh (the minimum should always be 1):'
  print *,minval(rmass_fictitious_inverse),maxval(rmass_fictitious_inverse)
  if (abs(minval(rmass_fictitious_inverse) - 1._CUSTOM_REAL) > 0.000001_CUSTOM_REAL) &
      stop 'error: the minimum valence in the mesh is not one!'
  print *

! invert that diagonal fictitious mass matrix once and for all
  do i = 1,NGLOB
    rmass_fictitious_inverse(i) = 1._CUSTOM_REAL / rmass_fictitious_inverse(i)
  enddo

! assign a fictitious value to the array for the test
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            displ(1,iglob) = cos(2.d0*PI*xstore(i,j,k,ispec))
            displ(2,iglob) = cos(PI*ystore(i,j,k,ispec))
            displ(3,iglob) = cos(0.5d0*PI*zstore(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo

! clear the gradient that we will need to average (this needs to be done every time you compute a new average)
  assembled_gradient_to_average(:) = 0._CUSTOM_REAL

! big loop over all the elements in the mesh to localize data
! from the global vectors to the local mesh
! using indirect addressing (contained in array ibool)
! and then to compute the elemental contribution
! to the acceleration vector of each element of the finite-element mesh
  do ispec = 1,NSPEC

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

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1l = tempx1l + dummyx_loc(l,j,k)*hp1
            tempy1l = tempy1l + dummyy_loc(l,j,k)*hp1
            tempz1l = tempz1l + dummyz_loc(l,j,k)*hp1

            hp2 = hprime_xx(j,l)
            tempx2l = tempx2l + dummyx_loc(i,l,k)*hp2
            tempy2l = tempy2l + dummyy_loc(i,l,k)*hp2
            tempz2l = tempz2l + dummyz_loc(i,l,k)*hp2

            hp3 = hprime_xx(k,l)
            tempx3l = tempx3l + dummyx_loc(i,j,l)*hp3
            tempy3l = tempy3l + dummyy_loc(i,j,l)*hp3
            tempz3l = tempz3l + dummyz_loc(i,j,l)*hp3
          enddo

!         compute derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl)- &
                xiyl*(etaxl*gammazl-etazl*gammaxl)+xizl*(etaxl*gammayl-etayl*gammaxl))

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! sum (i.e., assemble) the gradient that we will need to average
! here we use duxdxl as an example
          iglob = ibool(i,j,k,ispec)
          assembled_gradient_to_average(iglob) = assembled_gradient_to_average(iglob) + duxdxl

          enddo
        enddo
      enddo

  enddo   ! end of main loop on all the elements of the mesh

! compute the average by dividing the sum obtained by the sum obtained in the fictitious mass matrix,
! i.e. by the valence of each point in the mesh; since we have inverted the fictitious mass matrix
! once and for all above we multiply here instead of dividing (since on computers multiplying is always
! much faster than dividing)
  assembled_gradient_to_average(:) = assembled_gradient_to_average(:) * rmass_fictitious_inverse(:)

  print *,'Finished, array assembled_gradient_to_average now contains the average, which is continuous in the whole mesh'

  end program serial_specfem3D

