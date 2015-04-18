!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
! with this program; if not, read to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!! DK DK read binary unformatted arrays in Fortran format because they are faster
!! DK DK and smaller to use in the mesher, which is written in Fortran.

!! DK DK permutation of indices from Fortran order to C order is done automatically
!! DK DK (more precisely there is nothing to do) because these arrays are declared
!! DK DK with an inverted order of indices in the main calling program in C

  subroutine read_arrays_solver(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    kappav,muv,ibool,rmass,myrank,xstore,ystore,zstore)

  implicit none

  include "DATABASES_FOR_SOLVER/values_from_mesher_f90.h"

  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = 5
  integer, parameter :: NGLLZ = 5

  integer, parameter :: NGLL2 = NGLLX * NGLLX
  integer, parameter :: NGLL3 = NGLLX * NGLLX * NGLLX ! do not align for the regular C version for CPUs

! use single precision
  integer, parameter :: CUSTOM_REAL = 4

! only one region: crust_mantle
  integer, parameter :: iregion_code = 1

! unit to access the file
  integer, parameter :: IOUT = 56

  integer :: myrank,i,j,k,ispec

! arrays with jacobian matrix
  real(kind=CUSTOM_REAL), dimension(NGLL3*NSPEC) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,kappav,muv

  integer ibool(NGLL3*NSPEC)

! temporary arrays to convert from four indices: NGLLX,NGLLY,NGLLZ,NSPEC
! to only one index: NGLL3*NSPEC, and with padding to align memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: dummy_float_read
  integer dummy_int_read(NGLLX,NGLLY,NGLLZ,NSPEC)
! we can use the same temporary memory space to save space
  equivalence(dummy_float_read, dummy_int_read)

! mass matrix
  real(kind=CUSTOM_REAL) rmass(nglob)

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! processor identification
  character(len=250) prname

! create the name for the database of the current slide and region
  write(prname,"('DATABASES_FOR_SOLVER/proc',i6.6,'_reg',i1,'_')") myrank,iregion_code

  open(unit=IOUT,file=prname(1:len_trim(prname))//'database.dat',status='old',action='read',form='unformatted')

! read real numbers here
  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
! indices start at 1 in Fortran and 0 in C therefore we subtract 1 to all
! the indices in this formula to convert to a linear offset
          xix((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xiy((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xiz((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          etax((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          etay((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          etaz((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          gammax((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          gammay((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          gammaz((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          kappav((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  read(IOUT) dummy_float_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          muv((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_float_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

! read an integer here
  read(IOUT) dummy_int_read
  do ispec = 1,NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ibool((ispec-1)*NGLL3+(k-1)*NGLL2+(j-1)*NGLLX+i-1+1) = dummy_int_read(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

! read the mass matrix
  read(IOUT) rmass

! read the coordinates of the mesh
  read(IOUT) xstore
  read(IOUT) ystore
  read(IOUT) zstore

  close(IOUT)

  end subroutine read_arrays_solver

