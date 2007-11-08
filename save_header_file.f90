!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC, &
             UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,ATTENUATION,ANISOTROPY,NSTEP)

  implicit none

  include "constants.h"

  integer NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC,NSTEP

  logical ATTENUATION,ANISOTROPY

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX

  character(len=150) HEADER_FILE

! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')
  open(unit=IOUT,file=HEADER_FILE,status='unknown')
  write(IOUT,*)

  write(IOUT,*) '!'
  write(IOUT,*) '! this is the parameter file for static compilation of the solver'
  write(IOUT,*) '!'
  write(IOUT,*) '! mesh statistics:'
  write(IOUT,*) '! ---------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROC
  write(IOUT,*) '!'
  write(IOUT,*) '! number of ES nodes = ',real(NPROC)/8.
  write(IOUT,*) '! percentage of total 640 ES nodes = ',100.*(real(NPROC)/8.)/640.,' %'
  write(IOUT,*) '! total memory available on these ES nodes (Gb) = ',16.*real(NPROC)/8.

  write(IOUT,*) '!'
  write(IOUT,*) '! max points per processor = max vector length = ',NGLOB_AB
  write(IOUT,*) '! min vector length = ',NGLLSQUARE
  write(IOUT,*) '! min critical vector length = ',NGLLSQUARE_NDIM
  write(IOUT,*) '!'
  write(IOUT,*) '! on ES and SX-5, make sure "loopcnt=" parameter'
  write(IOUT,*) '! in Makefile is greater than ',NGLOB_AB
  write(IOUT,*) '!'

  write(IOUT,*) '! total elements per AB slice = ',NSPEC_AB
  write(IOUT,*) '! total points per AB slice = ',NGLOB_AB
  write(IOUT,*) '!'

  write(IOUT,*) '! total for full mesh:'
  write(IOUT,*) '! -------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '! ',NPROC*NSPEC_AB
  write(IOUT,*) '! approximate total number of points in entire mesh = '
  write(IOUT,*) '! ',dble(NPROC)*dble(NGLOB_AB)
! there are 3 DOFs in solid regions
  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '! ',3.d0*dble(NPROC)*dble(NGLOB_AB)
  write(IOUT,*) '!'

  write(IOUT,*) '! resolution of the mesh at the surface:'
  write(IOUT,*) '! -------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! spectral elements along X = ',NEX_XI
  write(IOUT,*) '! spectral elements along Y = ',NEX_ETA
  write(IOUT,*) '! GLL points along X = ',NEX_XI*(NGLLX-1) + 1
  write(IOUT,*) '! GLL points along Y = ',NEX_ETA*(NGLLY-1) + 1
  write(IOUT,*) '! average distance between points along X in m = ',sngl(UTM_X_MAX-UTM_X_MIN)/real(NEX_XI*(NGLLX-1))
  write(IOUT,*) '! average distance between points along Y in m = ',sngl(UTM_Y_MAX-UTM_Y_MIN)/real(NEX_ETA*(NGLLY-1))
  write(IOUT,*) '!'
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPEC_AB = ',NSPEC_AB
  write(IOUT,*) 'integer, parameter :: NGLOB_AB = ',NGLOB_AB
  write(IOUT,*)
  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'

! if attenuation is off, set dummy size of arrays to one
  if(ATTENUATION) then
    write(IOUT,*) 'integer, parameter :: NSPEC_ATTENUATION = NSPEC_AB'
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_ATTENUATION = 1'
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .false.'
  endif

! anisotropy
  if(ANISOTROPY) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .false.'
  endif

  write(IOUT,*)

  close(IOUT)

  end subroutine save_header_file

