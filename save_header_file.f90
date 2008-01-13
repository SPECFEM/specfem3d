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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC,NPROC_XI,NPROC_ETA,&
             UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,ATTENUATION,ANISOTROPY,NSTEP, &
             NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
             NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,SIMULATION_TYPE)

  implicit none

  include "constants.h"

  integer NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC,NSTEP,NPROC_XI,NPROC_ETA, &
             NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
             NPOIN2DMAX_XY,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,SIMULATION_TYPE

  logical ATTENUATION,ANISOTROPY

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX

  character(len=150) HEADER_FILE

! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')

! define maximum size for message buffers
  NPOIN2DMAX_XY = max(NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX)

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

  write(IOUT,*) 'integer, parameter ::  NPROC_VAL = ',NPROC
  write(IOUT,*) 'integer, parameter :: NPROC_XI_VAL = ', NPROC_XI
  write(IOUT,*) 'integer, parameter :: NPROC_ETA_VAL = ', NPROC_ETA

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
  write(IOUT,*) 'integer, parameter :: NSPEC_AB_VAL = ',NSPEC_AB
  write(IOUT,*) 'integer, parameter :: NGLOB_AB_VAL = ',NGLOB_AB
  write(IOUT,*)
  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'

! if attenuation is off, set dummy size of arrays to one
  if(ATTENUATION) then
    write(IOUT,*) 'integer, parameter :: NSPEC_ATTENUATION = ', NSPEC_AB
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_ATTENUATION = ', 1
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .false.'
  endif
  write(IOUT,*)

! anisotropy
  if(ANISOTROPY) then
    write(IOUT,*) 'integer, parameter :: NSPEC_ANISO = ',NSPEC_AB
    write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .true.'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_ANISO = ', 1
    write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .false.'
  endif

  write(IOUT,*)

! strain/attenuation
  if (ATTENUATION .and. SIMULATION_TYPE == 3) then
    write(IOUT,*) 'integer, parameter :: NSPEC_ATT_AND_KERNEL = ', NSPEC_AB
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_ATT_AND_KERNEL = ', 1
  endif

  ! adjoint
  if (SIMULATION_TYPE == 3) then
    write(IOUT,*) 'integer, parameter :: NSPEC_ADJOINT = ', NSPEC_AB
    write(IOUT,*) 'integer, parameter :: NGLOB_ADJOINT = ', NGLOB_AB
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_ADJOINT = ', 1
    write(IOUT,*) 'integer, parameter :: NGLOB_ADJOINT = ', 1
  endif

  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_XMIN_XMAX_VAL = ', NSPEC2DMAX_XMIN_XMAX
  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_YMIN_YMAX_VAL = ', NSPEC2DMAX_YMIN_YMAX
  write(IOUT,*) 'integer, parameter :: NSPEC2D_BOTTOM_VAL = ', NSPEC2D_BOTTOM
  write(IOUT,*) 'integer, parameter :: NSPEC2D_TOP_VAL = ', NSPEC2D_TOP
  write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_XMIN_XMAX_VAL = ', NPOIN2DMAX_XMIN_XMAX
  write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_YMIN_YMAX_VAL = ', NPOIN2DMAX_YMIN_YMAX
  write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_XY_VAL = ', NPOIN2DMAX_XY

  write(IOUT,*)

! Moho boundary
  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
    write(IOUT,*) 'integer, parameter :: NSPEC2D_MOHO_BOUN = ', NSPEC2D_BOTTOM
    write(IOUT,*) 'integer, parameter :: NSPEC_BOUN = ', NSPEC_AB
  else
    write(IOUT,*) 'integer, parameter :: NSPEC2D_MOHO_BOUN = ', 1
    write(IOUT,*) 'integer, parameter :: NSPEC_BOUN = ', 1
  endif

  close(IOUT)

  end subroutine save_header_file

