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

  subroutine save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
             ATTENUATION,ANISOTROPY,NSTEP,DT, &
             SIMULATION_TYPE,static_memory_size,nfaces_surface_glob_ext_mesh)

  implicit none

  include "constants.h"

! number of points per surface element
  integer, parameter :: NGLLSQUARE_NDIM = NGLLSQUARE * NDIM

  integer NSPEC_AB,NGLOB_AB,NPROC,NSTEP,SIMULATION_TYPE
           !  NPOIN2DMAX_XY,NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,

  logical ATTENUATION,ANISOTROPY

  double precision DT

  double precision :: static_memory_size

  character(len=256) HEADER_FILE

  integer :: nfaces_surface_glob_ext_mesh
  
! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')

! define maximum size for message buffers
  !NPOIN2DMAX_XY = max(NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX)

  open(unit=IOUT,file=HEADER_FILE,status='unknown')
  write(IOUT,*)

  write(IOUT,*) '!'
  write(IOUT,*) '! this is the parameter file for static compilation of the solver'
  write(IOUT,*) '!'
  write(IOUT,*) '! mesh statistics:'
  write(IOUT,*) '! ---------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK these statistics are now INCORRECT'
  write(IOUT,*) '! DK DK because the CUBIT + SCOTCH mesh has'
  write(IOUT,*) '! DK DK a different number of mesh elements and points in each slice'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '! DK DK'
  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROC
  write(IOUT,*) '!'
  write(IOUT,*) '! number of ES nodes = ',real(NPROC)/8.
  write(IOUT,*) '! percentage of total 640 ES nodes = ',100.*(real(NPROC)/8.)/640.,' %'
  write(IOUT,*) '! total memory available on these ES nodes (Gb) = ',16.*real(NPROC)/8.

! write(IOUT,*) 'integer, parameter ::  NPROC_VAL = ',NPROC
! write(IOUT,*) 'integer, parameter :: NPROC_XI_VAL = ', NPROC_XI
! write(IOUT,*) 'integer, parameter :: NPROC_ETA_VAL = ', NPROC_ETA

  write(IOUT,*) '!'
!  write(IOUT,*) '! max points per processor = max vector length = ',NGLOB_AB
  write(IOUT,*) '! min vector length = ',NGLLSQUARE
  write(IOUT,*) '! min critical vector length = ',NGLLSQUARE_NDIM
  write(IOUT,*) '!'
!  write(IOUT,*) '! on ES and SX-5, make sure "loopcnt=" parameter'
!  write(IOUT,*) '! in Makefile is greater than ',NGLOB_AB
!  write(IOUT,*) '!'

!  write(IOUT,*) '! total elements per AB slice = ',NSPEC_AB
!  write(IOUT,*) '! total points per AB slice = ',NGLOB_AB
  write(IOUT,*) '! not valid for external mesh files: total points per AB slice = ',NGLOB_AB
  write(IOUT,*) '! total elements per AB slice = (will be read in external file)'
  write(IOUT,*) '! total points per AB slice = (will be read in external file)'
  write(IOUT,*) '!'

  write(IOUT,*) '! total for full mesh:'
  write(IOUT,*) '! -------------------'
  write(IOUT,*) '!'
!  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
!  write(IOUT,*) '! ',NPROC*NSPEC_AB
!  write(IOUT,*) '! approximate total number of points in entire mesh = '
!  write(IOUT,*) '! ',dble(NPROC)*dble(NGLOB_AB)
! there are 3 DOFs in solid regions
!  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
!  write(IOUT,*) '! ',3.d0*dble(NPROC)*dble(NGLOB_AB)
!  write(IOUT,*) '!'

  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'
  write(IOUT,*) '! time step = ',DT
  write(IOUT,*) '!'

! if attenuation is off, set dummy size of arrays to one
! both parameters are obsolete for specfem3D
! they are only used in ampuero_implicit_ABC_specfem3D.f90
  write(IOUT,*) '! only needed for ampuero_implicit_ABC_specfem3D.f90 compilation: '
  write(IOUT,*) '! (uncomment next line) '
  if(ATTENUATION) then
    write(IOUT,*) '! integer, parameter :: NSPEC_ATTENUATION = ', NSPEC_AB
!    write(IOUT,*) '! logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) '! integer, parameter :: NSPEC_ATTENUATION = ', 1
!    write(IOUT,*) '! logical, parameter :: ATTENUATION_VAL = .false.'
  endif
  
  write(IOUT,*)

! anisotropy
  if(ANISOTROPY) then
    !stop 'ANISOTROPY not supported yet in the CUBIT + SCOTCH version because of arrays of constant size defined'
    !write(IOUT,*) 'integer, parameter :: NSPEC_ANISO = ',NSPEC_AB
    !write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .true.'
    write(IOUT,*) '! with anisotropy'
  else
    !write(IOUT,*) 'integer, parameter :: NSPEC_ANISO = ', 1
    !write(IOUT,*) 'logical, parameter :: ANISOTROPY_VAL = .false.'
    write(IOUT,*) '! no anisotropy'
  endif

  write(IOUT,*)
    
!! DK DK May 2009: removed all the things that are not supported in the CUBIT + SCOTCH version yet
!! DK DK May 2009: removed all the things that are not supported in the CUBIT + SCOTCH version yet
!! DK DK May 2009: removed all the things that are not supported in the CUBIT + SCOTCH version yet

  write(IOUT,*) '! approximate static memory needed by the solver:'
  write(IOUT,*) '! ----------------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! size of static arrays for the biggest slice = ',static_memory_size/1048576.d0,' MB'
  write(IOUT,*) '!                                             = ',static_memory_size/1073741824.d0,' GB'
  write(IOUT,*) '!'
  write(IOUT,*) '!   (should be below and typically equal to 80% of 1.5 GB = 1.2 GB on pangu'
  write(IOUT,*) '!    at Caltech, and below and typically equal to 85% of 2 GB = 1.7 GB'
  write(IOUT,*) '!    on Marenostrum in Barcelona)'
  write(IOUT,*) '!   (if significantly more, the job will not run by lack of memory)'
  write(IOUT,*) '!   (if significantly less, you waste a significant amount of memory)'
  write(IOUT,*) '!'

! strain/attenuation
  if (ATTENUATION .and. SIMULATION_TYPE == 3) then
!   write(IOUT,*) 'integer, parameter :: NSPEC_ATT_AND_KERNEL = ', NSPEC_AB
  else
!   write(IOUT,*) 'integer, parameter :: NSPEC_ATT_AND_KERNEL = ', 1
  endif

  ! adjoint
  if (SIMULATION_TYPE == 3) then
!   write(IOUT,*) 'integer, parameter :: NSPEC_ADJOINT = ', NSPEC_AB
!   write(IOUT,*) 'integer, parameter :: NGLOB_ADJOINT = ', NGLOB_AB
  else
!   write(IOUT,*) 'integer, parameter :: NSPEC_ADJOINT = ', 1
!   write(IOUT,*) 'integer, parameter :: NGLOB_ADJOINT = ', 1
  endif

  write(IOUT,*)

! write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_XMIN_XMAX_VAL = ', NSPEC2DMAX_XMIN_XMAX
! write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_YMIN_YMAX_VAL = ', NSPEC2DMAX_YMIN_YMAX
! write(IOUT,*) 'integer, parameter :: NSPEC2D_BOTTOM_VAL = ', NSPEC2D_BOTTOM
! write(IOUT,*) 'integer, parameter :: NSPEC2D_TOP_VAL = ', NSPEC2D_TOP
! write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_XMIN_XMAX_VAL = ', NPOIN2DMAX_XMIN_XMAX
! write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_YMIN_YMAX_VAL = ', NPOIN2DMAX_YMIN_YMAX
! write(IOUT,*) 'integer, parameter :: NPOIN2DMAX_XY_VAL = ', NPOIN2DMAX_XY

  write(IOUT,*)

! Moho boundary
!  if (SAVE_MOHO_MESH .and. SIMULATION_TYPE == 3) then
!   write(IOUT,*) 'integer, parameter :: NSPEC2D_MOHO_BOUN = ', NSPEC2D_BOTTOM
!   write(IOUT,*) 'integer, parameter :: NSPEC_BOUN = ', NSPEC_AB
!  else
!   write(IOUT,*) 'integer, parameter :: NSPEC2D_MOHO_BOUN = ', 1
!   write(IOUT,*) 'integer, parameter :: NSPEC_BOUN = ', 1
!  endif

  close(IOUT)


! copy number of surface elements in an include file for the movies
  if( nfaces_surface_glob_ext_mesh > 0 ) then

    call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/surface_from_mesher.h')

    open(unit=IOUT,file=HEADER_FILE,status='unknown')
    write(IOUT,*) '!'
    write(IOUT,*) '! this is the parameter file for static compilation for movie creation'
    write(IOUT,*) '!'
    write(IOUT,*) '! number of elements containing surface faces '
    write(IOUT,*) '! ---------------'
    write(IOUT,*)    
    write(IOUT,*) 'integer,parameter :: NSPEC_SURFACE_EXT_MESH = ',nfaces_surface_glob_ext_mesh
    write(IOUT,*)
    close(IOUT)
    
  endif

  end subroutine save_header_file

