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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
             ATTENUATION,ANISOTROPY,NSTEP,DT,STACEY_INSTEAD_OF_FREE_SURFACE, &
             SIMULATION_TYPE,memory_size,nfaces_surface_glob_ext_mesh)

  implicit none

  include "constants.h"

  integer NSPEC_AB,NGLOB_AB,NPROC,NSTEP,SIMULATION_TYPE

  logical ATTENUATION,ANISOTROPY
  logical STACEY_INSTEAD_OF_FREE_SURFACE, ABSORB_FREE_SURFACE_VAL

  double precision DT, memory_size

  character(len=256) HEADER_FILE

  integer nfaces_surface_glob_ext_mesh

  NAMELIST/MESHER/ABSORB_FREE_SURFACE_VAL

  if (STACEY_INSTEAD_OF_FREE_SURFACE) then
      ABSORB_FREE_SURFACE_VAL = .true.
  else
      ABSORB_FREE_SURFACE_VAL = .false.
  endif

! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', &
       OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))//'/values_from_mesher.h')

  open(unit=IOUT,file=HEADER_FILE,status='unknown')
  write(IOUT,*)
  write(IOUT,*) '!'
  write(IOUT,*) '! purely informative use'
  write(IOUT,*) '!'
  write(IOUT,*) '! mesh statistics:'
  write(IOUT,*) '! ---------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! note: '
  write(IOUT,*) '!    the values are only approximate and differ for different processes'
  write(IOUT,*) '!    because the CUBIT + SCOTCH mesh has'
  write(IOUT,*) '!    a different number of mesh elements and points in each slice'
  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROC
  write(IOUT,*) '!'
  write(IOUT,*) '! number of ES nodes = ',real(NPROC)/8.
  write(IOUT,*) '! percentage of total 640 ES nodes = ',100.*(real(NPROC)/8.)/640.,' %'
  write(IOUT,*) '! total memory available on these ES nodes (Gb) = ',16.*real(NPROC)/8.
  write(IOUT,*) '!'
  write(IOUT,*) '! min vector length = ',NGLLSQUARE
  write(IOUT,*) '! min critical vector length = ',NGLLSQUARE * NDIM
  write(IOUT,*) '!'
  write(IOUT,*) '! master process: total points per AB slice = ',NGLOB_AB
  write(IOUT,*) '! total elements per AB slice = (will be read in external file)'
  write(IOUT,*) '! total points per AB slice = (will be read in external file)'
  write(IOUT,*) '!'
  write(IOUT,*) '! total for full mesh:'
  write(IOUT,*) '! -------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'
  write(IOUT,*) '! time step = ',DT
  write(IOUT,*) '!'
  write(IOUT,*) '! attenuation uses:'
  if(ATTENUATION) then
    write(IOUT,*) '!  NSPEC_ATTENUATION = ', NSPEC_AB
  else
    write(IOUT,*) '!  NSPEC_ATTENUATION = ', 1
  endif
  write(IOUT,*) '! '
  write(IOUT,*) '! anisotropy uses:'
  if(ANISOTROPY) then
    write(IOUT,*) '!  NSPEC_ANISO = ',NSPEC_AB
  else
    write(IOUT,*) '!  NSPEC_ANISO = ', 1
  endif
  write(IOUT,*) '! '
  write(IOUT,*) '! adjoint uses:'
  if (SIMULATION_TYPE == 3) then
    write(IOUT,*) '!  NSPEC_ADJOINT = ', NSPEC_AB
    write(IOUT,*) '!  NGLOB_ADJOINT = ', NGLOB_AB
  else
    write(IOUT,*) '!  NSPEC_ADJOINT = ', 1
    write(IOUT,*) '!  NGLOB_ADJOINT = ', 1
  endif
  write(IOUT,*) '! '
  write(IOUT,*) '! approximate least memory needed by the solver:'
  write(IOUT,*) '! ----------------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! size of arrays for the largest slice = ',memory_size/1048576.d0,' MB'
  write(IOUT,*) '!                                      = ',memory_size/1073741824.d0,' GB'
  write(IOUT,*) '!'
  write(IOUT,*) '!   (should be below 90% or so of the amount of memory available per processor core'
  write(IOUT,*) '!   (if significantly more, the job will not run by lack of memory)'
  write(IOUT,*) '!   (if significantly less, you waste a significant amount of memory)'
  write(IOUT,*) '!'
  write(IOUT,*) '! check parameter to ensure the code has been compiled with the right values:'
  write(IOUT,NML=MESHER)
  write(IOUT,*)
  close(IOUT)

! copy number of surface elements in an include file for the movies
  if( nfaces_surface_glob_ext_mesh > 0 ) then

    call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', &
         OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))//'/surface_from_mesher.h')

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

