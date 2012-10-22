!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            November 2010
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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

  subroutine create_header_file

  implicit none

  include "constants.h"

  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY

! parameters to be computed based upon parameters above read from file
  integer NPROC

  integer NSPEC_AB, NGLOB_AB
   !   NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,

  double precision DT,HDUR_MOVIE

  logical ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=256) LOCAL_PATH,HEADER_FILE

! ************** PROGRAM STARTS HERE **************

  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', &
       OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH))//'/values_from_mesher.h')

  print *
  print *,'creating file ', trim(HEADER_FILE), ' to compile solver with correct values'

! read the parameter file
  call read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY)

! create include file for the solver
  call save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
             ATTENUATION,ANISOTROPY,NSTEP,DT, &
             SIMULATION_TYPE,0.d0,0)
  print *
  print *,'edit file OUTPUT_FILES/values_from_mesher.h to see some statistics about the mesh'
  print *
!! DK DK May 2009: removed this because now each slice of a CUBIT + SCOTCH mesh
!! DK DK May 2009: has a different number of spectral elements and therefore the
!! DK DK May 2009: value below should be the max() for all the slices
! print *,'on NEC SX, make sure "loopcnt=" parameter'
! print *,'in Makefile is greater than max vector length = ',NGLOB_AB

  print *
  print *,'done'
  print *

  end subroutine create_header_file

