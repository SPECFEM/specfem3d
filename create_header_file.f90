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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

  subroutine create_header_file

  implicit none

  include "constants.h"

  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
            NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

! parameters to be computed based upon parameters above read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER

  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,DT,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=150) LOCAL_PATH,MODEL,HEADER_FILE

! ************** PROGRAM STARTS HERE **************

  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')
  print *
  print *,'creating file ', trim(HEADER_FILE), ' to compile solver with correct values'

! read the parameter file
  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,USE_REGULAR_MESH)

! create include file for the solver
  call save_header_file(NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC,NPROC_XI,NPROC_ETA, &
             UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,ATTENUATION,ANISOTROPY,NSTEP, &
             NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
             NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,SIMULATION_TYPE)
  print *
  print *,'edit file OUTPUT_FILES/values_from_mesher.h to see some statistics about the mesh'
  print *
  print *,'on NEC SX, make sure "loopcnt=" parameter'
  print *,'in Makefile is greater than max vector length = ',NGLOB_AB
  print *
  print *, '********************************'
  print *,'number of time steps = ',NSTEP
  print *,'time-stepping of the solver will be: ',DT
  print *, '********************************'
  print *
  print *,'done'
  print *

  end subroutine create_header_file

