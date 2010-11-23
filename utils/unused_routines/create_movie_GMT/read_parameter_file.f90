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

  subroutine read_parameter_file(par_file,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO, &
        SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

  implicit none

  include "constants.h"

  character(len=*) :: par_file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
            NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, UTM_MAX
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,DT,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH

  character(len=256) LOCAL_PATH,MODEL,CMTSOLUTION

! local variables
  integer ios,icounter,isource,idummy,NEX_MAX

  double precision DEPTH_BLOCK_KM,RECORD_LENGTH_IN_SECONDS,hdur,minval_hdur

  character(len=256) dummystring

  integer, external :: err_occurred

  open(unit=IIN,file=trim(par_file),status='old',action='read')

  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) return

  call read_value_double_precision(LATITUDE_MIN, 'mesher.LATITUDE_MIN')
  if(err_occurred() /= 0) return
  call read_value_double_precision(LATITUDE_MAX, 'mesher.LATITUDE_MAX')
  if(err_occurred() /= 0) return
  call read_value_double_precision(LONGITUDE_MIN, 'mesher.LONGITUDE_MIN')
  if(err_occurred() /= 0) return
  call read_value_double_precision(LONGITUDE_MAX, 'mesher.LONGITUDE_MAX')
  if(err_occurred() /= 0) return
  call read_value_double_precision(DEPTH_BLOCK_KM, 'mesher.DEPTH_BLOCK_KM')
  if(err_occurred() /= 0) return
  call read_value_integer(UTM_PROJECTION_ZONE, 'mesher.UTM_PROJECTION_ZONE')
  if(err_occurred() /= 0) return
  call read_value_logical(SUPPRESS_UTM_PROJECTION, 'mesher.SUPPRESS_UTM_PROJECTION')
  if(err_occurred() /= 0) return

  call read_value_integer(NEX_XI, 'mesher.NEX_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(NEX_ETA, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) return
  call read_value_integer(NPROC_XI, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(NPROC_ETA, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) return

! convert model size to UTM coordinates and depth of mesh to meters
  call utm_geo(LONGITUDE_MIN,LATITUDE_MIN,UTM_X_MIN,UTM_Y_MIN,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
  call utm_geo(LONGITUDE_MAX,LATITUDE_MAX,UTM_X_MAX,UTM_Y_MAX,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

! check that parameters computed are consistent
  if(UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if(UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

! set time step and radial distribution of elements
! right distribution is determined based upon maximum value of NEX

  NEX_MAX = max(NEX_XI,NEX_ETA)
  UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM

  call read_value_string(MODEL, 'model.name')
  if(err_occurred() /= 0) return
  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) return
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) return
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_OLSEN_ATTENUATION, 'model.USE_OLSEN_ATTENUATION')
  if(err_occurred() /= 0) return

  if(dabs(DEPTH_BLOCK_KM) <= DEPTH_MOHO_SOCAL) &
    stop 'bottom of mesh must be deeper than deepest regional layer for Southern California'

! standard mesh for Southern California on Caltech cluster
  if (UTM_MAX/NEX_MAX >= 1.5) then

! time step in seconds
    DT                 = 0.011d0

! number of elements in the vertical direction
    NER_SEDIM          = 1
    NER_BASEMENT_SEDIM = 2
    NER_16_BASEMENT    = 2
    NER_MOHO_16        = 3
    NER_BOTTOM_MOHO    = 3

  else if(UTM_MAX/NEX_MAX >= 1.0) then

! time step in seconds
    DT                 = 0.009d0

! number of elements in the vertical direction
    NER_SEDIM          = 1
    NER_BASEMENT_SEDIM = 2
    NER_16_BASEMENT    = 3
    NER_MOHO_16        = 4
    NER_BOTTOM_MOHO    = 7

  else
    stop 'this value of NEX_MAX is not in the database, edit read_parameter_file.f90 and recompile'
  endif

! multiply parameters read by mesh scaling factor
  NER_BASEMENT_SEDIM = NER_BASEMENT_SEDIM * 2
  NER_16_BASEMENT = NER_16_BASEMENT * 2
  NER_MOHO_16 = NER_MOHO_16 * 2
  NER_BOTTOM_MOHO = NER_BOTTOM_MOHO * 4

  if(MODEL == 'SoCal') then

    BASEMENT_MAP             = .false.
    MOHO_MAP_LUPEI           = .false.
    HAUKSSON_REGIONAL_MODEL  = .false.
    HARVARD_3D_GOCAD_MODEL   = .false.
    THICKNESS_TAPER_BLOCK_HR = 1.d0
    THICKNESS_TAPER_BLOCK_MR = 1.d0
    IMPOSE_MINIMUM_VP_GOCAD  = .false.
    VP_MIN_GOCAD             = 1.d0
    VP_VS_RATIO_GOCAD_TOP    = 2.0d0
    VP_VS_RATIO_GOCAD_BOTTOM = 1.732d0
    ANISOTROPY               = .false.
    USE_REGULAR_MESH         = .false.

  else if(MODEL == 'Harvard_LA') then

    BASEMENT_MAP             = .true.
    MOHO_MAP_LUPEI           = .true.
    HAUKSSON_REGIONAL_MODEL  = .true.
    HARVARD_3D_GOCAD_MODEL   = .true.
    THICKNESS_TAPER_BLOCK_HR = 3000.d0
    THICKNESS_TAPER_BLOCK_MR = 15000.d0
    IMPOSE_MINIMUM_VP_GOCAD  = .true.
    VP_MIN_GOCAD             = 750.d0
    VP_VS_RATIO_GOCAD_TOP    = 2.0d0
    VP_VS_RATIO_GOCAD_BOTTOM = 1.732d0
    ANISOTROPY               = .false.
    USE_REGULAR_MESH         = .false.

  else if(MODEL == 'Min_Chen_anisotropy') then

    BASEMENT_MAP             = .false.
    MOHO_MAP_LUPEI           = .false.
    HAUKSSON_REGIONAL_MODEL  = .false.
    HARVARD_3D_GOCAD_MODEL   = .false.
    THICKNESS_TAPER_BLOCK_HR = 1.d0
    THICKNESS_TAPER_BLOCK_MR = 1.d0
    IMPOSE_MINIMUM_VP_GOCAD  = .false.
    VP_MIN_GOCAD             = 1.d0
    VP_VS_RATIO_GOCAD_TOP    = 2.0d0
    VP_VS_RATIO_GOCAD_BOTTOM = 1.732d0
    ANISOTROPY               = .true.
    USE_REGULAR_MESH         = .false.

  else
    stop 'model not implemented, edit read_parameter_file.f90 and recompile'
  endif

! check that Poisson's ratio is positive
  if(VP_VS_RATIO_GOCAD_TOP <= sqrt(2.d0) .or. VP_VS_RATIO_GOCAD_BOTTOM <= sqrt(2.d0)) &
      stop 'wrong value of Poisson''s ratio for Gocad Vs block'

  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) return
  call read_value_double_precision(RECORD_LENGTH_IN_SECONDS, 'solver.RECORD_LENGTH_IN_SECONDS')
  if(err_occurred() /= 0) return

! compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_SECONDS / (100.d0*DT)) + 1)
  if ( NOISE_TOMOGRAPHY /= 0 )   NSTEP = 2*NSTEP-1   ! time steps needs to be doubled, due to +/- branches

! compute the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')
  open(unit=1,file=CMTSOLUTION,iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file'
  icounter = 0
  do while(ios == 0)
    read(1,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(1)
  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) return
  call read_value_logical(CREATE_SHAKEMAP, 'solver.CREATE_SHAKEMAP')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_DISPLACEMENT, 'solver.SAVE_DISPLACEMENT')
  if(err_occurred() /= 0) return
  call read_value_logical(USE_HIGHRES_FOR_MOVIES, 'solver.USE_HIGHRES_FOR_MOVIES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) return
! computes a default hdur_movie that creates nice looking movies.
! Sets HDUR_MOVIE as the minimum period the mesh can resolve for Southern California model
  if(HDUR_MOVIE <=TINYVAL .and. (MODEL == 'Harvard_LA' .or. MODEL == 'SoCal')) &
  HDUR_MOVIE = max(384/NEX_XI*2.4,384/NEX_ETA*2.4)

! compute the minimum value of hdur in CMTSOLUTION file
  open(unit=1,file=CMTSOLUTION,status='old',action='read')
  minval_hdur = HUGEVAL
  do isource = 1,NSOURCES

! skip other information
    do idummy = 1,3
      read(1,"(a)") dummystring
    enddo

! read half duration and compute minimum
    read(1,"(a)") dummystring
    read(dummystring(15:len_trim(dummystring)),*) hdur
    minval_hdur = min(minval_hdur,hdur)

! skip other information
    do idummy = 1,9
      read(1,"(a)") dummystring
    enddo

  enddo
  close(1)

! one cannot use a Heaviside source for the movies
!  if((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. sqrt(minval_hdur**2 + HDUR_MOVIE**2) < TINYVAL) &
!    stop 'hdur too small for movie creation, movies do not make sense for Heaviside source'

  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) return
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) return
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) return

! close parameter file
  close(IIN)

  end subroutine read_parameter_file

