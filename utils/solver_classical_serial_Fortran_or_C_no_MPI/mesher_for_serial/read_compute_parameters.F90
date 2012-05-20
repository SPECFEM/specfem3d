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
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
         NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
         NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
         NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
         NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
         NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
         NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
         CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
         RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
         R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
         MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
         TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
         ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
         MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
         PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
         ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC, &
         NSPEC2D_XI, &
         NSPEC2D_ETA, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
         NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
         NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
         DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
         WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY)


  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
          NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP_KM,MOVIE_TOP,MOVIE_BOTTOM_KM,MOVIE_BOTTOM, &
          MOVIE_EAST_DEG,MOVIE_EAST,MOVIE_WEST_DEG,MOVIE_WEST,MOVIE_NORTH_DEG,MOVIE_NORTH,MOVIE_SOUTH_DEG,MOVIE_SOUTH

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! local variables
  integer NEX_MAX

  double precision RECORD_LENGTH_IN_MINUTES,ELEMENT_WIDTH

  integer, external :: err_occurred

! parameters to be computed based upon parameters above read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
      NSPEC2D_XI, &
      NSPEC2D_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
      NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
      NGLOB

  integer nblocks_xi,nblocks_eta

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                          DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval

! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical :: HONOR_1D_SPHERICAL_MOHO,CASE_3D

  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, padding, tmp_sum, tmp_sum_xi, tmp_sum_eta
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
              nb_lay_sb, nspec_sb, nglob_vol, nglob_surf, nglob_edge

  integer :: multiplication_factor

! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer :: lastdoubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
              normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  integer :: tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,nglob_edge_v,to_remove

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  call open_parameter_file

  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_integer(NCHUNKS, 'mesher.NCHUNKS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'NCHUNKS must be either 1, 2, 3 or 6'

  call read_value_double_precision(ANGULAR_WIDTH_XI_IN_DEGREES, 'mesher.ANGULAR_WIDTH_XI_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(ANGULAR_WIDTH_ETA_IN_DEGREES, 'mesher.ANGULAR_WIDTH_ETA_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(CENTER_LATITUDE_IN_DEGREES, 'mesher.CENTER_LATITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(CENTER_LONGITUDE_IN_DEGREES, 'mesher.CENTER_LONGITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(GAMMA_ROTATION_AZIMUTH, 'mesher.GAMMA_ROTATION_AZIMUTH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

!! DK DK this version of the mesher for the GPU + MPI solver is limited to one chunk for now
  if(NCHUNKS > 1) stop 'this version of the mesher for the GPU + MPI solver is limited to one chunk for now'

! this MUST be 90 degrees for two chunks or more to match geometrically
  if(NCHUNKS > 1 .and. abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_XI_IN_DEGREES must be 90 for more than one chunk'

! this can be any value in the case of two chunks
  if(NCHUNKS > 2 .and. abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_ETA_IN_DEGREES must be 90 for more than two chunks'

! include central cube or not
! use regular cubed sphere instead of cube for large distances
  if(NCHUNKS == 6) then
    INCLUDE_CENTRAL_CUBE = .true.
    INFLATE_CENTRAL_CUBE = .false.
  else
    INCLUDE_CENTRAL_CUBE = .false.
    INFLATE_CENTRAL_CUBE = .true.
  endif

! number of elements at the surface along the two sides of the first chunk
  call read_value_integer(NEX_XI_read, 'mesher.NEX_XI')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NEX_ETA_read, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NPROC_XI_read, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NPROC_ETA_read, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  if(.not. EMULATE_ONLY) then
    NEX_XI = NEX_XI_read
    NEX_ETA = NEX_ETA_read
    NPROC_XI = NPROC_XI_read
    NPROC_ETA = NPROC_ETA_read
  else
! this is used in UTILS/estimate_best_values_runs.f90 only, to estimate memory use
    NEX_ETA = NEX_XI
    NPROC_ETA = NPROC_XI
  endif

! define the velocity model
  call read_value_string(MODEL, 'model.name')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! use PREM as the 1D reference model by default
  REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM

! HONOR_1D_SPHERICAL_MOHO: honor PREM Moho or not: doing so drastically reduces
! the stability condition and therefore the time step, resulting in expensive
! calculations. If not, honor a fictitious Moho at the depth of 40 km
! in order to have even radial sampling from the d220 to the Earth surface.

! ONE_CRUST: in order to increase stability and therefore to allow cheaper
! simulations (larger time step), 1D models can be run with just one average crustal
! layer instead of two.

! CASE_3D : this flag allows the stretching of the elements in the crustal
! layers in the case of 3D models. The purpose of this stretching is to squeeze more
! GLL points per km in the upper part of the crust than in the lower part.
  HONOR_1D_SPHERICAL_MOHO = .false.
  ONE_CRUST = .false.
  CASE_3D = .false.

! default is no 3D model
  THREE_D_MODEL = 0

  if(MODEL == '1D_isotropic_prem') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_transversely_isotropic_prem') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_iasp91' .or. MODEL == '1D_1066a' .or. &
          MODEL == '1D_ak135' .or. MODEL == '1D_jp3d' .or. &
          MODEL == '1D_sea99') then
    if(MODEL == '1D_iasp91') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if(MODEL == '1D_1066a') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if(MODEL == '1D_ak135') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135
   else if(MODEL == '1D_jp3d') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
   else if(MODEL == '1D_sea99') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    else
      stop 'reference 1D Earth model unknown'
    endif
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_ref') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF

  else if(MODEL == '1D_ref_iso') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF

  else if(MODEL == '1D_isotropic_prem_onecrust') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == '1D_transversely_isotropic_prem_onecrust') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == '1D_iasp91_onecrust' .or. MODEL == '1D_1066a_onecrust' .or. MODEL == '1D_ak135_onecrust') then
    if(MODEL == '1D_iasp91_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if(MODEL == '1D_1066a_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if(MODEL == '1D_ak135_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135
    else
      stop 'reference 1D Earth model unknown'
    endif
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == 'transversely_isotropic_prem_plus_3D_crust_2.0') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else if(MODEL == 's20rts') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_S20RTS

  else if(MODEL == 'sea99_jp3d1994') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99_JP3D

  else if(MODEL == 'sea99') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99


  else if(MODEL == 'jp3d1994') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
    THREE_D_MODEL = THREE_D_MODEL_JP3D

  else if(MODEL == 's362ani') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  else if(MODEL == 's362iso') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  else if(MODEL == 's362wmani') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S362WMANI

  else if(MODEL == 's362ani_prem') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_S362ANI_PREM

  else if(MODEL == 's29ea') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S29EA

  else if(MODEL == '3D_attenuation') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .true.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else if(MODEL == '3D_anisotropic') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else
    stop 'model not implemented, edit read_compute_parameters.f90 and recompile'
  endif

! set time step, radial distribution of elements, and attenuation period range
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)

!----
!----  case prem_onecrust by default
!----
      if (SUPPRESS_CRUSTAL_MESH) then
        multiplication_factor=2
      else
        multiplication_factor=1
      endif

     ! element width =   0.5625000      degrees =    62.54715      km
      if(NEX_MAX*multiplication_factor <= 160) then
        DT                       = 0.252d0

        MIN_ATTENUATION_PERIOD   = 30
        MAX_ATTENUATION_PERIOD   = 1500

        NER_CRUST                = 1
        NER_80_MOHO              = 1
        NER_220_80               = 2
        NER_400_220              = 2
        NER_600_400              = 2
        NER_670_600              = 1
        NER_771_670              = 1
        NER_TOPDDOUBLEPRIME_771  = 15
        NER_CMB_TOPDDOUBLEPRIME  = 1
        NER_OUTER_CORE           = 16
        NER_TOP_CENTRAL_CUBE_ICB = 2
        R_CENTRAL_CUBE = 950000.d0

    ! element width =   0.3515625      degrees =    39.09196      km
      else if(NEX_MAX*multiplication_factor <= 256) then
        DT                       = 0.225d0

        MIN_ATTENUATION_PERIOD   = 20
        MAX_ATTENUATION_PERIOD   = 1000

        NER_CRUST                = 1
        NER_80_MOHO              = 1
        NER_220_80               = 2
        NER_400_220              = 3
        NER_600_400              = 3
        NER_670_600              = 1
        NER_771_670              = 1
        NER_TOPDDOUBLEPRIME_771  = 22
        NER_CMB_TOPDDOUBLEPRIME  = 2
        NER_OUTER_CORE           = 24
        NER_TOP_CENTRAL_CUBE_ICB = 3
        R_CENTRAL_CUBE = 965000.d0

    ! element width =   0.2812500      degrees =    31.27357      km
      else if(NEX_MAX*multiplication_factor <= 320) then
        DT                       = 0.16d0

        MIN_ATTENUATION_PERIOD   = 15
        MAX_ATTENUATION_PERIOD   = 750

        NER_CRUST                = 1
        NER_80_MOHO              = 1
        NER_220_80               = 3
        NER_400_220              = 4
        NER_600_400              = 4
        NER_670_600              = 1
        NER_771_670              = 2
        NER_TOPDDOUBLEPRIME_771  = 29
        NER_CMB_TOPDDOUBLEPRIME  = 2
        NER_OUTER_CORE           = 32
        NER_TOP_CENTRAL_CUBE_ICB = 4
        R_CENTRAL_CUBE = 940000.d0

    ! element width =   0.1875000      degrees =    20.84905      km
      else if(NEX_MAX*multiplication_factor <= 480) then
        DT                       = 0.11d0

        MIN_ATTENUATION_PERIOD   = 10
        MAX_ATTENUATION_PERIOD   = 500

        NER_CRUST                = 1
        NER_80_MOHO              = 2
        NER_220_80               = 4
        NER_400_220              = 5
        NER_600_400              = 6
        NER_670_600              = 2
        NER_771_670              = 2
        NER_TOPDDOUBLEPRIME_771  = 44
        NER_CMB_TOPDDOUBLEPRIME  = 3
        NER_OUTER_CORE           = 48
        NER_TOP_CENTRAL_CUBE_ICB = 5
        R_CENTRAL_CUBE = 988000.d0

    ! element width =   0.1757812      degrees =    19.54598      km
      else if(NEX_MAX*multiplication_factor <= 512) then
        DT                       = 0.1125d0

        MIN_ATTENUATION_PERIOD   = 9
        MAX_ATTENUATION_PERIOD   = 500

        NER_CRUST                = 1
        NER_80_MOHO              = 2
        NER_220_80               = 4
        NER_400_220              = 6
        NER_600_400              = 6
        NER_670_600              = 2
        NER_771_670              = 3
        NER_TOPDDOUBLEPRIME_771  = 47
        NER_CMB_TOPDDOUBLEPRIME  = 3
        NER_OUTER_CORE           = 51
        NER_TOP_CENTRAL_CUBE_ICB = 5
        R_CENTRAL_CUBE = 1010000.d0

    ! element width =   0.1406250      degrees =    15.63679      km
      else if(NEX_MAX*multiplication_factor <= 640) then
        DT                       = 0.09d0

        MIN_ATTENUATION_PERIOD   = 8
        MAX_ATTENUATION_PERIOD   = 400

        NER_CRUST                = 2
        NER_80_MOHO              = 3
        NER_220_80               = 5
        NER_400_220              = 7
        NER_600_400              = 8
        NER_670_600              = 3
        NER_771_670              = 3
        NER_TOPDDOUBLEPRIME_771  = 59
        NER_CMB_TOPDDOUBLEPRIME  = 4
        NER_OUTER_CORE           = 64
        NER_TOP_CENTRAL_CUBE_ICB = 6
        R_CENTRAL_CUBE = 1020000.d0

    ! element width =   0.1041667      degrees =    11.58280      km
      else if(NEX_MAX*multiplication_factor <= 864) then
        DT                       = 0.0667d0

        MIN_ATTENUATION_PERIOD   = 6
        MAX_ATTENUATION_PERIOD   = 300

        NER_CRUST                = 2
        NER_80_MOHO              = 4
        NER_220_80               = 6
        NER_400_220              = 10
        NER_600_400              = 10
        NER_670_600              = 3
        NER_771_670              = 4
        NER_TOPDDOUBLEPRIME_771  = 79
        NER_CMB_TOPDDOUBLEPRIME  = 5
        NER_OUTER_CORE           = 86
        NER_TOP_CENTRAL_CUBE_ICB = 9
        R_CENTRAL_CUBE = 990000.d0

    ! element width =   7.8125000E-02  degrees =    8.687103      km
      else if(NEX_MAX*multiplication_factor <= 1152) then
        DT                       = 0.05d0

        MIN_ATTENUATION_PERIOD   = 4
        MAX_ATTENUATION_PERIOD   = 200

        NER_CRUST                = 3
        NER_80_MOHO              = 6
        NER_220_80               = 8
        NER_400_220              = 13
        NER_600_400              = 13
        NER_670_600              = 4
        NER_771_670              = 6
        NER_TOPDDOUBLEPRIME_771  = 106
        NER_CMB_TOPDDOUBLEPRIME  = 7
        NER_OUTER_CORE           = 116
        NER_TOP_CENTRAL_CUBE_ICB = 12
        R_CENTRAL_CUBE = 985000.d0

    ! element width =   7.2115384E-02  degrees =    8.018865      km
      else if(NEX_MAX*multiplication_factor <= 1248) then
        DT                       = 0.0462d0

        MIN_ATTENUATION_PERIOD   = 4
        MAX_ATTENUATION_PERIOD   = 200

        NER_CRUST                = 3
        NER_80_MOHO              = 6
        NER_220_80               = 9
        NER_400_220              = 14
        NER_600_400              = 14
        NER_670_600              = 5
        NER_771_670              = 6
        NER_TOPDDOUBLEPRIME_771  = 114
        NER_CMB_TOPDDOUBLEPRIME  = 8
        NER_OUTER_CORE           = 124
        NER_TOP_CENTRAL_CUBE_ICB = 13
        R_CENTRAL_CUBE = 985000.d0

      else

! scale with respect to 1248 if above that limit
        DT                       = 0.0462d0 * 1248.d0 / (2.d0*NEX_MAX)

        MIN_ATTENUATION_PERIOD   = 4
        MAX_ATTENUATION_PERIOD   = 200

        NER_CRUST                = nint(3 * 2.d0*NEX_MAX / 1248.d0)
        NER_80_MOHO              = nint(6 * 2.d0*NEX_MAX / 1248.d0)
        NER_220_80               = nint(9 * 2.d0*NEX_MAX / 1248.d0)
        NER_400_220              = nint(14 * 2.d0*NEX_MAX / 1248.d0)
        NER_600_400              = nint(14 * 2.d0*NEX_MAX / 1248.d0)
        NER_670_600              = nint(5 * 2.d0*NEX_MAX / 1248.d0)
        NER_771_670              = nint(6 * 2.d0*NEX_MAX / 1248.d0)
        NER_TOPDDOUBLEPRIME_771  = nint(114 * 2.d0*NEX_MAX / 1248.d0)
        NER_CMB_TOPDDOUBLEPRIME  = nint(8 * 2.d0*NEX_MAX / 1248.d0)
        NER_OUTER_CORE           = nint(124 * 2.d0*NEX_MAX / 1248.d0)
        NER_TOP_CENTRAL_CUBE_ICB = nint(13 * 2.d0*NEX_MAX / 1248.d0)
        R_CENTRAL_CUBE = 985000.d0

!! removed this limit           else
!! removed this limit             stop 'problem with this value of NEX_MAX'
      endif

!----
!----  change some values in the case of regular PREM with two crustal layers or of 3D models
!----

! case of regular PREM with two crustal layers: change the time step for small meshes
! because of a different size of elements in the radial direction in the crust
    if (HONOR_1D_SPHERICAL_MOHO) then
      if (.not. ONE_CRUST) then
        ! case 1D + two crustal layers
        if (NER_CRUST<2) NER_CRUST=2
        if(NEX_MAX*multiplication_factor <= 160) then
          DT = 0.20d0
        else if(NEX_MAX*multiplication_factor <= 256) then
          DT = 0.20d0
        endif
      endif
    else
      ! case 3D
      if (NER_CRUST<2) NER_CRUST=2
      if(NEX_MAX*multiplication_factor <= 160) then
          DT = 0.15d0
      else if(NEX_MAX*multiplication_factor <= 256) then
          DT = 0.17d0
      else if(NEX_MAX*multiplication_factor <= 320) then
          DT = 0.155d0
      endif
    endif

    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
      DT = DT*0.20d0
    endif


    if( .not. ATTENUATION_RANGE_PREDEFINED ) then
       call auto_attenuation_periods(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, &
            MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)
    endif

    if(ANGULAR_WIDTH_XI_IN_DEGREES  < 90.0d0 .or. &
       ANGULAR_WIDTH_ETA_IN_DEGREES < 90.0d0 .or. &
       NEX_MAX > 1248) then

     call auto_ner(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, &
          NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
          NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
          R_CENTRAL_CUBE, CASE_3D)

     call auto_attenuation_periods(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, &
          MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)

     call auto_time_stepping(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, DT)

!! DK DK suppressed because this routine should not write anything to the screen
!    write(*,*)'##############################################################'
!    write(*,*)
!    write(*,*)' Auto Radial Meshing Code '
!    write(*,*)' Consult read_compute_parameters.f90 and auto_ner.f90 '
!    write(*,*)' This should only be invoked for chunks less than 90 degrees'
!    write(*,*)' and for chunks greater than 1248 elements wide'
!    write(*,*)
!    write(*,*)'CHUNK WIDTH:              ', ANGULAR_WIDTH_XI_IN_DEGREES
!    write(*,*)'NEX:                      ', NEX_MAX
!    write(*,*)'NER_CRUST:                ', NER_CRUST
!    write(*,*)'NER_80_MOHO:              ', NER_80_MOHO
!    write(*,*)'NER_220_80:               ', NER_220_80
!    write(*,*)'NER_400_220:              ', NER_400_220
!    write(*,*)'NER_600_400:              ', NER_600_400
!    write(*,*)'NER_670_600:              ', NER_670_600
!    write(*,*)'NER_771_670:              ', NER_771_670
!    write(*,*)'NER_TOPDDOUBLEPRIME_771:  ', NER_TOPDDOUBLEPRIME_771
!    write(*,*)'NER_CMB_TOPDDOUBLEPRIME:  ', NER_CMB_TOPDDOUBLEPRIME
!    write(*,*)'NER_OUTER_CORE:           ', NER_OUTER_CORE
!    write(*,*)'NER_TOP_CENTRAL_CUBE_ICB: ', NER_TOP_CENTRAL_CUBE_ICB
!    write(*,*)'R_CENTRAL_CUBE:           ', R_CENTRAL_CUBE
!    write(*,*)'multiplication factor:    ', multiplication_factor
!    write(*,*)
!    write(*,*)'DT:                       ',DT
!    write(*,*)'MIN_ATTENUATION_PERIOD    ',MIN_ATTENUATION_PERIOD
!    write(*,*)'MAX_ATTENUATION_PERIOD    ',MAX_ATTENUATION_PERIOD
!    write(*,*)
!    write(*,*)'##############################################################'

    if (HONOR_1D_SPHERICAL_MOHO) then
      if (.not. ONE_CRUST) then
        ! case 1D + two crustal layers
        if (NER_CRUST<2) NER_CRUST=2
      endif
    else
      ! case 3D
      if (NER_CRUST<2) NER_CRUST=2
    endif
  endif


! take a 5% safety margin on the maximum stable time step
! which was obtained by trial and error
!!!!!!!!  DT = DT * (1.d0 - 0.05d0)

  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ELLIPTICITY, 'model.ELLIPTICITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(GRAVITY, 'model.GRAVITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ROTATION, 'model.ROTATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 6) stop 'cannot have absorbing conditions in the full Earth'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 3) stop 'absorbing conditions not supported for three chunks yet'

  if(ATTENUATION_3D .and. .not. ATTENUATION) stop 'need ATTENUATION to use ATTENUATION_3D'

! radii in PREM or IASP91
! and normalized density at fluid-solid interface on fluid size for coupling
! ROCEAN: radius of the ocean (m)
! RMIDDLE_CRUST: radius of the middle crust (m)
! RMOHO: radius of the Moho (m)
! R80: radius of 80 km discontinuity (m)
! R120: radius of 120 km discontinuity (m) in IASP91
! R220: radius of 220 km discontinuity (m)
! R400: radius of 400 km discontinuity (m)
! R600: radius of 600 km 2nd order discontinuity (m)
! R670: radius of 670 km discontinuity (m)
! R771: radius of 771 km 2nd order discontinuity (m)
! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
! RCMB: radius of CMB (m)
! RICB: radius of ICB (m)

! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
  R120 = -1.d0

! value common to all models
  RHO_OCEANS = 1020.0 / RHOAV

  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then

! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0
    R80  = 6291000.d0
    R120 = 6251000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
! there is no d600 discontinuity in IASP91 therefore this value is useless
! but it needs to be there for compatibility with other subroutines
    R600 = R_EARTH - 600000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then

! our implementation of AK135 has not been checked carefully yet
! therefore let us doublecheck it carefully one day

! values below corrected by Ying Zhou <yingz@gps.caltech.edu>

! AK135 without the 300 meters of mud layer
   ROCEAN = 6368000.d0
   RMIDDLE_CRUST = 6361000.d0
   RMOHO  = 6353000.d0
   R80    = 6291000.d0
   R220   = 6161000.d0
   R400   = 5961000.d0
   R670   = 5711000.d0
   RTOPDDOUBLEPRIME = 3631000.d0
   RCMB   = 3479500.d0
   RICB   = 1217500.d0

! values for AK135 that are not discontinuities
   R600 = 5771000.d0
   R771 = 5611000.d0

   RHO_TOP_OC = 9914.5000 / RHOAV
   RHO_BOTTOM_OC = 12139.1000 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then

! values below corrected by Ying Zhou <yingz@gps.caltech.edu>

! 1066A
   RMOHO = 6360000.d0
   R400 = 5950000.d0
   R600 = 5781000.d0
   R670 = 5700000.d0
   RCMB = 3484300.d0
   RICB = 1229480.d0

! values for 1066A that are not discontinuities
   RTOPDDOUBLEPRIME = 3631000.d0
   R220 = 6161000.d0
   R771 = 5611000.d0
! RMIDDLE_CRUST used only for high resolution FFSW1C model, with 3 elements crust simulations
! mid_crust = 10 km
   RMIDDLE_CRUST = 6361000.d0
   R80 = 6291000.d0

! model 1066A has no oceans, therefore we use the radius of the Earth instead
   ROCEAN = R_EARTH

   RHO_TOP_OC = 9917.4500 / RHOAV
   RHO_BOTTOM_OC = 12160.6500 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then

! REF
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5961000.d0
    R600 = 5771000.d0
    R670 = 5721000.d0
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3479958.d0
    RICB = 1221491.d0

    RHO_TOP_OC = 9903.48 / RHOAV
    RHO_BOTTOM_OC = 12166.35 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D) then

! values below corrected by Min Chen <mchen@gps.caltech.edu>

! jp1d
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6359000.d0
    RMOHO = 6345000.d0
    R80 = 6291000.d0
    R220 = 6161000.d0
    R400 = 5949000.d0
    R600 = 5781000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0
    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then

! SEA1D without the 2 km of mud layer or the 3km water layer
   ROCEAN = 6371000.d0
   RMIDDLE_CRUST = 6361000.d0
   RMOHO  = 6346000.d0
   R80    = 6291000.d0
   R220   = 6161000.d0
   R400   = 5961000.d0
   R670   = 5711000.d0
   RTOPDDOUBLEPRIME = 3631000.d0
   RCMB   = 3485700.d0
   RICB   = 1217100.d0

! values for SEA1D that are not discontinuities
   R600 = 5771000.d0
   R771 = 5611000.d0

   RHO_TOP_OC = 9903.4384 / RHOAV
   RHO_BOTTOM_OC = 12166.5885 / RHOAV

  else

! PREM
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5971000.d0
    R600 = 5771000.d0
    R670 = 5701000.d0
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3480000.d0
    RICB = 1221000.d0

    RHO_TOP_OC = 9903.4384 / RHOAV
    RHO_BOTTOM_OC = 12166.5885 / RHOAV

  endif

! honor the PREM Moho or define a fictitious Moho in order to have even radial sampling
! from the d220 to the Earth surface
  if(HONOR_1D_SPHERICAL_MOHO) then
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO
  else
    RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2
  endif

  call read_value_double_precision(RECORD_LENGTH_IN_MINUTES, 'solver.RECORD_LENGTH_IN_MINUTES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(MOVIE_COARSE,'solver.MOVIE_COARSE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! computes a default hdur_movie that creates nice looking movies.
! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if(HDUR_MOVIE <= TINYVAL) &
    HDUR_MOVIE = 1.2d0*max(240.d0/NEX_XI*18.d0*ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA*18.d0*ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)

  call read_value_integer(MOVIE_VOLUME_TYPE, 'solver.MOVIE_VOLUME_TYPE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_TOP_KM, 'solver.MOVIE_TOP_KM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_BOTTOM_KM, 'solver.MOVIE_BOTTOM_KM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_WEST_DEG, 'solver.MOVIE_WEST_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_EAST_DEG, 'solver.MOVIE_EAST_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_NORTH_DEG, 'solver.MOVIE_NORTH_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(MOVIE_SOUTH_DEG, 'solver.MOVIE_SOUTH_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(MOVIE_START, 'solver.MOVIE_START')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(MOVIE_STOP, 'solver.MOVIE_STOP')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  MOVIE_TOP = (R_EARTH_KM-MOVIE_TOP_KM)/R_EARTH_KM
  MOVIE_BOTTOM = (R_EARTH_KM-MOVIE_BOTTOM_KM)/R_EARTH_KM
  MOVIE_EAST = MOVIE_EAST_DEG * DEGREES_TO_RADIANS
  MOVIE_WEST = MOVIE_WEST_DEG * DEGREES_TO_RADIANS
  MOVIE_NORTH = (90.0d0 - MOVIE_NORTH_DEG) * DEGREES_TO_RADIANS ! converting from latitude to colatitude
  MOVIE_SOUTH = (90.0d0 - MOVIE_SOUTH_DEG) * DEGREES_TO_RADIANS

  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NUMBER_OF_RUNS, 'solver.NUMBER_OF_RUNS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NUMBER_OF_THIS_RUN, 'solver.NUMBER_OF_THIS_RUN')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'solver.NTSTEP_BETWEEN_READ_ADJSRC')
  if(err_occurred() /= 0) return

  call read_value_logical(OUTPUT_SEISMOS_ASCII_TEXT, 'solver.OUTPUT_SEISMOS_ASCII_TEXT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(OUTPUT_SEISMOS_SAC_ALPHANUM, 'solver.OUTPUT_SEISMOS_SAC_ALPHANUM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(OUTPUT_SEISMOS_SAC_BINARY, 'solver.OUTPUT_SEISMOS_SAC_BINARY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ROTATE_SEISMOGRAMS_RT, 'solver.ROTATE_SEISMOGRAMS_RT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(WRITE_SEISMOGRAMS_BY_MASTER, 'solver.WRITE_SEISMOGRAMS_BY_MASTER')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(SAVE_ALL_SEISMOS_IN_ONE_FILE, 'solver.SAVE_ALL_SEISMOS_IN_ONE_FILE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(USE_BINARY_FOR_LARGE_FILE, 'solver.USE_BINARY_FOR_LARGE_FILE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_logical(RECEIVERS_CAN_BE_BURIED, 'solver.RECEIVERS_CAN_BE_BURIED')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')

  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! close parameter file
  call close_parameter_file
!--- check that parameters make sense

  if (OUTPUT_SEISMOS_SAC_ALPHANUM .and. (mod(NTSTEP_BETWEEN_OUTPUT_SEISMOS,5)/=0)) &
    stop 'if OUTPUT_SEISMOS_SAC_ALPHANUM = .true. then NTSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of 5, check the Par_file'

! subsets used to save seismograms must not be larger than the whole time series,
! otherwise we waste memory
  if(NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) then
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP
    if (OUTPUT_SEISMOS_SAC_ALPHANUM .and. (mod(NTSTEP_BETWEEN_OUTPUT_SEISMOS,5)/=0)) &
      stop 'if OUTPUT_SEISMOS_SAC_ALPHANUM = .true. then modified NTSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of 5'
  endif

! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) stop 'wrong size of CUSTOM_REAL for reals'

! check that the parameter file is correct
  if(NGNOD /= 27) stop 'number of control nodes must be 27'
  if(NGNOD == 27 .and. NGNOD2D /= 9) stop 'elements with 27 points should have NGNOD2D = 9'

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) stop 'number of SLS must be 3'

! check number of slices in each direction
  if(NCHUNKS < 1) stop 'must have at least one chunk'
#ifdef USE_MPI
  if(NPROC_XI < 2) stop 'NPROC_XI must be at least 2 for the MPI + GPU version'
  if(NPROC_ETA < 2) stop 'NPROC_ETA must be at least 2 for the MPI + GPU version'
#endif

! check number of chunks
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
     stop 'only one, two, three or six chunks can be meshed'

! check that the central cube can be included
  if(INCLUDE_CENTRAL_CUBE .and. NCHUNKS /= 6) stop 'need six chunks to include central cube'

! check that sphere can be cut into slices without getting negative Jacobian
  if(NEX_XI < 48) stop 'NEX_XI must be greater than 48 to cut the sphere into slices with positive Jacobian'
  if(NEX_ETA < 48) stop 'NEX_ETA must be greater than 48 to cut the sphere into slices with positive Jacobian'

! check that mesh can be coarsened in depth three or four times
  CUT_SUPERBRICK_XI=.false.
  CUT_SUPERBRICK_ETA=.false.

  if (SUPPRESS_CRUSTAL_MESH .and. .not. ADD_4TH_DOUBLING) then
    if(mod(NEX_XI,8) /= 0) stop 'NEX_XI must be a multiple of 8'
    if(mod(NEX_ETA,8) /= 0) stop 'NEX_ETA must be a multiple of 8'
    if(mod(NEX_XI/4,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 4*NPROC_XI'
    if(mod(NEX_ETA/4,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 4*NPROC_ETA'
    if(mod(NEX_XI/8,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/8,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  elseif (SUPPRESS_CRUSTAL_MESH .or. .not. ADD_4TH_DOUBLING) then
    if(mod(NEX_XI,16) /= 0) stop 'NEX_XI must be a multiple of 16'
    if(mod(NEX_ETA,16) /= 0) stop 'NEX_ETA must be a multiple of 16'
    if(mod(NEX_XI/8,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 8*NPROC_XI'
    if(mod(NEX_ETA/8,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 8*NPROC_ETA'
    if(mod(NEX_XI/16,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/16,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  else
    if(mod(NEX_XI,32) /= 0) stop 'NEX_XI must be a multiple of 32'
    if(mod(NEX_ETA,32) /= 0) stop 'NEX_ETA must be a multiple of 32'
    if(mod(NEX_XI/16,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 16*NPROC_XI'
    if(mod(NEX_ETA/16,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 16*NPROC_ETA'
    if(mod(NEX_XI/32,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/32,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  endif

! check that topology is correct if more than two chunks
  if(NCHUNKS > 2 .and. NEX_XI /= NEX_ETA) stop 'must have NEX_XI = NEX_ETA for more than two chunks'
  if(NCHUNKS > 2 .and. NPROC_XI /= NPROC_ETA) stop 'must have NPROC_XI = NPROC_ETA for more than two chunks'

! check that IASP91, AK135, 1066A, JP1D or SEA1D is isotropic
  if((REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91 .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135 .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) .and. TRANSVERSE_ISOTROPY) &
        stop 'models IASP91, AK135, 1066A, JP1D and SEA1D are currently isotropic'

  ELEMENT_WIDTH = ANGULAR_WIDTH_XI_IN_DEGREES/dble(NEX_MAX) * DEGREES_TO_RADIANS

!
!--- compute additional parameters
!

! number of elements horizontally in each slice (i.e. per processor)
! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

! total number of processors in the full Earth composed of the six chunks
  NPROCTOT = NCHUNKS * NPROC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  definition of general mesh parameters below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! find element below top of which we should implement the second doubling in the mantle
! locate element closest to optimal value
  distance_min = HUGEVAL
  do ielem = 2,NER_TOPDDOUBLEPRIME_771
    zval = RTOPDDOUBLEPRIME + ielem * (R771 - RTOPDDOUBLEPRIME) / dble(NER_TOPDDOUBLEPRIME_771)
    distance = abs(zval - (R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL))
    if(distance < distance_min) then
      elem_doubling_mantle = ielem
      distance_min = distance
      DEPTH_SECOND_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo

! find element below top of which we should implement the third doubling in the middle of the outer core
! locate element closest to optimal value
  distance_min = HUGEVAL
! start at element number 4 because we need at least two elements below for the fourth doubling
! implemented at the bottom of the outer core
  do ielem = 4,NER_OUTER_CORE
    zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
    distance = abs(zval - (R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL))
    if(distance < distance_min) then
      elem_doubling_middle_outer_core = ielem
      distance_min = distance
      DEPTH_THIRD_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo

  if (ADD_4TH_DOUBLING) then
! find element below top of which we should implement the fourth doubling in the middle of the outer core
! locate element closest to optimal value
    distance_min = HUGEVAL
! end two elements before the top because we need at least two elements above for the third doubling
! implemented in the middle of the outer core
    do ielem = 2,NER_OUTER_CORE-2
      zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
      distance = abs(zval - (R_EARTH - DEPTH_FOURTH_DOUBLING_OPTIMAL))
      if(distance < distance_min) then
        elem_doubling_bottom_outer_core = ielem
        distance_min = distance
        DEPTH_FOURTH_DOUBLING_REAL = R_EARTH - zval
      endif
    enddo
! make sure that the two doublings in the outer core are found in the right order
    if(elem_doubling_bottom_outer_core >= elem_doubling_middle_outer_core) &
                    stop 'error in location of the two doublings in the outer core'
  endif

  ratio_sampling_array(15) = 0

! define all the layers of the mesh
  if (.not. ADD_4TH_DOUBLING) then

    if (SUPPRESS_CRUSTAL_MESH) then

      ONE_CRUST = .false.
      OCEANS= .false.
      TOPOGRAPHY = .false.
      CRUSTAL = .false.

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 1

  ! now only one region
      ner( 1) = NER_CRUST + NER_80_MOHO
      ner( 2) = 0
      ner( 3) = 0

      ner( 4) = NER_220_80
      ner( 5) = NER_400_220
      ner( 6) = NER_600_400
      ner( 7) = NER_670_600
      ner( 8) = NER_771_670
      ner( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner(10) = elem_doubling_mantle
      ner(11) = NER_CMB_TOPDDOUBLEPRIME
      ner(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(13) = elem_doubling_middle_outer_core
      ner(14) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:9) = 1
      ratio_sampling_array(10:12) = 2
      ratio_sampling_array(13:14) = 4

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:3) = IFLAG_CRUST !!!!! IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(10) = .true.
      this_region_has_a_doubling(13) = .true.
      lastdoubling_layer = 13

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = R80

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80    !!!! now fictitious

      r_top(4) = R80
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = R80 / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious
      rmins(3) = R80 / R_EARTH    !!!! now fictitious

      rmaxs(4) = R80 / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:13) = RCMB / R_EARTH
      rmins(12:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    elseif (ONE_CRUST) then

      NUMBER_OF_MESH_LAYERS = 13
      layer_offset = 0

      ner( 1) = NER_CRUST
      ner( 2) = NER_80_MOHO
      ner( 3) = NER_220_80
      ner( 4) = NER_400_220
      ner( 5) = NER_600_400
      ner( 6) = NER_670_600
      ner( 7) = NER_771_670
      ner( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner( 9) = elem_doubling_mantle
      ner(10) = NER_CMB_TOPDDOUBLEPRIME
      ner(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(12) = elem_doubling_middle_outer_core
      ner(13) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      ratio_sampling_array(2:8) = 2
      ratio_sampling_array(9:11) = 4
      ratio_sampling_array(12:13) = 8

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1) = IFLAG_CRUST
      doubling_index(2) = IFLAG_80_MOHO
      doubling_index(3) = IFLAG_220_80
      doubling_index(4:6) = IFLAG_670_220
      doubling_index(7:10) = IFLAG_MANTLE_NORMAL
      doubling_index(11:12) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(13) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(2)  = .true.
      this_region_has_a_doubling(9)  = .true.
      this_region_has_a_doubling(12) = .true.
      lastdoubling_layer = 12

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

  !!!!!!!!!!! DK DK UGLY: beware, is there a bug when 3D crust crosses anisotropy in the mantle?
  !!!!!!!!!!! DK DK UGLY: i.e. if there is no thick crust there, some elements above the Moho
  !!!!!!!!!!! DK DK UGLY: should be anisotropic but anisotropy is currently only
  !!!!!!!!!!! DK DK UGLY: stored between d220 and MOHO to save memory? Clarify this one day.
  !!!!!!!!!!! DK DK UGLY: The Moho stretching and squishing that Jeroen added to V4.0
  !!!!!!!!!!! DK DK UGLY: should partly deal with this problem.

      r_top(1) = R_EARTH
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80

      r_top(3) = R80
      r_bottom(3) = R220

      r_top(4) = R220
      r_bottom(4) = R400

      r_top(5) = R400
      r_bottom(5) = R600

      r_top(6) = R600
      r_bottom(6) = R670

      r_top(7) = R670
      r_bottom(7) = R771

      r_top(8) = R771
      r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(9) = RTOPDDOUBLEPRIME

      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = RICB

      r_top(13) = RICB
      r_bottom(13) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(2) = R80 / R_EARTH

      rmaxs(3) = R80 / R_EARTH
      rmins(3) = R220 / R_EARTH

      rmaxs(4) = R220 / R_EARTH
      rmins(4) = R400 / R_EARTH

      rmaxs(5) = R400 / R_EARTH
      rmins(5) = R600 / R_EARTH

      rmaxs(6) = R600 / R_EARTH
      rmins(6) = R670 / R_EARTH

      rmaxs(7) = R670 / R_EARTH
      rmins(7) = R771 / R_EARTH

      rmaxs(8:9) = R771 / R_EARTH
      rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(10) = RCMB / R_EARTH

      rmaxs(11:12) = RCMB / R_EARTH
      rmins(11:12) = RICB / R_EARTH

      rmaxs(13) = RICB / R_EARTH
      rmins(13) = R_CENTRAL_CUBE / R_EARTH

    else

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 1
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER)<(R_EARTH-RMIDDLE_CRUST)) then
        ner( 1) = ceiling (NER_CRUST / 2.d0)
        ner( 2) = floor (NER_CRUST / 2.d0)
      else
        ner( 1) = floor (NER_CRUST / 2.d0)
        ner( 2) = ceiling (NER_CRUST / 2.d0)
      endif
      ner( 3) = NER_80_MOHO
      ner( 4) = NER_220_80
      ner( 5) = NER_400_220
      ner( 6) = NER_600_400
      ner( 7) = NER_670_600
      ner( 8) = NER_771_670
      ner( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner(10) = elem_doubling_mantle
      ner(11) = NER_CMB_TOPDDOUBLEPRIME
      ner(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(13) = elem_doubling_middle_outer_core
      ner(14) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      ratio_sampling_array(3:9) = 2
      ratio_sampling_array(10:12) = 4
      ratio_sampling_array(13:14) = 8

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:2) = IFLAG_CRUST
      doubling_index(3) = IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(3)  = .true.
      this_region_has_a_doubling(10) = .true.
      this_region_has_a_doubling(13) = .true.
      this_region_has_a_doubling(14) = .false.
      lastdoubling_layer = 13

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80

      r_top(4) = R80
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMIDDLE_CRUST / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R80 / R_EARTH

      rmaxs(4) = R80 / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:13) = RCMB / R_EARTH
      rmins(12:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    endif
  else
    if (SUPPRESS_CRUSTAL_MESH) then

      ONE_CRUST = .false.
      OCEANS= .false.
      TOPOGRAPHY = .false.
      CRUSTAL = .false.

      NUMBER_OF_MESH_LAYERS = 15
      layer_offset = 1

  ! now only one region
      ner( 1) = NER_CRUST + NER_80_MOHO
      ner( 2) = 0
      ner( 3) = 0

      ner( 4) = NER_220_80
      ner( 5) = NER_400_220
      ner( 6) = NER_600_400
      ner( 7) = NER_670_600
      ner( 8) = NER_771_670
      ner( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner(10) = elem_doubling_mantle
      ner(11) = NER_CMB_TOPDDOUBLEPRIME
      ner(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner(14) = elem_doubling_bottom_outer_core
      ner(15) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:9) = 1
      ratio_sampling_array(10:12) = 2
      ratio_sampling_array(13) = 4
      ratio_sampling_array(14:15) = 8

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:3) = IFLAG_CRUST !!!!! IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:14) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(15) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(10) = .true.
      this_region_has_a_doubling(13) = .true.
      this_region_has_a_doubling(14) = .true.
      lastdoubling_layer = 14

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = R80

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80    !!!! now fictitious

      r_top(4) = R80
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = R80 / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious
      rmins(3) = R80 / R_EARTH    !!!! now fictitious

      rmaxs(4) = R80 / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:14) = RCMB / R_EARTH
      rmins(12:14) = RICB / R_EARTH

      rmaxs(15) = RICB / R_EARTH
      rmins(15) = R_CENTRAL_CUBE / R_EARTH

    elseif (ONE_CRUST) then

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 0

      ner( 1) = NER_CRUST
      ner( 2) = NER_80_MOHO
      ner( 3) = NER_220_80
      ner( 4) = NER_400_220
      ner( 5) = NER_600_400
      ner( 6) = NER_670_600
      ner( 7) = NER_771_670
      ner( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner( 9) = elem_doubling_mantle
      ner(10) = NER_CMB_TOPDDOUBLEPRIME
      ner(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(12) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner(13) = elem_doubling_bottom_outer_core
      ner(14) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      ratio_sampling_array(2:8) = 2
      ratio_sampling_array(9:11) = 4
      ratio_sampling_array(12) = 8
      ratio_sampling_array(13:14) = 16

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1) = IFLAG_CRUST
      doubling_index(2) = IFLAG_80_MOHO
      doubling_index(3) = IFLAG_220_80
      doubling_index(4:6) = IFLAG_670_220
      doubling_index(7:10) = IFLAG_MANTLE_NORMAL
      doubling_index(11:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(2)  = .true.
      this_region_has_a_doubling(9)  = .true.
      this_region_has_a_doubling(12) = .true.
      this_region_has_a_doubling(13) = .true.
      lastdoubling_layer = 13

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

  !!!!!!!!!!! DK DK UGLY: beware, is there a bug when 3D crust crosses anisotropy in the mantle?
  !!!!!!!!!!! DK DK UGLY: i.e. if there is no thick crust there, some elements above the Moho
  !!!!!!!!!!! DK DK UGLY: should be anisotropic but anisotropy is currently only
  !!!!!!!!!!! DK DK UGLY: stored between d220 and MOHO to save memory? Clarify this one day.
  !!!!!!!!!!! DK DK UGLY: The Moho stretching and squishing that Jeroen added to V4.0
  !!!!!!!!!!! DK DK UGLY: should partly deal with this problem.

      r_top(1) = R_EARTH
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80

      r_top(3) = R80
      r_bottom(3) = R220

      r_top(4) = R220
      r_bottom(4) = R400

      r_top(5) = R400
      r_bottom(5) = R600

      r_top(6) = R600
      r_bottom(6) = R670

      r_top(7) = R670
      r_bottom(7) = R771

      r_top(8) = R771
      r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(9) = RTOPDDOUBLEPRIME

      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(2) = R80 / R_EARTH

      rmaxs(3) = R80 / R_EARTH
      rmins(3) = R220 / R_EARTH

      rmaxs(4) = R220 / R_EARTH
      rmins(4) = R400 / R_EARTH

      rmaxs(5) = R400 / R_EARTH
      rmins(5) = R600 / R_EARTH

      rmaxs(6) = R600 / R_EARTH
      rmins(6) = R670 / R_EARTH

      rmaxs(7) = R670 / R_EARTH
      rmins(7) = R771 / R_EARTH

      rmaxs(8:9) = R771 / R_EARTH
      rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(10) = RCMB / R_EARTH

      rmaxs(11:13) = RCMB / R_EARTH
      rmins(11:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    else

      NUMBER_OF_MESH_LAYERS = 15
      layer_offset = 1
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER)<(R_EARTH-RMIDDLE_CRUST)) then
        ner( 1) = ceiling (NER_CRUST / 2.d0)
        ner( 2) = floor (NER_CRUST / 2.d0)
      else
        ner( 1) = floor (NER_CRUST / 2.d0)
        ner( 2) = ceiling (NER_CRUST / 2.d0)
      endif
      ner( 3) = NER_80_MOHO
      ner( 4) = NER_220_80
      ner( 5) = NER_400_220
      ner( 6) = NER_600_400
      ner( 7) = NER_670_600
      ner( 8) = NER_771_670
      ner( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner(10) = elem_doubling_mantle
      ner(11) = NER_CMB_TOPDDOUBLEPRIME
      ner(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner(14) = elem_doubling_bottom_outer_core
      ner(15) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      ratio_sampling_array(3:9) = 2
      ratio_sampling_array(10:12) = 4
      ratio_sampling_array(13) = 8
      ratio_sampling_array(14:15) = 16

  ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:2) = IFLAG_CRUST
      doubling_index(3) = IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:14) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(15) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      this_region_has_a_doubling(3)  = .true.
      this_region_has_a_doubling(10) = .true.
      this_region_has_a_doubling(13) = .true.
      this_region_has_a_doubling(14) = .true.
      lastdoubling_layer = 14

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80

      r_top(4) = R80
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMIDDLE_CRUST / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R80 / R_EARTH

      rmaxs(4) = R80 / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:14) = RCMB / R_EARTH
      rmins(12:14) = RICB / R_EARTH

      rmaxs(15) = RICB / R_EARTH
      rmins(15) = R_CENTRAL_CUBE / R_EARTH
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of elements (NSPEC) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ratio_divide_central_cube = maxval(ratio_sampling_array)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! theoretical number of spectral elements in radial direction
do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
        if(iter_region == IREGION_CRUST_MANTLE) then
                ifirst_region = 1
                ilast_region = 10 + layer_offset
        else if(iter_region == IREGION_OUTER_CORE) then
                ifirst_region = 11 + layer_offset
                ilast_region = NUMBER_OF_MESH_LAYERS - 1
        else if(iter_region == IREGION_INNER_CORE) then
                ifirst_region = NUMBER_OF_MESH_LAYERS
                ilast_region = NUMBER_OF_MESH_LAYERS
        else
                stop 'incorrect region code detected'
        endif
        NSPEC1D_RADIAL(iter_region) = sum(ner(ifirst_region:ilast_region))
enddo

! difference of radial number of element for outer core if the superbrick is cut
  DIFF_NSPEC1D_RADIAL(:,:) = 0
  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC1D_RADIAL(2,1) = 1
      DIFF_NSPEC1D_RADIAL(3,1) = 2
      DIFF_NSPEC1D_RADIAL(4,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(2,2) = 2
      DIFF_NSPEC1D_RADIAL(3,2) = 1

      DIFF_NSPEC1D_RADIAL(1,3) = 1
      DIFF_NSPEC1D_RADIAL(3,3) = 1
      DIFF_NSPEC1D_RADIAL(4,3) = 2

      DIFF_NSPEC1D_RADIAL(1,4) = 2
      DIFF_NSPEC1D_RADIAL(2,4) = 1
      DIFF_NSPEC1D_RADIAL(4,4) = 1
    else
      DIFF_NSPEC1D_RADIAL(2,1) = 1
      DIFF_NSPEC1D_RADIAL(3,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(4,2) = 1
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC1D_RADIAL(3,1) = 1
      DIFF_NSPEC1D_RADIAL(4,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(2,2) = 1
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of surface elements for faces along XI and ETA

do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if(iter_region == IREGION_CRUST_MANTLE) then
        ifirst_region = 1
        ilast_region = 10 + layer_offset
    else if(iter_region == IREGION_OUTER_CORE) then
        ifirst_region = 11 + layer_offset
        ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if(iter_region == IREGION_INNER_CORE) then
        ifirst_region = NUMBER_OF_MESH_LAYERS
        ilast_region = NUMBER_OF_MESH_LAYERS
    else
        stop 'incorrect region code detected'
    endif
    tmp_sum_xi = 0
    tmp_sum_eta = 0
    tmp_sum_nglob2D_xi = 0
    tmp_sum_nglob2D_eta = 0
    do iter_layer = ifirst_region, ilast_region
        if (this_region_has_a_doubling(iter_layer)) then
            if (iter_region == IREGION_OUTER_CORE .and. iter_layer == lastdoubling_layer) then
              ! simple brick
              divider = 1
              nglob_surf = 6*NGLLX**2 - 7*NGLLX + 2
              nglob_edges_h = 2*(NGLLX-1)+1 + NGLLX
              ! minimum value to be safe
              nglob_edge_v = NGLLX-2
              nb_lay_sb = 2
              nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK
              nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK
            else
              ! double brick
              divider = 2
              if (ner(iter_layer) == 1) then
                nglob_surf = 6*NGLLX**2 - 8*NGLLX + 3
                nglob_edges_h = 4*(NGLLX-1)+1 + 2*(NGLLX-1)+1
                nglob_edge_v = NGLLX-2
                nb_lay_sb = 1
                nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK_1L
                nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK_1L
              else
                nglob_surf = 8*NGLLX**2 - 11*NGLLX + 4
                nglob_edges_h = 4*(NGLLX-1)+1 + 2*(NGLLX-1)+1
                nglob_edge_v = 2*(NGLLX-1)+1 -2
                nb_lay_sb = 2
                nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK
                nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK
                divider = 2
              endif
            endif
            doubling = 1
            to_remove = 1
        else
            if (iter_layer /= ifirst_region) then
              to_remove = 0
            else
              to_remove = 1
            endif
            ! dummy value
            divider = 1
            doubling = 0
            nb_lay_sb = 0
            nspec2D_xi_sb = 0
            nspec2D_eta_sb = 0
            nglob_surf = 0
            nglob_edges_h = 0
            nglob_edge_v = 0
        endif

        tmp_sum_xi = tmp_sum_xi + ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * (nspec2D_xi_sb/2))

        tmp_sum_eta = tmp_sum_eta + ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * (nspec2D_eta_sb/2))

        tmp_sum_nglob2D_xi = tmp_sum_nglob2D_xi + (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb))*NGLLX*NGLLX) - &
                ((((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - doubling*nb_lay_sb)) + &
                ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))*(ner(iter_layer) - to_remove - doubling*nb_lay_sb))*NGLLX) + &
                (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - to_remove - doubling*nb_lay_sb)) + &
                doubling * (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))/divider) * (nglob_surf-nglob_edges_h) - &
                ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))/divider -1) * nglob_edge_v)

        tmp_sum_nglob2D_eta = tmp_sum_nglob2D_eta + (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb))*NGLLX*NGLLX) - &
                ((((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - doubling*nb_lay_sb)) + &
                ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))* &
                   (ner(iter_layer) - to_remove - doubling*nb_lay_sb))*NGLLX) + &
                (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - to_remove - doubling*nb_lay_sb)) + &
                doubling * (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))/divider) * (nglob_surf-nglob_edges_h) - &
                ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))/divider -1) * nglob_edge_v)
    enddo
    NSPEC2D_XI(iter_region) = tmp_sum_xi
    NSPEC2D_ETA(iter_region) = tmp_sum_eta

    NGLOB2DMAX_YMIN_YMAX(iter_region) = tmp_sum_nglob2D_xi
    NGLOB2DMAX_XMIN_XMAX(iter_region) = tmp_sum_nglob2D_eta

    if (iter_region == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
        NSPEC2D_XI(iter_region) = NSPEC2D_XI(iter_region) + &
        ((NEX_PER_PROC_XI / ratio_divide_central_cube)*(NEX_XI / ratio_divide_central_cube))
        NSPEC2D_ETA(iter_region) = NSPEC2D_ETA(iter_region) + &
        ((NEX_PER_PROC_ETA / ratio_divide_central_cube)*(NEX_XI / ratio_divide_central_cube))

        NGLOB2DMAX_YMIN_YMAX(iter_region) = NGLOB2DMAX_YMIN_YMAX(iter_region) + &
        (((NEX_PER_PROC_XI / ratio_divide_central_cube)*(NGLLX-1)+1)*((NEX_XI / ratio_divide_central_cube)*(NGLLX-1)+1))

        NGLOB2DMAX_XMIN_XMAX(iter_region) = NGLOB2DMAX_XMIN_XMAX(iter_region) + &
        (((NEX_PER_PROC_ETA / ratio_divide_central_cube)*(NGLLX-1)+1)*((NEX_XI / ratio_divide_central_cube)*(NGLLX-1)+1))
    endif
enddo

! difference of number of surface elements along xi or eta for outer core if the superbrick is cut
  DIFF_NSPEC2D_XI(:,:) = 0
  DIFF_NSPEC2D_ETA(:,:) = 0
  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC2D_XI(2,1) = 2
      DIFF_NSPEC2D_XI(1,2) = 2
      DIFF_NSPEC2D_XI(2,3) = 2
      DIFF_NSPEC2D_XI(1,4) = 2

      DIFF_NSPEC2D_ETA(2,1) = 1
      DIFF_NSPEC2D_ETA(2,2) = 1
      DIFF_NSPEC2D_ETA(1,3) = 1
      DIFF_NSPEC2D_ETA(1,4) = 1
    else
      DIFF_NSPEC2D_ETA(2,1) = 1
      DIFF_NSPEC2D_ETA(1,2) = 1
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC2D_XI(2,1) = 2
      DIFF_NSPEC2D_XI(1,2) = 2
    endif
  endif
  DIFF_NSPEC2D_XI(:,:) = DIFF_NSPEC2D_XI(:,:) * (NEX_PER_PROC_XI / ratio_divide_central_cube)
  DIFF_NSPEC2D_ETA(:,:) = DIFF_NSPEC2D_ETA(:,:) * (NEX_PER_PROC_ETA / ratio_divide_central_cube)

! exact number of surface elements on the bottom and top boundaries

! in the crust and mantle
  NSPEC2D_TOP(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(1))*(NEX_ETA/ratio_sampling_array(1))/NPROC
  NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(10+layer_offset))*&
                                         (NEX_ETA/ratio_sampling_array(10+layer_offset))/NPROC

! in the outer core with mesh doubling
  if (ADD_4TH_DOUBLING) then
    NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/(ratio_divide_central_cube/4))*(NEX_ETA/(ratio_divide_central_cube/4))/NPROC
    NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  else
    NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/(ratio_divide_central_cube/2))*(NEX_ETA/(ratio_divide_central_cube/2))/NPROC
    NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  endif

! in the top of the inner core
  NSPEC2D_TOP(IREGION_INNER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  NSPEC2D_BOTTOM(IREGION_INNER_CORE) = NSPEC2D_TOP(IREGION_INNER_CORE)

! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX(:) = NSPEC2D_ETA(:)
  NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE) = NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE) + maxval(DIFF_NSPEC2D_ETA(:,:))
  NSPEC2DMAX_YMIN_YMAX(:) = NSPEC2D_XI(:)
  NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE) = NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE) + maxval(DIFF_NSPEC2D_XI(:,:))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of spectral elements in each region

do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if(iter_region == IREGION_CRUST_MANTLE) then
        ifirst_region = 1
        ilast_region = 10 + layer_offset
    else if(iter_region == IREGION_OUTER_CORE) then
        ifirst_region = 11 + layer_offset
        ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if(iter_region == IREGION_INNER_CORE) then
        ifirst_region = NUMBER_OF_MESH_LAYERS
        ilast_region = NUMBER_OF_MESH_LAYERS
    else
        stop 'incorrect region code detected'
    endif
    tmp_sum = 0;
    do iter_layer = ifirst_region, ilast_region
        if (this_region_has_a_doubling(iter_layer)) then
            if (ner(iter_layer) == 1) then
              nb_lay_sb = 1
              nspec_sb = NSPEC_SUPERBRICK_1L
            else
              nb_lay_sb = 2
              nspec_sb = NSPEC_DOUBLING_SUPERBRICK
            endif
            doubling = 1
        else
            doubling = 0
            nb_lay_sb = 0
            nspec_sb = 0
        endif
        tmp_sum = tmp_sum + (((NEX_XI / ratio_sampling_array(iter_layer)) * (NEX_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_XI / ratio_sampling_array(iter_layer)) * (NEX_ETA / ratio_sampling_array(iter_layer)) * &
                (nspec_sb/4))) / NPROC
    enddo
    NSPEC(iter_region) = tmp_sum
enddo

  if(INCLUDE_CENTRAL_CUBE) NSPEC(IREGION_INNER_CORE) = NSPEC(IREGION_INNER_CORE) + &
         (NEX_PER_PROC_XI / ratio_divide_central_cube) * &
         (NEX_PER_PROC_ETA / ratio_divide_central_cube) * &
         (NEX_XI / ratio_divide_central_cube)

  if(minval(NSPEC) <= 0) stop 'negative NSPEC, there is a problem somewhere'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of points (NGLOB) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! theoretical number of Gauss-Lobatto points in radial direction
  NGLOB1D_RADIAL(:) = NSPEC1D_RADIAL(:)*(NGLLZ-1)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NGLOB2DMAX_XMIN_XMAX(:) = NGLOB2DMAX_XMIN_XMAX(:) + 1
  NGLOB2DMAX_YMIN_YMAX(:) = NGLOB2DMAX_YMIN_YMAX(:) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of global points in each region

! initialize array
  NGLOB(:) = 0

! in the inner core (no doubling region + eventually central cube)
  if(INCLUDE_CENTRAL_CUBE) then
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/ratio_divide_central_cube) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/ratio_divide_central_cube) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB + NEX_XI / ratio_divide_central_cube)*(NGLLZ-1)+1)
  else
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/ratio_divide_central_cube) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/ratio_divide_central_cube) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB)*(NGLLZ-1)+1)
  endif

! in the crust-mantle and outercore
  do iter_region = IREGION_CRUST_MANTLE,IREGION_OUTER_CORE
      if(iter_region == IREGION_CRUST_MANTLE) then
            ifirst_region = 1
            ilast_region = 10 + layer_offset
      else if(iter_region == IREGION_OUTER_CORE) then
            ifirst_region = 11 + layer_offset
            ilast_region = NUMBER_OF_MESH_LAYERS - 1
      else
            stop 'incorrect region code detected'
      endif
      tmp_sum = 0;
      do iter_layer = ifirst_region, ilast_region
        nglob_int_surf_eta=0
        nglob_int_surf_xi=0
        nglob_ext_surf = 0
        nglob_center_edge = 0
        nglob_corner_edge = 0
        nglob_border_edge = 0
        if (this_region_has_a_doubling(iter_layer)) then
            if (iter_region == IREGION_OUTER_CORE .and. iter_layer == lastdoubling_layer .and. &
               (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
              doubling = 1
              normal_doubling = 0
              cut_doubling = 1
              nb_lay_sb = 2
              nglob_edge = 0
              nglob_surf = 0
              nglob_vol = 8*NGLLX**3 - 12*NGLLX**2 + 6*NGLLX - 1
              nglob_int_surf_eta = 6*NGLLX**2 - 7*NGLLX + 2
              nglob_int_surf_xi = 5*NGLLX**2 - 5*NGLLX + 1
              nglob_ext_surf = 4*NGLLX**2-4*NGLLX+1
              nglob_center_edge = 4*(NGLLX-1)+1
              nglob_corner_edge = 2*(NGLLX-1)+1
              nglob_border_edge = 3*(NGLLX-1)+1
            else
              if (ner(iter_layer) == 1) then
                nb_lay_sb = 1
                nglob_vol = 28*NGLLX**3 - 62*NGLLX**2 + 47*NGLLX - 12
                nglob_surf = 6*NGLLX**2-8*NGLLX+3
                nglob_edge = NGLLX
              else
                nb_lay_sb = 2
                nglob_vol = 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13
                nglob_surf = 8*NGLLX**2-11*NGLLX+4
                nglob_edge = 2*NGLLX-1
              endif
              doubling = 1
              normal_doubling = 1
              cut_doubling = 0
            endif
            padding = -1
        else
            doubling = 0
            normal_doubling = 0
            cut_doubling = 0
            padding = 0
            nb_lay_sb = 0
            nglob_vol = 0
            nglob_surf = 0
            nglob_edge = 0
        endif
        if (iter_layer == ilast_region) padding = padding +1
        nblocks_xi = NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)
        nblocks_eta = NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)

        tmp_sum = tmp_sum + &
        ((nblocks_xi)*(NGLLX-1)+1) * ((nblocks_eta)*(NGLLX-1)+1) * ((ner(iter_layer) - doubling*nb_lay_sb)*(NGLLX-1)+padding)+&
        normal_doubling * ((((nblocks_xi*nblocks_eta)/4)*nglob_vol) - &
        (((nblocks_eta/2-1)*nblocks_xi/2+(nblocks_xi/2-1)*nblocks_eta/2)*nglob_surf) + &
        ((nblocks_eta/2-1)*(nblocks_xi/2-1)*nglob_edge)) + &
        cut_doubling*(nglob_vol*(nblocks_xi*nblocks_eta) - &
            ( nblocks_eta*(int(nblocks_xi/2)*nglob_int_surf_xi + int((nblocks_xi-1)/2)*nglob_ext_surf) + &
              nblocks_xi*(int(nblocks_eta/2)*nglob_int_surf_eta + int((nblocks_eta-1)/2)*nglob_ext_surf)&
            ) + &
            ( int(nblocks_xi/2)*int(nblocks_eta/2)*nglob_center_edge + &
              int((nblocks_xi-1)/2)*int((nblocks_eta-1)/2)*nglob_corner_edge + &
              ((int(nblocks_eta/2)*int((nblocks_xi-1)/2))+(int((nblocks_eta-1)/2)*int(nblocks_xi/2)))*nglob_border_edge&
            ))
      enddo
      NGLOB(iter_region) = tmp_sum
  enddo

!!! example :
!!!                        nblocks_xi/2=5
!!!                  ____________________________________
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!! nblocks_eta/2=3  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______I______I______I______I______I
!!!
!!! NGLOB for this doubling layer = 3*5*Volume - ((3-1)*5+(5-1)*3)*Surface + (3-1)*(5-1)*Edge
!!!
!!! 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13 -> nb GLL points in a superbrick (Volume)
!!! 8*NGLLX**2-11*NGLLX+4 -> nb GLL points on a superbrick side (Surface)
!!! 2*NGLLX-1 -> nb GLL points on a corner edge of a superbrick (Edge)

!!! for the one layer superbrick :
!!! NGLOB = 28.NGLL^3 - 62.NGLL^2 + 47.NGLL - 12 (Volume)
!!! NGLOB = 6.NGLL^2 - 8.NGLL + 3 (Surface)
!!! NGLOB = NGLL (Edge)
!!!
!!! those results were obtained by using the script UTILS/doubling_brick/count_nglob_analytical.pl
!!! with an opendx file of the superbrick's geometry

!!! for the basic doubling bricks (two layers)
!!! NGLOB = 8.NGLL^3 - 12.NGLL^2 + 6.NGLL - 1 (VOLUME)
!!! NGLOB = 5.NGLL^2 - 5.NGLL + 1 (SURFACE 1)
!!! NGLOB = 6.NGLL^2 - 7.NGLL + 2 (SURFACE 2)

  end subroutine read_compute_parameters
