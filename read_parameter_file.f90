!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_parameter_file(LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCKS,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,STACEY_ABS_CONDITIONS, &
        SAVE_AVS_DX_MOVIE,SAVE_DISPLACEMENT,NMOVIE,HDUR_MIN_MOVIES)

  implicit none

  include "constants.h"

  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
            NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP,UTM_PROJECTION_ZONE
  integer NSOURCES,NMOVIE

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,DT
  double precision THICKNESS_TAPER_BLOCKS,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM
  double precision HDUR_MIN_MOVIES

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,STACEY_ABS_CONDITIONS
  logical SAVE_AVS_DX_MOVIE,SAVE_DISPLACEMENT

  character(len=150) LOCAL_PATH

  integer i
  double precision DEPTH_BLOCK_KM

! first 27 characters of each line in the file are a comment
  character(len=27) junk

  open(unit=IIN,file='DATA/Par_file',status='old')

! ignore header
  do i=1,11
    read(IIN,*)
  enddo

  read(IIN,1) junk,NEX_XI
  read(IIN,1) junk,NEX_ETA
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NPROC_XI
  read(IIN,1) junk,NPROC_ETA
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,DT
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSTEP
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,LAT_MIN
  read(IIN,2) junk,LAT_MAX
  read(IIN,2) junk,LONG_MIN
  read(IIN,2) junk,LONG_MAX
  read(IIN,2) junk,DEPTH_BLOCK_KM
  read(IIN,1) junk,UTM_PROJECTION_ZONE

! convert basin size to UTM coordinates and depth of mesh to meters
  call utm_geo(LONG_MIN,LAT_MIN,UTM_X_MIN,UTM_Y_MIN,UTM_PROJECTION_ZONE,ILONGLAT2UTM)
  call utm_geo(LONG_MAX,LAT_MAX,UTM_X_MAX,UTM_Y_MAX,UTM_PROJECTION_ZONE,ILONGLAT2UTM)

  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

  if(dabs(DEPTH_BLOCK_KM) <= DEPTH_MOHO_SOCAL) &
    stop 'bottom of mesh must be deeper than deepest regional layer'

! check that parameters computed are consistent
  if(UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if(UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,TOPOGRAPHY
  read(IIN,3) junk,BASEMENT_MAP
  read(IIN,3) junk,MOHO_MAP_LUPEI
  read(IIN,3) junk,HAUKSSON_REGIONAL_MODEL
  read(IIN,3) junk,HARVARD_3D_GOCAD_MODEL
  read(IIN,2) junk,THICKNESS_TAPER_BLOCKS
  read(IIN,3) junk,IMPOSE_MINIMUM_VP_GOCAD
  read(IIN,2) junk,VP_MIN_GOCAD
  read(IIN,2) junk,VP_VS_RATIO_GOCAD_TOP
  read(IIN,2) junk,VP_VS_RATIO_GOCAD_BOTTOM
  read(IIN,3) junk,ATTENUATION
  read(IIN,3) junk,OCEANS
  read(IIN,3) junk,STACEY_ABS_CONDITIONS

! check that Poisson's ratio is positive
  if(VP_VS_RATIO_GOCAD_TOP <= sqrt(2.) .or. VP_VS_RATIO_GOCAD_BOTTOM <= sqrt(2.)) &
      stop 'wrong value of Poisson''s ratio for Gocad Vs block'

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NER_SEDIM
  read(IIN,1) junk,NER_BASEMENT_SEDIM
  read(IIN,1) junk,NER_16_BASEMENT
  read(IIN,1) junk,NER_MOHO_16
  read(IIN,1) junk,NER_BOTTOM_MOHO
  read(IIN,*)
  read(IIN,*)
  read(IIN,4) junk,LOCAL_PATH

! multiply parameters read by mesh scaling factor
  NER_BASEMENT_SEDIM = NER_BASEMENT_SEDIM * 2
  NER_16_BASEMENT = NER_16_BASEMENT * 2
  NER_MOHO_16 = NER_MOHO_16 * 2
  NER_BOTTOM_MOHO = NER_BOTTOM_MOHO * 4

! ignore name of machine file (used by scripts but not by mesher nor solver)
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSEIS

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSOURCES

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,SAVE_AVS_DX_MOVIE
  read(IIN,3) junk,SAVE_DISPLACEMENT
  read(IIN,1) junk,NMOVIE
  read(IIN,2) junk,HDUR_MIN_MOVIES

! close parameter file
  close(IIN)

! formats

 1 format(a,i8)
 2 format(a,f12.5)
 3 format(a,l8)
 4 format(a,a)

  end subroutine read_parameter_file

