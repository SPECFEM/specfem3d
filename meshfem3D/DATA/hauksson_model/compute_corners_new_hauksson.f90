
  program jfdkfd

!! DK DK compute corners of new interpolated Hauksson model

  implicit none

  include "../../constants.h"

  ! UTM projection zone
  integer, parameter :: UTM_PROJECTION_ZONE = 11
  logical, parameter :: SUPPRESS_UTM_PROJECTION = .false.

  double precision, parameter :: ORIG_LONG_HAUKSSON = -121.d0,END_LONG_HAUKSSON = -114.d0
  double precision, parameter :: ORIG_LAT_HAUKSSON = 32.d0,END_LAT_HAUKSSON = 37.d0

  double precision utm_x,utm_y

  call utm_geo(ORIG_LONG_HAUKSSON,ORIG_LAT_HAUKSSON,utm_x,utm_y,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
  print *,'  double precision, parameter :: UTM_X_ORIG_HAUKSSON = ',utm_x
  print *,'  double precision, parameter :: UTM_Y_ORIG_HAUKSSON = ',utm_y

  call utm_geo(END_LONG_HAUKSSON,END_LAT_HAUKSSON,utm_x,utm_y,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
  print *,'  double precision, parameter :: UTM_X_END_HAUKSSON = ',utm_x
  print *,'  double precision, parameter :: UTM_Y_END_HAUKSSON = ',utm_y

  end

!! DK DK add UTM projection routine
  include "../../utm_geo.f90"

