
  subroutine interpolate_gocad_block_HR(vp_block_gocad_HR,vp_block_gocad_MR, &
      utm_x_eval,utm_y_eval,z_eval,rho_final,vp_final,vs_final,point_is_in_sediments, &
      VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
      IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
      vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL, MOHO_MAP_LUPEI)

  implicit none

  include "constants.h"
  include "constants_gocad.h"

  double precision vp_block_gocad_HR(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  double precision vp_block_gocad_MR(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)

  integer ix,iy,iz

  double precision utm_x_eval,utm_y_eval,z_eval
  double precision spacing_x,spacing_y,spacing_z
  double precision gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision v1,v2,v3,v4,v5,v6,v7,v8
  double precision vp_final,vs_final,rho_final,vp_vs_ratio
  double precision rho_ref_MR,vp_ref_MR,vs_ref_MR
  double precision THICKNESS_TAPER_BLOCK_HR, &
      VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical point_is_in_sediments,dummy_flag,IMPOSE_MINIMUM_VP_GOCAD

! for Hauksson's model
  integer doubling_index
  logical HAUKSSON_REGIONAL_MODEL,MOHO_MAP_LUPEI
  double precision, dimension(NLAYERS_HAUKSSON,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: vp_hauksson,vs_hauksson

! determine spacing and cell for linear interpolation
  spacing_x = (utm_x_eval - ORIG_X_GOCAD_HR) / SPACING_X_GOCAD_HR
  spacing_y = (utm_y_eval - ORIG_Y_GOCAD_HR) / SPACING_Y_GOCAD_HR
  spacing_z = (z_eval - ORIG_Z_GOCAD_HR) / SPACING_Z_GOCAD_HR

  ix = int(spacing_x)
  iy = int(spacing_y)
  iz = int(spacing_z)

  gamma_interp_x = spacing_x - dble(ix)
  gamma_interp_y = spacing_y - dble(iy)
  gamma_interp_z = spacing_z - dble(iz)

! this block is smaller than the grid, therefore just exit
! if the target point is outside of the block
  if(ix < 0 .or. ix > NX_GOCAD_HR-2 .or. iy < 0 .or. iy > NY_GOCAD_HR-2) return

! suppress edge effects in vertical direction
  if(iz < 0) then
    iz = 0
    gamma_interp_z = 0.d0
  endif
  if(iz > NZ_GOCAD_HR-2) then
    iz = NZ_GOCAD_HR-2
    gamma_interp_z = 1.d0
  endif

! define 8 corners of interpolation element
   v1 = vp_block_gocad_HR(ix,iy,iz)
   v2 = vp_block_gocad_HR(ix+1,iy,iz)
   v3 = vp_block_gocad_HR(ix+1,iy+1,iz)
   v4 = vp_block_gocad_HR(ix,iy+1,iz)

   v5 = vp_block_gocad_HR(ix,iy,iz+1)
   v6 = vp_block_gocad_HR(ix+1,iy,iz+1)
   v7 = vp_block_gocad_HR(ix+1,iy+1,iz+1)
   v8 = vp_block_gocad_HR(ix,iy+1,iz+1)

! check if element is defined (i.e. is in the sediments in Voxet)
! do nothing if element is undefined
! a P-velocity of 20 km/s is used to indicate fictitious elements
   if(v1 < 19000. .and. v2 < 19000. .and. &
      v3 < 19000. .and. v4 < 19000. .and. &
      v5 < 19000. .and. v6 < 19000. .and. &
      v7 < 19000. .and. v8 < 19000.) then

! set flag indicating whether point is in the sediments
         point_is_in_sediments = .true.

! use trilinear interpolation in cell to define Vp
         vp_final = &
           v1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z) + &
           v2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z) + &
           v3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z) + &
           v4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z) + &
           v5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z + &
           v6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z + &
           v7*gamma_interp_x*gamma_interp_y*gamma_interp_z + &
           v8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z

! impose minimum velocity if needed
         if(IMPOSE_MINIMUM_VP_GOCAD .and. vp_final < VP_MIN_GOCAD) vp_final = VP_MIN_GOCAD

! taper edges to make smooth transition between MR and HR blocks
! get value from edge of medium-resolution block
! then use linear interpolation from edge of the model
  if(TAPER_GOCAD_TRANSITIONS) then

! x = xmin
  if(utm_x_eval < ORIG_X_GOCAD_HR + THICKNESS_TAPER_BLOCK_HR) then
    gamma_interp_x = (utm_x_eval - ORIG_X_GOCAD_HR) / THICKNESS_TAPER_BLOCK_HR
    call interpolate_gocad_block_MR(vp_block_gocad_MR, &
              ORIG_X_GOCAD_HR,utm_y_eval,z_eval,rho_ref_MR,vp_ref_MR,vs_ref_MR,dummy_flag, &
              VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
              IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
              vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL,MOHO_MAP_LUPEI)
    vp_final = vp_ref_MR * (1. - gamma_interp_x) + vp_final * gamma_interp_x

! x = xmax
  else if(utm_x_eval > END_X_GOCAD_HR - THICKNESS_TAPER_BLOCK_HR) then
    gamma_interp_x = (utm_x_eval - (END_X_GOCAD_HR - THICKNESS_TAPER_BLOCK_HR)) / THICKNESS_TAPER_BLOCK_HR
    call interpolate_gocad_block_MR(vp_block_gocad_MR, &
              END_X_GOCAD_HR,utm_y_eval,z_eval,rho_ref_MR,vp_ref_MR,vs_ref_MR,dummy_flag, &
              VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
              IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
              vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL, MOHO_MAP_LUPEI)
    vp_final = vp_ref_MR * gamma_interp_x + vp_final * (1. - gamma_interp_x)

! y = ymin
  else if(utm_y_eval < ORIG_Y_GOCAD_HR + THICKNESS_TAPER_BLOCK_HR) then
    gamma_interp_y = (utm_y_eval - ORIG_Y_GOCAD_HR) / THICKNESS_TAPER_BLOCK_HR
    call interpolate_gocad_block_MR(vp_block_gocad_MR, &
              utm_x_eval,ORIG_Y_GOCAD_HR,z_eval,rho_ref_MR,vp_ref_MR,vs_ref_MR,dummy_flag, &
              VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
              IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
              vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL, MOHO_MAP_LUPEI)
    vp_final = vp_ref_MR * (1. - gamma_interp_y) + vp_final * gamma_interp_y

! y = ymax
  else if(utm_y_eval > END_Y_GOCAD_HR - THICKNESS_TAPER_BLOCK_HR) then
    gamma_interp_y = (utm_y_eval - (END_Y_GOCAD_HR - THICKNESS_TAPER_BLOCK_HR)) / THICKNESS_TAPER_BLOCK_HR
    call interpolate_gocad_block_MR(vp_block_gocad_MR, &
              utm_x_eval,END_Y_GOCAD_HR,z_eval,rho_ref_MR,vp_ref_MR,vs_ref_MR,dummy_flag, &
              VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
              IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
              vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL, MOHO_MAP_LUPEI)
    vp_final = vp_ref_MR * gamma_interp_y + vp_final * (1. - gamma_interp_y)

  endif

  endif

! use linear variation of vp/vs ratio with depth, between 0. and 8.5 km
         vp_vs_ratio = VP_VS_RATIO_GOCAD_BOTTOM + &
           (VP_VS_RATIO_GOCAD_TOP - VP_VS_RATIO_GOCAD_BOTTOM) * &
           (z_eval - (-8500.d0)) / (0.d0 - (-8500.d0))

! make sure ratio remains in interval
  if(vp_vs_ratio < VP_VS_RATIO_GOCAD_BOTTOM) vp_vs_ratio = VP_VS_RATIO_GOCAD_BOTTOM
  if(vp_vs_ratio > VP_VS_RATIO_GOCAD_TOP) vp_vs_ratio = VP_VS_RATIO_GOCAD_TOP

         vs_final = vp_final / vp_vs_ratio
         call compute_rho_estimate(rho_final,vp_final)

     endif

  end subroutine interpolate_gocad_block_HR

