
  subroutine interpolate_gocad_block_MR(vp_block_gocad_MR, &
      utm_x_eval,utm_y_eval,z_eval,rho_final,vp_final,vs_final,point_is_in_sediments, &
      VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
      IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_MR, &
      vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL,MOHO_MAP_LUPEI)

  implicit none

  include "constants.h"
  include "constants_gocad.h"

  double precision vp_block_gocad_MR(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)
  double precision VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM,THICKNESS_TAPER_BLOCK_MR

  integer ix,iy,iz

  double precision utm_x_eval,utm_y_eval,z_eval
  double precision spacing_x,spacing_y,spacing_z
  double precision gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision v1,v2,v3,v4,v5,v6,v7,v8
  double precision vp_final,vs_final,rho_final,vp_vs_ratio
  double precision xmesh,ymesh,zmesh,vs_dummy,rho_dummy

  logical point_is_in_sediments,IMPOSE_MINIMUM_VP_GOCAD

! for Hauksson's model
  integer doubling_index
  logical HAUKSSON_REGIONAL_MODEL,MOHO_MAP_LUPEI
  double precision vp_ref_hauksson
  double precision, dimension(NLAYERS_HAUKSSON,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: vp_hauksson,vs_hauksson

! determine spacing and cell for linear interpolation
  spacing_x = (utm_x_eval - ORIG_X_GOCAD_MR) / SPACING_X_GOCAD_MR
  spacing_y = (utm_y_eval - ORIG_Y_GOCAD_MR) / SPACING_Y_GOCAD_MR
  spacing_z = (z_eval - ORIG_Z_GOCAD_MR) / SPACING_Z_GOCAD_MR

  ix = int(spacing_x)
  iy = int(spacing_y)
  iz = int(spacing_z)

  gamma_interp_x = spacing_x - dble(ix)
  gamma_interp_y = spacing_y - dble(iy)
  gamma_interp_z = spacing_z - dble(iz)

! suppress edge effects for points outside of Gocad model
  if (ix < 0) then
    ix = 0
    gamma_interp_x = 0.d0
  endif
  if (ix > NX_GOCAD_MR-2) then
    ix = NX_GOCAD_MR-2
    gamma_interp_x = 1.d0
  endif

  if (iy < 0) then
    iy = 0
    gamma_interp_y = 0.d0
  endif
  if (iy > NY_GOCAD_MR-2) then
    iy = NY_GOCAD_MR-2
    gamma_interp_y = 1.d0
  endif

  if (iz < 0) then
    iz = 0
    gamma_interp_z = 0.d0
  endif
  if (iz > NZ_GOCAD_MR-2) then
    iz = NZ_GOCAD_MR-2
    gamma_interp_z = 1.d0
  endif

! define 8 corners of interpolation element
   v1 = vp_block_gocad_MR(ix,iy,iz)
   v2 = vp_block_gocad_MR(ix+1,iy,iz)
   v3 = vp_block_gocad_MR(ix+1,iy+1,iz)
   v4 = vp_block_gocad_MR(ix,iy+1,iz)

   v5 = vp_block_gocad_MR(ix,iy,iz+1)
   v6 = vp_block_gocad_MR(ix+1,iy,iz+1)
   v7 = vp_block_gocad_MR(ix+1,iy+1,iz+1)
   v8 = vp_block_gocad_MR(ix,iy+1,iz+1)

! check if element is defined (i.e. is in the sediments in Voxet)
! do nothing if element is undefined
! a P-velocity of 20 km/s is used to indicate fictitious elements
   if (v1 < 19000. .and. v2 < 19000. .and. &
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
         if (IMPOSE_MINIMUM_VP_GOCAD .and. vp_final < VP_MIN_GOCAD) vp_final = VP_MIN_GOCAD

! taper edges to make smooth transition between Hauksson and MR blocks
! get value from edge of medium-resolution block
! then use linear interpolation from edge of the model
  if (TAPER_GOCAD_TRANSITIONS) then

! x = xmin
  if (utm_x_eval < ORIG_X_GOCAD_MR + THICKNESS_TAPER_BLOCK_MR) then
    xmesh = ORIG_X_GOCAD_MR
    ymesh = utm_y_eval
    zmesh = z_eval
    if (HAUKSSON_REGIONAL_MODEL) then
      call hauksson_model(vp_hauksson,vs_hauksson,xmesh,ymesh,zmesh,vp_ref_hauksson,vs_dummy, MOHO_MAP_LUPEI)
    else
      call socal_model(doubling_index,rho_dummy,vp_ref_hauksson,vs_dummy)
    endif
    gamma_interp_x = (utm_x_eval - ORIG_X_GOCAD_MR) / THICKNESS_TAPER_BLOCK_MR
    vp_final = vp_ref_hauksson * (1. - gamma_interp_x) + vp_final * gamma_interp_x

! x = xmax
  else if (utm_x_eval > END_X_GOCAD_MR - THICKNESS_TAPER_BLOCK_MR) then
    xmesh = END_X_GOCAD_MR
    ymesh = utm_y_eval
    zmesh = z_eval
    if (HAUKSSON_REGIONAL_MODEL) then
      call hauksson_model(vp_hauksson,vs_hauksson,xmesh,ymesh,zmesh,vp_ref_hauksson,vs_dummy, MOHO_MAP_LUPEI)
    else
      call socal_model(doubling_index,rho_dummy,vp_ref_hauksson,vs_dummy)
    endif
    gamma_interp_x = (utm_x_eval - (END_X_GOCAD_MR - THICKNESS_TAPER_BLOCK_MR)) / THICKNESS_TAPER_BLOCK_MR
    vp_final = vp_ref_hauksson * gamma_interp_x + vp_final * (1. - gamma_interp_x)

! y = ymin
  else if (utm_y_eval < ORIG_Y_GOCAD_MR + THICKNESS_TAPER_BLOCK_MR) then
    xmesh = utm_x_eval
    ymesh = ORIG_Y_GOCAD_MR
    zmesh = z_eval
    if (HAUKSSON_REGIONAL_MODEL) then
      call hauksson_model(vp_hauksson,vs_hauksson,xmesh,ymesh,zmesh,vp_ref_hauksson,vs_dummy, MOHO_MAP_LUPEI)
    else
      call socal_model(doubling_index,rho_dummy,vp_ref_hauksson,vs_dummy)
    endif
    gamma_interp_y = (utm_y_eval - ORIG_Y_GOCAD_MR) / THICKNESS_TAPER_BLOCK_MR
    vp_final = vp_ref_hauksson * (1. - gamma_interp_y) + vp_final * gamma_interp_y

! y = ymax
  else if (utm_y_eval > END_Y_GOCAD_MR - THICKNESS_TAPER_BLOCK_MR) then
    xmesh = utm_x_eval
    ymesh = END_Y_GOCAD_MR
    zmesh = z_eval
    if (HAUKSSON_REGIONAL_MODEL) then
      call hauksson_model(vp_hauksson,vs_hauksson,xmesh,ymesh,zmesh,vp_ref_hauksson,vs_dummy, MOHO_MAP_LUPEI)
    else
      call socal_model(doubling_index,rho_dummy,vp_ref_hauksson,vs_dummy)
    endif
    gamma_interp_y = (utm_y_eval - (END_Y_GOCAD_MR - THICKNESS_TAPER_BLOCK_MR)) / THICKNESS_TAPER_BLOCK_MR
    vp_final = vp_ref_hauksson * gamma_interp_y + vp_final * (1. - gamma_interp_y)

  endif

  endif

! use linear variation of vp/vs ratio with depth, between 0. and 8.5 km
         vp_vs_ratio = VP_VS_RATIO_GOCAD_BOTTOM + &
           (VP_VS_RATIO_GOCAD_TOP - VP_VS_RATIO_GOCAD_BOTTOM) * &
           (z_eval - (-8500.d0)) / (0.d0 - (-8500.d0))

! make sure ratio remains in interval
  if (vp_vs_ratio < VP_VS_RATIO_GOCAD_BOTTOM) vp_vs_ratio = VP_VS_RATIO_GOCAD_BOTTOM
  if (vp_vs_ratio > VP_VS_RATIO_GOCAD_TOP) vp_vs_ratio = VP_VS_RATIO_GOCAD_TOP

         vs_final = vp_final / vp_vs_ratio
         call compute_rho_estimate(rho_final,vp_final)

     endif

  end subroutine interpolate_gocad_block_MR

