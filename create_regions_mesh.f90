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

  subroutine create_regions_mesh(xgrid,ygrid,zgrid,ibool,idoubling, &
           xstore,ystore,zstore,npx,npy,iproc_xi,iproc_eta,nspec, &
           volume_local,area_local_bottom,area_local_top, &
           NGLOB_AB,npointot, &
           NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM,NER, &
           NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
           HARVARD_3D_GOCAD_MODEL,NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI, &
           NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
           myrank,LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,UTM_PROJECTION_ZONE, &
           HAUKSSON_REGIONAL_MODEL,OCEANS, &
           VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
           IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,MOHO_MAP_LUPEI, &
           ANISOTROPY,SAVE_MESH_FILES,SUPPRESS_UTM_PROJECTION, &
           ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO,NX_TOPO,NY_TOPO,USE_REGULAR_MESH)

! create the different regions of the mesh

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer nspec

  integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,UTM_PROJECTION_ZONE
  integer NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM,NER

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI
  integer NSPEC2D_A_ETA,NSPEC2D_B_ETA
  integer NX_TOPO,NY_TOPO

  integer npx,npy
  integer npointot

  logical HARVARD_3D_GOCAD_MODEL,HAUKSSON_REGIONAL_MODEL
  logical OCEANS,IMPOSE_MINIMUM_VP_GOCAD,USE_REGULAR_MESH
  logical MOHO_MAP_LUPEI,ANISOTROPY,SAVE_MESH_FILES,SUPPRESS_UTM_PROJECTION

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM
  double precision horiz_size,vert_size,THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR
  double precision ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO

  character(len=150) LOCAL_PATH

! arrays with the mesh
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  double precision xstore_local(NGLLX,NGLLY,NGLLZ)
  double precision ystore_local(NGLLX,NGLLY,NGLLZ)
  double precision zstore_local(NGLLX,NGLLY,NGLLZ)

  double precision xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
  double precision ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)
  double precision zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA)

  double precision xmesh,ymesh,zmesh

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! use integer array to store topography values
  integer icornerlat,icornerlong
  double precision lat,long,elevation
  double precision long_corner,lat_corner,ratio_xi,ratio_eta
  integer itopo_bathy(NX_TOPO,NY_TOPO)

! auxiliary variables to generate the mesh
  integer ix,iy,iz,ir,ir1,ir2,dir
  integer ix1,ix2,dix,iy1,iy2,diy
  integer iax,iay,iar
  integer isubregion,nsubregions,doubling_index

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

! 2D shape functions and their derivatives
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer idoubling(nspec)

! for model density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: c11store,c12store,c13store,c14store,c15store,c16store,&
    c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store,c36store,c44store,c45store,c46store,&
    c55store,c56store,c66store

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

! boundary locator
  logical, dimension(:,:), allocatable :: iboun

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! mass matrix and bathymetry for ocean load
  integer ix_oceans,iy_oceans,iz_oceans,ispec_oceans
  integer ispec2D_ocean_bottom
  integer nglob_oceans
  double precision xval,yval
  double precision height_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight
  double precision area_local_bottom,area_local_top
  double precision volume_local

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: iglob,locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob,NGLOB_AB
  integer ieoff,ilocnum

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

! ---- Moho Vars here ------
! Moho boundary locator
  integer, dimension(:), allocatable :: ibelm_moho_top, ibelm_moho_bot
  logical, dimension(:), allocatable :: is_moho_top, is_moho_bot

! 2-D jacobian and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_moho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho

! number of elements on the boundaries
  integer nspec_moho_top, nspec_moho_bottom
! ---------------------------

! 2-D jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin,jacobian2D_xmax, &
    jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top

! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

! name of the database file
  character(len=150) prname

! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer i,j,k,ia,ispec,iglobnum,itype_element
  integer iproc_xi,iproc_eta

  double precision rho,vp,vs
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

! for the Harvard 3-D basin model
  double precision vp_block_gocad_MR(0:NX_GOCAD_MR-1,0:NY_GOCAD_MR-1,0:NZ_GOCAD_MR-1)
  double precision vp_block_gocad_HR(0:NX_GOCAD_HR-1,0:NY_GOCAD_HR-1,0:NZ_GOCAD_HR-1)
  integer irecord,nrecord,i_vp
  character(len=150) BASIN_MODEL_3D_MEDIUM_RES_FILE,BASIN_MODEL_3D_HIGH_RES_FILE

! for the harvard 3D salton sea model
  real :: vp_st_gocad(GOCAD_ST_NU,GOCAD_ST_NV,GOCAD_ST_NW)
  double precision :: umesh, vmesh, wmesh, vp_st, vs_st, rho_st

! for Hauksson's model
  double precision, dimension(NLAYERS_HAUKSSON,NGRID_NEW_HAUKSSON,NGRID_NEW_HAUKSSON) :: vp_hauksson,vs_hauksson
  integer ilayer
  character(len=150 ) HAUKSSON_REGIONAL_MODEL_FILE

! Stacey put back
! indices for Clayton-Engquist absorbing conditions
  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! flag indicating whether point is in the sediments
  logical point_is_in_sediments
  logical, dimension(:,:,:,:), allocatable :: flag_sediments
  logical, dimension(:), allocatable :: not_fully_in_bedrock

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

! **************

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX))
  allocate(yigll(NGLLY))
  allocate(zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX))
  allocate(wygll(NGLLY))
  allocate(wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ))
  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))

! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ))
  allocate(shape2D_y(NGNOD2D,NGLLX,NGLLZ))
  allocate(shape2D_bottom(NGNOD2D,NGLLX,NGLLY))
  allocate(shape2D_top(NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ))
  allocate(dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ))
  allocate(dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY))

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mustore(NGLLX,NGLLY,NGLLZ,nspec))

! array with anisotropy
  allocate(c11store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c12store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c13store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c14store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c15store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c16store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c22store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c23store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c24store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c25store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c26store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c33store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c34store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c35store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c36store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c44store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c45store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c46store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c55store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c56store(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(c66store(NGLLX,NGLLY,NGLLZ,nspec))

! Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,nspec))

! flag indicating whether point is in the sediments
  allocate(flag_sediments(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(not_fully_in_bedrock(nspec))

! boundary locator
  allocate(iboun(6,nspec))

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,nspec))

! boundary parameters locator
  allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX))
  allocate(ibelm_xmax(NSPEC2DMAX_XMIN_XMAX))
  allocate(ibelm_ymin(NSPEC2DMAX_YMIN_YMAX))
  allocate(ibelm_ymax(NSPEC2DMAX_YMIN_YMAX))
  allocate(ibelm_bottom(NSPEC2D_BOTTOM))
  allocate(ibelm_top(NSPEC2D_TOP))

! 2-D jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP))

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP))

! Moho boundary parameters, 2-D jacobians and normals
  if (SAVE_MOHO_MESH) then
    allocate(ibelm_moho_top(NSPEC2D_BOTTOM))
    allocate(ibelm_moho_bot(NSPEC2D_BOTTOM))
    allocate(is_moho_top(nspec))
    allocate(is_moho_bot(nspec))
    is_moho_top = .false.
    is_moho_bot = .false.
    nspec_moho_top = 0
    nspec_moho_bottom = 0
    allocate(jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_BOTTOM))
    allocate(normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM))
  endif


! Stacey put back
  allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX))
  allocate(nimax(2,NSPEC2DMAX_YMIN_YMAX))
  allocate(njmin(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(njmax(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX))

! MPI cut-planes parameters along xi and along eta
  allocate(iMPIcut_xi(2,nspec))
  allocate(iMPIcut_eta(2,nspec))

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! allocate memory for arrays
  allocate(iglob(npointot))
  allocate(locval(npointot))
  allocate(ifseg(npointot))
  allocate(xp(npointot))
  allocate(yp(npointot))
  allocate(zp(npointot))

!--- read Hauksson's model
  if(HAUKSSON_REGIONAL_MODEL) then
    call get_value_string(HAUKSSON_REGIONAL_MODEL_FILE, &
                          'model.HAUKSSON_REGIONAL_MODEL_FILE', &
                          'DATA/hauksson_model/hauksson_final_grid_smooth.dat')
!    call get_value_string(HAUKSSON_REGIONAL_MODEL_FILE, &
!                          'model.HAUKSSON_REGIONAL_MODEL_FILE', &
!                          'DATA/lin_model/lin_final_grid_smooth.dat')
    open(unit=14,file=HAUKSSON_REGIONAL_MODEL_FILE,status='old',action='read')
    do iy = 1,NGRID_NEW_HAUKSSON
      do ix = 1,NGRID_NEW_HAUKSSON
        read(14,*) (vp_hauksson(ilayer,ix,iy),ilayer=1,NLAYERS_HAUKSSON), &
                   (vs_hauksson(ilayer,ix,iy),ilayer=1,NLAYERS_HAUKSSON)
      enddo
    enddo
    close(14)
    vp_hauksson(:,:,:) = vp_hauksson(:,:,:) * 1000.d0
    vs_hauksson(:,:,:) = vs_hauksson(:,:,:) * 1000.d0
  endif

!--- read the Harvard 3-D basin model
  if(HARVARD_3D_GOCAD_MODEL) then

! read medium-resolution model

! initialize array to undefined values everywhere
  vp_block_gocad_MR(:,:,:) = 20000.

! read Vp from extracted text file
  call get_value_string(BASIN_MODEL_3D_MEDIUM_RES_FILE, &
                        'model.BASIN_MODEL_3D_MEDIUM_RES_FILE', &
                        'DATA/la_3D_block_harvard/la_3D_medium_res/LA_MR_voxet_extracted.txt')
  open(unit=27,file=BASIN_MODEL_3D_MEDIUM_RES_FILE,status='old',action='read')
  read(27,*) nrecord
  do irecord = 1,nrecord
    read(27,*) ix,iy,iz,i_vp
    if(ix<0 .or. ix>NX_GOCAD_MR-1 .or. iy<0 .or. iy>NY_GOCAD_MR-1 .or. iz<0 .or. iz>NZ_GOCAD_MR-1) &
      stop 'wrong array index read in Gocad medium-resolution file'
    vp_block_gocad_MR(ix,iy,iz) = dble(i_vp)
  enddo
  close(27)

! read high-resolution model

! initialize array to undefined values everywhere
  vp_block_gocad_HR(:,:,:) = 20000.

! read Vp from extracted text file
  call get_value_string(BASIN_MODEL_3D_HIGH_RES_FILE, &
                        'model.BASIN_MODEL_3D_HIGH_RES_FILE', &
                        'DATA/la_3D_block_harvard/la_3D_high_res/LA_HR_voxet_extracted.txt')
  open(unit=27,file=BASIN_MODEL_3D_HIGH_RES_FILE,status='old',action='read')
  read(27,*) nrecord
  do irecord = 1,nrecord
    read(27,*) ix,iy,iz,i_vp
    if(ix<0 .or. ix>NX_GOCAD_HR-1 .or. iy<0 .or. iy>NY_GOCAD_HR-1 .or. iz<0 .or. iz>NZ_GOCAD_HR-1) &
      stop 'wrong array index read in Gocad high-resolution file'
    vp_block_gocad_HR(ix,iy,iz) = dble(i_vp)
  enddo
  close(27)

! read Salton Trough model
  call read_salton_sea_model(vp_st_gocad)

  endif

!--- apply heuristic rule to modify doubling regions to balance angles

  if(APPLY_HEURISTIC_RULE .and. .not. USE_REGULAR_MESH) then

! define number of subregions affected by heuristic rule in doubling regions
  nsubregions = 8

  do isubregion = 1,nsubregions

! define shape of elements for heuristic
    call define_subregions_heuristic(myrank,isubregion,iaddx,iaddy,iaddz, &
              ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
              itype_element,npx,npy, &
              NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM)

! loop on all the mesh points in current subregion
  do ir = ir1,ir2,dir
    do iy = iy1,iy2,diy
      do ix = ix1,ix2,dix

! this heuristic rule is only valid for 8-node elements
! it would not work in the case of 27 nodes

!----
    if(itype_element == ITYPE_UNUSUAL_1) then

! side 1
      horiz_size = xgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) &
                 - xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      xgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
         xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * MAGIC_RATIO

      vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
                 - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
         zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

! side 2
      horiz_size = xgrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
                 - xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
      xgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
         xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + horiz_size * MAGIC_RATIO

      vert_size = zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) &
                 - zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
      zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
         zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + vert_size * MAGIC_RATIO / 0.50

!----
    else if(itype_element == ITYPE_UNUSUAL_1p) then

! side 1
      horiz_size = xgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) &
                 - xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      xgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
         xgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * (1. - MAGIC_RATIO)

      vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
                 - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
         zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

! side 2
      horiz_size = xgrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
                 - xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
      xgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
         xgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + horiz_size * (1. - MAGIC_RATIO)

      vert_size = zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) &
                 - zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4))
      zgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
         zgrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) + vert_size * MAGIC_RATIO / 0.50

!----
    else if(itype_element == ITYPE_UNUSUAL_4) then

! side 1
      horiz_size = ygrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
                 - ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
      ygrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
         ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + horiz_size * (1. - MAGIC_RATIO)

      vert_size = zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) &
                 - zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
      zgrid(ir+iar*iaddz(7),ix+iax*iaddx(7),iy+iay*iaddy(7)) = &
         zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + vert_size * MAGIC_RATIO / 0.50

! side 2
      horiz_size = ygrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) &
                 - ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      ygrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
         ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * (1. - MAGIC_RATIO)

      vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
                 - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      zgrid(ir+iar*iaddz(8),ix+iax*iaddx(8),iy+iay*iaddy(8)) = &
         zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

!----
    else if(itype_element == ITYPE_UNUSUAL_4p) then

! side 1
      horiz_size = ygrid(ir+iar*iaddz(3),ix+iax*iaddx(3),iy+iay*iaddy(3)) &
                 - ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
      ygrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
         ygrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + horiz_size * MAGIC_RATIO

      vert_size = zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) &
                 - zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2))
      zgrid(ir+iar*iaddz(6),ix+iax*iaddx(6),iy+iay*iaddy(6)) = &
         zgrid(ir+iar*iaddz(2),ix+iax*iaddx(2),iy+iay*iaddy(2)) + vert_size * MAGIC_RATIO / 0.50

! side 2
      horiz_size = ygrid(ir+iar*iaddz(4),ix+iax*iaddx(4),iy+iay*iaddy(4)) &
                 - ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      ygrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
         ygrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + horiz_size * MAGIC_RATIO

      vert_size = zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) &
                 - zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1))
      zgrid(ir+iar*iaddz(5),ix+iax*iaddx(5),iy+iay*iaddy(5)) = &
         zgrid(ir+iar*iaddz(1),ix+iax*iaddx(1),iy+iay*iaddy(1)) + vert_size * MAGIC_RATIO / 0.50

    endif

      enddo
    enddo
  enddo

  enddo

  endif

!---

! generate the elements in all the regions of the mesh
  ispec = 0

! define number of subregions in the mesh
  if(USE_REGULAR_MESH) then
    nsubregions = 2
  else
    if(NER_SEDIM > 1) then
      nsubregions = 30
    else
      nsubregions = 29
    endif
  endif

  do isubregion = 1,nsubregions

! define shape of elements
    call define_subregions(myrank,isubregion,iaddx,iaddy,iaddz, &
              ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
              doubling_index,npx,npy, &
              NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM,NER,USE_REGULAR_MESH)

! loop on all the mesh points in current subregion
  do ir = ir1,ir2,dir
    do iy = iy1,iy2,diy
      do ix = ix1,ix2,dix

!       loop over the NGNOD nodes
        do ia=1,NGNOD
          xelm(ia) = xgrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
          yelm(ia) = ygrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
          zelm(ia) = zgrid(ir+iar*iaddz(ia),ix+iax*iaddx(ia),iy+iay*iaddy(ia))
        enddo

! add one spectral element to the list and store its material number
        ispec = ispec + 1
        if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
        idoubling(ispec) = doubling_index

! assign Moho surface element
        if (SAVE_MOHO_MESH) then
        if (isubregion == 15 .and. ir == ir1) then
          nspec_moho_top = nspec_moho_top + 1
          if (nspec_moho_top > NSPEC2D_BOTTOM) call exit_mpi(myrank,"Error counting moho top elements")
          ibelm_moho_top(nspec_moho_top) = ispec
          call compute_jacobian_2D(myrank,nspec_moho_top,xelm(1:NGNOD2D),yelm(1:NGNOD2D),zelm(1:NGNOD2D), &
                     dershape2D_bottom,jacobian2D_moho,normal_moho,NGLLX,NGLLY,NSPEC2D_BOTTOM)
          is_moho_top(ispec) = .true.
        else if (isubregion == 28 .and. ir+dir > ir2) then
          nspec_moho_bottom = nspec_moho_bottom + 1
          if (nspec_moho_bottom > NSPEC2D_BOTTOM) call exit_mpi(myrank,"Error counting moho bottom elements")
          ibelm_moho_bot(nspec_moho_bottom) = ispec
          is_moho_bot(ispec) = .true.
        endif
        endif

! initialize flag indicating whether element is in sediments
  not_fully_in_bedrock(ispec) = .false.

! create mesh element
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

! compute mesh coordinates
       xmesh = ZERO
       ymesh = ZERO
       zmesh = ZERO
       do ia=1,NGNOD
         xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
         ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
         zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
       enddo

       xstore_local(i,j,k) = xmesh
       ystore_local(i,j,k) = ymesh
       zstore_local(i,j,k) = zmesh

! initialize flag indicating whether point is in the sediments
       point_is_in_sediments = .false.

       if(ANISOTROPY) then
          call aniso_model(doubling_index,zmesh,rho,vp,vs,c11,c12,c13,c14,c15,c16,&
               c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
       else
! get the regional model parameters
          if(HAUKSSON_REGIONAL_MODEL) then
! get density from socal model
             call socal_model(doubling_index,rho,vp,vs)
! get vp and vs from Hauksson
             call hauksson_model(vp_hauksson,vs_hauksson,xmesh,ymesh,zmesh,vp,vs,MOHO_MAP_LUPEI)
! if Moho map is used, then assume homogeneous medium below the Moho
! and use bottom layer of Hauksson's model in the halfspace
             if(MOHO_MAP_LUPEI .and. doubling_index == IFLAG_HALFSPACE_MOHO) &
                  call socal_model(IFLAG_HALFSPACE_MOHO,rho,vp,vs)
          else
             call socal_model(doubling_index,rho,vp,vs)
! include attenuation in first SoCal layer if needed
! uncomment line below to include attenuation in the 1D case
!        if(zmesh >= DEPTH_5p5km_SOCAL) point_is_in_sediments = .true.
          endif

! get the Harvard 3-D basin model
          if(HARVARD_3D_GOCAD_MODEL .and. &
               (doubling_index == IFLAG_ONE_LAYER_TOPOGRAPHY &
               .or. doubling_index == IFLAG_BASEMENT_TOPO) &
               .and. xmesh >= ORIG_X_GOCAD_MR &
               .and. xmesh <= END_X_GOCAD_MR &
               .and. ymesh >= ORIG_Y_GOCAD_MR &
               .and. ymesh <= END_Y_GOCAD_MR) then

! use medium-resolution model first
             call interpolate_gocad_block_MR(vp_block_gocad_MR, &
                  xmesh,ymesh,zmesh,rho,vp,vs,point_is_in_sediments, &
                  VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
                  IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_MR, &
                  vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL,&
                  MOHO_MAP_LUPEI)

! then superimpose high-resolution model
             if(xmesh >= ORIG_X_GOCAD_HR &
                  .and. xmesh <= END_X_GOCAD_HR &
                  .and. ymesh >= ORIG_Y_GOCAD_HR &
                  .and. ymesh <= END_Y_GOCAD_HR) &
                  call interpolate_gocad_block_HR(vp_block_gocad_HR,vp_block_gocad_MR,&
                  xmesh,ymesh,zmesh,rho,vp,vs,point_is_in_sediments, &
                  VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
                  IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR, &
                  vp_hauksson,vs_hauksson,doubling_index,HAUKSSON_REGIONAL_MODEL, &
                  MOHO_MAP_LUPEI)

          endif
! get the Harvard Salton Trough model
          if (HARVARD_3D_GOCAD_MODEL) then
            call vx_xyz2uvw(xmesh, ymesh, zmesh, umesh, vmesh, wmesh)
            if (umesh >= 0 .and. umesh <= GOCAD_ST_NU-1 .and. &
                  vmesh >= 0 .and. vmesh <=  GOCAD_ST_NV-1 .and. &
                  wmesh >= 0 .and. wmesh <= GOCAD_ST_NW-1) then
              call vx_xyz_interp(umesh,vmesh,wmesh, vp_st, vs_st, rho_st, vp_st_gocad)
              if (abs(vp_st - GOCAD_ST_NO_DATA_VALUE) > 1.0d-3) then
                vp = vp_st
                vs = vs_st
                rho = rho_st
              endif
            endif
          endif
       endif
! store flag indicating whether point is in the sediments
  flag_sediments(i,j,k,ispec) = point_is_in_sediments
  if(point_is_in_sediments) not_fully_in_bedrock(ispec) = .true.

! define elastic parameters in the model
! distinguish between single and double precision for reals
  if(ANISOTROPY) then

       if(CUSTOM_REAL == SIZE_REAL) then
         rhostore(i,j,k,ispec) = sngl(rho)
         kappastore(i,j,k,ispec) = sngl(rho*(vp*vp - 4.d0*vs*vs/3.d0))
         mustore(i,j,k,ispec) = sngl(rho*vs*vs)
         c11store(i,j,k,ispec) = sngl(c11)
         c12store(i,j,k,ispec) = sngl(c12)
         c13store(i,j,k,ispec) = sngl(c13)
         c14store(i,j,k,ispec) = sngl(c14)
         c15store(i,j,k,ispec) = sngl(c15)
         c16store(i,j,k,ispec) = sngl(c16)
         c22store(i,j,k,ispec) = sngl(c22)
         c23store(i,j,k,ispec) = sngl(c23)
         c24store(i,j,k,ispec) = sngl(c24)
         c25store(i,j,k,ispec) = sngl(c25)
         c26store(i,j,k,ispec) = sngl(c26)
         c33store(i,j,k,ispec) = sngl(c33)
         c34store(i,j,k,ispec) = sngl(c34)
         c35store(i,j,k,ispec) = sngl(c35)
         c36store(i,j,k,ispec) = sngl(c36)
         c44store(i,j,k,ispec) = sngl(c44)
         c45store(i,j,k,ispec) = sngl(c45)
         c46store(i,j,k,ispec) = sngl(c46)
         c55store(i,j,k,ispec) = sngl(c55)
         c56store(i,j,k,ispec) = sngl(c56)
         c66store(i,j,k,ispec) = sngl(c66)
! Stacey
         rho_vp(i,j,k,ispec) = sngl(rho*vp)
         rho_vs(i,j,k,ispec) = sngl(rho*vs)
      else
         rhostore(i,j,k,ispec) = rho
         kappastore(i,j,k,ispec) = rho*(vp*vp - 4.d0*vs*vs/3.d0)
         mustore(i,j,k,ispec) = rho*vs*vs
         c11store(i,j,k,ispec) = c11
         c12store(i,j,k,ispec) = c12
         c13store(i,j,k,ispec) = c13
         c14store(i,j,k,ispec) = c14
         c15store(i,j,k,ispec) = c15
         c16store(i,j,k,ispec) = c16
         c22store(i,j,k,ispec) = c22
         c23store(i,j,k,ispec) = c23
         c24store(i,j,k,ispec) = c24
         c25store(i,j,k,ispec) = c25
         c26store(i,j,k,ispec) = c26
         c33store(i,j,k,ispec) = c33
         c34store(i,j,k,ispec) = c34
         c35store(i,j,k,ispec) = c35
         c36store(i,j,k,ispec) = c36
         c44store(i,j,k,ispec) = c44
         c45store(i,j,k,ispec) = c45
         c46store(i,j,k,ispec) = c46
         c55store(i,j,k,ispec) = c55
         c56store(i,j,k,ispec) = c56
         c66store(i,j,k,ispec) = c66
! Stacey
         rho_vp(i,j,k,ispec) = rho*vp
         rho_vs(i,j,k,ispec) = rho*vs
      endif


   else
      if(CUSTOM_REAL == SIZE_REAL) then
         rhostore(i,j,k,ispec) = sngl(rho)
         kappastore(i,j,k,ispec) = sngl(rho*(vp*vp - 4.d0*vs*vs/3.d0))
         mustore(i,j,k,ispec) = sngl(rho*vs*vs)

! Stacey
         rho_vp(i,j,k,ispec) = sngl(rho*vp)
         rho_vs(i,j,k,ispec) = sngl(rho*vs)
      else
         rhostore(i,j,k,ispec) = rho
         kappastore(i,j,k,ispec) = rho*(vp*vp - 4.d0*vs*vs/3.d0)
         mustore(i,j,k,ispec) = rho*vs*vs

! Stacey
         rho_vp(i,j,k,ispec) = rho*vp
         rho_vs(i,j,k,ispec) = rho*vs
      endif
   endif

enddo
enddo
enddo

! detect mesh boundaries
  call get_flags_boundaries(nspec,iproc_xi,iproc_eta,ispec,doubling_index, &
        xstore_local,ystore_local,zstore_local, &
        iboun,iMPIcut_xi,iMPIcut_eta,NPROC_XI,NPROC_ETA, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK)

! compute coordinates and jacobian
        call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
               etaxstore,etaystore,etazstore, &
               gammaxstore,gammaystore,gammazstore,jacobianstore, &
               xstore,ystore,zstore, &
               xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

! end of loop on all the mesh points in current subregion
      enddo
    enddo
  enddo

! end of loop on all the subregions of the current region the mesh
  enddo

! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')
  if (SAVE_MOHO_MESH) then
    if (nspec_moho_top /= NSPEC2D_BOTTOM .or. nspec_moho_bottom /= NSPEC2D_BOTTOM) &
               call exit_mpi(myrank, "nspec_moho should equal NSPEC2D_BOTTOM")
  endif


  do ispec=1,nspec
  ieoff = NGLLCUBE*(ispec-1)
  ilocnum = 0
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ilocnum = ilocnum + 1
        xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
        yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
        zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
      enddo
    enddo
  enddo
  enddo

  call get_global(nspec,xp,yp,zp,iglob,locval,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

! put in classical format
  do ispec=1,nspec
  ieoff = NGLLCUBE*(ispec-1)
  ilocnum = 0
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ilocnum = ilocnum + 1
        ibool(i,j,k,ispec) = iglob(ilocnum+ieoff)
      enddo
    enddo
  enddo
  enddo

  if(minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_AB) &
    call exit_MPI(myrank,'incorrect global numbering')

! create a new indirect addressing array instead, to reduce cache misses
! in memory access in the solver
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mask_ibool(nglob))
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  inumber = 0
  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,k,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
        endif
      enddo
    enddo
  enddo
  enddo
  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

! creating mass matrix (will be fully assembled with MPI in the solver)
  allocate(rmass(nglob))
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        weight=wxgll(i)*wygll(j)*wzgll(k)
        iglobnum=ibool(i,j,k,ispec)

        jacobianl=jacobianstore(i,j,k,ispec)

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight)
    else
      rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
    endif

      enddo
    enddo
  enddo
  enddo

  call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
      dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
              jacobian2D_xmin,jacobian2D_xmax, &
              jacobian2D_ymin,jacobian2D_ymax, &
              jacobian2D_bottom,jacobian2D_top, &
              normal_xmin,normal_xmax, &
              normal_ymin,normal_ymax, &
              normal_bottom,normal_top, &
              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

! create MPI buffers
! arrays locval(npointot) and ifseg(npointot) used to save memory
  call get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC2D_A_ETA,NSPEC2D_B_ETA)
  call get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC2D_A_XI,NSPEC2D_B_XI)

! Stacey put back
  call get_absorb(myrank,prname,iboun,nspec, &
       nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! create AVS or DX mesh data for the slice, edges and faces
  if(SAVE_MESH_FILES) then
    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
  endif

! create ocean load mass matrix
  if(OCEANS) then

! adding ocean load mass matrix at ocean bottom
  nglob_oceans = nglob
  allocate(rmass_ocean_load(nglob_oceans))

! create ocean load mass matrix for degrees of freedom at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

! add contribution of the oceans for surface elements exactly at ocean bottom
    do ispec2D_ocean_bottom = 1,NSPEC2D_TOP

      ispec_oceans = ibelm_top(ispec2D_ocean_bottom)

      iz_oceans = NGLLZ

      do ix_oceans = 1,NGLLX
        do iy_oceans = 1,NGLLY

        iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! compute local height of oceans

! get coordinates of current point
          xval = xstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          yval = ystore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! project x and y in UTM back to long/lat since topo file is in long/lat
  call utm_geo(long,lat,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

! get coordinate of corner in bathy/topo model
    icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
    icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! avoid edge effects and extend with identical point if outside model
    if(icornerlong < 1) icornerlong = 1
    if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
    if(icornerlat < 1) icornerlat = 1
    if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! compute coordinates of corner
    long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
    lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! compute ratio for interpolation
    ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
    ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

! avoid edge effects
    if(ratio_xi < 0.) ratio_xi = 0.
    if(ratio_xi > 1.) ratio_xi = 1.
    if(ratio_eta < 0.) ratio_eta = 0.
    if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
    elevation = &
      itopo_bathy(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
      itopo_bathy(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
      itopo_bathy(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
      itopo_bathy(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

! suppress positive elevation, which means no oceans
    if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
      height_oceans = 0.d0
    else
      height_oceans = dabs(elevation)
    endif

! take into account inertia of water column
        weight = wxgll(ix_oceans)*wygll(iy_oceans)*dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_ocean_bottom)) &
                   * dble(RHO_OCEANS) * height_oceans

! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
        else
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
        endif

        enddo
      enddo

    enddo

! add regular mass matrix to ocean load contribution
  rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

! allocate dummy array if no oceans
    nglob_oceans = 1
    allocate(rmass_ocean_load(nglob_oceans))

  endif

! save the binary files
  call save_arrays_solver(flag_sediments,not_fully_in_bedrock,rho_vp,rho_vs,prname,xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore, &
            gammaxstore,gammaystore,gammazstore,jacobianstore, &
            xstore,ystore,zstore,kappastore,mustore, &
            ANISOTROPY, &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store,c36store, &
            c44store,c45store,c46store,c55store,c56store,c66store, &
            ibool,idoubling,rmass,rmass_ocean_load,nglob_oceans, &
            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
            jacobian2D_xmin,jacobian2D_xmax, &
            jacobian2D_ymin,jacobian2D_ymax, &
            jacobian2D_bottom,jacobian2D_top, &
            iMPIcut_xi,iMPIcut_eta,nspec,nglob, &
            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP,OCEANS)

! save Moho mesh arrays
  if (SAVE_MOHO_MESH) then
    open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='unknown',form='unformatted')
    ! total number of elements, corner points, all points
    write(27) NSPEC2D_BOTTOM
    write(27) (NEX_PER_PROC_XI/4 + 1) * (NEX_PER_PROC_ETA/4 + 1)
    write(27) (NEX_PER_PROC_XI/4 * (NGLLX - 1) + 1) * (NEX_PER_PROC_ETA/4 * (NGLLY - 1) + 1)
    write(27) ibelm_moho_top
    write(27) ibelm_moho_bot
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='unknown',form='unformatted')
    write(27) normal_moho
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='unknown',form='unformatted')
    write(27) is_moho_top
    write(27) is_moho_bot
    close(27)
  endif


  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          weight=wxgll(i)*wygll(j)*wzgll(k)
          jacobianl=jacobianstore(i,j,k,ispec)
          volume_local = volume_local + dble(jacobianl)*weight
        enddo
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_BOTTOM
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_TOP
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_top = area_local_top + dble(jacobian2D_top(i,j,ispec))*weight
      enddo
    enddo
  enddo

  end subroutine create_regions_mesh

!
!----
!

  subroutine create_regions_mesh_ext_mesh(ibool, &
           xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
           nnodes_ext_mesh,nelmnts_ext_mesh, &
           nodes_coords_ext_mesh,elmnts_ext_mesh,mat_ext_mesh,max_static_memory_size, &
           ninterface_ext_mesh,max_interface_size_ext_mesh, &
           my_neighbours_ext_mesh,my_nelmnts_neighbours_ext_mesh,my_interfaces_ext_mesh, &
           ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh)

! create the different regions of the mesh

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer nspec

  integer npointot

  character(len=150) LOCAL_PATH

! arrays with the mesh
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! static memory size needed by the solver
  double precision :: static_memory_size,max_static_memory_size 

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh
  integer, dimension(nelmnts_ext_mesh) :: mat_ext_mesh
  double precision, external :: materials_ext_mesh
  integer :: ninterface_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,ninterface_ext_mesh) :: my_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,ninterface_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: nibool_interfaces_ext_mesh

! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh,ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh_true
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer, dimension(:), allocatable :: ibool_interface_ext_mesh_dummy
  double precision, dimension(:), allocatable :: work_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! for model density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: iglob,locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob
  integer ieoff,ilocnum
  integer ier
  integer iinterface

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! ---------------------------

! name of the database file
  character(len=150) prname

  integer i,j,k,ia,ispec,iglobnum

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber


! **************

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX))
  allocate(yigll(NGLLY))
  allocate(zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX))
  allocate(wygll(NGLLY))
  allocate(wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ))
  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(mustore(NGLLX,NGLLY,NGLLZ,nspec))

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,nspec))

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! allocate memory for arrays
  allocate(iglob(npointot))
  allocate(locval(npointot))
  allocate(ifseg(npointot))
  allocate(xp(npointot))
  allocate(yp(npointot))
  allocate(zp(npointot))

!---

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  do ispec = 1, nspec
     !call get_xyzelm(xelm, yelm, zelm, ispec, elmnts_ext_mesh, nodes_coords_ext_mesh, nspec, nnodes_ext_mesh)
     do ia = 1,NGNOD
     xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
     yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
     zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
     enddo

     call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
          etaxstore,etaystore,etazstore, &
          gammaxstore,gammaystore,gammazstore,jacobianstore, &
          xstore,ystore,zstore, &
          xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

  enddo

! kappastore and mustore
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          kappastore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
               (materials_ext_mesh(2,mat_ext_mesh(ispec))*materials_ext_mesh(2,mat_ext_mesh(ispec)) - &
               4.d0*materials_ext_mesh(3,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))/3.d0)
          mustore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))*&
               materials_ext_mesh(3,mat_ext_mesh(ispec))
        enddo
      enddo
    enddo
  enddo

  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

  do ispec=1,nspec
  ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
  ilocnum = 0
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ilocnum = ilocnum + 1
        xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
        yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
        zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
      enddo
    enddo
  enddo
  enddo

  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot, &
       minval(nodes_coords_ext_mesh(1,:)),maxval(nodes_coords_ext_mesh(1,:)))

  deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'

!
!- we can create a new indirect addressing to reduce cache misses
!
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(mask_ibool(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  inumber = 0
  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,k,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
        endif
      enddo
    enddo
  enddo
  enddo

  deallocate(copy_ibool_ori,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(mask_ibool,stat=ier); if(ier /= 0) stop 'error in deallocate'

  allocate(xstore_dummy(nglob))
  allocate(ystore_dummy(nglob))
  allocate(zstore_dummy(nglob))
  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              xstore_dummy(iglobnum) = xstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              ystore_dummy(iglobnum) = ystore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              zstore_dummy(iglobnum) = zstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

! creating mass matrix (will be fully assembled with MPI in the solver)
  allocate(rmass(nglob))
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        weight=wxgll(i)*wygll(j)*wzgll(k)
        iglobnum=ibool(i,j,k,ispec)

        jacobianl=jacobianstore(i,j,k,ispec)

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(materials_ext_mesh(1,mat_ext_mesh(ispec))) * dble(jacobianl) * weight)
    else
      rmass(iglobnum) = rmass(iglobnum) + materials_ext_mesh(1,mat_ext_mesh(ispec)) * jacobianl * weight
    endif

      enddo
    enddo
  enddo
  enddo

  call prepare_assemble_MPI (nelmnts_ext_mesh,ibool, &
       elmnts_ext_mesh, ESIZE, &
       nglob, &
       ninterface_ext_mesh, max_interface_size_ext_mesh, &
       my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
       ibool_interfaces_ext_mesh, &
       nibool_interfaces_ext_mesh &
       )

! sort ibool comm buffers lexicographically
  allocate(nibool_interfaces_ext_mesh_true(ninterface_ext_mesh))

  do iinterface = 1, ninterface_ext_mesh

    allocate(xp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(yp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(zp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(locval(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ifseg(nibool_interfaces_ext_mesh(iinterface)))
    allocate(reorder_interface_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ibool_interface_ext_mesh_dummy(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ind_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ninseg_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(iwork_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(work_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))

    do ilocnum = 1, nibool_interfaces_ext_mesh(iinterface)
      xp(ilocnum) = xstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      yp(ilocnum) = ystore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      zp(ilocnum) = zstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
    enddo

    call sort_array_coordinates(nibool_interfaces_ext_mesh(iinterface),xp,yp,zp, &
         ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         reorder_interface_ext_mesh,locval,ifseg,nibool_interfaces_ext_mesh_true(iinterface), &
         ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh,work_ext_mesh)

    deallocate(xp)
    deallocate(yp)
    deallocate(zp)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(reorder_interface_ext_mesh)
    deallocate(ibool_interface_ext_mesh_dummy)
    deallocate(ind_ext_mesh)
    deallocate(ninseg_ext_mesh)
    deallocate(iwork_ext_mesh)
    deallocate(work_ext_mesh)

  enddo

! save the binary files
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IOUT,file=prname(1:len_trim(prname))//'external_mesh.bin',status='unknown',action='write',form='unformatted')
  write(IOUT) nspec
  write(IOUT) nglob

  write(IOUT) xixstore
  write(IOUT) xiystore
  write(IOUT) xizstore
  write(IOUT) etaxstore
  write(IOUT) etaystore
  write(IOUT) etazstore
  write(IOUT) gammaxstore
  write(IOUT) gammaystore
  write(IOUT) gammazstore

  write(IOUT) jacobianstore
  write(IOUT) kappastore
  write(IOUT) mustore

  write(IOUT) rmass

  write(IOUT) ibool

  write(IOUT) xstore_dummy
  write(IOUT) ystore_dummy
  write(IOUT) zstore_dummy

  write(IOUT) ninterface_ext_mesh
  write(IOUT) maxval(nibool_interfaces_ext_mesh)
  write(IOUT) my_neighbours_ext_mesh
  write(IOUT) nibool_interfaces_ext_mesh
  allocate(ibool_interfaces_ext_mesh_dummy(maxval(nibool_interfaces_ext_mesh),ninterface_ext_mesh))
  do i = 1, ninterface_ext_mesh
     ibool_interfaces_ext_mesh_dummy = ibool_interfaces_ext_mesh(1:maxval(nibool_interfaces_ext_mesh),:)
  enddo
  write(IOUT) ibool_interfaces_ext_mesh_dummy
  close(IOUT)

! compute the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),ninterface_ext_mesh,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)

  end subroutine create_regions_mesh_ext_mesh

!
!----
!

  double precision function materials_ext_mesh(i,j)

    implicit none

    integer :: i,j

    select case (j)
      case (1)
        select case (i)
          case (1)
            materials_ext_mesh = 2700.d0
          case (2)
            materials_ext_mesh = 3000.d0
          case (3)
            materials_ext_mesh = 1732.051d0
          case default
            call stop_all()
          end select
      case (2)
        select case (i)
          case (1)
            materials_ext_mesh = 2000.d0
          case (2)
            materials_ext_mesh = 900.d0
          case (3)
            materials_ext_mesh = 500.d0
          case default
            call stop_all()
          end select
      case default
        call stop_all()
    end select

  end function materials_ext_mesh

!
!----
!

subroutine prepare_assemble_MPI (nelmnts,ibool, &
     knods, ngnode, &
     npoin, &
     ninterface, max_interface_size, &
     my_nelmnts_neighbours, my_interfaces, &
     ibool_interfaces_asteroid, &
     nibool_interfaces_asteroid &
     )

  implicit none

  include 'constants.h'

  integer, intent(in)  :: nelmnts, npoin, ngnode
  integer, dimension(ngnode,nelmnts), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLY,NGLLZ,nelmnts), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours
  integer, dimension(6,max_interface_size,ninterface)  :: my_interfaces
  integer, dimension(NGLLX*NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_asteroid
  integer, dimension(ninterface)  :: &
       nibool_interfaces_asteroid

  integer  :: num_interface
  integer  :: ispec_interface

  logical, dimension(npoin)  :: mask_ibool_asteroid

  integer  :: ixmin, ixmax
  integer  :: iymin, iymax
  integer  :: izmin, izmax
  integer, dimension(ngnode)  :: n
  integer  :: e1, e2, e3, e4
  integer  :: type
  integer  :: ispec

  integer  :: k
  integer  :: npoin_interface_asteroid

  integer  :: ix,iy,iz


  ibool_interfaces_asteroid(:,:) = 0
  nibool_interfaces_asteroid(:) = 0

  do num_interface = 1, ninterface
     npoin_interface_asteroid = 0
     mask_ibool_asteroid(:) = .false.

     do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
        ispec = my_interfaces(1,ispec_interface,num_interface)
        type = my_interfaces(2,ispec_interface,num_interface)
        do k = 1, ngnode
           n(k) = knods(k,ispec)
        end do
        e1 = my_interfaces(3,ispec_interface,num_interface)
        e2 = my_interfaces(4,ispec_interface,num_interface)
        e3 = my_interfaces(5,ispec_interface,num_interface)
        e4 = my_interfaces(6,ispec_interface,num_interface)
        call get_edge(ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax)

        do iz = min(izmin,izmax), max(izmin,izmax)
           do iy = min(iymin,iymax), max(iymin,iymax)
              do ix = min(ixmin,ixmax), max(ixmin,ixmax)

                 if(.not. mask_ibool_asteroid(ibool(ix,iy,iz,ispec))) then
                    mask_ibool_asteroid(ibool(ix,iy,iz,ispec)) = .true.
                    npoin_interface_asteroid = npoin_interface_asteroid + 1
                    ibool_interfaces_asteroid(npoin_interface_asteroid,num_interface)=&
                         ibool(ix,iy,iz,ispec)
                 end if
              end do
           end do
        end do

     end do
     nibool_interfaces_asteroid(num_interface) = npoin_interface_asteroid


  end do

end subroutine prepare_assemble_MPI

!
!----
!

subroutine get_edge ( ngnode, n, type, e1, e2, e3, e4, ixmin, ixmax, iymin, iymax, izmin, izmax )

  implicit none

  include "constants.h"

  integer, intent(in)  :: ngnode
  integer, dimension(ngnode), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2, e3, e4
  integer, intent(out)  :: ixmin, ixmax, iymin, iymax, izmin, izmax

  integer, dimension(4) :: en
  integer :: valence, i

   if ( type == 1 ) then
     if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(5) ) then
        ixmin = 1
        ixmax = 1
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(6) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = 1
        iymax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(7) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(8) ) then
        ixmin = 1
        ixmax = 1
        iymin = NGLLY
        iymax = NGLLY
        izmin = NGLLZ
        izmax = NGLLZ
     end if
  else
     if ( type == 2 ) then
        if ( e1 ==  n(1) ) then
           ixmin = 1
           iymin = 1
           izmin = 1
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(2) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = 1
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if

        end if
        if ( e1 == n(3) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(4) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = 1
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(5) ) then
           ixmin = 1
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(1) ) then
              ixmax = 1
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(6) ) then
           ixmin = NGLLX
           iymin = 1
           izmin = NGLLZ
           if ( e2 == n(2) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = 1
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(7) ) then
           ixmin = NGLLX
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(3) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(8) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           end if
           if ( e2 == n(6) ) then
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           end if
        end if
        if ( e1 == n(8) ) then
           ixmin = 1
           iymin = NGLLY
           izmin = NGLLZ
           if ( e2 == n(4) ) then
              ixmax = 1
              iymax = NGLLY
              izmax = 1
           end if
           if ( e2 == n(5) ) then
              ixmax = 1
              iymax = 1
              izmax = NGLLZ
           end if
           if ( e2 == n(7) ) then
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           end if
        end if

     else
        if (type == 4) then
           en(1) = e1
           en(2) = e2
           en(3) = e3
           en(4) = e4

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = 1
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = 1
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(2)) then
                 valence = valence+1
              endif
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = NGLLX
              iymin = 1
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLZ
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(3)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = NGLLY
              izmin = 1
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(1)) then
                 valence = valence+1
              endif
              if ( en(i) == n(4)) then
                 valence = valence+1
              endif
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = 1
              ixmax = 1
              iymax = NGLLY
              izmax = NGLLZ
           endif

           valence = 0
           do i = 1, 4
              if ( en(i) == n(5)) then
                 valence = valence+1
              endif
              if ( en(i) == n(6)) then
                 valence = valence+1
              endif
              if ( en(i) == n(7)) then
                 valence = valence+1
              endif
              if ( en(i) == n(8)) then
                 valence = valence+1
              endif
           enddo
           if ( valence == 4 ) then
              ixmin = 1
              iymin = 1
              izmin = NGLLZ
              ixmax = NGLLX
              iymax = NGLLY
              izmax = NGLLZ
           endif

        else
           stop 'ERROR get_edge'
        endif

     end if
  end if

end subroutine get_edge
