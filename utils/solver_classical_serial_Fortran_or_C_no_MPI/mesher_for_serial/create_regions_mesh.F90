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

  subroutine create_regions_mesh(iregion_code,ibool,idoubling, &
           xstore,ystore,zstore,rmins,rmaxs, &
           iproc_xi,iproc_eta,ichunk,nspec,nspec_tiso, &
           volume_local,nspl,rspl,espl,espl2, &
           nglob_theor,npointot, &
           NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
           ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           NPROC_XI,NPROC_ETA,NSPEC2D_XI_FACE, &
!!!!           NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,NGLOB2DMAX_XY, &
           NSPEC2D_ETA_FACE,NGLOB2DMAX_XY, &
           myrank,LOCAL_PATH,OCEANS,ibathy_topo, &
           rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,&
           ATTENUATION,ATTENUATION_3D,SAVE_MESH_FILES, &
           NCHUNKS,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           ner,ratio_sampling_array,doubling_index,r_bottom,r_top,this_region_has_a_doubling,CASE_3D, &
           AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ipass,ratio_divide_central_cube, &
           CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
           NSTEP,DT,NPROCTOT,is_on_a_slice_edge)

! create the different regions of the mesh

  implicit none

  include "constants.h"

! standard include of the MPI library
#ifdef USE_MPI
  include 'mpif.h'
#endif

  integer :: NSTEP
  integer, save :: npoin2D_xi,npoin2D_eta
  double precision :: DT

! this to cut the doubling brick
!!!!  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer :: step_mult,offset_proc_xi,offset_proc_eta
  integer :: case_xi,case_eta,subblock_num

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ignod,ner_without_doubling,ispec_superbrick,ilayer,ilayer_loop,ix_elem,iy_elem,iz_elem, &
               ifirst_region,ilast_region,ratio_divide_central_cube
  integer, dimension(:), allocatable :: perm_layer

! mesh doubling superbrick
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick

  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

! aniso_mantle_model_variables
  type aniso_mantle_model_variables
    sequence
    double precision beta(14,34,37,73)
    double precision pro(47)
    integer npar1
  end type aniso_mantle_model_variables

  type (aniso_mantle_model_variables) AMM_V
! aniso_mantle_model_variables

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

  type (attenuation_model_variables) AM_V
! attenuation_model_variables

! model_1066a_variables
  type model_1066a_variables
    sequence
      double precision, dimension(NR_1066A) :: radius_1066a
      double precision, dimension(NR_1066A) :: density_1066a
      double precision, dimension(NR_1066A) :: vp_1066a
      double precision, dimension(NR_1066A) :: vs_1066a
      double precision, dimension(NR_1066A) :: Qkappa_1066a
      double precision, dimension(NR_1066A) :: Qmu_1066a
  end type model_1066a_variables

  type (model_1066a_variables) M1066a_V
! model_1066a_variables

! model_ak135_variables
  type model_ak135_variables
    sequence
    double precision, dimension(NR_AK135) :: radius_ak135
    double precision, dimension(NR_AK135) :: density_ak135
    double precision, dimension(NR_AK135) :: vp_ak135
    double precision, dimension(NR_AK135) :: vs_ak135
    double precision, dimension(NR_AK135) :: Qkappa_ak135
    double precision, dimension(NR_AK135) :: Qmu_ak135
  end type model_ak135_variables

 type (model_ak135_variables) Mak135_V
! model_ak135_variables

! model_ref_variables
  type model_ref_variables
    sequence
     double precision, dimension(NR_REF) :: radius_ref
     double precision, dimension(NR_REF) :: density_ref
     double precision, dimension(NR_REF) :: vpv_ref
     double precision, dimension(NR_REF) :: vph_ref
     double precision, dimension(NR_REF) :: vsv_ref
     double precision, dimension(NR_REF) :: vsh_ref
     double precision, dimension(NR_REF) :: eta_ref
     double precision, dimension(NR_REF) :: Qkappa_ref
     double precision, dimension(NR_REF) :: Qmu_ref
  end type model_ref_variables

 type (model_ref_variables) Mref_V
! model_ref_variables

! sea1d_model_variables
  type sea1d_model_variables
    sequence
     double precision, dimension(NR_SEA1D) :: radius_sea1d
     double precision, dimension(NR_SEA1D) :: density_sea1d
     double precision, dimension(NR_SEA1D) :: vp_sea1d
     double precision, dimension(NR_SEA1D) :: vs_sea1d
     double precision, dimension(NR_SEA1D) :: Qkappa_sea1d
     double precision, dimension(NR_SEA1D) :: Qmu_sea1d
  end type sea1d_model_variables

  type (sea1d_model_variables) SEA1DM_V
! sea1d_model_variables

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables

! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

! sea99_s_model_variables
  type sea99_s_model_variables
    sequence
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
 end type sea99_s_model_variables

 type (sea99_s_model_variables) SEA99M_V
! sea99_s_model_variables

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_tiso,nspec_stacey

  integer NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS,REFERENCE_1D_MODEL,THREE_D_MODEL

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer NPROC_XI,NPROC_ETA,NPROCTOT

  integer npointot

  logical ELLIPTICITY,TOPOGRAPHY,SAVE_MESH_FILES
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,OCEANS

  logical ATTENUATION,ATTENUATION_3D,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS

  double precision R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO, &
          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN

  character(len=150) LOCAL_PATH,errmsg

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! meshing parameters
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! to define the central cube in the inner core
  integer nx_central_cube,ny_central_cube,nz_central_cube
  double precision radius_cube
  double precision :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! auxiliary variables to generate the mesh
  integer ix,iy,iz

! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

! code for the four regions of the mesh
  integer iregion_code

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

! 2D shape functions and their derivatives
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  double precision, dimension(NGNOD) :: xelm,yelm,zelm,offset_x,offset_y,offset_z

  integer idoubling(nspec)

! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision rmin,rmax

! for model density and anisotropy
  integer nspec_ani
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

! boundary locator
  logical, dimension(:,:), allocatable :: iboun

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight
  double precision volume_local

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob,nglob_theor,ieoff,ilocnum,ier

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! mass matrix and bathymetry for ocean load
  integer ix_oceans,iy_oceans,iz_oceans,ispec_oceans
  integer ispec2D_top_crust
  integer nglob_oceans
  double precision xval,yval,zval,rval,thetaval,phival
  double precision lat,lon,colat
  double precision elevation,height_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  logical, dimension(:), allocatable :: mask_ibool2
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

! 2-D jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin,jacobian2D_xmax, &
    jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top

! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

! Stacey, indices for Clayton-Engquist absorbing conditions
  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! name of the database file
  character(len=150) prname

  integer i,j,k,ia,ispec,iglobnum,icolor
  integer iproc_xi,iproc_eta,ichunk

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

! attenuation
  double precision, dimension(:,:,:,:),   allocatable :: Qmu_store
  double precision, dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS)                  :: tau_s
  double precision  T_c_source

! **************
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  logical :: USE_ONE_LAYER_SB,CASE_3D
  integer :: nspec_sb

  integer NUMBER_OF_MESH_LAYERS,layer_shift,first_layer_aniso,last_layer_aniso,FIRST_ELT_NON_ANISO
  double precision, dimension(:,:), allocatable :: stretch_tab

  integer :: NGLOB2DMAX_XY

  integer :: nb_layer_above_aniso,FIRST_ELT_ABOVE_ANISO

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)
  real(kind=4) vercof(maxker)
  real(kind=4) vercofd(maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=80) kerstr
  character(len=40) varstr(maxker)

! now perform two passes in this part to be able to save memory
  integer :: ipass

  logical :: ACTUALLY_STORE_ARRAYS

! added for color permutation
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer
  integer, dimension(:), allocatable :: perm, first_elem_number_in_this_color,number_of_elements_in_this_color
  integer, dimension(:,:,:,:), allocatable :: temp_array_int
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: temp_array_real
  double precision, dimension(:,:,:,:), allocatable :: temp_array_dble
  integer, dimension(:), allocatable :: temp_array_1D_int

! added for inner elements/outer elements splitting of the mesh for the regular C version for a CPU
!!!! DK DK  integer :: inumber_in_new_list_after_perm,nspec_outer_min_global,nspec_outer_max_global
  integer :: nspec_outer_min_global,nspec_outer_max_global

! the height at which the central cube is cut
  integer :: nz_inf_limit

  logical, dimension(nspec) :: is_on_a_slice_edge

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

! Attenuation
  if(ATTENUATION .and. ATTENUATION_3D) then
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
    allocate(Qmu_store(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(tau_e_store(N_SLS,NGLLX,NGLLY,NGLLZ,nspec))
  else
    allocate(Qmu_store(1,1,1,1))
    allocate(tau_e_store(N_SLS,1,1,1,1))
    Qmu_store(1,1,1,1) = 0.0d0
    tau_e_store(:,1,1,1,1) = 0.0d0
  endif

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

! for anisotropy
  allocate(kappavstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(muvstore(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(kappahstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(muhstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec))

! Stacey
  if(NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey))
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey))

  nspec_ani = 1
  if((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
     (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) nspec_ani = nspec

  allocate(c11store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c12store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c13store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c14store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c15store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c16store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c22store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c23store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c24store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c25store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c26store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c33store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c34store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c35store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c36store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c44store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c45store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c46store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c55store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c56store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c66store(NGLLX,NGLLY,NGLLZ,nspec_ani))

! boundary locator
  allocate(iboun(6,nspec))

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

! Stacey
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

! define models 1066a and ak135 and ref
  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    call define_model_1066a(CRUSTAL, M1066a_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
    call define_model_ak135(CRUSTAL, Mak135_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
    call define_model_ref(Mref_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    call define_model_sea1d(CRUSTAL, SEA1DM_V)
  endif

!------------------------------------------------------------------------

! create the shape of the corner nodes of a regular mesh element
  call hex_nodes(iaddx,iaddy,iaddz)

! reference element has size one here, not two
  iaddx(:) = iaddx(:) / 2
  iaddy(:) = iaddy(:) / 2
  iaddz(:) = iaddz(:) / 2

  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
    layer_shift = 0
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
    layer_shift = 1
  endif

  if (.not. ADD_4TH_DOUBLING) NUMBER_OF_MESH_LAYERS = NUMBER_OF_MESH_LAYERS - 1

! define the first and last layers that define this region
  if(iregion_code == IREGION_CRUST_MANTLE) then
    ifirst_region = 1
    ilast_region = 10 + layer_shift

  else if(iregion_code == IREGION_OUTER_CORE) then
    ifirst_region = 11 + layer_shift
    ilast_region = NUMBER_OF_MESH_LAYERS - 1

  else if(iregion_code == IREGION_INNER_CORE) then
    ifirst_region = NUMBER_OF_MESH_LAYERS
    ilast_region = NUMBER_OF_MESH_LAYERS

  else
    call exit_MPI(myrank,'incorrect region code detected')

  endif

! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  if (ONE_CRUST) then
    first_layer_aniso=2
    last_layer_aniso=3
    nb_layer_above_aniso = 1
  else
    first_layer_aniso=3
    last_layer_aniso=4
    nb_layer_above_aniso = 2
  endif
  allocate (perm_layer(ifirst_region:ilast_region))
  perm_layer = (/ (i, i=ilast_region,ifirst_region,-1) /)
! if(iregion_code == IREGION_CRUST_MANTLE) then
!  cpt=3
!  perm_layer(1)=first_layer_aniso
!  perm_layer(2)=last_layer_aniso
!   do i = ilast_region,ifirst_region,-1
!    if (i/=first_layer_aniso .and. i/=last_layer_aniso) then
!       perm_layer(cpt) = i
!      cpt=cpt+1
!    endif
!   enddo
! endif

! initialize mesh arrays
  idoubling(:) = 0

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  if(ipass == 1) ibool(:,:,:,:) = 0

! initialize boundary arrays
  iboun(:,:) = .false.
  iMPIcut_xi(:,:) = .false.
  iMPIcut_eta(:,:) = .false.

! store and save the final arrays only in the second pass
! therefore in the first pass some arrays can be allocated with a dummy size
  if(ipass == 1) then

    ACTUALLY_STORE_ARRAYS = .false.

    allocate(xixstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(xiystore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(xizstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etaxstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etaystore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etazstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammaystore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammazstore(NGLLX,NGLLY,NGLLZ,1),stat=ier); if(ier /= 0) stop 'error in allocate'

  else

    ACTUALLY_STORE_ARRAYS = .true.

    allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'
    allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier); if(ier /= 0) stop 'error in allocate'

  endif

  if (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH) then
    allocate(stretch_tab(2,ner(1)))
    call stretching_function(r_top(1),r_bottom(1),ner(1),stretch_tab)
  endif

! generate and count all the elements in this region of the mesh
  ispec = 0

! loop on all the layers in this region of the mesh
  do ilayer_loop = ifirst_region,ilast_region

    ilayer = perm_layer(ilayer_loop)

! determine the radii that define the shell
  rmin = rmins(ilayer)
  rmax = rmaxs(ilayer)

  if(iregion_code == IREGION_CRUST_MANTLE .and. ilayer_loop==3) then
    FIRST_ELT_NON_ANISO = ispec+1
  endif
  if(iregion_code == IREGION_CRUST_MANTLE .and. ilayer_loop==(ilast_region-nb_layer_above_aniso+1)) then
    FIRST_ELT_ABOVE_ANISO = ispec+1
  endif

    ner_without_doubling = ner(ilayer)

! if there is a doubling at the top of this region, we implement it in the last two layers of elements
! and therefore we suppress two layers of regular elements here
    USE_ONE_LAYER_SB = .false.
    if(this_region_has_a_doubling(ilayer)) then
      if (ner(ilayer) == 1) then
        ner_without_doubling = ner_without_doubling - 1
        USE_ONE_LAYER_SB = .true.
      else
        ner_without_doubling = ner_without_doubling - 2
        USE_ONE_LAYER_SB = .false.
      endif
    endif

!----
!----   regular mesh elements
!----

! loop on all the elements
   do ix_elem = 1,NEX_PER_PROC_XI,ratio_sampling_array(ilayer)
   do iy_elem = 1,NEX_PER_PROC_ETA,ratio_sampling_array(ilayer)
   do iz_elem = 1,ner_without_doubling
! loop on all the corner nodes of this element
   do ignod = 1,NGNOD_EIGHT_CORNERS
! define topological coordinates of this mesh point
      offset_x(ignod) = (ix_elem - 1) + iaddx(ignod) * ratio_sampling_array(ilayer)
      offset_y(ignod) = (iy_elem - 1) + iaddy(ignod) * ratio_sampling_array(ilayer)
      if (ilayer == 1 .and. CASE_3D) then
        offset_z(ignod) = iaddz(ignod)
      else
        offset_z(ignod) = (iz_elem - 1) + iaddz(ignod)
      endif
   enddo
     call add_missing_nodes(offset_x,offset_y,offset_z)

! compute the actual position of all the grid points of that element
  if (ilayer == 1 .and. CASE_3D .and. .not. SUPPRESS_CRUSTAL_MESH) then
! crustal elements are stretched to be thinner in the upper crust than in lower crust in the 3D case
! max ratio between size of upper crust elements and lower crust elements is given by the param MAX_RATIO_STRETCHING
! to avoid stretching, set MAX_RATIO_STRETCHING = 1.0d  in constants.h
    call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               stretch_tab(1,ner_without_doubling-iz_elem+1),&
               stretch_tab(2,ner_without_doubling-iz_elem+1),1,ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
  else
     call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
  endif
! add one spectral element to the list
     ispec = ispec + 1
     if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')

! new get_flag_boundaries
! xmin & xmax
  if (ix_elem == 1) then
      iMPIcut_xi(1,ispec) = .true.
      if (iproc_xi == 0) iboun(1,ispec)= .true.
  endif
  if (ix_elem == (NEX_PER_PROC_XI-ratio_sampling_array(ilayer)+1)) then
      iMPIcut_xi(2,ispec) = .true.
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
  endif
! ymin & ymax
  if (iy_elem == 1) then
      iMPIcut_eta(1,ispec) = .true.
      if (iproc_eta == 0) iboun(3,ispec)= .true.
  endif
  if (iy_elem == (NEX_PER_PROC_ETA-ratio_sampling_array(ilayer)+1)) then
      iMPIcut_eta(2,ispec) = .true.
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
  endif
! zmin & zmax
  if (iz_elem == ner(ilayer) .and. ilayer == ifirst_region) then
      iboun(6,ispec)= .true.
  endif
  if (iz_elem == 1 .and. ilayer == ilast_region) then    ! defined if no doubling in this layer
      iboun(5,ispec)= .true.
  endif

! define the doubling flag of this element
     idoubling(ispec) = doubling_index(ilayer)

! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,rho_vp,rho_vs,&
           AMM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ACTUALLY_STORE_ARRAYS)

! end of loop on all the regular elements
  enddo
  enddo
  enddo
!----
!----   mesh doubling elements
!----
! If there is a doubling at the top of this region, let us add these elements.
! The superbrick implements a symmetric four-to-two doubling and therefore replaces
! a basic regular block of 2 x 2 = 4 elements.
! We have imposed that NEX be a multiple of 16 therefore we know that we can always create
! these 2 x 2 blocks because NEX_PER_PROC_XI / ratio_sampling_array(ilayer) and
! NEX_PER_PROC_ETA / ratio_sampling_array(ilayer) are always divisible by 2.
    if(this_region_has_a_doubling(ilayer)) then
      if (USE_ONE_LAYER_SB) then
        call define_superbrick_one_layer(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
        nspec_sb = NSPEC_SUPERBRICK_1L
        iz_elem = ner(ilayer)
        step_mult = 2
      else
        if(iregion_code==IREGION_OUTER_CORE .and. ilayer==ilast_region .and. (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
          nspec_sb = NSPEC_DOUBLING_BASICBRICK
          step_mult = 1
        else
          call define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
          nspec_sb = NSPEC_DOUBLING_SUPERBRICK
          step_mult = 2
        endif
! the doubling is implemented in the last two radial elements
! therefore we start one element before the last one
        iz_elem = ner(ilayer) - 1
      endif

! loop on all the elements in the 2 x 2 blocks
      do ix_elem = 1,NEX_PER_PROC_XI,step_mult*ratio_sampling_array(ilayer)
        do iy_elem = 1,NEX_PER_PROC_ETA,step_mult*ratio_sampling_array(ilayer)

          if (step_mult == 1) then
! for xi direction
            if (.not. CUT_SUPERBRICK_XI) then
              if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                case_xi = 1
              else
                case_xi = 2
              endif
            else
              if (offset_proc_xi == 0) then
                if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                  case_xi = 1
                else
                  case_xi = 2
                endif
              else
                if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))/=0) then
                  case_xi = 1
                else
                  case_xi = 2
                endif
              endif
            endif
! for eta direction
            if (.not. CUT_SUPERBRICK_ETA) then
              if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                case_eta = 1
              else
                case_eta = 2
              endif
            else
              if (offset_proc_eta == 0) then
                if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                  case_eta = 1
                else
                  case_eta = 2
                endif
              else
                if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))/=0) then
                  case_eta = 1
                else
                  case_eta = 2
                endif
              endif
            endif
! determine the current sub-block
            if (case_xi == 1) then
              if (case_eta == 1) then
                subblock_num = 1
              else
                subblock_num = 2
              endif
            else
              if (case_eta == 1) then
                subblock_num = 3
              else
                subblock_num = 4
              endif
            endif
! then define the geometry for this sub-block
            call define_basic_doubling_brick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb,subblock_num)
          endif
! loop on all the elements in the mesh doubling superbrick
          do ispec_superbrick = 1,nspec_sb
! loop on all the corner nodes of this element
            do ignod = 1,NGNOD_EIGHT_CORNERS

! define topological coordinates of this mesh point
              offset_x(ignod) = (ix_elem - 1) + &
         x_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
              offset_y(ignod) = (iy_elem - 1) + &
         y_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
              offset_z(ignod) = (iz_elem - 1) + &
         z_superbrick(ibool_superbrick(ignod,ispec_superbrick))

            enddo
! the rest of the 27 nodes are missing, therefore add them
     call add_missing_nodes(offset_x,offset_y,offset_z)

! compute the actual position of all the grid points of that element
     call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)

! add one spectral element to the list
     ispec = ispec + 1
     if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')

! new get_flag_boundaries
! xmin & xmax
  if (ix_elem == 1) then
      iMPIcut_xi(1,ispec) = iboun_sb(ispec_superbrick,1)
      if (iproc_xi == 0) iboun(1,ispec)= iboun_sb(ispec_superbrick,1)
  endif
  if (ix_elem == (NEX_PER_PROC_XI-step_mult*ratio_sampling_array(ilayer)+1)) then
      iMPIcut_xi(2,ispec) = iboun_sb(ispec_superbrick,2)
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= iboun_sb(ispec_superbrick,2)
  endif
!! ymin & ymax
  if (iy_elem == 1) then
      iMPIcut_eta(1,ispec) = iboun_sb(ispec_superbrick,3)
      if (iproc_eta == 0) iboun(3,ispec)= iboun_sb(ispec_superbrick,3)
  endif
  if (iy_elem == (NEX_PER_PROC_ETA-step_mult*ratio_sampling_array(ilayer)+1)) then
      iMPIcut_eta(2,ispec) = iboun_sb(ispec_superbrick,4)
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= iboun_sb(ispec_superbrick,4)
  endif
! zmax only
  if (ilayer==ifirst_region) then
    iboun(6,ispec)= iboun_sb(ispec_superbrick,6)
  endif
  if (ilayer==ilast_region .and. iz_elem==1) then
    iboun(5,ispec)= iboun_sb(ispec_superbrick,5)
  endif

! define the doubling flag of this element
     idoubling(ispec) = doubling_index(ilayer)

! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,rho_vp,rho_vs,&
           AMM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ACTUALLY_STORE_ARRAYS)

! end of loops on the mesh doubling elements
          enddo
        enddo
      enddo
    endif

! end of loop on all the layers of the mesh
  enddo

  if (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH) deallocate(stretch_tab)
  deallocate (perm_layer)

!---

! define central cube in inner core

  if(INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then

! create the shape of a regular mesh element in the inner core
  call hex_nodes(iaddx,iaddy,iaddz)

! define vertical slice in central cube on current processor
! we can assume that NEX_XI = NEX_ETA, otherwise central cube cannot be defined
  nx_central_cube = NEX_PER_PROC_XI / ratio_divide_central_cube
  ny_central_cube = NEX_PER_PROC_ETA / ratio_divide_central_cube
  nz_central_cube = NEX_XI / ratio_divide_central_cube

! size of the cube along cartesian axes before rotation
  radius_cube = (R_CENTRAL_CUBE / R_EARTH) / sqrt(3.d0)

! define spectral elements in central cube
  do iz = 0,2*nz_central_cube-2,2
    do iy = 0,2*ny_central_cube-2,2
      do ix = 0,2*nx_central_cube-2,2

!       radii that define the shell, we know that we are in the central cube
        rmin = 0.d0
        rmax = R_CENTRAL_CUBE / R_EARTH

!       loop over the NGNOD nodes
        do ia=1,NGNOD

! flat cubed sphere with correct mapping
          call compute_coord_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia), &
                  xgrid_central_cube,ygrid_central_cube,zgrid_central_cube, &
                  iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube,ny_central_cube,nz_central_cube,radius_cube)

          if(ichunk == CHUNK_AB) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + xgrid_central_cube
            zelm(ia) = + zgrid_central_cube

          else if(ichunk == CHUNK_AB_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - xgrid_central_cube
            zelm(ia) = - zgrid_central_cube

          else if(ichunk == CHUNK_AC) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - zgrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if(ichunk == CHUNK_AC_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + zgrid_central_cube
            zelm(ia) = - xgrid_central_cube

          else if(ichunk == CHUNK_BC) then
            xelm(ia) = - zgrid_central_cube
            yelm(ia) = + ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if(ichunk == CHUNK_BC_ANTIPODE) then
            xelm(ia) = + zgrid_central_cube
            yelm(ia) = - ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else
            call exit_MPI(myrank,'wrong chunk number in flat cubed sphere definition')
          endif

        enddo

! add one spectral element to the list
        ispec = ispec + 1
        if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in central cube creation')

! new get_flag_boundaries
! xmin & xmax
  if (ix == 0) then
      iMPIcut_xi(1,ispec) = .true.
      if (iproc_xi == 0) iboun(1,ispec)= .true.
  endif
  if (ix == 2*nx_central_cube-2) then
      iMPIcut_xi(2,ispec) = .true.
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
  endif
! ymin & ymax
  if (iy == 0) then
      iMPIcut_eta(1,ispec) = .true.
      if (iproc_eta == 0) iboun(3,ispec)= .true.
  endif
  if (iy == 2*ny_central_cube-2) then
      iMPIcut_eta(2,ispec) = .true.
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
  endif

! define the doubling flag of this element
! only two active central cubes, the four others are fictitious

! determine where we cut the central cube to share it between CHUNK_AB & CHUNK_AB_ANTIPODE
! in the case of mod(NPROC_XI,2)/=0, the cut is asymetric and the bigger part is for CHUNK_AB
  if (mod(NPROC_XI,2)/=0) then
    if (ichunk == CHUNK_AB) then
      nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*floor(NPROC_XI/2.d0)
    elseif (ichunk == CHUNK_AB_ANTIPODE) then
      nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*ceiling(NPROC_XI/2.d0)
    endif
  else
    nz_inf_limit = nz_central_cube
  endif

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    if(iz == nz_inf_limit) then
      idoubling(ispec) = IFLAG_BOTTOM_CENTRAL_CUBE
    else if(iz == 2*nz_central_cube-2) then
      idoubling(ispec) = IFLAG_TOP_CENTRAL_CUBE
    else if (iz > nz_inf_limit .and. iz < 2*nz_central_cube-2) then
      idoubling(ispec) = IFLAG_MIDDLE_CENTRAL_CUBE
    else
      idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
    endif
  else
    idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
  endif


! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,rho_vp,rho_vs,&
           AMM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ACTUALLY_STORE_ARRAYS)
      enddo
    enddo
  enddo

  endif    ! end of definition of central cube in inner core

!---

! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

! only create global addressing and the MPI buffers in the first pass
  if(ipass == 1) then

! allocate memory for arrays
  allocate(locval(npointot),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(ifseg(npointot),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(xp(npointot),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(yp(npointot),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(zp(npointot),stat=ier); if(ier /= 0) stop 'error in allocate'

  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

! we need to create a copy of the x, y and z arrays because sorting in get_global will swap
! these arrays and therefore destroy them
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

  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot)

  deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'

! check that number of points found equals theoretical value
  if(nglob /= nglob_theor) then
    write(errmsg,*) 'incorrect total number of points found: myrank,nglob,nglob_theor,ipass,iregion_code = ',&
      myrank,nglob,nglob_theor,ipass,iregion_code
    call exit_MPI(myrank,errmsg)
  endif

  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering')

! create a new indirect addressing to reduce cache misses in memory access in the solver
! this is *critical* to improve performance in the solver
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

  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering after sorting')

  allocate(mask_ibool2(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  mask_ibool2(:) = .false.
! create MPI buffers
! arrays locval(npointot) and ifseg(npointot) used to save memory
  call get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC2D_ETA_FACE,iregion_code,NGLOB2DMAX_XY,mask_ibool2,npoin2D_xi)
  call get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC2D_XI_FACE,iregion_code,NGLOB2DMAX_XY,mask_ibool2,npoin2D_eta)
!! DK DK only one chunk in current MPI+GPU version
! call get_MPI_1D_buffers(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool,idoubling, &
!                 xstore,ystore,zstore,ifseg,npointot, &
!                 NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion_code)

! now label all the elements that have at least one corner belonging
! to any of these buffers as elements that must contribute to the
! first step of the calculations (performed on the edges before starting
! the non-blocking communications); there is no need to examine the inside
! of the elements, checking their eight corners is sufficient
  do ispec = 1,nspec
    do k = 1,NGLLZ,NGLLZ-1
      do j  = 1,NGLLY,NGLLY-1
        do i = 1,NGLLX,NGLLX-1
          if(mask_ibool2(ibool(i,j,k,ispec))) then
            is_on_a_slice_edge(ispec) = .true.
            goto 888
          endif
        enddo
      enddo
    enddo
  888 continue
  enddo

! Stacey
  if(NCHUNKS /= 6) &
!!!!!!!!       call get_absorb(myrank,prname,iboun,nspec,nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
       call get_absorb(myrank,iboun,nspec,nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! create AVS or DX mesh data for the slices
  if(SAVE_MESH_FILES) then
    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
              rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,REFERENCE_1D_MODEL, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN,M1066a_V,Mak135_V,Mref_V,SEA1DM_V)
    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
  endif

  deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(mask_ibool2,stat=ier); if(ier /= 0) stop 'error in deallocate'

! only create mass matrix and save all the final arrays in the second pass
  else if(ipass == 2) then

! copy the theoretical number of points for the second pass
  nglob = nglob_theor

! count number of anisotropic elements in current region
! should be zero in all the regions except in the mantle
  nspec_tiso = count(idoubling(1:nspec) == IFLAG_220_80) + count(idoubling(1:nspec) == IFLAG_80_MOHO)

!!!! DM DM detection of the edges, coloring and permutation separately
  allocate(perm(nspec))
  if(USE_MESH_COLORING_INNER_OUTER) then

    stop 'USE_MESH_COLORING_INNER_OUTER should not be used in the serial case'

!! DK DK
!! DK DK the idea below of a different (regular) sorting to reduce cache misses did not work,
!! DK DK there is probably a small implementation problem somewhere in the mesher,
!! DK DK therefore I commented it out for now.
!! DK DK

!! DK DK for regular C version for CPUs: do not use colors but nonetheless put all the outer elements
!! DK DK first in order to be able to overlap non-blocking MPI communications with calculations
!   if(USE_REGULAR_C_CPU_VERSION) then

!   inumber_in_new_list_after_perm = 0

! first detect and list all the outer elements
!   do ispec = 1,nspec
!     if(is_on_a_slice_edge(ispec)) then
!       inumber_in_new_list_after_perm = inumber_in_new_list_after_perm + 1
!       perm(ispec) = inumber_in_new_list_after_perm
!     endif
!   enddo

! store the total number of outer elements found
!#ifdef USE_MPI
!   call MPI_ALLREDUCE(inumber_in_new_list_after_perm,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
!   call MPI_ALLREDUCE(inumber_in_new_list_after_perm,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
!#else
!   nspec_outer_min_global = inumber_in_new_list_after_perm
!   nspec_outer_max_global = inumber_in_new_list_after_perm
!#endif
!   if(myrank == 0) then
!     open(unit=99,file='../DATABASES_FOR_SOLVER/values_from_mesher_nspec_outer.h',status='unknown')
!     write(99,*) '#define NSPEC_OUTER ',nspec_outer_max_global
!     write(99,*) '// NSPEC_OUTER_min = ',nspec_outer_min_global
!     write(99,*) '// NSPEC_OUTER_max = ',nspec_outer_max_global
!     close(99)
!   endif

! just in case, test that we have detected outer elements
!   if(inumber_in_new_list_after_perm <= 0) stop 'fatal error: no outer elements detected!'

! then detect and list all the inner elements
!   do ispec = 1,nspec
!     if(.not. is_on_a_slice_edge(ispec)) then
!       inumber_in_new_list_after_perm = inumber_in_new_list_after_perm + 1
!       perm(ispec) = inumber_in_new_list_after_perm
!     endif
!   enddo

! test that all the elements have been used once and only once
!   if(inumber_in_new_list_after_perm /= nspec) stop 'fatal error: inumber_in_new_list_after_perm not equal to nspec'

!   else

! Il y a une routine get_perm_color que l'on ne peut pas enlever :
! c'est la routine de coloriage. Elle sert a grouper les elements en
! ensembles d'elements non-jointifs (pas de points globaux communs) pour
! garder accel() coherent lors de la somme des contributions en parallele sur la carte.
!
! De plus, a ce moment on separe aussi les elements externes des internes
! pour pouvoir faire des comms non boquantes avec recouvrement.
!
! This call is a bit expensive because it needs to build the adjacency table
! but it is necessary and cannot be removed.
    allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1))
    call get_perm_color(is_on_a_slice_edge,ibool,perm,nspec,nglob, &
      nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer,first_elem_number_in_this_color,myrank)

! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
    first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) = nspec + 1

    allocate(number_of_elements_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements))

! save mesh coloring
    open(unit=99,file=prname(1:len_trim(prname))//'number_of_elements_in_this_color.dat',status='unknown')

! number of colors for outer elements
    write(99,*) nb_colors_outer_elements

! number of colors for inner elements
    write(99,*) nb_colors_inner_elements

! number of elements in each color
    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
      number_of_elements_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) - first_elem_number_in_this_color(icolor)
      write(99,*) number_of_elements_in_this_color(icolor)
    enddo
    close(99)

! check that the sum of all the numbers of elements found in each color is equal
! to the total number of elements in the mesh
    if(sum(number_of_elements_in_this_color) /= nspec) then
      print *,'nspec = ',nspec
      print *,'total number of elements in all the colors of the mesh = ',sum(number_of_elements_in_this_color)
      stop 'incorrect total number of elements in all the colors of the mesh'
    endif

! check that the sum of all the numbers of elements found in each color for the outer elements is equal
! to the total number of outer elements found in the mesh
    if(sum(number_of_elements_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
      print *,'nspec_outer = ',nspec_outer
      print *,'total number of elements in all the colors of the mesh for outer elements = ',sum(number_of_elements_in_this_color)
      stop 'incorrect total number of elements in all the colors of the mesh for outer elements'
    endif

! store the total number of outer elements found
#ifdef USE_MPI
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
#else
    nspec_outer_min_global = nspec_outer
    nspec_outer_max_global = nspec_outer
#endif
    if(myrank == 0) then
      open(unit=99,file='../DATABASES_FOR_SOLVER/values_from_mesher_nspec_outer.h',status='unknown')
      write(99,*) '#define NSPEC_OUTER ',nspec_outer_max_global
      write(99,*) '// NSPEC_OUTER_min = ',nspec_outer_min_global
      write(99,*) '// NSPEC_OUTER_max = ',nspec_outer_max_global
      close(99)
    endif

    deallocate(first_elem_number_in_this_color)
    deallocate(number_of_elements_in_this_color)

!   endif  ! endif of the above section commented out because the idea did not work

  else

!   print *,'be careful, USE_MESH_COLORING_INNER_OUTER must always been set to .true. even for the regular C version for CPUs!'
!   print *,'generating identity permutation'
    do ispec = 1,nspec
      perm(ispec) = ispec
    enddo

  endif

  if (myrank == 0) then

! write a header file for the Fortran version of the solver
    open(unit=99,file='../DATABASES_FOR_SOLVER/values_from_mesher_f90.h',status='unknown')
    write(99,*) 'integer, parameter :: NSPEC = ',nspec
    write(99,*) 'integer, parameter :: NGLOB = ',nglob
    write(99,*) 'integer, parameter :: NSTEP = ',nstep
    write(99,*) 'real(kind=4), parameter :: deltat = ',DT
    write(99,*)
    write(99,*) '! element number of the source and of the station'
    write(99,*) '! after permutation of the elements by mesh coloring'
    write(99,*) '! and inner/outer set splitting in the mesher'
    write(99,*) '! (before permutation they are 1000 and NSPEC - 1000)'
    write(99,*) 'integer, parameter :: NSPEC_SOURCE = ',perm(1000)
    write(99,*) 'integer, parameter :: NSPEC_STATION = ',perm(NSPEC - 1000)
! save coordinates of the seismic source
!   write(99,*) xstore(2,2,2,10);
!   write(99,*) ystore(2,2,2,10);
!   write(99,*) zstore(2,2,2,10);

! save coordinates of the seismic station
!   write(99,*) xstore(2,2,2,nspec-10);
!   write(99,*) ystore(2,2,2,nspec-10);
!   write(99,*) zstore(2,2,2,nspec-10);
    close(99)

! write a header file for the C version of the solver
    open(unit=99,file='../DATABASES_FOR_SOLVER/values_from_mesher_C.h',status='unknown')
    write(99,*) '#define NSPEC ',nspec
    write(99,*) '#define NGLOB ',nglob
    write(99,*) '#define NSTEP ',nstep
!!!!!!! DK DK use 1000 time steps only for the scaling tests
!!!!!!!!!    write(99,*) '#define NSTEP 1000 // 5000 // ',nstep
! put an "f" at the end to force single precision
    write(99,"('#define deltat ',e18.10,'f')") DT
    write(99,*) '#define NGLOB2DMAX_XMIN_XMAX ',npoin2D_xi
    write(99,*) '#define NGLOB2DMAX_YMIN_YMAX ',npoin2D_eta
    write(99,*) '#define NGLOB2DMAX_ALL ',max(npoin2D_xi,npoin2D_eta)
    write(99,*) '#define NPROC_XI ',NPROC_XI
    write(99,*) '#define NPROC_ETA ',NPROC_ETA
    write(99,*)
    write(99,*) '// element and MPI slice number of the source and the station'
    write(99,*) '// after permutation of the elements by mesh coloring'
    write(99,*) '// and inner/outer set splitting in the mesher'
    write(99,*) '// (before permutation they are 1000 and NSPEC - 1000)'
    write(99,*) '#define RANK_SOURCE 0'
    write(99,*) '#define NSPEC_SOURCE ',perm(1000)
    write(99,*)
    write(99,*) '#define RANK_STATION (NPROC_XI*NPROC_ETA - 1)'
    write(99,*) '#define NSPEC_STATION ',perm(NSPEC - 1000)
    close(99)

  endif

  if (USE_MESH_COLORING_INNER_OUTER) then
    allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))
    call permute_elements_real(xixstore,temp_array_real,perm,nspec)
    call permute_elements_real(xiystore,temp_array_real,perm,nspec)
    call permute_elements_real(xizstore,temp_array_real,perm,nspec)
    call permute_elements_real(etaxstore,temp_array_real,perm,nspec)
    call permute_elements_real(etaystore,temp_array_real,perm,nspec)
    call permute_elements_real(etazstore,temp_array_real,perm,nspec)
    call permute_elements_real(gammaxstore,temp_array_real,perm,nspec)
    call permute_elements_real(gammaystore,temp_array_real,perm,nspec)
    call permute_elements_real(gammazstore,temp_array_real,perm,nspec)
    call permute_elements_real(kappavstore,temp_array_real,perm,nspec)
    call permute_elements_real(muvstore,temp_array_real,perm,nspec)
    call permute_elements_real(rhostore,temp_array_real,perm,nspec)
    deallocate(temp_array_real)

    allocate(temp_array_dble(NGLLX,NGLLY,NGLLZ,nspec))
    call permute_elements_dble(xstore,temp_array_dble,perm,nspec)
    call permute_elements_dble(ystore,temp_array_dble,perm,nspec)
    call permute_elements_dble(zstore,temp_array_dble,perm,nspec)
    deallocate(temp_array_dble)

    ! permutation of ibool
    allocate(temp_array_int(NGLLX,NGLLY,NGLLZ,nspec))
    call permute_elements_integer(ibool,temp_array_int,perm,nspec)
    deallocate(temp_array_int)

    ! permutation of idoubling
    allocate(temp_array_1D_int(nspec))
    temp_array_1D_int(:) = idoubling(:)
    do ispec = 1,nspec
      idoubling(perm(ispec)) = temp_array_1D_int(ispec)
    enddo
    deallocate(temp_array_1D_int)

    deallocate(perm)
  endif
!!!! DM DM coloring and permutation

! call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
!     dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
!     ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
!     nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
!             jacobian2D_xmin,jacobian2D_xmax, &
!             jacobian2D_ymin,jacobian2D_ymax, &
!             jacobian2D_bottom,jacobian2D_top, &
!             normal_xmin,normal_xmax, &
!             normal_ymin,normal_ymax, &
!             normal_bottom,normal_top, &
!             NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!             NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

! creating mass matrix in this slice (will be fully assembled in the solver)
  allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec

! suppress fictitious elements in central cube
  if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        weight = wxgll(i)*wygll(j)*wzgll(k)
        iglobnum = ibool(i,j,k,ispec)

! compute the jacobian
        xixl = xixstore(i,j,k,ispec)
        xiyl = xiystore(i,j,k,ispec)
        xizl = xizstore(i,j,k,ispec)
        etaxl = etaxstore(i,j,k,ispec)
        etayl = etaystore(i,j,k,ispec)
        etazl = etazstore(i,j,k,ispec)
        gammaxl = gammaxstore(i,j,k,ispec)
        gammayl = gammaystore(i,j,k,ispec)
        gammazl = gammazstore(i,j,k,ispec)

        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

! definition depends if region is fluid or solid
  if(iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight)
    else
      rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
    endif

! fluid in outer core
  else if(iregion_code == IREGION_OUTER_CORE) then

! no anisotropy in the fluid, use kappav

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)))
    else
      rmass(iglobnum) = rmass(iglobnum) + &
             jacobianl * weight * rhostore(i,j,k,ispec) / kappavstore(i,j,k,ispec)
    endif

  else
    call exit_MPI(myrank,'wrong region code')
  endif

      enddo
    enddo
  enddo
  enddo

! save the binary files
! save ocean load mass matrix as well if oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then

! adding ocean load mass matrix at the top of the crust for oceans
  nglob_oceans = nglob
  allocate(rmass_ocean_load(nglob_oceans))

! create ocean load mass matrix for degrees of freedom at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

! add contribution of the oceans
! for surface elements exactly at the top of the crust (ocean bottom)
    do ispec2D_top_crust = 1,NSPEC2D_TOP

      ispec_oceans = ibelm_top(ispec2D_top_crust)

      iz_oceans = NGLLZ

      do ix_oceans = 1,NGLLX
        do iy_oceans = 1,NGLLY

        iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! compute local height of oceans
        if(ISOTROPIC_3D_MANTLE) then

! get coordinates of current point
          xval = xstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          yval = ystore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          zval = zstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! map to latitude and longitude for bathymetry routine
          call xyz_2_rthetaphi_dble(xval,yval,zval,rval,thetaval,phival)
          call reduce(thetaval,phival)

! convert the geocentric colatitude to a geographic colatitude
          colat = PI/2.0d0 - datan(1.006760466d0*dcos(thetaval)/dmax1(TINYVAL,dsin(thetaval)))

! get geographic latitude and longitude in degrees
          lat = 90.0d0 - colat*180.0d0/PI
          lon = phival*180.0d0/PI
          elevation = 0.d0

! compute elevation at current point
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)

! non-dimensionalize the elevation, which is in meters
! and suppress positive elevation, which means no oceans
          if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation) / R_EARTH
          endif

        else
          height_oceans = THICKNESS_OCEANS_PREM
        endif

! take into account inertia of water column
        weight = wxgll(ix_oceans)*wygll(iy_oceans)*dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_top_crust)) &
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

  call save_arrays_solver(prname,xixstore,xiystore,xizstore, &
                          etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore, &
                          kappavstore,muvstore,ibool,rmass,nspec,nglob,myrank,NPROCTOT,xstore,ystore,zstore)

  deallocate(rmass,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(rmass_ocean_load,stat=ier); if(ier /= 0) stop 'error in deallocate'

! compute volume, bottom and top area of that part of the slice
  volume_local = ZERO

  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

! compute the jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          volume_local = volume_local + dble(jacobianl)*weight

        enddo
      enddo
    enddo
  enddo

  else
    stop 'there cannot be more than two passes in mesh creation'

  endif  ! end of test if first or second pass

! deallocate these arrays after each pass because they have a different size in each pass to save memory
  deallocate(xixstore,xiystore,xizstore,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(etaxstore,etaystore,etazstore,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(gammaxstore,gammaystore,gammazstore,stat=ier); if(ier /= 0) stop 'error in deallocate'

! deallocate arrays
  deallocate(rhostore,kappavstore,kappahstore)
  deallocate(muvstore,muhstore)
  deallocate(eta_anisostore)

  deallocate(c11store)
  deallocate(c12store)
  deallocate(c13store)
  deallocate(c14store)
  deallocate(c15store)
  deallocate(c16store)
  deallocate(c22store)
  deallocate(c23store)
  deallocate(c24store)
  deallocate(c25store)
  deallocate(c26store)
  deallocate(c33store)
  deallocate(c34store)
  deallocate(c35store)
  deallocate(c36store)
  deallocate(c44store)
  deallocate(c45store)
  deallocate(c46store)
  deallocate(c55store)
  deallocate(c56store)
  deallocate(c66store)

  deallocate(iboun)
  deallocate(xigll,yigll,zigll)
  deallocate(wxgll,wygll,wzgll)
  deallocate(shape3D,dershape3D)
  deallocate(shape2D_x,shape2D_y,shape2D_bottom,shape2D_top)
  deallocate(dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top)
  deallocate(ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax)
  deallocate(ibelm_bottom,ibelm_top)
  deallocate(jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax)
  deallocate(jacobian2D_bottom,jacobian2D_top)
  deallocate(normal_xmin,normal_xmax,normal_ymin,normal_ymax)
  deallocate(normal_bottom,normal_top)
  deallocate(iMPIcut_xi,iMPIcut_eta)

  deallocate(nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta)
  deallocate(rho_vp,rho_vs)

  deallocate(Qmu_store)
  deallocate(tau_e_store)

  end subroutine create_regions_mesh

