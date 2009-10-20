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
!
! United States and French Government Sponsorship Acknowledged.

module specfem_par

  implicit none

  include "constants.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

! standard include of the MPI library
!  include 'mpif.h'



! memory variables and standard linear solids for attenuation
  double precision, dimension(N_SLS) :: tau_mu_dble,tau_sigma_dble,beta_dble
  double precision factor_scale_dble,one_minus_sum_beta_dble
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: tau_mu,tau_sigma,beta
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: factor_scale,one_minus_sum_beta

  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: tauinv,factor_common, alphaval,betaval,gammaval
  integer iattenuation
  double precision scale_factor

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  integer :: NSPEC_ATTENUATION_AB
  integer, dimension(:,:,:,:),allocatable :: iflag_attenuation_store

! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: b_alphaval, b_betaval, b_gammaval
!! DK DK array not created yet for CUBIT
! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL,N_SLS) :: &
!            b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz
! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ATT_AND_KERNEL) ::  b_epsilondev_xx, &
!            b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz
! ADJOINT

! use integer array to store topography values
  integer NX_TOPO,NY_TOPO
  double precision ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  character(len=100) topo_file
  integer, dimension(:,:), allocatable :: itopo_bathy

! absorbing boundaries
!  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax
!  integer, dimension(:), allocatable :: ibelm_ymin,ibelm_ymax
!  integer, dimension(:), allocatable :: ibelm_bottom
!  integer, dimension(:), allocatable :: ibelm_top
!!  integer :: NSPEC2DMAX_XMIN_XMAX_ext,NSPEC2DMAX_YMIN_YMAX_ext
!  ! local indices i,j,k of all GLL points on xmin boundary in the element
!  integer,dimension(:,:,:,:),allocatable :: ibelm_gll_xmin,ibelm_gll_xmax, &
!                                          ibelm_gll_ymin,ibelm_gll_ymax, &
!                                          ibelm_gll_bottom,ibelm_gll_top  
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_xmin,jacobian2D_xmax
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_ymin,jacobian2D_ymax
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_bottom
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: jacobian2D_top
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_xmin,normal_xmax
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable  :: normal_ymin,normal_ymax
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable  :: normal_bottom
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable  :: normal_top

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorbing_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: absorbing_boundary_jacobian2D
  integer, dimension(:,:,:), allocatable :: absorbing_boundary_ijk
  integer, dimension(:), allocatable :: absorbing_boundary_ispec
  integer :: num_absorbing_boundary_faces

! free surface  
  integer :: nspec2D_top,ispec2D
  integer, dimension(:), allocatable :: ibelm_top
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable  :: normal_top
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: jacobian2D_top
  real(kind=CUSTOM_REAL) :: nx,ny,nz

!! DK DK array not created yet for CUBIT
! integer, dimension(NSPEC2D_TOP_VAL) :: ibelm_top
! real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_VAL) :: normal_top

!! DK DK array not created yet for CUBIT
! Moho mesh
! integer,dimension(NSPEC2D_MOHO_BOUN) :: ibelm_moho_top, ibelm_moho_bot
! real(CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO_BOUN) :: normal_moho
! integer :: nspec2D_moho

!! DK DK array not created yet for CUBIT
! buffers for send and receive between faces of the slices and the chunks
! real(kind=CUSTOM_REAL), dimension(NDIM,NPOIN2DMAX_XY_VAL) :: buffer_send_faces_vector,buffer_received_faces_vector

! -----------------

! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
        kappastore,mustore

! flag for sediments
  logical, dimension(:), allocatable :: not_fully_in_bedrock
  logical, dimension(:,:,:,:), allocatable :: flag_sediments

! Stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! local to global mapping
  integer, dimension(:), allocatable :: idoubling

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! additional mass matrix for ocean load
! ocean load mass matrix is always allocated statically even if no oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load
  logical, dimension(:), allocatable :: updated_dof_ocean_load
  real(kind=CUSTOM_REAL) additional_term,force_normal_comp

! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ,veloc,accel

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  real(kind=CUSTOM_REAL) hp1,hp2,hp3

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

! time scheme
  real(kind=CUSTOM_REAL) deltat,deltatover2,deltatsqover2

! ADJOINT
  real(kind=CUSTOM_REAL) b_additional_term,b_force_normal_comp
!! DK DK array not created yet for CUBIT
! real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_ADJOINT) :: b_displ, b_veloc, b_accel
! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ADJOINT) :: rho_kl, mu_kl, kappa_kl, &
!   rhop_kl, beta_kl, alpha_kl
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: absorb_xmin, absorb_xmax, &  
!       absorb_ymin, absorb_ymax, absorb_zmin ! for absorbing b.c.
!  integer reclen_xmin, reclen_xmax, reclen_ymin, reclen_ymax, reclen_zmin

  real(kind=CUSTOM_REAL) b_deltat, b_deltatover2, b_deltatsqover2
! ADJOINT

  integer l

! Moho kernel
! integer ispec2D_moho_top, ispec2D_moho_bot, k_top, k_bot, ispec_top, ispec_bot, iglob_top, iglob_bot
!! DK DK array not created yet for CUBIT
! real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ,NSPEC2D_MOHO_BOUN) :: dsdx_top, dsdx_bot, b_dsdx_top, b_dsdx_bot
! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_MOHO_BOUN) :: moho_kl
! real(kind=CUSTOM_REAL) :: kernel_moho_top, kernel_moho_bot

! --------

! parameters for the source
  integer it,isource
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  integer yr,jda,ho,mi
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays
  double precision, dimension(:,:,:), allocatable :: nu_source
!ADJOINT
  character(len=150) adj_source_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: adj_sourcearrays
!ADJOINT
  double precision sec,stf
  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: t_cmt,hdur,hdur_gaussian
  double precision, dimension(:), allocatable :: utm_x_source,utm_y_source
  double precision, external :: comp_source_time_function
  double precision :: t0

! receiver information
  character(len=150) rec_filename,filtered_rec_filename,dummystring
  integer nrec,nrec_local,nrec_tot_found,irec_local,ios
  integer, allocatable, dimension(:) :: islice_selected_rec,ispec_selected_rec,number_receiver_global
  double precision, allocatable, dimension(:) :: xi_receiver,eta_receiver,gamma_receiver
  double precision hlagrange
! ADJOINT
  integer nrec_simulation, nadj_rec_local
! source frechet derivatives
  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ), eps_s(NDIM,NDIM), eps_m_s(NDIM), stf_deltat
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store
! ADJOINT

! timing information for the stations
  double precision, allocatable, dimension(:,:,:) :: nu
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! seismograms
  double precision dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms_d,seismograms_v,seismograms_a
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: seismograms_eps

  integer i,j,k,ispec,irec,iglob

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hetar,hpxir,hpetar,hgammar,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hetar_store,hgammar_store

! 2-D addressing and buffers for summation between slices
! integer, dimension(NPOIN2DMAX_XMIN_XMAX_VAL) :: iboolleft_xi,iboolright_xi
! integer, dimension(NPOIN2DMAX_YMIN_YMAX_VAL) :: iboolleft_eta,iboolright_eta

! for addressing of the slices
! integer, dimension(0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL) :: addressing

! proc numbers for MPI
  integer :: myrank

! integer npoin2D_xi,npoin2D_eta

! integer iproc_xi,iproc_eta

! maximum of the norm of the displacement
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all
  integer:: Usolidnorm_index(1)
! ADJOINT
! real(kind=CUSTOM_REAL) b_Usolidnorm, b_Usolidnorm_all
! ADJOINT

! timer MPI
  double precision, external :: wtime
  integer :: ihours,iminutes,iseconds,int_tCPU, &
             ihours_remain,iminutes_remain,iseconds_remain,int_t_remain, &
             ihours_total,iminutes_total,iseconds_total,int_t_total
  double precision :: time_start,tCPU,t_remain,t_total


! parameters read from parameter file
  integer NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,HDUR_MOVIE

  logical TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=150) OUTPUT_FILES,LOCAL_PATH,prname,prname_Q

! parameters deduced from parameters read from file
  integer NPROC

  !integer :: NSPEC2D_BOTTOM
  !integer :: NSPEC2D_TOP
  
  integer :: NSPEC_AB, NGLOB_AB

! names of the data files for all the processors in MPI
  character(len=150) outputname

! Stacey conditions put back
  !integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,ispec2D
  !real(kind=CUSTOM_REAL) nx,ny,nz
  !integer, dimension(:,:),allocatable :: nimin,nimax,nkmin_eta
  !integer, dimension(:,:),allocatable :: njmin,njmax,nkmin_xi

! to save movie frames
  integer ipoin, nmovie_points, iloc, iorderi(NGNOD2D), iorderj(NGNOD2D)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
      store_val_x,store_val_y,store_val_z, &
      store_val_ux,store_val_uy,store_val_uz, &
      store_val_norm_displ,store_val_norm_veloc,store_val_norm_accel
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
      store_val_x_all,store_val_y_all,store_val_z_all, &
      store_val_ux_all,store_val_uy_all,store_val_uz_all

! to save full 3D snapshot of velocity
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dvxdxl,dvxdyl,dvxdzl,dvydxl,dvydyl,dvydzl,dvzdxl,dvzdyl,dvzdzl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable::  div, curl_x, curl_y, curl_z

! for assembling in case of external mesh
  integer :: ninterfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer, dimension(:), allocatable :: my_neighbours_ext_mesh
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar_ext_mesh
  integer, dimension(:), allocatable :: request_send_scalar_ext_mesh
  integer, dimension(:), allocatable :: request_recv_scalar_ext_mesh
  integer, dimension(:), allocatable :: request_send_vector_ext_mesh
  integer, dimension(:), allocatable :: request_recv_vector_ext_mesh

! for detecting surface receivers and source in case of external mesh
  integer, dimension(:), allocatable :: valence_external_mesh
  logical, dimension(:), allocatable :: iglob_is_surface_external_mesh
  logical, dimension(:), allocatable :: ispec_is_surface_external_mesh
  integer, dimension(:,:), allocatable :: buffer_send_scalar_i_ext_mesh
  integer, dimension(:,:), allocatable :: buffer_recv_scalar_i_ext_mesh
  integer :: nfaces_surface_external_mesh
  integer :: nfaces_surface_glob_ext_mesh
  integer,dimension(:),allocatable :: nfaces_perproc_surface_ext_mesh
  integer,dimension(:),allocatable :: faces_surface_offset_ext_mesh
  integer,dimension(:,:),allocatable :: faces_surface_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_y_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_z_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x_all_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_y_all_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_z_all_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uy_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uz_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux_all_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uy_all_external_mesh
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uz_all_external_mesh
  integer :: ii,jj,kk

! for communications overlapping
  logical, dimension(:), allocatable :: ispec_is_inner_ext_mesh
  logical, dimension(:), allocatable :: iglob_is_inner_ext_mesh
  integer :: iinterface

!  integer, dimension(:),allocatable :: spec_inner, spec_outer
!  integer :: nspec_inner,nspec_outer
  

!!!! NL NL REGOLITH : regolith layer for asteroid
!!$  double precision, external :: materials_ext_mesh
!!$  logical, dimension(:), allocatable :: ispec_is_regolith
!!$  real(kind=CUSTOM_REAL) :: weight, jacobianl
!!!! NL NL REGOLITH

  
end module
