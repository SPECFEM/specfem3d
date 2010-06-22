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

module constants

  include "constants.h"

end module constants

!=====================================================================

module specfem_par

! main parameter module for specfem simulations

  use constants
  
  implicit none

! attenuation  
  integer :: NSPEC_ATTENUATION_AB
  integer, dimension(:,:,:,:),allocatable :: iflag_attenuation_store

! use integer array to store topography values
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  character(len=100) :: topo_file
  integer, dimension(:,:), allocatable :: itopo_bathy

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ijk
  integer, dimension(:), allocatable :: abs_boundary_ispec
  integer :: num_abs_boundary_faces  

! free surface arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ijk
  integer, dimension(:), allocatable :: free_surface_ispec
  integer :: num_free_surface_faces

! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! material properties
  ! isotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore

! additional mass matrix for ocean load
! ocean load mass matrix is always allocated statically even if no oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! time scheme
  real(kind=CUSTOM_REAL) deltat,deltatover2,deltatsqover2

! time loop step
  integer :: it 

! parameters for the source
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays
  double precision, dimension(:,:,:), allocatable :: nu_source
  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: t_cmt,hdur,hdur_gaussian
  double precision, dimension(:), allocatable :: utm_x_source,utm_y_source
  double precision, external :: comp_source_time_function
  double precision :: t0
  real(kind=CUSTOM_REAL) :: stf_used_total
  integer :: NSOURCES
  
! receiver information
  character(len=256) :: rec_filename,filtered_rec_filename,dummystring
  integer :: nrec,nrec_local,nrec_tot_found
  integer :: nrec_simulation
  integer, allocatable, dimension(:) :: islice_selected_rec,ispec_selected_rec,number_receiver_global
  double precision, allocatable, dimension(:) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store

! timing information for the stations
  double precision, allocatable, dimension(:,:,:) :: nu
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! seismograms
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms_d,seismograms_v,seismograms_a

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

! proc numbers for MPI
  integer :: myrank

! timer MPI
  double precision, external :: wtime
  double precision :: time_start

! parameters read from parameter file
  integer :: NPROC_XI,NPROC_ETA
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE
  integer :: SIMULATION_TYPE

  double precision :: DT
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
            OCEANS,ABSORBING_CONDITIONS,ANISOTROPY
            
  logical :: SAVE_FORWARD,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical :: SUPPRESS_UTM_PROJECTION
  
  integer :: NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=256) OUTPUT_FILES,LOCAL_PATH,prname,prname_Q

! parameters deduced from parameters read from file
  integer :: NPROC
  integer :: NSPEC_AB, NGLOB_AB

! names of the data files for all the processors in MPI
  character(len=256) outputname

! for assembling in case of external mesh
  integer :: num_interfaces_ext_mesh
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
  logical, dimension(:), allocatable :: iglob_is_surface_external_mesh
  logical, dimension(:), allocatable :: ispec_is_surface_external_mesh

! MPI partition surfaces 
  logical, dimension(:), allocatable :: ispec_is_inner
  logical, dimension(:), allocatable :: iglob_is_inner

! maximum of the norm of the displacement
  real(kind=CUSTOM_REAL) Usolidnorm,Usolidnorm_all
  integer:: Usolidnorm_index(1)

! maximum speed in velocity model
  real(kind=CUSTOM_REAL):: model_speed_max

!!!! NL NL REGOLITH : regolith layer for asteroid
!!$  double precision, external :: materials_ext_mesh
!!$  logical, dimension(:), allocatable :: ispec_is_regolith
!!$  real(kind=CUSTOM_REAL) :: weight, jacobianl
!!!! NL NL REGOLITH


! ADJOINT parameters

  ! time scheme
  real(kind=CUSTOM_REAL) b_deltat, b_deltatover2, b_deltatsqover2

  ! absorbing stacey wavefield parts
  integer :: b_num_abs_boundary_faces

  ! Moho mesh
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC_BOUN,NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot

  ! adjoint sources
  character(len=256) adj_source_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: adj_sourcearrays
  integer :: nadj_rec_local
  ! adjoint source frechet derivatives
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Mxx_der,Myy_der,&
    Mzz_der,Mxy_der,Mxz_der,Myz_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: seismograms_eps  

  ! adjoint elements
  integer :: NSPEC_ADJOINT, NGLOB_ADJOINT

  ! norm of the backward displacement
   real(kind=CUSTOM_REAL) b_Usolidnorm, b_Usolidnorm_all

  
end module specfem_par


!=====================================================================

module specfem_par_elastic

! parameter module for elastic solver

  use constants,only: CUSTOM_REAL,N_SLS,NUM_REGIONS_ATTENUATION
  implicit none

! memory variables and standard linear solids for attenuation
  double precision, dimension(N_SLS) :: tau_mu_dble,tau_sigma_dble,beta_dble
  double precision factor_scale_dble,one_minus_sum_beta_dble
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: tau_mu,tau_sigma,beta
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: factor_scale,one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: &
    tauinv,factor_common, alphaval,betaval,gammaval
    
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ,veloc,accel
  
! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! Stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

  ! anisotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
            c11store,c12store,c13store,c14store,c15store,c16store,&
            c22store,c23store,c24store,c25store,c26store,c33store,&
            c34store,c35store,c36store,c44store,c45store,c46store,&
            c55store,c56store,c66store
  integer :: NSPEC_ANISO

! material flag
  logical, dimension(:), allocatable :: ispec_is_elastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_elastic
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic

  logical :: ELASTIC_SIMULATION


! ADJOINT elastic 

  ! (backward/reconstructed) wavefields
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ, b_veloc, b_accel

  ! backward attenuation arrays
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: &
    b_alphaval, b_betaval, b_gammaval
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz      
  integer:: NSPEC_ATT_AND_KERNEL

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl, &
    rhop_kl, beta_kl, alpha_kl

  ! topographic (Moho) kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:),allocatable :: &
    dsdx_top, dsdx_bot, b_dsdx_top, b_dsdx_bot
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: moho_kl
  integer :: ispec2D_moho_top,ispec2D_moho_bot

  ! absorbing stacey wavefield parts
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_absorb_field
  integer :: b_reclen_field

  ! for assembling backward field
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_mesh
  integer, dimension(:), allocatable :: b_request_send_vector_ext_mesh
  integer, dimension(:), allocatable :: b_request_recv_vector_ext_mesh
      
end module specfem_par_elastic

!=====================================================================

module specfem_par_acoustic

! parameter module for elastic solver

  use constants,only: CUSTOM_REAL
  implicit none

! potential
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic, &
                        potential_dot_acoustic,potential_dot_dot_acoustic

! density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore  

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic

! acoustic-elastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ijk
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer :: num_coupling_ac_el_faces

! material flag
  logical, dimension(:), allocatable :: ispec_is_acoustic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_acoustic
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  
  logical :: ACOUSTIC_SIMULATION

! ADJOINT acoustic

  ! (backward/reconstructed) wavefield potentials
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_acoustic, &
                        b_potential_dot_acoustic,b_potential_dot_dot_acoustic
  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl, &
    rhop_ac_kl, alpha_ac_kl

  ! absorbing stacey wavefield parts
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_absorb_potential
  integer :: b_reclen_potential

  ! for assembling backward field
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_recv_scalar_ext_mesh
  integer, dimension(:), allocatable :: b_request_send_scalar_ext_mesh
  integer, dimension(:), allocatable :: b_request_recv_scalar_ext_mesh

end module specfem_par_acoustic

!=====================================================================

module specfem_par_poroelastic

! parameter module for elastic solver

  use constants,only: CUSTOM_REAL
  implicit none

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_solid_poroelastic,&
    rmass_fluid_poroelastic

! material flag
  logical, dimension(:), allocatable :: ispec_is_poroelastic

  logical :: POROELASTIC_SIMULATION
  
end module specfem_par_poroelastic


!=====================================================================

module specfem_par_movie

! parameter module for movies/shakemovies

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD2D

  implicit none

! to save full 3D snapshot of velocity (movie volume
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: div, curl_x, curl_y, curl_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: velocity_x,velocity_y,velocity_z
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dvxdxl,dvxdyl,&
                                dvxdzl,dvydxl,dvydyl,dvydzl,dvzdxl,dvzdyl,dvzdzl

! shakemovies  
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

! movie volume
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  real(kind=CUSTOM_REAL) hp1,hp2,hp3

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

! for storing surface of external mesh
  integer,dimension(:),allocatable :: nfaces_perproc_surface_ext_mesh
  integer,dimension(:),allocatable :: faces_surface_offset_ext_mesh
  integer,dimension(:,:),allocatable :: faces_surface_ext_mesh
  integer,dimension(:),allocatable :: faces_surface_ext_mesh_ispec
  integer :: nfaces_surface_ext_mesh
  integer :: nfaces_surface_glob_ext_mesh
  ! face corner indices
  integer :: iorderi(NGNOD2D),iorderj(NGNOD2D)

! movie parameters
  double precision :: HDUR_MOVIE
  integer :: NTSTEP_BETWEEN_FRAMES  
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
            USE_HIGHRES_FOR_MOVIES

  logical :: MOVIE_SIMULATION  

end module specfem_par_movie

