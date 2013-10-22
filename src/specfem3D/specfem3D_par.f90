!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

! parameters deduced from parameters read from file
  integer :: NPROC
  integer :: NSPEC_AB, NGLOB_AB

! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! material properties
  ! isotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore

! density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore

! CUDA mesh pointer<->integer wrapper
  integer(kind=8) :: Mesh_pointer

! Global GPU toggle. Set in Par_file
  logical :: GPU_MODE

! use integer array to store topography values
  integer :: NX_TOPO,NY_TOPO
  integer, dimension(:,:), allocatable :: itopo_bathy

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ijk
  integer, dimension(:), allocatable :: abs_boundary_ispec
  integer :: num_abs_boundary_faces
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

! free surface arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ijk
  integer, dimension(:), allocatable :: free_surface_ispec
  integer :: num_free_surface_faces

! VM for new method
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: Veloc_dsm_boundary,Tract_dsm_boundary

! attenuation
  integer :: NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa
  character(len=256) prname_Q

! additional mass matrix for ocean load
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! time scheme
  real(kind=CUSTOM_REAL) deltat,deltatover2,deltatsqover2

! time loop step
  integer :: it

! VM for new  method
  integer :: it_dsm

! parameters for the source
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays
  double precision, dimension(:,:,:), allocatable :: nu_source
  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: tshift_src,hdur,hdur_gaussian,hdur_tiny
  double precision, dimension(:), allocatable :: utm_x_source,utm_y_source
  double precision, external :: comp_source_time_function
  double precision :: t0
  real(kind=CUSTOM_REAL) :: stf_used_total
  integer :: NSOURCES,nsources_local
  ! source encoding
  ! for acoustic sources: takes +/- 1 sign, depending on sign(Mxx)[ = sign(Myy) = sign(Mzz)
  ! since they have to equal in the acoustic setting]
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pm1_source_encoding

! receiver information
  character(len=256) :: rec_filename,filtered_rec_filename,dummystring
  integer :: nrec,nrec_local,nrec_tot_found
  integer :: nrec_simulation
  integer, dimension(:), allocatable :: islice_selected_rec,ispec_selected_rec
  integer, dimension(:), allocatable :: number_receiver_global
  double precision, dimension(:), allocatable :: xi_receiver,eta_receiver,gamma_receiver
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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hetar,hpxir,hpetar,hgammar,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hetar_store,hgammar_store

! proc numbers for MPI
  integer :: myrank, sizeprocs

! timer MPI
  double precision, external :: wtime
  double precision :: time_start

! parameters for a force source located exactly at a grid point
  logical :: USE_FORCE_POINT_SOURCE
  double precision, dimension(:), allocatable :: factor_force_source
  double precision, dimension(:), allocatable :: comp_dir_vect_source_E
  double precision, dimension(:), allocatable :: comp_dir_vect_source_N
  double precision, dimension(:), allocatable :: comp_dir_vect_source_Z_UP

! parameters
  integer :: SIMULATION_TYPE
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE
  integer :: IMODEL,NGNOD,NGNOD2D

  double precision :: DT,OLSEN_ATTENUATION_RATIO,f0_FOR_PML

  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
            APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,STACEY_ABSORBING_CONDITIONS,ANISOTROPY, &
            STACEY_INSTEAD_OF_FREE_SURFACE

  logical :: FULL_ATTENUATION_SOLID,PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  logical :: GRAVITY

  logical :: SAVE_FORWARD,SAVE_MESH_FILES

  logical :: USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION

  logical :: SUPPRESS_UTM_PROJECTION

  integer :: NTSTEP_BETWEEN_OUTPUT_INFO

! parameters read from mesh parameter file
  integer :: NPROC_XI,NPROC_ETA
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  character(len=256) OUTPUT_FILES,LOCAL_PATH,TOMOGRAPHY_PATH,prname,dsmname,TRAC_PATH

  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
             ADIOS_FOR_KERNELS

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
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh_s
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_ext_mesh_s
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh_w
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_ext_mesh_w
  integer, dimension(:), allocatable :: request_send_vector_ext_mesh_s
  integer, dimension(:), allocatable :: request_recv_vector_ext_mesh_s
  integer, dimension(:), allocatable :: request_send_vector_ext_mesh_w
  integer, dimension(:), allocatable :: request_recv_vector_ext_mesh_w

! for detecting surface receivers and source in case of external mesh
  logical, dimension(:), allocatable :: iglob_is_surface_external_mesh
  logical, dimension(:), allocatable :: ispec_is_surface_external_mesh

! MPI partition surfaces
  logical, dimension(:), allocatable :: ispec_is_inner

! maximum speed in velocity model
  real(kind=CUSTOM_REAL):: model_speed_max

  ! gravity
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: minus_deriv_gravity,minus_g
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

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

  ! length of reading blocks
  integer :: NTSTEP_BETWEEN_READ_ADJSRC

  ! parameter module for noise simulations
  integer :: irec_master_noise, NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sigma_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: noise_surface_movie
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
             normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

end module specfem_par

!=====================================================================

module specfem_par_elastic

! parameter module for elastic solver

  use constants,only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ

  implicit none

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: one_minus_sum_beta,one_minus_sum_beta_kappa
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common,factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_sigma
  real(kind=CUSTOM_REAL) :: min_resolved_period
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: &
    alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_trace,R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: epsilon_trace_over_3

! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_adj_coupling

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! Stacey
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassy,rmassz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

  ! anisotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
            c11store,c12store,c13store,c14store,c15store,c16store,&
            c22store,c23store,c24store,c25store,c26store,c33store,&
            c34store,c35store,c36store,c44store,c45store,c46store,&
            c55store,c56store,c66store
  integer :: NSPEC_ANISO

! for attenuation and/or kernel simulations
  integer :: NSPEC_STRAIN_ONLY
  logical :: COMPUTE_AND_STORE_STRAIN

! material flag
  logical, dimension(:), allocatable :: ispec_is_elastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_elastic
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic

! mesh coloring
  integer :: num_colors_outer_elastic,num_colors_inner_elastic
  integer, dimension(:), allocatable :: num_elem_colors_elastic
  integer :: nspec_elastic

  logical :: ELASTIC_SIMULATION

! ADJOINT elastic

  ! (backward/reconstructed) wavefields
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ, b_veloc, b_accel

  ! backward attenuation arrays
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: &
    b_alphaval, b_betaval, b_gammaval
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_epsilon_trace_over_3

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl

  ! anisotropic kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: cijkl_kl

  ! approximate hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_kl

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

! parameter module for acoustic solver

  use constants,only: CUSTOM_REAL
  implicit none

! potential
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic,potential_dot_acoustic, &
                                    potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_adj_coupling

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassz_acoustic

! acoustic-elastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ijk
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer :: num_coupling_ac_el_faces

! acoustic-poroelastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_po_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_po_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_ac_po_ijk
  integer, dimension(:), allocatable :: coupling_ac_po_ispec
  integer :: num_coupling_ac_po_faces

! material flag
  logical, dimension(:), allocatable :: ispec_is_acoustic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_acoustic
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic

! mesh coloring
  integer :: num_colors_outer_acoustic,num_colors_inner_acoustic
  integer, dimension(:), allocatable :: num_elem_colors_acoustic
  integer :: nspec_acoustic

  logical :: ACOUSTIC_SIMULATION

! ADJOINT acoustic

  ! (backward/reconstructed) wavefield potentials
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_acoustic, &
                        b_potential_dot_acoustic,b_potential_dot_dot_acoustic
  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl, &
    rhop_ac_kl, alpha_ac_kl

  ! approximate hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_ac_kl

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

! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accels_poroelastic,velocs_poroelastic,displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accelw_poroelastic,velocw_poroelastic,displw_poroelastic

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    epsilonsdev_xx,epsilonsdev_yy,epsilonsdev_xy,epsilonsdev_xz,epsilonsdev_yz, &
    epsilonwdev_xx,epsilonwdev_yy,epsilonwdev_xy,epsilonwdev_xz,epsilonwdev_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    epsilons_trace_over_3,epsilonw_trace_over_3

! material properties
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: mustore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: etastore,tortstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: phistore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: kappaarraystore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: permstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vpI,rho_vpII,rho_vsI

! elastic-poroelastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_el_po_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_el_po_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_el_po_ijk,coupling_po_el_ijk
  integer, dimension(:), allocatable :: coupling_el_po_ispec,coupling_po_el_ispec
  integer :: num_coupling_el_po_faces

! material flag
  logical, dimension(:), allocatable :: ispec_is_poroelastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_poroelastic
  integer :: num_phase_ispec_poroelastic,nspec_inner_poroelastic,nspec_outer_poroelastic

  logical :: POROELASTIC_SIMULATION

! ADJOINT poroelastic

  ! (backward/reconstructed) wavefields
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_epsilonsdev_xx,b_epsilonsdev_yy,b_epsilonsdev_xy,b_epsilonsdev_xz,b_epsilonsdev_yz, &
    b_epsilonwdev_xx,b_epsilonwdev_yy,b_epsilonwdev_xy,b_epsilonwdev_xz,b_epsilonwdev_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_epsilons_trace_over_3,b_epsilonw_trace_over_3

  ! adjoint kernels [primary kernels, density kernels, wavespeed kernels]
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
    C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, Bb_kl, Cb_kl, Mb_kl, mufrb_kl, &
    rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl

  ! absorbing stacey wavefield parts
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_absorb_fields,b_absorb_fieldw
  integer :: b_reclen_field_poro

  ! for assembling backward field
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_meshs
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_meshw
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_meshs
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_meshw
  integer, dimension(:), allocatable :: b_request_send_vector_ext_meshs
  integer, dimension(:), allocatable :: b_request_send_vector_ext_meshw
  integer, dimension(:), allocatable :: b_request_recv_vector_ext_meshs
  integer, dimension(:), allocatable :: b_request_recv_vector_ext_meshw


end module specfem_par_poroelastic

!=====================================================================

module specfem_par_movie

! parameter module for movies/shakemovies

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD2D_FOUR_CORNERS

  implicit none

! to save full 3D snapshot of velocity (movie volume
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: div, curl_x, curl_y, curl_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: velocity_x,velocity_y,velocity_z

! shakemovies and movie surface
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
  integer :: nfaces_surface_ext_mesh,nfaces_surface_ext_mesh_points
  integer :: nfaces_surface_glob_ext_mesh,nfaces_surface_glob_em_points
  ! face corner indices
  integer :: iorderi(NGNOD2D_FOUR_CORNERS),iorderj(NGNOD2D_FOUR_CORNERS)

! movie parameters
  double precision :: HDUR_MOVIE
  integer :: NTSTEP_BETWEEN_FRAMES,MOVIE_TYPE
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
            USE_HIGHRES_FOR_MOVIES
  logical :: MOVIE_SIMULATION

end module specfem_par_movie
