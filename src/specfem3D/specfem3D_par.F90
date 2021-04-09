!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


module specfem_par

! main parameter module for specfem simulations

  use constants

  use shared_parameters

  implicit none

!-----------------------------------------------------------------
! simulation
!-----------------------------------------------------------------

  ! number of spectral element and global points
  integer :: NSPEC_AB, NGLOB_AB

  ! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,ystore,zstore

  ! regular/irregular element shapes
  integer :: NSPEC_IRREGULAR
  integer, dimension(:), allocatable :: irregular_element_number
  real(kind=CUSTOM_REAL) :: xix_regular,jacobian_regular

  ! derivatives (of mapping to reference element)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    deriv_mapping

  ! material properties
  ! isotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kappastore,mustore

  ! density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore

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

  ! attenuation
  integer :: NSPEC_ATTENUATION_AB
  character(len=MAX_STRING_LEN) :: prname_Q

  ! additional mass matrix for ocean load
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

!-----------------------------------------------------------------
! coupling
!-----------------------------------------------------------------

  ! for couple with external code : DSM and AxiSEM (added by VM) for the moment
  integer :: it_dsm, it_fk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: Veloc_dsm_boundary, Tract_dsm_boundary
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: Veloc_axisem, Tract_axisem

  !! CD CD added this for RECIPROCITY_AND_KH_INTEGRAL
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: Displ_axisem_time, Tract_axisem_time
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: Displ_specfem_time, Tract_specfem_time


!-----------------------------------------------------------------
! time scheme
!-----------------------------------------------------------------

  ! time scheme
  real(kind=CUSTOM_REAL) :: deltat,deltatover2,deltatsqover2
  ! backward/reconstructed
  real(kind=CUSTOM_REAL) :: b_deltat, b_deltatover2, b_deltatsqover2

  ! LDDRK time scheme
  integer :: NSTAGE_TIME_SCHEME,istage
  integer :: NGLOB_AB_LDDRK,NSPEC_ATTENUATION_AB_LDDRK

  ! time loop step
  integer :: it

!! DK DK added this temporarily here to make SPECFEM3D and SPECFEM3D_GLOBE much more similar
!! DK DK in terms of the structure of their main time iteration loop; these are future features
!! DK DK that are missing in this code but implemented in the other and that could thus be cut and pasted one day
  integer :: it_begin,it_end

  ! UNDO_ATTENUATION_AND_OR_PML
  integer :: NSUBSET_ITERATIONS
  integer :: iteration_on_subset,it_of_this_subset
  integer :: it_subset_end

!-----------------------------------------------------------------
! sources
!-----------------------------------------------------------------

  ! parameters for the source
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays

  ! CMT moment tensors
  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  double precision, dimension(:,:,:), allocatable :: nu_source
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source

  double precision, dimension(:), allocatable :: tshift_src,hdur,hdur_Gaussian
  double precision, dimension(:), allocatable :: utm_x_source,utm_y_source
  double precision :: t0

  ! SAC time
  integer :: yr,mo,da,jda,ho,mi
  double precision :: sec

  ! source time function
  real(kind=CUSTOM_REAL) :: stf_used_total
  integer :: nsources_local
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: user_source_time_function

  ! source encoding
  ! for acoustic sources: takes +/- 1 sign, depending on sign(Mxx)
  ! [ = sign(Myy) = sign(Mzz) since they have to be equal in the acoustic setting]
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pm1_source_encoding

  ! parameters for a force source located exactly at a grid point
  double precision, dimension(:), allocatable :: factor_force_source
  double precision, dimension(:), allocatable :: comp_dir_vect_source_E
  double precision, dimension(:), allocatable :: comp_dir_vect_source_N
  double precision, dimension(:), allocatable :: comp_dir_vect_source_Z_UP

!-----------------------------------------------------------------
! receivers
!-----------------------------------------------------------------

  ! receiver information
  integer :: nrec,nrec_local
  integer :: nrec_tot_found
  character(len=MAX_STRING_LEN) :: rec_filename,filtered_rec_filename
  integer, dimension(:), allocatable :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store
  double precision, dimension(:,:,:), allocatable :: nu_rec

  ! Lagrange interpolators at receivers
  integer, dimension(:), allocatable, target :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable, target :: hxir_store,hetar_store,hgammar_store

  ! adjoint sources
  integer :: nadj_rec_local
  integer, dimension(:), pointer :: number_adjsources_global
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: hxir_adjstore,hetar_adjstore,hgammar_adjstore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_adjoint

  ! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

  ! seismograms
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms_d,seismograms_v,seismograms_a,seismograms_p
  integer :: nlength_seismogram

  integer :: seismo_offset,seismo_current
  ! adjoint seismograms
  integer :: it_adj_written

!-----------------------------------------------------------------
! GLL points & weights
!-----------------------------------------------------------------

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprime_yyT,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprime_zzT,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  ! arrays for Deville and force_vectorization
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  ! proc numbers for MPI
  integer :: sizeprocs
  character(len=MAX_STRING_LEN) :: prname

  ! timer MPI
  double precision :: time_start

  ! array for NB_RUN_ACOUSTIC_GPU > 1
  integer, dimension(:), allocatable :: run_number_of_the_source

  ! for assembling in case of external mesh
  integer :: num_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer, dimension(:), allocatable :: my_neighbors_ext_mesh
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
  integer :: nfaces_surface

  ! MPI partition surfaces
  logical, dimension(:), allocatable :: ispec_is_inner

  ! maximum speed in velocity model
  real(kind=CUSTOM_REAL):: model_speed_max

  ! gravity
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: minus_deriv_gravity,minus_g
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! for surface or volume integral on whole domain
  double precision, dimension(:), allocatable   :: integral_vol, integral_boun
  double precision, dimension(:,:), allocatable :: f_integrand_KH

!-----------------------------------------------------------------
! adjoint simulations
!-----------------------------------------------------------------

  ! absorbing stacey wavefield parts
  integer :: b_num_abs_boundary_faces
  logical :: SAVE_STACEY

  ! Moho mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC_BOUN,NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot

  ! adjoint source frechet derivatives
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: seismograms_eps

  ! adjoint elements
  integer :: NSPEC_ADJOINT, NGLOB_ADJOINT

!-----------------------------------------------------------------
! gravity
!-----------------------------------------------------------------

  ! for gravity integrals
  double precision, dimension(NTOTAL_OBSERVATION) :: x_observation,y_observation,z_observation, &
                                                g_x,g_y,g_z,G_xx,G_yy,G_zz,G_xy,G_xz,G_yz,temporary_array_for_sum

  ! force vectorization
#ifdef FORCE_VECTORIZATION
  logical, parameter :: FORCE_VECTORIZATION_VAL = .true.
#else
  logical, parameter :: FORCE_VECTORIZATION_VAL = .false.
#endif

#ifdef VTK_VIS
  ! VTK window mode, default is off
  logical :: VTK_MODE = .false.
#endif

!-----------------------------------------------------------------
! point search
!-----------------------------------------------------------------

  ! point search
  ! (i,j,k) indices of the control/anchor points of the element
  integer, dimension(:), allocatable :: anchor_iax,anchor_iay,anchor_iaz

  ! coordinates of element midpoints
  double precision, dimension(:,:), allocatable :: xyz_midpoints

!-----------------------------------------------------------------
! GPU
!-----------------------------------------------------------------

  ! CUDA mesh pointer to integer wrapper
  integer(kind=8) :: Mesh_pointer

  ! for dynamic rupture computations on GPU
  integer(kind=8) :: Fault_pointer

  ! ASDF
  ! asdf file handle
  integer :: current_asdf_handle

end module specfem_par

!=====================================================================

module specfem_par_elastic

! parameter module for elastic solver

  use constants, only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ

  implicit none

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: factor_common,factor_common_kappa
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_sigma
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL) :: min_resolved_period

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: R_trace,R_xx,R_yy,R_xy,R_xz,R_yz
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
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store, &
            c34store,c35store,c36store,c44store,c45store,c46store, &
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
  integer :: iglob_check_elastic

  ! ADJOINT elastic

  ! (backward/reconstructed) wavefields
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ, b_veloc, b_accel

  ! backward attenuation arrays
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: b_alphaval, b_betaval, b_gammaval
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy,b_epsilondev_xz,b_epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_epsilon_trace_over_3

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl

  ! anisotropic kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: cijkl_kl

  ! approximate Hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_kl, hess_rho_kl, hess_mu_kl, hess_kappa_kl

  ! topographic (Moho) kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:),allocatable :: dsdx_top, dsdx_bot, b_dsdx_top, b_dsdx_bot
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: moho_kl
  integer, dimension(:), allocatable :: ispec2D_moho_top,ispec2D_moho_bot

  ! absorbing stacey wavefield parts
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_absorb_field
  integer :: b_reclen_field

  ! for assembling backward field
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_mesh
  integer, dimension(:), allocatable :: b_request_send_vector_ext_mesh
  integer, dimension(:), allocatable :: b_request_recv_vector_ext_mesh

  ! LDDRK time scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_lddrk,veloc_lddrk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_trace_lddrk,R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk
  ! adjoint
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_lddrk, b_veloc_lddrk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    b_R_trace_lddrk,b_R_xx_lddrk,b_R_yy_lddrk,b_R_xy_lddrk,b_R_xz_lddrk,b_R_yz_lddrk


end module specfem_par_elastic

!=====================================================================

module specfem_par_acoustic

! parameter module for acoustic solver

  use constants, only: CUSTOM_REAL
  implicit none

  ! potential
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic,potential_dot_acoustic, &
                                    potential_dot_dot_acoustic
  !real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_adj_coupling ! not used yet

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
  integer :: iglob_check_acoustic

  ! ADJOINT acoustic

  ! (backward/reconstructed) wavefield potentials
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_acoustic, &
                        b_potential_dot_acoustic,b_potential_dot_dot_acoustic
  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl, &
    rhop_ac_kl, alpha_ac_kl

  ! approximate Hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_ac_kl, hess_rho_ac_kl, hess_kappa_ac_kl

  ! absorbing stacey wavefield parts
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_absorb_potential
  integer :: b_reclen_potential

  ! for assembling backward field
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_recv_scalar_ext_mesh
  integer, dimension(:), allocatable :: b_request_send_scalar_ext_mesh
  integer, dimension(:), allocatable :: b_request_recv_scalar_ext_mesh

  ! LDDRK time scheme
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_lddrk,potential_dot_acoustic_lddrk

end module specfem_par_acoustic

!=====================================================================

module specfem_par_poroelastic

! parameter module for elastic solver

  use constants, only: CUSTOM_REAL
  implicit none

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_solid_poroelastic, &
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
  integer :: iglob_check_poroelastic

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

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD2D_FOUR_CORNERS

  implicit none

  ! to save full 3D snapshot of velocity (movie volume
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: div, curl_x, curl_y, curl_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: velocity_x,velocity_y,velocity_z

  ! surface point locations
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_y
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_z
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_y_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_z_all

  ! movie data
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uz
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uy_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_uz_all

  ! shakemovie data
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_ux
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_uy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_uz
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_ux_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_uy_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shakemap_uz_all

  ! for storing surface of external mesh
  integer,dimension(:),allocatable :: nfaces_perproc_surface
  integer,dimension(:),allocatable :: faces_surface_offset
  integer,dimension(:,:),allocatable :: faces_surface_ibool
  integer,dimension(:),allocatable :: faces_surface_ispec
  integer :: nfaces_surface_points
  integer :: nfaces_surface_glob_ext_mesh,nfaces_surface_glob_points

  ! movie parameters
  logical :: MOVIE_SIMULATION

end module specfem_par_movie

!=====================================================================

module specfem_par_coupling

  use constants, only: CUSTOM_REAL

  implicit none

! added by Ping Tong (TP / Tong Ping) for the FK3D calculation

  ! FK elastic
  integer :: npt,nlayer,kpsv
  integer :: NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_FFT, NP_RESAMP, NPOW_FOR_INTERP
  integer :: NPTS_STORED, NPTS_INTERP

  integer, parameter :: NTIME_BETWEEN_FFT = 1  !! not used anymore

  integer,dimension(:),allocatable :: nbdglb
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: vxbd,vybd,vzbd,txxbd,txybd,txzbd,tyybd,tyzbd,tzzbd
  real(kind=CUSTOM_REAL) :: Z_REF_for_FK

  ! source
  real(kind=CUSTOM_REAL) :: p,xx0,yy0,zz0,ff0,tg,tt0,tmax_fk,df_fk
  real(kind=CUSTOM_REAL) :: phi_FK,theta_FK

  ! model
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: al_FK,be_FK,mu_FK,h_FK,tmp_for_interp
  complex(kind=8), dimension(:,:), allocatable :: VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t

  complex(kind=8), dimension(:), allocatable :: WKS_CMPLX_FOR_FFT
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: WKS_REAL_FOR_FFT

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
  logical :: stag
  real(kind=CUSTOM_REAL),dimension(:),allocatable  :: xx,yy,zz,xi1,xim,bdlambdamu
  ! normal
  real(kind=CUSTOM_REAL),dimension(:),allocatable  :: nmx,nmy,nmz

end module specfem_par_coupling

!=====================================================================

module specfem_par_noise

! parameter module for noise simulations

  use constants, only: CUSTOM_REAL

  implicit none

  ! NOISE_TOMOGRAPHY
  ! parameter module for noise simulations
  integer :: irec_main_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: noise_surface_movie
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

  ! noise strength kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sigma_kl

end module specfem_par_noise
