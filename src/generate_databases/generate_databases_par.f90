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


  module generate_databases_par

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer npointot

! local to global indexing array
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer :: myrank,sizeprocs,ier

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION

  integer :: NX_TOPO,NY_TOPO
  integer, dimension(:,:), allocatable :: itopo_bathy

! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

! parameters read from parameter file
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE
  integer :: NSOURCES,NGNOD,NGNOD2D,MOVIE_TYPE
  integer :: NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_READ_ADJSRC

  double precision :: DT,HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML

  logical :: ATTENUATION,USE_OLSEN_ATTENUATION,APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,SAVE_FORWARD,USE_FORCE_POINT_SOURCE
  logical :: ANISOTROPY,STACEY_ABSORBING_CONDITIONS,SAVE_MESH_FILES,STACEY_INSTEAD_OF_FREE_SURFACE
  logical :: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID
  logical :: USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES

  character(len=256) OUTPUT_FILES,LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH

  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
             ADIOS_FOR_KERNELS

! parameters deduced from parameters read from file
  integer :: NPROC

! memory size that will be needed by the solver
  double precision :: max_memory_size,max_memory_size_request

! this for all the regions
  integer NSPEC_AB,NGLOB_AB

  integer NSPEC2D_BOTTOM,NSPEC2D_TOP

  double precision min_elevation,max_elevation
  double precision min_elevation_all,max_elevation_all

! for Databases of external meshes
  double precision, dimension(:,:), allocatable :: nodes_coords_ext_mesh

  integer :: dummy_node
  integer :: dummy_elmnt

  integer :: ispec, inode, num_interface,ie,imat,iface,icorner
  integer :: nnodes_ext_mesh, nelmnts_ext_mesh
  integer  :: num_interfaces_ext_mesh
  integer  :: max_interface_size_ext_mesh
  integer  :: nmat_ext_mesh, nundefMat_ext_mesh
  integer, dimension(:), allocatable  :: my_neighbours_ext_mesh
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours_ext_mesh

  integer, dimension(:,:,:), allocatable  :: my_interfaces_ext_mesh
  integer, dimension(:,:), allocatable  :: ibool_interfaces_ext_mesh
  integer, dimension(:), allocatable  :: nibool_interfaces_ext_mesh

  integer, dimension(:,:), allocatable :: elmnts_ext_mesh
  integer, dimension(:,:), allocatable :: mat_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy

  character(len=256) prname

! boundaries and materials
  double precision, dimension(:,:), allocatable :: materials_ext_mesh

  integer :: ispec2D, boundary_number
  integer :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom_ext, nspec2D_top_ext

  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
              ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable :: nodes_ibelm_xmin,nodes_ibelm_xmax, &
              nodes_ibelm_ymin, nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top

  character (len=30), dimension(:,:), allocatable :: undef_mat_prop

! C-PML absorbing boundary conditions

  ! local number of C-PML spectral elements
  integer :: nspec_cpml

  ! global number of C-PML spectral elements
  integer :: nspec_cpml_tot

  ! C-PML spectral elements global indexing
  integer, dimension(:), allocatable :: CPML_to_spec

  ! C-PML regions
  integer, dimension(:), allocatable :: CPML_regions

  ! mask of C-PML elements for the global mesh
  logical, dimension(:), allocatable :: is_CPML

  ! thickness of C-PML layers in each direction
  real(kind=CUSTOM_REAL) :: CPML_width_x,CPML_width_y,CPML_width_z

  ! C-PML damping profile arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: d_store_x, d_store_y, d_store_z

  ! auxiliary parameters arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: K_store_x, K_store_y, K_store_z
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: alpha_store_x,alpha_store_y,alpha_store_z

  ! array recording the points on interface shared by PML and interior computational domain
  logical, dimension(:), allocatable :: mask_ibool_interior_domain
  integer :: nglob_interface_PML_acoustic,nglob_interface_PML_elastic
  integer, dimension(:), allocatable :: points_interface_PML_acoustic, points_interface_PML_elastic

! moho (optional)
  integer :: nspec2D_moho_ext
  integer, dimension(:), allocatable  :: ibelm_moho
  integer, dimension(:,:), allocatable  :: nodes_ibelm_moho

  integer :: nglob_total,nspec_total

  logical,dimension(:),allocatable :: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh
  integer :: nfaces_surface_ext_mesh,nfaces_surface_glob_ext_mesh

! flag for noise simulation
  integer :: NOISE_TOMOGRAPHY
  integer :: IMODEL

  end module generate_databases_par

!
!-------------------------------------------------------------------------------------------------
!

  module create_regions_mesh_ext_par

  include 'constants.h'

! global point coordinates
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy
  integer :: nglob_dummy

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision, dimension(:), allocatable :: xelm,yelm,zelm

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! for model density, kappa, mu
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore

! for poroelastic model
  real(kind=CUSTOM_REAL),dimension(:,:,:,:), allocatable :: etastore,phistore,tortstore
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: rhoarraystore,kappaarraystore,permstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vpI,rho_vpII,rho_vsI

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass,rmass_acoustic,&
                            rmass_solid_poroelastic,rmass_fluid_poroelastic

! mass matrix contributions
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassy,rmassz
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassz_acoustic
  integer :: nglob_xy

! ocean load
  integer :: NGLOB_OCEAN
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qkappa_attenuation_store,qmu_attenuation_store

! 2D shape functions and their derivatives, weights
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top
  double precision, dimension(:,:), allocatable :: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

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

! elastic-poroelastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_el_po_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_el_po_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_el_po_ijk,coupling_po_el_ijk
  integer, dimension(:), allocatable :: coupling_el_po_ispec,coupling_po_el_ispec
  integer :: num_coupling_el_po_faces

  ! Moho mesh
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot

! for stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! anisotropy
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
            c11store,c12store,c13store,c14store,c15store,c16store,&
            c22store,c23store,c24store,c25store,c26store,c33store,&
            c34store,c35store,c36store,c44store,c45store,c46store,&
            c55store,c56store,c66store

! material domain flags
  logical, dimension(:), allocatable :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic

! name of the database file
  character(len=256) prname

! inner/outer elements
  logical,dimension(:),allocatable :: ispec_is_inner
  integer :: nspec_inner_acoustic,nspec_outer_acoustic
  integer :: nspec_inner_elastic,nspec_outer_elastic
  integer :: nspec_inner_poroelastic,nspec_outer_poroelastic

  integer :: num_phase_ispec_acoustic
  integer,dimension(:,:),allocatable :: phase_ispec_inner_acoustic

  integer :: num_phase_ispec_elastic
  integer,dimension(:,:),allocatable :: phase_ispec_inner_elastic

  integer :: num_phase_ispec_poroelastic
  integer,dimension(:,:),allocatable :: phase_ispec_inner_poroelastic

  logical :: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION

  ! mesh coloring
  integer :: num_colors_outer_acoustic,num_colors_inner_acoustic
  integer, dimension(:), allocatable :: num_elem_colors_acoustic

  integer :: num_colors_outer_elastic,num_colors_inner_elastic
  integer, dimension(:), allocatable :: num_elem_colors_elastic

  end module create_regions_mesh_ext_par


