!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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


module constants_meshfem

  implicit none

  include "constants_meshfem3D.h"

end module constants_meshfem

!
!----------------------------------------------------------------------------------------------------
!

! main parameter module for xmeshfem3D mesher

module meshfem_par

  use constants

  use constants_meshfem

  use shared_parameters

  implicit none

  ! number of spectral elements in each block
  integer :: nspec = 0

  ! meshing parameters
  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid
  integer, dimension(:,:,:,:), allocatable :: ibool

  ! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

  ! proc numbers for MPI
  integer :: sizeprocs = 0

  ! mesh point steps for interfaces
  integer :: npx_element_steps,npy_element_steps

  ! for loop on all the slices
  integer, dimension(:,:), allocatable :: addressing

  ! addressing for all the slices
  integer, dimension(:), allocatable :: iproc_xi_slice,iproc_eta_slice
  integer :: iproc_xi_current,iproc_eta_current

  ! parameters read from mesh parameter file
  integer :: NEX_XI = 0
  integer :: NEX_ETA = 0
  integer :: NPROC_XI = 0
  integer :: NPROC_ETA = 0

  double precision :: UTM_X_MIN = 0.d0
  double precision :: UTM_X_MAX = 0.d0
  double precision :: UTM_Y_MIN = 0.d0
  double precision :: UTM_Y_MAX = 0.d0
  double precision :: Z_DEPTH_BLOCK = 0.d0
  double precision :: LATITUDE_MIN = 0.d0
  double precision :: LATITUDE_MAX = 0.d0
  double precision :: LONGITUDE_MIN = 0.d0
  double precision :: LONGITUDE_MAX = 0.d0

  logical :: USE_REGULAR_MESH = .false.

  ! Mesh files for visualization
  logical :: CREATE_ABAQUS_FILES = .false.
  logical :: CREATE_DX_FILES = .false.
  logical :: CREATE_VTK_FILES = .false.

  ! for Cubit postprocessing
  logical :: SAVE_MESH_AS_CUBIT = .false.

  ! CPML
  double precision :: THICKNESS_OF_X_PML = 0.d0
  double precision :: THICKNESS_OF_Y_PML = 0.d0
  double precision :: THICKNESS_OF_Z_PML = 0.d0
  logical :: ADD_PML_AS_EXTRA_MESH_LAYERS = .false.  ! for mesh extension with PML layers
  integer :: NUMBER_OF_PML_LAYERS_TO_ADD = 0
  logical, dimension(:), allocatable :: is_CPML
  integer, dimension(:), allocatable :: CPML_to_spec,CPML_regions
  integer :: nspec_CPML

  ! doublings parameters
  integer :: NDOUBLINGS = 0
  integer, dimension(:),allocatable :: ner_doublings

  ! parameters deduced from parameters read from file
  integer :: NEX_PER_PROC_XI = 0
  integer :: NEX_PER_PROC_ETA = 0
  integer :: NER = 0

  ! this for all the regions
  integer :: NSPEC_AB = 0
  integer :: NGLOB_AB = 0

  integer :: NSPEC2D_A_XI = 0
  integer :: NSPEC2D_B_XI = 0
  integer :: NSPEC2D_A_ETA = 0
  integer :: NSPEC2D_B_ETA = 0
  integer :: NSPEC2DMAX_XMIN_XMAX = 0
  integer :: NSPEC2DMAX_YMIN_YMAX = 0
  integer :: NSPEC2D_BOTTOM = 0
  integer :: NSPEC2D_TOP = 0

  !integer :: NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX - not needed...

  ! interfaces parameters
  integer :: number_of_interfaces = 0
  integer :: number_of_layers = 0
  integer :: max_npx_interface,max_npy_interface

  character(len=MAX_STRING_LEN) :: INTERFACES_FILE = 'interfaces.dat'

  integer, dimension(:), allocatable :: ner_layer

  ! cavity parameters
  character(len=MAX_STRING_LEN) :: CAVITY_FILE = 'dummy'

  ! subregions parameters
  integer :: NSUBREGIONS = 0
  !  definition of the different regions of the model in the mesh (nx,ny,nz)
  !  #1 #2 : nx_begining,nx_end
  !  #3 #4 : ny_begining,ny_end
  !  #5 #6 : nz_begining,nz_end
  !     #7 : material number
  integer, dimension(:,:), allocatable :: subregions

  ! material properties
  integer :: NMATERIALS = 0
  ! first dimension  : material_id
  ! second dimension : #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id ..
  double precision , dimension(:,:), allocatable :: material_properties
  ! tomography materials
  character(len=MAX_STRING_LEN), dimension(:,:), allocatable :: material_properties_undef

  logical :: BROADCAST_AFTER_READ = .false.

  ! name of the database file
  character(len=MAX_STRING_LEN) :: prname

  !! boundary of wavefield discontinuity, read from database file
  integer :: nb_wd

  !! boundary_to_ispec_wd(nb_wd)
  !! the element the boundary belongs to, read from database file
  !! each point on the boundary belongs to two sides of the boundary
  !! here the element must be on the inner side of the boundary
  integer, dimension(:), allocatable :: boundary_to_ispec_wd

  !! side_wd(nb_wd)
  !! integers specifying which side the boundary is in the element
  !! read from database file
  !! side_wd = 1--8: only one vertex is on the boundary
  !! side_wd = 9--20: only one edge is on the boundary
  !! side_wd = 21--26: one face is on the boundary
  integer, dimension(:), allocatable :: side_wd

end module meshfem_par

