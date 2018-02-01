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
!
! United States and French Government Sponsorship Acknowledged.

module constants_meshfem3D

  include "constants_meshfem3D.h"

end module constants_meshfem3D

!
!----------------------------------------------------------------------------------------------------
!

! main parameter module for xmeshfem3D mesher

module meshfem3D_par

  use constants

  use constants_meshfem3D

  use shared_parameters

  implicit none

! number of spectral elements in each block
  integer :: nspec,npointot

! meshing parameters
  double precision, dimension(:), allocatable :: rns

  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid

  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer :: myrank,sizeprocs

! mesh point steps for interfaces
  integer :: npx_element_steps,npy_element_steps

! for loop on all the slices
  integer, dimension(:,:), allocatable :: addressing

! addressing for all the slices
  integer, dimension(:), allocatable :: iproc_xi_slice,iproc_eta_slice
  integer :: iproc_xi_current,iproc_eta_current

! parameters read from mesh parameter file
  integer :: NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA

  double precision :: UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX
  double precision :: Z_DEPTH_BLOCK
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  logical :: USE_REGULAR_MESH

! Mesh files for visualization
  logical :: CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

! CPML
  double precision :: THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML
  logical, dimension(:), allocatable :: is_CPML
  integer, dimension(:), allocatable :: CPML_to_spec,CPML_regions
  integer :: nspec_CPML

! doublings parameters
  integer :: NDOUBLINGS
  integer, dimension(:),allocatable :: ner_doublings

! parameters deduced from parameters read from file
  integer :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer :: NER

! this for all the regions
  integer :: NSPEC_AB,NGLOB_AB
  integer :: NSPEC2D_A_XI,NSPEC2D_B_XI, &
             NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
             NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
             NSPEC2D_BOTTOM,NSPEC2D_TOP, &
             NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

! interfaces parameters
  integer :: number_of_interfaces,number_of_layers
  integer :: max_npx_interface,max_npy_interface

  character(len=MAX_STRING_LEN) :: INTERFACES_FILE

  integer, dimension(:), allocatable :: ner_layer

! cavity parameters
  character(len=MAX_STRING_LEN) :: CAVITY_FILE

! subregions parameters
  integer :: NSUBREGIONS
!  definition of the different regions of the model in the mesh (nx,ny,nz)
!  #1 #2 : nx_begining,nx_end
!  #3 #4 : ny_begining,ny_end
!  #5 #6 : nz_begining,nz_end
!     #7 : material number
  integer, dimension(:,:), pointer :: subregions

! material properties
  integer :: NMATERIALS
! first dimension  : material_id
! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision , dimension(:,:), pointer :: material_properties

  logical :: BROADCAST_AFTER_READ

end module meshfem3D_par


