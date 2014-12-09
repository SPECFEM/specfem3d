!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  ! a negative initial value is a convention that indicates that groups (i.e. sub-communicators, one per run) are off by default
  integer :: mygroup = -1

  ! create a copy of the original output file path, to which we may add a "run0001/", "run0002/", "run0003/" prefix later
  ! if NUMBER_OF_SIMULTANEOUS_RUNS > 1
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES_PATH = OUTPUT_FILES_PATH_BASE

  ! if doing simultaneous runs for the same mesh and model, see who should read the mesh and the model and broadcast it to others
  ! we put a default value here
  logical :: I_should_read_the_database = .true.

end module constants

!
!-------------------------------------------------------------------------------------------------
!

  module shared_input_parameters

! holds input parameters given in DATA/Par_file

  use constants,only: MAX_STRING_LEN

  implicit none

  ! parameters read from parameter file
  integer :: NPROC

  ! simulation parameters
  integer :: SIMULATION_TYPE
  integer :: NOISE_TOMOGRAPHY
  logical :: SAVE_FORWARD

  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION

  ! number of time steps
  integer :: NSTEP
  double precision :: DT

  integer :: NGNOD
!  character(len=MAX_STRING_LEN) :: MODEL

  character(len=MAX_STRING_LEN) :: SEP_MODEL_DIRECTORY

  ! physical parameters
  logical :: APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ATTENUATION,FULL_ATTENUATION_SOLID,ANISOTROPY
  logical :: GRAVITY

  character(len=MAX_STRING_LEN) :: TOMOGRAPHY_PATH

  ! attenuation
  logical :: USE_OLSEN_ATTENUATION
  double precision :: OLSEN_ATTENUATION_RATIO

  ! absorbing boundaries
  logical :: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE
  double precision :: f0_FOR_PML

  logical :: STACEY_ABSORBING_CONDITIONS,STACEY_INSTEAD_OF_FREE_SURFACE

  ! for simultaneous runs from the same batch job
  integer :: NUMBER_OF_SIMULTANEOUS_RUNS
  logical :: BROADCAST_SAME_MESH_AND_MODEL

  ! movies
  logical :: CREATE_SHAKEMAP
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  integer :: MOVIE_TYPE
  integer :: NTSTEP_BETWEEN_FRAMES
  double precision :: HDUR_MOVIE

  ! mesh
  logical :: SAVE_MESH_FILES
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! seismograms
  integer :: NTSTEP_BETWEEN_OUTPUT_INFO
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC
  logical :: SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY,SAVE_SEISMOGRAMS_ACCELERATION
  logical :: WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_SEISMOGRAMS,SU_FORMAT

  ! sources
  logical :: USE_FORCE_POINT_SOURCE
  logical :: USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION

  ! external code coupling (DSM)
  logical :: COUPLE_WITH_EXTERNAL_CODE
  integer :: EXTERNAL_CODE_TYPE
  character(len=MAX_STRING_LEN) :: TRACTION_PATH
  logical :: MESH_A_CHUNK_OF_THE_EARTH

  ! GPU simulations
  logical :: GPU_MODE

  ! adios file output
  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, ADIOS_FOR_KERNELS

  end module shared_input_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_compute_parameters

  ! parameters to be computed based upon parameters above read from file

  implicit none

  ! number of sources given in CMTSOLUTION file
  integer :: NSOURCES

  ! anchor points
  integer :: NGNOD2D

  ! model
  integer :: IMODEL


  end module shared_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_parameters

  use shared_input_parameters
  use shared_compute_parameters

  implicit none

  end module shared_parameters

