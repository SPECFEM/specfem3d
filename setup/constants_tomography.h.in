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

!
!--- user can modify parameters below
!

!!-----------------------------------------------------------
!!
!! Model parameterization
!!
!!-----------------------------------------------------------
! by default, this algorithm uses isotropic (alpha,beta) kernels to sum up
! if you prefer using isotropic kernels, set flags below accordingly

  ! if you prefer using isotropic kernels (bulk,bulk_beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ISO_KERNELS = .false.

  ! if you prefer isotropic  (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .true.

  ! region code
  ! (local version "_", globe version "_reg1_" for crust/mantle
  character(len=8),parameter :: REG = "_"

!!-----------------------------------------------------------
!!
!! Directory structure
!!
!!-----------------------------------------------------------

  ! directory containing kernels from current iteration
  character(len=MAX_STRING_LEN) :: INPUT_KERNELS_DIR = 'INPUT_GRADIENT/'

  ! directory containing kernels from previous iteration
  character(len=MAX_STRING_LEN) :: KERNEL_OLD_DIR = './KERNELS/OUTPUT_SUM.old'

  ! directory containing model database files
  character(len=MAX_STRING_LEN) :: INPUT_MODEL_DIR = 'INPUT_MODEL/'

  ! directory containing external_mesh.bin database files
  character(len=MAX_STRING_LEN) :: INPUT_DATABASES_DIR = 'topo/'

  ! destination for newly updated model database files
  character(len=MAX_STRING_LEN) :: OUTPUT_MODEL_DIR = 'OUTPUT_MODEL/'

  ! destination for log files containing minimum, maximum vales
  character(len=MAX_STRING_LEN) :: OUTPUT_STATISTICS_DIR = 'OUTPUT_MODEL/'

  ! file containing absolute or relative paths to kernel directories
  character(len=MAX_STRING_LEN) :: KERNEL_FILE_LIST = './kernels_list.txt'

  ! produce logs containing minimum, maximum vales and other statistics
  logical :: PRINT_STATISTICS_FILES = .false.

  ! maximum number of kernel paths
  integer, parameter :: MAX_KERNEL_PATHS = 10000

!!-----------------------------------------------------------
!!
!! Scaling laws
!!
!!-----------------------------------------------------------
  ! ignore rho kernel, but use density perturbations as a scaling of Vs perturbations
  logical, parameter :: USE_RHO_SCALING = .true.

  ! in case of rho scaling, specifies density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL),parameter :: RHO_SCALING = 0.33_CUSTOM_REAL

!!-----------------------------------------------------------
!!
!! Transversely isotropic (TI) model constraints
!!
!!-----------------------------------------------------------
  ! constraint on eta model
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MIN = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL),parameter :: LIMIT_ETA_MAX = 1.5_CUSTOM_REAL

!!-----------------------------------------------------------
!!
!! Approximate Hessian
!!
!!-----------------------------------------------------------
  ! 1 permille of maximum for inverting Hessian
  real(kind=CUSTOM_REAL),parameter :: THRESHOLD_HESS = 1.e-3

  ! sums all Hessians before inverting and preconditioning
  ! by default should be set to .true.
  logical, parameter :: USE_HESS_SUM = .true.

!!-----------------------------------------------------------
!!
!! Maximum kernel scaling
!!
!!-----------------------------------------------------------
! kernel values are maximum at very shallow depth (due to receivers at surface) which leads to strong
! model updates closest to the surface. scaling the kernel values, such that the maximum is taken slightly below
! the surface (between 1km - 5km) leads to a "more balanced" gradient, i.e., a better model update in deeper parts

  ! by default, sets maximum update in this depth range
  logical,parameter :: USE_DEPTH_RANGE_MAXIMUM = .false.

  ! depths
  ! top at 1km depth
  real(kind=CUSTOM_REAL),parameter :: R_TOP = -1000.0 ! shallow depth
  ! bottom at 5km depth
  real(kind=CUSTOM_REAL),parameter :: R_BOTTOM = -5000.0 ! deep depth

!!-----------------------------------------------------------
!!
!! Source mask
!!
!!-----------------------------------------------------------
  ! uses source mask to blend out source elements
  logical, parameter :: USE_SOURCE_MASK = .false.


!!-----------------------------------------------------------
!!
!! Deprecated
!!
!!-----------------------------------------------------------

  ! determines separate step length for each material parameter in nonlinear conjugate gradient model update
  logical,parameter :: USE_SEPARATE_CG_STEPLENGTHS = .false.

