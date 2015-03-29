!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
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


module inverse_problem_par

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN, &
    NGLLX,NGLLY,NGLLZ,IIN,IOUT, &
    FOUR_THIRDS,R_EARTH_KM,GAUSSALPHA,GAUSSBETA

  implicit none

  ! inverse_problem parameter settings
  include "constants_inverse_problem.h"

  ! mesh size
  integer :: NSPEC, NGLOB

  ! volume
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  integer, dimension(:,:,:,:),allocatable :: ibool

  ! model update length
  real(kind=CUSTOM_REAL) :: step_fac,step_length

  ! mpi process
  integer :: myrank,sizeprocs

  ! defaults
  ! directory which holds kernel files (alpha_kernel.bin,..)
  character(len=MAX_STRING_LEN) :: INPUT_KERNELS_DIR = 'INPUT_GRADIENT/'

  ! directory which holds model files (vp.bin,..)
  character(len=MAX_STRING_LEN) :: INPUT_MODEL_DIR = 'INPUT_MODEL/'

  ! directory which holds databases files (external_mesh.bin)
  character(len=MAX_STRING_LEN) :: INPUT_DATABASES_DIR = 'topo/'

  ! directory which holds new model files (vp_new.bin,..)
  character(len=MAX_STRING_LEN) :: OUTPUT_MODEL_DIR = 'OUTPUT_MODEL/'

  ! statistics
  ! set to true if you want to print out log files with statistics
  logical :: PRINT_STATISTICS_FILES = .false.
  ! directory where the statistics output files will be written
  character(len=MAX_STRING_LEN) :: OUTPUT_STATISTICS_DIR = 'OUTPUT_MODEL/'

end module inverse_problem_par

!
!-------------------------------------------------------------------------------------------------
!

module inverse_problem_kernels_iso

  use inverse_problem_par

  implicit none

  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_bulk,kernel_beta,kernel_rho

  ! gradients for model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dbulk,model_dbeta,model_drho

end module inverse_problem_kernels_iso


!
!-------------------------------------------------------------------------------------------------
!

module inverse_problem_kernels_tiso

  use inverse_problem_par

  implicit none

  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel_bulk,kernel_betav,kernel_betah,kernel_eta

  ! gradients for model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_dbulk,model_dbetah,model_dbetav,model_deta

end module inverse_problem_kernels_tiso


!
!-------------------------------------------------------------------------------------------------
!

module inverse_problem_kernels_tiso_cg

  use inverse_problem_par

  implicit none

  ! flags to determine wheter old gradients (model_dbulk,..) can be used or
  ! if update is based on old kernels only (kernel_bulk,..)
  logical :: USE_OLD_GRADIENT

  ! kernels from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        kernel_bulk_old,kernel_betav_old,kernel_betah_old,kernel_eta_old

  ! gradients for model updates from former iteration
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_dbulk_old,model_dbetah_old,model_dbetav_old,model_deta_old

end module inverse_problem_kernels_tiso_cg

!
!-------------------------------------------------------------------------------------------------
!

module inverse_problem_model_iso

  use inverse_problem_par

  implicit none

  ! isotropic model files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp,model_vs,model_rho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model_vp_new,model_vs_new,model_rho_new

end module inverse_problem_model_iso

!
!-------------------------------------------------------------------------------------------------
!

module inverse_problem_model_tiso

  use inverse_problem_par

  implicit none

  ! transverse isotropic model files
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_vpv,model_vph,model_vsv,model_vsh,model_eta,model_rho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
        model_vpv_new,model_vph_new,model_vsv_new,model_vsh_new,model_eta_new,model_rho_new

end module inverse_problem_model_tiso


