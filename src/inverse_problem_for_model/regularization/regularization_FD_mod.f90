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

module regularization_fd_mod

!! IMPORT SPECFEM VARIABLES -----------------------------------------------------------------------------------------------------
  use specfem_par, only: NGLLX, NGLLY, NGLLZ, NDIM, NSPEC_AB, &
       NGLOB_AB, ibool, xstore, ystore, zstore,  NUM_ITER, NGNOD, xigll, yigll, zigll, NPROC, HUGEVAL, &
       MIDX, MIDY, MIDZ


  !! IMPORT INVERSE_PROBLEM VARIABLES ---------------------------------------------------------------------------------------------
  use inverse_problem_par
  use projection_on_FD_grid

  implicit none

  real(kind=CUSTOM_REAL), private, dimension(:,:,:), allocatable ::      model_on_FD_grid
  real(kind=CUSTOM_REAL), private, dimension(:,:,:), allocatable :: diff_model_on_FD_grid
  real(kind=CUSTOM_REAL), private, dimension(:,:,:), allocatable :: diff2_model_on_FD_grid

contains


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  Setup regularization
!--------------------------------------------------------------------------------------------------------------------
  subroutine setup_FD_regularization(projection_fd, myrank)

    type(profd),                            intent(inout)  :: projection_fd
    integer,                                intent(in)     :: myrank

    integer :: ier

    call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)
    allocate(model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 218')
    allocate(diff_model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 219')
    allocate(diff2_model_on_FD_grid(nx_fd_proj, ny_fd_proj, nz_fd_proj),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 220')

  end subroutine setup_FD_regularization
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  compute FD laplacian of model
!--------------------------------------------------------------------------------------------------------------------
  subroutine gradient_FD_laplac(model_on_SEM_mesh, regul_penalty, gradient_regul_penalty, cost_penalty, projection_fd, myrank, ipar)

    real(kind=CUSTOM_REAL),   dimension(:,:,:,:), allocatable, intent(inout)  :: regul_penalty, gradient_regul_penalty
    type(profd),                                               intent(in)     :: projection_fd
    integer,                                                   intent(in)     :: myrank, ipar
    real(kind=CUSTOM_REAL),   dimension(:,:,:,:), allocatable, intent(in)     :: model_on_SEM_mesh
    real(kind=CUSTOM_REAL),                                    intent(inout)  :: cost_penalty
    character(len=100)                                                        :: name_file

    !! put model from SEM mesh to FD grid
    call Project_model_SEM2FD_grid(model_on_SEM_mesh, model_on_FD_grid, projection_fd, myrank)

    !! compute FD laplacian and FD bi_laplacian
    call Compute_bi_laplac_FD(cost_penalty)

    !! put laplacian from FD to SEM mesh
    call Project_model_FD_grid2SEM(regul_penalty, diff_model_on_FD_grid, myrank)

    !! put bi_laplacian from FD to SEM mesh
    call Project_model_FD_grid2SEM(gradient_regul_penalty, diff2_model_on_FD_grid, myrank)

    if (DEBUG_MODE .and. myrank == 0) then

       write(name_file,'(a9,i3.3,a4)') "Model_FD_",ipar,".bin"
       open(676,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx_fd_proj*ny_fd_proj*nz_fd_proj)
       write(676,rec=1) model_on_FD_grid
       close(676)

       write(name_file,'(a10,i3.3,a4)') "Laplac_FD_",ipar,".bin"
       open(676,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx_fd_proj*ny_fd_proj*nz_fd_proj)
       write(676,rec=1) diff_model_on_FD_grid
       close(676)

       write(name_file,'(a13,i3.3,a4)') "Bi_Laplac_FD_",ipar,".bin"
       open(676,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx_fd_proj*ny_fd_proj*nz_fd_proj)
       write(676,rec=1) diff2_model_on_FD_grid
       close(676)

    endif

  end subroutine gradient_FD_laplac
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  compute FD bi-laplacian of model
!--------------------------------------------------------------------------------------------------------------------
  subroutine Compute_bi_laplac_FD(cp)
    real(kind=CUSTOM_REAL), intent(inout) :: cp
    real(kind=CUSTOM_REAL)                :: cp_dummy

    call LaplacFD2D(model_on_FD_grid, diff_model_on_FD_grid, cp)
    call LaplacFD2D(diff_model_on_FD_grid, diff2_model_on_FD_grid, cp_dummy)

  end subroutine Compute_bi_laplac_FD
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  compute FD laplacian of model
!--------------------------------------------------------------------------------------------------------------------
  subroutine LaplacFD2D(m, lm, cp)

    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(in)    :: m
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable, intent(inout) :: lm
    real(kind=CUSTOM_REAL),                                intent(inout) :: cp
    integer                                                              :: i,j,k

    cp = 0._CUSTOM_REAL
    lm(:,:,:)=0._CUSTOM_REAL
    do k = 2, nz_fd_proj-1
       do j = 2, ny_fd_proj-1
          do i = 2,  nx_fd_proj-1
             lm(i,j,k) = -6._CUSTOM_REAL * m(i,j,k)    + &
                    m(i+1, j,   k  )  +  m(i-1, j,  k  ) + &
                    m(i,   j+1, k  )  +  m(i , j-1, k  ) + &
                    m(i,   j,   k+1)  +  m(i,  j,   k-1)
             cp = cp + lm(i,j,k)**2
          enddo
       enddo
    enddo

  end subroutine LaplacFD2D


!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module regularization_fd_mod
