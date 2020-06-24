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


subroutine read_parameters_tomo()

! reads in parameters needed (only step length for now...)

  use tomography_par

  implicit none
  integer :: ier
  character(len=MAX_STRING_LEN) :: s_step_fac,arg

  ! subjective step length to multiply to the gradient
  ! e.g. step_fac = 0.03

  call get_command_argument(1,s_step_fac)

  if (trim(s_step_fac) == '') then
    call tomo_usage()
  endif

  ! read in parameter information
  read(s_step_fac,*,iostat=ier) step_fac
  if (ier /= 0) then
    call tomo_usage()
  endif

  ! safety check
  if (abs(step_fac) < 1.e-15) then
    print *,'Error: step factor ',step_fac,' is too small and will lead to no update...'
    call exit_MPI(myrank,'Error step factor too small')
  endif

  ! optional arguments
  ! input directory which holds (summed) kernel files
  call get_command_argument(2,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      INPUT_KERNELS_DIR = trim(arg) // '/'
    else
      INPUT_KERNELS_DIR = trim(arg)
    endif
  endif

  ! output directory for new model files
  call get_command_argument(3,arg)
  if (len_trim(arg) > 0) then
    if (arg(len_trim(arg):len_trim(arg)) /= '/') then
      OUTPUT_MODEL_DIR = trim(arg) // '/'
    else
      OUTPUT_MODEL_DIR = trim(arg)
    endif

    ! sets output directory for statistics to same directory
    if (PRINT_STATISTICS_FILES) then
      OUTPUT_STATISTICS_DIR = trim(OUTPUT_MODEL_DIR)
    endif
  endif

  ! statistics
  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_step_fac',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(OUTPUT_STATISTICS_DIR)//'statistics_step_fac'
      print *,'Please make sure that directory '//trim(OUTPUT_STATISTICS_DIR)//' exists...'
      print *
      stop 'Error opening statistics file'
    endif
    write(IOUT,'(1e24.12)') step_fac
    close(IOUT)
  endif

  end subroutine read_parameters_tomo

!
!-------------------------------------------------------------------------------------------------
!

  subroutine tomo_usage()

  use tomography_par

  implicit none

  if (myrank == 0) then
    print *,'Usage: add_model step_factor [INPUT-KERNELS-DIR/] [OUTPUT-MODEL-DIR/]'
    print *
    print *,'with'
    print *,'  step_factor        - factor to scale gradient (e.g. 0.03 for 3 percent update)'
    print *,'  INPUT-KERNELS-DIR/ - (optional) directory which holds summed kernels (e.g. alpha_kernel.bin,..), default ' &
            // trim(INPUT_KERNELS_DIR)
    print *,'  OUTPUT-MODEL-DIR/  - (optional) directory which will hold new model files (e.g. vp_new.bin,..), default ' &
            // trim(OUTPUT_MODEL_DIR)
    print *
    print *,'Please rerun e.g. like: mpirun -np ',sizeprocs,' ./bin/xadd_model 0.03'
    print *
  endif
  call synchronize_all()
  call exit_MPI(myrank,'Error usage: add_model step_factor')

  end subroutine tomo_usage


