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

!##################################################################################################################################
!>
!!
!!
!!   main subroutine to launch the iterative Full Waveform Inversion.
!!
!!
!!
!!   -- the simulation is directly based on spedfem3D git devel version
!!   -- the MPI use directly the communicators defined in specfem
!!   -- additional module are used to perform an iterative full waveform inversion
!!   -- l-bfgs optimzation
!!
!!  V0.0.1
!!
!!  August 2016
!!  Vadim Monteiller, CNRS Marseille, France
!!
!#################################################################################################################################

subroutine inverse_problem_main()

  use specfem_par     !! module from specfem
  use inverse_problem_par
  use adjoint_source
  use input_output
  use regularization
  use inversion_scheme
  use projection_on_FD_grid
  use specfem_interface
  use fwi_iteration

  implicit none

  !! inversion global variables
  type(acqui), dimension(:), allocatable :: acqui_simu
!!!!!!!!!!!!!!!!!!  type(regul), dimension(:), allocatable :: regularization_fd
  type(inver)                            :: inversion_param
  type(profd)                            :: projection_fd

  integer                                :: ievent, iter_inverse, iter_frq
  logical                                :: finished
  character(len=MAX_LEN_STRING)          :: mode_running


  !!!##############################################################################################################################
  !!! ---------------------------------------------- INITIALIZE RUNTIME ----------------------------------------------------------
  !!!##############################################################################################################################

  !! get rank of my MPI process
  call world_rank(myrank)

  !! select which mode to run : only direct or FWI
  call get_mode_running(mode_running, inversion_param)

  !! open log files, distribute events if simultaneous simu, get acquisition file
  call SetUpInversion(inversion_param, myrank)

  !! intialize inversion (initialize_simulation read_parameter_file read_mesh_* prepare_time_run )
  call InitializeSpecfemForInversion()

  !! initialize input output (read inversion parameters input files, acquisition, get external model)
  call init_input_output_mod(inversion_param, acqui_simu, myrank)

  !! set up projection on fd grid if asked (for snapshot, movie, smoothing ...)
  if (PROJ_ON_FD) call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)

  !!!##############################################################################################################################
  !!! -------------------------------  different running mode : forward or FWI ----------------------------------------------------
  !!!##############################################################################################################################

  select case (trim(adjustl(mode_running)))

  case('forward') !--------------------------------------------------------------------------------------------------------

     if (myrank == 0) then
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) '--------------------------------------------------------------------------'
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)  '        Specfem inverse problem : Forward only problem ...  '
        write(INVERSE_LOG_FILE,*)  '        with generic family parameters  : ', trim(inversion_param%parameter_family_name)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
     endif

     do ievent = 1, acqui_simu(1)%nevent_tot
        call ComputeSismosPerEvent(ievent, acqui_simu, 0, inversion_param, myrank)
     enddo

     call synchronize_all_world()

     !! writing model in SEM mesh : (rho, vp, vs) or cijkl.
     if (inversion_param%output_model)  call WriteOuptutSEMmodel(inversion_param)

     !----------------------------------------------------------------------------------------------------------------

  case ('l-bfgs') !--------------------------------------------------------------------------------------------------------

     if (myrank == 0) then
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) '     --------------------------------------------------------------------------'
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)  '        Specfem inverse problem : L-BFGS FWI ...  '
        write(INVERSE_LOG_FILE,*)  '        with family parameters  : ', trim(inversion_param%parameter_family_name)
        write(INVERSE_LOG_FILE,*)  '        memory needed  ', inversion_param%max_history_bfgs
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
     endif

     !! initialize specifics arrays for optimization
     call AllocatememoryForFWI(inversion_param, acqui_simu(1)%nevent_tot)

     !! loop on frequencies
     do iter_frq = 1, inversion_param%Nifrq

        !! initialize (reset) optimization
        call InitializeOptimIteration(iter_frq, acqui_simu, inversion_param)

        !! loop on optimization iterations
        do iter_inverse = 0,  inversion_param%Niter

           call OneIterationOptim(iter_inverse, iter_frq, finished, acqui_simu, inversion_param) !!!!!! , regularization_fd)

           if (finished) exit

        enddo

        !! write model in disk
        call WriteOutputs(inversion_param)

     enddo

     !----------------------------------------------------------------------------------------------------------------
     if (myrank == 0) then
        if (iter_inverse >= inversion_param%Niter) then
            write(INVERSE_LOG_FILE,*)
            write(INVERSE_LOG_FILE,*) ' FWI STOP : maximum allowed iteration reached ', inversion_param%Niter
            write(INVERSE_LOG_FILE,*)
            iter_inverse = iter_inverse - 1
        endif
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) '  FWI STOPPED at iteration :', iter_inverse, ' / ', inversion_param%Niter
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) '--------------------------------------------------------------------------'
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
     endif

  case default

     if (myrank == 0) write(*,*) ' ERROR :', trim(mode_running),  ':  option not defined '

  end select


  !! SB
  close(IIDD)
  if (myrank == 0) then
     close(OUTPUT_ITERATION_FILE)
     close(INVERSE_LOG_FILE)
     close(OUTPUT_FWI_LOG)
  endif
  ! parce que sinon je dois mettre un rpint...
  call synchronize_all_world()

end subroutine inverse_problem_main

