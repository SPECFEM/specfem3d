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
! United States and French Government Sponsorship Acknowledged.

!##################################################################################################################################
!>
!!   main subroutine to launch the iterative Full Waveform Inversion.
!!
!!   -- the simulation is directly based on specfem3D git devel version
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

  integer                                :: ievent
  logical                                :: finished
  character(len=MAX_LEN_STRING)          :: mode_running
  ! timing
  double precision                       :: tCPU,tstart,tstart_begin
  integer                                :: ihours,iminutes,iseconds,int_tCPU
  double precision,             external :: wtime

  !!!##############################################################################################################################
  !!! ---------------------------------------------- INITIALIZE RUNTIME ----------------------------------------------------------
  !!!##############################################################################################################################

  ! get MPI starting time
  tstart_begin = wtime()

  !! get rank of my MPI process
  call world_rank(myrank)

  ! header output
  if (myrank == 0) call header_output(6)

  !! select which mode to run : only direct or FWI
  call get_mode_running(mode_running, inversion_param)

  !! open log files, distribute events if simultaneous simu, get acquisition file
  call SetUpInversion(inversion_param, myrank)

  !! intialize inversion (initialize_simulation read_parameter_file read_mesh_* prepare_time_run )
  call InitializeSpecfemForInversion(inversion_param)

  !! initialize input output (read inversion parameters input files, acquisition, get external model)
  call init_input_output_mod(inversion_param, acqui_simu, myrank)

  !! set up projection on fd grid if asked (for snapshot, movie, smoothing ...)
  if (PROJ_ON_FD) call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)

  ! elapsed time since beginning
  if (myrank == 0) then
    tCPU = wtime() - tstart_begin
    write(INVERSE_LOG_FILE,*)
    write(INVERSE_LOG_FILE,*) ' Elapsed time for setup in seconds = ',tCPU
    write(INVERSE_LOG_FILE,*)
    call flush_iunit(INVERSE_LOG_FILE)
  endif

  ! initializes iterators
  iter_frq     = 0   ! current frequency stages
  iter_inverse = 0   ! current inversion iteration (per stage)
  iter_wolfe   = 0   ! current sub-inversion for Wolfe conditions (per iteration)

  !!!##############################################################################################################################
  !!! -------------------------------  different running mode : forward or FWI ----------------------------------------------------
  !!!##############################################################################################################################

  select case (trim(adjustl(mode_running)))

  case('forward')

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*************************************************************************'
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' Specfem inverse problem : Forward only problem ...  '
      write(INVERSE_LOG_FILE,*) ' with generic family parameters  : ', trim(inversion_param%parameter_family_name)
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! output_solver.txt
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'Inverse problem for model:'
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'run mode              : forward modeling'
      write(IMAIN,*) 'total number of events: ',acqui_simu(1)%nevent_tot
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! loops over events
    do ievent = 1, acqui_simu(1)%nevent_tot
      ! timing
      tstart = wtime()

      call ComputeSismosPerEvent(ievent, acqui_simu, inversion_param, myrank)

      ! elapsed time since beginning
      if (myrank == 0) then
        tCPU = wtime() - tstart
        write(INVERSE_LOG_FILE,*) ' done event: ',ievent
        write(INVERSE_LOG_FILE,*) ' Elapsed time for computing event = ',tCPU
        write(INVERSE_LOG_FILE,*)
        call flush_iunit(INVERSE_LOG_FILE)
      endif

    enddo

    call synchronize_all_world()

    !! writing model in SEM mesh : (rho, vp, vs) or cijkl.
    if (inversion_param%output_model)  call WriteOutputSEMmodel(inversion_param)

  case ('l-bfgs')

    if (myrank == 0) then
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*************************************************************************'
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' Specfem inverse problem : L-BFGS FWI ...  '
      write(INVERSE_LOG_FILE,*) ' with family parameters  : ', trim(inversion_param%parameter_family_name)
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' L-BFGS history          : ', inversion_param%max_history_bfgs
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)
    endif

    ! output_solver.txt
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'Inverse problem for model:'
      write(IMAIN,*) '**************************'
      write(IMAIN,*) 'run mode                        : l-bfgs'
      write(IMAIN,*) 'total number of events          : ',acqui_simu(1)%nevent_tot
      write(IMAIN,*) 'total number of frequency bands : ',inversion_param%Nifrq
      write(IMAIN,*) 'maximum number of iterations    : ',inversion_param%Niter
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! initialize specifics arrays for optimization
    call AllocatememoryForFWI(inversion_param, acqui_simu(1)%nevent_tot)

    ! loop on frequencies
    do iter_frq = 1, inversion_param%Nifrq

      ! user output
      if (myrank == 0) then
        write(*,*) ' Frequency stage       : ',iter_frq,' / ',inversion_param%Nifrq
        call flush_iunit(6)
      endif

      ! initialize (reset) optimization
      call InitializeOptimIteration(acqui_simu, inversion_param)

      ! loop on optimization iterations
      ! note: we start with iter_inverse still at zero, as current iteration counter.
      do iter_inverse = 0, inversion_param%Niter
        ! timing
        tstart = wtime()

        ! user output
        if (myrank == 0) then
          write(*,*) '   running iteration : ',iter_inverse
          call flush_iunit(6)
        endif

        ! one step for optimization scheme iterations
        ! we will use the initial gradient computed in routine InitializeOptimIteration() for a first line search trial step;
        ! once the Wolfe's criteria are satisfied, the model and gradients are updated, and the iteration counter can increase.
        call OneIterationOptim(finished, acqui_simu, inversion_param) !!!!!! , regularization_fd)

        ! elapsed time
        if (myrank == 0) then
          tCPU = wtime() - tstart
          int_tCPU = int(tCPU)
          ihours = int_tCPU / 3600
          iminutes = (int_tCPU - 3600*ihours) / 60
          iseconds = int_tCPU - 3600*ihours - 60*iminutes
          write(INVERSE_LOG_FILE,*) 'done frequency band : ',iter_frq,' iteration : ',iter_inverse
          write(INVERSE_LOG_FILE,*) 'Elapsed time for computing iteration in seconds  = ',tCPU
          write(INVERSE_LOG_FILE,"(' Elapsed time for computing iteration in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
            ihours,iminutes,iseconds
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
        endif

        if (finished) exit

      enddo

      ! write model in disk
      call WriteOutputs(inversion_param)

      ! user output
      if (myrank == 0) then
        write(*,*)
        call flush_iunit(6)
      endif

    enddo

    ! user output
    if (myrank == 0) then
      if (iter_inverse >= inversion_param%Niter) then
        write(INVERSE_LOG_FILE,*)
        write(INVERSE_LOG_FILE,*) ' FWI STOP : maximum allowed iteration reached ', inversion_param%Niter
        write(INVERSE_LOG_FILE,*)
        iter_inverse = iter_inverse - 1
      endif
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) ' FWI STOPPED at iteration :', iter_inverse, ' / ', inversion_param%Niter
      write(INVERSE_LOG_FILE,*)
      write(INVERSE_LOG_FILE,*) '*************************************************************************'
      write(INVERSE_LOG_FILE,*)
      call flush_iunit(INVERSE_LOG_FILE)

      write(*,*)
      write(*,*) ' FWI STOPPED at iteration :', iter_inverse, ' / ', inversion_param%Niter
      write(*,*)
      call flush_iunit(6)
    endif

  case default

    if (myrank == 0) write(*,*) ' ERROR :', trim(mode_running),  ':  option not defined '

  end select

  ! closes debug files
  close(IIDD)

  ! closes log files
  if (myrank == 0) then
    ! elapsed time since beginning
    tCPU = wtime() - tstart_begin
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    ! log file
    write(INVERSE_LOG_FILE,*)
    write(INVERSE_LOG_FILE,*) 'Timing info:'
    write(INVERSE_LOG_FILE,*) 'Total elapsed time in seconds  = ',tCPU
    write(INVERSE_LOG_FILE,"(' Total elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") &
                              ihours,iminutes,iseconds
    write(INVERSE_LOG_FILE,*)
    write(INVERSE_LOG_FILE,*) 'End of the inversion'
    write(INVERSE_LOG_FILE,*)
    close(INVERSE_LOG_FILE)

    ! iteration logs
    if (.not. inversion_param%only_forward) then
      close(OUTPUT_ITERATION_FILE)
      close(OUTPUT_FWI_LOG)
    endif
  endif

  ! parce que sinon je dois mettre un rpint...
  call synchronize_all_world()

end subroutine inverse_problem_main
