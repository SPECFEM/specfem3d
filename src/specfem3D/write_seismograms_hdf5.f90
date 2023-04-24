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

! write seismograms to hdf5 files

  subroutine write_seismograms_h5()

! writes the seismograms with time shift

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  ! timing
  double precision, external :: wtime
  double precision :: write_time_begin,write_time
  logical :: do_save_seismograms

  ! note: there might be some confusion about adjoint simulations, i.e. adjoint receivers and adjoint sources:
  !       for pure adjoint simulations (SIMULATION_TYPE == 2), CMT source locations become adjoint receivers
  !       for recording seismograms; station locations become (possible) adjoint sources.
  !       thus,
  !       1. adjoint sources are located at the receiver positions given in STATIONS file,
  !          'nadj_rec_local' is the number of local adjoint sources, i.e. station positions acting as adjoint source
  !       2. adjoint "receivers" are located at the CMT source positions given in CMTSOLUTION file,
  !          'nrec_local' is the number of local adjoint receivers, i.e. source positions acting as receivers for
  !          recording adjoint seismograms.

  ! remember for pure adjoint runs (see setup_receivers routine):
  ! - nrec_local              -> between 1 to NSOURCES
  ! - number_receiver_global  -> local to global mapping for "adjoint" receivers, i.e. at source positions
  ! - ispec_selected_source   -> element containing source position, which becomes an adjoint "receiver"

  ! check if we need to save seismos
  if (SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) then
    do_save_seismograms = .false.
  else
    do_save_seismograms = .true.
  endif

  ! checks subsampling recurrence
  if (mod(it-1,NTSTEP_BETWEEN_OUTPUT_SAMPLE) == 0) then

    ! update position in seismograms
    seismo_current = seismo_current + 1

    ! check for edge effects
    if (seismo_current < 1 .or. seismo_current > nlength_seismogram) &
      call exit_mpi(myrank,'Error: seismo_current out of bounds in recording of seismograms')

    ! note: there might be some confusion about adjoint simulations, i.e. adjoint receivers and adjoint sources:
    !       for pure adjoint simulations (SIMULATION_TYPE == 2), CMT source locations become adjoint receivers
    !       for recording seismograms; station locations become (possible) adjoint sources.
    !       thus,
    !       1. adjoint sources are located at the receiver positions given in STATIONS file,
    !          'nadj_rec_local' is the number of local adjoint sources, i.e. station positions acting as adjoint source
    !       2. adjoint "receivers" are located at the CMT source positions given in CMTSOLUTION file,
    !          'nrec_local' is the number of local adjoint receivers, i.e. source positions acting as receivers for
    !          recording adjoint seismograms.

    ! remember for pure adjoint runs (see setup_receivers routine):
    ! - nrec_local              -> between 1 to NSOURCES
    ! - number_receiver_global  -> local to global mapping for "adjoint" receivers, i.e. at source positions
    ! - ispec_selected_source   -> element containing source position, which becomes an adjoint "receiver"

    ! gets seismogram values
    if (do_save_seismograms .and. nrec_local > 0) then
      ! on GPU
      if (GPU_MODE) then
        ! gets resulting array values onto CPU
        call compute_seismograms_cuda(Mesh_pointer,seismograms_d,seismograms_v,seismograms_a,seismograms_p, &
                                      seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS,it,it_end, &
                                      ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,USE_TRICK_FOR_BETTER_PRESSURE)
      else
        ! on CPU
        call compute_seismograms()
      endif

        ! strain seismograms
        if (SAVE_SEISMOGRAMS_STRAIN) call compute_seismograms_strain()

        ! computes source moment derivatives (only for adjoint simulations)
        if (SIMULATION_TYPE == 2) call compute_seismograms_moment_adjoint()
    endif

  endif ! NTSTEP_BETWEEN_OUTPUT_SAMPLE

  ! write the current or final seismograms
  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == it_end) then

    if (do_save_seismograms .and. .not. INVERSE_FWI_FULL_PROBLEM) then

      ! timing
      write_time_begin = wtime()

      ! checks if anything to do
      ! note: ASDF uses parallel hdf5 that defines the MPI communicator group that the solver is
      !       run with. this means every processor in the group is needed for write_seismograms
      if (nrec_local > 0 .or. (WRITE_SEISMOGRAMS_BY_MAIN .and. myrank == 0) .or. ASDF_FORMAT) then

        ! writes out seismogram files
        select case(SIMULATION_TYPE)
        case (1,3)
          ! forward/kernel simulations
          ! forward/backward wavefield
          ! forward & kernel simulations
          if (NIONOD > 0) then
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_seismograms_to_file_h5(1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_seismograms_to_file_h5(2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_seismograms_to_file_h5(3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_seismograms_to_file_h5(4)
          else
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_seismo_no_ioserv(1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_seismo_no_ioserv(2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_seismo_no_ioserv(3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_seismo_no_ioserv(4)
          endif
        case (2)
          ! adjoint simulations
          ! adjoint wavefield
          ! adjoint simulations
          !if (NIONOD > 0) then
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_adj_seismograms_to_file_h5(1) ! TODO: not implemented yet, currently dummy which out non-hdf5 format
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_adj_seismograms_to_file_h5(2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_adj_seismograms_to_file_h5(3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_adj_seismograms_to_file_h5(4)
          !else
          !  if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
          !    call write_seismo_no_ioserv(seismograms_d,1)
          !  if (SAVE_SEISMOGRAMS_VELOCITY) &
          !    call write_seismo_no_ioserv(seismograms_v,2)
          !  if (SAVE_SEISMOGRAMS_ACCELERATION) &
          !    call write_seismo_no_ioserv(seismograms_a,3)
          !  if (SAVE_SEISMOGRAMS_PRESSURE) &
          !    call write_seismo_no_ioserv(seismograms_p,4)
          !endif
        end select

        ! strain
        if (SAVE_SEISMOGRAMS_STRAIN) call write_seismograms_strain_to_file_h5() ! TODO: not implemented yet
      endif

      ! synchronizes processes (waits for all processes to finish writing)
      call synchronize_all()

      ! user output
      if (myrank == 0) then
        ! timing
        write_time = wtime() - write_time_begin
        ! output
        write(IMAIN,*)
        write(IMAIN,*) 'Total number of time steps done: ', it-it_begin+1
        if (WRITE_SEISMOGRAMS_BY_MAIN) then
          write(IMAIN,*) 'Writing the seismograms by main proc alone took ',sngl(write_time),' seconds'
        else
          write(IMAIN,*) 'Writing the seismograms in parallel took ',sngl(write_time),' seconds'
        endif
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    endif ! do_save_seismograms

    ! resets current seismogram position
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0

  endif


  end subroutine write_seismograms_h5



  subroutine write_seismograms_to_file_h5(istore)

  use constants

  use specfem_par, only: myrank,number_receiver_global,NPROC, &
          nrec,nrec_local,islice_selected_rec, &
          seismo_offset,seismo_current, &
          nlength_seismogram, NTSTEP_BETWEEN_OUTPUT_SEISMOS

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p

  implicit none

  integer :: istore, i

  ! local parameters
  character(len=4) component

  integer :: nrec_local_received,total_seismos,req
  integer :: io_tag_seismo, time_window

  real(kind=CUSTOM_REAL), dimension(NDIM,nlength_seismogram,nrec_local) :: seismograms
  seismograms(:,:,:) = 0.0

  ! saves displacement, velocity, acceleration, or pressure
  select case (istore)
  case (1)
    component = 'disp'
    io_tag_seismo = io_tag_seismo_body_disp
  case (2)
    component = 'velo'
    io_tag_seismo = io_tag_seismo_body_velo
  case (3)
    component = 'acce'
    io_tag_seismo = io_tag_seismo_body_acce
  case (4)
    component = 'pres'
    io_tag_seismo = io_tag_seismo_body_pres
  end select

  if (nrec_local > 0) then
    ! send seismograms (instead of one seismo)
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > nlength_seismogram) then
      time_window = nlength_seismogram
    else
      time_window = NTSTEP_BETWEEN_OUTPUT_SEISMOS
    endif

    select case (istore)
    case (1)
      do i = 1, time_window !
        seismograms(:,i,:) = seismograms_d(:,:,i)
      end do
      call isend_cr_inter(seismograms,NDIM*nrec_local*time_window,0,io_tag_seismo,req)
    case (2)
      do i = 1, time_window !
        seismograms(:,i,:) = seismograms_v(:,:,i)
      end do
      call isend_cr_inter(seismograms,NDIM*nrec_local*time_window,0,io_tag_seismo,req)
    case (3)
      do i = 1, time_window !
        seismograms(:,i,:) = seismograms_a(:,:,i)
      end do
      call isend_cr_inter(seismograms,NDIM*nrec_local*time_window,0,io_tag_seismo,req)
    case (4)
       do i = 1, time_window !
        seismograms(1,i,:) = seismograms_p(1,:,i)
      end do
      call isend_cr_inter(seismograms(1,:,:),1*nrec_local*time_window,0,io_tag_seismo,req)
    end select

  endif

  end subroutine write_seismograms_to_file_h5


  subroutine write_seismo_no_ioserv(istore)

  use specfem_par
  use phdf5_utils
  use io_server

  use constants

  use specfem_par, only: myrank,number_receiver_global,NPROC, &
          nrec,nrec_local,islice_selected_rec, &
          seismo_offset,seismo_current, &
          nlength_seismogram

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p


  implicit none

  character(len=4) :: component
  integer :: i,irec,irec_local,nrec_store,t_upper

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: all_seismograms
  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable :: tmp_seis_pre

  ! parameters for master collects seismograms
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
  integer :: nrec_local_received,total_seismos,receiver,sender
  integer :: iproc,ier
  integer,dimension(1) :: tmp_nrec_local_received,tmp_irec,tmp_nrec_local
  integer,dimension(:),allocatable:: islice_num_rec_local

  ! hdf5 vals
  character(len=64) :: fname_h5_base = "seismograms.h5"
  type(h5io) :: h5

  h5 = h5io()
  fname_h5_seismo = trim(OUTPUT_FILES)//fname_h5_base

  ! initialze hdf5
  call h5_init(h5, fname_h5_seismo)

  if ((seismo_h5_initialized .eqv. .false.) .and. (myrank == 0)) then
    call do_io_seismogram_init()
    seismo_h5_initialized = .true.
  endif

  ! checks: NB_RUNS_ACOUSTIC_GPU only works for pressure output, where we define the array
  !         seismograms_p(NDIM,nrec_store*NB_RUNS_ACOUSTIC_GPU,NSTEP_**)
  if (NB_RUNS_ACOUSTIC_GPU > 1 .and. istore /= 4) then
    print *,'Error: NB_RUNS_ACOUSTIC_GPU > 1 has invalid seismogram output type ',istore,' - component '//component//' choosen'
    print *,'Please output pressure only, other components are not implemented yet...'
    stop 'Invalid seismogram output for NB_RUNS_ACOUSTIC_GPU > 1'
  endif

  ! allocates single trace array
  allocate(one_seismogram(NDIM,nlength_seismogram),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2420')
  if (ier /= 0) stop 'error while allocating one temporary seismogram'
  one_seismogram(:,:) = 0._CUSTOM_REAL

  if (myrank == 0) then
    ! on the master, gather all the seismograms

    ! counts number of local receivers for each slice
    allocate(islice_num_rec_local(0:NPROC-1),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error allocating islice_num_rec_local')
    islice_num_rec_local(:) = 0
    do irec = 1,nrec
      iproc = islice_selected_rec(irec)
      islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
    enddo

    total_seismos = 0
    ! main process collects all traces
    if (myrank == 0) then
      ! explicit for pressure/acoustic runs and NB_RUNS_ACOUSTIC_GPU > 1
      nrec_store = nrec * NB_RUNS_ACOUSTIC_GPU
    else
      ! dummy for secondary processes
      nrec_store = 0
    endif

    allocate(all_seismograms(NDIM,nlength_seismogram,nrec_store),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating all_seismograms array')
    if (ier /= 0) stop 'error while allocating all_seismograms array'
    all_seismograms(:,:,:) = 0._CUSTOM_REAL

    allocate(tmp_seis_pre(nlength_seismogram,nrec_store),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating tmp_seis_pre array')
    if (ier /= 0) stop 'error while allocating tp_seis_pre array'
    tmp_seis_pre(:,:) = 0._CUSTOM_REAL

    ! loop on all the slices
    do iproc = 0,NPROC-1
      ! communicate only with processes which contain local receivers
      if (islice_num_rec_local(iproc) == 0) cycle

      ! receive except from proc 0, which is me and therefore I already have this value
      sender = iproc
      if (iproc /= 0) then
        call recv_i(tmp_nrec_local_received,1,sender,itag)
        nrec_local_received = tmp_nrec_local_received(1)
        if (nrec_local_received < 0) call exit_MPI(myrank,'error while receiving local number of receivers')
      else
        nrec_local_received = nrec_local * NB_RUNS_ACOUSTIC_GPU
      endif

      if (nrec_local_received > 0) then
        do irec_local = 1,nrec_local_received
          ! init trace
          one_seismogram(:,:) = 0._CUSTOM_REAL

          ! receive except from proc 0, which is myself and therefore I already have these values
          if (iproc == 0) then
            ! get global number of that receiver
            if (NB_RUNS_ACOUSTIC_GPU == 1) then
              irec = number_receiver_global(irec_local)
            else
              ! NB_RUNS_ACOUSTIC_GPU > 1
              ! if irec_local is a multiple of nrec_local then mod(irec_local,nrec_local) == 0 and
              ! the array access at number_receiver_global would be invalid;
              ! for those cases we want irec associated to irec_local == nrec_local
              if (mod(irec_local,nrec_local) == 0) then
                irec = number_receiver_global(nrec_local)
              else
                irec = number_receiver_global(mod(irec_local,nrec_local))
              endif
            endif

            ! gets trace
            select case (istore)
            case (1)
              ! displacement
              do i = 1,seismo_current
                one_seismogram(:,i) = seismograms_d(:,irec_local,i)
              enddo
            case (2)
              ! velocity
              do i = 1,seismo_current
                one_seismogram(:,i) = seismograms_v(:,irec_local,i)
              enddo
            case (3)
              ! acceleration
              do i = 1,seismo_current
                one_seismogram(:,i) = seismograms_a(:,irec_local,i)
              enddo
            case (4)
              ! pressure
              do i = 1,seismo_current
                one_seismogram(1,i) = seismograms_p(1,irec_local,i) ! single component
              enddo
            end select

          else
            call recv_i(tmp_irec,1,sender,itag)
            irec = tmp_irec(1)
            ! receives trace
            call recvv_cr(one_seismogram,NDIM*seismo_current,sender,itag)
          endif
          if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')

          total_seismos = total_seismos + 1

          ! stores this seismogram
          if (NB_RUNS_ACOUSTIC_GPU == 1) then
            all_seismograms(:,:,irec) = one_seismogram(:,:)
          else
            ! NB_RUNS_ACOUSTIC_GPU > 1
            ! fills all_seismgrams(..) array with all first run traces, then second run traces etc.
            if (irec_local <= islice_num_rec_local(sender)) then
              all_seismograms(:,:,irec) = one_seismogram(:,:)
            else
              ! irec_local > nrec_local in sender process
              irec = irec + int((irec_local-1)/islice_num_rec_local(sender)) * nrec
              all_seismograms(:,:,irec) = one_seismogram(:,:)
            endif
          endif
        enddo ! nrec_local_received
      endif ! if (nrec_local_received > 0)
    enddo ! NPROC-1

  else
    ! on the nodes, send the seismograms to the main
    receiver = 0
    tmp_nrec_local(1) = nrec_local * NB_RUNS_ACOUSTIC_GPU
    call send_i(tmp_nrec_local,1,receiver,itag)

    if (nrec_local > 0) then
      do irec_local = 1,nrec_local * NB_RUNS_ACOUSTIC_GPU
        ! get global number of that receiver
        if (NB_RUNS_ACOUSTIC_GPU == 1) then
          irec = number_receiver_global(irec_local)
        else
          ! NB_RUNS_ACOUSTIC_GPU > 1
          ! if irec_local is a multiple of nrec then mod(irec_local,nrec) == 0 and
          ! the array access at number_receiver_global would be invalid;
          ! for those cases we want irec associated to irec_local == nrec_store
          if (mod(irec_local,nrec_local) == 0) then
            irec = number_receiver_global(nrec_local)
          else
            irec = number_receiver_global(mod(irec_local,nrec_local))
          endif
        endif
        tmp_irec(1) = irec
        call send_i(tmp_irec,1,receiver,itag)

        ! sets trace
        one_seismogram(:,:) = 0._CUSTOM_REAL
        select case (istore)
          case (1)
            ! displacement
            do i = 1,seismo_current
              one_seismogram(:,i) = seismograms_d(:,irec_local,i)
            enddo
          case (2)
            ! velocity
            do i = 1,seismo_current
              one_seismogram(:,i) = seismograms_v(:,irec_local,i)
            enddo
          case (3)
            ! acceleration
            do i = 1,seismo_current
              one_seismogram(:,i) = seismograms_a(:,irec_local,i)
            enddo
          case (4)
            ! pressure
            do i = 1,seismo_current
              one_seismogram(1,i) = seismograms_p(1,irec_local,i) ! single component
            enddo
        end select

        ! sends seismogram of that receiver
        call sendv_cr(one_seismogram,NDIM*NSTEP,receiver,itag)
      enddo
    endif
  endif ! myrank

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  total number of receivers collected is ',total_seismos,' out of ',nrec * NB_RUNS_ACOUSTIC_GPU
    write(IMAIN,*)
    ! checks total
    if (total_seismos /= nrec * NB_RUNS_ACOUSTIC_GPU) call exit_MPI(myrank,'incorrect total number of receivers saved')
  endif

  if (myrank == 0) then
    ! check if the array length to be written > total timestep
    if (seismo_offset < 0) seismo_offset = 0
    if (seismo_offset+NTSTEP_BETWEEN_OUTPUT_SEISMOS > nlength_seismogram) then
      t_upper = nlength_seismogram - seismo_offset
    else
      t_upper = NTSTEP_BETWEEN_OUTPUT_SEISMOS
    endif

    ! writes out this seismogram
    call h5_open_file(h5)

    !call write_one_seismogram(one_seismogram,irec,component,istore)
    if (istore == 1) then
      component = 'disp'
      call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
        all_seismograms(:,1:t_upper,:), (/0, seismo_offset, 0/), .false.)

    else if (istore == 2) then
      component = 'velo'
      call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
           all_seismograms(:,1:t_upper,:), (/0, seismo_offset, 0/), .false.)
    else if (istore == 3) then
        component = 'acce'
      call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
           all_seismograms(:,1:t_upper,:), (/0, seismo_offset, 0/), .false.)
    else if (istore == 4) then
      component = 'pres'
      tmp_seis_pre = all_seismograms(1,1:t_upper,:)
      call h5_write_dataset_2d_r_collect_hyperslab(h5, component, &
           tmp_seis_pre, (/seismo_offset, 0/), .false.)
    else
      call exit_MPI(myrank,'wrong component to save for seismograms')
    endif

    call h5_close_file(h5)

    deallocate(all_seismograms)
    deallocate(tmp_seis_pre)

  endif ! myrank == 0

  deallocate(one_seismogram)

  end subroutine write_seismo_no_ioserv

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file_h5(istore)

  ! stop not implemented yet
  call exit_MPI(myrank,'write_adj_seismograms_to_file_h5 not implemented yet')

  end subroutine write_adj_seismograms_to_file_h5

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file_h5(myrank,seismograms_eps,number_receiver_global,nrec_local,it,DT,NSTEP,t0)

  ! stop not implemented yet
  call exit_MPI(myrank,'write_adj_seismograms2_to_file_h5 not implemented yet')

  end subroutine write_adj_seismograms2_to_file_h5


!=====================================================================

! write strain seismograms to text files

  subroutine write_seismograms_strain_to_file_h5()

  ! stop not implemented yet
  call exit_MPI(myrank,'write_seismograms_strain_to_file_h5 not implemented yet')

  end subroutine write_seismograms_strain_to_file_h5

