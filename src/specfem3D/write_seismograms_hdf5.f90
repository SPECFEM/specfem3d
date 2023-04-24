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

  endif ! myrank == 0
  deallocate(one_seismogram)
  deallocate(all_seismograms)
  deallocate(tmp_seis_pre)

  end subroutine write_seismo_no_ioserv

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file_h5(istore)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES,NB_RUNS_ACOUSTIC_GPU

  use specfem_par, only: myrank,number_receiver_global,nrec_local,it,DT, &
    WRITE_SEISMOGRAMS_BY_MAIN,SU_FORMAT,ASDF_FORMAT, &
    t0,seismo_offset,seismo_current, &
    NTSTEP_BETWEEN_OUTPUT_SAMPLE,nlength_seismogram

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p

  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: all_seismograms

  ! local parameters
  integer :: i,ier,irec,irec_local,nrec_store
  integer :: iorientation,isample,it_tmp,number_of_components
  real(kind=CUSTOM_REAL) :: time_t

  character(len=3) :: channel
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  ! checks
  if (NB_RUNS_ACOUSTIC_GPU > 1) &
    stop 'for NB_RUNS_ACOUSTIC_GPU > 1, adjoint wavefield seismograms output not implemented yet'
  if (WRITE_SEISMOGRAMS_BY_MAIN) &
    stop 'for WRITE_SEISMOGRAMS_BY_MAIN, adjoint wavefield seismograms output not implemented yet'
  if (ASDF_FORMAT) &
    stop 'for ASDF_FORMAT, adjoint wavefield seismograms output not implemented yet'

  ! saves displacement, velocity, acceleration, or pressurek
  if (istore == 1) then
    component = 'd'
  else if (istore == 2) then
    component = 'v'
  else if (istore == 3) then
    component = 'a'
  else if (istore == 4) then
    component = 'p'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  ! set number of traces in full seismogram array
  ! each process outputs its local traces
  nrec_store = nrec_local

  ! allocates array for all seismograms
  allocate(all_seismograms(NDIM,nlength_seismogram,nrec_store),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating all_seismograms array')
  if (ier /= 0) stop 'error while allocating all_seismograms array'
  all_seismograms(:,:,:) = 0._CUSTOM_REAL

  ! each process writes out its local receivers
  ! stores local seismograms
  do irec_local = 1,nrec_local
    ! sets trace
    select case (istore)
    case (1)
      ! displacement
      do i = 1,seismo_current
        all_seismograms(:,i,irec_local) = seismograms_d(:,irec_local,i)
      enddo
    case (2)
      ! velocity
      do i = 1,seismo_current
        all_seismograms(:,i,irec_local) = seismograms_v(:,irec_local,i)
      enddo
    case (3)
      ! acceleration
      do i = 1,seismo_current
        all_seismograms(:,i,irec_local) = seismograms_a(:,irec_local,i)
      enddo
    case (4)
      ! pressure
      do i = 1,seismo_current
        all_seismograms(1,i,irec_local) = seismograms_p(1,irec_local,i) ! single component
      enddo
    end select
  enddo

  ! writes file output
  if (SU_FORMAT) then
    ! SU_format file output
    ! write ONE binary file for all receivers (nrec_local) within one proc
    ! SU format, with 240-byte-header for each trace
    call write_output_SU(all_seismograms,nrec_store,istore)
  else
    ! ASCII ouput
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      ! see how many components we need to store: 1 for pressure, NDIM for a vector
      if (istore == 4) then ! this is for pressure
        number_of_components = 1
      else
        number_of_components = NDIM
      endif

      ! loop over each seismogram component
      do iorientation = 1,number_of_components

        ! gets channel name
        if (istore == 4) then
          ! this is for pressure
          call write_channel_name(istore,channel)
        else
          call write_channel_name(iorientation,channel)
        endif

        ! create the name of the seismogram file for each slice
        ! file name includes the name of the station, the network and the component
        write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,channel,component

        ! save seismograms in text format
        !
        ! note: subsampling the seismograms can result in a loss of accuracy when one later convolves
        !       the results with the source time function - know what you are doing.
        if (seismo_offset == 0) then
          ! open new file
          open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)), &
               status='unknown',action='write')
        else
          ! append to existing file
          open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)), &
               status='old',position='append',action='write')
        endif

        ! make sure we never write more than the maximum number of time steps
        do isample = 1,seismo_current
          ! current time
          ! current time increment
          it_tmp = seismo_offset + isample

          ! subtract half duration of the source to make sure travel time is correct
          time_t = real(dble((it_tmp-1)*NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0,kind=CUSTOM_REAL)

          ! distinguish between single and double precision for reals
          write(IOUT,*) time_t,all_seismograms(iorientation,isample,irec_local)
        enddo

        close(IOUT)

      enddo
    enddo
  endif ! SU_FORMAT

  deallocate(all_seismograms)



  end subroutine write_adj_seismograms_to_file_h5

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file_h5(myrank,seismograms_eps,number_receiver_global,nrec_local,it,DT,NSTEP,t0)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  implicit none
  integer :: myrank
  integer :: nrec_local,NSTEP,it
  integer, dimension(nrec_local) :: number_receiver_global
  ! note: seismograms here is still an array of size *,*,*,NSTEP
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local,NSTEP) :: seismograms_eps
  double precision :: t0,DT

  ! local parameters
  integer :: irec,irec_local
  integer :: idimval,jdimval,isample

  character(len=4) :: chn
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  component = 'd'

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do idimval = 1,NDIM
      do jdimval = idimval,NDIM

        ! strain channel name
        if (idimval == 1 .and. jdimval == 1) then
          chn = 'SNN'
        else if (idimval == 1 .and. jdimval == 2) then
          chn = 'SEN'
        else if (idimval == 1 .and. jdimval == 3) then
          chn = 'SEZ'
        else if (idimval == 2 .and. jdimval == 2) then
          chn = 'SEE'
        else if (idimval == 2 .and. jdimval == 3) then
          chn = 'SNZ'
        else if (idimval == 3 .and. jdimval == 3) then
          chn = 'SZZ'
        else
          call exit_MPI(myrank,'incorrect channel value')
        endif

        ! create the name of the seismogram file for each slice
        ! file name includes the name of the station, the network and the component
        write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,chn,component

        ! save seismograms in text format with no subsampling.
        ! Because we do not subsample the output, this can result in large files
        ! if the simulation uses many time steps. However, subsampling the output
        ! here would result in a loss of accuracy when one later convolves
        ! the results with the source time function
        open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)),status='unknown')

        ! make sure we never write more than the maximum number of time steps
        ! subtract half duration of the source to make sure travel time is correct
        do isample = 1,min(it,NSTEP)
          ! distinguish between single and double precision for reals
          write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ',seismograms_eps(jdimval,idimval,irec_local,isample)
        enddo

        close(IOUT)

      enddo ! jdimval
    enddo ! idimval
  enddo ! irec_local

  end subroutine write_adj_seismograms2_to_file_h5


!=====================================================================

! write strain seismograms to text files

  subroutine write_seismograms_strain_to_file_h5()

  use constants, only: myrank,CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  use specfem_par, only: number_receiver_global,nrec_local, &
    network_name,station_name,DT,NSTEP,t0, &
    seismo_offset,seismo_current,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    SUPPRESS_UTM_PROJECTION,WRITE_SEISMOGRAMS_BY_MAIN,SIMULATION_TYPE

  ! strain seismo
  use specfem_par, only: seismograms_eps

  implicit none

  ! local parameters
  integer :: irec,irec_local
  integer :: idimval,jdimval,isample,it_tmp,ier
  real(kind=CUSTOM_REAL) :: time_t

  character(len=4) :: chn
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

! strain
! for adjoint simulations, nrec_local is the number of "receivers" at which we store seismograms and the strain (epsilon)
! field. "receiver" locations for adjoint simulations are determined by the CMTSOLUTIONs locations, i.e., original sources.
! Thus, we flip this assignment for pure adjoint simulations, that is source locations becomes receiver locations,
! and receiver locations become "adjoint source" locations.

  ! safety check
  if (WRITE_SEISMOGRAMS_BY_MAIN) &
    call exit_MPI(myrank,'Error write_seismograms_strain_to_file() needs WRITE_SEISMOGRAMS_BY_MAIN turned off')

  ! checks if anything to do
  if (nrec_local <= 0 ) return

  ! using an ending .semd also for strain seismograms (since strain is derived from displacement)
  component = 'd'

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! reference frame convention:
    !   SPECFEM reference frame convention: x = East, y = North, z = Up
    !
    !   UTM convention: uses convention for N/E/Z
    !                   'longitude' parameters simply refer to the x axis
    !                   'latitude'  parameters simply refer to the y axis
    !
    ! SPECFEM3D Cartesian uses the traditional orientation code E/N/Z (East-West, North-South, Vertical)
    ! for three components when a UTM projection is used.
    !
    ! If the UTM conversion is suppressed, i.e. the flag SUPPRESS_UTM_PROJECTION is set to .true.,
    ! then the three components are labelled X/Y/Z according to the Cartesian reference frame.
    !
    !
    ! strain "seismogram":
    !   seismograms_eps(:,:) = eps_s(:,:)
    ! with
    !   symmetric strain tensor: eps_s(1,1) = eps_xx
    !                            eps_s(2,2) = eps_yy
    !                            eps_s(3,3) = eps_zz
    !                            eps_s(1,2) = eps_xy
    !                            eps_s(1,3) = eps_xz
    !                            eps_s(2,3) = eps_yz
    do idimval = 1,NDIM
      do jdimval = idimval,NDIM

        ! strain channel name
        if (SUPPRESS_UTM_PROJECTION) then
          ! no UTM
          ! uses Cartesian X/Y/Z direction to denote channel
          if (idimval == 1 .and. jdimval == 1) then
            ! eps(1,1) = eps_xx
            chn = 'SXX'
          else if (idimval == 1 .and. jdimval == 2) then
            ! eps(2,1) = eps_yx = eps_xy
            chn = 'SXY'
          else if (idimval == 1 .and. jdimval == 3) then
            ! eps(3,1) = eps_zx = eps_xz
            chn = 'SXZ'
          else if (idimval == 2 .and. jdimval == 2) then
            ! eps(2,2) = eps_yy
            chn = 'SYY'
          else if (idimval == 2 .and. jdimval == 3) then
            ! eps(3,2) = eps_zy = eps_yz
            chn = 'SYZ'
          else if (idimval == 3 .and. jdimval == 3) then
            ! eps(3,3) = eps_zz
            chn = 'SZZ'
          else
            call exit_MPI(myrank,'incorrect channel value')
          endif
        else
          ! UTM conversion
          ! uses convention for N/E/Z to denote channel
          if (idimval == 1 .and. jdimval == 1) then
            ! eps(1,1) = eps_xx
            chn = 'SEE'
          else if (idimval == 1 .and. jdimval == 2) then
            ! eps(2,1) = eps_yx = eps_xy
            chn = 'SEN'
          else if (idimval == 1 .and. jdimval == 3) then
            ! eps(3,1) = eps_zx = eps_xz
            chn = 'SEZ'
          else if (idimval == 2 .and. jdimval == 2) then
            ! eps(2,2) = eps_yy
            chn = 'SNN'
          else if (idimval == 2 .and. jdimval == 3) then
            ! eps(3,2) = eps_zy = eps_yz
            chn = 'SNZ'
          else if (idimval == 3 .and. jdimval == 3) then
            ! eps(3,3) = eps_zz
            chn = 'SZZ'
          else
            call exit_MPI(myrank,'incorrect channel value')
          endif
        endif

        ! name
        if (SIMULATION_TYPE == 2) then
          ! pure adjoint simulation
          ! create the name of the seismogram file for each source
          ! file name includes the name numbered by the source S0000*, channel and the component
          ! for example: /NT.S000001.SNN.semd
          write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,chn,component
        else
          ! strain seismograms will have channel-code S**
          !
          ! create the name of the strain seismogram file using the station name and network name
          ! using format: **net**.**sta**.channel
          ! for example: /IU.KONO.SNN.semd
          write(sisname,"('/',a,'.',a,'.',a3,'.sem',a1)") trim(network_name(irec)),trim(station_name(irec)),chn,component
        endif

        ! save seismograms in text format
        if (seismo_offset == 0) then
          ! open new file
          open(unit=IOUT,file=trim(OUTPUT_FILES)//sisname(1:len_trim(sisname)), &
               status='unknown',action='write',iostat=ier)
        else
          ! for it > NTSTEP_BETWEEN_OUTPUT_SEISMOS
          ! append to existing file
          open(unit=IOUT,file=trim(OUTPUT_FILES)//sisname(1:len_trim(sisname)), &
               status='old',position='append',action='write',iostat=ier)
        endif
        if (ier /= 0) call exit_mpi(myrank,'Error opening file: '//trim(OUTPUT_FILES)//trim(sisname))

        ! make sure we never write more than the maximum number of time steps
        do isample = 1,seismo_current
          ! current time
          it_tmp = seismo_offset + isample
          ! subtract onset time to make sure travel time is correct
          if (SIMULATION_TYPE == 3) then
            time_t = real(dble((NSTEP-it_tmp)*NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0,kind=CUSTOM_REAL)
          else
            time_t = real(dble((it_tmp-1)*NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0,kind=CUSTOM_REAL)
          endif

          ! output
          write(IOUT,*) time_t,seismograms_eps(jdimval,idimval,irec_local,isample)
        enddo

        close(IOUT)

      enddo ! jdimval
    enddo ! idimval
  enddo ! irec_local

  end subroutine write_seismograms_strain_to_file_h5

