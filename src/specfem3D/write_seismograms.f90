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

  subroutine write_seismograms()

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

  ! check if we need to save seismos
  if (SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) then
    do_save_seismograms = .false.
  else
    do_save_seismograms = .true.
  endif

  ! checks subsampling recurrence
  if (mod(it-1,subsamp_seismos) == 0) then

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
      ! seismograms displ/veloc/accel/pressure
      if (GPU_MODE) then
        ! on GPU
        ! gets resulting array values onto CPU
        call compute_seismograms_cuda(Mesh_pointer,seismograms_d,seismograms_v,seismograms_a,seismograms_p, &
                                      seismo_current,nlength_seismogram,it,it_end, &
                                      ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,USE_TRICK_FOR_BETTER_PRESSURE)
      else
        ! on CPU
        call compute_seismograms()
      endif

      ! strain "seismograms" (only for adjoint simulations)
      if (SIMULATION_TYPE == 2) call compute_seismograms_strain_adjoint()
    endif

  endif ! subsamp_seismos

  ! write the current or final seismograms
  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == it_end) then

    if (do_save_seismograms .and. .not. INVERSE_FWI_FULL_PROBLEM) then

      ! timing
      write_time_begin = wtime()

      ! checks if anything to do
      if (nrec_local > 0 .or. (WRITE_SEISMOGRAMS_BY_MAIN .and. myrank == 0) .or. ASDF_FORMAT) then

        ! writes out seismogram files
        select case(SIMULATION_TYPE)
        case (1,3)
          ! forward/kernel simulations
          ! forward/backward wavefield
          if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
            call write_seismograms_to_file(1)
          if (SAVE_SEISMOGRAMS_VELOCITY) &
            call write_seismograms_to_file(2)
          if (SAVE_SEISMOGRAMS_ACCELERATION) &
            call write_seismograms_to_file(3)
          if (SAVE_SEISMOGRAMS_PRESSURE) &
            call write_seismograms_to_file(4)
        case (2)
          ! adjoint simulations
          ! adjoint wavefield
          if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
            call write_adj_seismograms_to_file(1)
          if (SAVE_SEISMOGRAMS_VELOCITY) &
            call write_adj_seismograms_to_file(2)
          if (SAVE_SEISMOGRAMS_ACCELERATION) &
            call write_adj_seismograms_to_file(3)
          if (SAVE_SEISMOGRAMS_PRESSURE) &
            call write_adj_seismograms_to_file(4)
          ! seismograms (strain)
          if (it == it_end) call write_adj_seismograms2_to_file()
          ! updates adjoint time counter
          it_adj_written = it
        end select
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

  end subroutine write_seismograms


!================================================================

! write seismograms to text files

  subroutine write_seismograms_to_file(istore)

  use constants

  use specfem_par, only: myrank,number_receiver_global,NPROC, &
          nrec,nrec_local,islice_selected_rec, &
          seismo_offset,seismo_current, &
          ASDF_FORMAT,SU_FORMAT, &
          WRITE_SEISMOGRAMS_BY_MAIN,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_SEISMOGRAMS, &
          nlength_seismogram

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p

  implicit none

  integer,intent(in) :: istore

  ! local parameters
  integer :: irec,irec_local,nrec_store

  ! parameters for main collects seismograms
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: all_seismograms

  integer :: nrec_local_received,total_seismos,receiver,sender
  integer :: iproc,ier,i
  integer,dimension(1) :: tmp_nrec_local_received,tmp_irec,tmp_nrec_local
  integer,dimension(0:NPROC-1) :: islice_num_rec_local

  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  ! saves displacement, velocity, acceleration, or pressure
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

  ! checks: NB_RUNS_ACOUSTIC_GPU only works for pressure output, where we define the array
  !         seismograms_p(NDIM,nrec_store*NB_RUNS_ACOUSTIC_GPU,NSTEP_**)
  if (NB_RUNS_ACOUSTIC_GPU > 1 .and. istore /= 4) then
    print *,'Error: NB_RUNS_ACOUSTIC_GPU > 1 has invalid seismogram output type ',istore,' - component '//component//' choosen'
    print *,'Please output pressure only, other components are not implemented yet...'
    stop 'Invalid seismogram output for NB_RUNS_ACOUSTIC_GPU > 1'
  endif
  if (SU_FORMAT .and. ASDF_FORMAT) &
    stop 'Please choose either SU_FORMAT or ASDF_FORMAT, both outputs together are not implemented yet...'

  ! set number of traces in full seismogram array
  if (WRITE_SEISMOGRAMS_BY_MAIN) then
    ! main process collects all traces
    if (myrank == 0) then
      ! explicit for pressure/acoustic runs and NB_RUNS_ACOUSTIC_GPU > 1
      nrec_store = nrec * NB_RUNS_ACOUSTIC_GPU
    else
      ! dummy for secondary processes
      nrec_store = 0
    endif
  else
    ! each process outputs its local traces
    nrec_store = nrec_local * NB_RUNS_ACOUSTIC_GPU
  endif

  ! allocates single trace array
  allocate(one_seismogram(NDIM,nlength_seismogram),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2420')
  if (ier /= 0) stop 'error while allocating one temporary seismogram'
  one_seismogram(:,:) = 0._CUSTOM_REAL

  ! allocates array for all seismograms
  allocate(all_seismograms(NDIM,nlength_seismogram,nrec_store),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating all_seismograms array')
  if (ier /= 0) stop 'error while allocating all_seismograms array'
  all_seismograms(:,:,:) = 0._CUSTOM_REAL

  ! main process writes out all
  if (WRITE_SEISMOGRAMS_BY_MAIN) then
    ! only the main process does the writing of seismograms
    !
    ! collects seismograms from other processes
    if (myrank == 0) then
      ! on the main, gather all the seismograms
      ! counts number of local receivers for each slice
      islice_num_rec_local(:) = 0
      do irec = 1,nrec
        iproc = islice_selected_rec(irec)
        islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
      enddo
      total_seismos = 0

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
            if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while collecting global receiver number')

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
          call sendv_cr(one_seismogram,NDIM*seismo_current,receiver,itag)
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

  else
    ! each process writes out its local receivers
    ! stores local seismograms
    do irec_local = 1,nrec_local * NB_RUNS_ACOUSTIC_GPU
      ! init trace
      one_seismogram(:,:) = 0._CUSTOM_REAL

      ! sets trace
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

      ! stores this seismogram
      all_seismograms(:,:,irec_local) = one_seismogram(:,:)
    enddo
  endif ! WRITE_SEISMOGRAMS_BY_MAIN

  ! writes out seismograms
  if (.not. WRITE_SEISMOGRAMS_BY_MAIN) then
    ! all processes write their local seismograms themselves
    !
    ! file output
    if (SU_FORMAT) then
      ! SU_format
      ! write ONE binary file for all receivers (nrec_local) within one proc
      ! SU format, with 240-byte-header for each trace
      call write_output_SU(all_seismograms,nrec_store,istore)
    else
      ! ASCII output
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        ! note: saving in a single file while using multiple outputs when NTSTEP_BETWEEN_OUTPUT_SEISMOS is less than NSTEP
        !       will shuffle the entries and give unordered results. thus for now, we will only allow this when writing out
        !       the full arrays.
        write(sisname,'(A,I5.5)') '/all_seismograms_'//component//'_node_',myrank
        if (USE_BINARY_FOR_SEISMOGRAMS) then
          if (seismo_offset == 0) then
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown', &
                 form='unformatted',action='write')
          else
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old',position='append', &
                 form='unformatted',action='write')
          endif
        else
          if (seismo_offset == 0) then
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown', &
                 form='formatted',action='write')
          else
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old',position='append', &
                 form='formatted',action='write')
          endif
        endif
      endif

      ! writes one file per single trace
      ! loop on all the local receivers
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

        ! writes out this seismogram
        one_seismogram(:,:) = all_seismograms(:,:,irec_local)
        call write_one_seismogram(one_seismogram,irec_local,irec,component,istore)
      enddo ! nrec_local

      ! writes out ASDF container to the file
      if (ASDF_FORMAT) then
        call write_asdf()
        ! deallocate the container
        call close_asdf_data()
      endif

      ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)
    endif

  else
    ! main process writes all seismograms
    !
    ! only written out by main process - WRITE_SEISMOGRAMS_BY_MAIN
    if (myrank == 0) then
      ! file output
      if (SU_FORMAT) then
        ! SU_format
        ! write ONE binary file for all receivers (nrec) within this main proc
        ! SU format, with 240-byte-header for each trace
        call write_output_SU(all_seismograms,nrec_store,istore)
      else
        ! ASCII output
        ! create one large file instead of one small file per station to avoid file system overload
        if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
          write(sisname,'(A)') '/all_seismograms_'//component//'_main'
          if (USE_BINARY_FOR_SEISMOGRAMS) then
            if (seismo_offset == 0) then
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown', &
                   form='unformatted',action='write')
            else
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old',position='append', &
                   form='unformatted',action='write')
            endif
          else
            if (seismo_offset == 0) then
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown', &
                   form='formatted',action='write')
            else
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old',position='append', &
                   form='formatted',action='write')
            endif
          endif
        endif

        ! writes one file per single trace
        ! loop on all the receivers
        do irec_local = 1,nrec_store
          ! get global number of that receiver
          if (NB_RUNS_ACOUSTIC_GPU == 1) then
            irec = irec_local
          else
            ! NB_RUNS_ACOUSTIC_GPU > 1
            if (mod(irec_local,nrec) == 0) then
              irec = nrec
            else
              irec = mod(irec_local,nrec)
            endif
          endif

          ! single seismogram
          one_seismogram(:,:) = all_seismograms(:,:,irec_local)

          ! writes out this seismogram
          call write_one_seismogram(one_seismogram,irec_local,irec,component,istore)
        enddo
      endif ! SU_FORMAT

      write(IMAIN,*) 'Component: .sem'//component
      write(IMAIN,*) '  total number of receivers saved is ',nrec
      write(IMAIN,*)

      ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    endif ! myrank

    ! ASDF output
    if (ASDF_FORMAT) then
      ! ASDF format
      ! writes out ASDF container to the file
      call write_asdf()
      call synchronize_all()
      call close_asdf_data()
    endif

  endif ! WRITE_SEISMOGRAMS_BY_MAIN

  deallocate(one_seismogram)
  deallocate(all_seismograms)

  end subroutine write_seismograms_to_file

!=====================================================================

  subroutine write_one_seismogram(one_seismogram,irec_local,irec,component,istore)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,OUTPUT_FILES,NB_RUNS_ACOUSTIC_GPU

  use specfem_par, only: DT,t0,it, &
    NSTEP,ASDF_FORMAT,WRITE_SEISMOGRAMS_BY_MAIN, &
    SIMULATION_TYPE,station_name,network_name, &
    nrec,nrec_local,nlength_seismogram

  implicit none

  integer, intent(in) :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nlength_seismogram),intent(in) :: one_seismogram
  integer,intent(in) :: irec_local,irec
  character(len=1),intent(in) :: component

  ! local parameters
  integer :: iorientation,number_of_components
  integer :: length_station_name,length_network_name
  character(len=MAX_STRING_LEN) :: sisname,final_LOCAL_PATH,runname
  character(len=3) :: channel

  ! see how many components we need to store: 1 for pressure, NDIM for a vector
  if (istore == 4) then ! this is for pressure
    number_of_components = 1
  else
    number_of_components = NDIM
  endif

  ! loop over each seismogram component
  do iorientation = 1,number_of_components

    ! gets channel name
    if (istore == 4) then ! this is for pressure
      call write_channel_name(istore,channel)
    else
      call write_channel_name(iorientation,channel)
    endif

    ! create the name of the seismogram file for each slice
    ! file name includes the name of the station, the network and the component
    length_station_name = len_trim(station_name(irec))
    length_network_name = len_trim(network_name(irec))

    ! writes out **net**.**sta**.**BH**.sem* files
    write(sisname,"(a,'.',a,'.',a3,'.sem',a1)") network_name(irec)(1:length_network_name), &
                                                station_name(irec)(1:length_station_name),channel,component

    ! adds a run designator to trace name
    if (NB_RUNS_ACOUSTIC_GPU > 1) then
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        ! main process writes out all
        write(runname,"('NB_run',i1,'.')") int((irec_local-1)/nrec) + 1
        sisname = trim(runname) // trim(sisname)
      else
        ! each process writes out local receivers
        write(runname,"('NB_run',i1,'.')") int((irec_local-1)/nrec_local) + 1
        sisname = trim(runname) // trim(sisname)
      endif
    endif

    ! directory to store seismograms
    final_LOCAL_PATH = OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // '/'

    ! write/store seismogram
    if (ASDF_FORMAT) then
      ! ASDF storage
      call store_asdf_data(one_seismogram,irec_local,irec,channel,iorientation)
    else
      ! ASCII output format
      call write_output_ASCII_or_binary(one_seismogram, &
                                        NSTEP,it,SIMULATION_TYPE,DT,t0, &
                                        iorientation,sisname,final_LOCAL_PATH)
    endif

  enddo ! do iorientation

  end subroutine write_one_seismogram

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file(istore)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES,NB_RUNS_ACOUSTIC_GPU

  use specfem_par, only: myrank,number_receiver_global,nrec_local,DT, &
    WRITE_SEISMOGRAMS_BY_MAIN,SU_FORMAT,ASDF_FORMAT, &
    t0,seismo_offset,seismo_current, &
    it_adj_written,subsamp_seismos,nlength_seismogram

  ! seismograms
  use specfem_par, only: seismograms_d,seismograms_v,seismograms_a,seismograms_p

  implicit none

  integer,intent(in) :: istore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: all_seismograms

  ! local parameters
  integer :: i,ier,irec,irec_local,nrec_store
  integer :: iorientation,isample,number_of_components
  double precision :: time

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

  ! saves displacement, velocity, acceleration, or pressure
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
          !open new file
          open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)), &
               status='unknown',action='write')
        else
          !append to existing file
          open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)), &
               status='old',position='append',action='write')
        endif

        ! make sure we never write more than the maximum number of time steps
        do isample = 1,seismo_current
          ! current time
          ! subtract half duration of the source to make sure travel time is correct
          time = dble(isample-1)*DT*subsamp_seismos + dble(it_adj_written)*DT - t0

          ! distinguish between single and double precision for reals
          write(IOUT,*) real(time,kind=CUSTOM_REAL),' ',all_seismograms(iorientation,isample,irec_local)
        enddo

        close(IOUT)

      enddo
    enddo
  endif ! SU_FORMAT

  deallocate(all_seismograms)

  end subroutine write_adj_seismograms_to_file

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file()

  use constants, only: myrank,IMAIN,CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  use specfem_par, only: seismograms_eps, &
    number_receiver_global,nrec_local,DT,NSTEP,t0, &
    seismo_offset,seismo_current,subsamp_seismos, &
    SUPPRESS_UTM_PROJECTION

  implicit none

  ! local parameters
  integer :: irec,irec_local,total_length
  integer :: idimval,jdimval,isample
  real(kind=CUSTOM_REAL) :: time_t

  character(len=4) :: chn
  character(len=1) :: component
  character(len=MAX_STRING_LEN) :: sisname

  ! using an ending .semd also for strain seismograms (since strain is derived from displacement)
  component = 'd'

  ! full length of strain seismogram (routine will be called at the very end of the time stepping)
  total_length = seismo_offset + seismo_current

  ! this write routine here is called at the very end, thus seismo_offset should be equal to NSTEP/subsamp_seismos
  ! checks
  if (total_length /= NSTEP/subsamp_seismos) then
    print *,'Error: writing strain seismograms has wrong index lengths ',total_length,'should be ',NSTEP/subsamp_seismos
    call exit_MPI(myrank,'Error writing strain seismograms, has wrong index length')
  endif

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

        ! create the name of the seismogram file for each slice
        ! file name includes the name of the station, the network and the component
        write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,chn,component

        ! save seismograms in text format (with no subsampling).
        ! Because we do not subsample the output, this can result in large files
        ! if the simulation uses many time steps. However, subsampling the output
        ! here would result in a loss of accuracy when one later convolves
        ! the results with the source time function
        open(unit=IOUT,file=OUTPUT_FILES(1:len_trim(OUTPUT_FILES))//sisname(1:len_trim(sisname)),status='unknown')

        ! make sure we never write more than the maximum number of time steps
        ! subtract half duration of the source to make sure travel time is correct
        do isample = 1,total_length
          ! time
          ! distinguish between single and double precision for reals
          time_t = real(dble(isample-1)*DT*subsamp_seismos - t0,kind=CUSTOM_REAL)
          ! output
          write(IOUT,*) time_t,' ',seismograms_eps(jdimval,idimval,irec_local,isample)
        enddo

        close(IOUT)

      enddo ! jdimval
    enddo ! idimval
  enddo ! irec_local

  end subroutine write_adj_seismograms2_to_file

!=====================================================================

  subroutine write_channel_name(iorientation,channel)

  use specfem_par, only: DT,SUPPRESS_UTM_PROJECTION

  implicit none

  integer :: iorientation
  character(len=3) :: channel

  ! local parameters
  character(len=2) :: bic
  double precision:: sampling_rate

  ! gets band and instrument code
  sampling_rate = DT
  call band_instrument_code(sampling_rate,bic)

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

  ! sets channel name
  if (SUPPRESS_UTM_PROJECTION) then
    ! no UTM, pure Cartesian reference
    ! uses Cartesian X/Y/Z direction to denote channel
    select case (iorientation)
    case (1)
      channel = bic(1:2)//'X'
    case (2)
      channel = bic(1:2)//'Y'
    case (3)
      channel = bic(1:2)//'Z'
    case (4)
      channel = bic(1:2)//'P'  ! for pressure seismograms
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  else
    ! UTM conversion
    ! uses convention for N/E/Z to denote channel
    select case (iorientation)
    case (1)
      channel = bic(1:2)//'E'
    case (2)
      channel = bic(1:2)//'N'
    case (3)
      channel = bic(1:2)//'Z'
    case (4)
      channel = bic(1:2)//'P'  ! for pressure seismograms
    case default
      call exit_mpi(0,'error channel orientation value')
    end select

  endif

  end subroutine write_channel_name

!=====================================================================

  subroutine band_instrument_code(DT,bic)
  ! This subroutine is to choose the appropriate band and instrument codes for channel names of seismograms
  ! based on the IRIS convention (first two letters of channel codes, respectively,
  ! which were LH(Z/E/N) previously).
  ! For consistency with observed data, we now use the IRIS convention for band codes (first letter in channel codes) of
  ! SEM seismograms governed by their sampling rate.
  ! Instrument code (second letter in channel codes) is fixed to "X" which is assigned by IRIS for synthetic seismograms.
  ! See the manual for further explanations!
  ! Ebru Bozdag, November 2010
  implicit none
  double precision :: DT
  character(len=2) :: bic
  ! local parameter
  logical,parameter :: SUPPRESS_IRIS_CONVENTION = .false.

  ! see manual for ranges
  if (DT >= 1.0d0)  bic = 'LX'
  if (DT < 1.0d0 .and. DT > 0.1d0) bic = 'MX'
  if (DT <= 0.1d0 .and. DT > 0.0125d0) bic = 'BX'
  if (DT <= 0.0125d0 .and. DT > 0.004d0) bic = 'HX'
  if (DT <= 0.004d0 .and. DT > 0.001d0) bic = 'CX'
  if (DT <= 0.001d0) bic = 'FX'

  ! ignores IRIS convention, uses previous, constant band and instrument code
  if (SUPPRESS_IRIS_CONVENTION) then
    bic = 'BH'
  endif

 end subroutine band_instrument_code

