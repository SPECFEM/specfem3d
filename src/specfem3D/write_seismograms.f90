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

  ! update position in seismograms
  seismo_current = seismo_current + 1

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

  ! gets seismogram values
  if (GPU_MODE) then
    ! on GPU
    if (nrec_local > 0 .and. do_save_seismograms) then
      ! gets resulting array values onto CPU
      call compute_seismograms_cuda(Mesh_pointer,seismograms_d,seismograms_v,seismograms_a,seismograms_p, &
                                    seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS,it,it_end, &
                                    ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,USE_TRICK_FOR_BETTER_PRESSURE)
    endif
  else
    ! on CPU
    if (nrec_local > 0 .and. do_save_seismograms) then
      call compute_seismograms()
    endif
  endif

  ! write the current or final seismograms
  if (seismo_current == NTSTEP_BETWEEN_OUTPUT_SEISMOS .or. it == it_end) then

    if (do_save_seismograms .and. .not. INVERSE_FWI_FULL_PROBLEM) then

      ! timing
      write_time_begin = wtime()

      ! checks if anything to do
      if (nrec_local > 0 .or. (WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0) .or. ASDF_FORMAT) then

        ! writes out seismogram files
        select case(SIMULATION_TYPE)
        case (1,3)
          if (.not. SU_FORMAT) then
            ! forward & kernel simulations
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_seismograms_to_file(seismograms_d,1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_seismograms_to_file(seismograms_v,2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_seismograms_to_file(seismograms_a,3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_seismograms_to_file(seismograms_p,4)
          else
            ! SU_format
            ! write ONE binary file for all receivers (nrec_local) within one proc
            ! SU format, with 240-byte-header for each trace
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_output_SU(seismograms_d,1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_output_SU(seismograms_v,2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_output_SU(seismograms_a,3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_output_SU(seismograms_p,4)
          endif
        case (2)
          ! adjoint wavefield
          if (.not. SU_FORMAT) then
            ! adjoint simulations
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_adj_seismograms_to_file(seismograms_d,1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_adj_seismograms_to_file(seismograms_v,2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_adj_seismograms_to_file(seismograms_a,3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_adj_seismograms_to_file(seismograms_p,4)
          else
            ! SU_format
            ! write ONE binary file for all receivers (nrec_local) within one proc
            ! SU format, with 240-byte-header for each trace
            if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
              call write_output_SU(seismograms_d,1)
            if (SAVE_SEISMOGRAMS_VELOCITY) &
              call write_output_SU(seismograms_v,2)
            if (SAVE_SEISMOGRAMS_ACCELERATION) &
              call write_output_SU(seismograms_a,3)
            if (SAVE_SEISMOGRAMS_PRESSURE) &
              call write_output_SU(seismograms_p,4)
          endif
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
        write(IMAIN,*) 'Total number of time steps written: ', it-it_begin+1
        if (WRITE_SEISMOGRAMS_BY_MASTER) then
          write(IMAIN,*) 'Writing the seismograms by master proc alone took ',sngl(write_time),' seconds'
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

  subroutine write_seismograms_to_file(seismograms,istore)

  use constants

  use specfem_par, only: myrank,number_receiver_global,NPROC, &
          nrec,nrec_local,islice_selected_rec, &
          seismo_offset,seismo_current, &
          NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,ASDF_FORMAT, &
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_SEISMOGRAMS

  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms

  ! local parameters
  integer :: irec,irec_local

  character(len=1) component

  ! parameters for master collects seismograms
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
  integer :: nrec_local_received,total_seismos,receiver,sender
  integer :: iproc,ier,i
  integer,dimension(1) :: tmp_nrec_local_received,tmp_irec,tmp_nrec_local
  integer,dimension(0:NPROC-1) :: islice_num_rec_local

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

  allocate(one_seismogram(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2420')
  if (ier /= 0) stop 'error while allocating one temporary seismogram'

  ! write out seismograms: all processes write their local seismograms themselves
  if (.not. WRITE_SEISMOGRAMS_BY_MASTER) then

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

    if (ASDF_FORMAT) call init_asdf_data(nrec_local)

    ! loop on all the local receivers
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      ! writes out this seismogram
      do i = 1,seismo_current
        one_seismogram(:,i) = seismograms(:,irec_local,i)
      enddo

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

  else if (WRITE_SEISMOGRAMS_BY_MASTER .and. ASDF_FORMAT) then
    call init_asdf_data(nrec_local)
    call synchronize_all()

    do irec_local = 1, nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)
      one_seismogram(:,:) = seismograms(:,irec_local,:)
      call write_one_seismogram(one_seismogram,irec_local,irec,component,istore)
    enddo
    call write_asdf()
    call synchronize_all()
    call close_asdf_data()

! only the master process does the writing of seismograms and
! collects the data from all other processes
  else
    ! only written out by master process - WRITE_SEISMOGRAMS_BY_MASTER

    if (myrank == 0) then
      ! on the master, gather all the seismograms

      ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        write(sisname,'(A)') '/all_seismograms'
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
          nrec_local_received = nrec_local
        endif

        if (nrec_local_received > 0) then
          do irec_local = 1,nrec_local_received
            ! receive except from proc 0, which is myself and therefore I already have these values
            if (iproc == 0) then
              ! get global number of that receiver
              irec = number_receiver_global(irec_local)
              do i = 1,seismo_current
                one_seismogram(:,i) = seismograms(:,irec_local,i)
              enddo
            else
              call recv_i(tmp_irec,1,sender,itag)
              irec = tmp_irec(1)
              if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')

              call recvv_cr(one_seismogram,NDIM*seismo_current,sender,itag)
            endif

            total_seismos = total_seismos + 1

            ! writes out this seismogram
            call write_one_seismogram(one_seismogram,irec_local,irec,component,istore)

          enddo ! nrec_local_received
        endif ! if (nrec_local_received > 0)
      enddo ! NPROC-1

      write(IMAIN,*) 'Component: .sem'//component
      write(IMAIN,*) '  total number of receivers saved is ',total_seismos,' out of ',nrec
      write(IMAIN,*)

      if (total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

     ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    else
      ! on the nodes, send the seismograms to the master
      receiver = 0
      tmp_nrec_local(1) = nrec_local
      call send_i(tmp_nrec_local,1,receiver,itag)
      if (nrec_local > 0) then
        do irec_local = 1,nrec_local
          ! get global number of that receiver
          irec = number_receiver_global(irec_local)
          tmp_irec(1) = irec
          call send_i(tmp_irec,1,receiver,itag)

          ! sends seismogram of that receiver
          do i = 1,seismo_current
            one_seismogram(:,i) = seismograms(:,irec_local,i)
          enddo
          call sendv_cr(one_seismogram,NDIM*seismo_current,receiver,itag)
        enddo
      endif
    endif ! myrank

  endif ! of if (WRITE_SEISMOGRAMS_BY_MASTER)

  deallocate(one_seismogram)

  end subroutine write_seismograms_to_file

!=====================================================================

  subroutine write_one_seismogram(one_seismogram,irec_local,irec,component,istore)

  use constants

  use specfem_par, only: DT,t0,it, &
    NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,ASDF_FORMAT, &
    SIMULATION_TYPE,station_name,network_name

  implicit none

  integer, intent(in) :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: one_seismogram
  integer,intent(in) :: irec_local,irec
  character(len=1),intent(in) :: component

  ! local parameters
  integer :: iorientation,number_of_components
  integer :: length_station_name,length_network_name
  character(len=MAX_STRING_LEN) :: sisname,final_LOCAL_PATH
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

    ! directory to store seismograms
    final_LOCAL_PATH = OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // '/'


    if (ASDF_FORMAT) then
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

  subroutine write_adj_seismograms_to_file(seismograms,istore)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  use specfem_par, only: myrank,number_receiver_global,nrec_local,it,DT, &
    NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    t0,seismo_offset,it_adj_written

  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismograms

  ! local parameters
  integer :: irec,irec_local
  integer :: iorientation,isample,number_of_components

  character(len=3) :: channel
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
      if (istore == 4) then ! this is for pressure
        call write_channel_name(istore,channel)
      else
        call write_channel_name(iorientation,channel)
      endif

      ! create the name of the seismogram file for each slice
      ! file name includes the name of the station, the network and the component
      write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem',a1)") '/NT','S',irec,channel,component

      ! save seismograms in text format with no subsampling.
      ! Because we do not subsample the output, this can result in large files
      ! if the simulation uses many time steps. However, subsampling the output
      ! here would result in a loss of accuracy when one later convolves
      ! the results with the source time function
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
      ! subtract half duration of the source to make sure travel time is correct
      do isample = it_adj_written+1,min(it,NSTEP)
        ! distinguish between single and double precision for reals
        write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ', &
                      seismograms(iorientation,irec_local,isample-it_adj_written)
      enddo

      close(IOUT)

    enddo

  enddo

  end subroutine write_adj_seismograms_to_file

!=====================================================================

! write adjoint seismograms (strain) to text files

  subroutine write_adj_seismograms2_to_file(myrank,seismograms,number_receiver_global,nrec_local,it,DT,NSTEP,t0)

  use constants, only: CUSTOM_REAL,NDIM,MAX_STRING_LEN,IOUT,OUTPUT_FILES

  implicit none
  integer :: myrank
  integer :: nrec_local,NSTEP,it
  integer, dimension(nrec_local) :: number_receiver_global
  ! note: seismograms here is still an array of size *,*,*,NSTEP
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local,NSTEP) :: seismograms
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
          write(IOUT,*) real(dble(isample-1)*DT - t0,kind=CUSTOM_REAL),' ',seismograms(jdimval,idimval,irec_local,isample)
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

