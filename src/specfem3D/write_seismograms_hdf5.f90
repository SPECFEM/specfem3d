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

      ! writes out seismogram files
      select case(SIMULATION_TYPE)
      case (1,3)
        ! forward/kernel simulations
        ! forward/backward wavefield
        ! forward & kernel simulations
        if (NIONOD > 0) then
          if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
            call write_seismograms_to_file_h5(seismograms_d,1)
          if (SAVE_SEISMOGRAMS_VELOCITY) &
            call write_seismograms_to_file_h5(seismograms_v,2)
          if (SAVE_SEISMOGRAMS_ACCELERATION) &
            call write_seismograms_to_file_h5(seismograms_a,3)
          if (SAVE_SEISMOGRAMS_PRESSURE) &
            call write_seismograms_to_file_h5(seismograms_p,4)
        else
          if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
            call write_seismo_no_ioserv(seismograms_d,1)
          if (SAVE_SEISMOGRAMS_VELOCITY) &
            call write_seismo_no_ioserv(seismograms_v,2)
          if (SAVE_SEISMOGRAMS_ACCELERATION) &
            call write_seismo_no_ioserv(seismograms_a,3)
          if (SAVE_SEISMOGRAMS_PRESSURE) &
            call write_seismo_no_ioserv(seismograms_p,4)
        endif
      case (2)
        ! adjoint simulations
        ! adjoint wavefield
        ! adjoint simulations
        if (SAVE_SEISMOGRAMS_DISPLACEMENT) &
          call write_adj_seismograms_to_file_h5(seismograms_d,1)
        if (SAVE_SEISMOGRAMS_VELOCITY) &
          call write_adj_seismograms_to_file_h5(seismograms_v,2)
        if (SAVE_SEISMOGRAMS_ACCELERATION) &
          call write_adj_seismograms_to_file_h5(seismograms_a,3)
        if (SAVE_SEISMOGRAMS_PRESSURE) &
          call write_adj_seismograms_to_file_h5(seismograms_p,4)
        ! updates adjoint time counter
        it_adj_written = it
      end select

      ! synchronizes processes (waits for all processes to finish writing)
      call synchronize_all()

      ! user output
      if (myrank == 0) then
        ! timing
        write_time = wtime() - write_time_begin
        ! output
        write(IMAIN,*)
        write(IMAIN,*) 'Total number of time steps written: ', it-it_begin+1
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    endif ! do_save_seismograms

    ! resets current seismogram position
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0

  endif


  end subroutine write_seismograms_h5



  subroutine write_seismograms_to_file_h5(seismograms,istore)

  use constants

  use specfem_par, only: myrank,number_receiver_global,NPROC, &
          nrec,nrec_local,islice_selected_rec, &
          seismo_offset,seismo_current, &
          NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,ASDF_FORMAT, &
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_SEISMOGRAMS, NSTEP


  implicit none

  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms

  ! local parameters
  character(len=4) component

  integer :: nrec_local_received,total_seismos,req
  integer :: io_tag_seismo, time_window

  ! saves displacement, velocity, acceleration, or pressure
  if (istore == 1) then
    component = 'disp'
    io_tag_seismo = io_tag_seismo_body_disp
  else if (istore == 2) then
    component = 'velo'
    io_tag_seismo = io_tag_seismo_body_velo
  else if (istore == 3) then
    component = 'acce'
    io_tag_seismo = io_tag_seismo_body_acce
  else if (istore == 4) then
    component = 'pres'
    io_tag_seismo = io_tag_seismo_body_pres
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  if (nrec_local > 0) then
    ! send seismograms (instead of one seismo)
    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) then
      time_window = NSTEP
    else
      time_window = NTSTEP_BETWEEN_OUTPUT_SEISMOS
    endif

    if (istore /= 4) then
      call isend_cr_inter(seismograms(:,:,1:time_window),NDIM*nrec_local*time_window,0,io_tag_seismo,req)
    else
      call isend_cr_inter(seismograms(1,:,1:time_window),1*nrec_local*time_window,0,io_tag_seismo,req)
    endif
  endif

  end subroutine write_seismograms_to_file_h5


 subroutine write_seismo_no_ioserv(seismograms, istore)
    use specfem_par
    use phdf5_utils
    use io_server

    implicit none

    character(len=4) :: component
    integer :: t_upper, irec, irec_local, irec_glob

    integer :: istore
    real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable    :: seismograms_all_rec


    ! parameters for master collects seismograms
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram
    integer :: nrec_local_received,NPROCTOT,total_seismos,receiver,sender
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

    allocate(one_seismogram(NDIM,NSTEP),stat=ier)
    if (ier /= 0) stop 'error while allocating one temporary seismogram'

    if (myrank == 0) then
      ! on the master, gather all the seismograms
      total_seismos = 0

      ! loop on all the slices
      call world_size(NPROCTOT)

      ! counts number of local receivers for each slice
      allocate(islice_num_rec_local(0:NPROCTOT-1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'error allocating islice_num_rec_local')
      islice_num_rec_local(:) = 0
      do irec = 1,nrec
        iproc = islice_selected_rec(irec)
        islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
      enddo

      allocate(seismograms_all_rec(NDIM, sum(islice_num_rec_local(:)),NSTEP),stat=ier)
      if (ier /= 0) stop 'error while allocating all station dump (seismograms)'

      ! loop on all the slices
      do iproc = 0,NPROCTOT-1

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
              one_seismogram(:,:) = seismograms(:,irec_local,:)
            else
              call recv_i(tmp_irec,1,sender,itag)
              irec = tmp_irec(1)
              if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')

              call recvv_cr(one_seismogram,NDIM*NSTEP,sender,itag)
            endif

            ! gather all one_seismogram
            seismograms_all_rec(:,irec,:) = one_seismogram

            total_seismos = total_seismos + 1

          enddo ! nrec_local_received

        endif ! if (nrec_local_received > 0)
      enddo ! NPROCTOT-1
      deallocate(islice_num_rec_local)

      if (total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

    else

      if (nrec_local > 0) then
        ! on the nodes, send the seismograms to the master
        receiver = 0
        tmp_nrec_local(1) = nrec_local
        call send_i(tmp_nrec_local,1,receiver,itag)
        do irec_local = 1,nrec_local
          ! get global number of that receiver
          irec = number_receiver_global(irec_local)
          tmp_irec(1) = irec
          call send_i(tmp_irec,1,receiver,itag)

          ! sends seismogram of that receiver
          one_seismogram(:,:) = seismograms(:,irec_local,:)
          call sendv_cr(one_seismogram,NDIM*NSTEP,receiver,itag)
        enddo
      endif
    endif ! myrank

    if (myrank == 0) then
      ! check if the array length to be written > total timestep
      if (seismo_offset < 0) seismo_offset = 0
      if (seismo_offset+NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) then
        t_upper = NSTEP - seismo_offset
      else
        t_upper = NTSTEP_BETWEEN_OUTPUT_SEISMOS
      endif

      ! writes out this seismogram
      call h5_open_file(h5)

      !call write_one_seismogram(one_seismogram,irec,component,istore)
      if (istore == 1) then
        component = 'disp'
        call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
             seismograms_all_rec(:,:,1:t_upper), (/0, 0, seismo_offset/), .false.)

      else if (istore == 2) then
        component = 'velo'
        call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
             seismograms_all_rec(:,:,1:t_upper), (/0, 0, seismo_offset/), .false.)
      else if (istore == 3) then
          component = 'acce'
        call h5_write_dataset_3d_r_collect_hyperslab(h5, component, &
             seismograms_all_rec(:,:,1:t_upper), (/0, 0, seismo_offset/), .false.)
      else if (istore == 4) then
        component = 'pres'
        call h5_write_dataset_2d_r_collect_hyperslab(h5, component, &
             seismograms_all_rec(1,:,1:t_upper), (/0, seismo_offset/), .false.)
      else
        call exit_MPI(myrank,'wrong component to save for seismograms')
      endif

      call h5_close_file(h5)

      deallocate(seismograms_all_rec)
    endif ! myrank == 0
    deallocate(one_seismogram)

  end subroutine write_seismo_no_ioserv

!=====================================================================

! write adjoint seismograms (displacement) to text files

  subroutine write_adj_seismograms_to_file_h5(seismograms,istore)

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