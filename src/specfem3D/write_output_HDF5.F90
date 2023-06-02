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

! writes seismogram output in HDF5 file format

  subroutine write_output_HDF5(all_seismograms,nrec_store,istore)

  use specfem_par

#ifdef USE_HDF5
  use manager_hdf5
  use io_server_hdf5
#endif

  implicit none

  ! arguments
  integer,intent(in) :: istore,nrec_store
  real(kind=CUSTOM_REAL), dimension(NDIM,nlength_seismogram,nrec_store),intent(in) :: all_seismograms

#ifdef USE_HDF5
  ! local parameters
  integer :: ier
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tmp_seis_pre   ! for pressure output
  character(len=MAX_STRING_LEN) :: fname_h5_seismo
  character(len=5) :: comp

  ! flag for file initialization
  logical, save :: is_initialized = .false.

  ! hdf5 i/o server
  integer :: req
  integer :: io_tag_seismo

! HDF5 writes into a single file (by main process only)

  ! io server
  if (HDF5_IO_NODES > 0) then
    ! only io nodes do the file output
    ! collect seismograms on io nodes
    ! saves displacement, velocity, acceleration, or pressure
    select case (istore)
    case (1)
      io_tag_seismo = io_tag_seismo_body_disp
    case (2)
      io_tag_seismo = io_tag_seismo_body_velo
    case (3)
      io_tag_seismo = io_tag_seismo_body_acce
    case (4)
      io_tag_seismo = io_tag_seismo_body_pres
    end select

    ! send seismograms
    if (nrec_store > 0) then
      ! sends full array, all-components
      call isend_cr_inter(all_seismograms(:,:,:),NDIM*nrec_store*nlength_seismogram,0,io_tag_seismo,req)
    endif

    ! all done
    return
  else
    ! no io server
    ! main process writes out all
    continue
  endif

  ! safety check
  if (.not. WRITE_SEISMOGRAMS_BY_MAIN) &
    stop 'HDF5_FORMAT must have WRITE_SEISMOGRAMS_BY_MAIN set to .true.'

  ! initializes
  ! we only want to do this once
  if (.not. is_initialized) then
    ! creates and initializes seismogram file
    ! initialization needs to be done only once
    ! (first time & first component when this function gets called)
    call write_output_hdf5_seismogram_init()
    ! sets flag is done
    is_initialized = .true.
  endif

  ! only main process writes out
  if (myrank /= 0) return

  ! saves displacement, velocity, acceleration, or pressure
  select case (istore)
  case (1)
    comp = 'displ'
  case (2)
    comp = 'veloc'
  case (3)
    comp = 'accel'
  case (4)
    comp = 'press'
  case default
    stop 'Invalid seismogram component to save for seismograms, not implemented yet'
  end select

  fname_h5_seismo = trim(OUTPUT_FILES) // 'seismograms.h5'

  ! initialze hdf5
  call h5_initialize()

  ! note:
  ! - the total array length for a seismogram is NSTEP / NTSTEP_BETWEEN_OUTPUT_SAMPLE
  !
  ! - nlength_seismogram : is the length for writing out (intermediate) portions of the seismograms
  !                        with a constant length of
  !                        nlength_seismogram == NTSTEP_BETWEEN_OUTPUT_SEISMOS / NTSTEP_BETWEEN_OUTPUT_SAMPLE
  !
  ! - seismo_offset      : indicates from which time step (it) on the seismogram is written, starting with offset == 0
  !
  ! - seismo_current     : indicates the (variable) length of this current portion,
  !                        which can be the full nlength_seismogram length or just a remaining portion
  !                        (e.g. for a last portion in case NTSTEP_BETWEEN_OUTPUT_SEISMOS is not a multiple of NSTEP)

  ! writes out this seismogram
  call h5_open_file(fname_h5_seismo)

  select case (istore)
  case (1)
    ! displ
    call h5_write_dataset_collect_hyperslab(comp, all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)
  case (2)
    ! veloc
    call h5_write_dataset_collect_hyperslab(comp, all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)
  case (3)
    ! accel
    call h5_write_dataset_collect_hyperslab(comp, all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)
  case (4)
    ! pressure
    allocate(tmp_seis_pre(seismo_current,nrec_store),stat=ier)
    if (ier /= 0) stop 'error allocating tmp_seis_pre array'
    tmp_seis_pre(:,:) = all_seismograms(1,1:seismo_current,:)

    call h5_write_dataset_collect_hyperslab(comp, tmp_seis_pre, (/seismo_offset, 0/), .false.)

    deallocate(tmp_seis_pre)

  case default
    stop 'Invalid seismogram component to save for seismograms, not implemented yet'
  end select

  ! finish writing
  call h5_close_file()
  call h5_finalize()

#else
  ! no HDF5 compilation support
  ! to avoid compiler warnings
  real(kind=CUSTOM_REAL) :: dummy
  integer :: i_dummy

  i_dummy = istore
  if (nrec_store > 0 .and. nlength_seismogram > 0) dummy = all_seismograms(1,1,1)

  ! user output
  print *
  print *, "Error: HDF5 routine write_output_HDF5() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *
  stop 'Error HDF5 write_output_HDF5(): called without compilation support'

#endif

  end subroutine write_output_HDF5

!
!------------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5
! only with HDF5 compilation support

  subroutine write_output_hdf5_seismogram_init()

  use specfem_par
  use manager_hdf5

  implicit none

  ! local parameters
  integer :: i,irec,ier
  integer :: nlength_total_seismogram
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: time_array
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: val_array2d
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: val_array3d
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rec_coords
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: stations
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: networks
  character(len=3) :: channel
  character(len=3),dimension(NDIM) :: channels
  character(len=3),dimension(1) :: channels_press
  character(len=MAX_STRING_LEN) :: fname_h5_seismo

  ! only main process creates file
  if (myrank /= 0) return

  fname_h5_seismo = trim(OUTPUT_FILES) // 'seismograms.h5'

  ! user output
  if (myrank == 0 .and. IO_compute_task) then
    write(IMAIN,*) 'Creating seismograms in HDF5 file format'
    if (WRITE_SEISMOGRAMS_BY_MAIN) then
      write(IMAIN,*) '  writing waveforms by main...'
    else
      write(IMAIN,*) '  writing waveforms in parallel...'
    endif
    write(IMAIN,*) '  seismogram file: ',trim(fname_h5_seismo)
    call flush_IMAIN()
  endif

  ! initialze hdf5
  call h5_initialize()

  ! create file
  call h5_create_file(fname_h5_seismo)

  ! total length of seismograms (considering subsampling)
  nlength_total_seismogram = NSTEP / NTSTEP_BETWEEN_OUTPUT_SAMPLE

  ! create time dataset it = 1 ~ NSTEP
  allocate(time_array(nlength_total_seismogram),stat=ier)
  if (ier /= 0) stop 'Error allocating time_array'
  time_array(:) = 0.0_CUSTOM_REAL

  do i = 1, nlength_total_seismogram
    if (SIMULATION_TYPE == 1) then
      ! forward simulation
      ! distinguish between single and double precision for reals
      time_array(i) = real( dble((i-1) * NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0 ,kind=CUSTOM_REAL)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint simulation: backward/reconstructed wavefields
      ! distinguish between single and double precision for reals
      ! note: compare time_t with time used for source term
      time_array(i) = real( dble((NSTEP-i) * NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0 ,kind=CUSTOM_REAL)
    endif
  enddo

  ! time array
  call h5_write_dataset_no_group("time", time_array)

  ! free array
  deallocate(time_array)

  ! read out_list_stations.txt generated at locate_receivers.f90:431 here to write in the h5 file.
  allocate(stations(nrec), &
           networks(nrec), &
           rec_coords(nrec,3),stat=ier)
  if (ier /= 0) stop 'Error allocating station arrays'
  rec_coords(:,:) = 0.0_CUSTOM_REAL

  open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'output_list_stations.txt', &
       status='unknown',action='read',iostat=ier)
  if (ier /= 0) then
    call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'output_list_stations.txt')
  endif
  ! reads station infos
  do irec = 1,nrec
    read(IOUT_SU,*) stations(irec),networks(irec), rec_coords(irec, 1), rec_coords(irec, 2), rec_coords(irec, 3)
  enddo
  ! closes output file
  close(IOUT_SU)

  ! coordination
  call h5_write_dataset_no_group("coords", rec_coords)
  ! station name
  call h5_write_dataset_no_group("station", stations)
  ! network name
  call h5_write_dataset_no_group("network", networks)

  ! channels info
  if (SAVE_SEISMOGRAMS_DISPLACEMENT .or. SAVE_SEISMOGRAMS_VELOCITY .or. SAVE_SEISMOGRAMS_ACCELERATION) then
    do i = 1,NDIM
      call write_channel_name(i,channel)
      channels(i) = channel
    enddo
    call h5_write_dataset_no_group("channel", channels)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    ! this is for pressure
    call write_channel_name(4,channel)
    channels_press(1) = channel
    call h5_write_dataset_no_group("channel_press", channels_press)
  endif

  ! free arrays
  deallocate(stations,networks,rec_coords)

  ! prepare datasets for physical values
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    ! 3-component displacement
    allocate(val_array3d(NDIM,nlength_total_seismogram,nrec),stat=ier)
    call h5_create_dataset_gen("displ", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif

  if (SAVE_SEISMOGRAMS_VELOCITY) then
    ! 3-component velocity
    allocate(val_array3d(NDIM,nlength_total_seismogram,nrec),stat=ier)
    call h5_create_dataset_gen("veloc", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif

  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    ! 3-component acceleration
    allocate(val_array3d(NDIM,nlength_total_seismogram,nrec),stat=ier)
    call h5_create_dataset_gen("accel", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif

  if (SAVE_SEISMOGRAMS_PRESSURE) then
    ! single component pressure
    allocate(val_array2d(nlength_total_seismogram,nrec),stat=ier)
    call h5_create_dataset_gen("press", shape(val_array2d), 2, CUSTOM_REAL)
    deallocate(val_array2d)
  endif

  call h5_close_file()
  call h5_finalize()

  end subroutine write_output_hdf5_seismogram_init

#endif


!-------------------------------------------------------------------------------------------------
!
! HDF5 i/o server routines
!
!-------------------------------------------------------------------------------------------------

#ifdef USE_HDF5

  subroutine write_seismograms_io_hdf5(it_offset)

  use constants, only: MAX_STRING_LEN,OUTPUT_FILES

  use shared_parameters, only: NSTEP, SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY, &
    SAVE_SEISMOGRAMS_ACCELERATION, SAVE_SEISMOGRAMS_PRESSURE

  use specfem_par, only: nlength_seismogram

  use io_server_hdf5
  use manager_hdf5

  implicit none

  integer, intent(in) :: it_offset
  integer :: t_upper
  character(len=MAX_STRING_LEN) :: fname_h5_seismo
  character(len=5) comp

  ! initializes
  fname_h5_seismo = trim(OUTPUT_FILES) // 'seismograms.h5'

  ! initialze hdf5
  call h5_initialize()

  ! writes out to file
  call h5_open_file(fname_h5_seismo)

  ! check if the array length to be written > total timestep
  if (it_offset+nlength_seismogram > NSTEP) then
    t_upper = nlength_seismogram - (it_offset+nlength_seismogram - NSTEP)
  else
    t_upper = nlength_seismogram
  endif

  ! writes out this seismogram
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    comp = 'displ'
    call h5_write_dataset_collect_hyperslab(comp, seismo_displ(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    comp = 'veloc'
    call h5_write_dataset_collect_hyperslab(comp, seismo_veloc(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    comp = 'accel'
    call h5_write_dataset_collect_hyperslab(comp, seismo_accel(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    comp = 'press'
    call h5_write_dataset_collect_hyperslab(comp, seismo_press(1:t_upper,:), (/it_offset, 0/), .false.)
  endif

  call h5_close_file()
  call h5_finalize()

  end subroutine write_seismograms_io_hdf5

#endif
