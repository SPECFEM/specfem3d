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
#endif

  implicit none

  ! arguments
  integer,intent(in) :: istore,nrec_store
  real(kind=CUSTOM_REAL), dimension(NDIM,nlength_seismogram,nrec_store),intent(in) :: all_seismograms

#ifdef USE_HDF5
  ! local parameters
  integer :: ier
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tmp_seis_pre   ! for pressure output
  character(len=MAX_STRING_LEN) :: filename
  character(len=5) :: component

  ! flag for file initialization
  logical, save :: is_initialized = .false.

! HDF5 writes into a single file (by main process only)

  ! io server
  if (HDF5_IO_NNODES > 0) then
    ! only io nodes do the file output
    ! collect seismograms on io nodes
    stop 'HDF5 IO server not implemented yet for HDF5 file output'
  else
    ! no io server
    ! main process writes out all
    continue
  endif

  ! safety check
  if (.not. WRITE_SEISMOGRAMS_BY_MAIN) &
    stop 'HDF5_FORMAT must have WRITE_SEISMOGRAMS_BY_MAIN set to .true.'

  ! initializes
  filename = trim(OUTPUT_FILES) // 'seismograms.h5'

  ! we only want to do this once
  if (.not. is_initialized) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Creating seismograms in HDF5 file format'
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        write(IMAIN,*) '  writing waveforms by main...'
      else
        write(IMAIN,*) '  writing waveforms in parallel...'
      endif
      write(IMAIN,*) '  seismogram file: ',trim(filename)
      call flush_IMAIN()
    endif

    ! creates and initializes seismogram file
    ! initialization needs to be done only once
    ! (first time & first component when this function gets called)
    call write_output_hdf5_seismogram_init(filename)

    ! sets flag is done
    is_initialized = .true.
  endif

  ! only main process writes out
  if (myrank /= 0) return

  ! saves displacement, velocity, acceleration, or pressure
  select case (istore)
  case (1)
    component = 'displ'
  case (2)
    component = 'veloc'
  case (3)
    component = 'accel'
  case (4)
    component = 'press'
  case default
    stop 'Invalid seismogram component to save for seismograms, not implemented yet'
  end select

  ! initialze hdf5
  call h5_init()

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
  call h5_open_file(filename)

  select case (istore)
  case (1)
    ! displ
    call h5_write_dataset_collect_hyperslab(component, &
                                                 all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)

  case (2)
    ! veloc
    call h5_write_dataset_collect_hyperslab(component, &
                                                 all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)
  case (3)
    ! accel
    call h5_write_dataset_collect_hyperslab(component, &
                                                 all_seismograms(:,1:seismo_current,:), (/0, seismo_offset, 0/), .false.)
  case (4)
    ! pressure
    component = 'pres'
    allocate(tmp_seis_pre(seismo_current,nrec_store),stat=ier)
    if (ier /= 0) stop 'error allocating tmp_seis_pre array'
    tmp_seis_pre(:,:) = all_seismograms(1,1:seismo_current,:)
    call h5_write_dataset_collect_hyperslab(component, &
                                                 tmp_seis_pre, (/seismo_offset, 0/), .false.)
    deallocate(tmp_seis_pre)

  case default
    stop 'Invalid seismogram component to save for seismograms, not implemented yet'
  end select

  ! finish writing
  call h5_close_file()
  call h5_destructor()

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

  subroutine write_output_hdf5_seismogram_init(filename)

  use specfem_par
  use manager_hdf5

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename

  ! local parameters
  integer :: i,irec,ier
  integer :: nlength_total_seismogram
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: time_array
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: val_array2d
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: val_array3d
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rec_coords
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: stations
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: networks

  ! only main process creates file
  if (myrank /= 0) return

  ! initialze hdf5
  call h5_init()

  ! create file
  call h5_create_file(filename)

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
  call h5_close_dataset()

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
  call h5_close_dataset()

  ! station name
  call h5_write_dataset_no_group("station", stations)
  call h5_close_dataset()

  ! network name
  call h5_write_dataset_no_group("network", networks)
  call h5_close_dataset()

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
  call h5_destructor()

  end subroutine write_output_hdf5_seismogram_init

#endif
