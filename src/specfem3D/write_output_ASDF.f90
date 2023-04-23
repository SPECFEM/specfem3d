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

!> Initializes the data structure for ASDF
!! \param nrec_local The number of receivers on the local processor
  subroutine init_asdf_data(nrec_local)

  use constants, only: NDIM
  use specfem_par, only: myrank,NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS

  use asdf_data, only: asdf_container

  implicit none
  ! Parameters
  integer,intent(in) :: nrec_local

  ! Variables
  integer :: total_seismos_local, ier

  ! note: partial seismogram can have length equal to NTSTEP_BETWEEN_OUTPUT_SEISMOS / NTSTEP_BETWEEN_OUTPUT_SAMPLE
  !       ASDF records have the full (subsampled) length NSTEP / NTSTEP_BETWEEN_OUTPUT_SAMPLE
  !
  ! check that full trace length (subsampled) and partial trace length (subsampled) match
  ! that is, must have NTSTEP_BETWEEN_OUTPUT_SEISMOS == NSTEP setting for now.
  ! we don't support partial outputs to ASDF files yet...
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS /= NSTEP) then
    print *,'Error: ASDF trace lengths mismatch'
    print *,'Please check if the setting in Par_file has NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP'
    call exit_MPI(myrank,'error ASDF trace length mismatch')
  endif

  total_seismos_local = nrec_local * NDIM ! 3 components

  asdf_container%nrec_local = nrec_local

  allocate(asdf_container%receiver_name_array(nrec_local), &
           asdf_container%network_array(nrec_local), &
           asdf_container%component_array(total_seismos_local), &
           asdf_container%receiver_lat(nrec_local), &
           asdf_container%receiver_lo(nrec_local), &
           asdf_container%receiver_el(nrec_local), &
           asdf_container%receiver_dpt(nrec_local), &
           asdf_container%records(total_seismos_local),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2015')
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')

  ! initializes
  asdf_container%receiver_name_array(:) = ""
  asdf_container%network_array(:) = ""
  asdf_container%component_array(:) = ""

  asdf_container%receiver_lat(:) = 0.0
  asdf_container%receiver_lo(:) = 0.0
  asdf_container%receiver_el(:) = 0.0
  asdf_container%receiver_dpt(:) = 0.0

  end subroutine init_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Stores the records into the ASDF structure
!! \param seismogram_tmp The current seismogram to store
!! \param irec_local The local index of the receiver
!! \param irec The global index of the receiver
!! \param chn The broadband channel simulated
!! \param iorientation The recorded seismogram's orientation direction
  subroutine store_asdf_data(seismogram_tmp, irec_local, irec, chn, iorientation)

  use constants, only: CUSTOM_REAL,NDIM,myrank

  use specfem_par, only: NSTEP,NTSTEP_BETWEEN_OUTPUT_SAMPLE,nlength_seismogram, &
    stlat,stlon,stele,stbur,station_name,network_name

  use asdf_data, only: asdf_container

  implicit none

  ! Parameters
  character(len=3),intent(in) :: chn
  integer,intent(in) :: irec_local, irec, iorientation
  real(kind=CUSTOM_REAL),dimension(NDIM,nlength_seismogram),intent(in) :: seismogram_tmp

  ! local Variables
  integer :: length_station_name, length_network_name
  integer :: ier, i, index_increment
  integer :: seismo_current_used

  ! actual seismogram length
  seismo_current_used = ceiling(real(NSTEP) / NTSTEP_BETWEEN_OUTPUT_SAMPLE)

  index_increment = iorientation

  ! trace index
  i = (irec_local-1)*(NDIM) + (index_increment)

  length_station_name = len_trim(station_name(irec))
  length_network_name = len_trim(network_name(irec))

  asdf_container%receiver_name_array(irec_local) = station_name(irec)(1:length_station_name)
  asdf_container%network_array(irec_local) = network_name(irec)(1:length_network_name)
  asdf_container%component_array(i) = chn(1:3)

  asdf_container%receiver_lat(irec_local) = stlat(irec)
  asdf_container%receiver_lo(irec_local) = stlon(irec)
  asdf_container%receiver_el(irec_local) = stele(irec)
  asdf_container%receiver_dpt(irec_local) = stbur(irec)

  allocate(asdf_container%records(i)%record(seismo_current_used),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2017')
  if (ier /= 0) call exit_MPI (myrank, 'Allocating ASDF container failed.')

  ! initializes
  asdf_container%records(i)%record(:) = 0.0

  ! seismogram as real data
  asdf_container%records(i)%record(1:seismo_current_used) = real(seismogram_tmp(iorientation, 1:seismo_current_used),kind=4)

  end subroutine store_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Closes the ASDF data structure by deallocating all arrays
  subroutine close_asdf_data()

  use constants, only: NDIM
  use asdf_data, only: asdf_container

  implicit none

  ! local Variables
  integer :: i

  do i = 1, asdf_container%nrec_local * NDIM ! 3 components
    ! skips component if records was not allocated for this component
    ! (valid entries must have non-empty component name)
    if (len_trim(asdf_container%component_array(i)) == 0) cycle
    ! free memory
    deallocate(asdf_container%records(i)%record)
  enddo

  deallocate(asdf_container%receiver_name_array)
  deallocate(asdf_container%network_array)
  deallocate(asdf_container%component_array)

  end subroutine close_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to the file
  subroutine write_asdf()

  use constants, only: NDIM,itag,MAX_LENGTH_NETWORK_NAME,MAX_LENGTH_STATION_NAME,IMAIN,myrank

  ! for ASDF
  use constants, only: ASDF_OUTPUT_PROVENANCE,ASDF_MAX_STRING_LENGTH, &
    ASDF_MAX_PARFILE_LENGTH,ASDF_MAX_QUAKEML_LENGTH,ASDF_MAX_STATIONXML_LENGTH, &
    ASDF_MAX_CONSTANTS_LENGTH,ASDF_MAX_TIME_STRING_LENGTH

  use asdf_data, only: asdf_container

  use iso_c_binding, only: C_NULL_CHAR !,c_ptr
!  use iso_Fortran_env

  use specfem_par, only: seismo_offset,DT,NSTEP,NTSTEP_BETWEEN_OUTPUT_SAMPLE,OUTPUT_FILES,WRITE_SEISMOGRAMS_BY_MAIN

  implicit none

  !--- Character strings to be written to the ASDF file
  character(len=ASDF_MAX_QUAKEML_LENGTH) :: quakeml
  character(len=ASDF_MAX_STATIONXML_LENGTH) :: stationxml

  integer :: num_stations
  integer :: stationxml_length
  integer :: nsamples  ! constant, as in SPECFEM
  integer :: seismo_current_used
  double precision :: sampling_rate
  double precision :: startTime
  integer(kind=8) :: start_time

  ! Network names and station names are two different beast, as in SPECFEM
  ! network_names(i) is related to station_names(i)
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: networks_names
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: stations_names
  character(len=3), dimension(:), allocatable :: component_names

  ! data. dimension = nsamples * num_channels_per_station * num_stations

  !-- ASDF variables
  ! daniel: note that these are file pointers to hid_t which (at least in newer HDF5 versions) become int64
  !         the corresponding Fortran integer would be integer(kind=8).
  !         one might have to check if the compiler uses 64-bit or 32-bit pointers.
  integer(kind=8) :: file_id
  integer(kind=8), dimension(3) ::  data_ids
  !   These variables are used to know where further writes should be done.
  !   They have to be cleaned as soon as they become useless
  integer(kind=8) :: waveforms_grp  ! Group "/Waveforms/"
  integer(kind=8) :: station_grp
  integer(kind=8) :: stationxml_grp

  integer :: current_proc, sender, receiver
  real, dimension(:,:), allocatable :: one_seismogram

  !--- MPI variables
  integer :: mysize, comm
  !--- Loop variables
  integer :: i, j, k, l
  !--- Error variable
  integer :: ier

  !--- 'allgather' arrays. Variables that needs to be known by everyone in
  !    order to define ASDF groups and datasets or write them as attributes.
  integer, dimension(:), allocatable :: num_stations_gather
  integer :: max_num_stations_gather
  character(len=MAX_LENGTH_STATION_NAME), dimension(:,:), allocatable :: station_names_gather
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:,:), allocatable :: network_names_gather
  character(len=3), dimension(:,:), allocatable :: component_names_gather
  real, dimension(:,:), allocatable :: station_lats_gather, station_longs_gather, station_elevs_gather, station_depths_gather
  integer, dimension(:), allocatable :: displs, rcounts

  ! temporary name built from network, station and channel names.
  character(len=ASDF_MAX_STRING_LENGTH) :: waveform_name,group_name

  ! C/Fortran interop for C-allocated strings
  integer :: len_prov, len_constants, len_Parfile
  !type(c_ptr) :: cptr
  !character, pointer :: fptr(:)
  character(len=2048) :: fptr

  character, dimension(:), allocatable, TARGET :: provenance
  character(len=ASDF_MAX_CONSTANTS_LENGTH) :: sf_constants
  character(len=ASDF_MAX_PARFILE_LENGTH) :: sf_parfile

  ! Time variables
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: start_time_string, end_time_string, &
                                                cmt_start_time_string, pde_start_time_string

  ! alias MPI communicator
  call world_duplicate(comm)
  call world_size(mysize)

  ! actual seismogram length
  seismo_current_used = ceiling(real(NSTEP) / NTSTEP_BETWEEN_OUTPUT_SAMPLE)

  num_stations = asdf_container%nrec_local
  sampling_rate = 1.0/(DT*NTSTEP_BETWEEN_OUTPUT_SAMPLE)

  ! BS BS: The total number of samples to be written to the ASDF file should be NSTEP.
  !        seismo_current is equivalent to NTSTEP_BETWEEN_OUTPUT_SEISMOS (except when it==it_end), as only in this case
  !        the routine is called from write_seismograms()
  !nsamples = seismo_current * (NSTEP / NTSTEP_BETWEEN_OUTPUT_SEISMOS)
  nsamples = ceiling(real(NSTEP) / NTSTEP_BETWEEN_OUTPUT_SAMPLE)

  ! Calculate start_time
  call get_time_cmt(startTime, start_time_string, pde_start_time_string, cmt_start_time_string, end_time_string)
  start_time = startTime*(int(1000000000,kind=8)) ! convert to nanoseconds

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------

  ! Generate minimal QuakeML for SPECFEM3D
  call cmt_to_quakeml(quakeml, pde_start_time_string, cmt_start_time_string)

  if (ASDF_OUTPUT_PROVENANCE) then
    ! Generate specfem provenance string
    ! see: https://seismicdata.github.io/SEIS-PROV/index.html

    ! function interface has been removed in ASDF versions >= 1.0.0
    !call ASDF_generate_sf_provenance_f(trim(start_time_string)//C_NULL_CHAR, &
    !                                   trim(end_time_string)//C_NULL_CHAR, cptr, len_prov)
    !call c_f_pointer(cptr, fptr, [len_prov])

    ! creates simple provenance string
    ! following asdf-library routine generate_sf_provenance()
    ! header
    fptr = &
      '<?xml version="1.0" encoding="UTF-8"?>'// &
      '<prov:document xmlns:prov="http://www.w3.org/ns/prov#"'// &
      'xmlns:seis_prov="http://seisprov.org/seis_prov/0.1/#" xmlns:xsd="http://www.w3.org/2001/XMLSchema"'// &
      'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'

    ! software name ("random" software_id = 0123456)
    fptr = trim(fptr) // &
      '<prov:softwareAgent prov:id="seis_prov:sp000_sa_' // '0123456' // '">'// &
      '  <prov:label>SPECFEM3D_Cartesian</prov:label>'// &
      '  <seis_prov:software_name>SPECFEM3D_Cartesian</seis_prov:software_name>'// &
      '  <seis_prov:software_version>3.0.0</seis_prov:software_version>'// &
      '  <seis_prov:website>https://github.com/SPECFEM/specfem3d</seis_prov:website>'// &
      '</prov:softwareAgent>'

    ! constants (constants_id = abcdefg)
    fptr = trim(fptr) // &
      '<prov:entity prov:id="seis_prov:sp000_fi_' // 'abcdefg' // '">'// &
      '  <prov:label>File</prov:label>'// &
      '  <prov:type xsi:type="xsd:string">seis_prov:file</prov:type>'// &
      '  <seis_prov:filename>constants.h</seis_prov:filename>'// &
      '  <seis_prov:location>/AuxiliaryData/Files/constants_h</seis_prov:location>'// &
      '  <seis_prov:location_type>HDF5 Data Set</seis_prov:location_type>'// &
      '</prov:entity>'

    ! Par_file (parfile_id = zyxwvut)
    fptr = trim(fptr) // &
      '<prov:entity prov:id="seis_prov:sp000_fi_' // 'zyxwvut' // '">'// &
      '  <prov:label>File</prov:label>'// &
      '  <prov:type xsi:type="xsd:string">seis_prov:file</prov:type>'// &
      '  <seis_prov:filename>Parfile</seis_prov:filename>'// &
      '  <seis_prov:location>/AuxiliaryData/Files/Parfile</seis_prov:location>'// &
      '  <seis_prov:location_type>HDF5 Data Set</seis_prov:location_type>'// &
      '</prov:entity>'

    ! simulation (simulation_id = 9876543)
    fptr = trim(fptr) // &
      '<prov:activity prov:id="seis_prov:sp000_ws_' // '9876543' // '">'// &
      '  <prov:startTime>' // trim(start_time_string) // '</prov:startTime>'// &
      '  <prov:endTime>' // trim(end_time_string) // '</prov:endTime>'// &
      '  <prov:label>Waveform Simulation</prov:label>'// &
      '  <prov:type xsi:type="xsd:string">seis_prov:waveform_simulation</prov:type>'// &
      '</prov:activity>'

    ! association
    fptr = trim(fptr) // &
      '<prov:wasAssociatedWith>'// &
      '  <prov:agent prov:ref="seis_prov:sp000_sa_' // '0123456' // '"/>'// &
      '  <prov:activity prov:ref="seis_prov:sp000_ws_' // '9876543' // '"/>'// &
      '</prov:wasAssociatedWith>'

    ! closing
    fptr = trim(fptr) // &
      '</prov:document>'

    len_prov = len_trim(fptr)

    allocate(provenance(len_prov+1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2018')
    provenance(1:len_prov) = fptr(1:len_prov)
    provenance(len_prov+1) = C_NULL_CHAR
  endif

  allocate(networks_names(num_stations), &
           stations_names(num_stations), &
           component_names(num_stations*NDIM), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating arrays 2019')
  networks_names(:) = ""; stations_names(:) = ""; component_names(:) = ""

  !--------------------------------------------------------
  ! ASDF variables
  !--------------------------------------------------------
  ! Find how many stations are managed by each allgatheress
  allocate(num_stations_gather(mysize),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2022')
  num_stations_gather(:) = 0

  call all_gather_all_i(num_stations, num_stations_gather, mysize)

  ! find the largest number of stations per allgatheress
  max_num_stations_gather = maxval(num_stations_gather)

  allocate(displs(mysize), &
           rcounts(mysize),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2024')
  displs(:) = 0; rcounts(:) = 0

  ! Everyone should know about each and every station name and its coordinates
  allocate(station_names_gather(max_num_stations_gather,mysize), &
           network_names_gather(max_num_stations_gather,mysize), &
           station_lats_gather(max_num_stations_gather,mysize), &
           station_longs_gather(max_num_stations_gather,mysize), &
           station_elevs_gather(max_num_stations_gather,mysize), &
           station_depths_gather(max_num_stations_gather,mysize), &
           component_names_gather(max_num_stations_gather*NDIM, mysize),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2031')
  station_names_gather(:,:) = ""; network_names_gather(:,:) = ""; component_names_gather(:,:) = ""

  ! This needs to be done because asdf_data is a pointer
  do i = 1, num_stations
    write(networks_names(i), '(a)') asdf_container%network_array(i)
    write(stations_names(i), '(a)') asdf_container%receiver_name_array(i)
  enddo

  do i = 1, num_stations * NDIM
    write(component_names(i), '(a)') asdf_container%component_array(i)
  enddo

  ! The number of stations is not constant across processes
  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * MAX_LENGTH_STATION_NAME
    rcounts(i) = num_stations_gather(i) * MAX_LENGTH_STATION_NAME
  enddo

  call all_gather_all_ch(stations_names, &
                         num_stations * MAX_LENGTH_STATION_NAME, &
                         station_names_gather, &
                         rcounts, &
                         displs, &
                         max_num_stations_gather, &
                         MAX_LENGTH_STATION_NAME, &
                         mysize)

  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * MAX_LENGTH_NETWORK_NAME
    rcounts(i) = num_stations_gather(i) * MAX_LENGTH_NETWORK_NAME
  enddo

  call all_gather_all_ch(networks_names, &
                         num_stations * MAX_LENGTH_NETWORK_NAME, &
                         network_names_gather, &
                         rcounts, &
                         displs, &
                         max_num_stations_gather, &
                         MAX_LENGTH_NETWORK_NAME, &
                         mysize)

  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * NDIM
    rcounts(i) = num_stations_gather(i) * NDIM
  enddo

  call all_gather_all_ch(component_names, &
                         num_stations*NDIM*NDIM, &
                         component_names_gather, &
                         rcounts*NDIM, &
                         displs*NDIM, &
                         max_num_stations_gather*NDIM, &
                         NDIM, &
                         mysize)

  ! Now gather all the coordiante information for these stations
  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather
    rcounts(i) = num_stations_gather(i)
  enddo

  call all_gather_all_r(asdf_container%receiver_lat, &
                        num_stations, &
                        station_lats_gather, &
                        rcounts, &
                        displs, &
                        max_num_stations_gather, &
                        mysize)
  call all_gather_all_r(asdf_container%receiver_lo, &
                        num_stations, &
                        station_longs_gather, &
                        rcounts, &
                        displs, &
                        max_num_stations_gather, &
                        mysize)
  call all_gather_all_r(asdf_container%receiver_el, &
                        num_stations, &
                        station_elevs_gather, &
                        rcounts, &
                        displs, &
                        max_num_stations_gather, &
                        mysize)
  call all_gather_all_r(asdf_container%receiver_dpt, &
                        num_stations, &
                        station_depths_gather, &
                        rcounts, &
                        displs, &
                        max_num_stations_gather, &
                        mysize)

  deallocate(stations_names)
  deallocate(networks_names)
  deallocate(component_names)
  deallocate(displs)
  deallocate(rcounts)

  allocate(one_seismogram(NDIM,seismo_current_used),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2032')
  if (ier /= 0) call exit_MPI(myrank,'error allocating one_seismogram array')
  one_seismogram(:,:) = 0.0

  !--------------------------------------------------------
  ! write ASDF
  !--------------------------------------------------------

  ! we only want to do these steps one time
  if (seismo_offset == 0) then
    if (myrank == 0) then
      ! user output
      write(IMAIN,*) 'creating ASDF file: ',trim(OUTPUT_FILES) // "synthetic.h5"
      call flush_IMAIN()

      ! initialize HDF5
      call ASDF_initialize_hdf5_f(ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF initialize hdf5 failed')

      ! creates file for synthetics
      call ASDF_create_new_file_serial_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, file_id)

      ! checks file identifier
      if (file_id == 0) then
        print *,'Error: ASDF failed to create new file',trim(OUTPUT_FILES) // "synthetic.h5"
        print *,'Please check if your installation supports asdf-library calls... exiting'
        call exit_MPI(myrank,'Error creating ASDF file for synthetics')
      endif

      ! sets file attributes
      call ASDF_write_string_attribute_f(file_id, "file_format" // C_NULL_CHAR, &
                                         "ASDF" // C_NULL_CHAR, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF attribute file_format failed')

      call ASDF_write_string_attribute_f(file_id, "file_format_version" // C_NULL_CHAR, &
                                         "1.0.0" // C_NULL_CHAR, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF attribute file_format_version failed')

      ! sets event info in quakeML format
      call ASDF_write_quakeml_f(file_id, trim(quakeml) // C_NULL_CHAR, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF write quakeml failed')

      ! provenance
      if (ASDF_OUTPUT_PROVENANCE) then
        ! keeps track of all setup choices
        call ASDF_write_provenance_data_f(file_id, provenance(1:len_prov+1), ier)
        if (ier /= 0) call exit_MPI(myrank,'Error ASDF write provenance failed')

        call read_file("setup/constants.h", sf_constants, len_constants)
        call read_file("DATA/Par_file", sf_parfile, len_Parfile)

        call ASDF_write_auxiliary_data_f(file_id, trim(sf_constants) // C_NULL_CHAR, &
                                         trim(sf_parfile(1:len_Parfile)) // C_NULL_CHAR, ier)
        if (ier /= 0) call exit_MPI(myrank,'Error ASDF write auxiliary data failed')
      endif

      ! waveforms
      call ASDF_create_waveforms_group_f(file_id, waveforms_grp)

      do k = 1, mysize ! Need to set up metadata for all processes
        do j = 1, num_stations_gather(k) ! loop over number of stations on that process
          ! group name
          write(group_name, '(a)') trim(network_names_gather(j, k)) // "." // trim(station_names_gather(j, k))

          ! group for station
          call ASDF_create_stations_group_f(waveforms_grp, &
                                            trim(group_name) // C_NULL_CHAR, &
                                            station_grp)

          ! gathers station info
          call station_to_stationxml(station_names_gather(j,k), network_names_gather(j,k), &
                                     station_lats_gather(j,k), station_longs_gather(j,k), &
                                     station_elevs_gather(j,k), station_depths_gather(j,k), &
                                     start_time_string, stationxml)

          ! allocates memory block for stationXML info
          stationxml_length = len_trim(stationxml) + 1  ! plus one for the additional null character
          call ASDF_define_station_xml_f(station_grp, stationxml_length, stationxml_grp)

          ! writes stationXML
          call ASDF_write_station_xml_f(stationxml_grp, trim(stationxml)//C_NULL_CHAR, ier)
          if (ier /= 0) call exit_MPI(myrank,'Error ASDF write station XML failed')

          ! waveforms
          ! loop over each component
          do  i = 1,NDIM
            ! note: if only pressure is output, there is only a single valid component.
            !       the other component_name entries will be empty (""), which will lead to an error
            !       with the HDF5 output when two waveforms with the same name will be defined.
            !       we therefore skip components with an empty ("") entry.
            if (len_trim(component_names_gather(i+((j-1)*NDIM),k)) == 0) cycle

            ! Generate unique waveform name
            ! example: HT.LIT.S3.MXN__2008-01-06T05:14:19__2008-01-06T05:14:53
            write(waveform_name, '(a)') trim(network_names_gather(j,k)) // "." // &
                                        trim(station_names_gather(j,k)) // ".S3." // &
                                        trim(component_names_gather(i+((j-1)*NDIM),k)) //"__"// &
                                        trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"

            call ASDF_define_waveform_f(station_grp, &
                                        nsamples, start_time, sampling_rate, &
                                        "Shot_ASDF" // C_NULL_CHAR, &
                                        trim(waveform_name) // C_NULL_CHAR, &
                                        data_ids(i))

            call ASDF_close_dataset_f(data_ids(i), ier)
            if (ier /= 0) call exit_MPI(myrank,'Error ASDF close dataset waveform failed')
          enddo

          call ASDF_close_dataset_f(stationxml_grp, ier)
          if (ier /= 0) call exit_MPI(myrank,'Error ASDF close dataset failed')

          call ASDF_close_group_f(station_grp, ier)
          if (ier /= 0) call exit_MPI(myrank,'Error ASDF close station group failed')

        enddo
      enddo

      call ASDF_close_group_f(waveforms_grp, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF close waveforms group failed')

      call ASDF_close_file_f(file_id, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF close file failed')

      call ASDF_finalize_hdf5_f(ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF finalize hdf5 failed')

    endif ! (myrank == 0)
  endif ! (seismo_offset == 0)

  ! synchronizes MPI processes
  call synchronize_all()

  ! Now write waveforms
  if (WRITE_SEISMOGRAMS_BY_MAIN) then

    ! only main writes
    if (myrank == 0) then
      ! user output
      write(IMAIN,*) 'writing waveforms by main...'
      write(IMAIN,*)
      call flush_IMAIN()

      call ASDF_initialize_hdf5_f(ier);
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF initialize hdf5 seismogram failed')

      call ASDF_open_serial_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, file_id)
      ! checks file identifier
      if (file_id == 0) then
        print *,'Error: ASDF failed to open file',trim(OUTPUT_FILES) // "synthetic.h5"
        call exit_MPI(myrank,'Error opening ASDF file for synthetics')
      endif

      call ASDF_open_waveforms_group_f(file_id, waveforms_grp)
    endif

    do k = 1, mysize ! Need to write data from all processes

      current_proc = k - 1
      sender = current_proc
      receiver = 0 ! the main proc does all the writing

      do j = 1, num_stations_gather(k) ! loop over number of stations on that process

        l = (j-1)*(NDIM) ! Index of current receiver in asdf_container%records

        ! First get the information to the main proc
        if (current_proc == 0) then
          ! current_proc is main proc
          if (myrank == 0) then
            do i = 1, NDIM
              ! note: if only pressure is output, there is only a single valid component.
              !       the other component_name entries will be empty (""), which will lead to an error
              !       with the HDF5 output when two waveforms with the same name will be defined.
              !       we therefore skip components with an empty ("") entry.
              if (len_trim(component_names_gather(i+((j-1)*NDIM),k)) == 0) cycle

              !  write(*,*) j, l, l+i, size(asdf_container%records)
              one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current_used)
            enddo
          endif

        else
          ! current_proc is not main proc
          if (myrank == current_proc) then
            !one_seismogram(:,:) = seismograms(:,j,:)
            do i = 1, NDIM
              ! note: if only pressure is output, there is only a single valid component.
              !       the other component_name entries will be empty (""), which will lead to an error
              !       with the HDF5 output when two waveforms with the same name will be defined.
              !       we therefore skip components with an empty ("") entry.
              if (len_trim(component_names_gather(i+((j-1)*NDIM),k)) == 0) cycle

              one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current_used)
            enddo
            ! send (real) data
            call send_r(one_seismogram,NDIM*seismo_current_used,receiver,itag)

          else if (myrank == 0) then
            ! receive (real) data
            call recv_r(one_seismogram,NDIM*seismo_current_used,sender,itag)

          endif
        endif

        ! Now do the actual writing
        if (myrank == 0) then
          ! station group
          write(group_name, '(a)') trim(network_names_gather(j, k)) // "." // trim(station_names_gather(j, k))

          call ASDF_open_stations_group_f(waveforms_grp, &
                                          trim(group_name) // C_NULL_CHAR, &
                                          station_grp)

          ! loop over each component
          do i = 1,NDIM
            ! note: if only pressure is output, there is only a single valid component.
            !       the other component_name entries will be empty (""), which will lead to an error
            !       with the HDF5 output when two waveforms with the same name will be defined.
            !       we therefore skip components with an empty ("") entry.
            if (len_trim(component_names_gather(i+((j-1)*NDIM),k)) == 0) cycle

            ! Generate unique waveform name
            ! example: HT.LIT.S3.MXN__2008-01-06T05:14:19__2008-01-06T05:14:53
            write(waveform_name, '(a)') trim(network_names_gather(j,k)) // "." // &
                                        trim(station_names_gather(j,k)) // ".S3." // &
                                        trim(component_names_gather(i+((j-1)*NDIM),k)) //"__"// &
                                        trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"

            call ASDF_open_waveform_f(station_grp, &
                                      trim(waveform_name) // C_NULL_CHAR, &
                                      data_ids(i))

            ! writes (float) data
            call ASDF_write_partial_waveform_f(data_ids(i), &
                                               one_seismogram(i,1:seismo_current_used), 0, seismo_current_used, ier)
            if (ier /= 0) call exit_MPI(myrank,'Error ASDF write partial waveform failed')

            call ASDF_close_dataset_f(data_ids(i), ier)
            if (ier /= 0) call exit_MPI(myrank,'Error ASDF close dataset waveform failed')

          enddo

          call ASDF_close_group_f(station_grp, ier)
          if (ier /= 0) call exit_MPI(myrank,'Error ASDF close stations group waveform failed')

        endif
      enddo
    enddo

    ! closes file
    if (myrank == 0) then
      call ASDF_close_group_f(waveforms_grp, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF close group failed')

      call ASDF_close_file_f(file_id, ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF close file failed')

      call ASDF_finalize_hdf5_f(ier)
      if (ier /= 0) call exit_MPI(myrank,'Error ASDF finalize hdf5 failed')
    endif

  else
    ! write seismograms in parallel
    ! user output
    if (myrank == 0) then
      ! user output
      write(IMAIN,*) 'writing waveforms in parallel...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! initialize
    call ASDF_initialize_hdf5_f(ier);
    if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel initialize hdf5 failed')

    call ASDF_open_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, comm, file_id)
    ! checks file identifier
    if (file_id == 0) then
      print *,'Error rank ',myrank,': ASDF failed to open file',trim(OUTPUT_FILES) // "synthetic.h5"
      call exit_MPI(myrank,'Error opening ASDF file for synthetics')
    endif

    call ASDF_open_waveforms_group_f(file_id, waveforms_grp)

    do k = 1, mysize ! Need to open ASDF groups on all processes

      current_proc = k - 1

      do j = 1, num_stations_gather(k) ! loop over number of stations on that process
        ! station group
        write(group_name, '(a)') trim(network_names_gather(j, k)) // "." // trim(station_names_gather(j, k))

        call ASDF_open_stations_group_f(waveforms_grp, &
                                        trim(group_name) // C_NULL_CHAR, &
                                        station_grp)

        ! Index of current receiver in asdf_container%records
        l = (j-1)*(NDIM)

        ! loop over each component
        do i = 1,NDIM
          ! note: if only pressure is output, there is only a single valid component.
          !       the other component_name entries will be empty (""), which will lead to an error
          !       with the HDF5 output when two waveforms with the same name will be defined.
          !       we therefore skip components with an empty ("") entry.
          if (len_trim(component_names_gather(i+((j-1)*NDIM),k)) == 0) cycle

          ! Generate unique waveform name
          ! example: HT.LIT.S3.MXN__2008-01-06T05:14:19__2008-01-06T05:14:53
          write(waveform_name, '(a)') trim(network_names_gather(j,k)) // "." // &
                                      trim(station_names_gather(j,k)) // ".S3." // &
                                      trim(component_names_gather(i+((j-1)*NDIM),k)) //"__"// &
                                      trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"

          call ASDF_open_waveform_f(station_grp, &
                                    trim(waveform_name) // C_NULL_CHAR, &
                                    data_ids(i))

          if (myrank == current_proc) then
            one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current_used)

            ! writes (float) data
            call ASDF_write_partial_waveform_f(data_ids(i), &
                                               one_seismogram(i,1:seismo_current_used), 0, seismo_current_used, ier)
            if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel write partial waveform failed')
          endif

          call ASDF_close_dataset_f(data_ids(i), ier)
          if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel close dataset waveform failed')

        enddo

        call ASDF_close_group_f(station_grp, ier)
        if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel close group waveform failed')

      enddo
    enddo

    ! makes sure all are finished writing
    call synchronize_all()

    ! closes file
    call ASDF_close_group_f(waveforms_grp, ier)
    if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel close group failed')

    call ASDF_close_file_f(file_id, ier)
    if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel close file failed')

    call ASDF_finalize_hdf5_f(ier)
    if (ier /= 0) call exit_MPI(myrank,'Error ASDF parallel finalize hdf5 failed')

  endif ! WRITE_SEISMOGRAMS_BY_MAIN

  !--------------------------------------------------------
  ! Clean up
  !--------------------------------------------------------
  if (ASDF_OUTPUT_PROVENANCE) deallocate(provenance)
  deallocate(station_names_gather)
  deallocate(network_names_gather)
  deallocate(component_names_gather)
  deallocate(station_lats_gather)
  deallocate(station_longs_gather)
  deallocate(station_elevs_gather)
  deallocate(station_depths_gather)
  deallocate(num_stations_gather)
  deallocate(one_seismogram)

  ! all done
  call synchronize_all()

  end subroutine write_asdf

!
!-------------------------------------------------------------------------------------------------
!

!> Converts the CMT source file read by SPECFEM to a QuakeML file for the ASDF container
!! \param quakemlstring The QuakeML string to store in the ASDF file
!! \start_time_string The start date stored as a character string
  subroutine cmt_to_quakeml(quakemlstring, pde_start_time_string, cmt_start_time_string)

  use constants, only: ASDF_MAX_QUAKEML_LENGTH,ASDF_MAX_TIME_STRING_LENGTH

  use specfem_par, only: &
    hdur, &
    Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    USE_FORCE_POINT_SOURCE,utm_x_source,utm_y_source

  implicit none
  character(len=ASDF_MAX_QUAKEML_LENGTH) :: quakemlstring
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: pde_start_time_string
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: cmt_start_time_string

  ! local parameters
  character(len=13) :: cmt_lon_str, cmt_lat_str, cmt_depth_str, hdur_str
  character(len=13) :: pde_lat_str, pde_lon_str, pde_depth_str
  character(len=25) :: M0_str, mb_str, ms_str, Mw_str
  character(len=25) :: Mrr_str, Mtt_str, Mpp_str, Mrt_str, Mrp_str, Mtp_str
  character(len=20) :: event_name,type_name_str,type_str

  double precision :: M0,Mw,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  double precision :: cmt_lat,cmt_lon,cmt_depth
  double precision, external :: get_cmt_scalar_moment,get_cmt_moment_magnitude

  ! setup
  if (USE_FORCE_POINT_SOURCE) then
    ! force
    !(will need to find better infos, for now just artifical setup)
    !
    ! QuakeML:
    ! event types allowed: earthquake, explosion, crash, other event, thunder, slide, meteorite, ..
    ! event description types allowed: earthquake name, felt report, region name, ..
    !
    ! see: https://quake.ethz.ch/quakeml/docs/latest?action=AttachFile&do=get&target=QuakeML-BED.pdf
    !
    ! will use "controlled explosion" for force point source, sounds more exciting than "other event";
    ! and "earthquake name" for description, however, in principle it would be something else.
    type_str = "controlled explosion"  ! "force" doesn't work, unfortunately, "hammer" neither...
    type_name_str = "earthquake name"  ! "force name" doesn't work, "active source name" neither...
    event_name = "forcesolution1"
    cmt_lat = 0.d0
    cmt_lon = 0.d0
    cmt_depth = 0.d0
    Mrr = 1.d0
    Mtt = 0.d0
    Mpp = 0.d0
    Mrt = 0.d0
    Mrp = 0.d0
    Mtp = 0.d0
    M0 = 1.d0
    Mw = 1.d0
  else
    ! moment tensor
    type_str = "earthquake"
    type_name_str = "earthquake name"
    event_name = "cmtsolution1"
    ! coordinate system convention: x points east, y points north, z points vertical up
    cmt_lat = utm_y_source(1)
    cmt_lon = utm_x_source(1)
    cmt_depth = 0.d0    ! no info yet
    ! coordinate system convention: r -> z, theta -> -y, phi -> x
    !   Mrr =  Mzz
    !   Mtt =  Myy
    !   Mpp =  Mxx
    !   Mrt = -Myz
    !   Mrp =  Mxz
    !   Mtp = -Mxy
    Mrr = Myy(1)  !Mrr*1e-7
    Mtt = Myy(1)  !Mtt*1e-7
    Mpp = Mxx(1)  !Mpp*1e-7
    Mrt = -Myz(1) !Mrt*1e-7
    Mrp = Mxz(1)  !Mrp*1e-7
    Mtp = -Mxy(1) !Mtp*1e-7
    M0 = get_cmt_scalar_moment(Mxx(1),Myy(1),Mzz(1),Mxy(1),Mxz(1),Myz(1)) ! dyne-cm
    M0 = M0 * 1e-7 ! dyn-cm to N-m conversion
    Mw = get_cmt_moment_magnitude(Mxx(1),Myy(1),Mzz(1),Mxy(1),Mxz(1),Myz(1))
  endif

  ! Convert the CMT values to strings for the QuakeML string
  pde_start_time_string = cmt_start_time_string

  ! no infos yet stored for depth
  write(pde_lat_str, "(g12.5)") cmt_lat     !pde_lat
  write(pde_lon_str, "(g12.5)") cmt_lon     !pde_lon
  write(pde_depth_str, "(g12.5)") cmt_depth !pde_depth*1000 ! km to m conversion
  write(cmt_lat_str, "(g12.5)") cmt_lat     !cmt_lat
  write(cmt_lon_str, "(g12.5)") cmt_lon     !cmt_lon
  write(cmt_depth_str, "(g12.5)") cmt_depth !cmt_depth*1000 ! km to m conversion

  ! takes first source description
  write(hdur_str, "(g12.5)") hdur(1)
  write(M0_str, "(g12.5)") M0   !M0
  write(mb_str, "(g12.5)") 0.0  !mb
  write(ms_str, "(g12.5)") 0.0  !ms
  write(Mw_str, "(g12.5)") Mw   !Mw

  write(Mrr_str, "(g12.5)") Mrr
  write(Mtt_str, "(g12.5)") Mtt
  write(Mpp_str, "(g12.5)") Mpp
  write(Mrt_str, "(g12.5)") Mrt
  write(Mrp_str, "(g12.5)") Mrp
  write(Mtp_str, "(g12.5)") Mtp

  ! header version
  quakemlstring = &
    '<q:quakeml xmlns="http://quakeml.org/xmlns/bed/1.2" xmlns:q="http://quakeml.org/xmlns/quakeml/1.2">'

  ! event
  quakemlstring = trim(quakemlstring) // &
    '<eventParameters publicID="smi:local/'//trim(event_name)//'#eventPrm">'// &
    '<event publicID="smi:local/'//trim(event_name)//'#eventID">'

  ! info preferences
  quakemlstring = trim(quakemlstring) // &
    '<preferredOriginID>smi:local/'//trim(event_name)//'/origin#cmtorigin</preferredOriginID>'// &
    '<preferredMagnitudeID>smi:local/'//trim(event_name)//'/magnitude#moment_mag</preferredMagnitudeID>'// &
    '<preferredFocalMechanismID>smi:local/'//trim(event_name)//'/focal_mechanism</preferredFocalMechanismID>'

  ! event name
  quakemlstring = trim(quakemlstring) // &
    '<type>'//trim(type_str)//'</type>'// &
    '<typeCertainty>known</typeCertainty>'// &
    '<description>'// &
    '  <text>'//trim(event_name)//'</text>'// &
    '  <type>'//trim(type_name_str)//'</type>'// &
    '</description>'

  ! event origin
  quakemlstring = trim(quakemlstring) // &
    '<origin publicID="smi:local/'//trim(event_name)//'/origin#reforigin">'//&
    '  <time>'//&
    '    <value>'//trim(pde_start_time_string)//'</value>'//&
    '  </time>'//&
    '  <latitude>'//&
    '    <value>'//trim(pde_lat_str)//'</value>'//&
    '  </latitude>'//&
    '  <longitude>'//&
    '    <value>'//trim(pde_lon_str)//'</value>'//&
    '  </longitude>'//&
    '  <depth>'//&
    '    <value>'//trim(pde_depth_str)//'</value>'//&
    '  </depth>'//&
    '  <type>hypocenter</type>'//&
    '  <comment id="smi:local/'//trim(event_name)//'/comment#ref_origin">'//&
    '    <text>Hypocenter catalog: PDE</text>'//&
    '  </comment>'//&
    '</origin>'

  ! event cmtorigin
  quakemlstring = trim(quakemlstring) // &
    '<origin publicID="smi:local/'//trim(event_name)//'/origin#cmtorigin">'//&
    '  <time>'//&
    '    <value>'//trim(cmt_start_time_string)//'</value>'//&
    '  </time>'// &
    '  <latitude>'//&
    '    <value>'//trim(cmt_lat_str)//'</value>'//&
    '  </latitude>'//&
    '  <longitude>'//&
    '    <value>'//trim(cmt_lon_str)//'</value>'//&
    '  </longitude>'//&
    '  <depth>'//&
    '    <value>'//trim(cmt_depth_str)//'</value>'//&
    '  </depth>'//&
    '</origin>'

  ! event focal mechanism (empty for now...)
  quakemlstring = trim(quakemlstring) // &
    '<focalMechanism publicID="smi:local/'//trim(event_name)//'/focal_mechanism">'//&
    '<momentTensor publicID="smi:local/'//trim(event_name)//'/momenttensor">'//&
    '  <derivedOriginID>smi:local/'//trim(event_name)//'/origin#cmtorigin'//'</derivedOriginID>'//&
    '  <momentMagnitudeID>smi:local/'//trim(event_name)//'/magnitude#moment_mag'//'</momentMagnitudeID>'//&
    '  <scalarMoment>'//&
    '    <value>'//trim(M0_str)//'</value>'//&
    '  </scalarMoment>'//&
    '  <tensor>'//&
    '  <Mrr>'//&
    '    <value>'//trim(Mrr_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mrr>'//&
    '  <Mtt>'//&
    '    <value>'//trim(Mtt_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mtt>'//&
    '  <Mpp>'//&
    '    <value>'//trim(Mpp_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mpp>'//&
    '  <Mrt>'//&
    '    <value>'//trim(Mrt_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mrt>'//&
    '  <Mrp>'//&
    '    <value>'//trim(Mrp_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mrp>'//&
    '  <Mtp>'//&
    '    <value>'//trim(Mtp_str)//'</value>'//&
    '    <uncertainty>0</uncertainty>'//&
    '  </Mtp>'//&
    '  </tensor>'//&
    '  <sourceTimeFunction>'//&
    '    <type>triangle</type>'//&
    '    <duration>'//trim(hdur_str)//'</duration>'//&
    '  </sourceTimeFunction>'//&
    '</momentTensor>'//&
    '</focalMechanism>'

  ! event magnitudes
  quakemlstring = trim(quakemlstring) // &
    '<magnitude publicID="smi:local/'//trim(event_name)//'/magnitude#moment_mag">'//&
    '  <mag>'//&
    '    <value>'//trim(Mw_str)//'</value>'//&
    '  </mag>'//&
    '  <type>Mwc</type>'//&
    '</magnitude>'//&
    '<magnitude publicID="smi:local/'//trim(event_name)//'/magnitude#mb">'//&
    '  <mag>'//&
    '    <value>'//trim(mb_str)//'</value>'//&
    '  </mag>'//&
    '  <type>mb</type>'//&
    '</magnitude>'//&
    '<magnitude publicID="smi:local/'//trim(event_name)//'/magnitude#MS">'//&
    '  <mag>'//&
    '    <value>'//trim(ms_str)//'</value>'//&
    '  </mag>'//&
    '  <type>MS</type>'//&
    '</magnitude>'

  ! event finish
  quakemlstring = trim(quakemlstring) // &
    '</event>'//&
    '</eventParameters>'//&
    '</q:quakeml>'

  !debug
  !print *,'debug: quakeML:'
  !print *,'*****'//trim(quakemlstring)//'*****'

  end subroutine cmt_to_quakeml

!
!-------------------------------------------------------------------------------------------------
!

! Convert system time into a string
!! \param time system time
!! \param time_string a string representation of the input time
  subroutine convert_systime_to_string(time, time_string)

  use constants, only: ASDF_MAX_TIME_STRING_LENGTH

  implicit none

  double precision, intent(in) :: time
  character(len=ASDF_MAX_TIME_STRING_LENGTH), intent(out) :: time_string

  ! local parameters
  real :: fraction_sec
  integer :: iatime(9)
  character(len=4) :: yr
  character(len=2) :: mo, da, hr, minute
  character(len=15) :: str_second
  real :: real_sec
  integer :: stime,iyr,imo,ida,ihr,imin,isec,imillisec

  ! extract msec
  fraction_sec = time - int(time)

  ! function returns UTC time
  !call gmtime(int(time), iatime)
  !
  ! note: gmtime is not a standard Fortran function, but an extension.
  !       thus, it's implementation can differ from one compiler to another.
  !       we see problems with gmtime() on IBM xlf compilers.
  !       however, the C/C++ gmtime() function is C99 standard, thus we call here a wrapper function in param_reader.c
  !
  stime = int(time)
  call get_utctime_params(stime,iyr,imo,ida,ihr,imin,isec)

  ! see e.g.: https://gcc.gnu.org/onlinedocs/gcc-5.5.0/gfortran/GMTIME.html
  iatime(1) = isec
  iatime(2) = imin
  iatime(3) = ihr
  iatime(4) = ida
  iatime(5) = imo
  iatime(6) = iyr

  write(yr, "(I4.4)") iatime(6) + 1900
  write(mo, "(I2.2)") iatime(5) + 1
  write(da, "(I2.2)") iatime(4)
  write(hr, "(I2.2)") iatime(3)
  write(minute, "(I2.2)") iatime(2)

  real_sec = iatime(1) + fraction_sec
  imillisec = int(real_sec * 10**3 - int(real_sec) * 10**3)

  ! avoid zero-length specifier, leads to problems with different compilers
  !!write(str_second, "(I2.2,F0.4)") int(real_sec), real_sec-int(real_sec)
  ! for example: 2.891455 -> string 02.891
  write(str_second, "(I2.2,'.',I3.3)") int(real_sec), imillisec

  ! format example: 2018-01-31T16:40:02.89  - since ASDF_MAX_TIME_STRING_LENGTH is 22 characters
  time_string = trim(yr)//"-"//trim(mo)//"-"//trim(da)//"T"//&
                  trim(hr)//':'//trim(minute)//':'//trim(str_second)

  !debug
  !print *,'debug: time_string:',real_sec,int(real_sec),real_sec-int(real_sec),imillisec, &
  !        'str_second***',str_second,'***',trim(str_second),'***'
  !print *,'*****'//trim(time_string)//'*****'

  end subroutine convert_systime_to_string

!
!-------------------------------------------------------------------------------------------------
!

!> Uses the time in the CMTSOLUTION file to calculate the number of seconds since the epoch
!! \param starttime The start time of the simulation from the epoch
!! \param start_time_string A string for defining the waveform name start time
!! \param pde_start_time_string A string for defining the waveform name start time using PDE
!! \param cmt_start_time_string A string for defining the waveform name start time using CMT
!! \param end_time_string A string for defining the waveform name end time
  subroutine get_time_cmt(starttime, start_time_string, pde_start_time_string, cmt_start_time_string, end_time_string)

  use constants, only: ASDF_MAX_TIME_STRING_LENGTH

  use specfem_par, only: &
    NSTEP,DT,t0, &
    min_tshift_src_original, &
    yr_PDE,jda_PDE,ho_PDE,mi_PDE,sec_PDE !,hdur

  use specfem_par, only: USE_FORCE_POINT_SOURCE

  implicit none
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: start_time_string
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: end_time_string
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: cmt_start_time_string
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: pde_start_time_string
  double precision, intent(inout) :: starttime

  ! local parameters
  double precision :: trace_length_in_sec
  integer :: yr,year,jda,ho,mi
  double precision :: sec
  double precision :: pdetime, cmttime, endtime
  integer,dimension(8) :: values

  ! PDE time
  if (USE_FORCE_POINT_SOURCE) then
    ! force solution
    ! sets an artificial PDE time for SPECFEM seismograms:
    !   year 2000, January 1, time at 01:00:00 UTC (coordinated universal time)
    yr = 2000
    jda = 0
    ho = 1
    mi = 0
    sec = 0.d0
  else
    ! CMT solution
    ! takes PDE time from CMTSOLUTION header line
    yr = yr_PDE
    jda = jda_PDE
    ho = ho_PDE
    mi = mi_PDE
    sec = sec_PDE
  endif

  ! Calculates the start time since the epoch in seconds
  ! Reference:
  ! http://pubs.opengroup.org/onlinepubs/009695399/basedefs/xbd_chap04.html#tag_04_14
  year = yr - 1900
  pdetime = (year-70)*31536000.0d0 &
          + ((year-69)/4)*86400.0d0 &
          - ((year-1)/100)*86400.0d0 &
          + ((year+299)/400)*86400.0d0 &
          + (jda-1)*86400.0d0 + (ho)*3600.0d0 + (mi)*60.0d0 + sec
  call convert_systime_to_string(pdetime, pde_start_time_string)

  ! CMT centroid time
  ! (time shift as specified by description in CMTSOLUTION file; PDE time is in header line, CMT time can be slightly shifted)
  cmttime = pdetime + min_tshift_src_original
  call convert_systime_to_string(cmttime, cmt_start_time_string)

  ! trace start time
  ! (CMT time shifted by the onset time t0 of the simulation)
  starttime = cmttime - t0
  call convert_systime_to_string(starttime, start_time_string)

  ! Calculates the number of seconds to add to the start_time
  trace_length_in_sec = DT*NSTEP
  endtime = starttime + trace_length_in_sec
  call convert_systime_to_string(endtime, end_time_string)

  ! Calculates time in seconds since the epoch
  call date_and_time(VALUES=values)
  !time8()

  !write(*,*) fmtdate(values,'The CPU time used by this program is now %c seconds')

  end subroutine get_time_cmt

!
!-------------------------------------------------------------------------------------------------
!

!> Converts the Station information to a StationXML string
!! \param station_name The name of the station
!! \param network_name The name of the network
!! \param stationxmlstring The StationXML string that will be written to the ASDF file
  subroutine station_to_stationxml(station_name, network_name, latitude, longitude, elevation, &
                                   burial_depth, start_time_string, stationxmlstring)

  use constants, only: ASDF_MAX_STATIONXML_LENGTH,ASDF_MAX_TIME_STRING_LENGTH, &
    MAX_LENGTH_NETWORK_NAME,MAX_LENGTH_STATION_NAME

  implicit none
  character(len=ASDF_MAX_STATIONXML_LENGTH) :: stationxmlstring
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=ASDF_MAX_TIME_STRING_LENGTH) :: start_time_string

  character(len=15) :: station_lat, station_lon, station_ele, station_depth
  integer :: len_station_name, len_network_name, len_station_lat
  integer :: len_station_lon, len_station_depth, len_station_ele
  real, intent(in) :: latitude, longitude, elevation, burial_depth

  ! Convert double precision to character strings for the StationXML string

  write(station_lat, "(F11.4)") latitude
  write(station_lon, "(F11.4)") longitude
  write(station_ele, "(F8.1)") elevation
  write(station_depth, "(F8.1)") burial_depth

  ! print *, trim(station_lat), trim(station_lon), trim(station_depth), trim(station_ele)

  len_network_name = len(network_name)
  len_station_name = len(station_name)
  len_station_lat = len(trim(station_lat))
  len_station_lon = len(trim(station_lon))
  len_station_depth = len(trim(station_depth))
  len_station_ele = len(trim(station_ele))

  ! header
  stationxmlstring = &
    '<FDSNStationXML schemaVersion="1.0" xmlns="http://www.fdsn.org/xml/station/1">'

  ! info
  stationxmlstring = trim(stationxmlstring) // &
    '<Source>SPECFEM3D_Cartesian</Source>'//&
    '<Module>SPECFEM3D_Cartesian/asdf-library</Module>'//&
    '<ModuleURI>http://seismic-data.org</ModuleURI>'//&
    '<Created>'//trim(start_time_string)//'</Created>'

  ! station
  stationxmlstring = trim(stationxmlstring) // &
    '<Network code="'//trim(network_name(1:len(network_name)))//'">'//&
    '<Station code="'//trim(station_name(1:len(station_name)))//'">'//&
    '  <Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
    '  <Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
    '  <Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
    '  <Site>'//&
    '    <Name>N/A</Name>'//&
    '  </Site>'//&
    '  <CreationDate>'//trim(start_time_string)//'</CreationDate>'

  ! channels
  stationxmlstring = trim(stationxmlstring) // &
    '  <TotalNumberChannels>3</TotalNumberChannels>'//&
    '  <SelectedNumberChannels>3</SelectedNumberChannels>'//&
    '  <Channel locationCode="S3" code="MXN" startDate="'//trim(start_time_string)//'">'//&
    '    <Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
    '    <Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
    '    <Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
    '    <Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
    '    <Azimuth>0.0</Azimuth>'//&
    '    <Dip>0.0</Dip>'//&
    '  </Channel>'//&
    '  <Channel locationCode="S3" code="MXE" startDate="'//trim(start_time_string)//'">'//&
    '    <Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
    '    <Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
    '    <Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
    '    <Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
    '    <Azimuth>90.0</Azimuth>'//&
    '    <Dip>0.0</Dip>'//&
    '  </Channel>'//&
    '  <Channel locationCode="S3" code="MXZ" startDate="'//trim(start_time_string)//'">'//&
    '    <Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
    '    <Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
    '    <Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
    '    <Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
    '    <Azimuth>0.0</Azimuth>'//&
    '    <Dip>90.0</Dip>'//&
    '  </Channel>'

  ! finish
  stationxmlstring = trim(stationxmlstring) // &
    '</Station>'//&
    '</Network>'//&
    '</FDSNStationXML>'

  !debug
  !print *,'debug: station XML:'
  !print *,'*****'//trim(stationxmlstring)//'*****'

  end subroutine station_to_stationxml

!
!-------------------------------------------------------------------------------------------------
!

!> Reads an external file and stores it in filestring
!! \param filename The name of the file to read
!! \param filestring The string that the file is stored
  subroutine read_file(filename, filestring, filesize)

  use constants, only: IIN,myrank

  implicit none
  character(len=*) :: filestring
  character(len=*) :: filename
  integer,intent(out) :: filesize

  ! local parameters
  integer :: ier

  ! Get the size of the file using Fortran2003 feature
  open(unit=IIN, file=trim(filename), status='old', iostat=ier)
  if (ier /= 0) then
    print *,'Error: rank ',myrank,'failed opening file: ',trim(filename)
    call exit_MPI(myrank,'Error opening file in ASDF read_file() routine')
  endif
  inquire(unit=IIN, size=filesize)
  close(IIN)

  ! checks size
  if (len(filestring) < filesize) then
    print *,'Error: file ',trim(filename),' exceeds filesize ',filesize,'for string length limit ',len(filestring)
    call exit_MPI(myrank,'Error file size too large in ASDF read_file() routine')
  endif

  ! Read in the size of the file using direct access
  open(IIN, file=trim(filename), status='old', &
       recl=filesize, form='unformatted', access='direct')
  read(IIN, rec=1) filestring(1:filesize)
  close(IIN)

  end subroutine read_file
