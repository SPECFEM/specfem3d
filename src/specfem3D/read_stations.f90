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


  subroutine read_stations(rec_filename)

  use constants, only: MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,IIN,MAX_STRING_LEN

  use specfem_par, only: myrank,nrec,station_name,network_name,stlat,stlon,stele,stbur

  implicit none

  ! input receiver file name
  character(len=*),intent(in) :: rec_filename

  ! local parameters
  integer :: irec,ier
  character(len=MAX_STRING_LEN) :: line

  ! loop on all the stations to read the file
  if (myrank == 0) then

    ! opens STATIONS or STATIONS_ADJOINT file
    open(unit=IIN,file=trim(rec_filename),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(rec_filename))

    ! reads all stations
    do irec = 1,nrec
      read(IIN,"(a)",iostat=ier) line
      if (ier /= 0) call exit_mpi(myrank, 'Error reading station file '//trim(rec_filename))

      ! skip comment lines
      if (len_trim(line) > 0 .and. line(1:1) /= '#') then
        line = trim(line)
        read(line,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
        if (ier /= 0) call exit_mpi(myrank, 'Error reading station file line '//trim(line))
      endif
    enddo

    ! close receiver file
    close(IIN)

    ! find duplicate station names
    call read_stations_find_duplets()
  endif

  ! broadcast values to other slices
  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)

  end subroutine read_stations

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_stations_find_duplets()

  use constants, only: MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
                       DUPLETS_NREC_MINIMUM_FOR_HASH

  use specfem_par, only: myrank,nrec,station_name,network_name

  implicit none

  ! local parameters
  integer :: i,irec,ier

  integer, allocatable, dimension(:) :: station_duplet
  integer :: length_station_name,length_network_name

  ! hash table
  integer, allocatable, dimension(:) :: hash_table
  integer :: hash,hash_prob,hash_collisions

  ! In case that the same station and network name appear twice (or more times) in the STATIONS
  ! file, problems occur, as two (or more) seismograms are written (with mode
  ! "append") to a file with same name. The philosophy here is to accept multiple
  ! appearances and to just add a count to the station name in this case.
  allocate(station_duplet(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1970')
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating station_duplet array')
  station_duplet(:) = 0

  ! initializes the hash table
  if (nrec >= DUPLETS_NREC_MINIMUM_FOR_HASH) then
    ! makes hash table slightly larger than the actual number of stations to avoid many hash collisions
    allocate(hash_table(5*nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1971')
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating hash_table array')
    hash_table(:) = 0
    hash_collisions = 0
  endif

  do irec = 1,nrec
    ! checks if duplicate of another station read in so far
    if (nrec < DUPLETS_NREC_MINIMUM_FOR_HASH) then
      ! way 1:
      ! loops over all previous stations and checks station/network names
      ! this will get slow for very large STATIONS files (> 100,000 stations); scales approx. quadratic ~O(n**2)
      do i = 1,irec-1
        if ((station_name(irec) == station_name(i)) .and. (network_name(irec) == network_name(i))) then
          ! increases duplet count
          station_duplet(i) = station_duplet(i) + 1
          ! appends duplet number to station name
          if (len_trim(station_name(irec)) <= MAX_LENGTH_STATION_NAME-3) then
            write(station_name(irec),"(a,'_',i2.2)") trim(station_name(irec)),station_duplet(i)+1
          else
            call exit_MPI(myrank,'Please increase MAX_LENGTH_STATION_NAME by at least 3 to name station duplets')
          endif
        endif
      enddo

    else
      ! way 2:
      ! gets a hash value and checks in hash table for duplicates; scales approx. linearly ~O(n)
      hash = hashFunc(station_name(irec),network_name(irec))
      if (hash_table(hash) == 0) then
        ! stores station index
        hash_table(hash) = irec
      else
        ! found a duplicate hash
        ! check if name matches
        i = hash_table(hash)
        if ((station_name(irec) == station_name(i)) .and. (network_name(irec) == network_name(i))) then
          ! debug
          !print *,'debug: Duplicate found in hash table:', &
          !        irec,trim(station_name(irec)),trim(network_name(irec)), &
          !        ' - hash number',hash, &
          !        'return index',i,trim(station_name(i)),trim(network_name(i))
          ! increases duplet count
          station_duplet(i) = station_duplet(i) + 1
          ! appends duplet number to station name
          if (len_trim(station_name(irec)) <= MAX_LENGTH_STATION_NAME-3) then
            write(station_name(irec),"(a,'_',i2.2)") trim(station_name(irec)),station_duplet(i)+1
          else
            call exit_MPI(myrank,'Please increase MAX_LENGTH_STATION_NAME by at least 3 to name station duplets')
          endif
        else
          ! hash collision
          hash_collisions = hash_collisions + 1
          ! debug
          !print *,'debug: Collision found in hash table:', &
          !        irec,trim(station_name(irec)),trim(network_name(irec)), &
          !        ' - hash number',hash, &
          !        'return index',i,trim(station_name(i)),trim(network_name(i))
          ! put hash in next free slot (linear probing)
          hash_prob = hash
          do while (hash_table(hash_prob) /= 0)
            ! increases hash index
            hash_prob = mod(hash_prob + 1,size(hash_table))
            ! check if we reach again same hash, then table is full
            if (hash_prob == hash) then
              print *,'Error: Hash table is full, please consider a larger hash table!'
              call exit_MPI(myrank,'Please increase hash table size for station duplets search')
            endif
            ! check entry
            i = hash_table(hash_prob)
            if (i == 0) then
              ! stores station index in new free slot
              hash_table(hash_prob) = irec
              exit
            else if (i == irec) then
              ! station already stored, done adding
              exit
            else
              ! check entry, maybe station moved hash index due to previous collisions
              if ((station_name(irec) == station_name(i)) .and. (network_name(irec) == network_name(i))) then
                ! debug
                !print *,'debug: Duplicate found in hash table:', &
                !        irec,trim(station_name(irec)),trim(network_name(irec)), &
                !        ' - hash number by collision probing',hash_prob, &
                !        'return index',i,trim(station_name(i)),trim(network_name(i))
                ! increases duplet count
                station_duplet(i) = station_duplet(i) + 1
                ! appends duplet number to station name
                if (len_trim(station_name(irec)) <= MAX_LENGTH_STATION_NAME-3) then
                  write(station_name(irec),"(a,'_',i2.2)") trim(station_name(irec)),station_duplet(i)+1
                else
                  call exit_MPI(myrank,'Please increase MAX_LENGTH_STATION_NAME by at least 3 to name station duplets')
                endif
                ! done
                exit
              endif
            endif
          enddo
        endif

      endif
    endif

    ! checks name lengths
    length_station_name = len_trim(station_name(irec))
    length_network_name = len_trim(network_name(irec))

    ! check that length conforms to standard
    if (length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) then
      print *, 'Error: invalid station name ',trim(station_name(irec))
      call exit_MPI(myrank,'wrong length of station name')
    endif
    if (length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) then
      print *, 'Error: invalid network name ',trim(network_name(irec))
      call exit_MPI(myrank,'wrong length of network name')
    endif

  enddo

  ! user output
  if (sum(station_duplet) > 0) then
    print *
    print *,'Warning: found ',sum(station_duplet),' station duplets (having same network & station names)'
    do irec = 1,nrec
      if (station_duplet(irec) > 0) then
        print *,'  station_name  : ',station_name(irec),' network: ',network_name(irec)
        print *,'           irec : ',irec,' duplets: ',station_duplet(irec)
      endif
    enddo
    print *,'Please check your STATIONS file entries to avoid confusions...'
    print *
  endif

  ! debug: hash table info
  !if (nrec >= DUPLETS_NREC_MINIMUM_FOR_HASH .and. hash_collisions > 0) &
  !  print *,'debug: hash table collisions: ',hash_collisions

  ! free temporary array
  deallocate(station_duplet)
  if (nrec >= DUPLETS_NREC_MINIMUM_FOR_HASH) deallocate(hash_table)

  contains

    ! defines a simple hash function
    integer function hashFunc(sta_name,net_name)
      character(len=MAX_LENGTH_STATION_NAME), intent(in) :: sta_name
      character(len=MAX_LENGTH_NETWORK_NAME), intent(in) :: net_name
      ! local parameters
      integer :: i, sum, prime
      integer, parameter :: base = 31   ! seems to be a good number (arbitrary choice for a prime)
                                        ! smaller prime numbers lead to more collisions
      character(len=(MAX_LENGTH_STATION_NAME+MAX_LENGTH_NETWORK_NAME)) :: name

      ! full name, e.g., S00001DB
      name = sta_name // net_name

      sum = 0
      prime = 1
      do i = 1, len_trim(name)
        ! too simple, will get lots of hash collisions
        ! sum = sum + ichar(name(i:i))

        ! adding a bit more complexity, based on polynomial rolling
        sum = mod(sum + ichar(name(i:i)) * prime, size(hash_table))
        prime = mod(prime * base, size(hash_table))
      enddo

      hashFunc = mod(sum, size(hash_table))
    end function hashFunc


  end subroutine read_stations_find_duplets

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_stations_from_previous_run(is_done_stations)

  use constants, only: NDIM,IOUT_SU,IMAIN,OUTPUT_FILES,MAX_STRING_LEN

  use specfem_par, only: myrank, &
    nrec,station_name,network_name, &
    islice_selected_rec,ispec_selected_rec, &
    xi_receiver,eta_receiver,gamma_receiver, &
    nu_rec, stations_hashsum

  implicit none

  logical, intent(out) :: is_done_stations

  ! local parameters
  integer :: ier,irec

  ! file parameters
  logical :: station_file_exists

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  character(len=MAX_STRING_LEN) :: filename
  character(len=32) :: hashsum_read_from_file

  ! initializes
  is_done_stations = .false.

  ! only check for existing file if number of stations is big enough
  if (nrec < 1000) return

  filename = trim(OUTPUT_FILES)//'/stations_info.bin'

  ! checks if file with station infos located from previous run exists
  inquire(file=trim(filename),exist=station_file_exists)

  ! checks if something to do
  if (.not. station_file_exists) return

  ! main reads in available station information
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) 'reading station details:'
    write(IMAIN,*) '  from file: stations_info.bin'
    call flush_IMAIN()

    open(unit=IOUT_SU,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'error opening file '//trim(filename))

    ! temporary arrays
    allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1971')

    ! reads hash key
    read(IOUT_SU) hashsum_read_from_file

    ! user output
    write(IMAIN,*) '  hash key : ',hashsum_read_from_file
    write(IMAIN,*)
    call flush_IMAIN()

    ! get a hash checksum from STATIONS infos
    if (len_trim(stations_hashsum) == 0) then
      call create_stations_checksum(stations_hashsum)
    endif

    !debug
    !print *,'debug: hash from STATIONS = ',stations_hashsum
    !print *,'debug: hash from file     = ',hashsum_read_from_file,' match: ',(stations_hashsum==hashsum_read_from_file)

    ! checks hash keys
    if (stations_hashsum /= hashsum_read_from_file) then
      write(IMAIN,*) '  hash checksum does not match:'
      write(IMAIN,*) '    from STATIONS          file: ',stations_hashsum
      write(IMAIN,*) '    from stations_info.bin file: ',hashsum_read_from_file
      write(IMAIN,*) '  will fall back to locate stations again...'

      ! closes file
      close(IOUT_SU)

      ! return w/out infos
      is_done_stations = .false.
      return
    endif

    ! reads in station infos
    read(IOUT_SU) islice_selected_rec,ispec_selected_rec
    read(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
    read(IOUT_SU) x_found,y_found,z_found
    read(IOUT_SU) nu_rec
    close(IOUT_SU)

    ! write the locations of stations, so that we can load them and write them to SU headers later
    open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
         status='unknown',action='write',iostat=ier)
    if (ier /= 0) &
      call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'/output_list_stations.txt')

    do irec = 1,nrec
      write(IOUT_SU,*) station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
    enddo

    close(IOUT_SU)
    deallocate(x_found,y_found,z_found)
  endif

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)
  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)
  call bcast_all_dp(nu_rec,NDIM*NDIM*nrec)
  call synchronize_all()

  ! everything done
  is_done_stations = .true.

  end subroutine read_stations_from_previous_run

!
!-------------------------------------------------------------------------------------------
!

  subroutine write_stations_for_next_run(x_found,y_found,z_found)


  use constants, only: myrank,IMAIN,IOUT_SU,OUTPUT_FILES,MAX_STRING_LEN

  use specfem_par, only: nrec,islice_selected_rec,ispec_selected_rec, &
                         xi_receiver,eta_receiver,gamma_receiver,nu_rec, &
                         stations_hashsum
  implicit none

  double precision, dimension(nrec),intent(in) :: x_found,y_found,z_found

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  ! only output file if number of stations is big enough
  if (nrec < 1000) return

  ! only main process writes out all available station information
  if (myrank /= 0) return

  ! user output
  write(IMAIN,*) 'storing station details:'
  write(IMAIN,*) '  to file  : stations_info.bin'
  call flush_IMAIN()

  ! get checksum from STATIONS file infos
  ! could have been computed already in read_stations_from_previous run
  if (len_trim(stations_hashsum) == 0) then
    call create_stations_checksum(stations_hashsum)
  endif

  ! user output
  write(IMAIN,*) '  hash key : ',stations_hashsum
  write(IMAIN,*)
  call flush_IMAIN()

  ! creates file with all station location infos
  filename = trim(OUTPUT_FILES)//'/stations_info.bin'

  open(unit=IOUT_SU,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier == 0) then
    ! hash checksum for checking
    write(IOUT_SU) stations_hashsum
    ! stations informations
    write(IOUT_SU) islice_selected_rec,ispec_selected_rec
    write(IOUT_SU) xi_receiver,eta_receiver,gamma_receiver
    write(IOUT_SU) x_found,y_found,z_found
    write(IOUT_SU) nu_rec
    close(IOUT_SU)
  endif

  end subroutine write_stations_for_next_run

!
!-------------------------------------------------------------------------------------------
!

  subroutine create_stations_checksum(hashsum)

! creates a hash key to check if stations_info.bin file contains a valid number of stations
! compared to STATIONS file infos
!
! the checksum algorithm follows the principles of an MD5 hash algorithm,
! but its resulting hash differs to the exact MD5 one.

  use constants, only: MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME
  use specfem_par, only: nrec,station_name,network_name,stlat,stlon,stbur

  implicit none

  character(len=32),intent(out) :: hashsum         ! 128-bit hash key

  ! md5 parameters
  ! binary integer part of the sines of integers (Radians) as constants
  ! note: the initialization here with boz constants won't work with gfortran
  !integer(kind=4), dimension(64), parameter :: K =(/ &
  !    z'0xD76AA478', z'0xE8C7B756', z'0x242070DB', z'0xC1BDCEEE', &
  !    z'0xF57C0FAF', z'0x4787C62A', z'0xA8304613', z'0xFD469501', &
  !    z'0x698098D8', z'0x8B44F7AF', z'0xFFFF5BB1', z'0x895CD7BE', &
  !    z'0x6B901122', z'0xFD987193', z'0xA679438E', z'0x49B40821', &
  !    z'0xF61E2562', z'0xC040B340', z'0x265E5A51', z'0xE9B6C7AA', &
  !    z'0xD62F105D', z'0x02441453', z'0xD8A1E681', z'0xE7D3FBC8', &
  !    z'0x21E1CDE6', z'0xC33707D6', z'0xF4D50D87', z'0x455A14ED', &
  !    z'0xA9E3E905', z'0xFCEFA3F8', z'0x676F02D9', z'0x8D2A4C8A', &
  !    z'0xFFFA3942', z'0x8771F681', z'0x6D9D6122', z'0xFDE5380C', &
  !    z'0xA4BEEA44', z'0x4BDECFA9', z'0xF6BB4B60', z'0xBEBFBC70', &
  !    z'0x289B7EC6', z'0xEAA127FA', z'0xD4EF3085', z'0x04881D05', &
  !    z'0xD9D4D039', z'0xE6DB99E5', z'0x1FA27CF8', z'0xC4AC5665', &
  !    z'0xF4292244', z'0x432AFF97', z'0xAB9423A7', z'0xFC93A039', &
  !    z'0x655B59C3', z'0x8F0CCC92', z'0xFFEFF47D', z'0x85845DD1', &
  !    z'0x6FA87E4F', z'0xFE2CE6E0', z'0xA3014314', z'0x4E0811A1', &
  !    z'0xF7537E82', z'0xBD3AF235', z'0x2AD7D2BB', z'0xEB86D391' &
  !    /)
  ! or to be computed below
  integer(kind=4), dimension(64) :: K

  ! specifies the per-round shift amounts
  integer(kind=4), dimension(64), parameter :: S = (/ &
      7, 12, 17, 22, 7, 12, 17, 22, &
      7, 12, 17, 22, 7, 12, 17, 22, &
      5,  9, 14, 20, 5,  9, 14, 20, &
      5,  9, 14, 20, 5,  9, 14, 20, &
      4, 11, 16, 23, 4, 11, 16, 23, &
      4, 11, 16, 23, 4, 11, 16, 23, &
      6, 10, 15, 21, 6, 10, 15, 21, &
      6, 10, 15, 21, 6, 10, 15, 21 &
      /)

  integer(kind=4) :: a0, b0, c0, d0
  integer(kind=4) :: a, b, c, d, f, g, temp
  integer :: i, j, pos, str_len, bit_len, size_str4
  integer :: irec ! ier
  ! buffer 64-byte (16xint)
  character(len=8) :: wtmp
  integer(kind=4) :: w(16),i1,i2,i3,i4

  ! character string holding sta/net/x/y/z
  character(len=MAX_LENGTH_STATION_NAME + MAX_LENGTH_NETWORK_NAME + 3*12) :: single_station_line

  ! padded array with a multiple of 512-bit chunks (16 4-byte integer size), +2 to add a bit-1 and string length
  ! for the single_station_line length above of 76, the padded integer array should have size 32
  integer, parameter :: size_padded_int = 32
  integer(kind=4), dimension(size_padded_int) :: padded_int
  ! or dynamically for variable lengths of single_station_line - just in case
  !integer(kind=4), dimension(:), allocatable :: padded_int

  ! the original MD5 algorithm seems to use big endian initializations
  !   https://datatracker.ietf.org/doc/html/rfc1321
  ! thus, to match the exact MD5, one would need to be careful here - todo if one wants to match it in future.
  ! for now, we're happy with something similar :)
  !
  ! for an explanation of the MD5 algorithm, see also:
  !   https://github.com/Zunawe/md5-c

  ! endianness
  ! little-endian writes the lowest order byte first, big-endian writes it last.
  ! for example, given the integer value of 1, big endian uses 00 00 00 01
  ! and little endian uses 01 00 00 00 as bit representation
  !
  ! using transfer(1,'a') takes the bit-representation of integer value 1 (multi-byte, usually 32-bit)
  ! and interprets it as a character type (of 1-byte length).
  ! thus, looking at the first byte, it is either 0 or 1, respectively.
  ! finally, ichar(c) takes a character and returns its position number.
  !
  ! generally, an intel machine would use little endian ordering
  !logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)

  ! string as integer
  !allocate(padded_int((len(single_station_line) + 64 + 8) / 4),stat=ier)

  ! padded array with a multiple of 512-bit chunks (16 4-byte integer size), +2 to add a bit-1 and string length
  !if (mod(int(len(single_station_line)/4) + 2, 16) == 0) then
  !  size_padded_int = int(len(single_station_line)/4) + 2
  !else
  !  size_padded_int = ((int(len(single_station_line)/4) + 2)/16 + 1) * 16
  !endif
  !allocate(padded_int(size_padded_int),stat=ier)
  !if (ier /= 0) stop 'Error allocating padded_int array'
  !padded_int(:) = 0

  ! initializes MD5 hash values
  ! binary integer part of the sines of integers (Radians) as constants
  do i = 1,64
    K(i) = floor(2**32 * dabs(dsin((i-1) + 1.d0)))
  enddo

  ! magic numbers: the original MD5 algorithm seems to use big endian initializations
  !if (is_big_endian) then
  !  a0 = int(z'67452301',kind=4)
  !  b0 = int(z'EFCDAB89',kind=4)
  !  c0 = int(z'98BADCFE',kind=4)
  !  d0 = int(z'10325476',kind=4)
  !else
  !  a0 = int(z'01234567')
  !  b0 = int(z'89ABCDEF')
  !  c0 = int(z'FEDCBA98')
  !  d0 = int(z'76543210')
  !endif
  ! note: using gfortran versions <= 9 will return a compilation error for these boz-initializations:
  !       "Error: Arithmetic overflow converting INTEGER(16) to INTEGER(4) .."
  !       thus, instead of
  !          b0 = int(z'89ABCDEF')  or  b0 = int(z'89ABCDEF',kind=4)
  !       one could split it and use
  !          b0 = ior(ishft(int(z'89AB'), 16), int(z'CDEF'))
  !       or
  !          b0 = transfer(int(z'89ABCDEF',kind=8),b0)
  !       these hexadecimals all fit into a 32-bit representation, and therefore the transfer() function
  !       should return valid results.
  a0 = transfer(int(z'01234567',kind=8),a0)
  b0 = transfer(int(z'89ABCDEF',kind=8),b0)
  c0 = transfer(int(z'FEDCBA98',kind=8),c0)
  d0 = transfer(int(z'76543210',kind=8),d0)

  do irec = 1,nrec
    single_station_line(:) = ''

    ! creates string station info
    write(single_station_line,'(a,a,3F12.4)') station_name(irec),network_name(irec),stlat(irec),stlon(irec),stbur(irec)

    !to check:
    !single_station_line = 'The quick brown fox jumps over the lazy dog.'
    !single_station_line = 'Hello, World!'
    !debug
    !print *,'debug: station ',irec,len_trim(single_station_line),'***'//single_station_line//'***'

    ! Convert the string to padded integer array
    str_len = len_trim(single_station_line)
    bit_len = str_len * 8

    ! size as a multiple of 4 characters
    if (mod(str_len,4) == 0) then
      size_str4 = str_len / 4
    else
      size_str4 = (int(str_len / 4) + 1) * 4
    endif

    ! initializes padded integer array
    padded_int(:) = 0

    ! merges 4 characters into a hex and converts it to integer
    do j = 1,size_str4
      ! string index
      i = (j-1) * 4 + 1

      ! integer conversion
      ! converts 4 characters to a 4-byte integer
      i1 = 0; i2 = 0; i3 = 0; i4 = 0

      ! 4 char to int
      !if (is_big_endian) then
      !  ! the ordering here is big-endian
      !  if (i <= str_len) i4 = ichar(single_station_line(i:i))
      !  if (i+1 <= str_len) i3 = ichar(single_station_line(i+1:i+1))
      !  if (i+2 <= str_len) i2 = ichar(single_station_line(i+2:i+2))
      !  if (i+3 <= str_len) i1 = ichar(single_station_line(i+3:i+3))
      !else
        ! little endian
        if (i <= str_len) i1 = ichar(single_station_line(i:i))
        if (i+1 <= str_len) i2 = ichar(single_station_line(i+1:i+1))
        if (i+2 <= str_len) i3 = ichar(single_station_line(i+2:i+2))
        if (i+3 <= str_len) i4 = ichar(single_station_line(i+3:i+3))
      !endif

      ! concats 4 char values to a hex value
      write(wtmp,'(4(z2.2))') i1,i2,i3,i4

      ! reads hexadecimal value as integer
      read(wtmp,'(z8)') padded_int(j)

      !debug
      !print *,'debug: i ',i,size_str4,str_len,'wtmp ***'//wtmp//'*** line ***'//single_station_line(i:i+3)//'***', &
      !        'position',j,' int = ',padded_int(j)
    enddo

    ! Append the bit '1' to the message
    !if (is_big_endian) then
    !  ! 00000000 00000000 00000000 00000001
    !  padded_int(str_len / 4 + 1) = ibset(padded_int(str_len / 4 + 1), 8 * (4 - mod(str_len, 4)))
    !else
      ! 10000000 00000000 00000000 00000000
      padded_int(str_len / 4 + 1) = ibset(padded_int(str_len / 4 + 1), 8 * (mod(str_len, 4) + 1))
    !endif

    ! Append the original bit length to the message
    padded_int(size_padded_int) = bit_len

    !debug
    !print *,'debug: padded size ',size_padded_int,'int:',padded_int(:)

    ! initialize
    a = a0
    b = b0
    c = c0
    d = d0

    ! md5 loop
    do i = 1, size_padded_int,16
      ! takes 512-bit chunk, i.e., 16 4-byte integers
      do j = 1,16
        pos = (i-1) + j
        w(j) = padded_int(pos)
      enddo

      ! main loop
      do j = 1,64
        if (j >= 1 .and. j <= 16) then
          f = ior( iand(b,c) , iand(not(b), d)  )
          g = j
        else if (j >= 17 .and. j <= 32) then
          f = ior( iand(d , b) , iand(not(d) , c)   )
          g = mod(5*(j-1) + 1, 16) + 1
        else if (j >= 33 .and. j <= 48) then
          f = ieor(b, ieor(c, d))
          g = mod(3*(j-1) + 5, 16) + 1
        else if (j >= 49 .and. j <= 64) then
          f = ieor(c, ior(b , not(d)))
          g = mod(7*(j-1), 16) + 1
        endif

        ! Save current hash values
        !if (is_big_endian) then
        !  temp = a
        !  a = d
        !  d = c
        !  c = b
        !  b = b + leftrotate((a + f + K(j) + w(g)), S(j))
        !  a = temp
        !else
          temp = d
          d = c
          c = b
          b = b + leftrotate((a + f + K(j) + w(g)), S(j))
          a = temp
        !endif
      enddo
    enddo

    ! add this station to hash
    a0 = a0 + a
    b0 = b0 + b
    c0 = c0 + c
    d0 = d0 + d
  enddo

  ! Output hash
  hashsum = iachar2hex(a0, b0, c0, d0)

contains

  ! Bitwise left rotation function
  function leftrotate(x, n) result(rotated)
    integer(4), intent(in) :: x, n
    integer(4) :: rotated

    rotated = ieor(shiftl(x, n), shiftr(x, 32 - n))
  end function leftrotate

  ! Convert array of integers to hexadecimal string
  function iachar2hex(a0,b0,c0,d0) result(hex_string)
    integer(4), intent(in) :: a0,b0,c0,d0
    character(len=32) :: hex_string
    integer(4) :: a,b,c,d

    ! zero extend
    !a = zext(a0, 8)
    !b = zext(b0, 8)
    !c = zext(c0, 8)
    !d = zext(d0, 8)
    ! w/out extend
    a = a0
    b = b0
    c = c0
    d = d0

    ! digest is a 128-bit number written in little endian
    write(hex_string, '(4Z8.8)') a,b,c,d

  end function iachar2hex

  ! Zero extend function
  !function zext(x, width) result(zexted)
  !  integer(4), intent(in) :: x, width
  !  integer(4) :: zexted
  !  zexted = iand(x, ishft(1, width) - 1)
  !end function zext

  end subroutine create_stations_checksum
