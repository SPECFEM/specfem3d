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


!-------------------------------------------------------------------------------------------------
! Remove stations located outside of the mesh
!-------------------------------------------------------------------------------------------------

  subroutine station_filter(filename,filtered_filename,nrec_filtered)

  use constants, only: CUSTOM_REAL,MAX_LENGTH_NETWORK_NAME,MAX_LENGTH_STATION_NAME,MAX_STRING_LEN, &
    IIN,IOUT,IMAIN,ILONGLAT2UTM,IUTM2LONGLAT

  use specfem_par, only: SUPPRESS_UTM_PROJECTION,myrank,xstore,ystore

  implicit none

  ! input
  character(len=*),intent(in) :: filename,filtered_filename

  ! output
  integer,intent(out) :: nrec_filtered

  ! local
  integer :: nrec_all,ier
  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision :: minlat,minlon,maxlat,maxlon

  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name

  character(len=MAX_STRING_LEN) :: line
  character(len=32) :: fmt_str

  real(kind=CUSTOM_REAL):: minl,maxl,min_all,max_all
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  ! gets model dimensions
  minl = minval( xstore )
  maxl = maxval( xstore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LONGITUDE_MIN = min_all
  LONGITUDE_MAX = max_all

  minl = minval( ystore )
  maxl = maxval( ystore )
  call min_all_all_cr(minl,min_all)
  call max_all_all_cr(maxl,max_all)
  LATITUDE_MIN = min_all
  LATITUDE_MAX = max_all

  ! initialization
  nrec_all = 0
  nrec_filtered = 0

  if (myrank == 0) then

    ! counts number of stations in stations file, filter them and output the list of active stations in STATIONS_FILTERED file
    open(unit=IIN, file=trim(filename), status = 'old', iostat = ier)
    if (ier /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')

    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) line
      if (ier /= 0) exit

      if (len_trim(line) > 0 .and. line(1:1) /= '#') then
        nrec_all = nrec_all + 1

        line = trim(line)
        read(line, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y,ILONGLAT2UTM)

        ! counts stations within lon/lat region
        if (stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
            stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then
          nrec_filtered = nrec_filtered + 1

          ! check station name length
          if (len_trim(station_name) > 32) stop 'Station name length exceeds bounds of 32 characters'
          ! check network name length
          if (len_trim(network_name) > 10) stop 'Network name length exceeds bounds of 10 characters'

          ! with specific format
          ! fixing the format is needed for some compilers
          ! (say, cray Fortran would write same numbers as 2*1500.0 with write(..,*))
          ! however, we might loose location resolution if the range is not appropriate for the specifier (e,f,..).
          ! we thus try to estimate a good format string based on the maximum value of the numbers to output
          if (max(abs(stlat),abs(stlon),abs(stele),abs(stbur)) >= 1.d9) then
            ! uses exponential numbers format
            fmt_str = '(a32,1x,a10,4e24.10)'
          else if (max(abs(stlat),abs(stlon),abs(stele),abs(stbur)) >= 1.d-1) then
            ! uses float numbers format
            fmt_str = '(a32,1x,a10,4f24.12)'
          else
            ! uses exponential numbers format
            fmt_str = '(a32,1x,a10,4e24.10)'
          endif

          write(IOUT,fmt_str) trim(station_name),trim(network_name),stlat,stlon,stele,stbur

        endif
      endif
    enddo

    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec_all,' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec_all - nrec_filtered,' stations located outside the model'
    write(IMAIN,*)

    if (nrec_filtered < 1) then
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered
      write(IMAIN,*)
      write(IMAIN,*) '  check that stations in file '//trim(filename)//' are within'

      if (SUPPRESS_UTM_PROJECTION) then
        write(IMAIN,*) '    latitude min/max : ',LATITUDE_MIN,LATITUDE_MAX
        write(IMAIN,*) '    longitude min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
      else
        ! convert edge locations from UTM back to lat/lon
        call utm_geo(minlon,minlat,LONGITUDE_MIN,LATITUDE_MIN,IUTM2LONGLAT)
        call utm_geo(maxlon,maxlat,LONGITUDE_MAX,LATITUDE_MAX,IUTM2LONGLAT)
        write(IMAIN,*) '    longitude min/max: ',minlon,maxlon
        write(IMAIN,*) '    latitude min/max : ',minlat,maxlat
        write(IMAIN,*) '    UTM x min/max    : ',LONGITUDE_MIN,LONGITUDE_MAX
        write(IMAIN,*) '    UTM y min/max    : ',LATITUDE_MIN,LATITUDE_MAX
      endif

      write(IMAIN,*)
    endif

  endif ! myrank == 0

  ! broadcast to all processes
  call bcast_all_singlei(nrec_filtered)

  end subroutine station_filter
