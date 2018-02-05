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


!-------------------------------------------------------------------------------------------------
! Remove stations located outside of the mesh
!-------------------------------------------------------------------------------------------------

  subroutine station_filter(filename,filtered_filename,nfilter)

  use constants
  use specfem_par, only: SUPPRESS_UTM_PROJECTION,myrank,xstore,ystore

  implicit none

  ! input
  character(len=*) :: filename,filtered_filename

  ! output
  integer :: nfilter

  ! local
  integer,dimension(1) :: nrec, nrec_filtered
  integer :: ier
  double precision :: stlat,stlon,stele,stbur,stutm_x,stutm_y
  double precision :: minlat,minlon,maxlat,maxlon
  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name
  character(len=MAX_STRING_LEN) :: dummystring
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
  nrec = 0
  nrec_filtered = 0

  if (myrank == 0) then

    ! counts number of stations in stations file, filter them and output the list of active stations in STATIONS_FILTERED file
    open(unit=IIN, file=trim(filename), status = 'old', iostat = ier)
    if (ier /= 0) call exit_mpi(myrank, 'No file '//trim(filename)//', exit')
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    do while (ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier /= 0) exit

      if (len_trim(dummystring) > 0) then
        nrec(1) = nrec(1) + 1
        dummystring = trim(dummystring)
        read(dummystring, *) station_name, network_name, stlat, stlon, stele, stbur

        ! convert station location to UTM
        call utm_geo(stlon,stlat,stutm_x,stutm_y,ILONGLAT2UTM)

        ! counts stations within lon/lat region
        if (stutm_y >= LATITUDE_MIN .and. stutm_y <= LATITUDE_MAX .and. &
            stutm_x >= LONGITUDE_MIN .and. stutm_x <= LONGITUDE_MAX) then
          nrec_filtered(1) = nrec_filtered(1) + 1
          ! with specific format
          write(IOUT,'(a10,1x,a10,4e18.6)') &
                            trim(station_name),trim(network_name), &
                            sngl(stlat),sngl(stlon),sngl(stele),sngl(stbur)
        endif
      endif
    enddo

    close(IIN)
    close(IOUT)

    write(IMAIN,*)
    write(IMAIN,*) 'there are ',nrec(1),' stations in file ', trim(filename)
    write(IMAIN,*) 'saving ',nrec_filtered(1),' stations inside the model in file ', trim(filtered_filename)
    write(IMAIN,*) 'excluding ',nrec(1) - nrec_filtered(1),' stations located outside the model'
    write(IMAIN,*)

    if (nrec_filtered(1) < 1) then
      write(IMAIN,*) 'error filtered stations:'
      write(IMAIN,*) '  simulation needs at least 1 station but got ',nrec_filtered(1)
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
        write(IMAIN,*) '    UTM x min/max: ',LONGITUDE_MIN,LONGITUDE_MAX
        write(IMAIN,*) '    UTM y min/max : ',LATITUDE_MIN,LATITUDE_MAX
      endif

      write(IMAIN,*)
    endif

  endif ! myrank == 0

  call bcast_all_i(nrec,1)
  call bcast_all_i(nrec_filtered,1)

  nfilter = nrec_filtered(1)

  end subroutine station_filter
