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

  subroutine write_output_SU(seismograms,istore)

! writes out seismograms in SU (Seismic Unix) format

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! arguments
  integer,intent(in) :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local*NB_RUNS_ACOUSTIC_GPU,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: seismograms

  ! local parameters
  character(len=MAX_STRING_LEN) :: procname,final_LOCAL_PATH
  integer :: irec_local,irec,ioffset,ier

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source
  real :: dx

  ! arrays for Seismic Unix header
  integer, dimension(28) :: header1
  real(kind=4), dimension(30) :: header4
  integer(kind=2) :: header2(2),header3(2)

  character(len=1),parameter :: comp(4) = (/ 'd', 'v', 'a', 'p' /)

  allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2189')
  if (ier /= 0) stop 'error allocating arrays x_found y_found z_found'

  ! reads in station locations from output_list file
  open(unit=IIN_SU1,file=trim(OUTPUT_FILES)//'/output_list_stations.txt',status='old',iostat=ier)
  if (ier /= 0) stop 'error opening output_list_stations.txt file'

  do irec=1,nrec
   read(IIN_SU1,*) station_name(irec),network_name(irec),x_found(irec),y_found(irec),z_found(irec)
  enddo
  close(IIN_SU1)

  ! reads in source locations from output_list file
  open(unit=IIN_SU1,file=trim(OUTPUT_FILES)//'/output_list_sources.txt',status='old',iostat=ier)
  if (ier /= 0) stop 'error opening output_list_sources.txt file'

  read(IIN_SU1,*) x_found_source,y_found_source,z_found_source
  close(IIN_SU1)

  ! directory to store seismograms
  final_LOCAL_PATH = OUTPUT_FILES(1:len_trim(OUTPUT_FILES)) // '/'
  write(procname,"(i4)") myrank
  procname = adjustl(procname)

  select case(istore)
  case (1,2,3)
    ! open seismograms displacement, velocity, acceleration
    if (seismo_offset == 0) then
      open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'x_SU', &
           status='replace', access='stream', form='unformatted', action='write', iostat=ier)
      if (ier /= 0) stop 'error opening ***x_SU file'
      open(unit=IIN_SU2, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'y_SU', &
           status='replace', access='stream', form='unformatted', action='write', iostat=ier)
      if (ier /= 0) stop 'error opening ***y_SU file'
      open(unit=IIN_SU3, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'z_SU', &
           status='replace', access='stream', form='unformatted', action='write', iostat=ier)
      if (ier /= 0) stop 'error opening ***z_SU file'
    else
      open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'x_SU', &
           status='old', access='stream', form='unformatted', action='readwrite', iostat=ier)
      if (ier /= 0) stop 'error opening ***x_SU file'
      open(unit=IIN_SU2, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'y_SU', &
           status='old', access='stream', form='unformatted', action='readwrite', iostat=ier)
      if (ier /= 0) stop 'error opening ***y_SU file'
      open(unit=IIN_SU3, file=trim(final_LOCAL_PATH)//trim(procname)//'_'//comp(istore)//'z_SU', &
           status='old', access='stream', form='unformatted', action='readwrite', iostat=ier)
      if (ier /= 0) stop 'error opening ***z_SU file'
    endif
  case (4)
    ! open seismogram pressure
    if (seismo_offset == 0) then
      open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_p_SU', &
           status='replace', access='stream', form='unformatted', action='write', iostat=ier)
      if (ier /= 0) stop 'error opening ***_p_SU file'
    else
      open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_p_SU', &
           status='old', access='stream', form='unformatted', action='readwrite', iostat=ier)
      if (ier /= 0) stop 'error opening ***_p_SU file'
    endif
  case default
    stop 'Invalid istore value in writing output SU'
  end select

  ! receiver interval
  if (nrec > 1) then
    dx = SNGL(x_found(2)-x_found(1))
  else
    dx = 0.0
  endif

  do irec_local = 1,nrec_local*NB_RUNS_ACOUSTIC_GPU
    irec = number_receiver_global(mod(irec_local,nrec_local))

    if (seismo_offset == 0) then
      ! determines header
      call determine_SU_header(irec,dx,header1,header2,header3,header4, &
                               x_found(irec),y_found(irec),z_found(irec),x_found_source,y_found_source,z_found_source)
      ! writes section header
      ! position in bytes
      ioffset = 4*(irec_local-1)*(60+NSTEP) + 1
      select case (istore)
      case (1,2,3)
        write(IIN_SU1,pos=ioffset) header1,header2,header3,header4
        write(IIN_SU2,pos=ioffset) header1,header2,header3,header4
        write(IIN_SU3,pos=ioffset) header1,header2,header3,header4
      case (4)
        write(IIN_SU1,pos=ioffset) header1,header2,header3,header4
      end select
    endif

    ! writes seismos
    ! position in bytes
    ioffset = 4*(irec_local-1)*(60+NSTEP) + 4 * 60 + 4 * seismo_offset + 1
    select case (istore)
    case (1,2,3)
      write(IIN_SU1,pos=ioffset) seismograms(1,irec_local,1:seismo_current)
      write(IIN_SU2,pos=ioffset) seismograms(2,irec_local,1:seismo_current)
      write(IIN_SU3,pos=ioffset) seismograms(3,irec_local,1:seismo_current)
    case (4)
      write(IIN_SU1,pos=ioffset) seismograms(1,irec_local,1:seismo_current)
    end select
  enddo

  select case (istore)
  case (1,2,3)
    close(IIN_SU1)
    close(IIN_SU2)
    close(IIN_SU3)
  case (4)
    close(IIN_SU1)
  end select

  end subroutine write_output_SU

!
!------------------------------------------------------------------------------------------------------
!

  subroutine determine_SU_header(irec,dx,header1,header2,header3,header4, &
                                 x_found,y_found,z_found,x_found_source,y_found_source,z_found_source)


  use constants

  use specfem_par, only: nrec,NSTEP,DT

  implicit none

  integer :: irec
  double precision :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source
  real :: dx

! Arrays for Seismic Unix header
  integer, dimension(28) :: header1
  real(kind=4), dimension(30) :: header4
  integer(kind=2) :: header2(2),header3(2)

! info format, see:
! http://sepwww.stanford.edu/oldsep/cliner/files/suhelp/suhelp.html
! https://github.com/JohnWStockwellJr/SeisUnix/wiki/Seismic-Unix-data-format
! http://www.geophysik.uni-jena.de/Lehre/Master+of+Science/Reflexionsseismik/SEG_Y+Trace+Header+Format.html

  header1(:) = 0
  header2(:) = 0
  header3(:) = 0
  header4(:) = 0.0

  ! write SU headers (refer to Seismic Unix for details)
  header1(1)  = irec                          ! receiver ID
  header1(10) = NINT(x_found-x_found_source)  ! offset
  header1(11) = NINT(z_found)                 ! elevation

  header1(13) = NINT(z_found_source)                ! source location (depth)
  header1(19) = NINT(x_found_source)                ! source location (X)
  header1(20) = NINT(y_found_source)                ! source location (Y)

  header1(21) = NINT(x_found)           ! receiver location xr
  header1(22) = NINT(y_found)           ! receiver location zr

  if (nrec > 1) header4(18) = dx ! receiver interval

  ! time steps
  header2(1) = 0  ! dummy
  header2(2) = int(NSTEP, kind=2)

  ! time increment
  if (NINT(DT*1.0d6) < 65536) then
    header3(1) = NINT(DT*1.0d6, kind=2)  ! deltat (unit: 10^{-6} second)
  else if (NINT(DT*1.0d3) < 65536) then
    header3(1) = NINT(DT*1.0d3, kind=2)  ! deltat (unit: 10^{-3} second)
  else
    header3(1) = NINT(DT, kind=2)  ! deltat (unit: 10^{0} second)
  endif
  header3(2) = 0  ! dummy

  end subroutine determine_SU_header

