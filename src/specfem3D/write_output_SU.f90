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
  integer :: istore
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms

  ! local parameters
  character(len=MAX_STRING_LEN) :: procname,final_LOCAL_PATH
  integer :: irec_local,irec

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source

  real(kind=4),dimension(:),allocatable :: rtmpseis
  real :: dx
  integer :: ier

  ! arrays for Seismic Unix header
  integer, dimension(28) :: header1
  real(kind=4), dimension(30) :: header4
  integer(kind=2) :: header2(2),header3(2)

  allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
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

  allocate(rtmpseis(NSTEP),stat=ier)
  if (ier /= 0) stop 'error allocating rtmpseis array'

  ! deletes old files
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_dx_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_dy_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_dz_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_vx_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_vy_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_vz_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_ax_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_ay_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_az_SU')
  close(IIN_SU1,status='delete')
  open(unit=IIN_SU1,file=trim(final_LOCAL_PATH)//trim(procname)//'_p_SU')
  close(IIN_SU1,status='delete')

  if (istore == 1) then
    ! open seismograms displacement
    open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_dx_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_dx_SU file'
    open(unit=IIN_SU2, file=trim(final_LOCAL_PATH)//trim(procname)//'_dy_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_dy_SU file'
    open(unit=IIN_SU3, file=trim(final_LOCAL_PATH)//trim(procname)//'_dz_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_dz_SU file'
  else if (istore == 2) then
    ! open seismograms velocity
    open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_vx_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_vx_SU file'
    open(unit=IIN_SU2, file=trim(final_LOCAL_PATH)//trim(procname)//'_vy_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_vy_SU file'
    open(unit=IIN_SU3, file=trim(final_LOCAL_PATH)//trim(procname)//'_vz_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_vz_SU file'
  else if (istore == 3) then
    ! open seismograms acceleration
    open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_ax_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_ax_SU file'
    open(unit=IIN_SU2, file=trim(final_LOCAL_PATH)//trim(procname)//'_ay_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_ay_SU file'
    open(unit=IIN_SU3, file=trim(final_LOCAL_PATH)//trim(procname)//'_az_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_az_SU file'
  else
     ! open seismogram pressure
    open(unit=IIN_SU1, file=trim(final_LOCAL_PATH)//trim(procname)//'_p_SU', &
         status='unknown', access='stream', iostat=ier)
    if (ier /= 0) stop 'error opening ***_p_SU file'
  endif


    ! receiver interval
    if (nrec > 1) then
      dx = SNGL(x_found(2)-x_found(1))
    else
      dx = 0.0
    endif

    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)

      ! write section header
      call write_SU_header(irec,dx,header1,header2,header3,header4, &
                           x_found(irec),y_found(irec),z_found(irec),x_found_source,y_found_source,z_found_source)

      write(IIN_SU1,pos=4*(irec_local-1)*(60+NSTEP) + 1 ) header1,header2,header3,header4,seismograms(1,irec_local,:)

      if (istore /= 4) then
        write(IIN_SU2,pos=4*(irec_local-1)*(60+NSTEP) + 1 ) header1,header2,header3,header4,seismograms(2,irec_local,:)
        write(IIN_SU3,pos=4*(irec_local-1)*(60+NSTEP) + 1 ) header1,header2,header3,header4,seismograms(3,irec_local,:)
      endif

    enddo

  close(IIN_SU1)

  if (istore /= 4) then
    close(IIN_SU2)
    close(IIN_SU3)
  endif

  end subroutine write_output_SU

!
!------------------------------------------------------------------------------------------------------
!

  subroutine write_SU_header(irec,dx,header1,header2,header3,header4, &
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

  header1(:)=0
  header2(:)=0
  header3(:)=0
  header4(:)=0.0

  ! write SU headers (refer to Seismic Unix for details)
  header1(1) = irec                          ! receiver ID
  header1(10)= NINT(x_found-x_found_source)  ! offset
  header1(11)= NINT(z_found)

  header1(13)= NINT(z_found_source)                ! source location
  header1(19)= NINT(x_found_source)                ! source location
  header1(20)= NINT(y_found_source)                ! source location

  header1(21)= NINT(x_found)           ! receiver location xr
  header1(22)= NINT(y_found)           ! receiver location zr

  if (nrec > 1) header4(18)= dx ! receiver interval

  ! time steps
  header2(1)=0  ! dummy
  header2(2)=int(NSTEP, kind=2)

  ! time increment
  if (NINT(DT*1.0d6) < 65536) then
    header3(1)=NINT(DT*1.0d6, kind=2)  ! deltat (unit: 10^{-6} second)
  else if (NINT(DT*1.0d3) < 65536) then
    header3(1)=NINT(DT*1.0d3, kind=2)  ! deltat (unit: 10^{-3} second)
  else
    header3(1)=NINT(DT, kind=2)  ! deltat (unit: 10^{0} second)
  endif
  header3(2)=0  ! dummy

  end subroutine write_SU_header

