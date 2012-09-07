!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine write_output_SU()

! writes out seismograms in SU (Seismic Unix) format

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  character(len=256) procname,final_LOCAL_PATH
  integer :: irec_local,irec

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source

  real(kind=4),dimension(:),allocatable :: rtmpseis
  real :: dx
  integer :: i,ier

  allocate(x_found(nrec),y_found(nrec),z_found(nrec),stat=ier)
  if( ier /= 0 ) stop 'error allocating array x_found y_found z_found'

  ! reads in station locations from output_list file
  open(unit=IIN_SU1,file=trim(OUTPUT_FILES)//'/output_list_stations.txt',status='unknown',iostat=ier)
  if( ier /= 0 ) stop 'error opening output_list_stations.txt file'

  do irec=1,nrec
   read(IIN_SU1,*) x_found(irec),y_found(irec),z_found(irec)
  enddo
  close(IIN_SU1)

  ! reads in source locations from output_list file
  open(unit=IIN_SU1,file=trim(OUTPUT_FILES)//'/output_list_sources.txt',status='unknown',iostat=ier)
  if( ier /= 0 ) stop 'error opening output_list_sources.txt file'

  read(IIN_SU1,*) x_found_source,y_found_source,z_found_source
  close(IIN_SU1)

  ! directory to store seismograms
  if( USE_OUTPUT_FILES_PATH ) then
   final_LOCAL_PATH = OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)) // '/'
  else
   ! create full final local path
   final_LOCAL_PATH = trim(adjustl(LOCAL_PATH)) // '/'
  endif
  write(procname,"(i4)") myrank

  allocate(rtmpseis(NSTEP),stat=ier)
  if( ier /= 0 ) stop 'error allocating rtmpseis array'

  ! write seismograms (dx)
  open(unit=IOUT_SU, file=trim(adjustl(final_LOCAL_PATH))//trim(adjustl(procname))//'_dx_SU' ,&
      status='unknown', access='direct', recl=4, iostat=ier)

  if( ier /= 0 ) stop 'error opening ***_dx_SU file'

  ! receiver interval
  if( nrec > 1) then
    dx = SNGL(x_found(2)-x_found(1))
  else
    dx = 0.0
  endif


  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)

    ! write section header
    call write_SU_header(irec_local,irec,nrec,NSTEP,DT,dx, &
                          x_found(irec),y_found(irec),z_found(irec),x_found_source,y_found_source,z_found_source)

    ! convert trace to real 4-byte
    rtmpseis(1:NSTEP) = seismograms_d(1,irec_local,1:NSTEP)
    do i=1,NSTEP
      write(IOUT_SU,rec=irec_local*60+(irec_local-1)*NSTEP+i) rtmpseis(i)
    enddo

  enddo
  close(IOUT_SU)

  ! write seismograms (dy)
  open(unit=IOUT_SU, file=trim(adjustl(final_LOCAL_PATH))//trim(adjustl(procname))//'_dy_SU' ,&
      status='unknown', access='direct', recl=4, iostat=ier)

  if( ier /= 0 ) stop 'error opening ***_dy_SU file'

  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)

    ! write section header
    call write_SU_header(irec_local,irec,nrec,NSTEP,DT,dx, &
                         x_found(irec),y_found(irec),z_found(irec),x_found_source,y_found_source,z_found_source)

    ! convert trace to real 4-byte
    rtmpseis(1:NSTEP) = seismograms_d(2,irec_local,1:NSTEP)
    do i=1,NSTEP
      write(IOUT_SU,rec=irec_local*60+(irec_local-1)*NSTEP+i) rtmpseis(i)
    enddo
  enddo
  close(IOUT_SU)

  ! write seismograms (dz)
  open(unit=IOUT_SU, file=trim(adjustl(final_LOCAL_PATH))//trim(adjustl(procname))//'_dz_SU' ,&
      status='unknown', access='direct', recl=4, iostat=ier)

  if( ier /= 0 ) stop 'error opening ***_dz_SU file'

  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)

    ! write section header
    call write_SU_header(irec_local,irec,nrec,NSTEP,DT,dx, &
                         x_found(irec),y_found(irec),z_found(irec),x_found_source,y_found_source,z_found_source)

    ! convert trace to real 4-byte
    rtmpseis(1:NSTEP) = seismograms_d(3,irec_local,1:NSTEP)
    do i=1,NSTEP
      write(IOUT_SU,rec=irec_local*60+(irec_local-1)*NSTEP+i) rtmpseis(i)
    enddo

  enddo
  close(IOUT_SU)

  end subroutine write_output_SU

!
!------------------------------------------------------------------------------------------------------
!

  subroutine write_SU_header(irec_local,irec,nrec,NSTEP,DT,dx, &
                             x_found,y_found,z_found,x_found_source,y_found_source,z_found_source)

  implicit none

  include "constants.h"

  integer :: irec_local,irec,nrec
  integer :: NSTEP
  double precision :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source
  double precision :: DT
  real :: dx

  ! local parameters
  integer(kind=2) :: header2(2)

  ! write SU headers (refer to Seismic Unix for details)
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+1)  irec                          ! receiver ID
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+10) NINT(x_found-x_found_source)  ! offset
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+11) NINT(z_found)

  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+13) NINT(z_found_source)                ! source location
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+19) NINT(x_found_source)                ! source location
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+20) NINT(y_found_source)                ! source location

  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+21) NINT(x_found)           ! receiver location xr
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+22) NINT(y_found)           ! receiver location zr

  if (nrec>1) write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+48) dx ! receiver interval

  ! time steps
  header2(1)=0  ! dummy
  header2(2)=NSTEP
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+29) header2

  ! time increment
  if( NINT(DT*1.0d6) < 65536 ) then
    header2(1)=NINT(DT*1.0d6)  ! deltat (unit: 10^{-6} second)
  else if( NINT(DT*1.0d3) < 65536 ) then
    header2(1)=NINT(DT*1.0d3)  ! deltat (unit: 10^{-3} second)
  else
    header2(1)=NINT(DT)  ! deltat (unit: 10^{0} second)
  endif
  header2(2)=0  ! dummy
  write(IOUT_SU,rec=(irec_local-1)*60+(irec_local-1)*NSTEP+30) header2

  end subroutine write_SU_header

