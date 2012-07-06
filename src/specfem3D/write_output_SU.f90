!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  ! headers
  integer,parameter :: nheader=240      ! 240 bytes
  integer(kind=2) :: i2head(nheader/2)  ! 2-byte-integer
  integer(kind=4) :: i4head(nheader/4)  ! 4-byte-integer
  real(kind=4)    :: r4head(nheader/4)  ! 4-byte-real
  equivalence (i2head,i4head,r4head)    ! share the same 240-byte-memory

  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision :: x_found_source,y_found_source,z_found_source

  real(kind=4),dimension(:),allocatable :: rtmpseis
  integer :: ier

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
      form='unformatted', access='direct', recl=240+4*(NSTEP), iostat=ier)
  if( ier /= 0 ) stop 'error opening ***_dx_SU file'

  do irec_local = 1,nrec_local
   irec = number_receiver_global(irec_local)
   i4head(1)  =irec
   i4head(11) =z_found(irec)
   i4head(13) =z_found_source
   i4head(19) =x_found_source !utm_x_source(1)
   i4head(20) =y_found_source !utm_y_source(1)
   i4head(21) =x_found(irec)  !stutm_x(irec)
   i4head(22) =y_found(irec)  !stutm_y(irec)
   i2head(58) =NSTEP
   i2head(59) =DT*1.0d6

   ! convert to real 4-byte
   rtmpseis(1:NSTEP) = seismograms_d(1,irec_local,1:NSTEP)

   ! write record
   write(IOUT_SU,rec=irec_local) r4head, rtmpseis
  enddo
  close(IOUT_SU)

  ! write seismograms (dy)
  open(unit=IOUT_SU, file=trim(adjustl(final_LOCAL_PATH))//trim(adjustl(procname))//'_dy_SU' ,&
      form='unformatted', access='direct', recl=240+4*(NSTEP), iostat=ier)
  if( ier /= 0 ) stop 'error opening ***_dy_SU file'

  do irec_local = 1,nrec_local
   irec = number_receiver_global(irec_local)
   i4head(1)  =irec
   i4head(11) =z_found(irec)
   i4head(13) =z_found_source
   i4head(19) =x_found_source !utm_x_source(1)
   i4head(20) =y_found_source !utm_y_source(1)
   i4head(21) =x_found(irec)  !stutm_x(irec)
   i4head(22) =y_found(irec)  !stutm_y(irec)
   i2head(58) =NSTEP
   i2head(59) =DT*1.0d6

   ! convert to real 4-byte
   rtmpseis(1:NSTEP) = seismograms_d(2,irec_local,1:NSTEP)

   ! write record
   write(IOUT_SU,rec=irec_local) r4head, rtmpseis
  enddo
  close(IOUT_SU)

  ! write seismograms (dz)
  open(unit=IOUT_SU, file=trim(adjustl(final_LOCAL_PATH))//trim(adjustl(procname))//'_dz_SU' ,&
      form='unformatted', access='direct', recl=240+4*(NSTEP), iostat=ier)
  if( ier /= 0 ) stop 'error opening ***_dz_SU file'

  do irec_local = 1,nrec_local
   irec = number_receiver_global(irec_local)
   i4head(1)  =irec
   i4head(11) =z_found(irec)
   i4head(13) =z_found_source
   i4head(19) =x_found_source !utm_x_source(1)
   i4head(20) =y_found_source !utm_y_source(1)
   i4head(21) =x_found(irec)  !stutm_x(irec)
   i4head(22) =y_found(irec)  !stutm_y(irec)
   i2head(58) =NSTEP
   i2head(59) =DT*1.0d6

   ! convert to real 4-byte
   rtmpseis(1:NSTEP) = seismograms_d(3,irec_local,1:NSTEP)

   ! write record
   write(IOUT_SU,rec=irec_local) r4head, rtmpseis
  enddo
  close(IOUT_SU)

  end subroutine write_output_SU

