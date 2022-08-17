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

  subroutine read_external_source_time_function(isource,user_source_time_function,external_source_time_function_filename)

! reads in an external source time function file

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IO_STF

  use shared_parameters, only: NSTEP,NSTEP_STF,NSOURCES_STF,DT

  implicit none

  integer,intent(in) :: isource
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always OK:
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(inout) :: user_source_time_function

  character(len=MAX_STRING_LEN),intent(in) :: external_source_time_function_filename

  ! local variables below
  integer :: i,ier
  double precision :: dt_source
  character(len=256) :: line

  ! clear the array for that source
  user_source_time_function(:,isource) = 0._CUSTOM_REAL

  ! opens specified file
  open(IO_STF,file=trim(external_source_time_function_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open external source file: ',trim(external_source_time_function_filename)
    stop 'Error opening external source time function file'
  endif

  ! gets number of file entries
  i = 0
  do while (ier == 0)
    read(IO_STF,"(a256)",iostat=ier) line
    if (ier == 0) then
      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') &
        stop 'error in format of external_source_time_function_filename, no comments are allowed in it'

      ! increases counter
      i = i + 1
    endif
  enddo
  rewind(IO_STF)

  ! subtract one because the first line of the file contains the time step used
  i = i - 1

  if (i < 1) stop 'error: the number of time steps in external_source_time_function_filename is < 1'

  if (i > NSTEP_STF) then
    print *
    print *,'****************************************************************************************'
    print *,'Warning: external_source_time_function_filename contains more than NSTEP_STF time steps,'
    print *,'         only the first NSTEP_STF will be read, all the others will be ignored.'
    print *,'****************************************************************************************'
    print *
  endif

  ! checks number of time steps read
  if (i < NSTEP) then
    print *,'Problem when reading external source time file: ', trim(external_source_time_function_filename)
    print *,'  number of time steps in the simulation = ',NSTEP
    print *,'  number of time steps read from the source time function = ',i
    print *,'Please make sure that the number of time steps in the external source file read is greater or &
             &equal to the number of time steps in the simulation'
    stop 'Error invalid number of time steps in external source time file'
  endif

  ! read the time step used and check that it is the same as DT used for the code
  read(IO_STF,*) dt_source
  if (abs((dt_source - DT) / DT) > 1.d-3) stop 'error: the external source time file does not use the same time step as DT'

  ! read the source values
  do i = 1, NSTEP_STF
    read(IO_STF,*,iostat=ier) user_source_time_function(i,isource)
    if (ier /= 0) then
      print *,'Problem when reading external source time file: ', trim(external_source_time_function_filename)
      print *,'Please check, file format should be: #time #stf-value'
      stop 'Error reading external source time file with invalid format'
    endif
  enddo

  ! closes external STF file
  close(IO_STF)

  end subroutine read_external_source_time_function

