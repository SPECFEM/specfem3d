!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine read_external_stf(isource,user_source_time_function,external_stf_filename)

! reads in an external source time function file

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IOSTF

  use shared_parameters, only: NSTEP,DT,NSTEP_STF,NSOURCES_STF

  implicit none

  integer,intent(in) :: isource
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always rigth :
  !! in case of EXTERNAL_STF, it's equal to NSTEP,NSOURCES
  !! when .not. EXTERNAL_STF it' equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(inout) :: user_source_time_function

  character(len=MAX_STRING_LEN),intent(in) :: external_stf_filename

  ! local variables below
  integer :: i,ier
  double precision :: time_source,stf_val
  double precision :: time_source_start,time_source_end,dt_source
  double precision, parameter :: dt_tol = 1e-6
  character(len=256) :: line
  logical :: read_next_line

  ! initializes
  user_source_time_function(:,isource) = 0._CUSTOM_REAL

  ! opens specified file
  open(IOSTF,file=trim(external_stf_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open external source file: ',trim(external_stf_filename)
    stop 'Error opening external source time function file'
  endif

  ! debug
  !print *,'reading stf file: ',trim(external_stf_filename)

  ! gets number of file entries
  i = 0
  do while (ier == 0)
    read(IOSTF,"(a256)",iostat=ier) line
    !print *,i,'stf string: ',line
    if (ier == 0) then
      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! increases counter
      i = i + 1
    endif
  enddo
  rewind(IOSTF)

  ! checks number of entries
  if (i /= NSTEP) then
    print *,'Problem when reading external source time file: ', trim(external_stf_filename)
    print *,'  simulation number of time steps = ',NSTEP
    print *,'  source time function read number of entries = ',i
    print *,'Please make sure that the number of file entries correspond to the number of time steps in the simulation'
    stop 'Error invalid number of entries in external source time file'
  endif

  ! reads in external values
  do i = 1, NSTEP
    ! reads next valid line
    read_next_line = .true.
    do while (read_next_line)
      read(IOSTF,"(a256)",iostat=ier) line

      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! checks if done
      if (ier /= 0) then
        stop 'Error reading external source time file'
      else
        read_next_line = .false.
      endif
    enddo

    ! reads in values from line
    ! format: #time #stf-value
    read(line,*,iostat=ier) time_source, stf_val

    ! debug
    !print *,'stf: ',i,time_source,stf_val

    if (ier /= 0) then
      print *,'Problem when reading external source time file: ', trim(external_stf_filename)
      print *,'Please check, file format should be: #time #stf-value'
      stop 'Error reading external source time file with invalid format'
    endif

    ! sets STF
    user_source_time_function(i,isource) = stf_val

    ! to check time step size
    if (i == 1) then
      time_source_start = time_source
    else if (i == NSTEP) then
      time_source_end = time_source
    endif
  enddo

  ! checks if the time steps corresponds to the simulation time step
  dt_source =  (time_source_end - time_source_start) / dble(NSTEP-1)
  if (abs(dt_source - DT) > dt_tol ) then
    print *,'Error in time step in external source file ', trim(external_stf_filename)
    print *, ' simulation time step ', DT
    print *, ' source time function read time step ', dt_source
    stop 'Error invalid time step size in external STF file'
  endif

  ! closes external STF file
  close(IOSTF)

  end subroutine read_external_stf
