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

  subroutine read_external_source_time_function(isource,user_source_time_function,external_source_time_function_filename)

! reads in an external source time function file

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IO_STF

  use shared_parameters, only: NSTEP,NSTEP_STF,NSOURCES_STF

  implicit none

  integer,intent(in) :: isource
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always OK:
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(inout) :: user_source_time_function

  character(len=MAX_STRING_LEN),intent(in) :: external_source_time_function_filename

  ! local variables below
  integer :: i,ier,str_len
  real(kind=CUSTOM_REAL) :: stf_val
  character(len=256) :: line

  ! threshold value to issue a warning
  real(kind=CUSTOM_REAL),parameter :: SMALL_STF_VAL = 0.001

  ! debug timing
  !double precision :: tstart,tCPU
  !double precision, external :: wtime

  ! debug timing
  !tstart = wtime()

  ! saftey check
  if (NSTEP /= NSTEP_STF) then
    print *,'Error: external source time function should have NSTEP_STF = ',NSTEP_STF,' equal to NSTEP = ',NSTEP
    stop 'Error invalid number of NSTEP_STF'
  endif

  ! for STF file names like "***.bin" - ending with ".bin", we will read it in as a binary file instead of as an ASCII file
  ! to enhance the speed for reading in many STF files in case.

  ! check file ending
  ! ".bin" for binary source time function files
  str_len = len_trim(external_source_time_function_filename)
  if (str_len > 4) then
    if (external_source_time_function_filename(str_len-3:str_len) == ".bin") then
      ! binary source time function file
      ! format: can only have numbers, no comment lines
      !         also, it must have the exact NSTEP length (and values) for the simulation as no further checks will be done
      !
      ! open binary file
      open(IO_STF,file=trim(external_source_time_function_filename),form='unformatted',status='old',action='read',iostat=ier)
      if (ier /= 0) then
        print *,'Error: could not open external binary source file: ',trim(external_source_time_function_filename)
        print *,'       for source index: ',isource
        stop 'Error opening external binary source time function file'
      endif
      read(IO_STF) user_source_time_function(:,isource)
      close(IO_STF)

      ! debug timing
      !tCPU = wtime() - tstart
      !print *,'debug: binary read source ',isource,' - elapsed time: ',sngl(tCPU),'s'

      ! all done, no further checking (for speed)
      return
    endif
  endif

  ! clear the array for that source
  user_source_time_function(:,isource) = 0._CUSTOM_REAL

  ! ASCII source time function file
  ! format: the ASCII source time function file can use a format like
  !         # user comment
  !         ! another user comment
  !         stf_value1
  !         stf_value2
  !         ..
  !         stf_valueNSTEP
  !
  ! opens specified file
  open(IO_STF,file=trim(external_source_time_function_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open external ASCII source file: ',trim(external_source_time_function_filename)
    print *,'       for source index: ',isource
    stop 'Error opening external ASCII source time function file'
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
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle
      ! stop 'error in format of external_source_time_function_filename, no comments are allowed in it'

      ! increases counter
      i = i + 1
    endif
  enddo
  rewind(IO_STF)

  ! checks
  if (i < 1) then
    print *,'Error: External ASCII source time function file ',trim(trim(external_source_time_function_filename)), &
            'has no valid data;'
    print *,'       the number of time steps is < 1. Please check the file...'
    stop 'Error: the number of time steps in external_source_time_function_filename is < 1'
  endif

  if (i > NSTEP_STF) then
    print *
    print *,'****************************************************************************************'
    print *,'Warning: ',trim(external_source_time_function_filename), &
            ' contains ',i,' time steps'
    print *,'         only the first NSTEP_STF = ',NSTEP_STF,' will be read, all the others will be ignored.'
    print *,'****************************************************************************************'
    print *
  endif

  ! checks number of time steps read
  if (i < NSTEP_STF) then
    print *,'Problem when reading external ASCII source time file: ', trim(external_source_time_function_filename)
    print *,'  number of time steps in the simulation                  = ',NSTEP_STF
    print *,'  number of time steps read from the source time function = ',i
    print *,'Please make sure that the number of time steps in the external ASCII source file read is greater or &
             &equal to the number of time steps in the simulation'
    stop 'Error invalid number of time steps in external ASCII source time file'
  endif

  ! read the time step used and check that it is the same as DT used for the code
  ier = 0
  i = 0
  do while (ier == 0)
    read(IO_STF,"(a256)",iostat=ier) line
    if (ier == 0) then
      ! suppress leading white spaces, if any
      line = adjustl(line)

      ! skip empty/comment lines
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! reads the STF values
      read(line,*,iostat=ier) stf_val
      if (ier /= 0) then
        print *,'Problem when reading external ASCII source time file: ', trim(external_source_time_function_filename)
        print *,'Please check, file format should be: '
        print *,'  # DT-time-step-size'
        print *,'  # stf-value'
        print *,'  # ..'
        stop 'Error reading external ASCII source time file with invalid format'
      endif

      ! increases counter
      i = i + 1

      ! checks first STF value
      if (i == 1 .and. stf_val > SMALL_STF_VAL) then
        print *
        print *,'****************************************************************************************'
        print *,'Warning: ',trim(external_source_time_function_filename), &
                ' starts with an STF value ',stf_val
        print *,'         Onset values should be close to zero to avoid numerical, high-frequency oscillations artifacts'
        print *,'****************************************************************************************'
        print *
      endif

      ! stores source time function
      user_source_time_function(i,isource) = stf_val

      ! checks if all steps read
      if (i == NSTEP_STF) exit
    endif
  enddo

  ! closes external STF file
  close(IO_STF)

  ! debug timing
  !tCPU = wtime() - tstart
  !print *,'debug: ascii read source ',isource,' - elapsed time: ',sngl(tCPU),'s'

  end subroutine read_external_source_time_function

