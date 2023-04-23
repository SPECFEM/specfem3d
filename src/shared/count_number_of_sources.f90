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

  subroutine count_number_of_sources(NSOURCES,sources_filename)

! determines number of sources depending on number of lines in source file
! (only executed by main process)

  use constants, only: IIN,IIN_PAR,MAX_STRING_LEN,IN_DATA_FILES, &
    NLINES_PER_FORCESOLUTION_SOURCE,NLINES_PER_CMTSOLUTION_SOURCE, &
    mygroup

  use shared_parameters, only: USE_FORCE_POINT_SOURCE,USE_EXTERNAL_SOURCE_FILE, &
    HAS_FINITE_FAULT_SOURCE,NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer,intent(out) :: NSOURCES
  character(len=MAX_STRING_LEN),intent(in) :: sources_filename

  ! local variables
  integer :: icounter,ier
  integer :: nlines_per_source
  character(len=MAX_STRING_LEN) :: fault_filename
  character(len=MAX_STRING_LEN) :: path_to_add
  character(len=MAX_STRING_LEN) :: dummystring

  ! initializes
  NSOURCES = 0

  ! checks if finite fault source defined
  fault_filename = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file_faults'
  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    fault_filename = path_to_add(1:len_trim(path_to_add))//fault_filename(1:len_trim(fault_filename))
  endif

  open(unit=IIN_PAR,file=trim(fault_filename),status='old',iostat=ier)
  if (ier == 0) then
    HAS_FINITE_FAULT_SOURCE = .true.
    !write(IMAIN,*) 'provides finite faults'
    close(IIN_PAR)
  else
    HAS_FINITE_FAULT_SOURCE = .false.
  endif

  ! checks if anything to do, finite fault simulations ignore CMT and force sources
  if (HAS_FINITE_FAULT_SOURCE) return

  ! get the number of lines describing the sources
  if (USE_FORCE_POINT_SOURCE) then
    ! FORCESOLUTION
    nlines_per_source = NLINES_PER_FORCESOLUTION_SOURCE
  else
    ! CMTSOLUTION
    nlines_per_source = NLINES_PER_CMTSOLUTION_SOURCE
  endif

  ! number of lines for source description
  ! in case of USE_EXTERNAL_SOURCE_FILE we have to read one additional line per source (the name of external source file)
  if (USE_EXTERNAL_SOURCE_FILE) then
    nlines_per_source = nlines_per_source + 1
  endif

  ! gets number of point sources
  open(unit=IIN,file=trim(sources_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    if (HAS_FINITE_FAULT_SOURCE) then
      ! no need for source file
      return
    else
      print *,'Error opening source SOLUTION file: ',trim(sources_filename)
      stop 'Error opening source SOLUTION file'
    endif
  endif

  icounter = 0
  do while (ier == 0)
    read(IIN,"(a)",iostat=ier) dummystring
    if (ier == 0) icounter = icounter + 1
  enddo
  close(IIN)

  ! checks lines are a multiple
  if (mod(icounter,nlines_per_source) /= 0) then
    if (USE_FORCE_POINT_SOURCE) then
      print *,'Error: total number of lines in FORCESOLUTION file should be a multiple of ',nlines_per_source
      stop 'Error total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'
    else
      print *,'Error: total number of lines in CMTSOLUTION file should be a multiple of ',nlines_per_source
      stop 'Error total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
    endif
  endif

  ! number of sources in file
  NSOURCES = icounter / nlines_per_source

  ! checks if any
  if (NSOURCES < 1) then
    print *,'Error: ',trim(sources_filename),' has ',icounter,'lines, but need ',nlines_per_source, &
            'per source... ',NSOURCES
    stop 'Error need at least one source in CMTSOLUTION or FORCESOLUTION file'
  endif

  end subroutine count_number_of_sources

