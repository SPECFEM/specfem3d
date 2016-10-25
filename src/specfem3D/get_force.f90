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

  subroutine get_force(tshift_force,hdur,lat,long,depth,NSOURCES,min_tshift_force_original,factor_force_source, &
                      comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP,user_source_time_function)

  use constants, only: IIN,IN_DATA_FILES,MAX_STRING_LEN,TINYVAL,mygroup,CUSTOM_REAL
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,EXTERNAL_STF,NSTEP_STF,NSOURCES_STF

  implicit none

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES

  double precision, intent(out) :: min_tshift_force_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_force,hdur,lat,long,depth,factor_force_source
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_E
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_N
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_Z_UP
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always rigth :
  !! in case of EXTERNAL_STF, it's equal to NSTEP,NSOURCES
  !! when .not. EXTERNAL_STF it' equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(out) :: user_source_time_function

  ! local variables below
  integer :: isource,dummyval
  double precision :: t_shift(NSOURCES)
  double precision :: length
  character(len=7) :: dummy
  character(len=MAX_STRING_LEN) :: string
  character(len=MAX_STRING_LEN) :: FORCESOLUTION,path_to_add
  character(len=MAX_STRING_LEN) :: external_stf_filename
  integer :: ier

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_force(:) = 0.d0
  hdur(:) = 0.d0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0

!
!---- read info
!
  FORCESOLUTION = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'FORCESOLUTION'
! see if we are running several independent runs in parallel
! if so, add the right directory for that run
! (group numbers start at zero, but directory names start at run0001, thus we add one)
! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    FORCESOLUTION = path_to_add(1:len_trim(path_to_add))//FORCESOLUTION(1:len_trim(FORCESOLUTION))
  endif

  open(unit=IIN,file=trim(FORCESOLUTION),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(FORCESOLUTION)
    stop 'Error opening FORCESOLUTION file'
  endif

! read source number isource
  do isource = 1,NSOURCES

    read(IIN,"(a)") string
    ! skips empty lines
    do while (len_trim(string) == 0)
      read(IIN,"(a)") string
    enddo

    ! read header with event information
    read(string,"(a6,i4)") dummy,dummyval

    ! read time shift
    read(IIN,"(a)") string
    read(string(12:len_trim(string)),*) t_shift(isource)

    ! read f0 (stored in hdur() array for convenience, to use the same array as for CMTSOLUTION)
    ! Please be careful, if you meet an error in reading the file FORCESOLUTION,
    ! such as you still write "hdur:" instead of "f0:"
    ! Please change your file or do following change in the code, such as changing
    ! read(string(4:len_trim(string)),*) hdur(isource)
    ! to
    ! read(string(6:len_trim(string)),*) hdur(isource)
    read(IIN,"(a)") string
    read(string(4:len_trim(string)),*) hdur(isource)

    ! read latitude
    read(IIN,"(a)") string
    read(string(10:len_trim(string)),*) lat(isource)

    ! read longitude
    read(IIN,"(a)") string
    read(string(11:len_trim(string)),*) long(isource)

    ! read depth
    read(IIN,"(a)") string
    read(string(7:len_trim(string)),*) depth(isource)

    ! read magnitude
    read(IIN,"(a)") string
    read(string(21:len_trim(string)),*) factor_force_source(isource)

    ! read direction vector's East component
    read(IIN,"(a)") string
    read(string(29:len_trim(string)),*) comp_dir_vect_source_E(isource)

    ! read direction vector's North component
    read(IIN,"(a)") string
    read(string(29:len_trim(string)),*) comp_dir_vect_source_N(isource)

    ! read direction vector's vertical component
    read(IIN,"(a)") string
    read(string(32:len_trim(string)),*) comp_dir_vect_source_Z_UP(isource)

    ! reads USER EXTERNAL SOURCE if needed
    if (EXTERNAL_STF) then
      ! gets external STF file name
      read(IIN,"(a)") string
      external_stf_filename = trim(string)

      ! reads in stf values
      call read_external_stf(isource,user_source_time_function,external_stf_filename)
    endif

  enddo

  close(IIN)

  ! Sets tshift_force to zero to initiate the simulation!
  if (NSOURCES == 1) then
    tshift_force = 0.d0
    min_tshift_force_original = t_shift(1)
  else
    tshift_force(1:NSOURCES) = t_shift(1:NSOURCES)-minval(t_shift)
    min_tshift_force_original = minval(t_shift)
  endif

  do isource = 1,NSOURCES
    ! checks half-duration
    ! half-duration is the dominant frequency of the source
    ! point forces use a Ricker source time function
    ! null half-duration indicates a very low-frequency source
    ! (see constants.h: TINYVAL = 1.d-9 )
    if (hdur(isource) < TINYVAL) hdur(isource) = TINYVAL

    ! check (tilted) force source direction vector
    length = sqrt( comp_dir_vect_source_E(isource)**2 + &
                   comp_dir_vect_source_N(isource)**2 + &
                   comp_dir_vect_source_Z_UP(isource)**2)
    if (length < TINYVAL) then
      print *, 'normal length: ', length
      print *, 'isource: ',isource
      stop 'error set force point normal length, make sure all forces have a non null direction vector'
    endif
  enddo

  end subroutine get_force

