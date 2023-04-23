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


  subroutine get_force(FORCESOLUTION,tshift_src,hdur,lat,long,depth,NSOURCES, &
                       min_tshift_src_original,force_stf,factor_force_source, &
                       comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP, &
                       user_source_time_function)

  use constants, only: IIN,MAX_STRING_LEN,TINYVAL,CUSTOM_REAL
  use shared_parameters, only: USE_EXTERNAL_SOURCE_FILE,NSTEP_STF,NSOURCES_STF,DT,USE_RICKER_TIME_FUNCTION

  implicit none

!--- input or output arguments of the subroutine below

  character(len=MAX_STRING_LEN), intent(in) :: FORCESOLUTION
  integer, intent(in) :: NSOURCES
  integer, dimension(NSOURCES), intent(out) :: force_stf
  double precision, intent(out) :: min_tshift_src_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_src,hdur,lat,long,depth,factor_force_source
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_E
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_N
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_Z_UP
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always rigth :
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(out) :: user_source_time_function

  ! local variables below
  integer :: isource,ier,ipos,dummyval
  double precision :: t_shift(NSOURCES)
  double precision :: length
  character(len=7) :: dummy
  character(len=MAX_STRING_LEN) :: string
  character(len=MAX_STRING_LEN) :: external_source_time_function_filename

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0

  t_shift(:) = 0.d0
  tshift_src(:) = 0.d0
  hdur(:) = 0.d0

  force_stf(:) = 0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0

!
!---- read info
!
  open(unit=IIN,file=trim(FORCESOLUTION),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(FORCESOLUTION)
    stop 'Error opening FORCESOLUTION file'
  endif

! read source number isource
  do isource = 1,NSOURCES

    ! header
    read(IIN,"(a)") string
    ! skips empty lines
    do while (len_trim(string) == 0)
      read(IIN,"(a)") string
    enddo

    ! read header with event information
    ! format: FORCE  id
    ! as example: FORCE 001
    read(string,"(a6,i4)") dummy,dummyval  ! not used any further

    ! read time shift
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) t_shift(isource)
    else
      read(string(12:len_trim(string)),*) t_shift(isource)
    endif

    ! read hdur or f0 (stored in hdur() array for convenience, to use the same array as for CMTSOLUTION)
    ! Please be careful, if you meet an error in reading the file FORCESOLUTION,
    ! such as you still write "f0:" or "hdur:"
    ! Please change your file or do following change in the code, such as changing to
    ! read(string(4:len_trim(string)),*) hdur(isource)
    ! or
    ! read(string(6:len_trim(string)),*) hdur(isource)
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) hdur(isource)
    else
      read(string(10:len_trim(string)),*) hdur(isource)
    endif

    ! read latitude
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) lat(isource)
    else
      read(string(10:len_trim(string)),*) lat(isource)
    endif

    ! read longitude
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) long(isource)
    else
      read(string(11:len_trim(string)),*) long(isource)
    endif

    ! read depth
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) depth(isource)
    else
      read(string(7:len_trim(string)),*) depth(isource)
    endif

    ! source time function type
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) force_stf(isource)
    else
      read(string(22:len_trim(string)),*) force_stf(isource)
    endif

    ! read magnitude
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) factor_force_source(isource)
    else
      read(string(21:len_trim(string)),*) factor_force_source(isource)
    endif

    ! read direction vector's East component
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) comp_dir_vect_source_E(isource)
    else
      read(string(29:len_trim(string)),*) comp_dir_vect_source_E(isource)
    endif

    ! read direction vector's North component
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) comp_dir_vect_source_N(isource)
    else
      read(string(29:len_trim(string)),*) comp_dir_vect_source_N(isource)
    endif

    ! read direction vector's vertical component
    read(IIN,"(a)") string
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) comp_dir_vect_source_Z_UP(isource)
    else
      read(string(32:len_trim(string)),*) comp_dir_vect_source_Z_UP(isource)
    endif

    ! reads USER EXTERNAL SOURCE if needed
    if (USE_EXTERNAL_SOURCE_FILE) then
      ! gets external STF file name
      read(IIN,"(a)") string
      external_source_time_function_filename = trim(string)

      ! reads in stf values
      call read_external_source_time_function(isource,user_source_time_function,external_source_time_function_filename)
    endif

    ! checks Par_file flag to override type setting
    if (USE_RICKER_TIME_FUNCTION) force_stf(isource) = 1

    ! checks half-duration
    select  case(force_stf(isource))
    case (0)
      ! Gaussian
      ! null half-duration indicates a Dirac
      ! replace with very short Gaussian function
      if (hdur(isource) < 5. * DT ) hdur(isource) = 5. * DT
    case (1)
      ! Ricker source time function
      ! half-duration is the dominant frequency for the
      ! null half-duration indicates a very low-frequency source
      ! (see constants.h: TINYVAL = 1.d-9 )
      if (hdur(isource) < TINYVAL ) hdur(isource) = TINYVAL
    case (2)
      ! Step (Heaviside) source time function
      ! null half-duration indicates a Heaviside
      ! replace with very short error function
      if (hdur(isource) < 5. * DT ) hdur(isource) = 5. * DT
    case (3)
      ! Monochromatic source time function
      ! half-duration is the period
      ! (see constants.h: TINYVAL = 1.d-9 )
      if (hdur(isource) < TINYVAL ) then
        stop 'Error set force period, make sure all forces have a non-zero period'
      endif
    case (4)
      ! Gaussian by Meschede et al. (2011)
      ! null half-duration indicates a Dirac
      ! replace with very short Gaussian function
      if (hdur(isource) < 5. * DT ) hdur(isource) = 5. * DT
    case default
      stop 'unsupported source time function type (force_stf) value!'
    end select

  enddo

  close(IIN)

  ! Sets tshift_force to zero to initiate the simulation!
  if (NSOURCES == 1) then
    min_tshift_src_original = t_shift(1)
    tshift_src(1) = 0.d0
  else
    min_tshift_src_original = minval(t_shift)
    tshift_src(1:NSOURCES) = t_shift(1:NSOURCES) - min_tshift_src_original
  endif

  do isource = 1,NSOURCES
    ! checks half-duration
    ! half-duration is the dominant frequency of the source
    ! point forces use a Ricker source time function
    ! null half-duration indicates a very low-frequency source
    ! (see constants.h: TINYVAL = 1.d-9 )
    if (hdur(isource) < TINYVAL) hdur(isource) = TINYVAL

    ! check (tilted) force source direction vector
    length = sqrt( comp_dir_vect_source_E(isource)**2 &
                 + comp_dir_vect_source_N(isource)**2 &
                 + comp_dir_vect_source_Z_UP(isource)**2 )
    if (length < TINYVAL) then
      print *, 'normal length: ', length
      print *, 'isource: ',isource
      stop 'Error set force point normal length, make sure all forces have a non-zero direction vector'
    endif
  enddo

  end subroutine get_force

