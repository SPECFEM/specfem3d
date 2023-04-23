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

  subroutine get_cmt(CMTSOLUTION,yr,jda,mo,da,ho,mi,sec, &
                     tshift_src,hdur,lat,long,depth,moment_tensor, &
                     DT,NSOURCES,min_tshift_src_original,user_source_time_function)

  use constants, only: IIN,MAX_STRING_LEN,CUSTOM_REAL
  use shared_parameters, only: USE_EXTERNAL_SOURCE_FILE,NSTEP_STF,NSOURCES_STF,NOISE_TOMOGRAPHY

  implicit none

!--- input or output arguments of the subroutine below

  character(len=MAX_STRING_LEN), intent(in) :: CMTSOLUTION
  integer, intent(in) :: NSOURCES
  double precision, intent(in) :: DT

  ! PDE time info (from first header line)
  integer, intent(out) :: yr,jda,ho,mi,mo,da
  double precision, intent(out) :: sec

  double precision, dimension(NSOURCES), intent(out) :: tshift_src,hdur,lat,long,depth
  double precision, dimension(6,NSOURCES), intent(out) :: moment_tensor
  double precision, intent(out) :: min_tshift_src_original

  !! VM VM use NSTEP_STF, NSOURCES_STF which are always right:
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(out) :: user_source_time_function

  ! local variables below
  integer :: julian_day,isource
  integer :: i,itype,istart,iend,ier,ipos
  double precision :: t_shift(NSOURCES)
  character(len=256) :: string
  character(len=MAX_STRING_LEN) :: external_source_time_function_filename

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_src(:) = 0.d0
  hdur(:) = 0.d0
  moment_tensor(:,:) = 0.d0

  ! origin time
  yr = 0; da = 0
  ho = -1; mi = -1; sec = -1.d0

!
!---- read hypocenter info
!
  open(unit = IIN,file=trim(CMTSOLUTION),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(CMTSOLUTION)
    stop 'Error opening CMTSOLUTION file'
  endif

  ! read source number isource
  do isource = 1,NSOURCES

    ! gets header line
    read(IIN,"(a256)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading header line in source ',isource
      stop 'Error reading header line in station in CMTSOLUTION file'
    endif

    ! skips empty lines
    do while (len_trim(string) == 0)
      read(IIN,"(a256)",iostat=ier) string
      if (ier /= 0) then
        print *, 'Error reading header line in source ',isource
        stop 'Error reading header line in station in CMTSOLUTION file'
      endif
    enddo

    ! debug
    !print *,'line ----',string,'----'

    ! reads header line with event information (assumes fixed format)
    ! old line: read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec

    ! reads header line with event information (free format)
    ! gets rid of the first datasource qualifyer string which can have variable length, like:
    ! "PDE 2014 9 3 .."
    ! " PDEQ2014 9 3 .."
    ! " MLI   1971   1   1 .."
    ! note: globalcmt.org solutions might have missing spaces after datasource qualifier
    !
    ! reads in year,month,day,hour,minutes,seconds
    istart = 1
    do itype = 1,6
      ! determines where first number starts
      do i = istart,len_trim(string)
        if (is_numeric(string(i:i))) then
          istart = i
          exit
        endif
      enddo
      if (istart >= len_trim(string)) stop 'Error determining datasource length in header line in CMTSOLUTION file'
      if (istart <= 1) stop 'Error determining datasource length in header line in CMTSOLUTION file'

      ! determines end and length of number
      iend = istart
      do i = istart,len_trim(string)
        if (itype /= 6) then
          ! integer values
          if (.not. is_numeric(string(i:i))) then
            iend = i
            exit
          endif
        else
          ! seconds will have a digit number
          ! digit numbers, e.g. 39.60, can contain '.'
          if (.not. is_digit(string(i:i))) then
            iend = i
            exit
          endif
        endif
      enddo
      iend = iend-1
      if (iend >= len_trim(string)) stop 'Error determining number length in header line in CMTSOLUTION file'
      if (iend < istart) stop 'Error determining number with negative length in header line in CMTSOLUTION file'

      ! debug
      !print *,itype,'line ----',string(istart:iend),'----'

      ! reads in event time information
      ! in case of multiple sources, time refers to the first entry only
      if (isource == 1) then
        select case (itype)
        case (1)
          ! year (as integer value)
          read(string(istart:iend),*) yr
        case (2)
          ! month (as integer value)
          read(string(istart:iend),*) mo
        case (3)
          ! day (as integer value)
          read(string(istart:iend),*) da
        case (4)
          ! hour (as integer value)
          read(string(istart:iend),*) ho
        case (5)
          ! minutes (as integer value)
          read(string(istart:iend),*) mi
        case (6)
          ! seconds (as float value)
          read(string(istart:iend),*) sec
        end select
      endif

      ! advances string
      istart = iend + 1
    enddo

    ! checks time information
    if (yr <= 0 .or. yr > 3000) then
      print *, 'Error reading year: ',yr,' in source ',isource,'is invalid'
      stop 'Error reading year out of header line in CMTSOLUTION file'
    endif
    if (mo < 1 .or. mo > 12) then
      print *, 'Error reading month: ',mo,' in source ',isource,'is invalid'
      stop 'Error reading month out of header line in CMTSOLUTION file'
    endif
    if (da < 1 .or. da > 31) then
      print *, 'Error reading day: ',da,' in source ',isource,'is invalid'
      stop 'Error reading day of header line in CMTSOLUTION file'
    endif
    if (ho < 0 .or. ho > 24) then
      print *, 'Error reading hour: ',ho,' in source ',isource,'is invalid'
      stop 'Error reading hour of header line in CMTSOLUTION file'
    endif
    if (mi < 0 .or. mi > 59) then
      print *, 'Error reading minute: ',mi,' in source ',isource,'is invalid'
      stop 'Error reading minute of header line in CMTSOLUTION file'
    endif
    if (sec < 0.0 .or. sec >= 60.0) then
      print *, 'Error reading second: ',sec,' in source ',isource,'is invalid'
      stop 'Error reading second of header line in CMTSOLUTION file'
    endif

    ! gets julian day number
    jda = julian_day(yr,mo,da)

    ! ignore line with event name
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading event name in source ',isource
      stop 'Error reading event name in station in CMTSOLUTION file'
    endif

    ! read time shift
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading time shift in source ',isource
      stop 'Error reading time shift in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) t_shift(isource)
    else
      read(string(12:len_trim(string)),*) t_shift(isource)
    endif

    ! read half duration
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading half duration in source ',isource
      stop 'Error reading half duration in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) hdur(isource)
    else
      read(string(15:len_trim(string)),*) hdur(isource)
    endif

    ! read latitude
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading latitude in source ',isource
      stop 'Error reading latitude in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) lat(isource)
    else
      read(string(10:len_trim(string)),*) lat(isource)
    endif

    ! read longitude
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading longitude in source ',isource
      stop 'Error reading longitude in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) long(isource)
    else
      read(string(11:len_trim(string)),*) long(isource)
    endif

    ! read depth
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading depth in source ',isource
      stop 'Error reading depth in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) depth(isource)
    else
      read(string(7:len_trim(string)),*) depth(isource)
    endif

    ! seismic moment tensor
    ! CMTSOLUTION: components given in dyne-cm
    ! read Mrr
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mrr in source ',isource
      stop 'Error reading Mrr in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(1,isource)   ! Mrr
    else
      read(string(5:len_trim(string)),*) moment_tensor(1,isource)   ! Mrr
    endif

    ! read Mtt
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mtt in source ',isource
      stop 'Error reading Mtt in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(2,isource)   ! Mtt
    else
      read(string(5:len_trim(string)),*) moment_tensor(2,isource)   ! Mtt
    endif

    ! read Mpp
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mpp in source ',isource
      stop 'Error reading Mpp in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(3,isource)   ! Mpp
    else
      read(string(5:len_trim(string)),*) moment_tensor(3,isource)   ! Mpp
    endif

    ! read Mrt
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mrt in source ',isource
      stop 'Error reading Mrt in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(4,isource)   ! Mrt
    else
      read(string(5:len_trim(string)),*) moment_tensor(4,isource)   ! Mrt
    endif

    ! read Mrp
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mrp in source ',isource
      stop 'Error reading Mrp in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(5,isource)   ! Mrp
    else
      read(string(5:len_trim(string)),*) moment_tensor(5,isource)   ! Mrp
    endif

    ! read Mtp
    read(IIN,"(a)",iostat=ier) string
    if (ier /= 0) then
      print *, 'Error reading Mtp in source ',isource
      stop 'Error reading Mtp in station in CMTSOLUTION file'
    endif
    ipos = index(string,':')
    if (ipos > 1 .and. ipos < len_trim(string)) then
      read(string(ipos+1:len_trim(string)),*) moment_tensor(6,isource)   ! Mtp
    else
      read(string(5:len_trim(string)),*) moment_tensor(6,isource)   ! Mtp
    endif

    ! checks half-duration
    ! null half-duration indicates a Heaviside
    ! replace with very short error function
    if (hdur(isource) < 5.d0 * DT) hdur(isource) = 5.d0 * DT

    ! reads USER EXTERNAL SOURCE if needed
    if (USE_EXTERNAL_SOURCE_FILE) then
      ! gets external STF file name
      read(IIN,"(a)") string
      external_source_time_function_filename = trim(string)

      ! reads in stf values
      call read_external_source_time_function(isource,user_source_time_function,external_source_time_function_filename)
    endif

  enddo

  close(IIN)

  ! noise simulations don't use the CMTSOLUTION source but a noise-spectrum source defined in S_squared
  if (NOISE_TOMOGRAPHY /= 0) hdur(:) = 0.d0

  ! If we're using external stf, don't worry about hdur.
  if (USE_EXTERNAL_SOURCE_FILE) hdur(:) = 0.d0

  ! Sets tshift_src to zero to initiate the simulation!
  if (NSOURCES == 1) then
    min_tshift_src_original = t_shift(1)
    tshift_src(1) = 0.d0
  else
    min_tshift_src_original = minval(t_shift)
    tshift_src(1:NSOURCES) = t_shift(1:NSOURCES) - min_tshift_src_original
  endif

  ! scales the moment tensor to Newton.m
  !
  ! CMTSOLUTION file values are in dyne.cm
  ! (from Global CMT project, Dziewonski 1981, Ekstrom et al. 2012, moment-tensor elements are given in dyne-cm)
  !
  ! 1 dyne is 1 gram * 1 cm / (1 second)^2
  ! 1 Newton is 1 kg * 1 m / (1 second)^2
  ! thus 1 Newton = 100,000 dynes
  ! therefore 1 dyne.cm = 1e-7 Newton.m
  !
  moment_tensor(:,:) = moment_tensor(:,:) * 1.d-7

  contains

  !--------------------------------------------------------------

  logical function is_numeric(char)

  ! returns .true. if input character is a number

  implicit none
  character(len=1), intent(in) :: char

  is_numeric = .false.

  if (index('0123456789', char) /= 0) then
    is_numeric = .true.
  endif

  end function

  !--------------------------------------------------------------

  logical function is_digit(char)

  ! returns .true. if input character is a number or a '.'

  implicit none
  character(len=1), intent(in) :: char

  is_digit = .false.

  if (index('0123456789.', char) /= 0) then
    is_digit = .true.
  endif

  end function

  end subroutine get_cmt

!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! calculates scalar moment (M0) in dyne-cm

  implicit none

  double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  ! local parameters
  double precision :: scalar_moment,scaleM

  ! scalar moment:
  ! see equation (1.4) in P.G. Silver and T.H. Jordan, 1982,
  ! "Optimal estimation of scalar seismic moment",
  ! Geophys. J.R. astr. Soc., 70, 755 - 787
  !
  ! or see equation (5.91) in Dahlen & Tromp (1998)
  !
  ! moment tensor M is a symmetric 3x3 tensor, and has six independent components
  !
  ! the euclidean matrix norm is invariant under rotation.
  ! thus, input can be:
  !   Mxx,Myy,Mzz,Mxy,Mxz,Myz
  ! or
  !   Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  !
  ! euclidean (or Frobenius) norm of a matrix: M0**2 = sum( Mij**2)
  scalar_moment = Mxx**2 + Myy**2 + Mzz**2 + 2.d0 * ( Mxy**2 + Mxz**2 + Myz**2 )

  ! adds 1/2 to be coherent with double couple or point sources
  scalar_moment = dsqrt(0.5d0*scalar_moment)

  ! scale factor for the scalar moment in dyne-cm
  !
  ! Mxx,.. are given in in Newton.m
  !
  ! 1 dyne is 1 gram * 1 cm / (1 second)^2
  ! 1 Newton is 1 kg * 1 m / (1 second)^2
  ! thus 1 Newton = 100,000 dynes
  ! therefore 1 dyne.cm = 1e-7 Newton.m
  !       and 1 Newton.m = 1e7 dyne.cm
  !
  scaleM = 1.d7

  ! return value (in dyne-cm)
  get_cmt_scalar_moment = scalar_moment * scaleM

  end function

!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_cmt_moment_magnitude(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! calculates scalar moment (M0)

  implicit none

  double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  ! local parameters
  double precision :: M0,Mw
  double precision,external :: get_cmt_scalar_moment

  ! scalar moment (in dyne-cm)
  M0 = get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)

  ! moment magnitude by Hanks & Kanamori, 1979
  ! Mw = 2/3 log( M0 ) - 10.7       (dyne-cm)
  !
  ! alternative forms:
  ! Mw = 2/3 ( log( M0 ) - 16.1 )   (N-m) "moment magnitude" by Hanks & Kanamori(1979) or "energy magnitude" by Kanamori (1977)
  !
  ! Aki & Richards ("Quantitative Seismology",2002):
  ! Mw = 2/3 ( log( M0 ) - 9.1 )    (N-m)
  !
  ! conversion: dyne-cm = 10**-7 N-m
  !
  ! we follow here the USGS magnitude policy:
  ! "Another source of confusion is the form of the formula for converting from scalar moment M0 to moment magnitude, M.
  !  The preferred practice is to use M = (log Mo)/1.5-10.7, where Mo is in dyne-cm (dyne-cm=10-7 N-m),
  !  the definition given by Hanks and Kanamori in 1979. An alternate form in Hanks and Kanamori's paper, M=(log M0-16.1)/1.5,
  !  is sometimes used, with resulting confusion. These formulae look as if they should yield the same result, but the latter
  !  is equivalent to M = (log Mo)/1.5-10.7333. The resulting round-off error occasionally leads to differences of 0.1
  !  in the estimates of moment magnitude released by different groups.
  !  All USGS statements of moment magnitude should use M = (log Mo)/1.5 - 10.7 = 2/3 (log M0) - 10.7 (Hanks & Kanamori, 1979)
  !  for converting from scalar moment M0 to moment magnitude. (..)"
  ! see: http://earthquake.usgs.gov/aboutus/docs/020204mag_policy.php
  !      https://web.archive.org/web/20160428095841/http://earthquake.usgs.gov:80/aboutus/docs/020204mag_policy.php

  if (M0 > 0.0d0) then
    ! this is to ensure M0>0.0 inorder to avoid arithmetic error.
    Mw = 2.d0/3.d0 * log10( max(M0,tiny(M0)) ) - 10.7
  else
    Mw = 0.0d0
  endif

  ! return value
  get_cmt_moment_magnitude = Mw

  end function
