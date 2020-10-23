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

module passive_imaging_format_mod

  use interpolation_mod, only: si, sp, di, dp, cp, hp, deg2rad, rad2deg
  use rotations_mod, only: define_mesh_rotation_matrix, global_earth_coordinate, local_mesh_coordinate

  implicit none

  !*** Data gather declaration (ZNE, XYZ, ZRT, LQT)
  real(kind=sp), dimension(:,:), allocatable :: vx, vy, vz
  real(kind=sp), dimension(:,:), allocatable :: vn, ve
  real(kind=sp), dimension(:,:), allocatable :: vr, vt
  real(kind=sp), dimension(:,:), allocatable :: vl, vq

  !*** Data vector declaration (time, baz, inc,, gcarc and dist)
  real(kind=dp), dimension(:), allocatable :: tt
  real(kind=dp), dimension(:), allocatable :: baz
  real(kind=dp), dimension(:), allocatable :: inc
  real(kind=dp), dimension(:), allocatable :: gcarc
  real(kind=dp), dimension(:), allocatable :: dist
  real(kind=dp), dimension(:), allocatable :: stf

  !*** Some parameters for IO
  integer(kind=si),   private :: iunit=34, io_err, debug_level=0
  character(len=256), private :: line, keyword, keyval

  !*** Define type for station file header
  type hdr_type

     character(len=256) :: event_name        = 'undef'  ! name of event
     character(len=256) :: source_type       = 'undef'  ! type of source 'moment' or 'force'
     character(len=256) :: source_components = 'undef'  ! file with parameters of source
     character(len=256) :: modeling_tool     = 'undef'  ! kind of modeling (pointsource,
                                                        !    finitesource, Axisem, FK, PW)
     character(len=256) :: modeling_path     = 'undef'  ! repository of tractions or parameter file
                                                        !    or STF file for point source
     character(len=256) :: estimated_src     = 'undef'  ! repository of tractions or parameter file
     character(len=1)   :: data_type         = 'v'      ! kind of data (d=displacement, !
                                                        ! v=velocities, a=acceleration)
     character(len=3)   :: data_comp         = 'enz'    ! data coordinate system (geo 'zen'),
                                                        !    (local 'zxy'), (rotated 'zrt' or 'lqt')
     logical, dimension(3) :: is_comp        = (/ .true., .true., .true. /)  ! actual components available
     character(len=3)   :: coord_sys         = 'geo'    ! stations coordinate system 'geo' or 'car'

     real(kind=dp), dimension(3) :: mesh_origin = (/0._dp,0._dp,0._dp/) ! mesh top center geographic coordinates
              !(lat0, lon0, azi) (note: azi taken counterclockwise from east)

     real(kind=dp) :: otime  = 0._dp    ! event origin time wrt first sample
     real(kind=dp) :: tb     = 0._dp    ! time of data first sample
     real(kind=dp) :: te     = 0._dp    ! time of data last sample
     real(kind=dp) :: tl     = 0._dp    ! total time window
     real(kind=dp) :: dt     = 0._dp    ! time step
     real(kind=dp) :: fs     = 0._dp    ! sampling frequency

     integer(kind=si)  :: nt     = 0    ! number of time steps
     integer(kind=si)  :: nsta   = 0    ! number of stations

     logical :: is_pick     = .false.   ! flags : is picks
     logical :: is_window   = .false.   ! flags : is window

     real(kind=dp) :: tbef = 0._dp      ! time before pick for time window
     real(kind=dp) :: taft = 0._dp      ! time after  pick for time window

  end type hdr_type

  !*** Define type for one station information
  type station_type

     character(len=5) :: name = 'undef' ! name of station
     character(len=5) :: ntwk = 'undef' ! name of network

     real(kind=dp) :: x=0._dp,     y=0._dp,   z=0._dp    ! local Cartesian coordinates
     real(kind=dp) :: lat=0._dp, lon=0._dp, ele=0._dp    ! geographic coordinates

     real(kind=dp) :: tpick    ! time pick
     real(kind=dp) :: baz      ! baz
     real(kind=dp) :: gcarc    ! epicentral distance
     real(kind=dp) :: dist     ! distance to event

  end type station_type

  !*** Define type for one source information
  type source_type

     real(kind=dp), dimension(3) :: f = 0._dp   ! force vector  : fx, fy, fz
     real(kind=dp), dimension(6) :: m = 0._dp   ! moment tensor : mrr, mtt, mpp, mtp, mtr, mrt
     ! todo: rotate tensor to xyz local
     !  and add posibility o read it in local mesh

     real(kind=dp) :: inc=0._dp, baz=0._dp      ! event incidence and back-azimuth (for FK and PW)

     real(kind=dp) :: x=0._dp,     y=0._dp,   z=0._dp  ! local Cartesian coordinates
     real(kind=dp) :: lat=0._dp, lon=0._dp, ele=0._dp  ! local geographic coordinates

     integer(kind=si) :: yy=0,     mo=0,     dd=0, jday=0  ! event datetime : year, month, day, julian day
     real(kind=dp)    :: hh=0._dp, mi=0._dp, ss=0._dp      !                  hour, minute, decimal sec

     character(len=254) :: hdr    ! header of cmtfile
     character(len=254) :: name   ! name of event in cmt file

     real(kind=dp) :: tshift      ! time shift from cmtsolution
     real(kind=dp) :: hdur        ! half duration from cmtsolution

  end type source_type

  !*** Define type for the gather
  type gather

     type(hdr_type)    :: hdr
     type(source_type) :: source
     type(station_type), dimension(:), allocatable :: stations

  end type gather


contains

  !================================================================================
  ! Read passive imaging format gather file
  subroutine read_pif_header_file(filename,mygather)

    character(len=*), intent(in)  :: filename
    type(gather),     intent(out) :: mygather
    integer(kind=si) :: k

    integer :: ier

    write(*,*)'Read PIF-file header from ',trim(adjustl(filename)),' ...'
    open(iunit, file=trim(adjustl(filename)), status='old',action='read', iostat=io_err)

    if (io_err /= 0) then
       write(*,*)'PIF file: ',trim(adjustl(filename)),' does not exist!'
       stop 'PIF file does not exist!'
    endif

    do

       read(iunit,fmt='(a256)',iostat=io_err) line
       if (io_err < 0) exit
       if (trim(adjustl(line)) == 'STARTLIST') exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

       read(line,*) keyword, keyval

       !*** Ensure lower case
       keyword = lowcase(keyword)

       select case(trim(keyword))
       case('event_name')
          read(keyval,*) mygather%hdr%event_name
          write(*,*)'    event_name ',trim(adjustl(mygather%hdr%event_name))
       case('source_type')
          keyval  = lowcase(keyval)
          read(keyval,*) mygather%hdr%source_type
          write(*,*)'    source_type ',trim(adjustl(mygather%hdr%source_type))
       case('source_components')
          read(line,*)keyval, mygather%hdr%source_components
          write(*,*)'    receiver_component ',trim(adjustl(mygather%hdr%source_components))
       case('modeling_tool')
          read(line,*) keyval, mygather%hdr%modeling_tool, mygather%hdr%modeling_path
          mygather%hdr%modeling_tool  = lowcase(mygather%hdr%modeling_tool)
          write(*,*)'    modeling_tool ',trim(adjustl(mygather%hdr%modeling_tool)),' ', &
               trim(adjustl(mygather%hdr%modeling_path))
       case('data_components')
          read(line,*) keyval, mygather%hdr%data_type, mygather%hdr%data_comp
          mygather%hdr%data_type = lowcase(mygather%hdr%data_type)
          mygather%hdr%data_comp = lowcase(mygather%hdr%data_comp)
          write(*,*)'    data_components ',trim(adjustl(mygather%hdr%data_type)),' ', &
               trim(adjustl(mygather%hdr%data_comp))
       case('cartloc_mesh_origin')
          read(line,*) keyword, mygather%hdr%mesh_origin(1:3)
          write(*,*)'    cartloc_mesh_origin :',mygather%hdr%mesh_origin
       case('data_origin_time')
          read(keyval,*) mygather%hdr%otime
          write(*,*)'    data_origin_time ',mygather%hdr%otime
       case('number_of_station')
          read(keyval,*) mygather%hdr%nsta
          write(*,*)'    number_of_station ',mygather%hdr%nsta
       case('data_time_step')
          read(keyval,*) mygather%hdr%dt
          write(*,*)'    data_time_step ',mygather%hdr%dt
       case('data_sample_number')
          read(keyval,*) mygather%hdr%nt
          write(*,*)'    data_sample_number ',mygather%hdr%nt
       case('is_time_pick')
          read(keyval,*) mygather%hdr%is_pick
          write(*,*)'    is_time_pick ',mygather%hdr%is_pick
       case('time_window')
          read(line,*)      keyval, mygather%hdr%is_window, mygather%hdr%tbef, mygather%hdr%taft
          write(*,*)'    time_window ',mygather%hdr%is_window,mygather%hdr%tbef,mygather%hdr%taft
       case('station_coord_system')
          keyval = lowcase(keyval)
          read(keyval,*) mygather%hdr%coord_sys
          write(*,*)'    station_coord_system ',trim(adjustl(mygather%hdr%coord_sys))
       end select

    enddo

    !*** Get mesh matrix rotation
    call define_mesh_rotation_matrix(mygather%hdr%mesh_origin(1), &
                                     mygather%hdr%mesh_origin(2), &
                                     mygather%hdr%mesh_origin(3))

    !*** Now read stations
    k = 0
    if (.not. allocated(mygather%stations)) then
      allocate(mygather%stations(mygather%hdr%nsta),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 489')
    endif
    do
       read(iunit,fmt='(a256)',iostat=io_err) line
       if (trim(adjustl(line)) == 'ENDLIST') exit
       if (io_err < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

       k = k + 1
       if (k > mygather%hdr%nsta) stop 'Too much stations in PIF binary file..., check number_of_station'

       select case (mygather%hdr%coord_sys)
       case('car')

          if (mygather%hdr%is_pick) then
             read(line,*)mygather%stations(k)%name, &
                         mygather%stations(k)%ntwk, &
                         mygather%stations(k)%x, &
                         mygather%stations(k)%y, &
                         mygather%stations(k)%z, &
                         mygather%stations(k)%tpick
          else
             read(line,*)mygather%stations(k)%name, &
                         mygather%stations(k)%ntwk, &
                         mygather%stations(k)%x, &
                         mygather%stations(k)%y, &
                         mygather%stations(k)%z
          endif
          call global_earth_coordinate(mygather%stations(k)%x, &
                                       mygather%stations(k)%y, &
                                       mygather%stations(k)%z, &
                                       mygather%stations(k)%lat, &
                                       mygather%stations(k)%lon, &
                                       mygather%stations(k)%ele)
          if (debug_level > 2) then
             print *,k,mygather%stations(k)%x,mygather%stations(k)%y,mygather%stations(k)%z
             print *,k,mygather%stations(k)%lat,mygather%stations(k)%lon,mygather%stations(k)%ele
          endif

       case ('geo')

          if (mygather%hdr%is_pick) then
             read(line,*)mygather%stations(k)%name, &
                         mygather%stations(k)%ntwk, &
                         mygather%stations(k)%lat, &
                         mygather%stations(k)%lon, &
                         mygather%stations(k)%ele, &
                         mygather%stations(k)%tpick
          else
             read(line,*)mygather%stations(k)%name, &
                         mygather%stations(k)%ntwk, &
                         mygather%stations(k)%lat, &
                         mygather%stations(k)%lon, &
                         mygather%stations(k)%ele
          endif
          call local_mesh_coordinate(mygather%stations(k)%lat, &
                                     mygather%stations(k)%lon, &
                                     mygather%stations(k)%ele, &
                                     mygather%stations(k)%x, &
                                     mygather%stations(k)%y, &
                                     mygather%stations(k)%z)
          if (debug_level > 2) then
             print *,k,mygather%stations(k)%x,mygather%stations(k)%y,mygather%stations(k)%z
             print *,k,mygather%stations(k)%lat,mygather%stations(k)%lon,mygather%stations(k)%ele
          endif

       end select

       if (debug_level > 1) write(*,*)trim(adjustl(line))

    enddo

    if (k < mygather%hdr%nsta) stop 'Not enough stations in PIF binary file..., check number_of_station'
    close(iunit)
    write(*,*)'Done!'

    ! Read CMT file if needed
    select case (mygather%hdr%source_type)
    case('moment')
       call read_cmt_solution_file(mygather%hdr%source_components, mygather%source)
       ! Just in case, determine local coordinate (useful for local point source)
       call local_mesh_coordinate(mygather%source%lat, &
                                  mygather%source%lon, &
                                  mygather%source%ele, &
                                  mygather%source%x, &
                                  mygather%source%y, &
                                  mygather%source%z)
    case('force')
       call read_force_solution_file(mygather%hdr%source_components, mygather%source)
       ! Just in case, determine local coordinate (useful for local point source)
       call local_mesh_coordinate(mygather%source%lat, &
                                  mygather%source%lon, &
                                  mygather%source%ele, &
                                  mygather%source%x, &
                                  mygather%source%y, &
                                  mygather%source%z)
    case default
       write(*,*) 'WARNING : source_type undefined !'
    end select

!!!! Not here.... since only header
    ! Read source estimated time function
    !if (mygather%hdr%estimated_stf /= 'undef') then
    !    call read_binary_source_signature(filename,nt,stf)
    !endif

  end subroutine read_pif_header_file
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Read CMT solution file
  subroutine read_cmt_solution_file(filename,cmt)

    character(len=*),    intent(in) :: filename
    type(source_type),  intent(out) :: cmt

    write(*,*)'Read CMT solution file ...'
    print *,filename
    open(iunit, file=trim(adjustl(filename)), status='old',action='read', iostat=io_err)
    if (io_err /= 0) then
       write(*,*)'CMT solution file: ',trim(adjustl(filename)),' does not exist!'
       stop 'CMT solution file does not exist!'
    endif

    do
       read(iunit,fmt='(a256)',iostat=io_err) line
       if (io_err < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

       !*** Ensure lower case
       line    = lowcase(line)
       keyword = lowcase(keyword)
       keyval  = lowcase(keyval)

       read(line,*) keyword, keyval
       if (keyword(2:4) == 'pde') then
          cmt%hdr = trim(adjustl(line))
          ! to do : fill cmt header
          read(line(6:9),*)cmt%yy
          read(line(11:12),*)cmt%mo
          read(line(14:15),*)cmt%dd
          read(line(17:18),*)cmt%hh
          read(line(20:21),*)cmt%mi
          read(line(23:27),*)cmt%ss
          cycle
       endif

       select case(trim(keyword))
       case('event')
          read(line,*) keyword, cmt%name
       case('time')
          read(line,*) keyword, keyval, cmt%tshift
       case('half')
          read(line,*) keyword, keyval, cmt%hdur
       case('latitude:','latorUTM:')
          read(line,*) keyword, cmt%lat
       case('longitude:','lonorUTM:')
          read(line,*) keyword, cmt%lon
       case('depth:')
          read(line,*) keyword, cmt%ele
       case('mrr:','Mrr:')
          read(line,*) keyword, cmt%m(1)
       case('mtt:','Mtt:')
          read(line,*) keyword, cmt%m(2)
       case('mpp:','Mpp:')
          read(line,*) keyword, cmt%m(3)
       case('mrt:','Mrt:')
          read(line,*) keyword, cmt%m(6)
       case('mrp:','Mrp:')
          read(line,*) keyword, cmt%m(5)
       case('mtp:','Mtp:')
          read(line,*) keyword, cmt%m(4)
       end select

    enddo

    close(iunit)
    write(*,*)'Done!'

  end subroutine read_cmt_solution_file
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Read force solution file
  subroutine read_force_solution_file(filename,cmt)

    character(len=*),    intent(in) :: filename
    type(source_type),  intent(out) :: cmt

    write(*,*)'Read force solution file ...'
    open(iunit, file=trim(adjustl(filename)), status='old',action='read', iostat=io_err)
    if (io_err /= 0) then
       write(*,*)'Force solution file: ',trim(adjustl(filename)),' does not exist!'
       stop 'Force solution file does not exist!'
    endif

    do
       read(iunit,fmt='(a256)',iostat=io_err) line
       if (io_err < 0) exit
       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

       !*** Ensure lower case
       line    = lowcase(line)
       keyword = lowcase(keyword)
       keyval  = lowcase(keyval)

       read(line,*) keyword, keyval

       select case(trim(keyword))
       case('force_components')
          read(line,*) keyword, cmt%f(:)
       case('force_lat')
          read(line,*) keyword, cmt%lat
       case('force_lon')
          read(line,*) keyword, cmt%lon
       case('force_depth')
          read(line,*) keyword, cmt%ele
       case('force_year')
          read(line,*) keyword, cmt%yy
       case('force_month')
          read(line,*) keyword, cmt%mo
       case('force_day')
          read(line,*) keyword, cmt%dd
       case('force_hh')
          read(line,*) keyword, cmt%hh
       case('force_mm')
          read(line,*) keyword, cmt%mi
       case('froce_dsec')
          read(line,*) keyword, cmt%ss
       end select

    enddo

    close(iunit)
    write(*,*)'Done!'

  end subroutine read_force_solution_file
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Read one binary gather
  subroutine read_binary_data(filename,nrec,nt,data)

    character(len=*),                     intent(in)  :: filename
    integer(kind=si),                     intent(in)  :: nrec, nt
    integer(kind=si)                                  :: nsize, it
    real(kind=sp),    dimension(nrec,nt)              :: datas
    real(kind=cp),    dimension(nrec,nt), intent(out) :: data

    write(*,*)'Read binary data: '
    write(*,*)'    filename, nrec, nt = ',trim(adjustl(filename)),nrec,nt

    nsize = nrec * sp
    open(iunit,file=trim(adjustl(filename)),access='direct',recl=nsize,status='old')
    ! use a loop here because when large number of receiver integer overflow may happen
    do it = 1, nt
       read(iunit,rec=it)datas(:,it)
    enddo
    close(iunit)
    data = real(datas,kind=cp)

    write(*,*)'Done!'

  end subroutine read_binary_data
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Read one binary gather
  subroutine read_binary_source_signature(filename,nt,stf)

    character(len=*),                     intent(in)  :: filename
    integer(kind=si),                     intent(in)  :: nt
    integer(kind=si)                                  :: nsize
    real(kind=sp),    dimension(nt)              :: stfs
    real(kind=cp),    dimension(nt), intent(out) :: stf

    write(*,*)'Read binary source signature: '
    write(*,*)'    filename, nt = ',trim(adjustl(filename)),nt

    nsize = nt * sp

    open(iunit,file=trim(adjustl(filename)),access='direct',recl=nsize,status='old')
    read(iunit,rec=1)stfs(:)
    close(iunit)

    stf = real(stfs,kind=cp)

    write(*,*)'Done!'

  end subroutine read_binary_source_signature
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Read one binary gather
  subroutine write_binary_source_signature(filename,nt,stf)

    character(len=*),                     intent(in)  :: filename
    integer(kind=si),                     intent(in)  :: nt
    integer(kind=si)                                  :: nsize
    real(kind=sp),    dimension(nt)              :: stfs
    real(kind=cp),    dimension(nt), intent(in)  :: stf

    write(*,*)'Write binary source signature: '
    write(*,*)'    filename, nt = ',trim(adjustl(filename)),nt

    nsize = nt * sp
    stfs  = real(stf,kind=sp)

    open(iunit,file=trim(adjustl(filename)),access='direct',recl=nsize,status='replace')
    write(iunit,rec=1)stfs(:)
    close(iunit)

    write(*,*)'Done!'

  end subroutine write_binary_source_signature
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Write one binary gather
  subroutine write_binary_data(filename,nrec,nt,data)

    character(len=*),                     intent(in) :: filename
    integer(kind=si),                     intent(in) :: nrec, nt
    integer(kind=si)                                 :: nsize, it
    real(kind=cp), dimension(nrec,nt), intent(in) :: data

    write(*,*)'Write binary data: '
    write(*,*)'    filename, nrec, nt = ',trim(adjustl(filename)),nrec,nt

    nsize = nrec * sp
    open(iunit,file=trim(adjustl(filename)),access='direct',recl=nsize,status='replace')
    do it = 1, nt
       write(iunit,rec=it)real(data(:,it),kind=sp)
    enddo
    close(iunit)

    write(*,*)'Done!'

  end subroutine write_binary_data
  !--------------------------------------------------------------------------------

  ! ================================================================================
  ! Compute epicentral distance, baz and dit from stations/source informations
  !                   for now, only a spherical Earth of radius = 6371 km
  subroutine calc_delta_dist_baz(slati,sloni,rlati,rloni,gcarc,dist,baz)

    real(kind=dp), intent(in)  :: slati, sloni, rlati, rloni
    real(kind=dp), intent(out) :: gcarc, dist, baz

    real(kind=dp) :: slat, slon, rlat, rlon, den, num, cosarc
    real(kind=dp) :: r=6371000._dp

    !*** Pass in radian
    slat = deg2rad * slati
    slon = deg2rad * sloni
    rlat = deg2rad * rlati
    rlon = deg2rad * rloni

    !*** Compute distance in degree and in km
    cosarc = sin(rlat)*sin(slat) + cos(rlat)*cos(slat)*cos(rlon-slon)
    gcarc  = acos(cosarc)
    dist   = r * gcarc

    !*** Compute back-azimuth
    num = sin(rlon-slon) * cos(slat) * cos(rlat)
    den =     sin(slat) - cos(gcarc) * sin(rlat)
    baz= -(rad2deg * atan2(num,den)) - 360._dp
    if (baz >= 360._dp) then
       baz = baz -360._dp
    else if (baz < 0._dp) then
       baz = baz +360._dp
    else
       baz = baz
    endif

  end subroutine calc_delta_dist_baz
  !--------------------------------------------------------------------------------

  ! ================================================================================
  ! Compute offset and baz from stations/source informations
  subroutine calc_dist_baz_cart(src_x,src_y,sta_x,sta_y,dist,baz)

    real(kind=dp), intent(in)  :: src_x, src_y, sta_x, sta_y
    real(kind=dp), intent(out) :: dist, baz

    real(kind=dp) :: dx, dy

    !*** Compute distance in degree and in km
    dx   = src_x-sta_x
    dy   = src_y-sta_y
    dist = sqrt(dx*dx + dy*dy)

    !*** Compute back-azimuth
    baz  = rad2deg * atan2(dx,dy)
    if (baz >= 360._dp) then
       baz = baz -360._dp
    else if (baz < 0._dp) then
       baz = baz +360._dp
    else
       baz = baz
    endif

  end subroutine calc_dist_baz_cart
  !--------------------------------------------------------------------------------


  !================================================================================
  ! Make sur that everything is upper case...
  function upcase(strin) result(strout)

    character(len=*), intent(in) :: strin
    character(len=len(strin))    :: strout
    integer(kind=si)             :: i, j

    do i = 1,len(strin)
       j = iachar(strin(i:i))
       if (j >= iachar('a') .and. j <= iachar('z')) then
          strout(i:i) = achar(iachar(strin(i:i))-32)
       else
          strout(i:i) = strin(i:i)
       endif
    enddo

  end function upcase
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Make sur that everything is lower case...
  function lowcase(strin) result(strout)

    character(len=*), intent(in) :: strin
    character(len=len(strin))    :: strout
    integer(kind=si)             :: i, j

    do i = 1,len(strin)
       j = iachar(strin(i:i))
       if (j >= iachar('A') .and. j <= iachar('Z')) then
          strout(i:i) = achar(iachar(strin(i:i))+32)
       else
          strout(i:i) = strin(i:i)
       endif
    enddo

  end function lowcase
  !--------------------------------------------------------------------------------

  !================================================================================
  ! Give data comp in XX format (ex: vz)
  subroutine get_data_component(data_type, data_comp, ncomp, component)

    character(len=1), intent(in) :: data_type
    character(len=3), intent(in) :: data_comp

    character(len=2), dimension(3), intent(out) :: component
    integer(kind=si),               intent(out) :: ncomp

    integer(kind=si) :: i

    component(:) = '  '
    ncomp        = len(trim(adjustl(data_comp)))

    do i = 1, ncomp
       write(component(i),*)upcase(data_type),upcase(data_comp(i:i))
    enddo

  end subroutine get_data_component
  !--------------------------------------------------------------------------------


end module passive_imaging_format_mod
