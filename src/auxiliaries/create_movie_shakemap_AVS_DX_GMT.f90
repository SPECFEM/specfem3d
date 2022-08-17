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

!
!---  create a movie of the vertical component of surface displacement or velocity
!---  or a ShakeMap(R) (i.e. map of the maximum absolute value of the two horizontal components
!---  of the velocity vector) in AVS, OpenDX or GMT format
!
!
! AVS UCD format descriptions:
! https://lanl.github.io/LaGriT/pages/docs/read_avs.html
! http://people.sc.fsu.edu/~jburkardt/data/ucd/ucd.html
! http://www.hnware.de/rismo/dokumente/anwenderdoku/formate/avs_ucd.html

  program create_movie_shakemap

  use constants
  use shared_parameters

  implicit none

!-------------------------------------------------------------------------------------------------
! user parameters

! normalizes field display values
  logical, parameter :: NORMALIZE_OUTPUT = .false.

! threshold in percent of the maximum below which we cut the amplitude
  logical, parameter :: APPLY_THRESHOLD = .false.
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

! coefficient of power law used for non linear scaling
  logical, parameter :: NONLINEAR_SCALING = .false.
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.13_CUSTOM_REAL

! muting source region
  logical, parameter :: MUTE_SOURCE = .false.
  real(kind=CUSTOM_REAL), parameter :: RADIUS_TO_MUTE = 1000._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: X_SOURCE_EXT_MESH = -9023.021484375
  real(kind=CUSTOM_REAL), parameter :: Y_SOURCE_EXT_MESH = 6123.611328125
  real(kind=CUSTOM_REAL), parameter :: Z_SOURCE_EXT_MESH = 17.96331405639648

!-------------------------------------------------------------------------------------------------

  integer :: it,it1,it2,ivalue,nspectot_AVS_max,ispec
  integer :: iformat,nframes,iframe,inumber,inorm,iscaling_shake
  integer :: ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4

  logical :: USE_OPENDX,USE_AVS,USE_GMT,plot_shaking_map

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,display
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) :: vectorx,vectory,vectorz,vectornorm

  double precision :: min_field_current,max_field_current,max_absol

  character(len=MAX_STRING_LEN) :: outputname
  character(len=MAX_STRING_LEN) :: line

  integer :: ipoin

  ! GMT
  double precision :: lat,long

  ! for sorting routine
  integer :: npointot,ilocnum,nglob,i,j,ielm,ieoff,ispecloc
  integer, dimension(:), allocatable :: iglob,locval,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

  ! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

  integer :: ier

  ! order of points representing the 2D square element
  integer,dimension(NGNOD2D_FOUR_CORNERS_AVS_DX),parameter :: iorder = (/1,3,2,4/)

  integer :: NSPEC_SURFACE_EXT_MESH
  logical :: BROADCAST_AFTER_READ

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *

  print *
  print *,'reading parameter file'
  print *

  ! initializes
  myrank = 0
  BROADCAST_AFTER_READ = .false.

  ! read the parameter file
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! only one global array for movie data, but stored for all surfaces defined
  ! in file 'surface_from_mesher.h'
  open(unit=IIN,file=trim(OUTPUT_FILES)//'surface_from_mesher.h',status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'error opening file: ',trim(OUTPUT_FILES)//'surface_from_mesher.h'
    print *
    print *,'please run xgenerate_databases or xspecfem3D first to create this file, exiting now...'
    stop 'error opening moviedata header file'
  endif
  ! skips first few lines
  do i=1,6
    read(IIN,'(a)') line
  enddo
  ! line with info, e.g. "integer,parameter :: NSPEC_SURFACE_EXT_MESH = 23855"
  read(IIN,'(a)') line
  close(IIN)
  ! gets number from substring after = sign
  i = index(line,'=')
  if (i == 0) stop 'error reading in NSPEC_SURFACE_EXT_MESH from file OUTPUT_FILES/surface_from_mesher.h'

  read(line(i+1:len_trim(line)),*) NSPEC_SURFACE_EXT_MESH

  ! calculates number of total surface points
  if (USE_HIGHRES_FOR_MOVIES) then
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGLLSQUARE
  else
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGNOD2D_FOUR_CORNERS_AVS_DX
  endif
  print *,'  high-resolution: ',USE_HIGHRES_FOR_MOVIES
  print *,'  moviedata element surfaces: ',NSPEC_SURFACE_EXT_MESH
  print *,'  moviedata total elements all: ',ilocnum
  print *

  if (SAVE_DISPLACEMENT) then
    print *,'Displacement will be shown in movie'
  else
    print *,'Velocity will be shown in movie'
  endif
  print *

  if (MUTE_SOURCE) then
    print *,'Muting source region:'
    print *,'  radius = ',RADIUS_TO_MUTE
    print *,'  source location x/y/z = ',X_SOURCE_EXT_MESH,Y_SOURCE_EXT_MESH,Z_SOURCE_EXT_MESH
    print *
  endif

  ! user input
  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format'
  print *,'3 = create files in GMT xyz Ascii long/lat/Uz format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if (iformat < 1 .or. iformat > 3) stop 'exiting...'

  plot_shaking_map = .false.
  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *
  print *,'enter first time step of movie (e.g. 1, enter -1 for shaking map)'
  read(5,*) it1

  if (it1 < 1 .and. it1 /= -1) stop 'the first time step must be >= 1'

  if (it1 == -1) plot_shaking_map = .true.

  if (.not. plot_shaking_map) then
    print *,'enter last time step of movie (e.g. ',NSTEP,')'

    read(5,*) it2
    ! limits to maximum of NSTEP
    if (it2 > NSTEP) then
      it2 = NSTEP
    endif

    print *
    print *,'1 = define file names using frame number'
    print *,'2 = define file names using time step number'
    print *,'any other value = exit'
    print *
    print *,'enter value:'

    read(5,*) inumber
    if (inumber < 1 .or. inumber > 2) stop 'exiting...'

    print *
    print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'
    ! count number of movie frames
    nframes = 0
    do it = it1,it2
      if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
    enddo
  else
    ! only one frame if shaking map
    nframes = 1
    it1 = 1
    it2 = 1
  endif
  print *
  print *,'total number of frames will be ',nframes
  if (nframes == 0) stop 'null number of frames'

  iscaling_shake = 0
  if (plot_shaking_map) then
    print *
    print *,'norm to display in shaking map:'
    print *,'1=displacement  2=velocity  3=acceleration'
    print *
    read(5,*) inorm
    if (inorm < 1 .or. inorm > 3) stop 'incorrect value of inorm'
    print *
    print *,'apply non-linear scaling to shaking map:'
    print *,'1=non-linear  2=no scaling'
    print *
    read(5,*) iscaling_shake
    if (iscaling_shake < 1 .or. iscaling_shake > 2) stop 'incorrect value of iscaling_shake'
  else
    print *
    print *,'movie data:'
    if (SAVE_DISPLACEMENT) then
      print *,'1= norm of displacement  2=displacement x-comp 3=displacement y-comp 4=displacement z-comp'
    else
      print *,'1= norm of velocity  2=velocity x-comp 3=velocity y-comp 4=velocity z-comp'
    endif
    print *
    read(5,*) inorm
    if (inorm < 1 .or. inorm > 4) stop 'incorrect value of inorm'
  endif

! file format flags
  if (iformat == 1) then
    USE_OPENDX = .true.
    USE_AVS = .false.
    USE_GMT = .false.
  else if (iformat == 2) then
    USE_OPENDX = .false.
    USE_AVS = .true.
    USE_GMT = .false.
  else
    USE_OPENDX = .false.
    USE_AVS = .false.
    USE_GMT = .true.
  endif

  ! define the total number of elements at the surface
  if (USE_HIGHRES_FOR_MOVIES) then
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH * (NGLLX-1) * (NGLLY-1)
  else
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH
  endif

  ! maximum theoretical number of points at the surface
  npointot = NGNOD2D_FOUR_CORNERS_AVS_DX * nspectot_AVS_max

  ! allocate arrays for sorting routine
  allocate(iglob(npointot), &
           locval(npointot), &
           ifseg(npointot), &
           xp(npointot), &
           yp(npointot), &
           zp(npointot), &
           xp_save(npointot), &
           yp_save(npointot), &
           zp_save(npointot), &
           field_display(npointot), &
           mask_point(npointot), &
           ireorder(npointot),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for sorting routine'

  ! allocates data arrays
  allocate(store_val_x(ilocnum), &
           store_val_y(ilocnum), &
           store_val_z(ilocnum), &
           store_val_ux(ilocnum), &
           store_val_uy(ilocnum), &
           store_val_uz(ilocnum),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays for data arrays'

  if (USE_HIGHRES_FOR_MOVIES) then
    allocate(x(NGLLX,NGLLY), &
             y(NGLLX,NGLLY), &
             z(NGLLX,NGLLY), &
             display(NGLLX,NGLLY),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays for highres'
  endif

  ! user output
  print *
  print *,'there are a total of ',nspectot_AVS_max,' elements at the surface'
  print *
  print *
  if (APPLY_THRESHOLD .and. .not. plot_shaking_map) &
    print *,'Will apply a threshold to amplitude below ',100.*THRESHOLD,' %'
  if (NONLINEAR_SCALING .and. (.not. plot_shaking_map .or. iscaling_shake == 1)) &
    print *,'Will apply a non linear scaling with coef ',POWER_SCALING

  iframe = 0

! loop on all the time steps in the range entered
  do it = it1,it2

    ! check if time step corresponds to a movie frame
    if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0 .or. plot_shaking_map) then

      iframe = iframe + 1

      print *
      if (plot_shaking_map) then
        print *,'reading shaking map snapshot'
      else
        print *,'reading snapshot time step ',it,' out of ',NSTEP
      endif
      print *

      ! read all the elements from the same file
      if (plot_shaking_map) then
        write(outputname,"('/shakingdata')")
      else
        write(outputname,"('/moviedata',i6.6)") it
      endif
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname),status='old', &
            action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'error: ',trim(OUTPUT_FILES)//trim(outputname)
        stop 'error opening moviedata file'
      endif

      read(IOUT) store_val_x
      read(IOUT) store_val_y
      read(IOUT) store_val_z
      read(IOUT) store_val_ux
      read(IOUT) store_val_uy
      read(IOUT) store_val_uz
      close(IOUT)

      ! clear number of elements kept
      ispec = 0

      ! reset point number
      ipoin = 0

      do ispecloc = 1, NSPEC_SURFACE_EXT_MESH

        if (USE_HIGHRES_FOR_MOVIES) then
          ! assign the OpenDX "elements"
          do j = 1,NGLLY
            do i = 1,NGLLX
              ipoin = ipoin + 1

              ! x,y,z coordinates
              xcoord = store_val_x(ipoin)
              ycoord = store_val_y(ipoin)
              zcoord = store_val_z(ipoin)

              ! note:
              ! for shakemaps: ux = norm displacement, uy = norm velocity, uz = norm acceleration
              ! for movies: ux = velocity x-component, uy = velocity y-component, uz = velocity z-component
              vectorx = store_val_ux(ipoin)
              vectory = store_val_uy(ipoin)
              vectorz = store_val_uz(ipoin)

              x(i,j) = xcoord
              y(i,j) = ycoord
              z(i,j) = zcoord

              ! shakemap
              if (plot_shaking_map) then
                ! chooses norm
                if (inorm == 1) then
                  ! norm displacement
                  display(i,j) = vectorx
                else if (inorm == 2) then
                  ! norm velocity
                  display(i,j) = vectory
                else
                  ! norm acceleration
                  display(i,j) = vectorz
                endif
                !!!! NL NL mute value near source
                if (MUTE_SOURCE) then
                  if ( (sqrt(((x(i,j) - (X_SOURCE_EXT_MESH))**2 + &
                              (y(i,j) - (Y_SOURCE_EXT_MESH))**2 + &
                              (z(i,j) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) ) then
                    display(i,j) = 0.0
                  endif
                endif
              else
                ! movie
                if (inorm == 1) then
                  ! norm of velocity
                  vectornorm = sqrt(vectorz*vectorz + vectory*vectory + vectorx*vectorx)
                  display(i,j) = vectornorm
                else if (inorm == 2) then
                  ! velocity x-component
                  display(i,j) = vectorx
                else if (inorm == 3) then
                  ! velocity y-component
                  display(i,j) = vectory
                else if (inorm == 4) then
                  ! velocity z-component
                  display(i,j) = vectorz
                endif
              endif

            enddo
          enddo

          ! assign the values of the corners of the OpenDX "elements"
          ispec = ispec + 1
          ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)

          do j = 1,NGLLY-1
            do i = 1,NGLLX-1
              ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ielm+(i-1)+(j-1)*(NGLLX-1))
              do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX

                if (ilocnum == 1) then
                  xp(ieoff+ilocnum) = dble(x(i,j))
                  yp(ieoff+ilocnum) = dble(y(i,j))
                  zp(ieoff+ilocnum) = dble(z(i,j))
                  field_display(ieoff+ilocnum) = dble(display(i,j))
                else if (ilocnum == 2) then

                  ! accounts for different ordering of square points
                  xp(ieoff+ilocnum) = dble(x(i+1,j+1))
                  yp(ieoff+ilocnum) = dble(y(i+1,j+1))
                  zp(ieoff+ilocnum) = dble(z(i+1,j+1))
                  field_display(ieoff+ilocnum) = dble(display(i+1,j+1))

                else if (ilocnum == 3) then

                  ! accounts for different ordering of square points
                  xp(ieoff+ilocnum) = dble(x(i+1,j))
                  yp(ieoff+ilocnum) = dble(y(i+1,j))
                  zp(ieoff+ilocnum) = dble(z(i+1,j))
                  field_display(ieoff+ilocnum) = dble(display(i+1,j))

                else
                  xp(ieoff+ilocnum) = dble(x(i,j+1))
                  yp(ieoff+ilocnum) = dble(y(i,j+1))
                  zp(ieoff+ilocnum) = dble(z(i,j+1))
                  field_display(ieoff+ilocnum) = dble(display(i,j+1))
                endif

              enddo

            enddo
          enddo

        else
          ! low-resolution (only spectral element corners)
          ispec = ispec + 1
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)

          ! four points for each element
          do i = 1,NGNOD2D_FOUR_CORNERS_AVS_DX

            ! accounts for different ordering of square points
            ilocnum = iorder(i)

            ipoin = ipoin + 1

            xcoord = store_val_x(ipoin)
            ycoord = store_val_y(ipoin)
            zcoord = store_val_z(ipoin)

            vectorx = store_val_ux(ipoin)
            vectory = store_val_uy(ipoin)
            vectorz = store_val_uz(ipoin)


            xp(ilocnum+ieoff) = dble(xcoord)
            yp(ilocnum+ieoff) = dble(ycoord)
            zp(ilocnum+ieoff) = dble(zcoord)

            ! shakemap
            if (plot_shaking_map) then
              if (inorm == 1) then
                ! norm of displacement
                field_display(ilocnum+ieoff) = dble(vectorx)
              else if (inorm == 2) then
                ! norm of velocity
                field_display(ilocnum+ieoff) = dble(vectory)
              else
                ! norm of acceleration
                field_display(ilocnum+ieoff) = dble(vectorz)
              endif
              !!!! NL NL mute value near source
              if (MUTE_SOURCE) then
                if (sqrt(((dble(xcoord) - (X_SOURCE_EXT_MESH))**2 + &
                          (dble(ycoord) - (Y_SOURCE_EXT_MESH))**2 + &
                          (dble(zcoord) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) then
                  field_display(ilocnum+ieoff) = 0.0
                endif
              endif
            else
              ! movie
              if (inorm == 1) then
                ! norm of velocity
                field_display(ilocnum+ieoff) = sqrt(vectorz**2+vectory**2+vectorx**2)
              else if (inorm == 2) then
                ! velocity x-component
                field_display(ilocnum+ieoff) = vectorx
              else if (inorm == 3) then
                ! velocity y-component
                field_display(ilocnum+ieoff) = vectory
              else
                ! velocity z-component
                field_display(ilocnum+ieoff) = vectorz
              endif
            endif

          enddo
        endif ! USE_HIGHRES_FOR_MOVIES
      enddo ! NSPEC_SURFACE_EXT_MESH

      ! copy coordinate arrays since the sorting routine does not preserve them
      xp_save(:) = xp(:)
      yp_save(:) = yp(:)
      zp_save(:) = zp(:)

      ! sort the list based upon coordinates to get rid of multiples
      print *,'sorting list of points'
      call get_global(npointot,xp,yp,zp,iglob,locval,ifseg,nglob, &
                      dble(minval(store_val_x(:))),dble(maxval(store_val_x(:))))

      ! print total number of points found
      print *
      print *,'found a total of ',nglob,' points'
      print *,'initial number of points (with multiples) was ',npointot


      !  normalize and scale vector field

      ! compute min and max of data value to normalize
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))

      if (plot_shaking_map) then
        ! print minimum and maximum amplitude in current snapshot
        print *
        print *,'minimum amplitude in current snapshot after removal = ',min_field_current
        print *,'maximum amplitude in current snapshot after removal = ',max_field_current
        print *
      else
        ! print minimum and maximum amplitude in current snapshot
        print *
        print *,'minimum amplitude in current snapshot = ',min_field_current
        print *,'maximum amplitude in current snapshot = ',max_field_current
        print *
      endif

      ! apply scaling in all cases for movies
      if (.not. plot_shaking_map) then

        ! normalizes values
        if (NORMALIZE_OUTPUT) then
          ! make sure range is always symmetric and center is in zero
          ! this assumption works only for fields that can be negative
          ! would not work for norm of vector for instance
          ! (we would lose half of the color palette if no negative values)
          max_absol = max(abs(min_field_current),abs(max_field_current))
          min_field_current = - max_absol
          max_field_current = + max_absol

          ! normalize field to [0:1]
          if (abs(max_field_current - min_field_current) > TINYVAL) &
            field_display(:) = (field_display(:) - min_field_current) / (max_field_current - min_field_current)

          ! rescale to [-1,1]
          field_display(:) = 2.*field_display(:) - 1.

          ! apply threshold to normalized field
          if (APPLY_THRESHOLD) &
            where(abs(field_display(:)) <= THRESHOLD) field_display = 0.
        endif

        ! apply non linear scaling to normalized field if needed
        if (NONLINEAR_SCALING) then
          where(field_display(:) >= 0.)
            field_display = field_display ** POWER_SCALING
          elsewhere
            field_display = - abs(field_display) ** POWER_SCALING
          endwhere
        endif

        ! normalizes values
        if (NORMALIZE_OUTPUT) then
          ! map back to [0,1]
          field_display(:) = (field_display(:) + 1.) / 2.

          ! map field to [0:255] for AVS color scale
          field_display(:) = 255. * field_display(:)
        endif

      ! apply scaling only if selected for shaking map
      else if (NONLINEAR_SCALING .and. iscaling_shake == 1) then

        ! normalize field to [0:1]
        if (abs(max_field_current) > TINYVAL) &
          field_display(:) = field_display(:) / max_field_current

        ! apply non linear scaling to normalized field
        field_display = field_display ** POWER_SCALING

        ! map field to [0:255] for AVS color scale
        field_display(:) = 255. * field_display(:)

      endif

      !--- ****** create AVS file using sorted list ******

      if (.not. plot_shaking_map) then
        if (inumber == 1) then
          ivalue = iframe
        else
          ivalue = it
        endif
      endif

      ! create file name and open file
      if (plot_shaking_map) then

        if (USE_OPENDX) then
          write(outputname,"('/DX_shaking_map.dx')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
        else if (USE_AVS) then
          write(outputname,"('/AVS_shaking_map.inp')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
       else if (USE_GMT) then
          write(outputname,"('/gmt_shaking_map.xyz')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
        else
          stop 'wrong output format selected'
        endif

      else

        if (USE_OPENDX) then
          write(outputname,"('/DX_movie_',i6.6,'.dx')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
        else if (USE_AVS) then
          write(outputname,"('/AVS_movie_',i6.6,'.inp')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
       else if (USE_GMT) then
          write(outputname,"('/gmt_movie_',i6.6,'.xyz')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
        else
          stop 'wrong output format selected'
        endif

      endif


      if (USE_GMT) then

        ! output list of points
        mask_point = .false.
        do ispec = 1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX
            ibool_number = iglob(ilocnum+ieoff)
            if (.not. mask_point(ibool_number)) then
              call utm_geo(long,lat,xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff),IUTM2LONGLAT)
              write(11,*) sngl(long),sngl(lat),sngl(field_display(ilocnum+ieoff))
            endif
            mask_point(ibool_number) = .true.
          enddo
        enddo

      else

        ! output list of points
        mask_point = .false.
        ipoin = 0
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX
            ibool_number = iglob(ilocnum+ieoff)
            if (.not. mask_point(ibool_number)) then
              ipoin = ipoin + 1
              ireorder(ibool_number) = ipoin
              if (USE_OPENDX) then
                write(11,*) sngl(xp_save(ilocnum+ieoff)),sngl(yp_save(ilocnum+ieoff)),sngl(zp_save(ilocnum+ieoff))
              else if (USE_AVS) then
                write(11,*) ireorder(ibool_number),sngl(xp_save(ilocnum+ieoff)), &
                            sngl(yp_save(ilocnum+ieoff)),sngl(zp_save(ilocnum+ieoff))
              endif
            endif
            mask_point(ibool_number) = .true.
          enddo
        enddo

        if (USE_OPENDX) &
          write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspectot_AVS_max,' data follows'

        ! output list of elements
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          ibool_number1 = iglob(ieoff + 1)
          ibool_number2 = iglob(ieoff + 2)
          ibool_number3 = iglob(ieoff + 3)
          ibool_number4 = iglob(ieoff + 4)
          if (USE_OPENDX) then
            ! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
            write(11,"(i10,1x,i10,1x,i10,1x,i10)") ireorder(ibool_number1)-1, &
              ireorder(ibool_number4)-1,ireorder(ibool_number2)-1,ireorder(ibool_number3)-1
          else
            ! AVS UCD format
            write(11,"(i10,' 1 quad ',i10,1x,i10,1x,i10,1x,i10)") ispec,ireorder(ibool_number1), &
              ireorder(ibool_number4),ireorder(ibool_number2),ireorder(ibool_number3)
          endif
        enddo

        if (USE_OPENDX) then
          write(11,*) 'attribute "element type" string "quads"'
          write(11,*) 'attribute "ref" string "positions"'
          write(11,*) 'object 3 class array type float rank 0 items ',nglob,' data follows'
        else
          ! AVS UCD format
          ! dummy text for labels
          write(11,*) '1 1'
          write(11,*) 'a, b'
        endif

        ! output data values
        mask_point = .false.
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX
            ibool_number = iglob(ilocnum+ieoff)
            if (.not. mask_point(ibool_number)) then
              if (USE_OPENDX) then
                if (plot_shaking_map) then
                  write(11,*) sngl(field_display(ilocnum+ieoff))
                else
                  write(11,*) sngl(field_display(ilocnum+ieoff))
                endif
              else
                ! AVS UCD format
                if (plot_shaking_map) then
                  write(11,*) ireorder(ibool_number),sngl(field_display(ilocnum+ieoff))
                else
                  write(11,*) ireorder(ibool_number),sngl(field_display(ilocnum+ieoff))
                endif
              endif
            endif
            mask_point(ibool_number) = .true.
          enddo
        enddo

        ! define OpenDX field
        if (USE_OPENDX) then
          write(11,*) 'attribute "dep" string "positions"'
          write(11,*) 'object "irregular positions irregular connections" class field'
          write(11,*) 'component "positions" value 1'
          write(11,*) 'component "connections" value 2'
          write(11,*) 'component "data" value 3'
          write(11,*) 'end'
        endif

      ! end of test for GMT format
      endif

      close(11)

    ! end of loop and test on all the time steps for all the movie images
   endif
enddo ! it

  print *
  print *,'done creating movie or shaking map'
  print *
  if (USE_OPENDX) print *, 'DX files are stored in ', trim(OUTPUT_FILES), '/DX_*.dx'
  if (USE_AVS) print *, 'AVS files are stored in ', trim(OUTPUT_FILES), '/AVS_*.inp'
  if (USE_GMT) print *, 'GMT files are stored in ', trim(OUTPUT_FILES), '/gmt_*.xyz'
  print *

  deallocate(store_val_x)
  deallocate(store_val_y)
  deallocate(store_val_z)
  deallocate(store_val_ux)
  deallocate(store_val_uy)
  deallocate(store_val_uz)

  ! deallocate arrays for sorting routine
  deallocate(iglob,locval)
  deallocate(ifseg)
  deallocate(xp,yp,zp)
  deallocate(xp_save,yp_save,zp_save)
  deallocate(field_display)
  deallocate(mask_point)
  deallocate(ireorder)

  if (USE_HIGHRES_FOR_MOVIES) then
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(display)
  endif

  end program create_movie_shakemap

