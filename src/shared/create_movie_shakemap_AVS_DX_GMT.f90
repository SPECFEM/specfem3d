!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

!
!---  create a movie of the vertical component of surface displacement or velocity
!---  or a ShakeMap(R) (i.e. map of the maximum absolute value of the two horizontal components
!---  of the velocity vector) in AVS, OpenDX or GMT format
!

  program create_movie_shakemap

  implicit none

  include "constants.h"
  include "../../OUTPUT_FILES/surface_from_mesher.h"

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

!-------------------------------------------------------------------------------------------------

  integer it,it1,it2,ivalue,nspectot_AVS_max,ispec
  integer iformat,nframes,iframe,inumber,inorm,iscaling_shake
  integer ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4

  logical USE_OPENDX,USE_AVS,USE_GMT,plot_shaking_map

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,display
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) vectorx,vectory,vectorz

  double precision min_field_current,max_field_current,max_absol

  character(len=256) outputname

  integer ipoin

  ! GMT
  double precision lat,long

  ! for sorting routine
  integer npointot,ilocnum,nglob,i,j,ielm,ieoff,ispecloc
  integer, dimension(:), allocatable :: iglob,loc,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

  ! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

  ! parameters read from parameter file
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,NGNOD,NGNOD2D
  double precision DT
  double precision HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML
  logical ATTENUATION,USE_OLSEN_ATTENUATION, &
          APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE
  logical STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE
  logical ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  logical PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID
  character(len=256) OUTPUT_FILES,LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH
  integer NPROC
  integer ier
  integer MOVIE_TYPE,IMODEL

!--------------------------------------------
!!!! NL NL for external meshes
!--------------------------------------------
  ! muting source region
  logical, parameter :: MUTE_SOURCE = .false.
  real(kind=CUSTOM_REAL), parameter :: RADIUS_TO_MUTE = 1000._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: X_SOURCE_EXT_MESH = -9023.021484375
  real(kind=CUSTOM_REAL), parameter :: Y_SOURCE_EXT_MESH = 6123.611328125
  real(kind=CUSTOM_REAL), parameter :: Z_SOURCE_EXT_MESH = 17.96331405639648
!--------------------------------------------
!!!! NL NL

  ! order of points representing the 2D square element
  integer,dimension(NGNOD2D_FOUR_CORNERS_AVS_DX),parameter :: iorder = (/1,3,2,4/)


! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *

  print *
  print *,'reading parameter file'
  print *

  ! read the parameter file
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
        APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
        USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
        USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
        PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))

  ! only one global array for movie data, but stored for all surfaces defined
  ! in file 'surface_from_mesher.h'
  if(USE_HIGHRES_FOR_MOVIES) then
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGLLSQUARE
  else
     ilocnum = NSPEC_SURFACE_EXT_MESH*NGNOD2D_FOUR_CORNERS_AVS_DX
  endif
  print*,'  moviedata element surfaces: ',NSPEC_SURFACE_EXT_MESH
  print*,'  moviedata total elements all: ',ilocnum
  print *

  if(SAVE_DISPLACEMENT) then
    print *,'Vertical displacement will be shown in movie'
  else
    print *,'Vertical velocity will be shown in movie'
  endif
  print *


  ! user input
  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format'
  print *,'3 = create files in GMT xyz Ascii long/lat/Uz format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if(iformat < 1 .or. iformat > 3) stop 'exiting...'

  plot_shaking_map = .false.
  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *
  print *,'enter first time step of movie (e.g. 1, enter -1 for shaking map)'
  read(5,*) it1
  if(it1 == 0 ) it1 = 1
  if(it1 == -1) plot_shaking_map = .true.
  if(.not. plot_shaking_map) then
    print *,'enter last time step of movie (e.g. ',NSTEP,')'
    read(5,*) it2
    print *
    print *,'1 = define file names using frame number'
    print *,'2 = define file names using time step number'
    print *,'any other value = exit'
    print *
    print *,'enter value:'
    read(5,*) inumber
    if(inumber<1 .or. inumber>2) stop 'exiting...'
    print *
    print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'
    ! count number of movie frames
    nframes = 0
    do it = it1,it2
      if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
    enddo
  else
    ! only one frame if shaking map
    nframes = 1
    it1 = 1
    it2 = 1
  endif
  print *
  print *,'total number of frames will be ',nframes
  if(nframes == 0) stop 'null number of frames'

  iscaling_shake = 0
  if(plot_shaking_map) then
    print *
    print *,'norm to display in shaking map:'
    print *,'1=displacement  2=velocity  3=acceleration'
    print *
    read(5,*) inorm
    if(inorm < 1 .or. inorm > 3) stop 'incorrect value of inorm'
    print *
    print *,'apply non-linear scaling to shaking map:'
    print *,'1=non-linear  2=no scaling'
    print *
    read(5,*) iscaling_shake
    if(iscaling_shake < 1 .or. iscaling_shake > 2) stop 'incorrect value of iscaling_shake'
  else
    print *
    print *,'movie data:'
    print *,'1= norm of velocity  2=velocity x-comp 3=velocity y-comp 4=velocity z-comp'
    print *
    read(5,*) inorm
    if(inorm < 1 .or. inorm > 4) stop 'incorrect value of inorm'
  endif

! file format flags
  if(iformat == 1) then
    USE_OPENDX = .true.
    USE_AVS = .false.
    USE_GMT = .false.
  else if(iformat == 2) then
    USE_OPENDX = .false.
    USE_AVS = .true.
    USE_GMT = .false.
  else
    USE_OPENDX = .false.
    USE_AVS = .false.
    USE_GMT = .true.
  endif

  ! define the total number of elements at the surface
  if(USE_HIGHRES_FOR_MOVIES) then
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH * (NGLLX-1) * (NGLLY-1)
  else
     nspectot_AVS_max = NSPEC_SURFACE_EXT_MESH
  endif

  ! maximum theoretical number of points at the surface
  npointot = NGNOD2D_FOUR_CORNERS_AVS_DX * nspectot_AVS_max

  ! allocate arrays for sorting routine
  allocate(iglob(npointot),loc(npointot), &
          ifseg(npointot), &
          xp(npointot),yp(npointot),zp(npointot), &
          xp_save(npointot),yp_save(npointot),zp_save(npointot), &
          field_display(npointot), &
          mask_point(npointot), &
          ireorder(npointot),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for sorting routine'

  ! allocates data arrays
  allocate(store_val_x(ilocnum), &
          store_val_y(ilocnum), &
          store_val_z(ilocnum), &
          store_val_ux(ilocnum), &
          store_val_uy(ilocnum), &
          store_val_uz(ilocnum),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for data arrays'

  if(USE_HIGHRES_FOR_MOVIES) then
    allocate(x(NGLLX,NGLLY), &
            y(NGLLX,NGLLY), &
            z(NGLLX,NGLLY), &
            display(NGLLX,NGLLY),stat=ier)
    if( ier /= 0 ) stop 'error allocating arrays for highres'
  endif

  ! user output
  print *
  print *,'there are a total of ',nspectot_AVS_max,' elements at the surface'
  print *
  print *
  if(APPLY_THRESHOLD .and. .not. plot_shaking_map) &
    print *,'Will apply a threshold to amplitude below ',100.*THRESHOLD,' %'
  if(NONLINEAR_SCALING .and. (.not. plot_shaking_map .or. iscaling_shake == 1)) &
    print *,'Will apply a non linear scaling with coef ',POWER_SCALING


  iframe = 0

! loop on all the time steps in the range entered
  do it = it1,it2

    ! check if time step corresponds to a movie frame
    if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0 .or. plot_shaking_map) then

      iframe = iframe + 1

      print *
      if(plot_shaking_map) then
        print *,'reading shaking map snapshot'
      else
        print *,'reading snapshot time step ',it,' out of ',NSTEP
      endif
      print *

      ! read all the elements from the same file
      if(plot_shaking_map) then
        write(outputname,"('/shakingdata')")
      else
        write(outputname,"('/moviedata',i6.6)") it
      endif
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname),status='old', &
            action='read',form='unformatted',iostat=ier)
      if( ier /= 0 ) then
        print*,'error: ',trim(OUTPUT_FILES)//trim(outputname)
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

        if(USE_HIGHRES_FOR_MOVIES) then
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
              if(plot_shaking_map) then
                !!!! NL NL mute value near source
                if ( (sqrt(((x(i,j) - (X_SOURCE_EXT_MESH))**2 + &
                     (y(i,j) - (Y_SOURCE_EXT_MESH))**2 + &
                     (z(i,j) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) &
                     .and. MUTE_SOURCE) then

                  display(i,j) = 0.
                else
                  ! chooses norm
                  if(inorm == 1) then
                    ! norm displacement
                    display(i,j) = vectorx
                  else if(inorm == 2) then
                    ! norm velocity
                    display(i,j) = vectory
                  else
                    ! norm acceleration
                    display(i,j) = vectorz
                  endif
                endif
              else
                ! movie
                if(inorm == 1) then
                  ! norm of velocity
                  display(i,j) = sqrt(vectorz**2+vectory**2+vectorx**2)
                else if( inorm == 2 ) then
                  ! velocity x-component
                  display(i,j) = vectorx
                else if( inorm == 3 ) then
                  ! velocity y-component
                  display(i,j) = vectory
                else if( inorm == 4 ) then
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

                if(ilocnum == 1) then
                  xp(ieoff+ilocnum) = dble(x(i,j))
                  yp(ieoff+ilocnum) = dble(y(i,j))
                  zp(ieoff+ilocnum) = dble(z(i,j))
                  field_display(ieoff+ilocnum) = dble(display(i,j))
                else if(ilocnum == 2) then

                  ! accounts for different ordering of square points
                  xp(ieoff+ilocnum) = dble(x(i+1,j+1))
                  yp(ieoff+ilocnum) = dble(y(i+1,j+1))
                  zp(ieoff+ilocnum) = dble(z(i+1,j+1))
                  field_display(ieoff+ilocnum) = dble(display(i+1,j+1))

                else if(ilocnum == 3) then

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
            if(plot_shaking_map) then
              !!!! NL NL mute value near source
              if ( (sqrt(((dble(xcoord) - (X_SOURCE_EXT_MESH))**2 + &
                     (dble(ycoord) - (Y_SOURCE_EXT_MESH))**2 + &
                     (dble(zcoord) - (Z_SOURCE_EXT_MESH))**2)) < RADIUS_TO_MUTE) &
                     .and. MUTE_SOURCE) then
                  field_display(ilocnum+ieoff) = 0.
              else
                if(inorm == 1) then
                  ! norm of displacement
                  field_display(ilocnum+ieoff) = dble(vectorx)
                else if(inorm == 2) then
                  ! norm of velocity
                  field_display(ilocnum+ieoff) = dble(vectory)
                else
                  ! norm of acceleration
                  field_display(ilocnum+ieoff) = dble(vectorz)
                endif
              endif
            else
              ! movie
              if(inorm == 1) then
                ! norm of velocity
                field_display(ilocnum+ieoff) = sqrt(vectorz**2+vectory**2+vectorx**2)
              else if( inorm == 2 ) then
                ! velocity x-component
                field_display(ilocnum+ieoff) = vectorx
              else if( inorm == 3 ) then
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
      call get_global_AVS(nspectot_AVS_max,xp,yp,zp,iglob,loc,ifseg,nglob,npointot, &
           dble(minval(store_val_x(:))),dble(maxval(store_val_x(:))))

      ! print total number of points found
      print *
      print *,'found a total of ',nglob,' points'
      print *,'initial number of points (with multiples) was ',npointot


      !  normalize and scale vector field

      ! compute min and max of data value to normalize
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))

      if(plot_shaking_map) then
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
      if(.not. plot_shaking_map) then

        ! normalizes values
        if( NORMALIZE_OUTPUT ) then
          ! make sure range is always symmetric and center is in zero
          ! this assumption works only for fields that can be negative
          ! would not work for norm of vector for instance
          ! (we would lose half of the color palette if no negative values)
          max_absol = max(abs(min_field_current),abs(max_field_current))
          min_field_current = - max_absol
          max_field_current = + max_absol

          ! normalize field to [0:1]
          if( abs(max_field_current - min_field_current) > TINYVAL ) &
            field_display(:) = (field_display(:) - min_field_current) / (max_field_current - min_field_current)

          ! rescale to [-1,1]
          field_display(:) = 2.*field_display(:) - 1.

          ! apply threshold to normalized field
          if(APPLY_THRESHOLD) &
            where(abs(field_display(:)) <= THRESHOLD) field_display = 0.
        endif

        ! apply non linear scaling to normalized field if needed
        if(NONLINEAR_SCALING) then
          where(field_display(:) >= 0.)
            field_display = field_display ** POWER_SCALING
          elsewhere
            field_display = - abs(field_display) ** POWER_SCALING
          endwhere
        endif

        ! normalizes values
        if( NORMALIZE_OUTPUT ) then
          ! map back to [0,1]
          field_display(:) = (field_display(:) + 1.) / 2.

          ! map field to [0:255] for AVS color scale
          field_display(:) = 255. * field_display(:)
        endif

      ! apply scaling only if selected for shaking map
      else if(NONLINEAR_SCALING .and. iscaling_shake == 1) then

        ! normalize field to [0:1]
        if( abs(max_field_current) > TINYVAL ) &
          field_display(:) = field_display(:) / max_field_current

        ! apply non linear scaling to normalized field
        field_display = field_display ** POWER_SCALING

        ! map field to [0:255] for AVS color scale
        field_display(:) = 255. * field_display(:)

      endif

      !--- ****** create AVS file using sorted list ******

      if(.not. plot_shaking_map) then
        if(inumber == 1) then
          ivalue = iframe
        else
          ivalue = it
        endif
      endif

      ! create file name and open file
      if(plot_shaking_map) then

        if(USE_OPENDX) then
          write(outputname,"('/DX_shaking_map.dx')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
        else if(USE_AVS) then
          write(outputname,"('/AVS_shaking_map.inp')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
       else if(USE_GMT) then
          write(outputname,"('/gmt_shaking_map.xyz')")
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
        else
          stop 'wrong output format selected'
        endif

      else

        if(USE_OPENDX) then
          write(outputname,"('/DX_movie_',i6.6,'.dx')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
        else if(USE_AVS) then
          write(outputname,"('/AVS_movie_',i6.6,'.inp')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
          write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
       else if(USE_GMT) then
          write(outputname,"('/gmt_movie_',i6.6,'.xyz')") ivalue
          open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
        else
          stop 'wrong output format selected'
        endif

      endif


      if(USE_GMT) then

        ! output list of points
        mask_point = .false.
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          do ilocnum = 1,NGNOD2D_FOUR_CORNERS_AVS_DX
            ibool_number = iglob(ilocnum+ieoff)
            if(.not. mask_point(ibool_number)) then
              call utm_geo(long,lat,xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff), &
                      UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
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
            if(.not. mask_point(ibool_number)) then
              ipoin = ipoin + 1
              ireorder(ibool_number) = ipoin
              if(USE_OPENDX) then
                write(11,*) sngl(xp_save(ilocnum+ieoff)),sngl(yp_save(ilocnum+ieoff)),sngl(zp_save(ilocnum+ieoff))
              else if(USE_AVS) then
                write(11,*) ireorder(ibool_number),sngl(xp_save(ilocnum+ieoff)), &
                            sngl(yp_save(ilocnum+ieoff)),sngl(zp_save(ilocnum+ieoff))
              endif
            endif
            mask_point(ibool_number) = .true.
          enddo
        enddo

        if(USE_OPENDX) &
          write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspectot_AVS_max,' data follows'

        ! output list of elements
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
          ! four points for each element
          ibool_number1 = iglob(ieoff + 1)
          ibool_number2 = iglob(ieoff + 2)
          ibool_number3 = iglob(ieoff + 3)
          ibool_number4 = iglob(ieoff + 4)
          if(USE_OPENDX) then
            ! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
            write(11,"(i10,1x,i10,1x,i10,1x,i10)") ireorder(ibool_number1)-1, &
              ireorder(ibool_number4)-1,ireorder(ibool_number2)-1,ireorder(ibool_number3)-1
          else
            write(11,"(i10,' 1 quad ',i10,1x,i10,1x,i10,1x,i10)") ispec,ireorder(ibool_number1), &
              ireorder(ibool_number4),ireorder(ibool_number2),ireorder(ibool_number3)
          endif
        enddo

        if(USE_OPENDX) then
          write(11,*) 'attribute "element type" string "quads"'
          write(11,*) 'attribute "ref" string "positions"'
          write(11,*) 'object 3 class array type float rank 0 items ',nglob,' data follows'
        else
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
            if(.not. mask_point(ibool_number)) then
              if(USE_OPENDX) then
                if(plot_shaking_map) then
                  write(11,*) sngl(field_display(ilocnum+ieoff))
                else
                  write(11,*) sngl(field_display(ilocnum+ieoff))
                endif
              else
                if(plot_shaking_map) then
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
        if(USE_OPENDX) then
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
  if(USE_OPENDX) print *,'DX files are stored in ', trim(OUTPUT_FILES), '/DX_*.dx'
  if(USE_AVS) print *,'AVS files are stored in ', trim(OUTPUT_FILES), '/AVS_*.inp'
  if(USE_GMT) print *,'GMT files are stored in ', trim(OUTPUT_FILES), '/gmt_*.xyz'
  print *

  deallocate(store_val_x)
  deallocate(store_val_y)
  deallocate(store_val_z)
  deallocate(store_val_ux)
  deallocate(store_val_uy)
  deallocate(store_val_uz)

  ! deallocate arrays for sorting routine
  deallocate(iglob,loc)
  deallocate(ifseg)
  deallocate(xp,yp,zp)
  deallocate(xp_save,yp_save,zp_save)
  deallocate(field_display)
  deallocate(mask_point)
  deallocate(ireorder)

  if(USE_HIGHRES_FOR_MOVIES) then
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(display)
  endif

  end program create_movie_shakemap

!
!=====================================================================
!

  subroutine get_global_AVS(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  include "constants.h"

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision SMALLVALTOL

  integer npointot
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  integer nspec,nglob

  integer ispec,i,j,ier
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

  double precision UTM_X_MIN,UTM_X_MAX

! define geometrical tolerance based upon typical size of the model
    SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)
    print *, 'UTM_X_MAX', UTM_X_MAX
    print *, 'UTM_X_MIN', UTM_X_MIN
    print *, 'SMALLVALTOL', SMALLVALTOL

! dynamically allocate arrays
  allocate(ind(npointot), &
          ninseg(npointot), &
          iwork(npointot), &
          work(npointot),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays ind etc.'

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGNOD2D_FOUR_CORNERS_AVS_DX*(ispec-1)
    do ilocnum=1,NGNOD2D_FOUR_CORNERS_AVS_DX
      loc(ilocnum+ieoff)=ilocnum+ieoff
    enddo
  enddo

  ifseg(:)=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npointot

  do j=1,NDIM

! sort within each segment
  ioff=1
  do iseg=1,nseg
    if(j == 1) then
      call rank(xp(ioff),ind,ninseg(iseg))
    else if(j == 2) then
      call rank(yp(ioff),ind,ninseg(iseg))
    else
      call rank(zp(ioff),ind,ninseg(iseg))
    endif
    call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
    ioff=ioff+ninseg(iseg)
  enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
  if(j == 1) then
    do i=2,npointot
      if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else if(j == 2) then
    do i=2,npointot
      if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else
    do i=2,npointot
      if(dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  endif

! count up number of different segments
  nseg=0
  do i=1,npointot
    if(ifseg(i)) then
      nseg=nseg+1
      ninseg(nseg)=1
    else
      ninseg(nseg)=ninseg(nseg)+1
    endif
  enddo
  enddo

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot
    if(ifseg(i)) ig=ig+1
    iglob(loc(i))=ig
  enddo

  nglob=ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global_AVS

! -----------------------------------

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q<A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   goto 200
   ENDIF
   IND(I)=INDX
  goto 100
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all

