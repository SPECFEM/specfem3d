!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
!---  create a movie of vertical component of surface displacement or velocity
!---  in AVS or OpenDX format
!

program create_movie_GMT

  implicit none

  include 'constants.h'

  character(len=200) par_file, movie_data_prefix, start_frame, end_frame, &
       output_file_prefix
  integer ios1, ios2, nspectot_AVS_max
  ! threshold in percent of the maximum below which we cut the amplitude
  logical, parameter :: APPLY_THRESHOLD = .false.
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

  ! coefficient of power law used for non linear scaling
  logical, parameter :: NONLINEAR_SCALING = .false.
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.50_CUSTOM_REAL

  integer it,it1,it2,ivalue,ispec
  integer iformat,nframes,iframe,inumber,inorm,iscaling_shake
  integer ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4

  logical UNIQUE_FILE,plot_shaking_map

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,display
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord
  real(kind=CUSTOM_REAL) vectorx,vectory,vectorz

  double precision min_field_current,max_field_current,max_absol

  character(len=256) outputname

  integer iproc,ipoin

  ! GMT
  double precision lat,long,zscaling
  integer igmt

  ! for sorting routine
  integer npointot,ilocnum,nglob,i,j,ielm,ieoff,ispecloc
  integer, dimension(:), allocatable :: iglob,loc,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

  ! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
       store_val_x,store_val_y,store_val_z, &
       store_val_ux,store_val_uy,store_val_uz

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
            NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, UTM_MAX
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,DT,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH

  character(len=256) LOCAL_PATH,MODEL,CMTSOLUTION


  ! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
       NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
       NSPEC2D_BOTTOM,NSPEC2D_TOP, &
       NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

  double precision max_all_frames
  ! ************** PROGRAM STARTS HERE **************


  call getarg(1,movie_data_prefix)
  call getarg(2,par_file)
  call getarg(3,start_frame)
  call getarg(4,end_frame)
  call getarg(5,output_file_prefix)
  if (trim(movie_data_prefix) == '' .or. trim(par_file) == ''  &
       .or. trim(start_frame) == '' .or. trim(end_frame) == '' &
       .or. trim(output_file_prefix) == '') &
       stop 'Usage: xcreate_movie_GMT movie_data_prefix par_file start_frame end_frame output_file_prefix'

  read(start_frame, *,iostat=ios1) it1
  read(end_frame, *, iostat=ios2) it2
  if (ios1 /= 0 .or. ios2 /= 0) stop 'Error reading start_frame and end_frame'

! read the parameter file
  call read_parameter_file(par_file,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO, &
        SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,USE_REGULAR_MESH)

  print *
  print *,'There are ',NPROC,' slices numbered from 0 to ',NPROC-1
  print *

  print *, 'DT = ', DT , ' NSTEP = ', NSTEP
  print *

  if(SAVE_DISPLACEMENT) then
     print *,'Vertical displacement will be shown in movie'
  else
     print *,'Vertical velocity will be shown in movie'
  endif
  print *

  if(USE_HIGHRES_FOR_MOVIES) then
     print *, 'Movie is in high-resolution'
     ilocnum = NGLLSQUARE*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
  else
     print *, 'Movie is in low-resolution'
     ilocnum = NGNOD2D_AVS_DX*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
  endif
  allocate(store_val_x(ilocnum,0:NPROC-1))
  allocate(store_val_y(ilocnum,0:NPROC-1))
  allocate(store_val_z(ilocnum,0:NPROC-1))
  allocate(store_val_ux(ilocnum,0:NPROC-1))
  allocate(store_val_uy(ilocnum,0:NPROC-1))
  allocate(store_val_uz(ilocnum,0:NPROC-1))

  zscaling = 0.

  if(it1 == -1) then
     plot_shaking_map = .true.
     nframes = 1
     it1 = 1
     inorm = it2
     if(inorm < 1 .or. inorm > 3) stop 'incorrect value of inorm'
     it2 = 1
  else
     plot_shaking_map = .false.
     print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
     print *
     print *
     print *,'looping from ',it1,' to ',it2,' frame'
     ! count number of movie frames
     nframes = 0
     do it = it1,it2
        nframes = nframes + 1
     enddo
     print *
     print *,'total number of frames will be ',nframes
     if(nframes == 0) stop 'null number of frames'
     max_all_frames = -100.0
  endif


  ! define the total number of elements at the surface
  if(USE_HIGHRES_FOR_MOVIES) then
     nspectot_AVS_max = NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)
  else
     nspectot_AVS_max = NEX_XI * NEX_ETA
  endif

  print *
  print *,'there are a total of ',nspectot_AVS_max,' elements at the surface'
  print *

  ! maximum theoretical number of points at the surface
  npointot = NGNOD2D_AVS_DX * nspectot_AVS_max

  ! allocate arrays for sorting routine
  allocate(iglob(npointot),loc(npointot))
  allocate(ifseg(npointot))
  allocate(xp(npointot),yp(npointot),zp(npointot))
  allocate(xp_save(npointot),yp_save(npointot),zp_save(npointot))
  allocate(field_display(npointot))
  allocate(mask_point(npointot))
  allocate(ireorder(npointot))

  ! --------------------------------------

  if(USE_HIGHRES_FOR_MOVIES) then
     allocate(x(NGLLX,NGLLY))
     allocate(y(NGLLX,NGLLY))
     allocate(z(NGLLX,NGLLY))
     allocate(display(NGLLX,NGLLY))
  endif

  iframe = 0

  ! loop on all the time steps in the range entered
  do it = it1,it2

     ! check if time step corresponds to a movie frame
     iframe = iframe + 1
     ivalue = it * NTSTEP_BETWEEN_FRAMES
!     print *, 'ivalue = ' ,ivalue
     if(plot_shaking_map) then
        print *,'reading shaking map snapshot'
     else
        print *,'reading snapshot frame ',it,' out of ',NSTEP/NTSTEP_BETWEEN_FRAMES
     endif
     print *

     ! read all the elements from the same file
     if(plot_shaking_map) then
        write(outputname,"('/shakingdata')")
     else
        write(outputname,"('/moviedata',i6.6)") ivalue
     endif
     outputname = trim(movie_data_prefix) // trim(outputname)
     open(unit=IOUT,file=outputname,status='old',form='unformatted',iostat=ios1)
     if (ios1 /= 0) then
        print *, 'Error opening file ', trim(outputname)
        stop
     else
       print *, 'reading file ', trim(outputname)
     endif
     read(IOUT) store_val_x
     read(IOUT) store_val_y
     read(IOUT) store_val_z
     read(IOUT) store_val_ux
     read(IOUT) store_val_uy
     read(IOUT) store_val_uz
     close(IOUT)

     call utm_geo(long,lat,dble(minval(store_val_x)),dble(minval(store_val_y)), &
                   UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
!     print *
!     print *, 'min long = ', long, '  min y = ', lat
     call utm_geo(long,lat,dble(maxval(store_val_x)),dble(maxval(store_val_y)), &
                   UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
!     print *
!     print *, 'max long = ', long, '  max lat = ', lat

     ! clear number of elements kept
     ispec = 0

     ! read points for all the slices
     do iproc = 0,NPROC-1

        ! reset point number
        ipoin = 0

        do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA

           if(USE_HIGHRES_FOR_MOVIES) then
              ! assign the OpenDX "elements"

              do j = 1,NGLLY
                 do i = 1,NGLLX

                    ipoin = ipoin + 1

                    xcoord = store_val_x(ipoin,iproc)
                    ycoord = store_val_y(ipoin,iproc)
                    zcoord = store_val_z(ipoin,iproc)

                    ! amplify topography, otherwise too flat to see anything
                    zcoord = zcoord * zscaling

                    vectorx = store_val_ux(ipoin,iproc)
                    vectory = store_val_uy(ipoin,iproc)
                    vectorz = store_val_uz(ipoin,iproc)

                    x(i,j) = xcoord
                    y(i,j) = ycoord
                    z(i,j) = zcoord

                    if(plot_shaking_map) then
                       if(inorm == 1) then
                          display(i,j) = vectorx
                       else if(inorm == 2) then
                          display(i,j) = vectory
                       else
                          display(i,j) = vectorz
                       endif
                    else
                       display(i,j) = vectorz
                    endif

                 enddo
              enddo

              ! assign the values of the corners of the OpenDX "elements"
              ispec = ispec + 1
              ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)

              do j = 1,NGLLY-1
                 do i = 1,NGLLX-1
                    ieoff = NGNOD2D_AVS_DX*(ielm+(i-1)+(j-1)*(NGLLX-1))
                    do ilocnum = 1,NGNOD2D_AVS_DX

                       if(ilocnum == 1) then
                          xp(ieoff+ilocnum) = dble(x(i,j))
                          yp(ieoff+ilocnum) = dble(y(i,j))
                          zp(ieoff+ilocnum) = dble(z(i,j))
                          field_display(ieoff+ilocnum) = dble(display(i,j))
                       elseif(ilocnum == 2) then
                          xp(ieoff+ilocnum) = dble(x(i+1,j))
                          yp(ieoff+ilocnum) = dble(y(i+1,j))
                          zp(ieoff+ilocnum) = dble(z(i+1,j))
                          field_display(ieoff+ilocnum) = dble(display(i+1,j))
                       elseif(ilocnum == 3) then
                          xp(ieoff+ilocnum) = dble(x(i+1,j+1))
                          yp(ieoff+ilocnum) = dble(y(i+1,j+1))
                          zp(ieoff+ilocnum) = dble(z(i+1,j+1))
                          field_display(ieoff+ilocnum) = dble(display(i+1,j+1))
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

              ispec = ispec + 1
              ieoff = NGNOD2D_AVS_DX*(ispec-1)

              ! four points for each element
              do ilocnum = 1,NGNOD2D_AVS_DX

                 ipoin = ipoin + 1

                 xcoord = store_val_x(ipoin,iproc)
                 ycoord = store_val_y(ipoin,iproc)
                 zcoord = store_val_z(ipoin,iproc)

                 ! amplify topography, otherwise too flat to see anything
                 zcoord = zcoord * zscaling

                 vectorx = store_val_ux(ipoin,iproc)
                 vectory = store_val_uy(ipoin,iproc)
                 vectorz = store_val_uz(ipoin,iproc)

                 xp(ilocnum+ieoff) = dble(xcoord)
                 yp(ilocnum+ieoff) = dble(ycoord)
                 zp(ilocnum+ieoff) = dble(zcoord)

                 ! show vertical component of displacement or velocity in the movie
                 ! or show norm of vector if shaking map
                 ! for shaking map, norm of U stored in ux, V in uy and A in uz
                 if(plot_shaking_map) then
                    if(inorm == 1) then
                       field_display(ilocnum+ieoff) = dble(vectorx)
                    else if(inorm == 2) then
                       field_display(ilocnum+ieoff) = dble(vectory)
                    else
                       field_display(ilocnum+ieoff) = dble(vectorz)
                    endif
                 else
                    field_display(ilocnum+ieoff) = dble(vectorz)  ! only plotting z component for now.
                    !field_display(ilocnum+ieoff) = dble(vectorx * 0.87 - vectory * 0.5)  ! plot hroizontal movie
                 endif

              enddo

           endif

        enddo
     enddo

     ! copy coordinate arrays since the sorting routine does not preserve them
     xp_save(:) = xp(:)
     yp_save(:) = yp(:)
     zp_save(:) = zp(:)

     !--- sort the list based upon coordinates to get rid of multiples
     print *,'sorting list of points'
     call get_global_AVS(nspectot_AVS_max,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,UTM_X_MIN,UTM_X_MAX)

     !--- print total number of points found
     print *
     print *,'found a total of ',nglob,' points'
     print *,'initial number of points (with multiples) was ',npointot


     !--- normalize and scale vector field

     ! compute min and max of data value to normalize
     min_field_current = minval(field_display(:))
     max_field_current = maxval(field_display(:))


     if (max_field_current > max_all_frames) max_all_frames = max_field_current
     ! print minimum and maximum amplitude in current snapshot
     print *
     print *,'minimum amplitude in current snapshot = ',min_field_current
     print *,'maximum amplitude in current snapshot = ',max_field_current
     print *

     if(plot_shaking_map) then

        ! normalize field to [0:1]
        field_display(:) = field_display(:) / max_field_current

        ! apply non linear scaling to normalized field
        field_display = field_display ** POWER_SCALING

        ! map field to [0:255] for AVS color scale
        field_display(:) = 255. * field_display(:)

     endif

     !--- ****** create AVS file using sorted list ******


     ! create file name and open file
     if(plot_shaking_map) then
        write(outputname,"('/gmt_shaking_map.xyz')")
        open(unit=11,file=trim(output_file_prefix)//outputname,status='unknown',iostat=ios1)
     else
        write(outputname,"('/gmt_movie_',i6.6,'.xyz')") it
        outputname = trim(output_file_prefix) //trim(outputname)
        call system('\\rm -f ' // trim(outputname))
        open(unit=11,file=outputname,status='unknown')!,form='unformatted',access='direct',recl = 3 * 8,iostat=ios1)
     endif
     if (ios1 /= 0) stop 'Error opening output file'

     igmt = 1
     ! output list of points
     mask_point = .false.
     do ispec=1,nspectot_AVS_max
        ieoff = NGNOD2D_AVS_DX*(ispec-1)
        ! four points for each element
        do ilocnum = 1,NGNOD2D_AVS_DX
           ibool_number = iglob(ilocnum+ieoff)
           if(.not. mask_point(ibool_number)) then
              call utm_geo(long,lat,xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff), &
                   UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)
              if (plot_shaking_map) then
                 write(11,*) long,lat,field_display(ilocnum+ieoff)
              else
                 !write(11,rec=igmt) long,lat,field_display(ilocnum+ieoff)
                  write (11,'(3e17.6)') long, lat, field_display(ilocnum+ieoff)
                 igmt = igmt + 1
              endif
           endif
           mask_point(ibool_number) = .true.
        enddo
     enddo
     print *, 'total number of record is ', igmt - 1
     print *, ' '
     print *, ' '

  enddo

  ! step number for AVS multistep file
401 format('step',i1,' image',i1)
402 format('step',i2,' image',i2)
403 format('step',i3,' image',i3)
404 format('step',i4,' image',i4)

  print *,'GMT files are stored in OUTPUT_FILES/gmt_*.xyz'
  print *
  print *, 'The maximum absolute value of all frames is ', max_all_frames
  print *
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




end program create_movie_GMT

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

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

  double precision UTM_X_MIN,UTM_X_MAX

! define geometrical tolerance based upon typical size of the model
  SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)

! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGNOD2D_AVS_DX*(ispec-1)
    do ilocnum=1,NGNOD2D_AVS_DX
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

