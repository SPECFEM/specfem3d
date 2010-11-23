program convert_movie_real_to_double

  implicit none

  include 'constants.h'

  character(len=100) par_file, movie_data_prefix, start_frame, end_frame, &
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

  character(len=256) outputname,outputname1

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
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: store_val
  real*8, dimension(:,:), allocatable :: store_val_double

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_ETA,NEX_XI, &
             NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE
  integer NSOURCES

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS
  logical ANISOTROPY,SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT,USE_HIGHRES_FOR_MOVIES
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=256) LOCAL_PATH,clean_LOCAL_PATH,final_LOCAL_PATH,prname
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
  if (trim(movie_data_prefix) == '' .or. trim(par_file) == ''  &
       .or. trim(start_frame) == '' .or. trim(end_frame) == '' ) &
       stop 'Usage: xcreate_movie_GMT movie_data_prefix par_file start_frame end_frame output_file_prefix'

  read(start_frame, *,iostat=ios1) it1
  read(end_frame, *, iostat=ios2) it2
  if (ios1 /= 0 .or. ios2 /= 0) stop 'Error reading start_frame and end_frame'

! read the parameter file
  call read_parameter_file(par_file,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES, &
        SAVE_AVS_DX_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION,NTSTEP_BETWEEN_OUTPUT_INFO)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB)

  if(USE_HIGHRES_FOR_MOVIES) then
    ilocnum = NGLLSQUARE*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
  else
    ilocnum = NGNOD2D_AVS_DX*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
  endif
  allocate(store_val(ilocnum,0:NPROC-1))
  allocate(store_val_double(ilocnum,0:NPROC-1))
  print *, 'data = ', ilocnum*NPROC*8/1024/1024,'  MB'

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


 ! loop on all the time steps in the range entered
  do it = it1,it2

     ! check if time step corresponds to a movie frame
     iframe = iframe + 1
     ivalue = it * NTSTEP_BETWEEN_FRAMES
     print *, 'ivalue = ' ,ivalue
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
     outputname1 =  trim(outputname)//'.new'
     open(unit=IIN,file=outputname,status='old',form='unformatted',iostat=ios1)
     open(unit=IOUT,file=outputname1,form='unformatted',iostat=ios2)

     if (ios1 /= 0 .or. ios2 /= 0) then
        print *, 'Error opening file ', trim(outputname), ' and ', trim(outputname1)
        stop
     else
       print *, 'reading file ', trim(outputname)
       print *, 'writing file ', trim(outputname1)
     endif

     do i = 1, 6
       read(IIN,iostat=ios1) store_val_double
       if (ios1 /= 0) stop 'Error reading store_val_double file'
       store_val(:,:) = sngl(store_val_double(:,:))
       write(IOUT,iostat=ios1) store_val
       if (ios1 /= 0) stop 'Error writing store_val file'
     enddo

   close(IIN)
   close(IOUT)

 end do



 end program convert_movie_real_to_double
