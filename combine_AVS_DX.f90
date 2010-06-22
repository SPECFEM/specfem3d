!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

! combine AVS or DX global data files to check the mesh
! this is done in postprocessing after running the mesh generator
! can combine full mesh, edges only or faces only

  program combine_AVS_DX

  implicit none

  include "constants.h"

  integer iproc,nspec,npoin
  integer ispec
  integer iglob1,iglob2,iglob3,iglob4
  integer ipoin,numpoin,iglobpointoffset,ntotpoin,ntotspec
  integer numelem,iglobelemoffset,idoubling,maxdoubling,iformat,ivalue,icolor
  integer imaterial,imatprop
  integer nrec,ir
  integer ntotpoinAVS_DX,ntotspecAVS_DX

  real random_val
  integer ival_color
  integer, dimension(:), allocatable :: random_colors

  double precision xval,yval,zval
  double precision val_color

! for source location
  integer yr,jda,ho,mi
  double precision x_target_source,y_target_source,z_target_source
  double precision x_source_quad1,y_source_quad1,z_source_quad1
  double precision x_source_quad2,y_source_quad2,z_source_quad2
  double precision x_source_quad3,y_source_quad3,z_source_quad3
  double precision x_source_quad4,y_source_quad4,z_source_quad4
  double precision sec

  double precision, dimension(:), allocatable :: hdur,t_cmt,lat,long,depth
  double precision, dimension(:,:), allocatable :: moment_tensor

  logical USE_OPENDX

! for receiver location
  integer irec,ios
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name
  character(len=256) dummystring

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target

! processor identification
  character(len=256) prname

! small offset for source and receiver line in AVS_DX
! (small compared to normalized radius of the Earth)

! offset to represent source and receivers for model
  double precision, parameter :: small_offset = 2000.d0

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_XI,NEX_ETA, &
             NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  double precision zscaling

  character(len=256) OUTPUT_FILES,LOCAL_PATH,MODEL,filtered_rec_filename

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

  integer proc_p1,proc_p2

! ************** PROGRAM STARTS HERE **************

! only for old regular meshes!

  print *
  print *,'Recombining all AVS or DX files for slices'
  print *

  print *
  print *,'reading parameter file'
  print *

! read the parameter file
  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

  if(.not. SAVE_MESH_FILES) stop 'AVS or DX files were not saved by the mesher'

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  allocate(hdur(NSOURCES))
  allocate(t_cmt(NSOURCES))
  allocate(lat(NSOURCES))
  allocate(long(NSOURCES))
  allocate(depth(NSOURCES))
  allocate(moment_tensor(6,NSOURCES))

  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if(iformat<1 .or. iformat>2) stop 'exiting...'
  if(iformat == 1) then
    USE_OPENDX = .true.
  else
    USE_OPENDX = .false.
  endif

  print *
  print *,'1 = edges of all the slices only'
  print *,'2 = surface of the model only'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) ivalue
  if(ivalue<1 .or. ivalue>2) stop 'exiting...'

! apply scaling to topography if needed
  if(ivalue == 2) then
    print *
    print *,'scaling to apply to Z to amplify topography (1. to do nothing, 0. to get flat surface):'
    read(5,*) zscaling
  else
    zscaling = 1.d0
  endif

  print *
  print *,'1 = color by doubling flag'
  print *,'2 = by slice number'
  print *,'3 = by elevation of topography (for surface of model only)'
  print *,'4 = random color to show MPI slices'
  print *,'any other value=exit'
  print *
  print *,'enter value:'
  read(5,*) icolor
  if(icolor<1 .or. icolor >4) stop 'exiting...'
  if(icolor == 3 .and. ivalue /= 2) stop 'color by elevation of topography is for surface of model only'

  print *
  print *,'1 = material property by doubling flag'
  print *,'2 = by slice number'
  print *,'any other value=exit'
  print *
  print *,'enter value:'
  read(5,*) imaterial
  if(imaterial < 1 .or. imaterial > 2) stop 'exiting...'

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

! user can specify a range of processors here, enter 0 and -1 for all procs
  print *
  print *,'enter first proc (proc numbers start at 0) = '
  read(5,*) proc_p1
  if(proc_p1 < 0) proc_p1 = 0
  if(proc_p1 > NPROC-1) proc_p1 = NPROC-1

  print *,'enter last proc (enter -1 for all procs) = '
  read(5,*) proc_p2
  if(proc_p2 == -1) proc_p2 = NPROC-1
  if(proc_p2 < 0) proc_p2 = 0
  if(proc_p2 > NPROC-1) proc_p2 = NPROC-1

! set interval to maximum if user input is not correct
  if(proc_p1 <= 0) proc_p1 = 0
  if(proc_p2 < 0) proc_p2 = NPROC - 1

! set total number of points and elements to zero
  ntotpoin = 0
  ntotspec = 0

! initialize random colors
  allocate(random_colors(0:NPROC-1))
  do iproc=0,NPROC-1
    call random_number(random_val)
    ival_color = nint(random_val*NPROC)
    if(ival_color < 0) ival_color = 0
    if(ival_color > NPROC-1) ival_color = NPROC-1
    random_colors(iproc) = ival_color
  enddo

! loop on the selected range of processors
  do iproc = proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  if(ivalue == 1) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointsfaces.txt',status='old',action='read')
  else if(ivalue == 2) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='old',action='read')
  endif

  read(10,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'
  ntotpoin = ntotpoin + npoin
  close(10)

  if(ivalue == 1) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementsfaces.txt',status='old',action='read')
  else if(ivalue == 2) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementssurface.txt',status='old',action='read')
  endif

  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'
  ntotspec = ntotspec + nspec
  close(10)

  enddo

  print *
  print *,'There is a total of ',ntotspec,' elements in all the slices'
  print *,'There is a total of ',ntotpoin,' points in all the slices'
  print *

  ntotpoinAVS_DX = ntotpoin
  ntotspecAVS_DX = ntotspec

! use different name for surface and for slices
  if(USE_OPENDX) then
    open(unit=11,file=trim(OUTPUT_FILES)//'/DX_fullmesh.dx',status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ntotpoinAVS_DX,' data follows'
  else
    open(unit=11,file=trim(OUTPUT_FILES)//'/AVS_fullmesh.inp',status='unknown')
  endif

! write AVS or DX header with element data or point data
  if(.not. USE_OPENDX) then
    if(ivalue == 2 .and. icolor == 3) then
      write(11,*) ntotpoinAVS_DX,' ',ntotspecAVS_DX,' 1 0 0'
    else
      write(11,*) ntotpoinAVS_DX,' ',ntotspecAVS_DX,' 0 1 0'
    endif
  endif

! ************* generate points ******************

! set global point offset to zero
  iglobpointoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  if(ivalue == 1) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointsfaces.txt',status='old',action='read')
  else if(ivalue == 2) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='old',action='read')
  endif

  read(10,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'

! read local points in this slice and output global AVS or DX points
  do ipoin=1,npoin
      read(10,*) numpoin,xval,yval,zval
      if(numpoin /= ipoin) stop 'incorrect point number'
! write to AVS or DX global file with correct offset
      if(USE_OPENDX) then
        write(11,*) sngl(xval),' ',sngl(yval),' ',sngl(zval*zscaling)
      else
!!        write(11,*) numpoin + iglobpointoffset,' ',sngl(xval),' ',sngl(yval),' ',sngl(zval*zscaling)
!! XXX
 if(zval < 0.) then
        write(11,*) numpoin + iglobpointoffset,' ',sngl(xval),' ',sngl(yval),' ',sngl(zval*zscaling)
else
        write(11,*) numpoin + iglobpointoffset,' ',sngl(xval),' ',sngl(yval),' ',' 0'
endif
      endif

  enddo

  iglobpointoffset = iglobpointoffset + npoin

  close(10)

  enddo

! ************* generate elements ******************

! get source information for frequency for number of points per lambda
  print *,'reading source duration from the CMTSOLUTION file'
  call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,NSOURCES)

! set global element and point offsets to zero
  iglobpointoffset = 0
  iglobelemoffset = 0
  maxdoubling = -1

  if(USE_OPENDX) &
    write(11,*) 'object 2 class array type int rank 1 shape 4 items ',ntotspecAVS_DX,' data follows'

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  if(ivalue == 1) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementsfaces.txt',status='old',action='read')
    open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpointsfaces.txt',status='old',action='read')
  else if(ivalue == 2) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementssurface.txt',status='old',action='read')
    open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='old',action='read')
  endif

  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'

  read(12,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4
  if(numelem /= ispec) stop 'incorrect element number'
! compute max of the doubling flag
  maxdoubling = max(maxdoubling,idoubling)

! assign material property (which can be filtered later in AVS_DX)
  if(imaterial == 1) then
    imatprop = idoubling
  else if(imaterial == 2) then
    imatprop = iproc
  else
    stop 'invalid code for material property'
  endif

! write to AVS or DX global file with correct offset

! quadrangles (2-D)
      iglob1 = iglob1 + iglobpointoffset
      iglob2 = iglob2 + iglobpointoffset
      iglob3 = iglob3 + iglobpointoffset
      iglob4 = iglob4 + iglobpointoffset

! in the case of OpenDX, node numbers start at zero
! in the case of AVS, node numbers start at one
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
      if(USE_OPENDX) then
        write(11,"(i6,1x,i6,1x,i6,1x,i6)") iglob1-1,iglob4-1,iglob2-1,iglob3-1
      else
        write(11,"(i6,1x,i3,' quad ',i6,1x,i6,1x,i6,1x,i6)") numelem + iglobelemoffset,imatprop,iglob1,iglob2,iglob3,iglob4
      endif

  enddo

  iglobelemoffset = iglobelemoffset + nspec
  iglobpointoffset = iglobpointoffset + npoin

  close(10)
  close(12)

  enddo

! ************* generate data values ******************

! output AVS or DX header for data
  if(USE_OPENDX) then
    write(11,*) 'attribute "element type" string "quads"'
    write(11,*) 'attribute "ref" string "positions"'
    if(ivalue == 2 .and. icolor == 3) then
      write(11,*) 'object 3 class array type float rank 0 items ',ntotpoinAVS_DX,' data follows'
    else
      write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'
    endif
  else
    write(11,*) '1 1'
    write(11,*) 'Zcoord, meters'
  endif

!!!!
!!!! ###### element data in most cases
!!!!
  if(ivalue /= 2 .or. icolor /= 3) then

! set global element and point offsets to zero
  iglobelemoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  if(ivalue == 1) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementsfaces.txt',status='old',action='read')
  else if(ivalue == 2) then
    open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementssurface.txt',status='old',action='read')
  endif

  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4
      if(numelem /= ispec) stop 'incorrect element number'

! data is either the slice number or the mesh doubling region flag
      if(icolor == 1) then
        val_color = dble(idoubling)
      else if(icolor == 2) then
        val_color = dble(iproc)
      else if(icolor == 4) then
        val_color = dble(random_colors(iproc))
      else
        stop 'incorrect coloring code'
      endif

! write to AVS or DX global file with correct offset
      if(USE_OPENDX) then
        write(11,*) sngl(val_color)
      else
        write(11,*) numelem + iglobelemoffset,' ',sngl(val_color)
      endif
  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo

!!!!
!!!! ###### point data if surface colored according to topography
!!!!
  else

! set global point offset to zero
  iglobpointoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointssurface.txt',status='old',action='read')

  read(10,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'

! read local points in this slice and output global AVS or DX points
  do ipoin=1,npoin
      read(10,*) numpoin,xval,yval,zval
      if(numpoin /= ipoin) stop 'incorrect point number'
! write to AVS or DX global file with correct offset
      if(USE_OPENDX) then
        write(11,*) sngl(zval)
      else
        write(11,*) numpoin + iglobpointoffset,' ',sngl(zval)
      endif

  enddo

  iglobpointoffset = iglobpointoffset + npoin

  close(10)

  enddo

  endif     ! end test if element data or point data

! define OpenDX field
  if(USE_OPENDX) then
    if(ivalue == 2 .and. icolor == 3) then
      write(11,*) 'attribute "dep" string "positions"'
    else
      write(11,*) 'attribute "dep" string "connections"'
    endif
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'
  endif

  close(11)

  print *
  print *,'maximum value of doubling flag in all slices = ',maxdoubling
  print *

!
! create an AVS or DX file with the source and the receivers as well
!

  if(USE_OPENDX) then

    print *
    print *,'support for source and station file in OpenDX not added yet'
    print *
    stop 'warning: only partial support for OpenDX in current version (mesh ok, but no source)'

  else

!   get source information
    print *,'reading position of the source from the CMTSOLUTION file'
    call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,NSOURCES)

!   the point for the source is put at the surface for clarity (depth ignored)
!   even slightly above to superimpose to real surface
!   also save quadrangle for AVS or DX representation of epicenter

    z_target_source = 2000.
    z_source_quad1 = 2000.
    z_source_quad2 = 2000.
    z_source_quad3 = 2000.
    z_source_quad4 = 2000.

    call utm_geo(long,lat,x_target_source,y_target_source,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

    x_source_quad1 = x_target_source
    y_source_quad1 = y_target_source

    x_source_quad2 = x_target_source + small_offset
    y_source_quad2 = y_target_source

    x_source_quad3 = x_target_source + small_offset
    y_source_quad3 = y_target_source + small_offset

    x_source_quad4 = x_target_source
    y_source_quad4 = y_target_source + small_offset

    ntotpoinAVS_DX = 2
    ntotspecAVS_DX = 1

    print *
    print *,'reading position of the receivers from DATA/STATIONS_FILTERED file'
    call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_FILTERED')

! get number of stations from receiver file
    open(unit=11,file=filtered_rec_filename,iostat=ios,status='old',action='read')
    nrec = 0
    do while(ios == 0)
      read(11,"(a)",iostat=ios) dummystring
      if(ios == 0) nrec = nrec + 1
    enddo
    close(11)

    print *,'There are ',nrec,' three-component stations'
    print *
    if(nrec < 1) stop 'incorrect number of stations read - need at least one'

    allocate(station_name(nrec))
    allocate(network_name(nrec))
    allocate(stlat(nrec))
    allocate(stlon(nrec))
    allocate(stele(nrec))
    allocate(stbur(nrec))

    allocate(x_target(nrec))
    allocate(y_target(nrec))
    allocate(z_target(nrec))

! loop on all the stations
    open(unit=11,file=filtered_rec_filename,status='old',action='read')
    do irec=1,nrec
      read(11,*) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)

! points for the receivers are put at the surface for clarity (depth ignored)
      call utm_geo(stlon(irec),stlat(irec),x_target(irec),y_target(irec),UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

      z_target(irec) = 2000.

    enddo

    close(11)

! duplicate source to have right color normalization in AVS_DX
  ntotpoinAVS_DX = ntotpoinAVS_DX + 2*nrec + 1
  ntotspecAVS_DX = ntotspecAVS_DX + nrec + 1

  open(unit=11,file=trim(OUTPUT_FILES)//'/AVS_source_receivers.inp',status='unknown')

! write AVS or DX header with element data
  write(11,*) ntotpoinAVS_DX,' ',ntotspecAVS_DX,' 0 1 0'

! add source and receivers (small AVS or DX lines)
  write(11,*) '1 ',sngl(x_target_source),' ',sngl(y_target_source),' ',sngl(z_target_source)
  write(11,*) '2 ',sngl(x_target_source+small_offset),' ', &
    sngl(y_target_source+small_offset),' ',sngl(z_target_source+small_offset)
  write(11,*) '3 ',sngl(x_target_source+small_offset),' ', &
    sngl(y_target_source+small_offset),' ',sngl(z_target_source+small_offset)
  do ir=1,nrec
    write(11,*) 4+2*(ir-1),' ',sngl(x_target(ir)),' ',sngl(y_target(ir)),' ',sngl(z_target(ir))
    write(11,*) 4+2*(ir-1)+1,' ',sngl(x_target(ir)+small_offset),' ', &
      sngl(y_target(ir)+small_offset),' ',sngl(z_target(ir)+small_offset)
  enddo

! add source and receivers (small AVS or DX lines)
  write(11,*) '1 1 line 1 2'
  do ir=1,nrec
    write(11,*) ir+1,' 1 line ',4+2*(ir-1),' ',4+2*(ir-1)+1
  enddo
! duplicate source to have right color normalization in AVS_DX
  write(11,*) ir+1,' 1 line 1 3'

! output AVS or DX header for data
  write(11,*) '1 1'
  write(11,*) 'Zcoord, meters'

! add source and receiver data
  write(11,*) '1 1.'
  do ir=1,nrec
    write(11,*) ir+1,' 255.'
  enddo
! duplicate source to have right color normalization in AVS_DX
  write(11,*) ir+1,' 120.'

  close(11)

! create a file with the epicenter only, represented as a quadrangle

  open(unit=11,file=trim(OUTPUT_FILES)//'/AVS_epicenter.inp',status='unknown')

! write AVS or DX header with element data
  write(11,*) '4 1 0 1 0'

! add source and receivers (small AVS or DX lines)
  write(11,*) '1 ',sngl(x_source_quad1),' ',sngl(y_source_quad1),' ',sngl(z_source_quad1)
  write(11,*) '2 ',sngl(x_source_quad2),' ',sngl(y_source_quad2),' ',sngl(z_source_quad2)
  write(11,*) '3 ',sngl(x_source_quad3),' ',sngl(y_source_quad3),' ',sngl(z_source_quad3)
  write(11,*) '4 ',sngl(x_source_quad4),' ',sngl(y_source_quad4),' ',sngl(z_source_quad4)

! create a element for the source, some labels and element data
  write(11,*) '1 1 quad 1 2 3 4'
  write(11,*) '1 1'
  write(11,*) 'Zcoord, meters'
  write(11,*) '1 1.'

  close(11)

  endif

  end program combine_AVS_DX

