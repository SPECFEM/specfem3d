!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! combine mesh quality data files to check the mesh
! displays statistics on mesh quality
! and creates an AVS or DX file showing a given range of elements

!! DK DK
!! DK DK this routine could be improved:
!! DK DK
!! DK DK  - add full OpenDX support
!! DK DK  - debug aspect ratio
!! DK DK  - add mean in addition to min and max of ratios
!! DK DK

  program mesh_quality_AVS_DX

  implicit none

  include "constants.h"

  integer iproc,nspec,npoin
  integer ispec
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer ipoin,numpoin,iglobpointoffset,ntotpoin,ntotspec
  integer numelem,iglobelemoffset
  integer idoubling,iformat
  integer ntotpoinAVS_DX,ntotspecAVS_DX

  double precision xval,yval,zval

! processor identification
  character(len=150) prname

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_XI,NEX_ETA, &
             NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

! for all the regions
  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

! for quality of mesh
  logical, dimension(:), allocatable :: mask_ibool
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision skewness_AVS_DX_min,skewness_AVS_DX_max

! for stability and number of points per wavelength
  double precision stability,points_per_wavelength
  double precision stability_min,points_per_wavelength_min
  double precision stability_max,points_per_wavelength_max

! for histogram
  integer, parameter :: NCLASS = 20
  integer classes_skewness(0:NCLASS-1)
  integer iclass
  double precision current_percent,total_percent

  integer proc_p1,proc_p2

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all mesh quality files for slices'
  print *

  print *,'1 = create files in AVS UCD format'
  print *,'2 = create files in OpenDX format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
!!!!!! DK DK  read(5,*) iformat
!! DK DK impose AVS format for now, OpenDX format not implemented yet
  iformat = 1
  if(iformat<1 .or. iformat>2) stop 'exiting...'

! read range of skewness used for elements
  print *,'enter minimum skewness for AVS or DX (between 0. and 1.):'
  read(5,*) skewness_AVS_DX_min
  if(skewness_AVS_DX_min < 0.d0) skewness_AVS_DX_min = 0.d0
  if(skewness_AVS_DX_min > 0.99999d0) skewness_AVS_DX_min = 0.99999d0

  print *,'enter maximum skewness for AVS or DX (between 0. and 1.):'
  read(5,*) skewness_AVS_DX_max
  if(skewness_AVS_DX_max < 0.d0) skewness_AVS_DX_max = 0.d0
  if(skewness_AVS_DX_max > 0.99999d0) skewness_AVS_DX_max = 0.99999d0

  if(skewness_AVS_DX_min > skewness_AVS_DX_max) stop 'incorrect skewness range'

  print *
  print *,'reading parameter file'
  print *

! read the parameter file
  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

  if(.not. SAVE_MESH_FILES) stop 'AVS or DX files were not saved by the mesher'

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,USE_REGULAR_MESH)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  print *
  print *,'There are ',NPROC,' slices numbered from 0 to ',NPROC-1
  print *

! use all the slices to determine correct statistics
  proc_p1 = 0
  proc_p2 = NPROC - 1

! set total number of points and elements to zero
  ntotpoin = 0
  ntotspec = 0

! loop on the selected range of processors
  do iproc = proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old',action='read')
  read(10,*) npoin
  print *,'There are ',npoin,' global AVS or DX points in the slice'
  ntotpoin = ntotpoin + npoin
  close(10)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old',action='read')
  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'
  ntotspec = ntotspec + nspec
  close(10)

  enddo

  print *
  print *,'There is a total of ',ntotspec,' elements in all the slices'
  print *,'There is a total of ',ntotpoin,' points in all the slices'
  print *

! ************* compute min and max of skewness and ratios ******************

! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL
  stability_min = + HUGEVAL
  points_per_wavelength_min = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL
  stability_max = - HUGEVAL
  points_per_wavelength_max = - HUGEVAL

! set global element and point offsets to zero
  iglobelemoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

  print *,'Reading slice ',iproc

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old',action='read')

  read(10,*) nspec
  print *,'There are ',nspec,' AVS or DX elements in the slice'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength
      if(numelem /= ispec) stop 'incorrect element number'

! mulitply stability number by time step
      stability = stability * DT

! compute minimum and maximum of quality numbers
      equiangle_skewness_min = dmin1(equiangle_skewness_min,equiangle_skewness)
      edge_aspect_ratio_min = dmin1(edge_aspect_ratio_min,edge_aspect_ratio)
      diagonal_aspect_ratio_min = dmin1(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
      stability_min = dmin1(stability_min,stability)
      points_per_wavelength_min = dmin1(points_per_wavelength_min,points_per_wavelength)

      equiangle_skewness_max = dmax1(equiangle_skewness_max,equiangle_skewness)
      edge_aspect_ratio_max = dmax1(edge_aspect_ratio_max,edge_aspect_ratio)
      diagonal_aspect_ratio_max = dmax1(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
      stability_max = dmax1(stability_max,stability)
      points_per_wavelength_max = dmax1(points_per_wavelength_max,points_per_wavelength)

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions'
  print *
  print *,'equiangle skewness: 0. perfect  1. bad'
  print *,'skewness max deviation angle: 0. perfect  90. bad'
  print *,'skewness min mesh angle: 90. perfect  0. bad'
  print *,'edge aspect ratio: 1. perfect  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect  above 1. gives stretching factor'
  print *,'------------'

  print *
  print *,'equiangle skewness max = ',equiangle_skewness_max
  print *,'equiangle skewness min = ',equiangle_skewness_min
  print *
  print *,'skewness max deviation angle = ',90.*equiangle_skewness_max
  print *,'skewness min mesh angle = ',90.*(1. - equiangle_skewness_max)
  print *
  print *,'edge aspect ratio max = ',edge_aspect_ratio_max
  print *,'edge aspect ratio min = ',edge_aspect_ratio_min
  print *
  print *,'diagonal aspect ratio max = ',diagonal_aspect_ratio_max
  print *,'diagonal aspect ratio min = ',diagonal_aspect_ratio_min
  print *
  print *,'stability max = ',stability_max
  print *,'stability min = ',stability_min
  print *
  print *,'points per wavelength max = ',points_per_wavelength_max
  print *,'points per wavelength min = ',points_per_wavelength_min
  print *

! create statistics about mesh quality

  print *
  print *,'creating histogram and statistics of mesh quality - reading mesh data files'
  print *

! erase histogram of skewness
  classes_skewness(:) = 0

! erase number of elements belonging to skewness range for AVS_DX
  ntotspecAVS_DX = 0

! set global element and point offsets to zero
  iglobelemoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old',action='read')

  read(10,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength
      if(numelem /= ispec) stop 'incorrect element number'

! store skewness in histogram
    iclass = int(equiangle_skewness * dble(NCLASS))
    if(iclass < 0) iclass = 0
    if(iclass > NCLASS-1) iclass = NCLASS-1
    classes_skewness(iclass) = classes_skewness(iclass) + 1

! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. &
       equiangle_skewness <= skewness_AVS_DX_max) ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo

! create histogram of skewness and save in Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file=trim(OUTPUT_FILES)//'/mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_skewness(iclass))/dble(ntotspec)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewness(iclass),' ',sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

! create script for Gnuplot histogram file
  open(unit=14,file=trim(OUTPUT_FILES)//'/plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term x11'
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASS)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  print *
  print *,'total number of elements = ',ntotspec
  print *,'total percentage = ',total_percent,' %'
  print *

! display warning if maximum skewness is too high
  if(equiangle_skewness_max >= 0.80d0) then
    print *
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, mesh is bad (max skewness >= 0.80)'
    print *,'*********************************************'
    print *,'*********************************************'
    print *
  endif


! ************* create AVS or DX file with elements in a certain range of skewness

  print *
  print *,'creating AVS or DX file with subset of elements in skewness range'
  print *,'between ',skewness_AVS_DX_min,' and ',skewness_AVS_DX_max
  print *

! ************* count number of points without multiples  ******************

! set global element and point offsets to zero
  iglobpointoffset = 0
  iglobelemoffset = 0

! allocate flag to remove multiples
  allocate(mask_ibool(ntotpoin))
  mask_ibool(:) = .false.

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old',action='read')
  open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old',action='read')
  open(unit=14,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old',action='read')

  read(10,*) nspec
  read(12,*) npoin
  read(14,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
    read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
    if(numelem /= ispec) stop 'incorrect element number'

    read(14,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength
    if(numelem /= ispec) stop 'incorrect element number'

! check if element belongs to requested skewness range
! and flag all the points to remove multiples
      iglob1 = iglob1 + iglobpointoffset
      iglob2 = iglob2 + iglobpointoffset
      iglob3 = iglob3 + iglobpointoffset
      iglob4 = iglob4 + iglobpointoffset
      iglob5 = iglob5 + iglobpointoffset
      iglob6 = iglob6 + iglobpointoffset
      iglob7 = iglob7 + iglobpointoffset
      iglob8 = iglob8 + iglobpointoffset
      if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) then
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob5) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      endif

  enddo

  iglobelemoffset = iglobelemoffset + nspec
  iglobpointoffset = iglobpointoffset + npoin

  close(10)
  close(12)
  close(14)

  enddo

  if(ntotspecAVS_DX == 0) stop 'no elements in skewness range, no file created'

! count number of independent points
  ntotpoinAVS_DX = count(mask_ibool(:))

  open(unit=11,file='AVS_meshquality.inp',status='unknown')

! write AVS or DX header with element data
  write(11,*) ntotpoinAVS_DX,' ',ntotspecAVS_DX,' 0 1 0'

! ************* generate points ******************

! set global point offset to zero
  iglobpointoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old',action='read')
  read(10,*) npoin

! read local points in this slice and output global AVS or DX points
  do ipoin=1,npoin
      read(10,*) numpoin,xval,yval,zval
      if(numpoin /= ipoin) stop 'incorrect point number'
! write to AVS or DX global file with correct offset if point has been selected
      if(mask_ibool(numpoin + iglobpointoffset)) &
        write(11,*) numpoin + iglobpointoffset,' ',sngl(xval),' ',sngl(yval),' ',sngl(zval)

  enddo

  iglobpointoffset = iglobpointoffset + npoin

  close(10)

  enddo

! ************* generate elements ******************

! set global element and point offsets to zero
  iglobpointoffset = 0
  iglobelemoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelements.txt',status='old',action='read')
  open(unit=12,file=prname(1:len_trim(prname))//'AVS_DXpoints.txt',status='old',action='read')
  open(unit=14,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old',action='read')

  read(10,*) nspec
  read(12,*) npoin
  read(14,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
    read(10,*) numelem,idoubling,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
    if(numelem /= ispec) stop 'incorrect element number'

    read(14,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength
    if(numelem /= ispec) stop 'incorrect element number'

! write to AVS or DX global file with correct offset for hexahedra (3-D)
! check if element belongs to requested skewness range
      iglob1 = iglob1 + iglobpointoffset
      iglob2 = iglob2 + iglobpointoffset
      iglob3 = iglob3 + iglobpointoffset
      iglob4 = iglob4 + iglobpointoffset
      iglob5 = iglob5 + iglobpointoffset
      iglob6 = iglob6 + iglobpointoffset
      iglob7 = iglob7 + iglobpointoffset
      iglob8 = iglob8 + iglobpointoffset
      if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) &
        write(11,"(i6,' 1 hex ',i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
            numelem + iglobelemoffset,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  enddo

  iglobelemoffset = iglobelemoffset + nspec
  iglobpointoffset = iglobpointoffset + npoin

  close(10)
  close(12)
  close(14)

  enddo

! ************* generate element data values ******************

! output AVS or DX header for data
  write(11,*) '1 1'
  write(11,*) 'Zcoord, meters'

! set global element and point offsets to zero
  iglobelemoffset = 0

! loop on the selected range of processors
  do iproc=proc_p1,proc_p2

! create the name for the database of the current slide
  call create_serial_name_database(prname,iproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXmeshquality.txt',status='old',action='read')

  read(10,*) nspec

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec
      read(10,*) numelem,equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,points_per_wavelength
      if(numelem /= ispec) stop 'incorrect element number'

! write skewness data to AVS or DX global file with correct offset
! scale skewness to [0:255] for AVS or DX color palette
! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) &
      write(11,*) numelem + iglobelemoffset,' ',255.*sngl(equiangle_skewness)

  enddo

  iglobelemoffset = iglobelemoffset + nspec

  close(10)

  enddo

! close AVS or DX file
  close(11)

  print *
  print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX skewness range ',skewness_AVS_DX_min,skewness_AVS_DX_max
  print *

  end program mesh_quality_AVS_DX

