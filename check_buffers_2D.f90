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

! code to check that all the internal MPI buffers are okay along xi and eta
! we compare the coordinates of the points in the buffers

  program check_buffers_2D

  implicit none

  include "constants.h"

  integer ithisproc,iotherproc

  integer ipoin

  integer npoin2d_xi_save,npoin2d_xi_mesher,npoin2d_xi
  integer npoin2d_eta_save,npoin2d_eta_mesher,npoin2d_eta

! for addressing of the slices
  integer iproc_xi,iproc_eta,iproc
  integer iproc_read
  integer, dimension(:,:), allocatable :: addressing

  double precision diff

! 2-D addressing and buffers for summation between slices
  integer, dimension(:), allocatable :: iboolleft_xi,iboolright_xi, &
    iboolleft_eta,iboolright_eta

! coordinates of the points to compare
  double precision, dimension(:), allocatable :: xleft_xi,yleft_xi,zleft_xi, &
     xright_xi,yright_xi,zright_xi,xleft_eta,yleft_eta,zleft_eta, &
     xright_eta,yright_eta,zright_eta

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

  character(len=256) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

! now this is for all the regions
  integer NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB

! processor identification
  character(len=256) prname,prname_other

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers along xi and eta'
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
  print *,'There are ',NPROC_XI,' slices along xi'
  print *,'There are ',NPROC_ETA,' slices along eta'
  print *

! dynamic memory allocation for arrays
  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1))

! open file with global slice number addressing
  print *,'reading slice addressing'
  open(unit=34,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read')
  do iproc = 0,NPROC-1
      read(34,*) iproc_read,iproc_xi,iproc_eta
      if(iproc_read /= iproc) stop 'incorrect slice number read'
      addressing(iproc_xi,iproc_eta) = iproc
  enddo
  close(34)

! dynamic memory allocation for arrays
  allocate(iboolleft_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(iboolright_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(iboolleft_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(iboolright_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(xleft_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(yleft_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(zleft_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(xright_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(yright_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(zright_xi(NPOIN2DMAX_XMIN_XMAX))
  allocate(xleft_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(yleft_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(zleft_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(xright_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(yright_eta(NPOIN2DMAX_YMIN_YMAX))
  allocate(zright_eta(NPOIN2DMAX_YMIN_YMAX))

! double loop on NPROC_XI and NPROC_ETA
  do iproc_eta=0,NPROC_ETA-1

  print *,'checking row ',iproc_eta

  do iproc_xi=0,NPROC_XI-2

  print *,'checking slice ixi = ',iproc_xi,' in that row'

  ithisproc = addressing(iproc_xi,iproc_eta)
  iotherproc = addressing(iproc_xi+1,iproc_eta)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,LOCAL_PATH,NPROC,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_xi of this slice
  write(*,*) 'reading MPI buffer iboolright_xi slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 360  continue
  read(34,*) iboolright_xi(npoin2D_xi), &
              xright_xi(npoin2D_xi),yright_xi(npoin2D_xi),zright_xi(npoin2D_xi)
  if(iboolright_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 360
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolright_xi slice ',ithisproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolright_xi read'
  endif
  close(34)

! save to compare to other side
  npoin2D_xi_save = npoin2D_xi

! read iboolleft_xi of other slice
  write(*,*) 'reading MPI buffer iboolleft_xi slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 350  continue
  read(34,*) iboolleft_xi(npoin2D_xi), &
              xleft_xi(npoin2D_xi),yleft_xi(npoin2D_xi),zleft_xi(npoin2D_xi)
  if(iboolleft_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 350
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolleft_xi slice ',iotherproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolleft_xi read'
  endif
  close(34)

  if(npoin2D_xi_save == npoin2D_xi) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_xi
      diff = dmax1(dabs(xleft_xi(ipoin)-xright_xi(ipoin)), &
       dabs(yleft_xi(ipoin)-yright_xi(ipoin)),dabs(zleft_xi(ipoin)-zright_xi(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft_xi(ipoin),iboolright_xi(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo


! double loop on NPROC_XI and NPROC_ETA
  do iproc_xi=0,NPROC_XI-1

  print *,'checking row ',iproc_xi

  do iproc_eta=0,NPROC_ETA-2

  print *,'checking slice ieta = ',iproc_eta,' in that row'

  ithisproc = addressing(iproc_xi,iproc_eta)
  iotherproc = addressing(iproc_xi,iproc_eta+1)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,LOCAL_PATH,NPROC,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,LOCAL_PATH,NPROC,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_eta of this slice
  write(*,*) 'reading MPI buffer iboolright_eta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 460  continue
  read(34,*) iboolright_eta(npoin2D_eta), &
              xright_eta(npoin2D_eta),yright_eta(npoin2D_eta),zright_eta(npoin2D_eta)
  if(iboolright_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 460
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolright_eta slice ',ithisproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolright_eta read'
  endif
  close(34)

! save to compare to other side
  npoin2D_eta_save = npoin2D_eta

! read iboolleft_eta of other slice
  write(*,*) 'reading MPI buffer iboolleft_eta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 450  continue
  read(34,*) iboolleft_eta(npoin2D_eta), &
              xleft_eta(npoin2D_eta),yleft_eta(npoin2D_eta),zleft_eta(npoin2D_eta)
  if(iboolleft_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 450
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolleft_eta slice ',iotherproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolleft_eta read'
  endif
  close(34)

  if(npoin2D_eta_save == npoin2D_eta) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_eta
      diff = dmax1(dabs(xleft_eta(ipoin)-xright_eta(ipoin)), &
       dabs(yleft_eta(ipoin)-yright_eta(ipoin)),dabs(zleft_eta(ipoin)-zright_eta(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft_eta(ipoin),iboolright_eta(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_2D

