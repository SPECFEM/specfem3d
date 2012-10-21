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

! read a 2D or 3D CUBIT mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Pau, France, March 2009.

!! DK DK
!! DK DK this routine could be improved by computing the mean in addition to min and max of ratios
!! DK DK
!! DK DK also, the particular treatment for sub-blocks in Gocad files, e.g.,
!! DK DK if(IGNORE_OTHER_HEADERS .and. cubit_mesh_file == 'REGOLITE_only_no_fractures_2D_in_meters.inp' &
!! DK DK                  .and. i == 28429)
!! DK DK is not general at all and should be rewritten
!! DK DK

  program check_mesh_quality_CUBIT_Abaqus

  implicit none

  include "constants.h"

!------------------------------------------------------------------------------------------------
! EDIT YOUR PARAMETERS BELOW HERE

! number of points and of hex or quad elements
! number of points of a hex or quad element

! example: layered_halfspace
!                 Cubit -> File -> Export... Abacus (*.inp)
!                                         ( block ids: 1 2 3 ) volumes only
!                                         (optional: uncheck 'Export Using Cubit IDs' to have element IDs in increasing order)
  character(len=100), parameter :: cubit_mesh_file = 'examples/layered_halfspace/layered_halfspace_mesh.inp'
  integer, parameter :: NPOIN = 76819                    ! number of nodes
  integer, parameter :: NSPEC = 70200                    ! number of elements (only volumes, i.e. block ids 1,2,3 )
  integer, parameter :: NGNOD = 8                        ! hexahedral elements
  logical, parameter :: IGNORE_OTHER_HEADERS = .false.
  double precision, parameter :: delta_t = 0.005         ! arbitrary, initial guess
  double precision, parameter :: VP_MAX = 7500.d0        ! maximum vp in volume block id 3

!------------------------------------------------------------------------------------------------

  double precision, dimension(NPOIN) :: x,y,z

  integer, dimension(NGNOD,NSPEC) :: ibool

  integer :: i,ispec,iread,iformat,ispec_min_edge_length,ispec_max_edge_length, &
             ispec_begin,ispec_end,ispec_to_output,ier

  double precision :: xtmp,ytmp,ztmp
  integer :: n1,n2,n3,n4,n5,n6,n7,n8

! for quality of mesh
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision :: equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision :: equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision :: skewness_AVS_DX_min,skewness_AVS_DX_max,distance_min,distance_max
  double precision :: distmin,distmax

! for stability
  double precision :: stability,stability_min,stability_max,max_CFL_stability_limit

! for histogram
  integer, parameter :: NCLASS = 20
  integer classes_skewness(0:NCLASS-1)
  integer :: iclass
  double precision :: current_percent,total_percent
! to export elements that have a certain skewness range to OpenDX
  integer :: ntotspecAVS_DX
  logical :: USE_OPENDX

  character(len=256):: line

  if(NGNOD /= 8) then
    print *,'error: check_mesh_quality_CUBIT_Abaqus only supports NGNOD == 8 for now'
    stop 'thus if NGNOD == 27, just run the solver without checking the mesh with this program'
  endif

  print *
  print *,'1 = output elements above a certain skewness threshold in OpenDX format'
  print *,'2 = output a given element in OpenDX format'
  print *,'3 = do not output any OpenDX file'
  print *
  print *,'enter value:'
  read(5,*) iformat

  if(iformat < 1 .or. iformat > 3) stop 'exiting...'

  if(iformat == 1 .or. iformat == 2) then
    USE_OPENDX = .true.
  else
    USE_OPENDX = .false.
  endif

  if(USE_OPENDX) then

    if(iformat == 1) then

      ! read range of skewness used for elements
      print *,'enter minimum skewness for OpenDX (between 0. and 0.99):'
      read(5,*) skewness_AVS_DX_min
      if(skewness_AVS_DX_min < 0.d0) skewness_AVS_DX_min = 0.d0
      if(skewness_AVS_DX_min > 0.99999d0) skewness_AVS_DX_min = 0.99999d0

      !!!!!!!!  print *,'enter maximum skewness for OpenDX (between 0. and 1.):'
      !!!!!!!!!!!!!  read(5,*) skewness_AVS_DX_max
      skewness_AVS_DX_max = 0.99999d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(skewness_AVS_DX_max < 0.d0) skewness_AVS_DX_max = 0.d0
      if(skewness_AVS_DX_max > 0.99999d0) skewness_AVS_DX_max = 0.99999d0

      if(skewness_AVS_DX_min > skewness_AVS_DX_max) stop 'incorrect skewness range'

    else
      print *,'enter the element number to output in OpenDX format between 1 and ',NSPEC
      read(5,*) ispec_to_output
      if(ispec_to_output < 1 .or. ispec_to_output > NSPEC) stop 'incorrect element number to output'
    endif

  endif

! read the mesh
  print *
  print *,'start reading the CUBIT file: ',cubit_mesh_file(1:len_trim(cubit_mesh_file))
  print *,'  number of points: ',NPOIN
  print *,'  number of elements: ',NSPEC
  print *

  open(unit=10,file=cubit_mesh_file,status='unknown',action='read')

! skip the header:
!     *HEADING
!     cubit(M3D/examples/layered_halfspace/layered_halfspace_mesh.inp): 01/31/2011: 10
!     version: 12.2
  read(10,*)
  read(10,*)
  read(10,*)

! read the points / nodes section:
!   **
!   ********************************** P A R T S **********************************
!   *PART, NAME=Part-Default
!   **
!   ********************************** N O D E S **********************************
!   *NODE, NSET=ALLNODES
  iread = 0
  x(:) = 0.d0
  y(:) = 0.d0
  z(:) = 0.d0
  do i = 1,NPOIN

    ! reads in text line
    read(10,'(a256)',iostat=ier) line
    if(ier /= 0 ) then
      print *,'error read line:',i
      stop 'error read points'
    endif

    ! checks if line is a comment line (starts with *), and reads until it finds a non-comment line
    do while ( line(1:1) == "*" )
      ! skips comment line and goes to next line
      print*,'  comment:',trim(line)
      read(10,'(a256)',iostat=ier) line
      if(ier /= 0 ) then
        print *,'error read non-comment line:',i
        stop 'error read points'
      endif
    enddo

    ! gets node ID and position
    read(line,*,iostat=ier) iread,xtmp,ytmp,ztmp

    ! checks
    if(ier /= 0 ) then
      print *,'error point read:',iread,i
      print*, 'line: ',trim(line)
      stop 'error read points from current line'
    endif
    ! checks if out-of-range
    if(iread < 1 .or. iread > NPOIN ) then
      print *,'error at i,iread = ',i,iread
      stop 'wrong ID input for a point'
    endif

    ! stores locations
    x(iread) = xtmp
    y(iread) = ytmp
    z(iread) = ztmp

  enddo
  print*
  print*,'points read: ',iread
  print*

! skip the header
  !read(10,*)
  read(10,'(a256)',iostat=ier) line
  if( line(1:1) /= "*" ) then
    print*,'  new line: ',trim(line)
    print*,'  not a header line, check the number of points NPOIN specified'
    stop 'error reading elements'
  endif

! read the elements:
!   **
!   ********************************** E L E M E N T S ****************************
!   *ELEMENT, TYPE=C3D8R, ELSET=elastic_1
  iread = 0
  ibool(:,:) = 0
  do i = 1,NSPEC

    ! reads in element connectivity
    if(NGNOD == 4) then

      ! quadrangles

      !! DK DK ignore other headers for 2D mesh of Eros with fractures, which has multiple material sets
      if(IGNORE_OTHER_HEADERS .and. cubit_mesh_file == 'eros_complexe_2d_regolite_fractures_modifie_in_meters.inp' &
                 .and. i == 5709) read(10,*)

      if(IGNORE_OTHER_HEADERS .and. cubit_mesh_file == 'REGOLITE_only_no_fractures_2D_in_meters.inp' &
                 .and. i == 28429) read(10,*)

      ! reads in line
      read(10,'(a256)',iostat=ier) line
      if(ier /= 0 ) then
        print *,'error read:',iread
        stop 'error read elements line'
      endif

      ! checks if line is a comment line (starts with *), and reads until it finds a non-comment line
      do while ( line(1:1) == "*" )
        ! skips comment line and goes to next line
        print*,'  comment: ',trim(line)
        read(10,'(a256)',iostat=ier) line
        if(ier /= 0 ) then
          print *,'error read:',i
          stop 'error read non-comment elements line'
        endif
      enddo

      ! gets element connection nodes
      read(line,*,iostat=ier) iread,n1,n2,n3,n4

      ! checks
      if(ier /= 0 ) then
        print *,'error read:',iread
        stop 'error read elements'
      endif
      ! requires that elements are in increasing order
      if(iread /= i) then
        print *,'error at i,iread = ',i,iread
        stop 'wrong input ID for an element'
      endif

      ! stores element nodes
      ibool(1,iread) = n1
      ibool(2,iread) = n2
      ibool(3,iread) = n3
      ibool(4,iread) = n4

    else if(NGNOD == 8) then

      ! hexahedra

      if(IGNORE_OTHER_HEADERS .and. cubit_mesh_file == 'rego3d_70_disp.inp' &
                 .and. i == 252929) read(10,*)

      ! reads in line
      read(10,'(a256)',iostat=ier) line
      if(ier /= 0 ) then
        print *,'error read:',iread
        stop 'error read elements line'
      endif

      ! checks if line is a comment line (starts with *), and reads until it finds a non-comment line
      do while ( line(1:1) == "*" )
        ! skips comment line and goes to next line
        print*,'  comment: ',trim(line)
        read(10,'(a256)',iostat=ier) line
        if(ier /= 0 ) then
          print *,'error read:',i
          stop 'error read non-comment elements line'
        endif
      enddo

      ! gets element connection nodes
      read(line,*,iostat=ier) iread,n1,n2,n3,n4,n5,n6,n7,n8

      ! checks
      if(ier /= 0 ) then
        print *,'error element read:',i
        print *,'line: ',trim(line)
        stop 'error read elements connectivity'
      endif
      if( iread < 1 .or. iread > NSPEC ) then
        print *,'error at i,iread = ',i,iread
        stop 'wrong input ID for an element'
      endif

      ! if we analyze only the second layer of the mesh and ignore the first, shift iread
      ! so that it conforms with i
      if(cubit_mesh_file == 'rego3d_70_disp_bedrock_only.inp') iread = iread - 252928

      ! stores element nodes
      ibool(1,iread) = n1
      ibool(2,iread) = n2
      ibool(3,iread) = n3
      ibool(4,iread) = n4
      ibool(5,iread) = n5
      ibool(6,iread) = n6
      ibool(7,iread) = n7
      ibool(8,iread) = n8

    endif

  enddo
  close(10)
  print*
  print*,'elements read:',iread
  print*

  print *,'done reading the CUBIT file'
  print *

  print *,'start computing the minimum and maximum edge size'

! ************* compute min and max of skewness and ratios ******************

! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL
  stability_min = + HUGEVAL
  distance_min = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL
  stability_max = - HUGEVAL
  distance_max = - HUGEVAL

  ispec_min_edge_length = -1
  ispec_max_edge_length = -1

! loop on all the elements
  do ispec = 1,NSPEC

    if(mod(ispec,100000) == 0) print *,'processed ',ispec,' elements out of ',NSPEC

    if(NGNOD == 4) then
      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    else
      call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    endif

! store element number in which the edge of minimum or maximum length is located
    if(distmin < distance_min) ispec_min_edge_length = ispec
    if(distmax > distance_max) ispec_max_edge_length = ispec

! compute minimum and maximum of quality numbers
    equiangle_skewness_min = min(equiangle_skewness_min,equiangle_skewness)
    edge_aspect_ratio_min = min(edge_aspect_ratio_min,edge_aspect_ratio)
    diagonal_aspect_ratio_min = min(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
    stability_min = min(stability_min,stability)
    distance_min = min(distance_min,distmin)

    equiangle_skewness_max = max(equiangle_skewness_max,equiangle_skewness)
    edge_aspect_ratio_max = max(edge_aspect_ratio_max,edge_aspect_ratio)
    diagonal_aspect_ratio_max = max(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
    stability_max = max(stability_max,stability)
    distance_max = max(distance_max,distmax)

  enddo
  print *,'done processing ',NSPEC,' elements out of ',NSPEC

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions:'
  print *
  print *,'equiangle skewness: 0. perfect,  1. bad'
  print *,'skewness max deviation angle: 0. perfect,  90. bad'
  print *,'edge aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'------------'

  print *
  print *,'minimum length of an edge in the whole mesh (m) = ',distance_min,' in element ',ispec_min_edge_length
  print *
  print *,'maximum length of an edge in the whole mesh (m) = ',distance_max,' in element ',ispec_max_edge_length
  print *
  print *,'max equiangle skewness = ',equiangle_skewness_max
! print *,'min equiangle skewness = ',equiangle_skewness_min
  print *
  print *,'max deviation angle from a right angle (90 degrees) is therefore = ',90.*equiangle_skewness_max
  print *
  print *,'worst angle in the mesh is therefore ',90.*(1. - equiangle_skewness_max)
  print *,'or ',180. - 90.*(1. - equiangle_skewness_max),' degrees'
  print *
  print *,'max edge aspect ratio = ',edge_aspect_ratio_max
! print *,'min edge aspect ratio = ',edge_aspect_ratio_min
  print *
  print *,'max diagonal aspect ratio = ',diagonal_aspect_ratio_max
! print *,'min diagonal aspect ratio = ',diagonal_aspect_ratio_min
  print *
  print *,'max stability = ',stability_max
  print *,'computed using VP_MAX = ',VP_MAX
! print *,'min stability = ',stability_min

! max stability CFL value is different in 2D and in 3D
  if(NGNOD == 8) then
    max_CFL_stability_limit = 0.48d0
  else if(NGNOD == 4) then
    max_CFL_stability_limit = 0.68d0
  else
    stop 'NGNOD must be 4 or 8'
  endif

  if(stability_max >= max_CFL_stability_limit) then
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, that value is above the upper CFL limit of ',max_CFL_stability_limit
    print *,'therefore the run should be unstable'
    print *,'*********************************************'
    print *,'*********************************************'
  else
    print *,'that value is below the upper CFL limit of ',max_CFL_stability_limit
    print *,'therefore the run should be stable'
  endif
  print *

! create statistics about mesh quality
  print *,'creating histogram and statistics of mesh quality'

! erase histogram of skewness
  classes_skewness(:) = 0

! loop on all the elements
  do ispec = 1,NSPEC

    if(NGNOD == 4) then
      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    else
      call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    endif

! store skewness in histogram
    iclass = int(equiangle_skewness * dble(NCLASS))
    if(iclass < 0) iclass = 0
    if(iclass > NCLASS-1) iclass = NCLASS-1
    classes_skewness(iclass) = classes_skewness(iclass) + 1

  enddo

! create histogram of skewness and save in Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file='mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_skewness(iclass))/dble(NSPEC)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewness(iclass),' ',sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

! create script for Gnuplot histogram file
  open(unit=14,file='plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term x11'
  write(14,*) '#set term gif'
  write(14,*) '#set output "mesh_quality_histogram.gif"'
  write(14,*)
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASS)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  print *
  print *,'total number of elements = ',NSPEC
  print *

! display warning if maximum skewness is too high
  if(equiangle_skewness_max >= 0.75d0) then
    print *
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, mesh is bad (max skewness >= 0.75)'
    print *,'*********************************************'
    print *,'*********************************************'
    print *
  endif

  if(total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    print *,'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  endif

! ************* create OpenDX file with elements in a certain range of skewness

  if(USE_OPENDX) then

  print *
  if(iformat == 1) then
    print *,'creating OpenDX file with subset of elements in skewness range'
    print *,'between ',skewness_AVS_DX_min,' and ',skewness_AVS_DX_max
  else
    print *,'creating OpenDX file with element #',ispec_to_output
  endif
  print *

! ************* count number of elements in skewness range *************

! erase number of elements belonging to skewness range for AVS_DX
  ntotspecAVS_DX = 0

! loop on all the elements
  if(iformat == 1) then

  do ispec = 1,NSPEC

    if(NGNOD == 4) then
      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    else
      call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    endif

! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) &
        ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  else
! outputing a single element
    ntotspecAVS_DX = 1
  endif

  if(ntotspecAVS_DX == 0) then
    stop 'no elements in skewness range, no file created'
  else if(iformat == 1) then
    print *
    print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX skewness range ',skewness_AVS_DX_min,skewness_AVS_DX_max
    print *
  endif

  open(unit=11,file='DX_mesh_quality.dx',status='unknown')

! ************* generate points ******************

! write OpenDX header
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NPOIN,' data follows'

! write all the points
  do i = 1,NPOIN
    write(11,*) sngl(x(i)),sngl(y(i)),sngl(z(i))
  enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape ',NGNOD,' items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  if(iformat == 1) then
    ispec_begin = 1
    ispec_end = NSPEC
  else
    ispec_begin = ispec_to_output
    ispec_end = ispec_to_output
  endif

  do ispec = ispec_begin,ispec_end

    if(NGNOD == 4) then
      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    else
      call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    endif

! check if element needs to be output
    if(iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) then
! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
      if(NGNOD == 4) then
        write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool(1,ispec)-1, ibool(4,ispec)-1, ibool(2,ispec)-1, ibool(3,ispec)-1
      else
        write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool(4,ispec)-1, ibool(1,ispec)-1, ibool(8,ispec)-1, ibool(5,ispec)-1, &
            ibool(3,ispec)-1, ibool(2,ispec)-1, ibool(7,ispec)-1, ibool(6,ispec)-1
      endif
      if(iformat == 1) print *,'element ',ispec,' belongs to the range and has skewness = ',sngl(equiangle_skewness)
    endif

  enddo

! ************* generate element data values ******************

! output OpenDX header for data
  if(NGNOD == 4) then
    write(11,*) 'attribute "element type" string "quads"'
  else
    write(11,*) 'attribute "element type" string "cubes"'
  endif
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  do ispec = ispec_begin,ispec_end

    if(NGNOD == 4) then
      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    else
      call create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)
    endif

! check if element needs to be output
    if(iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) &
    write(11,*) sngl(equiangle_skewness)

  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

! close OpenDX file
  close(11)

  endif

  end program check_mesh_quality_CUBIT_Abaqus

!
!=====================================================================
!

! create mesh quality data for a given 3D spectral element

  subroutine create_mesh_quality_data_3D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)

  implicit none

  include "constants.h"

  integer :: iface,icorner,ispec,NSPEC,NPOIN,NGNOD,i

  double precision, dimension(NPOIN) :: x,y,z

  integer, dimension(NGNOD,NSPEC) :: ibool

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  double precision vectorA_x,vectorA_y,vectorA_z
  double precision vectorB_x,vectorB_y,vectorB_z
  double precision norm_A,norm_B,angle_vectors
  double precision distmin,distmax,dist,dist1,dist2,dist3,dist4
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio

! for stability
  double precision :: stability,VP_MAX,delta_t

! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision, dimension(NGLL_MAX_STABILITY) :: percent_GLL

! topology of faces of cube for skewness
  integer faces_topo(6,6)

! store the corners of this element for the skewness routine
  do i = 1,NGNOD
    xelm(i) = x(ibool(i,ispec))
    yelm(i) = y(ibool(i,ispec))
    zelm(i) = z(ibool(i,ispec))
  enddo

! define percentage of smallest distance between GLL points for NGLL points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(1) = 100.d0
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

! check that the degree is not above the threshold for list of percentages
  if(NGLLX > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

! define topology of faces of cube for skewness

! face 1
  faces_topo(1,1) = 1
  faces_topo(1,2) = 2
  faces_topo(1,3) = 6
  faces_topo(1,4) = 5

! face 2
  faces_topo(2,1) = 2
  faces_topo(2,2) = 3
  faces_topo(2,3) = 7
  faces_topo(2,4) = 6

! face 3
  faces_topo(3,1) = 4
  faces_topo(3,2) = 3
  faces_topo(3,3) = 7
  faces_topo(3,4) = 8

! face 4
  faces_topo(4,1) = 1
  faces_topo(4,2) = 5
  faces_topo(4,3) = 8
  faces_topo(4,4) = 4

! face 5
  faces_topo(5,1) = 1
  faces_topo(5,2) = 2
  faces_topo(5,3) = 3
  faces_topo(5,4) = 4

! face 6
  faces_topo(6,1) = 5
  faces_topo(6,2) = 6
  faces_topo(6,3) = 7
  faces_topo(6,4) = 8

! define wraparound for angles for skewness calculation
  faces_topo(:,5) = faces_topo(:,1)
  faces_topo(:,6) = faces_topo(:,2)

! compute equiangle skewness (as defined in Fluent/Gambit manual)
! and compute edge aspect ratio using the corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     equiangle_skewness = - HUGEVAL

     do iface = 1,6
       do icorner = 1,4

! first vector of angle
         vectorA_x = xelm(faces_topo(iface,icorner)) - xelm(faces_topo(iface,icorner+1))
         vectorA_y = yelm(faces_topo(iface,icorner)) - yelm(faces_topo(iface,icorner+1))
         vectorA_z = zelm(faces_topo(iface,icorner)) - zelm(faces_topo(iface,icorner+1))

! second vector of angle
         vectorB_x = xelm(faces_topo(iface,icorner+2)) - xelm(faces_topo(iface,icorner+1))
         vectorB_y = yelm(faces_topo(iface,icorner+2)) - yelm(faces_topo(iface,icorner+1))
         vectorB_z = zelm(faces_topo(iface,icorner+2)) - zelm(faces_topo(iface,icorner+1))

! norm of vectors A and B
         norm_A = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
         norm_B = sqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

! angle formed by the two vectors
         angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

! compute equiangle skewness
         equiangle_skewness = max(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

! compute min and max size of an edge
         dist = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)

         distmin = min(distmin,dist)
         distmax = max(distmax,dist)

       enddo
     enddo

! compute edge aspect ratio
   edge_aspect_ratio = distmax / distmin

   stability = delta_t * VP_MAX / (distmin * percent_GLL(NGLLX))

! compute diagonal aspect ratio
   dist1 = sqrt((xelm(1) - xelm(7))**2 + (yelm(1) - yelm(7))**2 + (zelm(1) - zelm(7))**2)
   dist2 = sqrt((xelm(2) - xelm(8))**2 + (yelm(2) - yelm(8))**2 + (zelm(2) - zelm(8))**2)
   dist3 = sqrt((xelm(3) - xelm(5))**2 + (yelm(3) - yelm(5))**2 + (zelm(3) - zelm(5))**2)
   dist4 = sqrt((xelm(4) - xelm(6))**2 + (yelm(4) - yelm(6))**2 + (zelm(4) - zelm(6))**2)
   diagonal_aspect_ratio = max(dist1,dist2,dist3,dist4) / min(dist1,dist2,dist3,dist4)

  end subroutine create_mesh_quality_data_3D

!
!=====================================================================
!

! create mesh quality data for a given 2D spectral element

  subroutine create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD,VP_MAX,delta_t, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,stability,distmin,distmax)

  implicit none

  include "constants.h"

  integer :: icorner,ispec,NSPEC,NPOIN,NGNOD,i

  double precision, dimension(NPOIN) :: x,y,z

  integer, dimension(NGNOD,NSPEC) :: ibool

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  double precision vectorA_x,vectorA_y,vectorA_z
  double precision vectorB_x,vectorB_y,vectorB_z
  double precision norm_A,norm_B,angle_vectors
  double precision distmin,distmax,dist,dist1,dist2
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio

! for stability
  double precision :: stability,VP_MAX,delta_t

! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLL_MAX_STABILITY = 15
  double precision, dimension(NGLL_MAX_STABILITY) :: percent_GLL

! topology of faces of cube for skewness
! only one face in 2D
  integer faces_topo(6)

! store the corners of this element for the skewness routine
  do i = 1,NGNOD
    xelm(i) = x(ibool(i,ispec))
    yelm(i) = y(ibool(i,ispec))
    zelm(i) = z(ibool(i,ispec))
  enddo

! define percentage of smallest distance between GLL points for NGLL points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

! check that the degree is not above the threshold for list of percentages
  if(NGLLX > NGLL_MAX_STABILITY) stop 'degree too high to compute stability value'

! define topology of faces of cube for skewness

! only one face in 2D
  faces_topo(1) = 1
  faces_topo(2) = 2
  faces_topo(3) = 3
  faces_topo(4) = 4

! define wraparound for angles for skewness calculation
  faces_topo(5) = faces_topo(1)
  faces_topo(6) = faces_topo(2)

! compute equiangle skewness (as defined in Fluent/Gambit manual)
! and compute edge aspect ratio using the corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     equiangle_skewness = - HUGEVAL

     do icorner = 1,4

! first vector of angle
       vectorA_x = xelm(faces_topo(icorner)) - xelm(faces_topo(icorner+1))
       vectorA_y = yelm(faces_topo(icorner)) - yelm(faces_topo(icorner+1))
       vectorA_z = zelm(faces_topo(icorner)) - zelm(faces_topo(icorner+1))

! second vector of angle
       vectorB_x = xelm(faces_topo(icorner+2)) - xelm(faces_topo(icorner+1))
       vectorB_y = yelm(faces_topo(icorner+2)) - yelm(faces_topo(icorner+1))
       vectorB_z = zelm(faces_topo(icorner+2)) - zelm(faces_topo(icorner+1))

! norm of vectors A and B
       norm_A = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
       norm_B = sqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

! angle formed by the two vectors
       angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

! compute equiangle skewness
       equiangle_skewness = max(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

! compute min and max size of an edge
       dist = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)

       distmin = min(distmin,dist)
       distmax = max(distmax,dist)

     enddo

! compute edge aspect ratio
   edge_aspect_ratio = distmax / distmin

   stability = delta_t * VP_MAX / (distmin * percent_GLL(NGLLX))

! compute diagonal aspect ratio
   dist1 = sqrt((xelm(1) - xelm(3))**2 + (yelm(1) - yelm(3))**2 + (zelm(1) - zelm(3))**2)
   dist2 = sqrt((xelm(2) - xelm(4))**2 + (yelm(2) - yelm(4))**2 + (zelm(2) - zelm(4))**2)
   diagonal_aspect_ratio = max(dist1,dist2) / min(dist1,dist2)

  end subroutine create_mesh_quality_data_2D

