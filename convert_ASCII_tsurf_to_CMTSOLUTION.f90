
! convert ASCII data based on tsurf to set of CMTSOLUTION files

  program convert_tsurf_to_CMTSOLUTION

  implicit none

  include "constants.h"

! constants provided by the user

! average slip length at each center (in meters)
!! DK DK mean value for SAF 1857 from Sieh (1978) is 4.90 m
  double precision, parameter :: AVERAGE_SLIP_LENGTH = 4.90d0

! fictitious constant value of mu = 1
! varies along the fault plane, therefore will be added later in the solver
  double precision, parameter :: mu = 1.d0

! rupture velocity in m/s
  double precision, parameter :: rupture_velocity = 2800.d0

! depth above which we suppress source elements that reach the surface
  double precision, parameter :: DEPTH_REMOVED_TOPO = 20.d0

! length of normal vectors for DX display
  double precision, parameter :: length_normal_display_DX = +20000.d0

 integer npoin    ! number of VRTX in GoCad file
 integer nspec    ! number of TRGL in GoCad file

 integer, dimension(:), allocatable :: iglob1,iglob2,iglob3
 integer, dimension(:), allocatable :: iglob1_copy,iglob2_copy,iglob3_copy
 double precision, dimension(:), allocatable :: x,y,z

 double precision horiz_dist_fault,time_shift,time_shift_min,time_shift_max

 integer ipoin,ispec,isource,NSOURCES,isourceshiftmin,isource_current
 integer iglob1_store,iglob2_store,iglob3_store,iglob_dummy

 logical exclude_source

 double precision x1,y1,z1,x2,y2,z2,x3,y3,z3,horizdistval,TOLERANCE
 double precision area_current,area_min,area_max,area_sum
 double precision nx,ny,nz,norm
 double precision x_center_triangle,y_center_triangle,z_center_triangle
 double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz,long,lat
 double precision moment_tensor(6)

! fault edges
 double precision x_begin,x_end,y_begin,y_end
 double precision xmin,xmax,ymin,ymax

! slip vector Cartesian components
 double precision ex,ey,ez,real_slip_length

 double precision, external :: area_triangle

! for parameter file
  integer NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,UTM_PROJECTION_ZONE,NSTEP
  double precision DEPTH_BLOCK_KM,LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,DT

! use integer array to store topography values
  integer icornerlat,icornerlong
  double precision elevation,long_corner,lat_corner,ratio_xi,ratio_eta
  integer itopo_bathy_basin(NX_TOPO,NY_TOPO)

! slip distribution from Sieh (1978) for San Andreas 1857
  integer, parameter :: NPOIN_SIEH_1978 = 2637
  double precision, dimension(NPOIN_SIEH_1978) :: xslip_sieh_1978,yslip_sieh_1978
  integer index_min_dist
  double precision dist_sieh,min_dist

  integer SIGN_NORMAL,event_number
  character(len=40) tsurf_file

! first 27 characters of each line in parameter file are a comment
  character(len=27) junk

  integer i

! read part of parameter file
  open(unit=33,file='DATA/Par_file',status='old')

! ignore header
  do i=1,11
    read(33,*)
  enddo

  read(33,1) junk,NEX_XI
  read(33,1) junk,NEX_ETA
  read(33,*)
  read(33,*)
  read(33,1) junk,NPROC_XI
  read(33,1) junk,NPROC_ETA
  read(33,*)
  read(33,*)
  read(33,2) junk,DT
  read(33,*)
  read(33,*)
  read(33,1) junk,NSTEP
  read(33,*)
  read(33,*)
  read(33,2) junk,LAT_MIN
  read(33,2) junk,LAT_MAX
  read(33,2) junk,LONG_MIN
  read(33,2) junk,LONG_MAX
  read(33,2) junk,DEPTH_BLOCK_KM
  read(33,1) junk,UTM_PROJECTION_ZONE

  close(33)

! formats
 1 format(a,i8)
 2 format(a,f12.5)

! compute tolerance in degrees to exclude sources located near model edges
! exclude typically 5 % of size of model
 TOLERANCE = dabs(LAT_MAX-LAT_MIN)*5.d0/100.d0

! ask user to select event
 print *,'select event:'
 print *
 print *,'  1 = 1857'
 print *,'  2 = Whittier'
 print *,'  3 = Southern San Andreas'
 print *

 read(5,*) event_number

 if(event_number == 1) then
  SIGN_NORMAL = +1
  tsurf_file = 'ASCII_1857_rupture_remeshed.dat'

 else if(event_number == 2) then
   SIGN_NORMAL = -1
   tsurf_file = 'ASCII_whittier_remeshed.dat'

 else if(event_number == 3) then
   SIGN_NORMAL = -1
   tsurf_file = 'ASCII_southern_san_andreas_remeshed.dat'

 else
   stop 'incorrect event selected'
 endif

 print *
 print *,'tsurf file name is ',tsurf_file
 print *

 open(unit=10,file=tsurf_file,status='old')
 read(10,*) npoin
 read(10,*) nspec

 print *
 print *,'number of points in data file is ',npoin
 print *,'number of elements in data file is ',nspec
 print *

 allocate(x(npoin))
 allocate(y(npoin))
 allocate(z(npoin))

 do ipoin = 1,npoin
   read(10,*) iglob_dummy,x(ipoin),y(ipoin),z(ipoin)
 enddo

 allocate(iglob1(nspec))
 allocate(iglob2(nspec))
 allocate(iglob3(nspec))

 allocate(iglob1_copy(nspec))
 allocate(iglob2_copy(nspec))
 allocate(iglob3_copy(nspec))

 do ispec = 1,nspec
   read(10,*) iglob1(ispec),iglob2(ispec),iglob3(ispec)
 enddo

 close(10)

!----

! read basin topography and bathymetry file
  write(*,*)
  write(*,*) 'reading topography file'
  call read_basin_topo_bathy_file(itopo_bathy_basin)
  write(*,*) 'regional topography file read ranges in m from ', &
        minval(itopo_bathy_basin),' to ',maxval(itopo_bathy_basin)
  write(*,*)

!----

! remove sources that are outside or above the model

  isource = 0

  do isource_current = 1,nspec

! get coordinates of corners of this triangle
  x1 = x(iglob1(isource_current))
  y1 = y(iglob1(isource_current))
  z1 = z(iglob1(isource_current))

  x2 = x(iglob2(isource_current))
  y2 = y(iglob2(isource_current))
  z2 = z(iglob2(isource_current))

  x3 = x(iglob3(isource_current))
  y3 = y(iglob3(isource_current))
  z3 = z(iglob3(isource_current))

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

! convert location of source
  call utm_geo(long,lat,x_center_triangle,y_center_triangle,UTM_PROJECTION_ZONE)

! determine if we need to exclude this source
  exclude_source = .false.

! exclude if horizontal coordinates are outside the block
  if(lat <= LAT_MIN + TOLERANCE .or. lat >= LAT_MAX - TOLERANCE .or. &
     long <= LONG_MIN + TOLERANCE .or. long >= LONG_MAX - TOLERANCE .or. &
     z_center_triangle <= - dabs(DEPTH_BLOCK_KM*1000.d0)) &
             exclude_source = .true.

! exclude if source is above topography

! get coordinate of corner in bathy/topo model
    icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
    icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! avoid edge effects and extend with identical point if outside model
    if(icornerlong < 1) icornerlong = 1
    if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
    if(icornerlat < 1) icornerlat = 1
    if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! compute coordinates of corner
    long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
    lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! compute ratio for interpolation
    ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
    ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

! avoid edge effects
    if(ratio_xi < 0.) ratio_xi = 0.
    if(ratio_xi > 1.) ratio_xi = 1.
    if(ratio_eta < 0.) ratio_eta = 0.
    if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
    elevation = &
      itopo_bathy_basin(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
      itopo_bathy_basin(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
      itopo_bathy_basin(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

  if(z_center_triangle > elevation - DEPTH_REMOVED_TOPO) &
    exclude_source = .true.

! store current point if source is kept
  if(.not. exclude_source) then
    isource = isource + 1
    iglob1_copy(isource) = iglob1(isource_current)
    iglob2_copy(isource) = iglob2(isource_current)
    iglob3_copy(isource) = iglob3(isource_current)
  endif

  enddo

! store total number of actual sources
  NSOURCES = isource

  print *
  print *,'keeping ',NSOURCES,' sources inside the model out of ',nspec
  print *,'excluding ',nspec-NSOURCES,' sources outside or above the model (', &
             sngl(100.d0*dble(nspec-NSOURCES)/dble(nspec)),' %)'
  print *

!----

! for min and max size of triangle
  area_min = + 10000000000.d0
  area_max = - 10000000000.d0
  area_sum = 0.d0

  time_shift_min = + 10000000000.d0
  time_shift_max = - 10000000000.d0
  isourceshiftmin = -1

  xmin = + 10000000000.d0
  xmax = - 10000000000.d0

  ymin = + 10000000000.d0
  ymax = - 10000000000.d0

! compute min and max of fault location

  do isource = 1,NSOURCES

! get coordinates of corners of this triangle

  x1 = x(iglob1_copy(isource))
  y1 = y(iglob1_copy(isource))
  z1 = z(iglob1_copy(isource))

  x2 = x(iglob2_copy(isource))
  y2 = y(iglob2_copy(isource))
  z2 = z(iglob2_copy(isource))

  x3 = x(iglob3_copy(isource))
  y3 = y(iglob3_copy(isource))
  z3 = z(iglob3_copy(isource))

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

  xmin = dmin1(xmin,x_center_triangle)
  xmax = dmax1(xmax,x_center_triangle)

  ymin = dmin1(ymin,y_center_triangle)
  ymax = dmax1(ymax,y_center_triangle)

  enddo

! fault edges
! rupture propagates from North-West to South-East
  x_begin = xmin
  x_end = xmax
  y_begin = ymax
  y_end = ymin

! length of the fault
  horiz_dist_fault = dsqrt((x_end-x_begin)**2 + (y_end-y_begin)**2)
  print *
  print *,'horizontal distance between XYmin and XYmax (km) = ',horiz_dist_fault/1000.d0
  print *

! compute min and max of time shift

  do isource = 1,NSOURCES

! get coordinates of corners of this triangle

  x1 = x(iglob1_copy(isource))
  y1 = y(iglob1_copy(isource))
  z1 = z(iglob1_copy(isource))

  x2 = x(iglob2_copy(isource))
  y2 = y(iglob2_copy(isource))
  z2 = z(iglob2_copy(isource))

  x3 = x(iglob3_copy(isource))
  y3 = y(iglob3_copy(isource))
  z3 = z(iglob3_copy(isource))

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

  time_shift = dsqrt((x_center_triangle-x_begin)**2 + (y_center_triangle-y_begin)**2) / rupture_velocity
! store source with minimum time shift for reference
  if(time_shift < time_shift_min) isourceshiftmin = isource
  time_shift_min = dmin1(time_shift_min,time_shift)
  time_shift_max = dmax1(time_shift_max,time_shift)

  enddo

!----

  print *
  print *,'minimum time shift detected for source ',isourceshiftmin
  print *,'automatically swapping it with source #1'
  print *

! making sure that minimum time shift is set for source #1
  if(isourceshiftmin /= 1) then
    iglob1_store = iglob1_copy(1)
    iglob2_store = iglob2_copy(1)
    iglob3_store = iglob3_copy(1)

    iglob1_copy(1) = iglob1_copy(isourceshiftmin)
    iglob2_copy(1) = iglob2_copy(isourceshiftmin)
    iglob3_copy(1) = iglob3_copy(isourceshiftmin)

    iglob1_copy(isourceshiftmin) = iglob1_store
    iglob2_copy(isourceshiftmin) = iglob2_store
    iglob3_copy(isourceshiftmin) = iglob3_store
  endif

!----

! read slip distribution from Sieh (1978) if San Andreas 1857
! convert horizontal distance from km to m
  if(event_number == 1) then
    print *,'reading slip ditribution from Sieh (1978)'
    open(unit=11,file='DATA/slip_sieh_1857_extracted.txt',status='old')
    do ipoin = 1,NPOIN_SIEH_1978
      read(11,*) xslip_sieh_1978(ipoin),yslip_sieh_1978(ipoin)
      xslip_sieh_1978(ipoin) = xslip_sieh_1978(ipoin) * 1000.d0
    enddo
    close(11)
  endif


!----

  time_shift_min = + 10000000000.d0
  time_shift_max = - 10000000000.d0
  isourceshiftmin = -1

! write result to CMTSOLUTION file
  open(unit=11,file='DATA/CMTSOLUTION',status='unknown')

!! DK DK this for seismic moment file
  open(unit=14,file='DATA/seismic_moment_calculation.dat',status='unknown')

  do isource = 1,NSOURCES

! get coordinates of corners of this triangle

  x1 = x(iglob1_copy(isource))
  y1 = y(iglob1_copy(isource))
  z1 = z(iglob1_copy(isource))

  x2 = x(iglob2_copy(isource))
  y2 = y(iglob2_copy(isource))
  z2 = z(iglob2_copy(isource))

  x3 = x(iglob3_copy(isource))
  y3 = y(iglob3_copy(isource))
  z3 = z(iglob3_copy(isource))

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

! compute area of triangle
  area_current = area_triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)

! compute normal vector to triangle
  nx = (y3-y2)*(z2-z1) - (z3-z2)*(y2-y1)
  ny = (z3-z2)*(x2-x1) - (x3-x2)*(z2-z1)
  nz = (x3-x2)*(y2-y1) - (y3-y2)*(x2-x1)

! normalize normal vector
  norm = dsqrt(nx**2 + ny**2 + nz**2)
!! DK DK invert normal or not
  nx = SIGN_NORMAL * nx / norm
  ny = SIGN_NORMAL * ny / norm
  nz = SIGN_NORMAL * nz / norm
!! DK DK fix problem of different normals for Southern San Andreas
  if(event_number == 3) then
    if(ny < 0.) then
      nx = - nx
      ny = - ny
      nz = - nz
    endif
  endif

! compute min, max and total area
  area_min = dmin1(area_min,area_current)
  area_max = dmax1(area_max,area_current)
  area_sum = area_sum + area_current

! compute slip vector: Whittier and San Andreas are right-lateral strike-slip
! i.e., rake = 180 degrees

! take a vector in the plane of the triangle, by using two corners
  ex = x2-x1
  ey = y2-y1
  ez = z2-z1
  horizdistval = dabs(ex)

! avoid null values when two points are aligned on the same vertical
! by choosing the edge with the longest projection along X
  if(dabs(x3-x1) > horizdistval) then
    ex = x3-x1
    ey = y3-y1
    ez = z3-z1
    horizdistval = dabs(ex)
  endif

  if(dabs(x3-x2) > horizdistval) then
    ex = x3-x2
    ey = y3-y2
    ez = z3-z2
    horizdistval = dabs(ex)
  endif

! because of fault orientation, slip should always be from East to West
!! DK DK UGLY check this, could be the other way around
  if(ex > 0.) then
    ex = - ex
    ey = - ey
    ez = - ez
  endif

! test of slip vector orientation for right-lateral strike-slip
!! DK DK UGLY check this, could be the other way around
!! DK DK UGLY fix this tomorrow  if(ex > 0. .or. ey < 0.) stop 'wrong orientation of slip vector'

! make sure slip is horizontal (for strike-slip)
  ez = 0.

! determine slip length from Sieh (1978) if SAF 1857, use constant otherwise
  if(event_number /= 1) then

    real_slip_length = AVERAGE_SLIP_LENGTH

  else

    min_dist = + 1000000000.d0
    index_min_dist = -1
! distance to origin of segment that ruptured
    horiz_dist_fault = dsqrt((x_center_triangle-x_begin)**2 + (y_center_triangle-y_begin)**2)
! loop on all the points to find the closest
    do ipoin = 1,NPOIN_SIEH_1978
! distance to current point in Sieh (1978) scanned curved
      dist_sieh = dabs(xslip_sieh_1978(ipoin) - horiz_dist_fault)
      if(dist_sieh < min_dist) then
        min_dist = dist_sieh
        index_min_dist = ipoin
      endif
    enddo
    real_slip_length = yslip_sieh_1978(index_min_dist)

  endif

! normalize slip vector and convert to real slip in meters
  norm = dsqrt(ex**2 + ey**2 + ez**2)
  ex = real_slip_length * ex / norm
  ey = real_slip_length * ey / norm
  ez = real_slip_length * ez / norm

! compute moment tensor
  Mxx = mu * area_current * (ex*nx + nx*ex)
  Myy = mu * area_current * (ey*ny + ny*ey)
  Mzz = mu * area_current * (ez*nz + nz*ez)

  Mxy = mu * area_current * (ex*ny + nx*ey)
  Mxz = mu * area_current * (ex*nz + nx*ez)
  Myz = mu * area_current * (ey*nz + ny*ez)

! fictitious info about event
  write(11,"(a)") 'PDE 2003  7  7 23 59 17.78  34.0745 -118.3792   6.4 4.2 4.2 FICTITIOUS'
  write(11,"(a)") 'event name:     9903873'

! time shift must be zero for first source, but then can be different from zero
! mimic rupture with constant rupture velocity
  time_shift = dsqrt((x_center_triangle-x_begin)**2 + (y_center_triangle-y_begin)**2) / rupture_velocity

! store source with minimum time shift for reference
  if(time_shift < time_shift_min) isourceshiftmin = isource
  time_shift_min = dmin1(time_shift_min,time_shift)
  time_shift_max = dmax1(time_shift_max,time_shift)

  if(isource == 1) then
    write(11,"('time shift:  0')")
  else
    write(11,"('time shift:  ',e)") time_shift
  endif

!! DK DK use reasonable constant value here
  write(11,"(a)") 'half duration:   2.0000'

! write location of source
  call utm_geo(long,lat,x_center_triangle,y_center_triangle,UTM_PROJECTION_ZONE)
  write(11,"('latitude: ',f)") lat
  write(11,"('longitude: ',f)") long

! write depth in km
  write(11,"('depth: ',f)") dabs(z_center_triangle)/1000.d0

! write the moment tensor
  moment_tensor(1) = + Mzz
  moment_tensor(2) = + Myy
  moment_tensor(3) = + Mxx
  moment_tensor(4) = - Myz
  moment_tensor(5) = + Mxz
  moment_tensor(6) = - Mxy

! scale the moment tensor
! CMTSOLUTION file values are in dyne.cm
! 1 dyne is 1 gram * 1 cm / (1 second)^2
! 1 Newton is 1 kg * 1 m / (1 second)^2
! thus 1 Newton = 100,000 dynes
! therefore 1 dyne.cm = 1e-7 Newton.m

  write(11,"('Mrr:  ',e)") moment_tensor(1) * 1.d7
  write(11,"('Mtt:  ',e)") moment_tensor(2) * 1.d7
  write(11,"('Mpp:  ',e)") moment_tensor(3) * 1.d7
  write(11,"('Mrt:  ',e)") moment_tensor(4) * 1.d7
  write(11,"('Mrp:  ',e)") moment_tensor(5) * 1.d7
  write(11,"('Mtp:  ',e)") moment_tensor(6) * 1.d7

!! DK DK write seismic moment for computation of total magnitude in solver
!! DK DK again, mu will be added in solver, because it depends on location
  write(14,*) mu * area_current * real_slip_length * 1.d7

  enddo  ! end of loop on all the sources

  close(11)

!! DK DK for file with seismic moment
  close(14)

  print *
  print *,'area of smallest triangle (m^2) = ',area_min
  print *,'area of largest triangle (m^2) = ',area_max
  print *,'ratio largest / smallest triangular surface = ',area_max/area_min
  print *,'total area of fault surface (m^2) = ',area_sum
  print *,'total area of fault surface (km^2) = ',area_sum/1.d6
  print *,'mean area of triangle (m^2) = ',area_sum/dble(NSOURCES)

  print *
  print *,'time shift min max (s) = ',time_shift_min,time_shift_max
  print *,'minimum time shift was detected for new source ',isourceshiftmin
  print *,'and was automatically set to exactly zero'
  if(isourceshiftmin /= 1) stop 'minimum time shift should be for first source'

  print *
  print *,'You need to set NSOURCES = ',NSOURCES,' in DATA/Par_file'
  print *

!====================================

! write DX file with normals, to check orientation

! write result to DX file
  open(unit=11,file='DX_normals.dx',status='unknown')

! write points
   write(11,*) 'object 1 class array type float rank 1 shape 3 items ',2*nspec,' data follows'

   do ispec = 1,nspec

  x1 = x(iglob1(ispec))
  y1 = y(iglob1(ispec))
  z1 = z(iglob1(ispec))

  x2 = x(iglob2(ispec))
  y2 = y(iglob2(ispec))
  z2 = z(iglob2(ispec))

  x3 = x(iglob3(ispec))
  y3 = y(iglob3(ispec))
  z3 = z(iglob3(ispec))

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

       write(11,*) sngl(x_center_triangle),sngl(y_center_triangle),sngl(z_center_triangle)

   enddo

! get coordinates of corners of this triangle

   do ispec = 1,nspec

  x1 = x(iglob1(ispec))
  y1 = y(iglob1(ispec))
  z1 = z(iglob1(ispec))

  x2 = x(iglob2(ispec))
  y2 = y(iglob2(ispec))
  z2 = z(iglob2(ispec))

  x3 = x(iglob3(ispec))
  y3 = y(iglob3(ispec))
  z3 = z(iglob3(ispec))

! compute area of triangle
  area_current = area_triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)

! compute center of triangle
  x_center_triangle = (x1 + x2 + x3) / 3.d0
  y_center_triangle = (y1 + y2 + y3) / 3.d0
  z_center_triangle = (z1 + z2 + z3) / 3.d0

! compute normal vector to triangle
  nx = (y3-y2)*(z2-z1) - (z3-z2)*(y2-y1)
  ny = (z3-z2)*(x2-x1) - (x3-x2)*(z2-z1)
  nz = (x3-x2)*(y2-y1) - (y3-y2)*(x2-x1)

! normalize normal vector
  norm = dsqrt(nx**2 + ny**2 + nz**2)
!! DK DK invert normal or not
  nx = SIGN_NORMAL * nx / norm
  ny = SIGN_NORMAL * ny / norm
  nz = SIGN_NORMAL * nz / norm
!! DK DK fix problem of different normals for Southern San Andreas
  if(event_number == 3) then
    if(ny < 0.) then
      nx = - nx
      ny = - ny
      nz = - nz
    endif
  endif

! compute min, max and total area
  area_min = dmin1(area_min,area_current)
  area_max = dmax1(area_max,area_current)
  area_sum = area_sum + area_current

  xmin = dmin1(xmin,x_center_triangle)
  xmax = dmax1(xmax,x_center_triangle)

  ymin = dmin1(ymin,y_center_triangle)
  ymax = dmax1(ymax,y_center_triangle)

       write(11,*) sngl(x_center_triangle+nx*length_normal_display_DX), &
                   sngl(y_center_triangle+ny*length_normal_display_DX), &
                   sngl(z_center_triangle+nz*length_normal_display_DX)
   enddo

! write elements
   write(11,*) 'object 2 class array type int rank 1 shape 2 items ',nspec,' data follows'

   do ispec = 1,nspec
       write(11,210) ispec-1,ispec+nspec-1
   enddo

 210     format(i6,1x,i6,1x,i6,1x,i6)

   write(11,*) 'attribute "element type" string "lines"'
   write(11,*) 'attribute "ref" string "positions"'
   write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! write data
   do ispec = 1,nspec

   x1 = x(iglob1(ispec))
   y1 = y(iglob1(ispec))
   z1 = z(iglob1(ispec))

   x2 = x(iglob2(ispec))
   y2 = y(iglob2(ispec))
   z2 = z(iglob2(ispec))

   x3 = x(iglob3(ispec))
   y3 = y(iglob3(ispec))
   z3 = z(iglob3(ispec))

! compute center of triangle
   x_center_triangle = (x1 + x2 + x3) / 3.d0
   y_center_triangle = (y1 + y2 + y3) / 3.d0
   z_center_triangle = (z1 + z2 + z3) / 3.d0

! convert location of source
   call utm_geo(long,lat,x_center_triangle,y_center_triangle,UTM_PROJECTION_ZONE)

! use different color depending on whether the source is inside the basin model
   if(lat <= LAT_MIN + TOLERANCE .or. lat >= LAT_MAX - TOLERANCE .or. &
           long <= LONG_MIN + TOLERANCE .or. long >= LONG_MAX - TOLERANCE .or. &
           z_center_triangle <= - dabs(DEPTH_BLOCK_KM*1000.d0)) then
       write(11,*) '200'
   else
       write(11,*) '0'
   endif

   enddo

   write(11,*) 'attribute "dep" string "connections"'
   write(11,*) 'object "irregular positions irregular connections" class field'
   write(11,*) 'component "positions" value 1'
   write(11,*) 'component "connections" value 2'
   write(11,*) 'component "data" value 3'
   write(11,*) 'end'

  close(11)

  deallocate(iglob1,iglob2,iglob3,x,y,z)
  deallocate(iglob1_copy,iglob2_copy,iglob3_copy)

  end program convert_tsurf_to_CMTSOLUTION

! ---------------

! compute area of a triangle given the coordinates of its corners
  double precision function area_triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3)

  implicit none

  double precision, intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3

  double precision theta,height,length1,length2,length3

! compute the length of the three sides
  length1 = dsqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
  length2 = dsqrt((x3-x1)**2 + (y3-y1)**2 + (z3-z1)**2)
  length3 = dsqrt((x3-x2)**2 + (y3-y2)**2 + (z3-z2)**2)

! compute the area
  theta = dacos((length1**2+length2**2-length3**2)/(2.0d0*length1*length2))
  height = length1*dsin(theta)
  area_triangle = 0.5d0*length2*height

  end function area_triangle

!! DK DK include routine to read topography
 include "read_basin_topo_bathy_file.f90"

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE)

! convert geodetic longitude and latitude to UTM, and back
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  implicit none

!
!-----CAMx v2.03
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Input/Output arguments:
!
!        rlon                  Longitude (deg, negative for West)
!        rlat                  Latitude (deg)
!        rx                    UTM easting (m)
!        ry                    UTM northing (m)
!        UTM_PROJECTION_ZONE  UTM zone
!

  integer UTM_PROJECTION_ZONE
  double precision rx,ry,rlon,rlat

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

  double precision, parameter :: degrad=PI/180., raddeg=180./PI
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr
  double precision f1,f2,f3,f4,rm,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

! save original parameters
  rlon_save = rlon
  rlat_save = rlat
  rx_save = rx
  ry_save = ry

! define parameters of reference ellipsoid
  e2=1.0-(semimin/semimaj)**2.0
  e4=e2*e2
  e6=e2*e4
  ep2=e2/(1.-e2)

  xx = rx
  yy = ry
!
!----- Set Zone parameters
!
  zone = dble(UTM_PROJECTION_ZONE)
  cm = zone*6.0 - 183.0
  cmr = cm*degrad

!
!---- UTM to Lat/Lon conversion
!

  xx = xx - east
  yy = yy - north
  e1 = sqrt(1. - e2)
  e1 = (1. - e1)/(1. + e1)
  rm = yy/scfa
  u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
  u = rm/(semimaj*u)

  f1 = 3.*e1/2. - 27.*e1**3./32.
  f1 = f1*sin(2.*u)
  f2 = 21.*e1**2/16. - 55.*e1**4/32.
  f2 = f2*sin(4.*u)
  f3 = 151.*e1**3./96.
  f3 = f3*sin(6.*u)
  rlat1 = u + f1 + f2 + f3
  dlat1 = rlat1*raddeg
  if (dlat1 >= 90. .or. dlat1 <= -90.) then
    dlat1 = dmin1(dlat1,dble(90.) )
    dlat1 = dmax1(dlat1,dble(-90.) )
    dlon = cm
  else
    c1 = ep2*cos(rlat1)**2.
    t1 = tan(rlat1)**2.
    f1 = 1. - e2*sin(rlat1)**2.
    rn1 = semimaj/sqrt(f1)
    r1 = semimaj*(1. - e2)/sqrt(f1**3)
    d = xx/(rn1*scfa)

    f1 = rn1*tan(rlat1)/r1
    f2 = d**2/2.
    f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
    f3 = f3*d**2*d**2/24.
    f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
    f4 = f4*(d**2)**3./720.
    rlat = rlat1 - f1*(f2 - f3 + f4)
    dlat = rlat*raddeg

    f1 = 1. + 2.*t1 + c1
    f1 = f1*d**2*d/6.
    f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
    f2 = f2*(d**2)**2*d/120.
    rlon = cmr + (d - f1 + f2)/cos(rlat1)
    dlon = rlon*raddeg
    if (dlon < -180.) dlon = dlon + 360.
    if (dlon > 180.) dlon = dlon - 360.
  endif

  rlon = dlon
  rlat = dlat
  rx = rx_save
  ry = ry_save

  end subroutine utm_geo

