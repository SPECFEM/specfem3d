
 program convert_northridge_CMT

 implicit none

  include "constants.h"

! UTM zone for L.A.
  integer, parameter :: UTM_PROJECTION_ZONE = 11

!! DK DK values from Wald 1996
  double precision, parameter :: strike=122.
  double precision, parameter :: dip=40.

!! DK DK coordinates of fault top center from Wald et al. 1996
 double precision, parameter :: long_top_fault_center = -118.55,lat_top_fault_center = 34.344,z_top_fault_center = -5000.

!! DK DK coordinates of hypocenter from Wald et al. 1996
 double precision, parameter :: long_hypocenter = -118.546,lat_hypocenter = 34.211,z_hypocenter = -17500.

!! DK DK rupture velocity from Wald et al. 1996
 double precision, parameter :: RUPTURE_VELOCITY = 3000.

!! DK DK number of fault patches in Wald et al. 1996
 integer, parameter :: NX_PATCH = 14
 integer, parameter :: NY_PATCH = NX_PATCH

!! DK DK total moment in Wald et al. 1996
 double precision, parameter :: scalar_moment_total = 1.3d26

!!! DK DK Northridge 1994 Harvard catalog CMT
!!Fault plane:  strike=278    dip=42   slip=65
!  strike=130.
!  dip=53.
!  rake=111.
!  scalar_moment = 1.18e+26    ! dyne/cm

 integer islip(NX_PATCH,NY_PATCH)

 double precision rake,anx,any,anz,dx,dy,dz,am0,scalar_moment_patch
 double precision am(3,3)

 integer ierr

 integer ix,iy,isum

  double precision xval,yval,zval,longval,latval
  double precision xi,eta

  double precision x1,y1,z1
  double precision x2,y2,z2
  double precision x3,y3,z3
  double precision x4,y4,z4

  double precision xold1,yold1,zold1
  double precision xold2,yold2,zold2
  double precision xold3,yold3,zold3
  double precision xold4,yold4,zold4

  double precision xveryold1,yveryold1,zveryold1
  double precision xveryold2,yveryold2,zveryold2
  double precision xveryold3,yveryold3,zveryold3
  double precision xveryold4,yveryold4,zveryold4

  double precision thetadip
  double precision thetastrike

  double precision x_top_fault_center,y_top_fault_center
  double precision x_hypocenter,y_hypocenter

  double precision time_shift,timeshift_min,distance_hypo

  integer ix_time_min,iy_time_min

  integer iglob1,iglob2,iglob3,iglob4

  print *
  print *,'creating file DATA/CMTSOLUTION'

!! DK DK convert coordinates of top center of fault
  call utm_geo(long_top_fault_center,lat_top_fault_center,x_top_fault_center,y_top_fault_center,UTM_PROJECTION_ZONE,ILONGLAT2UTM)

!! DK DK convert coordinates of hypocenter
  call utm_geo(long_hypocenter,lat_hypocenter,x_hypocenter,y_hypocenter,UTM_PROJECTION_ZONE,ILONGLAT2UTM)

!! DK DK define fault plane that ruptured once and for all
  xveryold1 = 0.
  xveryold2 = 0.
  xveryold3 = 0.
  xveryold4 = 0.

  yveryold1 = -9000.
  yveryold2 = +9000.
  yveryold3 = +9000.
  yveryold4 = -9000.

  zveryold1 = 0.
  zveryold2 = 0.
  zveryold3 = -24000.
  zveryold4 = -24000.

!! DK DK rotation to implement dip of fault
  thetadip = - (90. - dip) * PI / 180.
  xold1 =  cos(thetadip)*xveryold1 + sin(thetadip)*zveryold1
  zold1 = -sin(thetadip)*xveryold1 + cos(thetadip)*zveryold1
  yold1 = yveryold1

  xold2 =  cos(thetadip)*xveryold2 + sin(thetadip)*zveryold2
  zold2 = -sin(thetadip)*xveryold2 + cos(thetadip)*zveryold2
  yold2 = yveryold2

  xold3 =  cos(thetadip)*xveryold3 + sin(thetadip)*zveryold3
  zold3 = -sin(thetadip)*xveryold3 + cos(thetadip)*zveryold3
  yold3 = yveryold3

  xold4 =  cos(thetadip)*xveryold4 + sin(thetadip)*zveryold4
  zold4 = -sin(thetadip)*xveryold4 + cos(thetadip)*zveryold4
  yold4 = yveryold4

!! DK DK then rotation to implement strike of fault
  thetastrike = + strike * PI / 180.
  x1 =  cos(thetastrike)*xold1 + sin(thetastrike)*yold1   + x_top_fault_center
  y1 = -sin(thetastrike)*xold1 + cos(thetastrike)*yold1   + y_top_fault_center
  z1 = zold1     + z_top_fault_center

  x2 =  cos(thetastrike)*xold2 + sin(thetastrike)*yold2   + x_top_fault_center
  y2 = -sin(thetastrike)*xold2 + cos(thetastrike)*yold2   + y_top_fault_center
  z2 = zold2     + z_top_fault_center

  x3 =  cos(thetastrike)*xold3 + sin(thetastrike)*yold3   + x_top_fault_center
  y3 = -sin(thetastrike)*xold3 + cos(thetastrike)*yold3   + y_top_fault_center
  z3 = zold3     + z_top_fault_center

  x4 =  cos(thetastrike)*xold4 + sin(thetastrike)*yold4   + x_top_fault_center
  y4 = -sin(thetastrike)*xold4 + cos(thetastrike)*yold4   + y_top_fault_center
  z4 = zold4     + z_top_fault_center


! read slip value from Dave Wald's file and compute sum
  isum = 0
  do iy = NY_PATCH,1,-1
  do ix = 1,NX_PATCH
    read(*,*) islip(ix,iy)
    isum = isum + islip(ix,iy)
  enddo
  enddo


!! loop on all the patches to find minimum time shift
  timeshift_min = + 100000000.

! create an OpenDX file to check the slip map
  open(unit=15,file='DX_check_slip_map.dx',status='unknown')

! write points
  write(15,*) 'object 1 class array type float rank 1 shape 3 items ',NX_PATCH*NY_PATCH,' data follows'

  do iy = 1,NY_PATCH
  do ix = 1,NX_PATCH

!! DK DK compute longitude and depth of this point in the fault plane
!! DK DK use bilinear interpolation from the four corners of the fault
 xi = dble(ix-1)/dble(NX_PATCH-1)
 eta = dble(iy-1)/dble(NY_PATCH-1)

 xval = (1.-xi)*(1.-eta)*x4 + xi*(1.-eta)*x3 + xi*eta*x2 + (1.-xi)*eta*x1
 yval = (1.-xi)*(1.-eta)*y4 + xi*(1.-eta)*y3 + xi*eta*y2 + (1.-xi)*eta*y1
 zval = (1.-xi)*(1.-eta)*z4 + xi*(1.-eta)*z3 + xi*eta*z2 + (1.-xi)*eta*z1

! write to OpenDX file
 write(15,*) sngl(xval),sngl(yval),sngl(zval)

!! DK DK convert coordinates of current point back to long/lat
  call utm_geo(longval,latval,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT)

! compute distance from current point to hypocenter
  distance_hypo = dsqrt((xval-x_hypocenter)**2 + (yval-y_hypocenter)**2 + (zval-z_hypocenter)**2)

! time shift is distance to hypocenter divided by rupture velocity
!! DK DK need to manually move minimum value to first source in CMTSOLUTION then
 time_shift = distance_hypo / RUPTURE_VELOCITY

! compute and locate minimum time shift
  if(time_shift < timeshift_min) then
    timeshift_min = time_shift
    ix_time_min = ix
    iy_time_min = iy
  endif

  enddo
  enddo  ! end of loop on all the patches


!! write OpenDX elements
   write(15,*) 'object 2 class array type int rank 1 shape 4 items ',(NY_PATCH-1)*(NX_PATCH-1),' data follows'

   do iy = 1,NY_PATCH-1
   do ix = 1,NX_PATCH-1

     iglob1 = (iy-1)*NX_PATCH + ix
     iglob2 = (iy-1)*NX_PATCH + ix+1
     iglob3 = (iy+1-1)*NX_PATCH + ix+1
     iglob4 = (iy+1-1)*NX_PATCH + ix

! in the case of OpenDX, node numbers start at zero
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
     write(15,210) iglob1-1,iglob4-1,iglob2-1,iglob3-1
   enddo
   enddo

 210 format(i6,1x,i6,1x,i6,1x,i6)

   write(15,*) 'attribute "element type" string "quads"'
   write(15,*) 'attribute "ref" string "positions"'
   write(15,*) 'object 3 class array type float rank 0 items ',NX_PATCH*NY_PATCH,' data follows'

! write result to CMTSOLUTION file
  open(unit=11,file='DATA/CMTSOLUTION',status='unknown')

!! loop on all the patches to output CMTSOLUTION value
  do iy = 1,NY_PATCH
  do ix = 1,NX_PATCH

    rake=101.
    scalar_moment_patch = scalar_moment_total * dble(islip(ix,iy))/dble(isum)

!     compute Cartesian components of outward normal and slip
!     vectors from strike, dip and rake
!     arguments:
!     strike         strike angle in degrees (INPUT)
!     dip            dip angle in degrees (INPUT)
!     rake           rake angle in degrees (INPUT)
!     anx,any,anz    components of fault plane outward normal versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     dx,dy,dz       components of slip versor in the Aki-Richards
!                    Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
 ierr = 0
 call pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
 if(ierr /= 0) stop 'error in pl2nd conversion'


!     compute moment tensor Cartesian components (Harvard CMT convention)
!     from outward normal and slip vectors Cartesian components
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     am0            scalar seismic moment  (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
 ierr = 0
 am0 = scalar_moment_patch
 call nd2ha(anx,any,anz,dx,dy,dz,am0,am,ierr)
 if(ierr /= 0) stop 'error in nd2ha conversion'

!! DK DK compute longitude and depth of this point in the fault plane
!! DK DK use bilinear interpolation from the four corners of the fault
 xi = dble(ix-1)/dble(NX_PATCH-1)
 eta = dble(iy-1)/dble(NY_PATCH-1)

 xval = (1.-xi)*(1.-eta)*x4 + xi*(1.-eta)*x3 + xi*eta*x2 + (1.-xi)*eta*x1
 yval = (1.-xi)*(1.-eta)*y4 + xi*(1.-eta)*y3 + xi*eta*y2 + (1.-xi)*eta*y1
 zval = (1.-xi)*(1.-eta)*z4 + xi*(1.-eta)*z3 + xi*eta*z2 + (1.-xi)*eta*z1

!! DK DK convert coordinates of current point back to long/lat
  call utm_geo(longval,latval,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT)

! compute distance from current point to hypocenter
  distance_hypo = dsqrt((xval-x_hypocenter)**2 + (yval-y_hypocenter)**2 + (zval-z_hypocenter)**2)

! time shift is distance to hypocenter divided by rupture velocity
! subtract timeshift_min to make sure origin time is exactly zero
! even if rounded coordinates of hypocenter are slightly off the fault plane
 time_shift = distance_hypo / RUPTURE_VELOCITY - timeshift_min

!! write data value to OpenDX file
!! DK DK to visualize slip map
  write(15,*) islip(ix,iy)
!! DK DK to visualize time shift
!  write(15,*) time_shift

! fictitious info about event
  write(11,"(a)") 'PDE 2003  7  7 23 59 17.78  34.0745 -118.3792   6.4 4.2 4.2 FICTITIOUS'
  write(11,"(a)") 'event name:     9903873'

! time shift
  if(ix == ix_time_min .and. iy == iy_time_min) then
    write(11,"('time shift:   0')")
  else
    write(11,"('time shift:  ',e)") time_shift
  endif

  write(11,"(a)") 'half duration:    0.5'
  write(11,"('latitude: ',f)") latval
  write(11,"('longitude: ',f)") longval
  write(11,"('depth: ',f)") dabs(zval/1000.)

  write(11,"('Mrr:  ',e)") am(3,3)
  write(11,"('Mtt:  ',e)") am(1,1)
  write(11,"('Mpp:  ',e)") am(2,2)
  write(11,"('Mrt:  ',e)") am(1,3) !   ,am(3,1)
  write(11,"('Mrp:  ',e)") am(2,3) !   ,am(3,2)
  write(11,"('Mtp:  ',e)") am(1,2) !   ,am(2,1)

  enddo
  enddo  ! end of loop on all the patches

  close(11)

   write(15,*) 'attribute "dep" string "positions"'
   write(15,*) 'object "irregular positions irregular connections" class field'
   write(15,*) 'component "positions" value 1'
   write(15,*) 'component "connections" value 2'
   write(15,*) 'component "data" value 3'
   write(15,*) 'end'

  close(15)

 print *
 print *,' must set NSOURCES = ',NX_PATCH*NY_PATCH,' in DATA/Par_file'
 print *

 end program convert_northridge_CMT


!*******************************************************************************
!
! Gasperini P. and Vannucci G., FPSPACK: a package of simple Fortran subroutines
! to manage earthquake focal mechanism data
!
!*******************************************************************************
!     BASIC ROUTINES
!*******************************************************************************
  subroutine pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
!
!     compute Cartesian components of outward normal and slip
!     vectors from strike, dip and rake
!
!     usage:
!     call pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
!
!     arguments:
!     strike         strike angle in degrees (INPUT)
!     dip            dip angle in degrees (INPUT)
!     rake           rake angle in degrees (INPUT)
!     anx,any,anz    components of fault plane outward normal versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     dx,dy,dz       components of slip versor in the Aki-Richards
!                    Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision anx,any,anz,dx,dy,dz,strike,dip,rake,wstrik,wdip,wrake
  integer ierr
!
  call fpsset
!
  anx=c0
  any=c0
  anz=c0
  dx=c0
  dy=c0
  dz=c0
  ierr=0
  if(strike<amistr.or.strike>amastr) then
   write(io,'(1x,a,g10.4,a)') 'PL2ND: input STRIKE angle ',strike, &
     ' out of range'
   ierr=1
  endif
  if(dip<amidip.or.dip>amadip) then
   if(dip<amadip.and.dip>-ovrtol) then
      dip=amidip
   else if(dip>amidip.and.dip-amadip<ovrtol) then
      dip=amadip
   else
      write(io,'(1x,a,g10.4,a)') 'PL2ND: input DIP angle ',dip, &
        ' out of range'
      ierr=ierr+2
   endif
  endif
  if(rake<amirak.or.rake>amarak) then
   write(io,'(1x,a,g10.4,a)') 'PL2ND: input RAKE angle ',rake, &
     ' out of range'
   ierr=ierr+4
  endif
  if(ierr/=0) return
  wstrik=strike*dtor
  wdip=dip*dtor
  wrake=rake*dtor
!
  anx=-sin(wdip)*sin(wstrik)
  any=sin(wdip)*cos(wstrik)
  anz=-cos(wdip)
  dx=cos(wrake)*cos(wstrik)+cos(wdip)*sin(wrake)*sin(wstrik)
  dy=cos(wrake)*sin(wstrik)-cos(wdip)*sin(wrake)*cos(wstrik)
  dz=-sin(wdip)*sin(wrake)
  return
  end
!*******************************************************************************
  subroutine nd2pl(wanx,wany,wanz,wdx,wdy,wdz,phi,delta,alam,dipdir &
  ,ierr)
!
!     compute strike, dip, rake and dip directions from Cartesian
!     components of the outward normal and slip vectors
!
!     usage:
!     call nd2pl(anx,any,anz,dx,dy,dz,strike,dip,rake,dipdir,ierr)
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     strike         strike angle in degrees (OUTPUT)
!     dip            dip angle in degrees (OUTPUT)
!     rake           rake angle in degrees (OUTPUT)
!     dipdir         dip direction angle in degrees (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision wanx,wany,wanz,wdx,wdy,wdz,ang,anorm,anx,any,anz,dnorm &
  ,dx,dy,dz,wdelta,wphi,walam,phi,delta,alam,dipdir
  integer ierr
!
  call fpsset
!
  ierr=0
  call angle(wanx,wany,wanz,wdx,wdy,wdz,ang)
  if(abs(ang-c90)>orttol) then
   write(io,'(1x,a,g15.7,a)') 'ND2PL: input vectors not ' &
     //'perpendicular, angle=',ang
   ierr=1
  endif
  call norm(wanx,wany,wanz,anorm,anx,any,anz)
  call norm(wdx,wdy,wdz,dnorm,dx,dy,dz)
  if(anz>c0) then
   call invert(anx,any,anz)
   call invert(dx,dy,dz)
  endif
!
  if(anz==-c1) then
   wdelta=c0
   wphi=c0
   walam=atan2(-dy,dx)
  else
   wdelta=acos(-anz)
   wphi=atan2(-anx,any)
   walam=atan2(-dz/sin(wdelta),dx*cos(wphi)+dy*sin(wphi))
  endif
  phi=wphi/dtor
  delta=wdelta/dtor
  alam=walam/dtor
  phi=mod(phi+c360,c360)
  dipdir=phi+c90
  dipdir=mod(dipdir+c360,c360)
  return
  end
!*******************************************************************************
  subroutine ax2ca(trend,plunge,ax,ay,az,ierr)
!
!     compute cartesian components from trend and plunge
!
!     usage:
!     call ax2ca(trend,plunge,ax,ay,az,ierr)
!
!     arguments:
!     trend          clockwise angle from North in degrees (INPUT)
!     plunge         inclination angle in degrees (INPUT)
!     ax,ay,az       components of the axis direction downward versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input TREND angle out of range
!     2              input PLUNGE angle out of range
!     3              1+2
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision ax,ay,az,trend,plunge
  integer ierr
!
  call fpsset
!
  ax=c0
  ay=c0
  az=c0
  ierr=0
  if(trend<amitre.or.trend>amatre) then
   write(io,'(1x,a,g10.4,a)') 'AX2CA: input TREND angle ',trend, &
     ' out of range'
   ierr=1
  endif
  if(plunge<amiplu.or.plunge>amaplu) then
   if(plunge<amaplu.and.plunge>-ovrtol) then
      plunge=amiplu
   else if(plunge>amiplu.and.plunge-amaplu<ovrtol) then
      plunge=amaplu
   else
      write(io,'(1x,a,g10.4,a)') 'AX2CA: input PLUNGE angle ', &
        plunge,' out of range'
      ierr=ierr+2
   endif
  endif
  if(ierr/=0) return
  ax=cos(plunge*dtor)*cos(trend*dtor)
  ay=cos(plunge*dtor)*sin(trend*dtor)
  az=sin(plunge*dtor)
  return
  end
!*******************************************************************************
  subroutine ca2ax(wax,way,waz,trend,plunge,ierr)
!
!     compute trend and plunge from Cartesian components
!
!     usage:
!     call ca2ax(ax,ay,az,trend,plunge,ierr)
!
!     arguments:
!     ax,ay,az       components of axis direction vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     trend          clockwise angle from North in degrees (OUTPUT)
!     plunge         inclination angle in degrees (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision wax,way,waz,wnorm,ax,ay,az,trend,plunge
  integer ierr
!
  call fpsset
!
  ierr=0
  call norm(wax,way,waz,wnorm,ax,ay,az)
  if(az<c0) call invert(ax,ay,az)
  if(ay/=c0.or.ax/=c0) then
   trend=atan2(ay,ax)/dtor
  else
   trend=c0
  endif
  trend=mod(trend+c360,c360)
  plunge=asin(az)/dtor
  return
  end
!*******************************************************************************
  subroutine pt2nd(wpx,wpy,wpz,wtx,wty,wtz,anx,any,anz,dx,dy,dz &
  ,ierr)
!
!     compute Cartesian component of P and T versors
!     from outward normal and slip vectors
!
!     usage:
!     call pt2nd(px,py,pz,tx,ty,tz,anx,any,anz,dx,dy,dz,ierr)
!
!     arguments:
!     px,py,pz       components of P (maximum dilatation) axis vector
!                    in the Aki-Richards Cartesian coordinate system (INPUT)
!     tx,ty,tz       components of T (maximum tension) axis vector
!                    in the Aki-Richards Cartesian coordinate system (INPUT)
!     anx,any,anz    components of fault plane outward normal versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     dx,dy,dz       components of slip versor in the Aki-Richards
!                    Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision anx,any,anz,dx,dy,dz,wpx,wpy,wpz,wtx,wty,wtz,ang,pnorm &
  ,px,py,pz,tnorm,tx,ty,tz,amn
  integer ierr
!
  call fpsset
!
  anx=c0
  any=c0
  anz=c0
  dx=c0
  dy=c0
  dz=c0
  ierr=0
  call angle(wpx,wpy,wpz,wtx,wty,wtz,ang)
  if(abs(ang-c90)>orttol) then
   write(io,'(1x,a,g15.7,a)') 'PT2ND: input vectors not ' &
     //'perpendicular, angle=',ang
   ierr=1
  endif
  call norm(wpx,wpy,wpz,pnorm,px,py,pz)
  if(pz<c0) call invert(px,py,pz)
  call norm(wtx,wty,wtz,tnorm,tx,ty,tz)
  if(tz<c0) call invert(tx,ty,tz)
  anx=tx+px
  any=ty+py
  anz=tz+pz
  call norm(anx,any,anz,amn,anx,any,anz)
!
  dx=tx-px
  dy=ty-py
  dz=tz-pz
  call norm(dx,dy,dz,amn,dx,dy,dz)
  if(anz>c0) then
   call invert(anx,any,anz)
   call invert(dx,dy,dz)
  endif
  return
  end
!*******************************************************************************
  subroutine nd2pt(wanx,wany,wanz,wdx,wdy,wdz,px,py,pz,tx,ty,tz,bx &
  ,by,bz,ierr)
!
!     compute Cartesian component of P, T and B axes from outward normal
!     and slip vectors
!
!     usage:
!     call nd2pt(anx,any,anz,dx,dy,dz,px,py,pz,tx,ty,tz,bx,by,bz,ierr)
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     px,py,pz       components of downward P (maximum dilatation) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     tx,ty,tz       components of downward T (maximum tension) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     bx,by,bz       components of downward B (neutral) axis versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision wanx,wany,wanz,amn,anx,any,anz,wdx,wdy,wdz,amd,dx,dy,dz &
  ,ang,px,py,pz,tx,ty,tz,bx,by,bz,amp
  integer ierr
!
  call fpsset
!
  ierr=0
  call norm(wanx,wany,wanz,amn,anx,any,anz)
  call norm(wdx,wdy,wdz,amd,dx,dy,dz)
  call angle(anx,any,anz,dx,dy,dz,ang)
  if(abs(ang-c90)>orttol) then
   write(io,'(1x,a,g15.7,a)') 'ND2PT: input vectors not ' &
     //'perpendicular, angle=',ang
   ierr=1
  endif
  px=anx-dx
  py=any-dy
  pz=anz-dz
  call norm(px,py,pz,amp,px,py,pz)
  if(pz<c0) call invert(px,py,pz)
  tx=anx+dx
  ty=any+dy
  tz=anz+dz
  call norm(tx,ty,tz,amp,tx,ty,tz)
  if(tz<c0) call invert(tx,ty,tz)
  call vecpro(px,py,pz,tx,ty,tz,bx,by,bz)
  if(bz<c0) call invert(bx,by,bz)
  return
  end
!*******************************************************************************
  subroutine ar2pt(am,am0,am1,e,am0b,px,py,pz,tx,ty,tz,bx,by,bz &
  ,eta,ierr)
!
!     compute Cartesian components of deformation axes (P, T and B)
!     from moment tensor (Aki & Richards convention)
!
!     usage:
!     call ar2pt(am,am0,am1,e,am0b,px,py,pz,tx,ty,tz,bx,by,bz,eta,ierr)
!
!     arguments:
!     am             seismic moment tensor (3x3 matrix) (INPUT)
!     am0            scalar seismic moment of principal double couple
!                    (OUTPUT)
!     am1            scalar seismic moment of secondary double couple
!                    (OUTPUT)
!     e              isotropic component (trace of the input tensor)
!                    (OUTPUT)
!     am0b           scalar seismic moment of best double couple (OUTPUT)
!     px,py,pz       components of downward P (maximum dilatation) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     tx,ty,tz       components of downward T (maximum tension) axis versor
!                    in the Aki-Richards Cartesian coordinate system (OUTPUT)
!     bx,by,bz       components of downward B (neutral) axis versor in the
!                    Aki-Richards Cartesian coordinate system (OUTPUT)
!     eta            percentage of CLVD remainder (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input tensor not symmetrical: am(1,2).ne.am(2,1)
!     2              input tensor not symmetrical: am(1,3).ne.am(3,1)
!     3              input tensor not symmetrical: am(2,3).ne.am(3,2)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision am0,am1,e,px,py,pz,tx,ty,tz,bx,by,bz,am0b,eta,dum
  integer ierr,i,j,k
  double precision am(3,3),val(3),vec(3,3)
!
  call fpsset
!
  am0=c0
  am1=c0
  e=c0
  am0b=c0
  px=c0
  py=c0
  pz=c0
  tx=c0
  ty=c0
  tz=c0
  bx=c0
  by=c0
  bz=c0
  ierr=0
  if(abs(am(1,2)-am(2,1))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2PT: input tensor not' &
     //' symmetrical, m(1,2)=',am(1,2),' m(2,1)=',am(2,1)
   ierr=1
  endif
  if(abs(am(1,3)-am(3,1))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2PT: input tensor not' &
     //' symmetrical, m(1,3)=',am(1,3),' m(3,1)=',am(3,1)
   ierr=ierr+2
  endif
  if(abs(am(3,2)-am(2,3))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2PT: input tensor not' &
     //' symmetrical, m(2,3)=',am(2,3),' m(3,2)=',am(3,2)
   ierr=ierr+4
  endif
  if(ierr/=0) return
  call avec(am,val,vec)
  e=(val(1)+val(2)+val(3))/c3
!
!     compute deviatoric eigevalues
!
  do i=1,3
   val(i)=val(i)-e
  enddo
!
!     sort deviatoric eigenvalues (with isotropic component removed)
!     and eigenvectors by inverse order of eigenvalue modulus magnitude
!
  do 2 i=1,2
   do 3 j=i+1,3
      if(abs(val(i))<abs(val(j))) then
         dum=val(i)
         val(i)=val(j)
         val(j)=dum
         do 4 k=1,3
            dum=vec(k,i)
            vec(k,i)=vec(k,j)
            vec(k,j)=dum
  4              continue
      endif
  3        continue
  2     continue
  am0=val(1)
  eta=-val(3)/(c2*am0)
  am1=abs(val(3))
  am0b=(abs(val(1))+abs(val(2)))/c2
  if(am0<c0) then
   am0=-am0
   tx=vec(1,2)
   ty=vec(2,2)
   tz=vec(3,2)
   px=vec(1,1)
   py=vec(2,1)
   pz=vec(3,1)
   bx=vec(1,3)
   by=vec(2,3)
   bz=vec(3,3)
  else
   tx=vec(1,1)
   ty=vec(2,1)
   tz=vec(3,1)
   px=vec(1,2)
   py=vec(2,2)
   pz=vec(3,2)
   bx=vec(1,3)
   by=vec(2,3)
   bz=vec(3,3)
  endif
  return
  end
!*******************************************************************************
  subroutine nd2ar(anx,any,anz,dx,dy,dz,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Aki & Richards convention)
!     from outward normal and slip vectors Cartesian components
!
!     usage:
!     call nd2ar(anx,any,anz,dx,dy,dz,am0,am,ierr)
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     am0            scalar seismic moment  (INPPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision anx,any,anz,dx,dy,dz,am0,aam0,ang,anorm,wanx,wany,wanz &
  ,wdx,wdy,wdz
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do  2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  if(am0==c0) then
   aam0=c1
  else
   aam0=am0
  endif
  ierr=0
  call angle(anx,any,anz,dx,dy,dz,ang)
  if(abs(ang-c90)>orttol) then
   write(io,'(1x,a,g15.7,a)') 'ND2AR: input vectors not ' &
     //'perpendicular, angle=',ang
   ierr=1
  endif
  call norm(anx,any,anz,anorm,wanx,wany,wanz)
  call norm(dx,dy,dz,anorm,wdx,wdy,wdz)
  am(1,1)=aam0*c2*wdx*wanx
  am(1,2)=aam0*(wdx*wany+wdy*wanx)
  am(2,1)=am(1,2)
  am(1,3)=aam0*(wdx*wanz+wdz*wanx)
  am(3,1)=am(1,3)
  am(2,2)=aam0*c2*wdy*wany
  am(2,3)=aam0*(wdy*wanz+wdz*wany)
  am(3,2)=am(2,3)
  am(3,3)=aam0*c2*wdz*wanz
  return
  end
!*******************************************************************************
  subroutine ar2ha(am,amo,ierr)
!
!     transforms moment tensor component from Aki Richards to Harvard CMT
!     reference systems and viceversa
!
!     usage:
!     call ar2ha(am,amo,ierr)
!
!     arguments:
!     am             input seismic moment tensor (3x3 matrix) (INPUT)
!     amo            output seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input tensor not symmetrical: am(1,2).ne.am(2,1)
!     2              input tensor not symmetrical: am(1,3).ne.am(3,1)
!     3              input tensor not symmetrical: am(2,3).ne.am(3,2)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  integer ierr,i,j
  double precision am(3,3),amo(3,3)
!
  call fpsset
!
  ierr=0
  if(abs(am(1,2)-am(2,1))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2HA: input tensor not' &
     //' symmetrical, m(1,2)=',am(1,2),' m(2,1)=',am(2,1)
   ierr=1
  endif
  if(abs(am(1,3)-am(3,1))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2HA: input tensor not' &
     //' symmetrical, m(1,3)=',am(1,3),' m(3,1)=',am(3,1)
   ierr=ierr+2
  endif
  if(abs(am(3,2)-am(2,3))>tentol) then
   write(io,'(1x,a,g10.4,a,g10.4)') 'AR2HA: input tensor not' &
     //' symmetrical, m(2,3)=',am(2,3),' m(3,2)=',am(3,2)
   ierr=ierr+4
  endif
  if(ierr/=0) then
   do 1 i=1,3
      do 2 j=1,3
         amo(i,j)=c0
  2           continue
  1        continue
   return
  endif
  amo(1,1)=am(1,1)
  amo(1,2)=-am(1,2)
  amo(1,3)=am(1,3)
  amo(2,1)=-am(2,1)
  amo(2,2)=am(2,2)
  amo(2,3)=-am(2,3)
  amo(3,1)=am(3,1)
  amo(3,2)=-am(3,2)
  amo(3,3)=am(3,3)
  return
  end
!*******************************************************************************
!     COMPOSITE ROUTINES
!*******************************************************************************
  subroutine nd2ha(anx,any,anz,dx,dy,dz,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Harvard CMT convention)
!     from outward normal and slip vectors Cartesian components
!
!     usage:
!     call nd2ha(anx,any,anz,dx,dy,dz,am0,am,ierr)
!
!     arguments:
!     anx,any,anz    components of fault plane outward normal vector in the
!                    Aki-Richards Cartesian coordinate system (INPUT)
!     dx,dy,dz       components of slip vector in the Aki-Richards
!                    Cartesian coordinate system (INPUT)
!     am0            scalar seismic moment  (INPPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input vectors not perpendicular among each other
!     2              internal error
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision anx,any,anz,dx,dy,dz,am0
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do 2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  ierr=0
  call nd2ar(anx,any,anz,dx,dy,dz,am0,am,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'ND2HA: ierr=',ierr
   return
  endif
  call ar2ha(am,am,ierr)
  if(ierr/=0) then
   ierr=2
   write(io,'(1x,a,i3)') 'ND2HA: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pl2pl(strika,dipa,rakea,strikb,dipb,rakeb, &
  dipdib,ierr)
!
!     compute strike, dip and rake of a nodal plane
!     from strike, dip and rake of the other one
!
!
!     usage:
!     call pl2pl(strika,dipa,rakea,strikb,dipb,rakeb,dipdib,ierr)
!
!     arguments:
!     strika         strike angle in degrees of the first nodal plane (INPUT)
!     dipa           dip angle in degrees of the first nodal plane (INPUT)
!     rakea          rake angle in degrees of the first nodal plane (INPUT)
!     strikb         strike angle in degrees of the second nodal plane (OUTPUT)
!     dipb           dip angle in degrees of the second nodal plane (OUTPUT)
!     rakeb          rake angle in degrees of the second nodal plane (OUTPUT)
!     dipdib         dip direction in degrees of the second nodal plane (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8              internal error
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision strika,dipa,rakea,anx,any,anz,dx,dy,dz,strikb,dipb,rakeb, &
  dipdib
  integer ierr
!
  call fpsset
!
  call pl2nd(strika,dipa,rakea,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PL2PL: ierr=',ierr
   return
  endif
  call nd2pl(dx,dy,dz,anx,any,anz,strikb,dipb,rakeb,dipdib,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'PL2PL: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pl2pt(strike,dip,rake,trendp,plungp,trendt,plungt, &
  trendb,plungb,ierr)
!
!     compute trend and plunge of P, T and B axes
!     from strike, dip and rake of a nodal plane
!
!
!     usage:
!     call pl2pt(strike,dip,rake,trendp,plungp,trendt,plungt,trendb,plungb,ierr)
!
!     arguments:
!     strike         strike angle in degrees of the first nodal plane (INPUT)
!     dip            dip angle in degrees of the first nodal plane (INPUT)
!     rake           rake angle in degrees of the first nodal plane (INPUT)
!     trendp         trend of P axis (OUTPUT)
!     plungp         plunge or P axis (OUTPUT)
!     trendt         trend of T axis (OUTPUT)
!     plungt         plunge or T axis (OUTPUT)
!     trendb         trend of B axis (OUTPUT)
!     plungb         plunge or B axis (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8,9,10,11      internal error
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision strike,dip,rake,anx,any,anz,dx,dy,dz,px,py,pz,tx,ty,tz &
  ,bx,by,bz,trendp,plungp,trendt,plungt,trendb,plungb
  integer ierr
!
  call fpsset
!
  call pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PL2PT: ierr=',ierr
   return
  endif
  call nd2pt(dx,dy,dz,anx,any,anz,px,py,pz,tx,ty,tz,bx,by,bz,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'PL2PT: ierr=',ierr
  endif
  call ca2ax(px,py,pz,trendp,plungp,ierr)
  if(ierr/=0) then
   ierr=9
   write(io,'(1x,a,i3)') 'PL2PT: ierr=',ierr
  endif
  call ca2ax(tx,ty,tz,trendt,plungt,ierr)
  if(ierr/=0) then
   ierr=10
   write(io,'(1x,a,i3)') 'PL2PT: ierr=',ierr
  endif
  call ca2ax(bx,by,bz,trendb,plungb,ierr)
  if(ierr/=0) then
   ierr=11
   write(io,'(1x,a,i3)') 'PL2PT: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pt2pl(trendp,plungp,trendt,plungt,strika,dipa,rakea &
  ,dipdia,strikb,dipb,rakeb,dipdib,ierr)
!
!     compute strike dip and rake (and dip direction) of two nodal planes
!     from trend and plung of P and T axes
!
!     usage:
!     call pt2pl(trendp,plungp,trendt,plungt,strika,dipa,rakea
!    1,dipdia,strikb,dipb,rakeb,dipdib,ierr)
!
!     arguments:
!     trendp         trend of P axis in degrees (INPUT)
!     plungp         plunge of P axis in degrees (INPUT)
!     trendt         trend of P axis in degrees (INPUT)
!     plungt         plunge of P axis in degrees (INPUT)
!     strika         strike angle of first nodal plane in degrees (OUTPUT)
!     dipa           dip angle of first nodal plane in degrees (OUTPUT)
!     rakea          rake angle of first nodal plane in degrees (OUTPUT)
!     dipdia         dip direction angle of first nodal plane in degrees (OUTPUT)
!     strikb         strike angle of second nodal plane in degrees (OUTPUT)
!     dipb           dip angle of second nodal plane in degrees (OUTPUT)
!     rakeb          rake angle of second nodal plane in degrees (OUTPUT)
!     dipdib         dip direction angle of second nodal plane in degrees (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input TREND angle of P axis out of range
!     2              input PLUNGE angle P axis out of range
!     3              1+2
!     4              input TREND angle of P axis out of range
!     5              input PLUNGE angle P axis out of range
!     6              4+5
!     8,9,10         internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision trendp,plungp,trendt,plungt,px,py,pz,tx,ty,tz &
  ,anx,any,anz,dx,dy,dz,strika,dipa,rakea,dipdia,strikb,dipb &
  ,rakeb,dipdib
  integer ierr
!
  call fpsset
!
  call ax2ca(trendp,plungp,px,py,pz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PT2PL: ierr=',ierr
   return
  endif
  call ax2ca(trendt,plungt,tx,ty,tz,ierr)
  if(ierr/=0) then
   ierr=ierr+3
   write(io,'(1x,a,i3)') 'PT2PL: ierr=',ierr
   return
  endif
  call pt2nd(px,py,pz,tx,ty,tz,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'PT2PL: ierr=',ierr
   return
  endif
  call nd2pl(anx,any,anz,dx,dy,dz,strika,dipa,rakea,dipdia,ierr)
  if(ierr/=0) then
   ierr=9
   write(io,'(1x,a,i3)') 'PT2PL: ierr=',ierr
   return
  endif
  call nd2pl(dx,dy,dz,anx,any,anz,strikb,dipb,rakeb,dipdib,ierr)
  if(ierr/=0) then
   ierr=10
   write(io,'(1x,a,i3)') 'PT2PL: ierr=',ierr
   return
  endif
  return
  end
!
!*******************************************************************************
  subroutine ar2plp(am,am0,am1,e,am0b,phia,deltaa,alama,slipa, &
  phib,deltab,alamb,slipb,trendp,plungp,trendt,plungt,trendb, &
  plungb,eta,ierr)
!
!     compute planes and axes of principal double couple
!     from moment tensor (Aki & Richards convention)
!
!     usage:
!     call ar2plp(am,am0,am1,e,am0b,phia,deltaa,alama,slipa,
!    1phib,deltab,alamb,slipb,trendp,plungp,trendt,plungt,trendb,
!    2plungb,eta,ierr)
!
!     arguments:
!     am             seismic moment tensor (3x3 matrix) (INPUT)
!     am0            scalar seismic moment of principal double couple
!                    (OUTPUT)
!     am1            scalar seismic moment of secondary double couple
!                    (OUTPUT)
!     e              isotropic component (trace of the input tensor)
!                    (OUTPUT)
!     am0b           scalar seismic moment of best double couple (OUTPUT)
!     phia           strike angle in degrees of the first nodal plane (OUTPUT)
!     deltaa         dip angle in degrees of the first nodal plane (OUTPUT)
!     alama          rake angle in degrees of the first nodal plane (OUTPUT)
!     slipa          dip direction in degrees of the first nodal plane (OUTPUT)
!     phib           strike angle in degrees of the second nodal plane (OUTPUT)
!     deltab         dip angle in degrees of the second nodal plane (OUTPUT)
!     alamb          rake angle in degrees of the second nodal plane (OUTPUT)
!     slipb          dip direction in degrees of the second nodal plane (OUTPUT)
!     trendp         trend angle in degrees of the P axis (OUTPUT)
!     plungp         plunge angle in degrees of the P axis (OUTPUT)
!     trendt         trend angle in degrees of the T axis (OUTPUT)
!     plungt         plunge angle in degrees of the T axis (OUTPUT)
!     trendb         trend angle in degrees of the B axis (OUTPUT)
!     plungb         plunge angle in degrees of the B axis (OUTPUT)
!     eta            percentage of CLVD remainder (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input tensor not symmetrical: am(1,2).ne.am(2,1)
!     2              input tensor not symmetrical: am(1,3).ne.am(3,1)
!     3              input tensor not symmetrical: am(2,3).ne.am(3,2)
!     5,6,7,8,9,10   internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision am0,am1,e,am0b,phia,deltaa,alama,slipa,phib,deltab &
  ,alamb,slipb,trendp,plungp,trendt,plungt,trendb,plungb &
  ,px,py,pz,tx,ty,tz,anx,any,anz,dx,dy,dz,bx,by,bz,eta
  integer ierr
  double precision am(3,3)
!
  call fpsset
!
  am0=c0
  am1=c0
  e=c0
  am0b=c0
  phia=c0
  deltaa=c0
  alama=c0
  slipa=c0
  phib=c0
  deltab=c0
  alamb=c0
  slipb=c0
  trendp=c0
  plungp=c0
  trendt=c0
  plungt=c0
  trendb=c0
  plungb=c0
  ierr=0
  call ar2pt(am,am0,am1,e,am0b,px,py,pz,tx,ty,tz,bx,by,bz,eta,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call ca2ax(px,py,pz,trendp,plungp,ierr)
  if(ierr/=0) then
   ierr=5
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call ca2ax(tx,ty,tz,trendt,plungt,ierr)
  if(ierr/=0) then
   ierr=6
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call ca2ax(bx,by,bz,trendb,plungb,ierr)
  if(ierr/=0) then
   ierr=7
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call pt2nd(px,py,pz,tx,ty,tz,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call nd2pl(anx,any,anz,dx,dy,dz,phia,deltaa,alama,slipa,ierr)
  if(ierr/=0) then
   ierr=9
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  call nd2pl(dx,dy,dz,anx,any,anz,phib,deltab,alamb,slipb,ierr)
  if(ierr/=0) then
   ierr=10
   write(io,'(1x,a,i3)') 'AR2PLP: ierr=',ierr
   return
  endif
  return
  end
!*******************************************************************************
  subroutine ha2plp(am,am0,am1,e,am0b,strika,dipa,rakea,slipa, &
  strikb,dipb,rakeb,slipb,trendp,plungp,trendt,plungt,trendb, &
  plungb,eta,ierr)
!
!     compute planes and axes of principal double couple
!     from moment tensor (Harvard CMT convention)
!
!     usage:
!      subroutine ha2plp(am,am0,am1,e,am0b,strika,dipa,rakea,slipa,
!     1strikb,dipb,rakeb,slipb,trendp,plungp,trendt,plungt,trendb,
!     2plungb,eta,ierr)
!
!     arguments:
!     am             seismic moment tensor (3x3 matrix) (INPUT)
!     am0            scalar seismic moment of principal double couple
!                    (OUTPUT)
!     am1            scalar seismic moment of secondary double couple
!                    (OUTPUT)
!     e              isotropic component (trace of the input tensor)
!                    (OUTPUT)
!     am0b           scalar seismic moment of best double couple (OUTPUT)
!     strika         strike angle in degrees of the first nodal plane (OUTPUT)
!     dipa           dip angle in degrees of the first nodal plane (OUTPUT)
!     rakea          rake angle in degrees of the first nodal plane (OUTPUT)
!     slipa          dip direction in degrees of the first nodal plane (OUTPUT)
!     strikb         strike angle in degrees of the second nodal plane (OUTPUT)
!     dipb           dip angle in degrees of the second nodal plane (OUTPUT)
!     rakeb          rake angle in degrees of the second nodal plane (OUTPUT)
!     slipb          dip direction in degrees of the second nodal plane (OUTPUT)
!     trendp         trend angle in degrees of the P axis (OUTPUT)
!     plungp         plunge angle in degrees of the P axis (OUTPUT)
!     trendt         trend angle in degrees of the T axis (OUTPUT)
!     plungt         plunge angle in degrees of the T axis (OUTPUT)
!     trendb         trend angle in degrees of the B axis (OUTPUT)
!     plungb         plunge angle in degrees of the B axis (OUTPUT)
!     eta            percentage of CLVD remainder (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input tensor not symmetrical: am(1,2).ne.am(2,1)
!     2              input tensor not symmetrical: am(1,3).ne.am(3,1)
!     3              input tensor not symmetrical: am(2,3).ne.am(3,2)
!     5,6,7,8,9,10   internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision am0,am1,e,am0b,strika,dipa,rakea,slipa,strikb,dipb &
  ,rakeb,slipb,trendp,plungp,trendt,plungt,trendb,plungb,eta
  integer ierr
  double precision am(3,3),ama(3,3)
!
  call fpsset
!
  am0=c0
  am1=c0
  e=c0
  am0b=c0
  strika=c0
  dipa=c0
  rakea=c0
  slipa=c0
  strikb=c0
  dipb=c0
  rakeb=c0
  slipb=c0
  trendp=c0
  plungp=c0
  trendt=c0
  plungt=c0
  trendb=c0
  plungb=c0
  ierr=0
  call ar2ha(am,ama,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'HA2PLP: ierr=',ierr
   return
  endif
  call ar2plp(ama,am0,am1,e,am0b,strika,dipa,rakea,slipa, &
  strikb,dipb,rakeb,slipb,trendp,plungp,trendt,plungt,trendb, &
  plungb,eta,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'HA2PLP: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pl2ar(strike,dip,rake,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Aki & Richards convention)
!     from strike, dip and rake
!
!     usage:
!     call pl2ar(strike,dip,rake,am0,am,ierr)
!
!     arguments:
!     strike         strike angle in degrees  (INPUT)
!     dip            dip angle in degrees (INPUT)
!     rake           rake angle in degrees (INPUT)
!     am0            scalar seismic moment  (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8              internal error
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision anx,any,anz,dx,dy,dz,am0,strike,dip,rake
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do 2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  ierr=0
  call pl2nd(strike,dip,rake,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PL2AR: ierr=',ierr
   return
  endif
  call nd2ar(anx,any,anz,dx,dy,dz,am0,am,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'PL2AR: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pl2ha(strike,dip,rake,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Harvard CMT convention)
!     from strike, dip and rake
!
!     usage:
!     call pl2ha(strike,dip,rake,am0,am,ierr)
!
!     arguments:
!     strike         strike angle in degrees  (INPUT)
!     dip            dip angle in degrees (INPUT)
!     rake           rake angle in degrees (INPUT)
!     am0            scalar seismic moment  (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle out of range
!     2              input DIP angle out of range
!     4              input RAKE angle out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8,9            internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision strike,dip,rake,am0
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do 2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  ierr=0
  call pl2ar(strike,dip,rake,am0,am,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PL2HA: ierr=',ierr
   return
  endif
  call ar2ha(am,am,ierr)
  if(ierr/=0) then
   ierr=9
   write(io,'(1x,a,i3)') 'PL2HA: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pt2ar(trendp,plungp,trendt,plungt,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Aki & Richards convention)
!     from P and T axes
!
!     usage:
!     call pt2ar(trendp,plungp,trendt,plungt,am0,am)
!
!     arguments:
!     trendp         trend angle of P axis in degrees (INPUT)
!     plungp         plunge angle of P axis in degrees (INPUT)
!     trendt         trend angle of T axis in degrees (INPUT)
!     plungt         plunge angle of T axis in degrees (INPUT)
!     am0            scalar seismic moment  (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input TREND angle of P axis out of range
!     2              input PLUNGE angle P axis out of range
!     3              1+2
!     4              input TREND angle of P axis out of range
!     5              input PLUNGE angle P axis out of range
!     6              4+5
!     8,9            internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision trendp,plungp,trendt,plungt,am0,px,py,pz,tx,ty,tz &
  ,anx,any,anz,dx,dy,dz
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do 2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  ierr=0
  call ax2ca(trendp,plungp,px,py,pz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PT2AR: ierr=',ierr
   return
  endif
  call ax2ca(trendt,plungt,tx,ty,tz,ierr)
  if(ierr/=0) then
   ierr=ierr+3
   write(io,'(1x,a,i3)') 'PT2AR: ierr=',ierr
   return
  endif
  call pt2nd(px,py,pz,tx,ty,tz,anx,any,anz,dx,dy,dz,ierr)
  if(ierr/=0) then
   ierr=8
   write(io,'(1x,a,i3)') 'PT2AR: ierr=',ierr
   return
  endif
  call nd2ar(anx,any,anz,dx,dy,dz,am0,am,ierr)
  if(ierr/=0) then
   ierr=9
   write(io,'(1x,a,i3)') 'PT2AR: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
  subroutine pt2ha(trendp,plungp,trendt,plungt,am0,am,ierr)
!
!     compute moment tensor Cartesian components (Harvard CMT convention)
!     from P and T axes
!
!     usage:
!     call pt2ha(trendp,plungp,trendt,plungt,am0,am)
!
!     arguments:
!     trendp         trend angle of P axis in degrees (INPUT)
!     plungp         plunge angle of P axis in degrees (INPUT)
!     trendt         trend angle of T axis in degrees (INPUT)
!     plungt         plunge angle of T axis in degrees (INPUT)
!     am0            scalar seismic moment  (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input TREND angle of P axis out of range
!     2              input PLUNGE angle P axis out of range
!     3              1+2
!     4              input TREND angle of P axis out of range
!     5              input PLUNGE angle P axis out of range
!     6              4+5
!     8,9,10         internal errors
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision trendp,plungp,trendt,plungt,am0
  integer ierr,i,j
  double precision am(3,3)
!
  call fpsset
!
  do 1 i=1,3
   do 2 j=1,3
      am(i,j)=c0
  2        continue
  1     continue
  ierr=0
  call pt2ar(trendp,plungp,trendt,plungt,am0,am,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'PT2HA: ierr=',ierr
   return
  endif
  call ar2ha(am,am,ierr)
  if(ierr/=0) then
   ierr=10
   write(io,'(1x,a,i3)') 'PT2HA: ierr=',ierr
  endif
  return
  end
!*******************************************************************************
!     UTILITY ROUTINES
!*******************************************************************************
  subroutine avec(am,eval,evec)
!
!     compute eigenvalues and eigenvectors
!
!     usage:
!     utility routine for internal use only
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision dum
  integer i,j,k
  double precision am(3,3),eval(3),evec(3,3)
!
  call fpsset
!
!! DK DK routine not included  CALL EVCSF (3, AM, 3, EVAL, EVEC, 3)
  stop 'DK DK CALL EVCSF (3, AM, 3, EVAL, EVEC, 3) not included, error'
  do 2 i=1,2
   do 3 j=i+1,3
      if(abs(eval(i))<abs(eval(j))) then
         dum=eval(i)
         eval(i)=eval(j)
         eval(j)=dum
         do 4 k=1,3
            dum=evec(k,i)
            evec(k,i)=evec(k,j)
            evec(k,j)=dum
  4              continue
      endif
  3        continue
  2     continue
  return
  end
!*******************************************************************************
  subroutine angle(wax,way,waz,wbx,wby,wbz,ang)
!
!     compute the angle (in degrees) between two vectors
!
!     usage:
!     call angle(wax,way,waz,wbx,wby,wbz,ang)
!
!     arguments:
!     wax,way,waz    Cartesian component of first vector (INPUT)
!     wbx,wby,wbz    Cartesian component of second vector (INPUT)
!     ang            angle between the two vectors in degrees (OUTPUT)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision wax,way,waz,wbx,wby,wbz,ax,ay,az,bx,by,bz,ang &
  ,anorm,bnorm,prod
!
  call fpsset
!
  call norm(wax,way,waz,anorm,ax,ay,az)
  call norm(wbx,wby,wbz,bnorm,bx,by,bz)
  prod=ax*bx+ay*by+az*bz
  ang=acos(max(-c1,min(c1,prod)))/dtor
  return
  end
!*******************************************************************************
  subroutine angles(strika,dipa,rakea,strikb,dipb,rakeb, &
  anglep,angled,ierr)
!
!     compute the angle (in degrees) between nodal planes and slip vectors
!
!     usage:
!     call angles(strikea,dipa,rakea,strikeb,dipb,rakeb,anglep,angled,ierr)
!
!     arguments:
!     strika         strike angle in degrees of the first nodal plane (INPUT)
!     dipa           dip angle in degrees of the first nodal plane (INPUT)
!     rakea          rake angle in degrees of the first nodal plane (INPUT)
!     strikb         strike angle in degrees of the second nodal plane (INPUT)
!     dipb           dip angle in degrees of the second nodal plane (INPUT)
!     rakeb          rake angle in degrees of the second nodal plane (INPUT)
!     anglep         angle in degrees between planes
!     angled         angle in degrees between slip directions
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input STRIKE angle of first plane out of range
!     2              input DIP angle of first plane out of range
!     4              input RAKE angle of first plane out of range
!     3              1+2
!     5              1+4
!     7              1+2+4
!     8              input STRIKE angle of first plane out of range
!     9              input DIP angle of first plane out of range
!     11             input RAKE angle of first plane out of range
!     10             8+9
!     13             8+11
!     15             8+9+11
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision strika,dipa,rakea,strikb,dipb,rakeb,anglep,angled &
  ,anax,anay,anaz,anbx,anby,anbz,dax,day,daz,dbx,dby,dbz
  integer ierr
!
  call fpsset
!
  call pl2nd(strika,dipa,rakea,anax,anay,anaz,dax,day,daz,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'ANGLES: ierr=',ierr
   return
  endif
  call pl2nd(strikb,dipb,rakeb,anbx,anby,anbz,dbx,dby,dbz,ierr)
  if(ierr/=0) then
   ierr=ierr+8
   write(io,'(1x,a,i3)') 'ANGLES: ierr=',ierr
   return
  endif
  call angle(anax,anay,anaz,anbx,anby,anbz,anglep)
  call angle(dax,day,daz,dbx,dby,dbz,angled)
  return
  end
!*******************************************************************************
  subroutine anglea(trenda,plunga,trendb,plungb,ang,ierr)
!
!     usage:
!     call anglea(trenda,plunga,trendb,plungb,ang,ierr)
!
!     arguments:
!     trenda         clockwise angle from North in degrees of first axis (INPUT)
!     plunga         inclination angle in degrees of first axis (INPUT)
!     trendb         clockwise angle from North in degrees of second axis (INPUT)
!     plungb         inclination angle in degrees of second axis (INPUT)
!     ang            angle in degrees between the axes (OUTPUT)
!     ierr           error indicator (OUTPUT)
!
!     errors:
!     1              input TREND angle of first axis out of range
!     2              input PLUNGE angle of first axis out of range
!     3              1+2
!     5              input TREND angle of first axis out of range
!     6              input PLUNGE angle of first axis out of range
!     7              5+6
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision trenda,plunga,ax,ay,az,trendb,plungb,bx,by,bz,ang
  integer ierr
!
  call fpsset
!
  call ax2ca(trenda,plunga,ax,ay,az,ierr)
  if(ierr/=0) then
   write(io,'(1x,a,i3)') 'ANGLEA: ierr=',ierr
   return
  endif
  call ax2ca(trendb,plungb,bx,by,bz,ierr)
  if(ierr/=0) then
   ierr=ierr+4
   write(io,'(1x,a,i3)') 'ANGLEA: ierr=',ierr
   return
  endif
  call angle(ax,ay,az,bx,by,bz,ang)
  return
  end
!*******************************************************************************
  subroutine norm(wax,way,waz,anorm,ax,ay,az)
!
!     compute euclidean norm and versor components
!
!     usage:
!     call norm(wax,way,waz,anorm,ax,ay,az)
!
!     arguments:
!     wax,way,waz    Cartesian component of input vector (INPUT)
!     anorm          Euclidean norm of input vector (OUTPUT)
!     ax,ay,az       normalized Cartesian component of the vector (OUTPUT)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision wax,way,waz,anorm,ax,ay,az
!
  call fpsset
!
  anorm=sqrt(wax*wax+way*way+waz*waz)
  if(anorm==c0) return
  ax=wax/anorm
  ay=way/anorm
  az=waz/anorm
  return
  end
!*******************************************************************************
  subroutine vecpro(px,py,pz,tx,ty,tz,bx,by,bz)
!
!     compute vector products of two vectors
!
!     usage:
!     call vecpro(px,py,pz,tx,ty,tz,bx,by,bz)
!
!     arguments:
!
!     px,py,pz       Cartesian component of first vector (INPUT)
!     tx,ty,tz       Cartesian component of second vector (INPUT)
!     bx,by,bz       Cartesian component of vector product (OUTUT)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision px,py,pz,tx,ty,tz,bx,by,bz
!
  call fpsset
!
  bx=py*tz-pz*ty
  by=pz*tx-px*tz
  bz=px*ty-py*tx
  return
  end
!*******************************************************************************
  subroutine invert(ax,ay,az)
!
!     invert vector
!
!     usage:
!     utility routine for internal use only
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision ax,ay,az
!
  call fpsset
!
  ax=-ax
  ay=-ay
  az=-az
  return
  end
!*******************************************************************************
  subroutine hatens(amrr,amss,amee,amrs,amre,amse,am)
!
!     build the Harvard CMT tensor from independent components
!
!     usage:
!     call hatens(amrr,amss,amee,amrs,amre,amse,am)
!
!     arguments:
!     amss           Harvard CMT South-South component (INPUT)
!     amse           Harvard CMT South-East component (INPUT)
!     amrs           Harvard CMT Radial-South component (INPUT)
!     amee           Harvard CMT East-East component (INPUT)
!     amre           Harvard CMT Radial-East component (INPUT)
!     amrr           Harvard CMT Radial-Radial component (INPUT)
!     am             seismic moment tensor (3x3 matrix) (OUTPUT)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision amss,amse,amrs,amee,amre,amrr
  double precision am(3,3)
!
  call fpsset
!
  am(1,1)=amss
  am(1,2)=amse
  am(1,3)=amrs
  am(2,1)=amse
  am(2,2)=amee
  am(2,3)=amre
  am(3,1)=amrs
  am(3,2)=amre
  am(3,3)=amrr
  return
  end
!*******************************************************************************
  subroutine tensha(am,amrr,amss,amee,amrs,amre,amse)
!
!     give independent components from Harvard CMT tensor
!
!     usage:
!     call tensha(am,amrr,amss,amee,amrs,amre,amse)
!
!     arguments:
!     am             seismic moment tensor (3x3 matrix) (INPUT)
!     amss           Harvard CMT South-South component (OUTPUT)
!     amse           Harvard CMT South-East component (OUTPUT)
!     amrs           Harvard CMT Radial-South component (OUTPUT)
!     amee           Harvard CMT East-East component (OUTPUT)
!     amre           Harvard CMT Radial-East component (OUTPUT)
!     amrr           Harvard CMT Radial-Radial component (OUTPUT)
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  double precision amss,amse,amrs,amee,amre,amrr
  double precision am(3,3)
!
  call fpsset
!
  amrr=am(3,3)
  amss=am(1,1)
  amee=am(2,2)
  amrs=am(3,1)
  amre=am(3,2)
  amse=am(1,2)
  return
  end
!*******************************************************************************
  subroutine fpsset
!
!     define constants (i.e. input ranges and tolerances) used throughout the
!     package. It is called by every subroutines to setup constants
!
!     usage:
!     call fpsset
!
!     constants in fpscom common block:
!
!     amistr         strike lower limit
!     amastr         strike upper limit
!     amidip         dip lower limit
!     amadip         dip upper limit
!     amirak         rake lower limit
!     amarak         rake upper limit
!     amitre         trend lower limit
!     amatre         trend upper limit
!     amiplu         plunge lower limit
!     amaplu         plunge upper limit
!     orttol         orthogonality tolerance
!     ovrtol         dip overtaking tolerance
!     tentol         moment tensor symmetry tolerance
!     dtor       degree to radians
!     c360       360.
!     c90            90.
!     c0             0.
!     c1             1.
!     c2             2.
!     c3             3.
!     io             error messages file unit
!
       implicit none
!-------------------------------------------------------------------------------
  integer io
  double precision amistr,amastr,amidip,amadip,amirak,amarak,amitre,amatre &
  ,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2,c3
  common /fpscom/amistr,amastr,amidip,amadip,amirak,amarak,amitre &
  ,amatre,amiplu,amaplu,orttol,ovrtol,tentol,dtor,c360,c90,c0,c1,c2 &
  ,c3,io
!-------------------------------------------------------------------------------
  integer ifl
  save ifl
  data ifl/0/
  if(ifl==0) then
   amistr=-360.
   amastr=360.
   amidip=0.
   amadip=90.
   amirak=-360.
   amarak=360.
   amitre=-360.
   amatre=360.
   amiplu=0.
   amaplu=90.
   orttol=2.
   ovrtol=0.001
   tentol=0.0001
   dtor=0.017453292519943296
   c360=360.
   c90=90.
   c0=0.
   c1=1.
   c2=2.
   c3=3.
   io=6
   ifl=1
  endif
  return
  end

!=====================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

  subroutine utm_geo(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway)

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  implicit none

  include "constants.h"

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
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic
!

  integer UTM_PROJECTION_ZONE,iway
  double precision rx,ry,rlon,rlat

  double precision, parameter :: degrad=PI/180., raddeg=180./PI
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

  if(SUPPRESS_UTM_PROJECTION) then
    if (iway == ILONGLAT2UTM) then
      rx = rlon
      ry = rlat
    else
      rlon = rx
      rlat = ry
    endif
    return
  endif

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

  if (iway == IUTM2LONGLAT) then
    xx = rx
    yy = ry
  else
    dlon = rlon
    dlat = rlat
  endif
!
!----- Set Zone parameters
!
  zone = dble(UTM_PROJECTION_ZONE)
  cm = zone*6.0 - 183.0
  cmr = cm*degrad
!
!---- Lat/Lon to UTM conversion
!
  if (iway == ILONGLAT2UTM) then

  rlon = degrad*dlon
  rlat = degrad*dlat

  delam = dlon - cm
  if (delam < -180.) delam = delam + 360.
  if (delam > 180.) delam = delam - 360.
  delam = delam*degrad

  f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
  f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
  f2 = f2*sin(2.*rlat)
  f3 = 15.*e4/256.*45.*e6/1024.
  f3 = f3*sin(4.*rlat)
  f4 = 35.*e6/3072.
  f4 = f4*sin(6.*rlat)
  rm = semimaj*(f1 - f2 + f3 - f4)
  if (dlat == 90. .or. dlat == -90.) then
    xx = 0.
    yy = scfa*rm
  else
    rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
    t = tan(rlat)**2
    c = ep2*cos(rlat)**2
    a = cos(rlat)*delam

    f1 = (1. - t + c)*a**3/6.
    f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
    f2 = f2*a**5/120.
    xx = scfa*rn*(a + f1 + f2)
    f1 = a**2/2.
    f2 = 5. - t + 9.*c + 4.*c**2
    f2 = f2*a**4/24.
    f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
    f3 = f3*a**6/720.
    yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
  endif
  xx = xx + east
  yy = yy + north

!
!---- UTM to Lat/Lon conversion
!
  else

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
  endif

  if (iway == IUTM2LONGLAT) then
    rlon = dlon
    rlat = dlat
    rx = rx_save
    ry = ry_save
  else
    rx = xx
    ry = yy
    rlon = rlon_save
    rlat = rlat_save
  endif

  end subroutine utm_geo

