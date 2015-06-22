!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
! read and smooth crust2.0 model
! based on software routines provided with the crust2.0 model by Bassin et al.
!

  subroutine crustal_model(lat,lon,x,vp,vs,rho,moho,found_crust,CM_V)

  implicit none
  include "constants.h"

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

  double precision lat,lon,x,vp,vs,rho,moho
  logical found_crust

  double precision h_sed,h_uc
  double precision x3,x4,x5,x6,x7,scaleval
  double precision vps(NLAYERS_CRUST),vss(NLAYERS_CRUST),rhos(NLAYERS_CRUST),thicks(NLAYERS_CRUST)

  call crust(lat,lon,vps,vss,rhos,thicks,CM_V%abbreviation,CM_V%code,CM_V%thlr,CM_V%velocp,CM_V%velocs,CM_V%dens)

 x3 = (R_EARTH-thicks(3)*1000.0d0)/R_EARTH
 h_sed = thicks(3) + thicks(4)
 x4 = (R_EARTH-h_sed*1000.0d0)/R_EARTH
 h_uc = h_sed + thicks(5)
 x5 = (R_EARTH-h_uc*1000.0d0)/R_EARTH
 x6 = (R_EARTH-(h_uc+thicks(6))*1000.0d0)/R_EARTH
 x7 = (R_EARTH-(h_uc+thicks(6)+thicks(7))*1000.0d0)/R_EARTH

 found_crust = .true.
 if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST) then
   vp = vps(3)
   vs = vss(3)
   rho = rhos(3)
 else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST) then
   vp = vps(4)
   vs = vss(4)
   rho = rhos(4)
 else if(x > x5) then
   vp = vps(5)
   vs = vss(5)
   rho = rhos(5)
 else if(x > x6) then
   vp = vps(6)
   vs = vss(6)
   rho = rhos(6)
 else if(x > x7) then
   vp = vps(7)
   vs = vss(7)
   rho = rhos(7)
 else
   found_crust = .false.
 endif

 if (found_crust) then
!   non-dimensionalize
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vp = vp*1000.0d0/(R_EARTH*scaleval)
    vs = vs*1000.0d0/(R_EARTH*scaleval)
    rho = rho*1000.0d0/RHOAV
    moho = (h_uc+thicks(6)+thicks(7))*1000.0d0/R_EARTH
 endif

 end subroutine crustal_model

!---------------------------

  subroutine read_crustal_model(CM_V)

  implicit none
  include "constants.h"

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

! local variables
  integer i
  integer ila,icolat
  integer ikey

  double precision h_moho_min,h_moho_max

  character(len=150) CNtype2, CNtype2_key_modif

  call get_value_string(CNtype2, 'model.CNtype2', 'DATA/crust2.0/CNtype2.txt')
  call get_value_string(CNtype2_key_modif, 'model.CNtype2_key_modif', 'DATA/crust2.0/CNtype2_key_modif.txt')

  open(unit=1,file=CNtype2,status='old',action='read')
  do ila=1,NCAP_CRUST/2
    read(1,*) icolat,(CM_V%abbreviation(ila,i),i=1,NCAP_CRUST)
  enddo
  close(1)

  open(unit=1,file=CNtype2_key_modif,status='old',action='read')
  h_moho_min=HUGEVAL
  h_moho_max=-HUGEVAL
  do ikey=1,NKEYS_CRUST
    read (1,"(a2)") CM_V%code(ikey)
    read (1,*) (CM_V%velocp(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%velocs(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%dens(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%thlr(ikey,i),i=1,NLAYERS_CRUST-1),CM_V%thlr(ikey,NLAYERS_CRUST)
    if(CM_V%thlr(ikey,NLAYERS_CRUST) > h_moho_max) h_moho_max=CM_V%thlr(ikey,NLAYERS_CRUST)
    if(CM_V%thlr(ikey,NLAYERS_CRUST) < h_moho_min) h_moho_min=CM_V%thlr(ikey,NLAYERS_CRUST)
  enddo
  close(1)

  if(h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) &
    stop 'incorrect moho depths in read_3D_crustal_model'

  end subroutine read_crustal_model

!---------------------------

  subroutine crust(lat,lon,velp,vels,rho,thick,abbreviation,code,thlr,velocp,velocs,dens)

! crustal vp and vs in km/s, layer thickness in km
! crust2.0 is smoothed with a cap of size CAP using NTHETA points
! in the theta direction and NPHI in the phi direction.
! The cap is rotated to the North Pole.

  implicit none
  include "constants.h"

  integer, parameter :: NTHETA = 2
  integer, parameter :: NPHI = 10
  double precision, parameter :: CAP = 2.0d0*PI/180.0d0

! argument variables
  double precision lat,lon
  double precision rho(NLAYERS_CRUST),thick(NLAYERS_CRUST),velp(NLAYERS_CRUST),vels(NLAYERS_CRUST)
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
  character(len=2) code(NKEYS_CRUST),abbreviation(NCAP_CRUST/2,NCAP_CRUST)

! local variables
  integer i,j,k,icolat,ilon,ierr
  integer itheta,iphi,npoints
  double precision theta,phi,sint,cost,sinp,cosp,dtheta,dphi,cap_area,wght,total
  double precision r_rot,theta_rot,phi_rot
  double precision rotation_matrix(3,3),x(3),xc(3)
  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)
  double precision rhol(NLAYERS_CRUST),thickl(NLAYERS_CRUST),velpl(NLAYERS_CRUST),velsl(NLAYERS_CRUST)
  character(len=2) crustaltype

! get integer colatitude and longitude of crustal cap
! -90<lat<90 -180<lon<180
  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
    stop 'error in latitude/longitude range in crust'
  if(lat==90.0d0) lat=89.9999d0
  if(lat==-90.0d0) lat=-89.9999d0
  if(lon==180.0d0) lon=179.9999d0
  if(lon==-180.0d0) lon=-179.9999d0

  call icolat_ilon(lat,lon,icolat,ilon)
  crustaltype=abbreviation(icolat,ilon)
  call get_crust_structure(crustaltype,velp,vels,rho,thick, &
                    code,thlr,velocp,velocs,dens,ierr)

!  uncomment the following line to use crust2.0 as is, without smoothing
!
!  return

  theta = (90.0-lat)*PI/180.0
  phi = lon*PI/180.0

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

! set up rotation matrix to go from cap at North pole
! to cap around point of interest
  rotation_matrix(1,1) = cosp*cost
  rotation_matrix(1,2) = -sinp
  rotation_matrix(1,3) = cosp*sint
  rotation_matrix(2,1) = sinp*cost
  rotation_matrix(2,2) = cosp
  rotation_matrix(2,3) = sinp*sint
  rotation_matrix(3,1) = -sint
  rotation_matrix(3,2) = 0.0
  rotation_matrix(3,3) = cost

  dtheta = CAP/dble(NTHETA)
  dphi = 2.0*PI/dble(NPHI)
  cap_area = 2.0*PI*(1.0-cos(CAP))

! integrate over a cap at the North pole
  i = 0
  total = 0.0
  do itheta = 1,NTHETA

    theta = 0.5*dble(2*itheta-1)*CAP/dble(NTHETA)
    cost = cos(theta)
    sint = sin(theta)
    wght = sint*dtheta*dphi/cap_area

    do iphi = 1,NPHI

      i = i+1
!     get the weight associated with this integration point (same for all phi)
      weight(i) = wght
      total = total + weight(i)
      phi = dble(2*iphi-1)*PI/dble(NPHI)
      cosp = cos(phi)
      sinp = sin(phi)
!     x,y,z coordinates of integration point in cap at North pole
      xc(1) = sint*cosp
      xc(2) = sint*sinp
      xc(3) = cost
!     get x,y,z coordinates in cap around point of interest
      do j=1,3
        x(j) = 0.0
        do k=1,3
          x(j) = x(j)+rotation_matrix(j,k)*xc(k)
        enddo
      enddo
!     get latitude and longitude (degrees) of integration point
      call xyz_2_rthetaphi_dble(x(1),x(2),x(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      xlat(i) = (PI/2.0-theta_rot)*180.0/PI
      xlon(i) = phi_rot*180.0/PI
      if(xlon(i) > 180.0) xlon(i) = xlon(i)-360.0

    enddo

  enddo

  if(abs(total-1.0) > 0.001) stop 'error in cap integration for crust2.0'

  npoints = i

  do j=1,NLAYERS_CRUST
    rho(j)=0.0d0
    thick(j)=0.0d0
    velp(j)=0.0d0
    vels(j)=0.0d0
  enddo

  do i=1,npoints
    call icolat_ilon(xlat(i),xlon(i),icolat,ilon)
    crustaltype=abbreviation(icolat,ilon)
    call get_crust_structure(crustaltype,velpl,velsl,rhol,thickl, &
                    code,thlr,velocp,velocs,dens,ierr)
    if(ierr /= 0) stop 'error in routine get_crust_structure'
    do j=1,NLAYERS_CRUST
      rho(j)=rho(j)+weight(i)*rhol(j)
      thick(j)=thick(j)+weight(i)*thickl(j)
      velp(j)=velp(j)+weight(i)*velpl(j)
      vels(j)=vels(j)+weight(i)*velsl(j)
    enddo
  enddo

  end subroutine crust

!------------------------------------------------------

  subroutine icolat_ilon(xlat,xlon,icolat,ilon)

  implicit none


! argument variables
  double precision xlat,xlon
  integer icolat,ilon

  if(xlat > 90.0d0 .or. xlat < -90.0d0 .or. xlon > 180.0d0 .or. xlon < -180.0d0) &
    stop 'error in latitude/longitude range in icolat_ilon'
  icolat=int(1+((90.d0-xlat)/2.d0))
  if(icolat == 91) icolat=90
  ilon=int(1+((180.d0+xlon)/2.d0))
  if(ilon == 181) ilon=1

  if(icolat>90 .or. icolat<1) stop 'error in routine icolat_ilon'
  if(ilon<1 .or. ilon>180) stop 'error in routine icolat_ilon'

  end subroutine icolat_ilon

!---------------------------------------------------------------------

  subroutine get_crust_structure(type,vptyp,vstyp,rhtyp,thtp, &
               code,thlr,velocp,velocs,dens,ierr)

  implicit none
  include "constants.h"


! argument variables
  integer ierr
  double precision rhtyp(NLAYERS_CRUST),thtp(NLAYERS_CRUST)
  double precision vptyp(NLAYERS_CRUST),vstyp(NLAYERS_CRUST)
  character(len=2) type,code(NKEYS_CRUST)
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)

! local variables
  integer i,ikey

  ierr=1
  do ikey=1,NKEYS_CRUST
  if (code(ikey) == type) then
    do i=1,NLAYERS_CRUST
      vptyp(i)=velocp(ikey,i)
      vstyp(i)=velocs(ikey,i)
      rhtyp(i)=dens(ikey,i)
    enddo
    do i=1,NLAYERS_CRUST-1
      thtp(i)=thlr(ikey,i)
    enddo
!   get distance to Moho from the bottom of the ocean or the ice
    thtp(NLAYERS_CRUST)=thlr(ikey,NLAYERS_CRUST)-thtp(1)-thtp(2)
    ierr=0
  endif
  enddo

  end subroutine get_crust_structure

