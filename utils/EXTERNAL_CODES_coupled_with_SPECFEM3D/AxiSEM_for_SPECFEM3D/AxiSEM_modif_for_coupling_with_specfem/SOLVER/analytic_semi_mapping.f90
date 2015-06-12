!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!==================================
module analytic_semi_mapping
!==================================
!
!	10/01/2002: This module contains the 
! machinery necessary to describe analytically
! the transformation of the reference element
! into its deformed image in a semi-curved
! enveloppe. 
 
  use global_parameters 
  implicit none
  public :: map_semino, compute_partial_d_semino
  public :: map_semiso, compute_partial_d_semiso
  private

contains

!-----------------------------------------------------------------------------------------
real(kind=dp)    function map_semino(xil, etal, nodes_crd, idir)
! We are working in polar coordinates here: theta is the latitude. 

  real(kind=dp)     :: xil, etal
  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd
  integer           :: idir

  real(kind=dp)    :: atop, btop
  real(kind=dp)    :: thetabartop, dthetatop 
  real(kind=dp)    :: sbot, zbot, stop, ztop   ! careful: varname stop!!!
  real(kind=dp)    :: sbar, ds, slope, intersect

  map_semino = zero
  call compute_parameters_semino(nodes_crd, atop, btop, thetabartop, dthetatop)

  call compute_sz_xi_sline_no(sbot, zbot, xil, nodes_crd)
  call compute_sz_xi(stop, ztop, xil, atop, btop, thetabartop, dthetatop)

  sbar = half*(sbot+stop)
  ds = stop-sbot

  if (abs(ds)>1.d-10) then
     intersect = (zbot*stop-ztop*sbot)/ds   
     slope = (ztop-zbot)/ds
  end if   
  
  if (idir == 1) then
     map_semino = sbar+ds*etal*half
  elseif (idir == 2) then
     if (abs(ds)>1.d-10) then
     map_semino = slope*(sbar+half*ds*etal)+intersect 
     else
     map_semino = half*(zbot+ztop)+etal*(ztop-zbot)*half
     end if
  end if

end function map_semino
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp)    function map_semiso(xil,etal,nodes_crd,idir)
! We are working in polar coordinates here: theta is the latitude. 
 
  real(kind=dp)    :: xil, etal
  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd
  integer           :: idir
  
  real(kind=dp)    :: abot, bbot
  real(kind=dp)    :: thetabarbot, dthetabot
  real(kind=dp)    :: sbot, zbot, stop, ztop
  real(kind=dp)    :: sbar, ds, slope, intersect
! 
  map_semiso = zero
  call compute_parameters_semiso(nodes_crd, abot, bbot, thetabarbot, dthetabot)
  
  call compute_sz_xi_sline_so(stop, ztop, xil, nodes_crd)
  call compute_sz_xi(sbot, zbot, xil, abot, bbot, thetabarbot, dthetabot)
  
  sbar = half*(sbot+stop)
  ds = stop-sbot
  
  if (abs(ds)>1.d-10) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = (ztop-zbot)/ds
  end if
  
  if (idir == 1) then
     map_semiso = sbar+ds*etal*half
  elseif (idir == 2) then
     if (abs(ds)>1.d-10) then
     map_semiso = slope*(sbar+half*ds*etal)+intersect
     else
     map_semiso = half*(zbot+ztop)+etal*(ztop-zbot)*half
     end if
  end if

end function map_semiso
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_partial_d_semino(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, nodes_crd)

  real(kind=dp)   , intent(out) :: dsdxi, dzdxi, dsdeta, dzdeta
  real(kind=dp)   , intent(in)  :: xil, etal
  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd

  real(kind=dp)    :: atop, btop
  real(kind=dp)    :: thetabartop, dthetatop
  real(kind=dp)    :: sbot, zbot, stop, ztop
  real(kind=dp)    :: sbar, ds, dz, slope, intersect, sxieta
  real(kind=dp)    :: dsbotdxi, dzbotdxi
  real(kind=dp)    :: dstopdxi, dztopdxi
  real(kind=dp)    :: dsbardxi, ddsdxi
  real(kind=dp)    :: dzbardxi, ddzdxi
  real(kind=dp)    :: dslopedxi, dintersectdxi

  call compute_parameters_semino(nodes_crd, atop, btop, thetabartop, dthetatop)

  call compute_sz_xi_sline_no(sbot, zbot, xil, nodes_crd)
  call compute_sz_xi(stop, ztop, xil, atop, btop, thetabartop, dthetatop)

  sbar = half*(sbot+stop)
  ds = stop-sbot
  dz = ztop - zbot
  
  if (abs(ds)>1.d-10) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = dz/ds
  end if

  call compute_dsdxi_dzdxi_sline_no(dsbotdxi,dzbotdxi,nodes_crd)
  call compute_dsdxi_dzdxi(dstopdxi,dztopdxi,xil,atop,btop,thetabartop,dthetatop)

  dsbardxi = half*(dsbotdxi+dstopdxi)
  ddsdxi = (dstopdxi-dsbotdxi)
 
  dzbardxi = half*(dzbotdxi+dztopdxi)
  ddzdxi = (dztopdxi-dzbotdxi) 

  sxieta = sbar+ds*etal*half
  
  dsdxi = dsbardxi + half*etal*ddsdxi
  dsdeta = half*ds
  
  if (dabs(ds)>1.d-10) then
     dslopedxi     = (ddzdxi*ds-ddsdxi*dz)/ds**2
     dintersectdxi = ((dzbotdxi*stop-dztopdxi*sbot     &
                       +zbot*dstopdxi-ztop*dsbotdxi)*ds &
                      -ddsdxi*(zbot*stop-ztop*sbot))/ds**2
     dzdxi  = (dz/ds)*dsdxi+ sxieta*dslopedxi + dintersectdxi
     dzdeta = (dz/ds)*dsdeta  
  else
     write(6,*)'DS SMALL IN SEMIN ELEMENT!!!'; call flush(6)
     dzdxi = zero
     dzdeta = half*dz
  end if

end subroutine compute_partial_d_semino
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_partial_d_semiso(dsdxi, dzdxi, dsdeta, dzdeta, xil, etal, nodes_crd)

  real(kind=dp)   ,  intent(out) :: dsdxi, dzdxi, dsdeta, dzdeta
  real(kind=dp)   ,  intent(in)  :: xil, etal
  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd

  real(kind=dp)    :: abot,bbot
  real(kind=dp)    :: thetabarbot,dthetabot
  real(kind=dp)    :: sbot,zbot,stop,ztop
  real(kind=dp)    :: sbar,ds,dz,slope,intersect,sxieta
  real(kind=dp)    :: dsbotdxi,dzbotdxi
  real(kind=dp)    :: dstopdxi,dztopdxi
  real(kind=dp)    :: dsbardxi,ddsdxi
  real(kind=dp)    :: dzbardxi,ddzdxi
  real(kind=dp)    :: dslopedxi,dintersectdxi

  call compute_parameters_semiso(nodes_crd, abot, bbot, thetabarbot, dthetabot)
  call compute_sz_xi_sline_so(stop, ztop, xil, nodes_crd)
  call compute_sz_xi(sbot, zbot, xil, abot, bbot, thetabarbot, dthetabot)

  sbar = half*(sbot+stop); ds = stop-sbot ; dz = ztop - zbot
  if (abs(ds)>1.d-10) then
     intersect = (zbot*stop-ztop*sbot)/ds
     slope = dz/ds
  end if

  call compute_dsdxi_dzdxi(dsbotdxi, dzbotdxi, xil, abot, bbot, thetabarbot, dthetabot)
  call compute_dsdxi_dzdxi_sline_so(dstopdxi, dztopdxi, nodes_crd)

  dsbardxi = half*(dsbotdxi+dstopdxi)
  ddsdxi = (dstopdxi-dsbotdxi)
 
  dzbardxi = half*(dzbotdxi+dztopdxi)
  ddzdxi = (dztopdxi-dzbotdxi) 

  sxieta = sbar+ds*etal*half
  dsdxi  = dsbardxi + half*etal*ddsdxi
  dsdeta = half*ds
  
  if (dabs(ds)>1.d-10) then
    dslopedxi     = (ddzdxi*ds-ddsdxi*dz)/ds**2
    dintersectdxi = ((dzbotdxi*stop-dztopdxi*sbot     &
                       +zbot*dstopdxi-ztop*dsbotdxi)*ds &
                      -ddsdxi*(zbot*stop-ztop*sbot))/ds**2
    dzdxi  = (dz/ds)*dsdxi+ sxieta*dslopedxi + dintersectdxi
    dzdeta = (dz/ds)*dsdeta  
  else
    dzdxi  = zero
    dzdeta = half*dz
  end if

end subroutine compute_partial_d_semiso
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_sz_xi(s,z,xil,a,b,thetabar,dtheta)

  real(kind=dp)   , intent(out) :: s,z
  real(kind=dp)   , intent(in)  :: xil, a, b, thetabar, dtheta
  
  s = a*dcos(thetabar+xil*half*dtheta)
  z = b*dsin(thetabar+xil*half*dtheta)

end subroutine compute_sz_xi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_sz_xi_sline_no(s, z, xil, nodes_crd)

  real(kind=dp)   , intent(out) :: s,z
  real(kind=dp)   , intent(in)  :: xil,nodes_crd(8,2)
  real(kind=dp)    :: s1, z1, s3, z3
  
  s1 = nodes_crd(1,1)
  z1 = nodes_crd(1,2)
  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)

  s = half*((one+xil)*s3+(one-xil)*s1)
  z = half*((one+xil)*z3+(one-xil)*z1)

end subroutine compute_sz_xi_sline_no
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_sz_xi_sline_so(s, z, xil, nodes_crd)

  real(kind=dp)   , intent(out) :: s,z
  real(kind=dp)   , intent(in)  :: xil,nodes_crd(8,2)
  real(kind=dp)    :: s7, z7, s5, z5

  s7 = nodes_crd(7,1) ; z7 = nodes_crd(7,2)
  s5 = nodes_crd(5,1) ; z5 = nodes_crd(5,2)

  s = half*((one+xil)*s5+(one-xil)*s7)
  z = half*((one+xil)*z5+(one-xil)*z7)

end subroutine compute_sz_xi_sline_so
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_dsdxi_dzdxi(dsdxi, dzdxi, xil, a, b, thetabar, dtheta)
  
  real(kind=dp)   , intent(out) :: dsdxi, dzdxi
  real(kind=dp)   , intent(in)  :: xil, a, b, thetabar, dtheta 

  dsdxi =-a*half*dtheta*dsin(thetabar+xil*half*dtheta)
  dzdxi = b*half*dtheta*dcos(thetabar+xil*half*dtheta)  
 
end subroutine compute_dsdxi_dzdxi
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_dsdxi_dzdxi_sline_no(dsdxi,dzdxi,nodes_crd)

  real(kind=dp)   , intent(out) :: dsdxi, dzdxi
  real(kind=dp)   , intent(in)  :: nodes_crd(8,2)
  real(kind=dp)    :: s1, z1, s3, z3
  
  s1 = nodes_crd(1,1)
  z1 = nodes_crd(1,2)
  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)
  
  dsdxi = half*(s3-s1) 
  dzdxi = half*(z3-z1)

end subroutine compute_dsdxi_dzdxi_sline_no
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_dsdxi_dzdxi_sline_so(dsdxi, dzdxi, nodes_crd)
  
  real(kind=dp)   , intent(out) :: dsdxi, dzdxi
  real(kind=dp)   , intent(in)  :: nodes_crd(8,2)
  real(kind=dp)    :: s7, z7, s5, z5
  
  s7 = nodes_crd(7,1)
  z7 = nodes_crd(7,2)
  s5 = nodes_crd(5,1)
  z5 = nodes_crd(5,2)
  
  dsdxi = half*(s5-s7) ; dzdxi = half*(z5-z7)

end subroutine compute_dsdxi_dzdxi_sline_so
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_parameters_semino(nodes_crd, atop, btop, thetabartop, dthetatop)

  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd
  real(kind=dp)   ,intent(out) :: atop,btop
  real(kind=dp)   ,intent(out) :: thetabartop, dthetatop
  real(kind=dp)    :: s1, z1, s3, z3, s5, z5, s7, z7
  real(kind=dp)    :: theta5, theta7
  
  s1 = nodes_crd(1,1)
  z1 = nodes_crd(1,2)
  
  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)
  
  s5 = nodes_crd(5,1)
  z5 = nodes_crd(5,2)
  
  s7 = nodes_crd(7,1)
  z7 = nodes_crd(7,2)

  call compute_ab(atop, btop, s7, z7, s5, z5)
  call compute_theta(theta5, s5, z5, atop, btop) 
  call compute_theta(theta7, s7, z7, atop, btop) 

  thetabartop = half*(theta7+theta5)
  dthetatop = theta5 -theta7

end subroutine compute_parameters_semino
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_parameters_semiso(nodes_crd, abot, bbot, thetabarbot, dthetabot)

  real(kind=dp)   , dimension(8,2),intent(in) :: nodes_crd
  real(kind=dp)   ,intent(out) :: abot,bbot
  real(kind=dp)   ,intent(out) :: thetabarbot, dthetabot
  real(kind=dp)    :: s1, z1, s3, z3, s5, z5, s7, z7
  real(kind=dp)    :: theta1, theta3
  
  s1 = nodes_crd(1,1) 
  z1 = nodes_crd(1,2)
  
  s3 = nodes_crd(3,1)
  z3 = nodes_crd(3,2)
  
  s5 = nodes_crd(5,1)
  z5 = nodes_crd(5,2)
  
  s7 = nodes_crd(7,1)
  z7 = nodes_crd(7,2)

  call compute_ab(abot, bbot, s1, z1, s3, z3)
  call compute_theta(theta3, s3, z3, abot, bbot)
  call compute_theta(theta1, s1, z1, abot, bbot)

  thetabarbot = half*(theta1+theta3)
  dthetabot = theta3 -theta1

end subroutine compute_parameters_semiso
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_ab(a,b,s1,z1,s2,z2)

  real(kind=dp)   , intent(out) :: a,b
  real(kind=dp)   , intent(in) :: s1,z1,s2,z2

  a = dsqrt(dabs((s2**2*z1**2-z2**2*s1**2)/(z1**2-z2**2)))
  b = dsqrt(dabs((z1**2*s2**2-z2**2*s1**2)/(s2**2-s1**2))) 
end subroutine compute_ab
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_theta(theta,s,z,a,b)
!	This routine returns the latitude theta, given s and z.

  real(kind=dp)   , intent(out) :: theta
  real(kind=dp)   , intent(in) :: s,z,a,b

  if (s /= zero) then
     theta=datan(z*a/(s*b))  
  else
     if (z>0) theta=half*pi
     if (z<0) theta=-half*pi
  end if

end subroutine compute_theta
!-----------------------------------------------------------------------------------------


!=======================================
end module analytic_semi_mapping
!=======================================
