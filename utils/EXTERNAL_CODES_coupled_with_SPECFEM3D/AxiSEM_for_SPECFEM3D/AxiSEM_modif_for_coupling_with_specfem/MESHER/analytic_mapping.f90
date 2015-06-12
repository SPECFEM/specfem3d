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

module analytic_mapping
  
  use analytic_spheroid_mapping
  use analytic_semi_mapping
  use subpar_mapping
  use data_mesh

  implicit none

  public :: mapping_anal
 
  private

contains

!-----------------------------------------------------------------------------------------
real(kind=dp)    function mapping_anal(xi,eta,nodes_crd,iaxis,ielem0)
  
  integer           :: iaxis,ielem0
  real(kind=dp)     :: xi, eta, nodes_crd(8,2)
  if (eltypeg(ielem0) == 'curved') mapping_anal = map_spheroid(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'linear') mapping_anal = mapping_subpar(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'semino') mapping_anal = map_semino(xi,eta,nodes_crd,iaxis)
  if (eltypeg(ielem0) == 'semiso') mapping_anal = map_semiso(xi,eta,nodes_crd,iaxis)

end function mapping_anal
!-----------------------------------------------------------------------------------------

end module analytic_mapping
