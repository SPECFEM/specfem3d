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

module data_grid

  use global_parameters, only: sp, dp
  implicit none
  public
  
  integer           :: ns, nz
  real(kind=dp)     :: ri, ro
  real(kind=dp), dimension(:,:,:),allocatable :: crd_grdc, crd_grds
  real(kind=dp), dimension(:), allocatable    :: s_unif, z_unif
  real(kind=dp), dimension(:), allocatable    :: radius
  real(kind=dp), dimension(:), allocatable    :: ndeta
  
  ! Shrinking factor for axis elements
  real(kind=dp)                               :: axisfac
  ! Shrinking factor for fluid elements
  real(kind=dp)                               :: fluidfac
  
  !number of subdivisions for central square
  integer            :: ndivs 
  real(kind=dp)      :: lsq
  
  real(kind=dp)      :: router,lsq_fac
  integer            :: ngllcube
  
  real(kind=dp), dimension(:), allocatable :: rmax_el, rmin_el
  logical :: southern
end module data_grid
