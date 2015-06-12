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

module data_spec

  use global_parameters, only                   : sp, dp
  implicit none

  public 

  integer :: npol
  real(kind=dp)   , dimension(:),allocatable   :: xi_k, eta
  real(kind=dp)   , dimension(:),allocatable   :: dxi
  real(kind=dp)   , dimension (:), allocatable :: wt          !Quadrature weights
  real(kind=dp)   , dimension (:), allocatable :: wt_axial_k  !Quad. wgts for the 
                                                              !nonaxisymmetric components
end module data_spec
