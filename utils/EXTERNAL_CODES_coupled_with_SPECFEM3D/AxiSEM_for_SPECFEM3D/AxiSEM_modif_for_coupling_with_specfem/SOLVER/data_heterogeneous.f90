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

!===================
module data_heterogeneous
!===================

use global_parameters
implicit none
public 

! Heterogeneous region
   integer  :: num_het
   character(len=10), allocatable :: het_format(:),het_funct_type(:)
   character(len=200), allocatable :: het_file_discr(:), het_ani_discr(:), het_rel_discr(:)
   logical, allocatable :: rdep(:),grad(:)
   logical :: add_up, ani_hetero
   real(kind=dp)   , allocatable :: gradrdep1(:),gradrdep2(:)
   real(kind=dp)   , allocatable :: p_inv_dist(:), R_inv_dist(:)
   real(kind=dp)   , allocatable :: r_het1(:),r_het2(:),th_het1(:),th_het2(:)
   real(kind=dp)   , allocatable :: delta_rho(:), delta_vp(:), delta_vs(:)
   real(kind=dp)   , allocatable :: a_ica(:), b_ica(:), c_ica(:), fa_theta_ica(:), &
                                    fa_phi_ica(:), theta_slices(:)
   integer, allocatable :: inverseshape(:)
   integer :: num_slices
   real(kind=dp)    :: rhetmin, rhetmax, thhetmin, thhetmax



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!=======================
end module data_heterogeneous
!=======================
