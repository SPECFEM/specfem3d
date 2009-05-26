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
!
! United States and French Government Sponsorship Acknowledged.


! compute the approximate amount of static memory needed to run the solver

 subroutine memory_eval(NSPEC_AB,NGLOB_AB,max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh,static_memory_size)

  implicit none

  include "constants.h"

! input
!  logical, intent(in) :: ATTENUATION
  integer, intent(in) :: NSPEC_AB,NGLOB_AB
  integer, intent(in) :: max_nibool_interfaces_ext_mesh,ninterfaces_ext_mesh

! output
  double precision, intent(out) :: static_memory_size


  static_memory_size = 0.d0

! add size of each set of static arrays multiplied by the number of such arrays

! ibool,idoubling
  static_memory_size = static_memory_size + 2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(SIZE_INTEGER)

! xix,xiy,xiz,
! etax,etay,etaz,
! gammax,gammay,gammaz,jacobian
! kappavstore,muvstore
! flag_sediments,rho_vp,rho_vs  
  static_memory_size = static_memory_size + 15.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_AB*dble(CUSTOM_REAL)

! xstore,ystore,zstore,rmass,rmass_ocean_load
  static_memory_size = static_memory_size + 5.d0*NGLOB_AB*dble(CUSTOM_REAL)

! updated_dof_ocean_load,iglob_is_inner_ext_mesh 
  static_memory_size = static_memory_size + 2.d0*NGLOB_AB*dble(SIZE_LOGICAL)

! ispec_is_inner_ext_mesh 
  static_memory_size = static_memory_size + NSPEC_AB*dble(SIZE_LOGICAL)

! displ,veloc,accel
  static_memory_size = static_memory_size + 3.d0*dble(NDIM)*NGLOB_AB*dble(CUSTOM_REAL)

! my_neighbours_ext_mesh,nibool_interfaces_ext_mesh
  static_memory_size = static_memory_size + 2.d0*ninterfaces_ext_mesh*dble(SIZE_INTEGER)

! ibool_interfaces_ext_mesh
 static_memory_size = static_memory_size + max_nibool_interfaces_ext_mesh*ninterfaces_ext_mesh*dble(SIZE_INTEGER)

! buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh  
 static_memory_size = static_memory_size + 2.d0*dble(NDIM)*max_nibool_interfaces_ext_mesh*ninterfaces_ext_mesh*dble(CUSTOM_REAL)

! buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh 
 static_memory_size = static_memory_size + 2.d0*max_nibool_interfaces_ext_mesh*ninterfaces_ext_mesh*dble(CUSTOM_REAL)

! request_send_vector_ext_mesh,request_recv_vector_ext_mesh,request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh 
 static_memory_size = static_memory_size + 4.d0*ninterfaces_ext_mesh*dble(SIZE_INTEGER)


  end subroutine memory_eval



