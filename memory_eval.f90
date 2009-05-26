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

  integer, parameter :: SIZE_INTEGER = SIZE_REAL
  integer, parameter :: SIZE_LOGICAL = SIZE_REAL


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





!   subroutine memory_eval(ATTENUATION,NSPEC,nglob,static_memory_size)

!   implicit none

!   include "constants.h"

! ! input
!   logical, intent(in) :: ATTENUATION
!   integer, dimension(MAX_NUM_REGIONS), intent(in) :: NSPEC, nglob

! ! output
!   double precision, intent(out) :: static_memory_size

! ! variables
!   integer :: NSPEC_INNER_CORE_ATTENUATION,NUMBER_OF_MESH_LAYERS

!   if (ONE_CRUST) then
!     NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
!   else
!   NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
!   endif

!   if(ATTENUATION) then
!     NSPEC_INNER_CORE_ATTENUATION = NSPEC(IREGION_INNER_CORE)
!   else
!     NSPEC_INNER_CORE_ATTENUATION = 1
!   endif

! ! add size of each set of static arrays multiplied by the number of such arrays

!   static_memory_size = 0.d0

! ! R_memory_inner_core
!   static_memory_size = static_memory_size + 5.d0*dble(N_SLS)*dble(NGLLX)* &
!     dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)

! ! xix_outer_core,xiy_outer_core,xiz_outer_core,
! ! etax_outer_core,etay_outer_core,etaz_outer_core,
! ! gammax_outer_core,gammay_outer_core,gammaz_outer_core
! ! rhostore_outer_core,kappavstore_outer_core
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_OUTER_CORE)*11.d0*dble(CUSTOM_REAL)

! ! ibool_outer_core
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_OUTER_CORE)*dble(SIZE_INTEGER)

! ! xstore_outer_core, ystore_outer_core, zstore_outer_core, rmass_outer_core, displ_outer_core, veloc_outer_core, accel_outer_core
!   static_memory_size = static_memory_size + NGLOB(IREGION_OUTER_CORE)*7.d0*dble(CUSTOM_REAL)

! ! ibool_inner_core
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_INNER_CORE)*dble(SIZE_INTEGER)

! ! xix_inner_core,xiy_inner_core,xiz_inner_core,
! ! etax_inner_core,etay_inner_core,etaz_inner_core,
! ! gammax_inner_core,gammay_inner_core,gammaz_inner_core,
! ! rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core
!   static_memory_size = static_memory_size + dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC(IREGION_INNER_CORE)*12.d0*dble(CUSTOM_REAL)

! ! xstore_inner_core,ystore_inner_core,zstore_inner_core,rmass_inner_core
!   static_memory_size = static_memory_size + NGLOB(IREGION_INNER_CORE)*4.d0*dble(CUSTOM_REAL)

! ! displ_inner_core,veloc_inner_core,accel_inner_core
!   static_memory_size = static_memory_size + dble(NDIM)*NGLOB(IREGION_INNER_CORE)*3.d0*dble(CUSTOM_REAL)

!   end subroutine memory_eval
