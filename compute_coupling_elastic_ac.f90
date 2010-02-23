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

! for elastic solver

  subroutine compute_coupling_elastic_ac(NSPEC_AB,NGLOB_AB, &
                        ibool,accel,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! returns the updated acceleration array: accel                        

  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  
! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces) 
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)   

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! local parameters
  real(kind=CUSTOM_REAL) :: pressure
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  
  integer :: iface,igll,ispec,iglob
  integer :: i,j,k
  
! loops on all coupling faces
  do iface = 1,num_coupling_ac_el_faces

    ! gets corresponding spectral element 
    ! (note: can be either acoustic or elastic element, no need to specify since
    !           no material properties are needed for this coupling term)
    ispec = coupling_ac_el_ispec(iface)

    if( ispec_is_inner(ispec) .eqv. phase_is_inner ) then
    
      ! loops over common GLL points
      do igll = 1, NGLLSQUARE
        i = coupling_ac_el_ijk(1,igll,iface)
        j = coupling_ac_el_ijk(2,igll,iface)
        k = coupling_ac_el_ijk(3,igll,iface)
        
        ! gets global index of this common GLL point
        ! (note: should be the same as for corresponding i',j',k',ispec_elastic or ispec_elastic )
        iglob = ibool(i,j,k,ispec)
        
        ! acoustic pressure on global point
        pressure = - potential_dot_dot_acoustic(iglob)

        ! gets associated normal on GLL point
        ! (note convention: pointing outwards of acoustic element)
        nx = coupling_ac_el_normal(1,igll,iface)
        ny = coupling_ac_el_normal(2,igll,iface)
        nz = coupling_ac_el_normal(3,igll,iface)                   
        
        ! gets associated, weighted 2D jacobian 
        ! (note: should be the same for elastic and acoustic element)
        jacobianw = coupling_ac_el_jacobian2Dw(igll,iface)
        
        ! continuity of displacement and pressure on global point
        !
        ! note: newark time scheme together with definition of scalar potential: 
        !          pressure = - chi_dot_dot
        !          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
        !          pressure at time step [t + delta_t] 
        !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
        !          it means you have to calculate and update the acoustic pressure first before
        !          calculating this term...
        accel(1,iglob) = accel(1,iglob) + jacobianw*nx*pressure
        accel(2,iglob) = accel(2,iglob) + jacobianw*ny*pressure
        accel(3,iglob) = accel(3,iglob) + jacobianw*nz*pressure
        
      enddo ! igll

    endif
    
  enddo ! iface

end subroutine compute_coupling_elastic_ac

