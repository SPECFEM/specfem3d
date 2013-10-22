!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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

! for acoustic solver

  subroutine compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                        ibool,displ,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner,&
                        PML_CONDITIONS,spec_to_CPML,is_CPML,&
                        rmemory_coupling_ac_el_displ,&
                        SIMULATION_TYPE,backward_simulation)

! returns the updated pressure array: potential_dot_dot_acoustic
  use pml_par, only: NSPEC_CPML
  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  integer :: SIMULATION_TYPE
  logical :: backward_simulation

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces)
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2) :: rmemory_coupling_ac_el_displ

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! CPML
  logical :: PML_CONDITIONS
  integer :: spec_to_CPML(NSPEC_AB)
  logical :: is_CPML(NSPEC_AB)

! local parameters
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z,displ_n
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  integer :: iface,igll,ispec,iglob,ispec_CPML,i,j,k

! loops on all coupling faces
  do iface = 1,num_coupling_ac_el_faces

    ! gets corresponding elements
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
        ! (note: should be the same as for corresponding i',j',k',ispec_elastic or ispec_acoustic)
        iglob = ibool(i,j,k,ispec)

        ! elastic displacement on global point
        if(PML_CONDITIONS .and. NSPEC_CPML > 0)then
           if(.not. backward_simulation)then
               if(is_CPML(ispec))then
                  if(SIMULATION_TYPE == 1)then
                    ispec_CPML = spec_to_CPML(ispec)
                    call pml_compute_memory_variables_acoustic_elastic(ispec_CPML,iface,iglob,i,j,k,&
                                                      displ_x,displ_y,displ_z,displ,&
                                                      num_coupling_ac_el_faces,rmemory_coupling_ac_el_displ)
                  endif

                  if(SIMULATION_TYPE == 3)then
!left blank for change
                  endif

               else
                  displ_x = displ(1,iglob)
                  displ_y = displ(2,iglob)
                  displ_z = displ(3,iglob)
               endif
           else
              if(is_CPML(ispec))then
! left blank, since no operation needed
              else
                displ_x = displ(1,iglob)
                displ_y = displ(2,iglob)
                displ_z = displ(3,iglob)
              endif
           endif
        else
           displ_x = displ(1,iglob)
           displ_y = displ(2,iglob)
           displ_z = displ(3,iglob)
        endif

        ! gets associated normal on GLL point
        ! (note convention: pointing outwards of acoustic element)
        nx = coupling_ac_el_normal(1,igll,iface)
        ny = coupling_ac_el_normal(2,igll,iface)
        nz = coupling_ac_el_normal(3,igll,iface)

        ! calculates displacement component along normal
        ! (normal points outwards of acoustic element)
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz

        ! gets associated, weighted jacobian
        jacobianw = coupling_ac_el_jacobian2Dw(igll,iface)

        ! continuity of pressure and normal displacement on global point
        !
        ! note: Newmark time scheme together with definition of scalar potential:
        !          pressure = - chi_dot_dot
        !          requires that this coupling term uses the updated displacement at time step [t+delta_t],
        !          which is done at the very beginning of the time loop
        !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
        !          it also means you have to calculate and update this here first before
        !          calculating the coupling on the elastic side for the acceleration...
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + jacobianw*displ_n

      enddo ! igll

    endif

  enddo ! iface

end subroutine compute_coupling_acoustic_el
