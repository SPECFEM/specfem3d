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

! for acoustic solver

  subroutine compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                            potential_dot_dot_acoustic,potential_dot_acoustic, &
                            ibool,ispec_is_inner,phase_is_inner, &
                            abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic,&
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it,myrank,NGLOB_ADJOINT, &
                            b_potential_dot_dot_acoustic,b_reclen_potential, &
                            b_absorb_potential,b_num_abs_boundary_faces)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic,&
                                                 potential_dot_acoustic
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhostore,kappastore
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! absorbing boundary surface  
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces) 
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces) 

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,myrank,NGLOB_ADJOINT  
  integer:: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT) :: b_potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_potential
  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw 
  integer :: ispec,iglob,i,j,k,iface,igll
  !adjoint locals
  integer:: reclen1,reclen2

! adjoint simulations: 
  if (SIMULATION_TYPE == 3 .and. num_abs_boundary_faces > 0)  then
    read(IOABS_AC,rec=NSTEP-it+1) reclen1,b_absorb_potential,reclen2
    if (reclen1 /= b_reclen_potential .or. reclen1 /= reclen2) &
      call exit_mpi(myrank,'Error reading absorbing contribution b_absorb_potential')
  endif !adjoint
  
! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
      if( ispec_is_acoustic(ispec) ) then

        ! reference gll points on boundary face 
        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets global index
          iglob=ibool(i,j,k,ispec)

          ! determines bulk sound speed
          rhol = rhostore(i,j,k,ispec)
          cpl = sqrt( kappastore(i,j,k,ispec) / rhol )
             
          ! gets associated, weighted jacobian 
          jacobianw = abs_boundary_jacobian2Dw(igll,iface)
          
          ! Sommerfeld condition
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                              - potential_dot_acoustic(iglob) * jacobianw / cpl / rhol


          ! adjoint simulations          
          if (SIMULATION_TYPE == 3) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                                - b_absorb_potential(igll,iface)
          else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_potential(igll,iface) = potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
          endif !adjoint          
          
         enddo

      endif ! ispec_is_acoustic
    endif ! ispec_is_inner
  enddo ! num_abs_boundary_faces

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. num_abs_boundary_faces > 0 ) &
    write(IOABS_AC,rec=it) b_reclen_potential,b_absorb_potential,b_reclen_potential
  
  end subroutine compute_stacey_acoustic
