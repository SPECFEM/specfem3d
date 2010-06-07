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

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_elastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,myrank,SAVE_FORWARD, &
                        NSTEP,it,NGLOB_ADJOINT,b_accel, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface  
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces) 
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces) 
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces) 

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,myrank,NGLOB_ADJOINT  
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel
  logical:: SAVE_FORWARD
  
! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  !adjoint locals
  integer:: reclen1,reclen2
  
  
! adjoint simulations: 
  if (SIMULATION_TYPE == 3 .and. num_abs_boundary_faces > 0)  then
    read(IOABS,rec=NSTEP-it+1) reclen1,b_absorb_field,reclen2
    if (reclen1 /= b_reclen_field .or. reclen1 /= reclen2) &
      call exit_mpi(myrank,'Error reading absorbing contribution b_absorb_field')
  endif !adjoint
  

! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      if( ispec_is_elastic(ispec) ) then
      
        ! reference gll points on boundary face 
        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets velocity
          iglob=ibool(i,j,k,ispec)
          vx=veloc(1,iglob)
          vy=veloc(2,iglob)
          vz=veloc(3,iglob)

          ! gets associated normal
          nx = abs_boundary_normal(1,igll,iface)
          ny = abs_boundary_normal(2,igll,iface)
          nz = abs_boundary_normal(3,igll,iface)             

          ! velocity component in normal direction (normal points out of element)
          vn = vx*nx + vy*ny + vz*nz
             
          ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it 
          tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
          ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
          tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

          ! gets associated, weighted jacobian 
          jacobianw = abs_boundary_jacobian2Dw(igll,iface)
          
          ! adds stacey term (weak form)
          accel(1,iglob) = accel(1,iglob) - tx*jacobianw
          accel(2,iglob) = accel(2,iglob) - ty*jacobianw
          accel(3,iglob) = accel(3,iglob) - tz*jacobianw

          ! adjoint simulations
          if (SIMULATION_TYPE == 3) then
            b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
          else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
            b_absorb_field(1,igll,iface) = tx*jacobianw
            b_absorb_field(2,igll,iface) = ty*jacobianw
            b_absorb_field(3,igll,iface) = tz*jacobianw
          endif !adjoint

         enddo
      endif ! ispec_is_elastic
    endif ! ispec_is_inner    
  enddo

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. num_abs_boundary_faces > 0 ) &
    write(IOABS,rec=it) b_reclen_field,b_absorb_field,b_reclen_field
  
  end subroutine compute_stacey_elastic

