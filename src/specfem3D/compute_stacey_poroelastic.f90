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

! for poroelastic solver

! absorbing boundary terms for poroelastic media (type Stacey conditions)

  subroutine compute_stacey_poroelastic(NSPEC_AB,NGLOB_AB,accels,accelw, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        velocs,velocw,rho_vpI,rho_vpII,rho_vsI, &
                        rhoarraystore,phistore,tortstore, &
                        ispec_is_poroelastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        NSTEP,it,NGLOB_ADJOINT,b_accels,b_accelw, &
                        b_num_abs_boundary_faces,b_reclen_field_poro,b_absorb_fields, &
                        b_absorb_fieldw)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accels,accelw
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! type Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: velocs,velocw
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vpI,rho_vpII,rho_vsI
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: phistore,tortstore

  logical, dimension(NSPEC_AB) :: ispec_is_poroelastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,NGLOB_ADJOINT
  integer:: b_num_abs_boundary_faces,b_reclen_field_poro
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_fields, &
                                                                               b_absorb_fieldw
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accels,b_accelw
  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,vfx,vfy,vfz,nx,ny,nz,tx,ty,tz,tfx,tfy,tfz,vn,vfn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,rhol_bar,phil,tortl
  !integer:: reclen1,reclen2

  ! checks if anything to do
  if( num_abs_boundary_faces == 0 ) return

! adjoint simulations:
  if( SIMULATION_TYPE == 3 ) then
    ! reads in absorbing boundary array when first phase is running
    if( phase_is_inner .eqv. .false. ) then
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      ! uses fortran routine
      !read(IOABS,rec=NSTEP-it+1) reclen1,b_absorb_field,reclen2
      !if (reclen1 /= b_reclen_field .or. reclen1 /= reclen2) &
      !  call exit_mpi(0,'Error reading absorbing contribution b_absorb_field')
      ! uses c routine for faster reading
      call read_abs(0,b_absorb_fields,b_reclen_field_poro,NSTEP-it+1)
      call read_abs(0,b_absorb_fieldw,b_reclen_field_poro,NSTEP-it+1)
    endif
  endif !adjoint


     ! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
     do iface=1,num_abs_boundary_faces

        ispec = abs_boundary_ispec(iface)

        if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

           if( ispec_is_poroelastic(ispec) ) then

              ! reference gll points on boundary face
              do igll = 1,NGLLSQUARE

                 ! gets local indices for GLL point
                 i = abs_boundary_ijk(1,igll,iface)
                 j = abs_boundary_ijk(2,igll,iface)
                 k = abs_boundary_ijk(3,igll,iface)

                 ! gets velocity
                 iglob=ibool(i,j,k,ispec)
                 vx=velocs(1,iglob)
                 vy=velocs(2,iglob)
                 vz=velocs(3,iglob)
                 !
                 vfx=velocw(1,iglob)
                 vfy=velocw(2,iglob)
                 vfz=velocw(3,iglob)

                 ! gets properties
                 phil = phistore(i,j,k,ispec)
                 tortl = tortstore(i,j,k,ispec)
                 rhol_s = rhoarraystore(1,i,j,k,ispec)
                 rhol_f = rhoarraystore(2,i,j,k,ispec)
                 rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

                 ! gets associated normal
                 nx = abs_boundary_normal(1,igll,iface)
                 ny = abs_boundary_normal(2,igll,iface)
                 nz = abs_boundary_normal(3,igll,iface)

                 ! velocity component in normal direction (normal points out of element)
                 vn = vx*nx + vy*ny + vz*nz
                 vfn = vfx*nx + vfy*ny + vfz*nz

                 ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
                 tx = rho_vpI(i,j,k,ispec)*vn*nx + rho_vsI(i,j,k,ispec)*(vx-vn*nx)
                 ty = rho_vpI(i,j,k,ispec)*vn*ny + rho_vsI(i,j,k,ispec)*(vy-vn*ny)
                 tz = rho_vpI(i,j,k,ispec)*vn*nz + rho_vsI(i,j,k,ispec)*(vz-vn*nz)
                 !
                 tfx = rho_vpII(i,j,k,ispec)/(phil*rhol_bar)*tortl*rhol_f*vfn*nx - &
                       rho_vsI(i,j,k,ispec)/rhol_bar*rhol_f*(vx-vn*nx)
                 tfy = rho_vpII(i,j,k,ispec)/(phil*rhol_bar)*tortl*rhol_f*vfn*ny - &
                       rho_vsI(i,j,k,ispec)/rhol_bar*rhol_f*(vy-vn*ny)
                 tfz = rho_vpII(i,j,k,ispec)/(phil*rhol_bar)*tortl*rhol_f*vfn*nz - &
                       rho_vsI(i,j,k,ispec)/rhol_bar*rhol_f*(vz-vn*nz)

                 ! gets associated, weighted jacobian
                 jacobianw = abs_boundary_jacobian2Dw(igll,iface)

                 ! adds stacey term (weak form)
                 accels(1,iglob) = accels(1,iglob) - tx*jacobianw
                 accels(2,iglob) = accels(2,iglob) - ty*jacobianw
                 accels(3,iglob) = accels(3,iglob) - tz*jacobianw
                 !
                 accelw(1,iglob) = accelw(1,iglob) - tfx*jacobianw
                 accelw(2,iglob) = accelw(2,iglob) - tfy*jacobianw
                 accelw(3,iglob) = accelw(3,iglob) - tfz*jacobianw

                 ! adjoint simulations
                 if (SIMULATION_TYPE == 3) then
                    b_accels(:,iglob) = b_accels(:,iglob) - b_absorb_fields(:,igll,iface)
                    !
                    b_accelw(:,iglob) = b_accelw(:,iglob) - b_absorb_fieldw(:,igll,iface)
                 else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
                    b_absorb_fields(1,igll,iface) = tx*jacobianw
                    b_absorb_fields(2,igll,iface) = ty*jacobianw
                    b_absorb_fields(3,igll,iface) = tz*jacobianw
                    !
                    b_absorb_fieldw(1,igll,iface) = tfx*jacobianw
                    b_absorb_fieldw(2,igll,iface) = tfy*jacobianw
                    b_absorb_fieldw(3,igll,iface) = tfz*jacobianw
                 endif !adjoint

              enddo
           endif ! ispec_is_poroelastic
        endif ! ispec_is_inner
     enddo

  ! adjoint simulations: stores absorbed wavefield part
  if( SIMULATION_TYPE == 1 .and. SAVE_FORWARD ) then
    ! writes out absorbing boundary value only when second phase is running
    if( phase_is_inner .eqv. .true. ) then
      ! uses fortran routine
      !write(IOABS,rec=it) b_reclen_field,b_absorb_field,b_reclen_field
      ! uses c routine
      call write_abs(0,b_absorb_fields,b_reclen_field_poro,it)
      call write_abs(0,b_absorb_fieldw,b_reclen_field_poro,it)
    endif
  endif

  end subroutine compute_stacey_poroelastic

