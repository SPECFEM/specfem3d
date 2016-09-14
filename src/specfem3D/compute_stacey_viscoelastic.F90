!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,iphase, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_9.F90"
#endif

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

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
  integer:: it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_1.F90"
#endif

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        vx = veloc(1,iglob)
        vy = veloc(2,iglob)
        vz = veloc(3,iglob)

#ifdef DEBUG_COUPLED
  include "../../../add_to_compute_stacey_viscoelastic_2.F90"
#endif

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

#ifdef DEBUG_COUPLED
  include "../../../add_to_compute_stacey_viscoelastic_3.F90"
#endif

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! adds stacey term (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! adjoint simulations
        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          b_absorb_field(1,igll,iface) = tx*jacobianw
          b_absorb_field(2,igll,iface) = ty*jacobianw
          b_absorb_field(3,igll,iface) = tz*jacobianw
        endif !adjoint

#ifdef DEBUG_COUPLED
  include "../../../add_to_compute_stacey_viscoelastic_4.F90"
#endif

      enddo
    endif ! ispec_is_elastic
  enddo

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_5.F90"
#endif

  end subroutine compute_stacey_viscoelastic

!
!=====================================================================
!

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_backward(NSPEC_AB, &
                        ibool,iphase, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        ispec_is_elastic,SIMULATION_TYPE, &
                        NSTEP,it,NGLOB_ADJOINT,b_accel, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants
  use specfem_par, only: myrank

  implicit none

  integer :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,NGLOB_ADJOINT
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel

! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  ! reads in absorbing boundary array (when first phase is running)
  ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
  call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        ! adjoint simulations
        b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_6.F90"
#endif

!
!=====================================================================
!

! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                        Mesh_pointer)

  use constants

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_10.F90"
#endif

  implicit none

! communication overlap
  integer :: iphase

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_7.F90"
#endif

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array (when first phase is running)
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
  endif !adjoint

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_field)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_8.F90"
#endif

  end subroutine compute_stacey_viscoelastic_GPU

#ifdef DEBUG_COUPLED
    include "../../../add_to_compute_stacey_viscoelastic_11.F90"
#endif
