!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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
                            ibool,iphase, &
                            abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                            SIMULATION_TYPE,SAVE_FORWARD,it,b_reclen_potential, &
                            b_absorb_potential,b_num_abs_boundary_faces)

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic, &
                                                 potential_dot_acoustic
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhostore,kappastore
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE,it
  integer:: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_potential
  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw,absorbl
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer:: reclen1,reclen2

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
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
        absorbl = potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - absorbl

        ! adjoint simulations
        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          b_absorb_potential(igll,iface) = absorbl
        endif !adjoint

       enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
  endif

  end subroutine compute_stacey_acoustic
!
!=====================================================================
! for acoustic solver for back propagation wave field

  subroutine compute_stacey_acoustic_backward(NSPEC_AB, &
                            ibool,iphase, &
                            abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces,ispec_is_acoustic, &
                            SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                            b_potential_dot_dot_acoustic,b_reclen_potential, &
                            b_absorb_potential,b_num_abs_boundary_faces)

  use constants

  implicit none

  integer :: NSPEC_AB

! potentials
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,NGLOB_ADJOINT
  integer:: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT) :: b_potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_potential

! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer:: reclen1,reclen2

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,NSTEP-it+1)
  endif !adjoint

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets global index
        iglob=ibool(i,j,k,ispec)

        ! adjoint simulations
        if (SIMULATION_TYPE == 3) then
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                              - b_absorb_potential(igll,iface)
        endif !adjoint

      enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  end subroutine compute_stacey_acoustic_backward
!
!=====================================================================
! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                            b_reclen_potential,b_absorb_potential, &
                            b_num_abs_boundary_faces,Mesh_pointer)

  use constants

  implicit none

! potentials

! communication overlap
  integer :: iphase

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  integer:: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_potential
  logical:: SAVE_FORWARD

  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array (when first phase is running)
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,NSTEP-it+1)
  endif !adjoint

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  call compute_stacey_acoustic_cuda(Mesh_pointer,iphase,b_absorb_potential)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
  endif

  end subroutine compute_stacey_acoustic_GPU
