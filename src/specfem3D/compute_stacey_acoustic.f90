!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

  subroutine compute_stacey_acoustic_forward(NSPEC_AB,NGLOB_AB, &
                                             potential_dot_dot_acoustic,potential_dot_acoustic, &
                                             ibool,iphase, &
                                             abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                             num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                                             it,b_reclen_potential, &
                                             b_absorb_potential,b_num_abs_boundary_faces)

  use constants
  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! communication overlap
  integer,intent(in) :: iphase

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rhostore,kappastore
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer,intent(in) :: it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_potential

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
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_absorb_potential(igll,iface) = absorbl
        endif !adjoint

       enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  ! for kernel simulations: stores absorbed wavefield part
  if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
    ! writes out absorbing boundary value
    call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
  endif

  end subroutine compute_stacey_acoustic_forward
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

  integer,intent(in) :: NSPEC_AB

! potentials
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! communication overlap
  integer,intent(in) :: iphase

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer,intent(in) :: SIMULATION_TYPE
  integer,intent(in) :: NSTEP,it,NGLOB_ADJOINT
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT),intent(inout) :: b_potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces),intent(in) :: b_absorb_potential

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
!

  subroutine compute_stacey_acoustic_backward_undoatt(NSPEC_AB,NGLOB_AB, &
                                                      b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                      ibool,iphase, &
                                                      abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                                      num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic)

  use constants

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: b_potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: b_potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! communication overlap
  integer,intent(in) :: iphase

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rhostore,kappastore
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

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
        absorbl = b_potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
        b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - absorbl

       enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  end subroutine compute_stacey_acoustic_backward_undoatt

!
!=====================================================================
! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                                         NSTEP,it, &
                                         b_reclen_potential,b_absorb_potential, &
                                         b_num_abs_boundary_faces,Mesh_pointer, &
                                         FORWARD_OR_ADJOINT)

  use constants
  use specfem_par, only: SIMULATION_TYPE,SAVE_FORWARD,UNDO_ATTENUATION_AND_OR_PML

  implicit none

! potentials

! communication overlap
  integer,intent(in) :: iphase

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces

! adjoint simulations
  integer, intent(in) :: NSTEP,it
  integer, intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces), intent(inout) :: b_absorb_potential

  ! GPU_MODE variables
  integer(kind=8), intent(in) :: Mesh_pointer
  integer, intent(in) :: FORWARD_OR_ADJOINT

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! no need to store boundaries on disk
    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_acoustic_undoatt_cuda(Mesh_pointer,iphase,FORWARD_OR_ADJOINT)
  else
    ! adjoint simulations:
    if (SIMULATION_TYPE == 3) then
      ! reads in absorbing boundary array (when first phase is running)
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      call read_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,NSTEP-it+1)
    endif !adjoint

    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_acoustic_cuda(Mesh_pointer,iphase,b_absorb_potential,FORWARD_OR_ADJOINT)

    ! for kernel simulations: stores absorbed wavefield part
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! writes out absorbing boundary value
      call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
    endif
  endif

  end subroutine compute_stacey_acoustic_GPU
