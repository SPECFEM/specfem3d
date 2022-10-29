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


! for poroelastic solver

  subroutine compute_coupling_poroelastic_ac(NSPEC_AB,NGLOB_AB, &
                                             ibool,accels_poroelastic,accelw_poroelastic, &
                                             potential_dot_dot_acoustic, &
                                             num_coupling_ac_po_faces, &
                                             coupling_ac_po_ispec,coupling_ac_po_ijk, &
                                             coupling_ac_po_normal, &
                                             coupling_ac_po_jacobian2Dw, &
                                             rhoarraystore,phistore,tortstore, &
                                             iphase)

! returns the updated accelerations array: accels_poroelatsic & accelw_poroelastic

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accels_poroelastic,accelw_poroelastic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! acoustic-poroelastic coupling surface
  integer :: num_coupling_ac_po_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_po_normal(NDIM,NGLLSQUARE,num_coupling_ac_po_faces)
  real(kind=CUSTOM_REAL) :: coupling_ac_po_jacobian2Dw(NGLLSQUARE,num_coupling_ac_po_faces)
  integer :: coupling_ac_po_ijk(3,NGLLSQUARE,num_coupling_ac_po_faces)
  integer :: coupling_ac_po_ispec(num_coupling_ac_po_faces)

! properties
  real(kind=CUSTOM_REAL), dimension(2,NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhoarraystore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        phistore,tortstore

! communication overlap
  integer :: iphase

! local parameters
  real(kind=CUSTOM_REAL) :: pressure
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,phil,tortl,rhol_bar
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw,fac_s,fac_w

  integer :: iface,igll,ispec,iglob
  integer :: i,j,k

  ! only add these contributions in first pass
  if (iphase /= 1) return

! loops on all coupling faces
  do iface = 1,num_coupling_ac_po_faces

    ! gets corresponding spectral element, constructed such that it is a poroelastic
    ! element, since we need to have access to porous properties
    ispec = coupling_ac_po_ispec(iface)

    ! loops over common GLL points
    do igll = 1, NGLLSQUARE
      i = coupling_ac_po_ijk(1,igll,iface)
      j = coupling_ac_po_ijk(2,igll,iface)
      k = coupling_ac_po_ijk(3,igll,iface)

      ! gets global index of this common GLL point
      ! (note: should be the same as for corresponding i',j',k',ispec_poroelastic or ispec_acoustic )
      iglob = ibool(i,j,k,ispec)

      ! get poroelastic parameters
      phil = phistore(i,j,k,ispec)
      tortl = tortstore(i,j,k,ispec)
      rhol_s = rhoarraystore(1,i,j,k,ispec)
      rhol_f = rhoarraystore(2,i,j,k,ispec)
      rhol_bar = (1.0_CUSTOM_REAL - phil) * rhol_s + phil * rhol_f

      fac_s = 1.0_CUSTOM_REAL - phil/tortl
      fac_w = 1.0_CUSTOM_REAL - rhol_f/rhol_bar

      ! acoustic pressure on global point
      pressure = - potential_dot_dot_acoustic(iglob)

      ! gets associated normal on GLL point
      ! (note convention: pointing outwards of acoustic element)
      nx = coupling_ac_po_normal(1,igll,iface)
      ny = coupling_ac_po_normal(2,igll,iface)
      nz = coupling_ac_po_normal(3,igll,iface)

      ! gets associated, weighted 2D jacobian
      ! (note: should be the same for poroelastic and acoustic element)
      jacobianw = coupling_ac_po_jacobian2Dw(igll,iface)

      ! continuity of displacement and pressure on global point
      !
      ! note: Newmark time scheme together with definition of scalar potential:
      !          pressure = - chi_dot_dot
      !          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
      !          pressure at time step [t + delta_t]
      !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
      !          it means you have to calculate and update the acoustic pressure first before
      !          calculating this term...

      ! contribution to the solid phase
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + jacobianw * nx * pressure * fac_s
      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + jacobianw * ny * pressure * fac_s
      accels_poroelastic(3,iglob) = accels_poroelastic(3,iglob) + jacobianw * nz * pressure * fac_s
      ! contribution to the fluid phase
      accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + jacobianw * nx * pressure * fac_w
      accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + jacobianw * ny * pressure * fac_w
      accelw_poroelastic(3,iglob) = accelw_poroelastic(3,iglob) + jacobianw * nz * pressure * fac_w

    enddo ! igll

  enddo ! iface

  end subroutine compute_coupling_poroelastic_ac

