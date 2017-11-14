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

subroutine compute_interpolated_dva_viscoelast(displ,veloc,accel,NGLOB_AB, &
                                    ispec,NSPEC_AB,ibool, &
                                    hxir,hetar,hgammar, &
                                    dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

  use constants

  implicit none

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd

  integer,intent(in) :: ispec

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB),intent(in) :: displ,veloc,accel
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! receiver Lagrange interpolators
  double precision,dimension(NGLLX),intent(in) :: hxir
  double precision,dimension(NGLLY),intent(in) :: hetar
  double precision,dimension(NGLLZ),intent(in) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,iglob

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO
  vxd = ZERO
  vyd = ZERO
  vzd = ZERO
  axd = ZERO
  ayd = ZERO
  azd = ZERO

  ! interpolates seismograms at exact receiver locations
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)

        hlagrange = hxir(i) * hetar(j) * hgammar(k)

        ! displacement
        dxd = dxd + dble(displ(1,iglob))*hlagrange
        dyd = dyd + dble(displ(2,iglob))*hlagrange
        dzd = dzd + dble(displ(3,iglob))*hlagrange
        ! velocity
        vxd = vxd + dble(veloc(1,iglob))*hlagrange
        vyd = vyd + dble(veloc(2,iglob))*hlagrange
        vzd = vzd + dble(veloc(3,iglob))*hlagrange
        ! acceleration
        axd = axd + dble(accel(1,iglob))*hlagrange
        ayd = ayd + dble(accel(2,iglob))*hlagrange
        azd = azd + dble(accel(3,iglob))*hlagrange

      enddo
    enddo
  enddo

end subroutine compute_interpolated_dva_viscoelast

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                        potential_dot_dot_acoustic,potential_acoustic,NGLOB_AB, &
                        ispec,NSPEC_AB,ibool, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd,USE_TRICK_FOR_BETTER_PRESSURE)

! for acoustic elements
! returns displacement/velocity/acceleration/pressure (dxd,..,vxd,..,axd,..,pd) at receiver location

  use constants

  implicit none

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd

  integer,intent(in) :: ispec

  integer,intent(in) :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(in) :: displ_element,veloc_element,accel_element
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB),intent(in) :: potential_dot_dot_acoustic,potential_acoustic

  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  logical,intent(in) :: USE_TRICK_FOR_BETTER_PRESSURE

  ! Lagrange interpolators
  double precision,dimension(NGLLX),intent(in) :: hxir
  double precision,dimension(NGLLY),intent(in) :: hetar
  double precision,dimension(NGLLZ),intent(in) :: hgammar

! local parameters
  double precision :: hlagrange
  integer :: i,j,k,iglob

! perform the general interpolation using Lagrange polynomials
  dxd = ZERO
  dyd = ZERO
  dzd = ZERO

  vxd = ZERO
  vyd = ZERO
  vzd = ZERO

  axd = ZERO
  ayd = ZERO
  azd = ZERO

  pd  = ZERO

  ! interpolates seismograms at exact receiver locations
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange = hxir(i) * hetar(j) * hgammar(k)

        ! displacement
        dxd = dxd + hlagrange*displ_element(1,i,j,k)
        dyd = dyd + hlagrange*displ_element(2,i,j,k)
        dzd = dzd + hlagrange*displ_element(3,i,j,k)

        ! velocity
        vxd = vxd + hlagrange*veloc_element(1,i,j,k)
        vyd = vyd + hlagrange*veloc_element(2,i,j,k)
        vzd = vzd + hlagrange*veloc_element(3,i,j,k)

        ! acceleration
        axd = axd + hlagrange*accel_element(1,i,j,k)
        ayd = ayd + hlagrange*accel_element(2,i,j,k)
        azd = azd + hlagrange*accel_element(3,i,j,k)

        ! global index
        iglob = ibool(i,j,k,ispec)

        ! pressure
        if (USE_TRICK_FOR_BETTER_PRESSURE) then
          ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
          ! use the second derivative of the source for the source time function instead of the source itself,
          ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
          ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
          ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
          ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
          ! is accurate at second order and thus contains significantly less numerical noise.
          pd = pd - hlagrange*potential_acoustic(iglob)
          ! that trick is not implemented for the calculation of displacement, velocity nor acceleration seismograms
          ! in acoustic elements yet; to do so we would need to recompute them using the second integral in time of the
          ! current formulas in that case. Same remark for recording stations located in solid (elastic/viscoelastic) elements
          ! in the case of fluid/solid models when that trick is used; thus for now we erase these seismograms here just in case
          ! because they would be wrong
          dxd = ZERO
          dyd = ZERO
          dzd = ZERO
          vxd = ZERO
          vyd = ZERO
          vzd = ZERO
          axd = ZERO
          ayd = ZERO
          azd = ZERO
        else
          pd = pd - hlagrange*potential_dot_dot_acoustic(iglob)
        endif

      enddo
    enddo
  enddo

end subroutine compute_interpolated_dva_acoust

