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


subroutine compute_interpolated_dva(displ,veloc,accel,NGLOB_AB, &
                        ispec,NSPEC_AB,ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

  implicit none
  include 'constants.h'

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_AB) :: displ,veloc,accel
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

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

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! displacement
    dxd = dble(displ(1,iglob))
    dyd = dble(displ(2,iglob))
    dzd = dble(displ(3,iglob))
    ! velocity
    vxd = dble(veloc(1,iglob))
    vyd = dble(veloc(2,iglob))
    vzd = dble(veloc(3,iglob))
    ! acceleration
    axd = dble(accel(1,iglob))
    ayd = dble(accel(2,iglob))
    azd = dble(accel(3,iglob))

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

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

  endif

end subroutine compute_interpolated_dva

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_interpolated_dva_ac(displ_element,veloc_element,&
                        potential_dot_dot_acoustic,potential_dot_acoustic,&
                        potential_acoustic,NGLOB_AB, &
                        ispec,NSPEC_AB,ibool, &
                        xi_r,eta_r,gamma_r, &
                        hxir,hetar,hgammar, &
                        dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd)

! acoustic elements
! returns displacement/velocity/pressure (dxd,..,vxd,..,axd,.. ) at receiver location

  implicit none
  include 'constants.h'

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd

  integer :: ispec

  integer :: NSPEC_AB,NGLOB_AB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: displ_element,veloc_element
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLOB_AB) :: potential_acoustic

  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB):: ibool

  ! receiver information
  double precision :: xi_r,eta_r,gamma_r
  double precision,dimension(NGLLX) :: hxir
  double precision,dimension(NGLLY) :: hetar
  double precision,dimension(NGLLZ) :: hgammar

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

! takes closest GLL point only (no interpolation)
  if(FASTER_RECEIVERS_POINTS_ONLY) then

    ! displacement
    dxd = displ_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    dyd = displ_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    dzd = displ_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))
    ! velocity
    vxd = veloc_element(1,nint(xi_r),nint(eta_r),nint(gamma_r))
    vyd = veloc_element(2,nint(xi_r),nint(eta_r),nint(gamma_r))
    vzd = veloc_element(3,nint(xi_r),nint(eta_r),nint(gamma_r))

    ! global index
    iglob = ibool(nint(xi_r),nint(eta_r),nint(gamma_r),ispec)

    ! x component -> acoustic potential
    axd = potential_acoustic(iglob)
    ! y component -> first time derivative of potential
    ayd = potential_dot_acoustic(iglob)
    ! z component -> pressure
    azd = - potential_dot_dot_acoustic(iglob)

  else

! interpolates seismograms at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          hlagrange = hxir(i)*hetar(j)*hgammar(k)

          ! displacement
          dxd = dxd + hlagrange*displ_element(1,i,j,k)
          dyd = dyd + hlagrange*displ_element(2,i,j,k)
          dzd = dzd + hlagrange*displ_element(3,i,j,k)
          ! velocity
          vxd = vxd + hlagrange*veloc_element(1,i,j,k)
          vyd = vyd + hlagrange*veloc_element(2,i,j,k)
          vzd = vzd + hlagrange*veloc_element(3,i,j,k)

          ! x component -> acoustic potential
          axd = axd + hlagrange*potential_acoustic(iglob)
          ! y component -> first time derivative of potential
          ayd = ayd + hlagrange*potential_dot_acoustic(iglob)
          ! z component -> pressure
          azd = azd - hlagrange*potential_dot_dot_acoustic(iglob)

        enddo
      enddo
    enddo

  endif

end subroutine compute_interpolated_dva_ac
