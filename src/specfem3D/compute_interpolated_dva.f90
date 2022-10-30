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

  subroutine compute_interpolated_dva_viscoelast(displ,veloc,accel,NGLOB_AB, &
                                                 ispec,NSPEC_AB,ibool, &
                                                 hxir,hetar,hgammar, &
                                                 dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd)

! returns displacement/velocity/acceleration (dxd,..,vxd,..,axd,.. ) at receiver location

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ZERO,FOUR_THIRDS

  use specfem_par, only: SAVE_SEISMOGRAMS_PRESSURE,ANISOTROPY
  use specfem_par, only: kappastore,mustore

  use specfem_par_elastic, only: c11store,c12store,c13store,c14store,c15store, &
                                 c16store,c22store,c23store,c24store,c25store, &
                                 c26store,c33store,c34store,c35store,c36store

  implicit none

  double precision,intent(out) :: dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd

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
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_trace_over_3_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: pressure
  real(kind=CUSTOM_REAL) :: mul,kappal,lambdal,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz

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
  pd = ZERO

  ! interpolates seismograms at exact receiver locations
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec)

        hlagrange = hxir(i) * hetar(j) * hgammar(k)

        ! displacement
        dxd = dxd + dble(displ(1,iglob)) * hlagrange
        dyd = dyd + dble(displ(2,iglob)) * hlagrange
        dzd = dzd + dble(displ(3,iglob)) * hlagrange
        ! velocity
        vxd = vxd + dble(veloc(1,iglob)) * hlagrange
        vyd = vyd + dble(veloc(2,iglob)) * hlagrange
        vzd = vzd + dble(veloc(3,iglob)) * hlagrange
        ! acceleration
        axd = axd + dble(accel(1,iglob)) * hlagrange
        ayd = ayd + dble(accel(2,iglob)) * hlagrange
        azd = azd + dble(accel(3,iglob)) * hlagrange

      enddo
    enddo
  enddo

  if (SAVE_SEISMOGRAMS_PRESSURE) then
    ! for an elastic element:
    !
    ! isostatic stress or pressure: p = - 1/3 trace(sigma)
    ! (corresponds to hydrostatic pressure in fluids)
    !
    ! from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
    ! J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
    ! for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
    ! p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.
    !
    ! to compute pressure in 3D in an elastic solid, one uses: pressure = - trace(sigma) / 3
    !
    ! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
    !          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
    !
    ! sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
    ! sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
    ! sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
    !
    ! pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
    !
    ! this routines limits the pressure computations to: non-anisotropic, non-attenuation case
    ! todo for the future...

    ! strain for ispec element
    call compute_element_strain(ispec,NGLOB_AB,displ,epsilondev_loc,epsilondev_trace_over_3_loc)

    ! pressure within element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! single derivatives from strain elements
          !eps_trace_over_3_loc(i,j,k) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          !epsilondev_loc(1,i,j,k) = duxdxl - eps_trace_over_3_loc(i,j,k)
          !epsilondev_loc(2,i,j,k) = duydyl - eps_trace_over_3_loc(i,j,k)
          duxdxl = epsilondev_loc(1,i,j,k) + epsilondev_trace_over_3_loc(i,j,k)
          duydyl = epsilondev_loc(2,i,j,k) + epsilondev_trace_over_3_loc(i,j,k)
          duzdzl = 3.0_CUSTOM_REAL * epsilondev_trace_over_3_loc(i,j,k) - duxdxl - duydyl

          !epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          !epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          !epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
          duxdyl_plus_duydxl = 2.0_CUSTOM_REAL * epsilondev_loc(3,i,j,k)
          duzdxl_plus_duxdzl = 2.0_CUSTOM_REAL * epsilondev_loc(4,i,j,k)
          duzdyl_plus_duydzl = 2.0_CUSTOM_REAL * epsilondev_loc(5,i,j,k)

          ! computes either isotropic or anisotropic element stresses
          if (ANISOTROPY) then
            ! full anisotropic case, stress calculations
            c11 = c11store(i,j,k,ispec)
            c12 = c12store(i,j,k,ispec)
            c13 = c13store(i,j,k,ispec)
            c14 = c14store(i,j,k,ispec)
            c15 = c15store(i,j,k,ispec)
            c16 = c16store(i,j,k,ispec)
            c22 = c22store(i,j,k,ispec)
            c23 = c23store(i,j,k,ispec)
            c24 = c24store(i,j,k,ispec)
            c25 = c25store(i,j,k,ispec)
            c26 = c26store(i,j,k,ispec)
            c33 = c33store(i,j,k,ispec)
            c34 = c34store(i,j,k,ispec)
            c35 = c35store(i,j,k,ispec)
            c36 = c36store(i,j,k,ispec)

            sigma_xx = c11 * duxdxl + c16 * duxdyl_plus_duydxl + c12 * duydyl + &
                       c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * duzdzl
            sigma_yy = c12 * duxdxl + c26 * duxdyl_plus_duydxl + c22 * duydyl + &
                       c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * duzdzl
            sigma_zz = c13 * duxdxl + c36 * duxdyl_plus_duydxl + c23 * duydyl + &
                       c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * duzdzl
          else
            ! isotropic case
            kappal = kappastore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)

            lambdalplus2mul = kappal + FOUR_THIRDS * mul
            lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

            ! compute stress sigma
            sigma_xx = lambdalplus2mul * duxdxl + lambdal * (duydyl + duzdzl)
            sigma_yy = lambdalplus2mul * duydyl + lambdal * (duxdxl + duzdzl)
            sigma_zz = lambdalplus2mul * duzdzl + lambdal * (duxdxl + duydyl)
          endif ! ANISOTROPY

          ! pressure p = - 1/3 trace(sigma)
          pressure(i,j,k) = - (sigma_xx + sigma_yy + sigma_zz) / 3.0_CUSTOM_REAL
        enddo
      enddo
    enddo

    ! interpolates pressure at exact receiver locations
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          hlagrange = hxir(i) * hetar(j) * hgammar(k)
          ! pressure
          pd = pd + dble(pressure(i,j,k)) * hlagrange
        enddo
      enddo
    enddo

  endif

  end subroutine compute_interpolated_dva_viscoelast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_interpolated_dva_acoust(displ_element,veloc_element,accel_element, &
                                             potential_dot_dot_acoustic,potential_acoustic,NGLOB_AB, &
                                             ispec,NSPEC_AB,ibool, &
                                             hxir,hetar,hgammar, &
                                             dxd,dyd,dzd,vxd,vyd,vzd,axd,ayd,azd,pd, &
                                             USE_TRICK_FOR_BETTER_PRESSURE)

! for acoustic elements
! returns displacement/velocity/acceleration/pressure (dxd,..,vxd,..,axd,..,pd) at receiver location

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,ZERO

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
        dxd = dxd + hlagrange * displ_element(1,i,j,k)
        dyd = dyd + hlagrange * displ_element(2,i,j,k)
        dzd = dzd + hlagrange * displ_element(3,i,j,k)

        ! velocity
        vxd = vxd + hlagrange * veloc_element(1,i,j,k)
        vyd = vyd + hlagrange * veloc_element(2,i,j,k)
        vzd = vzd + hlagrange * veloc_element(3,i,j,k)

        ! acceleration
        axd = axd + hlagrange * accel_element(1,i,j,k)
        ayd = ayd + hlagrange * accel_element(2,i,j,k)
        azd = azd + hlagrange * accel_element(3,i,j,k)

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
          pd = pd - hlagrange * potential_acoustic(iglob)
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
          pd = pd - hlagrange * potential_dot_dot_acoustic(iglob)
        endif

      enddo
    enddo
  enddo

  end subroutine compute_interpolated_dva_acoust

