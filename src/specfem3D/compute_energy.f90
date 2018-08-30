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

  subroutine compute_energy()

! computes kinetic, potential and total energy
! in all the slices using an MPI reduction
! and output that to an energy file

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use pml_par

  implicit none

! local variables
  integer :: i,j,k,l,ispec,iglob,ispec_irreg

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  double precision :: epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz
  double precision :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz
  double precision :: vx,vy,vz,pressure

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul,rhol,rho_invl
  real(kind=CUSTOM_REAL) :: kappal

  double precision :: integration_weight
  double precision :: kinetic_energy,potential_energy
  double precision :: kinetic_energy_glob,potential_energy_glob,total_energy_glob

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer :: i_SLS

  ! local anisotropy parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local attenuation parameters
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val

! that trick (USE_TRICK_FOR_BETTER_PRESSURE) is not implemented for the calculation
! of displacement, velocity nor acceleration vectors in acoustic elements yet;
! to do so we would need to recompute them using the second integral in time of the
! current formulas in that case. And to compute kinetic energy, we need the velocity vector...
  if (USE_TRICK_FOR_BETTER_PRESSURE) &
    call exit_mpi(myrank,'USE_TRICK_FOR_BETTER_PRESSURE not implemented for OUTPUT_ENERGY, please turn one of them off')

  kinetic_energy = 0.d0
  potential_energy = 0.d0

! loop over spectral elements
  do ispec = 1,NSPEC_AB

! if element is a CPML then do not compute energy in it, since it is non physical;
! thus, we compute energy in the main domain only, without absorbing elements
    if (is_CPML(ispec)) cycle

    ispec_irreg = irregular_element_number(ispec)

    !---
    !--- elastic spectral element
    !---
    if (ispec_is_elastic(ispec)) then

      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob)
          enddo
        enddo
      enddo

      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            iglob = ibool(i,j,k,ispec)

            tempx1(i,j,k) = 0._CUSTOM_REAL
            tempx2(i,j,k) = 0._CUSTOM_REAL
            tempx3(i,j,k) = 0._CUSTOM_REAL

            tempy1(i,j,k) = 0._CUSTOM_REAL
            tempy2(i,j,k) = 0._CUSTOM_REAL
            tempy3(i,j,k) = 0._CUSTOM_REAL

            tempz1(i,j,k) = 0._CUSTOM_REAL
            tempz2(i,j,k) = 0._CUSTOM_REAL
            tempz3(i,j,k) = 0._CUSTOM_REAL

            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              hp1 = hprime_xx(i,l)
              tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1
              tempy1(i,j,k) = tempy1(i,j,k) + dummyy_loc(l,j,k)*hp1
              tempz1(i,j,k) = tempz1(i,j,k) + dummyz_loc(l,j,k)*hp1

              hp2 = hprime_yy(j,l)
              tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2
              tempy2(i,j,k) = tempy2(i,j,k) + dummyy_loc(i,l,k)*hp2
              tempz2(i,j,k) = tempz2(i,j,k) + dummyz_loc(i,l,k)*hp2

              hp3 = hprime_zz(k,l)
              tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
              tempy3(i,j,k) = tempy3(i,j,k) + dummyy_loc(i,j,l)*hp3
              tempz3(i,j,k) = tempz3(i,j,k) + dummyz_loc(i,j,l)*hp3
            enddo

            if (ispec_irreg /= 0 ) then !irregular element

              ! get derivatives of ux, uy and uz with respect to x, y and z
              xixl = xix(i,j,k,ispec_irreg)
              xiyl = xiy(i,j,k,ispec_irreg)
              xizl = xiz(i,j,k,ispec_irreg)
              etaxl = etax(i,j,k,ispec_irreg)
              etayl = etay(i,j,k,ispec_irreg)
              etazl = etaz(i,j,k,ispec_irreg)
              gammaxl = gammax(i,j,k,ispec_irreg)
              gammayl = gammay(i,j,k,ispec_irreg)
              gammazl = gammaz(i,j,k,ispec_irreg)

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

            else ! regular element

              duxdxl = xix_regular*tempx1(i,j,k)
              duxdyl = xix_regular*tempx2(i,j,k)
              duxdzl = xix_regular*tempx3(i,j,k)

              duydxl = xix_regular*tempy1(i,j,k)
              duydyl = xix_regular*tempy2(i,j,k)
              duydzl = xix_regular*tempy3(i,j,k)

              duzdxl = xix_regular*tempz1(i,j,k)
              duzdyl = xix_regular*tempz2(i,j,k)
              duzdzl = xix_regular*tempz3(i,j,k)

              jacobianl = jacobian_regular

            endif

            ! precompute some sums to save CPU time
            duxdxl_plus_duydyl = duxdxl + duydyl
            duxdxl_plus_duzdzl = duxdxl + duzdzl
            duydyl_plus_duzdzl = duydyl + duzdzl
            duxdyl_plus_duydxl = duxdyl + duydxl
            duzdxl_plus_duxdzl = duzdxl + duxdzl
            duzdyl_plus_duydzl = duzdyl + duydzl

            ! compute the strain
            epsilon_xx = duxdxl
            epsilon_yy = duydyl
            epsilon_zz = duzdzl
            epsilon_xy = 0.5d0 * duxdyl_plus_duydxl
            epsilon_xz = 0.5d0 * duzdxl_plus_duxdzl
            epsilon_yz = 0.5d0 * duzdyl_plus_duydzl

            kappal = kappastore(i,j,k,ispec)
            mul = mustore(i,j,k,ispec)
            rhol = rhostore(i,j,k,ispec)

            ! full anisotropic case, stress calculations
            if (ANISOTROPY) then
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
              c44 = c44store(i,j,k,ispec)
              c45 = c45store(i,j,k,ispec)
              c46 = c46store(i,j,k,ispec)
              c55 = c55store(i,j,k,ispec)
              c56 = c56store(i,j,k,ispec)
              c66 = c66store(i,j,k,ispec)

              sigma_xx = c11 * duxdxl + c16 * duxdyl_plus_duydxl + c12 * duydyl + &
                         c15 * duzdxl_plus_duxdzl + c14 * duzdyl_plus_duydzl + c13 * duzdzl
              sigma_yy = c12 * duxdxl + c26 * duxdyl_plus_duydxl + c22 * duydyl + &
                         c25 * duzdxl_plus_duxdzl + c24 * duzdyl_plus_duydzl + c23 * duzdzl
              sigma_zz = c13 * duxdxl + c36 * duxdyl_plus_duydxl + c23 * duydyl + &
                         c35 * duzdxl_plus_duxdzl + c34 * duzdyl_plus_duydzl + c33 * duzdzl
              sigma_xy = c16 * duxdxl + c66 * duxdyl_plus_duydxl + c26 * duydyl + &
                         c56 * duzdxl_plus_duxdzl + c46 * duzdyl_plus_duydzl + c36 * duzdzl
              sigma_xz = c15 * duxdxl + c56 * duxdyl_plus_duydxl + c25 * duydyl + &
                         c55 * duzdxl_plus_duxdzl + c45 * duzdyl_plus_duydzl + c35 * duzdzl
              sigma_yz = c14 * duxdxl + c46 * duxdyl_plus_duydxl + c24 * duydyl + &
                         c45 * duzdxl_plus_duxdzl + c44 * duzdyl_plus_duydzl + c34 * duzdzl

            else

              ! isotropic case
              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

              ! compute stress sigma
              sigma_xx = lambdalplus2mul * duxdxl + lambdal * duydyl_plus_duzdzl
              sigma_yy = lambdalplus2mul * duydyl + lambdal * duxdxl_plus_duzdzl
              sigma_zz = lambdalplus2mul * duzdzl + lambdal * duxdxl_plus_duydyl

              sigma_xy = mul * duxdyl_plus_duydxl
              sigma_xz = mul * duzdxl_plus_duxdzl
              sigma_yz = mul * duzdyl_plus_duydzl

            endif ! ANISOTROPY

            ! subtract memory variables if attenuation
            if (ATTENUATION) then
              do i_sls = 1,N_SLS
                R_xx_val = R_xx(i_sls,i,j,k,ispec)
                R_yy_val = R_yy(i_sls,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val
                sigma_yy = sigma_yy - R_yy_val
                sigma_zz = sigma_zz + R_xx_val + R_yy_val
                sigma_xy = sigma_xy - R_xy(i_sls,i,j,k,ispec)
                sigma_xz = sigma_xz - R_xz(i_sls,i,j,k,ispec)
                sigma_yz = sigma_yz - R_yz(i_sls,i,j,k,ispec)
              enddo
            endif

            integration_weight = wxgll(i)*wygll(j)*wzgll(k)*jacobianl

            ! velocity
            vx = veloc(1,iglob)
            vy = veloc(2,iglob)
            vz = veloc(3,iglob)

            ! compute kinetic energy  1/2 rho ||v||^2
            ! we will divide the total sum by 2 only once, at the end of this routine, to reduce compute time
            kinetic_energy = kinetic_energy + integration_weight * rhol*(vx**2 + vy**2 + vz**2)

            ! compute potential energy 1/2 sigma_ij epsilon_ij
            ! we will divide the total sum by 2 only once, at the end of this routine, to reduce compute time
            potential_energy = potential_energy + integration_weight * &
              (sigma_xx*epsilon_xx + sigma_yy*epsilon_yy + sigma_zz*epsilon_zz + &
               2.d0 * (sigma_xy*epsilon_xy + sigma_xz*epsilon_xz + sigma_yz*epsilon_yz))

          enddo ! of the triple loop on i,j,k
        enddo
      enddo

    !---
    !--- acoustic spectral element
    !---
    else if (ispec_is_acoustic(ispec)) then

      ! for the definition of potential energy in an acoustic fluid, see for instance
      ! equation (23) of M. Maess et al., Journal of Sound and Vibration 296 (2006) 264-276

      ! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
      ! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
      ! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
      ! This permits acoustic-elastic coupling based on a non-iterative time scheme.
      ! Displacement is then: u = grad(Chi) / rho
      ! velocity is then: v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
      ! and pressure is: p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).

      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = potential_dot_acoustic(iglob)
          enddo
        enddo
      enddo

      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            iglob = ibool(i,j,k,ispec)

            tempx1(i,j,k) = 0._CUSTOM_REAL
            tempx2(i,j,k) = 0._CUSTOM_REAL
            tempx3(i,j,k) = 0._CUSTOM_REAL

            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              hp1 = hprime_xx(i,l)
              tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1

              hp2 = hprime_yy(j,l)
              tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2

              hp3 = hprime_zz(k,l)
              tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
            enddo

            if (ispec_irreg /= 0 ) then !irregular element

              ! get derivatives of ux, uy and uz with respect to x, y and z
              xixl = xix(i,j,k,ispec_irreg)
              xiyl = xiy(i,j,k,ispec_irreg)
              xizl = xiz(i,j,k,ispec_irreg)
              etaxl = etax(i,j,k,ispec_irreg)
              etayl = etay(i,j,k,ispec_irreg)
              etazl = etaz(i,j,k,ispec_irreg)
              gammaxl = gammax(i,j,k,ispec_irreg)
              gammayl = gammay(i,j,k,ispec_irreg)
              gammazl = gammaz(i,j,k,ispec_irreg)

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

            else ! regular element

              duxdxl = xix_regular*tempx1(i,j,k)
              duxdyl = xix_regular*tempx2(i,j,k)
              duxdzl = xix_regular*tempx3(i,j,k)

              duydxl = xix_regular*tempy1(i,j,k)
              duydyl = xix_regular*tempy2(i,j,k)
              duydzl = xix_regular*tempy3(i,j,k)

              duzdxl = xix_regular*tempz1(i,j,k)
              duzdyl = xix_regular*tempz2(i,j,k)
              duzdzl = xix_regular*tempz3(i,j,k)

              jacobianl = jacobian_regular

            endif

            rhol = rhostore(i,j,k,ispec)
            rho_invl = 1._CUSTOM_REAL / rhol
            kappal = kappastore(i,j,k,ispec)

            ! velocity is v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
            vx = duxdxl * rho_invl
            vy = duxdyl * rho_invl
            vz = duxdzl * rho_invl

            ! pressure is p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi)
            pressure = - potential_dot_dot_acoustic(iglob)

            integration_weight = wxgll(i)*wygll(j)*wzgll(k)*jacobianl

            ! compute kinetic energy  1/2 rho ||v||^2
            ! we will divide the total sum by 2 only once, at the end of this routine, to reduce compute time
            kinetic_energy = kinetic_energy + integration_weight * rhol*(vx**2 + vy**2 + vz**2)

            ! compute potential energy 1/2 sigma_ij epsilon_ij
            ! we will divide the total sum by 2 only once, at the end of this routine, to reduce compute time
            potential_energy = potential_energy + integration_weight * pressure**2 / kappal

          enddo
        enddo
      enddo

    else

      call exit_MPI(myrank,'calculation of total energy implemented for acoustic and (visco)elastic elements only for now')

    endif

  enddo

! divide the total sum by 2 here because we have purposely not done it above in order to reduce compute time
  kinetic_energy = 0.5d0 * kinetic_energy
  potential_energy = 0.5d0 * potential_energy

! compute the total using a reduction between all the processors
  call sum_all_dp(kinetic_energy,kinetic_energy_glob)
  call sum_all_dp(potential_energy,potential_energy_glob)
  total_energy_glob = kinetic_energy_glob + potential_energy_glob

! write the total to disk from the master
  if (myrank == 0) write(IOUT_ENERGY,*) it,sngl(kinetic_energy_glob),sngl(potential_energy_glob),sngl(total_energy_glob)

  end subroutine compute_energy

