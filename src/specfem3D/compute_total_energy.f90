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


  subroutine compute_total_energy()

! computes kinetic, potential and total energy
! in all the slices using an MPI reduction
! and output that to an energy file

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use pml_par

  implicit none

! local variables
  integer :: i,j,k,l,ispec,iglob

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) :: epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz,epsilon_yx,epsilon_zx,epsilon_zy
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) :: vx,vy,vz,pressure

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul,rhol,cpl
  real(kind=CUSTOM_REAL) :: kappal

  real(kind=CUSTOM_REAL) :: integration_weight
  double precision :: kinetic_energy,potential_energy
  double precision :: kinetic_energy_glob,potential_energy_glob,total_energy_glob

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  kinetic_energy = 0.d0
  potential_energy = 0.d0

  if(ANISOTROPY .or. ATTENUATION) &
    call exit_MPI(myrank,'calculation of total energy currently implemented only for media with no anisotropy and no attenuation')

! loop over spectral elements
  do ispec = 1,NSPEC_AB

! if element is a CPML then do not compute energy in it, since it is non physical;
! thus, we compute energy in the main domain only, without absorbing elements
    if(PML_CONDITIONS) then
      ! do not merge this second line with the first using an ".and." statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if(is_CPML(ispec)) cycle
    endif

    !---
    !--- elastic spectral element
    !---
    if(ispec_is_elastic(ispec)) then

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

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1
            tempy1(i,j,k) = tempy1(i,j,k) + dummyy_loc(l,j,k)*hp1
            tempz1(i,j,k) = tempz1(i,j,k) + dummyz_loc(l,j,k)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2
            tempy2(i,j,k) = tempy2(i,j,k) + dummyy_loc(i,l,k)*hp2
            tempz2(i,j,k) = tempz2(i,j,k) + dummyz_loc(i,l,k)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
            tempy3(i,j,k) = tempy3(i,j,k) + dummyy_loc(i,j,l)*hp3
            tempz3(i,j,k) = tempz3(i,j,k) + dummyz_loc(i,j,l)*hp3
          enddo

              ! get derivatives of ux, uy and uz with respect to x, y and z
              xixl = xix(i,j,k,ispec)
              xiyl = xiy(i,j,k,ispec)
              xizl = xiz(i,j,k,ispec)
              etaxl = etax(i,j,k,ispec)
              etayl = etay(i,j,k,ispec)
              etazl = etaz(i,j,k,ispec)
              gammaxl = gammax(i,j,k,ispec)
              gammayl = gammay(i,j,k,ispec)
              gammazl = gammaz(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
              duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
              duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

              duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
              duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
              duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

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
              epsilon_xy = 0.5 * duxdyl_plus_duydxl
              epsilon_xz = 0.5 * duzdxl_plus_duxdzl
              epsilon_yz = 0.5 * duzdyl_plus_duydzl

              ! define symmetric components of epsilon
              epsilon_yx = epsilon_xy
              epsilon_zx = epsilon_xz
              epsilon_zy = epsilon_yz

              kappal = kappastore(i,j,k,ispec)
              mul = mustore(i,j,k,ispec)
              rhol = rhostore(i,j,k,ispec)

              ! isotropic case
              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2.*mul

              ! compute stress sigma
              sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
              sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
              sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

              sigma_xy = mul*duxdyl_plus_duydxl
              sigma_xz = mul*duzdxl_plus_duxdzl
              sigma_yz = mul*duzdyl_plus_duydzl

              ! define symmetric components of sigma
              sigma_yx = sigma_xy
              sigma_zx = sigma_xz
              sigma_zy = sigma_yz

              integration_weight = wxgll(i)*wygll(j)*wzgll(k)*jacobianl

              ! compute kinetic energy  1/2 rho ||v||^2
              kinetic_energy = kinetic_energy + integration_weight * rhol*(veloc(1,iglob)**2 + &
                                   veloc(2,iglob)**2 + veloc(3,iglob)**2) / 2.

              ! compute potential energy 1/2 sigma_ij epsilon_ij
              potential_energy = potential_energy + integration_weight * &
                (sigma_xx*epsilon_xx + sigma_xy*epsilon_xy + sigma_xz*epsilon_xz + &
                 sigma_yx*epsilon_yx + sigma_yy*epsilon_yy + sigma_yz*epsilon_yz + &
                 sigma_zx*epsilon_zx + sigma_zy*epsilon_zy + sigma_zz*epsilon_zz) / 2.

          enddo
        enddo
     enddo

    !---
    !--- acoustic spectral element
    !---
    else if(ispec_is_acoustic(ispec)) then

      ! for the definition of potential energy in an acoustic fluid, see for instance
      ! equation (23) of M. Maess et al., Journal of Sound and Vibration 296 (2006) 264-276

      ! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
      ! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
      ! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
      ! This permits acoustic-elastic coupling based on a non-iterative time scheme.
      ! Displacement is then: u = grad(Chi) / rho
      ! Velocity is then: v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
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

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1(i,j,k) = tempx1(i,j,k) + dummyx_loc(l,j,k)*hp1

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp2 = hprime_yy(j,l)
            tempx2(i,j,k) = tempx2(i,j,k) + dummyx_loc(i,l,k)*hp2

            !!! can merge these loops because NGLLX = NGLLY = NGLLZ
            hp3 = hprime_zz(k,l)
            tempx3(i,j,k) = tempx3(i,j,k) + dummyx_loc(i,j,l)*hp3
          enddo

              ! get derivatives of ux, uy and uz with respect to x, y and z
              xixl = xix(i,j,k,ispec)
              xiyl = xiy(i,j,k,ispec)
              xizl = xiz(i,j,k,ispec)
              etaxl = etax(i,j,k,ispec)
              etayl = etay(i,j,k,ispec)
              etazl = etaz(i,j,k,ispec)
              gammaxl = gammax(i,j,k,ispec)
              gammayl = gammay(i,j,k,ispec)
              gammazl = gammaz(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)

              duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
              duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
              duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

              rhol = rhostore(i,j,k,ispec)
              kappal = kappastore(i,j,k,ispec)
              cpl = sqrt(kappal / rhol)

              ! Velocity is v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
              vx = duxdxl / rhol
              vy = duxdyl / rhol
              vz = duxdzl / rhol

              ! pressure is p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi)
              pressure = - potential_dot_dot_acoustic(iglob)

              integration_weight = wxgll(i)*wygll(j)*wzgll(k)*jacobianl

              ! compute kinetic energy  1/2 rho ||v||^2
              kinetic_energy = kinetic_energy + integration_weight * rhol*(vx**2 + vy**2 + vz**2) / 2.

              ! compute potential energy 1/2 sigma_ij epsilon_ij
              potential_energy = potential_energy + integration_weight * pressure**2 / (2. * rhol * cpl**2)

          enddo
        enddo
     enddo

    else

      call exit_MPI(myrank,'calculation of total energy implemented for acoustic and (visco)elastic elements only for now')

    endif

  enddo

! compute the total using a reduction between all the processors
  call sum_all_dp(kinetic_energy,kinetic_energy_glob)
  call sum_all_dp(potential_energy,potential_energy_glob)
  total_energy_glob = kinetic_energy_glob + potential_energy_glob

! write the total to disk from the master
  if(myrank == 0) write(IOUT_ENERGY,*) it,sngl(kinetic_energy_glob),sngl(potential_energy_glob),sngl(total_energy_glob)

  end subroutine compute_total_energy
