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

  subroutine compute_forces_acoustic_NGLLnot5_generic_slow(iphase, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        backward_simulation)

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot

  use specfem_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,chi_elem,temp1,temp2,temp3,temp4, &
                         PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                         PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                         PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new, &
                         NGLOB_AB, &
                         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                         hprime_xx,hprime_yy,hprime_zz, &
                         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                         rhostore,jacobian,ibool, &
                         irregular_element_number,xix_regular,jacobian_regular

  use specfem_par_acoustic, only: nspec_inner_acoustic,nspec_outer_acoustic, &
                                   phase_ispec_inner_acoustic


  use pml_par, only: is_CPML, spec_to_CPML, potential_dot_dot_acoustic_CPML,rmemory_dpotential_dxl,rmemory_dpotential_dyl, &
                     rmemory_dpotential_dzl,rmemory_potential_acoustic, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  integer,intent(in) :: iphase
  logical,intent(in) :: backward_simulation

  ! local variables
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3

  integer :: ispec,ispec_irreg,iglob,i,j,k,l,ispec_p,num_elements

  ! CPML
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL) :: temp1l_old,temp2l_old,temp3l_old
  real(kind=CUSTOM_REAL) :: temp1l_new,temp2l_new,temp3l_new

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: rho_invl

  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

!
! all the loops below contain critical operations in which the solver spends 90% or so of its time;
! on modern machines, to get good performance it is thus crucial to make sure that most of these loops are fully vectorized
! by the compiler; we thus have several dedicated versions and also manually-unrolled loops, as well as small temporary arrays,
! whose goal is to ensure that compilers can vectorized them and make optimized use of the processor cache as well.
!

!
!-----------------
!

    ! gets value of the field inside the element and make it local
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
        enddo
      enddo
    enddo

!
!-----------------
!

  if (is_CPML(ispec) .and. .not. backward_simulation) then

    ispec_CPML = spec_to_CPML(ispec)
    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg == 0) jacobianl = jacobian_regular
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! derivative along x, y, z
          ! first double loop over GLL points to compute and store gradients
          temp1l = 0._CUSTOM_REAL
          temp2l = 0._CUSTOM_REAL
          temp3l = 0._CUSTOM_REAL

          ! we can merge the loops because NGLLX == NGLLY == NGLLZ
          do l = 1,NGLLX
            temp1l = temp1l + chi_elem(l,j,k)*hprime_xx(i,l)
            temp2l = temp2l + chi_elem(i,l,k)*hprime_yy(j,l)
            temp3l = temp3l + chi_elem(i,j,l)*hprime_zz(k,l)
          enddo

              temp1l_old = 0._CUSTOM_REAL
              temp2l_old = 0._CUSTOM_REAL
              temp3l_old = 0._CUSTOM_REAL

              temp1l_new = 0._CUSTOM_REAL
              temp2l_new = 0._CUSTOM_REAL
              temp3l_new = 0._CUSTOM_REAL

              ! we can merge these loops because NGLLX = NGLLY = NGLLZ
              do l=1,NGLLX
                hp1 = hprime_xx(i,l)
                iglob = ibool(l,j,k,ispec)
                temp1l_old = temp1l_old + PML_potential_acoustic_old(l,j,k,ispec_CPML)*hp1
                temp1l_new = temp1l_new + PML_potential_acoustic_new(l,j,k,ispec_CPML)*hp1

                hp2 = hprime_yy(j,l)
                iglob = ibool(i,l,k,ispec)
                temp2l_old = temp2l_old + PML_potential_acoustic_old(i,l,k,ispec_CPML)*hp2
                temp2l_new = temp2l_new + PML_potential_acoustic_new(i,l,k,ispec_CPML)*hp2

                hp3 = hprime_zz(k,l)
                iglob = ibool(i,j,l,ispec)
                temp3l_old = temp3l_old + PML_potential_acoustic_old(i,j,l,ispec_CPML)*hp3
                temp3l_new = temp3l_new + PML_potential_acoustic_new(i,j,l,ispec_CPML)*hp3
              enddo

              rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

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
                jacobianl = jacobian(i,j,k,ispec_irreg)

                ! derivatives of potential
                dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l
                dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l
                dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l

                PML_dpotential_dxl_old(i,j,k) = xixl*temp1l_old + etaxl*temp2l_old + gammaxl*temp3l_old
                PML_dpotential_dyl_old(i,j,k) = xiyl*temp1l_old + etayl*temp2l_old + gammayl*temp3l_old
                PML_dpotential_dzl_old(i,j,k) = xizl*temp1l_old + etazl*temp2l_old + gammazl*temp3l_old

                PML_dpotential_dxl_new(i,j,k) = xixl*temp1l_new + etaxl*temp2l_new + gammaxl*temp3l_new
                PML_dpotential_dyl_new(i,j,k) = xiyl*temp1l_new + etayl*temp2l_new + gammayl*temp3l_new
                PML_dpotential_dzl_new(i,j,k) = xizl*temp1l_new + etazl*temp2l_new + gammazl*temp3l_new

                ! for acoustic medium
                temp1(i,j,k) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
                temp2(i,j,k) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
                temp3(i,j,k) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

             else ! regular element

                dpotentialdxl = xix_regular*temp1l
                dpotentialdyl = xix_regular*temp2l
                dpotentialdzl = xix_regular*temp3l

                PML_dpotential_dxl_old(i,j,k) = xix_regular*temp1l_old
                PML_dpotential_dyl_old(i,j,k) = xix_regular*temp2l_old
                PML_dpotential_dzl_old(i,j,k) = xix_regular*temp3l_old

                PML_dpotential_dxl_new(i,j,k) = xix_regular*temp1l_new
                PML_dpotential_dyl_new(i,j,k) = xix_regular*temp2l_new
                PML_dpotential_dzl_new(i,j,k) = xix_regular*temp3l_new

                ! for acoustic medium
                temp1(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdxl
                temp2(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdyl
                temp3(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdzl

              endif

              ! stores derivatives of ux, uy and uz with respect to x, y and z
              PML_dpotential_dxl(i,j,k) = dpotentialdxl
              PML_dpotential_dyl(i,j,k) = dpotentialdyl
              PML_dpotential_dzl(i,j,k) = dpotentialdzl
        enddo
      enddo
    enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  else ! no PML in this element

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg == 0) jacobianl = jacobian_regular

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! derivative along x, y, z
          ! first double loop over GLL points to compute and store gradients
          temp1l = 0._CUSTOM_REAL
          temp2l = 0._CUSTOM_REAL
          temp3l = 0._CUSTOM_REAL

          ! we can merge the loops because NGLLX == NGLLY == NGLLZ
          do l = 1,NGLLX
            temp1l = temp1l + chi_elem(l,j,k)*hprime_xx(i,l)
            temp2l = temp2l + chi_elem(i,l,k)*hprime_yy(j,l)
            temp3l = temp3l + chi_elem(i,j,l)*hprime_zz(k,l)
          enddo

          ! density (reciproc)
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

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
            jacobianl = jacobian(i,j,k,ispec_irreg)

            ! derivatives of potential
            dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l
            dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l
            dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l

            temp1(i,j,k) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
            temp2(i,j,k) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
            temp3(i,j,k) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

          else ! regular element
            ! for acoustic medium
            temp1(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp1l
            temp2(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp2l
            temp3(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp3l

          endif

        enddo
      enddo
    enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  endif ! of if PML in this element

!
!-----------------
!

    if (is_CPML(ispec) .and. .not. backward_simulation) then
        ispec_CPML = spec_to_CPML(ispec)

        ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
        call pml_compute_memory_variables_acoustic(ispec,ispec_CPML,temp1,temp2,temp3, &
                                                   rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                                                   PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                                                   PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                                                   PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new)

        ! calculates contribution from each C-PML element to update acceleration
        call pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic, &
                                                     potential_dot_acoustic,rmemory_potential_acoustic)
    endif

!
!-----------------
!

    ! second double-loop over GLL to compute all the terms along the x,y,z directions and assemble the contributions
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          temp1l = 0._CUSTOM_REAL
          temp2l = 0._CUSTOM_REAL
          temp3l = 0._CUSTOM_REAL

          ! we can merge these loops because NGLLX = NGLLY = NGLLZ
          do l=1,NGLLX
            temp1l = temp1l + temp1(l,j,k) * hprimewgll_xx(l,i)
            temp2l = temp2l + temp2(i,l,k) * hprimewgll_yy(l,j)
            temp3l = temp3l + temp3(i,j,l) * hprimewgll_zz(l,k)
          enddo

          ! also add GLL integration weights
          temp4(i,j,k) = - ( wgllwgll_yz(j,k)*temp1l + wgllwgll_xz(i,k)*temp2l + wgllwgll_xy(i,j)*temp3l )

        enddo
      enddo
    enddo

!
!-----------------
!

    ! sum contributions from each element to the global values
    ! this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + temp4(i,j,k)
        enddo
      enddo
    enddo

!
!-----------------
!

    ! sum contributions from each element to the global values
    ! this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    if (is_CPML(ispec) .and. .not. backward_simulation) then
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_CPML(i,j,k)
           enddo
         enddo
       enddo
   endif

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_acoustic_NGLLnot5_generic_slow

