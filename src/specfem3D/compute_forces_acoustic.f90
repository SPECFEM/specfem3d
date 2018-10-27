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

  subroutine compute_forces_acoustic_generic_slow(iphase, &
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

! all the loops below contain critical operations in which the solver spends 90% or so of its time;
! on modern machines, to get good performance it is thus crucial to make sure that most of these loops are fully vectorized
! by the compiler; we thus have several dedicated versions and also manually-unrolled loops, as well as small temporary arrays,
! whose goal is to ensure that compilers can vectorized them and make optimized use of the processor cache as well.

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

              if (ispec_irreg /= 0 ) then ! irregular element

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

          if (ispec_irreg /= 0 ) then ! irregular element

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

  end subroutine compute_forces_acoustic_generic_slow

!
!=====================================================================
!

! for acoustic solver

  subroutine compute_forces_acoustic_fast_Deville(iphase, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        backward_simulation)

! computes forces for acoustic elements

  use specfem_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
                         chi_elem,chi_elem_old,chi_elem_new, &
                         temp1,temp2,temp3, &
                         temp1_old,temp2_old,temp3_old, &
                         temp1_new,temp2_new,temp3_new, &
                         newtemp1,newtemp2,newtemp3, &
                         PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                         PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                         PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new, &
                         NGLOB_AB, &
                         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                         hprime_xx,hprime_xxT, &
                         hprimewgll_xx,hprimewgll_xxT, &
                         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                         rhostore,jacobian,ibool, &
                         irregular_element_number,xix_regular,jacobian_regular,m1,m2

  use specfem_par_acoustic, only: nspec_inner_acoustic,nspec_outer_acoustic, &
                                   phase_ispec_inner_acoustic

  use pml_par, only: is_CPML, spec_to_CPML, potential_dot_dot_acoustic_CPML,rmemory_dpotential_dxl,rmemory_dpotential_dyl, &
                     rmemory_dpotential_dzl,rmemory_potential_acoustic, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  integer, intent(in) :: iphase
  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

  ! local variables
  integer :: ispec,ispec_irreg,iglob,i,j,k,ispec_p,num_elements

  ! CPML
  integer :: ispec_CPML

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

! all the loops below contain critical operations in which the solver spends 90% or so of its time; on modern machines,
! to get good performance it is thus crucial to make sure that most of these loops are fully vectorized by the compiler.
! In particular, PLEASE NEVER PUT ANY "IF" STATEMENT INSIDE THESE LOOPS;
! if you need to perform an "if" test, put it outside the loop
! and duplicate the content of the loop, using one for each case of the "if" test.

    ! gets value of the field inside the element and make it local
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
        enddo
      enddo
    enddo

    ! derivative along x, y, z
    ! first double loop over GLL points to compute and store gradients

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for temp1
    ! computes 2. matrix multiplication for temp2
    ! computes 3. matrix multiplication for temp3
  if (NGLLX == 5) then
    call mxm5_single(hprime_xx,m1,chi_elem,temp1,m2)
    call mxm5_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm5_single(chi_elem,m2,hprime_xxT,temp3,m1)
  else if (NGLLX == 6) then
    call mxm6_single(hprime_xx,m1,chi_elem,temp1,m2)
    call mxm6_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm6_single(chi_elem,m2,hprime_xxT,temp3,m1)
  else if (NGLLX == 7) then
    call mxm7_single(hprime_xx,m1,chi_elem,temp1,m2)
    call mxm7_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm7_single(chi_elem,m2,hprime_xxT,temp3,m1)
  endif

!
!-----------------
!

  if (is_CPML(ispec) .and. .not. backward_simulation) then

    ispec_CPML = spec_to_CPML(ispec)
    ispec_irreg = irregular_element_number(ispec)

    ! gets value of the field inside the element and make it local
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          chi_elem_old(i,j,k) = PML_potential_acoustic_old(i,j,k,ispec_CPML)
          chi_elem_new(i,j,k) = PML_potential_acoustic_new(i,j,k,ispec_CPML)
        enddo
      enddo
    enddo

        ! subroutines adapted from Deville, Fischer and Mund, High-order methods
        ! for incompressible fluid flow, Cambridge University Press (2002),
        ! pages 386 and 389 and Figure 8.3.1

        ! computes 1. matrix multiplication for temp1
        ! computes 2. matrix multiplication for temp2
        ! computes 3. matrix multiplication for temp3
    if (NGLLX == 5) then
        call mxm5_single_two_arrays_at_a_time(hprime_xx,m1,chi_elem_old,temp1_old,m2,chi_elem_new,temp1_new)
        call mxm5_3dmat_single_two_arrays_at_a_time(chi_elem_old,m1,hprime_xxT,m1,temp2_old,NGLLX,chi_elem_new,temp1_new)
        call mxm5_single_two_arrays_at_a_time(chi_elem_old,m2,hprime_xxT,temp3_old,m1,chi_elem_new,temp1_new)
    else if (NGLLX == 6) then
        call mxm6_single_two_arrays_at_a_time(hprime_xx,m1,chi_elem_old,temp1_old,m2,chi_elem_new,temp1_new)
        call mxm6_3dmat_single_two_arrays_at_a_time(chi_elem_old,m1,hprime_xxT,m1,temp2_old,NGLLX,chi_elem_new,temp1_new)
        call mxm6_single_two_arrays_at_a_time(chi_elem_old,m2,hprime_xxT,temp3_old,m1,chi_elem_new,temp1_new)
    else if (NGLLX == 7) then
        call mxm7_single_two_arrays_at_a_time(hprime_xx,m1,chi_elem_old,temp1_old,m2,chi_elem_new,temp1_new)
        call mxm7_3dmat_single_two_arrays_at_a_time(chi_elem_old,m1,hprime_xxT,m1,temp2_old,NGLLX,chi_elem_new,temp1_new)
        call mxm7_single_two_arrays_at_a_time(chi_elem_old,m2,hprime_xxT,temp3_old,m1,chi_elem_new,temp1_new)
    endif

    if (ispec_irreg /= 0) then ! irregular element

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

              ! reciprocal of density
              rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

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
                dpotentialdxl = xixl*temp1(i,j,k) + etaxl*temp2(i,j,k) + gammaxl*temp3(i,j,k)
                dpotentialdyl = xiyl*temp1(i,j,k) + etayl*temp2(i,j,k) + gammayl*temp3(i,j,k)
                dpotentialdzl = xizl*temp1(i,j,k) + etazl*temp2(i,j,k) + gammazl*temp3(i,j,k)

                PML_dpotential_dxl_old(i,j,k) = xixl*temp1_old(i,j,k) + etaxl*temp2_old(i,j,k) + gammaxl*temp3_old(i,j,k)
                PML_dpotential_dyl_old(i,j,k) = xiyl*temp1_old(i,j,k) + etayl*temp2_old(i,j,k) + gammayl*temp3_old(i,j,k)
                PML_dpotential_dzl_old(i,j,k) = xizl*temp1_old(i,j,k) + etazl*temp2_old(i,j,k) + gammazl*temp3_old(i,j,k)

                PML_dpotential_dxl_new(i,j,k) = xixl*temp1_new(i,j,k) + etaxl*temp2_new(i,j,k) + gammaxl*temp3_new(i,j,k)
                PML_dpotential_dyl_new(i,j,k) = xiyl*temp1_new(i,j,k) + etayl*temp2_new(i,j,k) + gammayl*temp3_new(i,j,k)
                PML_dpotential_dzl_new(i,j,k) = xizl*temp1_new(i,j,k) + etazl*temp2_new(i,j,k) + gammazl*temp3_new(i,j,k)

                ! for acoustic medium
                temp1(i,j,k) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
                temp2(i,j,k) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
                temp3(i,j,k) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

          ! stores derivatives of ux, uy and uz with respect to x, y and z
          PML_dpotential_dxl(i,j,k) = dpotentialdxl
          PML_dpotential_dyl(i,j,k) = dpotentialdyl
          PML_dpotential_dzl(i,j,k) = dpotentialdzl

        enddo
      enddo
    enddo

    else ! regular element

    jacobianl = jacobian_regular

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

              ! reciprocal of density
              rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

                dpotentialdxl = xix_regular*temp1(i,j,k)
                dpotentialdyl = xix_regular*temp2(i,j,k)
                dpotentialdzl = xix_regular*temp3(i,j,k)

                PML_dpotential_dxl_old(i,j,k) = xix_regular*temp1_old(i,j,k)
                PML_dpotential_dyl_old(i,j,k) = xix_regular*temp2_old(i,j,k)
                PML_dpotential_dzl_old(i,j,k) = xix_regular*temp3_old(i,j,k)

                PML_dpotential_dxl_new(i,j,k) = xix_regular*temp1_new(i,j,k)
                PML_dpotential_dyl_new(i,j,k) = xix_regular*temp2_new(i,j,k)
                PML_dpotential_dzl_new(i,j,k) = xix_regular*temp3_new(i,j,k)

                ! for acoustic medium
                temp1(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdxl
                temp2(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdyl
                temp3(i,j,k) = rho_invl * jacobianl * xix_regular * dpotentialdzl

          ! stores derivatives of ux, uy and uz with respect to x, y and z
          PML_dpotential_dxl(i,j,k) = dpotentialdxl
          PML_dpotential_dyl(i,j,k) = dpotentialdyl
          PML_dpotential_dzl(i,j,k) = dpotentialdzl

        enddo
      enddo
    enddo

    endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  else ! no PML in this element

    ispec_irreg = irregular_element_number(ispec)

    if (ispec_irreg /= 0 ) then ! irregular element

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! reciprocal of density
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

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
            dpotentialdxl = xixl*temp1(i,j,k) + etaxl*temp2(i,j,k) + gammaxl*temp3(i,j,k)
            dpotentialdyl = xiyl*temp1(i,j,k) + etayl*temp2(i,j,k) + gammayl*temp3(i,j,k)
            dpotentialdzl = xizl*temp1(i,j,k) + etazl*temp2(i,j,k) + gammazl*temp3(i,j,k)

            temp1(i,j,k) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
            temp2(i,j,k) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
            temp3(i,j,k) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

        enddo
      enddo
    enddo

  else ! regular element

    jacobianl = jacobian_regular

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          ! reciprocal of density
          rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec)

          ! for acoustic medium
          temp1(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp1(i,j,k)
          temp2(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp2(i,j,k)
          temp3(i,j,k) = rho_invl * jacobianl * xix_regular * xix_regular * temp3(i,j,k)

        enddo
      enddo
    enddo

  endif

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

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtemp1
    ! computes 2. matrix multiplication for newtemp2
    ! computes 3. matrix multiplication for newtemp3
  if (NGLLX == 5) then
    call mxm5_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    call mxm5_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    call mxm5_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
  else if (NGLLX == 6) then
    call mxm6_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    call mxm6_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    call mxm6_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
  else if (NGLLX == 7) then
    call mxm7_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    call mxm7_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    call mxm7_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
  endif

    ! sum contributions from each element to the global values
    ! this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! also add GLL integration weights
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
            (wgllwgll_yz(j,k)*newtemp1(i,j,k) + wgllwgll_xz(i,k)*newtemp2(i,j,k) + wgllwgll_xy(i,j)*newtemp3(i,j,k))
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

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices (5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications is in general slower
!
! please leave the routines here to help compilers inline the code

  subroutine mxm5_single(A,n1,B,C,n3)

! two-dimensional arrays (25,5)/(5,25)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_single

  !-------------

  subroutine mxm6_single(A,n1,B,C,n3)

! two-dimensional arrays (36,6)/(6,36)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j)
    enddo
  enddo

  end subroutine mxm6_single

  !-------------

  subroutine mxm7_single(A,n1,B,C,n3)

! two-dimensional arrays (49,7)/(7,49)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j) &
              + A(i,7) * B(7,j)
    enddo
  enddo

  end subroutine mxm7_single


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (5,5,5) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3dmat_single

  !-------------

  subroutine mxm6_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (6,6,6) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,6,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j)
      enddo
    enddo
  enddo

  end subroutine mxm6_3dmat_single

  !-------------

  subroutine mxm7_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (7,7,7) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,7,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j) &
                  + A(i,7,k) * B(7,j)
      enddo
    enddo
  enddo

  end subroutine mxm7_3dmat_single

  !-------------

  subroutine mxm5_single_two_arrays_at_a_time(A,n1,B,C,n3,B2,C2)

! two-dimensional arrays (25,5)/(5,25)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B,B2
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)
    enddo
  enddo

  end subroutine mxm5_single_two_arrays_at_a_time

  !-------------

  subroutine mxm6_single_two_arrays_at_a_time(A,n1,B,C,n3,B2,C2)

! two-dimensional arrays (36,6)/(6,36)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B,B2
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j) &
               + A(i,6) * B2(6,j)
    enddo
  enddo

  end subroutine mxm6_single_two_arrays_at_a_time

  !-------------

  subroutine mxm7_single_two_arrays_at_a_time(A,n1,B,C,n3,B2,C2)

! two-dimensional arrays (49,7)/(7,49)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B,B2
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j) &
              + A(i,7) * B(7,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j) &
               + A(i,6) * B2(6,j) &
               + A(i,7) * B2(7,j)
    enddo
  enddo

  end subroutine mxm7_single_two_arrays_at_a_time


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3dmat_single_two_arrays_at_a_time(A,n1,B,n2,C,n3,A2,C2)

! three-dimensional arrays (5,5,5) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A,A2
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3dmat_single_two_arrays_at_a_time

  !-------------

  subroutine mxm6_3dmat_single_two_arrays_at_a_time(A,n1,B,n2,C,n3,A2,C2)

! three-dimensional arrays (6,6,6) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,6,n3),intent(in) :: A,A2
  real(kind=CUSTOM_REAL),dimension(6,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j) &
                   + A2(i,6,k) * B(6,j)
      enddo
    enddo
  enddo

  end subroutine mxm6_3dmat_single_two_arrays_at_a_time

  !-------------

  subroutine mxm7_3dmat_single_two_arrays_at_a_time(A,n1,B,n2,C,n3,A2,C2)

! three-dimensional arrays (7,7,7) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,7,n3),intent(in) :: A,A2
  real(kind=CUSTOM_REAL),dimension(7,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C,C2

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j) &
                  + A(i,7,k) * B(7,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j) &
                   + A2(i,6,k) * B(6,j) &
                   + A2(i,7,k) * B(7,j)
      enddo
    enddo
  enddo

  end subroutine mxm7_3dmat_single_two_arrays_at_a_time

  end subroutine compute_forces_acoustic_fast_Deville

