!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


! for acoustic solver

  subroutine compute_forces_acoustic_Dev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic, &
                        phase_ispec_inner_acoustic,backward_simulation)

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot
!
  use specfem_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
                        m1,m2,NGLLCUBE,PML_CONDITIONS

  use pml_par, only: is_CPML, spec_to_CPML, NSPEC_CPML, &
                     potential_dot_dot_acoustic_CPML,rmemory_dpotential_dxl,rmemory_dpotential_dyl, &
                     rmemory_dpotential_dzl,rmemory_potential_acoustic, &
                     potential_acoustic_old,potential_acoustic_new

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: &
        rhostore,jacobian

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprimewgll_xx,hprimewgll_xxT

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  integer,intent(in) :: iphase
  integer,intent(in) :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(num_phase_ispec_acoustic,2),intent(in) :: phase_ispec_inner_acoustic

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: rho_invl

  integer :: ispec,iglob,ispec_p,num_elements

  ! CPML
  integer :: ispec_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem_old,chi_elem_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_old,tempx2_old,tempx3_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1_new,tempx2_new,tempx3_new
  ! derivatives of potential with respect to x, y and z
  ! in computation potential_acoustic at "n" time step is used
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl
  ! in computation of PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old
  ! we replace potential_acoustic with potential_acoustic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old
  ! we replace potential_acoustic at "n" time step with
  ! we replace potential_acoustic with potential_acoustic_old with potential_acoustic_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new


#ifdef FORCE_VECTORIZATION
! this will (purposely) give out-of-bound array accesses if run through range checking,
! thus use only for production runs with no bound checking
  integer :: ijk
#else
  integer :: i,j,k
#endif


  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! gets values for element
    DO_LOOP_IJK
      chi_elem(INDEX_IJK) = potential_acoustic(ibool(INDEX_IJK,ispec))
    ENDDO_LOOP_IJK

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for tempx1,..
    call mxm_single(hprime_xx,m1,chi_elem,tempx1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm_3dmat_single(chi_elem,m1,hprime_xxT,m1,tempx2,NGLLX)
    ! computes 3. matrix multiplication for tempx1,..
    call mxm_single(chi_elem,m2,hprime_xxT,tempx3,m1)

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
      ! gets values for element
        DO_LOOP_IJK
          chi_elem_old(INDEX_IJK) = potential_acoustic_old(ibool(INDEX_IJK,ispec))
          chi_elem_new(INDEX_IJK) = potential_acoustic_new(ibool(INDEX_IJK,ispec))
        ENDDO_LOOP_IJK

        ! subroutines adapted from Deville, Fischer and Mund, High-order methods
        ! for incompressible fluid flow, Cambridge University Press (2002),
        ! pages 386 and 389 and Figure 8.3.1

        ! computes 1. matrix multiplication for tempx1,..
        call mxm_single(hprime_xx,m1,chi_elem_old,tempx1_old,m2)
        ! computes 2. matrix multiplication for tempx2,..
        call mxm_3dmat_single(chi_elem_old,m1,hprime_xxT,m1,tempx2_old,NGLLX)
        ! computes 3. matrix multiplication for tempx1,..
        call mxm_single(chi_elem_old,m2,hprime_xxT,tempx3_old,m1)

        ! computes 1. matrix multiplication for tempx1,..
        call mxm_single(hprime_xx,m1,chi_elem_new,tempx1_new,m2)
        ! computes 2. matrix multiplication for tempx2,..
        call mxm_3dmat_single(chi_elem_new,m1,hprime_xxT,m1,tempx2_new,NGLLX)
        ! computes 3. matrix multiplication for tempx1,..
        call mxm_single(chi_elem_new,m2,hprime_xxT,tempx3_new,m1)
      endif ! is_CPML
    endif ! PML_CONDITIONS

    DO_LOOP_IJK
      ! get derivatives of potential with respect to x, y and z
      xixl = xix(INDEX_IJK,ispec)
      xiyl = xiy(INDEX_IJK,ispec)
      xizl = xiz(INDEX_IJK,ispec)
      etaxl = etax(INDEX_IJK,ispec)
      etayl = etay(INDEX_IJK,ispec)
      etazl = etaz(INDEX_IJK,ispec)
      gammaxl = gammax(INDEX_IJK,ispec)
      gammayl = gammay(INDEX_IJK,ispec)
      gammazl = gammaz(INDEX_IJK,ispec)
      jacobianl = jacobian(INDEX_IJK,ispec)

      ! derivatives of potential
      dpotentialdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
      dpotentialdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
      dpotentialdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

      ! stores derivatives of ux, uy and uz with respect to x, y and z
      if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
        ! do not merge this second line with the first using an " .and. " statement
        ! because array is_CPML() is unallocated when PML_CONDITIONS is false
        if (is_CPML(ispec)) then
          PML_dpotential_dxl(INDEX_IJK) = dpotentialdxl
          PML_dpotential_dyl(INDEX_IJK) = dpotentialdyl
          PML_dpotential_dzl(INDEX_IJK) = dpotentialdzl

          PML_dpotential_dxl_old(INDEX_IJK) = &
            xixl*tempx1_old(INDEX_IJK) + etaxl*tempx2_old(INDEX_IJK) + gammaxl*tempx3_old(INDEX_IJK)
          PML_dpotential_dyl_old(INDEX_IJK) = &
            xiyl*tempx1_old(INDEX_IJK) + etayl*tempx2_old(INDEX_IJK) + gammayl*tempx3_old(INDEX_IJK)
          PML_dpotential_dzl_old(INDEX_IJK) = &
            xizl*tempx1_old(INDEX_IJK) + etazl*tempx2_old(INDEX_IJK) + gammazl*tempx3_old(INDEX_IJK)

          PML_dpotential_dxl_new(INDEX_IJK) = &
            xixl*tempx1_new(INDEX_IJK) + etaxl*tempx2_new(INDEX_IJK) + gammaxl*tempx3_new(INDEX_IJK)
          PML_dpotential_dyl_new(INDEX_IJK) = &
            xiyl*tempx1_new(INDEX_IJK) + etayl*tempx2_new(INDEX_IJK) + gammayl*tempx3_new(INDEX_IJK)
          PML_dpotential_dzl_new(INDEX_IJK) = &
            xizl*tempx1_new(INDEX_IJK) + etazl*tempx2_new(INDEX_IJK) + gammazl*tempx3_new(INDEX_IJK)
        endif
      endif

      ! density (reciproc)
      rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)

      ! for acoustic medium
      ! also add GLL integration weights
      tempx1(INDEX_IJK) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
      tempx2(INDEX_IJK) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
      tempx3(INDEX_IJK) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
    ENDDO_LOOP_IJK

    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        ispec_CPML = spec_to_CPML(ispec)

        ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
        call pml_compute_memory_variables_acoustic(ispec,ispec_CPML, &
                                                   tempx1,tempx2,tempx3, &
                                                   rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                                                   PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                                                   PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                                                   PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new)

        ! calculates contribution from each C-PML element to update acceleration
        call pml_compute_accel_contribution_acoustic(ispec,ispec_CPML,potential_acoustic, &
                                                     potential_dot_acoustic,rmemory_potential_acoustic)
      endif
    endif ! PML_CONDITIONS

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    call mxm_single(hprimewgll_xxT,m1,tempx1,newtempx1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm_3dmat_single(tempx2,m1,hprimewgll_xx,m1,newtempx2,NGLLX)
    ! computes 3. matrix multiplication for newtempx3,..
    call mxm_single(tempx3,m2,hprimewgll_xx,newtempx3,m1)

    ! second double-loop over GLL to compute all the terms
#ifdef FORCE_VECTORIZATION
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
!DIR$ IVDEP
#endif
    DO_LOOP_IJK
      ! sum contributions from each element to the global values
      iglob = ibool(INDEX_IJK,ispec)
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                          - ( wgllwgll_yz_3D(INDEX_IJK)*newtempx1(INDEX_IJK) &
                                            + wgllwgll_xz_3D(INDEX_IJK)*newtempx2(INDEX_IJK) &
                                            + wgllwgll_xy_3D(INDEX_IJK)*newtempx3(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! updates potential_dot_dot_acoustic with contribution from each C-PML element
    if (PML_CONDITIONS .and. (.not. backward_simulation) .and. NSPEC_CPML > 0) then
      ! do not merge this second line with the first using an " .and. " statement
      ! because array is_CPML() is unallocated when PML_CONDITIONS is false
      if (is_CPML(ispec)) then
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_CPML(INDEX_IJK)
        ENDDO_LOOP_IJK
      endif
    endif ! PML_CONDITIONS

  enddo ! end of loop over all spectral elements

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

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
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm_single(A,n1,B,C,n3)

! 2-dimensional arrays e.g. (25,5)/(5,25)

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,NGLLX),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(NGLLX,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! matrix-matrix multiplication wrapper
  if (NGLLX == 5) then
    call mxm5_single(A,n1,B,C,n3)
  else if (NGLLX == 6) then
    call mxm6_single(A,n1,B,C,n3)
  else if (NGLLX == 7) then
    call mxm7_single(A,n1,B,C,n3)
  endif

  end subroutine mxm_single

  !-------------

  subroutine mxm5_single(A,n1,B,C,n3)

! 2-dimensional arrays (25,5)/(5,25)

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

! 2-dimensional arrays (36,6)/(6,36)

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

! 2-dimensional arrays (49,7)/(7,49)

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

  subroutine mxm_3dmat_single(A,n1,B,n2,C,n3)

! 3-dimensional arrays e.g. (5,5,5) for A and C

  use constants, only: CUSTOM_REAL,NGLLX

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,NGLLX,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(NGLLX,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! matrix-matrix multiplication wrapper
  if (NGLLX == 5) then
    call mxm5_3dmat_single(A,n1,B,n2,C,n3)
  else if (NGLLX == 6) then
    call mxm6_3dmat_single(A,n1,B,n2,C,n3)
  else if (NGLLX == 7) then
    call mxm7_3dmat_single(A,n1,B,n2,C,n3)
  endif

  end subroutine mxm_3dmat_single

  !-------------

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! 3-dimensional arrays (5,5,5) for A and C

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

! 3-dimensional arrays (6,6,6) for A and C

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

! 3-dimensional arrays (6,6,6) for A and C

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

  end subroutine compute_forces_acoustic_Dev

