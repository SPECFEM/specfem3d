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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


! for acoustic solver

  subroutine compute_forces_acoustic(iphase, &
                                     potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                     backward_simulation)

! computes forces for acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,m1,m2

  use specfem_par, only: NGLOB_AB, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,jacobianstore, &
                         hprime_xx,hprime_xxT, &
                         hprimewgll_xx,hprimewgll_xxT, &
                         hprime_yy,hprime_zz, &
                         hprimewgll_yy,hprimewgll_zz, &
                         rhostore,ibool, &
                         irregular_element_number,xix_regular,jacobian_regular

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D
  !or: use specfem_par, only: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  ! USE_DERIV_MAPPING_FUSION not used yet: use specfem_par, only: deriv_mapping

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  use specfem_par_acoustic, only: nspec_inner_acoustic,nspec_outer_acoustic, &
                                   phase_ispec_inner_acoustic

  use pml_par, only: is_CPML, NSPEC_CPML

  implicit none

  integer, intent(in) :: iphase

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_dot_dot_acoustic

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

  ! local variables

  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...
  !
  ! arrays for elemental computations inside a given spectral element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtemp1,newtemp2,newtemp3

  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: rho_invl

  integer :: i,j,k,l,ispec,ispec_irreg,iglob,ispec_p,num_elements
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

! note: we combine here previous generic_slow and fast_Deville versions of compute_forces_acoustic into a single routine.
!       Deville optimized matrix-matrix multiplications are used for NGLLX == 5, 6, 7, 8.
!       Generic loops are used for all other NGLL options.
!       Having a single routine to compute this stiffness term for acoustic elements will ease maintaining the code.

  ! number of elements to loop
  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

! openmp solver
!$OMP PARALLEL if (num_elements > 100) &
!$OMP DEFAULT(NONE) &
!$OMP SHARED( &
!$OMP num_elements,ibool, &
!$OMP iphase,phase_ispec_inner_acoustic, &
!$OMP irregular_element_number,jacobian_regular,xix_regular, &
!$OMP potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
!$OMP is_CPML,backward_simulation, &
!$OMP xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore,rhostore &
!$OMP ) &
!$OMP PRIVATE( &
!$OMP ispec_p,ispec,ispec_irreg,i,j,k,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#endif
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl,rho_invl, &
!$OMP dpotentialdxl,dpotentialdyl,dpotentialdzl, &
!$OMP temp1l,temp2l,temp3l, &
!$OMP chi_elem, &
!$OMP temp1,temp2,temp3, &
!$OMP newtemp1,newtemp2,newtemp3 &
!$OMP ) &
!$OMP FIRSTPRIVATE( &
!$OMP hprime_xx,hprime_xxT, hprimewgll_xxT, hprimewgll_xx, &
!$OMP hprime_yy,hprime_zz,hprimewgll_yy,hprimewgll_zz, &
!$OMP wgllwgll_yz_3D,wgllwgll_xz_3D,wgllwgll_xy_3D &
!$OMP )

  ! loop over spectral elements
!$OMP DO
  do ispec_p = 1,num_elements

    ! selects element contribution
    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! PML contributions will be computed after this main loop
    ! (this loop is intended for "pure" acoustic elements, otherwise performance might suffer)
    if (is_CPML(ispec) .and. .not. backward_simulation) cycle

! all the loops below contain critical operations in which the solver spends 90% or so of its time; on modern machines,
! to get good performance it is thus crucial to make sure that most of these loops are fully vectorized by the compiler.
! In particular, PLEASE NEVER PUT ANY "IF" STATEMENT INSIDE THESE LOOPS;
! if you need to perform an "if" test, put it outside the loop
! and duplicate the content of the loop, using one for each case of the "if" test.

    ! gets value of the field inside the element and make it local
    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          chi_elem(i,j,k) = potential_acoustic(iglob)
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
    select case (NGLLX)
    case (5)
      call mxm5_single(hprime_xx,m1,chi_elem,temp1,m2)
      call mxm5_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
      call mxm5_single(chi_elem,m2,hprime_xxT,temp3,m1)
    case (6)
      call mxm6_single(hprime_xx,m1,chi_elem,temp1,m2)
      call mxm6_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
      call mxm6_single(chi_elem,m2,hprime_xxT,temp3,m1)
    case (7)
      call mxm7_single(hprime_xx,m1,chi_elem,temp1,m2)
      call mxm7_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
      call mxm7_single(chi_elem,m2,hprime_xxT,temp3,m1)
    case (8)
      call mxm8_single(hprime_xx,m1,chi_elem,temp1,m2)
      call mxm8_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
      call mxm8_single(chi_elem,m2,hprime_xxT,temp3,m1)
    case default
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
            temp1(i,j,k) = temp1l
            temp2(i,j,k) = temp2l
            temp3(i,j,k) = temp3l
          enddo
        enddo
      enddo
    end select

    ! grad(potential)
    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg /= 0) then
      ! irregular element
      DO_LOOP_IJK
        ! reciprocal of density
        rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)

        ! get derivatives of ux, uy and uz with respect to x, y and z
        ! single arrays make it difficult for hardware pre-fetching...
        xixl = xixstore(INDEX_IJK,ispec_irreg)
        xiyl = xiystore(INDEX_IJK,ispec_irreg)
        xizl = xizstore(INDEX_IJK,ispec_irreg)
        etaxl = etaxstore(INDEX_IJK,ispec_irreg)
        etayl = etaystore(INDEX_IJK,ispec_irreg)
        etazl = etazstore(INDEX_IJK,ispec_irreg)
        gammaxl = gammaxstore(INDEX_IJK,ispec_irreg)
        gammayl = gammaystore(INDEX_IJK,ispec_irreg)
        gammazl = gammazstore(INDEX_IJK,ispec_irreg)
        jacobianl = jacobianstore(INDEX_IJK,ispec_irreg)
        ! using fused array instead
        !xixl = deriv_mapping(1,INDEX_IJK,ispec_irreg)
        !xiyl = deriv_mapping(2,INDEX_IJK,ispec_irreg)
        !xizl = deriv_mapping(3,INDEX_IJK,ispec_irreg)
        !etaxl = deriv_mapping(4,INDEX_IJK,ispec_irreg)
        !etayl = deriv_mapping(5,INDEX_IJK,ispec_irreg)
        !etazl = deriv_mapping(6,INDEX_IJK,ispec_irreg)
        !gammaxl = deriv_mapping(7,INDEX_IJK,ispec_irreg)
        !gammayl = deriv_mapping(8,INDEX_IJK,ispec_irreg)
        !gammazl = deriv_mapping(9,INDEX_IJK,ispec_irreg)
        !jacobianl = deriv_mapping(10,INDEX_IJK,ispec_irreg)

        ! derivatives of potential
        dpotentialdxl = xixl*temp1(INDEX_IJK) + etaxl*temp2(INDEX_IJK) + gammaxl*temp3(INDEX_IJK)
        dpotentialdyl = xiyl*temp1(INDEX_IJK) + etayl*temp2(INDEX_IJK) + gammayl*temp3(INDEX_IJK)
        dpotentialdzl = xizl*temp1(INDEX_IJK) + etazl*temp2(INDEX_IJK) + gammazl*temp3(INDEX_IJK)

        temp1(INDEX_IJK) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
        temp2(INDEX_IJK) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
        temp3(INDEX_IJK) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
      ENDDO_LOOP_IJK

    else
      ! regular element
      jacobianl = jacobian_regular
      DO_LOOP_IJK
        ! reciprocal of density
        rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)
        ! for acoustic medium
        temp1(INDEX_IJK) = rho_invl * jacobianl * xix_regular * xix_regular * temp1(INDEX_IJK)
        temp2(INDEX_IJK) = rho_invl * jacobianl * xix_regular * xix_regular * temp2(INDEX_IJK)
        temp3(INDEX_IJK) = rho_invl * jacobianl * xix_regular * xix_regular * temp3(INDEX_IJK)
      ENDDO_LOOP_IJK
    endif

    ! second double-loop over GLL to compute all the terms along the x,y,z directions and assemble the contributions

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtemp1
    ! computes 2. matrix multiplication for newtemp2
    ! computes 3. matrix multiplication for newtemp3
    select case (NGLLX)
    case (5)
      call mxm5_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm5_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm5_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (6)
      call mxm6_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm6_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm6_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (7)
      call mxm7_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm7_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm7_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (8)
      call mxm8_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm8_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm8_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case default
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
            ! to be compatible with matrix versions from above
            newtemp1(i,j,k) = temp1l
            newtemp2(i,j,k) = temp2l
            newtemp3(i,j,k) = temp3l
            ! or instead, also add GLL integration weights
            !temp4(i,j,k) = - ( wgllwgll_yz(j,k)*temp1l + wgllwgll_xz(i,k)*temp2l + wgllwgll_xy(i,j)*temp3l )
          enddo
        enddo
      enddo
    end select

    ! sum contributions from each element to the global values
    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! also add GLL integration weights
          ! alternative 1:
          !potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
          !     - ( wgllwgll_yz(j,k)*newtemp1(i,j,k) &
          !       + wgllwgll_xz(i,k)*newtemp2(i,j,k) &
          !       + wgllwgll_xy(i,j)*newtemp3(i,j,k))
          ! or alternative 2: using (i,j,k) 3D weight arrays (faster)
!$OMP ATOMIC
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
               - ( wgllwgll_yz_3D(i,j,k)*newtemp1(i,j,k) &
                 + wgllwgll_xz_3D(i,j,k)*newtemp2(i,j,k) &
                 + wgllwgll_xy_3D(i,j,k)*newtemp3(i,j,k))
        enddo
      enddo
    enddo

  enddo ! end of loop over all spectral elements
!$OMP ENDDO
!$OMP END PARALLEL

  ! adds contributions for PML elements
  if (NSPEC_CPML > 0 .and. .not. backward_simulation) then
    call compute_forces_acoustic_PML(iphase, &
                                     potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                     backward_simulation)
  endif

  end subroutine compute_forces_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_PML(iphase, &
                                         potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                         backward_simulation)

! computes forces for acoustic PML elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,m1,m2

  use specfem_par, only: NGLOB_AB, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore, &
                         hprime_xx,hprime_xxT, &
                         hprimewgll_xx,hprimewgll_xxT, &
                         hprime_yy,hprime_zz, &
                         hprimewgll_yy,hprimewgll_zz, &
                         ibool, &
                         irregular_element_number,xix_regular

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D
  !or: use specfem_par, only: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  ! USE_DERIV_MAPPING_FUSION not used yet: use specfem_par, only: deriv_mapping

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  use specfem_par_acoustic, only: nspec_inner_acoustic,nspec_outer_acoustic, &
                                   phase_ispec_inner_acoustic

  use pml_par, only: is_CPML, spec_to_CPML, &
                     rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                     rmemory_potential_acoustic, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  integer, intent(in) :: iphase

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_dot_dot_acoustic

  ! CPML adjoint
  logical,intent(in) :: backward_simulation

  ! local variables

  ! note: declaring arrays in this subroutine here will allocate them generally on the stack
  !       (intel by default; not for gfortran though, it always uses heap memory).
  !       stack memory access is faster, thus please let these declarations here for local element arrays...

  ! arrays for elemental computations inside a given spectral element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem_old,chi_elem_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtemp1,newtemp2,newtemp3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1_old,temp2_old,temp3_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1_new,temp2_new,temp3_new

  ! CPML
  integer :: ispec_CPML
  ! arrays for elemental computations in compute_forces() for PML elements
  ! stores C-PML contribution to update the second derivative of the potential to the global mesh
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: potential_dot_dot_acoustic_CPML
  ! derivatives of potential with respect to x, y and z
  ! in computation potential_acoustic at "n" time step is used
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl
  ! in computation of PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old
  ! we replace potential_acoustic with potential_acoustic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old
  ! we replace potential_acoustic at "n" time step with
  ! we replace potential_acoustic with potential_acoustic_old with potential_acoustic_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new

  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: temp1l_old,temp2l_old,temp3l_old
  real(kind=CUSTOM_REAL) :: temp1l_new,temp2l_new,temp3l_new

  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  integer :: i,j,k,l,ispec,ispec_irreg,iglob,ispec_p,num_elements
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

  ! number of elements to loop
  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

! openmp solver
!$OMP PARALLEL if (num_elements > 100) &
!$OMP DEFAULT(NONE) &
!$OMP SHARED( &
!$OMP num_elements,ibool, &
!$OMP iphase,phase_ispec_inner_acoustic, &
!$OMP irregular_element_number,xix_regular, &
!$OMP potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
!$OMP PML_potential_acoustic_old,PML_potential_acoustic_new, &
!$OMP rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl,rmemory_potential_acoustic, &
!$OMP spec_to_CPML,is_CPML,backward_simulation, &
!$OMP xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore &
!$OMP ) &
!$OMP PRIVATE( &
!$OMP ispec_p,ispec,ispec_irreg,ispec_CPML,i,j,k,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#endif
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
!$OMP hp1,hp2,hp3, &
!$OMP temp1l,temp2l,temp3l, &
!$OMP temp1l_old,temp2l_old,temp3l_old, &
!$OMP temp1l_new,temp2l_new,temp3l_new, &
!$OMP chi_elem,chi_elem_old,chi_elem_new, &
!$OMP temp1,temp2,temp3, &
!$OMP temp1_old,temp2_old,temp3_old, &
!$OMP temp1_new,temp2_new,temp3_new, &
!$OMP newtemp1,newtemp2,newtemp3, &
!$OMP potential_dot_dot_acoustic_CPML, &
!$OMP PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
!$OMP PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
!$OMP PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new &
!$OMP ) &
!$OMP FIRSTPRIVATE( &
!$OMP hprime_xx,hprime_xxT, hprimewgll_xxT, hprimewgll_xx, &
!$OMP hprime_yy,hprime_zz,hprimewgll_yy,hprimewgll_zz, &
!$OMP wgllwgll_yz_3D,wgllwgll_xz_3D,wgllwgll_xy_3D &
!$OMP )

  ! loop over spectral elements
!$OMP DO
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! selects element contribution only for PML elements
    if (.not. (is_CPML(ispec) .and. .not. backward_simulation)) cycle

! all the loops below contain critical operations in which the solver spends 90% or so of its time; on modern machines,
! to get good performance it is thus crucial to make sure that most of these loops are fully vectorized by the compiler.
! In particular, PLEASE NEVER PUT ANY "IF" STATEMENT INSIDE THESE LOOPS;
! if you need to perform an "if" test, put it outside the loop
! and duplicate the content of the loop, using one for each case of the "if" test.

    ! gets value of the field inside the element and make it local
    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          chi_elem(i,j,k) = potential_acoustic(iglob)
        enddo
      enddo
    enddo

    ! additional old and new arrays
    ispec_CPML = spec_to_CPML(ispec)

    ! gets value of the field inside the element and make it local
    DO_LOOP_IJK
      chi_elem_old(INDEX_IJK) = PML_potential_acoustic_old(INDEX_IJK,ispec_CPML)
      chi_elem_new(INDEX_IJK) = PML_potential_acoustic_new(INDEX_IJK,ispec_CPML)
    ENDDO_LOOP_IJK

    ! derivative along x, y, z
    ! first double loop over GLL points to compute and store gradients

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for temp1
    ! computes 2. matrix multiplication for temp2
    ! computes 3. matrix multiplication for temp3
    select case (NGLLX)
    case (5)
      call mxm5_3comp_singleA(hprime_xx,m1,chi_elem,chi_elem_old,chi_elem_new,temp1,temp1_old,temp1_new,m2)
      call mxm5_3comp_3dmat_single(chi_elem,chi_elem_old,chi_elem_new,m1,hprime_xxT,m1,temp2,temp2_old,temp2_new,m1)
      call mxm5_3comp_singleB(chi_elem,chi_elem_old,chi_elem_new,m2,hprime_xxT,temp3,temp3_old,temp3_new,m1)
    case (6)
      call mxm6_3comp_singleA(hprime_xx,m1,chi_elem,chi_elem_old,chi_elem_new,temp1,temp1_old,temp1_new,m2)
      call mxm6_3comp_3dmat_single(chi_elem,chi_elem_old,chi_elem_new,m1,hprime_xxT,m1,temp2,temp2_old,temp2_new,m1)
      call mxm6_3comp_singleB(chi_elem,chi_elem_old,chi_elem_new,m2,hprime_xxT,temp3,temp3_old,temp3_new,m1)
    case (7)
      call mxm7_3comp_singleA(hprime_xx,m1,chi_elem,chi_elem_old,chi_elem_new,temp1,temp1_old,temp1_new,m2)
      call mxm7_3comp_3dmat_single(chi_elem,chi_elem_old,chi_elem_new,m1,hprime_xxT,m1,temp2,temp2_old,temp2_new,m1)
      call mxm7_3comp_singleB(chi_elem,chi_elem_old,chi_elem_new,m2,hprime_xxT,temp3,temp3_old,temp3_new,m1)
    case (8)
      call mxm8_3comp_singleA(hprime_xx,m1,chi_elem,chi_elem_old,chi_elem_new,temp1,temp1_old,temp1_new,m2)
      call mxm8_3comp_3dmat_single(chi_elem,chi_elem_old,chi_elem_new,m1,hprime_xxT,m1,temp2,temp2_old,temp2_new,m1)
      call mxm8_3comp_singleB(chi_elem,chi_elem_old,chi_elem_new,m2,hprime_xxT,temp3,temp3_old,temp3_new,m1)
    case default
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            ! derivative along x, y, z
            temp1l = 0._CUSTOM_REAL
            temp2l = 0._CUSTOM_REAL
            temp3l = 0._CUSTOM_REAL

            temp1l_old = 0._CUSTOM_REAL
            temp2l_old = 0._CUSTOM_REAL
            temp3l_old = 0._CUSTOM_REAL

            temp1l_new = 0._CUSTOM_REAL
            temp2l_new = 0._CUSTOM_REAL
            temp3l_new = 0._CUSTOM_REAL
            ! we can merge these loops because NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              hp1 = hprime_xx(i,l)
              temp1l     = temp1l     + chi_elem(l,j,k)*hp1
              temp1l_old = temp1l_old + chi_elem_old(l,j,k)*hp1
              temp1l_new = temp1l_new + chi_elem_new(l,j,k)*hp1

              hp2 = hprime_yy(j,l)
              temp2l     = temp2l     + chi_elem(i,l,k)*hp2
              temp2l_old = temp2l_old + chi_elem_old(i,l,k)*hp2
              temp2l_new = temp2l_new + chi_elem_new(i,l,k)*hp2

              hp3 = hprime_zz(k,l)
              temp3l     = temp3l     + chi_elem(i,j,l)*hp3
              temp3l_old = temp3l_old + chi_elem_old(i,j,l)*hp3
              temp3l_new = temp3l_new + chi_elem_new(i,j,l)*hp3
            enddo
            temp1(i,j,k) = temp1l
            temp2(i,j,k) = temp2l
            temp3(i,j,k) = temp3l

            temp1_old(i,j,k) = temp1l_old
            temp2_old(i,j,k) = temp2l_old
            temp3_old(i,j,k) = temp3l_old

            temp1_new(i,j,k) = temp1l_new
            temp2_new(i,j,k) = temp2l_new
            temp3_new(i,j,k) = temp3l_new
          enddo
        enddo
      enddo
    end select

    ispec_irreg = irregular_element_number(ispec)
    if (ispec_irreg /= 0) then
      ! irregular element
      DO_LOOP_IJK
        ! get derivatives of ux, uy and uz with respect to x, y and z
        ! single arrays make it difficult for hardware pre-fetching...
        xixl = xixstore(INDEX_IJK,ispec_irreg)
        xiyl = xiystore(INDEX_IJK,ispec_irreg)
        xizl = xizstore(INDEX_IJK,ispec_irreg)
        etaxl = etaxstore(INDEX_IJK,ispec_irreg)
        etayl = etaystore(INDEX_IJK,ispec_irreg)
        etazl = etazstore(INDEX_IJK,ispec_irreg)
        gammaxl = gammaxstore(INDEX_IJK,ispec_irreg)
        gammayl = gammaystore(INDEX_IJK,ispec_irreg)
        gammazl = gammazstore(INDEX_IJK,ispec_irreg)
        ! using fused array
        !xixl = deriv_mapping(1,INDEX_IJK,ispec_irreg)
        !xiyl = deriv_mapping(2,INDEX_IJK,ispec_irreg)
        !xizl = deriv_mapping(3,INDEX_IJK,ispec_irreg)
        !etaxl = deriv_mapping(4,INDEX_IJK,ispec_irreg)
        !etayl = deriv_mapping(5,INDEX_IJK,ispec_irreg)
        !etazl = deriv_mapping(6,INDEX_IJK,ispec_irreg)
        !gammaxl = deriv_mapping(7,INDEX_IJK,ispec_irreg)
        !gammayl = deriv_mapping(8,INDEX_IJK,ispec_irreg)
        !gammazl = deriv_mapping(9,INDEX_IJK,ispec_irreg)

        ! following is not needed, will be computed later in pml_compute_memory_variables_acoustic() routine...
        !jacobianl = jacobianstore(INDEX_IJK,ispec_irreg)
        !! reciprocal of density
        !rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)
        !! derivatives of potential
        !dpotentialdxl = xixl*temp1(INDEX_IJK) + etaxl*temp2(INDEX_IJK) + gammaxl*temp3(INDEX_IJK)
        !dpotentialdyl = xiyl*temp1(INDEX_IJK) + etayl*temp2(INDEX_IJK) + gammayl*temp3(INDEX_IJK)
        !dpotentialdzl = xizl*temp1(INDEX_IJK) + etazl*temp2(INDEX_IJK) + gammazl*temp3(INDEX_IJK)
        !! for acoustic medium
        !temp1(INDEX_IJK) = rho_invl * jacobianl * (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
        !temp2(INDEX_IJK) = rho_invl * jacobianl * (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
        !temp3(INDEX_IJK) = rho_invl * jacobianl * (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

        ! stores derivatives of ux, uy and uz with respect to x, y and z
        PML_dpotential_dxl(INDEX_IJK) = &
          xixl*temp1(INDEX_IJK) + etaxl*temp2(INDEX_IJK) + gammaxl*temp3(INDEX_IJK)
        PML_dpotential_dyl(INDEX_IJK) = &
          xiyl*temp1(INDEX_IJK) + etayl*temp2(INDEX_IJK) + gammayl*temp3(INDEX_IJK)
        PML_dpotential_dzl(INDEX_IJK) = &
          xizl*temp1(INDEX_IJK) + etazl*temp2(INDEX_IJK) + gammazl*temp3(INDEX_IJK)

        PML_dpotential_dxl_old(INDEX_IJK) = &
          xixl*temp1_old(INDEX_IJK) + etaxl*temp2_old(INDEX_IJK) + gammaxl*temp3_old(INDEX_IJK)
        PML_dpotential_dyl_old(INDEX_IJK) = &
          xiyl*temp1_old(INDEX_IJK) + etayl*temp2_old(INDEX_IJK) + gammayl*temp3_old(INDEX_IJK)
        PML_dpotential_dzl_old(INDEX_IJK) = &
          xizl*temp1_old(INDEX_IJK) + etazl*temp2_old(INDEX_IJK) + gammazl*temp3_old(INDEX_IJK)

        PML_dpotential_dxl_new(INDEX_IJK) = &
          xixl*temp1_new(INDEX_IJK) + etaxl*temp2_new(INDEX_IJK) + gammaxl*temp3_new(INDEX_IJK)
        PML_dpotential_dyl_new(INDEX_IJK) = &
          xiyl*temp1_new(INDEX_IJK) + etayl*temp2_new(INDEX_IJK) + gammayl*temp3_new(INDEX_IJK)
        PML_dpotential_dzl_new(INDEX_IJK) = &
          xizl*temp1_new(INDEX_IJK) + etazl*temp2_new(INDEX_IJK) + gammazl*temp3_new(INDEX_IJK)
      ENDDO_LOOP_IJK

    else
      ! regular element
      DO_LOOP_IJK
        ! following is not needed, will be computed later in pml_compute_memory_variables_acoustic() routine...
        !! reciprocal of density
        !jacobianl = jacobian_regular
        !rho_invl = 1.0_CUSTOM_REAL / rhostore(INDEX_IJK,ispec)
        !dpotentialdxl = xix_regular*temp1(INDEX_IJK)
        !dpotentialdyl = xix_regular*temp2(INDEX_IJK)
        !dpotentialdzl = xix_regular*temp3(INDEX_IJK)
        !! for acoustic medium
        !temp1(INDEX_IJK) = rho_invl * jacobianl * xix_regular * dpotentialdxl
        !temp2(INDEX_IJK) = rho_invl * jacobianl * xix_regular * dpotentialdyl
        !temp3(INDEX_IJK) = rho_invl * jacobianl * xix_regular * dpotentialdzl

        ! stores derivatives of ux, uy and uz with respect to x, y and z
        PML_dpotential_dxl(INDEX_IJK) = xix_regular*temp1(INDEX_IJK)
        PML_dpotential_dyl(INDEX_IJK) = xix_regular*temp2(INDEX_IJK)
        PML_dpotential_dzl(INDEX_IJK) = xix_regular*temp3(INDEX_IJK)

        PML_dpotential_dxl_old(INDEX_IJK) = xix_regular*temp1_old(INDEX_IJK)
        PML_dpotential_dyl_old(INDEX_IJK) = xix_regular*temp2_old(INDEX_IJK)
        PML_dpotential_dzl_old(INDEX_IJK) = xix_regular*temp3_old(INDEX_IJK)

        PML_dpotential_dxl_new(INDEX_IJK) = xix_regular*temp1_new(INDEX_IJK)
        PML_dpotential_dyl_new(INDEX_IJK) = xix_regular*temp2_new(INDEX_IJK)
        PML_dpotential_dzl_new(INDEX_IJK) = xix_regular*temp3_new(INDEX_IJK)

      ENDDO_LOOP_IJK
    endif

    ! sets C-PML elastic memory variables to compute stress sigma and form dot product with test vector
    !
    ! note: this routine also updates temp1,temp2,temp3 arrays to calculate the "default" contribution
    !       in the second double-loop, done after this if (is_CPML ..) case
    call pml_compute_memory_variables_acoustic(ispec,ispec_CPML, &
                                               temp1,temp2,temp3, &
                                               rmemory_dpotential_dxl,rmemory_dpotential_dyl,rmemory_dpotential_dzl, &
                                               PML_dpotential_dxl,PML_dpotential_dyl,PML_dpotential_dzl, &
                                               PML_dpotential_dxl_old,PML_dpotential_dyl_old,PML_dpotential_dzl_old, &
                                               PML_dpotential_dxl_new,PML_dpotential_dyl_new,PML_dpotential_dzl_new)

    ! calculates contribution from each C-PML element to update potential
    call pml_compute_accel_contribution_acoustic(ispec,ispec_CPML, &
                                                 potential_acoustic,potential_dot_acoustic, &
                                                 rmemory_potential_acoustic, &
                                                 potential_dot_dot_acoustic_CPML)

    ! add contributions from each element to the global values
    ! this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    ! will be added at the end with final contributions...
    !do k = 1,NGLLZ
    !  do j = 1,NGLLY
    !    do i = 1,NGLLX
    !      iglob = ibool(i,j,k,ispec)
    !!$OMP ATOMIC
    !      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_CPML(i,j,k)
    !    enddo
    !  enddo
    !enddo

    ! second double-loop over GLL to compute all the terms along the x,y,z directions and assemble the contributions

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtemp1
    ! computes 2. matrix multiplication for newtemp2
    ! computes 3. matrix multiplication for newtemp3
    select case (NGLLX)
    case (5)
      call mxm5_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm5_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm5_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (6)
      call mxm6_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm6_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm6_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (7)
      call mxm7_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm7_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm7_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case (8)
      call mxm8_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
      call mxm8_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
      call mxm8_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    case default
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
            ! to be compatible with matrix versions from above
            newtemp1(i,j,k) = temp1l
            newtemp2(i,j,k) = temp2l
            newtemp3(i,j,k) = temp3l
            ! or instead, also add GLL integration weights
            !temp4(i,j,k) = - ( wgllwgll_yz(j,k)*temp1l + wgllwgll_xz(i,k)*temp2l + wgllwgll_xy(i,j)*temp3l )
          enddo
        enddo
      enddo
    end select

    ! sum contributions from each element to the global values
    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! adds GLL integration weights
          ! alternative 1:
          !potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
          !     - ( wgllwgll_yz(j,k)*newtemp1(i,j,k) &
          !       + wgllwgll_xz(i,k)*newtemp2(i,j,k) &
          !       + wgllwgll_xy(i,j)*newtemp3(i,j,k))
          !
          ! or alternative 2: using (i,j,k) 3D weight arrays (faster)
!$OMP ATOMIC
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
               - ( wgllwgll_yz_3D(i,j,k)*newtemp1(i,j,k) &
                 + wgllwgll_xz_3D(i,j,k)*newtemp2(i,j,k) &
                 + wgllwgll_xy_3D(i,j,k)*newtemp3(i,j,k) &
                 + potential_dot_dot_acoustic_CPML(i,j,k))
        enddo
      enddo
    enddo

  enddo ! end of loop over all spectral elements
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_forces_acoustic_PML




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
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ IVDEP
!DIR$ SIMD
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

  !-------------

  subroutine mxm8_single(A,n1,B,C,n3)

! two-dimensional arrays (64,8)/(8,64)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,8),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(8,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j) &
              + A(i,6) * B(6,j) &
              + A(i,7) * B(7,j) &
              + A(i,8) * B(8,j)
    enddo
  enddo

  end subroutine mxm8_single

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
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ IVDEP
!DIR$ SIMD
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

  subroutine mxm8_3dmat_single(A,n1,B,n2,C,n3)

! three-dimensional arrays (8,8,8) for A and C

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,8,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(8,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j) &
                  + A(i,6,k) * B(6,j) &
                  + A(i,7,k) * B(7,j) &
                  + A(i,8,k) * B(8,j)
      enddo
    enddo
  enddo

  end subroutine mxm8_3dmat_single

