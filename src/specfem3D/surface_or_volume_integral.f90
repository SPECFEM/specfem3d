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


  subroutine surface_or_volume_integral_on_whole_domain()

  use constants

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weightpt, jacobianpt

  double precision, dimension(3) :: integr_volloc, integr_bounloc
  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: f_integrand_bounloc
  double precision, dimension(3,NGLOB_AB) :: f_integrand_volloc

  real(kind=CUSTOM_REAL) :: xixpt,xiypt,xizpt,etaxpt,etaypt,etazpt,gammaxpt,gammaypt,gammazpt

  integer :: ier,iint,i,j,k,ispec,iglob,igll,iface,ispec_irreg

  ! NOTE : 'f_integrandloc' have to be defined at all 'iglob'
  ! (all the points of the surface, or all the point of the volume)

  if (RECIPROCITY_AND_KH_INTEGRAL) then

    allocate(f_integrand_KH(3,NGLLSQUARE*num_abs_boundary_faces), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1972')

    call integrand_for_computing_Kirchoff_Helmholtz_integral()

    do iint = 1,3
      f_integrand_bounloc(iint,:) = f_integrand_KH(iint,:)
    enddo

  endif

  do iint = 1,3
    integr_volloc(iint)  = 0.d0
    integr_bounloc(iint) = 0.d0
  enddo

  if ( (Surf_or_vol_integral == 2) .or. (Surf_or_vol_integral == 3) ) then

    allocate(integral_vol(3), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1973')

    ! calculates volume of all elements in mesh
    do ispec = 1, NSPEC_AB

      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg == 0) xixpt = dble(xix_regular)
      ! main computation
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            weightpt = wxgll(i)*wygll(j)*wzgll(k)

            if (ispec_irreg /= 0) then !irregular element
              ! compute the Jacobian
              xixpt    = xixstore(i,j,k,ispec_irreg)
              xiypt    = xiystore(i,j,k,ispec_irreg)
              xizpt    = xizstore(i,j,k,ispec_irreg)
              etaxpt   = etaxstore(i,j,k,ispec_irreg)
              etaypt   = etaystore(i,j,k,ispec_irreg)
              etazpt   = etazstore(i,j,k,ispec_irreg)
              gammaxpt = gammaxstore(i,j,k,ispec_irreg)
              gammaypt = gammaystore(i,j,k,ispec_irreg)
              gammazpt = gammazstore(i,j,k,ispec_irreg)

              ! do this in double precision for accuracy
              jacobianpt = 1.d0 / dble(xixpt*(etaypt*gammazpt-etazpt*gammaypt) &
                                      - xiypt*(etaxpt*gammazpt-etazpt*gammaxpt) &
                                      + xizpt*(etaxpt*gammaypt-etaypt*gammaxpt))

            else !regular element
              jacobianpt = 1.d0 / dble(xixpt*xixpt*xixpt)
            endif

            if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianpt <= ZERO) stop &
                                             'error: negative Jacobian found in volume integral calculation'

            iglob = ibool(i,j,k,ispec)

            integr_volloc(:)  = integr_volloc(:) + ( weightpt * jacobianpt * f_integrand_volloc(:,iglob) )

          enddo
        enddo
      enddo
    enddo

    call sum_all_1Darray_dp(integr_volloc, integral_vol, 3)

  endif

!
!---
!

  if ( (Surf_or_vol_integral == 1) .or. (Surf_or_vol_integral == 3) ) then

    allocate(integral_boun(3), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1974')

    ! calculates integral on all the surface of the whole domain

    if (num_free_surface_faces > 0) then

      do iface = 1, num_free_surface_faces

        ispec = free_surface_ispec(iface)

        do igll = 1, NGLLSQUARE

          ! gets local indices for GLL point
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

!!          integr_bounloc(:) = integr_bounloc(:) + (free_surface_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,iglob))
          integr_bounloc(:) = integr_bounloc(:) + (free_surface_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,igll*iface))

        enddo
      enddo
    endif

    if (num_abs_boundary_faces > 0) then

      do iface = 1,num_abs_boundary_faces

        ispec = abs_boundary_ispec(iface)

        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

!!          integr_bounloc(:) = integr_bounloc(:) + (abs_boundary_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,iglob))
          integr_bounloc(:) = integr_bounloc(:) + (abs_boundary_jacobian2Dw(igll,iface) * f_integrand_bounloc(:,igll*iface))

        enddo
      enddo
    endif

    call sum_all_1Darray_dp(integr_bounloc, integral_boun, 3)

  endif

  if (RECIPROCITY_AND_KH_INTEGRAL) deallocate(f_integrand_KH)

  end subroutine surface_or_volume_integral_on_whole_domain

!
!-------------------------------------------------------------------------------------------------
!

!--- CD CD : subroutine not validated yet

  subroutine integrand_for_computing_Kirchoff_Helmholtz_integral()

  use constants

  use specfem_par

  use shared_parameters

  implicit none

  ! local parameters
  integer :: itau

!!  integer :: i, iglob, iint

  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: convol_displsem_tractaxisem
  double precision, dimension(3,NGLLSQUARE*num_abs_boundary_faces) :: convol_tractsem_displaxisem

  convol_displsem_tractaxisem(:,:) = 0.d0
  convol_tractsem_displaxisem(:,:) = 0.d0

!---
!---

  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then

    if (old_DSM_coupling_from_Vadim) then

    else
      !! for 2D light version
    endif

  else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then

    call read_axisem_disp_file(num_abs_boundary_faces*NGLLSQUARE)
    call read_specfem_tract_file(NGLLSQUARE*num_abs_boundary_faces)
    call read_specfem_disp_file(NGLLSQUARE*num_abs_boundary_faces)

!!!-- ci-dessous : en cours de modification

    do itau = 1, it
      if ( it > itau .and. (it - itau <= NSTEP) ) then

        convol_displsem_tractaxisem(:,:) = convol_displsem_tractaxisem(:,:) + &
                                           Displ_specfem_time(:,:, itau) * Tract_axisem_time(:,:, it - itau) * dt

        convol_tractsem_displaxisem(:,:) = convol_tractsem_displaxisem(:,:) + &
                                           Tract_specfem_time(:,:, itau) * Displ_axisem_time(:,:, it - itau) * dt

      endif
    enddo

    f_integrand_KH(:,:) = convol_displsem_tractaxisem(:,:) - convol_tractsem_displaxisem(:,:)

  endif

!---
!---

  end subroutine integrand_for_computing_Kirchoff_Helmholtz_integral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_axisem_disp_file(nb)

  use constants

  use specfem_par, only: it, Displ_axisem_time, RECIPROCITY_AND_KH_INTEGRAL, IIN_displ_axisem

  implicit none

  ! argument
  integer nb

  ! local parameters
  real(kind=CUSTOM_REAL) :: Displ_axisem(3,nb)

  read(IIN_displ_axisem) Displ_axisem

  if (RECIPROCITY_AND_KH_INTEGRAL) Displ_axisem_time(:,:,it) = Displ_axisem(:,:)

  end subroutine read_axisem_disp_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_specfem_disp_file(nb)

  use constants

  use specfem_par, only: it, Displ_specfem_time, RECIPROCITY_AND_KH_INTEGRAL, num_abs_boundary_faces

  implicit none

  ! argument
  integer nb

  ! local parameters
  integer iface, igll
  real(kind=CUSTOM_REAL) :: Displ_specfem(3,nb)

  do iface = 1,num_abs_boundary_faces
    do igll = 1,NGLLSQUARE

      read(238) Displ_specfem(1,igll*iface), Displ_specfem(3,igll*iface), Displ_specfem(3,igll*iface)

    enddo
  enddo

  if (RECIPROCITY_AND_KH_INTEGRAL) Displ_specfem_time(:,:,it) = Displ_specfem(:,:)

  end subroutine read_specfem_disp_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_specfem_tract_file(nb)

  use constants

  use specfem_par, only: it, Tract_specfem_time, RECIPROCITY_AND_KH_INTEGRAL, num_abs_boundary_faces

  implicit none

  ! argument
  integer nb

  ! local parameters
  integer iface, igll
  real(kind=CUSTOM_REAL) :: Tract_specfem(3,nb)

  do iface = 1,num_abs_boundary_faces
    do igll = 1,NGLLSQUARE
      read(237) Tract_specfem(1,igll*iface), Tract_specfem(3,igll*iface), Tract_specfem(3,igll*iface)
    enddo
  enddo

  if (RECIPROCITY_AND_KH_INTEGRAL) Tract_specfem_time(:,:,it) = Tract_specfem(:,:)

  end subroutine read_specfem_tract_file
