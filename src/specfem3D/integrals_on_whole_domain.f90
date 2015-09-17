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
!
!-------------------------------------------------------------------------------------------------
!

  subroutine surface_or_volume_integral_on_whole_domain()

  use constants

  use specfem_par

  use shared_parameters, only: COUPLE_WITH_EXTERNAL_CODE

  implicit none

  ! local parameters
  double precision :: weightpt, jacobianpt
  double precision :: integr_volloc, integr_bounloc

  double precision, dimension(NGLOB_AB) :: f_integrandloc

  real(kind=CUSTOM_REAL) :: xixpt,xiypt,xizpt,etaxpt,etaypt,etazpt,gammaxpt,gammaypt,gammazpt

  integer :: ier,i,j,k,ispec,iglob,igll,iface

!-----------------------------------------------------------------------------

  allocate(f_integrand_KH(3,NGLOB_AB), stat=ier) !! Maybe to modify in the case of Surf_or_vol_integral == 1

!---

  ! NOTE : 'f_integrand' have to be defined at all 'iglob'
  ! (all the points of the surface, or all the point of the volume)

  if (COUPLE_WITH_EXTERNAL_CODE) then 
!!    call integrand_for_computing_Kirchoff_Helmholtz_integral(it)
    f_integrandloc(:) = f_integrand_KH(1,NGLOB_AB)
  endif 

  integr_volloc  = 0.d0
  integr_bounloc = 0.d0

  if ( (Surf_or_vol_integral == 2) .or. (Surf_or_vol_integral == 3) ) then

    ! calculates volume of all elements in mesh
    do ispec = 1, NSPEC_AB

      ! main computation
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            weightpt = wxgll(i)*wygll(j)*wzgll(k)

            ! compute the Jacobian
            xixpt    = xix(i,j,k,ispec)
            xiypt    = xiy(i,j,k,ispec)
            xizpt    = xiz(i,j,k,ispec)
            etaxpt   = etax(i,j,k,ispec)
            etaypt   = etay(i,j,k,ispec)
            etazpt   = etaz(i,j,k,ispec)
            gammaxpt = gammax(i,j,k,ispec)
            gammaypt = gammay(i,j,k,ispec)
            gammazpt = gammaz(i,j,k,ispec)

            ! do this in double precision for accuracy
            jacobianpt = 1.d0 / dble(xixpt*(etaypt*gammazpt-etazpt*gammaypt) &
                                    - xiypt*(etaxpt*gammazpt-etazpt*gammaxpt) &
                                    + xizpt*(etaxpt*gammaypt-etaypt*gammaxpt))

            if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianpt <= ZERO) stop &
                                             'error: negative Jacobian found in volume integral calculation'

            iglob = ibool(i,j,k,ispec)

            integr_volloc = integr_volloc + ( weightpt * jacobianpt * f_integrandloc(iglob) )

          enddo
        enddo
      enddo
    enddo

    call sum_all_dp(integr_volloc,integral_vol)

  endif

!
!---
!

  if ( (Surf_or_vol_integral == 1) .or. (Surf_or_vol_integral == 3) ) then

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

          integr_bounloc = integr_bounloc + ( free_surface_jacobian2Dw(igll,iface) * f_integrandloc(iglob) )

        enddo
      enddo
    endif

    if (num_abs_boundary_faces > 0) then

      do iface=1,num_abs_boundary_faces

        ispec = abs_boundary_ispec(iface)

        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

          integr_bounloc = integr_bounloc + ( abs_boundary_jacobian2Dw(igll,iface) * f_integrandloc(iglob) )

        enddo
      enddo
    endif

    call sum_all_dp(integr_bounloc,integral_boun)

  endif

  deallocate(f_integrand_KH)

  end subroutine surface_or_volume_integral_on_whole_domain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine integrand_for_computing_Kirchoff_Helmholtz_integral(ittmp)

  use constants

  use specfem_par

  use shared_parameters

  implicit none

  integer :: ittmp

  ! local parameters
  integer :: ier 
!!  integer :: i, iglob

  double precision, dimension(3,NGLOB_AB) :: convol_displsem_tractaxisem
  double precision, dimension(3,NGLOB_AB) :: convol_tractsem_displaxisem

  double precision, dimension(3,NGLOB_AB) :: Tract_axisem_inv
  double precision, dimension(3,NGLOB_AB) :: Displ_axisem_inv

!---

  allocate(f_integrand_KH(3,NGLOB_AB), stat=ier) !! Could maybe be allocated just on the surface points

!---
!---

  if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then

    if (old_DSM_coupling_from_Vadim) then

    else
      !! for 2D light version
    endif

  elseif (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then

    call read_axisem_disp_file(Displ_axisem, num_abs_boundary_faces*NGLLSQUARE)

    write(*,*) 'Read OKKKKKK'

!!!    Routine de convolution ??????

!    Tract_axisem_inv(:,:) = Tract_axisem(:,:) !! Lire les tractions Axisem en temps inverse
!    Displ_axisem_inv(:,:) = Displ_axisem(:,:) !! Lire le deplacement Axisem en temps inverse

!    convol_displsem_tractaxisem(:,:) = convol_displsem_tractaxisem(:,:) + displ(:,:)*Tract_axisem_inv(:,:) !! *dt ou equivalent ???
!    convol_tractsem_displaxisem(:,:) = convol_tractsem_displaxisem(:,:) + Traction(:,:)*Displ_axisem_inv(:,:) 
    !! ==> Traction specfem ???? + *dt ou equivalent ???

    !! une fois arrive au temps final : 
!    if (ittmp == NSTEP) f_integrand_KH(:,:) = convol_displsem_tractaxisem(:,:) - convol_tractsem_displaxisem(:,:)

  endif

!---
!---

  end subroutine integrand_for_computing_Kirchoff_Helmholtz_integral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_axisem_disp_file(Displ_axisem, nb)

  use constants

  implicit none

  integer nb
  real(kind=CUSTOM_REAL) :: Displ_axisem(3,nb)

  read(IIN_displ_axisem) Displ_axisem

  end subroutine read_axisem_disp_file

