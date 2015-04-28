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

  subroutine compute_integral_on_all_domain_surface_or_volume(f_integrand, integr_boun, integr_vol)

  use constants

  use specfem_par

  implicit none

  ! local parameters
  double precision :: integr_vol, integr_boun
  double precision :: f_integrand(NGLOB_AB)

  double precision :: weightloc, jacobianloc
!  double precision :: x_coord, y_coord, z_coord

  real(kind=CUSTOM_REAL) :: xixloc,xiyloc,xizloc,etaxloc,etayloc,etazloc,gammaxloc,gammayloc,gammazloc

  integer :: i,j,k,ispec,iglob,igll,iface

!  character(len=MAX_STRING_LEN) :: outputname

!-----------------------------------------------------------------------------
 
  integr_vol  = 0.d0
  integr_boun = 0.d0

  if ( (Surf_or_vol_integral == 2) .or. (Surf_or_vol_integral == 3) ) then

    ! calculates volume of all elements in mesh
    do ispec = 1, NSPEC_AB

      ! main computation
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            weightloc = wxgll(i)*wygll(j)*wzgll(k)

            ! compute the Jacobian
            xixloc    = xix(i,j,k,ispec)
            xiyloc    = xiy(i,j,k,ispec)
            xizloc    = xiz(i,j,k,ispec)
            etaxloc   = etax(i,j,k,ispec)
            etayloc   = etay(i,j,k,ispec)
            etazloc   = etaz(i,j,k,ispec)
            gammaxloc = gammax(i,j,k,ispec)
            gammayloc = gammay(i,j,k,ispec)
            gammazloc = gammaz(i,j,k,ispec)

            ! do this in double precision for accuracy
            jacobianloc = 1.d0 / dble(xixloc*(etayloc*gammazloc-etazloc*gammayloc) &
                                    - xiyloc*(etaxloc*gammazloc-etazloc*gammaxloc) &
                                    + xizloc*(etaxloc*gammayloc-etayloc*gammaxloc))

            if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianloc <= ZERO) stop & 
                                             'error: negative Jacobian found in volume integral calculation'

            iglob = ibool(i,j,k,ispec)

            integr_vol = integr_vol + ( weightloc * jacobianloc * f_integrand(iglob) )
  
          enddo
        enddo
      enddo
    enddo

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

          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

          integr_boun = integr_boun + ( free_surface_jacobian2Dw(igll,iface) * f_integrand(iglob) )

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

          integr_boun = integr_boun + ( abs_boundary_jacobian2Dw(igll,iface) * f_integrand(iglob) )

        enddo
      enddo
    endif

  endif

  end subroutine compute_integral_on_all_domain_surface_or_volume
