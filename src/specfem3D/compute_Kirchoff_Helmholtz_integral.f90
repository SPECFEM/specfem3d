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

  ! compute Kirchoff-Helmholtz integral on whole mesh

  subroutine compute_Kirchoff_Helmholtz_integral()

  use constants

  use specfem_par

  implicit none

  ! local parameters
  double precision :: integr_vol, integr_boun
  double precision :: f_integrand(NGLOB_AB)

  double precision :: weight, jacobianl
  double precision :: x_coord, y_coord, z_coord

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  integer :: i,j,k,ispec,iglob,igll,iface

  character(len=MAX_STRING_LEN) :: outputname


!  double precision :: xval,yval,zval
!  double precision :: xval_squared,yval_squared,zval_squared
!  double precision :: distance_squared,distance_cubed, &
!                      three_over_distance_squared,one_over_distance_cubed,three_over_distance_fifth_power
!  double precision :: common_multiplying_factor,common_mult_times_one_over,common_mult_times_three_over


!-----------------------------------------------------------------------------

   open(unit=157,file='testtemp_for_integr_comput',status='unknown',action='write') !!  a supprimer ensuite

!!! --- lecture du deplacement et traction en fonction de iglob (donc allant de 1 a NGLOB_AB)
!!!
!!! --- exemple simple pour test : ---------
 
  integr_vol  = 0.d0
  integr_boun = 0.d0

  do iglob = 1, NGLOB_AB

    x_coord = xstore(iglob)
    y_coord = ystore(iglob)
    z_coord = zstore(iglob)

    f_integrand(iglob) = dcos(x_coord) * dcos(y_coord) * dcos(z_coord)
!!    f_integrand(iglob) = x_coord**2.d0 * y_coord * z_coord**3.d0
!!    f_integrand(iglob) = 1.d0

  enddo

  write(157,*) 'min et max des X !!!!! = ', minval(xstore), maxval(xstore)
  write(157,*) 'min et max des Y !!!!! = ', minval(ystore), maxval(ystore)
  write(157,*) 'min et max des Z !!!!! = ', minval(zstore), maxval(zstore)

!!! -----------------------------------------

  ! calculates volume of all elements in mesh
  do ispec = 1, NSPEC_AB

    ! print information about number of elements done so far

    if (myrank == 0 .and. (mod(ispec,NSPEC_DISPLAY_INTERVAL) == 0 .or. ispec == 1 .or. ispec == NSPEC_AB)) then
       write(IMAIN,*) 'for Kirchoff-Helmholtz integral, ', ispec,' elements computed out of ', NSPEC_AB

       ! write time stamp file to give information about progression of the computation of Kirchoff-Helmholtz integral
       write(outputname,"('/timestamp_K_H_integral_calculation_ispec',i7.7,'_out_of_',i7.7)") ispec, NSPEC_AB

       ! timestamp file output
       open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')
       write(IOUT,*) ispec,' elements done for Kirchoff-Helmholtz integral calculation out of ',NSPEC_AB
       close(unit=IOUT)

    endif

    ! main computation

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the Jacobian
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! do this in double precision for accuracy
          jacobianl = 1.d0 / dble(xixl*(etayl*gammazl-etazl*gammayl) &
                                - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                                + xizl*(etaxl*gammayl-etayl*gammaxl))

          if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianl <= ZERO) stop & 
                                             'error: negative Jacobian found in K-H integral calculation'

          iglob = ibool(i,j,k,ispec)

          integr_vol = integr_vol + ( weight * jacobianl * f_integrand(iglob) )

!          common_multiplying_factor = jacobianl * weight * rhostore(i,j,k,ispec)

        enddo
      enddo
    enddo
  enddo

  write(157,*) ' '
  write(157,*) 'num_free, num_abs', num_free_surface_faces,num_abs_boundary_faces
  write(157,*) ' '

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

  write(157,*) ' ---------------------------- '
  write(157,*) ' '
  write(157,*) 'resultat de lintegrale volumique !!!!! = ',  integr_vol
  write(157,*) 'resultat de lintegrale surfacique !!!!! = ', integr_boun

  close(157)

  end subroutine compute_Kirchoff_Helmholtz_integral
