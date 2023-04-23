!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

module precond_mod

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, myrank, &
                         ibool, xstore, ystore, zstore


 !---------------------------------------------------------------------------------------------------------------------------------

  use inverse_problem_par

  implicit none

contains

  subroutine SetPrecond(inversion_param, current_gradient, hess_approxim, fwi_precond)

    type(inver),                                               intent(in)    :: inversion_param
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),              intent(inout) :: current_gradient, fwi_precond, hess_approxim
    real(kind=CUSTOM_REAL)                                                   :: taper, x,y,z
    real(kind=CUSTOM_REAL)                                                   :: a, z1, z2, dl
    real(kind=CUSTOM_REAL)                                                   :: nrme_coef_tmp, nrme_coef
    integer                                                                  :: i,j,k,ispec, iglob

    if (inversion_param%use_taper) then
       ! tapers model region

       ! debug output
       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) ' iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' ADD taper X:',  inversion_param%xmin_taper,   inversion_param%xmax_taper
          write(IIDD,*) ' ADD taper Y:',  inversion_param%ymin_taper,   inversion_param%ymax_taper
          write(IIDD,*) ' ADD taper Z:',  inversion_param%zmin_taper,   inversion_param%zmax_taper
          write(IIDD,*)
       endif

       do ispec = 1, NSPEC_ADJOINT
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX

                   iglob=ibool(i,j,k,ispec)

                   x=xstore(iglob)
                   y=ystore(iglob)
                   z=zstore(iglob)

!!$                   taperx = 1.
!!$                   tapery = 1.
!!$                   taperz = 1.
!!$
!!$   if (x < inversion_param%xmin_taper) taperx = cos_taper(inversion_param%xmin_taper, inversion_param%xmin_taper, x)
!!$   if (x > inversion_param%xmax_taper) taperx = cos_taper(inversion_param%xmax_taper, inversion_param%xmax, x)
!!$   if (y < inversion_param%ymin_taper) tapery = cos_taper(inversion_param%ymin_taper, inversion_param%xmin_taper, y)
!!$   if (y > inversion_param%ymax_taper) tapery = cos_taper(inversion_param%ymax_taper, inversion_param%xmax, y)
!!$   if (z < inversion_param%zmin_taper) taperz = cos_taper(inversion_param%zmin_taper, inversion_param%xmin_taper, z)
!!$   if (z > inversion_param%zmax_taper) taperz = cos_taper(inversion_param%zmax_taper, , z)

                   if (x > inversion_param%xmin_taper .and. x < inversion_param%xmax_taper .and. &
                       y > inversion_param%ymin_taper .and. y < inversion_param%ymax_taper .and. &
                       z > inversion_param%zmin_taper .and. z < inversion_param%zmax_taper ) then
                      taper = 1.
                   else
                      taper = 0.
                   endif

                   current_gradient(i,j,k,ispec,:) = taper *   current_gradient(i,j,k,ispec,:)

                enddo
             enddo
          enddo
       enddo
    endif

    !! only need to define the Z2 precond once.
    if (inversion_param%z2_precond .and. iter_inverse == 0) then
       ! depth Z-square preconditioner

       ! debug output
       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) ' iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' define Z*Z Precond :'
          write(IIDD,*)
       endif

       do ispec = 1, NSPEC_ADJOINT
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX

                   iglob = ibool(i,j,k,ispec)
                   z = zstore(iglob)
                   fwi_precond(i,j,k,ispec,:) = z*z

                enddo
             enddo
          enddo
       enddo
    endif

   if (inversion_param%z_precond .and. iter_inverse == 0) then
       ! depth Z preconditioner
       a  = inversion_param%aPrc
       z1 = inversion_param%zPrc1
       z2 = inversion_param%zPrc2
       dl = abs(z1 -z2)

       ! debug output
       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) ' iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' define Z Precond :'
          write(IIDD,*)
       endif

       do ispec = 1, NSPEC_ADJOINT
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX

                   iglob = ibool(i,j,k,ispec)
                   z = zstore(iglob)

                   if ( z >= z1) then
                       fwi_precond(i,j,k,ispec,:) = 1.e-8  !! small value to reduce perturbation at subsurface
                   else if ( z <= z1 .and. z >= z2) then
                       fwi_precond(i,j,k,ispec,:) = exp(- 0.5 * (a *( z - z1 - dl) / (0.5* dl )**2)  )
                   else
                        fwi_precond(i,j,k,ispec,:) = z / z2
                   endif
                enddo
             enddo
          enddo
       enddo
    endif

    if (inversion_param%shin_precond .and. iter_inverse == 0) then
       ! Shin preconditioner

       ! debug output
       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) ' iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' define Shin Precond :'
          write(IIDD,*)
       endif

       fwi_precond(:,:,:,:,:) = 1._CUSTOM_REAL / abs(hess_approxim(:,:,:,:,:))
    endif

    if (inversion_param%energy_precond .and. iter_inverse == 0) then
       ! energy preconditioner

       ! debug output
       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) ' iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' define energy Precond :'
          write(IIDD,*)
       endif

       ! normalisation of preconditionner
       nrme_coef_tmp = maxval(abs(hess_approxim(:,:,:,:,1)))
       call max_all_all_cr(nrme_coef_tmp, nrme_coef)

       do i = 1,inversion_param%NinvPar
          fwi_precond(:,:,:,:,i) = nrme_coef / abs(hess_approxim(:,:,:,:,1))
       enddo
    endif

  end subroutine SetPrecond

end module precond_mod
