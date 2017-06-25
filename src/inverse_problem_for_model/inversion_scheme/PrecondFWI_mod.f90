module precond_mod

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use specfem_par, only: CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC_ADJOINT, myrank, &
                         ibool, xstore, ystore, zstore


 !---------------------------------------------------------------------------------------------------------------------------------

  use inverse_problem_par

  implicit none

contains

  subroutine SetPrecond(iter_inverse, inversion_param, current_gradient, fwi_precond)

    type(inver),                                               intent(in)    :: inversion_param
    integer,                                                   intent(in)    :: iter_inverse
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(inout) :: current_gradient, fwi_precond
    real(kind=CUSTOM_REAL)                                                   :: tapper, x,y,z
    integer                                                                  :: i,j,k,ispec, iglob

    if (inversion_param%use_tapper) then

       if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) '       iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) ' ADD TAPPER X:',  inversion_param%xmin_tapper,   inversion_param%xmax_tapper
          write(IIDD,*) ' ADD TAPPER Y:',  inversion_param%ymin_tapper,   inversion_param%ymax_tapper
          write(IIDD,*) ' ADD TAPPER Z:',  inversion_param%zmin_tapper,   inversion_param%zmax_tapper
          write(IIDD,*)
          write(IIDD,*)
       endif

       do ispec=1, NSPEC_ADJOINT
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX

                   iglob=ibool(i,j,k,ispec)

                   x=xstore(iglob)
                   y=ystore(iglob)
                   z=zstore(iglob)

                   tapper = 1.

                   if (x > inversion_param%xmin_tapper .and. x < inversion_param%xmax_tapper .and. &
                       y > inversion_param%ymin_tapper .and. y < inversion_param%ymax_tapper .and. &
                       z > inversion_param%zmin_tapper .and. z < inversion_param%zmax_tapper ) then
                      tapper = 1.
                   else
                      tapper = 0.
                   endif

                   current_gradient(i,j,k,ispec,:) = tapper *   current_gradient(i,j,k,ispec,:)

                enddo
             enddo
          enddo
       enddo


    endif

    !! only need to define the Z2 precond once.
    if (inversion_param%z2_precond .and. iter_inverse == 0) then
        if (DEBUG_MODE) then
          write(IIDD,*)
          write(IIDD,*) '       iteration FWI : ', iter_inverse
          write(IIDD,*)
          write(IIDD,*) '             define Z*Z Precond :'
          write(IIDD,*)
          write(IIDD,*)
       endif

       do ispec=1, NSPEC_ADJOINT
          do k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX

                   iglob=ibool(i,j,k,ispec)
                   z=zstore(iglob)
                   fwi_precond(i,j,k,ispec,:) = z*z

                enddo
             enddo
          enddo
       enddo


    endif


  end subroutine SetPrecond

end module precond_mod
