!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!
!
!    This file contains a robust implementation of linear interpolation.
!    It is based on an example by Arjen Markus
!    http://flibs.sourceforge.net/robust_interp.f90
!    and has been adapted for descending order arrays.
!
module interpolation

    use global_parameters, only: sp, dp
    implicit none

    type interpolation_data
        logical                         :: useable = .false.
        integer                         :: extrapolation
        real(kind=dp), allocatable      :: x(:)
        real(kind=dp), allocatable      :: y(:)
    end type interpolation_data

    integer, parameter                  :: extrapolation_none     = 0
    integer, parameter                  :: extrapolation_constant = 1
    integer, parameter                  :: extrapolation_linear   = 2

contains

function interpolation_object( x, y, extrapolation )
    type(interpolation_data)        :: interpolation_object
    real(kind=dp), intent(in)       :: x(:)
    real(kind=dp), intent(in)       :: y(:)
    integer, intent(in)             :: extrapolation

    integer                         :: i
    integer                         :: ierr
    integer                         :: n
    logical                         :: success

    interpolation_object%useable = .false.

    if ( allocated(interpolation_object%x) ) then
        deallocate(interpolation_object%x )
    endif
    if ( allocated(interpolation_object%y) ) then
        deallocate(interpolation_object%y )
    endif

    !
    ! Set the extrapolation method
    !
    interpolation_object%extrapolation = extrapolation_none
    if ( extrapolation == extrapolation_constant .or. &
         extrapolation == extrapolation_linear        ) then
        interpolation_object%extrapolation = extrapolation
    endif

    !
    ! Enough data? If not, simply return
    !
    if ( size(x) < 2 .or. size(y) < size(x) ) then
        print *, 'ERROR: interpolation_object: Not enough data'
        print *, 'X: ', x
        print *, 'Y: ', y
        return
    endif

    !
    ! Data sorted?
    !
    success = .true.

    do i = 2,size(x)
        if ( x(i) > x(i-1) ) then
            print *, 'ERROR: interpolation_object: data not sorted'
            success = .false.
            exit
        endif
    enddo

    if ( .not. success ) then
        return
    endif

    !
    ! Copy the data
    !
    n = size(x)
    allocate( interpolation_object%x(n), &
              interpolation_object%y(n), stat = ierr )

    if ( ierr /= 0 ) then
        return
    endif

    !
    ! We allow array y to be larger than x, so take care of that
    !
    interpolation_object%x(1:n) = x(1:n)
    interpolation_object%y(1:n) = y(1:n)

    interpolation_object%useable = .true.


end function interpolation_object

subroutine interpolate( object, xp, estimate, success )

    type(interpolation_data)        :: object
    real(kind=dp), intent(in)       :: xp
    real(kind=dp), intent(out)      :: estimate
    logical, intent(out)            :: success

    integer                         :: i
    integer                         :: idx
    integer                         :: nd
    real                            :: dx
    real(kind=dp)                   :: eps=1d-6

    estimate = 0.0
    success  = .false.

    if ( .not. object%useable ) then
        print *, 'Unusable interpolation object'
        return
    endif

    !
    ! Check extrapolation
    !
    nd = size(object%x)

    if ( object%extrapolation == extrapolation_none ) then
        if ( xp > object%x(1)*(1+eps)  ) then
           print *, 'interpolation: x out of range (too large) and no extrapolation chosen'
           print *, 'xp:', xp, ', x(1): ', object%x(1)
           return
        end if
        if ( xp < object%x(nd)*(1-eps) ) then
           print *, 'interpolation: x out of range (too small) and no extrapolation chosen'
           print *, 'xp:', xp, ', x(nd): ', object%x(nd)
           return
        end if
    endif
    if ( object%extrapolation == extrapolation_constant ) then
        if ( xp > object%x(1)  ) then
            estimate = object%y(1)
            success = .true.
            return
        endif
        if ( xp < object%x(nd) ) then
            estimate = object%y(nd)
            success = .true.
            return
        endif
    endif

    !
    ! Search for the interval that contains xp
    ! (Linear extrapolation is taken care of
    ! automatically)
    !
    idx = nd - 1

    do i = 2,nd - 1
        if ( xp > object%x(i) ) then
            idx = i - 1
            exit
        endif
    enddo

    dx = object%x(idx+1) - object%x(idx)

    if ( dx /= 0.0 ) then
        estimate = object%y(idx) + &
            (xp - object%x(idx)) * (object%y(idx+1) - object%y(idx)) / dx
    else
        estimate = 0.5 * (object%y(idx+1) + object%y(idx))
    endif

    success  = .true.

end subroutine interpolate

!real function simple_interpolate( x, y, xp )
!
!    real, dimension(:), intent(in) :: x
!    real, dimension(:), intent(in)  :: y
!    real, intent(in)                :: xp
!
!    integer                         :: i
!    integer                         :: idx
!
!    !
!    ! Search for the interval that contains xp
!    !
!    idx = size(x) - 1
!
!    do i = 2,size(x)-1
!        if ( xp < x(i) ) then
!            idx = i - 1
!            exit
!        endif
!    enddo
!
!    simple_interpolate = y(idx) + (xp - x(idx)) * (y(idx+1) - y(idx)) / &
!                      (x(idx+1) -  x(idx))
!
!end function simple_interpolate

end module interpolation
