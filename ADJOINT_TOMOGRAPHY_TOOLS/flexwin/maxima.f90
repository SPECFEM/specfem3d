!
! $Id: $
!
! Given an array of floating point numbers, returns an array of 
! integers containing the indexes of the maxima of the array
!


  subroutine find_maxima(x,nx,maxima,nindex,w_level)
  implicit none

  integer, parameter :: UP = 1
  integer, parameter :: DOWN = -1
  integer, parameter :: FLAT = 0

  double precision, dimension (*), intent(in) :: x
  integer, intent(in)  :: nx
  double precision, intent(in)     :: w_level
   
  integer, dimension (*), intent(out) :: maxima
  integer,intent(out)  :: nindex

  integer :: i,k, prev_dir, next_dir
  integer, external :: direction

  k=0 


! set up the starting previous direction 
  prev_dir=direction(x(1),x(2), w_level)


! iterate through the points in the array 
  do i=2, nx-1
    !the next direction 
    next_dir=direction(x(i),x(i+1),w_level)
    if(next_dir.eq.DOWN .and. prev_dir .eq. UP) then
      k=k+1
      maxima(k)=i
    endif
!   set up prev_dir for next iternation 
    prev_dir = next_dir
  enddo

  nindex=k

end subroutine find_maxima


  integer function direction(a,b,w_level)
  implicit none

  integer, parameter :: UP = 1
  integer, parameter :: DOWN = -1
  integer, parameter :: FLAT = 0

  double precision, parameter :: EPS = 1.0e-12

  double precision, intent(in) :: a,b,w_level

!  write(*,*) 'Direction entered'

  if ( a < w_level .or. b < w_level ) then 
    direction = FLAT
  elseif ( abs(a-b) < EPS) then
    direction = FLAT
  elseif (b > a) then
    direction = UP
  else
    direction = DOWN
  endif

!  write(*,*) 'Direction =', direction

  end function direction
 

 

! given an array, the index of a local maximum, and a water level, give
! the index of the closest points in the array that dip below the water
! level

  subroutine bracket_maximum(x, nx, imax , w_level, i1, i2) 
  implicit none
  double precision, dimension(*), intent(in) :: x
  integer, intent(in) :: nx , imax
  double precision, intent(in) :: w_level

  integer, intent(out)  :: i1, i2
  integer  :: i

  i1=1
  i2=nx

  do i = imax, 1, -1
    if(x(i)<w_level) then
      i1=i
      exit
    endif
  enddo

  do i = imax, nx, +1
    if(x(i)<w_level) then
      i2=i
      exit
    endif
  enddo

  !write(*,*)'DEBUG: bracket returns :', i1, i2, nx

  end subroutine

!!$! given an array and the index of two end points, calculate the 
!!$! area under the curve by simple trapezium integration
!!$
!!$  function area_under_curve(x,dx,i1,i2)
!!$  double precision, dimension(*) :: x
!!$  double precision :: dx
!!$  integer :: i1, i2
!!$
!!$  integer :: i
!!$  double precision :: area
!!$
!!$  area=0.0
!!$
!!$  do i = i1, i2-1
!!$    area=area + dx*(x(i)+x(i+1))/2.0
!!$  enddo
!!$
!!$  ! compilation warning for this line, but do we use this?
!!$  area_under_curve = area
!!$
!!$  end function

  subroutine find_minima(x,nx,minima,nindex,w_level)
  implicit none

  integer, parameter :: UP = 1
  integer, parameter :: DOWN = -1
  integer, parameter :: FLAT = 0

  double precision, dimension (*), intent(in) :: x
  integer, intent(in)  :: nx
  double precision, intent(in)    :: w_level
   
  integer, dimension (*), intent(out) :: minima
  integer, intent(out) :: nindex

  integer :: i,k, prev_dir, next_dir
  integer, external :: direction
   
  k=0 


! set up the starting previous direction 
  prev_dir=direction(x(1),x(2), w_level)

! iterate through the points in the array 
  do i=2, nx-1
    !the next direction 
    next_dir=direction(x(i),x(i+1),w_level)
    if(next_dir.eq.UP .and. prev_dir .eq. DOWN) then
      k=k+1
      minima(k)=i
    endif
!   set up prev_dir for next iternation 
    prev_dir = next_dir
  enddo

  nindex=k

end subroutine find_minima



