

! these functions can be used directly in fortran

subroutine dread_ascfile_f(name,t0,dt,n,data)

  implicit none

  character(len=*) :: name
  real*8 :: t0, dt, data(*)
  integer :: n

  call dread_ascfile_c(trim(name)//char(0), t0, dt, n, data)

end subroutine dread_ascfile_f


subroutine dwrite_ascfile_f(name,t0,dt,n,data)

  implicit none

  character(len=*) :: name
  real*8 :: t0, dt, data(*)
  integer :: n

  call dwrite_ascfile_c(trim(name)//char(0), t0, dt, n, data)

end subroutine dwrite_ascfile_f

