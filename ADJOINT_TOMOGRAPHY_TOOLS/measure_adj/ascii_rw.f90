module ascii_rw

  ! This module replaces calls to previous libraries
  ! for reading and writing ascii files.
  ! All operations are in double precision.

  implicit none

contains

  !-----------------------------------------------------------------------------

  subroutine dwascii(fname,data,npt,b,dt)

    character(len=*),intent(in) :: fname
    integer, intent(in) :: npt
    double precision, intent(in) :: data(*)
    double precision, intent(in) :: b,dt
    integer :: ios,i

    open(9,file=trim(fname),iostat=ios,status='unknown')
    if (ios /= 0) then
       print *,'Error opening ascii file to write: ',trim(fname)
       stop
    endif
    do i = 1,npt
       write(9,*) b + dt*(i-1), data(i)
    enddo
    close(9)

  end subroutine dwascii

  !-----------------------------------------------------------------------------

  subroutine drascii(fname,data,npt,b,dt)

    character(len=*), intent(in) :: fname
    integer, intent(out) :: npt
    double precision, intent(out) :: b,dt
    double precision :: data(*)
    double precision :: dtemp,ttemp
    double precision, dimension(:), allocatable :: ti
    integer :: ios,i

    ! get number of lines in the input file
    npt=0
    open(10,file=trim(fname),iostat=ios,status='old')
    if (ios /= 0) then
       print *,'Error opening ascii file to read: ',trim(fname)
       stop
    endif
    do while (1==1)
       read(10,*,iostat=ios) ttemp, dtemp
       if (ios /= 0) exit
       npt=npt+1
    enddo
    close(10)

    allocate(ti(npt))
    open(10,file=trim(fname),iostat=ios,status='old')
    do i=1,npt
       read(10,*,iostat=ios) ttemp, dtemp
       ti(i) = ttemp
       data(i) = dtemp
    enddo
    close(10)

    ! assign b and dt -- assume time increment is constant
    b = ti(1)
    dt = (ti(npt) - ti(1)) / (npt-1)

    deallocate(ti)

  end subroutine drascii

  !-----------------------------------------------------------------------------

!!$! these functions can be used directly in fortran
!!$
!!$subroutine dread_ascfile_f(name,t0,dt,n,data)
!!$
!!$  implicit none
!!$
!!$  character(len=*) :: name
!!$  real*8 :: t0, dt, data(*)
!!$  integer :: n
!!$
!!$  call dread_ascfile_c(trim(name)//char(0), t0, dt, n, data)
!!$
!!$end subroutine dread_ascfile_f
!!$
!!$
!!$subroutine dwrite_ascfile_f(name,t0,dt,n,data)
!!$
!!$  implicit none
!!$
!!$  character(len=*) :: name
!!$  real*8 :: t0, dt, data(*)
!!$  integer :: n
!!$
!!$  call dwrite_ascfile_c(trim(name)//char(0), t0, dt, n, data)
!!$
!!$end subroutine dwrite_ascfile_f

end module ascii_rw
