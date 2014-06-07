
  program convert_points_and_weights

!! DK DK Dimitri Komatitsch, CNRS, Marseille, France, March 2014: added support for Gauss integration instead of Simpson

!! DK DK run this code using:    ./xconvert_points_and_weights < input_points_and_weights_in_old_format.txt

!! DK DK the tables of Gauss-Legendre points and weights in the input file are taken from Pavel Holoborodko
!! DK DK at http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/

  implicit none

  integer, parameter :: ns_max = 128
  real(kind=8), dimension(ns_max) :: xi_Gauss,weight_Gauss

  character(len=35) line_read
  character(len=40) line_to_write

  integer :: i,ns,ipower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DK DK for ns = 32, 64, 128, 256
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ipower = 5,8

  ns = 2**ipower

! read a line of header
  read(*,fmt="(a35)") line_read
  write(*,*)
  write(*,*) line_read
  write(*,*)

! read the points
  do i = 1,ns/2
    read(*,*) line_read
!   add a double precision Fortran symbol to the constant value read
    line_to_write = line_read(1:len_trim(line_read)) // 'd0'
!   the Gauss points are anti-symmetric
    write(*,*) '  xi_Gauss(',i,') = ',line_to_write(1:len_trim(line_to_write))
    write(*,*) '  xi_Gauss(',i + ns/2,') = - ',line_to_write(1:len_trim(line_to_write))
  enddo

! read a line of header
  read(*,fmt="(a35)") line_read
  write(*,*)
  write(*,*) line_read
  write(*,*)

! read the weights
  do i = 1,ns/2
    read(*,*) line_read
!   add a double precision Fortran symbol to the constant value read
    line_to_write = line_read(1:len_trim(line_read)) // 'd0'
    write(*,*) '  weight_Gauss(',i,') = ',line_to_write(1:len_trim(line_to_write))
!   the Gauss weights are symmetric
    write(*,*) '  weight_Gauss(',i + ns/2,') = ',line_to_write(1:len_trim(line_to_write))
  enddo

  enddo

  end program convert_points_and_weights

