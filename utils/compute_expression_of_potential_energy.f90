
  program expression_total_energy

! Dimitri Komatitsch, CNRS Marseille, France, April 2013

! wrote a small program to avoid risking having typos in the expression if I typed it manually

  implicit none

  integer :: i,j

  print *,'potential_energy = potential_energy + &'

! expression in 3D
! 1 = x, 2 = y and 3 = z
  do i = 1,3
    do j = 1,3
      print *,'sigma_',i,j,' * epsilon_',i,j,' + '
    enddo
  enddo

  end program expression_total_energy

