
  program expression_total_energy

! Dimitri Komatitsch, CNRS Marseille, France, April 2013

! wrote a small program to avoid risking having typos in the expression if I typed it manually

  implicit none

  integer :: i,j

! compute potential energy 1/2 sigma_ij epsilon_ij

  print *,'potential_energy = potential_energy + integration_weight * 0.5 * ( &'

! expression in 3D
! 1 = x, 2 = y and 3 = z
  do i = 1,3
    do j = 1,3
      if(i == 3 .and. j == 3) then
        print *,'sigma_',i,j,' * epsilon_',i,j,')'
      else
        print *,'sigma_',i,j,' * epsilon_',i,j,' + '
      endif
    enddo
  enddo

  end program expression_total_energy

