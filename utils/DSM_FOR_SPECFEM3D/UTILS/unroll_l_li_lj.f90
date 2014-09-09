
 program unroll_l_li_lj

!! DK DK writes a simple code to unroll some loops

 implicit none

  integer :: l(6),li(6),lj(6)

  integer :: i,j,ij

  ij = 0

  do i=1,3
    ij = ij + 1
    li(ij) = i + 1
    lj(ij) = 2
    l(ij) = 3 * (i-1) + 1
    do j=i+1,3
      ij = ij + 1
      l(ij) = l(ij-1) + 4
      li(ij) = li(ij-1) + 1
      lj(ij) = lj(ij-1) + 1
    enddo
  enddo

  do i = 1,6
    print *,'  integer, parameter :: l',i,' = ',l(i)
  enddo
  print *

  do i = 1,6
    print *,'  integer, parameter :: li',i,' = ',li(i)
  enddo
  print *

  do i = 1,6
    print *,'  integer, parameter :: lj',i,' = ',lj(i)
  enddo
  print *

 end program unroll_l_li_lj

