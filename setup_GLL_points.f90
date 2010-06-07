!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine setup_GLL_points()

  use specfem_par
  implicit none
  integer :: i,j

  if(myrank == 0) then
    write(IMAIN,*) '******************************************'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*) '******************************************'
    write(IMAIN,*)
  endif

! set up GLL points, weights and derivation matrices for reference element (between -1,1)
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
         hprime_xx,hprime_yy,hprime_zz, &
         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)

! define transpose of derivation matrix
  do j = 1,NGLLY
    do i = 1,NGLLX
      hprime_xxT(j,i) = hprime_xx(i,j)
      hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
    enddo
  enddo

! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX))
  allocate(hpxir(NGLLX))
  allocate(hetar(NGLLY))
  allocate(hpetar(NGLLY))
  allocate(hgammar(NGLLZ))
  allocate(hpgammar(NGLLZ))

! create name of database
  call create_name_database(prname,myrank,LOCAL_PATH)

  end subroutine