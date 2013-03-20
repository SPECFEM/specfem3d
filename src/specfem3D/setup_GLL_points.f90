!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
  integer :: i,j,k,ier

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

! define a 3D extension in order to be able to force vectorization in the compute_forces routines
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        wgllwgll_yz_3D(i,j,k) = wgllwgll_yz(j,k)
        wgllwgll_xz_3D(i,j,k) = wgllwgll_xz(i,k)
        wgllwgll_xy_3D(i,j,k) = wgllwgll_xy(i,j)
      enddo
    enddo
  enddo

! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX), &
          hpxir(NGLLX), &
          hetar(NGLLY), &
          hpetar(NGLLY), &
          hgammar(NGLLZ), &
          hpgammar(NGLLZ),stat=ier)
  if( ier /= 0 ) stop 'error allocating arrays for interpolators'

! create name of database
  call create_name_database(prname,myrank,LOCAL_PATH)

  end subroutine setup_GLL_points

