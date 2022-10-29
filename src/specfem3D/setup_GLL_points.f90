!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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


  subroutine setup_GLL_points()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none
  integer :: i,j,k,inum

  ! outputs total element numbers
  call sum_all_i(count(ispec_is_acoustic(:)),inum)
  if (myrank == 0) then
    write(IMAIN,*) 'total acoustic elements    :',inum
  endif
  call sum_all_i(count(ispec_is_elastic(:)),inum)
  if (myrank == 0) then
    write(IMAIN,*) 'total elastic elements     :',inum
  endif
  call sum_all_i(count(ispec_is_poroelastic(:)),inum)
  if (myrank == 0) then
    write(IMAIN,*) 'total poroelastic elements :',inum
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks Courant criteria on mesh
  call check_mesh_resolution(NSPEC_AB,NGLOB_AB, &
                             ibool,xstore,ystore,zstore, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore,mustore,rhostore, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             DT,model_speed_max,min_resolved_period)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '******************************************'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*) '******************************************'
    write(IMAIN,*)
    call flush_IMAIN()
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
      hprime_yyT(j,i) = hprime_yy(i,j)
      hprime_zzT(j,i) = hprime_zz(i,j)

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

  ! create name of database
  call create_name_database(prname,myrank,LOCAL_PATH)

  end subroutine setup_GLL_points

