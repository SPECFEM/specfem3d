!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine setup_GLL_points()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic

  implicit none
  integer :: i,j,k,ier,inum

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
    call flush_IMAIN()
  endif

  ! checks Courant criteria on mesh
  if (ELASTIC_SIMULATION) then
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)

  else if (POROELASTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2421')
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2422')
    rho_vp = 0.0_CUSTOM_REAL
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution_poro(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
                                    DT,model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  else if (ACOUSTIC_SIMULATION) then
    allocate(rho_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2423')
    if (ier /= 0) stop 'error allocating array rho_vp'
    allocate(rho_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2424')
    if (ier /= 0) stop 'error allocating array rho_vs'
    rho_vp = sqrt( kappastore / rhostore ) * rhostore
    rho_vs = 0.0_CUSTOM_REAL
    call check_mesh_resolution(myrank,NSPEC_AB,NGLOB_AB, &
                               ibool,xstore,ystore,zstore, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               DT,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
    deallocate(rho_vp,rho_vs)
  endif

  ! user output
  if (myrank == 0) then
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

  ! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2425')
  allocate(hpxir(NGLLX),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2426')
  allocate(hetar(NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2427')
  allocate(hpetar(NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2428')
  allocate(hgammar(NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2429')
  allocate(hpgammar(NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 2430')
  if (ier /= 0) stop 'error allocating arrays for interpolators'

  ! create name of database
  call create_name_database(prname,myrank,LOCAL_PATH)

  end subroutine setup_GLL_points

