!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

  subroutine get_MPI(myrank,nglob,nspec,ibool, &
                    nelmnts_ext_mesh,elmnts_ext_mesh, &
                    my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                    ibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh, &
                    num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                    my_neighbours_ext_mesh)

! sets up the MPI interface for communication between partitions
  use generate_databases_par, only: NPROC,NGNOD,NGLLX,NGLLY,NGLLZ,SMALLVAL_TOL,IMAIN
  use create_regions_mesh_ext_par

  implicit none

  integer :: myrank,nglob,nspec

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! external mesh, element indexing
  integer :: nelmnts_ext_mesh
  integer, dimension(NGNOD,nelmnts_ext_mesh) :: elmnts_ext_mesh

  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
    my_interfaces_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
    ibool_interfaces_ext_mesh


  !local parameters
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer, dimension(:), allocatable :: locval
  integer :: nibool_interfaces_ext_mesh_true

  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh, &
    ninseg_ext_mesh
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2

  ! assembly test
  integer :: i,j,k,ispec,iglob,countval,inum,ier
  integer :: max_nibool_interfaces_ext_mesh
  integer,dimension(:),allocatable :: test_flag
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr
  integer, dimension(:,:), allocatable :: ibool_interfaces_dummy

  ! gets global indices for points on MPI interfaces
  ! (defined by my_interfaces_ext_mesh) between different partitions
  ! and stores them in ibool_interfaces_ext_mesh & nibool_interfaces_ext_mesh
  ! (number of total points)
  call prepare_assemble_MPI( nelmnts_ext_mesh,elmnts_ext_mesh, &
                            ibool,nglob, &
                            num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
                            my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                            ibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,NGNOD )

  ! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  do iinterface = 1, num_interfaces_ext_mesh

    allocate(xp(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array xp'
    allocate(yp(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array yp'
    allocate(zp(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array zp'
    allocate(locval(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array locval'
    allocate(ifseg(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array ifseg'
    allocate(reorder_interface_ext_mesh(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array reorder_interface_ext_mesh'
    allocate(ninseg_ext_mesh(nibool_interfaces_ext_mesh(iinterface)),stat=ier)
    if (ier /= 0) stop 'error allocating array ninseg_ext_mesh'

    ! gets x,y,z coordinates of global points on MPI interface
    do ilocnum = 1, nibool_interfaces_ext_mesh(iinterface)
      xp(ilocnum) = xstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      yp(ilocnum) = ystore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      zp(ilocnum) = zstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
    enddo

    ! sorts (lexicographically?) ibool_interfaces_ext_mesh and updates value
    ! of total number of points nibool_interfaces_ext_mesh_true
    call sort_array_coordinates(nibool_interfaces_ext_mesh(iinterface),xp,yp,zp, &
         ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         reorder_interface_ext_mesh,locval,ifseg, &
         nibool_interfaces_ext_mesh_true, &
         ninseg_ext_mesh,SMALLVAL_TOL)

    ! checks that number of MPI points are still the same
    num_points1 = num_points1 + nibool_interfaces_ext_mesh(iinterface)
    num_points2 = num_points2 + nibool_interfaces_ext_mesh_true
    if (num_points1 /= num_points2) then
      write(*,*) 'error sorting MPI interface points:',myrank
      write(*,*) '   interface:',iinterface,num_points1,num_points2
      call exit_mpi(myrank,'error sorting MPI interface')
    endif

    ! cleanup temporary arrays
    deallocate(xp)
    deallocate(yp)
    deallocate(zp)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(reorder_interface_ext_mesh)
    deallocate(ninseg_ext_mesh)

  enddo

  ! outputs total number of MPI interface points
  call sum_all_i(num_points2,ilocnum)
  if (myrank == 0) then
    write(IMAIN,*) '     total MPI interface points: ',ilocnum
  endif

  ! checks with assembly of test fields
  allocate(test_flag(nglob),test_flag_cr(nglob),stat=ier)
  if (ier /= 0) stop 'error allocating array test_flag etc.'
  test_flag(:) = 0
  test_flag_cr(:) = 0._CUSTOM_REAL
  countval = 0
  do ispec = 1, nspec
    ! sets flags on global points
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,k,ispec)

          ! counts number of unique global points to set
          if (test_flag(iglob) == 0) countval = countval + 1

          ! sets identifier
          test_flag(iglob) = myrank + 1
          test_flag_cr(iglob) = myrank + 1.0
        enddo
      enddo
    enddo
  enddo
  call synchronize_all()

  ! collects contributions from different MPI partitions
  ! sets up MPI communications
  max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
  allocate(ibool_interfaces_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'error allocating array ibool_interfaces_dummy'

  countval = 0
  do iinterface = 1, num_interfaces_ext_mesh
     ibool_interfaces_dummy(:,iinterface) = &
      ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,iinterface)
     countval = countval + nibool_interfaces_ext_mesh(iinterface)
  enddo
  call synchronize_all()

  call sum_all_i(countval,iglob)
  if (myrank == 0) then
    if (iglob /= ilocnum) call exit_mpi(myrank,'error total global MPI interface points')
  endif

  ! adds contributions from different partitions to flag arrays
  ! integer arrays
  call assemble_MPI_scalar_i_blocking(NPROC,nglob,test_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_dummy, &
                        my_neighbours_ext_mesh)
  ! CUSTOM_REAL arrays
  call assemble_MPI_scalar_blocking(NPROC,nglob,test_flag_cr, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_dummy, &
                        my_neighbours_ext_mesh)

  ! checks number of interface points
  i = 0
  j = 0
  do iglob=1,nglob
    ! only counts flags with MPI contributions
    if (test_flag(iglob) > myrank+1) i = i + 1
    if (test_flag_cr(iglob) > myrank+1.0) j = j + 1
  enddo
  call sum_all_i(i,inum)
  call sum_all_i(j,iglob)
  if (myrank == 0) then
    write(IMAIN,*) '     total assembled MPI interface points:',inum
    if (inum /= iglob .or. inum > ilocnum) call exit_mpi(myrank,'error MPI assembly')
  endif

  end subroutine get_MPI

