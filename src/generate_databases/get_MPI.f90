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

  subroutine get_MPI_interface(nglob,nspec,ibool)

! sets up the MPI interface for communication between partitions
  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,SMALLVAL_TOL,IMAIN

  use generate_databases_par, only: NGNOD,NPROC

  ! MPI interfaces
  use generate_databases_par, only: num_interfaces_ext_mesh,my_neighbors_ext_mesh, &
    nibool_interfaces_ext_mesh,max_interface_size_ext_mesh,ibool_interfaces_ext_mesh

  ! external mesh, element indexing
  use generate_databases_par, only: nelmnts_ext_mesh,elmnts_ext_mesh, &
    my_nelmnts_neighbors_ext_mesh, my_interfaces_ext_mesh, &
    nodes_coords_ext_mesh

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nglob,nspec
  ! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! local parameters
  double precision, dimension(:), allocatable :: xp,yp,zp
  double precision :: x_min,x_max,x_min_all,x_max_all
  double precision :: SMALLVALTOL

  integer, dimension(:), allocatable :: locval
  integer :: nibool_interfaces_ext_mesh_true

  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh, ninseg_ext_mesh
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2

  ! assembly test
  integer :: i,j,k,ispec,iglob,countval,inum,ier
  integer :: max_nibool_interfaces_ext_mesh,num_interface_points
  integer,dimension(:),allocatable :: test_flag_i
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr
  integer, dimension(:,:), allocatable :: ibool_interfaces_dummy

  ! define sort geometrical tolerance based upon typical size of the model
  ! (compare with get_global() which uses same tolerance)
  !
  ! note: the tolerance value is important to adapt according to the mesh dimensions
  !       having the tolerance too small, e.g., for double precision runs, will lead to numerical artifacts
  !       when different partitions have slightly different node positions due to numerical round-offs.
  !       we adapt the tolerance to the dimensions of the mesh (using the x-direction as a typical dimension value),
  !       to avoid round-off problems for other kind of meshes
  !       (e.g., UTM meshes have often dimensions in ~ 10 km, local meshes for engineering ~ 1 m, ultrasonic ~ 1mm range).
  !
  !       one way to avoid this mesh size dependence would be to normalize the coordinates as done in the global version.
  !       again, that would need a typical dimension of the mesh geometry used. to check for the future...
  !
  ! min/max values in x-direction
  x_min = minval(nodes_coords_ext_mesh(1,:))
  x_max = maxval(nodes_coords_ext_mesh(1,:))
  call min_all_all_dp(x_min,x_min_all)
  call max_all_all_dp(x_max,x_max_all)

  ! define geometrical tolerance based upon typical size of the model
  ! (needs to decrease by a factor 0.1 to fix travis check TESTID=2)
  SMALLVALTOL = SMALLVAL_TOL * dabs(x_max_all - x_min_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     number of interfaces        : ',num_interfaces_ext_mesh
    if (num_interfaces_ext_mesh > 0) then
      write(IMAIN,*) '     creating MPI indexing       : x min/max = ',sngl(x_min_all),'/',sngl(x_max_all)
      write(IMAIN,*) '                                   tolerance = ',SMALLVALTOL
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! gets global indices for points on MPI interfaces
  ! (defined by my_interfaces_ext_mesh) between different partitions
  ! and stores them in ibool_interfaces_ext_mesh & nibool_interfaces_ext_mesh
  ! (number of total points)
  call prepare_assemble_MPI(nelmnts_ext_mesh,elmnts_ext_mesh, &
                            ibool,nglob, &
                            num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
                            my_nelmnts_neighbors_ext_mesh, my_interfaces_ext_mesh, &
                            ibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,NGNOD )

  ! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  do iinterface = 1, num_interfaces_ext_mesh
    ! number of points on this interface
    num_interface_points = nibool_interfaces_ext_mesh(iinterface)

    ! allocates sorting arrays
    allocate(xp(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 891')
    if (ier /= 0) stop 'error allocating array xp'
    allocate(yp(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 892')
    if (ier /= 0) stop 'error allocating array yp'
    allocate(zp(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 893')
    if (ier /= 0) stop 'error allocating array zp'
    allocate(locval(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 894')
    if (ier /= 0) stop 'error allocating array locval'
    allocate(ifseg(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 895')
    if (ier /= 0) stop 'error allocating array ifseg'
    allocate(reorder_interface_ext_mesh(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 896')
    if (ier /= 0) stop 'error allocating array reorder_interface_ext_mesh'
    allocate(ninseg_ext_mesh(num_interface_points),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 897')
    if (ier /= 0) stop 'error allocating array ninseg_ext_mesh'
    xp(:) = 0.0d0; yp(:) = 0.0d0; zp(:) = 0.0d0
    locval(:) = 0; ifseg(:) = .false.; reorder_interface_ext_mesh(:) = 0; ninseg_ext_mesh(:) = 0

    ! gets x,y,z coordinates of global points on MPI interface
    do ilocnum = 1, num_interface_points
      iglob = ibool_interfaces_ext_mesh(ilocnum,iinterface)
      ! we will use double precision locations to find/sort MPI points
      xp(ilocnum) = dble(xstore_unique(iglob))
      yp(ilocnum) = dble(ystore_unique(iglob))
      zp(ilocnum) = dble(zstore_unique(iglob))
    enddo

    ! sorts (lexicographically?) ibool_interfaces_ext_mesh and updates value
    ! of total number of points nibool_interfaces_ext_mesh_true
    call sort_array_coordinates(num_interface_points,xp,yp,zp, &
                                ibool_interfaces_ext_mesh(1:num_interface_points,iinterface), &
                                reorder_interface_ext_mesh,locval,ifseg, &
                                nibool_interfaces_ext_mesh_true, &
                                ninseg_ext_mesh,SMALLVALTOL)

    ! checks that number of MPI points are still the same
    num_points1 = num_points1 + num_interface_points
    num_points2 = num_points2 + nibool_interfaces_ext_mesh_true
    if (num_points1 /= num_points2) then
      print *,'Error: sorting MPI interface points:',myrank
      print *,'   interface:',iinterface,num_points1,num_points2
      print *,'   interface points: x min/max',minval(xp(:)),maxval(xp(:))
      print *,'                     y min/max',minval(yp(:)),maxval(yp(:))
      print *,'                     z min/max',minval(zp(:)),maxval(zp(:))
      print *,'   sort tolerance = ',SMALLVALTOL,'/ mesh dimensions = ',dabs(x_max_all - x_min_all)
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
  allocate(test_flag_i(nglob),test_flag_cr(nglob),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 898')
  if (ier /= 0) stop 'error allocating array test_flag etc.'
  test_flag_i(:) = 0
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
          if (test_flag_i(iglob) == 0) countval = countval + 1

          ! sets identifier
          test_flag_i(iglob) = myrank + 1
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
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 899')
  if (ier /= 0) stop 'error allocating array ibool_interfaces_dummy'
  ibool_interfaces_dummy(:,:) = 0

  countval = 0
  do iinterface = 1, num_interfaces_ext_mesh
    num_interface_points = nibool_interfaces_ext_mesh(iinterface)
    ! updates counts
    ibool_interfaces_dummy(:,iinterface) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,iinterface)
    countval = countval + num_interface_points
  enddo
  call synchronize_all()

  call sum_all_i(countval,iglob)
  if (myrank == 0) then
    if (iglob /= ilocnum) call exit_mpi(myrank,'error total global MPI interface points')
  endif

  ! adds contributions from different partitions to flag arrays
  ! integer arrays
  call assemble_MPI_scalar_i_blocking(NPROC,nglob,test_flag_i, &
                                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                      nibool_interfaces_ext_mesh,ibool_interfaces_dummy, &
                                      my_neighbors_ext_mesh)
  ! CUSTOM_REAL arrays
  call assemble_MPI_scalar_blocking(NPROC,nglob,test_flag_cr, &
                                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh,ibool_interfaces_dummy, &
                                    my_neighbors_ext_mesh)

  ! checks number of interface points
  i = 0
  j = 0
  do iglob = 1,nglob
    ! only counts flags with MPI contributions
    if (test_flag_i(iglob) > myrank+1) i = i + 1
    if (test_flag_cr(iglob) > myrank+1.0) j = j + 1
  enddo
  call sum_all_i(i,inum)
  call sum_all_i(j,iglob)
  if (myrank == 0) then
    write(IMAIN,*) '     total assembled MPI interface points:',inum
    if (inum /= iglob .or. inum > ilocnum) call exit_mpi(myrank,'error MPI assembly')
  endif

  end subroutine get_MPI_interface

