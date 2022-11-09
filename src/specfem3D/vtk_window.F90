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


! VTK visualization

#ifdef VTK_VIS

module vtk_window_par

  ! VTK module
  implicit none

  ! vtk buffer array
  real,dimension(:),allocatable :: vtkdata
  logical,dimension(:),allocatable :: vtkmask

  ! multi-mpi processes, gather data arrays on main
  real,dimension(:),allocatable :: vtkdata_all
  integer,dimension(:),allocatable :: vtkdata_points_all
  integer,dimension(:),allocatable :: vtkdata_offset_all
  integer :: vtkdata_numpoints_all

  ! single source location
  real :: vtkdata_source_x,vtkdata_source_y,vtkdata_source_z
  ! receiver locations
  real,dimension(:),allocatable :: vtkdata_recv_x,vtkdata_recv_y,vtkdata_recv_z

end module vtk_window_par

!
!-------------------------------------------------------------------------------------------------
!


  subroutine vtk_window_prepare()

  use specfem_par
  use specfem_par_elastic
  use vtk_window_par

  implicit none

  ! local parameters
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "preparing VTK runtime visualization"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! turns on VTK visualization mode
  VTK_MODE = .true.

  ! checks: only works with elastic wavefield (uses veloc for display)
  if (.not. ELASTIC_SIMULATION) then
    call exit_MPI(myrank,'Error: VTK_VIS only implemented for elastic simulations')
  endif

  ! initializes VTK window
  if (myrank == 0 ) call initialize_vtkwindow(GPU_MODE)
  call synchronize_all()

  ! note: the following ordering is intended:
  !       - source position is needed for clipping the volume.
  !       - free surface extends is needed for setting receiver glyph radius

  ! adds source sphere
  call vtk_window_prepare_source()

  ! mask
  allocate(vtkmask(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1738')
  if (ier /= 0) stop 'Error allocating arrays'
  vtkmask(:) = .false.

  ! free surface
  if (VTK_SHOW_FREESURFACE) call vtk_window_prepare_freesurface()

  ! adds receiver spheres
  call vtk_window_prepare_receivers()

  ! volume data
  if (VTK_SHOW_VOLUME) call vtk_window_prepare_volume()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  window will update every NTSTEP_BETWEEN_FRAMES = ",NTSTEP_BETWEEN_FRAMES
    write(IMAIN,*)
    write(IMAIN,*) "  VTK visualization preparation done"
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine vtk_window_prepare

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_prepare_source()

  use specfem_par
  use vtk_window_par

  implicit none

  ! local parameters
  integer :: ispec

  ! source location
  integer :: isource,ia
  double precision :: shape3D(NGNOD)
  double precision :: xil,etal,gammal
  double precision :: xmesh,ymesh,zmesh
  real(kind=CUSTOM_REAL),dimension(NGNOD) :: xelm,yelm,zelm

  ! adds source
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) "  VTK source sphere:"
    call flush_IMAIN()
  endif

  ! only consider single source location
  ! in case of multiple sources, estimate the "midpoint" of all source locations
  vtkdata_source_x = 0.0
  vtkdata_source_y = 0.0
  vtkdata_source_z = 0.0

  do isource = 1,NSOURCES
    ! spectral element id
    ispec = ispec_selected_source(isource)

    ! gets element anchor nodes
    if (myrank == islice_selected_source(isource)) then
      ! find the coordinates of the anchor (eight corner) nodes of the element
      call eval_shape3D_element_anchors(xelm,yelm,zelm,ispec,ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
    endif

    ! main collects corner locations
    if (islice_selected_source(isource) /= 0) then
      if (myrank == 0) then
        call recvv_cr(xelm,NGNOD,islice_selected_source(isource),0)
        call recvv_cr(yelm,NGNOD,islice_selected_source(isource),0)
        call recvv_cr(zelm,NGNOD,islice_selected_source(isource),0)
      else if (myrank == islice_selected_source(isource)) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif

    if (myrank == 0) then
      ! get the 3-D shape functions
      xil = xi_source(isource)
      etal = eta_source(isource)
      gammal = gamma_source(isource)
      call eval_shape3D_single(shape3D,xil,etal,gammal,NGNOD)

      ! interpolates source locations
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0
      do ia = 1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! stores location for VTK visualization
      vtkdata_source_x = vtkdata_source_x + sngl(xmesh) / dble(NSOURCES)
      vtkdata_source_y = vtkdata_source_y + sngl(ymesh) / dble(NSOURCES)
      vtkdata_source_z = vtkdata_source_z + sngl(zmesh) / dble(NSOURCES)
    endif
  enddo

  ! prepares source glyph
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) "    sphere location x/y/z : ",vtkdata_source_x,vtkdata_source_y,vtkdata_source_z
    call flush_IMAIN()

    ! creates source sphere
    call prepare_vtksource(vtkdata_source_x,vtkdata_source_y,vtkdata_source_z)
  endif
  call synchronize_all()

  end subroutine vtk_window_prepare_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_prepare_receivers()

  use specfem_par
  use vtk_window_par

  implicit none

  ! local parameters
  integer :: ispec,ier

  ! locations
  integer :: irec,ia,istart,iend
  double precision :: shape3D(NGNOD)
  double precision :: xil,etal,gammal
  double precision :: xmesh,ymesh,zmesh
  real(kind=CUSTOM_REAL),dimension(NGNOD) :: xelm,yelm,zelm
  character(len=MAX_LENGTH_STATION_NAME*nrec+1) :: stations_string

  ! adds source
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) "  VTK receiver spheres:"
    write(IMAIN,*) "    number of receivers : ",nrec
    call flush_IMAIN()
  endif

  ! array to hold receiver locations
  allocate(vtkdata_recv_x(nrec),vtkdata_recv_y(nrec),vtkdata_recv_z(nrec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1739')
  if (ier /= 0) stop 'Error allocating receiver arrays'
  vtkdata_recv_x(:) = 0.0
  vtkdata_recv_y(:) = 0.0
  vtkdata_recv_z(:) = 0.0

  do irec = 1,nrec
    ! spectral element id
    ispec = ispec_selected_rec(irec)

    ! gets element anchor nodes
    if (myrank == islice_selected_rec(irec)) then
      ! find the coordinates of the anchor (eight corner) nodes of the element
      call eval_shape3D_element_anchors(xelm,yelm,zelm,ispec,ibool,xstore,ystore,zstore,NSPEC_AB,NGLOB_AB)
    endif

    ! main collects corner locations
    if (islice_selected_rec(irec) /= 0) then
      if (myrank == 0) then
        call recvv_cr(xelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(yelm,NGNOD,islice_selected_rec(irec),0)
        call recvv_cr(zelm,NGNOD,islice_selected_rec(irec),0)
      else if (myrank == islice_selected_rec(irec)) then
        call sendv_cr(xelm,NGNOD,0,0)
        call sendv_cr(yelm,NGNOD,0,0)
        call sendv_cr(zelm,NGNOD,0,0)
      endif
    endif

    if (myrank == 0) then
      ! get the 3-D shape functions
      xil = xi_receiver(irec)
      etal = eta_receiver(irec)
      gammal = gamma_receiver(irec)
      call eval_shape3D_single(shape3D,xil,etal,gammal,NGNOD)

      ! interpolates receiver locations
      xmesh = 0.0
      ymesh = 0.0
      zmesh = 0.0
      do ia=1,NGNOD
        xmesh = xmesh + shape3D(ia)*xelm(ia)
        ymesh = ymesh + shape3D(ia)*yelm(ia)
        zmesh = zmesh + shape3D(ia)*zelm(ia)
      enddo

      ! stores location for VTK visualization
      vtkdata_recv_x(irec) = sngl(xmesh)
      vtkdata_recv_y(irec) = sngl(ymesh)
      vtkdata_recv_z(irec) = sngl(zmesh)
    endif
  enddo

  ! note: this is a workaround to make it more compatible between different compilers/platforms
  ! concatenates strings
  do irec = 1,nrec
    istart = (irec-1)*MAX_LENGTH_STATION_NAME + 1
    iend = irec*MAX_LENGTH_STATION_NAME
    stations_string(istart:iend) = station_name(irec)(1:MAX_LENGTH_STATION_NAME)
  enddo
  ! adds null termination character at the end (because C loves this)
  stations_string(len(stations_string):len(stations_string)) = char(0)

  ! prepares source glyph
  if (myrank == 0) then
    ! creates receiver spheres
    call prepare_vtkreceivers(nrec,vtkdata_recv_x,vtkdata_recv_y,vtkdata_recv_z,MAX_LENGTH_STATION_NAME,stations_string)
  endif
  call synchronize_all()

  end subroutine vtk_window_prepare_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_prepare_freesurface()

  use specfem_par
  use vtk_window_par

  implicit none

  ! local parameters
  integer :: iface,igll,i,j,k,iglob,ispec,inum,ier
  integer :: id1,id2,id3,id4
  integer :: ii,jj
  integer :: ispec_start,ispec_end
  real(kind=CUSTOM_REAL),dimension(1):: dummy
  integer,dimension(1):: dummy_i

  ! free surface points
  integer :: free_np,free_nspec
  real, dimension(:),allocatable :: free_x,free_y,free_z
  integer, dimension(:,:),allocatable :: free_conn
  integer, dimension(:),allocatable :: free_perm
  ! gather arrays for multi-mpi simulations
  real, dimension(:),allocatable :: free_x_all,free_y_all,free_z_all
  integer, dimension(:,:),allocatable :: free_conn_all
  integer, dimension(:),allocatable :: free_conn_offset_all,free_conn_nspec_all
  integer, dimension(:),allocatable :: free_points_all,free_offset_all
  integer :: free_np_all,free_nspec_all

  ! checks if anything to do
  if (.not. VTK_SHOW_FREESURFACE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  VTK free surface:"
    write(IMAIN,*) "    free surface faces    : ",num_free_surface_faces
  endif

  ! counts global free surface points
  vtkmask(:) = .false.

  ! determines number of global points on surface
  do iface = 1,num_free_surface_faces
    ispec = free_surface_ispec(iface)
    do igll = 1, NGLLSQUARE
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)
      ! hi/low resolution
      if (VTK_USE_HIRES) then
        ! all free surface points
        ! coordinates
        iglob = ibool(i,j,k,ispec)
        vtkmask(iglob) = .true.
      else
        ! only corner points
        if ((i == 1 .or. i == NGLLX) .and. &
            (j == 1 .or. j == NGLLY) .and. &
            (k == 1 .or. k == NGLLZ)) then
          ! coordinates
          iglob = ibool(i,j,k,ispec)
          vtkmask(iglob) = .true.
        endif
      endif
    enddo
  enddo

  ! loads free surface into data
  free_np = count(vtkmask(:))

  ! user output
  if (myrank == 0) write(IMAIN,*) "    loading surface points: ",free_np

  allocate(free_x(free_np),free_y(free_np),free_z(free_np),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1740')
  if (ier /= 0) stop 'Error allocating arrays'

  ! permutation array
  allocate(free_perm(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1741')
  if (ier /= 0) stop 'Error allocating arrays'

  free_perm(:) = 0
  inum = 0
  do iglob = 1,NGLOB_AB
    if (vtkmask(iglob) .eqv. .true.) then
      inum = inum + 1
      free_x(inum) = xstore(iglob)
      free_y(inum) = ystore(iglob)
      free_z(inum) = zstore(iglob)
      ! stores permutation
      free_perm(iglob) = inum
    endif
  enddo
  if (inum /= free_np) stop 'Error free_np count in loading free surface points'

  ! hi/low resolution
  if (VTK_USE_HIRES) then
    ! point connectivity
    free_nspec = num_free_surface_faces*(NGLLX-1)*(NGLLY-1)

    allocate(free_conn(4,free_nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1742')
    if (ier /= 0) stop 'Error allocating arrays'

    inum = 0
    free_conn(:,:) = -1
    do iface = 1,num_free_surface_faces
      ispec = free_surface_ispec(iface)
      do jj = 1,NGLLZ-1
        do ii = 1,NGLLY-1
          ! indices of corner points
          igll = (ii-1)*NGLLX + jj
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)
          id1 = free_perm(ibool(i,j,k,ispec))
          ! indices of corner points
          igll = ii*NGLLX + jj
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)
          id2 = free_perm(ibool(i,j,k,ispec))
          ! indices of corner points
          igll = ii*NGLLX + jj+1
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)
          id3 = free_perm(ibool(i,j,k,ispec))
          ! indices of corner points
          igll = (ii-1)*NGLLX + jj+1
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)
          id4 = free_perm(ibool(i,j,k,ispec))
          ! note: indices for vtk start at 0
          inum = inum+1
          free_conn(1,inum) = id1 - 1
          free_conn(2,inum) = id2 - 1
          free_conn(3,inum) = id3 - 1
          free_conn(4,inum) = id4 - 1
        enddo
      enddo
    enddo
  else
    ! point connectivity
    free_nspec = num_free_surface_faces

    allocate(free_conn(4,free_nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1743')
    if (ier /= 0) stop 'Error allocating arrays'

    inum = 0
    free_conn(:,:) = -1
    do iface = 1,num_free_surface_faces
      ispec = free_surface_ispec(iface)
      ! indices of corner points
      ii = 1
      jj = 1
      igll = (ii-1)*NGLLX + jj
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)
      id1 = free_perm(ibool(i,j,k,ispec))
      ! indices of corner points
      ii = NGLLX-1
      jj = 1
      igll = ii*NGLLX + jj
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)
      id2 = free_perm(ibool(i,j,k,ispec))
      ! indices of corner points
      ii = NGLLX-1
      jj = NGLLY-1
      igll = ii*NGLLX + jj+1
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)
      id3 = free_perm(ibool(i,j,k,ispec))
      ! indices of corner points
      ii = 1
      jj = NGLLY-1
      igll = (ii-1)*NGLLX + jj+1
      i = free_surface_ijk(1,igll,iface)
      j = free_surface_ijk(2,igll,iface)
      k = free_surface_ijk(3,igll,iface)
      id4 = free_perm(ibool(i,j,k,ispec))
      ! note: indices for vtk start at 0
      inum = inum+1
      free_conn(1,inum) = id1 - 1
      free_conn(2,inum) = id2 - 1
      free_conn(3,inum) = id3 - 1
      free_conn(4,inum) = id4 - 1
    enddo
  endif
  if (minval(free_conn(:,:)) < 0) stop 'Error vtk free surface point connectivity'


  ! gathers data from all MPI processes
  if (NPROC > 1) then
    ! multiple MPI processes

    ! user output
    !if (myrank == 0) print *,"    gathering all MPI infos... "

    ! number of volume points for all partitions together
    call sum_all_i(free_np,free_np_all)
    if (myrank == 0) write(IMAIN,*) "    all freesurface points: ",free_np_all

    ! gathers point infos
    allocate(free_points_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1744')
    if (ier /= 0) stop 'Error allocating arrays'

    free_points_all(:) = 0
    call gather_all_singlei(free_np,free_points_all,NPROC)

    ! array offsets
    allocate(free_offset_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1745')
    if (ier /= 0) stop 'Error allocating arrays'

    free_offset_all(1) = 0
    do i = 2, NPROC
      free_offset_all(i) = sum(free_points_all(1:i-1))
    enddo

    ! number of volume elements
    call sum_all_i(free_nspec,free_nspec_all)
    if (myrank == 0) write(IMAIN,*) "    all freesurface elements: ",free_nspec_all

    ! freesurface elements
    allocate(free_conn_nspec_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1746')
    if (ier /= 0) stop 'Error allocating arrays'

    free_conn_nspec_all(:) = 0
    call gather_all_singlei(4*free_nspec,free_conn_nspec_all,NPROC)

    ! array offsets
    allocate(free_conn_offset_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1747')
    if (ier /= 0) stop 'Error allocating arrays'

    free_conn_offset_all(1) = 0
    do i = 2, NPROC
      free_conn_offset_all(i) = sum(free_conn_nspec_all(1:i-1))
    enddo

    ! global data arrays (only needed on main process)
    if (myrank == 0) then
      ! gather locations
      allocate(free_x_all(free_np_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1748')
      allocate(free_y_all(free_np_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1749')
      allocate(free_z_all(free_np_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1750')
      if (ier /= 0) stop 'Error allocating free_x_all,... arrays'

      free_x_all(:) = 0.0
      free_y_all(:) = 0.0
      free_z_all(:) = 0.0

      ! connectivity
      allocate(free_conn_all(4,free_nspec_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1751')
      if (ier /= 0) stop 'Error allocating free_conn_all array'
      free_conn_all(:,:) = 0
    endif

    if (myrank == 0) then
      ! locations
      !if (myrank == 0) print *,"    locations..."
      call gatherv_all_cr(free_x,free_np, &
                          free_x_all,free_points_all,free_offset_all, &
                          free_np_all,NPROC)
      call gatherv_all_cr(free_y,free_np, &
                          free_y_all,free_points_all,free_offset_all, &
                          free_np_all,NPROC)
      call gatherv_all_cr(free_z,free_np, &
                          free_z_all,free_points_all,free_offset_all, &
                          free_np_all,NPROC)

      ! connectivity
      !if (myrank == 0) print *,"    connectivity..."
      call gatherv_all_i(free_conn,4*free_nspec, &
                         free_conn_all,free_conn_nspec_all,free_conn_offset_all, &
                         free_nspec_all,NPROC)

      ! shifts connectivity ids for all additional slices
      do i = 2, NPROC
        ! divides by 4 to get nspec numbers
        ispec_start = free_conn_offset_all(i)/4 + 1
        ispec_end = free_conn_offset_all(i)/4 + free_conn_nspec_all(i)/4
        do ispec = ispec_start,ispec_end
          free_conn_all(:,ispec) = free_conn_all(:,ispec) + free_offset_all(i)
        enddo
      enddo

      !if (myrank == 0) print *,"    preparing vtk field..."

      ! adds free surface to vtk window
      call prepare_vtkfreesurface(free_np_all,free_x_all,free_y_all,free_z_all, &
                                  free_nspec_all,free_conn_all)

    else
      ! all other process just send data
      ! locations
      call gatherv_all_cr(free_x,free_np, &
                          dummy,free_points_all,free_offset_all, &
                          1,NPROC)
      call gatherv_all_cr(free_y,free_np, &
                          dummy,free_points_all,free_offset_all, &
                          1,NPROC)
      call gatherv_all_cr(free_z,free_np, &
                          dummy,free_points_all,free_offset_all, &
                          1,NPROC)
      ! connectivity
      call gatherv_all_i(free_conn,4*free_nspec, &
                          dummy_i,free_conn_nspec_all,free_conn_offset_all, &
                          1,NPROC)
    endif
  else
    ! serial run
    ! creates vtk freesurface actor
    call prepare_vtkfreesurface(free_np,free_x,free_y,free_z, &
                                free_nspec,free_conn)

  endif

  ! frees memory
  deallocate(free_x,free_y,free_z)
  deallocate(free_conn,free_perm)
  if (NPROC > 1) then
    deallocate(free_conn_nspec_all,free_conn_offset_all)
    deallocate(free_points_all,free_offset_all)
    if (myrank == 0) deallocate(free_x_all,free_y_all,free_z_all,free_conn_all)
  endif

  call synchronize_all()

  end subroutine vtk_window_prepare_freesurface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_prepare_volume()

  use specfem_par
  use vtk_window_par

  implicit none

  ! local parameters
  integer :: i,j,k,iglob,ispec,inum,ier
  integer :: id1,id2,id3,id4,id5,id6,id7,id8
  integer :: NIT_res

  ! volume points
  integer :: vol_np,vol_nspec
  real, dimension(:),allocatable :: vol_x,vol_y,vol_z
  integer, dimension(:,:),allocatable :: vol_conn
  integer, dimension(:),allocatable :: vol_perm
  ! gather arrays for multi-mpi simulations
  real, dimension(:),allocatable :: vol_x_all,vol_y_all,vol_z_all
  integer, dimension(:,:),allocatable :: vol_conn_all
  integer, dimension(:),allocatable :: vol_conn_offset_all,vol_conn_nspec_all
  integer :: vol_nspec_all,ispec_start,ispec_end
  real(kind=CUSTOM_REAL),dimension(1):: dummy
  integer,dimension(1):: dummy_i

  ! checks if anything to do
  if (.not. VTK_SHOW_VOLUME) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  VTK volume:"
    write(IMAIN,*) "    spectral elements    : ",NSPEC_AB
    call flush_IMAIN()
  endif

  ! hi/low resolution increment
  if (VTK_USE_HIRES) then
    NIT_res = 1
  else
    NIT_res = NGLLX - 1
  endif

  ! sets new point mask
  vtkmask(:) = .false.
  do ispec = 1,NSPEC_AB
    ! hi/low resolution
    ! loops only over points
    do k = 1,NGLLZ,NIT_res
      do j = 1,NGLLY,NIT_res
        do i = 1,NGLLX,NIT_res
          iglob = ibool(i,j,k,ispec)
          ! sets mask
          vtkmask(iglob) = .true.
        enddo
      enddo
    enddo
  enddo
  vol_np = count(vtkmask(:))

  ! loads volume data arrays
  if (myrank == 0) write(IMAIN,*) "    loading volume points: ",vol_np

  allocate(vol_x(vol_np),vol_y(vol_np),vol_z(vol_np),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1752')
  if (ier /= 0) stop 'Error allocating arrays'

  ! permutation array
  allocate(vol_perm(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1753')
  if (ier /= 0) stop 'Error allocating arrays'

  vol_perm(:) = 0
  inum = 0
  do iglob = 1,NGLOB_AB
    if (vtkmask(iglob) .eqv. .true.) then
      inum = inum + 1
      vol_x(inum) = xstore(iglob)
      vol_y(inum) = ystore(iglob)
      vol_z(inum) = zstore(iglob)
      ! stores permutation
      vol_perm(iglob) = inum
    endif
  enddo

  ! hi/low resolution
  if (VTK_USE_HIRES) then
    ! point connectivity
    vol_nspec = NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)

    allocate(vol_conn(8,vol_nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1754')
    if (ier /= 0) stop 'Error allocating arrays'

    inum = 0
    vol_conn(:,:) = -1
    do ispec = 1,NSPEC_AB
      do k = 1, NGLLZ-1
        do j = 1, NGLLY-1
          do i = 1, NGLLX-1
            ! indices of corner points
            id1 = vol_perm(ibool(i,j,k,ispec))
            id2 = vol_perm(ibool(i+1,j,k,ispec))
            id3 = vol_perm(ibool(i+1,j+1,k,ispec))
            id4 = vol_perm(ibool(i,j+1,k,ispec))

            id5 = vol_perm(ibool(i,j,k+1,ispec))
            id6 = vol_perm(ibool(i+1,j,k+1,ispec))
            id7 = vol_perm(ibool(i+1,j+1,k+1,ispec))
            id8 = vol_perm(ibool(i,j+1,k+1,ispec))

            ! note: indices for vtk start at 0
            inum = inum+1
            vol_conn(1,inum) = id1 - 1
            vol_conn(2,inum) = id2 - 1
            vol_conn(3,inum) = id3 - 1
            vol_conn(4,inum) = id4 - 1
            vol_conn(5,inum) = id5 - 1
            vol_conn(6,inum) = id6 - 1
            vol_conn(7,inum) = id7 - 1
            vol_conn(8,inum) = id8 - 1
          enddo
        enddo
      enddo
    enddo
  else
    ! point connectivity
    vol_nspec = NSPEC_AB

    allocate(vol_conn(8,vol_nspec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1755')
    if (ier /= 0) stop 'Error allocating arrays'

    vol_conn(:,:) = -1
    do ispec = 1,NSPEC_AB
      ! indices of corner points
      id1 = vol_perm(ibool(1,1,1,ispec))
      id2 = vol_perm(ibool(NGLLX,1,1,ispec))
      id3 = vol_perm(ibool(NGLLX,NGLLY,1,ispec))
      id4 = vol_perm(ibool(1,NGLLY,1,ispec))

      id5 = vol_perm(ibool(1,1,NGLLZ,ispec))
      id6 = vol_perm(ibool(NGLLX,1,NGLLZ,ispec))
      id7 = vol_perm(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      id8 = vol_perm(ibool(1,NGLLY,NGLLZ,ispec))

      ! note: indices for vtk start at 0
      vol_conn(1,ispec) = id1 - 1
      vol_conn(2,ispec) = id2 - 1
      vol_conn(3,ispec) = id3 - 1
      vol_conn(4,ispec) = id4 - 1
      vol_conn(5,ispec) = id5 - 1
      vol_conn(6,ispec) = id6 - 1
      vol_conn(7,ispec) = id7 - 1
      vol_conn(8,ispec) = id8 - 1
    enddo
  endif
  if (minval(vol_conn(:,:)) < 0) stop 'Error vtk volume point connectivity'

  ! allocates local data array
  allocate(vtkdata(vol_np),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1756')
  if (ier /= 0) stop 'Error allocating arrays'

  vtkdata(:) = 0.0

  ! gathers data from all MPI processes
  if (NPROC > 1) then
    ! user output
    !if (myrank == 0) print *,"    gathering all MPI infos... "

    ! number of volume points for all partitions together
    call sum_all_i(vol_np,vtkdata_numpoints_all)
    if (myrank == 0) write(IMAIN,*) "    all volume points    : ",vtkdata_numpoints_all

    ! gathers point infos
    allocate(vtkdata_points_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1757')
    if (ier /= 0) stop 'Error allocating arrays'

    vtkdata_points_all(:) = 0
    call gather_all_singlei(vol_np,vtkdata_points_all,NPROC)

    ! array offsets
    allocate(vtkdata_offset_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1758')
    if (ier /= 0) stop 'Error allocating arrays'

    vtkdata_offset_all(1) = 0
    do i = 2, NPROC
      vtkdata_offset_all(i) = sum(vtkdata_points_all(1:i-1))
    enddo

    ! number of volume elements
    call sum_all_i(vol_nspec,vol_nspec_all)
    if (myrank == 0) write(IMAIN,*) "    all volume elements  : ",vol_nspec_all

    ! volume elements
    allocate(vol_conn_nspec_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1759')
    if (ier /= 0) stop 'Error allocating arrays'

    vol_conn_nspec_all(:) = 0
    call gather_all_singlei(8*vol_nspec,vol_conn_nspec_all,NPROC)

    ! array offsets
    allocate(vol_conn_offset_all(NPROC),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1760')
    if (ier /= 0) stop 'Error allocating arrays'

    vol_conn_offset_all(1) = 0
    do i = 2, NPROC
      vol_conn_offset_all(i) = sum(vol_conn_nspec_all(1:i-1))
    enddo

    ! global data arrays (only needed on main process)
    if (myrank == 0) then
      ! point data
      allocate(vtkdata_all(vtkdata_numpoints_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1761')
      if (ier /= 0) stop 'Error allocating vtkdata_all array'

      vtkdata_all(:) = 0.0

      ! gather locations
      allocate(vol_x_all(vtkdata_numpoints_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1762')
      allocate(vol_y_all(vtkdata_numpoints_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1763')
      allocate(vol_z_all(vtkdata_numpoints_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1764')
      if (ier /= 0) stop 'Error allocating vol_x_all,... arrays'

      vol_x_all(:) = 0.0
      vol_y_all(:) = 0.0
      vol_z_all(:) = 0.0

      ! connectivity
      allocate(vol_conn_all(8,vol_nspec_all),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 1765')
      if (ier /= 0) stop 'Error allocating vol_conn_all array'

      vol_conn_all(:,:) = 0

    endif

    if (myrank == 0) then
      ! locations
      !if (myrank == 0) print *,"    locations..."
      call gatherv_all_cr(vol_x,vol_np, &
                          vol_x_all,vtkdata_points_all,vtkdata_offset_all, &
                          vtkdata_numpoints_all,NPROC)
      call gatherv_all_cr(vol_y,vol_np, &
                          vol_y_all,vtkdata_points_all,vtkdata_offset_all, &
                          vtkdata_numpoints_all,NPROC)
      call gatherv_all_cr(vol_z,vol_np, &
                          vol_z_all,vtkdata_points_all,vtkdata_offset_all, &
                          vtkdata_numpoints_all,NPROC)

      ! connectivity
      !if (myrank == 0) print *,"    connectivity..."
      call gatherv_all_i(vol_conn,8*vol_nspec, &
                         vol_conn_all,vol_conn_nspec_all,vol_conn_offset_all, &
                         vol_nspec_all,NPROC)

      ! shifts connectivity ids for all additional slices
      do i = 2, NPROC
        ! divides by 8 to get nspec numbers
        ispec_start = vol_conn_offset_all(i)/8 + 1
        ispec_end = vol_conn_offset_all(i)/8 + vol_conn_nspec_all(i)/8
        do ispec = ispec_start,ispec_end
          vol_conn_all(:,ispec) = vol_conn_all(:,ispec) + vtkdata_offset_all(i)
        enddo
      enddo

      !if (myrank == 0) print *,"    preparing vtk field..."

      ! adds total volume wavefield to vtk window
      call prepare_vtkfield(vtkdata_numpoints_all,vol_x_all,vol_y_all,vol_z_all, &
                            vol_nspec_all,vol_conn_all)

    else
      ! all other process just send data
      ! locations
      call gatherv_all_cr(vol_x,vol_np, &
                          dummy,vtkdata_points_all,vtkdata_offset_all, &
                          1,NPROC)
      call gatherv_all_cr(vol_y,vol_np, &
                          dummy,vtkdata_points_all,vtkdata_offset_all, &
                          1,NPROC)
      call gatherv_all_cr(vol_z,vol_np, &
                          dummy,vtkdata_points_all,vtkdata_offset_all, &
                          1,NPROC)
      ! connectivity
      call gatherv_all_i(vol_conn,8*vol_nspec, &
                          dummy_i,vol_conn_nspec_all,vol_conn_offset_all, &
                          1,NPROC)

    endif

  else
    ! serial run
    !if (myrank == 0) print *,"    preparing vtk field..."

    ! adds volume wavefield to vtk window
    call prepare_vtkfield(vol_np,vol_x,vol_y,vol_z,vol_nspec,vol_conn)
  endif

  ! frees memory
  deallocate(vol_x,vol_y,vol_z)
  deallocate(vol_conn,vol_perm)
  if (NPROC > 1) then
    deallocate(vol_conn_nspec_all,vol_conn_offset_all)
    if (myrank == 0) deallocate(vol_x_all,vol_y_all,vol_z_all,vol_conn_all)
  endif
  call synchronize_all()

  end subroutine vtk_window_prepare_volume

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_update()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_movie
  use vtk_window_par

  implicit none

  ! local parameters
  real :: currenttime
  integer :: iglob,inum
  real(kind=CUSTOM_REAL),dimension(1):: dummy

  ! vtk rendering at frame interval
  if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

    ! user output
    !if (myrank == 0) print *,"  VTK rendering..."

    ! updates time
    currenttime = sngl((it-1)*DT-t0)

    ! transfers fields from GPU to host
    if (GPU_MODE) then
      !if (myrank == 0) print *,"  vtk: transfering velocity from gpu"
      call transfer_veloc_from_device(NDIM*NGLOB_AB,veloc,Mesh_pointer)
    endif

    ! updates wavefield
    !if (myrank == 0) print *,"  vtk: it = ",it," - updating velocity field"
    inum = 0
    vtkdata(:) = 0.0
    do iglob = 1,NGLOB_AB
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! stores norm of velocity vector
        vtkdata(inum) = sqrt(veloc(1,iglob)**2 + veloc(2,iglob)**2 + veloc(3,iglob)**2)
      endif
    enddo

    ! updates for multiple MPI process
    if (NPROC > 1) then
      if (myrank == 0) then
        ! gather data
        call gatherv_all_cr(vtkdata,size(vtkdata), &
                            vtkdata_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)
      else
        ! all other process just send data
        call gatherv_all_cr(vtkdata,size(vtkdata), &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
      endif

      ! updates vtk window
      if (myrank == 0) then
        call visualize_vtkdata(it,currenttime,vtkdata_all)
      endif

    else
      ! serial run
      ! updates vtk window
      call visualize_vtkdata(it,currenttime,vtkdata)
    endif

  endif

  end subroutine vtk_window_update

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vtk_window_cleanup(do_restart)

  use specfem_par
  use specfem_par_elastic
  use vtk_window_par

  implicit none

  logical,intent(inout) :: do_restart

  ! closes/cleans up vtk window
  if (myrank == 0) call finish_vtkwindow(do_restart)
  call bcast_all_singlel(do_restart)

  ! synchronize processes
  call synchronize_all()

  ! checks if we want to restart simulation
  if (do_restart) then
    ! re-run simulation
    ! reset wavefields
    call prepare_wavefields()

    ! reset seismograms
    seismograms_d(:,:,:) = 0._CUSTOM_REAL
    seismograms_v(:,:,:) = 0._CUSTOM_REAL
    seismograms_a(:,:,:) = 0._CUSTOM_REAL
    seismograms_p(:,:,:) = 0._CUSTOM_REAL

    ! clear memory variables if attenuation
    if (ATTENUATION) then
      epsilondev_trace(:,:,:,:) = 0._CUSTOM_REAL
      epsilondev_xx(:,:,:,:) = 0._CUSTOM_REAL
      epsilondev_yy(:,:,:,:) = 0._CUSTOM_REAL
      epsilondev_xy(:,:,:,:) = 0._CUSTOM_REAL
      epsilondev_xz(:,:,:,:) = 0._CUSTOM_REAL
      epsilondev_yz(:,:,:,:) = 0._CUSTOM_REAL
      R_trace(:,:,:,:,:) = 0._CUSTOM_REAL
      R_xx(:,:,:,:,:) = 0._CUSTOM_REAL
      R_yy(:,:,:,:,:) = 0._CUSTOM_REAL
      R_xy(:,:,:,:,:) = 0._CUSTOM_REAL
      R_xz(:,:,:,:,:) = 0._CUSTOM_REAL
      R_yz(:,:,:,:,:) = 0._CUSTOM_REAL
      if (FIX_UNDERFLOW_PROBLEM) then
        R_trace(:,:,:,:,:) = VERYSMALLVAL
        R_xx(:,:,:,:,:) = VERYSMALLVAL
        R_yy(:,:,:,:,:) = VERYSMALLVAL
        R_xy(:,:,:,:,:) = VERYSMALLVAL
        R_xz(:,:,:,:,:) = VERYSMALLVAL
        R_yz(:,:,:,:,:) = VERYSMALLVAL
      endif
    endif

    if (USE_LDDRK) stop 'VTK restarting is not supported yet for LDDRK'

    ! puts elastic initial fields onto GPU
    if (GPU_MODE) then
      if (ELASTIC_SIMULATION) then
        ! transfer forward and backward fields to device with initial values
        call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
      endif
      if (ATTENUATION) stop 'VTK restarting is not supported yet for attenuation on GPUs'
    endif

    ! return to restart simulation
    return
  endif

  ! frees memory
  deallocate(vtkdata,vtkmask)

  if (NPROC > 1) then
    deallocate(vtkdata_points_all,vtkdata_offset_all)
    if (myrank == 0) deallocate(vtkdata_all)
  endif

  end subroutine vtk_window_cleanup

#endif

