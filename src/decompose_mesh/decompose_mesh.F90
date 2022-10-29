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


  subroutine decompose_mesh()

! divides model into partitions using scotch library functions

  use constants, only: SCOTCH,METIS,PATOH,ROWS_PART,MAX_STRING_LEN

  use decompose_mesh_par

  implicit none
  ! local parameters
  integer, dimension(:), allocatable  :: num_material
  integer :: ispec,iproc,inum,ier
  integer :: max_neighbor   ! real maximum number of neighbors per element

  ! debug
  integer,dimension(:,:),allocatable :: tmp_elmnts
  double precision,dimension(:),allocatable :: xstore_dummy,ystore_dummy,zstore_dummy
  ! vtk output
  character(len=MAX_STRING_LEN) :: filename

  !-------------------------------------------------------------------------------------------

  ! VTK-file output for debugging
  logical,parameter :: DEBUG_VTK_OUTPUT = .false.

  !-------------------------------------------------------------------------------------------

  ! safety check
  if (LTS_MODE) then
    if (num_p_level < 1) stop 'Error decompose_mesh: LTS has no p-levels'
  else
    if (num_p_level /= 1) stop 'Error decompose_mesh: no-LTS mode has invalid p-levels'
  endif

  ! assumes indexing starts from 0 for partitioners
  elmnts(:,:) = elmnts(:,:) - 1

  ! determines maximum neighbors based on "ncommonnodes" common nodes
  allocate(xadj(1:nspec+1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 108')
  if (ier /= 0) stop 'Error allocating array xadj'
  allocate(adjncy(1:sup_neighbor*nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 109')
  if (ier /= 0) stop 'Error allocating array adjncy'
  allocate(nnodes_elmnts(1:nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 110')
  if (ier /= 0) stop 'Error allocating array nnodes_elmnts'
  allocate(nodes_elmnts(1:nsize*nnodes),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 111')
  if (ier /= 0) stop 'Error allocating array nodes_elmnts'
  xadj(:) = 0; adjncy(:) = 0; nnodes_elmnts(:) = 0; nodes_elmnts(:) = 0

  ! user output
  print *,'partitioning:'
  print *,'  number of partitions requested = ',nparts
  print *
  print *,'  array size xadj  : ', size(xadj)," memory: ",size(xadj) * 4 / 1024./1024.,"MB"
  print *,'  array size adjncy: ', size(adjncy)," memory: ",size(adjncy) * 4 / 1024./1024.,"MB"
  print *,'  sup_neighbor     : ',sup_neighbor
  print *

  ncommonnodes = 1

  if (PARTITIONING_TYPE /= PATOH .or. nparts == 1) then
    call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbor, elmnts, xadj, adjncy, nnodes_elmnts, &
                                nodes_elmnts, max_neighbor, ncommonnodes, NGNOD)

    ! user output
    print *,'  mesh2dual: max_neighbor = ',max_neighbor
    print *

!! DK DK Oct 2012: added this safety test
    if (max_neighbor > sup_neighbor) stop 'found max_neighbor > sup_neighbor in domain decomposition'
  endif

  nb_edges = xadj(nspec+1)

  ! allocates & initializes partioning of elements
  allocate(part(1:nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 112')
  if (ier /= 0) stop 'Error allocating array part'
  part(:) = -1

  ! initializes
  ! elements load array
  allocate(elmnts_load(1:nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 113')
  if (ier /= 0) stop 'Error allocating array elmnts_load'
!! DK DK Oct 2012: this should include CPML weights as well in the future
  ! uniform load by default
  elmnts_load(:) = ACOUSTIC_LOAD

  ! gets materials id associations
  allocate(num_material(1:nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 114')
  if (ier /= 0) stop 'Error allocating array num_material'
  ! note: num_material can be negative for tomographic material elements
  !       (which are counted then as elastic elements)
  num_material(:) = mat(1,:)

  ! in case of acoustic/elastic/poro simulation, assign different weights to elements accordingly
  call acoustic_elastic_poro_load(elmnts_load,nspec,count_def_mat,count_undef_mat, &
                                  num_material,mat_prop,undef_mat_prop,ATTENUATION)

  ! checks if anything to do
  if (nparts <= 1) then
    ! no scotch/metis/.. partitioning used for serial compilation without MPI support.
    ! puts all elements in a single partition
    part(:) = 0

  else
    ! graph partitioning (SCOTCH / METIS / PATOH / ..)
    ! choice depending on setting in constants.h
    ! todo: this could be moved to the Par_file if desired for more flexibility

    ! user output
    print *
    print *,'partitioning type: ', PARTITIONING_TYPE

    select case (PARTITIONING_TYPE)
    case (SCOTCH)
      print *,'  running SCOTCH partitioning'
      print *
      call partition_scotch()

    case (METIS)
      print *,'  running METIS partitioning'
      print *
      call partition_metis()

    case (PATOH)
      print *,'  running PATOH partitioning'
      print *
      call partition_patoh()

    case (ROWS_PART)
      print *,'  running rows partitioning'
      print *
      call partition_rows()

    case default
      print *,'Invalid PARTITIONING_TYPE ',PARTITIONING_TYPE,' choosen, please select 1 / 2 / 3 or 4'
      stop 'Error invalid PARTITIONING_TYPE'
    end select
    print *,'finished partitioning'
    print *
  endif

  ! debug: vtk output
  if (DEBUG_VTK_OUTPUT) then
    ! element partitioning
    filename = trim(LOCAL_PATH)//'/part_array_before'
    allocate(tmp_elmnts(NGNOD,nspec))
    tmp_elmnts(:,:) = elmnts(:,:)+1
    allocate(xstore_dummy(nnodes),ystore_dummy(nnodes),zstore_dummy(nnodes))
    xstore_dummy(:) = nodes_coords(1,:)
    ystore_dummy(:) = nodes_coords(2,:)
    zstore_dummy(:) = nodes_coords(3,:)

    call write_VTK_data_ngnod_elem_i(NSPEC,nnodes,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
                                     tmp_elmnts,part,filename)

    print *,'  written file: ',trim(filename)//'.vtk'
    print *
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
    deallocate(tmp_elmnts)
  endif

  ! vtk output for lts
  if (LTS_MODE .and. SAVE_MESH_FILES) then
    call lts_save_p_level_partitions()
  endif

  ! check
  if (minval(part(:)) < 0 .or. maxval(part(:)) > nparts) then
    print *,'Error: partitioning has invalid slice number! partitioning min/max ',minval(part(:)),maxval(part(:))
    stop 'Error partitioning has invalid slice number'
  endif


  ! note: at this point, the partitioning is defined by array part(1:nspec) which gives for each
  !       element the partition number it is assigned to.
  !
  !       the following routines modify the partitioning array part() by putting elements at a
  !       common interface into the same partition.
  !
  !       this likely break the load balancing created by the domain decomposer for high-performance computing.
  !
  !       todo in future: try to avoid the need of having coupled elements in the same partition? ...

  ! re-partitioning puts poroelastic-elastic coupled elements into same partition
  if (PORO_INTERFACE_REPARTITIONING) &
    call poro_elastic_repartitioning(nspec, nnodes, elmnts, &
                                     count_def_mat, num_material , mat_prop, &
                                     sup_neighbor, nsize, &
                                     nparts, part, NGNOD)

  ! re-partitioning transfers two coupled elements on fault side 1 and side 2 to the same partition
  if (ANY_FAULT) &
    call fault_repartition(nspec, nnodes, elmnts, nsize, nparts, part, NGNOD, nodes_coords)

  ! re-partitioning puts moho-surface coupled elements into same partition
  if (SAVE_MOHO_MESH) &
    call moho_surface_repartitioning(nspec, nnodes, elmnts, &
                                     sup_neighbor, nsize, nparts, part, &
                                     nspec2D_moho,ibelm_moho,nodes_ibelm_moho, NGNOD, NGNOD2D)

  ! final partitioning: vtk output
  if (SAVE_MESH_FILES) then
    ! element partitioning
    filename = trim(LOCAL_PATH)//'/part_array'
    allocate(tmp_elmnts(NGNOD,nspec))
    tmp_elmnts(:,:) = elmnts(:,:)+1
    allocate(xstore_dummy(nnodes),ystore_dummy(nnodes),zstore_dummy(nnodes))
    xstore_dummy(:) = nodes_coords(1,:)
    ystore_dummy(:) = nodes_coords(2,:)
    zstore_dummy(:) = nodes_coords(3,:)

    call write_VTK_data_ngnod_elem_i(NSPEC,nnodes,NGNOD,xstore_dummy,ystore_dummy,zstore_dummy, &
                                     tmp_elmnts,part,filename)

    print *,'  written file: ',trim(filename)//'.vtk'
    print *
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
    deallocate(tmp_elmnts)
  endif

  ! number of elements per partition
  call print_statistics_distribution()

  ! load distribution
  call print_statistics_load()

  ! local number of each element for each partition
  call build_glob2loc_elmnts(nspec, part, glob2loc_elmnts,nparts)

  ! local number of each node for each partition
  call build_glob2loc_nodes(nspec, nnodes,nsize, nnodes_elmnts, nodes_elmnts, part, &
                            glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, nparts)

  ! MPI interfaces
  ! acoustic/elastic/poroelastic boundaries will be split into different MPI partitions
  call build_interfaces(nspec, sup_neighbor, part, elmnts, &
                        xadj, adjncy, tab_interfaces, &
                        tab_size_interfaces, ninterfaces, &
                        nparts, NGNOD)

  ! obsolete: from when we wanted acoustic/elastic boundaries NOT to be separated into different MPI partitions
  ! call build_interfaces_no_ac_el_sep(nspec, sup_neighbor, part, elmnts, &
  !                           xadj, adjncy, tab_interfaces, &
  !                           tab_size_interfaces, ninterfaces, &
  !                           count_def_mat, mat_prop(3,:), mat(1,:), nparts, NGNOD)

  ! communication (edgecut) statistics
  if (LTS_MODE) then
    call lts_communication_load(ispec_p_refine,nspec)
  endif

  ! frees memory
  deallocate(elmnts_load)
  deallocate(num_material)

  deallocate(xadj,adjncy)
  deallocate(nnodes_elmnts,nodes_elmnts)

contains

    subroutine print_statistics_distribution()

    implicit none
    ! local parameters
    integer,dimension(:),allocatable :: num_pelem_part
    integer,dimension(:),allocatable :: min_pelem,max_pelem
    integer :: elem_min,elem_max
    real :: elem_balance
    integer :: i,p

    print *, 'element distribution:'
    elem_min = nspec
    elem_max = 0
    if (LTS_MODE) then
      ! local time stepping
      allocate(num_pelem_part(num_p_level), &
               min_pelem(num_p_level),max_pelem(num_p_level),stat=ier)
      if (ier /= 0) stop 'error allocating num_pelem_part'
      min_pelem(:) = NSPEC
      max_pelem(:) = 0
      do iproc = 0,nparts-1
        ! number of element in this process
        inum = count(part(:) == iproc )
        if (inum < elem_min) elem_min = inum
        if (inum > elem_max) elem_max = inum
        num_pelem_part(:) = 0
        do ispec = 1,nspec
          if (part(ispec) == iproc) then
            do i = 1,num_p_level
              p = p_level(i)
              if (ispec_p_refine(ispec) == p) num_pelem_part(i) = num_pelem_part(i) + 1
            enddo
          endif
        enddo
        ! min/max per p-level
        do i = 1,num_p_level
          if (num_pelem_part(i) < min_pelem(i)) min_pelem(i) = num_pelem_part(i)
          if (num_pelem_part(i) > max_pelem(i)) max_pelem(i) = num_pelem_part(i)
        enddo
        ! user output
        print *,'  partition ',iproc, '       has ',inum,' elements with p-elements = ',num_pelem_part(:)
      enddo
      print *
      ! partition balance
      do i = 1,num_p_level
        elem_min = min_pelem(i)
        elem_max = max_pelem(i)
        ! imbalance in percent
        elem_balance = 100.0 * ( elem_max - elem_min ) / elem_max
        ! user output
        print *,'  p-level: ',i,' p = ',p_level(i)
        print *,'    elements per partition: min/max   = ',elem_min,elem_max
        print *,'    elements per partition: imbalance = ',elem_balance,'%'
        print *,'                            (0% being totally balanced, 100% being unbalanced)'
        print *
      enddo
      print *
      ! frees temporary arrays
      deallocate(min_pelem,max_pelem)
      deallocate(num_pelem_part)
    else
      ! non-LTS
      do iproc = 0,nparts-1
        ! number of element in this process
        inum = count(part(:) == iproc )
        if (inum < elem_min) elem_min = inum
        if (inum > elem_max) elem_max = inum
        ! user output
        print *,'  partition ',iproc, '       has ',inum,' elements'
      enddo
      ! user output
      if (elem_max > 0) then
        ! imbalance in percent
        elem_balance = 100.0 * ( elem_max - elem_min ) / elem_max
      else
        stop 'error element count: partition has zero elements'
      endif
      print *,'  elements per partition: min/max   = ',elem_min,elem_max
      print *,'  elements per partition: imbalance = ',elem_balance,'%'
      print *,'                          (0% being totally balanced, 100% being unbalanced)'
      print *
    endif

    end subroutine print_statistics_distribution

    !--------------------------------------------------------------------

    subroutine print_statistics_load()

    implicit none
    ! local parameters
    integer :: load_min,load_max
    real :: load_balance

    print *, 'load distribution:'
    print *, '  element loads: min/max = ',minval(elmnts_load(:)),maxval(elmnts_load(:))
    print *

    load_min = nspec * maxval(elmnts_load(:))
    ! check integer size limit
    if (nspec > int(2147483646.0 / maxval(elmnts_load(:)))) then
      load_min = 2147483646
    endif
    load_max = 0
    do iproc = 0,nparts-1
      ! counts loads
      inum = 0
      do ispec = 1,nspec
        if (part(ispec) == iproc) then
          ! load per element
          inum = inum + elmnts_load(ispec)
        endif
      enddo
      if (inum < load_min) load_min = inum
      if (inum > load_max) load_max = inum
      ! user output
      print *,'  partition ',iproc, '       has ',inum,' load units'
    enddo
    if (load_max > 0) then
      ! imbalance in percent
      load_balance = 100.0 * ( load_max - load_min ) / load_max
    else
      stop 'error element count: partition has zero load'
    endif
    print *,'  load per partition: min/max   = ',load_min,load_max
    print *,'  load per partition: imbalance = ',load_balance,'%'
    print *,'                      (0% being totally balanced, 100% being unbalanced)'
    print *

    end subroutine print_statistics_load

  end subroutine decompose_mesh

