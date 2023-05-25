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


! Local Time Stepping

  module lts_generate_databases_par

  use constants, only: LTS_OVERLAP_REGION,LTS_SAFETY_MARGIN,LTS_SINGLE_P_LEVEL,LTS_TWO_P_LEVEL, &
                       LTS_DECREASE_DT,LTS_STABILITY_MARGIN_DT

  implicit none

  ! time steps
  double precision :: deltat_lts

  ! element p-refinement values (like 1 2 4 8 ..; p == 1 being coarsest, p == 8 finer local time step dt/p )
  integer,dimension(:),allocatable :: ispec_p_refine

  integer :: num_p_level
  integer,dimension(:),allocatable :: p_level
  integer, dimension(:),allocatable :: p_level_loops

  ! p-level stepping
  integer :: num_p_level_steps
  integer,dimension(:),allocatable :: p_level_steps

  integer,dimension(:),allocatable :: iglob_p_refine
  integer,dimension(:),allocatable :: p_level_iglob_start, p_level_iglob_end

  ! element flags
  logical,dimension(:,:),allocatable :: p_elem
  logical,dimension(:,:),allocatable :: boundary_elem

  end module lts_generate_databases_par

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_generate_databases()

! puts elements into separate refinement levels for local time stepping
! (see also routine in decompose_mesh.F90)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IMAIN,COURANT_SUGGESTED,myrank,itag

  use shared_parameters, only: NPROC,DT

  use generate_databases_par, only: ibool,NSPEC_AB,NGLOB_AB

  ! MPI interfaces
  use generate_databases_par, only: num_interfaces_ext_mesh,my_neighbors_ext_mesh, &
    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh

  use create_regions_mesh_ext_par, only: &
    xstore => xstore_unique, &
    ystore => ystore_unique, &
    zstore => zstore_unique, &
    nglob_unique, &
    kappastore,mustore,rhostore

  use fault_generate_databases, only: ANY_FAULT,ANY_FAULT_IN_THIS_PROC

  ! LTS module
  use lts_generate_databases_par

  implicit none

  ! local parameters
  ! ratio of minimum distance of GLL points, depending on NGLL
  double precision,dimension(15) :: percent_GLL

  ! time steps
  double precision :: dtmin,dtmax
  double precision :: dtmin_glob,dtmax_glob
  double precision :: dt_suggested

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: ispec_time_step

  ! p refinements ( local time step = dt/p )
  integer :: p
  integer :: p_min,p_max
  integer :: p_min_glob,p_max_glob
  integer :: p_fine,p_coarse

  integer :: max_p_l, p_elem_counter

  ! local arrays
  integer,dimension(:),allocatable :: tmp_p_level
  integer,dimension(:),allocatable :: tmp_i

  integer, dimension(:),allocatable :: p_level_relative

  ! p-level grid
  integer,dimension(:,:),allocatable :: p_level_grid
  integer,dimension(:),allocatable :: p_level_grid_Size

  integer :: i,j,k,ispec,iglob,ier
  integer :: iproc,inum,ii,jj
  integer :: step,l,ilevel
  logical :: found

  integer,dimension(:),allocatable :: num_ispec_level
  double precision :: total_work,total_work_lts,total_work_lts_b

  real(kind=CUSTOM_REAL) :: vpmin,vpmax,vsmin,vsmax
  real(kind=CUSTOM_REAL) :: poissonmin,poissonmax
  real(kind=CUSTOM_REAL) :: distance_min!,distance_max
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max

  logical :: has_vs_zero

  ! MPI communication
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh

  ! fault p-value re-assignments
  integer,dimension(:),allocatable :: tmp_iglob_p_refine
  integer :: iloop
  integer :: ichanged,ichanged_glob


  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     setting up elements for local time stepping'
    write(IMAIN,*) '     number of elements:',NSPEC_AB
    write(IMAIN,*) '     number of nodes   :',NGLOB_AB
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! safety check
  if (nglob_unique /= NGLOB_AB) stop 'Error nglob_unique not equal to NGLOB_AB'

  ! define percentage of smallest distance between GLL points for NGLL points
  ! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(1) = 0.d0
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

  ! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0
  if (NGLLX < 2 .or. NGLLX > 15) stop 'Error value of NGLLX not supported yet'

  ! allocates memory for p-refinement per element
  allocate(ispec_p_refine(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array ispec_p_refine'
  ispec_p_refine(:) = 0

  ! temporary arrays
  allocate(ispec_time_step(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array ispec_time_step'
  ispec_time_step(:) = 0.0

  allocate(iglob_p_refine(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array iglob_p_refine'
  iglob_p_refine(:) = 0

  ! estimated global, stable time step based on CFL condition
  dtmin = 1.e30
  dtmax = 0.0
  do ispec = 1, NSPEC_AB
    ! determines minimum/maximum velocities within this element
    call get_vpvs_minmax(vpmin,vpmax,vsmin,vsmax,poissonmin,poissonmax, &
                         ispec,has_vs_zero, &
                         NSPEC_AB,kappastore,mustore,rhostore)

    ! computes minimum and maximum size of this grid cell
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec, &
                             NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    ! minimum distance based on GLL point distribution
    distance_min = elemsize_min * percent_GLL(NGLLX)

    ! alternative:
    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
    !call get_GLL_minmaxdistance(distance_min,distance_max,ispec, &
    !                            NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore)

    ! estimated time step based on CFL condition
    dt_suggested = COURANT_SUGGESTED * distance_min / max( vpmax,vsmax )

    ! cut at a significant number of digits (2 digits)
    ! example: 0.0734815 -> lpow = (2 - (-1) = 3 -> 0.0730
    call get_timestep_limit_significant_digit_dp(dt_suggested)

    ! stores value
    ispec_time_step(ispec) = dt_suggested

    ! determines coarsest time step
    if (dtmax < dt_suggested) dtmax = dt_suggested

    ! determines finest time step
    if (dtmin > dt_suggested) dtmin = dt_suggested
  enddo

  ! gets min/max values over all processes
  call min_all_all_dp(dtmin,dtmin_glob)
  call max_all_all_dp(dtmax,dtmax_glob)
  dtmin = dtmin_glob
  dtmax = dtmax_glob

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     estimated time step min   = ',dtmin,' seconds'
    write(IMAIN,*) '     estimated time step max   = ',dtmax,' seconds'
    write(IMAIN,*) '     estimated time step ratio = ',dtmax/dtmin
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! users sets DT = .. in Par_file, we use this as a threshold limit for the minimum time step size
  ! (useful in case CFL is underestimating size)
  if (dtmin > DT) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '     DT time step size set by Par_file: DT = ',DT,' limits smallest time step size'
      if (dtmin / DT >= 2.0) then
        write(IMAIN,*) '     Please set DT to higher value    :      ',dtmin,' otherwise LTS will not be optimal'
      endif
      write(IMAIN,*)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     suggested minimum DT time step = ',dtmin
    call flush_IMAIN()
  endif

  ! we will use now the smallest time step estimate and put elements
  ! into bins with increasing time step size (power of 2 of smallest time step)

  ! setup p' refinements for each global point: local time step == dtmin * p'
  iglob_p_refine(:) = 1000000000
  do ispec = 1,NSPEC_AB

    ! gets estimated time step based on CFL condition
    dt_suggested = ispec_time_step(ispec)

    ! uses a user-defined safety margin for binning
    dt_suggested = ( 1.d0 - LTS_SAFETY_MARGIN ) * dt_suggested

    ! absorbing conditions need slightly smaller time steps
    ! if (ABSORBING_CONDITIONS) then
    !   ! uses a 60% safety margin for binning
    !   dt_suggested = ( 1.d0 - 0.3d0 ) * dt_suggested
    ! endif

    ! local time step level compared to global time step size
    p = floor( dt_suggested / dtmin )
    if (p < 1) p = 1
    ! debug
    !print *,'ispec',ispec,'time step=',dt_suggested,dtmin,'     p refinement = ',p

    ! debug
    ! uniform distribution for testing
    if (LTS_SINGLE_P_LEVEL) p = 1

    ! two-level distribution for testing
    if (LTS_TWO_P_LEVEL) then
      if (p > 1) p = 2
    endif

    ! enforces p' to be power of 2, find next smallest power of 2:
    ! example: 2 -> 2, 3 -> 2, 4 -> 4, 5 -> 4, .., 7 -> 4, 8 -> 8, 9 -> 8, ..
    if (p > 1) then
      i = 0
      do while ( 2**i <= p )
        i = i + 1
      enddo
      p = 2**(i-1)
    endif

    ! sets p' refinements on global points
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! sets only if smaller step needed
          if (p < iglob_p_refine(iglob)) iglob_p_refine(iglob) = p
        enddo
      enddo
    enddo
  enddo

  ! gets min/max over all processes
  p_min = minval(iglob_p_refine(:))
  p_max = maxval(iglob_p_refine(:))
  call min_all_all_i(p_min,p_min_glob)
  call max_all_all_i(p_max,p_max_glob)
  p_min = p_min_glob
  p_max = p_max_glob

  ! check
  if (p_min >= 1000000000) then
    print *,'Error: p minimum ',p_min,'on rank',myrank,'global',p_min_glob
    stop 'Error minimum p refinement must be 1, please check element distribution...'
  endif
  if (p_max >= 1000000000) then
    print *,'Error: p maximum ',p_max,'on rank',myrank,'global',p_max_glob
    stop 'Error maximum p refinement must be at least 1, please check element distribution...'
  endif
  call synchronize_all()

  ! coarsest time step for LTS (largest multiple p of smallest time step)
  deltat_lts = dtmin * dble(p_max)

  ! takes away percent of optimal time step for LTS stability
  ! LTS must decrease time step for stability depending on p depth
  ! to avoid this decrease, one can try to overlap the fine region by one/two elements
  if (LTS_DECREASE_DT) then
    ! sets new stable coarse time step for LTS
    deltat_lts = ( 1.d0 - LTS_STABILITY_MARGIN_DT ) * deltat_lts
    ! user output
    if (myrank == 0) write(IMAIN,*) '     decreases cfl time step by ',LTS_STABILITY_MARGIN_DT*100.0,'percent'
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     suggested global coarsest time step       = ',deltat_lts
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! re-orders p' -> p = 1 / p'
  ! note: 1 == coarsest, p_max == finest level such that local time step == dt / p
  do iglob = 1,NGLOB_AB
    p = iglob_p_refine(iglob)
    iglob_p_refine(iglob) = p_max / p
  enddo

  ! from here on, a higher p-value indicates a smaller local time step dt/p

  ! synchronizes global points shared with other processes
  if (ANY_FAULT) then
    ! re-assigns p-values on split nodes to have same p-value
    ! note: tohoku mesh has elements in 1 partition which only touches split nodes on edges.
    !       re-assigning p-values here makes sure that also after MPI synchronization,
    !       p-values stay the same on split nodes

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '     fault re-assigning p-values: '
      call flush_IMAIN()
    endif

    ! to store initial p-values
    allocate(tmp_iglob_p_refine(NGLOB_AB),stat=ier)
    if (ier /= 0) stop 'Error allocating array tmp_iglob_p_refine'
    tmp_iglob_p_refine(:) = 0

    ! we do this a few times to allow for changes
    iloop = 0
    ichanged = 1
    do while ( ichanged == 1 )

      ! counter
      iloop = iloop + 1

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '       re-assignment loop ',iloop
        call flush_IMAIN()
      endif

      ! stores initial p-values
      tmp_iglob_p_refine(:) = iglob_p_refine(:)

      ! makes sure split nodes have same p-value
      if (ANY_FAULT_IN_THIS_PROC) call lts_fault_assign_split_p_values()

      ! synchronizes global points shared with other processes
      if (NPROC > 1) then
        ! sets up MPI communications
        max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
        allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
        if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh_dummy'
        do i = 1, num_interfaces_ext_mesh
           ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
        enddo

        ! gets maximum p values on interfaces with other partitions
        call assemble_MPI_scalar_i_max(NPROC,NGLOB_AB,iglob_p_refine, &
                                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy, &
                                       my_neighbors_ext_mesh)

        ! frees temporary array
        deallocate(ibool_interfaces_ext_mesh_dummy)
      endif

      ! checks if any p-value changed
      ichanged = 0
      do iglob = 1,NGLOB_AB
        if (tmp_iglob_p_refine(iglob) /= iglob_p_refine(iglob)) then
          ichanged = 1
          exit
        endif
      enddo

      ! synchronize value with all others
      call max_all_all_i(ichanged,ichanged_glob)
      ichanged = ichanged_glob

      ! checks iloop counter
      if (iloop >= 10) then
        if (myrank == 0) then
          print *,'fault reassigning p-values exceeds ',iloop,'loops '
          print *,'exiting now, please check your mesh and setup...'
        endif
        call exit_mpi(myrank,'Error iloop counter exceeds tolerance in fault p-value re-assignment')
      endif
    enddo

    ! frees temporary array
    deallocate(tmp_iglob_p_refine)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  else
    ! only need to do this once
    if (NPROC > 1) then
      ! sets up MPI communications
      max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
      allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
      if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh_dummy'
      do i = 1, num_interfaces_ext_mesh
         ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
      enddo

      ! gets maximum p values on interfaces with other partitions
      call assemble_MPI_scalar_i_max(NPROC,NGLOB_AB,iglob_p_refine, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy, &
                                     my_neighbors_ext_mesh)

      ! frees temporary array
      deallocate(ibool_interfaces_ext_mesh_dummy)
    endif
  endif

  ! debug
  !do ispec = 1,NSPEC_AB
  !  do k = 1,NGLLZ
  !    do j = 1,NGLLY
  !      do i = 1,NGLLX
  !        iglob = ibool(i,j,k,ispec)
  !        print *,ispec,i,j,k,iglob,'iglob p value = ',iglob_p_refine(iglob)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! re-gets min/max p value on global nodes over all processes
  p_min = minval(iglob_p_refine(:))
  p_max = maxval(iglob_p_refine(:))
  call min_all_all_i(p_min,p_min_glob)
  call max_all_all_i(p_max,p_max_glob)
  p_min = p_min_glob
  p_max = p_max_glob

  ! sets p refinement for element
  ! (adds also elements touching global points with higher refinement)
  ispec_p_refine(:) = 0
  do ispec = 1,NSPEC_AB
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          p = iglob_p_refine(iglob)
          if (p > ispec_p_refine(ispec)) ispec_p_refine(ispec) = p
        enddo
      enddo
    enddo
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     p refinement of nodes   : min/max = ',p_min,p_max
    write(IMAIN,*) '     p refinement of elements: min/max = ',minval(ispec_p_refine(:)),maxval(ispec_p_refine(:))
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! adds coarser elements which overlap to finer p-level
  if (LTS_OVERLAP_REGION) call lts_add_overlap_elements_databases(iglob_p_refine,ispec_p_refine)

  ! temporary array to store different p values
  allocate(tmp_i(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating array tmp_i'
  tmp_i(:) = 0

  ! counts number of needed different levels of refinement
  num_p_level = 0
  do ispec = 1,NSPEC_AB
    p = ispec_p_refine(ispec)
    ! search in table
    found = .false.
    do i = 1,num_p_level
      if (p == tmp_i(i)) then
        found = .true.
        exit
      endif
    enddo
    ! inserts new p level
    if (.not. found) then
      ! increases levels
      num_p_level = num_p_level + 1
      ! adds new entry into sorted list
      tmp_i(num_p_level) = p
      ! sorts list in decreasing order (like in: 8 4 3 1)
      do i = 2, num_p_level
        j = i - 1
        p = tmp_i(i)
        do while (j > 1 .and. tmp_i(j) < p)
          tmp_i(j+1) = tmp_i(j)
          j = j - 1
        enddo
        if (j == 1 .and. tmp_i(j) < p) then
          tmp_i(j+1) = tmp_i(j)
          j = j - 1
        endif
        tmp_i(j+1) = p
      enddo
    endif
  enddo
  ! check
  if (num_p_level < 1) stop 'Error no p-level found'

  ! debug
  !print *,'process',myrank,'num_p_level = ',num_p_level
  !call synchronize_all()

  ! number of local time step refinements for each LTS domain
  ! note: may vary for different partitions/processes
  allocate(p_level(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level'
  p_level(:) = tmp_i(1:num_p_level)

  ! collect and assemble full array on master
  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1
        ! gets length of p_level array in other process
        call recv_singlei(inum,iproc,itag)

        ! gets p_level array from other process
        tmp_i(:) = 0
        call recv_i(tmp_i(1:inum),inum,iproc,itag)

        ! checks if we have all p values
        do i = 1,inum
          p = tmp_i(i)
          ! search in master table
          found = .false.
          do j = 1,num_p_level
            if (p == p_level(j)) then
              found = .true.
              exit
            endif
          enddo
          ! inserts new p level
          if (.not. found) then
            ! adds p value to list and increases number of levels
            ! make backup of list
            allocate(tmp_p_level(num_p_level),stat=ier)
            if (ier /= 0) stop 'Error allocating array tmp_p_level'
            tmp_p_level(:) = p_level(:)
            ! re-allocate new list
            num_p_level = num_p_level + 1
            deallocate(p_level)
            allocate(p_level(num_p_level),stat=ier)
            if (ier /= 0) stop 'Error allocating new array p_level'
            p_level(1:num_p_level-1) = tmp_p_level(:)
            deallocate(tmp_p_level)
            ! adds new entry into sorted list
            p_level(num_p_level) = p
            ! sorts list in decreasing order (like in: 8 4 3 1)
            do ii = 2, num_p_level
              jj = ii - 1
              p = p_level(ii)
              do while (jj > 1 .and. p_level(jj) < p)
                p_level(jj+1) = p_level(jj)
                jj = jj - 1
              enddo
              if (jj == 1 .and. p_level(jj) < p) then
                p_level(jj+1) = p_level(jj)
                jj = jj - 1
              endif
              p_level(jj+1) = p
            enddo
          endif
        enddo ! inum
      enddo ! NPROC-1
    else
      call send_singlei(num_p_level,0,itag)
      call send_i(p_level(1:num_p_level),num_p_level,0,itag)
    endif

    ! broadcast new p_level to all processes
    call bcast_all_singlei(num_p_level)

    ! uses temporary array to broadcast
    allocate(tmp_p_level(num_p_level),stat=ier)
    if (ier /= 0) stop 'Error allocating array tmp_p_level'

    if (myrank == 0) then
      tmp_p_level(:) = p_level(:)
    endif
    call bcast_all_i(tmp_p_level,num_p_level)

    ! re-set p_level array on all processes
    deallocate(p_level)
    allocate(p_level(num_p_level),stat=ier)
    if (ier /= 0) stop 'Error allocating new array p_level'
    p_level(:) = tmp_p_level(:)

    ! frees temporary array
    deallocate(tmp_p_level)
  endif ! NPROC
  ! frees temporary array
  deallocate(tmp_i)

  ! checks array
  if (minval(p_level(:)) > 1) stop 'Error coarsest p level starts with p > 1'
  if (maxval(p_level(:)) /= p_max) stop 'Error smallest p level has wrong maximum value'
  if (.not. p_max == p_level(1)) stop 'Error first p_level does not start with maximum p value'

  ! sets up relative depth (in power of 2 ) of fine region to enclosing coarse one
  ! like from 16 8 2 1 -> 2 2 1 1 which works with the p_level_grid stepping
  allocate(p_level_relative(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_relative'

  p_level_relative(:) = 1
  do i = 1,num_p_level-1
    p_fine = p_level(i)
    p_coarse = p_level(i+1)
    ! relative depth of fine region to coarse one
    ! (how many steps to reach coarser one)
    if (mod(p_fine,2*p_coarse) == 0) then
      p = p_fine / (2*p_coarse)
    else
      do while (mod(p_fine,2*p_coarse) /= 0 )
        p_fine = p_fine + 1
        p = p_fine / ( 2*p_coarse)
      enddo
    endif
    p_level_relative(i) = p
    if (myrank == 0) write(IMAIN,*) '     level',i,':  fine/coarse p refinement = ',p_fine,p_coarse,'relative = ',p
  enddo
  ! finest level must increase time steps by factor 2 with current p_level_grid stepping scheme
  if (num_p_level > 1) then
    p_level_relative(1) = 2 * p_level_relative(1)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '     number of p-levels    : ',num_p_level
    write(IMAIN,*) '     p-level array         : ',p_level(:)
    write(IMAIN,*) '     p-level relative array: ',p_level_relative(:)
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! setup of p-level grid stepping
  ! here you need to define max numbers of p-level steps according to the mesh

  ! calculates grid size
  allocate(p_level_grid_Size(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_grid_Size'

  if (num_p_level == 1) then
    p_level_grid_Size(1) = 1
  else if (num_p_level == 2) then
    p_level_grid_Size(1) = 1
    p_level_grid_Size(2) = 2
  else
    ! grid size
    p_level_grid_Size(1) = 1
    p_level_grid_Size(2) = 2
    do i = 3, num_p_level
      p_level_grid_Size(i) = 2 * p_level_grid_Size(i-1) * p_level_relative(i-1) + 1
    enddo
  endif

  !debug
  ! print grid size
  !print *,'p_level_grid_Size:'
  !do i=1, num_p_level
  !  print *,'debug: grid ', i, ' has size of ', p_level_grid_Size(i)
  !enddo
  !print *

  ! creates p_level_grid steps
  num_p_level_steps = p_level_grid_Size(num_p_level)

  allocate(p_level_grid(num_p_level_steps, num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_grid'

  p_level_grid(:,:) = 0
  if (num_p_level == 1) then
    p_level_grid(1,1) = 1
  else if (num_p_level == 2) then
    p_level_grid(1,1) = 1
    p_level_grid(1,2) = 1
    p_level_grid(2,2) = 2
  else
    ! Initialize first p_level_grid (no coarse grid)
    p_level_grid(1,1) = 1
    ! Initialize second p_level_grid (none fine grid, coarse grid)
    p_level_grid(1,2) = 1
    p_level_grid(2,2) = 2
    ! Initialize second grid (none fine grid, coarse grid)
    do i = 3, num_p_level
      l = 1
      do j = 1, 2*p_level_relative(i-1)
        do k = 1, p_level_grid_Size(i-1)
          p_level_grid(l,i) = p_level_grid(k,i-1)
          l = l + 1
        enddo
      enddo
      ! add last node the grid
      p_level_grid(l,i) = i
    enddo
  endif

  ! stores only relevant grid steps
  allocate(p_level_steps(num_p_level_steps),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_steps'
  allocate(p_level_loops(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_loops'

  p_level_steps(:) = p_level_grid(:,num_p_level)
  p_level_loops(:) = 0
  do l = 1,num_p_level
    ! counts number of level occurrences in p_level_steps array
    k = 0
    do step = 1, num_p_level_steps
      if (p_level_steps(step) == l) k = k + 1
    enddo
    ! checks occurrences
    if (k == 0) then
      print *,'Error occurrence: no level',l
      stop 'Error no occurrence in p_level_steps'
    endif
    ! checks dimensions
    if (k > p_level(l)) then
      print *,'Error occurrence: level',l,'occurs',k
      stop 'Error too many occurrences'
    endif
    p_level_loops(l) = p_level(l) / k
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     maximum p-level steps : ',num_p_level_steps
    write(IMAIN,*) '     p-level step array    : ',p_level_steps(:)
    write(IMAIN,*) '     p-level loops array   : ',p_level_loops(:)
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! frees temporary arrays
  deallocate(p_level_grid_Size)
  deallocate(p_level_grid)

  ! determines elements which belong to p-levels
  allocate(p_elem(NSPEC_AB,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array p_elem'
  allocate(boundary_elem(NSPEC_AB,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array boundary_elem'

  p_elem(:,:) = .false.
  boundary_elem(:,:) = .false.

  ! loops over refinement levels
  do ilevel = 1, num_p_level
    ! p refinement number in this specified p-level
    p = p_level(ilevel)

    ! sets element flags in this p-level
    ! p_elem is not duplicated over levels
    ! if one level is true, then all others are false
    do ispec = 1,NSPEC_AB
      max_p_l = 0
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            max_p_l = max(iglob_p_refine(iglob),max_p_l)
            ! if (iglob_p_refine(iglob) == p) p_elem(ispec,ilevel) = .true.
          enddo
        enddo
      enddo
      if (max_p_l == p) p_elem(ispec,ilevel) = .true.
    enddo
  enddo

  ! test p_elem
  do ispec = 1,NSPEC_AB
    ! test p_elem correctness
    p_elem_counter = 0
    do ilevel=1,num_p_level
      if (p_elem(ispec,ilevel) .eqv. .true.) p_elem_counter = p_elem_counter + 1
    enddo
    if (p_elem_counter /= 1) then
      print *, "ERROR: p_counter:",p_elem_counter
      call exit_mpi(myrank,"ASSERT(p_elem(ispec,:) only has 1 true) FAIL")
    endif
  enddo

  ! sets any boundary element flags with surrounding levels
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    do ispec = 1,NSPEC_AB
      ! set boundary elements
      ! boundaries are not doubled -- only one level owns the boundary
      if (p_elem(ispec,ilevel) .eqv. .true.) then
        ! if any node is in another p-level, this is a boundary
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              if (iglob_p_refine(iglob) < p) then
                boundary_elem(ispec,ilevel) = .true.
              else if (iglob_p_refine(iglob) > p) then
                call exit_mpi(myrank,"ASSERT(p-elem == .true. should if contains this p-level or coarser nodes)")
              endif
            enddo
          enddo
        enddo
      endif
    enddo
  enddo

  ! test boundary_elem
  do ispec = 1,NSPEC_AB
    ! test p_elem correctness
    p_elem_counter = 0
    do ilevel=1,num_p_level
      if (boundary_elem(ispec,ilevel) .eqv. .true.) p_elem_counter = p_elem_counter + 1
    enddo
    if (p_elem_counter > 1) then
      print *, "ERROR: p_counter:",p_elem_counter
      call exit_mpi(myrank,"ASSERT(boundary_elem(ispec,:) only has 1 true) FAIL")
    endif
  enddo

  ! counts number of elements for each level
  allocate(num_ispec_level(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating array num_ispec_level'
  num_ispec_level(:) = 0
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    num_ispec_level(ilevel) = count(ispec_p_refine(:) == p)
  enddo
  ! note: some partitions might not contain any element within a certain p-level
  ! collect on master
  if (NPROC > 1) then
    do ilevel = 1,num_p_level
      call sum_all_all_i(num_ispec_level(ilevel),inum)
      num_ispec_level(ilevel) = inum
    enddo
  endif
  if (minval(num_ispec_level(:)) == 0) stop 'Error p level without element found'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     p-level number of elements: ',num_ispec_level(:)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! theoretical speed-up values (valid only for all elements belonging to same domain type)
  total_work_lts = 0.0
  total_work = 0.0
  if (myrank == 0) then
    do ilevel = 1,num_p_level
      p = p_level(ilevel)
      total_work_lts = total_work_lts + num_ispec_level(ilevel) * p
    enddo
    ! work without lts
    total_work = sum( num_ispec_level(:)) * maxval(p_level(:))
  endif

  ! counts boundary elements
  num_ispec_level(:) = 0
  do ilevel = 1,num_p_level
    do ispec = 1,NSPEC_AB
      if (boundary_elem(ispec,ilevel)) then
        num_ispec_level(ilevel) = num_ispec_level(ilevel) + 1
      endif
    enddo
  enddo
  ! note: some partitions might not contain any element within a certain p-level
  ! collect on master
  if (NPROC > 1) then
    do ilevel = 1,num_p_level
      call sum_all_all_i(num_ispec_level(ilevel),inum)
      num_ispec_level(ilevel) = inum
    enddo
  endif
  ! work with additional coarse/fine boundary contributions
  total_work_lts_b = total_work_lts
  if (myrank == 0) then
    do ilevel = 1,num_p_level
      p = p_level(ilevel)
      total_work_lts_b = total_work_lts_b + num_ispec_level(ilevel) * p
    enddo
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     p-level boundary elements : ',num_ispec_level(:)
    write(IMAIN,*)
    write(IMAIN,*) '     theoretical speed-up value: ',sngl(total_work / total_work_lts),'(without boundary contributions)'
    write(IMAIN,*) '     theoretical speed-up value: ',sngl(total_work / total_work_lts_b),'(with boundary contributions)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! frees temporary array
  deallocate(num_ispec_level)

  ! re-order global nodes, puts nodes of same p-level together
  allocate(p_level_iglob_start(num_p_level), &
           p_level_iglob_end(num_p_level), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating array p_level_iglob_start,..'
  ! initializes start/end index arrays
  p_level_iglob_start(:) = 0
  p_level_iglob_end(:) = 0

  call lts_reorder_iglob_by_p_level()

  ! stores arrays in databases for solver
  call lts_save_databases()

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '     all done LTS'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! free memory
  deallocate(p_level)
  deallocate(p_level_relative)
  deallocate(p_level_loops)
  deallocate(p_level_steps)
  deallocate(ispec_p_refine)
  deallocate(iglob_p_refine)
  deallocate(p_elem)
  deallocate(boundary_elem)
  deallocate(p_level_iglob_start)
  deallocate(p_level_iglob_end)

  end subroutine lts_generate_databases

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_reorder_iglob_by_p_level()

! ## builds p_level_iglob_start,end ##
! re-orders ibool such that global nodes in a p-level are continuously listed

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,myrank

  use generate_databases_par, only: NSPEC_AB,NGLOB_AB, &
    ibool,num_interfaces_ext_mesh,ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh

  use create_regions_mesh_ext_par, only: &
    xstore => xstore_unique, &
    ystore => ystore_unique, &
    zstore => zstore_unique

  use fault_generate_databases, only: ANY_FAULT,ANY_FAULT_IN_THIS_PROC

  ! LTS module
  use lts_generate_databases_par, only: p_level_iglob_start,p_level_iglob_end,num_p_level,iglob_p_refine,p_level

  implicit none

  ! local parameters
  integer :: ispec, iglob, iglob_new, ip, p, ier, i,j,k, iinterface
  integer, dimension(:), allocatable :: iglob_touched
  integer, dimension(:,:,:,:), allocatable :: ibool_new
  integer, dimension(:), allocatable :: iglob_field_new

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: iglob_field_new_cr
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_orig,ystore_orig,zstore_orig

  integer, dimension(:), allocatable :: num_p
  integer, dimension(:), allocatable :: p_lookup

  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_new

  ! allocates temporary arrays
  allocate(num_p(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocate num_p,etc'
  num_p(:) = 0

  allocate(iglob_field_new(NGLOB_AB),iglob_field_new_cr(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocate iglob_p_field_new,etc'
  iglob_field_new(:) = 0
  iglob_field_new_cr(:) = 0.0

  allocate(xstore_orig(NGLOB_AB),ystore_orig(NGLOB_AB),zstore_orig(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocate iglob_p_field_new,etc'
  xstore_orig(:) = xstore(:)
  ystore_orig(:) = ystore(:)
  zstore_orig(:) = zstore(:)

  allocate(ibool_new(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool_new'
  ibool_new(:,:,:,:) = 0

  allocate(iglob_touched(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating iglob_touched'
  iglob_touched(:) = 0

  allocate(ibool_interfaces_ext_mesh_new(size(ibool_interfaces_ext_mesh,1),size(ibool_interfaces_ext_mesh,2)),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool_interfaces_ext_mesh_new'
  ibool_interfaces_ext_mesh_new(:,:) = 0

  allocate(p_lookup(maxval(p_level)),stat=ier)
  if (ier /= 0) stop 'Error allocating p_lookup'
  p_lookup(:) = 0

  ! 1. count number of elements in each p-level
  do ip = 1,num_p_level
    do iglob = 1,NGLOB_AB
      if (iglob_p_refine(iglob) == p_level(ip)) num_p(ip) = num_p(ip) + 1
    enddo
    ! lookup list to get p-level index (ilevel) from p-value (p)
    p_lookup(p_level(ip)) = ip
  enddo

  ! checks if every level has nodes
  !do ip = 1,num_p_level
  !  if (num_p(ip) == 0) then
  !    print *,'Error rank:',myrank,'level ',ip,' with p-value ',p_level(ip),'has no nodes!!! please check your mesh...'
  !    call exit_mpi(myrank,'Error setting up p-level nodes')
  !  endif
  !enddo

  ! checks total of p-nodes
  if (sum(num_p(:)) /= NGLOB_AB) then
    print *,'Error rank:',myrank,'has wrong number of p-nodes!!!',sum(num_p(:)),'instead of ',NGLOB_AB
    call exit_mpi(myrank,'Error setting up p-level nodes')
  endif

  ! use numbers of p-elements to build start-end array for dof
  ! note: for a level without p-nodes, start will be 1 bigger than end, e.g. start=1 end=0
  p_level_iglob_start(1) = 1
  p_level_iglob_end(1) = num_p(1)
  do ip = 2,num_p_level
    p_level_iglob_start(ip) = p_level_iglob_end(ip-1)+1
    p_level_iglob_end(ip) = p_level_iglob_start(ip) + num_p(ip)-1
  enddo

  ! checks if all nodes have been considered
  if (p_level_iglob_end(num_p_level) /= NGLOB_AB) then
    print *,'Error rank:',myrank,'p_level_iglob_end is invalid!!!',p_level_iglob_end(num_p_level),'instead of',NGLOB_AB
    print *,'p_level_iglob_end array:',p_level_iglob_end(:)
    call exit_mpi(myrank,'Error setting up p-level nodes')
  endif

  ! build reordering array "iglob_touched"
  num_p(:) = 0
  iglob_touched(:) = -1
  ibool_new(:,:,:,:) = -1
  do ispec = 1,NSPEC_AB
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          ! have not yet visited this dof (iglob)
          if (iglob_touched(iglob) == -1) then
            ! gets p-value for this node
            p = iglob_p_refine(iglob)
            ! gets p-level index
            ip = p_lookup(p)
            ! adds new iglob value
            ibool_new(i,j,k,ispec) = p_level_iglob_start(ip)+num_p(ip)
            num_p(ip) = num_p(ip) + 1
            ! mapping from old -> new
            iglob_touched(iglob) = ibool_new(i,j,k,ispec)
          else
            ! we have already visited this dof, but still need to update
            ! the ibool with the new iglob index
            ! (shared global node between different elements)
            ibool_new(i,j,k,ispec) = iglob_touched(iglob)
          endif
        enddo
      enddo
    enddo
  enddo

  ! checks if all nodes have been considered in new ibool array
  if (minval(ibool_new(:,:,:,:)) < 1 .or. maxval(ibool_new(:,:,:,:)) > NGLOB_AB) then
    print *,'Error rank:',myrank,'ibool_new array is invalid!!!',minval(ibool_new(:,:,:,:)),'is invalid entry'
    call exit_mpi(myrank,'Error setting up p-level nodes')
  endif

  ! re-checks total of p-nodes
  if (sum(num_p(:)) /= NGLOB_AB) then
    print *,'Error rank:',myrank,'has wrong number of p-nodes after reordering!!!',sum(num_p(:)),'instead of ',NGLOB_AB
    call exit_mpi(myrank,'Error setting up p-level nodes')
  endif

  ! use reordering array "iglob_touched" to reorder iglob arrays
  ! * iglob_p_refine
  ! * x,y,zstore
  iglob_field_new(:) = -1

  iglob_field_new(iglob_touched) = iglob_p_refine(1:NGLOB_AB)
  iglob_p_refine(:) = iglob_field_new(:)

  ! debug output
  !do i = 0,NPROC-1
  !  if (myrank == i) then
  !    print *,myrank,'iglob_p_refine:',iglob_p_refine(:)
  !  endif
  !  call synchronize_all()
  !enddo

  ! test: make sure iglob_touched mapping is 1-to-1.
  if (ANY(iglob_field_new(:) == -1)) stop 'Error: some iglobs not touched in reordered array!'

  iglob_field_new_cr(iglob_touched(:)) = xstore(1:NGLOB_AB)
  xstore(:) = iglob_field_new_cr(:)

  iglob_field_new_cr(iglob_touched(:)) = ystore(1:NGLOB_AB)
  ystore(:) = iglob_field_new_cr(:)

  iglob_field_new_cr(iglob_touched(:)) = zstore(1:NGLOB_AB)
  zstore(:) = iglob_field_new_cr(:)

  ! test x,y,z and new ibool
  do ispec = 1,NSPEC_AB
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          iglob_new = ibool_new(i,j,k,ispec)
          if (sqrt( (xstore(iglob_new) - xstore_orig(iglob))**2 + &
                    (ystore(iglob_new) - ystore_orig(iglob))**2 + &
                    (zstore(iglob_new) - zstore_orig(iglob))**2 ) > 1.e-4) then
            stop 'Error: (x,y,z) reordering didnt work!'
          endif
        enddo
      enddo
    enddo
  enddo

  ! copy new values into original arrays
  ibool(:,:,:,:) = ibool_new(:,:,:,:)

  ! fix MPI interface
  do iinterface = 1, num_interfaces_ext_mesh
    do i = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(i,iinterface)
      ibool_interfaces_ext_mesh_new(i,iinterface) = iglob_touched(iglob)
    enddo
  enddo
  ibool_interfaces_ext_mesh(:,:) = ibool_interfaces_ext_mesh_new(:,:)

  ! fix fault interfaces
  if (ANY_FAULT) then
    ! re-orders global values stored in ibulk1 & ibulk2 for fault split nodes
    if (ANY_FAULT_IN_THIS_PROC) call lts_fault_reorder_ibulk(iglob_touched)
  endif

  ! tests to make sure all arrays are correct
  if (ANY(iglob_touched(:) == -1)) stop 'Error: some iglobs not touched!'
  if (ANY(iglob_p_refine(:) < 1)) stop 'Error: some iglobs listed as p < 1!'
  if (ANY(ibool(:,:,:,:) < 1)) stop 'Error: some ibool still listed as -1!'

  deallocate(ibool_new)
  deallocate(iglob_field_new)
  deallocate(iglob_field_new_cr)
  deallocate(num_p)
  deallocate(iglob_touched)
  deallocate(xstore_orig,ystore_orig,zstore_orig)

  end subroutine lts_reorder_iglob_by_p_level

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_add_overlap_elements_databases(iglob_p_refine,ispec_p_refine)

  use constants, only: NGLLX,NGLLY,NGLLZ,IMAIN,myrank

  use generate_databases_par, only: &
    NSPEC_AB,NGLOB_AB,ibool

  implicit none

  integer,dimension(NGLOB_AB),intent(inout) :: iglob_p_refine
  integer,dimension(NSPEC_AB),intent(inout) :: ispec_p_refine

  ! local parameters
  integer :: ispec,iglob,i,j,k,p,icycle
  ! number of cycles to add overlap by one element layer
  integer, parameter :: NUMBER_OF_OVERLAP_CYCLES = 1

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     adding overlap region'
    write(IMAIN,*) '       number of overlap cycles = ', NUMBER_OF_OVERLAP_CYCLES
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! adds overlap by one element
  do icycle = 1,NUMBER_OF_OVERLAP_CYCLES
    ! flags all shared nodes according to the highest element p-level
    do ispec = 1,NSPEC_AB
      ! sets refinement on all points
      ! this enlarges the finer region by one element
      p = ispec_p_refine(ispec)
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            if (p > iglob_p_refine(iglob)) iglob_p_refine(iglob) = p
          enddo
        enddo
      enddo
    enddo
    ! re-sets p refinement for element due to overlap addition above
    ! (adds also elements touching global points with higher refinement)
    ispec_p_refine(:) = 0
    do ispec = 1,NSPEC_AB
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            p = iglob_p_refine(iglob)
            if (p > ispec_p_refine(ispec)) ispec_p_refine(ispec) = p
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine lts_add_overlap_elements_databases

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_save_databases()

  use constants, only: IOUT,IMAIN,MAX_STRING_LEN,myrank

  use shared_parameters, only: SAVE_MESH_FILES,LOCAL_PATH

  use generate_databases_par, only: NSPEC_AB,NGLOB_AB,ibool

  use create_regions_mesh_ext_par, only: &
    xstore => xstore_unique, &
    ystore => ystore_unique, &
    zstore => zstore_unique

  ! LTS module
  use lts_generate_databases_par

  implicit none

  integer :: ier

  ! output
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: prname

  ! VTK-file output for debugging
  logical,parameter :: DEBUG_VTK_OUTPUT = .false.

!#TODO: LTS database not stored yet in HDF5 / ADIOS2 format

  ! stores arrays in databases for solver
  call create_name_database(prname,myrank,LOCAL_PATH)

  filename = trim(prname) // 'lts.bin'

  open(unit=IOUT,file=trim(filename),status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening database proc******_lts.bin'

  write(IOUT) deltat_lts

  write(IOUT) num_p_level
  write(IOUT) p_level
  write(IOUT) p_level_loops

  write(IOUT) num_p_level_steps
  write(IOUT) p_level_steps

  write(IOUT) NGLOB_AB
  write(IOUT) iglob_p_refine

  write(IOUT) NSPEC_AB
  write(IOUT) ispec_p_refine

  write(IOUT) p_elem
  write(IOUT) boundary_elem

  write(IOUT) p_level_iglob_start
  write(IOUT) p_level_iglob_end

  close(IOUT)

  ! debug: for vtk output
  if (SAVE_MESH_FILES) then
    if (DEBUG_VTK_OUTPUT) then
      ! p-refinements
      filename = trim(prname) // 'ispec_p_refine'
      call write_VTK_data_elem_i(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,ispec_p_refine,filename)
      if (myrank == 0 ) then
        write(IMAIN,*) '     written file: ',trim(filename)//'.vtk'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif
  endif

  end subroutine lts_save_databases

!------------------------------------------------------------------------------------------------
!
! additional fault routines
!
!------------------------------------------------------------------------------------------------

  subroutine lts_fault_reorder_ibulk(iglob_touched)

  use generate_databases_par, only: NGLOB_AB

  use fault_generate_databases, only: fault_db_type,fault_db

  implicit none
  ! interface
  integer,dimension(NGLOB_AB),intent(in) :: iglob_touched

  ! local parameters
  ! fault
  type(fault_db_type) :: fdb
  integer :: iflt
  integer :: i,iglob1,iglob2

  ! loops over all faults
  do iflt = 1,size(fault_db)
    ! fault
    fdb = fault_db(iflt)
    ! re-assignes iglob values for p-ordered ibool,xstore_dummy,.. arrays
    do i = 1, fdb%nglob
      ! nodes on fault surface 1
      iglob1 = fdb%ibulk1(i)
      fdb%ibulk1(i) = iglob_touched(iglob1)
      ! nodes on fault surface 2
      iglob2 = fdb%ibulk2(i)
      fdb%ibulk2(i) = iglob_touched(iglob2)
    enddo
  enddo ! iflt

  end subroutine lts_fault_reorder_ibulk

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_fault_assign_split_p_values()

  use fault_generate_databases, only: fault_db_type,fault_db

  ! LTS module
  use lts_generate_databases_par, only: iglob_p_refine

  implicit none

  ! local parameters
  ! fault
  type(fault_db_type) :: fdb
  integer :: iflt
  integer :: i,j,iglob1,iglob2,p1,p2

  ! loops twice, note below:
  !JPA do it twice, to handle triple junctions (intersections between two faults)
  do i = 1,2
    ! loops over all faults
    do iflt = 1,size(fault_db)
      ! fault
      fdb = fault_db(iflt)
      ! re-assignes p-values on split nodes, takes maximum p-value for both
      do j = 1, fdb%nglob
        ! nodes on fault surface 1
        iglob1 = fdb%ibulk1(j)
        p1 = iglob_p_refine(iglob1)
        ! nodes on fault surface 2
        iglob2 = fdb%ibulk2(j)
        p2 = iglob_p_refine(iglob2)
        ! sets higher p-value for both split nodes
        iglob_p_refine(iglob1) = max(p1,p2)
        iglob_p_refine(iglob2) = max(p1,p2)
      enddo
    enddo ! iflt
  enddo

  end subroutine lts_fault_assign_split_p_values
