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

! Local Time Stepping (LTS)
!
! In case you use this local time stepping feature in your study, please reference this work:
!
! Rietmann, M., M. Grote, D. Peter, O. Schenk, 2017
! Newmark local time stepping on high-performance computing architectures,
! Journal of Comp. Physics, 334, p. 308-326.
! https://doi.org/10.1016/j.jcp.2016.11.012
!
! Rietmann, M., B. Ucar, D. Peter, O. Schenk, M. Grote, 2015.
! Load-balanced local time stepping for large-scale wave propagation,
! in: Parallel Distributed Processing Symposium (IPDPS), IEEE International, May 2015.
! https://doi.org/10.1109/IPDPS.2015.10


! LTS setup routines
!
! Authors: Max Rietmann, Daniel Peter

  subroutine lts_setup()

  use constants, only: NDIM,CUSTOM_REAL,VERYSMALLVAL,FIX_UNDERFLOW_PROBLEM,IMAIN,itag,myrank

  use shared_parameters, only: DT,NSTEP,NPROC,SIMULATION_TYPE,LTS_MODE,GPU_MODE,USE_LDDRK, &
    SAVE_SEISMOGRAMS_ACCELERATION,CREATE_SHAKEMAP

  use specfem_par, only: ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION, &
    deltat,NGLOB_AB

  ! MPI interfaces
  use specfem_par, only: num_interfaces_ext_mesh,ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
    my_neighbors_ext_mesh

  !use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  !use specfem_par_elastic, only: ELASTIC_SIMULATION
  !use specfem_par_poroelastic, only: POROELASTIC_SIMULATION
  !use specfem_par_movie

  use specfem_par_lts

  implicit none

  ! local parameters
  integer :: ilevel,p,m,step
  integer :: p_min,p_max,p_min_glob,p_max_glob,p_node
  integer :: i, ipoin, iinterface, iglob, ier
  real(kind=CUSTOM_REAL) :: duration
  integer, dimension(:), allocatable :: num_boundary_p, num_p_glob
  ! check
  integer, dimension(:,:), allocatable ::nibool_interface_p_refine_all
  integer, dimension(:), allocatable :: sendbuf,recvbuf
  double precision :: lts_speedup,lts_speedup_with_boundary
  double precision :: deltat_lts
  double precision :: memory_size
  integer :: is,ie,is0,ie0
  ! debugging
  integer :: iproc

  ! checks if anything to do
  if (.not. LTS_MODE) return

  ! safety checks
  ! checks simulation domain types
  if (ACOUSTIC_SIMULATION) call exit_MPI(myrank,'LTS routines only implemented for purely ELASTIC simulations')
  if (POROELASTIC_SIMULATION) call exit_MPI(myrank,'LTS routines only implemented for purely ELASTIC simulations')
  if (.not. ELASTIC_SIMULATION) call exit_MPI(myrank,'LTS routines only implemented for purely ELASTIC simulations')

  if (SIMULATION_TYPE /= 1) call exit_MPI(myrank,'LTS routines only implemented for forward simulations (SIMULATION_TYPE == 1)')
  if (USE_LDDRK) call exit_MPI(myrank,'LTS routines only implemented for Newark time scheme (USE_LDDRK == .false.)')
  if (GPU_MODE) call exit_MPI(myrank,'LTS routines only implemented for CPU-only runs (GPU_MODE == .false.)')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Local time stepping:'
    call flush_IMAIN()
  endif

  ! reads in local time stepping arrays from databases files
  call lts_read_databases()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of p-levels     : ',num_p_level
    write(IMAIN,*) '  p level                : ',p_level(:)
    write(IMAIN,*) '  p level loops          : ',p_level_loops(:)
    write(IMAIN,*)
    write(IMAIN,*) '  number of p-level steps: ',num_p_level_steps
    write(IMAIN,*) '  p level steps          : ',p_level_steps(:)
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! safety check
  if (num_p_level == 0) stop 'Error LTS mode needs at least a single p-level'

  ! we assume contiguous range of p-level nodes
  ! checks iglob range
  if (p_level_iglob_start(1) /= 1) call exit_mpi(myrank,"ASSERT: p-level iglob should start at 1")
  if (p_level_iglob_end(num_p_level) /= NGLOB_AB) call exit_mpi(myrank,"ASSERT: coarsest p-level must have all iglobs!")
  ! checks levels
  do ilevel = 1,num_p_level
    ! start index of current p-level
    is = p_level_iglob_start(ilevel)

    ! end index of current p-level
    ie = p_level_iglob_end(ilevel)

    ! start index of finest level
    is0 = p_level_iglob_start(1)

    ! end index of next finer p-level
    if (ilevel > 1) then
      ie0 = p_level_iglob_end(ilevel-1)
    endif

    ! checks
    if (ie < is) stop 'Error lts newmark update: end index of current level invalid'
    if (ie < is0) stop 'Error lts newmark update: start/end index invalid'
    if (ilevel > 1 .and. ie0 < is0) stop 'Error lts newmark update: end index of finer level invalid'
  enddo

  ! desired initial duration
  duration = sngl(NSTEP * DT)

  ! gets min/max over all processes
  p_min = minval(p_level(:))
  p_max = maxval(p_level(:))
  call min_all_all_i(p_min,p_min_glob)
  call max_all_all_i(p_max,p_max_glob)
  p_min = p_min_glob
  p_max = p_max_glob

  ! set "global" DT (coarsest step)
  ! DT as given in Par_file is the minimum step in smallest element
  deltat_lts = DT * dble(p_max)
  deltat = deltat_lts

  ! counts number of local time step evaluations
  lts_it_local = 0
  do step = 1,num_p_level_steps
    ilevel = p_level_steps(step)
    p = p_level(ilevel)
    do m = 1,p_level_loops(ilevel)
      lts_it_local = lts_it_local + 1
    enddo
  enddo
  NSTEP_LOCAL = lts_it_local

  ! gets theoretical speed-up values (valid only for all elements belonging to same domain type)
  call lts_get_theoretical_speedup(lts_speedup,lts_speedup_with_boundary)

  ! a map for p-level -> ilevel
  allocate(p_level_ilevel_map(p_level(1)),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating working LTS fields p_level_ilevel_map')
  p_level_ilevel_map(:) = 0
  do ilevel = 1,num_p_level
    p_level_ilevel_map(p_level(ilevel)) = ilevel
  enddo

  ! p_lookup maps p -> ilevel
  allocate(p_lookup(maxval(p_level)),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating p_lookup')
  p_lookup(:) = 0
  do ilevel = 1,num_p_level
    p_lookup(p_level(ilevel)) = ilevel
  enddo

  ! Find p-level information for MPI-boundaries.
  allocate(interface_p_refine_all(num_interfaces_ext_mesh,num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating interfrace_p_refine_all')
  interface_p_refine_all(:,:) = -1

  allocate(num_boundary_p(num_p_level),num_p_glob(num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating num_boundary_p, num_p_glob')
  num_boundary_p(:) = 0
  num_p_glob(:) = 0

  allocate(nibool_interface_p_refine_all(num_interfaces_ext_mesh,num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating nibool_interfrace_p_refine_all')
  nibool_interface_p_refine_all(:,:) = 0

  ! debug
  !do iproc = 0,NPROC-1
  !call synchronize_all()
  !if (myrank == iproc) then

  do iinterface = 1, num_interfaces_ext_mesh
    do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
      ! gets p value for this interface points
      iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
      p_node = iglob_p_refine(iglob)

      ! debug
      ! if (myrank == 1) print *, "p_node=", p_node

      ! for each interface, stores for each ilevel its corresponding p value if p-nodes on this interface where found
      interface_p_refine_all(iinterface,p_lookup(p_node)) = p_node
      num_boundary_p(p_lookup(p_node)) = num_boundary_p(p_lookup(p_node)) + 1

      ! counts p-nodes on interface
      nibool_interface_p_refine_all(iinterface,p_lookup(p_node)) = nibool_interface_p_refine_all(iinterface,p_lookup(p_node)) + 1
    enddo
    ! debug
    !print *,'rank ',myrank,'interface neighbor',my_neighbors_ext_mesh(iinterface)
    !print *,'  p_refine: ',interface_p_refine_all(iinterface,:)
    !print *,'  p_refine: ',nibool_interface_p_refine_all(iinterface,:)
  enddo

  ! debug
  !endif
  !enddo

  ! debug output
  !do iproc = 0,NPROC-1
  !  call synchronize_all()
  !  if (myrank == iproc) then
  !    print *,'rank',myrank,' num boundary p:'
  !    do ilevel = 1,num_p_level
  !      print *,'  level ',ilevel,' total number of p-nodes on MPI boundary = ',num_boundary_p(ilevel)
  !    enddo
  !  endif
  !enddo

  ! MPI setup
  if (NPROC > 1) then
    ! checks number match with other partitions
    allocate(sendbuf(num_p_level),recvbuf(num_p_level),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating temporary sendbuf/recvbuf arrays')

    do iinterface = 1, num_interfaces_ext_mesh
      ! collects nibools from other process
      sendbuf(:) = nibool_interface_p_refine_all(iinterface,:)
      recvbuf(:) = 0
      call sendrecv_all_i(sendbuf,num_p_level,my_neighbors_ext_mesh(iinterface),itag, &
                          recvbuf,num_p_level,my_neighbors_ext_mesh(iinterface),itag)
      ! checks if number of p-nodes on interface match
      do ilevel = 1,num_p_level
        if (sendbuf(ilevel) /= recvbuf(ilevel)) then
          print *,'Error p-node counts: rank',myrank,'neighbor:',my_neighbors_ext_mesh(iinterface)
          print *,'  send nibools:',sendbuf(:)
          print *,'  recv nibools:',recvbuf(:)
          print *,'rank',myrank,' num boundary p:'
          do i = 1,num_p_level
            print *,'  level ',i,' total number of p-nodes on MPI boundary = ',num_boundary_p(i)
          enddo
          call exit_MPI(myrank,'Error send & receive nibool buffers')
        endif
      enddo
    enddo
    ! frees temporary arrays
    deallocate(sendbuf,recvbuf)
  endif

  !debug
  ! num_p_glob(:) = 0
  ! do iglob = 1,NGLOB_AB
  !   p_node = iglob_p_refine(iglob)
  !   num_p_glob(p_lookup(p_node)) = num_p_glob(p_lookup(p_node)) + 1
  ! enddo
  ! if (myrank==0) print *, "Rank, p, #dof"
  ! do ilevel = 1,num_p_level
  !   if (num_boundary_p(ilevel) > 0) then
  !     print *, myrank, ",", p_level(ilevel), ",", num_p_glob(ilevel)
  !   endif
  ! enddo
  ! debug user output
  !call synchronize_all()
  !if (myrank == 0) then
  !  print *,'number of p-nodes on interfaces match'
  !endif

  ! checks if iglob values are ordered with respect to p_level_iglob_** start/end arrays
  do iproc = 0,NPROC-1
    if (myrank == iproc) then
      ! loops over all global points
      do iglob = 1,NGLOB_AB
        ! gets p value for this node
        p_node = iglob_p_refine(iglob)

        ! gets ilevel index for this p-value
        ilevel = p_lookup(p_node)

        ! start/end index for this p-level
        is = p_level_iglob_start(ilevel)
        ie = p_level_iglob_end(ilevel)

        ! checks if iglob is in start/end range
        if (iglob < is .or. iglob > ie) then
          print *,'Error rank:',myrank,'has iglob value ',iglob,'in level',ilevel,'out of range:',is,ie
          call exit_MPI(myrank,'Error iglob ordering for LTS invalid')
        endif
      enddo
      ! loops over all start/end index ranges
      do ilevel = 1,num_p_level
        ! start/end index for this p-level
        is = p_level_iglob_start(ilevel)
        ie = p_level_iglob_end(ilevel)
        ! checks that points in this range belong to ilevel
        do iglob = is,ie
          ! checks iglob value
          if (iglob < 1 .or. iglob > NGLOB_AB) then
            print *,'Error rank:',myrank,'has iglob value ',iglob,'in level',ilevel,'range:',is,ie
            call exit_MPI(myrank,'Error iglob invalid in LTS level')
          endif
          ! checks associated p-level
          p_node = iglob_p_refine(iglob)
          i = p_lookup(p_node)
          if (i /= ilevel) then
            print *,'Error rank:',myrank,'has iglob value ',iglob,'in level',i,'range:',is,ie,'p-value:',p_node,'ilevel:',ilevel
            call exit_MPI(myrank,'Error iglob level invalid in LTS level')
          endif
        enddo
      enddo
    endif
    call synchronize_all()
  enddo

  ! frees memory
  deallocate(nibool_interface_p_refine_all)
  deallocate(num_boundary_p)

  ! sets up p-level boundaries
  call lts_setup_level_boundary()

  ! keeps track of current local time step
  allocate(lts_current_m(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating lts_current_m array'
  lts_current_m(:) = 0

  ! initializes current pointers
  current_lts_elem => p_elem(:,1)                 ! current p-level elements
  current_lts_boundary_elem => boundary_elem(:,1) ! current boundary elements

  ! current p-level time
  current_lts_time = 0.d0

  ! user output
  if (myrank == 0) then
    ! wavefield sizes displ_p,veloc_p
    memory_size = dble(NDIM) * dble(NGLOB_AB) * dble(num_p_level) * dble(CUSTOM_REAL)
    write(IMAIN,*) '  array size for working fields displ_p, veloc_p: ', sngl(memory_size) / 1024./ 1024.,'MB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! working fields
  allocate(displ_p(NDIM,NGLOB_AB,num_p_level), &
           veloc_p(NDIM,NGLOB_AB,num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating working LTS fields displ_p,veloc_p')
  displ_p(:,:,:) = 0.0_CUSTOM_REAL
  veloc_p(:,:,:) = 0.0_CUSTOM_REAL
  ! put negligible initial value to avoid very slow underflow trapping
  if (FIX_UNDERFLOW_PROBLEM) displ_p(:,:,:) = VERYSMALLVAL

  ! collected acceleration
  ! only needed for seismograms and shakemaps
  if (SAVE_SEISMOGRAMS_ACCELERATION .or. CREATE_SHAKEMAP) then
    allocate(accel_collected(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating working LTS fields displ_p,veloc_p')
    accel_collected(:,:) = 0.0_CUSTOM_REAL
  endif

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  original time step                 : ',sngl(DT)
    write(IMAIN,*) '           number of time steps      : ',NSTEP
    write(IMAIN,*) '           total duration            : ',sngl(NSTEP * DT)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! re-sets simulation NSTEP
  NSTEP = ceiling( duration / deltat_lts )

  ! re-sets simulation global time step to be coarsest-level time step
  DT = deltat_lts

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  new LTS  time step                 : ',sngl(DT)
    write(IMAIN,*) '           number of time steps      : ',NSTEP
    write(IMAIN,*) '           total duration            : ',sngl(NSTEP * DT)
    write(IMAIN,*) '           number of local time steps: ',NSTEP_LOCAL
    write(IMAIN,*)
    write(IMAIN,*) '  theoretical speed-up value: ',sngl(lts_speedup),'(without boundary contributions)'
    write(IMAIN,*) '  theoretical speed-up value: ',sngl(lts_speedup_with_boundary),'(with boundary contributions)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  local time stepping setup done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine lts_setup

!
!------------------------------------------------------------------------------------------------
!

  subroutine lts_read_databases()

  use constants, only: MAX_STRING_LEN,IIN,IMAIN,myrank
  use shared_parameters, only: LOCAL_PATH

  use specfem_par, only: NGLOB_AB,NSPEC_AB

  use specfem_par_lts

  implicit none

  integer :: itmp,ier
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: prname

  ! reads in local time stepping arrays from databases files
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading in local time stepping arrays'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  call create_name_database(prname,myrank,LOCAL_PATH)
  filename = trim(prname) // 'lts.bin'

  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open database file : ',trim(filename)
    print *,'Please check if databases were generated with LTS_MODE set to .true. in Par_file'
    stop 'Error opening database proc******_lts.bin'
  endif

  ! suggested coarsest time step for LTS (largest multiple p of smallest time step)
  read(IIN) deltat_lts_suggested

  ! levels
  read(IIN) num_p_level

  ! allocates arrays
  allocate(p_level(num_p_level), &
           p_level_loops(num_p_level), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array p_level and p_level_loops')
  p_level(:) = 0
  p_level_loops(:) = 0

  read(IIN) p_level
  read(IIN) p_level_loops

  ! loops
  read(IIN) num_p_level_steps

  allocate(p_level_steps(num_p_level_steps),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array p_level_steps')
  p_level_steps(:) = 0

  read(IIN) p_level_steps

  ! p-refinement values on global points
  read(IIN) itmp !nglob
  if (itmp /= NGLOB_AB) call exit_MPI(myrank,'Error reading nglob in LTS databases file')

  allocate(iglob_p_refine(NGLOB_AB), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array iglob_p_refine')
  iglob_p_refine(:) = 0

  read(IIN) iglob_p_refine

  ! element flags
  read(IIN) itmp ! nspec
  if (itmp /= NSPEC_AB) call exit_MPI(myrank,'Error reading nspec in LTS databases file')

  allocate(ispec_p_refine(NSPEC_AB), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array ispec_p_refine')
  ispec_p_refine(:) = 0

  read(IIN) ispec_p_refine

  allocate(p_elem(NSPEC_AB,num_p_level), &
           boundary_elem(NSPEC_AB,num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array p_elem and boundary_elem')
  p_elem(:,:) = .false.
  boundary_elem(:,:) = .false.

  read(IIN) p_elem
  read(IIN) boundary_elem

  allocate(p_level_iglob_start(num_p_level), p_level_iglob_end(num_p_level), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array p_level_start/end')
  p_level_iglob_start(:) = 0
  p_level_iglob_end(:) = 0

  read(IIN) p_level_iglob_start
  read(IIN) p_level_iglob_end

  close(IIN)

  end subroutine lts_read_databases

!
!--------------------------------------------------------------------------------------------
!

  subroutine lts_get_theoretical_speedup(lts_speedup,lts_speedup_with_boundary)

  use constants, only: IMAIN,myrank
  use shared_parameters, only: NPROC
  use specfem_par, only: NSPEC_AB

  use specfem_par_lts

  implicit none

  double precision, intent(out) :: lts_speedup,lts_speedup_with_boundary

  ! local parameters
  double precision :: total_work,total_work_lts,total_work_lts_b
  integer,dimension(:),allocatable :: num_ispec_level
  integer :: ispec,p,ilevel,inum,ier

  ! theoretical speed-up values (valid only for all elements belonging to same domain type)
  total_work = 0.0
  total_work_lts = 0.0

  ! counts number of elements for each level
  allocate(num_ispec_level(num_p_level),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating array num_ispec_level')
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
  if (minval(num_ispec_level(:)) == 0) call exit_MPI(myrank,'Error p level without element found')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  p-level number of elements: ',num_ispec_level(:)
    call flush_IMAIN()
  endif

  ! computes total work
  if (myrank == 0) then
    ! "pure" speed-up, without taking account of additional coarse/fine boundary contributions
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
  if (NPROC > 1 ) then
    do ilevel = 1,num_p_level
      call sum_all_all_i(num_ispec_level(ilevel),inum)
      num_ispec_level(ilevel) = inum
    enddo
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  p-level boundary elements : ',num_ispec_level(:)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! work with additional coarse/fine boundary contributions
  total_work_lts_b = total_work_lts
  if (myrank == 0) then
    do ilevel = 1,num_p_level
      p = p_level(ilevel)
      total_work_lts_b = total_work_lts_b + num_ispec_level(ilevel) * p
    enddo
  endif

  ! speedup factors
  if (total_work_lts /= 0.d0) then
    lts_speedup = total_work / total_work_lts
  else
    lts_speedup = 0.d0
  endif

  ! speedup factor w/ boundary work
  if (total_work_lts_b /= 0.d0) then
    lts_speedup_with_boundary = total_work / total_work_lts_b
  else
    lts_speedup_with_boundary = 0.d0
  endif

  ! free temporary array
  deallocate(num_ispec_level)

  end subroutine lts_get_theoretical_speedup

!
!--------------------------------------------------------------------------------------------
!

  subroutine lts_setup_level_boundary()

  use constants, only: NDIM,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IMAIN,itag,myrank

  use shared_parameters, only: NPROC

  use specfem_par, only: NGLOB_AB,NSPEC_AB,ibool, &
    num_interfaces_ext_mesh,nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh, &
    my_neighbors_ext_mesh,request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
    max_nibool_interfaces_ext_mesh, &
    xstore,ystore,zstore

  use specfem_par_elastic, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  use specfem_par_lts, only: &
    p_level_ilevel_map, num_p_level, boundary_elem, iglob_p_refine, &
    num_p_level_boundary_ispec,num_p_level_boundary_nodes, p_level_boundary_node, &
    p_level_boundary_ispec, &
    p_level_boundary_ilevel_from, &
    num_p_level_coarser_to_update, p_level_coarser_to_update, &
    num_interface_p_refine_ibool, interface_p_refine_ibool, &
    p_level_iglob_start, p_level_iglob_end, &
    p_level, p_level_m_loops, current_lts_boundary_elem

  implicit none

  ! needed for qsort routine
  ! this drives the qsort subroutine
  !external cmp_function
  !integer :: cmp_function

  ! local parameters
  integer :: ilevel, iphase, ispec_p, ispec, num_elements, iglob
  integer :: jlevel, iglob_n, ipoin, iinterface, ipoin_n, ier, ipoin_old
  integer :: p,m
  integer :: i,j,k,icounter,inum_poin,inum_spec
  integer :: coarser_counter, coarser_counter_new, num_ipoin_coarser_neighbor
  integer :: is,ie
  integer :: p_node
  integer, dimension(:), allocatable :: interface_p_refine_ipoin
  integer, dimension(:), allocatable :: ipoin_coarser
  integer, dimension(:), allocatable :: recv_interface_p_refine_ipoin

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyz_send,xyz_recv
  integer :: num_points,num_points_recv

  integer, dimension(:,:), allocatable :: p_level_coarser_nodes
  integer, dimension(:,:), allocatable :: tmp_p_level_coarser_nodes
  integer, dimension(:), allocatable :: num_p_level_coarser_nodes

  logical, dimension(:), allocatable :: mask_ibool

  ! timing
  logical, parameter :: DEBUG_TIMING = .false.
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! p-level MPI interface testing
  logical, parameter :: TEST_ERROR = .true.

  integer :: mpi_cost
  integer, dimension(:), allocatable :: mpi_cost_gather

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  p-level boundary setup:'
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! get MPI starting time
  if (DEBUG_TIMING) time_start = wtime()

  ! boundary arrays
  allocate(p_level_boundary_ispec(NSPEC_AB,2,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays p_level_boundary_ispec,...'
  p_level_boundary_ispec(:,:,:) = 0

  allocate(p_level_boundary_node(NGLOB_AB,2,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays p_level_boundary_node,...'
  p_level_boundary_node(:,:,:) = 0

  allocate(p_level_boundary_ilevel_from(NGLOB_AB,2,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays p_level_boundary_ilevel_from,...'
  p_level_boundary_ilevel_from(:,:,:) = 0

  ! counters
  allocate(num_p_level_boundary_nodes(2,num_p_level), &
           num_p_level_boundary_ispec(2,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays num_p_level_boundary_nodes,...'
  num_p_level_boundary_nodes(:,:) = 0
  num_p_level_boundary_ispec(:,:) = 0

  ! temporary arrays
  allocate(mask_ibool(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays mask_ibool'
  mask_ibool(:) = .false.

  ! temporary arrays for nodes in coarser p-levels
  allocate(p_level_coarser_nodes(NGLOB_AB,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays p_level_coarser_nodes,...'
  p_level_coarser_nodes(:,:) = 0

  allocate(num_p_level_coarser_nodes(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays num_p_level_coarser_nodes,...'
  num_p_level_coarser_nodes(:) = 0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "    determining coarser p-level nodes"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! sets up coarser nodes
  num_p_level_coarser_nodes(:) = 0
  p_level_coarser_nodes(:,:) = 0

  do ilevel = 1,num_p_level

    ! current boundary elements
    current_lts_boundary_elem => boundary_elem(:,ilevel)

    do iphase = 1,2

      num_p_level_boundary_nodes(iphase,ilevel) = 0
      num_p_level_boundary_ispec(iphase,ilevel) = 0

      ! choses inner/outer elements
      if (iphase == 1) then
        num_elements = nspec_outer_elastic
      else
        num_elements = nspec_inner_elastic
      endif

      ! counters
      inum_spec = 0
      inum_poin = 0
      mask_ibool(:) = .false.

      do ispec_p = 1,num_elements

        ispec = phase_ispec_inner_elastic(ispec_p,iphase)

        ! only elements belonging to a p-level boundary
        if (current_lts_boundary_elem(ispec)) then
          ! boundary elements
          inum_spec = inum_spec + 1

          num_p_level_boundary_ispec(iphase,ilevel) = inum_spec
          p_level_boundary_ispec(inum_spec,iphase,ilevel) = ispec

          do k = 1,NGLLZ
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)

                ! associated p-level value
                p_node = iglob_p_refine(iglob)

                ! add global point only once
                if (.not. mask_ibool(iglob)) then
                  mask_ibool(iglob) = .true.
                  ! counter
                  inum_poin = inum_poin + 1

                  num_p_level_boundary_nodes(iphase,ilevel) = inum_poin
                  p_level_boundary_node(inum_poin,iphase,ilevel) = iglob

                  p_level_boundary_ilevel_from(inum_poin,iphase,ilevel) = p_level_ilevel_map(p_node)

                  ! adds node to coarser nodes
                  if (p_level_ilevel_map(p_node) /= ilevel) then
                    num_p_level_coarser_nodes(ilevel) = num_p_level_coarser_nodes(ilevel) + 1
                    p_level_coarser_nodes(num_p_level_coarser_nodes(ilevel),ilevel) = iglob
                  endif
                endif
              enddo
            enddo
          enddo
        endif
      enddo ! ispec
    enddo ! iphase
  enddo ! ilevel

  ! gets maximum number of coarser nodes
  if (myrank == 0) then
    write(IMAIN,*) "    maximum coarser nodes: ",maxval(num_p_level_coarser_nodes(:))
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! stores compressed array to save space
  allocate(tmp_p_level_coarser_nodes(maxval(num_p_level_coarser_nodes(:)),num_p_level-1),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays tmp_p_level_coarser_nodes'
  tmp_p_level_coarser_nodes(:,:) = 0

  do ilevel = 1,num_p_level-1
    do iglob_n = 1,num_p_level_coarser_nodes(ilevel)
      tmp_p_level_coarser_nodes(iglob_n,ilevel) = p_level_coarser_nodes(iglob_n,ilevel)
    enddo
  enddo

  deallocate(p_level_coarser_nodes)

  ! ----------------------------------- !
  ! --  Setup time-stepping arrays ---- !

  ! sets up all nodes from coarser contributions
  ! nodes in coarser p-levels to used in MPI update
  allocate(p_level_coarser_to_update(NGLOB_AB,num_p_level), &
           num_p_level_coarser_to_update(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays p_level_coarser_to_update,...'
  num_p_level_coarser_to_update(:) = 0
  p_level_coarser_to_update(:,:) = 0

  ! loops over finer p-levels
  do ilevel = 1,(num_p_level-1)
    ! checks if nodes have been taken multiple times from coarser p-boundary already
    mask_ibool(:) = .false.
    ! loops over all levels up to current one
    do jlevel = 1,ilevel
      do iglob_n = 1,num_p_level_coarser_nodes(jlevel)
        ! checks if node belongs to a coarser level
        iglob = tmp_p_level_coarser_nodes(iglob_n,jlevel)
        p_node = iglob_p_refine(iglob)

        if (p_level_ilevel_map(p_node) > ilevel .and. (.not. mask_ibool(iglob))) then
          ! adds node
          icounter = num_p_level_coarser_to_update(ilevel)
          icounter = icounter + 1

          num_p_level_coarser_to_update(ilevel) = icounter
          p_level_coarser_to_update(icounter,ilevel) = iglob

          mask_ibool(iglob) = .true.
        endif
      enddo
    enddo
  enddo

  ! frees temporary arrays
  deallocate(tmp_p_level_coarser_nodes)
  deallocate(num_p_level_coarser_nodes)
  deallocate(mask_ibool)

  ! user output
  if (DEBUG_TIMING) then
    if (myrank == 0 ) then
      tCPU = wtime() - time_start
      write(IMAIN,*) "    time in seconds = ", tCPU
      call flush_IMAIN()
    endif
    time_start = wtime()
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "    building MPI arrays"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! --- MPI setup for LTS ------ !
  allocate(num_interface_p_refine_ibool(num_interfaces_ext_mesh,num_p_level), &
           interface_p_refine_ibool(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh,num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating MPI arrays num_interface_p_refine_ibool,..'
  num_interface_p_refine_ibool(:,:) = 0
  interface_p_refine_ibool(:,:,:) = 0

  ! temporary arrays
  allocate(ipoin_coarser(max_nibool_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ipoin_coarser,...'
  ipoin_coarser(:) = 0

  ! coarsest level (ilevel == num_p_level) sends full boundary
  do iinterface = 1, num_interfaces_ext_mesh
    num_interface_p_refine_ibool(iinterface,num_p_level) = nibool_interfaces_ext_mesh(iinterface)
    interface_p_refine_ibool(:,iinterface,num_p_level) = ibool_interfaces_ext_mesh(:,iinterface)
  enddo

  if (TEST_ERROR) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      testing LTS initial MPI arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      checking initial interface locations"
      call flush_IMAIN()
    endif

    ! checks that locations of interace points are identical between MPI neighbors (in different partitions)
    ! (ilevel == num_p_level)

    ! loops over all interfaces (with other partitions)
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of interface points
      num_points = num_interface_p_refine_ibool(iinterface,num_p_level)
      ! checks locations of interface points
      allocate(xyz_send(NDIM,num_points),xyz_recv(NDIM,num_points),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating temporary array xyz_send')
      xyz_send(:,:) = 0._CUSTOM_REAL
      xyz_recv(:,:) = 1.e9_CUSTOM_REAL

      ! fills in locations
      do iglob_n = 1,num_points
        iglob = interface_p_refine_ibool(iglob_n,iinterface,num_p_level)
        ! checks index
        if (iglob < 1 .or. iglob > NGLOB_AB) then
          print *,'Error rank',myrank,' : iglob index ',iglob,' out of range ',NGLOB_AB
          call exit_mpi(myrank,'Error p-level MPI interface points range')
        endif
        ! fills send array
        xyz_send(1,iglob_n) = xstore(iglob)
        xyz_send(2,iglob_n) = ystore(iglob)
        xyz_send(3,iglob_n) = zstore(iglob)
      enddo

      ! sends interface point locations
      call isend_cr(xyz_send,NDIM*num_points, &
                    my_neighbors_ext_mesh(iinterface), &
                    itag,request_send_vector_ext_mesh(iinterface))
      ! gets the number of points which will be received from neighbor
      call get_count_i(my_neighbors_ext_mesh(iinterface),itag,num_points_recv)

      ! checks that sent/received number of points are equal
      if (NDIM*num_points /= num_points_recv) then
        print *,'Error rank',myrank,': number of interface points differ!!! ', &
                'ilevel:',num_p_level,'interface:',iinterface,'neighbor:',my_neighbors_ext_mesh(iinterface), &
                'num_points:',num_points,'num_points_recv:',num_points_recv,'per component',num_points_recv/NDIM
        call exit_MPI(myrank,'Error p-level MPI interface points')
      endif

      ! receives locations
      call irecv_cr(xyz_recv,NDIM*num_points, &
                    my_neighbors_ext_mesh(iinterface), &
                    itag,request_recv_vector_ext_mesh(iinterface))
      call wait_req(request_recv_vector_ext_mesh(iinterface))
      call wait_req(request_send_vector_ext_mesh(iinterface))

      ! checks that locations are identical
      do iglob_n = 1,num_points
        if (maxval(abs(xyz_send(:,iglob_n)-xyz_recv(:,iglob_n))) > 1.e-9) then
          print *,'Error rank',myrank,': level',num_p_level,'point:',iglob_n, &
                'send:',xyz_send(:,iglob_n),'recv:',xyz_recv(:,iglob_n), &
                'diff:',abs(xyz_send(:,iglob_n)-xyz_recv(:,iglob_n))
          call exit_mpi(myrank,"Error: (x,y,z) of MPI don't match")
        endif
      enddo
      deallocate(xyz_send)
      deallocate(xyz_recv)
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      test result okay"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (DEBUG_TIMING) time_start = wtime()
  endif ! TEST_ERROR

  ! temporary array of nodes on an MPI interface for refinement
  allocate(interface_p_refine_ipoin(max_nibool_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'Error allocating array interface_p_refine_ipoin'
  interface_p_refine_ipoin(:) = 0

  ! setup all levels but coarsest
  do ilevel = 1,num_p_level-1
    ! user output
    if (myrank == 0 ) then
      write(IMAIN,*) "    MPI interfaces for finer p-levels: level = ",ilevel," p-value = ",p_level(ilevel)
      call flush_IMAIN()
    endif
    if (DEBUG_TIMING) time_start = wtime()

    ! sets up MPI transfer interfaces for lts refinement
    do iinterface = 1, num_interfaces_ext_mesh

      ! determines all MPI interface nodes in a coarser p-level
      interface_p_refine_ipoin(:) = 0
      icounter = 0
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)

        ! checks if node belongs to finer/current p-levels
        if (iglob <= p_level_iglob_end(ilevel)) then
          ! if (iglob > -1) then
          ! num_interface_p_refine_ibool(iinterface,ilevel) = num_interface_p_refine_ibool(iinterface,ilevel) + 1
          ! interface_p_refine_ibool(num_interface_p_refine_ibool(iinterface,ilevel),iinterface,ilevel) = iglob
          icounter = icounter + 1
          interface_p_refine_ipoin(icounter) = ipoin
        else
          ! node belongs to coarser levels (on p-level boundary)
          if (count(p_level_coarser_to_update(1:num_p_level_coarser_to_update(ilevel),ilevel) == iglob) > 0) then
            icounter = icounter + 1
            interface_p_refine_ipoin(icounter) = ipoin
          endif
        endif
      enddo
      coarser_counter = icounter

      ! we now have the coarser points that are on the
      ! mpi-partition-boundary. The neighbor can have different
      ! coarser points, and we need to update all of them in the
      ! MPI-sync. We transfer our points to the neighbor, and
      ! receive "his" points, and compute their union, so that both
      ! sides send all needed points.
      if (ilevel < num_p_level) then
        ! finer p-levels
        ! send ipoin indices
        call isend_i(interface_p_refine_ipoin,coarser_counter, &
                     my_neighbors_ext_mesh(iinterface), &
                     itag,request_send_vector_ext_mesh(iinterface))

        ! gets number coarser points which will be received
        call get_count_i(my_neighbors_ext_mesh(iinterface),itag,num_ipoin_coarser_neighbor)

        ! receives ipoin indices
        allocate(recv_interface_p_refine_ipoin(num_ipoin_coarser_neighbor),stat=ier)
        if (ier /= 0) stop 'Error allocating array recv_interface_p_refine_ipoin'

        call irecv_i(recv_interface_p_refine_ipoin,num_ipoin_coarser_neighbor, &
                     my_neighbors_ext_mesh(iinterface), &
                     itag,request_recv_vector_ext_mesh(iinterface))

        call wait_req(request_recv_vector_ext_mesh(iinterface))

        ! compute the union of the two sets (this is done inefficiently (n^2))
        coarser_counter_new = coarser_counter
        ipoin_coarser(:) = 0
        ipoin_coarser(1:coarser_counter) = interface_p_refine_ipoin(1:coarser_counter)

        do ipoin_n = 1,num_ipoin_coarser_neighbor
          ipoin = recv_interface_p_refine_ipoin(ipoin_n)
          if (count(interface_p_refine_ipoin(1:coarser_counter) == ipoin) == 0) then
            coarser_counter_new = coarser_counter_new + 1
            ipoin_coarser(coarser_counter_new) = ipoin

            iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
            num_p_level_coarser_to_update(ilevel) = num_p_level_coarser_to_update(ilevel) + 1
            p_level_coarser_to_update(num_p_level_coarser_to_update(ilevel),ilevel) = iglob
            ! print *, "adding:", iglob, "myrank:",myrank,"ilevel",ilevel
          endif
        enddo
        ! frees temporary array
        deallocate(recv_interface_p_refine_ipoin)

        ! do ipoin_n = 1,coarser_counter
        !   ipoin = interface_p_refine_ipoin(ipoin_n)
        !   if (count(ipoin_coarser(1:num_ipoin_coarser_neighbor) == ipoin) == 0) then
        !     coarser_counter_new = coarser_counter_new + 1
        !     ipoin_coarser(coarser_counter_new) = ipoin
        !     iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
        !     ! num_p_level_coarser_to_update(ilevel) = num_p_level_coarser_to_update(ilevel) + 1
        !     ! p_level_coarser_to_update(num_p_level_coarser_to_update(ilevel),ilevel) = iglob
        !     print *, "adding:", iglob, "myrank:",myrank,"ilevel",ilevel
        !   endif
        ! enddo

        call wait_req(request_send_vector_ext_mesh(iinterface))

        ! call qsort(ipoin_coarser,coarser_counter_new,4,cmp_function)
        call Bubble_Sort(ipoin_coarser,coarser_counter_new)

        ! with list of ipoin indexes, build list of iglob indexes that we will use to populate MPI-send array
        ipoin = 0
        do ipoin_n = 1,coarser_counter_new
          ipoin_old = ipoin
          ipoin = ipoin_coarser(ipoin_n)
          if (ipoin <= ipoin_old) then
            print *, "Error: ipoin not sorted..."
          endif
          ! adds point to interface
          iglob = ibool_interfaces_ext_mesh(ipoin,iinterface)
          num_interface_p_refine_ibool(iinterface,ilevel) = num_interface_p_refine_ibool(iinterface,ilevel) + 1
          interface_p_refine_ibool(num_interface_p_refine_ibool(iinterface,ilevel),iinterface,ilevel) = iglob
        enddo

        ! tests MPI interface arrays
        if (TEST_ERROR) then
          ! check correctness of ipoin_coarser arrays
          ! sends ipoin
          call isend_i(ipoin_coarser,coarser_counter_new, &
                       my_neighbors_ext_mesh(iinterface), &
                       itag,request_send_vector_ext_mesh(iinterface))
          ! gets number of points which will be received by neighour
          call get_count_i(my_neighbors_ext_mesh(iinterface),itag,num_ipoin_coarser_neighbor)

          ! receives ipoin
          allocate(recv_interface_p_refine_ipoin(num_ipoin_coarser_neighbor),stat=ier)
          if (ier /= 0) stop 'Error allocating array recv_interace_p_refine_ipoin'

          call irecv_i(recv_interface_p_refine_ipoin,num_ipoin_coarser_neighbor,my_neighbors_ext_mesh(iinterface), &
                       itag,request_recv_vector_ext_mesh(iinterface))
          call wait_req(request_recv_vector_ext_mesh(iinterface))
          call wait_req(request_send_vector_ext_mesh(iinterface))

          ! compares sent/received arrays and makes sure they are identical for this p-level (ilevel)
          do ipoin_n = 1,coarser_counter_new
            if (ipoin_coarser(ipoin_n) /= recv_interface_p_refine_ipoin(ipoin_n)) then
              print *,'Error rank',myrank,': ilevel',ilevel,'interface',iinterface,'point',ipoin_n, &
                     'coarser point',ipoin_coarser(ipoin_n),recv_interface_p_refine_ipoin(ipoin_n)
              call exit_mpi(myrank,"ipoin_coarser don't match!")
            endif
          enddo
          ! frees temporary array
          deallocate(recv_interface_p_refine_ipoin)
        endif ! TEST_ERROR

      endif ! ilevel

      ! user output
      if (myrank == 0 ) then
        if (mod(iinterface,10) == 0) then
          write(IMAIN,*) "    ",iinterface," interfaces ", &
                         " out of ", num_interfaces_ext_mesh,"interfaces"
          if (DEBUG_TIMING) then
            tCPU = wtime() - time_start
            write(IMAIN,*) "    time in seconds = ", tCPU,"s"
          endif
          ! flushes file buffer for main output file (IMAIN)
          call flush_IMAIN()
        endif
      endif

    enddo ! interface

  enddo ! ilevel

  ! frees temporary array
  deallocate(interface_p_refine_ipoin)

  ! checks p_level_coarser_to_update
  do ilevel = 1,num_p_level
    ! gets start index of finest level and end index of current p-level
    is = p_level_iglob_start(1)
    ie = p_level_iglob_end(ilevel)
    if (ilevel < num_p_level) then
      ! considers contributions from coarser to finer p-levels
      do iglob_n = 1,num_p_level_coarser_to_update(ilevel)
        iglob = p_level_coarser_to_update(iglob_n,ilevel)
        ! checks
        if (iglob < ie) call exit_mpi(myrank,"ASSERT: coarser iglob should start in next coarser level")
        if (iglob > NGLOB_AB) call exit_mpi(myrank,"ASSERT: coarser iglob index is out of bounds!")
      enddo
    endif
  enddo

  ! tests MPI interface arrays
  if (TEST_ERROR) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      testing LTS MPI arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      checking iglob index range"
      call flush_IMAIN()
    endif

    ! debug p-level boundary elements
    do i = 0,NPROC-1
      if (myrank == i) then
        do ilevel = 1,num_p_level
          ! outputs boundary elements
          !print *,'myrank',myrank,'ilevel:',ilevel,'num boundary elements:',count(boundary_elem(:,ilevel) .eqv. .true.)
          !print *,'boundary elem:',boundary_elem(:,ilevel)
          ! outputs interface points
          do iinterface = 1,num_interfaces_ext_mesh
            do iglob_n = 1,num_interface_p_refine_ibool(iinterface,ilevel)
              iglob = interface_p_refine_ibool(iglob_n,iinterface,ilevel)
              ! checks global index
              if (iglob < 1 .or. iglob > NGLOB_AB) then
                print *,'Error rank',myrank,': ilevel',ilevel,'interface',iinterface,'point',iglob_n,'iglob:',iglob
                call exit_mpi(myrank,"Error: iglob index of MPI interface is invalid")
              endif
            enddo
          enddo
        enddo
      endif
      call synchronize_all()
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      checking interface locations"
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! checks that locations of interace points are identical between MPI neighbors (in different partitions)
    ! loops over all p-levels
    do ilevel = 1,num_p_level
      ! loops over all interfaces (with other partitions)
      do iinterface = 1, num_interfaces_ext_mesh
        ! number of interface points
        num_points = num_interface_p_refine_ibool(iinterface,ilevel)
        ! checks locations of interface points
        allocate(xyz_send(NDIM,num_points),xyz_recv(NDIM,num_points),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating temporary array xyz_send')
        xyz_send(:,:) = 0._CUSTOM_REAL
        xyz_recv(:,:) = 1.e9_CUSTOM_REAL
        ! fills in locations
        do iglob_n = 1,num_points
          iglob = interface_p_refine_ibool(iglob_n,iinterface,ilevel)
          ! checks index
          if (iglob < 1 .or. iglob > NGLOB_AB) then
            print *,'Error rank',myrank,' : iglob index ',iglob,' out of range ',NGLOB_AB
            call exit_mpi(myrank,'Error p-level MPI interface points range')
          endif
          ! fills send array
          xyz_send(1,iglob_n) = xstore(iglob)
          xyz_send(2,iglob_n) = ystore(iglob)
          xyz_send(3,iglob_n) = zstore(iglob)
        enddo
        ! sends interface point locations
        call isend_cr(xyz_send,NDIM*num_points, &
                      my_neighbors_ext_mesh(iinterface), &
                      itag,request_send_vector_ext_mesh(iinterface))
        ! gets the number of points which will be received from neighbor
        call get_count_i(my_neighbors_ext_mesh(iinterface),itag,num_points_recv)

        ! checks that sent/received number of points are equal
        if (NDIM*num_points /= num_points_recv) then
          print *,'Error rank',myrank,': number of interface points differ!!! ', &
                  'ilevel:',ilevel,'interface:',iinterface,'neighbor:',my_neighbors_ext_mesh(iinterface), &
                  'num_points:',num_points,'num_points_recv:',num_points_recv,'per component',num_points_recv/NDIM
          call exit_MPI(myrank,'Error p-level MPI interface points')
        endif
        ! receives locations
        call irecv_cr(xyz_recv,NDIM*num_points, &
                      my_neighbors_ext_mesh(iinterface), &
                      itag,request_recv_vector_ext_mesh(iinterface))
        call wait_req(request_recv_vector_ext_mesh(iinterface))
        call wait_req(request_send_vector_ext_mesh(iinterface))

        ! checks that locations are identical
        do iglob_n = 1,num_points
          if (maxval(abs(xyz_send(:,iglob_n)-xyz_recv(:,iglob_n))) > 1.e-9) then
            print *,'Error rank',myrank,': level',ilevel,'point:',iglob_n, &
                    'send:',xyz_send(:,iglob_n),'recv:',xyz_recv(:,iglob_n), &
                    'diff:',abs(xyz_send(:,iglob_n)-xyz_recv(:,iglob_n))
            call exit_mpi(myrank,"Error: (x,y,z) of MPI don't match")
          endif
        enddo
        deallocate(xyz_send)
        deallocate(xyz_recv)
      enddo
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "      test result okay"
      if (DEBUG_TIMING) then
        tCPU = wtime() - time_start
        write(IMAIN,*) "      time in seconds = ", tCPU
      endif
      call flush_IMAIN()
    endif
    call synchronize_all()

    time_start = wtime()
  endif ! TEST_ERROR

  ! sets up number of local time steps
  allocate(p_level_m_loops(num_p_level),stat=ier)
  if (ier /= 0) stop 'Error allocating p_level_m_loops array'

  ! example:
  !   num_p_level = 2 and p_level = [ 4  1 ]    -> p_level_m_loops = [ 4  1 ]
  !   num_p_level = 3 and p_level = [ 8  4  1 ] -> p_level_m_loops = [ 2  4  1 ]
  !
  ! note:
  !   p - integer fraction of coarsest time step, such that local time step = dt/p
  !   m - number of local loops before stepping into next level

  ! coarsest level only takes a single local time step
  p_level_m_loops(num_p_level) = 1
  ! finer levels
  do m = 1,(num_p_level-1)
    p_level_m_loops(m) = p_level(m)/p_level(m+1)
  enddo

  ! checks boundary node p-values
  do ilevel = 1,num_p_level
    ! current boundary elements
    current_lts_boundary_elem => boundary_elem(:,ilevel)

    ! p (dt/p) refinement number in this specified p-level
    p = p_level(ilevel)

    ! inner/outer elements
    do iphase = 1,2
      ! choses inner/outer elements
      if (iphase == 1) then
        num_elements = nspec_outer_elastic
      else
        num_elements = nspec_inner_elastic
      endif

      do ispec_p = 1,num_elements
        ! returns element id from stored element list
        ispec = phase_ispec_inner_elastic(ispec_p,iphase)

        ! only elements belonging to p-level boundary
        if (.not. current_lts_boundary_elem(ispec)) cycle

        ! checks element nodes
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool(i,j,k,ispec)
              ! checks if node belongs to this or coarser p-level
              p_node = iglob_p_refine(iglob)
              if (p_node > p) then
                ! coarser p-levels have smaller p values, finer p-levels have larger p values,
                ! such that local time step delta_lts==DT/p
                print *,'Error: boundary node p value is invalid: ',p_node,' on level',ilevel,'p',p
                print *,'       iglob ',iglob,'i/j/k/ispec',i,j,k,ispec
                call exit_mpi(myrank,"Error: invalid p-value node; Assert(boundary nodes are this level or coarser only)")
              endif
            enddo
          enddo
        enddo
      enddo ! ispec
    enddo ! iphase
  enddo ! ilevel

  ! computes true cost of MPI in terms of total elements sent
  mpi_cost = 0
  call lts_mpi_cost(num_p_level,mpi_cost)

  ! gathering all mpi
  allocate(mpi_cost_gather(NPROC),stat=ier)
  if (ier /= 0) stop 'Error allocating mpi_cost_gather'
  mpi_cost_gather(:) = 0

  call gather_all_singlei(mpi_cost,mpi_cost_gather,NPROC)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "    Communication cost: sum(mpi_cost)          = ",sum(mpi_cost_gather)
    write(IMAIN,*) "                        avg(mpi_cost per proc) = ",sngl(sum(mpi_cost_gather)/dble(NPROC))
    write(IMAIN,*) "    p-level boundary setup done"
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! frees temporary arrays
  deallocate(mpi_cost_gather)
  deallocate(ipoin_coarser)

  contains

  ! helper functions

  ! sorting routine
  subroutine Bubble_Sort(a,len)
    integer :: len
    integer, INTENT(inout), DIMENSION(len) :: a

    integer :: temp
    INTEGER :: i, j
    LOGICAL :: swapped = .true.

    DO j = len-1, 1, -1
      swapped = .false.
      DO i = 1, j
        if (a(i) > a(i+1)) then
          temp = a(i)
          a(i) = a(i+1)
          a(i+1) = temp
          swapped = .true.
        endif
      enddo
      if (.not. swapped) EXIT
    enddo
  end subroutine Bubble_Sort

  !
  !--------
  !

  ! comparison function for sorting integers
  integer function cmp_function(a1, a2)
    integer a1, a2
    cmp_function = a1 - a2
  end function cmp_function


  end subroutine lts_setup_level_boundary

!
!--------------------------------------------------------------
!

  recursive subroutine lts_mpi_cost(ilevel,mpi_cost)

  use constants, only: NDIM
  use specfem_par, only: num_interfaces_ext_mesh
  use specfem_par_lts, only: num_interface_p_refine_ibool,p_level_m_loops

  implicit none

  integer,intent(in) :: ilevel
  integer,intent(inout) :: mpi_cost

  ! local parameters
  integer :: m, iinterface

  ! computes true cost of MPI in terms of total elements sent
  do m = 1,p_level_m_loops(ilevel)
    ! recursive through all levels
    if (ilevel > 1) call lts_mpi_cost(ilevel-1, mpi_cost)

    do iinterface = 1, num_interfaces_ext_mesh
      mpi_cost = mpi_cost + NDIM * num_interface_p_refine_ibool(iinterface,ilevel)
    enddo
  enddo

  end subroutine lts_mpi_cost

!
!--------------------------------------------------------------
!

  subroutine lts_prepare_mass_matrices()

  use constants, only: CUSTOM_REAL,NDIM,myrank

  use shared_parameters, only: DT, &
    ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION, &
    STACEY_ABSORBING_CONDITIONS,PML_CONDITIONS

  use specfem_par, only: NGLOB_AB, NSPEC_AB, NSPEC_IRREGULAR, &
    ibool, rhostore, &
    jacobianstore,irregular_element_number,jacobian_regular,wxgll,wygll,wzgll

  ! absorbing boundary
  use specfem_par, only: num_abs_boundary_faces, &
    abs_boundary_ispec, abs_boundary_ijk, abs_boundary_normal, abs_boundary_jacobian2Dw

  use specfem_par_elastic, only: ispec_is_elastic,rho_vp,rho_vs,rmassx,rmassy,rmassz,rmass

  use specfem_par_lts, only: cmassxyz, rmassxyz, rmassxyz_mod

  implicit none
  integer :: ier

  ! safety checks
  if (ACOUSTIC_SIMULATION) &
    stop 'LTS mode for acoustic simulations and mass matrices not implemented yet'
  if (POROELASTIC_SIMULATION) &
    stop 'LTS mode for poroelastic simulations and mass matrices not implemented yet'
  if (PML_CONDITIONS) &
    stop 'LTS mode w/ PML boundaries not implemented yet'

  ! allocates additional LTS mass matrices
  allocate(rmassxyz(NDIM,NGLOB_AB),stat=ier)
  if ( ier /= 0 ) call exit_MPI(myrank,'error allocating temporary rmassxyz array')
  rmassxyz(:,:) = 0.0_CUSTOM_REAL

  allocate(rmassxyz_mod(NDIM,NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'error in allocate arrays rmassxyz_mod,.. ')
  rmassxyz_mod(:,:) = 0.0_CUSTOM_REAL

  allocate(cmassxyz(NDIM,NGLOB_AB), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'error in allocate arrays cmassx,.. ')
  cmassxyz(:,:) = 0.0_CUSTOM_REAL

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! since LTS modifies the time step size DT, we need to re-calculate the Stacey contributions to the mass matrix
    if (STACEY_ABSORBING_CONDITIONS) then
      ! note: in prepare_timerun_mass_matrices() routine, the rmassx,rmassy,rmassz will have the Stacey contributions already added.
      !       here we need the main mass matrix and the Stacey contributions separately, thus will re-create both.
      !       the Stacey contribution would need to be re-created anyway, in case the time step DT has been changed by LTS.
      !
      ! re-creates (main) mass matrix
      call define_mass_matrices_elastic(NGLOB_AB,NSPEC_AB,NSPEC_IRREGULAR,ibool,rhostore, &
                                        jacobianstore,irregular_element_number,jacobian_regular, &
                                        wxgll,wygll,wzgll,ispec_is_elastic, &
                                        rmass)

      ! re-creates Stacey contributions
      if (num_abs_boundary_faces > 0) then
        call define_mass_matrices_Stacey_elastic(NGLOB_AB,NSPEC_AB,DT,ibool,rho_vp,rho_vs, &
                                                 num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_ijk, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 ispec_is_elastic, &
                                                 rmassx, rmassy, rmassz)
      endif

      ! stores initial boundary contributions
      cmassxyz(1,:) = rmassx(:)
      cmassxyz(2,:) = rmassy(:)
      cmassxyz(3,:) = rmassz(:)

      ! modified mass matrices with contributions
      rmassxyz_mod(1,:) = rmass(:) + rmassx(:)
      rmassxyz_mod(2,:) = rmass(:) + rmassy(:)
      rmassxyz_mod(3,:) = rmass(:) + rmassz(:)

      ! re-sets mass matrices without contributions
      rmassx(:) = rmass(:)
      rmassy(:) = rmass(:)
      rmassz(:) = rmass(:)

    else
      ! no absorbing boundary contributions
      ! boundary contributions have zero value
      cmassxyz(:,:) = 0.0_CUSTOM_REAL
      ! default mass matrices
      rmassxyz_mod(1,:) = rmassx(:)
      rmassxyz_mod(2,:) = rmassy(:)
      rmassxyz_mod(3,:) = rmassz(:)
    endif
  endif

  end subroutine lts_prepare_mass_matrices

!
!--------------------------------------------------------------
!

  subroutine lts_prepare_mass_matrices_invert()

  use constants, only: CUSTOM_REAL,IMAIN,myrank

  use shared_parameters, only: NPROC, &
    ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION,STACEY_ABSORBING_CONDITIONS

  ! MPI
  use specfem_par, only: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,my_neighbors_ext_mesh

  use specfem_par, only: NGLOB_AB
  use specfem_par_elastic, only: rmassx,rmassy,rmassz

  use specfem_par_lts, only: cmassxyz, rmassxyz, rmassxyz_mod

  implicit none

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tmp_rmass
  real(kind=CUSTOM_REAL) :: rmass_min,rmass_max,rmass_min_glob,rmass_max_glob
  integer,dimension(1) :: rank_min,rank_max
  integer :: idim,ier

  ! safety checks
  if (ACOUSTIC_SIMULATION) &
    stop 'LTS mode for acoustic simulations and mass matrices not implemented yet'
  if (POROELASTIC_SIMULATION) &
    stop 'LTS mode for poroelastic simulations and mass matrices not implemented yet'

  ! elastic domains
  if (ELASTIC_SIMULATION) then
    ! assembles w/ MPI LTS mass matrices

    ! temporary array
    allocate(tmp_rmass(NGLOB_AB),stat=ier)
    if ( ier /= 0 ) call exit_MPI(myrank,'error allocating temporary tmp_rmass array')
    tmp_rmass(:) = 0.0_CUSTOM_REAL

    ! assembles "modified" mass matrix
    tmp_rmass(:) = rmassxyz_mod(1,:)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                       my_neighbors_ext_mesh)
    rmassxyz_mod(1,:) = tmp_rmass(:)

    tmp_rmass(:) = rmassxyz_mod(2,:)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                       my_neighbors_ext_mesh)
    rmassxyz_mod(2,:) = tmp_rmass(:)

    tmp_rmass(:) = rmassxyz_mod(3,:)
    call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                       my_neighbors_ext_mesh)
    rmassxyz_mod(3,:) = tmp_rmass(:)

    ! assembles absorbing boundary contributions
    if (STACEY_ABSORBING_CONDITIONS) then
      ! note: we will need cmassxyz in assembled form, but not inverted
      tmp_rmass(:) = cmassxyz(1,:)
      call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                         my_neighbors_ext_mesh)
      cmassxyz(1,:) = tmp_rmass(:)

      tmp_rmass(:) = cmassxyz(2,:)
      call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                         my_neighbors_ext_mesh)
      cmassxyz(2,:) = tmp_rmass(:)

      tmp_rmass(:) = cmassxyz(3,:)
      call assemble_MPI_scalar_blocking(NPROC,NGLOB_AB,tmp_rmass, &
                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                         nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                         my_neighbors_ext_mesh)
      cmassxyz(3,:) = tmp_rmass(:)
    endif

    ! frees temporary array
    deallocate(tmp_rmass)

    ! inverts modified mass matrices
    where(rmassxyz_mod <= 0.0_CUSTOM_REAL) rmassxyz_mod = 1.0_CUSTOM_REAL
    rmassxyz_mod(:,:) = 1.0_CUSTOM_REAL / rmassxyz_mod(:,:)

    ! stores final, inverted mass matrices of each component into single array
    ! note: when calling this routine, rmassx,rmassy,rmassz have been assembled and inverted.
    !       we will needs these inverted mass matrices and store them into rmassxyz
    rmassxyz(1,:) = rmassx(:)
    rmassxyz(2,:) = rmassy(:)
    rmassxyz(3,:) = rmassz(:)

    ! statistics
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  LTS mass matrix statistics: '
    endif

    ! mass matrix rmassxyz
    do idim = 1,3
      rmass_min = minval( rmassxyz(idim,:))
      rmass_max = maxval( rmassxyz(idim,:))
      rank_min = minloc( rmassxyz(idim,:))
      rank_max = maxloc( rmassxyz(idim,:))
      call max_all_cr(rmass_max,rmass_max_glob)
      call min_all_cr(rmass_min,rmass_min_glob)
      ! debug
      !print *,'rmassxyz min',rmass_min,rank_min,'max',rmass_max,rank_max,'rank',myrank
      ! user output
      if ( myrank == 0 ) then
        select case(idim)
        case (1)
          write(IMAIN,*) '  mass matrix rmassxyz    : x min/max   = ',rmass_min_glob,rmass_max_glob
        case (2)
          write(IMAIN,*) '  mass matrix rmassxyz    : y min/max   = ',rmass_min_glob,rmass_max_glob
        case (3)
          write(IMAIN,*) '  mass matrix rmassxyz    : z min/max   = ',rmass_min_glob,rmass_max_glob
        end select
        call flush_IMAIN()
      endif
    enddo
    if (myrank == 0) write(IMAIN,*)

    ! mass matrix rmassxyz_mod
    do idim = 1,3
      rmass_min = minval( rmassxyz_mod(idim,:))
      rmass_max = maxval( rmassxyz_mod(idim,:))
      rank_min = minloc( rmassxyz_mod(idim,:))
      rank_max = maxloc( rmassxyz_mod(idim,:))
      call max_all_cr(rmass_max,rmass_max_glob)
      call min_all_cr(rmass_min,rmass_min_glob)
      ! debug
      !print *,'rmassxyz_mod min',rmass_min,rank_min,'max',rmass_max,rank_max,'rank',myrank
      ! user output
      if (myrank == 0) then
        select case(idim)
        case (1)
          write(IMAIN,*) '  mass matrix rmassxyz_mod: x min/max   = ',rmass_min_glob,rmass_max_glob
        case (2)
          write(IMAIN,*) '  mass matrix rmassxyz_mod: y min/max   = ',rmass_min_glob,rmass_max_glob
        case (3)
          write(IMAIN,*) '  mass matrix rmassxyz_mod: z min/max   = ',rmass_min_glob,rmass_max_glob
        end select
        call flush_IMAIN()
      endif
    enddo
  endif

  end subroutine lts_prepare_mass_matrices_invert

!
!--------------------------------------------------------------
!

  subroutine lts_prepare_gpu()


!#TODO: LTS gpu setup
!       note that this routine will need to be called after prepare_constants_device()
!       as the Mesh_pointer will get created there, which is needed here.

  stop 'LTS setup on GPU not implemented yet'

!  use specfem_par, only: Mesh_pointer
!
!  use specfem_par_lts, only: &
!    num_interface_p_refine_boundary, interface_p_refine_boundary, max_nibool_interfaces_boundary, &
!    p_elem
!
!  implicit none
!
!  integer, dimension(:,:,:), allocatable :: element_list
!  integer, dimension(:,:), allocatable :: num_element_list
!  integer :: iispec
!
!  ! timing
!  logical, parameter :: DEBUG_TIMING = .false.
!  double precision, external :: wtime
!  double precision :: time_start,tCPU
!
!  ! sets up GPU transfers of LTS arrays
!
!  ! timing
!  if (DEBUG_TIMING) time_start = wtime()
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) "  loading LTS arrays"
!    call flush_IMAIN()
!  endif
!  call synchronize_all()
!
!  ! transfers LTS mass matrices
!  call prepare_mass_boundary_fields_lts_device(Mesh_pointer, rmassxyz, rmassxyz_mod, cmassxyz)

!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) "  transfering LTS p-level boundary arrays to GPU"
!    call flush_IMAIN()
!  endif
!  call synchronize_all()
!
!  ! allocates temporary arrays
!  allocate(element_list(NSPEC_AB,2,num_p_level), &
!           num_element_list(2,num_p_level),stat=ier)
!  if (ier /= 0) stop 'Error allocating arrays element_list,...'
!  element_list(:,:,:) = 0
!  num_element_list(:,:) = 0
!
!  ! builds list of elements belonging to p-level, without boundary elements
!  do ilevel = 1,num_p_level
!    do iphase = 1,2
!      if (iphase == 1) then
!        num_elements = nspec_outer_elastic
!      else
!        num_elements = nspec_inner_elastic
!      endif
!      ! loops over all inner/outer elements
!      iispec = 1
!      num_element_list(iphase,ilevel) = 0
!      do ispec=1,num_elements
!        ispec_p = phase_ispec_inner_elastic(ispec,iphase)
!
!        ! skips elements in other p-levels
!        if (.not. p_elem(ispec_p,ilevel)) cycle
!        ! skips elements belonging to p-level boundary
!        if (boundary_elem(ispec_p,ilevel) .eqv. .true.) cycle
!
!        ! adds element
!        element_list(iispec,iphase,ilevel) = ispec_p
!        num_element_list(iphase,ilevel) = iispec
!        iispec = iispec + 1
!      enddo
!    enddo ! iphase
!  enddo ! ilevel
!
!  ! sends p-element list to GPU
!  call transfer_element_list_to_device(Mesh_pointer,element_list,num_element_list)
!
!  ! transfer ibool and ilevel from CPU to GPU
!  call transfer_boundary_element_list_to_device(Mesh_pointer,p_level_boundary_node, &
!                                                p_level_boundary_ilevel_from,p_level_boundary_ispec)
!
!  ! sends list for coarser nodes to GPU
!  call setup_r_boundaries_time_stepping(Mesh_pointer,p_level_coarser_to_update)
!
!  ! sets up interface points
!  allocate(num_interface_p_refine_boundary(num_p_level),stat=ier)
!  if (ier /= 0) stop 'Error allocating arrays num_interface_p_refine_boundary'
!  num_interface_p_refine_boundary(:) = 0
!
!  ! counts all nodes on MPI interfaces between levels
!  max_nibool_interfaces_boundary = 0
!  do ilevel = 1,num_p_level-1
!    num_interface_p_refine_boundary(ilevel) = 0
!    do iinterface = 1, num_interfaces_ext_mesh
!      ! loops over all points on this interface
!      do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!        ! increases points counter
!        num_interface_p_refine_boundary(ilevel) = num_interface_p_refine_boundary(ilevel) + 1
!      enddo
!    enddo
!
!    ! gets maximum number of points on a refine interface
!    if (num_interface_p_refine_boundary(ilevel) > max_nibool_interfaces_boundary) then
!      max_nibool_interfaces_boundary = num_interface_p_refine_boundary(ilevel)
!    endif
!  enddo
!  ! checks
!  if (num_p_level > 1 .and. NPROC > 1 .and. max_nibool_interfaces_boundary == 0) &
!    stop 'Error LTS refinement without refinement nodes'
!
!  ! debug output
!  !do i = 0,NPROC
!  !  if (myrank == i) then
!  !    print *,'rank',myrank
!  !    do ilevel = 1,num_p_level-1
!  !      print *,'  level',ilevel,' num_interface_p_refine_boundary = ',num_interface_p_refine_boundary(ilevel)
!  !    enddo
!  !    print *,'  max_nibool_interfaces_ext_mesh = ',max_nibool_interfaces_ext_mesh,max_nibool_interfaces_boundary
!  !  endif
!  !  call synchronize_all()
!  !enddo
!
!  ! allocates refinement interfaces
!  allocate(interface_p_refine_boundary(max_nibool_interfaces_boundary,num_p_level),stat=ier)
!  if (ier /= 0) stop 'Error allocating arrays interface_p_refine_boundary'
!
!  ! Precompute list of iglob nodes for each level that will be sent over MPI
!  num_interface_p_refine_boundary(:) = 0
!  interface_p_refine_boundary(:,:) = 0
!  do ilevel = 1,num_p_level-1
!    do iinterface = 1, num_interfaces_ext_mesh
!      ! loops over all points on this interface
!      do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!        ! gets global index
!        iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
!
!        ! increases points counter
!        num_interface_p_refine_boundary(ilevel) = num_interface_p_refine_boundary(ilevel) + 1
!
!        ! checks
!        if (num_interface_p_refine_boundary(ilevel) > max_nibool_interfaces_boundary) then
!          print *,'Error rank: ',myrank,'level ',ilevel,'out of ',num_p_level, &
!                  ' num_interface_p_refine_boundary =',num_interface_p_refine_boundary(ilevel), &
!                  ' exceeds ',max_nibool_interfaces_ext_mesh,max_nibool_interfaces_boundary
!          stop 'Error num_interface_p_refine_boundary exceeds bounds'
!        endif
!
!        ! adds global index to interface
!        interface_p_refine_boundary(num_interface_p_refine_boundary(ilevel),ilevel) = iglob
!      enddo
!    enddo
!  enddo
!
!  ! sends interface lists to GPU
!  call setup_mpi_boundaries_lts(Mesh_pointer, &
!                                num_interface_p_refine_ibool, &
!                                interface_p_refine_ibool, &
!                                num_interface_p_refine_boundary, &
!                                interface_p_refine_boundary, &
!                                num_interfaces_ext_mesh, &
!                                max_nibool_interfaces_ext_mesh, &
!                                max_nibool_interfaces_boundary, &
!                                num_p_level)
!
!  ! frees temporary arrays
!  deallocate(element_list)
!  deallocate(num_element_list)
!
!  ! user output
!  call synchronize_all()
!  if (DEBUG_TIMING) then
!    ! timing
!    if (myrank == 0 ) then
!      tCPU = wtime() - time_start
!      write(IMAIN,*) "   time in seconds = ", tCPU
!      call flush_IMAIN()
!    endif
!  endif

  end subroutine lts_prepare_gpu
