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

! HDF5 I/O server
! (works with HDF5 file outputs)

! module for storing info. concering io

module io_server_hdf5

  use constants, only: CUSTOM_REAL

  use shared_parameters, only: HDF5_IO_NODES, &
    IO_storage_task, IO_compute_task

  implicit none

#ifdef USE_HDF5
! only with HDF5 support

  ! public routines
  public :: initialize_io_server
  public :: finalize_io_server

  public :: do_io_start_idle
  public :: pass_info_to_io
  public :: wait_all_send

  ! public parameters
  public :: NtoM
  public :: if_collect_server

  public :: nglob_this_io, nelm_this_io
  public :: nglob_all_server, nspec_all_server

  ! seismogram arrays
  public :: seismo_press
  public :: seismo_displ, seismo_veloc, seismo_accel

  ! movie arrays
  public :: surf_x,   surf_y,   surf_z, &
            surf_ux,  surf_uy,  surf_uz
  public :: surf_x_aug,   surf_y_aug,   surf_z_aug, &
            surf_ux_aug,  surf_uy_aug,  surf_uz_aug
  public :: shake_ux, shake_uy, shake_uz, &
            shake_ux_aug, shake_uy_aug, shake_uz_aug

  public :: vd_pres, vd_divglob, vd_div, &
            vd_curlx, vd_curly, vd_curlz, &
            vd_velox, vd_veloy, vd_veloz
  public :: vd_stressxx, vd_stressyy, vd_stresszz, &
            vd_stressxy, vd_stressxz, vd_stressyz

  public :: nglob_offset
  public :: id_proc_loc2glob, id_proc_glob2loc
  public :: dest_ioids, dest_ionod

  ! MPI message tags
  public :: io_tag_seismo_body_disp
  public :: io_tag_seismo_body_velo
  public :: io_tag_seismo_body_acce
  public :: io_tag_seismo_body_pres

  public :: io_tag_shake_ux, io_tag_shake_uy, io_tag_shake_uz

  public :: io_tag_surface_ux, io_tag_surface_uy, io_tag_surface_uz
  public :: io_tag_surface_x, io_tag_surface_y, io_tag_surface_z

  public :: io_tag_vol_curlx, io_tag_vol_curly, io_tag_vol_curlz
  public :: io_tag_vol_div, io_tag_vol_divglob
  public :: io_tag_vol_pres
  public :: io_tag_vol_velox, io_tag_vol_veloy, io_tag_vol_veloz
  public :: io_tag_vol_stressxx, io_tag_vol_stressyy, io_tag_vol_stresszz, &
            io_tag_vol_stressxy, io_tag_vol_stressxz, io_tag_vol_stressyz

  public :: io_tag_surface_coord_len
  public :: io_tag_surface_nfaces
  public :: io_tag_surface_offset

  public :: io_tag_vol_elmconn
  public :: io_tag_vol_ioid
  public :: io_tag_vol_nglob
  public :: io_tag_vol_nmsg
  public :: io_tag_vol_nodex, io_tag_vol_nodey, io_tag_vol_nodez
  public :: io_tag_vol_nspec
  public :: io_tag_vol_sendlist

  ! MPI requests
  public :: n_req_surf, n_req_vol
  public :: req_dump_surf, req_dump_vol

  public :: nproc_io
  public :: n_msg_vol_each_proc

  ! verbosity
  public :: VERBOSE

  private

  ! MPI tags for io server implementation
  !integer :: io_tag_seismo_nrec       = 10  ! unused
  integer :: io_tag_seismo_ids_rec    = 11
  integer :: io_tag_seismo_body_disp  = 12
  integer :: io_tag_seismo_body_velo  = 13
  integer :: io_tag_seismo_body_acce  = 14
  integer :: io_tag_seismo_body_pres  = 15
  integer :: io_tag_seismo_tzero      = 16
  integer :: io_tag_dt                = 17
  integer :: io_tag_nstep             = 18
  integer :: io_tag_seismo_length     = 19

  integer :: io_tag_surface_nfaces    = 20
  integer :: io_tag_surface_offset    = 21
  integer :: io_tag_surface_x         = 22
  integer :: io_tag_surface_y         = 23
  integer :: io_tag_surface_z         = 24
  integer :: io_tag_surface_ux        = 25
  integer :: io_tag_surface_uy        = 26
  integer :: io_tag_surface_uz        = 27
  integer :: io_tag_surface_coord_len = 28

  integer :: io_tag_shake_ux          = 30
  integer :: io_tag_shake_uy          = 31
  integer :: io_tag_shake_uz          = 32

  integer :: io_tag_vol_pres          = 40
  integer :: io_tag_vol_divglob       = 41
  integer :: io_tag_vol_div           = 42
  integer :: io_tag_vol_curlx         = 43
  integer :: io_tag_vol_curly         = 44
  integer :: io_tag_vol_curlz         = 45
  integer :: io_tag_vol_velox         = 46
  integer :: io_tag_vol_veloy         = 47
  integer :: io_tag_vol_veloz         = 48
  integer :: io_tag_vol_nmsg          = 49
  integer :: io_tag_vol_nspec         = 50
  integer :: io_tag_vol_nglob         = 51
  integer :: io_tag_vol_sendlist      = 52
  integer :: io_tag_vol_ioid          = 53
  integer :: io_tag_vol_elmconn       = 54
  integer :: io_tag_vol_nodex         = 55
  integer :: io_tag_vol_nodey         = 56
  integer :: io_tag_vol_nodez         = 57
  !integer :: io_tag_end               = 60   ! unused

  integer :: io_tag_num_recv          = 61
  integer :: io_tag_local_rec         = 62

  integer :: io_tag_vol_stressxx      = 63
  integer :: io_tag_vol_stressyy      = 64
  integer :: io_tag_vol_stresszz      = 65
  integer :: io_tag_vol_stressxy      = 66
  integer :: io_tag_vol_stressxz      = 67
  integer :: io_tag_vol_stressyz      = 68

  ! mpi_req dump
  integer :: n_req_surf = 0
  integer :: n_req_vol = 0
  integer, dimension(9) :: req_dump_vol
  integer, dimension(3) :: req_dump_surf

  ! responsible id of io node
  integer :: dest_ionod = 0

  ! movie arrays
  type vol_data_dump
    integer, dimension(:), allocatable :: req                       !=VAL_NOT_ASSIGNED should be
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: d1darr
  end type vol_data_dump

  ! dumps for volume data
  type(vol_data_dump) :: vd_pres, vd_divglob, vd_div, &
                         vd_curlx, vd_curly, vd_curlz, &
                         vd_velox, vd_veloy, vd_veloz

  type(vol_data_dump) :: vd_stressxx, vd_stressyy, vd_stresszz, &
                         vd_stressxy, vd_stressxz, vd_stressyz

  logical :: NtoM = .true. ! if this option set .false. SPECFEM3D takes N-to-1 I/O strategy,
                           ! i.e., volume data will be out to one single .h5 file.
                           ! This gives the output files great simplicity, however I/O performance will be
                           ! degraded for collective h5 I/O.

  logical :: if_collect_server = .false. ! need to be set .false. for NtoM = .true.

  ! number of computer nodes sending info to each io node
  integer :: nproc_io

  ! id for io node
  integer :: my_io_id

  ! local-global processor id relation
  integer, dimension(:), allocatable :: id_proc_glob2loc, id_proc_loc2glob

  ! io node id <-> proc id in compute nodes relation
  integer, dimension(:), allocatable :: dest_ioids

  integer :: nglob_all_server, nspec_all_server ! total number of elements and nodes for all procs
  integer :: nglob_this_io, nelm_this_io        ! number of nodes and elements in this IO group

  integer :: n_seismo_type = 0
  integer :: n_procs_with_rec = 0

  integer :: n_msg_seismo_each_proc = 1
  integer :: n_msg_surf_each_proc = 3
  integer :: n_msg_shake_each_proc = 3
  integer :: n_msg_vol_each_proc = 0

  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable   :: seismo_press
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable   :: seismo_displ, seismo_veloc, seismo_accel
  integer, dimension(:,:), allocatable                    :: id_rec_globs

  integer, parameter :: VAL_NOT_ASSIGNED = 9999999

  ! movie arrays for hdf5 output
  integer, dimension(:), allocatable :: nglob_offset
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_x,   surf_y,   surf_z, &
                                                       surf_ux,  surf_uy,  surf_uz
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_x_aug,   surf_y_aug,   surf_z_aug, &
                                                       surf_ux_aug,  surf_uy_aug,  surf_uz_aug
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: shake_ux, shake_uy, shake_uz, &
                                                       shake_ux_aug, shake_uy_aug, shake_uz_aug

  ! verbose output (for debugging)
  logical, parameter :: VERBOSE = .false.

#endif

contains


!-------------------------------------------------------------------------------------------------
!
! i/o server setup
!
!-------------------------------------------------------------------------------------------------

  subroutine initialize_io_server()

! initialization of IO server splits compute and io nodes
#ifdef USE_HDF5

  use constants, only: IMAIN,myrank,MAX_STRING_LEN
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,NPROC

  implicit none

  integer :: sizeval
  integer :: mpi_comm,split_comm,inter_comm
  integer :: key,io_start,comp_start
  ! test node name
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: node_names
  character(len=MAX_STRING_LEN) :: this_node_name
  integer :: node_len

  ! split comm into computation nodes and io node
  ! here we use the last HDF5_IO_NODES ranks (intra comm) as the io node
  ! for using one additional node, xspecfem3D need to be run with + HDF5_IO_NODES node
  ! thus for running mpirun, it should be like e.g.
  ! mpirun -n $((NPROC+HDF5_IO_NODES)) ./bin/xspecfem3D

  ! safety check
  if (HDF5_IO_NODES < 0) stop 'Invalid HDF5_IO_NODES, must be zero or positive'

  ! initializes
  call world_get_comm(mpi_comm)         ! mpi_comm == my_local_mpi_comm_world

  ! default MPI group, no i/o server
  call world_set_comm_inter(mpi_comm)   ! my_local_mpi_comm_inter = mpi_comm

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! select the task of this proc
  ! get the local mpi_size and rank
  call world_size(sizeval)
  call world_rank(myrank)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'HDF5 I/O server run:'
    write(IMAIN,*) '  number of dedicated HDF5 IO nodes = ',HDF5_IO_NODES
    write(IMAIN,*) '  total number of MPI processes     = ',sizeval
    write(IMAIN,*)
    write(IMAIN,*) '  separating subgroups for compute tasks with ',sizeval - HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*) '                       for io tasks      with ',HDF5_IO_NODES,'MPI processes'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks if solver run with a number of MPI processes equal to (NPROC + HDF5_IO_NODES)
  if (sizeval /= NPROC + HDF5_IO_NODES) then
    if (myrank == 0) then
      print *,"Error: HDF5 IO server needs ",HDF5_IO_NODES," additional process together with NPROC",NPROC,"processes"
      print *
      print *,"However, this simulation runs with ",sizeval,"MPI processes,"
      print *,"where instead it should run with ",NPROC + HDF5_IO_NODES,"MPI processes"
      print *
      print *,"Please check your call to 'mpirun -np .. xspecfem3D' and rerun with the correct number of processes"
      print *
    endif
    stop 'Invalid number of MPI processes for HDF5_IO_NODES > 0'
  endif

  ! to select io node and compute nodes on the same cluster node.
  allocate(node_names(0:sizeval-1))

  call world_get_processor_name(this_node_name,node_len)

  ! share the node name between procs
  call gather_all_all_single_ch(this_node_name,node_names,sizeval,MAX_STRING_LEN)

  call synchronize_all()

  ! select the task of this proc
  call select_io_node(node_names,myrank,sizeval,key,io_start)

  deallocate(node_names)

  ! split communicator into compute_comm and io_comm
  call world_comm_split(mpi_comm, key, myrank, split_comm)

  ! create inter communicator and set as my_local_mpi_comm_inter
  if (IO_storage_task) then
    comp_start = 0
    call world_create_intercomm(split_comm, 0, mpi_comm, comp_start, 1111, inter_comm)
  else
    !io_start = sizeval - HDF5_IO_NODES !+dest_ionod
    call world_create_intercomm(split_comm, 0, mpi_comm, io_start, 1111, inter_comm)
  endif

  ! sets new comm_world subgroup
  call world_set_comm(split_comm)            ! such that: my_local_mpi_comm_world == split_comm

  ! use inter_comm as my_local_mpi_comm_world for all send/recv
  call world_set_comm_inter(inter_comm)      ! such that: my_local_mpi_comm_inter == inter_comm

  ! exclude io nodes from the other compute nodes
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) NPROC = NPROC - HDF5_IO_NODES

  ! re-define myrank (within new comm_world subgroup)
  call world_rank(myrank)
  call world_size(sizeval)

  ! note: IO compute tasks cover the lower part of MPI processes, i.e., with ranks in
  !       the initial 0 to (sizeval-HDF5_IO_NODES)-1 range.
  !       the io nodes with IO_storage_task set to .true. are added at the end, with ranks in
  !       the upper range (sizeval-HDF5_IO_NODES) to (sizeval-1).
  !
  !       for user output, we only have the initial rank 0 as the main process opening the output_solver.txt file.
  !       after separating the tasks and creating new subgroup, we have a rank 0 process in the compute subgroup
  !       as well as a rank 0 process in the storage subgroup.
  !
  !       from here on, we only want the process with rank 0 in the compute group write output to IMAIN.
  !
  ! user output
  if (IO_compute_task) then
    if (myrank == 0) then
      write(IMAIN,*) '  new number of processes in compute group = ',sizeval
      write(IMAIN,*) '  new addtional processes in io group      = ',HDF5_IO_NODES
      write(IMAIN,*) '  MPI i/o server setup done'
      call flush_IMAIN()
    endif
  endif

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine initialize_io_server() called without HDF5 Support."
  print *, "       HDF5_IO_NODES > 0 requires HDF5 support and HDF5_ENABLED must be set to .true."
  print *
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'initialize_io_server() called without HDF5 compilation support'

#endif

  end subroutine initialize_io_server

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_io_server()

#ifdef USE_HDF5

  implicit none

  ! checks if anything to do
  if (HDF5_IO_NODES == 0) return

  ! finish MPI subgroup
  if (HDF5_IO_NODES > 0) then
    ! wait for all to finish
    call synchronize_inter()
    ! free subgroup
    call world_comm_free_inter()
  endif
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine finalize_io_server() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'finalize_io_server() called without HDF5 compilation support'

#endif

  end subroutine finalize_io_server

!
!-------------------------------------------------------------------------------------------------
!
!
#ifdef USE_HDF5

  subroutine select_io_node(node_names, myrank, sizeval, key, io_start)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer, intent(in)                         :: myrank, sizeval
  integer, intent(out)                        :: key, io_start
  character(len=MAX_STRING_LEN), dimension(0:sizeval-1), intent(in) :: node_names
  ! local parameters
  integer, dimension(sizeval) :: n_procs_on_node ! number of procs on each cluster node
  integer, dimension(:), allocatable :: n_ionode_on_cluster ! number of ionode on the cluster nodes
  integer :: i,j,c,n_cluster_node,my_cluster_id,n_rest_io,n_ionode,n_comp_node
  real(kind=CUSTOM_REAL) :: io_ratio ! dum
  character(len=MAX_STRING_LEN), dimension(sizeval) :: dump_node_names ! names of cluster nodes

  ! initialize
  my_cluster_id = -1
  n_cluster_node = 0
  n_procs_on_node(:) = 0
  dump_node_names(:) = "nan"

  ! BUG: nprocs goes wrong when 2 cluster nodes and 64 compute 8 io procs
  ! (only 62 compute nodes are assigned)

  ! cluster_node_nums = [n_1,...,n_i,...n_cn,-1,-1,...] ! n_i is the number of procs on
  ! each cluster node
  do i = 1, sizeval
    ! search the node name already found and registered in the dump_node_names array
    c = 0
    do j = 1, sizeval
      if (node_names(i-1) == dump_node_names(j)) c = j
    enddo

    ! if node_name[i-1] is not registered yet
    if (c == 0) then
      ! count the number of cluster node
      n_cluster_node = n_cluster_node + 1

      dump_node_names(n_cluster_node) = node_names(i-1) ! register the name
      n_procs_on_node(n_cluster_node) = 1 ! count up the number of procs on this cluster node

      c = n_cluster_node
    else ! if node_name[i] has already be found, count up the number of procs
      n_procs_on_node(c) = n_procs_on_node(c) + 1
    endif

    ! the id of cluster which this process belongs to
    if (i-1 == myrank) then
      my_cluster_id = c
    endif
  enddo

  ! warning when a cluster node has only one single proc.
  do i = 1, n_cluster_node
    if (n_procs_on_node(i) == 1) then
      print *
      print *, "***************************************************************************"
      print *, "IO server Warning:"
      print *, "  node name: " // trim(dump_node_names(i)) // " has only one procs."
      print *, "  this may lead a io performance issue by inter-clusternode communication."
      print *, "***************************************************************************"
      print *
    endif
  enddo

  ! select HDF5_IO_NODES of io nodes
  allocate(n_ionode_on_cluster(n_cluster_node))
  !! decide the number of io node on each cluster node
  ! at least one io node on each cluster node
  n_ionode_on_cluster(:) = 1

  ! check if the total number of io node > HDF5_IO_NODES
  if (sum(n_ionode_on_cluster) > HDF5_IO_NODES) then
    print *, "Error: HDF5_IO_NODES in Par_file is too small,"
    print *, "       at least one io node for each cluster node is necessary"
    stop 'Invalid HDF5_IO_NODES value too small'
  endif

  ! share the rest of ionodes based on the ratio of compute nodes
  n_rest_io = HDF5_IO_NODES - sum(n_ionode_on_cluster)
  if (n_rest_io /= 0) then
    ! add io nodes one by one to the cluster node where
    ! the io_node/total_node ratio is rowest
    do i = 1, n_rest_io
      io_ratio = 9999.0 !initialize at each i
      do j = 1, n_cluster_node
        if (io_ratio > real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))) then
          ! dump if largest
          io_ratio = real(n_ionode_on_cluster(j))/real(n_procs_on_node(j))
          ! destination of additional io node
          c = j
        endif
      enddo
      ! add one additional io node
      n_ionode_on_cluster(c) = n_ionode_on_cluster(c) + 1
    enddo
  endif

  !! choose the io node from the last rank of each cluster node
  n_ionode = 0
  do i = 1, n_cluster_node
    c = 0
    ! number of compute node on this cluster node
    n_comp_node = n_procs_on_node(i) - n_ionode_on_cluster(i)
    do j = 1, sizeval
      ! if this rank is on i cluster node
      if (dump_node_names(i) == node_names(j-1)) then
        c = c + 1
        if (n_comp_node < c) then
          ! j is io node
          if (j-1 == myrank) then
            ! check  if j is myrank
            IO_storage_task = .true. ! set io node flag
            IO_compute_task = .false.
            key          = 0
            my_io_id     = n_ionode ! id of io_node
            nproc_io     = n_comp_node/n_ionode_on_cluster(i) ! number of compute nodes which use this io node
            dest_ionod   = -1
            if (c-n_comp_node-1 < mod(n_comp_node,n_ionode_on_cluster(i))) nproc_io = nproc_io + 1
          endif

          ! rank of io_start
          if (n_ionode == 0) io_start = j-1

          n_ionode = n_ionode+1

        else
          ! j is compute node
          if (j-1 == myrank) then
            IO_storage_task = .false.
            IO_compute_task = .true.
            key          = 1
            ! set the destination of MPI communication
            dest_ionod   = mod(c-1,n_ionode_on_cluster(i)) + n_ionode
          endif
        endif
      endif
    enddo
  enddo

  ! debug
  if (VERBOSE) then
    if (myrank == 0) then
      print *, "io_server: n_procs_on_node", n_procs_on_node(:)
      print *, "io_server: n_ionode_on_cluster", n_ionode_on_cluster
      print *
    endif
    call flush_stdout()
    call synchronize_all()

    if (IO_storage_task) then
      print *, "io_server: io_task rank ", myrank, " nprocio ",nproc_io
      print *, "io_server: io_task rank ", myrank, " node ", trim(node_names(myrank)), " my_io_id", my_io_id
      print *
    endif
    call flush_stdout
    call synchronize_all()

    do i = 0, n_ionode-1
      if (.not. IO_storage_task) then
        if (dest_ionod == i) then
          print *, "io_server: rank ", myrank, " node ", trim(node_names(myrank)), " dest ", dest_ionod
        endif
      endif
      call flush_stdout()
      call synchronize_all()
    enddo
    call flush_stdout()
    call synchronize_all()

  endif

  deallocate(n_ionode_on_cluster)

  end subroutine select_io_node

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine do_io_start_idle()

#ifdef USE_HDF5

  use constants, only: myrank,my_status_size,my_status_source,my_status_tag

  use shared_parameters, only: NPROC,DT,NSTEP,NTSTEP_BETWEEN_FRAMES, &
    MOVIE_SURFACE,MOVIE_VOLUME,USE_HIGHRES_FOR_MOVIES,CREATE_SHAKEMAP, &
    MOVIE_VOLUME_STRESS, &
    NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS

  use specfem_par, only: nlength_seismogram

  implicit none

  integer :: ier

  ! vars seismo
  integer,dimension(0:NPROC-1) :: islice_num_rec_local
  integer                      :: status(my_status_size)
  integer                      :: rec_count_seismo, n_recv_msg_seismo, max_seismo_out
  integer                      :: it_offset, seismo_out_count
  !integer, dimension(1)        :: nrec_temp

  ! vars surface movie
  integer                       :: rec_count_surf, n_recv_msg_surf, surf_out_count, max_surf_out
  integer                       :: it_io
  integer, dimension(0:NPROC-1) :: nfaces_perproc, surface_offset
  integer                       :: num_nodes

  ! vars shakemap
  integer :: rec_count_shake, n_recv_msg_shake, shake_out_count, max_shake_out

  ! vars volumne movie
  integer                       :: rec_count_vol, n_recv_msg_vol, vol_out_count, max_vol_out
  ! storing the number of elements and GLL nodes
  integer, dimension(0:NPROC-1)           :: nelm_par_proc, nglob_par_proc, nglob_par_proc_offset
  integer, dimension(0:HDF5_IO_NODES-1)  :: nglob_par_io_offset, nelm_par_io_offset
  logical, dimension(6)                   :: val_type_mov ! true if movie file will be created,
                                                          !  (pressure, div_glob, div, curlxyz, velocity_xyz, stress)
  integer, dimension(1)                   :: tmp_1d_iarr
  double precision, dimension(1)          :: tmp_1d_darr
  integer :: tag, tag_src

  ! initialize counters
  rec_count_seismo  = 0
  n_recv_msg_seismo = 0
  it_offset         = 0
  seismo_out_count  = 0
  max_seismo_out    = 0

  rec_count_surf    = 0
  n_recv_msg_surf   = 0
  surf_out_count    = 0
  max_surf_out      = 0

  rec_count_shake   = 0
  n_recv_msg_shake  = 0
  shake_out_count   = 0
  max_shake_out     = 0

  rec_count_vol     = 0
  n_recv_msg_vol    = 0
  vol_out_count     = 0
  max_vol_out       = 0

  n_msg_vol_each_proc = 0

  ! prepare for receiving message from write_seismograms
  if (VERBOSE) then
    call flush_stdout()
    print *
    print *, "io_server: io node rank:", myrank," is waiting for the first message"
    print *
    call flush_stdout()
  endif

  ! obtain DT and NSTEP for each event
  call recv_i_inter(tmp_1d_iarr, 1, 0, io_tag_nstep)
  NSTEP = tmp_1d_iarr(1)

  call recv_dp_inter(tmp_1d_darr, 1, 0, io_tag_dt)
  DT = tmp_1d_darr(1)

  call recv_i_inter(tmp_1d_iarr, 1, 0, io_tag_seismo_length)
  nlength_seismogram = tmp_1d_iarr(1)

  !
  ! initialization seismo
  !

  if (myrank == 0) then
    ! get receiver info from compute nodes
    call get_receiver_info(islice_num_rec_local)

    ! initialize output file for seismo
    call write_output_hdf5_seismogram_init()

    ! count the number of procs having receivers (n_procs_with_rec)
    ! and the number of receivers on each procs (islice...)
    call count_nprocs_with_recs(islice_num_rec_local)

    ! check the seismo types to be saved
    call count_seismo_type()

    ! allocate temporal arrays for seismo signals
    call allocate_seismo_arrays(islice_num_rec_local)

    ! initialize receive count
    ! count the number of messages being sent
    n_recv_msg_seismo = n_procs_with_rec * n_msg_seismo_each_proc * n_seismo_type

    if (NTSTEP_BETWEEN_OUTPUT_SEISMOS < NSTEP) then
      max_seismo_out = int(NSTEP/NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      if (mod(NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS) /= 0) max_seismo_out = max_seismo_out+1
    else
      max_seismo_out = 1
    endif

    ! receive the global id of received
    call recv_id_rec(islice_num_rec_local)

    !
    ! initialize surface movie
    !
    if (MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
      call surf_mov_io_init_hdf5(nfaces_perproc, surface_offset)
      if (VERBOSE) print *, "io_server: surf move init done"

      if (.not. USE_HIGHRES_FOR_MOVIES) then
        num_nodes = size(surf_x)
      else
        num_nodes = size(surf_x_aug)
      endif

      ! surface movie
      if (MOVIE_SURFACE) then
        n_recv_msg_surf = n_msg_surf_each_proc * NPROC
        call write_xdmf_surface_header(num_nodes)
        if (VERBOSE) print *, "io_server: surface header done"
        max_surf_out = int(NSTEP/NTSTEP_BETWEEN_FRAMES)
      endif

      !
      ! initialize shakemap
      !
      if (CREATE_SHAKEMAP) then
        call shakemap_io_init_hdf5()
        if (VERBOSE) print *, "io_server: shakemap init done"
        n_recv_msg_shake = n_msg_shake_each_proc * NPROC
        max_shake_out = 1  ! only once at final timestep
      endif
    endif

  endif ! endif myrank == 0

  !
  ! initialize volume movie
  !
  if (MOVIE_VOLUME) then
    call movie_volume_io_init_hdf5(nelm_par_proc,nglob_par_proc,nglob_par_proc_offset,nglob_par_io_offset,nelm_par_io_offset)
    if (VERBOSE) print *, "io_server: movie volume init done"

    n_recv_msg_vol = n_msg_vol_each_proc*nproc_io
    max_vol_out    = int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ! initialize flags for the value types to be written out
    val_type_mov(:) = .false.

    ! allocate dumping arrays
    allocate(vd_pres%req(0:nproc_io-1),  vd_divglob%req(0:nproc_io-1), vd_div%req(0:nproc_io-1), &
             vd_curlx%req(0:nproc_io-1), vd_curly%req(0:nproc_io-1),   vd_curlz%req(0:nproc_io-1), &
             vd_velox%req(0:nproc_io-1), vd_veloy%req(0:nproc_io-1),   vd_veloz%req(0:nproc_io-1),stat=ier)
    if (ier /= 0) stop 'Error allocating vd_pres%req arrays'
    vd_pres%req(:)   = VAL_NOT_ASSIGNED
    vd_divglob%req(:)= VAL_NOT_ASSIGNED
    vd_div%req(:)    = VAL_NOT_ASSIGNED
    vd_curlx%req(:)  = VAL_NOT_ASSIGNED
    vd_curly%req(:)  = VAL_NOT_ASSIGNED
    vd_curlz%req(:)  = VAL_NOT_ASSIGNED
    vd_velox%req(:)  = VAL_NOT_ASSIGNED
    vd_veloy%req(:)  = VAL_NOT_ASSIGNED
    vd_veloz%req(:)  = VAL_NOT_ASSIGNED

    allocate(vd_pres%d1darr(nglob_par_io_offset(myrank)), &
             vd_divglob%d1darr(nglob_par_io_offset(myrank)), &
             vd_div%d1darr(nglob_par_io_offset(myrank)), &
             vd_curlx%d1darr(nglob_par_io_offset(myrank)), &
             vd_curly%d1darr(nglob_par_io_offset(myrank)), &
             vd_curlz%d1darr(nglob_par_io_offset(myrank)), &
             vd_velox%d1darr(nglob_par_io_offset(myrank)), &
             vd_veloy%d1darr(nglob_par_io_offset(myrank)), &
             vd_veloz%d1darr(nglob_par_io_offset(myrank)),stat=ier)
    if (ier /= 0) stop 'Error allocating vd_pres%d1darr arrays'
    vd_pres%d1darr(:) = 0.0; vd_divglob%d1darr(:) = 0.0; vd_div%d1darr(:) = 0.0
    vd_curlx%d1darr(:) = 0.0; vd_curly%d1darr(:) = 0.0; vd_curlz%d1darr(:) = 0.0
    vd_velox%d1darr(:) = 0.0; vd_veloy%d1darr(:) = 0.0; vd_veloz%d1darr(:) = 0.0

    if (MOVIE_VOLUME_STRESS) then
      ! allocate dumping arrays
      allocate(vd_stressxx%req(0:nproc_io-1), vd_stressyy%req(0:nproc_io-1), vd_stresszz%req(0:nproc_io-1), &
               vd_stressxy%req(0:nproc_io-1), vd_stressxz%req(0:nproc_io-1), vd_stressyz%req(0:nproc_io-1), stat=ier)
      if (ier /= 0) stop 'Error allocating vd_stress%req arrays'
      vd_stressxx%req(:) = VAL_NOT_ASSIGNED
      vd_stressyy%req(:) = VAL_NOT_ASSIGNED
      vd_stresszz%req(:) = VAL_NOT_ASSIGNED
      vd_stressxy%req(:) = VAL_NOT_ASSIGNED
      vd_stressxz%req(:) = VAL_NOT_ASSIGNED
      vd_stressyz%req(:) = VAL_NOT_ASSIGNED

      allocate(vd_stressxx%d1darr(nglob_par_io_offset(myrank)), &
               vd_stressyy%d1darr(nglob_par_io_offset(myrank)), &
               vd_stresszz%d1darr(nglob_par_io_offset(myrank)), &
               vd_stressxy%d1darr(nglob_par_io_offset(myrank)), &
               vd_stressxz%d1darr(nglob_par_io_offset(myrank)), &
               vd_stressyz%d1darr(nglob_par_io_offset(myrank)),stat=ier)
      if (ier /= 0) stop 'Error allocating vd_stress%d1darr arrays'
      vd_stressxx%d1darr(:) = 0.0; vd_stressyy%d1darr(:) = 0.0; vd_stresszz%d1darr(:) = 0.0
      vd_stressxy%d1darr(:) = 0.0; vd_stressxz%d1darr(:) = 0.0; vd_stressyz%d1darr(:) = 0.0
    endif

  endif ! if MOVIE_VOLUME

  if (VERBOSE .and. myrank == 0) then
    call flush_stdout()
    print *
    call flush_stdout()
  endif

  !
  ! idling loop
  !
  do while (seismo_out_count < max_seismo_out .or. &
            surf_out_count < max_surf_out .or. &
            shake_out_count < max_shake_out .or. &
            vol_out_count < max_vol_out)

    ! waiting for a MPI message
    call idle_mpi_io(status)

    tag = status(my_status_tag)
    tag_src = status(my_status_source)

    ! debug output
    if (VERBOSE) then
      print '(a,i2,2(a,i4),8(a,i2))'," io_server: msg: " , tag , " src/recv rank: ", tag_src, "/",myrank, &
              "  counters, seismo: " , rec_count_seismo, "/"      , n_recv_msg_seismo, &
              ", surf: " , rec_count_surf  , "/"      , n_recv_msg_surf, &
              ", shake: ", rec_count_shake , "/"      , n_recv_msg_shake, &
              ", vol: "  , rec_count_vol   , "/"      , n_recv_msg_vol
    endif

    !
    ! receive seismograms
    !
    if (tag == io_tag_seismo_body_disp .or. &
        tag == io_tag_seismo_body_velo .or. &
        tag == io_tag_seismo_body_acce .or. &
        tag == io_tag_seismo_body_pres) then
      call recv_seismo_data(status,islice_num_rec_local)
      rec_count_seismo = rec_count_seismo+1
    endif

    !
    ! receive surface movie data
    !
    if (MOVIE_SURFACE) then
      if (tag == io_tag_surface_ux .or. &
          tag == io_tag_surface_uy .or. &
          tag == io_tag_surface_uz) then
        call recv_surf_data(status, nfaces_perproc, surface_offset)
        rec_count_surf = rec_count_surf+1
      endif
    endif

    !
    ! receive shakemap data
    !
    if (CREATE_SHAKEMAP) then
      if (tag == io_tag_shake_ux .or. &
          tag == io_tag_shake_uy .or. &
          tag == io_tag_shake_uz) then
        call recv_shake_data(status, nfaces_perproc, surface_offset)
        rec_count_shake = rec_count_shake+1
      endif
    endif

    !
    ! receive volume movie data
    !
    if (MOVIE_VOLUME) then
      if (tag == io_tag_vol_pres .or. &
          tag == io_tag_vol_divglob .or. &
          tag == io_tag_vol_div .or. &
          tag == io_tag_vol_curlx .or. &
          tag == io_tag_vol_curly .or. &
          tag == io_tag_vol_curlz .or. &
          tag == io_tag_vol_velox .or. &
          tag == io_tag_vol_veloy .or. &
          tag == io_tag_vol_veloz .or. &
          tag == io_tag_vol_stressxx .or. &
          tag == io_tag_vol_stressyy .or. &
          tag == io_tag_vol_stresszz .or. &
          tag == io_tag_vol_stressxy .or. &
          tag == io_tag_vol_stressxz .or. &
          tag == io_tag_vol_stressyz &
          ) then
        it_io = NTSTEP_BETWEEN_FRAMES*(vol_out_count+1)
        call recv_volume_data(status, val_type_mov, nglob_par_proc_offset)
        rec_count_vol = rec_count_vol+1
        ! finish gathering the whole data at each time step
      endif
    endif

    !
    ! check if all data is collected then write
    !

    ! write seismo
    if (rec_count_seismo == n_recv_msg_seismo .and. myrank == 0) then
      ! only io process with rank 0 writes out
      if (VERBOSE) print *, "io_server: seismo out ",seismo_out_count,"rank",myrank

      it_offset = seismo_out_count * NTSTEP_BETWEEN_OUTPUT_SEISMOS ! calculate the offset of timestep
      call write_seismograms_io_hdf5(it_offset)

      rec_count_seismo = 0 ! reset the counter then wait for the messages of next iteration.
      seismo_out_count = seismo_out_count+1
      if (VERBOSE) print *, "io_server: seismo out done. at it_offset = ",it_offset
    endif

    ! write surf movie
    if (MOVIE_SURFACE .and. rec_count_surf == n_recv_msg_surf .and. myrank == 0) then
      ! only io process with rank 0 writes out
      if (VERBOSE) print *, "io_server: surf out ",surf_out_count,"rank",myrank

      it_io = NTSTEP_BETWEEN_FRAMES * (surf_out_count+1)
      call write_surf_io_hdf5(it_io)
      ! write out xdmf at each timestep
      if (myrank == 0) call write_xdmf_surface_body(it_io, num_nodes)

      rec_count_surf = 0 ! reset counter
      surf_out_count = surf_out_count+1
      if (VERBOSE) print *, "io_server: surface write done. at it_io = ", it_io
    endif

    ! write volume movie
    if (MOVIE_VOLUME .and. rec_count_vol == n_recv_msg_vol) then
      ! all io processes write out
      if (VERBOSE) print *, "io_server: vol out ",vol_out_count,"rank",myrank

      ! wait all vol data are reached (for all io processes)
      call wait_volume_recv()

      ! write dumped vol data
      it_io = NTSTEP_BETWEEN_FRAMES * (vol_out_count+1)
      call write_vol_data_io_hdf5(it_io,val_type_mov,nglob_par_io_offset)
      if (vol_out_count == 0) then
        ! create xdmf header file
        if (myrank == 0) then
          if (NtoM) then
            call write_xdmf_vol_io_mton_hdf5(val_type_mov,nglob_par_io_offset,nelm_par_io_offset)
          else
            call write_xdmf_vol_io_hdf5(val_type_mov)
          endif
        endif
      endif

      rec_count_vol = 0 ! reset counter
      vol_out_count = vol_out_count+1
      if (VERBOSE) print *, "io_server: volume write done. at it_io = ", it_io
    endif

    ! write shakemap
    if (CREATE_SHAKEMAP .and. rec_count_shake == n_recv_msg_shake .and. myrank == 0) then
      ! only io process with rank 0 writes out
      if (VERBOSE) print *, "io_server: shake out ",shake_out_count

      call write_shakemap_io_hdf5()
      ! write out xdmf (last timestep)
      call write_xdmf_shakemap(num_nodes)

      rec_count_shake = 0
      shake_out_count = shake_out_count+1
      if (VERBOSE) print *, "io_server: shakemap write done."
    endif
  enddo
  !
  !  end of idling loop
  !

  ! deallocate arrays
  call deallocate_io_arrays()

  if (VERBOSE .and. myrank == 0) then
    call flush_stdout()
    print *
    call flush_stdout()
  endif
  if (VERBOSE) print *,"io_server: rank",myrank,"do_io_start_idle done"
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine do_io_start_idle() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'do_io_start_idle() called without HDF5 compilation support'

#endif

  end subroutine do_io_start_idle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine pass_info_to_io()

#ifdef USE_HDF5

  use constants, only: NGLLX,NGLLY,NGLLZ,myrank,NB_RUNS_ACOUSTIC_GPU

  use shared_parameters, only: NPROC,DT,NSTEP, &
    MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,MOVIE_VOLUME_STRESS, &
    ACOUSTIC_SIMULATION,ELASTIC_SIMULATION,POROELASTIC_SIMULATION, &
    WRITE_SEISMOGRAMS_BY_MAIN

  use specfem_par, only: NSPEC_AB,NGLOB_AB, &
    nrec,nrec_local,number_receiver_global,t0, &
    nlength_seismogram, &
    xstore,ystore,zstore

  use specfem_par_movie, only: faces_surface_offset,nfaces_perproc_surface, &
    store_val_x_all,store_val_y_all,store_val_z_all

  implicit none
  integer :: nrec_store,irec,irec_local
  integer :: n_msg_vol_each_proc
  integer :: i_ionod,ier
  integer,dimension(:), allocatable :: tmp_irec
  integer, dimension(:,:), allocatable :: elm_conn_loc

  ! initialization of io node from compute node side
  n_msg_vol_each_proc = 0

  ! pass DT NSTEP of each event which may differ from DATA/Par_file
  if (myrank == 0) then
    do i_ionod = 0,HDF5_IO_NODES-1
      call send_i_inter((/NSTEP/), 1, i_ionod, io_tag_nstep)
      call send_dp_inter((/DT/), 1, i_ionod, io_tag_dt)
      call send_i_inter((/nlength_seismogram/), 1, i_ionod, io_tag_seismo_length)
      !debug
      if (VERBOSE) print *,"io_server: rank",myrank,"pass info. DT send to",i_ionod
    enddo
  endif

  ! send the receiver information to io node for the outputs of seismo signals
  ! seismo data is out only from rank=0 of io nodes
  if (myrank == 0) then
    ! send nrec
    call send_i_inter((/nrec/), 1, 0, io_tag_num_recv)
    ! send t0
    call send_dp_inter((/t0/), 1, 0, io_tag_seismo_tzero)
  endif

  ! set number of traces in full seismogram array
  if (WRITE_SEISMOGRAMS_BY_MAIN) then
    ! main process collects all traces
    if (myrank == 0) then
      ! explicit for pressure/acoustic runs and NB_RUNS_ACOUSTIC_GPU > 1
      nrec_store = nrec * NB_RUNS_ACOUSTIC_GPU
    else
      ! dummy for secondary processes
      nrec_store = 0
    endif
  else
    ! each process outputs its local traces
    nrec_store = nrec_local * NB_RUNS_ACOUSTIC_GPU
  endif
  ! send the number of local receiver to the io node
  call send_i_inter((/nrec_store/), 1, 0, io_tag_local_rec)
  !debug
  if (VERBOSE) print *,"io_server: rank",myrank,"pass info. nrec_store",nrec_store,"sent to",0

  if (nrec_store > 0) then
    ! temporary array with irec info
    allocate(tmp_irec(nrec_store),stat=ier)
    if (ier /= 0) stop 'Error allocating tmp_irec array'
    tmp_irec(:) = 0

    ! send global id of stations
    if (WRITE_SEISMOGRAMS_BY_MAIN) then
      ! main process collects all traces
      ! loop on all the receivers
      do irec_local = 1,nrec_store
        ! get global number of that receiver
        if (NB_RUNS_ACOUSTIC_GPU == 1) then
          irec = irec_local
        else
          ! NB_RUNS_ACOUSTIC_GPU > 1
          if (mod(irec_local,nrec) == 0) then
            irec = nrec
          else
            irec = mod(irec_local,nrec)
          endif
        endif
        tmp_irec(irec_local) = irec
      enddo
    else
      ! each process outputs its local traces
      ! send global ids of local receivers (integer array)
      do irec_local = 1,nrec_store    ! nrec_store == nrec_local * NB_RUNS_ACOUSTIC_GPU
        ! get global number of that receiver
        if (NB_RUNS_ACOUSTIC_GPU == 1) then
          irec = number_receiver_global(irec_local)
        else
          ! NB_RUNS_ACOUSTIC_GPU > 1
          ! if irec_local is a multiple of nrec then mod(irec_local,nrec) == 0 and
          ! the array access at number_receiver_global would be invalid;
          ! for those cases we want irec associated to irec_local == nrec_store
          if (mod(irec_local,nrec_local) == 0) then
            irec = number_receiver_global(nrec_local)
          else
            irec = number_receiver_global(mod(irec_local,nrec_local))
          endif
        endif
        tmp_irec(irec_local) = irec
      enddo
    endif
    call send_i_inter(tmp_irec,nrec_store,0,io_tag_seismo_ids_rec)
    ! free temporary array
    deallocate(tmp_irec)
  endif

  ! surface movie/vol/shakemap
  ! surface and shakemap data will be output from the first rank of io nodes
  if (myrank == 0) then
    if (MOVIE_SURFACE .or. CREATE_SHAKEMAP) then
      ! send nfaces_perproc_surface
      call send_i_inter(nfaces_perproc_surface,NPROC, 0, io_tag_surface_nfaces)
      ! send faces_surface_offset
      call send_i_inter(faces_surface_offset,NPROC, 0, io_tag_surface_offset)
      ! send size of store_val
      call send_i_inter((/size(store_val_x_all)/), 1, 0, io_tag_surface_coord_len)
      ! send store_val_x/y/z_all
      call sendv_cr_inter(store_val_x_all,size(store_val_x_all), 0, io_tag_surface_x)
      call sendv_cr_inter(store_val_y_all,size(store_val_y_all), 0, io_tag_surface_y)
      call sendv_cr_inter(store_val_z_all,size(store_val_z_all), 0, io_tag_surface_z)
      !debug
      if (VERBOSE) print *,"io_server: rank",myrank,"pass info. surface xyz size",size(store_val_x_all)
    endif

    if (MOVIE_VOLUME) then
      ! count number of messages for volume movie
      if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
        n_msg_vol_each_proc = n_msg_vol_each_proc+1 ! pressure
      endif
      if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
        if (ELASTIC_SIMULATION) n_msg_vol_each_proc = n_msg_vol_each_proc+1 ! div_glob
        n_msg_vol_each_proc = n_msg_vol_each_proc+4 ! div, curl_x, curl_y, curl_z
        if (MOVIE_VOLUME_STRESS) n_msg_vol_each_proc = n_msg_vol_each_proc+6 ! stress_xx,..
      endif
      if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
        n_msg_vol_each_proc = n_msg_vol_each_proc+3 ! velocity_x,velocity_y,velocity_z
      endif
      ! send the number of messages
      do i_ionod = 0,HDF5_IO_NODES-1
        call send_i_inter((/n_msg_vol_each_proc/),1,i_ionod,io_tag_vol_nmsg)
      enddo
    endif
  endif ! endif myrank == 0

  if (MOVIE_VOLUME) then
    ! send the compute node list to io node
    call send_i_inter((/0/),1,dest_ionod,io_tag_vol_sendlist)

    do i_ionod = 0,HDF5_IO_NODES-1
      ! send nspec and nglob in each process
      call send_i_inter((/NSPEC_AB/),1,i_ionod,io_tag_vol_nspec)
      call send_i_inter((/NGLOB_AB/),1,i_ionod,io_tag_vol_nglob)
      ! send the id of ionode which is the destination of this compute node
      call send_i_inter((/dest_ionod/), 1, i_ionod, io_tag_vol_ioid)
    enddo

    ! send mesh information to reconstract io node based mesh database on volume movie h5 file
    ! By this way, the required storage size will be larger but for taking some advantage
    ! to reduce the IO time.
    allocate(nglob_offset(0:NPROC-1), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array nglob_offset')
    if (ier /= 0) stop 'error allocating arrays for nglob_offset'
    nglob_offset(:) = 0

    ! gather nglobs of each proc
    call gather_all_all_singlei(NGLOB_AB,nglob_offset,NPROC)

    ! connectivity
    allocate(elm_conn_loc(9,NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)),stat=ier)
    if (ier /= 0) stop 'Error allocating elm_conn_loc array'
    elm_conn_loc(:,:) = 0

    ! create connectivity dataset (it's node base but element base)
    call get_conn_for_movie(elm_conn_loc, 0)

    ! send elm_conn_loc
    call send_i_inter(elm_conn_loc(1,1),size(elm_conn_loc),dest_ionod, io_tag_vol_elmconn)

    ! send xstore, ystore, zstore
    call sendv_cr_inter(xstore, size(xstore), dest_ionod, io_tag_vol_nodex)
    call sendv_cr_inter(ystore, size(ystore), dest_ionod, io_tag_vol_nodey)
    call sendv_cr_inter(zstore, size(zstore), dest_ionod, io_tag_vol_nodez)

    deallocate(nglob_offset, elm_conn_loc)

    !debug
    if (VERBOSE) print *,"io_server: rank",myrank,"pass info. movie sent to",dest_ionod
  endif

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine pass_info_to_io() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'pass_info_to_io() called without HDF5 compilation support'

#endif

  end subroutine pass_info_to_io

!
!-------------------------------------------------------------------------------------------------
!
#ifdef USE_HDF5

  subroutine recv_volume_data(status, val_type_mov, nglob_par_proc_offset)

  use constants, only: my_status_size,my_status_source,my_status_tag
  implicit none

  integer, intent(in)                  :: status(my_status_size)
  logical, dimension(6), intent(inout) :: val_type_mov
  integer, dimension(0:nproc_io-1), intent(in) :: nglob_par_proc_offset ! offset intra io node info

  integer :: sender_glob, sender_loc, tag, msg_size, iglo_sta, iglo_end

  sender_glob = status(my_status_source)
  sender_loc  = id_proc_glob2loc(sender_glob)

  tag         = status(my_status_tag)

  ! get message size
  call world_get_size_msg(status,msg_size)

  iglo_sta = nglob_par_proc_offset(sender_loc)+1
  iglo_end = iglo_sta + msg_size -1

  ! receive data and store to dump arrays
  !
  ! here necessary to check if file write has already finished before overwriting dump arrays
  !
  if (tag == io_tag_vol_pres) then
    val_type_mov(1) = .true.
    call irecvv_cr_inter(vd_pres%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_pres%req(sender_loc))
  ! div
  else if (tag == io_tag_vol_divglob) then
    val_type_mov(2) = .true.
    call irecvv_cr_inter(vd_divglob%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_divglob%req(sender_loc))
  else if (tag == io_tag_vol_div) then
    val_type_mov(3) = .true.
    call irecvv_cr_inter(vd_div%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_div%req(sender_loc))
  ! curl
  else if (tag == io_tag_vol_curlx) then
    val_type_mov(4) = .true.
    call irecvv_cr_inter(vd_curlx%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_curlx%req(sender_loc))
  else if (tag == io_tag_vol_curly) then
    val_type_mov(4) = .true.
    call irecvv_cr_inter(vd_curly%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_curly%req(sender_loc))
  else if (tag == io_tag_vol_curlz) then
    val_type_mov(4) = .true.
    call irecvv_cr_inter(vd_curlz%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_curlz%req(sender_loc))
  ! veloc
  else if (tag == io_tag_vol_velox) then
    val_type_mov(5) = .true.
    call irecvv_cr_inter(vd_velox%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_velox%req(sender_loc))
  else if (tag == io_tag_vol_veloy) then
    val_type_mov(5) = .true.
    call irecvv_cr_inter(vd_veloy%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_veloy%req(sender_loc))
  else if (tag == io_tag_vol_veloz) then
    val_type_mov(5) = .true.
    call irecvv_cr_inter(vd_veloz%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_veloz%req(sender_loc))
  ! stress
  else if (tag == io_tag_vol_stressxx) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stressxx%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stressxx%req(sender_loc))
  else if (tag == io_tag_vol_stressyy) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stressyy%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stressyy%req(sender_loc))
  else if (tag == io_tag_vol_stresszz) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stresszz%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stresszz%req(sender_loc))
  else if (tag == io_tag_vol_stressxy) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stressxy%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stressxy%req(sender_loc))
  else if (tag == io_tag_vol_stressxz) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stressxz%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stressxz%req(sender_loc))
  else if (tag == io_tag_vol_stressyz) then
    val_type_mov(6) = .true.
    call irecvv_cr_inter(vd_stressyz%d1darr(iglo_sta:iglo_end),msg_size,sender_glob,tag,vd_stressyz%req(sender_loc))
  endif

  end subroutine recv_volume_data

#endif

!-------------------------------------------------------------------------------------------------
!
! MPI communications
!
!-------------------------------------------------------------------------------------------------


  subroutine wait_all_send()

#ifdef USE_HDF5

  implicit none

  integer :: ireq

  ! surface movie
  if (n_req_surf /= 0) then
    ! wait till all mpi_isends are finished
    do ireq = 1,n_req_surf
      call wait_req(req_dump_surf(ireq))
    enddo
  endif

  ! volume movie
  if (n_req_vol /= 0) then
    ! wait till all mpi_isends are finished
    do ireq = 1,n_req_vol
      call wait_req(req_dump_vol(ireq))
    enddo
  endif

  n_req_surf = 0; n_req_vol = 0

  call synchronize_all()

#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 I/O server routine wait_all_send() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'wait_all_send() called without HDF5 compilation support'

#endif

  end subroutine wait_all_send



!-------------------------------------------------------------------------------
!
! HDF5 io routines (only available with HDF5 compilation support)
!
!-------------------------------------------------------------------------------

#if defined(USE_HDF5)
! only available with HDF5 compilation support


  subroutine idle_mpi_io(status)

! wait for an arrival of any MPI message

  use constants, only: my_status_size

  implicit none

  integer, intent(inout) :: status(my_status_size)

  call world_probe_any_inter(status)

  end subroutine idle_mpi_io

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_volume_recv()

  implicit none

  integer :: i

  do i = 0, nproc_io-1
    if (vd_pres%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_pres%req(i))
    if (vd_divglob%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_divglob%req(i))
    if (vd_div%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_div%req(i))
    if (vd_curlx%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_curlx%req(i))
    if (vd_curly%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_curly%req(i))
    if (vd_curlz%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_curlz%req(i))
    if (vd_velox%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_velox%req(i))
    if (vd_veloy%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_veloy%req(i))
    if (vd_veloz%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_veloz%req(i))
    if (vd_stressxx%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stressxx%req(i))
    if (vd_stressyy%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stressyy%req(i))
    if (vd_stresszz%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stresszz%req(i))
    if (vd_stressxy%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stressxy%req(i))
    if (vd_stressxz%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stressxz%req(i))
    if (vd_stressyz%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_stressyz%req(i))
  enddo

  end subroutine wait_volume_recv


!
!-------------------------------------------------------------------------------------------------
!

  subroutine count_nprocs_with_recs(islice_num_rec_local)

! counts number of local receivers for each slice

  use shared_parameters, only: NPROC

  implicit none

  integer, dimension(0:NPROC-1) :: islice_num_rec_local
  integer                       :: iproc

  do iproc = 0, NPROC-1
    if (islice_num_rec_local(iproc) > 0) &
      n_procs_with_rec = n_procs_with_rec + 1
  enddo

  end subroutine count_nprocs_with_recs

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_shake_data(status, nfaces_perproc, surface_offset)

  use constants, only: my_status_size,my_status_source,my_status_tag
  use shared_parameters, only: NPROC

  implicit none

  integer, dimension(0:NPROC-1), intent(in)         :: nfaces_perproc, surface_offset
  integer                                           :: ier, sender, tag, i
  integer, intent(in)                               :: status(my_status_size)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_array

  sender = status(my_status_source)
  tag    = status(my_status_tag)

  allocate(temp_array(nfaces_perproc(sender)), stat=ier)

  call recvv_cr_inter(temp_array,nfaces_perproc(sender),sender,tag)

  if (tag == io_tag_shake_ux) then
    do i = 1, size(temp_array)
      shake_ux(i+surface_offset(sender)) = temp_array(i)
    enddo
  else if (tag == io_tag_shake_uy) then
    do i = 1, size(temp_array)
      shake_uy(i+surface_offset(sender)) = temp_array(i)
    enddo
  else if (tag == io_tag_shake_uz) then
    do i = 1, size(temp_array)
      shake_uz(i+surface_offset(sender)) = temp_array(i)
    enddo
  endif

  deallocate(temp_array)

  end subroutine recv_shake_data

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_surf_data(status, nfaces_perproc, surface_offset)

  use constants, only: my_status_size,my_status_source,my_status_tag
  use shared_parameters, only: NPROC

  implicit none

  integer, intent(in)                               :: status(my_status_size)
  integer, dimension(0:NPROC-1), intent(in)         :: nfaces_perproc, surface_offset
  ! local parameters
  integer                                           :: ier, sender, tag, i
  integer                                           :: msg_size
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_array

  sender = status(my_status_source)
  tag    = status(my_status_tag)

  allocate(temp_array(nfaces_perproc(sender)), stat=ier)
  temp_array(:) = 0._CUSTOM_REAL

  call world_get_size_msg(status, msg_size)

  call recvv_cr_inter(temp_array,nfaces_perproc(sender),sender,tag)

  if (tag == io_tag_surface_ux) then
    do i = 1, size(temp_array)
      surf_ux(i+surface_offset(sender)) = temp_array(i)
    enddo
  else if (tag == io_tag_surface_uy) then
    do i = 1, size(temp_array)
      surf_uy(i+surface_offset(sender)) = temp_array(i)
    enddo
  else if (tag == io_tag_surface_uz) then
    do i = 1, size(temp_array)
      surf_uz(i+surface_offset(sender)) = temp_array(i)
    enddo
  endif

  deallocate(temp_array)

  end subroutine recv_surf_data


!
!-------------------------------------------------------------------------------------------------
!

!
! seismo
!
  subroutine get_receiver_info(islice_num_rec_local)

  use shared_parameters, only: NPROC
  use specfem_par, only: t0,nrec

  implicit none

  integer                      :: iproc
  integer, dimension(1)        :: tmp_1d_iarr
  integer,dimension(0:NPROC-1) :: islice_num_rec_local
  double precision, dimension(1) :: tmp_1d_darr

  call recv_i_inter(tmp_1d_iarr, 1, 0, io_tag_num_recv)
  nrec = tmp_1d_iarr(1)

  call recv_dp_inter(tmp_1d_darr, 1, 0, io_tag_seismo_tzero)
  t0 = tmp_1d_darr(1)

  do iproc = 0, NPROC-1
    call recv_i_inter(tmp_1d_iarr, 1, iproc, io_tag_local_rec)
    islice_num_rec_local(iproc) = tmp_1d_iarr(1)
  enddo

  end subroutine get_receiver_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine allocate_seismo_arrays(islice_num_rec_local)

  use constants, only: NDIM

  use shared_parameters, only: NPROC,NSTEP, &
    SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY, &
    SAVE_SEISMOGRAMS_ACCELERATION,SAVE_SEISMOGRAMS_PRESSURE

  use specfem_par, only: nlength_seismogram,nrec

  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local
  integer                                   :: ier, max_num_rec, nstep_temp

  if (nlength_seismogram >= NSTEP) then
    nstep_temp = NSTEP
  else
    nstep_temp = nlength_seismogram
  endif

  ! allocate id_rec_globs for storing global id of receivers
  max_num_rec = maxval(islice_num_rec_local)
  allocate(id_rec_globs(max_num_rec,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array id_rec_globs')
  if (ier /= 0) stop 'error allocating array id_rec_globs'
  ! initialize
  id_rec_globs(:,:) = 0

  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    allocate(seismo_displ(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seimo_disp')
    if (ier /= 0) stop 'error allocating array seismo_displ'
    seismo_displ(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    allocate(seismo_veloc(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_veloc')
    if (ier /= 0) stop 'error allocating array seismo_veloc'
    seismo_veloc(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    allocate(seismo_accel(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_accel')
    if (ier /= 0) stop 'error allocating array seismo_accel'
    seismo_accel(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    allocate(seismo_press(nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_press')
    if (ier /= 0) stop 'error allocating array seismo_press'
    seismo_press(:,:) = 0._CUSTOM_REAL
  endif

  end subroutine allocate_seismo_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine deallocate_io_arrays()

  use shared_parameters, only: MOVIE_VOLUME,USE_HIGHRES_FOR_MOVIES,MOVIE_VOLUME_STRESS

  implicit none

  ! seismo
  if (allocated(id_rec_globs)) deallocate(id_rec_globs)
  if (allocated(seismo_displ)) deallocate(seismo_displ)
  if (allocated(seismo_veloc)) deallocate(seismo_veloc)
  if (allocated(seismo_accel)) deallocate(seismo_accel)
  if (allocated(seismo_press)) deallocate(seismo_press)

  ! surface movie
  if (allocated(surf_x)) deallocate(surf_x,surf_y,surf_z)
  if (allocated(surf_ux)) deallocate(surf_ux,surf_uy,surf_uz)

  ! shake map
  if (allocated(shake_ux)) deallocate(shake_ux,shake_uy,shake_uz)

  if (USE_HIGHRES_FOR_MOVIES) then
    ! surface movie
    if (allocated(surf_x_aug)) deallocate(surf_x_aug,surf_y_aug,surf_z_aug)
    if (allocated(surf_ux_aug)) deallocate(surf_ux_aug,surf_uy_aug,surf_uz_aug)
    ! shake map
    if (allocated(shake_ux_aug)) deallocate(shake_ux_aug,shake_uy_aug,shake_uz_aug)
  endif

  if (MOVIE_VOLUME) then
    deallocate(vd_pres%req ,vd_divglob%req,   vd_div%req, &
               vd_curlx%req,  vd_curly%req, vd_curlz%req, &
               vd_velox%req,  vd_veloy%req, vd_veloz%req)
    deallocate(vd_pres%d1darr, &
               vd_divglob%d1darr, &
               vd_div%d1darr, &
               vd_curlx%d1darr, vd_curly%d1darr, vd_curlz%d1darr, &
               vd_velox%d1darr, vd_veloy%d1darr, vd_veloz%d1darr)
    if (MOVIE_VOLUME_STRESS) then
      deallocate(vd_stressxx%req,vd_stressyy%req,vd_stresszz%req, &
                 vd_stressxy%req,vd_stressxz%req,vd_stressyz%req)
      deallocate(vd_stressxx%d1darr,vd_stressyy%d1darr,vd_stresszz%d1darr, &
                 vd_stressxy%d1darr,vd_stressxz%d1darr,vd_stressyz%d1darr)
    endif
    deallocate(id_proc_loc2glob, id_proc_glob2loc, dest_ioids)
  endif

  end subroutine deallocate_io_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine count_seismo_type()

  use shared_parameters, only: SAVE_SEISMOGRAMS_DISPLACEMENT,SAVE_SEISMOGRAMS_VELOCITY, &
    SAVE_SEISMOGRAMS_ACCELERATION,SAVE_SEISMOGRAMS_PRESSURE

  implicit none

  integer :: n_type = 0

  if (SAVE_SEISMOGRAMS_DISPLACEMENT) n_type = n_type+1
  if (SAVE_SEISMOGRAMS_VELOCITY)     n_type = n_type+1
  if (SAVE_SEISMOGRAMS_ACCELERATION) n_type = n_type+1
  if (SAVE_SEISMOGRAMS_PRESSURE)     n_type = n_type+1

  n_seismo_type = n_type

  end subroutine count_seismo_type

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_id_rec(islice_num_rec_local)

  use shared_parameters, only: NPROC

  implicit none

  integer :: sender
  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local

  do sender = 0, NPROC-1
    if (islice_num_rec_local(sender) /= 0) then
      call recv_i_inter(id_rec_globs(:,sender), size(id_rec_globs(:,sender)), sender, io_tag_seismo_ids_rec)
    endif
  enddo

  end subroutine recv_id_rec

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_seismo_data(status, islice_num_rec_local)

  use constants, only: NDIM,my_status_size,my_status_source,my_status_tag
  use shared_parameters, only: NPROC

  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local
  integer, intent(in)                       :: status(my_status_size)
  ! local parameters
  integer :: sender, nrec_passed, irec_passed, tag, id_rec_glob, ier
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismo_temp
  integer                                               :: msg_size,time_window

  sender = status(my_status_source)
  tag    = status(my_status_tag)

  nrec_passed  = islice_num_rec_local(sender)

  call world_get_size_msg(status,msg_size)

  if (nrec_passed > 0) then
    ! get seismogram vector values
    ! allocate temp array size
    time_window = int(msg_size/NDIM/nrec_passed)
    allocate(seismo_temp(NDIM,time_window,nrec_passed),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_temp')
    if (ier /= 0) stop 'error allocating array seismo_temp'
    seismo_temp(:,:,:) = 0.0

    call recvv_cr_inter(seismo_temp, msg_size, sender, tag)

    ! set local array to the global array
    do irec_passed = 1,nrec_passed
      id_rec_glob = id_rec_globs(irec_passed,sender)

      if (tag == io_tag_seismo_body_disp) then
        ! displacement
        seismo_displ(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed)

      else if (tag == io_tag_seismo_body_velo) then
        ! velocity
        seismo_veloc(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed)

      else if (tag == io_tag_seismo_body_acce) then
        ! acceleration
        seismo_accel(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed)

      else if (tag == io_tag_seismo_body_pres) then
        ! pressure (single component only)
        seismo_press(:,id_rec_glob) = seismo_temp(1,:,irec_passed)
      endif
    enddo

    ! deallocate temp array
    deallocate(seismo_temp)

  endif

  end subroutine recv_seismo_data


! USE_HDF5
#endif

end module io_server_hdf5

