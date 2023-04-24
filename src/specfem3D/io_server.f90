! module for storing info. concering io
module io_server
  use specfem_par, only: CUSTOM_REAL, NPROC

  implicit none


  logical :: NtoM = .true. ! if this option set .false. SPECFEM3D takes N-to-1
                           ! I/O strategy i.e. volume data will be out to
                           ! one single .h5 file. This gives the output files
                           ! great simplicity however I/O performance will be
                           ! degredated for collective h5 I/O.
  logical :: if_collect_server = .false. ! need to be set .false. for NtoM = .true.

  integer :: nglob_all_server, nspec_all_server !  total number of elements and nodes for all procs
  integer :: nglob_this_io, nelm_this_io ! number of nodes and elements in this IO group

  integer :: n_msg_seismo_each_proc=1,n_seismo_type=0
  integer :: n_procs_with_rec
  integer :: n_msg_surf_each_proc=3,surf_offset
  integer :: n_msg_shake_each_proc=3
  integer :: n_msg_vol_each_proc=0

  integer                                           :: size_surf_array=0, surf_xdmf_pos
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surf_x,   surf_y,   surf_z,   &
                                                       surf_ux,  surf_uy,  surf_uz,  &
                                                       shake_ux, shake_uy, shake_uz, &
                                                       surf_x_aug,   surf_y_aug,   surf_z_aug,  &
                                                       surf_ux_aug,  surf_uy_aug,  surf_uz_aug, &
                                                       shake_ux_aug, shake_uy_aug, shake_uz_aug


  real(kind=CUSTOM_REAL), dimension(:,:),   allocatable   :: seismo_pres
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable   :: seismo_disp, seismo_velo, seismo_acce
  integer, dimension(:,:), allocatable                    :: id_rec_globs

  ! output file names
  character(len=64) :: fname_h5_seismo     = ""
  character(len=64) :: fname_h5_data_surf  = ""
  character(len=64) :: fname_h5_data_vol   = ""
  character(len=64) :: fname_h5_data_shake = ""

  character(len=64) :: fname_xdmf_surf     = ""
  character(len=64) :: fname_xdmf_vol      = ""
  character(len=64) :: fname_xdmf_vol_step = ""
  character(len=64) :: fname_xdmf_shake    = ""

  integer, parameter :: VAL_NOT_ASSIGNED = 9999999
  type vol_data_dump
    integer, dimension(:), allocatable :: req !=VAL_NOT_ASSIGNED should be
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: d1darr
  end type vol_data_dump

  ! dumps for volume data
  type(vol_data_dump) :: vd_pres, vd_divglob, vd_div,  &
                         vd_curlx, vd_curly, vd_curlz, &
                         vd_velox, vd_veloy, vd_veloz
  ! local-global processor id relation
  integer, dimension(:), allocatable :: id_proc_glob2loc, id_proc_loc2glob
  ! io node id <-> proc id in compute nodes relation
  integer, dimension(:), allocatable :: dest_ioids

  ! dummy variables
  integer :: dummy_int


contains
  function i2c(k) result(str)
  !   "Convert an integer to string."
      integer, intent(in) :: k
      character(len=20) str
      write (str, "(i20)") k
      str = adjustl(str)
  end function i2c

  function r2c(k) result(str)
  !   "Convert an real to string."
      real(kind=CUSTOM_REAL), intent(in) :: k
      character(len=20) str
      write (str, *) k
      str = adjustl(str)
  end function r2c

  ! reorder and expand the input array (for only corner nodes) to high res array (all GLL)
  subroutine recompose_for_hires(arr_in, arr_out)
    use specfem_par
    implicit none

    integer                                             :: nfaces_actual
    real(kind=CUSTOM_REAL), dimension(:), intent(in)    :: arr_in
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: arr_out
    integer :: i,j,k,c,factor_face_aug=(NGLLX-1)*(NGLLY-1),npoint_per_face=NGLLX*NGLLY,npoint_corner=4

    nfaces_actual=size(arr_in)/(NGLLX*NGLLY) ! expecting NGLLX==NGLLY==NGLLZ

    c=1
    do i=0, nfaces_actual-1
      do j=0,NGLLY-2 !y
        do k=0,NGLLX-2 !x
          arr_out(c  +j*(NGLLX-1)*4+k*4)=arr_in(i*npoint_per_face+1      +k+j*NGLLX)
          arr_out(c+1+j*(NGLLX-1)*4+k*4)=arr_in(i*npoint_per_face+2      +k+j*NGLLX)
          arr_out(c+2+j*(NGLLX-1)*4+k*4)=arr_in(i*npoint_per_face+2+NGLLX+k+j*NGLLX)
          arr_out(c+3+j*(NGLLX-1)*4+k*4)=arr_in(i*npoint_per_face+1+NGLLX+k+j*NGLLX)
        enddo
      enddo
      c=c+factor_face_aug*npoint_corner
    enddo

  end subroutine recompose_for_hires

end module io_server


subroutine do_io_start_idle()
  use my_mpi
  use specfem_par
  use io_server

  implicit none

  integer :: ier

  ! vars seismo
  integer,dimension(0:NPROC-1) :: islice_num_rec_local
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: rec_count_seismo=0, n_recv_msg_seismo=0, max_num_rec,idump,max_seismo_out=0
  integer                      :: it_offset=0, seismo_out_count=0
  integer, dimension(1)        :: nrec_temp

  ! vars surface movie
  integer                       :: rec_count_surf=0, n_recv_msg_surf=0,surf_out_count=0, it_io,max_surf_out=0
  integer, dimension(0:NPROC-1) :: nfaces_perproc, surface_offset
  integer                       :: num_nodes

  ! vars shakemap
  integer :: rec_count_shake=0, n_recv_msg_shake=0, shake_out_count=0,max_shake_out=0

  ! vars volumne movie
  integer                       :: rec_count_vol=0, n_recv_msg_vol=0, vol_out_count=0, max_vol_out=0
  ! storing the number of elements and gll nodes
  integer, dimension(0:NPROC-1)   :: nelm_par_proc, nglob_par_proc, nglob_par_proc_offset
  integer, dimension(0:NIONOD-1)  :: nglob_par_io_offset, nelm_par_io_offset
  logical, dimension(5)           :: val_type_mov ! true if movie file will be created, (pressure, div_glob, div, curlxyz, velocity_xyz)


  ! prepare for receiving message from write_seismograms
  print *, "io node rank:", myrank," is waiting for the first message"

  !
  ! initialization seismo
  !

  if (myrank == 0) then
    ! get receiver info from compute nodes
    call get_receiver_info(islice_num_rec_local)

    ! initialize output file for seismo
    call do_io_seismogram_init()

    ! count the number of procs having receivers (n_procs_with_rec)
    ! and the number of receivers on each procs (islice...)
    call count_nprocs_with_recs(islice_num_rec_local)

    ! check the seismo types to be saved
    call count_seismo_type()

    ! allocate temporal arrays for seismo signals
    call allocate_seismo_arrays(islice_num_rec_local)

    ! initialize receive count
    ! count the number of messages being sent
    n_recv_msg_seismo = n_procs_with_rec*n_msg_seismo_each_proc*n_seismo_type

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
      call surf_mov_init(nfaces_perproc, surface_offset)
      if (.not. USE_HIGHRES_FOR_MOVIES) then
        num_nodes = size(surf_x)
      else
        num_nodes = size(surf_x_aug)
      endif

      if (MOVIE_SURFACE) then
        n_recv_msg_surf = n_msg_surf_each_proc*NPROC
        print *, "surf move init done"
        call write_xdmf_surface_header(num_nodes, dummy_int)

        max_surf_out = int(NSTEP/NTSTEP_BETWEEN_FRAMES)
      endif
    !
    ! initialize shakemap
    !
      if (CREATE_SHAKEMAP) then
        call shakemap_init(nfaces_perproc, surface_offset)
        n_recv_msg_shake = n_msg_shake_each_proc*NPROC
        print *, "shakemap init done"
        max_shake_out = 1
      endif
    endif

  endif ! end if myrank == 0

  !
  ! initialize volume movie
  !
  if (MOVIE_VOLUME) then
    call movie_volume_init(nelm_par_proc,nglob_par_proc,nglob_par_proc_offset,nglob_par_io_offset,nelm_par_io_offset)
    print *, "movie volume init done"
    n_recv_msg_vol = n_msg_vol_each_proc*nproc_io
    max_vol_out    = int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ! initialize flags for the value types to be written out
    val_type_mov(:) = .false.

    ! allocate dumping arrays
    allocate(vd_pres%req(0:nproc_io-1),vd_divglob%req(0:nproc_io-1),   vd_div%req(0:nproc_io-1), &
            vd_curlx%req(0:nproc_io-1),  vd_curly%req(0:nproc_io-1), vd_curlz%req(0:nproc_io-1), &
            vd_velox%req(0:nproc_io-1),  vd_veloy%req(0:nproc_io-1), vd_veloz%req(0:nproc_io-1))
    vd_pres%req(:)   = VAL_NOT_ASSIGNED
    vd_divglob%req(:)= VAL_NOT_ASSIGNED
    vd_div%req(:)    = VAL_NOT_ASSIGNED
    vd_curlx%req(:)  = VAL_NOT_ASSIGNED
    vd_curly%req(:)  = VAL_NOT_ASSIGNED
    vd_curlz%req(:)  = VAL_NOT_ASSIGNED
    vd_velox%req(:)  = VAL_NOT_ASSIGNED
    vd_veloy%req(:)  = VAL_NOT_ASSIGNED
    vd_veloz%req(:)  = VAL_NOT_ASSIGNED
    allocate(   vd_pres%d1darr(nglob_par_io_offset(myrank)),&
             vd_divglob%d1darr(nglob_par_io_offset(myrank)),&
                 vd_div%d1darr(nglob_par_io_offset(myrank)),&
               vd_curlx%d1darr(nglob_par_io_offset(myrank)),&
               vd_curly%d1darr(nglob_par_io_offset(myrank)),&
               vd_curlz%d1darr(nglob_par_io_offset(myrank)),&
               vd_velox%d1darr(nglob_par_io_offset(myrank)),&
               vd_veloy%d1darr(nglob_par_io_offset(myrank)),&
               vd_veloz%d1darr(nglob_par_io_offset(myrank)))


  endif ! if MOVIE_VOLUME

 !
 ! idling loop
 !
  do while (seismo_out_count < max_seismo_out .or. &
            surf_out_count   < max_surf_out   .or. &
            shake_out_count  < max_shake_out  .or. &
            vol_out_count    < max_vol_out)
    ! waiting for a mpi message
    call idle_mpi_io(status)

    ! debug output
    !print *,                 "msg: " , status(MPI_TAG) , " send/recv rank: ", status(MPI_SOURCE), "/",myrank, &
    !          "  counters, seismo: " , rec_count_seismo, "/"      , n_recv_msg_seismo,  &
    !                      ", surf: " , rec_count_surf  , "/"      , n_recv_msg_surf,    &
    !                      ", shake: ", rec_count_shake , "/"      , n_recv_msg_shake,   &
    !                      ", vol: "  , rec_count_vol   , "/"      , n_recv_msg_vol

    !
    ! receive seismograms
    !

    if (status(MPI_TAG) == io_tag_seismo_body_disp .or. &
        status(MPI_TAG) == io_tag_seismo_body_velo .or. &
        status(MPI_TAG) == io_tag_seismo_body_acce .or. &
        status(MPI_TAG) == io_tag_seismo_body_pres      &
    ) then
      call recv_seismo_data(status,islice_num_rec_local,rec_count_seismo)
      rec_count_seismo = rec_count_seismo+1
    endif

    !
    ! receive surface movie data
    !
    if (MOVIE_SURFACE) then
      if (status(MPI_TAG) == io_tag_surface_ux .or. &
          status(MPI_TAG) == io_tag_surface_uy .or. &
          status(MPI_TAG) == io_tag_surface_uz      &
      ) then
        call recv_surf_data(status, nfaces_perproc, surface_offset)
        rec_count_surf = rec_count_surf+1
      endif
    endif

    !
    ! receive shakemap data
    !
    if (CREATE_SHAKEMAP) then
      if (status(MPI_TAG) == io_tag_shake_ux .or. &
          status(MPI_TAG) == io_tag_shake_uy .or. &
          status(MPI_TAG) == io_tag_shake_uz      &
      ) then
        call recv_shake_data(status, nfaces_perproc, surface_offset)
        rec_count_shake = rec_count_shake+1
      endif
    endif

    !
    ! receive volume movie data
    !
    if (MOVIE_VOLUME) then
      if ( status(MPI_TAG) == io_tag_vol_pres    .or. &
           status(MPI_TAG) == io_tag_vol_divglob .or. &
           status(MPI_TAG) == io_tag_vol_div     .or. &
           status(MPI_TAG) == io_tag_vol_curlx   .or. &
           status(MPI_TAG) == io_tag_vol_curly   .or. &
           status(MPI_TAG) == io_tag_vol_curlz   .or. &
           status(MPI_TAG) == io_tag_vol_velox   .or. &
           status(MPI_TAG) == io_tag_vol_veloy   .or. &
           status(MPI_TAG) == io_tag_vol_veloz        &
      ) then
        it_io = NTSTEP_BETWEEN_FRAMES*(vol_out_count+1)
        call recv_vol_data(status,rec_count_vol,it_io, val_type_mov, nglob_par_proc_offset)
        rec_count_vol = rec_count_vol+1
        ! finish gathering the whole data at each time step
      endif
    endif

    !
    ! check if all data is collected then write
    !

    ! write seismo
    if (rec_count_seismo == n_recv_msg_seismo .and. myrank==0) then
      it_offset        = seismo_out_count*NTSTEP_BETWEEN_OUTPUT_SEISMOS ! calculate the offset of timestep
      call write_seismograms_io(it_offset)
      rec_count_seismo = 0 ! reset the counter then wait for the messages of next iteration.
      seismo_out_count = seismo_out_count+1
    endif

    ! write surf movie
    if (MOVIE_SURFACE .and. rec_count_surf == n_recv_msg_surf .and. myrank==0) then
      it_io          = NTSTEP_BETWEEN_FRAMES*(surf_out_count+1)
      call write_surf_io(it_io)
      rec_count_surf = 0 ! reset counter
      surf_out_count = surf_out_count+1

      ! write out xdmf at each timestep
      call write_xdmf_surface_body(it_io, num_nodes, dummy_int)

      print *, "surface write done at it = ", it_io
    endif

    ! write volume movie
    if (MOVIE_VOLUME .and. rec_count_vol == n_recv_msg_vol) then
      ! wait all vol data are reached
      call wait_vol_recv()

      ! write dumped vol data
      call write_vol_data(it_io,val_type_mov,nglob_par_io_offset)

      rec_count_vol = 0 ! reset counter
      vol_out_count = vol_out_count+1
      if (vol_out_count==1) then
        ! create xdmf header file
        if (myrank == 0) then
          if (NtoM) then
            call write_xdmf_vol_mton(val_type_mov,nglob_par_io_offset,nelm_par_io_offset)
          else
            call write_xdmf_vol(val_type_mov)
          endif
        endif
      endif

      print *, "volume write done at it = ", it_io

    endif

    ! write shakemap
    if (CREATE_SHAKEMAP .and. rec_count_shake == n_recv_msg_shake .and. myrank==0) then
      call write_shake_io()
      rec_count_shake = 0
      shake_out_count = shake_out_count+1
      ! write out xdmf at each timestep
      call write_xdmf_shakemap(num_nodes)

      print *, "shakemap write done"
    endif

  enddo
  !
  !  end of idling loop
  !

  ! deallocate arrays
!  call deallocate_arrays()

end subroutine do_io_start_idle


!
! volume movie
!
subroutine movie_volume_init(nelm_par_proc,nglob_par_proc,nglob_par_proc_offset,nglob_par_io_offset,nelm_par_io_offset)
  use io_server
  use my_mpi
  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic
  use specfem_par_movie
  use phdf5_utils
  implicit none

  integer :: iproc, iproc2, count=0, id_glob, comm, info, sender, dump, ier

  integer, dimension(0:NPROC-1), intent(inout)  :: nelm_par_proc, nglob_par_proc ! storing the number of elements and gll nodes
  integer, dimension(0:nproc_io-1), intent(out) :: nglob_par_proc_offset ! offset intra io node info
  integer, dimension(0:NIONOD-1), intent(out)   :: nglob_par_io_offset   ! offset inter io node info
  integer, dimension(0:nproc_io-1)              :: nelm_par_proc_offset ! offset intra io node
  integer, dimension(0:NIONOD-1), intent(out)   :: nelm_par_io_offset   ! offset inter io node
  integer                                       :: status(MPI_STATUS_SIZE)

  ! arrays for storing elm_conn and xyzstore
  !integer                                           :: nglob_this_io=0, nelm_this_io=0
  integer, dimension(:,:), allocatable              :: elm_conn_tmp
  integer, dimension(:), allocatable                :: elm_conn_subtmp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_tmp, ystore_tmp, zstore_tmp
  integer :: nelm_tmp,nglo_tmp,ielm_sta,ielm_end,iglo_sta,iglo_end
  integer :: node_id_off, node_id_off_inter
  integer :: nelms_factor = (NGLLX-1)*(NGLLY-1)*(NGLLZ-1) ! divide one spectral element to
  ! small rectangles for visualization purpose
  integer, dimension(1)                             :: tmp_1d_arr


  ! make output file
  character(len=64) :: group_name
  character(len=64) :: dset_name
  character(len=64) :: proc_str
  type(h5io)        :: h5
  h5 = h5io()

  write(proc_str, "(i6.6)") myrank

  if (NtoM) then
    ! write out one h5 par 1 IO server process
    fname_h5_data_vol = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_volume_"//trim(proc_str)//".h5"
  else
    ! write one single h5 from all IO server
    fname_h5_data_vol = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_volume.h5"
  endif
  ! initialization of h5 file
  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)

  ! initialize h5 object
  call h5_init(h5, fname_h5_data_vol)
  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)

  ! get n_msg_vol_each_proc
  call recv_i_inter(tmp_1d_arr, 1, 0, io_tag_vol_nmsg)
  n_msg_vol_each_proc = tmp_1d_arr(1)

  ! make an array of local2global relation of sending compute node ids
  allocate(id_proc_loc2glob(0:nproc_io-1))
  ! make an array of global2local relation of sending compute node ids
  allocate(id_proc_glob2loc(0:NPROC-1))
  id_proc_glob2loc(:) = -999999

  allocate(dest_ioids(0:NPROC-1))

  ! make a sender list which communicate with this io node
  print *, "number of compute node: ", nproc_io, ", for io node rank: ", myrank
  do iproc = 0, nproc_io-1
    call mpi_probe(MPI_ANY_SOURCE, io_tag_vol_sendlist, my_local_mpi_comm_inter, status, ier)
    sender = status(MPI_SOURCE)
    call recv_i_inter((/dump/),1,sender,io_tag_vol_sendlist)
    id_proc_loc2glob(iproc) = sender
    id_proc_glob2loc(sender) = iproc
  enddo

  ! gather other informations for making a volume data output
  do iproc = 0, NPROC-1
    ! get nspec and nglob from each process
    call recv_i_inter(nelm_par_proc(iproc),  1, iproc, io_tag_vol_nspec) ! NSPEC_AB
    call recv_i_inter(nglob_par_proc(iproc), 1, iproc, io_tag_vol_nglob) ! NGLOB_AB
    ! array if io node ids of each compute procs
    call recv_i_inter(dest_ioids(iproc), 1, iproc, io_tag_vol_ioid)
  enddo

  ! calculate nglob and nelm for this io_node
  nglob_par_proc_offset(:) = 0
  nelm_par_proc_offset(:)  = 0

  ! count the number of nodes and elements in this IO group
  nglob_this_io = 0
  nelm_this_io  = 0
  do iproc=0, nproc_io-1
    nglob_this_io = nglob_this_io + nglob_par_proc(id_proc_loc2glob(iproc))
    nelm_this_io  = nelm_this_io  + nelm_par_proc (id_proc_loc2glob(iproc))
    if (iproc>0) then
      do iproc2=0, iproc-1
        nglob_par_proc_offset(iproc) = nglob_par_proc_offset(iproc)+nglob_par_proc(id_proc_loc2glob(iproc2))
        nelm_par_proc_offset(iproc)  = nelm_par_proc_offset(iproc) +nelm_par_proc (id_proc_loc2glob(iproc2))
      enddo
    endif
  enddo

  ! prepare intra io node offset array
  call gather_all_all_singlei((/nglob_this_io/), nglob_par_io_offset, NIONOD)
  call gather_all_all_singlei((/nelm_this_io/), nelm_par_io_offset, NIONOD)

  ! prepare arrays for storing elm_conn_loc and xyzstore
  allocate(xstore_tmp(nglob_this_io)) ! 1d array for this io_node
  allocate(ystore_tmp(nglob_this_io))
  allocate(zstore_tmp(nglob_this_io))
  allocate(elm_conn_tmp(9,nelm_this_io*nelms_factor))

  nspec_all_server = sum(nelm_par_proc(:))
  nglob_all_server = sum(nglob_par_proc(:))

  ! node id offset inter io nodes
  if (NtoM) then
    node_id_off_inter = 0 ! no offsets for inter io nodes
  else
    node_id_off_inter = sum(nglob_par_io_offset(0:myrank-1))
  endif

  do iproc=0, nproc_io-1
    nelm_tmp = nelm_par_proc(id_proc_loc2glob(iproc))
    nglo_tmp = nglob_par_proc(id_proc_loc2glob(iproc))
    ielm_sta = nelm_par_proc_offset(iproc)+1
    ielm_end = ielm_sta + nelm_tmp-1
    iglo_sta = nglob_par_proc_offset(iproc)+1
    iglo_end = iglo_sta + nglo_tmp-1
    !print*,"iorank, iproc, ielm_sta, ielm_end, nelm_tmp, nspec_all", myrank, iproc, ielm_sta, ielm_end, nelm_tmp, nspec_all_server
    !print*,"iorank, iproc, iglo_sta, iglo_end, nglo_tmp, nglob_all", myrank, iproc, iglo_sta, iglo_end, nglo_tmp, nglob_all_server

    allocate(elm_conn_subtmp(9*nelm_tmp*nelms_factor))

    ! receive elm_conn_loc
    call recv_i_inter(elm_conn_subtmp, &
        9*nelm_tmp*nelms_factor,id_proc_loc2glob(iproc),io_tag_vol_elmconn)

    elm_conn_tmp(:,(ielm_sta-1)*nelms_factor+1:ielm_end*nelms_factor) &
        = reshape(elm_conn_subtmp,(/9,nelm_tmp*nelms_factor/))

    ! node id offset intra io node
    node_id_off = 0
    do iproc2=0,iproc-1
      node_id_off = node_id_off + nglob_par_proc(id_proc_loc2glob(iproc2))
    enddo
    !print *, "iorank, iproc, off, off_inter", myrank, iproc, node_id_off, node_id_off_inter

    elm_conn_tmp(2:9,(ielm_sta-1)*nelms_factor+1:ielm_end*nelms_factor) &
    = elm_conn_tmp(2:9,(ielm_sta-1)*nelms_factor+1:ielm_end*nelms_factor) &
    + node_id_off + node_id_off_inter

    ! receive xstore, ystore, zstore
    call recvv_cr_inter(xstore_tmp(iglo_sta:iglo_end),nglo_tmp,id_proc_loc2glob(iproc),io_tag_vol_nodex)
    call recvv_cr_inter(ystore_tmp(iglo_sta:iglo_end),nglo_tmp,id_proc_loc2glob(iproc),io_tag_vol_nodey)
    call recvv_cr_inter(zstore_tmp(iglo_sta:iglo_end),nglo_tmp,id_proc_loc2glob(iproc),io_tag_vol_nodez)

    deallocate(elm_conn_subtmp)
  enddo

  ! write mesh database in movie_volue.h5
  group_name = "mesh"

  ! write mesh info. of each IO server into each movie_volume_*.h5 file
  if (NtoM) then

    ! create a hdf5 file
    call h5_create_file(h5)
    ! create coordinate dataset
    call h5_open_or_create_group(h5, group_name)

    dset_name = "elm_conn"
    call h5_create_dataset_gen_in_group(h5, dset_name, &
                (/9,nelm_this_io*nelms_factor/), 2, 1)
    dset_name = "x"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
    dset_name = "y"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
    dset_name = "z"
    call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)

    ! write data
    dset_name = "elm_conn"
    call h5_write_dataset_2d_i_collect_hyperslab_in_group(h5,dset_name,elm_conn_tmp,&
                (/0,0/),if_collect_server)
    dset_name = "x"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,xstore_tmp,&
                (/0/),if_collect_server)
    dset_name = "y"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,ystore_tmp,&
                (/0/),if_collect_server)
    dset_name = "z"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,zstore_tmp,&
                (/0/),if_collect_server)

    call h5_close_group(h5)
    call h5_close_file(h5)

  else ! for N-to-1 strategy
    if(myrank==0) then
      ! create a hdf5 file
      call h5_create_file(h5)
      ! create coordinate dataset
      call h5_open_or_create_group(h5, group_name)

      dset_name = "elm_conn"
      call h5_create_dataset_gen_in_group(h5, dset_name, &
                  (/9,nspec_all_server*nelms_factor/), 2, 1)
      dset_name = "x"
      call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
      dset_name = "y"
      call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
      dset_name = "z"
      call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)

      call h5_close_group(h5)
      call h5_close_file(h5)

    endif
    call synchronize_all()

    call h5_open_file_p(h5)
    call h5_open_group(h5,group_name)

    dset_name = "elm_conn"
    call h5_write_dataset_2d_i_collect_hyperslab_in_group(h5,dset_name,elm_conn_tmp,&
                (/0,sum(nelm_par_io_offset(0:myrank-1))*nelms_factor/),if_collect_server)
    dset_name = "x"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,xstore_tmp,&
                (/sum(nglob_par_io_offset(0:myrank-1))/),if_collect_server)
    dset_name = "y"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,ystore_tmp,&
                (/sum(nglob_par_io_offset(0:myrank-1))/),if_collect_server)
    dset_name = "z"
    call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5,dset_name,zstore_tmp,&
                (/sum(nglob_par_io_offset(0:myrank-1))/),if_collect_server)

    call h5_close_group(h5)
    call h5_close_file(h5)

  endif ! NtoM or Nto1

  call synchronize_all()

  deallocate(elm_conn_tmp)
  deallocate(xstore_tmp, ystore_tmp, zstore_tmp)
end subroutine movie_volume_init


subroutine recv_vol_data(status, rec_count_vol, it_io, val_type_mov, nglob_par_proc_offset)
  use io_server
  use specfem_par
  use my_mpi
  implicit none

  integer, intent(in)                  :: status(MPI_STATUS_SIZE)
  integer, intent(in)                  :: rec_count_vol,it_io
  logical, dimension(5), intent(inout) :: val_type_mov
  integer, dimension(0:nproc_io-1), intent(in) :: nglob_par_proc_offset ! offset intra io node info
  integer :: sender_glob, sender_loc, ier, tag, arrsize, msgsize, iglo_sta, iglo_end

  sender_glob = status(MPI_SOURCE)
  sender_loc  = id_proc_glob2loc(sender_glob)
  tag         = status(MPI_TAG)

  ! get message size
  call get_size_msg(status,msgsize)

  iglo_sta = nglob_par_proc_offset(sender_loc)+1
  iglo_end = iglo_sta + msgsize -1

  ! receive data and store to dump arrays
  !
  ! here necessary to check if file write has already finished before overwriting dump arrays
  !
  if (tag == io_tag_vol_pres) then
    val_type_mov(1) = .true.
    call irecvv_cr_inter(vd_pres%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_pres%req(sender_loc))
  elseif (tag == io_tag_vol_divglob) then
    val_type_mov(2) = .true.
    call irecvv_cr_inter(vd_divglob%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_divglob%req(sender_loc))
  elseif (tag == io_tag_vol_div) then
    val_type_mov(3) = .true.
    call irecvv_cr_inter(vd_div%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_div%req(sender_loc))
  elseif (tag == io_tag_vol_curlx) then
    val_type_mov(4) = .true.
    call irecvv_cr_inter(vd_curlx%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_curlx%req(sender_loc))
  elseif (tag == io_tag_vol_curly) then
    call irecvv_cr_inter(vd_curly%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_curly%req(sender_loc))
  elseif (tag == io_tag_vol_curlz) then
    call irecvv_cr_inter(vd_curlz%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_curlz%req(sender_loc))
  elseif (tag == io_tag_vol_velox) then
    val_type_mov(5) = .true.
    call irecvv_cr_inter(vd_velox%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_velox%req(sender_loc))
  elseif (tag == io_tag_vol_veloy) then
    call irecvv_cr_inter(vd_veloy%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_veloy%req(sender_loc))
  elseif (tag == io_tag_vol_veloz) then
    call irecvv_cr_inter(vd_veloz%d1darr(iglo_sta:iglo_end),msgsize,sender_glob,tag,vd_veloz%req(sender_loc))
  endif

end subroutine recv_vol_data


subroutine wait_vol_recv()
  use io_server
  use my_mpi
  use constants, only: nproc_io
  implicit none

  integer :: i

  do i = 0, nproc_io-1
    if (vd_pres%req(i)    /= VAL_NOT_ASSIGNED) call wait_req(vd_pres%req(i))
    if (vd_divglob%req(i) /= VAL_NOT_ASSIGNED) call wait_req(vd_divglob%req(i))
    if (vd_div%req(i)     /= VAL_NOT_ASSIGNED) call wait_req(vd_div%req(i))
    if (vd_curlx%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_curlx%req(i))
    if (vd_curly%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_curly%req(i))
    if (vd_curlz%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_curlz%req(i))
    if (vd_velox%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_velox%req(i))
    if (vd_veloy%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_veloy%req(i))
    if (vd_veloz%req(i)   /= VAL_NOT_ASSIGNED) call wait_req(vd_veloz%req(i))
  enddo

end subroutine wait_vol_recv


subroutine write_vol_data(it_io, val_type_mov,nglob_par_io_offset)
  use io_server
  use specfem_par
  use phdf5_utils
  implicit none

  integer, intent(in) :: it_io
  logical, dimension(5), intent(inout) :: val_type_mov
  integer :: i,j, num_max_type=5, comm, info, id_loc
  integer, dimension(0:NIONOD-1), intent(in) :: nglob_par_io_offset   ! offset inter io node info
  integer :: io_offset
  ! make output file
  character(len=10) :: tempstr
  character(len=64) :: dset_name
  character(len=64) :: group_name

  type(h5io) :: h5

  ! calculate total number of nglob
  !nglob_all = sum(nglob_par_io_offset(:))
  if (NtoM) then
    io_offset = 0 !
  else
    io_offset = sum(nglob_par_io_offset(0:myrank-1))
  endif

  ! initialization of h5 file
  h5 = h5io()
  call h5_init(h5, fname_h5_data_vol)

  ! create time group in h5
  write(tempstr, "(i6.6)") it_io
  group_name = "it_"//tempstr

  ! create dataset
  if (NtoM) then ! for NtoM strategy
      ! create it group
      call h5_open_file(h5)
      ! check if group_name exists
      call h5_open_or_create_group(h5, group_name)

      ! loop for each value type
      do j=1, num_max_type
        if (val_type_mov(j) .eqv. .true.) then
          if (j==1) then
            dset_name = "pressure"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
          elseif (j==2) then
            dset_name = "div_glob"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
          elseif (j==3) then
            dset_name = "div"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
          elseif (j==4) then
            dset_name = "curl_x"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
            dset_name = "curl_y"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
            dset_name = "curl_z"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
          elseif (j==5) then
            dset_name = "velo_x"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
            dset_name = "velo_y"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
            dset_name = "velo_z"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_this_io/), 1, CUSTOM_REAL)
          endif
        endif
      enddo

      call h5_close_group(h5)
      call h5_close_file(h5)

  else

    if (myrank==0) then
      ! create it group
      call h5_open_file(h5)
      ! check if group_name exists
      call h5_open_or_create_group(h5, group_name)

      ! loop for each value type
      do j=1, num_max_type
        if (val_type_mov(j) .eqv. .true.) then
          if (j==1) then
            dset_name = "pressure"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
          elseif (j==2) then
            dset_name = "div_glob"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
          elseif (j==3) then
            dset_name = "div"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
          elseif (j==4) then
            dset_name = "curl_x"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
            dset_name = "curl_y"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
            dset_name = "curl_z"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
          elseif (j==5) then
            dset_name = "velo_x"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
            dset_name = "velo_y"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
            dset_name = "velo_z"
            call h5_create_dataset_gen_in_group(h5, dset_name, (/nglob_all_server/), 1, CUSTOM_REAL)
          endif
        endif
      enddo

      call h5_close_group(h5)
      call h5_close_file(h5)
    endif ! if myrank==0
    call synchronize_all()

  endif ! Nto1

  if (NtoM) then
    call h5_open_file(h5)
  else
    call h5_open_file_p(h5)
  endif
  call h5_open_group(h5,group_name)

  ! loop for each value type
  do j=1, num_max_type
    if (val_type_mov(j) .eqv. .true.) then
      if (j==1) then
        dset_name = "pressure"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_pres%d1darr, (/io_offset/), if_collect_server)
      elseif (j==2) then
        dset_name = "div_glob"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_divglob%d1darr, (/io_offset/), if_collect_server)
      elseif (j==3) then
        dset_name = "div"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_div%d1darr, (/io_offset/), if_collect_server)
      elseif (j==4) then
        dset_name = "curl_x"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_curlx%d1darr, (/io_offset/), if_collect_server)
        dset_name = "curl_y"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_curly%d1darr, (/io_offset/), if_collect_server)
        dset_name = "curl_z"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_curlz%d1darr, (/io_offset/), if_collect_server)
      elseif (j==5) then
        dset_name = "velo_x"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_velox%d1darr, (/io_offset/), if_collect_server)
        dset_name = "velo_y"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_veloy%d1darr, (/io_offset/), if_collect_server)
        dset_name = "velo_z"
        call h5_write_dataset_1d_r_collect_hyperslab_in_group(h5, dset_name, &
          vd_veloz%d1darr, (/io_offset/), if_collect_server)
      endif

    endif
  enddo

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine write_vol_data


!
! shakemap
!
subroutine shakemap_init(nfaces_perproc, surface_offset)
  use io_server
  use phdf5_utils
  use specfem_par

  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: nfaces_perproc, surface_offset
  integer                                   :: ier, len_array_aug
  character(len=64)                         :: dset_name
  character(len=64)                         :: group_name

  type(h5io) :: h5
  h5 = h5io()

  fname_h5_data_shake = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/shakemap.h5"

  ! initialization of h5 file
  call h5_init(h5, fname_h5_data_shake)
  ! create a hdf5 file
  call h5_create_file(h5)

  ! information for computer node
  allocate(shake_ux(size_surf_array),stat=ier)
  allocate(shake_uy(size_surf_array),stat=ier)
  allocate(shake_uz(size_surf_array),stat=ier)
  if(USE_HIGHRES_FOR_MOVIES) then
    len_array_aug=size(surf_x_aug)
    allocate(shake_ux_aug(len_array_aug),stat=ier)
    allocate(shake_uy_aug(len_array_aug),stat=ier)
    allocate(shake_uz_aug(len_array_aug),stat=ier)
  endif

  ! write xyz coords in h5
  group_name = "surf_coord"
  call h5_create_group(h5, group_name)
  call h5_open_group(h5, group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "x"
    call h5_write_dataset_1d_d(h5, dset_name, surf_x)
    call h5_close_dataset(h5)
    dset_name = "y"
    call h5_write_dataset_1d_d(h5, dset_name, surf_y)
    call h5_close_dataset(h5)
    dset_name = "z"
    call h5_write_dataset_1d_d(h5, dset_name, surf_z)
    call h5_close_dataset(h5)
  else
    dset_name = "x"
    call h5_write_dataset_1d_d(h5, dset_name, surf_x_aug)
    call h5_close_dataset(h5)
    dset_name = "y"
    call h5_write_dataset_1d_d(h5, dset_name, surf_y_aug)
    call h5_close_dataset(h5)
    dset_name = "z"
    call h5_write_dataset_1d_d(h5, dset_name, surf_z_aug)
    call h5_close_dataset(h5)
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine shakemap_init


subroutine recv_shake_data(status, nfaces_perproc, surface_offset)
  use io_server
  use my_mpi
  use specfem_par
  implicit none

  integer, dimension(0:NPROC-1), intent(in)         :: nfaces_perproc, surface_offset
  integer                                           :: ier, sender, tag, i
  integer, intent(in)                               :: status(MPI_STATUS_SIZE)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_array

  sender = status(MPI_SOURCE)
  tag    = status(MPI_TAG)

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

  deallocate(temp_array, stat=ier)

end subroutine recv_shake_data


subroutine write_shake_io()
  use specfem_par
  use io_server
  use phdf5_utils

  implicit none

  character(len=64) :: dset_name
  character(len=64) :: group_name
  type(h5io)        :: h5
  h5 = h5io()

  ! continue opening hdf5 file till the end of write process
  call h5_init(h5, fname_h5_data_shake)
  call h5_open_file(h5)

  ! create a group for each io step
  group_name = "shakemap"
  call h5_create_group(h5, group_name)
  call h5_open_group(h5, group_name)
  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "shakemap_ux"
    call h5_write_dataset_1d_d(h5, dset_name, shake_ux)
    call h5_close_dataset(h5)
    dset_name = "shakemap_uy"
    call h5_write_dataset_1d_d(h5, dset_name, shake_uy)
    call h5_close_dataset(h5)
    dset_name = "shakemap_uz"
    call h5_write_dataset_1d_d(h5, dset_name, shake_uz)
    call h5_close_dataset(h5)
  else
    dset_name = "shakemap_ux"
    call recompose_for_hires(shake_ux, shake_ux_aug)
    call h5_write_dataset_1d_d(h5, dset_name, shake_ux_aug)
    call h5_close_dataset(h5)
    dset_name = "shakemap_uy"
     call recompose_for_hires(shake_uy, shake_uy_aug)
    call h5_write_dataset_1d_d(h5, dset_name, shake_uy_aug)
    call h5_close_dataset(h5)
    dset_name = "shakemap_uz"
    call recompose_for_hires(shake_uz, shake_uz_aug)
    call h5_write_dataset_1d_d(h5, dset_name, shake_uz_aug)
    call h5_close_dataset(h5)
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine write_shake_io

!
! surface movie
!

subroutine surf_mov_init(nfaces_perproc, surface_offset)
  use io_server
  use phdf5_utils
  use specfem_par

  implicit none

  integer, dimension(0:NPROC-1), intent(in)         :: nfaces_perproc, surface_offset
  integer                                           :: ier, nfaces_actual, nfaces_aug=(NGLLX-1)*(NGLLY-1),nnodes_per_face_aug=4
  integer                                           :: len_array_aug
  character(len=64)                                 :: dset_name
  character(len=64)                                 :: group_name
  integer, dimension(1)                             :: tmp_1d_arr

  type(h5io) :: h5
  h5 = h5io()

  fname_h5_data_surf = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_surface.h5"

  ! initialization of h5 file
  call h5_init(h5, fname_h5_data_surf)
  ! create a hdf5 file
  call h5_create_file(h5)

  ! information for computer node
  ! get nfaces_perproc_surface
  call recv_i_inter(nfaces_perproc,NPROC,0,io_tag_surface_nfaces)
  ! get faces_surface_offset
  call recv_i_inter(surface_offset,NPROC,0,io_tag_surface_offset)

  ! get xyz coordinates
  call recv_i_inter(tmp_1d_arr, 1, 0, io_tag_surface_coord_len)
  size_surf_array = tmp_1d_arr(1)
  !print *, "size surf array received: ", size_surf_array
  allocate(surf_x(size_surf_array),stat=ier)
  allocate(surf_y(size_surf_array),stat=ier)
  allocate(surf_z(size_surf_array),stat=ier)
  allocate(surf_ux(size_surf_array),stat=ier)
  allocate(surf_uy(size_surf_array),stat=ier)
  allocate(surf_uz(size_surf_array),stat=ier)

  if (USE_HIGHRES_FOR_MOVIES) then
    nfaces_actual = size_surf_array/(NGLLX*NGLLY)
    len_array_aug = nfaces_actual*nfaces_aug*nnodes_per_face_aug
    allocate(surf_x_aug(len_array_aug),stat=ier)
    allocate(surf_y_aug(len_array_aug),stat=ier)
    allocate(surf_z_aug(len_array_aug),stat=ier)
    allocate(surf_ux_aug(len_array_aug),stat=ier)
    allocate(surf_uy_aug(len_array_aug),stat=ier)
    allocate(surf_uz_aug(len_array_aug),stat=ier)
  endif

  ! x
  call recvv_cr_inter(surf_x, size_surf_array, 0, io_tag_surface_x)
  ! y
  call recvv_cr_inter(surf_y, size_surf_array, 0, io_tag_surface_y)
  ! z
  call recvv_cr_inter(surf_z, size_surf_array, 0, io_tag_surface_z)

  ! write xyz coords in h5
  group_name = "surf_coord"
  call h5_create_group(h5, group_name)
  call h5_open_group(h5, group_name)

  ! low resolution output
  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "x"
    call h5_write_dataset_1d_d(h5, dset_name, surf_x)
    call h5_close_dataset(h5)
    dset_name = "y"
    call h5_write_dataset_1d_d(h5, dset_name, surf_y)
    call h5_close_dataset(h5)
    dset_name = "z"
    call h5_write_dataset_1d_d(h5, dset_name, surf_z)
    call h5_close_dataset(h5)
  ! high resolution output
  else
    ! nfaces*25nodes => n*16faces*4
    dset_name = "x"
    call recompose_for_hires(surf_x,surf_x_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_x_aug)
    call h5_close_dataset(h5)
    dset_name = "y"
    call recompose_for_hires(surf_y,surf_y_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_y_aug)
    call h5_close_dataset(h5)
    dset_name = "z"
    call recompose_for_hires(surf_z,surf_z_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_z_aug)
    call h5_close_dataset(h5)
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)

end subroutine surf_mov_init


subroutine recv_surf_data(status, nfaces_perproc, surface_offset)
  use io_server
  use my_mpi
  use specfem_par
  implicit none

  integer, dimension(0:NPROC-1), intent(in)         :: nfaces_perproc, surface_offset
  integer                                           :: ier, sender, tag, i
  integer                                           :: msgsize
  integer, intent(in)                               :: status(MPI_STATUS_SIZE)
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: temp_array
  sender = status(MPI_SOURCE)
  tag    = status(MPI_TAG)

  allocate(temp_array(nfaces_perproc(sender)), stat=ier)
  temp_array(:) = 0._CUSTOM_REAL
  call get_size_msg(status, msgsize)

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

  deallocate(temp_array, stat=ier)

end subroutine recv_surf_data


subroutine write_surf_io(it_io)
  use io_server
  use specfem_par
  use phdf5_utils

  implicit none

  integer, intent(in)                               :: it_io
  !integer                                           :: ier, nfaces_actual, nfaces_aug=16,nnodes_per_face_aug=4
  character(len=64)                                 :: dset_name
  character(len=64)                                 :: group_name
  character(len=10)                                 :: tempstr
  type(h5io)                                        :: h5
  h5 = h5io()

  ! continue opening hdf5 file till the end of write process
  call h5_init(h5, fname_h5_data_surf)
  call h5_open_file(h5)
  ! create a group for each io step
  write(tempstr, "(i6.6)") it_io
  group_name = "it_"//tempstr
  call h5_create_group(h5, group_name)
  call h5_open_group(h5, group_name)

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    dset_name = "ux"
    call h5_write_dataset_1d_d(h5, dset_name, surf_ux)
    call h5_close_dataset(h5)
    dset_name = "uy"
    call h5_write_dataset_1d_d(h5, dset_name, surf_uy)
    call h5_close_dataset(h5)
    dset_name = "uz"
    call h5_write_dataset_1d_d(h5, dset_name, surf_uz)
    call h5_close_dataset(h5)

    surf_ux(:) = 0._CUSTOM_REAL
    surf_uy(:) = 0._CUSTOM_REAL
    surf_uz(:) = 0._CUSTOM_REAL
  else
    dset_name = "ux"
    call recompose_for_hires(surf_ux, surf_ux_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_ux_aug)
    call h5_close_dataset(h5)
    dset_name = "uy"
    call recompose_for_hires(surf_uy, surf_uy_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_uy_aug)
    call h5_close_dataset(h5)
    dset_name = "uz"
    call recompose_for_hires(surf_uz, surf_uz_aug)
    call h5_write_dataset_1d_d(h5, dset_name, surf_uz_aug)
    call h5_close_dataset(h5)

    surf_ux_aug(:) = 0._CUSTOM_REAL
    surf_uy_aug(:) = 0._CUSTOM_REAL
    surf_uz_aug(:) = 0._CUSTOM_REAL
  endif

  call h5_close_group(h5)
  call h5_close_file(h5)


end subroutine write_surf_io

!
! seismo
!

subroutine get_receiver_info(islice_num_rec_local)
  use specfem_par
  use my_mpi

  implicit none

  integer                      :: ier, iproc
  integer, dimension(1)        :: nrec_local_temp
  integer, dimension(1)        :: nrec_temp, nlength_seismogram_tmp
  integer,dimension(0:NPROC-1) :: islice_num_rec_local


  call recv_i_inter(nrec_temp, 1, 0, io_tag_num_recv)
  nrec = nrec_temp(1)

  call recv_dp_inter(t0, 1, 0, io_tag_seismo_tzero)

  call recv_i_inter(nlength_seismogram_tmp, 1, 0, io_tag_seismo_length)
  nlength_seismogram = nlength_seismogram_tmp(1)

  do iproc = 0, NPROC-1
    call recv_i_inter(nrec_local_temp, 1, iproc, io_tag_local_rec)
    islice_num_rec_local(iproc) = nrec_local_temp(1)
  enddo

end subroutine get_receiver_info


subroutine allocate_seismo_arrays(islice_num_rec_local)
  use specfem_par
  use io_server

  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local
  integer                                   :: ier, max_num_rec, nstep_temp

  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS >= nlength_seismogram) then
    nstep_temp = nlength_seismogram
  else
    nstep_temp = NTSTEP_BETWEEN_OUTPUT_SEISMOS
  endif

  ! allocate id_rec_globs for storing global id of receivers
  max_num_rec = maxval(islice_num_rec_local)
  allocate(id_rec_globs(max_num_rec,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array id_rec_globs')
  if (ier /= 0) stop 'error allocating array id_rec_globs'
  ! initialize
  id_rec_globs(:,:) = 0

  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    allocate(seismo_disp(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seimo_disp')
    if (ier /= 0) stop 'error allocating array seismo_disp'
    seismo_disp(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    allocate(seismo_velo(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_velo')
    if (ier /= 0) stop 'error allocating array seismo_velo'
    seismo_velo(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    allocate(seismo_acce(NDIM,nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_acce')
    if (ier /= 0) stop 'error allocating array seismo_acce'
    seismo_acce(:,:,:) = 0._CUSTOM_REAL
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    allocate(seismo_pres(nstep_temp,nrec),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_pres')
    if (ier /= 0) stop 'error allocating array seismo_pres'
    seismo_pres(:,:) = 0._CUSTOM_REAL
  endif

end subroutine allocate_seismo_arrays


subroutine deallocate_arrays()
  use specfem_par
  use io_server

  implicit none

  integer :: ier

  ! seismo
  deallocate(id_rec_globs,stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error deallocating array id_rec_globs')
  if (ier /= 0) stop 'error deallocating array id_rec_globs'

  if (allocated(seismo_disp)) then
    deallocate(seismo_disp,stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error deallocating array seimo_disp')
    if (ier /= 0) stop 'error deallocating array seismo_disp'
  endif
  if (allocated(seismo_velo)) then
    deallocate(seismo_velo,stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error deallocating array seismo_velo')
    if (ier /= 0) stop 'error deallocating array seismo_velo'
  endif
  if (allocated(seismo_acce)) then
    deallocate(seismo_acce,stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error deallocating array seismo_acce')
    if (ier /= 0) stop 'error deallocating array seismo_acce'
  endif
  if (allocated(seismo_pres)) then
    deallocate(seismo_pres,stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error deallocating array seismo_pres')
    if (ier /= 0) stop 'error deallocating array seismo_pres'
  endif

  ! surface movie
  deallocate(surf_x,stat=ier)
  deallocate(surf_y,stat=ier)
  deallocate(surf_z,stat=ier)

  deallocate(surf_ux,stat=ier)
  deallocate(surf_uy,stat=ier)
  deallocate(surf_uz,stat=ier)

  ! shake map
  deallocate(shake_ux,stat=ier)
  deallocate(shake_uy,stat=ier)
  deallocate(shake_uz,stat=ier)

  if (USE_HIGHRES_FOR_MOVIES) then
    ! surface movie
    deallocate(surf_x_aug,stat=ier)
    deallocate(surf_y_aug,stat=ier)
    deallocate(surf_z_aug,stat=ier)

    deallocate(surf_ux_aug,stat=ier)
    deallocate(surf_uy_aug,stat=ier)
    deallocate(surf_uz_aug,stat=ier)

    ! shake map
    deallocate(shake_ux_aug,stat=ier)
    deallocate(shake_uy_aug,stat=ier)
    deallocate(shake_uz_aug,stat=ier)
  endif

end subroutine deallocate_arrays


subroutine count_seismo_type()
  use specfem_par
  use io_server

  implicit none

  integer :: n_type = 0

  if (SAVE_SEISMOGRAMS_DISPLACEMENT) n_type = n_type+1
  if (SAVE_SEISMOGRAMS_VELOCITY)     n_type = n_type+1
  if (SAVE_SEISMOGRAMS_ACCELERATION) n_type = n_type+1
  if (SAVE_SEISMOGRAMS_PRESSURE)     n_type = n_type+1

  n_seismo_type = n_type

end subroutine count_seismo_type


subroutine recv_id_rec(islice_num_rec_local)
  use io_server
  use specfem_par
  implicit none

  integer :: sender
  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local

  do sender = 0, NPROC-1
    if (islice_num_rec_local(sender) /= 0) then
      call recv_i_inter(id_rec_globs(:,sender), size(id_rec_globs(:,sender)), sender, io_tag_seismo_ids_rec)
    endif
  enddo
end subroutine recv_id_rec


subroutine recv_seismo_data(status, islice_num_rec_local, rec_count_seismo)
  use my_mpi
  use specfem_par
  use io_server
  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local
  integer, intent(in)                       :: status(MPI_STATUS_SIZE), rec_count_seismo

  integer :: rec_id_glob, sender, nrec_passed, irec_passed, tag, irec, id_rec_glob, ier
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismo_temp
  integer                                               :: msg_size,time_window

  sender       = status(MPI_SOURCE)
  tag          = status(MPI_TAG)
  nrec_passed  = islice_num_rec_local(sender)
  call get_size_msg(status,msg_size)

  if (nrec_passed > 0) then
    ! get vector values i.e. disp, velo, acce
    if (tag /= io_tag_seismo_body_pres) then
      ! allocate temp array size
      time_window = int(msg_size/NDIM/nrec_passed)
      allocate(seismo_temp(NDIM,time_window,nrec_passed),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_temp')
      if (ier /= 0) stop 'error allocating array seismo_temp'
      seismo_temp(:,:,:) = 0.0

      call recvv_cr_inter(seismo_temp, msg_size, sender, tag)

    ! get scalar value i.e. pres
    else
      time_window = int(msg_size/1/nrec_passed)
      allocate(seismo_temp(1,time_window,nrec_passed),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_temp')
      if (ier /= 0) stop 'error allocating array seismo_temp'
      call recvv_cr_inter(seismo_temp, msg_size, sender, tag)
    endif

    ! set local array to the global array
    do irec_passed=1,nrec_passed
      id_rec_glob = id_rec_globs(irec_passed,sender)
      ! disp
      if (tag == io_tag_seismo_body_disp) then
        seismo_disp(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed)
      ! velo
      elseif (tag == io_tag_seismo_body_velo) then
        seismo_velo(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed) ! id_rec
      ! acce
      elseif (tag == io_tag_seismo_body_acce) then
        seismo_acce(:,:,id_rec_glob) = seismo_temp(:,:,irec_passed)
      ! pres
      elseif (tag == io_tag_seismo_body_pres) then
        seismo_pres(:,id_rec_glob) = seismo_temp(1,:,irec_passed)
      endif
    enddo

    ! deallocate temp array
    deallocate(seismo_temp,stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error deallocating array seismo_temp')
    if (ier /= 0) stop 'error allocating dearray seismo_temp'

  endif
end subroutine recv_seismo_data


! counts number of local receivers for each slice
subroutine count_nprocs_with_recs(islice_num_rec_local)
  use my_mpi
  use specfem_par, only: nrec,NPROC
  use io_server

  implicit none

  integer, dimension(0:NPROC-1) :: islice_num_rec_local
  integer                       :: irec, iproc

  do iproc = 0, NPROC-1
    if (islice_num_rec_local(iproc) > 0) &
      n_procs_with_rec = n_procs_with_rec+1
  enddo

end subroutine count_nprocs_with_recs


subroutine do_io_seismogram_init()
  use specfem_par
  use phdf5_utils
  use io_server

  implicit none

  ! hdf5 varianles
  character(len=64) :: fname_h5_base = "seismograms.h5"
  type(h5io)        :: h5

  ! mpi variables
  integer :: info, comm, error

  ! arrays
  integer                                                 :: i, irec
  real(kind=CUSTOM_REAL), dimension(nlength_seismogram)   :: time_array
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable     :: val_array2d
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable   :: val_array3d
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: stations
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: networks
  real(kind=CUSTOM_REAL), dimension(nrec,3)               :: rec_coords

  ! hdf5 utility
  h5 = h5io()
  fname_h5_seismo = trim(OUTPUT_FILES)//fname_h5_base
  ! initialze hdf5
  call h5_init(h5, fname_h5_seismo)

  ! create file
  call h5_create_file(h5)

  ! create time dataset it = 1 ~ NSTEP
  do i = 1, nlength_seismogram
    if (SIMULATION_TYPE == 1) then ! forward simulation ! distinguish between single and double precision for reals
      time_array(i) = real( dble(i-1)*DT - t0 ,kind=CUSTOM_REAL)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint simulation: backward/reconstructed wavefields
      ! distinguish between single and double precision for reals
      ! note: compare time_t with time used for source term
      time_array(i) = real( dble(nlength_seismogram-i)*DT - t0 ,kind=CUSTOM_REAL)
    endif
  enddo

  ! time array
  call h5_write_dataset_1d_d_no_group(h5, "time", time_array)
  call h5_close_dataset(h5)

  ! read out_list_stations.txt generated at locate_receivers.f90:431 here to write in the h5 file.
  open(unit=IOUT_SU,file=trim(OUTPUT_FILES)//'output_list_stations.txt', &
       status='unknown',action='read',iostat=error)
  if (error /= 0) &
    call exit_mpi(myrank,'error opening file '//trim(OUTPUT_FILES)//'output_list_stations.txt')
  ! writes station infos
  do irec=1,nrec
    read(IOUT_SU,*) stations(irec),networks(irec), rec_coords(irec, 1), rec_coords(irec, 2), rec_coords(irec, 3)
  enddo
  ! closes output file
  close(IOUT_SU)

  ! coordination
  call h5_write_dataset_2d_r_no_group(h5, "coords", rec_coords)
  call h5_close_dataset(h5)

  ! station name
  call h5_write_dataset_1d_c_no_group(h5, "station", stations)
  call h5_close_dataset(h5)

  ! network name
  call h5_write_dataset_1d_c_no_group(h5, "network", networks)
  call h5_close_dataset(h5)

  ! prepare datasets for physical values
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    allocate(val_array3d(NDIM,nlength_seismogram,nrec),stat=error)
    call h5_create_dataset_gen(h5, "disp", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    allocate(val_array3d(NDIM,nlength_seismogram,nrec),stat=error)
    call h5_create_dataset_gen(h5, "velo", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    allocate(val_array3d(NDIM,nlength_seismogram,nrec),stat=error)
    call h5_create_dataset_gen(h5, "acce", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    allocate(val_array2d(nlength_seismogram,nrec),stat=error)
    call h5_create_dataset_gen(h5, "pres", shape(val_array2d), 2, CUSTOM_REAL)
    deallocate(val_array2d)
  endif

  call h5_close_file(h5)

end subroutine do_io_seismogram_init


subroutine write_seismograms_io(it_offset)
  use specfem_par
  use io_server
  use phdf5_utils

  implicit none

  integer, intent(in) :: it_offset
  character(len=4) component
  integer :: t_upper

  ! hdf5 vals
  type(h5io) :: h5

  ! initialze hdf5
  call h5_init(h5, fname_h5_seismo)
  call h5_open_file(h5)

  ! check if the array length to be written > total timestep
  if (it_offset+NTSTEP_BETWEEN_OUTPUT_SEISMOS > nlength_seismogram) then
    t_upper = nlength_seismogram - it_offset
  else
    t_upper = NTSTEP_BETWEEN_OUTPUT_SEISMOS
  endif

  ! writes out this seismogram
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    component = 'disp'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_disp(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    component = 'velo'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_velo(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    component = 'acce'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_acce(:,1:t_upper,:), (/0, it_offset, 0/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    component = 'pres'
    call h5_write_dataset_2d_r_collect_hyperslab(h5, component, seismo_pres(1:t_upper,:), (/it_offset, 0/), .false.)
  endif

  call h5_close_file(h5)

end subroutine write_seismograms_io

!
! xdmf output routines
!
subroutine write_xdmf_surface_header(num_nodes, pos_line)
  use specfem_par
  use io_server
  implicit none
  integer :: num_elm, num_nodes
  integer, intent(out) :: pos_line ! useful for no ioserver mode

  num_elm = num_nodes/4

  ! writeout xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"

  open(unit=xdmf_surf, file=fname_xdmf_surf, recl=256)
  write(xdmf_surf,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_surf,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_surf,*) '<Xdmf Version="3.0">'
  write(xdmf_surf,*) '  <Domain Name="mesh">'
  write(xdmf_surf,*) '    <Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_surf,*) '    <Geometry GeometryType="X_Y_Z">'
  write(xdmf_surf,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        ./DATABASES_MPI/movie_surface.h5:/surf_coord/x'
  write(xdmf_surf,*) '      </DataItem>'
  write(xdmf_surf,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        ./DATABASES_MPI/movie_surface.h5:/surf_coord/y'
  write(xdmf_surf,*) '      </DataItem>'
  write(xdmf_surf,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '        ./DATABASES_MPI/movie_surface.h5:/surf_coord/z'
  write(xdmf_surf,*) '      </DataItem>'
  write(xdmf_surf,*) '    </Geometry>'

  write(xdmf_surf,*) '    <Grid Name="fensap" GridType="Collection" CollectionType="Temporal" >'
  write(xdmf_surf,*) '    </Grid>'

  write(xdmf_surf,*) '  </Domain>'
  write(xdmf_surf,*) '</Xdmf>'
  ! 20 lines

  ! position where the additional data will be inserted
  surf_xdmf_pos = 17

  close(xdmf_surf)

end subroutine write_xdmf_surface_header


subroutine write_xdmf_surface_body(it_io, num_nodes, pos_line_store)
  use specfem_par
  use io_server

  implicit none

  integer, intent(in)    :: it_io
  integer                :: i
  integer, intent(in)    :: num_nodes
  integer, intent(inout) :: pos_line_store ! useful for no ioserver mode

  character(len=20)  :: it_str
  character(len=20)  :: temp_str

  ! for no ioserver mode
  if (surf_xdmf_pos == 0) surf_xdmf_pos = pos_line_store

  ! redefinition for no ioserver case
  if (fname_xdmf_surf == '' ) fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"

  ! open xdmf file
  open(unit=xdmf_surf, file=fname_xdmf_surf, recl=256)

  ! skip lines till the position where we want to write new information
  do i = 1, surf_xdmf_pos
    read(xdmf_surf, *)
  enddo

  write(it_str, "(i6.6)") it_io
  write(xdmf_surf,*) '<Grid Name="surf_mov" GridType="Uniform">'
  write(xdmf_surf,*) '  <Time Value="'//trim(r2c(sngl((it_io-1)*DT-t0)))//'" />'
  write(xdmf_surf,*) '  <Topology Reference="/Xdmf/Domain/Topology" />'
  write(xdmf_surf,*) '  <Geometry Reference="/Xdmf/Domain/Geometry" />'
  write(xdmf_surf,*) '  <Attribute Name="ux" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      ./DATABASES_MPI/movie_surface.h5:/it_'//trim(it_str)//'/ux'
  write(xdmf_surf,*) '    </DataItem>'
  write(xdmf_surf,*) '  </Attribute>'
  write(xdmf_surf,*) '  <Attribute Name="uy" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      ./DATABASES_MPI/movie_surface.h5:/it_'//trim(it_str)//'/uy'
  write(xdmf_surf,*) '    </DataItem>'
  write(xdmf_surf,*) '  </Attribute>'
  write(xdmf_surf,*) '  <Attribute Name="uz" AttributeType="Scalar" Center="Node">'
  write(xdmf_surf,*) '     <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_surf,*) '      ./DATABASES_MPI/movie_surface.h5:/it_'//trim(it_str)//'/uz'
  write(xdmf_surf,*) '     </DataItem>'
  write(xdmf_surf,*) '  </Attribute>'
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Grid>'
  write(xdmf_surf,*) '</Domain>'
  write(xdmf_surf,*) '</Xdmf>'
  !
  surf_xdmf_pos = surf_xdmf_pos+20

  close(xdmf_surf)

  pos_line_store = surf_xdmf_pos

end subroutine write_xdmf_surface_body


subroutine write_xdmf_shakemap(num_nodes)
  use specfem_par
  use io_server
  implicit none
  integer :: num_elm, num_nodes

  num_elm = num_nodes/4

  ! writeout xdmf file for surface movie
  fname_xdmf_shake = trim(OUTPUT_FILES)//"/shakemap.xmf"

  open(unit=xdmf_shake, file=fname_xdmf_shake, recl=256)

  write(xdmf_shake,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_shake,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_shake,*) '<Xdmf Version="3.0">'
  write(xdmf_shake,*) '  <Domain Name="shakemap">'
  write(xdmf_shake,*) '  <Grid>'
  write(xdmf_shake,*) '    <Topology Name="topo" TopologyType="Quadrilateral" NumberOfElements="'//trim(i2c(num_elm))//'"/>'
  write(xdmf_shake,*) '    <Geometry GeometryType="X_Y_Z">'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="' &
                                                    //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/surf_coord/x'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                    //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/surf_coord/y'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                   //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/surf_coord/z'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '    </Geometry>'
  write(xdmf_shake,*) '    <Attribute Name="shake_ux" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                  //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/shakemap/shakemap_ux'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '    </Attribute>'
  write(xdmf_shake,*) '    <Attribute Name="shake_uy" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                 //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/shakemap/shakemap_uy'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '    </Attribute>'
  write(xdmf_shake,*) '    <Attribute Name="shake_uz" AttributeType="Scalar" Center="Node">'
  write(xdmf_shake,*) '      <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(i2c(num_nodes))//'">'
  write(xdmf_shake,*) '        ./DATABASES_MPI/shakemap.h5:/shakemap/shakemap_uz'
  write(xdmf_shake,*) '      </DataItem>'
  write(xdmf_shake,*) '    </Attribute>'
  write(xdmf_shake,*) '  </Grid>'
  write(xdmf_shake,*) '  </Domain>'
  write(xdmf_shake,*) '</Xdmf>'


  close(xdmf_shake)

end subroutine write_xdmf_shakemap


subroutine write_xdmf_vol(val_type_mov)
  use specfem_par
  use io_server
  implicit none
  logical, dimension(5), intent(inout) :: val_type_mov
  logical :: pressure_io, divglob_io, div_io, veloc_io, curl_io
  character(len=20)                         :: proc_str, it_str,nelm, nglo
  integer                                   :: iiout, nout, i, ii

  pressure_io = val_type_mov(1)
  divglob_io  = val_type_mov(2)
  div_io      = val_type_mov(3)
  curl_io     = val_type_mov(4)
  veloc_io    = val_type_mov(5)

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES)//"/movie_volume.xmf"

  open(unit=xdmf_vol, file=fname_xdmf_vol, recl=256)

  ! definition of topology and geometry
  ! refer only control nodes (8 or 27) as a coarse output
  ! data array need to be extracted from full data array on gll points
  write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
  write(xdmf_vol,*) '<Domain name="mesh">'
  ! loop for writing information of mesh partitions
  nelm=i2c(nspec_all_server*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
  nglo=i2c(nglob_all_server)

  write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm)//'">'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                         //trim(nelm)//' 9">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/elm_conn'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '</Topology>'
  write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/x'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/y'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                      //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
  write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume.h5:/mesh/z'
  write(xdmf_vol,*) '    </DataItem>'
  write(xdmf_vol,*) '</Geometry>'
  write(xdmf_vol,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'

  do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ii = i*NTSTEP_BETWEEN_FRAMES
    write(it_str, "(i6.6)") ii


    write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Uniform">'
    write(xdmf_vol,*) '  <Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'
    write(xdmf_vol,*) '  <Topology Reference="/Xdmf/Domain/Topology" />'
    write(xdmf_vol,*) '  <Geometry Reference="/Xdmf/Domain/Geometry" />'

    if (pressure_io) then
      write(xdmf_vol,*) '  <Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/pressure'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (divglob_io) then
      write(xdmf_vol,*) '  <Attribute Name="div_glob" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/div_glob'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (div_io) then
      write(xdmf_vol,*) '  <Attribute Name="div" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/div'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (veloc_io) then
      write(xdmf_vol,*) '  <Attribute Name="velo_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_x'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="velo_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_y'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="velo_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                         //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/velo_z'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    if (curl_io) then
      write(xdmf_vol,*) '  <Attribute Name="curl_x" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_x'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="curl_y" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_y'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
      write(xdmf_vol,*) '  <Attribute Name="curl_z" AttributeType="Scalar" Center="Node">'
      write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
      write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume.h5:/it_'//trim(it_str)//'/curl_z'
      write(xdmf_vol,*) '    </DataItem>'
      write(xdmf_vol,*) '  </Attribute>'
    endif

    write(xdmf_vol,*) '</Grid>'
  enddo

  write(xdmf_vol,*) '</Grid>'
  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

end subroutine write_xdmf_vol


subroutine write_xdmf_vol_mton(val_type_mov,nglob_par_io_offset,nelm_par_io_offset)
  use specfem_par
  use io_server
  implicit none

  logical, dimension(5), intent(inout)        :: val_type_mov
  integer, dimension(0:NIONOD-1), intent(in)  :: nglob_par_io_offset, nelm_par_io_offset
  logical :: pressure_io, divglob_io, div_io, veloc_io, curl_io
  character(len=20)                           :: proc_str, it_str,nelm, nglo
  integer                                     :: iiout, nout, i, ii, iproc

  pressure_io = val_type_mov(1)
  divglob_io  = val_type_mov(2)
  div_io      = val_type_mov(3)
  curl_io     = val_type_mov(4)
  veloc_io    = val_type_mov(5)

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES)//"/movie_volume.xmf"

  open(unit=xdmf_vol, file=fname_xdmf_vol, recl=256)

  ! definition of topology and geometry
  ! refer only control nodes (8 or 27) as a coarse output
  ! data array need to be extracted from full data array on gll points
  write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
  write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
  write(xdmf_vol,*) '<Domain>'
  write(xdmf_vol,*) '    <!-- mesh info -->'
  write(xdmf_vol,*) '    <Grid Name="mesh" GridType="Collection"  CollectionType="Spatial">'

  ! loop for writing information of mesh partitions
  do iproc=0,NIONOD-1

    ! loop for writing information of mesh partitions
    nelm=i2c(nelm_par_io_offset(iproc)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
    nglo=i2c(nglob_par_io_offset(iproc))
    write(proc_str, "(i6.6)") iproc

    write(xdmf_vol,*) '<Grid Name="mesh_'//trim(proc_str)//'">'
    write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm)//'">'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                           //trim(nelm)//' 9">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/mesh/elm_conn'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '</Topology>'
    write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/mesh/x'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/mesh/y'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/mesh/z'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '</Geometry>'
    write(xdmf_vol,*) '</Grid>'


  enddo ! loop for writing mesh info
  write(xdmf_vol,*) '</Grid>'


  ! write the data by
  ! loop for iteration
  !   loop for IONODES

  write(xdmf_vol,*) '<Grid Name="time_col" GridType="Collection" CollectionType="Temporal">'
  do i = 1, int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ii = i*NTSTEP_BETWEEN_FRAMES
    write(it_str, "(i6.6)") ii

    write(xdmf_vol,*) '<Grid Name="vol_mov" GridType="Collection" CollectionType="Spatial">'
    write(xdmf_vol,*) '  <Time Value="'//trim(r2c(sngl((ii-1)*DT-t0)))//'" />'


    do iproc=0, NIONOD-1
      write(proc_str, "(i6.6)") iproc
      nglo=i2c(nglob_par_io_offset(iproc))

      write(xdmf_vol, *)  '<Grid Name="data_'//trim(proc_str)//'" Type="Uniform">'
      write(xdmf_vol, *)  '    <Topology Reference="/Xdmf/Domain/Grid[@Name=''mesh'']/Grid[@Name=''mesh_'&
                                        //trim(proc_str)//''']/Topology" />'
      write(xdmf_vol, *)  '    <Geometry Reference="/Xdmf/Domain/Grid[@Name=''mesh'']/Grid[@Name=''mesh_'&
                                        //trim(proc_str)//''']/Geometry" />'

      if (pressure_io) then
        write(xdmf_vol,*) '  <Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/pressure'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
      endif

      if (divglob_io) then
        write(xdmf_vol,*) '  <Attribute Name="div_glob" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/div_glob'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
      endif

      if (div_io) then
        write(xdmf_vol,*) '  <Attribute Name="div" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/div'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
      endif

      if (veloc_io) then
        write(xdmf_vol,*) '  <Attribute Name="velo_x" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/velo_x'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
        write(xdmf_vol,*) '  <Attribute Name="velo_y" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/velo_y'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
        write(xdmf_vol,*) '  <Attribute Name="velo_z" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                           //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/velo_z'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
      endif

      if (curl_io) then
        write(xdmf_vol,*) '  <Attribute Name="curl_x" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/curl_x'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
        write(xdmf_vol,*) '  <Attribute Name="curl_y" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/curl_y'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
        write(xdmf_vol,*) '  <Attribute Name="curl_z" AttributeType="Scalar" Center="Node">'
        write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                       //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
        write(xdmf_vol,*) '      ./DATABASES_MPI/movie_volume_'//trim(proc_str)//'.h5:/it_'//trim(it_str)//'/curl_z'
        write(xdmf_vol,*) '    </DataItem>'
        write(xdmf_vol,*) '  </Attribute>'
      endif

      write(xdmf_vol,*) '</Grid>'

    enddo  ! loop proc

    write(xdmf_vol,*) '</Grid>'

  enddo ! loop time

  write(xdmf_vol,*) '</Grid>'
  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

end subroutine write_xdmf_vol_mton




subroutine pass_info_to_io()
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use constants, only: dest_ionod

  implicit none

  integer ::  n_msg_vol_each_proc = 0,irec,irec_local,i_ionod,ier
  integer,dimension(nrec_local) :: tmp_irec
  integer, dimension(9,NSPEC_AB*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)) :: elm_conn_loc

  ! initialization of io node from compute node side


  ! send the receiver information to io node for the outputs of seismo signals
  ! seismo data is out only from rank=0 of io nodes
  if (myrank == 0) then
    ! send nrec and nrec_local
    call send_i_inter((/nrec/), 1, 0, io_tag_num_recv)
    ! send t0
    call send_dp_inter((/t0/), 1, 0, io_tag_seismo_tzero)
    ! send nlength_seismogram
    call send_i_inter((/nlength_seismogram/), 1, 0, io_tag_seismo_length)
  endif
  ! send the number of local receiver to the io node
  call send_i_inter((/nrec_local/), 1, 0, io_tag_local_rec)

  ! send global id of stations
  if (nrec_local > 0) then
    ! send global ids of local receivers (integer array)
    do irec_local = 1,nrec_local
      ! get global number of that receiver
      irec = number_receiver_global(irec_local)
      tmp_irec(irec_local) = irec
    enddo

    call send_i_inter(tmp_irec,nrec_local,0,io_tag_seismo_ids_rec)
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
    endif
    if (MOVIE_VOLUME) then
      ! count number of messges for volume movie
      if (ACOUSTIC_SIMULATION .and. .not. ELASTIC_SIMULATION .and. .not. POROELASTIC_SIMULATION) then
          n_msg_vol_each_proc = n_msg_vol_each_proc+1 ! pressure
      endif
      if (ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
        if (ELASTIC_SIMULATION) n_msg_vol_each_proc = n_msg_vol_each_proc+1 ! div_glob
        n_msg_vol_each_proc = n_msg_vol_each_proc+4 ! div, curl_x, curl_y, curl_z
      endif
      if (ACOUSTIC_SIMULATION .or. ELASTIC_SIMULATION .or. POROELASTIC_SIMULATION) then
        n_msg_vol_each_proc = n_msg_vol_each_proc+3 ! velocity_x,velocity_y,velocity_z
      endif
      ! send the number of messages
      do i_ionod=0,NIONOD-1
        call send_i_inter((/n_msg_vol_each_proc/),1,i_ionod,io_tag_vol_nmsg)
      enddo
    endif
  endif ! end if myrank == 0

  if (MOVIE_VOLUME) then
    ! send the compute node list to io node
    call send_i_inter((/0/),1,dest_ionod,io_tag_vol_sendlist)

    do i_ionod=0, NIONOD-1
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
    ! gather nglobs of each proc
    call gather_all_all_singlei((/NGLOB_AB/),nglob_offset,NPROC)
    ! create connectivity dataset (it's node base but element base)
    call get_conn_for_movie(NSPEC_AB, elm_conn_loc, 0)!nglob_offset(myrank))
    ! send elm_conn_loc
    call send_i_inter(elm_conn_loc,size(elm_conn_loc),dest_ionod, io_tag_vol_elmconn)
    ! send xstore, ystore, zstore
    call sendv_cr_inter(xstore, size(xstore), dest_ionod, io_tag_vol_nodex)
    call sendv_cr_inter(ystore, size(ystore), dest_ionod, io_tag_vol_nodey)
    call sendv_cr_inter(zstore, size(zstore), dest_ionod, io_tag_vol_nodez)

    deallocate(nglob_offset, stat=ier)
  endif

end subroutine pass_info_to_io
