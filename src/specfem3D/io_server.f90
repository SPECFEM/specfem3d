! module for storing info. concering io
module io_server
  use specfem_par, only: CUSTOM_REAL, NPROC

  implicit none

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
    integer :: req=VAL_NOT_ASSIGNED
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: d1darr
  end type vol_data_dump

  ! dumps for volume data
  type(vol_data_dump), dimension(:), allocatable :: vd_pres, vd_divglob, vd_div,  &
                                                    vd_curlx, vd_curly, vd_curlz, &
                                                    vd_velox, vd_veloy, vd_veloz

  integer, dimension(:), allocatable:: id_proc_glob2loc, id_proc_loc2glob
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
    integer :: i,j,k,c,factor_face_aug=16,npoint_per_face=25,npoint_corner=4

    nfaces_actual=size(arr_in)/(NGLLX*NGLLY)

    c=1
    do i=0, nfaces_actual-1
      do j=0,3 !y
        do k=0,3 !x
          arr_out(c  +j*4*4+k*4)=arr_in(i*npoint_per_face+1+k+j*NGLLX)
          arr_out(c+1+j*4*4+k*4)=arr_in(i*npoint_per_face+2+k+j*NGLLX)
          arr_out(c+2+j*4*4+k*4)=arr_in(i*npoint_per_face+7+k+j*NGLLX)
          arr_out(c+3+j*4*4+k*4)=arr_in(i*npoint_per_face+6+k+j*NGLLX)
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

  ! vars shakemap
  integer :: rec_count_shake=0, n_recv_msg_shake=0, shake_out_count=0,max_shake_out=0

  ! vars volumne movie
  integer                       :: rec_count_vol=0, n_recv_msg_vol=0, vol_out_count=0, max_vol_out=0
  integer, dimension(0:NPROC-1) :: nelm_par_proc, nglob_par_proc ! storing the number of elements and gll nodes
  logical, dimension(5)         :: val_type_mov ! true if movie file will be created, (pressure, div_glob, div, curlxyz, velocity_xyz)

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
      if (MOVIE_SURFACE) then
        n_recv_msg_surf = n_msg_surf_each_proc*NPROC
        print *, "surf move init done"
        call write_xdmf_surface_header()

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
    call movie_volume_init(nelm_par_proc,nglob_par_proc)
    print *, "movie volume init done"
    n_recv_msg_vol = n_msg_vol_each_proc*nproc_io
    max_vol_out    = int(NSTEP/NTSTEP_BETWEEN_FRAMES)

    ! initialize flags for the value types to be written out
    val_type_mov(:) = .false.

    ! allocate dumping arrays
    allocate(vd_pres(0:nproc_io-1),vd_divglob(0:nproc_io-1),   vd_div(0:nproc_io-1), &
            vd_curlx(0:nproc_io-1),  vd_curly(0:nproc_io-1), vd_curlz(0:nproc_io-1), &
            vd_velox(0:nproc_io-1),  vd_veloy(0:nproc_io-1), vd_veloz(0:nproc_io-1))

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
!    print *,                 "msg: " , status(MPI_TAG) , " send/recv rank: ", status(MPI_SOURCE), "/",myrank, &
!              "  counters, seismo: " , rec_count_seismo, "/"      , n_recv_msg_seismo,  &
!                          ", surf: " , rec_count_surf  , "/"      , n_recv_msg_surf,    &
!                          ", shake: ", rec_count_shake , "/"      , n_recv_msg_shake,   &
!                          ", vol: "  , rec_count_vol   , "/"      , n_recv_msg_vol

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
        call recv_vol_data(status,rec_count_vol,it_io, val_type_mov)
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
      call write_xdmf_surface_body(it_io)

      print *, "surface write done at it = ", it_io
    endif

    ! write volume movie
    if (MOVIE_VOLUME .and. rec_count_vol  == n_recv_msg_vol) then
      ! wait all vol data are reached
      call wait_vol_recv()

      ! write dumped vol data
      call write_vol_data(it_io,val_type_mov)

      rec_count_vol = 0 ! reset counter
      vol_out_count = vol_out_count+1
      if (vol_out_count==1) then
        ! create xdmf header file
        if (myrank == 0) then
          call write_xdmf_vol_header(nelm_par_proc,nglob_par_proc)
        endif
      endif

      if (myrank == 0) then
        call write_xdmf_vol_body_header(it_io)
        call write_xdmf_vol_body(it_io, nelm_par_proc, nglob_par_proc, val_type_mov)
        call write_xdmf_vol_body_close()
      endif
      print *, "volume write done at it = ", it_io

    endif

    ! write shakemap
    if (CREATE_SHAKEMAP .and. rec_count_shake == n_recv_msg_shake .and. myrank==0) then
      call write_shake_io()
      rec_count_shake = 0
      shake_out_count = shake_out_count+1
      ! write out xdmf at each timestep
      call write_xdmf_shakemap()

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
subroutine movie_volume_init(nelm_par_proc,nglob_par_proc)
  use io_server
  use specfem_par
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_acoustic
  use specfem_par_movie
  use phdf5_utils
  implicit none

  integer :: iproc, count=0, id_glob, comm, info

  integer, dimension(0:NPROC-1), intent(inout) :: nelm_par_proc, nglob_par_proc ! storing the number of elements and gll nodes

  ! make output file
  character(len=64) :: group_name
  character(len=64) :: dset_name
  character(len=5)  :: ioidstr
  type(h5io)        :: h5
  h5 = h5io()

  write(ioidstr, "(i5.5)") myrank
  fname_h5_data_vol = LOCAL_PATH(1:len_trim(LOCAL_PATH))//"/movie_volume_"//ioidstr//".h5"

  ! initialization of h5 file
  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)

  ! initialize h5 object
  call h5_init(h5, fname_h5_data_vol)
  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)

  ! create a hdf5 file
  call h5_create_file(h5)
  call h5_close_file(h5)

  ! get n_msg_vol_each_proc
  call recv_i_inter(n_msg_vol_each_proc, 1, 0, io_tag_vol_nmsg)

  ! make an array of local2global relation of sending compute node ids
  allocate(id_proc_loc2glob(0:nproc_io-1))
  ! make an array of global2local relation of sending compute node ids
  allocate(id_proc_glob2loc(0:NPROC-1))
  id_proc_glob2loc(:) = -999999

  do iproc = 0, NPROC-1
    if(mod(iproc,NIONOD)==myrank) then
      id_proc_loc2glob(count) = iproc
      id_proc_glob2loc(iproc) = count
      count                   = count+1
    endif
  enddo

  ! get nspec and nglob from each process
  if (myrank == 0) then
    do iproc = 0, NPROC-1
      call recv_i_inter(nelm_par_proc(iproc),  1, iproc, io_tag_vol_nspec) ! NSPEC_AB
      call recv_i_inter(nglob_par_proc(iproc), 1, iproc, io_tag_vol_nglob) ! NGLOB_AB
    enddo
  endif

end subroutine movie_volume_init


subroutine recv_vol_data(status, rec_count_vol, it_io, val_type_mov)
  use io_server
  use specfem_par
  use my_mpi
  implicit none

  integer, intent(in)                  :: status(MPI_STATUS_SIZE)
  integer, intent(in)                  :: rec_count_vol,it_io
  logical, dimension(5), intent(inout) :: val_type_mov
  integer :: sender_glob, sender_loc, ier, tag, arrsize, msgsize
  logical :: if_aloc

  ! flag for allocating dump arrays only at the initial timestep
  if (it_io == NTSTEP_BETWEEN_FRAMES) then
    if_aloc = .true.
  else
    if_aloc = .false.
  endif

  sender_glob = status(MPI_SOURCE)
  sender_loc  = id_proc_glob2loc(sender_glob)
  tag         = status(MPI_TAG)

  ! get message size
  call get_size_msg(status,msgsize)

  ! receive data and store to dump arrays
  !
  ! here necessary to check if file write has already finished before overwriting dump arrays
  !
  if (tag == io_tag_vol_pres) then
    val_type_mov(1) = .true.
    if(if_aloc) allocate(vd_pres(sender_loc)%d1darr(msgsize),stat=ier)
    vd_pres(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_pres(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_pres(sender_loc)%req)
  elseif (tag == io_tag_vol_divglob) then
    val_type_mov(2) = .true.
  if(if_aloc) allocate(vd_divglob(sender_loc)%d1darr(msgsize),stat=ier)
    vd_divglob(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_divglob(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_divglob(sender_loc)%req)
  elseif (tag == io_tag_vol_div) then
    val_type_mov(3) = .true.
    if(if_aloc) allocate(vd_div(sender_loc)%d1darr(msgsize),stat=ier)
    vd_div(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    ! #BUG: cannot receive the array when using single proc
    call irecvv_cr_inter(vd_div(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_div(sender_loc)%req)
  elseif (tag == io_tag_vol_curlx) then
    val_type_mov(4) = .true.
    if(if_aloc) allocate(vd_curlx(sender_loc)%d1darr(msgsize),stat=ier)
    vd_curlx(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_curlx(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_curlx(sender_loc)%req)
  elseif (tag == io_tag_vol_curly) then
    if(if_aloc) allocate(vd_curly(sender_loc)%d1darr(msgsize),stat=ier)
    vd_curly(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_curly(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_curly(sender_loc)%req)
  elseif (tag == io_tag_vol_curlz) then
    if(if_aloc) allocate(vd_curlz(sender_loc)%d1darr(msgsize),stat=ier)
    vd_curlz(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_curlz(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_curlz(sender_loc)%req)
  elseif (tag == io_tag_vol_velox) then
    val_type_mov(5) = .true.
    if(if_aloc) allocate(vd_velox(sender_loc)%d1darr(msgsize),stat=ier)
    vd_velox(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_velox(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_velox(sender_loc)%req)
  elseif (tag == io_tag_vol_veloy) then
    if(if_aloc) allocate(vd_veloy(sender_loc)%d1darr(msgsize),stat=ier)
    vd_veloy(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_veloy(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_veloy(sender_loc)%req)
  elseif (tag == io_tag_vol_veloz) then
    if(if_aloc) allocate(vd_veloz(sender_loc)%d1darr(msgsize),stat=ier)
    vd_veloz(sender_loc)%d1darr(:) = 0._CUSTOM_REAL
    call irecvv_cr_inter(vd_veloz(sender_loc)%d1darr,msgsize,sender_glob,tag,vd_veloz(sender_loc)%req)
  endif

end subroutine recv_vol_data


subroutine wait_vol_recv()
  use io_server
  use my_mpi
  use constants, only: nproc_io
  implicit none

  integer :: i

  do i = 0, nproc_io-1
    if (vd_pres(i)%req    /= VAL_NOT_ASSIGNED) call wait_req(vd_pres(i)%req)
    if (vd_divglob(i)%req /= VAL_NOT_ASSIGNED) call wait_req(vd_divglob(i)%req)
    if (vd_div(i)%req     /= VAL_NOT_ASSIGNED) call wait_req(vd_div(i)%req)
    if (vd_curlx(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_curlx(i)%req)
    if (vd_curly(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_curly(i)%req)
    if (vd_curlz(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_curlz(i)%req)
    if (vd_velox(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_velox(i)%req)
    if (vd_veloy(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_veloy(i)%req)
    if (vd_veloz(i)%req   /= VAL_NOT_ASSIGNED) call wait_req(vd_veloz(i)%req)
  enddo

end subroutine wait_vol_recv


subroutine write_vol_data(it_io, val_type_mov)
  use io_server
  use specfem_par
  use phdf5_utils
  implicit none

  integer, intent(in) :: it_io
  logical, dimension(5), intent(inout) :: val_type_mov
  integer :: i,j, num_max_type=5, comm, info, id_loc

  ! make output file
  character(len=10) :: tempstr
  character(len=64) :: dset_name
  character(len=64) :: group_name

  type(h5io) :: h5
  h5 = h5io()

  ! get mpi parameters
  call world_get_comm(comm)
  call get_info_null(info)

  ! initialization of h5 file
  call h5_init(h5, fname_h5_data_vol)
  call h5_set_mpi_info(h5, comm, info, myrank, NPROC)

  ! create a hdf5 file
  call h5_open_file(h5)

  ! create time group in h5
  write(tempstr, "(i6.6)") it_io
  group_name = "it_"//tempstr
  call h5_create_group(h5, group_name)
  call h5_open_group(h5, group_name)

  ! loop to write the volume data for each process
  do i = 0, NPROC-1

    id_loc = id_proc_glob2loc(i)

    if (id_loc >= 0) then ! if i (proc) is assigned to this io node

      ! create or open a processor subgroup
      write(tempstr, "(i6.6)") i
      group_name = "proc_"//tempstr
      call h5_create_subgroup(h5, group_name)
      call h5_open_subgroup(h5, group_name)

      ! loop for each value type
      do j=1, num_max_type
        if (val_type_mov(j) .eqv. .true.) then
          if (j==1) then
            dset_name = "pressure"
            call h5_write_dataset_1d_d(h5, dset_name, vd_pres(id_loc)%d1darr)
            call h5_close_dataset(h5)

          elseif (j==2) then
            dset_name = "div_glob"
            call h5_write_dataset_1d_d(h5, dset_name, vd_divglob(id_loc)%d1darr)
            call h5_close_dataset(h5)

          elseif (j==3) then
            dset_name = "div"
            call h5_write_dataset_1d_d(h5, dset_name, vd_div(id_loc)%d1darr)
            call h5_close_dataset(h5)

          elseif (j==4) then
            dset_name = "curl_x"
            call h5_write_dataset_1d_d(h5, dset_name, vd_curlx(id_loc)%d1darr)
            call h5_close_dataset(h5)

            dset_name = "curl_y"
            call h5_write_dataset_1d_d(h5, dset_name, vd_curly(id_loc)%d1darr)
             call h5_close_dataset(h5)

            dset_name = "curl_z"
            call h5_write_dataset_1d_d(h5, dset_name, vd_curlz(id_loc)%d1darr)
            call h5_close_dataset(h5)

          else
            dset_name = "velo_x"
            call h5_write_dataset_1d_d(h5, dset_name, vd_velox(id_loc)%d1darr)
            call h5_close_dataset(h5)

            dset_name = "velo_y"
            call h5_write_dataset_1d_d(h5, dset_name, vd_veloy(id_loc)%d1darr)
            call h5_close_dataset(h5)

            dset_name = "velo_z"
            call h5_write_dataset_1d_d(h5, dset_name, vd_veloz(id_loc)%d1darr)
            call h5_close_dataset(h5)
          endif

        endif
      enddo

    call h5_close_subgroup(h5)

    endif ! end if i (proc) is assigned to this io node
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
  integer                                           :: ier, nfaces_actual, nfaces_aug=16,nnodes_per_face_aug=4
  integer                                           :: len_array_aug
  character(len=64)                                 :: dset_name
  character(len=64)                                 :: group_name

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
  call recv_i_inter(size_surf_array, 1, 0, io_tag_surface_coord_len)
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
  integer                                           :: ier, nfaces_actual, nfaces_aug=16,nnodes_per_face_aug=4
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

  integer                      :: ier, iproc, nrec_local_temp
  integer, dimension(1)        :: nrec_temp
  integer,dimension(0:NPROC-1) :: islice_num_rec_local


  call recv_i_inter(nrec_temp, 1, 0, io_tag_num_recv)
  nrec = nrec_temp(1)

  call recv_dp_inter(t0, 1, 0, io_tag_seismo_tzero)

  do iproc = 0, NPROC-1
    call recv_i_inter(nrec_local_temp, 1, iproc, io_tag_local_rec)
    islice_num_rec_local(iproc) = nrec_local_temp
  enddo

end subroutine get_receiver_info


subroutine allocate_seismo_arrays(islice_num_rec_local)
  use specfem_par
  use io_server

  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: islice_num_rec_local
  integer                                   :: ier, max_num_rec, nstep_temp

  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS >= NSTEP) then
    nstep_temp = NSTEP
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
    allocate(seismo_disp(NDIM,nrec,nstep_temp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seimo_disp')
    if (ier /= 0) stop 'error allocating array seismo_disp'
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    allocate(seismo_velo(NDIM,nrec,nstep_temp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_velo')
    if (ier /= 0) stop 'error allocating array seismo_velo'
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    allocate(seismo_acce(NDIM,nrec,nstep_temp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_acce')
    if (ier /= 0) stop 'error allocating array seismo_acce'
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    allocate(seismo_pres(nrec,nstep_temp),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_pres')
    if (ier /= 0) stop 'error allocating array seismo_pres'
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
      allocate(seismo_temp(NDIM,nrec_passed,time_window),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_temp')
      if (ier /= 0) stop 'error allocating array seismo_temp'

      call recvv_cr_inter(seismo_temp, msg_size, sender, tag)

    ! get scalar value i.e. pres
    else
      time_window = int(msg_size/1/nrec_passed)
      allocate(seismo_temp(1,nrec_passed,time_window),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array seismo_temp')
      if (ier /= 0) stop 'error allocating array seismo_temp'
      call recvv_cr_inter(seismo_temp, msg_size, sender, tag)
    endif

    ! set local array to the global array
    do irec_passed=1,nrec_passed
      id_rec_glob = id_rec_globs(irec_passed,sender)
      ! disp
      if (tag == io_tag_seismo_body_disp) then
        seismo_disp(:,id_rec_glob,:) = seismo_temp(:,irec_passed,:)
      ! velo
      elseif (tag == io_tag_seismo_body_velo) then
        seismo_velo(:,id_rec_glob,:) = seismo_temp(:,irec_passed,:) ! id_rec
      ! acce
      elseif (tag == io_tag_seismo_body_acce) then
        seismo_acce(:,id_rec_glob,:) = seismo_temp(:,irec_passed,:)
      ! pres
      elseif (tag == io_tag_seismo_body_pres) then
        seismo_pres(id_rec_glob,:) = seismo_temp(1,irec_passed,:)
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
  real(kind=CUSTOM_REAL), dimension(NSTEP)                :: time_array
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
  do i = 1, NSTEP
    if (SIMULATION_TYPE == 1) then ! forward simulation ! distinguish between single and double precision for reals
      time_array(i) = real( dble(i-1)*DT - t0 ,kind=CUSTOM_REAL)
    else if (SIMULATION_TYPE == 3) then
      ! adjoint simulation: backward/reconstructed wavefields
      ! distinguish between single and double precision for reals
      ! note: compare time_t with time used for source term
      time_array(i) = real( dble(NSTEP-i)*DT - t0 ,kind=CUSTOM_REAL)
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
    allocate(val_array3d(NDIM,nrec,NSTEP),stat=error)
    call h5_create_dataset_gen(h5, "disp", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    allocate(val_array3d(NDIM,nrec,NSTEP),stat=error)
    call h5_create_dataset_gen(h5, "velo", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    allocate(val_array3d(NDIM,nrec,NSTEP),stat=error)
    call h5_create_dataset_gen(h5, "acce", shape(val_array3d), 3, CUSTOM_REAL)
    deallocate(val_array3d)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    allocate(val_array2d(nrec,NSTEP),stat=error)
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
  if (it_offset+NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) then
    t_upper = NSTEP - it_offset
  else
    t_upper = NTSTEP_BETWEEN_OUTPUT_SEISMOS
  endif

  ! writes out this seismogram
  if (SAVE_SEISMOGRAMS_DISPLACEMENT) then
    component = 'disp'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_disp(:,:,1:t_upper), (/0, 0, it_offset/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_VELOCITY) then
    component = 'velo'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_velo(:,:,1:t_upper), (/0, 0, it_offset/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_ACCELERATION) then
    component = 'acce'
    call h5_write_dataset_3d_r_collect_hyperslab(h5, component, seismo_acce(:,:,1:t_upper), (/0, 0, it_offset/), .false.)
  endif
  if (SAVE_SEISMOGRAMS_PRESSURE) then
    component = 'pres'
    call h5_write_dataset_2d_r_collect_hyperslab(h5, component, seismo_pres(:,1:t_upper), (/0, it_offset/), .false.)
  endif

  call h5_close_file(h5)

end subroutine write_seismograms_io

!
! xdmf output routines
!
subroutine write_xdmf_surface_header()
  use specfem_par
  use io_server
  implicit none
  integer :: num_elm, num_nodes

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    num_nodes = size(surf_x)
  else
    num_nodes = size(surf_x_aug)
  endif

  num_elm = num_nodes/4

  ! writeout xdmf file for surface movie
  fname_xdmf_surf = trim(OUTPUT_FILES)//"/movie_surface.xmf"

  open(unit=xdmf_surf, file=fname_xdmf_surf)
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


subroutine write_xdmf_surface_body(it_io)
  use specfem_par
  use io_server

  implicit none

  integer, intent(in) :: it_io
  integer             :: i

  character(len=20)  :: it_str
  character(len=20)  :: temp_str

  integer :: num_nodes

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    num_nodes = size(surf_x)
  else
    num_nodes = size(surf_x_aug)
  endif

  ! create a group for each io step

  ! open xdmf file
  open(unit=xdmf_surf, file=fname_xdmf_surf)

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

end subroutine write_xdmf_surface_body


subroutine write_xdmf_shakemap()
  use specfem_par
  use io_server
  implicit none
  integer :: num_elm, num_nodes

  if (.not. USE_HIGHRES_FOR_MOVIES) then
    num_nodes = size(surf_x)
  else
    num_nodes = size(surf_x_aug)
  endif

  num_elm = num_nodes/4

  ! writeout xdmf file for surface movie
  fname_xdmf_shake = trim(OUTPUT_FILES)//"/shakemap.xmf"

  open(unit=xdmf_shake, file=fname_xdmf_shake)

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


subroutine write_xdmf_vol_header(nelm_par_proc,nglob_par_proc)
  use specfem_par
  use io_server
  implicit none

  integer, dimension(0:NPROC-1), intent(in) :: nelm_par_proc, nglob_par_proc
  character(len=20)                         :: proc_str, it_str,nelm, nglo
  integer                                   :: iproc, iiout, nout

  ! writeout xdmf file for volume movie
  fname_xdmf_vol = trim(OUTPUT_FILES)//"/movie_volume.xmf"

  open(unit=xdmf_vol, file=fname_xdmf_vol)

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
  do iproc=0,NPROC-1
    nelm=i2c(nelm_par_proc(iproc)*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1))
    nglo=i2c(nglob_par_proc(iproc))
    write(proc_str, "(i6.6)") iproc

    write(xdmf_vol,*) '<Grid Name="mesh_'//trim(proc_str)//'">'
    write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm)//'">'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                           //trim(nelm)//' 9">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/external_mesh.h5:/proc_'&
                           //trim(proc_str)//'/spec_elm_conn_xdmf'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '</Topology>'
    write(xdmf_vol,*) '<Geometry GeometryType="X_Y_Z">'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/external_mesh.h5:/proc_'//trim(proc_str)//'/xstore_dummy'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/external_mesh.h5:/proc_'//trim(proc_str)//'/ystore_dummy'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                        //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
    write(xdmf_vol,*) '       ./DATABASES_MPI/external_mesh.h5:/proc_'//trim(proc_str)//'/zstore_dummy'
    write(xdmf_vol,*) '    </DataItem>'
    write(xdmf_vol,*) '</Geometry>'
    write(xdmf_vol,*) '</Grid>'
  enddo

  ! loop for writing xml includes of timestep information
  nout = int(NSTEP/NTSTEP_BETWEEN_FRAMES) +1
  write(xdmf_vol,*) '</Grid>' ! close mesh info
  write(xdmf_vol,*) '<!-- time series data -->'
  write(xdmf_vol,*) '<Grid Name="results" GridType="Collection" CollectionType="Temporal">'
    do iiout = 1,nout-1
      write(it_str, "(i6.6)") iiout*NTSTEP_BETWEEN_FRAMES
      write(xdmf_vol,*) '    <xi:include href="it_'//trim(it_str)//'.xmf" />'
    enddo
  write(xdmf_vol,*) '</Grid>'

  write(xdmf_vol,*) '</Domain>'
  write(xdmf_vol,*) '</Xdmf>'

  close(xdmf_vol)

end subroutine write_xdmf_vol_header


subroutine write_xdmf_vol_body(it_io,nelm_par_proc, nglob_par_proc, val_type_mov)
  use specfem_par
  use io_server
  implicit none

  integer, intent(in)                       :: it_io
  integer, dimension(0:NPROC-1), intent(in) :: nelm_par_proc, nglob_par_proc
  logical, dimension(5), intent(in)         :: val_type_mov
  character(len=20) :: it_str, proc_str, type_str, type_str1, type_str2, nglo
  integer           :: itype,iproc,idionod
  character(len=5)  :: ioidstr

  ! writeout xdmf file for volume movie
  write(it_str, "(i6.6)") it_io

  open(unit=xdmf_vol_step, file=fname_xdmf_vol_step, position="append", action="write")


  do iproc=0, NPROC-1
    write(proc_str, "(i6.6)") iproc

    idionod = mod(iproc,NIONOD)
    write(ioidstr, "(i5.5)") idionod

    nglo=i2c(nglob_par_proc(iproc))

    write(xdmf_vol_step, *)  '<Grid Name="data_'//trim(proc_str)//'" Type="Uniform">'
    write(xdmf_vol_step, *)  '    <Topology Reference="/Xdmf/Domain/Grid[@Name=''mesh'']/Grid[@Name=''mesh_'&
                                      //trim(proc_str)//''']/Topology" />'
    write(xdmf_vol_step, *)  '    <Geometry Reference="/Xdmf/Domain/Grid[@Name=''mesh'']/Grid[@Name=''mesh_'&
                                      //trim(proc_str)//''']/Geometry" />'

    do itype=1,5
      if (val_type_mov(itype)) then

        if (itype < 4) then

          ! write pressure
          if (itype == 1) then
             type_str = "pressure"
          ! write div_glob
          elseif (itype == 2) then
             type_str = "div_glob"
          ! write div
          elseif (itype == 3) then
            type_str = "div"
          endif

          write(xdmf_vol_step, *)  '    <Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Node">'
          write(xdmf_vol_step, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                              //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
          write(xdmf_vol_step, *)  '            ./DATABASES_MPI/movie_volume_'//ioidstr//'.h5:/it_'&
                                                 //trim(it_str)//'/proc_'//trim(proc_str)//'/'//trim(type_str)
          write(xdmf_vol_step, *)  '        </DataItem>'
          write(xdmf_vol_step, *)  '    </Attribute>'

        else  ! curl or velocity
          ! write curl xyz
          if (val_type_mov(itype) .and. itype == 4) then
            type_str  = "curl_x"
            type_str1 = "curl_y"
            type_str2 = "curl_z"
          else if (val_type_mov(itype) .and. itype == 5) then
            ! write velocity xyz
            type_str  = "velo_x"
            type_str1 = "velo_y"
            type_str2 = "velo_z"
          endif
          ! x
          write(xdmf_vol_step, *)  '    <Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Node">'
          write(xdmf_vol_step, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
          write(xdmf_vol_step, *)  '            ./DATABASES_MPI/movie_volume_'//ioidstr//'.h5:/it_'&
                                                //trim(it_str)//'/proc_'//trim(proc_str)//'/'//trim(type_str)
          write(xdmf_vol_step, *)  '        </DataItem>'
          write(xdmf_vol_step, *)  '    </Attribute>'
          ! y
          write(xdmf_vol_step, *)  '    <Attribute Name="'//trim(type_str1)//'" AttributeType="Scalar" Center="Node">'
          write(xdmf_vol_step, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
          write(xdmf_vol_step, *)  '            ./DATABASES_MPI/movie_volume_'//ioidstr//'.h5:/it_'&
                                                //trim(it_str)//'/proc_'//trim(proc_str)//'/'//trim(type_str1)
          write(xdmf_vol_step, *)  '        </DataItem>'
          write(xdmf_vol_step, *)  '    </Attribute>'
          ! z
          write(xdmf_vol_step, *)  '    <Attribute Name="'//trim(type_str2)//'" AttributeType="Scalar" Center="Node">'
          write(xdmf_vol_step, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                                //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nglo)//'">'
          write(xdmf_vol_step, *)  '            ./DATABASES_MPI/movie_volume_'//ioidstr//'.h5:/it_'&
                                                //trim(it_str)//'/proc_'//trim(proc_str)//'/'//trim(type_str2)
          write(xdmf_vol_step, *)  '        </DataItem>'
          write(xdmf_vol_step, *)  '    </Attribute>'

        endif

      endif ! if vol_type_mov == true
    enddo
    write(xdmf_vol_step, *)  '</Grid>'
  enddo

  close(xdmf_vol_step)

end subroutine write_xdmf_vol_body

subroutine write_xdmf_vol_body_header(it_io)
  use specfem_par
  use io_server
  implicit none
  integer, intent(in) :: it_io
  character(len=20)   :: it_str

  write(it_str, "(i6.6)") it_io
  fname_xdmf_vol_step = trim(OUTPUT_FILES)//"it_"//trim(it_str)//".xmf"

  open(unit=xdmf_vol_step, file=fname_xdmf_vol_step)
  write(xdmf_vol_step,*) '<Grid Name="result"  GridType="Collection"  CollectionType="Spatial">'
  write(xdmf_vol_step,*) '<Time Value="'//trim(r2c(sngl((it_io-1)*DT-t0)))//'" />'

  close(xdmf_vol_step)
end subroutine write_xdmf_vol_body_header


subroutine write_xdmf_vol_body_close()
  use specfem_par
  use io_server
  implicit none

  open(unit=xdmf_vol_step, file=fname_xdmf_vol_step, position="append", action="write")
  write(xdmf_vol_step, *) '</Grid>'
  close(xdmf_vol_step)
end subroutine write_xdmf_vol_body_close


subroutine pass_info_to_io()
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use specfem_par_movie
  use constants, only: dest_ionod

  implicit none

  integer ::  n_msg_vol_each_proc = 0,irec,irec_local,i_ionod
  integer,dimension(nrec_local) :: tmp_irec

  ! initialization of io node from compute node side


  ! send the receiver information to io node for the outputs of seismo signals
  ! seismo data is out only from rank=0 of io nodes
  if (myrank == 0) then
    ! send nrec and nrec_local
    call send_i_inter((/nrec/), 1, 0, io_tag_num_recv)
    ! send t0
    call send_dp_inter((/t0/), 1, 0, io_tag_seismo_tzero)
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
      do i_ionod=0,NIONOD-1
        call send_i_inter((/n_msg_vol_each_proc/),1,i_ionod,io_tag_vol_nmsg)
      enddo
    endif
  endif ! end if myrank == 0

  if (MOVIE_VOLUME) then
    ! send nspec and nglob in each process
    call send_i_inter((/NSPEC_AB/),1,0,io_tag_vol_nspec)
    call send_i_inter((/NGLOB_AB/),1,0,io_tag_vol_nglob)
  endif

end subroutine pass_info_to_io
