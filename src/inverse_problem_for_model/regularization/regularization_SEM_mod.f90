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

module regularization

  !! IMPORT SPECFEM VARIABLES ------------------------------------------------------------------------------------------------------
  use specfem_par, only: prname, LOCAL_PATH, myrank, IIN, NDIM, IMAIN, NGNOD, NSPEC_AB, &
       NGLOB_AB, NGNOD2D, NGNOD, NGLLX, NGLLY, NGLLZ, GAUSSALPHA, GAUSSBETA, xigll, &
       NPROC, SAVE_MOHO_MESH, NUMBER_OF_SIMULTANEOUS_RUNS, ibool, xstore, ystore, zstore, &
       MAX_STRING_LEN, num_interfaces_ext_mesh_sp => num_interfaces_ext_mesh, &
       max_nibool_interfaces_ext_mesh_sp => max_nibool_interfaces_ext_mesh, &
       nibool_interfaces_ext_mesh_sp => nibool_interfaces_ext_mesh, &
       ibool_interfaces_ext_mesh_sp => ibool_interfaces_ext_mesh, &
       my_neighbors_ext_mesh_sp => my_neighbors_ext_mesh

  !! IMPORT INVERSE_PROBLEM VARIABLES ----------------------------------------------------------------------------------------------
  use inverse_problem_par

  implicit none
  integer, public                                                            :: mem_used_for_reg
  integer, public                                                            :: Nb_iglob_on_faces

  type, private :: list_iglob
     integer                                                                 :: Niglob
     integer,                               dimension(:),       allocatable  :: list_of_iglob
  end type list_iglob

  type, private :: comm_fd
     integer                                                                 :: ns, nr, ibegin
     integer,                              dimension(:),        allocatable  :: iglob_to_send
     real(kind=CUSTOM_REAL),               dimension(:),        allocatable  :: array_to_send
     real(kind=CUSTOM_REAL),               dimension(:),        allocatable  :: array_to_recv
  end type comm_fd

  !! LOCAL VARIABLES ---------------------------------------------------------------------------------------------------------------
  type(comm_fd),                 private, dimension(:),         allocatable  :: struct_comm
  integer,                       private, dimension(:),         allocatable  :: request_send_scalar,request_recv_scalar

  !! TEMPORARY LOCAL VARIABLES -----------------------------------------------------------------------------------------------------
  !! for reading mesh partition (same as used in generate_databases)
  integer,                       private                                     :: boundary_number !!!!!! ,MAX_NEIG
!!!!!!  integer,                       private                                     :: MAX_POINTS_ALLOWED=1000
  integer,                       private                                     :: nspec2D_xmin, nspec2D_xmax
  integer,                       private                                     :: nspec2D_ymin, nspec2D_ymax
  integer,                       private                                     :: nspec2D_bottom, nspec2D_top
  integer,                       private                                     :: nspec2D_bottom_ext, nspec2D_top_ext
  integer,                       private                                     :: max_interface_size_ext_mesh
  integer,                       private                                     :: num_interfaces_ext_mesh
  integer,                       private                                     :: nmat_ext_mesh, nundefMat_ext_mesh
  integer,                       private, dimension(:),         allocatable  :: my_nelmnts_neighbors_ext_mesh
  integer,                       private, dimension(:),         allocatable  :: my_neighbors_ext_mesh
  integer,                       private, dimension(:),         allocatable  :: nibool_interfaces_ext_mesh
  integer,                       private, dimension(:),         allocatable  :: ibelm_xmin, ibelm_xmax
  integer,                       private, dimension(:),         allocatable  :: ibelm_ymin, ibelm_ymax
  integer,                       private, dimension(:),         allocatable  :: ibelm_bottom, ibelm_top
  integer,                       private, dimension(:,:),       allocatable  :: elmnts_ext_mesh
  integer,                       private, dimension(:,:),       allocatable  :: mat_ext_mesh
  integer,                       private, dimension(:,:),       allocatable  :: ibool_interfaces_ext_mesh
  integer,                       private, dimension(:,:),       allocatable  :: nodes_ibelm_xmin, nodes_ibelm_xmax
  integer,                       private, dimension(:,:),       allocatable  :: nodes_ibelm_ymin, nodes_ibelm_ymax
  integer,                       private, dimension(:,:),       allocatable  :: nodes_ibelm_bottom, nodes_ibelm_top
  integer,                       private, dimension(:,:,:),     allocatable  :: my_interfaces_ext_mesh
  double precision,              private, dimension(:,:),       allocatable  :: nodes_coords_ext_mesh
  double precision,              private, dimension(:,:),       allocatable  :: mat_prop
  character(len=MAX_STRING_LEN), private, dimension(:,:),       allocatable  :: undef_mat_prop

  !! MPI communication
!!!!!!  integer,                       private, dimension(:),         allocatable  :: indx_send
  integer,                       private, dimension(:),         allocatable  :: indx_recv

  !! auxiliary arrays to build MPI communication
!!!!!!  integer,                       private, dimension(:),         allocatable  :: ispec_to_send
!!!!!!  integer,                       private, dimension(:,:,:,:),   allocatable  :: iglob_to_send, iglob_to_recv
!!!!!!  integer,                       private, dimension(:,:,:,:),   allocatable  :: iglob_store
!!!!!!  real(kind=CUSTOM_REAL),        private, dimension(:,:,:,:),   allocatable  :: xgrid, ygrid, zgrid
!!!!!!  real(kind=CUSTOM_REAL),        private, dimension(:,:,:,:),   allocatable  :: xgrid_to_send, ygrid_to_send, zgrid_to_send
!!!!!!  real(kind=CUSTOM_REAL),        private, dimension(:,:,:,:),   allocatable  :: xgrid_to_recv, ygrid_to_recv, zgrid_to_recv
!!!!!!  real(kind=CUSTOM_REAL),        private                                     :: radius_disk
!!!!!!  type(list_iglob),              private, dimension(:),         allocatable  :: list_iglob_neig
!!!!!!  type(list_iglob),              private, dimension(:),         allocatable  :: list_iglob_to_send

  !! for derivatives in gll
  integer,                       private,              parameter             :: NGLLd = 10
  double precision,              private,              parameter             :: ONE_EIGHTH = 0.125d0, ONE = 1.d0, ZERO = 1.d0
  double precision,              private, dimension(NGLLd)                   :: gll_points, wgll_points
  double precision,              private, dimension(NGLLd,NGLLd)             :: hlagrange_prime
  double precision,              private, dimension(NGLLX,NGLLd)             :: hlagrange
  double precision,              private, dimension(NGLLd,NGLLX)             :: hlagrange_old
  double precision,              private, dimension(8,NGLLd,NGLLd,NGLLd)     :: shape_function
  double precision,              private, dimension(:,:,:,:,:),allocatable   :: dershape_function
  double precision,              private, dimension(8)                       :: xnodelm, ynodelm, znodelm
  double precision,              private, dimension(NGLLX,NGLLX,NGLLX)       :: field_initial
  double precision,              private, dimension(NGLLd,NGLLd,NGLLd)       :: dfdx, dfdy, dfdz
  double precision,              private, dimension(NGLLd,NGLLd,NGLLd)       :: xstore_interp, ystore_interp, zstore_interp
  double precision,              private, dimension(NGLLd,NGLLd,NGLLd)       :: field_interpolated, Laplacian_of_field_interpolated
  double precision,              private, dimension(:,:,:,:,:),allocatable   :: Jacobian_shape_function

  !! for checking
  real(kind=CUSTOM_REAL),        private                                      :: lambda

  !! memory evaluation
  integer,                       private                                     :: mem_tmp, max_memory_used

  !---------------------------------------------------------------------------------------------------------------------------------

  public  :: initialize_regularization, compute_laplacian, compute_laplacian2_vp_vs_rho
!!                                       DK DK unused for now
  private :: read_partition_files, & !!  find_all_glls_in_disk, prepare_MPI_regularization, find_list_of_element_adj, &
!! DK DK unused for now   add_or_not_one_more, count_gll_in_face, store_list_gll_in_face, iglob_not_used, testing_regularization, &
!! DK DK unused for now   send_recv_blocking, check_regularization_on_mesh, check_field_accuracy, compute_cosine_field, &
             send_recv_blocking, check_regularization_on_mesh, compute_cosine_field, &
!! DK DK unused for now   compute_first_derivatives_lagrange, compute_laplac_lagrange, compute_and_store_FD_derivatives_matrix, &
             compute_first_derivatives_lagrange, compute_laplac_lagrange, &
!! DK DK unused for now   compute_laplacian_FD, line_vandermonde, compute_M_transpose_times_M, compute_inv_MtM_times_Mt, &
             compute_laplacian_FD, &
!! DK DK unused for now   matrix_times_vector, inverse
             matrix_times_vector

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  set up MPI commnication for regularization and define derivatives matrices
!---------------------------------------------------------------------------------

!!!!!!!!!!!!  subroutine initialize_regularization(regularization_fd)
  subroutine initialize_regularization()

    implicit none

!!!!!!!!!!!!    type(regul), dimension(:), allocatable, intent(inout) :: regularization_fd

    integer :: num_interface
!!!!!!!!!!!!    integer :: ispec, inode , ie, imat, iface, icorner

    if (DEBUG_MODE) then
       write(IIDD,*) ' reading partition files ...'
    endif
    mem_used_for_reg = 0
    max_memory_used = 0
    ! read the origial mesh partition
    call read_partition_files()

    ! prepare arrays for MPI comm
    !!  NOT USE IT FOR NOW
    !! call prepare_MPI_regularization(regularization_fd)

    ! compute and store Derivative Matrices for each gll
    !! NOT USE IT FOR NOW
    !! call compute_and_store_FD_derivatives_matrix(regularization_fd)

    !! TO DO we may be can avoid to store all matrices by computing them each time we need
    !! in case of memory issues (because it can be require lot of space)

    call setup_interpolation_for_higher_degree()

    !!! check the accuray of FD derivatives on cosine function
!!!!!!!!!!!!!!!    call check_regularization_on_mesh(regularization_fd)
    call check_regularization_on_mesh()

    if (DEBUG_MODE) then
       !call testing_regularization(regularization_fd)
       write(IIDD,*) ' reading partition files : passed  '
       write(IIDD,*)
       write(IIDD,'(a26, i8, a17, i8)') &
       'Number of MPI interfaces  ', num_interfaces_ext_mesh, '  Maximum element ', max_interface_size_ext_mesh
       do  num_interface = 1, num_interfaces_ext_mesh
          write(IIDD,'(5x, a14, i6, a20, i8)') &
          'id neighbors ', my_neighbors_ext_mesh(num_interface), ' number of elements ', &
          my_nelmnts_neighbors_ext_mesh(num_interface)
!!$          do ie = 1, my_nelmnts_neighbors_ext_mesh(num_interface)
!!$             write(IIDD,'(a7, 2i8)') 'TOSEND ', my_interfaces_ext_mesh(1,ie,num_interface), &
!!$          my_interfaces_ext_mesh(2,ie,num_interface)
!!$          enddo
       enddo

       write(IIDD,*)
       write(IIDD,*)
       write(IIDD,*) '               max memory  : ' ,max_memory_used, max_memory_used / 1000000, ' Mo'
       write(IIDD,*) '           current memory  : ' ,mem_used_for_reg , mem_used_for_reg / 1000000, ' Mo'
       write(IIDD,*) '         memory for grid  : ' ,CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*NSPEC_AB, &
            CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*NSPEC_AB / 1000000, ' Mo'
       write(IIDD,*) '         memory for field : ' ,CUSTOM_REAL*NGLOB_AB, CUSTOM_REAL*NGLOB_AB / 1000000 , ' Mo '
       write(IIDD,*) ' memory used for structure : ', mem_tmp, &
            mem_tmp/1000000, ' Mo'

    endif

    max_memory_used =  mem_used_for_reg

  end subroutine initialize_regularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH MANIPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  set up temporary MPI communication for regularization
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

!  subroutine prepare_MPI_regularization(regularization_fd)

!    use my_mpi

!    integer :: i,j,k, ki, iglob
!    integer :: iglob_to_find, irank_to_look
!    integer :: ispec, num_interface, ie
!    integer :: nspec_to_send, nspec_sended, nspec_recevd
!    integer :: nspec_sended_dummy(1), nspec_recevd_dummy(1)
!    integer :: sendtag0=101, sendtag1=100
!    integer :: recvtag0=100, recvtag1=101
!    integer :: ier, igll, iglob_index, nused_iglob
!    integer :: ifac_loc, nb_gll_in_face
!    integer :: irank, irank0, irank1, Npts
!    integer :: in_index, out_index, igll_local, out_index_to_allocate
!    integer, dimension(:,:),   allocatable :: xcomm, xcomm_dummy, EtoE, EtoMPI, list_of_gll_in_face
!    integer, dimension(:),     allocatable :: nspec_send, nspec_recv, Nv, list_of_gll, rank_of_gll
!    integer, dimension(:),     allocatable :: tmp_l1, tmp_i1, used_iglob, number_points_to_receive_from_neig

!    type(regul), dimension(:), allocatable, intent(inout)  :: regularization_fd
!    real(kind=CUSTOM_REAL)                                 :: distance_min_glob,distance_max_glob
!    real(kind=CUSTOM_REAL)                                 :: elemsize_min_glob,elemsize_max_glob
!    real(kind=CUSTOM_REAL)                                 :: x_min_glob,x_max_glob
!    real(kind=CUSTOM_REAL)                                 :: y_min_glob,y_max_glob
!    real(kind=CUSTOM_REAL)                                 :: z_min_glob,z_max_glob
!    integer,                 dimension(NGNOD)              :: iaddx,iaddy,iaddz

!    real(kind=CUSTOM_REAL)                                 :: dist_sq_min, dist_sq
!    real(kind=CUSTOM_REAL)                                 :: xp,yp,zp, xs,ys,zs
!    logical                                                :: not_used_yet, find_iglob


!    integer, dimension(:), allocatable :: list_of_gll_to_send
!    integer, dimension(:), allocatable :: list_of_gll_to_recv
!    integer, dimension(:), allocatable :: Npts_recv_dummy, sendtag, recvtag
!    integer, dimension(1)              :: Npts_send_dummy


!!!! SEND RECV OVERLAP ---------------------------------------------------------------------------------------------------------
!!! send the neighbor elements grid to neighbor MPI slice
!!!
!    allocate(xcomm(0:NPROC-1,0:NPROC-1), xcomm_dummy(0:NPROC-1,0:NPROC-1), nspec_send(0:NPROC-1), nspec_recv(0:NPROC-1))
!    mem_used_for_reg =  mem_used_for_reg + 2*4*NPROC*NPROC + 2* NPROC

!    xcomm(:,:)=0
!    xcomm_dummy(:,:)=0
!    nspec_send(:)=0
!    nspec_recv(:)=0
!    do num_interface = 1, num_interfaces_ext_mesh
!       xcomm_dummy(myrank,my_neighbors_ext_mesh(num_interface))=1
!    enddo
!    call MPI_ALLREDUCE(xcomm_dummy,xcomm,NPROC*NPROC,MPI_INTEGER,MPI_SUM,my_local_mpi_comm_world,ier)

!    xcomm_dummy(:,:)=0
!    do irank0=0,NPROC-1
!       do irank1=irank0+1,NPROC-1
!          if (xcomm(irank0,irank1) > 0) then

!             if (myrank == irank0) then
!                do num_interface = 1, num_interfaces_ext_mesh
!                   if (irank1 == my_neighbors_ext_mesh(num_interface)) then

!                      nspec_to_send = my_nelmnts_neighbors_ext_mesh(num_interface)
!                      nspec_sended  = 0
!                      allocate(ispec_to_send(nspec_to_send))
!                      mem_used_for_reg = mem_used_for_reg + 4*nspec_to_send
!                      max_memory_used = max(mem_used_for_reg, max_memory_used)

!                      ispec_to_send(:)=0
!                      !! retrieve the elements to send
!                      do ie = 1, nspec_to_send
!                         ispec =  my_interfaces_ext_mesh(1,ie,num_interface)
!                         call add_or_not_one_more(ispec_to_send, ispec, nspec_to_send, nspec_sended)
!                      enddo

!                      nspec_sended_dummy(1)=nspec_sended
!                      nspec_send(irank1)=nspec_sended

!                      call send_i_t(nspec_sended_dummy, 1, irank1)
!                      call recv_i_t(nspec_recevd_dummy, 1, irank1)
!                      nspec_recevd=nspec_recevd_dummy(1)
!                      nspec_recv(irank1)=nspec_recevd

!                      deallocate(ispec_to_send)
!                      mem_used_for_reg = mem_used_for_reg - 4*nspec_to_send


!                   endif
!                enddo
!             endif

!             if (myrank == irank1) then
!                do num_interface = 1, num_interfaces_ext_mesh
!                   if (irank0 == my_neighbors_ext_mesh(num_interface)) then

!                      nspec_to_send = my_nelmnts_neighbors_ext_mesh(num_interface)
!                      nspec_sended  = 0
!                      allocate(ispec_to_send(nspec_to_send))
!                      mem_used_for_reg =  mem_used_for_reg + 4*nspec_to_send
!                      max_memory_used = max( mem_used_for_reg, max_memory_used)

!                      ispec_to_send(:)=0
!                      !! retrieve the elements to send
!                      do ie = 1, nspec_to_send
!                         ispec =  my_interfaces_ext_mesh(1,ie,num_interface)
!                         call add_or_not_one_more(ispec_to_send, ispec, nspec_to_send, nspec_sended)
!                      enddo

!                      nspec_sended_dummy(1)=nspec_sended
!                      nspec_send(irank0)=nspec_sended

!                      call recv_i_t(nspec_recevd_dummy, 1, irank0)
!                      call send_i_t(nspec_sended_dummy, 1, irank0)
!                      nspec_recevd=nspec_recevd_dummy(1)
!                      nspec_recv(irank0)=nspec_recevd

!                      deallocate(ispec_to_send)
!                      mem_used_for_reg = mem_used_for_reg - 4*nspec_to_send

!                   endif
!                enddo
!             endif

!          endif

!          call synchronize_all()

!       enddo
!    enddo


!    !! MPI comm
!    nspec_sended=sum(nspec_send)
!    nspec_recevd=sum(nspec_recv)

!    allocate(indx_recv(0:NPROC))
!    allocate(indx_send(0:NPROC))

!    mem_used_for_reg = mem_used_for_reg + 2* 4 * (NPROC+1)
!    max_memory_used = max(mem_used_for_reg, max_memory_used)

!    indx_recv(:)=0;indx_send(:)=0
!    indx_recv(0)=1;indx_send(0)=1
!    do irank=1,NPROC
!       indx_recv(irank)=indx_recv(irank-1) +  nspec_recv(irank-1)
!       indx_send(irank)=indx_send(irank-1) +  nspec_send(irank-1)
!    enddo

!    if (DEBUG_MODE) then
!       write(IIDD,*)
!       write(IIDD,*) 'MYRANK ', myrank
!       write(IIDD,*)
!       write(IIDD,*) ' communication : MPI rank, nspec to recv, nspec to send'
!       write(IIDD,*)
!       do irank=0, NPROC-1
!          write(IIDD,'(5i10)')  irank, nspec_recv(irank), nspec_send(irank), indx_recv(irank), indx_send(irank)
!       enddo
!       write(IIDD,*)
!       write(IIDD,*)
!       write(IIDD,*)
!    endif


!    allocate(xgrid(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!    allocate(ygrid(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!    allocate(zgrid(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!    allocate(iglob_store(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!    mem_used_for_reg = mem_used_for_reg + CUSTOM_REAL * NGLLX*NGLLY*NGLLZ*nspec_recevd*4 +&
!                                                    4 * NGLLX*NGLLY*NGLLZ*nspec_recevd

!    do irank0=0,NPROC-1
!       do irank1=irank0+1,NPROC-1

!          if (xcomm(irank0,irank1) > 0) then

!             if (myrank == irank0) then

!                do num_interface = 1, num_interfaces_ext_mesh
!                   if (irank1 == my_neighbors_ext_mesh(num_interface)) then

!                      !! retrieve the elements to send
!                      nspec_sended  = nspec_send(irank1)
!                      nspec_recevd  = nspec_recv(irank1)
!                      nspec_to_send = my_nelmnts_neighbors_ext_mesh(num_interface)

!                      allocate(ispec_to_send(nspec_to_send))

!                      allocate(iglob_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                      allocate(iglob_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!                      allocate(xgrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                      allocate(ygrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                      allocate(zgrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))

!                      allocate(xgrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!                      allocate(ygrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!                      allocate(zgrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!                      mem_tmp = 4*nspec_to_send + 4*NGLLX*NGLLY*NGLLZ*nspec_sended + &
!                                4*NGLLX*NGLLY*NGLLZ*nspec_recevd + &
!                                3*CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*nspec_sended + &
!                                3*CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*nspec_recevd

!                      mem_used_for_reg = mem_used_for_reg + mem_tmp
!                      max_memory_used = max( mem_used_for_reg, max_memory_used)

!                      nspec_sended=0
!                      ispec_to_send(:)=0
!                      do ie = 1, nspec_to_send
!                         ispec =  my_interfaces_ext_mesh(1,ie,num_interface)

!                         call add_or_not_one_more(ispec_to_send, ispec, nspec_to_send, nspec_sended)
!                      enddo
!                      !------------------------------------------------------

!                      !! store positions to send
!                      do ie = 1, nspec_sended
!                         ispec=ispec_to_send(ie)
!                         do k=1,NGLLZ
!                            do j=1,NGLLY
!                               do i=1,NGLLX
!                                  iglob=ibool(i,j,k,ispec)
!                                  iglob_to_send(i,j,k,ie) =  iglob
!                                  xgrid_to_send(i,j,k,ie) =  xstore(iglob)
!                                  ygrid_to_send(i,j,k,ie) =  ystore(iglob)
!                                  zgrid_to_send(i,j,k,ie) =  zstore(iglob)
!                               enddo
!                            enddo
!                         enddo
!                      enddo
!                      !------------------------------------------------------
!                      nspec_sended_dummy(1)=nspec_sended

!                      call send_i(iglob_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank1, sendtag1)
!                      call recv_i(iglob_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank1, recvtag1)

!                      call sendv_cr(xgrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank1, sendtag1)
!                      call recvv_cr(xgrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank1, recvtag1)

!                      call sendv_cr(ygrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank1, sendtag1)
!                      call recvv_cr(ygrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank1, recvtag1)

!                      call sendv_cr(zgrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank1, sendtag1)
!                      call recvv_cr(zgrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank1, recvtag1)

!                      !! store
!                      iglob_store(:,:,:,indx_recv(irank1):indx_recv(irank1+1)-1)=iglob_to_recv(:,:,:,:)
!                      xgrid(:,:,:,indx_recv(irank1):indx_recv(irank1+1)-1)=xgrid_to_recv(:,:,:,:)
!                      ygrid(:,:,:,indx_recv(irank1):indx_recv(irank1+1)-1)=ygrid_to_recv(:,:,:,:)
!                      zgrid(:,:,:,indx_recv(irank1):indx_recv(irank1+1)-1)=zgrid_to_recv(:,:,:,:)

!                      !----------------- deallocate -------------------------
!                      deallocate(ispec_to_send, iglob_to_send, iglob_to_recv)
!                      deallocate(xgrid_to_send, ygrid_to_send, zgrid_to_send)
!                      deallocate(xgrid_to_recv, ygrid_to_recv, zgrid_to_recv)
!                      mem_used_for_reg = mem_used_for_reg - mem_tmp

!                  endif
!               enddo

!            endif


!            if (myrank == irank1) then

!               do num_interface = 1, num_interfaces_ext_mesh
!                  if (irank0 == my_neighbors_ext_mesh(num_interface)) then

!                     !! retrieve the elements to send
!                     nspec_sended  = nspec_send(irank0)
!                     nspec_recevd  = nspec_recv(irank0)
!                     nspec_to_send =  my_nelmnts_neighbors_ext_mesh(num_interface)

!                     allocate(ispec_to_send(nspec_to_send))

!                     allocate(iglob_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                     allocate(iglob_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!                     allocate(xgrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                     allocate(ygrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))
!                     allocate(zgrid_to_send(NGLLX,NGLLY,NGLLZ,nspec_sended))

!                     allocate(xgrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!                     allocate(ygrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))
!                     allocate(zgrid_to_recv(NGLLX,NGLLY,NGLLZ,nspec_recevd))

!                     mem_tmp = 4*nspec_to_send + 4*NGLLX*NGLLY*NGLLZ*nspec_sended + &
!                               4*NGLLX*NGLLY*NGLLZ*nspec_recevd+&
!                               3*CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*nspec_sended+&
!                               3*CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*nspec_recevd
!                     mem_used_for_reg = mem_used_for_reg + mem_tmp
!                     max_memory_used = max( mem_used_for_reg, max_memory_used)

!                     nspec_sended=0
!                     ispec_to_send(:)=0
!                     do ie = 1, nspec_to_send
!                        ispec =  my_interfaces_ext_mesh(1,ie,num_interface)
!                        call add_or_not_one_more(ispec_to_send, ispec, nspec_to_send, nspec_sended)
!                     enddo
!                     !----------------------------------------------------------------

!                     !! store position to send
!                     do ie = 1, nspec_sended
!                        ispec=ispec_to_send(ie)
!                        do k=1,NGLLZ
!                           do j=1,NGLLY
!                              do i=1,NGLLX
!                                 iglob=ibool(i,j,k,ispec)
!                                 iglob_to_send(i,j,k,ie) =  iglob
!                                 xgrid_to_send(i,j,k,ie) =  xstore(iglob)
!                                 ygrid_to_send(i,j,k,ie) =  ystore(iglob)
!                                 zgrid_to_send(i,j,k,ie) =  zstore(iglob)
!                              enddo
!                           enddo
!                        enddo
!                     enddo
!                     !-----------------------------------------------------------------
!                     nspec_sended_dummy(1)=nspec_sended

!                     call recv_i(iglob_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank0 ,recvtag0)
!                     call send_i(iglob_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank0, sendtag0)

!                     call recvv_cr(xgrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank0, recvtag0)
!                     call sendv_cr(xgrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank0, sendtag0)

!                     call recvv_cr(ygrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank0, recvtag0)
!                     call sendv_cr(ygrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank0, sendtag0)

!                     call recvv_cr(zgrid_to_recv, NGLLX*NGLLY*NGLLZ*nspec_recevd, irank0, recvtag0)
!                     call sendv_cr(zgrid_to_send, NGLLX*NGLLY*NGLLZ*nspec_sended, irank0, sendtag0)

!                     !! store
!                     iglob_store(:,:,:,indx_recv(irank0):indx_recv(irank0+1)-1)=iglob_to_recv(:,:,:,:)

!                     xgrid(:,:,:,indx_recv(irank0):indx_recv(irank0+1)-1)=xgrid_to_recv(:,:,:,:)
!                     ygrid(:,:,:,indx_recv(irank0):indx_recv(irank0+1)-1)=ygrid_to_recv(:,:,:,:)
!                     zgrid(:,:,:,indx_recv(irank0):indx_recv(irank0+1)-1)=zgrid_to_recv(:,:,:,:)

!                     !----------------- deallocate -------------------------
!                     deallocate(ispec_to_send, iglob_to_send, iglob_to_recv)
!                     deallocate(xgrid_to_send, ygrid_to_send, zgrid_to_send)
!                     deallocate(xgrid_to_recv, ygrid_to_recv, zgrid_to_recv)
!                     mem_used_for_reg = mem_used_for_reg - mem_tmp

!                  endif
!               enddo

!            endif

!         endif

!         call synchronize_all()

!      enddo
!   enddo

!   deallocate(xcomm_dummy)
!   mem_used_for_reg = mem_used_for_reg - 4*NPROC*NPROC
!!!! END SEND RECV OVERLAP ---------------------------------------------------------------------------------------------------------


!   ! get mesh properties
!   call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
!   call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
!        x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
!        elemsize_min_glob,elemsize_max_glob, &
!        distance_min_glob,distance_max_glob)

!   !! for neighour search
!   dist_sq_min =  5* elemsize_max_glob**2
!   !! radius for disk
!   radius_disk= 3. * distance_min_glob
!   if (DEBUG_MODE) write(IIDD,*) 'RADIUS DISK :',radius_disk


!   !!  count max neighbor for each elements -----------------------
!   allocate(Nv(NSPEC_AB))

!   mem_used_for_reg = mem_used_for_reg + 4*NSPEC_AB
!   max_memory_used = max(mem_used_for_reg, max_memory_used)

!   Nv(:)=0
!   do ie = 1, NSPEC_AB

!      iglob = ibool(3,3,3,ie)
!      xp = xstore(iglob)
!      yp = ystore(iglob)
!      zp = zstore(iglob)

!      ! regular elements
!      do ispec = 1, NSPEC_AB

!         iglob = ibool(3,3,3,ispec)
!         xs = xstore(iglob)
!         ys = ystore(iglob)
!         zs = zstore(iglob)

!         dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2
!         if (dist_sq < dist_sq_min) Nv(ie)=Nv(ie)+1
!      enddo

!      ! additional element from neighbor MPI slice
!      do ispec = 1, nspec_recevd

!         xs=xgrid(3,3,3,ispec)
!         ys=ygrid(3,3,3,ispec)
!         zs=zgrid(3,3,3,ispec)

!         dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2
!         if (dist_sq < dist_sq_min) Nv(ie)=Nv(ie)+1

!      enddo

!   enddo

!   MAX_NEIG=maxval(Nv)
!   deallocate(Nv)
!   mem_used_for_reg = mem_used_for_reg - 4*NSPEC_AB
!   !! end of count max neighbor for each elements -----------------------

!   if (DEBUG_MODE) write(IIDD,*) 'MAX NEIG ', MAX_NEIG

!   !! search for neighbor elements -------------------------------------
!   allocate(EtoE(NSPEC_AB, MAX_NEIG),EtoMPI(NSPEC_AB,MAX_NEIG))
!   mem_used_for_reg = mem_used_for_reg + 2 *4 * NSPEC_AB * MAX_NEIG
!   max_memory_used = max( mem_used_for_reg, max_memory_used)

!   EtoE(:,:)=-1
!   EtoMPI(:,:)=-1
!   do ie = 1, NSPEC_AB

!      iglob = ibool(3,3,3,ie)
!      xp = xstore(iglob)
!      yp = ystore(iglob)
!      zp = zstore(iglob)
!      ki=0

!      !! look element in my MPI slice
!      do ispec = 1, NSPEC_AB

!         iglob = ibool(3,3,3,ispec)
!         xs = xstore(iglob)
!         ys = ystore(iglob)
!         zs = zstore(iglob)

!         dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2

!         if (dist_sq < dist_sq_min) then
!            ki = ki + 1
!            EtoE(ie, ki)   = ispec
!            EtoMPI(ie,ki) = myrank
!         endif
!      enddo

!      !! look additional element from neighbor MPI slice
!      do irank = 0, NPROC-1
!         do ispec = indx_recv(irank), indx_recv(irank+1)-1

!            xs=xgrid(3,3,3,ispec)
!            ys=ygrid(3,3,3,ispec)
!            zs=zgrid(3,3,3,ispec)

!            dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2
!            if (dist_sq < dist_sq_min) then
!               ki = ki + 1
!               EtoE(ie, ki)  = ispec
!               EtoMPI(ie,ki) = irank
!               !if (DEBUG_MODE) write(IIDD,*) 'EtoE on MPI slice boundary', ie
!            endif
!         enddo
!      enddo

!   enddo
!   if (DEBUG_MODE) write(IIDD,*) 'STORED EtoE '
!   !! end search for neighborg elements -------------------------------------

!   !! a way to run over face gll
!   call count_gll_in_face(nb_gll_in_face)
!   allocate(list_of_gll_in_face(4,nb_gll_in_face))

!   mem_used_for_reg = mem_used_for_reg + 4 * 4 * nb_gll_in_face
!   max_memory_used = max( mem_used_for_reg, max_memory_used)

!   call store_list_gll_in_face(list_of_gll_in_face)

!   !! count number of GLL iglob in all faces----------------
!   allocate(used_iglob(NGLOB_AB))

!   mem_used_for_reg = mem_used_for_reg + 4 * NGLOB_AB
!   max_memory_used = max( mem_used_for_reg, max_memory_used)

!   nused_iglob=0
!   used_iglob(:)=-1
!   do ie = 1, NSPEC_AB
!      !! loop over GLL in faces
!      do igll = 1, nb_gll_in_face

!         i = list_of_gll_in_face(1,igll)
!         j = list_of_gll_in_face(2,igll)
!         k = list_of_gll_in_face(3,igll)
!         ifac_loc = list_of_gll_in_face(4,igll)  !! not needed anylore

!         iglob = ibool(i,j,k,ie)

!         call  add_or_not_one_more(used_iglob, iglob, NGLOB_AB, nused_iglob)
!      enddo
!   enddo
!   !!--------------------------------------

!   used_iglob(:)=-1
!   iglob_index=0

!   allocate(list_of_gll(MAX_POINTS_ALLOWED), rank_of_gll(MAX_POINTS_ALLOWED))
!   allocate(tmp_l1(MAX_NEIG),tmp_i1(MAX_NEIG))

!   mem_used_for_reg = mem_used_for_reg + &
!        2 * 4 * MAX_POINTS_ALLOWED + 2 * 4 * MAX_NEIG


!   !!-- structure to save
!   Nb_iglob_on_faces=nused_iglob
!   allocate(regularization_fd(nused_iglob))
!   mem_used_for_reg = mem_used_for_reg + 4 * 4 * nused_iglob
!   mem_tmp=0
!   !!
!   allocate(number_points_to_receive_from_neig(0:NPROC-1))
!   number_points_to_receive_from_neig(:)=0

!   mem_used_for_reg = mem_used_for_reg + 4 * NPROC
!   max_memory_used = max(mem_used_for_reg, max_memory_used)

!   !! deternine point in disk for each iglob------------------------------
!   do ie =1, NSPEC_AB

!      !! loop over GLL in faces
!      do igll = 1, nb_gll_in_face

!         i = list_of_gll_in_face(1,igll)
!         j = list_of_gll_in_face(2,igll)
!         k = list_of_gll_in_face(3,igll)
!         ifac_loc = list_of_gll_in_face(4,igll) !! not needed anylore

!         iglob = ibool(i,j,k,ie)

!         !! look if we already used this GLL : iglob
!         call iglob_not_used(used_iglob, iglob, nused_iglob, iglob_index, not_used_yet)

!         if (not_used_yet) then

!            xs = xstore(iglob)
!            ys = ystore(iglob)
!            zs = zstore(iglob)

!            !! find all GLL in disk
!            tmp_l1(:)=EtoE(ie,:)    !! list of elements to look
!            tmp_i1(:)=EtoMPI(ie,:)  !! MPI slices where are those elements
!            call find_all_glls_in_disk(xs, ys, zs, tmp_l1, tmp_i1, list_of_gll, rank_of_gll, in_index, out_index)

!            ! store it in structure
!            regularization_fd(iglob_index)%iglob=iglob      !! store the main GLL point where we want to compute FD derivatives
!            regularization_fd(iglob_index)%nReg=in_index    !! numbers of GLL points in my MPI slice used for FD
!            regularization_fd(iglob_index)%nNei=out_index   !! numbres of GLL points out of my MPI slice needed for FD
!            regularization_fd(iglob_index)%MaxOrder=2     !! max order FD derivative computation
!            !! store indices of GLL points needed for FD derivatives in iglob
!            allocate(regularization_fd(iglob_index)%iglob_regular_point_to_use(in_index))

!            if (out_index == 0) then
!               out_index_to_allocate=1  !! just to avoid to allocate 0 sized array (to do : check if it is necessary)
!            else
!                out_index_to_allocate=out_index
!            endif

!            allocate(regularization_fd(iglob_index)%iglob_neighbo_point_to_use(out_index_to_allocate))
!            allocate(regularization_fd(iglob_index)%iglob_neighbo_rank_slice(out_index_to_allocate))

!            !! initialozation
!            regularization_fd(iglob_index)%iglob_regular_point_to_use(:)=-1
!            regularization_fd(iglob_index)%iglob_neighbo_point_to_use(:)=-1
!            regularization_fd(iglob_index)%iglob_neighbo_rank_slice(:)=-1

!            !! memory usage
!            mem_tmp = mem_tmp+ 4*(in_index + out_index)
!            mem_used_for_reg = mem_used_for_reg + 4*(in_index + 2*out_index_to_allocate)
!            max_memory_used = max( mem_used_for_reg, max_memory_used)

!            igll_local=0
!            in_index=0
!            out_index=0
!            do
!               igll_local=igll_local+1
!               if ( list_of_gll(igll_local) < 0) exit

!               !! store GLL point in my MPI slice
!               if ( rank_of_gll(igll_local) == myrank) then
!                  in_index=in_index+1
!                  regularization_fd(iglob_index)%iglob_regular_point_to_use(in_index)=list_of_gll(igll_local)
!               endif

!               !! store GLL point from neigbhour MPI slices
!               if ( rank_of_gll(igll_local) /= myrank .and. rank_of_gll(igll_local) > -1) then
!                  out_index=out_index+1
!                  regularization_fd(iglob_index)%iglob_neighbo_point_to_use(out_index)=list_of_gll(igll_local)
!                  regularization_fd(iglob_index)%iglob_neighbo_rank_slice(out_index)=rank_of_gll(igll_local)
!                  number_points_to_receive_from_neig(rank_of_gll(igll_local)) = &
!                       number_points_to_receive_from_neig(rank_of_gll(igll_local)) + 1
!               endif
!            enddo

!         endif
!      enddo

!   enddo

!   !! determine list of iglob to receive from different MPI slices -----------------------------------
!   allocate(list_iglob_neig(0:NPROC-1))
!   mem_used_for_reg = mem_used_for_reg  + 4 * NPROC
!   max_memory_used = max(mem_used_for_reg, max_memory_used)

!   deallocate(list_of_gll)
!   mem_used_for_reg = mem_used_for_reg  - 4 * MAX_POINTS_ALLOWED

!   do irank = 0, NPROC-1
!      list_iglob_neig(irank)%Niglob=0
!      Npts = number_points_to_receive_from_neig(irank)
!      if (Npts > 0) then
!         allocate(list_of_gll(Npts))
!         list_of_gll(:)=-1
!         out_index=0
!         do iglob_index=1, nused_iglob
!            do igll=1, regularization_fd(iglob_index)%nNei
!               if (regularization_fd(iglob_index)%iglob_neighbo_rank_slice(igll) == irank) then
!                  iglob = regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
!                  call add_or_not_one_more(list_of_gll, iglob, Npts, out_index)
!               endif
!            enddo
!         enddo
!         !! store list of GLL to receive from irank
!         list_iglob_neig(irank)%Niglob=out_index
!         allocate(list_iglob_neig(irank)%list_of_iglob(out_index))
!         list_iglob_neig(irank)%list_of_iglob(:)=list_of_gll(1:out_index)
!         deallocate(list_of_gll)
!         mem_used_for_reg = mem_used_for_reg  + 4 * out_index
!      endif
!   enddo

!   allocate(list_iglob_to_send(0:NPROC-1))
!   allocate(Npts_recv_dummy(0:NPROC-1))
!   allocate(request_recv_scalar(0:NPROC-1),request_send_scalar(0:NPROC-1))
!   allocate(sendtag(0:NPROC-1),recvtag(0:NPROC-1))

!   mem_used_for_reg = mem_used_for_reg + 6*4*NPROC
!   max_memory_used = max(mem_used_for_reg, max_memory_used)

!   Npts_recv_dummy(:)=0
!   do irank=0,NPROC-1
!      sendtag(irank)=irank
!      recvtag(irank)=irank
!   enddo

!   !! determine list of iglob to send to different MPI slices ---------------------------------------
!   do irank0 = 0, NPROC-1
!      do irank1 = irank0+1, NPROC-1

!         if (xcomm(irank0,irank1) > 0) then

!             if (myrank == irank0) then

!                 Npts_send_dummy(1)=list_iglob_neig(irank1)%Niglob

!                 if (DEBUG_MODE) write(IIDD,*) myrank, 'send ', Npts_send_dummy(1), ' to ', irank1

!                 call send_i(Npts_send_dummy,1,irank1,sendtag1)
!                 call recv_i(Npts_recv_dummy(irank1),1,irank1,recvtag1)

!                 if (DEBUG_MODE)  write(IIDD,*) myrank, 'received ', Npts_recv_dummy(irank1), ' from ', irank1

!                 allocate(list_of_gll_to_send(Npts_send_dummy(1)))
!                 allocate(list_of_gll_to_recv(Npts_recv_dummy(irank1)))

!                 list_of_gll_to_send(:)=list_iglob_neig(irank1)%list_of_iglob(:)

!                 call send_i(list_of_gll_to_send, Npts_send_dummy(1), irank1, sendtag1)
!                 call recv_i(list_of_gll_to_recv, Npts_recv_dummy(irank1), irank1, recvtag1)

!                 allocate(list_iglob_to_send(irank1)%list_of_iglob(Npts_recv_dummy(irank1)))

!                 list_iglob_to_send(irank1)%Niglob=Npts_recv_dummy(irank1)

!                 list_iglob_to_send(irank1)%list_of_iglob(:)=list_of_gll_to_recv(:)  !! irank0 recoit liste de irank1

!                 deallocate(list_of_gll_to_send)
!                 deallocate(list_of_gll_to_recv)
!             endif

!             if (myrank == irank1) then

!                call recv_i(Npts_recv_dummy(irank0),1,irank0,recvtag0)

!                if (DEBUG_MODE)  write(IIDD,*) myrank, 'recevied ', Npts_recv_dummy(irank0), ' from ', irank0

!                Npts_send_dummy(1)=list_iglob_neig(irank0)%Niglob

!                if (DEBUG_MODE) write(IIDD,*) myrank, 'send ', Npts_send_dummy(1), ' to ', irank0

!                call send_i(Npts_send_dummy, 1, irank0, sendtag0)

!                allocate(list_of_gll_to_send(Npts_send_dummy(1)))
!                allocate(list_of_gll_to_recv(Npts_recv_dummy(irank0)))

!                call recv_i(list_of_gll_to_recv, Npts_recv_dummy(irank0), irank0, recvtag0)


!! pb avec cette allocation !!!
!                allocate(list_iglob_to_send(irank0)%list_of_iglob(Npts_recv_dummy(irank0)))

!                list_iglob_to_send(irank0)%Niglob=Npts_recv_dummy(irank0)
!                list_iglob_to_send(irank0)%list_of_iglob(:)=list_of_gll_to_recv(:)    !! irank1 recoit liste de irank0

!                list_of_gll_to_send(:)=list_iglob_neig(irank0)%list_of_iglob(:)
!                call send_i(list_of_gll_to_send, Npts_send_dummy(1), irank0, sendtag0)

!                deallocate(list_of_gll_to_send)
!                deallocate(list_of_gll_to_recv)
!             endif

!          endif
!          call synchronize_all()
!       enddo
!    enddo


!    !! modification of regularization_fd in order to be able to find directly the index in received buffer
!    indx_recv(0)=0
!    do irank=1,NPROC
!       indx_recv(irank)=indx_recv(irank-1) + list_iglob_neig(irank-1)%Niglob
!    enddo

!    do iglob_index = 1, nused_iglob
!       if ( regularization_fd(iglob_index)%nNei > 0) then
!          !if (DEBUG_MODE) write(IIDD,*) iglob_index
!          do igll = 1, regularization_fd(iglob_index)%nNei
!             iglob_to_find =  regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
!             irank_to_look =  regularization_fd(iglob_index)%iglob_neighbo_rank_slice(igll)
!             find_iglob=.false.
!             do i=1, list_iglob_neig(irank_to_look)%Niglob
!                if (iglob_to_find == list_iglob_neig(irank_to_look)%list_of_iglob(i)) then
!                   find_iglob=.true.
!                   exit
!                endif
!             enddo
!             if (find_iglob) then
!                regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)=i+indx_recv(irank_to_look)
!             else
!                write(*,*) ' ABORT:', myrank, ' not find needed iglob ',iglob_to_find, ' in rank ',irank_to_look
!                stop
!             endif
!          enddo
!       endif
!    enddo

!    !! finalize structure for communication -------------------------------------------------------
!    allocate(struct_comm(0:NPROC-1))
!    do irank = 0, NPROC-1
!       struct_comm(irank)%ns=list_iglob_to_send(irank)%Niglob
!       struct_comm(irank)%nr=list_iglob_neig(irank)%Niglob
!       struct_comm(irank)%ibegin=indx_recv(irank)
!       if (list_iglob_to_send(irank)%Niglob > 0) then  ! send to irank
!          allocate(struct_comm(irank)%array_to_send(struct_comm(irank)%ns))
!          allocate(struct_comm(irank)%iglob_to_send(struct_comm(irank)%ns))
!          struct_comm(irank)%iglob_to_send(:)=list_iglob_to_send(irank)%list_of_iglob(:)
!       endif
!       if (list_iglob_neig(irank)%Niglob > 0) then ! receive from irank
!          allocate(struct_comm(irank)%array_to_recv(list_iglob_neig(irank)%Niglob))
!       endif
!    enddo

!    !! deallocate all temporary arrays ---------
!    deallocate(EtoE)
!    deallocate(EtoMPI)
!    mem_used_for_reg = mem_used_for_reg - 2 *4 * NSPEC_AB * MAX_NEIG

!   !! -----------------------------------------------------------------------------------------------

!   if (DEBUG_MODE) then
!      write(IIDD,*) 'NGLOB_AB ', NGLOB_AB
!      write(IIDD,*)
!      iglob_index=990
!      write(IIDD, *)
!      write(IIDD, *) ' iglob ',   regularization_fd(iglob_index)%iglob
!      write(IIDD, *) ' Niglobs ', regularization_fd(iglob_index)%nReg, regularization_fd(iglob_index)%nNei
!      do igll = 1, regularization_fd(iglob_index)%nReg
!         write(IIDD,*) regularization_fd(iglob_index)%iglob_regular_point_to_use(igll)
!      enddo
!      do igll = 1, regularization_fd(iglob_index)%nNei
!         write(IIDD,*) regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
!      enddo

!      write(IIDD,*)
!      write(IIDD,*)
!      write(IIDD,*) '               max memory  : ' ,max_memory_used, max_memory_used / 1000000, ' Mo'
!      write(IIDD,*) '           current memory  : ' ,mem_used_for_reg , mem_used_for_reg / 1000000, ' Mo'
!      write(IIDD,*) '         memory for grid  : ' ,CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*NSPEC_AB, &
!           CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*NSPEC_AB / 1000000, ' Mo'
!      write(IIDD,*) '         memory for field : ' ,CUSTOM_REAL*NGLOB_AB, CUSTOM_REAL*NGLOB_AB / 1000000 , ' Mo '
!      write(IIDD,*) ' memory used for structure : ', mem_tmp, &
!           mem_tmp/1000000, ' Mo'
!   endif
!   !! end of  deternine point in disk for each iglob------------------------------


!  end subroutine prepare_MPI_regularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  find all points involved in DF derivative around a central GLL point
!---------------------------------------------------------------------------------
!##################################################################################################################################

!  subroutine find_all_glls_in_disk(xp, yp, zp, list_of_elements, rank_of_elements, list_of_gll, rank_of_gll, in_index, out_index)

!    implicit none

!    real(kind=CUSTOM_REAL)                 :: xp, yp, zp
!    real(kind=CUSTOM_REAL)                 :: xs, ys, zs

!    integer, dimension(:),   allocatable   :: list_of_elements, rank_of_elements
!    integer, dimension(:),   allocatable   :: list_of_gll, rank_of_gll

!    real(kind=CUSTOM_REAL)                 :: dist_sq
!    integer                                :: in_index, out_index
!    integer                                :: irank, current_index
!    integer                                :: ispec, ie, iglob
!    integer                                :: i, j, k

!    current_index=0
!    in_index=0
!    out_index=0
!    list_of_gll(:)=-1
!    rank_of_gll(:)=-1
!    do ie=1, MAX_NEIG

!       ispec = list_of_elements(ie)

!       if (ispec <= 0) return

!       irank = rank_of_elements(ie)

!       if (irank /= myrank) then
!          do k=1,NGLLZ
!             do j=1,NGLLY
!                do i=1,NGLLX

!                   xs=xgrid(i,j,k,ispec)
!                   ys=ygrid(i,j,k,ispec)
!                   zs=zgrid(i,j,k,ispec)

!                   dist_sq = sqrt((xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2)

!                   if (dist_sq <= radius_disk) then
!                      current_index = current_index + 1
!                      out_index=out_index+1
!                      list_of_gll(current_index)=iglob_store(i,j,k,ispec)
!                      rank_of_gll(current_index)=irank
!                   endif

!                enddo
!             enddo
!          enddo

!       else

!          do k=1,NGLLZ
!             do j=1,NGLLY
!                do i=1,NGLLX

!                   iglob = ibool(i,j,k,ispec)
!                   xs = xstore(iglob)
!                   ys = ystore(iglob)
!                   zs = zstore(iglob)

!                   dist_sq = sqrt((xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2)

!                   if (dist_sq <= radius_disk) then
!                      in_index=in_index+1
!                      current_index = current_index + 1
!                      list_of_gll(current_index)=iglob
!                      rank_of_gll(current_index)=myrank
!                   endif

!                enddo
!             enddo
!          enddo

!       endif

!    enddo

!  end subroutine find_all_glls_in_disk

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  find list of neighbor elements for a given GLL point
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine find_list_of_element_adj(xp,yp,zp, list_of_all_elements, list_of_elements, ne_stored)

!   implicit none

!   real(kind=CUSTOM_REAL)             :: xp, yp, zp
!   real(kind=CUSTOM_REAL)             :: xs, ys, zs
!   integer                            :: ne_stored
!   integer, dimension(:), allocatable :: list_of_all_elements, list_of_elements
!   integer                            :: i,j,k,ie, ispec, iglob
!   real(kind=CUSTOM_REAL)             :: dist_sq, xtol=1e-6

!   do ie=1, MAX_NEIG

!      ispec = list_of_all_elements(ie)

!      if (ispec <= 0) return

!      if (ispec > NSPEC_AB) then

!         do k=1,NGLLZ,NGLLZ-1
!            do j=1,NGLLY,NGLLY-1
!               do i=1,NGLLX,NGLLX-1

!                  xs=xgrid(i,j,k,ispec-NSPEC_AB)
!                  ys=ygrid(i,j,k,ispec-NSPEC_AB)
!                  zs=zgrid(i,j,k,ispec-NSPEC_AB)

!                  dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2

!                  if (dist_sq <= xtol) then
!                     call add_or_not_one_more(list_of_elements, ispec, MAX_NEIG, ne_stored)
!                  endif
!               enddo
!            enddo
!         enddo

!      else

!         do k=1,NGLLZ,NGLLZ-1
!            do j=1,NGLLY,NGLLY-1
!               do i=1,NGLLX,NGLLX-1

!                  iglob = ibool(i,j,k,ispec)
!                  xs = xstore(iglob)
!                  ys = ystore(iglob)
!                  zs = zstore(iglob)

!                  dist_sq = (xp - xs)**2 + (yp - ys)**2 + (zp - zs)**2

!                  if (dist_sq <= xtol) call add_or_not_one_more(list_of_elements, ispec, MAX_NEIG, ne_stored)

!               enddo
!            enddo
!         enddo

!      endif

!   enddo

! end subroutine find_list_of_element_adj

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  if 'is' not appears in list 'it' then add it
!---------------------------------------------------------------------------------
!##################################################################################################################################

!  subroutine add_or_not_one_more(it, is, n, k)

!    implicit none

!    integer                            :: is, n, k
!    integer, dimension(:), allocatable :: it
!    integer i

!    do i=1, n
!       if (it(i) == is ) return
!    enddo

!    k = k + 1
!    it(k) = is

!  end subroutine add_or_not_one_more

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  count the number of GLL on element faces
!---------------------------------------------------------------------------------
!##################################################################################################################################

!  subroutine count_gll_in_face(nb_gll_in_face)

!    implicit none

!    integer                              :: nb_gll_in_face
!    integer                              :: i, j, k, ifac_loc

!    nb_gll_in_face=0
!    do k=1, NGLLZ
!       do j=1,NGLLY
!          do i=1,NGLLX

!             ifac_loc=0
!             !! warning priority depends on this order
!             if (i == 1)     ifac_loc=5
!             if (i == NGLLX) ifac_loc=6
!             if (j == 1)     ifac_loc=3
!             if (j == NGLLY) ifac_loc=4
!             if (k == 1)     ifac_loc=1
!             if (k == NGLLZ) ifac_loc=2
!             if (ifac_loc > 0) nb_gll_in_face=nb_gll_in_face+1

!          enddo
!       enddo
!    enddo

!  end subroutine count_gll_in_face

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  make list of all GLL on element faces
!---------------------------------------------------------------------------------
!##################################################################################################################################

!  subroutine store_list_gll_in_face(list_of_gll_in_face)

!    implicit none

!    integer, dimension(:,:), allocatable :: list_of_gll_in_face
!    integer                              :: i, j, k, kgll, ifac_loc

!    !! just to be able to loop over GLL on surface but only
!    !! once per GLL (avoid duplication in corners and edges)
!    kgll=0
!    do k=1, NGLLZ
!       do j=1,NGLLY
!          do i=1,NGLLX

!             ifac_loc=0
!             !! warning priority depends on this order
!             if (i == 1)     ifac_loc=5
!             if (i == NGLLX) ifac_loc=6
!             if (j == 1)     ifac_loc=3
!             if (j == NGLLY) ifac_loc=4
!             if (k == 1)     ifac_loc=1
!             if (k == NGLLZ) ifac_loc=2

!             if (ifac_loc > 0) then
!                kgll = kgll + 1
!                list_of_gll_in_face(1,kgll)=i
!                list_of_gll_in_face(2,kgll)=j
!                list_of_gll_in_face(3,kgll)=k
!                list_of_gll_in_face(4,kgll)=ifac_loc
!             endif

!          enddo
!       enddo
!    enddo

!  end subroutine store_list_gll_in_face

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  test if 'is' is in list 'it', add it if not appears and return logical flag nuy (true if no appears)
!---------------------------------------------------------------------------------
!##################################################################################################################################

! subroutine iglob_not_used(it, is, n, k, nuy)

!    implicit none

!    integer                            :: is, n, k
!    integer, dimension(:), allocatable :: it
!    integer i
!    logical :: nuy
!    nuy=.false.
!    do i=1, n
!       if (it(i) == is ) return
!    enddo
!    k = k + 1
!    it(k) = is
!    nuy=.true.

!  end subroutine iglob_not_used

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!##################################################################################################################################
  !! VM VM :: COPIED / PASTED subroutine from genetrate_databases/read_partition_files.f90
  !! TOFIX : find a way to use the orignal file to avoid this copy/paste

  subroutine read_partition_files

! reads in proc***_Databases files

  implicit none

  integer :: ispec, inode, num_interface, ie, imat
  integer :: num_xmin,num_xmax,num_ymin,num_ymax,num_top,num_bottom,num
  integer :: num_cpml
  integer :: num_moho
  integer :: i,j
  integer :: ier
  integer :: ispec2D
  integer :: dummy_node
  integer :: dummy_elmnt

  integer :: nnodes_ext_mesh, nelmnts_ext_mesh

  ! C-PML absorbing boundary conditions

  ! local number of C-PML spectral elements
  integer :: nspec_cpml

  ! global number of C-PML spectral elements
  integer :: nspec_cpml_tot

  ! C-PML spectral elements global indexing
  integer, dimension(:), allocatable :: CPML_to_spec

  ! C-PML regions
  integer, dimension(:), allocatable :: CPML_regions

  ! mask of C-PML elements for the global mesh
  logical, dimension(:), allocatable :: is_CPML

  ! moho (optional)
  integer :: nspec2D_moho_ext
  integer, dimension(:), allocatable  :: ibelm_moho
  integer, dimension(:,:), allocatable  :: nodes_ibelm_moho

  character(len=MAX_STRING_LEN) :: LOCAL_PATH_FOR_READING, prname_for_reading
  character(len=8)              :: path_to_add

  !! all group will read the partition file
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
     write(path_to_add,"('run',i4.4,'/')") 1
     LOCAL_PATH_FOR_READING=trim(path_to_add)//trim(LOCAL_PATH(9:len_trim(LOCAL_PATH)))
  else
     LOCAL_PATH_FOR_READING=LOCAL_PATH
  endif

! read databases about external mesh simulation
! global node coordinates
  call create_name_database(prname_for_reading,myrank,LOCAL_PATH_FOR_READING)
  open(unit=IIN,file=prname_for_reading(1:len_trim(prname_for_reading))//'Database', &
        status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'rank ',myrank,'- Error opening file: ', &
         prname_for_reading(1:len_trim(prname_for_reading))//'Database'
    print *,'please make sure file exists'
    call exit_mpi(myrank,'Error opening database file')
  endif
  read(IIN) nnodes_ext_mesh

  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 164')
  if (ier /= 0) stop 'Error allocating array nodes_coords_ext_mesh'

  do inode = 1, nnodes_ext_mesh
     read(IIN) dummy_node, nodes_coords_ext_mesh(1,inode), nodes_coords_ext_mesh(2,inode), &
                nodes_coords_ext_mesh(3,inode)
  enddo

  call sum_all_i(nnodes_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'external mesh points : ',num
  endif
  call synchronize_all()

! read physical properties of the materials
! added poroelastic properties and filled with 0 the last 10 entries for elastic/acoustic
  read(IIN) nmat_ext_mesh, nundefMat_ext_mesh

  allocate(mat_prop(17,nmat_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 165')
  if (ier /= 0) stop 'Error allocating array mat_prop'
  mat_prop(:,:) = 0.d0

  do imat = 1, nmat_ext_mesh
     ! (visco)elastic or acoustic format:
     ! #(1) rho   #(2) vp  #(3) vs #(4) Q_kappa  #(5) Q_mu  #(6) anisotropy_flag  #(7) material_domain_id
     ! and remaining entries are filled with zeros.
     ! Q_kappa is not stored next to Q_mu for historical reasons, because it was added later.
     !
     ! poroelastic format:  rhos,rhof,phi,tort,eta,0,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
     read(IIN) mat_prop(1,imat),  mat_prop(2,imat),  mat_prop(3,imat), &
               mat_prop(4,imat),  mat_prop(5,imat),  mat_prop(6,imat), &
               mat_prop(7,imat),  mat_prop(8,imat),  mat_prop(9,imat), &
               mat_prop(10,imat), mat_prop(11,imat), mat_prop(12,imat), &
               mat_prop(13,imat), mat_prop(14,imat), mat_prop(15,imat), &
               mat_prop(16,imat), mat_prop(17,imat)
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'defined materials    : ',nmat_ext_mesh
  endif
  call synchronize_all()

  allocate(undef_mat_prop(6,nundefMat_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 166')
  if (ier /= 0) stop 'Error allocating array undef_mat_prop'
  do imat = 1, nundefMat_ext_mesh
     ! format example tomography:
     ! #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
     ! e.g.: -1 tomography elastic tomography_model.xyz 0 2
     ! format example interface:
     ! e.g.: -1 interface 14 15 0 2
     read(IIN) undef_mat_prop(1,imat),undef_mat_prop(2,imat),undef_mat_prop(3,imat),undef_mat_prop(4,imat), &
          undef_mat_prop(5,imat), undef_mat_prop(6,imat)
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'undefined materials  : ',nundefMat_ext_mesh
  endif
  call synchronize_all()

! element indexing
  read(IIN) nelmnts_ext_mesh
  allocate(elmnts_ext_mesh(NGNOD,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 167')
  if (ier /= 0) stop 'Error allocating array elmnts_ext_mesh'
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 168')
  if (ier /= 0) stop 'Error allocating array mat_ext_mesh'

  ! reads in material association for each spectral element and corner node indices
  do ispec = 1, nelmnts_ext_mesh
     ! format:
     ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
     ! or
     ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id27
     read(IIN) dummy_elmnt,mat_ext_mesh(1,ispec),mat_ext_mesh(2,ispec),(elmnts_ext_mesh(j,ispec),j=1,NGNOD)

     ! check debug
     if (dummy_elmnt /= ispec) stop 'Error ispec order in materials file'

  enddo
  NSPEC_AB = nelmnts_ext_mesh

  call sum_all_i(nspec_ab,num)
  if (myrank == 0) then
    write(IMAIN,*) 'total number of spectral elements: ',num
  endif
  call synchronize_all()

! reads absorbing/free-surface boundaries
  read(IIN) boundary_number ,nspec2D_xmin
  if (boundary_number /= 1) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_xmax
  if (boundary_number /= 2) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_ymin
  if (boundary_number /= 3) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_ymax
  if (boundary_number /= 4) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_bottom_ext
  if (boundary_number /= 5) stop "Error : invalid database file"

  read(IIN) boundary_number ,nspec2D_top_ext
  if (boundary_number /= 6) stop "Error : invalid database file"

  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext

  allocate(ibelm_xmin(nspec2D_xmin),nodes_ibelm_xmin(NGNOD2D,nspec2D_xmin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 169')
  if (ier /= 0) stop 'Error allocating array ibelm_xmin etc.'
  do ispec2D = 1,nspec2D_xmin
     read(IIN) ibelm_xmin(ispec2D),(nodes_ibelm_xmin(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_xmax(nspec2D_xmax),nodes_ibelm_xmax(NGNOD2D,nspec2D_xmax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 170')
  if (ier /= 0) stop 'Error allocating array ibelm_xmax etc.'
  do ispec2D = 1,nspec2D_xmax
     read(IIN) ibelm_xmax(ispec2D),(nodes_ibelm_xmax(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_ymin(nspec2D_ymin),nodes_ibelm_ymin(NGNOD2D,nspec2D_ymin),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 171')
  if (ier /= 0) stop 'Error allocating array ibelm_ymin'
  do ispec2D = 1,nspec2D_ymin
     read(IIN) ibelm_ymin(ispec2D),(nodes_ibelm_ymin(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_ymax(nspec2D_ymax),nodes_ibelm_ymax(NGNOD2D,nspec2D_ymax),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 172')
  if (ier /= 0) stop 'Error allocating array ibelm_ymax etc.'
  do ispec2D = 1,nspec2D_ymax
     read(IIN) ibelm_ymax(ispec2D),(nodes_ibelm_ymax(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_bottom(nspec2D_bottom_ext),nodes_ibelm_bottom(NGNOD2D,nspec2D_bottom_ext),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 173')
  if (ier /= 0) stop 'Error allocating array ibelm_bottom etc.'
  do ispec2D = 1,nspec2D_bottom_ext
     read(IIN) ibelm_bottom(ispec2D),(nodes_ibelm_bottom(j,ispec2D),j=1,NGNOD2D)
  enddo

  allocate(ibelm_top(nspec2D_top_ext),nodes_ibelm_top(NGNOD2D,nspec2D_top_ext),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 174')
  if (ier /= 0) stop 'Error allocating array ibelm_top etc.'
  do ispec2D = 1,nspec2D_top_ext
     read(IIN) ibelm_top(ispec2D),(nodes_ibelm_top(j,ispec2D),j=1,NGNOD2D)
  enddo

  call sum_all_i(nspec2D_xmin,num_xmin)
  call sum_all_i(nspec2D_xmax,num_xmax)
  call sum_all_i(nspec2D_ymin,num_ymin)
  call sum_all_i(nspec2D_ymax,num_ymax)
  call sum_all_i(nspec2D_top_ext,num_top)
  call sum_all_i(nspec2D_bottom_ext,num_bottom)

  if (myrank == 0) then
    write(IMAIN,*) 'absorbing boundaries: '
    write(IMAIN,*) '  xmin,xmax : ',num_xmin,num_xmax
    write(IMAIN,*) '  ymin,ymax : ',num_ymin,num_ymax
    write(IMAIN,*) '  bottom,top: ',num_bottom,num_top
  endif
  call synchronize_all()

  ! reads number of C-PML elements in the global mesh
  nspec_cpml_tot = 0
  nspec_cpml = 0

  read(IIN) nspec_cpml_tot
  if (myrank == 0) then
     write(IMAIN,*) 'total number of C-PML elements in the global mesh: ',nspec_cpml_tot
  endif
  call synchronize_all()

  ! reads mask of C-PML elements for all elements in this partition
  allocate(is_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 177')
  if (ier /= 0) stop 'Error allocating array is_CPML'
  is_CPML(:) = .false.

  if (nspec_cpml_tot > 0) then
     ! reads number of C-PML elements in this partition
     read(IIN) nspec_cpml

     if (myrank == 0) then
        write(IMAIN,*) '  number of C-PML spectral elements in this partition: ',nspec_cpml
     endif
     call synchronize_all()

     call sum_all_i(nspec_cpml,num_cpml)

     ! checks that the sum of C-PML elements over all partitions is correct
     if (myrank == 0 .and. nspec_cpml_tot /= num_cpml) stop 'Error while summing C-PML elements over all partitions'

     ! reads C-PML regions and C-PML spectral elements global indexing
     allocate(CPML_to_spec(nspec_cpml),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 175')
     if (ier /= 0) stop 'Error allocating array CPML_to_spec'
     allocate(CPML_regions(nspec_cpml),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 176')
     if (ier /= 0) stop 'Error allocating array CPML_regions'

     do i=1,nspec_cpml
        ! #id_cpml_regions = 1 : X_surface C-PML
        ! #id_cpml_regions = 2 : Y_surface C-PML
        ! #id_cpml_regions = 3 : Z_surface C-PML
        ! #id_cpml_regions = 4 : XY_edge C-PML
        ! #id_cpml_regions = 5 : XZ_edge C-PML
        ! #id_cpml_regions = 6 : YZ_edge C-PML
        ! #id_cpml_regions = 7 : XYZ_corner C-PML
        !
        ! format: #id_cpml_element #id_cpml_regions
        read(IIN) CPML_to_spec(i), CPML_regions(i)
     enddo

     do i=1,NSPEC_AB
        read(IIN) is_CPML(i)
     enddo
  endif

! MPI interfaces between different partitions
  if (NPROC > 1) then
    ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
    read(IIN) num_interfaces_ext_mesh, max_interface_size_ext_mesh
  else
    num_interfaces_ext_mesh = 0
    max_interface_size_ext_mesh = 0
  endif

  ! allocates interfaces
  allocate(my_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 178')
  if (ier /= 0) stop 'Error allocating array my_neighbors_ext_mesh'
  allocate(my_nelmnts_neighbors_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 179')
  if (ier /= 0) stop 'Error allocating array my_nelmnts_neighbors_ext_mesh'
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 180')
  if (ier /= 0) stop 'Error allocating array my_interfaces_ext_mesh'
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 181')
  if (ier /= 0) stop 'Error allocating array ibool_interfaces_ext_mesh'
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 182')
  if (ier /= 0) stop 'Error allocating array nibool_interfaces_ext_mesh'

  ! loops over MPI interfaces with other partitions
  do num_interface = 1, num_interfaces_ext_mesh
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    read(IIN) my_neighbors_ext_mesh(num_interface), my_nelmnts_neighbors_ext_mesh(num_interface)

    ! loops over interface elements
    do ie = 1, my_nelmnts_neighbors_ext_mesh(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5)...
      !
      ! interface types:
      !     1  -  corner point only
      !     2  -  element edge
      !     4  -  element face
      read(IIN) my_interfaces_ext_mesh(1,ie,num_interface), my_interfaces_ext_mesh(2,ie,num_interface), &
                  my_interfaces_ext_mesh(3,ie,num_interface), my_interfaces_ext_mesh(4,ie,num_interface), &
                  my_interfaces_ext_mesh(5,ie,num_interface), my_interfaces_ext_mesh(6,ie,num_interface)
    enddo
  enddo

  call sum_all_i(num_interfaces_ext_mesh,num)
  if (myrank == 0) then
    write(IMAIN,*) 'number of MPI partition interfaces: ',num
  endif
  call synchronize_all()

  ! optional moho
  if (SAVE_MOHO_MESH) then
     ! checks if additional line exists
     read(IIN,iostat=ier) boundary_number,nspec2D_moho_ext
     if (ier /= 0) then
        ! no moho informations given
        nspec2D_moho_ext = 0
        boundary_number = 7
     endif
     if (boundary_number /= 7) stop "Error invalid database file"

     ! checks total number of elements
     call sum_all_i(nspec2D_moho_ext,num_moho)
     if (num_moho == 0) call exit_mpi(myrank,'Error no moho mesh in database')

     ! reads in element informations
     allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 183')
     if (ier /= 0) stop 'Error allocating array ibelm_moho etc.'
     do ispec2D = 1,nspec2D_moho_ext
        ! format: #element_id #node_id1 #node_id2 #node_id3 #node_id4
        read(IIN) ibelm_moho(ispec2D),(nodes_ibelm_moho(j,ispec2D),j=1,NGNOD2D)
     enddo

     ! user output
     if (myrank == 0) then
        write(IMAIN,*) 'moho surfaces: ',num_moho
     endif
     call synchronize_all()
  else
     ! allocate dummy array
     nspec2D_moho_ext = 0
     allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(NGNOD2D,nspec2D_moho_ext),stat=ier)
     if (ier /= 0) call exit_MPI_without_rank('error allocating array 184')
     if (ier /= 0) stop 'Error allocating dumy array ibelm_moho etc.'
  endif

  close(IIN)

  end subroutine read_partition_files



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!!!!#####################################################################################################################
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  check the accuracy of numerical laplacian by comparison with analytical cosine functions
!---------------------------------------------------------------------------------
!##################################################################################################################################
  !! check the accuracy of laplacian using cosinre function define on sepcfem
  !! mesh: the wave length of cosine is :  2x, 4x, 8x the biggest element size
  !!

!!!!!!!!!!  subroutine check_regularization_on_mesh(regularization_fd)
  subroutine check_regularization_on_mesh()

!!!!!!!!!!    type(regul), dimension(:), allocatable, intent(in)         :: regularization_fd
    real(kind=CUSTOM_REAL)                                     :: distance_min_glob,distance_max_glob
    real(kind=CUSTOM_REAL)                                     :: elemsize_min_glob,elemsize_max_glob
    real(kind=CUSTOM_REAL)                                     :: x_min_glob,x_max_glob
    real(kind=CUSTOM_REAL)                                     :: y_min_glob,y_max_glob
    real(kind=CUSTOM_REAL)                                     :: z_min_glob,z_max_glob
    integer,                 dimension(NGNOD)                  :: iaddx,iaddy,iaddz

    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable    :: field
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable    :: laplacian_of_field
    real(kind=CUSTOM_REAL)                                     :: value_to_test(5)
    real(kind=CUSTOM_REAL)                                     :: length
    integer                                                    :: itest, Nb_test

  integer :: ier

    value_to_test(1)=2.
    value_to_test(2)=4.
    value_to_test(3)=8.
    value_to_test(4)=16.
    value_to_test(5)=32.

    ! get mesh properties
    call usual_hex_nodes(NGNOD,iaddx,iaddy,iaddz)
    call check_mesh_distances(myrank,NSPEC_AB,NGLOB_AB,ibool,xstore,ystore,zstore, &
         x_min_glob,x_max_glob,y_min_glob,y_max_glob,z_min_glob,z_max_glob, &
         elemsize_min_glob,elemsize_max_glob, &
         distance_min_glob,distance_max_glob)

    allocate(field(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 185')
    allocate(laplacian_of_field(NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 186')

    if (myrank == 0) then
       write(INVERSE_LOG_FILE,*)
       write(INVERSE_LOG_FILE,*) ' PERFORMING ACCURACY TEST FOR FD REGULARIZATION '
       write(INVERSE_LOG_FILE,*)
    endif

    Nb_test=5
    do itest=1, Nb_test
       length = value_to_test(itest) * elemsize_max_glob
       lambda = 2* 3.1459265359/length

       if (myrank == 0) write(INVERSE_LOG_FILE,*) ' check cosine function numerical laplacian :  wavelength  :', length

       call compute_cosine_field(lambda, field, laplacian_of_field)


!!$       current_model_vp(:,:,:,:) = field(:,:,:,:)
!!$       current_model_vs(:,:,:,:) = field(:,:,:,:)
!!$       current_model_rh(:,:,:,:) = field(:,:,:,:)


       !call compute_laplacian2_vp_vs_rho(Lap1_vp, Lap1_vs, Lap1_rh, Lap2_vp, Lap2_vs, Lap2_rh, &
       !regularization_fd)

       !! penalaty : norm L2 of gradient
!!$       call compute_grad_laplacian_vp_vs_rho(nGrad_vp, nGrad_vs, nGrad_rh, Lap1_vp, Lap1_vs, Lap1_rh) !!!!! , regularization_fd )
!!$
!!$       call check_field_accuracy(itest, Lap1_vp, laplacian_of_field)

    enddo

  end subroutine check_regularization_on_mesh

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  check the accuracy of numerical laplacian
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine check_field_accuracy(itest, numerical_laplacian_of_field, laplacian_of_field)

!   implicit none

!   integer,                                                    intent(in) :: itest
!   real(kind=CUSTOM_REAL),  dimension(:,:,:,:),   allocatable, intent(in) :: laplacian_of_field
!   real(kind=CUSTOM_REAL),  dimension(:,:,:,:),   allocatable, intent(in) :: numerical_laplacian_of_field

!   real(kind=CUSTOM_REAL)                                                 :: L2_error, L1_error, Linf_error, max_val

    !! for DEBUG_MODE ----- can be removed --------------------------------------------------------------
!   character(len=256)                                                     :: path_file, name_file
    !! --------------------------------------------------------------------------------------------------

!   L2_error=0.
!   L1_error=0.
!   Linf_error=0.

!   max_val = maxval(abs (laplacian_of_field))

!   L2_error = L2_error + sum( numerical_laplacian_of_field(:,:,:,:) - laplacian_of_field(:,:,:,:) )**2
!   L1_error = L1_error + sum(abs( numerical_laplacian_of_field(:,:,:,:) - laplacian_of_field(:,:,:,:) ))
!   Linf_error = maxval(abs( numerical_laplacian_of_field(:,:,:,:) - laplacian_of_field(:,:,:,:) ))
!   !! todo manque les comms MPI pour les erreurs

!   if (myrank == 0) then
!      call flush(INVERSE_LOG_FILE)
!      write(INVERSE_LOG_FILE,*) 'ABSOLUTE ERRORS : L1, L2, Linf :', L1_error, L2_error, Linf_error
!      !write(INVERSE_LOG_FILE,*) 'RELATIVE ERRORS : L1, L2, Linf :', L1_error/ max_val, L2_error/ max_val, &
!      !     Linf_error/ max_val
!   endif

!   if (DEBUG_MODE) then !! ------ can be removed: from here  ----

!      !allocate(field_vtk(NGLLX, NGLLY, NGLLZ, NSPEC_AB))

!      !! display analytical laplacian
!      !do ispec = 1, NSPEC_AB
!      !   do k =1, NGLLZ
!      !      do j =1, NGLLY
!      !         do i=1, NGLLX
!      !            iglob = ibool(i,j,k,ispec)
!      !            field_vtk(i,j,k,ispec) = laplacian_of_field(iglob)
!      !         enddo
!      !      enddo
!      !   enddo
!      !enddo

!      path_file='OUTPUT_FILES/DATABASES_MPI/proc'
!      write(name_file,'(i6.6,a8,i2.2,a4)') myrank, '_laplac_',itest,'.bin'
!      path_file=(trim(path_file))//trim(name_file)
!      open(888,file=trim(path_file),form='unformatted')
!      write(888) laplacian_of_field
!      close(888)

!      !! display numerical laplacian
!      !do ispec = 1, NSPEC_AB
!      !   do k =1, NGLLZ
!      !      do j =1, NGLLY
!      !         do i=1, NGLLX
!      !            iglob = ibool(i,j,k,ispec)
!      !            field_vtk(i,j,k,ispec) = numerical_laplacian_of_field(iglob)
!      !         enddo
!      !      enddo
!      !   enddo
!      !enddo

!      path_file='OUTPUT_FILES/DATABASES_MPI/proc'
!      write(name_file,'(i6.6,a11,i2.2,a4)') myrank, '_munlaplac_',itest,'.bin'
!      path_file=(trim(path_file))//trim(name_file)
!      open(888,file=trim(path_file),form='unformatted')
!      write(888) numerical_laplacian_of_field
!      close(888)

!   endif  !! -------------------- to here ------------------------------

! end subroutine check_field_accuracy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute analytical cosine function in order to compare with the numerical one
!---------------------------------------------------------------------------------
!##################################################################################################################################
  subroutine compute_cosine_field(lambda, field, laplacian_of_field)

    implicit none
    real(kind=CUSTOM_REAL),                                     intent(in)    :: lambda
    real(kind=CUSTOM_REAL),  dimension(:,:,:,:),   allocatable, intent(inout) :: laplacian_of_field
    real(kind=CUSTOM_REAL),  dimension(:,:,:,:),   allocatable, intent(inout) :: field

    integer :: i, j, k, ispec, iglob

    !if (myrank==0)write (*,*) lambda
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1, NGLLX
                iglob = ibool(i,j,k,ispec)

                field(i,j,k,ispec) =  cos( lambda * xstore(iglob) ) * cos( lambda* ystore(iglob) ) * cos( lambda * zstore(iglob) ) &
                     / (3 * lambda * lambda)

!!$                x =  - sin( lambda * xstore(iglob) ) * &
!!$                       cos( lambda * ystore(iglob) ) * &
!!$                       cos( lambda * zstore(iglob) ) / (3 * lambda)
!!$
!!$                y = - cos( lambda * xstore(iglob) ) * &
!!$                      sin( lambda * ystore(iglob) ) * &
!!$                      cos( lambda * zstore(iglob) ) / (3 * lambda)
!!$                z = - cos( lambda * xstore(iglob) ) * &
!!$                      cos( lambda * ystore(iglob) ) * &
!!$                      sin( lambda * zstore(iglob) ) / (3 * lambda)

                laplacian_of_field(i,j,k,ispec) = &
                                         cos( lambda * xstore(iglob) ) * &
                                         cos( lambda * ystore(iglob) ) * &
                                         cos( lambda * zstore(iglob) )
!!$
!!$                laplacian_of_field(i,j,k,ispec) = sqrt(x**2 + y**2 +z**2)
                !field(i,j,k,ispec) = cos( lambda * zstore(iglob) ) / (lambda * lambda*lambda*lambda)
                !laplacian_of_field(i,j,k,ispec) = cos( lambda * zstore(iglob) ) * (lambda * lambda)
                !field(i,j,k,ispec) = ( -0.5 * ( ( zstore(iglob) + 30000.)  * lambda )**2)
                !field(i,j,k, ispec) = ((zstore(iglob) / 1000 ) ** 5 )/ (5*4*3*2)
                !laplacian_of_field(i,j,k,ispec) = zstore(iglob) /1000

             enddo
          enddo
       enddo
    enddo

  end subroutine compute_cosine_field

!###########################################################################################################################
!###########################################################################################################################
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute laplacian of a given field in all gll
!---------------------------------------------------------------------------------
!##################################################################################################################################
  subroutine compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)

    implicit none

    type(regul),             dimension(:),   allocatable, intent(in)    :: regularization_fd
    real(kind=CUSTOM_REAL),  dimension(:,:), allocatable, intent(in)    :: field
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable, intent(inout) :: numerical_laplacian_of_field

    real(kind=CUSTOM_REAL),  dimension(:),   allocatable                :: field_to_derivate
    real(kind=CUSTOM_REAL),  dimension(:,:), allocatable                :: Laplac_boundary, LapF

    integer :: ier

    !! the order below is important do not change it.
    allocate(Laplac_boundary(NDIM, Nb_iglob_on_faces), field_to_derivate(NGLOB_AB), LapF(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 187')

    !! 1/ compute FD derivatives
    call compute_laplacian_FD(Laplac_boundary, field, regularization_fd)

   !! 2/ compute deivatives inside elements by lagrange interpolation
    field_to_derivate(:)=field(1,:)
    call compute_laplac_lagrange(LapF, field_to_derivate)

    !! 3/ replace boundary elements
    numerical_laplacian_of_field(:) = LapF(1,:)
!!$    do iglob_index=1,  Nb_iglob_on_faces
!!$       iglob = regularization_fd(iglob_index)%iglob
!!$       numerical_laplacian_of_field(iglob) = Laplac_boundary(1, iglob_index)
!!$    enddo

    deallocate(Laplac_boundary, field_to_derivate, LapF)

  end subroutine compute_laplacian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute gradient norm and laplacian of a given field at all GLL points
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_grad_laplacian(field, laplacian_of_field, norm_grad_of_field) !!!!!!!!!!!!! , regularization_fd)

    implicit none

!!!!!!!!!!!!!    type(regul),             dimension(:),   allocatable, intent(in)    :: regularization_fd
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable, intent(in)    :: field
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable, intent(inout) :: laplacian_of_field
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable, intent(inout) :: norm_grad_of_field

    real(kind=CUSTOM_REAL),  dimension(:),   allocatable                :: field_to_derivate
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable                :: Laplac_boundary, nGrad_boundary
    real(kind=CUSTOM_REAL),  dimension(:),   allocatable                :: nGrad, LapF

  integer :: ier

    allocate(Laplac_boundary(Nb_iglob_on_faces), nGrad_boundary(Nb_iglob_on_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 188')
    allocate(field_to_derivate(NGLOB_AB), LapF(NGLOB_AB), nGrad(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 189')

    !! the order below is important, do not change it.

    !! 1/ compute FD derivatives
    !! NOT USE FOR NOW call compute_gradient_laplacian_FD(nGrad_boundary, Laplac_boundary, field, regularization_fd)

   !! 2/ compute deivatives inside elements by lagrange interpolation
    field_to_derivate(:)=field(:)
    call compute_grad_laplac_lagrange(nGrad, LapF, field_to_derivate)

    !! 3/ replace boundary elements
    laplacian_of_field(:) = LapF(:)
    norm_grad_of_field(:) = nGrad(:)
!!$    do iglob_index=1,  Nb_iglob_on_faces
!!$       iglob = regularization_fd(iglob_index)%iglob
!!$       laplacian_of_field(iglob) = Laplac_boundary(iglob_index)
!!$       norm_grad_of_field(iglob) = nGrad_boundary(iglob_index)
!!$    enddo

    deallocate(Laplac_boundary, nGrad_boundary, field_to_derivate, LapF, nGrad)

  end subroutine compute_grad_laplacian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------------------
! compute laplacian of laplacian of vp, vs and rho curent fields
!---------------------------------------------------------------------------------

  subroutine compute_grad_laplacian_vp_vs_rho(nGrad_vp, nGrad_vs, nGrad_rh, Lap1_vp, Lap1_vs, Lap1_rh) !!!!!!! , regularization_fd)

    implicit none

!!!!!!!!!    type(regul),            dimension(:),          allocatable, intent(in)    :: regularization_fd
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),    allocatable, intent(inout) :: nGrad_vp, nGrad_vs, nGrad_rh
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),    allocatable, intent(inout) :: Lap1_vp, Lap1_vs, Lap1_rh

    real(kind=CUSTOM_REAL),  dimension(:),         allocatable                :: field_to_derivate
    real(kind=CUSTOM_REAL),  dimension(:,:),       allocatable                :: field
    real(kind=CUSTOM_REAL),  dimension(:),         allocatable                :: laplacian_of_field, norm_grad_of_field
    real(kind=CUSTOM_REAL),  dimension(:),         allocatable                :: valence
    real(kind=CUSTOM_REAL)                                                    :: penalty
    integer                                                                   :: i,j,k,ispec,iglob

  integer :: ier

    allocate(field(NDIM,NGLOB_AB), field_to_derivate(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 190')
    allocate(laplacian_of_field(NGLOB_AB),  norm_grad_of_field(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 191')
    allocate(valence(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 192')


    field(:,:)=0._CUSTOM_REAL
    valence(:)=0._CUSTOM_REAL
    penalty=0._CUSTOM_REAL

    !! mean values of field on element boundary
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
!!$                field(1,iglob)=field(1,iglob) + current_model_vp(i,j,k,ispec)
!!$                field(2,iglob)=field(2,iglob) + current_model_vs(i,j,k,ispec)
!!$                field(3,iglob)=field(3,iglob) + current_model_rh(i,j,k,ispec)
                valence(iglob)=valence(iglob) + 1._CUSTOM_REAL
             enddo
          enddo
       enddo
    enddo

    field(1,:)=field(1,:)/valence(:)
    field(2,:)=field(2,:)/valence(:)
    field(3,:)=field(3,:)/valence(:)

    field_to_derivate(:)=field(1,:)
    call compute_grad_laplacian(field_to_derivate, laplacian_of_field, norm_grad_of_field) !!!!!!!!! , regularization_fd)
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                nGrad_vp(i,j,k,ispec)=norm_grad_of_field(iglob)
                Lap1_vp(i,j,k,ispec)=laplacian_of_field(iglob)
             enddo
          enddo
       enddo
    enddo

    field_to_derivate(:)=field(2,:)
    call compute_grad_laplacian(field_to_derivate, laplacian_of_field, norm_grad_of_field) !!!!!!!!! , regularization_fd)
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                nGrad_vs(i,j,k,ispec)=norm_grad_of_field(iglob)
                Lap1_vs(i,j,k,ispec)=laplacian_of_field(iglob)
             enddo
          enddo
       enddo
    enddo

    field_to_derivate(:)=field(3,:)
    call compute_grad_laplacian(field_to_derivate, laplacian_of_field, norm_grad_of_field) !!!!!!!!! , regularization_fd)
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                nGrad_rh(i,j,k,ispec)=norm_grad_of_field(iglob)
                Lap1_rh(i,j,k,ispec)=laplacian_of_field(iglob)
             enddo
          enddo
       enddo
    enddo

  end subroutine compute_grad_laplacian_vp_vs_rho


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------------------
! compute laplacian of laplacian of vp, vs and rho curent fields
!---------------------------------------------------------------------------------

  subroutine compute_laplacian2_vp_vs_rho(Lap1_vp, Lap1_vs, Lap1_rh, Lap2_vp, Lap2_vs, Lap2_rh) !!!!!! , regularization_fd)

    implicit none

!!!!!!    type(regul),            dimension(:),          allocatable, intent(in)    :: regularization_fd
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),    allocatable, intent(inout) :: Lap1_vp, Lap1_vs, Lap1_rh
    real(kind=CUSTOM_REAL), dimension(:,:,:,:),    allocatable, intent(inout) :: Lap2_vp, Lap2_vs, Lap2_rh

    real(kind=CUSTOM_REAL), dimension(:,:,:,:),    allocatable                :: field_to_derivate
    real(kind=CUSTOM_REAL),  dimension(:,:),       allocatable                :: field
    real(kind=CUSTOM_REAL),  dimension(:),         allocatable                :: numerical_laplacian_of_field
    real(kind=CUSTOM_REAL),  dimension(:),         allocatable                :: valence
    real(kind=CUSTOM_REAL)                                                    :: penalty
    integer                                                                   :: i,j,k,ispec,iglob

  integer :: ier

    allocate(field(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 193')
    allocate(numerical_laplacian_of_field(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 194')
    allocate(valence(NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 195')
    allocate(field_to_derivate(NGLLX,NGLLX,NGLLX,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 196')

    field(:,:)=0._CUSTOM_REAL
    valence(:)=0._CUSTOM_REAL
    penalty=0._CUSTOM_REAL

    !! mean values of field on element boundary
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
!!$                field(1,iglob)=field(1,iglob) + current_model_vp(i,j,k,ispec)
!!$                field(2,iglob)=field(2,iglob) + current_model_vs(i,j,k,ispec)
!!$                field(3,iglob)=field(3,iglob) + current_model_rh(i,j,k,ispec)
                valence(iglob)=valence(iglob) + 1._CUSTOM_REAL
             enddo
          enddo
       enddo
    enddo

    field(1,:)=field(1,:)/valence(:)
    field(2,:)=field(2,:)/valence(:)
    field(3,:)=field(3,:)/valence(:)

    !! VP
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                field_to_derivate(i,j,k,ispec)=field(1,iglob)
             enddo
          enddo
       enddo
    enddo
    call  compute_derivatives_with_interpolation(field_to_derivate, lap1_vp, Lap2_vp)

   !! VS
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                field_to_derivate(i,j,k,ispec)=field(2,iglob)
             enddo
          enddo
       enddo
    enddo
    call  compute_derivatives_with_interpolation(field_to_derivate, lap1_vs, Lap2_vs)

    !! rho
    do ispec=1,NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                field_to_derivate(i,j,k,ispec)=field(3,iglob)
             enddo
          enddo
       enddo
    enddo
    call  compute_derivatives_with_interpolation(field_to_derivate, lap1_rh, Lap2_rh)






!!$    !! VP  ------------------------------------------------------------------------
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap1_vp(i,j,k,ispec)= numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$
!!$    field(1,:)=numerical_laplacian_of_field(:)
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap2_vp(i,j,k,ispec)=numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    !! VS -----------------------------------------------------------------------------
!!$    field(1,:)=field(2,:)
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap1_vs(i,j,k,ispec)=numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    field(1,:)=numerical_laplacian_of_field(:)
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap2_vs(i,j,k,ispec)=numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$
!!$    !! RHO --------------------------------------------------------------------------
!!$    field(1,:)=field(3,:)
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap1_rh(i,j,k,ispec)=numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    field(1,:)=numerical_laplacian_of_field(:)
!!$    call compute_laplacian(field, numerical_laplacian_of_field, regularization_fd)
!!$
!!$    do ispec=1,NSPEC_AB
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                Lap2_rh(i,j,k,ispec)=numerical_laplacian_of_field(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo


    deallocate(field, numerical_laplacian_of_field, valence, field_to_derivate)

  end subroutine compute_laplacian2_vp_vs_rho

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute mean value at edge of element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_mean_values_on_edge(field)

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),    allocatable, intent(inout)      :: field
    real(kind=CUSTOM_REAL), dimension(:,:),          allocatable                     :: field_wksp, valence
    integer                                                       :: ispec, iglob, i, j, k

    integer :: ier

    allocate(field_wksp(NDIM,NGLOB_AB), valence(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 197')

    field_wksp(:,:) = 0.
    valence(:,:)=0.

    !! 1/ compute average value on boundary of elements
    do ispec=1, NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)

                field_wksp(1,iglob) = field_wksp(1,iglob) + field(1,i,j,k,ispec)
                field_wksp(2,iglob) = field_wksp(2,iglob) + field(2,i,j,k,ispec)
                field_wksp(3,iglob) = field_wksp(3,iglob) + field(3,i,j,k,ispec)

                valence(1,iglob) = valence(1,iglob) + 1.
                valence(2,iglob) = valence(2,iglob) + 1.
                valence(3,iglob) = valence(3,iglob) + 1.

             enddo
          enddo
       enddo
    enddo

    field_wksp(:,:) = field_wksp(:,:) /  valence(:,:)

    valence(:,:) = 1.

    !! 2/ compute average values of field_to_derivate at the edge of MPI slices -------------------
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,valence, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,field_wksp, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    field_wksp(:,:) = field_wksp(:,:) /  valence(:,:)

    do ispec = 1, NSPEC_AB
       do k = 1,NGLLZ
          do j = 1,NGLLY
             do i = 1,NGLLX
                iglob = ibool(i,j,k,ispec)

                field(1,i,j,k,ispec) = field_wksp(1,iglob)
                field(2,i,j,k,ispec) = field_wksp(2,iglob)
                field(3,i,j,k,ispec) = field_wksp(3,iglob)

             enddo
          enddo
       enddo
    enddo

    deallocate(field_wksp, valence)

  end subroutine compute_mean_values_on_edge


!###################################################################################################################################
!######################## SUBORUITNES RELATED TO GLL SHAPE and INTERPOLATION #######################################################
!######################## SUBORUITNES RELATED TO GLL SHAPE and INTERPOLATION #######################################################
!######################## SUBORUITNES RELATED TO GLL SHAPE and INTERPOLATION #######################################################
!######################## SUBORUITNES RELATED TO GLL SHAPE and INTERPOLATION #######################################################
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute lagange and GLL for higher interpolation
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine setup_interpolation_for_higher_degree()

    integer                              :: i1, i2, ier
    double precision                     :: lagrange_deriv_GLL
    double precision,dimension(NGLLX)    :: h,hprime
    double precision,dimension(NGLLd)    :: h1,hprime1

    call zwgljd(gll_points, wgll_points, NGLLd,GAUSSALPHA,GAUSSBETA)

    do i1=1,NGLLd
       call lagrange_any(gll_points(i1),NGLLX,xigll,h,hprime)
       hlagrange(:,i1)=h(:)
    enddo

    do i1=1,NGLLd
       do i2=1,NGLLd
          hlagrange_prime(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,gll_points,NGLLd)
       enddo
    enddo

    do i1=1,NGLLX
       call lagrange_any(xigll(i1), NGLLd, gll_points, h1, hprime1)
       hlagrange_old(:,i1)=h1(:)
    enddo

    ! allocates dershape_function array (uses dynamic allocation to avoid stack issues with newer gcc compilers)
    allocate(dershape_function(3,8,NGLLd,NGLLd,NGLLd),stat=ier)
    if (ier /= 0) stop 'Error allocating dershape_function array'

    call get_shape3D_genric(NGLLd, gll_points, shape_function, dershape_function)

  end subroutine setup_interpolation_for_higher_degree

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute lagrange interpolation in new GLL points
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine interpolation_in_new_gll(field_in_element, field_in_higher_degree)

    double precision, dimension(NGLLX,NGLLX,NGLLX), intent(in)    :: field_in_element
    double precision, dimension(NGLLd,NGLLd,NGLLd), intent(inout) :: field_in_higher_degree

    integer :: i, j, k
    integer :: ih, jh, kh

    field_in_higher_degree(:,:,:) = 0.
    do kh=1,NGLLd
       do jh=1,NGLLd
          do ih=1,NGLLd

             do k=1,NGLLX
                do j=1,NGLLX
                   do i=1,NGLLX
                      field_in_higher_degree(ih,jh,kh) = field_in_higher_degree(ih,jh,kh) + &
                           hlagrange(i,ih)*hlagrange(j,jh)*hlagrange(k,kh)*field_in_element(i,j,k)
                   enddo
                enddo
             enddo

          enddo
       enddo
    enddo

  end subroutine interpolation_in_new_gll

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute lagrange interpolation in new GLL points
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine interpolation_in_old_gll(field_in_higher_degree, field_in_element)

    double precision, dimension(NGLLX,NGLLX,NGLLX), intent(inout) :: field_in_element
    double precision, dimension(NGLLd,NGLLd,NGLLd), intent(in)    :: field_in_higher_degree

    integer :: i, j, k
    integer :: ih, jh, kh

    field_in_element(:,:,:) = 0.
    do k=1,NGLLX
       do j=1,NGLLX
          do i=1,NGLLX
             do kh=1,NGLLd
                do jh=1,NGLLd
                   do ih=1,NGLLd
                      field_in_element(i,j,k) = field_in_element(i,j,k) + &
                           hlagrange_old(ih,i)*hlagrange_old(jh,j)*hlagrange_old(kh,k)*field_in_higher_degree(ih,jh,kh)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine interpolation_in_old_gll

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute shape function for a given 8 node element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine get_shape3D_genric(NGLL, gll_coord, shape_func, dershape_func)

    integer,                                          intent(in)    :: NGLL
    double precision, dimension(NGLL),                intent(in)    :: gll_coord
    double precision, dimension(8,NGLL,NGLL,NGLL),    intent(inout) :: shape_func
    double precision, dimension(3, 8,NGLL,NGLL,NGLL), intent(inout) :: dershape_func

    double precision                                                :: xi, eta, gamma
    double precision                                                :: ra1, ra2, rb1, rb2, rc1, rc2
    integer                                                         :: i,j,k

    do i=1,NGLL
       do j=1,NGLL
          do k=1,NGLL

             xi =    gll_coord(i)
             eta =   gll_coord(j)
             gamma = gll_coord(k)

             ra1 = one + xi
             ra2 = one - xi

             rb1 = one + eta
             rb2 = one - eta

             rc1 = one + gamma
             rc2 = one - gamma

             shape_func(1,i,j,k) = ONE_EIGHTH*ra2*rb2*rc2
             shape_func(2,i,j,k) = ONE_EIGHTH*ra1*rb2*rc2
             shape_func(3,i,j,k) = ONE_EIGHTH*ra1*rb1*rc2
             shape_func(4,i,j,k) = ONE_EIGHTH*ra2*rb1*rc2
             shape_func(5,i,j,k) = ONE_EIGHTH*ra2*rb2*rc1
             shape_func(6,i,j,k) = ONE_EIGHTH*ra1*rb2*rc1
             shape_func(7,i,j,k) = ONE_EIGHTH*ra1*rb1*rc1
             shape_func(8,i,j,k) = ONE_EIGHTH*ra2*rb1*rc1

             dershape_func(1,1,i,j,k) = - ONE_EIGHTH*rb2*rc2
             dershape_func(1,2,i,j,k) =   ONE_EIGHTH*rb2*rc2
             dershape_func(1,3,i,j,k) =   ONE_EIGHTH*rb1*rc2
             dershape_func(1,4,i,j,k) = - ONE_EIGHTH*rb1*rc2
             dershape_func(1,5,i,j,k) = - ONE_EIGHTH*rb2*rc1
             dershape_func(1,6,i,j,k) =   ONE_EIGHTH*rb2*rc1
             dershape_func(1,7,i,j,k) =   ONE_EIGHTH*rb1*rc1
             dershape_func(1,8,i,j,k) = - ONE_EIGHTH*rb1*rc1

             dershape_func(2,1,i,j,k) = - ONE_EIGHTH*ra2*rc2
             dershape_func(2,2,i,j,k) = - ONE_EIGHTH*ra1*rc2
             dershape_func(2,3,i,j,k) =   ONE_EIGHTH*ra1*rc2
             dershape_func(2,4,i,j,k) =   ONE_EIGHTH*ra2*rc2
             dershape_func(2,5,i,j,k) = - ONE_EIGHTH*ra2*rc1
             dershape_func(2,6,i,j,k) = - ONE_EIGHTH*ra1*rc1
             dershape_func(2,7,i,j,k) =   ONE_EIGHTH*ra1*rc1
             dershape_func(2,8,i,j,k) =   ONE_EIGHTH*ra2*rc1

             dershape_func(3,1,i,j,k) = - ONE_EIGHTH*ra2*rb2
             dershape_func(3,2,i,j,k) = - ONE_EIGHTH*ra1*rb2
             dershape_func(3,3,i,j,k) = - ONE_EIGHTH*ra1*rb1
             dershape_func(3,4,i,j,k) = - ONE_EIGHTH*ra2*rb1
             dershape_func(3,5,i,j,k) =   ONE_EIGHTH*ra2*rb2
             dershape_func(3,6,i,j,k) =   ONE_EIGHTH*ra1*rb2
             dershape_func(3,7,i,j,k) =   ONE_EIGHTH*ra1*rb1
             dershape_func(3,8,i,j,k) =   ONE_EIGHTH*ra2*rb1

          enddo
       enddo
    enddo

  end subroutine get_shape3D_genric

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute shape function for a given 8 node element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_jacobian_one_element(xnode, ynode, znode, Jacobian_shape)

    double precision,  dimension(NGNOD),                        intent(in)     :: xnode, ynode, znode
    double precision,  dimension(NGLLd,NGLLd,NGLLd,NDIM, NDIM), intent(inout)  :: Jacobian_shape

    integer          ::  ia, i,j,k
    double precision ::  xxi, xeta, xgamma
    double precision ::  yxi, yeta, ygamma
    double precision ::  zxi, zeta, zgamma
    double precision ::  xix, etax, gammax
    double precision ::  xiy, etay, gammay
    double precision ::  xiz, etaz, gammaz
    double precision ::  xmesh, ymesh, zmesh
    double precision :: jacobian

    do k=1,NGLLd
       do j=1,NGLLd
          do i=1,NGLLd

             xxi = ZERO
             xeta = ZERO
             xgamma = ZERO
             yxi = ZERO
             yeta = ZERO
             ygamma = ZERO
             zxi = ZERO
             zeta = ZERO
             zgamma = ZERO
             xmesh = ZERO
             ymesh = ZERO
             zmesh = ZERO

             do ia=1,NGNOD
                xxi    = xxi    + dershape_function(1,ia,i,j,k)*xnode(ia)
                xeta   = xeta   + dershape_function(2,ia,i,j,k)*xnode(ia)
                xgamma = xgamma + dershape_function(3,ia,i,j,k)*xnode(ia)

                yxi    = yxi    + dershape_function(1,ia,i,j,k)*ynode(ia)
                yeta   = yeta   + dershape_function(2,ia,i,j,k)*ynode(ia)
                ygamma = ygamma + dershape_function(3,ia,i,j,k)*ynode(ia)

                zxi    = zxi    + dershape_function(1,ia,i,j,k)*znode(ia)
                zeta   = zeta   + dershape_function(2,ia,i,j,k)*znode(ia)
                zgamma = zgamma + dershape_function(3,ia,i,j,k)*znode(ia)

                xmesh = xmesh + shape_function(ia,i,j,k)*xnode(ia)
                ymesh = ymesh + shape_function(ia,i,j,k)*ynode(ia)
                zmesh = zmesh + shape_function(ia,i,j,k)*znode(ia)
             enddo

             jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
                       xeta*(yxi*zgamma-ygamma*zxi) + &
                       xgamma*(yxi*zeta-yeta*zxi)

             ! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
             if (jacobian <= ZERO) call exit_MPI(myrank,'Error negative or null 3D Jacobian found')

             !     invert the relation (Fletcher p. 50 vol. 2)
             xix = (yeta*zgamma-ygamma*zeta) / jacobian
             xiy = (xgamma*zeta-xeta*zgamma) / jacobian
             xiz = (xeta*ygamma-xgamma*yeta) / jacobian
             etax = (ygamma*zxi-yxi*zgamma) / jacobian
             etay = (xxi*zgamma-xgamma*zxi) / jacobian
             etaz = (xgamma*yxi-xxi*ygamma) / jacobian
             gammax = (yxi*zeta-yeta*zxi) / jacobian
             gammay = (xeta*zxi-xxi*zeta) / jacobian
             gammaz = (xxi*yeta-xeta*yxi) / jacobian

             !     compute and store the jacobian for the solver
             jacobian = 1. / (xix*(etay*gammaz-etaz*gammay) &
                  -xiy*(etax*gammaz-etaz*gammax) &
                  +xiz*(etax*gammay-etay*gammax))

             !     save the derivatives and the jacobian

             ! distinguish between single and double precision for reals
             Jacobian_shape(i,j,k,1,1) = xix
             Jacobian_shape(i,j,k,1,2) = xiy
             Jacobian_shape(i,j,k,1,3) = xiz
             Jacobian_shape(i,j,k,2,1) = etax
             Jacobian_shape(i,j,k,2,2) = etay
             Jacobian_shape(i,j,k,2,3) = etaz
             Jacobian_shape(i,j,k,3,1) = gammax
             Jacobian_shape(i,j,k,3,2) = gammay
             Jacobian_shape(i,j,k,3,3) = gammaz

!!$             jacobianstore(i,j,k) = jacobian
             xstore_interp(i,j,k) = xmesh
             ystore_interp(i,j,k) = ymesh
             zstore_interp(i,j,k) = zmesh
!!$
          enddo
       enddo
    enddo

  end subroutine compute_jacobian_one_element

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute 2 and 4 derivatives in one element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_derivatives_with_interpolation(field_to_derivate, laplacian_of_field, double_laplacian_of_field)

    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(in)    ::  field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(inout) ::  laplacian_of_field, double_laplacian_of_field

    integer :: ispec, iglob, ier

    ! allocates jacobian array for element (uses dynamic allocation to avoid stack issues with gcc compilers version >= 10)
    if (.not. allocated(Jacobian_shape_function)) then
      allocate(Jacobian_shape_function(NGLLd,NGLLd,NGLLd,3,3),stat=ier)
      if (ier /= 0) stop 'Error allocating Jacobian_shape_function array'
    endif

    do ispec = 1, NSPEC_AB

       !! store  locals
       field_initial(:,:,:) = field_to_derivate(:,:,:,ispec)

       !! store elements control points
       iglob=ibool(1,1,1,ispec)
       xnodelm(1)=xstore(iglob)
       ynodelm(1)=ystore(iglob)
       znodelm(1)=zstore(iglob)

       iglob=ibool(NGLLX,1,1,ispec)
       xnodelm(2)=xstore(iglob)
       ynodelm(2)=ystore(iglob)
       znodelm(2)=zstore(iglob)

       iglob=ibool(NGLLX,NGLLY,1,ispec)
       xnodelm(3)=xstore(iglob)
       ynodelm(3)=ystore(iglob)
       znodelm(3)=zstore(iglob)

       iglob=ibool(1,NGLLY,1,ispec)
       xnodelm(4)=xstore(iglob)
       ynodelm(4)=ystore(iglob)
       znodelm(4)=zstore(iglob)

       iglob=ibool(1,1,NGLLZ,ispec)
       xnodelm(5)=xstore(iglob)
       ynodelm(5)=ystore(iglob)
       znodelm(5)=zstore(iglob)

       iglob=ibool(NGLLX,1,NGLLZ,ispec)
       xnodelm(6)=xstore(iglob)
       ynodelm(6)=ystore(iglob)
       znodelm(6)=zstore(iglob)

       iglob=ibool(NGLLX,NGLLY,NGLLZ,ispec)
       xnodelm(7)=xstore(iglob)
       ynodelm(7)=ystore(iglob)
       znodelm(7)=zstore(iglob)

       iglob=ibool(1,NGLLY,NGLLZ,ispec)
       xnodelm(8)=xstore(iglob)
       ynodelm(8)=ystore(iglob)
       znodelm(8)=zstore(iglob)

       call compute_jacobian_one_element(xnodelm, ynodelm, znodelm, Jacobian_shape_function)
       call interpolation_in_new_gll(field_initial, field_interpolated)
       !call compute_function_to_test_in_one_element(field_interpolated)

       call compute_laplacian_lagrange_element(Laplacian_of_field_interpolated, field_interpolated)
       call interpolation_in_old_gll(Laplacian_of_field_interpolated, field_initial)
       laplacian_of_field(:,:,:,ispec) = field_initial(:,:,:)

       field_interpolated(:,:,:) = Laplacian_of_field_interpolated(:,:,:)
       call compute_laplacian_lagrange_element(Laplacian_of_field_interpolated, field_interpolated)
       call interpolation_in_old_gll(Laplacian_of_field_interpolated, field_initial)
       double_laplacian_of_field(:,:,:,ispec) = field_initial(:,:,:)

    enddo

  end subroutine compute_derivatives_with_interpolation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute laplacian in 1 element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_laplacian_lagrange_element(Laplacian_of_field_in_element, field_in_element)

    double precision, dimension(NGLLd,NGLLd,NGLLd), intent(in)    :: field_in_element
    double precision, dimension(NGLLd,NGLLd,NGLLd), intent(inout) :: Laplacian_of_field_in_element

    integer          :: i,j,k,l
    double precision :: hp1, hp2, hp3
    double precision :: tempx1l,  tempx2l,  tempx3l
    double precision :: tempy1l,  tempy2l,  tempy3l
    double precision :: tempz1l,  tempz2l,  tempz3l
    double precision :: xixl, xiyl, xizl
    double precision :: etaxl, etayl, etazl
    double precision :: gammaxl, gammayl, gammazl

    do k=1,NGLLd
       do j=1, NGLLd
          do i=1,NGLLd

             tempx1l = 0.d0
             tempx2l = 0.d0
             tempx3l = 0.d0
             do l=1,NGLLd

                hp1 = hlagrange_prime(i,l)
                tempx1l = tempx1l + field_in_element(l,j,k)*hp1

                hp2 = hlagrange_prime(j,l)
                tempx2l = tempx2l + field_in_element(i,l,k)*hp2

                hp3 = hlagrange_prime(k,l)
                tempx3l = tempx3l + field_in_element(i,j,l)*hp3

             enddo

             xixl    = Jacobian_shape_function(i,j,k,1,1)
             xiyl    = Jacobian_shape_function(i,j,k,1,2)
             xizl    = Jacobian_shape_function(i,j,k,1,3)
             etaxl   = Jacobian_shape_function(i,j,k,2,1)
             etayl   = Jacobian_shape_function(i,j,k,2,2)
             etazl   = Jacobian_shape_function(i,j,k,2,3)
             gammaxl = Jacobian_shape_function(i,j,k,3,1)
             gammayl = Jacobian_shape_function(i,j,k,3,2)
             gammazl = Jacobian_shape_function(i,j,k,3,3)

             dfdx(i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
             dfdy(i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
             dfdz(i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l


          enddo
       enddo
    enddo

    do k=1,NGLLd
       do j=1, NGLLd
          do i=1,NGLLd

             tempx1l = 0.d0
             tempx2l = 0.d0
             tempx3l = 0.d0
             tempy1l = 0.d0
             tempy2l = 0.d0
             tempy3l = 0.d0
             tempz1l = 0.d0
             tempz2l = 0.d0
             tempz3l = 0.d0
             do l=1,NGLLd

                hp1 = hlagrange_prime(i,l)
                tempx1l = tempx1l + dfdx(l,j,k)*hp1
                tempy1l = tempy1l + dfdy(l,j,k)*hp1
                tempz1l = tempz1l + dfdz(l,j,k)*hp1

                hp2 = hlagrange_prime(j,l)
                tempx2l = tempx2l + dfdx(i,l,k)*hp2
                tempy2l = tempy2l + dfdy(i,l,k)*hp2
                tempz2l = tempz2l + dfdz(i,l,k)*hp2

                hp3 = hlagrange_prime(k,l)
                tempx3l = tempx3l + dfdx(i,j,l)*hp3
                tempy3l = tempy3l + dfdy(i,j,l)*hp3
                tempz3l = tempz3l + dfdz(i,j,l)*hp3

             enddo

             xixl    = Jacobian_shape_function(i,j,k,1,1)
             xiyl    = Jacobian_shape_function(i,j,k,1,2)
             xizl    = Jacobian_shape_function(i,j,k,1,3)
             etaxl   = Jacobian_shape_function(i,j,k,2,1)
             etayl   = Jacobian_shape_function(i,j,k,2,2)
             etazl   = Jacobian_shape_function(i,j,k,2,3)
             gammaxl = Jacobian_shape_function(i,j,k,3,1)
             gammayl = Jacobian_shape_function(i,j,k,3,2)
             gammazl = Jacobian_shape_function(i,j,k,3,3)

             Laplacian_of_field_in_element(i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l + &
                                                    xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l + &
                                                    xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l
          enddo
       enddo
    enddo

  end subroutine compute_laplacian_lagrange_element

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute laplacian in 1 element
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_function_to_test_in_one_element(field_to_define)

    double precision, dimension(NGLLd,NGLLd,NGLLd), intent(inout)    :: field_to_define
    double precision :: x,y,z
    integer i,j,k

    do k=1,NGLLd
       do j=1, NGLLd
          do i=1,NGLLd
             x = xstore_interp(i,j,k)
             y = ystore_interp(i,j,k)
             z = zstore_interp(i,j,k)
             field_to_define(i,j,k)= cos( lambda * x ) * cos( lambda* y ) * cos( lambda * z ) &
                  / (3 * lambda * lambda * lambda * lambda)
          enddo
       enddo
    enddo


  end subroutine compute_function_to_test_in_one_element

!########################################################################################################################
!################################################# UNIT TESTING #########################################################
!########################################################################################################################

!---------------------------------------------------------------------------------
!  just for testing and debugging
!---------------------------------------------------------------------------------

! subroutine testing_regularization(regularization_fd)

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

!   implicit none

!   type(regul), dimension(:), allocatable, intent(inout)   :: regularization_fd

!   ! locals
!   real(kind=CUSTOM_REAL), dimension(:),       allocatable :: field_to_derivate
!   real(kind=CUSTOM_REAL), dimension(:,:),     allocatable :: field_glob, Laplac_boundary
!   real(kind=CUSTOM_REAL), dimension(:,:),     allocatable :: check_glob

!   real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: field_gll
!   real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: check_gll

!   integer                                                 :: ig, i, j, k, ie
!   integer                                                 :: iglob, iglob_index
!   real(kind=CUSTOM_REAL)                                  :: omega, lambda
!   character(len=256)                                      :: name_file, path_file

!   !! compute Finite Difference Derivative Matrix and store in structure --------------------------------
!   !!call compute_and_store_FD_derivatives_matrix(regularization_fd)

!   !! allocate arrays for testing ---------------------------
!   allocate(field_gll(NGLLX, NGLLY, NGLLZ, NSPEC_AB))
!   allocate(check_gll(NGLLX, NGLLY, NGLLZ, NSPEC_AB))

!   allocate(field_glob(NDIM,NGLOB_AB))
!   allocate(check_glob(NDIM,NGLOB_AB))

!! DK DK not clear where this value of 50000 comes from
!   lambda=50000
!   omega = 2*3.1459265359/lambda

!   !! 1/ define fields and valence to process -----------------
!   do ig=1,NGLOB_AB
!      field_glob(:,ig)=cos(omega*xstore(ig))*cos(omega*ystore(ig))*cos(omega*zstore(ig))
!   enddo

!   !! in: field_glob(NDIM,NGLOB_AB)    out: check_glob(NDIM,NGLOB_AB)
!   allocate(Laplac_boundary(NDIM, Nb_iglob_on_faces), field_to_derivate(NGLOB_AB))

!   !! 1/ compute FD derivatives
!   call compute_laplacian_FD(Laplac_boundary, field_glob, regularization_fd)
!   !! 2/ compute deivatives inside elements by lagrange interpolation
!   field_to_derivate(:)=field_glob(1,:)
!   call compute_laplac_lagrange(check_glob, field_to_derivate)
!   !! 3/ replace boundary elements  !! for debug
!   do iglob_index=1,  Nb_iglob_on_faces
!      iglob = regularization_fd(iglob_index)%iglob
!      check_glob(1, iglob) = Laplac_boundary(1, iglob_index)
!   enddo

!   deallocate(Laplac_boundary, field_to_derivate)

!   ! call compute_first_derivatives_lagrange(check_glob, field_to_derivate)  !! to test

!   if (DEBUG_MODE)  write(IIDD,*)

!   !! save in disk ----------------------------------
!   ! store it for vtk visualization
!   do ie = 1, NSPEC_AB
!      do k=1,NGLLZ
!         do j=1,NGLLY
!            do i=1,NGLLX
!               ig=ibool(i,j,k,ie)
!               check_gll(i,j,k,ie)=check_glob(1,ig)
!            enddo
!         enddo
!      enddo
!   enddo

!   !! store in all GLL just for vtk visualization --------------------------------------
!   do ie=1,NSPEC_AB
!      do k=1,NGLLZ
!         do j=1,NGLLY
!            do i=1,NGLLX
!               ig=ibool(i,j,k,ie)
!               field_gll(i,j,k,ie)   = field_glob(1,ig)
!!$                valence_gll(i,j,k,ie) = valence_glob(1,ig)
!            enddo
!         enddo
!      enddo
!   enddo
!   path_file='OUTPUT_FILES/DATABASES_MPI/proc'
!   write(name_file,'(i6.6,a21)') myrank, '_model_test_field.bin'
!   path_file=(trim(path_file))//trim(name_file)
!   open(888,file=trim(path_file),form='unformatted')
!   write(888) field_gll
!   close(888)

!!$    path_file='OUTPUT_FILES/DATABASES_MPI/proc'
!!$    write(name_file,'(i6.6,a16)') myrank, '_valence_gll.bin'
!!$    path_file=(trim(path_file))//trim(name_file)
!!$    open(888,file=trim(path_file),form='unformatted')
!!$    write(888) valence_gll
!!$    close(888)

!   path_file='OUTPUT_FILES/DATABASES_MPI/proc'
!   write(name_file,'(i6.6,a14)') myrank, '_check_gll.bin'
!   path_file=(trim(path_file))//trim(name_file)
!   open(888,file=trim(path_file),form='unformatted')
!   write(888) check_gll
!   close(888)

!   deallocate(field_gll, check_gll, field_glob, check_glob, Laplac_boundary)

! end subroutine testing_regularization

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! send / recv of buffer for overlap of field in MPI slices with neighbor MPI slices
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine send_recv_blocking(buffer_to_send, buffer_to_recv)

    implicit none

    real(kind=CUSTOM_REAL), dimension(:), allocatable ::  buffer_to_send, buffer_to_recv

    integer :: irank, igll, ishift
    integer :: itag=0

    do irank = 0, NPROC-1
       if (struct_comm(irank)%ns > 0) then
          !! store array to send
           do igll = 1, struct_comm(irank)%ns
              struct_comm(irank)%array_to_send(igll)=buffer_to_send(struct_comm(irank)%iglob_to_send(igll))
           enddo
            call isend_cr(struct_comm(irank)%array_to_send,struct_comm(irank)%ns, irank, itag,  request_send_scalar(irank))
        endif
        if (struct_comm(irank)%nr > 0) then  !!
           call irecv_cr(struct_comm(irank)%array_to_recv, struct_comm(irank)%nr, irank, itag, request_recv_scalar(irank))
        endif
    enddo

    do irank = 0, NPROC - 1
       if (struct_comm(irank)%nr > 0) call wait_req(request_recv_scalar(irank))
    enddo

    do irank = 0, NPROC - 1
       if (struct_comm(irank)%ns > 0) call wait_req(request_send_scalar(irank))
    enddo

    do irank = 0, NPROC-1
       ishift=struct_comm(irank)%ibegin
       if (struct_comm(irank)%nr > 0) then
          do igll = 1, struct_comm(irank)%nr
             buffer_to_recv(igll + ishift) = struct_comm(irank)%array_to_recv(igll)
          enddo
       endif
    enddo

  end subroutine send_recv_blocking

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! for each GLL in element compute first derivatives using lagrange interpolation and derivatives
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_first_derivatives_lagrange(Df, field_to_derivate)

    use specfem_par, only: xixstore, xiystore, xizstore, &
                           etaxstore, etaystore, etazstore, &
                           gammaxstore, gammaystore, gammazstore, &
                           hprime_xx, irregular_element_number, &
                           xix_regular

    implicit none
    real(kind=CUSTOM_REAL), dimension(:),         allocatable, intent(in)     :: field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable, intent(inout)  :: Df

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ)              :: dummyloc
    real(kind=CUSTOM_REAL)                                              :: hp1, hp2, hp3
    real(kind=CUSTOM_REAL)                                              :: tempx1l, tempx2l, tempx3l
    real(kind=CUSTOM_REAL)                                              :: xixl,xiyl,xizl
    real(kind=CUSTOM_REAL)                                              :: etaxl,etayl,etazl
    real(kind=CUSTOM_REAL)                                              :: gammaxl,gammayl,gammazl
    integer                                                             :: ispec, ispec_irreg, iglob, i, j, k, l


    do ispec =1, NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                dummyloc(i,j,k)=field_to_derivate(iglob)
             enddo
          enddo
       enddo

       ispec_irreg = irregular_element_number(ispec)

       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX

                iglob = ibool(i,j,k,ispec)

                tempx1l = 0.
                tempx2l = 0.
                tempx3l = 0.

                do l=1,NGLLX

                   hp1 = hprime_xx(i,l)
                   tempx1l = tempx1l + dummyloc(l,j,k)*hp1

                   hp2 = hprime_xx(j,l)
                   tempx2l = tempx2l + dummyloc(i,l,k)*hp2

                   hp3 = hprime_xx(k,l)
                   tempx3l = tempx3l + dummyloc(i,j,l)*hp3

                enddo
                if (ispec_irreg /= 0) then
                  !irregular element
                  xixl = xixstore(i,j,k,ispec_irreg)
                  xiyl = xiystore(i,j,k,ispec_irreg)
                  xizl = xizstore(i,j,k,ispec_irreg)
                  etaxl = etaxstore(i,j,k,ispec_irreg)
                  etayl = etaystore(i,j,k,ispec_irreg)
                  etazl = etazstore(i,j,k,ispec_irreg)
                  gammaxl = gammaxstore(i,j,k,ispec_irreg)
                  gammayl = gammaystore(i,j,k,ispec_irreg)
                  gammazl = gammazstore(i,j,k,ispec_irreg)

                  Df(1,i,j,k,ispec) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                  Df(2,i,j,k,ispec) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                  Df(3,i,j,k,ispec) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

                else
                  !regular element
                  Df(1,i,j,k,ispec) = xix_regular * tempx1l
                  Df(2,i,j,k,ispec) = xix_regular * tempx2l
                  Df(3,i,j,k,ispec) = xix_regular * tempx3l

                endif

             enddo
          enddo
       enddo

    enddo

  end subroutine compute_first_derivatives_lagrange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! for each GLL in element compute laplacian using lagrange interpolation and derivatives
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_laplac_lagrange(Lapf, field_to_derivate)

    use specfem_par, only: xixstore, xiystore, xizstore, etaxstore, etaystore, etazstore, &
                           gammaxstore, gammaystore, gammazstore, hprime_xx, irregular_element_number, &
                           xix_regular


    implicit none
    real(kind=CUSTOM_REAL), dimension(:),     allocatable, intent(in)     :: field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:,:),   allocatable, intent(inout)  :: Lapf

    double precision, dimension(NGLLX, NGLLY, NGLLZ)              :: F
    double precision, dimension(NDIM, NGLLX, NGLLY, NGLLZ)        :: DF
    double precision                                              :: hp1, hp2, hp3
    double precision                                              :: tempx1l, tempx2l, tempx3l
    double precision                                              :: xixl,xiyl,xizl
    double precision                                              :: etaxl,etayl,etazl
    double precision                                              :: gammaxl,gammayl,gammazl
    double precision                                              :: coef_norm
    integer                                                       :: ispec, ispec_irreg, iglob, i, j, k, l
    integer :: ispec_to_debug = 100

    if (DEBUG_MODE) write(IIDD,*)

    do ispec =1, NSPEC_AB

       !! first derivatives
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)
                F(i,j,k)=field_to_derivate(iglob)
             enddo
          enddo
       enddo

       coef_norm = maxval(abs(F(:,:,:)))
       !coef_norm=1.
       F(:,:,:)=F(:,:,:)/coef_norm
       ispec_irreg = irregular_element_number(ispec)

       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX

                iglob = ibool(i,j,k,ispec)

                tempx1l = 0.
                tempx2l = 0.
                tempx3l = 0.

                do l=1,NGLLX

                   hp1 = hprime_xx(i,l)
                   tempx1l = tempx1l + F(l,j,k)*hp1

                   hp2 = hprime_xx(j,l)
                   tempx2l = tempx2l + F(i,l,k)*hp2

                   hp3 = hprime_xx(k,l)
                   tempx3l = tempx3l + F(i,j,l)*hp3

                enddo

                if (ispec_irreg /= 0) then
                  !irregular element
                  xixl = xixstore(i,j,k,ispec_irreg)
                  xiyl = xiystore(i,j,k,ispec_irreg)
                  xizl = xizstore(i,j,k,ispec_irreg)
                  etaxl = etaxstore(i,j,k,ispec_irreg)
                  etayl = etaystore(i,j,k,ispec_irreg)
                  etazl = etazstore(i,j,k,ispec_irreg)
                  gammaxl = gammaxstore(i,j,k,ispec_irreg)
                  gammayl = gammaystore(i,j,k,ispec_irreg)
                  gammazl = gammazstore(i,j,k,ispec_irreg)

                  dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                  dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                  dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

                else
                  !regular element
                  dF(1,i,j,k) = xix_regular * tempx1l
                  dF(2,i,j,k) = xix_regular * tempx2l
                  dF(3,i,j,k) = xix_regular * tempx3l

                endif

                if (DEBUG_MODE .and. ispec == ispec_to_debug .and. i == 1 .and. j == 1) write(IIDD,*)   &
                     k,zstore(iglob), dF(3,i,j,k)*coef_norm,  F(i,j,k)*coef_norm
                !if (DEBUG_MODE) write(IIDD,*)   dF(1,i,j,k),  dF(2,i,j,k), dF(3,i,j,k),  F(i,j,k)
                !LapF(1,iglob) = dF(3,i,j,k)
             enddo
          enddo
       enddo

       dF(:,:,:,:) = dF(:,:,:,:) * coef_norm

       !! deivate again df/dx / dx
       F(:,:,:)=dF(1,:,:,:)

       coef_norm = maxval(abs(F(:,:,:)))
       !coef_norm=1.
       F(:,:,:)=F(:,:,:)/coef_norm

       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX

                iglob = ibool(i,j,k,ispec)

                tempx1l = 0.
                tempx2l = 0.
                tempx3l = 0.

                do l=1,NGLLX

                   hp1 = hprime_xx(i,l)
                   tempx1l = tempx1l + F(l,j,k)*hp1

                   hp2 = hprime_xx(j,l)
                   tempx2l = tempx2l + F(i,l,k)*hp2

                   hp3 = hprime_xx(k,l)
                   tempx3l = tempx3l + F(i,j,l)*hp3

                enddo
                if (ispec_irreg /= 0) then
                  !irregular element
                  xixl = xixstore(i,j,k,ispec_irreg)
                  xiyl = xiystore(i,j,k,ispec_irreg)
                  xizl = xizstore(i,j,k,ispec_irreg)
                  etaxl = etaxstore(i,j,k,ispec_irreg)
                  etayl = etaystore(i,j,k,ispec_irreg)
                  etazl = etazstore(i,j,k,ispec_irreg)
                  gammaxl = gammaxstore(i,j,k,ispec_irreg)
                  gammayl = gammaystore(i,j,k,ispec_irreg)
                  gammazl = gammazstore(i,j,k,ispec_irreg)

                  dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                  !dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                  !dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

                else
                  !regular element
                  dF(1,i,j,k) = xix_regular * tempx1l
                  !dF(2,i,j,k) = xix_regular * tempx2l
                  !dF(3,i,j,k) = xix_regular * tempx3l

                endif

                dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                !dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                !dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l
                !if (DEBUG_MODE) write(IIDD,*) 'x', dF(1,i,j,k)
                LapF(1,iglob) = dF(1,i,j,k) * coef_norm

             enddo
          enddo
       enddo

       if (DEBUG_MODE .and. ispec == ispec_to_debug) write(IIDD,*)
       !! deivate again df/dy / dy
       F(:,:,:)=dF(2,:,:,:)

       coef_norm = maxval(abs(F(:,:,:)))
       !coef_norm=1.
       F(:,:,:)=F(:,:,:)/coef_norm

       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX

                iglob = ibool(i,j,k,ispec)

                tempx1l = 0.
                tempx2l = 0.
                tempx3l = 0.

                do l=1,NGLLX

                   hp1 = hprime_xx(i,l)
                   tempx1l = tempx1l + F(l,j,k)*hp1

                   hp2 = hprime_xx(j,l)
                   tempx2l = tempx2l + F(i,l,k)*hp2

                   hp3 = hprime_xx(k,l)
                   tempx3l = tempx3l + F(i,j,l)*hp3

                enddo

                if (ispec_irreg /= 0) then
                  !irregular element
                  xixl = xixstore(i,j,k,ispec_irreg)
                  xiyl = xiystore(i,j,k,ispec_irreg)
                  xizl = xizstore(i,j,k,ispec_irreg)
                  etaxl = etaxstore(i,j,k,ispec_irreg)
                  etayl = etaystore(i,j,k,ispec_irreg)
                  etazl = etazstore(i,j,k,ispec_irreg)
                  gammaxl = gammaxstore(i,j,k,ispec_irreg)
                  gammayl = gammaystore(i,j,k,ispec_irreg)
                  gammazl = gammazstore(i,j,k,ispec_irreg)

                  !dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                  dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                  !dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

                else
                  !regular element

                  !dF(1,i,j,k) = xix_regular * tempx1l
                  dF(2,i,j,k) = xix_regular * tempx2l
                  !dF(3,i,j,k) = xix_regular * tempx3l

                endif

                !if (DEBUG_MODE) write(IIDD,*) 'y', dF(2,i,j,k)
                LapF(1,iglob) = LapF(1,iglob) + dF(2,i,j,k) * coef_norm
             enddo
          enddo
       enddo

       !! deivate again df/dz / dz
       F(:,:,:)=dF(3,:,:,:)

       coef_norm = maxval(abs(F(:,:,:)))
       !coef_norm=1.
       F(:,:,:)=F(:,:,:)/coef_norm

       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX

                iglob = ibool(i,j,k,ispec)

                tempx1l = 0.
                tempx2l = 0.
                tempx3l = 0.

                do l=1,NGLLX

                   hp1 = hprime_xx(i,l)
                   tempx1l = tempx1l + F(l,j,k)*hp1

                   hp2 = hprime_xx(j,l)
                   tempx2l = tempx2l + F(i,l,k)*hp2

                   hp3 = hprime_xx(k,l)
                   tempx3l = tempx3l + F(i,j,l)*hp3

                enddo
                if (ispec_irreg /= 0) then
                  !irregular element
                  xixl = xixstore(i,j,k,ispec_irreg)
                  xiyl = xiystore(i,j,k,ispec_irreg)
                  xizl = xizstore(i,j,k,ispec_irreg)
                  etaxl = etaxstore(i,j,k,ispec_irreg)
                  etayl = etaystore(i,j,k,ispec_irreg)
                  etazl = etazstore(i,j,k,ispec_irreg)
                  gammaxl = gammaxstore(i,j,k,ispec_irreg)
                  gammayl = gammaystore(i,j,k,ispec_irreg)
                  gammazl = gammazstore(i,j,k,ispec_irreg)

                  !dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                  !dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                  dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

                else
                  !regular element

                  !dF(1,i,j,k) = xix_regular * tempx1l
                  !dF(2,i,j,k) = xix_regular * tempx2l
                  dF(3,i,j,k) = xix_regular * tempx3l

                endif

                !if (DEBUG_MODE) write(IIDD,*) 'z', dF(3,i,j,k)
                LapF(1,iglob) = LapF(1,iglob) + dF(3,i,j,k) * coef_norm
                if (DEBUG_MODE .and. ispec == ispec_to_debug .and. i == 1 .and. j == 1) write(IIDD,*)   &
                     k, zstore(iglob), dF(3,i,j,k)*coef_norm,  F(i,j,k)*coef_norm
             enddo
          enddo
       enddo

       !!!!!!!!!!!!!!!!! HERE BEGINS DEBUG
!!$        do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                F(i,j,k)=LapF(1,iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       coef_norm = maxval(abs(F(:,:,:)))
!!$       F(:,:,:)=F(:,:,:)/coef_norm
!!$
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$                dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l  ! df/dx
!!$                dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l  ! df/dy
!!$                dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l  ! df/dz
!!$                if (DEBUG_MODE .and. ispec==1 .and. i==1 .and. j==1) write(IIDD,*)   k, dF(3,i,j,k),  F(i,j,k)
!!$                LapF(1,iglob) = dF(3,i,j,k) * coef_norm
!!$             enddo
!!$          enddo
!!$       enddo

            !! deivate again df/dx / dx
!!$       F(:,:,:)=dF(1,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$                dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
!!$                !dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
!!$                !dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l
!!$                !if (DEBUG_MODE) write(IIDD,*) 'x', dF(1,i,j,k)
!!$                LapF(1,iglob) = dF(1,i,j,k)
!!$             enddo
!!$          enddo
!!$       enddo

              !! deivate again df/dy / dy
!!$       F(:,:,:)=dF(2,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$               ! dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
!!$                dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
!!$               ! dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l
!!$                !if (DEBUG_MODE) write(IIDD,*) 'y', dF(2,i,j,k)
!!$                LapF(1,iglob) = LapF(1,iglob) + dF(2,i,j,k)
!!$             enddo
!!$          enddo
!!$       enddo

      !! deivate again df/dz / dz
!!$       F(:,:,:)=dF(3,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$               ! dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
!!$               ! dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
!!$                dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l
!!$                !if (DEBUG_MODE) write(IIDD,*) 'z', dF(3,i,j,k)
!!$                LapF(1,iglob) = LapF(1,iglob) + dF(3,i,j,k)
!!$             enddo
!!$          enddo
!!$       enddo

       !!!!!!!!!!!!!!!!! END OF DEBUG


    enddo

  end subroutine compute_laplac_lagrange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! for each GLL in element compute gradient norm and laplacian using lagrange interpolation and derivatives
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_grad_laplac_lagrange(nGrad, Lapf, field_to_derivate)

    implicit none

    real(kind=CUSTOM_REAL), dimension(:),     allocatable, intent(inout)  :: field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:),     allocatable, intent(inout)  :: nGrad, Lapf

    integer                                                       :: ispec, iglob, i, j, k

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),    allocatable        :: Derivatives_of_field
    real(kind=CUSTOM_REAL), dimension(:,:),          allocatable        :: field_to_derivate_wks, Fwks

    integer :: ier

    allocate(Derivatives_of_field(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 198')
    allocate(field_to_derivate_wks(NDIM,NGLOB_AB), Fwks(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 199')

    call compute_first_derivatives_lagrange(Derivatives_of_field, field_to_derivate)

    call compute_mean_values_on_edge(Derivatives_of_field)

    ! store norm L2 of gradient
    do ispec=1, NSPEC_AB
       do k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                nGrad(iglob) = sqrt(Derivatives_of_field(1,i,j,k,ispec)**2 + Derivatives_of_field(2,i,j,k,ispec)**2 +&
                     Derivatives_of_field(3,i,j,k,ispec)**2)
                field_to_derivate_wks(1,iglob)=Derivatives_of_field(1,i,j,k,ispec)
                field_to_derivate_wks(2,iglob)=Derivatives_of_field(2,i,j,k,ispec)
                field_to_derivate_wks(3,iglob)=Derivatives_of_field(3,i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    !! Df/dx 2nd derivatives
    field_to_derivate(:)=field_to_derivate_wks(1,:)
    call compute_first_derivatives_lagrange(Derivatives_of_field, field_to_derivate)
    !! store result
    do ispec=1, NSPEC_AB
       do k=1, NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                !LapF(iglob) =  Derivatives_of_field(1,i,j,k,ispec)
                Fwks(1,iglob) =  Derivatives_of_field(1,i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    !! Df/dy 2nd derivatives
    field_to_derivate(:)=field_to_derivate_wks(2,:)
    call compute_first_derivatives_lagrange(Derivatives_of_field, field_to_derivate)
    !! store result
    do ispec=1, NSPEC_AB
       do k=1, NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                !LapF(iglob) =   LapF(iglob) + Derivatives_of_field(2,i,j,k,ispec)
                !LapF(iglob) =  Derivatives_of_field(2,i,j,k,ispec)
                Fwks(2,iglob) =  Derivatives_of_field(2,i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo


    !! Df/dz 2nd derivatives
    field_to_derivate(:)=field_to_derivate_wks(3,:)
    call compute_first_derivatives_lagrange(Derivatives_of_field, field_to_derivate)
    !! store result
    do ispec=1, NSPEC_AB
       do k=1, NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob = ibool(i,j,k,ispec)
                !LapF(iglob) =   LapF(iglob) + Derivatives_of_field(3,i,j,k,ispec)
                !LapF(iglob) = Derivatives_of_field(3,i,j,k,ispec)
                Fwks(3,iglob) =  Derivatives_of_field(3,i,j,k,ispec)
             enddo
          enddo
       enddo
    enddo

    LapF(:) = Fwks(1,:) + Fwks(2,:) + Fwks(3,:)

    deallocate(Derivatives_of_field)
    deallocate(field_to_derivate_wks,Fwks)

!!$    do ispec =1, NSPEC_AB
!!$
!!$       !! ----------------------------------------------- 1st derivatives -----------------------
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$                iglob=ibool(i,j,k,ispec)
!!$                F(i,j,k)=field_to_derivate(iglob)
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$                dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l  ! df/dx
!!$                dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l  ! df/dy
!!$                dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l  ! df/dz
!!$
!!$                nGrad(iglob) =  sqrt(dF(1,i,j,k)**2 + dF(2,i,j,k)**2 + dF(3,i,j,k)**2)
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       !! --------------------------------- 2nd derivatives ----------------------------------
!!$       F(:,:,:)=dF(1,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$                dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
!!$                LapF(iglob) = dF(1,i,j,k)
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       F(:,:,:)=dF(2,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$
!!$                dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
!!$                LapF(iglob) = LapF(iglob) + dF(2,i,j,k)
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$       F(:,:,:)=dF(3,:,:,:)
!!$       do k=1,NGLLZ
!!$          do j=1,NGLLY
!!$             do i=1,NGLLX
!!$
!!$                iglob = ibool(i,j,k,ispec)
!!$
!!$                tempx1l = 0.
!!$                tempx2l = 0.
!!$                tempx3l = 0.
!!$
!!$                do l=1,NGLLX
!!$
!!$                   hp1 = hprime_xx(i,l)
!!$                   tempx1l = tempx1l + F(l,j,k)*hp1
!!$
!!$                   hp2 = hprime_xx(j,l)
!!$                   tempx2l = tempx2l + F(i,l,k)*hp2
!!$
!!$                   hp3 = hprime_xx(k,l)
!!$                   tempx3l = tempx3l + F(i,j,l)*hp3
!!$
!!$                enddo
!!$
!!$                xixl = xixstore(i,j,k,ispec)
!!$                xiyl = xiystore(i,j,k,ispec)
!!$                xizl = xizstore(i,j,k,ispec)
!!$                etaxl = etaxstore(i,j,k,ispec)
!!$                etayl = etaystore(i,j,k,ispec)
!!$                etazl = etazstore(i,j,k,ispec)
!!$                gammaxl = gammaxstore(i,j,k,ispec)
!!$                gammayl = gammaystore(i,j,k,ispec)
!!$                gammazl = gammazstore(i,j,k,ispec)
!!$
!!$                dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l
!!$                LapF(iglob) = LapF(iglob) + dF(3,i,j,k)
!!$
!!$             enddo
!!$          enddo
!!$       enddo
!!$
!!$    enddo

  end subroutine compute_grad_laplac_lagrange

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! for each GLL in element faces compute and store FD derivative matrix
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine compute_and_store_FD_derivatives_matrix(regularization_fd)

!   implicit none
!   type(regul), dimension(:), allocatable, intent(inout)  :: regularization_fd

!   ! locals
!   real(kind=CUSTOM_REAL), dimension(:),       allocatable :: xstore_recv, ystore_recv, zstore_recv
!   double precision,       dimension(:),       allocatable :: line_vander
!   double precision,       dimension(:,:),     allocatable :: VanderMonde, MtM, invMtM, FDm_dp, MtM_tmp

!   integer                                                 :: iglob_index, ig, igll, iglob
!   integer                                                 :: nline, ncolu, ip, norder
!   real(kind=CUSTOM_REAL)                                  :: xp, yp, zp
!   real(kind=CUSTOM_REAL)                                  :: xs, ys, zs

!   !! communication of grid from neighbor MPI slices
!   allocate(xstore_recv(indx_recv(NPROC)), ystore_recv(indx_recv(NPROC)), zstore_recv(indx_recv(NPROC)))
!   call send_recv_blocking(xstore, xstore_recv)
!   call send_recv_blocking(ystore, ystore_recv)
!   call send_recv_blocking(zstore, zstore_recv)

!   !! for order 2
!   allocate(MtM(10,10), invMtM(10,10), MtM_tmp(10,10))

!   do iglob_index = 1, Nb_iglob_on_faces

!      ! position of point
!      ig= regularization_fd(iglob_index)%iglob
!      xp=xstore(ig)
!      yp=ystore(ig)
!      zp=zstore(ig)

!      nline=10
!      ncolu=regularization_fd(iglob_index)%nReg+regularization_fd(iglob_index)%nNei
!      allocate(VanderMonde(ncolu, nline), line_vander(nline), FDm_dp(nline, ncolu))
!      VanderMonde(:,:)=0.

!      ip=0
!      do igll = 1, regularization_fd(iglob_index)%nReg
!         iglob=regularization_fd(iglob_index)%iglob_regular_point_to_use(igll)
!         xs=xstore(iglob)
!         ys=ystore(iglob)
!         zs=zstore(iglob)
!         norder=regularization_fd(iglob_index)%MaxOrder
!         call line_vandermonde(line_vander, xp,yp,zp, xs, ys, zs, norder)
!         ip=ip+1
!         VanderMonde(ip,:)=line_vander(:)
!      enddo

!      do igll = 1, regularization_fd(iglob_index)%nNei
!         iglob = regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
!         xs = xstore_recv(iglob)
!         ys = ystore_recv(iglob)
!         zs = zstore_recv(iglob)
!         norder=regularization_fd(iglob_index)%MaxOrder
!         call line_vandermonde(line_vander, xp,yp,zp, xs, ys, zs, norder)
!         ip=ip+1
!         VanderMonde(ip,:)=line_vander(:)
!      enddo

!      call compute_M_transpose_times_M(MtM, VanderMonde, nline, ncolu)
!      invMtM(:,:)=0.d0
!      MtM_tmp(:,:)=MtM(:,:)
!      call inverse(MtM,invMtM,nline)

!      call compute_inv_MtM_times_Mt(Fdm_dp, invMtM, VanderMonde, nline, ncolu)
!      allocate(regularization_fd(iglob_index)%Deriv_FD_Matrix(nline, ncolu))

!      mem_used_for_reg =  mem_used_for_reg + CUSTOM_REAL*nline*ncolu
!      mem_tmp =  mem_tmp +  CUSTOM_REAL*nline*ncolu
!      max_memory_used = max(max_memory_used, mem_used_for_reg)

!      !! store derivatice matrix
!      regularization_fd(iglob_index)%Deriv_FD_Matrix(:,:)=FDm_dp(:,:)

!      deallocate(VanderMonde, line_vander, FDm_dp)
!   enddo
!   deallocate(xstore_recv, ystore_recv, zstore_recv)

! end subroutine compute_and_store_FD_derivatives_matrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute laplacian of field_input on elements faces using FD
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_laplacian_FD(Dfb, field_input, regularization_fd)


    implicit none
    type(regul),            dimension(:),   allocatable, intent(in)     :: regularization_fd
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable, intent(in)     :: field_input
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable, intent(inout)  :: Dfb

    real(kind=CUSTOM_REAL), dimension(:,:), allocatable                 :: valence, field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:),   allocatable                 :: field_to_send, field_overlap
    real(kind=CUSTOM_REAL), dimension(:),   allocatable                 :: result_df, Values
    integer                                                             :: iglob, iglob_index, igll, idim, ip
    integer                                                             :: nline, ncolu

    integer :: ier

    nline=10

    allocate(valence(NDIM,NGLOB_AB),  field_to_derivate(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 200')
    allocate(field_to_send(NGLOB_AB), field_overlap(indx_recv(NPROC)), result_df(nline),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 201')
    valence(:,:) = 1.
    field_to_derivate(:,:) = field_input(:,:)

    !! 1/ compute average values of field_to_derivate at the edge of MPI slices -------------------
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,valence, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,field_to_derivate, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    !! divide by valence
    field_to_derivate(:,:) = field_to_derivate(:,:) / valence(:,:)

    do idim = 1, NDIM
       !! 2/ communicate overlap of MPI slices
       field_to_send(:)=field_to_derivate(idim,:)
       call send_recv_blocking(field_to_send, field_overlap)

       !! 3/ process derivatives in edges by FD
       do iglob_index=1, Nb_iglob_on_faces
          ncolu=regularization_fd(iglob_index)%nReg+regularization_fd(iglob_index)%nNei
          allocate(Values(ncolu),stat=ier)
          if (ier /= 0) call exit_MPI_without_rank('error allocating array 202')
          ip=0
          do igll=1, regularization_fd(iglob_index)%nReg
             ip = ip + 1
             iglob=regularization_fd(iglob_index)%iglob_regular_point_to_use(igll)
             Values(ip)=field_to_send(iglob)
          enddo


          do igll = 1,  regularization_fd(iglob_index)%nNei
             ip=ip+1
             iglob = regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
             Values(ip) = field_overlap(iglob)
          enddo

          !! 4 apply derivative matrix
          call matrix_times_vector(result_df, regularization_fd(iglob_index)%Deriv_FD_Matrix, Values, nline, ncolu)

          deallocate(Values)

          Dfb(idim,iglob_index) = result_df(5) + result_df(6) + result_df(7)

       enddo


    enddo

    deallocate(valence, field_to_send, field_overlap, result_df, field_to_derivate)

  end subroutine compute_laplacian_FD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
!  compute norm gradient and laplacian of field_input on elements faces using FD
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine compute_gradient_laplacian_FD(nGrad, Laplac, field_input, regularization_fd)

    implicit none
    type(regul),            dimension(:),   allocatable, intent(in)     :: regularization_fd
    real(kind=CUSTOM_REAL), dimension(:),   allocatable, intent(in)     :: field_input
    real(kind=CUSTOM_REAL), dimension(:),   allocatable, intent(inout)  :: Laplac, nGrad

    real(kind=CUSTOM_REAL), dimension(:,:), allocatable                 :: valence, field_to_derivate
    real(kind=CUSTOM_REAL), dimension(:),   allocatable                 :: field_to_send, field_overlap
    real(kind=CUSTOM_REAL), dimension(:),   allocatable                 :: result_df, Values
    integer                                                             :: iglob, iglob_index, igll, ip
    integer                                                             :: nline, ncolu

    integer :: ier

    nline=10
    allocate(valence(NDIM,NGLOB_AB),  field_to_derivate(NDIM,NGLOB_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 203')
    allocate(field_to_send(NGLOB_AB), field_overlap(indx_recv(NPROC)), result_df(nline),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 204')
    valence(:,:) = 1.

    !! need to duplicate in order to use the already build subroutine from sepcfem package : assembel_MPI..
    field_to_derivate(1,:) = field_input(:)
    field_to_derivate(2,:) = field_input(:)
    field_to_derivate(3,:) = field_input(:)

    !! 1/ compute average values of field_to_derivate at the edge of MPI slices -------------------
    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,valence, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    call assemble_MPI_vector_blocking(NPROC,NGLOB_AB,field_to_derivate, &
                                      num_interfaces_ext_mesh_sp,max_nibool_interfaces_ext_mesh_sp, &
                                      nibool_interfaces_ext_mesh_sp,ibool_interfaces_ext_mesh_sp, &
                                      my_neighbors_ext_mesh_sp)

    !! divide by valence
    field_to_derivate(1,:) = field_to_derivate(1,:) / valence(1,:)


    !! 2/ communicate overlap of MPI slices
    field_to_send(:)=field_to_derivate(1,:)
    call send_recv_blocking(field_to_send, field_overlap)

    !! 3/ process derivatives in edges by FD
    do iglob_index=1, Nb_iglob_on_faces
       ncolu=regularization_fd(iglob_index)%nReg+regularization_fd(iglob_index)%nNei
       allocate(Values(ncolu),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 205')
       ip=0
       do igll=1, regularization_fd(iglob_index)%nReg
          ip = ip + 1
          iglob=regularization_fd(iglob_index)%iglob_regular_point_to_use(igll)
          Values(ip)=field_to_send(iglob)
       enddo


       do igll = 1,  regularization_fd(iglob_index)%nNei
          ip=ip+1
          iglob = regularization_fd(iglob_index)%iglob_neighbo_point_to_use(igll)
          Values(ip) = field_overlap(iglob)
       enddo

          !! 4 apply derivative matrix
       call matrix_times_vector(result_df, regularization_fd(iglob_index)%Deriv_FD_Matrix, Values, nline, ncolu)

       deallocate(Values)

       Laplac(iglob_index) = result_df(5) + result_df(6) + result_df(7)
       nGrad(iglob_index) = result_df(2)**2 + result_df(3)**2 + result_df(4)**2

    enddo

    deallocate(valence, field_to_send, field_overlap, result_df, field_to_derivate)

  end subroutine compute_gradient_laplacian_FD

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATHEMATICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! compute line of vandermonde matrix : limited devloppement
! f(u,v,w) =  f(x,y,z) + (u-x)*dxf(x,y,z) + (v-y)*dyf(x,y,z) + ...
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine line_vandermonde(l, x,y,z, u,v,w, n)

!   implicit none

!   integer, intent(in) :: n
!   real(kind=CUSTOM_REAL), intent(in) :: x,y,z, u,v,w
!   double precision,   dimension(:), allocatable :: l

!   integer :: nmax

!   if (n /= 2 ) then
!      stop 'ABORT: Vandermonde matrix defined only for order = 2 !!'
!   endif

!   nmax = (n+1)*(n+2)*(n+3) / 6

!   l(1)=1
!   l(2)=u-x
!   l(3)=v-y
!   l(4)=w-z
!   l(5)=0.5*(u-x)**2
!   l(6)=0.5*(v-y)**2
!   l(7)=0.5*(w-z)**2
!   l(8)=(u-x)*(v-y)
!   l(9)=(u-x)*(w-z)
!   l(10)=(v-y)*(w-z)

! end subroutine line_vandermonde

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! compute matrix product tranpose(V)*V and store in Mt
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine  compute_M_transpose_times_M(Mt,V,n,m)

!   implicit none
!   double precision, dimension(:,:), allocatable, intent(in)    :: V
!   double precision, dimension(:,:), allocatable, intent(inout) :: Mt
!   integer,                                       intent(in)    :: n,m

!   integer                                                      :: i, j, k

!   Mt(:,:)=0.

!   do i=1,n
!      do j=1,n
!         do k=1,m
!            Mt(i,j) = Mt(i,j) + V(k,i)*V(k,j)
!         enddo
!      enddo
!   enddo

! end subroutine compute_M_transpose_times_M

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! compute matrix product iMtM*tranpose(V) and store in DFm
!---------------------------------------------------------------------------------
!##################################################################################################################################

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine compute_inv_MtM_times_Mt(FDm, iMtM, V, n, m)

!   implicit none
!   integer,                                       intent(in)    :: n, m
!   double precision, dimension(:,:), allocatable, intent(in)    :: iMtM, V
!   double precision, dimension(:,:), allocatable, intent(inout) :: FDm

!   integer                                                      :: i,j,k

!   FDm(:,:)=0.d0

!   do i = 1, n
!      do j = 1, m
!         do k=1, n
!            FDm(i,j) = FDm(i,j) + iMtM(i,k)*V(j,k)
!         enddo
!      enddo
!   enddo

! end subroutine compute_inv_MtM_times_Mt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!---------------------------------------------------------------------------------
! matrix vector mulitplication r = A*b
!---------------------------------------------------------------------------------
!##################################################################################################################################

  subroutine matrix_times_vector(r, A, b, n, m)

    implicit none
    integer,                                              intent(in)    :: n,m
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable,  intent(in)    :: A
    real(kind=CUSTOM_REAL), dimension(:),   allocatable,  intent(in)    :: b
    real(kind=CUSTOM_REAL), dimension(:),   allocatable,  intent(inout) :: r

    integer                                                             :: i,k


    do i=1,n
       r(i)=0.
       do k=1,m
          r(i)=r(i)+A(i,k)*b(k)
       enddo
    enddo

  end subroutine matrix_times_vector

!!########################################################################################################
!!  inverse matrix program taken in web page of
!!  Alexander Godunov : http://ww2.odu.edu/~agodunov/
!!  http://ww2.odu.edu/~agodunov/computing/programs/index.html
!!  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
!!
!! codes from book :
!!                   Introductory Computational Physics
!!                   Andi Klein and Alexander Godunov
!!                   Cambridge University Press 2006
!!                   ISBN 0-521-82862-7
!!
!!
!!

!! DK DK commented out because this subroutine is currently unused, and gfortran -Wall gives an error in case of unused functions

! subroutine inverse(a,c,n)
!   !============================================================
!   ! Inverse matrix
!   ! Method: Based on Doolittle LU factorization for Ax=b
!   ! Alex G. December 2009
!   !-----------------------------------------------------------
!   ! input ...
!   ! a(n,n) - array of coefficients for matrix A
!   ! n      - dimension
!   ! output ...
!   ! c(n,n) - inverse matrix of A
!   ! comments ...
!   ! the original matrix a(n,n) will be destroyed
!   ! during the calculation
!   !===========================================================

!   !! little modif by VM

!   implicit none
!   integer,                                       intent(in)     :: n
!   double precision, dimension(:,:), allocatable, intent(inout)  :: a,c     ! a(n,n), c(n,n)
!   double precision, dimension(:,:), allocatable                 :: L, U
!   double precision, dimension(:),   allocatable                 :: b, d, x
!   double precision                                              :: coeff
!   integer                                                       :: i, j, k

!   ! step 0: initialization for matrices L and U and b
!   allocate(L(n,n), U(n,n))
!   allocate(b(n), d(n), x(n))

!   L=0.0
!   U=0.0
!   b=0.0

!   ! step 1: forward elimination
!   do k=1, n-1
!      do i=k+1,n
!         coeff=a(i,k)/a(k,k)
!         L(i,k) = coeff
!         do j=k+1,n
!            a(i,j) = a(i,j)-coeff*a(k,j)
!         enddo
!      enddo
!   enddo

!   ! Step 2: prepare L and U matrices
!   ! L matrix is a matrix of the elimination coefficient
!   ! + the diagonal elements are 1.0
!   do i=1,n
!      L(i,i) = 1.0
!   enddo
!   ! U matrix is the upper triangular part of A
!   do j=1,n
!      do i=1,j
!         U(i,j) = a(i,j)
!      enddo
!   enddo

!   ! Step 3: compute columns of the inverse matrix C
!   do k=1,n
!      b(k)=1.0
!      d(1) = b(1)
!      ! Step 3a: Solve Ld=b using the forward substitution
!      do i=2,n
!         d(i)=b(i)
!         do j=1,i-1
!            d(i) = d(i) - L(i,j)*d(j)
!         enddo
!      enddo
!      ! Step 3b: Solve Ux=d using the back substitution
!      x(n)=d(n)/U(n,n)
!      do i = n-1,1,-1
!         x(i) = d(i)
!         do j=n,i+1,-1
!            x(i)=x(i)-U(i,j)*x(j)
!         enddo
!         x(i) = x(i)/u(i,i)
!      enddo
!      ! Step 3c: fill the solutions x(n) into column k of C
!      do i=1,n
!         c(i,k) = x(i)
!      enddo
!      b(k)=0.0
!   enddo

!   deallocate(L,U, b, d, x)

! end subroutine inverse

!!!!!!!!!!!!!!!! DEBUG subroutine !!!!!!!!!!!!!!!!!!
  subroutine write_in_disk_this(f)

    integer :: ier

    real(kind=CUSTOM_REAL), dimension(:), allocatable :: f
    integer i,j,k,ispec, iglob
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dd
    character(len=256)                                                     :: path_file, name_file
    integer itest
    itest=1
    allocate(dd(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 206')
    do ispec=1,nspec_ab
       do k=1,ngllz
          do j=1,nglly
             do i=1,ngllx
                iglob = ibool(i,j,k,ispec)
                dd(i,j,k,ispec)=f(iglob)
             enddo
          enddo
       enddo
    enddo

    path_file='OUTPUT_FILES/DATABASES_MPI/proc'
    write(name_file,'(i6.6,a8,i2.2,a4)') myrank, '_lalala_',itest,'.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) dd
    close(888)
    deallocate(dd)

  end subroutine write_in_disk_this




!!!========================== NEW WAY TO DERIVATE ===================

  subroutine compute_lapalacian_of_field(field, laplacian_of_field)

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),      allocatable, intent(in)    :: field
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),      allocatable, intent(inout) :: laplacian_of_field

  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),    allocatable                :: field_wkstmp
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),    allocatable                :: derivative_of_field
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),    allocatable                :: second_derivative_of_field

  integer :: ier

  allocate(field_wkstmp(3,NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 207')
  allocate(derivative_of_field(3,NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 208')
  allocate(second_derivative_of_field(3,NGLLX, NGLLY, NGLLZ, NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 209')

  field_wkstmp(1,:,:,:,:) = field(:,:,:,:)
  field_wkstmp(2,:,:,:,:) = field(:,:,:,:)
  field_wkstmp(3,:,:,:,:) = field(:,:,:,:)
  call compute_mean_values_on_edge(field_wkstmp)
  call compute_derivative_with_lagrange_polynomials(derivative_of_field, field_wkstmp)

  call compute_mean_values_on_edge(derivative_of_field)
  call compute_2nd_derivative_with_lagrange_polynomials(second_derivative_of_field, derivative_of_field)

  call compute_mean_values_on_edge(second_derivative_of_field)

  laplacian_of_field(:,:,:,:) = second_derivative_of_field(1,:,:,:,:) + &
                                second_derivative_of_field(2,:,:,:,:) + &
                                second_derivative_of_field(3,:,:,:,:)

  deallocate(field_wkstmp, derivative_of_field, second_derivative_of_field)

  end subroutine compute_lapalacian_of_field

!===================

  subroutine compute_bi_laplacian_of_field(field, laplacian_of_field, bi_laplacian_of_field)

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),      allocatable, intent(in)    :: field
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),      allocatable, intent(inout) :: laplacian_of_field, bi_laplacian_of_field

  call compute_lapalacian_of_field(field, laplacian_of_field)
  call compute_lapalacian_of_field(laplacian_of_field, bi_laplacian_of_field)

  end subroutine compute_bi_laplacian_of_field


!!===============
!!===============
!!===============

  subroutine compute_derivative_with_lagrange_polynomials(derivative_of_field, field_to_derivate)

  use specfem_par, only: xixstore, xiystore, xizstore, etaxstore, etaystore, etazstore, &
                         gammaxstore, gammaystore, gammazstore, hprime_xx, irregular_element_number, &
                         xix_regular

  implicit none
  !! size (NDIM,NGLLX, NGLLY, NGLLZ, NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),     allocatable, intent(in)     :: field_to_derivate
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),     allocatable, intent(inout)  :: derivative_of_field


  double precision, dimension(NGLLX, NGLLY, NGLLZ)              :: F
  double precision, dimension(NDIM, NGLLX, NGLLY, NGLLZ)        :: DF
  double precision                                              :: hp1, hp2, hp3
  double precision                                              :: tempx1l, tempx2l, tempx3l
  double precision                                              :: xixl,xiyl,xizl
  double precision                                              :: etaxl,etayl,etazl
  double precision                                              :: gammaxl,gammayl,gammazl
  integer                                                       :: ispec, ispec_irreg, i, j, k, l


  do ispec =1, NSPEC_AB

!! store field in local array
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          !iglob=ibool(i,j,k,ispec)
          F(i,j,k)=field_to_derivate(1,i,j,k,ispec)
        enddo
      enddo
    enddo
    ispec_irreg = irregular_element_number(ispec)
    !! derivative of field based on lagrange polynomials
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          !iglob = ibool(i,j,k,ispec)

          tempx1l = 0.
          tempx2l = 0.
          tempx3l = 0.

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1l = tempx1l + F(l,j,k)*hp1

            hp2 = hprime_xx(j,l)
            tempx2l = tempx2l + F(i,l,k)*hp2

            hp3 = hprime_xx(k,l)
            tempx3l = tempx3l + F(i,j,l)*hp3
          enddo

          if (ispec_irreg /= 0) then
            !irregular element
            xixl = xixstore(i,j,k,ispec_irreg)
            xiyl = xiystore(i,j,k,ispec_irreg)
            xizl = xizstore(i,j,k,ispec_irreg)
            etaxl = etaxstore(i,j,k,ispec_irreg)
            etayl = etaystore(i,j,k,ispec_irreg)
            etazl = etazstore(i,j,k,ispec_irreg)
            gammaxl = gammaxstore(i,j,k,ispec_irreg)
            gammayl = gammaystore(i,j,k,ispec_irreg)
            gammazl = gammazstore(i,j,k,ispec_irreg)

            dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
            dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
            dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          else
            !regular element
            dF(1,i,j,k) = xix_regular * tempx1l
            dF(2,i,j,k) = xix_regular * tempx2l
            dF(3,i,j,k) = xix_regular * tempx3l
          endif

        enddo
      enddo
    enddo

    !! get derivative of field from  local array
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          derivative_of_field(1,i,j,k,ispec)= dF(1,i,j,k)
          derivative_of_field(2,i,j,k,ispec)= dF(2,i,j,k)
          derivative_of_field(3,i,j,k,ispec)= dF(3,i,j,k)
        enddo
      enddo
    enddo

  enddo

  end subroutine compute_derivative_with_lagrange_polynomials

!!============

  subroutine compute_2nd_derivative_with_lagrange_polynomials(derivative_of_field, field_to_derivate)

   use specfem_par, only: xixstore, xiystore, xizstore, etaxstore, etaystore, etazstore, &
                          gammaxstore, gammaystore, gammazstore, hprime_xx, irregular_element_number, &
                          xix_regular

   implicit none
   !! size (NDIM, NGLLX, NGLLY, NGLLZ, NSPEC_AB)
   real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),     allocatable, intent(in)     :: field_to_derivate
   real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),     allocatable, intent(inout)  :: derivative_of_field


   double precision, dimension(NDIM, NGLLX, NGLLY, NGLLZ)        :: F
   double precision, dimension(NDIM, NGLLX, NGLLY, NGLLZ)        :: DF
   double precision                                              :: hp1, hp2, hp3
   double precision                                              :: tempx1l, tempx2l, tempx3l
   double precision                                              :: tempy1l, tempy2l, tempy3l
   double precision                                              :: tempz1l, tempz2l, tempz3l
   double precision                                              :: xixl,xiyl,xizl
   double precision                                              :: etaxl,etayl,etazl
   double precision                                              :: gammaxl,gammayl,gammazl
   integer                                                       :: ispec, ispec_irreg,i, j, k, l

   do ispec =1, NSPEC_AB

      !! store field in local array
      do k=1,NGLLZ
         do j=1,NGLLY
            do i=1,NGLLX
               F(1,i,j,k)=field_to_derivate(1,i,j,k,ispec)
               F(2,i,j,k)=field_to_derivate(2,i,j,k,ispec)
               F(3,i,j,k)=field_to_derivate(3,i,j,k,ispec)
            enddo
         enddo
      enddo
      ispec_irreg = irregular_element_number(ispec)

      !! derivative of field based on lagrange polynomials
      do k=1,NGLLZ
         do j=1,NGLLY
            do i=1,NGLLX

               !iglob = ibool(i,j,k,ispec)

               tempx1l = 0.
               tempx2l = 0.
               tempx3l = 0.

               tempy1l = 0.
               tempy2l = 0.
               tempy3l = 0.

               tempz1l = 0.
               tempz2l = 0.
               tempz3l = 0

               do l=1,NGLLX

                  hp1 = hprime_xx(i,l)
                  tempx1l = tempx1l + F(1,l,j,k)*hp1
                  tempy1l = tempy1l + F(2,l,j,k)*hp1
                  tempz1l = tempz1l + F(3,l,j,k)*hp1

                  hp2 = hprime_xx(j,l)
                  tempx2l = tempx2l + F(1,i,l,k)*hp2
                  tempy2l = tempy2l + F(2,i,l,k)*hp2
                  tempz2l = tempz2l + F(3,i,l,k)*hp2

                  hp3 = hprime_xx(k,l)
                  tempx3l = tempx3l + F(1,i,j,l)*hp3
                  tempy3l = tempy3l + F(2,i,j,l)*hp3
                  tempz3l = tempz3l + F(3,i,j,l)*hp3


               enddo
               if (ispec_irreg /= 0) then !irregular element

                 xixl = xixstore(i,j,k,ispec_irreg)
                 xiyl = xiystore(i,j,k,ispec_irreg)
                 xizl = xizstore(i,j,k,ispec_irreg)
                 etaxl = etaxstore(i,j,k,ispec_irreg)
                 etayl = etaystore(i,j,k,ispec_irreg)
                 etazl = etazstore(i,j,k,ispec_irreg)
                 gammaxl = gammaxstore(i,j,k,ispec_irreg)
                 gammayl = gammaystore(i,j,k,ispec_irreg)
                 gammazl = gammazstore(i,j,k,ispec_irreg)

                 dF(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
                 dF(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
                 dF(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

               else !regular element

                 dF(1,i,j,k) = xix_regular * tempx1l
                 dF(2,i,j,k) = xix_regular * tempx2l
                 dF(3,i,j,k) = xix_regular * tempx3l

               endif

          enddo
       enddo
    enddo

    !! get derivative of field from  local array
    do k=1,NGLLZ
       do j=1,NGLLY
          do i=1,NGLLX
             derivative_of_field(1,i,j,k,ispec)= dF(1,i,j,k)
             derivative_of_field(2,i,j,k,ispec)= dF(2,i,j,k)
             derivative_of_field(3,i,j,k,ispec)= dF(3,i,j,k)
          enddo
       enddo
    enddo

 enddo

  end subroutine compute_2nd_derivative_with_lagrange_polynomials

!!
!! variable damping regularization for trying to kill suprious variations close to the point sources

  subroutine compute_spatial_damping_for_source_singularities(acqui_simu, inversion_param, spatial_damping)

  type(inver),                                                    intent(inout) :: inversion_param
  type(acqui),            dimension(:),       allocatable,        intent(inout) :: acqui_simu
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable,        intent(inout) :: spatial_damping
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                       :: spatial_damping_tmp
  integer                                                                       :: isrc, ievent, iglob, ispec, i, j, k
  real(kind=CUSTOM_REAL)                                                        :: xgll, ygll, zgll
  real(kind=CUSTOM_REAL)                                                        :: distance_from_source, value_of_damping

  integer :: ier

  do ispec = 1, NSPEC_AB

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob=ibool(i,j,k,ispec)
          xgll=xstore(iglob)
          ygll=ystore(iglob)
          zgll=zstore(iglob)

          !! compute distance form sources
          do ievent = 1, acqui_simu(1)%nevent_tot
            if (trim(acqui_simu(ievent)%source_type) == 'moment' .or. trim(acqui_simu(ievent)%source_type) == 'force' .or. &
                trim(acqui_simu(ievent)%source_type) == 'shot' ) then

              do isrc = 1, acqui_simu(ievent)%nsources_local
                distance_from_source = sqrt(  (acqui_simu(ievent)%Xs(isrc) - xgll)**2 + &
                                              (acqui_simu(ievent)%Ys(isrc) - ygll)**2 + &
                                              (acqui_simu(ievent)%Zs(isrc) - zgll)**2)

                value_of_damping = inversion_param%min_damp + &
                                   (inversion_param%max_damp - inversion_param%min_damp )*&
                      exp(-0.5 * (distance_from_source/(inversion_param%distance_from_source/3.))**2 )

                spatial_damping(i,j,k,ispec) = max(spatial_damping(i,j,k,ispec), value_of_damping)
                !write(*,*) inversion_param%max_damp,  inversion_param%min_damp,
                !inversion_param%distance_from_source, value_of_damping
              enddo
            endif
          enddo

        enddo
      enddo
    enddo
  enddo

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
    allocate(spatial_damping_tmp(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 210')
    spatial_damping_tmp(:,:,:,:)=spatial_damping(:,:,:,:)

    call max_all_all_cr_for_simulatenous_runs(spatial_damping_tmp(1,1,1,1), spatial_damping(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

    deallocate(spatial_damping_tmp)
  endif

  end subroutine compute_spatial_damping_for_source_singularities

end module regularization
