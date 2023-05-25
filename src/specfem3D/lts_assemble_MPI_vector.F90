!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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


! MPI assembly routines for LTS

  subroutine assemble_MPI_vector_async_send_lts(NPROC,NGLOB_AB,array_val,ilevel, &
                                                buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                                num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                                my_neighbors_ext_mesh, &
                                                request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! sends data on LTS level

  use constants, only: NDIM,CUSTOM_REAL,itag,ASSEMBLE_MPI_OFF

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_elastic, only: nspec_outer_elastic

  ! LTS
  use specfem_par_lts, only: num_interface_p_refine_ibool, interface_p_refine_ibool, num_p_level

  ! GPU
  !use specfem_par_lts, only: num_interface_p_refine_boundary

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer,intent(in) :: ilevel

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  ! local parameters
  integer :: ipoin,iinterface,iglob
  integer :: num_buffer_points
  ! testing
  !logical, parameter :: TEST_GPU = .false.
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: test_buffer
  !integer :: loc(3)
  !integer :: ier

  ! here we have to assemble all the contributions between partitions using MPI

  ! debug: no mpi
  if (ASSEMBLE_MPI_OFF) return

  ! assemble only if more than one partition
  if (NPROC == 1) return

  ! sends MPI buffers
  if (.not. GPU_MODE) then
    ! on CPU
    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)

      if (num_buffer_points > 0) then
        ! fills buffer
        do ipoin = 1, num_buffer_points
          iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
          buffer_send_vector_ext_mesh(:,ipoin,iinterface) = array_val(:,iglob)
        enddo
      endif
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)

      if (num_buffer_points > 0) then
        call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
                      NDIM*num_buffer_points, &
                      my_neighbors_ext_mesh(iinterface), &
                      itag, &
                      request_send_vector_ext_mesh(iinterface))
        call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
                      NDIM*num_buffer_points, &
                      my_neighbors_ext_mesh(iinterface), &
                      itag, &
                      request_recv_vector_ext_mesh(iinterface))
      endif
    enddo

  else ! GPU_MODE

    if (ilevel == num_p_level) then
      ! coarse level does standard xfer
      call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)

    else
      ! ilevel < num_p_level (fine levels)

      !#TODO: LTS on GPU stacey contribution
      stop 'LTS on GPU w/ assemble MPI send ilevel contribution not implemented yet'

      ! transfers boundary region to host asynchronously. The
      ! MPI-send is done from within compute_forces_viscoelastic_cuda,
      ! once the inner element kernels are launched, and the
      ! memcpy has finished. see compute_forces_viscoelastic_cuda:1655

      !if (TEST_GPU) then
      !  allocate(test_buffer(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
      !  if (ier /= 0) stop 'Error allocating test buffer'
      !  test_buffer(:,:,:) = 0
      !  call test_boundary_transfer_lts(Mesh_pointer,ilevel,maxval(num_interface_p_refine_ibool),test_buffer)
      !  call copy_accelfield_from_gpu(Mesh_pointer,array_val)
      !  ! zero out for test
      !  buffer_send_vector_ext_mesh = 0
      !  do iinterface = 1, num_interfaces_ext_mesh
      !    if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
      !      do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
      !        iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
      !        buffer_send_vector_ext_mesh(:,ipoin,iinterface) = array_val(:,iglob)
      !      enddo
      !    endif
      !  enddo
      !  if (maxval(abs(test_buffer-buffer_send_vector_ext_mesh)) > 1e-10) then
      !    loc = maxloc(abs(test_buffer-buffer_send_vector_ext_mesh))
      !    print *, "test_buffer vs. buffer_send_vector_ext_mesh=",maxval(abs(test_buffer-buffer_send_vector_ext_mesh)), &
      !            "ilevel=",ilevel, "@", loc, "... ", test_buffer(loc(1),loc(2),loc(3)), &
      !            "vs.", buffer_send_vector_ext_mesh(loc(1),loc(2),loc(3))
      !  endif
      !  deallocate(test_buffer)
      !endif
      !
      ! only transfer nodes on ilevel:
      !call transfer_reduced_boundary_from_device_async_lts(Mesh_pointer,ilevel,num_interface_p_refine_boundary(ilevel))
      ! or not optimized, transfers full acceleration field:
      !call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)

    endif ! ilevel < num_p_level

  endif ! GPU


  end subroutine assemble_MPI_vector_async_send_lts

!
!-------------------------------------------------------------------------------------------------
!

! MPI routines that only recv if *this* p-level has elements on the MPI boundary

  subroutine assemble_MPI_vector_async_recv_lts(NPROC,NGLOB_AB,array_val,ilevel, &
                                                buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                                max_nibool_interfaces_ext_mesh, &
                                                nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                                request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                                my_neighbors_ext_mesh)

! waits for send/receiver to be completed and assembles contributions

  use constants, only: NDIM,CUSTOM_REAL,ASSEMBLE_MPI_OFF

  use specfem_par, only: FAULT_SIMULATION

  use specfem_par, only: GPU_MODE, Mesh_pointer, buffer_send_vector_ext_mesh

  use specfem_par_lts, only: num_interface_p_refine_ibool, interface_p_refine_ibool, num_p_level

  ! GPU
  !use specfem_par_lts, only: num_interface_p_refine_boundary, interface_p_refine_boundary

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer,intent(in) :: ilevel

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
    buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! local parameters
  integer :: ipoin,iinterface,iglob
  integer :: num_buffer_points
  ! GPU
  !real(kind=CUSTOM_REAL), dimension(:), allocatable :: reduced_buffer_recv_vector_ext_mesh, reduced_buffer_send_vector_ext_mesh
  !integer :: num_refine,index,ier
  ! testing
  !logical, parameter :: TEST_GPU = .false.
  !integer :: loc(3)
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh_test

  ! here we have to assemble all the contributions between partitions using MPI

  ! debug: no mpi
  if (ASSEMBLE_MPI_OFF) return

  ! assemble only if more than one partition
  if (NPROC == 1) return

  ! fault ruptures
  if (FAULT_SIMULATION) then
    ! receives MPI buffers with ordered assembly
    call assemble_MPI_vector_async_w_ord_lts(NPROC,NGLOB_AB,array_val,ilevel, &
                                             buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                             max_nibool_interfaces_ext_mesh, &
                                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                             request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                             my_neighbors_ext_mesh)
    ! all done
    return
  endif

  ! receives and assembles values
  if (.not. GPU_MODE) then
    ! on CPU
    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        call wait_req(request_recv_vector_ext_mesh(iinterface))
      endif
    enddo

    ! adding all contributions
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        do ipoin = 1, num_buffer_points
          iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
          array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
        enddo
      endif
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        call wait_req(request_send_vector_ext_mesh(iinterface))
      endif
    enddo

  else
    ! on GPU
    ! while Inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (ilevel == num_p_level) then
      ! default transfers
      ! coarse level needs full array
      call sync_copy_from_device(Mesh_pointer,2,buffer_send_vector_ext_mesh)

      ! sends MPI buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                                         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, &
                                         my_neighbors_ext_mesh, &
                                         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! transfers MPI buffers onto GPU, includes mpi-wait
      call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                       request_recv_vector_ext_mesh)

      ! assemble values on GPU to accel field
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,array_val, Mesh_pointer, &
                                          buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                          1)

    else

      !#TODO: LTS on GPU stacey contribution
      stop 'LTS on GPU w/ assemble MPI receive ilevel contribution not implemented yet'

!
!      !daniel: this avoids calling the Fortran vector send from CUDA routine
!      ! wait for asynchronous copy to finish
!
!      num_refine = num_interface_p_refine_boundary(ilevel)
!
!      allocate(reduced_buffer_send_vector_ext_mesh(3*num_refine), &
!               reduced_buffer_recv_vector_ext_mesh(3*num_refine),stat=ier)
!      if (ier /= 0) call exit_mpi(myrank,"Error allocating reduced_buffer_send_vector")
!
!      call sync_copy_reduced_from_device(Mesh_pointer,2,reduced_buffer_send_vector_ext_mesh,num_refine)
!
!      ! when testing, zero out unused portions
!      if (TEST_GPU) buffer_send_vector_ext_mesh = 0
!
!      ! fill reduced buffer into full buffer for sending
!      ! we could just calculate index+offsets in the send, but this is easier for now
!      index = 0
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!            index = index + 1
!            buffer_send_vector_ext_mesh(1,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!            index = index + 1
!            buffer_send_vector_ext_mesh(2,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!            index = index + 1
!            buffer_send_vector_ext_mesh(3,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!          enddo
!        endif
!      enddo
!      ! checks
!      if (3*num_refine /= index) call exit_mpi(myrank,"ASSERT FAIL: reduced boundary interface count incorrect")
!
!      if (TEST_GPU) then
!        ! get accel boundary
!        call copy_accelfield_from_gpu(Mesh_pointer,array_val)
!        allocate(buffer_send_vector_ext_mesh_test(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
!        buffer_send_vector_ext_mesh_test = 0
!        ! build dummy send vector
!        index = 0
!        do iinterface = 1, num_interfaces_ext_mesh
!          if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!            do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!              iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
!
!              index = index + 1
!              if (interface_p_refine_boundary(index,ilevel) /= iglob) call exit_mpi(myrank,"iglob's don't match!")
!
!              buffer_send_vector_ext_mesh_test(:,ipoin,iinterface) = array_val(:,iglob)
!            enddo
!          endif
!        enddo
!        ! compare two results
!
!        if (maxval(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test)) > 1e-10) then
!          loc = maxloc(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test))
!          print *, "diff send_vector:",maxval(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test)), &
!                ":: ",buffer_send_vector_ext_mesh(loc(1),loc(2),loc(3)),buffer_send_vector_ext_mesh_test(loc(1),loc(2),loc(3)), &
!                "ilevel=",ilevel
!        endif
!
!      endif
!
!      ! sends MPI buffers
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
!                        NDIM*num_interface_p_refine_ibool(iinterface,ilevel), &
!                        my_neighbors_ext_mesh(iinterface), &
!                        itag, &
!                        request_send_vector_ext_mesh(iinterface))
!          call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
!                        NDIM*num_interface_p_refine_ibool(iinterface,ilevel), &
!                        my_neighbors_ext_mesh(iinterface), &
!                        itag, &
!                        request_recv_vector_ext_mesh(iinterface))
!        endif
!      enddo
!
!      ! wait for communications completion (recv)
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0) then
!          call wait_req(request_recv_vector_ext_mesh(iinterface))
!        endif
!      enddo
!
!      ! fill full buffer into reduced buffer for transfer to GPU
!      ! we could just calculate index+offsets in the send, but this is easier for now
!      index = 0
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(1,ipoin,iinterface)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(2,ipoin,iinterface)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(3,ipoin,iinterface)
!          enddo
!        endif
!      enddo
!
!      ! async sends boundary to device (on memory copy stream)
!      call transfer_reduced_boundary_to_device_async_lts(Mesh_pointer,reduced_buffer_recv_vector_ext_mesh,num_refine)
!
!      ! assemble values on GPU to accel field
!      call assemble_reduced_mpi_device_lts(Mesh_pointer,ilevel,num_refine)
!
!      ! wait for sends to complete
!      do iinterface = 1, num_interfaces_ext_mesh
!        call wait_req(request_send_vector_ext_mesh(iinterface))
!      enddo
!
!      deallocate(reduced_buffer_send_vector_ext_mesh)
!      deallocate(reduced_buffer_recv_vector_ext_mesh)

    endif ! ilevel == num_p_level
  endif ! GPU

  end subroutine assemble_MPI_vector_async_recv_lts



!
!-------------------------------------------------------------------------------------------------
!

! MPI routines that only recv if *this* p-level has elements on the MPI boundary

  subroutine assemble_MPI_vector_async_w_ord_lts(NPROC,NGLOB_AB,array_val,ilevel, &
                                                 buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                                 max_nibool_interfaces_ext_mesh, &
                                                 nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                                 request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                                 my_neighbors_ext_mesh)

! waits for send/receiver to be completed and assembles contributions

! The goal of this version is to avoid different round-off errors in different processors.
! The contribution of each processor is added following the order of its rank.
! This guarantees that the sums are done in the same order on all processors.
!
! NOTE: this version assumes that the interfaces are ordered by increasing rank of the neighbor.
! That is currently done so in subroutine write_interfaces_database in decompose_mesh_SCOTCH/part_decompose_mesh_SCOTCH.f90
! A safety test could be added here.
!
! October 2012 - Surendra Somala and Jean-Paul Ampuero - Caltech Seismolab

  use constants, only: NDIM,CUSTOM_REAL,ASSEMBLE_MPI_OFF,myrank

  use specfem_par, only: GPU_MODE, Mesh_pointer, buffer_send_vector_ext_mesh

  use specfem_par_lts, only: num_interface_p_refine_ibool, interface_p_refine_ibool, num_p_level

  ! GPU
  !use specfem_par_lts, only: num_interface_p_refine_boundary, interface_p_refine_boundary

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB
  integer,intent(in) :: ilevel

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: array_val

  integer,intent(in) :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(inout) :: &
    buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh),intent(in) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),intent(in) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(inout) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh),intent(in) :: my_neighbors_ext_mesh

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: mybuffer
  integer :: ipoin,iinterface,iglob
  logical :: need_add_my_contrib
  integer :: num_buffer_points
  ! GPU
  !real(kind=CUSTOM_REAL), dimension(:), allocatable :: reduced_buffer_recv_vector_ext_mesh, reduced_buffer_send_vector_ext_mesh
  !integer :: num_refine,index,ier
  ! testing
  !logical, parameter :: TEST_GPU = .false.
  !integer :: loc(3)
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh_test

  ! here we have to assemble all the contributions between partitions using MPI

  ! debug: no mpi
  if (ASSEMBLE_MPI_OFF) return

  ! assemble only if more than one partition
  if (NPROC == 1) return

  ! receives and assembles values
  if (.not. GPU_MODE) then
    ! on CPU
    ! move interface values of array_val to local buffers
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        do ipoin = 1, num_buffer_points
          iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
          mybuffer(:,ipoin,iinterface) = array_val(:,iglob)
          ! set them to zero right away to avoid counting it more than once during assembly:
          ! buffers of higher rank get zeros on nodes shared with current buffer
          array_val(:,iglob) = 0.0_CUSTOM_REAL
        enddo
      endif
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        call wait_req(request_recv_vector_ext_mesh(iinterface))
      endif
    enddo

    ! note: assumes that the ranks in my_neighbors***(:) array are ordered in increasing values
    !       in this case, the contributions are added in the same order on both partition sides

    ! adding all contributions in order of processor rank
    need_add_my_contrib = .true.
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        if (need_add_my_contrib .and. myrank < my_neighbors_ext_mesh(iinterface)) call add_my_contrib_lts()
        do ipoin = 1, num_buffer_points
          iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
          array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
        enddo
      endif
    enddo
    if (need_add_my_contrib) call add_my_contrib_lts()

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      ! number of points on this interface
      num_buffer_points = num_interface_p_refine_ibool(iinterface,ilevel)
      if (num_buffer_points > 0) then
        call wait_req(request_send_vector_ext_mesh(iinterface))
      endif
    enddo

  else
    ! on GPU
    ! while Inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (ilevel == num_p_level) then
      ! default transfers
      ! coarse level needs full array
      call sync_copy_from_device(Mesh_pointer,2,buffer_send_vector_ext_mesh)

      ! sends MPI buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                                         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, &
                                         my_neighbors_ext_mesh, &
                                         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! transfers MPI buffers onto GPU, includes mpi-wait
      call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                       request_recv_vector_ext_mesh)

      ! assemble values on GPU to accel field
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,array_val, Mesh_pointer, &
                                          buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                          1)

    else

      !#TODO: LTS on GPU stacey contribution
      stop 'LTS on GPU w/ assemble MPI receive ilevel contribution not implemented yet'

!
!      !daniel: this avoids calling the Fortran vector send from CUDA routine
!      ! wait for asynchronous copy to finish
!
!      num_refine = num_interface_p_refine_boundary(ilevel)
!
!      allocate(reduced_buffer_send_vector_ext_mesh(3*num_refine), &
!               reduced_buffer_recv_vector_ext_mesh(3*num_refine),stat=ier)
!      if (ier /= 0) call exit_mpi(myrank,"Error allocating reduced_buffer_send_vector")
!
!      call sync_copy_reduced_from_device(Mesh_pointer,2,reduced_buffer_send_vector_ext_mesh,num_refine)
!
!      ! when testing, zero out unused portions
!      if (TEST_GPU) buffer_send_vector_ext_mesh = 0
!
!      ! fill reduced buffer into full buffer for sending
!      ! we could just calculate index+offsets in the send, but this is easier for now
!      index = 0
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!            index = index + 1
!            buffer_send_vector_ext_mesh(1,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!            index = index + 1
!            buffer_send_vector_ext_mesh(2,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!            index = index + 1
!            buffer_send_vector_ext_mesh(3,ipoin,iinterface) = reduced_buffer_send_vector_ext_mesh(index)
!          enddo
!        endif
!      enddo
!      ! checks
!      if (3*num_refine /= index) call exit_mpi(myrank,"ASSERT FAIL: reduced boundary interface count incorrect")
!
!      if (TEST_GPU) then
!        ! get accel boundary
!        call copy_accelfield_from_gpu(Mesh_pointer,array_val)
!        allocate(buffer_send_vector_ext_mesh_test(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
!        buffer_send_vector_ext_mesh_test = 0
!        ! build dummy send vector
!        index = 0
!        do iinterface = 1, num_interfaces_ext_mesh
!          if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!            do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!              iglob = interface_p_refine_ibool(ipoin,iinterface,ilevel)
!
!              index = index + 1
!              if (interface_p_refine_boundary(index,ilevel) /= iglob) call exit_mpi(myrank,"iglob's don't match!")
!
!              buffer_send_vector_ext_mesh_test(:,ipoin,iinterface) = array_val(:,iglob)
!            enddo
!          endif
!        enddo
!        ! compare two results
!
!        if (maxval(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test)) > 1e-10) then
!          loc = maxloc(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test))
!          print *, "diff send_vector:",maxval(abs(buffer_send_vector_ext_mesh - buffer_send_vector_ext_mesh_test)), &
!                ":: ",buffer_send_vector_ext_mesh(loc(1),loc(2),loc(3)),buffer_send_vector_ext_mesh_test(loc(1),loc(2),loc(3)), &
!                "ilevel=",ilevel
!        endif
!
!      endif
!
!      ! sends MPI buffers
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
!                        NDIM*num_interface_p_refine_ibool(iinterface,ilevel), &
!                        my_neighbors_ext_mesh(iinterface), &
!                        itag, &
!                        request_send_vector_ext_mesh(iinterface))
!          call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
!                        NDIM*num_interface_p_refine_ibool(iinterface,ilevel), &
!                        my_neighbors_ext_mesh(iinterface), &
!                        itag, &
!                        request_recv_vector_ext_mesh(iinterface))
!        endif
!      enddo
!
!      ! wait for communications completion (recv)
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0) then
!          call wait_req(request_recv_vector_ext_mesh(iinterface))
!        endif
!      enddo
!
!      ! fill full buffer into reduced buffer for transfer to GPU
!      ! we could just calculate index+offsets in the send, but this is easier for now
!      index = 0
!      do iinterface = 1, num_interfaces_ext_mesh
!        if (num_interface_p_refine_ibool(iinterface,ilevel) > 0 ) then
!          do ipoin = 1, num_interface_p_refine_ibool(iinterface,ilevel)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(1,ipoin,iinterface)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(2,ipoin,iinterface)
!            index = index + 1
!            reduced_buffer_recv_vector_ext_mesh(index) = buffer_recv_vector_ext_mesh(3,ipoin,iinterface)
!          enddo
!        endif
!      enddo
!
!      ! async sends boundary to device (on memory copy stream)
!      call transfer_reduced_boundary_to_device_async_lts(Mesh_pointer,reduced_buffer_recv_vector_ext_mesh,num_refine)
!
!      ! assemble values on GPU to accel field
!      call assemble_reduced_mpi_device_lts(Mesh_pointer,ilevel,num_refine)
!
!      ! wait for sends to complete
!      do iinterface = 1, num_interfaces_ext_mesh
!        call wait_req(request_send_vector_ext_mesh(iinterface))
!      enddo
!
!      deallocate(reduced_buffer_send_vector_ext_mesh)
!      deallocate(reduced_buffer_recv_vector_ext_mesh)
    endif ! ilevel == num_p_level
  endif ! GPU

  contains

    subroutine add_my_contrib_lts()

      integer :: my_iinterface,my_ipoin,my_num_buffer_points,my_iglob

      do my_iinterface = 1, num_interfaces_ext_mesh
        ! number of points on this interface
        my_num_buffer_points = num_interface_p_refine_ibool(my_iinterface,ilevel)
        if (my_num_buffer_points > 0) then
          do my_ipoin = 1, my_num_buffer_points
            my_iglob = interface_p_refine_ibool(my_ipoin,my_iinterface,ilevel)
            array_val(:,my_iglob) = array_val(:,my_iglob) + mybuffer(:,my_ipoin,my_iinterface)
          enddo
        endif
      enddo
      need_add_my_contrib = .false.

    end subroutine add_my_contrib_lts

  end subroutine assemble_MPI_vector_async_w_ord_lts
