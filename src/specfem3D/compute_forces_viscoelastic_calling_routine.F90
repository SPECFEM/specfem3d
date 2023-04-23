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

! elastic solver

  subroutine compute_forces_viscoelastic_calling()


  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par

  ! fault simulations
  use constants, only: FAULT_SYNCHRONIZE_DISPL_VELOC,FAULT_SYNCHRONIZE_ACCEL
  use fault_solver_dynamic, only: bc_dynflt_set3d_all,SIMULATION_TYPE_DYN,fault_output_synchronize_GPU,NT_RECORD_LENGTH
  use fault_solver_kinematic, only: bc_kinflt_set_all,SIMULATION_TYPE_KIN

  implicit none

  integer :: iphase
  integer :: iface,ispec,igll,i,j,k,iglob,ispec_CPML
  logical :: backward_simulation

  ! debug timing
  double precision, external :: wtime
  double precision :: t_start,tCPU
  logical, parameter :: DO_TIMING = .false.

  ! GPU
  if (GPU_MODE) then
    ! checks if for kernel simulation with both, forward & backward fields
    if (SIMULATION_TYPE == 3 &
        .and. .not. UNDO_ATTENUATION_AND_OR_PML &
        .and. .not. (ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION)) then
      ! runs with the additionally optimized GPU routine
      ! (combines forward/backward fields in main compute_kernel_acoustic)
      call compute_forces_viscoelastic_GPU_calling()
      ! all done
      return
    endif
  endif

  ! forward fields
  backward_simulation = .false.

  ! saftey check
  if (GPU_MODE .and. PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! kbai added the following two synchronizations to ensure that the displacement and velocity values
  ! at nodes on MPI interfaces stay equal on all processors that share the node.
  ! Do this only for dynamic rupture simulations
  if (FAULT_SIMULATION) then
    ! only needed for parallel simulations
    if (FAULT_SYNCHRONIZE_DISPL_VELOC .and. SIMULATION_TYPE_DYN .and. NPROC > 1) then
      ! GPU needs to transfer fields to CPU first - todo: in future, limit this to the MPI boundary buffers only
      if (GPU_MODE) then
        ! transfers displacement & velocity to the CPU
        call transfer_displ_from_device(NDIM*NGLOB_AB, displ, Mesh_pointer)
        call transfer_veloc_from_device(NDIM*NGLOB_AB, veloc, Mesh_pointer)
      endif
      ! displacement
      call synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,displ, &
                                               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                               my_neighbors_ext_mesh)
      ! velocity
      call synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,veloc, &
                                               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                               my_neighbors_ext_mesh)
      ! GPU: copies "corrected" wavefields from CPU to GPU
      if (GPU_MODE) then
        ! transfers displacement & velocity back to the GPU
        call transfer_displ_to_device(NDIM*NGLOB_AB,displ, Mesh_pointer)
        call transfer_veloc_to_device(NDIM*NGLOB_AB,veloc, Mesh_pointer)
      endif
    endif
  endif

! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! debug timing
    if (DO_TIMING .and. myrank == 0 .and. iphase == 2) then
      t_start = wtime()
    endif

    ! elastic term
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_forces_viscoelastic(iphase,deltat, &
                                       displ,veloc,accel, &
                                       alphaval,betaval,gammaval, &
                                       R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                       R_trace_lddrk, &
                                       R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                       epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                       epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                                       backward_simulation)
    else
      ! on GPU
      call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, deltat, &
                                            nspec_outer_elastic, &
                                            nspec_inner_elastic, &
                                            COMPUTE_AND_STORE_STRAIN,ATTENUATION,1) ! 1 == forward
    endif

    ! debug timing
    if (DO_TIMING .and. myrank == 0 .and. iphase == 2) then
      tCPU = wtime() - t_start
      print *,'timing: compute_forces_viscoelastic elapsed time ',tCPU,'s'
    endif

    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (GPU_MODE .and. iphase == 2) then
      !daniel: todo - this avoids calling the Fortran vector send from CUDA routine
      ! wait for asynchronous copy to finish
      call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)

      ! sends MPI buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                                         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, &
                                         my_neighbors_ext_mesh, &
                                         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! transfers MPI buffers onto GPU
      call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                       request_recv_vector_ext_mesh)
    endif ! inner elements

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        if (.not. GPU_MODE) then
          ! on CPU
          call compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel, &
                                                   ibool,iphase, &
                                                   abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                   abs_boundary_ijk,abs_boundary_ispec, &
                                                   num_abs_boundary_faces, &
                                                   veloc,rho_vp,rho_vs, &
                                                   ispec_is_elastic, &
                                                   it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
        else
          ! on GPU
          call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                               NSTEP,it, &
                                               b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                                               Mesh_pointer,1) ! 1 == forward
        endif
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          if (.not. GPU_MODE) then
            ! on CPU
            if (SIMULATION_TYPE == 1) then
              ! forward definition: pressure = - potential_dot_dot
              call compute_coupling_viscoelastic_ac(NSPEC_AB,NGLOB_AB, &
                                                    ibool,accel,potential_dot_dot_acoustic, &
                                                    num_coupling_ac_el_faces, &
                                                    coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                                    coupling_ac_el_normal, &
                                                    coupling_ac_el_jacobian2Dw, &
                                                    iphase, &
                                                    PML_CONDITIONS, &
                                                    SIMULATION_TYPE,backward_simulation, &
                                                    potential_acoustic,potential_dot_acoustic)


            else
              ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
              ! adjoint definition: pressure^\dagger = potential^\dagger
              call compute_coupling_viscoelastic_ac(NSPEC_AB,NGLOB_AB, &
                                                    ibool,accel,potential_acoustic, &
                                                    num_coupling_ac_el_faces, &
                                                    coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                                    coupling_ac_el_normal, &
                                                    coupling_ac_el_jacobian2Dw, &
                                                    iphase, &
                                                    PML_CONDITIONS, &
                                                    SIMULATION_TYPE,backward_simulation, &
                                                    potential_acoustic,potential_dot_acoustic)

            endif
          else
            ! on GPU
            call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,1) ! 1 == forward
          endif
        endif ! num_coupling_ac_el_faces
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          ! on CPU
          if (SIMULATION_TYPE == 1) then
            call compute_coupling_viscoelastic_po(iphase)
          else
            stop 'poroelastic/elastic coupling for adjoint/kernel simulations not implemented yet'
          endif
        else
          ! on GPU
          if (SIMULATION_TYPE == 1) then
            ! note:
            ! these routines are not implemented as CUDA kernels, we just transfer the fields
            ! from the GPU to the CPU and vice versa

            ! transfers displacement & acceleration to the CPU
            call transfer_displ_from_device(NDIM*NGLOB_AB,displ, Mesh_pointer)
            call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

            call compute_coupling_viscoelastic_po(iphase)

            ! transfers acceleration back to GPU
            call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
          else
            stop 'for kernel simulations, poroelastic/elastic coupling on GPUs not implemented yet'
          endif
        endif
      endif


      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      ! adds source term (single-force/moment-tensor solution)
      if (.not. GPU_MODE) then
        ! on CPU
        call compute_add_sources_viscoelastic()
      else
        ! on GPU
        call compute_add_sources_viscoelastic_GPU()
      endif
    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      if (.not. GPU_MODE) then
        ! on CPU
        ! sends accel values to corresponding MPI interface neighbors
        call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,accel, &
                                            buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            my_neighbors_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh)
      else
        ! on GPU
        ! transfers boundary region to host asynchronously. The
        ! MPI-send is done from within compute_forces_viscoelastic_cuda,
        ! once the inner element kernels are launched, and the
        ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655
        call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)
      endif
    else
      ! waits for send/receive requests to be completed and assembles values
      if (.not. GPU_MODE) then
        ! on CPU
        ! receives MPI buffers
        call assemble_MPI_vector_async_recv(NPROC,NGLOB_AB,accel, &
                                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            my_neighbors_ext_mesh)
      else
        ! on GPU
        ! waits for send/receive requests to be completed and assembles values
        call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,accel, Mesh_pointer, &
                                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            1)
      endif
    endif
  enddo

  ! Fault boundary term B*tau is added to the assembled forces
  ! which at this point are stored in the array 'accel'
  if (FAULT_SIMULATION) then
    ! makes sure all processes have same accel values on MPI boundary points
    if (FAULT_SYNCHRONIZE_ACCEL .and. NPROC > 1) then
      ! note: the CPU version uses an "ordered" assembly stage, where the MPI buffer contributions are added
      !       in the same order on all processes. this should assure that accel is the same on all MPI boundary points.
      !
      !       the GPU version doesn't ensure this ordered assembly. thus, one might have to enforce the same values
      !       on the boundary points with this following setup.
      !
      ! using the fault_examples/tpv5 as a test, this is not needed for accuracy.
      ! however, this should to be tested further with more cases...
      !
      ! GPU needs to transfer fields to CPU first - todo: in future, limit this to the MPI boundary buffers only
      if (GPU_MODE) then
        ! transfers to the CPU
        call transfer_accel_from_device(NDIM*NGLOB_AB, accel, Mesh_pointer)
      endif
      ! acceleration
      call synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,accel, &
                                               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                               my_neighbors_ext_mesh)
      ! GPU: copies "corrected" wavefields from CPU to GPU
      if (GPU_MODE) then
        ! transfers back to the GPU
        call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
      endif
    endif

    ! computes fault force contribution
    if (.not. GPU_MODE) then
      ! on CPU
      ! dynamic rupture
      if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(accel,veloc,displ)
      ! kinematic rupture
      if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(accel,veloc,displ)
    else
      ! on GPU
      ! GPU fault solver
      call fault_solver_gpu(Mesh_pointer,Fault_pointer,deltat,it)
      ! output results every 500 steps
      if ((mod(it,NT_RECORD_LENGTH) == 0 .or. it == it_end) .and. it /= 0) call fault_output_synchronize_GPU(it)
    endif
  endif

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  if (.not. GPU_MODE) then
    ! on CPU
    do iglob = 1,NGLOB_AB
      accel(1,iglob) = accel(1,iglob)*rmassx(iglob)
      accel(2,iglob) = accel(2,iglob)*rmassy(iglob)
      accel(3,iglob) = accel(3,iglob)*rmassz(iglob)
    enddo
  else
    ! on GPU
    call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,APPROXIMATE_OCEAN_LOAD,1) ! 1 == forward
  endif

  ! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_coupling_ocean(NSPEC_AB,NGLOB_AB, &
                                  ibool,rmassx,rmassy,rmassz, &
                                  rmass_ocean_load,accel, &
                                  free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                  num_free_surface_faces)
    else
      ! on GPU
      call compute_coupling_ocean_cuda(Mesh_pointer,1) ! 1 == forward
    endif
  endif

  ! PML
  ! impose Dirichlet conditions on the outer edges of the C-PML layers
  if (PML_CONDITIONS) then
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
!!! It is better to move this into do iphase=1,2 loop
      if (ispec_is_elastic(ispec) .and. is_CPML(ispec)) then
        ! reference GLL points on boundary face
        ispec_CPML = spec_to_CPML(ispec)
        do igll = 1,NGLLSQUARE
          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          iglob = ibool(i,j,k,ispec)

          accel(:,iglob) = 0._CUSTOM_REAL
          veloc(:,iglob) = 0._CUSTOM_REAL
          displ(:,iglob) = 0._CUSTOM_REAL
          PML_displ_old(:,i,j,k,ispec_CPML) = 0._CUSTOM_REAL
          PML_displ_new(:,i,j,k,ispec_CPML) = 0._CUSTOM_REAL

        enddo
      endif ! ispec_is_elastic
!!!        endif
    enddo
  endif

! updates velocities
! Newmark finite-difference time scheme with elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t))
!
! where
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!   chi_dot_dot is acoustic (fluid) potential ( dotted twice with respect to time)
!
! corrector:
!   updates the velocity term which requires a(t+delta)
  if (.not. GPU_MODE) then
    ! on CPU
    if (USE_LDDRK) then
      call update_veloc_elastic_lddrk()
    else
      veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
    endif
  else
    ! on GPU
    ! only call in case of ocean load, otherwise already done in kernel 3 a
    if (APPROXIMATE_OCEAN_LOAD) then
      call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2,1) ! 1 == forward
    endif
  endif

  ! PML
  if (PML_CONDITIONS) then
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      if (nglob_interface_PML_elastic > 0) then
        call save_field_on_pml_interface(displ,veloc,accel,nglob_interface_PML_elastic, &
                                         b_PML_field,b_reclen_PML_field)
      endif
    endif
  endif

  end subroutine compute_forces_viscoelastic_calling

!
!=====================================================================
!

! elastic solver for backward/reconstructed wavefields

  subroutine compute_forces_viscoelastic_backward_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  integer :: iphase
  logical :: backward_simulation

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling compute_forces_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! GPU
  if (GPU_MODE) then
    ! checks if for kernel simulation with both, forward & backward fields
    if (SIMULATION_TYPE == 3 &
        .and. .not. UNDO_ATTENUATION_AND_OR_PML &
        .and. .not. (ELASTIC_SIMULATION .and. ACOUSTIC_SIMULATION)) then
      ! runs with the additionally optimized GPU routine
      ! (combines forward/backward fields in main compute_kernel_acoustic)
      ! all done in compute_forces_acoustic_GPU_calling()
      return
    endif
  endif

  ! backward fields
  backward_simulation = .true.

  ! saftey check
  if (GPU_MODE .and. PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! check
  if (FAULT_SIMULATION) then
    stop 'FAULT_SIMULATION (of TYPE_DYN and TYPE_KIN) for backward simulation routine not implemented yet'
  endif

  ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! elastic term
    if (.not. GPU_MODE) then
      ! on CPU
      ! adjoint simulations: backward/reconstructed wavefield
      call compute_forces_viscoelastic(iphase,b_deltat, &
                                       b_displ,b_veloc,b_accel, &
                                       b_alphaval,b_betaval,b_gammaval, &
                                       b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                                       b_R_trace_lddrk, &
                                       b_R_xx_lddrk,b_R_yy_lddrk,b_R_xy_lddrk,b_R_xz_lddrk,b_R_yz_lddrk, &
                                       b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                                       b_epsilondev_xz,b_epsilondev_yz,b_epsilon_trace_over_3, &
                                       backward_simulation)
    else
      ! on GPU
      call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, b_deltat, &
                                            nspec_outer_elastic, &
                                            nspec_inner_elastic, &
                                            COMPUTE_AND_STORE_STRAIN,ATTENUATION,3) ! 3 == backward
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        if (.not. GPU_MODE) then
          ! on CPU
          if (UNDO_ATTENUATION_AND_OR_PML) then
            call compute_stacey_viscoelastic_backward_undoatt(NSPEC_AB,NGLOB_AB,b_accel,b_veloc, &
                                                              ibool,iphase, &
                                                              abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                              abs_boundary_ijk,abs_boundary_ispec, &
                                                              num_abs_boundary_faces, &
                                                              rho_vp,rho_vs,ispec_is_elastic)
          else
            call compute_stacey_viscoelastic_backward(NSPEC_AB, &
                                                      ibool,iphase, &
                                                      abs_boundary_ijk,abs_boundary_ispec, &
                                                      num_abs_boundary_faces, &
                                                      ispec_is_elastic,SIMULATION_TYPE, &
                                                      NSTEP,it,NGLOB_ADJOINT,b_accel, &
                                                      b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
          endif
        else
          ! on GPU
          call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                               NSTEP,it, &
                                               b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                                               Mesh_pointer,3) ! 3 == backward
        endif
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          if (.not. GPU_MODE) then
            ! on CPU
            ! backward simulations
            call compute_coupling_viscoelastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                                                  ibool,b_accel,b_potential_dot_dot_acoustic, &
                                                  num_coupling_ac_el_faces, &
                                                  coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                                  coupling_ac_el_normal, &
                                                  coupling_ac_el_jacobian2Dw, &
                                                  iphase, &
                                                  PML_CONDITIONS, &
                                                  SIMULATION_TYPE,backward_simulation, &
                                                  potential_acoustic,potential_dot_acoustic)
          else
            ! on GPU
            call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,3) ! 3 == backward
          endif
        endif ! num_coupling_ac_el_faces
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        stop 'poroelastic-elastic coupling error, for backward routine not implemented yet'
      endif

      ! adds source term (single-force/moment-tensor solution)
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      if (.not. GPU_MODE) then
        ! on CPU
        call compute_add_sources_viscoelastic_backward()
      else
        ! on GPU
        call compute_add_sources_viscoelastic_backward_GPU()
      endif

    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends accel values to corresponding MPI interface neighbors
      if (.not. GPU_MODE) then
        ! on CPU
        ! adjoint simulations
        call assemble_MPI_vector_async_send(NPROC,NGLOB_ADJOINT,b_accel, &
                                            b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            my_neighbors_ext_mesh, &
                                            b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      else
        ! on GPU
        call transfer_boun_accel_from_device(Mesh_pointer, &
                                             b_buffer_send_vector_ext_mesh, &
                                             3) ! 3 == adjoint b_accel

        call assemble_MPI_vector_send_cuda(NPROC, &
                                           b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif
    else
      ! waits for send/receive requests to be completed and assembles values
      if (.not. GPU_MODE) then
        ! on CPU
        ! adjoint simulations
        call assemble_MPI_vector_async_recv(NPROC,NGLOB_ADJOINT,b_accel, &
                                            b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh, &
                                            my_neighbors_ext_mesh)
      else
        ! on GPU
        call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,b_accel, Mesh_pointer, &
                                            b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh, &
                                            3)
      endif
    endif
  enddo

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  if (.not. GPU_MODE) then
    ! on CPU
    ! adjoint simulations
    b_accel(1,:) = b_accel(1,:)*rmassx(:)
    b_accel(2,:) = b_accel(2,:)*rmassy(:)
    b_accel(3,:) = b_accel(3,:)*rmassz(:)
  else
    ! on GPU
    call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,APPROXIMATE_OCEAN_LOAD,3) ! 3 == backward
  endif

  ! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_coupling_ocean_backward(NSPEC_AB,NGLOB_AB, &
                                           ibool,rmassx,rmassy,rmassz, &
                                           rmass_ocean_load, &
                                           free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                           num_free_surface_faces, &
                                           SIMULATION_TYPE, &
                                           NGLOB_ADJOINT,b_accel)
    else
      ! on GPU
      call compute_coupling_ocean_cuda(Mesh_pointer,3) ! 3 == backward
    endif
  endif


! updates velocities
! Newmark finite-difference time scheme with elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t))
!
! where
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!   chi_dot_dot is acoustic (fluid) potential ( dotted twice with respect to time)
!
! corrector:
!   updates the velocity term which requires a(t+delta)
  if (.not. GPU_MODE) then
    ! on CPU
    ! adjoint simulations
    b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)
  else
    ! on GPU
    ! only call in case of ocean load, otherwise already done in kernel 3 a
    if (APPROXIMATE_OCEAN_LOAD) then
      call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2,3) ! 3 == backward
    endif
  endif

  end subroutine compute_forces_viscoelastic_backward_calling

!
!-------------------------------------------------------------------------------------------------
!

! elastic solver

  subroutine compute_forces_viscoelastic_GPU_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par

  use fault_solver_dynamic, only: bc_dynflt_set3d_all,fault_output_synchronize_GPU,NT_RECORD_LENGTH
  use fault_solver_kinematic, only: bc_kinflt_set_all

  implicit none

  integer:: iphase

  ! check
  if (PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! elastic term
    ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
    call compute_forces_viscoelastic_cuda(Mesh_pointer, iphase, deltat, &
                                          nspec_outer_elastic, &
                                          nspec_inner_elastic, &
                                          COMPUTE_AND_STORE_STRAIN,ATTENUATION,0) ! 0 == both

    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (iphase == 2) then
      !daniel: todo - this avoids calling the Fortran vector send from CUDA routine
      ! wait for asynchronous copy to finish
      call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)

      ! sends MPI buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                                         buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, &
                                         my_neighbors_ext_mesh, &
                                         request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

      ! transfers MPI buffers onto GPU
      call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
                                       num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                       request_recv_vector_ext_mesh)
    endif ! inner elements

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                             NSTEP,it, &
                                             b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                                             Mesh_pointer,0)
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,1) ! 1 == forward
          if (SIMULATION_TYPE == 3) &
            call compute_coupling_el_ac_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,3) ! 3 == backward
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        ! note:
        ! these routines are not implemented as CUDA kernels, we just transfer the fields
        ! from the GPU to the CPU and vice versa

        ! transfers displacement & acceleration to the CPU
        call transfer_displ_from_device(NDIM*NGLOB_AB,displ, Mesh_pointer)
        call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

        call compute_coupling_viscoelastic_po(iphase)

        ! transfers acceleration back to GPU
        call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
      endif

      ! adds source term (single-force/moment-tensor solution)
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      call compute_add_sources_viscoelastic_GPU()

    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends accel values to corresponding MPI interface neighbors

      ! transfers boundary region to host asynchronously. The
      ! MPI-send is done from within compute_forces_viscoelastic_cuda,
      ! once the inner element kernels are launched, and the
      ! memcpy has finished. see compute_forces_viscoelastic_cuda: ~ line 1655
      call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call transfer_boun_accel_from_device(Mesh_pointer, &
                                             b_buffer_send_vector_ext_mesh, &
                                             3) ! 3 == adjoint b_accel

        call assemble_MPI_vector_send_cuda(NPROC, &
                                           b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
      endif !adjoint

    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,accel, Mesh_pointer, &
                                          buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                          1)
      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call assemble_MPI_vector_write_cuda(NPROC,NGLOB_AB,b_accel, Mesh_pointer, &
                                            b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh, &
                                            3)
      endif !adjoint
    endif

  enddo

  ! Fault boundary term B*tau is added to the assembled forces
  ! which at this point are stored in the array 'accel'
  if (FAULT_SIMULATION) then
    ! GPU fault solver
    call fault_solver_gpu(Mesh_pointer,Fault_pointer,deltat,it)
    ! output results every 500 steps
    if ((mod(it,NT_RECORD_LENGTH) == 0 .or. it == it_end) .and. it /= 0) call fault_output_synchronize_GPU(it)
  endif

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,APPROXIMATE_OCEAN_LOAD,1) ! 1 == forward
  if (SIMULATION_TYPE == 3) call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,APPROXIMATE_OCEAN_LOAD,3) ! 3 == backward

  ! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    call compute_coupling_ocean_cuda(Mesh_pointer,1) ! 1 == forward
    if (SIMULATION_TYPE == 3) call compute_coupling_ocean_cuda(Mesh_pointer,3) ! 3 == backward

    ! updates velocities
    ! Newmark finite-difference time scheme with elastic domains:
    ! (see e.g. Hughes, 1987; Chaljub et al., 2003)
    !
    ! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
    ! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
    ! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t))
    !
    ! where
    !   u, v, a are displacement,velocity & acceleration
    !   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
    !   f denotes a source term (acoustic/elastic)
    !   chi_dot_dot is acoustic (fluid) potential ( dotted twice with respect to time)
    !
    ! corrector:
    ! updates the velocity term which requires a(t+delta)
    ! GPU_MODE: this is handled in 'kernel_3' at the same time as accel*rmass
    call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2,1) ! 1 == forward
    if (SIMULATION_TYPE == 3) call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2,3) ! 3 == backward
  endif

  end subroutine compute_forces_viscoelastic_GPU_calling

