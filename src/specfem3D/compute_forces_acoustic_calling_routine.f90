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

! acoustic solver

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(Chi) / rho
! Velocity is then:
!     v = grad(Chi_dot) / rho
! (Chi_dot being the time derivative of Chi)
! and pressure is:
!     p = - Chi_dot_dot
! (Chi_dot_dot being the time second derivative of Chi).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore Chi_dot_dot
! is continuous, therefore Chi is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u = grad(Chi) would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.

  subroutine compute_forces_acoustic_forward_calling()

! acoustic domains for forward or adjoint simulations (SIMULATION_TYPE == 1 or 2 )

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  use pml_par, only: is_CPML,spec_to_CPML,nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  ! local parameters
  integer:: iface,ispec,iglob,igll,i,j,k,ispec_CPML
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase
  logical :: backward_simulation

  ! debug timing
  double precision, external :: wtime
  double precision :: t_start,tCPU
  logical, parameter :: DO_TIMING = .false.

  ! GPU
  if (GPU_MODE) then
    ! checks if for kernel simulation with both, forward & backward fields
    if (SIMULATION_TYPE == 3 .and. .not. UNDO_ATTENUATION_AND_OR_PML) then
      ! runs with the additionally optimized GPU routine
      ! (combines forward/backward fields in main compute_kernel_acoustic)
      call compute_forces_acoustic_GPU_calling()
      ! all done
      return
    endif
  endif

  ! for forward fields only
  backward_simulation = .false.

  ! saftey check
  if (GPU_MODE .and. PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                     backward_simulation)
  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NGLOB_AB_LDDRK,potential_acoustic_lddrk,potential_dot_acoustic_lddrk)
  endif

  ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

    !debug timing
    if (DO_TIMING .and. myrank == 0 .and. iphase == 2) then
      t_start = wtime()
    endif

    ! acoustic pressure term
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_forces_acoustic(iphase,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                   backward_simulation)
    else
      ! on GPU
      call compute_forces_acoustic_cuda(Mesh_pointer, iphase, nspec_outer_acoustic, nspec_inner_acoustic, 1) ! 1 == forward
    endif

    ! debug timing
    if (DO_TIMING .and. myrank == 0 .and. iphase == 2) then
      tCPU = wtime() - t_start
      print *,'timing: compute_forces_acoustic elapsed time ',tCPU,'s'
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        if (.not. GPU_MODE) then
          ! on CPU
          call compute_stacey_acoustic_forward(NSPEC_AB,NGLOB_AB, &
                                               potential_dot_dot_acoustic,potential_dot_acoustic, &
                                               ibool,iphase, &
                                               abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                               num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                                               it,b_reclen_potential, &
                                               b_absorb_potential,b_num_abs_boundary_faces)
        else
          ! on GPU
          call compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                                           SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                           b_reclen_potential,b_absorb_potential, &
                                           b_num_abs_boundary_faces,Mesh_pointer,1) ! 1 == forward
        endif
      endif

      ! elastic coupling
      if (ELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          ! on CPU
          if (SIMULATION_TYPE == 1) then
            ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                                              ibool,displ,potential_dot_dot_acoustic, &
                                              num_coupling_ac_el_faces, &
                                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                              coupling_ac_el_normal, &
                                              coupling_ac_el_jacobian2Dw, &
                                              iphase, &
                                              PML_CONDITIONS,SIMULATION_TYPE,backward_simulation)
          else
            ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
            ! adjoint definition: \partial_t^2 \bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                                              ibool,accel,potential_dot_dot_acoustic, &
                                              num_coupling_ac_el_faces, &
                                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                              coupling_ac_el_normal, &
                                              coupling_ac_el_jacobian2Dw, &
                                              iphase, &
                                              PML_CONDITIONS,SIMULATION_TYPE,backward_simulation)
          endif
        else
          ! on GPU
          call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,1) ! 1 == forward
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION ) then
        if (.not. GPU_MODE) then
          ! on CPU
          if (SIMULATION_TYPE == 1) then
            call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                                              ibool,displs_poroelastic,displw_poroelastic, &
                                              potential_dot_dot_acoustic, &
                                              num_coupling_ac_po_faces, &
                                              coupling_ac_po_ispec,coupling_ac_po_ijk, &
                                              coupling_ac_po_normal, &
                                              coupling_ac_po_jacobian2Dw, &
                                              iphase)
          else
            stop 'poroelastic/acoustic coupling for adjoint/kernel simulations not implemented yet'
          endif
        else
          ! on GPU
          stop 'poroelastic/acoustic coupling on GPUs not implemented yet'
        endif
      endif

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      if (.not. GPU_MODE) then
        ! on CPU
        call compute_add_sources_acoustic()
      else
        ! on GPU
        call compute_add_sources_acoustic_GPU()
      endif
    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      if (.not. GPU_MODE) then
        ! on CPU
        ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
        call assemble_MPI_scalar_async_send(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                                            buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            my_neighbors_ext_mesh, &
                                            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      else
        ! on GPU
        ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           buffer_send_scalar_ext_mesh, &
                                           1) ! -- 1 == fwd accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                                           buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      endif
    else
      if (.not. GPU_MODE) then
        ! on CPU
        ! waits for send/receive requests to be completed and assembles values
        call assemble_MPI_scalar_async_recv(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                                            buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      else
        ! on GPU
        ! waits for send/receive requests to be completed and assembles values
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                                            Mesh_pointer, &
                                            buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                           ! nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                                            1)
      endif
    endif
  enddo

  ! divides pressure with mass matrix
  if (.not. GPU_MODE) then
    ! on CPU
    potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:)

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

!! DK DK this paragraph seems to be from Zhinan or from ChangHua:
! However, enforcing explicitly potential_dot_dot_acoustic, potential_dot_acoustic, potential_acoustic
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation.

! impose Dirichlet conditions for the potential (i.e. Neumann for displacement) on the outer edges of the C-PML layers
    if (PML_CONDITIONS .and. SET_NEUMANN_RATHER_THAN_DIRICHLET_FOR_FLUID_PMLs) then
      do iface = 1,num_abs_boundary_faces
        ispec = abs_boundary_ispec(iface)
  !!! It is better to move this into do iphase=1,2 loop
          if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
            ! reference GLL points on boundary face
            ispec_CPML = spec_to_CPML(ispec)
            do igll = 1,NGLLSQUARE
              ! gets local indices for GLL point
              i = abs_boundary_ijk(1,igll,iface)
              j = abs_boundary_ijk(2,igll,iface)
              k = abs_boundary_ijk(3,igll,iface)

              iglob=ibool(i,j,k,ispec)

              potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
              potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
              potential_acoustic(iglob) = 0._CUSTOM_REAL
              if (ELASTIC_SIMULATION) then
                PML_potential_acoustic_old(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
                PML_potential_acoustic_new(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
              endif
            enddo
          endif ! ispec_is_acoustic
  !!!   endif
      enddo
    endif

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the chi_dot term which requires chi_dot_dot(t+delta)
    ! corrector
    if (USE_LDDRK) then
      call update_potential_dot_acoustic_lddrk()
    else
      potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2*potential_dot_dot_acoustic(:)
    endif
  else
    ! on GPU
    call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2,1) ! 1 == forward
  endif

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                     backward_simulation)
  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NGLOB_AB_LDDRK,potential_acoustic_lddrk,potential_dot_acoustic_lddrk)
  endif

  ! coupling - not used yet
  !if (SIMULATION_TYPE /= 1) then
  !  potential_acoustic_adj_coupling(:) = potential_acoustic(:) &
  !                                      + deltat * potential_dot_acoustic(:) &
  !                                      + deltatsqover2 * potential_dot_dot_acoustic(:)
  !endif

  ! PML
  if (PML_CONDITIONS) then
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      if (nglob_interface_PML_acoustic > 0) then
        call save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)
      endif
    endif
  endif

  end subroutine compute_forces_acoustic_forward_calling

!
!=====================================================================

! acoustic solver

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(Chi) / rho
! Velocity is then:
!     v = grad(Chi_dot) / rho
! (Chi_dot being the time derivative of Chi)
! and pressure is:
!     p = - Chi_dot_dot
! (Chi_dot_dot being the time second derivative of Chi).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore Chi_dot_dot
! is continuous, therefore Chi is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u = grad(Chi) would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.


  subroutine compute_forces_acoustic_backward_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer :: iphase
  logical :: backward_simulation

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling compute_forces_acoustic_backward() with wrong SIMULATION_TYPE')

  ! GPU
  if (GPU_MODE) then
    ! checks if for kernel simulation with both, forward & backward fields
    if (SIMULATION_TYPE == 3 .and. .not. UNDO_ATTENUATION_AND_OR_PML) then
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

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NGLOB_ADJOINT,b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                                     backward_simulation)
  if (USE_LDDRK) then
    stop 'LDDRK for backward fields not implemented yet'
    !call acoustic_enforce_free_surface_lddrk(NGLOB_AB_LDDRK,b_potential_acoustic_lddrk,b_potential_dot_acoustic_lddrk)
  endif

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! acoustic pressure term
    if (.not. GPU_MODE) then
      ! on CPU
    call compute_forces_acoustic(iphase,b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                                 backward_simulation)
    else
      ! on GPU
      call compute_forces_acoustic_cuda(Mesh_pointer, iphase, nspec_outer_acoustic, nspec_inner_acoustic, 3) ! 3 == backward
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        if (.not. GPU_MODE) then
          if (UNDO_ATTENUATION_AND_OR_PML) then
            call compute_stacey_acoustic_backward_undoatt(NSPEC_AB,NGLOB_AB, &
                                                          b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                          ibool,iphase, &
                                                          abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                                          num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic)
          else
            call compute_stacey_acoustic_backward(NSPEC_AB, &
                                                  ibool,iphase, &
                                                  abs_boundary_ijk,abs_boundary_ispec, &
                                                  num_abs_boundary_faces,ispec_is_acoustic, &
                                                  SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                                                  b_potential_dot_dot_acoustic,b_reclen_potential, &
                                                  b_absorb_potential,b_num_abs_boundary_faces)
          endif
        else
          ! on GPU
          call compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                                           SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                           b_reclen_potential,b_absorb_potential, &
                                           b_num_abs_boundary_faces,Mesh_pointer,3) ! 3 == backward
        endif
      endif

      ! elastic coupling
      if (ELASTIC_SIMULATION) then
        if (.not. GPU_MODE) then
          ! on CPU
          ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
          call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                                            ibool,b_displ,b_potential_dot_dot_acoustic, &
                                            num_coupling_ac_el_faces, &
                                            coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                            coupling_ac_el_normal, &
                                            coupling_ac_el_jacobian2Dw, &
                                            iphase, &
                                            PML_CONDITIONS,SIMULATION_TYPE,backward_simulation)
        else
          ! on GPU
          call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,3) ! 3 == backward
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION ) then
        if (num_coupling_ac_po_faces > 0) then
            stop 'coupling acoustic-poroelastic domains in backward routine not implemented yet...'
        endif
      endif

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      if (.not. GPU_MODE) then
        ! on CPU
        call compute_add_sources_acoustic_backward()
      else
        ! on GPU
        call compute_add_sources_acoustic_backward_GPU()
      endif

    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      if (.not. GPU_MODE) then
        ! on CPU
        ! sends b_potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
        ! adjoint simulations
        call assemble_MPI_scalar_async_send(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                                            b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            my_neighbors_ext_mesh, &
                                            b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
      else
        ! on GPU
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_buffer_send_scalar_ext_mesh, &
                                           3) ! -- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                                           b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
      endif
    else
      if (.not. GPU_MODE) then
        ! on CPU
        ! adjoint simulations
        call assemble_MPI_scalar_async_recv(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                                            b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
      else
        ! on GPU
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_potential_dot_dot_acoustic, &
                                            Mesh_pointer, &
                                            b_buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                           ! nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh, &
                                            3) ! 3 == backward
      endif
    endif
  enddo

  ! divides pressure with mass matrix
  if (.not. GPU_MODE) then
    ! on CPU
    ! adjoint simulations
    b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_acoustic(:)

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the chi_dot term which requires chi_dot_dot(t+delta)

    ! corrector
    if (USE_LDDRK) then
      stop 'LDDRK scheme for backward propagation not implemented yet'
    else
      ! adjoint simulations
      b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2*b_potential_dot_dot_acoustic(:)
    endif
  else
    ! on GPU
    call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2,3) ! 3 == backward
  endif

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NGLOB_ADJOINT,b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                                     backward_simulation)
  if (USE_LDDRK) then
    stop 'LDDRK for backward fields not implemented yet'
    !call acoustic_enforce_free_surface_lddrk(NGLOB_AB_LDDRK,b_potential_acoustic_lddrk,b_potential_dot_acoustic_lddrk)
  endif

  end subroutine compute_forces_acoustic_backward_calling

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_forces_acoustic_GPU_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: iphase

  ! safety check
  if (.not. GPU_MODE) return
  ! check
  if (PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,1)
  if (SIMULATION_TYPE == 3) call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,3)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase,nspec_outer_acoustic, nspec_inner_acoustic,0) ! 0 == both

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
         call compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                                          SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                          b_reclen_potential,b_absorb_potential, &
                                          b_num_abs_boundary_faces,Mesh_pointer,0) ! 0 == both
      endif

      ! elastic coupling
      if (ELASTIC_SIMULATION) then
        call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,1) ! 1 == forward
        if (SIMULATION_TYPE == 3) call compute_coupling_ac_el_cuda(Mesh_pointer,iphase,num_coupling_ac_el_faces,3)
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION ) then
        ! on CPU
        !if (SIMULATION_TYPE == 1) then
        !  call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
        !                                    ibool,displs_poroelastic,displw_poroelastic, &
        !                                    potential_dot_dot_acoustic, &
        !                                    num_coupling_ac_po_faces, &
        !                                    coupling_ac_po_ispec,coupling_ac_po_ijk, &
        !                                    coupling_ac_po_normal, &
        !                                    coupling_ac_po_jacobian2Dw, &
        !                                    iphase)
        !else
        !  stop 'poroelastic/acoustic coupling for adjoint/kernel simulations not implemented yet'
        !endif
        ! on GPU
        stop 'poroelastic/acoustic coupling on GPUs not implemented yet'
      endif

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      call compute_add_sources_acoustic_GPU()          ! forward/adjoint sources
      call compute_add_sources_acoustic_backward_GPU() ! backward sources
    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call transfer_boun_pot_from_device(Mesh_pointer, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! -- 1 == fwd accel

      call assemble_MPI_scalar_send_cuda(NPROC, &
                                         buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                         num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                         nibool_interfaces_ext_mesh, &
                                         my_neighbors_ext_mesh, &
                                         request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_buffer_send_scalar_ext_mesh, &
                                           3) ! -- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                                           b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

      endif

    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                                          Mesh_pointer, &
                                          buffer_recv_scalar_ext_mesh, &
                                          num_interfaces_ext_mesh, &
                                          max_nibool_interfaces_ext_mesh, &
                                         ! nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                          request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                                          1)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_potential_dot_dot_acoustic, &
                                            Mesh_pointer, &
                                            b_buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                           ! nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh, &
                                            3)
      endif
    endif
  enddo

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T
! CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic
! u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with
!   respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
! updates the chi_dot term which requires chi_dot_dot(t+delta)

  call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2,0) ! 0 == both

! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,1)
  if (SIMULATION_TYPE == 3) call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,3)

  end subroutine compute_forces_acoustic_GPU_calling

!
!-------------------------------------------------------------------------------------------------
!

  subroutine acoustic_enforce_free_surface(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                           backward_simulation)

  use constants
  use specfem_par, only: Mesh_pointer,ibool,free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
    STACEY_INSTEAD_OF_FREE_SURFACE,BOTTOM_FREE_SURFACE,GPU_MODE
  use specfem_par_acoustic, only: ispec_is_acoustic

  implicit none

  integer,intent(in) :: NGLOB_AB

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  logical, intent(in) :: backward_simulation

  ! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE .and. .not. BOTTOM_FREE_SURFACE) return

  ! enforce potentials to be zero at surface
  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1, num_free_surface_faces
      ispec = free_surface_ispec(iface)
      if (ispec_is_acoustic(ispec)) then
        do igll = 1, NGLLSQUARE
          i = free_surface_ijk(1,igll,iface)
          j = free_surface_ijk(2,igll,iface)
          k = free_surface_ijk(3,igll,iface)
          iglob = ibool(i,j,k,ispec)

          ! sets potentials to zero
          potential_acoustic(iglob)         = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob)     = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif
    enddo
  else
    ! on GPU
    if (backward_simulation) then
      ! backward arrays
      call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,3) ! 3 == backward
    else
      ! forward arrays
      call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE,1) ! 1 == forward
    endif
  endif

  end subroutine acoustic_enforce_free_surface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine acoustic_enforce_free_surface_lddrk(NGLOB_AB_LDDRK,potential_acoustic_lddrk,potential_dot_acoustic_lddrk)

  use constants
  use specfem_par, only: ibool,free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
    STACEY_INSTEAD_OF_FREE_SURFACE,BOTTOM_FREE_SURFACE,GPU_MODE
  use specfem_par_acoustic, only: ispec_is_acoustic

  implicit none

  integer,intent(in) :: NGLOB_AB_LDDRK

  ! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB_LDDRK),intent(inout) :: &
            potential_acoustic_lddrk,potential_dot_acoustic_lddrk

  ! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! safety check
  if (GPU_MODE) stop 'for GPU_MODE, enforce free surface with LDDRK not implemented yet'

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE .and. .not. BOTTOM_FREE_SURFACE) return

  ! enforce potentials to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        ! sets potentials to zero
        potential_acoustic_lddrk(iglob)     = 0._CUSTOM_REAL
        potential_dot_acoustic_lddrk(iglob) = 0._CUSTOM_REAL
      enddo
    endif

  enddo

  end subroutine acoustic_enforce_free_surface_lddrk

