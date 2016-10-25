!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! elastic solver

subroutine compute_forces_viscoelastic_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par
  use fault_solver_dynamic, only: bc_dynflt_set3d_all,SIMULATION_TYPE_DYN
  use fault_solver_kinematic, only: bc_kinflt_set_all,SIMULATION_TYPE_KIN

  implicit none

  integer:: iphase
  integer:: iface,ispec,igll,i,j,k,iglob,ispec_CPML

  ! kbai added the following two synchronizations to ensure that the displacement and velocity values
  ! at nodes on MPI interfaces stay equal on all processors that share the node.
  ! Do this only for dynamic rupture simulations
  if (SIMULATION_TYPE_DYN) then
    call synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,displ, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     my_neighbours_ext_mesh,myrank)
    call synchronize_MPI_vector_blocking_ord(NPROC,NGLOB_AB,veloc, &
                                     num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     my_neighbours_ext_mesh,myrank)
  endif


! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

! elastic term
    call compute_forces_viscoelastic(iphase,NSPEC_AB,NGLOB_AB, &
                        displ,veloc,accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,deltat,PML_CONDITIONS, &
                        one_minus_sum_beta,factor_common, &
                        one_minus_sum_beta_kappa,factor_common_kappa, &
                        alphaval,betaval,gammaval, &
                        NSPEC_ATTENUATION_AB, &
                        R_trace,R_xx,R_yy,R_xy,R_xz,R_yz, &
                        NSPEC_ATTENUATION_AB_LDDRK,R_trace_lddrk, &
                        R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                        epsilondev_trace,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                        epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT, &
                        is_moho_top,is_moho_bot, &
                        dsdx_top,dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic,.false.)

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                         ibool,iphase, &
                         abs_boundary_normal,abs_boundary_jacobian2Dw, &
                         abs_boundary_ijk,abs_boundary_ispec, &
                         num_abs_boundary_faces, &
                         veloc,rho_vp,rho_vs, &
                         ispec_is_elastic,SIMULATION_TYPE,SAVE_FORWARD, &
                         it, &
                         b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
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
                         SIMULATION_TYPE,.false., &
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
                                SIMULATION_TYPE,.false., &
                                potential_acoustic,potential_dot_acoustic)

          endif
        endif ! num_coupling_ac_el_faces
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        call compute_coupling_viscoelastic_po(NSPEC_AB,NGLOB_AB,ibool, &
                          displs_poroelastic,displw_poroelastic, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          kappaarraystore,rhoarraystore,mustore, &
                          phistore,tortstore,jacobian, &
                          displ,accel,kappastore, &
                          ANISOTROPY,NSPEC_ANISO, &
                          c11store,c12store,c13store,c14store,c15store,c16store, &
                          c22store,c23store,c24store,c25store,c26store,c33store, &
                          c34store,c35store,c36store,c44store,c45store,c46store, &
                          c55store,c56store,c66store, &
                          SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                          num_coupling_el_po_faces, &
                          coupling_el_po_ispec,coupling_po_el_ispec, &
                          coupling_el_po_ijk,coupling_po_el_ijk, &
                          coupling_el_po_normal, &
                          coupling_el_po_jacobian2Dw, &
                          iphase)
      endif

    !! CD CD
#ifndef DEBUG_COUPLED
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      ! adds source term (single-force/moment-tensor solution)
      call compute_add_sources_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                                            ibool, &
                                            NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                                            sourcearrays, &
                                            ispec_is_elastic,SIMULATION_TYPE,NSTEP, &
                                            nrec,islice_selected_rec,ispec_selected_rec, &
                                            nadj_rec_local,adj_sourcearrays, &
                                            NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY)
#endif
    !! CD CD
    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
       ! sends accel values to corresponding MPI interface neighbors
       call assemble_MPI_vector_async_send(NPROC,NGLOB_AB,accel, &
               buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
               num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
               my_neighbours_ext_mesh, &
               request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_vector_async_w_ord(NPROC,NGLOB_AB,accel, &
                            buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                            max_nibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                            my_neighbours_ext_mesh,myrank)
    endif
  enddo

!Percy , Fault boundary term B*tau is added to the assembled forces
!        which at this point are stored in the array 'accel'
  if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(accel,veloc,displ)
  if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(accel,veloc,displ)

 ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  accel(1,:) = accel(1,:)*rmassx(:)
  accel(2,:) = accel(2,:)*rmassy(:)
  accel(3,:) = accel(3,:)*rmassz(:)

! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    call compute_coupling_ocean(NSPEC_AB,NGLOB_AB, &
                                ibool,rmassx,rmassy,rmassz, &
                                rmass_ocean_load,accel, &
                                free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                num_free_surface_faces)
  endif

  ! C-PML boundary
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

          iglob=ibool(i,j,k,ispec)

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
  if (USE_LDDRK) then
    call update_veloc_elastic_lddrk()
  else
    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
  endif

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
  use pml_par
  use fault_solver_dynamic, only: bc_dynflt_set3d_all
  use fault_solver_kinematic, only: bc_kinflt_set_all

  implicit none

  integer:: iphase

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling compute_forces_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2

! elastic term
    ! adjoint simulations: backward/reconstructed wavefield
    call compute_forces_viscoelastic(iphase,NSPEC_AB,NGLOB_AB, &
                        b_displ,b_veloc,b_accel, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        kappastore,mustore,jacobian,ibool, &
                        ATTENUATION,deltat,PML_CONDITIONS, &
                        one_minus_sum_beta,factor_common, &
                        one_minus_sum_beta_kappa,factor_common_kappa, &
                        b_alphaval,b_betaval,b_gammaval, &
                        NSPEC_ATTENUATION_AB, &
                        b_R_trace,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                        NSPEC_ATTENUATION_AB_LDDRK,b_R_trace_lddrk, &
                        b_R_xx_lddrk,b_R_yy_lddrk,b_R_xy_lddrk,b_R_xz_lddrk,b_R_yz_lddrk, &
                        b_epsilondev_trace,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                        b_epsilondev_xz,b_epsilondev_yz,b_epsilon_trace_over_3, &
                        ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        SIMULATION_TYPE,COMPUTE_AND_STORE_STRAIN,NSPEC_STRAIN_ONLY, &
                        NSPEC_BOUN,NSPEC2D_MOHO,NSPEC_ADJOINT, &
                        is_moho_top,is_moho_bot, &
                        b_dsdx_top,b_dsdx_bot, &
                        ispec2D_moho_top,ispec2D_moho_bot, &
                        num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic, &
                        phase_ispec_inner_elastic,.true.)

    ! computes additional contributions
    if (iphase == 1) then
      ! adds elastic absorbing boundary term to acceleration (Stacey conditions)
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_viscoelastic_backward(NSPEC_AB, &
                         ibool,iphase, &
                         abs_boundary_ijk,abs_boundary_ispec, &
                         num_abs_boundary_faces, &
                         ispec_is_elastic,SIMULATION_TYPE, &
                         NSTEP,it,NGLOB_ADJOINT,b_accel, &
                         b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          ! backward simulations
          call compute_coupling_viscoelastic_ac(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        ibool,b_accel,b_potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        iphase, &
                        PML_CONDITIONS, &
                        SIMULATION_TYPE,.true., &
                        potential_acoustic,potential_dot_acoustic)

        endif ! num_coupling_ac_el_faces
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION) then
        stop 'poroelastic-elastic coupling error'
      endif

      ! adds source term (single-force/moment-tensor solution)
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      call compute_add_sources_viscoelastic_backward( NSPEC_AB,NGLOB_AB, &
                          ibool, &
                          NSOURCES,myrank,it,islice_selected_source,ispec_selected_source, &
                          sourcearrays, &
                          ispec_is_elastic,SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                          b_accel,NOISE_TOMOGRAPHY)
    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends accel values to corresponding MPI interface neighbors
      ! adjoint simulations
      call assemble_MPI_vector_async_send(NPROC,NGLOB_ADJOINT,b_accel, &
                 b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                 num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                 nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                 my_neighbours_ext_mesh, &
                 b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh)
    else
      ! waits for send/receive requests to be completed and assembles values
      ! adjoint simulations
      call assemble_MPI_vector_async_w_ord(NPROC,NGLOB_ADJOINT,b_accel, &
                             b_buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
                             max_nibool_interfaces_ext_mesh, &
                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                             b_request_send_vector_ext_mesh,b_request_recv_vector_ext_mesh, &
                             my_neighbours_ext_mesh,myrank)
    endif
  enddo

  ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
  ! adjoint simulations
  b_accel(1,:) = b_accel(1,:)*rmassx(:)
  b_accel(2,:) = b_accel(2,:)*rmassy(:)
  b_accel(3,:) = b_accel(3,:)*rmassz(:)

! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    call compute_coupling_ocean_backward(NSPEC_AB,NGLOB_AB, &
                                     ibool,rmassx,rmassy,rmassz, &
                                     rmass_ocean_load, &
                                     free_surface_normal,free_surface_ijk,free_surface_ispec, &
                                     num_free_surface_faces, &
                                     SIMULATION_TYPE, &
                                     NGLOB_ADJOINT,b_accel)
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
  ! adjoint simulations
  b_veloc(:,:) = b_veloc(:,:) + b_deltatover2*b_accel(:,:)

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
  use fault_solver_dynamic, only: bc_dynflt_set3d_all,SIMULATION_TYPE_DYN,synchronize_GPU
  use fault_solver_kinematic, only: bc_kinflt_set_all,SIMULATION_TYPE_KIN

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
                                          COMPUTE_AND_STORE_STRAIN,ATTENUATION,ANISOTROPY)

    ! while inner elements compute "Kernel_2", we wait for MPI to
    ! finish and transfer the boundary terms to the device asynchronously
    if (iphase == 2) then
      !daniel: todo - this avoids calling the fortran vector send from CUDA routine
      ! wait for asynchronous copy to finish
      call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_ext_mesh)

      ! sends MPI buffers
      call assemble_MPI_vector_send_cuda(NPROC, &
                  buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                  num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                  nibool_interfaces_ext_mesh, &
                  my_neighbours_ext_mesh, &
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
                                             SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                             b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                                             Mesh_pointer)
      endif

      ! acoustic coupling
      if (ACOUSTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          call compute_coupling_el_ac_cuda(Mesh_pointer,iphase, &
                                           num_coupling_ac_el_faces)
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

        call compute_coupling_viscoelastic_po(NSPEC_AB,NGLOB_AB,ibool, &
                          displs_poroelastic,displw_poroelastic, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          kappaarraystore,rhoarraystore,mustore, &
                          phistore,tortstore,jacobian, &
                          displ,accel,kappastore, &
                          ANISOTROPY,NSPEC_ANISO, &
                          c11store,c12store,c13store,c14store,c15store,c16store, &
                          c22store,c23store,c24store,c25store,c26store,c33store, &
                          c34store,c35store,c36store,c44store,c45store,c46store, &
                          c55store,c56store,c66store, &
                          SIMULATION_TYPE,NGLOB_ADJOINT,NSPEC_ADJOINT, &
                          num_coupling_el_po_faces, &
                          coupling_el_po_ispec,coupling_po_el_ispec, &
                          coupling_el_po_ijk,coupling_po_el_ijk, &
                          coupling_el_po_normal, &
                          coupling_el_po_jacobian2Dw, &
                          iphase)

        ! transfers acceleration back to GPU
        call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
      endif

      ! adds source term (single-force/moment-tensor solution)
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      call compute_add_sources_viscoelastic_GPU(NSPEC_AB, &
                          NSOURCES,myrank,it, &
                          ispec_is_elastic,SIMULATION_TYPE,NSTEP, &
                          nrec,islice_selected_rec,ispec_selected_rec, &
                          nadj_rec_local,adj_sourcearrays, &
                          NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                          Mesh_pointer)
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
        call transfer_boun_accel_from_device(Mesh_pointer, b_accel, &
                        b_buffer_send_vector_ext_mesh, &
                        3) ! 3 == adjoint b_accel
        call assemble_MPI_vector_send_cuda(NPROC, &
                        b_buffer_send_vector_ext_mesh,b_buffer_recv_vector_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh, &
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

  !Percy , Fault boundary term B*tau is added to the assembled forces
  !        which at this point are stored in the array 'accel'
  if (SIMULATION_TYPE_DYN .or. SIMULATION_TYPE_KIN) then
    ! transfers wavefields to the CPU
    ! call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel, Mesh_pointer)
    ! will remove later if GPU fault solver is fully tested

    ! adds dynamic source
    ! if (SIMULATION_TYPE_DYN) call bc_dynflt_set3d_all(accel,veloc,displ)
    ! if (SIMULATION_TYPE_KIN) call bc_kinflt_set_all(accel,veloc,displ)
    call fault_solver_gpu(Mesh_pointer,Fault_pointer,deltat,myrank,it)  ! GPU fault solver
    !call transfer_boundary_from_device_a(Mesh_pointer,nspec_outer_elastic)
    ! transfer data from mp->d_boundary to mp->h_boundary
    !call sync_copy_from_device(Mesh_pointer,2,buffer_send_vector_ext_mesh)
    ! transfer data from mp->h_boundary to send_buffer
    !call assemble_MPI_vector_send_cuda(NPROC, &
    !              buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
    !              num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
    !              nibool_interfaces_ext_mesh, &
    !              my_neighbours_ext_mesh, &
    !              request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

    ! transfers MPI buffers onto GPU
    !call transfer_boundary_to_device(NPROC,Mesh_pointer,buffer_recv_vector_ext_mesh, &
    !              num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
    !              request_recv_vector_ext_mesh)
    ! waits for send/receive requests to be completed and assembles values
    !call synchronize_MPI_vector_write_cuda(NPROC,NGLOB_AB,accel, Mesh_pointer, &
    !                  buffer_recv_vector_ext_mesh,num_interfaces_ext_mesh, &
    !                  max_nibool_interfaces_ext_mesh, &
    !                  nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
    !                  request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
    !                  1)

    if (mod(it,500) == 0 .and. it /= 0) call synchronize_GPU(it)  ! output results every 500 steps

    ! transfers acceleration back to GPU
    ! call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
    ! will remove later if GPU fault solver is fully tested
  endif

 ! multiplies with inverse of mass matrix (note: rmass has been inverted already)
 call kernel_3_a_cuda(Mesh_pointer,deltatover2,b_deltatover2,APPROXIMATE_OCEAN_LOAD)

! updates acceleration with ocean load term
  if (APPROXIMATE_OCEAN_LOAD) then
    call compute_coupling_ocean_cuda(Mesh_pointer)

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
    call kernel_3_b_cuda(Mesh_pointer,deltatover2,b_deltatover2)
  endif

end subroutine compute_forces_viscoelastic_GPU_calling

