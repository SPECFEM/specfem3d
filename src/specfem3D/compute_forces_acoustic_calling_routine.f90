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

! acoustic solver

! in case of an acoustic medium, a scalar related to grad(density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(minus_int_int_pressure) / rho
! Velocity is then:
!     v = grad(minus_int_pressure) / rho
! (minus_int_pressure being the time derivative of minus_int_int_pressure)
! and pressure is:
!     p = - minus_pressure
! (minus_pressure being the time second derivative of minus_int_int_pressure).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore minus_pressure
! is continuous, therefore minus_int_int_pressure is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.


subroutine compute_forces_acoustic()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par,only: is_CPML,spec_to_CPML,nglob_interface_PML_acoustic,&
                    b_PML_minus_int_int_pressure,b_reclen_PML_minus_int_int_pressure,&
                    PML_minus_int_int_pressure_old,PML_minus_int_int_pressure_new

  implicit none

  ! local parameters
  integer:: iphase,iface,ispec,iglob,igll,i,j,k,ispec_CPML
  logical:: phase_is_inner

  ! enforces free surface (zeroes scalars at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure,minus_int_pressure,minus_pressure, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure_lddrk,minus_int_pressure_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  endif

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if (iphase == 1) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! acoustic pressure term
    call compute_forces_acoustic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
                        minus_int_int_pressure,minus_int_pressure,minus_pressure, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic,.false.)

    ! Stacey absorbing boundary conditions
    if (STACEY_ABSORBING_CONDITIONS) then
      call compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                        minus_pressure,minus_int_pressure, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                        SIMULATION_TYPE,SAVE_FORWARD,it,b_reclen_minus_int_int_pressure, &
                        b_absorb_minus_int_int_pressure,b_num_abs_boundary_faces)
    endif

    ! elastic coupling
    if (ELASTIC_SIMULATION) then
      if (num_coupling_ac_el_faces > 0) then
        if (SIMULATION_TYPE == 1) then
          ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
          call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,displ,minus_pressure, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner,&
                              PML_CONDITIONS,SIMULATION_TYPE,.false.)

        else
          ! handles adjoint runs coupling between adjoint pressure and adjoint elastic wavefield
          ! adjoint definition: \partial_t^2 \bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger
          call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,accel,minus_pressure, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner,&
                              PML_CONDITIONS,SIMULATION_TYPE,.false.)
        endif
      endif
    endif

! poroelastic coupling
    if (POROELASTIC_SIMULATION )  then
      if (num_coupling_ac_po_faces > 0) then
        if (SIMULATION_TYPE == 1) then
          call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                        ibool,displs_poroelastic,displw_poroelastic, &
                        minus_pressure, &
                        num_coupling_ac_po_faces, &
                        coupling_ac_po_ispec,coupling_ac_po_ijk, &
                        coupling_ac_po_normal, &
                        coupling_ac_po_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
        else
          stop 'not implemented yet'
        endif
      endif
    endif

    ! sources
    call compute_add_sources_acoustic(NSPEC_AB,NGLOB_AB,minus_pressure, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        hdur,hdur_gaussian,tshift_src,dt,t0, &
                        sourcearrays,kappastore,ispec_is_acoustic,&
                        SIMULATION_TYPE,NSTEP, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,NTSTEP_BETWEEN_READ_ADJSRC)

    ! assemble all the contributions between slices using MPI
    if (phase_is_inner .eqv. .false.) then
      ! sends minus_pressure values to corresponding MPI interface neighbors (non-blocking)
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_AB,minus_pressure, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_AB,minus_pressure, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    endif !phase_is_inner

  enddo

  ! divides pressure with mass matrix
  minus_pressure(:) = minus_pressure(:) * rmass_acoustic(:)

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the pressure
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the pressure.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

! However, enforcing explicitly minus_pressure, minus_int_pressure, minus_int_int_pressure
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation

! C-PML boundary
  if (PML_CONDITIONS) then
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
!!! It is better to move this into do iphase=1,2 loop
!!!   if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
        if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
          ! reference gll points on boundary face
          ispec_CPML = spec_to_CPML(ispec)
          do igll = 1,NGLLSQUARE
            ! gets local indices for GLL point
            i = abs_boundary_ijk(1,igll,iface)
            j = abs_boundary_ijk(2,igll,iface)
            k = abs_boundary_ijk(3,igll,iface)

            iglob=ibool(i,j,k,ispec)

            minus_pressure(iglob) = 0._CUSTOM_REAL
            minus_int_pressure(iglob) = 0._CUSTOM_REAL
            minus_int_int_pressure(iglob) = 0._CUSTOM_REAL
            if (ELASTIC_SIMULATION) then
              PML_minus_int_int_pressure_old(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
              PML_minus_int_int_pressure_new(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
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
! minus_int_int_pressure(t+delta_t) = minus_int_int_pressure(t) + delta_t minus_int_pressure(t) + 1/2 delta_t**2 minus_pressure(t)
! minus_int_pressure(t+delta_t) = minus_int_pressure(t) + 1/2 delta_t minus_pressure(t) + 1/2 DELTA_T minus_pressure( T + DELTA_T )
! minus_pressure(t+delta_t) = 1/M_acoustic( -K_acoustic minus_int_int_pressure(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   minus_int_int_pressure, minus_int_pressure, minus_pressure are acoustic (fluid) scalars ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the minus_int_pressure term which requires minus_pressure(t+delta)
  ! corrector
  if (USE_LDDRK) then
    call update_minus_int_pressure_lddrk()
  else
    minus_int_pressure(:) = minus_int_pressure(:) + deltatover2*minus_pressure(:)
  endif

! enforces free surface (zeroes scalars at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure,minus_int_pressure,minus_pressure, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure_lddrk,minus_int_pressure_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  endif

  if (SIMULATION_TYPE /= 1) then
    minus_int_int_pressure_adj_coupling(:) = minus_int_int_pressure(:) &
                            + deltat * minus_int_pressure(:) &
                            + deltatsqover2 * minus_pressure(:)
  endif

  if (PML_CONDITIONS) then
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      if (nglob_interface_PML_acoustic > 0) then
        call save_pressure_scalar_on_pml_interface(minus_int_int_pressure,minus_int_pressure,minus_pressure,&
                    nglob_interface_PML_acoustic,b_PML_minus_int_int_pressure,b_reclen_PML_minus_int_int_pressure)
      endif
    endif
  endif

end subroutine compute_forces_acoustic

!
!=====================================================================

! acoustic solver

! in case of an acoustic medium, a scalar related to grad(density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(minus_int_int_pressure) / rho
! Velocity is then:
!     v = grad(minus_int_pressure) / rho
! (minus_int_pressure being the time derivative of minus_int_int_pressure)
! and pressure is:
!     p = - minus_pressure
! (minus_pressure being the time second derivative of minus_int_int_pressure).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore minus_pressure
! is continuous, therefore minus_int_int_pressure is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.


subroutine compute_forces_acoustic_backward()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: iphase
  logical:: phase_is_inner

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling compute_forces_acoustic_backward() with wrong SIMULATION_TYPE')

  ! adjoint simulations
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      b_minus_int_int_pressure,b_minus_int_pressure,b_minus_pressure, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if (iphase == 1) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! adjoint simulations
    call compute_forces_acoustic_noDev(iphase,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                      b_minus_int_int_pressure,b_minus_int_pressure,b_minus_pressure, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      hprime_xx,hprime_yy,hprime_zz, &
                      hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                      wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                      rhostore,jacobian,ibool, &
                      num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                      phase_ispec_inner_acoustic,.true.)

    ! Stacey absorbing boundary conditions
    if (STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic_backward(NSPEC_AB, &
                            ibool,ispec_is_inner,phase_is_inner, &
                            abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces,ispec_is_acoustic,&
                            SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                            b_minus_pressure,b_reclen_minus_int_int_pressure, &
                            b_absorb_minus_int_int_pressure,b_num_abs_boundary_faces)
    endif

    ! elastic coupling
    if (ELASTIC_SIMULATION) then
      if (num_coupling_ac_el_faces > 0) then
        ! adjoint/kernel simulations
        call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                          ibool,b_displ,b_minus_pressure, &
                          num_coupling_ac_el_faces, &
                          coupling_ac_el_ispec,coupling_ac_el_ijk, &
                          coupling_ac_el_normal, &
                          coupling_ac_el_jacobian2Dw, &
                          ispec_is_inner,phase_is_inner,&
                          PML_CONDITIONS,SIMULATION_TYPE,.true.)
      endif
    endif

! poroelastic coupling
    if (POROELASTIC_SIMULATION )  then
      if (num_coupling_ac_po_faces > 0) then
          stop 'coupling acoustic-poroelastic domains not implemented yet...'
      endif
    endif

    ! sources
    call compute_add_sources_acoustic_backward(NSPEC_AB, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  sourcearrays,kappastore,ispec_is_acoustic,&
                                  SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                                  b_minus_pressure)

    ! assemble all the contributions between slices using MPI
    if (phase_is_inner .eqv. .false.) then
      ! sends b_minus_pressure values to corresponding MPI interface neighbors (non-blocking)
      ! adjoint simulations
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_ADJOINT,b_minus_pressure, &
                      b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                      my_neighbours_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

    else
      ! adjoint simulations
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_ADJOINT,b_minus_pressure, &
                      b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
    endif !phase_is_inner

  enddo

  ! divides pressure with mass matrix
  ! adjoint simulations
  b_minus_pressure(:) = b_minus_pressure(:) * rmass_acoustic(:)

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! minus_int_int_pressure(t+delta_t) = minus_int_int_pressure(t) + delta_t minus_int_pressure(t) + 1/2 delta_t**2 minus_pressure(t)
! minus_int_pressure(t+delta_t) = minus_int_pressure(t) + 1/2 delta_t minus_pressure(t) + 1/2 DELTA_T minus_pressure( T + DELTA_T )
! minus_pressure(t+delta_t) = 1/M_acoustic( -K_acoustic minus_int_int_pressure(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   minus_int_int_pressure, minus_int_pressure, minus_pressure are acoustic (fluid) scalars ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the minus_int_pressure term which requires minus_pressure(t+delta)

  ! corrector
  ! adjoint simulations
  b_minus_int_pressure(:) = b_minus_int_pressure(:) + b_deltatover2*b_minus_pressure(:)

! enforces free surface (zeroes scalars at free surface)
  ! adjoint simulations
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      b_minus_int_int_pressure,b_minus_int_pressure,b_minus_pressure, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

end subroutine compute_forces_acoustic_backward
!
!-------------------------------------------------------------------------------------------------
!
subroutine compute_forces_acoustic_GPU()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: iphase
  logical:: phase_is_inner

  ! enforces free surface (zeroes scalars at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if (iphase == 1) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                      nspec_outer_acoustic, nspec_inner_acoustic)

    ! ! Stacey absorbing boundary conditions
    if (STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic_GPU(phase_is_inner,num_abs_boundary_faces,&
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                            b_reclen_minus_int_int_pressure,b_absorb_minus_int_int_pressure, &
                            b_num_abs_boundary_faces,Mesh_pointer)
    endif

    ! elastic coupling
    if (ELASTIC_SIMULATION) then
      if (num_coupling_ac_el_faces > 0) then
        ! on GPU
        call compute_coupling_ac_el_cuda(Mesh_pointer,phase_is_inner, &
                                         num_coupling_ac_el_faces)
      endif
    endif

! poroelastic coupling
    if (POROELASTIC_SIMULATION )  then
      if (num_coupling_ac_po_faces > 0) then
        if (SIMULATION_TYPE == 1) then
          call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                        ibool,displs_poroelastic,displw_poroelastic, &
                        minus_pressure, &
                        num_coupling_ac_po_faces, &
                        coupling_ac_po_ispec,coupling_ac_po_ijk, &
                        coupling_ac_po_normal, &
                        coupling_ac_po_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
        else
          stop 'not implemented yet'
        endif
        if (SIMULATION_TYPE == 3) &
          stop 'not implemented yet'
      endif
    endif

    ! sources
    call compute_add_sources_acoustic_GPU(NSPEC_AB,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,&
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  ispec_is_acoustic,SIMULATION_TYPE,NSTEP, &
                                  nrec,islice_selected_rec,ispec_selected_rec, &
                                  nadj_rec_local,adj_sourcearrays, &
                                  NTSTEP_BETWEEN_READ_ADJSRC,Mesh_pointer)

    ! assemble all the contributions between slices using MPI
    if (phase_is_inner .eqv. .false.) then
      ! sends minus_pressure values to corresponding MPI interface neighbors (non-blocking)
      call transfer_boun_pot_from_device(Mesh_pointer, &
                                         minus_pressure, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! <-- 1 == fwd accel
      call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_minus_pressure, &
                                           b_buffer_send_scalar_ext_mesh,&
                                           3) ! <-- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                          b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh,&
                          my_neighbours_ext_mesh, &
                          b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

      endif

    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,minus_pressure, &
                        Mesh_pointer,&
                        buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        1)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_minus_pressure, &
                        Mesh_pointer, &
                        b_buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh, &
                        3)
      endif
    endif !phase_is_inner

  enddo

 ! divides pressure with mass matrix
  call kernel_3_a_acoustic_cuda(Mesh_pointer)

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! minus_int_int_pressure(t+delta_t) = minus_int_int_pressure(t) + delta_t minus_int_pressure(t) + 1/2 delta_t**2 minus_pressure(t)
! minus_int_pressure(t+delta_t) = minus_int_pressure(t) + 1/2 delta_t minus_pressure(t) + 1/2 DELTA_T minus_pressure( T + DELTA_T )
! minus_pressure(t+delta_t) = 1/M_acoustic( -K_acoustic minus_int_int_pressure(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   minus_int_int_pressure, minus_int_pressure, minus_pressure are acoustic (fluid) scalars ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
! updates the minus_int_pressure term which requires minus_pressure(t+delta)
  call kernel_3_b_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2)

! enforces free surface (zeroes scalars at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

end subroutine compute_forces_acoustic_GPU
!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure,minus_int_pressure,minus_pressure, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE

! acoustic scalars
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        minus_int_int_pressure,minus_int_pressure,minus_pressure

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! free surface
  integer :: num_free_surface_faces
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE ) return

! enforce scalars to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        ! sets scalars to zero
        minus_int_int_pressure(iglob)         = 0._CUSTOM_REAL
        minus_int_pressure(iglob)     = 0._CUSTOM_REAL
        minus_pressure(iglob) = 0._CUSTOM_REAL
      enddo
    endif

  enddo

end subroutine acoustic_enforce_free_surface

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        minus_int_int_pressure_lddrk,minus_int_pressure_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB_LDDRK
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE

! acoustic scalars
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB_LDDRK) :: &
            minus_int_int_pressure_lddrk,minus_int_pressure_lddrk

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! free surface
  integer :: num_free_surface_faces
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE ) return

! enforce scalars to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        ! sets scalars to zero
        minus_int_int_pressure_lddrk(iglob)         = 0._CUSTOM_REAL
        minus_int_pressure_lddrk(iglob)     = 0._CUSTOM_REAL
      enddo
    endif

  enddo

end subroutine acoustic_enforce_free_surface_lddrk

