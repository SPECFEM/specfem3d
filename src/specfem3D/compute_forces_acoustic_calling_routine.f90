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


subroutine compute_forces_acoustic()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par,only: spec_to_CPML,is_CPML,rmemory_coupling_ac_el_displ,nglob_interface_PML_acoustic,&
                    b_PML_potential,b_reclen_PML_potential,potential_dot_dot_acoustic_old,potential_acoustic_old
  implicit none

  ! local parameters
  integer:: iphase,iface,ispec,iglob,igll,i,j,k
  logical:: phase_is_inner

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! acoustic pressure term
    if(USE_DEVILLE_PRODUCTS) then
      ! uses Deville (2002) optimizations
      call compute_forces_acoustic_Dev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic,.false.)
    else
      call compute_forces_acoustic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic,.false.)
    endif

    ! ! Stacey absorbing boundary conditions
    if(STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                         potential_dot_dot_acoustic,potential_dot_acoustic, &
                         ibool,ispec_is_inner,phase_is_inner, &
                         abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                         num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                         SIMULATION_TYPE,SAVE_FORWARD,it,b_reclen_potential, &
                         b_absorb_potential,b_num_abs_boundary_faces)
    endif

    ! elastic coupling
    if(ELASTIC_SIMULATION ) then
      if( num_coupling_ac_el_faces > 0 ) then
        if( SIMULATION_TYPE == 1 ) then
          ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
          call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,displ,potential_dot_dot_acoustic, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner,&
                              PML_CONDITIONS,spec_to_CPML,is_CPML,&
                              rmemory_coupling_ac_el_displ,&
                              SIMULATION_TYPE,.false.)

        else
          ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
          ! adjoint definition: \partial_t^2 \bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger
          call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,-accel,potential_dot_dot_acoustic, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner,&
                              PML_CONDITIONS,spec_to_CPML,is_CPML,&
                              rmemory_coupling_ac_el_displ,&
                              SIMULATION_TYPE,.false.)
        endif
      endif
    endif

! poroelastic coupling
    if(POROELASTIC_SIMULATION )  then
      if( num_coupling_ac_po_faces > 0 ) then
        if( SIMULATION_TYPE == 1 ) then
          call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                        ibool,displs_poroelastic,displw_poroelastic, &
                        potential_dot_dot_acoustic, &
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
    call compute_add_sources_acoustic(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        hdur,hdur_gaussian,tshift_src,dt,t0, &
                        sourcearrays,kappastore,ispec_is_acoustic,&
                        SIMULATION_TYPE,NSTEP, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,NTSTEP_BETWEEN_READ_ADJSRC)

    ! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    endif !phase_is_inner

  enddo

  ! divides pressure with mass matrix
  potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:)

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

! However, enforcing explicitly potential_dot_dot_acoustic, potential_dot_acoustic, potential_acoustic
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation

! C-PML boundary
  if(PML_CONDITIONS)then
    do iface=1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
!!! It is better to move this into do iphase=1,2 loop
!!!   if(ispec_is_inner(ispec) .eqv. phase_is_inner) then
        if(ispec_is_acoustic(ispec) .and. is_CPML(ispec) ) then
          ! reference gll points on boundary face
          do igll = 1,NGLLSQUARE
            ! gets local indices for GLL point
            i = abs_boundary_ijk(1,igll,iface)
            j = abs_boundary_ijk(2,igll,iface)
            k = abs_boundary_ijk(3,igll,iface)

            iglob=ibool(i,j,k,ispec)

            potential_dot_dot_acoustic(iglob) = 0.0
            potential_dot_acoustic(iglob) = 0.0
            potential_acoustic(iglob) = 0.0
            if(ELASTIC_SIMULATION ) then
              potential_dot_dot_acoustic_old(iglob) = 0.0
              potential_acoustic_old(iglob) = 0.0
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
  potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2*potential_dot_dot_acoustic(:)

! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  if( SIMULATION_TYPE /= 1 ) then
    potential_acoustic_adj_coupling(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
  endif

  if(PML_CONDITIONS)then
    if(SIMULATION_TYPE == 1 .and. SAVE_FORWARD)then
      if(nglob_interface_PML_acoustic > 0)then
        call save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,&
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)
      endif
    endif
  endif

end subroutine compute_forces_acoustic

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


subroutine compute_forces_acoustic_bpwf()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par,only: spec_to_CPML,is_CPML,rmemory_coupling_ac_el_displ
  implicit none

  ! local parameters
  integer:: iphase
  logical:: phase_is_inner

  ! checks
  if( SIMULATION_TYPE /= 3 ) &
    call exit_MPI(myrank,'error calling compute_forces_acoustic_bpwf() with wrong SIMULATION_TYPE')

  ! adjoint simulations
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! adjoint simulations
    if(USE_DEVILLE_PRODUCTS) then
      ! uses Deville (2002) optimizations
      call compute_forces_acoustic_Dev(iphase,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                      wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                      rhostore,jacobian,ibool, &
                      num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                      phase_ispec_inner_acoustic,.true.)
    else
      call compute_forces_acoustic_noDev(iphase,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      hprime_xx,hprime_yy,hprime_zz, &
                      hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                      wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                      rhostore,jacobian,ibool, &
                      num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                      phase_ispec_inner_acoustic,.true.)
    endif

    ! ! Stacey absorbing boundary conditions
    if(STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic_bpwf(NSPEC_AB, &
                            ibool,ispec_is_inner,phase_is_inner, &
                            abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces,ispec_is_acoustic,&
                            SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                            b_potential_dot_dot_acoustic,b_reclen_potential, &
                            b_absorb_potential,b_num_abs_boundary_faces)
    endif

    ! elastic coupling
    if(ELASTIC_SIMULATION ) then
      if( num_coupling_ac_el_faces > 0 ) then
        ! adjoint/kernel simulations
        call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                          ibool,b_displ,b_potential_dot_dot_acoustic, &
                          num_coupling_ac_el_faces, &
                          coupling_ac_el_ispec,coupling_ac_el_ijk, &
                          coupling_ac_el_normal, &
                          coupling_ac_el_jacobian2Dw, &
                          ispec_is_inner,phase_is_inner,&
                          PML_CONDITIONS,spec_to_CPML,is_CPML,&
                          rmemory_coupling_ac_el_displ,&
                          SIMULATION_TYPE,.true.)
      endif
    endif

! poroelastic coupling
    if(POROELASTIC_SIMULATION )  then
      if( num_coupling_ac_po_faces > 0 ) then
          stop 'coupling acoustic-poroelastic domains not implemented yet...'
      endif
    endif

    ! sources
    call compute_add_sources_acoustic_bpwf(NSPEC_AB, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  hdur,hdur_gaussian,tshift_src,dt,t0, &
                                  sourcearrays,kappastore,ispec_is_acoustic,&
                                  SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                                  b_potential_dot_dot_acoustic)

    ! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends b_potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      ! adjoint simulations
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                      b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                      my_neighbours_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

    else
      ! adjoint simulations
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                      b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
    endif !phase_is_inner

  enddo

  ! divides pressure with mass matrix
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
  ! adjoint simulations
  b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2*b_potential_dot_dot_acoustic(:)

! enforces free surface (zeroes potentials at free surface)
  ! adjoint simulations
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

end subroutine compute_forces_acoustic_bpwf
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

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                      nspec_outer_acoustic, nspec_inner_acoustic)

    ! ! Stacey absorbing boundary conditions
    if(STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic_GPU(phase_is_inner,num_abs_boundary_faces,&
                            SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                            b_reclen_potential,b_absorb_potential, &
                            b_num_abs_boundary_faces,Mesh_pointer)
    endif

    ! elastic coupling
    if(ELASTIC_SIMULATION ) then
      if( num_coupling_ac_el_faces > 0 ) then
        ! on GPU
        call compute_coupling_ac_el_cuda(Mesh_pointer,phase_is_inner, &
                                         num_coupling_ac_el_faces)
      endif
    endif

! poroelastic coupling
    if(POROELASTIC_SIMULATION )  then
      if( num_coupling_ac_po_faces > 0 ) then
        if( SIMULATION_TYPE == 1 ) then
          call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                        ibool,displs_poroelastic,displw_poroelastic, &
                        potential_dot_dot_acoustic, &
                        num_coupling_ac_po_faces, &
                        coupling_ac_po_ispec,coupling_ac_po_ijk, &
                        coupling_ac_po_normal, &
                        coupling_ac_po_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
        else
          stop 'not implemented yet'
        endif
        if( SIMULATION_TYPE == 3 ) &
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
    if( phase_is_inner .eqv. .false. ) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call transfer_boun_pot_from_device(Mesh_pointer, &
                                         potential_dot_dot_acoustic, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! <-- 1 == fwd accel
      call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_potential_dot_dot_acoustic, &
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
      call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        Mesh_pointer,&
                        buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        1)

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_potential_dot_dot_acoustic, &
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
! updates the chi_dot term which requires chi_dot_dot(t+delta)
  call kernel_3_b_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2)

! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

end subroutine compute_forces_acoustic_GPU
!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE

! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! free surface
  integer :: num_free_surface_faces
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if( STACEY_INSTEAD_OF_FREE_SURFACE ) return

! enforce potentials to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if( ispec_is_acoustic(ispec) ) then

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

end subroutine acoustic_enforce_free_surface

