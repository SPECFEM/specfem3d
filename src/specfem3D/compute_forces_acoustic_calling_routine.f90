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

  implicit none

  ! local parameters
  integer:: iphase
  logical:: phase_is_inner

  ! enforces free surface (zeroes potentials at free surface)
  if(.NOT. GPU_MODE) then
    ! on CPU
    call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,ABSORB_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

    ! adjoint simulations
    if( SIMULATION_TYPE == 3 ) &
      call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,ABSORB_INSTEAD_OF_FREE_SURFACE, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  else
    ! on GPU
    call acoustic_enforce_free_surf_cuda(Mesh_pointer,ABSORB_INSTEAD_OF_FREE_SURFACE)
  endif

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    !first for points on MPI interfaces, thus outer elements
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    ! acoustic pressure term
    if(.NOT. GPU_MODE) then
      ! on CPU
      if(USE_DEVILLE_PRODUCTS) then
        ! uses Deville (2002) optimizations
        call compute_forces_acoustic_Dev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic)
      else
        call compute_forces_acoustic_noDev(iphase,NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool,deltat, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic)
      endif

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        if(USE_DEVILLE_PRODUCTS) then
          ! uses Deville (2002) optimizations
          call compute_forces_acoustic_Dev(iphase,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        b_potential_acoustic,b_potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
                        wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic)
        else
          call compute_forces_acoustic_noDev(iphase,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool,deltat, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic)
        endif
      endif

    else
      ! on GPU
      ! includes code for SIMULATION_TYPE==3
      call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                        nspec_outer_acoustic, nspec_inner_acoustic)
    endif

    ! ! Stacey absorbing boundary conditions
    if(STACEY_ABSORBING_CONDITIONS) then
       call compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                         potential_dot_dot_acoustic,potential_dot_acoustic, &
                         ibool,ispec_is_inner,phase_is_inner, &
                         abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                         num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                         SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it,NGLOB_ADJOINT, &
                         b_potential_dot_dot_acoustic,b_reclen_potential, &
                         b_absorb_potential,b_num_abs_boundary_faces, &
                         GPU_MODE,Mesh_pointer)
    endif

    ! elastic coupling
    if(ELASTIC_SIMULATION ) then
      if( num_coupling_ac_el_faces > 0 ) then
        if( .NOT. GPU_MODE ) then
          if( SIMULATION_TYPE == 1 ) then
            ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,displ,potential_dot_dot_acoustic, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner)
          else
            ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
            ! adjoint definition: \partial_t^2 \bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                              ibool,-accel_adj_coupling,potential_dot_dot_acoustic, &
                              num_coupling_ac_el_faces, &
                              coupling_ac_el_ispec,coupling_ac_el_ijk, &
                              coupling_ac_el_normal, &
                              coupling_ac_el_jacobian2Dw, &
                              ispec_is_inner,phase_is_inner)
          endif
          ! adjoint/kernel simulations
          if( SIMULATION_TYPE == 3 ) &
            call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                            ibool,b_displ,b_potential_dot_dot_acoustic, &
                            num_coupling_ac_el_faces, &
                            coupling_ac_el_ispec,coupling_ac_el_ijk, &
                            coupling_ac_el_normal, &
                            coupling_ac_el_jacobian2Dw, &
                            ispec_is_inner,phase_is_inner)

        else
          ! on GPU
          call compute_coupling_ac_el_cuda(Mesh_pointer,phase_is_inner, &
                                              num_coupling_ac_el_faces)
        endif ! GPU_MODE
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
    call compute_add_sources_acoustic(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source, &
                        hdur,hdur_gaussian,tshift_src,dt,t0, &
                        sourcearrays,kappastore,ispec_is_acoustic,&
                        SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,b_potential_dot_dot_acoustic, &
                        NTSTEP_BETWEEN_READ_ADJSRC, &
                        GPU_MODE, Mesh_pointer)

    ! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      if(.NOT. GPU_MODE) then
        call assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      else
        ! on GPU
        call transfer_boun_pot_from_device(NGLOB_AB, Mesh_pointer, &
                                            potential_dot_dot_acoustic, &
                                            buffer_send_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh, &
                                            max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh, &
                                            ibool_interfaces_ext_mesh, &
                                            1) ! <-- 1 == fwd accel
        call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      endif

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        if(.NOT. GPU_MODE) then
          call assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                        b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
        else
          ! on GPU
          call transfer_boun_pot_from_device(NGLOB_AB, Mesh_pointer, &
                                                  b_potential_dot_dot_acoustic, &
                                                  b_buffer_send_scalar_ext_mesh,&
                                                  num_interfaces_ext_mesh, &
                                                  max_nibool_interfaces_ext_mesh, &
                                                  nibool_interfaces_ext_mesh, &
                                                  ibool_interfaces_ext_mesh, &
                                                  3) ! <-- 3 == adjoint b_accel

          call assemble_MPI_scalar_send_cuda(NPROC, &
                          b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh,&
                          my_neighbours_ext_mesh, &
                          b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

        endif
      endif

    else

      ! waits for send/receive requests to be completed and assembles values
      if(.NOT. GPU_MODE) then
        call assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      else
        ! on GPU
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        Mesh_pointer,&
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        1)
      endif

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        if(.NOT. GPU_MODE) then
          call assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                        b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
        else
          ! on GPU
          call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_potential_dot_dot_acoustic, &
                        Mesh_pointer, &
                        b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        3)
        endif
      endif
    endif !phase_is_inner

  enddo

  if(.NOT. GPU_MODE) then
    ! divides pressure with mass matrix
    potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:)

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) &
      b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_acoustic(:)
  else
    ! on GPU
    call kernel_3_a_acoustic_cuda(Mesh_pointer,NGLOB_AB)
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
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the chi_dot term which requires chi_dot_dot(t+delta)
 if( .NOT. GPU_MODE ) then
    ! corrector
    potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2*potential_dot_dot_acoustic(:)

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) &
      b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2*b_potential_dot_dot_acoustic(:)
  else
    ! on GPU
    call kernel_3_b_acoustic_cuda(Mesh_pointer,NGLOB_AB,deltatover2,b_deltatover2)
  endif

! enforces free surface (zeroes potentials at free surface)
  if(.NOT. GPU_MODE) then
    ! on CPU
    call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,ABSORB_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

    if( SIMULATION_TYPE /= 1 ) then
      potential_acoustic_adj_coupling(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
    endif

    ! adjoint simulations
    if (SIMULATION_TYPE == 3) &
      call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,ABSORB_INSTEAD_OF_FREE_SURFACE, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  else
    ! on GPU
    call acoustic_enforce_free_surf_cuda(Mesh_pointer,ABSORB_INSTEAD_OF_FREE_SURFACE)
  endif

end subroutine compute_forces_acoustic

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,ABSORB_INSTEAD_OF_FREE_SURFACE, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  implicit none
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB
  logical :: ABSORB_INSTEAD_OF_FREE_SURFACE

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
  if( ABSORB_INSTEAD_OF_FREE_SURFACE ) return

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

