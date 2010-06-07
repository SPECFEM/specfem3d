!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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
  use PML_par
  use PML_par_acoustic
  implicit none
  ! local parameters
  integer:: iphase
  logical:: phase_is_inner
  
  
! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  ! adjoint simulations
  if( SIMULATION_TYPE == 3 ) &
    call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  

  if(PML) call PML_acoustic_enforce_free_srfc(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic, &
                        num_PML_ispec,PML_ispec,&
                        chi1,chi2,chi2_t,chi3,chi4,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,&
                        chi3_dot_dot,chi4_dot_dot)             

! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase=1,2
  
    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

! acoustic pressure term
    call compute_forces_acoustic_pot( iphase, NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic )
                    
    ! adjoint simulations
    if( SIMULATION_TYPE == 3 ) &
      call compute_forces_acoustic_pot( iphase, NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        b_potential_acoustic,b_potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic,&
                        phase_ispec_inner_acoustic )

    
    if(PML) then
      call compute_forces_acoustic_PML(NSPEC_AB,NGLOB_AB, &
                        ibool,ispec_is_inner,phase_is_inner, &                        
                        rhostore,ispec_is_acoustic,potential_acoustic, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian,&
                        wxgll,wygll,wzgll,&
                        PML_damping_dprime,num_PML_ispec,&
                        PML_ispec,PML_normal,&
                        chi1_dot_dot,chi2_t_dot_dot,&
                        chi3_dot_dot,chi4_dot_dot)

      ! couples potential_dot_dot with PML interface contributions
      call PML_acoustic_interface_coupling(phase_is_inner,NSPEC_AB,NGLOB_AB,&
                        potential_dot_dot_acoustic,&
                        ibool,ispec_is_inner,ispec_is_acoustic,&
                        num_PML_ispec,PML_ispec,iglob_is_PML_interface,&
                        chi1_dot_dot,chi3_dot_dot,chi4_dot_dot)
      
    endif ! PML

! absorbing boundaries
    if(ABSORBING_CONDITIONS) then
      if( PML .and. PML_USE_SOMMERFELD ) then
        ! adds a Sommerfeld condition on the domain's absorbing boundaries
        call PML_acoustic_abs_boundaries(phase_is_inner,NSPEC_AB,NGLOB_AB,&
                        abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        kappastore,ibool,ispec_is_inner, &
                        rhostore,ispec_is_acoustic,&
                        potential_dot_acoustic,potential_dot_dot_acoustic,&
                        num_PML_ispec,PML_ispec,ispec_is_PML_inum,&
                        chi1_dot,chi2_t,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi3_dot_dot,chi4_dot_dot)      
      else
        call compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                        potential_dot_dot_acoustic,potential_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it,myrank,NGLOB_ADJOINT, &
                        b_potential_dot_dot_acoustic,b_reclen_potential, &
                        b_absorb_potential,b_num_abs_boundary_faces)    
      endif
    endif
    
! elastic coupling
    if(ELASTIC_SIMULATION ) then
      call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                        ibool,displ,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) &
        call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                        ibool,b_displ,b_potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)
    endif

! poroelastic coupling 
! not implemented yet
    !if(POROELASTIC_SIMULATION ) &
    !  call compute_coupling_acoustic_poro()
    
! sources
    call compute_add_sources_acoustic(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source, &
                        hdur,hdur_gaussian,t_cmt,dt,t0, &
                        sourcearrays,kappastore,ispec_is_acoustic,&
                        SIMULATION_TYPE,NSTEP,NGLOB_ADJOINT, &
                        nrec,islice_selected_rec,ispec_selected_rec, &
                        nadj_rec_local,adj_sourcearrays,b_potential_dot_dot_acoustic )

! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) &  
        call assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                        b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) &  
      call assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                        b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

    endif


  enddo

  ! divides pressure with mass matrix 
  potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:)

  ! adjoint simulations  
  if (SIMULATION_TYPE == 3) &
    b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_acoustic(:)


  if(PML) then
    ! divides local contributions with mass term
    call PML_acoustic_mass_update(NSPEC_AB,NGLOB_AB,&
                        ispec_is_acoustic,rmass_acoustic,ibool,&
                        num_PML_ispec,PML_ispec,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot)

    ! Newark time scheme corrector terms
    call PML_acoustic_time_corrector(NSPEC_AB,ispec_is_acoustic,deltatover2,&
                        num_PML_ispec,PML_ispec,PML_damping_d,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot)                        
  endif


! update velocity
! note: Newark finite-difference time scheme with acoustic domains:
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
  potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2*potential_dot_dot_acoustic(:)

  ! adjoint simulations  
  if (SIMULATION_TYPE == 3) &
    b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + deltatover2*b_potential_dot_dot_acoustic(:)

  ! updates potential_dot_acoustic and potential_dot_dot_acoustic inside PML region for plotting seismograms/movies
  if(PML) call PML_acoustic_update_potentials(NGLOB_AB,NSPEC_AB, &
                        ibool,ispec_is_acoustic, &
                        potential_dot_acoustic,potential_dot_dot_acoustic,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC,&
                        num_PML_ispec,PML_ispec,iglob_is_PML_interface,&
                        PML_mask_ibool,PML_damping_d,&
                        chi1,chi2,chi2_t,chi3,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi3_dot_dot,chi4_dot_dot)


! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  ! adjoint simulations  
  if (SIMULATION_TYPE == 3) &
    call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT, &
                        b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
                        

  if(PML) call PML_acoustic_enforce_free_srfc(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces, &
                        ispec_is_acoustic, &
                        num_PML_ispec,PML_ispec,&
                        chi1,chi2,chi2_t,chi3,chi4,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,&
                        chi3_dot_dot,chi4_dot_dot)             


end subroutine compute_forces_acoustic


!
!-------------------------------------------------------------------------------------------------
!


subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  implicit none 
  include 'constants.h'

  integer :: NSPEC_AB,NGLOB_AB

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

