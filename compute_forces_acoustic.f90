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
  
  ! time marching potentials
  if(PML) call PML_acoustic_time_march(NSPEC_AB,NGLOB_AB,ibool,&
                        potential_acoustic,potential_dot_acoustic,&
                        deltat,deltatsqover2,deltatover2,&
                        num_PML_ispec,PML_ispec,PML_damping_d,&
                        chi1,chi2,chi2_t,chi3,chi4,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot,&
                        iglob_is_PML_interface,PML_mask_ibool,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,&
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC,&
                        ispec_is_acoustic)
  
! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool, &
                        free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces, &
                        ispec_is_acoustic)

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

! distinguishes two runs: for points on MPI interfaces, and points within the partitions
  do iphase=1,2
  
    !first for points on MPI interfaces
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

! acoustic pressure term
    call acoustic_pressure( phase_is_inner, NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        ispec_is_inner,ispec_is_acoustic)
                    
    
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
      
    endif

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
        call acoustic_absorbing_boundaries(NSPEC_AB,NGLOB_AB, &
                        potential_dot_dot_acoustic,potential_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        rhostore,kappastore, &
                        ispec_is_acoustic)    

      endif
    endif
    
! elastic coupling
    if(ELASTIC_SIMULATION ) &
      call acoustic_coupling_elastic(NSPEC_AB,NGLOB_AB, &
                        ibool,displ,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! poroelastic coupling
    if(POROELASTIC_SIMULATION ) &
      call acoustic_coupling_poroelastic()
    
! sources
    call acoustic_sources(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                        xi_source,eta_source,gamma_source, &
                        hdur,hdur_gaussian,t_cmt,dt,stf,t0, &
                        sourcearrays,kappastore, &
                        ispec_is_acoustic)

! assemble all the contributions between slices using MPI
    if( phase_is_inner .eqv. .false. ) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    else
      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh,&
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    endif


  enddo

  ! divides pressure with mass matrix 
  potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:)

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
                        ibool, &
                        free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces, &
                        ispec_is_acoustic)

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

subroutine acoustic_pressure( phase_is_inner, NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_dot_acoustic, &
                        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        rhostore,jacobian,ibool, &
                        ispec_is_inner,ispec_is_acoustic )

! compute forces for the acoustic elements
!
! note that pressure is defined as:
!     p = - Chi_dot_dot  
!
  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,TINYVAL_SNGL
  use PML_par,only:PML,ispec_is_PML_inum
  implicit none
  !include "constants.h"
  integer :: NSPEC_AB,NGLOB_AB

! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_dot_acoustic

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        rhostore,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local variables
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL) temp1l,temp2l,temp3l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) rho_invl
  
  integer :: ispec,iglob,i,j,k,l

! loop over spectral elements
  do ispec = 1,NSPEC_AB

    if ( (ispec_is_inner(ispec) .eqv. phase_is_inner) ) then

      ! only elements outside PML, inside "regular" domain
      if( PML ) then
        if( ispec_is_PML_inum(ispec) > 0 ) then
         cycle
        endif
      endif
      
      if( ispec_is_acoustic(ispec) ) then

        ! gets values for element
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
            enddo
          enddo
        enddo
        ! would check if anything to do, but might lower accuracy of computation
        !if( maxval( abs( chi_elem ) ) < TINYVAL_SNGL ) cycle

        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX

              ! density (reciproc)
              rho_invl = 1.0_CUSTOM_REAL / rhostore(i,j,k,ispec) 
              
              ! derivative along x, y, z
              ! first double loop over GLL points to compute and store gradients
              ! we can merge the loops because NGLLX == NGLLY == NGLLZ
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l = 1,NGLLX
                !temp1l = temp1l + potential_acoustic(ibool(l,j,k,ispec))*hprime_xx(i,l)
                !temp2l = temp2l + potential_acoustic(ibool(i,l,k,ispec))*hprime_yy(j,l)
                !temp3l = temp3l + potential_acoustic(ibool(i,j,l,ispec))*hprime_zz(k,l)
                temp1l = temp1l + chi_elem(l,j,k)*hprime_xx(i,l)
                temp2l = temp2l + chi_elem(i,l,k)*hprime_yy(j,l)
                temp3l = temp3l + chi_elem(i,j,l)*hprime_zz(k,l)
              enddo 

              ! get derivatives of potential with respect to x, y and z
              xixl = xix(i,j,k,ispec)
              xiyl = xiy(i,j,k,ispec)
              xizl = xiz(i,j,k,ispec)
              etaxl = etax(i,j,k,ispec)
              etayl = etay(i,j,k,ispec)
              etazl = etaz(i,j,k,ispec)
              gammaxl = gammax(i,j,k,ispec)
              gammayl = gammay(i,j,k,ispec)
              gammazl = gammaz(i,j,k,ispec)
              jacobianl = jacobian(i,j,k,ispec)

              ! derivatives of potential
              dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l
              dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l
              dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l

              ! for acoustic medium
              ! also add GLL integration weights
              temp1(i,j,k) = rho_invl * wgllwgll_yz(j,k) * jacobianl* &
                            (xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
              temp2(i,j,k) = rho_invl * wgllwgll_xz(i,k) * jacobianl* &
                            (etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
              temp3(i,j,k) = rho_invl * wgllwgll_xy(i,j) * jacobianl* &
                            (gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)
            enddo
          enddo
        enddo

        ! second double-loop over GLL to compute all the terms
        do k = 1,NGLLZ
          do j = 1,NGLLZ
            do i = 1,NGLLX

              ! along x,y,z direction
              ! and assemble the contributions
              !!! can merge these loops because NGLLX = NGLLY = NGLLZ   
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l=1,NGLLX
                temp1l = temp1l + temp1(l,j,k) * hprimewgll_xx(l,i)
                temp2l = temp2l + temp2(i,l,k) * hprimewgll_yy(l,j)
                temp3l = temp3l + temp3(i,j,l) * hprimewgll_zz(l,k)
              enddo

              ! sum contributions from each element to the global values              
              iglob = ibool(i,j,k,ispec)
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                  - ( temp1l + temp2l + temp3l )

            enddo
          enddo 
        enddo

      endif ! end of test if acoustic element
    endif ! ispec_is_inner
    
  enddo ! end of loop over all spectral elements

end subroutine acoustic_pressure


!
!-------------------------------------------------------------------------------------------------
!


subroutine acoustic_absorbing_boundaries(NSPEC_AB,NGLOB_AB, &
                            potential_dot_dot_acoustic,potential_dot_acoustic, &
                            ibool,ispec_is_inner,phase_is_inner, &
                            abs_boundary_jacobian2Dw, &
                            abs_boundary_ijk,abs_boundary_ispec, &
                            num_abs_boundary_faces, &
                            rhostore,kappastore, &
                            ispec_is_acoustic)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic,&
                                                 potential_dot_acoustic
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhostore,kappastore
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! absorbing boundary surface  
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces) 
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces) 

! local parameters
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw !weight,jacobianl
  integer :: ispec,iglob,i,j,k,iface,igll
  
! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
      if( ispec_is_acoustic(ispec) ) then

        ! reference gll points on boundary face 
        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets global index
          iglob=ibool(i,j,k,ispec)

          ! determines bulk sound speed
          rhol = rhostore(i,j,k,ispec)
          cpl = sqrt( kappastore(i,j,k,ispec) / rhol )
             
          ! gets associated, weighted jacobian 
          jacobianw = abs_boundary_jacobian2Dw(igll,iface)
          
          ! Sommerfeld condition
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                              - potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
          
         enddo

      endif ! ispec_is_acoustic
    endif ! ispec_is_inner
  enddo ! num_abs_boundary_faces
  
end subroutine acoustic_absorbing_boundaries

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_coupling_elastic(NSPEC_AB,NGLOB_AB, &
                        ibool,displ,potential_dot_dot_acoustic, &
                        num_coupling_ac_el_faces, &
                        coupling_ac_el_ispec,coupling_ac_el_ijk, &
                        coupling_ac_el_normal, &
                        coupling_ac_el_jacobian2Dw, &
                        ispec_is_inner,phase_is_inner)

! returns the updated pressure array: potential_dot_dot_acoustic 
                        
  implicit none
  include 'constants.h'
  
  integer :: NSPEC_AB,NGLOB_AB

! displacement and pressure
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: displ
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic
  
! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! acoustic-elastic coupling surface
  integer :: num_coupling_ac_el_faces
  real(kind=CUSTOM_REAL) :: coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces) 
  real(kind=CUSTOM_REAL) :: coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces) 
  integer :: coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces)
  integer :: coupling_ac_el_ispec(num_coupling_ac_el_faces)   

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! local parameters
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z,displ_n
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  
  integer :: iface,igll,ispec,iglob
  integer :: i,j,k
  
! loops on all coupling faces
  do iface = 1,num_coupling_ac_el_faces

    ! gets corresponding elements
    ! (note: can be either acoustic or elastic element, no need to specify since
    !           no material properties are needed for this coupling term)
    ispec = coupling_ac_el_ispec(iface)

    if( ispec_is_inner(ispec) .eqv. phase_is_inner ) then

      ! loops over common GLL points
      do igll = 1, NGLLSQUARE
        i = coupling_ac_el_ijk(1,igll,iface)
        j = coupling_ac_el_ijk(2,igll,iface)
        k = coupling_ac_el_ijk(3,igll,iface)
        
        ! gets global index of this common GLL point
        ! (note: should be the same as for corresponding i',j',k',ispec_elastic or ispec_acoustic)
        iglob = ibool(i,j,k,ispec)
        
        ! elastic displacement on global point
        displ_x = displ(1,iglob)
        displ_y = displ(2,iglob)
        displ_z = displ(3,iglob)

        ! gets associated normal on GLL point
        ! (note convention: pointing outwards of acoustic element)
        nx = coupling_ac_el_normal(1,igll,iface)
        ny = coupling_ac_el_normal(2,igll,iface)
        nz = coupling_ac_el_normal(3,igll,iface)                   

        ! calculates displacement component along normal
        ! (normal points outwards of acoustic element)
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz    
        
        ! gets associated, weighted jacobian
        jacobianw = coupling_ac_el_jacobian2Dw(igll,iface)
        
        ! continuity of pressure and normal displacement on global point
        !
        ! note: newark time scheme together with definition of scalar potential: 
        !          pressure = - chi_dot_dot
        !          requires that this coupling term uses the updated displacement at time step [t+delta_t],
        !          which is done at the very beginning of the time loop
        !          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
        !          it also means you have to calculate and update this here first before
        !          calculating the coupling on the elastic side for the acceleration...      
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + jacobianw*displ_n
        
      enddo ! igll

    endif

  enddo ! iface
   
end subroutine acoustic_coupling_elastic

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_coupling_poroelastic()
  implicit none
 
  stop 'not yet implemented'
  
end subroutine acoustic_coupling_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_sources(NSPEC_AB,NGLOB_AB,potential_dot_dot_acoustic, &
                                  ibool,ispec_is_inner,phase_is_inner, &
                                  NSOURCES,myrank,it,islice_selected_source,ispec_selected_source,&
                                  xi_source,eta_source,gamma_source, &
                                  hdur,hdur_gaussian,t_cmt,dt,stf,t0, &
                                  sourcearrays,kappastore, &
                                  ispec_is_acoustic)

  implicit none

  include "constants.h"

  integer :: NSPEC_AB,NGLOB_AB

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: kappastore

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

! source
  integer :: NSOURCES,myrank,it
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NSOURCES) :: hdur,hdur_gaussian,t_cmt 
  double precision :: dt

  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrays 

  double precision, external :: comp_source_time_function 

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic
  
! local parameters
  double precision :: t0,f0
  double precision :: stf 
  real(kind=CUSTOM_REAL) stf_used 
  integer :: isource,iglob,ispec,i,j,k

! adds acoustic sources
  do isource = 1,NSOURCES

    !   add the source (only if this proc carries the source)
    if(myrank == islice_selected_source(isource)) then

      ispec = ispec_selected_source(isource)

      if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
  
        if( ispec_is_acoustic(ispec) ) then

          if(USE_FORCE_POINT_SOURCE) then

            ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
            iglob = ibool(nint(xi_source(isource)), &
                           nint(eta_source(isource)), &
                           nint(gamma_source(isource)), &
                           ispec)
             
            f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format
            t0 = 1.2d0/f0

            if (it == 1 .and. myrank == 0) then
              print *,'using a source of dominant frequency ',f0
              print *,'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
              print *,'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
            endif

            ! gaussian source time function
            !stf_used = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

            ! we use nu_source(:,3) here because we want a source normal to the surface.
            ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
            stf_used = 1.d10 * ( 1.d0-2.d0*PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0) ) * &
                        exp( -PI*PI*f0*f0*(dble(it-1)*DT-t0)*(dble(it-1)*DT-t0) )

            ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
            ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
            ! to add minus the source to Chi_dot_dot to get plus the source in pressure:
            
            ! acoustic source for pressure gets divided by kappa
            stf_used = stf_used / kappastore(nint(xi_source(isource)), &
                                             nint(eta_source(isource)), &
                                             nint(gamma_source(isource)),ispec)            
            
            ! source contribution
            potential_dot_dot_acoustic(iglob) = &
                        potential_dot_dot_acoustic(iglob) - stf_used
             
          else   

            ! gaussian source time 
            stf = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

            ! distinguishes between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              stf_used = sngl(stf)
            else
              stf_used = stf
            endif

            ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
            ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
            ! to add minus the source to Chi_dot_dot to get plus the source in pressure

            !     add source array
            do k=1,NGLLZ
              do j=1,NGLLY
                 do i=1,NGLLX
                    ! adds source contribution
                    ! note: acoustic source for pressure gets divided by kappa
                    iglob = ibool(i,j,k,ispec)
                    potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                            - sourcearrays(isource,1,i,j,k) * stf_used / kappastore(i,j,k,ispec)                          
                 enddo
              enddo
            enddo

          endif ! USE_FORCE_POINT_SOURCE
        endif ! ispec_is_elastic
      endif ! ispec_is_inner     
    endif ! myrank
  
  enddo ! NSOURCES

end subroutine acoustic_sources

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool, &
                        free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces, &
                        ispec_is_acoustic)
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

