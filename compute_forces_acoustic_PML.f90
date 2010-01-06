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

subroutine compute_forces_acoustic_PML(NSPEC_AB,NGLOB_AB, &
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

  use constants,only: NGLLX,NGLLY,NGLLZ,NDIM,TINYVAL_SNGL,CUSTOM_REAL
  implicit none

  integer :: NSPEC_AB,NGLOB_AB

  ! potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_acoustic 
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          PML_damping_dprime
  integer,dimension(num_PML_ispec):: PML_ispec
  real(kind=CUSTOM_REAL),dimension(NDIM,num_PML_ispec):: PML_normal
          
  ! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner
  
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rhostore 

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLY) :: wygll
  double precision, dimension(NGLLZ) :: wzgll

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1_n,temp2_n,temp3_n
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1_p,temp2_p,temp3_p
  real(kind=CUSTOM_REAL) :: rho_invl 
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp3l
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl  
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: dpotentialdxl_n,dpotentialdyl_n,dpotentialdzl_n
  real(kind=CUSTOM_REAL) :: dpotentialdxl_p,dpotentialdyl_p,dpotentialdzl_p  
  real(kind=CUSTOM_REAL) :: nx,ny,nz,grad_n,dprime,weights  
  integer :: ispec,iglob,i,j,k,l,ispecPML 
   
  ! loops over all PML elements           
  do ispecPML=1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)
    
    ! checks with MPI interface flag
    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
      ! only acoustic part
      if( ispec_is_acoustic(ispec) ) then

        ! gets values for element
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              chi_elem(i,j,k) = potential_acoustic(ibool(i,j,k,ispec))
            enddo
          enddo
        enddo
        ! checks if anything to do
        if( maxval( abs( chi_elem ) ) < TINYVAL_SNGL ) cycle

        ! PML normal 
        nx = PML_normal(1,ispecPML)
        ny = PML_normal(2,ispecPML)
        nz = PML_normal(3,ispecPML)

        ! calculates terms:
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
              ! \npartial_i \chi
              dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l
              dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l
              dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l

              ! splits derivatives of potential into normal and parallel components
              ! dpotential normal to PML plane
              ! \hat{n} \partial_n \chi
              grad_n = dpotentialdxl*nx + dpotentialdyl*ny + dpotentialdzl*nz              
              dpotentialdxl_n = nx * grad_n
              dpotentialdyl_n = ny * grad_n
              dpotentialdzl_n = nz * grad_n              
              
              
              ! dpotential parallel to plane                            
              ! \nabla^{parallel} \chi
              dpotentialdxl_p = dpotentialdxl - dpotentialdxl_n
              dpotentialdyl_p = dpotentialdyl - dpotentialdyl_n
              dpotentialdzl_p = dpotentialdzl - dpotentialdzl_n
              
              ! normal incidence term: ( 1/rho J \hat{n} \partial_n \chi )
              ! (note: we can add two weights at this point to save some computations )
              temp1_n(i,j,k) = rho_invl * jacobianl * dpotentialdxl_n 
              temp2_n(i,j,k) = rho_invl * jacobianl * dpotentialdyl_n  
              temp3_n(i,j,k) = rho_invl * jacobianl * dpotentialdzl_n 
                            
              ! parallel incidence 1/rho J \nabla^{parallel} \chi
              temp1_p(i,j,k) = rho_invl * jacobianl * dpotentialdxl_p  
              temp2_p(i,j,k) = rho_invl * jacobianl * dpotentialdyl_p 
              temp3_p(i,j,k) = rho_invl * jacobianl * dpotentialdzl_p 
            enddo
          enddo
        enddo

        ! second double-loop over GLL to compute all the terms
        do k = 1,NGLLZ
          do j = 1,NGLLZ
            do i = 1,NGLLX              
              
              iglob = ibool(i,j,k,ispec)
              
              ! 1. split term:
              !-----------------
              ! normal derivative of w dotted with normal dpotential
              ! ( \hat{n} \nabla_n w ) \cdot ( 1/rho \hat{n} \nabla_n \chi )
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l=1,NGLLX
                ! derivatives
                xixl = xix(l,j,k,ispec)
                xiyl = xiy(l,j,k,ispec)
                xizl = xiz(l,j,k,ispec)
                ! note: hprimewgll_xx(l,i) = hprime_xx(l,i)*wxgll(l)
                !          don't confuse order of indices in hprime_xx: they are l and i 
                !           -> lagrangian (hprime) function i evaluated at point xi_{ l }
                temp1l = temp1l + hprimewgll_xx(l,i)   &
                                  *(nx*temp1_n(l,j,k)+ny*temp2_n(l,j,k)+nz*temp3_n(l,j,k)) &
                                  *(nx*xixl+ny*xiyl+nz*xizl)

                etaxl = etax(i,l,k,ispec)
                etayl = etay(i,l,k,ispec)
                etazl = etaz(i,l,k,ispec)                                  

                temp2l = temp2l + hprimewgll_yy(l,j)  &
                                  *(nx*temp1_n(i,l,k)+ny*temp2_n(i,l,k)+nz*temp3_n(i,l,k)) &
                                  *(nx*etaxl+ny*etayl+nz*etazl)

                gammaxl = gammax(i,j,l,ispec)
                gammayl = gammay(i,j,l,ispec)
                gammazl = gammaz(i,j,l,ispec)                                  

                temp3l = temp3l + hprimewgll_zz(l,k)  &
                                  *(nx*temp1_n(i,j,l)+ny*temp2_n(i,j,l)+nz*temp3_n(i,j,l)) &
                                  *(nx*gammaxl+ny*gammayl+nz*gammazl)
              enddo
              temp1l = temp1l * wgllwgll_yz(j,k)      
              temp2l = temp2l * wgllwgll_xz(i,k)      
              temp3l = temp3l * wgllwgll_xy(i,j)      
              
              chi1_dot_dot(i,j,k,ispecPML) = - (temp1l + temp2l + temp3l)

              ! 2. split term:
              !-----------------
              ! dprime times normal w dotted with normal dpotential
              ! w dprime \hat{n} \cdot ( 1/rho \hat{n} \nabla_n \chi )

              weights = wxgll(i)*wygll(j)*wzgll(k)
              
              temp1l = nx*temp1_n(i,j,k)*weights
              temp2l = ny*temp2_n(i,j,k)*weights
              temp3l = nz*temp3_n(i,j,k)*weights

              dprime = PML_damping_dprime(i,j,k,ispecPML)
              
              ! contribution has negative sign?
              chi2_t_dot_dot(i,j,k,ispecPML) = - dprime*(temp1l + temp2l + temp3l )


              ! 3. split term:
              !-----------------
              ! parallel derivative of w dotted with normal dpotential
              ! ( \nabla^{parallel} w ) \cdot ( 1/rho \hat{n} \nabla_n \chi )      
              ! and
              ! normal derivative of w dotted with parallel dpotential
              ! ( \hat{n} \nabla_n w ) \cdot ( 1/rho \nabla_{parallel} \chi )                    
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l=1,NGLLX
                ! derivatives
                xixl = xix(l,j,k,ispec)
                xiyl = xiy(l,j,k,ispec)
                xizl = xiz(l,j,k,ispec)
                etaxl = etax(i,l,k,ispec)
                etayl = etay(i,l,k,ispec)
                etazl = etaz(i,l,k,ispec)
                gammaxl = gammax(i,j,l,ispec)
                gammayl = gammay(i,j,l,ispec)
                gammazl = gammaz(i,j,l,ispec)

                ! normal derivative of w dotted with parallel dpotential
                temp1l = temp1l + hprimewgll_xx(l,i)  &
                        *(nx*temp1_p(l,j,k)+ny*temp2_p(l,j,k)+nz*temp3_p(l,j,k)) &
                        *(nx*xixl+ny*xiyl+nz*xizl)
                                  
                temp2l = temp2l + hprimewgll_yy(l,j)  &
                        *(nx*temp1_p(i,l,k)+ny*temp2_p(i,l,k)+nz*temp3_p(i,l,k)) &
                        *(nx*etaxl+ny*etayl+nz*etazl)
                                  
                temp3l = temp3l + hprimewgll_zz(l,k)  &
                        *(nx*temp1_p(i,j,l)+ny*temp2_p(i,j,l)+nz*temp3_p(i,j,l)) &
                        *(nx*gammaxl+ny*gammayl+nz*gammazl)


                ! parallel derivative of w dotted with normal dpotential
                temp1l = temp1l + hprimewgll_xx(l,i)  &
                        *( (xixl - nx*(nx*xixl+ny*xiyl+nz*xizl))*temp1_n(l,j,k) &
                          +(xiyl - ny*(nx*xixl+ny*xiyl+nz*xizl))*temp2_n(l,j,k) & 
                          +(xizl - nz*(nx*xixl+ny*xiyl+nz*xizl))*temp3_n(l,j,k) )

                temp2l = temp2l + hprimewgll_yy(l,j)  &
                        *( (etaxl - nx*(nx*etaxl+ny*etayl+nz*etazl))*temp1_n(i,l,k) &
                          +(etayl - ny*(nx*etaxl+ny*etayl+nz*etazl))*temp2_n(i,l,k) & 
                          +(etazl - nz*(nx*etaxl+ny*etayl+nz*etazl))*temp3_n(i,l,k) )

                temp3l = temp3l + hprimewgll_zz(l,k)  &
                        *( (gammaxl - nx*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp1_n(i,j,l) &
                          +(gammayl - ny*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp2_n(i,j,l) & 
                          +(gammazl - nz*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp3_n(i,j,l) )
              enddo
              temp1l = temp1l * wgllwgll_yz(j,k)      
              temp2l = temp2l * wgllwgll_xz(i,k)      
              temp3l = temp3l * wgllwgll_xy(i,j)      
                         
              chi3_dot_dot(i,j,k,ispecPML) = - (temp1l + temp2l + temp3l)
              

              ! 4. split term:
              !-----------------
              ! parallel derivative of w dotted with parallel dpotential
              ! ( \nabla_{parallel} w ) \cdot ( 1/rho \nabla_{parallel} \chi )
              temp1l = 0._CUSTOM_REAL
              temp2l = 0._CUSTOM_REAL
              temp3l = 0._CUSTOM_REAL
              do l=1,NGLLX
                ! derivatives
                xixl = xix(l,j,k,ispec)
                xiyl = xiy(l,j,k,ispec)
                xizl = xiz(l,j,k,ispec)
                etaxl = etax(i,l,k,ispec)
                etayl = etay(i,l,k,ispec)
                etazl = etaz(i,l,k,ispec)
                gammaxl = gammax(i,j,l,ispec)
                gammayl = gammay(i,j,l,ispec)
                gammazl = gammaz(i,j,l,ispec)

                temp1l = temp1l + hprimewgll_xx(l,i) &
                        *( (xixl - nx*(nx*xixl+ny*xiyl+nz*xizl))*temp1_p(l,j,k) &
                          +(xiyl - ny*(nx*xixl+ny*xiyl+nz*xizl))*temp2_p(l,j,k) & 
                          +(xizl - nz*(nx*xixl+ny*xiyl+nz*xizl))*temp3_p(l,j,k) )

                temp2l = temp2l + hprimewgll_yy(l,j)  &
                        *( (etaxl - nx*(nx*etaxl+ny*etayl+nz*etazl))*temp1_p(i,l,k) &
                          +(etayl - ny*(nx*etaxl+ny*etayl+nz*etazl))*temp2_p(i,l,k) & 
                          +(etazl - nz*(nx*etaxl+ny*etayl+nz*etazl))*temp3_p(i,l,k) )

                temp3l = temp3l + hprimewgll_zz(l,k)  &
                        *( (gammaxl - nx*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp1_p(i,j,l) &
                          +(gammayl - ny*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp2_p(i,j,l) & 
                          +(gammazl - nz*(nx*gammaxl+ny*gammayl+nz*gammazl))*temp3_p(i,j,l) )
              enddo
              temp1l = temp1l * wgllwgll_yz(j,k)      
              temp2l = temp2l * wgllwgll_xz(i,k)      
              temp3l = temp3l * wgllwgll_xy(i,j)      
                         
              chi4_dot_dot(i,j,k,ispecPML) = - (temp1l + temp2l + temp3l)

            enddo
          enddo 
        enddo

        ! note: the surface integral expressions would be needed for points on a free surface
        !
        ! BUT at the free surface: potentials are set to zero (zero pressure condition), 
        ! thus the additional surface term contributions would be zeored again.
        
      endif ! ispec_is_acoustic
    endif ! ispec_is_inner
  enddo ! num_PML_ispec
  
end subroutine compute_forces_acoustic_PML

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_abs_boundaries(phase_is_inner,NSPEC_AB,NGLOB_AB,&
                        abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        kappastore,ibool,ispec_is_inner, &
                        rhostore,ispec_is_acoustic,&
                        potential_dot_acoustic,potential_dot_dot_acoustic,&
                        num_PML_ispec,PML_ispec,ispec_is_PML_inum,&
                        chi1_dot,chi2_t,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi3_dot_dot,chi4_dot_dot)

! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)

  use constants,only: NGLLX,NGLLY,NGLLZ,NGLLSQUARE,CUSTOM_REAL
  implicit none

  integer :: NSPEC_AB,NGLOB_AB

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot,chi2_t,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi3_dot_dot,chi4_dot_dot
  integer,dimension(num_PML_ispec):: PML_ispec
  integer,dimension(NSPEC_AB):: ispec_is_PML_inum  

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
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw,temp
  integer :: ispec,iglob,i,j,k,iface,igll,ispecPML
  
  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then
    
      if( ispec_is_acoustic(ispec) .and. ispec_is_PML_inum(ispec) > 0 ) then
      
        do ispecPML=1,num_PML_ispec
        
          if( PML_ispec(ispecPML) == ispec) then

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

              temp = jacobianw / cpl / rhol
              
              ! Sommerfeld condition
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                  - potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
              ! split-potentials
              chi1_dot_dot(i,j,k,ispecPML) = chi1_dot_dot(i,j,k,ispecPML) - chi1_dot(i,j,k,ispecPML) * temp
              chi3_dot_dot(i,j,k,ispecPML) = chi3_dot_dot(i,j,k,ispecPML) - chi3_dot(i,j,k,ispecPML) * temp
              chi4_dot_dot(i,j,k,ispecPML) = chi4_dot_dot(i,j,k,ispecPML) - chi4_dot(i,j,k,ispecPML) * temp
              
              ! chi2 potential?
              chi2_t_dot(i,j,k,ispecPML) = chi2_t_dot(i,j,k,ispecPML) - chi2_t(i,j,k,ispecPML) * temp              
              
            enddo
          endif
        enddo
      endif ! ispec_is_acoustic
    endif ! ispec_is_inner
  enddo ! num_abs_boundary_faces
  
end subroutine PML_acoustic_abs_boundaries

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_interface_coupling(phase_is_inner,NSPEC_AB,NGLOB_AB,&
                        potential_dot_dot_acoustic,&
                        ibool,ispec_is_inner,ispec_is_acoustic,&
                        num_PML_ispec,PML_ispec,iglob_is_PML_interface,&
                        chi1_dot_dot,chi3_dot_dot,chi4_dot_dot)

! couples potential_dot_dot with PML interface contributions

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL 
  implicit none

  integer :: NGLOB_AB,NSPEC_AB
  
  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi3_dot_dot,chi4_dot_dot
  integer,dimension(num_PML_ispec):: PML_ispec
  integer,dimension(NGLOB_AB):: iglob_is_PML_interface
  
  ! potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic    
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  
  !local parameters
  integer :: iglob,ispecPML,i,j,k,ispec

  ! experimental:
  ! updates with the contribution from potential_dot_dot_acoustic on split potentials and vice versa
  
  do ispecPML = 1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)

    if( ispec_is_inner(ispec) .eqv. phase_is_inner ) then
    
      ! acoustic potentials
      if( ispec_is_acoustic(ispec) ) then 
    
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX          
              iglob = ibool(i,j,k,ispec)
              
              ! sums contributions to PML potentials on interface points    
              if( iglob_is_PML_interface(iglob) > 0 ) then   
   
                ! this would be the contribution to the potential_dot_dot array
                ! note: on PML interface, damping coefficient d should to be zero
                !           as well as dprime (-> no chi2 contribution)
                
                potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                              + chi1_dot_dot(i,j,k,ispecPML) &
                                              + chi3_dot_dot(i,j,k,ispecPML) &
                                              + chi4_dot_dot(i,j,k,ispecPML) 

              endif ! interface iglob
            enddo
          enddo
        enddo

      endif ! ispec_is_acoustic      
    endif ! ispec_is_inner    
  enddo ! ispecPML

                        
end subroutine PML_acoustic_interface_coupling


!
!-------------------------------------------------------------------------------------------------
!


subroutine PML_acoustic_mass_update(NSPEC_AB,NGLOB_AB,&
                        ispec_is_acoustic,rmass_acoustic,ibool,&
                        num_PML_ispec,PML_ispec,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot)

! updates split-potentials with local mass in PML region

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL                 
  implicit none  
  integer :: NSPEC_AB,NGLOB_AB

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  integer,dimension(num_PML_ispec):: PML_ispec
  
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: rmass_acoustic
  
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  !local parameters
  real(kind=CUSTOM_REAL):: mass
  integer :: ispec,ispecPML,i,j,k,iglob

  ! updates the dot_dot potentials for the PML
  do ispecPML = 1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)    
    
    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)

            ! global mass ( sum over elements included)
            mass = rmass_acoustic(iglob)
            
            chi1_dot_dot(i,j,k,ispecPML)    = chi1_dot_dot(i,j,k,ispecPML) * mass
            chi2_t_dot_dot(i,j,k,ispecPML)  = chi2_t_dot_dot(i,j,k,ispecPML) * mass
            chi3_dot_dot(i,j,k,ispecPML)    = chi3_dot_dot(i,j,k,ispecPML) * mass
            chi4_dot_dot(i,j,k,ispecPML)    = chi4_dot_dot(i,j,k,ispecPML) * mass
            
          enddo
        enddo
      enddo
    endif
  enddo

end subroutine PML_acoustic_mass_update

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_time_march(NSPEC_AB,NGLOB_AB,ibool,&
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


! time marching scheme - updates with corrector terms
!
! note that the value of d changes according to the distance to the boundary,
! thus instead of updating the whole arrays chi1(:) this scheme updates every given,single value chi1,...

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL  
  implicit none

  integer :: NSPEC_AB,NGLOB_AB
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

  ! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_acoustic  
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_acoustic  
  
  real(kind=CUSTOM_REAL):: deltat,deltatsqover2,deltatover2

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1,chi2,chi2_t,chi3,chi4
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          PML_damping_d

  integer,dimension(num_PML_ispec):: PML_ispec
  integer,dimension(NGLOB_AB) :: iglob_is_PML_interface    
  logical,dimension(NGLOB_AB) :: PML_mask_ibool
  
  ! MPI communication
  integer :: NPROC
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  !local parameters
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: contributions,contributions_dot
  real(kind=CUSTOM_REAL):: d
  integer :: ispec,ispecPML,i,j,k,iglob

  ! updates local points in PML
  allocate(contributions_dot(NGLOB_AB))
  allocate(contributions(NGLOB_AB))  
  contributions_dot(:) = 0._CUSTOM_REAL
  contributions(:) = 0._CUSTOM_REAL

  do ispecPML = 1,num_PML_ispec
    
    ispec = PML_ispec(ispecPML)    

    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX

            ! updates split-potential arrays
            d = PML_damping_d(i,j,k,ispecPML)

            call PML_acoustic_time_march_s(chi1(i,j,k,ispecPML),chi2(i,j,k,ispecPML),&
                      chi2_t(i,j,k,ispecPML),chi3(i,j,k,ispecPML),chi4(i,j,k,ispecPML), &
                      chi1_dot(i,j,k,ispecPML),chi2_t_dot(i,j,k,ispecPML),&
                      chi3_dot(i,j,k,ispecPML),chi4_dot(i,j,k,ispecPML), &
                      chi1_dot_dot(i,j,k,ispecPML),chi2_t_dot_dot(i,j,k,ispecPML),&
                      chi3_dot_dot(i,j,k,ispecPML),chi4_dot_dot(i,j,k,ispecPML), &
                      deltat,deltatsqover2,deltatover2,d)

            ! adds new contributions
            iglob = ibool(i,j,k,ispec)
            if( iglob_is_PML_interface(iglob) > 0 ) then  
                ! on interface points, the time marched global potential from the regular domains applies
                contributions(iglob) = 0._CUSTOM_REAL
                contributions_dot(iglob) = 0._CUSTOM_REAL                
            else
              contributions(iglob) = contributions(iglob) &
                                      + chi1(i,j,k,ispecPML) &
                                      + chi2(i,j,k,ispecPML) &
                                      + chi3(i,j,k,ispecPML) &
                                      + chi4(i,j,k,ispecPML) 

              contributions_dot(iglob) = contributions_dot(iglob) &
                                      + chi1_dot(i,j,k,ispecPML) - d*chi1(i,j,k,ispecPML) &
                                      + chi2_t(i,j,k,ispecPML) - d*chi2(i,j,k,ispecPML) &
                                      + chi3_dot(i,j,k,ispecPML) - d*chi3(i,j,k,ispecPML) &
                                      + chi4_dot(i,j,k,ispecPML) 
            endif
          enddo
        enddo
      enddo
    endif
  enddo

  ! assembles contributions from different MPI processes
  call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,contributions, &
                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                    my_neighbours_ext_mesh)
  call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,contributions_dot, &
                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                    my_neighbours_ext_mesh)

  ! separates contributions from regular domain
  PML_mask_ibool = .false.  

  !do ispec = 1,NSPEC_AB    
  do ispecPML = 1,num_PML_ispec
    
    ispec = PML_ispec(ispecPML)    
      
    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            
            if( PML_mask_ibool(iglob) .eqv. .false. ) then
              ! on interface points, leave contribution from regular domain

              ! inside PML region, split potentials determine the global acoustic potential  
              if( iglob_is_PML_interface(iglob) == 0 ) then  
                potential_acoustic(iglob) = contributions(iglob) 
                potential_dot_acoustic(iglob) = contributions_dot(iglob)     
              endif
                
              PML_mask_ibool(iglob) = .true.
            endif
          enddo
        enddo
      enddo
    endif
  enddo
 
end subroutine PML_acoustic_time_march

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_time_march_s(chi1,chi2,chi2_t,chi3,chi4, &
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot, &
                        chi1_dot_dot,chi2_t_dot_dot, &
                        chi3_dot_dot,chi4_dot_dot, &
                        deltat,deltatsqover2,deltatover2,d)

! time marching scheme
!
! note that the value of d changes according to the distance to the boundary,
! thus instead of updating the whole arrays chi1(:) this scheme updates every given,single value chi1,...
  use constants,only: CUSTOM_REAL
  implicit none
  real(kind=CUSTOM_REAL):: chi1,chi2,chi2_t,chi3,chi4
  real(kind=CUSTOM_REAL):: chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL):: chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL):: deltat,deltatsqover2,deltatover2,d
  !local parameters
  real(kind=CUSTOM_REAL):: fac1,fac2,fac3,fac4
  
  ! pre-computes some factors
  fac1 = 1._CUSTOM_REAL/(1.0_CUSTOM_REAL + deltatover2*d)
  fac2 = 1._CUSTOM_REAL/(d + 1.0_CUSTOM_REAL/deltatover2)
  fac3 = 1._CUSTOM_REAL/(2.0_CUSTOM_REAL + deltat*d)
  fac4 = deltatsqover2*d*d - deltat*d
    
  ! first term: chi1(t+deltat) update
  chi1            = chi1 + deltat*chi1_dot + deltatsqover2*chi1_dot_dot &
                    + fac4*chi1 - deltat*deltat*d*chi1_dot 
                
  ! chi1_dot predictor                      
  chi1_dot        = fac1 * chi1_dot - d*fac2 * chi1_dot + fac2 * chi1_dot_dot
  chi1_dot_dot    = 0._CUSTOM_REAL

  ! second term: chi2  
  ! note that it uses chi2_t at time ( t )  
  chi2            = 2.0*fac3 * chi2 - deltat*d*fac3 * chi2 + deltat*fac3 * chi2_t
            
  ! temporary chi2_t(t+deltat) update  
  chi2_t          = chi2_t + deltat*chi2_t_dot + deltatsqover2*chi2_t_dot_dot &
                    + fac4*chi2_t - deltat*deltat*d*chi2_t_dot
            
  ! chi2 - corrector using updated chi2_t(t+deltat)
  chi2            = chi2 + deltat*fac3 * chi2_t
  
  ! temporary chi2_t_dot - predictor  
  chi2_t_dot      = fac1 * chi2_t_dot - d*fac2 * chi2_t_dot + fac2 * chi2_t_dot_dot
  chi2_t_dot_dot  = 0._CUSTOM_REAL  
  
  ! third term: chi3 (t+deltat) update  
  chi3            = chi3 + deltat*chi3_dot + deltatsqover2*chi3_dot_dot &
                    + fac4*chi3 - deltatsqover2*d*chi3_dot            
  chi3_dot        = chi3_dot + deltatover2*chi3_dot_dot
  chi3_dot_dot    = 0._CUSTOM_REAL
    
  ! fourth term: chi4 (t+deltat) update  
  chi4            = chi4 + deltat*chi4_dot + deltatsqover2*chi4_dot_dot  
  chi4_dot        = chi4_dot + deltatover2*chi4_dot_dot
  chi4_dot_dot    = 0._CUSTOM_REAL
  
end subroutine PML_acoustic_time_march_s


!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_time_corrector(NSPEC_AB,ispec_is_acoustic,deltatover2,&
                        num_PML_ispec,PML_ispec,PML_damping_d,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot)

! time marching scheme - updates with corrector terms
!
! note that the value of d changes according to the distance to the boundary,
! thus instead of updating the whole arrays chi1(:) this scheme updates every given,single value chi1,...

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL 
  implicit none  
  
  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: PML_damping_d

  integer,dimension(num_PML_ispec):: PML_ispec

  real(kind=CUSTOM_REAL):: deltatover2

  integer :: NSPEC_AB
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  !local parameters
  real(kind=CUSTOM_REAL):: d
  integer :: ispec,ispecPML,i,j,k

  ! updates "velocity" potentials in PML with corrector terms
  do ispecPML = 1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)   

    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then 
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
          
            ! time marches chi_dot,.. potentials
            d = PML_damping_d(i,j,k,ispecPML)
              
            call PML_acoustic_time_corrector_s(chi1_dot(i,j,k,ispecPML),chi2_t_dot(i,j,k,ispecPML), &
                      chi3_dot(i,j,k,ispecPML),chi4_dot(i,j,k,ispecPML), &
                      chi1_dot_dot(i,j,k,ispecPML),chi2_t_dot_dot(i,j,k,ispecPML), &
                      chi3_dot_dot(i,j,k,ispecPML),chi4_dot_dot(i,j,k,ispecPML), &
                      deltatover2,d)
          enddo
        enddo
      enddo
    endif
  enddo    

  
end subroutine PML_acoustic_time_corrector

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_time_corrector_s(chi1_dot,chi2_t_dot,chi3_dot,chi4_dot, &
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot, &
                        deltatover2,d)

! time marching scheme - updates with corrector terms
!
! note that the value of d changes according to the distance to the boundary,
! thus instead of updating the whole arrays chi1(:) this scheme updates every given,single value chi1,...
  use constants,only: CUSTOM_REAL
  implicit none
  real(kind=CUSTOM_REAL):: chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL):: chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL):: deltatover2,d
  real(kind=CUSTOM_REAL):: fac1
  
  fac1 = 1.0_CUSTOM_REAL/(d + 1.0_CUSTOM_REAL/deltatover2)

  ! first term:
  chi1_dot = chi1_dot + fac1*chi1_dot_dot
  
  ! second term:
  chi2_t_dot = chi2_t_dot + fac1*chi2_t_dot_dot

  ! third term:
  chi3_dot = chi3_dot + deltatover2*chi3_dot_dot
  
  ! fourth term:
  chi4_dot = chi4_dot + deltatover2*chi4_dot_dot

end subroutine PML_acoustic_time_corrector_s


!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_acoustic_enforce_free_srfc(NSPEC_AB,NGLOB_AB, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces, &
                        ispec_is_acoustic, &
                        num_PML_ispec,PML_ispec,&
                        chi1,chi2,chi2_t,chi3,chi4,&
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot,&
                        chi1_dot_dot,chi2_t_dot_dot,&
                        chi3_dot_dot,chi4_dot_dot)
                      
  use constants,only: NGLLX,NGLLY,NGLLZ,NGLLSQUARE,CUSTOM_REAL
  implicit none 

  integer :: NSPEC_AB,NGLOB_AB

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1,chi2,chi2_t,chi3,chi4
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot
  integer,dimension(num_PML_ispec):: PML_ispec
  
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
  integer :: iface,igll,i,j,k,ispec,iglob,ispecPML

  ! enforce potentials to be zero at surface 
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if( ispec_is_acoustic(ispec) ) then 
      
      do ispecPML=1,num_PML_ispec
        if( PML_ispec(ispecPML) == ispec ) then
      
          do igll = 1, NGLLSQUARE
            i = free_surface_ijk(1,igll,iface)
            j = free_surface_ijk(2,igll,iface)
            k = free_surface_ijk(3,igll,iface)
            iglob = ibool(i,j,k,ispec)

            ! sets potentials to zero
            potential_acoustic(iglob)         = 0._CUSTOM_REAL
            potential_dot_acoustic(iglob)     = 0._CUSTOM_REAL
            potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
            
            ! sets PML potentials to zero 
            chi1(i,j,k,ispecPML) = 0._CUSTOM_REAL  
            chi1_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi1_dot_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            
            chi2(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi2_t(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi2_t_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi2_t_dot_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            
            chi3(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi3_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi3_dot_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            
            chi4(i,j,k,ispecPML) = 0._CUSTOM_REAL  
            chi4_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
            chi4_dot_dot(i,j,k,ispecPML) = 0._CUSTOM_REAL
          enddo
        endif
      enddo
    endif
    
  enddo

end subroutine PML_acoustic_enforce_free_srfc



!
!-------------------------------------------------------------------------------------------------
!


subroutine PML_acoustic_update_potentials(NGLOB_AB,NSPEC_AB, &
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

! updates potential_dot_acoustic and potential_dot_dot_acoustic inside PML region

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  implicit none
  
  integer :: NGLOB_AB,NSPEC_AB

  ! split-potentials
  integer :: num_PML_ispec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1,chi2,chi2_t,chi3
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          chi1_dot_dot,chi3_dot_dot,chi4_dot_dot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,num_PML_ispec):: &
          PML_damping_d
  integer,dimension(num_PML_ispec):: PML_ispec
  integer,dimension(NGLOB_AB):: iglob_is_PML_interface
  logical,dimension(NGLOB_AB):: PML_mask_ibool
  
  
  ! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_acoustic  
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: potential_dot_dot_acoustic  
  
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

  ! MPI communication
  integer :: NPROC
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh

  !local parameters
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: contributions_dot_dot,contributions_dot
  real(kind=CUSTOM_REAL):: d
  integer :: ispec,ispecPML,i,j,k,iglob

  allocate(contributions_dot_dot(NGLOB_AB),contributions_dot(NGLOB_AB))
  contributions_dot_dot = 0._CUSTOM_REAL
  contributions_dot = 0._CUSTOM_REAL

  ! updates the potential_dot & potential_dot_dot_acoustic array inside the PML
  do ispecPML = 1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)    
    
    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then
    
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)

            ! for points inside PML region
            if( iglob_is_PML_interface(iglob) == 0 ) then
              
              ! damping coefficient                
              d = PML_damping_d(i,j,k,ispecPML)

              ! inside PML region: at this stage, this is only needed for seismogram/plotting output
              !                                afterwards potential_dot_dot, resp. chi1_dot_dot,.. get reset to zero

              ! potential_dot: note that we defined 
              !   chi1_dot = (\partial_t + d) chi1 
              !   chi2_t = (\partial_t + d) chi2
              !   chi3_dot = (\partial_t + d) chi3
              !   chi4_dot = \partial_t chi4
              ! where \partial_t is the time derivative, thus \partial_t (chi1+chi2+chi3+chi4) equals
              contributions_dot(iglob) = contributions_dot(iglob) &
                                            + chi1_dot(i,j,k,ispecPML) - d*chi1(i,j,k,ispecPML) &
                                            + chi2_t(i,j,k,ispecPML) - d*chi2(i,j,k,ispecPML) &
                                            + chi3_dot(i,j,k,ispecPML) - d*chi3(i,j,k,ispecPML) &
                                            + chi4_dot(i,j,k,ispecPML)
                            
              ! potential_dot_dot: note that we defined 
              !   chi1_dot_dot = (\partial_t + d)**2 chi1 
              !   chi2_t_dot = (\partial_t + d)**2 chi2
              !   chi3_dot = \partial_t (\partial_t + d) chi3
              !   chi4_dot = \partial_t**2 chi4
              ! where \partial_t is the time derivative, thus \partial_t**2 (chi1+chi2+chi3+chi4) equals  
              contributions_dot_dot(iglob) = contributions_dot_dot(iglob) &
                + chi1_dot_dot(i,j,k,ispecPML) - 2.0*d*chi1_dot(i,j,k,ispecPML) + d*d*chi1(i,j,k,ispecPML) &
                + chi2_t_dot(i,j,k,ispecPML) - 2.0*d*chi2_t(i,j,k,ispecPML) + d*d*chi2(i,j,k,ispecPML) &
                + chi3_dot_dot(i,j,k,ispecPML) - d*chi3_dot(i,j,k,ispecPML) + d*d*chi3(i,j,k,ispecPML) &
                + chi4_dot_dot(i,j,k,ispecPML)
            
            endif
          
          enddo
        enddo
      enddo
    endif
  enddo

  ! assembles contributions from different MPI processes
  call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,contributions_dot, &
                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                    my_neighbours_ext_mesh)
  call assemble_MPI_scalar_ext_mesh(NPROC,NGLOB_AB,contributions_dot_dot, &
                    num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                    my_neighbours_ext_mesh)

  ! updates the potential_dot & potential_dot_dot_acoustic array inside the PML
  PML_mask_ibool = .false.
  do ispecPML = 1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)    
    
    ! acoustic potentials
    if( ispec_is_acoustic(ispec) ) then
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)

            if( PML_mask_ibool(iglob) .eqv. .false. ) then
              ! for points inside PML region
              if( iglob_is_PML_interface(iglob) == 0 ) then
                potential_dot_acoustic(iglob) = contributions_dot(iglob)
                potential_dot_dot_acoustic(iglob) = contributions_dot(iglob)                
              endif
              PML_mask_ibool(iglob) = .true.
            endif
          enddo
        enddo
      enddo
    endif
  enddo

end subroutine PML_acoustic_update_potentials

