!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
!> Read elastic information of the background model, define precomputable
!! matrices for mass, stiffness, boundary terms, pointwise derivatives.
module def_precomp_terms

  use global_parameters
  use data_mesh
  use data_spec
  use data_source, only: src_type
  use data_io, only: verbose, coupling
  use data_proc

  use get_mesh, only: compute_coordinates_mesh
  use utlity
  use analytic_mapping

  implicit none

  public :: read_model_compute_terms
  private
contains

!-----------------------------------------------------------------------------------------
!> Wrapper routine to contain globally defined large matrices that are not
!! used in the time loop to this module (e.g. rho, lambda, mu).
!! Also fills up Q with values (which is used in the time loop)
subroutine read_model_compute_terms

  use get_model
  use attenuation, only: prepare_attenuation
  use commun, only: barrier
  use data_matr, only: Q_mu, Q_kappa, M_w_fl, M0_w_fl, M1chi_fl, M2chi_fl, M4chi_fl, &
                          bdry_matr
  use coupling_mod, only: lambda_cp,mu_cp,rho_cp   !! SB coupling


  real(kind=dp), dimension(:,:,:),allocatable :: rho, lambda, mu, massmat_kwts2
  real(kind=dp), dimension(:,:,:),allocatable :: xi_ani, phi_ani, eta_ani
  real(kind=dp), dimension(:,:,:),allocatable :: fa_ani_theta, fa_ani_phi

  if (lpr .and. verbose > 0) write(*,'(a)') &
            '  ::::::::: BACKGROUND MODEL & PRECOMPUTED MATRICES:::::::'

  if (lpr .and. verbose > 1) write(*,'(a)') '    allocate elastic fields....'

  allocate(rho(0:npol,0:npol,1:nelem),massmat_kwts2(0:npol,0:npol,1:nelem))
  allocate(lambda(0:npol,0:npol,1:nelem),mu(0:npol,0:npol,1:nelem))
  allocate(xi_ani(0:npol,0:npol,1:nelem))
  allocate(phi_ani(0:npol,0:npol,1:nelem))
  allocate(eta_ani(0:npol,0:npol,1:nelem))
  allocate(fa_ani_theta(0:npol,0:npol,1:nelem))
  allocate(fa_ani_phi(0:npol,0:npol,1:nelem))

  !!! SB coupling
  if (coupling) then
     allocate(rho_cp(0:npol,0:npol,1:nelem))
     allocate(lambda_cp(0:npol,0:npol,1:nelem),mu_cp(0:npol,0:npol,1:nelem))
  endif
  !!! SB

  if (anel_true) then
    allocate(Q_mu(1:nel_solid))
    allocate(Q_kappa(1:nel_solid))
  endif

  ! load velocity/density model  (velocities in m/s, density in kg/m^3 )
  if (lpr .and. verbose > 1) write(*,*) '   define background model....'

  if (lpr .and. verbose > 1) write(*,*) '   model is anisotropic....'
  if (anel_true) then
    if (lpr .and. verbose > 1) write(*,*)'   ....and anelastic...'
    call read_model(rho, lambda, mu, xi_ani, phi_ani, eta_ani, fa_ani_theta, &
                        fa_ani_phi, Q_mu, Q_kappa)
  else
    call read_model(rho, lambda, mu, xi_ani, phi_ani, eta_ani, fa_ani_theta, fa_ani_phi)
  endif

  !!! SB
  if (coupling) then
     rho_cp=rho
     lambda_cp=lambda
     mu_cp=mu
  endif
  !!! SB coupling

  if (lpr .and. verbose > 1) write(*,*) '   define mass matrix....'
  call def_mass_matrix_k(rho, lambda, mu, massmat_kwts2)

  if (do_mesh_tests) then
     if (lpr .and. verbose > 1) write(*,*) '   compute mass of the earth model....'
     call compute_mass_earth(rho)
  endif

  if (lpr .and. verbose > 1) write(*,*) &
        '   define precomputed matrices for pointwise derivatives...'
  call compute_pointwisederiv_matrices

  if (do_mesh_tests) then
     if (lpr .and. verbose > 1) &
         write(*,*)'   test pointwise derivatives & Laplacian in solid....'
     call test_pntwsdrvtvs_solid
     if (lpr .and. verbose > 1) &
         write(*,*)'   test pointwise derivatives & Laplacian in fluid....'
     if (have_fluid) call test_pntwsdrvtvs_fluid
  endif

  if (anel_true) then
     if (lpr .and. verbose > 1) write(*, '(/,a,/)') '    preparing ATTENUATION model'
     ! this needs to be done before def_solid_stiffness_terms, as it calculates
     ! the unrelaxed moduli from the ones at reference frequency
     call prepare_attenuation(lambda, mu)
     if (lpr .and. verbose > 1) write(*, '(/,a,/)') '    done preparing ATTENUATION model'
  endif

  if (lpr .and. verbose > 1) write(*,*) '   define solid stiffness terms....'
  call def_solid_stiffness_terms(lambda, mu, massmat_kwts2, xi_ani, phi_ani, &
                                 eta_ani, fa_ani_theta, fa_ani_phi)

  deallocate(xi_ani, phi_ani, eta_ani, fa_ani_theta, fa_ani_phi)
  deallocate(lambda,mu)

  if (have_fluid) then
     if (lpr .and. verbose > 1) write(*,*) '   define fluid stiffness terms....'
     call def_fluid_stiffness_terms(rho, massmat_kwts2)

     if (lpr .and. verbose > 1) write(*,*) '   define solid-fluid boundary terms....'
     call def_solid_fluid_boundary_terms
  endif

  if (lpr .and. verbose > 1) write(*,*) '   ...defined all precomputed arrays'
  deallocate(rho, massmat_kwts2)

  if (lpr .and. verbose > 1) write(*,*) '   ...deallocated unnecessary elastic arrays'

  if (lpr .and. verbose > 0) &
     write(*,*) ' :::::::DONE BACKGROUND MODEL & PRECOMPUTED MATRICES:::::'
  call flush(6)

end subroutine read_model_compute_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_pointwisederiv_matrices
! < The 4 necessary global matrices due to pointwise derivatives d/ds and d/dz:
!! dzdeta/J, dzdxi/J, dsdeta/J, dsdxi/J (J: Jacobian).
!! These are known during the time loop if the strain is computed on-the-fly.
!! This is convenient to avoid recomputing these mapping derivatives at each
!! dumping stage and to avoid knowing the grid itself during the time loop
!! (hence the additional 4 global variables in exchange for at least 2 for the
!! mesh, i.e. only slightly more memory intensive).

  use data_pointwise

  integer          :: iel,inode,ipol,jpol
  real(kind=dp)    :: dsdxi,dzdxi,dsdeta,dzdeta
  real(kind=dp)    :: local_crd_nodes(8,2)

  ! fluid pointwise derivatives
  allocate(DsDeta_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DzDeta_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DsDxi_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(DzDxi_over_J_flu(0:npol,0:npol,1:nel_fluid))
  allocate(inv_s_fluid(0:npol,0:npol,1:nel_fluid))

  ! solid pointwise derivatives
  allocate(DsDeta_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DzDeta_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DsDxi_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(DzDxi_over_J_sol(0:npol,0:npol,1:nel_solid))
  allocate(inv_s_solid(0:npol,0:npol,1:nel_solid))

  ! Solid region
  do iel = 1, nel_solid
    do inode = 1, 8
      call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                                    local_crd_nodes(inode,2),ielsolid(iel),inode)
    enddo
    if (.not. axis_solid(iel)) then ! non-axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDeta_over_J_sol(ipol,jpol,iel) = -dsdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDeta_over_J_sol(ipol,jpol,iel) = dzdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDxi_over_J_sol(ipol,jpol,iel) = dsdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDxi_over_J_sol(ipol,jpol,iel) = -dzdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         inv_s_solid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielsolid(iel))

       enddo
    enddo
   else ! axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                  xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         DsDeta_over_J_sol(ipol,jpol,iel) = -dsdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDeta_over_J_sol(ipol,jpol,iel) = dzdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DsDxi_over_J_sol(ipol,jpol,iel) = dsdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))
         DzDxi_over_J_sol(ipol,jpol,iel) = -dzdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(iel))

         if (ipol > 0) then
            inv_s_solid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielsolid(iel))
         else
            inv_s_solid(ipol,jpol,iel) = one
         endif

        enddo
    enddo
    endif !axis
  enddo

  if (verbose > 1) then
     write(69,*)'Pointwise derivative precomputed terms in solid:'
     write(69,8)'  min/max DsDeta/J [1/m]:',minval(DsDeta_over_J_sol), &
          maxval(DsDeta_over_J_sol)
     write(69,8)'  min/max DzDeta/J [1/m]:',minval(DzDeta_over_J_sol), &
          maxval(DzDeta_over_J_sol)
     write(69,8)'  min/max DsDxi/J  [1/m]:',minval(DsDxi_over_J_sol), &
          maxval(DsDxi_over_J_sol)
     write(69,8)'  min/max DzDxi/J  [1/m]:',minval(DzDxi_over_J_sol), &
          maxval(DzDxi_over_J_sol)
     write(69,*)
  endif

  ! Fluid region
  do iel = 1, nel_fluid
    do inode = 1, 8
      call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                                  local_crd_nodes(inode,2),ielfluid(iel),inode)
    enddo
    if (.not. axis_fluid(iel)) then ! non-axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDeta_over_J_flu(ipol,jpol,iel) = -dsdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDeta_over_J_flu(ipol,jpol,iel) = dzdeta / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDxi_over_J_flu(ipol,jpol,iel) = dsdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDxi_over_J_flu(ipol,jpol,iel) = -dzdxi / &
                  jacobian(eta(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))

         inv_s_fluid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielfluid(iel))
      enddo
    enddo
   else ! axial elements
    do ipol=0,npol
       do jpol=0,npol
         call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                  xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))

         DsDeta_over_J_flu(ipol,jpol,iel) = -dsdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDeta_over_J_flu(ipol,jpol,iel) = dzdeta / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DsDxi_over_J_flu(ipol,jpol,iel) = dsdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))
         DzDxi_over_J_flu(ipol,jpol,iel) = -dzdxi / &
                  jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,ielfluid(iel))

         if (ipol > 0) then
            inv_s_fluid(ipol,jpol,iel) = one/scoord(ipol,jpol,ielfluid(iel))
         else
            inv_s_fluid(ipol,jpol,iel) = one
         endif

       enddo
    enddo
    endif !axis
  enddo

  if (verbose > 1) then
     write(69,*)'Pointwise derivative precomputed terms in fluid:'
     write(69,8)'  min/max DsDeta/J [1/m]:',minval(DsDeta_over_J_flu), &
          maxval(DsDeta_over_J_flu)
     write(69,8)'  min/max DzDeta/J [1/m]:',minval(DzDeta_over_J_flu), &
          maxval(DzDeta_over_J_flu)
     write(69,8)'  min/max DsDxi/J  [1/m]:',minval(DsDxi_over_J_flu), &
          maxval(DsDxi_over_J_flu)
     write(69,8)'  min/max DzDxi/J  [1/m]:',minval(DzDxi_over_J_flu), &
          maxval(DzDxi_over_J_flu)
     write(69,*)
  endif

8 format(a25,2(1pe14.4))

end subroutine compute_pointwisederiv_matrices
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_pntwsdrvtvs_solid
  ! < Test pointwise derivatives & axisymmetric Laplacian in solid region

  use data_io
  use pointwise_derivatives


  real(kind=realkind),allocatable :: tmpsolfield(:,:,:)
  real(kind=realkind),allocatable :: tmpsolfieldcomp(:,:,:,:)
  real(kind=realkind),allocatable :: tmpsolfielddiff(:,:,:,:)
  real(kind=realkind)             :: meandiff(2)
  real(kind=dp)   ,allocatable    :: elderiv(:,:)
  real(kind=dp)                   :: s,z,r,theta
  integer                         :: iel,ipol,jpol,iarr(3)
  character(len=16)               :: fmt1

  allocate(tmpsolfield(0:npol,0:npol,1:nel_solid))
  allocate(tmpsolfieldcomp(0:npol,0:npol,1:nel_solid,1:3))
  allocate(tmpsolfielddiff(0:npol,0:npol,1:nel_solid,1:3))
  allocate(elderiv(0:npol,0:npol))

  ! Test derivatives: define scalar field inside solid
  do iel=1,nel_solid
     do jpol=0,npol
        do ipol=0,npol
           ! make tmpsolfield equal to z s^3 + s z^3
           call compute_coordinates(s,z,r,theta,ielsolid(iel),ipol,jpol)
           tmpsolfield(ipol,jpol,iel)= z*s**3+s*z**3 + (s+z)*router
        enddo
     enddo
  enddo

  call axisym_gradient_solid(tmpsolfield,tmpsolfieldcomp(:,:,:,1:2))

  meandiff = zero

  open(unit=34,file=infopath(1:lfinfo)// '/pointwise_deriv_sol_num.dat'//appmynum)
  open(unit=36,file=infopath(1:lfinfo)// '/pointwise_deriv_reldiff.dat'//appmynum)

  do iel=1,nel_solid
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielsolid(iel),ipol,jpol)

           tmpsolfielddiff(ipol,jpol,iel,1)=absreldiff(&
                                   tmpsolfieldcomp(ipol,jpol,iel,1), &
                                   real(3.d0*z*s**2+z**3+router,kind=realkind))
           tmpsolfielddiff(ipol,jpol,iel,2)=absreldiff(&
                                   tmpsolfieldcomp(ipol,jpol,iel,2), &
                                   real(3.d0*s*z**2+s**3+router,kind=realkind))

           meandiff(1)=meandiff(1)+tmpsolfielddiff(ipol,jpol,iel,1)
           meandiff(2)=meandiff(2)+tmpsolfielddiff(ipol,jpol,iel,2)

           if (ipol == npol/2 .and. jpol == npol/2) then
              write(34,12)s,z,tmpsolfieldcomp(npol/2,npol/2,iel,1), &
                   tmpsolfieldcomp(npol/2,npol/2,iel,2)
              write(36,13)r/1000.,theta*180./pi, &
                          absreldiff(tmpsolfieldcomp(npol/2,npol/2,iel,1), &
                          real(3.d0*z*s**2+z**3+router,kind=realkind)), &
                          absreldiff(tmpsolfieldcomp(npol/2,npol/2,iel,2), &
                          real(3.d0*s*z**2+s**3+router,kind=realkind))
           endif
        enddo
     enddo
  enddo
  close(34)
  close(36)


  if (verbose > 1) then
     write(69,*)
     write(69,*)'/_/_/_/_/_/_/SOLID pointwise deriv with f=zs^3+sz^3/_/_/_/_/_/_/'
     write(69,9)'  mean error df/ds          :',meandiff(1)/ &
                                               real((npol)**2*nel_solid)
     write(69,8)'  min/max error df/ds       :',minval(tmpsolfielddiff(:,:,:,1)), &
                                               maxval(tmpsolfielddiff(:,:,:,1))
     iarr = maxloc(tmpsolfielddiff(:,:,:,1))
     write(69,8)'  r[m],theta max error df/ds:', &
                        rcoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3))), &
                        thetacoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3)))*180./pi
     write(69,7)'  elem type,coarsening ,axis :',eltype(ielsolid(iarr(3))), &
                                                coarsing(ielsolid(iarr(3))), &
                                                axis_solid(iarr(3))

     fmt1 = "(K(1pe13.4))"
     write(fmt1(2:2),'(i1.1)') npol+1
     write(69,*)' Analytical derivatives df/ds across that element (-- > xi):'
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielsolid(iarr(3)),ipol,jpol)
           elderiv(ipol,jpol)=3.d0*z*s**2+z**3+router
        enddo
        write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' Numerical derivatives df/ds across that element (-- > xi):'
     do jpol=0,npol
        write(69,fmt1)(tmpsolfieldcomp(ipol,jpol,iarr(3),1),ipol=0,npol)
     enddo
     write(69,*)

     write(69,9)'  mean error df/dz          :',meandiff(2)/ &
                                                real((npol)**2*nel_solid)
     write(69,8)'  min/max error df/dz       :',minval(tmpsolfielddiff(:,:,:,2)), &
                                                maxval(tmpsolfielddiff(:,:,:,2))
     iarr = maxloc(tmpsolfielddiff(:,:,:,2))
     write(69,8)'  r[m],theta max error df/dz:', &
                         rcoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3))), &
                         thetacoord(iarr(1)-1,iarr(2)-1,ielsolid(iarr(3)))*180./pi
     write(69,7)'  elem type,coarsening,axis :',eltype(ielsolid(iarr(3))), &
                                               coarsing(ielsolid(iarr(3))), &
                                               axis_solid(iarr(3))
     write(69,*)' Analytical derivatives df/dz across that element (-- > xi):'
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielsolid(iarr(3)),ipol,jpol)
                                    elderiv(ipol,jpol)=3.d0*s*z**2+s**3+router
        enddo
        write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' Numerical derivatives df/dz across that element (-- > xi):'
     do jpol=0,npol
        write(69,fmt1)(tmpsolfieldcomp(ipol,jpol,iarr(3),2),ipol=0,npol)
     enddo
     write(69,*)'/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'

     write(69,*)
  endif

7  format(a30,a10,2(l4))
8  format(a30,2(1pe14.4))
9  format(a30,1pe14.4)
12 format(4(1pe16.7))
13 format(2(1pe15.5),2(1pe16.7))

  deallocate(tmpsolfield)
  deallocate(tmpsolfieldcomp)
  deallocate(tmpsolfielddiff)
  deallocate(elderiv)

end subroutine test_pntwsdrvtvs_solid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine test_pntwsdrvtvs_fluid
  ! < Test pointwise derivatives & axisymmetric Laplacian in fluid

  use data_io
  use pointwise_derivatives


  real(kind=realkind),allocatable :: tmpflufield(:,:,:)
  real(kind=realkind),allocatable :: tmpflufieldcomp(:,:,:,:)
  real(kind=realkind),allocatable :: tmpflufielddiff(:,:,:,:)
  real(kind=realkind)             :: meandiff(2)
  real(kind=dp)   ,allocatable    :: elderiv(:,:)
  real(kind=dp)                   :: s,z,r,theta
  integer                         :: iel,ipol,jpol,iarr(3)
  character(len=16)               :: fmt1

  allocate(tmpflufield(0:npol,0:npol,1:nel_fluid))
  allocate(tmpflufieldcomp(0:npol,0:npol,1:nel_fluid,1:3))
  allocate(tmpflufielddiff(0:npol,0:npol,1:nel_fluid,1:3))
  allocate(elderiv(0:npol,0:npol))

  ! Test derivatives: define scalar field inside fluid
  do iel=1,nel_fluid
     do jpol=0,npol
        do ipol=0,npol
           ! make tmpflufield equal to z s^3 + s z^3
           call compute_coordinates(s,z,r,theta,ielfluid(iel),ipol,jpol)
           tmpflufield(ipol,jpol,iel)= z*s**3+s*z**3 + (s+z)*router
        enddo
     enddo
  enddo

  call axisym_gradient_fluid(tmpflufield,tmpflufieldcomp(:,:,:,1:2))

  meandiff = zero

  open(unit=34,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_flu_num.dat'//appmynum)
  open(unit=36,file=infopath(1:lfinfo)//&
       '/pointwise_deriv_reldiff.dat'//appmynum)
  do iel=1,nel_fluid
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielfluid(iel),ipol,jpol)

           tmpflufielddiff(ipol,jpol,iel,1)=absreldiff(&
                                   tmpflufieldcomp(ipol,jpol,iel,1), &
                                   real(3.d0*z*s**2+z**3+router,kind=realkind))
           tmpflufielddiff(ipol,jpol,iel,2)=absreldiff(&
                                   tmpflufieldcomp(ipol,jpol,iel,2), &
                                   real(3.d0*s*z**2+s**3+router,kind=realkind))

           meandiff(1)=meandiff(1)+tmpflufielddiff(ipol,jpol,iel,1)
           meandiff(2)=meandiff(2)+tmpflufielddiff(ipol,jpol,iel,2)

           if (ipol == npol/2 .and. jpol == npol/2) then
              write(34,12)s,z,tmpflufieldcomp(npol/2,npol/2,iel,1), &
                   tmpflufieldcomp(npol/2,npol/2,iel,2)
              write(36,13)r/1000.,theta*180./pi, &
                          absreldiff(tmpflufieldcomp(npol/2,npol/2,iel,1), &
                          real(3.d0*z*s**2+z**3+router,kind=realkind)), &
                          absreldiff(tmpflufieldcomp(npol/2,npol/2,iel,2), &
                          real(3.d0*s*z**2+s**3+router,kind=realkind))
           endif
        enddo
     enddo
  enddo
  close(34)
  close(36)

  if (verbose > 1) then
     write(69,*)
     write(69,*)'/_/_/_/_/_/_/FLUID pointwise deriv with f=zs^3+sz^3/_/_/_/_/_/_/'
     write(69,9)'  mean error df/ds          :',meandiff(1)/ &
                                               real((npol)**2*nel_fluid)
     write(69,8)'  min/max error df/ds       :',minval(tmpflufielddiff(:,:,:,1)), &
                                               maxval(tmpflufielddiff(:,:,:,1))
     iarr = maxloc(tmpflufielddiff(:,:,:,1))
     write(69,8)'  r[m],theta max error df/ds:', &
                        rcoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3))), &
                        thetacoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3)))*180./pi
     write(69,7)'  elem type,coarsening ,axis :',eltype(ielfluid(iarr(3))), &
                                                coarsing(ielfluid(iarr(3))), &
                                                axis_fluid(iarr(3))

     fmt1 = "(K(1pe13.4))"
     write(fmt1(2:2),'(i1.1)') npol+1
     write(69,*)' Analytical derivatives df/ds across that element (-- > xi):'
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielfluid(iarr(3)),ipol,jpol)
           elderiv(ipol,jpol)=3.d0*z*s**2+z**3+router
        enddo
        write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' Numerical derivatives df/ds across that element (-- > xi):'
     do jpol=0,npol
        write(69,fmt1)(tmpflufieldcomp(ipol,jpol,iarr(3),1),ipol=0,npol)
     enddo
     write(69,*)

     write(69,9)'  mean error df/dz          :',meandiff(2)/ &
                                                real((npol)**2*nel_fluid)
     write(69,8)'  min/max error df/dz       :',minval(tmpflufielddiff(:,:,:,2)), &
                                                maxval(tmpflufielddiff(:,:,:,2))
     iarr = maxloc(tmpflufielddiff(:,:,:,2))
     write(69,8)'  r[m],theta max error df/dz:', &
                         rcoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3))), &
                         thetacoord(iarr(1)-1,iarr(2)-1,ielfluid(iarr(3)))*180./pi
     write(69,7)'  elem type,coarsening,axis :',eltype(ielfluid(iarr(3))), &
                                               coarsing(ielfluid(iarr(3))), &
                                               axis_fluid(iarr(3))
     write(69,*)' Analytical derivatives df/dz across that element (-- > xi):'
     do jpol=0,npol
        do ipol=0,npol
           call compute_coordinates(s,z,r,theta,ielfluid(iarr(3)),ipol,jpol)
                                    elderiv(ipol,jpol)=3.d0*s*z**2+s**3+router
        enddo
        write(69,fmt1)(elderiv(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' Numerical derivatives df/dz across that element (-- > xi):'
     do jpol=0,npol
        write(69,fmt1)(tmpflufieldcomp(ipol,jpol,iarr(3),2),ipol=0,npol)
     enddo
     write(69,*)'/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/'

     write(69,*)
  endif

7  format(a30,a10,2(l4))
8  format(a30,2(1pe14.4))
9  format(a30,1pe14.4)
12 format(4(1pe16.7))
13 format(2(1pe15.5),2(1pe16.7))

  deallocate(tmpflufield)
  deallocate(tmpflufieldcomp)
  deallocate(tmpflufielddiff)
  deallocate(elderiv)

end subroutine test_pntwsdrvtvs_fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_mass_matrix_k(rho,lambda,mu,massmat_kwts2)
! < This routine computes and stores the coefficients of the diagonal
!! mass matrix, when a weighted Gauss-Lobatto quadrature for axial elements.
!! It is built here with a factor of volume equal to s * ds * dz, as required
!! by our approach.  we also define in this routine the mass matrix weighted
!! by 1/s^2, as required by some components of the Laplacian of the
!! fields. Note the special contribution arising in the case of an
!! axial element.
!! massmat_k    : The actual mass term:    \sigma_I \sigma_J J^{IJ} s^{IJ}
!! massmat_kwts2: Paper 2, table 1, term A=\sigma_I \sigma_J J^{IJ} / s^{IJ}
!! jacob: just defined locally for check on extrema further below....

  use data_io, only: need_fluid_displ, dump_energy
  use commun, only: pdistsum_solid_1D, pdistsum_fluid
  use data_pointwise, only: inv_rho_fluid
  use data_matr, only: set_mass_matrices, unassem_mass_rho_solid, &
                               unassem_mass_lam_fluid

  use data_mesh, only: npol, nelem, nel_solid, nel_fluid

  real(kind=dp), dimension(0:,0:,:),intent(in)  :: rho, lambda, mu
  real(kind=dp), dimension(0:npol,0:npol,nelem),intent(out) :: massmat_kwts2

  real(kind=realkind), allocatable :: inv_mass_rho(:,:,:), inv_mass_fluid(:,:,:)
  real(kind=dp), allocatable       :: massmat_k(:,:,:)   ! < Mass matrix
  real(kind=dp), allocatable       :: jacob (:,:,:)      ! < jacobian array
  real(kind=realkind), allocatable :: drdxi(:,:,:,:)     ! < min/max derivs

  real(kind=dp)     :: local_crd_nodes(8,2),s,z,r,theta
  integer           :: iel,inode,iarr(3),ipol,jpol
  character(len=16) :: fmt1
  real(kind=dp)     :: dsdxi,dzdeta,dzdxi,dsdeta

  allocate(massmat_k(0:npol,0:npol,1:nelem),jacob(0:npol,0:npol,1:nelem))
  allocate(drdxi(0:npol,0:npol,1:nelem,1:4))
  allocate(inv_rho_fluid(0:npol,0:npol,1:nel_fluid))

  allocate(inv_mass_rho(0:npol,0:npol,1:nel_solid))
  allocate(inv_mass_fluid(0:npol,0:npol,1:nel_fluid))

  massmat_k(:,:,:) = zero; massmat_kwts2(:,:,:) = zero

  do iel = 1, nelem
     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1), &
             local_crd_nodes(inode,2),iel,inode)
     enddo

     ! computing global arrays for the respective coordinate mapping derivatives
     ! only needed for min/max write statements below!
     if ( axis(iel) ) then
        do ipol  = 0, npol
           do jpol = 0, npol
              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   xi_k(ipol),eta(jpol),local_crd_nodes,iel)
              drdxi(ipol,jpol,iel,1)=dsdxi; drdxi(ipol,jpol,iel,2)=dzdxi
              drdxi(ipol,jpol,iel,3)=dsdeta; drdxi(ipol,jpol,iel,4)=dzdeta
           enddo
        enddo
     else
        do ipol  = 0, npol
           do jpol = 0, npol
              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,iel)
              drdxi(ipol,jpol,iel,1)=dsdxi; drdxi(ipol,jpol,iel,2)=dzdxi
              drdxi(ipol,jpol,iel,3)=dsdeta; drdxi(ipol,jpol,iel,4)=dzdeta
           enddo
        enddo
     endif

     ! ::::::::::::::::non-axial elements::::::::::::::::
     if (.not. axis(iel)) then
        do ipol  = 0, npol
           do jpol = 0, npol
              massmat_k(ipol,jpol,iel) = &
                   jacobian(eta(ipol),eta(jpol),local_crd_nodes,iel) &
                   *scoord(ipol,jpol,iel)*wt(ipol)*wt(jpol)
              massmat_kwts2(ipol,jpol,iel) =  &
                   jacobian(eta(ipol),eta(jpol),local_crd_nodes,iel) &
                   *scoord(ipol,jpol,iel)**(-1)*wt(ipol)*wt(jpol)
              jacob(ipol,jpol,iel) = jacobian(eta(ipol),eta(jpol), &
                                              local_crd_nodes,iel)
           enddo
        enddo

     ! ::::::::::::::::axial elements::::::::::::::::
     else if (axis(iel)) then
        do ipol  = 0, npol ! Be careful here !!!!
           do jpol = 0, npol
              massmat_k(ipol,jpol,iel) = &
                    jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,iel) &
                    *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                    local_crd_nodes,iel)*wt_axial_k(ipol)*wt(jpol)
              jacob(ipol,jpol,iel) = jacobian(xi_k(ipol),eta(jpol), &
                                              local_crd_nodes,iel)
           enddo
        enddo
        do ipol = 1, npol
           do jpol = 0, npol
              massmat_kwts2(ipol,jpol,iel) = &
                   jacobian(xi_k(ipol),eta(jpol),local_crd_nodes,iel) / &
                   ( scoord(ipol,jpol,iel) * ( one+xi_k(ipol) ) ) * &
                   wt_axial_k(ipol)*wt(jpol)
           enddo
        enddo
        do jpol = 0, npol
           massmat_kwts2(0,jpol,iel) = &
                jacobian(xi_k(0),eta(jpol),local_crd_nodes,iel)* &
                (s_over_oneplusxi_axis(xi_k(0),eta(jpol), &
                local_crd_nodes,iel))**(-1) * &
                wt_axial_k(0)*wt(jpol)
        enddo

        ! check consistency of s/(1+xi) definitions...
        ! axial element but off-axis
        do ipol = 1, npol
           do jpol = 0, npol
              if (.not. dblreldiff_small( scoord(ipol,jpol,iel)/ &
                                        ( one+xi_k(ipol) ), &
                                       s_over_oneplusxi_axis(xi_k(ipol), &
                                        eta(jpol),local_crd_nodes,iel)) ) then
                 write(*,*)procstrg,'PROBLEM: 2 definitions of s/(1+xi) differ'
                 write(*,*)procstrg,'scoord/(1+xi)=',scoord(ipol,jpol,iel)/ &
                                                     ( one+xi_k(ipol) )
                 write(*,*)procstrg,'s_over_onexi=', &
                                     s_over_oneplusxi_axis(xi_k(ipol), &
                                     eta(jpol),local_crd_nodes,iel)
                 write(*,*)procstrg,'iel,ipol,jpol:',iel,ipol,jpol
                 write(*,*)procstrg,'s,r,theta:',scoord(ipol,jpol,iel), &
                           rcoord(ipol,jpol,iel),thetacoord(ipol,jpol,iel)
                 stop
              endif
           enddo
        enddo

     endif ! axial?
  enddo ! iel

  if (lpr .and. verbose > 1) write(*,*) '   solid mass matrix...'
  ! Solid inverse mass term
  do iel=1,nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
              inv_mass_rho(ipol,jpol,iel)= rho(ipol,jpol,ielsolid(iel))*&
                                           massmat_k(ipol,jpol,ielsolid(iel))
        enddo
     enddo
  enddo

  ! unassembled mass matrix in solid for energy.
  if (dump_energy) then
    allocate(unassem_mass_rho_solid(0:npol,0:npol,nel_solid))
    unassem_mass_rho_solid = inv_mass_rho
    if (src_type(1) == 'dipole') &
       unassem_mass_rho_solid = two * unassem_mass_rho_solid
  endif


  ! Exchange boundary information
  if (lpr .and. verbose > 1) write(*,*) '   assemble solid mass matrix...'
  call pdistsum_solid_1D(inv_mass_rho)

  if (lpr .and. verbose > 1) write(*,*) '   compute inverse solid mass matrix...'
  do iel=1,nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           if ( inv_mass_rho(ipol,jpol,iel) /= zero) then
              inv_mass_rho(ipol,jpol,iel) = one / inv_mass_rho(ipol,jpol,iel)
           else
              write(*,*)procstrg,'WARNING: solid mass term zero!', &
                         ipol,jpol,iel,ielsolid(iel)
              inv_mass_rho(ipol,jpol,iel) = one / rho(ipol,jpol,ielsolid(iel))
           endif
        enddo
     enddo
  enddo

  if (src_type(1) == 'dipole') inv_mass_rho = half * inv_mass_rho

  ! Fluid inverse mass term
  if (lpr .and. verbose > 1) write(*,*) '   fluid mass matrix...'
  do iel=1,nel_fluid
     ! check if fluid element is really fluid throughout
     if (maxval(mu(:,:,ielfluid(iel))) > zero) then
        call compute_coordinates(s,z,r,theta,ielfluid(iel), &
                                 int(npol/2),int(npol/2))
        fmt1 = '(A,A,4(1PE11.7))'
        write(*,*)
        write(*,*)procstrg,'!!!!!!!!!!!!!!!!  PROBLEM  !!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)procstrg,'Have a non-zero mu in the fluid region!'
        write(*,*)procstrg,'Element: ', ielfluid(iel)
        write(*,fmt1)procstrg,' Mu, Vs, location r[km], theta[deg]:', &
                  maxval(mu(:,:,ielfluid(iel))), &
                  maxval(sqrt(mu(:,:,ielfluid(iel))/rho(:,:,ielfluid(iel)))), &
                  r/1000., theta*180./pi
        call flush(6)
        stop
     endif

     do ipol = 0, npol
        do jpol = 0, npol
           ! since mu=0 only consider lambda here for the
           ! bulk modulus/incompressibility kappa = lambda + 2/3 mu
           ! in the ideal fluid that harbors the outer core...
           inv_mass_fluid(ipol,jpol,iel)=massmat_k(ipol,jpol,ielfluid(iel))/&
                                         lambda(ipol,jpol,ielfluid(iel))

           ! inverse density inside the fluid, needed to calculate displacement in fluid
           if (need_fluid_displ) &
              inv_rho_fluid(ipol,jpol,iel) = one / rho(ipol,jpol,ielfluid(iel))
        enddo
     enddo
  enddo


  ! unassembled mass matrix in fluid for the energy.
  if (dump_energy) then
    allocate(unassem_mass_lam_fluid(0:npol,0:npol,nel_fluid))
    unassem_mass_lam_fluid = inv_mass_fluid
  endif


  ! Exchange boundary information
  if (lpr .and. verbose > 1) write(*,*) '   assemble fluid mass matrix...'
  call pdistsum_fluid(inv_mass_fluid)

  if (lpr .and. verbose > 1) write(*,*) '   compute inverse fluid mass matrix...'
  do iel=1,nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           if (inv_mass_fluid(ipol,jpol,iel) /= zero) then
              inv_mass_fluid(ipol,jpol,iel) = one / inv_mass_fluid(ipol,jpol,iel)
           else
              write(*,*)procstrg,'WARNING: Fluid mass term zero!', &
                        ipol,jpol,iel,ielfluid(iel); call flush(6)
              inv_mass_fluid(ipol,jpol,iel)=lambda(ipol,jpol,ielfluid(iel))
           endif
        enddo
     enddo
  enddo


  ! Call routine in data_matr to actually set the values
  call set_mass_matrices(npol, nel_solid, nel_fluid, inv_mass_rho, inv_mass_fluid)

  ! In the remainder: document min/max values and locations for velocities,
  !    density, Jacobian, mass terms, GLL points. Lagrange derivatives etc.
  if (verbose > 1) then
     fmt1 = "(a18,K(f9.4))"
     write(fmt1(6:6),'(i1.1)') npol+1
     write(69,*)
     write(69,*)'-+-+-+-+-+-+-+-+-+-+ Integration weights +-+-+-+-+-+-+-+-+-+-+-+'
     write(69,fmt1)'GLJ sigma_I ax   :',(wt_axial_k(ipol),ipol=0,npol)
     write(69,fmt1)'GLL sigma_J nonax:',(wt(ipol),ipol=0,npol)

     fmt1 = "(a11,K(f9.4))"
     write(fmt1(6:6),'(i1.1)') npol+1
     write(69,*)
     write(69,*)'-+-+-+- Gauss-Lobatto-Legendre pts, Lagrange derivatives +-+-+-+'
     write(69,fmt1)'GLL eta  :',(eta(ipol),ipol=0,npol)
     write(69,*)
     write(69,*)' Lagrange derivatives eta: \partial_\eta (l_j(\eta_q)), -- > li_i'
     do jpol=0,npol
        write(69,fmt1)'',(G2(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' transpose derivatives eta: -- > eta_p'
     do jpol=0,npol
        write(69,fmt1)'',(G2T(ipol,jpol),ipol=0,npol)
     enddo

     write(69,*)
     write(69,*)'-+-+-+ Gauss-Lobatto-Jacobi (0,1) pts, Lagrange derivatives +-+-'
     write(69,fmt1)'GLJ xi_ax:',(xi_k(ipol),ipol=0,npol)
     write(69,*)
     write(69,*)' Lagrange derivatives G1 xi : \partial_\xi (l_i(\xi_p)), -- > li_i'
     do jpol=0,npol
        write(69,fmt1)'',(G1(ipol,jpol),ipol=0,npol)
     enddo
     write(69,*)' transpose G1T derivatives xi -- > xi_p'
     do jpol=0,npol
        write(69,fmt1)'',(G1T(ipol,jpol),ipol=0,npol)
     enddo

     write(69,*)' Lagrange derivatives axial vector: \partial_\xi (l_i(\xi_0)) '
     write(69,fmt1)'',(G0(ipol),ipol=0,npol)
     write(69,*)' '

     write(69,*)'+-+-+-+-+-+- Coordinate mapping/derivative extrema +-+-+-+-+-+--'
     iarr = minloc(drdxi(:,:,:,1))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th dsdxi [m,m,deg]        :',minval(drdxi(:,:,:,1)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(drdxi(:,:,:,1))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th dsdxi [m,m,deg]        :',maxval(drdxi(:,:,:,1)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     iarr = minloc(drdxi(:,:,:,2))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th dzdxi [m,m,deg]        :',minval(drdxi(:,:,:,2)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(drdxi(:,:,:,2))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th dzdxi [m,m,deg]        :',maxval(drdxi(:,:,:,2)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     iarr = minloc(drdxi(:,:,:,3))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th dsdeta [m,m,deg]       :',minval(drdxi(:,:,:,3)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(drdxi(:,:,:,3))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th dsdeta [m,m,deg]       :',maxval(drdxi(:,:,:,3)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     iarr = minloc(drdxi(:,:,:,4))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th dzdeta [m,m,deg]       :',minval(drdxi(:,:,:,4)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(drdxi(:,:,:,4))
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th dzdeta [m,m,deg]       :',maxval(drdxi(:,:,:,4)), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     write(69,*)' '
     write(69,*)'+-+-+-+-+-+-+-+- Jacobian & mass term extrema -+-+-+-+-+-+-+-+-+'
     write(69,*)'   Jacobian  = dsdxi dzdeta - dsdeta dzdxi'
     write(69,*)'   (mass term)_{IJ} = Jacobian_{IJ} \sigma_I \sigma_J s_{IJ}))'
     iarr = minloc(jacob)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th Jacobian [m,m,deg]        :',minval(jacob), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(jacob)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th Jacobian [m,m,deg]        :',maxval(jacob), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     iarr = minloc(massmat_k)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th mass term [m^3,m,deg]     :',minval(massmat_k), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta

     iarr = maxloc(massmat_k)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th mass term [m^3,m,deg]     :',maxval(massmat_k), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = minloc(inv_mass_rho)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Min,r,th sol. invmass [1/kg,m,deg] :',minval(inv_mass_rho), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     iarr = maxloc(inv_mass_rho)
     theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
     write(69,9)'Max,r,th sol. invmass [1/kg,m,deg] :',maxval(inv_mass_rho), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     if (have_fluid) then
        iarr = minloc(inv_mass_fluid)
        theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
        write(69,9)'Min,r,th flu. invmass [N/m^5,m,deg]:',minval(inv_mass_fluid), &
             rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
        iarr = maxloc(inv_mass_fluid)
        theta = thetacoord(iarr(1)-1,iarr(2)-1,iarr(3)); theta = theta*180./pi
        write(69,9)'Max,r,th flu. invmass [N/m^5,m,deg]:',maxval(inv_mass_fluid), &
                                   rcoord(iarr(1)-1,iarr(2)-1,iarr(3)),theta
     endif


     write(69,*)' '
     write(69,*)'+-+-+-+-+-+-+-+-+-+-+- Elastic extrema -+-+-+-+-+-+-+-+-+-+-+-+-'
     iarr = minloc(sqrt((lambda+2*mu)/rho))
     write(69,8)'Min,r P-vel [m/s,m]   :',minval(dsqrt((lambda+2*mu)/rho)), &
                                        rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
     iarr = minloc(sqrt(mu/rho))
     write(69,8)'Min,r S-vel [m/s,m]   :',minval(dsqrt(mu/rho)), &
                                         rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
     iarr = maxloc(sqrt((lambda+2*mu)/rho))
     write(69,8)'Max,r P-vel [m/s,m]   :',maxval(dsqrt((lambda+2*mu)/rho)), &
                                         rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
     iarr = maxloc(sqrt(mu/rho))
     write(69,8)'Max,r S-vel [m/s,m]   :',maxval(dsqrt(mu/rho)), &
                                         rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
     iarr = minloc(rho)
     write(69,8)'Min,r rho [kg/m^3,m]  :',minval(rho), &
                                         rcoord(iarr(1)-1,iarr(2)-1,iarr(3))
     iarr = maxloc(rho)
     write(69,8)'Max,r rho [kg/m^3,m]  :',maxval(rho), &
                                          rcoord(iarr(1)-1,iarr(2)-1,iarr(3))/1.d3
     write(69,8)'Min/max mu [N/m^2]    :',minval(mu),maxval(mu)
     write(69,8)'Min/max lambda [N/m^2]:',minval(lambda),maxval(lambda)

     write(69,*)''
  endif

8 format(a25,2(1pe12.4))
9 format(a37,3(1pe12.4))

  deallocate(massmat_k,jacob)
  deallocate(drdxi)

end subroutine def_mass_matrix_k
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_mass_earth(rho)
! < A straight computation of the mass of the sphere and that of its
!! solid and fluid sub-volumes. This is the same as computing the volume
!! (see def_grid.f90), but with a multiplicative density factor, i.e. the
!! actual mass matrix. The comparison to the exact value is merely
!! done as an indicator and does not cause any error. One should keep an eye
!! on these values when generating any new kind of background model for which
!! one then needs to dig up the total mass...

  use def_grid, only: massmatrix
  use background_models, only: velocity
  use commun, only: psum
  use data_io, only: infopath,lfinfo

  real(kind=dp)   , intent(in)     :: rho(0:npol,0:npol,nelem)
  integer                          :: iel,ipol,jpol,idom,iidom,idisc
  integer                          :: count_solid,count_fluid,count_sic
  real(kind=dp)                    :: mass_glob,mass_solid,mass_fluid,mass_sic
  real(kind=dp)                    :: mass_glob_num,mass_solid_num
  real(kind=dp)                    :: mass_fluid_num,mass_sic_num
  real(kind=dp)                    :: mass_layer(ndisc)
  real(kind=dp)                    :: r,density,vs,dr
  real(kind=realkind), allocatable :: massmat(:,:,:),massmat_solid(:,:,:)
  real(kind=realkind), allocatable :: massmat_fluid(:,:,:)
  character(len=100)               :: modelstring

  allocate(massmat(0:npol,0:npol,1:nelem))
  allocate(massmat_solid(0:npol,0:npol,1:nel_solid))
  allocate(massmat_fluid(0:npol,0:npol,1:nel_fluid))

  mass_fluid = zero
  mass_glob = zero
  mass_solid = zero
  mass_sic = zero
  count_fluid = 0
  count_solid = 0
  count_sic = 0

  modelstring=bkgrdmodel!(1:(lfbkgrdmodel))
  ! actual masses, calculated by summation over layers spaced 100 meters
  do iel=1,int(router/100.)+1
     r = router-(real(iel)-.5)*100.
     dr = (router-(real(iel)-1.)*100.)**3 - (router-(real(iel)*100.))**3
     idom=10*ndisc
     do iidom=1,ndisc-1
        if (r <= discont(iidom) .and. r > discont(iidom+1) ) then
           idom=iidom
           exit
        endif
        if (idom == iidom) exit
     enddo
     if (r <= discont(ndisc)) idom=ndisc
     density=velocity(r,'rho',idom,modelstring,lfbkgrdmodel)
     mass_glob = mass_glob +  density*dr
     vs=velocity(r,'v_s',idom,modelstring,lfbkgrdmodel)
     if (vs < 0.1d0) then
        mass_fluid = mass_fluid + density*dr
     else
        mass_solid = mass_solid + density*dr
     endif
     ! Solid inner core only
     if (idom == ndisc) mass_sic = mass_sic + density*dr
  enddo
  mass_glob  = 4.d0/3.d0*pi*mass_glob
  mass_fluid = 4.d0/3.d0*pi*mass_fluid
  mass_solid = 4.d0/3.d0*pi*mass_solid
  mass_sic = 4.d0/3.d0*pi*mass_sic

  mass_glob_num=zero; mass_solid_num=zero; mass_fluid_num=zero
  mass_sic_num=zero
  mass_layer = 0

  ! numerically computed masses
  call massmatrix(massmat,nelem,'total')
  call massmatrix(massmat_solid,nel_solid,'solid')
  call massmatrix(massmat_fluid,nel_fluid,'fluid')

  do iel = 1, nelem
     do ipol = 0, npol
        do jpol = 0, npol
           mass_glob_num = mass_glob_num + &
                           rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)
           ! solid inner core only
           if (rcoord(npol/2,npol/2,iel) < discont(ndisc)) &
                mass_sic_num = mass_sic_num + &
                               rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)

           ! compute mass for each layer between discontinuities
           do idisc = 1,ndisc-1
              if (rcoord(npol/2,npol/2,iel) <= discont(idisc) .and. &
                  rcoord(npol/2,npol/2,iel) > discont(idisc+1) ) then
                 mass_layer(idisc) = mass_layer(idisc) + &
                                     rho(ipol,jpol,iel)*massmat(ipol,jpol,iel)
              endif
           enddo

        enddo
     enddo
  enddo
  mass_layer(ndisc) = mass_sic_num
  mass_glob_num=2.d0*pi*mass_glob_num
  mass_glob_num=psum(real(mass_glob_num,kind=realkind))
  mass_sic_num=2.d0*pi*mass_sic_num
  mass_sic_num=psum(real(mass_sic_num,kind=realkind))
  mass_layer=2.d0*pi*mass_layer
  do idisc = 1,ndisc-1
     mass_layer(idisc) = psum(real(mass_layer(idisc),kind=realkind))
  enddo

  do iel = 1, nel_solid
     do ipol = 0, npol
        do jpol = 0, npol
           mass_solid_num = mass_solid_num + rho(ipol,jpol,ielsolid(iel))* &
                                             massmat_solid(ipol,jpol,iel)
        enddo
     enddo
  enddo
  mass_solid_num=2.d0*pi*mass_solid_num
  mass_solid_num=psum(real(mass_solid_num,kind=realkind))

  do iel = 1, nel_fluid
     do ipol = 0, npol
        do jpol = 0, npol
           mass_fluid_num = mass_fluid_num + rho(ipol,jpol,ielfluid(iel))* &
                                             massmat_fluid(ipol,jpol,iel)
        enddo
     enddo
  enddo
  mass_fluid_num=2.d0*pi*mass_fluid_num
  mass_fluid_num=psum(real(mass_fluid_num,kind=realkind))

  if (lpr) then
     write(*,*)'  Calculated masses for earth model: ', &
                                          bkgrdmodel(1:lfbkgrdmodel)
     write(*,10)'  Total mass (real,num,diff)    :',mass_glob,mass_glob_num, &
                                         abs(mass_glob-mass_glob_num)/mass_glob
     write(*,11)'  Sum of layers (real,num,diff) :',sum(mass_layer)
     write(*,10)'    Solid mass (real,num,diff)    :',mass_solid, &
                       mass_solid_num,abs(mass_solid-mass_solid_num)/mass_solid
     if (have_fluid) then
        write(*,10)'    Fluid mass (real,num,diff)    :',mass_fluid, &
                    mass_fluid_num,abs(mass_fluid-mass_fluid_num)/mass_fluid
     endif
     write(*,10)'    Innercore mass (real,num,diff):',mass_sic, &
                       mass_sic_num,abs(mass_sic-mass_sic_num)/mass_sic
     write(*,*) '   Mass of the Earth             :   5.97 x 10^24 kg'
10   format(a35,2(1pe14.5),1pe12.2)
11   format(a35,1pe14.5)
     write(*,*)
  endif

! write out total masses for each layer
  if (lpr) then
     open(unit=9876,file=infopath(1:lfinfo)//'/mass_kg_per_discont_layer.dat')
     do idisc=1,ndisc
        write(9876,*)discont(idisc),mass_layer(idisc)
     enddo
     close(9876)
  endif

  deallocate(massmat)
  deallocate(massmat_solid)
  deallocate(massmat_fluid)

end subroutine compute_mass_earth
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_solid_stiffness_terms(lambda, mu, massmat_kwts2, xi_ani, phi_ani, &
                                     eta_ani, fa_ani_theta, fa_ani_phi)
! < This routine is a merged version to minimize global work
!! array definitions. The terms alpha_wt_k etc. are now
!! merely elemental arrays, and defined on the fly when
!! computing the global, final precomputable matrices
!! for the solid stiffness term. The loop is over solid elements only.
!! Adding optional arguments for anisotropy, MvD

  use attenuation, only: att_coarse_grained
  use data_mesh, only: npol, nel_solid, nelem
  use data_matr, only: M11s, M21s, M41s, M12s, M22s, M42s, M11z, M21z, M41z, M32s, &
                       M13s, M33s, M43s, &
                       M_1, M_2, M_3, M_4, M_5, M_6, M_7, M_8, &
                       M_w1, M_w2, M_w3, M_w4, M_w5, &
                       M0_w1, M0_w2, M0_w3, M0_w4, M0_w5, M0_w6, M0_w7, M0_w8, M0_w9, M0_w10, &
                       M1phi, M2phi, M4phi, &
                       Y_cg4,  V_s_eta_cg4, V_s_xi_cg4, V_z_eta_cg4, V_z_xi_cg4, &
                       Y,  V_s_eta, V_s_xi, V_z_eta, V_z_xi, &
                       Y0, V0_s_eta, V0_s_xi, V0_z_eta, V0_z_xi

  real(kind=dp), dimension(0:,0:,:), intent(in) :: lambda,mu
  real(kind=dp), dimension(0:,0:,:), intent(in) :: massmat_kwts2
  real(kind=dp), dimension(0:,0:,:), intent(in), optional :: xi_ani, phi_ani, eta_ani, &
                                                             fa_ani_theta, fa_ani_phi

  real(kind=dp) :: local_crd_nodes(8,2)
  integer       :: ielem, ipol, jpol, inode
  real(kind=dp) :: dsdxi, dzdeta, dzdxi, dsdeta

  real(kind=dp) :: alpha_wt_k(0:npol,0:npol)
  real(kind=dp) :: beta_wt_k(0:npol,0:npol)
  real(kind=dp) :: gamma_wt_k(0:npol,0:npol)
  real(kind=dp) :: delta_wt_k(0:npol,0:npol)
  real(kind=dp) :: epsil_wt_k(0:npol,0:npol)
  real(kind=dp) :: zeta_wt_k(0:npol,0:npol)

  real(kind=dp) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp) :: M_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp) :: M_z_xi_wt_k(0:npol,0:npol)
  real(kind=dp) :: M_z_eta_wt_k(0:npol,0:npol)
  real(kind=dp) :: M_s_eta_wt_k(0:npol,0:npol)

  ! non-diagfac
  real(kind=dp), allocatable :: non_diag_fact(:,:)

  ! Allocate Global Stiffness Arrays depending on source type:
  allocate(M11s(0:npol, 0:npol, nel_solid))
  allocate(M21s(0:npol, 0:npol, nel_solid))
  allocate(M41s(0:npol, 0:npol, nel_solid))
  allocate(M12s(0:npol, 0:npol, nel_solid))
  allocate(M22s(0:npol, 0:npol, nel_solid))
  allocate(M42s(0:npol, 0:npol, nel_solid))
  allocate(M11z(0:npol, 0:npol, nel_solid))
  allocate(M21z(0:npol, 0:npol, nel_solid))
  allocate(M41z(0:npol, 0:npol, nel_solid))
  allocate(M32s(0:npol, 0:npol, nel_solid))

  allocate(M_1(0:npol,0:npol,nel_solid))
  allocate(M_2(0:npol,0:npol,nel_solid))
  allocate(M_3(0:npol,0:npol,nel_solid))
  allocate(M_4(0:npol,0:npol,nel_solid))


  allocate(M_w1(0:npol,0:npol,nel_solid))

  allocate(M0_w1(0:npol,nel_solid))
  allocate(M0_w2(0:npol,nel_solid))
  allocate(M0_w3(0:npol,nel_solid))


  select case (src_type(1))

  case ('dipole')

      allocate(M13s(0:npol,0:npol,nel_solid))
      allocate(M33s(0:npol,0:npol,nel_solid))
      allocate(M43s(0:npol,0:npol,nel_solid))

      allocate(M_5(0:npol,0:npol,nel_solid))
      allocate(M_6(0:npol,0:npol,nel_solid))
      allocate(M_7(0:npol,0:npol,nel_solid))
      allocate(M_8(0:npol,0:npol,nel_solid))

      allocate(M_w2(0:npol,0:npol,nel_solid))
      allocate(M_w3(0:npol,0:npol,nel_solid))

      allocate(M0_w4(0:npol,nel_solid))
      allocate(M0_w5(0:npol,nel_solid))
      allocate(M0_w6(0:npol,nel_solid))
      allocate(M0_w7(0:npol,nel_solid))
      allocate(M0_w8(0:npol,nel_solid))
      allocate(M0_w9(0:npol,nel_solid))
      allocate(M0_w10(0:npol,nel_solid))

  case ('quadpole')

      allocate(M1phi(0:npol,0:npol,nel_solid))
      allocate(M2phi(0:npol,0:npol,nel_solid))
      allocate(M4phi(0:npol,0:npol,nel_solid))

      allocate(M_5(0:npol,0:npol,nel_solid))
      allocate(M_6(0:npol,0:npol,nel_solid))
      allocate(M_7(0:npol,0:npol,nel_solid))
      allocate(M_8(0:npol,0:npol,nel_solid))

      allocate(M_w2(0:npol,0:npol,nel_solid))
      allocate(M_w3(0:npol,0:npol,nel_solid))
      allocate(M_w4(0:npol,0:npol,nel_solid))
      allocate(M_w5(0:npol,0:npol,nel_solid))

      allocate(M0_w4(0:npol,nel_solid))
      allocate(M0_w5(0:npol,nel_solid))
      allocate(M0_w6(0:npol,nel_solid))

  end select

  ! Allocate Global Anelastic Stiffness Arrays
  if (anel_true) then
     allocate(Y(0:npol,0:npol,nel_solid))
     allocate(V_s_eta(0:npol,0:npol,nel_solid))
     allocate(V_s_xi(0:npol,0:npol,nel_solid))
     allocate(V_z_eta(0:npol,0:npol,nel_solid))
     allocate(V_z_xi(0:npol,0:npol,nel_solid))

     allocate(Y0(0:npol,nel_solid))
     allocate(V0_s_eta(0:npol,nel_solid))
     allocate(V0_s_xi(0:npol,nel_solid))
     allocate(V0_z_eta(0:npol,nel_solid))
     allocate(V0_z_xi(0:npol,nel_solid))

     Y = 0
     V_s_eta = 0
     V_s_xi = 0
     V_z_eta = 0
     V_z_xi = 0

     Y0 = 0
     V0_s_eta = 0
     V0_s_xi = 0
     V0_z_eta = 0
     V0_z_xi = 0
  endif

  ! NON DIAG FAC
  allocate(non_diag_fact(0:npol,1:nel_solid))
  non_diag_fact(:,:) = zero
  do ielem = 1, nel_solid
     if ( axis_solid(ielem) ) then
        do inode = 1, 8
           call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                local_crd_nodes(inode,2),ielsolid(ielem),inode)
        enddo
        do jpol = 0,npol
           non_diag_fact(jpol,ielem) = wt_axial_k(0)*wt(jpol)* &
              jacobian(xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))/&
              s_over_oneplusxi_axis(xi_k(0),eta(jpol), &
                                    local_crd_nodes,ielsolid(ielem))
        enddo
     endif
  enddo

  ! SOLID STIFFNESS TERMS

  do ielem=1, nel_solid

     alpha_wt_k(:,:) = zero
     beta_wt_k(:,:)  = zero
     gamma_wt_k(:,:) = zero
     delta_wt_k(:,:) = zero
     epsil_wt_k(:,:) = zero
     zeta_wt_k(:,:)  = zero

     Ms_z_eta_s_xi_wt_k(:,:) = zero
     Ms_z_eta_s_eta_wt_k(:,:) = zero
     Ms_z_xi_s_eta_wt_k(:,:) = zero
     Ms_z_xi_s_xi_wt_k(:,:) = zero

     M_s_xi_wt_k(:,:) = zero
     M_z_xi_wt_k(:,:) = zero
     M_z_eta_wt_k(:,:) = zero
     M_s_eta_wt_k(:,:) = zero

     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1), &
             local_crd_nodes(inode,2),ielsolid(ielem),inode)
     enddo

     ! ::::::::::::::::non-axial elements::::::::::::::::
     if (.not. axis_solid(ielem) ) then
        do ipol  = 0, npol
           do jpol = 0, npol

              call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                   eta(ipol),eta(jpol),local_crd_nodes,ielsolid(ielem))

              alpha_wt_k(ipol,jpol) = alpha(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              beta_wt_k(ipol,jpol) = beta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              gamma_wt_k(ipol,jpol) = gamma1(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              delta_wt_k(ipol,jpol) = delta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              epsil_wt_k(ipol,jpol) = epsilon1(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              zeta_wt_k(ipol,jpol) = zeta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_eta_s_xi_wt_k(ipol,jpol) = Ms_z_eta_s_xi(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_eta_s_eta_wt_k(ipol,jpol) = Ms_z_eta_s_eta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_xi_s_eta_wt_k(ipol,jpol) = Ms_z_xi_s_eta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              Ms_z_xi_s_xi_wt_k(ipol,jpol) = Ms_z_xi_s_xi(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                   *wt(ipol)*wt(jpol)
              M_s_xi_wt_k(ipol,jpol)  =  dsdxi * wt(ipol) * wt(jpol)
              M_z_xi_wt_k(ipol,jpol)  = -dzdxi * wt(ipol) * wt(jpol)
              M_z_eta_wt_k(ipol,jpol) =  dzdeta * wt(ipol) * wt(jpol)
              M_s_eta_wt_k(ipol,jpol) = -dsdeta * wt(ipol) * wt(jpol)

              if (anel_true) then
                 Y(ipol,jpol,ielem) = wt(ipol) * wt(jpol) &
                    * jacobian(eta(ipol), eta(jpol), local_crd_nodes, ielsolid(ielem))

                 V_s_eta(ipol,jpol,ielem) &
                    = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                        * M_s_eta_wt_k(ipol,jpol)
                 V_s_xi(ipol,jpol,ielem) &
                    = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                        * M_s_xi_wt_k(ipol,jpol)
                 V_z_eta(ipol,jpol,ielem) &
                    = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                        * M_z_eta_wt_k(ipol,jpol)
                 V_z_xi(ipol,jpol,ielem) &
                    = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                        * M_z_xi_wt_k(ipol,jpol)

              endif
           enddo
        enddo

     ! ::::::::::::::::axial elements::::::::::::::::
     else if ( axis_solid(ielem) ) then
        do jpol  = 0, npol
           do ipol = 0, npol
              alpha_wt_k(ipol,jpol) = alphak(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              beta_wt_k(ipol,jpol) = betak(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              gamma_wt_k(ipol,jpol) = gammak(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              delta_wt_k(ipol,jpol) = deltak(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              epsil_wt_k(ipol,jpol) = epsilonk(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              zeta_wt_k(ipol,jpol) = zetak(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem)) &
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)

              Ms_z_eta_s_xi_wt_k(ipol,jpol) = Ms_z_eta_s_xi_k(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              Ms_z_eta_s_eta_wt_k(ipol,jpol) = Ms_z_eta_s_eta_k(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)
              Ms_z_xi_s_eta_wt_k(ipol,jpol) = Ms_z_xi_s_eta_k(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)

              Ms_z_xi_s_xi_wt_k(ipol,jpol) = Ms_z_xi_s_xi_k(xi_k(ipol), &
                 eta(jpol),local_crd_nodes,ielsolid(ielem))&
                 *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                 local_crd_nodes,ielsolid(ielem))*wt_axial_k(ipol)*wt(jpol)


              if (ipol > 0) then
                 call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
                      xi_k(ipol),eta(jpol),local_crd_nodes,ielsolid(ielem))
                 M_s_xi_wt_k(ipol,jpol)  = dsdxi*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_z_xi_wt_k(ipol,jpol)  = -dzdxi*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_z_eta_wt_k(ipol,jpol) = dzdeta*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))
                 M_s_eta_wt_k(ipol,jpol) = -dsdeta*wt_axial_k(ipol) &
                      *wt(jpol)/(one+xi_k(ipol))

                 if (anel_true) then
                    Y(ipol,jpol,ielem) = wt_axial_k(ipol) * wt(jpol) / (1 + xi_k(ipol))&
                       * jacobian(xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(ielem))

                    V_s_eta(ipol,jpol,ielem) &
                       = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                           * M_s_eta_wt_k(ipol,jpol)
                    V_s_xi(ipol,jpol,ielem) &
                       = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                           * M_s_xi_wt_k(ipol,jpol)
                    V_z_eta(ipol,jpol,ielem) &
                       = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                           * M_z_eta_wt_k(ipol,jpol)
                    V_z_xi(ipol,jpol,ielem) &
                       = mapping(eta(ipol), eta(jpol), local_crd_nodes,1,ielsolid(ielem)) &
                           * M_z_xi_wt_k(ipol,jpol)
                 endif
              else
                 M_s_xi_wt_k(ipol,jpol) = zero
                 M_z_xi_wt_k(ipol,jpol) = zero
                 M_z_eta_wt_k(ipol,jpol) = zero
                 M_s_eta_wt_k(ipol,jpol) = zero

                 if (anel_true) then
                    dsdxi = 0
                    dsdeta = 0
                    dzdeta = 0
                    dzdxi = 0
                    Y(ipol,jpol,ielem) = 0
                    V_s_eta(ipol,jpol,ielem) = 0
                    V_s_xi(ipol,jpol,ielem) = 0
                    V_z_eta(ipol,jpol,ielem) = 0
                    V_z_xi(ipol,jpol,ielem) = 0

                    Y0(jpol,ielem) = wt_axial_k(ipol) * wt(jpol) &
                       * jacobian(xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(ielem))

                    call compute_partial_derivatives(dsdxi, dzdxi, dsdeta, dzdeta, &
                         xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(ielem))

                    !V0_s_eta(jpol,ielem) = wt_axial_k(ipol) * wt(jpol) * dsdxi * (-dsdeta)
                    V0_s_eta(jpol,ielem) = 0 ! dsdeta = 0 at the axis
                    V0_s_xi(jpol,ielem)  = wt_axial_k(ipol) * wt(jpol) * dsdxi * dsdxi
                    V0_z_eta(jpol,ielem) = wt_axial_k(ipol) * wt(jpol) * dsdxi * dzdeta
                    V0_z_xi(jpol,ielem)  = wt_axial_k(ipol) * wt(jpol) * dsdxi * (-dzdxi)
                 endif
              endif
           enddo
        enddo
     endif ! axis_solid or not

     do jpol=0,npol
       select case (src_type(1))
       case ('monopole')
          call compute_monopole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                     lambda,mu,xi_ani,phi_ani,eta_ani, &
                                     fa_ani_theta, fa_ani_phi, &
                                     massmat_kwts2, &
                                     non_diag_fact,alpha_wt_k,beta_wt_k, &
                                     gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                     zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                     M_z_eta_wt_k,M_s_eta_wt_k, &
                                     Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                     Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

       case('dipole')
          call compute_dipole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                     lambda,mu,xi_ani,phi_ani,eta_ani, &
                                     fa_ani_theta, fa_ani_phi, &
                                     massmat_kwts2, &
                                     non_diag_fact,alpha_wt_k,beta_wt_k, &
                                     gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                     zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                     M_z_eta_wt_k,M_s_eta_wt_k, &
                                     Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                     Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

       case('quadpole')
          call compute_quadrupole_stiff_terms(ielem,jpol, &
                                     lambda,mu,xi_ani,phi_ani,eta_ani, &
                                     fa_ani_theta, fa_ani_phi, &
                                     massmat_kwts2, &
                                     non_diag_fact,alpha_wt_k,beta_wt_k, &
                                     gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                     zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                     M_z_eta_wt_k,M_s_eta_wt_k, &
                                     Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                     Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)
      end select

     enddo  ! jpol
  enddo ! solid elements

  if (anel_true .and. att_coarse_grained) then
     allocate(Y_cg4(1:4,nel_solid))
     allocate(V_s_eta_cg4(1:4,nel_solid))
     allocate(V_s_xi_cg4(1:4,nel_solid))
     allocate(V_z_eta_cg4(1:4,nel_solid))
     allocate(V_z_xi_cg4(1:4,nel_solid))

     Y_cg4(1,:) = Y(1,1,:)
     Y_cg4(2,:) = Y(1,3,:)
     Y_cg4(3,:) = Y(3,1,:)
     Y_cg4(4,:) = Y(3,3,:)

     V_s_eta_cg4(1,:) = V_s_eta(1,1,:)
     V_s_eta_cg4(2,:) = V_s_eta(1,3,:)
     V_s_eta_cg4(3,:) = V_s_eta(3,1,:)
     V_s_eta_cg4(4,:) = V_s_eta(3,3,:)

     V_s_xi_cg4(1,:)  = V_s_xi(1,1,:)
     V_s_xi_cg4(2,:)  = V_s_xi(1,3,:)
     V_s_xi_cg4(3,:)  = V_s_xi(3,1,:)
     V_s_xi_cg4(4,:)  = V_s_xi(3,3,:)

     V_z_eta_cg4(1,:) = V_z_eta(1,1,:)
     V_z_eta_cg4(2,:) = V_z_eta(1,3,:)
     V_z_eta_cg4(3,:) = V_z_eta(3,1,:)
     V_z_eta_cg4(4,:) = V_z_eta(3,3,:)

     V_z_xi_cg4(1,:)  = V_z_xi(1,1,:)
     V_z_xi_cg4(2,:)  = V_z_xi(1,3,:)
     V_z_xi_cg4(3,:)  = V_z_xi(3,1,:)
     V_z_xi_cg4(4,:)  = V_z_xi(3,3,:)

     deallocate(Y)
     deallocate(V_s_eta)
     deallocate(V_s_xi)
     deallocate(V_z_eta)
     deallocate(V_z_xi)

     deallocate(Y0)
     deallocate(V0_s_eta)
     deallocate(V0_s_xi)
     deallocate(V0_z_eta)
     deallocate(V0_z_xi)

  endif

  if (verbose > 1) then
     write(69,*) ' '
     write(69,*) 'Min/max M11s [kg/s^2]:', minval(M11s),maxval(M11s)
     write(69,*) ' '
  endif

  deallocate(non_diag_fact)

end subroutine def_solid_stiffness_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_monopole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                       lambda,mu,xi_ani,phi_ani,eta_ani, &
                                       fa_ani_theta, fa_ani_phi, &
                                       massmat_kwts2, &
                                       non_diag_fact,alpha_wt_k,beta_wt_k, &
                                       gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                       zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                       M_z_eta_wt_k,M_s_eta_wt_k, &
                                       Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                       Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

  use data_mesh, only: npol, nelem, nel_solid, nel_fluid
  use data_matr, only: M0_w1, M0_w2, M0_w3, &
                       M11s, M21s, M41s, M12s, M22s, M32s, M42s, M11z, M21z, M41z, &
                       M_1, M_2, M_3, M_4, M_w1, M0_w1

  integer, intent(in) :: ielem, jpol

  real(kind=dp), intent(in) :: lambda(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: mu(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: xi_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: phi_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: eta_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: fa_ani_theta(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: fa_ani_phi(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

  real(kind=dp), intent(in) :: non_diag_fact(0:npol,nel_solid)
  real(kind=dp), intent(in) :: local_crd_nodes(8,2)

  real(kind=dp), intent(in) :: alpha_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: beta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: gamma_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: delta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: epsil_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: zeta_wt_k(0:npol,0:npol)

  real(kind=dp), intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

  real(kind=dp), intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

  integer          :: ipol
  real(kind=dp)    :: dsdxi,dzdeta,dzdxi,dsdeta
  real(kind=dp)    :: fa_ani_thetal, fa_ani_phil
  real(kind=dp)    :: C11, C22, C33, C12, C13, C23, C15, C25, C35, C44, C46, C55, C66, Ctmp
  real(kind=dp)    :: lambdal, mul, xil, phil, etal

  ! Clumsy to initialize inside a loop... but hey, the easiest way in this setup.
  if ( ielem == 1 .and. jpol == 0 ) then
     M0_w1(0:npol,1:nel_solid) = zero
     M0_w2(0:npol,1:nel_solid) = zero
     M0_w3(0:npol,1:nel_solid) = zero
  endif

  do ipol=0, npol
     fa_ani_thetal = fa_ani_theta(ipol,jpol,ielsolid(ielem))
     fa_ani_phil = fa_ani_phi(ipol,jpol,ielsolid(ielem))

     lambdal = lambda(ipol,jpol,ielsolid(ielem))
     mul = mu(ipol,jpol,ielsolid(ielem))
     xil = xi_ani(ipol,jpol,ielsolid(ielem))
     phil = phi_ani(ipol,jpol,ielsolid(ielem))
     etal = eta_ani(ipol,jpol,ielsolid(ielem))

     C11 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 1, 1)
     C12 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 2, 2)
     C13 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 3)
     C15 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 1)
     C22 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 2, 2)
     C23 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 3)
     C25 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 1)
     C33 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 3)
     C35 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 1)
     C44 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 2, 3)
     C46 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 1, 2)
     C55 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 1, 3, 1)
     C66 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 2, 1, 2)

     ! Test for the components that should be zero:
     if (do_mesh_tests) then
        if ( ielem == 1 .and. jpol == 0 .and. ipol == 0 ) then
           if (lpr) write(*,*) &
                ' Test for the components of c_ijkl that should be zero in anisotropic case'
        endif
        Ctmp = zero
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 3, 3, 1))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 1, 1, 2))

        if (Ctmp > smallval_sngl) then
           write(*,*) procstrg, ' ERROR: some stiffness term that should be zero '
           write(*,*) procstrg, '        is not: in compute_monopole_stiff_terms()'
           stop
        endif
     endif

     M11s(ipol,jpol,ielem) = C11 * delta_wt_k(ipol,jpol) &
                           + C15 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C15 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * alpha_wt_k(ipol,jpol)

     M21s(ipol,jpol,ielem) = C11 * zeta_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C55 * gamma_wt_k(ipol,jpol)

     M41s(ipol,jpol,ielem) = C11 * epsil_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C55 * beta_wt_k(ipol,jpol)

     M12s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C55 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M22s(ipol,jpol,ielem) = C15 * zeta_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C35 * gamma_wt_k(ipol,jpol)

     M32s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M42s(ipol,jpol,ielem) = C15 * epsil_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C35 * beta_wt_k(ipol,jpol)

     M11z(ipol,jpol,ielem) = C55 * delta_wt_k(ipol,jpol) &
                           + C35 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C33 * alpha_wt_k(ipol,jpol)

     M21z(ipol,jpol,ielem) = C55 * zeta_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C33 * gamma_wt_k(ipol,jpol)

     M41z(ipol,jpol,ielem) = C55 * epsil_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C33 * beta_wt_k(ipol,jpol)


     M_1(ipol,jpol,ielem) = C12 * M_z_eta_wt_k(ipol,jpol) + C25  * M_s_eta_wt_k(ipol,jpol)
     M_2(ipol,jpol,ielem) = C12 * M_z_xi_wt_k(ipol,jpol)  + C25  * M_s_xi_wt_k(ipol,jpol)
     M_3(ipol,jpol,ielem) = C23 * M_s_eta_wt_k(ipol,jpol) + C25  * M_z_eta_wt_k(ipol,jpol)
     M_4(ipol,jpol,ielem) = C23 * M_s_xi_wt_k(ipol,jpol)  + C25  * M_z_xi_wt_k(ipol,jpol)

     M_w1(ipol,jpol,ielem) = C22 * massmat_kwts2(ipol,jpol,ielsolid(ielem))

     if (axis_solid(ielem)) M_w1(0,jpol,ielem) = zero

     if (axis_solid(ielem) .and. ipol == 0) then
        call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
             xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))

        M0_w1(jpol,ielem) = (2 * C12 + C22) * non_diag_fact(jpol,ielem)
        M0_w2(jpol,ielem) = C25 * non_diag_fact(jpol,ielem)
        M0_w3(jpol,ielem) = C23 * dsdxi * wt_axial_k(0) * wt(jpol) &
                          + C25 * dzdxi * wt_axial_k(0) * wt(jpol)
     endif
  enddo ! ipol

end subroutine compute_monopole_stiff_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_dipole_stiff_terms(ielem,jpol,local_crd_nodes, &
                                      lambda,mu,xi_ani,phi_ani,eta_ani, &
                                      fa_ani_theta, fa_ani_phi, &
                                      massmat_kwts2, &
                                      non_diag_fact,alpha_wt_k,beta_wt_k, &
                                      gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                      zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                      M_z_eta_wt_k,M_s_eta_wt_k, &
                                      Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                      Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

  use data_matr

  integer, intent(in) :: ielem,jpol

  real(kind=dp), intent(in) :: lambda(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: mu(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: xi_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: phi_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: eta_ani(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: fa_ani_theta(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: fa_ani_phi(0:npol,0:npol,nelem)
  real(kind=dp), intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

  real(kind=dp), intent(in) :: non_diag_fact(0:npol,nel_solid)
  real(kind=dp), intent(in) :: local_crd_nodes(8,2)

  real(kind=dp), intent(in) :: alpha_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: beta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: gamma_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: delta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: epsil_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: zeta_wt_k(0:npol,0:npol)

  real(kind=dp), intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

  real(kind=dp), intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
  real(kind=dp), intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

  integer          :: ipol
  real(kind=dp)    :: dsdxi, dzdeta, dzdxi, dsdeta
  real(kind=dp)    :: fa_ani_thetal, fa_ani_phil
  real(kind=dp)    :: C11, C22, C33, C12, C13, C23, C15, C25, C35, C44, C46, C55, C66, Ctmp
  real(kind=dp)    :: lambdal, mul, xil, phil, etal

  if ( ielem == 1 .and. jpol == 0 ) then
     M0_w1(0:npol,1:nel_solid) = zero
     M0_w2(0:npol,1:nel_solid) = zero
     M0_w3(0:npol,1:nel_solid) = zero
     M0_w4(0:npol,1:nel_solid) = zero
     M0_w5(0:npol,1:nel_solid) = zero
     M0_w6(0:npol,1:nel_solid) = zero
     M0_w7(0:npol,1:nel_solid) = zero
     M0_w8(0:npol,1:nel_solid) = zero
     M0_w9(0:npol,1:nel_solid) = zero
     M0_w10(0:npol,1:nel_solid) = zero
  endif

  do ipol=0,npol
     fa_ani_thetal = fa_ani_theta(ipol,jpol,ielsolid(ielem))
     fa_ani_phil = fa_ani_phi(ipol,jpol,ielsolid(ielem))

     lambdal = lambda(ipol,jpol,ielsolid(ielem))
     mul = mu(ipol,jpol,ielsolid(ielem))
     xil = xi_ani(ipol,jpol,ielsolid(ielem))
     phil = phi_ani(ipol,jpol,ielsolid(ielem))
     etal = eta_ani(ipol,jpol,ielsolid(ielem))

     C11 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 1, 1)
     C12 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 2, 2)
     C13 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 3)
     C15 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 1)
     C22 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 2, 2)
     C23 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 3)
     C25 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 1)
     C33 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 3)
     C35 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 1)
     C44 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 2, 3)
     C46 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 1, 2)
     C55 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 1, 3, 1)
     C66 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 2, 1, 2)

     ! Test for the components that should be zero:
     if (do_mesh_tests) then
        if ( ielem == 1 .and. jpol == 0 .and. ipol == 0 ) then
           if (lpr) write(*,*) ' Test for the components of c_ijkl that should be zero in anisotropic case'
        endif
        Ctmp = zero
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 3, 3, 1))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 1, 1, 2))

        if (Ctmp > smallval_sngl) then
           write(*,*) procstrg, ' ERROR: some stiffness term that should be zero '
           write(*,*) procstrg, '        is not: in compute_dipole_stiff_terms()'
           stop
        endif
     endif

     ! + and - components
     M11s(ipol,jpol,ielem) = (C11 + C66) * delta_wt_k(ipol,jpol) &
                           + (C15 + C46) * (Ms_z_eta_s_xi_wt_k(ipol,jpol) &
                                            +  Ms_z_xi_s_eta_wt_k(ipol,jpol)) &
                           + (C55 + C44) * alpha_wt_k(ipol,jpol)

     M21s(ipol,jpol,ielem) = (C11 + C66) * zeta_wt_k(ipol,jpol) &
                           + (C15 + C46) * two * Ms_z_eta_s_eta_wt_k(ipol,jpol) &
                           + (C55 + C44) * gamma_wt_k(ipol,jpol)

     M41s(ipol,jpol,ielem) = (C11 + C66) * epsil_wt_k(ipol,jpol) &
                           + (C15 + C46) * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + (C55 + C44) * beta_wt_k(ipol,jpol)


     M12s(ipol,jpol,ielem) = (C11 - C66) * delta_wt_k(ipol,jpol) &
                           + (C15 - C46) * (Ms_z_eta_s_xi_wt_k(ipol,jpol) &
                                            +  Ms_z_xi_s_eta_wt_k(ipol,jpol)) &
                           + (C55 - C44) * alpha_wt_k(ipol,jpol)

     M22s(ipol,jpol,ielem) = (C11 - C66) * zeta_wt_k(ipol,jpol) &
                           + (C15 - C46) * two * Ms_z_eta_s_eta_wt_k(ipol,jpol) &
                           + (C55 - C44) * gamma_wt_k(ipol,jpol)

     M42s(ipol,jpol,ielem) = (C11 - C66) * epsil_wt_k(ipol,jpol) &
                           + (C15 - C46) * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + (C55 - C44) * beta_wt_k(ipol,jpol)


     M13s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C55 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M32s(ipol,jpol,ielem) = C15 * zeta_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C35 * gamma_wt_k(ipol,jpol)

     M33s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M43s(ipol,jpol,ielem) = C15 * epsil_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C35 * beta_wt_k(ipol,jpol)


     ! z component
     M11z(ipol,jpol,ielem) = C55 * delta_wt_k(ipol,jpol) &
                           + C35 * (Ms_z_eta_s_xi_wt_k(ipol,jpol) &
                                    +  Ms_z_xi_s_eta_wt_k(ipol,jpol)) &
                           + C33 * alpha_wt_k(ipol,jpol)

     M21z(ipol,jpol,ielem) = C55 * zeta_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol) &
                           + C33 * gamma_wt_k(ipol,jpol)

     M41z(ipol,jpol,ielem) = C55 * epsil_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C33 * beta_wt_k(ipol,jpol)

     ! D^x_y terms
     M_1(ipol,jpol,ielem) = (C12 + C66) * two * M_z_eta_wt_k(ipol,jpol) &
                          + (C25 + C46) * two * M_s_eta_wt_k(ipol,jpol)
     M_2(ipol,jpol,ielem) = (C12 + C66) * two * M_z_xi_wt_k(ipol,jpol)  &
                          + (C25 + C46) * two * M_s_xi_wt_k(ipol,jpol)

     M_3(ipol,jpol,ielem) = C46 * M_z_eta_wt_k(ipol,jpol) &
                          + C44 * M_s_eta_wt_k(ipol,jpol)
     M_4(ipol,jpol,ielem) = C46 * M_z_xi_wt_k(ipol,jpol)  &
                          + C44 * M_s_xi_wt_k(ipol,jpol)

     M_5(ipol,jpol,ielem) = (C12 - C66) * two * M_z_eta_wt_k(ipol,jpol) &
                          + (C25 - C46) * two * M_s_eta_wt_k(ipol,jpol)
     M_6(ipol,jpol,ielem) = (C12 - C66) * two * M_z_xi_wt_k(ipol,jpol)  &
                          + (C25 - C46) * two * M_s_xi_wt_k(ipol,jpol)

     M_7(ipol,jpol,ielem) = C25 * two * M_z_eta_wt_k(ipol,jpol) &
                          + C23 * two * M_s_eta_wt_k(ipol,jpol)
     M_8(ipol,jpol,ielem) = C25 * two * M_z_xi_wt_k(ipol,jpol)  &
                          + C23 * two * M_s_xi_wt_k(ipol,jpol)

     ! 2nd order terms
     M_w1(ipol,jpol,ielem) = four * (C22 + C66) * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w2(ipol,jpol,ielem) = two * C46 * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w3(ipol,jpol,ielem) = C44 * massmat_kwts2(ipol,jpol,ielsolid(ielem))

     if (axis_solid(ielem) .and. ipol == 0) then

        call compute_partial_derivatives(dsdxi,dzdxi,dsdeta,dzdeta, &
             xi_k(0),eta(jpol),local_crd_nodes,ielsolid(ielem))

        M0_w1(jpol,ielem) = (C12 + C66) * two * non_diag_fact(jpol,ielem)
        M0_w2(jpol,ielem) = (C12 + C66) * two * dzdxi * wt_axial_k(0) * wt(jpol)

        M0_w3(jpol,ielem) = C46 * non_diag_fact(jpol,ielem)
        M0_w4(jpol,ielem) = C46 * dzdxi * wt_axial_k(0) * wt(jpol)

        M0_w6(jpol,ielem) = (C25 + C46) * two * dsdxi * wt_axial_k(0) * wt(jpol)

        M0_w7(jpol,ielem) = C44 * non_diag_fact(jpol,ielem)
        M0_w8(jpol,ielem) = C44 * dsdxi * wt_axial_k(0) * wt(jpol)

        M0_w9(jpol,ielem) = (C12 + C22) * four * non_diag_fact(jpol,ielem)

        M0_w10(jpol,ielem) = (two * C25 + C46) * non_diag_fact(jpol,ielem)

     endif
  enddo

  if ( axis_solid(ielem) ) then

     M11s(0,jpol,ielem) = zero
     M12s(0,jpol,ielem) = zero
     M42s(0,jpol,ielem) = zero
     M32s(0,jpol,ielem) = zero
     M43s(0,jpol,ielem) = zero
     M11z(0,jpol,ielem) = zero

     M_1(0,jpol,ielem) = zero
     M_2(0,jpol,ielem) = zero
     M_3(0,jpol,ielem) = zero
     M_4(0,jpol,ielem) = zero
     M_5(0,jpol,ielem) = zero
     M_6(0,jpol,ielem) = zero
     M_7(0,jpol,ielem) = zero
     M_8(0,jpol,ielem) = zero

     M_w1(0,jpol,ielem) = zero
     M_w3(0,jpol,ielem) = zero

  endif

end subroutine compute_dipole_stiff_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_quadrupole_stiff_terms(ielem,jpol, &
                                      lambda,mu,xi_ani,phi_ani,eta_ani, &
                                      fa_ani_theta, fa_ani_phi, &
                                      massmat_kwts2, &
                                      non_diag_fact,alpha_wt_k,beta_wt_k, &
                                      gamma_wt_k,delta_wt_k,epsil_wt_k, &
                                      zeta_wt_k,M_s_xi_wt_k,M_z_xi_wt_k, &
                                      M_z_eta_wt_k,M_s_eta_wt_k, &
                                      Ms_z_eta_s_xi_wt_k,Ms_z_eta_s_eta_wt_k, &
                                      Ms_z_xi_s_eta_wt_k,Ms_z_xi_s_xi_wt_k)

  use data_matr

  integer, intent(in)          :: ielem, jpol

  real(kind=dp)   , intent(in) :: lambda(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: mu(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: xi_ani(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: phi_ani(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: eta_ani(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: fa_ani_theta(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: fa_ani_phi(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: massmat_kwts2(0:npol,0:npol,nelem)

  real(kind=dp)   , intent(in) :: non_diag_fact(0:npol,nel_solid)

  real(kind=dp)   , intent(in) :: alpha_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: beta_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: gamma_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: delta_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: epsil_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: zeta_wt_k(0:npol,0:npol)

  real(kind=dp)   , intent(in) :: Ms_z_eta_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: Ms_z_eta_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: Ms_z_xi_s_eta_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: Ms_z_xi_s_xi_wt_k(0:npol,0:npol)

  real(kind=dp)   , intent(in) :: M_s_xi_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: M_z_xi_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: M_z_eta_wt_k(0:npol,0:npol)
  real(kind=dp)   , intent(in) :: M_s_eta_wt_k(0:npol,0:npol)

  integer          :: ipol
  real(kind=dp)    :: fa_ani_thetal, fa_ani_phil
  real(kind=dp)    :: C11, C22, C33, C12, C13, C23, C15, C25, C35, C44, C46, C55, C66, Ctmp
  real(kind=dp)    :: lambdal, mul, xil, phil, etal

  if ( ielem == 1 .and. jpol == 0 ) then
     M0_w1(0:npol,1:nel_solid) = zero
     M0_w2(0:npol,1:nel_solid) = zero
     M0_w3(0:npol,1:nel_solid) = zero
     M0_w4(0:npol,1:nel_solid) = zero
     M0_w5(0:npol,1:nel_solid) = zero
     M0_w6(0:npol,1:nel_solid) = zero
  endif

  do ipol = 0, npol
     fa_ani_thetal = fa_ani_theta(ipol,jpol,ielsolid(ielem))
     fa_ani_phil = fa_ani_phi(ipol,jpol,ielsolid(ielem))

     lambdal = lambda(ipol,jpol,ielsolid(ielem))
     mul = mu(ipol,jpol,ielsolid(ielem))
     xil = xi_ani(ipol,jpol,ielsolid(ielem))
     phil = phi_ani(ipol,jpol,ielsolid(ielem))
     etal = eta_ani(ipol,jpol,ielsolid(ielem))

     C11 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 1, 1)
     C12 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 2, 2)
     C13 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 3)
     C15 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 1, 3, 1)
     C22 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 2, 2)
     C23 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 3)
     C25 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 2, 3, 1)
     C33 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 3)
     C35 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 3, 3, 1)
     C44 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 2, 3)
     C46 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 2, 3, 1, 2)
     C55 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 3, 1, 3, 1)
     C66 = c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, fa_ani_phil, 1, 2, 1, 2)

    ! Test for the components that should be zero:
     if (do_mesh_tests) then
        if ( ielem == 1 .and. jpol == 0 .and. ipol == 0 ) then
           if (lpr) write(*,*) &
               ' Test for the components of c_ijkl that should be zero in anisotropic case'
        endif
        Ctmp = zero
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 1, 1, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 2, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 2, 3))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 3, 1, 2))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 2, 3, 3, 1))
        Ctmp = Ctmp + dabs(c_ijkl_ani(lambdal, mul, xil, phil, etal, fa_ani_thetal, &
                                      fa_ani_phil, 3, 1, 1, 2))

        if (Ctmp > smallval_sngl) then
           write(*,*)procstrg,' ERROR: some stiffness term that should be zero '
           write(*,*)procstrg,'        is not: in compute_quadrupole_stiff_terms()'
           stop
        endif
     endif

     M11s(ipol,jpol,ielem) = C11 * delta_wt_k(ipol,jpol) &
                           + C15 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C15 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * alpha_wt_k(ipol,jpol)

     M21s(ipol,jpol,ielem) = C11 * zeta_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C55 * gamma_wt_k(ipol,jpol)

     M41s(ipol,jpol,ielem) = C11 * epsil_wt_k(ipol,jpol) &
                           + C15 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C55 * beta_wt_k(ipol,jpol)


     M12s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C55 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M22s(ipol,jpol,ielem) = C15 * zeta_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C35 * gamma_wt_k(ipol,jpol)

     M32s(ipol,jpol,ielem) = C15 * delta_wt_k(ipol,jpol) &
                           + C13 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C55 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * alpha_wt_k(ipol,jpol)

     M42s(ipol,jpol,ielem) = C15 * epsil_wt_k(ipol,jpol) &
                           + (C13 + C55) * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C35 * beta_wt_k(ipol,jpol)


     M11z(ipol,jpol,ielem) = C55 * delta_wt_k(ipol,jpol) &
                           + C35 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                           + C35 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                           + C33 * alpha_wt_k(ipol,jpol)

     M21z(ipol,jpol,ielem) = C55 * zeta_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                           + C33 * gamma_wt_k(ipol,jpol)

     M41z(ipol,jpol,ielem) = C55 * epsil_wt_k(ipol,jpol) &
                           + C35 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                           + C33 * beta_wt_k(ipol,jpol)


     M1phi(ipol,jpol,ielem) = C66 * delta_wt_k(ipol,jpol) &
                            + C46 * Ms_z_eta_s_xi_wt_k(ipol,jpol)&
                            + C46 * Ms_z_xi_s_eta_wt_k(ipol,jpol)&
                            + C44 * alpha_wt_k(ipol,jpol)

     M2phi(ipol,jpol,ielem) = C66 * zeta_wt_k(ipol,jpol) &
                            + C46 * two * Ms_z_eta_s_eta_wt_k(ipol,jpol)&
                            + C44 * gamma_wt_k(ipol,jpol)

     M4phi(ipol,jpol,ielem) = C66 * epsil_wt_k(ipol,jpol) &
                            + C46 * two * Ms_z_xi_s_xi_wt_k(ipol,jpol)&
                            + C44 * beta_wt_k(ipol,jpol)

     M_1(ipol,jpol,ielem) = C12 * M_z_eta_wt_k(ipol,jpol) + C25 * M_s_eta_wt_k(ipol,jpol)
     M_2(ipol,jpol,ielem) = C12 * M_z_xi_wt_k(ipol,jpol)  + C25 * M_s_xi_wt_k(ipol,jpol)
     M_3(ipol,jpol,ielem) = C23 * M_s_eta_wt_k(ipol,jpol) + C25 * M_z_eta_wt_k(ipol,jpol)
     M_4(ipol,jpol,ielem) = C23 * M_s_xi_wt_k(ipol,jpol)  + C25 * M_z_xi_wt_k(ipol,jpol)

     M_5(ipol,jpol,ielem) = C66 * M_z_eta_wt_k(ipol,jpol) + C46 * M_s_eta_wt_k(ipol,jpol)
     M_6(ipol,jpol,ielem) = C66 * M_z_xi_wt_k(ipol,jpol)  + C46 * M_s_xi_wt_k(ipol,jpol)
     M_7(ipol,jpol,ielem) = C44 * M_s_eta_wt_k(ipol,jpol) + C46 * M_z_eta_wt_k(ipol,jpol)
     M_8(ipol,jpol,ielem) = C44 * M_s_xi_wt_k(ipol,jpol)  + C46 * M_z_xi_wt_k(ipol,jpol)

     M_w1(ipol,jpol,ielem) = (C22 + 4.d0 * C66)   * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w2(ipol,jpol,ielem) = - 2.d0 * (C22 + C66) * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w3(ipol,jpol,ielem) = 2.d0 * C46           * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w4(ipol,jpol,ielem) = (4.d0 * C22 + C66)   * massmat_kwts2(ipol,jpol,ielsolid(ielem))
     M_w5(ipol,jpol,ielem) = 4.d0 * C44           * massmat_kwts2(ipol,jpol,ielsolid(ielem))

     if (axis_solid(ielem) .and. ipol == 0) then
        M0_w1(jpol,ielem) = (2.d0 * C12 + C22 + 4.d0 * C66) * non_diag_fact(jpol,ielem)
        M0_w2(jpol,ielem) = - 2.d0 * (C12 + C22)            * non_diag_fact(jpol,ielem)
        M0_w3(jpol,ielem) = (C25 + 4.d0 * C46)              * non_diag_fact(jpol,ielem)
        M0_w4(jpol,ielem) = (4.d0 * C22 - C66)              * non_diag_fact(jpol,ielem)
        M0_w5(jpol,ielem) = - 2.d0 * C25                    * non_diag_fact(jpol,ielem)
        M0_w6(jpol,ielem) = 4.d0 * C44                      * non_diag_fact(jpol,ielem)
     endif

   enddo ! ipol

   if (axis_solid(ielem)) then
      M_1(0,jpol,ielem) = zero
      M_2(0,jpol,ielem) = zero
      M_3(0,jpol,ielem) = zero
      M_4(0,jpol,ielem) = zero
      M_5(0,jpol,ielem) = zero
      M_6(0,jpol,ielem) = zero
      M_7(0,jpol,ielem) = zero
      M_8(0,jpol,ielem) = zero

      M_w1(0,jpol,ielem) = zero
      M_w2(0,jpol,ielem) = zero
      M_w3(0,jpol,ielem) = zero
      M_w4(0,jpol,ielem) = zero
      M_w5(0,jpol,ielem) = zero
   endif

end subroutine compute_quadrupole_stiff_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
real(kind=dp) function c_ijkl_ani(lambda, mu, xi_ani, phi_ani, eta_ani, &
                                  theta_fa, phi_fa, i, j, k, l)
! < returns the stiffness tensor as defined in Nolet(2008), Eq. (16.2)
!! i, j, k and l should be in [1,3]
!
! MvD [Anisotropy Notes, p. 13.4]


  real(kind=dp), intent(in) :: lambda, mu, xi_ani, phi_ani, eta_ani
  real(kind=dp), intent(in) :: theta_fa, phi_fa
  integer, intent(in)       :: i, j, k, l
  real(kind=dp), dimension(1:3, 1:3) :: deltaf
  real(kind=dp), dimension(1:3) :: s

  deltaf = zero
  deltaf(1,1) = one
  deltaf(2,2) = one
  deltaf(3,3) = one

  s(1) = dcos(phi_fa) * dsin(theta_fa)
  s(2) = dsin(phi_fa) * dsin(theta_fa)
  s(3) = dcos(theta_fa)

  c_ijkl_ani = zero

  ! isotropic part:
  c_ijkl_ani = c_ijkl_ani + lambda * deltaf(i,j) * deltaf(k,l)

  c_ijkl_ani = c_ijkl_ani + mu * (deltaf(i,k) * deltaf(j,l) + deltaf(i,l) * deltaf(j,k))


  ! anisotropic part:
  ! in xi, phi, eta

  c_ijkl_ani = c_ijkl_ani &
      + ((eta_ani - one) * lambda + two * eta_ani * mu * (one - one / xi_ani)) &
          * (deltaf(i,j) * s(k) * s(l) + deltaf(k,l) * s(i) * s(j))

  c_ijkl_ani = c_ijkl_ani &
      + mu * (one / xi_ani - one) &
          * (deltaf(i,k) * s(j) * s(l) + deltaf(i,l) * s(j) * s(k) + &
             deltaf(j,k) * s(i) * s(l) + deltaf(j,l) * s(i) * s(k))

  c_ijkl_ani = c_ijkl_ani &
      + ((one - two * eta_ani + phi_ani) * (lambda + two * mu) &
              + (4. * eta_ani - 4.) * mu / xi_ani) &
          * (s(i) * s(j) * s(k) * s(l))

end function c_ijkl_ani
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_fluid_stiffness_terms(rho,massmat_kwts2)
! < Fluid precomputed matrices definitions for all sources.
!! Note that in this routine terms alpha etc. are scalars
!! (as opposed to the solid case of being elemental arrays).

  use data_matr
  real(kind=dp)   , intent(in)  :: rho(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in)  :: massmat_kwts2(0:npol,0:npol,nelem)

  real(kind=dp)   , allocatable :: non_diag_fact(:,:)
  real(kind=dp)                 :: local_crd_nodes(8,2)
  real(kind=dp)                 :: alpha_wt_k,beta_wt_k,gamma_wt_k
  real(kind=dp)                 :: delta_wt_k,epsil_wt_k,zeta_wt_k
  integer                       :: iel,ipol,jpol,inode

  allocate(M1chi_fl(0:npol,0:npol,nel_fluid))
  allocate(M2chi_fl(0:npol,0:npol,nel_fluid))
  allocate(M4chi_fl(0:npol,0:npol,nel_fluid))
  allocate(M_w_fl(0:npol,0:npol,nel_fluid))
  allocate(M0_w_fl(0:npol,nel_fluid))

  allocate(non_diag_fact(0:npol,1:nel_fluid))

  do iel=1,nel_fluid

     do inode = 1, 8
        call compute_coordinates_mesh(local_crd_nodes(inode,1), &
             local_crd_nodes(inode,2),ielfluid(iel),inode)
     enddo

     do jpol=0,npol
        do ipol=0, npol
           ! ::::::::::::::::non-axial elements::::::::::::::::
           if (.not. axis_fluid(iel) ) then

              alpha_wt_k = alpha(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              beta_wt_k = beta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              gamma_wt_k = gamma1(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              delta_wt_k = delta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              epsil_wt_k = epsilon1(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)
              zeta_wt_k = zeta(eta(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *wt(ipol)*wt(jpol)

           ! ::::::::::::::::axial elements::::::::::::::::
           else if (axis_fluid(iel) ) then
              alpha_wt_k = alphak(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              beta_wt_k = betak(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel))&
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              gamma_wt_k = gammak(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              delta_wt_k = deltak(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              epsil_wt_k = epsilonk(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel))&
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
              zeta_wt_k = zetak(xi_k(ipol), &
                   eta(jpol),local_crd_nodes,ielfluid(iel)) &
                   *s_over_oneplusxi_axis(xi_k(ipol),eta(jpol), &
                   local_crd_nodes,ielfluid(iel))*wt_axial_k(ipol)*wt(jpol)
           endif

           M1chi_fl(ipol,jpol,iel)=(delta_wt_k + alpha_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
           M2chi_fl(ipol,jpol,iel)=(zeta_wt_k + gamma_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
           M4chi_fl(ipol,jpol,iel)=(epsil_wt_k + beta_wt_k) / &
                rho(ipol,jpol,ielfluid(iel))
        enddo
     enddo
  enddo

  !  2nd order term for multipoles
  if (src_type(1) == 'dipole' .or. src_type(1) == 'quadpole') then

     non_diag_fact(:,:) = zero
     do iel=1,nel_fluid

        ! ::::::::::::::::all elements::::::::::::::::
        M_w_fl(:,:,iel) = massmat_kwts2(:,:,ielfluid(iel))/rho(:,:,ielfluid(iel))

        ! ::::::::::::::::axial elements::::::::::::::::
        if ( axis_fluid(iel) ) then
           do inode = 1, 8
              call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                   local_crd_nodes(inode,2),ielfluid(iel),inode)
           enddo
           do jpol = 0,npol
              non_diag_fact(jpol,iel) = wt_axial_k(0)*wt(jpol)* &
                   jacobian(xi_k(0),eta(jpol),local_crd_nodes,ielfluid(iel))/&
                   s_over_oneplusxi_axis(xi_k(0), &
                   eta(jpol),local_crd_nodes,ielfluid(iel))
           enddo

           ! axial masking of main term
           M_w_fl(0,0:npol,iel)=zero

           ! additional axial term
           M0_w_fl(:,iel)=non_diag_fact(:,iel)/rho(0,:,ielfluid(iel))
        endif ! axis
     enddo
     if (src_type(1) == 'quadpole') then
        M_w_fl=four*M_w_fl
        M0_w_fl=four*M0_w_fl
     endif
  endif ! multipole

  if (verbose > 1) then
    write(69,*) ' '
    write(69,*) 'Min/max M1chi_fl [m^4/kg]:',minval(M1chi_fl),maxval(M1chi_fl)
  endif

  deallocate(non_diag_fact)

end subroutine def_fluid_stiffness_terms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_solid_fluid_boundary_terms
! < Defines the 1-d vector-array bdry_matr which acts as the diagonal matrix
!! to accomodate the exchange of fields across solid-fluid boundaries
!! in both directions. Take note of the sign conventions in accordance with
!! those used in the time loop.

  use commun, only: psum
  use data_io
  use data_mesh, only: npol, nel_bdry
  use data_matr

  real(kind=dp)                :: local_crd_nodes(8,2)
  real(kind=dp)                :: s,z,r,theta,rf,thetaf
  real(kind=dp)                :: theta1,theta2,r1,r2,delta_th,bdry_sum
  integer                      :: iel,ielglob,ipol,inode,idom
  integer                      :: count_lower_disc,count_upper_disc

  allocate(bdry_matr(0:npol,nel_bdry,2))
  allocate(solflubdry_radius(nel_bdry))

  bdry_sum = zero

  count_lower_disc = 0
  count_upper_disc = 0


  ! check if proc has boundary elements
  if (have_bdry_elem) then

     do iel=1,nel_bdry

        ! Map from boundary to global indexing. Choice of the solid side is random...
        ielglob = ielsolid(bdry_solid_el(iel))

        ! closest to axis
        call compute_coordinates(s, z, r1, theta1, ielglob, 0, bdry_jpol_solid(iel))

        call compute_coordinates(s, z, rf, thetaf, ielfluid(bdry_fluid_el(iel)), &
                                 0, bdry_jpol_fluid(iel))

        ! test if the mapping of solid element & jpol numbers agrees for solid & fluid
        if (abs( (rf-r1) /r1 ) > 1.e-5 .or. abs((thetaf-theta1)) > 1.e-3) then
           write(*,*)
           write(*,*)procstrg,'Problem with boundary term mapping near axis!'
           write(*,*)procstrg,'radius,theta solid index:',r1/1.d3,theta1/pi*180.
           write(*,*)procstrg,'radius,theta fluid index:',rf/1.d3,thetaf/pi*180.
           write(*,*)procstrg,'Possible reason: doubling layer directly on the solid side of'
           write(*,*)procstrg,'                 solid/fluid boundary. Check your mesh!'
           write(*,*)procstrg,'                 see ticket 26'
           stop
        endif

        ! furthest away from axis
        call compute_coordinates(s, z, r2, theta2, ielglob, npol, bdry_jpol_solid(iel))
        call compute_coordinates(s, z, rf, thetaf, ielfluid(bdry_fluid_el(iel)), &
                                 npol, bdry_jpol_fluid(iel))

        ! test if the mapping of solid element & jpol numbers agrees for solid & fluid
        if (abs( (rf-r2) /r2 ) > 1.e-5 .or. abs((thetaf-theta2)) > 1.e-3) then
           write(*,*)
           write(*,*)procstrg,'Problem with boundary term mapping far axis!'
           write(*,*)procstrg,'radius,theta solid index:',r2/1.d3,theta2/pi*180.
           write(*,*)procstrg,'radius,theta fluid index:',rf/1.d3,thetaf/pi*180.
           stop
        endif

        if ( abs(r1-r2) > min_distance_dim) then
           write(*,*)
           write(*,*)procstrg,'Problem with S/F boundary element',ielglob
           write(*,*)procstrg,'radii at min./max theta are not equal!'
           write(*,*)procstrg,'r1,r2 [km],theta [deg]:', &
                               r1/1000.,r2/1000.,theta1*180./pi
           stop
        endif
        !1/2 comes from d_th=1/2 (th_2 - th_1)d_xi; abs() takes care of southern elems
        delta_th = half * abs(theta2 - theta1)

        ! ::::::::::::::::axial elements::::::::::::::::
        if ( axis(ielglob) ) then

           if (abs(sin(theta1)) * two * pi * r1 > min_distance_dim) then
              write(*,*)
              write(*,*)procstrg,'Problem with axial S/F boundary element',ielglob
              write(*,*)procstrg,'Min theta is not exactly on the axis'
              write(*,*)procstrg,'r [km],theta [deg]:',r1/1000.,theta1*180./pi
              stop
           endif

           do inode = 1, 8
              call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                                            local_crd_nodes(inode,2), ielglob, inode)
           enddo

           ! I>0 (off the axis)
           do ipol=1, npol
              call compute_coordinates(s, z, r, theta, ielglob, ipol, &
                                       bdry_jpol_solid(iel))
              if (abs(r - r1) > min_distance_dim) then
                 write(*,*)
                 write(*,*)procstrg,'Problem with axial S/F boundary element', &
                           ielglob
                 write(*,*)procstrg,'radius at ipol=',ipol,'different from ipol=0'
                 write(*,*)procstrg,'r,r1 [km],theta [deg]:',r/1000.,r1/1000., &
                                                              theta*180./pi
                 stop
              endif

              bdry_matr(ipol,iel,1) = delta_th * wt_axial_k(ipol) * dsin(theta) &
                                        / (one + xi_k(ipol)) * dsin(theta)
              bdry_matr(ipol,iel,2) = delta_th * wt_axial_k(ipol) * dsin(theta) &
                                        / (one + xi_k(ipol)) * dcos(theta)

              ! Test: Integrated over the whole boundary, the boundary term san sin/cos
              ! equals two: \int_0^pi sin(\theta) d\theta = two, i.e. 4 if 2 boundaries
              bdry_sum = bdry_sum + delta_th * wt_axial_k(ipol) * dsin(theta) / &
                                    (one + xi_k(ipol))
           enddo

           ! I=0 axis
           bdry_sum = bdry_sum + 1/r * delta_th * wt_axial_k(0) &
                            * s_over_oneplusxi_axis(xi_k(0), &
                                                    eta(bdry_jpol_solid(iel)), &
                                                    local_crd_nodes,ielglob)

           ! I=0 axis itself
           bdry_matr(0,iel,1) = zero ! note sin(0) = 0

           ! need factor 1/r to compensate for using s/(1+xi) rather than sin theta/(1+xi)
           ! note cos(0) = 1
           bdry_matr(0,iel,2)= 1/r * delta_th * wt_axial_k(0) &
                              * s_over_oneplusxi_axis(xi_k(0), &
                                                      eta(bdry_jpol_solid(iel)), &
                                                      local_crd_nodes,ielglob)

           if (verbose > 1) then
              write(69,*)
              write(69,11)'S/F axial r[km], theta[deg]   :',r1/1000.,theta1*180/pi
              write(69,10)'S/F axial element (bdry,glob.):',iel,ielglob
              write(69,9)'S/F deltath, sigma_0, s/(1+xi):',delta_th*180/pi, &
                                                        wt_axial_k(0),delta_th*r
           endif

10         format(a33,2(i8))
11         format(a33,2(1pe14.5))
9          format(a33,3(1pe14.5))

           ! testing some algebra...
           if (abs(s_over_oneplusxi_axis(xi_k(0),eta(bdry_jpol_solid(iel)), &
               local_crd_nodes,ielglob)-delta_th*r) > min_distance_dim) then
              write(*,*)
              write(*,*)procstrg, &
                        'Problem with some axialgebra/definitions, elem:',ielglob
              write(*,*)procstrg,'r [km],theta [deg]:',r1/1000.,theta1*180/pi
              write(*,*)procstrg,'s_0 / (1+xi_0)  =', &
                                            s_over_oneplusxi_axis(xi_k(0), &
                                            eta(bdry_jpol_solid(iel)), &
                                            local_crd_nodes,ielglob)
              write(*,*)procstrg,'1/2 theta2 r_sf =',delta_th*r
              write(*,*)procstrg,'...are not the same :('
              write(*,*)procstrg,'xi,eta,eltype', &
                         xi_k(0),eta(bdry_jpol_solid(iel)),eltype(ielglob)
              write(*,*)procstrg,'theta1,theta2',theta1*180/pi,theta2*180/pi
              write(*,*)procstrg,'s,z min:', &
                        minval(local_crd_nodes(:,1)),minval(local_crd_nodes(:,2))
              stop
           endif

           ! ::::::::::::::::non-axial elements::::::::::::::::
           else
              do ipol=0, npol

                 call compute_coordinates(s, z, r, theta, ielglob, ipol, &
                                          bdry_jpol_solid(iel))

                 if (abs(r - r1) > min_distance_dim) then
                    write(*,*)
                    write(*,*)procstrg, &
                              'Problem with non-axial S/F boundary element',ielglob
                    write(*,*)procstrg,'radius at ipol=',ipol,'different from ipol=0'
                    write(*,12)procstrg,'r,r1 [km],theta [deg]:',r/1000.,r1/1000., &
                                                                 theta*180./pi
                    stop
                 endif

                 bdry_matr(ipol,iel,1) = delta_th * wt(ipol) * dsin(theta) * dsin(theta)
                 bdry_matr(ipol,iel,2) = delta_th * wt(ipol) * dsin(theta) * dcos(theta)

              ! Test: Integrated over the whole boundary, the boundary term san sin/cos
              ! equals two: \int_0^pi sin(\theta) d\theta = two, i.e. 4 if 2 boundaries
              bdry_sum = bdry_sum + delta_th * wt(ipol) * dsin(theta)

           enddo
        endif ! ax/nonax

        ! Define the term such that B >( < ) 0 of solid above(below) fluid
        do idom=1,ndisc-1
           if (.not. solid_domain(idom) ) then
              ! run a check to make sure radius is either discontinuity
              if ( abs(r1-discont(idom)) > min_distance_dim .and. &
                   abs(r1-discont(idom+1)) > min_distance_dim ) then
                 write(*,*)
                 write(*,*)procstrg, &
                           'Problem: S/F boundary radius is not one of the'
                 write(*,*)procstrg, &
                           '         two discontinuities bounding the fluid!!'
                 write(*,*)procstrg,'   r,elem(loc,glob):',r1,iel,ielglob
                 write(*,*)procstrg,'Upper/lower discont:', &
                                                  discont(idom),discont(idom+1)
                 stop
              endif
              ! if current radius=bottom radius of fluid layer, set negative
              if (abs(r1-discont(idom+1)) < min_distance_dim) then
                 bdry_matr(:,iel,:) = -bdry_matr(:,iel,:)
                 count_lower_disc = count_lower_disc+1
              else ! element is in upper radius of fluid layer, keep positive
                 count_upper_disc=count_upper_disc+1
              endif
           endif
        enddo ! idom

        ! Factor r^2 stemming from the integration over spherical domain
        bdry_matr(0:npol,iel,1:2) = bdry_matr(0:npol,iel,1:2) * r * r
        solflubdry_radius(iel) = r

     enddo ! elements along solid-fluid boundary

  endif ! have_bdry_elem

  bdry_sum = psum(real(bdry_sum,kind=realkind))

  ! yet another check....see if # elements above fluid is multiple of # below
  ! or the same (this is the case for no coarsening layer in the fluid)
  if ((count_upper_disc /= count_lower_disc) &
        .and. mod(count_upper_disc,2*count_lower_disc) /= 0) then
     write(*,*)procstrg, &
               'Problem: Number of elements found to be at discont above fluid'
     write(*,*)procstrg, &
               '   is not an even multiple of or the same as elements found to be below fluid'
     write(*,*)procstrg, &
               '   check doubling layers'
     write(*,*)procstrg,'# elems above fluid:',count_upper_disc
     write(*,*)procstrg,'# elems below fluid:',count_lower_disc
     stop
  endif

12 format(a25,3(1pe14.6))

  if (verbose > 1) then
     write(69,*)
     write(69,*) '# bdry elems above fluid:',count_upper_disc
     write(69,*) '# bdry elems below fluid:',count_lower_disc
     write(69,*) 'Min/max of boundary term [m^2]:', &
                                               minval(bdry_matr),maxval(bdry_matr)
     write(69,*) 'Min location bdry term  :',minloc(bdry_matr)
     write(69,*) 'Max location bdry term  :',maxloc(bdry_matr)
     write(69,*) 'Integr. bdry term, diff :',bdry_sum,dbleabsreldiff(bdry_sum,four)
     write(69,*)
  endif

  if (.not. dblreldiff_small(bdry_sum,four) ) then
     if (lpr) then
        write(*,*)'WARNING: boundary term not all that precise!'
        write(*,*)' Term should equal four for 2 boundaries (i.e. int (sin) )'
        write(*,*)' Actual numerical value:',bdry_sum
     endif
     ! deactivating the stop, because this prevents simulation with other then 2
     ! solid/fluid boundaries
     !stop
  endif

  if (diagfiles) then

      if (verbose > 1) then
         write(69,*)
         write(69,*)'saving boundary matrix with solid radius/colatitude into ', &
                     'boundary_term_sol_'//appmynum//'.dat'
         write(69,*)'saving boundary matrix with fluid radius/colatitude into ', &
                     'boundary_term_flu_'//appmynum//'.dat'
         write(69,*)
      endif

      ! output boundary precomputable matrix with radius [k]m and colatitude [deg]
      open(unit=500+mynum,file=infopath(1:lfinfo)//'/boundary_term_sol'&
                               //appmynum//'.dat')
      open(unit=400+mynum,file=infopath(1:lfinfo)//'/boundary_term_flu'&
                               //appmynum//'.dat')

      do iel=1,nel_bdry
         ielglob=ielsolid(bdry_solid_el(iel))
         do ipol=0,npol
            write(500+mynum,15)rcoord(ipol,bdry_jpol_solid(iel),ielglob)/1.d3, &
                               thetacoord(ipol,bdry_jpol_solid(iel),ielglob)/pi*180., &
                               bdry_matr(ipol,iel,1),bdry_matr(ipol,iel,2)

            write(400+mynum,15)rcoord(ipol,bdry_jpol_fluid(iel), &
                                       ielfluid(bdry_fluid_el(iel)))/1.d3, &
                                thetacoord(ipol,bdry_jpol_fluid(iel), &
                                           ielfluid(bdry_fluid_el(iel)))/pi*180., &
                                bdry_matr(ipol,iel,1),bdry_matr(ipol,iel,2)
         enddo
      enddo
      close(500+mynum)
      close(400+mynum)
  endif

15 format(4(1pe12.4))

end subroutine def_solid_fluid_boundary_terms
!-----------------------------------------------------------------------------------------

end module def_precomp_terms
!=========================================================================================
