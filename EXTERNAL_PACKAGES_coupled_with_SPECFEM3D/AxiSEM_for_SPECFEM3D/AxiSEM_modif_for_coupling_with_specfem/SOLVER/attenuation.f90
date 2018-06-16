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
module attenuation
! < Variables and routines for viscoelastic wave propagation

  use global_parameters, only: realkind, third, sp, dp
  use data_io, only: verbose
  use data_time
  implicit none

  private
  public :: prepare_attenuation
  public :: n_sls_attenuation
  public :: dump_memory_vars
  public :: time_step_memvars
  public :: att_coarse_grained

  real(kind=dp), allocatable       :: y_j_attenuation(:)
  real(kind=dp), allocatable       :: w_j_attenuation(:), exp_w_j_deltat(:)
  real(kind=dp), allocatable       :: ts_fac_t(:), ts_fac_tm1(:)
  integer                          :: n_sls_attenuation
  logical                          :: do_corr_lowq, dump_memory_vars = .false.
  logical                          :: att_coarse_grained = .true.

  real(kind=realkind), allocatable :: src_dev_tm1_glob(:,:,:,:)
  real(kind=realkind), allocatable :: src_tr_tm1_glob(:,:,:)

  real(kind=realkind), allocatable :: src_dev_tm1_glob_cg4(:,:,:)
  real(kind=realkind), allocatable :: src_tr_tm1_glob_cg4(:,:)

contains

!-----------------------------------------------------------------------------------------
!> Wrapper routine to avoid if statements in the time loop
subroutine time_step_memvars(memvar, memvar_cg, disp, cg)
  use data_mesh, only: npol

  real(kind=realkind), intent(in)                  :: disp(*)
  real(kind=realkind), intent(inout), allocatable  :: memvar(:,:,:,:,:)
  real(kind=realkind), intent(inout), allocatable  :: memvar_cg(:,:,:,:)
  logical, intent(in)                              :: cg

  if (cg) then
     call time_step_memvars_cg4(memvar_cg, disp)
  else
     if (npol == 4) then
        call time_step_memvars_4(memvar, disp)
     else
        call time_step_memvars_generic(memvar, disp)
     endif
  endif


end subroutine time_step_memvars
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> analytical time integration of memory variables (linear interpolation for
!! the strain), coarse grained version
!! MvD, attenutation notes, p 13.2
subroutine time_step_memvars_cg4(memvar, disp)

  use data_time, only: deltat
  use data_matr, only: Q_mu, Q_kappa
  use data_matr, only: delta_mu_cg4, delta_kappa_cg4
  use data_mesh, only: axis_solid, nel_solid
  use data_source, only: src_type

  real(kind=realkind), intent(inout)    :: memvar(1:4,6,n_sls_attenuation,nel_solid)
  real(kind=realkind), intent(in)       :: disp(0:4,0:4,nel_solid,3)

  integer               :: iel, j
  real(kind=dp)         :: yp_j_mu(n_sls_attenuation)
  real(kind=dp)         :: yp_j_kappa(n_sls_attenuation)
  real(kind=dp)         :: a_j_mu(n_sls_attenuation)
  real(kind=dp)         :: a_j_kappa(n_sls_attenuation)

  real(kind=realkind)   :: grad_t_cg4(1:4,6)

  real(kind=realkind)   :: trace_grad_t(1:4)
  real(kind=realkind)   :: src_tr_t(1:4)
  real(kind=realkind)   :: src_tr_tm1(1:4)
  real(kind=realkind)   :: src_dev_t(1:4,6)
  real(kind=realkind)   :: src_dev_tm1(1:4,6)

  real(kind=realkind)   :: src_tr_buf(1:4)
  real(kind=realkind)   :: src_dev_buf(1:4,6)

  real(kind=realkind)   :: Q_mu_last, Q_kappa_last

  Q_mu_last = -1
  Q_kappa_last = -1

  do iel=1, nel_solid
     call compute_strain_att_el_cg4(disp(:,:,iel,:), grad_t_cg4, iel)

     ! compute local coefficients y_j for kappa and mu (only if different from
     ! previous element)
     if (Q_mu(iel) /= Q_mu_last) then
        Q_mu_last = Q_mu(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        else
           yp_j_mu = y_j_attenuation / Q_mu(iel)
        endif
        a_j_mu = yp_j_mu / sum(yp_j_mu)
     endif

     if (Q_kappa(iel) /= Q_kappa_last) then
        Q_kappa_last = Q_kappa(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
        else
           yp_j_kappa = y_j_attenuation / Q_kappa(iel)
        endif
        a_j_kappa = yp_j_kappa / sum(yp_j_kappa)
     endif

     trace_grad_t(:) = sum(grad_t_cg4(:,1:3), dim=2)

     ! analytical time stepping, monopole/isotropic hardcoded

     ! compute new source terms (excluding the weighting)
     src_tr_t(:) = 0
     src_dev_t(:,:) = 0

     src_tr_t(:) = delta_kappa_cg4(:,iel) * trace_grad_t(:)

     src_dev_t(:,1) = delta_mu_cg4(:,iel) * 2 * (grad_t_cg4(:,1) - trace_grad_t(:) * third)
     src_dev_t(:,2) = delta_mu_cg4(:,iel) * 2 * (grad_t_cg4(:,2) - trace_grad_t(:) * third)
     src_dev_t(:,3) = delta_mu_cg4(:,iel) * 2 * (grad_t_cg4(:,3) - trace_grad_t(:) * third)
     src_dev_t(:,5) = delta_mu_cg4(:,iel) * grad_t_cg4(:,5)

     if (src_type(1) /= 'monopole') then
        src_dev_t(:,4) = delta_mu_cg4(:,iel) * grad_t_cg4(:,4)
        src_dev_t(:,6) = delta_mu_cg4(:,iel) * grad_t_cg4(:,6)
     endif

     ! load old source terms
     src_tr_tm1(:) = src_tr_tm1_glob_cg4(:,iel)
     src_dev_tm1(:,:) = src_dev_tm1_glob_cg4(:,:,iel)

     do j=1, n_sls_attenuation
        ! do the timestep
        src_tr_buf(:) = ts_fac_t(j) * a_j_kappa(j) * src_tr_t(:) &
                        + ts_fac_tm1(j) * a_j_kappa(j) * src_tr_tm1(:)

        src_dev_buf(:,1:3) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,1:3) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,1:3)
        src_dev_buf(:,5) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,5) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,5)

        memvar(:,1,j,iel) = exp_w_j_deltat(j) * memvar(:,1,j,iel) &
                        + src_dev_buf(:,1) + src_tr_buf(:)
        memvar(:,2,j,iel) = exp_w_j_deltat(j) * memvar(:,2,j,iel) &
                        + src_dev_buf(:,2) + src_tr_buf(:)
        memvar(:,3,j,iel) = exp_w_j_deltat(j) * memvar(:,3,j,iel) &
                        + src_dev_buf(:,3) + src_tr_buf(:)
        memvar(:,5,j,iel) = exp_w_j_deltat(j) * memvar(:,5,j,iel) &
                        + src_dev_buf(:,5)

        if (src_type(1) /= 'monopole') then
           src_dev_buf(:,4) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,4) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,4)
           src_dev_buf(:,6) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,6) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,6)

           memvar(:,4,j,iel) = exp_w_j_deltat(j) * memvar(:,4,j,iel) &
                        + src_dev_buf(:,4)
           memvar(:,6,j,iel) = exp_w_j_deltat(j) * memvar(:,6,j,iel) &
                        + src_dev_buf(:,6)
        endif
     enddo

     ! save srcs for next iteration
     src_tr_tm1_glob_cg4(:,iel) = src_tr_t(:)
     src_dev_tm1_glob_cg4(:,:,iel) = src_dev_t(:,:)

  enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> analytical time integration of memory variables (linear interpolation for
!! the strain)
!! MvD, attenutation notes, p 13.2
!! explicitly optimized for npol = 4
subroutine time_step_memvars_4(memvar, disp)

  use data_time, only: deltat
  use data_matr, only: Q_mu, Q_kappa
  use data_matr, only: delta_mu, delta_kappa
  use data_mesh, only: axis_solid, nel_solid!, npol
  use data_source, only: src_type

  real(kind=realkind), intent(inout)    :: memvar(0:4,0:4,6,n_sls_attenuation,nel_solid)
  real(kind=realkind), intent(in)       :: disp(0:4,0:4,nel_solid,3)

  integer               :: iel, j
  real(kind=dp)         :: yp_j_mu(n_sls_attenuation)
  real(kind=dp)         :: yp_j_kappa(n_sls_attenuation)
  real(kind=dp)         :: a_j_mu(n_sls_attenuation)
  real(kind=dp)         :: a_j_kappa(n_sls_attenuation)

  real(kind=realkind)   :: grad_t(0:4,0:4,6)

  real(kind=realkind)   :: trace_grad_t(0:4,0:4)
  real(kind=realkind)   :: src_tr_t(0:4,0:4)
  real(kind=realkind)   :: src_tr_tm1(0:4,0:4)
  real(kind=realkind)   :: src_dev_t(0:4,0:4,6)
  real(kind=realkind)   :: src_dev_tm1(0:4,0:4,6)

  real(kind=realkind)   :: src_tr_buf(0:4,0:4)
  real(kind=realkind)   :: src_dev_buf(0:4,0:4,6)

  real(kind=realkind)   :: Q_mu_last, Q_kappa_last

  Q_mu_last = -1
  Q_kappa_last = -1

  do iel=1, nel_solid
     call compute_strain_att_el_4(disp(:,:,iel,:), grad_t, iel)

     ! compute local coefficients y_j for kappa and mu (only if different from
     ! previous element)
     if (Q_mu(iel) /= Q_mu_last) then
        Q_mu_last = Q_mu(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        else
           yp_j_mu = y_j_attenuation / Q_mu(iel)
        endif
        a_j_mu = yp_j_mu / sum(yp_j_mu)
     endif

     if (Q_kappa(iel) /= Q_kappa_last) then
        Q_kappa_last = Q_kappa(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
        else
           yp_j_kappa = y_j_attenuation / Q_kappa(iel)
        endif
        a_j_kappa = yp_j_kappa / sum(yp_j_kappa)
     endif

     trace_grad_t(:,:) = sum(grad_t(:,:,1:3), dim=3)

     ! analytical time stepping, monopole/isotropic hardcoded

     ! compute new source terms
     src_tr_t(:,:) = delta_kappa(:,:,iel) * trace_grad_t(:,:)
     src_dev_t(:,:,1) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,1) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,2) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,2) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,3) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,3) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,5) = delta_mu(:,:,iel) * grad_t(:,:,5)

     if (src_type(1) /= 'monopole') then
        src_dev_t(:,:,4) = delta_mu(:,:,iel) * grad_t(:,:,4)
        src_dev_t(:,:,6) = delta_mu(:,:,iel) * grad_t(:,:,6)
     endif

     ! load old source terms
     src_tr_tm1(:,:) = src_tr_tm1_glob(:,:,iel)
     src_dev_tm1(:,:,:) = src_dev_tm1_glob(:,:,:,iel)


     do j=1, n_sls_attenuation
        ! do the timestep
        src_tr_buf(:,:) = ts_fac_t(j) * a_j_kappa(j) * src_tr_t(:,:) &
                        + ts_fac_tm1(j) * a_j_kappa(j) * src_tr_tm1(:,:)

        src_dev_buf(:,:,1:3) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,1:3) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,1:3)
        src_dev_buf(:,:,5) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,5) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,5)

        memvar(:,:,1,j,iel) = exp_w_j_deltat(j) * memvar(:,:,1,j,iel) &
                        + src_dev_buf(:,:,1) + src_tr_buf(:,:)
        memvar(:,:,2,j,iel) = exp_w_j_deltat(j) * memvar(:,:,2,j,iel) &
                        + src_dev_buf(:,:,2) + src_tr_buf(:,:)
        memvar(:,:,3,j,iel) = exp_w_j_deltat(j) * memvar(:,:,3,j,iel) &
                        + src_dev_buf(:,:,3) + src_tr_buf(:,:)
        memvar(:,:,5,j,iel) = exp_w_j_deltat(j) * memvar(:,:,5,j,iel) &
                        + src_dev_buf(:,:,5)

        ! maybe a bit uggly with the if inside the loop, but keeping it
        ! for now
        if (src_type(1) /= 'monopole') then

           src_dev_buf(:,:,4) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,4) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,4)
           src_dev_buf(:,:,6) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,6) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,6)

           memvar(:,:,4,j,iel) = exp_w_j_deltat(j) * memvar(:,:,4,j,iel) &
                        + src_dev_buf(:,:,4)
           memvar(:,:,6,j,iel) = exp_w_j_deltat(j) * memvar(:,:,6,j,iel) &
                        + src_dev_buf(:,:,6)

        endif
     enddo

     ! save srcs for next iteration
     src_tr_tm1_glob(:,:,iel) = src_tr_t(:,:)
     src_dev_tm1_glob(:,:,:,iel) = src_dev_t(:,:,:)
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> analytical time integration of memory variables (linear interpolation for
!! the strain)
!! MvD, attenutation notes, p 13.2
subroutine time_step_memvars_generic(memvar, disp)

  use data_time, only: deltat
  use data_matr, only: Q_mu, Q_kappa
  use data_matr, only: delta_mu, delta_kappa
  use data_mesh, only: axis_solid, nel_solid, npol
  use data_source, only: src_type

  real(kind=realkind), intent(inout)    :: memvar(0:npol,0:npol,6,n_sls_attenuation,nel_solid)
  real(kind=realkind), intent(in)       :: disp(0:npol,0:npol,nel_solid,3)

  integer               :: iel, j
  real(kind=dp)         :: yp_j_mu(n_sls_attenuation)
  real(kind=dp)         :: yp_j_kappa(n_sls_attenuation)
  real(kind=dp)         :: a_j_mu(n_sls_attenuation)
  real(kind=dp)         :: a_j_kappa(n_sls_attenuation)

  real(kind=realkind)   :: grad_t(0:npol,0:npol,6)

  real(kind=realkind)   :: trace_grad_t(0:npol,0:npol)
  real(kind=realkind)   :: src_tr_t(0:npol,0:npol)
  real(kind=realkind)   :: src_tr_tm1(0:npol,0:npol)
  real(kind=realkind)   :: src_dev_t(0:npol,0:npol,6)
  real(kind=realkind)   :: src_dev_tm1(0:npol,0:npol,6)

  real(kind=realkind)   :: src_tr_buf(0:npol,0:npol)
  real(kind=realkind)   :: src_dev_buf(0:npol,0:npol,6)

  real(kind=realkind)   :: Q_mu_last, Q_kappa_last

  Q_mu_last = -1
  Q_kappa_last = -1

  do iel=1, nel_solid
     call compute_strain_att_el(disp(:,:,iel,:), grad_t, iel)

     ! compute local coefficients y_j for kappa and mu (only if different from
     ! previous element)
     if (Q_mu(iel) /= Q_mu_last) then
        Q_mu_last = Q_mu(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        else
           yp_j_mu = y_j_attenuation / Q_mu(iel)
        endif
        a_j_mu = yp_j_mu / sum(yp_j_mu)
     endif

     if (Q_kappa(iel) /= Q_kappa_last) then
        Q_kappa_last = Q_kappa(iel)
        if (do_corr_lowq) then
           call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
        else
           yp_j_kappa = y_j_attenuation / Q_kappa(iel)
        endif
        a_j_kappa = yp_j_kappa / sum(yp_j_kappa)
     endif

     trace_grad_t(:,:) = sum(grad_t(:,:,1:3), dim=3)

     ! analytical time stepping, monopole/isotropic hardcoded

     ! compute new source terms
     src_tr_t(:,:) = delta_kappa(:,:,iel) * trace_grad_t(:,:)
     src_dev_t(:,:,1) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,1) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,2) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,2) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,3) = delta_mu(:,:,iel) * 2 * (grad_t(:,:,3) - trace_grad_t(:,:) * third)
     src_dev_t(:,:,5) = delta_mu(:,:,iel) * grad_t(:,:,5)

     if (src_type(1) /= 'monopole') then
        src_dev_t(:,:,4) = delta_mu(:,:,iel) * grad_t(:,:,4)
        src_dev_t(:,:,6) = delta_mu(:,:,iel) * grad_t(:,:,6)
     endif

     ! load old source terms
     src_tr_tm1(:,:) = src_tr_tm1_glob(:,:,iel)
     src_dev_tm1(:,:,:) = src_dev_tm1_glob(:,:,:,iel)

     do j=1, n_sls_attenuation
        ! do the timestep
        src_tr_buf(:,:) = ts_fac_t(j) * a_j_kappa(j) * src_tr_t(:,:) &
                        + ts_fac_tm1(j) * a_j_kappa(j) * src_tr_tm1(:,:)

        src_dev_buf(:,:,1:3) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,1:3) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,1:3)
        src_dev_buf(:,:,5) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,5) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,5)

        memvar(:,:,1,j,iel) = exp_w_j_deltat(j) * memvar(:,:,1,j,iel) &
                        + src_dev_buf(:,:,1) + src_tr_buf(:,:)
        memvar(:,:,2,j,iel) = exp_w_j_deltat(j) * memvar(:,:,2,j,iel) &
                        + src_dev_buf(:,:,2) + src_tr_buf(:,:)
        memvar(:,:,3,j,iel) = exp_w_j_deltat(j) * memvar(:,:,3,j,iel) &
                        + src_dev_buf(:,:,3) + src_tr_buf(:,:)
        memvar(:,:,5,j,iel) = exp_w_j_deltat(j) * memvar(:,:,5,j,iel) &
                        + src_dev_buf(:,:,5)

        ! maybe a bit uggly with the if inside the loop, but keeping it
        ! for now
        if (src_type(1) /= 'monopole') then

           src_dev_buf(:,:,4) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,4) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,4)
           src_dev_buf(:,:,6) = &
                    ts_fac_t(j) * a_j_mu(j) * src_dev_t(:,:,6) &
                        + ts_fac_tm1(j) * a_j_mu(j) * src_dev_tm1(:,:,6)

           memvar(:,:,4,j,iel) = exp_w_j_deltat(j) * memvar(:,:,4,j,iel) &
                        + src_dev_buf(:,:,4)
           memvar(:,:,6,j,iel) = exp_w_j_deltat(j) * memvar(:,:,6,j,iel) &
                        + src_dev_buf(:,:,6)

        endif
     enddo

     ! save srcs for next iteration
     src_tr_tm1_glob(:,:,iel) = src_tr_t(:,:)
     src_dev_tm1_glob(:,:,:,iel) = src_dev_t(:,:,:)
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> compute strain in Voigt notation for a single element
!! (i.e. E1 = E11, E2 = E22, E3 = E33, E4 = 2E23, E5 = 2E31, E6 = 2E12)
!! coarse grained version, so only computed at GLL points (1,1), (1,3), (3,1) and (3,3)
subroutine compute_strain_att_el_cg4(u, grad_u, iel)

  use data_source, only: src_type
  use pointwise_derivatives, only: axisym_gradient_solid_el_cg4
  use pointwise_derivatives, only: f_over_s_solid_el_cg4


  real(kind=realkind), intent(in)   :: u(0:,0:,:)
  real(kind=realkind), intent(out)  :: grad_u(1:4,6)
  integer, intent(in)               :: iel

  real(kind=realkind)               :: grad_buff1(1:4,2)
  real(kind=realkind)               :: grad_buff2(1:4,2)

  grad_u(:,:) = 0

  ! s,z components, identical for all source types..........................
  if (src_type(1) == 'dipole') then
     call axisym_gradient_solid_el_cg4(u(:,:,1) + u(:,:,2), grad_buff1, iel)
  else
     ! 1: dsus, 2: dzus
     call axisym_gradient_solid_el_cg4(u(:,:,1), grad_buff1, iel)
  endif

  ! 1:dsuz, 2:dzuz
  call axisym_gradient_solid_el_cg4(u(:,:,3), grad_buff2, iel)

  grad_u(:,1) = grad_buff1(:,1)  ! dsus
  grad_u(:,3) = grad_buff2(:,2)  ! dzuz

  ! dsuz + dzus (factor of 2 from voigt notation)
  grad_u(:,5) = grad_buff1(:,2) + grad_buff2(:,1)

  ! Components involving phi....................................................

  if (src_type(1) == 'monopole') then
     ! us / s
     grad_u(:,2) = f_over_s_solid_el_cg4(u(:,:,1), iel)

  else if (src_type(1) == 'dipole') then
     ! 2 u- / s
     grad_u(:,2) = 2 * f_over_s_solid_el_cg4(u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el_cg4(u(:,:,1) - u(:,:,2), grad_buff1, iel)

     ! -uz/s - dzup
     grad_u(:,4) = - f_over_s_solid_el_cg4(u(:,:,3), iel) - grad_buff1(:,2)
     ! -2 u-/s - dsup
     grad_u(:,6) = - grad_u(:,2) - grad_buff1(:,1)

  else if (src_type(1) == 'quadpole') then
     ! (us - 2 up) / s
     grad_u(:,2) = f_over_s_solid_el_cg4(u(:,:,1) - 2 * u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el_cg4(u(:,:,2), grad_buff1, iel)

     ! -2 uz/s - dzup
     grad_u(:,4) = - 2 * f_over_s_solid_el_cg4(u(:,:,3), iel) - grad_buff1(:,2)
     ! (up - 2 us) /s - dsup
     grad_u(:,6) = f_over_s_solid_el_cg4(u(:,:,2) - 2 * u(:,:,1), iel) - grad_buff1(:,1)
  endif

end subroutine compute_strain_att_el_cg4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> compute strain in Voigt notation for a single element
!! (i.e. E1 = E11, E2 = E22, E3 = E33, E4 = 2E23, E5 = 2E31, E6 = 2E12)
!! explicitly for npol = 4
subroutine compute_strain_att_el_4(u, grad_u, iel)

  use data_source, only: src_type
  use pointwise_derivatives, only: axisym_gradient_solid_el_4
  use pointwise_derivatives, only: f_over_s_solid_el_4


  real(kind=realkind), intent(in)   :: u(0:,0:,:)
  real(kind=realkind), intent(out)  :: grad_u(0:4,0:4,6)
  integer, intent(in)               :: iel

  real(kind=realkind)               :: grad_buff1(0:4,0:4,2)
  real(kind=realkind)               :: grad_buff2(0:4,0:4,2)

  grad_u(:,:,:) = 0

  ! s,z components, identical for all source types..........................
  if (src_type(1) == 'dipole') then
     call axisym_gradient_solid_el_4(u(:,:,1) + u(:,:,2), grad_buff1, iel)
  else
     ! 1: dsus, 2: dzus
     call axisym_gradient_solid_el_4(u(:,:,1), grad_buff1, iel)
  endif

  ! 1:dsuz, 2:dzuz
  call axisym_gradient_solid_el_4(u(:,:,3), grad_buff2, iel)

  grad_u(:,:,1) = grad_buff1(:,:,1)  ! dsus
  grad_u(:,:,3) = grad_buff2(:,:,2)  ! dzuz

  ! dsuz + dzus (factor of 2 from voigt notation)
  grad_u(:,:,5) = grad_buff1(:,:,2) + grad_buff2(:,:,1)

  ! Components involving phi....................................................

  if (src_type(1) == 'monopole') then
     ! us / s
     grad_u(:,:,2) = f_over_s_solid_el_4(u(:,:,1), iel)

  else if (src_type(1) == 'dipole') then
     ! 2 u- / s
     grad_u(:,:,2) = 2 * f_over_s_solid_el_4(u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el_4(u(:,:,1) - u(:,:,2), grad_buff1, iel)
     ! -uz/s - dzup
     grad_u(:,:,4) = - f_over_s_solid_el_4(u(:,:,3), iel) - grad_buff1(:,:,2)
     ! -2 u-/s - dsup
     grad_u(:,:,6) = - grad_u(:,:,2) - grad_buff1(:,:,1)

  else if (src_type(1) == 'quadpole') then
     ! (us - 2 up) / s
     grad_u(:,:,2) = f_over_s_solid_el_4(u(:,:,1) - 2 * u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el_4(u(:,:,2), grad_buff1, iel)

     ! -2 uz/s - dzup
     grad_u(:,:,4) = - 2 * f_over_s_solid_el_4(u(:,:,3), iel) - grad_buff1(:,:,2)
     ! (up - 2 us) /s - dsup
     grad_u(:,:,6) = f_over_s_solid_el_4(u(:,:,2) - 2 * u(:,:,1), iel) - grad_buff1(:,:,1)

  endif

end subroutine compute_strain_att_el_4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> compute strain in Voigt notation for a single element
!! (i.e. E1 = E11, E2 = E22, E3 = E33, E4 = 2E23, E5 = 2E31, E6 = 2E12)
subroutine compute_strain_att_el(u, grad_u, iel)

  use data_source, only: src_type
  use pointwise_derivatives, only: axisym_gradient_solid_el
  use pointwise_derivatives, only: f_over_s_solid_el
  use data_mesh, only: npol


  real(kind=realkind), intent(in)   :: u(0:,0:,:)
  real(kind=realkind), intent(out)  :: grad_u(0:npol,0:npol,6)
  integer, intent(in)               :: iel

  real(kind=realkind)               :: grad_buff1(0:npol,0:npol,2)
  real(kind=realkind)               :: grad_buff2(0:npol,0:npol,2)

  grad_u(:,:,:) = 0

  ! s,z components, identical for all source types..........................
  if (src_type(1) == 'dipole') then
     call axisym_gradient_solid_el(u(:,:,1) + u(:,:,2), grad_buff1, iel)
  else
     ! 1: dsus, 2: dzus
     call axisym_gradient_solid_el(u(:,:,1), grad_buff1, iel)
  endif

  ! 1:dsuz, 2:dzuz
  call axisym_gradient_solid_el(u(:,:,3), grad_buff2, iel)

  grad_u(:,:,1) = grad_buff1(:,:,1)  ! dsus
  grad_u(:,:,3) = grad_buff2(:,:,2)  ! dzuz

  ! dsuz + dzus (factor of 2 from voigt notation)
  grad_u(:,:,5) = grad_buff1(:,:,2) + grad_buff2(:,:,1)

  ! Components involving phi....................................................

  if (src_type(1) == 'monopole') then
     ! us / s
     grad_u(:,:,2) = f_over_s_solid_el(u(:,:,1), iel)

  else if (src_type(1) == 'dipole') then
     ! 2 u- / s
     grad_u(:,:,2) = 2 * f_over_s_solid_el(u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el(u(:,:,1) - u(:,:,2), grad_buff1, iel)
     ! -uz/s - dzup
     grad_u(:,:,4) = - f_over_s_solid_el(u(:,:,3), iel) - grad_buff1(:,:,2)
     ! -2 u-/s - dsup
     grad_u(:,:,6) = - grad_u(:,:,2) - grad_buff1(:,:,1)

  else if (src_type(1) == 'quadpole') then
     ! (us - 2 up) / s
     grad_u(:,:,2) = f_over_s_solid_el(u(:,:,1) - 2 * u(:,:,2), iel)

     ! 1:dsup, 2:dzup
     call axisym_gradient_solid_el(u(:,:,2), grad_buff1, iel)

     ! -2 uz/s - dzup
     grad_u(:,:,4) = - 2 * f_over_s_solid_el(u(:,:,3), iel) - grad_buff1(:,:,2)
     ! (up - 2 us) /s - dsup
     grad_u(:,:,6) = f_over_s_solid_el(u(:,:,2) - 2 * u(:,:,1), iel) - grad_buff1(:,:,1)

  endif

end subroutine compute_strain_att_el
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> read 'inparam_attenuation' file and compute precomputable terms
subroutine prepare_attenuation(lambda, mu)

  use data_io, only: infopath, lfinfo, diagfiles
  use data_proc, only: lpr, mynum
  use data_time, only: deltat
  use global_parameters, only: pi
  use data_matr, only: Q_mu, Q_kappa
  use data_matr, only: delta_mu, delta_kappa
  use data_matr, only: delta_mu_cg4, delta_kappa_cg4
  use get_mesh, only: compute_coordinates_mesh
  use data_mesh, only: axis_solid, npol, nelem, ielsolid, nel_solid
  use data_spec, only: eta, xi_k, wt, wt_axial_k
  use utlity, only: scoord
  use analytic_mapping, only: compute_partial_derivatives, jacobian
  use data_source, only: nelsrc, ielsrc
  use data_pointwise
  use commun, only: broadcast_int, broadcast_log, &
                                  broadcast_char, broadcast_dble, barrier

  real(kind=dp), intent(inout)   :: lambda(0:,0:,1:)
  real(kind=dp), intent(inout)   :: mu(0:,0:,1:)

  real(kind=dp)                  :: mu_w1(0:npol,0:npol)

  real(kind=dp)                  :: delta_mu_0(0:npol,0:npol)
  real(kind=dp)                  :: kappa_w1(0:npol,0:npol)
  real(kind=dp)                  :: delta_kappa_0(0:npol,0:npol)

  real(kind=dp)                  :: kappa_fac, mu_fac

  real(kind=dp)                  :: f_min, f_max, w_1, w_0
  real(kind=dp)                  :: qpl_w_ref, qpl_alpha
  integer                        :: nfsamp, max_it, i, iel, j
  real(kind=dp)                  :: Tw, Ty, d
  logical                        :: fixfreq, freq_weight
  real(kind=dp), allocatable     :: w_samp(:), q_fit(:), chil(:)
  real(kind=dp), allocatable     :: yp_j_mu(:)
  real(kind=dp), allocatable     :: yp_j_kappa(:)

  real(kind=dp)                  :: local_crd_nodes(8,2)
  real(kind=dp)                  :: gamma_w_l(0:npol,0:npol)
  integer                        :: inode, ipol, jpol
  real(kind=dp)                  :: dsdxi, dzdxi, dsdeta, dzdeta
  real(kind=dp)                  :: weights_cg(0:npol,0:npol)

  integer                        :: iinparam_advanced=500, ioerr
  character(len=256)             :: line
  character(len=256)             :: keyword, keyvalue

  ! Default values
  n_sls_attenuation = 5
  f_min = 0.001
  f_max = 1.0
  w_0 = 1.0
  qpl_w_ref = 1
  qpl_alpha = 0
  do_corr_lowq = .true.
  nfsamp = 100
  max_it = 100000
  Tw = 0.1
  Ty = 0.1
  d = 0.99995
  fixfreq = .false.
  freq_weight = .true.
  dump_memory_vars = .false.
  att_coarse_grained = .true.

  call barrier ! only for nicer output

  if (mynum == 0) then
     keyword = ' '
     keyvalue = ' '

     if (verbose > 1) write(*, '(A)', advance='no') &
            '   Reading attenuation parameters from inparam_advanced...'
     open(unit=iinparam_advanced, file='inparam_advanced', status='old', action='read',  iostat=ioerr)
     if (ioerr /= 0) stop 'Check input file ''inparam_advanced''! Is it still there?'

     do
        read(iinparam_advanced,fmt='(a256)',iostat=ioerr) line
        if (ioerr < 0) exit
        if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
        read(line,*) keyword, keyvalue

        select case(keyword)
        case('NR_LIN_SOLIDS')
            read(keyvalue,*) n_sls_attenuation

        case('F_MIN')
            read(keyvalue,*) f_min

        case('F_MAX')
            read(keyvalue,*) f_max

        case('F_REFERENCE')
            read(keyvalue,*) w_0

        case('QPL_F_REFERENCE')
            read(keyvalue,*) qpl_w_ref

        case('QPL_ALPHA')
            read(keyvalue,*) qpl_alpha

        case('SMALL_Q_CORRECTION')
            read(keyvalue,*) do_corr_lowq

        case('NR_F_SAMPLE')
            read(keyvalue,*) nfsamp

        case('FREQ_WEIGHT')
            read(keyvalue,*) freq_weight

        case('MAXINT_SA')
            read(keyvalue,*) max_it

        case('TSTART_SR')
            read(keyvalue,*) Tw

        case('TSTART_AMP')
            read(keyvalue,*) Ty

        case('T_DECAY')
            read(keyvalue,*) d

        case('FIX_FREQ')
            read(keyvalue,*) fixfreq

        case('DUMP_VTK')
            read(keyvalue,*) dump_memory_vars

        case('COARSE_GRAINED')
            read(keyvalue,*) att_coarse_grained

        end select

     enddo
  endif ! mynum

  ! broadcast values to other processors
  call broadcast_int(n_sls_attenuation, 0)
  call broadcast_dble(f_min, 0)
  call broadcast_dble(f_max, 0)
  call broadcast_dble(w_0, 0)
  call broadcast_dble(qpl_w_ref, 0)
  call broadcast_dble(qpl_alpha, 0)
  call broadcast_log(do_corr_lowq, 0)
  call broadcast_int(nfsamp, 0)
  call broadcast_log(freq_weight, 0)
  call broadcast_int(max_it, 0)
  call broadcast_dble(Tw, 0)
  call broadcast_dble(Ty, 0)
  call broadcast_dble(d, 0)
  call broadcast_log(fixfreq, 0)
  call broadcast_log(dump_memory_vars, 0)
  call broadcast_log(att_coarse_grained, 0)

  if (lpr .and. verbose > 1) print *, 'done'

  w_0 = w_0 * (2 * pi)
  if (lpr .and. verbose > 1) print '(a,f6.3)', '       w_0 = ', w_0

  w_1 = dsqrt(f_min * f_max) * (2 * pi)
  if (lpr .and. verbose > 1) print '(a,f6.3)', '       w_1 = ', w_1

  if (f_min > f_max) then
     print *, "ERROR: minimum frequency larger then maximum frequency"
     print *, "       in inparam_attenuation:2-3"
     stop 2
  endif

  qpl_w_ref = qpl_w_ref * (2 * pi)

  allocate(w_samp(nfsamp))
  allocate(q_fit(nfsamp))
  allocate(chil(max_it))

  allocate(w_j_attenuation(1:n_sls_attenuation))
  allocate(exp_w_j_deltat(1:n_sls_attenuation))
  allocate(y_j_attenuation(1:n_sls_attenuation))

  allocate(yp_j_mu(1:n_sls_attenuation))
  allocate(yp_j_kappa(1:n_sls_attenuation))

  allocate(ts_fac_t(1:n_sls_attenuation))
  allocate(ts_fac_tm1(1:n_sls_attenuation))


  if (lpr .and. verbose > 1) print *, &
        '  inverting for standard linear solid parameters...'

  call invert_linear_solids(Q=1.d0, f_min=f_min, f_max=f_max, N=n_sls_attenuation, &
                            nfsamp=nfsamp, max_it=max_it, Tw=Tw, Ty=Ty, d=d, &
                            fixfreq=fixfreq, verbose=.false., exact=.false., &
                            freq_weight=freq_weight, w_ref=qpl_w_ref, alpha=qpl_alpha, &
                            w_j=w_j_attenuation, y_j=y_j_attenuation, w=w_samp, &
                            q_fit=q_fit, chil=chil)

  if (lpr .and. verbose > 1) print *, '  ...done'

  ! prefactors for the exact time stepping (att notes p 13.3)
  do j=1, n_sls_attenuation
     exp_w_j_deltat(j) = dexp(-w_j_attenuation(j) * deltat)
     ts_fac_tm1(j) = ((1 - exp_w_j_deltat(j)) / (w_j_attenuation(j) * deltat) &
                      - exp_w_j_deltat(j))
     ts_fac_t(j) = ((exp_w_j_deltat(j) - 1) / (w_j_attenuation(j) * deltat) + 1)
  enddo

  if (lpr) then
     if (verbose > 1) then
        print *, '  ...log-l2 misfit    : ', chil(max_it)
        print *, '  ...frequencies      : ', w_j_attenuation / (2. * pi)
        print *, '  ...coefficients y_j : ', y_j_attenuation
        print *, '  ...coarse grained   : ', att_coarse_grained

     endif
     if (diagfiles) then
        print *, '  ...writing fitted Q to file...'
        open(unit=165, file=infopath(1:lfinfo)//'/attenuation_q_fitted', status='replace')
        write(165,*) (w_samp(i), q_fit(i), char(10), i=1,nfsamp)
        close(unit=165)

        if (verbose > 1) print *, '  ...writing convergence of chi to file...'
        open(unit=166, file=infopath(1:lfinfo)//'/attenuation_convergence', status='replace')
        write(166,*) (chil(i), char(10), i=1,max_it)
        close(unit=166)
     endif
  endif

  if (lpr .and. verbose > 1) print *, '  ...calculating relaxed moduli...'

  if (att_coarse_grained) then

     if ( npol /= 4 ) then
        write(*,*) "ERROR: coarse grained memvar only implemented for npol = 4, but is ", npol
        stop 2
     endif

     allocate(delta_mu_cg4(1:4,nel_solid))
     allocate(delta_kappa_cg4(1:4,nel_solid))
     allocate(src_dev_tm1_glob_cg4(1:4,6,nel_solid))
     allocate(src_tr_tm1_glob_cg4(1:4,nel_solid))
     src_dev_tm1_glob_cg4 = 0
     src_tr_tm1_glob_cg4 = 0
  else
     allocate(delta_mu(0:npol,0:npol,nel_solid))
     allocate(delta_kappa(0:npol,0:npol,nel_solid))
     allocate(src_dev_tm1_glob(0:npol,0:npol,6,nel_solid))
     allocate(src_tr_tm1_glob(0:npol,0:npol,nel_solid))
     src_dev_tm1_glob = 0
     src_tr_tm1_glob = 0
  endif


  do iel=1, nel_solid

     if (att_coarse_grained) then
        !weighting for coarse grained memory vars (hard coded for polynomial order 4)

        do inode = 1, 8
           call compute_coordinates_mesh(local_crd_nodes(inode,1), &
                                         local_crd_nodes(inode,2), ielsolid(iel), inode)
        enddo

        if (.not. axis_solid(iel)) then ! non-axial elements
           do ipol=0, npol
              do jpol=0, npol
                 gamma_w_l(ipol, jpol) = wt(ipol) * wt(jpol) &
                       * jacobian(eta(ipol), eta(jpol), local_crd_nodes, ielsolid(iel)) &
                       * scoord(ipol,jpol,ielsolid(iel))
              enddo
           enddo
        else   ! axial elements
           do ipol=1, npol
              do jpol=0, npol
                 gamma_w_l(ipol, jpol) = wt_axial_k(ipol) * wt(jpol) / (1 + xi_k(ipol)) &
                       * jacobian(xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(iel)) &
                       * scoord(ipol,jpol,ielsolid(iel))
              enddo
           enddo
           ! extra axis terms
           ipol = 0
           do jpol=0, npol
              call compute_partial_derivatives(dsdxi, dzdxi, dsdeta, dzdeta, &
                      xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(iel))
              gamma_w_l(ipol, jpol) = wt_axial_k(ipol) * wt(jpol) &
                       * jacobian(xi_k(ipol), eta(jpol), local_crd_nodes, ielsolid(iel)) &
                       * dsdxi
           enddo
        endif

        weights_cg(:,:) = 0
        !! cg with 4 points out of 25
        weights_cg(1,1) = (   gamma_w_l(0,0) + gamma_w_l(0,1) &
                            + gamma_w_l(1,0) + gamma_w_l(1,1) &
                            + 0.5 * (  gamma_w_l(0,2) + gamma_w_l(1,2) &
                                     + gamma_w_l(2,0) + gamma_w_l(2,1)) &
                            + 0.25 * gamma_w_l(2,2) ) &
                          / gamma_w_l(1,1)

        weights_cg(1,3) = (   gamma_w_l(0,3) + gamma_w_l(0,4) &
                            + gamma_w_l(1,3) + gamma_w_l(1,4) &
                            + 0.5 * (  gamma_w_l(0,2) + gamma_w_l(1,2) &
                                     + gamma_w_l(2,3) + gamma_w_l(2,4)) &
                            + 0.25 * gamma_w_l(2,2) ) &
                          / gamma_w_l(1,3)

        weights_cg(3,1) = (   gamma_w_l(3,0) + gamma_w_l(3,1) &
                            + gamma_w_l(4,0) + gamma_w_l(4,1) &
                            + 0.5 * (  gamma_w_l(2,0) + gamma_w_l(2,1) &
                                     + gamma_w_l(3,2) + gamma_w_l(4,2)) &
                            + 0.25 * gamma_w_l(2,2) ) &
                          / gamma_w_l(3,1)

        weights_cg(3,3) = (   gamma_w_l(3,3) + gamma_w_l(3,4) &
                            + gamma_w_l(4,3) + gamma_w_l(4,4) &
                            + 0.5 * (  gamma_w_l(2,3) + gamma_w_l(2,4) &
                                     + gamma_w_l(3,2) + gamma_w_l(4,2)) &
                            + 0.25 * gamma_w_l(2,2) ) &
                          / gamma_w_l(3,3)

     endif ! att_coarse_grained

     if (do_corr_lowq) then
        call fast_correct(y_j_attenuation / Q_mu(iel), yp_j_mu)
        call fast_correct(y_j_attenuation / Q_kappa(iel), yp_j_kappa)
     else
       yp_j_mu = y_j_attenuation / Q_mu(iel)
       yp_j_kappa = y_j_attenuation / Q_kappa(iel)
     endif

     mu_fac = 0
     do i=1, n_sls_attenuation
        mu_fac = mu_fac + yp_j_mu(i) * w_j_attenuation(i)**2 &
                            / (w_1**2 + w_j_attenuation(i)**2)
     enddo
     mu_fac = mu_fac / sum(yp_j_mu)

     kappa_fac = 0
     do i=1, n_sls_attenuation
        kappa_fac = kappa_fac + yp_j_kappa(i) * w_j_attenuation(i)**2 &
                            / (w_1**2 + w_j_attenuation(i)**2)
     enddo
     kappa_fac = kappa_fac / sum(yp_j_kappa)

     ! compute moduli at central frequency w_1
     mu_w1(:,:) =  mu(:,:,ielsolid(iel)) * (1 + 2. / (pi * Q_mu(iel)) * log(w_1 / w_0))
     kappa_w1(:,:) =  (lambda(:,:,ielsolid(iel)) + 2.d0 / 3.d0 * mu(:,:,ielsolid(iel))) &
                         * (1 + 2. / (pi * Q_kappa(iel)) * log(w_1 / w_0))

     ! delta moduli
     delta_mu_0(:,:) = mu_w1(:,:) / (1.d0 / sum(yp_j_mu) + 1 - mu_fac)
     delta_kappa_0(:,:) = kappa_w1(:,:) / (1.d0 / sum(yp_j_kappa) + 1 - kappa_fac)

     if (att_coarse_grained) then
        ! compute unrelaxed moduli
        mu(:,:,ielsolid(iel)) = mu_w1(:,:) + weights_cg(:,:) * delta_mu_0(:,:) * mu_fac
        lambda(:,:,ielsolid(iel)) = kappa_w1(:,:) &
                                       + weights_cg(:,:) * delta_kappa_0(:,:) * kappa_fac &
                                       - 2.d0 / 3.d0 * mu(:,:,ielsolid(iel))

        ! weighted delta moduli
        delta_mu_cg4(1,iel) = weights_cg(1,1) * delta_mu_0(1,1)
        delta_mu_cg4(2,iel) = weights_cg(1,3) * delta_mu_0(1,3)
        delta_mu_cg4(3,iel) = weights_cg(3,1) * delta_mu_0(3,1)
        delta_mu_cg4(4,iel) = weights_cg(3,3) * delta_mu_0(3,3)

        delta_kappa_cg4(1,iel) = weights_cg(1,1) * delta_kappa_0(1,1)
        delta_kappa_cg4(2,iel) = weights_cg(1,3) * delta_kappa_0(1,3)
        delta_kappa_cg4(3,iel) = weights_cg(3,1) * delta_kappa_0(3,1)
        delta_kappa_cg4(4,iel) = weights_cg(3,3) * delta_kappa_0(3,3)
     else
        ! compute unrelaxed moduli
        mu(:,:,ielsolid(iel)) = mu_w1(:,:) + delta_mu_0(:,:) * mu_fac
        lambda(:,:,ielsolid(iel)) = kappa_w1(:,:) + delta_kappa_0(:,:) * kappa_fac &
                                       - 2.d0 / 3.d0 * mu(:,:,ielsolid(iel))

        ! delta moduli
        delta_mu(:,:,iel) = delta_mu_0(:,:)
        delta_kappa(:,:,iel) = delta_kappa_0(:,:)
     endif

  enddo


  if (att_coarse_grained) then
     allocate(DsDeta_over_J_sol_cg4(1:4,1:nel_solid))
     allocate(DzDeta_over_J_sol_cg4(1:4,1:nel_solid))
     allocate(DsDxi_over_J_sol_cg4(1:4,1:nel_solid))
     allocate(DzDxi_over_J_sol_cg4(1:4,1:nel_solid))

     DzDeta_over_J_sol_cg4(1,:) = DzDeta_over_J_sol(1,1,:)
     DzDeta_over_J_sol_cg4(2,:) = DzDeta_over_J_sol(1,3,:)
     DzDeta_over_J_sol_cg4(3,:) = DzDeta_over_J_sol(3,1,:)
     DzDeta_over_J_sol_cg4(4,:) = DzDeta_over_J_sol(3,3,:)

     DzDxi_over_J_sol_cg4(1,:) = DzDxi_over_J_sol(1,1,:)
     DzDxi_over_J_sol_cg4(2,:) = DzDxi_over_J_sol(1,3,:)
     DzDxi_over_J_sol_cg4(3,:) = DzDxi_over_J_sol(3,1,:)
     DzDxi_over_J_sol_cg4(4,:) = DzDxi_over_J_sol(3,3,:)

     DsDeta_over_J_sol_cg4(1,:) = DsDeta_over_J_sol(1,1,:)
     DsDeta_over_J_sol_cg4(2,:) = DsDeta_over_J_sol(1,3,:)
     DsDeta_over_J_sol_cg4(3,:) = DsDeta_over_J_sol(3,1,:)
     DsDeta_over_J_sol_cg4(4,:) = DsDeta_over_J_sol(3,3,:)

     DsDxi_over_J_sol_cg4(1,:) = DsDxi_over_J_sol(1,1,:)
     DsDxi_over_J_sol_cg4(2,:) = DsDxi_over_J_sol(1,3,:)
     DsDxi_over_J_sol_cg4(3,:) = DsDxi_over_J_sol(3,1,:)
     DsDxi_over_J_sol_cg4(4,:) = DsDxi_over_J_sol(3,3,:)
  endif

  if (lpr .and. verbose > 1) print *, '  ...DONE'

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> compute Q after (Emmerich & Korn, inverse of eq 21)
!! linearized version (exact = false) is eq 22 in E&K
pure subroutine q_linear_solid(y_j, w_j, w, exact, Qls)

  real(kind=dp)   , intent(in)    :: y_j(:), w_j(:), w(:)
  real(kind=dp)   , intent(out)   :: Qls(size(w))
  integer                         :: j

  logical, optional, intent(in)           :: exact
  !f2py logical, optional, intent(in)     :: exact = 0
  logical                                 :: exact_loc

  real(kind=dp)                   :: Qls_denom(size(w))

  if (present(exact)) then
      exact_loc = exact
  else
      exact_loc = .false.
  endif

  Qls = 1
  if (exact_loc) then
     do j=1, size(y_j)
        Qls = Qls + y_j(j) *  w**2 / (w**2 + w_j(j)**2)
     enddo
  endif

  Qls_denom = 0
  do j=1, size(y_j)
     Qls_denom = Qls_denom + y_j(j) * w * w_j(j) / (w**2 + w_j(j)**2)
  enddo

  Qls = Qls / Qls_denom
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> computes a first order correction to the linearized coefficients:
!! yp_j_corrected = y_j * delta_j
!! MvD Attenuation Notes, p. 17.3 bottom
pure subroutine fast_correct(y_j, yp_j)

  real(kind=dp), intent(in)    :: y_j(:)
  real(kind=dp), intent(out)   :: yp_j(size(y_j))

  real(kind=dp)                :: dy_j(size(y_j))
  integer                      :: k

  dy_j(1) = 1 + .5 * y_j(1)

  do k=2, size(y_j)
     dy_j(k) = dy_j(k-1) + (dy_j(k-1) - .5) * y_j(k-1) + .5 * y_j(k)
  enddo

  yp_j = y_j * dy_j

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> returns l2 misfit between Q_target and fitted Q using standard linear solids
pure subroutine l2_error(Q_target, Qls, weights, lse)

  real(kind=dp), intent(in)       :: Q_target(:), Qls(:), weights(:)

  real(kind=dp), intent(out)      :: lse
  integer                         :: nfsamp, i

  lse = 0
  nfsamp = size(Qls)

  ! log-l2 norm
  do i=1, nfsamp
     lse = lse + (log(Q_target(i) / Qls(i)))**2 * weights(i)
  enddo

  lse = lse / float(nfsamp)
  lse = dsqrt(lse)
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Inverts for constant Q, minimizing the L2 error for 1/Q using a simulated annealing
!! approach (varying peak frequencies and amplitudes).
subroutine invert_linear_solids(Q, f_min, f_max, N, nfsamp, max_it, Tw, Ty, d, &
                                fixfreq, verbose, exact, freq_weight, w_ref, alpha, w_j, &
                                y_j, w, q_fit, chil)

  ! Parameters:
  ! Q:              clear
  ! f_min, fmax:    frequency band (in Hz)
  ! N:              number of standard linear solids
  ! nfsamp:         number of sampling frequencies for computation of the misfit (log
  !                   spaced in freqeuncy band)
  ! max_it:         number of iterations
  ! Tw:             starting temperature for the frequencies
  ! Ty:             starting temperature for the amplitudes
  ! d:              temperature decay
  ! fixfreq:        use log spaced peak frequencies (fixed)
  ! verbose:        clear
  ! exact:          use exact relation for Q and the coefficients (Emmerich & Korn, eq
  !                   21). If false, linearized version is used (large Q approximation, eq
  !                   22).
  ! freq_weight:    use frequency weighting to ensure better fit at high frequencies
  ! w_ref:          reference angular frequency for power law Q
  ! alpha:          exponent for power law Q
  !
  ! returns:
  ! w_j:            relaxation frequencies, equals 1/tau_sigma in zener
  !                   formulation
  ! y_j:            coefficients of the linear solids, (Emmerich & Korn, eq 23 and 24)
  ! w:              sampling frequencies at which Q(w) is minimized
  ! q_fit:          resulting q(w) at these frequencies
  ! chil:           error as a function of iteration to check convergence,
  !                   Note that this version uses log-l2 norm!
  !
  use data_proc, only: lpr, mynum

  real(kind=dp), intent(in)                 :: Q, f_min, f_max
  integer, intent(in)                       :: N, nfsamp, max_it

  real(kind=dp), optional, intent(in)       :: Tw, Ty, d, w_ref, alpha
  !f2py real(kind=dp), optional, intent(in) :: Tw = .1, Ty = .1, d = .99995
  !f2py real(kind=dp), optional, intent(in) :: w_ref = 1., alpha = 0
  real(kind=dp)                             :: Tw_loc = .1, Ty_loc = .1
  real(kind=dp)                             :: d_loc = .99995, w_ref_loc = 1.
  real(kind=dp)                             :: alpha_loc = 0

  logical, optional, intent(in)             :: fixfreq, verbose, exact, freq_weight
  !f2py logical, optional, intent(in)       :: fixfreq = 0, verbose = 0
  !f2py logical, optional, intent(in)       :: exact = 0, freq_weight = 0
  logical                                   :: fixfreq_loc = .false., verbose_loc = .false.
  logical                                   :: exact_loc = .false., freq_weight_loc = .true.

  real(kind=dp), intent(out)    :: w_j(N)
  real(kind=dp), intent(out)    :: y_j(N)
  real(kind=dp), intent(out)    :: w(nfsamp)
  real(kind=dp), intent(out)    :: q_fit(nfsamp)
  real(kind=dp), intent(out)    :: chil(max_it)

  real(kind=dp)     :: w_j_test(N)
  real(kind=dp)     :: y_j_test(N)
  real(kind=dp)     :: expo
  real(kind=dp)     :: chi, chi_test
  real(kind=dp)     :: randnr
  real(kind=dp)     :: Q_target(nfsamp)
  real(kind=dp)     :: weights(nfsamp)

  integer           :: j, it, last_it_print

  ! set default values
  if (present(Tw)) Tw_loc = Tw
  if (present(Ty)) Ty_loc = Ty
  if (present(d)) d_loc = d
  if (present(w_ref)) w_ref_loc = w_ref
  if (present(alpha)) alpha_loc = alpha

  if (present(fixfreq)) fixfreq_loc = fixfreq
  if (present(verbose)) verbose_loc = verbose
  if (present(exact)) exact_loc = exact
  if (present(freq_weight)) freq_weight_loc = freq_weight

  if (.not. lpr) verbose_loc = .false.

  ! Set the starting test frequencies equally spaced in log frequency
  if (N > 1) then
     expo = (log10(f_max) - log10(f_min)) / (N - 1.d0)
     do j=1, N
        w_j_test(j) = 2 * pi * 10**(log10(f_min) + (j - 1) * expo)
     enddo
  else
     w_j_test(1) = (f_max * f_min)**.5 * 2 * pi
  endif

  if (verbose_loc) print *, w_j_test

  ! Set the sampling frequencies equally spaced in log frequency
  expo = (log10(f_max) - log10(f_min)) / (nfsamp - 1.d0)
  do j=1, nfsamp
     w(j) = 2 * pi * 10**(log10(f_min) + (j - 1) * expo)
  enddo

  if (verbose_loc) print *, w

  ! compute target Q from power law
  Q_target(:) = Q * (w / w_ref_loc) ** alpha_loc

  ! compute weights for linear frequency weighting
  if (freq_weight_loc) then
     weights = w / sum(w) * nfsamp
  else
     weights(:) = 1
  endif

  ! initial weights y_j based on an empirical guess
  y_j_test = 1.d0 / Q * 1.5
  if (verbose_loc) print *, y_j_test

  ! initial Q(omega)
  call q_linear_solid(y_j=y_j_test, w_j=w_j_test, w=w, exact=exact_loc, Qls=q_fit)

  if (verbose_loc) print *, q_fit

  ! initial chi
  call l2_error(Q_target=Q_target, Qls=q_fit, weights=weights, lse=chi)
  if (verbose_loc) print *, 'initital chi: ', chi

  y_j(:) = y_j_test(:)
  w_j(:) = w_j_test(:)

  last_it_print = -1
  ! actuall simulated annealing loop:
  do it=1, max_it
     do j=1, N
        call random_number(randnr)
        if (.not. fixfreq_loc) &
           w_j_test(j) = w_j(j) * (1.0 + (0.5 - randnr) * Tw_loc)

        call random_number(randnr)
        y_j_test(j) = y_j(j) * (1.0 + (0.5 - randnr) * Ty_loc)
     enddo

     ! compute Q with test parameters
     call q_linear_solid(y_j=y_j_test, w_j=w_j_test, w=w, exact=exact_loc, Qls=q_fit)

     ! compute new misfit and new temperature
     call l2_error(Q_target=Q_target, Qls=q_fit, weights=weights, lse=chi_test)
     Tw_loc = Tw_loc * d_loc
     Ty_loc = Ty_loc * d_loc

     ! check if the tested parameters are better, if so, update
     if (chi_test < chi) then
        y_j(:) = y_j_test(:)
        w_j(:) = w_j_test(:)
        chi = chi_test
     endif

     chil(it) = chi
  enddo

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
