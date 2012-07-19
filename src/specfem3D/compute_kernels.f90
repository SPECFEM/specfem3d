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

  subroutine compute_kernels()

! kernel calculations
! see e.g. Tromp et al. (2005) for elastic calculation
! and Morency et al. (2009) for poroelastic calculation

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  implicit none

  ! elastic simulations
  if( ELASTIC_SIMULATION ) then
    call compute_kernels_el()
  endif

  ! elastic simulations
  if( ACOUSTIC_SIMULATION ) then
    call compute_kernels_ac()
  endif

  ! poroelastic simulations
  if( POROELASTIC_SIMULATION ) then
    call compute_kernels_po()
  endif

  ! computes an approximative hessian for preconditioning kernels
  if ( APPROXIMATE_HESS_KL ) then
    call compute_kernels_hessian()
  endif

  end subroutine compute_kernels


!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_el()

! kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par
  use specfem_par_elastic

  implicit none
  ! local parameters
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc,b_epsilondev_loc

  if( .not. GPU_MODE ) then
    ! updates kernels on CPU
    do ispec = 1, NSPEC_AB

      ! elastic domains
      if( ispec_is_elastic(ispec) ) then

         do k = 1, NGLLZ
            do j = 1, NGLLY
               do i = 1, NGLLX
                  iglob = ibool(i,j,k,ispec)

                  ! isotropic kernels
                  ! note: takes displacement from backward/reconstructed (forward) field b_displ
                  !          and acceleration from adjoint field accel (containing adjoint sources)
                  !
                  ! note: : time integral summation uses deltat
                  !
                  ! compare with Tromp et al. (2005), eq. (14), which takes adjoint displacement
                  ! and forward acceleration, that is the symmetric form of what is calculated here
                  ! however, this kernel expression is symmetric with regards
                  ! to interchange adjoint - forward field
                  rho_kl(i,j,k,ispec) =  rho_kl(i,j,k,ispec) &
                       + deltat * dot_product(accel(:,iglob), b_displ(:,iglob))

                  ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
                  ! note: multiplication with 2*mu(x) will be done after the time loop
                  epsilondev_loc(1) = epsilondev_xx(i,j,k,ispec)
                  epsilondev_loc(2) = epsilondev_yy(i,j,k,ispec)
                  epsilondev_loc(3) = epsilondev_xy(i,j,k,ispec)
                  epsilondev_loc(4) = epsilondev_xz(i,j,k,ispec)
                  epsilondev_loc(5) = epsilondev_yz(i,j,k,ispec)

                  b_epsilondev_loc(1) = b_epsilondev_xx(i,j,k,ispec)
                  b_epsilondev_loc(2) = b_epsilondev_yy(i,j,k,ispec)
                  b_epsilondev_loc(3) = b_epsilondev_xy(i,j,k,ispec)
                  b_epsilondev_loc(4) = b_epsilondev_xz(i,j,k,ispec)
                  b_epsilondev_loc(5) = b_epsilondev_yz(i,j,k,ispec)

                  mu_kl(i,j,k,ispec) =  mu_kl(i,j,k,ispec) &
                       + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
                       + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
                       + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                       epsilondev_loc(5)*b_epsilondev_loc(5)) )

                  ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
                  ! note: multiplication with kappa(x) will be done after the time loop
                  kappa_kl(i,j,k,ispec) = kappa_kl(i,j,k,ispec) &
                       + deltat * (9 * epsilon_trace_over_3(i,j,k,ispec) &
                       * b_epsilon_trace_over_3(i,j,k,ispec))

               enddo
            enddo
         enddo
      endif !ispec_is_elastic

    enddo

  else
    ! updates kernels on GPU
    call compute_kernels_elastic_cuda(Mesh_pointer,deltat)
  endif


  ! moho kernel
  if( SAVE_MOHO_MESH ) then
    if( GPU_MODE ) then
      call transfer_accel_from_device(NDIM*NGLOB_AB,accel,Mesh_pointer)
      call transfer_b_displ_from_device(NDIM*NGLOB_AB,b_displ,Mesh_pointer)
    endif
    ! updates on CPU
    call compute_boundary_kernel()
  endif

  ! for noise simulations --- source strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    call compute_kernels_strength_noise(NGLLSQUARE*num_free_surface_faces,ibool, &
                        sigma_kl,displ,deltat,it, &
                        normal_x_noise,normal_y_noise,normal_z_noise, &
                        noise_surface_movie, &
                        NSPEC_AB,NGLOB_AB, &
                        num_free_surface_faces,free_surface_ispec,free_surface_ijk,&
                        GPU_MODE,Mesh_pointer)
  endif

  end subroutine compute_kernels_el

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par
  use specfem_par_acoustic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: b_displ_elm,accel_elm
  real(kind=CUSTOM_REAL) :: kappal,rhol
  integer :: i,j,k,ispec,iglob

  ! updates kernels on GPU
  if(GPU_MODE) then

    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_acoustic_cuda(Mesh_pointer,deltat)

    ! kernels are done
    return
  endif

  ! updates kernels
  do ispec = 1, NSPEC_AB

    ! acoustic domains
    if( ispec_is_acoustic(ispec) ) then

      ! backward fields: displacement vector
      call compute_gradient(ispec,NSPEC_ADJOINT,NGLOB_ADJOINT, &
                      b_potential_acoustic, b_displ_elm,&
                      hprime_xx,hprime_yy,hprime_zz, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      ibool,rhostore,GRAVITY)
      ! adjoint fields: acceleration vector
      ! new expression (\partial_t^2\bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger)
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                      potential_acoustic, accel_elm,&
                      hprime_xx,hprime_yy,hprime_zz, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      ibool,rhostore,GRAVITY)

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! new expression
            ! density kernel
            rhol = rhostore(i,j,k,ispec)
            rho_ac_kl(i,j,k,ispec) =  rho_ac_kl(i,j,k,ispec) &
                      + deltat * rhol * dot_product(accel_elm(:,i,j,k), b_displ_elm(:,i,j,k))

            ! bulk modulus kernel
            kappal = 1._CUSTOM_REAL / kappastore(i,j,k,ispec)
            kappa_ac_kl(i,j,k,ispec) = kappa_ac_kl(i,j,k,ispec) &
                                  + deltat * kappal  &
                                  * potential_acoustic(iglob) &
                                  * b_potential_dot_dot_acoustic(iglob)

          enddo
        enddo
      enddo
    endif ! ispec_is_acoustic

  enddo

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use specfem_par
  use specfem_par_poroelastic

  implicit none
  ! local parameters
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL), dimension(5) :: epsilonsdev_loc,b_epsilonsdev_loc

  if( .not. GPU_MODE ) then
    ! updates kernels on CPU
    do ispec = 1, NSPEC_AB

      ! poroelastic domains
      if( ispec_is_poroelastic(ispec) ) then

         do k = 1, NGLLZ
            do j = 1, NGLLY
               do i = 1, NGLLX
                  iglob = ibool(i,j,k,ispec)

                  ! 8 isotropic kernels
                  ! note: takes displacement from backward/reconstructed
                  ! (forward) field b_displ
                  !          and acceleration from adjoint field accel
                  !          (containing adjoint sources)
                  !
                  ! note: : time integral summation uses deltat
                  ! 3 density kernels, see Morency et al. (2009),
                  ! equations (39)-(41)
                  rhot_kl(i,j,k,ispec) =  rhot_kl(i,j,k,ispec) &
                       + deltat * dot_product(accels_poroelastic(:,iglob), b_displs_poroelastic(:,iglob))

                  rhof_kl(i,j,k,ispec) =  rhof_kl(i,j,k,ispec) &
                       + deltat * ( dot_product(accelw_poroelastic(:,iglob), b_displs_poroelastic(:,iglob)) &
                       + dot_product(accels_poroelastic(:,iglob), b_displw_poroelastic(:,iglob)) )

                  sm_kl(i,j,k,ispec) =  sm_kl(i,j,k,ispec) &
                       + deltat * dot_product(accelw_poroelastic(:,iglob), b_displw_poroelastic(:,iglob))

                  ! kernel for viscous damping, see e.g. Morency et al. (2009),
                  ! equation (42)
                  eta_kl(i,j,k,ispec) =  eta_kl(i,j,k,ispec) &
                       + deltat * dot_product(velocw_poroelastic(:,iglob), b_displw_poroelastic(:,iglob))

                  ! kernel for frame shear modulus, see e.g. Morency et al. (2009),
                  ! equation (46)
                  ! note: multiplication with 2*mufr(x) will be done after the
                  ! time loop
                  epsilonsdev_loc(1) = epsilonsdev_xx(i,j,k,ispec)
                  epsilonsdev_loc(2) = epsilonsdev_yy(i,j,k,ispec)
                  epsilonsdev_loc(3) = epsilonsdev_xy(i,j,k,ispec)
                  epsilonsdev_loc(4) = epsilonsdev_xz(i,j,k,ispec)
                  epsilonsdev_loc(5) = epsilonsdev_yz(i,j,k,ispec)

                  b_epsilonsdev_loc(1) = b_epsilonsdev_xx(i,j,k,ispec)
                  b_epsilonsdev_loc(2) = b_epsilonsdev_yy(i,j,k,ispec)
                  b_epsilonsdev_loc(3) = b_epsilonsdev_xy(i,j,k,ispec)
                  b_epsilonsdev_loc(4) = b_epsilonsdev_xz(i,j,k,ispec)
                  b_epsilonsdev_loc(5) = b_epsilonsdev_yz(i,j,k,ispec)

                  mufr_kl(i,j,k,ispec) =  mufr_kl(i,j,k,ispec) &
                       + deltat * (epsilonsdev_loc(1)*b_epsilonsdev_loc(1) + epsilonsdev_loc(2)*b_epsilonsdev_loc(2) &
                       + (epsilonsdev_loc(1)+epsilonsdev_loc(2)) * (b_epsilonsdev_loc(1)+b_epsilonsdev_loc(2)) &
                       + 2 * (epsilonsdev_loc(3)*b_epsilonsdev_loc(3) + epsilonsdev_loc(4)*b_epsilonsdev_loc(4) + &
                       epsilonsdev_loc(5)*b_epsilonsdev_loc(5)) )

                  ! 3 kernels for the bulk moduli, see e.g. Morency et al. (2009),
                  ! equations (43)-(45)
                  ! note: multiplication by respective moduli will be done after the
                  ! time loop
                  B_kl(i,j,k,ispec) = B_kl(i,j,k,ispec) &
                       + deltat * (9 * epsilons_trace_over_3(i,j,k,ispec) &
                       * b_epsilons_trace_over_3(i,j,k,ispec))

                  C_kl(i,j,k,ispec) = C_kl(i,j,k,ispec) &
                       + deltat * (9 * epsilons_trace_over_3(i,j,k,ispec) &
                       * b_epsilonw_trace_over_3(i,j,k,ispec) &
                       + 9 * epsilonw_trace_over_3(i,j,k,ispec) &
                       * b_epsilons_trace_over_3(i,j,k,ispec) )

                  M_kl(i,j,k,ispec) = M_kl(i,j,k,ispec) &
                       + deltat * (9 * epsilonw_trace_over_3(i,j,k,ispec) &
                       * b_epsilonw_trace_over_3(i,j,k,ispec))
               enddo
            enddo
         enddo
      endif !ispec_is_poroelastic

    enddo

  else
    ! TO DO
    ! updates kernels on GPU
    !call compute_kernels_poroelastic_cuda(Mesh_pointer,deltat)
  endif

  end subroutine compute_kernels_po
!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_hessian()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ):: b_accel_elm,accel_elm
  integer :: i,j,k,ispec,iglob

  ! updates kernels on GPU
  if(GPU_MODE) then

    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_hess_cuda(Mesh_pointer,deltat, &
                                  ELASTIC_SIMULATION,ACOUSTIC_SIMULATION)

    ! done on GPU
    return
  endif

  ! loops over all elements
  do ispec = 1, NSPEC_AB

    ! acoustic domains
    if( ispec_is_acoustic(ispec) ) then

      ! adjoint fields: acceleration vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                      potential_dot_dot_acoustic, accel_elm,&
                      hprime_xx,hprime_yy,hprime_zz, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      ibool,rhostore,GRAVITY)

      ! adjoint fields: acceleration vector
      call compute_gradient(ispec,NSPEC_AB,NGLOB_AB, &
                      b_potential_dot_dot_acoustic, b_accel_elm,&
                      hprime_xx,hprime_yy,hprime_zz, &
                      xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                      ibool,rhostore,GRAVITY)

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! approximates hessian
            ! term with adjoint acceleration and backward/reconstructed acceleration
            hess_ac_kl(i,j,k,ispec) =  hess_ac_kl(i,j,k,ispec) &
               + deltat * dot_product(accel_elm(:,i,j,k), b_accel_elm(:,i,j,k))

          enddo
        enddo
      enddo
    endif

    ! elastic domains
    if( ispec_is_elastic(ispec) ) then
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)

            ! approximates hessian
            ! term with adjoint acceleration and backward/reconstructed acceleration
            hess_kl(i,j,k,ispec) =  hess_kl(i,j,k,ispec) &
               + deltat * dot_product(accel(:,iglob), b_accel(:,iglob))

          enddo
        enddo
      enddo
    endif

  enddo

  end subroutine compute_kernels_hessian


