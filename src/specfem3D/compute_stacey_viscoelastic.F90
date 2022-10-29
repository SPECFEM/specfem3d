!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_forward(NSPEC_AB,NGLOB_AB,accel, &
                                                 ibool,iphase, &
                                                 abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                 abs_boundary_ijk,abs_boundary_ispec, &
                                                 num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                                                 ispec_is_elastic, &
                                                 it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  ! Kirchoff-Helmholtz integrals
  use specfem_par_elastic, only: displ

  ! wavefield injection
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! injecting boundary wavefield
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then
    ! adds boundary contribution from injected wavefield
    call compute_coupled_injection_contribution(NSPEC_AB,NGLOB_AB,accel, &
                                                ibool,iphase, &
                                                abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                abs_boundary_ijk,abs_boundary_ispec, &
                                                num_abs_boundary_faces,rho_vp,rho_vs, &
                                                ispec_is_elastic, &
                                                it)
  endif ! COUPLE_WITH_INJECTION_TECHNIQUE

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,igll,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,jacobianw)
!$OMP DO
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        vx = veloc(1,iglob)
        vy = veloc(2,iglob)
        vz = veloc(3,iglob)

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz

        ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! adds stacey term (weak form)
!$OMP ATOMIC
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
!$OMP ATOMIC
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
!$OMP ATOMIC
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! for kernel simulations
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_absorb_field(1,igll,iface) = tx*jacobianw
          b_absorb_field(2,igll,iface) = ty*jacobianw
          b_absorb_field(3,igll,iface) = tz*jacobianw
        endif

      enddo
    endif ! ispec_is_elastic
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  ! for kernel simulations: stores absorbed wavefield part
  if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
    ! adds boundary injection contribution to stacey contribution before saving to disk
    ! this avoids storing the boundary injection arrays as a separate file
    ! for kernel simulations to reconstruct forward wavefields.
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      b_absorb_field(:,:,:) = b_absorb_field(:,:,:) + b_boundary_injection_field(:,:,:)
    endif

    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  ! Kirchoff-Helmholtz integrals
  !! CD CD added this
  if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      if (ispec_is_elastic(ispec)) then
        ! reference GLL points on boundary face
        do igll = 1,NGLLSQUARE
          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)
          iglob = ibool(i,j,k,ispec)
          write(237) b_absorb_field(1,igll,iface), b_absorb_field(2,igll,iface), b_absorb_field(3,igll,iface)
          write(238) displ(1,iglob), displ(2,iglob), displ(3,iglob)
        enddo
      endif
    enddo
  endif

  ! wavefield injection
  ! not used and implemented yet..
  !if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
  !  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
  !     if (iphase == 1) it_dsm = it_dsm + 1
  !     !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  !  endif
  !endif

  end subroutine compute_stacey_viscoelastic_forward

!
!=====================================================================
!

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_backward(NSPEC_AB, &
                                                  ibool,iphase, &
                                                  abs_boundary_ijk,abs_boundary_ispec, &
                                                  num_abs_boundary_faces, &
                                                  ispec_is_elastic,SIMULATION_TYPE, &
                                                  NSTEP,it,NGLOB_ADJOINT,b_accel, &
                                                  b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants
  use specfem_par, only: myrank

  implicit none

  integer,intent(in) :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: SIMULATION_TYPE
  integer,intent(in) :: NSTEP,it,NGLOB_ADJOINT
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT),intent(inout) :: b_accel

  ! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  ! reads in absorbing boundary array (when first phase is running)
  ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
  call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        ! adjoint simulations
        b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward

!
!=====================================================================
!

  subroutine compute_stacey_viscoelastic_backward_undoatt(NSPEC_AB,NGLOB_AB,b_accel,b_veloc, &
                                                          ibool,iphase, &
                                                          abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                          abs_boundary_ijk,abs_boundary_ispec, &
                                                          num_abs_boundary_faces, &
                                                          rho_vp,rho_vs,ispec_is_elastic)

  use constants
  use specfem_par, only: myrank,SIMULATION_TYPE

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: b_accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(in) :: b_veloc

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_backward() with wrong SIMULATION_TYPE')

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Engquist)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_elastic(ispec)) then
      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets velocity
        iglob = ibool(i,j,k,ispec)

        vx = b_veloc(1,iglob)
        vy = b_veloc(2,iglob)
        vz = b_veloc(3,iglob)

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz

        ! stacey term: velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! adds stacey term (weak form)
        b_accel(1,iglob) = b_accel(1,iglob) - tx*jacobianw
        b_accel(2,iglob) = b_accel(2,iglob) - ty*jacobianw
        b_accel(3,iglob) = b_accel(3,iglob) - tz*jacobianw
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_stacey_viscoelastic_backward_undoatt

!
!=====================================================================
!

! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                                             NSTEP,it, &
                                             b_num_abs_boundary_faces,b_reclen_field,b_absorb_field,Mesh_pointer, &
                                             FORWARD_OR_ADJOINT)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,UNDO_ATTENUATION_AND_OR_PML

  ! wavefield injection
  use specfem_par, only: NSPEC_AB,NGLOB_AB,ibool, &
                         abs_boundary_normal,abs_boundary_jacobian2Dw, &
                         abs_boundary_ijk,abs_boundary_ispec
  use specfem_par_elastic, only: accel,rho_vp,rho_vs,ispec_is_elastic
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none

  ! communication overlap
  integer,intent(in) :: iphase

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces

  ! adjoint simulations
  integer,intent(in) :: NSTEP,it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  ! GPU_MODE variables
  integer(kind=8),intent(in) :: Mesh_pointer
  integer, intent(in) :: FORWARD_OR_ADJOINT

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! wavefield injection
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then
    ! note: wavefield injection arrays are only available on CPU
    !
    !       as a quick work-around, we transfer here the velocity and acceleration arrays between GPU-CPU and vice versa,
    !       until the full injection contribution will be implemented on the GPU side as CUDA kernels.
    !       this is slowing down the simulation a bit.
    !
    ! transfers acceleration to the CPU
    call transfer_accel_from_device(NDIM*NGLOB_AB,accel, Mesh_pointer)

    ! adds boundary contribution from injected wavefield
    call compute_coupled_injection_contribution(NSPEC_AB,NGLOB_AB,accel, &
                                                ibool,iphase, &
                                                abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                abs_boundary_ijk,abs_boundary_ispec, &
                                                num_abs_boundary_faces,rho_vp,rho_vs, &
                                                ispec_is_elastic, &
                                                it)

    ! transfers updated acceleration field back to the GPU
    call transfer_accel_to_device(NDIM*NGLOB_AB,accel, Mesh_pointer)
  endif ! COUPLE_WITH_INJECTION_TECHNIQUE

  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! no need to store boundaries on disk
    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_viscoelastic_undoatt_cuda(Mesh_pointer,iphase,FORWARD_OR_ADJOINT)
  else
    ! adjoint simulations:
    if (SIMULATION_TYPE == 3) then
      ! reads in absorbing boundary array (when first phase is running)
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
    endif !adjoint

    call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_field,FORWARD_OR_ADJOINT)

    ! for kernel simulations: stores absorbed wavefield part
    if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
      ! adds boundary injection contribution to stacey contribution before saving to disk
      ! this avoids storing the boundary injection arrays as a separate file
      ! for kernel simulations to reconstruct forward wavefields.
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        b_absorb_field(:,:,:) = b_absorb_field(:,:,:) + b_boundary_injection_field(:,:,:)
      endif

      ! writes out absorbing boundary value
      call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
    endif
  endif

  ! wavefield injection
  ! not used and implemented yet..
  !if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
  !  if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
  !     if (iphase == 1) it_dsm = it_dsm + 1
  !     !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  !  endif
  !endif

  end subroutine compute_stacey_viscoelastic_GPU


!=============================================================================
!
! For coupling with external code
!
!=============================================================================

  subroutine compute_coupled_injection_contribution(NSPEC_AB,NGLOB_AB,accel, &
                                                    ibool,iphase, &
                                                    abs_boundary_normal,abs_boundary_jacobian2Dw, &
                                                    abs_boundary_ijk,abs_boundary_ispec, &
                                                    num_abs_boundary_faces,rho_vp,rho_vs, &
                                                    ispec_is_elastic, &
                                                    it)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  ! boundary coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE,RECIPROCITY_AND_KH_INTEGRAL
  use specfem_par_coupling, only: it_dsm, it_fk, &
    Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, Tract_axisem_time
  ! FK3D calculation
  use specfem_par_coupling, only: npt, nbdglb, &
    VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, NP_RESAMP, &
    vx_FK, vy_FK, vz_FK, tx_FK, ty_FK, tz_FK
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_field

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  ! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_elastic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: it

  ! local parameters
  real(kind=CUSTOM_REAL) :: vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll
  ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  ! FK surface
  integer :: ipt,ixglob
  integer :: ii, kk, iim1, iip1, iip2
  double precision :: cs(4), w
  real(kind=CUSTOM_REAL) ::  cs_single(4) !vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
  !! CD modif. : begin (implemented by VM) !! For coupling with DSM
  integer :: kaxisem, ip
  integer,dimension(NGLLSQUARE,num_abs_boundary_faces) :: ipt_table

!! comment from Vadim Monteiller, Feb 2017:

! txxbd est calcule dans add_to_compute_stacey_viscoelastic_1.F90 qui lui-meme appelle une subroutine qui se trouve
! dans add_to_compute_stacey_viscoelastic_11.F90. En fait je ne peux pas directement stocker txxbd en memoire
! (sinon il faut 500Go de ram, mais j'ai pas ca). Donc tous les 100 pas de temps je fais une fft
! et je prends la partie du sismo qui m'interesse. C'est a ce moment qu txxbd est rempli. C'est le call suivant qui fait ca:
!
! call store_next_FK_solution( VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f, &
!                    WKS_CMPLX_FOR_FFT, WKS_REAL_FOR_FFT, NF_FOR_STORING, &
!                    NF_FOR_FFT, NTIME_BETWEEN_FFT, NPOW_FOR_FFT, &
!                    vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, npt, it, deltat)

! je prends la solution en frequence : VX_f, VY_f, VZ_f, TX_f, TY_f, TZ_f
! et je la sors en temps : vxbd, vybd, vzbd, txxbd, tyybd, tzzbd, pour les 100 pas de temps suivants.

! la subroutine store_next_FK_solution est dans add_to_compute_stacey_viscoelastic_11.F90 ligne 882

! on stocke directement la traction (on fait toujours ca avec DSM et AxiSEM aussi),
! ca evite de stocker 6 composantes de stress, surtout qu'on a des problemes de memoire.

  ! safety checks
  if (.not. COUPLE_WITH_INJECTION_TECHNIQUE) return

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! only for forward wavefield
  if (SIMULATION_TYPE /= 1) return

  ! injecting boundary
  !
  ! note: Following Tong et al. (2014, GJI, 197 (1); https://doi.org/10.1093/gji/ggt508) or
  !       Tong et al. (2014, GRL, 41; https://doi.org/10.1002/2014GL061644), coupling the wavefield
  !       by injecting the displacement and stresses at the mesh boundary can be merged together
  !       with the Clayton-Engquist absorbing boundary by:
  !            ( T_total - T_FK ) * n = - rho alpha [ n * d\dt(u_total - u_FK)]n - rho beta [ t * d/dt(u_total - u_FK)]t
  !       see equation (1) in Tong et al.'s GJI paper; equivalently using velocities v:
  !            ( T_total - T_FK ) * n = - rho alpha [ n * (v_total - v_FK)]n - rho beta [ t * (v_total - v_FK)]t
  !        in this case, the velocity rather than the displacement is needed and the Clayton-Engquist formulation modified.
  !        this allows to absorb the outgoing wavefield components, while still injecting an external wavefield.
  !
  !       -> option A: (original implementation) use the above modified Clayton-Engquist formulation to inject the wavefield
  !                    and absorb outgoing components
  !
  !       Changing the Clayton-Engquist boundary expression however would require to inject the wavefield with a
  !       different formulation. thus, instead of using the modified equation (1) from above, an equivalent way
  !       would be to inject the displacement and stress from the wavefield before the absorbing boundary gets computed,
  !       separating absorbing boundary implementations from the wavefield injection.
  !
  !       -> option B: 1. take only T_FK as a right-hand-side term to update acceleration
  !                    2. apply the absorbing boundary stresses based on the injected velocity:
  !                         T_abs_inj * n =  - rho alpha [ n * (- v_FK)]n - rho beta [ t * (- v_FK)]t
  !                    3. add the final stress as an additional right-hand-side term due to the injection:
  !                         T_final * n = (- T_FK + T_abs_inj) * n
  !                    this might help to separate absorbing boundary and wavefield injection, and could also be slightly
  !                    faster in case this coupling is not used as we avoid additional if-statements within the loops.
  !
  !       we now implement option B, as this also allows an easier integration for GPU simulations.
  !       thus, coupling with wavefield injection now is supported for both GPU and kernel simulations.
  !
  !       however, coupling boundary wavefields still requires to have Stacey absorbing boundaries set as well.
  !       this is due to the fact that we only call this routine within the Stacey routine,
  !       and also because we store the coupling contribution together with the Stacey ones for reconstructing the wavefield
  !       in kernels simulations.

  ! gets velocity & stress for boundary points
  select case(INJECTION_TECHNIQUE_TYPE)
  case (INJECTION_TECHNIQUE_IS_DSM)
    ! DSM coupling
    if (old_DSM_coupling_from_Vadim) then
      if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
        call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
      endif
    else
      !! MODIFS DE NOBU 2D
    endif

  case (INJECTION_TECHNIQUE_IS_AXISEM)
    ! AxiSEM coupling
    call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)

    !! CD CD add this
    if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)

  case (INJECTION_TECHNIQUE_IS_FK)
    ! FK coupling
    !! find indices
    ! example: np_resamp = 2 and it = 1,2,3,4,5,6, ..
    !          --> ii = 1,1,2,2,3,3,..
    ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
    ! example: --> kk = 1,2,1,2,1,2,,..
    kk = it - (ii-1) * NP_RESAMP

    w = dble(kk-1) / dble(NP_RESAMP)

    ! Cubic spline values
    cs(4) = w*w*w/6.d0
    cs(1) = 1.d0/6.d0+w*(w-1.d0)/2.d0-cs(4)
    cs(3) = w+cs(1)-2.d0*cs(4)
    cs(2) = 1.d0-cs(1)-cs(3)-cs(4)

    cs_single(:) = sngl(cs(:))

    iim1 = ii-1
    iip1 = ii+1
    iip2 = ii+2

    do ip = 1, npt
       vx_FK(ip) = cs_single(1)* VX_t(ip,iim1) + cs_single(2)* VX_t(ip,ii) + cs_single(3)* VX_t(ip,iip1) + &
            cs_single(4)* VX_t(ip,iip2)
       vy_FK(ip) = cs_single(1)* VY_t(ip,iim1) + cs_single(2)* VY_t(ip,ii) + cs_single(3)* VY_t(ip,iip1) + &
            cs_single(4)* VY_t(ip,iip2)
       vz_FK(ip) = cs_single(1)* VZ_t(ip,iim1) + cs_single(2)* VZ_t(ip,ii) + cs_single(3)* VZ_t(ip,iip1) + &
            cs_single(4)* VZ_t(ip,iip2)
       tx_FK(ip) = cs_single(1)* TX_t(ip,iim1) + cs_single(2)* TX_t(ip,ii) + cs_single(3)* TX_t(ip,iip1) + &
            cs_single(4)* TX_t(ip,iip2)
       ty_FK(ip) = cs_single(1)* TY_t(ip,iim1) + cs_single(2)* TY_t(ip,ii) + cs_single(3)* TY_t(ip,iip1) + &
            cs_single(4)* TY_t(ip,iip2)
       tz_FK(ip) = cs_single(1)* TZ_t(ip,iim1) + cs_single(2)* TZ_t(ip,ii) + cs_single(3)* TZ_t(ip,iip1) + &
            cs_single(4)* TZ_t(ip,iip2)
    enddo
    it_fk = it_fk+1

    ! prepares ipt table
    ! (needed also for openmp threads accessing correct ipt index)
    ipt = 0
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      if (ispec_is_elastic(ispec)) then
        do igll = 1,NGLLSQUARE
          ! counter
          ipt = ipt + 1
          ipt_table(igll,iface) = ipt
        enddo
      endif
    enddo
  end select

  ! wavefield injection
  ! option B: with injection of velocity/stress before the absorbing boundary condition
  !           as separate terms to the right-hand-side
  do iface = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)
    if (ispec_is_elastic(ispec)) then
      ! GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! global index
        iglob = ibool(i,j,k,ispec)

        ! gets injected velocity & stress
        ! velocity and stresses would be subtracted from total, therefore we add a minus sign when getting the values
        select case(INJECTION_TECHNIQUE_TYPE)
        case (INJECTION_TECHNIQUE_IS_DSM)                 !! To verify for NOBU version
          ! velocity
          vx = - Veloc_dsm_boundary(1,it_dsm,igll,iface)
          vy = - Veloc_dsm_boundary(2,it_dsm,igll,iface)
          vz = - Veloc_dsm_boundary(3,it_dsm,igll,iface)
          ! stress
          tx = - Tract_dsm_boundary(1,it_dsm,igll,iface)
          ty = - Tract_dsm_boundary(2,it_dsm,igll,iface)
          tz = - Tract_dsm_boundary(3,it_dsm,igll,iface)
        case (INJECTION_TECHNIQUE_IS_AXISEM)              !! VM VM add AxiSEM
          kaxisem = igll + NGLLSQUARE*(iface - 1)
          ! velocity
          vx = - Veloc_axisem(1,kaxisem)
          vy = - Veloc_axisem(2,kaxisem)
          vz = - Veloc_axisem(3,kaxisem)
          ! stress
          tx = - Tract_axisem(1,kaxisem)
          ty = - Tract_axisem(2,kaxisem)
          tz = - Tract_axisem(3,kaxisem)
        case (INJECTION_TECHNIQUE_IS_FK)
          ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          ! point index using table lookup
          !ipt = ipt + 1
          ipt = ipt_table(igll,iface)
          !! DEBUGVM pour eviter de stocker pour profiler la vitesse de FK
          !vx_FK = vxbd(it_fk,ipt)
          !vy_FK = vybd(it_fk,ipt)
          !vz_FK = vzbd(it_fk,ipt)
          ! sanity check, make sure we are at the right point
          ixglob = nbdglb(ipt)
          !write(*,'(3(i10,1x),3x,  3i3, 3x, i10)') ipt, ixglob, iglob, i,j,k,ispec
          if (iglob /= ixglob) stop 'wrong boundary index for FK coupling'
          ! velocity
          vx = - vx_FK(ipt)
          vy = - vy_FK(ipt)
          vz = - vz_FK(ipt)
          ! stress
          tx = - tx_FK(ipt)
          ty = - ty_FK(ipt)
          tz = - tz_FK(ipt)
        end select

        ! computes absorbing boundary for injected velocity
        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! gets associated normal
        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        ! velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz

        ! adds stacey term to injected stresses:
        ! velocity vector component * vp * rho in normal direction + vs * rho component tangential to it
        tx = tx + rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(vx-vn*nx)
        ty = ty + rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(vy-vn*ny)
        tz = tz + rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(vz-vn*nz)

        ! adds final stress term for injected wavefield (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! for kernel simulations: stores contribution to buffer array and add it to stacey buffer before saving to disk
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_boundary_injection_field(1,igll,iface) = tx*jacobianw
          b_boundary_injection_field(2,igll,iface) = ty*jacobianw
          b_boundary_injection_field(3,igll,iface) = tz*jacobianw
        endif
      enddo
    endif ! ispec_is_elastic
  enddo

  end subroutine compute_coupled_injection_contribution

!=============================================================================

  subroutine read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)

  use constants

  implicit none

  integer igll,it_dsm
  integer iface,num_abs_boundary_faces,i,j
  real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  real(kind=CUSTOM_REAL) :: dsm_boundary_tmp(3,Ntime_step_dsm,NGLLX,NGLLY)

  it_dsm = 1

  do iface = 1,num_abs_boundary_faces
    igll = 0
    do j = 1,NGLLY
      do i = 1,NGLLX
        igll = igll + 1
        read(IIN_veloc_dsm) dsm_boundary_tmp(:,:,i,j)
        Veloc_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
        read(IIN_tract_dsm) dsm_boundary_tmp(:,:,i,j)
        Tract_dsm_boundary(:,:,igll,iface) = dsm_boundary_tmp(:,:,i,j)
      enddo
    enddo
  enddo

  end subroutine read_dsm_file

!=============================================================================

  subroutine read_axisem_file(Veloc_axisem,Tract_axisem,nb)

  use constants

  implicit none

  integer :: nb
  real(kind=CUSTOM_REAL) :: Veloc_axisem(3,nb)
  real(kind=CUSTOM_REAL) :: Tract_axisem(3,nb)

  read(IIN_veloc_dsm) Veloc_axisem, Tract_axisem

  end subroutine read_axisem_file

