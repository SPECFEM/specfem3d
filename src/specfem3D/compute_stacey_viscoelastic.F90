!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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


  use specfem_par, only: SAVE_STACEY
  use specfem_par_elastic, only: displ

  ! boundary coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE,RECIPROCITY_AND_KH_INTEGRAL
  use specfem_par, only: it_dsm, it_fk, Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, Tract_axisem_time
  ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  use specfem_par_coupling, only: npt,nbdglb, &
     VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, NP_RESAMP, &
     vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK

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


! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
! FK surface
  integer :: ipt,ixglob
  integer :: ii, kk, iim1, iip1, iip2
  double precision :: cs(4), w
  real(kind=CUSTOM_REAL) ::  cs_single(4) !vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
! *********************************************************************************

  !! CD modif. : begin (implemented by VM) !! For coupling with DSM
  integer :: kaxisem, ip
  integer,dimension(NGLLSQUARE,num_abs_boundary_faces) :: ipt_table

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! injecting boundary
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    ! note: only applies boundary condition in first phase
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
      if (old_DSM_coupling_from_Vadim) then
        if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
          call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
        endif
      else
        !! MODIFS DE NOBU 2D
      endif

    else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
      call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)

      !! CD CD add this
      if (RECIPROCITY_AND_KH_INTEGRAL) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)

    else if ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
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

      do ip=1, npt
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
      it_fk=it_fk+1

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
    endif
  endif ! COUPLE_WITH_INJECTION_TECHNIQUE

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,igll,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,jacobianw, &
!$OMP kaxisem,ixglob,ipt)
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

        !! CD CD !! For coupling with EXTERNAL CODE
        if (COUPLE_WITH_INJECTION_TECHNIQUE) then
          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then  !! To verify for NOBU version
            vx = vx - Veloc_dsm_boundary(1,it_dsm,igll,iface)
            vy = vy - Veloc_dsm_boundary(2,it_dsm,igll,iface)
            vz = vz - Veloc_dsm_boundary(3,it_dsm,igll,iface)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then !! VM VM add AxiSEM
            kaxisem = igll + NGLLSQUARE*(iface - 1)
            vx = vx - Veloc_axisem(1,kaxisem)
            vy = vy - Veloc_axisem(2,kaxisem)
            vz = vz - Veloc_axisem(3,kaxisem)
          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
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

            vx = vx - vx_FK(ipt)
            vy = vy - vy_FK(ipt)
            vz = vz - vz_FK(ipt)
          endif
        endif

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

        !! CD CD !! For coupling with DSM
        if (COUPLE_WITH_INJECTION_TECHNIQUE) then
          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then    !! To verify for NOBU version
            tx = tx - Tract_dsm_boundary(1,it_dsm,igll,iface)
            ty = ty - Tract_dsm_boundary(2,it_dsm,igll,iface)
            tz = tz - Tract_dsm_boundary(3,it_dsm,igll,iface)
          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            tx = tx - Tract_axisem(1,kaxisem)
            ty = ty - Tract_axisem(2,kaxisem)
            tz = tz - Tract_axisem(3,kaxisem)
          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_FK) then
            ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
            tx = tx - tx_FK(ipt)
            ty = ty - ty_FK(ipt)
            tz = tz - tz_FK(ipt)
          endif
        endif

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! adds stacey term (weak form)
!$OMP ATOMIC
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
!$OMP ATOMIC
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
!$OMP ATOMIC
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! adjoint simulations
        if (SAVE_STACEY) then
          b_absorb_field(1,igll,iface) = tx*jacobianw
          b_absorb_field(2,igll,iface) = ty*jacobianw
          b_absorb_field(3,igll,iface) = tz*jacobianw
        endif !adjoint

      enddo
    endif ! ispec_is_elastic
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

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

  ! adjoint simulations: stores absorbed wavefield part
  if (SAVE_STACEY) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
    if (iphase == 1) it_dsm = it_dsm + 1
    !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  endif

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

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
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

!=============================================================================
!
  !! VM VM & CD CD !! For coupling with external code
!
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
    do j=1,NGLLY
      do i=1,NGLLX
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

  integer nb
  real(kind=CUSTOM_REAL) :: Veloc_axisem(3,nb)
  real(kind=CUSTOM_REAL) :: Tract_axisem(3,nb)

  read(IIN_veloc_dsm) Veloc_axisem, Tract_axisem

  end subroutine read_axisem_file

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

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
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
                                             SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                                             b_num_abs_boundary_faces,b_reclen_field,b_absorb_field,Mesh_pointer, &
                                             FORWARD_OR_ADJOINT)

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,UNDO_ATTENUATION_AND_OR_PML

  implicit none

! communication overlap
  integer,intent(in) :: iphase

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces

! adjoint simulations
  integer,intent(in) :: SIMULATION_TYPE
  integer,intent(in) :: NSTEP,it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_field

  logical,intent(in) :: SAVE_FORWARD

  ! GPU_MODE variables
  integer(kind=8),intent(in) :: Mesh_pointer
  integer, intent(in) :: FORWARD_OR_ADJOINT

  !! For coupling with DSM

!! DK DK beware:: these two arrays are automatic arrays, not subroutine arguments, thus allocated every time
!! DK DK beware:: this routine is called, and erased when the routine exits; that is very strange
! real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
! real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  integer :: it_dsm

  it_dsm = 1 !! initialize dsm iterator
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    if (old_DSM_coupling_from_Vadim) then
      if (iphase == 1) then
        if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
          stop 'DK DK old_DSM_coupling_from_Vadim support is discontinued'
!         call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
        endif
      endif
    else
      !! MODIFS DE NOBU 2D
    endif
  endif

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

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

    ! adjoint simulations: stores absorbed wavefield part
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! writes out absorbing boundary value
      call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
    endif
  endif

  !! CD CD : begin
  !! For coupling with DSM
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then !! To verify for NOBU version
    if (iphase == 1) then
      it_dsm = it_dsm + 1
    endif
  endif
  !! CD CD : end

  end subroutine compute_stacey_viscoelastic_GPU

