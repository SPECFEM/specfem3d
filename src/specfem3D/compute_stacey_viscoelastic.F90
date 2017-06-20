!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic(NSPEC_AB,NGLOB_AB,accel, &
                        ibool,iphase, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces,veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        it,b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

  use specfem_par, only: it_dsm, it_fk, Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem, Tract_axisem_time

  use specfem_par_elastic, only: displ

  use shared_parameters, only: COUPLE_WITH_EXTERNAL_CODE, &
                  EXTERNAL_CODE_TYPE,EXTERNAL_CODE_IS_DSM,EXTERNAL_CODE_IS_AXISEM,EXTERNAL_CODE_IS_FK, &
                  old_DSM_coupling_from_Vadim,CUT_SOLUTION_FOR_VISU,SAVE_RUN_BOUN_FOR_KH_INTEGRAL,Ntime_step_dsm

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  use specfem_par_elastic, only: npt,nbdglb, NTIME_BETWEEN_FFT, &
     VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, NP_RESAMP, &
     vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
! *********************************************************************************

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

! Stacey conditions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: veloc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: rho_vp,rho_vs

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL) :: abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

! local parameters
  real(kind=CUSTOM_REAL) vx,vy,vz,nx,ny,nz,tx,ty,tz,vn,jacobianw
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


  if (COUPLE_WITH_EXTERNAL_CODE) then
     ipt = 0
    if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then

      if (old_DSM_coupling_from_Vadim) then
        if (iphase == 1) then
          if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
            call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
          endif
        endif
      else
        !! MODIFS DE NOBU 2D
      endif

   else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then

      if (iphase == 1) then
        call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)

        !! CD CD add this
        if (CUT_SOLUTION_FOR_VISU) Tract_axisem_time(:,:,it) = Tract_axisem(:,:)
     endif

  else if ( EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_FK) then
     if (iphase == 1) then

        !! find indices
        ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
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

     endif
  endif

endif

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

        vx = veloc(1,iglob)
        vy = veloc(2,iglob)
        vz = veloc(3,iglob)

          !! CD CD !! For coupling with EXTERNAL CODE
          if (COUPLE_WITH_EXTERNAL_CODE) then

            if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then  !! To verify for NOBU version
              vx = vx - Veloc_dsm_boundary(1,it_dsm,igll,iface)
              vy = vy - Veloc_dsm_boundary(2,it_dsm,igll,iface)
              vz = vz - Veloc_dsm_boundary(3,it_dsm,igll,iface)

            else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then !! VM VM add AxiSEM
                kaxisem = igll + NGLLSQUARE*(iface - 1)
                vx = vx - Veloc_axisem(1,kaxisem)
                vy = vy - Veloc_axisem(2,kaxisem)
                vz = vz - Veloc_axisem(3,kaxisem)
            endif

          endif

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          if (COUPLE_WITH_EXTERNAL_CODE .and. EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_FK) then
            ipt = ipt + 1

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
! *********************************************************************************

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
          if (COUPLE_WITH_EXTERNAL_CODE) then

            if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then    !! To verify for NOBU version
              tx = tx - Tract_dsm_boundary(1,it_dsm,igll,iface)
              ty = ty - Tract_dsm_boundary(2,it_dsm,igll,iface)
              tz = tz - Tract_dsm_boundary(3,it_dsm,igll,iface)

            else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then
                tx = tx - Tract_axisem(1,kaxisem)
                ty = ty - Tract_axisem(2,kaxisem)
                tz = tz - Tract_axisem(3,kaxisem)
            endif

          endif

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          if (COUPLE_WITH_EXTERNAL_CODE .and. EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_FK) then
            tx = tx - tx_FK(ipt)
            ty = ty - ty_FK(ipt)
            tz = tz - tz_FK(ipt)
          endif
! *********************************************************************************

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! adds stacey term (weak form)
        accel(1,iglob) = accel(1,iglob) - tx*jacobianw
        accel(2,iglob) = accel(2,iglob) - ty*jacobianw
        accel(3,iglob) = accel(3,iglob) - tz*jacobianw

        ! adjoint simulations
        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          b_absorb_field(1,igll,iface) = tx*jacobianw
          b_absorb_field(2,igll,iface) = ty*jacobianw
          b_absorb_field(3,igll,iface) = tz*jacobianw
        endif !adjoint

          !! CD CD added this
          if (SAVE_RUN_BOUN_FOR_KH_INTEGRAL) then
              write(237) b_absorb_field(1,igll,iface), b_absorb_field(2,igll,iface), b_absorb_field(3,igll,iface)
              write(238) displ(1,iglob), displ(2,iglob), displ(3,iglob)
          endif

      enddo
    endif ! ispec_is_elastic
  enddo

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  if (COUPLE_WITH_EXTERNAL_CODE) then !! To verify for NOBU version
    if (iphase == 1) it_dsm = it_dsm + 1
    !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here
  endif

  end subroutine compute_stacey_viscoelastic

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

  integer :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  integer :: iphase

  logical, dimension(NSPEC_AB) :: ispec_is_elastic

! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it,NGLOB_ADJOINT
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_ADJOINT):: b_accel

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

  use shared_parameters, only: Ntime_step_dsm,IIN_veloc_dsm,IIN_tract_dsm

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

  use shared_parameters, only: IIN_veloc_dsm

  implicit none

  integer nb
  real(kind=CUSTOM_REAL) :: Veloc_axisem(3,nb)
  real(kind=CUSTOM_REAL) :: Tract_axisem(3,nb)

  read(IIN_veloc_dsm) Veloc_axisem, Tract_axisem

  end subroutine read_axisem_file

!
!=====================================================================
!

! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(iphase,num_abs_boundary_faces, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field,Mesh_pointer)

  use constants

  use shared_parameters, only: COUPLE_WITH_EXTERNAL_CODE,old_DSM_coupling_from_Vadim,Ntime_step_dsm

  implicit none

! communication overlap
  integer :: iphase

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD

  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

  !! For coupling with DSM

!! DK DK beware:: these two arrays are automatic arrays, not subroutine arguments, thus allocated every time
!! DK DK beware:: this routine is called, and erased when the routine exits; that is very strange
! real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
! real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  integer :: it_dsm

  if (COUPLE_WITH_EXTERNAL_CODE) then
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

! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array (when first phase is running)
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
  endif !adjoint

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,iphase,b_absorb_field)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value
    call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
  endif

  !! CD CD : begin
  !! For coupling with DSM
  if (COUPLE_WITH_EXTERNAL_CODE) then !! To verify for NOBU version
    if (iphase == 1) then
      it_dsm = it_dsm + 1
    endif
  endif
  !! CD CD : end

  end subroutine compute_stacey_viscoelastic_GPU

!
!=====================================================================
!

! *********************************************************************************
! added by Ping Tong (TP / Tong Ping) for the FK3D calculation

!! count the number of point in the mesh partition boundary : npt
subroutine nbound(NSPEC_AB, num_abs_boundary_faces, abs_boundary_ispec, ispec_is_elastic, npt)

  use constants

  implicit none

  integer,                                    intent(inout)     :: npt
  integer,                                    intent(in)        :: NSPEC_AB, num_abs_boundary_faces
  ! elastic domain flag
  logical, dimension(NSPEC_AB),               intent(in)        :: ispec_is_elastic
  ! absorbing boundary surface
  integer, dimension(num_abs_boundary_faces), intent(in)        :: abs_boundary_ispec
  ! local parameters
  integer                                                       :: ispec, iface

  npt = 0
  do iface = 1, num_abs_boundary_faces
     ispec = abs_boundary_ispec(iface)
     if ( ispec_is_elastic(ispec) ) then
        ! reference GLL points on boundary face
        npt = npt + NGLLSQUARE
     endif
  enddo

end subroutine nbound

!==========================================================

subroutine FK3D(myrank, NSPEC_AB, ibool, abs_boundary_ijk, abs_boundary_normal, &
           abs_boundary_ispec, num_abs_boundary_faces, ispec_is_elastic, kpsv, nlayer, nstep, npt, nbdglb, &
           p, phi, xx0, yy0, zz0, tg, tt0, al_FK, be_FK, mu_FK, h_FK, deltat, &
           NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK)

  use constants

  use specfem_par, only: xstore,ystore,zstore,kappastore,mustore,rhostore
  use specfem_par_elastic, only: xx, yy, zz, xi1, xim, bdlambdamu, &
                                   nmx, nmy, nmz

  implicit none


  integer              :: NSPEC_AB,kpsv,nlayer,npt,nstep,ipt,iphase,myrank
  integer              :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP
  ! global index
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool
  integer, dimension(npt) :: nbdglb

  ! source
  real(kind=CUSTOM_REAL) :: p,phi,xx0,yy0,zz0,tt0
  real(kind=CUSTOM_REAL) :: DF_FK

  ! model
  real(kind=CUSTOM_REAL),dimension(nlayer) :: al_FK,be_FK,mu_FK,h_FK

  real(kind=CUSTOM_REAL) :: rhotmp,kappatmp,mutmp,xi,deltat,tg
  logical, dimension(NSPEC_AB) :: ispec_is_elastic
  logical ::  phase_is_inner

  ! absorbing boundary surface
  integer :: num_abs_boundary_faces
  integer :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer :: abs_boundary_ispec(num_abs_boundary_faces)
  real(kind=CUSTOM_REAL),dimension(3,NGLLSQUARE,num_abs_boundary_faces) :: abs_boundary_normal


  ! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll


  ! absorbs absorbing-boundary surface using Stacey condition (Clayton and Enquist)
  if (npt > 0) then
     allocate(xx(npt),yy(npt),zz(npt),xi1(npt),xim(npt),bdlambdamu(npt),nmx(npt),nmy(npt),nmz(npt))
  else
     allocate(xx(1),yy(1),zz(1),xi1(1),xim(1),bdlambdamu(1),nmx(1),nmy(1),nmz(1))
  endif

  nbdglb(:) = 0
  ipt = 0

  do iphase=1,2
     ! first for points on MPI interfaces
     if ( iphase == 1 ) then
        phase_is_inner = .false.
     else
        phase_is_inner = .true.
     endif

     do iface = 1,num_abs_boundary_faces
        ispec = abs_boundary_ispec(iface)
        if (iphase == 1) then
           if ( ispec_is_elastic(ispec) ) then

              ! reference GLL points on boundary face
              do igll = 1,NGLLSQUARE

                 ! gets local indices for GLL point
                 i = abs_boundary_ijk(1,igll,iface)
                 j = abs_boundary_ijk(2,igll,iface)
                 k = abs_boundary_ijk(3,igll,iface)

                 iglob = ibool(i,j,k,ispec)
                 ipt = ipt + 1
                 nbdglb(ipt) = iglob

                 xx(ipt) = xstore(iglob)
                 yy(ipt) = ystore(iglob)
                 zz(ipt) = zstore(iglob)

                 nmx(ipt) = abs_boundary_normal(1,igll,iface)
                 nmy(ipt) = abs_boundary_normal(2,igll,iface)
                 nmz(ipt) = abs_boundary_normal(3,igll,iface)


                 rhotmp   = rhostore(i,j,k,ispec)
                 kappatmp = kappastore(i,j,k,ispec)
                 mutmp    = mustore(i,j,k,ispec)

                 xi       = mutmp/(kappatmp+4.0*mutmp/3.0)
                 xi1(ipt) = 1.0-2.0*xi
                 xim(ipt) = (1.0-xi)*mutmp
                 bdlambdamu(ipt) = (3.0*kappatmp-2.0*mutmp)/(6.0*kappatmp+2.0*mutmp)

              enddo

           endif ! ispec_is_elastic
        endif ! ispec_is_inner

     enddo

  enddo


  call FK(al_FK, be_FK, mu_FK, h_FK, nlayer, tg, p, phi, xx0, yy0, zz0, tt0, deltat, nstep, npt, &
       kpsv,NF_FOR_STORING, NPOW_FOR_FFT,  NP_RESAMP, DF_FK, myrank)

  deallocate(xx, yy, zz, xi1, xim, bdlambdamu, nmx, nmy, nmz)

end subroutine FK3D

!==========================================================

  subroutine FK(al, be, mu, h, nlayer, Tg, p, phi, x0, y0, z0, t0, dt, npts, np, &
                kpsv,NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP, DF_FK, myrank)

    use constants
    use specfem_par_elastic, only: VX_t, VY_t, VZ_t, TX_t, TY_t, TZ_t, &
                                     xx, yy, zz, xi1, xim, bdlambdamu, &
                                     nmx, nmy, nmz, &
                                     NPTS_STORED, NPTS_INTERP

    implicit none

    integer,                parameter                         :: CUSTOM_CMPLX=8
    real(kind=CUSTOM_REAL), parameter                         :: zign_neg=-1.
    ! input and output
    integer,                                     intent(in)   :: nlayer, np, npts, kpsv, myrank
    integer                                                   :: NF_FOR_STORING, NPOW_FOR_FFT, NP_RESAMP

    ! model
    real(kind=CUSTOM_REAL),  dimension(nlayer),  intent(in)   :: al(nlayer),be(nlayer),mu(nlayer),H(nlayer)

    ! source
    real(kind=CUSTOM_REAL),                      intent(in)    :: dt, p, phi, x0, y0, z0, Tg, t0, DF_FK

    real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: fvec, dtmp
    complex(kind=CUSTOM_CMPLX), dimension(:,:), allocatable    :: coeff, field_f
    real(kind=CUSTOM_REAL),     dimension(:,:), allocatable    :: field
    complex(kind=CUSTOM_CMPLX), dimension(:),   allocatable    :: tmp_f1, tmp_f2, tmp_f3
    real(kind=CUSTOM_REAL),     dimension(:),   allocatable    :: tmp_t1, tmp_t2, tmp_t3, tmp_it1

    complex(kind=CUSTOM_CMPLX)                                 :: C_3,stf_coeff,a,b,c,d,delta_mat,N_mat(4,4),dx_f,dz_f,txz_f,tzz_f
    real(kind=CUSTOM_REAL)                                     :: epsil,dt_fk
    real(kind=CUSTOM_REAL)                                     :: sigma_rr,sigma_rt,sigma_rz,sigma_tt,sigma_tz,sigma_zz
    real(kind=CUSTOM_REAL)                                     ::  Txx_tmp, Txy_tmp, Txz_tmp, Tyy_tmp, Tyz_tmp, Tzz_tmp
    real(kind=CUSTOM_REAL)                                     :: df,om,tdelay,eta_p,eta_s,fmax,C_1
    integer                                                    :: lnpts2,npts2,nf,nf2,nn,ii,ip,j,nvar,lpts
    integer                                                    :: npoints2
    integer                                                    :: ier
    logical                                                    :: comp_stress, pout

    epsil = 1.0e-7
    comp_stress=.true.
    nvar=5
    pout =.false.

    fmax=1/(2*dt)    ! Nyquist frequency of specfem time serie

    !! new way to do time domain resampling
    df    = df_fk
    nf2   = NF_FOR_STORING+1   ! number of positive frequency sample points
    nf    = 2*NF_FOR_STORING   ! number of total frequencies after symetrisation
    npts2 = nf                 ! number of samples in time serie

    !! VM VM recompute new values for new way to do
    lnpts2=ceiling(log(npts2*1.)/log(2.))
    npts2=2**lnpts2
    NPOW_FOR_FFT = lnpts2

    dt_fk=1./(df*(npts2-1))
    if (myrank == 0) write(*,*) ' DT_FK ', dt_fk, ' df ' , df, 'npts2 ',npts2, 'lnpts2 ', lnpts2

    !! number of points for resmpled vector
    npoints2 = NP_RESAMP*(npts2-1)+1

    if (myrank == 0) then
       print *
       print *, ' entering the FK synthetics program .... '
       print *, '   model = ', nlayer
       print *, '   source = ',x0, y0, z0
       print *, '   azimuth ', phi
       print *, '   number of points used for FFT = ', npts2
       print *, '   total time length = ',t0,dt_fk,t0+(npts2-1)*dt_fk
       print *,  '  Number of samples stored for FK solution ',  NF_FOR_STORING
    endif


    allocate(fvec(nf2),stat=ier)
    fvec=0.
    do ii = 1, nf2
       fvec(ii)=(ii-1)*df
    enddo


    allocate(coeff(2,nf2),stat=ier)
    if (ier /= 0) stop 'error while allocating'

    allocate(field_f(nf,nvar),stat=ier)
    if (ier /= 0) stop 'error while allocating'

    allocate(field(npts2,nvar),dtmp(npts),stat=ier)
    if (ier /= 0) stop 'error while allocating'

    !! allocate debug vectors
    allocate(tmp_f1(npts2), tmp_f2(npts2), tmp_f3(npts2),stat=ier)

    if (ier /= 0) stop 'error while allocating'
    allocate(tmp_t1(npts2), tmp_t2(npts2), tmp_t3(npts2),stat=ier)

    if (ier /= 0) stop 'error while allocating'
    allocate(tmp_it1(npoints2))

    NPTS_STORED=npts2
    NPTS_INTERP=npoints2


    tmp_t1(:)=0.
    tmp_t2(:)=0.
    tmp_t3(:)=0.

    tmp_f1(:)=(0.,0.)
    tmp_f2(:)=(0.,0.)
    tmp_f3(:)=(0.,0.)

    nn=int(-t0/dt) ! what if this is not an integer number?

  if (myrank == 0) print *, '   starting from ',nn,' points before time 0'

  if (kpsv == 1) then

     ! for C_3=i sin(inc) (u=[sin(inc), cos(inc)])
     C_3=cmplx(0,1.)*p*al(nlayer)  ! amp. of incoming P in the bot. layer
     eta_p=sqrt(1./al(nlayer)**2-p**2) ! vertical slowness for lower layer
     if (myrank == 0) print *, 'Incoming P : C_3,  p, eta = ', C_3, p, eta_p

     N_mat(:,:) =(0.0,0.0)

     ! find out the wave coefficients in the bottom layer for all freqs -------------------------------
     do ii = 1, nf2
        om=2*pi*fvec(ii)
        ! propagation matrix
        call compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)),N_mat) !total-thickness=sum(h)
        a=N_mat(3,2); b=N_mat(3,4); c=N_mat(4,2); d=N_mat(4,4)
        delta_mat=a*d-b*c
        coeff(1,ii)=-(d*N_mat(3,3)-b*N_mat(4,3))/delta_mat*C_3
        coeff(2,ii)=-(-c*N_mat(3,3)+a*N_mat(4,3))/delta_mat*C_3
     enddo

     ! loop over all data points -------------------------------------------------
     do ip = 1, np  ! maybe this can be run faster by shifting t for diff. x of fixed z


        field_f=0.
        tdelay=p*(xx(ip)-x0)*cos(phi)+p*(yy(ip)-y0)*sin(phi)+eta_p*(0-z0)

        do ii = 1, nf2

           om=2*pi*fvec(ii)               !! pulsation
           stf_coeff=exp(-(om*tg/2)**2)   !! apodization window
           stf_coeff=stf_coeff*exp(cmplx(0,-1)*om*tdelay)

           !! zz(ip) is the height of point with respec to the lower layer
           call compute_N_rayleigh(al,be,mu,H,nlayer,om,p,zz(ip),N_mat)

           dx_f=N_mat(1,2)*coeff(1,ii)+N_mat(1,4)*coeff(2,ii)+N_mat(1,3)*C_3  ! y_1
           dz_f=N_mat(2,2)*coeff(1,ii)+N_mat(2,4)*coeff(2,ii)+N_mat(2,3)*C_3  ! y_3
           field_f(ii,1)=stf_coeff*dx_f*cmplx(0,-1)*cmplx(0,om)               ! (i om)u_x
           field_f(ii,2)=stf_coeff*dz_f*cmplx(0,om)                           ! (i om)u_z

           if (comp_stress) then
              txz_f=N_mat(3,2)*coeff(1,ii)+N_mat(3,4)*coeff(2,ii)+N_mat(3,3)*C_3 ! tilde{y}_4
              tzz_f=N_mat(4,2)*coeff(1,ii)+N_mat(4,4)*coeff(2,ii)+N_mat(4,3)*C_3 ! tilde{y}_6
              field_f(ii,3)=stf_coeff*om*p*(xi1(ip)*tzz_f-4*xim(ip)*dx_f) ! T_xx
              field_f(ii,4)=stf_coeff*om*p*txz_f*cmplx(0,-1)              ! T_xz
              field_f(ii,5)=stf_coeff*om*p*tzz_f                          ! T_zz
           endif

        enddo

        ! pad negative f, and convert to time series
        do ii = 2, nf2-1
           field_f(nf+2-ii,:) = conjg(field_f(ii,:))
        enddo

        !! inverse fast fourier transform
        field=0.
        do j = 1, nvar
           call FFTinv(lnpts2,field_f(:,j),zign_neg,dt,field(:,j))
           ! wrap around to start from t0: here one has to be careful if t0/dt is not
           ! exactly an integer, assume nn > 0
           if (nn > 0) then
              dtmp(1:nn)=field(npts2-nn+1:npts2,j)
              field(nn+1:npts2,j)=field(1:npts2-nn,j)
              field(1:nn,j)=dtmp(1:nn)
           else if (nn < 0) then
              dtmp(1:nn)=field(1:nn,j)
              field(1:npts-nn,j)=field(nn+1:npts,j)
              field(npts-nn+1:npts,j)=dtmp(1:nn)
           endif
        enddo


        !! store undersampled version of velocity  FK solution
        tmp_t1(:)=field(:,1)*cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vx_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:)=field(:,1)*sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vy_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:)=field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vz_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        !! compute traction
        do lpts = 1, NF_FOR_STORING

           sigma_rr = field(lpts,3)
           sigma_rt = 0.0
           sigma_rz = field(lpts,4)
           sigma_zz = field(lpts,5)
           sigma_tt = bdlambdamu(ip)*(sigma_rr+sigma_zz)
           sigma_tz = 0.0

           Txx_tmp = sigma_rr*cos(phi)*cos(phi)+sigma_tt*sin(phi)*sin(phi)
           Txy_tmp = cos(phi)*sin(phi)*(sigma_rr-sigma_tt)
           Txz_tmp = sigma_rz*cos(phi)
           Tyy_tmp = sigma_rr*sin(phi)*sin(phi)+sigma_tt*cos(phi)*cos(phi)
           Tyz_tmp = sigma_rz*sin(phi)
           Tzz_tmp = sigma_zz

           !! store directly the traction
           Tx_t(ip,lpts) = Txx_tmp*nmx(ip) +  Txy_tmp*nmy(ip) +  Txz_tmp*nmz(ip)
           Ty_t(ip,lpts) = Txy_tmp*nmx(ip) +  Tyy_tmp*nmy(ip) +  Tyz_tmp*nmz(ip)
           Tz_t(ip,lpts) = Txz_tmp*nmx(ip) +  Tyz_tmp*nmy(ip) +  Tzz_tmp*nmz(ip)

        enddo

        !! store undersamped version of tractions FK solution
        tmp_t1(1:NF_FOR_STORING)=Tx_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tx_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING)=Ty_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Ty_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING)=Tz_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tz_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

     enddo

  else if (kpsv == 2) then

     ! for C_2= sin(inc) (u=[cos(inc), sin(inc)])
     C_1= p*be(nlayer)  ! amp. of incoming S in the bot. layer
     eta_s=sqrt(1./be(nlayer)**2-p**2) ! vertical slowness for lower layer
     if (myrank == 0) print *, 'Incoming S :  C_1,  p, eta = ', C_1, p, eta_s

     N_mat(:,:) =(0.0,0.0)

     ! find out the wave coefficients in the bottom layer for all freqs
     do ii = 1, nf2
        om=2*pi*fvec(ii)
        ! propagation matrix
        !if (ii == nf2) pout = .true.
        call compute_N_rayleigh(al,be,mu,H,nlayer,om,p,sum(H(1:nlayer-1)),N_mat) !total-thickness=sum(h)
        a=N_mat(3,2); b=N_mat(3,4); c=N_mat(4,2); d=N_mat(4,4)
        delta_mat=a*d-b*c
        coeff(1,ii)=-(d*N_mat(3,1)-b*N_mat(4,1))/delta_mat*C_1
        coeff(2,ii)=-(-c*N_mat(3,1)+a*N_mat(4,1))/delta_mat*C_1
     enddo

     ! loop over all data points
     do ip = 1, np  ! maybe this can be run faster by shifting t for diff. x of fixed z

        field_f=0.
        tdelay=p*(xx(ip)-x0)*cos(phi)+p*(yy(ip)-y0)*sin(phi)+eta_s*(0-z0)

        do ii = 1, nf2

           om=2*pi*fvec(ii)
           stf_coeff=exp(-(om*tg/2)**2)*exp(cmplx(0,-1)*om*tdelay)

           ! z is the height of position with respect to the lowest layer interface.
           call compute_N_rayleigh(al,be,mu,H,nlayer,om,p,zz(ip),N_mat)

           dx_f=N_mat(1,2)*coeff(1,ii)+N_mat(1,4)*coeff(2,ii)+N_mat(1,1)*C_1  ! y_1
           dz_f=N_mat(2,2)*coeff(1,ii)+N_mat(2,4)*coeff(2,ii)+N_mat(2,1)*C_1  ! y_3
           field_f(ii,1)=stf_coeff*dx_f*cmplx(0,-1)*cmplx(0,om)  ! (i om)u_x(1.20)
           field_f(ii,2)=stf_coeff*dz_f*cmplx(0,om)              ! (i om)u_z

           if (comp_stress) then
              txz_f=N_mat(3,2)*coeff(1,ii)+N_mat(3,4)*coeff(2,ii)+N_mat(3,1)*C_1 ! tilde{y}_4
              tzz_f=N_mat(4,2)*coeff(1,ii)+N_mat(4,4)*coeff(2,ii)+N_mat(4,1)*C_1 ! tilde{y}_6
              field_f(ii,3)=stf_coeff*om*p*(xi1(ip)*tzz_f-4*xim(ip)*dx_f) !T_xx
              field_f(ii,4)=stf_coeff*om*p*txz_f*cmplx(0,-1) ! T_xz
              field_f(ii,5)=stf_coeff*om*p*tzz_f  ! T_zz
           endif

        enddo

        ! pad negative f, and convert to time series
        do ii = 2, nf2-1
           field_f(nf+2-ii,:) = conjg(field_f(ii,:))
        enddo

        field=0.
        do j = 1, nvar
           call FFTinv(lnpts2,field_f(:,j),zign_neg,dt,field(:,j))
           ! wrap around to start from t0: here one has to be careful if t0/dt is not
           ! exactly an integer, assume nn > 0
           if (nn > 0) then
              dtmp(1:nn)=field(npts2-nn+1:npts2,j)
              field(nn+1:npts2,j)=field(1:npts2-nn,j)
              field(1:nn,j)=dtmp(1:nn)
           else if (nn < 0) then
              dtmp(1:nn)=field(1:nn,j)
              field(1:npts-nn,j)=field(nn+1:npts,j)
              field(npts-nn+1:npts,j)=dtmp(1:nn)
           endif
        enddo

        !! store undersampled version of velocity  FK solution
        tmp_t1(:)=field(:,1)*cos(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vx_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:)=field(:,1)*sin(phi)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vy_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(:)=field(:,2)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        vz_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        !! compute traction
        do lpts = 1, NF_FOR_STORING

           sigma_rr = field(lpts,3)
           sigma_rt = 0.0
           sigma_rz = field(lpts,4)
           sigma_zz = field(lpts,5)
           sigma_tt = bdlambdamu(ip)*(sigma_rr+sigma_zz)
           sigma_tz = 0.0

           Txx_tmp = sigma_rr*cos(phi)*cos(phi)+sigma_tt*sin(phi)*sin(phi)
           Txy_tmp = cos(phi)*sin(phi)*(sigma_rr-sigma_tt)
           Txz_tmp = sigma_rz*cos(phi)
           Tyy_tmp = sigma_rr*sin(phi)*sin(phi)+sigma_tt*cos(phi)*cos(phi)
           Tyz_tmp = sigma_rz*sin(phi)
           Tzz_tmp = sigma_zz

           !! store directly the traction
           Tx_t(ip,lpts) = Txx_tmp*nmx(ip) +  Txy_tmp*nmy(ip) +  Txz_tmp*nmz(ip)
           Ty_t(ip,lpts) = Txy_tmp*nmx(ip) +  Tyy_tmp*nmy(ip) +  Tyz_tmp*nmz(ip)
           Tz_t(ip,lpts) = Txz_tmp*nmx(ip) +  Tyz_tmp*nmy(ip) +  Tzz_tmp*nmz(ip)

        enddo

        !! store undersamped version of tractions FK solution
        tmp_t1(1:NF_FOR_STORING)=Tx_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tx_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING)=Ty_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Ty_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

        tmp_t1(1:NF_FOR_STORING)=Tz_t(ip,1:NF_FOR_STORING)
        call compute_spline_coef_to_store(tmp_t1, npts2, tmp_t2)
        Tz_t(ip,1:NF_FOR_STORING)=tmp_t2(1:NF_FOR_STORING)

     enddo
  endif

  deallocate(fvec,coeff, field_f, field, dtmp)
  deallocate(tmp_f1, tmp_f2, tmp_f3, tmp_t1, tmp_t2, tmp_t3)

end subroutine FK

!==========================================================

subroutine compute_N_rayleigh(alpha, beta, mu, H, nlayer, om, p, ht, Nmat)

  ! assumes that ht = 0 is the bottom interface

  use constants
  implicit none

  ! precision for complex
  integer, parameter :: CUSTOM_CMPLX = 8

  ! input
  integer,                                       intent(in)   :: nlayer
  real(kind=CUSTOM_REAL),     dimension(nlayer), intent(in)   :: alpha, beta, mu, H
  real(kind=CUSTOM_REAL),                        intent(in)   :: om,p,ht

  ! output
  complex(kind=CUSTOM_CMPLX), dimension(4,4),    intent(inout) :: Nmat(4,4)

  ! local vars
  integer                                         :: i, j, ilayer
  complex(kind=CUSTOM_CMPLX), dimension(nlayer)   :: eta_al, eta_be, nu_al, nu_be
  complex(kind=CUSTOM_CMPLX), dimension(4,4)      :: Emat, Gmat, Player
  complex(kind=CUSTOM_CMPLX)                      :: ca, sa, xa, ya, cb, sb, xb, yb, g1, mul, c1, c2
  real(kind=CUSTOM_REAL),     dimension(nlayer)   :: gamma0, gamma1
  real(kind=CUSTOM_REAL),     dimension(nlayer)   :: hh

  if (nlayer < 1) stop 'nlayer has to be larger than or equal to 1'

  do i=1,nlayer
     eta_al(i)=-cmplx(0,1)*sqrt(1.0/alpha(i)+p)*sqrt(1.0/alpha(i)-p) ! i*vertical slowness, purely imaginary
     eta_be(i)=-cmplx(0,1)*sqrt(1.0/beta(i)+p)*sqrt(1.0/beta(i)-p)
     nu_al(i)=om*eta_al(i)
     nu_be(i)=om*eta_be(i) ! i * vertical wavenumber
     gamma0(i)=2*beta(i)**2*p**2
     gamma1(i)=1-1/gamma0(i)
  enddo

  ! note Emat is not omega dependent
  Emat(1,1) =  eta_be(nlayer)/p
  Emat(1,2) = -Emat(1,1)
  Emat(1,3) = 1
  Emat(1,4) = 1
  Emat(2,1) = 1
  Emat(2,2) = 1
  Emat(2,3) =  eta_al(nlayer)/p
  Emat(2,4) = -Emat(2,3)

  Emat(3,1) = 2*mu(nlayer)*gamma1(nlayer)
  Emat(3,2) = Emat(3,1)
  Emat(3,3) = 2*mu(nlayer)*eta_al(nlayer)/p
  Emat(3,4) = -Emat(3,3)
  Emat(4,1) = 2*mu(nlayer)*eta_be(nlayer)/p
  Emat(4,2) = -Emat(4,1)
  Emat(4,3) = Emat(3,1)
  Emat(4,4) = Emat(3,1)

  if (ht > sum(h(1:nlayer-1))) stop 'Z point is located in the air above the surface rather than in the solid!'

  ! figure out the location z with respect to layer stack
  if (ht <= 0) then ! in lower half space
     Gmat = 0.
     Gmat(1,1)=exp(nu_be(nlayer)*ht)
     Gmat(2,2)=exp(-nu_be(nlayer)*ht)
     Gmat(3,3)=exp(nu_al(nlayer)*ht)
     Gmat(4,4)=exp(-nu_al(nlayer)*ht)
     Nmat = matmul(Emat,Gmat)
  else ! in layers
     hh=H
     ilayer=nlayer
     do j = nlayer-1 , 1 , -1
        if (ht <= sum(H(j:nlayer-1))) then
           ilayer=j; exit
        endif
     enddo
     hh(ilayer+1:nlayer-1)=H(ilayer+1:nlayer-1)
     hh(ilayer)=ht-sum(H(ilayer+1:nlayer-1))
     if (hh(ilayer) < 0) stop 'Error setting layer thickness'

     ! compute propagation matrices
     Nmat = Emat

     do j=nlayer-1, ilayer, -1
        c1=nu_al(j)*hh(j)
        ca=(exp(c1)+exp(-c1))/2; sa=(exp(c1)-exp(-c1))/2
        xa=eta_al(j)*sa/p; ya=p*sa/eta_al(j)
        c2=nu_be(j)*hh(j)
        cb=(exp(c2)+exp(-c2))/2; sb=(exp(c2)-exp(-c2))/2
        xb=eta_be(j)*sb/p; yb=p*sb/eta_be(j)
        g1=gamma1(j); mul=mu(j)

        Player(1,1) = ca-g1*cb
        Player(1,2) = xb-g1*ya
        Player(1,3) = (ya-xb)/(2*mul)
        Player(1,4) = (cb-ca)/(2*mul)
        Player(2,1) = xa-g1*yb
        Player(2,2) = cb-g1*ca
        Player(2,3) = (ca-cb)/(2*mul)
        Player(2,4) = (yb-xa)/(2*mul)
        Player(3,1) = 2*mul*(xa-g1**2*yb)
        Player(3,2) = 2*mul*g1*(cb-ca)
        Player(3,3) = ca-g1*cb
        Player(3,4) = g1*yb-xa
        Player(4,1) = 2*mul*g1*(ca-cb)
        Player(4,2) = 2*mul*(xb-g1**2*ya)
        Player(4,3) = g1*ya-xb
        Player(4,4) = cb-g1*ca

        ! if (pout) print *,'j,player',j,player

        Nmat=gamma0(j)*matmul(Player,Nmat)
     enddo

  endif

end subroutine compute_N_rayleigh

!==========================================================

  subroutine FFT(n,xi,zign,dtt)

! Fourier transform
! This inputs AND outputs a complex function.
! The convention is FFT --> e^(-iwt)
! numerical factor for Plancherel theorem: planch_fac = dble(NPT * dt * dt)

  implicit none

      integer, parameter :: CUSTOM_REAL = 4, CUSTOM_CMPLX = 8

      integer :: n
      integer :: lblock,k,FK,jh,ii,istart
      integer :: l,iblock,nblock,i,lbhalf,j,lx

      complex(kind=CUSTOM_CMPLX),dimension(*) :: xi
      complex(kind=CUSTOM_CMPLX) :: wk, hold, q

      real(kind=CUSTOM_REAL) :: zign,flx,inv_of_flx,v,dtt,pi

!! DK DK here is the hardwired maximum size of the array
!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
      real(kind=CUSTOM_REAL) :: m(30)

      pi = acos(-1.0)

!! DK DK added this sanity check
      if (n > 30) stop 'error: the FK FTT routine has an hardwired maximum of 30 levels'
!! DK DK in any case the line below imposes a maximum of 31, otherwise the integer 2**n will overflow

      lx = 2**n

!! DK DK Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
!! DK DK Aug 2016: this should be moved to the calling program and precomputed once and for all
      do i = 1,n
        m(i) = 2**(n-i)
      enddo

      do l = 1,n

      nblock = 2**(l-1)
      lblock = lx/nblock
      lbhalf = lblock/2

      k = 0

      do iblock = 1,nblock

      FK = k
      flx = lx

      v = zign*2.*PI*FK/flx         ! Fourier convention
     ! - sign: MATLAB convention: forward e^{-i om t}
     ! + sign: engineering convention: forward e^{i om t}
      wk = cmplx(cos(v),-sin(v))   ! sign change to -sin(v) or sin(v)
      istart = lblock*(iblock-1)

      do i = 1,lbhalf
        j  = istart+i
        jh = j+lbhalf
        q = xi(jh)*wk
        xi(jh) = xi(j)-q
        xi(j)  = xi(j)+q
      enddo

      do i = 2,n
        ii = i
        if (k < m(i)) goto 4
        k = k-m(i)
      enddo

    4 k = k+m(ii)

      enddo
      enddo

      k = 0

      do j = 1,lx
        if (k < j) goto 5

        hold = xi(j)
        xi(j) = xi(k+1)
        xi(k+1) = hold

    5   do i = 1,n
          ii = i
          if (k < m(i)) goto 7
          k = k-m(i)
        enddo

    7   k = k+m(ii)
      enddo

      ! final steps deal with dt factors
      if (zign > 0.) then      ! FORWARD FFT

        xi(1:lx) = xi(1:lx) * dtt    ! multiplication by dt

      else                     ! REVERSE FFT

        flx = flx*dtt
        inv_of_flx = 1._CUSTOM_REAL / flx

!! DK DK Aug 2016: changed to multiplication by the precomputed inverse to make the routine faster
!       xi(1:lx) = xi(1:lx) / flx         ! division by dt
        xi(1:lx) = xi(1:lx) * inv_of_flx  ! division by dt

      endif

  end subroutine FFT

!==========================================================

  subroutine FFTinv(npow,s,zign,dtt,r)

! inverse Fourier transform -- calls FFT

  implicit none

      integer, parameter :: CUSTOM_CMPLX = 8,CUSTOM_REAL = 4

      integer :: npow, nsmp, nhalf
      real(kind=CUSTOM_REAL)  :: dtt,zign

      complex(kind=CUSTOM_CMPLX), intent(in) :: s(*)
      real(kind=CUSTOM_REAL), intent(out) :: r(*)   ! note that this is real, not double precision

      nsmp = 2**npow
      nhalf = nsmp/2

      call rspec(s,nhalf)   ! restructuring
      call FFT(npow,s,zign,dtt)    ! Fourier transform

      r(1:nsmp) = real(s(1:nsmp))     ! take the real part

  end subroutine FFTinv

!==========================================================

  subroutine rspec(s,np2)

  implicit none

      integer, parameter :: CUSTOM_CMPLX = 8

      integer :: np2,n,n1,i

      complex(kind=CUSTOM_CMPLX) :: s(*)

      n = 2*np2
      n1 = np2+1

      s(n1) = 0.0
      s(1)  = cmplx(real(s(1)),0.0)

      do i = 1,np2
        s(np2+i) = conjg(s(np2+2-i))
      enddo

  end subroutine rspec


!-------------------------------------------------------------------

   subroutine find_size_of_working_arrays(deltat, tmax_fk, NF_FOR_STORING, &
    NF_FOR_FFT, NPOW_FOR_INTERP, np_resampling, DF_FK)

      use constants

      real(kind=CUSTOM_REAL),intent(inout) :: tmax_fk
      real(kind=CUSTOM_REAL),intent(inout) :: DF_FK, deltat
      integer,               intent(inout) :: NF_FOR_STORING, NF_FOR_FFT, NPOW_FOR_INTERP, np_resampling

      real(kind=CUSTOM_REAL)               :: df, dt_min_fk, Frq_ech_Fk

      !! sampling frequency to store fk solution
      Frq_ech_Fk = 10._CUSTOM_REAL  !! WM WM TO DO PUT THIS IN PARAMETER

      !! sampling time step to store fk solution
      dt_min_fk = 1. /  Frq_ech_Fk

      !!  compute resampling rate
      np_resampling  = floor(dt_min_fk / deltat)

      !! update dt for fk with respect to integer np_resampling
      dt_min_fk = np_resampling * deltat  !! this is the time step sampling for FK storage

      !! compute number of time steps to store
      NF_FOR_STORING  = ceiling( tmax_fk / dt_min_fk)

      !! in power of two
      NF_FOR_STORING  =   ceiling(log(real(NF_FOR_STORING))/log(2.))

      !! multiply by 2 in order to do an inverse FFT
      NF_FOR_FFT      =   2** (NF_FOR_STORING+1)

      NPOW_FOR_INTERP =   NF_FOR_STORING+1
      NF_FOR_STORING  =   2** NF_FOR_STORING

      !! now we have this new time window
      tmax_fk = dt_min_fk * (NF_FOR_FFT - 1)

      !! step in frequency for fk
      df = 1. / tmax_fk

      DF_FK = df

    end subroutine find_size_of_working_arrays


!! #################  INTERPOLATION ROUTINES IN TIME DOMAIN ######################################

    !! compute and store spline coefficients
subroutine compute_spline_coef_to_store(Sig, npts, spline_coeff)


  use constants
  implicit none
  integer,                                             intent(in)     :: npts
  real(kind=CUSTOM_REAL), dimension(npts),             intent(in)     :: Sig
  real(kind=CUSTOM_REAL), dimension(npts),             intent(inout)  :: spline_coeff


  !! computation in double precision
  double precision                                                    :: error=1.d-24
  double precision                                                    :: z1, zn, sumc
  double precision, dimension(:), allocatable                         :: c
  integer                                                             :: i, n_init

  allocate(c(npts))

  ! Compute pole value
  z1 = dsqrt(3.d0)-2.d0
  c(:)=dble(Sig(:)) * (1.d0-z1) *( 1.d0 - 1.d0/ z1)

  ! Initialisation causal filter
  n_init = ceiling(log(error)/log(abs(z1)))
  sumc = c(1)
  zn = z1
  do i = 1,n_init
     sumc = sumc + zn*c(i)
     zn = zn*z1
  enddo
  c(1) = sumc

  ! Causal filter
  do i = 2,npts
     c(i) = c(i) + z1* c(i-1)
  enddo

  ! Initialisation anti-causal filter
  c(npts) = ( z1 / (z1-1.d0) ) *c(npts)
  do i = npts-1,1,-1
     c(i) = z1*(c(i+1)-c(i))
  enddo

  !! store spline coeff in CUSTOM_REAL precision
  spline_coeff(:)=c(:)

  deallocate(c)

end subroutine compute_spline_coef_to_store

