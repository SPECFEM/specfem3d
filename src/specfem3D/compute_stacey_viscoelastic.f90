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
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        veloc,rho_vp,rho_vs, &
                        ispec_is_elastic,SIMULATION_TYPE,SAVE_FORWARD, &
                        it, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants

  use specfem_par, only: it_dsm, Veloc_dsm_boundary, Tract_dsm_boundary, Veloc_axisem, Tract_axisem

  use shared_parameters, only: COUPLE_WITH_EXTERNAL_CODE, EXTERNAL_CODE_TYPE

  implicit none

  integer :: NSPEC_AB,NGLOB_AB

! acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

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

  !! CD modif. : begin (implemented by VM) !! For coupling with DSM

  integer :: kaxisem

  if (COUPLE_WITH_EXTERNAL_CODE) then

    if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_DSM) then

      if (old_DSM_coupling_from_Vadim) then
        if (phase_is_inner .eqv. .false.) then
          if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
            call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
          endif
        endif
      else
        !! MODIFS DE NOBU 2D
      endif

    else if (EXTERNAL_CODE_TYPE == EXTERNAL_CODE_IS_AXISEM) then

      if (phase_is_inner .eqv. .false.) then
        call read_axisem_file(Veloc_axisem,Tract_axisem,num_abs_boundary_faces*NGLLSQUARE)
      endif

    endif

  endif

  !! CD modif. : end

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      if (ispec_is_elastic(ispec)) then

        ! reference gll points on boundary face
        do igll = 1,NGLLSQUARE

          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets velocity
          iglob=ibool(i,j,k,ispec)
          vx=veloc(1,iglob)
          vy=veloc(2,iglob)
          vz=veloc(3,iglob)

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
          !! CD CD

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

          !! CD CD

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

        enddo
      endif ! ispec_is_elastic
    endif ! ispec_is_inner
  enddo

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value only when second phase is running
    if (phase_is_inner .eqv. .true.) call write_abs(IOABS,b_absorb_field,b_reclen_field,it)

  endif

  !! CD CD !! For coupling with DSM
  if (COUPLE_WITH_EXTERNAL_CODE) then !! To verify for NOBU version

    if (phase_is_inner .eqv. .true.) it_dsm = it_dsm + 1

  !! TODO: maybe call integrand_for_computing_Kirchoff_Helmholtz_integral here

  endif
  !! CD CD

  end subroutine compute_stacey_viscoelastic
!
!=====================================================================

! for elastic solver

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_bpwf(NSPEC_AB, &
                        ibool,ispec_is_inner,phase_is_inner, &
                        abs_boundary_ijk,abs_boundary_ispec, &
                        num_abs_boundary_faces, &
                        ispec_is_elastic,SIMULATION_TYPE, &
                        NSTEP,it,NGLOB_ADJOINT,b_accel, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field)

  use constants
  use specfem_par,only: myrank

  implicit none

  integer :: NSPEC_AB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! communication overlap
  logical, dimension(NSPEC_AB) :: ispec_is_inner
  logical :: phase_is_inner

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
    call exit_MPI(myrank,'error calling routine compute_stacey_viscoelastic_bpwf() with wrong SIMULATION_TYPE')

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  ! reads in absorbing boundary array when first phase is running
  if (phase_is_inner .eqv. .false.) then
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
  endif

  ! absorbs absorbing-boundary surface using Stacey condition (Clayton & Enquist)
  do iface=1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_inner(ispec) .eqv. phase_is_inner) then

      if (ispec_is_elastic(ispec)) then
        ! reference gll points on boundary face
        do igll = 1,NGLLSQUARE
          ! gets local indices for GLL point
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)

          ! gets velocity
          iglob=ibool(i,j,k,ispec)

          ! adjoint simulations
          b_accel(:,iglob) = b_accel(:,iglob) - b_absorb_field(:,igll,iface)
        enddo
      endif ! ispec_is_elastic
    endif ! ispec_is_inner
  enddo

  end subroutine compute_stacey_viscoelastic_bpwf

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
  do iface=1,num_abs_boundary_faces

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
!=============================================================================

!
!=====================================================================
! for elastic solver on GPU

! absorbing boundary term for elastic media (Stacey conditions)

  subroutine compute_stacey_viscoelastic_GPU(phase_is_inner,num_abs_boundary_faces, &
                        SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                        b_num_abs_boundary_faces,b_reclen_field,b_absorb_field, &
                        Mesh_pointer,it_dsm,Veloc_dsm_boundary,Tract_dsm_boundary,COUPLE_WITH_EXTERNAL_CODE)

  use constants

  implicit none

! communication overlap
  logical :: phase_is_inner

! absorbing boundary surface
  integer :: num_abs_boundary_faces

! adjoint simulations
  integer:: SIMULATION_TYPE
  integer:: NSTEP,it
  integer:: b_num_abs_boundary_faces,b_reclen_field
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLSQUARE,b_num_abs_boundary_faces):: b_absorb_field

  logical:: SAVE_FORWARD
  logical:: COUPLE_WITH_EXTERNAL_CODE

  ! GPU_MODE variables
  integer(kind=8) :: Mesh_pointer

  !! CD CD : begin !! For coupling with DSM
  real(kind=CUSTOM_REAL) :: Veloc_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)
  real(kind=CUSTOM_REAL) :: Tract_dsm_boundary(3,Ntime_step_dsm,NGLLSQUARE,num_abs_boundary_faces)

  integer :: it_dsm

  if (COUPLE_WITH_EXTERNAL_CODE) then
    if (old_DSM_coupling_from_Vadim) then
      if (phase_is_inner .eqv. .false.) then
        if (mod(it_dsm,Ntime_step_dsm+1) == 0 .or. it == 1) then
          call read_dsm_file(Veloc_dsm_boundary,Tract_dsm_boundary,num_abs_boundary_faces,it_dsm)
        endif
      endif
    else
      !! MODIFS DE NOBU 2D
    endif
  endif

  !! CD CD : end

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array when first phase is running
    if (phase_is_inner .eqv. .false.) then
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      call read_abs(IOABS,b_absorb_field,b_reclen_field,NSTEP-it+1)
    endif
  endif !adjoint

  call compute_stacey_viscoelastic_cuda(Mesh_pointer,phase_is_inner,b_absorb_field)

  ! adjoint simulations: stores absorbed wavefield part
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! writes out absorbing boundary value only when second phase is running
    if (phase_is_inner .eqv. .true.) then
      call write_abs(IOABS,b_absorb_field,b_reclen_field,it)
    endif
  endif

  !! CD CD : begin !! For coupling with DSM
  if (COUPLE_WITH_EXTERNAL_CODE) then !! To verify for NOBU version
    if (phase_is_inner .eqv. .true.) then
      it_dsm = it_dsm + 1
    endif
  endif
  !! CD CD : end

  end subroutine compute_stacey_viscoelastic_GPU

