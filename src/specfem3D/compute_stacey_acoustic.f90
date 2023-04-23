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

! for acoustic solver

  subroutine compute_stacey_acoustic_forward(NSPEC_AB,NGLOB_AB, &
                                             potential_dot_dot_acoustic,potential_dot_acoustic, &
                                             ibool,iphase, &
                                             abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                             num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                                             it,b_reclen_potential, &
                                             b_absorb_potential,b_num_abs_boundary_faces)

  use constants
  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  ! wavefield injection
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_potential

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

  ! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

  ! communication overlap
  integer,intent(in) :: iphase

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rhostore,kappastore
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

  ! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

  ! adjoint simulations
  integer,intent(in) :: it
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces),intent(inout) :: b_absorb_potential

  ! local parameters
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw,absorbl
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer:: reclen1,reclen2

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! injecting boundary wavefield
  if (COUPLE_WITH_INJECTION_TECHNIQUE .and. SIMULATION_TYPE == 1) then
    ! adds boundary contribution from injected wavefield
    call compute_coupled_injection_contribution_ac(NGLOB_AB,potential_dot_dot_acoustic,iphase,it)
  endif

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets global index
        iglob = ibool(i,j,k,ispec)

        ! determines bulk sound speed
        rhol = rhostore(i,j,k,ispec)
        cpl = sqrt( kappastore(i,j,k,ispec) / rhol )

        ! gets associated, weighted jacobian
        jacobianw = abs_boundary_jacobian2Dw(igll,iface)

        ! Sommerfeld condition
        absorbl = potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - absorbl

        ! adjoint simulations
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_absorb_potential(igll,iface) = absorbl
        endif !adjoint

       enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  ! for kernel simulations: stores absorbed wavefield part
  if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
    ! adds boundary injection contribution to stacey contribution before saving to disk
    ! this avoids storing the boundary injection arrays as a separate file
    ! for kernel simulations to reconstruct forward wavefields.
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      b_absorb_potential(:,:) = b_absorb_potential(:,:) + b_boundary_injection_potential(:,:)
    endif

    ! writes out absorbing boundary value
    call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
  endif

  end subroutine compute_stacey_acoustic_forward
!
!=====================================================================
! for acoustic solver for back propagation wave field

  subroutine compute_stacey_acoustic_backward(NSPEC_AB, &
                                              ibool,iphase, &
                                              abs_boundary_ijk,abs_boundary_ispec, &
                                              num_abs_boundary_faces,ispec_is_acoustic, &
                                              SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                                              b_potential_dot_dot_acoustic,b_reclen_potential, &
                                              b_absorb_potential,b_num_abs_boundary_faces)

  use constants

  implicit none

  integer,intent(in) :: NSPEC_AB

! potentials
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! communication overlap
  integer,intent(in) :: iphase

  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

! adjoint simulations
  integer,intent(in) :: SIMULATION_TYPE
  integer,intent(in) :: NSTEP,it,NGLOB_ADJOINT
  integer,intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLOB_ADJOINT),intent(inout) :: b_potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces),intent(in) :: b_absorb_potential

! local parameters
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer:: reclen1,reclen2

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! adjoint simulations:
  if (SIMULATION_TYPE == 3) then
    ! reads in absorbing boundary array
    ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
    call read_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,NSTEP-it+1)
  endif !adjoint

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! gets global index
        iglob=ibool(i,j,k,ispec)

        ! adjoint simulations
        if (SIMULATION_TYPE == 3) then
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                              - b_absorb_potential(igll,iface)
        endif !adjoint

      enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  end subroutine compute_stacey_acoustic_backward

!
!=====================================================================
!

  subroutine compute_stacey_acoustic_backward_undoatt(NSPEC_AB,NGLOB_AB, &
                                                      b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                      ibool,iphase, &
                                                      abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                                                      num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic)

  use constants

  implicit none

  integer,intent(in) :: NSPEC_AB,NGLOB_AB

! potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(in) :: b_potential_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: b_potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: ibool

! communication overlap
  integer,intent(in) :: iphase

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB),intent(in) :: rhostore,kappastore
  logical, dimension(NSPEC_AB),intent(in) :: ispec_is_acoustic

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces
  real(kind=CUSTOM_REAL),intent(in) :: abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces)
  integer,intent(in) :: abs_boundary_ispec(num_abs_boundary_faces)

! local parameters
  real(kind=CUSTOM_REAL) :: rhol,cpl,jacobianw,absorbl
  integer :: ispec,iglob,i,j,k,iface,igll
  !integer:: reclen1,reclen2

  ! only add these contributions in first pass
  if (iphase /= 1) return

  ! checks if anything to do
  if (num_abs_boundary_faces == 0) return

  ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
  do iface = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      ! reference GLL points on boundary face
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
        absorbl = b_potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
        b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - absorbl

       enddo
    endif ! ispec_is_acoustic
  enddo ! num_abs_boundary_faces

  end subroutine compute_stacey_acoustic_backward_undoatt

!
!=====================================================================
! for acoustic solver on GPU

  subroutine compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                                         NSTEP,it, &
                                         b_reclen_potential,b_absorb_potential, &
                                         b_num_abs_boundary_faces,Mesh_pointer, &
                                         FORWARD_OR_ADJOINT)

  use constants

  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,UNDO_ATTENUATION_AND_OR_PML

  use specfem_par, only: SIMULATION_TYPE,SAVE_FORWARD

  ! wavefield injection
  use specfem_par, only: NGLOB_AB
  use specfem_par_acoustic, only: potential_dot_dot_acoustic
  ! boundary injection wavefield parts for saving together with b_absorb_field
  use specfem_par_coupling, only: b_boundary_injection_potential

  implicit none

! potentials

! communication overlap
  integer,intent(in) :: iphase

! absorbing boundary surface
  integer,intent(in) :: num_abs_boundary_faces

! adjoint simulations
  integer, intent(in) :: NSTEP,it
  integer, intent(in) :: b_num_abs_boundary_faces,b_reclen_potential
  real(kind=CUSTOM_REAL),dimension(NGLLSQUARE,b_num_abs_boundary_faces), intent(inout) :: b_absorb_potential

  ! GPU_MODE variables
  integer(kind=8), intent(in) :: Mesh_pointer
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
    ! transfers potential to the CPU
    call transfer_dot_dot_from_device(NGLOB_AB, potential_dot_dot_acoustic, Mesh_pointer)

    ! adds boundary contribution from injected wavefield
    call compute_coupled_injection_contribution_ac(NGLOB_AB,potential_dot_dot_acoustic,iphase,it)

    ! transfers updated potential field back to the GPU
    call transfer_dot_dot_to_device(NGLOB_AB, potential_dot_dot_acoustic, Mesh_pointer)
  endif

  if (UNDO_ATTENUATION_AND_OR_PML) then
    ! no need to store boundaries on disk
    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_acoustic_undoatt_cuda(Mesh_pointer,iphase,FORWARD_OR_ADJOINT)
  else
    ! adjoint simulations:
    if (SIMULATION_TYPE == 3) then
      ! reads in absorbing boundary array (when first phase is running)
      ! note: the index NSTEP-it+1 is valid if b_displ is read in after the Newmark scheme
      call read_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,NSTEP-it+1)
    endif !adjoint

    ! absorbs absorbing-boundary surface using Sommerfeld condition (vanishing field in the outer-space)
    call compute_stacey_acoustic_cuda(Mesh_pointer,iphase,b_absorb_potential,FORWARD_OR_ADJOINT)

    ! for kernel simulations: stores absorbed wavefield part
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! adds boundary injection contribution to stacey contribution before saving to disk
      ! this avoids storing the boundary injection arrays as a separate file
      ! for kernel simulations to reconstruct forward wavefields.
      if (COUPLE_WITH_INJECTION_TECHNIQUE) then
        b_absorb_potential(:,:) = b_absorb_potential(:,:) + b_boundary_injection_potential(:,:)
      endif

      ! writes out absorbing boundary value
      call write_abs(IOABS_AC,b_absorb_potential,b_reclen_potential,it)
    endif
  endif

  end subroutine compute_stacey_acoustic_GPU

!=============================================================================
!
! For coupling with external code
!
!=============================================================================

  subroutine compute_coupled_injection_contribution_ac(NGLOB_AB,potential_dot_dot_acoustic,iphase,it)

  use constants

  use specfem_par, only: SAVE_STACEY,SIMULATION_TYPE

  use specfem_par, only: abs_boundary_normal,abs_boundary_jacobian2Dw, &
    abs_boundary_ijk,abs_boundary_ispec, &
    num_abs_boundary_faces

  use specfem_par, only: ibool, rhostore, kappastore
  use specfem_par_acoustic, only: ispec_is_acoustic

  ! boundary coupling
  use shared_parameters, only: COUPLE_WITH_INJECTION_TECHNIQUE,INJECTION_TECHNIQUE_TYPE
  ! FK3D calculation
  use specfem_par_coupling, only: type_kpsv_fk, ipt_table, NP_RESAMP, Veloc_FK, Tract_FK
  ! boundary injection wavefield parts for saving together with b_absorb_potential
  use specfem_par_coupling, only: b_boundary_injection_potential

  implicit none

  integer,intent(in) :: NGLOB_AB

  ! acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB),intent(inout) :: potential_dot_dot_acoustic

  ! communication overlap
  integer,intent(in) :: iphase

  ! adjoint simulations
  integer,intent(in) :: it

  ! local parameters
  real(kind=CUSTOM_REAL) :: nx,ny,nz,jacobianw
  integer :: ispec,iglob,i,j,k,iface,igll
  ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
  ! FK surface
  integer :: ipt, ii, kk, iim1, iip1, iip2
  real(kind=CUSTOM_REAL) :: cs1,cs2,cs3,cs4,w
  real(kind=CUSTOM_REAL) :: vx_FK,vy_FK,vz_FK,tx_FK,ty_FK,tz_FK
  real(kind=CUSTOM_REAL) :: vx,vy,vz,tx,ty,tz,vn,tn
  real(kind=CUSTOM_REAL) :: rhol,kappal,cpl,absorbl

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
    ! safety stop
    stop 'DSM coupling for acoustic domains is not implemented yet'

  case (INJECTION_TECHNIQUE_IS_AXISEM)
    ! AxiSEM coupling
    ! safety stop
    stop 'AxiSEM coupling for acoustic domains is not implemented yet'

  case (INJECTION_TECHNIQUE_IS_FK)
    ! FK coupling
    ! safety check only P-wave incident
    if (type_kpsv_fk /= 1) return

    !! find indices
    ! example: np_resamp = 2 and it = 1,2,3,4,5,6, ..
    !          --> ii = 1,1,2,2,3,3,..
    ii = floor( real(it + NP_RESAMP - 1) / real( NP_RESAMP))
    ! example: --> kk = 1,2,1,2,1,2,,..
    kk = it - (ii-1) * NP_RESAMP
    w = dble(kk-1) / dble(NP_RESAMP)

    ! Cubic spline values
    cs4 = w*w*w/6.d0
    cs1 = 1.d0/6.d0+w*(w-1.d0)/2.d0-cs4
    cs3 = w+cs1-2.d0*cs4
    cs2 = 1.d0-cs1-cs3-cs4

    ! interpolation indices
    iim1 = ii-1
    iip1 = ii+1
    iip2 = ii+2
  end select

  ! wavefield injection
  ! note: here, we inject the wavefield by coupling between acoustic/elastic domains.
  !       at the moment, only the FK solution works and this analytical wavefield solution is available only for elastic domains.
  !       thus, instead of velocity/stress as in the elastic boundary injection, we use the displacement solution to inject
  !       as separate terms to the right-hand-side. The Stacey absorbing contribution will be computed as usual.
  !       this might need a better formulation in future...

  do iface = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(iface)
    if (ispec_is_acoustic(ispec)) then
      ! GLL points on boundary face
      do igll = 1,NGLLSQUARE
        ! gets injected velocity & stress
        ! velocity and stresses would be subtracted from total, therefore we add a minus sign when getting the values
        select case(INJECTION_TECHNIQUE_TYPE)
        case (INJECTION_TECHNIQUE_IS_DSM)                 !! To verify for NOBU version
          ! not implemented yet
          return
        case (INJECTION_TECHNIQUE_IS_AXISEM)              !! VM VM add AxiSEM
          ! not implemented yet
          return
        case (INJECTION_TECHNIQUE_IS_FK)
          ! added by Ping Tong (TP / Tong Ping) for the FK3D calculation
          ! point index using table lookup
          ipt = ipt_table(igll,iface)

          ! interpolates velocity/stress
          vx_FK = cs1 * Veloc_FK(1,ipt,iim1) + cs2 * Veloc_FK(1,ipt,ii) + cs3 * Veloc_FK(1,ipt,iip1) + cs4 * Veloc_FK(1,ipt,iip2)
          vy_FK = cs1 * Veloc_FK(2,ipt,iim1) + cs2 * Veloc_FK(2,ipt,ii) + cs3 * Veloc_FK(2,ipt,iip1) + cs4 * Veloc_FK(2,ipt,iip2)
          vz_FK = cs1 * Veloc_FK(3,ipt,iim1) + cs2 * Veloc_FK(3,ipt,ii) + cs3 * Veloc_FK(3,ipt,iip1) + cs4 * Veloc_FK(3,ipt,iip2)

          tx_FK = cs1 * Tract_FK(1,ipt,iim1) + cs2 * Tract_FK(1,ipt,ii) + cs3 * Tract_FK(1,ipt,iip1) + cs4 * Tract_FK(1,ipt,iip2)
          ty_FK = cs1 * Tract_FK(2,ipt,iim1) + cs2 * Tract_FK(2,ipt,ii) + cs3 * Tract_FK(2,ipt,iip1) + cs4 * Tract_FK(2,ipt,iip2)
          tz_FK = cs1 * Tract_FK(3,ipt,iim1) + cs2 * Tract_FK(3,ipt,ii) + cs3 * Tract_FK(3,ipt,iip1) + cs4 * Tract_FK(3,ipt,iip2)

          ! velocity
          vx = - vx_FK
          vy = - vy_FK
          vz = - vz_FK
          ! stress
          tx = - tx_FK
          ty = - ty_FK
          tz = - tz_FK
        end select

        ! gets local indices for GLL point
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        ! global index
        iglob = ibool(i,j,k,ispec)

        ! determines bulk sound speed
        rhol = rhostore(i,j,k,ispec)
        kappal = kappastore(i,j,k,ispec)
        cpl = sqrt( kappal / rhol )

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
        ! velocity vector component * vp * rho in normal direction (+ vs * rho component tangential to it for elastic case)
        tx = tx + rhol * cpl * vn * nx
        ty = ty + rhol * cpl * vn * ny
        tz = tz + rhol * cpl * vn * nz

        ! stress component in normal direction (normal points out of element)
        tn = tx*nx + ty*ny + tz*nz

        ! total contribution (traction and Stacey/Sommerfeld absorbing)
        !
        ! note: the additional factor 1/(cpl**2 * rhol**2) is ad-hoc and might need a proper derivation.
        !       here, it is added to first undo the factor (rhol*cpl) used to add the velocity contribution to the traction together
        !       with another (rhol*cpl) used in the expressions of the Sommerfeld condition, where
        !          absorbl = potential_dot_acoustic(iglob) * jacobianw / cpl / rhol
        absorbl = tn * jacobianw / (cpl**2 * rhol**2)
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - absorbl

        ! for kernel simulations: stores contribution to buffer array and add it to stacey buffer before saving to disk
        if (SAVE_STACEY .and. SIMULATION_TYPE == 1) then
          b_boundary_injection_potential(igll,iface) = - absorbl
        endif
      enddo
    endif
  enddo

  end subroutine compute_coupled_injection_contribution_ac

