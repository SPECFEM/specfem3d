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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine pml_allocate_arrays()

  use pml_par
  use specfem_par, only: NSPEC_AB,PML_CONDITIONS,SIMULATION_TYPE,SAVE_FORWARD,NSTEP,myrank,prname !
  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION,num_coupling_ac_el_faces
  use specfem_par_elastic, only: ELASTIC_SIMULATION

  implicit none

  ! local parameters
  integer :: ier,b_nglob_interface_PML_elastic,b_nglob_interface_PML_acoustic
  integer(kind=8) :: filesize

  ! checks PML flag
  if (.not. PML_CONDITIONS) return

  ! slices without PML elements
  if (NSPEC_CPML == 0) then
    ! dummy allocation with a size of 1 for all the PML arrays that have not yet been allocated
    ! in order to be able to use these arrays as arguments in subroutine calls
    call pml_allocate_arrays_dummy()
    return
  endif

  ! C-PML spectral elements local indexing
  allocate(spec_to_CPML(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'error allocating array spec_to_CPML'

  ! C-PML element type array: 1 = face, 2 = edge, 3 = corner
  allocate(CPML_type(NSPEC_CPML),stat=ier)
  if (ier /= 0) stop 'error allocating array CPML_type'

  if (ELASTIC_SIMULATION) then
    ! store the displ field at n-1 time step
    allocate(PML_displ_old(3,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating PML_displ_old array'
    ! store the displ field at n time step
    allocate(PML_displ_new(3,NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating PML_displ_new array'
    if (ier /= 0) stop 'error allocating displ_new array'
    ! stores derivatives of ux, uy and uz with respect to x, y and z
    allocate(PML_dux_dxl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dxl array'
    allocate(PML_dux_dyl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dyl array'
    allocate(PML_dux_dzl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dzl array'
    allocate(PML_duy_dxl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dxl array'
    allocate(PML_duy_dyl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dyl array'
    allocate(PML_duy_dzl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dzl array'
    allocate(PML_duz_dxl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dxl array'
    allocate(PML_duz_dyl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dyl array'
    allocate(PML_duz_dzl(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dzl array'

    allocate(PML_dux_dxl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dxl_old array'
    allocate(PML_dux_dyl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dyl_old array'
    allocate(PML_dux_dzl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dzl_old array'
    allocate(PML_duy_dxl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dxl_old array'
    allocate(PML_duy_dyl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dyl_old array'
    allocate(PML_duy_dzl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dzl_old array'
    allocate(PML_duz_dxl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dxl_old array'
    allocate(PML_duz_dyl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dyl_old array'
    allocate(PML_duz_dzl_old(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dzl_old array'

    allocate(PML_dux_dxl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dxl_new array'
    allocate(PML_dux_dyl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dyl_new array'
    allocate(PML_dux_dzl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_dux_dzl_new array'
    allocate(PML_duy_dxl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dxl_new array'
    allocate(PML_duy_dyl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dyl_new array'
    allocate(PML_duy_dzl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duy_dzl_new array'
    allocate(PML_duz_dxl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dxl_new array'
    allocate(PML_duz_dyl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dyl_new array'
    allocate(PML_duz_dzl_new(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating PML_duz_dzl_new array'

    ! stores C-PML memory variables
    allocate(rmemory_dux_dxl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dxl_x array'
    allocate(rmemory_dux_dyl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dyl_x array'
    allocate(rmemory_dux_dzl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dzl_x array'
    allocate(rmemory_duy_dxl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dxl_x array'
    allocate(rmemory_duy_dyl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dyl_x array'
    allocate(rmemory_duz_dxl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dxl_x array'
    allocate(rmemory_duz_dzl_x(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dzl_x array'

    allocate(rmemory_dux_dxl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dxl_y array'
    allocate(rmemory_dux_dyl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dyl_y array'
    allocate(rmemory_duy_dxl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dxl_y array'
    allocate(rmemory_duy_dyl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dyl_y array'
    allocate(rmemory_duy_dzl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dzl_y array'
    allocate(rmemory_duz_dyl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dyl_y array'
    allocate(rmemory_duz_dzl_y(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dzl_y array'

    allocate(rmemory_dux_dxl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dxl_z array'
    allocate(rmemory_dux_dzl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dux_dzl_z array'
    allocate(rmemory_duy_dyl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dyl_z array'
    allocate(rmemory_duy_dzl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duy_dzl_z array'
    allocate(rmemory_duz_dxl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dxl_z array'
    allocate(rmemory_duz_dyl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dyl_z array'
    allocate(rmemory_duz_dzl_z(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_duz_dzl_z array'

    ! stores C-PML memory variables needed for displacement
    allocate(rmemory_displ_elastic(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_displ_elastic array'

    ! stores C-PML contribution to update acceleration to the global mesh
    allocate(accel_elastic_CPML(NDIM,NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating accel_elastic_CPML array'
  endif

  if (ACOUSTIC_SIMULATION) then
    ! store the potential acoustic field at n-1 time step for CMPL
    allocate(PML_potential_acoustic_old(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating PML_potential_acoustic_old array'
    ! store the potential acoustic field at n time step for CMPL
    allocate(PML_potential_acoustic_new(NGLLX,NGLLY,NGLLZ,NSPEC_CPML),stat=ier)
    if (ier /= 0) stop 'error allocating PML_potential_acoustic_new array'

    ! store the potential acoustic field at n-1 time step
!    allocate(potential_dot_dot_acoustic_old(NGLOB_AB),stat=ier)
!    if (ier /= 0) stop 'error allocating potential_dot_dot_acoustic_old array'

    ! stores C-PML memory variables
    allocate(rmemory_dpotential_dxl(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dpotential_dxl array'
    allocate(rmemory_dpotential_dyl(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dpotential_dyl array'
    allocate(rmemory_dpotential_dzl(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_dpotential_dzl array'

    ! stores C-PML memory variables needed for potential
    allocate(rmemory_potential_acoustic(NGLLX,NGLLY,NGLLZ,NSPEC_CPML,3),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_potential_acoustic array'

    ! stores C-PML contribution to update the second derivative of the potential to the global mesh
    allocate(potential_dot_dot_acoustic_CPML(NGLLX,NGLLY,NGLLZ),stat=ier)
    if (ier /= 0) stop 'error allocating potential_dot_dot_acoustic_CPML array'
  endif

  ! stores C-PML contribution on elastic/acoustic interface
  if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
    allocate(rmemory_coupling_ac_el_displ(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_coupling_ac_el_displ array'
    allocate(rmemory_coupling_el_ac_potential_dot_dot(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2),stat=ier)
    if (ier /= 0) stop 'error allocating rmemory_coupling_el_ac_potential_dot_dot array'
    if (SIMULATION_TYPE == 3) then
      allocate(rmemory_coupling_el_ac_potential(3,NGLLX,NGLLY,NGLLZ,num_coupling_ac_el_faces,2),stat=ier)
      if (ier /= 0) stop 'error allocating rmemory_coupling_el_ac_potential array'
    endif
  endif

  ! initializes arrays
  spec_to_CPML(:) = 0
  CPML_type(:) = 0

  if (ELASTIC_SIMULATION) then
    PML_displ_old(:,:,:,:,:) = 0._CUSTOM_REAL
    PML_displ_new(:,:,:,:,:) = 0._CUSTOM_REAL

    PML_dux_dxl(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dyl(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dzl(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dxl(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dyl(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dzl(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dxl(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dyl(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dzl(:,:,:) = 0._CUSTOM_REAL

    PML_dux_dxl_old(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dyl_old(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dzl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dxl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dyl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dzl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dxl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dyl_old(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dzl_old(:,:,:) = 0._CUSTOM_REAL

    PML_dux_dxl_new(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dyl_new(:,:,:) = 0._CUSTOM_REAL
    PML_dux_dzl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dxl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dyl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duy_dzl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dxl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dyl_new(:,:,:) = 0._CUSTOM_REAL
    PML_duz_dzl_new(:,:,:) = 0._CUSTOM_REAL

    rmemory_dux_dxl_x(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dux_dyl_x(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dux_dzl_x(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dxl_x(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dyl_x(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dxl_x(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dzl_x(:,:,:,:) = 0._CUSTOM_REAL

    rmemory_dux_dxl_y(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dux_dyl_y(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dxl_y(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dyl_y(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dzl_y(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dyl_y(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dzl_y(:,:,:,:) = 0._CUSTOM_REAL

    rmemory_dux_dxl_z(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dux_dzl_z(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dyl_z(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duy_dzl_z(:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dxl_z(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dyl_z(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_duz_dzl_z(:,:,:,:,:) = 0._CUSTOM_REAL

    rmemory_displ_elastic(:,:,:,:,:,:) = 0._CUSTOM_REAL
    accel_elastic_CPML(:,:,:,:) = 0._CUSTOM_REAL
  endif

  if (ACOUSTIC_SIMULATION) then
    PML_potential_acoustic_old(:,:,:,:) = 0._CUSTOM_REAL
    PML_potential_acoustic_new(:,:,:,:) = 0._CUSTOM_REAL

!    potential_dot_dot_acoustic_old(:) = 0._CUSTOM_REAL

    rmemory_dpotential_dxl(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dpotential_dyl(:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_dpotential_dzl(:,:,:,:,:) = 0._CUSTOM_REAL

    rmemory_potential_acoustic(:,:,:,:,:) = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_CPML(:,:,:) = 0._CUSTOM_REAL
  endif

  if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
    rmemory_coupling_ac_el_displ(:,:,:,:,:,:) = 0._CUSTOM_REAL
    rmemory_coupling_el_ac_potential_dot_dot(:,:,:,:,:,:) = 0._CUSTOM_REAL
  endif

  if (SIMULATION_TYPE == 3) then
    if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
      rmemory_coupling_el_ac_potential(:,:,:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

! fields on PML interface will be reconstructed for adjoint simulations
! using snapshot files of wavefields

! elastic domains
  if (ELASTIC_SIMULATION) then
    ! opens absorbing wavefield saved/to-be-saved by forward simulations
    if (nglob_interface_PML_elastic > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

      b_nglob_interface_PML_elastic = nglob_interface_PML_elastic

      ! allocates wavefield
      allocate(b_PML_field(9,b_nglob_interface_PML_elastic),stat=ier)
      if (ier /= 0) stop 'error allocating array b_PML_field'

      ! size of single record
      b_reclen_PML_field = CUSTOM_REAL * 9 * nglob_interface_PML_elastic

      ! check integer size limit: size of b_reclen_PML_field must fit onto an 4-byte integer
      if (nglob_interface_PML_elastic > 2147483646 / (CUSTOM_REAL * 9)) then
        print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_PML_field
        print *,'  ',CUSTOM_REAL, NDIM, 9, nglob_interface_PML_elastic
        print *,'bit size fortran: ',bit_size(b_reclen_PML_field)
        call exit_MPI(myrank,"error b_reclen_PML_field integer limit")
      endif

      ! total file size
      filesize = b_reclen_PML_field
      filesize = filesize*NSTEP

      if (SIMULATION_TYPE == 3) then
        call open_file_abs_r(0,trim(prname)//'absorb_PML_field.bin', &
                            len_trim(trim(prname)//'absorb_PML_field.bin'), &
                            filesize)

      else
        call open_file_abs_w(0,trim(prname)//'absorb_PML_field.bin', &
                            len_trim(trim(prname)//'absorb_PML_field.bin'), &
                            filesize)
      endif
    else
      ! needs dummy array
      allocate(b_PML_field(9,1),stat=ier)
      if (ier /= 0) stop 'error allocating array b_PML_field'
    endif
  endif

! acoustic domains
  if (ACOUSTIC_SIMULATION) then
    ! opens absorbing wavefield saved/to-be-saved by forward simulations
    if (nglob_interface_PML_acoustic > 0 .and. (SIMULATION_TYPE == 3 .or. &
        (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

      b_nglob_interface_PML_acoustic = nglob_interface_PML_acoustic

      ! allocates wavefield
      allocate(b_PML_potential(3,b_nglob_interface_PML_acoustic),stat=ier)
      if (ier /= 0) stop 'error allocating array b_PML_potential'

      ! size of single record
      b_reclen_PML_potential = CUSTOM_REAL * nglob_interface_PML_acoustic

      ! check integer size limit: size of b_reclen_PML_field must fit onto an 4-byte integer
      if (nglob_interface_PML_acoustic > 2147483646 / (CUSTOM_REAL)) then
        print *,'reclen needed exceeds integer 4-byte limit: ',b_reclen_PML_potential
        print *,'  ',CUSTOM_REAL, nglob_interface_PML_acoustic
        print *,'bit size fortran: ',bit_size(b_reclen_PML_potential)
        call exit_MPI(myrank,"error b_reclen_PML_potential integer limit")
      endif

      ! total file size
      filesize = b_reclen_PML_potential
      filesize = filesize*NSTEP

      if (SIMULATION_TYPE == 3) then
        call open_file_abs_r(1,trim(prname)//'absorb_PML_potential.bin', &
                            len_trim(trim(prname)//'absorb_PML_potential.bin'), &
                            filesize)

      else
        call open_file_abs_w(1,trim(prname)//'absorb_PML_potential.bin', &
                            len_trim(trim(prname)//'absorb_PML_potential.bin'), &
                            filesize)
      endif
    else
      ! needs dummy array
      allocate(b_PML_potential(3,1),stat=ier)
      if (ier /= 0) stop 'error allocating array b_PML_potential'
    endif
  endif

  end subroutine pml_allocate_arrays

!=====================================================================

  subroutine pml_allocate_arrays_dummy()

  ! dummy allocation with a size of 1 for all the PML arrays that have not yet been allocated
  ! in order to be able to use these arrays as arguments in subroutine calls

  use specfem_par, only: SIMULATION_TYPE
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION

  use pml_par

  implicit none

  if (.not. allocated(spec_to_CPML)) allocate(spec_to_CPML(1))
  if (.not. allocated(CPML_type)) allocate(CPML_type(1))

  if (ELASTIC_SIMULATION) then
    if (.not. allocated(PML_displ_old)) allocate(PML_displ_old(3,1,1,1,1))
    if (.not. allocated(PML_displ_new)) allocate(PML_displ_new(3,1,1,1,1))

    if (.not. allocated(PML_dux_dxl)) allocate(PML_dux_dxl(1,1,1))
    if (.not. allocated(PML_dux_dyl)) allocate(PML_dux_dyl(1,1,1))
    if (.not. allocated(PML_dux_dzl)) allocate(PML_dux_dzl(1,1,1))
    if (.not. allocated(PML_duy_dxl)) allocate(PML_duy_dxl(1,1,1))
    if (.not. allocated(PML_duy_dyl)) allocate(PML_duy_dyl(1,1,1))
    if (.not. allocated(PML_duy_dzl)) allocate(PML_duy_dzl(1,1,1))
    if (.not. allocated(PML_duz_dxl)) allocate(PML_duz_dxl(1,1,1))
    if (.not. allocated(PML_duz_dyl)) allocate(PML_duz_dyl(1,1,1))
    if (.not. allocated(PML_duz_dzl)) allocate(PML_duz_dzl(1,1,1))
    if (.not. allocated(PML_dux_dxl_old)) allocate(PML_dux_dxl_old(1,1,1))
    if (.not. allocated(PML_dux_dyl_old)) allocate(PML_dux_dyl_old(1,1,1))
    if (.not. allocated(PML_dux_dzl_old)) allocate(PML_dux_dzl_old(1,1,1))
    if (.not. allocated(PML_duy_dxl_old)) allocate(PML_duy_dxl_old(1,1,1))
    if (.not. allocated(PML_duy_dyl_old)) allocate(PML_duy_dyl_old(1,1,1))
    if (.not. allocated(PML_duy_dzl_old)) allocate(PML_duy_dzl_old(1,1,1))
    if (.not. allocated(PML_duz_dxl_old)) allocate(PML_duz_dxl_old(1,1,1))
    if (.not. allocated(PML_duz_dyl_old)) allocate(PML_duz_dyl_old(1,1,1))
    if (.not. allocated(PML_duz_dzl_old)) allocate(PML_duz_dzl_old(1,1,1))
    if (.not. allocated(PML_dux_dxl_new)) allocate(PML_dux_dxl_new(1,1,1))
    if (.not. allocated(PML_dux_dyl_new)) allocate(PML_dux_dyl_new(1,1,1))
    if (.not. allocated(PML_dux_dzl_new)) allocate(PML_dux_dzl_new(1,1,1))
    if (.not. allocated(PML_duy_dxl_new)) allocate(PML_duy_dxl_new(1,1,1))
    if (.not. allocated(PML_duy_dyl_new)) allocate(PML_duy_dyl_new(1,1,1))
    if (.not. allocated(PML_duy_dzl_new)) allocate(PML_duy_dzl_new(1,1,1))
    if (.not. allocated(PML_duz_dxl_new)) allocate(PML_duz_dxl_new(1,1,1))
    if (.not. allocated(PML_duz_dyl_new)) allocate(PML_duz_dyl_new(1,1,1))
    if (.not. allocated(PML_duz_dzl_new)) allocate(PML_duz_dzl_new(1,1,1))
    if (.not. allocated(rmemory_dux_dxl_x)) allocate(rmemory_dux_dxl_x(1,1,1,1,3))
    if (.not. allocated(rmemory_dux_dyl_x)) allocate(rmemory_dux_dyl_x(1,1,1,1,3))
    if (.not. allocated(rmemory_dux_dzl_x)) allocate(rmemory_dux_dzl_x(1,1,1,1,3))
    if (.not. allocated(rmemory_duy_dxl_x)) allocate(rmemory_duy_dxl_x(1,1,1,1))
    if (.not. allocated(rmemory_duy_dyl_x)) allocate(rmemory_duy_dyl_x(1,1,1,1))
    if (.not. allocated(rmemory_duz_dxl_x)) allocate(rmemory_duz_dxl_x(1,1,1,1))
    if (.not. allocated(rmemory_duz_dzl_x)) allocate(rmemory_duz_dzl_x(1,1,1,1))
    if (.not. allocated(rmemory_dux_dxl_y)) allocate(rmemory_dux_dxl_y(1,1,1,1))
    if (.not. allocated(rmemory_dux_dyl_y)) allocate(rmemory_dux_dyl_y(1,1,1,1))
    if (.not. allocated(rmemory_duy_dxl_y)) allocate(rmemory_duy_dxl_y(1,1,1,1,3))
    if (.not. allocated(rmemory_duy_dyl_y)) allocate(rmemory_duy_dyl_y(1,1,1,1,3))
    if (.not. allocated(rmemory_duy_dzl_y)) allocate(rmemory_duy_dzl_y(1,1,1,1,3))
    if (.not. allocated(rmemory_duz_dyl_y)) allocate(rmemory_duz_dyl_y(1,1,1,1))
    if (.not. allocated(rmemory_duz_dzl_y)) allocate(rmemory_duz_dzl_y(1,1,1,1))
    if (.not. allocated(rmemory_dux_dxl_z)) allocate(rmemory_dux_dxl_z(1,1,1,1))
    if (.not. allocated(rmemory_dux_dzl_z)) allocate(rmemory_dux_dzl_z(1,1,1,1))
    if (.not. allocated(rmemory_duy_dyl_z)) allocate(rmemory_duy_dyl_z(1,1,1,1))
    if (.not. allocated(rmemory_duy_dzl_z)) allocate(rmemory_duy_dzl_z(1,1,1,1))
    if (.not. allocated(rmemory_duz_dxl_z)) allocate(rmemory_duz_dxl_z(1,1,1,1,3))
    if (.not. allocated(rmemory_duz_dyl_z)) allocate(rmemory_duz_dyl_z(1,1,1,1,3))
    if (.not. allocated(rmemory_duz_dzl_z)) allocate(rmemory_duz_dzl_z(1,1,1,1,3))
    if (.not. allocated(rmemory_displ_elastic)) allocate(rmemory_displ_elastic(1,1,1,1,1,3))
    if (.not. allocated(accel_elastic_CPML)) allocate(accel_elastic_CPML(1,1,1,1))

    ! allocates wavefield
    if (.not. allocated(b_PML_field)) allocate(b_PML_field(9,1))
  endif

  if (ACOUSTIC_SIMULATION) then
    if (.not. allocated(PML_potential_acoustic_old)) allocate(PML_potential_acoustic_old(1,1,1,1))
    if (.not. allocated(PML_potential_acoustic_new)) allocate(PML_potential_acoustic_new(1,1,1,1))

    if (.not. allocated(rmemory_dpotential_dxl)) allocate(rmemory_dpotential_dxl(1,1,1,1,3))
    if (.not. allocated(rmemory_dpotential_dyl)) allocate(rmemory_dpotential_dyl(1,1,1,1,3))
    if (.not. allocated(rmemory_dpotential_dzl)) allocate(rmemory_dpotential_dzl(1,1,1,1,3))
    if (.not. allocated(rmemory_potential_acoustic)) allocate(rmemory_potential_acoustic(1,1,1,1,3))
    if (.not. allocated(potential_dot_dot_acoustic_CPML)) allocate(potential_dot_dot_acoustic_CPML(1,1,1))

    ! allocates wavefield
    if (.not. allocated(b_PML_potential)) allocate(b_PML_potential(3,1))
  endif

  if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
    if (.not. allocated(rmemory_coupling_ac_el_displ)) allocate(rmemory_coupling_ac_el_displ(3,1,1,1,1,2))
    if (.not. allocated(rmemory_coupling_el_ac_potential_dot_dot)) allocate(rmemory_coupling_el_ac_potential_dot_dot(3,1,1,1,1,2))
    if (SIMULATION_TYPE == 3) then
      if (.not. allocated(rmemory_coupling_el_ac_potential)) allocate(rmemory_coupling_el_ac_potential(3,1,1,1,1,2))
    endif
  endif

  end subroutine pml_allocate_arrays_dummy

!=====================================================================

  subroutine pml_cleanup()

! deallocates C_PML arrays

  use specfem_par, only: SIMULATION_TYPE
  use specfem_par_acoustic, only: ACOUSTIC_SIMULATION
  use specfem_par_elastic, only: ELASTIC_SIMULATION

  use pml_par

  implicit none

  deallocate(is_CPML)

  if (NSPEC_CPML > 0) then
    deallocate(CPML_regions)
    deallocate(CPML_to_spec)

    deallocate(d_store_x)
    deallocate(d_store_y)
    deallocate(d_store_z)
    deallocate(k_store_x)
    deallocate(k_store_y)
    deallocate(k_store_z)
    deallocate(alpha_store_x)
    deallocate(alpha_store_y)
    deallocate(alpha_store_z)
  endif

  deallocate(spec_to_CPML)
  deallocate(CPML_type)

  if (ELASTIC_SIMULATION) then
    deallocate(PML_displ_old)
    deallocate(PML_displ_new)

    deallocate(PML_dux_dxl)
    deallocate(PML_dux_dyl)
    deallocate(PML_dux_dzl)
    deallocate(PML_duy_dxl)
    deallocate(PML_duy_dyl)
    deallocate(PML_duy_dzl)
    deallocate(PML_duz_dxl)
    deallocate(PML_duz_dyl)
    deallocate(PML_duz_dzl)
    deallocate(PML_dux_dxl_old)
    deallocate(PML_dux_dyl_old)
    deallocate(PML_dux_dzl_old)
    deallocate(PML_duy_dxl_old)
    deallocate(PML_duy_dyl_old)
    deallocate(PML_duy_dzl_old)
    deallocate(PML_duz_dxl_old)
    deallocate(PML_duz_dyl_old)
    deallocate(PML_duz_dzl_old)
    deallocate(PML_dux_dxl_new)
    deallocate(PML_dux_dyl_new)
    deallocate(PML_dux_dzl_new)
    deallocate(PML_duy_dxl_new)
    deallocate(PML_duy_dyl_new)
    deallocate(PML_duy_dzl_new)
    deallocate(PML_duz_dxl_new)
    deallocate(PML_duz_dyl_new)
    deallocate(PML_duz_dzl_new)
    deallocate(rmemory_dux_dxl_x)
    deallocate(rmemory_dux_dyl_x)
    deallocate(rmemory_dux_dzl_x)
    deallocate(rmemory_duy_dxl_x)
    deallocate(rmemory_duy_dyl_x)
    deallocate(rmemory_duz_dxl_x)
    deallocate(rmemory_duz_dzl_x)
    deallocate(rmemory_dux_dxl_y)
    deallocate(rmemory_dux_dyl_y)
    deallocate(rmemory_duy_dxl_y)
    deallocate(rmemory_duy_dyl_y)
    deallocate(rmemory_duy_dzl_y)
    deallocate(rmemory_duz_dyl_y)
    deallocate(rmemory_duz_dzl_y)
    deallocate(rmemory_dux_dxl_z)
    deallocate(rmemory_dux_dzl_z)
    deallocate(rmemory_duy_dyl_z)
    deallocate(rmemory_duy_dzl_z)
    deallocate(rmemory_duz_dxl_z)
    deallocate(rmemory_duz_dyl_z)
    deallocate(rmemory_duz_dzl_z)
    deallocate(rmemory_displ_elastic)
    deallocate(accel_elastic_CPML)
  endif

  if (ACOUSTIC_SIMULATION) then
    deallocate(PML_potential_acoustic_old)
    deallocate(PML_potential_acoustic_new)

    deallocate(rmemory_dpotential_dxl)
    deallocate(rmemory_dpotential_dyl)
    deallocate(rmemory_dpotential_dzl)
    deallocate(rmemory_potential_acoustic)
    deallocate(potential_dot_dot_acoustic_CPML)
  endif

  if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
    deallocate(rmemory_coupling_ac_el_displ)
    deallocate(rmemory_coupling_el_ac_potential_dot_dot)
    if (SIMULATION_TYPE == 3) deallocate(rmemory_coupling_el_ac_potential)
  endif

  end subroutine pml_cleanup
