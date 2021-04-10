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

module model_ipati_adios_mod

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,IMAIN,FOUR_THIRDS

contains

!-----------------------------------------------------------------------------
!
! IPATI
!
! based on given rho and vp structure for GLL files
!
!------------------------------------------------------------------------------

  subroutine model_ipati_adios(myrank,nspec,LOCAL_PATH)

  use adios_read_mod
  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier

  ! ---------------------------------------------------------------------------
  ! note: vp not vs structure is available
  ! (as is often the case in exploration seismology),
  ! scaling factor
  real, parameter :: SCALING_FACTOR = 1.0/1.8
  ! ---------------------------------------------------------------------------

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using IPATI model (ADIOS) from: ',trim(LOCAL_PATH)
    write(IMAIN,*) '     scaling factor: ',SCALING_FACTOR
  endif

  ! density
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 572')
  if (ier /= 0) stop 'error allocating array rho_read'
  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 573')
  if (ier /= 0) stop 'error allocating array vp_read'
  ! vs scaled from vp
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 574')
  if (ier /= 0) stop 'error allocating array vs_read'

  call read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! isotropic model parameters
  rhostore    = rho_read
  kappastore  = rhostore * ( vp_read * vp_read - FOUR_THIRDS * vs_read * vs_read )
  mustore     = rhostore * vs_read * vs_read
  rho_vp = rhostore * vp_read
  rho_vs = rhostore * vs_read

  ! free memory
  deallocate( rho_read,vp_read,vs_read)

  end subroutine model_ipati_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_ipati_water_adios(myrank,nspec,LOCAL_PATH)

  use create_regions_mesh_ext_par
  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ispec,ier

  ! -----------------------------------------------------------------------------

  ! note: vp not vs structure is available (as is often the case in exploration seismology),
  ! scaling factor
  real, parameter :: SCALING_FACTOR = 1.0/1.8

  ! -----------------------------------------------------------------------------

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using IPATI_WATER model (ADIOS) from: ',trim(LOCAL_PATH)
    write(IMAIN,*) '     scaling factor: ',SCALING_FACTOR
  endif

  ! density
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 575')
  if (ier /= 0) stop 'error allocating array rho_read'
  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 576')
  if (ier /= 0) stop 'error allocating array vp_read'
  ! vs scaled from vp
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 577')
  if (ier /= 0) stop 'error allocating array vs_read'

  call read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! overwrites only elastic elements
  do ispec=1,nspec
    ! assumes water layer with acoustic elements are set properly
    ! only overwrites elastic elements
    if (ispec_is_elastic(ispec)) then
      ! isotropic model parameters
      rhostore(:,:,:,ispec) = rho_read(:,:,:,ispec)
      kappastore(:,:,:,ispec) = rhostore(:,:,:,ispec) * ( vp_read(:,:,:,ispec) * vp_read(:,:,:,ispec) &
                                    - FOUR_THIRDS * vs_read(:,:,:,ispec) * vs_read(:,:,:,ispec) )
      mustore(:,:,:,ispec) = rhostore(:,:,:,ispec) * vs_read(:,:,:,ispec) * vs_read(:,:,:,ispec)
      rho_vp(:,:,:,ispec) = rhostore(:,:,:,ispec) * vp_read(:,:,:,ispec)
      rho_vs(:,:,:,ispec) = rhostore(:,:,:,ispec) * vs_read(:,:,:,ispec)
    endif
  enddo

  ! free memory
  deallocate( rho_read,vp_read,vs_read)

  end subroutine model_ipati_water_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_vp_rho_adios (myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  use adios_read_mod
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH
  real, dimension(:,:,:,:), intent(inout) :: vp_read,rho_read

  ! ADIOS stuffs
  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle, sel
  integer(kind=8), dimension(1) :: start, count_ad
  integer :: local_dim_rho, local_dim_vp
  integer :: ier
  integer :: comm

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = LOCAL_PATH(1:len_trim(LOCAL_PATH)) //"/model_values.bp"

  ! gets MPI communicator
  comm = comm_adios

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, ADIOS_VERBOSITY, ier)
  call adios_read_open_file (handle, database_name, 0, comm, ier)
  if (ier /= 0) call abort_mpi()

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  call adios_get_scalar(handle, "rho/local_dim", local_dim_rho, ier)
  call adios_get_scalar(handle, "vp/local_dim", local_dim_vp, ier)

  start(1) = local_dim_rho * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox(sel, 1, start, count_ad)
  call adios_schedule_read(handle, sel, "rho/array", 0, 1, &
                           rho_read, ier)
  call adios_schedule_read(handle, sel, "vp/array", 0, 1, &
                           vp_read, ier)

  !---------------------------------------.
  ! Perform read and close the adios file |
  !---------------------------------------'
  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  end subroutine read_model_vp_rho_adios

end module model_ipati_adios_mod
