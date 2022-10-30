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


!-----------------------------------------------------------------------------
!
! IPATI
!
! based on given rho and vp structure for GLL files
!
!------------------------------------------------------------------------------

  subroutine model_ipati_adios(myrank,nspec,LOCAL_PATH)

  use constants, only: NGLLX,NGLLY,NGLLZ,IMAIN,FOUR_THIRDS
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
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 572')
  if (ier /= 0) stop 'error allocating array rho_read'
  rho_read(:,:,:,:) = 0.0
  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 573')
  if (ier /= 0) stop 'error allocating array vp_read'
  vp_read(:,:,:,:) = 0.0
  ! vs scaled from vp
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 574')
  if (ier /= 0) stop 'error allocating array vs_read'
  vs_read(:,:,:,:) = 0.0

  call read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! isotropic model parameters
  rhostore   = rho_read
  kappastore = rhostore * ( vp_read * vp_read - FOUR_THIRDS * vs_read * vs_read )
  mustore    = rhostore * vs_read * vs_read
  rho_vp = rhostore * vp_read
  rho_vs = rhostore * vs_read

  ! free memory
  deallocate(rho_read,vp_read,vs_read)

  end subroutine model_ipati_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_ipati_water_adios(myrank,nspec,LOCAL_PATH)

  use constants, only: NGLLX,NGLLY,NGLLZ,IMAIN,FOUR_THIRDS
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
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 575')
  if (ier /= 0) stop 'error allocating array rho_read'
  rho_read(:,:,:,:) = 0.0
  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 576')
  if (ier /= 0) stop 'error allocating array vp_read'
  vp_read(:,:,:,:) = 0.0
  ! vs scaled from vp
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 577')
  if (ier /= 0) stop 'error allocating array vs_read'
  vs_read(:,:,:,:) = 0.0

  call read_model_vp_rho_adios(myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  ! scaling
  vs_read = vp_read * SCALING_FACTOR

  ! overwrites only elastic elements
  do ispec = 1,nspec
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
  deallocate(rho_read,vp_read,vs_read)

  end subroutine model_ipati_water_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_vp_rho_adios (myrank, nspec, LOCAL_PATH, rho_read, vp_read)

  use constants, only: NGLLX,NGLLY,NGLLZ,IMAIN
  use create_regions_mesh_ext_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH
  real, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(inout) :: vp_read,rho_read

  ! ADIOS stuffs
  character(len=MAX_STRING_LEN) :: database_name

  !-------------------------------------.
  ! Open ADIOS Database file, read mode |
  !-------------------------------------'
  database_name = get_adios_filename(trim(LOCAL_PATH) // "/model_values")

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  model file: ',trim(database_name)
#if defined(USE_ADIOS)
    write(IMAIN,*) '  using ADIOS1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  using ADIOS2 file format'
#endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initiate new group
  call init_adios_group(myadios_group,"ModelIPATIReader")

  ! setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  ! assumes GLL type array size (NGLLX,NGLLY,NGLLZ,nspec)
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, "rho", rho_read)
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, "vp", vp_read)

  !---------------------------------------.
  ! Perform read and close the adios file |
  !---------------------------------------'
  ! explicit read
  ! in principle, already done in read routine, otherwise executed when closing
  call read_adios_perform(myadios_file)

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"ModelIPATIReader")

  end subroutine read_model_vp_rho_adios
