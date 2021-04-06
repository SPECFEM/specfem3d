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

!==============================================================================
!> Read rho, vp, vs from model_values.bp
!
! based on modified GLL mesh output from mesher
! used for iterative inversion procedures
!
! \param myrank rank of the MPI process
! \param nspec  number of spectral elements in the model
! \param LOCAL_PATH path where the '.bp' file is located
subroutine model_gll_adios(myrank,nspec,LOCAL_PATH)

  use adios_read_mod
  use adios_manager_mod, only: comm_adios,ADIOS_VERBOSITY

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,ATTENUATION

  use create_regions_mesh_ext_par

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier

  ! ADIOS stuffs
  character(len=MAX_STRING_LEN) :: database_name
  integer(kind=8) :: handle, sel
  integer(kind=8), dimension(1) :: start, count_ad
  integer :: local_dim_rho, local_dim_vp, local_dim_vs
  integer :: comm

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model (ADIOS) from: ',trim(LOCAL_PATH)
  endif

  ! density
  allocate( rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 617')
  if (ier /= 0) stop 'error allocating array rho_read'
  ! vp
  allocate( vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 618')
  if (ier /= 0) stop 'error allocating array vp_read'
  ! vs
  allocate( vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 619')
  if (ier /= 0) stop 'error allocating array vs_read'

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
  call adios_get_scalar(handle, "vs/local_dim", local_dim_vs, ier)

  start(1) = local_dim_rho * myrank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox(sel, 1, start, count_ad)

  ! wave speeds
  call adios_schedule_read(handle, sel, "rho/array", 0, 1, &
                           rho_read, ier)
  call adios_schedule_read(handle, sel, "vp/array", 0, 1, &
                           vp_read, ier)
  call adios_schedule_read(handle, sel, "vp/array", 0, 1, &
                           vp_read, ier)

  ! attenuation
  if (ATTENUATION) then
    ! shear attenuation
    call adios_schedule_read(handle, sel, "qmu/array", 0, 1, &
                             qmu_attenuation_store, ier)
    ! bulk attenuation
    call adios_schedule_read(handle, sel, "qkappa/array", 0, 1, &
                             qkappa_attenuation_store, ier)
  endif

  !---------------------------------------.
  ! Perform read and close the adios file |
  !---------------------------------------'
  call adios_perform_reads(handle, ier)
  if (ier /= 0) call abort_mpi()

  call adios_read_close(handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error adios read finalize'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where density structure is not given
  !!! modify according to your desire

  !  rho_read = 1000.0
  !  where ( mustore > 100.0 )  &
  !           rho_read = (1.6612 * (vp_read / 1000.0)     &
  !                      -0.4720 * (vp_read / 1000.0)**2  &
  !                      +0.0671 * (vp_read / 1000.0)**3  &
  !                      -0.0043 * (vp_read / 1000.0)**4  &
  !                      +0.000106*(vp_read / 1000.0)**5)*1000.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! in cases where shear wavespeed structure is not given
  !!! modify according to your desire

  !   vs_read = 0.0
  !   where ( mustore > 100.0 )       vs_read = vp_read / sqrt(3.0)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! update arrays that will be saved and used in the solver xspecfem3D
  !!! the following part is neccessary if you uncommented something above

  rhostore(:,:,:,:) = rho_read(:,:,:,:)
  kappastore(:,:,:,:) = rhostore(:,:,:,:) * ( vp_read(:,:,:,:) * vp_read(:,:,:,:) &
                          - FOUR_THIRDS * vs_read(:,:,:,:) * vs_read(:,:,:,:) )
  mustore(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:) * vs_read(:,:,:,:)
  rho_vp(:,:,:,:) = rhostore(:,:,:,:) * vp_read(:,:,:,:)
  rho_vs(:,:,:,:) = rhostore(:,:,:,:) * vs_read(:,:,:,:)

  ! free memory
  deallocate( rho_read,vp_read,vs_read)

  end subroutine model_gll_adios
