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

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,FOUR_THIRDS,IMAIN,ATTENUATION

  use create_regions_mesh_ext_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer, intent(in) :: myrank,nspec
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  real, dimension(:,:,:,:),allocatable :: vp_read,vs_read,rho_read
  integer :: ier

  ! ADIOS stuffs
  character(len=MAX_STRING_LEN) :: database_name
  !integer(kind=8) :: sel
  !integer(kind=8), dimension(1) :: start, count
  !integer :: local_dim_rho

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     using GLL model (ADIOS) from: ',trim(LOCAL_PATH)
  endif

  ! density
  allocate(rho_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 617')
  if (ier /= 0) stop 'error allocating array rho_read'
  rho_read(:,:,:,:) = 0.0
  ! vp
  allocate(vp_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 618')
  if (ier /= 0) stop 'error allocating array vp_read'
  vp_read(:,:,:,:) = 0.0
  ! vs
  allocate(vs_read(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 619')
  if (ier /= 0) stop 'error allocating array vs_read'
  vs_read(:,:,:,:) = 0.0

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
  call init_adios_group(myadios_group,"ModelGLLReader")

  ! setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,database_name)

  !------------------------.
  ! Get the 'chunks' sizes |
  !------------------------'
  !call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "rho", local_dim_rho)
  !call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "vp", local_dim_vp)
  !call read_adios_scalar_local_dim(myadios_file, myadios_group, myrank, "vs", local_dim_vs)
  ! assuming local_dim is the same for rho,vp,vs,.. arrays
  !start(1) = local_dim_rho * myrank
  !count(1) = NGLLX * NGLLY * NGLLZ * nspec
  !call set_selection_boundingbox(sel, start, count)
  !call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "rho/array", rho_read)
  !call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "vp/array", vp_read)
  !call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "vs/array", vs_read)

  ! wave speeds
  ! assumes GLL type array size (NGLLX,NGLLY,NGLLZ,nspec)
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, "rho", rho_read)
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, "vp", vp_read)
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, "vs", vs_read)

  ! attenuation
  if (ATTENUATION) then
    ! shear attenuation
    !call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "qmu/array", qmu_attenuation_store)
    call read_adios_array(myadios_file, myadios_group, myrank, nspec, "qmu", qmu_attenuation_store)

    ! bulk attenuation
    !call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "qkappa/array", qkappa_attenuation_store)
    call read_adios_array(myadios_file, myadios_group, myrank, nspec, "qkappa", qkappa_attenuation_store)
  endif

  !---------------------------------------.
  ! Perform read and close the adios file |
  !---------------------------------------'
  ! explicit read
  ! in principle, already done in read routine, otherwise executed when closing
  call read_adios_perform(myadios_file)

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"ModelGLLReader")

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
  deallocate(rho_read,vp_read,vs_read)

  end subroutine model_gll_adios
