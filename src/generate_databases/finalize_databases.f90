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
!

  subroutine finalize_databases()

! finalizes generate databases

  use generate_databases_par

  implicit none

  ! local parameters
  ! timing
  double precision, external :: wtime

  ! outputs partitioning info
  call finalize_repartition_info()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! copy number of elements and points in an include file for the solver
  if (myrank == 0) then
    call save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
                          ATTENUATION,ANISOTROPY,NSTEP,DT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                          SIMULATION_TYPE,max_memory_size,nfaces_surface_glob_ext_mesh)
  endif

  ! elapsed time since beginning of mesh generation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
  endif

  ! close main output file
  if (myrank == 0) then
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine finalize_databases

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_repartition_info()

! this routine is similar to routine print_statistics_load() for xdecompose_mesh

  use generate_databases_par

  use create_regions_mesh_ext_par, only: nspec_irregular,ispec_is_acoustic,ispec_is_elastic

  implicit none

  ! local parameters
  integer :: nspec_total
  ! this can overflow if more than 2 Gigapoints in the whole mesh, thus replaced with double precision version
  integer(kind=8) :: nglob_total
  double precision :: nglob_l,nglob_total_db

  ! loads
  ! same integer loads as in part_decompose_mesh
  integer, parameter :: ACOUSTIC_LOAD = 10     ! is in reality 1.0
  integer, parameter :: ELASTIC_LOAD = 41      ! is in reality 4.1
  integer, parameter :: VISCOELASTIC_LOAD = 59 ! is in reality 5.9
  integer, parameter :: POROELASTIC_LOAD = 81  ! is in reality 8.1

  integer, dimension(:), allocatable  :: elmnts_load
  integer, dimension(:), allocatable  :: gather_loads
  integer :: load_min,load_max,load_per_proc
  integer :: iproc,ispec,ier
  real :: load_balance

  ! elements load array
  allocate(elmnts_load(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 113')
  if (ier /= 0) stop 'Error allocating array elmnts_load'
  !! DK DK Oct 2012: this should include CPML weights as well in the future
  ! uniform load by default
  elmnts_load(:) = ACOUSTIC_LOAD

  ! determine loads based on element type
  do ispec = 1,NSPEC_AB
    if (ispec_is_acoustic(ispec)) then
      ! acoustic
      elmnts_load(ispec) = ACOUSTIC_LOAD
    else if (ispec_is_elastic(ispec)) then
      ! elastic
      if (ATTENUATION) then
        elmnts_load(ispec) = VISCOELASTIC_LOAD
      else
        elmnts_load(ispec) = ELASTIC_LOAD
      endif
    else
      ! poroelastic
      elmnts_load(ispec) = POROELASTIC_LOAD
    endif
  enddo

  ! stats
  ! check integer size limit
  load_min = NSPEC_AB * maxval(elmnts_load(:))
  if (NSPEC_AB > int(2147483646.0 / maxval(elmnts_load(:)))) then
    load_min = 2147483646
  endif
  load_max = 0

  ! determine load per process
  allocate(gather_loads(0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 113')
  if (ier /= 0) stop 'Error allocating array elmnts_load'
  gather_loads(:) = 0

  ! counts loads in this slice
  load_per_proc = 0
  do ispec = 1,NSPEC_AB
    ! load per element
    load_per_proc = load_per_proc + elmnts_load(ispec)
  enddo

  ! collect loads from all processes on master
  call gather_all_singlei(load_per_proc,gather_loads,NPROC)

  ! print number of points and elements in the mesh
  call sum_all_i(NSPEC_AB,nspec_total)

  ! this can overflow if more than 2 Gigapoints in the whole mesh, thus replaced with double precision version
  nglob_l = dble(NGLOB_AB)
  call sum_all_dp(nglob_l,nglob_total_db)

  ! user output
  if (myrank == 0) then
    ! converts to integer*8
    ! note: only the main process has the total sum in nglob_total_db and a valid value;
    !       this conversion could lead to compiler errors if done by other processes.
    nglob_total = int(nglob_total_db,kind=8)

    ! load-balancing
    ! determine min/max loads over all processes
    load_min = minval(gather_loads(:))
    load_max = maxval(gather_loads(:))
    ! computes load-balance
    if (load_max > 0) then
      ! imbalance in percent
      load_balance = 100.0 * (load_max - load_min ) / load_max
    else
      stop 'Error load-balance: partitioning has zero load'
    endif

    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements:'
    write(IMAIN,*) '-----------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'load distribution:'
    write(IMAIN,*) '  element loads: min/max = ',load_min,load_max
    write(IMAIN,*)
    do iproc = 0,NPROC-1
      load_per_proc = gather_loads(iproc)
      write(IMAIN,*) '  partition ',iproc, '       has ',load_per_proc,' load units'
    enddo
    write(IMAIN,*)
    write(IMAIN,*) '  load per partition: min/max   = ',load_min,load_max
    write(IMAIN,*) '  load per partition: imbalance = ',load_balance,'%'
    write(IMAIN,*)'                      (0% being totally balanced, 100% being unbalanced)'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in mesh slice 0: ',NSPEC_AB
    write(IMAIN,*) 'total number of   regular elements in mesh slice 0: ',NSPEC_AB - nspec_irregular
    write(IMAIN,*) 'total number of irregular elements in mesh slice 0: ',nspec_irregular
    write(IMAIN,*) 'total number of points in mesh slice 0: ',NGLOB_AB
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in entire mesh: ',nspec_total
    write(IMAIN,*) 'approximate total number of points in entire mesh (with duplicates on MPI edges): ',nglob_total
    write(IMAIN,*) 'approximate total number of DOFs   in entire mesh (with duplicates on MPI edges): ',nglob_total*NDIM
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine finalize_repartition_info
