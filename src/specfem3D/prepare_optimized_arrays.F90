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


! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


  subroutine prepare_optimized_arrays()

! optimizes array memory layout to increase computational efficiency

  use constants, only: myrank,IMAIN

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing optimized arrays"
    call flush_IMAIN()
  endif

#ifdef USE_OPENMP
  ! prepares arrays for OpenMP
  call prepare_timerun_OpenMP()
#endif

  ! prepares element flags
  call prepare_irregular_elements()

  ! prepare fused array for computational kernel
  call prepare_fused_array()

  ! a memory bandwidth benchmark
  call prepare_bandwidth_test()

  end subroutine prepare_optimized_arrays


!
!-------------------------------------------------------------------------------------------------
!

! OpenMP version uses "special" compute_forces_viscoelastic routine
! we need to set num_elem_colors_elastic arrays

#ifdef USE_OPENMP
  subroutine prepare_timerun_OpenMP()

  use specfem_par
  use specfem_par_elastic

  implicit none

  ! local parameters
  integer :: ier
  integer :: NUM_THREADS
  integer :: OMP_GET_MAX_THREADS


! the old OpenMP implementation for compute_forces_viscoelastic is in utils/unused_routines/:
! older_please_do_not_use_anymore_partial_OpenMP_port/older_not_maintained_compute_forces_viscoelastic_Dev_openmp.f90

  ! gets number of openMP threads
  ! will be determined by environment setting OMP_NUM_THREADS
  !
  ! for example, run executable with:
  ! OMP_NUM_THREADS=4 mpirun -np 2 ./bin/xspecfem3D
  !
  NUM_THREADS = OMP_GET_MAX_THREADS()

  ! output info
  if (myrank == 0) then
    write(IMAIN,*) '  OpenMP:'
    write(IMAIN,*) '    using',NUM_THREADS,' OpenMP threads'
    call flush_IMAIN()
  endif

  ! OpenMP for elastic simulation only supported yet
  if (ELASTIC_SIMULATION) then
    ! safety check
    if (SAVE_MOHO_MESH) then
      print *,'Using the OpenMP feature together with setting SAVE_MOHO_MESH in Par_file is not supported yet!'
      print *,'Please consider adding it as a contribution to this code. For now, re-compile the code without OpenMP feature.'
      call exit_mpi(myrank,'Invalid SAVE_MOHO_MESH setting in OpenMP version')
    endif

    ! below color arrays is probably not needed for cpu version, still left here, just to be sure.
    !
    ! set num_elem_colors array in case no mesh coloring is used
    if (.not. USE_MESH_COLORING_GPU) then
      ! deallocate dummy array
      if (allocated(num_elem_colors_elastic)) deallocate(num_elem_colors_elastic)

      ! loads with corresonding values
      num_colors_outer_elastic = 1
      num_colors_inner_elastic = 1
      allocate(num_elem_colors_elastic(num_colors_outer_elastic + num_colors_inner_elastic),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 2155')
      if (ier /= 0) stop 'error allocating num_elem_colors_elastic array'

      ! sets to all elements in inner/outer phase
      num_elem_colors_elastic(1) = nspec_outer_elastic
      num_elem_colors_elastic(2) = nspec_inner_elastic
    endif

  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_OpenMP
#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_irregular_elements()

! checks indexing of irregular element array

  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,ispec_irreg,ier,num_irreg
  logical, dimension(:), allocatable :: ispec_is_irregular

  ! allocates flags for irregular elements
  allocate(ispec_is_irregular(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_mpi(myrank, 'Error allocating ispec_is_irregular array')
  ! initializes
  ispec_is_irregular(:) = .true.

  ! checks each element
  do ispec = 1,NSPEC_AB
    ! gets irregular number
    ispec_irreg = irregular_element_number(ispec)

    ! sets irregular flag
    if (ispec_irreg == 0) then
      ! regular shaped (cube)
      ispec_is_irregular(ispec) = .false.
    else
      ! irregular shaped (not a perfect cube)
      ispec_is_irregular(ispec) = .true.
    endif

    ! checks index range
    if (ispec_irreg < 0 .or. ispec_irreg > NSPEC_AB) then
      print *,'Error rank',myrank,': invalid irregular index ',ispec_irreg,'in irregular_element_number'
      call exit_mpi(myrank, 'Error invalid index in irregular_element_number array')
    endif
  enddo

  ! number of irregular elements in this partition
  num_irreg = count(ispec_is_irregular(:) .eqv. .true.)

  ! checks
  if (num_irreg /= NSPEC_IRREGULAR) then
    print *,'Error rank',myrank,': invalid number of irregular elements ',num_irreg,'should be ',NSPEC_IRREGULAR
    call exit_mpi(myrank, 'Error invalid number of irregular elements')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  number of regular shaped elements  : ",NSPEC_AB - NSPEC_IRREGULAR
    write(IMAIN,*)"  number of irregular shaped elements: ",NSPEC_IRREGULAR
    call flush_IMAIN()
  endif

  ! frees temporary array
  deallocate(ispec_is_irregular)

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_irregular_elements


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_fused_array()

! prepare fused array for computational kernel
!
! note: fusing separate arrays for xi/eta/gamma into a single one increases efficiency for hardware pre-fetching

  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,ier,ispec_irreg

  ! feature unused so far, will lead to higher memory usage, but might improve performance on some hardware
  logical, parameter :: USE_DERIV_MAPPING_FUSION = .false.

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! allocates fused array
  if (USE_DERIV_MAPPING_FUSION .and. NSPEC_IRREGULAR > 0) then
    allocate(deriv_mapping(10,NGLLX,NGLLY,NGLLZ,NSPEC_IRREGULAR),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'Error allocating array deriv_mapping')

    ! fused array of mapping matrix
    ! (mapping from reference element back to physical element)
    ! d(xi)/d(x)
    do ispec = 1,NSPEC_AB
      ispec_irreg = irregular_element_number(ispec)
      if (ispec_irreg /= 0) then
        DO_LOOP_IJK
          deriv_mapping(1,INDEX_IJK,ispec_irreg) = xixstore(INDEX_IJK,ispec_irreg)
          deriv_mapping(2,INDEX_IJK,ispec_irreg) = xiystore(INDEX_IJK,ispec_irreg)
          deriv_mapping(3,INDEX_IJK,ispec_irreg) = xizstore(INDEX_IJK,ispec_irreg)

          deriv_mapping(4,INDEX_IJK,ispec_irreg) = etaxstore(INDEX_IJK,ispec_irreg)
          deriv_mapping(5,INDEX_IJK,ispec_irreg) = etaystore(INDEX_IJK,ispec_irreg)
          deriv_mapping(6,INDEX_IJK,ispec_irreg) = etazstore(INDEX_IJK,ispec_irreg)

          deriv_mapping(7,INDEX_IJK,ispec_irreg) = gammaxstore(INDEX_IJK,ispec_irreg)
          deriv_mapping(8,INDEX_IJK,ispec_irreg) = gammaystore(INDEX_IJK,ispec_irreg)
          deriv_mapping(9,INDEX_IJK,ispec_irreg) = gammazstore(INDEX_IJK,ispec_irreg)

          deriv_mapping(10,INDEX_IJK,ispec_irreg) = jacobianstore(INDEX_IJK,ispec_irreg)
        ENDDO_LOOP_IJK
      endif
    enddo
  else
    ! dummy
    allocate(deriv_mapping(1,1,1,1,1),stat=ier)
    if (ier /= 0) call exit_mpi(myrank,'Error allocating array deriv_mapping')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  fused array done"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_fused_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_bandwidth_test()

! outputs memory bandwidth performance to know more about the memory bandwidth. this should indicate how fast we can go,
! since our routines are memory-bound...
!
! motivated by: STEAM benchmarks
! http://www.cs.virginia.edu/stream/ref.html

  use constants, only: CUSTOM_REAL,NDIM,IMAIN,myrank
  use specfem_par, only: deltat,FORCE_VECTORIZATION_VAL,NGLOB_AB

  implicit none

  ! local parameters
  integer :: i,ier
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: displ_tmp,veloc_tmp,accel_tmp

  ! timing
  integer :: k
  double precision :: t_s,t_e
  double precision :: t_min,t_max,t_avg
  double precision :: t_min_force,t_max_force,t_avg_force
  double precision :: mem
  double precision, external :: wtime

  ! repeats test
  integer, parameter :: NTIMES = 15

  ! note: we mimick dimension of elastic arrays displ/veloc/accel which will be used later in the code

  ! stores original values
  allocate(displ_tmp(NDIM,NGLOB_AB),veloc_tmp(NDIM,NGLOB_AB),accel_tmp(NDIM,NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays displ_tmp,...'

  ! stream benchmark uses initial values
  displ_tmp(:,:) = 3.0_CUSTOM_REAL
  veloc_tmp(:,:) = 2.0_CUSTOM_REAL
  accel_tmp(:,:) = 1.0_CUSTOM_REAL

  ! synchronizes processes
  call synchronize_all()

  ! timing for Newmark time scheme update

  ! repeats timing exercise for more reliable measure
  t_avg = 0.d0
  t_min = 1.e30
  t_max = 0.d0
  do k = 1,NTIMES
    ! timing
    t_s = wtime()

! openmp solver
!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(NGLOB_AB,displ_tmp,veloc_tmp,accel_tmp,deltat) &
!$OMP PRIVATE(i)
!$OMP DO
    do i = 1,NGLOB_AB
      accel_tmp(1,i) = veloc_tmp(1,i) + deltat * displ_tmp(1,i)
      accel_tmp(2,i) = veloc_tmp(2,i) + deltat * displ_tmp(2,i)
      accel_tmp(3,i) = veloc_tmp(3,i) + deltat * displ_tmp(3,i)
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

    ! synchronizes processes
    call synchronize_all()

    ! timing
    t_e = wtime()

    ! for average time, skips first 5 measurements
    if (k > 5) then
      t_avg = t_avg + (t_e - t_s)
      ! min/max
      if (t_min > (t_e - t_s)) t_min = t_e - t_s
      if (t_max < (t_e - t_s)) t_max = t_e - t_s
    endif
  enddo
  ! takes average
  t_avg = t_avg / dble(NTIMES - 5)

  ! again, for forcing 1D array accessing
  if (FORCE_VECTORIZATION_VAL) then
    ! synchronizes processes
    call synchronize_all()

    ! repeats timing exercise for more reliable measure
    t_avg_force = 0.d0
    t_min_force = 1.e30
    t_max_force = 0.d0
    do k = 1,NTIMES
      ! timing
      t_s = wtime()

! openmp solver
!$OMP PARALLEL &
!$OMP DEFAULT(NONE) &
!$OMP SHARED(NGLOB_AB,displ_tmp,veloc_tmp,accel_tmp,deltat) &
!$OMP PRIVATE(i)
!$OMP DO
      do i = 1,NGLOB_AB * NDIM
        ! following STREAM benchmark: TRIAD a = b + fac * c
        ! we use here the wavefield allocated earlier, since this test could lead to different result
        ! when newly allocated arrays are used

        ! counts:
        ! + 2 FLOP
        !
        ! + 3 * 1 float = 3 * 4 BYTE
        accel_tmp(i,1) = veloc_tmp(i,1) + deltat * displ_tmp(i,1)
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! synchronizes processes
      call synchronize_all()

      ! timing
      t_e = wtime()

      ! for average time, skips first 5 measurements
      if (k > 5) then
        t_avg_force = t_avg_force + (t_e - t_s)
        ! min/max
        if (t_min_force > (t_e - t_s)) t_min_force = t_e - t_s
        if (t_max_force < (t_e - t_s)) t_max_force = t_e - t_s
      endif
    enddo
    ! takes average
    t_avg_force = t_avg_force / dble(NTIMES - 5)
  endif ! FORCE_VECTORIZATION_VAL

! total counts:
!  2 FLOP / 12 BYTE = 0.166 FLOP/BYTE
!
! total memory accesses:
!  NGLOB * NDIM * 3 * 4 BYTE
!
! memory bandwidth
! NVIDIA K20: theoretical 250 GB/s
!      single-precision peak performance: 3.95 TFlop/s -> corner arithmetic intensity = 3950 / 250 ~ 15.8 flop/byte
!      hand-count performance: 0.166 flop/byte ~ 1% of the peak performance ~ 2.5 GB/s
! Intel(R) Xeon(R) Haswell CPU (E5-2680 v3 @ 2.50GHz):
!      http://ark.intel.com/products/81908/Intel-Xeon-Processor-E5-2680-v3-30M-Cache-2_50-GHz
!      Max Memory Bandwidth 68 GB/s

  ! total memory accesses (in MB)
  mem = NGLOB_AB * NDIM * 3.d0 * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0

  ! output bandwidth
  if (myrank == 0) then
    write(IMAIN,*) "  bandwidth test (STREAM TRIAD): "
    write(IMAIN,*) "     memory accesses = ",sngl(mem),'MB'
    write(IMAIN,*) "     timing  min/max = ",sngl(t_min),'s / ',sngl(t_max),'s'
    write(IMAIN,*) "     timing      avg = ",sngl(t_avg),'s'
    if (abs(t_avg) > 0.0) then
      write(IMAIN,*) "     bandwidth       = ",sngl(mem / t_avg / 1024.d0),'GB/s'
    else
      write(IMAIN,*) "     bandwidth       = ",'n/a'
    endif

    if (FORCE_VECTORIZATION_VAL) then
      write(IMAIN,*) "     with force_vectorization:"
      write(IMAIN,*) "     timing  min/max = ",sngl(t_min_force),'s / ',sngl(t_max_force),'s'
      write(IMAIN,*) "     timing      avg = ",sngl(t_avg_force),'s'
      if (abs(t_avg) > 0.0) then
        write(IMAIN,*) "     bandwidth       = ",sngl(mem / t_avg_force / 1024.d0),'GB/s'
      else
        write(IMAIN,*) "     bandwidth       = ",'n/a'
      endif
      ! warning/suggestion to switch, if forcing vectorization slows down
      if (t_avg < 0.9*t_avg_force) then
        write(IMAIN,*) "****"
        write(IMAIN,*) "**** force vectorization slows down performance                           ****"
        write(IMAIN,*) "**** Please consider turning FORCE_VECTORIZATION off, and re-compile code ****"
        write(IMAIN,*) "****"
      endif
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! free temporary arrays
  deallocate(displ_tmp,veloc_tmp,accel_tmp)

  ! todo: would be interesting for GPU bandwidth as well...

  end subroutine prepare_bandwidth_test
