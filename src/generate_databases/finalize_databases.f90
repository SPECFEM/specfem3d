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
!

  subroutine finalize_databases()

! checks user input parameters

  use generate_databases_par
  use create_regions_mesh_ext_par, only: nspec_irregular

  implicit none

  ! local parameters
  integer :: nspec_total
  ! this can overflow if more than 2 Gigapoints in the whole mesh, thus replaced with double precision version
  integer(kind=8) :: nglob_total
  double precision :: nglob_l,nglob_total_db

  ! timing
  double precision, external :: wtime

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

    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements:'
    write(IMAIN,*) '-----------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in mesh slice 0: ',NSPEC_AB
    write(IMAIN,*) 'total number of   regular elements in mesh slice 0: ',NSPEC_AB - nspec_irregular
    write(IMAIN,*) 'total number of irregular elements in mesh slice 0: ',nspec_irregular
    write(IMAIN,*) 'total number of points in mesh slice 0: ',NGLOB_AB
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in entire mesh: ',nspec_total
    write(IMAIN,*) 'approximate total number of points in entire mesh (with duplicates on MPI edges): ',nglob_total
    write(IMAIN,*) 'approximate total number of DOFs in entire mesh (with duplicates on MPI edges): ',nglob_total*NDIM
    write(IMAIN,*)
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

