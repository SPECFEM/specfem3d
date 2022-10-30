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

  program project_and_combine_vol_data_on_regular_grid

! combines the database files on several slices, project it on a regular grid
! and saves it in a binary file that can be read with Paraview in RAW format
!
! works for external, unregular meshes

  use constants
  use shared_parameters
  use specfem_par, only: xigll, yigll, zigll, wxgll, wygll, wzgll

  use projection_on_FD_grid

  implicit none

  integer, parameter :: NARGS = 3

  ! data must be of dimension: (NGLLX,NGLLY,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: data
  integer :: i, ier
  character(len=MAX_STRING_LEN) :: arg(9), indir, outdir
  character(len=MAX_STRING_LEN) :: prname, prname_lp, data_filename
  character(len=MAX_STRING_LEN*2) :: local_data_file
  logical :: BROADCAST_AFTER_READ
  integer :: sizeprocs, NSPEC_IRREGULAR

  type(profd)  :: projection_fd
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::      model_on_FD_grid

  ! MPI initialization
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    print *
    print *,'Projecting volumetric data on a regular grid'
    print *
  endif

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
      print *,'USAGE:  mpirun -np NPROC bin/xproject_and_combine_vol_data_on_regular_grid data_filename input_dir output_dir'
      stop 'Please check command line arguments'
    endif
  endif

  ! reads in arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i))
  enddo

  data_filename = arg(1)
  indir = arg(2)
  outdir = arg(3)

  ! Get dimensions of current model, stored in proc******_external_mesh.bin
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'
  open(unit=27,file=prname_lp(1:len_trim(prname_lp))//'external_mesh.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  read(27) NSPEC_AB
  read(27) NGLOB_AB

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1102')
  if (ier /= 0) stop 'error allocating array ibool'
  allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1103')
  if (ier /= 0) stop 'error allocating array xstore etc.'

  read(27) NSPEC_IRREGULAR
  read(27) ibool
  read(27) xstore
  read(27) ystore
  read(27) zstore
  close(27)

  ! Get regular grid properties, and precomputes all interpolation coefficients
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  call compute_interpolation_coeff_FD_SEM(projection_fd, myrank)

  allocate(model_on_FD_grid(projection_fd%nx, projection_fd%ny, projection_fd%nz),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1104')

  if (myrank == 0) then
    print *, 'Grid size is : ',projection_fd%nx, projection_fd%ny, projection_fd%nz
  endif

  ! Get data to project
  allocate(data(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1106')
  if (ier /= 0) stop 'error allocating double precision data array'

  ! data file
  write(prname,'(a,i6.6,a)') trim(indir)//'proc',myrank,'_'
  local_data_file = trim(prname) // trim(data_filename) // '.bin'
  open(unit = 28,file = trim(local_data_file),status='old',action='read',form ='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening ',trim(local_data_file)
    stop
  endif

  ! Read (either SP or DP) floating point numbers.
  read(28) data
  close(28)

  ! put data from SEM mesh to a regular grid
  call Project_model_SEM2FD_grid(data, model_on_FD_grid, projection_fd, myrank)

  ! Write output on a Fortran binary file
  if (myrank == 0) then
    open(unit = 28,file = trim(outdir)//trim(data_filename) // '_projected.bin',status='unknown', &
          action='write',form ='unformatted',iostat=ier)
    write(28) model_on_FD_grid
    close(28)
  endif

  call finalize_mpi()

  print *, 'Done writing ' // trim(outdir)//trim(data_filename) // '_projected.bin'

  end program project_and_combine_vol_data_on_regular_grid
