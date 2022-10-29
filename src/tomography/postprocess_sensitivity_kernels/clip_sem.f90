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

! XCLIP_SEM
!
! USAGE
!   mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL KERNEL_NAMES INPUT_FILE OUTPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   MIN_VAL                - threshold below which array values are clipped
!   MAX_VAL                - threshold above which array values are clipped
!   KERNEL_NAMES           - one or more kernel names separated by commas
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which clipped array are written
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, reads kernels from INPUT_DIR, applies
!   thresholds, and writes the resulting clipped kernels to OUTPUT_DIR.
!
!   KERNEL_NAMES is comma-delimited list of kernel names,
!   e.g. 'alphav_kernel,alphah_kernel'
!
!   Files written to OUTPUT_DIR have the suffix 'clip' appended,
!   e.g. proc***alpha_kernel.bin becomes proc***alpha_kernel_clip.bin
!
!   This program's primary use case is to clip kernels. It can be used though on
!   any scalar field of dimension (NGLLX,NGLLY,NGLLZ,NSPEC).
!
!   This is a parrallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.


program clip_sem

  use postprocess_par, only: MAX_STRING_LEN,IIN,IOUT, &
    myrank,sizeprocs,NGLOB,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,MAX_KERNEL_NAMES

  use shared_parameters

  implicit none

  integer, parameter :: NARGS = 5

  character(len=MAX_STRING_LEN) :: input_dir,output_dir,filename
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  integer :: ier, iker,nker,i,j,k,ispec

  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_names
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited, kernel_name

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sem_array

  double precision :: min_val, max_val

  logical :: BROADCAST_AFTER_READ


  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then
    write(*,*) 'Running XCLIP_SEM'
    write(*,*)
  endif
  call synchronize_all()

  ! check command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
    endif
    stop 'Error wrong number of arguments'
    call synchronize_all()
  endif

  ! allocates array
  allocate(kernel_names(MAX_KERNEL_NAMES),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_names array'
  kernel_names(:) = ''

  ! parse command line arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo
  call synchronize_all()

  read(arg(1),*) min_val
  read(arg(2),*) max_val
  read(arg(3),'(a)') kernel_names_comma_delimited
  read(arg(4),'(a)') input_dir
  read(arg(5),'(a)') output_dir

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  ! read simulation parameters
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

  ! checks number of MPI processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *
      print *,'Expected number of MPI processes: ', NPROC
      print *,'Actual number of MPI processes: ', sizeprocs
      print *
    endif
    call synchronize_all()
    stop 'Error wrong number of MPI processes'
  endif
  call synchronize_all()

  ! read mesh dimensions
  write(filename,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=IIN,file=trim(filename), &
          status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open external mesh file '
    print *,'path: ',trim(filename)
    stop 'Error reading external mesh file'
  endif

  read(IIN) NSPEC
  read(IIN) NGLOB
  close(IIN)
  call synchronize_all()

  allocate(sem_array(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 977')

  ! clip kernels
  do iker=1,nker

      kernel_name = trim(kernel_names(iker))
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_name)//'.bin'

      ! read array
      open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(filename)
        stop 'File not found'
      endif
      read(IIN) sem_array
      close(IIN)

     if (myrank == 0) then
        write(*,*) 'clipping array: ',trim(kernel_names(iker))
        write(*,*) '  min/max values = ',min_val,max_val
     endif

     ! apply thresholds
      do ispec=1,NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              if (sem_array(i,j,k,ispec) < min_val) sem_array(i,j,k,ispec) = min_val
              if (sem_array(i,j,k,ispec) > max_val) sem_array(i,j,k,ispec) = max_val
            enddo
          enddo
        enddo
      enddo

      ! write clipped array
      kernel_name = trim(kernel_names(iker))//'_clip'
      write(filename,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'.bin'

      open(IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  error opening file: ',trim(filename)
        stop 'Error opening file'
      endif
      write(IOUT) sem_array
      close(IOUT)

  enddo


  if (myrank == 0) write(*,*) 'done clipping all arrays, see directory: ', trim(output_dir)
  deallocate(sem_array)
  call finalize_mpi()

end program clip_sem

