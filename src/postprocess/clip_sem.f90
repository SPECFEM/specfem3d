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

! XCLIP_SEM
!
! USAGE
!   mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL INPUT_FILE OUTPUT_DIR KERNEL_NAMES
!
!
! COMMAND LINE ARGUMENTS
!   MIN_VAL                - threshold below which array values are clipped
!   MAX_VAL                - threshold above which array values are clipped
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which clipped array are written
!   KERNEL_NAMES         - one or more material parameter names separated by commas
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, reads kernels from INPUT_DIR, applies 
!   thresholds, and writes the resulting clipped kernels to OUTPUT_DIR.
!
!   KERNEL_NAMES is comma-delimited list of material names, 
!   e.g. 'alpha_kernel,beta_kernel,rho_kernel'
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

  use postprocess_par,only: MAX_STRING_LEN,IIN,IOUT, &
    myrank,sizeprocs,NGLOB,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL, &
    MAX_KERNEL_NAMES

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN) :: input_dir,output_dir,filename
  character(len=MAX_STRING_LEN) :: arg(5)
  integer :: ier

  ! machinery for tokenizing comma-delimited list of material names
  character(len=255) :: strtok
  character(len=1) :: delimiter
  integer :: imat,nmat,i,j,k,ispec

  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: mat

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sem_array

  double precision :: min_val, max_val

  logical :: BROADCAST_AFTER_READ

  ! ============ program starts here =====================

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! parse command line arguments
  do i = 1, 5
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 1 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL INPUT_FILE OUTPUT_DIR KERNEL_NAMES'
      print *, ''
      stop 'Please check command line arguments'
      endif
    endif
  enddo

  read(arg(1),*) min_val
  read(arg(2),*) max_val
  read(arg(3),'(a)') input_dir
  read(arg(4),'(a)') output_dir
  read(arg(5),'(a)') kernel_names_comma_delimited

  ! parse kernel names
  delimiter = ','
  imat = 1
  kernel_names(imat) = trim(strtok(kernel_names_comma_delimited, delimiter))
  do while (kernel_names(imat) /= char(0))
     imat = imat + 1
     kernel_names(imat) = trim(strtok(char(0), delimiter))
  enddo
  nmat = imat-1

  ! print status update
  if (myrank==0) then
    write(*,*) 'Running XCLIP_SEM'
  endif
  call synchronize_all()

  ! read simulation parameters
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks number of MPI processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print*,''
      print*,'Expected number of MPI processes: ', NPROC
      print*,'Actual number of MPI processes: ', sizeprocs
      print*,''
    endif
    call synchronize_all()
    stop 'Error wrong number of MPI processes'
  endif
  call synchronize_all()

  ! read mesh dimensions
  write(filename,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=27,file=trim(filename),&
          status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error: could not open external mesh file '
    print*,'path: ',trim(filename)
    stop 'Error reading external mesh file'
  endif

  read(27) NSPEC
  read(27) NGLOB
  close(27)
  call synchronize_all()

  allocate(sem_array(NGLLX,NGLLY,NGLLZ,NSPEC))

  ! clip kernels
  do imat=1,nmat

      mat = trim(kernel_names(imat))
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(mat)//'.bin'

      ! read array
      open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(filename)
        stop 'File not found'
      endif
      read(IIN) sem_array
      close(IIN)

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
      if (myrank==0) write(*,*) 'writing array: ',trim(kernel_names(imat))
      mat = trim(kernel_names(imat))//'_clip'
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(mat)//'.bin'

      open(IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  error opening file: ',trim(filename)
        stop 'Error opening file'
      endif
      write(IOUT) sem_array
      close(IOUT)

  enddo


  if (myrank==0) write(*,*) 'done clipping all arrays, see directory', trim(output_dir)
  deallocate(sem_array)
  call finalize_mpi()

end program clip_sem

