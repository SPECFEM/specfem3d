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

! XCOMBINE_SEM
!
! USAGE
!   mpirun -np NPROC bin/xcombine_sem INPUT_FILE OUTPUT_DIR KERNELL_NAMES
!
!
! COMMAND LINE ARGUMENTS
!   INPUT_FILE             - text file containing list of kernel directories
!   OUTPUT_PATH            - directory to which summed kernels are written
!   KERNEL_NAMES           - one or more material parameter names separated by
!   commas
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, sums kernels from directories specified in
!   INPUT_FILE.
!   Writes the resulting sums to OUTPUT_DIR.
!
!   INPUT_FILE is a text file containing a list of absolute or relative paths to
!   kernel direcotires, one directoy per line.
!
!   KERNEL_NAMES is comma-delimited list of kernel names, 
!   e.g.'alpha_kernel,beta_kernel,rho_kernel'.



program combine_sem

  use postprocess_par,only: MAX_STRING_LEN,MAX_KERNEL_PATHS,IIN, &
    myrank,sizeprocs,NGLOB,NSPEC

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN) :: paths_list(MAX_KERNEL_PATHS), material_names(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: sline,prname_lp,output_dir,input_file,material_names_comma_delimited
  character(len=MAX_STRING_LEN) :: arg(3)
  character(len=255) :: strtok
  character(len=1) :: delimiter
  integer :: nker,nmat
  integer :: i,ier,imat

  logical :: BROADCAST_AFTER_READ

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  do i = 1, 3
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 1 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xcombine_sem INPUT_FILE OUTPUT_DIR MATERIAL_NAMES'
      stop ' Reenter command line options'
      endif
    endif
  enddo

  ! gets arguments
  read(arg(1),'(a)') input_file
  read(arg(2),'(a)') output_dir
  read(arg(3),'(a)') material_names_comma_delimited

  ! tokenize comma-delimited list of material names
  delimiter = ','
  imat = 1
  material_names(imat) = trim(strtok(material_names_comma_delimited, delimiter))
  do while (material_names(imat) /= char(0))
     imat = imat + 1
     material_names(imat) = trim(strtok(char(0), delimiter))
  enddo
  nmat = imat-1

  if (myrank==0) then
    write(*,*) 'Running COMBINE_SEM:'
    write(*,*)
  endif
  call synchronize_all()

  ! reads in event list
  nker=0
  open(unit = IIN, file = trim(input_file), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(input_file), myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     nker = nker+1
     if (nker > MAX_KERNEL_PATHS) stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     paths_list(nker) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print*,''
      print*,'Error: run xcombine_sem with the same number of MPI processes '
      print*,'       as specified in Par_file by NPROC when slices were created'
      print*,''
      print*,'for example: mpirun -np ',NPROC,' ./xcombine_sem ...'
      print*,''
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()

  ! reads mesh file
  !
  ! needs to get array dimensions

  ! opens external mesh file
  write(prname_lp,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=27,file=trim(prname_lp),&
          status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error: could not open database '
    print*,'path: ',trim(prname_lp)
    stop 'Error reading external mesh file'
  endif

  ! gets number of elements and global points for this partition
  read(27) NSPEC
  read(27) NGLOB

  close(27)

  ! user output
  if (myrank == 0) then
    print*,'summing arrays in:'
    print*,paths_list(1:nker)
    print*
  endif

  ! synchronizes
  call synchronize_all()

  do imat=1,nmat
      call combine_sem_array(material_names(imat),paths_list,output_dir,nker)
  enddo


  if (myrank==0) write(*,*) 'done writing all arrays, see directory', output_dir

  ! stop all the processes, and exit
  call finalize_mpi()

end program combine_sem

!
!-------------------------------------------------------------------------------------------------
!

subroutine combine_sem_array(material_name,paths_list,output_dir,nker)

  use postprocess_par

  implicit none

  character(len=MAX_STRING_LEN) :: material_name,paths_list(MAX_KERNEL_PATHS),output_dir
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN*2) :: filename
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: array,sum_arrays
  double precision :: norm,norm_sum
  integer :: iker,ier

  ! initializes arrays
  allocate(array(NGLLX,NGLLY,NGLLZ,NSPEC), &
           sum_arrays(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating array'

 ! loop over array paths
  sum_arrays = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if (myrank==0) then
      write(*,*) 'reading in array for: ',trim(material_name)
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! read array
    array = 0._CUSTOM_REAL
    write(filename,'(a,i6.6,a)') trim(paths_list(iker)) //'/proc',myrank,'_'//trim(material_name)//'.bin'

    open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  array not found: ',trim(filename)
      stop 'Error array file not found'
    endif
    read(IIN) array
    close(IIN)

    norm = sum( array * array )
    call sum_all_dp(norm, norm_sum)
    if (myrank == 0) then
      print*,'  norm array: ',sqrt(norm_sum)
      print*
    endif

    ! add to sum
    sum_arrays = sum_arrays + array

  enddo

  ! stores summed arrays
  if (myrank==0) write(*,*) 'writing out summed array for: ',trim(material_name)

  write(filename,'(a,i6.6,a)') trim(output_dir)//'/'//'proc',myrank,'_'//trim(material_name)//'.bin'

  open(IOUT,file=trim(filename),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error array not written:',trim(filename)
    stop 'Error array write'
  endif
  write(IOUT) sum_arrays
  close(IOUT)

  if (myrank==0) write(*,*)

  ! frees memory
  deallocate(array,sum_arrays)

end subroutine combine_sem_array


