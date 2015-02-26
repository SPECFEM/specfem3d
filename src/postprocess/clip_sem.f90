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
!   mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL INPUT_FILE OUTPUT_DIR MATERIAL_NAMES
!
!
! COMMAND LINE ARGUMENTS
!   MIN_VAL                - threshold below which array values are clipped
!   MAX_VAL                - threshold above which array values are clipped
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which clipped array are written
!   MATERIAL_NAMES         - one or more material parameter names separated by commas
!
!
! DESCRIPTION
!   MATERIAL_NAMES is comma-delimited list of material names, e.g.'alpha_kernel,beta_kernel'.
!
!   This program works on any  scalar field of appropriate dimension,
!   i.e. (NGLLX,NGLLY,NGLLZ,NSPEC). Its primary use case though is to clip kernels.
!
!   This is an embarassingly-parallel program.


program clip_sem

  use postprocess_par,only: MAX_STRING_LEN,IIN,IOUT, &
    myrank,sizeprocs,NGLOB,NSPEC,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL

  use shared_parameters

  implicit none

  character(len=MAX_STRING_LEN) :: input_dir,output_dir,filename
  character(len=MAX_STRING_LEN) :: arg(5)
  integer :: i,ier

  ! machinery for tokenizing comma-delimited list of material names
  character(len=255) :: strtok
  character(len=1) :: delimiter
  integer :: imat,nmat

  ! given that there are a maximum of 21 elastic moduli plus density,
  ! it is unlikely there will ever be a need for more than 22 names
  integer,parameter :: MAX_MATERIAL_NAMES = 22
  character(len=MAX_STRING_LEN) :: material_names(MAX_MATERIAL_NAMES)
  character(len=MAX_STRING_LEN) :: material_names_comma_delimited
  character(len=MAX_STRING_LEN) :: mat

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sem_array

  double precision :: min_val, max_val

  logical :: BROADCAST_AFTER_READ

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  do i = 1, 5
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 1 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL INPUT_FILE OUTPUT_DIR MATERIAL_NAMES'
      stop ' Reenter command line options'
      endif
    endif
  enddo

  ! gets arguments
  read(arg(1),*) min_val
  read(arg(2),*) max_val
  read(arg(3),'(a)') input_dir
  read(arg(4),'(a)') output_dir
  read(arg(5),'(a)') material_names_comma_delimited

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
    write(*,*) 'Running XCLIP_SEM'
  endif
  call synchronize_all()

  ! needs local_path for mesh files
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print*,''
      print*,'Error: run xclip_sem with the same number of MPI processes '
      print*,'       as specified in Par_file by NPROC when slices were created'
      print*,''
      print*,'for example: mpirun -np ',NPROC,' ./xlcip_sem ...'
      print*,''
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()

  ! opens external mesh file
  write(filename,'(a,i6.6,a)') trim(LOCAL_PATH)//'/proc',myrank,'_'//'external_mesh.bin'
  open(unit=27,file=trim(filename),&
          status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error: could not open database '
    print*,'path: ',trim(filename)
    stop 'Error reading external mesh file'
  endif

  read(27) NSPEC
  read(27) NGLOB
  close(27)


  allocate(sem_array(NGLLX,NGLLY,NGLLZ,NSPEC))

  do imat=1,nmat

      mat = trim(material_names(imat))
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(mat)//'.bin'

      ! read array
      open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(filename)
        stop 'File not found'
      endif
      read(IIN) sem_array
      close(IIN)

      call clip_sem_array(sem_array,min_val,max_val)

      ! write clipped array
      mat = trim(material_names(imat))//'_clip'
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

  ! stop all the processes, and exit
  call finalize_mpi()

end program clip_sem

!
!-------------------------------------------------------------------------------------------------
!

subroutine clip_sem_array(sem_array,min_val,max_val)

  use postprocess_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: sem_array
  double precision :: min_val, max_val
  integer :: i,j,k,ispec

  ! apply threshold
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

end subroutine clip_sem_array


!
!-------------------------------------------------------------------------------------------------
!

! The following utility function was made freely available by the Fortran Wiki:
! http://fortranwiki.org/fortran/show/strtok
!
character(len=255) function strtok (source_string, delimiters)

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in DELIMITERS.
!
!             then, if the returned value is not equal to char(0), keep calling until it is
!             with SOURCE_STRING set to char(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning char(0).
!
!     Input:  source_string =   Source string to tokenize.
!             delimiters    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures
!     input string length limited to set size
!     function returns fixed 255 character length
!     length of returned string not given

!     PARAMETERS:
      character(len=*),intent(in)  :: source_string
      character(len=*),intent(in)  :: delimiters

!     SAVED VALUES:
      character(len=255),save :: saved_string
      integer,save :: isaved_start  ! points to beginning of unprocessed data
      integer,save :: isource_len   ! length of original input string

!     LOCAL VALUES:
      integer :: ibegin        ! beginning of token to return
      integer :: ifinish       ! end of token to return

      ! initialize stored copy of input string and pointer into input string on first call
      if (source_string(1:1) /= char(0)) then
          isaved_start = 1                 ! beginning of unprocessed data
          saved_string = source_string     ! save input string from first call in series
          isource_len = LEN(saved_string)  ! length of input string from first call
      endif

      ibegin = isaved_start

      do
         if ( (ibegin <= isource_len) .AND. (index(delimiters,saved_string(ibegin:ibegin)) /= 0)) then
             ibegin = ibegin + 1
         else
             exit
         endif
      enddo

      if (ibegin > isource_len) then
          strtok = char(0)
          RETURN
      endif

      ifinish = ibegin

      do
         if ((ifinish <= isource_len) .AND.  (index(delimiters,saved_string(ifinish:ifinish)) == 0)) then
             ifinish = ifinish + 1
         else
             exit
         endif
      enddo

      !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
      strtok = saved_string(ibegin:ifinish-1)
      isaved_start = ifinish

end function strtok




