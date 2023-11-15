program combine_signal_files_into_a_single_one


  use hdf5

  implicit none

  integer :: number_of_files ! number of .semp files to be combined and binarized
  integer :: number_of_rows  ! number of time steps
  integer, parameter :: number_of_columns = 2

  integer :: line_start
  integer :: line_end

  double precision, allocatable, dimension(:,:,:) :: big_array
  real, allocatable, dimension(:) :: xx, yy, zz
  integer :: pos

  integer, parameter :: strlen=200

  character(len=strlen) :: filename_to_open, strdump
  character(len=strlen) :: out_h5_name
  character(len=strlen) :: list_file_name

  integer :: ifile, iline

  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dspace_id

  integer(HSIZE_T), dimension(3) :: dims ! Dataset dimensions
  integer :: rank = 3                    ! Dataset rank

  integer :: error  ! Error flag

  ! get commandline arguments
  character(len=strlen) :: argv
  integer :: stat
  call mygetarg(1, argv)
    list_file_name = trim(adjustl(argv))
  call mygetarg(2, argv)
    out_h5_name = trim(adjustl(argv))
  call mygetarg(3, argv)
    call str2int(argv,number_of_files,stat)
  call mygetarg(4, argv)
    call str2int(argv,line_start,stat)
  call mygetarg(5, argv)
    call str2int(argv,line_end,stat)
  call mygetarg(6, argv)
    call str2int(argv,number_of_rows,stat)


  allocate(big_array(number_of_files,line_start:line_end,number_of_columns))
  dims = (/number_of_files,line_end-line_start+1,number_of_columns/)


  open(unit=24,file=list_file_name,status='old',action='read')

  do ifile = 1, number_of_files
    ! read the name of the file to open
    read(24,'(a)') filename_to_open
    write(*,*) filename_to_open

    !write(*,*) trim(adjustl(filename_to_open))
    open(unit=23,file=filename_to_open,status='old',action='read')
      do iline = 1,number_of_rows
        if (iline >= line_start .and. iline <= line_end) then
          read(23,*) big_array(ifile,iline,1),big_array(ifile,iline,2)
        else
          ! read the whole line and ignore it
          read(23,*)
        endif
      enddo
    close(23)

  enddo

  close(24)



  write(*,*) "loading signals end."
  write(*,*) "exporting as a hdf5 file..."


  call h5open_f (error)
  ! Create a new file using default properties.
  call h5fcreate_f(out_h5_name, H5F_ACC_TRUNC_F, file_id, error)
  ! Create the dataspace
  call h5screate_simple_f(rank,dims,dspace_id, error)

  ! Create and write dataset using default properties.
  call h5dcreate_f(file_id, 'time_pressure', H5T_NATIVE_DOUBLE, dspace_id, &
                   dset_id, error, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                   H5P_DEFAULT_F)

  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, big_array, dims, error)


  ! Terminate access to the file.
  call h5dclose_f(dset_id, error)

  call h5sclose_f(dspace_id, error)

  call h5fclose_f(file_id, error)

  call h5close_f(error)

  contains

    subroutine str2int(str,int,stat)
      implicit none
      ! Arguments
      character(len=*),intent(in) :: str
      integer,intent(out)         :: int
      integer,intent(out)         :: stat

      read(str,*,iostat=stat)  int
    end subroutine str2int

    subroutine mygetarg(i, argc)
      implicit none
      integer, intent(in) :: i
      character(len=*), intent(out) :: argc

      call get_command_argument(0, argc)
      if (argc == "") then
        call get_command_argument(i + 1, argc)
      else
        call get_command_argument(i, argc)
      endif
    end subroutine

    integer function myiargc() result(oresult)
      implicit none
      integer :: iargc
      character(len=8) :: argc
      oresult = iargc()
      call get_command_argument(0, argc)
      if (argc == "") then
        oresult = oresult - 1
      endif
    end function

end program combine_signal_files_into_a_single_one

