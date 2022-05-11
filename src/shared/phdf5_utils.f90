!
! phdf5 file I/O routines
!


module phdf5_utils ! class-like module
    use constants, only: CUSTOM_REAL, MAX_LENGTH_NETWORK_NAME, MAX_LENGTH_STATION_NAME
    use hdf5
    implicit none

    private
    public :: h5io, h5_init, h5_destructor, &
        h5_check_collective, &
        h5_create_file, h5_open_file, h5_close_file, &
        h5_create_group, h5_open_group, h5_close_group, &
        h5_create_subgroup, h5_open_subgroup, h5_open_or_create_group, h5_close_subgroup, &
        h5_open_dataset, h5_open_dataset2, h5_close_dataset, &
        h5_write_dataset_1d_i, h5_write_dataset_1d_d, h5_write_dataset_1d_c,&
        h5_write_dataset_1d_i_no_group, h5_write_dataset_1d_d_no_group, h5_write_dataset_1d_c_no_group,&
        h5_write_dataset_2d_d, h5_write_dataset_2d_r, h5_write_dataset_2d_i, h5_write_dataset_2d_c, &
        h5_write_dataset_2d_r_no_group, &
        h5_write_dataset_4d_r, &
        h5_add_attribute_i, &
        h5_set_mpi_info, &
        h5_set_buffer_size, &
        h5_set_sieve_buffer_size, &
        h5_create_file_p, h5_create_file_p_collect, &
        h5_open_file_p, h5_open_file_p_collect, &
        h5_close_file_p, &
        h5_create_file_prop_list, h5_close_prop_list,h5_close_prop_list_nocheck, &
        h5_open_group_prop_list, h5_create_group_prop_list, h5_close_group_prop_list,&
        h5_create_dataset_prop_list, &
        h5_create_dataset_gen, h5_create_dataset_gen_in_group, &
        h5_create_group_p, h5_open_group_p, &
        h5_read_dataset_p_scalar_i, h5_read_dataset_p_scalar_r,&
        h5_read_dataset_p_1d_i, h5_read_dataset_p_1d_r, h5_read_dataset_p_1d_l,&
        h5_read_dataset_p_2d_i, h5_read_dataset_p_2d_r, h5_read_dataset_p_2d_c, h5_read_dataset_p_2d_d,&
        h5_read_dataset_p_3d_i, h5_read_dataset_p_3d_r, &
        h5_read_dataset_p_4d_i, h5_read_dataset_p_4d_r, &
        h5_read_dataset_p_5d_r, &
        h5_read_attribute_p, &
        h5_read_dataset_scalar_i_collect_hyperslab, &
        h5_read_dataset_scalar_r_collect_hyperslab, &
        h5_read_dataset_1d_l_collect_hyperslab, &
        h5_read_dataset_1d_i_collect_hyperslab, &
        h5_read_dataset_1d_r_collect_hyperslab, &
        h5_read_dataset_1d_r_collect_hyperslab_in_group, &
        h5_read_dataset_2d_i_collect_hyperslab, &
        h5_read_dataset_2d_i_collect_hyperslab_in_group, &
        h5_read_dataset_2d_r_collect_hyperslab, &
        h5_read_dataset_2d_d_collect_hyperslab, &
        h5_read_dataset_3d_i_collect_hyperslab, &
        h5_read_dataset_3d_r_collect_hyperslab, &
        h5_read_dataset_4d_i_collect_hyperslab, &
        h5_read_dataset_4d_r_collect_hyperslab, &
        h5_read_dataset_5d_r_collect_hyperslab, &
        h5_set_group_name, &
        h5_write_dataset_p_1d_i, h5_write_dataset_p_1d_ia, h5_write_dataset_p_1d_r, h5_write_dataset_p_1d_l, &
        h5_write_dataset_p_2d_i, h5_write_dataset_p_2d_ia, h5_write_dataset_p_2d_r, h5_write_dataset_p_2d_d,&
        h5_write_dataset_p_3d_i, h5_write_dataset_p_3d_r, &
        h5_write_dataset_p_4d_i, h5_write_dataset_p_4d_r, &
        h5_write_dataset_p_5d_r, &
        bool_array2integer, int_array2bool, &
        h5_gather_dsetsize, create_dataset_collect, &
        h5_write_dataset_1d_to_2d_r_collect_hyperslab, h5_write_dataset_2d_to_3d_r_collect_hyperslab, &
        h5_write_dataset_1d_l_collect_hyperslab, &
        h5_write_dataset_1d_i_collect_hyperslab, &
        h5_write_dataset_1d_r_collect_hyperslab, &
        h5_write_dataset_1d_r_collect_hyperslab_in_group, &
        h5_write_dataset_2d_i_collect_hyperslab, &
        h5_write_dataset_2d_i_collect_hyperslab_in_group, &
        h5_write_dataset_2d_r_collect_hyperslab, &
        h5_write_dataset_2d_d_collect_hyperslab, &
        h5_write_dataset_3d_i_collect_hyperslab, &
        h5_write_dataset_3d_r_collect_hyperslab, &
        h5_write_dataset_4d_i_collect_hyperslab, &
        h5_write_dataset_4d_r_collect_hyperslab, &
        h5_write_dataset_5d_r_collect_hyperslab, &
        write_attenuation_file_in_h5, read_attenuation_file_in_h5, &
        write_checkmesh_data_h5, write_checkmesh_xdmf_h5



    ! class-wide private variables
    character(len=256) :: file_path ! file path

    character(len=256) :: store_group_name !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! groupname for debug
    integer(HID_T) :: file_id, group_id, parent_group_id, dataset_id!, attribute_id
    integer(HID_T) :: mem_dspace_id, file_dspace_id ! for collective IO
    integer :: error

    ! parallel process
    integer(HID_T) :: plist_id, fplist_id, gplist_id

    ! string array
    integer(SIZE_T)                                    :: str_len = 16
    integer(HID_T)                                     :: str_type


    ! mpi info
    integer :: this_info, this_comm, this_rank, total_proc

    ! dummy array for generate 0 length dataset
    integer, dimension(1) :: dummy_1d_array = (/0/)

    type(c_ptr) :: f_ptr

    ! default constructor
    type h5io
    end type h5io

contains
    !
    ! serial write routines
    !
    subroutine h5_init(this,fpath_in)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: fpath_in

        ! set the target file path
        file_path = fpath_in
        call h5open_f(error)

        ! prepare string array type
        if (MAX_LENGTH_STATION_NAME >= MAX_LENGTH_NETWORK_NAME) then
            str_len = MAX_LENGTH_STATION_NAME
        else
            str_len = MAX_LENGTH_NETWORK_NAME
        endif

        call h5tcopy_f(H5T_FORTRAN_S1, str_type, error)
        call h5tset_size_f(str_type, str_len, error)
    end subroutine h5_init


    subroutine h5_destructor(this)
        type(h5io), intent(in) :: this
        call h5close_f(error)
    end subroutine h5_destructor


    subroutine h5_check_collective(this,dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer :: io_mode
        call h5pget_mpio_actual_io_mode_f(plist_id,io_mode,error)
        !if (io_mode == H5D_MPIO_NO_COLLECTIVE_F) print*, &
        !    "collective read/write not possible for dataset: ", dataset_name

    end subroutine h5_check_collective


    subroutine h5_create_file(this)
        type(h5io), intent(in) :: this

        call h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, error)
        if (error /= 0) then
            print *,'Error file open:', file_path
            print *
            print *,'check if path exists:', file_path
            stop 'Error file open h5'
        endif
    end subroutine h5_create_file


    subroutine h5_open_file(this)
        type(h5io), intent(in) :: this
        call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error)
        if (error /= 0) then
            print *, 'Error file open', file_path
            stop 'Error file open h5'
        endif
    end subroutine h5_open_file


    subroutine h5_close_file(this)
        type(h5io), intent(in) :: this
        call h5fclose_f(file_id, error)
        if (error /= 0) write(*,*) 'error while closing a file.'
    end subroutine h5_close_file

    subroutine h5_close_file_p(this)
        type(h5io), intent(in) :: this
        call h5_close_file_prop_list(this)
        call h5fclose_f(file_id, error)
        if (error /= 0) write(*,*) 'error while closing a file.'
    end subroutine h5_close_file_p


    subroutine h5_create_group(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name

        call h5gcreate_f(file_id, trim(group_name), group_id, error)
        if (error /= 0) write(*,*) 'error while creating a group, ', group_name
        call h5gclose_f(group_id, error)
        if (error /= 0) write(*,*) 'error while closing a group, ' , group_name
    end subroutine h5_create_group


    subroutine h5_create_subgroup(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        integer(HID_T) :: temp_group_id

        call h5gcreate_f(group_id, trim(group_name), temp_group_id, error)
        if (error /= 0) write(*,*) 'error while creating a subgroup, ', group_name
        call h5gclose_f(temp_group_id, error)
        if (error /= 0) write(*,*) 'error while closing a subgroup, ' , group_name
    end subroutine h5_create_subgroup


    subroutine h5_open_group(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        call h5gopen_f(file_id, trim(group_name), group_id, error) ! group open
        if (error /= 0) write(*,*) 'hdf5 open group failed for, ', group_name

        store_group_name = group_name
    end subroutine h5_open_group


    subroutine h5_open_subgroup(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        integer(HID_T) :: temp_group_id

        call h5gopen_f(group_id, trim(group_name), temp_group_id, error) ! group open
        if (error /= 0) write(*,*) 'hdf5 open subgroup failed for, ', group_name

        ! put the group id of the first level becomes parent group
        ! only while the group at second level exists
        parent_group_id = group_id
        group_id        = temp_group_id
    end subroutine h5_open_subgroup

    subroutine h5_open_or_create_group(this, group_name)
        ! Check if group with the given name exists. Create it if it doesn't,
        ! open it if it does.
        implicit none
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        ! Variable for checking if a group exists or not
        logical :: group_exists

        call h5lexists_f(file_id, group_name, group_exists, error)
        if (group_exists) then
            call h5gopen_f(file_id, group_name, group_id, error)
        else
            call h5gcreate_f(file_id, group_name, group_id, error)
        end if
    end subroutine h5_open_or_create_group


    subroutine h5_close_group(this)
        type(h5io), intent(in) :: this
        call h5gclose_f(group_id, error) ! group open
        if (error /= 0) write(*,*) 'hdf5 close group failed'
    end subroutine h5_close_group


    subroutine h5_close_subgroup(this)
        type(h5io), intent(in) :: this
        call h5gclose_f(group_id, error) ! group open
        if (error /= 0) write(*,*) 'hdf5 close group failed'

        group_id = parent_group_id
    end subroutine h5_close_subgroup



    subroutine h5_set_group_name(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        store_group_name = group_name
    end subroutine h5_set_group_name


    ! open dataset. group need to be opened before.
    subroutine h5_open_dataset(this, dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name

        call h5dopen_f(group_id, trim(dataset_name), dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 open dataset failed for, ', dataset_name
    end subroutine h5_open_dataset


    ! open dataset without open group
    subroutine h5_open_dataset2(this, dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name

        call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 open dataset failed for, ', dataset_name
    end subroutine h5_open_dataset2


    subroutine h5_close_dataset(this)
        type(h5io), intent(in) :: this
        call h5dclose_f(dataset_id, error) ! group open
        if (error /= 0) write(*,*) 'hdf5 close dataset failed'
    end subroutine h5_close_dataset


    ! dataset writer for 1d integer array
    subroutine h5_write_dataset_1d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 1
        integer(HSIZE_T), dimension(1)    :: dim
        dim = size(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_1d_i


    subroutine h5_write_dataset_1d_i_no_group(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 1
        integer(HSIZE_T), dimension(1)    :: dim
        dim = size(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name

        call h5_close_dataset(this)
    end subroutine h5_write_dataset_1d_i_no_group



    ! dataset writer for 1d custom real array
    subroutine h5_write_dataset_1d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 1
        integer(HSIZE_T), dimension(1)    :: dim
        dim = size(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        if (CUSTOM_REAL == 4) then
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        else
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        endif
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_1d_d


    ! dataset writer for 1d custom real array without grouping
    subroutine h5_write_dataset_1d_d_no_group(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 1
        integer(HSIZE_T), dimension(1)    :: dim
        dim = size(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        if (CUSTOM_REAL == 4) then
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        else
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        endif
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_1d_d_no_group


    subroutine h5_write_dataset_1d_c(this, dataset_name, data_in)
        type(h5io), intent(in)                            :: this
        character(len=*), intent(in)                      :: dataset_name
        character(len=*), dimension(:), intent(in)        :: data_in
        integer(HID_T)                                    :: dspace_id ! dataspace id is local.
        character(len=str_len), dimension(:), allocatable, target ::  data
        integer                                           :: rank = 1
        integer(HSIZE_T), dimension(1)                    :: dim
        integer                                           :: i

        dim = shape(data_in)
        allocate(data(dim(1)),stat=error)
        ! fill blancs after each string to have the same length of strings
        do i = 1, dim(1)
            data(i) = data_in(i)//repeat(' ',str_len-len(data_in(i)))
        enddo

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(group_id, trim(dataset_name), str_type, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, str_type, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name

        deallocate(data, stat=error)

    end subroutine h5_write_dataset_1d_c


    subroutine h5_write_dataset_1d_c_no_group(this, dataset_name, data_in)
        type(h5io), intent(in)                            :: this
        character(len=*), intent(in)                      :: dataset_name
        character(len=*), dimension(:), intent(in)        :: data_in
        integer(HID_T)                                    :: dspace_id ! dataspace id is local.
        character(len=str_len), dimension(:), allocatable, target ::  data
        integer                                           :: rank = 1
        integer(HSIZE_T), dimension(1)                    :: dim
        integer                                           :: i

        dim = shape(data_in)
        allocate(data(dim(1)),stat=error)
        ! fill blancs after each string to have the same length of strings
        do i = 1, dim(1)
            data(i) = data_in(i)//repeat(' ',str_len-len(data_in(i)))
        enddo

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(file_id, trim(dataset_name), str_type, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, str_type, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name

        deallocate(data, stat=error)

    end subroutine h5_write_dataset_1d_c_no_group


    ! dataset writer for 2d integer array
    subroutine h5_write_dataset_2d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 2
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_2d_i


    ! dataset writer for 2d double array
    subroutine h5_write_dataset_2d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        double precision, dimension(:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 2
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_2d_d


    ! dataset writer for 2d double array
    subroutine h5_write_dataset_2d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 2
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        if(CUSTOM_REAL == 4) then
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        else
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        endif
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_2d_r


    ! dataset writer for 2d double array
    subroutine h5_write_dataset_2d_r_no_group(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 2
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        if(CUSTOM_REAL == 4) then
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        else
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        endif
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_2d_r_no_group


    subroutine h5_write_dataset_4d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 4
        integer(HSIZE_T), dimension(4)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        if(CUSTOM_REAL == 4) then
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        else
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
            if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
        endif
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name

    end subroutine h5_write_dataset_4d_r


    ! dataset writer for 2d character array
    subroutine h5_write_dataset_2d_c(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        character(len=*), dimension(:,:), intent(in) :: data
        integer(HID_T)                    :: dspace_id ! dataspace id is local.
        integer                           :: rank = 2
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
        call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_CHARACTER, dspace_id, dataset_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
        call h5dwrite_f(dataset_id, H5T_NATIVE_CHARACTER, data, dim, error)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', dataset_name
    end subroutine h5_write_dataset_2d_c


    ! set attribute to a dataset
    subroutine h5_add_attribute_i(this, attribute_name, data)
        type(h5io), intent(in) :: this

        integer(HID_T)                             :: aspace_id, atype_id, attr_id ! attribute id
        integer                                    :: rank = 1 ! attribute rank
        integer, dimension(:), intent(in)          :: data
        integer(HSIZE_T)                           :: taglen
        character(len=*), intent(in)               :: attribute_name
        integer(HSIZE_T), dimension(1)             :: dim

        dim    = shape(data)
        taglen = len(trim(attribute_name))

        call h5screate_simple_f(rank, dim, aspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 screate failed for attribute, ', attribute_name
        call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error) ! for a string tag of attribute value
        if (error /= 0) write(*,*) 'hdf5 tcopy failed for attribute, ', attribute_name
        call h5tset_size_f(atype_id, taglen, error)
        if (error /= 0) write(*,*) 'hdf5 set_size failed for attribute, ', attribute_name
        ! here the attribute is written on the current openning dataset
        call h5acreate_f(dataset_id, trim(attribute_name), atype_id, aspace_id, attr_id, error)
        if (error /= 0) write(*,*) 'hdf5 acreate failed for attribute, ', attribute_name
        call h5awrite_f(attr_id, atype_id, data, dim, error) ! write
        if (error /= 0) write(*,*) 'hdf5 awrite failed for attribute, ', attribute_name

        call h5sclose_f(aspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for, ', attribute_name
        call h5aclose_f(attr_id,    error)
        if (error /= 0) write(*,*) 'hdf5 aclose failed for, ', attribute_name
        call h5tclose_f(atype_id,   error)
        if (error /= 0) write(*,*) 'hdf5 tclose failed for, ', attribute_name
    end subroutine h5_add_attribute_i


    !
    ! parallel routines
    !
    subroutine h5_set_mpi_info(this, comm, info, rank, nproc)
        type(h5io), intent(in) :: this
        integer, intent(in) :: comm, info, rank, nproc

        this_comm  = comm
        this_info  = info
        this_rank  = rank
        total_proc = nproc
    end subroutine h5_set_mpi_info

    subroutine h5_set_sieve_buffer_size(this)
        type(h5io), intent(in) :: this
        integer(hsize_t)       :: buf_size  = 1*1024*1024
        integer(hsize_t)       :: alig_size = 1*1024*1024
        call h5pset_sieve_buf_size_f(fplist_id, buf_size, error) ! buf_size may vary depending on machiens
        call h5pset_alignment_f(fplist_id, buf_size, alig_size, error)
    end subroutine h5_set_sieve_buffer_size


    subroutine h5_set_buffer_size(this)
        type(h5io), intent(in) :: this
        integer(hsize_t)       :: buf_size = 1*1024*1024
        call h5pset_buffer_f(plist_id, buf_size, error)
    end subroutine h5_set_buffer_size


    subroutine h5_create_file_prop_list(this, if_collective)
        type(h5io), intent(in) :: this
        logical, intent(in) :: if_collective

        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_FILE_ACCESS_F, fplist_id, error)
        if (error /= 0) write(*,*) 'hdf5 create plist failed.'

        call h5pset_all_coll_metadata_ops_f(fplist_id, if_collective, error)
        call h5pset_coll_metadata_write_f(fplist_id, if_collective, error)

        call h5pset_fapl_mpio_f(fplist_id, this_comm, this_info, error)
        if (error /= 0) write(*,*) 'hdf5 set_fapl failed.'
    end subroutine h5_create_file_prop_list


    subroutine h5_open_group_prop_list(this)
        type(h5io), intent(in) :: this
        logical, parameter :: if_collective = .true.
        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_GROUP_ACCESS_F, gplist_id, error)
        if (error /= 0) write(*,*) 'hdf5 create group plist failed.'
        call h5pset_all_coll_metadata_ops_f(gplist_id, if_collective, error)
        call h5pset_coll_metadata_write_f(gplist_id, if_collective, error)

        if (error /= 0) write(*,*) 'hdf5 group set_all_coll_metadata_ops failed.'
    end subroutine h5_open_group_prop_list


    subroutine h5_create_group_prop_list(this)
        type(h5io), intent(in) :: this
        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_GROUP_CREATE_F, gplist_id, error)
        if (error /= 0) write(*,*) 'hdf5 create group plist failed.'
    end subroutine h5_create_group_prop_list


    subroutine h5_create_dataset_prop_list(this, if_collect)
        type(h5io), intent(in) :: this
        logical, intent(in)    :: if_collect
        ! Setup file access property list with parallel I/O access.
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        if (error /= 0) write(*,*) 'hdf5 create dataset plist failed.'

        ! use larger buffer size from the default
        call h5_set_buffer_size(this)
        if (if_collect) then
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        else
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset set_dxpl_mpio failed.'
    end subroutine h5_create_dataset_prop_list


    subroutine h5_close_prop_list(this,dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        !print *,"DEBUG closing dataset prop: ", dataset_name
        ! check if the latest dataset transfer has been done in collective mode
        call h5_check_collective(this, dataset_name)
        call h5pclose_f(plist_id, error) ! property list can be closed soon
        if (error /= 0) write(*,*) 'hdf5 close_prop_list failed.'
    end subroutine h5_close_prop_list


    subroutine h5_close_prop_list_nocheck(this,dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        !print *,"DEBUG closing dataset prop: ", dataset_name
        ! check if the latest dataset transfer has been done in collective mode
        call h5pclose_f(plist_id, error) ! property list can be closed soon
        if (error /= 0) write(*,*) 'hdf5 close_prop_list failed.'
    end subroutine h5_close_prop_list_nocheck



    subroutine h5_close_file_prop_list(this)
        type(h5io), intent(in) :: this
        call h5pclose_f(fplist_id, error) ! property list can be closed soon
        if (error /= 0) write(*,*) 'hdf5 close_file_prop_list failed.'
    end subroutine h5_close_file_prop_list


    subroutine h5_close_group_prop_list(this)
        type(h5io), intent(in) :: this
        call h5pclose_f(gplist_id, error) ! property list can be closed soon
        if (error /= 0) write(*,*) 'hdf5 close_prop_list failed.'
    end subroutine h5_close_group_prop_list


    subroutine h5_create_file_p(this)
        type(h5io), intent(in) :: this
        call h5_create_file_prop_list(this, .false.)
        call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, error, creation_prp=H5P_DEFAULT_F, access_prp=fplist_id)
        if (error /= 0) write(*,*) 'hdf5 create file p failed.'
    end subroutine h5_create_file_p


    subroutine h5_create_file_p_collect(this)
        type(h5io), intent(in) :: this
        call h5_create_file_prop_list(this, .true.)
        call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, error, creation_prp=H5P_DEFAULT_F, access_prp=fplist_id)
        if (error /= 0) write(*,*) 'hdf5 create file p failed.'
    end subroutine h5_create_file_p_collect


    subroutine h5_open_file_p(this)
        type(h5io), intent(in) :: this
        call h5_create_file_prop_list(this,.false.)
        call h5fopen_f(file_path, H5F_ACC_RDWR_F, file_id, error, access_prp=fplist_id)
        if (error /= 0) write(*,*) 'hdf5 open file p failed.'
    end subroutine h5_open_file_p


    subroutine h5_open_file_p_collect(this)
        type(h5io), intent(in) :: this
        call h5_create_file_prop_list(this, .true.)
        call h5_set_sieve_buffer_size(this)
        call h5fopen_f(file_path, H5F_ACC_RDWR_F, file_id, error, access_prp=fplist_id)
        if (error /= 0) write(*,*) 'hdf5 open file p failed.'
    end subroutine h5_open_file_p_collect


    subroutine h5_create_dataset_gen(this, dataset_name, dim_in, rank, dtype_id)
        type(h5io), intent(in)       :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:), intent(in)  :: dim_in

        integer(HSIZE_T), dimension(:), allocatable :: dim
        integer, intent(in)                :: dtype_id ! 1:int, 4:real4, 8:real8,
        integer, intent(in)                :: rank

        integer(HID_T)                     :: dspace_id
        integer :: i
        logical :: if_chunk = .true.

        dim = dim_in ! convert data type

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name

        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
        ! chunk size setting
        !do i = 1, rank
        !    if (dim(i) <= 0) then
        !        if_chunk = .false.
        !        print *, "dataset not chunk set: ", dataset_name
        !    endif
        !enddo
        !if (if_chunk) call h5pset_chunk_f(plist_id,rank,dim,error)


        if (dtype_id == 0) then ! bool uses integer
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 1) then ! integer
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 2) then ! character
            call h5dcreate_f(file_id, trim(dataset_name), str_type, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 4) then  !  real
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 8) then ! double
            call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else
            print *, "specified dtype_id is not implemented yet for hdf5 io. aborting..."
            stop
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name

        call h5_close_prop_list_nocheck(this,dataset_name)

        call h5dclose_f(dataset_id,error)
        if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name

    end subroutine h5_create_dataset_gen


    subroutine h5_create_dataset_gen_in_group(this, dataset_name, dim_in, rank, dtype_id)
        type(h5io), intent(in)       :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:), intent(in)  :: dim_in

        integer(HSIZE_T), dimension(:), allocatable :: dim
        integer, intent(in)                :: dtype_id ! 1:int, 4:real4, 8:real8,
        integer, intent(in)                :: rank

        integer(HID_T)                     :: dspace_id

        dim = dim_in ! convert data type

        call h5screate_simple_f(rank, dim, dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name

        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if (dtype_id == 1) then ! integer
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 2) then ! character
            call h5dcreate_f(group_id, trim(dataset_name), str_type, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 4) then  !  real
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else if (dtype_id == 8) then ! double
            call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error, &
                                     dcpl_id=plist_id)
        else
            print *, "specified dtype_id is not implemented yet for hdf5 io. aborting..."
            stop
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name

        call h5_close_prop_list_nocheck(this,dataset_name)

        call h5dclose_f(dataset_id,error)
        if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
        call h5sclose_f(dspace_id, error)
        if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name

    end subroutine h5_create_dataset_gen_in_group



    subroutine h5_create_group_p(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name

        call h5_create_group_prop_list(this)
        call h5gcreate_f(file_id, trim(group_name), group_id, error, gcpl_id=gplist_id)
        if (error /= 0) write(*,*) 'error while creating a group, ', group_name
        call h5gclose_f(group_id, error)
        if (error /= 0) write(*,*) 'error while closing a group, ' , group_name
        call h5_close_group_prop_list(this)
    end subroutine h5_create_group_p


    subroutine h5_open_group_p(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
        call h5_open_group_prop_list(this)
        call h5gopen_f(file_id, trim(group_name), group_id, error, gapl_id=gplist_id)
        if (error /= 0) write(*,*) 'hdf5 open group failed for, ', group_name
        call h5_close_group_prop_list(this)
        store_group_name = group_name
    end subroutine h5_open_group_p


    subroutine h5_read_dataset_p_scalar_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, intent(out) :: data
        integer, dimension(1) :: rdata
        integer(HSIZE_T), dimension(1)    :: dim
        dim = shape(rdata)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, rdata, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)

        data = rdata(1)
    end subroutine h5_read_dataset_p_scalar_i


    subroutine h5_read_dataset_p_scalar_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), intent(out) :: data
        real(kind=CUSTOM_REAL), dimension(1) :: rdata
        integer(HSIZE_T), dimension(1)    :: dim
        dim = shape(rdata)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, rdata, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, rdata, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)

        data = rdata(1)
    end subroutine h5_read_dataset_p_scalar_r


    subroutine h5_read_dataset_p_1d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:), intent(inout) :: data
        integer(HSIZE_T), dimension(1)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_1d_i


    subroutine h5_read_dataset_p_1d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: data
        integer(HSIZE_T), dimension(1)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_1d_r


    subroutine h5_read_dataset_p_1d_l(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        logical, dimension(:), intent(inout) :: data
        integer, dimension(:), allocatable :: ldata
        integer(HSIZE_T), dimension(1)    :: dim
        integer :: lsize
        dim = shape(data)
        lsize = size(data)
        allocate(ldata(lsize),stat=error)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, ldata, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)

        call int_array2bool(this, ldata, data)

        deallocate(ldata, stat=error)
    end subroutine h5_read_dataset_p_1d_l


    subroutine h5_read_dataset_p_2d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_2d_i


    subroutine h5_read_dataset_p_2d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_2d_d


    subroutine h5_read_dataset_p_2d_c(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        character(len=*), dimension(:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_CHARACTER, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_2d_c


    subroutine h5_read_dataset_p_2d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(2)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_2d_r


    subroutine h5_read_dataset_p_3d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(3)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_3d_i


    subroutine h5_read_dataset_p_3d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(3)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_3d_r


    subroutine h5_read_dataset_p_4d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(4)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_4d_r


    subroutine h5_read_dataset_p_5d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(5)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        if (CUSTOM_REAL == 4) then
          call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
        else
          call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
        endif
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_5d_r


    subroutine h5_read_attribute_p(this, attribute_name, dataset_name, data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: attribute_name
        character(len=*), intent(in)         :: dataset_name
        integer(HID_T)                       :: attr_id
        integer, dimension(:), intent(inout) :: data
        integer(HSIZE_T), dimension(1)       :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        call h5aopen_f(dataset_id, attribute_name, attr_id, error)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, dim, error)
        call h5aclose_f(attr_id,error)

        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_attribute_p


    subroutine h5_read_dataset_p_4d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:,:,:), intent(inout) :: data
        integer(HSIZE_T), dimension(4)    :: dim
        dim = shape(data)

        call h5dopen_f(group_id, dataset_name, dataset_id, error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dataset_id, error)
    end subroutine h5_read_dataset_p_4d_i


    subroutine h5_read_dataset_scalar_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, intent(inout), target                              :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        !dim = shape(data)
        dim = (/1/)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_scalar_i_collect_hyperslab


    subroutine h5_read_dataset_scalar_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), intent(inout), target               :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        !dim = shape(data)
        dim = (/1/)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if(CUSTOM_REAL==4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_scalar_r_collect_hyperslab



    subroutine h5_read_dataset_1d_l_collect_hyperslab(this, dataset_name, data_out, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:), allocatable, target                  :: data
        logical, dimension(:), intent(inout), target                :: data_out
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data_out)
        offset = offset_in ! convert data type
        allocate(data(dim(1)))

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

        call int_array2bool(this, data, data_out)
        deallocate(data)
    end subroutine h5_read_dataset_1d_l_collect_hyperslab


    subroutine h5_read_dataset_1d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:), intent(inout), target                :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_1d_i_collect_hyperslab


    ! store local 1d array to global 1d array
    subroutine h5_read_dataset_1d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(inout), target :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_1d_r_collect_hyperslab




    ! store local 1d array to global 1d array
    subroutine h5_read_dataset_1d_r_collect_hyperslab_in_group(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(inout), target :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_1d_r_collect_hyperslab_in_group


    subroutine h5_read_dataset_2d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:), intent(inout), target              :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_2d_i_collect_hyperslab


    subroutine h5_read_dataset_2d_i_collect_hyperslab_in_group(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:), intent(inout), target              :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_2d_i_collect_hyperslab_in_group


    ! store local 2d array to global 2d array
    subroutine h5_read_dataset_2d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(inout), target  :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)! nrec
        count(2)  = dim(2) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_2d_r_collect_hyperslab


    subroutine h5_read_dataset_2d_d_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        double precision, dimension(:,:), intent(inout), target        :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)! nrec
        count(2)  = dim(2) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_read_dataset_2d_d_collect_hyperslab


    subroutine h5_read_dataset_3d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:,:), intent(inout), target            :: data
        integer                                                     :: rank = 3
        integer(HSIZE_T), dimension(3)                              :: dim
        integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial


        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_read_dataset_3d_i_collect_hyperslab


    subroutine h5_read_dataset_3d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:), intent(inout), target  :: data
        integer                                                     :: rank = 3
        integer(HSIZE_T), dimension(3)                              :: dim
        integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial


        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_read_dataset_3d_r_collect_hyperslab


    subroutine h5_read_dataset_4d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:,:,:), intent(inout), target          :: data
        integer                                                     :: rank = 4
        integer(HSIZE_T), dimension(4)                              :: dim
        integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial
        count(4) = dim(4) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_read_dataset_4d_i_collect_hyperslab




    subroutine h5_read_dataset_4d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout), target  :: data
        integer                                                     :: rank = 4
        integer(HSIZE_T), dimension(4)                              :: dim
        integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial
        count(4) = dim(4) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_read_dataset_4d_r_collect_hyperslab


    subroutine h5_read_dataset_5d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout), target  :: data
        integer                                                     :: rank = 5
        integer(HSIZE_T), dimension(5)                              :: dim
        integer(HSIZE_T), dimension(5)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(5)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial
        count(4) = dim(4) ! NSTEP partial
        count(5) = dim(5) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_read_dataset_5d_r_collect_hyperslab


    !
    ! parallel write routines
    !


    !
    ! collective writers
    !

    subroutine h5_write_dataset_p_1d_i(this, dataset_name, data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: dataset_name
        integer, dimension(:), intent(in), target    :: data
        integer                              :: rank = 1
        integer(HSIZE_T), dimension(1)       :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! add dummy 0 for no_element array
        if (size(data) == 0) then
            dim = 1
        endif

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)


        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        if (size(data) == 0) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, dummy_1d_array, dim, error,xfer_prp=plist_id)
        else
            f_ptr = c_loc(data(1))
            call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                            xfer_prp=plist_id)
       endif
       if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_1d_i


    ! write 1d integer and attribute
    subroutine h5_write_dataset_p_1d_ia(this, dataset_name, data, attribute_name, attr_data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: dataset_name
        character(len=*), intent(in)         :: attribute_name
        integer, dimension(:), intent(in), target    :: data
        integer, intent(in)                  :: attr_data
        integer                              :: rank = 1
        integer(HSIZE_T), dimension(1)       :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! add dummy 0 for no_element array
        if (size(data) == 0) then
            dim = 1
        endif

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        if (size(data) == 0) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, dummy_1d_array, dim, error,xfer_prp=plist_id)
        else
            f_ptr = c_loc(data(1))
            call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                            xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        ! add attribute
        call h5_add_attribute_i(this, attribute_name, (/attr_data/))

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_1d_ia


    subroutine h5_write_dataset_p_1d_r(this, dataset_name, data)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name


        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, CUSTOM_REAL, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1))
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                            xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                           xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_1d_r


    subroutine h5_write_dataset_p_2d_i(this, dataset_name, data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: dataset_name
        integer, dimension(:,:), intent(in), target    :: data
        integer                              :: rank = 2
        integer(HSIZE_T), dimension(2)       :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                        xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_2d_i


    subroutine h5_write_dataset_p_2d_ia(this, dataset_name, data, attribute_name, attr_data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: dataset_name
        character(len=*), intent(in)         :: attribute_name
        integer, dimension(:,:), intent(in), target    :: data
        integer, intent(in)                  :: attr_data
        integer                              :: rank = 2
        integer(HSIZE_T), dimension(2)       :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                        xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        ! add attribute
        call h5_add_attribute_i(this, attribute_name, (/attr_data/))

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_2d_ia


    subroutine h5_write_dataset_p_2d_r(this, dataset_name, data)
        type(h5io), intent(in)                                        :: this
        character(len=*), intent(in)                                  :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target    :: data
        integer                                                       :: rank = 2
        integer(HSIZE_T), dimension(2)                                :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, CUSTOM_REAL, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1))
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                            xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                           xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_2d_r


    subroutine h5_write_dataset_p_2d_d(this, dataset_name, data)
        type(h5io), intent(in)                                        :: this
        character(len=*), intent(in)                                  :: dataset_name
        double precision, dimension(:,:), intent(in), target    :: data
        integer                                                       :: rank = 2
        integer(HSIZE_T), dimension(2)                                :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 8, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                       xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_2d_d


    subroutine h5_write_dataset_p_3d_i(this, dataset_name, data)
        type(h5io), intent(in)                           :: this
        character(len=*), intent(in)                     :: dataset_name
        integer, dimension(:,:,:), intent(in), target    :: data
        integer                                          :: rank = 3
        integer(HSIZE_T), dimension(3)                   :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1,1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                        xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_3d_i


    subroutine h5_write_dataset_p_3d_r(this, dataset_name, data)
        type(h5io), intent(in)                                        :: this
        character(len=*), intent(in)                                  :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:), intent(in), target  :: data
        integer                                                       :: rank = 3
        integer(HSIZE_T), dimension(3)                                :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, CUSTOM_REAL, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1,1))
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                            xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                           xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_3d_r


    subroutine h5_write_dataset_p_4d_i(this, dataset_name, data)
        type(h5io), intent(in)                           :: this
        character(len=*), intent(in)                     :: dataset_name
        integer, dimension(:,:,:,:), intent(in), target  :: data
        integer                                          :: rank = 4
        integer(HSIZE_T), dimension(4)                   :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1,1,1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                        xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_4d_i


    subroutine h5_write_dataset_p_4d_r(this, dataset_name, data)
        type(h5io), intent(in)                                         :: this
        character(len=*), intent(in)                                   :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in), target :: data
        integer                                                        :: rank = 4
        integer(HSIZE_T), dimension(4)                                 :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, CUSTOM_REAL, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1,1,1))
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                            xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                           xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_4d_r


    subroutine h5_write_dataset_p_5d_r(this, dataset_name, data)
        type(h5io), intent(in)                                            :: this
        character(len=*), intent(in)                                      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in), target  :: data
        integer                                                           :: rank = 5
        integer(HSIZE_T), dimension(5)                                    :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, CUSTOM_REAL, file_id)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


        ! write array using fortran pointer
        f_ptr = c_loc(data(1,1,1,1,1))
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                            xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                           xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
    end subroutine h5_write_dataset_p_5d_r


    subroutine h5_write_dataset_p_1d_l(this, dataset_name, data_in)
        ! writer for logical array
        ! logical array will be converted to integer array before write
        type(h5io), intent(in)                     :: this
        character(len=*), intent(in)               :: dataset_name
        logical, dimension(:), intent(in)          :: data_in
        integer, dimension(:), allocatable, target :: data

        integer(HID_T)                     :: dspace_id ! dataspace id is local.
        integer                            :: rank = 1
        integer(HSIZE_T), dimension(1)     :: dim

        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        dim = shape(data_in)

        ! create dataset
        call create_dataset_collect(this, dataset_name, dim, rank, 1, file_id)

        allocate(data(size(data_in)), stat=error)
        call bool_array2integer(this, data_in, data)

        write(tempstr, "(i6.6)") this_rank
        group_name = gname_proc_head // trim(tempstr)

        !! create datasets of all processes for independent write
        call h5_open_group(this, group_name)
        call h5_open_dataset(this,trim(dataset_name))

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ! write array using fortran pointer
        f_ptr = c_loc(data(1))
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                        xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5_close_dataset(this)
        call h5_close_group(this)
        deallocate(data, stat=error)
    end subroutine h5_write_dataset_p_1d_l


    ! from 1d local array to 2d global array
    subroutine h5_write_dataset_1d_to_2d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(2), intent(in)                           :: offset_in ! the position where the datablock is inserted
        integer(HSSIZE_T), dimension(2)                             :: offset
        logical                                                     :: if_collect

        dim = size(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! size of data array inserted.
        count(1)  = dim(1)
        count(2)  = 1

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_1d_to_2d_r_collect_hyperslab


    ! store local 2d array to global 3d array
    subroutine h5_write_dataset_2d_to_3d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target  :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = 3! NDIM
        count(2)  = dim(2) ! NSTEP partial
        count(3)  = 1 ! nrec

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_2d_to_3d_r_collect_hyperslab


    subroutine h5_write_dataset_1d_l_collect_hyperslab(this, dataset_name, data_in, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        logical, dimension(:), intent(in), target                   :: data_in
        integer, dimension(:), allocatable, target                  :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical, intent(in)                                         :: if_collect

        dim = shape(data_in)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! if_collect == .false., this function gather all data from procs then write in a file at once

        ! select a place where data is inserted.
        count(1)  = dim(1)
        allocate(data(dim(1)))

        ! convert logical array to  integer array
        call bool_array2integer(this, data_in, data)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

        deallocate(data)

    end subroutine h5_write_dataset_1d_l_collect_hyperslab



    ! store local 1d array to global 1d array
    subroutine h5_write_dataset_1d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:), intent(in), target                   :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_1d_i_collect_hyperslab


    subroutine h5_write_dataset_1d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_1d_r_collect_hyperslab


    ! store local 1d array to global 1d array
    subroutine h5_write_dataset_1d_r_collect_hyperslab_in_group(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
        integer                                                     :: rank = 1
        integer(HSIZE_T), dimension(1)                              :: dim
        integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_1d_r_collect_hyperslab_in_group


    ! store local 2d array to global 2d array
    subroutine h5_write_dataset_2d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:), intent(in), target  :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_2d_i_collect_hyperslab


    subroutine h5_write_dataset_2d_i_collect_hyperslab_in_group(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:), intent(in), target  :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_2d_i_collect_hyperslab_in_group


    ! store local 2d array to global 2d array
    subroutine h5_write_dataset_2d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target  :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)! nrec
        count(2)  = dim(2) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_2d_r_collect_hyperslab


    subroutine h5_write_dataset_2d_d_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        double precision, dimension(:,:), intent(in), target        :: data
        integer                                                     :: rank = 2
        integer(HSIZE_T), dimension(2)                              :: dim
        integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)! nrec
        count(2)  = dim(2) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_2d_d_collect_hyperslab


    subroutine h5_write_dataset_3d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                        :: this
        character(len=*), intent(in)                  :: dataset_name
        integer, dimension(:,:,:), intent(in), target :: data
        integer                                       :: rank = 3
        integer(HSIZE_T), dimension(3)                :: dim
        integer(HSIZE_T), dimension(3)                :: count ! size of hyperslab
        integer, dimension(:), intent(in)             :: offset_in
        integer(HSSIZE_T), dimension(3)               :: offset ! the position where the datablock is inserted
        logical                                       :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_write_dataset_3d_i_collect_hyperslab

    ! store local 3d array to global 3d array
    subroutine h5_write_dataset_3d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:), intent(in), target  :: data
        integer                                                     :: rank = 3
        integer(HSIZE_T), dimension(3)                              :: dim
        integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1) = dim(1) ! NDIM
        count(2) = dim(2) ! nrec
        count(3) = dim(3) ! NSTEP partial

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)

    end subroutine h5_write_dataset_3d_r_collect_hyperslab


    subroutine h5_write_dataset_4d_i_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                      :: this
        character(len=*), intent(in)                                :: dataset_name
        integer, dimension(:,:,:,:), intent(in), target             :: data
        integer                                                     :: rank = 4
        integer(HSIZE_T), dimension(4)                              :: dim
        integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
        integer, dimension(:), intent(in)                           :: offset_in
        integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted
        logical                                                     :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)
        count(3)  = dim(3)
        count(4)  = dim(4)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                        file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_4d_i_collect_hyperslab


    subroutine h5_write_dataset_4d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                         :: this
        character(len=*), intent(in)                                   :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in), target :: data
        integer                                                        :: rank = 4
        integer(HSIZE_T), dimension(4)                                 :: dim
        integer(HSIZE_T), dimension(4)                                 :: count ! size of hyperslab
        integer, dimension(:), intent(in)                              :: offset_in
        integer(HSSIZE_T), dimension(4)                                :: offset ! the position where the datablock is inserted
        logical                                                        :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)
        count(3)  = dim(3)
        count(4)  = dim(4)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_4d_r_collect_hyperslab


    subroutine h5_write_dataset_5d_r_collect_hyperslab(this, dataset_name, data, offset_in, if_collect)
        type(h5io), intent(in)                                           :: this
        character(len=*), intent(in)                                     :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in), target :: data
        integer                                                          :: rank = 5
        integer(HSIZE_T), dimension(5)                                   :: dim
        integer(HSIZE_T), dimension(5)                                   :: count ! size of hyperslab
        integer, dimension(:), intent(in)                                :: offset_in
        integer(HSSIZE_T), dimension(5)                                  :: offset ! the position where the datablock is inserted
        logical                                                          :: if_collect

        dim = shape(data)
        offset = offset_in ! convert data type

        ! open dataset
        call h5_open_dataset2(this,trim(dataset_name))

        ! select a place where data is inserted.
        count(1)  = dim(1)
        count(2)  = dim(2)
        count(3)  = dim(3)
        count(4)  = dim(4)
        count(5)  = dim(5)

        ! select hyperslab in the file
        call h5screate_simple_f(rank,count, mem_dspace_id, error)
        call h5dget_space_f(dataset_id, file_dspace_id, error)
        call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)

        call h5_create_dataset_prop_list(this,if_collect)

        ! write array using fortran pointer
        if (CUSTOM_REAL == 4) then
            call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        else
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                            file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
        endif
        if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name

        call h5_close_prop_list(this,dataset_name)
        call h5sclose_f(mem_dspace_id, error)
        call h5sclose_f(file_dspace_id, error)
        call h5_close_dataset(this)
    end subroutine h5_write_dataset_5d_r_collect_hyperslab


    !
    ! other utilities
    !
    subroutine bool_array2integer(this, boolarray, intarray)
        type(h5io), intent(in) :: this
        logical, dimension(:), intent(in) :: boolarray
        integer, dimension(:), intent(out) :: intarray
        integer :: array_size, i

        array_size = size(boolarray)
        intarray(:) = 0
        do i = 1, array_size
            if (boolarray(i)) then
                intarray(i) = 1
            endif
        enddo
    end subroutine bool_array2integer


    subroutine int_array2bool(this,intarray, boolarray)
        type(h5io), intent(in)             :: this
        integer, dimension(:), intent(in)  :: intarray
        logical, dimension(:), intent(out) :: boolarray
        integer                            :: array_size,i

        array_size = size(intarray)
        boolarray(:) = .false.

        do i = 1, array_size
            if (intarray(i) /=0) then
                boolarray(i) = .true.
            endif
        enddo
    end subroutine int_array2bool


    subroutine h5_gather_dsetsize(this, dimin, data_rank, all_dim)
        use mpi
        type(h5io), intent(in)                      :: this
        integer, intent(in)                         :: data_rank
        integer(HSIZE_T), dimension(:), intent(in)  :: dimin
        integer, dimension(data_rank)               :: dim
        integer, dimension(data_rank,0:total_proc-1), intent(out)      :: all_dim

        dim = dimin ! convert integer 8(HSIZE_T) to 4
        if (data_rank==1) then
           call gather_all_all_singlei(dim(1),all_dim,total_proc)
            if (this_rank == 0) then
!                print *, "rank: ", this_rank &
!                 ,"dimin", dimin, ", kind: ", kind(dimin) &
!                 ,"dim", dim, ", kind: ", kind(dim) &
!                 ,"shape dimin", shape(dimin) &
!                 ,"data_rank", data_rank &
!                 ,"all dim", all_dim &
!                 ,"shape all dim", shape(all_dim) &
!                 ,"kind all_dim", kind(all_dim(1,0)) &
!                 ,"total proc " , total_proc &
!                 ,"kind mpi_integer ", kind(MPI_INTEGER)
            endif
        else
!            call gather_all_i(dim, data_rank, all_dim, data_rank, total_proc)
            call gather_all_all_i(dim, data_rank, all_dim, data_rank, total_proc)

        endif

        call synchronize_all()
    end subroutine


    subroutine create_dataset_collect(this, dataset_name, dim, data_rank, dtype_id, base_id)
        type(h5io), intent(in)                      :: this
        character(len=*), intent(in)                :: dataset_name
        integer(HSIZE_T), dimension(:), intent(in)  :: dim
        integer, intent(in)                         :: dtype_id ! 1:int, 4:real4, 8:real8,
        integer, intent(in)                         :: data_rank
        integer(HID_T), intent(in)                  :: base_id ! base (file_id or group_id) of dataset creation

        integer, dimension(data_rank,0:total_proc-1)  :: all_dim
        integer(HSIZE_T), dimension(data_rank)        :: dim_h5
        integer                                       :: iproc
        character(len=128)                            :: group_and_dataset_name
        integer(HID_T)                                :: dspace_id
        character(len=10) :: tempstr
        character(len=5)  :: gname_proc_head = "proc_"
        character(len=64) :: group_name

        ! gather dataset dimension of other processors
        call h5_gather_dsetsize(this, dim, data_rank, all_dim)
        ! create all datasets of all process by rank0

        !print *, dataset_name, ": debug returned all dim", all_dim
        do iproc = 0, total_proc -1
            write(tempstr, "(i6.6)") iproc
            group_name = gname_proc_head // trim(tempstr)
            group_and_dataset_name = trim(group_name) //"/"// dataset_name

            dim_h5 = all_dim(:,iproc)
            call h5screate_simple_f(data_rank, dim_h5, dspace_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name

            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

            ! This is required for this data pattern
            call H5Pset_alloc_time_f(plist_id, H5D_ALLOC_TIME_EARLY_F, error)

            if (dtype_id == 1) then
                call h5dcreate_f(base_id, trim(group_and_dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error, &
                                         dcpl_id=plist_id)
            else if (dtype_id == 4) then
                call h5dcreate_f(base_id, trim(group_and_dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error, &
                                         dcpl_id=plist_id)
            else
                call h5dcreate_f(base_id, trim(group_and_dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error, &
                                         dcpl_id=plist_id)
            endif
            if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name

            call h5_close_prop_list(this,dataset_name)
            call h5dclose_f(dataset_id,error)
            if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
            call h5sclose_f(dspace_id, error)
            if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name

        enddo

        call synchronize_all()
    end subroutine create_dataset_collect


!
! higher level utilities
!

    subroutine write_attenuation_file_in_h5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)
        use shared_parameters, only: NPROC, LOCAL_PATH
        use constants, only: CUSTOM_REAL,myrank,N_SLS,NGLLX,NGLLY,NGLLZ
        implicit none

        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor_kappa

        integer, dimension(4) :: dims
        integer :: nspec

        ! offset arrays
        integer, dimension(0:NPROC-1) :: offset_nspec

        ! hdf5 valiables
        character(len=64) :: filename, dset_name, tempstr
        type(h5io)        :: h5

        dims = shape(scale_factor)
        nspec=dims(4)

        ! prepare offset arrays
        call gather_all_all_singlei(nspec,offset_nspec,NPROC) ! n spec in each proc

        tempstr = "/external_mesh.h5"
        filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)
        h5 = h5io()

        ! initialize h5 object
        call h5_init(h5, filename)

        ! prepare dataset
        if(myrank == 0) then
            call h5_open_file(h5)
            dset_name = "scale_factor"
            call h5_create_dataset_gen(h5, dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 4, CUSTOM_REAL)
            dset_name = "scale_factor_kappa"
            call h5_create_dataset_gen(h5, dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 4, CUSTOM_REAL)
            dset_name = "factor_common"
            call h5_create_dataset_gen(h5, dset_name, (/N_SLS,NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 5, CUSTOM_REAL)
            dset_name = "factor_common_kappa"
            call h5_create_dataset_gen(h5, dset_name, (/N_SLS,NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 5, CUSTOM_REAL)
            call h5_close_file(h5)
        endif

        call synchronize_all()

        ! write
        call h5_open_file_p_collect(h5)
        dset_name = "scale_factor"
        call h5_write_dataset_4d_r_collect_hyperslab(h5,dset_name,&
                            scale_factor,(/0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
        dset_name = "scale_factor_kappa"
        call h5_write_dataset_4d_r_collect_hyperslab(h5,dset_name,&
                            scale_factor_kappa,(/0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
        dset_name = "factor_common"
        call h5_write_dataset_5d_r_collect_hyperslab(h5,dset_name,&
                            factor_common,(/0,0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
        dset_name = "factor_common_kappa"
        call h5_write_dataset_5d_r_collect_hyperslab(h5,dset_name,&
                            factor_common_kappa,(/0,0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)

        call h5_close_file_p(h5)

    end subroutine write_attenuation_file_in_h5


    subroutine read_attenuation_file_in_h5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)
        use shared_parameters, only: NPROC, LOCAL_PATH
        use constants, only: CUSTOM_REAL,myrank

        implicit none

        ! offset array
        integer, dimension(0:NPROC-1) :: offset_nspec

        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
        real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor_kappa

        ! hdf5 valiables
        character(len=64) :: fname, prname, gname, dset_name, tempstr
        type(h5io)        :: h5

        tempstr = "/external_mesh.h5"
        fname = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)
        h5 = h5io()

        ! initialize h5 object
        call h5_init(h5, fname)


        ! initialize h5 object
        call h5_init(h5, fname)
        call h5_open_file_p_collect(h5)
        ! read offset array
        call h5_read_dataset_1d_i_collect_hyperslab(h5,"offset_nspec",offset_nspec, (/0/), .true.)

        call h5_read_dataset_4d_r_collect_hyperslab(h5, "scale_factor", scale_factor, &
                                                    (/0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
        call h5_read_dataset_4d_r_collect_hyperslab(h5, "scale_factor_kappa", scale_factor_kappa,&
                                                    (/0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
        call h5_read_dataset_5d_r_collect_hyperslab(h5, "factor_common", factor_common, &
                                                    (/0,0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
        call h5_read_dataset_5d_r_collect_hyperslab(h5, "factor_common_kappa", factor_common_kappa, &
                                                    (/0,0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)

        call h5_close_file_p(h5)

    end subroutine read_attenuation_file_in_h5


    subroutine write_checkmesh_data_h5(dset_name,dump_array)
        use shared_parameters, only: LOCAL_PATH, NPROC
        use constants, only: CUSTOM_REAL, MAX_STRING_LEN, myrank

        implicit none

        real(kind=CUSTOM_REAL),dimension(:), intent(in) :: dump_array
        character(len=MAX_STRING_LEN), intent(in)       :: dset_name
        character(len=MAX_STRING_LEN)                   :: filename
        integer, dimension(0:NPROC-1)                   :: offset

        ! mpi variables
        integer :: info, comm

        ! hdf5 valiables
        character(len=64) :: tempstr
        type(h5io)        :: h5

        ! saves mesh file external_mesh.h5
        tempstr = "/external_mesh.h5"
        filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

        h5 = h5io()
        ! get mpi parameters
        call world_get_comm(comm)
        call get_info_null(info)

        ! initialize h5 object
        call h5_init(h5, filename)
        call h5_set_mpi_info(h5, comm, info, myrank, NPROC)

        ! get offset info
        call gather_all_all_singlei(size(dump_array), offset, NPROC)

        ! make dataset
        if(myrank==0) then
            call h5_open_file(h5)
            call h5_create_dataset_gen(h5, dset_name, (/sum(offset(:))/), 1, CUSTOM_REAL)
            call h5_close_file(h5)
        endif
        call synchronize_all()

        ! open file
        call h5_open_file_p_collect(h5)
        call h5_write_dataset_1d_r_collect_hyperslab(h5, dset_name, dump_array, (/sum(offset(0:myrank-1))/),.true.)
        call h5_close_file_p(h5)
        call h5_destructor(h5)

    end subroutine write_checkmesh_data_h5


    subroutine write_checkmesh_xdmf_h5(NSPEC_AB)
        use shared_parameters
        use constants
        implicit none

        integer, intent(in)                             :: NSPEC_AB
        character(len=64)                               :: fname_xdmf_checkmesh    = ""
        integer, dimension(0:NPROC-1)                   :: nelms, nnodes
        character(len=20)                               :: proc_str,type_str,nelm,nnode
        integer                                         :: iproc

        ! gather number of elements in each proc
        call gather_all_singlei(NSPEC_AB, nelms, NPROC)
        ! count and gather the number of controle nodes in each proc
        call gather_all_singlei(n_control_node, nnodes, NPROC)

        if (myrank==0) then

            ! writeout xdmf file for surface movie
            fname_xdmf_checkmesh = trim(OUTPUT_FILES)//"/checkmesh.xmf"

            open(unit=xdmf_vol, file=fname_xdmf_checkmesh, recl=512)

           ! definition of topology and geometry
            ! refer only control nodes (8 or 27) as a coarse output
            ! data array need to be extracted from full data array on gll points
            write(xdmf_vol,'(a)') '<?xml version="1.0" ?>'
            write(xdmf_vol,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
            write(xdmf_vol,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
            write(xdmf_vol,*) '<Domain>'
            write(xdmf_vol,*) '    <!-- mesh info -->'
            !write(xdmf_vol,*) '    <Grid Name="mesh" GridType="Collection"  CollectionType="Spatial">'
                nelm  = i2c(sum(nelms(:)))
                nnode = i2c(sum(nnodes(:)))

                write(xdmf_vol,*) '<Grid Name="mesh">'
                write(xdmf_vol,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm)//'">'
                write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                                       //trim(nelm)//' '//trim(i2c(8+1))//'">'
                write(xdmf_vol,*) '       ./DATABASES_MPI/Database.h5:/elm_conn_xdmf'
                write(xdmf_vol,*) '    </DataItem>'
                write(xdmf_vol,*) '</Topology>'
                write(xdmf_vol,*) '<Geometry GeometryType="XYZ">'
                write(xdmf_vol,*) '    <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                    //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nnode)//' 3">'
                write(xdmf_vol,*) '       ./DATABASES_MPI/Database.h5:/nodes_coords'
                write(xdmf_vol,*) '    </DataItem>'
                write(xdmf_vol,*) '</Geometry>'

                type_str = "res_Courant_number"
                write(xdmf_vol, *)  '    <Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Cell">'
                write(xdmf_vol, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                               //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nelm)//'">'
                write(xdmf_vol, *)  '            ./DATABASES_MPI/external_mesh.h5:/'//trim(type_str)
                write(xdmf_vol, *)  '        </DataItem>'
                write(xdmf_vol, *)  '    </Attribute>'

                type_str = "res_minimum_period"
                write(xdmf_vol, *)  '    <Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Cell">'
                write(xdmf_vol, *)  '        <DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                               //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nelm)//'">'
                write(xdmf_vol, *)  '            ./DATABASES_MPI/external_mesh.h5:/'//trim(type_str)
                write(xdmf_vol, *)  '        </DataItem>'
                write(xdmf_vol, *)  '    </Attribute>'

                write(xdmf_vol,*) '</Grid>'

            !write(xdmf_vol,*) '</Grid>'
            write(xdmf_vol,*) '</Domain>'
            write(xdmf_vol,*) '</Xdmf>'

            close(xdmf_vol)
        endif

    end subroutine write_checkmesh_xdmf_h5

    function i2c(k) result(str)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        character(len=20) str
        write (str, "(i20)") k
        str = adjustl(str)
    end function i2c

end module phdf5_utils


