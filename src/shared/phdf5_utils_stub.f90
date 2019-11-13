module phdf5_utils
    use constants, only: CUSTOM_REAL
    implicit none

    private
    public :: h5io, h5_init, h5_destructor, &
         h5_create_file, h5_open_file, h5_close_file, &
         h5_create_group, h5_open_group, h5_close_group, &
         h5_open_dataset, h5_close_dataset, &
         h5_write_dataset_1d_i, h5_write_dataset_1d_d, &
         h5_write_dataset_2d_d, h5_write_dataset_2d_i, h5_write_dataset_2d_c, &
         h5_add_attribute_i, &
         h5_set_mpi_info, h5_create_file_p, h5_open_file_p,  &
         h5_create_file_prop_list, h5_close_prop_list, h5_open_group_prop_list, h5_create_group_prop_list, &
         h5_create_dataset_prop_list, &
         h5_create_group_p, h5_open_group_p, &
         h5_open_dataset_p, &
         h5_read_dataset_p_1d_i, &
         h5_read_dataset_p_2d_i, h5_read_dataset_p_2d_d, h5_read_dataset_p_2d_c, &
         h5_read_attribute_p, &
         h5_set_group_name, &
         h5_write_dataset_p_1d_i, h5_write_dataset_p_1d_r, h5_write_dataset_p_1d_l, &
         h5_write_dataset_p_2d_i, h5_write_dataset_p_2d_r, &
         h5_write_dataset_p_3d_i, h5_write_dataset_p_3d_r, &
         h5_write_dataset_p_4d_i, h5_write_dataset_p_4d_r, &
         h5_write_dataset_p_5d_r, &
         bool_array2integer, int_array2bool, &
         h5_gather_dsetsize, h5_create_dataset_in_main_proc, h5_create_dataset_setter, &
         write_attenuation_file_in_h5, read_attenuation_file_in_h5


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
    end subroutine h5_init


    subroutine h5_destructor(this)
        type(h5io), intent(in) :: this
    end subroutine h5_destructor


    subroutine h5_create_file(this)
        type(h5io), intent(in) :: this

    end subroutine h5_create_file


    subroutine h5_open_file(this)
        type(h5io), intent(in) :: this
    end subroutine h5_open_file


    subroutine h5_close_file(this)
        type(h5io), intent(in) :: this
    end subroutine h5_close_file


    subroutine h5_create_group(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name

    end subroutine h5_create_group


    subroutine h5_open_group(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
    end subroutine h5_open_group


    subroutine h5_close_group(this)
        type(h5io), intent(in) :: this
    end subroutine h5_close_group


    subroutine h5_set_group_name(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
    end subroutine h5_set_group_name

    ! open dataset. group need to be opened before.
    subroutine h5_open_dataset(this, dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
    end subroutine h5_open_dataset


    subroutine h5_close_dataset(this)
        type(h5io), intent(in) :: this
    end subroutine h5_close_dataset


    ! dataset writer for 1d integer array
    subroutine h5_write_dataset_1d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:), intent(in) :: data
    end subroutine h5_write_dataset_1d_i


    ! dataset writer for 1d custom real array
    subroutine h5_write_dataset_1d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
    end subroutine h5_write_dataset_1d_d


    ! dataset writer for 2d integer array
    subroutine h5_write_dataset_2d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:,:), intent(in) :: data
    end subroutine h5_write_dataset_2d_i


    ! dataset writer for 2d double array
    subroutine h5_write_dataset_2d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        double precision, dimension(:,:), intent(in) :: data
    end subroutine h5_write_dataset_2d_d

    ! dataset writer for 2d character array
    subroutine h5_write_dataset_2d_c(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        character(len=*), dimension(:,:), intent(in) :: data
    end subroutine h5_write_dataset_2d_c


    ! set attribute to a dataset
    subroutine h5_add_attribute_i(this, attribute_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: attribute_name
        integer, dimension(:), intent(in)          :: data
    end subroutine h5_add_attribute_i


    !
    ! parallel routines
    !
    subroutine h5_set_mpi_info(this, comm, info, rank, nproc)
        type(h5io), intent(in) :: this
        integer, intent(in) :: comm, info, rank, nproc
    end subroutine h5_set_mpi_info

    subroutine h5_create_file_prop_list(this)
        type(h5io), intent(in) :: this
    end subroutine h5_create_file_prop_list


    subroutine h5_open_group_prop_list(this)
        type(h5io), intent(in) :: this
    end subroutine h5_open_group_prop_list


    subroutine h5_create_group_prop_list(this)
        type(h5io), intent(in) :: this
    end subroutine h5_create_group_prop_list


    subroutine h5_create_dataset_prop_list(this)
        type(h5io), intent(in) :: this
    end subroutine h5_create_dataset_prop_list


    subroutine h5_close_prop_list(this)
        type(h5io), intent(in) :: this
    end subroutine h5_close_prop_list


    subroutine h5_create_file_p(this)
        type(h5io), intent(in) :: this
    end subroutine h5_create_file_p


    subroutine h5_open_file_p(this)
        type(h5io), intent(in) :: this
    end subroutine h5_open_file_p


    subroutine h5_create_group_p(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
    end subroutine h5_create_group_p


    subroutine h5_open_group_p(this, group_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: group_name
    end subroutine h5_open_group_p


    subroutine h5_read_dataset_p_1d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:), intent(inout) :: data
    end subroutine h5_read_dataset_p_1d_i


    subroutine h5_read_dataset_p_2d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:), intent(inout) :: data
    end subroutine h5_read_dataset_p_2d_i


    subroutine h5_read_dataset_p_2d_d(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(:,:), intent(inout) :: data
    end subroutine h5_read_dataset_p_2d_d


    subroutine h5_read_dataset_p_2d_c(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in) :: dataset_name
        character(len=*), dimension(:,:), intent(inout) :: data
    end subroutine h5_read_dataset_p_2d_c


    subroutine h5_read_attribute_p(this, attribute_name, dataset_name, data)
        type(h5io), intent(in)               :: this
        character(len=*), intent(in)         :: attribute_name
        character(len=*), intent(in)         :: dataset_name
        integer, dimension(:), intent(inout) :: data
    end subroutine h5_read_attribute_p

    !
    ! parallel write routines
    !

    !
    ! independent writers
    !


    subroutine h5_create_dataset_setter(this, dataset_name, dim_in, rank, dtype_id)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:), intent(in) :: dim_in
        integer, intent(in)               :: rank
        integer, intent(in)               :: dtype_id ! 1: integer 4: real4 8:real8
    end subroutine h5_create_dataset_setter


    subroutine h5_open_dataset_p(this, dataset_name)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
    end subroutine h5_open_dataset_p


    ! dataset writer for 1d integer array
    subroutine h5_write_dataset_p_1d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:), intent(in) :: data
    end subroutine h5_write_dataset_p_1d_i


    subroutine h5_write_dataset_p_1d_r(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
    end subroutine h5_write_dataset_p_1d_r


    subroutine h5_write_dataset_p_2d_i(this, dataset_name, data)
        type(h5io), intent(in) :: this
        character(len=*), intent(in)      :: dataset_name
        integer, dimension(:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_2d_i


    subroutine h5_write_dataset_p_2d_r(this, dataset_name, data)
        type(h5io), intent(in)                             :: this
        character(len=*), intent(in)                       :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_2d_r


    subroutine h5_write_dataset_p_3d_i(this, dataset_name, data)
        type(h5io), intent(in)                               :: this
        character(len=*), intent(in)                         :: dataset_name
        integer, dimension(:,:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_3d_i


    subroutine h5_write_dataset_p_3d_r(this, dataset_name, data)
        type(h5io), intent(in)                :: this
        character(len=*), intent(in)          :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_3d_r


    subroutine h5_write_dataset_p_4d_i(this, dataset_name, data)
        type(h5io), intent(in)                  :: this
        character(len=*), intent(in)            :: dataset_name
        integer, dimension(:,:,:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_4d_i


    subroutine h5_write_dataset_p_4d_r(this, dataset_name, data)
        type(h5io), intent(in)                                 :: this
        character(len=*), intent(in)                           :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_4d_r


    subroutine h5_write_dataset_p_5d_r(this, dataset_name, data)
        type(h5io), intent(in)                                   :: this
        character(len=*), intent(in)                             :: dataset_name
        real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in) :: data
    end subroutine h5_write_dataset_p_5d_r


    subroutine h5_write_dataset_p_1d_l(this, dataset_name, data)
        type(h5io), intent(in)             :: this
        character(len=*), intent(in)       :: dataset_name
        logical, dimension(:), intent(in)  :: data
    end subroutine h5_write_dataset_p_1d_l


    !
    ! collective writers
    !


    !
    ! other utilities
    !
    subroutine bool_array2integer(this, boolarray, intarray)
        type(h5io), intent(in) :: this
        logical, dimension(:), intent(in) :: boolarray
        integer, dimension(:), intent(out) :: intarray
    end subroutine bool_array2integer


    subroutine int_array2bool(this,intarray, boolarray)
        type(h5io), intent(in) :: this
        integer, dimension(:), intent(in) :: intarray
        logical, dimension(:), intent(out) :: boolarray
    end subroutine int_array2bool


    subroutine h5_gather_dsetsize(this, dimin, data_rank, all_dim)
        type(h5io), intent(in) :: this
        integer, intent(in) :: data_rank
        integer, dimension(:), intent(in)  :: dimin
        integer, dimension(:,:) :: all_dim
    end subroutine


    subroutine h5_create_dataset_in_main_proc(this, dataset_name, dim, data_rank, dtype_id)
        type(h5io), intent(in)             :: this
        character(len=*), intent(in)       :: dataset_name
        integer, dimension(:), intent(in)  :: dim
        integer, intent(in)                :: dtype_id ! 1:int, 4:real4, 8:real8,
        integer, intent(in)                :: data_rank
    end subroutine

!
! higher level utilities
!
    subroutine write_attenuation_file_in_h5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)
        real, allocatable, dimension(:,:,:,:,:) :: factor_common
        real, allocatable, dimension(:,:,:,:)   :: scale_factor
        real, allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
        real, allocatable, dimension(:,:,:,:)   :: scale_factor_kappa
    end subroutine write_attenuation_file_in_h5


    subroutine read_attenuation_file_in_h5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)
        real, allocatable, dimension(:,:,:,:,:) :: factor_common
        real, allocatable, dimension(:,:,:,:)   :: scale_factor
        real, allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
        real, allocatable, dimension(:,:,:,:)   :: scale_factor_kappa
    end subroutine read_attenuation_file_in_h5


end module phdf5_utils