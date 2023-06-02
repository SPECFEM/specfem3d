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

!-------------------------------------------------------------------------------
!> Tools for parallel HDF5 support,
!> contains wrapper subroutines for common hdf5 calls
!
! The idea for these HDF5 tools is to convert all database files to a HDF5-format
! for better scalability to very large simulations.
! It is complementary to ADIOS2 support in case those libraries are not installed/wanted.
! These HDF5 files should avoid bottlenecks on HPC clusters due to the possible higher number of file outputs
! in our simulations (database files, seismograms, visualization snapshots, etc.).

! note: HDF5 library calls use a format like "h5**do_something***_f()"
!       our own wrapper functions in this module will rather use something like "h5_***do_something_***()",
!       and higher-level routines called in other Fortran routines like "***do_something_h5***()",
!       to better distinguish between library functions and wrappers.
!
! version info:
! - initial version, April 2023:
!   Masaru Nagaso (main developer, all initial hdf5 functionality)
!   Daniel Peter (just adding Masaru-san's original implementation, sometimes renaming things to more specfem-coherent ways)
!
!-------------------------------------------------------------------------------

!
! parallel HDF5 file I/O routines
!

module manager_hdf5

! class-like module for HDF5 routines

  use constants, only: CUSTOM_REAL, MAX_STRING_LEN

#if defined(USE_HDF5)
  use hdf5
#endif

  implicit none

  private

  ! public routines
  ! also to act as empty stubs (for missing HDF5 compilation support)
  public :: write_attenuation_file_hdf5, read_attenuation_file_hdf5
  public :: write_checkmesh_data_hdf5, write_checkmesh_xdmf_hdf5

  ! functional interface
  public :: h5_initialize

#if defined(USE_HDF5)
  ! only w/ HDF5 compilation

  public :: &
      h5_finalize

  public :: &
      h5_create_file, &
      h5_open_file, &
      h5_close_file, &
      h5_create_group, &
      h5_open_group, &
      h5_close_group, &
      h5_create_subgroup, &
      h5_open_subgroup, &
      h5_open_or_create_group, &
      h5_close_subgroup

  public :: &
      h5_open_dataset, &
      h5_open_dataset2, &
      h5_close_dataset

  public :: &
      h5_add_attribute_i, &
      h5_set_mpi_info, &
      h5_set_buffer_size, &
      h5_set_sieve_buffer_size, &
      h5_set_group_name, &
      h5_gather_dsetsize

  public :: &
      h5_create_file_p, &
      h5_create_file_p_collect, &
      h5_open_file_p, &
      h5_open_file_p_collect, &
      h5_close_file_p, &
      h5_create_file_prop_list, &
      h5_close_prop_list, &
      h5_close_prop_list_nocheck, &
      h5_open_group_prop_list, &
      h5_create_group_prop_list, &
      h5_close_group_prop_list, &
      h5_create_dataset_prop_list, &
      h5_create_dataset_gen, &
      h5_create_dataset_gen_in_group, &
      h5_check_dataset_exists, &
      h5_create_group_p, &
      h5_open_group_p

  public :: &
      h5_read_attribute_p, &
      h5_read_dataset_p, &
      h5_read_dataset_p_scalar, &
      h5_read_dataset_scalar_collect_hyperslab, &
      h5_read_dataset_collect_hyperslab_in_group, &
      h5_read_dataset_collect_hyperslab

  public :: &
      i2c

  public :: &
      h5_write_dataset, &
      h5_write_dataset_no_group, &
      h5_write_dataset_p, &
      h5_write_dataset_p_a, &
      h5_write_dataset_1d_to_2d_r_collect_hyperslab, &
      h5_write_dataset_2d_to_3d_r_collect_hyperslab, &
      h5_write_dataset_collect_hyperslab_in_group, &
      h5_write_dataset_collect_hyperslab

  ! private routines

  private :: &
    h5_check_arr_dim, &
    h5_check_collective, &
    h5_check_group

  private :: &
    create_dataset_collect, &
    bool_array2integer, &
    int_array2bool

  ! generic interface to read dataset
  interface h5_read_dataset_p
    module procedure h5_read_dataset_p_1d_l ! logical
    module procedure h5_read_dataset_p_1d_i ! integer
    module procedure h5_read_dataset_p_2d_i
    module procedure h5_read_dataset_p_3d_i
    module procedure h5_read_dataset_p_4d_i
    module procedure h5_read_dataset_p_1d_r ! real
    module procedure h5_read_dataset_p_2d_r
    module procedure h5_read_dataset_p_3d_r
    module procedure h5_read_dataset_p_4d_r
    module procedure h5_read_dataset_p_5d_r
    module procedure h5_read_dataset_p_2d_d ! double
    module procedure h5_read_dataset_p_2d_c ! char
  end interface h5_read_dataset_p

  ! generic interface to read scalar dataset
  interface h5_read_dataset_p_scalar
    module procedure h5_read_dataset_p_scalar_i
    module procedure h5_read_dataset_p_scalar_r
  end interface h5_read_dataset_p_scalar

  ! generic interface to read dataset in collective mode
  interface h5_read_dataset_scalar_collect_hyperslab
    module procedure h5_read_dataset_scalar_i_collect_hyperslab
    module procedure h5_read_dataset_scalar_r_collect_hyperslab
  end interface h5_read_dataset_scalar_collect_hyperslab

  ! generic interface to read dataset in collective mode
  interface h5_read_dataset_collect_hyperslab_in_group
    module procedure h5_read_dataset_1d_r_collect_hyperslab_in_group
    module procedure h5_read_dataset_2d_i_collect_hyperslab_in_group
  end interface h5_read_dataset_collect_hyperslab_in_group

  ! generic interface to read dataset in collective mode
  interface h5_read_dataset_collect_hyperslab
    module procedure h5_read_dataset_1d_l_collect_hyperslab ! logical
    module procedure h5_read_dataset_1d_i_collect_hyperslab ! integer
    module procedure h5_read_dataset_2d_i_collect_hyperslab
    module procedure h5_read_dataset_3d_i_collect_hyperslab
    module procedure h5_read_dataset_4d_i_collect_hyperslab
    module procedure h5_read_dataset_1d_r_collect_hyperslab ! real
    module procedure h5_read_dataset_2d_r_collect_hyperslab
    module procedure h5_read_dataset_3d_r_collect_hyperslab
    module procedure h5_read_dataset_4d_r_collect_hyperslab
    module procedure h5_read_dataset_5d_r_collect_hyperslab
    module procedure h5_read_dataset_2d_d_collect_hyperslab ! double
  end interface h5_read_dataset_collect_hyperslab

  ! generic interface to write dataset
  interface h5_write_dataset
    module procedure h5_write_dataset_1d_i ! integer
    module procedure h5_write_dataset_2d_i
    module procedure h5_write_dataset_1d_d ! double
    module procedure h5_write_dataset_2d_d
    module procedure h5_write_dataset_1d_c ! char
    module procedure h5_write_dataset_2d_c
    module procedure h5_write_dataset_2d_r ! real
    module procedure h5_write_dataset_4d_r
  end interface h5_write_dataset

  interface h5_write_dataset_no_group
    module procedure h5_write_dataset_1d_i_no_group
    module procedure h5_write_dataset_1d_d_no_group
    module procedure h5_write_dataset_1d_c_no_group
    module procedure h5_write_dataset_2d_r_no_group
  end interface h5_write_dataset_no_group

  ! generic interface to write dataset
  interface h5_write_dataset_p
    module procedure h5_write_dataset_p_1d_l  ! logical
    module procedure h5_write_dataset_p_1d_i  ! integer
    module procedure h5_write_dataset_p_2d_i
    module procedure h5_write_dataset_p_3d_i
    module procedure h5_write_dataset_p_4d_i
    module procedure h5_write_dataset_p_1d_r  ! real
    module procedure h5_write_dataset_p_2d_r
    module procedure h5_write_dataset_p_3d_r
    module procedure h5_write_dataset_p_4d_r
    module procedure h5_write_dataset_p_5d_r
    module procedure h5_write_dataset_p_2d_d  ! double
  end interface h5_write_dataset_p

  interface h5_write_dataset_p_a
    module procedure h5_write_dataset_p_1d_ia
    module procedure h5_write_dataset_p_2d_ia
  end interface h5_write_dataset_p_a

  ! generic interface to write dataset in collective mode
  interface h5_write_dataset_collect_hyperslab_in_group
    module procedure h5_write_dataset_2d_i_collect_hyperslab_in_group
    module procedure h5_write_dataset_1d_r_collect_hyperslab_in_group
  end interface h5_write_dataset_collect_hyperslab_in_group

  ! generic interface to write dataset in collective mode
  interface h5_write_dataset_collect_hyperslab
    module procedure h5_write_dataset_1d_l_collect_hyperslab ! logical
    module procedure h5_write_dataset_1d_i_collect_hyperslab ! integer
    module procedure h5_write_dataset_2d_i_collect_hyperslab
    module procedure h5_write_dataset_3d_i_collect_hyperslab
    module procedure h5_write_dataset_4d_i_collect_hyperslab
    module procedure h5_write_dataset_1d_r_collect_hyperslab ! real
    module procedure h5_write_dataset_2d_r_collect_hyperslab
    module procedure h5_write_dataset_3d_r_collect_hyperslab
    module procedure h5_write_dataset_4d_r_collect_hyperslab
    module procedure h5_write_dataset_5d_r_collect_hyperslab
    module procedure h5_write_dataset_2d_d_collect_hyperslab ! double
  end interface h5_write_dataset_collect_hyperslab

  ! object-oriented interface
  ! (Fortran 2003 standard style)
  !
  ! usage example:
  !   subroutine **
  !     ..
  !     use manager_hdf5, only: h5io
  !     type(h5io) :: h5
  !     ..
  !     h5 = h5io()               ! calls the object constructor, i.e. h5io_constructor()
  !     call h5%open(fname)       ! object function call
  !     ! or in a single call:
  !     !   h5 = h5io(fname)
  !     call h5%write("x",store_val_x_all)
  !     call h5%close()
  !     ..
  !   end subroutine              ! h5 goes out of scope, will call the object destructor, i.e., h5io_destructor()
  !
  type, public :: h5io
    private
    logical :: is_initialized = .false.
    integer(HID_T) :: f_id = -1, g_id = -1, d_id = -1
    character(len=MAX_STRING_LEN) :: filename = ""
  contains
    procedure :: open => h5io_open_file
    procedure :: close => h5io_close_file
    generic   :: write => h5io_write_dataset_i, h5io_write_dataset_r
    procedure, private :: h5io_write_dataset_i
    procedure, private :: h5io_write_dataset_r
    final :: h5io_destructor
  end type h5io
  ! object constructor
  ! called by: h5 = h5io()
  interface h5io
    module procedure :: h5io_constructor
    module procedure :: h5io_constructor_from_file
  end interface h5io

  ! module parameters
  ! ids
  integer(HID_T) :: file_id, group_id, parent_group_id, dataset_id
  integer(HID_T) :: mem_dspace_id, file_dspace_id                     ! for collective IO

  ! parallel process
  integer(HID_T) :: plist_id, fplist_id, gplist_id

  ! string array
  integer(SIZE_T) :: str_len = 16
  integer(HID_T) :: str_type

  type(c_ptr) :: f_ptr

  ! MPI info
  integer :: this_rank, this_info, this_comm, total_proc

  ! io unit
  integer, parameter :: xdmf_mesh = 190

  ! mesh nodes
  public :: xdmf_mesh_nnodes
  integer :: xdmf_mesh_nnodes

  ! function call errors
  integer :: error

  ! class-wide private variables
  character(len=256) :: file_path
  character(len=256) :: store_group_name ! groupname for debug

  public :: name_database_hdf5

  ! store HDF5 output file name
  character(len=MAX_STRING_LEN) :: name_database_hdf5

#endif

contains

!-------------------------------------------------------------------------------
!
! public HDF5 wrapper routines (also available without hdf5 compilation support)
!
!-------------------------------------------------------------------------------

  subroutine h5_initialize()

#if defined(USE_HDF5)
    use constants, only: MAX_LENGTH_NETWORK_NAME, MAX_LENGTH_STATION_NAME
#endif

    implicit none

#if defined(USE_HDF5)
    ! initialize Fortran interface
    call h5open_f(error)
    call check_error()

    ! prepare string array type
    if (MAX_LENGTH_STATION_NAME >= MAX_LENGTH_NETWORK_NAME) then
      str_len = MAX_LENGTH_STATION_NAME
    else
      str_len = MAX_LENGTH_NETWORK_NAME
    endif

    call h5tcopy_f(H5T_Fortran_S1, str_type, error)
    call check_error()

    call h5tset_size_f(str_type, str_len, error)
    call check_error()
#else
  ! no HDF5 compilation support

  ! compilation without HDF5 support
  print *
  print *, "Error: HDF5 routine h5_initialize() called without HDF5 Support."
  print *, "To enable HDF5 support, reconfigure with --with-hdf5 flag."
  print *

  ! safety stop
  stop 'Error HDF5 manager: h5_initialize() intitialization called without compilation support'

#endif

  end subroutine h5_initialize

!-------------------------------------------------------------------------------
!
! higher level utilities
!
!-------------------------------------------------------------------------------

  subroutine write_attenuation_file_hdf5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)

#if defined(USE_HDF5)
    use shared_parameters, only: NPROC, LOCAL_PATH
    use constants, only: myrank,N_SLS,NGLLX,NGLLY,NGLLZ
#endif

    implicit none

    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor_kappa

#if defined(USE_HDF5)
    integer, dimension(4) :: dims
    integer :: nspec

    ! offset arrays
    integer, dimension(0:NPROC-1) :: offset_nspec

    ! hdf5 valiables
    character(len=64) :: filename, dset_name, tempstr

    dims = shape(scale_factor)
    nspec = dims(4)

    ! prepare offset arrays
    call gather_all_all_singlei(nspec,offset_nspec,NPROC) ! n spec in each proc

    tempstr = "/external_mesh.h5"
    filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

    ! initialize h5 object
    call h5_initialize()

    ! prepare dataset
    if (myrank == 0) then
      call h5_open_file(filename)
      dset_name = "scale_factor"
      call h5_create_dataset_gen(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 4, CUSTOM_REAL)
      dset_name = "scale_factor_kappa"
      call h5_create_dataset_gen(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 4, CUSTOM_REAL)
      dset_name = "factor_common"
      call h5_create_dataset_gen(dset_name, (/N_SLS,NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 5, CUSTOM_REAL)
      dset_name = "factor_common_kappa"
      call h5_create_dataset_gen(dset_name, (/N_SLS,NGLLX,NGLLY,NGLLZ,sum(offset_nspec)/), 5, CUSTOM_REAL)
      call h5_close_file()
    endif

    call synchronize_all()

    ! write
    call h5_open_file_p_collect(filename)
    dset_name = "scale_factor"
    call h5_write_dataset_4d_r_collect_hyperslab(dset_name, &
                                                 scale_factor,(/0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
    dset_name = "scale_factor_kappa"
    call h5_write_dataset_4d_r_collect_hyperslab(dset_name, &
                                                 scale_factor_kappa,(/0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
    dset_name = "factor_common"
    call h5_write_dataset_5d_r_collect_hyperslab(dset_name, &
                                                 factor_common,(/0,0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)
    dset_name = "factor_common_kappa"
    call h5_write_dataset_5d_r_collect_hyperslab(dset_name, &
                                                 factor_common_kappa,(/0,0,0,0,sum(offset_nspec(0:myrank-1))/),.true.)

    call h5_close_file_p()
    call h5_finalize()
#else
  ! no compilation support
  ! to avoid compiler warnings
  real(kind=CUSTOM_REAL) :: dummy

  dummy = factor_common(1,1,1,1,1)
  dummy = factor_common_kappa(1,1,1,1,1)
  dummy = scale_factor(1,1,1,1)
  dummy = scale_factor_kappa(1,1,1,1)

  ! safety stop
  stop 'Error HDF5 manager: write_attenuation_file_hdf5() called without compilation support'
#endif

  end subroutine write_attenuation_file_hdf5

!
!-------------------------------------------------------------------------------
!

  subroutine read_attenuation_file_hdf5(factor_common, scale_factor, factor_common_kappa, scale_factor_kappa)

#if defined(USE_HDF5)
    use shared_parameters, only: NPROC, LOCAL_PATH
    use constants, only: myrank
#endif

    implicit none

    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: factor_common_kappa
    real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:)   :: scale_factor_kappa

#if defined(USE_HDF5)
    ! offset array
    integer, dimension(0:NPROC-1) :: offset_nspec

    ! hdf5 valiables
    character(len=64) :: fname, tempstr

    tempstr = "/external_mesh.h5"
    fname = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

    ! initialize h5 object
    call h5_initialize()

    ! open file
    call h5_open_file_p_collect(fname)
    ! read offset array
    call h5_read_dataset_1d_i_collect_hyperslab("offset_nspec",offset_nspec, (/0/), .true.)

    call h5_read_dataset_4d_r_collect_hyperslab("scale_factor", scale_factor, &
                                                (/0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
    call h5_read_dataset_4d_r_collect_hyperslab("scale_factor_kappa", scale_factor_kappa, &
                                                (/0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
    call h5_read_dataset_5d_r_collect_hyperslab("factor_common", factor_common, &
                                                (/0,0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)
    call h5_read_dataset_5d_r_collect_hyperslab("factor_common_kappa", factor_common_kappa, &
                                                (/0,0,0,0,sum(offset_nspec(0:myrank-1))/), .true.)

    call h5_close_file_p()
    call h5_finalize()
#else
  ! no compilation support
  ! to avoid compiler warnings
  real(kind=CUSTOM_REAL) :: dummy

  dummy = factor_common(1,1,1,1,1)
  dummy = factor_common_kappa(1,1,1,1,1)
  dummy = scale_factor(1,1,1,1)
  dummy = scale_factor_kappa(1,1,1,1)

  ! safety stop
  stop 'Error HDF5 manager: read_attenuation_file_hdf5() called without compilation support'
#endif

  end subroutine read_attenuation_file_hdf5

!
!-------------------------------------------------------------------------------
!

  subroutine write_checkmesh_data_hdf5(dset_name,dump_array)

#if defined(USE_HDF5)
    use shared_parameters, only: LOCAL_PATH, NPROC
    use constants, only: myrank
#endif

    implicit none

    real(kind=CUSTOM_REAL),dimension(:), intent(in) :: dump_array
    character(len=MAX_STRING_LEN), intent(in)       :: dset_name

#if defined(USE_HDF5)
    character(len=MAX_STRING_LEN)                   :: filename
    integer, dimension(0:NPROC-1)                   :: offset

    ! MPI variables
    integer :: info, comm

    ! hdf5 valiables
    character(len=64) :: tempstr

    ! flag if dataset exists
    logical           :: exists = .false.

    ! saves mesh file external_mesh.h5
    tempstr = "/external_mesh.h5"
    filename = LOCAL_PATH(1:len_trim(LOCAL_PATH))//trim(tempstr)

    ! get MPI parameters
    call world_get_comm(comm)
    call world_get_info_null(info)

    ! initialize h5 object
    call h5_initialize()

    call h5_set_mpi_info(comm, info, myrank, NPROC)

    ! get offset info
    call gather_all_all_singlei(size(dump_array), offset, NPROC)

    ! make dataset
    if (myrank == 0) then
      call h5_open_file(filename)
      ! check if dataset exists
      call h5_check_dataset_exists(dset_name, exists)
      if (.not. exists) then
        call h5_create_dataset_gen(dset_name, (/sum(offset(:))/), 1, CUSTOM_REAL)
      endif
      call h5_close_file()
    endif
    call synchronize_all()

    ! open file
    call h5_open_file_p_collect(filename)
    call h5_write_dataset_1d_r_collect_hyperslab(dset_name, dump_array, (/sum(offset(0:myrank-1))/),.true.)
    call h5_close_file_p()

    call h5_finalize()
#else
  ! no compilation support
  ! to avoid compiler warnings
  real(kind=CUSTOM_REAL) :: dummy
  character(len=1) :: c_dummy

  dummy = dump_array(1)
  c_dummy = dset_name(1:1)

  ! safety stop
  stop 'Error HDF5 manager: write_checkmesh_data_hdf5() called without compilation support'
#endif

  end subroutine write_checkmesh_data_hdf5

!
!-------------------------------------------------------------------------------
!

  subroutine write_checkmesh_xdmf_hdf5(NSPEC_AB)

#if defined(USE_HDF5)
    use constants, only: myrank
    use shared_parameters
#endif
    implicit none

    integer, intent(in)                             :: NSPEC_AB

#if defined(USE_HDF5)
    character(len=MAX_STRING_LEN)                   :: fname_xdmf_checkmesh,fname_h5_database_xdmf,fname_h5_extmesh_xdmf
    integer, dimension(0:NPROC-1)                   :: nelms, nnodes
    character(len=20)                               :: type_str,nelm_str,nnode_str

    ! gather number of elements in each proc
    call gather_all_singlei(NSPEC_AB, nelms, NPROC)

    ! count and gather the number of control nodes in each proc
    call gather_all_singlei(xdmf_mesh_nnodes, nnodes, NPROC)

    if (myrank == 0) then
      ! writeout xdmf file for surface movie
      fname_xdmf_checkmesh = trim(LOCAL_PATH)//"/checkmesh.xmf"
      fname_h5_database_xdmf = "./Database.h5"      ! relative to checkmesh.xmf file
      fname_h5_extmesh_xdmf = "./external_mesh.h5"

      open(unit=xdmf_mesh, file=trim(fname_xdmf_checkmesh), recl=512)

      ! definition of topology and geometry
      ! refer only control nodes (8 or 27) as a coarse output
      ! data array need to be extracted from full data array on GLL points
      write(xdmf_mesh,'(a)') '<?xml version="1.0" ?>'
      write(xdmf_mesh,*) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(xdmf_mesh,*) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="3.0">'
      write(xdmf_mesh,*) '<Domain>'
      write(xdmf_mesh,*) '<!-- mesh info -->'
      !write(xdmf_mesh,*) '<Grid Name="mesh" GridType="Collection"  CollectionType="Spatial">'

      nelm_str  = i2c(sum(nelms(:)))
      nnode_str = i2c(sum(nnodes(:)))

      write(xdmf_mesh,*) '<Grid Name="mesh">'
      write(xdmf_mesh,*) '<Topology TopologyType="Mixed" NumberOfElements="'//trim(nelm_str)//'">'
      write(xdmf_mesh,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Int" Precision="4" Dimensions="'&
                             //trim(nelm_str)//' '//trim(i2c(8+1))//'">'
      write(xdmf_mesh,*) '       '//trim(fname_h5_database_xdmf)//':/elm_conn_xdmf'
      write(xdmf_mesh,*) '</DataItem>'
      write(xdmf_mesh,*) '</Topology>'
      write(xdmf_mesh,*) '<Geometry GeometryType="XYZ">'
      write(xdmf_mesh,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                          //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nnode_str)//' 3">'
      write(xdmf_mesh,*) '       '//trim(fname_h5_database_xdmf)//':/nodes_coords'
      write(xdmf_mesh,*) '</DataItem>'
      write(xdmf_mesh,*) '</Geometry>'

      type_str = "res_Courant_number"
      write(xdmf_mesh,*) '<Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Cell">'
      write(xdmf_mesh,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nelm_str)//'">'
      write(xdmf_mesh,*) '       '//trim(fname_h5_extmesh_xdmf)//':/'//trim(type_str)
      write(xdmf_mesh,*) '</DataItem>'
      write(xdmf_mesh,*) '</Attribute>'

      type_str = "res_minimum_period"
      write(xdmf_mesh,*) '<Attribute Name="'//trim(type_str)//'" AttributeType="Scalar" Center="Cell">'
      write(xdmf_mesh,*) '<DataItem ItemType="Uniform" Format="HDF" NumberType="Float" Precision="'&
                                     //trim(i2c(CUSTOM_REAL))//'" Dimensions="'//trim(nelm_str)//'">'
      write(xdmf_mesh,*) '       '//trim(fname_h5_extmesh_xdmf)//':/'//trim(type_str)
      write(xdmf_mesh,*) '</DataItem>'
      write(xdmf_mesh,*) '</Attribute>'

      !write(xdmf_mesh,*) '</Grid>'
      write(xdmf_mesh,*) '</Grid>'
      write(xdmf_mesh,*) '</Domain>'
      write(xdmf_mesh,*) '</Xdmf>'

      close(xdmf_mesh)
    endif

#else
  ! no compilation support
  ! to avoid compiler warnings
  integer :: dummy

  dummy = NSPEC_AB

  ! safety stop
  stop 'Error HDF5 manager: write_checkmesh_xdmf_hdf5() called without compilation support'
#endif

  end subroutine write_checkmesh_xdmf_hdf5


!-------------------------------------------------------------------------------
!
! HDF5 wrapper routines (only available with HDF5 compilation support)
!
!-------------------------------------------------------------------------------

#if defined(USE_HDF5)
! only available with HDF5 compilation support


  function i2c(k) result(str)
  ! "Convert an integer to string."
    implicit none
    integer, intent(in) :: k
    character(len=20) str
    write (str, "(i20)") k
    str = adjustl(str)
  end function i2c

!
!-------------------------------------------------------------------------------
!

  subroutine h5_finalize()
    implicit none
    ! close Fortran interface
    call h5close_f(error)
    call check_error()
  end subroutine h5_finalize

!
!-------------------------------------------------------------------------------
!

  subroutine h5_check_collective(dataset_name)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer :: io_mode
    call h5pget_mpio_actual_io_mode_f(plist_id,io_mode,error)
    if (error /= 0) write(*,*) 'hdf5 get_mpio_actual_io_mode failed for ',trim(dataset_name)
    !if (io_mode == H5D_MPIO_NO_COLLECTIVE_F) print *, &
    !    "collective read/write not possible for dataset: ", dataset_name
    call check_error()
  end subroutine h5_check_collective

!
!-------------------------------------------------------------------------------
!

  subroutine h5_check_arr_dim(dim)
    implicit none
    integer(kind=HSIZE_T), dimension(:), intent(in) :: dim
    integer :: i
    ! if one of the dimension is 0, cancel space_id and file_id
    ! loop all elements in dim
    do i = 1, size(dim)
      if (dim(i) == 0) then
        call h5sselect_none_f(mem_dspace_id, error)
        call check_error()
        call h5sselect_none_f(file_dspace_id, error)
        call check_error()
      endif
    enddo
  end subroutine h5_check_arr_dim

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_file(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, error)
    if (error /= 0) then
      print *,'Error file create: ', trim(file_path)
      print *
      print *,'check if path exists: ', trim(file_path)
      stop 'Error file create h5'
    endif
  end subroutine h5_create_file

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_file(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error)
    if (error /= 0) then
      print *, 'Error file open: ', trim(file_path)
      stop 'Error file open h5'
    endif
  end subroutine h5_open_file

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_file()
    implicit none
    call h5fclose_f(file_id, error)
    if (error /= 0) write(*,*) 'error while closing a file.'
    call check_error()
  end subroutine h5_close_file

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_file_p()
    implicit none
    call h5_close_file_prop_list()
    call h5fclose_f(file_id, error)
    if (error /= 0) write(*,*) 'error while closing a file.'
    call check_error()
  end subroutine h5_close_file_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_group(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    call h5gcreate_f(file_id, trim(group_name), group_id, error)
    if (error /= 0) write(*,*) 'error while creating a group, ', group_name
    call check_error()
    call h5gclose_f(group_id, error)
    if (error /= 0) write(*,*) 'error while closing a group, ' , group_name
    call check_error()
  end subroutine h5_create_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_subgroup(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    integer(HID_T) :: temp_group_id

    call h5gcreate_f(group_id, trim(group_name), temp_group_id, error)
    if (error /= 0) write(*,*) 'error while creating a subgroup, ', group_name
    call check_error()
    call h5gclose_f(temp_group_id, error)
    if (error /= 0) write(*,*) 'error while closing a subgroup, ' , group_name
    call check_error()
  end subroutine h5_create_subgroup

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_group(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    ! open group
    call h5gopen_f(file_id, trim(group_name), group_id, error) ! group open
    if (error /= 0) write(*,*) 'hdf5 open group failed for ', trim(group_name)
    call check_error()
    ! stores name
    store_group_name = group_name
  end subroutine h5_open_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_subgroup(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    integer(HID_T) :: temp_group_id
    ! open subgroup
    call h5gopen_f(group_id, trim(group_name), temp_group_id, error) ! group open
    if (error /= 0) write(*,*) 'hdf5 open subgroup failed for ', trim(group_name)
    call check_error()
    ! put the group id of the first level becomes parent group
    ! only while the group at second level exists
    parent_group_id = group_id
    group_id        = temp_group_id
  end subroutine h5_open_subgroup

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_or_create_group(group_name)
    ! Check if group with the given name exists. Create it if it doesn't,
    ! open it if it does.
    implicit none
    character(len=*), intent(in) :: group_name
    ! Variable for checking if a group exists or not
    logical :: group_exists
    ! check group
    call h5_check_group(group_name, group_exists)
    ! open or create
    if (group_exists) then
      call h5gopen_f(file_id, group_name, group_id, error)
    else
      call h5gcreate_f(file_id, group_name, group_id, error)
    endif
    call check_error()
  end subroutine h5_open_or_create_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_check_group(group_name,group_exists)
    ! Check if group with the given name exists. Create it if it doesn't,
    ! open it if it does.
    implicit none
    character(len=*), intent(in) :: group_name
    logical, intent(out) :: group_exists
    ! Variable for checking if a group exists or not
    group_exists = .false.
    ! check
    call h5lexists_f(file_id, group_name, group_exists, error)
    call check_error()
  end subroutine h5_check_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_group()
    implicit none
    call h5gclose_f(group_id, error) ! group open
    if (error /= 0) write(*,*) 'hdf5 close group failed'
    call check_error()
  end subroutine h5_close_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_subgroup()
    implicit none
    call h5gclose_f(group_id, error) ! group open
    if (error /= 0) write(*,*) 'hdf5 close group failed'
    call check_error()
    ! stores group
    group_id = parent_group_id
  end subroutine h5_close_subgroup

!
!-------------------------------------------------------------------------------
!

  subroutine h5_set_group_name(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    store_group_name = group_name
  end subroutine h5_set_group_name

!
!-------------------------------------------------------------------------------
!

  ! open dataset. group need to be opened before.
  subroutine h5_open_dataset(dataset_name)
    implicit none
    character(len=*), intent(in) :: dataset_name

    call h5dopen_f(group_id, trim(dataset_name), dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 open dataset failed for ', dataset_name
    call check_error()
  end subroutine h5_open_dataset

!
!-------------------------------------------------------------------------------
!

  ! open dataset without open group
  subroutine h5_open_dataset2(dataset_name)
    implicit none
    character(len=*), intent(in) :: dataset_name

    call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 open dataset2 failed for ', dataset_name
    call check_error()
  end subroutine h5_open_dataset2

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_dataset()
    implicit none
    call h5dclose_f(dataset_id, error) ! group open
    if (error /= 0) write(*,*) 'hdf5 close dataset failed'
    call check_error()
  end subroutine h5_close_dataset

!
!-------------------------------------------------------------------------------
!

  subroutine h5_check_dataset_exists(dataset_name, exists)
    implicit none
    character(len=*), intent(in) :: dataset_name
    logical, intent(out) :: exists

    call h5lexists_f(file_id, dataset_name, exists, error)
    call check_error()
  end subroutine h5_check_dataset_exists


!-------------------------------------------------------------------------------
!
! serial write routines
!
!-------------------------------------------------------------------------------

  ! dataset writer for 1d integer array
  subroutine h5_write_dataset_1d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    integer, dimension(:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 1
    integer(HSIZE_T), dimension(1)    :: dim

    dim = size(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_1d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_1d_i_no_group(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    integer, dimension(:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 1
    integer(HSIZE_T), dimension(1)    :: dim

    dim = size(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
    ! closes dataset
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_i_no_group

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 1d custom real array
  subroutine h5_write_dataset_1d_d(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 1
    integer(HSIZE_T), dimension(1)    :: dim

    dim = size(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    else
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    endif
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_1d_d

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 1d custom real array without grouping
  subroutine h5_write_dataset_1d_d_no_group(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 1
    integer(HSIZE_T), dimension(1)    :: dim

    dim = size(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    else
      call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    endif
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
    ! closes dataset
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_d_no_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_1d_c(dataset_name, data_in)
    implicit none
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
    call check_error()
    call h5dcreate_f(group_id, trim(dataset_name), str_type, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, str_type, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
    deallocate(data, stat=error)

  end subroutine h5_write_dataset_1d_c

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_1d_c_no_group(dataset_name, data_in)
    implicit none
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
    call check_error()
    call h5dcreate_f(file_id, trim(dataset_name), str_type, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, str_type, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
    ! closes dataset
    call h5_close_dataset()
    deallocate(data, stat=error)
  end subroutine h5_write_dataset_1d_c_no_group

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 2d integer array
  subroutine h5_write_dataset_2d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    integer, dimension(:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 2
    integer(HSIZE_T), dimension(2)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_INTEGER, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_2d_i

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 2d double array
  subroutine h5_write_dataset_2d_d(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    double precision, dimension(:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 2
    integer(HSIZE_T), dimension(2)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_2d_d

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 2d double array
  subroutine h5_write_dataset_2d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 2
    integer(HSIZE_T), dimension(2)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    else
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    endif
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_2d_r

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 2d double array
  subroutine h5_write_dataset_2d_r_no_group(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 2
    integer(HSIZE_T), dimension(2)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    else
      call h5dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    endif
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
    ! closes dataset
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_r_no_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_4d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 4
    integer(HSIZE_T), dimension(4)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_REAL, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    else
      call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
      call check_error()
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error)
      if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
      call check_error()
    endif
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_4d_r

!
!-------------------------------------------------------------------------------
!

  ! dataset writer for 2d character array
  subroutine h5_write_dataset_2d_c(dataset_name, data)
    implicit none
    character(len=*), intent(in)      :: dataset_name
    character(len=*), dimension(:,:), intent(in) :: data
    integer(HID_T)                    :: dspace_id ! dataspace id is local.
    integer                           :: rank = 2
    integer(HSIZE_T), dimension(2)    :: dim

    dim = shape(data)

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5dcreate_f(group_id, trim(dataset_name), H5T_NATIVE_CHARACTER, dspace_id, dataset_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()
    call h5dwrite_f(dataset_id, H5T_NATIVE_CHARACTER, data, dim, error)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', dataset_name
    call check_error()
  end subroutine h5_write_dataset_2d_c

!
!-------------------------------------------------------------------------------
!

  ! set attribute to a dataset
  subroutine h5_add_attribute_i(attribute_name, data)
    implicit none
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
    call check_error()
    call h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error) ! for a string tag of attribute value
    if (error /= 0) write(*,*) 'hdf5 tcopy failed for attribute, ', attribute_name
    call check_error()
    call h5tset_size_f(atype_id, taglen, error)
    if (error /= 0) write(*,*) 'hdf5 set_size failed for attribute, ', attribute_name
    call check_error()
    ! here the attribute is written on the current openning dataset
    call h5acreate_f(dataset_id, trim(attribute_name), atype_id, aspace_id, attr_id, error)
    if (error /= 0) write(*,*) 'hdf5 acreate failed for attribute, ', attribute_name
    call check_error()
    call h5awrite_f(attr_id, atype_id, data, dim, error) ! write
    if (error /= 0) write(*,*) 'hdf5 awrite failed for attribute, ', attribute_name
    call check_error()
    call h5sclose_f(aspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace closing failed for ', attribute_name
    call check_error()
    call h5aclose_f(attr_id,    error)
    if (error /= 0) write(*,*) 'hdf5 aclose failed for ', attribute_name
    call check_error()
    call h5tclose_f(atype_id,   error)
    if (error /= 0) write(*,*) 'hdf5 tclose failed for ', attribute_name
    call check_error()
  end subroutine h5_add_attribute_i


!-------------------------------------------------------------------------------
!
! parallel routines
!
!-------------------------------------------------------------------------------

  subroutine h5_set_mpi_info(comm, info, rank, nproc)
    implicit none
    integer, intent(in) :: comm, info, rank, nproc

    this_comm  = comm
    this_info  = info
    this_rank  = rank
    total_proc = nproc
  end subroutine h5_set_mpi_info

!
!-------------------------------------------------------------------------------
!

  subroutine h5_set_sieve_buffer_size()
    implicit none
    integer(hsize_t)       :: buf_size  = 1*1024*1024   ! 1 MB
    integer(hsize_t)       :: alig_size = 1*1024*1024   ! 1 MB
    call h5pset_sieve_buf_size_f(fplist_id, buf_size, error) ! buf_size may vary depending on machiens
    call check_error()
    call h5pset_alignment_f(fplist_id, buf_size, alig_size, error)
    call check_error()
  end subroutine h5_set_sieve_buffer_size

!
!-------------------------------------------------------------------------------
!

  subroutine h5_set_buffer_size()
    implicit none
    integer(hsize_t)       :: buf_size = 1*1024*1024    ! 1 MB
    call h5pset_buffer_f(plist_id, buf_size, error)
    call check_error()
  end subroutine h5_set_buffer_size

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_file_prop_list(if_collective)
    implicit none
    logical, intent(in) :: if_collective

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, fplist_id, error)
    if (error /= 0) write(*,*) 'hdf5 create plist failed.'
    call check_error()
    call h5pset_all_coll_metadata_ops_f(fplist_id, if_collective, error)
    call check_error()
    call h5pset_coll_metadata_write_f(fplist_id, if_collective, error)
    call check_error()
    call h5pset_fapl_mpio_f(fplist_id, this_comm, this_info, error)
    if (error /= 0) write(*,*) 'hdf5 set_fapl failed.'
    call check_error()
  end subroutine h5_create_file_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_group_prop_list()
    implicit none
    logical, parameter :: if_collective = .true.
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_GROUP_ACCESS_F, gplist_id, error)
    if (error /= 0) write(*,*) 'hdf5 create group plist failed.'
    call check_error()
    call h5pset_all_coll_metadata_ops_f(gplist_id, if_collective, error)
    call check_error()
    call h5pset_coll_metadata_write_f(gplist_id, if_collective, error)
    if (error /= 0) write(*,*) 'hdf5 group set_all_coll_metadata_ops failed.'
    call check_error()
  end subroutine h5_open_group_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_group_prop_list()
    implicit none
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_GROUP_CREATE_F, gplist_id, error)
    if (error /= 0) write(*,*) 'hdf5 create group plist failed.'
    call check_error()
  end subroutine h5_create_group_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_dataset_prop_list(if_collective)
    implicit none
    logical, intent(in) :: if_collective
    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error /= 0) write(*,*) 'hdf5 create dataset plist failed.'
    call check_error()
    ! use larger buffer size from the default
    call h5_set_buffer_size()
    if (if_collective) then
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    else
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset set_dxpl_mpio failed.'
    call check_error()
  end subroutine h5_create_dataset_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_prop_list(dataset_name)
    implicit none
    character(len=*), intent(in) :: dataset_name
    !print *,"DEBUG closing dataset prop: ", dataset_name
    ! check if the latest dataset transfer has been done in collective mode
    call h5_check_collective(dataset_name)
    call h5pclose_f(plist_id, error) ! property list can be closed soon
    if (error /= 0) write(*,*) 'hdf5 close_prop_list failed for ',trim(dataset_name)
    call check_error()
  end subroutine h5_close_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_prop_list_nocheck(dataset_name)
    implicit none
    character(len=*), intent(in) :: dataset_name
    !print *,"DEBUG closing dataset prop: ", dataset_name
    ! check if the latest dataset transfer has been done in collective mode
    call h5pclose_f(plist_id, error) ! property list can be closed soon
    if (error /= 0) write(*,*) 'hdf5 close_prop_list failed for ',trim(dataset_name)
    call check_error()
  end subroutine h5_close_prop_list_nocheck

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_file_prop_list()
    implicit none
    call h5pclose_f(fplist_id, error) ! property list can be closed soon
    if (error /= 0) write(*,*) 'hdf5 close_file_prop_list failed.'
    call check_error()
  end subroutine h5_close_file_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_close_group_prop_list()
    implicit none
    call h5pclose_f(gplist_id, error) ! property list can be closed soon
    if (error /= 0) write(*,*) 'hdf5 close_prop_list failed.'
    call check_error()
  end subroutine h5_close_group_prop_list

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_file_p(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5_create_file_prop_list(.false.)
    call h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, error, creation_prp=H5P_DEFAULT_F, access_prp=fplist_id)
    if (error /= 0) write(*,*) 'hdf5 create file p failed.'
    call check_error()
  end subroutine h5_create_file_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_file_p_collect(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5_create_file_prop_list(.true.)
    call h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, error, creation_prp=H5P_DEFAULT_F, access_prp=fplist_id)
    if (error /= 0) write(*,*) 'hdf5 create file p failed.'
    call check_error()
  end subroutine h5_create_file_p_collect

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_file_p(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5_create_file_prop_list(.false.)
    call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error, access_prp=fplist_id)
    if (error /= 0) write(*,*) 'hdf5 open file p failed.'
    call check_error()
  end subroutine h5_open_file_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_file_p_collect(fpath_in)
    implicit none
    character(len=*), intent(in) :: fpath_in

    ! set the target file path
    file_path = fpath_in

    call h5_create_file_prop_list(.true.)
    call h5_set_sieve_buffer_size()

    call h5fopen_f(trim(file_path), H5F_ACC_RDWR_F, file_id, error, access_prp=fplist_id)
    if (error /= 0) write(*,*) 'hdf5 open file p failed.'
    call check_error()
  end subroutine h5_open_file_p_collect

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_dataset_gen(dataset_name, dim_in, rank, dtype_id)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:), intent(in)  :: dim_in

    integer(HSIZE_T), dimension(size(dim_in)) :: dim
    integer, intent(in)                :: dtype_id ! 1:int, 4:real4, 8:real8,
    integer, intent(in)                :: rank

    integer(HID_T)                     :: dspace_id
    !logical :: if_chunk = .true.
    !integer :: i

    dim = dim_in ! convert data type

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call check_error()

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
      stop 'Invalid dtype_id, not implemented yet in h5_create_dataset_gen() routine'
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()

    call h5_close_prop_list_nocheck(dataset_name)

    call h5dclose_f(dataset_id,error)
    if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name
    call check_error()

  end subroutine h5_create_dataset_gen

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_dataset_gen_in_group(dataset_name, dim_in, rank, dtype_id)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:), intent(in)  :: dim_in

    integer(HSIZE_T), dimension(size(dim_in)) :: dim
    integer, intent(in)                :: dtype_id ! 1:int, 4:real4, 8:real8,
    integer, intent(in)                :: rank

    integer(HID_T)                     :: dspace_id

    dim = dim_in ! convert data type

    call h5screate_simple_f(rank, dim, dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
    call check_error()
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call check_error()
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
      print *, "Error: specified dtype_id is not implemented yet for hdf5 io. aborting..."
      stop 'Invalid dtype_id, not implemented yet in h5_create_dataset_gen_in_group()'
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset create failed for ', dataset_name
    call check_error()

    call h5_close_prop_list_nocheck(dataset_name)

    call h5dclose_f(dataset_id,error)
    if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
    call check_error()
    call h5sclose_f(dspace_id, error)
    if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name
    call check_error()
  end subroutine h5_create_dataset_gen_in_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_create_group_p(group_name)
    implicit none
    character(len=*), intent(in) :: group_name

    call h5_create_group_prop_list()
    call h5gcreate_f(file_id, trim(group_name), group_id, error, gcpl_id=gplist_id)
    if (error /= 0) write(*,*) 'error while creating a group, ', group_name
    call check_error()
    call h5gclose_f(group_id, error)
    if (error /= 0) write(*,*) 'error while closing a group, ' , group_name
    call check_error()
    call h5_close_group_prop_list()
  end subroutine h5_create_group_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_open_group_p(group_name)
    implicit none
    character(len=*), intent(in) :: group_name
    call h5_open_group_prop_list()
    call h5gopen_f(file_id, trim(group_name), group_id, error, gapl_id=gplist_id)
    if (error /= 0) write(*,*) 'hdf5 open group failed for ', group_name
    call check_error()
    call h5_close_group_prop_list()
    store_group_name = group_name
  end subroutine h5_open_group_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_scalar_i(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, intent(out) :: data
    integer, dimension(1) :: rdata
    integer(HSIZE_T), dimension(1)    :: dim
    dim = shape(rdata)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, rdata, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
    data = rdata(1)
  end subroutine h5_read_dataset_p_scalar_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_scalar_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), intent(out) :: data
    real(kind=CUSTOM_REAL), dimension(1) :: rdata
    integer(HSIZE_T), dimension(1)    :: dim
    dim = shape(rdata)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, rdata, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, rdata, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
    data = rdata(1)
  end subroutine h5_read_dataset_p_scalar_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_1d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:), intent(inout) :: data
    integer(HSIZE_T), dimension(1)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_1d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_1d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(inout) :: data
    integer(HSIZE_T), dimension(1)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_1d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_1d_l(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    logical, dimension(:), intent(inout) :: data
    integer, dimension(:), allocatable :: ldata
    integer(HSIZE_T), dimension(1)    :: dim
    integer :: lsize
    dim = shape(data)
    lsize = size(data)
    allocate(ldata(lsize),stat=error)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, ldata, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
    call int_array2bool(ldata, data)

    deallocate(ldata, stat=error)
  end subroutine h5_read_dataset_p_1d_l

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_2d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(2)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_2d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_2d_d(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    double precision, dimension(:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(2)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_2d_d

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_2d_c(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    character(len=*), dimension(:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(2)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_CHARACTER, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_2d_c

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_2d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(2)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_2d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_3d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:,:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(3)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_3d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_3d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(3)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_3d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_4d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(4)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_4d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_5d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(5)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, xfer_prp=plist_id)
    endif
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_5d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_attribute_p(attribute_name, dataset_name, data)
    implicit none
    character(len=*), intent(in)         :: attribute_name
    character(len=*), intent(in)         :: dataset_name
    integer(HID_T)                       :: attr_id
    integer, dimension(:), intent(inout) :: data
    integer(HSIZE_T), dimension(1)       :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5aopen_f(dataset_id, attribute_name, attr_id, error)
    call check_error()
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, data, dim, error)
    call check_error()
    call h5aclose_f(attr_id,error)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_attribute_p

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_p_4d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in) :: dataset_name
    integer, dimension(:,:,:,:), intent(inout) :: data
    integer(HSIZE_T), dimension(4)    :: dim
    dim = shape(data)

    call h5dopen_f(group_id, dataset_name, dataset_id, error)
    call check_error()
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, xfer_prp=plist_id)
    call check_error()
    call h5pclose_f(plist_id, error)
    call check_error()
    call h5dclose_f(dataset_id, error)
    call check_error()
  end subroutine h5_read_dataset_p_4d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_scalar_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, intent(inout), target                              :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    !dim = shape(data)
    dim = (/1/)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_scalar_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_scalar_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), intent(inout), target               :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    !dim = shape(data)
    dim = (/1/)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_scalar_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_1d_l_collect_hyperslab(dataset_name, data_out, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    logical, dimension(:), intent(inout), target                :: data_out
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer, dimension(:), allocatable, target                  :: data
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data_out)
    offset = offset_in ! convert data type
    allocate(data(dim(1)))

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    f_ptr = c_loc(data(1))
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                   file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()

    call int_array2bool(data, data_out)
    deallocate(data)
  end subroutine h5_read_dataset_1d_l_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_1d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:), intent(inout), target                :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    call h5_check_arr_dim(dim)

    ! write array using Fortran pointer
    !call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    ! use F2003 API
    f_ptr = c_loc(data(1))
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_1d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 1d array to global 1d array
  subroutine h5_read_dataset_1d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(inout), target :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_1d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 1d array to global 1d array
  subroutine h5_read_dataset_1d_r_collect_hyperslab_in_group(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(inout), target :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_1d_r_collect_hyperslab_in_group

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_2d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:), intent(inout), target              :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)
    count(2)  = dim(2)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1))

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_2d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_2d_i_collect_hyperslab_in_group(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:), intent(inout), target              :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)
    count(2)  = dim(2)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1))

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_2d_i_collect_hyperslab_in_group

!
!-------------------------------------------------------------------------------
!

  ! store local 2d array to global 2d array
  subroutine h5_read_dataset_2d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(inout), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)! nrec
    count(2)  = dim(2) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_2d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_2d_d_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    double precision, dimension(:,:), intent(inout), target        :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1)  = dim(1)! nrec
    count(2)  = dim(2) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1))

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_2d_d_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_3d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:,:), intent(inout), target            :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 3
    integer(HSIZE_T), dimension(3)                              :: dim
    integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1,1))

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_3d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_3d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:), intent(inout), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 3
    integer(HSIZE_T), dimension(3)                              :: dim
    integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1,1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_3d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_4d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:,:,:), intent(inout), target          :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 4
    integer(HSIZE_T), dimension(4)                              :: dim
    integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial
    count(4) = dim(4) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1,1,1))

    ! write array using Fortran pointer
    call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                   file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_4d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_4d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 4
    integer(HSIZE_T), dimension(4)                              :: dim
    integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial
    count(4) = dim(4) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1,1,1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_4d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_read_dataset_5d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(inout), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 5
    integer(HSIZE_T), dimension(5)                              :: dim
    integer(HSIZE_T), dimension(5)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(5)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial
    count(4) = dim(4) ! NSTEP partial
    count(5) = dim(5) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    f_ptr = c_loc(data(1,1,1,1,1))

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dread_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_read_dataset_5d_r_collect_hyperslab

!-------------------------------------------------------------------------------
!
! parallel write routines
!
!-------------------------------------------------------------------------------

!
! collective writers
!

  subroutine h5_write_dataset_p_1d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)         :: dataset_name
    integer, dimension(:), intent(in), target    :: data
    integer                              :: rank = 1
    integer(HSIZE_T), dimension(1)       :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    ! dummy array for generate 0 length dataset
    integer, dimension(1) :: dummy_1d_array = (/0/)

    dim = shape(data)

    ! add dummy 0 for no_element array
    if (size(data) == 0) then
      dim = 1
    endif

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    if (size(data) == 0) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, dummy_1d_array, dim, error,xfer_prp=plist_id)
    else
      f_ptr = c_loc(data(1))
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_1d_i

!
!-------------------------------------------------------------------------------
!

  ! write 1d integer and attribute
  subroutine h5_write_dataset_p_1d_ia(dataset_name, data, attribute_name, attr_data)
    implicit none
    character(len=*), intent(in)         :: dataset_name
    character(len=*), intent(in)         :: attribute_name
    integer, dimension(:), intent(in), target    :: data
    integer, intent(in)                  :: attr_data
    integer                              :: rank = 1
    integer(HSIZE_T), dimension(1)       :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    ! dummy array for generate 0 length dataset
    integer, dimension(1) :: dummy_1d_array = (/0/)

    dim = shape(data)

    ! add dummy 0 for no_element array
    if (size(data) == 0) then
        dim = 1
    endif

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    if (size(data) == 0) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, dummy_1d_array, dim, error,xfer_prp=plist_id)
    else
      f_ptr = c_loc(data(1))
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    ! add attribute
    call h5_add_attribute_i(attribute_name, (/attr_data/))

    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_1d_ia

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_1d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name


    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, CUSTOM_REAL, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1))
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                      xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_1d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_2d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)         :: dataset_name
    integer, dimension(:,:), intent(in), target    :: data
    integer                              :: rank = 2
    integer(HSIZE_T), dimension(2)       :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()

    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_2d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_2d_ia(dataset_name, data, attribute_name, attr_data)
    implicit none
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
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()

    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    ! add attribute
    call h5_add_attribute_i(attribute_name, (/attr_data/))

    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_2d_ia

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_2d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                  :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target    :: data
    integer                                                       :: rank = 2
    integer(HSIZE_T), dimension(2)                                :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, CUSTOM_REAL, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1))
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                      xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_2d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_2d_d(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                  :: dataset_name
    double precision, dimension(:,:), intent(in), target    :: data
    integer                                                       :: rank = 2
    integer(HSIZE_T), dimension(2)                                :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 8, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                   xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_2d_d

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_3d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)                     :: dataset_name
    integer, dimension(:,:,:), intent(in), target    :: data
    integer                                          :: rank = 3
    integer(HSIZE_T), dimension(3)                   :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1,1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_3d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_3d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                  :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:), intent(in), target  :: data
    integer                                                       :: rank = 3
    integer(HSIZE_T), dimension(3)                                :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, CUSTOM_REAL, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()

    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1,1))
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                      xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_3d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_4d_i(dataset_name, data)
    implicit none
    character(len=*), intent(in)                     :: dataset_name
    integer, dimension(:,:,:,:), intent(in), target  :: data
    integer                                          :: rank = 4
    integer(HSIZE_T), dimension(4)                   :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1,1,1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_4d_i

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_4d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                   :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in), target :: data
    integer                                                        :: rank = 4
    integer(HSIZE_T), dimension(4)                                 :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, CUSTOM_REAL, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1,1,1))
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                      xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_4d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_5d_r(dataset_name, data)
    implicit none
    character(len=*), intent(in)                                      :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in), target  :: data
    integer                                                           :: rank = 5
    integer(HSIZE_T), dimension(5)                                    :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, CUSTOM_REAL, file_id)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()

    ! write array using Fortran pointer
    f_ptr = c_loc(data(1,1,1,1,1))
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, f_ptr, error, &
                      xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                      xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
  end subroutine h5_write_dataset_p_5d_r

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_p_1d_l(dataset_name, data_in)
    ! writer for logical array
    ! logical array will be converted to integer array before write
    implicit none
    character(len=*), intent(in)               :: dataset_name
    logical, dimension(:), intent(in)          :: data_in
    integer, dimension(:), allocatable, target :: data

    integer                            :: rank = 1
    integer(HSIZE_T), dimension(1)     :: dim

    character(len=10) :: tempstr
    character(len=5)  :: gname_proc_head = "proc_"
    character(len=64) :: group_name

    dim = shape(data_in)

    ! create dataset
    call create_dataset_collect(dataset_name, dim, rank, 1, file_id)

    allocate(data(size(data_in)), stat=error)
    call bool_array2integer(data_in, data)

    write(tempstr, "(i6.6)") this_rank
    group_name = gname_proc_head // trim(tempstr)

    !! create datasets of all processes for independent write
    call h5_open_group(group_name)
    call h5_open_dataset(trim(dataset_name))

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call check_error()
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    call check_error()
    ! write array using Fortran pointer
    f_ptr = c_loc(data(1))
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, f_ptr, error, &
                    xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5_close_dataset()
    call h5_close_group()
    deallocate(data)
  end subroutine h5_write_dataset_p_1d_l

!
!-------------------------------------------------------------------------------
!

  ! from 1d local array to 2d global array
  subroutine h5_write_dataset_1d_to_2d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
    integer, dimension(2), intent(in)                           :: offset_in ! the position where the datablock is inserted
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset

    dim = size(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! size of data array inserted.
    count(1) = dim(1)
    count(2) = 1

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_to_2d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 2d array to global 3d array
  subroutine h5_write_dataset_2d_to_3d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = 3! NDIM
    count(2) = dim(2) ! NSTEP partial
    count(3) = 1 ! nrec

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_to_3d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_1d_l_collect_hyperslab(dataset_name, data_in, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    logical, dimension(:), intent(in), target                   :: data_in
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer, dimension(:), allocatable, target                  :: data
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data_in)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! if_collective == .false., this function gather all data from procs then write in a file at once

    ! select a place where data is inserted.
    count(1) = dim(1)
    allocate(data(dim(1)))

    ! convert logical array to  integer array
    call bool_array2integer(data_in, data)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    !call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    ! F2003 API
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
    deallocate(data)

  end subroutine h5_write_dataset_1d_l_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 1d array to global 1d array
  subroutine h5_write_dataset_1d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:), intent(in), target                   :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    call h5_check_arr_dim(dim)
    ! write array using Fortran pointer
    !call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    ! use F2003 API
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_1d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      !call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1)),error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      !call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 1d array to global 1d array
  subroutine h5_write_dataset_1d_r_collect_hyperslab_in_group(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:), intent(in), target    :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 1
    integer(HSIZE_T), dimension(1)                              :: dim
    integer(HSIZE_T), dimension(1)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(1)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      !call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      !call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_1d_r_collect_hyperslab_in_group

!
!-------------------------------------------------------------------------------
!

  ! store local 2d array to global 2d array
  subroutine h5_write_dataset_2d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:), intent(in), target                 :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)
    count(2) = dim(2)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    !call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1,1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_2d_i_collect_hyperslab_in_group(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:), intent(in), target                 :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)
    count(2) = dim(2)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    !call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1,1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_i_collect_hyperslab_in_group

!
!-------------------------------------------------------------------------------
!

  ! store local 2d array to global 2d array
  subroutine h5_write_dataset_2d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! nrec
    count(2) = dim(2) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      !call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      !call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dim, error, &
      !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_2d_d_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    double precision, dimension(:,:), intent(in), target        :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 2
    integer(HSIZE_T), dimension(2)                              :: dim
    integer(HSIZE_T), dimension(2)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(2)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! nrec
    count(2) = dim(2) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1,1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_2d_d_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_3d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                  :: dataset_name
    integer, dimension(:,:,:), intent(in), target :: data
    integer, dimension(:), intent(in)             :: offset_in
    logical, intent(in)                           :: if_collective
    ! local parameters
    integer                                       :: rank = 3
    integer(HSIZE_T), dimension(3)                :: dim
    integer(HSIZE_T), dimension(3)                :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(3)               :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    !call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, data, dim, error, &
    !                file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    ! F2003 API
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1,1,1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_3d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  ! store local 3d array to global 3d array
  subroutine h5_write_dataset_3d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:), intent(in), target  :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 3
    integer(HSIZE_T), dimension(3)                              :: dim
    integer(HSIZE_T), dimension(3)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(3)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1) ! NDIM
    count(2) = dim(2) ! nrec
    count(3) = dim(3) ! NSTEP partial

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_3d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_4d_i_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                :: dataset_name
    integer, dimension(:,:,:,:), intent(in), target             :: data
    integer, dimension(:), intent(in)                           :: offset_in
    logical, intent(in)                                         :: if_collective
    ! local parameters
    integer                                                     :: rank = 4
    integer(HSIZE_T), dimension(4)                              :: dim
    integer(HSIZE_T), dimension(4)                              :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(4)                             :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)
    count(2) = dim(2)
    count(3) = dim(3)
    count(4) = dim(4)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, c_loc(data(1,1,1,1)), error, &
                    file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_4d_i_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_4d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                   :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(in), target :: data
    integer, dimension(:), intent(in)                              :: offset_in
    logical, intent(in)                                            :: if_collective
    ! local parameters
    integer                                                        :: rank = 4
    integer(HSIZE_T), dimension(4)                                 :: dim
    integer(HSIZE_T), dimension(4)                                 :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(4)                                :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)
    count(2) = dim(2)
    count(3) = dim(3)
    count(4) = dim(4)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1,1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1,1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_4d_r_collect_hyperslab

!
!-------------------------------------------------------------------------------
!

  subroutine h5_write_dataset_5d_r_collect_hyperslab(dataset_name, data, offset_in, if_collective)
    implicit none
    character(len=*), intent(in)                                     :: dataset_name
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in), target :: data
    integer, dimension(:), intent(in)                                :: offset_in
    logical, intent(in)                                              :: if_collective
    ! local parameters
    integer                                                          :: rank = 5
    integer(HSIZE_T), dimension(5)                                   :: dim
    integer(HSIZE_T), dimension(5)                                   :: count ! size of hyperslab
    integer(HSSIZE_T), dimension(5)                                  :: offset ! the position where the datablock is inserted

    dim = shape(data)
    offset = offset_in ! convert data type

    ! open dataset
    call h5_open_dataset2(trim(dataset_name))

    ! select a place where data is inserted.
    count(1) = dim(1)
    count(2) = dim(2)
    count(3) = dim(3)
    count(4) = dim(4)
    count(5) = dim(5)

    ! select hyperslab in the file
    call h5screate_simple_f(rank,count, mem_dspace_id, error)
    call check_error()
    call h5dget_space_f(dataset_id, file_dspace_id, error)
    call check_error()
    call h5sselect_hyperslab_f(file_dspace_id, H5S_SELECT_SET_F, offset, count, error)
    call check_error()
    call h5_create_dataset_prop_list(if_collective)

    ! write array using Fortran pointer
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dataset_id, H5T_NATIVE_REAL, c_loc(data(1,1,1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    else
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, c_loc(data(1,1,1,1,1)), error, &
                      file_space_id=file_dspace_id, mem_space_id=mem_dspace_id, xfer_prp=plist_id)
    endif
    if (error /= 0) write(*,*) 'hdf5 dataset write failed for ', dataset_name
    call check_error()
    call h5_close_prop_list(dataset_name)
    call h5sclose_f(mem_dspace_id, error)
    call check_error()
    call h5sclose_f(file_dspace_id, error)
    call check_error()
    call h5_close_dataset()
  end subroutine h5_write_dataset_5d_r_collect_hyperslab

!-------------------------------------------------------------------------------
!
! other utilities
!
!-------------------------------------------------------------------------------

  subroutine check_error()
    implicit none
    if (error /= 0) then
      print *,'Error: HDF5 routine call got error ',error
      stop 'HDF5 error'
    endif
  end subroutine check_error

!
!-------------------------------------------------------------------------------
!

  subroutine bool_array2integer(boolarray, intarray)
    implicit none
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

!
!-------------------------------------------------------------------------------
!

  subroutine int_array2bool(intarray, boolarray)
    implicit none
    integer, dimension(:), intent(in)  :: intarray
    logical, dimension(:), intent(out) :: boolarray
    integer                            :: array_size,i

    array_size = size(intarray)
    boolarray(:) = .false.

    do i = 1, array_size
      if (intarray(i) /= 0) then
        boolarray(i) = .true.
      endif
    enddo
  end subroutine int_array2bool

!
!-------------------------------------------------------------------------------
!

  subroutine h5_gather_dsetsize(dimin, data_rank, all_dim)

    implicit none
    integer, intent(in)                         :: data_rank
    integer(HSIZE_T), dimension(:), intent(in)  :: dimin
    integer, dimension(data_rank)               :: dim
    integer, dimension(data_rank,0:total_proc-1), intent(out)      :: all_dim

    dim = dimin ! convert integer 8(HSIZE_T) to 4

    if (data_rank == 1) then
      call gather_all_all_singlei(dim(1),all_dim,total_proc)
      !debug
      !if (this_rank == 0) then
      !  print *, "rank: ", this_rank &
      !         ,"dimin", dimin, ", kind: ", kind(dimin) &
      !         ,"dim", dim, ", kind: ", kind(dim) &
      !         ,"shape dimin", shape(dimin) &
      !         ,"data_rank", data_rank &
      !         ,"all dim", all_dim &
      !         ,"shape all dim", shape(all_dim) &
      !         ,"kind all_dim", kind(all_dim(1,0)) &
      !         ,"total proc " , total_proc &
      !         ,"kind mpi_integer ", kind(MPI_INTEGER)
      !endif
    else
      ! gather only on main rank
      !call gather_all_i(dim, data_rank, all_dim, data_rank, total_proc)
      ! gather on all processes
      call gather_all_all_i(dim, data_rank, all_dim, data_rank, total_proc)
    endif

    call synchronize_all()
  end subroutine

!
!-------------------------------------------------------------------------------
!

  subroutine create_dataset_collect(dataset_name, dim, data_rank, dtype_id, base_id)
    implicit none
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
    call h5_gather_dsetsize(dim, data_rank, all_dim)
    ! create all datasets of all process by rank0

    !print *, dataset_name, ": debug returned all dim", all_dim
    do iproc = 0, total_proc -1
      write(tempstr, "(i6.6)") iproc
      group_name = gname_proc_head // trim(tempstr)
      group_and_dataset_name = trim(group_name) //"/"// dataset_name

      dim_h5 = all_dim(:,iproc)

      call h5screate_simple_f(data_rank, dim_h5, dspace_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataspace create failed for ', dataset_name
      call check_error()
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call check_error()
      ! This is required for this data pattern
      call H5Pset_alloc_time_f(plist_id, H5D_ALLOC_TIME_EARLY_F, error)
      call check_error()
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
      call check_error()
      call h5_close_prop_list(dataset_name)
      call h5dclose_f(dataset_id,error)
      if (error /= 0) write(*,*) 'hdf5 dataset close failed for ', dataset_name
      call check_error()
      call h5sclose_f(dspace_id, error)
      if (error /= 0) write(*,*) 'hdf5 dataspace close failed for ', dataset_name
      call check_error()
    enddo

    call synchronize_all()
  end subroutine create_dataset_collect

!-------------------------------------------------------------------------------
!
! object-oriented interface
!
!-------------------------------------------------------------------------------
! Fortran 2003 standard
!
! this interface is mainly for testing purposes, and for checking if compilers by now support these features :(
!
! in the past, we had issues with compilers due to keywords like class().
! here, we use strictly Fortran 2003 features, like polymorphism with 'class()', object destructor procedure by 'final ::',
! and type binding by a pass statement like 'procedure :: . => ..'
!
! for a status report on compiler support, see Fortran+2003+status on:
!   https://Fortranwiki.org/
! Cray, GNU, IBM and Intel compilers should well support these features.

  ! constructors

  function h5io_constructor() result(this)
    implicit none
    type(h5io) :: this
    call h5_initialize()
    this%is_initialized = .true.
  end function h5io_constructor

  function h5io_constructor_from_file(filename) result(this)
    implicit none
    type(h5io) :: this
    character(*), intent(in) :: filename
    call h5_initialize()
    this%is_initialized = .true.
    call this%open(filename)
    this%filename = trim(filename)
  end function h5io_constructor_from_file

  ! destructor

  subroutine h5io_destructor(this)
    implicit none
    type(h5io) :: this
    call h5_finalize()
    this%is_initialized = .false.
  end subroutine h5io_destructor

  ! object procedures

  ! note: with the pass statement "procedure :: open => h5io_open_file", the routine must have at least one argument and
  !       the argument 'this' must be defined as polymorphic, i.e., by class() and not type() as above

  subroutine h5io_open_file(this,filename,status)
    implicit none
    class(h5io) :: this
    character(*), intent(in) :: filename
    character(*), intent(in), optional :: status
    ! checks if initialized
    if (.not. this%is_initialized) then
      call h5_initialize()
      this%is_initialized = .true.
    endif
    ! only main process creates file if needed
    if (present(status)) then
      if (status == 'new') call h5_create_file(filename)
    endif
    ! opens file
    call h5_open_file(filename)
    this%f_id = file_id
    this%filename = trim(filename)
  end subroutine h5io_open_file

  subroutine h5io_close_file(this)
    implicit none
    class(h5io) :: this
    call h5_close_file()
    this%f_id = file_id
  end subroutine h5io_close_file

  subroutine h5io_write_dataset_i(this,dname,data_i,group)
    implicit none
    class(h5io) :: this
    character(len=*), intent(in) :: dname
    integer, dimension(:), intent(in) :: data_i
    character(len=*), intent(in), optional :: group
    ! writes either with or without group
    if (present(group)) then
      call h5_open_or_create_group(group)
      call h5_write_dataset(dname,data_i)
    else
      call h5_write_dataset_no_group(dname,data_i)
    endif
    this%d_id = dataset_id
  end subroutine h5io_write_dataset_i

  subroutine h5io_write_dataset_r(this,dname,data_r,group)
    implicit none
    class(h5io) :: this
    character(len=*), intent(in) :: dname
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: data_r
    character(len=*), intent(in), optional :: group
    ! writes either with or without group
    if (present(group)) then
      call h5_open_or_create_group(group)
      call h5_write_dataset(dname,data_r)
    else
      call h5_write_dataset_no_group(dname,data_r)
    endif
    this%d_id = dataset_id
  end subroutine h5io_write_dataset_r

#endif

end module manager_hdf5

!-----------------------------------------------------

#if defined(USE_HDF5)

! test function for object-oriented interface

  subroutine test_io_hdf5()
    use constants, only: myrank
    use manager_hdf5, only: h5io
    implicit none
    type(h5io) :: h5
    integer :: store_x(10)
    integer :: i

    ! serial test by main process only
    if (myrank /= 0) return

    ! initialize
    store_x(:) = (/(i,i=1,10)/)

    ! hdf5
    ! calls the object constructor, i.e. h5io_constructor()
    h5 = h5io()

    ! open file
    call h5%open("tmp_test.h5",status='new')  ! object function call

    ! write out data
    call h5%write("x",store_x)

    ! close file
    call h5%close()

    ! h5 goes out of scope, will call the object destructor, i.e., h5io_destructor()
  end subroutine test_io_hdf5

#endif
