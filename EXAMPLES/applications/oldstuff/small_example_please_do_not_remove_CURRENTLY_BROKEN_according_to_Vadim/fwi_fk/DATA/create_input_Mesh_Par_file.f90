module create_input_Mesh_Par_file

  implicit none

  type :: mesh_par

     double precision    :: LATITUDE_MIN                    = 0.d0
     double precision    :: LATITUDE_MAX                    = 32000.d0
     double precision    :: LONGITUDE_MIN                   = 0.d0
     double precision    :: LONGITUDE_MAX                   = 32000.d0
     double precision    :: DEPTH_BLOCK_KM                  = 5.d0
     integer             :: UTM_PROJECTION_ZONE             = 11
     logical             :: SUPPRESS_UTM_PROJECTION         = .true.

     character(len=100)  :: INTERFACES_FILE                 = interfaces.dat
     character(len=100)  :: CAVITY_FILE                     = no_cavity.dat
     integer             :: NEX_XI                          = 32
     integer             :: NEX_ETA                         = 32
     integer             :: NPROC_XI                        = 2
     integer             :: NPROC_ETA                       = 2
     logical             :: USE_REGULAR_MESH                = .true.
     integer             :: NDOUBLINGS                      = 0
     integer             :: NZ_DOUBLING_1                   = 11
     integer             :: NZ_DOUBLING_2                   = 0
     logical             :: CREATE_ABAQUS_FILES             = .false.
     logical             :: CREATE_DX_FILES                 = .false.
     logical             :: CREATE_VTK_FILES                = .false.
     character(len=100)  :: LOCAL_PATH                      = ./DATABASES_MPI
     double precision    :: THICKNESS_OF_X_PML              = 12.3d0
     double precision    :: THICKNESS_OF_Y_PML              = 12.3d0
     double precision    :: THICKNESS_OF_Z_PML              = 12.3d0
     integer             :: NMATERIALS                      = 4
     integer             :: NREGIONS                        = 2

     real, dimension(4,6)                 :: default_materials
     real, dimension(:,:), allocatable    :: custom_materials

     integer, dimension(2,7)              :: default_regions
     integer, dimension(:,:), allocatable :: custom_regions

  end type mesh_par


  character(len=100)                            :: line_to_read, line_to_write
  character(len=1), parameter                   :: separator='='
  character(len=256)                            :: reference_file_to_modify, file_modifed

  character(len=32), dimension(:), allocatable  :: FIELD_to_MODIF

!!$
!!$  !!
!!$  reference_file_to_modify='Mesh_Par_file_ref'
!!$  file_modifed='Mesh_Par_file_new'
!!$
!!$  open(10,file=trim(reference_file_to_modify))
!!$  open(11,file=trim(file_modifed))
!!$
!!$
!!$  open(12,file='input_field_to_modif.txt')
!!$
!!$


end module create_input_Mesh_Par_file
