  module global_parameters


    character(len=256) working_axisem_dir,input_field_name
    character(len=256) input_veloc_name(3), output_veloc_name(3), input_stress_name(6), output_stress_name(6)

    character(len=256) input_displ_name(3), output_displ_name(3), input_deriv_name(9), output_deriv_name(9)

    character(len=256) output_field_name,input_point_file,input_point_file_cart
    integer, parameter :: CUSTOM_REAL=8,SINGLE_REAL=4
    real(kind=CUSTOM_REAL) lat_src,lon_src,lat_mesh,lon_mesh,azi_rot,dtt
    integer nbproc,nel,ntime
    integer ibeg,iend
    integer, allocatable :: nb_stored(:)
    real(kind=CUSTOM_REAL), allocatable :: scoor(:,:,:),zcoor(:,:,:)
    real(kind=CUSTOM_REAL), allocatable :: depth_ele(:)

    !! mesh specifications
    integer, parameter :: NGNOD=8,NGLLX=5,NGLLY=5,NGLLS=5,NGLLZ=5

    !! field arrays
    real(kind=SINGLE_REAL), allocatable :: data_read(:,:,:),data_rec(:,:),stress_rec(:,:),stress_to_write(:,:)
    real(kind=SINGLE_REAL), allocatable :: strain_rec(:,:)

    real(kind=SINGLE_REAL), allocatable :: deriv_rec(:,:) !! CD CD

    real(kind=SINGLE_REAL), allocatable :: data_tmpKH_rec(:,:,:), deriv_tmpKH_rec(:,:,:) !! CD CD

    !! work arrays
    real(kind=SINGLE_REAL), allocatable :: data_reduce(:,:),stress_reduce(:,:),deriv_reduce(:,:)

    !! input box point
    integer nbrec
    real(kind=CUSTOM_REAL), allocatable :: reciever_cyl(:,:),reciever_geogr(:,:),reciever_sph(:,:)
    real(kind=CUSTOM_REAL), allocatable :: reciever_interp_value(:),xi_rec(:),eta_rec(:)
    integer, allocatable :: rec2elm(:),ele_candidate(:,:),ele_seen(:)
    logical, allocatable :: up(:)

    !!
    real(kind=CUSTOM_REAL), parameter :: pi=3.1415926535898_CUSTOM_REAL
    real(kind=CUSTOM_REAL), parameter :: pick_tresh=0.1
    integer ipick

    !! rotation matrix with recpect to the source
    real(kind=CUSTOM_REAL) rot_mat(3,3),trans_rot_mat(3,3)

    !! rotation matrix with respcet to the mesh
    real(kind=CUSTOM_REAL) rot_mat_mesh(3,3),trans_rot_mat_mesh(3,3)
    real(kind=CUSTOM_REAL) mat(3,3),tmat(3,3)
    real(kind=CUSTOM_REAL) rot_azi_chunk(3,3),trans_rot_azi_chunk(3,3)

    !! info about simulation
    integer nsim
    character(len=100), allocatable :: simdir(:)
    character(len=10),  allocatable :: src_type(:,:)

    !! moment tensor
    real(kind=SINGLE_REAL)  :: Mij(6)
    real, allocatable, dimension(:)     :: magnitude

    !! post process
    real(kind=SINGLE_REAL), allocatable ::  f1(:),f2(:),phi(:)

    !! mpi
    integer irecmin, irecmax

    !! if the expand_2D_3D and reformat are used for KH integral and reciprocity
    logical :: recip_KH_integral
    integer Xk_force

  end module global_parameters
