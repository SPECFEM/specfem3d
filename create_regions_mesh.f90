!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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


  subroutine create_regions_mesh_ext_mesh(ibool, &
                    xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
                    nnodes_ext_mesh,nelmnts_ext_mesh, &
                    nodes_coords_ext_mesh, elmnts_ext_mesh, &
                    max_static_memory_size, mat_ext_mesh, materials_ext_mesh, &
                    nmat_ext_mesh, undef_mat_prop, nundefMat_ext_mesh, &
                    ninterface_ext_mesh, max_interface_size_ext_mesh, &
                    my_neighbours_ext_mesh, my_nelmnts_neighbours_ext_mesh, &
                    my_interfaces_ext_mesh, &
                    ibool_interfaces_ext_mesh, nibool_interfaces_ext_mesh, &
                    nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, &
                    NSPEC2D_BOTTOM, NSPEC2D_TOP,&
                    ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
                    nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax,&
                    nodes_ibelm_bottom,nodes_ibelm_top, &
                    SAVE_MESH_FILES,nglob)

! create the different regions of the mesh

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer :: nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: npointot

! proc numbers for MPI
  integer :: myrank

  character(len=150) :: LOCAL_PATH

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! static memory size needed by the solver
  double precision :: max_static_memory_size 
  
  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh

!pll
  integer :: nmat_ext_mesh,nundefMat_ext_mesh 
  double precision, dimension(5,nmat_ext_mesh) :: materials_ext_mesh  
  character (len=30), dimension(5,nundefMat_ext_mesh):: undef_mat_prop
  
!  double precision, external :: materials_ext_mesh
  integer :: ninterface_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,ninterface_ext_mesh) :: my_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,ninterface_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(ninterface_ext_mesh) :: nibool_interfaces_ext_mesh

! absorbing boundaries
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
!  integer  :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer, dimension(nspec2D_xmin)  :: ibelm_xmin  
  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top
  ! node indices of boundary faces
  integer, dimension(4,nspec2D_xmin)  :: nodes_ibelm_xmin  
  integer, dimension(4,nspec2D_xmax)  :: nodes_ibelm_xmax
  integer, dimension(4,nspec2D_ymin)  :: nodes_ibelm_ymin
  integer, dimension(4,nspec2D_ymax)  :: nodes_ibelm_ymax
  integer, dimension(4,NSPEC2D_BOTTOM)  :: nodes_ibelm_bottom
  integer, dimension(4,NSPEC2D_TOP)  :: nodes_ibelm_top

  logical :: SAVE_MESH_FILES
  integer :: nglob

! local parameters
!-----------------------    

! for MPI buffers
!  integer, dimension(:), allocatable :: reorder_interface_ext_mesh,ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh
!  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh_true
  !integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
!  integer, dimension(:), allocatable :: ibool_interface_ext_mesh_dummy
!  double precision, dimension(:), allocatable :: work_ext_mesh
  
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision, dimension(:), allocatable :: xelm,yelm,zelm

! static memory size needed by the solver
  double precision :: static_memory_size

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! for model density
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore,vpstore,vsstore 
! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! attenuation 
  integer, dimension(:,:,:,:), allocatable :: iflag_attenuation_store

! 2D shape functions and their derivatives, weights
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top
  double precision, dimension(:,:), allocatable :: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

! absorbing boundaries
! pll 
!  logical, dimension(:,:),allocatable :: iboun  
!  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xcoord_iboun,ycoord_iboun,zcoord_iboun
! 2-D jacobians and normals
!  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
!       jacobian2D_xmin,jacobian2D_xmax, &
!       jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_top
  
!  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
!    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_top

!  ! local indices i,j,k of all GLL points on xmin boundary in the element
!  integer,dimension(:,:,:,:),allocatable :: ibelm_gll_xmin,ibelm_gll_xmax, &
!                                          ibelm_gll_ymin,ibelm_gll_ymax, &
!                                          ibelm_gll_bottom,ibelm_gll_top
!  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorbing_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: absorbing_boundary_jacobian2D
  integer, dimension(:,:,:), allocatable :: absorbing_boundary_ijk
  integer, dimension(:), allocatable :: absorbing_boundary_ispec
  integer :: num_absorbing_boundary_faces

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! variables for creating array ibool (some arrays also used for AVS or DX files)
!  integer, dimension(:), allocatable :: locval !,iglob
!  logical, dimension(:), allocatable :: ifseg
!  double precision, dimension(:), allocatable :: xp,yp,zp

!  integer :: ilocnum,ier,iinterface !,ieoff
  integer, dimension(:), allocatable :: elem_flag
  integer :: ier
  integer :: i,j,k,ispec,iglobnum
!  integer  :: ispec2D

! name of the database file
  character(len=150) prname
  character(len=150) prname_file

! mask to sort ibool
!  integer, dimension(:), allocatable :: mask_ibool
!  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
!  integer :: inumber
  
! memory test
!  logical,dimension(:),allocatable :: test_mem 


! for vtk output
!  integer,dimension(:),allocatable :: itest_flag

! For Piero Basini :
! integer :: doubling_value_found_for_Piero
!   double precision :: xmesh,ymesh,zmesh
!   double precision :: rho,vp,vs

!   integer,dimension(nspec) ::  idoubling
!   integer :: doubling_value_found_for_Piero
!   integer, parameter :: NUMBER_OF_STATIONS = 6
!   double precision, parameter :: RADIUS_TO_EXCLUDE = 250.d0
!   double precision, dimension(NUMBER_OF_STATIONS) :: utm_x_station,utm_y_station

!   logical :: is_around_a_station
!   integer :: istation

! ! store bedrock values
!   integer ::  icornerlat,icornerlong
!   double precision ::  lat,long,elevation_bedrock
!   double precision ::  lat_corner,long_corner,ratio_xi,ratio_eta
  
! ! size of topography and bathymetry file for Piero Basini's model
!   integer, parameter :: NX_TOPO = 787, NY_TOPO = 793
!   double precision, parameter :: ORIG_LAT_TOPO = -102352.d0
!   double precision, parameter :: ORIG_LONG_TOPO = 729806.d0
!   character(len=150), parameter :: TOPO_FILE = 'DATA/piero_model/dem_EV_UTM_regular_250_reordered.dat'
! ! for Piero Basini's model this is the resolution in meters of the topo file
!   double precision, parameter :: DEGREES_PER_CELL_TOPO = 250.d0

!real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ibedrock


! **************


  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) 
    write(IMAIN,*) '  ...allocating arrays '
  endif

! tests memory availability (including some small buffer of 10*1024 byte)
!  allocate( test_mem(int(max_static_memory_size)+10*1024),stat=ier)
!  if(ier /= 0) then
!    write(IMAIN,*) 'error: try to increase the available process stack size by'
!    write(IMAIN,*) '       ulimit -s **** '
!    call exit_MPI(myrank,'not enough memory to allocate arrays')
!  endif
!  test_mem(:) = .true.
!  deallocate( test_mem, stat=ier) 
!  if(ier /= 0) call exit_MPI(myrank,'error to allocate arrays')
!  call sync_all()

          
  allocate( xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),stat=ier)

  allocate( iflag_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), &
          dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)

! pll 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
          shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
          shape2D_bottom(NGNOD2D,NGLLX,NGLLY), &
          shape2D_top(NGNOD2D,NGLLX,NGLLY), stat=ier)
  
  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
          dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
          dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY), &
          dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  
  allocate(wgllwgll_xy(NGLLX,NGLLY), &
          wgllwgll_xz(NGLLX,NGLLZ), &
          wgllwgll_yz(NGLLY,NGLLZ),stat=ier)  
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! pll Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
          kappastore(NGLLX,NGLLY,NGLLZ,nspec), &
          mustore(NGLLX,NGLLY,NGLLZ,nspec), &
          vpstore(NGLLX,NGLLY,NGLLZ,nspec), &
          vsstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier) !pll
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec), &
          xiystore(NGLLX,NGLLY,NGLLZ,nspec), &
          xizstore(NGLLX,NGLLY,NGLLZ,nspec), &
          etaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
          etaystore(NGLLX,NGLLY,NGLLZ,nspec), &
          etazstore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammaystore(NGLLX,NGLLY,NGLLZ,nspec), &
          gammazstore(NGLLX,NGLLY,NGLLZ,nspec), &
          jacobianstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! allocates arrays for Stacey boundaries
!  allocate( nimin(2,NSPEC2DMAX_YMIN_YMAX),nimax(2,NSPEC2DMAX_YMIN_YMAX), &
!          njmin(2,NSPEC2DMAX_XMIN_XMAX),njmax(2,NSPEC2DMAX_XMIN_XMAX), &
!          nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX),nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX),stat=ier)
!  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

!  ! local indices i,j,k of all GLL points on xmin boundary in the element
!  allocate(ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin),ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax), &
!            ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin),ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax), &
!            ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom),ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top),stat=ier)          
!  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')  

!  ! pll 2-D jacobians and normals
!  allocate(jacobian2D_xmin(NGLLY,NGLLZ,nspec2D_xmin),jacobian2D_xmax(NGLLY,NGLLZ,nspec2D_xmax), &
!          jacobian2D_ymin(NGLLX,NGLLZ,nspec2D_ymin),jacobian2D_ymax(NGLLX,NGLLZ,nspec2D_ymax), &
!          jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM),jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
!  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

!  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,nspec2D_xmin),normal_xmax(NDIM,NGLLY,NGLLZ,nspec2D_xmax), &
!          normal_ymin(NDIM,NGLLX,NGLLZ,nspec2D_ymin),normal_ymax(NDIM,NGLLX,NGLLZ,nspec2D_ymax), &
!          normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM),normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
!  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! free surface
  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP),&
          normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! absorbing boundary 
  ! absorbing faces
  num_absorbing_boundary_faces = nspec2D_xmin + nspec2D_xmax + nspec2D_ymin + nspec2D_ymax + nspec2D_bottom
  ! free surface also absorbs
  if( ABSORB_FREE_SURFACE ) num_absorbing_boundary_faces = num_absorbing_boundary_faces + nspec2D_top

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( absorbing_boundary_ispec(num_absorbing_boundary_faces), &
           absorbing_boundary_ijk(3,NGLLSQUARE,num_absorbing_boundary_faces), &
           absorbing_boundary_jacobian2D(NGLLSQUARE,num_absorbing_boundary_faces), &
           absorbing_boundary_normal(NDIM,NGLLSQUARE,num_absorbing_boundary_faces),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')


! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
! returns jacobianstore,xixstore,...gammazstore
! and GLL-point locations in xstore,ystore,zstore
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up jacobian '
  endif
  
  call create_regions_mesh_ext_mesh_setup_jacobian(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                      myrank,shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                      dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                      wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                      xstore,ystore,zstore,nspec,xelm,yelm,zelm, &
                      nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh, &
                      xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                      gammaxstore,gammaystore,gammazstore, &
                      jacobianstore)

! sets material velocities
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...determining kappa and mu parameters'
  endif

  call create_regions_mesh_ext_mesh_determine_kappamu(nspec,mat_ext_mesh,nelmnts_ext_mesh,&
                        materials_ext_mesh,nmat_ext_mesh,&
                        undef_mat_prop,nundefMat_ext_mesh,&
                        rhostore,kappastore,mustore,vpstore,vsstore,&
                        iflag_attenuation_store,rho_vp,rho_vs)

! creates ibool index array for projection from local to global points
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...indexing global points'
  endif

  call create_regions_mesh_ext_mesh_setup_global_indexing(ibool, &
           xstore,ystore,zstore,nspec,nglob,npointot, &
           nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! unique global point locations
  allocate(xstore_dummy(nglob), &
          ystore_dummy(nglob), &
          zstore_dummy(nglob),stat=ier) 
  if(ier /= 0) stop 'error in allocate'  
  do ispec = 1, nspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              iglobnum = ibool(i,j,k,ispec)
              xstore_dummy(iglobnum) = xstore(i,j,k,ispec)
              ystore_dummy(iglobnum) = ystore(i,j,k,ispec)
              zstore_dummy(iglobnum) = zstore(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo  

! creating mass matrix (will be fully assembled with MPI in the solver)
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
  endif

  allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

  call create_regions_mesh_ext_mesh_create_mass_matrix(nglob,rmass,&
                  nspec,wxgll,wygll,wzgll,ibool,jacobianstore,rhostore)
          
! sets up absorbing/free surface boundaries  
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries '
  endif

  call create_regions_mesh_ext_mesh_setup_absorbing_bound(myrank,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                            wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                            normal_top,jacobian2D_top, &
                            absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                            absorbing_boundary_ijk,absorbing_boundary_ispec, &
                            num_absorbing_boundary_faces)

! sets up MPI interfaces between partitions
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...preparing MPI interfaces '
  endif
       
  call create_regions_mesh_ext_mesh_prepare_MPI_interfaces(nglob,nspec,ibool, &
                                    nelmnts_ext_mesh,elmnts_ext_mesh, &
                                    my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                                    ibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh, &
                                    ninterface_ext_mesh,max_interface_size_ext_mesh, &
                                    xstore_dummy,ystore_dummy,zstore_dummy)

! saves the binary files
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...saving databases'
  endif

  call create_name_database(prname,myrank,LOCAL_PATH)
  call save_arrays_solver_ext_mesh(nspec,nglob, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore, &
                        jacobianstore, rho_vp,rho_vs,iflag_attenuation_store, &
                        kappastore,mustore,rmass,ibool, &
                        xstore_dummy,ystore_dummy,zstore_dummy, &
                        NSPEC2D_TOP,ibelm_top,normal_top,jacobian2D_top, &
                        absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                        absorbing_boundary_ijk,absorbing_boundary_ispec, &
                        num_absorbing_boundary_faces, &
                        ninterface_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname,SAVE_MESH_FILES)

! computes the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),ninterface_ext_mesh,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)


! checks the mesh, stability and resolved period 
  call sync_all()
  call check_mesh_resolution(myrank,nspec,nglob,ibool,&
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            kappastore,mustore,rho_vp,rho_vs, &
                            -1.0d0 )

! VTK file output
  if( SAVE_MESH_FILES ) then
    ! saves material flag assigned for each spectral element into a vtk file 
    prname_file = prname(1:len_trim(prname))//'material_flag'
    allocate(elem_flag(nspec))
    elem_flag(:) = mat_ext_mesh(1,:)
    call save_arrays_solver_ext_mesh_elem_vtk(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            elem_flag,prname_file)
    deallocate(elem_flag)
    
    ! saves attenuation flag assigned on each gll point into a vtk file 
    prname_file = prname(1:len_trim(prname))//'attenuation_flag'
    call save_arrays_solver_ext_mesh_glldata_vtk(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            iflag_attenuation_store,prname_file)

    !daniel
    !plotting abs boundaries
    !  allocate(itest_flag(nspec))
    !  itest_flag(:) = 0
    !  do ispec=1,nspec
    !    if( iboun(1,ispec) ) itest_flag(ispec) = 1
    !  enddo
    !  prname_file = prname(1:len_trim(prname))//'iboundary1_flag'
    !  call save_arrays_solver_ext_mesh_elem_vtk(nspec,nglob, &
    !            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
    !            itest_flag,prname_file)
    !  deallocate(itest_flag)
  endif  

! AVS/DX file output
! create AVS or DX mesh data for the slice, edges and faces
!  if(SAVE_MESH_FILES) then
! check: no idoubling
!    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
!    call write_AVS_DX_mesh_quality_data(prname,nspec,xstore,ystore,zstore, &
!                   kappastore,mustore,rhostore)
! check: no iMPIcut_xi,iMPIcut_eta,idoubling
!    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
!              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
! check: no idoubling
!    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
!              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
!  endif

! cleanup
  deallocate(xixstore,xiystore,xizstore,&
            etaxstore,etaystore,etazstore,&
            gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore,iflag_attenuation_store)
  deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  deallocate(kappastore,mustore,rho_vp,rho_vs)

  end subroutine create_regions_mesh_ext_mesh

!
!----
!

subroutine create_regions_mesh_ext_mesh_setup_jacobian(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                      myrank,shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                      dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                      wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                      xstore,ystore,zstore,nspec,xelm,yelm,zelm, &
                      nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh, &
                      xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                      gammaxstore,gammaystore,gammazstore,&
                      jacobianstore)

  implicit none

  include 'constants.h'

! number of spectral elements in each block
  integer :: nspec

! Gauss-Lobatto-Legendre points and weights of integration
  double precision :: xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ),wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

! 3D shape functions and their derivatives
  double precision :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision :: dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

! 2D shape functions and their derivatives
  double precision :: shape2D_x(NGNOD2D,NGLLY,NGLLZ),shape2D_y(NGNOD2D,NGLLX,NGLLZ),&
                  shape2D_bottom(NGNOD2D,NGLLX,NGLLY),shape2D_top(NGNOD2D,NGLLX,NGLLY)
  double precision :: dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ),dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ),&
              dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY),dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  double precision,dimension(NGNOD) :: xelm,yelm,zelm

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xixstore,xiystore,xizstore, &
                        etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                        jacobianstore

! proc numbers for MPI
  integer :: myrank

  integer :: ispec,ia,i,j,k

!  integer :: ielm
!  logical :: inorder
  
! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! pll get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! 2D weights
  do j=1,NGLLY
    do i=1,NGLLX
      wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo
  do k=1,NGLLZ
    do i=1,NGLLX
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo
  do k=1,NGLLZ
    do j=1,NGLLY
      wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

! point locations
  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  do ispec = 1, nspec
    !call get_xyzelm(xelm, yelm, zelm, ispec, elmnts_ext_mesh, nodes_coords_ext_mesh, nspec, nnodes_ext_mesh)
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo


!daniel
!    ! do we have to test CUBIT order - or will 3D jacobian be defined?
!
!    ! bottom - top?
!    ! point 1 (0,0,0) vs point 5 (0,0,1)
!    inorder = .true.
!    if( nodes_coords(3,elmnts(1,num_elmnt)) > nodes_coords(3,elmnts(5,num_elmnt)) ) then
!      print*,num_elmnt,'z1-5 :',nodes_coords(3,elmnts(1,num_elmnt)),nodes_coords(3,elmnts(5,num_elmnt))
!      inorder = .false.
!    endif
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(1,num_elmnt)
!      elmnts(1,num_elmnt) = elmnts(5,num_elmnt)
!      elmnts(5,num_elmnt) = ielm
!
!      ! assumes to switch the others as well
!      ielm = elmnts(2,num_elmnt)
!      elmnts(2,num_elmnt) = elmnts(6,num_elmnt)
!      elmnts(6,num_elmnt) = ielm
!
!      ielm = elmnts(3,num_elmnt)
!      elmnts(3,num_elmnt) = elmnts(7,num_elmnt)
!      elmnts(7,num_elmnt) = ielm
!
!      ielm = elmnts(4,num_elmnt)
!      elmnts(4,num_elmnt) = elmnts(8,num_elmnt)
!      elmnts(8,num_elmnt) = ielm
!
!    endif
!    ! makes sure bottom - top is o.k.
!    ! point 2 (0,1,0) vs point 6 (0,1,1)
!    inorder = .true.
!    if( nodes_coords(3,elmnts(2,num_elmnt)) > nodes_coords(3,elmnts(6,num_elmnt)) ) then
!      print*,num_elmnt,'z2-6 :',nodes_coords(3,elmnts(2,num_elmnt)),nodes_coords(3,elmnts(6,num_elmnt))
!      inorder = .false.
!    endif    
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(2,num_elmnt)
!      elmnts(2,num_elmnt) = elmnts(6,num_elmnt)
!      elmnts(6,num_elmnt) = ielm
!    endif
!    
!    ! point 3 (1,1,0) vs point 7 (1,1,1)
!    inorder = .true.
!    if( nodes_coords(3,elmnts(3,num_elmnt)) > nodes_coords(3,elmnts(7,num_elmnt)) ) then
!      print*,num_elmnt,'z3-7 :',nodes_coords(3,elmnts(3,num_elmnt)),nodes_coords(3,elmnts(7,num_elmnt))
!      inorder = .false.
!    endif    
!    if( inorder .eqv. .false. ) then    
!      ielm = elmnts(3,num_elmnt)
!      elmnts(3,num_elmnt) = elmnts(7,num_elmnt)
!      elmnts(7,num_elmnt) = ielm
!    endif
!    
!    ! point 4 (1,0,0) vs point 8 (1,0,1)
!    inorder = .true.
!    if( nodes_coords(3,elmnts(4,num_elmnt)) > nodes_coords(3,elmnts(8,num_elmnt)) ) then
!      print*,num_elmnt,'z4-8 :',nodes_coords(3,elmnts(4,num_elmnt)),nodes_coords(3,elmnts(8,num_elmnt))
!      inorder = .false.
!    endif    
!    if( inorder .eqv. .false. ) then    
!      ielm = elmnts(4,num_elmnt)
!      elmnts(4,num_elmnt) = elmnts(8,num_elmnt)
!      elmnts(8,num_elmnt) = ielm
!    endif
!
!    ! clock-wise order?
!    ! point 1 (0,0,0) vs point 3 (1,1,0)
!    inorder = .true.
!    if( nodes_coords(1,elmnts(1,num_elmnt)) > nodes_coords(1,elmnts(3,num_elmnt)) ) then
!      print*,num_elmnt,'x1-3 :',nodes_coords(1,elmnts(1,num_elmnt)),nodes_coords(1,elmnts(3,num_elmnt))
!      inorder = .false.
!    endif
!    if( nodes_coords(2,elmnts(1,num_elmnt)) > nodes_coords(2,elmnts(3,num_elmnt)) ) then
!      print*,num_elmnt,'y1-3 :',nodes_coords(2,elmnts(1,num_elmnt)),nodes_coords(2,elmnts(3,num_elmnt))
!      inorder = .false.
!    endif
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(1,num_elmnt)
!      elmnts(1,num_elmnt) = elmnts(3,num_elmnt)
!      elmnts(3,num_elmnt) = ielm
!    endif
!
!    ! point 2 (0,1,0) vs point 4 (1,0,0)
!    inorder = .true.
!    if( nodes_coords(1,elmnts(2,num_elmnt)) > nodes_coords(1,elmnts(4,num_elmnt)) ) then
!      print*,num_elmnt,'x2-4 :',nodes_coords(1,elmnts(2,num_elmnt)),nodes_coords(1,elmnts(4,num_elmnt))
!      inorder = .false.
!    endif
!    if( nodes_coords(2,elmnts(2,num_elmnt)) < nodes_coords(2,elmnts(4,num_elmnt)) ) then
!      print*,num_elmnt,'y2-4 :',nodes_coords(2,elmnts(2,num_elmnt)),nodes_coords(2,elmnts(4,num_elmnt))
!      inorder = .false.
!    endif
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(2,num_elmnt)
!      elmnts(2,num_elmnt) = elmnts(4,num_elmnt)
!      elmnts(4,num_elmnt) = ielm
!    endif
!
!    ! point 5 (0,0,1) vs point 7 (1,1,1)
!    inorder = .true.
!    if( nodes_coords(1,elmnts(5,num_elmnt)) > nodes_coords(1,elmnts(7,num_elmnt)) ) then
!      print*,num_elmnt,'x5-7 :',nodes_coords(1,elmnts(5,num_elmnt)),nodes_coords(1,elmnts(7,num_elmnt))
!      inorder = .false.
!    endif
!    if( nodes_coords(2,elmnts(5,num_elmnt)) > nodes_coords(2,elmnts(7,num_elmnt)) ) then
!      print*,num_elmnt,'y5-7 :',nodes_coords(2,elmnts(5,num_elmnt)),nodes_coords(2,elmnts(7,num_elmnt))
!      inorder = .false.
!    endif
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(5,num_elmnt)
!      elmnts(5,num_elmnt) = elmnts(7,num_elmnt)
!      elmnts(7,num_elmnt) = ielm
!    endif
!
!    ! point 6 (0,1,1) vs point 8 (1,0,1)
!    inorder = .true.
!    if( nodes_coords(1,elmnts(6,num_elmnt)) > nodes_coords(1,elmnts(8,num_elmnt)) ) then
!      print*,num_elmnt,'x6-8 :',nodes_coords(1,elmnts(6,num_elmnt)),nodes_coords(1,elmnts(8,num_elmnt))
!      inorder = .false.
!    endif
!    if( nodes_coords(2,elmnts(6,num_elmnt)) < nodes_coords(2,elmnts(8,num_elmnt)) ) then
!      print*,num_elmnt,'y6-8 :',nodes_coords(2,elmnts(6,num_elmnt)),nodes_coords(2,elmnts(8,num_elmnt))
!      inorder = .false.
!    endif
!    if( inorder .eqv. .false. ) then
!      ielm = elmnts(6,num_elmnt)
!      elmnts(6,num_elmnt) = elmnts(8,num_elmnt)
!      elmnts(8,num_elmnt) = ielm
!    endif
!
! or    
!    if( .false. ) then
!      ! trys to order points in increasing z direction first, then y and x
!      inorder = .false.
!      do while (inorder .eqv. .false.)
!        inorder = .true.       
!        do i=1,8              
!          ! If z needs to be swapped, do so 
!          if (nodes_coords(3,elmnts(i,num_elmnt)) > nodes_coords(3,elmnts(i+1,num_elmnt)) )then
!            i_temp = elmnts(i,num_elmnt)
!            elmnts(i,num_elmnt) = elmnts(i+1,num_elmnt)
!            elmnts(i+1,num_elmnt) = i_temp
!            inorder = .false.
!            exit
!          endif         
!          ! Check Equilivant Points and swap those on Y
!          if (nodes_coords(3,elmnts(i,num_elmnt)) == nodes_coords(3,elmnts(i+1,num_elmnt))) then
!            if (nodes_coords(2,elmnts(i,num_elmnt)) > nodes_coords(2,elmnts(i+1,num_elmnt)) ) then
!              i_temp = elmnts(i,num_elmnt)
!              elmnts(i,num_elmnt) = elmnts(i+1,num_elmnt)
!              elmnts(i+1,num_elmnt) = i_temp
!              inorder = .false.
!              exit
!            endif
!          endif
!          ! Check Equilivant Points and swap those on X
!          if (nodes_coords(3,elmnts(i,num_elmnt)) == nodes_coords(3,elmnts(i+1,num_elmnt))) then
!            if (nodes_coords(2,elmnts(i,num_elmnt)) == nodes_coords(2,elmnts(i+1,num_elmnt)) ) then
!              if (nodes_coords(1,elmnts(i,num_elmnt)) > nodes_coords(1,elmnts(i+1,num_elmnt)) )then
!                i_temp = elmnts(i,num_elmnt)
!                elmnts(i,num_elmnt) = elmnts(i+1,num_elmnt)
!                elmnts(i+1,num_elmnt) = i_temp
!                inorder = .false.
!                exit
!              endif
!            endif 
!          endif
!        enddo
!      enddo    
!      ! respect anti-clockwise ordering bottom face
!      i_temp = elmnts(3,num_elmnt)
!      elmnts(3,num_elmnt) = elmnts(4,num_elmnt)   
!      elmnts(4,num_elmnt) = i_temp
!      ! respect anti-clockwise ordering top face
!      i_temp = elmnts(7,num_elmnt)
!      elmnts(7,num_elmnt) = elmnts(8,num_elmnt)   
!      elmnts(8,num_elmnt) = i_temp        
!      if( nodes_coords(1,elmnts(1,num_elmnt)) > nodes_coords(1,elmnts(2,num_elmnt)) ) then
!        print*,'elem:',num_elmnt
!        stop 'error sorting x'
!      endif
!      if( nodes_coords(2,elmnts(1,num_elmnt)) > nodes_coords(2,elmnts(4,num_elmnt)) ) then
!        print*,'elem:',num_elmnt
!        stop 'error sorting y'
!     endif
!      if( nodes_coords(3,elmnts(1,num_elmnt)) > nodes_coords(3,elmnts(5,num_elmnt)) ) then
!        print*,'elem:',num_elmnt
!        stop 'error sorting z'
!      endif
!    endif

    call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
                      etaxstore,etaystore,etazstore, &
                      gammaxstore,gammaystore,gammazstore,jacobianstore, &
                      xstore,ystore,zstore, &
                      xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

  enddo

end subroutine create_regions_mesh_ext_mesh_setup_jacobian

!
!----
!

subroutine create_regions_mesh_ext_mesh_determine_kappamu(nspec,mat_ext_mesh,nelmnts_ext_mesh,&
                        materials_ext_mesh,nmat_ext_mesh,&
                        undef_mat_prop,nundefMat_ext_mesh,&
                        rhostore,kappastore,mustore,vpstore,vsstore,&
                        iflag_attenuation_store,rho_vp,rho_vs)

  implicit none

  include 'constants.h'

! number of spectral elements in each block
  integer :: nspec

! external mesh
  integer :: nelmnts_ext_mesh
  integer :: nmat_ext_mesh,nundefMat_ext_mesh 

  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh
  double precision, dimension(5,nmat_ext_mesh) :: materials_ext_mesh  
  character (len=30), dimension(5,nundefMat_ext_mesh):: undef_mat_prop

! for model density
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore, &
                                        kappastore,mustore,vpstore,vsstore 
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rho_vp,rho_vs

! attenuation 
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: iflag_attenuation_store

! local parameters
  integer :: ispec,i,j,k,iundef
  integer  :: iflag,flag_below,flag_above

! !  Piero, read bedrock file
!  allocate(ibedrock(NX_TOPO_ANT,NY_TOPO_ANT))              
!  if(myrank == 0) then
!      call read_bedrock_file(ibedrock)
!  !    write(IMAIN,*)
!  !    write(IMAIN,*) 'regional bedrock file read ranges in m from ',minval(ibedrock),' to ',maxval(ibedrock)
!  !    write(IMAIN,*)
!   endif
!  ! broadcast the information read on the master to the nodes
!  ! call MPI_BCAST(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT,MPI_REAL,0,MPI_COMM_WORLD,ier)
! call bcast_all_cr(ibedrock,NX_TOPO_ANT*NY_TOPO_ANT)
  
! kappastore and mustore
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
! check if the material is known or unknown
           if (mat_ext_mesh(1,ispec) > 0) then
              ! materials_ext_mesh format: #index1 = rho #index2 = vp #index3 = vs #index4 = Q_flag #index5 = 0 
              rhostore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(1,ispec))
              vpstore(i,j,k,ispec) = materials_ext_mesh(2,mat_ext_mesh(1,ispec))
              vsstore(i,j,k,ispec) = materials_ext_mesh(3,mat_ext_mesh(1,ispec))
              iflag_attenuation_store(i,j,k,ispec) = materials_ext_mesh(4,mat_ext_mesh(1,ispec))                            
              !change for piero :
              !if(mat_ext_mesh(1,ispec) == 1) then
              !   iflag_attenuation_store(i,j,k,ispec) = 1
              !else
              !   iflag_attenuation_store(i,j,k,ispec) = 2
              !endif
           else if (mat_ext_mesh(2,ispec) == 1) then
              stop 'material: interface not implemented yet'
              
              do iundef = 1,nundefMat_ext_mesh 
                 if(trim(undef_mat_prop(2,iundef)) == 'interface') then
                    read(undef_mat_prop(4,iundef),'(1i3)') flag_below
                    read(undef_mat_prop(5,iundef),'(1i3)') flag_above
                 endif
              enddo
              !call interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)
              iflag = 1
              rhostore(i,j,k,ispec) = materials_ext_mesh(1,iflag)
              vpstore(i,j,k,ispec) = materials_ext_mesh(2,iflag)
              vsstore(i,j,k,ispec) = materials_ext_mesh(3,iflag)
              iflag_attenuation_store(i,j,k,ispec) = materials_ext_mesh(4,iflag)
              !change for piero :
              !  if(iflag == 1) then
              !     iflag_attenuation_store(i,j,k,ispec) = 1
              !  else
              !     iflag_attenuation_store(i,j,k,ispec) = 2
              !  endif
             else
              stop 'material: tomography not implemented yet'
             ! call tomography()
           end if

           kappastore(i,j,k,ispec) = rhostore(i,j,k,ispec)* &
                ( vpstore(i,j,k,ispec)*vpstore(i,j,k,ispec) &
                - FOUR_THIRDS*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec) )
                
           mustore(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec)
           
           ! Stacey, a completer par la suite  
           rho_vp(i,j,k,ispec) = rhostore(i,j,k,ispec)*vpstore(i,j,k,ispec)
           rho_vs(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)
           !end pll

        enddo
      enddo
    enddo
    !print*,myrank,'ispec:',ispec,'rho:',rhostore(1,1,1,ispec),'vp:',vpstore(1,1,1,ispec),'vs:',vsstore(1,1,1,ispec)    
  enddo

! !! DK DK store the position of the six stations to be able to
! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!    utm_x_station(1) =  783500.6250000d0
!    utm_y_station(1) = -11828.7519531d0

!    utm_x_station(2) =  853644.5000000d0
!    utm_y_station(2) = -114.0138092d0

!    utm_x_station(3) = 863406.0000000d0
!    utm_y_station(3) = -53736.1640625d0

!    utm_x_station(4) =   823398.8125000d0
!    utm_y_station(4) = 29847.4511719d0

!    utm_x_station(5) = 863545.3750000d0
!    utm_y_station(5) = 19669.6621094d0

!    utm_x_station(6) =  817099.3750000d0
!    utm_y_station(6) = -24430.2871094d0

!  print*,myrank,'après store the position of the six stations'
!  call flush(6)

!  print*, myrank,minval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


! print*, myrank,maxval(nodes_coords_ext_mesh(1,:))
!  call flush(6)


!  do ispec = 1, nspec

!     zmesh = zstore(2,2,2,ispec)

!    ! if(doubling_index == IFLAG_ONE_LAYER_TOPOGRAPHY) then
!     if(any(ibelm_top == ispec)) then
!     doubling_value_found_for_Piero = IFLAG_ONE_LAYER_TOPOGRAPHY
       
!     else if(zmesh < Z_23p4km) then
!        doubling_value_found_for_Piero = IFLAG_MANTLE_BELOW_23p4km
       
!     else if(zmesh < Z_14km) then
!        doubling_value_found_for_Piero = IFLAG_14km_to_23p4km
       
!     else
!        doubling_value_found_for_Piero = IFLAG_BEDROCK_down_to_14km
!     endif
!    idoubling(ispec) = doubling_value_found_for_Piero

!     do k = 1, NGLLZ
!       do j = 1, NGLLY
!         do i = 1, NGLLX

           
!            if(idoubling(ispec) == IFLAG_ONE_LAYER_TOPOGRAPHY .or. &
!               idoubling(ispec) == IFLAG_BEDROCK_down_to_14km) then
              
!               ! since we have suppressed UTM projection for Piero Basini, UTMx is the same as long
!               ! and UTMy is the same as lat
!               long = xstore(i,j,k,ispec)
!               lat = ystore(i,j,k,ispec)
              
!               ! get coordinate of corner in model
!               icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
!               icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1
              
!               ! avoid edge effects and extend with identical point if outside model
!               if(icornerlong < 1) icornerlong = 1
!               if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
!               if(icornerlat < 1) icornerlat = 1
!               if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1
              
!               ! compute coordinates of corner
!               long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
!               lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO
                   
!               ! compute ratio for interpolation
!               ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
!               ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO
                   
!               ! avoid edge effects
!               if(ratio_xi < 0.) ratio_xi = 0.
!               if(ratio_xi > 1.) ratio_xi = 1.
!               if(ratio_eta < 0.) ratio_eta = 0.
!               if(ratio_eta > 1.) ratio_eta = 1.
                   
!               ! interpolate elevation at current point
!               elevation_bedrock = &
!                    ibedrock(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!                    ibedrock(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!                    ibedrock(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta
                   
!               !! DK DK exclude circles around each station to make sure they are on the bedrock
!               !! DK DK and not in the ice
!               is_around_a_station = .false.
!               do istation = 1,NUMBER_OF_STATIONS
!                  if(sqrt((long - utm_x_station(istation))**2 + (lat - utm_y_station(istation))**2) < RADIUS_TO_EXCLUDE) then
!                     is_around_a_station = .true.
!                     exit
!                  endif
!               enddo
              
!               ! define elastic parameters in the model
              
!               ! we are above the bedrock interface i.e. in the ice, and not too close to a station
!               if(zmesh >= elevation_bedrock .and. .not. is_around_a_station) then
!                  vp = 3800.d0
!                  vs = 1900.d0
!                  rho = 900.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_ICE
                 
!                  ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!               else
!                  vp = 5800.d0
!                  vs = 3200.d0
!                  rho = 2600.d0
!                  iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
!               endif
              
!            else if(idoubling(ispec) == IFLAG_14km_to_23p4km) then
!               vp = 6800.d0
!               vs = 3900.d0
!               rho = 2900.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            else if(idoubling(ispec) == IFLAG_MANTLE_BELOW_23p4km) then
!               vp = 8100.d0
!               vs = 4480.d0
!               rho = 3380.d0
!               iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
              
!            endif
           
!                 !pll  8/06
!                     if(CUSTOM_REAL == SIZE_REAL) then
!                        rhostore(i,j,k,ispec) = sngl(rho)
!                        vpstore(i,j,k,ispec) = sngl(vp)
!                        vsstore(i,j,k,ispec) = sngl(vs)
!                     else
!                        rhostore(i,j,k,ispec) = rho
!                        vpstore(i,j,k,ispec) = vp
!                        vsstore(i,j,k,ispec) = vs
!                     end if
                
!                 kappastore(i,j,k,ispec) = rhostore(i,j,k,ispec)*(vpstore(i,j,k,ispec)*vpstore(i,j,k,ispec) - &
!                      4.d0*vsstore(i,j,k,ispec)*vsstore(i,j,k,ispec)/3.d0)
!                 mustore(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)*&
!                      vsstore(i,j,k,ispec)
           
!                 ! Stacey, a completer par la suite  
!                 rho_vp(i,j,k,ispec) = rhostore(i,j,k,ispec)*vpstore(i,j,k,ispec)
!                 rho_vs(i,j,k,ispec) = rhostore(i,j,k,ispec)*vsstore(i,j,k,ispec)
!                 !end pll
                
!                 !      kappastore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                 !       (materials_ext_mesh(2,mat_ext_mesh(ispec))*materials_ext_mesh(2,mat_ext_mesh(ispec)) - &
!                 !        4.d0*materials_ext_mesh(3,mat_ext_mesh(ispec))*materials_ext_mesh(3,mat_ext_mesh(ispec))/3.d0)
!                 !      mustore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(ispec))* &
!                                                         materials_ext_mesh(3,mat_ext_mesh(ispec))*&
!                 !  x    materials_ext_mesh(3,mat_ext_mesh(ispec))
!              enddo
!           enddo
!        enddo
!     enddo

end subroutine create_regions_mesh_ext_mesh_determine_kappamu

!
!----
!

subroutine create_regions_mesh_ext_mesh_setup_global_indexing(ibool, &
                            xstore,ystore,zstore,nspec,nglob,npointot, &
                            nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! creates global indexing array ibool

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer :: nspec,nglob,npointot,myrank

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

! local parameters
! variables for creating array ibool 
  double precision, dimension(:), allocatable :: xp,yp,zp
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg

  integer :: ieoff,ilocnum,ier
  integer :: i,j,k,ispec

! allocate memory for arrays
  allocate(locval(npointot), &
          ifseg(npointot), &
          xp(npointot), &
          yp(npointot), &
          zp(npointot),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! creates temporary global point arrays
  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

  do ispec=1,nspec
    ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
    ilocnum = 0
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ilocnum = ilocnum + 1
          xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
          yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
          zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

! gets ibool indexing from local (GLL points) to global points
  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot, &
       minval(nodes_coords_ext_mesh(1,:)),maxval(nodes_coords_ext_mesh(1,:)))

!- we can create a new indirect addressing to reduce cache misses
  call get_global_indirect_addressing(nspec,nglob,ibool)

!cleanup
  deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'

end subroutine create_regions_mesh_ext_mesh_setup_global_indexing

!
!----
!

subroutine create_regions_mesh_ext_mesh_create_mass_matrix(nglob,rmass,&
          nspec,wxgll,wygll,wzgll,ibool,jacobianstore,rhostore)

! returns precomputed mass matrix in rmass array

  implicit none

  include 'constants.h'

! number of spectral elements in each block
  integer :: nglob,nspec

! mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass

! Gauss-Lobatto-Legendre weights of integration
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: jacobianstore,rhostore

! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,i,j,k,iglobnum

! creates mass matrix  
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          weight=wxgll(i)*wygll(j)*wzgll(k)
          iglobnum=ibool(i,j,k,ispec)

          jacobianl=jacobianstore(i,j,k,ispec)

! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglobnum) = rmass(iglobnum) + &
                sngl((dble(rhostore(i,j,k,ispec)))  * dble(jacobianl) * weight)
          else
             rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
          endif

        enddo
      enddo
    enddo
  enddo  


end subroutine create_regions_mesh_ext_mesh_create_mass_matrix

!
!----
!

subroutine create_regions_mesh_ext_mesh_setup_absorbing_bound(myrank,nspec,nglob,&
                            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                            wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top, &
                            normal_top,jacobian2D_top, &
                            absorbing_boundary_normal,absorbing_boundary_jacobian2D, &
                            absorbing_boundary_ijk,absorbing_boundary_ispec, &
                            num_absorbing_boundary_faces)

! determines absorbing boundaries/free-surface, 2D jacobians, face normals for Stacey conditions
  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer :: myrank,nspec,nglob

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
!  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
! global point locations          
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)

! 2D shape functions derivatives and weights
  double precision :: dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ),dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
          dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY),dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz

! data from the external mesh
  integer :: nnodes_ext_mesh !,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
!  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! absorbing boundaries
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
  ! element indices containing a boundary
  integer, dimension(nspec2D_xmin)  :: ibelm_xmin  
  integer, dimension(nspec2D_xmax)  :: ibelm_xmax
  integer, dimension(nspec2D_ymin)  :: ibelm_ymin
  integer, dimension(nspec2D_ymax)  :: ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM)  :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP)  :: ibelm_top

  ! corner node indices of boundary faces coming from CUBIT
  integer, dimension(4,nspec2D_xmin)  :: nodes_ibelm_xmin  
  integer, dimension(4,nspec2D_xmax)  :: nodes_ibelm_xmax
  integer, dimension(4,nspec2D_ymin)  :: nodes_ibelm_ymin
  integer, dimension(4,nspec2D_ymax)  :: nodes_ibelm_ymax
  integer, dimension(4,NSPEC2D_BOTTOM)  :: nodes_ibelm_bottom
  integer, dimension(4,NSPEC2D_TOP)  :: nodes_ibelm_top

  ! local indices i,j,k of all GLL points on an absorbing boundary in the element, 
  ! defines all gll points located on the absorbing surfaces
!  integer :: ibelm_gll_xmin(3,NGLLY,NGLLZ,nspec2D_xmin),ibelm_gll_xmax(3,NGLLY,NGLLZ,nspec2D_xmax), &
!            ibelm_gll_ymin(3,NGLLX,NGLLZ,nspec2D_ymin),ibelm_gll_ymax(3,NGLLX,NGLLZ,nspec2D_ymax), &
!            ibelm_gll_bottom(3,NGLLY,NGLLY,nspec2D_bottom),ibelm_gll_top(3,NGLLY,NGLLY,nspec2D_top)

! overlap indices for elements at corners and edges with more than one aborbing boundary face
!  integer  :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
!  integer :: nimin(2,NSPEC2DMAX_YMIN_YMAX),nimax(2,NSPEC2DMAX_YMIN_YMAX), &
!            njmin(2,NSPEC2DMAX_XMIN_XMAX),njmax(2,NSPEC2DMAX_XMIN_XMAX), &
!            nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX),nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX)

  ! 2-D jacobians and normals
!  real(kind=CUSTOM_REAL) :: jacobian2D_xmin(NGLLY,NGLLZ,nspec2D_xmin),&
!                jacobian2D_xmax(NGLLY,NGLLZ,nspec2D_xmax), &
!                 jacobian2D_ymin(NGLLX,NGLLZ,nspec2D_ymin), &
!                 jacobian2D_ymax(NGLLX,NGLLZ,nspec2D_ymax),&
!                 jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM),&
  real(kind=CUSTOM_REAL):: jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  integer :: num_absorbing_boundary_faces
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_normal
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_jacobian2D
  integer, dimension(3,NGLLSQUARE,num_absorbing_boundary_faces) :: absorbing_boundary_ijk
  integer, dimension(num_absorbing_boundary_faces) :: absorbing_boundary_ispec    
  
  ! normals for all GLL points on boundaries
!  real(kind=CUSTOM_REAL) :: normal_xmin(NDIM,NGLLY,NGLLZ,nspec2D_xmin),&
!           normal_xmax(NDIM,NGLLY,NGLLZ,nspec2D_xmax), &
!           normal_ymin(NDIM,NGLLX,NGLLZ,nspec2D_ymin), &
!           normal_ymax(NDIM,NGLLX,NGLLZ,nspec2D_ymax), &
!           normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM),
  real(kind=CUSTOM_REAL) :: normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)  
          
! local parameters
! pll 
  logical, dimension(:,:),allocatable :: iboun  

  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL) :: jacobian2D_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  integer:: ijk_face(3,NGLLX,NGLLY)
  
  ! corner locations for faces
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xcoord_iboun,ycoord_iboun,zcoord_iboun
  
  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  integer  :: ispec,ispec2D,icorner,ier,iabs,iface,igll,i,j
  
! allocate temporary flag array
  allocate(iboun(6,nspec), &
          xcoord_iboun(NGNOD2D,6,nspec), &
          ycoord_iboun(NGNOD2D,6,nspec), &
          zcoord_iboun(NGNOD2D,6,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  
! sets flag in array iboun for elements with an absorbing boundary faces
  iboun(:,:) = .false. 

! abs face counter  
  iabs = 0
  
  ! xmin   
  do ispec2D = 1, nspec2D_xmin 
    ! sets element 
    ispec = ibelm_xmin(ispec2D)
     
    !if(myrank == 0 ) print*,'xmin:',ispec2D,ispec
    
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmin(icorner,ispec2D))
      !print*,'corner look:',icorner,xcoord(icorner),ycoord(icorner),zcoord(icorner)
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    iboun(iface,ispec) = .true. 

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)    
    
    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2D_face,normal_face,NGLLX,NGLLZ)                              

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (-1,0,0)    
    !if( myrank == 0 ) then
    !if( abs(normal_face(1,1,1) + 1.0 ) > 0.1 ) then
    !  print*,'error normal xmin',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !  stop
    !endif    
    !if( abs(xstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) - 0.0) > 0.1 ) &
    !  print*,'error element xmin:',ispec,xstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))
    !endif
                            
    ! sets face infos
    iabs = iabs + 1
    absorbing_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
        absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo        

  enddo ! nspec2D_xmin
 
  ! xmax
  do ispec2D = 1, nspec2D_xmax 
    ! sets element
    ispec = ibelm_xmax(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_xmax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_xmax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_xmax(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)
    
    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2D_face,normal_face,NGLLX,NGLLZ)                              

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (1,0,0)    
    !if( abs(normal_face(1,1,1) - 1.0 ) > 0.1 ) then
    !  print*,'error normal xmax',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !endif    
    !if( abs(xstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) - 134000.0) > 0.1 ) &
    !  print*,'error element xmax:',ispec,xstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))

    ! sets face infos
    iabs = iabs + 1
    absorbing_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
        absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo            
    
  enddo

  ! ymin
  do ispec2D = 1, nspec2D_ymin 
    ! sets element 
    ispec = ibelm_ymin(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymin(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymin(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymin(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2D_face,normal_face,NGLLY,NGLLZ)                              

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (0,-1,0)    
    !if( abs(normal_face(2,1,1) + 1.0 ) > 0.1 ) then
    !  print*,'error normal ymin',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !endif    
    !if( abs(ystore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) - 0.0) > 0.1 ) &
    !  print*,'error element ymin:',ispec,ystore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))

    ! sets face infos
    iabs = iabs + 1
    absorbing_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLY
        igll = igll+1
        absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
        absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo        
                                  
  enddo

  ! ymax
  do ispec2D = 1, nspec2D_ymax 
    ! sets element 
    ispec = ibelm_ymax(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_ymax(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_ymax(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_ymax(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLY,NGLLZ)                              

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2D_face,normal_face,NGLLY,NGLLZ) 

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLZ
      do i=1,NGLLY
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (0,1,0)    
    !if( abs(normal_face(2,1,1) - 1.0 ) > 0.1 ) then
    !  print*,'error normal ymax',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !endif    
    !if( abs(ystore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) - 134000.0) > 0.1 ) &
    !  print*,'error element ymax:',ispec,ystore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))

    ! sets face infos
    iabs = iabs + 1
    absorbing_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
        absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo
    
  enddo
  
  ! bottom
  do ispec2D = 1, NSPEC2D_BOTTOM
    ! sets element 
    ispec = ibelm_bottom(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_bottom(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_bottom(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_bottom(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )   
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)
    
    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2D_face,normal_face,NGLLX,NGLLY) 

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (0,0,-1)    
    !if( abs(normal_face(3,1,1) + 1.0 ) > 0.1 ) then
    !  print*,'error normal bottom',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !endif    
    !if( abs(zstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) + 60000.0) > 0.1 ) &
    !  print*,'error element bottom:',ispec,zstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))

    ! sets face infos
    iabs = iabs + 1
    absorbing_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
        absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo    
    
  enddo
  
  ! top
  do ispec2D = 1, NSPEC2D_TOP
    ! sets element 
    ispec = ibelm_top(ispec2D)
     
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_top(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_top(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_top(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface )
    iboun(iface,ispec) = .true. 
                              
    ! ijk indices of GLL points on face
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2D_face,normal_face,NGLLX,NGLLY) 

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob, &
                                      xstore_dummy,ystore_dummy,zstore_dummy, &
                                      normal_face(:,i,j) )
      enddo
    enddo

    !daniel
    ! checks: layered halfspace  normals
    ! for boundary on xmin, outward direction must be (0,0,1)    
    !if( abs(normal_face(3,1,1) - 1.0 ) > 0.1 ) then
    !  print*,'error normal top',myrank,ispec
    !  print*,sngl(normal_face(:,1,1))
    !endif    
    !if( abs(zstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec)) - 0.0) > 0.1 ) &
    !  print*,'error element top:',ispec,zstore_dummy(ibool(ijk_face(1,2,2),ijk_face(2,2,2),ijk_face(3,2,2),ispec))

    ! store for free surface
    jacobian2D_top(:,:,ispec2D) = jacobian2D_face(:,:)
    normal_top(:,:,:,ispec2D) = normal_face(:,:,:)  

    ! store for absorbing boundaries
    if( ABSORB_FREE_SURFACE ) then
      ! sets face infos
      iabs = iabs + 1
      absorbing_boundary_ispec(iabs) = ispec      
      
      ! gll points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j=1,NGLLY
        do i=1,NGLLX
          igll = igll+1
          absorbing_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
          absorbing_boundary_jacobian2D(igll,iabs) = jacobian2D_face(i,j)
          absorbing_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
        enddo
      enddo
    endif
  enddo
  
  if( iabs /= num_absorbing_boundary_faces ) then
    print*,'error number of absorbing faces:',iabs,num_absorbing_boundary_faces
    stop 'error number of absorbing faces'
  endif

  call sum_all_i(num_absorbing_boundary_faces,iabs)
  if( myrank == 0 ) then
    write(IMAIN,*) '     absorbing boundary:'
    write(IMAIN,*) '     total number of faces = ',iabs
    if( ABSORB_FREE_SURFACE ) then
    write(IMAIN,*) 'absorbing boundary includes free surface'
    endif
  endif

!obsolete...
! calculates 2D jacobians and normals for each GLL point on face
!  call get_jacobian_boundaries(myrank,iboun,nspec, &
!                   xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&  
!                   dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
!                   wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                   
!                   ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
!                   xcoord_iboun,ycoord_iboun,zcoord_iboun, &
!                   nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
!                   jacobian2D_xmin,jacobian2D_xmax, &
!                   jacobian2D_ymin,jacobian2D_ymax, &
!                   jacobian2D_bottom,jacobian2D_top, &
!                   normal_xmin,normal_xmax, &
!                   normal_ymin,normal_ymax, &
!                   normal_bottom,normal_top, &
!                   NSPEC2D_BOTTOM,NSPEC2D_TOP)
! obsolete... arrays not used anymore...  
! Stacey put back
!  call get_absorb_ext_mesh(myrank,iboun,nspec, &
!       nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
!       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

end subroutine create_regions_mesh_ext_mesh_setup_absorbing_bound

!
!----
!

subroutine create_regions_mesh_ext_mesh_prepare_MPI_interfaces(nglob,nspec,ibool, &
                                    nelmnts_ext_mesh,elmnts_ext_mesh, &
                                    my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                                    ibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh, &
                                    ninterface_ext_mesh,max_interface_size_ext_mesh, &
                                    xstore_dummy,ystore_dummy,zstore_dummy)

! sets up the MPI interface for communication between partitions

  implicit none

  include "constants.h"

  integer :: nglob,nspec

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  
  integer :: nelmnts_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh
  
  integer :: ninterface_ext_mesh,max_interface_size_ext_mesh
  
  integer, dimension(ninterface_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,ninterface_ext_mesh) :: my_interfaces_ext_mesh
  
  integer, dimension(ninterface_ext_mesh) :: nibool_interfaces_ext_mesh  
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,ninterface_ext_mesh) :: ibool_interfaces_ext_mesh
  
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  
!local parameters
  double precision, dimension(:), allocatable :: xp,yp,zp
  double precision, dimension(:), allocatable :: work_ext_mesh

  integer, dimension(:), allocatable :: locval !,iglob
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh_true

! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh,ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh
  integer, dimension(:), allocatable :: ibool_interface_ext_mesh_dummy

  logical, dimension(:), allocatable :: ifseg

  integer :: iinterface,ilocnum
  

! get global indices for MPI interfaces between different partitions
  call prepare_assemble_MPI (nelmnts_ext_mesh,ibool, &
                            elmnts_ext_mesh, ESIZE, &
                            nglob, &
                            ninterface_ext_mesh, max_interface_size_ext_mesh, &
                            my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                            ibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh &
                            )

  allocate(nibool_interfaces_ext_mesh_true(ninterface_ext_mesh))

! sort ibool comm buffers lexicographically  
  do iinterface = 1, ninterface_ext_mesh

    allocate(xp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(yp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(zp(nibool_interfaces_ext_mesh(iinterface)))
    allocate(locval(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ifseg(nibool_interfaces_ext_mesh(iinterface)))
    allocate(reorder_interface_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ibool_interface_ext_mesh_dummy(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ind_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(ninseg_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(iwork_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))
    allocate(work_ext_mesh(nibool_interfaces_ext_mesh(iinterface)))

    do ilocnum = 1, nibool_interfaces_ext_mesh(iinterface)
      xp(ilocnum) = xstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      yp(ilocnum) = ystore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      zp(ilocnum) = zstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
    enddo

    call sort_array_coordinates(nibool_interfaces_ext_mesh(iinterface),xp,yp,zp, &
         ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         reorder_interface_ext_mesh,locval,ifseg,nibool_interfaces_ext_mesh_true(iinterface), &
         ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh,work_ext_mesh)

    deallocate(xp)
    deallocate(yp)
    deallocate(zp)
    deallocate(locval)
    deallocate(ifseg)
    deallocate(reorder_interface_ext_mesh)
    deallocate(ibool_interface_ext_mesh_dummy)
    deallocate(ind_ext_mesh)
    deallocate(ninseg_ext_mesh)
    deallocate(iwork_ext_mesh)
    deallocate(work_ext_mesh)

  enddo

end subroutine create_regions_mesh_ext_mesh_prepare_MPI_interfaces

!pll
! subroutine interface(iflag,flag_below,flag_above,ispec,nspec,i,j,k,xstore,ystore,zstore,ibedrock)

! implicit none

! include "constants.h"

! integer :: iflag,flag_below,flag_above
! integer :: ispec,nspec
! integer :: i,j,k
! double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
! real(kind=CUSTOM_REAL), dimension(NX_TOPO_ANT,NY_TOPO_ANT) :: ibedrock
! integer, parameter :: NUMBER_OF_STATIONS = 1
! double precision, parameter :: RADIUS_TO_EXCLUDE = 250.d0
! double precision, dimension(NUMBER_OF_STATIONS) :: utm_x_station,utm_y_station

! !-------------------

! !for Piero
! logical :: is_around_a_station
! integer :: istation

! ! store bedrock values
! integer ::  icornerlat,icornerlong
! double precision ::  lat,long,elevation_bedrock
! double precision ::  lat_corner,long_corner,ratio_xi,ratio_eta


! !! DK DK store the position of the six stations to be able to
! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!    utm_x_station(1) =  783500.6250000d0
!    utm_y_station(1) = -11828.7519531d0

!    utm_x_station(2) =  853644.5000000d0
!    utm_y_station(2) = -114.0138092d0

!    utm_x_station(3) = 863406.0000000d0
!    utm_y_station(3) = -53736.1640625d0

!    utm_x_station(4) =   823398.8125000d0
!    utm_y_station(4) = 29847.4511719d0

!    utm_x_station(5) = 863545.3750000d0
!    utm_y_station(5) = 19669.6621094d0

!    utm_x_station(6) =  817099.3750000d0
!    utm_y_station(6) = -24430.2871094d0

! ! since we have suppressed UTM projection for Piero Basini, UTMx is the same as long
! ! and UTMy is the same as lat
!     long = xstore(i,j,k,ispec)
!     lat =  ystore(i,j,k,ispec)

! ! get coordinate of corner in model
!     icornerlong = int((long - ORIG_LONG_TOPO_ANT) / DEGREES_PER_CELL_TOPO_ANT) + 1
!     icornerlat = int((lat - ORIG_LAT_TOPO_ANT) / DEGREES_PER_CELL_TOPO_ANT) + 1

! ! avoid edge effects and extend with identical point if outside model
!     if(icornerlong < 1) icornerlong = 1
!     if(icornerlong > NX_TOPO_ANT-1) icornerlong = NX_TOPO_ANT-1
!     if(icornerlat < 1) icornerlat = 1
!     if(icornerlat > NY_TOPO_ANT-1) icornerlat = NY_TOPO_ANT-1

! ! compute coordinates of corner
!     long_corner = ORIG_LONG_TOPO_ANT + (icornerlong-1)*DEGREES_PER_CELL_TOPO_ANT
!     lat_corner = ORIG_LAT_TOPO_ANT + (icornerlat-1)*DEGREES_PER_CELL_TOPO_ANT

! ! compute ratio for interpolation
!     ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO_ANT
!     ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO_ANT

! ! avoid edge effects
!     if(ratio_xi < 0.) ratio_xi = 0.
!     if(ratio_xi > 1.) ratio_xi = 1.
!     if(ratio_eta < 0.) ratio_eta = 0.
!     if(ratio_eta > 1.) ratio_eta = 1.

! ! interpolate elevation at current point
!     elevation_bedrock = &
!       ibedrock(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!       ibedrock(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!       ibedrock(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!       ibedrock(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

! !! DK DK exclude circles around each station to make sure they are on the bedrock
! !! DK DK and not in the ice
!   is_around_a_station = .false.
!   do istation = 1,NUMBER_OF_STATIONS
!     if(sqrt((xstore(i,j,k,ispec) - utm_x_station(istation))**2 + (ystore(i,j,k,ispec) - &
!          utm_y_station(istation))**2) < RADIUS_TO_EXCLUDE) then
!       is_around_a_station = .true.
!       exit
!     endif
!   enddo

! ! we are above the bedrock interface i.e. in the ice, and not too close to a station
!   if(zstore(i,j,k,ispec) >= elevation_bedrock .and. .not. is_around_a_station) then
!      iflag = flag_above
!      !iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_ICE
!      ! we are below the bedrock interface i.e. in the bedrock, or close to a station
!   else
!      iflag = flag_below
!      !iflag_attenuation_store(i,j,k,ispec) = IATTENUATION_BEDROCK
!   endif
    

! end subroutine interface
