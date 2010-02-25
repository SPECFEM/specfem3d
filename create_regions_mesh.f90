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


module create_regions_mesh_ext_par

  include 'constants.h'
  
! global point coordinates
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ystore_dummy
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: zstore_dummy

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  double precision, dimension(:), allocatable :: xelm,yelm,zelm

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! for model density, kappa, mu
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappastore,mustore  
  
! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass,rmass_acoustic,&
                            rmass_solid_poroelastic,rmass_fluid_poroelastic

! ocean load
  integer :: NGLOB_OCEAN
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! attenuation 
  integer, dimension(:,:,:,:), allocatable :: iflag_attenuation_store

! 2D shape functions and their derivatives, weights
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top
  double precision, dimension(:,:), allocatable :: wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ijk
  integer, dimension(:), allocatable :: abs_boundary_ispec
  integer :: num_abs_boundary_faces

! free surface arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ijk
  integer, dimension(:), allocatable :: free_surface_ispec
  integer :: num_free_surface_faces

! acoustic-elastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ijk
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer :: num_coupling_ac_el_faces

! for stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

! anisotropy
  integer :: NSPEC_ANISO
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
            c11store,c12store,c13store,c14store,c15store,c16store,&
            c22store,c23store,c24store,c25store,c26store,c33store,&
            c34store,c35store,c36store,c44store,c45store,c46store,&
            c55store,c56store,c66store

! material domain flags
  logical, dimension(:), allocatable :: ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic

! name of the database file
  character(len=256) prname
  
end module create_regions_mesh_ext_par

!
!-------------------------------------------------------------------------------------------------
!

! main routine

subroutine create_regions_mesh_ext(ibool, &
                        xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
                        nnodes_ext_mesh,nelmnts_ext_mesh, &
                        nodes_coords_ext_mesh, elmnts_ext_mesh, &
                        max_static_memory_size, mat_ext_mesh, materials_ext_mesh, &
                        nmat_ext_mesh, undef_mat_prop, nundefMat_ext_mesh, &
                        num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh, my_nelmnts_neighbours_ext_mesh, &
                        my_interfaces_ext_mesh, &
                        ibool_interfaces_ext_mesh, nibool_interfaces_ext_mesh, &
                        nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, &
                        NSPEC2D_BOTTOM, NSPEC2D_TOP,&
                        ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
                        nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax,&
                        nodes_ibelm_bottom,nodes_ibelm_top, &
                        SAVE_MESH_FILES,nglob, &
                        ANISOTROPY,NPROC,OCEANS, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! create the different regions of the mesh

  use create_regions_mesh_ext_par  
  implicit none
  !include "constants.h"

! number of spectral elements in each block
  integer :: nspec

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: npointot

! proc numbers for MPI
  integer :: myrank
  integer :: NPROC

  character(len=256) :: LOCAL_PATH

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! static memory size needed by the solver
  double precision :: max_static_memory_size 
  
  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh

!pll
  integer :: nmat_ext_mesh,nundefMat_ext_mesh 
  double precision, dimension(6,nmat_ext_mesh) :: materials_ext_mesh  
  character (len=30), dimension(6,nundefMat_ext_mesh):: undef_mat_prop
  
!  double precision, external :: materials_ext_mesh

! MPI communication
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: my_interfaces_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

! absorbing boundaries
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP
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

  integer :: nglob

  logical :: SAVE_MESH_FILES
  logical :: ANISOTROPY
  logical :: OCEANS

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy
  
! local parameters
! static memory size needed by the solver
  double precision :: static_memory_size
  real(kind=CUSTOM_REAL) :: model_speed_max
  
! for vtk output
!  character(len=256) prname_file
!  integer,dimension(:),allocatable :: itest_flag
!  integer, dimension(:), allocatable :: elem_flag

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
!real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ibedrock

! initializes arrays
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) 
    write(IMAIN,*) '  ...allocating arrays '
  endif
  call crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                        nspec2D_bottom,nspec2D_top,ANISOTROPY)
  

! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
! returns jacobianstore,xixstore,...gammazstore
! and GLL-point locations in xstore,ystore,zstore
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up jacobian '
  endif  
  call crm_ext_setup_jacobian(myrank, &                      
                        xstore,ystore,zstore,nspec, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,&
                        elmnts_ext_mesh,nelmnts_ext_mesh)
  
! sets material velocities
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...determining velocity model'
  endif
  call crm_ext_determine_velocity(nspec,&
                        mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY)
  
! creates ibool index array for projection from local to global points
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...indexing global points'
  endif
  call crm_ext_setup_indexing(ibool, &
                        xstore,ystore,zstore,nspec,nglob,npointot, &
                        nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! sets up MPI interfaces between partitions
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...preparing MPI interfaces '
  endif
  call crm_ext_prepare_MPI(myrank,nglob,nspec,ibool, &
                        nelmnts_ext_mesh,elmnts_ext_mesh, &
                        my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                        ibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC)

! creates mass matrix 
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
  endif
  call crm_ext_create_mass_matrix(nglob,nspec,ibool)
  
! sets up absorbing/free surface boundaries  
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries '
  endif
  call crm_ext_setup_abs_boundary(myrank,nspec,nglob,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)
    
! sets up acoustic-elastic coupling surfaces
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...detecting acoustic-elastic surfaces '
  endif
  call crm_ext_detect_ac_el_surface(myrank, &
                        nspec,nglob,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)

! creates ocean load mass matrix 
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating ocean load mass matrix '
  endif
  call crm_ext_create_ocean_load_mass(nglob,nspec,ibool,OCEANS,&
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)


! saves the binary mesh files
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...saving databases'
  endif
  !call create_name_database(prname,myrank,LOCAL_PATH)
  call save_arrays_solver_ext_mesh(nspec,nglob, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore, &
                        jacobianstore, rho_vp,rho_vs,iflag_attenuation_store, &
                        rhostore,kappastore,mustore, &
                        rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                        OCEANS,rmass_ocean_load,NGLOB_OCEAN,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec,num_abs_boundary_faces, &
                        free_surface_normal,free_surface_jacobian2Dw, &
                        free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
                        coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                        coupling_ac_el_ijk,coupling_ac_el_ispec,num_coupling_ac_el_faces, &                        
                        num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname,SAVE_MESH_FILES,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

! computes the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)

! checks the mesh, stability and resolved period 
  call sync_all()
  call check_mesh_resolution(myrank,nspec,nglob,ibool,&
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            kappastore,mustore,rho_vp,rho_vs, &
                            -1.0d0, model_speed_max )

! VTK file output
!  if( SAVE_MESH_FILES ) then
!    ! saves material flag assigned for each spectral element into a vtk file 
!    prname_file = prname(1:len_trim(prname))//'material_flag'
!    allocate(elem_flag(nspec))
!    elem_flag(:) = mat_ext_mesh(1,:)
!    call write_VTK_data_elem_i(nspec,nglob, &
!            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!            elem_flag,prname_file)
!    deallocate(elem_flag)
!    
!    !plotting abs boundaries
!    !  allocate(itest_flag(nspec))
!    !  itest_flag(:) = 0
!    !  do ispec=1,nspec
!    !    if( iboun(1,ispec) ) itest_flag(ispec) = 1
!    !  enddo
!    !  prname_file = prname(1:len_trim(prname))//'iboundary1_flag'
!    !  call write_VTK_data_elem_i(nspec,nglob, &
!    !            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!    !            itest_flag,prname_file)
!    !  deallocate(itest_flag)
!  endif  

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
  if( .not. SAVE_MOHO_MESH ) deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  deallocate(xixstore,xiystore,xizstore,&
              etaxstore,etaystore,etazstore,&
              gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore,iflag_attenuation_store)  
  deallocate(kappastore,mustore,rho_vp,rho_vs)

end subroutine create_regions_mesh_ext

!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                        nspec2D_bottom,nspec2D_top,ANISOTROPY)

  use create_regions_mesh_ext_par  
  implicit none

  integer :: nspec,myrank
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            nspec2D_bottom,nspec2D_top

  character(len=256) :: LOCAL_PATH
            
  logical :: ANISOTROPY

! local parameters  
  integer :: ier

! memory test
!  logical,dimension(:),allocatable :: test_mem 
!  
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
          mustore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier) !pll
          !vpstore(NGLLX,NGLLY,NGLLZ,nspec), &
          !vsstore(NGLLX,NGLLY,NGLLZ,nspec),          
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

! absorbing boundary 
  ! absorbing faces
  num_abs_boundary_faces = nspec2D_xmin + nspec2D_xmax + nspec2D_ymin + nspec2D_ymax + nspec2D_bottom
  ! adds faces of free surface if it also absorbs
  if( ABSORB_FREE_SURFACE ) num_abs_boundary_faces = num_abs_boundary_faces + nspec2D_top

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( abs_boundary_ispec(num_abs_boundary_faces), &
           abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! free surface
  ! free surface faces
  if( ABSORB_FREE_SURFACE ) then
    ! no free surface - uses a dummy size
    num_free_surface_faces = 1
  else
    num_free_surface_faces = nspec2D_top  
  endif

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( free_surface_ispec(num_free_surface_faces), &
           free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
           free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
           free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with anisotropy
  if( ANISOTROPY ) then
    NSPEC_ANISO = nspec
  else
    NSPEC_ANISO = 1
  endif
  allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), &
          c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! material flags
  allocate( ispec_is_acoustic(nspec), &
           ispec_is_elastic(nspec), &
           ispec_is_poroelastic(nspec), stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  
end subroutine crm_ext_allocate_arrays


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_jacobian(myrank, &                      
                        xstore,ystore,zstore,nspec, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,&
                        elmnts_ext_mesh,nelmnts_ext_mesh)

  use create_regions_mesh_ext_par
  implicit none

! number of spectral elements in each block
  integer :: nspec

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh

! proc numbers for MPI
  integer :: myrank

! local parameters
  integer :: ispec,ia,i,j,k

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

    ! CUBIT should provide a mesh ordering such that the 3D jacobian is defined
    ! (otherwise mesh would be degenerated)
    call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
                      etaxstore,etaystore,etazstore, &
                      gammaxstore,gammaystore,gammazstore,jacobianstore, &
                      xstore,ystore,zstore, &
                      xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)

  enddo

end subroutine crm_ext_setup_jacobian


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_determine_velocity(nspec,&
                        mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY)

  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: nspec

! external mesh
  integer :: nelmnts_ext_mesh
  integer :: nmat_ext_mesh,nundefMat_ext_mesh 

  integer, dimension(2,nelmnts_ext_mesh) :: mat_ext_mesh
  double precision, dimension(6,nmat_ext_mesh) :: materials_ext_mesh  
  character (len=30), dimension(6,nundefMat_ext_mesh):: undef_mat_prop

! anisotropy
  logical :: ANISOTROPY

! local parameters
  real(kind=CUSTOM_REAL) :: vp,vs,rho  
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  
  integer :: ispec,i,j,k,iundef,iflag_atten
  integer :: iflag,flag_below,flag_above
  integer :: iflag_aniso,idomain_id
  
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

  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

! material properties on all GLL points: taken from material values defined for 
! each spectral element in input mesh
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
! check if the material is known or unknown
           if (mat_ext_mesh(1,ispec) > 0) then
              ! density    
              ! materials_ext_mesh format: #index1 = rho #index2 = vp #index3 = vs #index4 = Q_flag #index5 = 0 
              !rhostore(i,j,k,ispec) = materials_ext_mesh(1,mat_ext_mesh(1,ispec))
              rho = materials_ext_mesh(1,mat_ext_mesh(1,ispec))
              
              ! isotropic values: vp, vs              
              !vpstore(i,j,k,ispec) = materials_ext_mesh(2,mat_ext_mesh(1,ispec))
              !vsstore(i,j,k,ispec) = materials_ext_mesh(3,mat_ext_mesh(1,ispec))
              vp = materials_ext_mesh(2,mat_ext_mesh(1,ispec))
              vs = materials_ext_mesh(3,mat_ext_mesh(1,ispec))

              ! attenuation
              iflag_atten = materials_ext_mesh(4,mat_ext_mesh(1,ispec))                            
              !change for piero :
              !if(mat_ext_mesh(1,ispec) == 1) then
              !   iflag_attenuation_store(i,j,k,ispec) = 1
              !else
              !   iflag_attenuation_store(i,j,k,ispec) = 2
              !endif
              
              ! anisotropy
              iflag_aniso = materials_ext_mesh(5,mat_ext_mesh(1,ispec))
              
              ! material domain_id
              idomain_id = materials_ext_mesh(6,mat_ext_mesh(1,ispec))
              
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
              !rhostore(i,j,k,ispec) = materials_ext_mesh(1,iflag)
              !vpstore(i,j,k,ispec) = materials_ext_mesh(2,iflag)
              !vsstore(i,j,k,ispec) = materials_ext_mesh(3,iflag)
              rho = materials_ext_mesh(1,iflag)
              vp = materials_ext_mesh(2,iflag)
              vs = materials_ext_mesh(3,iflag)
              iflag_atten = materials_ext_mesh(4,iflag)
              !change for piero :
              !  if(iflag == 1) then
              !     iflag_attenuation_store(i,j,k,ispec) = 1
              !  else
              !     iflag_attenuation_store(i,j,k,ispec) = 2
              !  endif
              iflag_aniso = materials_ext_mesh(5,iflag)
              idomain_id = materials_ext_mesh(6,iflag)
             else
              stop 'material: tomography not implemented yet'
             ! call tomography()
           end if

! adds anisotropic perturbation to vp, vs
           if( ANISOTROPY ) then
             call aniso_model(iflag_aniso,rho,vp,vs,c11,c12,c13,c14,c15,c16, &
                     c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66) 

             c11store(i,j,k,ispec) = c11
             c12store(i,j,k,ispec) = c12
             c13store(i,j,k,ispec) = c13
             c14store(i,j,k,ispec) = c14
             c15store(i,j,k,ispec) = c15
             c16store(i,j,k,ispec) = c16
             c22store(i,j,k,ispec) = c22
             c23store(i,j,k,ispec) = c23
             c24store(i,j,k,ispec) = c24
             c25store(i,j,k,ispec) = c25
             c26store(i,j,k,ispec) = c26
             c33store(i,j,k,ispec) = c33
             c34store(i,j,k,ispec) = c34
             c35store(i,j,k,ispec) = c35
             c36store(i,j,k,ispec) = c36
             c44store(i,j,k,ispec) = c44
             c45store(i,j,k,ispec) = c45
             c46store(i,j,k,ispec) = c46
             c55store(i,j,k,ispec) = c55
             c56store(i,j,k,ispec) = c56
             c66store(i,j,k,ispec) = c66
                     
           endif
! density
           rhostore(i,j,k,ispec) = rho
          
! kappa, mu
           kappastore(i,j,k,ispec) = rho*( vp*vp - FOUR_THIRDS*vs*vs )                
           mustore(i,j,k,ispec) = rho*vs*vs

! attenuation
           iflag_attenuation_store(i,j,k,ispec) = iflag_atten
           ! Stacey, a completer par la suite  
           rho_vp(i,j,k,ispec) = rho*vp
           rho_vs(i,j,k,ispec) = rho*vs
           !end pll

! material domain
           !print*,'velocity model:',ispec,idomain_id           
           if( idomain_id == IDOMAIN_ACOUSTIC ) then
             ispec_is_acoustic(ispec) = .true.            
           else if( idomain_id == IDOMAIN_ELASTIC ) then
             ispec_is_elastic(ispec) = .true.
           else if( idomain_id == IDOMAIN_POROELASTIC ) then
             stop 'poroelastic material domain not implemented yet'
             ispec_is_poroelastic(ispec) = .true.
           else
             stop 'error material domain index'
           endif
        enddo
      enddo
    enddo
    !print*,myrank,'ispec:',ispec,'rho:',rhostore(1,1,1,ispec),'vp:',vpstore(1,1,1,ispec),'vs:',vsstore(1,1,1,ispec)    
  enddo

! checks material domains
  do ispec=1,nspec
    if( (ispec_is_acoustic(ispec) .eqv. .false.) &
          .and. (ispec_is_elastic(ispec) .eqv. .false.) &
          .and. (ispec_is_poroelastic(ispec) .eqv. .false.) ) then
      print*,'error material domain not assigned to element:',ispec
      print*,'acoustic: ',ispec_is_acoustic(ispec)
      print*,'elastic: ',ispec_is_elastic(ispec)
      print*,'poroelastic: ',ispec_is_poroelastic(ispec)      
      stop 'error material domain index element'
    endif
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

end subroutine crm_ext_determine_velocity


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_indexing(ibool, &
                            xstore,ystore,zstore,nspec,nglob,npointot, &
                            nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! creates global indexing array ibool

  use create_regions_mesh_ext_par 
  implicit none

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
  integer :: i,j,k,ispec,iglobnum

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

end subroutine crm_ext_setup_indexing

!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_prepare_MPI(myrank,nglob,nspec,ibool, &
                                    nelmnts_ext_mesh,elmnts_ext_mesh, &
                                    my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                                    ibool_interfaces_ext_mesh, &
                                    nibool_interfaces_ext_mesh, &
                                    num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                                    my_neighbours_ext_mesh,NPROC)

! sets up the MPI interface for communication between partitions

  use create_regions_mesh_ext_par 
  implicit none

  integer :: myrank,nglob,nspec,NPROC

! global indexing
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! external mesh, element indexing  
  integer :: nelmnts_ext_mesh
  integer, dimension(ESIZE,nelmnts_ext_mesh) :: elmnts_ext_mesh
  
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  
  integer, dimension(num_interfaces_ext_mesh) :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: my_interfaces_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh  
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh


  !integer :: nnodes_ext_mesh
  !double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh  
  
!local parameters
  double precision, dimension(:), allocatable :: xp,yp,zp
  double precision, dimension(:), allocatable :: work_ext_mesh

  integer, dimension(:), allocatable :: locval
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh_true

  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface_ext_mesh,ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh
  integer, dimension(:), allocatable :: ibool_interface_ext_mesh_dummy
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2 

  ! assembly test
  integer :: i,j,k,ispec,iglob,count,inum
  integer :: max_nibool_interfaces_ext_mesh
  integer,dimension(:),allocatable :: test_flag
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr
  integer, dimension(:,:), allocatable :: ibool_interfaces_dummy  

! gets global indices for points on MPI interfaces (defined by my_interfaces_ext_mesh) between different partitions
! and stores them in ibool_interfaces_ext_mesh & nibool_interfaces_ext_mesh (number of total points)
  call prepare_assemble_MPI( nelmnts_ext_mesh,elmnts_ext_mesh, &
                            ibool,nglob,ESIZE, &
                            num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
                            my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                            ibool_interfaces_ext_mesh, &
                            nibool_interfaces_ext_mesh )

  allocate(nibool_interfaces_ext_mesh_true(num_interfaces_ext_mesh))

! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  do iinterface = 1, num_interfaces_ext_mesh

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

    ! gets x,y,z coordinates of global points on MPI interface
    do ilocnum = 1, nibool_interfaces_ext_mesh(iinterface)
      xp(ilocnum) = xstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      yp(ilocnum) = ystore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
      zp(ilocnum) = zstore_dummy(ibool_interfaces_ext_mesh(ilocnum,iinterface))
    enddo

    ! sorts (lexicographically?) ibool_interfaces_ext_mesh and updates value
    ! of total number of points nibool_interfaces_ext_mesh_true(iinterface)
    call sort_array_coordinates(nibool_interfaces_ext_mesh(iinterface),xp,yp,zp, &
         ibool_interfaces_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
         reorder_interface_ext_mesh,locval,ifseg,nibool_interfaces_ext_mesh_true(iinterface), &
         ind_ext_mesh,ninseg_ext_mesh,iwork_ext_mesh,work_ext_mesh)

    ! checks that number of MPI points are still the same
    num_points1 = num_points1 + nibool_interfaces_ext_mesh(iinterface)
    num_points2 = num_points2 + nibool_interfaces_ext_mesh_true(iinterface)    
    if( num_points1 /= num_points2 ) then
      write(*,*) 'error sorting MPI interface points:',myrank
      write(*,*) '   interface:',iinterface,num_points1,num_points2
      call exit_mpi(myrank,'error sorting MPI interface')
    endif
    !write(*,*) myrank,'intfc',iinterface,num_points2,nibool_interfaces_ext_mesh_true(iinterface)
    
    ! cleanup temporary arrays
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

  ! cleanup
  deallocate(nibool_interfaces_ext_mesh_true)

  ! outputs total number of MPI interface points
  call sum_all_i(num_points2,ilocnum)  
  if( myrank == 0 ) then
    write(IMAIN,*) '     total MPI interface points: ',ilocnum  
  endif
  
! checks with assembly of test fields
  allocate(test_flag(nglob),test_flag_cr(nglob))
  test_flag(:) = 0
  test_flag_cr(:) = 0._CUSTOM_REAL
  count = 0
  do ispec = 1, nspec    
    ! sets flags on global points
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,k,ispec)         
          
          ! counts number of unique global points to set
          if( test_flag(iglob) == 0 ) count = count+1
          
          ! sets identifier
          test_flag(iglob) = myrank + 1 
          test_flag_cr(iglob) = myrank + 1.0
        enddo
      enddo
    enddo
  enddo
  call sync_all()

  ! collects contributions from different MPI partitions
  ! sets up MPI communications
  max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
  allocate(ibool_interfaces_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh))
  
  count = 0
  do iinterface = 1, num_interfaces_ext_mesh
     ibool_interfaces_dummy(:,iinterface) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,iinterface)
     count = count + nibool_interfaces_ext_mesh(iinterface)
     !write(*,*) myrank,'interfaces ',iinterface,nibool_interfaces_ext_mesh(iinterface),max_nibool_interfaces_ext_mesh
  enddo
  call sync_all()
  
  call sum_all_i(count,iglob)
  if( myrank == 0 ) then
    if( iglob /= ilocnum ) call exit_mpi(myrank,'error total global MPI interface points')
  endif
  
  ! adds contributions from different partitions to flag arrays
  ! integer arrays
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,test_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_dummy,&
                        my_neighbours_ext_mesh)
  ! custom_real arrays
  call assemble_MPI_scalar_ext_mesh(NPROC,nglob,test_flag_cr, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_dummy, &
                        my_neighbours_ext_mesh)

  ! checks number of interface points
  i = 0
  j = 0
  do iglob=1,nglob
    ! only counts flags with MPI contributions
    if( test_flag(iglob) > myrank+1 ) i = i + 1
    if( test_flag_cr(iglob) > myrank+1.0) j = j + 1
  enddo  
  call sum_all_i(i,inum)
  call sum_all_i(j,iglob)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total assembled MPI interface points:',inum
    if( inum /= iglob .or. inum > ilocnum ) call exit_mpi(myrank,'error MPI assembly')
  endif
  
end subroutine crm_ext_prepare_MPI


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_create_mass_matrix(nglob,nspec,ibool)

! returns precomputed mass matrix in rmass array
  
  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: nspec
  integer :: nglob
  
! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: jacobianl
  integer :: ispec,i,j,k,iglob,ier

! allocates memory
  allocate(rmass(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(rmass_acoustic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(rmass_solid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'
  allocate(rmass_fluid_poroelastic(nglob),stat=ier); if(ier /= 0) stop 'error in allocate'

! creates mass matrix  
  rmass(:) = 0._CUSTOM_REAL
  rmass_acoustic(:) = 0._CUSTOM_REAL
  rmass_solid_poroelastic(:) = 0._CUSTOM_REAL
  rmass_fluid_poroelastic(:) = 0._CUSTOM_REAL
  
  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)

          weight = wxgll(i)*wygll(j)*wzgll(k)
          jacobianl = jacobianstore(i,j,k,ispec)

! acoustic mass matrix
          if( ispec_is_acoustic(ispec) ) then
            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                    sngl( dble(jacobianl) * weight / dble(kappastore(i,j,k,ispec)) )
            else
               rmass_acoustic(iglob) = rmass_acoustic(iglob) + &
                    jacobianl * weight / kappastore(i,j,k,ispec)
            endif
          endif

! elastic mass matrix
          if( ispec_is_elastic(ispec) ) then
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass(iglob) = rmass(iglob) + &
                    sngl( dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) )
            else
               rmass(iglob) = rmass(iglob) + &
                    jacobianl * weight * rhostore(i,j,k,ispec)
            endif
          endif
          
! poroelastic mass matrices
          if( ispec_is_poroelastic(ispec) ) then
            
            stop 'poroelastic mass matrices not implemented yet'
            
            !rho_solid = density(1,kmato(ispec))
            !rho_fluid = density(2,kmato(ispec))
            !phi = porosity(kmato(ispec))
            !tort = tortuosity(kmato(ispec))
            !rho_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f          
            !
            !if(CUSTOM_REAL == SIZE_REAL) then            
            !  ! for the solid mass matrix
            !  rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
            !      sngl( dble(jacobianl) * weight * dble(rho_bar - phi*rho_fluid/tort) )
            !  
            !  ! for the fluid mass matrix
            !  rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
            !      sngl( dble(jacobianl) * weight * dble(rho_bar*rho_fluid*tort - &
            !                                  phi*rho_fluid*rho_fluid)/dble(rho_bar*phi) )            
            !else
            !  rmass_solid_poroelastic(iglob) = rmass_solid_poroelastic(iglob) + &
            !      jacobianl * weight * (rho_bar - phi*rho_fluid/tort)
            !  
            !  rmass_fluid_poroelastic(iglob) = rmass_fluid_poroelastic(iglob) + &
            !      jacobianl * weight * (rho_bar*rho_fluid*tort - &
            !                                  phi*rho_fluid*rho_fluid) / (rho_bar*phi) 
            !endif
          endif
          
        enddo
      enddo
    enddo
  enddo ! nspec  

end subroutine crm_ext_create_mass_matrix


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_abs_boundary(myrank,nspec,nglob,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)

! determines absorbing boundaries/free-surface, 2D jacobians, face normals for Stacey conditions

  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec,nglob

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! data from the external mesh
  integer :: nnodes_ext_mesh 
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

! absorbing boundaries (as defined in CUBIT)
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
  
! local parameters
  logical, dimension(:,:),allocatable :: iboun   ! pll 

  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  integer:: ijk_face(3,NGLLX,NGLLY)
  
  ! corner locations for faces
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xcoord_iboun,ycoord_iboun,zcoord_iboun
  
  ! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  integer  :: ispec,ispec2D,icorner,ier,iabs,iface,igll,i,j,igllfree,ifree
  
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)                              

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

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)                              

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

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ)                              

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

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLZ
      do i=1,NGLLY
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLY,NGLLZ) 

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

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY) 

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

    ! sets face infos
    iabs = iabs + 1
    abs_boundary_ispec(iabs) = ispec      
    
    ! gll points -- assuming NGLLX = NGLLY = NGLLZ
    igll = 0
    do j=1,NGLLY
      do i=1,NGLLX
        igll = igll+1
        abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
        abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
        abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
      enddo
    enddo    
    
  enddo
  
  ! top 
  ! free surface face counter
  ifree = 0
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
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLY) 

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

    ! stores surface infos
    if( .not. ABSORB_FREE_SURFACE ) then
      ! store for free surface
      !jacobian2D_top(:,:,ispec2D) = jacobian2Dw_face(:,:)
      !normal_top(:,:,:,ispec2D) = normal_face(:,:,:)  

      ! sets face infos
      ifree = ifree + 1
      free_surface_ispec(ifree) = ispec      
      
      ! gll points -- assuming NGLLX = NGLLY = NGLLZ
      igllfree = 0
      do j=1,NGLLY
        do i=1,NGLLX
          igllfree = igllfree+1
          free_surface_ijk(:,igllfree,ifree) = ijk_face(:,i,j)
          free_surface_jacobian2Dw(igllfree,ifree) = jacobian2Dw_face(i,j)
          free_surface_normal(:,igllfree,ifree) = normal_face(:,i,j)  
        enddo
      enddo        
    else
      ! adds face infos to absorbing boundary surface
      iabs = iabs + 1
      abs_boundary_ispec(iabs) = ispec      
      
      ! gll points -- assuming NGLLX = NGLLY = NGLLZ
      igll = 0
      do j=1,NGLLY
        do i=1,NGLLX
          igll = igll+1
          abs_boundary_ijk(:,igll,iabs) = ijk_face(:,i,j)
          abs_boundary_jacobian2Dw(igll,iabs) = jacobian2Dw_face(i,j)
          abs_boundary_normal(:,igll,iabs) = normal_face(:,i,j)  
        enddo
      enddo
      
      ! resets free surface 
      ifree = 1
      free_surface_ispec(:) = 0
      free_surface_ijk(:,:,:) = 0
      free_surface_jacobian2Dw(:,:) = 0.0
      free_surface_normal(:,:,:) = 0.0
    endif
  enddo
  
! checks counters  
  if( ifree /= num_free_surface_faces ) then  
    print*,'error number of free surface faces:',ifree,num_free_surface_faces
    stop 'error number of free surface faces'
  endif
  
  if( iabs /= num_abs_boundary_faces ) then
    print*,'error number of absorbing faces:',iabs,num_abs_boundary_faces
    stop 'error number of absorbing faces'
  endif

  call sum_all_i(num_abs_boundary_faces,iabs)
  if( myrank == 0 ) then
    write(IMAIN,*) '     absorbing boundary:'
    write(IMAIN,*) '     total number of faces = ',iabs
    if( ABSORB_FREE_SURFACE ) then
    write(IMAIN,*) '     absorbing boundary includes free surface'
    endif
  endif

end subroutine crm_ext_setup_abs_boundary


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_detect_ac_el_surface(myrank, &
                        nspec,nglob,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)
                            
! determines coupling surface for acoustic-elastic domains

  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: myrank,nspec,nglob,NPROC

! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! MPI communication
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: my_neighbours_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
            ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

! local parameters
  ! (assumes NGLLX=NGLLY=NGLLZ)
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: tmp_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: tmp_jacobian2Dw
  integer :: ijk_face(3,NGLLX,NGLLY)
  integer,dimension(:,:,:),allocatable :: tmp_ijk
  integer,dimension(:),allocatable :: tmp_ispec

  integer,dimension(NGNOD2D) :: iglob_corners_ref !,iglob_corners
  integer :: ispec,i,j,k,igll,ier,iglob
  integer :: inum,iface_ref,icorner,iglob_midpoint ! iface,ispec_neighbor
  integer :: count_elastic,count_acoustic
  
  ! mpi interface communication
  integer, dimension(:), allocatable :: elastic_flag,acoustic_flag,test_flag
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy
  integer :: max_nibool_interfaces_ext_mesh
  logical, dimension(:), allocatable :: mask_ibool
  
  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top  
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)               
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))   ! top  

  
  ! test vtk output
  !integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: gll_data
  !character(len=256):: prname_file
  
! allocates temporary arrays  
  allocate(tmp_normal(NDIM,NGLLSQUARE,nspec*6))
  allocate(tmp_jacobian2Dw(NGLLSQUARE,nspec*6))  
  allocate(tmp_ijk(3,NGLLSQUARE,nspec*6))
  allocate(tmp_ispec(nspec*6))
  tmp_ispec(:) = 0
  tmp_ijk(:,:,:) = 0
  tmp_normal(:,:,:) = 0.0
  tmp_jacobian2Dw(:,:) = 0.0
  
  ! sets flags for acoustic / elastic on global points
  allocate(elastic_flag(nglob),stat=ier)
  allocate(acoustic_flag(nglob),stat=ier)  
  allocate(test_flag(nglob),stat=ier)  
  allocate(mask_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocate flag array'  
  elastic_flag(:) = 0
  acoustic_flag(:) = 0
  test_flag(:) = 0
  count_elastic = 0
  count_acoustic = 0
  do ispec = 1, nspec
    ! counts elements
    if( ispec_is_elastic(ispec) ) count_elastic = count_elastic + 1
    if( ispec_is_acoustic(ispec) ) count_acoustic = count_acoustic + 1
    
    ! sets flags on global points
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,k,ispec)         
          ! sets elastic flag
          if( ispec_is_elastic(ispec) ) elastic_flag(iglob) =  myrank+1
          ! sets acoustic flag
          if( ispec_is_acoustic(ispec) ) acoustic_flag(iglob) =  myrank+1
          ! sets test flag
          test_flag(iglob) = myrank+1
        enddo
      enddo
    enddo
  enddo
  call sum_all_i(count_acoustic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total acoustic elements:',inum
  endif   
  call sum_all_i(count_elastic,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     total elastic elements :',inum
  endif   

  ! collects contributions from different MPI partitions
  ! sets up MPI communications
  max_nibool_interfaces_ext_mesh = maxval( nibool_interfaces_ext_mesh(:) )
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'  
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,i) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,i)
  enddo  
  ! sums elastic flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,elastic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)
  ! sums acoustic flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,acoustic_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)

  ! sums test flags
  call assemble_MPI_scalar_i_ext_mesh(NPROC,nglob,test_flag, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh_dummy,&
                        my_neighbours_ext_mesh)

  ! loops over all element faces and 
  ! counts number of coupling faces between acoustic and elastic elements
  mask_ibool(:) = .false.
  inum = 0    
  do ispec=1,nspec

    ! loops over each face
    do iface_ref= 1, 6      

      ! takes indices of corners of reference face
      do icorner = 1,NGNOD2D
        i = iface_all_corner_ijk(1,icorner,iface_ref)
        j = iface_all_corner_ijk(2,icorner,iface_ref)
        k = iface_all_corner_ijk(3,icorner,iface_ref)
        ! global reference indices
        iglob_corners_ref(icorner) = ibool(i,j,k,ispec)

        ! reference corner coordinates
        xcoord(icorner) = xstore_dummy(iglob_corners_ref(icorner))
        ycoord(icorner) = ystore_dummy(iglob_corners_ref(icorner))
        zcoord(icorner) = zstore_dummy(iglob_corners_ref(icorner))                  
      enddo
      
      ! checks if face has acoustic side
      if( acoustic_flag( iglob_corners_ref(1) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(2) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(3) ) >= 1 .and. &
         acoustic_flag( iglob_corners_ref(4) ) >= 1) then        
        ! checks if face is has an elastic side 
        if( elastic_flag( iglob_corners_ref(1) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(2) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(3) ) >= 1 .and. &
           elastic_flag( iglob_corners_ref(4) ) >= 1) then

          ! reference midpoint on face (used to avoid redundant face counting)
          i = iface_all_midpointijk(1,iface_ref)
          j = iface_all_midpointijk(2,iface_ref)
          k = iface_all_midpointijk(3,iface_ref)      
          iglob_midpoint = ibool(i,j,k,ispec)

          ! checks if points on this face are masked already
          if( .not. mask_ibool(iglob_midpoint) ) then

            ! gets face GLL points i,j,k indices from element face
            call get_element_face_gll_indices(iface_ref,ijk_face,NGLLX,NGLLY)
            
            ! takes each element face only once, if it lies on an MPI interface
            ! note: this is not exactly load balanced
            !          lowest rank process collects as many faces as possible, second lowest as so forth
            if( (test_flag(iglob_midpoint) == myrank+1) .or. &
               (test_flag(iglob_midpoint) > 2*(myrank+1)) ) then
            
              ! gets face GLL 2Djacobian, weighted from element face
              call get_jacobian_boundary_face(myrank,nspec, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                        ispec,iface_ref,jacobian2Dw_face,normal_face,NGLLX,NGLLY)

              ! normal convention: points away from acoustic, reference element
              !                                switch normal direction if necessary
              do j=1,NGLLY
                do i=1,NGLLX
                    ! directs normals such that they point outwards of element
                    call get_element_face_normal(ispec,iface_ref,xcoord,ycoord,zcoord, &
                                                ibool,nspec,nglob, &
                                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                                normal_face(:,i,j) )
                    ! makes sure that it always points away from acoustic element, 
                    ! otherwise switch direction
                    if( ispec_is_elastic(ispec) ) normal_face(:,i,j) = - normal_face(:,i,j)
                enddo
              enddo

              ! stores informations about this face
              inum = inum + 1
              tmp_ispec(inum) = ispec
              igll = 0
              do j=1,NGLLY
                do i=1,NGLLX
                  ! adds all gll points on this face
                  igll = igll + 1
                  
                  ! do we need to store local i,j,k,ispec info? or only global indices iglob?
                  tmp_ijk(:,igll,inum) = ijk_face(:,i,j)
                  
                  ! stores weighted jacobian and normals
                  tmp_jacobian2Dw(igll,inum) = jacobian2Dw_face(i,j)
                  tmp_normal(:,igll,inum) = normal_face(:,i,j)
                  
                  ! masks global points ( to avoid redundant counting of faces)
                  iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec)
                  mask_ibool(iglob) = .true.
                enddo
              enddo
            else
              ! assumes to be already collected by lower rank process, masks face points
              do j=1,NGLLY
                do i=1,NGLLX
                  iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec)
                  mask_ibool(iglob) = .true. 
                enddo
              enddo
            endif ! test_flag
          endif ! mask_ibool          
        endif ! elastic_flag
      endif ! acoustic_flag
    enddo ! iface_ref
  enddo ! ispec
    
! stores completed coupling face informations  
! 
! note: no need to store material parameters on these coupling points 
!          for acoustic-elastic interface
  num_coupling_ac_el_faces = inum
  allocate(coupling_ac_el_normal(NDIM,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_jacobian2Dw(NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ijk(3,NGLLSQUARE,num_coupling_ac_el_faces))
  allocate(coupling_ac_el_ispec(num_coupling_ac_el_faces))
  do inum = 1,num_coupling_ac_el_faces
    coupling_ac_el_normal(:,:,inum) = tmp_normal(:,:,inum)
    coupling_ac_el_jacobian2Dw(:,inum) = tmp_jacobian2Dw(:,inum)
    coupling_ac_el_ijk(:,:,inum) = tmp_ijk(:,:,inum)
    coupling_ac_el_ispec(inum) = tmp_ispec(inum)    
  enddo

! user output
  call sum_all_i(num_coupling_ac_el_faces,inum)
  if( myrank == 0 ) then
    write(IMAIN,*) '     acoustic-elastic coupling:'
    write(IMAIN,*) '     total number of faces = ',inum
  endif  

end subroutine crm_ext_detect_ac_el_surface

!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_create_ocean_load_mass(nglob,nspec,ibool,OCEANS,&
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! returns precomputed mass matrix in rmass array
  
  use create_regions_mesh_ext_par 
  implicit none

! number of spectral elements in each block
  integer :: nspec
  integer :: nglob
  
! arrays with the mesh global indices
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  logical :: OCEANS

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

  
! local parameters
  double precision :: weight
  double precision :: xval,yval,long,lat,elevation
  double precision :: height_oceans
  double precision :: long_corner,lat_corner,ratio_xi,ratio_eta
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D,igll,iglobnum
  integer :: icornerlong,icornerlat

! creates ocean load mass matrix
  if(OCEANS) then

    ! adding ocean load mass matrix at ocean bottom
    NGLOB_OCEAN = nglob
    allocate(rmass_ocean_load(NGLOB_OCEAN))

    ! create ocean load mass matrix for degrees of freedom at ocean bottom
    rmass_ocean_load(:) = 0._CUSTOM_REAL

    ! add contribution of the oceans for surface elements exactly at ocean bottom
    do ispec2D = 1,num_free_surface_faces

      ispec_oceans = free_surface_ispec(ispec2D)

      ! only adds contribution if top boundary is elastic, no need to add this approximate calculation
      ! if top is already acoustic/poroelastic
      if( ispec_is_elastic(ispec_oceans) ) then

        do igll=1,NGLLSQUARE
          ix_oceans = free_surface_ijk(1,igll,ispec2D)
          iy_oceans = free_surface_ijk(1,igll,ispec2D)
          iz_oceans = free_surface_ijk(1,igll,ispec2D)
        
          iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

          ! compute local height of oceans

          ! get coordinates of current point
          xval = xstore_dummy(iglobnum)
          yval = ystore_dummy(iglobnum)

          ! project x and y in UTM back to long/lat since topo file is in long/lat
          call utm_geo(long,lat,xval,yval,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

          ! get coordinate of corner in bathy/topo model
          icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
          icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

          ! avoid edge effects and extend with identical point if outside model
          if(icornerlong < 1) icornerlong = 1
          if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
          if(icornerlat < 1) icornerlat = 1
          if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

          ! compute coordinates of corner
          long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
          lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

          ! compute ratio for interpolation
          ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
          ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

          ! avoid edge effects
          if(ratio_xi < 0.) ratio_xi = 0.
          if(ratio_xi > 1.) ratio_xi = 1.
          if(ratio_eta < 0.) ratio_eta = 0.
          if(ratio_eta > 1.) ratio_eta = 1.

          ! interpolate elevation at current point
          elevation = &
                itopo_bathy(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
                itopo_bathy(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
                itopo_bathy(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
                itopo_bathy(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

          ! suppress positive elevation, which means no oceans
          if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation)
          endif

          ! take into account inertia of water column
          weight = dble( free_surface_jacobian2Dw(igll,ispec2D)) &
                   * dble(RHO_OCEANS) * height_oceans

          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
          else
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
          endif

        enddo ! igll
      endif ! ispec_is_elastic
    enddo ! num_free_surface_faces

    ! add regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

    ! allocate dummy array if no oceans
    NGLOB_OCEAN = 1
    allocate(rmass_ocean_load(NGLOB_OCEAN))

  endif

end subroutine crm_ext_create_ocean_load_mass


!
!-------------------------------------------------------------------------------------------------
!

subroutine create_regions_mesh_save_moho( myrank,nglob,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )
  
  use create_regions_mesh_ext_par  
  implicit none

  integer :: nspec2D_moho_ext
  integer, dimension(nspec2D_moho_ext) :: ibelm_moho
  integer, dimension(4,nspec2D_moho_ext) :: nodes_ibelm_moho

  integer :: myrank,nglob,nspec

  ! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! local parameters        
  ! Moho mesh
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot  
  
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord    
  real(kind=CUSTOM_REAL) :: jacobian2Dw_face(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: normal_face(NDIM,NGLLX,NGLLY)
  real(kind=CUSTOM_REAL),dimension(NDIM):: normal
  integer :: ijk_face(3,NGLLX,NGLLY)  

  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: iglob_normals
  integer,dimension(:),allocatable:: iglob_is_surface

  integer :: imoho_bot,imoho_top
  integer :: ispec2D,ispec,icorner,iface,i,j,k,igll,iglob
  integer :: iglob_midpoint,idirect,counter
  integer :: imoho_top_all,imoho_bot_all,imoho_all
  
  ! corners indices of reference cube faces
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/))   ! xmin
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
             reshape( (/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/))   ! xmax
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
             reshape( (/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/))   ! ymin
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
             reshape( (/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/))   ! ymax
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
             reshape( (/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/))  ! bottom
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
             reshape( (/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/))   ! top  
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
             reshape( (/ iface1_corner_ijk,iface2_corner_ijk, &
                 iface3_corner_ijk,iface4_corner_ijk, &
                 iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/))   ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)               
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))   ! top  
  
  ! temporary arrays for passing information
  allocate(iglob_is_surface(nglob))
  allocate(iglob_normals(NDIM,nglob))
  iglob_is_surface = 0
  iglob_normals = 0._CUSTOM_REAL
  
  ! loops over given moho surface elements
  do ispec2D=1, nspec2D_moho_ext

    ! gets element id
    ispec = ibelm_moho(ispec2D)
           
    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    ! (note: uses point locations rather than point indices to find the element face,
    !            because the indices refer no more to the newly indexed ibool array )
    do icorner=1,NGNOD2D
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_moho(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_moho(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_moho(icorner,ispec2D))
    enddo
    
    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)    
    
    ! weighted jacobian and normal                          
    call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)                              

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

    ! stores information on global points on moho surface
    igll = 0
    do j=1,NGLLY    
      do i=1,NGLLX
        iglob = ibool(ijk_face(1,i,j),ijk_face(2,i,j),ijk_face(3,i,j),ispec) 
        ! sets flag
        iglob_is_surface(iglob) = ispec2D
        ! sets normals
        iglob_normals(:,iglob) = normal_face(:,i,j)
      enddo
    enddo
  enddo
  
  ! stores moho elements
  NSPEC2D_MOHO = nspec2D_moho_ext
  
  allocate(ibelm_moho_bot(NSPEC2D_MOHO))
  allocate(ibelm_moho_top(NSPEC2D_MOHO))
  allocate(normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO))
  allocate(ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO))
  ibelm_moho_bot = 0
  ibelm_moho_top = 0
  
  ! element flags 
  allocate(is_moho_top(nspec))
  allocate(is_moho_bot(nspec))
  is_moho_top = .false.
  is_moho_bot = .false.  

  ! finds spectral elements with moho surface
  imoho_top = 0
  imoho_bot = 0
  do ispec=1,nspec
  
    ! loops over each face
    do iface = 1,6      
      ! checks if corners of face on surface
      counter = 0
      do icorner = 1,NGNOD2D
        i = iface_all_corner_ijk(1,icorner,iface)
        j = iface_all_corner_ijk(2,icorner,iface)
        k = iface_all_corner_ijk(3,icorner,iface)
        iglob = ibool(i,j,k,ispec)
  
        ! checks if point on surface  
        if( iglob_is_surface(iglob) > 0 ) then
          counter = counter+1
  
          ! reference corner coordinates
          xcoord(icorner) = xstore_dummy(iglob)
          ycoord(icorner) = ystore_dummy(iglob)
          zcoord(icorner) = zstore_dummy(iglob)
        endif
      enddo

      ! stores moho informations    
      if( counter == NGNOD2D ) then

        ! gets face GLL points i,j,k indices from element face
        call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

        ! re-computes face infos  
        ! weighted jacobian and normal                          
        call get_jacobian_boundary_face(myrank,nspec, & 
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob,&
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,&                                          
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ)                              

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

        ! takes normal stored temporary on a face midpoint
        i = iface_all_midpointijk(1,iface)
        j = iface_all_midpointijk(2,iface)
        k = iface_all_midpointijk(3,iface)      
        iglob_midpoint = ibool(i,j,k,ispec)
        normal(:) = iglob_normals(:,iglob_midpoint)
        
        ! determines whether normal points into element or not (top/bottom distinction)
        call get_element_face_normal_idirect(ispec,iface,xcoord,ycoord,zcoord, &
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              normal,idirect )        

        ! takes moho surface element id given by id on midpoint
        ispec2D = iglob_is_surface(iglob_midpoint)

        ! sets face infos for bottom (normal points away from element)
        if( idirect == 1 ) then
          
          ! checks validity
          if( is_moho_bot( ispec) .eqv. .true. ) then
            print*,'error: moho surface geometry bottom'
            print*,'  does not allow for mulitple element faces in kernel computation'
            call exit_mpi(myrank,'error moho bottom elements')
          endif
          
          imoho_bot = imoho_bot + 1        
          is_moho_bot(ispec) = .true.
          ibelm_moho_bot(ispec2D) = ispec

          ! stores on surface gll points (assuming NGLLX = NGLLY = NGLLZ)
          igll = 0
          do j=1,NGLLZ
            do i=1,NGLLX
              igll = igll+1
              ijk_moho_bot(:,igll,ispec2D) = ijk_face(:,i,j)
              normal_moho_bot(:,igll,ispec2D) = normal_face(:,i,j)  
            enddo
          enddo 
                 
        ! sets face infos for top element  
        else if( idirect == 2 ) then

          ! checks validity
          if( is_moho_top( ispec) .eqv. .true. ) then
            print*,'error: moho surface geometry top'
            print*,'  does not allow for mulitple element faces kernel computation'
            call exit_mpi(myrank,'error moho top elements')
          endif

          imoho_top = imoho_top + 1        
          is_moho_top(ispec) = .true.
          ibelm_moho_top(ispec2D) = ispec

          ! gll points 
          igll = 0
          do j=1,NGLLZ
            do i=1,NGLLX
              igll = igll+1
              ijk_moho_top(:,igll,ispec) = ijk_face(:,i,j)
              ! note: top elements have normal pointing into element
              normal_moho_top(:,igll,ispec) = - normal_face(:,i,j)  
            enddo
          enddo                
        endif
    
      endif ! counter
      
    enddo ! iface
    
    ! checks validity of top/bottom distinction
    if( is_moho_top(ispec) .and. is_moho_bot(ispec) ) then
      print*,'error: moho surface elements confusing'
      print*,'  element:',ispec,'has top and bottom surface'
      call exit_mpi(myrank,'error moho surface element')
    endif
    
  enddo ! ispec2D
  
  ! note: surface e.g. could be at the free-surface and have no top elements etc...
  ! user output
  call sum_all_i( imoho_top, imoho_top_all )
  call sum_all_i( imoho_bot, imoho_bot_all )
  call sum_all_i( NSPEC2D_MOHO, imoho_all )
  if( myrank == 0 ) then
    write(IMAIN,*) '********'
    write(IMAIN,*) 'Moho surface:'
    write(IMAIN,*) '    total surface elements: ',imoho_all
    write(IMAIN,*) '    top elements   :',imoho_top_all
    write(IMAIN,*) '    bottom elements:',imoho_bot_all
    write(IMAIN,*) '********'
  endif

  ! saves moho files: total number of elements, corner points, all points
  open(unit=27,file=prname(1:len_trim(prname))//'ibelm_moho.bin',status='unknown',form='unformatted')
  write(27) NSPEC2D_MOHO
  write(27) ibelm_moho_top
  write(27) ibelm_moho_bot
  write(27) ijk_moho_top
  write(27) ijk_moho_bot
  close(27)
  open(unit=27,file=prname(1:len_trim(prname))//'normal_moho.bin',status='unknown',form='unformatted')
  write(27) normal_moho_top
  write(27) normal_moho_bot
  close(27)
  open(unit=27,file=prname(1:len_trim(prname))//'is_moho.bin',status='unknown',form='unformatted')
  write(27) is_moho_top
  write(27) is_moho_bot
  close(27)
  
end subroutine create_regions_mesh_save_moho

!
!-------------------------------------------------------------------------------------------------
!

subroutine bubble_sort( arr, ndim )

! sorts values in array arr[ndim] in increasing order

  implicit none
  
  integer :: ndim
  integer :: arr(ndim)
  
  logical :: swapped
  integer :: j,tmp
  
  swapped = .true.
  do while( swapped )
    swapped = .false.
    do j = 1, ndim-1
      if( arr(j+1) < arr(j) ) then
        tmp = arr(j) 
        arr(j) = arr(j+1)
        arr(j+1) = tmp
        swapped = .true.        
      endif
    enddo
  enddo
  
end subroutine bubble_sort

!
!-------------------------------------------------------------------------------------------------
!


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
