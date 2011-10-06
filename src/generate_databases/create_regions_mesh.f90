!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 0
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

! for poroelastic model
  real(kind=CUSTOM_REAL),dimension(:,:,:,:), allocatable :: etastore,phistore,tortstore
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: rhoarraystore,kappaarraystore,permstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vpI,rho_vpII,rho_vsI

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass,rmass_acoustic,&
                            rmass_solid_poroelastic,rmass_fluid_poroelastic

! ocean load
  integer :: NGLOB_OCEAN
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qmu_attenuation_store

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

! acoustic-poroelastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_po_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_po_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_ac_po_ijk
  integer, dimension(:), allocatable :: coupling_ac_po_ispec
  integer :: num_coupling_ac_po_faces

! elastic-poroelastic coupling surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_el_po_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_el_po_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: coupling_el_po_ijk
  integer, dimension(:), allocatable :: coupling_el_po_ispec
  integer :: num_coupling_el_po_faces

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
                        ANISOTROPY,NPROC,OCEANS,TOPOGRAPHY, &
                        ATTENUATION,USE_OLSEN_ATTENUATION, &
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

! material properties
  integer :: nmat_ext_mesh,nundefMat_ext_mesh
  double precision, dimension(16,nmat_ext_mesh) :: materials_ext_mesh
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
  logical :: OCEANS,TOPOGRAPHY
  logical :: ATTENUATION,USE_OLSEN_ATTENUATION

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  integer, dimension(NX_TOPO,NY_TOPO) :: itopo_bathy

! local parameters
! static memory size needed by the solver
  double precision :: static_memory_size
  real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period

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
  call get_MPI(myrank,nglob,nspec,ibool, &
                        nelmnts_ext_mesh,elmnts_ext_mesh, &
                        my_nelmnts_neighbours_ext_mesh, my_interfaces_ext_mesh, &
                        ibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh,&
                        my_neighbours_ext_mesh,NPROC)

! sets material velocities
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...determining velocity model'
  endif
  call get_model(myrank,nspec,ibool,mat_ext_mesh,nelmnts_ext_mesh, &
                        materials_ext_mesh,nmat_ext_mesh, &
                        undef_mat_prop,nundefMat_ext_mesh, &
                        ANISOTROPY,LOCAL_PATH)

! sets up absorbing/free surface boundaries
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries '
  endif
  call get_absorbing_boundary(myrank,nspec,nglob,ibool, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh, &
                            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                            nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                            nodes_ibelm_bottom,nodes_ibelm_top, &
                            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                            nspec2D_bottom,nspec2D_top)

! sets up acoustic-elastic-poroelastic coupling surfaces
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...detecting acoustic-elastic-poroelastic surfaces '
  endif
  call get_coupling_surfaces(myrank, &
                        nspec,nglob,ibool,NPROC, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                        my_neighbours_ext_mesh)

! creates mass matrix
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
  endif
  call create_mass_matrices(nglob,nspec,ibool)

! creates ocean load mass matrix
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating ocean load mass matrix '
  endif
  call create_mass_matrices_ocean_load(nglob,nspec,ibool,OCEANS,TOPOGRAPHY, &
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
                        jacobianstore, rho_vp,rho_vs,qmu_attenuation_store, &
                        rhostore,kappastore,mustore, &
                        rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore, &
                        rho_vpI,rho_vpII,rho_vsI, &
                        rmass,rmass_acoustic,rmass_solid_poroelastic,rmass_fluid_poroelastic, &
                        OCEANS,rmass_ocean_load,NGLOB_OCEAN,ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                        abs_boundary_normal,abs_boundary_jacobian2Dw, &
                        abs_boundary_ijk,abs_boundary_ispec,num_abs_boundary_faces, &
                        free_surface_normal,free_surface_jacobian2Dw, &
                        free_surface_ijk,free_surface_ispec,num_free_surface_faces, &
                        coupling_ac_el_normal,coupling_ac_el_jacobian2Dw, &
                        coupling_ac_el_ijk,coupling_ac_el_ispec,num_coupling_ac_el_faces, &
                        coupling_ac_po_normal,coupling_ac_po_jacobian2Dw, &
                        coupling_ac_po_ijk,coupling_ac_po_ispec,num_coupling_ac_po_faces, &
                        coupling_el_po_normal,coupling_el_po_jacobian2Dw, &
                        coupling_el_po_ijk,coupling_el_po_ispec,num_coupling_el_po_faces, &
                        num_interfaces_ext_mesh,my_neighbours_ext_mesh,nibool_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                        prname,SAVE_MESH_FILES,ANISOTROPY,NSPEC_ANISO, &
                        c11store,c12store,c13store,c14store,c15store,c16store, &
                        c22store,c23store,c24store,c25store,c26store,c33store, &
                        c34store,c35store,c36store,c44store,c45store,c46store, &
                        c55store,c56store,c66store, &
                        ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

! computes the approximate amount of static memory needed to run the solver
  call memory_eval(nspec,nglob,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh, &
                  OCEANS,static_memory_size)
  call max_all_dp(static_memory_size, max_static_memory_size)

! checks the mesh, stability and resolved period
  call sync_all()
!chris: check for poro: At the moment cpI & cpII are for eta=0
  call check_mesh_resolution(myrank,nspec,nglob,ibool,&
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            kappastore,mustore,rho_vp,rho_vs, &
                            -1.0d0, model_speed_max,min_resolved_period, &
                            phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI )

! saves binary mesh files for attenuation
  if( ATTENUATION ) then
    call get_attenuation_model(myrank,nspec,USE_OLSEN_ATTENUATION, &
                          mustore,rho_vs,qmu_attenuation_store, &
                          ispec_is_elastic,min_resolved_period,prname)
  endif


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

! cleanup
  if( .not. SAVE_MOHO_MESH ) deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  deallocate(xixstore,xiystore,xizstore,&
              etaxstore,etaystore,etazstore,&
              gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore,qmu_attenuation_store)
  deallocate(kappastore,mustore,rho_vp,rho_vs)
  deallocate(rho_vpI,rho_vpII,rho_vsI)
  deallocate(rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore)

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
  if( ier /= 0 ) stop 'error allocating array xelm etc.'

  allocate( qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ),stat=ier)
  if( ier /= 0 ) stop 'error allocating array xigll etc.'

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ),stat=ier)
  if( ier /= 0 ) stop 'error allocating array wxgll etc.'

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), &
          dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
  if( ier /= 0 ) stop 'error allocating array shape3D etc.'

! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
          shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
          shape2D_bottom(NGNOD2D,NGLLX,NGLLY), &
          shape2D_top(NGNOD2D,NGLLX,NGLLY),stat=ier)
  if( ier /= 0 ) stop 'error allocating array shape2D_x etc.'

  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
          dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
          dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY), &
          dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if( ier /= 0 ) stop 'error allocating array dershape2D_x etc.'

  allocate(wgllwgll_xy(NGLLX,NGLLY), &
          wgllwgll_xz(NGLLX,NGLLZ), &
          wgllwgll_yz(NGLLY,NGLLZ),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
          kappastore(NGLLX,NGLLY,NGLLZ,nspec), &
          mustore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
          !vpstore(NGLLX,NGLLY,NGLLZ,nspec), &
          !vsstore(NGLLX,NGLLY,NGLLZ,nspec),
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with poroelastic model
  allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,nspec), &
          kappaarraystore(3,NGLLX,NGLLY,NGLLZ,nspec), &
          etastore(NGLLX,NGLLY,NGLLZ,nspec), &
          tortstore(NGLLX,NGLLY,NGLLZ,nspec), &
          phistore(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vpI(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vpII(NGLLX,NGLLY,NGLLZ,nspec), &
          rho_vsI(NGLLX,NGLLY,NGLLZ,nspec), &
          permstore(6,NGLLX,NGLLY,NGLLZ,nspec), stat=ier)
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

  ! free surface faces
  num_free_surface_faces = nspec2D_top

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

! get the 2-D shape functions
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
  integer :: ispec2D,ispec,icorner,iface,i,j,k,igll,iglob,ier
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
  allocate(iglob_is_surface(nglob), &
          iglob_normals(NDIM,nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating array iglob_is_surface'

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

  allocate(ibelm_moho_bot(NSPEC2D_MOHO), &
          ibelm_moho_top(NSPEC2D_MOHO), &
          normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
          normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO), &
          ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO), &
          ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if( ier /= 0 ) stop 'error allocating ibelm_moho_bot'

  ibelm_moho_bot = 0
  ibelm_moho_top = 0

  ! element flags
  allocate(is_moho_top(nspec), &
          is_moho_bot(nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating is_moho_top'
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
