!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

  subroutine create_regions_mesh()

! create the different regions of the mesh

  use generate_databases_par, only: nspec => NSPEC_AB,nglob => NGLOB_AB, &
      ibool,xstore,ystore,zstore, &
      npointot,myrank,LOCAL_PATH, &
      nnodes_ext_mesh,nelmnts_ext_mesh, &
      nodes_coords_ext_mesh, elmnts_ext_mesh, &
      max_memory_size,num_interfaces_ext_mesh, max_interface_size_ext_mesh, &
      my_neighbors_ext_mesh, my_nelmnts_neighbors_ext_mesh, &
      my_interfaces_ext_mesh, &
      ibool_interfaces_ext_mesh, nibool_interfaces_ext_mesh, &
      STACEY_ABSORBING_CONDITIONS, nspec2D_xmin, nspec2D_xmax, &
      nspec2D_ymin, nspec2D_ymax, &
      NSPEC2D_BOTTOM, NSPEC2D_TOP, &
      ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
      nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
      nodes_ibelm_bottom,nodes_ibelm_top, &
      SAVE_MESH_FILES,PML_CONDITIONS, &
      ANISOTROPY,NPROC,APPROXIMATE_OCEAN_LOAD,OLSEN_ATTENUATION_RATIO, &
      ATTENUATION,USE_OLSEN_ATTENUATION, &
      nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
      ADIOS_FOR_MESH,IMAIN,SAVE_MOHO_MESH,ATTENUATION_f0_REFERENCE

  use create_regions_mesh_ext_par
  use fault_generate_databases, only: fault_read_input,fault_setup, &
                          fault_save_arrays,fault_save_arrays_test, &
                          nnodes_coords_open,nodes_coords_open,ANY_FAULT_IN_THIS_PROC, &
                          ANY_FAULT, PARALLEL_FAULT

  implicit none

! local parameters
! memory size needed by the solver
  double precision :: memory_size
  real(kind=CUSTOM_REAL) :: model_speed_max,min_resolved_period

! initializes arrays
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...allocating arrays '
    call flush_IMAIN()
  endif
  call crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                               nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                               nspec2D_bottom,nspec2D_top,ANISOTROPY, &
                               nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh,ANY_FAULT_IN_THIS_PROC)

 ! if faults exist this reads nodes_coords_open
  call fault_read_input(prname,myrank)

! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
! returns jacobianstore,xixstore,...gammazstore
! and GLL-point locations in xstore,ystore,zstore
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up jacobian '
    call flush_IMAIN()
  endif
  if (ANY_FAULT_IN_THIS_PROC) then
   ! compute jacobians with fault open and *store needed for ibool.
    call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                nodes_coords_open, nnodes_coords_open,elmnts_ext_mesh,nelmnts_ext_mesh)
  else ! with no fault
    call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh)
  endif


! creates ibool index array for projection from local to global points
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...indexing global points'
    call flush_IMAIN()
  endif
  if (ANY_FAULT_IN_THIS_PROC) then
    call crm_ext_setup_indexing(ibool, &
                                xstore,ystore,zstore,nspec,nglob,npointot, &
                                nnodes_coords_open,nodes_coords_open,myrank)
  else ! with no fault
    call crm_ext_setup_indexing(ibool, &
                                xstore,ystore,zstore,nspec,nglob,npointot, &
                                nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)
  endif

  if (ANY_FAULT) then
    ! recalculate *store with faults closed
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ... resetting up jacobian in fault domains'
      call flush_IMAIN()
    endif
    if (ANY_FAULT_IN_THIS_PROC) then
      call crm_ext_setup_jacobian(myrank,xstore,ystore,zstore,nspec, &
                                  nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                  elmnts_ext_mesh,nelmnts_ext_mesh)
    endif
    ! at this point (xyz)store_dummy are still open
    if (.not. PARALLEL_FAULT) then
      call fault_setup (ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                        xstore,ystore,zstore,nspec,nglob,myrank)
    endif
    ! this closes (xyz)store_dummy
  endif


! sets up MPI interfaces between partitions
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...preparing MPI interfaces '
    call flush_IMAIN()
  endif
  call get_MPI(myrank,nglob_dummy,nspec,ibool, &
               nelmnts_ext_mesh,elmnts_ext_mesh, &
               my_nelmnts_neighbors_ext_mesh, my_interfaces_ext_mesh, &
               ibool_interfaces_ext_mesh, &
               nibool_interfaces_ext_mesh, &
               num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
               my_neighbors_ext_mesh)

! setting up parallel fault
  if (PARALLEL_FAULT .and. ANY_FAULT) then
    call synchronize_all()
    !at this point (xyz)store_dummy are still open
    call fault_setup (ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                      xstore,ystore,zstore,nspec,nglob,myrank)
   ! this closes (xyz)store_dummy
  endif

! sets up absorbing/free surface boundaries
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up absorbing boundaries'
    call flush_IMAIN()
  endif
  call get_absorbing_boundary(myrank,nspec,ibool, &
                              nodes_coords_ext_mesh,nnodes_ext_mesh, &
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                              nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                              nodes_ibelm_bottom,nodes_ibelm_top, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              nspec2D_bottom,nspec2D_top)

! sets up mesh surface
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...setting up mesh surface'
    call flush_IMAIN()
  endif
  call crm_setup_mesh_surface()

! sets up up Moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...setting up Moho surface'
      call flush_IMAIN()
    endif
    call crm_setup_moho(myrank,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )
  endif

! sets material velocities
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...determining velocity model'
    call flush_IMAIN()
  endif
  call get_model(myrank)

! sets up acoustic-elastic-poroelastic coupling surfaces
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...detecting acoustic-elastic-poroelastic surfaces '
    call flush_IMAIN()
  endif
  call get_coupling_surfaces(myrank, &
                             nspec,ibool,NPROC, &
                             nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                             num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                             my_neighbors_ext_mesh)

! locates inner and outer elements
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...element inner/outer separation '
    call flush_IMAIN()
  endif
  call crm_setup_inner_outer_elemnts(myrank,nspec, &
                                     num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                                     nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                     ibool,SAVE_MESH_FILES)

! colors mesh if requested
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...element mesh coloring '
    call flush_IMAIN()
  endif
  call setup_color_perm(myrank,nspec,nglob,ibool,ANISOTROPY,SAVE_MESH_FILES)

! overwrites material parameters from external binary files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...external binary models '
    call flush_IMAIN()
  endif
  call get_model_binaries(myrank,nspec,LOCAL_PATH)

! calculates damping profiles and auxiliary coefficients on all C-PML points
  if (PML_CONDITIONS) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...creating C-PML damping profiles '
      call flush_IMAIN()
    endif
    call pml_set_local_dampingcoeff(myrank,xstore_dummy,ystore_dummy,zstore_dummy)
  endif

! creates mass matrix
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...creating mass matrix '
    call flush_IMAIN()
  endif
  call create_mass_matrices(nglob_dummy,nspec,ibool,PML_CONDITIONS,STACEY_ABSORBING_CONDITIONS)

! saves the binary mesh files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...saving databases'
    call flush_IMAIN()
  endif
  !call create_name_database(prname,myrank,LOCAL_PATH)
  if (ADIOS_FOR_MESH) then
    call save_arrays_solver_ext_mesh_adios(nspec,nglob,APPROXIMATE_OCEAN_LOAD, &
                                           ibool,num_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           max_interface_size_ext_mesh, &
                                           ibool_interfaces_ext_mesh,SAVE_MESH_FILES,ANISOTROPY)
  else
    call save_arrays_solver_ext_mesh(nspec,nglob_dummy,APPROXIMATE_OCEAN_LOAD,ibool, &
                                     num_interfaces_ext_mesh,my_neighbors_ext_mesh,nibool_interfaces_ext_mesh, &
                                     max_interface_size_ext_mesh,ibool_interfaces_ext_mesh, &
                                     SAVE_MESH_FILES,ANISOTROPY)
  endif

! saves faults
  if (ANY_FAULT) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving fault databases'
      call flush_IMAIN()
    endif
    !  call fault_save_arrays_test(prname)  ! for debugging
    call fault_save_arrays(prname)
  endif

! saves moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving Moho surfaces'
      call flush_IMAIN()
    endif
    call crm_save_moho()
  endif

! computes the approximate amount of memory needed to run the solver
  call synchronize_all()
  call memory_eval(nspec,nglob_dummy,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh, &
                   APPROXIMATE_OCEAN_LOAD,memory_size)
  call max_all_dp(memory_size, max_memory_size)

! checks the mesh, stability and resolved period
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  ...checking mesh resolution'
    call flush_IMAIN()
  endif

  if (POROELASTIC_SIMULATION) then
    !chris: check for poro: At the moment cpI & cpII are for eta=0
    call check_mesh_resolution_poro(myrank,nspec,nglob_dummy,ibool, &
                                    xstore_dummy,ystore_dummy,zstore_dummy, &
                                    -1.0d0, model_speed_max,min_resolved_period, &
                                    phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                                    LOCAL_PATH,SAVE_MESH_FILES )
  else
    call check_mesh_resolution(myrank,nspec,nglob_dummy, &
                               ibool,xstore_dummy,ystore_dummy,zstore_dummy, &
                               kappastore,mustore,rho_vp,rho_vs, &
                               -1.0d0,model_speed_max,min_resolved_period, &
                               LOCAL_PATH,SAVE_MESH_FILES)
  endif

! saves binary mesh files for attenuation
  if (ATTENUATION) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  ...saving attenuation databases'
      call flush_IMAIN()
    endif
    call get_attenuation_model(myrank,nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore,rho_vs,kappastore,rho_vp,qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)
  endif

  ! cleanup
  deallocate(xixstore,xiystore,xizstore, &
             etaxstore,etaystore,etazstore, &
             gammaxstore,gammaystore,gammazstore)
  deallocate(jacobianstore)
  deallocate(qkappa_attenuation_store,qmu_attenuation_store)
  deallocate(kappastore,mustore,rho_vp,rho_vs)
  deallocate(rho_vpI,rho_vpII,rho_vsI)
  deallocate(rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore)

  if (.not. SAVE_MOHO_MESH) then
    deallocate(xstore_dummy,ystore_dummy,zstore_dummy)
  endif

  if (ACOUSTIC_SIMULATION) then
    deallocate(rmass_acoustic)
  endif

  if (ELASTIC_SIMULATION) then
    deallocate(rmass)
  endif

  if (POROELASTIC_SIMULATION) then
    deallocate(rmass_solid_poroelastic,rmass_fluid_poroelastic)
  endif

  if (STACEY_ABSORBING_CONDITIONS) then
     if (ELASTIC_SIMULATION) then
       deallocate(rmassx,rmassy,rmassz)
     endif
     if (ACOUSTIC_SIMULATION) then
       deallocate(rmassz_acoustic)
     endif
  endif

end subroutine create_regions_mesh

!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_allocate_arrays(nspec,LOCAL_PATH,myrank, &
                        nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                        nspec2D_bottom,nspec2D_top,ANISOTROPY, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh, &
                        elmnts_ext_mesh,nelmnts_ext_mesh,ANY_FAULT_IN_THIS_PROC)

  use generate_databases_par, only: STACEY_INSTEAD_OF_FREE_SURFACE,PML_INSTEAD_OF_FREE_SURFACE, &
    NGNOD,NGNOD2D,NDIM,NDIM2D,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,BOTTOM_FREE_SURFACE,nspec_irregular

  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec,myrank
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            nspec2D_bottom,nspec2D_top

  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  logical :: ANISOTROPY,ANY_FAULT_IN_THIS_PROC

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(NGNOD,nelmnts_ext_mesh) :: elmnts_ext_mesh


! local parameters
  integer :: ier,ispec,ia
  logical :: any_regular_elem
  real    :: cube_edge_size_squared
  real, dimension(NGNOD) :: xelm_real,yelm_real,zelm_real

  allocate(xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),stat=ier)
  if (ier /= 0) stop 'error allocating array xelm etc.'

! attenuation
  allocate(qkappa_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  allocate(qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ),stat=ier)
  if (ier /= 0) stop 'error allocating array xigll etc.'

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ),stat=ier)
  if (ier /= 0) stop 'error allocating array wxgll etc.'

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), &
           dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) stop 'error allocating array shape3D etc.'

! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
           shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
           shape2D_bottom(NGNOD2D,NGLLX,NGLLY), &
           shape2D_top(NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'error allocating array shape2D_x etc.'

  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
           dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
           dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY), &
           dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'error allocating array dershape2D_x etc.'

  allocate(wgllwgll_xy(NGLLX,NGLLY), &
           wgllwgll_xz(NGLLX,NGLLZ), &
           wgllwgll_yz(NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec), &
           rho_vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
           kappastore(NGLLX,NGLLY,NGLLZ,nspec), &
           mustore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  if (POROELASTIC_SIMULATION) then
    NSPEC_PORO = nspec
  else
    NSPEC_PORO = 1
  endif

! array with poroelastic model
  allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           etastore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           phistore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_PORO), &
           permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_PORO), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! get the number of regular and irregular elements
  any_regular_elem = .false.
  allocate(irregular_element_number(nspec))
  irregular_element_number(:) = 0
  nspec_irregular = nspec

  do ispec = 1, nspec
    do ia = 1,NGNOD
      xelm_real(ia) = real(nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec)))
      yelm_real(ia) = real(nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec)))
      zelm_real(ia) = real(nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec)))
    enddo

    ! checks if element is regular (is a cube)
    call check_element_regularity(xelm_real,yelm_real,zelm_real,any_regular_elem,cube_edge_size_squared, &
                                  nspec_irregular,ispec,nspec,irregular_element_number,ANY_FAULT_IN_THIS_PROC)
  enddo

  if (nspec_irregular == 0) then
    allocate(xixstore(1,1,1,1), &
             xiystore(1,1,1,1), &
             xizstore(1,1,1,1), &
             etaxstore(1,1,1,1), &
             etaystore(1,1,1,1), &
             etazstore(1,1,1,1), &
             gammaxstore(1,1,1,1), &
             gammaystore(1,1,1,1), &
             gammazstore(1,1,1,1), &
             jacobianstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  else
    allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             xiystore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             xizstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             etaxstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             etaystore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             etazstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             gammaxstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             gammaystore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             gammazstore(NGLLX,NGLLY,NGLLZ,nspec_irregular), &
             jacobianstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
   if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  endif

! absorbing boundary
  ! absorbing faces
  num_abs_boundary_faces = nspec2D_xmin + nspec2D_xmax + nspec2D_ymin + nspec2D_ymax + nspec2D_bottom
  ! adds faces of free surface if it also absorbs
  if (STACEY_INSTEAD_OF_FREE_SURFACE .or. PML_INSTEAD_OF_FREE_SURFACE) then
     num_abs_boundary_faces = num_abs_boundary_faces + nspec2D_top
  endif

  ! subtract bottom surface from Stacey condition when free surface
  if (BOTTOM_FREE_SURFACE)  num_abs_boundary_faces = num_abs_boundary_faces -  nspec2D_bottom

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate(abs_boundary_ispec(num_abs_boundary_faces), &
           abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces), &
           abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  ! free surface faces
  num_free_surface_faces = nspec2D_top

   ! add bottom surface to free surface condition
  if (BOTTOM_FREE_SURFACE) then
     num_free_surface_faces = num_free_surface_faces + nspec2D_bottom
     if (STACEY_INSTEAD_OF_FREE_SURFACE)  num_free_surface_faces = num_free_surface_faces - nspec2D_top
  endif
  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  allocate( free_surface_ispec(num_free_surface_faces), &
            free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces), &
            free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces), &
            free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! array with anisotropy
  if (ANISOTROPY) then
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
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  ! material flags
  allocate(ispec_is_acoustic(nspec), &
           ispec_is_elastic(nspec), &
           ispec_is_poroelastic(nspec), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  ! initializes Moho surface
  NSPEC2D_MOHO = 0

end subroutine crm_ext_allocate_arrays


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_jacobian(myrank, &
                        xstore,ystore,zstore,nspec, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh, &
                        elmnts_ext_mesh,nelmnts_ext_mesh)

  use generate_databases_par, only: NGNOD,NGNOD2D,NDIM,NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA,CUSTOM_REAL
  use create_regions_mesh_ext_par

  implicit none

! number of spectral elements in each block
  integer :: nspec

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(NGNOD,nelmnts_ext_mesh) :: elmnts_ext_mesh

! proc numbers for MPI
  integer :: myrank

! local parameters
  integer :: ispec,ia,i,j,k,ispec_irreg
  logical :: any_regular_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix_reg,xiy_reg,xiz_reg,etax_reg,etay_reg,etaz_reg, &
                                                          gammax_reg,gammay_reg,gammaz_reg,jacobian_reg

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)

! get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)

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

  any_regular_element = .false.

  do ispec = 1, nspec
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo

    ! CUBIT should provide a mesh ordering such that the 3D jacobian is defined
    ! (otherwise mesh would be degenerated)
    ispec_irreg = irregular_element_number(ispec)
    !irregular_element_number is 0 only if the element is regular
    if (ispec_irreg /= 0 ) then
      call calc_jacobian(myrank,xixstore(:,:,:,ispec_irreg),xiystore(:,:,:,ispec_irreg),xizstore(:,:,:,ispec_irreg), &
                         etaxstore(:,:,:,ispec_irreg),etaystore(:,:,:,ispec_irreg),etazstore(:,:,:,ispec_irreg), &
                         gammaxstore(:,:,:,ispec_irreg),gammaystore(:,:,:,ispec_irreg),gammazstore(:,:,:,ispec_irreg), &
                         jacobianstore(:,:,:,ispec_irreg),xelm,yelm,zelm,dershape3D)
    else
      any_regular_elem = .true.
    endif
    call calc_coords(xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec), &
                     xelm,yelm,zelm,shape3D)

  enddo
  ! get xix derivative and jacobian on a regular element
  if (any_regular_elem) then

    !find a regular element
    ispec = 1
    do while (irregular_element_number(ispec) /= 0)
      ispec = ispec + 1
    enddo
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo
    call calc_jacobian(myrank,xix_reg,xiy_reg,xiz_reg, &
                       etax_reg,etay_reg,etaz_reg, &
                       gammax_reg,gammay_reg,gammaz_reg, &
                       jacobian_reg,xelm,yelm,zelm,dershape3D)
    xix_regular = xix_reg(1,1,1)
    jacobian_regular  = jacobian_reg(1,1,1)

  endif

end subroutine crm_ext_setup_jacobian


!
!-------------------------------------------------------------------------------------------------
!

subroutine crm_ext_setup_indexing(ibool, &
                            xstore,ystore,zstore,nspec,nglob,npointot, &
                            nnodes_ext_mesh,nodes_coords_ext_mesh,myrank)

! creates global indexing array ibool

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,NDIM
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
  double precision :: x_min,x_max

! allocate memory for arrays
  allocate(locval(npointot), &
           ifseg(npointot), &
           xp(npointot), &
           yp(npointot), &
           zp(npointot),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

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

  ! min/max values in x-direction
  x_min = minval(nodes_coords_ext_mesh(1,:))
  x_max = maxval(nodes_coords_ext_mesh(1,:))

! gets ibool indexing from local (GLL points) to global points
  call get_global(npointot,xp,yp,zp,ibool,locval,ifseg,nglob,x_min,x_max)

! we can create a new indirect addressing to reduce cache misses
  call get_global_indirect_addressing(nspec,nglob,ibool)

! cleanup
  deallocate(xp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if (ier /= 0) stop 'error in deallocate'

! unique global point locations
  nglob_dummy = nglob
  allocate(xstore_dummy(nglob_dummy), &
           ystore_dummy(nglob_dummy), &
           zstore_dummy(nglob_dummy),stat=ier)
  if (ier /= 0) stop 'error in allocate'
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

  subroutine crm_setup_moho( myrank,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )

  use generate_databases_par, only: NGNOD2D,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,SIZE_REAL,IMAIN, &
    NDIM,NGLLSQUARE,NGNOD2D_FOUR_CORNERS
  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec2D_moho_ext
  integer, dimension(nspec2D_moho_ext) :: ibelm_moho
  integer, dimension(NGNOD2D,nspec2D_moho_ext) :: nodes_ibelm_moho

  integer :: myrank,nspec

  ! data from the external mesh
  integer :: nnodes_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! local parameters
  real(kind=CUSTOM_REAL),dimension(NGNOD2D_FOUR_CORNERS) :: xcoord,ycoord,zcoord
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
  allocate(iglob_is_surface(nglob_dummy), &
          iglob_normals(NDIM,nglob_dummy),stat=ier)
  if (ier /= 0) stop 'error allocating array iglob_is_surface'

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
    do icorner=1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_moho(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_moho(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_moho(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                            ibool,nspec,nglob_dummy, &
                            xstore_dummy,ystore_dummy,zstore_dummy, &
                            iface)

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j=1,NGLLY
      do i=1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
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
  if (ier /= 0) stop 'error allocating ibelm_moho_bot'

  ibelm_moho_bot = 0
  ibelm_moho_top = 0

  ! element flags
  allocate(is_moho_top(nspec), &
          is_moho_bot(nspec),stat=ier)
  if (ier /= 0) stop 'error allocating is_moho_top'
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
      do icorner = 1,NGNOD2D_FOUR_CORNERS
        i = iface_all_corner_ijk(1,icorner,iface)
        j = iface_all_corner_ijk(2,icorner,iface)
        k = iface_all_corner_ijk(3,icorner,iface)
        iglob = ibool(i,j,k,ispec)

        ! checks if point on surface
        if (iglob_is_surface(iglob) > 0) then
          counter = counter+1

          ! reference corner coordinates
          xcoord(icorner) = xstore_dummy(iglob)
          ycoord(icorner) = ystore_dummy(iglob)
          zcoord(icorner) = zstore_dummy(iglob)
        endif
      enddo

      ! stores moho informations
      if (counter == NGNOD2D_FOUR_CORNERS) then

        ! gets face GLL points i,j,k indices from element face
        call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

        ! re-computes face infos
        ! weighted jacobian and normal
        call get_jacobian_boundary_face(myrank,nspec, &
              xstore_dummy,ystore_dummy,zstore_dummy,ibool,nglob_dummy, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
              ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

        ! normal convention: points away from element
        ! switch normal direction if necessary
        do j=1,NGLLZ
          do i=1,NGLLX
            call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                      ibool,nspec,nglob_dummy, &
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
                              ibool,nspec,nglob_dummy, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              normal,idirect )

        ! takes moho surface element id given by id on midpoint
        ispec2D = iglob_is_surface(iglob_midpoint)

        ! sets face infos for bottom (normal points away from element)
        if (idirect == 1) then

          ! checks validity
          if (is_moho_bot( ispec) .eqv. .true.) then
            print *,'error: moho surface geometry bottom'
            print *,'  does not allow for mulitple element faces in kernel computation'
            call exit_mpi(myrank,'error moho bottom elements')
          endif

          imoho_bot = imoho_bot + 1
          is_moho_bot(ispec) = .true.
          ibelm_moho_bot(ispec2D) = ispec

          ! stores on surface GLL points (assuming NGLLX = NGLLY = NGLLZ)
          igll = 0
          do j=1,NGLLZ
            do i=1,NGLLX
              igll = igll+1
              ijk_moho_bot(:,igll,ispec2D) = ijk_face(:,i,j)
              normal_moho_bot(:,igll,ispec2D) = normal_face(:,i,j)
            enddo
          enddo

        ! sets face infos for top element
        else if (idirect == 2) then

          ! checks validity
          if (is_moho_top( ispec) .eqv. .true.) then
            print *,'error: moho surface geometry top'
            print *,'  does not allow for mulitple element faces kernel computation'
            call exit_mpi(myrank,'error moho top elements')
          endif

          imoho_top = imoho_top + 1
          is_moho_top(ispec) = .true.
          ibelm_moho_top(ispec2D) = ispec

          ! GLL points
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
    if (is_moho_top(ispec) .and. is_moho_bot(ispec)) then
      print *,'error: moho surface elements confusing'
      print *,'  element:',ispec,'has top and bottom surface'
      call exit_mpi(myrank,'error moho surface element')
    endif

  enddo ! ispec2D

  ! note: surface e.g. could be at the free-surface and have no top elements etc...
  ! user output
  call sum_all_i( imoho_top, imoho_top_all )
  call sum_all_i( imoho_bot, imoho_bot_all )
  call sum_all_i( NSPEC2D_MOHO, imoho_all )
  if (myrank == 0) then
    write(IMAIN,*) '     Moho surface:'
    write(IMAIN,*) '     total surface elements: ',imoho_all
    write(IMAIN,*) '     top elements   :',imoho_top_all
    write(IMAIN,*) '     bottom elements:',imoho_bot_all
    call flush_IMAIN()
  endif

  deallocate(iglob_is_surface)
  deallocate(iglob_normals)

  end subroutine crm_setup_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_save_moho()

  use generate_databases_par, only: ADIOS_FOR_MESH,IOUT
  use create_regions_mesh_ext_par

  implicit none

  ! local parameters
  integer :: ier

  if (ADIOS_FOR_MESH) then
    call crm_save_moho_adios()
  else
    ! saves moho files: total number of elements, corner points, all points
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ibelm_moho.bin', &
          status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening ibelm_moho.bin file'
    write(IOUT) NSPEC2D_MOHO
    write(IOUT) ibelm_moho_top
    write(IOUT) ibelm_moho_bot
    write(IOUT) ijk_moho_top
    write(IOUT) ijk_moho_bot
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'normal_moho.bin', &
          status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening normal_moho.bin file'
    write(IOUT) normal_moho_top
    write(IOUT) normal_moho_bot
    close(IOUT)
    open(unit=IOUT,file=prname(1:len_trim(prname))//'is_moho.bin', &
      status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'error opening is_moho.bin file'
    write(IOUT) is_moho_top
    write(IOUT) is_moho_bot
    close(IOUT)
  endif

  end subroutine crm_save_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_inner_outer_elemnts(myrank,nspec, &
                                  num_interfaces_ext_mesh,max_interface_size_ext_mesh, &
                                  nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                  ibool,SAVE_MESH_FILES)

! locates inner and outer elements

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,IMAIN
  use create_regions_mesh_ext_par

  implicit none

  integer :: myrank,nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! MPI interfaces
  integer :: num_interfaces_ext_mesh,max_interface_size_ext_mesh
  integer, dimension(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh) :: &
    ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh

  logical :: SAVE_MESH_FILES

  ! local parameters
  integer :: i,j,k,ispec,iglob
  integer :: iinterface,ier
  integer :: ispec_inner,ispec_outer
  real :: percentage_edge
  character(len=MAX_STRING_LEN) :: filename
  logical,dimension(:),allocatable :: iglob_is_inner

  logical,parameter :: DEBUG = .false.

  ! allocates arrays
  allocate(ispec_is_inner(nspec),stat=ier)
  if (ier /= 0) stop 'error allocating array ispec_is_inner'

  ! temporary array
  allocate(iglob_is_inner(nglob_dummy),stat=ier)
  if (ier /= 0) stop 'error allocating temporary array  iglob_is_inner'

  ! initialize flags
  ispec_is_inner(:) = .true.
  iglob_is_inner(:) = .true.
  do iinterface = 1, num_interfaces_ext_mesh
    do i = 1, nibool_interfaces_ext_mesh(iinterface)
      iglob = ibool_interfaces_ext_mesh(i,iinterface)
      iglob_is_inner(iglob) = .false.
    enddo
  enddo

  ! determines flags for inner elements (purely inside the partition)
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          ispec_is_inner(ispec) = ( iglob_is_inner(iglob) .and. ispec_is_inner(ispec) )
        enddo
      enddo
    enddo
  enddo

  ! frees temporary array
  deallocate( iglob_is_inner )

  if (SAVE_MESH_FILES .and. DEBUG) then
    filename = prname(1:len_trim(prname))//'ispec_is_inner'
    call write_VTK_data_elem_l(nspec,nglob_dummy, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        ispec_is_inner,filename)
  endif


  ! sets up elements for loops in acoustic simulations
  nspec_inner_acoustic = 0
  nspec_outer_acoustic = 0
  if (ACOUSTIC_SIMULATION) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_acoustic = max(nspec_inner_acoustic,nspec_outer_acoustic)
    if (num_phase_ispec_acoustic < 0) stop 'error acoustic simulation: num_phase_ispec_acoustic is < zero'

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_acoustic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_acoustic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_acoustic = 0
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0) stop 'error allocating dummy array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0
  endif

  ! sets up elements for loops in acoustic simulations
  nspec_inner_elastic = 0
  nspec_outer_elastic = 0
  if (ELASTIC_SIMULATION) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_elastic = nspec_inner_elastic + 1
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_elastic = max(nspec_inner_elastic,nspec_outer_elastic)
    if (num_phase_ispec_elastic < 0) stop 'error elastic simulation: num_phase_ispec_elastic is < zero'

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_elastic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_elastic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_elastic = 0
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0) stop 'error allocating dummy array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0
  endif

  ! sets up elements for loops in poroelastic simulations
  nspec_inner_poroelastic = 0
  nspec_outer_poroelastic = 0
  if (POROELASTIC_SIMULATION) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_poroelastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_poroelastic = nspec_inner_poroelastic + 1
        else
          nspec_outer_poroelastic = nspec_outer_poroelastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_poroelastic = max(nspec_inner_poroelastic,nspec_outer_poroelastic)
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0) stop 'error allocating array phase_ispec_inner_poroelastic'
    nspec_inner_poroelastic = 0
    nspec_outer_poroelastic = 0
    do ispec = 1, nspec
      if (ispec_is_poroelastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_poroelastic = nspec_inner_poroelastic + 1
          phase_ispec_inner_poroelastic(nspec_inner_poroelastic,2) = ispec
        else
          nspec_outer_poroelastic = nspec_outer_poroelastic + 1
          phase_ispec_inner_poroelastic(nspec_outer_poroelastic,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_poroelastic = 0
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0) stop 'error allocating dummy array phase_ispec_inner_poroelastic'
    phase_ispec_inner_poroelastic(:,:) = 0
  endif

  ! user output
  if (myrank == 0) then
    percentage_edge = 100.*count(ispec_is_inner(:))/real(nspec)
    write(IMAIN,*) '     for overlapping of communications with calculations:'
    write(IMAIN,*) '     percentage of   edge elements ',100. -percentage_edge,'%'
    write(IMAIN,*) '     percentage of volume elements ',percentage_edge,'%'
    call flush_IMAIN()
  endif

  end subroutine crm_setup_inner_outer_elemnts

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_mesh_surface()

! gets number of surface elements
! (for receiver detection, movie outputs)

  use generate_databases_par, only: NPROC,NSPEC_AB,NGLOB_AB,ibool, &
    ispec_is_surface_external_mesh,iglob_is_surface_external_mesh,nfaces_surface,nfaces_surface_glob_ext_mesh, &
    max_nibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
    num_interfaces_ext_mesh,my_neighbors_ext_mesh, &
    MOVIE_TYPE,NSPEC2D_TOP

  use create_regions_mesh_ext_par

  implicit none

  ! local parameters
  integer :: i,ier
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy

  ! allocates mesh surface arrays
  allocate(ispec_is_surface_external_mesh(NSPEC_AB), &
           iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if (ier /= 0) stop 'error allocating array'
  nfaces_surface = 0

  ! collects MPI interfaces for detection
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh)
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) stop 'error allocating array'

  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,:) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,:)
  enddo
  call synchronize_all()

  call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                      ispec_is_surface_external_mesh, &
                      iglob_is_surface_external_mesh, &
                      nfaces_surface, &
                      num_interfaces_ext_mesh, &
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh, &
                      my_neighbors_ext_mesh, &
                      ibool_interfaces_ext_mesh_dummy )

  ! takes number of faces for top, free surface only
  if (MOVIE_TYPE == 1) then
    nfaces_surface = NSPEC2D_TOP
  endif

! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface,nfaces_surface_glob_ext_mesh)

  deallocate(ibool_interfaces_ext_mesh_dummy)

  end subroutine crm_setup_mesh_surface

