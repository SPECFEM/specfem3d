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

  use constants, only: myrank,PARALLEL_FAULT,IMAIN

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION, &
    STACEY_ABSORBING_CONDITIONS,SAVE_MESH_FILES,PML_CONDITIONS, &
    ANISOTROPY,APPROXIMATE_OCEAN_LOAD,OLSEN_ATTENUATION_RATIO, &
    ATTENUATION,USE_OLSEN_ATTENUATION, &
    ADIOS_FOR_MESH,SAVE_MOHO_MESH,ATTENUATION_f0_REFERENCE, &
    LOCAL_PATH

  use generate_databases_par, only: nspec => NSPEC_AB,nglob => NGLOB_AB, &
      ibool,xstore,ystore,zstore, &
      npointot, &
      nnodes_ext_mesh,nelmnts_ext_mesh, &
      nodes_coords_ext_mesh, elmnts_ext_mesh, &
      max_memory_size,num_interfaces_ext_mesh, &
      nibool_interfaces_ext_mesh, &
      nspec2D_xmin, nspec2D_xmax, &
      nspec2D_ymin, nspec2D_ymax, &
      NSPEC2D_BOTTOM, NSPEC2D_TOP, &
      ibelm_xmin, ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top, &
      nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
      nodes_ibelm_bottom,nodes_ibelm_top, &
      nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho

  use create_regions_mesh_ext_par

  use fault_generate_databases, only: fault_read_input,fault_setup, &
                          fault_save_arrays,fault_save_arrays_test, &
                          nnodes_coords_open,nodes_coords_open,ANY_FAULT_IN_THIS_PROC, &
                          ANY_FAULT

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
  call crm_ext_allocate_arrays(nspec,LOCAL_PATH, &
                               nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                               nspec2D_bottom,nspec2D_top,ANISOTROPY, &
                               nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh,ANY_FAULT_IN_THIS_PROC)

  ! if faults exist this reads nodes_coords_open
  call fault_read_input(prname)

  ! fills location and weights for Gauss-Lobatto-Legendre points, shape and derivations,
  ! returns jacobianstore,xixstore,...gammazstore
  ! and GLL-point locations in xstore,ystore,zstore
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up jacobian '
    call flush_IMAIN()
  endif
  if (ANY_FAULT_IN_THIS_PROC) then
   ! compute jacobians with fault open and *store needed for ibool.
    call crm_ext_setup_jacobian(xstore,ystore,zstore,nspec, &
                                nodes_coords_open, nnodes_coords_open,elmnts_ext_mesh,nelmnts_ext_mesh)
  else ! with no fault
    call crm_ext_setup_jacobian(xstore,ystore,zstore,nspec, &
                                nodes_coords_ext_mesh,nnodes_ext_mesh,elmnts_ext_mesh,nelmnts_ext_mesh)
  endif


  ! creates ibool index array for projection from local to global points
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
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
      write(IMAIN,*)
      write(IMAIN,*) '  ...resetting up jacobian in fault domains'
      call flush_IMAIN()
    endif
    if (ANY_FAULT_IN_THIS_PROC) then
      call crm_ext_setup_jacobian(xstore,ystore,zstore,nspec, &
                                  nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                  elmnts_ext_mesh,nelmnts_ext_mesh)
    endif
    ! at this point (xyz)store_unique are still open
    if (.not. PARALLEL_FAULT) then
      call fault_setup(ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                       xstore,ystore,zstore,nspec,nglob)
    endif
    ! this closes (xyz)store_unique
  endif


  ! sets up MPI interfaces between partitions
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...preparing MPI interfaces '
    call flush_IMAIN()
  endif
  call get_MPI_interface(nglob_unique,nspec,ibool)

  ! setting up parallel fault
  if (PARALLEL_FAULT .and. ANY_FAULT) then
    call synchronize_all()
    !at this point (xyz)store_unique are still open
    call fault_setup(ibool,nnodes_ext_mesh,nodes_coords_ext_mesh, &
                     xstore,ystore,zstore,nspec,nglob)
   ! this closes (xyz)store_unique
  endif

  ! sets up absorbing/free surface boundaries
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up absorbing boundaries'
    call flush_IMAIN()
  endif
  call get_absorbing_boundary(nspec,ibool, &
                              nodes_coords_ext_mesh,nnodes_ext_mesh, &
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                              nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                              nodes_ibelm_bottom,nodes_ibelm_top, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              nspec2D_bottom,nspec2D_top)

  ! sets up mesh surface
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up mesh surface'
    call flush_IMAIN()
  endif
  call crm_setup_mesh_surface()

  ! sets up up Moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...setting up Moho surface'
      call flush_IMAIN()
    endif
    call crm_setup_moho(nspec,nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )
  endif

  ! sets material velocities
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...determining velocity model'
    call flush_IMAIN()
  endif
  call get_model()

  ! sets up acoustic-elastic-poroelastic coupling surfaces
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...detecting acoustic-elastic-poroelastic surfaces '
    call flush_IMAIN()
  endif
  call get_coupling_surfaces(nspec,ibool)

  ! locates inner and outer elements
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...element inner/outer separation '
    call flush_IMAIN()
  endif
  call crm_setup_inner_outer_elemnts(nspec,ibool,SAVE_MESH_FILES)

  ! colors mesh if requested
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...element mesh coloring '
    call flush_IMAIN()
  endif
  call setup_color_perm(nspec,nglob,ibool,ANISOTROPY,SAVE_MESH_FILES)

  ! overwrites material parameters from external binary files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...external binary models '
    call flush_IMAIN()
  endif
  call get_model_binaries(nspec,LOCAL_PATH)

  ! calculates damping profiles and auxiliary coefficients on all C-PML points
  if (PML_CONDITIONS) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating C-PML damping profiles '
      call flush_IMAIN()
    endif
    call pml_set_local_dampingcoeff(xstore_unique,ystore_unique,zstore_unique)
  endif

  ! creates mass matrix
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...creating mass matrix '
    call flush_IMAIN()
  endif
  call create_mass_matrices(nglob_unique,nspec,ibool,PML_CONDITIONS,STACEY_ABSORBING_CONDITIONS)

  ! saves the binary mesh files
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...saving databases'
    call flush_IMAIN()
  endif
  !call create_name_database(prname,myrank,LOCAL_PATH)
  if (ADIOS_FOR_MESH) then
    call save_arrays_solver_ext_mesh_adios(nspec,ibool)
  else
    call save_arrays_solver_ext_mesh(nspec,ibool)
  endif

  ! saves faults
  if (ANY_FAULT) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving fault databases'
      call flush_IMAIN()
    endif
    ! for debugging
    !call fault_save_arrays_test(prname)

    ! saves fault arrays
    call fault_save_arrays(prname)
  endif

  ! saves moho surface
  if (SAVE_MOHO_MESH) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving Moho surfaces'
      call flush_IMAIN()
    endif
    call crm_save_moho()
  endif
  call synchronize_all()

  ! computes the approximate amount of memory needed to run the solver
  call memory_eval(nspec,nglob_unique,maxval(nibool_interfaces_ext_mesh),num_interfaces_ext_mesh, &
                   APPROXIMATE_OCEAN_LOAD,memory_size)

  call max_all_dp(memory_size, max_memory_size)

  ! checks the mesh, stability and resolved period
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...checking mesh resolution'
    call flush_IMAIN()
  endif
  call check_mesh_resolution(nspec,nglob_unique, &
                             ibool,xstore_unique,ystore_unique,zstore_unique, &
                             ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
                             kappastore,mustore,rhostore, &
                             phistore,tortstore,rhoarraystore,rho_vpI,rho_vpII,rho_vsI, &
                             -1.0d0,model_speed_max,min_resolved_period)

  ! saves binary mesh files for attenuation
  if (ATTENUATION) then
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...saving attenuation databases'
      call flush_IMAIN()
    endif
    call get_attenuation_model(nspec,USE_OLSEN_ATTENUATION,OLSEN_ATTENUATION_RATIO, &
                               mustore,rho_vs,kappastore,rho_vp, &
                               qkappa_attenuation_store,qmu_attenuation_store, &
                               ispec_is_elastic,min_resolved_period,prname,ATTENUATION_f0_REFERENCE)
  endif

  ! synchronizes processes, making sure everybody has finished
  call synchronize_all()

  ! cleanup
  deallocate(xixstore,xiystore,xizstore, &
             etaxstore,etaystore,etazstore, &
             gammaxstore,gammaystore,gammazstore,jacobianstore)
  deallocate(qkappa_attenuation_store,qmu_attenuation_store)
  deallocate(kappastore,mustore,rhostore,rho_vp,rho_vs)
  deallocate(rho_vpI,rho_vpII,rho_vsI)
  deallocate(rhoarraystore,kappaarraystore,etastore,phistore,tortstore,permstore)

  if (.not. SAVE_MOHO_MESH) deallocate(xstore_unique,ystore_unique,zstore_unique)

  if (ACOUSTIC_SIMULATION) deallocate(rmass_acoustic)
  if (ELASTIC_SIMULATION) deallocate(rmass)
  if (POROELASTIC_SIMULATION) deallocate(rmass_solid_poroelastic,rmass_fluid_poroelastic)

  if (STACEY_ABSORBING_CONDITIONS) then
    if (ELASTIC_SIMULATION) deallocate(rmassx,rmassy,rmassz)
    if (ACOUSTIC_SIMULATION) deallocate(rmassz_acoustic)
  endif

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'mesh regions done'
    call flush_IMAIN()
  endif

  end subroutine create_regions_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_ext_allocate_arrays(nspec,LOCAL_PATH, &
                                     nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                                     nspec2D_bottom,nspec2D_top,ANISOTROPY, &
                                     nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                     elmnts_ext_mesh,nelmnts_ext_mesh,ANY_FAULT_IN_THIS_PROC)

  use constants, only: myrank,NDIM,NDIM2D,NGLLX,NGLLY,NGLLZ,NGLLSQUARE,IMAIN,GAUSSALPHA,GAUSSBETA, &
    DO_IRREGULAR_ELEMENT_SEPARATION

  use generate_databases_par, only: STACEY_INSTEAD_OF_FREE_SURFACE,PML_INSTEAD_OF_FREE_SURFACE,BOTTOM_FREE_SURFACE, &
    NGNOD,NGNOD2D

  use create_regions_mesh_ext_par

  implicit none

  integer :: nspec
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
             nspec2D_bottom,nspec2D_top

  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  logical :: ANISOTROPY,ANY_FAULT_IN_THIS_PROC

  ! data from the external mesh
  integer :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh) :: nodes_coords_ext_mesh
  integer, dimension(NGNOD,nelmnts_ext_mesh) :: elmnts_ext_mesh

  ! local parameters
  integer :: ier,ispec,ia,i,j,k
  logical :: any_regular_elem
  double precision :: cube_edge_size_squared
  real, dimension(NGNOD) :: xelm_real,yelm_real,zelm_real

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    NGLLX   = ',NGLLX
    write(IMAIN,*) '    NGLLY   = ',NGLLY
    write(IMAIN,*) '    NGLLZ   = ',NGLLZ
    write(IMAIN,*) '    NGNOD   = ',NGNOD
    write(IMAIN,*) '    NGNOD2D = ',NGNOD2D
    write(IMAIN,*)
    write(IMAIN,*) '    main process setup: '
    write(IMAIN,*) '    nspec           = ',nspec
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! attenuation
  allocate(qkappa_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 724')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  qkappa_attenuation_store(:,:,:,:) = 9999.9_CUSTOM_REAL

  allocate(qmu_attenuation_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 725')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  qmu_attenuation_store(:,:,:,:) = 9999.9_CUSTOM_REAL

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,LOCAL_PATH)

  ! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 726')
  if (ier /= 0) stop 'error allocating array xigll etc.'
  xigll(:) = 0.d0; yigll(:) = 0.d0; zigll(:) = 0.d0

  ! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 727')
  if (ier /= 0) stop 'error allocating array wxgll etc.'
  wxgll(:) = 0.0; wygll(:) = 0.d0; wzgll(:) = 0.d0

  ! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 728')
  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 729')
  if (ier /= 0) stop 'error allocating array shape3D etc.'
  shape3D(:,:,:,:) = 0.d0; dershape3D(:,:,:,:,:) = 0.d0

  ! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 730')
  allocate(shape2D_y(NGNOD2D,NGLLX,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 731')
  allocate(shape2D_bottom(NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 732')
  allocate(shape2D_top(NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 733')
  if (ier /= 0) stop 'error allocating array shape2D_x etc.'
  shape2D_x(:,:,:) = 0.d0; shape2D_y(:,:,:) = 0.d0
  shape2D_bottom(:,:,:) = 0.d0; shape2D_top(:,:,:) = 0.d0

  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 734')
  allocate(dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 735')
  allocate(dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 736')
  allocate(dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 737')
  if (ier /= 0) stop 'error allocating array dershape2D_x etc.'
  dershape2D_x(:,:,:,:) = 0.d0; dershape2D_y(:,:,:,:) = 0.d0
  dershape2D_bottom(:,:,:,:) = 0.d0; dershape2D_top(:,:,:,:) = 0.d0

  allocate(wgllwgll_xy(NGLLX,NGLLY),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 738')
  allocate(wgllwgll_xz(NGLLX,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 739')
  allocate(wgllwgll_yz(NGLLY,NGLLZ),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 740')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  wgllwgll_xy(:,:) = 0.d0; wgllwgll_xz(:,:) = 0.d0; wgllwgll_yz(:,:) = 0.d0

  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX,NGLLY,NGLLZ)

  ! get the 2-D shape functions
  call get_shape2D(shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ,NGNOD,NGNOD2D)
  call get_shape2D(shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)
  call get_shape2D(shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY,NGNOD,NGNOD2D)

  ! 2D weights
  do j = 1,NGLLY
    do i = 1,NGLLX
      wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo
  do k = 1,NGLLZ
    do i = 1,NGLLX
      wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo
  do k = 1,NGLLZ
    do j = 1,NGLLY
      wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

  ! model parameters
  ! (full arrays needed for reading in acoustic/elastic/poroelastic velocity models)

  ! acoustic/elastic array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 743')
  allocate(kappastore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 744')
  allocate(mustore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 745')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  rhostore(:,:,:,:) = 0.0_CUSTOM_REAL; kappastore(:,:,:,:) = 0.0_CUSTOM_REAL; mustore(:,:,:,:) = 0.0_CUSTOM_REAL

  ! Stacey
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 741')
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 742')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  rho_vp(:,:,:,:) = 0.0_CUSTOM_REAL; rho_vs(:,:,:,:) = 0.0_CUSTOM_REAL

!EB EB April 2018 : we should find a way to know if there are any poroelastic
! element before these costly allocations. The commented part does notwork
! because POROELASTIC_SIMULATION flag is set after (in get_model)
!  if (POROELASTIC_SIMULATION) then
!    NSPEC_PORO = nspec
!  else
!    NSPEC_PORO = 1
!  endif

  ! array with poroelastic model
  NSPEC_PORO = nspec
  allocate(rhoarraystore(2,NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 746')
  allocate(kappaarraystore(3,NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 747')
  allocate(etastore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 748')
  allocate(tortstore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 749')
  allocate(phistore(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 750')
  allocate(rho_vpI(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 751')
  allocate(rho_vpII(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 752')
  allocate(rho_vsI(NGLLX,NGLLY,NGLLZ,NSPEC_PORO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 753')
  allocate(permstore(6,NGLLX,NGLLY,NGLLZ,NSPEC_PORO), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 754')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  ! mesh arrays
  ! get the number of regular and irregular elements
  any_regular_elem = .false.

  allocate(irregular_element_number(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 755')
  irregular_element_number(:) = 0

  ! distinguishes between regular and irregular elements (based on cube shape)
  if (DO_IRREGULAR_ELEMENT_SEPARATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    separating regular/irregular element shapes'
      call flush_IMAIN()
    endif

    ! checks each element shape
    irregular_element_number(:) = 0
    nspec_irregular = nspec
    cube_edge_size_squared = 1.e9

    do ispec = 1, nspec
      ! gets element corner positions
      do ia = 1,NGNOD
        xelm_real(ia) = real(nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec)))
        yelm_real(ia) = real(nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec)))
        zelm_real(ia) = real(nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec)))
      enddo
      ! checks if element is regular (is a cube)
      call check_element_regularity(xelm_real,yelm_real,zelm_real,any_regular_elem,cube_edge_size_squared, &
                                    nspec_irregular,ispec,nspec,irregular_element_number,ANY_FAULT_IN_THIS_PROC)
    enddo
  else
    ! default case of previous versions: don't separate and assume each element to be irregular in shape
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    assuming general irregular element shapes'
      call flush_IMAIN()
    endif
    ! assigns each element as being irregular (value /= 0)
    nspec_irregular = nspec
    do ispec = 1,nspec
      irregular_element_number(ispec) = ispec
    enddo
  endif
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    nspec regular   = ',nspec - nspec_irregular
    write(IMAIN,*) '    nspec irregular = ',nspec_irregular
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (nspec_irregular > 0) then
    allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 766')
    allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 767')
    allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 768')
    allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 769')
    allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 770')
    allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 771')
    allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 772')
    allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 773')
    allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 774')
    allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,nspec_irregular),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 775')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  else
    ! dummy arrays
    allocate(xixstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 756')
    allocate(xiystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 757')
    allocate(xizstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 758')
    allocate(etaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 759')
    allocate(etaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 760')
    allocate(etazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 761')
    allocate(gammaxstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 762')
    allocate(gammaystore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 763')
    allocate(gammazstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 764')
    allocate(jacobianstore(1,1,1,1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 765')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  endif
  xixstore(:,:,:,:) = 0.0_CUSTOM_REAL; xiystore(:,:,:,:) = 0.0_CUSTOM_REAL; xizstore(:,:,:,:) = 0.0_CUSTOM_REAL
  etaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; etaystore(:,:,:,:) = 0.0_CUSTOM_REAL; etazstore(:,:,:,:) = 0.0_CUSTOM_REAL
  gammaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; gammaystore(:,:,:,:) = 0.0_CUSTOM_REAL; gammazstore(:,:,:,:) = 0.0_CUSTOM_REAL
  jacobianstore(:,:,:,:) = 0.0_CUSTOM_REAL

  ! absorbing boundary
  ! absorbing faces
  num_abs_boundary_faces = nspec2D_xmin + nspec2D_xmax + nspec2D_ymin + nspec2D_ymax + nspec2D_bottom
  ! adds faces of free surface if it also absorbs
  if (STACEY_INSTEAD_OF_FREE_SURFACE .or. PML_INSTEAD_OF_FREE_SURFACE) then
     num_abs_boundary_faces = num_abs_boundary_faces + nspec2D_top
  endif

  ! subtract bottom surface from Stacey condition when free surface
  if (BOTTOM_FREE_SURFACE)  num_abs_boundary_faces = num_abs_boundary_faces -  nspec2D_bottom

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    absorbing boundary faces:'
    write(IMAIN,*) '      num_abs_boundary_faces = ',num_abs_boundary_faces
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  if (num_abs_boundary_faces > 0) then
    allocate(abs_boundary_ispec(num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 776')
    allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 777')
    allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 778')
    allocate(abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 779')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  else
    ! dummy allocations
    allocate(abs_boundary_ispec(1),abs_boundary_ijk(1,1,1), &
             abs_boundary_jacobian2Dw(1,1),abs_boundary_normal(1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy arrays'
  endif
  abs_boundary_ispec(:) = 0; abs_boundary_ijk(:,:,:) = 0
  abs_boundary_jacobian2Dw(:,:) = 0.0; abs_boundary_normal(:,:,:) = 0.0

  ! free surface faces
  num_free_surface_faces = nspec2D_top

   ! add bottom surface to free surface condition
  if (BOTTOM_FREE_SURFACE) then
    num_free_surface_faces = num_free_surface_faces + nspec2D_bottom
    if (STACEY_INSTEAD_OF_FREE_SURFACE)  num_free_surface_faces = num_free_surface_faces - nspec2D_top
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    free surface faces:'
    write(IMAIN,*) '      num_free_surface_faces = ',num_free_surface_faces
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates arrays to store info for each face (assumes NGLLX=NGLLY=NGLLZ)
  if (num_free_surface_faces > 0) then
    allocate(free_surface_ispec(num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 780')
    allocate(free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 781')
    allocate(free_surface_jacobian2Dw(NGLLSQUARE,num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 782')
    allocate(free_surface_normal(NDIM,NGLLSQUARE,num_free_surface_faces),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 783')
    if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  else
    ! dummy allocation
    allocate(free_surface_ispec(1),free_surface_ijk(1,1,1), &
             free_surface_jacobian2Dw(1,1),free_surface_normal(1,1,1),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy arrays'
  endif
  free_surface_ispec(:) = 0; free_surface_ijk(:,:,:) = 0
  free_surface_jacobian2Dw(:,:) = 0.0; free_surface_normal(:,:,:) = 0.0

  ! array with anisotropy
  if (ANISOTROPY) then
    NSPEC_ANISO = nspec
  else
    NSPEC_ANISO = 1
  endif
  allocate(c11store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 784')
  allocate(c12store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 785')
  allocate(c13store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 786')
  allocate(c14store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 787')
  allocate(c15store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 788')
  allocate(c16store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 789')
  allocate(c22store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 790')
  allocate(c23store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 791')
  allocate(c24store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 792')
  allocate(c25store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 793')
  allocate(c26store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 794')
  allocate(c33store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 795')
  allocate(c34store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 796')
  allocate(c35store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 797')
  allocate(c36store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 798')
  allocate(c44store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 799')
  allocate(c45store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 800')
  allocate(c46store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 801')
  allocate(c55store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 802')
  allocate(c56store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 803')
  allocate(c66store(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 804')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  c11store(:,:,:,:) = 0.0_CUSTOM_REAL; c12store(:,:,:,:) = 0.0_CUSTOM_REAL; c13store(:,:,:,:) = 0.0_CUSTOM_REAL
  c14store(:,:,:,:) = 0.0_CUSTOM_REAL; c15store(:,:,:,:) = 0.0_CUSTOM_REAL; c16store(:,:,:,:) = 0.0_CUSTOM_REAL
  c22store(:,:,:,:) = 0.0_CUSTOM_REAL; c23store(:,:,:,:) = 0.0_CUSTOM_REAL; c24store(:,:,:,:) = 0.0_CUSTOM_REAL
  c25store(:,:,:,:) = 0.0_CUSTOM_REAL; c26store(:,:,:,:) = 0.0_CUSTOM_REAL; c33store(:,:,:,:) = 0.0_CUSTOM_REAL
  c34store(:,:,:,:) = 0.0_CUSTOM_REAL; c35store(:,:,:,:) = 0.0_CUSTOM_REAL; c36store(:,:,:,:) = 0.0_CUSTOM_REAL
  c44store(:,:,:,:) = 0.0_CUSTOM_REAL; c45store(:,:,:,:) = 0.0_CUSTOM_REAL; c46store(:,:,:,:) = 0.0_CUSTOM_REAL
  c55store(:,:,:,:) = 0.0_CUSTOM_REAL; c56store(:,:,:,:) = 0.0_CUSTOM_REAL; c66store(:,:,:,:) = 0.0_CUSTOM_REAL

  ! material flags
  allocate(ispec_is_acoustic(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 805')
  allocate(ispec_is_elastic(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 806')
  allocate(ispec_is_poroelastic(nspec), stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 807')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  ispec_is_acoustic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! initializes Moho surface
  NSPEC2D_MOHO = 0

  end subroutine crm_ext_allocate_arrays


!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_ext_setup_jacobian(xstore,ystore,zstore,nspec, &
                                    nodes_coords_ext_mesh,nnodes_ext_mesh, &
                                    elmnts_ext_mesh,nelmnts_ext_mesh)

  use constants, only: myrank,NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES
  use generate_databases_par, only: NGNOD

  use create_regions_mesh_ext_par, only: shape3D,dershape3D, &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,jacobianstore, &
    irregular_element_number,xix_regular,jacobian_regular

  implicit none

! number of spectral elements in each block
  integer,intent(in) :: nspec

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: xstore,ystore,zstore

! data from the external mesh
  integer,intent(in) :: nnodes_ext_mesh,nelmnts_ext_mesh
  double precision, dimension(NDIM,nnodes_ext_mesh),intent(in) :: nodes_coords_ext_mesh
  integer, dimension(NGNOD,nelmnts_ext_mesh),intent(in) :: elmnts_ext_mesh

! local parameters
  integer :: ispec,ia,ispec_irreg
  logical :: any_regular_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix_reg,xiy_reg,xiz_reg,etax_reg,etay_reg,etaz_reg, &
                                                          gammax_reg,gammay_reg,gammaz_reg,jacobian_reg
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision,parameter :: threshold_zero = 1.e-25

  ! debug
  logical, parameter :: DEBUG_ELEMENT = .false.
  character(len=MAX_STRING_LEN) :: filename

  ! point locations
  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  any_regular_elem = .false.

  do ispec = 1, nspec
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo

    ! CUBIT should provide a mesh ordering such that the 3D jacobian is defined
    ! (otherwise mesh would be degenerated)
    ispec_irreg = irregular_element_number(ispec)

    ! irregular_element_number is 0 only if the element is regular
    if (ispec_irreg /= 0) then
      ! note: we pass arrays by x**(1,1,1,ispec_irreg) to avoid possible array copies and since Fortran uses pointers,
      !       this will point to the first entry for the ispec element
      call calc_jacobian(myrank,xixstore(1,1,1,ispec_irreg),xiystore(1,1,1,ispec_irreg),xizstore(1,1,1,ispec_irreg), &
                         etaxstore(1,1,1,ispec_irreg),etaystore(1,1,1,ispec_irreg),etazstore(1,1,1,ispec_irreg), &
                         gammaxstore(1,1,1,ispec_irreg),gammaystore(1,1,1,ispec_irreg),gammazstore(1,1,1,ispec_irreg), &
                         jacobianstore(1,1,1,ispec_irreg),xelm,yelm,zelm,dershape3D)
    else
      ! sets flag for regular elements
      any_regular_elem = .true.
    endif
    call calc_coords(xstore(1,1,1,ispec),ystore(1,1,1,ispec),zstore(1,1,1,ispec), &
                     xelm,yelm,zelm,shape3D)

    !debug
    if (DEBUG_ELEMENT) then
      if (myrank == 0 .and. ispec == 1) then
        write(filename,'(a,i6.6,a)') trim(OUTPUT_FILES)//'/proc',myrank,'_debug_element'
        call write_VTK_data_points_elem(NGNOD,xelm,yelm,zelm,dble(jacobianstore(1,1,1,ispec_irreg)),filename)
        print *,'  written out:',trim(filename)
      endif
    endif
  enddo

  ! get xix derivative and jacobian on a regular element
  if (any_regular_elem) then
    ! find a regular element
    ispec = 1
    do while (irregular_element_number(ispec) /= 0)
      ispec = ispec + 1
    enddo

    ! checks index
    if (ispec < 1 .or. ispec > nspec) then
      print *,'Error: element ',ispec,' - regular element index out of bounds ',nspec
      print *,'  regular elements indices:'
      do ispec = 1,nspec
        if (irregular_element_number(ispec) == 0) then
          print *,'  ispec index ',ispec,' has ispec_irreg = ',irregular_element_number(ispec)
        endif
      enddo
      call exit_mpi(myrank,'Invalid regular element index')
    endif

    ! gets corner positions of regular element
    do ia = 1,NGNOD
      xelm(ia) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(ia,ispec))
      yelm(ia) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(ia,ispec))
      zelm(ia) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(ia,ispec))
    enddo

    ! jacobian and derivatives of mapping
    call calc_jacobian(myrank,xix_reg,xiy_reg,xiz_reg, &
                       etax_reg,etay_reg,etaz_reg, &
                       gammax_reg,gammay_reg,gammaz_reg, &
                       jacobian_reg,xelm,yelm,zelm,dershape3D)

    ! debug
    !print *,'debug: rank ',myrank,' regular elements indices:',any_regular_elem
    !do ispec = 1,nspec
    !  if (irregular_element_number(ispec) == 0) then
    !    print *,'debug: rank',myrank,'  ispec index ',ispec,' has ispec_irreg = ',irregular_element_number(ispec)
    !  endif
    !enddo
    ! only xix == etay == gammaz are non-zero for regular elements
    ! debug
    !print *,'debug: xix    ',xix_reg(1,1,1),xiy_reg(1,1,1),xiz_reg(1,1,1)
    !print *,'debug: etax   ',etax_reg(1,1,1),etay_reg(1,1,1),etaz_reg(1,1,1)
    !print *,'debug: gammax ',gammax_reg(1,1,1),gammay_reg(1,1,1),gammaz_reg(1,1,1)

    ! check
    if (abs(xix_reg(1,1,1) - etay_reg(1,1,1)) > threshold_zero) then
      print *,'Error: element ',ispec,' - regular element should have xix == etay ',xix_reg(1,1,1),etay_reg(1,1,1)
      call exit_mpi(myrank,'Invalid regular element xix/etay')
    endif
    if (abs(xix_reg(1,1,1) - gammaz_reg(1,1,1)) > threshold_zero) then
      print *,'Error: element ',ispec,' - regular element should have xix == gammaz ',xix_reg(1,1,1),gammaz_reg(1,1,1)
      call exit_mpi(myrank,'Invalid regular element xix/gammaz')
    endif

    ! saves regular values
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

  use generate_databases_par, only: NGLLX,NGLLY,NGLLZ,NDIM,IMAIN
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
  double precision :: x_min,x_max,x_min_all,x_max_all

! allocate memory for arrays
  allocate(locval(npointot),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 808')
  locval = 0

  allocate(ifseg(npointot),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 809')
  ifseg = .false.

  allocate(xp(npointot),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 810')
  xp = 0.d0

  allocate(yp(npointot),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 811')
  yp = 0.d0

  allocate(zp(npointot),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 812')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  zp = 0.d0

  ! creates temporary global point arrays
  do ispec = 1,nspec
    ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
    ilocnum = 0
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
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
  call min_all_all_dp(x_min,x_min_all)
  call max_all_all_dp(x_max,x_max_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     creating ibool indexing     : x min/max = ',sngl(x_min_all),'/',sngl(x_max_all)
    call flush_IMAIN()
  endif

  ! gets ibool indexing from local (GLL points) to global points
  call get_global(npointot,xp,yp,zp,ibool,locval,ifseg,nglob,x_min_all,x_max_all)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     creating indirect addressing: nglob = ',nglob
    call flush_IMAIN()
  endif

  ! we can create a new indirect addressing to reduce cache misses
  call get_global_indirect_addressing(nspec,nglob,ibool)

  ! cleanup
  deallocate(xp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(yp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(zp,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(locval,stat=ier); if (ier /= 0) stop 'error in deallocate'
  deallocate(ifseg,stat=ier); if (ier /= 0) stop 'error in deallocate'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     creating unique point locations'
    call flush_IMAIN()
  endif

  ! unique global point locations
  nglob_unique = nglob              ! note: nglob => NGLOB_AB

  ! coordinates for global points
  allocate(xstore_unique(nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 813')
  xstore_unique(:) = 0._CUSTOM_REAL

  allocate(ystore_unique(nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 814')
  ystore_unique(:) = 0._CUSTOM_REAL

  allocate(zstore_unique(nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 815')
  if (ier /= 0) stop 'error in allocate'
  zstore_unique(:) = 0._CUSTOM_REAL

  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglobnum = ibool(i,j,k,ispec)
          xstore_unique(iglobnum) = real(xstore(i,j,k,ispec),kind=CUSTOM_REAL)
          ystore_unique(iglobnum) = real(ystore(i,j,k,ispec),kind=CUSTOM_REAL)
          zstore_unique(iglobnum) = real(zstore(i,j,k,ispec),kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo

  end subroutine crm_ext_setup_indexing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_moho(nspec,nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                            nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,IMAIN,NDIM,NGLLSQUARE,NGNOD2D_FOUR_CORNERS
  use generate_databases_par, only: NGNOD2D
  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nspec

  integer :: nspec2D_moho_ext
  integer, dimension(nspec2D_moho_ext) :: ibelm_moho
  integer, dimension(NGNOD2D,nspec2D_moho_ext) :: nodes_ibelm_moho

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
  allocate(iglob_is_surface(nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 816')
  allocate(iglob_normals(NDIM,nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 817')
  if (ier /= 0) stop 'error allocating array iglob_is_surface'

  iglob_is_surface(:) = 0
  iglob_normals(:,:) = 0._CUSTOM_REAL

  ! loops over given moho surface elements
  do ispec2D = 1,nspec2D_moho_ext

    ! gets element id
    ispec = ibelm_moho(ispec2D)

    ! looks for i,j,k indices of GLL points on boundary face
    ! determines element face by given CUBIT corners
    ! (note: uses point locations rather than point indices to find the element face,
    !            because the indices refer no more to the newly indexed ibool array )
    do icorner = 1,NGNOD2D_FOUR_CORNERS
      xcoord(icorner) = nodes_coords_ext_mesh(1,nodes_ibelm_moho(icorner,ispec2D))
      ycoord(icorner) = nodes_coords_ext_mesh(2,nodes_ibelm_moho(icorner,ispec2D))
      zcoord(icorner) = nodes_coords_ext_mesh(3,nodes_ibelm_moho(icorner,ispec2D))
    enddo

    ! sets face id of reference element associated with this face
    call get_element_face_id(ispec,xcoord,ycoord,zcoord, &
                             ibool,nspec,nglob_unique, &
                             xstore_unique,ystore_unique,zstore_unique, &
                             iface)

    ! ijk indices of GLL points for face id
    call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLZ)

    ! weighted jacobian and normal
    call get_jacobian_boundary_face(nspec, &
                                    xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                    dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                    ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

    ! normal convention: points away from element
    ! switch normal direction if necessary
    do j = 1,NGLLY
      do i = 1,NGLLX
          call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                       ibool,nspec,nglob_unique, &
                                       xstore_unique,ystore_unique,zstore_unique, &
                                       normal_face(:,i,j))
      enddo
    enddo

    ! stores information on global points on moho surface
    igll = 0
    do j = 1,NGLLY
      do i = 1,NGLLX
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

  allocate(ibelm_moho_bot(NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 818')
  allocate(ibelm_moho_top(NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 819')
  allocate(normal_moho_top(NDIM,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 820')
  allocate(normal_moho_bot(NDIM,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 821')
  allocate(ijk_moho_bot(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 822')
  allocate(ijk_moho_top(3,NGLLSQUARE,NSPEC2D_MOHO),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 823')
  if (ier /= 0) stop 'error allocating ibelm_moho_bot'

  ibelm_moho_bot = 0
  ibelm_moho_top = 0

  ! element flags
  allocate(is_moho_top(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 824')
  allocate(is_moho_bot(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 825')
  if (ier /= 0) stop 'error allocating is_moho_top'
  is_moho_top = .false.
  is_moho_bot = .false.

  ! finds spectral elements with moho surface
  imoho_top = 0
  imoho_bot = 0
  do ispec = 1,nspec
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
          xcoord(icorner) = xstore_unique(iglob)
          ycoord(icorner) = ystore_unique(iglob)
          zcoord(icorner) = zstore_unique(iglob)
        endif
      enddo

      ! stores moho informations
      if (counter == NGNOD2D_FOUR_CORNERS) then

        ! gets face GLL points i,j,k indices from element face
        call get_element_face_gll_indices(iface,ijk_face,NGLLX,NGLLY)

        ! re-computes face infos
        ! weighted jacobian and normal
        call get_jacobian_boundary_face(nspec, &
                                        xstore_unique,ystore_unique,zstore_unique,ibool,nglob_unique, &
                                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                                        wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                        ispec,iface,jacobian2Dw_face,normal_face,NGLLX,NGLLZ,NGNOD2D)

        ! normal convention: points away from element
        ! switch normal direction if necessary
        do j = 1,NGLLZ
          do i = 1,NGLLX
            call get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                         ibool,nspec,nglob_unique, &
                                         xstore_unique,ystore_unique,zstore_unique, &
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
                                             ibool,nspec,nglob_unique, &
                                             xstore_unique,ystore_unique,zstore_unique, &
                                             normal,idirect)

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

  ! note: surface e.g. could be at the free-surface and have no top elements etc.
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

  subroutine crm_setup_inner_outer_elemnts(nspec,ibool,SAVE_MESH_FILES)

! locates inner and outer elements

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,IMAIN

  use shared_parameters, only: ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, POROELASTIC_SIMULATION

  ! MPI interfaces
  use generate_databases_par, only: num_interfaces_ext_mesh, &
    nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh

  use create_regions_mesh_ext_par

  implicit none

  integer,intent(in) :: nspec
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  logical,intent(in) :: SAVE_MESH_FILES

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
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 826')
  if (ier /= 0) stop 'error allocating array ispec_is_inner'

  ! temporary array
  allocate(iglob_is_inner(nglob_unique),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 827')
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
    call write_VTK_data_elem_l(nspec,nglob_unique, &
                               xstore_unique,ystore_unique,zstore_unique,ibool, &
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 828')
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 829')
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 830')
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 831')
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 832')
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
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 833')
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
  allocate(ispec_is_surface_external_mesh(NSPEC_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 834')
  allocate(iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 835')
  if (ier /= 0) stop 'error allocating array'
  ispec_is_surface_external_mesh(:) = .false.; iglob_is_surface_external_mesh(:) = .false.
  nfaces_surface = 0

  ! collects MPI interfaces for detection
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh)
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 836')
  if (ier /= 0) stop 'error allocating array'
  ibool_interfaces_ext_mesh_dummy(:,:) = 0

  do i = 1,num_interfaces_ext_mesh
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

