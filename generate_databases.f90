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
!
! United States and French Government Sponsorship Acknowledged.
!
!=============================================================================!
!                                                                             !
!  generate_databases produces a spectral element grid                        !
!  for a local or regional model.                                             !
!  The mesher uses the UTM projection                                         !
!                                                                             !
!=============================================================================!
!
! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiPoKoTr04,
! author = {Qinya Liu and Jascha Polet and Dimitri Komatitsch and Jeroen Tromp},
! title = {Spectral-element moment tensor inversions for earthquakes in {S}outhern {C}alifornia},
! journal={Bull. Seismol. Soc. Am.},
! year = {2004},
! volume = {94},
! pages = {1748-1761},
! number = {5},
! doi = {10.1785/012004038}}
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoLiTrSuStSh04,
! author={Dimitri Komatitsch and Qinya Liu and Jeroen Tromp and Peter S\"{u}ss
!   and Christiane Stidham and John H. Shaw},
! year=2004,
! title={Simulations of Ground Motion in the {L}os {A}ngeles {B}asin
!   based upon the Spectral-Element Method},
! journal={Bull. Seism. Soc. Am.},
! volume=94,
! number=1,
! pages={187-206}}
!
! and/or another article from http://web.univ-pau.fr/~dkomati1/publications.html
!
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
!
! Reference frame - convention:
! ----------------------------
!
! The code uses the following convention for the reference frame:
!
!  - X axis is East
!  - Y axis is North
!  - Z axis is up
!
! Note that this convention is different from both the Aki-Richards convention
! and the Harvard CMT convention.
!
! Let us recall that the Aki-Richards convention is:
!
!  - X axis is North
!  - Y axis is East
!  - Z axis is down
!
! and that the Harvard CMT convention is:
!
!  - X axis is South
!  - Y axis is East
!  - Z axis is up
!
! To report bugs or suggest improvements to the code, please send an email
! to Jeroen Tromp <jtromp AT princeton.edu> and/or use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! MPI v. 2.0 "SESAME" (Spectral ElementS on Any MEsh), Fall 2009:
! Dimitri Komatitsch, Nicolas Le Goff, Roland Martin and Pieyre Le Loher, University of Pau, France,
! Jeroen Tromp and the Princeton group of developers, Princeton University, USA,
! and Emanuele Casarotti, INGV Roma, Italy:
!  support for CUBIT meshes decomposed by SCOTCH, METIS or ZOLTAN;
!  much faster solver using Michel Deville's inlined matrix products.
!
! MPI v. 1.4 Dimitri Komatitsch, University of Pau, Qinya Liu and others, Caltech, September 2006:
!  better adjoint and kernel calculations, faster and better I/Os
!  on very large systems, many small improvements and bug fixes
!
! MPI v. 1.3 Dimitri Komatitsch, University of Pau, and Qinya Liu, Caltech, July 2005:
!  serial version, regular mesh, adjoint and kernel calculations, ParaView support
!
! MPI v. 1.2 Min Chen and Dimitri Komatitsch, Caltech, July 2004:
!  full anisotropy, volume movie
!
! MPI v. 1.1 Dimitri Komatitsch, Caltech, October 2002: Zhu's Moho map, scaling
!  of Vs with depth, Hauksson's regional model, attenuation, oceans, movies
!
! MPI v. 1.0 Dimitri Komatitsch, Caltech, May 2002: first MPI version based on global code

!
!-------------------------------------------------------------------------------------------------
!

  module generate_databases_par

  implicit none

  include "constants.h"

! number of spectral elements in each block
  integer nspec,npointot

! local to global indexing array
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer :: myrank,sizeprocs,ier

! use integer array to store topography values
  integer :: UTM_PROJECTION_ZONE
  logical :: SUPPRESS_UTM_PROJECTION
  integer :: NX_TOPO,NY_TOPO
  double precision :: ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  character(len=100) :: topo_file
  integer, dimension(:,:), allocatable :: itopo_bathy
  
! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

! parameters read from parameter file
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE
  integer :: NSOURCES

  double precision :: DT,HDUR_MOVIE

  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS, SAVE_FORWARD
  logical :: ANISOTROPY,ABSORBING_CONDITIONS,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical :: MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES
  integer :: NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_READ_ADJSRC

  character(len=256) OUTPUT_FILES,LOCAL_PATH

! parameters deduced from parameters read from file
  integer :: NPROC

! static memory size that will be needed by the solver
  double precision :: max_static_memory_size,max_static_memory_size_request

! this for all the regions
  integer NSPEC_AB,NGLOB_AB
  
  integer NSPEC2D_BOTTOM,NSPEC2D_TOP
  
  double precision min_elevation,max_elevation
  double precision min_elevation_all,max_elevation_all

! for Databases of external meshes
  character(len=256) prname
  integer :: dummy_node
  integer :: dummy_elmnt
  integer :: ispec, inode, num_interface,ie,imat,iface,icorner
  integer :: nnodes_ext_mesh, nelmnts_ext_mesh
  integer  :: num_interfaces_ext_mesh
  integer  :: max_interface_size_ext_mesh
  integer  :: nmat_ext_mesh, nundefMat_ext_mesh   
  integer, dimension(:), allocatable  :: my_neighbours_ext_mesh
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(:,:,:), allocatable  :: my_interfaces_ext_mesh
  integer, dimension(:,:), allocatable  :: ibool_interfaces_ext_mesh
  integer, dimension(:), allocatable  :: nibool_interfaces_ext_mesh
  double precision, dimension(:,:), allocatable :: nodes_coords_ext_mesh
  integer, dimension(:,:), allocatable :: elmnts_ext_mesh
  integer, dimension(:,:), allocatable :: mat_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_dummy

! boundaries and materials
  integer  :: ispec2D, boundary_number
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom_ext, nspec2D_top_ext
  character (len=30), dimension(:,:), allocatable :: undef_mat_prop   
  integer, dimension(:), allocatable  :: ibelm_xmin,ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top
  integer, dimension(:,:), allocatable  :: nodes_ibelm_xmin,nodes_ibelm_xmax, &
              nodes_ibelm_ymin, nodes_ibelm_ymax, nodes_ibelm_bottom, nodes_ibelm_top
  double precision, dimension(:,:), allocatable :: materials_ext_mesh 

! moho (optional)  
  integer :: nspec2D_moho_ext
  integer, dimension(:), allocatable  :: ibelm_moho
  integer, dimension(:,:), allocatable  :: nodes_ibelm_moho
    
! number of points per spectral element
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

  integer :: nglob,nglob_total,nspec_total

  integer,dimension(:),allocatable :: ispec_is_surface_external_mesh,iglob_is_surface_external_mesh
  integer :: nfaces_surface_ext_mesh,nfaces_surface_glob_ext_mesh
  
  end module generate_databases_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine generate_databases

  use generate_databases_par
  implicit none
  
! sizeprocs returns number of processes started (should be equal to NPROC).
! myrank is the rank of each process, between 0 and NPROC-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

! get MPI starting time
  time_start = wtime()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '******************************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher - f90 version ***'
    write(IMAIN,*) '******************************************'
    write(IMAIN,*)
  endif

! read the parameter file
  call gd_read_parameters()
      
! makes sure processes are synchronized  
  call sync_all()
  
! reads topography and bathymetry file
  call gd_read_topography()
  
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! reads Databases files
  call gd_read_partition_files()

! external mesh creation
  call gd_setup_mesh()

! finalize mesher
  call gd_finalize()
  
  end subroutine generate_databases
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_read_parameters

! reads and checks user input parameters

  use generate_databases_par
  implicit none

! reads DATA/Par_file 
  call read_parameter_file( NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        OCEANS,ANISOTROPY,ABSORBING_CONDITIONS, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC)

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) then
    if( myrank == 0 ) then
      write(IMAIN,*) 'error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'error: number of processors actually run on: ',sizeprocs    
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

! there would be a problem with absorbing boundaries for different NGLLX,NGLLY,NGLLZ values
! just to be sure for now..
  if( ABSORBING_CONDITIONS ) then
    if( NGLLX /= NGLLY .and. NGLLY /= NGLLZ ) &
      call exit_MPI(myrank,'must have NGLLX = NGLLY = NGLLZ for external meshes')
  endif

! info about external mesh simulation
! nlegoff -- should be put in compute_parameters and read_parameter_file for clarity
! chris -- once the steps in decompose_mesh_SCOTCH are integrated into generate_database.f90,
! NPROC will be known

  if(myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'
    write(IMAIN,*)
    write(IMAIN,*) 'NGLLX = ',NGLLX
    write(IMAIN,*) 'NGLLY = ',NGLLY
    write(IMAIN,*) 'NGLLZ = ',NGLLZ

    write(IMAIN,*)
    write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD,' control nodes'
    write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D,' control nodes'
    write(IMAIN,*)
  endif

! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) &
    call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

  if(NGNOD /= 8) call exit_MPI(myrank,'number of control nodes must be 8')
  if(NGNOD2D /= 4) call exit_MPI(myrank,'elements with 8 points should have NGNOD2D = 4')

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) call exit_MPI(myrank,'number of SLS must be 3')

  ! exclusive movie flags
  if( EXTERNAL_MESH_MOVIE_SURFACE .or. EXTERNAL_MESH_CREATE_SHAKEMAP ) then  
    MOVIE_SURFACE = .false.
    CREATE_SHAKEMAP = .false.
  endif


  if(myrank == 0) then
! chris: I am not sure if we should suppress the following. topography should appear in the external mesh
! leave it for now

    write(IMAIN,*)
    if(SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'suppressing UTM projection'
    else
      write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
    endif

    write(IMAIN,*)
    if(ATTENUATION) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if(USE_OLSEN_ATTENUATION) then
        write(IMAIN,*) 'using Olsen''s attenuation'
      else
        write(IMAIN,*) 'not using Olsen''s attenuation'
      endif
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    if(ANISOTROPY) then
      write(IMAIN,*) 'incorporating anisotropy'
    else
      write(IMAIN,*) 'no anisotropy'
    endif

    write(IMAIN,*)
    if(OCEANS) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)

  endif

  end subroutine gd_read_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_read_topography

! reads in topography files

  use generate_databases_par
  implicit none

  allocate(itopo_bathy(NX_TOPO,NY_TOPO))

  if(OCEANS) then

! for Southern California
    NX_TOPO = NX_TOPO_SOCAL
    NY_TOPO = NY_TOPO_SOCAL
    ORIG_LAT_TOPO = ORIG_LAT_TOPO_SOCAL
    ORIG_LONG_TOPO = ORIG_LONG_TOPO_SOCAL
    DEGREES_PER_CELL_TOPO = DEGREES_PER_CELL_TOPO_SOCAL
    topo_file = TOPO_FILE_SOCAL

    call read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO,topo_file)

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'regional topography file read ranges in m from ',minval(itopo_bathy),' to ',maxval(itopo_bathy)
      write(IMAIN,*)
    endif
  endif

!! read basement map
!  if(BASEMENT_MAP) then
!    call get_value_string(BASEMENT_MAP_FILE,'model.BASEMENT_MAP_FILE','DATA/la_basement/reggridbase2_filtered_ascii.dat')
!    open(unit=55,file=BASEMENT_MAP_FILE,status='old',action='read')
!    do ix=1,NX_BASEMENT
!      do iy=1,NY_BASEMENT
!        read(55,*) iz_basement
!        z_basement(ix,iy) = dble(iz_basement)
!      enddo
!    enddo
!    close(55)
!  endif

  end subroutine gd_read_topography
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_read_partition_files

! reads in proc***_Databases files

  use generate_databases_par
  implicit none

  integer :: num_xmin,num_xmax,num_ymin,num_ymax,num_top,num_bottom,num
  integer :: num_moho
  integer :: j
  character(len=128) :: line
  
! read databases about external mesh simulation
! global node coordinates
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IIN,file=prname(1:len_trim(prname))//'Database',status='old',action='read',form='formatted',iostat=ier)
  if( ier /= 0 ) then
    write(IMAIN,*) 'error opening file: ',prname(1:len_trim(prname))//'Database'
    write(IMAIN,*) 'make sure file exists'
    call exit_mpi(myrank,'error opening database file')
  endif
  read(IIN,*) nnodes_ext_mesh
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh))
  do inode = 1, nnodes_ext_mesh
     read(IIN,*) dummy_node, nodes_coords_ext_mesh(1,inode), nodes_coords_ext_mesh(2,inode), &
                nodes_coords_ext_mesh(3,inode)
  enddo

  call sum_all_i(nnodes_ext_mesh,num)
  if(myrank == 0) then
    write(IMAIN,*) '  external mesh points: ',num
  endif
  call sync_all()

! read materials' physical properties
  read(IIN,*) nmat_ext_mesh, nundefMat_ext_mesh
  allocate(materials_ext_mesh(6,nmat_ext_mesh))
  do imat = 1, nmat_ext_mesh
     ! format:        #(1) rho   #(2) vp  #(3) vs  #(4) Q_flag  #(5) anisotropy_flag  #(6) material_domain_id 
     read(IIN,*) materials_ext_mesh(1,imat),  materials_ext_mesh(2,imat),  materials_ext_mesh(3,imat), &
          materials_ext_mesh(4,imat),  materials_ext_mesh(5,imat), materials_ext_mesh(6,imat)
     
     ! output
     !print*,'materials:',materials_ext_mesh(1,imat),  materials_ext_mesh(2,imat),  materials_ext_mesh(3,imat), &
     !     materials_ext_mesh(4,imat),  materials_ext_mesh(5,imat), materials_ext_mesh(6,imat)
  end do

  if(myrank == 0) then
    write(IMAIN,*) '  defined materials: ',nmat_ext_mesh
  endif
  call sync_all()

  allocate(undef_mat_prop(6,nundefMat_ext_mesh))
  do imat = 1, nundefMat_ext_mesh
     ! format example tomography: 
     ! -1 tomography elastic tomography_model.xyz 1 2              
     ! format example interface: 
     ! -1 interface 14 15 1 2              
     read(IIN,*) undef_mat_prop(1,imat),undef_mat_prop(2,imat),undef_mat_prop(3,imat),undef_mat_prop(4,imat), &
          undef_mat_prop(5,imat), undef_mat_prop(6,imat)

     ! output debug
     !print*,'undefined materials:'
     !print*,undef_mat_prop(:,imat)     
  end do

  if(myrank == 0) then
    write(IMAIN,*) '  undefined materials: ',nundefMat_ext_mesh
  endif
  call sync_all()

! element indexing
  read(IIN,*) nelmnts_ext_mesh
  allocate(elmnts_ext_mesh(esize,nelmnts_ext_mesh))
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh))
  
  ! reads in material association for each spectral element and corner node indices
  do ispec = 1, nelmnts_ext_mesh
     ! format:
     ! # ispec_local # material_index_1 # material_index_2 # corner_id1 # corner_id2 # ... # corner_id8
     read(IIN,*) dummy_elmnt, mat_ext_mesh(1,ispec),mat_ext_mesh(2,ispec), &
          elmnts_ext_mesh(1,ispec), elmnts_ext_mesh(2,ispec), elmnts_ext_mesh(3,ispec), elmnts_ext_mesh(4,ispec), &
          elmnts_ext_mesh(5,ispec), elmnts_ext_mesh(6,ispec), elmnts_ext_mesh(7,ispec), elmnts_ext_mesh(8,ispec)

     ! check debug     
     if( dummy_elmnt /= ispec) stop "error ispec order in materials file"

  enddo
  NSPEC_AB = nelmnts_ext_mesh

  call sum_all_i(nspec_ab,num)
  if(myrank == 0) then
    write(IMAIN,*) '  spectral elements: ',num
  endif
  call sync_all()


! read boundaries
  read(IIN,*) boundary_number ,nspec2D_xmin
  if(boundary_number /= 1) stop "Error : invalid database file"
  read(IIN,*) boundary_number ,nspec2D_xmax
  if(boundary_number /= 2) stop "Error : invalid database file"
  read(IIN,*) boundary_number ,nspec2D_ymin
  if(boundary_number /= 3) stop "Error : invalid database file"
  read(IIN,*) boundary_number ,nspec2D_ymax
  if(boundary_number /= 4) stop "Error : invalid database file"
  read(IIN,*) boundary_number ,nspec2D_bottom_ext
  if(boundary_number /= 5) stop "Error : invalid database file"
  read(IIN,*) boundary_number ,nspec2D_top_ext
  if(boundary_number /= 6) stop "Error : invalid database file"

  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext

  allocate(ibelm_xmin(nspec2D_xmin),nodes_ibelm_xmin(4,nspec2D_xmin))
  do ispec2D = 1,nspec2D_xmin
     read(IIN,*) ibelm_xmin(ispec2D),(nodes_ibelm_xmin(j,ispec2D),j=1,4)
  end do

  allocate(ibelm_xmax(nspec2D_xmax),nodes_ibelm_xmax(4,nspec2D_xmax))
  do ispec2D = 1,nspec2D_xmax
     read(IIN,*) ibelm_xmax(ispec2D),(nodes_ibelm_xmax(j,ispec2D),j=1,4)
  end do

  allocate(ibelm_ymin(nspec2D_ymin),nodes_ibelm_ymin(4,nspec2D_ymin))
  do ispec2D = 1,nspec2D_ymin
     read(IIN,*) ibelm_ymin(ispec2D),(nodes_ibelm_ymin(j,ispec2D),j=1,4)
  end do

  allocate(ibelm_ymax(nspec2D_ymax),nodes_ibelm_ymax(4,nspec2D_ymax))
  do ispec2D = 1,nspec2D_ymax
     read(IIN,*) ibelm_ymax(ispec2D),(nodes_ibelm_ymax(j,ispec2D),j=1,4)
  end do

  allocate(ibelm_bottom(nspec2D_bottom_ext),nodes_ibelm_bottom(4,nspec2D_bottom_ext))
  do ispec2D = 1,nspec2D_bottom_ext
     read(IIN,*) ibelm_bottom(ispec2D),(nodes_ibelm_bottom(j,ispec2D),j=1,4)
  end do

  allocate(ibelm_top(nspec2D_top_ext),nodes_ibelm_top(4,nspec2D_top_ext))
  do ispec2D = 1,nspec2D_top_ext
     read(IIN,*) ibelm_top(ispec2D),(nodes_ibelm_top(j,ispec2D),j=1,4)
  end do

  call sum_all_i(nspec2D_xmin,num_xmin)
  call sum_all_i(nspec2D_xmax,num_xmax)
  call sum_all_i(nspec2D_ymin,num_ymin)
  call sum_all_i(nspec2D_ymax,num_ymax)
  call sum_all_i(nspec2D_top_ext,num_top)
  call sum_all_i(nspec2D_bottom_ext,num_bottom)
  
  if(myrank == 0) then
    write(IMAIN,*) '  absorbing boundaries: '
    write(IMAIN,*) '    xmin,xmax: ',num_xmin,num_xmax
    write(IMAIN,*) '    ymin,ymax: ',num_ymin,num_ymax
    write(IMAIN,*) '    bottom,top: ',num_bottom,num_top
  endif
  call sync_all()

! MPI interfaces between different partitions
  ! format: #number_of_MPI_interfaces  #maximum_number_of_elements_on_each_interface
  read(IIN,*) num_interfaces_ext_mesh, max_interface_size_ext_mesh

  ! allocates interfaces
  allocate(my_neighbours_ext_mesh(num_interfaces_ext_mesh))
  allocate(my_nelmnts_neighbours_ext_mesh(num_interfaces_ext_mesh))
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,num_interfaces_ext_mesh))
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,num_interfaces_ext_mesh))
  allocate(nibool_interfaces_ext_mesh(num_interfaces_ext_mesh))

  ! loops over MPI interfaces with other partitions
  do num_interface = 1, num_interfaces_ext_mesh
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    read(IIN,*) my_neighbours_ext_mesh(num_interface), my_nelmnts_neighbours_ext_mesh(num_interface)
    
    ! loops over interface elements
    do ie = 1, my_nelmnts_neighbours_ext_mesh(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2 #(5)...
      !
      ! interface types: 
      !     1  -  corner point only
      !     2  -  element edge
      !     4  -  element face
      read(IIN,*) my_interfaces_ext_mesh(1,ie,num_interface), my_interfaces_ext_mesh(2,ie,num_interface), &
                  my_interfaces_ext_mesh(3,ie,num_interface), my_interfaces_ext_mesh(4,ie,num_interface), &
                  my_interfaces_ext_mesh(5,ie,num_interface), my_interfaces_ext_mesh(6,ie,num_interface)
    enddo
  enddo
  
  call sum_all_i(num_interfaces_ext_mesh,num)  
  if(myrank == 0) then
    write(IMAIN,*) '  number of MPI partition interfaces: ',num
  endif
  call sync_all()

  ! optional moho
  if( SAVE_MOHO_MESH ) then
    ! checks if additional line exists
    read(IIN,'(a128)',iostat=ier) line 
    if( ier /= 0 ) then 
      ! no moho informations given
      nspec2D_moho_ext = 0
      boundary_number = 7
    else
      ! tries to read in number of moho elements
      read(line,*,iostat=ier) boundary_number ,nspec2D_moho_ext
      if( ier /= 0 ) call exit_mpi(myrank,'error reading moho mesh in database')
    endif    
    if(boundary_number /= 7) stop "Error : invalid database file"

    ! checks total number of elements  
    call sum_all_i(nspec2D_moho_ext,num_moho)
    if( num_moho == 0 ) call exit_mpi(myrank,'error no moho mesh in database')
    
    ! reads in element informations
    allocate(ibelm_moho(nspec2D_moho_ext),nodes_ibelm_moho(4,nspec2D_moho_ext))
    do ispec2D = 1,nspec2D_moho_ext
      ! format: #element_id #node_id1 #node_id2 #node_id3 #node_id4
      read(IIN,*) ibelm_moho(ispec2D),(nodes_ibelm_moho(j,ispec2D),j=1,4)
    end do
  
    ! user output
    if(myrank == 0) then
      write(IMAIN,*) '  moho surfaces: ',num_moho
    endif    
    call sync_all()
  endif
  
  close(IIN)
  
  end subroutine gd_read_partition_files

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_setup_mesh

! mesh creation for static solver

  use generate_databases_par
  implicit none

! assign theoretical number of elements
  nspec = NSPEC_AB

! compute maximum number of points
  npointot = nspec * NGLLCUBE

! use dynamic allocation to allocate memory for arrays
!  allocate(idoubling(nspec))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier) 
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  call memory_eval_mesher(myrank,nspec,npointot,nnodes_ext_mesh,&
                        nelmnts_ext_mesh,nmat_ext_mesh,num_interfaces_ext_mesh, &
                        max_interface_size_ext_mesh,nspec2D_xmin,nspec2D_xmax,&
                        nspec2D_ymin,nspec2D_ymax,nspec2D_bottom,nspec2D_top,&
                        max_static_memory_size_request)
                            
  max_static_memory_size = max_static_memory_size_request    

! make sure everybody is synchronized
  call sync_all()

! main working routine to create all the regions of the mesh
  if(myrank == 0) then
    write(IMAIN,*) 'create regions: '
  endif  
  call create_regions_mesh_ext(ibool, &
                        xstore, ystore, zstore, nspec, npointot, myrank, LOCAL_PATH, &
                        nnodes_ext_mesh, nelmnts_ext_mesh, &
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
                        nodes_ibelm_xmin,nodes_ibelm_xmax,nodes_ibelm_ymin,nodes_ibelm_ymax, &
                        nodes_ibelm_bottom,nodes_ibelm_top, &
                        SAVE_MESH_FILES,nglob, &
                        ANISOTROPY,NPROC,OCEANS, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,NX_TOPO,NY_TOPO, &
                        ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO, &
                        itopo_bathy)

! Moho boundary parameters, 2-D jacobians and normals
  if( SAVE_MOHO_MESH ) then
    call create_regions_mesh_save_moho(myrank,nglob,nspec, &
                        nspec2D_moho_ext,ibelm_moho,nodes_ibelm_moho, &
                        nodes_coords_ext_mesh,nnodes_ext_mesh,ibool )    
  endif

! defines global number of nodes in model
  NGLOB_AB = nglob

! print min and max of topography included
  min_elevation = HUGEVAL
  max_elevation = -HUGEVAL
  do iface = 1,nspec2D_top_ext
     do icorner = 1,NGNOD2D
        inode = nodes_ibelm_top(icorner,iface)
        if (nodes_coords_ext_mesh(3,inode) < min_elevation) then
           min_elevation = nodes_coords_ext_mesh(3,inode)
        end if
        if (nodes_coords_ext_mesh(3,inode) > max_elevation) then
           max_elevation = nodes_coords_ext_mesh(3,inode) 
        end if
     end do
  end do

! compute the maximum of the maxima for all the slices using an MPI reduction
  call min_all_dp(min_elevation,min_elevation_all)
  call max_all_dp(max_elevation,max_elevation_all)
  
  if(myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'min and max of topography included in mesh in m is ',min_elevation_all,' ',max_elevation_all
     write(IMAIN,*)
  endif

! clean-up
  deallocate(xstore,ystore,zstore)

! make sure everybody is synchronized
  call sync_all()

  end subroutine gd_setup_mesh
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_finalize

! checks user input parameters

  use generate_databases_par
  implicit none

  integer :: i
  
! print number of points and elements in the mesh
  call sum_all_i(NGLOB_AB,nglob_total)
  call sum_all_i(NSPEC_AB,nspec_total)
  call sync_all()  
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements:'
    write(IMAIN,*) '-----------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in each slice: ',NSPEC_AB
    write(IMAIN,*) 'total number of points in each slice: ',NGLOB_AB
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in entire mesh: ',nspec_total     ! NSPEC_AB*NPROC
    write(IMAIN,*) 'total number of points in entire mesh: ',nglob_total        !NGLOB_AB*NPROC
    write(IMAIN,*) 'total number of DOFs in entire mesh: ',nglob_total*NDIM   !NGLOB_AB*NPROC*NDIM
    write(IMAIN,*)
    write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if(CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
  endif
  
! gets number of surface elements (for movie outputs)
  allocate( ispec_is_surface_external_mesh(NSPEC_AB), &
           iglob_is_surface_external_mesh(NGLOB_AB),stat=ier)  
  if( ier /= 0 ) stop 'error allocating array'  
  max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh)
  allocate(ibool_interfaces_ext_mesh_dummy(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh),stat=ier)
  if( ier /= 0 ) stop 'error allocating array'  
  do i = 1, num_interfaces_ext_mesh
     ibool_interfaces_ext_mesh_dummy(:,:) = ibool_interfaces_ext_mesh(1:max_nibool_interfaces_ext_mesh,:)
  enddo
  call sync_all()  
  call detect_surface(NPROC,NGLOB_AB,NSPEC_AB,ibool, &
                        ispec_is_surface_external_mesh, &
                        iglob_is_surface_external_mesh, &
                        nfaces_surface_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh, &
                        ibool_interfaces_ext_mesh_dummy )

  deallocate(ibool)
  deallocate(ispec_is_surface_external_mesh)
  deallocate(iglob_is_surface_external_mesh)
  deallocate(ibool_interfaces_ext_mesh_dummy)

  ! takes number of faces for top, free surface only
  if( MOVIE_SURFACE .or. CREATE_SHAKEMAP ) then
    nfaces_surface_ext_mesh = NSPEC2D_TOP
  endif
  
! number of surface faces for all partitions together
  call sum_all_i(nfaces_surface_ext_mesh,nfaces_surface_glob_ext_mesh)

  
! copy number of elements and points in an include file for the solver
  if( myrank == 0 ) then
    call save_header_file(NSPEC_AB,NGLOB_AB,NPROC, &
               ATTENUATION,ANISOTROPY,NSTEP,DT, &
               SIMULATION_TYPE,max_static_memory_size,nfaces_surface_glob_ext_mesh)
  endif 
  
! elapsed time since beginning of mesh generation
  if(myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
  endif

! close main output file
  if(myrank == 0) then
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call sync_all()
  
  end subroutine gd_finalize
