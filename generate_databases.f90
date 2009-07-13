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

  subroutine generate_databases

  implicit none

  include "constants.h"

!=============================================================================!
!                                                                             !
!  generate_databases produces a spectral element grid for a local or regional model.  !
!  The mesher uses the UTM projection                                         !
!                                                                             !
!=============================================================================!
!
! If you use this code for your own research, please cite some of these articles:
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
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! If you use the kernel capabilities of the code, please cite
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
! to Jeroen Tromp <jtromp AT caltech.edu> and/or use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! MPI v. 1.4 Dimitri Komatitsch, University of Pau, Qinya Liu and others, Caltech, September 2006:
!  better adjoint and kernel calculations, faster and better I/Os
!  on very large systems, many small improvements and bug fixes
! MPI v. 1.3 Dimitri Komatitsch, University of Pau, and Qinya Liu, Caltech, July 2005:
!  serial version, regular mesh, adjoint and kernel calculations, ParaView support
! MPI v. 1.2 Min Chen and Dimitri Komatitsch, Caltech, July 2004:
!  full anisotropy, volume movie
! MPI v. 1.1 Dimitri Komatitsch, Caltech, October 2002: Zhu's Moho map, scaling
!  of Vs with depth, Hauksson's regional model, attenuation, oceans, movies
! MPI v. 1.0 Dimitri Komatitsch, Caltech, May 2002: first MPI version
!                        based on global code

! number of spectral elements in each block
  integer nspec,npointot

! auxiliary variables to generate the mesh
  integer ix,iy

! parameters needed to store the radii of the grid points
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier

! check area and volume of the final mesh
  double precision area_local_bottom,area_total_bottom
  double precision area_local_top,area_total_top
  double precision volume_local,volume_total

! use integer array to store topography values
  integer NX_TOPO,NY_TOPO
  double precision ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  character(len=100) topo_file
  integer, dimension(:,:), allocatable :: itopo_bathy

! use integer array to store Moho depth
  integer imoho_depth(NX_MOHO,NY_MOHO)

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_XI,NEX_ETA, &
             NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX
  double precision Z_DEPTH_BLOCK
  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

! static memory size that will be needed by the solver
  double precision :: static_memory_size

! this for all the regions
  integer NSPEC_AB,NGLOB_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

  double precision min_elevation,max_elevation
  double precision min_elevation_all,max_elevation_all

! for tapered basement map
  integer iz_basement
  double precision z_basement(NX_BASEMENT,NY_BASEMENT)
  character(len=150) BASEMENT_MAP_FILE

! to filter list of stations
!   integer nrec,nrec_filtered
!   double precision stlat,stlon,stele,stbur
!   character(len=MAX_LENGTH_STATION_NAME) station_name
!   character(len=MAX_LENGTH_NETWORK_NAME) network_name
!   character(len=150) rec_filename!,filtered_rec_filename

! for Databases of external meshes
  character(len=150) prname
  integer :: dummy_node
  integer :: dummy_elmnt
  integer :: ispec, inode, num_interface, ie,imat !pll
  integer :: nnodes_ext_mesh, nelmnts_ext_mesh
  integer  :: ninterface_ext_mesh
  integer  :: max_interface_size_ext_mesh
  integer  :: nmat_ext_mesh, nundefMat_ext_mesh   !pll
  integer, dimension(:), allocatable  :: my_neighbours_ext_mesh
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours_ext_mesh
  integer, dimension(:,:,:), allocatable  :: my_interfaces_ext_mesh
  integer, dimension(:,:), allocatable  :: ibool_interfaces_ext_mesh
  integer, dimension(:), allocatable  :: nibool_interfaces_ext_mesh
  double precision, dimension(:,:), allocatable :: nodes_coords_ext_mesh
  integer, dimension(:,:), allocatable :: elmnts_ext_mesh
  integer, dimension(:,:), allocatable :: mat_ext_mesh

  ! pll
  double precision, dimension(:,:), allocatable :: materials_ext_mesh
  integer, dimension(:), allocatable  :: ibelm_xmin,ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top  
  integer  :: ispec2D, boundary_number
  integer  :: nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, nspec2D_bottom_ext, nspec2D_top_ext
  character (len=30), dimension(:,:), allocatable :: undef_mat_prop

! ************** PROGRAM STARTS HERE **************

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
  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT,NER_MOHO_16,NER_BOTTOM_MOHO, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,DT, &
        ATTENUATION,USE_OLSEN_ATTENUATION,HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
        OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL,ANISOTROPY, &
        BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS, &
        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

  if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) then
    stop 'must have NPROC_XI = NPROC_ETA = 1 for a serial run'
  endif

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,USE_REGULAR_MESH)

! info about external mesh simulation
! nlegoff -- should be put in compute_parameters and read_parameter_file for clarity
  NPROC = sizeprocs

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) call exit_MPI(myrank,'wrong number of MPI processes')

  if(myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta'
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
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

  if(NGNOD /= 8) call exit_MPI(myrank,'number of control nodes must be 8')
  if(NGNOD2D /= 4) call exit_MPI(myrank,'elements with 8 points should have NGNOD2D = 4')

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) call exit_MPI(myrank,'number of SLS must be 3')

  if(myrank == 0) then

  write(IMAIN,*) 'region selected:'
  write(IMAIN,*)
  write(IMAIN,*) 'latitude min = ',LATITUDE_MIN
  write(IMAIN,*) 'latitude max = ',LATITUDE_MAX
  write(IMAIN,*)
  write(IMAIN,*) 'longitude min = ',LONGITUDE_MIN
  write(IMAIN,*) 'longitude max = ',LONGITUDE_MAX
  write(IMAIN,*)
  write(IMAIN,*) 'this is mapped to UTM in region ',UTM_PROJECTION_ZONE
  write(IMAIN,*)
  write(IMAIN,*) 'UTM X min = ',UTM_X_MIN
  write(IMAIN,*) 'UTM X max = ',UTM_X_MAX
  write(IMAIN,*)
  write(IMAIN,*) 'UTM Y min = ',UTM_Y_MIN
  write(IMAIN,*) 'UTM Y max = ',UTM_Y_MAX
  write(IMAIN,*)
  write(IMAIN,*) 'UTM size of model along X is ',(UTM_X_MAX-UTM_X_MIN)/1000.,' km'
  write(IMAIN,*) 'UTM size of model along Y is ',(UTM_Y_MAX-UTM_Y_MIN)/1000.,' km'
  write(IMAIN,*)
  write(IMAIN,*) 'Bottom of the mesh is at a depth of ',dabs(Z_DEPTH_BLOCK)/1000.,' km'
  write(IMAIN,*)


  write(IMAIN,*)
  if(TOPOGRAPHY) then
    write(IMAIN,*) 'incorporating surface topography'
  else
    write(IMAIN,*) 'no surface topography'
  endif

  write(IMAIN,*)
  if(SUPPRESS_UTM_PROJECTION) then
    write(IMAIN,*) 'suppressing UTM projection'
  else
    write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
  endif

  write(IMAIN,*)
  if(HARVARD_3D_GOCAD_MODEL) then
    write(IMAIN,*) 'incorporating 3-D lateral variations'
  else
    write(IMAIN,*) 'no 3-D lateral variations'
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
  if(OCEANS) then
    write(IMAIN,*) 'incorporating the oceans using equivalent load'
  else
    write(IMAIN,*) 'no oceans'
  endif

  write(IMAIN,*)

  endif

! read topography and bathymetry file
  if(TOPOGRAPHY .or. OCEANS) then

! for Southern California
    NX_TOPO = NX_TOPO_SOCAL
    NY_TOPO = NY_TOPO_SOCAL
    ORIG_LAT_TOPO = ORIG_LAT_TOPO_SOCAL
    ORIG_LONG_TOPO = ORIG_LONG_TOPO_SOCAL
    DEGREES_PER_CELL_TOPO = DEGREES_PER_CELL_TOPO_SOCAL
    topo_file = TOPO_FILE_SOCAL

    allocate(itopo_bathy(NX_TOPO,NY_TOPO))

    call read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO,topo_file)

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'regional topography file read ranges in m from ',minval(itopo_bathy),' to ',maxval(itopo_bathy)
      write(IMAIN,*)
    endif

  endif

! read Moho map
  if(MOHO_MAP_LUPEI) then
    call read_moho_map(imoho_depth)
    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'regional Moho depth read ranges in m from ',minval(imoho_depth),' to ',maxval(imoho_depth)
      write(IMAIN,*)
    endif
  endif

! read basement map
  if(BASEMENT_MAP) then
    call get_value_string(BASEMENT_MAP_FILE,'model.BASEMENT_MAP_FILE','DATA/la_basement/reggridbase2_filtered_ascii.dat')
    open(unit=55,file=BASEMENT_MAP_FILE,status='old',action='read')
    do ix=1,NX_BASEMENT
      do iy=1,NY_BASEMENT
        read(55,*) iz_basement
        z_basement(ix,iy) = dble(iz_basement)
      enddo
    enddo
    close(55)
  endif

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! volume of bottom and top area of the slice
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

! read databases about external mesh simulation
  
  call create_name_database(prname,myrank,LOCAL_PATH)
  open(unit=IIN,file=prname(1:len_trim(prname))//'Database',status='old',action='read',form='formatted')
  read(IIN,*) nnodes_ext_mesh
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh))
  do inode = 1, nnodes_ext_mesh
     read(IIN,*) dummy_node, nodes_coords_ext_mesh(1,inode), nodes_coords_ext_mesh(2,inode), nodes_coords_ext_mesh(3,inode)
  enddo
  
! read materials' physical properties    
  read(IIN,*) nmat_ext_mesh, nundefMat_ext_mesh 
  allocate(materials_ext_mesh(5,nmat_ext_mesh))
  allocate(undef_mat_prop(5,nundefMat_ext_mesh))
  do imat = 1, nmat_ext_mesh
     read(IIN,*) materials_ext_mesh(1,imat),  materials_ext_mesh(2,imat),  materials_ext_mesh(3,imat), &
          materials_ext_mesh(4,imat),  materials_ext_mesh(5,imat)
  end do
  
  do imat = 1, nundefMat_ext_mesh
     read(IIN,*) undef_mat_prop(1,imat),undef_mat_prop(2,imat),undef_mat_prop(3,imat),undef_mat_prop(4,imat), &
          undef_mat_prop(5,imat)
  end do

  read(IIN,*) nelmnts_ext_mesh
  allocate(elmnts_ext_mesh(esize,nelmnts_ext_mesh))
  allocate(mat_ext_mesh(2,nelmnts_ext_mesh))
  do ispec = 1, nelmnts_ext_mesh
     read(IIN,*) dummy_elmnt, mat_ext_mesh(1,ispec),mat_ext_mesh(2,ispec), &
          elmnts_ext_mesh(1,ispec), elmnts_ext_mesh(2,ispec), elmnts_ext_mesh(3,ispec), elmnts_ext_mesh(4,ispec), &
          elmnts_ext_mesh(5,ispec), elmnts_ext_mesh(6,ispec), elmnts_ext_mesh(7,ispec), elmnts_ext_mesh(8,ispec)
  enddo
  NSPEC_AB = nelmnts_ext_mesh
  
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
  NSPEC2DMAX_XMIN_XMAX = max(nspec2D_xmin,nspec2D_xmax)
  NSPEC2DMAX_YMIN_YMAX = max(nspec2D_ymin,nspec2D_ymax)
  NSPEC2D_BOTTOM = nspec2D_bottom_ext
  NSPEC2D_TOP = nspec2D_top_ext  
  
  allocate(ibelm_xmin(nspec2D_xmin))
  do ispec2D = 1,nspec2D_xmin
     read(IIN,*) ibelm_xmin(ispec2D)
  end do
  
  allocate(ibelm_xmax(nspec2D_xmax))
  do ispec2D = 1,nspec2D_xmax
     read(IIN,*) ibelm_xmax(ispec2D)
  end do
  
  allocate(ibelm_ymin(nspec2D_ymin))
  do ispec2D = 1,nspec2D_ymin
     read(IIN,*) ibelm_ymin(ispec2D)
  end do
  
  allocate(ibelm_ymax(nspec2D_ymax))
  do ispec2D = 1,nspec2D_ymax
     read(IIN,*) ibelm_ymax(ispec2D)
  end do
  
  allocate(ibelm_bottom(nspec2D_bottom_ext))
  do ispec2D = 1,nspec2D_bottom_ext
     read(IIN,*) ibelm_bottom(ispec2D)
  end do

  allocate(ibelm_top(nspec2D_top_ext))
  do ispec2D = 1,nspec2D_top_ext
     read(IIN,*) ibelm_top(ispec2D)
  end do
  
  read(IIN,*) ninterface_ext_mesh, max_interface_size_ext_mesh
  allocate(my_neighbours_ext_mesh(ninterface_ext_mesh))
  allocate(my_nelmnts_neighbours_ext_mesh(ninterface_ext_mesh))
  allocate(my_interfaces_ext_mesh(6,max_interface_size_ext_mesh,ninterface_ext_mesh))
  allocate(ibool_interfaces_ext_mesh(NGLLX*NGLLX*max_interface_size_ext_mesh,ninterface_ext_mesh))
  allocate(nibool_interfaces_ext_mesh(ninterface_ext_mesh))
  do num_interface = 1, ninterface_ext_mesh
     read(IIN,*) my_neighbours_ext_mesh(num_interface), my_nelmnts_neighbours_ext_mesh(num_interface)
     do ie = 1, my_nelmnts_neighbours_ext_mesh(num_interface)
        read(IIN,*) my_interfaces_ext_mesh(1,ie,num_interface), my_interfaces_ext_mesh(2,ie,num_interface), &
             my_interfaces_ext_mesh(3,ie,num_interface), my_interfaces_ext_mesh(4,ie,num_interface), &
             my_interfaces_ext_mesh(5,ie,num_interface), my_interfaces_ext_mesh(6,ie,num_interface)
     enddo
  enddo
  
  close(IIN)


! assign theoretical number of elements
  nspec = NSPEC_AB

! compute maximum number of points
  npointot = nspec * NGLLCUBE

! make sure everybody is synchronized
  call sync_all()

! use dynamic allocation to allocate memory for arrays
  allocate(idoubling(nspec))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)

! exit if there is not enough memory to allocate all the arrays
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create all the regions of the mesh
  call create_regions_mesh_ext_mesh(ibool, &
       xstore,ystore,zstore,nspec,npointot,myrank,LOCAL_PATH, &
       nnodes_ext_mesh,nelmnts_ext_mesh, &
       nodes_coords_ext_mesh,elmnts_ext_mesh,static_memory_size,mat_ext_mesh,materials_ext_mesh, &
       nmat_ext_mesh,undef_mat_prop,nundefMat_ext_mesh,ninterface_ext_mesh,max_interface_size_ext_mesh, &
       my_neighbours_ext_mesh,my_nelmnts_neighbours_ext_mesh,my_interfaces_ext_mesh, &
       ibool_interfaces_ext_mesh,nibool_interfaces_ext_mesh, &
       nspec2D_xmin, nspec2D_xmax, nspec2D_ymin, nspec2D_ymax, NSPEC2D_BOTTOM, NSPEC2D_TOP,&
       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
       ibelm_xmin,ibelm_xmax, ibelm_ymin, ibelm_ymax, ibelm_bottom, ibelm_top)
  

! print min and max of topography included
  if(TOPOGRAPHY) then

! compute the maximum of the maxima for all the slices using an MPI reduction
      call min_all_dp(min_elevation,min_elevation_all)
      call max_all_dp(max_elevation,max_elevation_all)

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'min and max of topography included in mesh in m is ',min_elevation_all,' ',max_elevation_all
      write(IMAIN,*)
    endif
  endif


! use MPI reduction to compute total area and volume
  area_total_bottom   = ZERO
  area_total_top   = ZERO
  call sum_all_dp(area_local_bottom,area_total_bottom)
  call sum_all_dp(area_local_top,area_total_top)
  call sum_all_dp(volume_local,volume_total)

  if(myrank == 0) then

!   check volume, and bottom and top area

      write(IMAIN,*)
      write(IMAIN,*) '   calculated top area: ',area_total_top

! compare to exact theoretical value
    if(.not. TOPOGRAPHY) &
          write(IMAIN,*) '            exact area: ',(UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)

      write(IMAIN,*)
      write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

! compare to exact theoretical value (bottom is always flat)
      write(IMAIN,*) '            exact area: ',(UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)

  endif

! make sure everybody is synchronized
  call sync_all()

  if(myrank == 0) then
! check volume
      write(IMAIN,*)
      write(IMAIN,*) 'calculated volume: ',volume_total
! take the central cube into account
   if(.not. TOPOGRAPHY) &
      write(IMAIN,*) '     exact volume: ', &
        (UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)*dabs(Z_DEPTH_BLOCK)

  endif

!--- print number of points and elements in the mesh

  if(myrank == 0) then

  write(IMAIN,*)
  write(IMAIN,*) 'Repartition of elements:'
  write(IMAIN,*) '-----------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements in each slice: ',NSPEC_AB
  write(IMAIN,*)
  write(IMAIN,*) 'total number of points in each slice: ',NGLOB_AB

  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements in entire mesh: ',NSPEC_AB*NPROC
  write(IMAIN,*) 'total number of points in entire mesh: ',NGLOB_AB*NPROC
  write(IMAIN,*) 'total number of DOFs in entire mesh: ',NGLOB_AB*NPROC*NDIM
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

! copy number of elements and points in an include file for the solver
  call save_header_file(NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC, &
             UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,ATTENUATION,ANISOTROPY,NSTEP, &
             NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,SIMULATION_TYPE,static_memory_size)

!  call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS')
!  call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_FILTERED')

! get total number of stations
! open(unit=IIN,file=rec_filename,iostat=ios,status='old',action='read')
! nrec = 0
! do while(ios == 0)
!   read(IIN,"(a)",iostat=ios) dummystring
!   if(ios == 0) nrec = nrec + 1
! enddo
! close(IIN)

! filter list of stations, only retain stations that are in the model
!  nrec_filtered = 0
!  open(unit=IIN,file=rec_filename,status='old',action='read')
!  do irec = 1,nrec
!    read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
!    if((stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
!         .or. USE_EXTERNAL_MESH) &
!      nrec_filtered = nrec_filtered + 1
!  enddo
!  close(IIN)

!  write(IMAIN,*)
!  write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(rec_filename)
!  write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_rec_filename)
!  write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
!  write(IMAIN,*)

!  if(nrec_filtered < 1) call exit_MPI(myrank,'need at least one station in the model')

!  if(nrec < 1) call exit_MPI(myrank,'need at least one station in the model')

!  open(unit=IIN,file=rec_filename,status='old',action='read')
!  open(unit=IOUT,file=filtered_rec_filename,status='unknown')

!  do irec = 1,nrec
!    read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
!    if((stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
!         .or. USE_EXTERNAL_MESH) &
!      write(IOUT,*) station_name(1:len_trim(station_name)),' ',network_name(1:len_trim(network_name)),' ', &
!              sngl(stlat),' ',sngl(stlon), ' ', sngl(stele), ' ', sngl(stbur)
!  enddo

!  close(IIN)
!  close(IOUT)

  endif   ! end of section executed by main process only

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

  end subroutine generate_databases

