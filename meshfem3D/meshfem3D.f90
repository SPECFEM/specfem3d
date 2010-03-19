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

  subroutine meshfem3D

  use readParFile
  use createRegMesh
  implicit none
  include "constants.h"
 
!=============================================================================!
!                                                                             !
!  meshfem3D produces a spectral element grid for a local or regional model.  !
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

! number of spectral elements in each block
  integer nspec,npointot

! meshing parameters
  double precision, dimension(:), allocatable :: rns

! auxiliary variables to generate the mesh
  integer ix,iy,ir

  double precision xin,etan,rn
  double precision x_current,y_current,z_top,z_bot

  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid

  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier

! check area and volume of the final mesh
!   double precision area_local_bottom,area_total_bottom
!   double precision area_local_top,area_total_top
!   double precision volume_local,volume_total

  integer iprocnum,npx,npy

! for loop on all the slices
  integer iproc_xi,iproc_eta
  integer, dimension(:,:), allocatable :: addressing

! use integer array to store topography values
  integer icornerlat,icornerlong,NX_TOPO,NY_TOPO
  double precision lat,long,elevation,ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO
  double precision long_corner,lat_corner,ratio_xi,ratio_eta
  character(len=100) topo_file
  integer, dimension(:,:), allocatable :: itopo_bathy

! use integer array to store Moho depth
  integer imoho_depth(NX_MOHO,NY_MOHO)

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

! addressing for all the slices
  integer, dimension(:), allocatable :: iproc_xi_slice,iproc_eta_slice

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_XI,NEX_ETA, &
             NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX
  double precision Z_DEPTH_BLOCK,Z_BASEMENT_SURFACE,Z_DEPTH_MOHO
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX!,HDUR_MOVIE

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,&!ATTENUATION,USE_OLSEN_ATTENUATION, &
          IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI!,ABSORBING_CONDITIONS,SAVE_FORWARD
!  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION

  logical SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
!  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

! this for all the regions
  integer NSPEC_AB,NGLOB_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX

  double precision min_elevation,max_elevation
  double precision min_elevation_all,max_elevation_all

! for tapered basement map
  integer icorner_x,icorner_y
  integer iz_basement
  double precision x_corner,y_corner
  double precision z_basement(NX_BASEMENT,NY_BASEMENT)
  character(len=150) BASEMENT_MAP_FILE

! interfaces parameters
  logical SUPPRESS_UTM_PROJECTION_BOTTOM,SUPPRESS_UTM_PROJECTION_TOP
  integer ipoint_current,ilayer,interface_current
  integer number_of_interfaces,number_of_layers
  integer max_npx_interface,max_npy_interface
  integer npx_interface_bottom,npy_interface_bottom
  integer npx_interface_top,npy_interface_top
  double precision z_interface_bottom,z_interface_top
  double precision orig_x_interface_bottom,orig_y_interface_bottom
  double precision orig_x_interface_top,orig_y_interface_top
  double precision spacing_x_interface_bottom,spacing_y_interface_bottom
  double precision spacing_x_interface_top,spacing_y_interface_top
  character(len=50) INTERFACES_FILE
  integer, dimension(:), allocatable :: ner_layer
  double precision, dimension(:,:),allocatable :: interface_bottom,interface_top

! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma

! subregions parameters
  integer NSUBREGIONS
!  definition of the different regions of the model in the mesh (nx,ny,nz) 
!  #1 #2 : nx_begining,nx_end  
!  #3 #4 : ny_begining,ny_end 
!  #5 #6 : nz_begining,nz_end   
!     #7 : material number 
  integer, dimension(:,:), pointer :: subregions

! material properties
  integer NMATERIALS
! first dimension  : material_id
! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision , dimension(:,:), pointer :: material_properties

! to filter list of stations
!   integer irec,nrec,nrec_filtered,ios
!   double precision stlat,stlon,stele,stbur
!   character(len=MAX_LENGTH_STATION_NAME) station_name
!   character(len=MAX_LENGTH_NETWORK_NAME) network_name
!   character(len=150) rec_filename,filtered_rec_filename,dummystring

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

 ! nullify(subregions,material_properties) 
  call read_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE, &
        LOCAL_PATH,SUPPRESS_UTM_PROJECTION,&
        INTERFACES_FILE,NSUBREGIONS,subregions,NMATERIALS,material_properties, &
        USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) then
    stop 'must have NPROC_XI = NPROC_ETA = 1 for a serial run'
  endif

! get interface data from external file to count the spectral elements along Z
  if(myrank == 0) then
     write(IMAIN,*) 'Reading interface data from file DATA/',INTERFACES_FILE(1:len_trim(INTERFACES_FILE)), &
          ' to count the spectral elements'
  end if

  open(unit=IIN,file='DATA/'//INTERFACES_FILE,status='old')

  max_npx_interface  = -1
  max_npy_interface  = -1

! read number of interfaces
  call read_value_integer(DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES')
  if(number_of_interfaces < 1) stop 'not enough interfaces (minimum is 1, for topography)'

! loop on all the interfaces
  do interface_current = 1,number_of_interfaces
     
     call read_interface_parameters(SUPPRESS_UTM_PROJECTION_BOTTOM,npx_interface_bottom,npy_interface_bottom,&
             orig_x_interface_bottom,orig_y_interface_bottom,spacing_x_interface_bottom,spacing_y_interface_bottom)
     max_npx_interface = max(npx_interface_bottom,max_npx_interface)
     max_npy_interface = max(npy_interface_bottom,max_npy_interface)

     if((max_npx_interface < 2) .or.(max_npy_interface < 2)) stop 'not enough interface points (minimum is 2x2)'
     
! skip all the points describing this interface
     do ipoint_current = 1,npx_interface_bottom*npy_interface_bottom
        read(IIN,*)
     enddo
     
  enddo

  ! define number of layers
  number_of_layers = number_of_interfaces! - 1
  allocate(ner_layer(number_of_layers))

! loop on all the layers
  do ilayer = 1,number_of_layers

! read number of spectral elements in vertical direction in this layer
     call read_value_integer(DONT_IGNORE_JUNK,ner_layer(ilayer),'NER_LAYER')
     if(ner_layer(ilayer) < 1) stop 'not enough spectral elements along Z in layer (minimum is 1)'

  enddo

  close(IIN)

! compute total number of spectral elements in vertical direction
  NER = sum(ner_layer)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,&
      USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) call exit_MPI(myrank,'wrong number of MPI processes')

! dynamic allocation of mesh arrays
  allocate(rns(0:2*NER))

  allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA))
  allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA))
  allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA))

  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(iproc_xi_slice(0:NPROC-1))
  allocate(iproc_eta_slice(0:NPROC-1))

! clear arrays
  xgrid(:,:,:) = 0.
  ygrid(:,:,:) = 0.
  zgrid(:,:,:) = 0.

  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

! create global slice addressing for solver
  if(myrank == 0) then
    write(IMAIN,*) 'creating global slice addressing'
    write(IMAIN,*)
 endif
    do iproc_eta=0,NPROC_ETA-1
      do iproc_xi=0,NPROC_XI-1
        iprocnum = iproc_eta * NPROC_XI + iproc_xi
        iproc_xi_slice(iprocnum) = iproc_xi
        iproc_eta_slice(iprocnum) = iproc_eta
        addressing(iproc_xi,iproc_eta) = iprocnum
      enddo
    enddo

  if (myrank == 0) then
    write(IMAIN,*) 'Spatial distribution of slice numbers:'
    do iproc_eta = NPROC_ETA-1, 0, -1
      do iproc_xi = 0, NPROC_XI-1, 1
        write(IMAIN,'(i5)',advance='no') addressing(iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(a1)',advance='yes') ' '
    enddo
  endif

  if(myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta'
    write(IMAIN,*) 'There are ',NER,' elements along Z'
    do ilayer = 1,number_of_layers
       write(IMAIN,*) 'There are ',ner_layer(ilayer),' spectral elements along Z in layer ',ilayer
    end do
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

! check that number of slices is at least 1 in each direction
  if(NPROC_XI < 1) call exit_MPI(myrank,'NPROC_XI must be greater than 1')
  if(NPROC_ETA < 1) call exit_MPI(myrank,'NPROC_ETA must be greater than 1')

! check that mesh can be cut into the right number of slices
! also check that mesh can be coarsened in depth twice (block size multiple of 8)
  if(USE_REGULAR_MESH) then
    if(mod(NEX_XI,NPROC_XI) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of NPROC_XI for a regular mesh')
    if(mod(NEX_ETA,NPROC_ETA) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of NPROC_ETA for a regular mesh')
  else
    !if(mod(NEX_XI,8) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8 for a non-regular mesh')
    !if(mod(NEX_ETA,8) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8 for a non-regular mesh')

    !if(mod(NEX_XI/8,NPROC_XI) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8*NPROC_XI for a non-regular mesh')
    !if(mod(NEX_ETA/8,NPROC_ETA) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8*NPROC_ETA for a non-regular mesh')
  endif

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
!   if(TOPOGRAPHY) then
!     write(IMAIN,*) 'incorporating surface topography'
!   else
!     write(IMAIN,*) 'no surface topography'
!   endif

  write(IMAIN,*)
  if(SUPPRESS_UTM_PROJECTION) then
    write(IMAIN,*) 'suppressing UTM projection'
  else
    write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
  endif

  write(IMAIN,*)
!   if(HARVARD_3D_GOCAD_MODEL) then
!     write(IMAIN,*) 'incorporating 3-D lateral variations'
!   else
!     write(IMAIN,*) 'no 3-D lateral variations'
!   endif

!   write(IMAIN,*)
!   if(ATTENUATION) then
!     write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
!     if(USE_OLSEN_ATTENUATION) then
!       write(IMAIN,*) 'using Olsen''s attenuation'
!     else
!       write(IMAIN,*) 'not using Olsen''s attenuation'
!     endif
!   else
!     write(IMAIN,*) 'no attenuation'
!   endif

!   write(IMAIN,*)
!   if(OCEANS) then
!     write(IMAIN,*) 'incorporating the oceans using equivalent load'
!   else
!     write(IMAIN,*) 'no oceans'
!   endif

  write(IMAIN,*)

  endif

! read topography and bathymetry file
!  if(TOPOGRAPHY .or. OCEANS) then
!   if(TOPOGRAPHY) then

! ! for Southern California
!     NX_TOPO = NX_TOPO_SOCAL
!     NY_TOPO = NY_TOPO_SOCAL
!     ORIG_LAT_TOPO = ORIG_LAT_TOPO_SOCAL
!     ORIG_LONG_TOPO = ORIG_LONG_TOPO_SOCAL
!     DEGREES_PER_CELL_TOPO = DEGREES_PER_CELL_TOPO_SOCAL
!     topo_file = TOPO_FILE_SOCAL

!     allocate(itopo_bathy(NX_TOPO,NY_TOPO))

!     call read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO,topo_file)

!     if(myrank == 0) then
!       write(IMAIN,*)
!       write(IMAIN,*) 'regional topography file read ranges in m from ',minval(itopo_bathy),' to ',maxval(itopo_bathy)
!       write(IMAIN,*)
!     endif

!   endif

! ! read Moho map
!   if(MOHO_MAP_LUPEI) then
!     call read_moho_map(imoho_depth)
!     if(myrank == 0) then
!       write(IMAIN,*)
!       write(IMAIN,*) 'regional Moho depth read ranges in m from ',minval(imoho_depth),' to ',maxval(imoho_depth)
!       write(IMAIN,*)
!     endif
!   endif

! ! read basement map
!   if(BASEMENT_MAP) then
!     call get_value_string(BASEMENT_MAP_FILE,'model.BASEMENT_MAP_FILE','DATA/la_basement/reggridbase2_filtered_ascii.dat')
!     open(unit=55,file=BASEMENT_MAP_FILE,status='old',action='read')
!     do ix=1,NX_BASEMENT
!       do iy=1,NY_BASEMENT
!         read(55,*) iz_basement
!         z_basement(ix,iy) = dble(iz_basement)
!       enddo
!     enddo
!     close(55)
!   endif

! get addressing for this process
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

! number of elements in each slice
  npx = 2*NEX_PER_PROC_XI
  npy = 2*NEX_PER_PROC_ETA
  ner_layer(:) = 2 * ner_layer(:)
  min_elevation = +HUGEVAL
  max_elevation = -HUGEVAL

  if(myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Reading interface data from file DATA/',INTERFACES_FILE(1:len_trim(INTERFACES_FILE))
     write(IMAIN,*)
  end if

     open(unit=IIN,file='DATA/'//INTERFACES_FILE,status='old')

     allocate(interface_bottom(max_npx_interface,max_npy_interface))
     allocate(interface_top(max_npx_interface,max_npy_interface))

     ! read number of interfaces
     call read_value_integer(DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES')

!      ! read bottom interface
!      call read_interface_parameters(SUPPRESS_UTM_PROJECTION_BOTTOM,npx_interface_bottom,npy_interface_bottom,&
!           orig_x_interface_bottom,orig_y_interface_bottom,spacing_x_interface_bottom,spacing_y_interface_bottom) 

!      ! loop on all the points describing this interface
!      do ix=1,npx_interface_bottom
!         do iy=1,npy_interface_bottom
!            call read_value_double_precision(DONT_IGNORE_JUNK,interface_bottom(ix,iy),'Z_INTERFACE_BOTTOM') 
!         enddo
!      enddo
     SUPPRESS_UTM_PROJECTION_BOTTOM = SUPPRESS_UTM_PROJECTION
     npx_interface_bottom = 2
     npy_interface_bottom = 2
     orig_x_interface_bottom = UTM_X_MIN
     orig_y_interface_bottom = UTM_Y_MIN
     spacing_x_interface_bottom = UTM_X_MAX - UTM_X_MIN
     spacing_y_interface_bottom = UTM_Y_MAX - UTM_Y_MIN
     interface_bottom(:,:) = - dabs(Z_DEPTH_BLOCK)

     ! loop on all the layers
     do ilayer = 1,number_of_layers

        ! read top interface
        call read_interface_parameters(SUPPRESS_UTM_PROJECTION_TOP,npx_interface_top,npy_interface_top,&
             orig_x_interface_top,orig_y_interface_top,spacing_x_interface_top,spacing_y_interface_top) 

        !npoints_interface_top = npx_interface_top * npy_interface
        ! loop on all the points describing this interface
        do ix=1,npx_interface_top
           do iy=1,npy_interface_top
              call read_value_double_precision(DONT_IGNORE_JUNK,interface_top(ix,iy),'Z_INTERFACE_TOP') 
           enddo
        enddo        

        ! compute the offset of this layer in terms of number of spectral elements below along Z
        if(ilayer > 1) then
           ioffset = sum(ner_layer(1:ilayer-1))
        else
           ioffset = 0
        endif

        !--- definition of the mesh

        do iy=0,npy
           do ix=0,npx

!   define the mesh points on the top and the bottom
              xin=dble(ix)/dble(npx)
              x_current = UTM_X_MIN + (dble(iproc_xi)+xin)*(UTM_X_MAX-UTM_X_MIN)/dble(NPROC_XI)

              etan=dble(iy)/dble(npy)
              y_current = UTM_Y_MIN + (dble(iproc_eta)+etan)*(UTM_Y_MAX-UTM_Y_MIN)/dble(NPROC_ETA)

! get bottom interface value 
! project x and y in UTM back to long/lat since topo file is in long/lat
              call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION_BOTTOM)
                 
! get coordinate of corner in bathy/topo model
              icornerlong = int((long - orig_x_interface_bottom) / spacing_x_interface_bottom) + 1
              icornerlat = int((lat - orig_y_interface_bottom) / spacing_x_interface_bottom) + 1

! avoid edge effects and extend with identical point if outside model
              if(icornerlong < 1) icornerlong = 1
              if(icornerlong > npx_interface_bottom-1) icornerlong = npx_interface_bottom-1
              if(icornerlat < 1) icornerlat = 1
              if(icornerlat > npy_interface_bottom-1) icornerlat = npy_interface_bottom-1

! compute coordinates of corner
              long_corner = orig_x_interface_bottom + (icornerlong-1)*spacing_x_interface_bottom
              lat_corner = orig_y_interface_bottom + (icornerlat-1)*spacing_y_interface_bottom

! compute ratio for interpolation
              ratio_xi = (long - long_corner) / spacing_x_interface_bottom
              ratio_eta = (lat - lat_corner) / spacing_y_interface_bottom

! avoid edge effects
              if(ratio_xi < 0.) ratio_xi = 0.
              if(ratio_xi > 1.) ratio_xi = 1.
              if(ratio_eta < 0.) ratio_eta = 0.
              if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
              z_interface_bottom = &
                   interface_bottom(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
                   interface_bottom(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
                   interface_bottom(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
                   interface_bottom(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

! get top interface value
! project x and y in UTM back to long/lat since topo file is in long/lat
              call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION_TOP)

! get coordinate of corner in bathy/topo model
              icornerlong = int((long - orig_x_interface_top) / spacing_x_interface_top) + 1
              icornerlat = int((lat - orig_y_interface_top) / spacing_x_interface_top) + 1

! avoid edge effects and extend with identical point if outside model
              if(icornerlong < 1) icornerlong = 1
              if(icornerlong > npx_interface_top-1) icornerlong = npx_interface_top-1
              if(icornerlat < 1) icornerlat = 1
              if(icornerlat > npy_interface_top-1) icornerlat = npy_interface_top-1

! compute coordinates of corner
              long_corner = orig_x_interface_top + (icornerlong-1)*spacing_x_interface_top
              lat_corner = orig_y_interface_top + (icornerlat-1)*spacing_y_interface_top

! compute ratio for interpolation
              ratio_xi = (long - long_corner) / spacing_x_interface_top
              ratio_eta = (lat - lat_corner) / spacing_y_interface_top

! avoid edge effects
              if(ratio_xi < 0.) ratio_xi = 0.
              if(ratio_xi > 1.) ratio_xi = 1.
              if(ratio_eta < 0.) ratio_eta = 0.
              if(ratio_eta > 1.) ratio_eta = 1.

! interpolate elevation at current point
              z_interface_top = &
                   interface_top(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
                   interface_top(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
                   interface_top(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
                   interface_top(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

              do ir = 0,ner_layer(ilayer)
              ! linear interpolation between bottom and top
                 gamma = dble(ir) / dble(ner_layer(ilayer))

              ! coordinates of the grid points
                 xgrid(ir + ioffset,ix,iy) = x_current
                 ygrid(ir + ioffset,ix,iy) = y_current
                 zgrid(ir + ioffset,ix,iy) = gamma*z_interface_top + (1.d0 - gamma)*z_interface_bottom
              enddo

           end do
        end do

        ! the top interface becomes the bottom interface before switching to the next layer
        SUPPRESS_UTM_PROJECTION_BOTTOM = SUPPRESS_UTM_PROJECTION_TOP
        npx_interface_bottom = npx_interface_top
        npy_interface_bottom = npy_interface_top
        orig_x_interface_bottom = orig_x_interface_top
        orig_y_interface_bottom = orig_y_interface_top
        spacing_x_interface_bottom = spacing_x_interface_top
        spacing_y_interface_bottom = spacing_y_interface_top
        interface_bottom(:,:) = interface_top(:,:)

     enddo

     close(IIN_INTERFACES)


!!! FIN PLL

! ! fill the region between the cutoff depth and the free surface
!   do iy=0,npy
!   do ix=0,npx

! !   define the mesh points on the top and the bottom

!     xin=dble(ix)/dble(npx)
!     x_current = UTM_X_MIN + (dble(iproc_xi)+xin)*(UTM_X_MAX-UTM_X_MIN)/dble(NPROC_XI)

!     etan=dble(iy)/dble(npy)
!     y_current = UTM_Y_MIN + (dble(iproc_eta)+etan)*(UTM_Y_MAX-UTM_Y_MIN)/dble(NPROC_ETA)

! ! define model between topography surface and fictitious bottom
!     if(TOPOGRAPHY) then

! ! project x and y in UTM back to long/lat since topo file is in long/lat
!   call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

! ! get coordinate of corner in bathy/topo model
!     icornerlong = int((long - ORIG_LONG_TOPO) / DEGREES_PER_CELL_TOPO) + 1
!     icornerlat = int((lat - ORIG_LAT_TOPO) / DEGREES_PER_CELL_TOPO) + 1

! ! avoid edge effects and extend with identical point if outside model
!     if(icornerlong < 1) icornerlong = 1
!     if(icornerlong > NX_TOPO-1) icornerlong = NX_TOPO-1
!     if(icornerlat < 1) icornerlat = 1
!     if(icornerlat > NY_TOPO-1) icornerlat = NY_TOPO-1

! ! compute coordinates of corner
!     long_corner = ORIG_LONG_TOPO + (icornerlong-1)*DEGREES_PER_CELL_TOPO
!     lat_corner = ORIG_LAT_TOPO + (icornerlat-1)*DEGREES_PER_CELL_TOPO

! ! compute ratio for interpolation
!     ratio_xi = (long - long_corner) / DEGREES_PER_CELL_TOPO
!     ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_TOPO

! ! avoid edge effects
!     if(ratio_xi < 0.) ratio_xi = 0.
!     if(ratio_xi > 1.) ratio_xi = 1.
!     if(ratio_eta < 0.) ratio_eta = 0.
!     if(ratio_eta > 1.) ratio_eta = 1.

! ! interpolate elevation at current point
!     elevation = &
!       itopo_bathy(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!       itopo_bathy(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!       itopo_bathy(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!       itopo_bathy(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta

!     else

!       elevation = 0.d0

!     endif

!     z_top = Z_SURFACE + elevation
!     z_bot = - dabs(Z_DEPTH_BLOCK)

! ! compute global min and max of elevation
!   min_elevation = dmin1(min_elevation,elevation)
!   max_elevation = dmax1(max_elevation,elevation)

! ! create vertical point distribution at current horizontal point
!   if(BASEMENT_MAP) then

! ! get coordinate of corner in bathy/topo model
!     icorner_x = int((x_current - ORIG_X_BASEMENT) / SPACING_X_BASEMENT) + 1
!     icorner_y = int((y_current - ORIG_Y_BASEMENT) / SPACING_Y_BASEMENT) + 1

! ! avoid edge effects and extend with identical point if outside model
!     if(icorner_x < 1) icorner_x = 1
!     if(icorner_x > NX_BASEMENT-1) icorner_x = NX_BASEMENT-1
!     if(icorner_y < 1) icorner_y = 1
!     if(icorner_y > NY_BASEMENT-1) icorner_y = NY_BASEMENT-1

! ! compute coordinates of corner
!     x_corner = ORIG_X_BASEMENT + (icorner_x-1)*SPACING_X_BASEMENT
!     y_corner = ORIG_Y_BASEMENT + (icorner_y-1)*SPACING_Y_BASEMENT

! ! compute ratio for interpolation
!     ratio_xi = (x_current - x_corner) / SPACING_X_BASEMENT
!     ratio_eta = (y_current - y_corner) / SPACING_Y_BASEMENT

! ! avoid edge effects
!     if(ratio_xi < 0.) ratio_xi = 0.
!     if(ratio_xi > 1.) ratio_xi = 1.
!     if(ratio_eta < 0.) ratio_eta = 0.
!     if(ratio_eta > 1.) ratio_eta = 1.

! ! interpolate basement surface at current point
!     Z_BASEMENT_SURFACE = &
!       z_basement(icorner_x,icorner_y)*(1.-ratio_xi)*(1.-ratio_eta) + &
!       z_basement(icorner_x+1,icorner_y)*ratio_xi*(1.-ratio_eta) + &
!       z_basement(icorner_x+1,icorner_y+1)*ratio_xi*ratio_eta + &
!       z_basement(icorner_x,icorner_y+1)*(1.-ratio_xi)*ratio_eta

!   else
!     Z_BASEMENT_SURFACE = DEPTH_5p5km_SOCAL
!   endif

! ! honor Lupei Zhu's Moho map
!   if(MOHO_MAP_LUPEI) then

! ! project x and y in UTM back to long/lat since topo file is in long/lat
!     call utm_geo(long,lat,x_current,y_current,UTM_PROJECTION_ZONE,IUTM2LONGLAT,SUPPRESS_UTM_PROJECTION)

! ! get coordinate of corner in Moho map
!     icornerlong = int((long - ORIG_LONG_MOHO) / DEGREES_PER_CELL_MOHO) + 1
!     icornerlat = int((lat - ORIG_LAT_MOHO) / DEGREES_PER_CELL_MOHO) + 1

! ! avoid edge effects and extend with identical point if outside model
!     if(icornerlong < 1) icornerlong = 1
!     if(icornerlong > NX_MOHO-1) icornerlong = NX_MOHO-1
!     if(icornerlat < 1) icornerlat = 1
!     if(icornerlat > NY_MOHO-1) icornerlat = NY_MOHO-1

! ! compute coordinates of corner
!     long_corner = ORIG_LONG_MOHO + (icornerlong-1)*DEGREES_PER_CELL_MOHO
!     lat_corner = ORIG_LAT_MOHO + (icornerlat-1)*DEGREES_PER_CELL_MOHO

! ! compute ratio for interpolation
!     ratio_xi = (long - long_corner) / DEGREES_PER_CELL_MOHO
!     ratio_eta = (lat - lat_corner) / DEGREES_PER_CELL_MOHO

! ! avoid edge effects
!     if(ratio_xi < 0.) ratio_xi = 0.
!     if(ratio_xi > 1.) ratio_xi = 1.
!     if(ratio_eta < 0.) ratio_eta = 0.
!     if(ratio_eta > 1.) ratio_eta = 1.

! ! interpolate Moho depth at current point
!     Z_DEPTH_MOHO = &
!      - (imoho_depth(icornerlong,icornerlat)*(1.-ratio_xi)*(1.-ratio_eta) + &
!         imoho_depth(icornerlong+1,icornerlat)*ratio_xi*(1.-ratio_eta) + &
!         imoho_depth(icornerlong+1,icornerlat+1)*ratio_xi*ratio_eta + &
!         imoho_depth(icornerlong,icornerlat+1)*(1.-ratio_xi)*ratio_eta)

!   else
!     Z_DEPTH_MOHO = DEPTH_MOHO_SOCAL
!   endif

! ! define vertical spacing of the mesh in case of a non-regular mesh with mesh doublings
!   if(.not. USE_REGULAR_MESH) call mesh_vertical(myrank,rns,NER,NER_BOTTOM_MOHO,NER_MOHO_16, &
!                      NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
! !! DK DK UGLY modif z_top by Emmanuel Chaljub here
! !! DK DK UGLY modif Manu removed                     z_top, &
!                      Z_DEPTH_BLOCK,Z_BASEMENT_SURFACE,Z_DEPTH_MOHO,MOHO_MAP_LUPEI)

! !   fill the volume
!     do ir = 0,2*NER
!       if(USE_REGULAR_MESH) then
!         rn = dble(ir) / dble(2*NER)
!       else
!         rn = rns(ir)
!       endif
!       xgrid(ir,ix,iy) = x_current
!       ygrid(ir,ix,iy) = y_current
!       zgrid(ir,ix,iy) = z_bot*(ONE-rn) + z_top*rn
!     enddo

!   enddo
!   enddo

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! volume of bottom and top area of the slice
!   volume_local = ZERO
!   area_local_bottom = ZERO
!   area_local_top = ZERO

! assign theoretical number of elements
  nspec = NSPEC_AB

! compute maximum number of points
  npointot = nspec * NGLLCUBE

! make sure everybody is synchronized
  call sync_all()

! use dynamic allocation to allocate memory for arrays
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)

! exit if there is not enough memory to allocate all the arrays
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

! create all the regions of the mesh
!   call create_regions_mesh(xgrid,ygrid,zgrid,ibool,idoubling, &
!          xstore,ystore,zstore,npx,npy, &
!          iproc_xi,iproc_eta,addressing,nspec, &
!          volume_local,area_local_bottom,area_local_top, &
!          NGLOB_AB,npointot, &
!          NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM,NER, &
!          NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
!          NSPEC2DMAX_XMIN_XMAX, &
!          NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!          HARVARD_3D_GOCAD_MODEL,NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI, &
!          NSPEC2D_A_ETA,NSPEC2D_B_ETA,myrank,LOCAL_PATH, &
!          UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,UTM_PROJECTION_ZONE, &
!          HAUKSSON_REGIONAL_MODEL,OCEANS, &
!          VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM, &
!          IMPOSE_MINIMUM_VP_GOCAD,THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,MOHO_MAP_LUPEI, &
!          ANISOTROPY,SAVE_MESH_FILES,SUPPRESS_UTM_PROJECTION, &
!          ORIG_LAT_TOPO,ORIG_LONG_TOPO,DEGREES_PER_CELL_TOPO,NX_TOPO,NY_TOPO,USE_REGULAR_MESH)

! call create_regions_mesh(xgrid,ygrid,zgrid,ibool,idoubling, &
!            xstore,ystore,zstore,npx,npy,iproc_xi,iproc_eta,addressing,nspec, &
!            NGLOB_AB,npointot, &
!            NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM,NER, &
!            NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
!            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!            HARVARD_3D_GOCAD_MODEL,NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI, &
!            NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
!            myrank,LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
!            HAUKSSON_REGIONAL_MODEL,USE_REGULAR_MESH)

call create_regions_mesh(xgrid,ygrid,zgrid,ibool, &
           xstore,ystore,zstore,npx,npy,iproc_xi,iproc_eta,addressing,nspec, &
           NGLOB_AB,npointot, &
           NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
           NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI, &
           NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
           NSUBREGIONS,subregions,number_of_layers,ner_layer,NMATERIALS,material_properties, &
           myrank,LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,&
           USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

! ! print min and max of topography included
!   if(TOPOGRAPHY) then

! ! compute the maximum of the maxima for all the slices using an MPI reduction
!       call min_all_dp(min_elevation,min_elevation_all)
!       call max_all_dp(max_elevation,max_elevation_all)

!     if(myrank == 0) then
!       write(IMAIN,*)
!       write(IMAIN,*) 'min and max of topography included in mesh in m is ',min_elevation_all,' ',max_elevation_all
!       write(IMAIN,*)
!     endif
!   endif


! ! use MPI reduction to compute total area and volume
!   area_total_bottom   = ZERO
!   area_total_top   = ZERO
!   call sum_all_dp(area_local_bottom,area_total_bottom)
!   call sum_all_dp(area_local_top,area_total_top)
!   call sum_all_dp(volume_local,volume_total)

  if(myrank == 0) then

! !   check volume, and bottom and top area

!       write(IMAIN,*)
!       write(IMAIN,*) '   calculated top area: ',area_total_top

! ! compare to exact theoretical value
!     if(.not. TOPOGRAPHY) &
!           write(IMAIN,*) '            exact area: ',(UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)

!       write(IMAIN,*)
!       write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

! compare to exact theoretical value (bottom is always flat)
      write(IMAIN,*) '            exact area: ',(UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)

  endif

! make sure everybody is synchronized
  call sync_all()

!   if(myrank == 0) then
! ! check volume
!       write(IMAIN,*)
!       write(IMAIN,*) 'calculated volume: ',volume_total
! ! take the central cube into account
!    if(.not. TOPOGRAPHY) &
!       write(IMAIN,*) '     exact volume: ', &
!         (UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)*dabs(Z_DEPTH_BLOCK)

!   endif

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

! write information about precision used for floating-point operations
  if(CUSTOM_REAL == SIZE_REAL) then
    write(IMAIN,*) 'using single precision for the calculations'
  else
    write(IMAIN,*) 'using double precision for the calculations'
  endif
  write(IMAIN,*)
  write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
  write(IMAIN,*)

! ! copy number of elements and points in an include file for the solver
!   call save_header_file(NSPEC_AB,NGLOB_AB,NEX_XI,NEX_ETA,NPROC,NPROC_XI,NPROC_ETA, &
!              UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,ATTENUATION,ANISOTROPY,NSTEP, &
!              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!              NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,SIMULATION_TYPE)

!   call get_value_string(rec_filename, 'solver.STATIONS', 'DATA/STATIONS')
!   call get_value_string(filtered_rec_filename, 'solver.STATIONS_FILTERED', 'DATA/STATIONS_FILTERED')

! ! get total number of stations
!   open(unit=IIN,file=rec_filename,iostat=ios,status='old',action='read')
!   nrec = 0
!   do while(ios == 0)
!     read(IIN,"(a)",iostat=ios) dummystring
!     if(ios == 0) nrec = nrec + 1
!   enddo
!   close(IIN)

! ! filter list of stations, only retain stations that are in the model
!   nrec_filtered = 0
!   open(unit=IIN,file=rec_filename,status='old',action='read')
!   do irec = 1,nrec
!     read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
!     if(stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
!       nrec_filtered = nrec_filtered + 1
!   enddo
!   close(IIN)

!   write(IMAIN,*)
!   write(IMAIN,*) 'there are ',nrec,' stations in file ', trim(rec_filename)
!   write(IMAIN,*) 'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_rec_filename)
!   write(IMAIN,*) 'excluding ',nrec - nrec_filtered,' stations located outside the model'
!   write(IMAIN,*)

!   if(nrec_filtered < 1) call exit_MPI(myrank,'need at least one station in the model')

!   open(unit=IIN,file=rec_filename,status='old',action='read')
!   open(unit=IOUT,file=filtered_rec_filename,status='unknown')

!   do irec = 1,nrec
!     read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
!     if(stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
!       write(IOUT,*) station_name(1:len_trim(station_name)),' ',network_name(1:len_trim(network_name)),' ', &
!               sngl(stlat),' ',sngl(stlon), ' ', sngl(stele), ' ', sngl(stbur)
!   enddo

!   close(IIN)
!   close(IOUT)

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

  end subroutine meshfem3D

