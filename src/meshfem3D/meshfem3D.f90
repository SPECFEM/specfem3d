!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
  include "constants_meshfem3D.h"

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
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
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
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
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
! MPI v. 2.1, July 2012:
! Max Rietmann, Peter Messmer, Daniel Peter, Dimitri Komatitsch, Joseph Charles, Zhinan Xie:
! support for CUDA GPUs, better CFL stability for the Stacey absorbing conditions.
!
! MPI v. 2.0, November 2010:
! Dimitri Komatitsch, Nicolas Le Goff, Roland Martin and Pieyre Le Loher, University of Pau, France,
! Daniel Peter, Jeroen Tromp and the Princeton group of developers, Princeton University, USA,
! and Emanuele Casarotti, INGV Roma, Italy:
!  support for CUBIT meshes decomposed by SCOTCH;
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
! MPI v. 1.0 Dimitri Komatitsch, Caltech, USA, May 2002: first MPI version based on global code
!
! Dimitri Komatitsch, IPG Paris, France, December 1996: first 3-D solver for the CM-5 Connection Machine,
!    parallelized on 128 processors using Connection Machine Fortran
!

! number of spectral elements in each block
  integer nspec,npointot

! meshing parameters
  double precision, dimension(:), allocatable :: rns

! auxiliary variables to generate the mesh
  integer ix,iy,ir

  double precision xin,etan
  double precision x_current,y_current

  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid

  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier
  integer iprocnum,npx,npy

! for loop on all the slices
  integer iproc_xi,iproc_eta
  integer, dimension(:,:), allocatable :: addressing

! use integer array to store topography values
  integer icornerlat,icornerlong
  double precision lat,long
  double precision long_corner,lat_corner,ratio_xi,ratio_eta

! timer MPI
  double precision, external :: wtime
  double precision time_start,tCPU

! addressing for all the slices
  integer, dimension(:), allocatable :: iproc_xi_slice,iproc_eta_slice

! parameters read from mesh parameter file
  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX
  double precision Z_DEPTH_BLOCK
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  logical SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH

! Mesh files for visualization
  logical CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

! doublings parameters
  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

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

! interfaces parameters
  logical SUPPRESS_UTM_PROJECTION_BOTTOM,SUPPRESS_UTM_PROJECTION_TOP
  integer ilayer,interface_current ! ipoint_current
  integer number_of_interfaces,number_of_layers
  integer max_npx_interface,max_npy_interface
  integer npx_interface_bottom,npy_interface_bottom
  integer npx_interface_top,npy_interface_top
  double precision z_interface_bottom,z_interface_top
  double precision orig_x_interface_bottom,orig_y_interface_bottom
  double precision orig_x_interface_top,orig_y_interface_top
  double precision spacing_x_interface_bottom,spacing_y_interface_bottom
  double precision spacing_x_interface_top,spacing_y_interface_top
  character(len=50) INTERFACES_FILE,interface_top_file
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

  ! parameters read from parameter file
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,SIMULATION_TYPE
  integer NSOURCES,NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO,NGNOD,NGNOD2D
  double precision DT
  double precision HDUR_MOVIE,OLSEN_ATTENUATION_RATIO,f0_FOR_PML
  logical ATTENUATION,USE_OLSEN_ATTENUATION, &
          APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,USE_FORCE_POINT_SOURCE
  logical STACEY_ABSORBING_CONDITIONS,SAVE_FORWARD,STACEY_INSTEAD_OF_FREE_SURFACE
  logical ANISOTROPY,SAVE_MESH_FILES,USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION
  logical PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE,FULL_ATTENUATION_SOLID
  integer MOVIE_TYPE,IMODEL
  character(len=256) OUTPUT_FILES,LOCAL_PATH,TOMOGRAPHY_PATH,TRAC_PATH
  logical :: ADIOS_ENABLED, ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, &
             ADIOS_FOR_FORWARD_ARRAYS, ADIOS_FOR_KERNELS

! ************** PROGRAM STARTS HERE **************

! sizeprocs returns number of processes started (should be equal to NPROC).
! myrank is the rank of each process, between 0 and NPROC-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', OUTPUT_FILES_PATH(1:len_trim(OUTPUT_FILES_PATH)))

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
  call read_parameter_file(NPROC,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,DT,NGNOD,NGNOD2D, &
                        UTM_PROJECTION_ZONE,SUPPRESS_UTM_PROJECTION,TOMOGRAPHY_PATH, &
                        ATTENUATION,USE_OLSEN_ATTENUATION,LOCAL_PATH,NSOURCES, &
                        APPROXIMATE_OCEAN_LOAD,TOPOGRAPHY,ANISOTROPY,STACEY_ABSORBING_CONDITIONS,MOVIE_TYPE, &
                        MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
                        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,HDUR_MOVIE, &
                        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
                        NTSTEP_BETWEEN_OUTPUT_INFO,SIMULATION_TYPE,SAVE_FORWARD, &
                        NTSTEP_BETWEEN_READ_ADJSRC,NOISE_TOMOGRAPHY, &
                        USE_FORCE_POINT_SOURCE,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        USE_RICKER_TIME_FUNCTION,OLSEN_ATTENUATION_RATIO,PML_CONDITIONS, &
                        PML_INSTEAD_OF_FREE_SURFACE,f0_FOR_PML,IMODEL,FULL_ATTENUATION_SOLID,TRAC_PATH)

  call read_adios_parameters(ADIOS_ENABLED, ADIOS_FOR_DATABASES,       &
                             ADIOS_FOR_FORWARD_ARRAYS, ADIOS_FOR_MESH, &
                             ADIOS_FOR_KERNELS)

! read the mesh parameter file
! nullify(subregions,material_properties)
  call read_mesh_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
                          UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
                          NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE, &
                          LOCAL_PATH,SUPPRESS_UTM_PROJECTION,&
                          INTERFACES_FILE,NSUBREGIONS,subregions,NMATERIALS,material_properties, &
                          CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                          USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) &
    stop 'error: must have NPROC_XI = NPROC_ETA = 1 for a serial run'

! get interface data from external file to count the spectral elements along Z
  if(myrank == 0) then
     write(IMAIN,*) 'Reading interface data from file ',&
      MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH))//INTERFACES_FILE(1:len_trim(INTERFACES_FILE)), &
      ' to count the spectral elements'
  endif

  open(unit=IIN,file=MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH))//INTERFACES_FILE, &
      status='old')

  max_npx_interface  = -1
  max_npy_interface  = -1

! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES')
  if(number_of_interfaces < 1) stop 'error: not enough interfaces (minimum is 1, for topography)'

! loop on all the interfaces
  do interface_current = 1,number_of_interfaces
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_BOTTOM,interface_top_file, &
          npx_interface_bottom,npy_interface_bottom,&
          orig_x_interface_bottom,orig_y_interface_bottom,spacing_x_interface_bottom,spacing_y_interface_bottom)

    max_npx_interface = max(npx_interface_bottom,max_npx_interface)
    max_npy_interface = max(npy_interface_bottom,max_npy_interface)

    if((max_npx_interface < 2) .or.(max_npy_interface < 2)) stop 'not enough interface points (minimum is 2x2)'

  enddo

  ! define number of layers
  number_of_layers = number_of_interfaces! - 1
  allocate(ner_layer(number_of_layers),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ner_layer'

! loop on all the layers
  do ilayer = 1,number_of_layers

! read number of spectral elements in vertical direction in this layer
    call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,ner_layer(ilayer),'NER_LAYER')
    if(ner_layer(ilayer) < 1) stop 'not enough spectral elements along Z in layer (minimum is 1)'

  enddo

  close(IIN)

! compute total number of spectral elements in vertical direction
  NER = sum(ner_layer)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                        NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                        NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
                        NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB,&
                        USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) then
    if( myrank == 0 ) then
      write(IMAIN,*) 'error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'error: number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'error meshfem3D: number of processors supposed to run on: ',NPROC
      print*, 'error meshfem3D: number of MPI processors actually run on: ',sizeprocs
      print*
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

! dynamic allocation of mesh arrays
  allocate(rns(0:2*NER),stat=ier)
  if( ier /= 0 ) stop 'error allocating array rns'

  allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if( ier /= 0 ) stop 'error allocating array xgrid'
  allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ygrid'
  allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'not enough memory to allocate arrays')

  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1),stat=ier)
  if( ier /= 0 ) stop 'error allocating array addressing'
  allocate(iproc_xi_slice(0:NPROC-1),stat=ier)
  if( ier /= 0 ) stop 'error allocating array iproc_xi_slice'
  allocate(iproc_eta_slice(0:NPROC-1),stat=ier)
  if( ier /= 0 ) stop 'error allocating array iproc_eta_slice'

! clear arrays
  xgrid(:,:,:) = 0.d0
  ygrid(:,:,:) = 0.d0
  zgrid(:,:,:) = 0.d0

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
    enddo
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices'

    write(IMAIN,*)
    write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD_EIGHT_CORNERS,' control nodes'
    write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D_FOUR_CORNERS,' control nodes'
    write(IMAIN,*) 'Beware! Curvature (i.e. HEX27 elements) is not handled by our internal mesher'
    write(IMAIN,*)
  endif

  ! check that the constants.h file is correct
  if(NGNOD /= 8) call exit_MPI(myrank,'volume elements should have 8 control nodes in our internal mesher')
  if(NGNOD2D /= 4) call exit_MPI(myrank,'surface elements should have 4 control nodes in our internal mesher')

  ! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

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
  endif

  ! checks that nex is dividable by 8
  if(mod(NEX_XI,8) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8')
  if(mod(NEX_ETA,8) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8')

  ! checks that nex_per_proc is dividable by 8
  if(mod(NEX_PER_PROC_XI,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 8')
  if(mod(NEX_PER_PROC_ETA,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 8')


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
    if(SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'suppressing UTM projection'
    else
      write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
    endif
    write(IMAIN,*)
  endif

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
    write(IMAIN,*) 'Reading interface data from file ', &
          MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH)) &
          //INTERFACES_FILE(1:len_trim(INTERFACES_FILE))
    write(IMAIN,*)
  endif

  open(unit=IIN,file=MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH)) &
          //INTERFACES_FILE,status='old')

  allocate(interface_bottom(max_npx_interface,max_npy_interface),stat=ier)
  if( ier /= 0 ) stop 'error allocating array interface_bottom'
  allocate(interface_top(max_npx_interface,max_npy_interface),stat=ier)
  if( ier /= 0 ) stop 'error allocating array interface_top'

  ! read number of interfaces
  call read_value_integer_mesh(IIN,DONT_IGNORE_JUNK,number_of_interfaces,'NINTERFACES')

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
    call read_interface_parameters(IIN,SUPPRESS_UTM_PROJECTION_TOP,interface_top_file,&
         npx_interface_top,npy_interface_top,&
         orig_x_interface_top,orig_y_interface_top,spacing_x_interface_top,spacing_y_interface_top)

    !npoints_interface_top = npx_interface_top * npy_interface
    ! loop on all the points describing this interface
    open(unit=45,file=MF_IN_DATA_FILES_PATH(1:len_trim(MF_IN_DATA_FILES_PATH)) &
         //interface_top_file,status='old')
    do iy=1,npy_interface_top
      do ix=1,npx_interface_top
        call read_value_dble_precision_mesh(45,DONT_IGNORE_JUNK,interface_top(ix,iy),'Z_INTERFACE_TOP')
      enddo
    enddo
    close(45)

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
        icornerlat = int((lat - orig_y_interface_bottom) / spacing_y_interface_bottom) + 1

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
        icornerlat = int((lat - orig_y_interface_top) / spacing_y_interface_top) + 1

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

      enddo
    enddo

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

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
  endif

! assign theoretical number of elements
  nspec = NSPEC_AB

! compute maximum number of points
  npointot = nspec * NGLLCUBE_M

! make sure everybody is synchronized
  call sync_all()

! use dynamic allocation to allocate memory for arrays
  allocate(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ibool'
  allocate(xstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating array xstore'
  allocate(ystore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating array ystore'
  allocate(zstore(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
  ! exit if there is not enough memory to allocate all the arrays
  if(ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  call create_regions_mesh(xgrid,ygrid,zgrid,ibool, &
                         xstore,ystore,zstore,iproc_xi,iproc_eta,addressing,nspec, &
                         NGLOB_AB,npointot, &
                         NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NER, &
                         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                         NPROC_XI,NPROC_ETA, &
                         NSUBREGIONS,subregions,number_of_layers,ner_layer,NMATERIALS,material_properties, &
                         myrank, sizeprocs, &
                         LOCAL_PATH,UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK,&
                         CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                        USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
                        ADIOS_ENABLED, ADIOS_FOR_DATABASES)

  if(myrank == 0) then
! compare to exact theoretical value (bottom is always flat)
    write(IMAIN,*) '            exact area: ',(UTM_Y_MAX-UTM_Y_MIN)*(UTM_X_MAX-UTM_X_MIN)
  endif

! make sure everybody is synchronized
  call sync_all()

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

