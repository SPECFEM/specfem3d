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
!
! United States and French Government Sponsorship Acknowledged.
!

  subroutine meshfem3D()

  use meshfem3D_par
  use chunk_earth_mod

  implicit none

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
! journal = {communications in Computational Physics},
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
! journal = {communications in Computational Physics},
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
! To report bugs or suggest improvements to the code, please use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! MPI v. 3.0, December 2014: many developers.
! Convolutional PML, LDDRK time scheme, bulk attenuation support, simultaneous MPI runs,
! ADIOS file I/O support, coupling with external codes, new seismogram names,
! Deville routines for additional GLL degrees, tomography tools, unit/regression test framework,
! improved CUDA GPUs performance, additonal GEOCUBIT support, better make compilation,
! git versioning system.
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

! local parameters
  integer :: iprocnum
  integer :: iproc_xi,iproc_eta
  integer :: ier

  ! interface parameters
  integer :: ilayer

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU

! ************** PROGRAM STARTS HERE **************

  ! sizeprocs returns number of processes started (should be equal to NPROC).
  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) then
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_meshfem3D.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open output file :',trim(OUTPUT_FILES)//'/output_meshfem3D.txt'
      stop 'Error opening output file'
    endif
  endif

  ! get MPI starting time
  time_start = wtime()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '******************************************'
    write(IMAIN,*) '*** Specfem3D MPI meshfem3D - f90 version ***'
    write(IMAIN,*) '******************************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read the parameter file (DATA/Par_file)
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(myrank,BROADCAST_AFTER_READ)

  ! make sure everybody is synchronized
  call synchronize_all()

  ! if meshing a chunk of the Earth, call a specific internal mesher designed specifically for that
  ! CD CD change this to have also the possibility to use a chunk without coupling
  if (MESH_A_CHUNK_OF_THE_EARTH) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'creating chunk of the Earth mesh'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    !! VM VM : new way to mesh and store mesh in geocubit format
    if (myrank == 0) then  !! serial mesh and use decompose_mesh after
       call mesh_chunk_earth()
       write(*,*) 'Done creating a chunk of the earth Mesh (HEX8 elements), see directory MESH/'
    endif
    ! make sure everybody is synchronized
    call synchronize_all()
    return
    !call bcast_input_param_to_all()
    !call read_mesh_parameter_file()


    !! VM VM old way to create  MESH_A_CHUNK_OF_THE_EARTH, but still some routines that will
    !! move in the new way thus not remove for now
    if (NGNOD == 8) then
      ! creates mesh in MESH/
      call earth_chunk_HEX8_Mesher(NGNOD)
      ! done with mesher
      stop 'Done creating a chunk of the earth Mesh (HEX8 elements), see directory MESH/'

    else if (NGNOD == 27) then

      ! creates mesh in MESH/
      call earth_chunk_HEX27_Mesher(NGNOD)
      ! done with mesher
      stop 'Done creating a chunk of the earth Mesh (HEX27 elements), see directory MESH/'

    else
      stop 'Bad number of nodes per hexahedron: NGNOD must be equal to 8 or 27'
    endif

    ! make sure everybody is synchronized
    call synchronize_all()
  endif

  ! read the mesh parameter file (Data/meshfem3D_files/Mesh_Par_file)
  ! nullify(subregions,material_properties)
  call read_mesh_parameter_file()

  ! get interface data from external file to count the spectral elements along Z
  call get_interfaces_mesh_count()

  ! compute total number of spectral elements in vertical direction
  NER = sum(ner_layer)

  ! checks if regions and vertical layers from interfaces file match
  if (maxval(subregions(:,6)) /= NER) then
    print *,'Error invalid total number of element layers in vertical direction!'
    print *,'from interface file, total layers = ',NER
    print *,'should be equal to maximum layer NZ_END specified in regions:', maxval(subregions(:,6))
    stop 'Error invalid total number of vertical layers'
  endif

  ! checks irregular grid entries
  if (.not. USE_REGULAR_MESH) then
    if (maxval(ner_doublings(:)) == NER) then
      print *,'Error invalid doubling layer NZ index too close to surface layer ',NER
      print *,'Please decrease maximum doubling layer index NZ_DOUBLING'
      stop 'Error invalid doubling layer index too close to surface layer'
    endif
  endif


  ! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                          NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC_AB,NSPEC2D_A_XI,NSPEC2D_B_XI, &
                          NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
                          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                          NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NGLOB_AB, &
                          USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      write(IMAIN,*) 'Error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'Error: number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'Error meshfem3D: number of processors supposed to run on: ',NPROC
      print *, 'Error meshfem3D: number of MPI processors actually run on: ',sizeprocs
      print *
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

  ! dynamic allocation of mesh arrays
  allocate(rns(0:2*NER),stat=ier)
  if (ier /= 0) stop 'Error allocating array rns'

  allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) stop 'Error allocating array xgrid'
  allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) stop 'Error allocating array ygrid'
  allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')

  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1),stat=ier)
  if (ier /= 0) stop 'Error allocating array addressing'
  allocate(iproc_xi_slice(0:NPROC-1),stat=ier)
  if (ier /= 0) stop 'Error allocating array iproc_xi_slice'
  allocate(iproc_eta_slice(0:NPROC-1),stat=ier)
  if (ier /= 0) stop 'Error allocating array iproc_eta_slice'

  ! clear arrays
  xgrid(:,:,:) = 0.d0
  ygrid(:,:,:) = 0.d0
  zgrid(:,:,:) = 0.d0

  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

  ! create global slice addressing for solver
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Creating global slice addressing'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  do iproc_eta = 0,NPROC_ETA-1
    do iproc_xi = 0,NPROC_XI-1
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
    call flush_IMAIN()
  endif

  if (myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta'
    write(IMAIN,*) 'There are ',NER,' elements along Z'
    write(IMAIN,*)
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
    call flush_IMAIN()
  endif

  ! check that the constants.h file is correct
  if (NGNOD /= 8) call exit_MPI(myrank,'volume elements should have 8 control nodes in our internal mesher')
  if (NGNOD2D /= 4) call exit_MPI(myrank,'surface elements should have 4 control nodes in our internal mesher')

  ! check that reals are either 4 or 8 bytes
  if (CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

  ! for the number of standard linear solids for attenuation
  if (N_SLS /= 3) call exit_MPI(myrank,'number of SLS must be 3')

  ! check that number of slices is at least 1 in each direction
  if (NPROC_XI < 1) call exit_MPI(myrank,'NPROC_XI must be greater than 1')
  if (NPROC_ETA < 1) call exit_MPI(myrank,'NPROC_ETA must be greater than 1')

  ! check that mesh can be cut into the right number of slices
  if (mod(NEX_XI,NPROC_XI) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of NPROC_XI for a regular mesh')
  if (mod(NEX_ETA,NPROC_ETA) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of NPROC_ETA for a regular mesh')

  ! also check that mesh can be coarsened in depth twice (block size multiple of 8)
  ! i.e. check that NEX is divisible by 8 and that NEX_PER_PROC is divisible by 8
  ! This is not required for a regular mesh
  if (.not. USE_REGULAR_MESH) then
    if (mod(NEX_XI,8) /= 0) call exit_MPI(myrank,'NEX_XI must be a multiple of 8')
    if (mod(NEX_ETA,8) /= 0) call exit_MPI(myrank,'NEX_ETA must be a multiple of 8')

    if (mod(NEX_PER_PROC_XI,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 8')
    if (mod(NEX_PER_PROC_ETA,8) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 8')

    if (mod(NEX_PER_PROC_XI, 2**NDOUBLINGS * 2) /= 0 ) &
      call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 2 * 2**NDOUBLINGS')
    if (mod(NEX_PER_PROC_ETA, 2**NDOUBLINGS * 2) /= 0 ) &
      call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 2 * 2**NDOUBLINGS')
  endif

  if (myrank == 0) then
    write(IMAIN,*) 'region selected:'
    write(IMAIN,*)
    write(IMAIN,*) 'latitude min = ',LATITUDE_MIN
    write(IMAIN,*) 'latitude max = ',LATITUDE_MAX
    write(IMAIN,*)
    write(IMAIN,*) 'longitude min = ',LONGITUDE_MIN
    write(IMAIN,*) 'longitude max = ',LONGITUDE_MAX
    write(IMAIN,*)
    if (SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'this is given directly as UTM'
    else
      write(IMAIN,*) 'this is mapped to UTM in region ',UTM_PROJECTION_ZONE
    endif
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
    if (SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,*) 'suppressing UTM projection'
    else
      write(IMAIN,*) 'using UTM projection in region ',UTM_PROJECTION_ZONE
    endif
     if (PML_CONDITIONS) then
       write(IMAIN,*)
       write(IMAIN,*) 'PML thickness in X direction = ',THICKNESS_OF_X_PML,'m'
       write(IMAIN,*) 'PML thickness in Y direction = ',THICKNESS_OF_Y_PML,'m'
       write(IMAIN,*) 'PML thickness in Z direction = ',THICKNESS_OF_Z_PML,'m'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get addressing for this process
  iproc_xi_current = iproc_xi_slice(myrank)
  iproc_eta_current = iproc_eta_slice(myrank)

  ! number of elements in each slice
  npx_element_steps = 2*NEX_PER_PROC_XI
  npy_element_steps = 2*NEX_PER_PROC_ETA
  ner_layer(:) = 2 * ner_layer(:)

  !min_elevation = +HUGEVAL
  !max_elevation = -HUGEVAL

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'Creating interfaces'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! creates mesh interfaces
  call create_interfaces_mesh()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'Creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! creates mesh element points
  call create_meshfem_mesh()

  ! make sure everybody is synchronized
  call synchronize_all()

  !--- print number of points and elements in the mesh
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Repartition of elements:'
    write(IMAIN,*) '-----------------------'
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in mesh slice 0: ',NSPEC_AB
    write(IMAIN,*) 'total number of points in mesh slice 0: ',NGLOB_AB
    write(IMAIN,*)
    write(IMAIN,*) 'total number of elements in entire mesh: ',NSPEC_AB*NPROC
    ! the float() statements here are for the case of more than 2 Gigapoints per mesh, in which
    ! case and integer(kind=4) counter would overflow and display an incorrect negative value;
    ! converting it to float fixes the problem (but then prints some extra decimals equal to zero).
    ! Another option would be to declare the sum as integer(kind=8) and then print it.
    write(IMAIN,*) 'approximate total number of points in entire mesh (with duplicates on MPI edges): ', &
                                             dble(NGLOB_AB)*dble(NPROC)
    write(IMAIN,*) 'approximate total number of DOFs in entire mesh (with duplicates on MPI edges): ', &
                                             dble(NGLOB_AB)*dble(NPROC*NDIM)
    write(IMAIN,*)
    ! write information about precision used for floating-point operations
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)
    call flush_IMAIN()
  endif   ! end of section executed by main process only

  ! elapsed time since beginning of mesh generation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    call flush_IMAIN()

    ! close main output file
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine meshfem3D

!=====================================================================

  subroutine earth_chunk_HEX8_Mesher(NGNOD)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM, R_EARTH, PI, ZERO, TINYVAL, &
    old_DSM_coupling_from_Vadim, INJECTION_TECHNIQUE_IS_AXISEM, INJECTION_TECHNIQUE_IS_DSM

  use shared_parameters, only: INJECTION_TECHNIQUE_TYPE

  implicit none

!==============================================================================================!
!                                                                                              !
!  Singular option of meshfem3D : MESH OF A GLOBE EARTH CHUNK FOR THE INTERFACE DSM-SPECFEM3D  !
!  Case of 8 nodes per element (HEX8)                                                          !
!                                                                                              !
!  VM, February 2013                                                                           !
!  Integrated in meshfem3d by CD, September 2014                                               !
!                                                                                              !
!  WARNING : A local convention is used for the mapping of                                     !
!            the cubic sphere (to complete)                                                    !
!                                                                                              !
!==============================================================================================!

!
!--- Parameters
!

  integer, parameter :: myrank = 0
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)

  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0

  logical, parameter ::  RUN_BENCHMARK = .false.

!
!--- Other
!

  integer NGNOD

  integer  nel_lat, nel_lon, nel_depth, NX, NY, NZ, Ndepth, nglob, kglob, ilocnum, ieoff, npointot
  integer ilat, ilon, ispec, iz, i, j, k, nspec, ia, izshift, index_mat
  integer ispec2Dxmin, ispec2Dxmax, ispec2Dymin, ispec2Dymax, ispec2Dzmin, ispec2Dzmax
  integer ilayer_current, ilayer
  integer nlat_dsm, nlon_dsm

  integer iaddx(NGNOD), iaddy(NGNOD), iaddz(NGNOD)

  integer, allocatable :: inum_loc(:,:,:,:), iglob(:), loc(:), current_layer(:)

  double precision ratio_eta, ratio_xi
  double precision ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD, Z_DEPTH_BLOCK, UTM_X_MIN, UTM_X_MAX
  double precision lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision deg2rad
  double precision x, y, z, px, py, pz, z_bottom

  double precision rotation_matrix(3,3)
  double precision zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  double precision xelm(NGNOD), yelm(NGNOD), zelm(NGNOD)
  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)

  !! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  !! GLL points and weights of integration
  double precision xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ), wxgll(NGLLX), wygll(NGLLY), wzgll(NGLLZ)

  double precision, allocatable :: xp(:), yp(:), zp(:), xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  double precision, allocatable :: lon_zmin(:,:), lat_zmin(:,:)
  double precision, dimension(:,:), allocatable :: ProfForGemini


  !! For new outputs (list of ggl on boundary, spherical or Cartesian)
  !! AND for coupling with AxiSEM
  integer ::  istore_for_new_outputs
  integer ::   updown(NGLLZ)
  double precision , dimension(NGLLX,NGLLY,NGLLZ) ::  longitud, latitud, radius

  logical test

  logical, allocatable :: ifseg(:)
  logical, dimension(:,:), allocatable :: iboun ! boundary locator

  character(len=100) line
  character(len=250) model1D_file

  character(len=10), parameter :: MESH = "./MESH/"

!! Unused
!
! integer iii, jjj, kkk,
! double precision long, lati, x_bot, y_bot, z_bot,rayon,ratio_depth,
! double precision theta, colat, vector_ori(3), vector_rotated(3)

!
!--- WARNING, CONVENTION : (lon,lat) -> (xi,eta)
!---                       (k = 6 with -z for the mapping of the cubic sphere, cf Chevrot 2012)
!---                       We define the mesh of a chunk of the earth in the cubic sphere
!

  deg2rad = 3.141592653589793d0/180.d0

  open(49, file=trim(MESH)//'output_mesher_chunk_HEX8.txt')

  if (RUN_BENCHMARK) then
!
!--- Parameters fixed inside the subroutine for the moment
!
     ANGULAR_WIDTH_ETA_RAD = 10.d0 * deg2rad ! latitude 2.5
     ANGULAR_WIDTH_XI_RAD  = 20.d0 * deg2rad ! longitude 2.0

!    Chunk center
     lat_center_chunk      = 0.d0  ! 42.35d0 ! 42.5d0 !* deg2rad
     lon_center_chunk      = 60.d0 ! 1.3d0   ! 1.2d0  !* deg2rad

!    Azimuth
     chunk_azi             = 0.d0 !90.d0 !80.d0 !10.d0  !* deg2rad

!    Depth
     chunk_depth           = 1000.d0 * 1000.d0 ! 250.d0 * 1000.d0

!    Number of elements
     nel_lat               = 20 ! 120  ! 15
     nel_lon               = 40 ! 96   ! 15
     nel_depth             = 20 ! 100  ! 10

  else

     open(10, file=trim(MESH)//'ParFileMeshChunk')

     read(10,'(a)') line
     read(10,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
     read(10,'(a)') line
     read(10,*) lon_center_chunk, lat_center_chunk, chunk_azi
     read(10,'(a)') line
     read(10,*) chunk_depth
     read(10,'(a)') line
     read(10,*) nel_lon,nel_lat, nel_depth
     read(10,'(a)') line
     read(10,'(a)') model1D_file

     model1D_file = 'MESH/'//trim(model1D_file)

     close(10)

     ANGULAR_WIDTH_XI_RAD  = deg2rad * ANGULAR_WIDTH_XI_RAD
     ANGULAR_WIDTH_ETA_RAD = deg2rad * ANGULAR_WIDTH_ETA_RAD
     chunk_depth           = chunk_depth * 1000.d0

  endif

  NX = nel_lon
  NY = nel_lat
  NZ = nel_depth

!
!===========================================================================
!
!--- TO DO : the reference chunk must be always symmetric (EW) and (NS)
!

  nlon_dsm = (ngllx - 1) * NX + 1
  nlat_dsm = (nglly - 1) * NY + 1
  nglob    = (nel_lat + 1) * (nel_lon + 1) * (nel_depth + 1)
  nspec    = nel_lat * nel_lon * nel_depth
  npointot = 8 * nspec

  allocate(xp(npointot), yp(npointot), zp(npointot))
  allocate(iglob(npointot), loc(npointot))
  allocate(ifseg(npointot))
  allocate(ProfForGemini(0:NZ-1,3))
  allocate(current_layer(0:NZ-1))
  allocate(inum_loc(2,2,2,nspec))
  allocate(xgrid(2,2,2,nspec), ygrid(2,2,2,nspec), zgrid(2,2,2,nspec))
  allocate(lon_zmin(nlon_dsm,nlat_dsm), lat_zmin(nlon_dsm,nlat_dsm))
  allocate(iboun(6,nspec)) ! boundary locator

  iboun(:,:) = .false.

  iaddx(1) = 0
  iaddy(1) = 0
  iaddz(1) = 0

  iaddx(2) = 1
  iaddy(2) = 0
  iaddz(2) = 0

  iaddx(3) = 1
  iaddy(3) = 1
  iaddz(3) = 0

  iaddx(4) = 0
  iaddy(4) = 1
  iaddz(4) = 0

  iaddx(5) = 0
  iaddy(5) = 0
  iaddz(5) = 1

  iaddx(6) = 1
  iaddy(6) = 0
  iaddz(6) = 1

  iaddx(7) = 1
  iaddy(7) = 1
  iaddz(7) = 1

  iaddx(8) = 0
  iaddy(8) = 1
  iaddz(8) = 1

!
!===========================================================================
!
!--- set up coordinates of the Gauss-Lobatto-Legendre points
!
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

!
!--- if number of points is odd, the middle abscissa is exactly zero
!
  if (mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
  if (mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
  if (mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

!
!--- get the 3-D shape functions
!
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)

!
!--- rotation matrix to switch to the geographical coordinates
!--- call euler_angles(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
!--- new rotation matrix
!

  call compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

!
!--- call ReadIasp91(vpv,vsv,density,zlayer,nlayer)
!

  call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

!
!--- calculation of the vertical discretization of layers
!

  Z_DEPTH_BLOCK = chunk_depth/1000.d0 !!!! switch to km

  call CalGridProf(ProfForGemini,current_layer,zlayer,nlayer,NZ,Z_DEPTH_BLOCK)

!
!===========================================================================
!
!--- GRID OF THE MESH
!
  izshift     = 0
  ispec       = 0
  kglob       = 0
  Ndepth      = 0
  ispec2Dxmin = 0
  ispec2Dxmax = 0
  ispec2Dymin = 0
  ispec2Dymax = 0
  ispec2Dzmin = 0
  ispec2Dzmax = 0

!
!--- Interface file  DSM-SPECFEM3D
!

  open(27, file = trim(MESH)//'.recdepth')  ! receptors on the vertical
  open(28, file = trim(MESH)//'stxmin')
  write(28,*) nlat_dsm          ! face xmin

  open(29, file = trim(MESH)//'stxmax')
  write(29,*) nlat_dsm ! face xmax

  open(30, file = trim(MESH)//'stymin')
  write(30,*) nlon_dsm ! face ymin

  open(31, file = trim(MESH)//'stymax')
  write(31,*) nlon_dsm ! face ymax

  open(38, file = trim(MESH)//'IgXmin')
  open(39, file = trim(MESH)//'IgXmax')
  open(40, file = trim(MESH)//'IgYmin')
  open(41, file = trim(MESH)//'IgYmax')
  open(42, file = trim(MESH)//'IgZmin')

! MESH for SPECFEM3D
  open(86, file = trim(MESH)//'nummaterial_velocity_file')
  open(87, file = trim(MESH)//'materials_file')

! open(88, file = 'model_1D.in')
  open(88, file = trim(MESH)//'OrigRepSpecfm')
  write(88,*) lon_center_chunk, lat_center_chunk
  write(88,*) chunk_azi, ANGULAR_WIDTH_XI_RAD/deg2rad, ANGULAR_WIDTH_ETA_RAD/deg2rad
  close(88)

  open(89, file = trim(MESH)//'flags_boundary.txt')
  open(90, file = trim(MESH)//'Nb_ielm_faces.txt')

!
!--- for new output mesh files and VM coupling with AxiSEM
!

  open(91, file = trim(MESH)//'list_ggl_boundary_spherical.txt')
  open(92, file = trim(MESH)//'list_ggl_boundary_Cartesian.txt')

! open(32,file='gll_zmin')
! open(125,file='ggl_elemts')

!
!-- Loop on the grid of the spectral elements
!

  ilayer    = 0
  index_mat = 0

  do iz = 0, nel_depth - 1

     ilayer_current=current_layer(iz) - 1 ! Caution between pickets and intervals !!

     if (iz /= 0) then
        if (current_layer(iz-1) /= current_layer(iz)) then
           izshift   = izshift + 1 ! point is repeated on the interface for DSM
           index_mat = index_mat - 1
           write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
                '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
        endif

     else
!       We write the first material
        index_mat = index_mat - 1
        write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
             '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
     endif

     do ilat=0,nel_lat-1
        do ilon=0,nel_lon-1

           ispec = ispec + 1
           ! material file
           write(87 ,*) ispec,index_mat

           istore_for_new_outputs = 0

           ! get boundary

           ! on boundary 1: x=xmin
           if (ilon == 0) then

              iboun(1,ispec)=.true.
              ispec2Dxmin=ispec2Dxmin+1
              write(89,*) ispec,ispec2Dxmin,1

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 2: xmax
           if (ilon == nel_lon-1) then

              iboun(2,ispec)=.true.
              ispec2Dxmax=ispec2Dxmax+1
              !write(*,*) '------ TOZ',ispec,ilon
              write(89,*) ispec,ispec2Dxmax,2

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 3: ymin
           if (ilat == 0) then

              iboun(3,ispec)=.true.
              ispec2Dymin=ispec2Dymin+1
              write(89,*) ispec,ispec2Dymin,3

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 4: ymax
           if (ilat == nel_lat-1) then

              iboun(4,ispec) =.true.
              ispec2Dymax=ispec2Dymax+1
              write(89,*) ispec,ispec2Dymax,4

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 5: bottom
           if (iz == 0) then

              iboun(5,ispec)=.true.
              ispec2Dzmin=ispec2Dzmin+1
              write(89,*) ispec,ispec2Dzmin,5

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 6: top
           if (iz == nel_depth-1) then
              ispec2Dzmax= ispec2Dzmax+1
              iboun(6,ispec)=.true.
           endif

          ! 8 vertices of the element ispec
           do ia=1,NGNOD

              i=iaddx(ia)
              j=iaddy(ia)
              k=iaddz(ia)

              z = 1000d0*ProfForGemini(iz,1+k)

              ! longitude
              ratio_xi = (dble(ilon+i)) / dble(NX)
              x = 2.d0*ratio_xi-1.d0
              x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

              ! latitude
              ratio_eta = (dble(ilat+j)) / dble(NY)
              y = 2.d0*ratio_eta-1.d0
              y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

              pz=  z/dsqrt(1.d0 + y*y + x*x) !(=r/s)
              px= pz * x !(tan(xi) * r/s)
              py= pz * y !(tan(eta) * r/s)

              ! old version
              xgrid(i+1,j+1,k+1,ispec) = px !px
              ygrid(i+1,j+1,k+1,ispec) = py !py
              zgrid(i+1,j+1,k+1,ispec) = pz

              xelm(ia)=xgrid(i+1,j+1,k+1,ispec)
              yelm(ia)=ygrid(i+1,j+1,k+1,ispec)
              zelm(ia)=zgrid(i+1,j+1,k+1,ispec)

           enddo

           ! INTERFACE FOR DSM ------

           ! Vertical receptors

           if (ilat == 0 .and. ilon == 0) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              call write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth,updown)
           endif

          ! Write two files giving Spherical coordinate on ALL the GLL points on the surface of the 3D chunk for the new DSM
          ! coupling (light version using 2D chunk)
          !
          ! CAUTION : will be also used later as INTERFACE for the VM coupling with AxiSEM
          !
          ! (must be after write_gllz_points to know the value of ilayer)

          if ( ( ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM .and. (.not. old_DSM_coupling_from_Vadim) ) .or. &
                 ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM                                        ) ) .and. &
               ( istore_for_new_outputs > 0 ) ) then


            call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
            call write_all_chunk_surface_GLL_in_spherical_and_Cartesian_coords(xstore,ystore,zstore, &
                                                          deg2rad,ilayer,iboun,ispec,nspec,longitud, &
                                                          latitud,radius,rotation_matrix,updown)

          endif

           ! Horizontal receptors

           ! stxmin
           if (ilon == 0 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilat == nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

           endif
            if (ilon == 0) call write_Igm_file(38,ispec2Dxmin,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stxmax
           if (ilon == nel_lon - 1 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
               if (ilat == nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilon == nel_lon-1)  call write_Igm_file(39,ispec2Dxmax,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stymin
           if (ilat == 0 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon == nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
               call write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat == 0) call write_Igm_file(40,ispec2Dymin,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)
           ! stymax
           if (ilat == nel_lat-1 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon == nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat == nel_lat-1) call write_Igm_file(41,ispec2Dymax,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)

           ! stzmin
           if (iz == 0) then ! pas besoin du test comme precedemment car je stocke tout dans des tableaux et c'est pas
                             ! grave si on recrit les memes choses
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
               call write_Igm_file(42,ispec2Dzmin,NGLLX,NGLLY,ilon,ilat,0,ilayer_current)
              call store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix, &
                   lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)
           endif

        enddo
     enddo
  enddo
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)

  ! ecriture des profondeurs de calcul pour DSM
  call write_recdepth_dsm(Ndepth,R_EARTH,MESH)
  ! ecriture de stzmin
  call write_stzmin(lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,MESH)
  !

  z_bottom = minval(zgrid(:,:,:,:))
  zgrid(:,:,:,:) = zgrid(:,:,:,:) - z_bottom
  UTM_X_MIN=minval(xgrid)
  UTM_X_MAX=maxval(xgrid)
 ! modele 1D
  open(88,file=trim(MESH)//'model_1D.in')
  write(88,*) nlayer,4
  do i=1,nlayer
     write(88,*) zlayer(i)
     write(88,'(4f20.10)') vpv(i,:)
     write(88,'(4f20.10)') vsv(i,:)
     write(88,'(4f20.10)') density(i,:)
  enddo
  write(88,*)  z_bottom
  write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
  close(88)

  !---------------- NUMEROTATION DES POINTS DE LA GRILLE ----

  ! on stocke tous les points de tous les elements
  do ispec=1,nspec

     ieoff = 8 * (ispec - 1)
     ilocnum = 0
     do k=1,2
        do j=1,2
           do i=1,2

              ilocnum = ilocnum + 1
              xp(ilocnum + ieoff)= xgrid(i,j,k,ispec)
              yp(ilocnum + ieoff)= ygrid(i,j,k,ispec)
              zp(ilocnum + ieoff)= zgrid(i,j,k,ispec)

           enddo
        enddo
     enddo
  enddo

  ! on identifie les points semblables et on les numerote
  call getglob_for_chunk(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD,UTM_X_MIN,UTM_X_MAX)

  deallocate(xp,yp,zp)
  allocate(xp(nglob),yp(nglob),zp(nglob))

  ! on ne stocke que les points de la grille et leur numeros
  do ispec=1,nspec
     ieoff = 8 * (ispec - 1)
     ilocnum = 0
     do k=1,2
        do j=1,2
           do i=1,2
              ilocnum=ilocnum+1
              inum_loc(i,j,k,ispec) = iglob(ilocnum+ieoff)
              xp(iglob(ilocnum+ieoff)) = xgrid(i,j,k,ispec)
              yp(iglob(ilocnum+ieoff)) = ygrid(i,j,k,ispec)
              zp(iglob(ilocnum+ieoff)) = zgrid(i,j,k,ispec)
           enddo
        enddo
     enddo
  enddo

!---------------------------------------------------------------------

  write(90,*)  ispec2Dxmin
  write(90,*)  ispec2Dxmax
  write(90,*)  ispec2Dymin
  write(90,*)  ispec2Dymax
  write(90,*)  ispec2Dzmin
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  close(32)
  close(37)
  close(38)
  close(39)
  close(40)
  close(41)
  close(42)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(86)
  close(87)
  close(88)
  close(89)
  close(90)
  close(91)
  close(92)

  ! -------------------------------- SAUVEGARDE DES MESH FILES -----------

  open(27,file=trim(MESH)//'nodes_coords_file')
  write(27,*) nglob ! nb de sommets
  do kglob=1,nglob
     write(27,'(i14,3x,3(f20.5,1x))') kglob,xp(kglob),yp(kglob),zp(kglob)
  enddo
  close(27)

  open(27,file=trim(MESH)//'mesh_file')
  write(27,*) nspec
  do ispec=1,nspec
     write(27,'(9i15)')  ispec,inum_loc(1,1,1,ispec),inum_loc(2,1,1,ispec), &
          inum_loc(2,2,1,ispec),inum_loc(1,2,1,ispec), &
          inum_loc(1,1,2,ispec),inum_loc(2,1,2,ispec), &
          inum_loc(2,2,2,ispec),inum_loc(1,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmin')
  write(27,*)  ispec2Dxmin
  do ispec=1,nspec
     if (iboun(1,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(1,2,1,ispec), &
          inum_loc(1,2,2,ispec),inum_loc(1,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmax')
  write(27,*) ispec2Dxmax
  do ispec=1,nspec
     if (iboun(2,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(2,1,1,ispec),inum_loc(2,2,1,ispec), &
          inum_loc(2,2,2,ispec),inum_loc(2,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymin')
  write(27,*) ispec2Dymin
  do ispec=1,nspec
     if (iboun(3,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(2,1,1,ispec), &
          inum_loc(2,1,2,ispec),inum_loc(1,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymax')
  write(27,*) ispec2Dymax
  do ispec=1,nspec
     if (iboun(4,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,2,1,ispec),inum_loc(2,2,1,ispec), &
          inum_loc(2,2,2,ispec),inum_loc(1,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_bottom')
  write(27,*) ispec2Dzmin
  do ispec=1,nspec
     if (iboun(5,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,1,ispec),inum_loc(1,2,1,ispec), &
          inum_loc(2,2,1,ispec),inum_loc(2,1,1,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'free_surface')
  write(27,*) ispec2Dzmax
  do ispec=1,nspec
     if (iboun(6,ispec)) write(27,'(5(i10,1x))') ispec,inum_loc(1,1,2,ispec),inum_loc(1,2,2,ispec), &
          inum_loc(2,2,2,ispec),inum_loc(2,1,2,ispec)
  enddo
  close(27)

  close(49)

  end subroutine earth_chunk_HEX8_Mesher

!=======================================================================================================

  subroutine earth_chunk_HEX27_Mesher(NGNOD)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM, R_EARTH, PI, ZERO, TINYVAL, &
    old_DSM_coupling_from_Vadim, INJECTION_TECHNIQUE_IS_AXISEM, INJECTION_TECHNIQUE_IS_DSM

  use shared_parameters, only: INJECTION_TECHNIQUE_TYPE

  implicit none

!==============================================================================================!
!                                                                                              !
!  Singular option of meshfem3D : MESH OF A GLOBE EARTH CHUNK FOR THE INTERFACE DSM-SPECFEM3D  !
!  Case of 27 nodes per element (HEX27)                                                        !
!                                                                                              !
!  Integrated in meshfem3d by CD, October 2014                                                 !
!                                                                                              !
!  WARNING : A local convention is used for the mapping of                                     !
!            the cubic sphere (to complete)                                                    !
!                                                                                              !
!==============================================================================================!

!
!--- Parameters
!

  integer, parameter :: myrank = 0
  integer, parameter :: nlayer = 12 !! (number of layer in the model iasp91, or ak135, or prem (one more layer than the model)

  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0

  logical, parameter ::  RUN_BENCHMARK = .false.

!
!--- Other
!

  integer NGNOD

  integer  nel_lat, nel_lon, nel_depth, NX, NY, NZ, Ndepth, nglob, kglob, ilocnum, ieoff, npointot
  integer ilat, ilon, ispec, iz, i, j, k, nspec, ia, izshift, index_mat
  integer ispec2Dxmin, ispec2Dxmax, ispec2Dymin, ispec2Dymax, ispec2Dzmin, ispec2Dzmax
  integer ilayer_current, ilayer
  integer nlat_dsm, nlon_dsm

  integer iaddx(NGNOD), iaddy(NGNOD), iaddz(NGNOD)

  integer, allocatable :: inum_loc(:,:,:,:), iglob(:), loc(:), current_layer(:)

  double precision ratio_eta, ratio_xi
  double precision ANGULAR_WIDTH_ETA_RAD, ANGULAR_WIDTH_XI_RAD, Z_DEPTH_BLOCK, UTM_X_MIN, UTM_X_MAX
  double precision lat_center_chunk, lon_center_chunk, chunk_depth, chunk_azi
  double precision deg2rad
  double precision x, y, z, px, py, pz, z_bottom

  double precision rotation_matrix(3,3)
  double precision zlayer(nlayer), vpv(nlayer,4), vsv(nlayer,4), density(nlayer,4)
  double precision xelm(NGNOD), yelm(NGNOD), zelm(NGNOD)
  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)

  !! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  !! GLL points and weights of integration
  double precision xigll(NGLLX), yigll(NGLLY), zigll(NGLLZ), wxgll(NGLLX), wygll(NGLLY), wzgll(NGLLZ)

  double precision, allocatable :: xp(:), yp(:), zp(:), xgrid(:,:,:,:), ygrid(:,:,:,:), zgrid(:,:,:,:)
  double precision, allocatable :: lon_zmin(:,:), lat_zmin(:,:)
  double precision, dimension(:,:), allocatable :: ProfForGemini

  integer ::  istore_for_new_outputs
  integer ::   updown(NGLLZ)
  double precision , dimension(NGLLX,NGLLY,NGLLZ) ::  longitud, latitud, radius

  logical test

  logical, allocatable :: ifseg(:)
  logical, dimension(:,:), allocatable :: iboun ! boundary locator

  character(len=100) line
  character(len=250) model1D_file

  character(len=10), parameter :: MESH = "./MESH/"

!
!--- WARNING, CONVENTION : (lon,lat) -> (xi,eta)
!---                       (k = 6 with -z for the mapping of the cubic sphere, cf Chevrot 2012)
!---                       We define the mesh of a chunk of the earth in the cubic sphere
!

  deg2rad = 3.141592653589793d0/180.d0

  open(49, file=trim(MESH)//'output_mesher_chunk_HEX27.txt')

  if (RUN_BENCHMARK) then
!
!--- Parameters fixed inside the subroutine for the moment
!
     ANGULAR_WIDTH_ETA_RAD = 10.d0 * deg2rad ! latitude 2.5
     ANGULAR_WIDTH_XI_RAD  = 20.d0 * deg2rad ! longitude 2.0

!    Chunk center
     lat_center_chunk      = 0.d0  ! 42.35d0 ! 42.5d0 !* deg2rad
     lon_center_chunk      = 60.d0 ! 1.3d0   ! 1.2d0  !* deg2rad

!    Azimuth
     chunk_azi             = 0.d0 !90.d0 !80.d0 !10.d0  !* deg2rad

!    Depth
     chunk_depth           = 1000.d0 * 1000.d0 ! 250.d0 * 1000.d0

!    Number of elements
     nel_lat               = 20 ! 120  ! 15
     nel_lon               = 40 ! 96   ! 15
     nel_depth             = 20 ! 100  ! 10

  else

     open(10, file=trim(MESH)//'ParFileMeshChunk')

     read(10,'(a)') line
     read(10,*) ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
     read(10,'(a)') line
     read(10,*) lon_center_chunk, lat_center_chunk, chunk_azi
     read(10,'(a)') line
     read(10,*) chunk_depth
     read(10,'(a)') line
     read(10,*) nel_lon,nel_lat, nel_depth
     read(10,'(a)') line
     read(10,'(a)') model1D_file

     model1D_file = 'MESH/'//trim(model1D_file)

     close(10)

     ANGULAR_WIDTH_XI_RAD  = deg2rad * ANGULAR_WIDTH_XI_RAD
     ANGULAR_WIDTH_ETA_RAD = deg2rad * ANGULAR_WIDTH_ETA_RAD
     chunk_depth           = chunk_depth * 1000.d0

  endif

  NX = nel_lon
  NY = nel_lat
  NZ = nel_depth

!
!===========================================================================
!
!--- TO DO : the reference chunk must be always symmetric (EW) and (NS)
!

  nlon_dsm = (ngllx - 1) * NX + 1
  nlat_dsm = (nglly - 1) * NY + 1
  nglob    = (2*nel_lat + 1) * (2*nel_lon + 1) * (2*nel_depth + 1)
  nspec    = nel_lat * nel_lon * nel_depth
  npointot = 27 * nspec

  allocate(xp(npointot), yp(npointot), zp(npointot))
  allocate(iglob(npointot), loc(npointot))
  allocate(ifseg(npointot))
  allocate(ProfForGemini(0:NZ-1,3))
  allocate(current_layer(0:NZ-1))
  allocate(inum_loc(3,3,3,nspec))
  allocate(xgrid(3,3,3,nspec), ygrid(3,3,3,nspec), zgrid(3,3,3,nspec))
  allocate(lon_zmin(nlon_dsm,nlat_dsm), lat_zmin(nlon_dsm,nlat_dsm))
  allocate(iboun(6,nspec)) ! boundary locator

  iboun(:,:) = .false.


!! MODIF HEX27 LA, cf call hex_nodes-----------------------------

 ! corner nodes

  iaddx(1) = 0
  iaddy(1) = 0
  iaddz(1) = 0

  iaddx(2) = 2
  iaddy(2) = 0
  iaddz(2) = 0

  iaddx(3) = 2
  iaddy(3) = 2
  iaddz(3) = 0

  iaddx(4) = 0
  iaddy(4) = 2
  iaddz(4) = 0

  iaddx(5) = 0
  iaddy(5) = 0
  iaddz(5) = 2

  iaddx(6) = 2
  iaddy(6) = 0
  iaddz(6) = 2

  iaddx(7) = 2
  iaddy(7) = 2
  iaddz(7) = 2

  iaddx(8) = 0
  iaddy(8) = 2
  iaddz(8) = 2

! midside nodes (nodes located in the middle of an edge)

  iaddx(9) = 1
  iaddy(9) = 0
  iaddz(9) = 0

  iaddx(10) = 2
  iaddy(10) = 1
  iaddz(10) = 0

  iaddx(11) = 1
  iaddy(11) = 2
  iaddz(11) = 0

  iaddx(12) = 0
  iaddy(12) = 1
  iaddz(12) = 0

  iaddx(13) = 0
  iaddy(13) = 0
  iaddz(13) = 1

  iaddx(14) = 2
  iaddy(14) = 0
  iaddz(14) = 1

  iaddx(15) = 2
  iaddy(15) = 2
  iaddz(15) = 1

  iaddx(16) = 0
  iaddy(16) = 2
  iaddz(16) = 1

  iaddx(17) = 1
  iaddy(17) = 0
  iaddz(17) = 2

  iaddx(18) = 2
  iaddy(18) = 1
  iaddz(18) = 2

  iaddx(19) = 1
  iaddy(19) = 2
  iaddz(19) = 2

  iaddx(20) = 0
  iaddy(20) = 1
  iaddz(20) = 2

! side center nodes (nodes located in the middle of a face)

  iaddx(21) = 1
  iaddy(21) = 1
  iaddz(21) = 0

  iaddx(22) = 1
  iaddy(22) = 0
  iaddz(22) = 1

  iaddx(23) = 2
  iaddy(23) = 1
  iaddz(23) = 1

  iaddx(24) = 1
  iaddy(24) = 2
  iaddz(24) = 1

  iaddx(25) = 0
  iaddy(25) = 1
  iaddz(25) = 1

  iaddx(26) = 1
  iaddy(26) = 1
  iaddz(26) = 2

! center node (barycenter of the eight corners)

  iaddx(27) = 1
  iaddy(27) = 1
  iaddz(27) = 1

!! --------------------------------------------
!
!===========================================================================
!
!--- set up coordinates of the Gauss-Lobatto-Legendre points
!
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

!
!--- if number of points is odd, the middle abscissa is exactly zero
!
  if (mod(NGLLX,2) /= 0) xigll((NGLLX - 1)/2 + 1) = ZERO
  if (mod(NGLLY,2) /= 0) yigll((NGLLY - 1)/2 + 1) = ZERO
  if (mod(NGLLZ,2) /= 0) zigll((NGLLZ - 1)/2 + 1) = ZERO

!
!--- get the 3-D shape functions
!
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll,NGNOD)

!
!--- rotation matrix to switch to the geographical coordinates
!--- call euler_angles(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)
!--- new rotation matrix
!

  call compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

!
!--- call ReadIasp91(vpv,vsv,density,zlayer,nlayer)
!

  call Read_dsm_model(model1D_file,vpv,vsv,density,zlayer,nlayer)

!
!--- calculation of the vertical discretization of layers
!

  Z_DEPTH_BLOCK = chunk_depth/1000.d0 !!!! switch to km

  call CalGridProf(ProfForGemini,current_layer,zlayer,nlayer,NZ,Z_DEPTH_BLOCK)

!
!===========================================================================
!
!--- GRID OF THE MESH
!
  izshift     = 0
  ispec       = 0
  kglob       = 0
  Ndepth      = 0
  ispec2Dxmin = 0
  ispec2Dxmax = 0
  ispec2Dymin = 0
  ispec2Dymax = 0
  ispec2Dzmin = 0
  ispec2Dzmax = 0

!
!--- Interface file  DSM-SPECFEM3D
!

  open(27, file = trim(MESH)//'.recdepth')  ! receptors on the vertical
  open(28, file = trim(MESH)//'stxmin')
  write(28,*) nlat_dsm          ! face xmin

  open(29, file = trim(MESH)//'stxmax')
  write(29,*) nlat_dsm ! face xmax

  open(30, file = trim(MESH)//'stymin')
  write(30,*) nlon_dsm ! face ymin

  open(31, file = trim(MESH)//'stymax')
  write(31,*) nlon_dsm ! face ymax

  open(38, file = trim(MESH)//'IgXmin')
  open(39, file = trim(MESH)//'IgXmax')
  open(40, file = trim(MESH)//'IgYmin')
  open(41, file = trim(MESH)//'IgYmax')
  open(42, file = trim(MESH)//'IgZmin')

! MESH for SPECFEM3D
  open(86, file = trim(MESH)//'nummaterial_velocity_file')
  open(87, file = trim(MESH)//'materials_file')

! open(88, file = 'model_1D.in')
  open(88, file = trim(MESH)//'OrigRepSpecfm')
  write(88,*) lon_center_chunk, lat_center_chunk
  write(88,*) chunk_azi, ANGULAR_WIDTH_XI_RAD/deg2rad, ANGULAR_WIDTH_ETA_RAD/deg2rad
  close(88)

  open(89, file = trim(MESH)//'flags_boundary.txt')
  open(90, file = trim(MESH)//'Nb_ielm_faces.txt')

!
!--- for new output mesh files and VM coupling with AxiSEM
!

  open(91, file = trim(MESH)//'list_ggl_boundary_spherical.txt')
  open(92, file = trim(MESH)//'list_ggl_boundary_Cartesian.txt')

!
!-- Loop on the grid of the spectral elements
!

  ilayer    = 0
  index_mat = 0

  do iz = 0, nel_depth - 1

     ilayer_current=current_layer(iz) - 1 ! Caution between pickets and intervals !!

     if (iz /= 0) then
        if (current_layer(iz-1) /= current_layer(iz)) then
           izshift   = izshift + 1 ! point is repeated on the interface for DSM
           index_mat = index_mat - 1
           write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
                '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
        endif

     else
!       We write the first material
        index_mat = index_mat - 1
        write(86,'(a1,2x,i10,2x,a10,2x,a7,2x,a20,2x,a1)') &
             '2', index_mat, 'tomography', 'elastic', 'tomography_model.xyz', '1'
     endif

     do ilat=0,nel_lat-1
        do ilon=0,nel_lon-1

           ispec = ispec + 1
           ! material file
           write(87 ,*) ispec,index_mat

           istore_for_new_outputs = 0

           ! get boundary

           ! on boundary 1: x=xmin
           if (ilon == 0) then

              iboun(1,ispec)=.true.
              ispec2Dxmin=ispec2Dxmin+1
              write(89,*) ispec,ispec2Dxmin,1

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 2: xmax
           if (ilon == nel_lon-1) then

              iboun(2,ispec)=.true.
              ispec2Dxmax=ispec2Dxmax+1
              !write(*,*) '------ TOZ',ispec,ilon
              write(89,*) ispec,ispec2Dxmax,2

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 3: ymin
           if (ilat == 0) then

              iboun(3,ispec)=.true.
              ispec2Dymin=ispec2Dymin+1
              write(89,*) ispec,ispec2Dymin,3

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 4: ymax
           if (ilat == nel_lat-1) then

              iboun(4,ispec) =.true.
              ispec2Dymax=ispec2Dymax+1
              write(89,*) ispec,ispec2Dymax,4

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 5: bottom
           if (iz == 0) then

              iboun(5,ispec)=.true.
              ispec2Dzmin=ispec2Dzmin+1
              write(89,*) ispec,ispec2Dzmin,5

              istore_for_new_outputs = istore_for_new_outputs + 1

           endif

           ! on boundary 6: top
           if (iz == nel_depth-1) then
              ispec2Dzmax= ispec2Dzmax+1
              iboun(6,ispec)=.true.
           endif

           do ia=1,NGNOD

!! MODIF HEX27 LA -----------------------------

              i=iaddx(ia)
              j=iaddy(ia)
              k=iaddz(ia)

              SELECT CASE (k)
                CASE(0)
                  z = 1000d0*ProfForGemini(iz,1)
                CASE(1)
                  z = 1000d0*ProfForGemini(iz,3)
                CASE(2)
                  z = 1000d0*ProfForGemini(iz,2)
              END SELECT

              ! longitude
              ratio_xi = (dble(ilon) + dble(i)/2.d0 ) / dble(NX)
              x = 2.d0*ratio_xi-1.d0
              x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

              ! latitude
              ratio_eta = (dble(ilat) + dble(j)/2.d0 ) / dble(NY)
              y = 2.d0*ratio_eta-1.d0
              y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

              ! mapping cubic sphere (k=6, Chevrot at al 2012, avec -z)

              pz= z/dsqrt(1.d0 + y*y + x*x) !(=r/s)
              px= pz * x !(tan(xi) * r/s)
              py= pz * y !(tan(eta) * r/s)

              ! old version
              xgrid(i+1,j+1,k+1,ispec) = px !px
              ygrid(i+1,j+1,k+1,ispec) = py !py
              zgrid(i+1,j+1,k+1,ispec) = pz

              xelm(ia)=xgrid(i+1,j+1,k+1,ispec)
              yelm(ia)=ygrid(i+1,j+1,k+1,ispec)
              zelm(ia)=zgrid(i+1,j+1,k+1,ispec)

           enddo

           ! INTERFACE FOR DSM ------

           ! Vertical receptors

           if (ilat == 0 .and. ilon == 0) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              call write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth,updown)
           endif

          ! Write two files giving Spherical coordinate on ALL the GLL points on the surface of the 3D chunk for the new DSM
          ! coupling (light version using 2D chunk)
          !
          ! CAUTION : will be also used later for the VM coupling with AxiSEM
          !
          ! (must be after write_gllz_points to know the value of ilayer)

          if ( ( ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM .and. (.not. old_DSM_coupling_from_Vadim) ) .or. &
                 ( INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM                                        ) ) .and. &
               ( istore_for_new_outputs > 0 ) ) then


            call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
            call write_all_chunk_surface_GLL_in_spherical_and_Cartesian_coords(xstore,ystore,zstore, &
                                                          deg2rad,ilayer,iboun,ispec,nspec,longitud, &
                                                          latitud,radius,rotation_matrix,updown)

          endif

           ! Horizontal receptors

           ! stxmin
           if (ilon == 0 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilat == nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

           endif
            if (ilon == 0) call write_Igm_file(38,ispec2Dxmin,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stxmax
           if (ilon == nel_lon - 1 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
               if (ilat == nel_lat-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call  write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilon == nel_lon-1)  call write_Igm_file(39,ispec2Dxmax,NGLLY,NGLLZ,ilat,iz,izshift,ilayer_current)

           ! stymin
           if (ilat == 0 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon == nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
               call write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat == 0) call write_Igm_file(40,ispec2Dymin,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)
           ! stymax
           if (ilat == nel_lat-1 .and. iz == nel_depth-1) then
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)
              if (ilon == nel_lon-1) then ! This test is for add the last GLL point
                 test=.true.
              else
                 test=.false.
              endif
              call write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)
           endif
           if (ilat == nel_lat-1) call write_Igm_file(41,ispec2Dymax,NGLLX,NGLLZ,ilon,iz,izshift,ilayer_current)

           ! stzmin
           if (iz == 0) then ! pas besoin du test comme precedemment car je stocke tout dans des tableaux et c'est pas
                             ! grave si on recrit les memes choses
              call calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)

              call write_Igm_file(42,ispec2Dzmin,NGLLX,NGLLY,ilon,ilat,0,ilayer_current)

              call store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix, &
                   lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)
           endif

        enddo
     enddo
  enddo
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)

  ! ecriture des profondeurs de calcul pour DSM
  call write_recdepth_dsm(Ndepth,R_EARTH,MESH)
  ! ecriture de stzmin
  call write_stzmin(lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,MESH)
  !

  z_bottom = minval(zgrid(:,:,:,:))
  zgrid(:,:,:,:) = zgrid(:,:,:,:) - z_bottom
  UTM_X_MIN=minval(xgrid)
  UTM_X_MAX=maxval(xgrid)
 ! modele 1D
  open(88,file=trim(MESH)//'model_1D.in')
  write(88,*) nlayer,4
  do i=1,nlayer
     write(88,*) zlayer(i)
     write(88,'(4f20.10)') vpv(i,:)
     write(88,'(4f20.10)') vsv(i,:)
     write(88,'(4f20.10)') density(i,:)
  enddo
  write(88,*)  z_bottom
  write(88,*)  lon_center_chunk,  lat_center_chunk,  chunk_azi
  close(88)

  !---------------- NUMEROTATION DES POINTS DE LA GRILLE ----

  ! on stocke tous les points de tous les elements
  do ispec=1,nspec

     ieoff = 27 * (ispec - 1)
     ilocnum = 0

     do k=1,3
        do j=1,3
           do i=1,3

              ilocnum = ilocnum + 1
              xp(ilocnum + ieoff)= xgrid(i,j,k,ispec)
              yp(ilocnum + ieoff)= ygrid(i,j,k,ispec)
              zp(ilocnum + ieoff)= zgrid(i,j,k,ispec)

           enddo
        enddo
     enddo
  enddo

  ! on identifie les points semblables et on les numerote
  call getglob_for_chunk(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD,UTM_X_MIN,UTM_X_MAX)

  deallocate(xp,yp,zp)
  allocate(xp(nglob),yp(nglob),zp(nglob))

!! MODIF HEX27 LA -----------------------------

  ! on ne stocke que les points de la grille et leur numeros
  do ispec=1,nspec

     ieoff = 27 * (ispec - 1)
     ilocnum = 0

     do k=1,3
        do j=1,3
           do i=1,3

              ilocnum                  = ilocnum + 1
              inum_loc(i,j,k,ispec)    = iglob(ilocnum+ieoff)
              xp(iglob(ilocnum+ieoff)) = xgrid(i,j,k,ispec)
              yp(iglob(ilocnum+ieoff)) = ygrid(i,j,k,ispec)
              zp(iglob(ilocnum+ieoff)) = zgrid(i,j,k,ispec)

           enddo
        enddo
     enddo
  enddo

!---------------------------------------------------------------------

  write(90,*)  ispec2Dxmin
  write(90,*)  ispec2Dxmax
  write(90,*)  ispec2Dymin
  write(90,*)  ispec2Dymax
  write(90,*)  ispec2Dzmin
  close(27)
  close(28)
  close(29)
  close(30)
  close(31)
  close(32)
  close(37)
  close(38)
  close(39)
  close(40)
  close(41)
  close(42)
  close(81)
  close(82)
  close(83)
  close(84)
  close(85)
  close(86)
  close(87)
  close(88)
  close(89)
  close(90)
  close(91)
  close(92)

  ! -------------------------------- SAUVEGARDE DES MESH FILES -----------

!! MODIF HEX27 LA -----------------------------

  open(27,file=trim(MESH)//'nodes_coords_file')
  write(27,*) nglob ! nb de sommets
  do kglob=1,nglob
     write(27,'(i14,3x,3(f20.5,1x))') kglob,xp(kglob),yp(kglob),zp(kglob)
  enddo
  close(27)

  open(27,file=trim(MESH)//'mesh_file')
  write(27,*) nspec
  do ispec=1,nspec
    write(27,'(28i15)') ispec, &
                        inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), inum_loc(3,3,1,ispec), inum_loc(1,3,1,ispec), &
                        inum_loc(1,1,3,ispec), inum_loc(3,1,3,ispec), inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                        inum_loc(2,1,1,ispec), inum_loc(3,2,1,ispec), inum_loc(2,3,1,ispec), inum_loc(1,2,1,ispec), &
                        inum_loc(1,1,2,ispec), inum_loc(3,1,2,ispec), inum_loc(3,3,2,ispec), inum_loc(1,3,2,ispec), &
                        inum_loc(2,1,3,ispec), inum_loc(3,2,3,ispec), inum_loc(2,3,3,ispec), inum_loc(1,2,3,ispec), &
                        inum_loc(2,2,1,ispec), inum_loc(2,1,2,ispec), inum_loc(3,2,2,ispec), inum_loc(2,3,2,ispec), &
                        inum_loc(1,2,2,ispec), inum_loc(2,2,3,ispec), inum_loc(2,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmin')
  write(27,*)  ispec2Dxmin
  do ispec=1,nspec
     if (iboun(1,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(1,3,1,ispec), &
                                                  inum_loc(1,3,3,ispec), inum_loc(1,1,3,ispec), &
                                                  inum_loc(1,2,1,ispec), inum_loc(1,3,2,ispec), &
                                                  inum_loc(1,2,3,ispec), inum_loc(1,1,2,ispec), &
                                                  inum_loc(1,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_xmax')
  write(27,*) ispec2Dxmax
  do ispec=1,nspec
     if (iboun(2,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(3,1,1,ispec), inum_loc(3,3,1,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(3,1,3,ispec), &
                                                  inum_loc(3,2,1,ispec), inum_loc(3,3,2,ispec), &
                                                  inum_loc(3,2,3,ispec), inum_loc(3,1,2,ispec), &
                                                  inum_loc(3,2,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymin')
  write(27,*) ispec2Dymin
  do ispec=1,nspec
     if (iboun(3,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), &
                                                  inum_loc(3,1,3,ispec), inum_loc(1,1,3,ispec), &
                                                  inum_loc(2,1,1,ispec), inum_loc(3,1,2,ispec), &
                                                  inum_loc(2,1,3,ispec), inum_loc(1,1,2,ispec), &
                                                  inum_loc(2,1,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_ymax')
  write(27,*) ispec2Dymax
  do ispec=1,nspec
     if (iboun(4,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,3,1,ispec), inum_loc(3,3,1,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                                                  inum_loc(2,3,1,ispec), inum_loc(3,3,2,ispec), &
                                                  inum_loc(2,3,3,ispec), inum_loc(1,3,2,ispec), &
                                                  inum_loc(2,3,2,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'absorbing_surface_file_bottom')
  write(27,*) ispec2Dzmin
  do ispec=1,nspec
     if (iboun(5,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,1,ispec), inum_loc(3,1,1,ispec), &
                                                  inum_loc(3,3,1,ispec), inum_loc(1,3,1,ispec), &
                                                  inum_loc(2,1,1,ispec), inum_loc(3,2,1,ispec), &
                                                  inum_loc(2,3,1,ispec), inum_loc(1,2,1,ispec), &
                                                  inum_loc(2,2,1,ispec)
  enddo
  close(27)

  open(27,file=trim(MESH)//'free_surface')
  write(27,*) ispec2Dzmax
  do ispec=1,nspec
     if (iboun(6,ispec)) write(27,'(10(i10,1x))') ispec, &
                                                  inum_loc(1,1,3,ispec), inum_loc(3,1,3,ispec), &
                                                  inum_loc(3,3,3,ispec), inum_loc(1,3,3,ispec), &
                                                  inum_loc(2,1,3,ispec), inum_loc(3,2,3,ispec), &
                                                  inum_loc(2,3,3,ispec), inum_loc(1,2,3,ispec), &
                                                  inum_loc(2,2,3,ispec)
  enddo
  close(27)

  close(49)

  end subroutine earth_chunk_HEX27_Mesher

!=======================================================================================================!

  subroutine earth_chunk_ReadIasp91(vp,vs,rho,rb,n)

  implicit none

  integer i,j,n,iunit,nlay,nco(n),ifanis
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n)
  real fref,vph,vsh,qm,qk,eta

  character(len=80) text
  character(len=2) cnlay
  character(len=11) format_to_use

  do i=1,n
     !qm(i)=0.d0
        !qk(i)=0.d0
     rb(i)=0.d0
     !iflso(i)=0
     nco(i)=0
     do j=1,4
        rho(i,j)=0.d0
        vp(i,j)=0.d0
        !vph(i,j)=0.d0
        vs(i,j)=0.d0
        !vsh(i,j)=0.d0
        !eta(i,j)=0.d0
     enddo
  enddo
  iunit=26
  open(unit=iunit,file='iasp91',status='old')

1 read(iunit,'(a72)') text
  if (text(1:1) == '#') then
     goto 1
  endif
  backspace iunit


  read(iunit,'(i2)') nlay                ! Number of layers

  write(cnlay,'(i2)') nlay
  format_to_use='('//cnlay//'i2)'                 ! Number of polynomial
  read(iunit,format_to_use) (nco(i),i=1,nlay)     ! coefficients for each layer

  read(iunit,*) fref               ! reference frequency of Qs in Hertz
  read(iunit,*) ifanis             ! Transversal isotropic? 1=y, else=n
  read(iunit,'(1x/1x/)')


  do i = 1, nlay

     read(iunit,*) rb(i),rho(i,1),vp(i,1),vph,vs(i,1),vsh,qm,qk,eta
     do j = 2, nco(i)
        read(iunit,*) rho(i,j),vp(i,j),vph,vs(i,j),vsh,eta
     enddo
     read(iunit,'(1x)')
  enddo
  i = nlay+1
  read(iunit,*) rb(i)
  j = 1
  rho(i,j) =  rho(i-1,j)
  vp(i,j) = vp(i-1,j)
  vs(i,j) = vs(i-1,j)

  end subroutine earth_chunk_ReadIasp91

!
!===========================================================================
!

  subroutine Read_dsm_model(model_file,vp,vs,rho,rb,n)

  implicit none

  integer i,n,iunit,nco(n)
  double precision vp(n,4),vs(n,4),rho(n,4),rb(n),eta(4),vrmin,vrmax
  real vph(4),vsh(4),qm,qk
  integer nzone

  character(len=250) model_file

  rb    = 0.d0
  rho   = 0.d0
  vp    = 0.d0
  vs    = 0.d0
  nco   = 0
  iunit = 26

  open(unit=iunit,file=trim(model_file),status='old',action='read')

  read(iunit,*) nzone

  do i=1, nzone
     read(iunit,*) vrmin, vrmax, &
          rho(i,1), rho(i,2), rho(i,3), rho(i,4), &
          vp(i,1), vp(i,2), vp(i,3), vp(i,4), &
          vph(1), vph(2), vph(3), vph(4), &
          vs(i,1), vs(i,2), vs(i,3), vs(i,4), &
          vsh(1), vsh(2), vsh(3), vsh(4), &
          eta(1), eta(2), eta(3), eta(4), &
          qm, qk
          rb(i)=vrmin
  enddo

  i        = nzone+1
  rb(i)    = vrmax
  vp(i,:)  = vp(i-1,:)
  vs(i,:)  = vs(i-1,:)
  rho(i,:) = rho(i-1,:)

  close(iunit)

  end subroutine Read_dsm_model

!
!=======================================================================================================
!
! compute the Euler angles and the associated rotation matrix

  subroutine euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  implicit none

  double precision rotation_matrix(3,3)
  double precision CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  double precision alpha,beta,gamma
  double precision sina,cosa,sinb,cosb,sing,cosg

  double precision DEGREES_TO_RADIANS

  DEGREES_TO_RADIANS = 3.141592653589793d0/180.d0


! compute colatitude and longitude and convert to radians
  alpha = CENTER_LONGITUDE_IN_DEGREES * DEGREES_TO_RADIANS
  beta = (90.0d0 - CENTER_LATITUDE_IN_DEGREES) * DEGREES_TO_RADIANS
  gamma = GAMMA_ROTATION_AZIMUTH * DEGREES_TO_RADIANS

  sina = dsin(alpha)
  cosa = dcos(alpha)
  sinb = dsin(beta)
  cosb = dcos(beta)
  sing = dsin(gamma)
  cosg = dcos(gamma)

! define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
  rotation_matrix(1,3) = sinb*cosa
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
  rotation_matrix(2,3) = sinb*sina
  rotation_matrix(3,1) = -cosg*sinb
  rotation_matrix(3,2) = sing*sinb
  rotation_matrix(3,3) = cosb

  end subroutine euler_angles

!=======================================================================================================

  subroutine write_gllz_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,current_layer,nel_depth,ilayer,iz,Ndepth,updown)

  implicit none

  integer NGLLX,NGLLY,NGLLZ,nel_depth,iz,Ndepth
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision profondeur
  integer current_layer(0:nel_depth-1),ilayer,k
  integer updown(NGLLZ) !! will be also used for VM coupling with AxiSEM

  updown(:) = 0
  if (ilayer == current_layer(iz)) then

    do k=2,NGLLZ
      profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
      write(27,*) profondeur/1000., ilayer-1,1
      Ndepth = Ndepth + 1
      updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM
    enddo

  else ! new layer

     k=1
     profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
     if (ilayer == 0) then
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1

        updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM

     else
        ilayer =  current_layer(iz)
        write(27,*)  profondeur/1000., ilayer-1,-1
        Ndepth=Ndepth+1

        updown(k) = -1 !! for new output mesh files and VM coupling with AxiSEM

     endif
     do k=2,NGLLZ ! on duplique le dernier point
        profondeur = dsqrt(xstore(1,1,k)**2 + ystore(1,1,k)**2 + (zstore(1,1,k) )**2 )
        write(27,*)  profondeur/1000., ilayer-1,1
        Ndepth=Ndepth+1

        updown(k) = 0 !! for new output mesh files and VM coupling with AxiSEM

     enddo


  endif

  end subroutine write_gllz_points


!=======================================================================================================
!
!=======================================================================================================
!
! Useless for the moment, we will maybe need it later
!
!!$! To have one and only file, who give the Spherical coordinate on ALL the GLL points
!!$! on the surface of the 3D chunk, for the new DSM coupling (light version using 2D chunk)
!!$!
!!$
!!$  subroutine Cartesian_product_to_r_theta_phi_on_chunk_surface_GLL(MESH,deg2rad)
!!$
!!$  use constants, only: R_EARTH_KM
!!$
!!$  implicit none
!!$
!!$  character(len=10)  :: MESH
!!$  double precision   :: deg2rad
!!$  integer            :: np_r, np_xmin, np_xmax, np_ymin, np_ymax, np_zmin, recflag1, recflag2, i, j, np_surf, ios
!!$  double precision   :: rec_val, xmin_val1, xmin_val2, xmax_val1, xmax_val2, ymin_val1, ymin_val2
!!$  double precision   :: ymax_val1, ymax_val2, zmin_val1, zmin_val2, zmin_fix, x, y ,z, R, R_m, latrad, lgrad
!!$
!!$  open(unit=10,file=trim(MESH)//'recdepth',action='read',status='unknown',iostat=ios)
!!$  open(unit=11,file=trim(MESH)//'stxmin',action='read',status='unknown',iostat=ios)
!!$  open(unit=12,file=trim(MESH)//'stxmax',action='read',status='unknown',iostat=ios)
!!$  open(unit=13,file=trim(MESH)//'stymin',action='read',status='unknown',iostat=ios)
!!$  open(unit=14,file=trim(MESH)//'stymax',action='read',status='unknown',iostat=ios)
!!$  open(unit=15,file=trim(MESH)//'stzmin',action='read',status='unknown',iostat=ios)
!!$
!!$  open(unit=20,file=trim(MESH)//'chunk_surface_GLL_r_theta_phi.out',status='unknown',iostat=ios)
!!$
!!$  read(10,*) np_r
!!$
!!$  read(11,*) np_xmin
!!$  read(12,*) np_xmax
!!$  read(13,*) np_ymin
!!$  read(14,*) np_ymax
!!$  read(15,*) np_zmin
!!$
!!$  np_surf = np_r*(np_xmin + np_xmax + np_ymin + np_ymax) + np_zmin
!!$
!!$  write(20,*) np_surf
!!$
!!$  do i=1,np_r
!!$
!!$    rewind(11)
!!$    read(11,*)
!!$    rewind(12)
!!$    read(12,*)
!!$    rewind(13)
!!$    read(13,*)
!!$    rewind(14)
!!$    read(14,*)
!!$
!!$    read(10,*) rec_val, recflag1, recflag2
!!$
!!$    R = dabs(R_EARTH_KM - rec_val) !! kM

!!$    do j=1,np_xmin
!!$
!!$      read(11,*) xmin_val1, xmin_val2
!!$      write(20,*) R, xmin_val1, xmin_val2
!!$
!!$ enddo
!!$
!!$    do j=1,np_xmax
!!$
!!$      read(12,*) xmax_val1, xmax_val2
!!$      write(20,*) R, xmax_val1, xmax_val2
!!$
!!$    enddo
!!$
!!$    do j=1,np_ymin
!!$
!!$      read(13,*) ymin_val1, ymin_val2
!!$      write(20,*) R, ymin_val1, ymin_val2
!!$
!!$    enddo
!!$
!!$    do j=1,np_ymax
!!$
!!$      read(14,*) ymax_val1, ymax_val2
!!$      write(20,*) R, ymax_val1, ymax_val2
!!$
!!$ enddo
!!$
!!$ if (i == np_r) zmin_fix = rec_val !! maximal depth
!!$
!!$  enddo
!!$
!!$  rewind(15)v
!!$  read(15,*)
!!$
!!$  R = dabs(R_EARTH_KM - zmin_fix) !! kM
!!$
!!$  do j=1,np_zmin
!!$
!!$    read(15,*) zmin_val1, zmin_val2
!!$    write(20,*) R, zmin_val1, zmin_val2
!!$
!!$  enddo
!!$
!!$  close(10)
!!$  close(11)
!!$  close(12)
!!$  close(13)
!!$  close(14)
!!$  close(15)
!!$  close(20)
!!$
!!$  end subroutine Cartesian_product_to_r_theta_phi_on_chunk_surface_GLL
!
!=======================================================================================================
!
!=======================================================================================================
!
!! Used (among other) for VM coupling with AxiSEM

  subroutine write_all_chunk_surface_GLL_in_spherical_and_Cartesian_coords(xstore,ystore,zstore, &
                                                      deg2rad,ilayer,iboun,ispec,nspec,longitud, &
                                                      latitud,radius,rotation_matrix,updown)

  use constants, only: NGLLX, NGLLY, NGLLZ, &
    INJECTION_TECHNIQUE_IS_DSM, INJECTION_TECHNIQUE_IS_AXISEM

  use shared_parameters, only: INJECTION_TECHNIQUE_TYPE

  implicit none

  integer ispec, nspec, ilayer
  integer i, j, k, imin, imax, jmin, jmax, kmin, kmax
  integer updown(NGLLZ)
  logical :: iboun(6,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: longitud, latitud, radius
  double precision rotation_matrix(3,3)
  double precision deg2rad

!
!---- CF 'earth_chunk_HEX8_Mesher' and 'earth_chunk_HEX27_Mesher' to see the name of files whose units are 91 and 92
!

1000 format(3f30.10)

!-- all GLL points in geographical coordinates

  call Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)

!
!-- xmin ----
!

  if (iboun(1,ispec)) then

    imin = 1
    imax = 1
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          ! CF 'earth_chunk_HEX8_Mesher' and 'earth_chunk_HEX27_Mesher' to see files whose units are 91 and 92

          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,1,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- xmax ----
!

  if (iboun(2,ispec)) then

    imin = NGLLX
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,2,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- ymin ----
!

  if (iboun(3,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = 1
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,3,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- ymax ----
!

  if (iboun(4,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = NGLLY
    jmax = NGLLY
    kmin = 1
    kmax = NGLLZ

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,4,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

!
!-- zmin ----
!

  if (iboun(5,ispec)) then

    imin = 1
    imax = NGLLX
    jmin = 1
    jmax = NGLLY
    kmin = 1
    kmax = 1

    do k=kmin,kmax
      do j=jmin,jmax
        do i=imin,imax

          if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) then
            write(92,1000) xstore(i,j,k), ystore(i,j,k), zstore(i,j,k)

          else if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) then
            write(92,'(3f25.10,i10,6i3)') xstore(i,j,k),ystore(i,j,k),zstore(i,j,k),ispec,i,j,k,5,ilayer,updown(k)

          endif

          write(91,1000) radius(i,j,k), latitud(i,j,k), longitud(i,j,k)

        enddo
      enddo
    enddo

  endif

  end subroutine write_all_chunk_surface_GLL_in_spherical_and_Cartesian_coords

!
!=======================================================================================================
!

  subroutine Cartesian2spheric(xstore,ystore,zstore,rotation_matrix,longitud,latitud,radius,deg2rad)

  use constants, only: NGLLX, NGLLY, NGLLZ, NDIM

  implicit none

  integer i, j, igll, jgll, kgll

  double precision xstore(NGLLX,NGLLY,NGLLZ), ystore(NGLLX,NGLLY,NGLLZ), zstore(NGLLX,NGLLY,NGLLZ)
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: longitud, latitud, radius
  double precision rotation_matrix(3,3)
  double precision vector_ori(3), vector_rotated(3)
  double precision rayon, x, y, z, long, lati, deg2rad

!
!----
!

  do kgll=1,NGLLZ
    do jgll=1,NGLLY
      do igll=1,NGLLX

        vector_ori(1) = xstore(igll,jgll,kgll)
        vector_ori(2) = ystore(igll,jgll,kgll)
        vector_ori(3) = zstore(igll,jgll,kgll)

        do i = 1,NDIM

          vector_rotated(i) = 0.d0

          do j = 1,NDIM

            vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)

          enddo
        enddo

        x     = vector_rotated(1)
        y     = vector_rotated(2)
        z     = vector_rotated(3)
        rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

        long  = datan2(y,x)
        lati  = dasin(z/rayon)

        longitud(igll,jgll,kgll) = long/deg2rad
        latitud(igll,jgll,kgll)  = lati/deg2rad
        radius(igll,jgll,kgll)   = rayon/1000.d0

      enddo
    enddo
  enddo

  end subroutine Cartesian2spheric


!=======================================================================================================
!
!=======================================================================================================

  subroutine write_recdepth_dsm(Ndepth,R_EARTH,MESH)

  implicit none

  integer Ndepth,i
  double precision R_EARTH,prof
  double precision, allocatable :: z(:)
  integer, allocatable :: zindex(:),ziflag(:)
  integer ilayer,flag
  character(len=10) MESH

  open(27,file=trim(MESH)//'.recdepth')
  allocate(zindex(Ndepth),ziflag(Ndepth))
  allocate(z(Ndepth))

  do i=1,Ndepth
     read(27,*) prof,ilayer,flag
     z(Ndepth-i+1)=R_EARTH/1000.d0-prof
     zindex(Ndepth-i+1)=ilayer
     ziflag(Ndepth-i+1)=flag
  enddo
  close(27)

  open(27,file=trim(MESH)//'recdepth')
  write(27,*) Ndepth
  i=1
  write(27,*) z(i),zindex(i),ziflag(i)
  do i=2,Ndepth-1
     if (ziflag(i-1) == -1 ) then
        write(27,*) z(i),zindex(i),-1
     else
         write(27,*) z(i),zindex(i),1
     endif
  enddo
  i=Ndepth
  write(27,*) z(i),zindex(i),ziflag(i)
  close(27)

end subroutine write_recdepth_dsm

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(1,jgll,NGLLZ)
     vector_ori(2)=ystore(1,jgll,NGLLZ)
     vector_ori(3)=zstore(1,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      write(28,*) long/deg2rad,lati/deg2rad !,rayon/1000

  enddo

  end subroutine write_stxmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stxmax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLY_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLY_eff = NGLLY
  else
     NGLLY_eff = NGLLY - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY_eff
     vector_ori(1)=xstore(NGLLX,jgll,NGLLZ)
     vector_ori(2)=ystore(NGLLX,jgll,NGLLZ)
     vector_ori(3) =zstore(NGLLX,jgll,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      write(29,*) long/deg2rad,lati/deg2rad !,rayon/1000

  enddo

  end subroutine write_stxmax

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymin(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

   if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,1,NGLLZ)
     vector_ori(2)=ystore(jgll,1,NGLLZ)
     vector_ori(3) =zstore(jgll,1,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      write(30,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stymax(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix,test)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,jgll,i,j,NGLLX_eff
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  logical test

  if (test) then
     NGLLX_eff = NGLLX
  else
     NGLLX_eff = NGLLX - 1
  endif

  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLX_eff
     vector_ori(1)=xstore(jgll,NGLLY,NGLLZ)
     vector_ori(2)=ystore(jgll,NGLLY,NGLLZ)
     vector_ori(3) =zstore(jgll,NGLLY,NGLLZ)

     do i = 1,NDIM
        vector_rotated(i) = 0.d0
        do j = 1,NDIM
           vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
        enddo
     enddo
     x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
     rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

      long=atan2(y,x)
      lati=asin(z/rayon)

      write(31,*) long/deg2rad,lati/deg2rad !,rayon/1000
  enddo

  end subroutine write_stymax

!=======================================================================================================
!
!=======================================================================================================

  subroutine store_zmin_points(xstore,ystore,zstore,NGLLX,NGLLY,NGLLZ,rotation_matrix, &
             lon_zmin,lat_zmin,nlon_dsm,nlat_dsm,ilon,ilat)

  implicit none

  integer NDIM,NGLLX,NGLLY,NGLLZ,igll,jgll,i,j
  integer ilon,ilat,iglob,jglob,nlat_dsm,nlon_dsm
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision rayon,x,y,z,deg2rad,long,lati
  double precision lon_zmin(nlon_dsm,nlat_dsm),lat_zmin(nlon_dsm,nlat_dsm)


  deg2rad=3.141592653589793d0/180.d0
  NDIM=3

  do jgll=1,NGLLY
     do igll=1,NGLLX
        vector_ori(1)=xstore(igll,jgll,1)
        vector_ori(2)=ystore(igll,jgll,1)
        vector_ori(3) =zstore(igll,jgll,1)

        do i = 1,NDIM
           vector_rotated(i) = 0.d0
           do j = 1,NDIM
              vector_rotated(i) = vector_rotated(i) + rotation_matrix(i,j)*vector_ori(j)
           enddo
        enddo
        x=vector_rotated(1);y=vector_rotated(2);z=vector_rotated(3)
        rayon = dsqrt(vector_rotated(1)**2 + vector_rotated(2)**2 + vector_rotated(3)**2)

        long=atan2(y,x)
        lati=asin(z/rayon)

        iglob=(ilon)*(NGLLX-1)+igll
        jglob=(ilat)*(NGLLY-1)+jgll
        lon_zmin(iglob,jglob)= long/deg2rad
        lat_zmin(iglob,jglob)= lati/deg2rad

     enddo
  enddo

  end subroutine store_zmin_points

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_stzmin(x,y,nx,ny,MESH)

  implicit none

  integer i,j,nx,ny
  double precision x(nx,ny),y(nx,ny)
  character(len=10) MESH

  open(27,file=trim(MESH)//'stzmin')
  write(27,*) nx*ny
  do j=1,ny
     do i=1,nx
        write(27,*) x(i,j),y(i,j)
     enddo
  enddo
  close(27)

  end subroutine write_stzmin

!=======================================================================================================
!
!=======================================================================================================

  subroutine write_Igm_file(iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il)

  implicit none

  integer iunit,ispec2D,NGLL1,NGLL2,ie,je,js,il
  integer i,j
  do j=1,NGLL2
     do i=1,NGLL1
        write(iunit,*) i,j,ispec2D,(NGLL1-1)*ie+i,(NGLL2-1)*je+j+js,il
     enddo
  enddo

  end subroutine write_Igm_file

!=======================================================================================================
!
!=======================================================================================================

  subroutine  compute_rotation_matrix(rotation_matrix, lon_center_chunk,lat_center_chunk, chunk_azi)

  implicit none

  double precision rotation_matrix(3,3),lon_center_chunk,lat_center_chunk, chunk_azi
  double precision R0(3,3),R1(3,3),R2(3,3),axe_rotation(3),R00(3,3)

  ! je met le chunk en 0,0
  axe_rotation(1)=0.d0; axe_rotation(2)=1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R00,axe_rotation,90.d0)  ! je ramene le chunk en (0,0)
  ! rotation de l'azimuth du chunk
  axe_rotation(1)=1.d0; axe_rotation(2)=0.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R0,axe_rotation,90.-chunk_azi)
  ! on met le chunk a la bonne latitude
  axe_rotation(1)=0.d0; axe_rotation(2)=-1.d0; axe_rotation(3)=0.d0
  call rotation_matrix_axe(R1,axe_rotation,lat_center_chunk)
  ! on met le chunk a la bonne longitude
  axe_rotation(1)=0.d0; axe_rotation(2)=0.d0; axe_rotation(3)=1.d0
  call rotation_matrix_axe(R2,axe_rotation, lon_center_chunk)
  ! rotation resultante
  call compose4matrix(rotation_matrix,R00,R0,R1,R2)

  end subroutine compute_rotation_matrix

!=======================================================================================================
!
!   ROUTINES POUR FAIRE DES ROTATIONS 3D ET DIVERS CHANGEMENTS DE REPERES
!
! Vadim Monteiller Mars 2013
!
!-------------------------------------------------------------------------------
! matrice de rotation 3D d'axe "axe" et d'angle theta (en degres)
! cette matrice est en complexe
!
!=======================================================================================================
!
  subroutine rotation_matrix_axe(R,axe,theta)

  implicit none

  double precision axe(3),theta,pi,deg2rad
  double precision R(3,3)
  double precision c,s,ux,uy,uz,norme_axe

  pi=3.1415926535897932d0
  deg2rad = pi / 180.d0
  ! on normalise l'axe
  norme_axe=dsqrt(axe(1)**2 + axe(2)**2 + axe(3)**2)

  ! composantes de l'axe
  ux=axe(1)/norme_axe
  uy=axe(2)/norme_axe
  uz=axe(3)/norme_axe

  ! on calcule le cos et sin
  c=dcos(deg2rad * theta);s=dsin(deg2rad * theta)

  ! matrice de rotation complexe
  R(1,1)=(ux**2 + (1.d0-ux**2)*c)
  R(1,2)=(ux*uy*(1.d0-c)-uz*s)
  R(1,3)=(ux*uz*(1.d0-c)+uy*s)

  R(2,1)=(ux*uy*(1.d0-c)+uz*s)
  R(2,2)=(uy**2+(1.d0-uy**2)*c)
  R(2,3)=(uy*uz*(1.d0-c)-ux*s)

  R(3,1)=(ux*uz*(1.d0-c)-uy*s)
  R(3,2)=(uy*uz*(1.d0-c)+ux*s)
  R(3,3)=(uz**2+(1.d0-uz**2)*c)

  !write(49,*) ' MATRICE ROTATION '
  !write(49,*) R(1,:)
  !write(49,*) R(2,:)
  !write(49,*) R(3,:)
  !write(49,*)

  end subroutine rotation_matrix_axe

!=======================================================================================================
!
! R=R2*R1*R0
!
!=======================================================================================================

  subroutine compose4matrix(R,R00,R0,R1,R2)

  implicit none

  double precision R(3,3),R0(3,3),R1(3,3),R2(3,3),R00(3,3),Rtmp(3,3)
  integer i,j,k


  R(:,:)=0.d0
  ! multiplication R=R0*R00
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R0(i,k)*R00(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R1*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R1(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  ! multiplication R=R2*R
  Rtmp=R
  R(:,:)=0.d0
  do j=1,3
     do i=1,3
        do k=1,3
           R(i,j)=R(i,j) + R2(i,k)*Rtmp(k,j)
        enddo
     enddo
  enddo

  end subroutine compose4matrix

!------------------------------------------------------------------------------
! rotation pour passer d'un repere local a un autre
!===========================================================================!
!

subroutine Lyfnd(r,rb,n,i)

  implicit none

  integer i,n
  double precision r,rb(n)

  i=1
  do while (r > rb(i) )
     i = i + 1
  enddo
  i = i - 1

end subroutine Lyfnd

function IsNewLayer(x,r,n)
  implicit none
  integer IsNewLayer,n,i
  double precision x,r(n)
  IsNewLayer = 0
  ! ce test fonctionne que si les mailles sont suffisament petites !! ATTENTION
  do i = 1, n-1
     if (abs(x-r(i)) < 1.d-10) then
        IsNewLayer = 1
        return
     endif
  enddo
end function IsNewLayer


subroutine StorePoint(z,k,zc)

  implicit none

  integer k
  double precision z(*),zc

  if (k == 0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k) == zc) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePoint

subroutine StorePointZ(z,k,zc,NoInter)

  implicit none

  integer k
  double precision z(*),zc
  logical NoInter

  if (k == 0) then
     k = k + 1
     z(k) = zc
     return
  else
     if (z(k) == zc .and. NoInter) then
        return
     else
        k = k + 1
        z(k) = zc
     endif
  endif
end subroutine StorePointZ

 subroutine CalGridProf(ProfForGemini,Niveau_elm,zlayer,nlayer,NEX_GAMMA,Z_DEPTH_BLOCK)

  implicit none
  integer NEX_GAMMA,nlayer,nbbloc(100000),Niveau_elm(0:NEX_GAMMA-1)
  double precision ProfForGemini(0:NEX_GAMMA-1,3),zlayer(nlayer)
  double precision Z_DEPTH_BLOCK,zpoint(100000),zz(100000)
  double precision epsillon
  integer nb, n, i,j,k,ilayer,ilay,nd,niveau
  double precision p, pas, longeur
  logical test

  epsillon=1d-3
   nbbloc(:)=0
   ! point de depart
   zpoint(1)=zlayer(nlayer) - Z_DEPTH_BLOCK
   write(*,*) zlayer(nlayer) ,  Z_DEPTH_BLOCK
   !! niveau de depart
   call FindLayer_for_earth_chunk_mesh(ilayer,zlayer,zpoint(1),nlayer)
   write(*,*) '              INITIALISATION calcul du niveau de depart : '
   write(*,*)
   write(*,*) 'zlayer : ', zlayer
   write(*,*) 'premier point : '   , zpoint(1),ilayer
    write(*,*)

  !! on compte le nombre d'elements par niveau
  i = 1
  k = ilayer - 1
  nb = 0
  do while (zpoint(i) < zlayer(nlayer))
    i = i + 1
    k = k + 1
    zpoint(i) = zlayer(k)
  enddo

  nb = i
  nd = i-1
  longeur = zlayer(nlayer) - zpoint(1)


  do i=1,nb-1

     pas = zpoint(i+1) - zpoint(i)
     p = NEX_GAMMA * pas / longeur

    if (p < 0.8d0) then
        n = 1
    else
        n = max(int(p),2)
    endif

    nbbloc(i)=n

  enddo

  do j=1,nb-1
    write(*,*) j,nbbloc(j)
  enddo

  !! on elimine les blocs en trop
   write(*,*) 'SUM ',sum(nbbloc)

   nb = sum(nbbloc)

   do while (nb > NEX_GAMMA)

      k  =  1
      test = .true.

    do  while (test)

         j =  maxval(nbbloc)
         ! on cherche l'indice du max

         if (j == nbbloc(k)) then
            nbbloc(k ) = nbbloc(k) -1
            test = .false.
         endif

         k = k + 1

      enddo

      nb = sum(nbbloc)
      write(*,*) 'nb, ',nb,NEX_GAMMA
   enddo

  longeur = zlayer(nlayer) - zpoint(1)
  k=1
  zz(k)=zpoint(1)
  do i=1,nd
     pas = (zpoint(i+1) - zpoint(i)) / nbbloc(i)
     write(*,*) i,nbbloc(i),pas
     do while (zz(k) < zpoint(i+1) - epsillon)
        k = k + 1
        zz(k) = zz(k-1) + pas
        write(*,*) zz(k), zpoint(i+1)
     enddo
  enddo

   do ilay=1,NEX_GAMMA

      ProfForGemini(ilay-1,1)  =  zz(ilay)
      ProfForGemini(ilay-1,2)  =  zz(ilay+1)
      ProfForGemini(ilay-1,3)  = 0.5d0 * (zz(ilay) + zz(ilay+1))

      call FindLayer_for_earth_chunk_mesh(niveau,zlayer, ProfForGemini(ilay-1,3),nlayer)
      Niveau_elm(ilay-1)=niveau
      write(*,'(i5,2f15.3,i10)') ilay,zz(ilay),zz(ilay+1),niveau
   enddo

 end subroutine CalGridProf

 subroutine  FindLayer_for_earth_chunk_mesh(i,z,r,n)
   implicit none
   integer i,n
   double precision z(n),r

   if (r > z(n) .or. r < z(1)) then
    write(*,*) 'STOP :: point ouside grid'
    stop
   endif
   i = 1
   do while (r > z(i))
     i = i + 1
   enddo


 end subroutine FindLayer_for_earth_chunk_mesh


!! VM VM add this for Axisem coupling

subroutine  find_layer_in_axisem_model(i,u,r,z,n)

   implicit none

   integer i,n,u(5)
   double precision z(n),r(5)

   if (r(3) > z(n) .or. r(3) < z(1)) then
      write(*,*) 'STOP :: point ouside grid'
      stop
   endif
   i = 1
   do while (r(3) > z(i))
      i = i + 1
   enddo

   u(:)=0

end subroutine find_layer_in_axisem_model


!=====================================================================

  subroutine getglob_for_chunk(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD,UTM_X_MIN,UTM_X_MAX)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  integer NGNOD

  integer npointot
  integer nspec,nglob
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  double precision UTM_X_MIN,UTM_X_MAX

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! geometry tolerance parameter to calculate number of independent grid points
! small value for double precision and to avoid sensitivity to roundoff
  double precision SMALLVALTOL

! define geometrical tolerance based upon typical size of the model
  SMALLVALTOL = 1.d-10 * dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) dabs(UTM_X_MAX - UTM_X_MIN)
  write(*,*) ' SMALLVALTOL  ',SMALLVALTOL
! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers (!! VM changed NGLLCUBE (as in Specfem3D Basin Version 1.1) to NGNOD !!)
  do ispec=1,nspec
    ieoff = NGNOD * (ispec - 1)
    do ilocnum = 1,NGNOD
      loc(ilocnum + ieoff) = ilocnum + ieoff
    enddo
  enddo

  ifseg(:)  = .false.
  nseg      = 1
  ifseg(1)  = .true.
  ninseg(1) = npointot

  do j=1,3 !,NDIM

! sort within each segment
  ioff = 1

  do iseg=1,nseg

    if (j == 1) then
      call rank(xp(ioff),ind,ninseg(iseg))
    else if (j == 2) then
      call rank(yp(ioff),ind,ninseg(iseg))
    else
      call rank(zp(ioff),ind,ninseg(iseg))
    endif

    call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))

    ioff = ioff + ninseg(iseg)

  enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
  if (j == 1) then

    do i=2,npointot
      if (dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  else if (j == 2) then

    do i=2,npointot
      if (dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  else

    do i=2,npointot
      if (dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i) = .true.
    enddo

  endif

! count up number of different segments
  nseg = 0

  do i=1,npointot
    if (ifseg(i)) then
      nseg = nseg + 1
      ninseg(nseg) = 1
    else
      ninseg(nseg) = ninseg(nseg) + 1
    endif
  enddo
  enddo

! assign global node numbers (now sorted lexicographically)
  ig = 0

  do i=1,npointot
    if (ifseg(i)) ig = ig + 1

    iglob(loc(i)) = ig
  enddo

  nglob = ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine getglob_for_chunk


! -----------------------------------

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 continue
   if (l > 1) then
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   endif
   i=l
   j=l+l
  200    continue
   if (J <= IR) then
      if (J < IR) then
         if (A(IND(j)) < A(IND(j+1))) j=j+1
      endif
      if (q < A(IND(j))) then
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   goto 200
   endif
   IND(I)=INDX
  goto 100
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all

!----------------

subroutine calc_gll_points(xelm,yelm,zelm,xstore,ystore,zstore,shape3D,NGNOD,NGLLX,NGLLY,NGLLZ)

  implicit none
  integer NGNOD,NGLLX,NGLLY,NGLLZ
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision xstore(NGLLX,NGLLY,NGLLZ),ystore(NGLLX,NGLLY,NGLLZ),zstore(NGLLX,NGLLY,NGLLZ)
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD),xmesh,ymesh,zmesh
  double precision, parameter :: ZERO = 0.d0
  integer ia,i,j,k

  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX
           ! compute mesh coordinates
           xmesh = ZERO
           ymesh = ZERO
           zmesh = ZERO
           do ia=1,NGNOD
              xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
              ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
              zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
           enddo
           xstore(i,j,k) = xmesh
           ystore(i,j,k) = ymesh
           zstore(i,j,k) = zmesh
        enddo
     enddo
  enddo

end subroutine calc_gll_points

