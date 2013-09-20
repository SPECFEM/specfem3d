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
! MPI v. 1.0 Dimitri Komatitsch, Caltech, May 2002: first MPI version based on global code

!
!-------------------------------------------------------------------------------------------------
!


  subroutine generate_databases

  use adios_manager_mod
  use mpi
  use generate_databases_par

  implicit none

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
    call flush_IMAIN()
  endif

! read the parameter file
  call gd_read_parameters()

! reads topography and bathymetry file
  call gd_read_topography()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) 'creating mesh in the model'
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Initialize ADIOS I/O
  if (ADIOS_ENABLED) then
    call adios_setup()
  endif

! reads Databases files
  if (ADIOS_FOR_DATABASES) then
    call read_partition_files_adios()
  else
    call read_partition_files()
  endif

! external mesh creation
  call setup_mesh()

! finalize mesher
  call finalize_databases()

  if (ADIOS_ENABLED) then
    call adios_cleanup()
  endif
  end subroutine generate_databases

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_read_parameters

! reads and checks user input parameters

  use generate_databases_par
  implicit none

! reads Par_file
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

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROC) then
    if( myrank == 0 ) then
      write(IMAIN,*) 'error: number of processors supposed to run on: ',NPROC
      write(IMAIN,*) 'error: number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'error generate_databases: number of processors supposed to run on: ',NPROC
      print*, 'error generate_databases: number of MPI processors actually run on: ',sizeprocs
      print*
    endif
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

! there would be a problem with absorbing boundaries for different NGLLX,NGLLY,NGLLZ values
! just to be sure for now..
  if(STACEY_ABSORBING_CONDITIONS) then
    if( NGLLX /= NGLLY .and. NGLLY /= NGLLZ ) &
      call exit_MPI(myrank,'must have NGLLX = NGLLY = NGLLZ for external meshes')
  endif

! info about external mesh simulation
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
    write(IMAIN,*) 'Beware! Curvature (i.e. HEX27 elements) is not handled by our internal mesher'
    write(IMAIN,*)

! check that the constants.h file is correct
    if ( NGNOD /= 8 .and. NGNOD /= 27 ) then
       stop 'elements should have 8 or 27 control nodes, please modify NGNOD in Par_file'
    endif

    write(IMAIN,'(a)',advance='no') ' velocity model: '
    select case(IMODEL)
    case( IMODEL_DEFAULT )
    write(IMAIN,'(a)',advance='yes') '  default '
    case( IMODEL_GLL )
    write(IMAIN,'(a)',advance='yes') '  gll'
    case( IMODEL_1D_PREM )
    write(IMAIN,'(a)',advance='yes') '  1d_prem'
    case( IMODEL_1D_CASCADIA )
    write(IMAIN,'(a)',advance='yes') '  1d_cascadia'
    case( IMODEL_1D_SOCAL )
    write(IMAIN,'(a)',advance='yes') '  1d_socal'
    case( IMODEL_SALTON_TROUGH )
    write(IMAIN,'(a)',advance='yes') '  salton_trough'
    case( IMODEL_TOMO )
    write(IMAIN,'(a)',advance='yes') '  tomo'
    case( IMODEL_USER_EXTERNAL )
    write(IMAIN,'(a)',advance='yes') '  external'
    case( IMODEL_IPATI )
    write(IMAIN,'(a)',advance='yes') '  ipati'
    case( IMODEL_IPATI_WATER )
    write(IMAIN,'(a)',advance='yes') '  ipati_water'
    end select

    write(IMAIN,*)
  endif

! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) &
    call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) call exit_MPI(myrank,'number of SLS must be 3')

  ! for noise simulations, we need to save movies at the surface (where the noise is generated)
  ! and thus we force MOVIE_SURFACE to be .true., in order to use variables defined for surface movies later
  if ( NOISE_TOMOGRAPHY /= 0 ) then
    MOVIE_TYPE = 1
    MOVIE_SURFACE = .true.
    USE_HIGHRES_FOR_MOVIES = .true.     ! we need to save surface movie everywhere, i.e. at all GLL points on the surface
  endif

  if(myrank == 0) then

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
    if(APPROXIMATE_OCEAN_LOAD) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
      if(TOPOGRAPHY) write(IMAIN,*) ' with elevation from topography file'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)
    if(STACEY_ABSORBING_CONDITIONS) then
      write(IMAIN,*) 'incorporating Stacey absorbing conditions'
    else
      if(PML_CONDITIONS) then
        write(IMAIN,*) 'incorporating absorbing conditions of perfectly matched layer'
      else
        write(IMAIN,*) 'no absorbing condition'
      endif
    endif

    write(IMAIN,*)
    if(USE_FORCE_POINT_SOURCE) then
       write(IMAIN,*) 'using a FORCESOLUTION source instead of a CMTSOLUTION source'
    else
       write(IMAIN,*) 'using a CMTSOLUTION source'
       write(IMAIN,*)
    endif

    write(IMAIN,*)
    if(USE_RICKER_TIME_FUNCTION) then
       write(IMAIN,*) 'using a Ricker source time function'
    else
       if(USE_FORCE_POINT_SOURCE) then
          write(IMAIN,*) 'using a quasi-Heaviside source time function'
          write(IMAIN,*)
       else
          write(IMAIN,*) 'using a Gaussian source time function'
          write(IMAIN,*)
       endif
    endif
    call flush_IMAIN()
  endif

! makes sure processes are synchronized
  call sync_all()

  end subroutine gd_read_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gd_read_topography

! reads in topography files

  use generate_databases_par
  implicit none

  if( APPROXIMATE_OCEAN_LOAD .and. TOPOGRAPHY ) then

    ! values given in constants.h
    NX_TOPO = NX_TOPO_FILE
    NY_TOPO = NY_TOPO_FILE
    allocate(itopo_bathy(NX_TOPO,NY_TOPO),stat=ier)
    if( ier /= 0 ) stop 'error allocating array itopo_bathy'

    call read_topo_bathy_file(itopo_bathy,NX_TOPO,NY_TOPO)

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'regional topography file read ranges in m from ',minval(itopo_bathy),' to ',maxval(itopo_bathy)
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  else
    NX_TOPO = 1
    NY_TOPO = 1
    allocate(itopo_bathy(NX_TOPO,NY_TOPO),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array itopo_bathy'

  endif

  end subroutine gd_read_topography

