!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
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

!=============================================================================!
!                                                                             !
!  meshfem3D produces a spectral element grid for a local or regional model.  !
!  The mesher uses the UTM projection                                         !
!                                                                             !
!=============================================================================!
!
! Please find in the header of specfem3D.F90 further code informations.
!
! ************** PROGRAM STARTS HERE **************

  program xmeshfem3D

  use meshfem_par
  use chunk_earth_mod

  implicit none

  include 'version.fh'

  ! local parameters
  integer :: iprocnum
  integer :: iproc_xi,iproc_eta
  integer :: ier
  integer(kind=8) :: nglob_total
  ! interface parameters
  integer :: ilayer

  ! timer MPI
  double precision, external :: wtime
  double precision :: time_start,tCPU
  character(len=3) :: str_unit

  ! MPI initialization
  call init_mpi()

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
    write(IMAIN,*) 'Running Git package version of the code: ', git_package_version
    write(IMAIN,*) 'which is Git ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (myrank == 0) then
    write(IMAIN,*) 'Reading parameters from ',IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read the parameter file (DATA/Par_file)
  BROADCAST_AFTER_READ = .true.
  call read_parameter_file(BROADCAST_AFTER_READ)

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

    ! stop program
    call finalize_mpi()
    stop

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
  call synchronize_all()

  ! dynamic allocation of mesh arrays
  allocate(xgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1353')
  if (ier /= 0) stop 'Error allocating array xgrid'
  xgrid(:,:,:) = 0.d0

  allocate(ygrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1354')
  if (ier /= 0) stop 'Error allocating array ygrid'
  ygrid(:,:,:) = 0.d0

  allocate(zgrid(0:2*NER,0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1355')
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate arrays')
  zgrid(:,:,:) = 0.d0

  allocate(addressing(0:NPROC_XI-1,0:NPROC_ETA-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1356')
  if (ier /= 0) stop 'Error allocating array addressing'
  addressing(:,:) = 0

  allocate(iproc_xi_slice(0:NPROC-1), &
           iproc_eta_slice(0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1358')
  if (ier /= 0) stop 'Error allocating array iproc_eta_slice'
  ! clear arrays
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
    write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD,' control nodes'
    write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D,' control nodes'
    write(IMAIN,*) 'Beware! Curvature (i.e. HEX27 elements) is not handled by our internal mesher'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check that the constants.h file is correct
  if (NGNOD /= 8 .and. NGNOD /= 27) &
    call exit_MPI(myrank,'Error must have set NGNOD == 8 or NGNOD == 27')
  if (NGNOD2D /= 4 .and. NGNOD2D /= 9) &
    call exit_MPI(myrank,'Error must have set NGNOD2D == 4 or NGNOD2D == 9')

  if (NGLLX_M == 2 .and. NGLLY_M == 2 .and. NGLLZ_M == 2) then
    if (NGNOD /= 8) &
      call exit_MPI(myrank,'With NGLLX_M == 2, volume elements should have NGNOD == 8 control nodes in our internal mesher')
    if (NGNOD2D /= 4) &
      call exit_MPI(myrank,'With NGLLX_M == 2, surface elements should have NGNOD2D == 4 control nodes in our internal mesher')
  endif

  if (NGNOD == 27 .and. (NGLLX_M < 3 .or. NGLLY_M < 3 .or. NGLLZ_M < 3)) &
    call exit_MPI(myrank,'NGNOD = 27 control nodes needs at least NGLLX_M == NGLLY_M == NGLLZ_M >= 3 in our internal mesher')
  if (NGNOD2D == 9 .and. (NGLLX_M < 3 .or. NGLLY_M < 3 .or. NGLLZ_M < 3)) &
    call exit_MPI(myrank,'NGNOD2D = 9 control nodes needs at least NGLLX_M == NGLLY_M == NGLLZ_M >= 3 in our internal mesher')

!  if (.not. USE_REGULAR_MESH .and. NGLLX_M >= 3) &
!    call exit_MPI(myrank,'NGLLX_M == NGLLY_M == NGLLZ_M >= 3 only supported with USE_REGULAR_MESH = .true. at the moment')

  ! check that reals are either 4 or 8 bytes
  if (CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) call exit_MPI(myrank,'wrong size of CUSTOM_REAL for reals')

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
      if (SUPPRESS_UTM_PROJECTION) then
        ! no UTM, thickness given in m
        str_unit = '(m)'
      else
        ! UTM, thickness given in degree
        str_unit = 'deg'
      endif
      write(IMAIN,*)
      write(IMAIN,*) 'PML thickness in X direction = ',sngl(THICKNESS_OF_X_PML),str_unit
      write(IMAIN,*) 'PML thickness in Y direction = ',sngl(THICKNESS_OF_Y_PML),str_unit
      write(IMAIN,*) 'PML thickness in Z direction = ',sngl(THICKNESS_OF_Z_PML),str_unit
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

  ! to avoid overflow for large meshes
  nglob_total = int( dble(NGLOB_AB)*dble(NPROC) ,kind=8)

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
    write(IMAIN,*) 'approximate total number of points in entire mesh (with duplicates on MPI edges): ',nglob_total
    write(IMAIN,*) 'approximate total number of DOFs in entire mesh (with duplicates on MPI edges): ',nglob_total*NDIM
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

  call finalize_mpi()

  end program xmeshfem3D

