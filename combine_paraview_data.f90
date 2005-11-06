!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program combine_paraview_data

! puts the output of SPECFEM3D in ParaView format.
! see http://www.paraview.org for details

! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

  implicit none

  include "constants.h"

  integer i,j,k,ispec, ios, it
  integer iproc, proc1, proc2, num_node, node_list(300)
  integer np, ne, npp, nee, npoint, nelement, njunk, njunk2, n1, n2, n3, n4, n5, n6, n7, n8
  integer numpoin, iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8, iglob

  real x, y, z

  logical, dimension(:), allocatable :: mask_ibool

  integer, dimension(:), allocatable :: num_ibool
  integer, dimension(:,:,:,:), allocatable :: ibool

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore, ystore, zstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: data,dat

  character(len=150) :: sline, arg(5), filename, indir, outdir, prname
  character(len=150) :: mesh_file, local_point_file, local_element_file, local_file, local_data_file, local_ibool_file
  logical HIGH_RESOLUTION_MESH
  integer ires

! parameters read from parameter file
  integer NER_SEDIM,NER_BASEMENT_SEDIM,NER_16_BASEMENT, &
             NER_MOHO_16,NER_BOTTOM_MOHO,NEX_XI,NEX_ETA, &
             NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE,SIMULATION_TYPE
  integer NSOURCES

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision DT,LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX,HDUR_MOVIE
  double precision THICKNESS_TAPER_BLOCK_HR,THICKNESS_TAPER_BLOCK_MR,VP_MIN_GOCAD,VP_VS_RATIO_GOCAD_TOP,VP_VS_RATIO_GOCAD_BOTTOM

  logical HARVARD_3D_GOCAD_MODEL,TOPOGRAPHY,ATTENUATION,USE_OLSEN_ATTENUATION, &
          OCEANS,IMPOSE_MINIMUM_VP_GOCAD,HAUKSSON_REGIONAL_MODEL, &
          BASEMENT_MAP,MOHO_MAP_LUPEI,ABSORBING_CONDITIONS,SAVE_FORWARD
  logical ANISOTROPY,SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION
  logical MOVIE_SURFACE,MOVIE_VOLUME,CREATE_SHAKEMAP,SAVE_DISPLACEMENT, &
          USE_HIGHRES_FOR_MOVIES,SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  integer NTSTEP_BETWEEN_FRAMES,NTSTEP_BETWEEN_OUTPUT_INFO

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER

  integer nspec,NSPEC2D_A_XI,NSPEC2D_B_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,nglob

  print *
  print *,'Recombining ParaView data for slices'
  print *

  do i = 1, 6
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print *, 'Usage: xcombine_data start_slice end_slice filename input_dir output_dir high/low-resolution'
      print *, '    or xcombine_data slice_list filename input_dir output_dir high/low-resolution'
      print *, ' possible filenames are '
      print *, '   rho_vp, rho_vs, kappastore, mustore etc'
      print *, '   that are stored in the local directory as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec)  '
      print *, '   in filename.bin'
      print *, ' files have been collected in input_dir, output mesh file goes to output_dir '
      print *, ' give 0 for low resolution and 1 for high resolution'
      stop ' Reenter command line options'
    endif
  enddo

  print *
  print *,'reading parameter file'
  print *

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
        NTSTEP_BETWEEN_FRAMES,USE_HIGHRES_FOR_MOVIES,,HDUR_MOVIE, &
        SAVE_MESH_FILES,PRINT_SOURCE_TIME_FUNCTION, &
        NTSTEP_BETWEEN_OUTPUT_INFO,SUPPRESS_UTM_PROJECTION,MODEL,USE_REGULAR_MESH,SIMULATION_TYPE,SAVE_FORWARD)

! compute other parameters based upon values read
  call compute_parameters(NER,NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER_BOTTOM_MOHO,NER_MOHO_16,NER_16_BASEMENT,NER_BASEMENT_SEDIM,NER_SEDIM, &
      nspec,NSPEC2D_A_XI,NSPEC2D_B_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,nglob,USE_REGULAR_MESH)

  allocate(mask_ibool(nglob))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(data(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(dat(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(nglob))
  allocate(ystore(nglob))
  allocate(zstore(nglob))
  allocate(num_ibool(nglob))

! get slice list
  if (trim(arg(6)) == '') then
    num_node = 0
    open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ios)
    if (ios /= 0) then
      print *,'Error opening ',trim(arg(1))
      stop
    endif
    do while ( 1 == 1)
      read(20,'(a)',iostat=ios) sline
      if (ios /= 0) exit
      read(sline,*,iostat=ios) njunk
      if (ios /= 0) exit
      num_node = num_node + 1
      node_list(num_node) = njunk
    enddo
    close(20)
    filename = arg(2)
    indir= arg(3)
    outdir = arg(4)
    read(arg(5),*) ires
  else
    read(arg(1),*) proc1
    read(arg(2),*) proc2
    do iproc = proc1, proc2
      node_list(iproc - proc1 + 1) = iproc
    enddo
    num_node = proc2 - proc1 + 1
    filename = arg(3)
    indir = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  endif

  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
  else
    HIGH_RESOLUTION_MESH = .true.
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)

  ! open paraview output mesh file
  mesh_file = trim(outdir) // '/' // trim(filename)//'.mesh'
  call open_file(trim(mesh_file)//char(0))

  np = 0

  ! write point and scalar information

  do it = 1, num_node

    iproc = node_list(it)

    print *, ' '
    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a)') trim(indir)//'/proc',iproc,'_'

  ! data file
    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 27,file = trim(local_data_file),status='old', iostat = ios,form ='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(27) data
    close(27)
    print *, trim(local_data_file)

    if (CUSTOM_REAL == SIZE_DOUBLE) then
      dat = sngl(data)
    else
      dat = data
    endif

  ! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old', iostat = ios, form='unformatted')
    if (ios /= 0) then
      print *,'Error opening ',trim(local_data_file)
      stop
    endif
    read(28) ibool
    close(28)
    print *, trim(local_ibool_file)

     mask_ibool(:) = .false.
     numpoin = 0

    if (.not. HIGH_RESOLUTION_MESH) then

      local_point_file = trim(prname) // 'AVS_DXpoints.txt'
      open(unit = 25, file = trim(local_point_file), status = 'old', iostat = ios)
      if (ios /= 0) then
        print *,'Error opening ',trim(local_point_file)
        stop
      endif
      read(25,*) npoint

      if (it == 1) then
        npp = npoint * num_node
        call write_integer(npp)
      endif

      do ispec=1,nspec
        iglob1=ibool(1,1,1,ispec)
        iglob2=ibool(NGLLX,1,1,ispec)
        iglob3=ibool(NGLLX,NGLLY,1,ispec)
        iglob4=ibool(1,NGLLY,1,ispec)
        iglob5=ibool(1,1,NGLLZ,ispec)
        iglob6=ibool(NGLLX,1,NGLLZ,ispec)
        iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8=ibool(1,NGLLY,NGLLZ,ispec)

        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,1,1,ispec))
          mask_ibool(iglob1) = .true.
        endif
        if(.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,1,1,ispec))
          mask_ibool(iglob2) = .true.
        endif
        if(.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,NGLLY,1,ispec))
          mask_ibool(iglob3) = .true.
        endif
        if(.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,NGLLY,1,ispec))
          mask_ibool(iglob4) = .true.
        endif
        if(.not. mask_ibool(iglob5)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,1,NGLLZ,ispec))
          mask_ibool(iglob5) = .true.
        endif
        if(.not. mask_ibool(iglob6)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,1,NGLLZ,ispec))
          mask_ibool(iglob6) = .true.
        endif
        if(.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,NGLLY,NGLLZ,ispec))
          mask_ibool(iglob7) = .true.
        endif
        if(.not. mask_ibool(iglob8)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,NGLLY,NGLLZ,ispec))
          mask_ibool(iglob8) = .true.
        endif
      enddo ! ispec
      close(25)

    else  ! high resolution

      if (it == 1) then
        npp = nglob * num_node
        npoint = nglob
        call write_integer(npp)
      endif

      local_file = trim(prname)//'x.bin'
      open(unit = 27,file = trim(prname)//'x.bin',status='old', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print *,'Error opening ',trim(local_file)
        stop
      endif
      read(27) xstore
      close(27)
      local_file = trim(prname)//'y.bin'
      open(unit = 27,file = trim(prname)//'y.bin',status='old', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print *,'Error opening ',trim(local_file)
        stop
      endif
      read(27) ystore
      close(27)
      local_file = trim(prname)//'z.bin'
      open(unit = 27,file = trim(prname)//'z.bin',status='old', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print *,'Error opening ',trim(local_file)
        stop
      endif
      read(27) zstore
      close(27)

      do ispec=1,nspec
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              if(.not. mask_ibool(iglob)) then
                numpoin = numpoin + 1
                if (CUSTOM_REAL == SIZE_DOUBLE) then
                  x = sngl(xstore(iglob))
                  y = sngl(ystore(iglob))
                  z = sngl(zstore(iglob))
                else
                  x = xstore(iglob)
                  y = ystore(iglob)
                  z = zstore(iglob)
                endif
                call write_real(x)
                call write_real(y)
                call write_real(z)
                call write_real(dat(i,j,k,ispec))
                mask_ibool(iglob) = .true.
              endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec
    endif

    if (numpoin /= npoint) stop 'Error: number of points are not consistent'
    np = np + npoint

  enddo  ! all slices for points

 if (np /=  npp) stop 'Error: Number of total points are not consistent'
 print *, 'Total number of points: ', np
 print *, ' '


 ne = 0
! write element information
 do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a)') trim(indir)//'/proc',iproc,'_'

    np = npoint * (it-1)

    if (.not. HIGH_RESOLUTION_MESH) then

      local_element_file = trim(prname) // 'AVS_DXelements.txt'
      open(unit = 26, file = trim(local_element_file), status = 'old', iostat = ios)
      if (ios /= 0) then
        print *,'Error opening ',trim(local_element_file)
        stop
      endif
      print *, trim(local_element_file)

      read(26, *) nelement
      if (it == 1) then
        nee = nelement * num_node
        call write_integer(nee)
      end if

      do i = 1, nelement
        read(26,*) njunk, njunk2, n1, n2, n3, n4, n5, n6, n7, n8
        n1 = n1+np-1
        n2 = n2+np-1
        n3 = n3+np-1
        n4 = n4+np-1
        n5 = n5+np-1
        n6 = n6+np-1
        n7 = n7+np-1
        n8 = n8+np-1
        call write_integer(n1)
        call write_integer(n2)
        call write_integer(n3)
        call write_integer(n4)
        call write_integer(n5)
        call write_integer(n6)
        call write_integer(n7)
        call write_integer(n8)
      enddo
      close(26)

    else ! high resolution mesh

      if (it == 1) then
        nelement = nspec * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
        nee = nelement * num_node
        call write_integer(nee)
      endif

      numpoin = 0
      mask_ibool = .false.
      do ispec=1,nspec
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              if(.not. mask_ibool(iglob)) then
                numpoin = numpoin + 1
                num_ibool(iglob) = numpoin
                mask_ibool(iglob) = .true.
              endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec

      do ispec = 1, nspec
        do k = 1, NGLLZ-1
          do j = 1, NGLLY-1
            do i = 1, NGLLX-1
              iglob1 = ibool(i,j,k,ispec)
              iglob2 = ibool(i+1,j,k,ispec)
              iglob3 = ibool(i+1,j+1,k,ispec)
              iglob4 = ibool(i,j+1,k,ispec)
              iglob5 = ibool(i,j,k+1,ispec)
              iglob6 = ibool(i+1,j,k+1,ispec)
              iglob7 = ibool(i+1,j+1,k+1,ispec)
              iglob8 = ibool(i,j+1,k+1,ispec)
              n1 = num_ibool(iglob1)+np-1
              n2 = num_ibool(iglob2)+np-1
              n3 = num_ibool(iglob3)+np-1
              n4 = num_ibool(iglob4)+np-1
              n5 = num_ibool(iglob5)+np-1
              n6 = num_ibool(iglob6)+np-1
              n7 = num_ibool(iglob7)+np-1
              n8 = num_ibool(iglob8)+np-1
              call write_integer(n1)
              call write_integer(n2)
              call write_integer(n3)
              call write_integer(n4)
              call write_integer(n5)
              call write_integer(n6)
              call write_integer(n7)
              call write_integer(n8)
            enddo
          enddo
        enddo
      enddo

    endif
    ne = ne + nelement

  enddo ! num_node
  if (ne /= nee) stop 'Number of total elements are not consistent'
  print *, 'Total number of elements: ', ne

  call close_file()

  print *, 'Done writing '//trim(mesh_file)

  end program combine_paraview_data

