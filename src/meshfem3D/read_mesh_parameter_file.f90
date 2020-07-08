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

  subroutine read_mesh_parameter_file()

  use meshfem3D_par, only: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
    UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
    NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE, &
    LOCAL_PATH,SUPPRESS_UTM_PROJECTION, &
    INTERFACES_FILE,CAVITY_FILE,NSUBREGIONS,subregions,NMATERIALS,material_properties, &
    CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
    USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
    THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML, &
    myrank,sizeprocs,NUMBER_OF_MATERIAL_PROPERTIES

  use constants, only: IIN,MF_IN_DATA_FILES,MAX_STRING_LEN,IMAIN, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC, &
    ILONGLAT2UTM,IGNORE_JUNK,DONT_IGNORE_JUNK

  implicit none

! local variables
  integer :: NEX_MAX
  double precision :: UTM_MAX
  double precision :: DEPTH_BLOCK_KM

  integer :: ier

  integer :: ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer :: iz_beg_region,iz_end_region,imaterial_number

  double precision :: rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag

  integer :: ireg,imat,ndef,nundef
  integer :: mat_id,domain_id
  integer :: ner_value,ner_max,idoub
  logical :: found
  character(len=MAX_STRING_LEN) :: filename

  ! Mesh Parameter file
  filename = MF_IN_DATA_FILES(1:len_trim(MF_IN_DATA_FILES)) // 'Mesh_Par_file'

  if (myrank == 0) then
    write(IMAIN,*) 'Reading mesh parameters from file ',trim(filename)
    call flush_IMAIN()
  endif

  ! note: please be careful that the order of reading in parameters here
  !       must match the order of appearance in the Mesh_Par_file
  !
  ! open parameter file Mesh_Par_file
  call open_parameter_file_mesh(filename)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  input parameters...'
    call flush_IMAIN()
  endif

  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MIN, 'LATITUDE_MIN', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LATITUDE_MIN'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MAX, 'LATITUDE_MAX', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LATITUDE_MAX'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MIN, 'LONGITUDE_MIN', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LONGITUDE_MIN'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MAX, 'LONGITUDE_MAX', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LONGITUDE_MAX'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM, 'DEPTH_BLOCK_KM', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter DEPTH_BLOCK_KM'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,UTM_PROJECTION_ZONE, 'UTM_PROJECTION_ZONE', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter UTM_PROJECTION_ZONE'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,SUPPRESS_UTM_PROJECTION, 'SUPPRESS_UTM_PROJECTION', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter SUPPRESS_UTM_PROJECTION'
  call read_value_string_mesh(IIN,IGNORE_JUNK,INTERFACES_FILE, 'INTERFACES_FILE', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter INTERFACES_FILE'
  call read_value_string_mesh(IIN,IGNORE_JUNK,CAVITY_FILE, 'CAVITY_FILE', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CAVITY_FILE'

  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_XI, 'NEX_XI', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NEX_XI'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_ETA, 'NEX_ETA', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NEX_ETA'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_XI, 'NPROC_XI', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NPROC_XI'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_ETA, 'NPROC_ETA', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NPROC_ETA'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  doubling layers...'
    call flush_IMAIN()
  endif

  call read_value_logical_mesh(IIN,IGNORE_JUNK,USE_REGULAR_MESH, 'USE_REGULAR_MESH', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter USE_REGULAR_MESH'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NDOUBLINGS, 'NDOUBLINGS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NDOUBLINGS'

  ! checks value
  if (USE_REGULAR_MESH) then
    ! sets NDOUBLINGS to zero for regular grids
    NDOUBLINGS = 0
  else
    ! irregular grid with doubling layer
    if (NDOUBLINGS < 1) stop 'Error parameter NDOUBLINGS is invalid! value must be at least 1 for irregular grids'
  endif

  ! allocate doubling array
  allocate(ner_doublings(NDOUBLINGS),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1318')
  if (ier /= 0) stop 'Error allocating ner_doublings array'
  ner_doublings(:) = 0

  do idoub = 1,NDOUBLINGS
    call read_value_doubling_integer_mesh(IIN,DONT_IGNORE_JUNK, ner_value, 'NZ_DOUBLING', ier)
    if (ier /= 0) then
      print *,'Error reading doubling entry for doubling layer: ',idoub
      print *,'Please check NDOUBLINGS value and corresponding lines with NZ_DOUBLING entries'
      stop 'Error reading Mesh parameter NZ_DOUBLING'
    endif
    ner_doublings(idoub) = ner_value
  enddo
  ! jump over unused doubling entries to reach lines with visualization parameters below
  call read_value_doubling_skip_mesh(IIN,DONT_IGNORE_JUNK,'NZ_DOUBLING', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter after NDOUBLINGS'

  ! only fix 2 doubling layer entries
  !call read_value_integer_mesh(IIN,IGNORE_JUNK, ner_doublings(1), 'NZ_DOUBLING_1', ier)
  !if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUBLING_1'
  !call read_value_integer_mesh(IIN,IGNORE_JUNK, ner_doublings(2), 'NZ_DOUBLING_2', ier)
  !if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUBLING_2'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  visualization...'
    call flush_IMAIN()
  endif

  ! visualization file output
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_ABAQUS_FILES, 'CREATE_ABAQUS_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_ABAQUS_FILES'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_DX_FILES, 'CREATE_DX_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_DX_FILES'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_VTK_FILES, 'CREATE_VTK_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_VTK_FILES'

  ! file in which we store the databases
  call read_value_string_mesh(IIN,IGNORE_JUNK,LOCAL_PATH, 'LOCAL_PATH', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LOCAL_PATH'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  CPML...'
    call flush_IMAIN()
  endif

  ! CPML thickness
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_X_PML, 'THICKNESS_OF_X_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_X_PML'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Y_PML, 'THICKNESS_OF_Y_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Y_PML'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Z_PML, 'THICKNESS_OF_Z_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Z_PML'

  if (THICKNESS_OF_X_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_X_PML'
  if (THICKNESS_OF_Y_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_Y_PML'
  if (THICKNESS_OF_Z_PML < 0.0) stop 'Error invalid negative value THICKNESS_OF_Z_PML'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  domain materials...'
    call flush_IMAIN()
  endif

  ! read number of materials
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NMATERIALS, 'NMATERIALS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NMATERIALS'

  ! read materials properties
  allocate(material_properties(NMATERIALS,NUMBER_OF_MATERIAL_PROPERTIES),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1319')
  if (ier /= 0) stop 'Error allocation of material_properties'
  material_properties(:,:) = 0.d0

  do imat = 1,NMATERIALS
    call read_material_parameters(IIN,material_properties,imat,NMATERIALS,ier)
    if (ier /= 0) then
      print *,'Error reading material ',imat,' out of ',NMATERIALS
      stop 'Error reading materials in Mesh_Par_file'
    endif
    ! user output
    domain_id = material_properties(imat,7)
    mat_id = material_properties(imat,8)
    if (myrank == 0) then
      select case(domain_id)
      case(IDOMAIN_ACOUSTIC)
        write(IMAIN,*) '    material ',mat_id,' acoustic'
      case(IDOMAIN_ELASTIC)
        write(IMAIN,*) '    material ',mat_id,' elastic'
      case (IDOMAIN_POROELASTIC)
        write(IMAIN,*) '    material ',mat_id,' poroelastic'
      end select
      call flush_IMAIN()
    endif
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  domain regions...'
    call flush_IMAIN()
  endif

  ! read number of subregions
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NSUBREGIONS, 'NSUBREGIONS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NSUBREGIONS'

  ! read subregions properties
  allocate(subregions(NSUBREGIONS,7),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1320')
  if (ier /= 0) stop 'Error allocation of subregions'
  subregions(:,:) = 0
  do ireg = 1,NSUBREGIONS
    call read_region_parameters(IIN,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region, &
                                iz_beg_region,iz_end_region,imaterial_number,ier)
    if (ier /= 0) then
      print *,'Error reading region ',ireg,' out of ',NSUBREGIONS
      stop 'Error reading regions in Mesh_Par_file'
    endif
    ! check for negative values: if ix or iy == -1, it means we take full range
    if (ix_end_region == -1) ix_end_region = NEX_XI
    if (iy_end_region == -1) iy_end_region = NEX_ETA

    ! stores region
    subregions(ireg,1) = ix_beg_region
    subregions(ireg,2) = ix_end_region
    subregions(ireg,3) = iy_beg_region
    subregions(ireg,4) = iy_end_region
    subregions(ireg,5) = iz_beg_region
    subregions(ireg,6) = iz_end_region
    subregions(ireg,7) = imaterial_number

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    region ',ireg,' with material ',imaterial_number
      write(IMAIN,*) '      nex_xi  begin/end = ',ix_beg_region,ix_end_region
      write(IMAIN,*) '      nex_eta begin/end = ',iy_beg_region,iy_end_region
      write(IMAIN,*) '      nz      begin/end = ',iz_beg_region,iz_end_region
      call flush_IMAIN()
    endif
  enddo

  ! close parameter file
  call close_parameter_file_mesh()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  reading Mesh_Par_file done successfully'
    write(IMAIN,*)
    write(IMAIN,*) '  checking mesh setup...'
    call flush_IMAIN()
  endif

  ! convert model size to UTM coordinates and depth of mesh to meters
  call utm_geo(LONGITUDE_MIN,LATITUDE_MIN,UTM_X_MIN,UTM_Y_MIN,ILONGLAT2UTM)
  call utm_geo(LONGITUDE_MAX,LATITUDE_MAX,UTM_X_MAX,UTM_Y_MAX,ILONGLAT2UTM)

  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

  ! check that parameters computed are consistent
  if (UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if (UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

  ! set time step and radial distribution of elements
  ! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)
  UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM

  !------------------------------------
  ! Mesh_Par_file parameter check
  !------------------------------------

  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do imat = 1,NMATERIALS
    mat_id = material_properties(imat,8)
    ! checks material id (must be positive or negative)
    if (mat_id == 0) stop 'Error incorrect material ID 0 found'
    ! counters
    if (mat_id > 0) ndef = ndef + 1
    if (mat_id < 0) nundef = nundef + 1
  enddo

  ! checks given material properties
  do imat = 1,NMATERIALS
    ! material properties
    rho = material_properties(imat,1)
    vp = material_properties(imat,2)
    vs = material_properties(imat,3)
    Q_Kappa = material_properties(imat,4)
    Q_mu = material_properties(imat,5)
    anisotropy_flag = material_properties(imat,6)
    domain_id = material_properties(imat,7)
    mat_id = material_properties(imat,8)

    ! checks material parameters
    if (rho <= 0.d0 .or. vp <= 0.d0 .or. vs < 0.d0) stop 'Error material with negative value of velocity or density found'
    if (Q_Kappa < 0.d0) stop 'Error material with negative Q_Kappa value found'
    if (Q_mu < 0.d0) stop 'Error material with negative Q_mu value found'

    ! checks domain id (1 = acoustic / 2 = elastic / 3 = poroelastic)
    select case (domain_id)
    case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC)
      continue
    case default
      stop 'Error invalid domain ID in materials'
    end select

    ! checks that material id doesn't exceed bounds of defined/undefined materials
    if (mat_id > 0 .and. mat_id > ndef) stop 'Error material ID bigger than total number of defined materials'
    if (mat_id < 0 .and. abs(mat_id) > nundef) stop 'Error negative material ID exceeds total number of undefined materials'
  enddo

  ! checks subregion ranges
  ner_max = 0
  do ireg = 1,NSUBREGIONS
    ix_beg_region = subregions(ireg,1)
    ix_end_region = subregions(ireg,2)
    iy_beg_region = subregions(ireg,3)
    iy_end_region = subregions(ireg,4)
    iz_beg_region = subregions(ireg,5)
    iz_end_region = subregions(ireg,6)
    imaterial_number = subregions(ireg,7)

    ! xi-direction range
    if (ix_beg_region < 1) stop 'XI coordinate of region negative!'
    if (ix_end_region > NEX_XI) stop 'XI coordinate of region too high!'

    ! eta-direction range
    if (iy_beg_region < 1) stop 'ETA coordinate of region negative!'
    if (iy_end_region > NEX_ETA) stop 'ETA coordinate of region too high!'

    ! depth range
    if (iz_beg_region < 1 .or. iz_end_region < 1) stop 'Z coordinate of region negative!'
    if (iz_end_region > ner_max) ner_max = iz_end_region

    if (imaterial_number == 0) stop 'Material ID of region zero!'
    ! searches material in given material_properties
    found = .false.
    do imat = 1,NMATERIALS
      mat_id = material_properties(imat,8)
      if (imaterial_number == mat_id) then
        found = .true.
        exit
      endif
    enddo
    if (.not. found) then
      print *,'Error: material id ',imaterial_number,' given in region ',ireg,' not found in materials section'
      stop 'Material ID of region not matching any given material'
    endif
  enddo
  ! checks if full range is covers
  if (minval(subregions(:,1)) > 1) stop 'Error minimum XI coordinate of regions must start at 1'
  if (maxval(subregions(:,2)) < NEX_XI) stop 'Error maximum XI coordinate of regions must end at NEX_XI'

  if (minval(subregions(:,3)) > 1) stop 'Error minimum ETA coordinate of regions must start at 1'
  if (maxval(subregions(:,4)) < NEX_ETA) stop 'Error maximum ETA coordinate of regions must end at NEX_ETA'

  if (minval(subregions(:,5)) > 1) stop 'Error minimum Z coordinate of regions must start at 1'
  ! maximum Z-layers must match with interface which is not read in here yet...
  !if (maxval(subregions(:,6)) < NER) stop 'Error maximum Z coordinate of regions must end at NER'

  ! checks doubling layer selection
  if (.not. USE_REGULAR_MESH) then
    do idoub = 1,NDOUBLINGS
      if (ner_doublings(idoub) < 1) then
        print *,'Error doubling layer ',idoub,' has invalid NZ value ',ner_doublings(idoub)
        stop 'Error doubling layer value must be at least 1'
      else if (ner_doublings(idoub) > ner_max) then
        print *,'Error doubling layer ',idoub,' with NZ value ',ner_doublings(idoub),'exceeds regions layers ',ner_max
        stop 'Error invalid doubling layer number, NZ exceeds regions NZ_END specification'
      endif
    enddo
  endif

  ! sorts doubling layer entries
  ! note: creating mesh layers starts from bottom, i.e. doubling builds layers from bottom to top,
  !       and a smaller NZ entry means the layer is closer to the bottom
  if (NDOUBLINGS > 1) then
    ! sorts with decreasing order, e.g. (14,10) instead of (10,14)
    ! switches doubling layer entries to have first doubling number with biggest and last entry with smallest value
    call bubble_sort_decreasing(NDOUBLINGS,ner_doublings)

    ! checks entries
    do idoub = 1,NDOUBLINGS-1
      ! checks if there are duplicate entries
      if (ner_doublings(idoub) == ner_doublings(idoub+1)) then
        print *,'Error doubling layer indices are invalid: ',ner_doublings(:)
        stop 'Error doubling layer entries contain duplicate values, please use unique layering'
      endif

      ! checks that entries are at least 2 layers apart
      ! (we alternate between regular and doubling subregions)
      if (ner_doublings(idoub) - ner_doublings(idoub+1) < 2 ) then
        print *,'Error doubling layer indices are too close: ',ner_doublings(idoub),' and ',ner_doublings(idoub+1)
        stop 'Error doubling layer entries must be at least 2 layers apart, please use larger layer spacings'
      endif
    enddo
  endif

  ! checks number of processes
  if (sizeprocs == 1 .and. (NPROC_XI /= 1 .or. NPROC_ETA /= 1)) &
    stop 'Error: must have NPROC_XI = NPROC_ETA = 1 for a serial run'

  ! make sure everybody is synchronized
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  all okay'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_mesh_parameter_file

!
!-------------------------------------------------------------------------------
!

  subroutine bubble_sort_decreasing(N,array)

! bubble sort routine
! orders integer array in decreasing order

  implicit none

  integer,intent(in) :: N
  integer,dimension(N),intent(inout) :: array

  ! local parameters
  integer :: i,j
  integer :: temp

  do i = 1,N
    do j = N, i+1, -1
      ! switches elements to have decreasing order
      if (array(j-1) < array(j)) then
        temp = array(j-1)
        array(j-1) = array(j)
        array(j) = temp
      endif
    enddo
  enddo

  end subroutine bubble_sort_decreasing
