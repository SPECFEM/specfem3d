!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
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

module readParFile

contains

  subroutine read_mesh_parameter_file(LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX, &
                                      UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
                                      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE, &
                                      LOCAL_PATH,SUPPRESS_UTM_PROJECTION,&
                                      INTERFACES_FILE,NSUBREGIONS,subregions,NMATERIALS,material_properties,&
                                      CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES, &
                                      USE_REGULAR_MESH,NDOUBLINGS,ner_doublings, &
                                      THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML)

  use constants

  implicit none

  double precision,intent(out) :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX
  double precision,intent(out) :: UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK

  integer,intent(out) :: NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA

  integer,intent(out) :: UTM_PROJECTION_ZONE
  logical,intent(out) :: SUPPRESS_UTM_PROJECTION

  character(len=MAX_STRING_LEN),intent(out) :: LOCAL_PATH
  character(len=MAX_STRING_LEN),intent(out) :: INTERFACES_FILE

  ! subregions parameters
  integer,intent(out) :: NSUBREGIONS
  !  definition of the different regions of the model in the mesh (nx,ny,nz)
  !  #1 #2 : nx_begining,nx_end
  !  #3 #4 : ny_begining,ny_end
  !  #5 #6 : nz_begining,nz_end
  !     #7 : material number
  integer, dimension(:,:), pointer,intent(out) :: subregions

  ! material properties
  integer,intent(out) :: NMATERIALS
  ! first dimension  : local number of material
  ! (note: material_id can be positive/negative for defined/undefined materials)
  ! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id #material_id
  double precision, dimension(:,:), pointer,intent(out) :: material_properties

  logical,intent(out) :: CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

  logical,intent(out) :: USE_REGULAR_MESH
  integer,intent(out) :: NDOUBLINGS
  integer, dimension(2),intent(out) :: ner_doublings

  ! CPML
  double precision, intent(out) :: THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML

! local variables
  integer :: NEX_MAX
  double precision :: UTM_MAX
  double precision :: DEPTH_BLOCK_KM
  integer :: ier
  integer :: ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer :: iz_beg_region,iz_end_region,imaterial_number
  double precision :: rho,vp,vs,Q_flag,anisotropy_flag
  integer :: ireg,imat,idoubl,ndef,nundef
  integer :: mat_id,domain_id
  logical :: found

! open parameter file Mesh_Par_file
  call open_parameter_file_mesh()

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
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_XI, 'NEX_XI', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NEX_XI'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_ETA, 'NEX_ETA', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NEX_ETA'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_XI, 'NPROC_XI', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NPROC_XI'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_ETA, 'NPROC_ETA', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NPROC_ETA'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,USE_REGULAR_MESH, 'USE_REGULAR_MESH', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter USE_REGULAR_MESH'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NDOUBLINGS, 'NDOUBLINGS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NDOUBLINGS'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,ner_doublings(1), 'NZ_DOUGLING_1', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUGLING_1'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,ner_doublings(2), 'NZ_DOUGLING_2', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NZ_DOUGLING_2'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_ABAQUS_FILES, 'CREATE_ABAQUS_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_ABAQUS_FILES'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_DX_FILES, 'CREATE_DX_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_DX_FILES'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_VTK_FILES, 'CREATE_VTK_FILES', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter CREATE_VTK_FILES'

  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_X_PML, 'THICKNESS_OF_X_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_X_PML'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Y_PML, 'THICKNESS_OF_Y_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Y_PML'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,THICKNESS_OF_Z_PML, 'THICKNESS_OF_Z_PML', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter THICKNESS_OF_Z_PML'

  ! file in which we store the databases
  call read_value_string_mesh(IIN,IGNORE_JUNK,LOCAL_PATH, 'LOCAL_PATH', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter LOCAL_PATH'

  ! read number of materials
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NMATERIALS, 'NMATERIALS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NMATERIALS'
  ! read materials properties
  allocate(material_properties(NMATERIALS,7),stat=ier)
  if (ier /= 0) stop 'Error allocation of material_properties'
  material_properties(:,:) = 0.d0
  do imat = 1,NMATERIALS
    call read_material_parameters(IIN,mat_id,rho,vp,vs,Q_flag,anisotropy_flag,domain_id,ier)
    if (ier /= 0) stop 'Error reading materials in Mesh_Par_file'
    ! stores material
    material_properties(imat,1) = rho
    material_properties(imat,2) = vp
    material_properties(imat,3) = vs
    material_properties(imat,4) = Q_flag
    material_properties(imat,5) = anisotropy_flag
    material_properties(imat,6) = domain_id
    material_properties(imat,7) = mat_id
  enddo

  ! read number of subregions
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NSUBREGIONS, 'NSUBREGIONS', ier)
  if (ier /= 0) stop 'Error reading Mesh parameter NSUBREGIONS'
  ! read subregions properties
  allocate(subregions(NSUBREGIONS,7),stat=ier)
  if (ier /= 0) stop 'Error allocation of subregions'
  subregions(:,:) = 0
  do ireg =1,NSUBREGIONS
    call read_region_parameters(IIN,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
                                iz_beg_region,iz_end_region,imaterial_number,ier)
    if (ier /= 0) stop 'Error reading regions in Mesh_Par_file'
    ! stores region
    subregions(ireg,1) = ix_beg_region
    subregions(ireg,2) = ix_end_region
    subregions(ireg,3) = iy_beg_region
    subregions(ireg,4) = iy_end_region
    subregions(ireg,5) = iz_beg_region
    subregions(ireg,6) = iz_end_region
    subregions(ireg,7) = imaterial_number
  enddo

  ! close parameter file
  call close_parameter_file_mesh()

  ! convert model size to UTM coordinates and depth of mesh to meters
  call utm_geo(LONGITUDE_MIN,LATITUDE_MIN,UTM_X_MIN,UTM_Y_MIN,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
  call utm_geo(LONGITUDE_MAX,LATITUDE_MAX,UTM_X_MAX,UTM_Y_MAX,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

  ! check that parameters computed are consistent
  if (UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if (UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

  ! set time step and radial distribution of elements
  ! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)
  UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM

  ! switches number of doubling layers to have first doubling number with bigger value
  ! (doubling counts layers from bottom to top)
  if (NDOUBLINGS == 2) then
    if (ner_doublings(1) < ner_doublings(2)) then
      idoubl = ner_doublings(1)
      ner_doublings(1) = ner_doublings(2)
      ner_doublings(2) = idoubl
    endif
  endif

  !------------------------------------
  ! Mesh_Par_file parameter check
  !------------------------------------

  ! counts defined/undefined materials
  ndef = 0
  nundef = 0
  do imat = 1,NMATERIALS
    mat_id = material_properties(imat,7)
    ! checks material id (must be positiv or negative)
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
    Q_flag = material_properties(imat,4)
    anisotropy_flag = material_properties(imat,5)
    domain_id = material_properties(imat,6)
    mat_id = material_properties(imat,7)

    ! checks material parameters
    if (rho <= 0.d0 .or. vp <= 0.d0 .or. vs < 0.d0) stop 'Error material with negative value of velocity or density found'
    if (Q_flag < 0.d0) stop 'Error material with negative Q value found'

    ! checks domain id (1 = acoustic / 2 = elastic / 3 = poroelastic)
    select case (domain_id)
    case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC)
      continue
    case (IDOMAIN_POROELASTIC)
      stop 'Error materials for poro-elastic domains not implemented yet'
    case default
      stop 'Error invalid domain ID in materials'
    end select

    ! checks that material id doesn't exceed bounds of defined/undefined materials
    if (mat_id > 0 .and. mat_id > ndef) stop 'Error material ID bigger than total number of defined materials'
    if (mat_id < 0 .and. abs(mat_id) > nundef) stop 'Error negative material ID exceeds total number of undefined materials'
  enddo

  ! checks subregion ranges
  do ireg =1,NSUBREGIONS
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
    if (iz_beg_region < 1) stop 'Z coordinate of region negative!'

    if (imaterial_number == 0) stop 'Material ID of region zero!'
    ! searches material in given material_properties
    found = .false.
    do imat = 1,NMATERIALS
      mat_id = material_properties(imat,7)
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

  end subroutine read_mesh_parameter_file

end module readParFile
