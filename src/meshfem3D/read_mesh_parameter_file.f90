!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
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
                                USE_REGULAR_MESH,NDOUBLINGS,ner_doublings)

  implicit none

  include "constants.h"

  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,UTM_PROJECTION_ZONE

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, UTM_MAX
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  logical SUPPRESS_UTM_PROJECTION,USE_REGULAR_MESH
  logical CREATE_ABAQUS_FILES,CREATE_DX_FILES,CREATE_VTK_FILES

  integer NDOUBLINGS
  integer, dimension(2) :: ner_doublings

  character(len=256) LOCAL_PATH
  character(len=50) INTERFACES_FILE

! local variables
  integer NEX_MAX

  double precision DEPTH_BLOCK_KM

  integer ierr

! subregions parameters
  integer NSUBREGIONS
  integer ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer iz_beg_region,iz_end_region,imaterial_number
!  definition of the different regions of the model in the mesh (nx,ny,nz)
!  #1 #2 : nx_begining,nx_end
!  #3 #4 : ny_begining,ny_end
!  #5 #6 : nz_begining,nz_end
!     #7 : material number
  integer, dimension(:,:), pointer :: subregions

! material properties
  integer NMATERIALS
  double precision rho,vp,vs,Q_flag,anisotropy_flag,domain_id
! first dimension  : material_id
! second dimension : #rho  #vp  #vs  #Q_flag  #anisotropy_flag #domain_id
  double precision, dimension(:,:), pointer :: material_properties

  integer i,ireg,imat,idoubl

! open parameter file Mesh_Par_file
  call open_parameter_file_mesh()

  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MIN, 'LATITUDE_MIN', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 1'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LATITUDE_MAX, 'LATITUDE_MAX', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 2'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MIN, 'LONGITUDE_MIN', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 3'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,LONGITUDE_MAX, 'LONGITUDE_MAX', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 4'
  call read_value_dble_precision_mesh(IIN,IGNORE_JUNK,DEPTH_BLOCK_KM, 'DEPTH_BLOCK_KM', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 5'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,UTM_PROJECTION_ZONE, 'UTM_PROJECTION_ZONE', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 6'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,SUPPRESS_UTM_PROJECTION, 'SUPPRESS_UTM_PROJECTION', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 7'

  call read_value_string_mesh(IIN,IGNORE_JUNK,INTERFACES_FILE, 'INTERFACES_FILE', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 8'

  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_XI, 'NEX_XI', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 9'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NEX_ETA, 'NEX_ETA', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 10'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_XI, 'NPROC_XI', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 11'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NPROC_ETA, 'NPROC_ETA', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 12'

! convert model size to UTM coordinates and depth of mesh to meters
  call utm_geo(LONGITUDE_MIN,LATITUDE_MIN,UTM_X_MIN,UTM_Y_MIN,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)
  call utm_geo(LONGITUDE_MAX,LATITUDE_MAX,UTM_X_MAX,UTM_Y_MAX,UTM_PROJECTION_ZONE,ILONGLAT2UTM,SUPPRESS_UTM_PROJECTION)

  Z_DEPTH_BLOCK = - dabs(DEPTH_BLOCK_KM) * 1000.d0

! check that parameters computed are consistent
  if(UTM_X_MIN >= UTM_X_MAX) stop 'horizontal dimension of UTM block incorrect'
  if(UTM_Y_MIN >= UTM_Y_MAX) stop 'vertical dimension of UTM block incorrect'

! set time step and radial distribution of elements
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)
  UTM_MAX = max(UTM_Y_MAX-UTM_Y_MIN, UTM_X_MAX-UTM_X_MIN)/1000.0 ! in KM

  call read_value_logical_mesh(IIN,IGNORE_JUNK,USE_REGULAR_MESH, 'USE_REGULAR_MESH', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 13'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NDOUBLINGS, 'NDOUBLINGS', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 14'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,ner_doublings(1), 'NZ_DOUGLING_1', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 15'
  call read_value_integer_mesh(IIN,IGNORE_JUNK,ner_doublings(2), 'NZ_DOUGLING_2', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 16'

  if(ner_doublings(1) < ner_doublings(2) .and. NDOUBLINGS == 2) then
    idoubl = ner_doublings(1)
    ner_doublings(1) = ner_doublings(2)
    ner_doublings(2) = idoubl
  endif

  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_ABAQUS_FILES, 'CREATE_ABAQUS_FILES', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 17'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_DX_FILES, 'CREATE_DX_FILES', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 18'
  call read_value_logical_mesh(IIN,IGNORE_JUNK,CREATE_VTK_FILES, 'CREATE_VTK_FILES', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 19'

! file in which we store the databases
  call read_value_string_mesh(IIN,IGNORE_JUNK,LOCAL_PATH, 'LOCAL_PATH', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 20'

! read number of materials
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NMATERIALS, 'NMATERIALS', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 21'

! read materials properties
  allocate(material_properties(NMATERIALS,6),stat=ierr)
  if (ierr /= 0) stop 'Allocation error of material_properties'

  do imat =1,NMATERIALS
     call read_material_parameters(IIN,i,rho,vp,vs,Q_flag,anisotropy_flag,domain_id,ierr)
     if (i /= imat) stop "Incorrect material ID"
     if(rho <= 0.d0 .or. vp <= 0.d0 .or. vs < 0.d0) stop 'negative value of velocity or density'
     material_properties(imat,1) = rho
     material_properties(imat,2) = vp
     material_properties(imat,3) = vs
     material_properties(imat,4) = Q_flag
     material_properties(imat,5) = anisotropy_flag
     material_properties(imat,6) = domain_id
  enddo

! read number of subregions
  call read_value_integer_mesh(IIN,IGNORE_JUNK,NSUBREGIONS, 'NSUBREGIONS', ierr)
  if (ierr /= 0) stop 'Error reading Mesh parameter 22'

! read subregions properties
  allocate(subregions(NSUBREGIONS,7),stat=ierr)
  if (ierr /= 0) stop 'Allocation error of subregions'
  do ireg =1,NSUBREGIONS
    call read_region_parameters(IIN,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region,&
         iz_beg_region,iz_end_region,imaterial_number,ierr)
    if (ix_beg_region < 1) stop 'XI coordinate of region negative!'
    if (ix_end_region > NEX_XI) stop 'XI coordinate of region too high!'
    if (iy_beg_region < 1) stop 'ETA coordinate of region negative!'
    if (iy_end_region > NEX_ETA) stop 'ETA coordinate of region too high!'
    if (iz_beg_region < 1) stop 'Z coordinate of region negative!'

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

  end subroutine read_mesh_parameter_file

end module readParFile
