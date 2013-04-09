
  program convert_mesh_to_CPML

! Dimitri Komatitsch, CNRS Marseille, France, April 2013

! convert outer layers of an existing CUBIT mesh file stored in SPECFEM3D_Cartesian format to CPML layers,
! i.e., create a 'absorbing_cpml_file' file for that existing mesh

  implicit none

! this is for HEX8; the code below also works for HEX27, it then just uses the first 8 points of an element
! to determine if it belongs to a CPML layer
  integer, parameter :: NGNOD = 8

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: ipoin_read,ispec_loop,iflag
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,number_of_CPML_elements

  real, dimension(:), allocatable :: x,y,z

  logical, dimension(:), allocatable :: is_X_CPML,is_Y_CPML,is_Z_CPML

  real :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax,limit

  logical :: ALSO_ADD_ON_THE_TOP_SURFACE

  real :: THICKNESS_OF_X_PML, THICKNESS_OF_Y_PML, THICKNESS_OF_Z_PML

! to make sure coordinate roundoff problems do not occur, use a tolerance of 1.5%
  real, parameter :: SMALL_PERCENTAGE_TOLERANCE = 1.015

! flags for the seven CPML regions
  integer, parameter :: CPML_X_ONLY = 1
  integer, parameter :: CPML_Y_ONLY = 2
  integer, parameter :: CPML_Z_ONLY = 3
  integer, parameter :: CPML_XY_ONLY = 4
  integer, parameter :: CPML_XZ_ONLY = 5
  integer, parameter :: CPML_YZ_ONLY = 6
  integer, parameter :: CPML_XYZ = 7

  print *
  print *,'IMPORTANT: it is your responsibility to make sure that in the input CUBIT (or similar) mesh'
  print *,'that this code will read in SPECFEM3D format from files "nodes_coords_file" and "mesh_file"'
  print *,'you have created layers of elements that constitute a layer of constant thickness aligned'
  print *,'with the coordinate grid axes (X, Y and/or Z), so that this code can assign CPML flags to them.'
  print *,'This code does NOT check that (because it cannot, in any easy way).'
  print *,'The mesh inside these CPML layers does not need to be structured nor regular, any non-structured'
  print *,'mesh is fine as long as it has flat PML inner and outer faces, parallel to the axes, and thus'
  print *,'of a constant thickness.'
  print *,'The thickness can be different between the X, Y and Z sides. But for X it must not vary,'
  print *,'for Y it must not vary, and for Z it must not vary.'
  print *,'If you do not know the exact thickness, you can use a slightly LARGER value'
  print *,'in this code (say 2% to 5% more) and this code will fix that and will adjust it;'
  print *,'never use a SMALLER value otherwise this code will miss some CPML elements.'
  print *
  print *,'Note: in a future release we will remove the constraint of having CPML layers aligned with the'
  print *,'coordinate axes; we will allow for meshes that are titled by any constant angle in the horizontal plane.'
  print *,'However this is not implemented yet.'
  print *

  print *,'1 = use a free surface at the top of the mesh (most classical option)'
  print *,'2 = use a CPML absorbing layer at the top of the mesh (less classical option)'
  print *,'3 = exit'
  read(*,*) iflag
  if(iflag /= 1 .and. iflag /= 2) stop 'exiting...'
  if(iflag == 1) then
    ALSO_ADD_ON_THE_TOP_SURFACE = .false.
  else
    ALSO_ADD_ON_THE_TOP_SURFACE = .true.
  endif
  print *

! open SPECFEM3D_Cartesian mesh file to read the points
    open(unit=23,file='nodes_coords_file',status='old',action='read')
    read(23,*) npoin
    allocate(x(npoin))
    allocate(y(npoin))
    allocate(z(npoin))
    do ipoin = 1,npoin
      read(23,*) ipoin_read,xread,yread,zread
      x(ipoin_read) = xread
      y(ipoin_read) = yread
      z(ipoin_read) = zread
    enddo
    close(23)

! compute the min and max values of each coordinate
   xmin = minval(x)
   xmax = maxval(x)

   ymin = minval(y)
   ymax = maxval(y)

   zmin = minval(z)
   zmax = maxval(z)

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Ymin and Ymax of the mesh read = ',ymin,ymax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

  print *,'What is the exact thickness of the PML layer that you want'
  print *,'on the Xmin and Xmax faces of your mesh? (it needs to correspond exactly'
  print *,'to the flat layers you created in your input CUBIT mesh, as mentioned in'
  print *,'the comment printed above; if you think you have roundoff issues or very'
  print *,'slightly varying thickness, give 2% or 5% more here, but never less'
  read(*,*) THICKNESS_OF_X_PML
  if(THICKNESS_OF_X_PML <= 0) stop 'negative thickness is not allowed; exiting...'
  if(THICKNESS_OF_X_PML > 0.30*(xmax - xmin)) &
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  print *

  print *,'What is the exact thickness of the PML layer that you want'
  print *,'on the Ymin and Ymax faces of your mesh?'
  read(*,*) THICKNESS_OF_Y_PML
  if(THICKNESS_OF_Y_PML <= 0) stop 'negative thickness is not allowed; exiting...'
  if(THICKNESS_OF_Y_PML > 0.30*(ymax - ymin)) &
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  print *

  print *,'What is the exact thickness of the PML layer that you want'
  print *,'on the Zmin and Zmax faces of your mesh?'
  read(*,*) THICKNESS_OF_Z_PML
  if(THICKNESS_OF_Z_PML <= 0) stop 'negative thickness is not allowed; exiting...'
  if(THICKNESS_OF_Z_PML > 0.30*(zmax - zmin)) &
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  print *

! ************* read mesh elements and generate CPML flags *************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  allocate(is_X_CPML(nspec))
  allocate(is_Y_CPML(nspec))
  allocate(is_Z_CPML(nspec))

  is_X_CPML(:) = .false.
  is_Y_CPML(:) = .false.
  is_Z_CPML(:) = .false.

! loop on the whole mesh
  do ispec_loop = 1,nspec

    read(23,*) ispec,i1,i2,i3,i4,i5,i6,i7,i8

! Xmin CPML
    limit = xmin + THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
    if(x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit .and. &
       x(i5) < limit .and. x(i6) < limit .and. x(i7) < limit .and. x(i8) < limit) is_X_CPML(ispec) = .true.

! Xmax CPML
    limit = xmax - THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
    if(x(i1) > limit .and. x(i2) > limit .and. x(i3) > limit .and. x(i4) > limit .and. &
       x(i5) > limit .and. x(i6) > limit .and. x(i7) > limit .and. x(i8) > limit) is_X_CPML(ispec) = .true.

! Ymin CPML
    limit = ymin + THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
    if(y(i1) < limit .and. y(i2) < limit .and. y(i3) < limit .and. y(i4) < limit .and. &
       y(i5) < limit .and. y(i6) < limit .and. y(i7) < limit .and. y(i8) < limit) is_Y_CPML(ispec) = .true.

! Ymax CPML
    limit = ymax - THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
    if(y(i1) > limit .and. y(i2) > limit .and. y(i3) > limit .and. y(i4) > limit .and. &
       y(i5) > limit .and. y(i6) > limit .and. y(i7) > limit .and. y(i8) > limit) is_Y_CPML(ispec) = .true.

! Zmin CPML
    limit = zmin + THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
    if(z(i1) < limit .and. z(i2) < limit .and. z(i3) < limit .and. z(i4) < limit .and. &
       z(i5) < limit .and. z(i6) < limit .and. z(i7) < limit .and. z(i8) < limit) is_Z_CPML(ispec) = .true.

! Zmax CPML
  if(ALSO_ADD_ON_THE_TOP_SURFACE) then
    limit = zmax - THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
    if(z(i1) > limit .and. z(i2) > limit .and. z(i3) > limit .and. z(i4) > limit .and. &
       z(i5) > limit .and. z(i6) > limit .and. z(i7) > limit .and. z(i8) > limit) is_Z_CPML(ispec) = .true.
  endif

  enddo

  close(23)

  print *,'Total number of elements in the mesh = ',nspec
  print *
  print *,'Found ',count(is_X_CPML),' X_CPML elements'
  print *,'Found ',count(is_Y_CPML),' Y_CPML elements'
  print *,'Found ',count(is_Z_CPML),' Z_CPML elements'
  if(ALSO_ADD_ON_THE_TOP_SURFACE) print *,'    (also converted the top surface from free surface to CPML)'
  print *

  number_of_CPML_elements = 0
  do ispec=1,nspec
    if(is_X_CPML(ispec) .or. is_Y_CPML(ispec) .or. is_Z_CPML(ispec)) &
          number_of_CPML_elements = number_of_CPML_elements + 1
  enddo
  print *,'Created a total of ',number_of_CPML_elements,' unique CPML elements'
  print *,'   (i.e., ',100.*number_of_CPML_elements/real(nspec),'% of the mesh)'

! ************* generate the CPML database file *************

  open(unit=23,file='absorbing_cpml_file',status='unknown',action='write')

! write the total number of unique CPML elements
  write(23,*) number_of_CPML_elements

! write the CPML flag for each CPML element
  do ispec=1,nspec
    if(is_X_CPML(ispec) .and. is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,CPML_XYZ

    else if(is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,CPML_YZ_ONLY

    else if(is_X_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,CPML_XZ_ONLY

    else if(is_X_CPML(ispec) .and. is_Y_CPML(ispec)) then
      write(23,*) ispec,CPML_XY_ONLY

    else if(is_Z_CPML(ispec)) then
      write(23,*) ispec,CPML_Z_ONLY

    else if(is_Y_CPML(ispec)) then
      write(23,*) ispec,CPML_Y_ONLY

    else if(is_X_CPML(ispec)) then
      write(23,*) ispec,CPML_X_ONLY
    endif

  enddo
  close(23)

  print *
  print *,'CPML absorbing layer file "absorbing_cpml_file" has been successfully created'
  print *

  end program convert_mesh_to_CPML

