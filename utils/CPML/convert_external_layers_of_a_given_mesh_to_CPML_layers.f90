
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
  integer :: ipoin_read,ispec_read
  integer :: i1,i2,i3,i4,i5,i6,i7,i8,number_of_CPML_elements

  real, dimension(:), allocatable :: x,y,z

  integer, dimension(:,:), allocatable :: ibool

  logical, dimension(:), allocatable :: is_X_CPML,is_Y_CPML,is_Z_CPML

  real :: xread,yread,zread,xmin,xmax,ymin,ymax,zmin,zmax,limit

  logical, parameter :: ALSO_ADD_ON_THE_TOP_SURFACE = .false.

  real, parameter :: THICKNESS_OF_X_PML = 3722.219 * 3 !! DK DK because we want 3 CPML elements
  real, parameter :: THICKNESS_OF_Y_PML = THICKNESS_OF_X_PML
  real, parameter :: THICKNESS_OF_Z_PML = 3750. * 3 !! DK DK because we want 3 CPML elements

  real, parameter :: SMALL_PERCENTAGE_TOLERANCE = 1.03 !! to make sure coordinate roundoff problems do not occur

!! DK DK dire dans manuel et ici en comment: it is your responsibility to have flat mesh elements of the right thickness aligned
!! with the coordinate axes in your input CUBIT (or similar) mesh; this code does NOT check that

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

! ************* generate elements ******************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  allocate(ibool(NGNOD,nspec))

  allocate(is_X_CPML(nspec))
  allocate(is_Y_CPML(nspec))
  allocate(is_Z_CPML(nspec))

! read local elements in this slice and output global DX elements
  do ispec=1,nspec
    read(23,*) ispec_read,i1,i2,i3,i4,i5,i6,i7,i8
    ibool(1,ispec_read) = i1
    ibool(2,ispec_read) = i2
    ibool(3,ispec_read) = i3
    ibool(4,ispec_read) = i4
    ibool(5,ispec_read) = i5
    ibool(6,ispec_read) = i6
    ibool(7,ispec_read) = i7
    ibool(8,ispec_read) = i8
  enddo
  close(23)

! ************* generate CPML flags *************

  is_X_CPML(:) = .false.
  is_Y_CPML(:) = .false.
  is_Z_CPML(:) = .false.

! loop on the whole mesh
  do ispec=1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)
    i5 = ibool(5,ispec)
    i6 = ibool(6,ispec)
    i7 = ibool(7,ispec)
    i8 = ibool(8,ispec)

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

! write the encoded CPML flag for each CPML element
  do ispec=1,nspec
    if(is_X_CPML(ispec) .and. is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,' 7'

    else if(is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,' 6'

    else if(is_X_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(23,*) ispec,' 5'

    else if(is_X_CPML(ispec) .and. is_Y_CPML(ispec)) then
      write(23,*) ispec,' 4'

    else if(is_Z_CPML(ispec)) then
      write(23,*) ispec,' 3'

    else if(is_Y_CPML(ispec)) then
      write(23,*) ispec,' 2'

    else if(is_X_CPML(ispec)) then
      write(23,*) ispec,' 1'
    endif

  enddo
  close(23)

  end program convert_mesh_to_CPML

