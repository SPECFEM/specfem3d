
  program convert_CUBIT_SPECFEM_to_DX

! Dimitri Komatitsch, CNRS Marseille, France, March 2013

! convert CUBIT files that are in SPECFEM3D_Cartesian format to OpenDX format for visualization

  implicit none

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: iread,idummy
  integer :: i1,i2,i3,i4,i5,i6,i7,i8

  real, dimension(:), allocatable :: x,y,z

  real :: xread,yread,zread
  real :: val_color

! open SPECFEM3D_Cartesian mesh file to read the points
    open(unit=23,file='nodes_coords_file',status='old',action='read')
    read(23,*) npoin
    allocate(x(npoin))
    allocate(y(npoin))
    allocate(z(npoin))
    do ipoin = 1,npoin
      read(23,*) iread,xread,yread,zread
      x(iread) = xread
      y(iread) = yread
      z(iread) = zread
    enddo
    close(23)

! write DX header with element data
    open(unit=11,file='DX_fullmesh.dx',status='unknown',action='write')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

! read local points in this slice and output global DX points
  do ipoin=1,npoin
    write(11,*) x(ipoin),y(ipoin),z(ipoin)
  enddo

! ************* generate elements ******************

! open SPECFEM3D_Cartesian topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',nspec,' data follows'

! read local elements in this slice and output global DX elements
  do ispec=1,nspec
    read(23,*) idummy,i1,i2,i3,i4,i5,i6,i7,i8

! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS or SPECFEM
! and point numbers start at 0 rather than 1
    write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") i4-1,i1-1,i8-1,i5-1,i3-1,i2-1,i7-1,i6-1
  enddo

  close(23)

! ************* generate element data values ******************

! output DX header for data
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! read local elements in this slice and output global DX elements
  do ispec=1,nspec
    val_color = 1.0 ! dummy uniform color
    write(11,*) val_color
  enddo

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

  close(11)

  end program convert_CUBIT_SPECFEM_to_DX

