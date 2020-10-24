!=========================================================================================
module xdmf
  implicit none
  private
  public :: dump_mesh_xdmf
  public :: dump_mesh_data_xdmf

contains

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_xdmf(filename, npoints)
  character(len=*), intent(in)      :: filename
  integer, intent(in)               :: npoints

  integer                           :: iinput_xdmf, iinput_heavy_data
  character(len=512)                :: filename_np

  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 732) npoints, npoints, trim(filename_np), npoints, trim(filename_np)
  close(iinput_xdmf)

732 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/&
    '<Grid GridType="Uniform">',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="', i10, '">',/&
    '</Topology>',/&
    '<Geometry GeometryType="X_Y_Z">',/&
    '<DataItem DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '            ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '            ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</Geometry>',/&
    '</Grid>',/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dump_mesh_data_xdmf(filename, varname, npoints, nsnap)
  character(len=*), intent(in)      :: filename, varname
  integer, intent(in)               :: npoints, nsnap

  integer                           :: iinput_xdmf, iinput_heavy_data
  integer                           :: i
  character(len=512)                :: filename_np


  ! relative filename for xdmf content
  filename_np = trim(filename(index(filename, '/', back=.true.)+1:))

  ! XML Data
  open(newunit=iinput_xdmf, file=trim(filename)//'.xdmf')
  write(iinput_xdmf, 733) npoints, npoints, trim(filename_np), npoints, trim(filename_np)

  do i=1, nsnap
     ! create new snapshot in the temporal collection
     write(iinput_xdmf, 7341) dble(i), npoints, "'", "'"

     ! write attribute
     write(iinput_xdmf, 7342) varname, npoints, i-1, npoints, nsnap, npoints, &
                              trim(filename_np), trim(varname)

     write(iinput_xdmf, 7343)
  enddo

  ! finish xdmf file
  write(iinput_xdmf, 736)
  close(iinput_xdmf)

733 format(&
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="points" ItemType="function" function="join($0, $1)" Dimensions="', i10, ' 2">',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_S',/&
    '</DataItem>',/&
    '<DataItem Name="points" DataType="Float" Precision="8" Dimensions="', i10, '" Format="HDF">',/&
    '    ', A, ':/Mesh/mesh_Z',/&
    '</DataItem>',/&
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)

7341 format(&
    '<Grid Name="grid" GridType="Uniform">',/&
    '<Time Value="',F8.2,'" />',/&
    '<Topology TopologyType="Polyvertex" NumberOfElements="',i10,'">',/&
    '</Topology>',/&
    '<Geometry GeometryType="XY">',/&
    '<DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '</Geometry>')

7342 format(&
    '<Attribute Name="', A,'" AttributeType="Scalar" Center="Node">',/&
    '<DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '<DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '</DataItem>',/&
    '<DataItem DataType="Float" Precision="8" Dimensions="', i10, i10, '" Format="HDF">',/&
    '                    ', A, ':/', A, /&
    '</DataItem>',/,/&
    '</DataItem>',/&
    '</Attribute>')

7343 format(&
    '</Grid>',/)

736 format(&
    '</Grid>',/,/&
    '</Domain>',/&
    '</Xdmf>')

end subroutine
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================


program test_inversion_mesh
  use xdmf
  implicit none

  call dump_mesh_data_xdmf('axisem_output.nc4', 'Snapshots/straintrace', 6384, 32)

end program

