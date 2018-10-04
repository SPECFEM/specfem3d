
  program mesh_a_Y_shaped_crack

  implicit none

! points defining the different regions to mesh
  double precision, dimension(31) :: xnpgeo,znpgeo

  integer :: ip1,ip2,ip3,ip4,nx,nz,ispec,ier

  xnpgeo(01) = 0.d0
  znpgeo(01) = 0.d0

  xnpgeo(02) = 18.6
  znpgeo(02) = 0.d0

  xnpgeo(03) = 50
  znpgeo(03) = 0.d0

  xnpgeo(04) = 0.d0
  znpgeo(04) = 15.d0

  xnpgeo(05) = 18.6
  znpgeo(05) = 15.d0

  xnpgeo(06) = 50
  znpgeo(06) = 15.d0

  xnpgeo(07) = 0.d0
  znpgeo(07) = 19.87

  xnpgeo(08) = 18.6
  znpgeo(08) = 19.87

  xnpgeo(09) = 50
  znpgeo(09) = 19.87

  xnpgeo(10) = 0
  znpgeo(10) = 25.37

  xnpgeo(11) = 18.6
  znpgeo(11) = 25.37

  xnpgeo(12) = 23.36
  znpgeo(12) = 23.53

  xnpgeo(13) = 50
  znpgeo(13) = 23.53

  xnpgeo(14) = 0
  znpgeo(14) = 29.46

  xnpgeo(15) = 23.95
  znpgeo(15) = 29.46

  xnpgeo(16) = 50
  znpgeo(16) = 29.46

  xnpgeo(17) = 0
  znpgeo(17) = 35

  xnpgeo(18) = 26.697
  znpgeo(18) = 35

  xnpgeo(19) = 50
  znpgeo(19) = 35

  xnpgeo(20) = 0
  znpgeo(20) = 43

  xnpgeo(21) = 30.85
  znpgeo(21) = 43

  xnpgeo(22) = 50
  znpgeo(22) = 43

  xnpgeo(23) = 0
  znpgeo(23) = 50

  xnpgeo(24) = 30.85
  znpgeo(24) = 50

  xnpgeo(25) = 50
  znpgeo(25) = 50

  xnpgeo(26) = 0
  znpgeo(26) = 65

  xnpgeo(27) = 30.85
  znpgeo(27) = 65

  xnpgeo(28) = 50
  znpgeo(28) = 65

  xnpgeo(29) = 0
  znpgeo(29) = 85

  xnpgeo(30) = 30.85
  znpgeo(30) = 85

  xnpgeo(31) = 50
  znpgeo(31) = 85

! convert to millimiters
  xnpgeo(:) = xnpgeo(:) / 1000.d0
  znpgeo(:) = znpgeo(:) / 1000.d0

  ! user output
  write(*,*)
  write(*,*) 'Saving the mesh in Gnuplot format...'
  write(*,*)

  open(unit=20,file='gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) stop 'Error opening gnuplot file for writing: gridfile.gnu'

! generate the mesh for all the regions
  ispec = 0

! region 01
  ip1 = 1
  ip2 = 2
  ip3 = 5
  ip4 = 4
  nx = 32
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 02
  ip1 = 2
  ip2 = 3
  ip3 = 6
  ip4 = 5
  nx = 32
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 03
  ip1 = 4
  ip2 = 5
  ip3 = 8
  ip4 = 7
  nx = 32
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 04
  ip1 = 5
  ip2 = 6
  ip3 = 9
  ip4 = 8
  nx = 32
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 05
  ip1 = 7
  ip2 = 8
  ip3 = 11
  ip4 = 10
  nx = 32
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 06
  ip1 = 8
  ip2 = 9
  ip3 = 13
  ip4 = 12
  nx = 32
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 07
  ip1 = 8
  ip2 = 12
  ip3 = 15
  ip4 = 11
  nx = 6
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 08
  ip1 = 10
  ip2 = 11
  ip3 = 15
  ip4 = 14
  nx = 32
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 09
  ip1 = 12
  ip2 = 13
  ip3 = 16
  ip4 = 15
  nx = 32
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 10
  ip1 = 14
  ip2 = 15
  ip3 = 18
  ip4 = 17
  nx = 32
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 11
  ip1 = 15
  ip2 = 16
  ip3 = 19
  ip4 = 18
  nx = 32
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 12
  ip1 = 17
  ip2 = 18
  ip3 = 21
  ip4 = 20
  nx = 32
  nz = 10
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 13
  ip1 = 18
  ip2 = 19
  ip3 = 22
  ip4 = 21
  nx = 32
  nz = 10
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 14
  ip1 = 20
  ip2 = 21
  ip3 = 24
  ip4 = 23
  nx = 32
  nz = 9
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 15
  ip1 = 21
  ip2 = 22
  ip3 = 25
  ip4 = 24
  nx = 32
  nz = 9
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 16
  ip1 = 23
  ip2 = 24
  ip3 = 27
  ip4 = 26
  nx = 32
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 17
  ip1 = 24
  ip2 = 25
  ip3 = 28
  ip4 = 27
  nx = 32
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 18
  ip1 = 26
  ip2 = 27
  ip3 = 30
  ip4 = 29
  nx = 32
  nz = 26
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

! region 19
  ip1 = 27
  ip2 = 28
  ip3 = 31
  ip4 = 30
  nx = 32
  nz = 26
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

  print *,'total number of mesh elements generated = ',ispec
  if (ispec /= (32 + 32) * (26*2 + 19*3) + 6*7) stop 'error in total number of mesh elements created'

  close(20)

  ! create a Gnuplot script to display the grid
  open(unit=20,file='plot_gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) stop 'Error saving plotgnu file'

  write(20,*) '#set term wxt'
  write(20,*)
  write(20,*) 'set term postscript portrait color solid "Helvetica" 22'
  write(20,*) 'set output "maillage_du_crack_OK_avec_petit_mailleur_ecrit_par_Dimitri.ps"'
  write(20,*)

  ! use same unit length on both X and Y axes
  write(20,*) 'set size ratio -1'

  ! size of our model
  write(20,*) 'set xrange [0:0.050]'
  write(20,*) 'set yrange [0:0.085]'

  ! draw rectangles showing the water and steel layers
  write(20,*) 'set object 1 rect from 0,0.065 to 0.050,0.085 fc rgb "#99FFFF" back'
  write(20,*) 'set object 2 rect from 0,0.035 to 0.050,0.050 fc rgb "#888888" back'
  write(20,*) 'set object 3 rect from 0,0 to 0.050,0.015 fc rgb "#888888" back'

  write(20,*) 'plot "gridfile.gnu" title "" w l lc black'
  write(20,*) 'pause -1 "Hit any key..."'
  close(20)

  write(*,*) 'Mesh saved in Gnuplot format...'
  write(*,*)

  end program mesh_a_Y_shaped_crack

!---------

  subroutine generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,xnpgeo,znpgeo)

  implicit none

  integer :: ip1,ip2,ip3,ip4,nx,nz,ispec

! points defining the different regions to mesh
  double precision, dimension(31) :: xnpgeo,znpgeo

  double precision :: ratio_ix,ratio_iz,ratio_ixplus1,ratio_izplus1
  double precision :: x1,z1,x2,z2,x3,z3,x4,z4

  integer :: ix,iz

  double precision, dimension(:,:), allocatable :: x,z

  allocate(x(0:nx,0:nz))
  allocate(z(0:nx,0:nz))

  do iz = 0,nz-1
    do ix = 0,nx-1

! generate one more mesh element
      ispec = ispec + 1

      ratio_ix = ix / dble(nx)
      ratio_iz = iz / dble(nz)

      ratio_ixplus1 = (ix+1) / dble(nx)
      ratio_izplus1 = (iz+1) / dble(nz)

! point 1
      call interpolate_bilinear(x1,z1,ip1,ip2,ip3,ip4,xnpgeo,znpgeo,ratio_ix,ratio_iz)

! point 2
      call interpolate_bilinear(x2,z2,ip1,ip2,ip3,ip4,xnpgeo,znpgeo,ratio_ixplus1,ratio_iz)

! point 3
      call interpolate_bilinear(x3,z3,ip1,ip2,ip3,ip4,xnpgeo,znpgeo,ratio_ixplus1,ratio_izplus1)

! point 4
      call interpolate_bilinear(x4,z4,ip1,ip2,ip3,ip4,xnpgeo,znpgeo,ratio_ix,ratio_izplus1)

! save the points created
      x(ix,iz) = x1
      z(ix,iz) = z1

      x(ix+1,iz) = x2
      z(ix+1,iz) = z2

      x(ix+1,iz+1) = x3
      z(ix+1,iz+1) = z3

      x(ix,iz+1) = x4
      z(ix,iz+1) = z4

    enddo
  enddo

  call save_gnuplot_file(nx,nz,x,z)

  deallocate(x)
  deallocate(z)

  end subroutine generate_region

!---------

  subroutine interpolate_bilinear(x,z,ip1,ip2,ip3,ip4,xnpgeo,znpgeo,gamma_interp_x,gamma_interp_z)

  implicit none

  integer :: ip1,ip2,ip3,ip4

! points defining the different regions to mesh
  double precision, dimension(31) :: xnpgeo,znpgeo

  double precision :: gamma_interp_x,gamma_interp_z
  double precision :: x,z
  double precision :: val1,val2,val3,val4

  ! interpolation rule
  val1 = xnpgeo(ip1)
  val2 = xnpgeo(ip2)
  val3 = xnpgeo(ip3)
  val4 = xnpgeo(ip4)
  x =  val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_z        + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_z

  val1 = znpgeo(ip1)
  val2 = znpgeo(ip2)
  val3 = znpgeo(ip3)
  val4 = znpgeo(ip4)
  z =  val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_z        + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_z

  end subroutine interpolate_bilinear

!-------------------------

!!! this is src/meshfem2D/save_gnuplot_file.f90 from SPECFEM2D to display the mesh in Gnuplot format

  subroutine save_gnuplot_file(nx,nz,x,z)

! creates a Gnuplot file that displays the grid

  implicit none

  integer :: nx,nz
  double precision, dimension(0:nx,0:nz) :: x,z

  ! local parameters
  integer :: ier,ili,icol

  ! draw horizontal lines of the grid
  do ili=0,nz
    do icol=0,nx-1
       write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(20,*) sngl(x(icol+1,ili)),sngl(z(icol+1,ili))
       write(20,10)
    enddo
  enddo

  ! draw vertical lines of the grid
  do icol=0,nx
    do ili=0,nz-1
       write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(20,*) sngl(x(icol,ili+1)),sngl(z(icol,ili+1))
       write(20,10)
    enddo
  enddo

10   format('')

  end subroutine save_gnuplot_file

