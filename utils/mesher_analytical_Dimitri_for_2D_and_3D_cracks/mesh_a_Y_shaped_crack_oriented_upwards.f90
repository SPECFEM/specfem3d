
  program mesh_a_Y_shaped_crack

  implicit none

! spatial dimension of the mesh
  integer, parameter :: NDIM = 2

! total number of spectral elements to create
!! DK DK half of the elements in the NX direction
  integer, parameter :: nx_base = 48 ! 32 ! DK DK I have densified by 1.5
!! DK DK number of elements in the vertical direction in the water layer
  integer, parameter :: nz_water = 39 ! 26 ! DK DK I have densified by 1.5
  integer, parameter :: ngnod = 4
  integer, parameter :: nspec = (nx_base + nx_base) * (19 + 19 + 7 + 6 + 7 + 6 + 3 + 19 + nz_water) + 6*7
  integer, parameter :: NGLOB_MAX = nspec * ngnod

! number of mesh elements in the third direction of the mesh
! (direction in which we extend the mesh identical to itself to convert it from 2D to 3D)
  integer, parameter :: ny = nx_base + nx_base
  double precision, parameter :: Ymin = 0.d0, Ymax = 0.050d0
  double precision, parameter :: deltaY = (Ymax - Ymin) / dble(ny)

! points defining the different regions to mesh
  double precision, dimension(31) :: xnpgeo,znpgeo

! introduce split nodes on the crack or not (if not, the crack becomes invisible i.e. it does not exist any more
! and the layers that contain it become homogeneous and continuous; this is useful to provide a reference solution to compare to)
  logical, parameter :: INTRODUCE_SPLIT_NODES_ON_THE_CRACK_CLEAN = .false.

  double precision, dimension(NGLOB_MAX) :: x,z
  integer, dimension(nspec) :: imat
  integer, dimension(ngnod,nspec) :: ibool

  integer :: ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,material_number,region_number,ipoin3D,ispec3D,iy
!!!!!!!!!!!!  integer :: ispec_free_surface_found
  integer :: ispec_absorbing_surface_found
  double precision :: y
  double precision :: geomtol
  double precision :: big_X_offset
  double precision :: angle_crack_to_horizontal,x1,z1,x2,z2,x3,z3

! additional arrays needed to sort the mesh points and remove the multiples
  integer, dimension(:), allocatable :: locval,ind,ninseg,iglob,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,work

  integer :: nseg,ioff,iseg,ig,i,j,nglob,ia
  integer :: ieoff,ilocnum

  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist

! flags for Stacey absorbing boundaries
  integer, parameter :: IBOTTOM = 1
  integer, parameter :: IRIGHT = 2
  integer, parameter :: ITOP = 3
  integer, parameter :: ILEFT = 4

! very large value
  double precision, parameter :: HUGEVAL = 1.d+30

! relative mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 1.d-10

  double precision, parameter :: PI = 3.141592653589793d0

  xnpgeo(01) = 0.d0
  znpgeo(01) = 0.d0

  xnpgeo(02) = 21.28
  znpgeo(02) = 0.d0

  xnpgeo(03) = 50
  znpgeo(03) = 0.d0

  xnpgeo(04) = 0.d0
  znpgeo(04) = 15.d0

  xnpgeo(05) = 21.28
  znpgeo(05) = 15.d0

  xnpgeo(06) = 50
  znpgeo(06) = 15.d0

  xnpgeo(07) = 0.d0
  znpgeo(07) = 28.53

  xnpgeo(08) = 21.28
  znpgeo(08) = 28.53

  xnpgeo(09) = 50
  znpgeo(09) = 28.53

  xnpgeo(10) = 0
  znpgeo(10) = 35.

  xnpgeo(11) = 24.75
  znpgeo(11) = 35.

  xnpgeo(12) = 50.
  znpgeo(12) = 35.

  xnpgeo(13) = 0
  znpgeo(13) = 38.82

  xnpgeo(14) = 27.09
  znpgeo(14) = 38.82

  xnpgeo(15) = 50.
  znpgeo(15) = 38.82

  xnpgeo(16) = 0
  znpgeo(16) = 44.7

  xnpgeo(17) = 27.7
  znpgeo(17) = 44.7

  xnpgeo(18) = 31.87
  znpgeo(18) = 42.75

  xnpgeo(19) = 50
  znpgeo(19) = 42.75

  xnpgeo(20) = 0
  znpgeo(20) = 47.89

  xnpgeo(21) = 30.35
  znpgeo(21) = 47.89

  xnpgeo(22) = 50
  znpgeo(22) = 47.89

  xnpgeo(23) = 0
  znpgeo(23) = 50

  xnpgeo(24) = 30.35
  znpgeo(24) = 50

  xnpgeo(25) = 50
  znpgeo(25) = 50

  xnpgeo(26) = 0
  znpgeo(26) = 65

  xnpgeo(27) = 30.35
  znpgeo(27) = 65

  xnpgeo(28) = 50
  znpgeo(28) = 65

  xnpgeo(29) = 0
  znpgeo(29) = 85

  xnpgeo(30) = 30.35
  znpgeo(30) = 85

  xnpgeo(31) = 50
  znpgeo(31) = 85

! convert to millimiters
  xnpgeo(:) = xnpgeo(:) / 1000.d0
  znpgeo(:) = znpgeo(:) / 1000.d0

! generate the mesh for all the regions
  ispec = 0
  ipoin = 0

! define a huge offset that is larger than the (Xmax - Xmin) horizontal size of the model
! in order to very easily create split nodes for the crack using a trick
! (the trick is that I shift the X coordinate of one side of the crack by a huge value,
! then perform the sorting of the mesh coordinates, and thus the two points are kept because one
! one is shifted from the other, and after sorting I shift the second point back to its original location).
! By doing so I automatically get split nodes for the two sides of the crack
  if (INTRODUCE_SPLIT_NODES_ON_THE_CRACK_CLEAN) then
    big_X_offset = 1000.d0 * (maxval(xnpgeo) - minval(xnpgeo))
  else
    big_X_offset = 0.d0
  endif

  print *
  x1 = xnpgeo(9)
  z1 = znpgeo(9)
  x2 = xnpgeo(14)
  z2 = znpgeo(14)
  x3 = xnpgeo(8)
  z3 = znpgeo(8)
  angle_crack_to_horizontal = acos(( (x1-x3)*(x2-x3) + (z1-z3)*(z2-z3) ) / &
              ( sqrt((x1-x3)**2 + (z1-z3)**2) * sqrt((x2-x3)**2 + (z2-z3)**2) )) * 180.d0 / PI
  print *,'The crack makes an angle of ',angle_crack_to_horizontal,' degrees with respect to the horizontal'
  print *,'The crack makes an angle of ',90.d0 - angle_crack_to_horizontal,' degrees with respect to the vertical'
  print *

  x1 = xnpgeo(8)
  z1 = znpgeo(8)
  x2 = (xnpgeo(17) + xnpgeo(18)) / 2.d0
  z2 = (znpgeo(17) + znpgeo(18)) / 2.d0
  print *,'The crack has an approximate total length of ',sqrt((x1-x2)**2 + (z1-z2)**2),' m'
  print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! this is the bottom layer (labeled Layer 4 in the InkScape figure)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  material_number = 5

! region 01
  region_number = 1
  ip1 = 1
  ip2 = 2
  ip3 = 5
  ip4 = 4
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 02
  region_number = 2
  ip1 = 2
  ip2 = 3
  ip3 = 6
  ip4 = 5
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! this is Layer 3 in the InkScape figure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  material_number = 4

! region 03
  region_number = 3
  ip1 = 4
  ip2 = 5
  ip3 = 8
  ip4 = 7
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 04
  region_number = 4
  ip1 = 5
  ip2 = 6
  ip3 = 9
  ip4 = 8
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 05
  region_number = 5
  ip1 = 7
  ip2 = 8
  ip3 = 11
  ip4 = 10
  nx = nx_base
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 06
  region_number = 6
  ip1 = 8
  ip2 = 9
  ip3 = 12
  ip4 = 11
  nx = nx_base
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! this is Layer 2 in the InkScape figure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  material_number = 3

! region 07
  region_number = 7
  ip1 = 10
  ip2 = 11
  ip3 = 14
  ip4 = 13
  nx = nx_base
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 08
  region_number = 8
  ip1 = 11
  ip2 = 12
  ip3 = 15
  ip4 = 14
  nx = nx_base
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 09
  region_number = 9
  ip1 = 13
  ip2 = 14
  ip3 = 17
  ip4 = 16
  nx = nx_base
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 10
  region_number = 10
  ip1 = 14
  ip2 = 18
  ip3 = 21
  ip4 = 17
  nx = 6
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 11
  region_number = 11
  ip1 = 14
  ip2 = 15
  ip3 = 19
  ip4 = 18
  nx = nx_base
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 12
  region_number = 12
  ip1 = 16
  ip2 = 17
  ip3 = 21
  ip4 = 20
  nx = nx_base
  nz = 6
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 13
  region_number = 13
  ip1 = 18
  ip2 = 19
  ip3 = 22
  ip4 = 21
  nx = nx_base
  nz = 7
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 14
  region_number = 14
  ip1 = 20
  ip2 = 21
  ip3 = 24
  ip4 = 23
  nx = nx_base
  nz = 3
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 15
  region_number = 15
  ip1 = 21
  ip2 = 22
  ip3 = 25
  ip4 = 24
  nx = nx_base
  nz = 3
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! this is Layer 1 in the InkScape figure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  material_number = 2

! region 16
  region_number = 16
  ip1 = 23
  ip2 = 24
  ip3 = 27
  ip4 = 26
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 17
  region_number = 17
  ip1 = 24
  ip2 = 25
  ip3 = 28
  ip4 = 27
  nx = nx_base
  nz = 19
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! this is the top layer (the water layer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  material_number = 1

! region 18
  region_number = 18
  ip1 = 26
  ip2 = 27
  ip3 = 30
  ip4 = 29
  nx = nx_base
  nz = nz_water
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

! region 19
  region_number = 19
  ip1 = 27
  ip2 = 28
  ip3 = 31
  ip4 = 30
  nx = nx_base
  nz = nz_water
  call generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                                        region_number,big_X_offset)

  print *,'total number of mesh elements generated = ',ispec
  if (ispec /= nspec) stop 'error in total number of mesh elements created'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----  create global mesh numbering, removing the multiples

  allocate(locval(NGLOB_MAX))
  allocate(ind(NGLOB_MAX))
  allocate(ninseg(NGLOB_MAX))
  allocate(iglob(NGLOB_MAX))
  allocate(ifseg(NGLOB_MAX))
  allocate(xp(NGLOB_MAX))
  allocate(yp(NGLOB_MAX))
  allocate(work(NGLOB_MAX))
  allocate(iwork(NGLOB_MAX))

! copy coordinates of the grid points (since sorting will destroy the array order, and we will need to access the original)
  xp(:) = x(:)
  yp(:) = z(:)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Establish initial pointers
  do ispec = 1,nspec
    ieoff = ngnod*(ispec -1)
    do ia = 1,ngnod
      locval (ia+ieoff) = ia+ieoff
    enddo
  enddo

! set up a local geometric tolerance

  xtypdist = +HUGEVAL

  do ispec = 1,nspec

    xminval = +HUGEVAL
    yminval = +HUGEVAL
    xmaxval = -HUGEVAL
    ymaxval = -HUGEVAL
    ieoff = ngnod*(ispec-1)
    do ilocnum = 1,ngnod
      xmaxval = max(xp(ieoff+ilocnum),xmaxval)
      xminval = min(xp(ieoff+ilocnum),xminval)
      ymaxval = max(yp(ieoff+ilocnum),ymaxval)
      yminval = min(yp(ieoff+ilocnum),yminval)
    enddo

! compute the minimum typical "size" of an element in the mesh
    xtypdist = min(xtypdist,xmaxval-xminval)
    xtypdist = min(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol = SMALLVALTOL * xtypdist

! print *,'DK DK xtol = ',xtol

  ifseg(:) = .false.
  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = NGLOB_MAX

  do j = 1,NDIM
!  Sort within each segment
    ioff = 1
    do iseg = 1,nseg
      if (j == 1) then
        call rank (xp(ioff),ind,ninseg(iseg))
      else
        call rank (yp(ioff),ind,ninseg(iseg))
      endif
      call swap(xp(ioff),work,ind,ninseg(iseg))
      call swap(yp(ioff),work,ind,ninseg(iseg))
      call iswap(locval(ioff),iwork,ind,ninseg(iseg))
      ioff = ioff + ninseg(iseg)
    enddo
!  Check for jumps in current coordinate
    if (j == 1) then
      do i = 2,NGLOB_MAX
        if (abs(xp(i)-xp(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else
      do i = 2,NGLOB_MAX
        if (abs(yp(i)-yp(i-1)) > xtol) ifseg(i) = .true.
      enddo
    endif
!  Count up number of different segments
    nseg = 0
    do i = 1,NGLOB_MAX
      if (ifseg(i)) then
        nseg = nseg + 1
        ninseg(nseg) = 1
      else
        ninseg(nseg) = ninseg(nseg) + 1
      endif
    enddo
  enddo
!
!  Assign global node numbers (now sorted lexicographically!)
!
  ig = 0
  do i = 1,NGLOB_MAX
   if (ifseg(i)) then
     ig = ig+1
!! DK DK sep 2018: added this to save the list of grid points without multiples
     x(ig) = xp(i)
     z(ig) = yp(i)
!! DK DK sep 2018: added this to save the list of grid points without multiples

!! DK DK restore the correct physical location of the split nodes, i.e. remove the shift in X that I have introduced
!! DK DK purposely in order to create split nodes automatically in the mesh point sorting routine and point multiple removal
     if (INTRODUCE_SPLIT_NODES_ON_THE_CRACK_CLEAN) then
       if (x(ig) < -SMALLVALTOL) x(ig) = x(ig) + big_X_offset
       if (x(ig) > big_X_offset) x(ig) = x(ig) - big_X_offset
     endif

   endif
   iglob(locval(i)) = ig
  enddo

  nglob = ig

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! get result in my format
  do ispec = 1,nspec
    ieoff = ngnod*(ispec - 1)
    ilocnum = 0
    do ia = 1,ngnod
        ilocnum = ilocnum + 1
        ibool(ia,ispec) = iglob(ilocnum + ieoff)
    enddo
  enddo

  deallocate(locval)
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iglob)
  deallocate(ifseg)
  deallocate(xp)
  deallocate(yp)
  deallocate(work)
  deallocate(iwork)

! check the numbering obtained
  if (minval(ibool) /= 1 .or. maxval(ibool) /= nglob) stop 'Error while generating global numbering'

  print *,'total number of unique mesh points generated = ',nglob
  if (nglob > NGLOB_MAX) stop 'error in total number of unique mesh points created'

! save the mesh files in the format needed by SPECFEM2D for external meshes
  open(unit=20,file='DATA/mesh_file',status='unknown')
  write(20,*) nspec
  do ispec = 1,nspec
    write(20,*) (ibool(ia,ispec),ia=1,ngnod)
  enddo
  close(20)

  open(unit=20,file='DATA/nodes_coords_file',status='unknown')
  write(20,*) nglob
  do ipoin = 1,nglob
    write(20,*) x(ipoin),z(ipoin)
  enddo
  close(20)

  open(unit=20,file='DATA/materials_file',status='unknown')
  do ispec = 1,nspec
    write(20,*) imat(ispec)
  enddo
  close(20)

! geometrical tolerance to detect the edges of the model (1/1000th of the total model size)
  geomtol = 0.085d0 / 1000.d0

  open(unit=20,file='DATA/free_surface_file',status='unknown')
  write(20,*) '0'  ! create an empty file (file listing zero elements)
! write(20,*) nx_base + nx_base ! elements at the top free surface
! ispec_free_surface_found = 0
! do ispec = 1,nspec
! check if the element is located in the water
! and if the top edge of this element belongs to the acoustic free surface, which is located at Z = 0.085
!   if (imat(ispec) == 1 .and. z(ibool(3,ispec)) > 0.085d0 - geomtol .and. z(ibool(4,ispec)) > 0.085d0 - geomtol) then
!     ispec_free_surface_found = ispec_free_surface_found + 1
! ! 'acoustic_surface' contains 1/ element number, 2/ number of nodes that form the free surface,
! ! 3/ first node on the free surface, 4/ second node on the free surface, if relevant (if 2/ is equal to 2)
!     write(20,*) ispec,' 2 ',ibool(3,ispec),ibool(4,ispec)
!   endif
! enddo
! if (ispec_free_surface_found /= nx_base + nx_base) stop 'error in detecting acoustic free surface'
  close(20)

  open(unit=20,file='DATA/axial_elements_file',status='unknown')
  write(20,*) '0'  ! create an empty file (file listing zero elements)
  close(20)

  open(unit=20,file='DATA/absorbing_surface_file',status='unknown')
! write(20,*) '0'  ! create an empty file (file listing zero elements)

! count the number of absorbing edges
  ispec_absorbing_surface_found = 0
  do ispec = 1,nspec
    if (x(ibool(1,ispec)) < geomtol .and. x(ibool(4,ispec)) < geomtol) &
                      ispec_absorbing_surface_found = ispec_absorbing_surface_found + 1
    if (x(ibool(2,ispec)) > 0.050d0 - geomtol .and. x(ibool(3,ispec)) > 0.050d0 - geomtol) &
                      ispec_absorbing_surface_found = ispec_absorbing_surface_found + 1

    if (z(ibool(1,ispec)) < geomtol .and. z(ibool(2,ispec)) < geomtol) &
                      ispec_absorbing_surface_found = ispec_absorbing_surface_found + 1
    if (z(ibool(3,ispec)) > 0.085d0 - geomtol .and. z(ibool(4,ispec)) > 0.085d0 - geomtol) &
                      ispec_absorbing_surface_found = ispec_absorbing_surface_found + 1
  enddo
  write(20,*) ispec_absorbing_surface_found

! write the absorbing edges
! ! 'absorbing_surface' contains 1/ element number, 2/ number of nodes that form the free surface,
! ! 3/ first node on the free surface, 4/ second node on the free surface, if relevant (if 2/ is equal to 2)
  do ispec = 1,nspec
    if (x(ibool(1,ispec)) < geomtol .and. x(ibool(4,ispec)) < geomtol) then
      write(20,*) ispec,' 2 ',ibool(1,ispec),ibool(4,ispec),ILEFT
    endif

    if (x(ibool(2,ispec)) > 0.050d0 - geomtol .and. x(ibool(3,ispec)) > 0.050d0 - geomtol) then
      write(20,*) ispec,' 2 ',ibool(2,ispec),ibool(3,ispec),IRIGHT
    endif

    if (z(ibool(1,ispec)) < geomtol .and. z(ibool(2,ispec)) < geomtol) then
      write(20,*) ispec,' 2 ',ibool(1,ispec),ibool(2,ispec),IBOTTOM
    endif

    if (z(ibool(3,ispec)) > 0.085d0 - geomtol .and. z(ibool(4,ispec)) > 0.085d0 - geomtol) then
      write(20,*) ispec,' 2 ',ibool(3,ispec),ibool(4,ispec),ITOP
    endif
  enddo

  close(20)

  open(unit=20,file='DATA/Surf_acforcing_Bottom_enforcing_mesh',status='unknown')
  write(20,*) '0'  ! create an empty file (file listing zero elements)
  close(20)

  open(unit=20,file='DATA/absorbing_cpml_file',status='unknown')
  write(20,*) '0'  ! create an empty file (file listing zero elements)
  close(20)

  open(unit=20,file='DATA/courbe_eros_nodes',status='unknown')
  write(20,*) '0'  ! create an empty file (file listing zero elements)
  close(20)


! save the mesh in Gnuplot format for verification
  call save_gnuplot_file(x,z,nglob,ibool,ngnod,nspec)

  print *,'done creating and saving the 2D mesh'

!! DK DK
!! DK DK also save a 3D mesh, by extending the mesh identical to itself in the third direction
!! DK DK

  print *,'creating and saving the 3D mesh...'
  print *,'the 3D mesh contains ',nglob * (ny + 1),' unique points, and ',nspec * ny,' elements'

! create SPECFEM3D_Cartesian mesh file for the points
  open(unit=20,file='DATA_3D/nodes_coords_file',status='unknown',action='write')
  write(20,*) nglob * (ny + 1)
  ipoin3D = 0
  do iy = 1,ny + 1
    y = (iy - 1) * deltaY
    do ipoin = 1,nglob
      ipoin3D = ipoin3D + 1
      write(20,*) ipoin3D,x(ipoin),y,z(ipoin)
    enddo
  enddo
  close(20)

  open(unit=20,file='DATA_3D/mesh_file',status='unknown')
  write(20,*) nspec * ny
  ispec3D = 0
  do iy = 1,ny
    do ispec = 1,nspec
      ispec3D = ispec3D + 1
! the 3D element is created by using the 2D element, and the 2D element of the next plane in the Y direction,
! which is nglob grid points further in the list of points, due to the fact that to create the 3D mesh
! we extend of the 2D mesh (which contains nglob points) identical to itself in the third direction.
! To get the point order right, see file numbering_convention_8_nodes.png
      write(20,*) ispec3D,ibool(1,ispec) + (iy-1) * nglob,ibool(2,ispec) + (iy-1) * nglob, &
                          ibool(2,ispec) + iy * nglob,ibool(1,ispec) + iy * nglob, &
                          ibool(4,ispec) + (iy-1) * nglob,ibool(3,ispec) + (iy-1) * nglob, &
                          ibool(3,ispec) + iy * nglob,ibool(4,ispec) + iy * nglob
    enddo
  enddo
  close(20)

  open(unit=20,file='DATA_3D/materials_file',status='unknown')
  ispec3D = 0
  do iy = 1,ny
    do ispec = 1,nspec
      ispec3D = ispec3D + 1
      write(20,*) ispec3D,imat(ispec)
    enddo
  enddo
  close(20)

  print *,'done creating and saving the 3D mesh'

  end program mesh_a_Y_shaped_crack

!---------

  subroutine generate_region(ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,xnpgeo,znpgeo,material_number,x,z,imat,nspec,NGLOB_MAX, &
                    region_number,big_X_offset)

  implicit none

  integer :: ip1,ip2,ip3,ip4,nx,nz,ispec,ipoin,material_number,region_number
  double precision :: big_X_offset

! points defining the different regions to mesh
  double precision, dimension(31) :: xnpgeo,znpgeo

  double precision :: ratio_ix,ratio_iz,ratio_ixplus1,ratio_izplus1
  double precision :: x1,z1,x2,z2,x3,z3,x4,z4
  double precision :: xbackup1,zbackup1
  double precision :: xbackup2
  double precision :: big_offset_local

  integer :: ix,iz

  integer :: nspec,NGLOB_MAX
  double precision, dimension(NGLOB_MAX) :: x,z
  integer, dimension(nspec) :: imat

! introduce split nodes on the crack or not (if not, the crack becomes invisible i.e. it does not exist any more
! and the layers that contain it become homogeneous and continuous; this is useful to provide a reference solution to compare to)
  logical, parameter :: INTRODUCE_SPLIT_NODES_ON_THE_CRACK_UGLY = .true.

  if (INTRODUCE_SPLIT_NODES_ON_THE_CRACK_UGLY) then
!   big_offset_local = 0.03d0 * (50.d0/1000.d0/64.d0) !! 3% of average size of an element
    big_offset_local = 0.0003d0 * (50.d0/1000.d0/64.d0) !! 0.03% of average size of an element
  else
    big_offset_local = 0.d0
  endif

!! DK DK YYYYYYYYYYYYYYYYYYYYYY
  if (INTRODUCE_SPLIT_NODES_ON_THE_CRACK_UGLY) then
    if (region_number == 5) then
      xbackup1 = xnpgeo(ip3)
      xnpgeo(ip3) = xnpgeo(ip3) - big_offset_local
    endif

  if (region_number == 7) then
    xbackup1 = xnpgeo(ip2)
    xnpgeo(ip2) = xnpgeo(ip2) - big_offset_local
    xbackup2 = xnpgeo(ip3)
    xnpgeo(ip3) = xnpgeo(ip3) - big_offset_local
  endif

    if (region_number == 9) then
      xbackup1 = xnpgeo(ip2)
      xnpgeo(ip2) = xnpgeo(ip2) - big_offset_local
    endif

    if (region_number == 10) then
      xbackup1 = xnpgeo(ip1)
      zbackup1 = znpgeo(ip1)
      xnpgeo(ip1) = xnpgeo(ip1) + big_offset_local
      znpgeo(ip1) = znpgeo(ip1) + big_offset_local
    endif
  endif
!! DK DK YYYYYYYYYYYYYYYYYYYYYY

  do iz = 0,nz-1
    do ix = 0,nx-1

! generate one more mesh element
      ispec = ispec + 1

! store the material number
      imat(ispec) = material_number

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

! shift the points to introduce split nodes automatically if needed
      if (region_number == 5 .or. region_number == 7 .or. region_number == 9) then
        ! right edge of the region left of the crack, shift the two points that are on the IRIGHT edge of the element,
        ! which corresponds to the crack
        if (ix == nx-1) then
          x2 = x2 + big_X_offset
          x3 = x3 + big_X_offset
        endif
      endif

      if (region_number == 11) then
        ! left edge of the region right of the crack, shift the two points that are on the ILEFT edge of the element,
        ! which corresponds to the crack
        if (ix == 0) then
          ! use a negative sign this time, because the center point of the Y of the crack needs to be split twice,
          ! once to the right and the second time to the left, since that particular point is shared between three branches
          ! of the crack instead of shared between only two side
          x1 = x1 - big_X_offset
          x4 = x4 - big_X_offset
        endif
      endif

! save the points created
      x(ipoin+1) = x1
      z(ipoin+1) = z1

      x(ipoin+2) = x2
      z(ipoin+2) = z2

      x(ipoin+3) = x3
      z(ipoin+3) = z3

      x(ipoin+4) = x4
      z(ipoin+4) = z4

      ipoin = ipoin + 4

    enddo
  enddo

!! DK DK YYYYYYYYYYYYYYYYYYYYYY
! restore the right geometrical position of the nodes, after having introduced the splitting above
  if (INTRODUCE_SPLIT_NODES_ON_THE_CRACK_UGLY) then
    if (region_number == 5) then
      xnpgeo(ip3) = xbackup1
    endif

    if (region_number == 7) then
      xnpgeo(ip2) = xbackup1
      xnpgeo(ip3) = xbackup2
    endif

    if (region_number == 9) then
      xnpgeo(ip2) = xbackup1
    endif

    if (region_number == 10) then
      xnpgeo(ip1) = xbackup1
      znpgeo(ip1) = zbackup1
    endif
  endif
!! DK DK YYYYYYYYYYYYYYYYYYYYYY

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

!!! this is adapted from src/meshfem2D/save_gnuplot_file.f90 of SPECFEM2D to display the mesh in Gnuplot format

  subroutine save_gnuplot_file(x,z,nglob,ibool,ngnod,nspec)

! creates a Gnuplot file that displays the grid

  implicit none

  integer :: ngnod,nspec,nglob
  integer, dimension(ngnod,nspec) :: ibool
  double precision, dimension(nglob) :: x,z

  ! local parameters
  integer :: ispec,ier

  open(unit=20,file='gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening gnuplot file for writing: gridfile.gnu'

  do ispec = 1,nspec

    ! draw horizontal lines of the grid
    write(20,*) sngl(x(ibool(1,ispec))),sngl(z(ibool(1,ispec)))
    write(20,*) sngl(x(ibool(2,ispec))),sngl(z(ibool(2,ispec)))
    write(20,10)

    write(20,*) sngl(x(ibool(3,ispec))),sngl(z(ibool(3,ispec)))
    write(20,*) sngl(x(ibool(4,ispec))),sngl(z(ibool(4,ispec)))
    write(20,10)

    ! draw vertical lines of the grid
    write(20,*) sngl(x(ibool(1,ispec))),sngl(z(ibool(1,ispec)))
    write(20,*) sngl(x(ibool(4,ispec))),sngl(z(ibool(4,ispec)))
    write(20,10)

    write(20,*) sngl(x(ibool(2,ispec))),sngl(z(ibool(2,ispec)))
    write(20,*) sngl(x(ibool(3,ispec))),sngl(z(ibool(3,ispec)))
    write(20,10)

  enddo

  close(20)

 10 format('')

  ! create a Gnuplot script to display the grid
  open(unit=20,file='plot_gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error saving plot_gridfile.gnu file'

  write(20,*) '#set term wxt'
  write(20,*) 'set term x11'
  write(20,*) 'set term qt'
  write(20,*)
  write(20,*) '#set term postscript portrait color solid "Helvetica" 22'
  write(20,100) '#set output "maillage_du_crack_newer_upwards_OK_avec_petit_mailleur_ecrit_par_Dimitri.ps"'
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

 100 format(a200)

  end subroutine save_gnuplot_file

!
!-----------------------------------------------------------------------
!

  subroutine rank(A,IND,N)
!
! Use Heap Sort (p 233 Numerical Recipes)
!
  implicit none

  integer N
  double precision A(N)
  integer IND(N)

  integer i,j,l,ir,indx
  double precision q

  do j = 1,N
   IND(j)=j
  enddo

  if (n == 1) return
  L=n/2+1
  ir=n
  100 continue
   if (l > 1) then
     l=l-1
     indx=ind(l)
     q=a(indx)
   ELSE
     indx=ind(ir)
     q=a(indx)
     ind(ir)=ind(1)
     ir=ir-1
     if (ir == 1) then
       ind(1)=indx
       return
     endif
   endif
   i=l
   j=l+l
  200 continue
   if (J <= IR) then
      if (J < IR) then
         if (A(IND(j)) < A(IND(j+1))) j=j+1
      endif
      if (q < A(IND(j))) then
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   goto 200
   endif
   IND(I)=INDX
  goto 100

  end subroutine rank

!-----------------------------------------------------------------------

  subroutine swap(a,w,ind,n)
!
! Use IND to sort array A (p 233 Numerical Recipes)
!
  implicit none

  integer n
  double precision A(N),W(N)
  integer IND(N)

  integer j

  W(:) = A(:)

  do j = 1,N
    A(j) = W(ind(j))
  enddo

  end subroutine swap

!-----------------------------------------------------------------------

  subroutine iswap(a,w,ind,n)
!
! Use IND to sort array A
!
  implicit none

  integer n
  integer A(N),W(N),IND(N)

  integer j

  W(:) = A(:)

  do j = 1,N
    A(j) = W(ind(j))
  enddo

  end subroutine iswap

