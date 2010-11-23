
  program decimate_mesh

! mesh decimation, by Nicolas Le Goff, 2008

! cuts each hexahedron of a given mesh in 8 hexahedra (2 * 2 * 2)
! recursively in order to increase mesh density by a factor of 2

  implicit none

  include './constants.h'

  integer :: nelmnts_ext_mesh
  integer, dimension(:,:), allocatable  :: elmnts_ext_mesh
  integer, dimension(:,:), allocatable  :: elmnts_ext_mesh_sub
  integer, dimension(:), allocatable  :: mat_ext_mesh
  integer, dimension(:), allocatable  :: mat_ext_mesh_sub

  integer :: nnodes_ext_mesh
  real, dimension(:,:), allocatable  :: nodes_coords_ext_mesh
  real, dimension(:,:), allocatable  :: nodes_coords_ext_mesh_sub

  real, dimension(NDIM,NSUB+1,NSUB+1,NSUB+1)  :: temporary_nodes
  integer, dimension(NSUB+1,NSUB+1,NSUB+1)  :: temporary_nodes_lookup

  integer, dimension(:), allocatable  :: xadj
  integer, dimension(:), allocatable  :: adjncy
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts

  integer  :: ispec, inode, ispec_neighbours, ispec_neighbours_sub
  integer  :: nnodes_ext_mesh_sub
  integer  :: i, j, k
  integer  :: ix, iy, iz
  integer :: idim

  real :: xtol

  real :: xminval,yminval,xmaxval,ymaxval,xtypdist,zminval,zmaxval

  open(unit=98, file='./mesh', status='old', form='formatted')
  read(98,*) nelmnts_ext_mesh
  allocate(elmnts_ext_mesh(ESIZE,nelmnts_ext_mesh))
   do ispec = 1, nelmnts_ext_mesh
     read(98,*) elmnts_ext_mesh(1,ispec), elmnts_ext_mesh(2,ispec), elmnts_ext_mesh(3,ispec), elmnts_ext_mesh(4,ispec), &
          elmnts_ext_mesh(5,ispec), elmnts_ext_mesh(6,ispec), elmnts_ext_mesh(7,ispec), elmnts_ext_mesh(8,ispec)
  end do
  close(98)

  open(unit=98, file='./mat', status='old', form='formatted')
  allocate(mat_ext_mesh(nelmnts_ext_mesh))
   do ispec = 1, nelmnts_ext_mesh
     read(98,*) mat_ext_mesh(ispec)
  end do
  close(98)

  open(unit=98, file='./nodes_coords', status='old', form='formatted')
  read(98,*) nnodes_ext_mesh
  allocate(nodes_coords_ext_mesh(NDIM,nnodes_ext_mesh))
  do inode = 1, nnodes_ext_mesh
     read(98,*) nodes_coords_ext_mesh(1,inode), nodes_coords_ext_mesh(2,inode),nodes_coords_ext_mesh(3,inode)
  end do
  close(98)


! check that there really are 8 nodes per element.
  do ispec = 1, nelmnts_ext_mesh
    do inode = 1, ESIZE
      do ix = inode+1, ESIZE
         if (elmnts_ext_mesh(inode,ispec) == elmnts_ext_mesh(ix,ispec)) then
            stop 'ERRORERROR'
         endif
      enddo

   enddo
enddo

! set up local geometric tolerances
  xtypdist=+HUGEVAL

  do ispec=1,nelmnts_ext_mesh

  xminval=+HUGEVAL
  yminval=+HUGEVAL
  zminval=+HUGEVAL
  xmaxval=-HUGEVAL
  ymaxval=-HUGEVAL
  zmaxval=-HUGEVAL

  do inode = 1, 8
     xmaxval=max(nodes_coords_ext_mesh(1,elmnts_ext_mesh(inode,ispec)),xmaxval)
     xminval=min(nodes_coords_ext_mesh(1,elmnts_ext_mesh(inode,ispec)),xminval)
     ymaxval=max(nodes_coords_ext_mesh(2,elmnts_ext_mesh(inode,ispec)),ymaxval)
     yminval=min(nodes_coords_ext_mesh(2,elmnts_ext_mesh(inode,ispec)),yminval)
     zmaxval=max(nodes_coords_ext_mesh(3,elmnts_ext_mesh(inode,ispec)),zmaxval)
     zminval=min(nodes_coords_ext_mesh(3,elmnts_ext_mesh(inode,ispec)),zminval)
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = min(xtypdist,xmaxval-xminval)
  xtypdist = min(xtypdist,ymaxval-yminval)
  xtypdist = min(xtypdist,zmaxval-zminval)
  !xtypdist = min(xtypdist,sqrt((xmaxval-xminval)**2 + (ymaxval-yminval)**2 + (zmaxval-zminval)**2))

  enddo

! define a tolerance, small with respect to the minimum size
  xtol=smallval_tol*xtypdist*1.d7

  print *, 'xtypdist' , xtypdist
  print *, 'facteur de tolerance XTOL = ', xtol

  print *, 'xmin', minval(nodes_coords_ext_mesh(1,:))
  print *, 'xmax', maxval(nodes_coords_ext_mesh(1,:))
  print *, 'ymin', minval(nodes_coords_ext_mesh(2,:))
  print *, 'ymax', maxval(nodes_coords_ext_mesh(2,:))
  print *, 'zmin', minval(nodes_coords_ext_mesh(3,:))
  print *, 'zmax', maxval(nodes_coords_ext_mesh(3,:))



! we build the graph
    elmnts_ext_mesh(:,:) = elmnts_ext_mesh(:,:) - 1

    allocate(xadj(1:nelmnts_ext_mesh+1))
    allocate(adjncy(1:MAX_NEIGHBOURS*nelmnts_ext_mesh))
    allocate(nnodes_elmnts(1:nnodes_ext_mesh))
    allocate(nodes_elmnts(1:NSIZE*nnodes_ext_mesh))

    call mesh2dual_ncommonnodes(nelmnts_ext_mesh, nnodes_ext_mesh, elmnts_ext_mesh, xadj, adjncy, nnodes_elmnts, nodes_elmnts,1)

        print *, 'ZZZZ'

    elmnts_ext_mesh(:,:) = elmnts_ext_mesh(:,:) + 1
    adjncy(:) = adjncy(:) + 1
    xadj(:) = xadj(:) + 1

    allocate(elmnts_ext_mesh_sub(ESIZE,nelmnts_ext_mesh*NSUB*NSUB*NSUB))
    allocate(nodes_coords_ext_mesh_sub(NDIM,ESIZE*nelmnts_ext_mesh*(NSUB+1)*(NSUB+1)*(NSUB+1)))
    allocate(mat_ext_mesh_sub(nelmnts_ext_mesh*NSUB*NSUB*NSUB))

    nnodes_ext_mesh_sub = 0

    do ispec = 1, nelmnts_ext_mesh

      do ix = 1, NSUB+1

        temporary_nodes(1,ix,1,1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(2,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,1,1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(2,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(3,ix,1,1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(2,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (ix-1)

        temporary_nodes(1,ix,NSUB+1,1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,NSUB+1,1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(3,ix,NSUB+1,1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (ix-1)

        temporary_nodes(1,ix,1,NSUB+1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,1,NSUB+1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(3,ix,1,NSUB+1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (ix-1)

        temporary_nodes(1,ix,NSUB+1,NSUB+1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(8,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(8,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,NSUB+1,NSUB+1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(8,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(8,ispec))) &
             / real(NSUB))  * (ix-1)
        temporary_nodes(3,ix,NSUB+1,NSUB+1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(8,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(8,ispec))) &
             / real(NSUB))  * (ix-1)

      enddo

      do iy = 1, NSUB+1

        temporary_nodes(1,1,iy,1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(4,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(2,1,iy,1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(4,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(3,1,iy,1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(4,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iy-1)

        temporary_nodes(1,NSUB+1,iy,1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(2,NSUB+1,iy,1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(3,NSUB+1,iy,1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(3,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iy-1)

        temporary_nodes(1,1,iy,NSUB+1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(2,1,iy,NSUB+1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(3,1,iy,NSUB+1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(5,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(5,ispec))) &
             / real(NSUB))  * (iy-1)

        temporary_nodes(1,NSUB+1,iy,NSUB+1) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(6,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(6,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(2,NSUB+1,iy,NSUB+1) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(6,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(6,ispec))) &
             / real(NSUB))  * (iy-1)
        temporary_nodes(3,NSUB+1,iy,NSUB+1) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(6,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(6,ispec))) &
             / real(NSUB))  * (iy-1)

      enddo

      do iz = 1, NSUB+1

        temporary_nodes(1,1,1,iz) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(5,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(2,1,1,iz) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(5,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(3,1,1,iz) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(5,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(1,ispec))) &
             / real(NSUB))  * (iz-1)

        temporary_nodes(1,NSUB+1,1,iz) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(2,NSUB+1,1,iz) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(3,NSUB+1,1,iz) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(2,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(6,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(2,ispec))) &
             / real(NSUB))  * (iz-1)

        temporary_nodes(1,1,NSUB+1,iz) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(2,1,NSUB+1,iz) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(3,1,NSUB+1,iz) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(4,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(8,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(4,ispec))) &
             / real(NSUB))  * (iz-1)

        temporary_nodes(1,NSUB+1,NSUB+1,iz) = nodes_coords_ext_mesh(1,elmnts_ext_mesh(3,ispec)) + &
             ( (nodes_coords_ext_mesh(1,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(1,elmnts_ext_mesh(3,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(2,NSUB+1,NSUB+1,iz) = nodes_coords_ext_mesh(2,elmnts_ext_mesh(3,ispec)) + &
             ( (nodes_coords_ext_mesh(2,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(2,elmnts_ext_mesh(3,ispec))) &
             / real(NSUB))  * (iz-1)
        temporary_nodes(3,NSUB+1,NSUB+1,iz) = nodes_coords_ext_mesh(3,elmnts_ext_mesh(3,ispec)) + &
             ( (nodes_coords_ext_mesh(3,elmnts_ext_mesh(7,ispec)) - nodes_coords_ext_mesh(3,elmnts_ext_mesh(3,ispec))) &
             / real(NSUB))  * (iz-1)

      enddo

      ix = 1
      do iy = 2, NSUB
      do iz = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,ix,1,iz) + &
             ((temporary_nodes(idim,ix,NSUB+1,iz)-temporary_nodes(idim,ix,1,iz)) &
             / real(NSUB))  * (iy-1)) &
             + (temporary_nodes(idim,ix,iy,1) + ((temporary_nodes(idim,ix,iy,NSUB+1)-temporary_nodes(idim,ix,iy,1)) &
             / real(NSUB))  * (iz-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      ix = NSUB+1
      do iy = 2, NSUB
      do iz = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,ix,1,iz) + &
             ((temporary_nodes(idim,ix,NSUB+1,iz)-temporary_nodes(idim,ix,1,iz)) &
             / real(NSUB))  * (iy-1)) &
             + (temporary_nodes(idim,ix,iy,1) + ((temporary_nodes(idim,ix,iy,NSUB+1)-temporary_nodes(idim,ix,iy,1)) &
             / real(NSUB))  * (iz-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      iy = 1
      do ix = 2, NSUB
      do iz = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,1,iy,iz) + &
             ((temporary_nodes(idim,NSUB+1,iy,iz)-temporary_nodes(idim,1,iy,iz)) &
             / real(NSUB))  * (ix-1)) &
             + (temporary_nodes(idim,ix,iy,1) + ((temporary_nodes(idim,ix,iy,NSUB+1)-temporary_nodes(idim,ix,iy,1)) &
             / real(NSUB))  * (iz-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      iy = NSUB+1
      do ix = 2, NSUB
      do iz = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,1,iy,iz) + &
             ((temporary_nodes(idim,NSUB+1,iy,iz)-temporary_nodes(idim,1,iy,iz)) &
             / real(NSUB))  * (ix-1)) &
             + (temporary_nodes(idim,ix,iy,1) + ((temporary_nodes(idim,ix,iy,NSUB+1)-temporary_nodes(idim,ix,iy,1)) &
             / real(NSUB))  * (iz-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      iz = 1
      do ix = 2, NSUB
      do iy = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,1,iy,iz) + &
             ((temporary_nodes(idim,NSUB+1,iy,iz)-temporary_nodes(idim,1,iy,iz)) &
             / real(NSUB))  * (ix-1)) &
             + (temporary_nodes(idim,ix,1,iz) + ((temporary_nodes(idim,ix,NSUB+1,iz)-temporary_nodes(idim,ix,1,iz)) &
             / real(NSUB))  * (iy-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      iz = NSUB+1
      do ix = 2, NSUB
      do iy = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,1,iy,iz) + &
             ((temporary_nodes(idim,NSUB+1,iy,iz)-temporary_nodes(idim,1,iy,iz)) &
             / real(NSUB))  * (ix-1)) &
             + (temporary_nodes(idim,ix,1,iz) + ((temporary_nodes(idim,ix,NSUB+1,iz)-temporary_nodes(idim,ix,1,iz)) &
             / real(NSUB))  * (iy-1))) &
             * 1./2.

      enddo
      enddo
      enddo

      do ix = 2, NSUB
      do iy = 2, NSUB
      do iz = 2, NSUB
      do idim = 1,NDIM
        temporary_nodes(idim,ix,iy,iz) = ((temporary_nodes(idim,1,iy,iz) + &
             ((temporary_nodes(idim,NSUB+1,iy,iz)-temporary_nodes(idim,1,iy,iz)) &
             / real(NSUB))  * (ix-1)) &
             + (temporary_nodes(idim,ix,1,iz) + ((temporary_nodes(idim,ix,NSUB+1,iz)-temporary_nodes(idim,ix,1,iz)) &
             / real(NSUB))  * (iy-1)) &
             + (temporary_nodes(idim,ix,iy,1) + ((temporary_nodes(idim,ix,iy,NSUB+1)-temporary_nodes(idim,ix,iy,1)) &
             / real(NSUB))  * (iz-1))) &
             * 1./3.

      enddo
      enddo
      enddo
      enddo

      temporary_nodes_lookup(:,:,:) = 0

      do ispec_neighbours = xadj(ispec), xadj(ispec+1)-1
        if ( adjncy(ispec_neighbours) < ispec ) then
          do ispec_neighbours_sub = (adjncy(ispec_neighbours)-1)*NSUB*NSUB*NSUB + 1, adjncy(ispec_neighbours)*NSUB*NSUB*NSUB

            do ix = 1, NSUB+1
            do iy = 1, NSUB+1
            do iz = 1, NSUB+1
              do inode = 1, ESIZE
                if ( sqrt( &
                  (temporary_nodes(1,ix,iy,iz)-nodes_coords_ext_mesh_sub(1,elmnts_ext_mesh_sub(inode,ispec_neighbours_sub)))**2 + &
                  (temporary_nodes(2,ix,iy,iz)-nodes_coords_ext_mesh_sub(2,elmnts_ext_mesh_sub(inode,ispec_neighbours_sub)))**2 + &
                  (temporary_nodes(3,ix,iy,iz)-nodes_coords_ext_mesh_sub(3,elmnts_ext_mesh_sub(inode,ispec_neighbours_sub)))**2 ) &
                     < xtol ) then
                  temporary_nodes_lookup(ix,iy,iz) = elmnts_ext_mesh_sub(inode,ispec_neighbours_sub)
                end if

              enddo
            enddo
            enddo
            enddo
          enddo
        end if
      enddo

      do ix = 1, NSUB+1
      do iy = 1, NSUB+1
      do iz = 1, NSUB+1
        if (temporary_nodes_lookup(ix,iy,iz) == 0 ) then
           nnodes_ext_mesh_sub = nnodes_ext_mesh_sub + 1
           temporary_nodes_lookup(ix,iy,iz) = nnodes_ext_mesh_sub
           nodes_coords_ext_mesh_sub(1,nnodes_ext_mesh_sub) = temporary_nodes(1,ix,iy,iz)
           nodes_coords_ext_mesh_sub(2,nnodes_ext_mesh_sub) = temporary_nodes(2,ix,iy,iz)
           nodes_coords_ext_mesh_sub(3,nnodes_ext_mesh_sub) = temporary_nodes(3,ix,iy,iz)
        end if
      enddo
      enddo
      enddo

     do ix = 1, NSUB
     do iy = 1, NSUB
     do iz = 1, NSUB
        elmnts_ext_mesh_sub(1,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix,iy,iz)
        elmnts_ext_mesh_sub(2,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix+1,iy,iz)
        elmnts_ext_mesh_sub(3,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix+1,iy+1,iz)
        elmnts_ext_mesh_sub(4,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix,iy+1,iz)
        elmnts_ext_mesh_sub(5,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix,iy,iz+1)
        elmnts_ext_mesh_sub(6,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix+1,iy,iz+1)
        elmnts_ext_mesh_sub(7,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix+1,iy+1,iz+1)
        elmnts_ext_mesh_sub(8,(ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = temporary_nodes_lookup(ix,iy+1,iz+1)

        mat_ext_mesh_sub((ispec-1)*NSUB*NSUB*NSUB+(ix-1)*NSUB*NSUB+(iy-1)*NSUB+iz) = mat_ext_mesh(ispec)

     enddo
     enddo
     enddo

    enddo

! check that there really are 8 nodes per element.
  do ispec = 1, nelmnts_ext_mesh*NSUB*NSUB*NSUB
    do inode = 1, ESIZE
      do ix = inode+1, ESIZE
         if (elmnts_ext_mesh_sub(inode,ispec) == elmnts_ext_mesh_sub(ix,ispec)) then
            stop 'ERRORERROR'
         endif
      enddo

   enddo
enddo


  print *, 'xmin', minval(nodes_coords_ext_mesh_sub(1,:))
  print *, 'xmax', maxval(nodes_coords_ext_mesh_sub(1,:))
  print *, 'ymin', minval(nodes_coords_ext_mesh_sub(2,:))
  print *, 'ymax', maxval(nodes_coords_ext_mesh_sub(2,:))
  print *, 'zmin', minval(nodes_coords_ext_mesh_sub(3,:))
  print *, 'zmax', maxval(nodes_coords_ext_mesh_sub(3,:))


  open(unit=99, file='./mesh_sub', status='unknown', form='formatted')
  write(99,*) nelmnts_ext_mesh*NSUB*NSUB*NSUB
  do ispec = 1, nelmnts_ext_mesh*NSUB*NSUB*NSUB
     write(99,*) elmnts_ext_mesh_sub(1,ispec), elmnts_ext_mesh_sub(2,ispec), elmnts_ext_mesh_sub(3,ispec), &
          elmnts_ext_mesh_sub(4,ispec), elmnts_ext_mesh_sub(5,ispec), elmnts_ext_mesh_sub(6,ispec), &
          elmnts_ext_mesh_sub(7,ispec), elmnts_ext_mesh_sub(8,ispec)
  end do
  close(99)

  open(unit=99, file='./mat_sub', status='unknown', form='formatted')
  do ispec = 1, nelmnts_ext_mesh*NSUB*NSUB*NSUB
     write(99,*) mat_ext_mesh_sub(ispec)
  end do
  close(99)


  open(unit=99, file='./nodes_coords_sub', status='unknown', form='formatted')
  write(99,*) nnodes_ext_mesh_sub
  do inode = 1, nnodes_ext_mesh_sub
     write(99,*) nodes_coords_ext_mesh_sub(1,inode), nodes_coords_ext_mesh_sub(2,inode), nodes_coords_ext_mesh_sub(3,inode)
  end do
  close(99)

  end program decimate_mesh


  !-----------------------------------------------
  ! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
  !-----------------------------------------------
  subroutine mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts, ncommonnodes)

  include './constants.h'

    integer, intent(in)  :: nelmnts
    integer, intent(in)  :: nnodes
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts
    integer, dimension(0:nelmnts)  :: xadj
    integer, dimension(0:MAX_NEIGHBOURS*nelmnts-1)  :: adjncy
    integer, dimension(0:nnodes-1)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1)  :: nodes_elmnts
    integer, intent(in)  :: ncommonnodes

    integer  :: i, j, k, l, m, nb_edges
    logical  ::  is_neighbour
    integer  :: num_node, n
    integer  :: elem_base, elem_target
    integer  :: connectivity

        print *, 'RRRRRRRRRR'

    !allocate(xadj(0:nelmnts))
    xadj(:) = 0
    !allocate(adjncy(0:MAX_NEIGHBOURS*nelmnts-1))
    adjncy(:) = 0
    !allocate(nnodes_elmnts(0:nnodes-1))
    nnodes_elmnts(:) = 0
    !allocate(nodes_elmnts(0:nsize*nnodes-1))
    nodes_elmnts(:) = 0

    nb_edges = 0


    ! list of elements per node
    do i = 0, esize*nelmnts-1

       nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
       nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1

    end do

    print *, 'nnodes_elmnts'

    ! checking which elements are neighbours ('ncommonnodes' criteria)
    do j = 0, nnodes-1
       do k = 0, nnodes_elmnts(j)-1
          do l = k+1, nnodes_elmnts(j)-1

             connectivity = 0
             elem_base = nodes_elmnts(k+j*nsize)
             elem_target = nodes_elmnts(l+j*nsize)
             do n = 1, esize
                num_node = elmnts(esize*elem_base+n-1)
                do m = 0, nnodes_elmnts(num_node)-1
                   if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                      connectivity = connectivity + 1
                   end if
                end do
             end do

             if ( connectivity >=  ncommonnodes) then

                is_neighbour = .false.

                do m = 0, xadj(nodes_elmnts(k+j*nsize))
                   if ( .not.is_neighbour ) then
                      if ( adjncy(nodes_elmnts(k+j*nsize)*MAX_NEIGHBOURS+m) == nodes_elmnts(l+j*nsize) ) then
                         is_neighbour = .true.


                      end if
                   end if
                end do
                if ( .not.is_neighbour ) then
                   adjncy(nodes_elmnts(k+j*nsize)*MAX_NEIGHBOURS+xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)
                   xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                   adjncy(nodes_elmnts(l+j*nsize)*MAX_NEIGHBOURS+xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)
                   xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                end if
             end if
          end do
       end do
    end do

    ! making adjacency arrays compact (to be used for partitioning)
    do i = 0, nelmnts-1
       k = xadj(i)
       xadj(i) = nb_edges
       do j = 0, k-1
          adjncy(nb_edges) = adjncy(i*MAX_NEIGHBOURS+j)
          nb_edges = nb_edges + 1
       end do
    end do

    xadj(nelmnts) = nb_edges


  end subroutine mesh2dual_ncommonnodes

