
  program decompose_mesh_SCOTCH

  use part_decompose_mesh_SCOTCH
  implicit none

  include './constants_decompose_mesh_SCOTCH.h'
  include "./scotchf.h"

  integer, parameter :: nparts = 8

  integer(long) :: nspec
  integer, dimension(:,:), allocatable  :: elmnts
  integer, dimension(:), allocatable  :: mat
  integer, dimension(:), allocatable  :: part

  integer :: nnodes
  double precision, dimension(:,:), allocatable  :: nodes_coords

  integer, dimension(:), allocatable  :: xadj
  integer, dimension(:), allocatable  :: adjncy
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts

  integer, dimension(:), allocatable  :: vwgt
  integer, dimension(:), allocatable  :: adjwgt
  integer, dimension(5)  :: metis_options

  integer, dimension(:), pointer  :: glob2loc_elmnts
  integer, dimension(:), pointer  :: glob2loc_nodes_nparts
  integer, dimension(:), pointer  :: glob2loc_nodes_parts
  integer, dimension(:), pointer  :: glob2loc_nodes

  integer, dimension(:), pointer  :: tab_size_interfaces, tab_interfaces
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces
  integer  ::  ninterfaces
  integer  :: my_ninterface

  integer  :: nb_materials
  double precision, dimension(:), allocatable :: cs
  integer, dimension(:), allocatable :: num_material

  integer(long)  :: nsize  ! Max number of elements that contain the same node.
  integer  :: edgecut
  integer  :: nb_edges

  integer  :: ispec, inode
  integer  :: wgtflag
  integer  :: num_start
  integer  :: ngnod
  integer  :: max_neighbour   ! Real maximum number of neighbours per element
  integer(long)  :: sup_neighbour   ! Majoration of the maximum number of neighbours per element

  integer  :: ipart, nnodes_loc, nspec_loc
  character(len=256)  :: prname

  logical, dimension(:), allocatable :: mask_nodes_elmnts
  integer, dimension(:), allocatable :: used_nodes_elmnts

!!!! NL NL for SCOTCH partitioner
 double precision, dimension(SCOTCH_GRAPHDIM)  :: scotchgraph
 double precision, dimension(SCOTCH_STRATDIM)  :: scotchstrat
 character(len=256), parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'
 integer  :: ierr
!!!! NL NL


  ngnod = esize

! read the elements, material and nodes files
  open(unit=98, file='../model_asteroid_subdivide/mesh', status='old', form='formatted')
  read(98,*) nspec
  allocate(elmnts(esize,nspec))
   do ispec = 1, nspec
     read(98,*) elmnts(1,ispec), elmnts(2,ispec), elmnts(3,ispec), elmnts(4,ispec), &
          elmnts(5,ispec), elmnts(6,ispec), elmnts(7,ispec), elmnts(8,ispec)
  end do
  close(98)

  open(unit=98, file='../model_asteroid_subdivide/mat', status='old', form='formatted')
  allocate(mat(nspec))
   do ispec = 1, nspec
     read(98,*) mat(ispec)
  end do
  close(98)

  open(unit=98, file='../model_asteroid_subdivide/nodes_coords', status='old', form='formatted')
  read(98,*) nnodes
  allocate(nodes_coords(3,nnodes))
  do inode = 1, nnodes
     read(98,*) nodes_coords(1,inode), nodes_coords(2,inode), nodes_coords(3,inode)
  end do
  close(98)

  allocate(mask_nodes_elmnts(nnodes))
  allocate(used_nodes_elmnts(nnodes))
  mask_nodes_elmnts(:) = .false.
  used_nodes_elmnts(:) = 0
  do ispec = 1, nspec
    do inode = 1, ESIZE
      mask_nodes_elmnts(elmnts(inode,ispec)) = .true.
      used_nodes_elmnts(elmnts(inode,ispec)) = used_nodes_elmnts(elmnts(inode,ispec)) + 1
    enddo
  enddo
  nsize = maxval(used_nodes_elmnts(:))
  sup_neighbour = ngnod * nsize - (ngnod + (ngnod/2 - 1)*nfaces)
  print*, 'nsize = ',nsize, 'sup_neighbour = ', sup_neighbour

  do inode = 1, nnodes
    if (.not. mask_nodes_elmnts(inode)) then
      stop 'ERROR : nodes not used.'
    endif
  enddo
!   if (maxval(used_nodes_elmnts(:))>nsize) then
!     stop 'ERROR : increase nsize or modify the mesh.'
!   endif

  elmnts(:,:) = elmnts(:,:) - 1

  allocate(xadj(1:nspec+1))
  allocate(adjncy(1:sup_neighbour*nspec))
  allocate(nnodes_elmnts(1:nnodes))
  allocate(nodes_elmnts(1:nsize*nnodes))

  call mesh2dual_ncommonnodes(nspec, nnodes, nsize, sup_neighbour, elmnts, xadj, adjncy, nnodes_elmnts, &
       nodes_elmnts, max_neighbour, 1)
  print*, 'max_neighbour = ',max_neighbour


 ! elmnts(:,:) = elmnts(:,:) + 1
 ! adjncy(:) = adjncy(:) + 1
 ! xadj(:) = xadj(:) + 1
!  allocate(vwgt(0:nspec-1))
  nb_edges = xadj(nspec+1)
!   allocate(adjwgt(0:nb_edges-1))
!   vwgt(:) = 1
!   adjwgt(:) = 1

!   metis_options(1) = 0
!   metis_options(2) = 3
!   metis_options(3) = 1
!   metis_options(4) = 1
!   metis_options(5) = 0


!   num_start = 0
!   wgtflag = 0

  allocate(part(1:nspec))


! Old metis partitioning
!   call METIS_PartGraphRecursive(nspec, xadj(1), adjncy(1), vwgt(0), adjwgt(0), wgtflag, num_start, nparts, &
!        metis_options, edgecut, part(1));

! SCOTCH partitioning
    call scotchfstratinit (scotchstrat(1), ierr)
     if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot initialize strat'
    endif

    call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), ierr)
     if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot build strat'
    endif

    call scotchfgraphinit (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot initialize graph'
    endif

    call scotchfgraphbuild (scotchgraph (1), 0, nspec, xadj (1), xadj (1), &
         xadj (1), xadj (1), nb_edges, adjncy (1), adjncy (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot build graph'
    endif

    call scotchfgraphcheck (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Invalid check'
    endif

    call scotchfgraphpart (scotchgraph (1), nparts, scotchstrat(1),part(1),ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot part graph'
    endif

    call scotchfgraphexit (scotchgraph (1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot destroy graph'
    endif

    call scotchfstratexit (scotchstrat(1), ierr)
    if (ierr /= 0) then
       stop 'ERROR : MAIN : Cannot destroy strat'
    endif


! local number of each element for each partition
  call Construct_glob2loc_elmnts(nspec, part, nparts, glob2loc_elmnts)

! local number of each node for each partition
  call Construct_glob2loc_nodes(nspec, nnodes,nsize, nnodes_elmnts, nodes_elmnts, part, nparts, &
       glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes)

  nb_materials = 1
  allocate(cs(nb_materials))
  allocate(num_material(nspec))
  cs(:) = 1000.d0
  num_material(:) = 1

  call Construct_interfaces(nspec, nparts, sup_neighbour, part, elmnts, xadj, adjncy, tab_interfaces, &
             tab_size_interfaces, ninterfaces, nb_materials, cs, num_material)

  allocate(my_interfaces(0:ninterfaces-1))
  allocate(my_nb_interfaces(0:ninterfaces-1))



  do ipart = 0, nparts-1

     !write(prname, "('/Database',i5.5)") ipart
     write(prname, "(i6.6,'_Database')") ipart
     open(unit=15,file='./OUTPUT_FILES/proc'//prname,status='unknown', action='write', form='formatted')

     call write_glob2loc_nodes_database(15, ipart, nnodes_loc, nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, nnodes, 1)
     call write_partition_database(15, ipart, nspec_loc, nspec, elmnts, glob2loc_elmnts, glob2loc_nodes_nparts, &
          glob2loc_nodes_parts, glob2loc_nodes, part, num_material, ngnod, 1)

     write(15,*) nnodes_loc
     call write_glob2loc_nodes_database(15, ipart, nnodes_loc, nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, nnodes, 2)
     write(15,*) nspec_loc
     call write_partition_database(15, ipart, nspec_loc, nspec, elmnts, glob2loc_elmnts, glob2loc_nodes_nparts, &
          glob2loc_nodes_parts, glob2loc_nodes, part, num_material, ngnod, 2)

     call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, nparts, ipart, ninterfaces, &
          my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, 1)
     write(15,*) my_ninterface, maxval(my_nb_interfaces)
     call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, nparts, ipart, ninterfaces, &
          my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, 2)

     close(15)

  enddo

  end program decompose_mesh_SCOTCH

