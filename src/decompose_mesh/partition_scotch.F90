!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! SCOTCH is a project carried out within the Satanas team of the Laboratoire Bordelais de Recherche en Informatique (LaBRI).
! It is part of the ScAlApplix project of INRIA Bordeaux - Sud-Ouest.
!
! Its purpose is to apply graph theory, with a divide and conquer approach, to scientific computing problems
! such as graph and mesh partitioning, static mapping, and sparse matrix ordering, in application domains ranging
! from structural mechanics to operating systems or bio-chemistry.
!
! The SCOTCH distribution is a set of programs and libraries which implement the static mapping and
! sparse matrix reordering algorithms developed within the SCOTCH project.
!
! http://www.labri.fr/perso/pelegrin/scotch/
! https://gitlab.inria.fr/scotch/scotch
!
! reference:
! F. Pellegrini and J. Roman,
! SCOTCH: A Software Package for Static Mapping by Dual Recursive Bipartitioning of Process and Architecture Graphs.
! Proceedings of HPCN'96, Brussels, Belgium. LNCS 1067, pages 493-498. Springer, April 1996.
! https://link.springer.com/chapter/10.1007/3-540-61142-8_588
!
!
! support for SCOTCH:
!  By default, SCOTCH is chosen as partitioning algorithm for the xdecompose_mesh tool
!  and compilation support should have been added automatically by running the configure script.
!
!  You can check the Makefile in the SPECFEM3D/ root directory:
!    - should have added the path to the SCOTCH library files and compiler flag -DUSE_SCOTCH
!
!  To select SCOTCH partioning, in Par_file choose:
!    PARTITIONING_TYPE = 1

module scotch_par

  implicit none

#if defined(USE_SCOTCH)
  include 'scotchf.h'

  ! SCOTCH
  double precision, dimension(SCOTCH_GRAPHDIM)  :: scotchgraph
  double precision, dimension(SCOTCH_STRATDIM)  :: scotchstrat
!!!!!! character(len=*), parameter :: scotch_strategy='b{job=t,map=t,poli=S,sep=h{pass=30}}'

#endif

end module scotch_par

!
!----------------------------------------------------------------------------------------------
!

  subroutine partition_scotch()

#if defined(USE_SCOTCH)

  use constants, only: BALANCE_P_LEVELS

  use decompose_mesh_par, only: NSPEC, nparts, part, nb_edges, elmnts_load, xadj, adjncy, &
    LTS_MODE

  use scotch_par

  implicit none

  ! local parameters
  integer :: ier

  ! checks if partitioning needed
  if (nparts == 1) return

  ! SCOTCH partitioning
  print *,'SCOTCH partitioning'

  ! SCOTCH partitioning

  ! we use default strategy for partitioning, thus omit specifing explicit strategy.

  ! workflow preferred by F. Pellegrini (SCOTCH):
  !!This comes from the fact that, in version 5.1.8, the name
  !!for the "recursive bisection" method has changed from "b"
  !!("bipartitioning") to "r" ("recursive").
  !!
  !!As a general rule, do not try to set up strategies by
  !!yourself. The default strategy in Scotch will most probably
  !!provide better results. To use it, just call:
  !!
  !!SCOTCHFstratInit (),
  !!
  !!and use this "empty" strategy in the mapping routine
  !!(consequently, no call to SCOTCHFstratGraphMap () is
  !!required).
  !!
  !!This will make you independent from further changes
  !!(improvements) in the strategy syntax.
  !!And you should see an improvement in performance, too,
  !!as your hand-made strategy did not make use of the
  !!multi-level framework.
  call scotchfstratinit (scotchstrat(1), ier)
   if (ier /= 0) then
     stop 'Scotch Error : MAIN : Cannot initialize strategy'
  endif

  ! resets SCOTCH random number generator to produce deterministic partitions
  call scotchfrandomReset()

  !call scotchfstratgraphmap (scotchstrat(1), trim(scotch_strategy), ier)
  ! if (ier /= 0) then
  !   stop 'Scotch Error : MAIN : Cannot build strategy'
  !endif

  call scotchfgraphinit (scotchgraph (1), ier)
  if (ier /= 0) then
     stop 'Scotch Error : MAIN : Cannot initialize graph'
  endif

  ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
  ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
  !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
  !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
  !                    #(7) number_of_arcs                    #(8) adjacency_array
  !                    #(9) arc_load_array (optional)      #(10) ierror
  !
  !
  ! note: for graph partitioning, each spectral-element is represented by a single vertex.
  !       xadj(:)  : number of neighbors connected to this single vertex (and listed in adjncy array)
  !       adjncy(:): every vertex has neighbors connected to it and specified in this adjacency array.
  !       load(:)  : a non-zero integer value representing the associated load of the vertex (= element)


  if (LTS_MODE .and. BALANCE_P_LEVELS) then
    ! LTS mode partitioning for each p-level
    call lts_partition_scotch()

  else
    ! over-all scotch partitioning
    call scotchfgraphbuild (scotchgraph (1), 0, nspec, &
                            xadj (1), xadj (1), &
                            elmnts_load (1), xadj (1), &
                            nb_edges, adjncy (1), &
                            adjncy (1), ier)
    if (ier /= 0) then
      stop 'Scotch ERROR : MAIN : Cannot build graph'
    endif

    call scotchfgraphcheck (scotchgraph(1), ier)
    if (ier /= 0) then
      stop 'Scotch ERROR : MAIN : Invalid check'
    endif

    call scotchfgraphpart (scotchgraph(1), nparts, scotchstrat(1),part(1),ier)
    if (ier /= 0) then
      stop 'Scotch ERROR : MAIN : Cannot part graph'
    endif

  endif ! LTS_MODE

  call scotchfgraphexit (scotchgraph(1), ier)
  if (ier /= 0) then
    stop 'Scotch ERROR : MAIN : Cannot destroy graph'
  endif

  call scotchfstratexit (scotchstrat(1), ier)
  if (ier /= 0) then
    stop 'Scotch ERROR : MAIN : Cannot destroy strategy'
  endif

#else

  ! compiled without support
  print *,'Scotch Error: Please re-compile with -DUSE_SCOTCH to use SCOTCH partitioning.'
  stop 'No SCOTCH support'

#endif

  end subroutine partition_scotch

!
!----------------------------------------------------------------------------------------------
!

#if defined(USE_SCOTCH)

  subroutine lts_partition_scotch()

  ! note: SCOTCH is unable to handle constrained partitioning as Patoh.
  !       the goal would be to have about the same number of p-elements in each partition (due to MPI scheme).
  !       as an ugly work-around, we partition for each p-level again, using newly assigned element weights

  use decompose_mesh_par, only: NSPEC, nparts, part, nb_edges, elmnts_load, xadj, adjncy
  ! p-levels
  use constants, only: BALANCE_P_LEVELS_PARTIAL,P_LEVEL_PARTIAL_SUBSET_MINIMUM
  use decompose_mesh_par, only: p_level, num_p_level, ispec_p_refine, num_ispec_level

  use scotch_par

  implicit none

  ! local parameters
  integer :: ier

  ! scotch p-level partitioning
  integer,dimension(:),allocatable :: part_tmp,elmnts_load_tmp
  integer,dimension(:),allocatable :: ispec_global,ispec_local
  integer,dimension(:),allocatable :: xadj_tmp,adjncy_tmp
  integer :: nb_edges_tmp
  integer :: i,j,k,ispec,iproc

  integer :: ilevel,nspec_p
  integer :: nparts_partial
  integer :: p,inum

  ! allocates temporary arrays
  allocate(part_tmp(1:nspec), &
           ispec_global(1:nspec), &
           ispec_local(1:nspec), &
           elmnts_load_tmp(1:nspec), &
           xadj_tmp(1:nspec+1), &
           adjncy_tmp(1:nb_edges+1), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating array part_tmp,...'

  ! partitions each p-level separately
  do ilevel = 1,num_p_level
    p = p_level(ilevel)
    nspec_p = num_ispec_level(ilevel)

    ! element lookup tables
    ispec_global(:) = 0
    ispec_local(:) = 0
    inum = 0
    do ispec = 1,nspec
      if (ispec_p_refine(ispec) == p) then
        inum = inum + 1
        if (inum > nspec_p) stop 'Error index out of bounds of nspec_p'
        ispec_global(inum) = ispec
        ispec_local(ispec) = inum
      endif
    enddo
    if (inum /= nspec_p) stop 'Error number of p-vertices not equals to number of p-elements'

    ! builds up local arrays only for p-elements
    xadj_tmp(:) = 0
    adjncy_tmp(:) = -1
    elmnts_load_tmp(:) = 0
    nb_edges_tmp = 1
    do inum = 1,nspec_p
      ispec = ispec_global(inum)

      ! element load
      elmnts_load_tmp(inum) = elmnts_load(ispec)

      ! number of neighbors
      ! note: xadj() values start from 0; adjncy() values start from 0 to nspec-1
      !       we shift these by +1 since arrays in this subroutine are defined between 1 to ..
      !
      !       xadj(i)  -> adjncy( . ) : values in xadj point to the first index in adjncy() for vertex i
      !                   adjncy() holds all indices of neighbor vertices(=elements)
      ! note: for xadj_tmp and adjncy_tmp we start indexing at 1
      xadj_tmp(inum) = nb_edges_tmp
      do j = xadj(ispec),xadj(ispec+1)-1
        ! gets neighbor index
        k = adjncy(j+1) + 1
        ! checks bounds
        if (k < 1 .or. k > nspec) stop 'Error adjncy index'
        ! adds vertex if it belongs to p-elements
        if (ispec_p_refine(k) == p) then
          ! we store the local index between 1 and nspec_p
          i = ispec_local( k )
          adjncy_tmp(nb_edges_tmp) = i
          nb_edges_tmp = nb_edges_tmp + 1
        endif
      enddo
    enddo
    ! last entry for xadj for a contiguous range of indices ("compact edge array")
    xadj_tmp(nspec_p+1) = nb_edges_tmp
    nb_edges_tmp = nb_edges_tmp - 1

    ! debug
    !print *,'xadj:',xadj(ispec_global(1):ispec_global(1)+1),xadj(nspec-1:nspec+1)
    !print *,'xadj:',xadj(nspec+1)
    !print *,'adjncy:'
    !do j = xadj(ispec_global(1)),xadj(ispec_global(1)+1)
    !  print *,j-xadj(ispec_global(1))+1,j,adjncy(j+1)+1,ispec_local(adjncy(j+1)+1)
    !enddo
    !print *
    !print *,'temporary:',ispec_global(1)
    !print *,'xadj:',xadj_tmp(1:2),xadj_tmp(nspec_p-1:nspec_p+2)
    !print *,'xadj:',xadj_tmp(nspec_p+1)
    !print *,'adjncy:'
    !do j = xadj_tmp(1),xadj_tmp(2)
    !  print *,j-xadj_tmp(1)+1,j,adjncy_tmp(j),ispec_global(adjncy_tmp(j))
    !enddo

    ! checks ranges
    if (minval(adjncy_tmp(1:nb_edges_tmp)) < 1 .or. maxval(adjncy_tmp(1:nb_edges_tmp)) > nspec_p) then
      print *,'Error adjncy bounds invalid', minval(adjncy_tmp(1:nb_edges_tmp)),maxval(adjncy_tmp(1:nb_edges_tmp))
      stop 'Error adjncy bounds invalid'
    endif

    ! sets up p-element partitioning
    part_tmp(:) = -1

    ! user output
    print *,'  p-level partitioning: p-level =',ilevel,'p =',p,'elements =',nspec_p !,'edges =',nb_edges_tmp

    ! scotch partitioning
    ! arguments: #(1) graph_structure       #(2)baseval (either 0/1)    #(3)vertnbr (number_of_vertices)
    !            #(4)verttab (adjacency_index_array)            #(5)vendtab (adjacency_end_index_array (optional))
    !            #(6)velotab (vertex_load_array (optional))     #(7)vlbltab (vertex_label_array)
    !            #(7)edgenbr (number_of_arcs)                   #(8)edgetab (adjacency_array)
    !            #(9)edlotab (arc_load_array (optional))        #(10)ierror
    !
    ! Since, in Fortran, there is no null reference, passing the scotchfgraphbuild routine
    ! a reference equal to verttab in the velotab or vlbltab fields makes them be considered as missing arrays.
    ! The same holds for edlotab when it is passed a reference equal to edgetab.
    ! Setting vendtab to refer to one cell after verttab yields the same result,
    ! as it is the exact semantics of a compact vertex array.

    call scotchfgraphbuild (scotchgraph(1), 1, nspec_p, &
                            xadj_tmp(1), xadj_tmp(1), &
                            elmnts_load_tmp(1), xadj_tmp(1), &
                            nb_edges_tmp, adjncy_tmp(1), &
                            adjncy_tmp(1), ier)
    if (ier /= 0) stop 'ERROR : MAIN : Cannot build graph'

    ! checks scotch graph
    call scotchfgraphcheck (scotchgraph(1), ier)
    if (ier /= 0) stop 'Scotch ERROR : MAIN : Invalid check'

    ! computes partition based on scotch graph
    if (BALANCE_P_LEVELS_PARTIAL) then
      ! partitions level using only a subset of nparts processes
      nparts_partial = int( nspec_p / P_LEVEL_PARTIAL_SUBSET_MINIMUM )
      if (nparts_partial > nparts) nparts_partial = nparts
      if (nparts_partial < 1) nparts_partial = 1

      ! user output
      print *,'    partial partitioning: nparts_partial =',nparts_partial, &
        'using subset maximum',P_LEVEL_PARTIAL_SUBSET_MINIMUM

      ! partitions among subset nparts_partial processes
      call scotchfgraphpart (scotchgraph(1), nparts_partial, scotchstrat(1),part_tmp(1),ier)
      if (ier /= 0) stop 'Scotch ERROR : MAIN : Cannot part graph'

    else

      ! partitions among all nparts processes
      call scotchfgraphpart (scotchgraph(1), nparts, scotchstrat(1),part_tmp(1),ier)
      if (ier /= 0) stop 'Scotch ERROR : MAIN : Cannot part graph'

    endif

    ! frees scotch graph for subsequent calls of scotch again
    call scotchfgraphfree (scotchgraph(1),ier)
    if (ier /= 0) stop 'Scotch ERROR : MAIN : Cannot free graph'

    ! stitch partitioning together
    do inum = 1,nspec_p
      ispec = ispec_global(inum)
      part(ispec) = part_tmp(inum)
    enddo

    ! user output: counts number of p-elements in each partition
    do iproc = 0,nparts-1
      inum = 0
      do ispec = 1,nspec
        if (ispec_p_refine(ispec) == p) then
          if (part(ispec) == iproc) inum = inum + 1
        endif
      enddo
      ! user output
      print *,'    partition ',iproc, '       has ',inum,' p-elements'
    enddo

  enddo ! num_p_level
  print *

  ! reassemble partitions trying to optimize for communication costs
  call remap_partitions(part)

  ! frees memory
  deallocate(part_tmp,ispec_global,ispec_local,elmnts_load_tmp,xadj_tmp,adjncy_tmp)

  end subroutine lts_partition_scotch

#endif

!
!----------------------------------------------------------------------------------------------
!

  subroutine remap_partitions(part)

  use constants, only: SCOTCH_P_REMAP

  use decompose_mesh_par, only: nspec,num_p_level, ispec_p_refine, p_level, nparts

  implicit none
  integer, dimension(nspec) :: part

  integer :: current_partition
  logical, dimension(:), allocatable :: partition_available
  integer, dimension(:), allocatable :: part_remap, part_remap_inverse
  integer :: ipart, unused_count, best_choice, ilevel
  integer :: p, ispec

  ! Try to reconnect the individually partitioned levels in a
  ! semi-optimal way.  Similar to traveling salesman problem, which is
  ! known to be NP-complete. We use a simple greedy algorithm, which
  ! should be sufficient.

  if (SCOTCH_P_REMAP) then
      allocate(partition_available(0:nparts-1))
      allocate(part_remap(0:nparts-1))
      allocate(part_remap_inverse(0:nparts-1))

      do ilevel = 1,num_p_level-1
        part_remap(:) = -1
        p = p_level(ilevel)
        partition_available(:) = .true.

        ! cycle through partitions. last partition is forced to take last remaining choice
        do current_partition = 0,nparts-2
          call find_best_choice(nparts,current_partition,part,ilevel,partition_available,best_choice)
          partition_available(best_choice) = .false.
          part_remap(current_partition) = best_choice
        enddo ! current_partition

        ! get unused choice
        unused_count = 0
        do ipart = 0,nparts-1
          if (partition_available(ipart)) then
            part_remap(nparts-1) = ipart
            partition_available(ipart) = .false.
            ! sanity check
            unused_count = unused_count + 1
          endif
        enddo

        ! sanity check
        if (unused_count /= 1) then
          print *, "Error: unused_count should be 1, was = ", unused_count
          stop 'Error invalid unused_count'
        endif

        ! sanity check
        if (minval(part_remap) == -1) then
          print *, "Error: SCOTCH-P reassignment mapping not complete"
          print *, "part_remap=", part_remap
          stop 'Error invalid part_remap'
        endif

        ! build inverse map
        do ipart = 0,nparts-1
          part_remap_inverse(part_remap(ipart)) = ipart
        enddo

        ! remap parts, 1 level up. easier to do this as we go up the levels
        do ispec=1,nspec
          if (ispec_p_refine(ispec) == p_level(ilevel + 1)) then
            part(ispec) = part_remap_inverse(part(ispec))
          endif
        enddo

      enddo ! ilevel

    endif

  end subroutine remap_partitions

!
!----------------------------------------------------------------------------------------------
!

  subroutine find_best_choice(nparts,current_partition,part,ilevel,partition_available,best_choice)

  use decompose_mesh_par, only: nspec,xadj,adjncy,ispec_p_refine,p_level

  implicit none
  integer :: nparts, current_partition, ilevel, best_choice, ispec, num_neighbors, ipart, neighbor
  integer, dimension(nspec) :: part
  logical, dimension(0:nspec-1) :: partition_available
  integer :: number_neighbors(0:nparts-1)
  integer :: iconn

  ! dummy to test compilation
  ! best_choice = current_partition
  number_neighbors(:) = 0

  do ispec = 1,nspec
    if (ispec_p_refine(ispec) == p_level(ilevel)) then
      if (part(ispec) == current_partition) then
        num_neighbors = xadj(ispec+1)-xadj(ispec)

        do iconn = 1,num_neighbors
          neighbor = adjncy(xadj(ispec) + iconn)+1
          ! only care if neighbor is 1 level higher
          if (ispec_p_refine(neighbor) == p_level(ilevel+1)) then
            ! neighbor's communication cost added to its partition
            number_neighbors(part(neighbor)) = number_neighbors(part(neighbor)) + 1
          endif ! neighbor ilevel+1
        enddo ! iconn

      endif ! level
    endif ! partition
  enddo ! ispec

  ! print *, "number_neighbors=",number_neighbors
  do ipart = 0,nparts-1
    best_choice = maxloc(number_neighbors,1)-1
    ! print *, "best_choice:", best_choice

    if (partition_available(best_choice)) then
      exit
    endif

    ! remove this choice
    ! print *, "choice not available"
    number_neighbors(best_choice) = -1
  enddo

  ! print *, "best_choice=", best_choice, "for ilevel", ilevel, "for current_partition", current_partition
  end subroutine find_best_choice
