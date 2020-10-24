!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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

! PaToH is a sequential, multilevel, hypergraph partitioning tool that can be used to solve various
! combinatorial scientific computing problems that could be modeled as hypergraph partitioning problem,
! including sparse matrix partitioning, ordering, and load balancing for parallel processing.
!
! PaToH can be downloaded here (accessed Jan 2019):
! https://www.cc.gatech.edu/~umit/software.html
! https://www.cc.gatech.edu/~umit/PaToH/manual.pdf

! reference:
! Cevdet Aykanata, B. Barla Cambazoglu, Bora Ucar, 2008:
! Multi-level direct K-way hypergraph partitioning with multiple constraints and fixed vertices,
! Journal of Parallel and Distributed Computing, Volume 68, Issue 5, May 2008, Pages 609-625.
! https://www.sciencedirect.com/science/article/pii/S0743731507001724
!
!
! support for PATOH:
!  If you want to use PATOH instead of SCOTCH, then the program xdecompose_mesh must be compiled with PATOH support:
!
!    - Add the path to the PATOH library files and add compiler flag -DUSE_PATOH in Makefile in the SPECFEM3D/ directory
!    - re-compile xdecompose_mesh: > make xdecompose_mesh
!
!  To select PATOH instead of SCOTCH, in Par_file choose:
!    PARTITIONING_TYPE = 3


  subroutine partition_patoh()

  use part_decompose_mesh, only: mesh2hyp_ncommonnodes
  use decompose_mesh_par

  implicit none

  ! local parameters
  integer, dimension(:), allocatable :: xpins
  integer, dimension(:), allocatable :: pins
  integer :: max_neighbor

  integer :: bcuttype
  integer :: newlvlids
  integer :: useFixCells
  integer :: cut
  integer, dimension(:), allocatable :: belmts_mc_loads
  integer, dimension(:), allocatable :: nwghts
  integer :: ier

  ! checks if partitioning needed
  if (nparts == 1) return

  ! PATOH partitioning
  print *,'PATOH partitioning'

  ! determines maximum neighbors based on 1 common node
  print *, 'patoh: nnodes ', nnodes, ' nelems ', nspec

  allocate(xpins(1:nnodes+1),stat=ier)
  if (ier /= 0) stop 'Error allocating array xpins'
  allocate(pins(1:sup_neighbor*nspec+1),stat=ier)
  if (ier /= 0) stop 'Error allocating array pins'
  xpins(:) = -1
  pins(:) = -1

  print *, 'patoh: calling mesh2hyp_ncommonnodes in partition_patoh'

  call mesh2hyp_ncommonnodes(nspec, nnodes, nsize, sup_neighbor, elmnts, xadj, adjncy, &
                             xpins, pins, &
                             nnodes_elmnts, nodes_elmnts, max_neighbor, ncommonnodes, NGNOD)


  ! initializes
  useFixCells = 0 ! for moving higher p-elements into single partition

  ! cut type, either:
  ! 1 == "Connectivity-1" metric
  ! 2 == for cutnet metric
  bcuttype = 1
  print *,'patoh: cut type: ',bcuttype
  print *

  if (LTS_MODE) then
    call plevel_partitioning()
  else
    ! no-LTS mode
    ! single level constraints
    newlvlids = 1

    ! weighting according to element loads
    allocate(belmts_mc_loads(nspec))
    ! weighting differs for each element (acoustic/elastic/..)
    belmts_mc_loads(:) = elmnts_load(:)

    ! weight nets
    allocate(nwghts(nnodes))
    nwghts(:) = 1
  endif

#ifndef USE_PATOH
  print *, 'PATOH Error: called without PATOH Support. To enable PATOH support, please re-compile with flag -DUSE_PATOH'
  stop 'PATOH not enabled'
#else
  print *, 'patoh: calling patoh with nconst = ', newlvlids, ' constraints'
  !print *, 'patoh: part(1)=', part(1), 'part(2)=',part(2)
  ! patoh partitioning
  call wrap_patoh(nspec, nnodes, nwghts, xpins, pins, belmts_mc_loads, nparts, newlvlids, part, useFixCells, bcuttype, cut)
#endif

  print *,'patoh: Reports cut of ', cut
  print *,'patoh: done'
  print *

  deallocate(xpins,pins)
  deallocate(nwghts,belmts_mc_loads)

  contains

    !----------------------------------------------------------------------------------------------

    subroutine plevel_partitioning()

    use constants, only: PLEVEL_ISLAND

    implicit none
    ! local parameters
    integer :: bmaxplevel,bnumconst,bminel,bmaxel
    integer :: bunitw
    integer,dimension(:),allocatable :: bnumplevel,bplevelids
    integer :: j,ielem,num_elements,ispec,p,lvlid

    ! sets level ids for levels with p-elements
    if (PLEVEL_ISLAND) then
      ! single constraint partitioning
      ! user output
      print *,'p-level island partitioning'

      ! single constraint
      newlvlids = 1
    else
      ! multi-constraint
      ! user output
      print *,'p-level multi-constraint partitioning'

      ! multi-constraint
      newlvlids = num_p_level
    endif

    ! gets maximum p-value of all elements
    bmaxplevel = 0
    do ispec = 1,nspec
      ! checks if p values are valid
      if (ispec_p_refine(ispec) < 1) stop 'ERROR: MAIN: patoh level <= 0'
      ! stores maximum p value
      if (ispec_p_refine(ispec) > bmaxplevel) bmaxplevel = ispec_p_refine(ispec)
    enddo

    allocate(bnumplevel(1:bmaxplevel), stat=ier)
    if (ier /= 0) stop 'ERROR: MAIN: patoh cannot allocate  bnumplevel'
    bnumplevel(:) = 0

    ! counts number of elements for each p-value
    do ispec = 1,nspec
      bnumplevel(ispec_p_refine(ispec)) = bnumplevel(ispec_p_refine(ispec))+1
    enddo

    ! debug
    !print *, 'patoh: elements per levels:',bnumplevel(:)

    !TODO Bora: for the time being we partition all (noneempty) levels among all
    !      processes. Later on we may want to combine some levels in case there are
    !      only a few elems at a level.

    ! counts number of non-zero p-value levels
    bnumconst = 0
    do p = 1,bmaxplevel
      if (bnumplevel(p) > 0) bnumconst = bnumconst + 1
    enddo
    print *, 'patoh: number of p-levels constraints = ', bnumconst

    ! gets minimum/maximum element load values
    bminel = elmnts_load(1)
    bmaxel = elmnts_load(1)
    do ispec = 2, nspec
      if (elmnts_load(ispec) < bminel) then
        bminel = elmnts_load(ispec)
      else if (elmnts_load(ispec) > bmaxel) then
        bmaxel = elmnts_load(ispec)
      endif
    enddo
    print *, 'patoh: element load min = ', bminel, ' max = ', bmaxel

    ! determines if unit weighting (constant element load for all elements)
    bunitw = 0
    if (bminel == bmaxel) bunitw = 1
    print *, 'patoh: unit weighting = ', bunitw

    allocate(bplevelids(1:bmaxplevel), stat=ier)
    if (ier /= 0) stop 'ERROR: MAIN: patoh cannot allocate  bplevelids'
    bplevelids(:) = 0
    ! sets level ids
    lvlid = 1
    do p = 1, bmaxplevel
      if (bnumplevel(p) /= 0) then
        bplevelids(p) = lvlid
        lvlid = lvlid + 1
      endif
    enddo
    lvlid = lvlid - 1 ! Bora: now it is the number of constraints
    if (lvlid /= newlvlids) stop 'ERROR: patoh invalid number of level ids'

    ! element loads
    allocate(belmts_mc_loads(1:newlvlids*nspec), stat = ier)
    if (ier /= 0) stop 'ERROR: MAIN: patoh cannot allocate belmts_mc_loads'
    belmts_mc_loads(:) = 0

    if (PLEVEL_ISLAND) then
      ! element weigths according to p-value
      do ispec = 1,nspec
        belmts_mc_loads(ispec) = ispec_p_refine(ispec) ! * elmnts_load(ispec)
      enddo
    else
      ! sets multi-constrained element loads
      if (bunitw /= 0) then
        ! non-uniform weighting according to p-level
        do ispec = 1,nspec
          ! unit weight for all elements within their p-level set
          belmts_mc_loads((ispec-1) * newlvlids + bplevelids(ispec_p_refine(ispec))) = 1
        enddo
      else
        ! single p-value for all elements
        ! weighting differs for each element according to acoustic/elastic/.. element type
        do ispec = 1,nspec
          ! adds p-value weighted load
          belmts_mc_loads((ispec-1) * newlvlids + bplevelids(ispec_p_refine(ispec))) = elmnts_load(ispec)
        enddo
      endif
    endif

    ! weight nets by p-level. Fine levels generate more
    ! communication between nodes than coarser levels.
    allocate(nwghts(nnodes))
    nwghts(:) = -1
    do j = 1, nnodes
      num_elements = xpins(j+1) - xpins(j)
      do ielem = 1,num_elements
        ! (ielem starts at 1 b/c Fortran 1-index)
        ! +1 for Fortran
        ispec = pins(xpins(j)+ielem)+1
        p = ispec_p_refine(ispec)
        ! set net weight to maximum p of contained cells
        if (p > nwghts(j)) then
          nwghts(j) = p
        endif
      enddo
    enddo
    if (minval(nwghts) <= 0) then
      print *,"Patoh ERROR: net weights > 0:"
      print *,"number of minval=", count(nwghts == minval(nwghts))
      stop 'Error net weights > 0'
    endif

    ! we are trying to ensure that all elements with p>0 are on the same
    ! partition, so we will fix p>1 to partition 0
    useFixCells = 0
    if (PLEVEL_ISLAND) then
      ! single-constraint partitioning
      useFixCells = 1
      part(:) = -1
      do ispec = 1,nspec
        p = ispec_p_refine(ispec)
        if (p > 1) part(ispec) = 0
      enddo
    endif
    print *, 'patoh: use Fix cells  = ', useFixCells

    deallocate(bnumplevel,bplevelids)

    end subroutine plevel_partitioning

  end subroutine partition_patoh
