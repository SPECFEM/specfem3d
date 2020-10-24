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

! METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes,
! and producing fill reducing orderings for sparse matrices.
!
! The algorithms implemented in METIS are based on the multilevel recursive-bisection,
! multilevel k-way, and multi-constraint partitioning schemes.
!
! developed by George Karypis' Lab
! http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
!
! reference:
! Karypis, G. & Kumar, V. (1999). "A fast and high quality multilevel scheme for partitioning irregular graphs".
! SIAM Journal on Scientific Computing. 20 (1): 359.
! https://epubs.siam.org/doi/10.1137/S1064827595287997
!
!
! support for METIS:
!  If you want to use METIS instead of SCOTCH, then the program xdecompose_mesh must be compiled with METIS support:
!
!    - Add the path to the METIS library files and add compiler flag -DUSE_METIS in Makefile in the SPECFEM3D/ directory
!    - re-compile xdecompose_mesh: > make xdecompose_mesh
!
!  To select METIS instead of SCOTCH, in Par_file choose:
!    PARTITIONING_TYPE = 2

  subroutine partition_metis()

  use decompose_mesh_par

  implicit none

  integer :: edgecut
  integer :: nconstraints
  real(4), dimension(:), allocatable :: ubvec
  integer, dimension(:,:), allocatable :: vwgt
  integer, dimension(:), allocatable :: adjwgt

  ! checks if partitioning needed
  if (nparts == 1) return

  ! METIS partitioning
  print *,'METIS partitioning'

  ! initializes
  edgecut = 0

  if (LTS_MODE) then
    ! LTS partitioning
    call plevel_partioning()
  else
    ! no-LTS partitioning
    nconstraints = 1

    ! single edge weights only
    ! safety check
    if (num_p_level /= 1) stop 'Error Metis partioning: no-LTS call has multiple p-levels'

    allocate(vwgt(1,nspec)) ! vertices weights
    allocate(adjwgt(1))     ! dummy edges weights
    allocate(ubvec(1))      ! dummy maximum load imbalance
    adjwgt(:) = 0
    ubvec(:) = 0
    vwgt(1,:) = elmnts_load(:)
  endif

  ! METIS partioning
  print *, "Metis: nconstraints == ",nconstraints

  call wrap_metis(nspec,nconstraints,xadj,adjncy,vwgt,adjwgt,nparts,ubvec,edgecut,part)

  ! frees arrays
  deallocate(vwgt,adjwgt,ubvec)

  print *
  print *,'total edgecut created by METIS = ',edgecut
  print *

  ! check for error
  if (edgecut < 0) stop 'METIS Error: negative edgecut'

  contains

    !----------------------------------------------------------------------------------------------

    subroutine plevel_partioning()

    use constants, only: PLEVEL_ISLAND

    implicit none

    ! LTS
    integer :: iconn,num_neighbors,vtx_id,p_neighbor,p_me,neighbor
    integer :: ispec,p,pmax,ilevel
    integer,dimension(:),allocatable :: ilevel_from_p

    ! METIS partitioning

    ! constraints
    if (PLEVEL_ISLAND) then
      ! single constraint partitioning
      ! user output
      print *,'  p-level island partitioning'

      ! constraints
      nconstraints = 1
    else
      ! multi-constraint
      ! user output
      print *,'  p-level multi-constraint partitioning'

      ! multi-constraint
      nconstraints = num_p_level
    endif

    ! lookup table for p values
    pmax = maxval(p_level) ! maximum p values
    allocate(ilevel_from_p(pmax))
    ! for example: set {1,2,8) -> num_p_level = 3
    !              has ilevel_from_p = (3,2,0,0,0,0,0,1)
    ilevel_from_p(:) = 0
    do ilevel = 1,num_p_level
      p = p_level(ilevel)
      ilevel_from_p(p) = ilevel
    enddo
    ! debug
    !print *,'p_level       = ',p_level
    !print *,'ilevel_from_p = ',ilevel_from_p

    print *, "  size adjncy = ",size(adjncy),"  minval adjncy = ",minval(adjncy)
    print *, "  minval xadj = ",minval(xadj)
    print *, "  maxval xadj = ",maxval(xadj)


    if (PLEVEL_ISLAND) then
      ! single constraint partitioning
      allocate(vwgt(1,nspec))           ! vertices weights
    else
      ! multi-constraint
      allocate(vwgt(num_p_level,nspec)) ! vertices weights
    endif
    allocate(adjwgt(maxval(xadj)))    ! edges weights
    allocate(ubvec(num_p_level))      ! maximum load imbalance

    ! sets graph vertices weights
    vwgt(:,:) = 0
    ! hard coded for now...
    do ispec = 1,nspec
      ! p refinement value
      p = ispec_p_refine(ispec)

      if (PLEVEL_ISLAND) then
        ! single constraint
        ! weighting according to p value
        vwgt(1,ispec) = p
      else
        ! multi-constraint
        ! determines level index
        ! example: p = 1 -> ilevel = 0+1 = 1; p = 2 -> ilevel = 1+1 = 2; p = 4 -> ilevel = 2+1 = 3
        !ilevel = int(log(real(p))/log(2.0))+1
        ! this will not work if we have p values which skip a level, e.g., a set like {1,2,8} -> num_p_level = 3
        ! thus we use the lookup table
        ilevel = ilevel_from_p(p)
        if (ilevel < 1 .or. ilevel > num_p_level) then
          print *,'Error: p = ',p,'ilevel = ',ilevel,'num_p_level = ',num_p_level
          stop 'Error invalid ilevel for vwgt found'
        endif
        ! weighting according to p value
        vwgt(ilevel,ispec) = p
      endif
    enddo
    deallocate(ilevel_from_p)

    ! sets maximum load imbalance constraint
    do ilevel = 1,num_p_level
      ubvec(ilevel) = 1.02     ! 1.02 indicates a maximum load imbalance of 2%
    enddo

    ! weight graph edges by p-level. Fine levels generate more
    ! communication between elements than coarser levels.
    if (PLEVEL_ISLAND) then
      adjwgt(:) = 1
    else
      adjwgt(:) = 0
    endif
    do ispec = 1,nspec
      num_neighbors = xadj(ispec+1)-xadj(ispec)
      ! loops over connections to other elements
      do iconn = 1,num_neighbors
        ! vertex id
        vtx_id = xadj(ispec) + iconn

        ! get neighbor (+1 for C->Fortran)
        neighbor = adjncy(vtx_id)+1

        ! p-values
        p_neighbor = ispec_p_refine(neighbor)
        p_me = ispec_p_refine(ispec)

        if (PLEVEL_ISLAND) then
          ! set very high cut weight to avoid cutting p-levels
          if (p_me > 1 .or. p_neighbor > 1) then
            adjwgt(vtx_id) = 1000
          endif
        else
          ! initially just set edge to my value
          adjwgt(vtx_id) = p_me
          ! check neighbor
          if (adjwgt(vtx_id) < p_neighbor) then
            ! neighbor is costlier
            adjwgt(vtx_id) = p_neighbor
          endif
        endif
      enddo
    enddo
    if (minval(adjwgt) <= 0) then
      print *, "Metis Error: vertex weights must be > 0:"
      print *, "number of minval=", count(adjwgt == minval(adjwgt))
      stop
    endif

    end subroutine plevel_partioning

  end subroutine partition_metis
