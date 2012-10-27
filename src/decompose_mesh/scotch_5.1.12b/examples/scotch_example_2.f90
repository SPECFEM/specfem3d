!***********************************************************
!*                                                        **
!*   NAME       : scotch_example_2.f90                    **
!*                                                        **
!*   AUTHOR     : Francois PELLEGRINI                     **
!*                Frederic COUDERC                        **
!*                                                        **
!*   FUNCTION   : FORTRAN testbed for the LibSCOTCH       **
!*                library routines.                       **
!*                                                        **
!*   DATES      : # Version 5.1  : from : 24 jul 2010     **
!*                                 to     24 jul 2010     **
!*                                                        **
!*   NOTES      : # This program is to be compiled with   **
!*                  the 64-bit version of the libScotch.  **
!*                  It can be deduced from the fact that  **
!*                  SCOTCH_Num values are declared as     **
!*                  "integer*8" quantities.               **
!*                                                        **
!*                # To be compiled with : gfortran        **
!*                  scotch_example_2.f90                  **
!*                  -L/scotch/lib/path/ -lscotch          **
!*                  -lscotcherr -lz -lpthread -lrt        **
!*                                                        **
!234567*****************************************************

PROGRAM Scotch_test

!   USE IFPOSIX

    implicit none

    include "scotchf.h"

    real*8, dimension(Scotch_graphdim) :: Scotchgraph
    real*8, dimension(Scotch_stratdim) :: Scotchstrat

    integer :: ierr, file_u
    integer*8 :: nnt, nelt, i, j, k
    integer*8 :: base, flag, nparts

    integer*8, dimension(2) :: indxtab

    integer*8 :: baseval, vertnbr, vertidx, vendidx, veloidx
    integer*8 :: vlblidx, edgenbr, edgeidx, edloidx

    real*8, dimension(:,:), allocatable :: coordinates
    real*8, dimension(:,:), allocatable :: connectivity
    real*8, dimension(:,:), allocatable :: neighbors

    integer*8, dimension(:), allocatable :: part

    write(6,*) 'Starting'

    call ScotchFstratInit ( Scotchstrat(1), ierr )
    if (ierr /= 0) then
      write(6,*) 'Scotch error : cannot initialize graph'
      STOP
    end if

    call ScotchFgraphInit ( Scotchgraph(1), ierr )
    if (ierr /= 0) then
      write(6,*) 'Scotch error : cannot initialize strat'
      STOP
    end if

    open(10,iostat=ierr,file="bump.grf")

!   call PXFFILENO (10, file_u, ierr)
    file_u = FNUM (10)

    base = 1
    flag = 0
    call ScotchFgraphLoad ( Scotchgraph(1), file_u, base, flag, ierr )
    if (ierr /= 0) then
      write(6,*) 'Scotch error : cannot load graph'
      STOP
    end if

    close(10)

    call ScotchFgraphData ( Scotchgraph(1), indxtab(1), baseval, &
      vertnbr, vertidx, vendidx, veloidx, vlblidx, edgenbr, &
      edgeidx, edloidx )

    write(6,*) baseval, vertnbr, vertidx, vendidx, veloidx, vlblidx, &
      edgenbr, edgeidx, edloidx

    write(6,*) indxtab(vertidx)
    write(6,*) indxtab(vendidx)
    write(6,*) indxtab(veloidx)
    write(6,*) indxtab(vlblidx)
    write(6,*) indxtab(edgeidx)
    write(6,*) indxtab(edloidx)

    call ScotchFgraphCheck ( Scotchgraph(1), ierr )
    if (ierr /= 0) then
      write(6,*) 'ERROR Scotch : Invalid check'
      STOP
    end if

    open (10,iostat=ierr,file="test.grf")

!   call PXFFILENO (10, file_u, ierr)
    file_u = FNUM (10)

    call ScotchFgraphSave ( Scotchgraph(1), file_u, ierr )
    if (ierr /= 0) then
      write(6,*) 'ERROR Scotch : Invalid save '
      STOP
    end if

    close(10)

    nparts = 8
    call ScotchFstratGraphMapChoose ( Scotchstrat(1), &
      SCOTCH_STRATQUALITY, nparts, 0.01, ierr )
    if (ierr /= 0) then
      write(6,*) 'ERROR Scotch : Cannot initialize partitioning&
        &strategy'
      STOP
    end if

    allocate(part(1:vertnbr + 1))

    call ScotchFgraphPart ( Scotchgraph(1), nparts, Scotchstrat(1), &
      part(1), ierr )
    if (ierr /= 0) then
      write(6,*) 'ERROR Scotch : Cannot partition graph'
      STOP
    end if

    do i=1,vertnbr + 1
        write(10,*) i,part(i)
    end do

    call ScotchFgraphExit ( Scotchgraph(1) )

    call ScotchFstratExit ( Scotchstrat(1) )

    write(6,*) 'Complete'

END PROGRAM Scotch_test
