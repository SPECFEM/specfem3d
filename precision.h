!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 1
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology October 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
! solver in single or double precision depending on the machine
!
!  ALSO CHANGE FILE  constants.h ACCORDINGLY
!
! uncomment this to run in single precision
! integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
! uncomment this to run in double precision
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_DOUBLE_PRECISION

