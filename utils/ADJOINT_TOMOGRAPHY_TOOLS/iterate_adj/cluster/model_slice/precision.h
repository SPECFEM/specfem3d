!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology July 2005
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! precision.h.  Generated from precision.h.in by configure.

!
! solver in single or double precision depending on the machine
!
! set to MPI_REAL to run in single precision
! set to MPI_DOUBLE_PRECISION to run in double precision
!
! ALSO CHANGE FILE constants.h ACCORDINGLY
!
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
  integer, parameter :: CUSTOM_MPI_2REAL = MPI_2REAL
