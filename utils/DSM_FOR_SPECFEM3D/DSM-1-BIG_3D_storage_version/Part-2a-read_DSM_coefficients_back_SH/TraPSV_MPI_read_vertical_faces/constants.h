
! very slow debug mode or not: open many ASCII text files from each node to print debug information.
! *ALWAYS* set this flag to false and recompile the whole code before starting big production runs
  logical, parameter :: SLOW_DEBUG_MODE = .false.  !!! .true.
!
! use timer
  logical, parameter :: USE_TIMER = .true.
!
!  maximum azimuthal order
!!  integer, parameter :: maxlmax = 35001
! write coefficient expansion in disk every maxlmax_g elements
!  integer, parameter :: maxlmax_g = 1000
