! Number of slices for mesh partitioning
integer, parameter  :: nparts = 4

! Useful kind types
integer ,parameter :: short = SELECTED_INT_KIND(4), long = SELECTED_INT_KIND(18)

! Number of nodes per elements.
integer, parameter  :: ESIZE = 8

! Number of faces per element.
integer, parameter  :: nfaces = 6

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9
