
! variables for numerical integrations using Simpson's or Gauss' quadrature rule for mass lumping

!! DK DK number of integration points in the sub-grid used for mass lumping
!! DK DK is USE_GAUSS_INTEGRATION below is false, you can use any even value you want for ns.
!! DK DK is USE_GAUSS_INTEGRATION below is true, you can use ns = 32, 64 or 128.
!! DK DK
!! DK DK in the DSM code, Simpson integration uses points going from 0 to 2*ns
!! DK DK but Gauss integration only uses points going from 1 to ns.
!! DK DK Thus, to use roughly the same number of points (with a difference of 1)
!! DK DK when using Gauss you should use twice the value of ns compared to the same run with Simpson (e.g. use 256 instead of 128)
!! DK DK
! integer, parameter :: ns = 128 ! use this for Simpson
  integer, parameter :: ns = 32 ! use this for Gauss

!! DK DK use Gauss integration instead, or keep using the older Simpson's integration routine
  logical, parameter :: USE_GAUSS_INTEGRATION = .true.

