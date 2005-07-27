
!! DK DK the 21 elastic coefficients and density for copper (cubic system, three parameters)

! rotation angle for vertical free surface
!! DK DK UGLY  double precision, parameter :: ANGLE_MESH = -45.d0
!! DK DK UGLY  double precision, parameter :: ANGLE_MESH = 0.d0
  double precision, parameter :: ANGLE_MESH = -45.d0
  double precision, parameter :: ANGLE_ROTATE = ANGLE_MESH * PI / 180.d0

!! DK DK
!! DK DK copper crystal (cubic symmetry)
!! DK DK
  double precision, parameter :: rho_copper = 8920.d0

!! DK DK regular c_ijkl with no rotation
  double precision, parameter :: c11_copper = 169.d9
  double precision, parameter :: c12_copper = 122.d9
  double precision, parameter :: c13_copper = c12_copper
  double precision, parameter :: c14_copper = 0.d0
  double precision, parameter :: c15_copper = 0.d0
  double precision, parameter :: c16_copper = 0.d0

  double precision, parameter :: c22_copper = c11_copper
  double precision, parameter :: c23_copper = c12_copper
  double precision, parameter :: c24_copper = 0.d0
  double precision, parameter :: c25_copper = 0.d0
  double precision, parameter :: c26_copper = 0.d0

  double precision, parameter :: c33_copper = c11_copper
  double precision, parameter :: c34_copper = 0.d0
  double precision, parameter :: c35_copper = 0.d0
  double precision, parameter :: c36_copper = 0.d0

  double precision, parameter :: c44_copper = 75.3d9
  double precision, parameter :: c45_copper = 0.d0
  double precision, parameter :: c46_copper = 0.d0

  double precision, parameter :: c55_copper = c44_copper
  double precision, parameter :: c56_copper = 0.d0

  double precision, parameter :: c66_copper = c44_copper

!! DK DK new c_ijkl rotated by 45 degrees using Helbig's rotation routine
!  double precision, parameter :: c11_copper =     0.22080000000000000000d12
!  double precision, parameter :: c12_copper =     0.70200000000000000000d11
!  double precision, parameter :: c13_copper =     0.12200000000000000000d12
!  double precision, parameter :: c14_copper =     0.d0
!  double precision, parameter :: c15_copper =     0.d0
!  double precision, parameter :: c16_copper =     0.d0
!  double precision, parameter :: c22_copper =     0.22080000000000000000d12
!  double precision, parameter :: c23_copper =     0.12200000000000000000d12
!  double precision, parameter :: c24_copper =     0.d0
!  double precision, parameter :: c25_copper =     0.d0
!  double precision, parameter :: c26_copper =     0.d0
!  double precision, parameter :: c33_copper =     0.16900000000000000000d12
!  double precision, parameter :: c34_copper =     0.d0
!  double precision, parameter :: c35_copper =     0.d0
!  double precision, parameter :: c36_copper =     0.d0
!  double precision, parameter :: c44_copper =     0.75300000000000000000d11
!  double precision, parameter :: c45_copper =     0.d0
!  double precision, parameter :: c46_copper =     0.d0
!  double precision, parameter :: c55_copper =     0.75300000000000000000d11
!  double precision, parameter :: c56_copper =     0.d0
!  double precision, parameter :: c66_copper =     0.23500000000000000000d11


!! DK DK
!! DK DK isotropic equivalent for copper for tests with anisotropy off
!! DK DK
!!!!!!!!!!!!!! XXXXXXXXXXXXXX YYYYYYYYYYYYYY UGLY  double precision, parameter :: cp_copper = 3570.d0
!!!!!!!!!!!!!! XXXXXXXXXXXXXX YYYYYYYYYYYYYY UGLY  double precision, parameter :: cs_copper = cp_copper / 2.d0

!!! DK DK example of Mesaverde clay shale (transversely anisotropic, hexagonal)
!!! DK DK taken from Carcione et al., Geophysics, vol. 57,  p. 1595 (1992)
!!! DK DK and Komatitsch, Barnes and Tromp, Geophysics (2000)
!!! DK DK in standard units (N.m-2)
!
!  double precision, parameter :: rho_copper = 2590.d0
!
!  double precision, parameter :: c11_copper = 66.6e9
!  double precision, parameter :: c22_copper = c11_copper
!  double precision, parameter :: c33_copper = 39.9e9
!  double precision, parameter :: c44_copper = 10.9e9
!  double precision, parameter :: c55_copper = c44_copper
!  double precision, parameter :: c13_copper = 39.4e9
!  double precision, parameter :: c12_copper = 19.7e9
!  double precision, parameter :: c66_copper = - (c12_copper - c11_copper)/2.d0
!  double precision, parameter :: c23_copper = c13_copper
!  double precision, parameter :: c14_copper = 0.d0
!  double precision, parameter :: c24_copper = 0.d0
!  double precision, parameter :: c34_copper = 0.d0
!  double precision, parameter :: c15_copper = 0.d0
!  double precision, parameter :: c25_copper = 0.d0
!  double precision, parameter :: c35_copper = 0.d0
!  double precision, parameter :: c45_copper = 0.d0
!  double precision, parameter :: c16_copper = 0.d0
!  double precision, parameter :: c26_copper = 0.d0
!  double precision, parameter :: c36_copper = 0.d0
!  double precision, parameter :: c46_copper = 0.d0
!  double precision, parameter :: c56_copper = 0.d0

! source time function:  1 = Gaussian   2 = Ricker
  integer, parameter :: SOURCE_TIME_FUNCTION = 1

! source dominant frequency (not used if Heaviside)
  double precision, parameter :: SOURCE_DOMINANT_FREQ = 1.6d6

! scaling factor for point source
  double precision, parameter :: FACTOR_SOURCE = 1.d6

! save postscript snapshots of the simulation
  logical, parameter :: POSTSCRIPT_SNAPSHOTS = .true.

! interval at which we output polarograms
  integer, parameter :: ISAMP_POLAROGRAMS = 5

