module mtadj_variables

  use mtadj_constants

  ! ------------------------------------------------------------
  ! parameter file
  integer :: iker  ! kernel type
  integer :: itap  ! taper type
  character(len=150) :: meas_dir ! output dir for measurements
  character(len=150) :: adj_dir ! output dir for adj src

  ! preprocessing
  logical :: BANDPASS  ! band-pass option
  real :: tshort, tlong ! band-pass range
  real :: fstart, fend

  ! measurements
  real :: BEFORE_SHIFT   ! tshift_cc max control
  real :: wtr ! water level for estimating f_right through amp. of s(f)
  ! multi-taper measurements
  real :: npi
  integer :: ntaper ! ntaper = 2*npi
  real :: wtr_mtm ! water level for bottom of MTM
  real :: MIN_DT_SIGMA,MIN_DlnA_SIGMA ! minimum dt,dlnA standard error

  ! compute error ????

  ! window selection after measurements
  logical :: SELECT_WINDOW
  ! selection criteria
  integer :: ncycle_in_window
  ! dtau_i < min(T_i/dt_fac,dt_max_scale*tshift), err_dt_i < T_i/err_fac
  real :: dt_fac, dt_max_scale, err_fac
  !  real ::  after_quality,after_tshift ! ??????

  ! write adjoint source
  logical :: INCLUDE_ERROR, BANDPASS_ADJ
  real :: b_adj, dt_adj
  integer :: npts_adj   !interpolation of adjoint source

  ! -----------------------------------------------------
  ! global parameters

  ! data and syn array
  real, dimension(NDIM) :: data, syn
  real :: b, dt
  integer :: npts

  ! windowed data and syn
  real, dimension(NPT) :: dataw, synw
  integer :: nstart, nend, nlen

  ! cross-correlation measurements
  real :: tshift_cc, dlnA
  real :: sigma_tshift_cc, sigma_dlnA_cc ! error estimates

  ! freq-dependent measurements
  real,dimension(NPT,NMAX_TAPER) :: tas

  ! left and right of frequency range to output measurements
  integer :: i_left, i_right, idf_fd
  real :: f_left, f_right, df_fd ! freq-dep indep df spacing
  ! maximum power index
  integer :: i_pmax

  ! complex transfer function from freq-dep measurements
  complex, dimension(NPT) :: trans_fdm
  real, dimension(NPT) :: dtau_fdm, dlnA_fdm ! dtau and dlnA
  real, dimension(NPT) :: sigma_dtau_fdm, sigma_dlnA_fdm ! error bars

  ! adjoint source
  integer :: ipwr_t, ipwr_w

end module mtadj_variables
