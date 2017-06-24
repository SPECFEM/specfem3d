module grid3d_variables

  use grid3d_constants

  ! GRID3D.PAR

  ! derivative files
  character(len=150) :: cmt_file, new_cmt_file
  real :: dmoment

  ! data selection
  character(len=150) :: flexwin_out_file
  logical :: weigh_data_files, read_weight
  real :: comp_z_weight, comp_r_weight, comp_t_weight, &
       az_exp_weight, &
       pnl_dist_weight, rayleigh_dist_weight, love_dist_weight

  ! grid search scheme
  logical :: station_correction,global_search
  integer :: ncalc
  real ::  tshift_max,s_strike,e_strike,d_strike,s_dip,e_dip,d_dip, &
       s_rake,e_rake,d_rake,s_mw,e_mw,d_mw

  ! misc
  logical ::write_new_cmt


  ! global search interval
  real :: t_strike,t_dip,t_rake,t_mw

  ! flexwin output variables
  integer :: nfiles, nwin_total, nwins(NRECMAX)

  ! data weights array
  real :: data_weights(NWINMAX)

  ! data and syn files/arrays
  character(len=150) :: data_file, syn_file
  real, dimension(NDATAMAX) :: data, syn
  real, dimension(NM,NDATAMAX) :: dsyn

  ! mij and misfit arrays (these are the largest arrays)
  real, dimension(NM,NMEC_MAX) :: mij
  real, dimension(NMEC_MAX) :: misfit

end module grid3d_variables
