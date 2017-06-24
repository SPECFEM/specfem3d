program grid3d_flexwin

  use grid3d_sub
  implicit none

  character(len=150) par_file

  ! read and print parameters
  par_file = 'GRID3D.PAR'
  print *
  print *, 'Reading grid3d parameters ...'
  call set_parameters(par_file)

 ! read data/syn files and setup data weights
  print *, 'Compute weights associated with each window ...'
  call setup_data_weights

  ! compute waveform differences to determine the best mechanism
  print *, 'Grid search over strike, dip and rake for minimum misfit value ...'
  call grid_search


end program grid3d_flexwin
