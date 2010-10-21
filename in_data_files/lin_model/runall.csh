
pgf90 -fast -o xdens densify_hauksson_ori.f90
xdens < lin_new_format.dat > lin_denser.dat

pgf90 -fast -o xregrid regrid_hauksson_regular.f90 
xregrid < lin_denser.dat > lin_final_grid_raw.dat

pgf90 -fast -o xsmooth smooth_final_hauksson.f90
./xsmooth < lin_final_grid_raw.dat > lin_final_grid_smooth.dat
