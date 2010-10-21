
pgf90 -fast -o xdens densify_hauksson_ori.f90
xdens < hauksson_new_format.dat > hauksson_denser.dat

pgf90 -fast -o xregrid regrid_hauksson_regular.f90 
xregrid < hauksson_denser.dat > hauksson_final_grid_raw.dat

pgf90 -fast -o xshow show_hauksson_regrid_dx.f90
./xshow < hauksson_final_grid_raw.dat > final_regular_raw.dx 

pgf90 -fast -o xsmooth smooth_final_hauksson.f90
./xsmooth < hauksson_final_grid_raw.dat > hauksson_final_grid_smooth.dat

pgf90 -fast -o xshow show_hauksson_regrid_dx.f90
./xshow < hauksson_final_grid_smooth.dat > final_regular_smooth.dx 

