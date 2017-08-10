
# fill the fortran compiler name :
FC=ifort


$FC create_uniform_stations.f90 -o create_uniform_stations.x
$FC create_movie_slice.f90 -o create_movie_slice.x
$FC filtre_traces.f90 filters.f90 -o filtre_traces.x
$FC create_interface.f90 -o create_interface.x
$FC stalta.f90 -o stalta.x
$FC pick_list_stalta.f90 -o pick_list_stalta.x
$FC apply_tapper_on_interface.f90 -o apply_tapper_on_interface.x
#$FC create_model.f90 -o create_model.x
