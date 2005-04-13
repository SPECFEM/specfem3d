#=====================================================================
#
#          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 2
#          --------------------------------------------------
#
#                 Dimitri Komatitsch and Jeroen Tromp
#    Seismological Laboratory - California Institute of Technology
#         (c) California Institute of Technology July 2004
#
#    A signed non-commercial agreement is required to use this program.
#   Please check http://www.gps.caltech.edu/research/jtromp for details.
#           Free for non-commercial academic research ONLY.
#      This program is distributed WITHOUT ANY WARRANTY whatsoever.
#      Do not redistribute this program without written permission.
#
#=====================================================================
#
# Copyright July 2004, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship Acknowledged.
#
# Any commercial use must be negotiated with the Office of Technology
# Transfer at the California Institute of Technology. This software may be
# subject to U.S. export control laws and regulations. By accepting
# this software, the user agrees to comply with all applicable U.S. export laws
# and regulations, including the International Traffic and Arms Regulations,
# 22 C.F.R. 120-130 and the Export Administration Regulations,
# 15 C.F.R. 730-744. User has the responsibility to obtain export licenses,
# or other export authority as may be required before exporting such
# information to foreign countries or providing access to foreign nationals.
# In no event shall the California Institute of Technology be liable to any
# party for direct, indirect, special, incidental or consequential damages,
# including lost profits, arising out of the use of this software and its
# documentation, even if the California Institute of Technology has been
# advised of the possibility of such damage.
#
# The California Institute of Technology specifically disclaims any
# warranties, including the implied warranties or merchantability and fitness
# for a particular purpose. The software and documentation provided hereunder
# is on an "as is" basis, and the California Institute of Technology has no
# obligations to provide maintenance, support, updates, enhancements or
# modifications.
#

O = obj

################ PC Linux #################
#
# Beowulf Portland pgf90
#
F90 = pgf90
MPIF90 = mpif90 #/home/local/mpich/bin/mpif90
FLAGS_CHECK = -fast -Mnobounds -Mneginfo -Mdclchk -Mstandard -Knoieee
FLAGS_NO_CHECK = -fast -Mnobounds -Mneginfo -Mdclchk -Munroll=c:6 -Mstandard -Knoieee
MPI_FLAGS = 

#
# Intel ifort Fortran90 for Linux
#
#F90 = ifort
#MPIF90 = mpif90
#
# Caltech cluster  (Hrothgar)
#
#FLAGS_NO_CHECK = -O3 -tpp6 -xK  -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -check nobounds  #-std95 -e95
#
# more recent machines (e.g. Pangu)
#
#FLAGS_NO_CHECK = -fast -tpp7 -xN  -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds # -e95
#
# debug with range checking
#
#FLAGS_NO_CHECK = -O0 -static -e95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check bounds
#FLAGS_CHECK = $(FLAGS_NO_CHECK)
#MPI_FLAGS = -Vaxlib

#
# g95 (free f95 compiler from http://www.g95.org, still under development, but works)
#
#F90 = g95
#MPIF90 = g95
#FLAGS_CHECK = -O
#FLAGS_NO_CHECK = -O
#MPI_FLAGS =

#
# AbSoft
#
#F90 = f90
#MPIF90 = mpif90
#FLAGS_CHECK = -W132 -s -O2 -cpu:p7 -v -YDEALLOC=ALL
#FLAGS_NO_CHECK = $(FLAGS_CHECK)
#MPI_FLAGS =   

#
# NAG compiler for Linux
#
#F90 = f95
#MPIF90 = mpif90
#FLAGS_CHECK = -O -u -strict95 -C=all
#FLAGS_NO_CHECK = -O -u -strict95
#MPI_FLAGS = 

#
# Lahey f90
#
#F90 = lf95
#MPIF90 = mpif90
#FLAGS_CHECK = --warn --wo --tpp --f95 --dal -O --chk
#FLAGS_NO_CHECK = --warn --wo --tpp --f95 --dal -O
#MPI_FLAGS = 

################ SGI Irix #################
##
##  CAUTION: always define setenv TRAP_FPE OFF on SGI before compiling
##
#F90 = f90
#MPIF90 = f90
#FLAGS_NO_CHECK = -ansi -u -64 -O3 -OPT:Olimit=0 -OPT:roundoff=3 -OPT:IEEE_arithmetic=3 -r10000 -mips4
#FLAGS_CHECK = $(FLAGS_NO_CHECK) -check_bounds
#MPI_FLAGS = -lmpi -lfastm -lfpe

################## Compaq Dec Alpha #################
#F90 = f90
#MPIF90 = f90
#FLAGS_NO_CHECK = -fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow
#FLAGS_CHECK = $(FLAGS_NO_CHECK) -check bounds
#MPI_FLAGS = -lfmpi -lmpi

################## Earth Simulator and NEC SX-5 ##################
#F90 = f90
#MPIF90 = f90
#FLAGS_CHECK = -C hopt -R2 -Wf" -L nostdout noinclist mrgmsg noeject -msg b -pvctl loopcnt=5000000 expand=6 fullmsg vecthreshold=20 -s" -pi auto line=100 exp=swap_all,rank
#FLAGS_NO_CHECK = $(FLAGS_CHECK)
#MPI_FLAGS =

######## IBM SP or Power 4 - file names need to be *.f instead of *.f90 ######
#F90 = mpxlf_r
#MPIF90 = mpxlf_r
#FLAGS_CHECK = -O3 -Q -qarch=auto -qcache=auto -qtune=auto -qlanglvl=95pure -qmaxmem=65536 -qflag=L:L -qhalt=L -qsuffix=f=f90
#FLAGS_CHECK = -q64 -O3 -qarch=pwr4 -qlanglvl=95pure -qflag=L:L -qhalt=L -qsuffix=f=f90
# use this on IDRIS machines, www.idris.fr
#FLAGS_CHECK = -q64 -O4 -qfree=f90 -qsuffix=f=f90
#FLAGS_NO_CHECK = $(FLAGS_CHECK)
#MPI_FLAGS = 

baksave:
	cp *f90 *h README_SPECFEM3D_BASIN DATA/Par_file* Makefile go_mesher go_solver mymachines bak

meshfem3D: constants.h \
       $O/meshfem3D.o \
       $O/create_regions_mesh.o \
       $O/calc_jacobian.o \
       $O/gll_library.o \
       $O/get_jacobian_boundaries.o \
       $O/get_absorb.o \
       $O/get_flags_boundaries.o \
       $O/get_MPI_cutplanes_xi.o \
       $O/get_MPI_cutplanes_eta.o \
       $O/get_global.o \
       $O/write_AVS_DX_global_faces_data.o \
       $O/write_AVS_DX_surface_data.o \
       $O/write_AVS_DX_global_data.o \
       $O/write_AVS_DX_mesh_quality_data.o \
       $O/create_name_database.o \
       $O/define_subregions_basin.o \
       $O/define_subregions_heuristic.o \
       $O/get_shape3D.o \
       $O/get_shape2D.o \
       $O/hex_nodes.o \
       $O/lagrange_poly.o \
       $O/mesh_vertical.o \
       $O/numerical_recipes.o \
       $O/interpolate_gocad_block_MR.o \
       $O/interpolate_gocad_block_HR.o \
       $O/salton_trough_gocad.o \
       $O/socal_model.o \
       $O/aniso_model.o \
       $O/compute_rho_estimate.o \
       $O/hauksson_model.o \
       $O/save_arrays.o \
       $O/save_header_file.o \
       $O/read_basin_topo_bathy_file.o \
       $O/read_moho_map.o \
       $O/exit_mpi.o \
       $O/read_parameter_file.o \
       $O/read_value_parameters.o \
       $O/utm_geo.o \
       $O/compute_parameters.o
## use MPI here
	${MPIF90} $(FLAGS_CHECK) -o xmeshfem3D \
       $O/meshfem3D.o \
       $O/create_regions_mesh.o \
       $O/calc_jacobian.o \
       $O/gll_library.o \
       $O/get_jacobian_boundaries.o \
       $O/get_absorb.o \
       $O/get_flags_boundaries.o \
       $O/get_MPI_cutplanes_xi.o \
       $O/get_MPI_cutplanes_eta.o \
       $O/get_global.o \
       $O/write_AVS_DX_global_faces_data.o \
       $O/write_AVS_DX_surface_data.o \
       $O/write_AVS_DX_global_data.o \
       $O/write_AVS_DX_mesh_quality_data.o \
       $O/create_name_database.o \
       $O/define_subregions_basin.o \
       $O/define_subregions_heuristic.o \
       $O/get_shape3D.o \
       $O/get_shape2D.o \
       $O/hex_nodes.o \
       $O/lagrange_poly.o \
       $O/mesh_vertical.o \
       $O/numerical_recipes.o \
       $O/interpolate_gocad_block_MR.o \
       $O/interpolate_gocad_block_HR.o \
       $O/salton_trough_gocad.o \
       $O/socal_model.o \
       $O/aniso_model.o \
       $O/compute_rho_estimate.o \
       $O/hauksson_model.o \
       $O/save_arrays.o \
       $O/save_header_file.o \
       $O/read_basin_topo_bathy_file.o \
       $O/read_moho_map.o \
       $O/exit_mpi.o \
       $O/read_parameter_file.o \
       $O/read_value_parameters.o \
       $O/utm_geo.o \
       $O/compute_parameters.o $(MPI_FLAGS)

# solver also depends on values from mesher
specfem3D: constants.h OUTPUT_FILES/values_from_mesher.h \
       $O/specfem3D.o \
       $O/read_arrays_solver.o \
       $O/calc_jacobian.o \
       $O/gll_library.o \
       $O/get_cmt.o \
       $O/numerical_recipes.o \
       $O/write_seismograms.o \
       $O/read_parameter_file.o \
       $O/read_value_parameters.o \
       $O/utm_geo.o \
       $O/compute_parameters.o \
       $O/locate_source.o \
       $O/locate_receivers.o \
       $O/get_global.o \
       $O/comp_source_time_function.o \
       $O/recompute_jacobian.o \
       $O/hex_nodes.o \
       $O/lagrange_poly.o \
       $O/read_basin_topo_bathy_file.o \
       $O/exit_mpi.o \
       $O/get_shape3D.o \
       $O/create_name_database.o \
       $O/read_arrays_buffers_solver.o \
       $O/define_derivation_matrices.o \
       $O/compute_arrays_source.o \
       $O/get_attenuation_model.o \
       $O/assemble_MPI_vector.o \
       $O/assemble_MPI_scalar.o
## use MPI here
	${MPIF90} $(FLAGS_NO_CHECK) -o xspecfem3D \
       $O/specfem3D.o \
       $O/read_arrays_solver.o \
       $O/calc_jacobian.o \
       $O/gll_library.o \
       $O/get_cmt.o \
       $O/numerical_recipes.o \
       $O/write_seismograms.o \
       $O/read_parameter_file.o \
       $O/read_value_parameters.o \
       $O/utm_geo.o \
       $O/compute_parameters.o \
       $O/locate_source.o \
       $O/locate_receivers.o \
       $O/get_global.o \
       $O/comp_source_time_function.o \
       $O/recompute_jacobian.o \
       $O/hex_nodes.o \
       $O/lagrange_poly.o \
       $O/read_basin_topo_bathy_file.o \
       $O/exit_mpi.o \
       $O/get_shape3D.o \
       $O/create_name_database.o \
       $O/read_arrays_buffers_solver.o \
       $O/define_derivation_matrices.o \
       $O/compute_arrays_source.o \
       $O/get_attenuation_model.o \
       $O/assemble_MPI_vector.o \
       $O/assemble_MPI_scalar.o $(MPI_FLAGS)

convolve_source_timefunction: $O/convolve_source_timefunction.o
	${F90} $(FLAGS_CHECK) -o xconvolve_source_timefunction $O/convolve_source_timefunction.o

create_header_file: $O/create_header_file.o $O/read_parameter_file.o \
     $O/compute_parameters.o $O/save_header_file.o $O/utm_geo.o $O/read_value_parameters.o
	${F90} $(FLAGS_CHECK) -o xcreate_header_file $O/create_header_file.o \
     $O/read_parameter_file.o $O/compute_parameters.o $O/save_header_file.o $O/utm_geo.o $O/read_value_parameters.o

create_movie_AVS_DX: $O/create_movie_AVS_DX.o $O/read_parameter_file.o \
     $O/compute_parameters.o $O/utm_geo.o $O/read_value_parameters.o
	${F90} $(FLAGS_CHECK) -o xcreate_movie_AVS_DX $O/create_movie_AVS_DX.o \
     $O/read_parameter_file.o $O/compute_parameters.o $O/utm_geo.o $O/read_value_parameters.o

combine_AVS_DX: constants.h $O/combine_AVS_DX.o $O/get_cmt.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o
	${F90} $(FLAGS_CHECK) -o xcombine_AVS_DX $O/combine_AVS_DX.o $O/get_cmt.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o

check_mesh_quality_AVS_DX: constants.h $O/check_mesh_quality_AVS_DX.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o
	${F90} $(FLAGS_CHECK) -o xcheck_mesh_quality_AVS_DX $O/check_mesh_quality_AVS_DX.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o

check_buffers_2D: constants.h $O/check_buffers_2D.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o
	${F90} $(FLAGS_CHECK) -o xcheck_buffers_2D $O/check_buffers_2D.o \
       $O/read_parameter_file.o $O/compute_parameters.o $O/create_serial_name_database.o $O/utm_geo.o $O/read_value_parameters.o

combine_paraview_data: constants.h $O/combine_paraview_data.o $O/write_c_binary.o
	${F90} $(FLAGS_CHECK) -o xcombine_paraview_data  $O/combine_paraview_data.o $O/write_c_binary.o

clean:
	rm -f $O/*.o *.o *.gnu OUTPUT_FILES/timestamp* OUTPUT_FILES/starttime*txt work.pc* xmeshfem3D xspecfem3D xcombine_AVS_DX xcheck_mesh_quality_AVS_DX xcheck_buffers_2D xconvolve_source_timefunction xcreate_header_file xcreate_movie_AVS_DX xcombine_paraview_data
####
#### rule for each .o file below
####

###
### optimized flags and dependence on values from mesher here
###

### use MPI here
$O/specfem3D.o: constants.h OUTPUT_FILES/values_from_mesher.h specfem3D.f90
	${MPIF90} $(FLAGS_NO_CHECK) -c -o $O/specfem3D.o specfem3D.f90

### use MPI here
$O/assemble_MPI_vector.o: constants.h OUTPUT_FILES/values_from_mesher.h assemble_MPI_vector.f90
	${MPIF90} $(FLAGS_NO_CHECK) -c -o $O/assemble_MPI_vector.o assemble_MPI_vector.f90

### use MPI here
$O/assemble_MPI_scalar.o: constants.h OUTPUT_FILES/values_from_mesher.h assemble_MPI_scalar.f90
	${MPIF90} $(FLAGS_NO_CHECK) -c -o $O/assemble_MPI_scalar.o assemble_MPI_scalar.f90

$O/read_arrays_solver.o: constants.h OUTPUT_FILES/values_from_mesher.h read_arrays_solver.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_arrays_solver.o read_arrays_solver.f90

###
### regular compilation options here
###

$O/convolve_source_timefunction.o: convolve_source_timefunction.f90
	${F90} $(FLAGS_CHECK) -c -o $O/convolve_source_timefunction.o convolve_source_timefunction.f90

$O/create_header_file.o: create_header_file.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_header_file.o create_header_file.f90

$O/combine_AVS_DX.o: constants.h combine_AVS_DX.f90
	${F90} $(FLAGS_CHECK) -c -o $O/combine_AVS_DX.o combine_AVS_DX.f90

$O/save_header_file.o: constants.h save_header_file.f90
	${F90} $(FLAGS_CHECK) -c -o $O/save_header_file.o save_header_file.f90

$O/check_mesh_quality_AVS_DX.o: constants.h check_mesh_quality_AVS_DX.f90
	${F90} $(FLAGS_CHECK) -c -o $O/check_mesh_quality_AVS_DX.o check_mesh_quality_AVS_DX.f90

$O/check_buffers_2D.o: constants.h check_buffers_2D.f90
	${F90} $(FLAGS_CHECK) -c -o $O/check_buffers_2D.o check_buffers_2D.f90

### use MPI here
$O/locate_source.o: constants.h locate_source.f90
	${MPIF90} $(FLAGS_CHECK) -c -o $O/locate_source.o locate_source.f90

### use MPI here
$O/locate_receivers.o: constants.h locate_receivers.f90
	${MPIF90} $(FLAGS_CHECK) -c -o $O/locate_receivers.o locate_receivers.f90

## use MPI here
$O/exit_mpi.o: constants.h exit_mpi.f90
	${MPIF90} $(FLAGS_CHECK) -c -o $O/exit_mpi.o exit_mpi.f90

$O/read_parameter_file.o: constants.h read_parameter_file.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_parameter_file.o read_parameter_file.f90

$O/read_value_parameters.o: constants.h read_value_parameters.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_value_parameters.o read_value_parameters.f90

$O/utm_geo.o: constants.h utm_geo.f90
	${F90} $(FLAGS_CHECK) -c -o $O/utm_geo.o utm_geo.f90

$O/compute_parameters.o: constants.h compute_parameters.f90
	${F90} $(FLAGS_CHECK) -c -o $O/compute_parameters.o compute_parameters.f90

$O/calc_jacobian.o: constants.h calc_jacobian.f90
	${F90} $(FLAGS_CHECK) -c -o $O/calc_jacobian.o calc_jacobian.f90

$O/gll_library.o: constants.h gll_library.f90
	${F90} $(FLAGS_CHECK) -c -o $O/gll_library.o gll_library.f90

$O/get_absorb.o: constants.h get_absorb.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_absorb.o get_absorb.f90

$O/get_jacobian_boundaries.o: constants.h get_jacobian_boundaries.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_jacobian_boundaries.o get_jacobian_boundaries.f90

$O/get_flags_boundaries.o: constants.h get_flags_boundaries.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_flags_boundaries.o get_flags_boundaries.f90

$O/get_MPI_cutplanes_xi.o: constants.h get_MPI_cutplanes_xi.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_MPI_cutplanes_xi.o get_MPI_cutplanes_xi.f90

$O/get_MPI_cutplanes_eta.o: constants.h get_MPI_cutplanes_eta.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_MPI_cutplanes_eta.o get_MPI_cutplanes_eta.f90

$O/get_cmt.o: constants.h get_cmt.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_cmt.o get_cmt.f90

$O/create_movie_AVS_DX.o: constants.h create_movie_AVS_DX.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_movie_AVS_DX.o create_movie_AVS_DX.f90

$O/get_global.o: constants.h get_global.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_global.o get_global.f90

$O/write_AVS_DX_global_faces_data.o: constants.h write_AVS_DX_global_faces_data.f90
	${F90} $(FLAGS_CHECK) -c -o $O/write_AVS_DX_global_faces_data.o write_AVS_DX_global_faces_data.f90

$O/write_AVS_DX_surface_data.o: constants.h write_AVS_DX_surface_data.f90
	${F90} $(FLAGS_CHECK) -c -o $O/write_AVS_DX_surface_data.o write_AVS_DX_surface_data.f90

$O/write_AVS_DX_global_data.o: constants.h write_AVS_DX_global_data.f90
	${F90} $(FLAGS_CHECK) -c -o $O/write_AVS_DX_global_data.o write_AVS_DX_global_data.f90

$O/write_AVS_DX_mesh_quality_data.o: constants.h write_AVS_DX_mesh_quality_data.f90
	${F90} $(FLAGS_CHECK) -c -o $O/write_AVS_DX_mesh_quality_data.o write_AVS_DX_mesh_quality_data.f90

$O/get_shape3D.o: constants.h get_shape3D.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_shape3D.o get_shape3D.f90

$O/get_shape2D.o: constants.h get_shape2D.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_shape2D.o get_shape2D.f90

$O/hex_nodes.o: constants.h hex_nodes.f90
	${F90} $(FLAGS_CHECK) -c -o $O/hex_nodes.o hex_nodes.f90

$O/mesh_vertical.o: constants.h mesh_vertical.f90
	${F90} $(FLAGS_CHECK) -c -o $O/mesh_vertical.o mesh_vertical.f90

## use MPI here
$O/meshfem3D.o: constants.h meshfem3D.f90
	${MPIF90} $(FLAGS_CHECK) -c -o $O/meshfem3D.o meshfem3D.f90

$O/numerical_recipes.o: constants.h numerical_recipes.f90
	${F90} $(FLAGS_CHECK) -c -o $O/numerical_recipes.o numerical_recipes.f90

$O/interpolate_gocad_block_MR.o: constants.h interpolate_gocad_block_MR.f90
	${F90} $(FLAGS_CHECK) -c -o $O/interpolate_gocad_block_MR.o interpolate_gocad_block_MR.f90

$O/interpolate_gocad_block_HR.o: constants.h interpolate_gocad_block_HR.f90
	${F90} $(FLAGS_CHECK) -c -o $O/interpolate_gocad_block_HR.o interpolate_gocad_block_HR.f90

$O/salton_trough_gocad.o: constants.h salton_trough_gocad.f90
	${F90} $(FLAGS_CHECK) -c -o $O/salton_trough_gocad.o salton_trough_gocad.f90

$O/socal_model.o: constants.h socal_model.f90
	${F90} $(FLAGS_CHECK) -c -o $O/socal_model.o socal_model.f90

$O/aniso_model.o: constants.h aniso_model.f90
	${F90} $(FLAGS_CHECK) -c -o $O/aniso_model.o aniso_model.f90

$O/compute_rho_estimate.o: constants.h compute_rho_estimate.f90
	${F90} $(FLAGS_CHECK) -c -o $O/compute_rho_estimate.o compute_rho_estimate.f90

$O/hauksson_model.o: constants.h hauksson_model.f90
	${F90} $(FLAGS_CHECK) -c -o $O/hauksson_model.o hauksson_model.f90

$O/save_arrays.o: constants.h save_arrays.f90
	${F90} $(FLAGS_CHECK) -c -o $O/save_arrays.o save_arrays.f90

$O/comp_source_time_function.o: constants.h comp_source_time_function.f90
	${F90} $(FLAGS_CHECK) -c -o $O/comp_source_time_function.o comp_source_time_function.f90

$O/read_basin_topo_bathy_file.o: constants.h read_basin_topo_bathy_file.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_basin_topo_bathy_file.o read_basin_topo_bathy_file.f90

$O/read_moho_map.o: constants.h read_moho_map.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_moho_map.o read_moho_map.f90

$O/write_seismograms.o: constants.h write_seismograms.f90
	${F90} $(FLAGS_CHECK) -c -o $O/write_seismograms.o write_seismograms.f90

$O/lagrange_poly.o: constants.h lagrange_poly.f90
	${F90} $(FLAGS_CHECK) -c -o $O/lagrange_poly.o lagrange_poly.f90

$O/recompute_jacobian.o: constants.h recompute_jacobian.f90
	${F90} $(FLAGS_CHECK) -c -o $O/recompute_jacobian.o recompute_jacobian.f90

$O/create_regions_mesh.o: constants.h create_regions_mesh.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_regions_mesh.o create_regions_mesh.f90

$O/create_name_database.o: constants.h create_name_database.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_name_database.o create_name_database.f90

$O/create_serial_name_database.o: constants.h create_serial_name_database.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_serial_name_database.o create_serial_name_database.f90

$O/define_subregions_basin.o: constants.h define_subregions_basin.f90
	${F90} $(FLAGS_CHECK) -c -o $O/define_subregions_basin.o define_subregions_basin.f90

$O/define_subregions_heuristic.o: constants.h define_subregions_heuristic.f90
	${F90} $(FLAGS_CHECK) -c -o $O/define_subregions_heuristic.o define_subregions_heuristic.f90

$O/read_arrays_buffers_solver.o: constants.h read_arrays_buffers_solver.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_arrays_buffers_solver.o read_arrays_buffers_solver.f90

$O/define_derivation_matrices.o: constants.h define_derivation_matrices.f90
	${F90} $(FLAGS_CHECK) -c -o $O/define_derivation_matrices.o define_derivation_matrices.f90

$O/compute_arrays_source.o: constants.h compute_arrays_source.f90
	${F90} $(FLAGS_CHECK) -c -o $O/compute_arrays_source.o compute_arrays_source.f90

$O/get_attenuation_model.o: constants.h get_attenuation_model.f90
	${F90} $(FLAGS_CHECK) -c -o $O/get_attenuation_model.o get_attenuation_model.f90

$O/combine_paraview_data.o: constants.h combine_paraview_data.f90
	${F90} $(FLAGS_CHECK) -c -o $O/combine_paraview_data.o combine_paraview_data.f90

$O/write_c_binary.o: write_c_binary.c
	cc -c -o $O/write_c_binary.o write_c_binary.c


