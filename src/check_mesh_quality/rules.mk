#=====================================================================
#
#                         S p e c f e m 3 D
#                         -----------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                              CNRS, France
#                       and Princeton University, USA
#                 (there are currently many more authors!)
#                           (c) October 2017
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================


## compilation directories
S := ${S_TOP}/src/check_mesh_quality
$(check_mesh_quality_OBJECTS): S := ${S_TOP}/src/check_mesh_quality

#######################################

####
#### targets
####

check_mesh_quality_TARGETS = \
	$E/xcheck_mesh_quality \
	$E/xconvert_skewness_to_angle \
	$(EMPTY_MACRO)

check_mesh_quality_OBJECTS = \
	$O/check_mesh_quality.check.o \
	$O/convert_skewness_to_angle.check.o \
	$(EMPTY_MACRO)


#######################################


####
#### rules for executables
####


check: $(check_mesh_quality_TARGETS)
check_mesh: $(check_mesh_quality_TARGETS)


check_mesh_quality: xcheck_mesh_quality
xcheck_mesh_quality: $E/xcheck_mesh_quality

convert_skewness_to_angle: xconvert_skewness_to_angle
xconvert_skewness_to_angle: $E/xconvert_skewness_to_angle

xcheck_mesh_quality_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/count_number_of_sources.shared.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

$E/xcheck_mesh_quality: $O/check_mesh_quality.check.o $(xcheck_mesh_quality_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xcheck_mesh_quality"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""

$E/xconvert_skewness_to_angle: $O/convert_skewness_to_angle.check.o $O/shared_par.shared_module.o
	@echo ""
	@echo "building xconvert_skewness_to_angle"
	@echo ""
	${FCLINK} -o  $E/xconvert_skewness_to_angle $O/convert_skewness_to_angle.check.o $O/shared_par.shared_module.o
	@echo ""

#######################################

####
#### rule to build each .o file below
####

$O/%.check.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<



