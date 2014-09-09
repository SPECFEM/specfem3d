#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  2 . 1
#               ---------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, July 2012
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
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
S := ${S_TOP}/src/check_mesh_quality_CUBIT_Abaqus
$(check_mesh_quality_CUBIT_Abaqus_OBJECTS): S := ${S_TOP}/src/check_mesh_quality_CUBIT_Abaqus

#######################################

####
#### targets
####

check_mesh_quality_CUBIT_Abaqus_TARGETS = \
	$E/xcheck_mesh_quality_CUBIT_Abaqus \
	$E/xconvert_skewness_to_angle \
	$E/xmultiply_CUBIT_Abaqus_mesh_by_1000 \
	$(EMPTY_MACRO)

check_mesh_quality_CUBIT_Abaqus_OBJECTS = \
	$O/check_mesh_quality_CUBIT_Abaqus.check.o \
	$O/convert_skewness_to_angle.check.o \
	$O/unused_mod.shared_module.o \
	$O/multiply_CUBIT_Abaqus_mesh_by_1000.check.o \
	$(EMPTY_MACRO)


#######################################


####
#### rules for executables
####

check: $(check_mesh_quality_CUBIT_Abaqus_TARGETS)


check_mesh_quality_CUBIT_Abaqus: xcheck_mesh_quality_CUBIT_Abaqus
xcheck_mesh_quality_CUBIT_Abaqus: $E/xcheck_mesh_quality_CUBIT_Abaqus

convert_skewness_to_angle: xconvert_skewness_to_angle
xconvert_skewness_to_angle: $E/xconvert_skewness_to_angle

multiply_CUBIT_Abaqus_mesh_by_1000: xmultiply_CUBIT_Abaqus_mesh_by_1000
xmultiply_CUBIT_Abaqus_mesh_by_1000: $E/xmultiply_CUBIT_Abaqus_mesh_by_1000

# rules for the pure Fortran version

$E/xcheck_mesh_quality_CUBIT_Abaqus: $O/check_mesh_quality_CUBIT_Abaqus.check.o
	${FCLINK} -o  $E/xcheck_mesh_quality_CUBIT_Abaqus $O/check_mesh_quality_CUBIT_Abaqus.check.o

$E/xconvert_skewness_to_angle: $O/convert_skewness_to_angle.check.o
	${FCLINK} -o  $E/xconvert_skewness_to_angle $O/convert_skewness_to_angle.check.o

$E/xmultiply_CUBIT_Abaqus_mesh_by_1000: $O/multiply_CUBIT_Abaqus_mesh_by_1000.check.o
	${FCLINK} -o  $E/xmultiply_CUBIT_Abaqus_mesh_by_1000 $O/multiply_CUBIT_Abaqus_mesh_by_1000.check.o


#######################################

####
#### rule to build each .o file below
####

$O/%.check.o: $S/%.f90 $O/constants_mod.shared_module.o  $O/unused_mod.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<



