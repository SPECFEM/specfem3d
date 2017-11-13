#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  3 . 0
#               ---------------------------------------
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

DIR = specfem3D

# The rest of this file is generic
#######################################

####
#### targets
####

default:
	$(MAKE) -C ../.. $(DIR)

all:
	$(MAKE) -C ../.. all

clean:
	$(MAKE) -C ../.. CLEAN=$(DIR) clean

cleanall:
	$(MAKE) -C ../.. clean

backup:
	mkdir -p bak
	cp *f90 *h Makefile bak

bak: backup

.PHONY: default all clean cleanall backup bak


