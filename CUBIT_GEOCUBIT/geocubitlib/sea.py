#############################################################################
# sea.py                                                    
# this file is part of GEOCUBIT                                             #
#                                                                           #
# Created by Emanuele Casarotti                                             #
# Copyright (c) 2008 Istituto Nazionale di Geofisica e Vulcanologia         #
#                                                                           #
#############################################################################
#                                                                           #
# GEOCUBIT is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# GEOCUBIT is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with GEOCUBIT.  If not, see <http://www.gnu.org/licenses/>.         #
#                                                                           #
#############################################################################

def adjust_sea_layers(zvertex,sealevel,bathymetry,cfg):
    if cfg.seaup: 
        if sealevel:
            vertex=max(cfg.sea_level,vertex-cfg.sea_threshold)    
        elif bathymetry:
            #if zvertex > cfg.sea_threshold: zvertex=max(zvertex,cfg.sea_level)+cfg.sea_threshold #move node below the topography of sea_threshold
            vertex=vertex
    else:
        if sealevel and zvertex < cfg.sea_level:
            zvertex=cfg.sea_level
        elif bathymetry:
            if zvertex > cfg.sea_threshold: zvertex=max(zvertex,cfg.sea_level)+cfg.sea_threshold #move node below the topography of sea_threshold
    return zvertex
    


