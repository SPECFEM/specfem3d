#############################################################################
# read_parameter_cfg.py                                                    
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
import os

def readcfg(filename=None,importmenu=False,mpiflag=False):
    """
    read the configuration file, filename is defined in the command line arguments (see menu.py)
    """
    if importmenu:
        import menu as menu
        cfgname=menu.cfg_name
        id_proc=menu.id_proc
        create_plane=menu.create_plane
        menusurface=menu.surface
        single=menu.single
    elif filename:
        cfgname=filename
        id_proc=0
        menu=False
        create_plane=False
        menusurface=False
        single=False
    else:
        print 'error: no configuration file'
        import sys
        sys.exit()
    #
    from utilities import geo2utm,get_cubit_version #here I can use pyproj but I prefere to include a function in pure python in order to avoid an additional installation
    #
    #
    import ConfigParser
    config = ConfigParser.ConfigParser()
    #
    #
    def converter(s):
        if s == 'True':
            value=True
        elif s == 'False':
            value=False
        elif s == 'None':
            value= None
        else:
            if s.count(',') != 0:
                value=s.split(',')
                while value.count(''):                
                    value.remove('')
            else:
                value=s
            try:
                if type(value).__name__ == 'str':
                    if str(value).count('.') != 0:
                        value=float(value)
                    else:
                        value=int(value)
                else:
                    if str(value).count('.') != 0:
                        value=map(float,value)
                    else:
                        value=map(int,value)
            except:
                pass
        return value
    #
    def section_dict(section):
        dict_o = {}
        options = config.options(section)
        for option in options:
            try:
                value=converter(config.get(section, option))
                dict_o[option] = value
            except:
                dict_o[option] = None
        return dict_o
    
    class attrdict(dict):
        def __init__(self, *args, **kwargs):
            dict.__init__(self, *args, **kwargs)
            self.__dict__ = self
        def __str__(self):
            names=[]
            values=[]
            for name,value in self.items():
                names.append(name)
                values.append(value)
            print names,values
            a=zip(names,values)
            a.sort()
            arc=''
            for o in a:
                if o[0][0] != arc: print
                print o[0],' -> ',o[1]
                arc=o[0][0]
            return '____'
    
    #
    dcfg={}
    #
    #CONSTANTS
    dcfg['osystem'] = 'linux'
    dcfg['debug_cfg']=False
    dcfg['version_cubit']=get_cubit_version()
    dcfg['checkbound']=False
    dcfg['top_partitioner'] = 10000
    dcfg['tres']=0.3 #if n is the vertical component of the normal at a surface pointing horizontally, when -tres < n < tres then the surface is vertical
    dcfg['precision'] =0.02 #precision for the boundary check (0.02 m)
    #
    #INIT
    dcfg['debug']=True 
    dcfg['cubit_info']="on"
    dcfg['echo_info']="on"
    dcfg['jou_info']="on"
    dcfg['jer_info']="on"
    dcfg['monitored_cpu']=0
    dcfg['parallel_import']=True
    dcfg['save_geometry_cubit']=True
    dcfg['save_surface_cubit']=False
    dcfg['save_geometry_paraview'] = False #not implemented 
    dcfg['save_geometry_ACIS'] = False #not implemented
    dcfg['export_exodus_mesh']=False
    dcfg['manual_adj']                  =False
    dcfg['play_adj']                    =False
    dcfg['no_adj']                      =False
    dcfg['nx']=False
    dcfg['ny']=False
    dcfg['nstep']=False
    dcfg['localdir_is_globaldir']=True
    dcfg['refinement_depth']=[]
    dcfg['scratchdir']=None
    dcfg['map_meshing_type']='regularmap'
    dcfg['4sideparallel']=True
    dcfg["outlinebasin_curve"]=False
    dcfg["transition_curve"]=False
    dcfg["faulttrace_curve"]=False
    dcfg['geological_imprint']=False
    dcfg['number_processor_xi']=1
    dcfg['number_processor_eta']=1
    dcfg['filename']=None
    dcfg['actual_vertical_interval_top_layer']=1
    dcfg['coarsening_top_layer']=False
    dcfg['refineinsidevol']=False
    dcfg['sea']=False
    dcfg['seaup']=False
    dcfg['sea_level']=False
    dcfg['sea_threshold']=False
    dcfg['subduction']=True
    dcfg['subduction_thres']=500
    dcfg['debugsurface']=False #if true it creates only the surface not the lofted volumes
    dcfg['lat_orientation']=False
    dcfg['irregulargridded_surf']=False
    dcfg['chktop']=False
    dcfg['smoothing']=False
    dcfg['ntripl']=0
    dcfg['debug_geometry']=False
    
    if float(dcfg['version_cubit']) >= 13.1:
        dcfg['volumecreation_method']=None
    else:
        dcfg['volumecreation_method']='loft'
    
    
    
    dcfg['nsurf'] = None
    if cfgname:
        config.read(cfgname)
        sections=['cubit.options','simulation.cpu_parameters','geometry.surfaces','geometry.volumes','geometry.volumes.layercake','geometry.volumes.flatcake','geometry.volumes.partitioner','geometry.partitioner','meshing']
        
        
        for section in sections:
            try:
                d=section_dict(section)
                dcfg.update(d)
            except:
                pass
        
        if dcfg['nsurf']:
           surface_name=[]
           num_x=[]
           num_y=[]
           xstep=[]
           ystep=[]
           step=[]
           directionx=[]
           directiony=[]
           unit2=[]
           surf_type=[]
           delimiter=[]
           nsurf=int(dcfg['nsurf'])
           for i in range(1,nsurf+1):
               section='surface'+str(i)+'.parameters'
               d=section_dict(section)
               surface_name.append(d['name'])
               surf_type.append(d['surf_type'])
               unit2.append(d['unit_surf'])
               if d['surf_type'] == 'regular_grid':
                   xstep.append(d['step_along_x'])
                   ystep.append(d['step_along_y'])
                   num_x.append(d['number_point_along_x'])
                   num_y.append(d['number_point_along_y'])
               elif d['surf_type'] == 'skin':
                   step.append(d['step'])
                   try:
                      delimiter.append(d['delimiter'])
                   except:
                      pass
                   directionx.append(d['directionx'])
                   directiony.append(d['directiony'])
           dcfg['surface_name']=surface_name
           dcfg['num_x']=num_x       
           dcfg['num_y']=num_y       
           dcfg['xstep']=xstep       
           dcfg['ystep']=ystep       
           dcfg['step']=step        
           dcfg['directionx']=directionx  
           dcfg['directiony']=directiony  
           dcfg['unit2']=unit2       
           dcfg['surf_type']=surf_type   
           dcfg['delimiter']=delimiter  
           
        try:
            tres=0
            xmin,ymin=geo2utm(dcfg['longitude_min'],dcfg['latitude_min'],dcfg['unit'])
            xmax,ymax=geo2utm(dcfg['longitude_max'],dcfg['latitude_max'],dcfg['unit'])
            dcfg['xmin']=xmin
            dcfg['ymin']=ymin
            dcfg['xmax']=xmax
            dcfg['ymax']=ymax
            x1,y1=geo2utm(dcfg['longitude_min'],dcfg['latitude_min'],dcfg['unit'])
            x2,y2=geo2utm(dcfg['longitude_max'],dcfg['latitude_min'],dcfg['unit'])
            x3,y3=geo2utm(dcfg['longitude_max'],dcfg['latitude_max'],dcfg['unit'])
            x4,y4=geo2utm(dcfg['longitude_min'],dcfg['latitude_max'],dcfg['unit'])
            dcfg['x1_box']=x1
            dcfg['y1_box']=y1
            dcfg['x2_box']=x2
            dcfg['y2_box']=y2
            dcfg['x3_box']=x3
            dcfg['y3_box']=y3
            dcfg['x4_box']=x4
            dcfg['y4_box']=y4
            dcfg['tres_boundarydetection']=tres
        except:
            pass
            
        if dcfg['sea']:
            if not dcfg['sea_level']: dcfg['sea_level']=0
            if not dcfg['sea_threshold']: dcfg['sea_threshold']=-200
            dcfg['actual_vertical_interval_top_layer']=1
            dcfg['coarsening_top_layer']=True
        
    
    dcfg['optionsea']={'sea':dcfg['sea'],'seaup':dcfg['seaup'],'sealevel':dcfg['sea_level'],'seathres':dcfg['sea_threshold']}
    cfg=attrdict(dcfg)
    
    if menu:
        try:
            if cfg.working_dir[-1] == '/': cfg.working_dir=cfg.working_dir[:-1]
            if cfg.working_dir[0] != '/': cfg.working_dir='./'+cfg.working_dir
        except:
            cfg.working_dir=os.getcwd()
        
        try:
            if cfg.output_dir[-1] == '/': cfg.output_dir=cfg.output_dir[:-1]
            if cfg.output_dir[0] != '/': cfg.output_dir='./'+cfg.output_dir
        except:
            cfg.output_dir=os.getcwd()
        
        try:
            if cfg.SPECFEM3D_output_dir[-1] == '/': cfg.SPECFEM3D_output_dir=cfg.SPECFEM3D_output_dir[:-1]
            if cfg.SPECFEM3D_output_dir[0] != '/': cfg.SPECFEM3D_output_dir='./'+cfg.SPECFEM3D_output_dir
        except:
            cfg.SPECFEM3D_output_dir=os.getcwd()
        
        cfg.single=single
        
        if menusurface:
            cfg.nsurf=1
            cfg.name=[menu.surface_name]
            cfg.num_x=[menu.num_x]
            cfg.num_y=[menu.num_x]
            cfg.unit=[menu.unit]
            cfg.surf_type=[menu.surf_type]
            try:
                cfg.delimiter=[menu.delimiter]
            except:
                cfg.delimiter=[' ']
            cfg.directionx=[menu.directionx]
            cfg.directiony=[menu.directiony]
    else:
        cfg.SPECFEM3D_output_dir=os.getcwd()
        
        
    if not cfg.number_processor_eta and cfg.nodes:
        cfg.number_processor_xi,cfg.number_processor_eta=split(cfg.nodes)
        
    if isinstance(cfg.filename,str): cfg.filename=[cfg.filename]
    
    try:
        cfg.nproc_eta=cfg.number_processor_eta
        cfg.nproc_xi=cfg.number_processor_xi
        cfg.cpuy=cfg.number_processor_eta
        cfg.cpux=cfg.number_processor_xi
    except:
        pass
    
    if create_plane:
        cfg.x1=map(float,menu.x1.split(','))
        cfg.x2=map(float,menu.x2.split(','))
        cfg.x3=map(float,menu.x3.split(','))
        cfg.x4=map(float,menu.x4.split(','))
        cfg.unit=menu.unit
    #
    if menu:
        cfg.id_proc=menu.id_proc
    else:
        cfg.id_proc=id_proc
    #
    try:
        if isinstance(cfg.tripl,int): cfg.tripl=[cfg.tripl]
    except:
        pass
        
    try:
        if isinstance(cfg.iv_interval,int): cfg.iv_interval=[cfg.iv_interval]
    except:
        pass
        
    try:
        if isinstance(cfg.refinement_depth,int): cfg.refinement_depth=[cfg.refinement_depth]
    except:
        pass
    
    return cfg


class getparameter(dict): 
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self.__dict__ = self


def split(x):
    import math
    c=int(math.sqrt(x))
    while math.fmod(x,c):
        c=c+1
    return c,x/c