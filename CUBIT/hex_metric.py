#############################################################################
# hex_metric.py                                                    
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
##mesh=SEM_metric_3D
##mesh.check_metric()
##print mesh
#
#
import math

try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass

class SEM_metric_3D(object):
    def __init__(self,list_hex=None,volume=None):
        super(SEM_metric_3D, self).__init__()
        self.list_hex=list_hex
        self.volume=volume
        self.skew_hystogram=None
        self.max_skewness=None
        self.max_angle=None
        self.min_angle=None
        self.min_edge_length=None
        self.cubit_skew=0
        self.max_valence=None
        self.node_with_max_valence=None
        self.valence_threshold=7
        self.resolution=.05
        self.resolution_edge_length=.1
        self.nbin=10
        self.spatial_unit=1
        self.gll=0.17267316464601141
        self.dt=None
    def __repr__(self):
        if self.max_skewness is not None:
            self.skew_hystogram=self.hyst(0,self.max_skewness,self.skew_hyst)
            print '-'*70
            print 'SKEWNESS'
            print 
            if len(self.hex_max_skewness) <= 30:
                print 'max = ',self.max_skewness,' in hexes ',self.hex_max_skewness
                print '(angle -> minimun =', self.min_angle, ' maximun =', self.max_angle,')'
            else:
                print 'max = ',self.max_skewness,' in ', len(self.hex_max_skewness), ' hexes '
                print '(angle -> minimun =', self.min_angle, ' maximun =', self.max_angle,')'
            print 
            print 'skew hystogram'
            print 
            tot=0
            for i in self.skew_hystogram.values():
                tot=tot+len(i)
            #k=self.skew_hystogram.keys()
            #k.sort()
            factor=self.max_skewness/self.nbin
            for i in range(0,self.nbin+1):
                if self.skew_hystogram.has_key(i):
                    if (i+1)*factor <= 1:
                        print i,' [',i*factor,'->',(i+1)*factor,'[ : ',len(self.skew_hystogram[i]),'/',tot,' hexes (',len(self.skew_hystogram[i])/float(tot)*100.,'%)'
                else:
                    if (i+1)*factor <= 1:
                        print i,' [',i*factor,'->',(i+1)*factor,'[ : ',0,'/',tot,' hexes (0%)'
            print
        ###############################################
        if self.min_edge_length is not None:
            self.edgemin_hystogram=self.hyst(self.min_edge_length,self.max_edge_length,self.edgemin_hyst)
            self.edgemax_hystogram=self.hyst(self.min_edge_length,self.max_edge_length,self.edgemax_hyst)
            print '-'*70
            print 'edge length'
            print 
            if len(self.hex_min_edge_length) <= 30:
                print 'minimum edge length: ', self.min_edge_length, ' in hexes ',  self.hex_min_edge_length
            else:
                print 'minimum edge length: ', self.min_edge_length, ' in ',  len(self.hex_min_edge_length), ' hexes.'
            if len(self.hex_max_edge_length) <= 30:
                print 'maximum edge length: ', self.max_edge_length, ' in hexes ',  self.hex_max_edge_length              
            else:                                                                                                                        
                print 'maximum edge length: ', self.max_edge_length, ' in ',  len(self.hex_max_edge_length), ' hexes.'    
            print
            print 'edge length hystogram'
            print 
            factor=(self.max_edge_length-self.min_edge_length)/self.nbin
            print 'minimum edge length'
            tot=0
            for i in self.edgemin_hystogram.values():
                tot=tot+len(i)
            #k=self.edgemin_hystogram.keys()
            #k.sort()
            for i in range(0,self.nbin+1):
                if self.edgemin_hystogram.has_key(i):
                    print i,' [',i*factor+self.min_edge_length,'->',(i+1)*factor+self.min_edge_length,'[ : ',len(self.edgemin_hystogram[i]),'/',tot,' hexes (',len(self.edgemin_hystogram[i])/float(tot)*100.,'%)'
                else:
                    print i,' [',i*factor+self.min_edge_length,'->',(i+1)*factor+self.min_edge_length,'[ : ',0,'/',tot,' hexes (0%)'
            print 
            print 'maximum edge length'
            tot=0
            for i in self.edgemax_hystogram.values():
                tot=tot+len(i)
            #k=self.edgemax_hystogram.keys()
            #k.sort()
            for i in range(0,self.nbin+1):
                if self.edgemax_hystogram.has_key(i):
                    print i,' [',i*factor+self.min_edge_length,'->',(i+1)*factor+self.min_edge_length,'[ : ',len(self.edgemax_hystogram[i]),'/',tot,' hexes (',len(self.edgemax_hystogram[i])/float(tot)*100.,'%)'
                else:
                    print i,' [',i*factor+self.min_edge_length,'->',(i+1)*factor+self.min_edge_length,'[ : ',0,'/',tot,' hexes (0%)'
        if self.dt is not None:
            print '-'*70
            print
            print 'STABILITY'
            print
            print 'time step < ',self.dt,'s, for velocity = ',self.velocity
            print 
        try:
            return str(len(self.list_hex))+' hexes checked'
        except:
            return 'please performe a metric check: ex ... object.check_metric()'
    #
    #
    def invert_dict(self,d):
         inv = {}
         for k,v in d.iteritems():
             keys = inv.setdefault(v, [])
             keys.append(k)
         return inv
    def hex_metric(self,h):
        nodes=cubit.get_connectivity('Hex',h)
        equiangle_skewness=None
        min_angle=float('infinity')
        edge_length_min=float('infinity')
        max_angle=None
        edge_length_max=None
        #loopnode=[(i,j) for i in range(8) for j in range(i+1,8)]
        if len(nodes) == 4:
            faces=[[0,1,2,3]]
        elif len(nodes) == 8:
            faces=[[0,1,2,3],[4,5,6,7],[1,5,6,2],[0,4,7,3],[2,6,7,3],[1,5,4,0]]
        else:
            print 'bad definition of nodes'
            return None,None,None
        x=[]
        y=[]
        z=[]
        for i in nodes:
            n=cubit.get_nodal_coordinates(i)
            x.append(n[0])
            y.append(n[1])
            z.append(n[2])
        for face in faces:
            for i in range(-1,3):
                vx1=x[face[i-1]]-x[face[i]]
                vy1=y[face[i-1]]-y[face[i]]
                vz1=z[face[i-1]]-z[face[i]]
                #
                vx2=x[face[i+1]]-x[face[i]]
                vy2=y[face[i+1]]-y[face[i]]
                vz2=z[face[i+1]]-z[face[i]]
                #
                norm1=math.sqrt(vx1*vx1+vy1*vy1+vz1*vz1)
                norm2=math.sqrt(vx2*vx2+vy2*vy2+vz2*vz2)
                #
                if norm1 == 0 or norm2 == 0:
                    print 'degenerated mesh, 0 length edge'
                    import sys
                    sys.exit()
                angle=math.acos((vx1*vx2+vy1*vy2+vz1*vz2)/(norm1*norm2))
                #
                equiangle_skewness=max(equiangle_skewness,math.fabs(2.*angle-math.pi)/math.pi)
                min_angle=min(min_angle,angle*180./math.pi)
                max_angle=max(max_angle,angle*180./math.pi)
                edge_length_min=min(edge_length_min,norm1)
                edge_length_max=max(edge_length_max,norm1)
        return round(equiangle_skewness,5),min_angle,max_angle,round(edge_length_min,5),round(edge_length_max,5)
    def show(self,minvalue,maxvalue,dic):
        hexlist=[]
        for k in dic.keys():
            if minvalue <= dic[k] <= maxvalue: hexlist.extend([k])
        #command = "draw hex "+str(hexlist)
        #command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        #cubit.cmd(command)
        return hexlist
    def hyst(self,minvalue,maxvalue,dic):
        if maxvalue != minvalue:
            factor=(maxvalue-minvalue)/self.nbin
        else:
            factor=1
        dic_new={}
        for k in dic.keys():
            dic_new[k]=int((dic[k]-minvalue)/factor)
        inv_dic=self.invert_dict(dic_new)
        return inv_dic                                                                                          
    #
    #
    #
    #
    #
    ############################################################################################
    def pick_hex(self,list_hex=None,volume=None):
        if self.cubit_skew > 0:
            command = "del group skew_top"
            cubit.cmd(command)
            command = "group 'skew_top' add quality volume all skew low "+str(self.cubit_skew)
            cubit.silent_cmd(command)
            group=cubit.get_id_from_name("skew_top")
            self.list_hex=cubit.get_group_hexes(group)
        elif list_hex is not None:
            self.list_hex=list_hex
        elif volume is not None:
            command = "group 'hextmp' add hex in volume "+str(volume)    
            cubit.silent_cmd(command)                                           
            group=cubit.get_id_from_name("hextmp")                       
            self.list_hex=cubit.get_group_hexes(group)
            command = "del group hextmp"
            cubit.silent_cmd(command)   
        elif self.volume is not None:
            command = "group 'hextmp' add hex in volume "+str(self.volume) 
            cubit.silent_cmd(command)                                        
            group=cubit.get_id_from_name("hextmp")                    
            self.list_hex=cubit.get_group_hexes(group)
            command = "del group hextmp"
            cubit.silent_cmd(command)                 
        elif list_hex is None and self.list_hex is None:
            self.list_hex=cubit.parse_cubit_list('hex','all')
        print 'list_hex: ',len(self.list_hex),' hexes'
    #
    #
    ############################################################################################
    #SKEWNESS
    def check_skew(self,list_hex=None,volume=None):
        self.pick_hex(list_hex=list_hex,volume=volume)
        global_max_skew=None
        global_max_angle=None
        global_min_angle=float('infinity')
        tmp=int(1./self.resolution)+1
        skew_hyst={}
        for i in range(0,tmp):
            if i*self.resolution < 1:
                skew_hyst[i]=[]
        h_global_max_skew=[]
        h_global_max_angle=[]
        h_global_min_angle=[]
        for h in self.list_hex:
            equiangle_skewness,min_angle,max_angle,tmp,tmp2=self.hex_metric(h)
            skew_hyst[h]=equiangle_skewness
            #
            if equiangle_skewness > global_max_skew: 
                global_max_skew=equiangle_skewness
                h_global_max_shew=[]
                h_global_max_shew.append(h)
                global_min_angle=min_angle
                global_max_angle=max_angle
            elif equiangle_skewness == global_max_skew:
                h_global_max_shew.append(h)
            s=h_global_max_shew
        self.skew_hyst=skew_hyst
        self.max_skewness=global_max_skew
        self.max_angle=global_max_angle
        self.min_angle=global_min_angle
        self.hex_max_skewness=s
    def show_skew(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_skewness,self.skew_hyst)
        elif low is not None and high is None:
            hexlist=self.show(0,high,self.skew_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.skew_hyst)
        else:
            hexlist=self.hex_min_edge_length
        command = "draw hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.silent_cmd(command)
    def list_skew(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_skewness,self.skew_hyst)
        elif low is not None and high is None:
            hexlist=self.show(0,high,self.skew_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.skew_hyst)
        else:
            hexlist=self.hex_min_edge_length
        return hexlist
    def highlight_skew(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_skewness,self.skew_hyst)
        elif low is not None and high is None:
            hexlist=self.show(0,high,self.skew_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.skew_hyst)
        else:
            hexlist=self.hex_min_edge_length
        command = "highlight hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.silent_cmd(command)
    #
    #
    #
    #############################################################################################
    #HEX VALENCE
    def check_valence(self):
        #usage: hex_valence()
        #
        list_hex=cubit.parse_cubit_list('hex','all')
        list_node=cubit.parse_cubit_list('node','all')
        lookup=dict(zip(list_node,[0]*len(list_node)))
        for h in list_hex:
            n=cubit.get_connectivity('Hex',h)
            for i in n:
                if lookup.has_key(i): lookup[i]=lookup[i]+1
        inv_lookup=self.invert_dict(lookup)
        #
        vmax=max(inv_lookup.keys())
        nmax=inv_lookup[max(inv_lookup.keys())]
        print 'max hex valence ',vmax,' at node ',nmax
        print '_____'
        for v in inv_lookup.keys():
            print ('valence %2i - %9i hexes '% (v,len(inv_lookup[v])))
        self.valence_summary=inv_lookup
        self.max_valence=vmax
        self.node_with_max_valence=nmax
    #
    #
    ###########################################################################################
    #EDGE LENGTH
    def check_edge_length(self,list_hex=None,volume=None):
        self.pick_hex(list_hex=list_hex,volume=volume)
        global_min_edge_length=float('infinity')
        h_global_min_edge_length=[]
        global_max_edge_length=None
        h_global_max_edge_length=[]
        edgemin_hyst={}
        edgemax_hyst={}     
        for h in self.list_hex:
            tmp,tmp2,tmp3,edge_length_min,edge_length_max=self.hex_metric(h)
            edgemax_hyst[h]=edge_length_max
            edgemin_hyst[h]=edge_length_min
            if edge_length_min < global_min_edge_length:
                global_min_edge_length=edge_length_min
                h_global_min_edge_length=[]
                h_global_min_edge_length.append(h)
            elif edge_length_min == global_min_edge_length:
                h_global_min_edge_length.append(h)
            if edge_length_max > global_max_edge_length:
                global_max_edge_length=edge_length_max
                h_global_max_edge_length=[]
                h_global_max_edge_length.append(h)      
            elif edge_length_max == global_max_edge_length:    
                h_global_max_edge_length.append(h)
        self.min_edge_length=global_min_edge_length
        self.hex_min_edge_length=h_global_min_edge_length
        self.max_edge_length=global_max_edge_length
        self.hex_max_edge_length=h_global_max_edge_length
        self.edgemin_hyst=edgemin_hyst
        self.edgemax_hyst=edgemax_hyst        
    def show_edgemin(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemin_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemin_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemin_hyst)
        else:
            hexlist=self.hex_min_edge_length
        command = "draw hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.cmd(command)
    def show_edgemax(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemax_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemax_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemax_hyst)
        else:
            hexlist=self.hex_max_edge_length
        command = "draw hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.cmd(command)
    def list_edgemax(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemax_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemax_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemax_hyst)
        else:
            hexlist=self.hex_max_edge_length        
        return hexlist
    def list_edgemin(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemin_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemin_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemin_hyst)
        else:
            hexlist=self.hex_min_edge_length
        return hexlist
    def highlight_edgemin(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemin_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemin_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemin_hyst)
        else:
            hexlist=self.hex_min_edge_length
        command = "highlight hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.cmd(command)
    def highlight_edgemax(self,low=None,high=None):
        if low is None and high is not None:
            hexlist=self.show(low,self.max_edge_length,self.edgemax_hyst)
        elif low is not None and high is None:
            hexlist=self.show(self.min_edge_length,high,self.edgemax_hyst)
        if low is not None and high is not None:
            hexlist=self.show(low,high,self.edgemax_hyst)
        else:
            hexlist=self.hex_max_edge_length
        command = "highlight hex "+str(hexlist)
        command = command.replace("["," ").replace("]"," ").replace("("," ").replace(")"," ")
        cubit.cmd(command)    #
    #
    #
    #
    #
    ############################################################################################
    def check_metric(self,list_hex=None,volume=None):
        self.pick_hex(list_hex=list_hex,volume=volume)
        skew_hyst={}
        edata={}
        edgemin_hyst={}
        edgemax_hyst={}
        edgemaxdata={}
        edgemindata={}
        #
        global_min_edge_length=float('infinity')
        h_global_min_edge_length=[]
        global_max_edge_length=None
        h_global_max_edge_length=[]
        #
        global_max_skew=None
        global_max_angle=None
        global_min_angle=float('infinity')
        h_global_max_skew=[]
        h_global_max_angle=[]
        h_global_min_angle=[]
        #
        s=0.
        for h in self.list_hex:
            equiangle_skewness,min_angle,max_angle,edge_length_min,edge_length_max=self.hex_metric(h)
            skew_hyst[h]=equiangle_skewness
            edgemax_hyst[h]=edge_length_max
            edgemin_hyst[h]=edge_length_min
            if edge_length_min < global_min_edge_length:
                global_min_edge_length=edge_length_min
                h_global_min_edge_length=[]
                h_global_min_edge_length.append(h)
            elif edge_length_min == global_min_edge_length:
                h_global_min_edge_length.append(h)
            if edge_length_max > global_max_edge_length:
                global_max_edge_length=edge_length_max
                h_global_max_edge_length=[]
                h_global_max_edge_length.append(h)     
            elif edge_length_max == global_max_edge_length:   
                h_global_max_edge_length.append(h)
            if equiangle_skewness > global_max_skew: 
                global_max_skew=equiangle_skewness
                h_global_max_shew=[]
                h_global_max_shew.append(h)
                global_min_angle=min_angle
                global_max_angle=max_angle
            elif equiangle_skewness == global_max_skew:
                h_global_max_shew.append(h)
            s=h_global_max_shew
        #
        self.max_skewness=global_max_skew    
        self.max_angle=global_max_angle      
        self.min_angle=global_min_angle      
        self.hex_max_skewness=s
        self.skew_hyst=skew_hyst
        #                
        self.min_edge_length=global_min_edge_length
        self.hex_min_edge_length=h_global_min_edge_length
        self.max_edge_length=global_max_edge_length
        self.hex_max_edge_length=h_global_max_edge_length
        self.edgemin_hyst=edgemin_hyst
        self.edgemax_hyst=edgemax_hyst
        #
        
class SEM_stability_3D(SEM_metric_3D):
    def __init__(self,list_hex=None,volume=None):
        super(SEM_metric_3D, self).__init__()
        self.list_hex=list_hex
        self.volume=volume
        self.Ngll_per_wavelength=5
        self.gllcoeff=.17
        self.maxgllcoeff=0.5-.17
        self.Cmax=.3
        self.nbin=10
        self.period_hyst=None
        self.dt=None
        self.cubit_skew =0
        #
        #
        #
    def __repr__(self):
        print 'check mesh stability'
        #
    def check_simulation_parameter(self,list_hex=None,volume=None,tomofile=None,vp_static=None,vs_static=None):
        self.pick_hex(list_hex=list_hex,volume=volume)
        timestep_hyst={}
        stability_hyst={}
        global_min_timestep=float(1)
        global_max_periodresolved=None
        if tomofile: 
            self.read_tomo(tomofile)
        #
        #
        for ind,h in enumerate(self.list_hex):
            if ind%10000==0: print 'hex checked: '+str(int(float(ind)/len(self.list_hex)*100))+'%'
            dt_tmp,pmax_tmp,_,_=self.hex_simulation_parameter(h,vp_static=vp_static,vs_static=vs_static)
            timestep_hyst[h]=dt_tmp
            stability_hyst[h]=pmax_tmp
            if dt_tmp < global_min_timestep:
                global_min_timestep=dt_tmp
            if pmax_tmp > global_max_periodresolved:
                global_max_periodresolved=pmax_tmp   
        #
        self.period=global_max_periodresolved    
        self.dt=global_min_timestep      
        self.dt_hyst=timestep_hyst
        self.period_hyst=stability_hyst
        #                
    def read_tomo(self,tomofile=None):
        if tomofile:
            print 'reading tomography file ',tomofile
            import numpy
            #xtomo,ytomo,ztomo,vp,vs,rho=numpy.loadtxt(tomofile,skiprows=4)
            print 'tomography file loaded'
            tf=open(tomofile,'r')
            orig_x, orig_y, orig_z, end_x, end_y, end_z        =map(float,tf.readline().split())
            spacing_x, spacing_y, spacing_z                    =map(float,tf.readline().split())
            nx, ny, nz                                         =map(int,tf.readline().split())
            vp_min, vp_max, vs_min, vs_max, rho_min, rho_max   =map(float,tf.readline().split())
            #
            ind=0
            import sys
            xtomo,ytomo,ztomo,vp,vs=[],[],[],[],[]
            while ind<nx*ny*nz:
                ind=ind+1
                if ind%100000==0: sys.stdout.write("reading progress: %i/%i   \r" % (ind,(nx*ny*nz)) )
                x,y,z,v1,v2,_=map(float,tf.readline().split())
                xtomo.append(x)
                ytomo.append(y)
                ztomo.append(z)
                vp.append(v1)
                vs.append(v2)
            tf.close()
            print 'tomography file loaded'
            #
            self.orig_x, self.orig_y, self.orig_z, self.end_x, self.end_y, self.end_z        = orig_x, orig_y, orig_z, end_x, end_y, end_z        
            self.spacing_x, self.spacing_y, self.spacing_z                    = spacing_x, spacing_y, spacing_z                    
            self.nx, self.ny, self.nz                                         = nx, ny, nz                                         
            self.vp_min, self.vp_max, self.vs_min, self.vs_max = vp_min, vp_max, vs_min, vs_max
            self.xtomo,self.ytomo,self.ztomo,self.vp,self.vs=xtomo,ytomo,ztomo,vp,vs
        else:
            print 'no tomofile!!!!'
            #
            #
            #
    def tomo(self,x,y,z,vp_static=None,vs_static=None):
        if not vp_static:
            spac_x = (x - self.orig_x) / self.spacing_x
            spac_y = (y - self.orig_y) / self.spacing_y
            spac_z = (z - self.orig_z) / self.spacing_z
            #
            ix = int(spac_x)
            iy = int(spac_y)
            iz = int(spac_z)
            #
            gamma_interp_x = spac_x - float(ix)
            gamma_interp_y = spac_y - float(iy)
            #
            NX=self.nx
            NY=self.ny
            NZ=self.nz
            if(ix < 0):
                ix = 0
                gamma_interp_x = 0.
            if(ix > NX-2):         
                ix = NX-2          
                gamma_interp_x = 1.
            if(iy < 0):            
                iy = 0             
                gamma_interp_y = 0.
            if(iy > NY-2):         
                iy = NY-2          
                gamma_interp_y = 1.
            if(iz < 0):            
                 iz = 0
            if(iz > NZ-2):
                 iz = NZ-2
            #
            p0 = ix+iy*NX+iz*(NX*NY)
            p1 = (ix+1)+iy*NX+iz*(NX*NY)
            p2 = (ix+1)+(iy+1)*NX+iz*(NX*NY)
            p3 = ix+(iy+1)*NX+iz*(NX*NY)
            p4 = ix+iy*NX+(iz+1)*(NX*NY)
            p5 = (ix+1)+iy*NX+(iz+1)*(NX*NY)
            p6 = (ix+1)+(iy+1)*NX+(iz+1)*(NX*NY)
            p7 = ix+(iy+1)*NX+(iz+1)*(NX*NY)
            #
            if self.ztomo[p4]==self.ztomo[p0]:
                gamma_interp_z1 = 1
            else:
                gamma_interp_z1 = (z-self.ztomo[p0])/(self.ztomo[p4]-self.ztomo[p0])
            if(gamma_interp_z1 > 1.): gamma_interp_z1 = 1.
            if(gamma_interp_z1 < 0.): gamma_interp_z1 = 0.
            #
            if self.ztomo[p5]==self.ztomo[p1]:
                gamma_interp_z2 = 1
            else:
                gamma_interp_z2 = (z-self.ztomo[p1])/(self.ztomo[p5]-self.ztomo[p1])
            if(gamma_interp_z2 > 1.): gamma_interp_z2 = 1.
            if(gamma_interp_z2 < 0.): gamma_interp_z2 = 0.
            #
            if self.ztomo[p6]==self.ztomo[p2]:
                gamma_interp_z3 = 1
            else:
                gamma_interp_z3 = (z-self.ztomo[p2])/(self.ztomo[p6]-self.ztomo[p2])
            if(gamma_interp_z3 > 1.): gamma_interp_z3 = 1.
            if(gamma_interp_z3 < 0.): gamma_interp_z3 = 0.
            #
            if self.ztomo[p7]==self.ztomo[p3]:
                gamma_interp_z4 = 1
            else:
                gamma_interp_z4 = (z-self.ztomo[p3])/(self.ztomo[p7]-self.ztomo[p3])
            if(gamma_interp_z4 > 1.): gamma_interp_z4 = 1.
            if(gamma_interp_z4 < 0.): gamma_interp_z4 = 0.
            #
            gamma_interp_z5 = 1. - gamma_interp_z1
            gamma_interp_z6 = 1. - gamma_interp_z2
            gamma_interp_z7 = 1. - gamma_interp_z3
            gamma_interp_z8 = 1. - gamma_interp_z4
            #
            vp1 = self.vp[p0]
            vp2 = self.vp[p1]
            vp3 = self.vp[p2]
            vp4 = self.vp[p3]
            vp5 = self.vp[p4]
            vp6 = self.vp[p5]
            vp7 = self.vp[p6]
            vp8 = self.vp[p7]
            #       [  ]
            vs1 = self.vs[p0]
            vs2 = self.vs[p1]
            vs3 = self.vs[p2]
            vs4 = self.vs[p3]
            vs5 = self.vs[p4]
            vs6 = self.vs[p5]
            vs7 = self.vs[p6]
            vs8 = self.vs[p7]
            #
            vp_final = vp1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + \
               vp2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + \
               vp3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + \
               vp4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + \
               vp5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + \
               vp6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + \
               vp7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + \
               vp8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
            #
            vs_final = vs1*(1.-gamma_interp_x)*(1.-gamma_interp_y)*(1.-gamma_interp_z1) + \
               vs2*gamma_interp_x*(1.-gamma_interp_y)*(1.-gamma_interp_z2) + \
               vs3*gamma_interp_x*gamma_interp_y*(1.-gamma_interp_z3) + \
               vs4*(1.-gamma_interp_x)*gamma_interp_y*(1.-gamma_interp_z4) + \
               vs5*(1.-gamma_interp_x)*(1.-gamma_interp_y)*gamma_interp_z1 + \
               vs6*gamma_interp_x*(1.-gamma_interp_y)*gamma_interp_z2 + \
               vs7*gamma_interp_x*gamma_interp_y*gamma_interp_z3 + \
               vs8*(1.-gamma_interp_x)*gamma_interp_y*gamma_interp_z4
            #
            if(vp_final < self.vp_min): vp_final = self.vp_min
            if(vs_final < self.vs_min): vs_final = self.vs_min
            if(vp_final > self.vp_max): vp_final = self.vp_max
            if(vs_final > self.vs_max): vs_final = self.vs_max
            return vp_final,vs_final
        else:
            return vp_static,vs_static
    #
    def hex_simulation_parameter(self,h,vp_static=None,vs_static=None):
        nodes=cubit.get_connectivity('Hex',h)
        id_nodes=[0,1,2,3,4,5,6,7]
        id_faces=[[1,3,4],[0,2,5],[1,3,6],[0,2,7],[5,7],[4,6],[5,7],[4,6]]
        dt=[]
        pmax=[]
        vmax=[]
        vmin=[]
        for i in id_nodes[:-1]:
            for j in id_nodes[i+1:]:
                x1,y1,z1=cubit.get_nodal_coordinates(nodes[i])
                x2,y2,z2=cubit.get_nodal_coordinates(nodes[j])
                nvp1,nvs1=self.tomo(x1,y1,z1,vp_static,vs_static)
                nvp2,nvs2=self.tomo(x2,y2,z2,vp_static,vs_static)
                d=math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                #pmax_tmp=d*self.gllcoeff/min(nvp1,nvp2,nvs1,nvs2)*self.Ngll_per_wavelength
                pmax_tmp=d*(.5-self.gllcoeff)/min(nvp1,nvp2,nvs1,nvs2)*self.Ngll_per_wavelength #more conservative.....
                dt_tmp=self.Cmax*d*self.gllcoeff/max(nvp1,nvp2,nvs1,nvs2)
                dt.append(dt_tmp)
                pmax.append(pmax_tmp)
                vmax.append(max(nvp1,nvp2,nvs1,nvs2))
                vmin.append(min(nvp1,nvp2,nvs1,nvs2))
        return min(dt),max(pmax),min(vmin),max(vmax)
    #
    def group_period(self):
        tot=0.
        if self.period_hyst is not None:
            pmin=min(self.period_hyst.values())
            period=self.hyst(pmin,self.period,self.period_hyst)
            factor=(self.period-pmin)/self.nbin
            for i in range(0,self.nbin+1):
                if period.has_key(i):
                    txt='group "period_%.1e_%.1e" add hex ' %(pmin+factor*i,pmin+factor*(i+1))
                    txt=txt+' '.join(str(hh) for hh in period[i])
                    cubit.cmd(txt)
    def group_timestep(self):
         tot=0.
         if self.dt_hyst is not None:
             dtmax=max(self.dt_hyst.values())
             dt=self.hyst(self.dt,dtmax,self.dt_hyst)
             factor=(dtmax-self.dt)/self.nbin
             for i in range(0,self.nbin+1):
                 if dt.has_key(i):
                     txt='group "timestep_%.1e_%.1e" add hex ' %(self.dt+factor*i,self.dt+factor*(i+1))
                     txt=txt+' '.join(str(hh) for hh in dt[i])
                     cubit.cmd(txt)


#cubit.cmd('brick x 10000')
#cubit.cmd('mesh vol 1')
#cubit.cmd('refine hex in node in surf 1')
#cubit.cmd('refine hex in node in surf 3')
#vp=1000
#vs=600
#mesh=SEM_stability_3D()
#mesh.check_simulation_parameter(vp_static=vp,vs_static=vs)
#mesh.group_timestep()
#mesh.group_period()