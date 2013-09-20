#############################################################################
# quality_log.py                                                    
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
try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print 'error importing cubit, check if cubit is installed'
        pass


def quality_log(tqfile=None):
    """
    creation of the quality parameter file
    """
    import start as start
    #
    #
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    #
    from hex_metric import SEM_metric_3D
    #
    lvol=cubit.parse_cubit_list('volume','all')
    if len(lvol)!=0:
        cubit.cmd('quality vol all allmetric                                      ')
    else:
        cubit.cmd('quality hex in block all allmetric                                      ')
    cubit.cmd('list model                                                      ')
    #
    toclose=True
    if isinstance(tqfile,file): 
        totstat_file=tqfile
    elif isinstance(tqfile,str):
        totstat_file=open(tqfile+'_cubitquality_skewness_proc_'+str(iproc)+'.log','w')
    else:
        import sys
        totstat_file=sys.stdout
        toclose=False
    
    mesh=SEM_metric_3D()
    mesh.check_metric()
    
    if mesh.max_skewness is not None:
        mesh.skew_hystogram=mesh.hyst(0,mesh.max_skewness,mesh.skew_hyst)
        totstat_file.write('-'*70+'\n')
        
        
        
        
        
        totstat_file.write('='*70+'\n')
        totstat_file.write('SKEWNESS'+'\n')
        totstat_file.write('='*70+'\n') 
        if len(mesh.hex_max_skewness) <= 30:
            totstat_file.write('max = '+str(mesh.max_skewness)+' in hexes '+str(mesh.hex_max_skewness)+'\n')
            totstat_file.write('(angle -> minimun ='+str(mesh.min_angle)+ ' maximun ='+str(mesh.max_angle)+')'+'\n')
        else:
            totstat_file.write('max = '+str(mesh.max_skewness)+' in '+str(len(mesh.hex_max_skewness))+' hexes '+'\n')
            totstat_file.write('(angle -> minimun ='+str(mesh.min_angle)+' maximun ='+str(mesh.max_angle)+')'+'\n')
        totstat_file.write('-'*70+'\n')
        totstat_file.write('skew hystogram')
        totstat_file.write('-'*70+'\n')
        tot=0
        for i in mesh.skew_hystogram.values():
            tot=tot+len(i)
        #k=mesh.skew_hystogram.keys()
        #k.sort()
        factor=mesh.max_skewness/mesh.nbin
        for i in range(0,mesh.nbin+1):
            if mesh.skew_hystogram.has_key(i):
                if (i+1)*factor <= 1:
                    totstat_file.write(str(i)+' ['+str(i*factor)+'->'+str((i+1)*factor)+'[ : '+str(len(mesh.skew_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.skew_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                if (i+1)*factor <= 1:
                    totstat_file.write(str(i)+' ['+str(i*factor)+'->'+str((i+1)*factor)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
        totstat_file.write('-'*70+'\n')
    ###############################################
    if mesh.min_edge_length is not None:
        mesh.edgemin_hystogram=mesh.hyst(mesh.min_edge_length,mesh.max_edge_length,mesh.edgemin_hyst)
        mesh.edgemax_hystogram=mesh.hyst(mesh.min_edge_length,mesh.max_edge_length,mesh.edgemax_hyst)
        totstat_file.write('='*70+'\n')
        totstat_file.write('edge length')
        totstat_file.write('='*70+'\n')
        if len(mesh.hex_min_edge_length) <= 30:
            totstat_file.write('minimum edge length: '+str(mesh.min_edge_length)+ ' in hexes '+str(mesh.hex_min_edge_length)+'\n')
        else:
            totstat_file.write('minimum edge length: '+str(mesh.min_edge_length)+ ' in '+str(len(mesh.hex_min_edge_length))+ ' hexes.'+'\n')
        if len(mesh.hex_max_edge_length) <= 30:
            totstat_file.write('maximum edge length: '+str(mesh.max_edge_length)+ ' in hexes '+str(mesh.hex_max_edge_length)+'\n')              
        else:                                                                                                                        
            totstat_file.write('maximum edge length: '+str(mesh.max_edge_length)+' in '+str(len(mesh.hex_max_edge_length))+ ' hexes.'+'\n')    
        totstat_file.write('-'*70+'\n')
        totstat_file.write('edge length hystogram')
        totstat_file.write('-'*70+'\n')
        factor=(mesh.max_edge_length-mesh.min_edge_length)/mesh.nbin
        totstat_file.write('minimum edge length'+'\n')
        tot=0
        for i in mesh.edgemin_hystogram.values():
            tot=tot+len(i)
        #k=mesh.edgemin_hystogram.keys()
        #k.sort()
        for i in range(0,mesh.nbin+1):
            if mesh.edgemin_hystogram.has_key(i):
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : '+str(len(mesh.edgemin_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.edgemin_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
        totstat_file.write('-'*70+'\n')
        totstat_file.write('maximum edge length')
        tot=0
        for i in mesh.edgemax_hystogram.values():
            tot=tot+len(i)
        #k=mesh.edgemax_hystogram.keys()
        #k.sort()
        for i in range(0,mesh.nbin+1):
            if mesh.edgemax_hystogram.has_key(i):
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : '+str(len(mesh.edgemax_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.edgemax_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
    try:
        if mesh.dt is not None:
            totstat_file.write('='*70+'\n')
            totstat_file.write('STABILITY')
            totstat_file.write('='*70+'\n')
            totstat_file.write('time step < '+str(mesh.dt)+'s, for velocity = '+str(mesh.velocity)+'\n')
    except:
        pass
    
    if toclose: 
        totstat_file.close()
    
    print 'max specfem3d skewness: ',mesh.max_skewness
    print 'min edge length: ',mesh.min_edge_length
    return mesh.max_skewness,mesh.min_edge_length
    
    
def quality_log(tqfile=None):
    """
    creation of the quality parameter file
    """
    import start as start
    #
    #
    mpiflag,iproc,numproc,mpi   = start.start_mpi()
    #
    #
    from hex_metric import SEM_metric_3D
    #
    lvol=cubit.parse_cubit_list('volume','all')
    if len(lvol)!=0:
        cubit.cmd('quality vol all allmetric                                      ')
    else:
        cubit.cmd('quality hex in block all allmetric                                      ')
    cubit.cmd('list model                                                      ')
    #
    toclose=True
    if isinstance(tqfile,file): 
        totstat_file=tqfile
    elif isinstance(tqfile,str):
        totstat_file=open(tqfile+'_cubitquality_skewness_proc_'+str(iproc)+'.log','w')
    else:
        import sys
        totstat_file=sys.stdout
        toclose=False

    mesh=SEM_metric_3D()
    mesh.check_metric()

    if mesh.max_skewness is not None:
        mesh.skew_hystogram=mesh.hyst(0,mesh.max_skewness,mesh.skew_hyst)
        totstat_file.write('-'*70+'\n')





        totstat_file.write('='*70+'\n')
        totstat_file.write('SKEWNESS'+'\n')
        totstat_file.write('='*70+'\n') 
        if len(mesh.hex_max_skewness) <= 30:
            totstat_file.write('max = '+str(mesh.max_skewness)+' in hexes '+str(mesh.hex_max_skewness)+'\n')
            totstat_file.write('(angle -> minimun ='+str(mesh.min_angle)+ ' maximun ='+str(mesh.max_angle)+')'+'\n')
        else:
            totstat_file.write('max = '+str(mesh.max_skewness)+' in '+str(len(mesh.hex_max_skewness))+' hexes '+'\n')
            totstat_file.write('(angle -> minimun ='+str(mesh.min_angle)+' maximun ='+str(mesh.max_angle)+')'+'\n')
        totstat_file.write('-'*70+'\n')
        totstat_file.write('skew hystogram')
        totstat_file.write('-'*70+'\n')
        tot=0
        for i in mesh.skew_hystogram.values():
            tot=tot+len(i)
        #k=mesh.skew_hystogram.keys()
        #k.sort()
        factor=mesh.max_skewness/mesh.nbin
        for i in range(0,mesh.nbin+1):
            if mesh.skew_hystogram.has_key(i):
                if (i+1)*factor <= 1:
                    totstat_file.write(str(i)+' ['+str(i*factor)+'->'+str((i+1)*factor)+'[ : '+str(len(mesh.skew_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.skew_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                if (i+1)*factor <= 1:
                    totstat_file.write(str(i)+' ['+str(i*factor)+'->'+str((i+1)*factor)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
        totstat_file.write('-'*70+'\n')
    ###############################################
    if mesh.min_edge_length is not None:
        mesh.edgemin_hystogram=mesh.hyst(mesh.min_edge_length,mesh.max_edge_length,mesh.edgemin_hyst)
        mesh.edgemax_hystogram=mesh.hyst(mesh.min_edge_length,mesh.max_edge_length,mesh.edgemax_hyst)
        totstat_file.write('='*70+'\n')
        totstat_file.write('edge length')
        totstat_file.write('='*70+'\n')
        if len(mesh.hex_min_edge_length) <= 30:
            totstat_file.write('minimum edge length: '+str(mesh.min_edge_length)+ ' in hexes '+str(mesh.hex_min_edge_length)+'\n')
        else:
            totstat_file.write('minimum edge length: '+str(mesh.min_edge_length)+ ' in '+str(len(mesh.hex_min_edge_length))+ ' hexes.'+'\n')
        if len(mesh.hex_max_edge_length) <= 30:
            totstat_file.write('maximum edge length: '+str(mesh.max_edge_length)+ ' in hexes '+str(mesh.hex_max_edge_length)+'\n')              
        else:                                                                                                                        
            totstat_file.write('maximum edge length: '+str(mesh.max_edge_length)+' in '+str(len(mesh.hex_max_edge_length))+ ' hexes.'+'\n')    
        totstat_file.write('-'*70+'\n')
        totstat_file.write('edge length hystogram')
        totstat_file.write('-'*70+'\n')
        factor=(mesh.max_edge_length-mesh.min_edge_length)/mesh.nbin
        totstat_file.write('minimum edge length'+'\n')
        tot=0
        for i in mesh.edgemin_hystogram.values():
            tot=tot+len(i)
        #k=mesh.edgemin_hystogram.keys()
        #k.sort()
        for i in range(0,mesh.nbin+1):
            if mesh.edgemin_hystogram.has_key(i):
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : '+str(len(mesh.edgemin_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.edgemin_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
        totstat_file.write('-'*70+'\n')
        totstat_file.write('maximum edge length')
        tot=0
        for i in mesh.edgemax_hystogram.values():
            tot=tot+len(i)
        #k=mesh.edgemax_hystogram.keys()
        #k.sort()
        for i in range(0,mesh.nbin+1):
            if mesh.edgemax_hystogram.has_key(i):
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : '+str(len(mesh.edgemax_hystogram[i]))+'/'+str(tot)+' hexes ('+str(len(mesh.edgemax_hystogram[i])/float(tot)*100.)+'%)'+'\n')
            else:
                totstat_file.write(str(i)+' ['+str(i*factor+mesh.min_edge_length)+'->'+str((i+1)*factor+mesh.min_edge_length)+'[ : 0/'+str(tot)+' hexes (0%)'+'\n')
    try:
        if mesh.dt is not None:
            totstat_file.write('='*70+'\n')
            totstat_file.write('STABILITY')
            totstat_file.write('='*70+'\n')
            totstat_file.write('time step < '+str(mesh.dt)+'s, for velocity = '+str(mesh.velocity)+'\n')
    except:
        pass

    if toclose: 
        totstat_file.close()

    print 'max specfem3d skewness: ',mesh.max_skewness
    print 'min edge length: ',mesh.min_edge_length
    return mesh.max_skewness,mesh.min_edge_length