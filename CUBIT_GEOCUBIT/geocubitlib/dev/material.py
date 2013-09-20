#############################################################################
# material.py                                                    
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



class material():
    def __init__(self,name,vp,vs,rho,q,ani):
        self.name=name
        self.vp=vp
        self.vs=vs
        self.rho=rho
        self.qp=q
        self.qs=ani
    def __repr__(self):
        print self.name
        print 'vp',self.vp
        print 'vs',self.vs
        print 'rho',self.rho
        print 'Qp',self.qp
        print 'Qs',self.qs
        return
    def __str__(self):
        txt=self.name+' vp '+str(self.vp)+ ' vs '+str(self.vs)+ ' rho '+str(self.rho)+ ' Qp '+str(self.qp)+ ' Qs '+str(self.qs)
        return txt


def read_material_archive(archivefile='material_archive.dat'):
    archive=open(archivefile,'r')
    archive_dict={}
    for record in archive:
        m=material(*record.split())
        archive_dict[m.name]=m
    return archive_dict


def block_description(archivefile='material_archive.dat'):
    try:
        import initializing as initializing
        cubit   = initializing.initializing_cubit()
    except:
        pass
    import menu as menu
    archive=read_material_archive(archivefile=menu.material_file)
    b_assign=menu.material_assignement
    print b_assign
    if len(b_assign) != 0:
        for block in b_assign:
            mat=archive[block[1]]
            print mat
            cubit.cmd("block "+block[0]+" attribute count 6")
            cubit.cmd("block "+block[0]+" name             '"+ mat.name   +'"')
            cubit.cmd("block "+block[0]+" attribute index 1 "+ '1'            )
            cubit.cmd("block "+block[0]+" attribute index 2 "+ str(mat.vp)    )
            cubit.cmd("block "+block[0]+" attribute index 3 "+ str(mat.vs)    )
            cubit.cmd("block "+block[0]+" attribute index 4 "+ str(mat.rho)   )
            cubit.cmd("block "+block[0]+" attribute index 5 "+ str(mat.qp)    )
            cubit.cmd("block "+block[0]+" attribute index 5 "+ str(mat.qs)    )
    
    