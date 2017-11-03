# -*-coding:Latin-1 -*
from math import *
from numpy import *



# ################################# functions ################################
# ################################# functions ################################
# ################################# functions ################################
# ################################# functions ################################

#----------------------------------------------------------------------------
# compute 3D rotation matrix with respect to axis for rotation_angle in degree
#
def rotation_with_respect_to_axis(axis, rotation_angle):
    pi=3.1415926535897932
    deg2rad = pi / 180.
    # on normalise l'axe
    norme_axis=sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)

    # normalization
    ux=axis[0]/norme_axis
    uy=axis[1]/norme_axis
    uz=axis[2]/norme_axis

    # compute cos and sin
    c=cos(deg2rad * rotation_angle)
    s=sin(deg2rad * rotation_angle)
    R=zeros((3,3))

    # rotation matrix
    R[0,0]=(ux**2 + (1.-ux**2)*c)
    R[0,1]=(ux*uy*(1.-c)-uz*s)
    R[0,2]=(ux*uz*(1.-c)+uy*s)

    R[1,0]=(ux*uy*(1.-c)+uz*s)
    R[1,1]=(uy**2+(1.-uy**2)*c)
    R[1,2]=(uy*uz*(1.-c)-ux*s)

    R[2,0]=(ux*uz*(1.-c)-uy*s)
    R[2,1]=(uy*uz*(1.-c)+ux*s)
    R[2,2]=(uz**2+(1.-uz**2)*c)

    return R
#
#----------------------------------------------------------------------------
# compute 3D rotation matrix to change coordiante spheric - cartesian
# in order to switch specfem local coordiante - axisem global coordiantes
#
def general_matrix_rotation(lon, lat, azi):
    # angles in radian
    cla=(90.-lat) * 3.1415926535897932 / 180.
    lo=lon * 3.1415926535897932 / 180.

    # compute rotation matrix without azimuth rotation of chunk
    R=zeros((3,3))

    R[0,0]=-sin(lo)
    R[1,0]=cos(lo)
    R[2,0]=0.

    R[0,1]=-cos(cla)*cos(lo)
    R[1,1]=-cos(cla)*sin(lo)
    R[2,1]=sin(cla)

    R[0,2]=sin(cla)*cos(lo)
    R[1,2]=sin(cla)*sin(lo)
    R[2,2]=cos(cla)

    # add azimuth chunk rotation contribution
    axis=zeros(3)
    axis=[cos(lo)*sin(cla), sin(lo)*sin(cla), cos(cla)]
    R1=rotation_with_respect_to_axis(axis, -azi)
    rot=zeros((3,3))
    # compute rot=R1*R
    i=0
    j=0
    while j<3:
        i=0
        while i<3:
            k=0
            rot[i,j]=0.
            while k<3:
                rot[i,j] = rot[i,j] + R1[i,k]*R[k,j]
                k+=1
            i+=1
        j+=1

    return rot
#----------------------------------------------------------------------
def Cartesian2spheric(Rotation, point):
    deg2rad= pi / 180.
    Spoint=zeros(3)
    j=0
    while j<3:
        i=0
        while i<3:
            Spoint[j]=Spoint[j]+Rotation[j,i]*point[i]
            i+=1
        j+=1
    # compute geographic coordinate
    x=Spoint[0]
    y=Spoint[1]
    z=Spoint[2]
    radius=sqrt(Spoint[0]**2 + Spoint[1]**2 + Spoint[2]**2)
    longitude=atan2(y,x)
    latitude=asin(z/radius)
    Spoint[0]=longitude/deg2rad
    Spoint[1]=latitude/deg2rad
    Spoint[2]=radius/1000.

    return Spoint

def Cartesian2spheric1(Rotation, point):
    deg2rad= pi / 180.
    Spoint=zeros(3)
    j=0
    while j<3:
        i=0
        while i<3:
            Spoint[j]=Spoint[j]+Rotation[j,i]*point[i]
            i+=1
        j+=1
    # compute geographic coordinate
    #x=Spoint[0]
    #y=Spoint[1]
    #z=Spoint[2]
    #radius=sqrt(Spoint[0]**2 + Spoint[1]**2 + Spoint[2]**2)
    #longitude=atan2(y,x)
    #latitude=asin(z/radius)
    #Spoint[0]=longitude/deg2rad
    #Spoint[1]=latitude/deg2rad
    #Spoint[2]=radius/1000.

    return Spoint

########################################################## MAIN #############################################################
########################################################## MAIN #############################################################
########################################################## MAIN #############################################################
########################################################## MAIN #############################################################
########################################################## MAIN #############################################################
########################################################## MAIN #############################################################

# constants
RADIUS_EARTH=6371000.

# 1/ --- READING CHUNK PARAMETER --------
# hardcoded file that contains position of chunk
chunk_parameter_file="MESH/ParFileMeshChunk"

fid=open(chunk_parameter_file,"r")
list_of_line=fid.readlines()
fid.close()

print "\n reading chunk parameter in "+ chunk_parameter_file+ " file \n"

# get third line that define chunk position
line=list_of_line[3].split()
chunk_longitue_in_degree=float(line[0])
chunk_latitude_in_degree=float(line[1])
chunk_azimuth_in_degree=float(line[2])


# 2/ ------ COMPUTING ROTATION MATRIX -------------
print"\n computing rotation matrix for chunk :"
print "    long :"+str(chunk_longitue_in_degree),"  lat :"+str(chunk_latitude_in_degree)+"   azimuth :"+str(chunk_azimuth_in_degree)+"\n\n"
Rotation_matrix=general_matrix_rotation(chunk_longitue_in_degree, chunk_latitude_in_degree, chunk_azimuth_in_degree)
print "\n rotation matrix : \n"
print Rotation_matrix



# 3/ ------ CARTESIAN TO SPHERIC TRANSFORMATION ---------
# hardcoded file contains gll points
list_file="normals.txt"
fid=open(list_file,"r")
list_of_ggl_point=fid.readlines()
fid.close()
# open file to write gll points in spherical coordinates
fid=open("normals.sph",'w')
gll_point=zeros(3)
gll_normal=zeros(3)
for line in list_of_ggl_point:
    one_line_list=line.split()
    gll_point[0]=float(one_line_list[0])
    gll_point[1]=float(one_line_list[1])
    gll_point[2]=float(one_line_list[2]) + RADIUS_EARTH
    gll_spheric=Cartesian2spheric(Rotation_matrix, gll_point)
    gll_normal[0]=float(one_line_list[3])
    gll_normal[1]=float(one_line_list[4])
    gll_normal[2]=float(one_line_list[5])
    gll_normal_spheric=Cartesian2spheric(Rotation_matrix, gll_normal)
    fid.write("{0}   {1}   {2}   {3}   {4}   {5}\n".format(gll_spheric[0],gll_spheric[1],gll_spheric[2], gll_normal_spheric[0], gll_normal_spheric[1], gll_normal_spheric[2]))
    #fid.write("{0}   {1}   {2}   {3}   {4}   {5}\n".format(gll_point[0],gll_point[1],gll_point[2], gll_normal[0], gll_normal[1], gll_normal[2]))
fid.close()
