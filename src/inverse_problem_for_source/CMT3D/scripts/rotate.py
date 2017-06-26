from numpy import *
from math import *
from seismo_const import *


def rotate_matrix(lon,lat):
  R=zeros((3,3))
  th=(90-lat)*PI/180; phi=lon*PI/180
  R[0,0]=sin(th)*cos(phi)
  R[0,1]=cos(th)*cos(phi)
  R[0,2]=-sin(phi)
  R[1,0]=sin(th)*sin(phi)
  R[1,1]=cos(th)*sin(phi)
  R[1,2]=cos(phi)
  R[2,0]=cos(th)
  R[2,1]=-sin(th)
  R[2,2]=0

  return R

def rotate_cmt(lon,lat,mrr,mtt,mpp,mrt,mrp,mtp,local2global=True):

  Min=zeros((3,3)); Mout=zeros((3,3)); R=zeros((3,3))
  R=rotate_matrix(lon,lat)
  Min[0,0]=mrr;      Min[0,1]=mrt;      Min[0,2]=mrp
  Min[1,0]=Min[0,1]; Min[1,1]=mtt;      Min[1,2]=mtp
  Min[2,0]=Min[0,2]; Min[2,1]=Min[1,2]; Min[2,2]=mpp
  if local2global:
    Mout=dot(R,  dot(Min,R.transpose()))
  else:
    Mout=dot(R.transpose(), dot(Min, R))
  mrr1=Mout[0,0];  mtt1=Mout[1,1];   mpp1=Mout[2,2]
  mrt1=Mout[0,1];  mrp1=Mout[0,2];   mtp1=Mout[1,2]

  return mrr1,mtt1,mpp1,mrt1,mrp1,mtp1


def rotate_loc(lon,lat,locx,locy,locz,local2global=True):
  lin=array([locx,locy,locz]); lout=zeros((3)); R=zeros((3,3))
  R=rotate_matrix(lon,lat)
  if local2global:
    lout=dot(R, lin)
  else:
    lout=dot(R.transpose(),lin)

  return lout[0],lout[1],lout[2]

def sph2xyz(lon,lat):
  th=(90-lat)*PI/180;
  phi=lon*PI/180;
  X=sin(th)*cos(phi)
  Y=sin(th)*sin(phi)
  Z=cos(th)
  return X,Y,Z

def xyz2sph(x,y,z):
  r=sqrt(x**2+y**2+z**2);
  th=acos(z/r);
  if abs(th)<1.e-5 or abs(th-PI)<1.e-5:
    phi=0
  else:
    phi=atan2(y/sin(th),x/sin(th))
  lon=phi*180/PI; lat=90-th*180/PI
  return lon, lat, r

