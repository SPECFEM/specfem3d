from __future__ import print_function

import sys
import check_jump_or_not
import re
import os
import subprocess
def setparameter(stress, depth):
    with open('./Par_file_faults','rw') as f1:
        content1 = f1.readlines()
    f2 =  open('./Par_file_faults_new','w')
    s1 = "TAG1"
    s2 = "TAG2"
    s3 = "TAG3"
    for x in content1:
        newx = x
        if( s1 in x):
            newx = ""
            str_ = x.split(',')
            str_[7] = 'lz='+str(depth)
            newx = ','.join(str_)
            newx = newx + '/ #TAG1\n'
        if( s2 in x):
            newx = ""
            str_ = x.split(',')
            str_[3] = str(stress)
            newx = ','.join(str_)
        if( s3 in x):
            newx = ""
            str_ = x.split(',')
            str_[7] = 'lz='+str(depth)
            newx = ','.join(str_)
            newx = newx + '/ #TAG3\n'
        f2.write(newx)
    f1.close()
    f2.close()
    p = subprocess.Popen("mv Par_file_faults_new Par_file_faults", stdout = subprocess.PIPE, shell = True)
    p.communicate()


def run_mpi_job():
    p = subprocess.Popen("mpirun -np 90 -perhost 3 -f ~/cpunodes ../bin/xspecfem3D",stdout = subprocess.PIPE,shell = True)
    print("mpi job running")
    p.communicate()

def bisect_search(Tmax,Tmin,h):
    if(Tmax < Tmin):
        return bisect_search(Tmin,Tmax,h)
    if(abs(Tmax - Tmin)<1e5):
        return (Tmax,Tmin)
    Tmid = 0.5*(Tmax + Tmin)
    print('now testing T = '+str(Tmid))
    pf = test(Tmid,h)
    print('test results:'+ str(pf))
    if(pf):
        Tmax = Tmid
    else:
        Tmin = Tmid
    (a,b) = bisect_search(Tmax,Tmin,h)
    return (a,b)

def test(T,h):
    setparameter(T,h)
    run_mpi_job()
    f = open('logfile','a')
    pf,Vmax,Dmax,D1max,supershear = check_jump_or_not.checkjump(80000)
    f.write(str(h)+'\t'+str(T)+'\t'+str(Vmax)+'\t'+str(Dmax)+'\t'+str(pf)+'\t'+str(D1max)+'\t'+str(supershear)+'\n')
    return pf



def main():

    (a,b) = bisect_search(68e6,74e6,90e3)
    (a,b) = bisect_search(67e6,73e6,85e3)
    (a,b) = bisect_search(66e6,72e6,80e3)
    (a,b) = bisect_search(65e6,71e6,75e3)
    (a,b) = bisect_search(64e6,70e6,70e3)
    (a,b) = bisect_search(63e6,69e6,65e3)
    (a,b) = bisect_search(62e6,68e6,60e3)

    print(a,b)


if __name__ == "__main__":
    main()

