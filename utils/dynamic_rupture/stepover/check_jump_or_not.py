from __future__ import print_function

import numpy as np
import sys

def checkjump(time):
    rvalue = False
    prefix = '../OUTPUT_FILES/'
    with open(prefix+'Snapshot0_F2.bin','r') as f1:
        for i in range(1,9):
            header = np.fromfile(f1,dtype=np.int32,count=1)
            length = header/4
            data = np.fromfile(f1,dtype=np.float32,count=length)
            tail = np.fromfile(f1,dtype=np.int32,count=1)
        if(i==8):
            print("shear stress",max(data))

    with open(prefix+'Snapshot%s_F1.bin'%str(time),'r') as f:
        for i in range(1,9):
            header = np.fromfile(f,dtype=np.int32,count=1)
            length = header/4
            data = np.fromfile(f,dtype=np.float32,count=length)
            tail = np.fromfile(f,dtype=np.int32,count=1)
            if(i == 2):
                print("stepover distance:",max(data))
            if( i == 4):
                print("max displacement fault 1",max(abs(data)))
                Dmax1 = max(abs(data))
            if( i == 6):
                print("slip veloc on 1",max(abs(data)))

    with open(prefix+'Snapshot%s_F2.bin'%str(time),'r') as f:
        for i in range(1,9):
            header = np.fromfile(f,dtype=np.int32,count=1)
            length = header/4
            data = np.fromfile(f,dtype=np.float32,count=length)
            tail = np.fromfile(f,dtype=np.int32,count=1)
            if(i == 1):
                print("x range",max(data),'-',min(data))
            if(i == 2):
                print("stepover distance:",max(data))
            if( i == 4 or i == 5):
                print("max displacement",max(abs(data)))
                if(i == 4):
                    Dmax = max(abs(data))
                    if(max(abs(data))>2.0):
                        rvalue = True
            if( i == 6):
                Vmax = max(abs(data))
                if(max(abs(data)) > 0.001):
                    print("jump happens!",max(abs(data)))
                else:
                    print("jump fails!",max(abs(data)))
#this part for checking if the rupture on the first fault goes supershear!
    with open(prefix+'Snapshot10000_F1.bin','r') as f:
         for i in range(1,9):
            header = np.fromfile(f,dtype=np.int32,count=1)
            length = header/4
            data = np.fromfile(f,dtype=np.float32,count=length)
            tail = np.fromfile(f,dtype=np.int32,count=1)
            if(i == 1):
                A = (data > 15000.0)
            if( i== 4):
                Dmax = max(abs(data[A]))
                if(Dmax>0.1):
                    supershear = True
                else:
                    supershear = False


    return rvalue,Dmax,Vmax,Dmax1,supershear
def main():
    print(checkjump(20000))
if __name__ == "__main__":
    main()

