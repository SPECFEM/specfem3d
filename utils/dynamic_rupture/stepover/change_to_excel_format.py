from __future__ import print_function

import sys
import os
import numpy as np
import pandas as pd

def convert_to_excel():
    data = np.genfromtxt('./logfile',dtype=None)
    f = open('excelinput.txt','w')
    for x in data:
        depth = 0.0005*(1e5-x[0])
        stress = x[1]
        Dmax = x[2]
        Vmax = x[3]
        if x[4]:
            PF = 'S'
        else:
            PF = 'F'

        f.write('%f, , ,%f,%f,%s,%f\n'%(depth,stress,Vmax,PF,Dmax))
    f.close()

convert_to_excel()

