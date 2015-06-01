#!/usr/bin/env python
# File: mk_regions.py
# Created on: Sun May 27 21:56:23 2012
# Last Change: Tue May 29 16:13:20 2012
# Purpose of script: <+INSERT+>
# Author: Steven Boada

import sys
import os
import numpy as np

#def mk_regions(coords):
def mk_regions():

    for coords in os.listdir('.'):
        if coords.endswith('.txt'):
            d = np.loadtxt(coords,dtype=str)
            f1=open(coords.rstrip('txt')+'reg','wt')
            f1.writelines('# Region file formart: DS9 version 4.1\n')
            f1.writelines('# Filename: '+coords.rstrip('txt')+'reg\n')
    #f1.writelines('global color=green\n')
    
            for i in range(len(d)):
                f1.writelines('fk5;circle('+d[i][1]+',')
                #f1.writelines(d[i][2]+'0.00055)')
                f1.writelines(d[i][2]+',2")')
                f1.writelines('# width=2 text={'+d[i][0]+'} ')
                f1.writelines('tag={'+coords.rstrip('_coords.txt')+'}\n')

            f1.close()
if __name__=='__main__':
    #mk_regions(sys.argv[1])
    mk_regions()
