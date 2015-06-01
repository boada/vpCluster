#!/usr/bin/env python
# File: mk_lists.py
# Created on: Wed 16 May 2012 01:38:20 PM CDT
# Last Change: Sat May 26 14:39:24 2012
# Purpose of script: <+INSERT+>
# Author: Steven Boada

import pyfits as pyf
import numpy as np
import os
import sys

def mk_lists():
    # Read these from a param file
    std_star = 'G191B2B'
    std_star2 = 'HD84937' 
    exp = '1200'
    Ist_science = 'vp0078'


    # Open the required lists for building
    f1 = open('bias.list','wt')
    f2 = open('flat_n1.list','wt')
    f3 = open('arc_n1.list','wt')
    f4 = open('science_n1.list','wt')
    #f5 = open('sky_n1.list','wt') -- DO NOT NEED

    #Run through all of the fits images
    #os.chdir("./")
    a = os.listdir(".")
    a=np.sort(a)
    for index, files in enumerate(a):
        if files.endswith(".fits"):
            hdulist=pyf.open(files)
            obj = hdulist[0].header['object']
            typ = hdulist[0].header['imagetyp']
            print files,typ
            files = files.rstrip('.fits')
            # Add the image to the right list
            if typ =='object': # science
                if std_star in obj or std_star2 in obj:
                    print 'star'
                    f4.writelines(files+' '+files+' '+files+'\n')
                else:
                    f4.writelines(files+' '+files+' '+files+'\n')
                    if files == Ist_science:
                        f5.writelines(files+' '+files+' '+files+' '+exp+' '+\
                                exp+'\n')
                    # Commented the Sky part
                    #else:
                    #    f5.writelines(files+' '+a[index-1].rstrip('.fits')+' '+\
                    #            a[index+1].rstrip('.fits')+' '+exp+' '+exp+'\n')

            elif typ =='comp': # arc
                f3.writelines(files+'\n')
            elif typ =='zero': #bias
                f1.writelines(files+'\n')
            elif typ =='flat': #flat
                f2.writelines(files+'\n')

        else:
            pass

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    #f5.close()

if __name__ =="__main__":
    mk_lists()

