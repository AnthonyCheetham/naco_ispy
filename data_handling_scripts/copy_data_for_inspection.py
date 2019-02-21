#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:03:27 2018

@author: cheetham
"""

import numpy as np
import glob
import os
import argparse
import shutil


parser = argparse.ArgumentParser(description='Copy relevant files to one location for easy inspection of results')

parser.add_argument('--data_dir', action="store", dest="data_dir",default="/data/NACO/Science/", 
                    help='Location of the main NACO data directory to search for files')
parser.add_argument('--date_pattern', action="store", dest="date_pattern",  default="201*", 
                    help='Find files matching this date')
parser.add_argument('--target_pattern', action="store", dest="target_pattern",  default="*", 
                    help='Find files matching this target name')
parser.add_argument('--output_dir', action="store", dest="output_dir",  default="/home/spectro/cheetham/tmp/inspection/", 
                    help='The directory to save the files')

args = parser.parse_args()

data_dir = args.data_dir
date_pattern = args.date_pattern
target_pattern = args.target_pattern
output_dir = args.output_dir

# Find the relevant directories
dirs = glob.glob(data_dir+date_pattern+'/'+target_pattern+'/ADI/')
print('Found '+str(len(dirs))+' matching datasets')

# Loop through them
for d in dirs:
    
    # Get a directory ready in the output
    targname = d.split('/')[-3]
    dir_out = output_dir+targname+'/'
    dir_exists=os.access(dir_out, os.F_OK)    
    if not dir_exists:
        os.makedirs(dir_out)

    
    # Copy the files
    for f in ['master_cube_PCA.fits','flux.fits','GRAPHIC_PCA/pca_multimodes.fits',
              'GRAPHIC_PCA/snr_map.fits','GRAPHIC_PCA/contrast.txt','parallactic_angle.txt']:
        
        fname = d+f
        file_exists=os.access(fname, os.F_OK)    
        if file_exists:
            fname_out = dir_out+f.replace('GRAPHIC_PCA/','')
            shutil.copy2(fname,fname_out)
        else:
            print('Didnt find '+f+' for target '+targname)
            print('  '+fname)
            