#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:01:31 2017

@author: cheetham
"""

import glob,shutil
import numpy as np

data_dir = '/data/NACO/Science/'
save_dir = '/data/NACO/contrasts/'

# Look for the data in here
contrast_files = glob.glob(data_dir+'*/*/ADI/GRAPHIC_PCA/contrast.txt')

# Loop through and copy them
for contrast_file in contrast_files:
    
    targ_name = contrast_file.split('/')[-4]
    date = contrast_file.split('/')[-5]
    
    date = date.replace('-','')
    
    new_name = save_dir+targ_name+'_'+date+'_contrast.txt'
    
    # Copy the file
    # print('Copying:'+contrast_file)
    # shutil.copy2(contrast_file,new_name)

    # Load the file and write it out again
    con_data = np.loadtxt(contrast_file)
    hdr = "Separation (arcsec) Contrast (mag)  using plate scale of 0.027"

    # if con_data.shape[0] < con_data.shape[1] :
	con_data = con_data.T

    np.savetxt(new_name,con_data,header=hdr)

