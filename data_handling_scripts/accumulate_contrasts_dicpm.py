#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:01:31 2017

@author: cheetham
"""

import glob,shutil
import astropy.io.fits as pf
import numpy as np

data_dir = '/data/NACO/Science/'
save_dir = '/data/NACO/dicpm_contrasts/'

# Look for the data in here
contrast_files = glob.glob(data_dir+'*/*/ADI/DICPM/contrast_table*.fits')

# Loop through and copy them
for contrast_file in contrast_files:
    
    targ_name = contrast_file.split('/')[-4]
    date = contrast_file.split('/')[-5]
    
    date = date.replace('-','')
    
    new_name = save_dir+targ_name+'_'+date+'_contrast.txt'

    # Load the file and write it out again

    all_con_data = pf.getdata(contrast_file)
    seps = all_con_data['sep (pix)'] * 0.027
    cons = -2.5*np.log10(all_con_data['contrast_50'])
    hdr = "Separation (arcsec) Contrast (mag)  using plate scale of 0.027"
    con_data = np.array([seps,cons])

    # Remove all the NaNs
    good_ix = np.isfinite(con_data[1])
	con_data = con_data.T
    con_data = con_data[good_ix]

    np.savetxt(new_name,con_data,header=hdr)

