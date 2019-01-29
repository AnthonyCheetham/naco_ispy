#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 13:13:25 2019

Copy all of the data cubes from the nacodrs into a format suitable for Andre
to run RDI.
He needs:
    - 50 x 50 stamps of each image
    - parallactic angles in a second image extension

@author: cheetham
"""

import numpy as np
import astropy.io.fits as pf
import glob
import os
import subprocess


data_dir = '/data/NACO/Science/'
data_dir = '/Users/cheetham/data/naco_data/GTO/Science/'

dirs = glob.glob(data_dir+'*/*/ADI/')

save_dir = '/Users/cheetham/tmp/rdi_cubes/'

output_radius = 25 # so a 50x50 image

# Loop through directories
for wdir in dirs:
    
    # Check if there is one of each file
    cube_file = wdir+'master_cube_PCA.fits'
    cube_exists = os.access(cube_file, os.F_OK)
    parang_file = wdir+'parallactic_angle.txt'
    parang_exists = os.access(parang_file, os.F_OK)
    
    if cube_exists and parang_exists:
        # Load them both
        im,hdr = pf.getdata(cube_file,header=True)
        parangs = np.loadtxt(parang_file)
        
        # Cut the image down
        cen = [im.shape[1]/2,im.shape[2]/2]
        im = im[:,cen[0]-output_radius:cen[0]+output_radius,
                cen[1]-output_radius:cen[1]+output_radius]
        
        targname = hdr['ESO OBS TARG NAME']
        targdate = hdr['DATE-OBS'].split('T')[0].replace('-','-')
        
        # Save them
        outname = save_dir+targname+'_'+targdate+'.fits'
        try:
            pf.writeto(outname,im,header=hdr,overwrite=True)
        except:
            pf.writeto(outname,im,header=hdr,clobber=True)
        # Add the parangs as a second extension
        header2 = pf.Header()
        header2['UNIT'] = 'degrees'
        hdulist = pf.append(outname,parangs,header2)
        
        print(wdir)
    else:
        print('One or more cubes missing for '+wdir)

# Now pack it all into a tarball
print 'Saving tarball'
cmd = 'tar -czvf rdi_cubes.tar.gz rdi_cubes'
os.chdir(save_dir+'../')
subprocess.call(cmd,shell=True)