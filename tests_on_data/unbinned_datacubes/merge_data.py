#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 16:08:04 2017

@author: cheetham
"""
# Merge all of the files into one

import numpy as np
import glob
import argparse
import graphic_nompi_lib_330 as graphic_nompi_lib
import os
import astropy.io.fits as pf

#os.chdir('/Users/cheetham/data/naco_data/GTO/Science/2017-05-19/HIP65426/Targ/')

parser = argparse.ArgumentParser(description='Read the parangs from the RDB files and save them out as txt files')
parser.add_argument('--output_dir', action="store", dest="output_dir",  default='./unbinned/',
                    help='Output directory')
parser.add_argument('--parang_input_dir', action="store", dest="parang_input_dir",  default='./cube-info/',
                    help='Input directory for parallactic angle rdb file')
parser.add_argument('--info_prefix', action="store", dest="info_prefix",  default='all_info_r_',
                    help='Filename prefix for info files')
parser.add_argument('--fits_prefix', action="store", dest="fits_prefix",  default='r_',
                    help='Filename prefix for fits cubes')

args = parser.parse_args()

parang_input_dir = args.parang_input_dir
output_dir = args.output_dir
info_prefix = args.info_prefix
fits_prefix = args.fits_prefix

# Check that output_dir exists
if not os.access(output_dir,os.F_OK):
    os.mkdir(output_dir)

###########

# Find the rdb files
rdb_files = glob.glob(parang_input_dir+info_prefix+'*.rdb')
rdb_files = np.sort(rdb_files)

parang_array = []

# Loop through them
for rdb_file in rdb_files:
    
    # Read it
    data_info = graphic_nompi_lib.read_rdb(rdb_file)
    
    # Add to the array
    parangs = data_info['paralactic_angle']
    
    parang_array.extend(parangs)
    
# Write it out
parang_out_name = output_dir+ 'parangs.txt'
np.savetxt(parang_out_name,parang_array)

###########

# Find the fits files
fits_files = glob.glob(fits_prefix+'*.fits')
fits_files = np.sort(fits_files)

up_to_ix = 0

for ix,fits_file in enumerate(fits_files):
    
    # Read it
    if ix == 0:
        cube,hdr = pf.getdata(fits_file,header=True)
        
        # Instead of constantly resizing the array, guess how big it will be and then
        # remove the extra frames at the end
        n_frames_est = np.int(len(fits_files)*cube.shape[0]*1.2)
        cube_out = np.zeros((n_frames_est,cube.shape[1],cube.shape[2]))
    else:
        cube = pf.getdata(fits_file)

    # Quick check that it will fit in the array (could make it adapt the array size if needed)
    if (up_to_ix+cube.shape[0]) > cube_out.shape[0]:
        raise Exception('The pre-made output array is too small. Fix merge_data to increase the size')
        
    cube_out[up_to_ix:up_to_ix+cube.shape[0]] = cube
    
    up_to_ix += cube.shape[0]
        
# Remove the unneeded frames
cube_out = cube_out[:up_to_ix]

# Write it out
fits_save_name = output_dir+'master_cube_PCA.fits'
pf.writeto(fits_save_name,cube_out,header=hdr,clobber=True)
