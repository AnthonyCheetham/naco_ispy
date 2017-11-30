#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 18:04:41 2017

@author: cheetham
"""

# One-time script to update the contrast curves and SNR maps,
# so that they use the median absolute deviation instead of the standard deviation


import graphic_contrast_lib
import astropy.io.fits as pf
import numpy as np

# Assume PCA and fake-planet injection have been run. 

# Options
smooth_image_length=1.25
median_filter_length=20
fwhm = 4.
plate_scale = 0.027
n_radii = 200

# Filenames
wdir = 'ADI/'
save_dir = 'GRAPHIC_PCA/'
psf_file = wdir+'flux.fits'
image_file = wdir+'GRAPHIC_PCA/smart_annular_pca_derot.fits'
cube_file = wdir+'master_cube_PCA.fits'
parang_file = wdir+'parallactic_angle.txt'


# Set up the names of the output files
fp_name = wdir+save_dir+'fake_planets.fits'
fp_pca_name = wdir+save_dir+'fake_planets_pca.fits'
fp_derot_name = wdir+save_dir+'fake_planets_derot.fits'
throughput_file = wdir+save_dir+'throughput.txt'
all_throughput_file = wdir+save_dir+'all_throughputs.txt'
contrast_im_file = wdir+save_dir+'contrast_im.fits'
contrast_file = wdir+save_dir+'contrast.txt'
snr_map_file = wdir+save_dir+'snr_map.fits'
noise_file = wdir+save_dir+'noise.txt'


##############

# Load the psf image, pca subtracted image and data cubes
# (and their headers)
psf_frame,psf_header = pf.getdata(psf_file,header=True)
image,header = pf.getdata(image_file,header=True)
cube,cube_header = pf.getdata(cube_file,header=True)

# Load the parallactic angles
parangs_deg = np.loadtxt(parang_file)
parangs_rad = parangs_deg*np.pi/180.

pca_r_min = header['HIERARCH GC PCA RMIN']
pca_r_max = header['HIERARCH GC PCA RMAX']
r_max = pca_r_max

# Correct the PSF for ND filters etc
# This needs to be properly implemented for NICI, SCEXAO etc.
if 'INSTRUME' in psf_header.keys():
    instrument = psf_header['INSTRUME'].strip()
else:
    raise Exception('Unknown instrument in GRAPHIC_contrast_pca')
    
if instrument == 'NAOS+CONICA':
    flux_factor = graphic_contrast_lib.correct_flux_frame_naco(cube_header,psf_header)
elif instrument == 'SPHERE':
    flux_factor = graphic_contrast_lib.correct_flux_frame_sphere(cube_header,psf_header)
else:
    raise Exception('Unknown instrument in GRAPHIC_contrast_pca: '+instrument)
    
psf_frame *= flux_factor
print('Multiplying the flux frame by:'+str(flux_factor))

# Apply some cosmetics to the image and flux frame before measuring the contrast
graphic_contrast_lib.prepare_detection_image(image_file,save_name=contrast_im_file,
                         smooth_image_length=smooth_image_length,
                         median_filter_length=median_filter_length)
# Apply the same smoothing to the flux frame (but keep in memory rather than saving it)
# However, DON'T apply the median filter, since that would reduce the flux of the PSF much
# more than the image, due to the high SNR
psf_frame = graphic_contrast_lib.prepare_detection_image(psf_frame,
                         smooth_image_length=smooth_image_length)

# Calculate the contrast
graphic_contrast_lib.contrast_curve(contrast_im_file,psf_frame,
       r_min=pca_r_min,r_max=r_max,fwhm=fwhm,plate_scale=plate_scale,
       self_subtraction_file=throughput_file,save_contrast = contrast_file,
       n_radii=n_radii,save_noise = noise_file,mad=True,robust_sigma=False)
       
# Make a signal-to-noise map
graphic_contrast_lib.snr_map(contrast_im_file,noise_file,
                  remove_planet=False,save_name=snr_map_file)
