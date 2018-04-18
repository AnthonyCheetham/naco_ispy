#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:59:18 2017

@author: cheetham
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import os
import graphic_companion_fitting_lib
import graphic_contrast_lib
from scipy import interpolate

# Script to run on a directory to calculate the rough astrometry and photometry of
# a companion in NACO-ISPY data

wdir = os.getcwd()+os.sep

image_name = 'smart_annular_pca_derot.fits'
psf_name = 'flux.fits'
tp_name = 'throughput.txt'

plate_scale = 0.0272 # arcsec/pix

print('Using plate scale of: '+str(plate_scale*1000)+' mas/pix')

# Make the image, PSF and throughput global so they can be accessed by the onclick function
global im
global psf
global tp_interp

# Load the reduced image file
if os.access(wdir+image_name,os.F_OK):
    im,im_hdr = pf.getdata(wdir+image_name,header=True)
elif os.access(wdir+'GRAPHIC_PCA'+os.sep+image_name,os.F_OK):
    im,im_hdr = pf.getdata(wdir+'GRAPHIC_PCA'+os.sep+image_name,header=True)
elif os.access(wdir+'ADI'+os.sep+'GRAPHIC_PCA'+os.sep+image_name,os.F_OK):
    im,im_hdr = pf.getdata(wdir+'ADI'+os.sep+'GRAPHIC_PCA'+os.sep+image_name,header=True)
else:
    raise IOError('Reduced image not found')
    
#Try loading the psf frame
if os.access(wdir+psf_name,os.F_OK):
    psf,psf_hdr = pf.getdata(wdir+psf_name,header=True)
elif os.access(wdir+'ADI'+os.sep+psf_name,os.F_OK):
    psf,psf_hdr = pf.getdata(wdir+'ADI'+os.sep+psf_name,header=True)
elif os.access(wdir+'..'+os.sep+'ADI'+os.sep+psf_name,os.F_OK):
    psf,psf_hdr = pf.getdata(wdir+'..'+os.sep+'ADI'+os.sep+psf_name,header=True)
else:
    raise IOError('PSF image not found')
    
# Load the throughput
if os.access(wdir+tp_name,os.F_OK):
    throughput_sep,throughput = np.loadtxt(wdir+tp_name)
elif os.access(wdir+'GRAPHIC_PCA'+os.sep+tp_name,os.F_OK):
    throughput_sep,throughput = np.loadtxt(wdir+'GRAPHIC_PCA'+os.sep+tp_name)
elif os.access(wdir+'ADI'+os.sep+'GRAPHIC_PCA'+os.sep+tp_name,os.F_OK):
    throughput_sep,throughput = np.loadtxt(wdir+'ADI'+os.sep+'GRAPHIC_PCA'+os.sep+tp_name)
else:
    raise IOError('Throughput file not found')
    
# Set up an interpolation function for the throughput
# Use the end values for points outside the range
tp_interp = interpolate.interp1d(throughput_sep,throughput,kind='linear',
                     bounds_error=False,fill_value=(throughput[0],throughput[-1]))


# Rescale the flux of the PSF so it is in the same units as the file
flux_factor = graphic_contrast_lib.correct_flux_frame_naco(im_hdr,psf_hdr)
psf *= flux_factor
print('Multiplying the flux frame by:'+str(flux_factor))

# Function to run after clicking on a part of the image
def onclick(event):
    # Handle clicks outside the image first
    if event.xdata is None:
        return
    coords = [event.xdata,event.ydata]
    coords = np.round(coords).astype(int) # round
    
    cutout_rad = 6
    
    # Cut down the psf frame
    psf_cut = psf[psf.shape[0]/2-cutout_rad:psf.shape[0]/2+cutout_rad,
                  psf.shape[1]/2-cutout_rad:psf.shape[1]/2+cutout_rad]
    
    # Cut out a small part of the image
    xmin = coords[1] - cutout_rad
    xmax = coords[1] + cutout_rad
    ymin = coords[0] - cutout_rad
    ymax = coords[0] + cutout_rad
    cutout = im[xmin:xmax,ymin:ymax]
    
    # Do the fitting
    flux_guess = im[coords[1],coords[0]] / np.max(psf_cut)
    initial_guess = [cutout_rad,cutout_rad,flux_guess]
#    results = companion_fitting.companion_mcmc(im,psf,initial_guess,image_err=1.,n_walkers=50.,n_iterations=1e3,
#                   threads=3,plot=False,burn_in=200)
    
    result = graphic_companion_fitting_lib.simple_companion_leastsq(cutout,psf_cut,initial_guess,
             image_err=1.,threads=3,
             plot=False,method='Powell',fit_to_stdev=False)
    
    # Fix any wrapping in the fitting
    result.x[0] = result.x[0] % (2*cutout_rad)
    if result.x[0] > cutout_rad:
        result.x[0]-=(2*cutout_rad)
    result.x[1] = result.x[1] % (2*cutout_rad)
    if result.x[1] > cutout_rad:
        result.x[1]-=(2*cutout_rad)
        
    # Turn into physical values
    # the results of the fitting are the shift in x and y needed for the psf frame
    # So we need to add cutout_rad (since the psf is shifted that much from zero)
    # and add xmin (or ymin), since that is the distance from the edge of the frame 
    # to the edge of the cutout
    xpos = result.x[1]+ymin+cutout_rad
    ypos = result.x[0]+xmin+cutout_rad
    
    xdist = xpos - im.shape[0]/2
    ydist = ypos - im.shape[0]/2
    
    sep_pix = np.sqrt(xdist**2+ydist**2) 
    sep = sep_pix * plate_scale
    pa = (np.arctan2(-xdist,ydist)*180./np.pi + 720) % 360
    
    # Use the throughput to calculate the contrast
    tp = tp_interp(sep_pix)
    flux_ratio = result.x[2]/tp
    contrast = -2.5*np.log10(flux_ratio)

    # Print the results
    print('  ')
    print('  XPOS:'+str(xdist)+'YPOS:'+str(ydist))
    print('  Sep:'+str(sep),'PA:'+str(pa),'Contrast:'+str(contrast))
    
    # Make a residual image for display
    resids = graphic_companion_fitting_lib.simple_companion_resids(result.x, cutout, psf_cut)
    
    cmax = np.max(np.abs(cutout))
    
    # Display it
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(121)
    implt = ax2.imshow(cutout,vmin=-cmax,vmax=cmax,origin='lower')
    ax2.set_title('Cutout of image')
    fig2.colorbar(implt)
    
    ax3 = fig2.add_subplot(122)
    resplt = ax3.imshow(resids,vmin=-cmax,vmax=cmax,origin='lower')
    ax3.set_title('Residuals after fit')
    fig2.colorbar(resplt)


    plt.show()

print('Click on a companion to fit to its position')
print('Exit all plot windows to finish')

# Set up the plot
fig = plt.figure(1)
ax = fig.add_subplot(111)

# Plot the image
ax.imshow(im,origin='lower')

cid = fig.canvas.mpl_connect('button_press_event',onclick)

plt.show()

#while True:
#    plt.pause(0.5)
fig.canvas.mpl_disconnect(cid)