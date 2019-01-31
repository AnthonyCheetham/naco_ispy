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
from astropy.table import Table


data_dir = '/data/NACO/Science/'
save_dir = '/home/spectro/cheetham/tmp/rdi_cubes/'

data_dir = '/Users/cheetham/data/naco_data/GTO/Science/'
save_dir = '/Users/cheetham/tmp/rdi_cubes/'

dirs = glob.glob(data_dir+'*/*/ADI/')

output_radius = 25 # so a 50x50 image

info_table=Table(names=['MJD','Airmass','Seeing','r0','t0','WindSpeed','Humidity'])

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

        # Remove any frames that are entirely NaNs:
        good_frames = np.count_nonzero(~np.isnan(im),axis=(1,2)) > 1
        
        # Save them
        outname = save_dir+targname+'_'+targdate+'.fits'
#        pf.writeto(outname,im[good_frames],header=hdr,overwrite=True)
        hdu1 = pf.PrimaryHDU(data=im[good_frames],header=hdr)

        # Add the parangs as a second extension
        header2 = pf.Header()
        header2['UNIT'] = 'degrees'
#        hdulist = pf.append(outname,parangs[good_frames],header2)
        hdu2 = pf.ImageHDU(data=parangs[good_frames],header=header2,name='Parangs')

        # Also get some info from the headers of each of the files
        individual_files = glob.glob(wdir+'../Targ/nomed*.fits')
        for frame_ix,f in enumerate(individual_files):
            if good_frames[frame_ix]:

                head =  pf.getheader(f)

                mjd = head['MJD-OBS']
                airmass = head['AIRMASS']
                seeing = head['HIERARCH ESO TEL AMBI FWHM START']
                r0 = head['HIERARCH ESO AOS RTC DET DST R0MEAN']
                t0 = head['HIERARCH ESO TEL AMBI TAU0']*1000
                windspeed = head['HIERARCH ESO TEL AMBI WINDSP']
                humidity = head['HIERARCH ESO TEL AMBI RHUM']

                info_table.add_row([mjd,airmass,seeing,r0,t0,windspeed,humidity])

        # Now add it to the file
#        pf.append(outname,info_table)
        hdu3 = pf.BinTableHDU(data=info_table,name='SequenceInfo')
        
        hdulist = pf.HDUList(hdus=[hdu1,hdu2,hdu3])
        hdulist.writeto(outname,overwrite=True)
        
        print(wdir)
    else:
        print('One or more cubes missing for '+wdir)

# Now pack it all into a tarball
print 'Saving tarball'
cmd = 'tar -czvf rdi_cubes.tar.gz rdi_cubes'
os.chdir(save_dir+'../')
subprocess.call(cmd,shell=True)