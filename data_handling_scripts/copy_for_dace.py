# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 10:35:43 2017

This program copies the final products for a dataset into a directory and
format that is easily accessible by DACE

We want:
    - high_contrast.fits
    - non_saturated.fits
    - detection_limits.rdb
    - signal_to_noise_map.fits

Designed specifically for NACO-ISPY

@author: cheetham
"""

import numpy as np
import glob,os,sys,argparse
import astropy.io.fits as pf
from astropy import table
import naco_ispy


parser = argparse.ArgumentParser(description='This program copies all processed NACO ISPY data to the right format and directory for DACE.')

parser.add_argument('-dry_run', dest="dry_run",action='store_const',const=True,
                    default=False, help='Dont actually copy the files, but print the commands it will do')   
parser.add_argument('--num', action="store",  dest="num", type=int, default=None,
                    help='Maximum number of datasets to copy')
parser.add_argument('--skip',action="store",  dest="skip", type=int, default=0,
                    help='Skip this many datasets before starting the copying')

# Get the input arguments
args = parser.parse_args()
num = args.num
skip = args.skip
dry_run = args.dry_run

# Define the directories etc here
data_folder = '/data/NACO/'
db_filename = '/data/NACO/obs_table.dat'
# data_folder='/Users/cheetham/data/naco_data/GTO/'
# db_filename='/Users/cheetham/data/naco_data/GTO/obs_table.dat'



###############

###############

def add_dace_headers(hdr,targ_row):
    ''' Add the necessary keywords to the header of a file so that DACE
    knows how to catalogue it.
    
    All of the information is already in headers, but we need to rearrange it 
    to save DACE some effort'''

    # Fix the name of some association members
    name = hdr['HIERARCH ESO OBS TARG NAME']
    for association in ['Sco','Oph','Hya','Sgr','Cen','Tau','Aur','CrA']:
        name = name.replace(association,' '+association)
    
    hdr['HIERARCH DACE INSTRUMENT NAME'] = 'NACO'
    hdr['HIERARCH DACE OBS PROGRAM'] = 'NACO-ISPY'  #hdr['HIERARCH ESO OBS PROG ID']
    hdr['HIERARCH DACE OBS DATE'] = targ_row['Date']
    hdr['HIERARCH DACE OBS RJD'] = np.round(hdr['MJD-OBS'] + 0.5,2) # convert from MJD to Reduced JD
    hdr['HIERARCH DACE OBJECT NAME'] = name
    hdr['HIERARCH DACE OBJECT RA-DEG'] = hdr['RA']
    hdr['HIERARCH DACE OBJECT DEC-DEG'] = hdr['Dec']
    hdr['HIERARCH DACE IMAGING CALIBRATION'] = ''
    hdr['HIERARCH DACE IMAGING FILTER'] = 'Lp'

    # AGPM or SatPSF?
    if targ_row['AGPM'] == 'True':
        imaging_type = 'AGPM'
    else:
        imaging_type = 'Saturated PSF'
    hdr['HIERARCH DACE IMAGING TYPE'] = imaging_type

    # Extra keywords for assessing data quality:
    hdr['HIERARCH DACE IMAGING FIELD ROTATION'] = np.round(targ_row['FieldRotation'],1)
    hdr['HIERARCH DACE IMAGING SEEING R0'] = targ_row['r0']
    hdr['HIERARCH DACE IMAGING SEEING T0'] = targ_row['t0']
    hdr['HIERARCH DACE IMAGING DIT'] = targ_row['ExpTime']
    
    return hdr
    
###############
    
###############

def copy_dataset_for_dace(dir_in, targ_row, flux_file = 'flux_cube.fits',
#          adi_file = 'GRAPHIC_PCA/smart_annular_pca_derot.fits',
        adi_file = 'GRAPHIC_PCA/pca_multimodes.fits',
        contrast_file = 'GRAPHIC_PCA/contrast.txt',
        snr_file = 'GRAPHIC_PCA/snr_map.fits',
        dir_out = ''):
    ''' Copy the final products for a single dataset to a directory accessible 
    by DACE
    targ_row is the row from the database containing the target
    '''
    
    # Check that the files are there    
    flux_exists = os.access(dir_in+flux_file, os.F_OK)
    adi_exists = os.access(dir_in+adi_file, os.F_OK)
    contrast_exists = os.access(dir_in+contrast_file, os.F_OK)
    snr_exists = os.access(dir_in+snr_file, os.F_OK)
    
    # Check that the output directory exists and create it
    if not os.access(dir_out, os.F_OK):
        os.makedirs(dir_out)
    
    # Now go through them one by one
    if flux_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+flux_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row)
        
        # Save it out
        pf.writeto(dir_out+'non_saturated.fits',im,header=hdr,clobber=True,
                   output_verify='silentfix')
    
    if adi_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+adi_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row)
        
        # Save it out
        pf.writeto(dir_out+'high_contrast.fits',im,header=hdr,clobber=True,
                   output_verify='silentfix')
                   
    if snr_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+snr_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row)
        
        # Save it out
        pf.writeto(dir_out+'signal_to_noise_map.fits',im,header=hdr,clobber=True,
                   output_verify='silentfix')
        
    
    if contrast_exists:
        
        # Load it
        sep,con = np.loadtxt(dir_in+contrast_file)
        
        # Remove any NaNs
        sep[np.isnan(con) == False]
        con[np.isnan(con) == False]
        
        # If there are no values left, print an error
        if len(con) == 0:
            print('    Contrast curve is only NaNs! Not copied...')
        else:
            
            # convert to rdb
            data={'Separation':sep,'Contrast':con}
            tab = table.Table(data,names=['Separation','Contrast'])
            
            tab.write(dir_out+'detection_limits.rdb',format='rdb')
    print('    '+str(flux_exists)+' '+str(adi_exists)+' '+str(snr_exists)+' '+str(contrast_exists))

###############

# Load the obs_db 
# First, load the target database
obs_db = naco_ispy.databases.obs_table(filename=db_filename, data_folder=data_folder)

# Now loop through all of the entries
data = obs_db.data

for targ_row in data[skip:num]:

    # Check if it has been processed
    processed = np.bool(targ_row['ADIProcessed'])

    # If it has been processed, then copy it over
    if processed:

        # Where do we want to save it
        band = targ_row['Band']
        if band =='L_prime':
            band = 'Lp'
        else:
            continue
        dir_out = data_folder + 'Final_products'+os.sep+'NACO'+os.sep + targ_row['TargetName'] +os.sep+ targ_row['Date']+os.sep+band+os.sep

        dir_in = targ_row['Location']+'ADI'+os.sep


        if dry_run:
            # print('  Would copy from : '+targ_row['Location'])
            print('target:'+dir_out)
        else:
            print('  Copying from: '+targ_row['Location'])
            copy_dataset_for_dace(dir_in,targ_row,dir_out=dir_out)

# Find all of the relevant files
#star='hd15115/'
#dir_in = '/Users/cheetham/tmp/'+star
#dir_out = '/Users/cheetham/data/naco_data/GTO/DACE/'+star
# test it
#hdr=copy_dataset_for_dace(dir_in,dir_out=dir_out)#,adi_file = 'GRAPHIC_PCA/smart_annular_pca_derot.fits')

# dir_in = '/data/NACO/Science/2015-12-15/HD17925/ADI/'
# dir_out = '/data/NACO/Final_products/NACO/HD17925/2015-12-15/Lp/'
# hdr=copy_dataset_for_dace(dir_in,dir_out=dir_out)

# dir_in = '/data/NACO/Science/2016-11-08/HD15115/ADI/'
# dir_out = '/data/NACO/Final_products/NACO/HD15115/2016-11-08/Lp/'
# hdr=copy_dataset_for_dace(dir_in,dir_out=dir_out)

# dir_in = '/data/NACO/Science/2016-05-02/HD179218/ADI/'
# dir_out = '/data/NACO/Final_products/NACO/HD179218/2016-05-02/Lp/'
# hdr=copy_dataset_for_dace(dir_in,dir_out=dir_out)

