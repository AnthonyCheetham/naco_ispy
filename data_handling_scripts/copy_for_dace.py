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
import glob


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
save_folder = '/data/NACO/NACO-ISPYDRS/DRS-1.0/reduced/'
# data_folder='/Users/cheetham/data/naco_data/GTO/'
# db_filename='/Users/cheetham/data/naco_data/GTO/obs_table.dat'



###############

###############

def add_dace_headers(hdr,targ_row,postproc_algorithm=None):
    ''' Add the necessary keywords to the header of a file so that DACE
    knows how to catalogue it.
    
    All of the information is already in headers, but we need to rearrange it 
    to save DACE some effort'''

    # Fix the name of some association members
    name = hdr['HIERARCH ESO OBS TARG NAME']
    for association in ['Sco','Oph','Hya','Sgr','Cen','Tau','Aur','CrA','Lup']:
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

    # Add postprocessing algorithm if neeeded
    if not (postproc_algorithm is None):
        hdr['HIERARCH DACE ALGORITHM'] = postproc_algorithm
    
    return hdr
    
###############
    
###############

def copy_dataset_for_dace_flux(dir_in, targ_row, flux_file = 'flux_cube.fits',
        dir_out = '', archive_name = '', postproc_algorithm='gpca'):
    ''' Copy the flux frame for a single dataset to a directory accessible 
    by DACE
    targ_row is the row from the database containing the target
    '''
    
    # Check that the files are there
    flux_exists = os.access(dir_in+flux_file, os.F_OK)

    nonsat_filename = dir_out + archive_name+'_ns.fits' # Only need one
    
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
        pf.writeto(nonsat_filename,im,header=hdr,overwrite=True,
                   output_verify='silentfix')

def copy_dataset_for_dace_pca(dir_in, targ_row, flux_file = 'flux_cube.fits',
        adi_file = 'GRAPHIC_PCA/pca_multimodes.fits',
        contrast_file = 'GRAPHIC_PCA/contrast.txt',
        snr_file = 'GRAPHIC_PCA/snr_map.fits',
        dir_out = '',archive_name = '',
        postproc_algorithm='gpca'):
    ''' Copy the final products for a single dataset to a directory accessible 
    by DACE
    targ_row is the row from the database containing the target
    '''
    
    # Check that the files are there
    adi_exists = os.access(dir_in+adi_file, os.F_OK)
    contrast_exists = os.access(dir_in+contrast_file, os.F_OK)
    snr_exists = os.access(dir_in+snr_file, os.F_OK)

    highcontrast_filename = dir_out + archive_name+'_'+postproc_algorithm+'_hc.fits'
    snr_filename =          dir_out + archive_name+'_'+postproc_algorithm+'_snr.fits'
    lims_filename =         dir_out + archive_name+'_'+postproc_algorithm+'_dl.rdb'
    
    # Check that the output directory exists and create it
    if not os.access(dir_out, os.F_OK):
        os.makedirs(dir_out)

    # Now go through them one by one
    if adi_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+adi_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row,postproc_algorithm='PCA')
        
        # Save it out
        pf.writeto(highcontrast_filename,im,header=hdr,overwrite=True,
                   output_verify='silentfix')
                       
    if snr_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+snr_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row,postproc_algorithm='PCA')
        
        # Save it out
        pf.writeto(snr_filename,im,header=hdr,overwrite=True,
                   output_verify='silentfix')
        
    
    if contrast_exists:
        # Load it
        sep,con = np.loadtxt(dir_in+contrast_file)
        
        # Remove any NaNs
        sep = sep[np.isnan(con) == False]
        con = con[np.isnan(con) == False]
        
        # If there are no values left, print an error
        if len(con) == 0:
            print('    Contrast curve is only NaNs! Not copied...')
        else:
            
            # convert to rdb
            data={'Separation':sep,'Contrast':con}
            tab = table.Table(data,names=['Separation','Contrast'])
            
            tab.write(lims_filename,format='rdb',overwrite=True)


def copy_dataset_for_dace_cadi(dir_in, targ_row,adi_file = 'cADI/smart_adi_derot.fits',archive_name='',
    dir_out = '',postproc_algorithm='cadi'):
    ''' Copy the cADI final product for a single dataset to a directory accessible 
    by DACE
    targ_row is the row from the database containing the target
    '''
    
    # Check that the files are there
    adi_exists = os.access(dir_in+adi_file, os.F_OK)

    highcontrast_filename = dir_out + archive_name+'_'+postproc_algorithm+'_hc.fits'
    
    # Check that the output directory exists and create it
    if not os.access(dir_out, os.F_OK):
        os.makedirs(dir_out)

    # Now go through them one by one
    if adi_exists:
        # Load it
        im,hdr = pf.getdata(dir_in+adi_file,header=True)
        
        # Fix the header
        hdr = add_dace_headers(hdr,targ_row,postproc_algorithm='cADI')
        
        # Save it out
        pf.writeto(highcontrast_filename,im,header=hdr,overwrite=True,
                   output_verify='silentfix')

def copy_dataset_for_dace_trap(dir_in, targ_row, dir_out = '',archive_name='',postproc_algorithm='trap',header=None):
    ''' Copy the TRAP/DICPM final products for a single dataset to a directory accessible 
    by DACE
    targ_row is the row from the database containing the target
    '''

    # We have to find the files in case they have a different name
    adi_file = glob.glob(dir_in+'DICPM'+os.sep+'detection*.fits')
    snr_file = glob.glob(dir_in+'DICPM'+os.sep+'norm_detection*.fits')
    contrast_file = glob.glob(dir_in+'DICPM'+os.sep+'contrast_table*.fits')
    
    # Check that the files are there
    adi_exists = len(adi_file) > 0
    snr_exists = len(snr_file) > 0
    contrast_exists = len(contrast_file) > 0

    highcontrast_filename = dir_out + archive_name+'_'+postproc_algorithm+'_hc.fits'
    snr_filename =          dir_out + archive_name+'_'+postproc_algorithm+'_snr.fits'
    lims_filename =         dir_out + archive_name+'_'+postproc_algorithm+'_dl.rdb'
    
    # Check that the output directory exists and create it
    if not os.access(dir_out, os.F_OK):
        os.makedirs(dir_out)

    # Now go through them one by one
    if adi_exists:
        # Load it
        im = pf.getdata(adi_file[0])[0] # First frame is the best-fit contrast
        
        # Fix the header
        header = add_dace_headers(header,targ_row,postproc_algorithm='TRAP')
        
        # Save it out
        pf.writeto(highcontrast_filename,im,header=header,overwrite=True,
                   output_verify='silentfix')

                       
    if snr_exists:
        # Load it
        im = pf.getdata(snr_file[0])
        
        # Fix the header
        header = add_dace_headers(header,targ_row,postproc_algorithm='TRAP')
        
        # Save it out
        pf.writeto(snr_filename,im,header=header,overwrite=True,
                   output_verify='silentfix')
        
    
    if contrast_exists:
        # Load it
        hdulist = pf.open(contrast_file[0])
        tab = hdulist[1].data

        sep = tab['sep (mas)']/1000. # in arcsec
        con_ratio = tab['contrast_50'] # contrast ratio
        con = -2.5*np.log10(con_ratio) # contrast (mag)
        
        # Remove any NaNs
        sep = sep[np.isnan(con) == False]
        con = con[np.isnan(con) == False]
        
        # If there are no values left, print an error
        if len(con) == 0:
            print('    Contrast curve is only NaNs! Not copied...')
        else:
            
            # convert to rdb
            data={'Separation':sep,'Contrast':con}
            tab = table.Table(data,names=['Separation','Contrast'])
            
            tab.write(lims_filename,format='rdb',overwrite=True)

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

        dir_out = save_folder + targ_row['Date'] + os.sep

        dir_in = targ_row['Location']+'ADI'+os.sep


        if dry_run:
            # print('  Would copy from : '+targ_row['Location'])
            print('target:'+dir_out)
        else:
            print('  Copying from: '+targ_row['Location'])
            # Get the archive name from the master_cube file
            if os.access(dir_in+'master_cube_PCA.fits',os.F_OK):
                header = pf.getheader(dir_in+'master_cube_PCA.fits')
            else:
                # Otherwise skip this one...
                print('Didnt find cleaned cube for {0}. Skipping.'.format(dir_in))
                continue
            archive_name = header['ARCFILE'].replace('.fits','')

            # Flux frame
            copy_dataset_for_dace_flux(dir_in,targ_row,dir_out=dir_out,archive_name=archive_name)

            # GRAPHIC_PCA
            copy_dataset_for_dace_pca(dir_in,targ_row,dir_out=dir_out,archive_name=archive_name)

            # TRAP
            copy_dataset_for_dace_trap(dir_in,targ_row,dir_out=dir_out,archive_name=archive_name,header=header)

            # cADI
            copy_dataset_for_dace_cadi(dir_in,targ_row,dir_out=dir_out,archive_name=archive_name)

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

