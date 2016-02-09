# -*- coding: utf-8 -*-
"""
This module contains all of the tools useful for sorting the data
"""
import numpy as np
import sys,os,glob,time,astropy
import astropy.io.fits as pyfits
from astropy.time import Time

# Some constants
nexpo_limit=2 # If nexpo > 2, this indicates that it is a target observation. otherwise, sky
obstime_limits=[0.1,0.5] # all target/sky observations have exp times in this range

def detect_filetype(hdr,get_folder_string=False):
    ''' Works out what kind of file it is based on the header'''
    type_flag=hdr['HIERARCH ESO DPR TYPE']
    expt=hdr['EXPTIME'] # exposure time.

    try:
        nexpo=hdr['HIERARCH ESO SEQ NEXPO']
    except:
        nexpo=0
    
    
    # Now format all of these strings      
    if type_flag=='OBJECT':
        # We need to work out which of the "OBJECT" frames are skies, flux 
        #  frames and actual target observations.
        # For the GTO, we can use the number of exposures. Skies are a single cube for the AGPM.
        # For non-AGPM observations, we don't need separate skies, so label them all as Targ.
        if nexpo > nexpo_limit:
            obstype='Target'
            folder_string='Raw'
        elif (expt < obstime_limits[1]) and (expt > obstime_limits[0]):
            obstype='Sky'
            folder_string='Sky'
        else:
            obstype='Flux'
            folder_string='Flux'
    elif type_flag=='SKY':
        obstype='Sky'
        folder_string='Sky'
    elif type_flag=='FLAT,LAMP' or type_flag=='FLAT,SKY':
        obstype='Flat'
        folder_string='Flats'
    
    # We don't actually use any of the following types, but I thought we might as well
    # put them somewhere
    elif type_flag=='STD':
        obstype='Std'
        folder_string='STD'
    elif type_flag=='DARK':
        obstype='Dark'
        folder_string='Dark'
    else:
        # Put all of the unknown file types into a single folder to make it easy
        print 'Unrecognised DPR type:',type_flag
        obstype='Unknown'
        folder_string='Uncategorized'
    
    if get_folder_string:
        return obstype,folder_string
    else:
        return obstype
    

def organise_data(folder,prefix='NACO_CORO',suffix='.fits',dry_run=True,
              save_dir='/Users/cheetham/data/data_archive/naco_watch_folder/',
              flux_exp_cutoff=0.15):
    """ 
    Looks at all fits files in a folder and moves them to a consistent 
    directory structure.
    
    dry_run: if this is set to True, no files will be moved, but their 
        suggested locations will be printed
    flux_exp_cutoff: The smallest

    Calib/[date]/[type]/Raw/ : the location of raw calibration data. The current 
        types are: flats.
        
    Science/[date]/[target_name]/Raw/ : the location of the raw science data
    Science/[date]/[target_name]/Sky/ : the location of the skies (for AGPM data)
    Science/[date]/[target_name]/Flux/ : the location of the flux frames
    """

    # Find the files    
    files=glob.glob(folder+prefix+'*'+suffix)
    
    if save_dir[-1] != '/':
        save_dir=save_dir+'/'
    
    if len(files) ==0:
        print 'No files found!'
    else:
        print len(files),'files to be sorted'
    
    for filename in files:
        
        hdr=pyfits.getheader(filename)
        
        # Now find the relevant information
        try:
            target_name=hdr['HIERARCH ESO OBS TARG NAME'] # this will fail for cal files
        except:
            target_name=''
        mjd=hdr['MJD-OBS']
        science_flag=hdr['HIERARCH ESO DPR CATG']
                
        obstype,folder_string=detect_filetype(hdr,get_folder_string=True)
        
        # Get the date in a nice format and convert so that observations made in
        #  the same night have the same date.
        # Chile is UTC-3 all year round (no daylight savings), so doing
        #  (time-9hrs) moves midnight Chilean time to midday UTC. Then, the raw date
        #  (i.e. ignoring the time) will give the date of the *night*
        
        date_time=Time(mjd-9./24.,format='mjd')
        date,dummy=date_time.iso.split(' ')
        
        if science_flag=='SCIENCE':
            location=save_dir+'Science/'+date+'/'+target_name+'/'+folder_string+'/'
        elif science_flag=='CALIB':            
            location=save_dir+'Calib/'+date+'/'+folder_string+'/'
        elif science_flag=='ACQUISITION':
            location=save_dir+'Science/'+date+'/'+target_name+'/'
        else:
            print('Unrecognised science_flag:',science_flag)
            location=''
        
        if not os.path.isdir(location):
            os.makedirs(location)
            print 'Making new folder'
        
        # Now move the file:
        outname=filename.split('/')[-1] # this line is not particularly robust
        if dry_run:
            print("Moving to:",location+outname)
        else:
            try:
                os.rename(filename,location+outname)
            except:
                print("Failed to move:",filename,'to',location+outname)
        
    return hdr
    
# As a test, run it on the watch folder
folder='/Users/cheetham/data/data_archive/naco_watch_folder/New/'
folder='/Users/cheetham/data/naco_data/naco_51Eri/Raw/data_with_raw_calibs/'

hdr=organise_data(folder,dry_run=False,save_dir='/Users/cheetham/data/naco_data/naco_51Eri/Raw/',prefix='NACO.')