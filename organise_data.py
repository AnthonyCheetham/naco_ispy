# -*- coding: utf-8 -*-
"""
This module contains all of the tools useful for sorting the data.
Note that the important detect_filetype module is in monitor.py

The two main modules are a data sorter, and a data un-sorter.
The unsorter should be useful in case you want to change the file structure or
to undo any changes.

For both modules, you can set dry_run=True to make it print any proposed
file moves without actually moving them
"""
import os,glob
import astropy.io.fits as pyfits
from astropy.time import Time
from monitor import detect_filetype

###############

def unorganise_data(folder,prefix='NACO_',suffix='.fits',dry_run=True,
                    save_dir='/Users/cheetham/data/naco_data/GTO/New/'):
                        
    ''' 
    This undoes the effects of organise_data, an puts all of the files into 
    a single directory. This is most useful for testing changes and fixing bugs
    in organise_data
    '''

    # Add a slash in case it is missing    
    if folder[-1]!='/':
        folder=folder+'/'
        
    files=glob.glob(folder+'*/*/*/*/'+prefix+'*'+suffix)
    calib_files=glob.glob(folder+'*/*/*/'+prefix+'*'+suffix)
    files.extend(calib_files)
    
    print len(files),'Files found'
    
    # Loop through the files and move them
    for filename in files:
        
        new_name=filename.rsplit('/',1)[1]
        new_name=save_dir+new_name
        
        if dry_run:
            print 'Moving to:',new_name
        else:
            try:
                os.rename(filename,new_name)
            except:
                print("Failed to move:",filename,'to',new_name)

###############

###############    

def organise_data(folder,prefix='NACO_',suffix='.fits',dry_run=True,
              save_dir='/Users/cheetham/data/naco_data/GTO/New/'):
    """ 
    Looks at all fits files in a folder and moves them to a consistent 
    directory structure.
    
    dry_run: if this is set to True, no files will be moved, but their 
        suggested locations will be printed

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
        # MJD is in UTC, so split by midday by subtracting 12hrs and then 
        # taking the *raw date* .
        
        date_time=Time(mjd-0.5,format='mjd')
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
    
###############
    
# As a test, run it on the watch folder
#folder='/Users/cheetham/data/data_archive/naco_watch_folder/New/'
#folder='/Users/cheetham/data/naco_data/GTO/New/'
#save_dir='/Users/cheetham/data/naco_data/GTO/'

#hdr=unorganise_data(save_dir,dry_run=False,save_dir=folder,prefix='NACO.',suffix='.fits')
#hdr=organise_data(folder,dry_run=False,save_dir=save_dir,prefix='NACO',suffix='.fits')
