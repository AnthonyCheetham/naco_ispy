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
import os,glob,shutil
import astropy.io.fits as pyfits
from astropy.time import Time
from monitor import detect_filetype
import numpy as np
import header

###############

def unorganise_data(folder,prefix='NACO',suffix='.fits',dry_run=True,
                    save_dir='./New/',silent=False):
                        
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
    
    if not silent:
        print(str(len(files))+' files found')
    
    # Loop through the files and move them
    for filename in files:
        
        new_name=filename.rsplit('/',1)[1]
        new_name=save_dir+new_name
        
        if dry_run:
            if not silent:
                print('Moving to: '+new_name)
        else:
            try:
                os.rename(filename,new_name)
            except:
                if not silent:
                    print("Failed to move: "+filename+' to '+new_name)

###############

###############    

def organise_data(folder,prefix='NACO',suffix='.fits',dry_run=True,
              save_dir='./',silent=False):
    """ 
    Looks at all fits files in a folder and moves them to a consistent 
    directory structure.
    
    dry_run: if this is set to True, no files will be moved, but their 
        suggested locations will be printed

    Calib/[date]/[type]/Raw/ : the location of raw calibration data. The current 
        types are: flats.
        
    Science/[date]/[target_name]/Targ/ : the location of the target data
    Science/[date]/[target_name]/Sky/ : the location of the skies (for AGPM data)
    Science/[date]/[target_name]/Flux/ : the location of the flux frames
    """

    # Find the files    
    files=glob.glob(folder+prefix+'*'+suffix)
    files.sort()
    
    if save_dir[-1] != '/':
        save_dir=save_dir+'/'
    
    if len(files) ==0:
        if not silent:
            print('No files to move!')
    else:
        if not silent:
            print(str(len(files))+' files to be sorted')
    
    for filename in files:
        
        hdr=pyfits.getheader(filename)
        
        # Now find the relevant information
        try:
            # this will fail for cal files
            target_name=hdr['HIERARCH ESO OBS TARG NAME'] 
#            target_name=hdr['HIERARCH ESO OBS NAME']# this was changed for astcals
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

        # We have to handle the astrometric calibrators first since they can 
        #  look like target obs
        if 'AstCal' in obstype:
            location=save_dir+'Calib/'+date+'/'+folder_string+'/'
        elif science_flag=='SCIENCE':
            location=save_dir+'Science/'+date+'/'+target_name+'/'+folder_string+'/'
        elif science_flag=='CALIB':            
            location=save_dir+'Calib/'+date+'/'+folder_string+'/'
        elif science_flag=='ACQUISITION':
            location=save_dir+'Science/'+date+'/'+target_name+'/Acq/'
        else:
            if not silent:
                print('Unrecognised science_flag: '+science_flag)
            location=''
        
        if not os.path.isdir(location):
            os.makedirs(location)
            if not silent:
                print('Making new folder')
        
        # Now move the file:
        outname=filename.split('/')[-1] # this line is not particularly robust
        if dry_run:
            if not silent:
                print("Moving to:",location+outname)
        else:
            try:
                os.rename(filename,location+outname)
            except:
                if not silent:
                    print("Failed to move: "+filename+' to '+location+outname)

###############

###############   

def copy_from_archive(archive_folder,save_folder,filename_list='copied_files.txt',
                      prefix='NACO',suffix='.fits',dry_run=True,silent=False):
    ''' Script to be ran that copies files from the archive machine to the working
    directory where the data reduction takes place.
    It tracks the files that have been copied, so that it can be run at any time.
    '''
    
    # Load the list of already copied files
    filename_list_exists=os.access(filename_list,os.F_OK)
    if not filename_list_exists:
        copied_files=[]
    else:
        copied_files=list(np.loadtxt(filename_list,dtype=str))
        
    ncopied=0
    
    # Find the files in the archive
    all_files=glob.glob(archive_folder+'*/'+prefix+'*'+suffix)
    
    # Loop through the files to work out what to do with them
    for filename in all_files:
        
        # Remove the directory path from the filename
        short_fname=filename.split(os.sep)[-1]
        
        # Check if it is in the list of already copied files
        if short_fname in copied_files:
            pass
        else:
            # Move it, and add it to the list
            if dry_run:
                if not silent:
                    print('Copying '+short_fname+' to '+save_folder)
            else:
                shutil.copy2(filename,save_folder)
                copied_files.append(short_fname)
                ncopied+=1
    

    # Save the results
    np.savetxt(filename_list,copied_files,fmt='%50s')

    if not silent:
        print(str(ncopied)+" files successfully copied from archive")

###############

############### 

def check_consistency(all_info,setup_keys,n_min,n_max,obstype='targ',silent=False):
    ''' This program goes through an observation sequence and checks that all 
    of the instrument setup parameters are the same throughout. It will 
    return True if the setup is consistent, and False if one of the params
    changes (e.g. DIT, filter, camera, ND, subarray size).
    Inputs:
        all_info    : list of dictionaries containing the parameters for each file.
        setup_keys  : the dictionary keys to check for consistency.
        n_min,n_max : the minimum and maximum number of files required for a good sequence.
        obstype     : the type of sequence. (Used for printing where the problem was).
        '''

    is_ok=True # innocent until proven guilty
    
    # Check 1 : Number of files
    n_files=len(all_info)
    if (n_files < n_min) or (n_files > n_max):
        is_ok=False
        if not silent:
            print('  Incorrect number of '+obstype+' files')
    
    # Check 2 : Setup consistency
    for key in setup_keys:
        setups=np.unique([info[key] for info in all_info])
        if len(setups) > 1:
            if not silent:
                print('    Different '+obstype+' "'+key+'" values found')
            is_ok=False
        
    return is_ok

###############

###############              

def find_consistent_flux(flux_files,all_info,setup_keys,n_min,n_max,silent=False,
                         dry_run=True,other_flux_ok = -1):
    ''' Find if there is a consistent sequence of flux frames that can be used.
    Since often the DIT will be changed after the OB starts, there may be one
    or two extra files with an incorrect DIT. This should flag them and move
    them away so the sequence is consistent and ready to be reduced.
    '''
    
    is_ok=False # guilty until proven innocent
    
    # Loop over files:
    setups=[]
    for info in all_info:
            
        # Turn the instrument parameters into a set so we can check which files match
        arr=[]
        for key in setup_keys:
            arr.append(info[key])
        setups.append(set(arr))
    
    # Now loop again and check how many match each one
    for test_setup in setups:
        matching_files = [ix for ix, this_setup in enumerate(setups) if this_setup == test_setup]
        other_files = [ix for ix, this_setup in enumerate(setups) if this_setup != test_setup]
        n_matching = len(matching_files)
        
        if ((n_matching >= n_min) & (n_matching <= n_max)) or (n_matching == other_flux_ok):
            if not silent:
                print('  Consistent setup found!')
                
            is_ok=True
            
            # Just assume there is only one "good" sequence and move the remaining files away
            for bad_ix in other_files:
                old_name=flux_files[bad_ix]
                new_name=old_name.replace('/Flux/','/Flux/bad_files/')
                
                # check if the directory exists
                bad_files_dir=old_name.rsplit(os.sep,1)[0]+os.sep+'bad_files'+os.sep
                if not os.access(bad_files_dir,os.F_OK):
                    os.mkdir(bad_files_dir)
                    
                # Move the bad file
                if not dry_run:
                    os.rename(flux_files[bad_ix],new_name)
            
            break
    if is_ok == False and (not silent):
        print("  Didn't find consistent setup!")
        
    return is_ok
    
###############

############### 
            
         
def check_directory_consistency(folder,setup_keys,silent=False,dry_run=True,
                n_flux_min=6,n_flux_max=6,n_targ_min=10,n_targ_max=1e6,
                other_flux_ok = 3):
    ''' Check a directory to see if the observing sequence is consistent.
    This is just a wrapper that combines check_consistency and 
    header.get_info_from_files

    other_flux_ok: An alternative number of flux frames that is acceptable. 
    For example, as long as a full dither sequence is completed, any multiple
    of n_dither is ok. (3,6,9...)
    
    silent : if False, will print lots of diagnostic info at each step
    dry_run: if True, will not move or rename any files.
    
    '''
    if not silent:    
        print(folder)
    
    ########
    # Check the flux frames
    ########
    flux_files=glob.glob(folder+'Flux/NACO*.fits')

    # Get the info for these files
    flux_info=header.get_info_from_files(flux_files,header.target_keys,wdir=folder)
    
    flux_ok=check_consistency(flux_info,setup_keys,n_flux_min,n_flux_max,
                obstype='flux',silent=silent)
    
    # For the flux, we can check if there is a nice setup
    if not flux_ok:
        flux_ok=find_consistent_flux(flux_files,flux_info,setup_keys,n_flux_min,
                 n_flux_max,silent=silent,dry_run=dry_run,other_flux_ok = other_flux_ok)

    ########
    # Check the targ/sky frames
    ########
    targ_files=glob.glob(folder+'Targ/NACO*.fits')
    sky_files=glob.glob(folder+'Sky/NACO*.fits')
    targ_files+=sky_files

    # Get the info for these files
    targ_info=header.get_info_from_files(targ_files,
                header.target_keys,wdir=folder)
    
    targ_ok=check_consistency(targ_info,setup_keys,n_targ_min,n_targ_max,
                obstype='targ',silent=silent)

    # Print the results (temp)
    if not silent:
        print('  Final status: '+str(flux_ok)+' '+str(targ_ok))
    
    return flux_ok,targ_ok    

###############
    
# As a test, run it on the watch folder
#folder='/Users/cheetham/data/data_archive/naco_watch_folder/New/'
#folder='/Users/cheetham/data/naco_data/GTO/New/'
#save_dir='/Users/cheetham/data/naco_data/GTO/'

#hdr=unorganise_data(save_dir,dry_run=False,save_dir=folder,prefix='NACO.',suffix='.fits')
#hdr=organise_data(folder,dry_run=False,save_dir=save_dir,prefix='NACO',suffix='.fits')

#archive_folder='/Users/cheetham/data/data_archive/naco_archive/'
#save_folder='/Users/cheetham/data/data_archive/naco_watch_folder/New/'
#
#test=copy_from_archive(archive_folder,save_folder,dry_run=False)
