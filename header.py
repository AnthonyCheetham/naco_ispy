# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:32:19 2015

Module that reads the headers of NACO files and performs various operations

Todo:
 - Program to make a summary file

@author: cheetham
"""

import astropy.io.fits as pf
import glob

#####################
target_keys={'tname':'HIERARCH ESO OBS NAME',
      'datetime':'DATE-OBS',
      'RA':'RA',
      'Dec':'DEC',
      'obstype':'HIERARCH ESO DPR TYPE',
      'filter1':'HIERARCH ESO INS OPTI6 ID',
      'filter2':'HIERARCH ESO INS OPTI5 ID',
      'dit':'HIERARCH ESO DET DIT',
      'ndit':'HIERARCH ESO DET NDIT',
      'NAXIS':'NAXIS',
      'ax1':'NAXIS1',
      'ax3':'NAXIS3',
      'nexpo':'HIERARCH ESO SEQ NEXPO',
      'rotstart':'HIERARCH ESO ADA ABSROT START',
      'rotend':'HIERARCH ESO ADA ABSROT END',
      'pastart':'HIERARCH ESO TEL PARANG START',
      'paend':'HIERARCH ESO TEL PARANG END',
      'alt':'ESO TEL ALT',
      'camera':'HIERARCH ESO INS OPTI7 ID',
      'nd':'HIERARCH ESO INS OPTI3 ID'
      }

#####################

#####################

def get_info_from_files(files,keys_dict,wdir=''):
    ''' Loop through a list of files, store the value of each header keyword.
    wdir is used only to truncate the folder name from the file name
    
    In the future it may be possible to use fitsio to read the headers of 
    .fits.Z and .fits.gz files, but for now only .fits is supported'''
#    if extn=='.fits.Z' or extn=='.fits.gz':
#        import fitsio
#        reader=fitsio.read_header
#        hierarch='' # fitsio removes the HIERARCH part of the header keyword
#    else:
#        import astropy.io.fits as pyfits
#        reader=pyfits.getheader

    all_info=[]
    # Loop through files
    for fname in files:
        
        this_info={}
        # Load the file
        hdr = pf.getheader(fname,ext=0)
        try:
            hdr2 = pf.getheader(fname,ext=1)
            if hdr2['EXTNAME'].strip() != 'MVCO': # We don't want the AO covariance matrix header...
                hdr = hdr2+hdr
        except:
            pass
        
        # Loop through the header keys
        for name in keys_dict.keys():
            try:
                result=hdr[keys_dict[name]]
            except:
                result=0
            if result=='':
                result=0
            this_info[name]=result
        
        all_info.append(this_info)

    # Fix a few potential entries
    for ix,info in enumerate(all_info):
        # Datetime > date + time
        try:
            obsdate,obstime=info['datetime'].split('T')
            info['obsdate']=str(obsdate)
            info['obstime']=str(obstime)
        except:
            pass
        info['ix']=ix
        info['fname']=files[ix].replace(wdir,'')
         # Filter
        try:
            if info['filter1'].lower().startswith('empty'):
                info['filter']=info['filter2']
            else:
                info['filter']=info['filter1']
        except:
            pass

        # Parallactic angle        
        try:
            info['parang']=(info['rotstart']+info['rotend'])/2.+info['alt']-(180.-(info['pastart']+info['paend'])/2.)
        except:
            info['parang']=0.
            pass
        
    return all_info

#####################

#####################      

def make_header_file(wdir,prefix='NACO',extn='.fits',save_name='header.txt',header_keys=target_keys):
    ''' This program makes a file summarizing the information contained in the
    headers of each cube.
    wdir = working directory
    prefix = filename prefix
    extn = filename extension
    We want:
        file name, target name, time, RA, Dec, type, filter, AX1, NDIT, DIT,
        PosAng, Camera
    '''
    
    # Find what files are in the directory
    data_files=sorted(glob.glob(wdir+prefix+'*'+extn))
    
    nfiles=len(data_files)
    print '  Files found:',nfiles

    # Use get_info_from_files to read through all of the headers and get the right values    
    all_info=get_info_from_files(data_files,header_keys,wdir=wdir)
    
    line0='IX       TARGET          Date          Time         RA        Dec    ObsType   Filter   ND_Filter  AX1   AX3    DIT  NDIT  ParAng Camera NEXPO     File Name\n'
    
    # Now print everything to the file
    with open(save_name,'w') as myf:
        
        line_template ="{ix:4}  {tname:17}  {obsdate:10}  {obstime:10}  {RA:15.9}  "
        line_template+="{Dec:15.9}  {obstype:8}  {filter:8}  {nd:7}  {ax1:4}  "
        line_template+="{ax3:4}  {dit:5.2}  {ndit:3}  {parang:8.3}  {camera:3}  "
        line_template+="{nexpo:2}  {fname:37}\n"
    
        myf.write(line0)
        for info in all_info:
            printline=line_template.format(**info)
            myf.write(printline)

#####################

#####################

def make_calib_header_file(wdir,prefix='NACO',extn='.fits',save_name='header.txt'):
    ''' This program makes a file summarizing the information contained in the
    headers of each cube.
    wdir = working directory
    prefix = filename prefix
    extn = filename extension
    We want:
        file name, target name, time, RA, Dec, type, filter, AX1, NDIT, DIT,
        PosAng, Camera
    '''
    
    # Find what files are in the directory
    data_files=sorted(glob.glob(wdir+prefix+'*'+extn))
    
    nfiles=len(data_files)
    print '  Files found:',nfiles
    
    all_info=get_info_from_files(data_files,target_keys,wdir=wdir)

    line0='IX       Date          Time         ObsType             Filter   ND_Filter  AX1   AX3  DIT  NDIT  Cam NEXPO     File Name\n'

    # Now print everything to the file
    with open(save_name,'w') as myf:
    
        myf.write(line0)

        line_template="{ix:4d}  {obsdate:10s}  {obstime:10s}  {obstype:20s}  {filter:8s}  "
        line_template+="{nd:7s}  {ax1:4d}  {ax3:4d}  {dit:5.2f}  {ndit:3d}  "
        line_template+="{camera:3s}  {nexpo:2d}  {fname:37}\n"        
        
        for info in all_info:
            printline=line_template.format(**info)
            myf.write(printline)