#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 11:41:33 2018

Script to remove intermediate products from directories 

@author: cheetham
"""

import naco_ispy,subprocess,os,argparse

parser = argparse.ArgumentParser(description='This program queues up all unprocessed NACO ISPY data for analysis.')
parser.add_argument('-flux', dest="flux",action='store_const',const=True,
                    default=False, help='Remove flux products')
parser.add_argument('-targ', dest="targ",action='store_const',const=True,
                    default=False, help='Remove targ products')    
#parser.add_argument('-pca', dest="pca",action='store_const',const=True,
#                    default=False, help='Remove PCA products')    
parser.add_argument('-dry_run', dest="dry_run",action='store_const',const=True,
                    default=False, help='Dont actually queue the analysis, but print the commands it will do')    
parser.add_argument('--num', action="store",  dest="num", type=int, default=None,
                    help='Maximum number of datasets to process')
parser.add_argument('--skip',action="store",  dest="skip", type=int, default=0,
                    help='Skip this many datasets before starting the queing')

# Get the input arguments
args = parser.parse_args()
num = args.num
skip = args.skip

data_folder = '/data/NACO/'
db_filename = '/data/NACO/obs_table.dat'
#data_folder='/Users/cheetham/data/naco_data/GTO/'
#db_filename='/Users/cheetham/data/naco_data/GTO/obs_table.dat'

scripts_directory = os.path.expanduser('~/code/naco_ispy/processing_scripts/')
        
dry_run = args.dry_run

# First, load the target database
obs_db = naco_ispy.databases.obs_table(filename=db_filename, data_folder=data_folder)

# Loop through the targets in the database
for targ_ix,targ_row in enumerate(obs_db.data[skip:num]):
    
    main_dir = targ_row['Location'] 
    os.chdir(main_dir)
    cmd = 'clean_graphic' # script to run to clean the directories
#    cmd = 'test_script'
    
    # Now go through each possible subdirectory to clean
    if dry_run:
        print('Cleaning directory:'+targ_row['Location'])
    else:
        
        if args.flux:
            os.chdir('Flux')
            subprocess.call(cmd,shell=True)
            os.chdir(main_dir)
        
        if args.targ:
            try:
                os.chdir('Sky')
                subprocess.call(cmd,shell=True)
                os.chdir(main_dir)
            except:
                print('Sky directory not accessible')
            
            # Do Sky as well if the target has AGPM
            if str(targ_row['AGPM']) == 'True':
                os.chdir('Sky')
                subprocess.call(cmd,shell=True)
                os.chdir(main_dir)
                
        # Put PCA here:
#        if
