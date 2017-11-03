# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:41:32 2016

Program to run through the calibrations folders and queue all data for analysis

It isn't yet smart enough to check which ones are done already

@author: cheetham
"""

import naco_ispy,subprocess,os,argparse,glob

parser = argparse.ArgumentParser(description='This program queues up all unprocessed NACO ISPY calibration data for analysis.')
parser.add_argument('-dry_run', dest="dry_run",action='store_const',const=True,
                    default=False, help='Dont actually queue the analysis, but print the commands it will do')    
parser.add_argument('--num', action="store",  dest="num", type=int, default=-1,
                    help='Maximum number of datasets to process')

# Get the input arguments
args = parser.parse_args()
num = args.num

data_folder = '/data/NACO/'
# db_filename = '/data/NACO/calib_table.dat'
# data_folder='/Users/cheetham/data/naco_data/GTO/'
#db_filename='/Users/cheetham/data/data_archive/GTO/obs_table.dat'
        
dry_run = args.dry_run

# First, load the target database
# calib_db = naco_ispy.databases.calib_table(filename=db_filename, data_folder=data_folder)

scripts_directory = os.path.expanduser('~/code/naco_ispy/processing_scripts/')

# Instead of using the database, use glob to find all folders
all_folders = glob.glob(data_folder+'Calib/*/')

# Loop through the targets in the database
for targ_ix,targ_folder in enumerate(all_folders[0:num]):

    # Check what we want to process    
    process_script = scripts_directory+'naco_calibrations.slurm'
    
    # The command to run:
    cmd = "echo 'bash "+process_script+"' | at -q b now"
        
    # Change to the right directory
    os.chdir(targ_folder)
    
    if dry_run:
        print('Queueing analysis for '+targ_folder)
        print('  '+cmd)
    else:
        # Execute the processing command
        subprocess.call(cmd,shell=True)
        
       