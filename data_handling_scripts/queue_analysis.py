# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:41:32 2016

Program to run through the target database and queue the unprocessed flux
data for analysis

@author: cheetham
"""

import naco_ispy,subprocess,os,argparse

parser = argparse.ArgumentParser(description='This program queues up all unprocessed NACO ISPY data for analysis.')
parser.add_argument('-flux', dest="flux",action='store_const',const=True,
                    default=False, help='Process flux data')
parser.add_argument('-targ', dest="targ",action='store_const',const=True,
                    default=False, help='Process target (agpm or saturated psf) data')    
parser.add_argument('-pca', dest="pca",action='store_const',const=True,
                    default=False, help='Run PCA on the datasets')    
parser.add_argument('-dry_run', dest="dry_run",action='store_const',const=True,
                    default=False, help='Dont actually queue the analysis, but print the commands it will do')    
parser.add_argument('--num', action="store",  dest="num", type=int, default=None,
                    help='Maximum number of datasets to process')
parser.add_argument('--skip',action="store",  dest="skip", type=int, default=0,
                    help='Skip this many datasets before starting the queing')
parser.add_argument('-reprocess',action="store_const",  dest="reprocess", const=True, default=False,
                    help='Ignore whether datasets have already been processed, and launch the reduction anyway')

# Get the input arguments
args = parser.parse_args()
num = args.num
skip = args.skip
reprocess = args.reprocess

data_folder = '/data/NACO/'
db_filename = '/data/NACO/obs_table.dat'
#data_folder='/Users/cheetham/data/naco_data/GTO/'
#db_filename='/Users/cheetham/data/data_archive/GTO/obs_table.dat'
        
dry_run = args.dry_run

# First, load the target database
obs_db = naco_ispy.databases.obs_table(filename=db_filename, data_folder=data_folder)

# Loop through the targets in the database
for targ_ix,targ_row in enumerate(obs_db.data[skip:num]):

    # Check if the dataset is consistent
    consistent = targ_row['ConsistentSequence']

    # Check which type of processing we want to launch
    if args.flux:
        process_script = '/home/spectro/cheetham/code/naco_ispy/processing_scripts/flux.slurm'
        processed = targ_row['FluxProcessed']
    elif args.targ:
        processed = targ_row['TargProcessed']
        if str(targ_row['AGPM']) == 'True':
            process_script = '/home/spectro/cheetham/code/naco_ispy/processing_scripts/agpm_graphic.slurm'
        else:
            process_script = '/home/spectro/cheetham/code/naco_ispy/processing_scripts/saturated_psf_graphic.slurm'
    elif args.pca:
        # Check that it is ready
        targ_processed = targ_row['TargProcessed']
        flux_processed = targ_row['FluxProcessed']
        
        if str(targ_processed) == 'True' and str(flux_processed) == 'True':
            
            # Check if it has already been processed    
            processed = targ_row['ADIProcessed']

            # Also, if the targ and flux have already been processed, we dont need to worry about consistency
            consistent = 'True'
            
            if str(targ_row['AGPM']) == 'True':
                process_script = '/home/spectro/cheetham/code/naco_ispy/processing_scripts/pca_graphic_agpm.slurm'
            else:
                process_script = '/home/spectro/cheetham/code/naco_ispy/processing_scripts/pca_graphic_satpsf.slurm'
        else:
            consistent = 'False' # Just to make sure it doesnt get processed
            continue
    
    # Overwrite the processed status if we want to reprocess
    if reprocess:
        processed = False
    
    # The command to run:
    cmd = "echo 'bash "+process_script+"' | at -q b now"
    if (str(processed) == 'False') and (str(consistent) == 'True'):
        
        # Change to the right directory
        os.chdir(targ_row['Location'])
        
        if dry_run:
            print('Queueing analysis for '+targ_row['Location'])
            print('  '+cmd)
        else:
            # Execute the processing command
            subprocess.call(cmd,shell=True)
            
           