# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:41:32 2016

Program to run through the target database and queue the unprocessed flux
data for analysis

@author: cheetham
"""

import naco_ispy,subprocess,os,argparse
import mpia_queuer

parser = argparse.ArgumentParser(description='This program queues up all unprocessed NACO ISPY data for analysis.')
parser.add_argument('-flux', dest="flux",action='store_const',const=True,
                    default=False, help='Process flux data')
parser.add_argument('-targ', dest="targ",action='store_const',const=True,
                    default=False, help='Process target (agpm or saturated psf) data')    
parser.add_argument('-pca', dest="pca",action='store_const',const=True,
                    default=False, help='Run PCA on the datasets')    
parser.add_argument('-dicpm', dest="dicpm",action='store_const',const=True,
                    default=False, help='Run DICPM on the datasets')    
parser.add_argument('-dry_run', dest="dry_run",action='store_const',const=True,
                    default=False, help='Dont actually queue the analysis, but print the commands it will do')    
parser.add_argument('--num', action="store",  dest="num", type=int, default=None,
                    help='Maximum number of datasets to process')
parser.add_argument('--skip',action="store",  dest="skip", type=int, default=0,
                    help='Skip this many datasets before starting the queing')
parser.add_argument('-reprocess',action="store_const",  dest="reprocess", const=True, default=False,
                    help='Ignore whether datasets have already been processed, and launch the reduction anyway')

parser.add_argument('--queue_script', dest="queue_script",action='store',type=str,default=None,
                    help='Queue this script instead of one of the default options.')


# Get the input arguments
args = parser.parse_args()
num = args.num
skip = args.skip
reprocess = args.reprocess
queue_script = args.queue_script

data_folder = '/data/beegfs/astro-storage/groups/henning/cheetham/naco-ispy-shared/'
db_filename = '/data/beegfs/astro-storage/groups/henning/cheetham/naco-ispy-shared/obs_table.dat'
#data_folder='/Users/cheetham/data/naco_data/GTO/'
#db_filename='/Users/cheetham/data/naco_data/GTO/obs_table.dat'

scripts_directory = os.path.expanduser('~/code/naco_ispy/processing_scripts/')
        
dry_run = args.dry_run

# First, load the target database
obs_db = naco_ispy.databases.obs_table(filename=db_filename, data_folder=data_folder)

# Loop through the targets in the database
for targ_ix,targ_row in enumerate(obs_db.data[skip:num]):

    # Check if the dataset is consistent
    consistent = targ_row['ConsistentSequence']

    # Check which type of processing we want to launch
    if args.flux:
        process_script = scripts_directory+'flux.slurm'
        processed = targ_row['FluxProcessed']
    elif args.targ:
        processed = targ_row['TargProcessed']
        if str(targ_row['AGPM']) == 'True':
            process_script = scripts_directory+'agpm_graphic.slurm'
        else:
            process_script = scripts_directory+'saturated_psf_graphic.slurm'
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
                process_script = scripts_directory+'pca_graphic_agpm.slurm'
            else:
                process_script = scripts_directory+'pca_graphic_satpsf.slurm'
        else:
            consistent = 'False' # Just to make sure it doesnt get processed
            continue

    elif args.dicpm:
        # Process using Matthias' DICPM pipeline
        # Check that it is ready
        targ_processed = targ_row['TargProcessed']
        flux_processed = targ_row['FluxProcessed']
        
        if str(targ_processed) == 'True' and str(flux_processed) == 'True':
            
            # If the targ and flux have already been processed, we dont need to worry about consistency
            consistent = 'True'
            
            process_script = scripts_directory+'dicpm_reduction.slurm'
            processed = 'False'
        else:
            consistent = 'False' # Just to make sure it doesnt get processed
            processed = 'True' # Just to make sure it doesnt get processed
            continue

    elif args.queue_script:
        # Queue the input script instead
        ###########
        # This code should be temporary, and can be used to run once-off scripts
        process_script = queue_script
        processed=False
        
        # Currently set up for contrast estimation, so need the data to be ADIProcessed
        if str(targ_row['ADIProcessed']) != 'True':
            continue
        
        ###########
        
    # Overwrite the processed status if we want to reprocess
    if reprocess:
        processed = False
    
    # The command to run:
    if (str(processed) == 'False') and (str(consistent) == 'True'):
        
        job = [targ_row['Location'],process_script]
        
        if dry_run:
            print('Queueing analysis for '+targ_row['Location'])
            print('  '+process_script)
        else:
            # Save it to the queue file since we're using our own system
            # Get the existing list, append the new job and rewrite it
            all_jobs = mpia_queuer.read_file(mpia_queuer.queue_file)
            all_jobs.append(job)
            mpia_queuer.write_file(mpia_queuer.queue_file,all_jobs)
            
            
           