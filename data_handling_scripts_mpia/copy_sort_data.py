#!/usr/bin/python
# -*- coding: utf-8 -*-
# This script copies the data from the archive machine to the NACO-DRS.
# It then runs the sorting code on the data

import naco_ispy


#########
## Define the directories
#########

## THESE ARE FOR ANTHONY'S LAPTOP, used for testing
# archive_folder = '/Users/cheetham/data/data_archive/naco_archive/'
# new_data_folder = '/Users/cheetham/data/data_archive/GTO/New/'
# reduction_folder = '/Users/cheetham/data/data_archive/GTO/'
# copied_files_list = '/Users/cheetham/data/data_archive/GTO/copied_files.txt'

## The folders on astro-node4
# archive_folder = '/obs/insdata/NACO/raw/'
new_data_folder = '/data/beegfs/astro-storage/groups/henning/cheetham/naco-ispy-shared/raw_data/'
reduction_folder = '/data/beegfs/astro-storage/groups/henning/cheetham/naco-ispy-shared/'
# copied_files_list = '/data/NACO/copied_files.txt'

#########
## Copy the data from the archive to the new data directory
#########
# Don't have data download set up at the moment.
# naco_ispy.organise_data.copy_from_archive(archive_folder,new_data_folder,dry_run=False,filename_list=copied_files_list)


#########
## Sort the data by date, target and type
#########

naco_ispy.organise_data.organise_data(new_data_folder,save_dir=reduction_folder,dry_run=False)

