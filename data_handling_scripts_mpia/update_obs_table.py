# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 14:18:56 2016

This will update the observations table when run on the NACO DRS

@author: cheetham
"""

import naco_ispy

data_folder = '/data/beegfs/astro-storage/groups/henning/rlau/naco-ispy-shared/'
db_filename = '/data/beegfs/astro-storage/groups/henning/rlau/naco-ispy-shared/obs_table.dat'

# Make the object 
obs_db = naco_ispy.databases.obs_table(filename=db_filename, data_folder=data_folder)

# Run the search function to find the relevant files and update the table
# Update previous entries if possible
obs_db.search(dry_run=False, fast_update = True)

# Save the updated table
obs_db.save()