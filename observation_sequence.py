# -*- coding: utf-8 -*-
""" 
This defines the observation_sequence object, which is created by the pipeline.


"""

class observation_sequence(object):
    """ 
    This is the basic class for the pipeline, and holds all of the information
    related to a single dataset for a target.
    Basic properties include lists of all the related files and their locations.
    """
    
    def __init__(self,target_name):
        
        self.target_name=target_name
        
        
    def find_calibrations(self,calib_database):
        """ 
        Queries the calibration database to find relevant calibration files
        (flat fields and astrometric calibrations)
        """
        pass

    def make_analysis_script(self):
        """ 
        Saves an analysis script that will run the pipeline on this dataset"""
        
    def update_database(self,obs_database):
        """ 
        Updates the observations database with the properties of this observation sequence.
        
        Input: 
            - database: a database object
        """
        pass
    