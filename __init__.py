# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 10:36:26 2015

Automated program that searches directories for new data, determines the 
appropriate way to process them and schedules their analysis.

Also hosts a calibration and an observation database that keeps track of what
has already been processed and what calibrations are available

@author: cheetham
"""

import numpy as np
import scipy
import monitor
import databases
import organise_data
import header