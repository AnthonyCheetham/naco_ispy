# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 15:29:56 2016

This will display the observations table in your browser.
The paths are set for my laptop.

@author: cheetham
"""

import naco_ispy

#Load it
obs_table=naco_ispy.databases.obs_table(filename='/Users/cheetham/data/data_archive/GTO/obs_table.dat')

# display it
obs_table.data.show_in_browser(jsviewer=True)