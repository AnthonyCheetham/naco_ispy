# naco-ispy
Python tools for the NACO-ISPY survey

This repository contains two main functions:

1. The automated reduction pipeline. This consists of two submodules: a data organiser and a database tool.

2. The real-time analysis code. This is an infinite loop that you can run on the offline machine at Paranal, 
  and it will plot important information about the data as it arrives.

## Real-time analysis code:

  This is intended to be run on the offline machine at Paranal, and monitors a folder for incoming data.
  You should exit and relaunch the program when you change stars, so that the plots are all reset.
  Unfortunately it will not plot anything until it finds a sky frame, since it can't find the peak or estimate the background level without a sky. For non-AGPM frames, it will 
wait until it has at least 2 frames.
  It will ignore flux frames if their exposure times are NOT between 0.1-0.5s 

To run it, put monitor.py somewhere on the NACO offline machine, then run
	> python
    > import monitor
    > monitor.run_and_process(folder='/path/to/the/data/')

###Advanced Options:
*   prefix (default='NACO') + suffix (default='.fits'): To find the files, it looks for prefix+'*'+suffix. So if your files are named differently, change these
*   pause_interval (default=2): The number of seconds to wait between updates.
*   crop_size (default=500): The number of pixels to consider. 1024x1024 images are slow to process, so this will cut the images down (around the centre) before processing.
*   new_only (default=True): This switch allows you to look at all files in the directory, or just those that are added after the program is started.