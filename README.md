# naco-ispy
Python tools for the NACO-ISPY survey

This repository contains two main functions:

1. The automated reduction pipeline. This consists of two submodules: a data organiser and a database tool.

2. The real-time analysis code. This is an infinite loop that you can run on the offline machine at Paranal, 
  and it will plot important information about the data as it arrives.

## Real-time analysis code:

  This is intended to be run on the offline machine at Paranal, and monitors a folder for incoming data.
  It will detect changes in the OB name and reset the plots so that you only need to run it once.
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




Copyright (C) 2016 Anthony Cheetham and the NACO-ISPY team

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.