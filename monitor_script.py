# -*- coding: utf-8 -*-
"""
This script makes running the real-time monitor as simple as possible
"""

import monitor,datetime

current_time=datetime.datetime.today()

# What was the date at the beginning of the night?
datestr='{0:4d}-{1:02d}-{2:02d}' # yyyy-mm-dd
if current_time.hour <12: # So midday in Chile is where the date swaps.
    # it is after midnight but before midday so take away a day
    date=datestr.format(current_time.year,current_time.month,current_time.day-1)
else:
    # it is after midday so the date is correct
    date=datestr.format(current_time.year,current_time.month,current_time.day)

# Where is the data?
folder='/data-ut1/raw/'+date+'/'
# Run the monitor
monitor.run_and_process(folder=folder)