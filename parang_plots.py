# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:49:38 2015

Some parallactic angle plots (useful for deciding on which targets to observe and when)
I expect that the OPS will make this obsolete, but I need something for the first run

@author: cheetham
"""
from numpy import sin, cos, tan, arctan2, pi, deg2rad, rad2deg
import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.time import Time
import astropy.time
import astropy.units as u
from astroquery.simbad import Simbad

# Some params about the target (can get some through query simbad in astropy)
ra_deg = 125.71 #deg
dec_deg = 1.86 # deg

date_start='2015-12-18'
tstart_utc=8.00+00./60. # in hrs
total_time=1.0+00./60. # in hrs

# plot params
deltat=30.# seconds
###########

# Some params about Paranal
geolat_deg=-24.6270
geolat_rad=geolat_deg*np.pi/180.
geolon_rad=-70.4040*np.pi/180.
location=(geolon_rad*u.rad,geolat_rad*u.rad)

#Make Time objects out of the times
t_hr=np.int(tstart_utc)
t_min=np.int((tstart_utc-t_hr)*60.)
t_sec=np.int((tstart_utc-t_hr-t_min/60.)*60.*60.)
t_start=date_start+'T'+str(t_hr).zfill(2)+':'+str(t_min).zfill(2)+':'+str(t_sec).zfill(2)

obs_start=Time(t_start,format='isot',scale='utc',location=location) # UTC

# explicitly set the number of leap-seconds to zero (since astropy wont predict the future)
obs_start.delta_ut1_utc = 0.

parangs=[]
times=[]

nt=np.int(total_time*60*60./deltat)
for ix in range(nt):
    times.append(ix*deltat/60./60.)
    current_time=astropy.time.TimeDelta(ix*deltat,format='sec')+obs_start
    
    current_time.delta_ut1_utc = 0.
    
    # convert to sidereal time
    lst=current_time.sidereal_time('apparent')
    
    if ix ==0:
        print 'Starting LST:',lst
    if ix ==nt-1:
        print 'Ending LST:',lst
        
    # convert to hour angle
    ha=lst.hourangle*15.-ra_deg # degrees
    
    # Calculate the parallactic angle
    d2r=np.pi/180.
    r2d=180./np.pi
    f1 = cos(geolat_rad) * sin(d2r*ha)
    f2 = sin(geolat_rad) * cos(d2r*dec_deg) - cos(geolat_rad) * sin(d2r*dec_deg) * cos(d2r*ha)
    pa = -r2d*arctan2(-f1,f2)
    if dec_deg > geolat_deg:
        pa = ((pa + 360) % 360)
    
    parangs.append(pa)

delta_pa=np.round(np.max(parangs)-np.min(parangs),decimals=1)
print 'Field Rotation:',delta_pa

# Now the plot
plt.clf()
plt.plot(times,parangs)
plt.xlabel('Time since start (hrs)')
plt.ylabel('Parallactic angle')
plt.title(u'$\Delta$PA='+str(delta_pa))
plt.tight_layout()