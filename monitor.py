# -*- coding: utf-8 -*-
"""
NACO AGPM/Saturated PSF real-time statistics module.

  The main program here is run_and_process. Everything else is defined so the 
code is easier to follow. It is intended to be run on the offline machine at 
Paranal, and monitors a folder for incoming data.
  You should exit and relaunch the program when you change stars, so that the plots
are all reset.
  Unfortunately it will idle until it finds a sky frame, since it can't find the 
peak or estimate the background level without a sky. For non-AGPM frames, it will 
wait until it has at least 2 frames.
  It will ignore flux frames if their exposure times are NOT between 0.1-0.5s 

Known bugs:
- If you try to exit while matplotlib is thinking, it won't exit properly and you may 
  have to close the terminal window. The only way to fix this is to change the plots 
  to an interactive GUI, which I might do later. For now, I've intentionally added 
  a 2s pause to reduce the chances of this happening.
  
  An example of how to run it from the directory you put the code:
    import monitor
    monitor.run_and_process(folder='/path/to/the/data/')
  
"""

#Ideas:
#    - Just plot last ~10 cubes? No, it is more useful to show everything. Can always
#        rerun the program for every new target.
#    - Organise data by target and then plot all data for a certain target?
#    - Plot the standard deviation of the flux in the donut (for the agpm)?
#    - Plot the agpm centering?
#
#Problems:
#  - Peak stellar flux + background doesnt work with dithered data.
#  - Infinite loops don't play nicely with matplotlib. Sometimes ctrl+c doesn't work.
#      Very hard to reproduce, but might be related to exception generated during plt.pause
#  
#Solved Problems:
# - SOLVED: Need to do sky subtraction to check peak flux, since the centre of the agpm 
#     can be saturated while the actual light that we care about is not.
# - SOLVED: Infinite loops don't work with the default MacOSX backend for matplotlib.
#     Have to use (e.g.) ipython --matplotlib='tk'


# Before we start, change the backend to Qt4Agg, since the MacOSX default doesnt work.
#  The combination of the infinite loop and the plt.tight_layout() call (as well as the 
#  plt.show() and plt.pause() calls) causes problems with the macosx backend
import matplotlib as mpl
mpl.use('TkAgg') # if this doesn't work, try the next line instead
#mpl.use('QT4Agg')

import numpy as np

import matplotlib.pyplot as plt
import glob,time,datetime
import astropy.io.fits as pyfits
from astropy.time import Time
from matplotlib.dates import DateFormatter
import scipy.signal as signal
#from Tkinter import *
#import tkMessageBox
#import pdb

plt.interactive(True)

# Here are all of the hard-coded numbers in case we need to change any
#nonlinear_limit=18000.
#saturate_limit=22000.
#minimum_flux=-7000. # The actual "zero" point
nonlinear_limit=11000.
saturate_limit=15000.
minimum_flux=0. # This should reduce some confusion

nexpo_limit=2 # If nexpo > 2, this indicates that it is a target observation. otherwise, sky
obstime_limits=[0.1,0.5] # all target/sky observations have exp times in this range. Anything outside is a flux frame.

smooth_dist=4 # FWHM of gaussian used to smooth the images before measuring the 
              # background and finding the centre

def detect_filetype(hdr,get_folder_string=False):
    ''' Works out what kind of file it is based on the header.
    This function is necessarily complicated and hard to read, since there
    are so many cases it has to cover.'''
    
    type_flag=hdr['HIERARCH ESO DPR TYPE']
#    type_cat = hdr['HIERARCH ESO DPR CATG']
    expt=hdr['EXPTIME'] # exposure time.
    agpm=hdr['HIERARCH ESO INS OPTI1 ID'] # this is AGPM if it is used
    try:
        targ_name=hdr['HIERARCH ESO OBS NAME']
    except:
        targ_name='NoName'
    naxis=hdr['NAXIS']
    naxis1=hdr['NAXIS1']

    try:
        nexpo=hdr['HIERARCH ESO SEQ NEXPO']
    except:
        nexpo=0

    # Now format all of these strings
    if 'astcal' in targ_name.lower():
    # Astrometric calibrators are an annoying case that we have to deal with first
    # For now, assume they have "AstCal" in their target names

        obstype='AstCal'
        folder_string='AstCal'
    
    elif type_flag=='OBJECT':
        # We need to work out which of the "OBJECT" frames are skies, flux 
        #  frames and actual target observations.
        #
        # For the AGPM skies, we can use the number of exposures. Skies are a single cube for the AGPM.
        # There are no separate skies for non-AGPM observations, so label them all as Targ.
        # 
        # For the flux, the only way to guess is the exposure time (or possibly the ND?)
    
        # Handle the AGPM and non-AGPM cases differently
        if agpm=='AGPM':
            if nexpo > nexpo_limit:
                obstype='Target_AGPM'
                folder_string='Targ'
            elif (expt < obstime_limits[1]) and (expt > obstime_limits[0]) and (naxis1 >512):
                obstype='Sky'
                folder_string='Sky'
            else:
                obstype='Flux'
                folder_string='Flux'
        else:
            if ((expt < obstime_limits[1]) and (expt > obstime_limits[0])):
                obstype='Target_saturated'
                folder_string='Targ'
            else:
                obstype='Flux'
                folder_string='Flux'
                
    elif type_flag=='SKY':
        obstype='Sky'
        folder_string='Sky'
    elif 'FLAT' in type_flag:
        obstype='Flat'
        folder_string='Flats'
    elif type_flag.lower=='psf-calibrator':
        obstype='Flux'
        folder_string='Flux'
    
    # We don't actually use any of the following types, but I thought we might as well
    # put them somewhere
    elif type_flag=='STD':
        obstype='Std'
        folder_string='STD'
    elif 'DARK' in type_flag:
        obstype='Dark'
        folder_string='Dark'
    else:
        # Put all of the unknown file types into a single folder to make it easy
        print 'Unrecognised DPR type:',type_flag
        obstype='Unknown'
        folder_string='Uncategorized'
        
    # But if it has NAXIS3=0, it is really an acquisition!
    if naxis==2 and (obstype != 'Flat') and (obstype !='Dark'):
        folder_string='Acq_'+folder_string
        obstype='Acq'
        
    
    if get_folder_string:
        return obstype,folder_string
    else:
        return obstype
    


###################

###################

def diagnostic_plots(axes,capture_time,peakcounts,bgflux,parangs,clean_im):
    ''' Wrapper function for the diagnostic plots for the real-time monitor'''

    # Clear the plots
    ax1,ax2,ax3,ax4=axes.flatten()
    ax1.cla()
    ax2.cla()
    ax3.cla()
    ax4.cla()

    # Work out the order of the data, just in case it is not in chronological order    
    order=np.argsort(capture_time)
    
    t_lims=[np.min(capture_time),np.max(capture_time)]
    
    # Plot 1: The peak flux
    ax1.cla()
    ax1.plot_date(capture_time[order],peakcounts[order],'x',label='Peak flux')
    ax1.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax1.set_title('Peak Stellar Flux (or peak around agpm donut)')
    ax1.set_xlabel('Time')
    ax1.set_ylim(np.min([0,np.min(bgflux)]),1.2*np.max([np.max(bgflux),saturate_limit])) # force the plot to start at zero so it is easier to read
    # plot the nonlinear and saturation regimes
    ax1.plot(t_lims,[nonlinear_limit,nonlinear_limit],'r')
    ax1.plot(t_lims,[saturate_limit,saturate_limit],'k')
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)
    
    # Plot 2: Background flux
    ax2.cla()
    ax2.plot_date(capture_time[order],bgflux[order],'x',label='Background flux')
    ax2.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax2.set_title('Background Flux')
    ax2.set_xlabel('Time')
    ax2.ticklabel_format(axis='y',useOffset=False)
#    ax2.set_ylim(np.min([0,np.min(bgflux)]),1.2*np.max([np.max(bgflux),saturate_limit])) # force the plot to start at zero so it is easier to read
    # plot the nonlinear and saturation regimes
#    ax2.plot(t_lims,[nonlinear_limit,nonlinear_limit],'r')
#    ax2.plot(t_lims,[saturate_limit,saturate_limit],'k')
    for tick in ax2.get_xticklabels():
        tick.set_rotation(45)
    
    # Plot 3: Parallactic angle
    ax3.cla()
    ax3.plot_date(capture_time[order],parangs[order],label='Parallactic angle')
    ax3.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax3.set_title('Parallactic Angle')
    ax3.set_xlabel('Time')
    for tick in ax3.get_xticklabels():
        tick.set_rotation(45)
    
    # plot 4: FWHM of image... Need to fit these first
    # For now, just plot the image (should be x and y position in the future)
    ax4.cla()
    try:
        ax4.imshow(clean_im)
        ax4.colorbar()
    except:
        pass
    ax4.set_title('Clean image')
#    ax4.set_title('Clean image (will be psf width in the future)')

###################

###################

def quick_clean(im,sky,crop_size):
    ''' Does some quick data cosmetics so it can be used in the real-time analysis plots'''
    
    image_size=np.min([im.shape[0],crop_size])
    # crop the image so we don't have to deal with the region outside the agpm
    im=im[im.shape[0]/2-image_size/2:im.shape[0]/2+image_size/2,
          im.shape[1]/2-image_size/2:im.shape[1]/2+image_size/2]
             
    # change it so that the zero point is actually zero
    im-=minimum_flux

    # sky subtract and return
    return im-sky

###################

###################

def check_data(head,window=None):
    ''' This program is a place to put any warnings that the data are bad
        First use is for the Full_Uszd mask, since we never want to use it
        but some datasets seem to have it even though the OB didn't ask for it
    '''
    warning = False
    
    if head['ESO INS OPTI3 NAME'] == 'Full_Uszd':

        # If we want to make a text box pop up, use this code:
        #centre screen message
#        window.geometry("1x1+"+str(window.winfo_screenwidth()/2)+"+"+str(window.winfo_screenheight()/2))
#        tkMessageBox.showwarning(title="NACO-ISPY monitor", message="WARNING! Full_Uszd mask is inserted. Ask the night astronomer to send a new preset")

        # If we want to print the message to the console, use this code:
        print('WARNING: Full_Uszd mask is inserted. Ask the night astronomer to send a new preset.')
        warning = True
        
    return warning
    

###################

###################


def run_and_process(folder='./',prefix='NACO',suffix='.fits',
                    pause_interval=2.,crop_size=500,new_only=True):
    '''
    Run this program on a directory of NACO data to display some important
    information in real-time as data is added. Currently this program plots the
    peak flux, background flux and parallactic angle of each frame as a function 
    of time. 
    Options:
        - folder: the path to the folder to monitor
        - prefix: the name of the naco data (e.g. NACO_CORO_SCI)
        - suffix: the file type (e.g. .fits)
        - pause_interval: the delay between updates for the plots
        - crop_size: the number of pixels to consider. This is used to crop 
            the image and speed up the processing. This is taken from the _centre_
            so be careful not to crop the star out!
        - new_only: if True, it will ignore files that already exist in a folder
             and only display the statistics of new files.
        '''
    print 'Monitoring folder:',folder
    print 'Press ctrl+c to exit'
    
    # Make sure the / is included in the filename
    if folder[-1]!='/':
        folder=folder+'/'
        
    # Set up all of the arrays that will hold the information
    known_files=np.array([])
    capture_time=np.array([])
    peakcounts=np.array([])
    bgflux=np.array([])
    parangs=np.array([])
    target_names=np.array([])
    
    # Set up the plots
    fig,axes=plt.subplots(2,2,num=0)
    
    fig.canvas.set_window_title('Summary of data')
    # Set up some arrays:
    skies={} # a dictionary to contain all of the skies
    clean_im=0

    # Set up error window    
#    window = Tk()
#    window.wm_withdraw()
    window = None

    # Begin the main loop
    repeats=0
    first_plot=True
    while True:
        try:
            # Find the data in the folder
            files=glob.glob(folder+prefix+'*'+suffix)
            
            # Remove any acquisition files
            acq_files = glob.glob(folder+'*ACQ*'+suffix)
            files = set(files)
            acq_files = set(acq_files)
            files = list(files-acq_files)
            
            
            nfiles=len(files)
            
            if new_only and repeats==0:
                known_files=files
        
            # Now find which files are new
            new_files=list(set(files)-set(known_files))
            n_new=len(new_files)
            if nfiles ==0 and repeats ==0:
                print 'No files found'
                time.sleep(pause_interval)
            elif n_new >0:
                pass

            # Sort them so they are in filename order 
            #  (which should also correspond to the time they were made)
            new_files=sorted(new_files)
                
            # Go through them and see what to do with them
            for f in new_files:
                head=pyfits.getheader(f)
                exptime=np.round(head['EXPTIME'],decimals=3) # to the nearest ms to avoid mismatches
                
                # Classify the file (we only care about sky and target for now)
                obstype=detect_filetype(head)
                
                # Print any warnings about the data or observing strategy
                warnings = check_data(head,window)
                
                
                # If it is a saturated psf, we can make a dodgy sky by combining all of the data
                if obstype=='Target_saturated':
                    # Work out if we already have a sky for this target
                    if skies.has_key(str(exptime)):
                        sky=skies[str(exptime)]
                        nsky=skies['n_'+str(exptime)]                        
                    else:
                        sky=0
                        nsky=0
                        skies['last4_'+str(exptime)]=[]
                    
                    im=pyfits.getdata(f)[0]                    
                    this_sky=quick_clean(im,0,crop_size)
                    
                    # Update the existing sky estimate
#                    sky=(nsky*sky+this_sky)/(nsky+1) # the old way that has self-subtraction
                    skies['last4_'+str(exptime)].append(this_sky) # track the last 4                    
                    skies['n_'+str(exptime)]=nsky+1
                    
                    # if we have more than 4, pop the first one and continue
                    if (nsky+1) >4:
                        skies['last4_'+str(exptime)].pop(0)
                        skies['n_'+str(exptime)]=nsky
                    skies[str(exptime)]=np.median(skies['last4_'+str(exptime)],axis=0)
                
                if obstype=='Sky':
                    # If it is a sky, update the master sky frame (for that exposure time)
                    im=pyfits.getdata(f)[0]
                    this_sky=quick_clean(im,0,crop_size)
                    
                    skies[str(exptime)]=this_sky
                    
                if obstype=='Target_AGPM' or obstype=='Target_saturated':                    
                    
                    # sky subtract
                    if skies.has_key(str(exptime)):
                        sky=skies[str(exptime)]
                    else:
                        # if the sky doesnt exist yet, skip this file for now 
                        #  and come back to it
                        files.remove(f)
                        continue

                    # We don't want to sky subtract the first frame with itself...
                    if obstype=='Target_saturated' and skies['n_'+str(exptime)]==1:
                        files.remove(f)
                        # To avoid problems with the case of only 1 file (where it 
                        #  make a sky from 2 copies of itself), reset the sky until 
                        #  we have another file
                        if len(files)==1:
#                            print 'deleting sky'
                            skies.pop(str(exptime))
                            skies.pop('n_'+str(exptime))

                        continue
                    
                    im=pyfits.getdata(f)[0]                
                    im=quick_clean(im,0,crop_size)
                    
                    clean_im=im-sky
                    
                    #  measure the background level
                    bg=np.median(sky)
                    bgflux=np.append(bgflux,bg)
                
                    # Save the observing time
                    t=Time(head['MJD-OBS'],format='mjd')
                    capture_time=np.append(capture_time,t.datetime)
                
                    # Measure the peak flux
                    # Pixel distance map
                    npix=im.shape[1]
                    xarr=np.arange(0,npix)-npix/2
                    xx,yy=np.meshgrid(xarr,xarr)
                    pix_dist_map=np.sqrt(xx**2+yy**2)
                    # Smooth the image for centering
                    circ_ap=np.zeros((npix,npix))
                    circ_ap[pix_dist_map<(smooth_dist/2)]=1
                    convol_sz=np.int(np.ceil(smooth_dist)+3)
                    circ_ap=circ_ap[npix/2-convol_sz/2:npix/2+convol_sz/2,
                                    npix/2-convol_sz/2:npix/2+convol_sz/2]
                    smooth_image=signal.fftconvolve(clean_im,circ_ap,mode='same')
                    
                    mx=np.where(smooth_image ==np.max(smooth_image))
                    peak_flux=im[mx[0][0],mx[1][0]]
#                    pdb.set_trace()
                    
                    peakcounts=np.append(peakcounts,peak_flux)
                    
                    # the parang (just use the rough value in the header...)
                    parang=head['HIERARCH ESO ADA POSANG']
                    parang = ((parang + 360) % 360)
                    parangs=np.append(parangs,parang)
                    
                    # Find the target name
                    target_name=head['HIERARCH ESO OBS NAME']
                    target_names=np.append(target_names,target_name)
                    last_target_name=target_name
                
                
            # Find the order that the data was taken in, by sorting the observation times
            if len(capture_time) >0:
                display_sz=80
                cropped_im=clean_im[mx[0][0]-display_sz/2:mx[0][0]+display_sz/2,
                                    mx[1][0]-display_sz/2:mx[1][0]+display_sz/2]
                                    
                # Remove all data from previous targets and only plot the current one
                target_ix=target_names==last_target_name
#                pdb.set_trace()
                
                diagnostic_plots(axes,capture_time[target_ix],peakcounts[target_ix],
                                 bgflux[target_ix],parangs[target_ix],cropped_im)
                if first_plot==True:
                    plt.tight_layout()
                    first_plot=False
            
            known_files=files

            plt.pause(0.05) # this gives python some time to make the plot
            time.sleep(pause_interval) # we cant use plt.pause because it catches 
            #             the KeyboardInterrupt and makes it hard to exit 
        except KeyboardInterrupt:
            break
        repeats+=1
   
###################

###################

# This code should run if you run it directly, e.g. python monitor.py
# It should get the correct date and monitor the correct folder on the offline 
#  machine at Paranal.
if __name__ == "__main__":
    current_time=datetime.datetime.today()

    # What was the date at the beginning of the night?
    datestr='{0:4d}-{1:02d}-{2:02d}' # yyyy-mm-dd
    if current_time.hour <12: # So midday in Chile is where the date swaps.
        # it is after midnight but before midday so take away a day
        delt=datetime.timedelta(1)
        current_time-=delt
        date=datestr.format(current_time.year,current_time.month,current_time.day)
    else:
        # it is after midday so the date is correct
        date=datestr.format(current_time.year,current_time.month,current_time.day)
    
    # Where is the data?
    folder='/data-ut1/raw/'+date+'/'

    # Run the monitor
    run_and_process(folder=folder)