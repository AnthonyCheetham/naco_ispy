# -*- coding: utf-8 -*-
"""

This module contains all of the database related code for the NACO-GTO pipeline
automator

Purpose of Target Database:
    - Keep track of each target that has been observed, some relevant info and its status
    - Include the seeing, elevation, V+L mag, exp time, saturation level (?), psf width,
        parallactic angle change, a manual quality flag, field contamination, location of
        data and end products, reduction date, pipeline version(?).
    - Idea: Run at the end of the data reduction, to add the entry.
        - Need some kind of GUI to interact and edit manually...
    
Purpose of Calib Database:
    - Keep track of each calibration file (flats, true north offsets?)

TODO:
    - Define database fields properly somewhere easy to access

"""

import numpy as np
import os,glob,time,string
import astropy.io.fits as pf
from astropy.time import Time
from astropy.table import Table
import datetime
from monitor import detect_filetype
import organise_data

########################

def read_rdb(rdb_file,h=0,comment=None):
	"""
	Reads an rdb file

	Input:
	 rdb_file: filename
	 h: line number of header
	 comment: comment char, for lines to be ignored

	Output: content of the file in form of a list
	"""

	f = open(rdb_file,'r');
	data = f.readlines()
	f.close()

	# take the second line to define the list keys.
	key = string.split(data[h][:-1],'\t')
	data_list = {}
	for i in range(len(key)): data_list[key[i]] = []

	for line in data[h+2:]:
		if not line[0]==comment or line[0:2]=='--':
			qq = string.split(line[:-1],'\t')
			for i in range(len(key)):
				try: value = float(qq[i])
				except ValueError: value = qq[i]
				data_list[key[i]].append(value)

	return data_list

        
########################
        
class obs_table(object):
    ''' This class represents a locally stored astropy table which contains all
    of the relevant information about the data it finds (such as the target name,
    date, type, status, field rotation etc.
    
    If you want to change from astropy tables to a database, it should be easy
    to change the __init__ , create, save and add_entry functions. The rest 
    should work.

    filename : the name of the stored table
    table_format : the format of the stored table
    data_folder : the root directory of the data folder.
        It is assumed that the data are stored in data_folder+'Science/'
    '''
    
    ########################
    
    def __init__(self,filename='./obs_table.dat',table_format='ascii.csv',
                 data_folder='/Users/cheetham/data/naco_data/GTO/'):
        
        self.filename=filename
        self.table_format=table_format
        
        # Define the columns and column names
        self.columns=('TargetName','Band','AGPM','AnalysisStatus','ConsistentSequence','FieldRotation',
                 'r0','t0','Location','Date','Time','WindSpeed','ExpTime','SaturationLevel',
                 'PsfXWidth','PsfYWidth','Vmag','Kmag','PsfReference')
        f=np.float64
        self.dtypes=('S40','S20',np.bool_,'S40',np.bool_,f,
                f,f,'S200','S12','S14',f,f,f,
                f,f,f,f,np.bool_)
                
        if data_folder[-1] != os.sep:
            data_folder=data_folder+os.sep
        self.data_folder=data_folder+'Science'+os.sep
        
        try:
            data=Table.read(filename,format=table_format)
            self.data=data
        except:
            print 'Observations Table not found!',filename
    
    ########################
    
    def create(self):
        ''' Create a new table with all of the right columns for the database.'''

        # create the table
        data=Table(names=self.columns,dtype=self.dtypes,masked=False)
        self.data=data
        
    ########################
    def save(self):
        ''' Save changes to the table to disk'''
        print 'Saving table as',self.filename
        self.data.write(self.filename,format=self.table_format)
        
    ########################
        
    def add_entry(self,obs_sequence):
        ''' A wrapper that adds the observing sequence data to the table of observations.
        Assumes that all of the data is stored in a dictionary as obs_sequence.dict'''
        
        # Remove any previous entry based on a duplication of the Location value
        self.data.remove_rows(self.data['Location'] == obs_sequence['Location'])
        
        # Add the new row
        self.data.add_row(obs_sequence)
        
        # save it to disk (actually we can do this at the end)
#        print 'Saving table as',self.filename
#        self.data.write(self.filename,format=self.table_format)
        
    ########################

    def get_header_info(self,targ_row,head):
        ''' Looks through a header to find relevant info for the table.
        targ_row is a dictionary corresponding to the table row
        '''

        # AGPM?              
        obstype=detect_filetype(head)
        if obstype=='Target_AGPM':
            targ_row['AGPM']=True
        elif obstype=='Target_saturated':
            targ_row['AGPM']=False
        else:
            raise ValueError('The file chosen by the database is not agpm or saturated!')
            
        # Filter:
        filt1=head['HIERARCH ESO INS OPTI6 ID']
        if filt1=='L_prime':
            filt=filt1
        else:
            filt=head['HIERARCH ESO INS OPTI5 ID']
        targ_row['Band']=filt
        
        # Exposure Time
        targ_row['ExpTime']=head['HIERARCH ESO DET DIT']

        # Target Name        
        try:
            target_name=head['HIERARCH ESO OBS TARG NAME'] # this will fail for cal files
        except:
            target_name=''
        targ_row['TargetName']=target_name
        
        # And the easy one
        date=Time(head['MJD-OBS']-0.5,format='mjd').iso.split(' ')
        targ_row['Date']=date[0]
        targ_row['Time']=date[1]
        
        return targ_row
        
    ########################
        
    def get_sequence_info(self,targ_row,targ_files,read_every_n=5):
        ''' Get header information from a series of files'''
        
        # Set up the arrays first
        parangs=np.array([])
        r0=np.array([])
        t0=np.array([])
        windspeed=np.array([])

        # To save time, only loop over every n files
        # But make sure we read the first and last
        loop_files=targ_files[::read_every_n]
        loop_files.append(targ_files[-1])

        # Loop over files
        for f in loop_files:
            hdulist=pf.open(f)
            head=hdulist[0].header
            
            # Append the values to the arrays defined earlier
            try:
                parangs=np.append(parangs,(head['HIERARCH ESO ADA POSANG'] + 360) % 360)
            except:
                print('Error finding ESO ADA POSANG header in: '+f)
                parangs=0
            r0=np.append(r0,head['HIERARCH ESO AOS RTC DET DST R0MEAN'])
            t0=np.append(t0,head['HIERARCH ESO AOS RTC DET DST T0MEAN'])
            windspeed=np.append(windspeed,head['HIERARCH ESO TEL AMBI WINDSP'])
        
        # Average and save
        # We need to take care with the parangs in case they wrap.
        #   We can calculate the rotation using parangs and parangs+180 and 
        #   take the min, which works since FieldRotation<180deg
        parangs180=(parangs+180.) % 360.
        targ_row['FieldRotation']=np.min([np.max(parangs)-np.min(parangs),
                                    np.max(parangs180)-np.min(parangs180)])
        targ_row['r0']=np.median(r0)
        targ_row['t0']=np.median(t0)
        targ_row['WindSpeed']=np.median(windspeed)
        
        return targ_row
            

    ########################

    def get_rbd_info(self,targ_row,rdb_files):
        ''' Read other relevant information from the rdb files generated by
        the GRAPHIC pipeline
        '''

        psf_xwidths=np.array([])
        psf_ywidths=np.array([])
        
        # Loop over files        
        for rf in rdb_files:
            
            rdb_info=read_rdb(rf)
            
            psf_xwidths=np.append(psf_xwidths,np.median(np.abs(rdb_info['psf_fit_width_x'])))
            psf_ywidths=np.append(psf_ywidths,np.median(np.abs(rdb_info['psf_fit_width_y'])))


        targ_row['PsfXWidth']=np.median(psf_xwidths)
        targ_row['PsfYWidth']=np.median(psf_ywidths)
        return targ_row
        
    ########################
        
    def search(self,read_every_n=5,silent=False,dry_run=True):
        ''' Search self.data_folder for relevant files and update the table 
        with relevant information
        '''
        
        # Find the list of directories
        targ_directories=glob.glob(self.data_folder+'*/*/')
        
        # Big loop over directories
        for targ_dir in targ_directories:
            
            # make a new row
            targ_row={}
            
            targ_files=sorted(glob.glob(targ_dir+'Targ/NACO*.fits'))
            if len(targ_files) ==0:
                print('Warning: No files found in directory:',targ_dir+'Targ/')
                
            # Take the middle file to estimate some params          
            head=pf.getheader(targ_files[len(targ_files)/2])
            targ_row=self.get_header_info(targ_row,head)
            
            # Now get the parameters we need the whole sequence for
            targ_row=self.get_sequence_info(targ_row,targ_files,read_every_n=read_every_n)
            
            # Work out if the data has been processed based on whether any rdb files exist
            rdb_files=sorted(glob.glob(targ_dir+'Targ/cube-info/all_info_*.rdb'))
            if len(rdb_files) ==0:
                targ_row['AnalysisStatus']='NotStarted'
            else:
                targ_row['AnalysisStatus']='InProgress'
                # Get the info from the rdb files
                targ_row=self.get_rbd_info(targ_row,rdb_files)
                
            # ACTUALLY, we should get the X and Y widths from the processed psf files, if they exist.
            
            # And update the AnalysisStatus if the final_products have been added
            # TODO
            
            # Miscellaneous things
            targ_row['Location']=targ_dir
            
            # Check that the sequence is consistent (i.e. all files use the same setup)
            setup_keys=['filter','camera','nd','dit','ax1']
            flux_ok,targ_ok=organise_data.check_directory_consistency(targ_dir,setup_keys,
                      silent=silent,dry_run=dry_run)
            targ_row['ConsistentSequence']= flux_ok & targ_ok
            
            # How do I get these?
#            targ_row['SaturationLevel']=
#            targ_row['Vmag']=
#            targ_row['Kmag']=
            targ_row['PsfReference']=True # Default value. Otherwise, set it to False manually
            
            # Update the table
            self.add_entry(targ_row)
    
    
########################
        
########################
        
class calib_table(obs_table):
    ''' An alternative to the calibrations database object. This is just a locally
    stored astropy table. This is a sub-class of the obs_table class.
    Differences to the obs_table:
        - Columns are different
        - 
    '''
    
    def __init__(self,filename='./calib_table.dat',table_format='ascii.csv',
                 data_folder='/Users/cheetham/data/naco_data/GTO/'):
                
        self.filename=filename
        self.table_format=table_format
        
        # Define the columns and column names
        self.columns=('Band','Type','QualityFlag','Location',
                 'Date','NAXIS1','NAXIS2')

        self.dtypes=('S10','S10',np.bool_,'S200',Time,np.int,np.int)
        
        try:
            data=Table.read(filename,format=table_format)
            self.data=data
        except:
            print 'Calibrations Table not found!',filename
            
        
        if data_folder[-1] != os.sep:
            data_folder=data_folder+os.sep
        self.data_folder=data_folder+'Calib'+os.sep
            
    
    ########################
    
    def search(self,read_every_n=5):
        ''' Search self.data_folder for relevant files and update the table 
        with relevant information
        '''
        
        # Go through the reduced cal files one-by-one
        # Find the flats
        all_flats=sorted(glob.glob(self.data_folder+'*/*/Flat*.fits'))
        
        for flat_file in all_flats:
            
            calib_row={}
            
            # Read the header
            with pf.open(flat_file) as f:
                header=f[0].header

            calib_row['Location']=flat_file
                
            #Get the relevant params from the header            
            calib_row['Band']=header['Filter']
            calib_row['Type']=header['Flat Type']
            calib_row['QualityFlag']=header['Flat Quality']
            calib_row['Date']=header['MJD-OBS']
            calib_row['NAXIS1']=header['NAXIS1']
            calib_row['NAXIS2']=header['NAXIS2']
            
            # How do we get this?
#            calib_row['CalibID']
            
            # Update the table
            self.data.add_row(calib_row)
        
        # Find the astrometric calibrations
        print 'Astrometric calibrations not implemented'
        

# Maybe astropy tables are the best way to go?
#db=obs_table(data_folder='/Users/cheetham/data/data_archive/GTO/')
#db.create()
#######
#db.search(silent=True)
##db.search(silent=True)
#db.data.show_in_browser(jsviewer=True)
#db.save()

#cal_db=calib_table()
#cal_db.create()
##
#cal_db.search()
#cal_db.save()
#