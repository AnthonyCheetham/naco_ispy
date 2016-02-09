# -*- coding: utf-8 -*-
"""

This module contains all of the database related code for the NACO-GTO pipeline
automator

TODO:
    - Encrypt passwords
    - Update entry method
    - Define database fields properly somewhere easy to access

"""

import numpy as np
import os,glob,time,psycopg2
import astropy.io.fits as pyfits
import astropy
from astropy.table import Table
import datetime

class obs_database(object):
    """ The observations database object
    """
    
    def __init__(self):
        ''' 
        Connects to the database and adds a conn object
        '''
        self.conn=psycopg2.connect(user='naco_pipeline',host='localhost',dbname='naco_pipeline',
              password='naco')
        
    def create(self):
        '''
        Create the database, with all of the right columns.
        This is all just hard-coded for now, but maybe it could be passed
        as an option if required
        '''
        
        # Command to create table
        create_command='''CREATE TABLE "Observations"
        ("ObsID" bigserial NOT NULL, "TargetID" bigint, "Band" text, "AGPM" boolean, 
          "AnalysisStatus" integer, "FieldRotation" double precision, 
           "r0" double precision, "t0" double precision, "Location" text, "Date" date,
          "WindSpeed" double precision, "ExpTime" double precision, 
          "SaturationLevel" double precision, "PsfXWidth" double precision, 
          "PsfYWidth" double precision, "Vmag" double precision, "Kmag" double precision,
          "PsfReference?" boolean, 
          CONSTRAINT "ObsID" PRIMARY KEY ("ObsID")) WITH ( OIDS=FALSE );'''
          
        cur=self.conn.cursor()
        cur.execute(create_command)
        cur.close()
        self.conn.commit()
        
class obs_table(object):
    ''' An alternative to the observations database object. This is just a locally
    stored astropy table.
    '''
    def __init__(self,filename='./obs_table.dat',table_format='ascii.csv'):
        
        self.filename=filename
        self.table_format=table_format
        
        # Define the columns and column names
        self.columns=('ObsID','TargetID','Band','AGPM','AnalysisStatus','FieldRotation',
                 'r0','t0','Location','Date','WindSpeed','ExpTime','SaturationLevel',
                 'PsfXWidth','PsfYWidth','Vmag','Kmag','PsfReference?')
        f=np.float64
        self.dtypes=(np.int_,np.int_,'S5',np.bool_,np.int_,f,
                f,f,'S200',datetime.datetime,f,f,f,
                f,f,f,f,np.bool_)    
        
        try:
            data=Table.read(filename,format=table_format)
            self.data=data
        except:
            print 'Observations Table not found!',filename
    
    def create(self):
        ''' Create a new table with all of the right columns for the database.'''

        # create the table
        data=Table(names=self.columns,dtype=self.dtypes,masked=False)
        self.data=data

        # save it to disk        
        print 'Saving table as',self.filename
        self.data.write(self.filename,format=self.table_format)
        
    def add_entry(self,obs_sequence):
        ''' A wrapper that adds the observing sequence data to the table of observations.
        Assumes that all of the data is stored in a dictionary as obs_sequence.dict'''
        
        self.data.add_row(obs_sequence.dict)
        # save it to disk        
        print 'Saving table as',self.filename
        self.data.write(self.filename,format=self.table_format)
        
class calib_table(obs_table):
    ''' An alternative to the calibrations database object. This is just a locally
    stored astropy table. This is a sub-class of the obs_table class
    '''
    
    def __init__(self,filename='./calib_table.dat',table_format='ascii.csv'):
                
        self.filename=filename
        self.table_format=table_format
        
        # Define the columns and column names
        self.columns=('CalibID','Band','Type','AnalysisStatus','Location',
                 'Date')
        f=np.float64
        self.dtypes=(np.int_,np.int_,'S5',np.bool_,np.int_,f,
                f,f,'S200',datetime.datetime,f,f,f,
                f,f,f,f,np.bool_)    
        
        try:
            data=Table.read(filename,format=table_format)
            self.data=data
        except:
            print 'Calibrations Table not found!',filename
        
        

# Maybe astropy tables are the best way to go?
db=obs_table()
db.create()
#
## try adding a value:
#entry={'ObsID':0,'TargetID':1,'Band':'Lp','AGPM':True,'AnalysisStatus':0,'FieldRotation':50.,
#       'Date':datetime.datetime(2015,11,20)}
#db.data.add_row(entry)
#x=db.data
##print x['Date']

        
#drop_command='DROP TABLE "Observations"'

# Some code so I can test this
#db=obs_database()
#db.create()
        
#try:
#    if not conn.closed:
#        conn.close()
#except:
#    pass
        
# Create a connection to the default database        
#conn=psycopg2.connect(user='naco_pipeline',host='localhost',dbname='naco_pipeline',
#                      password='naco')
# Create a cursor
#cur = conn.cursor()

# First delete the table
#cur.execute(drop_command)

# Then create a new one
#cur.execute(command)
#db.create()

# Make a random table and add some data
#command="CREATE TABLE test (id serial PRIMARY KEY, num integer, data varchar);"
#cur.execute(command)
#cur.execute("INSERT INTO test (num, data) VALUES (%s, %s)",(100, "abc'def"))
#cur.fetchone()

#conn.commit()
#cur.close()
#
#conn.close()
