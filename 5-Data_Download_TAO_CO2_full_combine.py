#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 13:51:04 2020

@author: npittman

Kind of duplicated the last script, but these are the suggested links. Nice to try this method
And combine data where possible. (155W)
"""
import requests                     
import xarray as xr                
import numpy as np                   
import matplotlib.pyplot as plt    
import os
from bs4 import BeautifulSoup
import multiprocessing
import sys
import pandas as pd
from tarfile import TarFile
import gzip
import subprocess
import calendar
import shutil
import glob

def downloader(url):
    path='datasets/tao/tao_co2/full_timeseries/downloaded/'
        
    def exists(fileloc):
        if str(fileloc[-4:])=='.pdf':
            try:
                open(fileloc)
                return True
            except:
                return False
        #print(fileloc)
        try:
            pd.read_csv(fileloc,header=50)
            return True
        except:
            return False

    if not os.path.isdir(path):
        print('Creating directory: ',path)
        os.makedirs(path)    
    
    fileloc=path+url.split('/')[-1][:-3]+'csv' #Rename to csv
    r = requests.get(url)#,timeout=s20)
    with open(fileloc, 'wb') as f:
        f.write(r.content)

    #Ensure that the file actually downloaded, this can fail sometimes for some reason, maybe too soon.
    #time.sleep(1) #Can fail sometimes so maybe this is a fix
    if (r.status_code==200):
        print('Downloaded: ',fileloc)

    else:
        print('Download failed:', fileloc,'status:',r.status_code)
                        
    #dat.to_csv(fl,index=False)
    return fileloc


def add_19972005_to_155w():
    """
    Bit of a hack script to add the 1997-2010 data to the 155W mooring. 
    Obtained from other sources. See the other TAO downloadscript.
    """
    t155=pd.read_csv('datasets/tao/tao_co2/full_timeseries/downloaded/TAO155W.csv',comment='#',delimiter='\t')

    fs=['datasets/tao/tao_co2/pieces/downloaded/TAO155W_0N/TAO155W_0N_Jun2005_Aug2008.csv',
        'datasets/tao/tao_co2/pieces/downloaded/TAO155W_0N/TAO155W_0N_Nov1997_Jun2005.csv']
    t155['delta_pCO2']=t155.pCO2_sw-t155.pCO2_air
    all_data=t155
    for i in [0,1]:
        dat=pd.read_csv(fs[i],header=0)[1:]
        try:
            datet=pd.to_datetime(dat.Date+' '+dat.Time,format='%m/%d/%Y %H:%M')
        except:
            dd=dat['Date'].astype(str).str.zfill(8)
            datet=pd.to_datetime(dd +' '+dat.Time,format='%d%m%Y %H:%M:%S')
        try:
            sst=dat.SST
            sss=dat.SSS
            pco2_air=dat.fCO2_Air_sat.astype(float)
            pco2_sw=dat.fCO2_SW_sat.astype(float)
            xco2=dat.xCO2_Air_dry.astype(float)
            data=pd.DataFrame({'datetime_utc':datet,
                      'SST':sst, 
                      'SSS':sss,
                      'pCO2_sw':pco2_sw, 
                      'pCO2_air':pco2_air,
                      'xCO2_air':xco2,
                      'pH_sw':np.nan,
                      'delta_pCO2':pco2_sw-pco2_air})
        except:
            print("Vars don't exist in 1997 data")
            data=pd.DataFrame({'datetime_utc':datet,
              'SST':np.nan, 
              'SSS':np.nan,
              'pCO2_sw':np.nan, 
              'pCO2_air':np.nan,
              'xCO2_air':np.nan,
              'pH_sw':np.nan,
              'delta_pCO2':dat.Delta_pCO2})
    
        all_data=all_data.append(data)
        
    all_data['datetime_utc']=pd.to_datetime(all_data.datetime_utc)
    all_data=all_data.sort_values('datetime_utc')
    all_data.to_csv('datasets/tao/tao_co2/full_timeseries/downloaded/TAO155W_full.csv',index=False)
    print('Combined TAO155W.csv with extra data')

def downloadrun(flocs=[]):
    for url in url_list:
        flocs.append(downloader(url))
    return flocs


url_list=['https://www.pmel.noaa.gov/co2/timeseries/TAO110W.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO125W.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO140W.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO155W.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO165E.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO8S165E.txt',
          'https://www.pmel.noaa.gov/co2/timeseries/TAO170W.txt']

#Can turn download on / off here
print('Check here because some downloading and processings can be commented on and off')
flocs=downloadrun() #Can be commented out after run
add_19972005_to_155w() #Can be commented out after run
flocs=glob.glob('datasets/tao/tao_co2/full_timeseries/downloaded/*')

for fl in flocs:
    dat=pd.read_csv(fl,delimiter='\t',comment='#')
    if len(dat.columns)==1:
        dat=pd.read_csv(fl,delimiter=',',comment='#')
    dat=dat.rename(columns={'datetime_utc':'Date'})
    try:
        check=dat.delta_pCO2
    except:
        dat['delta_pCO2']=dat.pCO2_sw-dat.pCO2_air
        
    if dat.shape[1]==1:
        dat=pd.read_csv(fl,delimiter=',',comment='#')
        dat=dat.rename(columns={'datetime_utc':'Date'})
   
    dat['Date']=pd.to_datetime(dat['Date'])
    
    #plt.plot_date(dat.Date,dat.xCO2_air)
    #plt.title( fl.split('/')[1][:-4]+' xCO2')
    #plt.show()
    #plt.plot_date(dat.Date,dat.pCO2_sw-dat.pCO2_air,c='r')
    #plt.title( fl.split('/')[1][:-4]+' pCO2sw-pCO2atm')
    #plt.axhline(0,c='k')
    #plt.show()

    #plt.plot_date(dat.Date,dat.delta_pCO2)
    #plt.plot_date(dat.Date,dat.pCO2_sw-dat.pCO2_air,c='r',ms=2)
    #plt.title( fl.split('/')[1][:-4]+' pCO2sw-pCO2atm')
    #plt.axhline(0,c='k')
    #plt.show()
    
    path='datasets/tao/tao_co2/'
    if not os.path.isdir(path):
        print('Creating directory: ',path)
        os.makedirs(path)    
    dat.to_csv(path+'/'+fl.split('/')[-1],index=False)
    print('saved to: '+path+fl.split('/')[-1])
