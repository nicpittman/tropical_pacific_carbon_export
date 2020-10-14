#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:31:01 2019

@author: npittman
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
import lxml

def downloader(urls,mooring=None, meta=False):
    """" Function designed to download TaoTriton data from: """
    
    path='datasets/tao/tao_co2/pieces/downloaded/'
    meta_path=path+'meta/'
    if meta==True:
        path=meta_path
        
    def exists(fileloc):
        if str(fileloc[-4:])=='.pdf':
            try:
                open(fileloc)
                return True
            except:
                return False
        #print(fileloc)
        try:
            pd.read_csv(fileloc,header=3)
            return True
        except:
            return False

    if not os.path.isdir(path):
        print('Creating directory: ',path)
        os.makedirs(path)    
    if not os.path.isdir(meta_path):
        print('Creating directory: ',meta_path)
        os.makedirs(meta_path)    
        
    if mooring != None:
        for moor in mooring:
            if not os.path.isdir(path+moor):
                print('Creating directory: ',path+moor)
                os.makedirs(path+moor)    
        
    file_locations=[]
    
    for url in urls:
        while True:
        #Download the files to their file name in the directory we just created.
        #ie: seawifs_data/S2000001.L3m_DAY_RRS_Rrs_443_9km.nc
            try:
                if mooring ==None:
                    fileloc=path+url.split('/')[-1]
                else:
                    if str(url.split('/')[-1][7])!='_':
                        if str(url.split('/')[-1][3])!='_':
                            moorpath=url.split('/')[-1][0:7]+'_'+url.split('/')[-1][7:9]
                            fileloc=path+moorpath+'/'+url.split('/')[-1]
                        else:
                            moorpath=url.split('/')[-1][0:3]+url.split('/')[-1][4:11]
                            fileloc=path+moorpath+'/'+url.split('/')[-1]
                    else:
                        moorpath=url.split('/')[-1][0:10]
                        fileloc=path+moorpath+'/'+url.split('/')[-1]
            except:
               # print('something broke at:',url)
                continue
            if fileloc[-3:]=='txt':
                fileloc=fileloc[0:-3]+'csv'
            print(url)
            if exists(fileloc):
                print('Exists: ',fileloc)
                file_locations.append(fileloc)
                break
            r = requests.get(url)#,timeout=s20)
            with open(fileloc, 'wb') as f:
                f.write(r.content)
    
            #Ensure that the file actually downloaded, this can fail sometimes for some reason, maybe too soon.
            #time.sleep(1) #Can fail sometimes so maybe this is a fix
            if (r.status_code==200) & (exists(fileloc)==True):
                print('Downloaded: ',fileloc)
                file_locations.append(fileloc)
                break
            else:
                print('Download failed:', fileloc,'status:',r.status_code)
                        
    return file_locations


def gather_links(url_list):
    urls=[]
    metadata_urls=[]
    
    for url in url_list: #Loop through the different wavelengths
        print("Gathering URLS for TaoTRITON")
        print(url)
        r = requests.get(url)
        soup = BeautifulSoup(r.content,features="lxml")
        for tag in soup.find_all('td'): #Loop through years
            tags=str(tag).split('"')
            try:
                this_tag=tags[1] #Get a year's URL
                if (str(this_tag[-4:])=='.csv') or (str(this_tag[-4:])=='.txt'):
                    print('DATA: '+this_tag)
                    urls.append(url+this_tag)
                elif str(this_tag[-8:])=='Meta.pdf':
                    print('META: '+this_tag) 
                    metadata_urls.append(url+this_tag)
            except:
                pass
    return urls,metadata_urls



url_list=['https://www.nodc.noaa.gov/archive/arc0061/0113238/5.5/data/0-data/',	#TAO165E_0N
          'https://www.nodc.noaa.gov/archive/arc0063/0117073/4.4/data/0-data/', #TAO165E_8S
          'https://www.nodc.noaa.gov/archive/arc0051/0100078/8.8/data/0-data/', #TAO170W_0N
          'https://www.nodc.noaa.gov/archive/arc0051/0100084/6.6/data/0-data/', #TAO155W_0N
          'https://www.nodc.noaa.gov/archive/arc0051/0100077/7.7/data/0-data/', #TAO140W_0N
          'https://www.nodc.noaa.gov/archive/arc0051/0100076/8.8/data/0-data/', #TAO125W_0N
          'https://www.nodc.noaa.gov/archive/arc0061/0112885/6.6/data/0-data/', #TAO110W_0N
          'https://www.nodc.noaa.gov/archive/arc0051/0100079/1.1/data/0-data/TAO170W_2S_Aug07_Aug08/', #TAO170W_2S_Aug07_Aug08 1
          'https://www.nodc.noaa.gov/archive/arc0051/0100079/1.1/data/0-data/TAO170W_2S_Jun98_Nov04/'] #TAO170W_2S_Jun98_Nov04 2
mooring=['TAO165E_0N','TAO165E_8S','TAO170W_0N','TAO155W_0N','TAO140W_0N','TAO125W_0N','TAO110W_0N','TAO170W_2N','TAO170W_2S']

urls,metadata=gather_links(url_list) #Can comment out onece complete.
downloader(urls,mooring) #Can comment out onece complete.
downloader(metadata,meta=True) #Can comment out onece complete.

data=[]
#Clean up and save into the pieces folder.
for x, folder in enumerate(glob.glob('datasets/tao/tao_co2/pieces/downloaded/T*')):
    print('\n'+folder)
    files=glob.glob(folder+'/*')
    indexes=[]
    data.append([])
    #Remove files including LOG or QC from our list of files.
    for file in files:
        if 'Log' in file:
            index=files.index(file)
            indexes.append(index)    
    for i in reversed(indexes):
        files.pop(i)

    
    for floc in files:
        #print(floc)
        check_head=0
        while True:
            #print(check_head)
            try:
                dat=pd.read_csv(floc,header=check_head)
                if ('mooring' in str.lower(dat.columns[0])) and ('lat' in str.lower(dat.columns[1])) :
                    print('SUCCESS\t',floc)
                    #print(dat.head(3))
                    data[x].append(dat)
                    break
                else:
                    check_head+=1
                    if check_head>10:
                        print('FAIL\t',floc)
                        break
            except:
                if check_head>10:
                    print('FAIL\t',floc)
                    break
                else: 
                    check_head+=1
                    #print('continue?')
                    continue
for dat in data:
        if dat==[]:
            print('Empty Array')
            continue
        dat1=pd.concat(dat,sort=True)
        try:
            dat1.Date=pd.to_datetime(dat1.Date,errors='coerce')
            dat1=dat1.sort_values('Date')
            #dat1 = dat1.dropna(subset=['Date'])
        except:
            print('Aborted')
            pass
        
        try:
            fname=dat1.iloc[5].Mooring
        except:
            fname=dat1.iloc[5]['Mooring Name']

        print(fname)
        dat1.to_csv('datasets/tao/tao_co2/pieces/'+str(fname)+'.csv')
        print('Saved:',fname)
        print(dat1.Date.head(5))

        #Possible Plotting        
        try:
            try:
                plt.plot(dat[1:].Date.astype(np.datetime64),dat[1:].SST.astype(float))
                plt.show()
            except:
                plt.plot(dat[1:].Date.astype(np.datetime64),dat[1:]['SST (C)'].astype(float))
                plt.show()
        except:
            pass
    
