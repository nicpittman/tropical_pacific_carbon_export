#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:56:43 2020

Cut out the NPP for each mooring

Has been updated to include the TPCA chlorophyll. Downloadable 
from http://dx.doi.org/10.25914/5dccbd3b64bdc (Pittman et al., 2019. JGR Oceans)

Chlorophyll will also need to be downloaded from https://oceandata.sci.gsfc.nasa.gov/ 
This script has been designed for chlor_a products for viirs, modis, meris and seawifs.

@author: N. Pittman 2020
"""

import xarray as xr
import os
import requests
import numpy as np



def Download_TPCA():
    '''
    Function to download TPCA MODIS and SeaWIFS from nci.org.
    
    Possible to also use DAPPS through xarray and save the files
    rather than using requests.
    '''
    path_sw='datasets/chl/tpca/seawifs/'
    if not os.path.isdir(path_sw):
        print('Creating directory: ',path_sw)
        os.makedirs(path_sw)
       
    path_mod='datasets/chl/tpca/modis/'
    if not os.path.isdir(path_mod):
        print('Creating directory: ',path_mod)
        os.makedirs(path_mod)
    
    
    tpca_link=['http://dapds00.nci.org.au/thredds/fileServer/ks32/CLEX_Data/TPCA_reprocessing/v2019_01/']
    sensors=['SeaWiFS/tpca_seawifs_','MODIS-Aqua/tpca_modis_aqua_'] #and then year
    sens=['sw','mod']           
    #Download SeaWiFS files from the above array, spaced by each year.
    
    for i in range(0,2): #To do SeaWiFS and then MODIS
        for yr in np.arange(1997,2020):
            if i==0: #SW
                sensor=tpca_link[0]+sensors[0]+str(yr)+'.nc'
                path=path_sw
            elif i==1: #MODIS
                sensor=tpca_link[0]+sensors[1]+str(yr)+'.nc'
                path=path_mod
        
            #Start the download
            try:
                r = requests.get(sensor)#,timeout=s20)
                fileloc=path+sensors[0].split('/')[1]+str(yr)+'.nc'
                if r.status_code!=404:
                    with open(fileloc, 'wb') as f:
                        f.write(r.content)
                    print('Downloaded: ' + sens[i] + str(yr))
                else:
                    print(i,str(r.status_code))
            except KeyboardInterrupt:
                import sys
                sys.exit()
            except:
                print(str(yr)+ sens[i]+'  Unavailable')
            pass

def process_npp_moorings():
    buff=0.5 #in degrees.
    l=0
    lats=[l+buff,l-buff]

    lns=[165,190,205,220,235,250]
    lons=[[lns[0]-buff,lns[0]+buff],
          [lns[1]-buff,lns[1]+buff],
          [lns[2]-buff,lns[2]+buff],
          [lns[3]-buff,lns[3]+buff],
          [lns[4]-buff,lns[4]+buff],
          [lns[5]-buff,lns[5]+buff]]
    #lons=[[164.5,165.5],[189.5,190.5],[204.5,205.5],[219.5,220.5],[234.5,235.5],[249.5,250.5]] #One degree over the moorings.
    #Might attempt some different sizes in the future.

    mooring_sites=['165E','170W','155W','140W','125W','110W']

    npp_models=['datasets/npp_satellite/vgpm_mod_nc/*',
                'datasets/npp_satellite/cbpm_mod_nc/*',
                'datasets/npp_satellite/eppley_mod_nc/*',
                
                'datasets/npp_satellite/cafe_mod_nc/*',     
                'datasets/npp_satellite/cafe_sw_nc/*',   
                
                'datasets/npp_satellite/vgbm_sw_nc/*',
                'datasets/npp_satellite/cbpm_sw_nc/*',
                'datasets/npp_satellite/eppley_sw_nc/*',  
                
                'datasets/npp_satellite/vgpm_viirs_nc/*', 
                'datasets/npp_satellite/cbpm_viirs_nc/*', 
                'datasets/npp_satellite/eppley_viirs_nc/*']
    out_path='processed/npp_mooring_timeseries/'

    for i,x in enumerate(npp_models):
        try:
            dat=xr.open_mfdataset(x,concat_dim='time',combine='nested')
            for ii,ll in enumerate(lons):
                datslice=dat.sel(lat=slice(lats[0],lats[1]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).npp.to_series()

                dataset=npp_models[i]
                mooring_site=mooring_sites[ii]
                series_name= dataset.split('/')[2]+'_'+mooring_site+'.csv'
                
                print('saved to: '+out_path+series_name)        
                if not os.path.isdir(out_path):
                    print('Creating directory: ',out_path)
                    os.makedirs(out_path)   
                datslice.to_csv(out_path+series_name)
                
            
        except:
            print('Skipped: '+x)
            pass
        
process_npp_moorings()    
#Depending where data stored probably want to comment the breakpoint out.
print('CHECK HERE - There is a breakpoint installed.')
print('If you have TPCA or NASA chlor_a downloaded, please ignore this, and then update the file paths accordingly')
print('If not, TPCA will be downloaded automatically now into datasets/chl/tpca/seawifs/')
print('NASA Chlor_a is a little harder, you will need an account with NASA')
print('https://oceancolor.gsfc.nasa.gov/data/download_methods/')
print('A shell script in datasets/chl/download_NASA_chlora.sh is provided, you will need to follow the instructions at:')
print('https://oceancolor.gsfc.nasa.gov/data/download_methods/ to get the auth cookie working')

Download_TPCA() #Downloads Pittman 2019 data

print('TPCA Downloaded, make sure chlor_a is also available')

#sw_tpca='/g/data/ua8/ocean_color/TPCA_reprocessing/SeaWiFS/*nc'
#mod_tpca='/g/data/ua8/ocean_color/TPCA_reprocessing/MODIS-Aqua/*nc'
sw_tpca='datasets/chl/tpca/seawifs/*nc'
mod_tpca='datasets/chl/tpca/modis/*nc'

buff=0.5 #in degrees.
l=0
lats=[l+buff,l-buff]

lns=[165,190,205,220,235,250]
lons=[[lns[0]-buff,lns[0]+buff],
      [lns[1]-buff,lns[1]+buff],
      [lns[2]-buff,lns[2]+buff],
      [lns[3]-buff,lns[3]+buff],
      [lns[4]-buff,lns[4]+buff],
      [lns[5]-buff,lns[5]+buff]]
 
mooring_sites=['165E','170W','155W','140W','125W','110W']
out_path='processed/npp_mooring_timeseries/'
   

#paths=['/g/data/ua8/ocean_color/tropics/','/day/9km/chlor_a/*nc']
paths='datasets/chl/chlor_a/'

chl_models=['seawifs','modis','meris','viirs']

#This uses the downloaded versions in datasets/chl/download shell script. If the data is stored locally, you can indicate it like the commented paths. You might also need to change the path definitions below as well (commented line).

for i,x in enumerate(chl_models):
        path = paths+x+'/*.nc'
        #path=paths[0]+x+paths[1]
        dat=xr.open_mfdataset(path,concat_dim='time',combine='nested',engine='h5netcdf')
        for ii,ll in enumerate(lons):
            datslice=dat.sel(lat=slice(lats[0],lats[1]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).chlor_a.to_series()

            dataset=chl_models[i]
            mooring_site=mooring_sites[ii]
            series_name= dataset+'_chlor_a_'+mooring_site+'.csv'
            
            print('saved to: '+out_path+series_name)        
            datslice.to_csv(out_path+series_name)
                       
            if i==0:
                chlor_dat=xr.open_mfdataset(sw_tpca,concat_dim='time',combine='nested')
                chlslice=chlor_dat.sel(lat=slice(lats[0],lats[1]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).chl_tpca.to_series()
                chlslice.to_csv(out_path+'seawifs_chl_tpca_'+mooring_site+'.csv') 
                
                chlor_dat=xr.open_mfdataset(mod_tpca,concat_dim='time',combine='nested')
                chlslice=chlor_dat.sel(lat=slice(lats[0],lats[1]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).chl_tpca.to_series()
                chlslice.to_csv(out_path+'modis_chl_tpca_'+mooring_site+'.csv') 
                
            #Will save this 5 times but whatever... Sorry for the inefficiency.
        
