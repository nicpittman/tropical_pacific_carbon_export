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
        dat=xr.open_mfdataset(x,concat_dim='time')
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
    
    
#Depending where data stored probably want to comment the breakpoint out.
print('CHECK HERE - There is a breakpoint installed here but you may need to reproduce this if the tpca and chlor_a files are not in processed/npp_mooring_timeseries/')
import sys
sys.exit()

#This does the samea as above, but for chlorophyll datasets.
#A hack method but this should copy-paste to where the data lives. + top 35 lines.
        
#This runs on Gadi super computer. As the data is stored here. 
#Or you will need to change this and dow','/day/9km/chlor_a/*nc']

sw_tpca='/g/data/ua8/ocean_color/TPCA_reprocessing/SeaWiFS/*nc'
mod_tpca='/g/data/ua8/ocean_color/TPCA_reprocessing/MODIS-Aqua/*nc'

#sw_tpca='datasets/tpca/modis/*nc' #Actually the second one but I only have seawifs stored on gadi. #'datasets/tpca/seawifs/*nc'
#mod_tpca='datasets/tpca/modis/*nc'nload all of the chlorophyll data locally.
paths=['/g/data/ua8/ocean_color/tropics/']

chl_models=['seawifs','modis','meris','viirs']

for i,x in enumerate(chl_models):

        path=paths[0]+x+paths[1]
        dat=xr.open_mfdataset(path,concat_dim='time',combine='nested')
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
                
            #Will save this 5 times but whatever haha.
        
