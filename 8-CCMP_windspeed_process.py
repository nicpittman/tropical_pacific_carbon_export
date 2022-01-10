#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:11:42 2021

@author: npittman
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import requests
import os
import sys
import xesmf as xe


glob=True
download=False

#http://data.remss.com/ccmp/v02.0/Y2000/M02/CCMP_Wind_Analysis_200002_V02.0_L3.5_RSS.nc

years=np.arange(1987,2020)
months=np.arange(1,13)
raw_link=['http://data.remss.com/ccmp/v02.0/Y','/M','/CCMP_Wind_Analysis_','_V02.0_L3.5_RSS.nc']
links=[]
for y in years:
    y=str(y)
    for m in months:
        if m<10:
            m='0'+str(m)
        else:
            m=str(m)
        links.append(raw_link[0]+y+raw_link[1]+m+raw_link[2]+y+m+raw_link[3])

# %%

if download==True:
    print('Downloading global CCMP windspeed data from www.data.remss.com/ccmp/v02.0 ')
    #print(links)
    if not os.path.isdir('datasets/ws_ccmp'):
         print('Creating directory: ','datasets/ws_ccmp')
         os.makedirs('datasets')  
    
    def check_exists(fp):
        try:
            d=xr.open_dataset(fp)
            return True
        except:
            return False
        
    for l in links:
        fileloc='datasets/ws_ccmp/'+l.split('/')[-1]
        if check_exists(fileloc)==False:
          
            try: 
                print('Downloading: '+l)
                r = requests.get(l,timeout=20)
                if r.status_code!=404:
                    with open(fileloc, 'wb') as f:
                        f.write(r.content)
                        print('Saved to: '+fileloc)
                else: 
                    print(str(r.status_code))
            except KeyboardInterrupt:
                sys.exit()
            except:
                print(fileloc+' is unavailable')

    
# %%
#Need to shuffle the landschutzer data a little bit.
landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
landschutzer=xr.open_dataset(landsch_fp)
landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')
land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
land_pac.to_netcdf('processed/flux/landshutzer.nc',engine='h5netcdf',mode='w')

#Process and clean windspeed (CCMP)


# # THIS NEEDS TO BE RUN ONCE BUT CAN be memory intensive

w_ccmp=xr.open_mfdataset('datasets/ws_ccmp/*.nc') #Downloaded manually
w_ccmp['time']=w_ccmp.time.astype('datetime64[M]')

suffix=''
if glob==True:
    suffix='_global'
else:
    w_ccmp=w_ccmp.sel(longitude=slice(120,290),latitude=slice(-20,20))
    
w_ccmp=w_ccmp.rename({'latitude':'lat','longitude':'lon'})

#w_ccmp['windspeed']=np.sqrt((w_ccmp.uwnd**2)+(w_ccmp.vwnd**2))
try:
    print('saving..')
    w_ccmp.to_netcdf(f'datasets/CCMP_windspeed{suffix}.nc')
    
    print('saved')
except:
    pass

# %% Regrid the windspeed

ws=xr.open_dataset(f'datasets/CCMP_windspeed{suffix}.nc')
#ws_ccmp=np.sqrt((w_ccmp.uwnd**2)+(w_ccmp.vwnd**2))   
landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
landschutzer=xr.open_dataset(landsch_fp)
landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')
land_pac=landschutzer.fgco2_smoothed.to_dataset(name='co2')
#land_pac=land_pac.sel(lat=slice(-20,20))

lat_rad=land_pac.lat.diff(dim='lat').mean().values/2
lon_rad=land_pac.lon.diff(dim='lon').mean().values/2
landlats=np.linspace(land_pac.lat.min().values-lat_rad,land_pac.lat.max().values+lat_rad,len(land_pac.lat)+1)
landlons=np.linspace(land_pac.lon.min().values-lon_rad,land_pac.lon.max().values+lon_rad,len(land_pac.lon)+1)
    
lat_rad=ws.lat.diff(dim='lat').mean().values/2
lon_rad=ws.lon.diff(dim='lon').mean().values/2
lats=np.linspace(ws.lat.min().values-lat_rad,ws.lat.max().values+lat_rad,len(ws.lat)+1)
lons=np.linspace(ws.lon.min().values-lon_rad,ws.lon.max().values+lon_rad,len(ws.lon)+1)
    
ws.coords['lon_b']=lons
ws.coords['lat_b']=lats
ws['lat_b']=lats#(['lonb','latb'],x[1])
ws['lon_b']=lons#(['lonb','latb'],x[0])

land_pac.coords['lon_b']=landlons
land_pac.coords['lat_b']=landlats
land_pac['lat_b']=landlats#(['lonb','latb'],x[1])
land_pac['lon_b']=landlons##(['lonb','latb'],x[0])
print('Regridding WS')
regridder = xe.Regridder(ws, land_pac, 'conservative')
wsreg=regridder(ws)
wsreg.to_netcdf(f'processed/CCMP_ws_1deg{suffix}.nc',mode='w') #engine='h5netcdf'
    
    
