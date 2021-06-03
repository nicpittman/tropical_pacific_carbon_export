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

w_ccmp_a=xr.open_mfdataset('datasets/ws_ccmp/*.nc') #Downloaded manually
w_ccmp_a['time']=w_ccmp_a.time.astype('datetime64[M]')
w_ccmp=w_ccmp_a.sel(latitude=slice(-20,20))

# w_ccmp_b=xr.open_mfdataset('datasets/CCMP_winds.nc') #Bulk ErDap download
# dt=w_ccmp_b.indexes['time'].to_datetimeindex()
# w_ccmp_b['time']=dt

# w_ccmp=xr.merge([w_ccmp_b,w_ccmp_a])

w_ccmp=w_ccmp.sel(longitude=slice(120,290),latitude=slice(-20,20))
ws_ccmp=np.sqrt((w_ccmp.uwnd**2)+(w_ccmp.vwnd**2))
ws_ccmp=ws_ccmp.rename({'latitude':'lat','longitude':'lon'})
try:
    ws_ccmp.to_netcdf('datasets/CCMP_windspeed.nc')
    print('saved')
except:
    pass

# %% Regrid the windspeed

ws=xr.open_dataarray('datasets/CCMP_windspeed.nc')

ws=ws.to_dataset(name='ws')
   
landschutzer=xr.open_dataset('processed/flux/landshutzer.nc')
land_pac=landschutzer.fgco2_smoothed.to_dataset(name='co2')
land_pac=land_pac.sel(lat=slice(-20,20))

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
wsreg.to_netcdf('processed/CCMP_ws_1deg.nc',mode='w') #engine='h5netcdf'
    
    
