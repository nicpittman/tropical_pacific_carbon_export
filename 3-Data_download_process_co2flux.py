#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:05:38 2020

Similar to 2-Data_Cut_mooring_npp_timeseries
But is focused for the JMA CO2 flux data and 

@author: npittman
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import requests
import os
from zipfile import ZipFile


'''
Iida et al. (2015; http://www.data.jma.go.jp/gmd/kaiyou/english/co2_flux/co2_flux_data_en.html; JMA_MLR). ETH_SOM-
FFN was estimated by a combined two-step neural network (a
self-organizing map and a feed-forward network) approach with SST,
SSS, mixed layer depth, and chlorophyll-a. JMA_MLR was estimated
using multiple linear regression of SST, SSS, and chlorophyll-a.
'''

def process_jams():
        
    #Download Co2 flux data
    #Want to download Japan
    
    #But also: https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc
    #Download this one manually.
    #(https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html)
    
    path='datasets/co2/JMA_co2/'
    if not os.path.isdir(path):
        print('Creating directory: ',path)
        os.makedirs(path)    
            
    jamslinks=['http://www.data.jma.go.jp/gmd/kaiyou/data/english/co2_flux/grid/JMA_co2map_','.ZIP']
    #Download files from the above array, spaced by each year.
    for i in np.arange(1990,2019):
        link=jamslinks[0]+str(i)+jamslinks[1]
        r = requests.get(link)#,timeout=s20)
        fileloc=path+'JMA_co2'+str(i)+'.zip'
        with open(fileloc, 'wb') as f:
            f.write(r.content)
        
        with ZipFile(fileloc, 'r') as zipObj:
        # Extract all the contents of zip file in current directory
            zipObj.extractall(path)
        print(i)
        os.remove(fileloc)
        
    
    #Process it into a single file
    yearlist=np.arange(1990,2019)
    dat=xr.open_mfdataset('datasets/co2/JMA_co2/*nc',decode_times=False,concat_dim='Year')
    dat['Year']=yearlist
    dat=dat.rename({'time':'Month'})
    datasets=[]
    
    for y in yearlist:
        yearly=dat.sel(Year=y)
        for m in np.arange(0,12):
            month_box=yearly.sel(Month=m)
            month=month_box.Month.values+1
            if len(str(month))==1:
                month='0'+str(month)
                
            month_box=month_box.assign_coords(time=(np.datetime64(str(month_box.Year.values)+'-'+str(month))))
            month_box=month_box.drop(['Month','Year'])
            datasets.append(month_box)
            
    new_flux=xr.concat(datasets,dim='time')
    
    
    
    out_path='processed/flux/'      
    if not os.path.isdir(out_path):
        print('Creating directory: ',out_path)
        os.makedirs(out_path)   
    new_flux.to_netcdf(out_path+'jma_flux.nc')
    return True

def process_yasanaka():
	"""
	Save as .nc
	Open Yasunaka (2019) data from binary to xarray and save to labelled netcdf. 
	Data obtained from personal communication. I am unable to share this. 

	3 variables: co2flux, error, pco2
	35 years 1981 - 2015
	12 months 0 - 11
	40 lats -19.5 - 19.5
	200 lons 100 - 300E

	@author: Nic.Pittman@utas.edu.au

	Yasunaka, S., Kouketsu, S., Strutton, P. G., Sutton, A. J., Murata, A., Nakaoka, S., & Nojiri, Y. (2019). Spatio-temporal variability of surface water pCO2 and nutrients in the tropical Pacific from 1981 to 2015. Deep Sea Research Part II: Topical Studies in Oceanography, 169, 104680.
	"""
	fps=['datasets/co2/yasanaka_co2/co2flux_1807.bin',                #File paths - Array for the 3 variables.
	     'datasets/co2/yasanaka_co2/error_1807.bin',
	     'datasets/co2/yasanaka_co2/pco2_1807.bin']
	     
	arrs=[]
	year=np.arange(1981,2016)            #Setting years
	month=np.arange(0,12)                #Setting months
	lat=np.arange(-19.5,20.5)            #Setting lats
	lon=np.arange(100.5,300.5)  
	print('Processing')
	for fp in fps:
	    data = np.array('B')                        #Make a numpy bytes array
	    f = open(fp, "r")                           #Open the file as bytes
	    a = np.fromfile(f, dtype=np.float32)        #Convert into an array of 3360000
	    count=0                                     #Set count to run through bytes
	    data=np.ndarray([35,12,40,200])             #Set the data array size 
	    for iy in range(0,35):                      #Loop through them all to count the bytes in order
	       for im in range(0,12):
	           for la in range(0,40):
	               data[iy,im,la,:]=a[count:count+200]  #Read and append the bytes to the data array
	               count=count+200                      #Cycle our counts up
	    piece=xr.DataArray(data,dims=['year','month','lat','lon'])
	    piece.name=fp.split('/')[1].split('_')[0]
	    arrs.append(piece) #Make an xarray so we can save as labeled netcdf.

		       
	  
	out=xr.merge(arrs)
	out['year']=year         #Setting years
	out['month']=month       #Setting months
	out['lat']=lat        #Setting lats
	out['lon']=lon       #Setting lons
	out['flux_masked']=out.co2flux.where(out.error<0.9)
	#out['var']=['co2flux','error','pco2s']
	#print(out)                                 #Check it works
	#out.flux_masked.sel(year=2000).mean(dim='month').plot() #Plot it to check ok
	
	#Merge and convert dates
	print('Merging Dates')
	datasets=[]
	month_box=out
	for y in month_box.year.values:
	    for m in month_box.month.values:
	        box=month_box.sel(year=y,month=m)
	        m+=1
	        if len(str(m))==1:
	            m=str(0)+str(m)
	        box=box.assign_coords(time=(np.datetime64(str(y)+'-'+str(m))))
	        box=box.drop(['month','year'])
	        datasets.append(box)
	out1=xr.concat(datasets,dim='time')
	#Convert the month and years into a numpy datetime. Best way besides just looping through surely...

	out1.to_netcdf('datasets/co2/yasanaka_co2/Yasunaka_pCO2_flux.nc')           #Save to labelled netcdf. 
	print('Finished')
	return True
print('')
print('')
print('')
print('------------------- NOTICE --------------------')
print('You will need to download the Landschutzer product manually. You can use curl for example.')
print('source: https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html and can be downloaded like:')
print('curl https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc --output spco2_MPI_SOM-FFN_v2018.nc')
print('You will then need to use NCO to convert this variable name so they can be opened. Ie:')
print('ncrename -v date,t MPI-SOM_FFN_SOCCOMv2018.nc')
print('')
print('If this is your first time running, #process_jams() and process_yasanaka() should be working but these can be commented out if rerunning this analysis')
print('Yasanaka Data was obtained from private communication. Please email her for this dataset- Paper DOI: https://doi.org/10.1016/j.dsr2.2019.104680')
print('---------------------------------------------')

#Uncomment this to preprocess the JAMS and Yasanaka data.
process_jams()
process_yasanaka()
out_path='processed/flux/'    
new_flux=xr.open_mfdataset(out_path+'jma_flux.nc')

#Include the Landschutzer products:
#Need to use NCO to convert if getting the MissingDimensionsError
chunkz={'time':10}
landschutzer=xr.open_mfdataset('datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc',chunks=chunkz)

landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False)		#EPIC 1 line fix for the dateline problem.

yasanaka=xr.open_mfdataset('datasets/co2/yasanaka_co2/Yasunaka_pCO2_flux.nc')

#Cutout 

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

dat=new_flux

moorings=[]
dat=xr.open_mfdataset(out_path+'jma_flux.nc')

land_moorings=[]
yasanaka_moorings=[]

for ii,ll in enumerate(lons):
    datslice=dat.sel(lat=slice(lats[1],lats[0]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon'])
    datslice=datslice.assign_coords(Mooring=mooring_sites[ii])
    moorings.append(datslice)
    
    yasanaka_datslice=yasanaka.sel(lat=slice(lats[1],lats[0]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).assign_coords(Mooring=mooring_sites[ii])
    yasanaka_moorings.append(yasanaka_datslice)
    land_datslice=landschutzer.sel(lat=slice(lats[1],lats[0]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).assign_coords(Mooring=mooring_sites[ii])
    
    land_moorings.append(xr.merge([land_datslice.spco2_raw,
                        land_datslice.aco2,
                        land_datslice.fgco2_raw,
                        land_datslice.fgco2_smoothed,
                        land_datslice.sol,
                        land_datslice.kw]))
    
yasanaka_mooring_flux=xr.concat(yasanaka_moorings,dim='Mooring')    
yasanaka_mooring_flux.to_netcdf(out_path+'yasanaka_mooring_co2_flux.nc')

mooring_fluxes=xr.concat(moorings,dim='Mooring')
mooring_fluxes.to_netcdf(out_path+'JMA_mooring_co2_flux.nc') #Renamed this filepath to JMA_*


land_mooring_flux=xr.concat(land_moorings,dim='Mooring')
land_mooring_flux.to_netcdf(out_path+'landsch_mooring_co2_flux.nc')
