#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:51:54 2020
A collection of cleanup functions that should just run. 
Just keep them in one place and clean up the file structure a bit.

Expects that the entire pipeline up until now has been completed. 
Hopefully this all works because will be hard to debug!!
Things might be out of order so need to check this.

@author: npittman
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import xesmf as xe
import cbsyst as cb
import os
from carbon_math import *


    
def find_enso_events(threshold=0.5):
    '''
    A function to pull ENSO data from our datasets/indexes/meiv2.csv
    save events (months) stronger than threshold (0.5 by default)
    'processed/indexes/el_nino_events.csv'
    'processed/indexes/la_nina_events.csv'
    Returns
    -------
    None.

    '''
    #enso=pd.read_csv('datasets/indexes/meiv2.csv',index_col='Year')
    enso=pd.read_csv('datasets/indexes/meiv2.csv',index_col=0,header=None)
    enso_flat=enso.stack()
    enso_dates=pd.date_range('1979','2020-07-01',freq='M')- pd.offsets.MonthBegin(1) #Probably want to check this is correct if updating.
    enso_timeseries=pd.DataFrame({'Date':enso_dates,'mei':enso_flat})
    
    #Check if we are in or out of an event so far
    el_event=False
    la_event=False
    el_startdate=''
    la_startdate=''
    elnino=pd.DataFrame()
    lanina=pd.DataFrame()
    
    for i,today in enumerate(enso_timeseries.Date):
        val=enso_timeseries.mei.iloc[i]
        if val>=threshold:
            if el_event==False:  #And we havent yet entered an event
                el_startdate=today
                el_event=True
            else:
                pass
                #Dont need to do anything because it will get caught later
        else:
            if el_event==True:
                elnino=elnino.append({'start':el_startdate.to_datetime64(),
                                      'end':enso_timeseries.Date.iloc[i-1],
                                      'mei':enso_timeseries.mei.iloc[i-1]},ignore_index=True)
                el_event=False
    
    for i,today in enumerate(enso_timeseries.Date):
        val=enso_timeseries.mei.iloc[i]
        if val<=-threshold:
            if la_event==False:  #And we havent yet entered an event
                la_startdate=today
                la_event=True
            else:
                pass
                #Dont need to do anything because it will get caught later
        else:
            if la_event==True:
                lanina=lanina.append({'start':la_startdate.to_datetime64(),
                                      'end':enso_timeseries.Date.iloc[i-1],
                                      'mei':enso_timeseries.mei.iloc[i-1]},ignore_index=True)
                la_event=False
    
    print(elnino)
    print(lanina)
    elnino.to_csv('processed/indexes/el_nino_events.csv')
    lanina.to_csv('processed/indexes/la_nina_events.csv')


def combine_csvs_to_nc():
    '''
    Combine all our data into daily, weekly, monthly or all data files. 
    This should be done at the end of 7ab already so may not be essential to run this one.
    Previously known as data day average.
    '''
    moorings=['110W','125W','140W','155W','170W','165E']
    mooring_int=[110,125,140,155,170,195]
    aavg_a=[]
    davg_a=[]
    wavg_a=[]
    mavg_a=[]
    
    for mooring in moorings:
        fp='processed/combined_dataset/'+mooring+'_combined.csv'
        dat=pd.read_csv(fp,index_col=False)
        #print(dat)
        
        dat['Date']=dat.Date.astype(np.datetime64)
        dat.set_index(pd.DatetimeIndex(dat.Date),inplace=True)
        
        alld = dat.to_xarray()#.drop('Unnamed: 0')
        davg = dat.resample('D').mean().to_xarray()#.drop('Unnamed: 0') #Day average
        wavg = dat.resample('W').mean().to_xarray()#.drop('Unnamed: 0') #Week average
        mavg = dat.resample('M').mean().to_xarray()#.drop('Unnamed: 0') #Month average
    
        aavg_a.append(alld)
        davg_a.append(davg)
        wavg_a.append(mavg)
        mavg_a.append(wavg)
        
        #plt.scatter(wavg.co2flux_gmyr,wavg.mod_vgbm)
        #plt.scatter(wavg.co2flux_gmyr,wavg.mod_cpbm)
        #plt.scatter(wavg.co2flux_gmyr,wavg.mod_cafe)
        #plt.scatter(wavg.co2flux_gmyr,wavg.mod_eppley)
        
        
    all_data=xr.concat(aavg_a,dim='Mooring')
    daily=xr.concat(davg_a,dim='Mooring')
    weekly=xr.concat(wavg_a,dim='Mooring')
    monthly=xr.concat(mavg_a,dim='Mooring')
    all_data.coords['Mooring']=moorings
    daily.coords['Mooring']=mooring_int
    weekly.coords['Mooring']=mooring_int
    monthly.coords['Mooring']=mooring_int
    
    fp='processed/combined_dataset/'
    all_data.to_netcdf(fp+'all_data.nc',engine='h5netcdf')
    daily.to_netcdf(fp+'day_data.nc')
    weekly.to_netcdf(fp+'week_data.nc')
    monthly.to_netcdf(fp+'month_data.nc')
    
def npp_csvs_to_nc():
    '''
    Cleanup function to combine the csvs for each mooring (A heap of files) into a single 4d xarray netcdf.    
    '''
    
    moorings=['110W','125W','140W','155W','170W','165E']
    #moorings=['155W','170W','165E']
    #moorings=['155W']
    
    flux_holder=[]
    chl_holder=[]
    for i, mooring_name in enumerate(moorings):
        
        #Primary Productivity Models and Time Series
        mod_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgpm_mod_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','mod_vgpm'],index_col=0).to_xarray()
        mod_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_mod_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','mod_cbpm'],index_col=0).to_xarray()
        mod_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_mod_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','mod_eppley'],index_col=0).to_xarray()
        mod_cafe=pd.read_csv('processed/npp_mooring_timeseries/cafe_mod_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','mod_cafe'],index_col=0).to_xarray()
        
        sw_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgbm_sw_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','sw_vgbm'],index_col=0).to_xarray()
        sw_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_sw_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','sw_cbpm'],index_col=0).to_xarray()
        sw_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_sw_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','sw_eppley'],index_col=0).to_xarray()
        sw_cafe=pd.read_csv('processed/npp_mooring_timeseries/cafe_sw_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','sw_cafe'],index_col=0).to_xarray()
        
        viirs_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgpm_viirs_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','viirs_vgpm'],index_col=0).to_xarray()
        viirs_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_viirs_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','viirs_cbpm'],index_col=0).to_xarray()
        viirs_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_viirs_nc_'+mooring_name+'.csv',skiprows=1,names=['Date','viirs_eppley'],index_col=0).to_xarray()
          
        #Chl Models
        modis_tpca=pd.read_csv('processed/npp_mooring_timeseries/modis_chl_tpca_'+mooring_name+'.csv',skiprows=1,names=['Date','mod_tpca'],index_col=0).to_xarray()
        seawifs_tpca=pd.read_csv('processed/npp_mooring_timeseries/seawifs_chl_tpca_'+mooring_name+'.csv',skiprows=1,names=['Date','sw_tpca'],index_col=0).to_xarray()
        viirs_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/viirs_chlor_a_'+mooring_name+'.csv',skiprows=1,names=['Date','viirs_chlora'],index_col=0).to_xarray()
        modis_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/modis_chlor_a_'+mooring_name+'.csv',skiprows=1,names=['Date','modis_chlora'],index_col=0).to_xarray()
        seawifs_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/seawifs_chlor_a_'+mooring_name+'.csv',skiprows=1,names=['Date','seawifs_chlora'],index_col=0).to_xarray()
        meris_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/meris_chlor_a_'+mooring_name+'.csv',skiprows=1,names=['Date','meris_chlora'],index_col=0).to_xarray()
    
        combined_flux=xr.merge([mod_vgpm,mod_cbpm,mod_eppley,mod_cafe,
                           sw_vgpm,sw_cbpm,sw_eppley,sw_cafe,
                           viirs_vgpm,viirs_cbpm,viirs_eppley,
                           modis_tpca,seawifs_tpca,
                           viirs_chlor_a,modis_chlor_a,seawifs_chlor_a,meris_chlor_a])
        
        combined_flux=combined_flux.assign_coords(Mooring=mooring_name)
        #combined_chl=combined_chl.assign_coords(Mooring=mooring_name)
        
        flux_holder.append(combined_flux)
        #chl_holder.append(combined_chl)
      
    flux=xr.concat(flux_holder,dim='Mooring')    
    flux['Date']=flux.Date.astype('datetime64[D]')
    flux=flux.resample(Date='M').mean()
    flux['Date']=flux.Date.astype('datetime64[M]')
    flux=flux.rename({'sw_vgbm':'sw_vgpm'})
    flux.to_netcdf('processed/flux/npp.nc')
    
    
def add_cafe_and_sst(fp='processed/combined_dataset/month_data_exports.nc'):
    #Bit of a hacky way to modify our data netcdf with the sw and cafe product.

    
    dat=xr.open_mfdataset(fp)
    npp=xr.open_mfdataset('processed/flux/npp.nc')
    dat['Date']=dat['Date'].astype('datetime64[M]')
    npp['Date']=npp['Date'].astype('datetime64[M]')
    npp['Mooring']=dat.Mooring.values
    dat['sw_cafe']=npp.sw_cafe.T
    
    
    cafe=dat[['sw_cafe','mod_cafe']]
    mean = cafe.to_array(dim='new').mean('new')
    dat=dat.assign(cafe=mean)
    
    
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
    
    lonz=[110, 125, 140, 155, 170, 195]
    mooring_sites=['165E','170W','155W','140W','125W','110W']
    
    sst=['datasets/sst/sst.mnmean.nc']
    
    sst=xr.open_mfdataset(sst)
    datslice=pd.DataFrame()
    dats=[]
    for ii,ll in enumerate(lons):
        mooring=mooring_sites[ii]
        oursst1=sst.sel(lat=slice(lats[0],lats[1]),lon=slice(ll[0],ll[1])).mean(dim=['lat','lon']).sst
        oursst1.name=mooring
        dats.append(oursst1)
    
    
    out=xr.concat(dats,dim='Mooring')
    out['Mooring']=lonz[::-1]
    out=out.rename({'time':'Date'})
    
    
    dat['reySST']=out    
    dat.load()
    dat.close()

    os.remove(fp)
    dat.to_netcdf(fp,mode='w',engine='h5netcdf')
    


def create_npp_avgs(): #This function will regrid New production models to the Landschutzer grid.
    landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
    landschutzer=xr.open_dataset(landsch_fp)
    landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')
    land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
    land_pac.to_netcdf('processed/flux/landshutzer.nc',engine='h5netcdf',mode='w')
    land_pac=land_pac.fgco2_smoothed.to_dataset()
    
    
    vgpm=['datasets/npp_satellite/vgpm_mod_nc/*',
            'datasets/npp_satellite/vgbm_sw_nc/*',
            'datasets/npp_satellite/vgpm_viirs_nc/*']
             
    cbpm=['datasets/npp_satellite/cbpm_mod_nc/*',
            'datasets/npp_satellite/cbpm_sw_nc/*',
            'datasets/npp_satellite/cbpm_viirs_nc/*']
            
    eppley=['datasets/npp_satellite/eppley_mod_nc/*',
            'datasets/npp_satellite/eppley_sw_nc/*',  
            'datasets/npp_satellite/eppley_viirs_nc/*']
    cafe=['datasets/npp_satellite/cafe_mod_nc/*','datasets/npp_satellite/cafe_sw_nc/*']
    title=['vgpm','cbpm','eppley','cafe']
    npps=[[],[],[],[]]
    for ix,npp_models in enumerate([vgpm,cbpm,eppley,cafe]):
        model_holder=[]
        for i, model_fp in enumerate(npp_models):
            #if i>=6:
            #    break
            model=xr.open_mfdataset(model_fp,concat_dim='time',combine='nested')
            name=model_fp.split('/')[-2][:-3]
            model=model.rename({'npp':name})
            print(name+' '+str(model.nbytes/1e9)+' GB')
            model_holder.append(model)
        npps[ix]=xr.merge(model_holder)
    
        
        _, index = np.unique(npps[ix]['lon'], return_index=True)
        npps[ix]=npps[ix].isel(lon=index)
        mean = npps[ix].to_array(dim='new').mean('new')
        npp=npps[ix].assign(avg_npp=mean).avg_npp
        
        #USE XESMF to regrid
        #Used Bilinear previously however we need to use conservative method as below
        # regridder = xe.Regridder(npp, land_pac, 'bilinear',reuse_weights=True)
        # npp1d=regridder(npp)
        
        # npp1=npp1d.resample(time='M').mean(dim='time') 
        # npp1['time']=npp1.time.astype('datetime64[M]')#_rg #First day of month
        
        # npp1.to_netcdf('datasets/npp_satellite/avg_npp_rg_'+title[ix]+'.nc')
        
        #Need to add boundary coordinates to perform conservative regridding.
        mod=npp.where(npp>0).to_dataset()
        mod=mod.sortby(['lat'], ascending=True) #This is backwards and was doing everything upside down. 
        lat_rad=land_pac.lat.diff(dim='lat').mean().values/2
        lon_rad=land_pac.lon.diff(dim='lon').mean().values/2
        landlats=np.linspace(land_pac.lat.min().values-lat_rad,land_pac.lat.max().values+lat_rad,len(land_pac.lat)+1)
        landlons=np.linspace(land_pac.lon.min().values-lon_rad,land_pac.lon.max().values+lon_rad,len(land_pac.lon)+1)
        
        lat_rad=mod.lat.diff(dim='lat').mean().values/2
        lon_rad=mod.lon.diff(dim='lon').mean().values/2
        lats=np.linspace(mod.lat.min().values-lat_rad,mod.lat.max().values+lat_rad,len(mod.lat)+1)
        lons=np.linspace(mod.lon.min().values-lon_rad,mod.lon.max().values+lon_rad,len(mod.lon)+1)
        
        x=np.meshgrid(lats,lons)
        #mod=mod.expand_dims({'latb':lats})#,'lonb':lons})
        #mod.coords('lon_b')=lats#=(['lonb','latb'],x[0])
        mod.coords['lon_b']=lons
        mod.coords['lat_b']=lats
        mod['lat_b']=lats#(['lonb','latb'],x[1])
        mod['lon_b']=lons#(['lonb','latb'],x[0])
        mod['lat']=mod.lat
        
        land_pac.coords['lon_b']=landlons
        land_pac.coords['lat_b']=landlats
        land_pac['lat_b']=landlats#(['lonb','latb'],x[1])
        land_pac['lon_b']=landlons##(['lonb','latb'],x[0])
        
        
        # Regridding
        regridder = xe.Regridder(mod, land_pac, 'conservative',reuse_weights=True)
        npp1d=regridder(mod)
        
        npp1=npp1d.resample(time='M').mean(dim='time') 
        npp1['time']=npp1.time.astype('datetime64[M]')#_rg #First day of month
        
        npp1.to_netcdf('processed/flux/avg_npp_rg_'+title[ix]+'.nc',engine='h5netcdf')
    
    
        
    
def convert_tpca_to_month():
    #Process TPCA into monthly so we can regrid it (too big before and this is the resolution we want)
    #sw=xr.open_mfdataset('/g/data/ua8/ocean_color/TPCA_reprocessing/SeaWiFS/*nc',combine='nested',concat_dim='time')# You will need to modify this file path#datasets/tpca/seawifs/*nc')
    sw=xr.open_mfdataset('datasets/chl/tpca/seawifs/*nc',combine='nested',concat_dim='time')
    sw=sw.rename(chl_tpca='sw_tpca')
    sw=sw.resample(time='M').mean(dim='time') 
    sw.to_netcdf('datasets/tpca/sw_month.nc',engine='h5netcdf',mode='w')
    print('saved monthly SW TPCA')

    #mod=xr.open_mfdataset('/g/data/ua8/ocean_color/TPCA_reprocessing/MODIS-Aqua/*nc',combine='nested',concat_dim='time')# Modify this file path #datasets/tpca/modis/*nc')
    mod=xr.open_mfdataset('datasets/chl/tpca/modis/*nc',combine='nested',concat_dim='time')#
    mod=mod.rename(chl_tpca='mod_tpca')
    mod=mod.resample(time='M').mean(dim='time') 
    mod.to_netcdf('datasets/tpca/mod_month.nc',engine='h5netcdf',mode='w')
    print('saved monthly modis TPCA')


def regrid_tpca():
    sw=xr.open_dataset('datasets/tpca/sw_month.nc')
    mod=xr.open_dataset('datasets/tpca/mod_month.nc')
    
    tpca=sw
    tpca=tpca.merge(mod)
    tpca = tpca.to_array(dim='tpca').mean('tpca')
    #tpca.to_netcdf('datasets/tpca/tpca.nc')
    #tpca=tpca.resample(time='M').mean(dim='time') 
    
    landschutzer=xr.open_dataset('processed/flux/landshutzer.nc')
    land_pac=landschutzer.fgco2_smoothed.to_dataset(name='co2')
    land_pac=land_pac.sel(lat=slice(-10,10))


    mod=tpca.to_dataset(name='tpca')
    mod=mod.sortby(['lat'], ascending=True) #This is backwards and was doing everything upside down. 
    lat_rad=land_pac.lat.diff(dim='lat').mean().values/2
    lon_rad=land_pac.lon.diff(dim='lon').mean().values/2
    landlats=np.linspace(land_pac.lat.min().values-lat_rad,land_pac.lat.max().values+lat_rad,len(land_pac.lat)+1)
    landlons=np.linspace(land_pac.lon.min().values-lon_rad,land_pac.lon.max().values+lon_rad,len(land_pac.lon)+1)
    
    lat_rad=mod.lat.diff(dim='lat').mean().values/2
    lon_rad=mod.lon.diff(dim='lon').mean().values/2
    lats=np.linspace(mod.lat.min().values-lat_rad,mod.lat.max().values+lat_rad,len(mod.lat)+1)
    lons=np.linspace(mod.lon.min().values-lon_rad,mod.lon.max().values+lon_rad,len(mod.lon)+1)
    
    mod.coords['lon_b']=lons
    mod.coords['lat_b']=lats
    mod['lat_b']=lats#(['lonb','latb'],x[1])
    mod['lon_b']=lons#(['lonb','latb'],x[0])
    mod['lat']=mod.lat
    
    land_pac.coords['lon_b']=landlons
    land_pac.coords['lat_b']=landlats
    land_pac['lat_b']=landlats#(['lonb','latb'],x[1])
    land_pac['lon_b']=landlons##(['lonb','latb'],x[0])
    print('Regridding TPCA')
    regridder = xe.Regridder(mod, land_pac, 'conservative')
    chl=regridder(mod)
    chl.to_netcdf('datasets/tpca/tpca.nc',engine='h5netcdf',mode='w')
 

def cut_sst_moorings():
    dat = xr.open_dataset('datasets/sst/sst.mnmean.nc')

    lns=[165,190,205,220,235,250]
    ln_names=[110, 125, 140, 155, 170, 195][::-1]
    dat=dat.sel(lat=0.5,method='nearest')
    data=[]
    for i,x in enumerate(lns):
        mooring=dat.sel(lon=x,method='nearest')
        #mooring.sst.plot()
        d=mooring.sst.to_dataframe()
        d=d.drop(columns=['lat','lon'])
        d.columns.name = str(ln_names[i])
        data.append(d.to_xarray())
        
    d=xr.concat(data,dim=ln_names)
    d=d.rename({'concat_dim':'Mooring'})
    d.to_netcdf('processed/indexes/sst.nc',engine='h5netcdf',mode='w')
    return True

def make_earth_grid_m2():
    boxlo,boxla=np.array(np.meshgrid(np.arange(-179.5,179.5,1),np.arange(-89.5,89.5,1)))
    actual_grid=np.cos(np.radians(abs(boxla)))*(111.1*111.1*1000*1000)
    grid_nc=xr.DataArray(actual_grid,coords={'lat':boxla[:,1],'lon':boxlo[1,:]},dims=['lat','lon'])
    lat_size=110567 #in m
    grid_nc['m2']=grid_nc#*lat_size
    grid_nc=grid_nc['m2']
    grid_nc.to_netcdf('processed/earth_m2.nc',engine='h5netcdf',mode='w')
    return True
    
	
def carbon_uatm_to_grams(plotter=0):
    
	co2=xr.open_dataset('processed/flux/landshutzer.nc')
	co2['time']=co2.time.astype('datetime64[M]')
	ratios=xr.open_mfdataset('processed/flux/fratios.nc').laws2011b
	#ratio=f_ratios.laws2011a #laws2000,laws2011a,laws2011b,henson2011

	npp=(xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc').avg_npp/1000*365)

	sst = xr.open_dataset('datasets/sst/sst.mnmean.nc')
	sst= sst.assign_coords(lon=(sst.lon % 360)).roll(lon=(sst.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
	sst=sst.sel(lon=slice(120,290),lat=slice(20,-20)).sst
	sst=sst.sel(time=slice(npp.time.min().values,npp.time.max().values))

	sst=sst.reindex(lat=sst.lat[::-1]) #Lat indexes are backwards


	co2=co2.sel(time=slice(npp.time.min().values,npp.time.max().values))
	sst=sst.sel(time=slice(co2.time.min().values,co2.time.max().values))


	rolling_months=13
	pco2=co2.spco2_smoothed
	co2g=co2.spco2_smoothed#moles_to_carbon(co2.spco2_smoothed)
	co2_rolling=pco2.rolling(time=rolling_months,center=True).mean(dim='time')
	sst_rolling=sst.rolling(time=rolling_months,center=True).mean(dim='time')


	#Equation 1
	pco2Tmean=pco2* np.exp(0.0423*(sst_rolling-sst)) #Equation 1
	pco2Tmean_min=pco2Tmean.rolling(time=rolling_months,center=True).min(dim='time')
	pco2Tmean_max=pco2Tmean.rolling(time=rolling_months,center=True).max(dim='time')

	#Calculate into mg /m3 here.

	#Equation 2
	pco2Tobs=co2_rolling*np.exp(0.0423*(sst_rolling-sst))
	pco2Tobs_min=pco2Tobs.rolling(time=rolling_months,center=True).min(dim='time')
	pco2Tobs_max=pco2Tobs.rolling(time=rolling_months,center=True).max(dim='time')

	#Equation 3
	deltapCO2bio=pco2Tmean_max-pco2Tmean_min
	
        #if plotter==1:
        #    deltapCO2bio.mean(dim='time').plot(),plt.title('deltapco2bio'),plt.show()
	#Equation 4
	deltapCO2temp=pco2Tobs_max-pco2Tobs_min
        
        #if plotter==1:
        #    deltapCO2temp.mean(dim='time').plot(),plt.title('deltapco2temp'),plt.show()
	#Equation 5
	diff=(deltapCO2temp-deltapCO2bio)
	ratio=(deltapCO2temp/deltapCO2bio)

	iida=xr.open_dataset('processed/flux/jma_flux.nc')
	iida=iida.sel(lon=slice(120,290),lat=slice(-20,20))
	iida=iida.sel(time=slice(co2.time.min().values,co2.time.max().values))
	dic=iida.dic.values

	pco2monthvals=deltapCO2bio.values
	pco2monthvals=pco2Tmean.values
	seasurf=sst.values
	flat_pco2=pco2monthvals.reshape(pco2monthvals.shape[0]*pco2monthvals.shape[1]*pco2monthvals.shape[2])
	flat_sst=seasurf.reshape(pco2monthvals.shape[0]*pco2monthvals.shape[1]*pco2monthvals.shape[2])
	flat_dic=dic.reshape(pco2monthvals.shape[0]*pco2monthvals.shape[1]*pco2monthvals.shape[2])
	volumes=np.array([])

	#Flatten and rebuild the array with bespoke sst and dic. 

	#Annoying that I have to calculate it this way. Obviously quite slow but allows bespoke pixel calcultaion. 
	for i,dat in enumerate(flat_pco2):
	    vol=cb.CBsys(DIC=flat_dic[i],pCO2=flat_pco2[i],T_in=flat_sst[i])#T_in=28)#
	    co=(vol.CO2/1000000)*1000
	    #if ~np.isnan(co):
	    #    print(co)
	    #print(i,co)
	    volumes=np.append(volumes,co)
	answer=xr.DataArray(volumes.reshape((pco2monthvals.shape[0],pco2monthvals.shape[1],pco2monthvals.shape[2])),coords=deltapCO2bio.coords)
	answer.to_netcdf('processed/flux/pco2grams.nc',engine='h5netcdf',mode='w')
	return True



#Calculate different f/ep ratios
def calculate_exports_add_to_mooring():
    fp='processed/combined_dataset/month_data.nc'
    dat=xr.open_mfdataset(fp)
    
    ds=dat[['modis_tpca','meris_chlor_a','seawifs_tpca']]
    mean = ds.to_array(dim='new').mean('new')
    ds=ds.assign(chl=mean)
    dat=dat.assign(chl=mean)
    
    ds=dat[['mod_vgpm','sw_vgpm','viirs_vgpm']]
    vgpmmean = ds.to_array(dim='new').mean('new')
    dat=dat.assign(vgpmmean=vgpmmean)
    
    ds1=dat[['mod_cbpm','sw_cbpm','viirs_cbpm']]
    cbpmmean = ds1.to_array(dim='new').mean('new')
    dat=dat.assign(cbpmmean=cbpmmean)
    
    
    
    dat['Date']=dat.Date.astype('datetime64[M]')
    
    sst=xr.open_mfdataset('processed/indexes/sst.nc')
    sst=sst.rename({'time':'Date'})
    sst['Date']=sst.Date.astype('datetime64[M]')
    dat= dat.assign(sst_rey=sst.sst)
    #lns=[165,190,205,220,235,250]
    
    
    
        
    #HAD A BUG HERE. Basically, dat.avgnpp was used to calculate f-ratios. Which is wrong.
    #We are Using CAFE so it must all be CAFE. This should fix it.
    #If changing the npp model used, modifications will need to be made here.
    
    #dat=ds
    zeu1= 34*dat.chl**-0.39#lee 2007
    zeu2=38*dat.chl**-0.428#Morel 1989
    pe_dunne=-0.0101*dat.sst_rey+0.0582*np.log(dat.cafe/zeu1)+0.419
    pe_dunne2=-0.0101*dat.sst_rey+0.0582*np.log(dat.cafe/zeu2)+0.419
    pe_dunne3=0.0081*dat.sst_rey+0.0668*np.log(dat.chl/zeu2)+0.426
    f_ratio=(0.62-(0.02*dat.sst_rey))
    th_e_ratio=(0.23*np.exp(-0.08*dat.sst_rey))
    laws2011a=((0.5857-0.0165*dat.sst_rey)*dat.cafe)/(51.7+dat.cafe)
    
    laws2011b=0.04756*(0.78-((0.43*dat.sst_rey)/30))*dat.cafe**0.307 #avgnpp
    #laws2011b1=0.04756*(0.78-((0.43*dat.sst_rey)/30))*dat.cafe**0.307
    
    dat= dat.assign(zeu_lee=zeu1)
    dat= dat.assign(zeu_morel89=zeu2)
    dat= dat.assign(dunne_zeu1=pe_dunne)
    dat= dat.assign(dunne_zeu2=pe_dunne2)
    dat= dat.assign(f_ratio=f_ratio)
    dat= dat.assign(thE_ratio=th_e_ratio)
    dat= dat.assign(laws2011a=laws2011a)
    dat= dat.assign(laws2011b=laws2011b)
    #dat= dat.assign(laws2011b_vgpm=laws2011b1)
    
    dat.to_netcdf('processed/combined_dataset/month_data_exports.nc',engine='h5netcdf',mode='w')
    
    #pe_dunne.plot()
    #plt.suptitle('Dunne 2005')
    #plt.show()
    
    #pe_dunne2.plot()
    #plt.suptitle('Dunne2 2005')
    #plt.show()
    #pe_dunne3.plot()
    #plt.suptitle('Dunne3_chl 2005')
    #plt.show()
    #f_ratio.plot(cmap='viridis',vmin=0)
    #plt.suptitle('Laws 2000')
    #plt.show()
    #laws2011a.plot()
    #plt.suptitle('laws2011a')
    #plt.show()
    #laws2011b.plot()
    #plt.suptitle('laws2011b')
    #plt.show()
    #plt.contourf(laws2011b1.Date.astype(np.datetime64).values,laws2011b1.Mooring,laws2011b1,levels=np.arange(0.05,0.25,0.025))
    #plt.colorbar()
    #plt.suptitle('laws2011b_vgpm')
    #plt.show()
    
    #th_e_ratio.plot()
    #plt.suptitle('Henson 2011')
    #plt.show()



# %% Prepare Chl data to calc Euz
def calc_euc():
    #this was calculated but hasn't really been used. This file is included in processed/flux/zeu.nc but isn't used in any cases.
    landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
    landschutzer=xr.open_dataset(landsch_fp)
    landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
    land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
    #land_pac.to_netcdf('processed/fluxmaps/landshutzer.nc')
    land_pac=land_pac.fgco2_smoothed
    
    
    chl_mod=xr.open_mfdataset('datasets/tpca/modis/*.nc',combine='nested')
    chl_sw=xr.open_mfdataset('datasets/tpca/seawifs/*.nc',combine='nested')
    #d=xr.concat([chl_sw,chl_mod],dim='model')
    #d['model']=['sw','mod']
    #tpca_chl = d.mean(dim='model')
    #month_tpca_chl=tpca_chl.resample({'time':'M'}).mean()
    #month_tpca_chl['time']=month_tpca_chl.time.astype('datetime64[M]')
    # d.to_dataset().to_array(dim='new').mean('new').mean('model')
    
    #zeu= 34*month_tpca_chl.chl_tpca**-0.39#lee 2007
    zeu= 34*chl_mod.chl_tpca**-0.39#lee 2007
    #USE XESMF to regrid
    regridder = xe.Regridder(zeu, land_pac, 'bilinear',reuse_weights=True)
    zeu1d=regridder(zeu)
    
    #Make a new variable of the average NPP. 
    zeu1=zeu1d.resample(time='M').mean(dim='time') 
    zeu1['time']=zeu1.time.astype('datetime64[M]')#_rg #First day of month
    zeu1.to_netcdf('processed/flux/zeu_mod.nc')
    
    #zeu= 34*month_tpca_chl.chl_tpca**-0.39#lee 2007
    zeu= 34*chl_sw.chl_tpca**-0.39#lee 2007
    #USE XESMF to regrid
    regridder = xe.Regridder(zeu, land_pac, 'bilinear',reuse_weights=True)
    zeu1d=regridder(zeu)
    
    #Make a new variable of the average NPP. 
    zeu1=zeu1d.resample(time='M').mean(dim='time') 
    zeu1['time']=zeu1.time.astype('datetime64[M]')#_rg #First day of month
    zeu1.to_netcdf('processed/flux/zeu_sw.nc')
    
    zeu_mod=xr.open_dataset('processed/flux/zeu_mod.nc')
    zeu_sw=xr.open_dataset('processed/flux/zeu_sw.nc')
    zeu=xr.concat([zeu_sw.chl_tpca,zeu_mod.chl_tpca],dim=['mod','sw']).mean(dim='concat_dim')
    zeu.to_netcdf('processed/flux/zeu.nc')


def make_fratio_nc():
    def convert_trim_fratios():
        '''
        This data is sourced from Tim DeVris, 2017
        https://tdevries.eri.ucsb.edu/models-and-data-products/
        DeVries, T., and Weber, T. (2017). The export and fate of organic matter in the ocean: New constraints from combining satellite and oceanographic tracer observations: EXPORT AND FATE OF MARINE ORGANIC MATTER. Global Biogeochem. Cycles 31, 535–555.
    
        This function will return the regridded eqpac TRIM ef ratio.
    
        '''
        trim=xr.open_dataset('datasets/SIMPLE_TRIM_output.nc')
        
        #ratio of sinking particle flux at the base of the euphotic zone to the NPP at each grid point)
        # This loops through to see the difference between each
        # Tim said that he just averages the 
        # for i in range(0,12):
        #     ver=trim.sel(version=i)
        #     efratio1=ver.NPP/ver.FPOCex
        #     efratio2=ver.FPOCex/ver.NPP
        #     efratio2.T.plot.contourf(levels=np.arange(0,0.425,0.025),cmap='viridis')
        #     plt.suptitle(str(i))
        #     plt.show()
        #     print(i)
        #     print((((ver.FPOCex/1000)*12)*ver.Area).sum().values/1e15) #Global integrated carbon removal
        # #trim=trim.set_coords('version')
        
        ratio=(trim.FPOCex/trim.NPP).mean(dim='version')
        ratio.name='avg'
        std=(trim.FPOCex/trim.NPP).std(dim='version')
        std.name='stdev'
        ratio=xr.merge([ratio,std])
        ourset=trim.mean(dim='version')
        print((((ourset.FPOCex/1000)*12)*ourset.Area).sum().values/1e15) #Global integrated carbon removal
        #sel(version=5)
        ratio['latitude']= trim.LAT.values[0][0]
        ratio['longitude']= trim.LON.values[0][:,0]
        
        ratio=ratio.rename({'latitude':'lat','longitude':'lon'})
        ratio=ratio.where((ratio>0)&(ratio<100))
        
        eqpac_trim=ratio.sel(lat=slice(-21,21),lon=slice(119,301))
        eqpac_trim['avg']=eqpac_trim.avg.T
        eqpac_trim['stdev']=eqpac_trim.stdev.T
        
        #eqpac_trim.stdev.plot.contourf(levels=np.arange(0,0.425,0.025))
         
        landschutzer=xr.open_dataset('processed/flux/landshutzer.nc')
        land_pac=landschutzer.fgco2_smoothed
        land_pac=land_pac.sel(lat=slice(-20,20))
        
        regridder = xe.Regridder(eqpac_trim, land_pac, 'bilinear')
        eqpac_trim_grid=regridder(eqpac_trim)
        return eqpac_trim_grid

    landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
    landschutzer=xr.open_dataset(landsch_fp)
    landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
    land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
    #land_pac.to_netcdf('processed/fluxmaps/landshutzer.nc')
    land_pac=moles_to_carbon(land_pac.fgco2_smoothed)/365
    
    npp=xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc') #MAKE SURE THIS IS CORRECT MODEL
    zeu=xr.open_dataset('processed/flux/zeu.nc').chl_tpca
    

    sst = xr.open_dataset('datasets/sst/sst.mnmean.nc')
    sst= sst.assign_coords(lon=(sst.lon % 360)).roll(lon=(sst.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
    sst=sst.sel(lon=slice(120,290),lat=slice(20,-20))
    sst=sst.sel(time=slice(npp.time.min().values,npp.time.max().values))
    land_pac=land_pac.sel(time=slice(npp.time.min().values,npp.time.max().values))
    land_pac['time']=land_pac.time.astype('datetime64[M]')
    
    #Make sure using the correct npp model to calculate these.
    laws2011a=((0.5857-0.0165*sst.sst)*npp.avg_npp)/(51.7+npp.avg_npp)
    laws2011b=0.04756*(0.78-((0.43*sst.sst)/30))*npp.avg_npp**0.307
    laws2000=(0.62-(0.02*sst.sst))
    henson2011=(0.23*np.exp(-0.08*sst.sst))
    pe_dunne=-0.0101*sst.sst+0.0582*np.log((npp.avg_npp/12)/zeu)+0.419
    pe_dunne.name='dunne2005'
    laws2011a.name='laws2011a'
    laws2011b.name='laws2011b'
    laws2000.name='laws2000'
    henson2011.name='henson2011'
    trim=convert_trim_fratios() #Calculate devris TRIM model.
    trim_av=trim.avg
    trim_std=trim.stdev
    trim_av.name='trim'
    trim_std.name='trim_std'
    
    ratios=xr.merge([laws2011a,laws2011b,laws2000,henson2011,pe_dunne,trim_av,trim_std])
    ratios=ratios.where((ratios>0)&(ratios<100))
    ratios.to_netcdf('processed/flux/fratios.nc',mode='w',engine='h5netcdf')
    
    #grid=xr.open_dataarray('processed/earth_size.nc')
    
    #((npp.avg_npp*laws2011a).mean(dim='time')*grid.m2).plot()
    #plt.show()
    
    
def save_landschutzer_2018_seamask():
    seamask=xr.open_dataset('datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc').seamask
    seamask.to_netcdf('processed/seamask.nc')
# %% RUN FUNCS HERE.

# print('Regridding the NPP models')
# #create_npp_avgs() #Regrid the NPP models.
# print("Running TPCA to month calc")
# convert_tpca_to_month()
# print('Regridding TPCA - xesmf')
# regrid_tpca()
# print('Calculate when ENSOs occured +- 0.5 MEI')
# find_enso_events()
# print('Convert our primary productivity moorings to netcdf')
# npp_csvs_to_nc()
# #combine_csvs_to_nc() #This one should already have been run. 
print('Working out the SST for each mooring')
cut_sst_moorings()
# print('Calculate earth size per pixel')
# make_earth_grid_m2()
# print('Convert uatm carbon to grams of carbon')
# #carbon_uatm_to_grams(plotter=0) #This one uses lots of memory and time. might need to qsub it.
print('Calculate f-ratio maps')
make_fratio_nc()
print('Adding Cafe and SST')
add_cafe_and_sst(fp='processed/combined_dataset/month_data.nc')
print('Combining export data to processed/combined_dataset/month_data_exports.nc')

calculate_exports_add_to_mooring()

print('Save landschutzer 2018 seamask')
save_landschutzer_2018_seamask() #Just to make sure that it stays in the folder as the 2020 version doesnt have this.

#print('Calculating euphotic depth from chl - lee 2007')
#calc_euc() # This one isn't actually used in any of the analysis. Uses too much memory anyway. File is processed and stored in processed/flux/zeu.nc if desired.
