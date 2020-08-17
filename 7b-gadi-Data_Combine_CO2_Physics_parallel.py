#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 14:10:46 2020

Run through + compile data for NP ch2.
Datasets:
    TaoTriton CO2 / pco2 / delta pco2 / xco2?
    TaoTriton Physics SST, Windspeed, 
    Satellite Chl
    Primary Productivity
    Climate Indicies: MEI, SOI, PDO

Compiles into one at 1:3 hour intervals. 
8 Day NPP where possible
Climate indicies monthly

some ocean physics missing still.

Moorings across pacific Ocean. 

Steps through windspeed as in daily files, easiest. 
Keeps hourly pco2 data
8 Day NPP
@author: nicpittman
"""

import xarray as xr
import glob
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pandas as pd
from carbon_math import carbon_flux

#start_time = datetime.now()

def moles_to_carbon(moles):
    #In: moles Carbon/Co2
    #Out: g/C

    molarmassC=12
    grams=moles*molarmassC
    return grams

def nearest_ind(items, pivot,label=None):
    try:
        find=items.NPP
    except:
        try:
            find=items.chl
        except:
            find=pd.DataFrame(eval('items.'+str(label)).values)

    item_date=items.Date
   # print(item_date)
   # print(type(item_date))
    try:
        time_diff = (np.abs([(pd.to_datetime(date) - pivot).days for date in item_date.values])).flatten()
    except:
        #Weird label and cut the start of the array:
        time_diff = (np.abs([(pd.to_datetime(date) - pivot).days for date in item_date.values[1:]])).flatten()


    if time_diff.min()<=5:
        return find.iloc[time_diff.argmin(0)]
    else:
        return np.nan

#Tao 155W
#Loading and combing physical and chemical variables for each timeseries (Starting with 155W.)
#CO2

moorings=['110W','125W','140W','155W','170W','165E']

#for mooring_name in moorings:
def process(mooring_name):
    print("Starting: "+mooring_name)
    #for mooring_name in moorings:
    #def process(mooring_name):
    print("Starting: "+mooring_name)
    #for mooring_name in moorings:
    if mooring_name == '155W':
         #Preprocessed this one check other scripts
        mooring_path='datasets/tao/tao_co2/TAO'+mooring_name+'_full.csv'
    else:
        #General filepath
        mooring_path='datasets/tao/tao_co2/TAO'+mooring_name+'.csv'
    
    mooring=pd.read_csv(mooring_path,comment='#') #Hourly to 3 hourly timesteps
    mooring['Date'] = pd.to_datetime(mooring['Date'])  
    
    #Climate Indicies
    soifp='datasets/indexes/soi.csv'
    pdofp='datasets/indexes/pdo.csv'
    meifp='datasets/indexes/meiv2.csv'
    emifp='datasets/indexes/emi.csv'
    soi=pd.read_csv(soifp)
    pdo=pd.read_csv(pdofp,header=1)
    mei=pd.read_csv(meifp)# If using original (not v2) need to use this. ,header=1)
    emi=pd.read_csv(emifp)
    #Primary Productivity Models and Time Series
    mod_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgpm_mod_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    mod_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_mod_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    mod_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_mod_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    mod_cafe=pd.read_csv('processed/npp_mooring_timeseries/cafe_mod_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    sw_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgbm_sw_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    sw_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_sw_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    try:
        sw_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_sw_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    except:
        sw_eppley=np.nan
    viirs_vgpm=pd.read_csv('processed/npp_mooring_timeseries/vgpm_viirs_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    viirs_cbpm=pd.read_csv('processed/npp_mooring_timeseries/cbpm_viirs_nc_'+mooring_name+'.csv',names=['Date','NPP'])
    viirs_eppley=pd.read_csv('processed/npp_mooring_timeseries/eppley_viirs_nc_'+mooring_name+'.csv',names=['Date','NPP'])
   
    #Chl Models
    modis_tpca=pd.read_csv('processed/npp_mooring_timeseries/modis_chl_tpca_'+mooring_name+'.csv',names=['Date','chl'])
    seawifs_tpca=pd.read_csv('processed/npp_mooring_timeseries/seawifs_chl_tpca_'+mooring_name+'.csv',names=['Date','chl'])
    viirs_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/viirs_chlor_a_'+mooring_name+'.csv',names=['Date','chl'])
    modis_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/modis_chlor_a_'+mooring_name+'.csv',names=['Date','chl'])
    seawifs_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/seawifs_chlor_a_'+mooring_name+'.csv',names=['Date','chl'])
    meris_chlor_a=pd.read_csv('processed/npp_mooring_timeseries/meris_chlor_a_'+mooring_name+'.csv',names=['Date','chl'])
    
    #Physical Variables
    start_day=mooring.iloc[0].Date-np.timedelta64(3,'M') #Well start minus three months stops needless processing for some of the moorings
    windspeed=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').WS_401.sel(time=slice(start_day,'2021'))
    winddir=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').WD_410.sel(time=slice(start_day,'2021'))
    thermocline=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/iso0n'+mooring_name.lower()+'_dy.cdf').ISO_6
    airtemp=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').AT_21.sel(time=slice(start_day,'2021'))
    sst2=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').T_25.sel(time=slice(start_day,'2021'))
    #Silence warning: 
    # RuntimeWarning: invalid value encountered in greater
    # sst2=sst2.where(sst2.values>0)
    with np.errstate(invalid='ignore'):
        sst2=sst2.where(sst2.values>0) #Quality control above 0C
        
    wu=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').WU_422.sel(time=slice(start_day,'2021'))
    wv=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/met0n'+mooring_name.lower()+'_dy.cdf').WV_423.sel(time=slice(start_day,'2021'))
    sss2=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/sss0n'+mooring_name.lower()+'_dy.cdf').S_41.sel(time=slice(start_day,'2021'))
    ssd= xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/ssd0n'+mooring_name.lower()+'_dy.cdf').STH_71.sel(time=slice(start_day,'2021'))
    lat=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/pos0n'+mooring_name.lower()+'_dy.cdf')['LAT_500']
    lon=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/pos0n'+mooring_name.lower()+'_dy.cdf')['LON_502']
    precip=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/rain0n'+mooring_name.lower()+'_dy.cdf').RN_485.sel(time=slice(start_day,'2021'))
    #Still want to use this but shape makes it hard (~30 depths - probably ok to remain standalone)
    temps=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/t0n'+mooring_name.lower()+'_dy.cdf')
    

    #Open Flux Maps
    co2flux_JMA=xr.open_mfdataset('processed/flux/JMA_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
    
    co2flux_landshutz=xr.open_dataset('processed/flux/landsch_mooring_co2_flux.nc').sel(Mooring=mooring_name)
    #co2flux_yasanaka=xr.open_dataset('processed/flux/yasanaka_mooring_co2_flux.nc').sel(Mooring=mooring_name)
    #Dates are set as the 15th of the month, but with the backfill d1 method we need to change that to the first day of each month.
    co2flux_landshutz['time']=co2flux_landshutz.time-np.timedelta64(14,'D') 
    
    #temps_monthly=temps.resample(time='M').mean()
    #temps_selected=temps_monthly.T_20.sel(depth=slice(0,300))
    #temps_selected=temps_selected.interpolate_na(dim='depth',method='linear')
    #temps_selected=temps_selected.interpolate_na(dim='time',method='linear').T
    #temps_selected['depth']=temps_selected['depth']*-1
    
    #plt.figure(figsize=(20,4))
    #x=plt.contourf(temps_selected.time.values.astype(np.datetime64),temps_selected.depth.values,temps_selected.squeeze(),levels=np.arange(10,32,3),cmap='RdYlBu_r')
    #plt.colorbar(x)
    #plt.contour(temps_selected.time.values.astype(np.datetime64),temps_selected.depth.values,temps_selected.squeeze(),levels=[20],colors='k',linewidths=0.5)
    #plt.title(mooring_name+': Upper ocean temperatures')
    #plt.ylim([-250,-0])
    #plt.xlim([np.datetime64('1998'),np.datetime64('2020')])
    #plt.show()
    
    #Setup
    final_mooring=pd.DataFrame()
    
    for i,dt in enumerate(np.arange(np.datetime64('1997-01-01'),
                                    np.datetime64('2020-01-01'))):
        wspd=windspeed[windspeed.time.astype('datetime64[D]')==dt] 
#    for i,wspd in enumerate(windspeed): #In daily files.
        #Getting a day cut out and info.
        #d0=np.datetime64(wspd.time.values)-np.timedelta64(12, 'h')
        #d1=d0+np.timedelta64(1, 'D')#-np.timedelta64(12, 'h')
        d0=dt
        d1=d0+np.timedelta64(1, 'D')
        #print(d0,d1)
        month=str(d0)[5:7].lstrip('0')
        year=str(d0)[0:4]
        day=str(d0)[8:10]
        holder=mooring[(mooring.Date>np.datetime64(d0))&(mooring.Date<np.datetime64(d1))]
        holder=holder.dropna(subset=['Date'])

        #if len(holder)>0: #This is where we are losing the extra data.
            #These here will save us a little time.
            #Calculate Physical variables and Climate Modes
        try:#This keeps breaking on the last days.. So lets just end it where it breaks as much as possible.
            holder=holder.copy()
    
            try:
                la=lat.sel(time=slice(d0,d1)).values.flatten()[0] 
                lo=lon.sel(time=slice(d0,d1)).values.flatten()[0] 
            except:
                lo,la=np.nan,np.nan
            if len(holder)==0:
                holder=holder.append(pd.Series({'mooring_lat':la}),ignore_index=True)
                holder['mooring_lon']=lo
                holder['Date']=d0
                
            else:
                holder['mooring_lat']=la#=holder.append({'mooring_lat':la},ignore_index=True)
                holder['mooring_lon']=lo#holder=holder.append({'mooring_lon':lo},ignore_index=True)
            try:
                atemp=airtemp.sel(time=slice(d0,d1)).values.flatten()[0] 
                holder['airtemp']=atemp
            except:
                pass
            try:
                sst2a=sst2.sel(time=slice(d0,d1)).values.flatten()[0] 
                holder['sst2']=sst2a   
            except:
                pass
            try:
                ssd2a=ssd.sel(time=slice(d0,d1)).values.flatten()[0]    
                holder['ssd']=ssd2a
            except:
                pass
            try:
                prec=precip.sel(time=slice(d0,d1)).values.flatten()[0]    
                holder['precip']=prec
            except: #One of these broke in one permutation
                pass
            
            if wspd.values.size!=0:
                ws=wspd.values.flatten()[0]
                wus=wu.sel(time=slice(d0,d1)).values.flatten()[0]
                wvs=wv.sel(time=slice(d0,d1)).values.flatten()[0]
                wdir=winddir.sel(time=slice(d0,d1)).values.flatten()[0] 
                                  #Windspeed
            else:
                ws,wus,wvs,wdir=np.nan,np.nan,np.nan,np.nan
                
            try:
                tcd=thermocline.sel(time=slice(d0,d1)).values.flatten()[0]  #Theromcline Depth
            except:
                tcd=np.nan
                
    
            sss2a=sss2.sel(time=slice(d0,d1)).values.flatten()
            if sss2a.size>0:
                holder['sss2']=sss2a[0]
            else:
                holder['sss2']=np.nan
    
            holder['windspeed']=ws
            holder['thermoclinedepth']=tcd
            holder['winddirection']=wdir
            holder['wu']=wus
            holder['wv']=wvs
     
            try:
                vsoi=soi[soi.Year==int(year)][month].values[0]              #SOI
                holder['soi']=vsoi
            except:
                holder['soi']=np.nan
            try:   
                vmei=mei[mei.Year==int(year)][month].values[0]              #MEI #Migh need YEAR (Caps) for original mei
                holder['mei']=vmei
            except:
                holder['mei']=np.nan
            try:
                vpdo=pdo[pdo.Year==int(year)][month].values[0]              #PDO
                holder['pdo']=vpdo
            except:
                holder['pdo']=np.nan
            try:
                holder['emi']=emi[emi.Date.astype('datetime64[M]')==d0.astype('datetime64[M]')].EMI.values[0]
            except:
                holder['emi']=np.nan
            #Append our desired data into this.        
            holder['mod_vgpm']=nearest_ind(mod_vgpm,d0)                              #Modis VGPM      
            holder['mod_cbpm']=nearest_ind(mod_cbpm,d0)                              #Modis CPBM
            holder['mod_eppley']=nearest_ind(mod_eppley,d0)                          #Modis Eppley
            holder['mod_cafe']=nearest_ind(mod_cafe, d0)                             #Modis Cafe           
            holder['sw_vgpm']=nearest_ind(sw_vgpm,d0)                                #SeaWiFS VGPM
            holder['sw_cbpm']=nearest_ind(sw_cbpm,d0)                                #SeaWiFS CPBM
            holder['sw_eppley']=nearest_ind(sw_eppley,d0)                            #SeaWiFS Eppley         
            holder['viirs_vgpm']=nearest_ind(viirs_vgpm,d0)                          #SeaWiFS VGPM
            holder['viirs_cbpm']=nearest_ind(viirs_cbpm,d0)                          #SeaWiFS CPBM
            holder['viirs_eppley']=nearest_ind(viirs_eppley,d0)                      #SeaWiFS Eppley
            
            holder['modis_tpca']=nearest_ind(modis_tpca,d0)
            holder['seawifs_tpca']=nearest_ind(seawifs_tpca,d0)
            holder['viirs_chlor_a']=nearest_ind(viirs_chlor_a,d0)
            holder['modis_chlor_a']=nearest_ind(modis_chlor_a,d0)
            holder['seawifs_chlor_a']=nearest_ind(seawifs_chlor_a,d0)
            holder['meris_chlor_a']=nearest_ind(meris_chlor_a,d0)
                 
            #Calculate our own CO2 flux estimates
            try:
                holder['co2flux']=carbon_flux(holder.SSS,holder.SST,holder.windspeed,holder.pCO2_sw,holder.pCO2_air)[0]
                holder['co2flux_gmyr']=moles_to_carbon(holder.co2flux)
            except:
                pass
            
            try:
                holder['co2flux2']=carbon_flux(holder.sss2,holder.sst2,holder.windspeed,holder.pCO2_sw,holder.pCO2_air)[0]
                holder['co2flux2_gmyr']=moles_to_carbon(holder.co2flux2)
            except:
                pass
            
            
            holder['avgnpp']=holder[['viirs_eppley','viirs_cbpm','viirs_vgpm','sw_eppley','sw_cbpm','sw_vgpm','mod_cafe','mod_eppley','mod_cbpm','mod_vgpm']].mean(axis=1)
    
            #Use pad and d1 so we get the months as is
            try:
                holder['co2flux3_JMA']=co2flux_JMA.sel(Date=d1,method='pad').flux.values
                holder['pH']=co2flux_JMA.sel(Date=d1,method='pad').pH.values
                holder['dic']=co2flux_JMA.sel(Date=d1,method='pad').dic.values        
                holder['co2flux3_JMA_gmyr']=moles_to_carbon(holder.co2flux3_JMA)
            except:
                holder['co2flux3_JMA']=np.nan
                holder['pH']=np.nan
                holder['dic']=np.nan
                holder['co2flux3_JMA_gmyr']=np.nan
            
            try:
                holder['co2flux4_landshutz']=co2flux_landshutz.sel(time=d1,method='pad').fgco2_raw.values
                holder['co2flux4_land_gmyr']=moles_to_carbon(holder.co2flux4_landshutz)
                holder['pco2atm_landshutz']=co2flux_landshutz.sel(time=d1,method='pad').atm_co2.values
                holder['pco2sw_landshutz']=co2flux_landshutz.sel(time=d1,method='pad').spco2_raw.values
                holder['sol_landshutz']=co2flux_landshutz.sel(time=d1,method='pad').sol.values
                holder['kw_landshutz']=co2flux_landshutz.sel(time=d1,method='pad').kw.values
            except:
                holder['co2flux4_landshutz']=np.nan
                holder['co2flux4_land_gmyr']=np.nan
                holder['pco2atm_landshutz']=np.nan
                holder['pco2sw_landshutz']=np.nan
                holder['sol_landshutz']=np.nan
                holder['kw_landshutz']=np.nan
            
            
            try:
                holder['co2flux5_yasanaka']=co2flux_yasanaka.sel(time=d1,method='pad').flux_masked.values
                holder['co2flux5_yasanaka_gmyr']=moles_to_carbon(holder.co2flux5_yasanaka)
                holder['co2flux5_yasanaka_unmasked']=co2flux_yasanaka.sel(time=d1,method='pad').co2flux.values
                holder['pco2sw_yasanaka']=co2flux_yasanaka.sel(time=d1,method='pad').pco2.values
            except:
                holder['co2flux5_yasanaka']=np.nan
                holder['co2flux5_yasanaka_gmyr']=np.nan
                holder['co2flux5_yasanaka_unmasked']=np.nan
                holder['pco2sw_yasanaka']=np.nan
            
       
            #print(holder)
        
        except KeyboardInterrupt:
            print('interrupted!')
            import sys
            sys.exit()
            break
        except Exception as e:
            print(e)
            pass
        #Bit of printing so we know whats going on.
        final_mooring=final_mooring.append(holder,sort=True)
        if i%50==0:
    
            print(holder.head(1))
            #days_since=str(final_mooring.iloc[-1].Date-pd.datetime(int(year),int(month),int(day))).split(' ')[0]
            #Track time taken
            #print("Cycle at %s seconds." % (str(datetime.now()- start_time)))
            #start_time = datetime.now()
                
        #else: #holder =0 - no co2 data
            #was breaking because some data didn't exist during 2019 but not a problem
            #because the CO2 timeseries doesn't go that long (yet?)
            #try:
            #    print('Broken and skipped. Days since: '+days_since)
            #except:
        #    print('Skipped',d0,d1)
        #    continue

    #And save it!
    #Add filepath magic!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    path='processed/combined_dataset/'+mooring_name+'_combined.csv'
    final_mooring.to_csv(path,index=False)
    print('Finished: '+ mooring_name)    
    print('Saved to: '+path)





#Chuck everything above into a single element function, here mooring loc
#Then this is the parrel code below.
#from joblib import Parallel, delayed
#import multiprocessing
#from multiprocessing import Process
#moorings=['155W','170W','165E']
#num_cores=multiprocessing.cpu_count()-1 #-1 so we don't lock up. Up to 6 cores I guess, one for each mooring.
#res=Parallel(n_jobs=num_cores)(delayed(process)(mooring) for mooring in moorings)





#Chuck everything above into a single element function, here mooring loc
#Then this is the parrel code below.
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import Process
moorings=['110W','125W','140W','155W','170W','165E']
num_cores=6#multiprocessing.cpu_count()-1 #-1 so we don't lock up. Up to 6 cores I guess, one for each mooring.
res=Parallel(n_jobs=num_cores)(delayed(process)(mooring) for mooring in moorings)



#Calculate Day, Week and Monthly averages for the dataset
#Save the data in binary .nc

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
    wavg_a.append(wavg)
    mavg_a.append(mavg)
    
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
all_data.to_netcdf(fp+'all_data.nc',engine='h5netcdf') #May or may not need engine, just that there is a type issue otherwise.
daily.to_netcdf(fp+'day_data.nc')
weekly.to_netcdf(fp+'week_data.nc')
monthly.to_netcdf(fp+'month_data.nc')


#plt.contourf(daily.Mooring.values*-1,
#             daily.Date.astype(np.datetime64).values,
#             daily.co2flux2.interpolate_na(dim='Date').interpolate_na(dim='Mooring').T)

#plt.colorbar()
#plt.show()
#daily.co2flux2.plot()

