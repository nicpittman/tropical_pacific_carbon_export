#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:52:00 2020

@author: Nic Pittman

Ewp
Takes a little while to run.

This code is a pretty average (poorly named variables and reuse of dat, and different names for the moorings).
All of the values are calculated on the fly and printed. Not saved anywhere. 
Quite inefficient
Recommended to turn warnings off to save the data and put in text. 
Could have turned this into a function. Have refractored where possible but this script is provided as is.

Requirements:
    processed/combined_dataset/month_data_exports.nc
    processed/flux/pco2grams.nc
    
Produces:
    figs/Figure5a_ENSO_seasonality.png
    
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
from matplotlib.dates import MonthLocator, DateFormatter 
from matplotlib.ticker import FuncFormatter

#%% Seasonal Comparison  
#Air sea flux minus new production

#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

plt.figure(figsize=(12,16))

month_fmt = DateFormatter('%b')
def m_fmt(x, pos=None):
    return month_fmt(x)#[0]


ax = plt.subplot(6,2,5)
cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

lns=[165,190,205,220,235,250]
moors=[110, 125, 140, 155, 170, 195]
moorings=['110W','125W','140W','155W','170W','165E'][::-1]


fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    asfdiff=((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((asfdiff.month-1).astype('datetime64[M]'),asfdiff,c=cols[i],linewidth=3,label=moorings[i])
  
plt.title('e) Air-Sea Flux - New Production seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()

ax = plt.subplot(6,2,6)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()
for i, mooring_name in enumerate(moorings):

    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
   
    dat['Date'].astype('datetime64[M]')
    
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.select_model)]
    
    nino1=info[info.mei>0.5].mean()
    nina1=info[info.mei<-0.5].mean()
    neutral1=info[(info.mei<0.5)&(info.mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1.co2-nino1.select_model,
        'La Nina':nina1.co2-nina1.select_model,
        'Neutral':neutral1.co2-neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':nino1.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'El Nino NP':nino1.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name)
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()

plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('f)  Air-Sea Flux - New Production ENSO',loc='left')
print(final_mooring_enso_avgs)
ENSOAVG=final_mooring_enso_avgs


# %% New Production


ax = plt.subplot(6,2,1)
cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])

plt.title('a) New Production seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()


ax = plt.subplot(6,2,2)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()
for i, mooring_name in enumerate(moorings):
    print(mooring_name)
   
    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
    
    dat['Date'].astype('datetime64[M]')
    
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.co2)]
    
    nino1=info[info.mei>0.5].mean()
    nina1=info[info.mei<-0.5].mean()
    neutral1=info[(info.mei<0.5)&(info.mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1.select_model,
        'La Nina':nina1.select_model,
        'Neutral':neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':nino1.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'El Nino NP':nino1.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name)
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()

plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('b)  New Production ENSO',loc='left')

# %%    PCO2 section
   
plt.subplot(627)


cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4']#[::-1]


pco2_month = xr.open_dataarray('processed/flux/pco2grams.nc')
pco2_month = (pco2_month*12*50)#/30
for i, mooring_name in enumerate(moors):
    
    pco2_m=pco2_month.sel(lat=0,lon=lns[i],method='nearest')
    
    
    pco2=pco2_m.groupby('time.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(pco2_m.groupby('time.month').mean()-pco2_m.groupby('time.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])

plt.title('g) pCO2t seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$')
plt.grid()

#ENSO BREAKDOWN
plt.subplot(628)
#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()


pco2_month=pco2_month.rename({'time':'Date'})
lns=[165,190,205,220,235,250]
moors=[110, 125, 140, 155, 170, 195]
moorings=['110W','125W','140W','155W','170W','165E'][::-1]
for i, mooring_name in enumerate(moors):
    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    #try:
    dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name))
    #except:
    #    dat=xr.open_mfdataset(fp).sel(Mooring=195)
   
    pco2_m=pco2_month.sel(lat=0,lon=lns[i],method='nearest')
    dat['pco2t']=pco2_m
#for i, mooring_name in enumerate(moorings):
    print(mooring_name)
  
    dat['Date'].astype('datetime64[M]')
    
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.select_model)]
    
    nino1=info[info.mei>0.5].mean()
    nina1=info[info.mei<-0.5].mean()
    neutral1=info[(info.mei<0.5)&(info.mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1.co2-nino1.select_model,
        'La Nina':nina1.co2-nina1.select_model,
        'Neutral':neutral1.co2-neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':nino1.pco2t,
        'La Nina CO2':nina1.pco2t,
        'Neutral CO2':neutral1.pco2t,
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=moorings[i])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=moorings[i])
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso_avgs.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()
plt.ylabel('gC m$^{-2}$')
plt.title('h) pCO2t ENSO',loc='left')



# %% How about the SST
plt.subplot(629)


cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=xx.sst_rey.groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])


plt.title('i) SST seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
#plt.axhline(0,c='gray')
#plt.text(np.datetime64('1970-07-15'),0.025,'Air-Sea flux Dominated',fontsize=11)
#plt.text(np.datetime64('1970-07-15'),-0.02,'Biology Dominated',fontsize=11)
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.xlabel('Month')
plt.ylabel('Degree C')
plt.grid()

#ENSO BREAKDOWN


plt.subplot(6,2,10)
#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()

ty='month' #Actually month though need to fix this.
fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp).sel(Mooring=int(110))
dat['Date'].astype('datetime64[M]')
dat=dat.rename({'Date':'time'})


moors=[110, 125, 140, 155, 170, 195][::-1]

for i, mooring_name in enumerate(lns):
    fp='processed/combined_dataset/month_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(moors[i]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)

    d=dat.sel(Date=slice(pco2_month.Date.min().values,pco2_month.Date.max().values))
    
    pco2_m=d.sst_rey#seasonaltrend_sstobs[i]
    mei=d.mei
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring


    nino1=pco2_m[mei>0.5].mean()
    nina1=pco2_m[mei<-0.5].mean()
    neutral1=pco2_m[(mei<0.5)&(mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1,
        'La Nina':nina1,
        'Neutral':neutral1
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    
    ensoavgs=pd.Series(x,name=moorings[i])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()

plt.xlabel('Mooring')
plt.ylabel('Degree C')
plt.title('j) SST ENSO',loc='left')



# %% Air Sea Flux

ax = plt.subplot(6,2,3)

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.co2flux4_land_gmyr/365).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])

   
   
plt.title('c) Air-Sea Flux seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()


ax = plt.subplot(6,2,4)


means=pd.DataFrame()

final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()
for i, mooring_name in enumerate(moorings):
    print(mooring_name)

    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
   
    dat['Date'].astype('datetime64[M]')
    
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.co2)]
    
    nino1=info[info.mei>0.5].mean()
    nina1=info[info.mei<-0.5].mean()
    neutral1=info[(info.mei<0.5)&(info.mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1.co2,
        'La Nina':nina1.co2,
        'Neutral':neutral1.co2
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':nino1.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'El Nino NP':nino1.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name)
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()

plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('d)  Air-Sea Flux ENSO',loc='left')


#%% How about the windspeed

plt.subplot(6,2,11)

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.windspeed).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[::-1][i])


plt.title('k) Windspeed seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
plt.legend(ncol=2)
plt.xlabel('Month')
plt.ylabel('Windspeed (m/s)')
plt.grid()

#ENSO BREAKDOWN

plt.subplot(6,2,12)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()

ty='month' #Actually month though need to fix this.
fp='processed/combined_dataset/month_data_exports.nc'

lns=[165,190,205,220,235,250]
for i, mooring_name in enumerate(lns):
    fp='processed/combined_dataset/month_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(moors[i]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
        
    #plt.subplot(6,1,i+1)
    #pco2_m=pco2_month.sel(lat=0,lon=mooring_name,method='nearest')

    d=dat.sel(Date=slice(pco2_month.Date.min().values,pco2_month.Date.max().values))
    
    pco2_m=d.windspeed#seasonaltrend_sstobs[i]
    mei=d.mei
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring


    nino1=pco2_m[mei>0.5].mean()
    nina1=pco2_m[mei<-0.5].mean()
    neutral1=pco2_m[(mei<0.5)&(mei>-0.5)].mean()
    
    #pd.Serie
    x ={'El Nino':nino1,
        'La Nina':nina1,
        'Neutral':neutral1
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    
    ensoavgs=pd.Series(x,name=moorings[i])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    if 'Neutral' in x[0]:
        c='darkgoldenrod'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
            
    plt.plot(x[1].index,x[1],c=c,label=x[0],linewidth=3)
plt.grid()
plt.legend(final_mooring_enso.columns)

plt.xlabel('Mooring')
plt.ylabel('Windspeeed (m/s)')
plt.title('l) Windspeed ENSO',loc='left')
plt.tight_layout()

plt.savefig('figs/Figure5a_ENSO_seasonality.png',dpi=200)
plt.show()

print(ENSOAVG.mean())
