#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:52:00 2020

@author: Nic Pittman

Ewp
Takes a little while to run.

This code is a pretty average (poorly named variables and reuse of dat, and different names for the moorings).
All of the values are calculated on the fly and printed. Not saved anywhere. 
Quite inefficient, but provided as is.
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


xl0=0.0
yl0=0.18

xl1=0.0
yl1=0.18

xl2=-0.09
yl2=0


lanina=pd.read_csv('processed/indexes/la_nina_events.csv')
cp_nino=pd.read_csv('processed/indexes/cp_events.csv')
ep_nino=pd.read_csv('processed/indexes/ep_events.csv')


fp='processed/combined_dataset/month_data_exports.nc'
info=xr.open_mfdataset(fp).sel(Mooring=195).to_dataframe()


#Process EP, CP and Nino events.
nina=pd.DataFrame()
ep=pd.DataFrame()
cp=pd.DataFrame()
for i in lanina.iterrows(): nina=nina.append(info[slice(i[1].start,i[1].end)])
for i in ep_nino.iterrows(): ep=ep.append(info[slice(i[1].start,i[1].end)])
for i in cp_nino.iterrows(): cp=cp.append(info[slice(i[1].start,i[1].end)])
nina_dates=nina.index
ep_dates=ep.index
cp_dates=cp.index
    

#%% Seasonal Comparison  
#Air sea flux minus new production

#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

plt.figure(figsize=(12,16))

month_fmt = DateFormatter('%b')
def m_fmt(x, pos=None):
    return month_fmt(x)#[0]


ax = plt.subplot(6,2,6)
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
  
plt.ylim([xl2,yl2])
plt.title('f) Air-sea flux - new production seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()

ax = plt.subplot(6,2,5)

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
    

    
    
   
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
 #   nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    
    #pd.Serie
    x ={#'El Nino':nino1.co2-nino1.select_model,
        'La Nina':nina1.co2-nina1.select_model,
        'Neutral':neutral1.co2-neutral1.select_model,
        'CP Nino':modoki.co2-modoki.select_model,
 #       'Cold Modoki':nmodoki.co2-nmodoki.select_model,
        'EP Nino':EP.co2-EP.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':EP.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'Modoki CO2':modoki.co2,
        'El Nino NP':EP.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model,
        'Modoki NP':modoki.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name[0:3]+'\xb0'+mooring_name[3:4]) #Add degrees to the lab
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    if 'EP' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='blue'
    elif 'CP' in x[0]:
        c='darkred'
        ls='--'

            
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()

plt.ylim([xl2,yl2])
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('e)  Air-sea flux - new production ENSO',loc='left')
#print(final_mooring_enso_avgs)
print('AIR SEA FLUX MINUS NEW PROD')
print(final_mooring_enso.mean())
ENSOAVG=final_mooring_enso_avgs


# %% New Production


ax = plt.subplot(6,2,2)
cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    newprod=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),newprod,c=cols[i],linewidth=3,label=moorings[i])

plt.title('b) New production seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()
plt.ylim(xl1,yl1)


ax = plt.subplot(6,2,1)

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
            
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
   # nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    #pd.Serie
    x ={#'El Nino':nino1.select_model,
        'La Nina':nina1.select_model,
        'Neutral':neutral1.select_model,
        'Modoki':modoki.select_model,
    #    'Cold Modoki':nmodoki.select_model,
        'EP':EP.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={'El Nino CO2':EP.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'Modoki CO2':modoki.co2,
        'El Nino NP':EP.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model,
        #'Modoki NP':modoki.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name[0:3]+'\xb0'+mooring_name[3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'Modoki' in x[0]:
        c='darkred'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
            
            
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()

plt.ylim([xl0,yl0])
print('NEW PRODUCTION')
print(final_mooring_enso.mean())
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('a)  New production ENSO',loc='left')

# %%    PCO2 section
   
plt.subplot(628)


cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4']#[::-1]


pco2_month = xr.open_dataarray('processed/flux/pco2grams.nc')
pco2_month = (pco2_month*12*50)#/30

LandSch_co2flux_data=xr.open_mfdataset('processed/flux/landsch_mooring_co2_flux.nc').rename({'time':'Date'})
LandSch_co2flux_data=LandSch_co2flux_data.rename({'Date':'time'})
mooringz_lon=['165E', '170W', '155W', '140W', '125W', '110W']

for i, mooring_name in enumerate(mooringz_lon):
    pco2_m=LandSch_co2flux_data.dco2.sel(Mooring=mooring_name)#.plot()
    
    #pco2_m=pco2_month.sel(lat=0,lon=lns[i],method='nearest')
    
    
    pco2=pco2_m.groupby('time.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(pco2_m.groupby('time.month').mean()-pco2_m.groupby('time.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])
plt.ylim([0,110])
plt.title('h) \u0394pCO$_{2}$ seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('μatm')
plt.grid()

#ENSO BREAKDOWN
plt.subplot(627)
#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()


LandSch_co2flux_data=LandSch_co2flux_data.rename({'time':'Date'})

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
    pco2_m=LandSch_co2flux_data.dco2.sel(Mooring=moorings[i])
    pco2_m['Date']=pco2_m['Date'].astype('datetime64[M]')
    dat['pco2t']=pco2_m
#for i, mooring_name in enumerate(moorings):
    print(mooring_name)
  
    dat['Date'].astype('datetime64[M]')
    
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.select_model)]
    
   
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
   # nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    #pd.Serie
    x ={'El Nino':nino1.co2-nino1.select_model,
        'La Nina':nina1.co2-nina1.select_model,
        'Neutral':neutral1.co2-neutral1.select_model,
        'Modoki':modoki.co2-modoki.select_model,
    #    'Cold Modoki':nmodoki.co2-nmodoki.select_model,
        'EP':EP.co2-EP.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={#'El Nino CO2':nino1.pco2t,
        'La Nina CO2':nina1.pco2t,
        'Neutral CO2':neutral1.pco2t,
        'Modoki CO2':modoki.pco2t,
    #    'Cold Modoki CO2':nmodoki.pco2t,
        
        'EP CO2':EP.pco2t
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=moorings[i][0:3]+'\xb0'+moorings[i][3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=moorings[i])
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso_avgs.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'Modoki' in x[0]:
        c='darkred'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
            
            
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()
plt.ylabel('μatm')
plt.ylim([0,110])
plt.title('g) \u0394pCO$_{2}$ ENSO',loc='left')
print('pCO2')
print(final_mooring_enso_avgs.mean())


# %% How about the SST
plt.subplot(6,2,10)


cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=xx.sst_rey.groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])

plt.ylim([23,31])
plt.title('j) SST seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
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


plt.subplot(6,2,9)
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

    #d=dat.sel(Date=slice(pco2_month.Date.min().values,pco2_month.Date.max().values))
    
    #pco2_m=d.sst_rey#seasonaltrend_sstobs[i]
    #mei=d.mei
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.co2)]
       
   
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
   # nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    
    #pd.Serie
    x ={#'El Nino':nino1.sst_rey,
        'La Nina':nina1.sst_rey,
        'Neutral':neutral1.sst_rey,
        'Modoki':modoki.sst_rey,
#        'Cold Modoki':nmodoki.sst_rey,
        'EP':EP.sst_rey
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    
    ensoavgs=pd.Series(x,name=moorings[i][0:3]+'\xb0'+moorings[i][3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'Modoki' in x[0]:
        c='darkred'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
            
            
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()

plt.xlabel('Mooring')
plt.ylabel('Degree C')
plt.ylim([23,31])
plt.title('i) SST ENSO',loc='left')
print('SST')
print(final_mooring_enso.mean())


# %% Air Sea Flux

ax = plt.subplot(6,2,4)

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.co2flux4_land_gmyr/365).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[i])

   
plt.ylim([xl1,yl1])  
plt.title('d) Air-sea flux seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()


ax = plt.subplot(6,2,3)


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
    
        
   
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
  #  nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    #pd.Serie
    x ={#'El Nino':nino1.co2,
        'La Nina':nina1.co2,
        'Neutral':neutral1.co2,
        'Modoki':modoki.co2,
  #      'Cold Modoki':nmodoki.co2,
        'EP':EP.co2
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    

    avgs ={#'El Nino CO2':nino1.co2,
        'La Nina CO2':nina1.co2,
        'Neutral CO2':neutral1.co2,
        'El Nino NP':EP.select_model,
        'La Nina NP':nina1.select_model,
        'Neutral NP':neutral1.select_model
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    ensoavgs=pd.Series(x,name=mooring_name[0:3]+'\xb0'+mooring_name[3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    avgz=pd.Series(avgs,name=mooring_name)
    final_mooring_enso_avgs=final_mooring_enso_avgs.append(avgz)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    elif 'Nino' in x[0]:
        c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'Modoki' in x[0]:
        c='darkred'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
            
            
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()

plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('c)  Air-sea flux ENSO',loc='left')
plt.ylim([xl0,yl0])
print('AIR SEA FLUX')
print(final_mooring_enso.mean())

#%% How about the windspeed

plt.subplot(6,2,12)

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.windspeed).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=3,label=moorings[::-1][i])

plt.ylim([3,7])
plt.title('l) Windspeed seasonality',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
plt.legend(ncol=2)
plt.xlabel('Month')
plt.ylabel('Windspeed (m s$^{-1}$)')
plt.grid()

#ENSO BREAKDOWN

plt.subplot(6,2,11)

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

    d=dat.sel(Date=slice(pco2_month.time.min().values,pco2_month.time.max().values))
    
    #pco2_m=d.windspeed#seasonaltrend_sstobs[i]
    #mei=d.mei
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.laws2011a*info.cafe/1000 #Was cbpmmean
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.co2)]
        
    
    nino1=info[info.index.isin(ep_dates)]#info[info.mei>0.5]#.mean()
    modoki=info[info.index.isin(cp_dates)]#info[info.emi>0.5]#.mean()
    #nmodoki=info[info.emi<-0.5]
    nina1=info[info.index.isin(nina_dates)]#info[info.mei<-0.5]#.mean()
    neutral1=info[~info.index.isin(cp_dates)]#info[(info.mei<0.5)&(info.mei>-0.5)]#.mean()
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral1=neutral1[~neutral1.index.isin(nina_dates)]
    
    EP=nino1.mean() #Remove if it is in MODOKI
    modoki=modoki.mean()
    #nmodoki=nmodoki.mean()
    nina1=nina1.mean()
    neutral1=neutral1.mean()
    nino1=nino1.mean()
    
    #pd.Serie
    x ={#'El Nino':nino1.windspeed,
        'La Nina':nina1.windspeed,
 #       'Cold Modoki':nmodoki.windspeed,
        'Neutral':neutral1.windspeed,
        'CP El Nino':modoki.windspeed,
        'EP El Nino':EP.windspeed
       # 'All time':info.co2.mean()-info.select_model.mean()
        }
    
    ensoavgs=pd.Series(x,name=moorings[i][0:3]+'\xb0'+moorings[i][3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    print(x[0])
    print(x[0]=='CP El Nino')
    ls='-'
    if 'Neutral' in x[0]:
        c='black'
    #elif 'Nino' in x[0]:
    #    c='darkred'
    elif 'Nina' in x[0]:
        c='royalblue'
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
    elif 'CP' in x[0]:
        c='darkred'
        ls='--'
            
    
    plt.plot(x[1].index,x[1],c=c,linestyle=ls,label=x[0],linewidth=3,alpha=0.8)
plt.grid()
plt.legend(final_mooring_enso.columns)
plt.ylim([3,7])
plt.xlabel('Mooring')
plt.ylabel('Windspeeed (m s$^{-1}$)')
plt.title('k) Windspeed ENSO',loc='left')
print('WINDSPEED')
print(final_mooring_enso.mean())
plt.tight_layout()

plt.savefig('figs/Figure5a_ENSO_seasonality.png',dpi=200)
plt.savefig('figs/Figure5a_ENSO_seasonality.eps')
plt.show()

print(ENSOAVG.mean())




# %% Check lag between NP and ASF in the east


fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

xx=dat.sel(Mooring=110)

asf=((xx.co2flux4_land_gmyr/365)).groupby('Date.month').mean()
np=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()
year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
plt.xcorr(asf,np,maxlags=6)