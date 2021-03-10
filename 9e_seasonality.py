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


#%% Seasonal Comparison  
#Chlorophyll

#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

plt.figure(figsize=(10,8))

month_fmt = DateFormatter('%b')
def m_fmt(x, pos=None):
    return month_fmt(x)#[0]

lw=1.5
ax = plt.subplot(4,2,3)

#cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

lns=[165,190,205,220,235,250]
moors=[110, 125, 140, 155, 170, 195]
moorings=['110$^\circ$W','125$^\circ$W','140$^\circ$W','155$^\circ$W','170$^\circ$W','165$^\circ$E']#[::-1]
#ZGot a bit confused here. Everything goes East (110W) to West (165E). Some mistakes should be fixed now.

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
precip= xr.open_dataset('datasets/precip.mon.mean.enhanced.nc').sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).precip


for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    #asfdiff=((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
    asfdiff=xx.chl.groupby('Date.month').mean()  
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((asfdiff.month-1).astype('datetime64[M]'),asfdiff,c=cols[i],linewidth=lw,label=moorings[i])
  
plt.ylim([0.05,0.35])
plt.title('c) TPCA Chlorophyll' ,loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('mg Chl m$^{-3}$ day$^{-1}$')
plt.grid()


#print(final_mooring_enso_avgs)
#print('AIR SEA FLUX MINUS NEW PROD')
#print(final_mooring_enso.mean())
#ENSOAVG=final_mooring_enso_avgs


# %% New Production


ax = plt.subplot(4,2,4)
#cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
#cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    newprod=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),newprod,c=cols[i],linewidth=lw,label=moorings[i])

plt.title('d) New production',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()
plt.ylim(xl1,yl1)

# %%    PCO2 section
   
plt.subplot(4,2,6)


#cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
#cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4']#[::-1]


pco2_month = xr.open_dataarray('processed/flux/pco2grams.nc')
pco2_month = (pco2_month*12*50)#/30

LandSch_co2flux_data=xr.open_mfdataset('processed/flux/landsch_mooring_co2_flux.nc').rename({'time':'Date'})
LandSch_co2flux_data=LandSch_co2flux_data.rename({'Date':'time'})
mooringz_lon=['165E', '170W', '155W', '140W', '125W', '110W'][::-1]

for i, mooring_name in enumerate(mooringz_lon):
    pco2_m=LandSch_co2flux_data.dco2.sel(Mooring=mooring_name)#.plot()
    
    #pco2_m=pco2_month.sel(lat=0,lon=lns[i],method='nearest')
    
    
    pco2=pco2_m.groupby('time.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(pco2_m.groupby('time.month').mean()-pco2_m.groupby('time.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=lw,label=moorings[i])
plt.ylim([0,110])
plt.title('f) \u0394pCO$_{2}$',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('μatm')
plt.grid()



# %% How about the SST
plt.subplot(4,2,1)


#cols=['#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4']
#cols=['#d73027','#fc8d59','#fee090','#d1e5f0','#91bfdb','#4575b4'][::-1]

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
    
    pco2=xx.sst_rey.groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=lw,label=moorings[i])

plt.ylim([23,31])
plt.title('a) SST',loc='left')#'Seasonal cycle in '+t,loc='left')
#plt.axhline(0,c='gray')
#plt.text(np.datetime64('1970-07-15'),0.025,'Air-Sea flux Dominated',fontsize=11)
#plt.text(np.datetime64('1970-07-15'),-0.02,'Biology Dominated',fontsize=11)
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.xlabel('Month')
plt.ylabel('Degree C')
plt.grid()



#%% How about the windspeed
lns=[165,190,205,220,235,250][::-1]
moors=[110, 125, 140, 155, 170, 195]


plt.subplot(4,2,2)
ws1d=xr.open_dataset('processed/ws_1deg.nc').windspeed
fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
for i, mooring in enumerate(moors):
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.windspeed).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())

    

    sat_ws=ws1d.sel(lat=0,lon=lns[i],method='nearest')
    sat_ws_group=sat_ws.groupby('time.month').mean()
    sat_year=(sat_ws.groupby('time.month').mean()-sat_ws.groupby('time.month').mean().mean())
 
    #plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linestyle='--',linewidth=lw,label=moorings[i])
    plt.plot((sat_year.month-1).astype('datetime64[M]'),sat_ws_group,c=cols[i],linewidth=lw,label=moorings[i])




  
#plt.ylim([3,7])
plt.title('b) Wind speed',loc='left')#'Seasonal cycle in '+t,loc='left')

plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))

plt.xlabel('Month')
plt.ylabel('m s$^{-1}$')
plt.grid()


# %% Precip
ax = plt.subplot(4,2,5)

lns=[165,190,205,220,235,250][::-1]
moors=[110, 125, 140, 155, 170, 195]



fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(moors):
    #print(moors, i,mooring)
    xx=dat.sel(Mooring=mooring)
    
    pco2=(xx.precip).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    sat_precip=precip.sel(lat=0,lon=lns[i],method='nearest')
    sat_precip_group=sat_precip.groupby('time.month').mean()
    sat_year=(sat_precip.groupby('time.month').mean()-sat_precip.groupby('time.month').mean().mean())
    #print('Satellite' + str(sat_precip_group.lon.values))
    #plt.plot((year.month-1).astype('datetime64[M]'),pco2*24,c=cols[i],linewidth=lw,label=moorings[i],linestyle='--')
    plt.plot((sat_year.month-1).astype('datetime64[M]'),sat_precip_group,c=cols[i],linewidth=lw,label=moorings[i])
    #print('Mooring'+ str(xx.Mooring.values))
   
#plt.ylim([xl1,yl1])  
plt.title('e) Precipitation',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('mm day$^{-1}$')
plt.grid()
#plt.legend(ncol=3,loc='upper left')


# %% Air Sea Flux

ax = plt.subplot(4,2,7)

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

for i, mooring in enumerate(dat.Mooring.values):
    xx=dat.sel(Mooring=mooring)
   
    pco2=(xx.co2flux4_land_gmyr/365).groupby('Date.month').mean()#((xx.co2flux4_land_gmyr/365)-(xx.laws2011a*xx.cafe/1000)).groupby('Date.month').mean()
      
    year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
    plt.plot((year.month-1).astype('datetime64[M]'),pco2,c=cols[i],linewidth=lw,label=moorings[i])
    print(i,mooring,xx.Mooring.values,moorings[i])
   
plt.ylim([xl1,yl1])  
plt.title('g) Air-sea CO$_{2}$ flux',loc='left')#'Seasonal cycle in '+t,loc='left')
plt.gca().xaxis.set_major_locator(MonthLocator())
plt.gca().xaxis.set_major_formatter(FuncFormatter(m_fmt))
#plt.legend()
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.grid()


# %%
plt.tight_layout()
plt.legend(ncol=3,fontsize=12,bbox_to_anchor=(1.19,0.5), loc="center left")#, borderaxespad=0)

plt.savefig('figs/Figure5.png',dpi=200)

try:
    plt.savefig('figs/Figure5.jpeg',dpi=300)
except:
    pass
plt.savefig('figs/vector/Figure5.eps')
plt.savefig('figs/vector/Figure5.pdf')
plt.show()


# %% Check lag between NP and ASF in the east


fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

xx=dat.sel(Mooring=110)

asf=((xx.co2flux4_land_gmyr/365)).groupby('Date.month').mean()
npr=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()
year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
plt.xcorr(asf,npr,maxlags=6)