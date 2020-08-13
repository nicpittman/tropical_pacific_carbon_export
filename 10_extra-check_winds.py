#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 17:07:42 2020

@author: npittman
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import linregress
def trends(x,y,c):  
    from scipy.stats import linregress
    mean=np.nanmean(y)
    std=np.nanstd(y)*1

    mask=~np.isnan(x)
    x=x[mask]
    y=y[mask]
    mask1=~np.isnan(y)
    x=x[mask1]
    y=y[mask1]

    #x_n=np.arange(0,len(x))
    # x1=np.arange(np.datetime64(x[0],'M'),np.datetime64(x[-1],'M')+np.timedelta64(1,'M'))
    #x1=trd.index.values.astype('datetime64[D]')
    #x1=x.values.astype('datetime64[D]')
    #pd.to_numeric(x1)
    slope, intercept, r_value, p_value,std_err = linregress(x,y)
    mn=min(x)
    mx=max(x)
    x1=np.linspace(mn,mx,len(x))
    y1=slope*x1+intercept
    
    plt.plot(x1,y1,c,linestyle='--',linewidth=2.5)  
    #ax.text(x1[-1]-(x1[-1]*0.1),y1[-1]-(y1[-1]*0.1),'R2='+str(np.round(r_value**2,3)))
    print( slope, intercept, r_value**2,p_value,std_err)
    return slope, intercept, r_value,p_value,std_err
    
plt.figure(figsize=(20,10))
plt.subplot(321)
fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
dat.windspeed.groupby('Date.month').mean().plot()
plt.title('Windspeed seasonal hovmoller')
#plt.show()

plt.subplot(343)
plt.scatter(dat.delta_pCO2,(dat.co2flux4_land_gmyr/365)*1000,label='delta pCO2 (uatm)',alpha=0.45)
plt.xlabel('pCO2 (uatm')
plt.ylabel('CO2 flux (mgC/day)')
plt.title('pCO2 vs CO2 flux')
trends(dat.delta_pCO2.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='b')

plt.subplot(344)
plt.scatter(dat.delta_pCO2,dat.laws2011a*dat.cafe,label='delta pCO2 (uatm)',alpha=0.45)
plt.ylabel('New Production (mgC/day)')
plt.xlabel('pCO2 (uatm')
plt.title('pCO2 vs New Production')
trends((dat.laws2011a*dat.cafe).values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='b')


plt.subplot(323)
plt.scatter(dat.windspeed,dat.laws2011a*dat.cafe,label='NP (mgC/day)',alpha=0.9)
plt.scatter(dat.windspeed,(dat.co2flux4_land_gmyr/365)*1000,label='CO2 flux (mgC/day)',alpha=0.5)
plt.scatter(dat.windspeed,dat.delta_pCO2,label='delta pCO2 (uatm)',alpha=0.45)
#plt.legend()
trends(dat.windspeed.values.flatten(),(dat.laws2011a*dat.cafe).values.flatten(),c='b')
trends(dat.windspeed.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='orange')
trends(dat.windspeed.values.flatten(),((dat.delta_pCO2)).values.flatten(),c='darkgreen')
plt.xlabel('Windspeed m/s')
plt.ylabel('NP and CO2 flux (mgC/day)')
plt.title('Windspeed')
#plt.show()

plt.subplot(324)
plt.scatter(dat.winddirection,dat.laws2011a*dat.cafe,label='NP (mgC/day)',alpha=0.9)
plt.scatter(dat.winddirection,(dat.co2flux4_land_gmyr/365)*1000,label='CO2 flux (mgC/day)',alpha=0.5)
plt.scatter(dat.winddirection,dat.delta_pCO2,label='delta pCO2 (uatm)',alpha=0.45)
#plt.legend()
trends(dat.winddirection.values.flatten(),(dat.laws2011a*dat.cafe).values.flatten(),c='b')
trends(dat.winddirection.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='orange')
trends(dat.winddirection.values.flatten(),((dat.delta_pCO2)).values.flatten(),c='darkgreen')
plt.xlabel('wind direction Degrees')
plt.ylabel('NP and CO2 flux (mgC/day)')
plt.title('Wind Direction)')

plt.subplot(325)
plt.scatter(dat.wu,dat.laws2011a*dat.cafe,label='NP (mgC/day)',alpha=0.9)
plt.scatter(dat.wu,(dat.co2flux4_land_gmyr/365)*1000,label='CO2 flux (mgC/day)',alpha=0.5)
plt.scatter(dat.wu,dat.delta_pCO2,label='delta pCO2 (uatm)',alpha=0.45)
#plt.legend()
trends(dat.wu.values.flatten(),(dat.laws2011a*dat.cafe).values.flatten(),c='b')
trends(dat.wu.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='orange')
trends(dat.wu.values.flatten(),((dat.delta_pCO2)).values.flatten(),c='darkgreen')
plt.xlabel('wu speed (ms)')
plt.title('Zonal WU speed')
plt.ylabel('NP and CO2 flux (mgC/day)')
#plt.show()


plt.subplot(326)
plt.scatter(dat.wv,dat.laws2011a*dat.cafe,label='NP (mgC/day)',alpha=0.9)
plt.scatter(dat.wv,(dat.co2flux4_land_gmyr/365)*1000,label='CO2 flux (mgC/day)',alpha=0.5)
plt.scatter(dat.wv,dat.delta_pCO2,label='delta pCO2 (uatm)',alpha=0.45)
plt.legend()
trends(dat.wv.values.flatten(),(dat.laws2011a*dat.cafe).values.flatten(),c='b')
trends(dat.wv.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='orange')
trends(dat.wv.values.flatten(),((dat.delta_pCO2)).values.flatten(),c='darkgreen')
plt.xlabel('wv speed (ms)')
plt.title('Meridional WV speed')
plt.ylabel('NP and CO2 flux (mgC/day)')

plt.tight_layout()
plt.show()

#Seasonality

m_ws=dat.windspeed.groupby('Date.month').mean()
m_np=(dat.laws2011a*dat.cafe).groupby('Date.month').mean()
m_co2=((dat.co2flux4_land_gmyr/365)*1000).groupby('Date.month').mean()
m_pco2=dat.delta_pCO2.groupby('Date.month').mean()

for i,e in enumerate(m_pco2):
    f=m_ws.sel(Mooring=e.Mooring)
    plt.scatter(e,f,label=e.Mooring.values)
    print(i,e)
plt.legend()
#plt.show()


# x=(dat.laws2011a*dat.cafe).values.flatten()
# y=dat.windspeed.values.flatten()
# mask=~np.isnan(x)
# x=x[mask]
# y=y[mask]
# mask1=~np.isnan(y)
# x=x[mask1]
# y=y[mask1]
# slope, intercept, r_value, p_value,std_err = linregress(x,y)
# mn=min(x)
# mx=max(x)
# x1=np.linspace(mn,mx,len(x))
# y1=slope*x1+intercept
# #plt.show()
# plt.plot(y1,x1)