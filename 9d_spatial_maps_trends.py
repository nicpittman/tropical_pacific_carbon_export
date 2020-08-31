#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:40:28 2020
@author: Nic Pittman

This code will reproduce Figure 4 in Pittman et al., 2021. 

Trends and pvalues are calculated on the fly and not saved anywhere, however could be done easily. 
regridded data is required for this process

This results in a slower script but works well. All of the processing occurs in the main function.
Easy to call modified version of this figure.

Produces mean, trend and pval (Stipples) for the following:
    
    figs/Figure4_Spatial_map_update_'+ratio.name+'.png
    
    air-sea flux
    new production 
    difference is calculated here
    SST
    TPCA chlorophyll (regridded) 
    carbon (as processed into grams)
    
Requires: 
        datasets/co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc
        processed/seamask.nc
        processed/flux/fratios.nc
    
        processed/flux/avg_npp_rg_cafe.nc'
        processed/flux/tpca.nc
        datasets/sst/sst.mnmean.nc
        processed/flux/pco2grams.nc
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
from mpl_toolkits.basemap import Basemap
from scipy.stats import linregress

    
def plot_basemap():
    m = Basemap(llcrnrlon=120.,llcrnrlat=-14.5,urcrnrlon=290,urcrnrlat=14.51,
                resolution='l',projection='merc',fix_aspect=False)
    m.drawcoastlines()
    m.fillcontinents()
    # draw parallels     # labels = [left,right,top,bottom]
    m.drawparallels(np.arange(-20,21,10),labels=[1,0,1,1],fontsize=12,latmax=20)
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],fontsize=12)
    return m



def plot_basemap_row(fig,axn,hovmol,units,title,units_tr,levs=None,levs_trend=None,trend_conversion=None,sb1=7,sb2=2,cmap='viridis',cmaptr='RdBu_r'):
    '''
    Create a plotting function to make it repeatable and nicer
    colormaps should either be viridis or RdBu_r
    axis (number) will be 1,3,5,7 (plots both avg and trend at once)
     
    Unfortunately this function does the processing of mean, trends and pvals on the fly.
    Could save these if needed, but not provided here. 
    '''
    fr=0.03
    fs=12
    ms=10
    startday=np.datetime64('2000-01-01')
    
    if title.endswith('pCO2t'):
        endday=np.datetime64('2016-12-01') 
        print(title)
    elif title.endswith('chlorophyll'):
        endday=np.datetime64('2017-12-01')
    else:
        endday=np.datetime64('2020-01-01') 
        
    ax1=fig.add_subplot(sb1,sb2,axn)
    m=plot_basemap()

    lo,la=np.meshgrid(hovmol.lon.values,hovmol.lat.values)
    lo1,la1=m(lo,la)
    
    if type(levs_trend)==type(None):
        f=m.contourf(lo1,la1,hovmol.mean(dim='time'),cmap=cmap) #11 colors
    else:
        if title=='TPCA Chlorophyll':
            f=m.contourf(lo1,la1,hovmol.mean(dim='time'),extend='max',cmap=cmap,levels=levs) #11 colors
        else:
            f=m.contourf(lo1,la1,hovmol.mean(dim='time'),cmap=cmap,levels=levs) #11 colors

    ax1.axhline(0,c='k',linestyle=':')

    moorings=[165,190,205,220,235,250]
    for x in moorings:
        x1,y1=m(x,0)
        ax1.plot(x1,y1,marker='x',c='k',markersize=ms)
    
    if title=='SST':
        
        lev=28.5#29.2 #rather than 28.5
        early_sst=hovmol.sel(time=slice('1997-01-01','2002-01-01')).mean(dim='time')#.where(co2.seamask==1)
        m.contour(lo1,la1,early_sst,levels=[lev],linestyles='dotted',colors='k')
       
        late_sst=hovmol.sel(time=slice('2015-01-01','2020-01-01')).mean(dim='time')#.where(co2.seamask==1)
        m.contour(lo1,la1,late_sst,levels=[lev],linestyles='solid',colors='k')


    cb=plt.colorbar(f,ax=ax1,fraction=fr)
    cb.set_label(units,fontsize=fs)
    cb.ax.tick_params(labelsize=fs-1)
    ax1.set_title(chr(ord('`')+axn)+') Average: '+title,fontsize=fs)
    ax1.tick_params(labelsize=fs)

    #Trends
    hovmol=hovmol.where(hovmol!=-0.9999,np.nan)
    hm=hovmol.interpolate_na(dim='time').sel(time=slice(startday,endday))
    months=hm.time
    
    dt_dates=pd.to_numeric(months.values.astype('datetime64[D]'))
    num_dates=dt_dates
    hm['time']=num_dates



    #This will calculate the per pixel trends and pvalues

    time=hm.time.values
    xx=np.concatenate(hm.T)
    #print(xx.shape)
    tr=[]
    pv=[]
    for i in range(xx.shape[0]):
        #print(xx[i,:])
        stat=linregress(time,xx[i,:])
        #print(stat)
        tr.append(stat.slope*365)
        pv.append(stat.pvalue)
        
    tr=np.array(tr).reshape(len(hm.lon),len(hm.lat)).T
    pv=np.array(pv).reshape(len(hm.lon),len(hm.lat)).T
    
    hh=hm.copy()
    hh=hh.drop('time')
    hh['trend']=(['lat','lon'],tr)
    hh['pval']=(['lat','lon'],pv)
    ##hh now contains a flat xarray with trend and pvalue, could save this for each variable. 
    
    if type(trend_conversion)!=type(None):
        hh['trend']=hh['trend']*trend_conversion
   
    
    #print(hh.sel(lat=0,lon=165,method='nearest').mean(dim='time'))    
    ax2=plt.subplot(sb1,sb2,axn+1)
    
    m1=plot_basemap()
    
    if type(levs_trend)==type(None):
        f=m1.contourf(lo1,la1,hh.trend,cmap=cmaptr,extend='both')
    else:
        f=m1.contourf(lo1,la1,hh.trend,cmap=cmaptr,extend='both',levels=levs_trend) #11 colors 
    cb=plt.colorbar(f,ax=ax2,extend='both',fraction=fr)
        
    m1.contourf(lo1,la1,hh.pval,colors='none',hatches=['.'],levels=[0,0.05])
    #m1.contour(lo1,la1,hh.pval,colors='k',levels=[0.05])
    
    cb.set_label(units_tr,fontsize=fs)
    cb.ax.tick_params(labelsize=fs)
    ax2.axhline(0,c='k',linestyle=':')
    
    ax2.set_title(chr(ord('`')+axn+1)+') Trends: '+title,fontsize=fs)
    for x in moorings:
        x1,y1=m(x,0)
        ax2.plot(x1,y1,marker='x',c='k',markersize=ms)
        
        ax2.tick_params(labelsize=fs)

    return hh


#Functions above make plotting easy.
# # Code begins
    
# %%Load data in
         
#landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc'



seamask=xr.open_dataset('processed/seamask.nc') #Because 2020 version doesn't have it.
seamask= seamask.assign_coords(lon=(seamask.lon % 360)).roll(lon=(seamask.dims['lon']),roll_coords=False).sortby('lon')	

#It would be preferable to use the 2020 version,
# landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc'
#However it doesn't include seamask so we are going to need both.... (Unless I save the seamask)
landschutzer=xr.open_dataset(landsch_fp)
landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
land_pac=land_pac.fgco2_smoothed

f_ratios=xr.open_mfdataset('processed/flux/fratios.nc')
ratio=f_ratios.laws2011a#laws2000#laws2000,laws2011a,laws2011b,henson2011

npp1=xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc')
avg_npp=(npp1.avg_npp/1000)*ratio

land=moles_to_carbon(land_pac)/365  #LANDSCHUTZ

land['time']=land.time.astype('datetime64[M]')

diff=land-avg_npp
diff1=diff.where((diff<0.1)|(diff<-0.1),np.nan)


# Need to combine the chlorophyll products, takes a bit of memory.
chl=xr.open_dataset('processed/flux/tpca.nc').tpca#'sw_month.nc')
#mod=xr.open_dataset('datasets/tpca/mod_month.nc')
#sw['time']=sw.time.astype('datetime64[M]')
#mod['time']=mod.time.astype('datetime64[M]')
#tpca=sw
#tpca=tpca.merge(mod)
#chl = tpca.to_array(dim='tpca').mean('tpca')

#SST
sst = xr.open_dataset('datasets/sst/sst.mnmean.nc')
sst= sst.assign_coords(lon=(sst.lon % 360)).roll(lon=(sst.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
sst=sst.sel(lon=slice(120,290),lat=slice(20,-20)).sst
sst=sst.where(seamask.seamask==1)

pCO2 = xr.open_dataarray('processed/flux/pco2grams.nc')
integratedpCO2 = (pCO2*12*50)
#monthlyPCO2=integratedpCO2.diff('time',1)/30

# %% Prepare Figure 

fig=plt.figure(figsize=(19*2/2.54,23*2/2.54))#(figsize=(30,15))
sb1=7
sb2=2

plot_basemap_row(fig,axn=1,
                 hovmol=avg_npp.sel(lat=slice(-15,15)),
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='New Production',
                 units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',
                 levs=np.arange(0,0.26,0.025),
                 levs_trend=np.arange(-2,2.1,0.25),
                 trend_conversion=1000,
                 cmap='viridis')

plot_basemap_row(fig,axn=3,
                 hovmol=land,
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='Air-Sea CO2 flux',
                 units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',
                 levs=np.arange(-0.12,0.13,0.02),
                 levs_trend=np.arange(-2,2.1,0.5),
                 trend_conversion=1000,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=5,
                 hovmol=diff1,
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='CO2 flux - New Production',
                 units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',                 
                 levs=np.arange(-0.12,0.13,0.02),
                 
                 levs_trend=np.arange(-2,2.1,0.5),
                 trend_conversion=1000,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=7,
                 hovmol=sst,
                 units='Degrees C',
                 title='SST',
                 units_tr='Degrees C year$^{-1}$',
                 levs=np.arange(20,32,1),
                 
                 levs_trend=np.arange(-0.06,0.07,0.01),
                 #trend_conversion=1000,
                 cmap='viridis')


plot_basemap_row(fig,axn=9,
                 hovmol=chl,
                 units='mg chl m$^{-3}$ day$^{-1}$',
                 title='TPCA chlorophyll',
                 units_tr='ug chl m$^{-3}$ day$^{-1}$ year$^{-1}$',
                 levs=np.arange(0,0.65,0.05),
                 
                 levs_trend=np.arange(-4,4.1,1),
                 trend_conversion=1000,
                 cmap='viridis')

plot_basemap_row(fig,axn=11,
                  hovmol=integratedpCO2,#monthlyPCO2*1000,
                  units='gC m$^{-2}$',
                  title='pCO2t',
                  units_tr='mgC m$^{-2}$ year$^{-1}$',
                  levs=np.arange(5.5,9.5,0.5),
                  levs_trend=np.arange(0,100,10),
                  trend_conversion=1000,
                  cmap='viridis',
                  cmaptr='Reds')



plt.tight_layout()
plt.savefig('figs/Figure4_Spatial_map_update_'+ratio.name+'.png',dpi=100)
plt.show()