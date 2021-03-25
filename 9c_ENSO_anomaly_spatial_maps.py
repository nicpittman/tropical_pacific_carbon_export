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
    m = Basemap(llcrnrlon=120.,llcrnrlat=-15,urcrnrlon=290,urcrnrlat=15.01,
                resolution='l',projection='merc',fix_aspect=False)
    m.drawcoastlines()
    m.fillcontinents()
    # draw parallels     # labels = [left,right,top,bottom]
    m.drawparallels(np.arange(-20,21,10),labels=[1,0,1,1],fontsize=12,latmax=20)
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],fontsize=12)
    return m



def plot_basemap_row(fig,axn,hovmol,mean,units,title,levs=None,levs_trend=None,trend_conversion=None,sb1=7,sb2=3,cmap='viridis',cmaptr='RdBu_r',wu=None,wv=None):
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
    
    if type(levs)==type(None):
        f=m.contourf(lo1,la1,hovmol.mean(dim='time')-mean.mean(dim='time'),cmap=cmap,extend='both') #11 colors
        #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
        for c in f.collections:
            c.set_edgecolor("face")
    else:
        if title=='TPCA Chlorophyll':
            f=m.contourf(lo1,la1,hovmol.mean(dim='time')-mean.mean(dim='time'),extend='both',cmap=cmap,levels=levs) #11 colors
                #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
            for c in f.collections:
                c.set_edgecolor("face")
        else:
            
            f=m.contourf(lo1,la1,hovmol.mean(dim='time')-mean.mean(dim='time'),cmap=cmap,levels=levs,extend='both') #11 colors
            #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
            for c in f.collections:
                c.set_edgecolor("face")
    ax1.axhline(0,c='k',linestyle=':')

    moorings=[165,190,205,220,235,250]
    for x in moorings:
        x1,y1=m(x,0)
        ax1.plot(x1,y1,marker='x',c='k',markersize=ms)
    
    if title=='SST':
        
        
        meansst=hovmol.mean(dim='time')#.where(co2.seamask==1)
        m.contour(lo1,la1,meansst,levels=[28.5],linestyles='solid',colors='k')
        m.contour(lo1,la1,meansst,levels=[25],linestyles='dashed',colors='k')
        
    
    #wu['lon'],wu['lat']=m(lo,la,wu.lon.values,wu.lat.values)
    #No windspeed vectors now
    #if title=='Wind speed':
    #      skip=(slice(None,None,4),slice(None,None,4)) #2 for NCEP 2
    #      m.quiver(lo1[skip],la1[skip],wu.mean(dim='time')[skip]/2,wv.mean(dim='time')[skip]/2,scale=90,headwidth=4.5)#,minshaft=2)


    cb=plt.colorbar(f,ax=ax1,fraction=fr,extend='both')
    cb.set_label(units,fontsize=fs)
    cb.ax.tick_params(labelsize=fs-1)

    if axn==1:
        name='EP Events'
    elif axn==2:
        name='CP Events'
    elif axn==3:
        name='La Nina Events'
    if axn<=3:
        ax1.set_title(name+'\n'+chr(ord('`')+axn)+') Anomaly: '+title,fontsize=fs)
    else:
        ax1.set_title(chr(ord('`')+axn)+') Anomaly: '+title,fontsize=fs)
            
    ax1.tick_params(labelsize=fs)

    #Trends
    
    #hovmol=hovmol.where(hovmol!=-0.9999,np.nan)
    #hm=hovmol.interpolate_na(dim='time').sel(time=slice(startday,endday))
    #months=hm.time
    
    #dt_dates=pd.to_numeric(months.values.astype('datetime64[D]'))
    #num_dates=dt_dates
    #hm['time']=num_dates

    
 


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
land_pac['time']=land_pac.time.astype('datetime64[M]')
land_pac_all=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))

land_pac=land_pac.fgco2_smoothed

atmco2=land_pac_all.atm_co2
dco2=land_pac_all.dco2
pco2=land_pac_all.spco2_smoothed
kw=land_pac_all.kw

f_ratios=xr.open_mfdataset('processed/flux/fratios.nc')
ratio=f_ratios.laws2011a#laws2000#laws2000,laws2011a,laws2011b,henson2011

npp1=xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc')
avg_npp=(npp1.avg_npp/1000)*ratio

land=moles_to_carbon(land_pac)/365  #LANDSCHUTZ


diff=land-avg_npp
diff1=diff.where((diff<0.1)|(diff<-0.1),np.nan)


# Need to combine the chlorophyll products, takes a bit of memory.
chl=xr.open_dataset('processed/flux/tpca.nc').tpca#'sw_month.nc')

#mod=xr.open_dataset('datasets/tpca/mod_month.nc')
chl['time']=chl.time.astype('datetime64[M]')
#mod['time']=mod.time.astype('datetime64[M]')
#tpca=sw
#tpca=tpca.merge(mod)
#chl = tpca.to_array(dim='tpca').mean('tpca')

#SST
sst = xr.open_dataset('datasets/sst/sst.mnmean.nc')
sst= sst.assign_coords(lon=(sst.lon % 360)).roll(lon=(sst.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
sst=sst.sel(lon=slice(120,290),lat=slice(20,-20)).sst
sst=sst.where(seamask.seamask==1)

pCO2 = xr.open_dataarray('processed/flux/pco2grams.nc') #_norm
integratedpCO2 = (pCO2*12*50)

#wu=xr.open_dataset('datasets/uwnd.mon.mean.nc').sel(level=1000,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
#wv=xr.open_dataset('datasets/vwnd.mon.mean.nc').sel(level=1000,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd
wu=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
wv=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd
dco2['time']=dco2.time.astype('datetime64[M]')

ws=np.sqrt((wu**2)+(wv**2))



precip= xr.open_dataset('datasets/precip.mon.mean.enhanced.nc').sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).precip

# # THIS NEEDS TO BE RUN ONCE BUT CAN be memory intensive

# w_ccmp_a=xr.open_mfdataset('datasets/ws_ccmp/*.nc') #Downloaded manually
# w_ccmp_a['time']=w_ccmp_a.time.astype('datetime64[M]')
# w_ccmp_a=w_ccmp_a.sel(latitude=slice(-20,20))

# w_ccmp_b=xr.open_mfdataset('datasets/CCMP_winds.nc') #Bulk ErDap download
# dt=w_ccmp_b.indexes['time'].to_datetimeindex()
# w_ccmp_b['time']=dt

# w_ccmp=xr.merge([w_ccmp_b,w_ccmp_a])


# w_ccmp=w_ccmp.sel(longitude=slice(120,290),latitude=slice(-20,20))
# ws_ccmp=np.sqrt((w_ccmp.uwnd**2)+(w_ccmp.vwnd**2))
# ws_ccmp=ws_ccmp.rename({'latitude':'lat','longitude':'lon'})
# try:
#     ws_ccmp.to_netcdf('datasets/CCMP_windspeed.nc')
#     print('saved')
# except:
#     pass

ws_ccmp=xr.open_dataarray('datasets/CCMP_windspeed.nc')
ws_ccmp=xr.open_dataarray('processed/CCMP_ws_1deg.nc')

# %% Prepare Figure 


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
ep_dates=ep.index[4:]
cp_dates=cp.index
#all_dates=chl.time
all_dates=info.index#[36:] #2000 - 2020
old_all_dates=all_dates
neutral=all_dates.drop(cp_dates).drop(nina_dates).drop(ep_dates)
all_dates=neutral
fig=plt.figure(figsize=(19*2/2.54,23*2/2.54))#(figsize=(30,15))
sb1=7
sb2=3


sst_range=np.arange(-2.5,2.75,0.25)
ws_range=np.arange(-2,2.2,0.2)
chl_range=np.arange(-0.1,0.11,0.01)
npp_range=np.arange(-0.04,0.0425,0.0025)
precip_range=np.arange(-6,7,1)
dco2_range=np.arange(-50,55,5)
co2_range=np.arange(-0.04,0.045,0.005)
#%% EP

plot_basemap_row(fig,axn=1,
                 hovmol=sst.sel(time=ep_dates,method='nearest'),
                 mean=sst.sel(time=all_dates,method='nearest'),
                 units='Degrees C',
                 title='SST',
                 levs=sst_range,
             
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=4,
                 hovmol=ws_ccmp.sel(time=ep_dates,method='nearest'),
                 mean=ws_ccmp.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed',              
                 levs=ws_range,
                 
             
                 cmap='RdBu_r')

plot_basemap_row(fig,axn=7,
                 hovmol=chl.sel(time=ep_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$ day$^{-1}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 
                 trend_conversion=1000,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=10,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=ep_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='New production',
  
                 levs=npp_range,
        
                 trend_conversion=1000,
                 cmap='RdBu_r')

plot_basemap_row(fig,axn=13,
                 hovmol=precip.sel(time=ep_dates,method='nearest'),
                 mean=precip.sel(time=all_dates,method='nearest'),
                 units='mm day$^{-1}$',
                 title='Precipitation',
       
                 levs=precip_range,

                 cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=16,
                  hovmol=dco2.sel(time=ep_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
     
                  levs=dco2_range,#(200,1200,10),#(5.5,9.5,0.5),
    
                  cmap='RdBu_r')

plot_basemap_row(fig,axn=19,
                  hovmol=land.sel(time=ep_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='gC m$_{-2}$ day$^{-1}$',
                  title='Air-sea CO$_{2}$ flux',
                  
                  levs=co2_range,

                  cmap='RdBu_r')


# %% CP

plot_basemap_row(fig,axn=2,
                 hovmol=sst.sel(time=cp_dates,method='nearest'),
                 mean=sst.sel(time=all_dates,method='nearest'),
                 units='Degrees C',
                 title='SST',
                 levs=sst_range,#np.arange(20,32,1),
                 
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=5,
                 hovmol=ws_ccmp.sel(time=cp_dates,method='nearest'),
                 mean=ws_ccmp.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed',             
                 levs=ws_range,#p.arange(0,11,1),
               
                 cmap='RdBu_r',
                 wu=wu,wv=wv)

plot_basemap_row(fig,axn=8,
                 hovmol=chl.sel(time=cp_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$ day$^{-1}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=11,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=cp_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='New production',
                 levs=npp_range,

                 cmap='RdBu_r')

plot_basemap_row(fig,axn=14,
                 hovmol=precip.sel(time=cp_dates,method='nearest'),
                 mean=precip.sel(time=all_dates,method='nearest'),
                 units='mm day$^{-1}$',
                 title='Precipitation',
                 levs=precip_range,

                 cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=17,
                  hovmol=dco2.sel(time=cp_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
                  levs=dco2_range,

                  cmap='RdBu_r')

plot_basemap_row(fig,axn=20,
                  hovmol=land.sel(time=cp_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='gC m$_{-2}$ day$^{-1}$',
                  title='Air-sea CO$_{2}$ flux',
                  levs=co2_range,
                  
                  cmap='RdBu_r')


#%% NINA


plot_basemap_row(fig,axn=3,
                 hovmol=sst.sel(time=nina_dates,method='nearest'),
                 mean=sst.sel(time=all_dates,method='nearest'),
                 units='Degrees C',
                 title='SST',
                 levs=sst_range,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=6,
                 hovmol=ws_ccmp.sel(time=nina_dates,method='nearest'),
                 mean=ws_ccmp.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed',
                            
                 levs=ws_range,
                 cmap='RdBu_r')

plot_basemap_row(fig,axn=9,
                 hovmol=chl.sel(time=nina_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$ day$^{-1}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=12,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=nina_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='New production',
                 levs=npp_range,

                 cmap='RdBu_r')

plot_basemap_row(fig,axn=15,
                 hovmol=precip.sel(time=nina_dates,method='nearest'),
                 mean=precip.sel(time=all_dates,method='nearest'),
                 units='mm day$^{-1}$',
                 title='Precipitation',
                 
                 levs=precip_range,

                 cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=18,
                  hovmol=dco2.sel(time=nina_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
                  levs=dco2_range,
 
                  cmap='RdBu_r')

plot_basemap_row(fig,axn=21,
                  hovmol=land.sel(time=nina_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='gC m$_{-2}$ day$^{-1}$',
                  title='Air-sea CO$_{2}$ flux',
                  levs=co2_range,

                  cmap='RdBu_r')






plt.tight_layout()
plt.savefig('figs/Figure4_ENSO_Anomaly_spatial_map_.png',dpi=200)
# plt.savefig('figs/vector/Figure3_Spatial_map_'+ratio.name+etype+'.eps')
# plt.savefig('figs/vector/Figure3_Spatial_map_'+ratio.name+etype+'.pdf')

# try:
#     plt.savefig('figs/Figure3_Spatial_map_'+ratio.name+'.jpeg',dpi=300)
# except:
#     pass
# plt.show()