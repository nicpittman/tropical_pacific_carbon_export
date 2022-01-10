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
from scipy.stats import ttest_ind, ttest_rel
#from windspharm.xarray import VectorWind
import matplotlib

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
   def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
       self.oom = order
       self.fformat = fformat
       matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
   def _set_order_of_magnitude(self):
       self.orderOfMagnitude = self.oom
   def _set_format(self, vmin=None, vmax=None):
       self.format = self.fformat
       if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format




def plot_basemap():
    m = Basemap(llcrnrlon=120.,llcrnrlat=-15,urcrnrlon=290,urcrnrlat=15.01,
                resolution='l',projection='merc',fix_aspect=False)
    m.drawcoastlines()
    m.fillcontinents()
    # draw parallels     # labels = [left,right,top,bottom]
    m.drawparallels(np.arange(-20,21,10),labels=[1,0,1,1],fontsize=12,latmax=20)
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],fontsize=12)
    return m



def plot_basemap_row(fig,
                     axn,
                     hovmol,
                     mean,
                     units,
                     title,
                     levs=None,

                     trend_conversion=None,
                     sb1=7,
                     sb2=3,
                     cmap='viridis',
                     cmaptr='RdBu_r',
                     wu=None,wv=None,
                     wu_all=None,wv_all=None):
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
        
    hovmol=hovmol.sel(time=slice(startday,endday))
    if wu is not None:
        wu=wu.sel(time=slice(startday,endday))
        wv=wv.sel(time=slice(startday,endday))
        wu_all=wu_all.sel(time=slice(startday,endday))
        wv_all=wv_all.sel(time=slice(startday,endday))
            
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
        
    
    
    #LETS do a t-test to see if they where there is significant differnces
    
    anom=np.concatenate(hovmol.T)
    mean=np.concatenate(mean.T)
     
    pv=[]
    for i in range(anom.shape[0]):
        #print(xx[i,:])
        stat=ttest_ind(anom[i],mean[i],nan_policy='omit')#linregress(time,xx[i,:])
        #print(stat)
        #tr.append(stat.slope*365)
        pv.append(stat.pvalue)
    
    pv=np.array(pv).reshape(len(hovmol.lon),len(hovmol.lat)).T
    
    hh=hovmol.copy()
    hh=hh.drop('time')
    hh['pval']=(['lat','lon'],pv)

    ####No windspeed vectors now
    #pass

  

    # if title=='Wind divergence':  
    #     windFmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    #     windFmt.set_powerlimits((0, 0))
    #     cnt=m.contourf(lo1,la1,hh.pval,colors='none',hatches=['.'],levels=[0,0.05])
    #     #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
    #     for c in cnt.collections:
    #         c.set_edgecolor("face")
    # elif title=='Wind divergence and direction': 

    if title=='Wind speed and direction':

        lo2,la2=np.meshgrid(wu.lon.values,wu.lat.values)
        lo2a,la2a=m(lo2,la2)
       
        skip=(slice(None,None,5),slice(None,None,5)) #2 for NCEP 2
        qu=m.quiver(lo2a[skip],
                    la2a[skip],
                    (wu.mean(dim='time')-wu_all.mean(dim='time'))[skip],
                    (wv.mean(dim='time')-wv_all.mean(dim='time'))[skip],
                    scale=17,headaxislength=4,headlength=5,headwidth=5)
        #x,y=m(-10,150)
        #ax1.quiverkey(qu,label='Wind direction m/s',labelpos='S',U=1,X=150,Y=-10)
        if axn!=6:
            plt.quiverkey(qu,1.3,1.028,U=1,label='Wind speed 5m s$^{-1}$')
        
    else:
        #windFmt=None
        cnt=m.contourf(lo1,la1,hh.pval,colors='none',hatches=['.'],levels=[0,0.05])
        #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
        for c in cnt.collections:
            c.set_edgecolor("face")


        
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
avg_npp=(npp1.avg_npp/12)*ratio

land=(land_pac*1000)/365  #LANDSCHUTZ


diff=land-avg_npp
diff1=diff.where((diff<0.1)|(diff<-0.1),np.nan)


# Need to combine the chlorophyll products, takes a bit of memory.
chl=xr.open_dataset('processed/flux/tpca.nc').tpca#'sw_month.nc')

#mod=xr.open_dataset('datasets/tpca/mod_month.nc')
chl['time']=chl.time.astype('datetime64[M]')
chl=chl.interpolate_na(dim='time')
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

#NCEP2 winds
#wu=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10).uwnd
#wv=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10).vwnd
#ws=np.sqrt((wu**2)+(wv**2))
#ws=ws.sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01'))

ws_ccmp=xr.open_dataset('processed/CCMP_ws_1deg_global.nc')
wu=ws_ccmp.uwnd
wv=ws_ccmp.vwnd

# %% Test Horizontal Divergence
#w = VectorWind(wu, wv)
#spd = w.magnitude()
#divergence = w.divergence().sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01'))
#div.mean(dim='time').plot()

wu=wu.sel(lat=slice(-20,20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01'))
wv=wv.sel(lat=slice(-20,20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01'))



# # THIS NEEDS TO BE RUN ONCE BUT CAN be memory intensive. In own file. 

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

#ws_ccmp=xr.open_dataarray('processed/CCMP_ws_1deg.nc')

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
all_dates=info.index[8:]#[36:] #2000 - 2020
old_all_dates=all_dates
neutral=all_dates.drop(cp_dates).drop(nina_dates).drop(ep_dates)
all_dates=neutral
fig=plt.figure(figsize=(19*2/2.54,23*2/2.54))#(figsize=(30,15))
sb1=7
sb2=3


sst_range=np.arange(-2.5,2.75,0.25)
ws_range=np.arange(-2,2.2,0.2)
chl_range=np.arange(-0.1,0.11,0.01)
npp_range=np.arange(-2.75,3,0.25)#None#np.arange(-0.04,0.0425,0.0025)
precip_range=np.arange(-6,7,1)
dco2_range=np.arange(-45,50,5)
co2_range=npp_range#None#np.arange(-0.04,0.045,0.005)



#%% EP




plot_basemap_row(fig,axn=1,
                 hovmol=sst.sel(time=ep_dates,method='nearest'),
                 mean=sst.sel(time=all_dates,method='nearest'),
                 units='Degrees C',
                 title='SST',
                 levs=sst_range,
             
                 cmap='RdBu_r')



plot_basemap_row(fig,axn=4,
                 hovmol=ws_ccmp.wspd.sel(time=ep_dates,method='nearest'),
                 mean=ws_ccmp.wspd.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed and direction',              
                 levs=ws_range,   
                 wu=wu.sel(time=ep_dates,method='nearest'),
                 wv=wv.sel(time=ep_dates,method='nearest'),
                 wv_all=wv.sel(time=all_dates,method='nearest'),
                 wu_all=wu.sel(time=all_dates,method='nearest'),
                 cmap='RdBu_r')

# plot_basemap_row(fig,axn=7,
#                  hovmol=divergence.sel(time=ep_dates,method='nearest'),
#                  mean=divergence.sel(time=all_dates,method='nearest'),
#                  units='m s$^{-1}$',
#                  title='Wind divergence',              
#                  levs=np.arange(-6*10**-6,6.1*10**-6,0.5*10**-6),
#                  wu=wu.sel(time=ep_dates,method='nearest'),
#                  wv=wv.sel(time=ep_dates,method='nearest'),
#                  wv_all=wv.sel(time=all_dates,method='nearest'),
#                  wu_all=wu.sel(time=all_dates,method='nearest'),
#                  cmap='RdBu_r')


plot_basemap_row(fig,axn=7,
                 hovmol=chl.sel(time=ep_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 
                 trend_conversion=1000,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=10,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=ep_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='mmol C m$^{-2}$ day$^{-1}$',
                 title='New production',
  
                 levs=npp_range,
        
                 #trend_conversion=1000,
                 cmap='RdBu_r')

# plot_basemap_row(fig,axn=13,
#                  hovmol=precip.sel(time=ep_dates,method='nearest'),
#                  mean=precip.sel(time=all_dates,method='nearest'),
#                  units='mm day$^{-1}$',
#                  title='Precipitation',
       
#                  levs=precip_range,

#                  cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=13,
                  hovmol=dco2.sel(time=ep_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
     
                  levs=dco2_range,#(200,1200,10),#(5.5,9.5,0.5),
    
                  cmap='RdBu_r')

plot_basemap_row(fig,axn=16,
                  hovmol=land.sel(time=ep_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='mmol C m$^{-2}$ day$^{-1}$',
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
                 hovmol=ws_ccmp.wspd.sel(time=cp_dates,method='nearest'),
                 mean=ws_ccmp.wspd.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed and direction',             
                 levs=ws_range,#p.arange(0,11,1),
                 wu=wu.sel(time=cp_dates,method='nearest'),
                 wv=wv.sel(time=cp_dates,method='nearest'),
                 wv_all=wv.sel(time=all_dates,method='nearest'),
                 wu_all=wu.sel(time=all_dates,method='nearest'),
                 cmap='RdBu_r')


# plot_basemap_row(fig,axn=8,
#                  hovmol=divergence.sel(time=cp_dates,method='nearest'),
#                  mean=divergence.sel(time=all_dates,method='nearest'),
#                  units='m s$^{-1}$',
#                  title='Wind divergence',             
#                  levs=np.arange(-6*10**-6,6.1*10**-6,0.5*10**-6),##p.arange(0,11,1),
#                  wu=wu.sel(time=cp_dates,method='nearest'),
#                  wv=wv.sel(time=cp_dates,method='nearest'),
#                  wv_all=wv.sel(time=all_dates,method='nearest'),
#                  wu_all=wu.sel(time=all_dates,method='nearest'),
#                  cmap='RdBu_r')


plot_basemap_row(fig,axn=8,
                 hovmol=chl.sel(time=cp_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=11,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=cp_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='mmol C m$^{-2}$ day$^{-1}$',
                 title='New production',
                 levs=npp_range,

                 cmap='RdBu_r')

# plot_basemap_row(fig,axn=14,
#                  hovmol=precip.sel(time=cp_dates,method='nearest'),
#                  mean=precip.sel(time=all_dates,method='nearest'),
#                  units='mm day$^{-1}$',
#                  title='Precipitation',
#                  levs=precip_range,

#                  cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=14,
                  hovmol=dco2.sel(time=cp_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
                  levs=dco2_range,

                  cmap='RdBu_r')

plot_basemap_row(fig,axn=17,
                  hovmol=land.sel(time=cp_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='mmol C m$^{-2}$ day$^{-1}$',
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
                 hovmol=ws_ccmp.wspd.sel(time=nina_dates,method='nearest'),
                 mean=ws_ccmp.wspd.sel(time=all_dates,method='nearest'),
                 units='m s$^{-1}$',
                 title='Wind speed and direction',
                 wu=wu.sel(time=nina_dates,method='nearest'),
                 wv=wv.sel(time=nina_dates,method='nearest'),
                 wv_all=wv.sel(time=all_dates,method='nearest'),
                 wu_all=wu.sel(time=all_dates,method='nearest'),
                 levs=ws_range,
                 cmap='RdBu_r')


# plot_basemap_row(fig,axn=9,
#                  hovmol=divergence.sel(time=nina_dates,method='nearest'),
#                  mean=divergence.sel(time=all_dates,method='nearest'),
#                  units='m s$^{-1}$',
#                  title='Wind divergence',
#                  wu=wu.sel(time=nina_dates,method='nearest'),
#                  wv=wv.sel(time=nina_dates,method='nearest'),
#                  wv_all=wv.sel(time=all_dates,method='nearest'),
#                  wu_all=wu.sel(time=all_dates,method='nearest'),
#                  levs=np.arange(-6*10**-6,6.1*10**-6,0.5*10**-6),
#                  cmap='RdBu_r')


plot_basemap_row(fig,axn=9,
                 hovmol=chl.sel(time=nina_dates,method='nearest'),
                 mean=chl.sel(time=all_dates,method='nearest'),
                 units='mg chl m$^{-3}$',
                 title='TPCA chlorophyll',
                 levs=chl_range,
                 cmap='RdBu_r')


plot_basemap_row(fig,axn=12,
                 hovmol=avg_npp.sel(lat=slice(-15,15)).sel(time=nina_dates,method='nearest'),
                 mean=avg_npp.sel(lat=slice(-15,15)).sel(time=all_dates,method='nearest'),
                 units='mmol C m$^{-2}$ day$^{-1}$',
                 title='New production',
                 levs=npp_range,

                 cmap='RdBu_r')

# plot_basemap_row(fig,axn=15,
#                  hovmol=precip.sel(time=nina_dates,method='nearest'),
#                  mean=precip.sel(time=all_dates,method='nearest'),
#                  units='mm day$^{-1}$',
#                  title='Precipitation',
                 
#                  levs=precip_range,

#                  cmap='RdBu_r')


#Delta pCO2
h=plot_basemap_row(fig,axn=15,
                  hovmol=dco2.sel(time=nina_dates,method='nearest'),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  mean=dco2.sel(time=all_dates,method='nearest'),
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
                  levs=dco2_range,
 
                  cmap='RdBu_r')

plot_basemap_row(fig,axn=18,
                  hovmol=land.sel(time=nina_dates,method='nearest'),
                  mean=land.sel(time=all_dates,method='nearest'),
                  units='mmol C m$^{-2}$ day$^{-1}$',
                  title='Air-sea CO$_{2}$ flux',
                  levs=co2_range,

                  cmap='RdBu_r')





plt.tight_layout()
plt.savefig('figs/Figure4.png',dpi=300)
plt.savefig('figs/vector/Figure4.eps')
plt.savefig('figs/vector/Figure4.pdf')
plt.show()

try:
    plt.savefig('figs/Figure4.jpeg',dpi=300)
except:
    pass
plt.show()


# Check correlation between new prod and sst
# %%
ev=[ep_dates,cp_dates,nina_dates,info.index,neutral]
for e in ev:
    sst_corr=sst.sel(time=e,method='nearest').sel(time=slice(np.datetime64('1997-09-01'),np.datetime64('2020-01-01')))#.mean(dim='time')-sst.sel(time=all_dates,method='nearest').mean(dim='time')
    avg_npp_corr=avg_npp.sel(time=e,method='nearest').sel(lat=slice(-15,15))#.mean(dim='time')-avg_npp.sel(time=all_dates,method='nearest').mean(dim='time')).sel(lat=slice(-15,15))

    startday=np.datetime64('2000-01-01')
    endday=np.datetime64('2020-01-01') 
        
    sst_corr=sst_corr.sel(time=slice(startday,endday))
    avg_npp_corr=avg_npp_corr.sel(time=slice(startday,endday))
    
    c=xr.corr(sst_corr,avg_npp_corr,dim='time').mean().values
    print(c)