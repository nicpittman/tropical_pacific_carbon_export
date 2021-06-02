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
import matplotlib.patches as patches
from matplotlib.patches import Polygon
    
def plot_basemap():
    m = Basemap(llcrnrlon=120.,llcrnrlat=-15,urcrnrlon=290,urcrnrlat=15.01,
                resolution='l',projection='merc',fix_aspect=False)
    m.drawcoastlines()
    m.fillcontinents()
    # draw parallels     # labels = [left,right,top,bottom]
    m.drawparallels(np.arange(-20,21,10),labels=[1,0,1,1],fontsize=12,latmax=20)
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],fontsize=12)
    return m



def plot_basemap_row(fig,axn,hovmol,units,title,units_tr,levs=None,levs_trend=None,trend_conversion=None,sb1=7,sb2=2,cmap='viridis',cmaptr='RdBu_r',wu=None,wv=None):
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
        #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
        for c in f.collections:
            c.set_edgecolor("face")
    else:
        if title=='TPCA Chlorophyll':
            f=m.contourf(lo1,la1,hovmol.mean(dim='time'),extend='max',cmap=cmap,levels=levs) #11 colors
                #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
            for c in f.collections:
                c.set_edgecolor("face")
        else:
            f=m.contourf(lo1,la1,hovmol.mean(dim='time'),cmap=cmap,levels=levs) #11 colors
            #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
            for c in f.collections:
                c.set_edgecolor("face")
    ax1.axhline(0,c='k',linestyle=':')

    moorings=[165,190,205,220,235,250]
    for x in moorings:
        x1,y1=m(x,0)
        ax1.plot(x1,y1,marker='x',c='k',markersize=ms)
    
    if title=='SST':
        
        lev=28.5#29.2 #rather than 28.5
        early_sst=hovmol.sel(time=slice('1997-01-01','2002-01-01')).mean(dim='time')#.where(co2.seamask==1)
        late_sst=hovmol.sel(time=slice('2015-01-01','2020-01-01')).mean(dim='time')#.where(co2.seamask==1)
       
        m.contour(lo1,la1,early_sst,levels=[lev],linestyles='dotted',colors='k')
        
        m.contour(lo1,la1,late_sst,levels=[lev],linestyles='solid',colors='k')
        m.contour(lo1,la1,hovmol.mean(dim='time'),levels=[25],linestyles='dashed',colors='k')
        
    if title=='Air-sea CO$_{2}$ flux':
        lo_rect,la_rect=m(165,-15)
    
        def draw_screen_poly( lats, lons, m):
            #https://stackoverflow.com/questions/12251189/how-to-draw-rectangles-on-a-basemap
            x, y = m( lons, lats )
            xy = zip(x,y)
            poly = Polygon(list(xy), facecolor='none',edgecolor='k', linewidth=2)
            plt.gca().add_patch(poly)
        
        lats = [ -14.5, 14.5, 14.5, -14.5 ]
        lons = [ 165, 165, 180, 180 ]
        draw_screen_poly( lats, lons, m )
        
        lats = [ -14.5, 14.5, 14.5, -14.5 ]
        lons = [ 205, 205, 220, 220 ]
        draw_screen_poly( lats, lons, m )
        
        lats = [ -14.5, 14.5, 14.5, -14.5 ]
        lons = [ 235, 235, 250, 250 ]
        draw_screen_poly( lats, lons, m )
        
        
    #wu['lon'],wu['lat']=m(lo,la,wu.lon.values,wu.lat.values)
    #No windspeed vectors now
    #if title=='Wind speed':
    #      skip=(slice(None,None,4),slice(None,None,4)) #2 for NCEP 2
    #      m.quiver(lo1[skip],la1[skip],wu.mean(dim='time')[skip]/2,wv.mean(dim='time')[skip]/2,scale=90,headwidth=4.5)#,minshaft=2)


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
        #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
        for c in f.collections:
            c.set_edgecolor("face")
    else:
        f=m1.contourf(lo1,la1,hh.trend,cmap=cmaptr,extend='both',levels=levs_trend) #11 colors 
        #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
        for c in f.collections:
            c.set_edgecolor("face")
    cb=plt.colorbar(f,ax=ax2,extend='both',fraction=fr)
        
    cnt=m1.contourf(lo1,la1,hh.pval,colors='none',hatches=['.'],levels=[0,0.05])
    #Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
    for c in cnt.collections:
        c.set_edgecolor("face")
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

#Old NCEP2 winds
#wu=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
#wv=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd
dco2['time']=dco2.time.astype('datetime64[M]')

#ws=np.sqrt((wu**2)+(wv**2))


precip= xr.open_dataset('datasets/precip.mon.mean.enhanced.nc').sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).precip


# import seaborn as sns
# precip1= xr.open_dataset('processed/prec_1deg.nc').sel(lat=slice(-20,20),lon=slice(120,290),time=slice('1997-07-01','2019-12-31')).precip
# dco21=dco2.sel(time=slice('1997-07-01','2020-01-01'))
# prec=precip1.values.reshape(270*40*170)
# dc=dco21.values.reshape(270*40*170)

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
all_dates=info.index[36:] #2000 - 2020

#So we can select which time range. Select nina/ep.cp/ep dates 
ensodates=all_dates
etype=''


fig=plt.figure(figsize=(19*2/2.54,23*2/2.54))#(figsize=(30,15))
sb1=7
sb2=2


plot_basemap_row(fig,axn=1,
                 hovmol=sst.sel(time=ensodates),
                 units='Degrees C',
                 title='SST',
                 units_tr='Degrees C year$^{-1}$',
                 levs=np.arange(20,32,1),
                 
                 levs_trend=np.arange(-0.06,0.07,0.01),
                 #trend_conversion=1000,
                 cmap='viridis')


#wind #_ccmp
if etype=='':
    ws_ccmp_d=ensodates[:-8]
plot_basemap_row(fig,axn=3,
                 hovmol=ws_ccmp.sel(time=ws_ccmp_d),
                 units='m s$^{-1}$',
                 title='Wind speed',
                 units_tr='m s$^{-1}$ year$^{-1}$',                 
                 levs=np.arange(0,11,1),
                 levs_trend=np.arange(-0.1,0.125,0.025),
                 
                 #levs_trend=np.arange(-0.15,0.175,0.025),
                 #trend_conversion=1000,
                 cmap='viridis')

if etype=='':
    chl_d=chl.time#ensodates[:-5] #Ok so this is actually from 1997 but doesn't change anything except fills in missing trend due to strange start years
plot_basemap_row(fig,axn=5,
                 hovmol=chl.sel(time=chl_d),
                 units='mg chl m$^{-3}$',
                 title='TPCA chlorophyll',
                 units_tr='ug chl m$^{-3}$ year$^{-1}$',
                 levs=np.arange(0,0.65,0.05),
                 
                 levs_trend=np.arange(-4,4.1,1),
                 trend_conversion=1000,
                 cmap='viridis')


plot_basemap_row(fig,axn=7,
                 hovmol=avg_npp.sel(lat=slice(-15,15),time=ensodates),
                 units='gC m$^{-2}$ day$^{-1}$',
                 title='New production',
                 units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',
                 levs=np.arange(0,0.26,0.025),
                 levs_trend=np.arange(-2,2.1,0.25),
                 trend_conversion=1000,
                 cmap='viridis')



plot_basemap_row(fig,axn=9,
                 hovmol=precip.sel(time=ensodates),
                 units='mm day$^{-1}$',
                 title='Precipitation',
                 units_tr='mm day$^{-1}$ year$^{-1}$',
                 levs=np.arange(0,13,1),
                 levs_trend=np.arange(-0.08,0.085,0.005),
                 #trend_conversion=1000,
                 cmap='viridis')




#Delta pCO2
h=plot_basemap_row(fig,axn=11,
                  hovmol=dco2.sel(time=ensodates),#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
                  units='μatm',
                  title='\u0394pCO$_{2}$',#'pCO21',
                  units_tr='μatm year$^{-1}$',
                  levs=np.arange(-15,121,10),#(200,1200,10),#(5.5,9.5,0.5),
                  levs_trend=np.arange(-2.5,2.6,0.1),
                  trend_conversion=1,#1000,
                  cmap='viridis',
                  cmaptr='RdBu_r')#'Reds')


CO2_tr=plot_basemap_row(fig,axn=13,
                  hovmol=land.sel(time=ensodates),
                  units='gC m$^{-2}$ day$^{-1}$',
                  title='Air-sea CO$_{2}$ flux',
                  units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',
                  levs=np.arange(-0.14,0.15,0.02),
                  levs_trend=np.arange(-2,2.1,0.5),
                  trend_conversion=1000,
                  cmap='RdBu_r')

#Save the CO2 trends so we can superimpose onto plot later.
CO2_tr.to_netcdf('processed/results/CO2f_trend_pval_alltime.nc')

# #Ocean pCO2
# plot_basemap_row(fig,axn=13,
#                   hovmol=pco2,#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
#                   units='ppm',
#                   title='ocean pco2',#'pCO21',
#                   units_tr='ppm year$^{-1}$',
#                   levs=np.arange(300,550,10),#(200,1200,10),#(5.5,9.5,0.5),
#                   levs_trend=np.arange(-4,4,0.01),
#                   trend_conversion=1,#1000,
#                   cmap='viridis',
#                   cmaptr='RdBu_r')#'Reds')



# #Gas transfer velocity
# plot_basemap_row(fig,axn=11,
#                   hovmol=kw,#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
#                   units='ppm',
#                   title='\u0394pCO$_{2}$',#'pCO21',
#                   units_tr='ppm year$^{-1}$',
#                   levs=np.arange(500,3000,10),#(200,1200,10),#(5.5,9.5,0.5),
#                   levs_trend=np.arange(-30,30,1),
#                   trend_conversion=1,#1000,
#                   cmap='viridis',
#                   cmaptr='RdBu_r')#'Reds')

# #Atmospheric pCO2
# plot_basemap_row(fig,axn=11,
#                   hovmol=atm pco2,#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
#                   units='ppm',
#                   title='atm co2',#'pCO21',
#                   units_tr='ppm year$^{-1}$',
#                   levs=np.arange(300,500,10),#(200,1200,10),#(5.5,9.5,0.5),
#                   levs_trend=np.arange(4,6,0.01),
#                   trend_conversion=1,#1000,
#                   cmap='viridis',
#                   cmaptr='RdBu_r')#'Reds')


# plot_basemap_row(fig,axn=11,
#                   hovmol=integratedpCO2,#monthlyPCO2*1000,
#                   units='gC m$^{-2}$',
#                   title='pCO2t',
#                   units_tr='mgC m$^{-2}$ year$^{-1}$',
#                   levs=np.arange(5.5,9.5,0.5),
#                   levs_trend=np.arange(0,100,10),
#                   trend_conversion=1000,
#                   cmap='viridis',
#                   cmaptr='Reds')


# #F-ratio - Need to comment out integratedpCO2 above.
# plot_basemap_row(fig,axn=11,
#                   hovmol=ratio,#npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
#                   units='gC m$^{-2}$',
#                   title='fratio',#'pCO21',
#                   units_tr='ppm year$^{-1}$',
#                   levs=np.arange(0.05,0.35,0.01),#(200,1200,10),#(5.5,9.5,0.5),
#                   levs_trend=np.arange(-0.1,0.1,0.001),
#                   trend_conversion=100,#1000,
#                   cmap='viridis',
#                   cmaptr='RdBu_r')

# #CAFE NPP mean
# hh=plot_basemap_row(fig,axn=11,
#                   hovmol=npp1.avg_npp,#dco2,#integratedpCO2,#monthlyPCO2*1000,
#                   units='gC m$^{-2}$',
#                   title='NPP Average',
#                   units_tr='ppm year$^{-1}$',
#                   levs=np.arange(200,1200,10),#(5.5,9.5,0.5),
#                   levs_trend=np.arange(-15,15,1),
#                   trend_conversion=1,#1000,
#                   cmap='viridis',
#                   cmaptr='RdBu_r')



plt.tight_layout()
plt.savefig('figs/Figure3'+etype+'.png',dpi=300)
plt.savefig('figs/vector/Figure3'+etype+'.eps')
plt.savefig('figs/vector/Figure3'+etype+'.pdf')

try:
    plt.savefig('figs/Figure3.jpeg',dpi=300)
except:
    pass
plt.show()
