#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:19:00 2020

@author: npittman
"""

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap
from scipy.signal import detrend as detrend
import xscale
import numpy as np
import cartopy.crs as ccrs
import geopandas

import geopandas as gpd
from shapely.geometry import LineString
from shapely.ops import split
from shapely.affinity import translate


def plot_basemap():
    m = Basemap(llcrnrlon=120.,llcrnrlat=-14.5,urcrnrlon=290,urcrnrlat=14.51,
                resolution='l',projection='merc',fix_aspect=False)
    m.drawcoastlines()
    m.fillcontinents()
    # draw parallels     # labels = [left,right,top,bottom]
    m.drawparallels(np.arange(-20,21,10),labels=[1,0,1,1],fontsize=12,latmax=20)
    m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],fontsize=12)
    return m


def plot_basemap_row(fig,axn,hovmol,hovmol2,units,title,units_tr,levs=None,levs_trend=None,trend_conversion=None,sb1=7,sb2=2,cmap='viridis',cmaptr='RdBu_r'):
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

    endday=np.datetime64('2019-12-01')
    
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


    cb=plt.colorbar(f,ax=ax1,fraction=fr)
    cb.set_label(units,fontsize=fs)
    cb.ax.tick_params(labelsize=fs-1)
    ax1.set_title(chr(ord('`')+axn)+') Average: '+title,fontsize=fs)
    ax1.tick_params(labelsize=fs)

    #Trends
    hovmol=hovmol.where(hovmol!=-0.9999,np.nan)
    hm=hovmol.interpolate_na(dim='time').sel(time=slice(startday,endday))
    months=hm.time
    
    hovmol2=hovmol2.where(hovmol2!=-0.9999,np.nan)
    hm2=hovmol2.interpolate_na(dim='time').sel(time=slice(startday,endday))
    months2=hm2.time
    
    dt_dates=pd.to_numeric(months.values.astype('datetime64[D]'))
    num_dates=dt_dates
    hm['time']=num_dates

    print(hm)
    print(hm2)


    #This will calculate the per pixel trends and pvalues

    time=hm.time.values
    xx=np.concatenate(hm.T)
    yy=np.concatenate(hm2.T)
    #print(xx.shape)
    tr=[]
    pv=[]
    for i in range(xx.shape[0]):
        #print(xx[i,:])
        stat=linregress(yy[i,:],xx[i,:])
        #print(stat)
        tr.append(stat.slope)
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
    
    
    

seamask=xr.open_dataset('processed/seamask.nc') #Because 2020 version doesn't have it.
seamask= seamask.assign_coords(lon=(seamask.lon % 360)).roll(lon=(seamask.dims['lon']),roll_coords=False).sortby('lon')	


#landsch_fp='processed/flux/landshutzer.nc'
#landschutzer=xr.open_dataset(landsch_fp)


#landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc'
landschutzer=xr.open_dataset(landsch_fp)
landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
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

land['time']=land.time.astype('datetime64[M]')

diff=land-avg_npp
diff1=diff.where((diff<0.1)|(diff<-0.1),np.nan)


sst = xr.open_dataset('datasets/sst/sst.mnmean.nc')
sst= sst.assign_coords(lon=(sst.lon % 360)).roll(lon=(sst.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
sst=sst.sel(lon=slice(120,290),lat=slice(20,-20)).sst
sst=sst.where(seamask.seamask==1)

#wu=xr.open_dataset('datasets/uwnd.mon.mean.nc')
#wv=xr.open_dataset('datasets/vwnd.mon.mean.nc')


startday=np.datetime64('2000-01-01')
endday=np.datetime64('2019-12-01')


newprod=avg_npp.sel(lat=slice(-15,15)).interpolate_na(dim='time').sel(time=slice(startday,endday))
co2=land.sel(lat=slice(-15,15)).interpolate_na(dim='time').sel(time=slice(startday,endday))
pco2_intrp=pco2.sel(lat=slice(-15,15)).interpolate_na(dim='time').sel(time=slice(startday,endday))
pco2_intrp['time']=pco2_intrp.time.astype('datetime64[M]')
kw['time']=kw.time.astype('datetime64[M]')

dco2['time']=dco2.time.astype('datetime64[M]')
pco2=pco2_intrp

kw1=kw.sel(lat=slice(-15,15)).interpolate_na(dim='time').sel(time=slice(startday,endday))


sst=sst.sel(lat=slice(15,-15)).interpolate_na(dim='time').sel(time=slice(startday,endday))
      


#plt.contourf(lo,la,detrend(pco2_intrp.fillna(-999)).mean(axis=0),levels=np.arange(0,500,1)),plt.colorbar()
#plt.show()

#plt.contourf(lo,la,detrend(newprod.fillna(-999)).mean(axis=0),levels=np.arange(0,1000,1)),plt.colorbar()
#plt.show()

#Detrend the variables. 
newprod_dt=xscale.signal.fitting.detrend(newprod,dim='time',type='linear') #I don't know why these methods are different. But seem to work.
#newprod_dt1=newprod.copy(data=detrend(newprod.fillna(-999)))
co2_dt=xscale.signal.fitting.detrend(co2.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.
pco2_dt=xscale.signal.fitting.detrend(pco2_intrp.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.
kw_dt=xscale.signal.fitting.detrend(kw1.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.
dco2_dt=xscale.signal.fitting.detrend(dco2.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.
sst_dt=xscale.signal.fitting.detrend(sst.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.





# co2_dt=co2.copy(data=detrend(co2.fillna(-999)))
# pco2_dt=pco2_intrp.copy(data=detrend(pco2_intrp.fillna(-999)))
# kw_dt=kw1.copy(data=detrend(kw1.fillna(-999)))
# sst_dt=sst.copy(data=detrend(sst.fillna(-999)))
# dco2_dt=dco2.copy(data=detrend(dco2.fillna(-999)))


# newprod=xscale.signal.fitting.detrend(newprod,dim='time',type='linear') #I don't know why these methods are different. But seem to work.
# pco2_intrp=pco2_intrp.copy(data=detrend(pco2_intrp.fillna(-999)))
# co2=co2.copy(data=detrend(co2.fillna(-999)))

# %%
import xesmf as xe
#wu=xr.open_dataset('datasets/uwnd.mon.mean.nc').sel(level=1000,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
#wv=xr.open_dataset('datasets/vwnd.mon.mean.nc').sel(level=1000,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd
wu=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
wv=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd
precip= xr.open_dataset('datasets/precip.mon.mean.enhanced.nc').sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).precip

ws=np.sqrt((wu**2)+(wv**2))



ws=ws.to_dataset(name='windspeed')
ws=ws.sortby(['lat'], ascending=True) #This is backwards and was doing everything upside down. 
lat_rad=land_pac.lat.diff(dim='lat').mean().values/2
lon_rad=land_pac.lon.diff(dim='lon').mean().values/2
landlats=np.linspace(land_pac.lat.min().values-lat_rad,land_pac.lat.max().values+lat_rad,len(land_pac.lat)+1)
landlons=np.linspace(land_pac.lon.min().values-lon_rad,land_pac.lon.max().values+lon_rad,len(land_pac.lon)+1)

lat_rad=ws.lat.diff(dim='lat').mean().values/2
lon_rad=ws.lon.diff(dim='lon').mean().values/2
lats=np.linspace(ws.lat.min().values-lat_rad,ws.lat.max().values+lat_rad,len(ws.lat)+1)
lons=np.linspace(ws.lon.min().values-lon_rad,ws.lon.max().values+lon_rad,len(ws.lon)+1)

x=np.meshgrid(lats,lons)
#mod=mod.expand_dims({'latb':lats})#,'lonb':lons})
#mod.coords('lon_b')=lats#=(['lonb','latb'],x[0])
ws.coords['lon_b']=lons
ws.coords['lat_b']=lats
ws['lat_b']=lats#(['lonb','latb'],x[1])
ws['lon_b']=lons#(['lonb','latb'],x[0])
ws['lat']=ws.lat
land_pac=land_pac_all
land_pac.coords['lon_b']=landlons
land_pac.coords['lat_b']=landlats
land_pac['lat_b']=landlats#(['lonb','latb'],x[1])
land_pac['lon_b']=landlons##(['lonb','latb'],x[0])


regridder = xe.Regridder(ws, land_pac, 'conservative',reuse_weights=True)
ws_1d=regridder(ws)

ws_1d=ws_1d.resample(time='M').mean(dim='time') 
ws_1d['time']=ws_1d.time.astype('datetime64[M]')#_rg #First day of month

try:
    ws_1d.to_netcdf('processed/ws_1deg.nc')
except:
    pass
ws1d=xr.open_dataset('processed/ws_1deg.nc')
#ws_dt=ws1d.windspeed.copy(data=detrend(ws1d.windspeed.fillna(-999)))
ws_dt=xscale.signal.fitting.detrend(ws1d.windspeed.chunk(None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.



prec=precip.to_dataset(name='precip')
prec=prec.sortby(['lat'], ascending=True) #This is backwards and was doing everything upside down. 

lat_rad=prec.lat.diff(dim='lat').mean().values/2
lon_rad=prec.lon.diff(dim='lon').mean().values/2
lats=np.linspace(prec.lat.min().values-lat_rad,prec.lat.max().values+lat_rad,len(prec.lat)+1)
lons=np.linspace(prec.lon.min().values-lon_rad,prec.lon.max().values+lon_rad,len(prec.lon)+1)

x=np.meshgrid(lats,lons)
#mod=mod.expand_dims({'latb':lats})#,'lonb':lons})
#mod.coords('lon_b')=lats#=(['lonb','latb'],x[0])
ws.coords['lon_b']=lons
ws.coords['lat_b']=lats
ws['lat_b']=lats#(['lonb','latb'],x[1])
ws['lon_b']=lons#(['lonb','latb'],x[0])
ws['lat']=ws.lat
land_pac=land_pac_all
land_pac.coords['lon_b']=landlons
land_pac.coords['lat_b']=landlats
land_pac['lat_b']=landlats#(['lonb','latb'],x[1])
land_pac['lon_b']=landlons##(['lonb','latb'],x[0])


x=np.meshgrid(lats,lons)
#mod=mod.expand_dims({'latb':lats})#,'lonb':lons})
#mod.coords('lon_b')=lats#=(['lonb','latb'],x[0])
prec.coords['lon_b']=lons
prec.coords['lat_b']=lats
prec['lat_b']=lats#(['lonb','latb'],x[1])
prec['lon_b']=lons#(['lonb','latb'],x[0])
prec['lat']=prec.lat




regridder = xe.Regridder(precip, land_pac, 'bilinear',reuse_weights=True)
prec_1d=regridder(prec)

prec_1d=prec_1d.resample(time='M').mean(dim='time') 
prec_1d['time']=prec_1d.time.astype('datetime64[M]')#_rg #First day of month

try:
    prec_1d.to_netcdf('processed/prec_1deg.nc')
except:
    pass
prec_1d=xr.open_dataset('processed/prec_1deg.nc')
#prec_dt=prec_1d.precip.copy(data=detrend(prec_1d.precip.fillna(-999)))
prec_dt=xscale.signal.fitting.detrend(prec_1d.precip.chunk(chunks=None),dim='time',type='linear') #I don't know why these methods are different. But seem to work.



# %This will calculate the per pixel trends and pvalues

def get_spatial_trends(var1,var2):
    #dt_dates=pd.to_numeric(var.time.values.astype('datetime64[D]'))
    #num_dates=dt_dates
    #var['time']=num_dates
    #time=var.time.values
    xx=np.concatenate(var1.T)
    yy=np.concatenate(var2.T)
    tr=[]
    pv=[]
    for i in range(xx.shape[0]):
        #print(xx[i,:])
        stat=linregress(yy[i,:],xx[i,:])
        #print(stat)
        tr.append(stat.slope)
        pv.append(stat.pvalue)
        
    tr=np.array(tr).reshape(len(var1.lon),len(var1.lat)).T
    pv=np.array(pv).reshape(len(var1.lon),len(var1.lat)).T
    
    hh=var1.copy()
    hh=hh.drop('time')
    hh['trend']=(['lat','lon'],tr)
    hh['pval']=(['lat','lon'],pv)
    return hh


def shift_geom(shift, gdataframe, plotQ=False, extent=[-60,120, -50,40]):
    # shift: expect positive values for good results
    # crs: crs="EPSG:4326" (not anything else); this handles world geometries
    # this code is adapted from answer found in SO
    # will be credited here: ??? https://stackoverflow.com/questions/65100498/changing-longitude-coordinates-of-shape-file-for-plotting/65119696#65119696
    shift -= 180
    moved_geom = []
    splitted_geom = []
    border = LineString([(shift,90),(shift,-90)])

    for row in gdataframe["geometry"]:
        splitted_geom.append(split(row, border))
    for element in splitted_geom:
        items = list(element)
        for item in items:
            minx, miny, maxx, maxy = item.bounds
            if minx >= shift:
                moved_geom.append(translate(item, xoff=-180-shift))
            else:
                moved_geom.append(translate(item, xoff=180-shift))

    # got `moved_geom` as the moved geometry
    # must specify CRS here
    moved_geom_gdf = gpd.GeoDataFrame({"geometry": moved_geom}, crs="EPSG:4326")

    # can change crs here
    if plotQ:
        #fig1, ax1 = plt.subplots(figsize=[8,6])
        fig = plt.figure( figsize=(8,8) )
        ax1 = fig.add_subplot( projection=ccrs.PlateCarree() )
        moved_geom_gdf.plot( ax=ax1, figsize=(8,5), scheme='quantiles', cmap='tab20')
        ax1.set_extent(extent, crs=ccrs.PlateCarree())
        # ax1.gridlines( draw_labels=True ) # careful, misleading labels
        plt.show()

    return moved_geom_gdf

wideseas = geopandas.GeoDataFrame.from_file('datasets/longhurst/')
wideseas180 = shift_geom(180, wideseas, False, extent=[-60,120, -50,70])
wideseas180.to_file('longhurst_fixed')
# %% FIGURE WINDSPEED
#pCO22
import matplotlib.patches as patches
from matplotlib.patches import Polygon
moorings=[165,190,205,220,235,250]
ms=20
def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), facecolor=None, alpha=0.9,edgecolor='k',linewidth=4,fill=False)
    plt.gca().add_patch(poly)


def plotpolys(m):
    lats_boxa=[-10,10,10,-10]
    lons_boxa=[150,150,210,210]
    lats_boxb=[-10,10,10,-10]
    lons_boxb=[210,210,245,245]
    lats_boxc=[10,0,0,10]
    lons_boxc=[245,245,280,280]
    lats_boxd=[-10,0,0,-10]
    lons_boxd=[245,245,280,280]
    lats_boxe=[10,14,14,10]
    lons_boxe=[150,150,260,260]
    lats_boxf=[-14,-10,-10,-14]
    lons_boxf=[150,150,240,240]

    draw_screen_poly(lats_boxa,lons_boxa,m)
    draw_screen_poly(lats_boxb,lons_boxb,m)
    draw_screen_poly(lats_boxc,lons_boxc,m)
    draw_screen_poly(lats_boxd,lons_boxd,m)
    draw_screen_poly(lats_boxe,lons_boxe,m)
    draw_screen_poly(lats_boxf,lons_boxf,m)


ws1d=xr.open_dataset('processed/ws_1deg.nc').windspeed
corr=xr.corr(pco2,ws1d,dim='time')
#corr=xr.corr(co2,newprod,dim='time')

fig=plt.figure(figsize=(20,15))
prec_dt=prec_1d.precip.copy(data=detrend(prec_1d.precip.fillna(-999)))
ws1d_dt=ws1d.copy(data=detrend(ws1d.fillna(-999)))


# NEW PRODUCTION AND WINDSPEED
corr=xr.corr(newprod,prec_1d.precip,dim='time')
#pv=get_spatial_trends(newprod,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(8,2,1)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
#plotpolys(m)
plt.title('Correlation between new production and windspeed',fontsize=16)



#CO2 and windspeed
corr=xr.corr(co2,prec_1d.precip,dim='time')
#pv=get_spatial_trends(co2,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(8,2,3)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.title('Correlation between CO2 flux and windspeed',fontsize=16)

#New Production and CO2 flux
corr=xr.corr(newprod,dco2,dim='time')
#pv=get_spatial_trends(newprod,sst).pval#.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
span=max(abs(corr.min()),corr.max())
ax1=fig.add_subplot(8,2,5)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
       #ax2.tick_params(labelsize=fs)
plt.title('Correlation between new production and dCO2',fontsize=16)


#DCO2 and Windspeed
corr=xr.corr(dco2.sel(lat=slice(-15,15)),prec_1d.precip.sel(lat=slice(-15,15)),dim='time')
#pv=get_spatial_trends(co2,newprod).pval
ax2=fig.add_subplot(8,2,7)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between dCO2 and windspeed',fontsize=16)



#New Production and CO2 flux
corr=xr.corr(newprod,prec_1d.precip,dim='time')
#pv=get_spatial_trends(newprod,sst).pval#.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
span=max(abs(corr.min()),corr.max())
ax1=fig.add_subplot(8,2,9)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
       #ax2.tick_params(labelsize=fs)
plt.title('Correlation between new production and precipitation',fontsize=16)


#DCO2 and Windspeed
corr=xr.corr(ws1d,prec_1d.precip,dim='time')
#pv=get_spatial_trends(co2,newprod).pval
ax2=fig.add_subplot(8,2,11)
m=plot_basemap()
lo,la=np.meshgrid(ws1d.lon.values,ws1d.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between dCO2 and precipitation ', fontsize=16)


#CHECK
# plt.plot(ws1d.sel(lat=slice(-15,15)).mean(dim='time'),newprod.mean(dim='time'))



# NEW PRODUCTION AND WINDSPEED
corr=xr.corr(newprod_dt,ws1d_dt,dim='time')
#pv=get_spatial_trends(newprod,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(8,2,2)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
#plotpolys(m)
plt.title('Correlation between new production and windspeed DT',fontsize=16)



#CO2 and windspeed
corr=xr.corr(co2_dt,ws1d_dt,dim='time')
#pv=get_spatial_trends(co2,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(8,2,4)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.title('Correlation between CO2 flux and windspeed DT',fontsize=16)

#New Production and CO2 flux
corr=xr.corr(newprod_dt,dco2_dt,dim='time')
#pv=get_spatial_trends(newprod,sst).pval#.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
span=max(abs(corr.min()),corr.max())
ax1=fig.add_subplot(8,2,6)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
       #ax2.tick_params(labelsize=fs)
plt.title('Correlation between new production and dCO2 DT',fontsize=16)


#DCO2 and Windspeed
corr=xr.corr(dco2_dt.sel(lat=slice(-15,15)),ws1d_dt.sel(lat=slice(-15,15)),dim='time')
#pv=get_spatial_trends(co2,newprod).pval
ax2=fig.add_subplot(8,2,8)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between dCO2 and windspeed DT', fontsize=16)


#New Production and CO2 flux
corr=xr.corr(newprod_dt,prec_dt,dim='time')
#pv=get_spatial_trends(newprod,sst).pval#.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
span=max(abs(corr.min()),corr.max())
ax1=fig.add_subplot(8,2,10)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
       #ax2.tick_params(labelsize=fs)
plt.title('Correlation between new production and precipitation Detrended',fontsize=16)


#DCO2 and Windspeed
corr=xr.corr(ws1d_dt,prec_dt,dim='time')
#pv=get_spatial_trends(co2,newprod).pval
ax2=fig.add_subplot(8,2,12)
m=plot_basemap()
lo,la=np.meshgrid(ws1d.lon.values,ws1d.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
for x in moorings:
       x1,y1=m(x,0)
       m.plot(x1,y1,marker='x',c='k',markersize=ms)
       
plt.colorbar()
#plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between dCO2 and precipitation Detrended ', fontsize=16)



plt.tight_layout()
plt.show()




import sys
sys.exit()

# %% Detrended 


fig=plt.figure(figsize=(20,12))


corr=xr.corr(newprod_dt,ws_dt,dim='time')
#pv=get_spatial_trends(newprod,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(4,2,1)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
plotpolys(m)
plt.title('Correlation between new production and windspeed DETRENDED',fontsize=16)


corr=xr.corr(newprod_dt,sst_dt,dim='time')
#pv=get_spatial_trends(newprod,sst).pval#.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
span=max(abs(corr.min()),corr.max())
ax1=fig.add_subplot(4,2,2)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between new production and sst',fontsize=16)


corr=xr.corr(co2_dt,kw_dt,dim='time')
#pv=get_spatial_trends(co2,ws1d.sel(time=slice(startday,endday),lat=slice(-15,15))).pval
ax2=fig.add_subplot(4,2,3)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between CO2 flux and windspeed',fontsize=16)



corr=xr.corr(co2_dt,sst_dt,dim='time')
#pv=get_spatial_trends(co2,sst).pval
ax1=fig.add_subplot(4,2,4)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(corr.min()),corr.max()).round(2)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#),levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])
plotpolys(m)
m.readshapefile('datasets/longhurst/Longhurst_world_v4_2010',name='longhurst')
plt.title('Correlation between CO2 and SST',fontsize=16)

corr=xr.corr(co2_dt,newprod_dt,dim='time')
#pv=get_spatial_trends(co2,newprod).pval
ax2=fig.add_subplot(4,2,5)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between newprod and CO2',fontsize=16)


corr=xr.corr(pco2_dt,newprod_dt,dim='time')
#pv=get_spatial_trends(pco2,newprod.sel(time=slice(startday,endday-1))).pval
ax2=fig.add_subplot(4,2,6)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between newprod and pco2',fontsize=16)

corr=xr.corr(dco2_dt,newprod_dt,dim='time')
#pv=get_spatial_trends(dco2.sel(lat=slice(-15,15)).sel(time=slice(startday,endday)),newprod).pval
ax2=fig.add_subplot(4,2,7)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between newprod and dco2',fontsize=16)

corr=xr.corr(dco2_dt.sel(lat=slice(-15,15)),ws_dt.sel(lat=slice(-15,15)),dim='time')
#pv=get_spatial_trends(dco2.sel(lat=slice(-15,15),time=slice(startday,endday)),ws1d.sel(lat=slice(-15,15),time=slice(startday,endday))).pval
ax2=fig.add_subplot(4,2,8)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.25))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plotpolys(m)
#m.contourf(lo1,la1,pv,colors='none',hatches=['.'],levels=[0,0.05])

plt.title('Correlation between dCO2 and windspeed',fontsize=16)




plt.tight_layout()
plt.show()


import sys
sys.exit()

# %%
#pCO22
covariance=xr.cov(co2,pco2_intrp,dim='time')
corr=xr.corr(co2,pco2_intrp,dim='time')
#corr=xr.corr(co2,newprod,dim='time')

fig=plt.figure(figsize=(20,12))
ax1=fig.add_subplot(3,2,1)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)

span=max(abs(covariance.min()),covariance.max()).round(3)
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))#,levels=np.arange(-0.001,0.0012,0.0002))
#m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-0.0015,0.0016,0.0001),extend='max')#span+(span/20),span/20))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plt.title('Covariance between pCO2 and new production',fontsize=16)

ax2=fig.add_subplot(3,2,2)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between pCO2 and new production',fontsize=16)


covariance=xr.cov(newprod,sst,dim='time')
corr=xr.corr(newprod,sst,dim='time')
span=max(abs(covariance.min()),covariance.max())
ax1=fig.add_subplot(3,2,3)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(covariance.min()),covariance.max()).round(2)
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
plt.title('Covariance between new production and SST',fontsize=16)

ax2=fig.add_subplot(3,2,4)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between new production and SST',fontsize=16)

kw1=kw.sel(lat=slice(-15,15)).interpolate_na(dim='time').sel(time=slice(startday,endday))

covariance=xr.cov(newprod,kw1,dim='time')
corr=xr.corr(newprod,kw1,dim='time')

ax1=fig.add_subplot(3,2,5)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(covariance.min()),covariance.max()).round(2)
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
plt.title('Covariance between new prod and solubility',fontsize=16)


ax2=fig.add_subplot(3,2,6)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plt.title('Correlation between new prod and solubility',fontsize=16)


plt.tight_layout()
plt.show()

import sys
sys.exit()
# %%
#CO2 FLUXXXXX
covariance=xr.cov(newprod,co2,dim='time')
corr=xr.corr(newprod,co2,dim='time')
#corr=xr.corr(co2,newprod,dim='time')

fig=plt.figure(figsize=(20,12))
ax1=fig.add_subplot(3,2,1)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)

span=max(abs(covariance.min()),covariance.max()).round(3)
#m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))#,levels=np.arange(-0.001,0.0012,0.0002))
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-0.0015,0.0016,0.0001),extend='max')#span+(span/20),span/20))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plt.title('Covariance between air-sea CO2 flux and new production',fontsize=16)

ax2=fig.add_subplot(3,2,2)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between air-sea CO2 flux and new production',fontsize=16)


covariance=xr.cov(newprod,sst,dim='time')
corr=xr.corr(newprod,sst,dim='time')
span=max(abs(covariance.min()),covariance.max())
ax1=fig.add_subplot(3,2,3)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(covariance.min()),covariance.max()).round(2)
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
plt.title('Covariance between new production and SST',fontsize=16)

ax2=fig.add_subplot(3,2,4)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between new production and SST',fontsize=16)



covariance=xr.cov(co2,sst,dim='time')
corr=xr.corr(co2,sst,dim='time')

ax1=fig.add_subplot(3,2,5)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
span=max(abs(covariance.min()),covariance.max()).round(2)
m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-span,span+(span/20),span/20))
plt.colorbar()
plt.title('Covariance between CO2 and SST',fontsize=16)


ax2=fig.add_subplot(3,2,6)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plt.title('Correlation between CO2 and SST',fontsize=16)


plt.tight_layout()
plt.show()



# %% Trying different idea; Maps / corrleations during CP, EP and Nina events
lanina=pd.read_csv('processed/indexes/la_nina_events.csv')
cp_nino=pd.read_csv('processed/indexes/cp_events.csv')
ep_nino=pd.read_csv('processed/indexes/ep_events.csv')


fp='processed/combined_dataset/month_data_exports.nc'
info=xr.open_mfdataset(fp).sel(Mooring=195).to_dataframe()

nina=pd.DataFrame()
ep=pd.DataFrame()
cp=pd.DataFrame()
for i in lanina.iterrows(): nina=nina.append(info[slice(i[1].start,i[1].end)])
for i in ep_nino.iterrows(): ep=ep.append(info[slice(i[1].start,i[1].end)])
for i in cp_nino.iterrows(): cp=cp.append(info[slice(i[1].start,i[1].end)])
nina_dates=nina.index
ep_dates=ep.index
cp_dates=cp.index


newprod.sel(time=cp_dates).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()
co2.sel(time=cp_dates).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()
newprod.sel(time=ep_dates[13:]).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()
co2.sel(time=ep_dates[13:]).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()
newprod.sel(time=nina_dates[18:]).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()
co2.sel(time=nina_dates[18:]).mean(dim='time').plot(vmin=0,vmax=0.15)
plt.show()


# %%
###
"""
Correlations SPlit by ENSO


"""
#newprod_backup=newprod

#newprod=avg_npp*f_ratios.trim
#newprod=avg_npp*f_ratios.dunne2005
#newprod=avg_npp*f_ratios.trim
newprod=avg_npp*f_ratios.laws2011a


#either new prod, co2 or do2
dco2['time']=dco2.time.astype('datetime64[M]')
cp_corr=xr.corr(newprod.sel(time=cp_dates),dco2.sel(time=cp_dates),dim='time')
ep_corr=xr.corr(newprod.sel(time=ep_dates[13:]),dco2.sel(time=ep_dates[13:]),dim='time')
nina_corr=xr.corr(newprod.sel(time=nina_dates[18:]),dco2.sel(time=nina_dates[18:]),dim='time')





fig=plt.figure(figsize=(20,12))


corr=xr.corr(dco2,newprod,dim='time')
ax2=fig.add_subplot(5,2,1)
m=plot_basemap()
lo,la=np.meshgrid(newprod.lon.values,newprod.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between air-sea flux and new production (All time) Laws2011a',fontsize=16)


corr_dt=xr.corr(co2_dt,newprod_dt,dim='time')
ax2=fig.add_subplot(5,2,3)
m=plot_basemap()
lo,la=np.meshgrid(corr_dt.lon.values,corr_dt.lat.values)
lo1,la1=m(lo,la)
m.contourf(lo1,la1,corr_dt,cmap='bwr',levels=np.arange(-1,1.1,0.1))#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between air-sea flux and new production (All time: Detrended)',fontsize=16)



ax1=fig.add_subplot(5,2,5)
m=plot_basemap()
lo,la=np.meshgrid(cp_corr.lon.values,cp_corr.lat.values)
lo1,la1=m(lo,la)
span=np.arange(-1,1.1,0.1)
    
m.contourf(lo1,la1,cp_corr,cmap='bwr',levels=span)#,levels=np.arange(-0.001,0.0012,0.0002))
#m.contourf(lo1,la1,covariance,cmap='bwr',levels=np.arange(-0.0015,0.0016,0.0001),extend='max')#span+(span/20),span/20))#,levels=np.arange(-0.001,0.0012,0.0002))

plt.colorbar()
plt.title('Correlation between CO2 flux and new production: CP El Nino',fontsize=16)


#corr=xr.corr(co2,newprod,dim='time')
ax2=fig.add_subplot(5,2,7)
m=plot_basemap()
lo,la=np.meshgrid(ep_corr.lon.values,ep_corr.lat.values)
lo1,la1=m(lo,la)

m.contourf(lo1,la1,ep_corr,cmap='bwr',levels=span)
plt.colorbar()
plt.title('Correlation between CO2 flux and new production: EP El Nino',fontsize=16)

corr=xr.corr(co2,newprod,dim='time')
ax2=fig.add_subplot(5,2,9)
m=plot_basemap()
lo,la=np.meshgrid(nina_corr.lon.values,nina_corr.lat.values)
lo1,la1=m(lo,la)

m.contourf(lo1,la1,nina_corr,cmap='bwr',levels=span)#,levels=np.arange(-0.001,0.0012,0.0002))
plt.colorbar()
plt.title('Correlation between CO2 flux and new production: La Nina',fontsize=16)
plt.tight_layout()
plt.show()
# %%
#covariance.plot(vmin=-0.001,vmax=0.001,cmap='bwr')

#corr.plot()
#corr.show()



# fig=plt.figure(figsize=(19*2/2.54,23*2/2.54))#(figsize=(30,15))
# sb1=1
# sb2=1

# plot_basemap_row(fig,axn=1,
#                  hovmol=avg_npp.sel(lat=slice(-15,15)),
#                  hovmol2=land.sel(lat=slice(-15,15)),
#                  units='gC m$^{-2}$ day$^{-1}$',
#                  title='New production',
#                  units_tr='mgC m$^{-2}$ day$^{-1}$ year$^{-1}$',
#                  levs=np.arange(0,0.26,0.025),
#                  levs_trend=np.arange(-200,210,5),
#                  trend_conversion=1000,
#                  cmap='viridis')    
    
# def trends(x,y,c):  
#     from scipy.stats import linregress
#     mean=np.nanmean(y)
#     std=np.nanstd(y)*1

#     mask=~np.isnan(x)
#     x=x[mask]
#     y=y[mask]
#     mask1=~np.isnan(y)
#     x=x[mask1]
#     y=y[mask1]

#     #x_n=np.arange(0,len(x))
#     # x1=np.arange(np.datetime64(x[0],'M'),np.datetime64(x[-1],'M')+np.timedelta64(1,'M'))
#     #x1=trd.index.values.astype('datetime64[D]')
#     #x1=x.values.astype('datetime64[D]')
#     #pd.to_numeric(x1)
#     slope, intercept, r_value, p_value,std_err = linregress(x,y)
#     mn=min(x)
#     mx=max(x)
#     x1=np.linspace(mn,mx,len(x))
#     y1=slope*x1+intercept
    
#     plt.plot(x1,y1,c,linestyle='--',linewidth=2.5)  
#     #ax.text(x1[-1]-(x1[-1]*0.1),y1[-1]-(y1[-1]*0.1),'R2='+str(np.round(r_value**2,3)))
#     print( slope, intercept, r_value**2,p_value,std_err)
#     return slope, intercept, r_value,p_value,std_err
    
# # plt.figure(figsize=(20,10))
# plt.subplot(321)
# fp='processed/combined_dataset/month_data_exports.nc'
# dat=xr.open_mfdataset(fp)
# dat.windspeed.groupby('Date.month').mean().plot()
# plt.title('Windspeed seasonal hovmoller')
# #plt.show()

# plt.subplot(343)
# plt.scatter(dat.delta_pCO2,(dat.co2flux4_land_gmyr/365)*1000,label='delta pCO2 (uatm)',alpha=0.45)
# plt.xlabel('pCO2 (uatm')
# plt.ylabel('CO2 flux (mgC/day)')
# plt.title('pCO2 vs CO2 flux')
# print('pCO2 vs CO2 flux')
# trends(dat.delta_pCO2.values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='b')

# plt.subplot(344)
# plt.scatter(dat.delta_pCO2,dat.laws2011a*dat.cafe,label='delta pCO2 (uatm)',alpha=0.45)
# plt.ylabel('New Production (mgC/day)')
# plt.xlabel('pCO2 (uatm')
# plt.title('pCO2 vs New Production')
# print('pCO2 vs New Production')
# trends((dat.laws2011a*dat.cafe).values.flatten(),((dat.co2flux4_land_gmyr/365)*1000).values.flatten(),c='b')
