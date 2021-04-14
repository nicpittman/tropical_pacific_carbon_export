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
ep_dates=ep.index[4:]
cp_dates=cp.index



# %% Load Useful files


seamask=xr.open_dataset('processed/seamask.nc') #Because 2020 version doesn't have it.
seamask= seamask.assign_coords(lon=(seamask.lon % 360)).roll(lon=(seamask.dims['lon']),roll_coords=False).sortby('lon')	
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




#startday=np.datetime64('2000-01-01')
#endday=np.datetime64('2019-12-01')

#wu=xr.open_dataset('datasets/uwnd.mon.mean.nc').sel(level=1000,lat=0,lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
#wv=xr.open_dataset('datasets/vwnd.mon.mean.nc').sel(level=1000,lat=0,lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd

wu=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).uwnd
wv=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).vwnd

ws_ncep2=np.sqrt((wu**2)+(wv**2))

#CHeck line 164 depending if using NCEP2 or windspeed
ws=xr.open_dataarray('datasets/CCMP_windspeed.nc')

#wind=uw.sel(lat=)
precip= xr.open_dataset('datasets/precip.mon.mean.enhanced.nc').sel(lat=slice(20,-20),lon=slice(120,290),time=slice('1997-07-01','2020-01-01')).precip


newprod=avg_npp.sel(lat=slice(-15,15))#.interpolate_na(dim='time').sel(time=slice(startday,endday))
co2=land.sel(lat=slice(-15,15))#.interpolate_na(dim='time').sel(time=slice(startday,endday))
pco2=pco2.sel(lat=slice(-15,15))#.interpolate_na(dim='time').sel(time=slice(startday,endday))
pco2['time']=pco2.time.astype('datetime64[M]')
kw['time']=kw.time.astype('datetime64[M]')
dco2['time']=dco2.time.astype('datetime64[M]')

#pco2=pco2_intrp

kw1=kw.sel(lat=slice(-15,15))#.interpolate_na(dim='time').sel(time=slice(startday,endday))
sst=sst.sel(lat=slice(15,-15))#.interpolate_na(dim='time').sel(time=slice(startday,endday))

chl=xr.open_dataset('processed/flux/tpca.nc').tpca#'sw_month.nc')
chl['time']=chl.time.astype('datetime64[M]')




# %%
#Test windspeed 
minlat=-2
maxlat=2
#enso=pd.read_csv('datasets/indexes/meiv2.csv',index_col='Year')
enso=pd.read_csv('datasets/indexes/meiv2.csv',index_col=0,header=None)
enso=enso.iloc[3:] #Just so Both EMI and MEI start in 1981-01-01
enso_flat=enso.stack()
enso_dates=pd.date_range('1982','2020-07-01',freq='M')- pd.offsets.MonthBegin(1) #Probably want to check this is correct if updating.
enso_timeseries=pd.DataFrame({'Date':enso_dates,'mei':enso_flat})
enso_timeseries=enso_timeseries.where((enso_timeseries.Date>np.datetime64('1997-01-01'))&(enso_timeseries.Date<np.datetime64('2000-01-01'))).dropna()


plt.figure(figsize=(30,10))
plt.subplot(151)
wu1=xr.open_dataset('datasets/uwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290)).uwnd
d=wu1.sel(lat=slice(maxlat,minlat),time=slice(np.datetime64('1997-01-01'),np.datetime64('2000-01-01'))).mean(dim='lat')
d.plot.contourf(vmin=-6,vmax=6)
gc=plt.gca()
gc.invert_yaxis()
positions = (140, 160, 180,200,220,240,260)
labels = ("140$^\circ$E", "160$^\circ$E", "180$^\circ$",'160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W')
plt.xticks(positions, labels)
plt.title('zonal (uwnd) speed')
plt.xlim(150,265)

plt.subplot(152)
wv1=xr.open_dataset('datasets/vwnd.10m.mon.mean.nc').sel(level=10,lat=slice(20,-20),lon=slice(120,290)).vwnd
d=wv1.sel(lat=slice(maxlat,minlat),time=slice(np.datetime64('1997-01-01'),np.datetime64('2000-01-01'))).mean(dim='lat')
d.plot.contourf(vmin=-6,vmax=6)
gc=plt.gca()
gc.invert_yaxis()
positions = (140, 160, 180,200,220,240,260)
labels = ("140$^\circ$E", "160$^\circ$E", "180$^\circ$",'160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W')
plt.xticks(positions, labels)
plt.xlim(150,265)

plt.title('meridional (vwnd) speed')


plt.subplot(153)
nws=ws.sel(time=slice(np.datetime64('1997-01-01'),np.datetime64('2000-01-01')),lat=slice(minlat,maxlat)).mean(dim='lat')
nws.plot.contourf()
gc=plt.gca()
gc.invert_yaxis()
positions = (140, 160, 180,200,220,240,260)
labels = ("140$^\circ$E", "160$^\circ$E", "180$^\circ$",'160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W')
plt.xticks(positions, labels)
plt.xlim(150,265)

plt.title('Wind Speed (m/s)')


plt.subplot(154)
nws=ws.sel(time=slice(np.datetime64('1997-01-01'),np.datetime64('2000-01-01')),lat=slice(minlat,maxlat)).mean(dim='lat')
(nws-nws.mean(dim='time')).plot.contourf(cmap='bwr')
gc=plt.gca()
gc.invert_yaxis()
positions = (140, 160, 180,200,220,240,260)
labels = ("140$^\circ$E", "160$^\circ$E", "180$^\circ$",'160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W')
plt.xticks(positions, labels)
plt.title('Wind Speed Anomaly (m/s, wind-wind mean)')
plt.xlim(150,265)


plt.subplot(155)
plt.title('MEI')
plt.plot(enso_timeseries.mei,enso_timeseries.Date,linewidth=5,c='k')
gc=plt.gca()
gc.invert_yaxis()

plt.tight_layout()
# %% Plot them
latbnd=1
co2_eq=co2.sel(lat=slice(-latbnd,latbnd)).mean(dim='lat')
co2_ep=co2_eq.sel(time=ep_dates).mean(dim='time')
co2_cp=co2_eq.sel(time=cp_dates).mean(dim='time')
co2_nina=co2_eq.sel(time=nina_dates).mean(dim='time')
co2_neutral=co2_eq.drop_sel(time=ep_dates.append(cp_dates).append(nina_dates)).mean(dim='time')

dco2_eq=dco2.sel(lat=slice(-latbnd,latbnd)).mean(dim='lat')
dco2_ep=dco2_eq.sel(time=ep_dates).mean(dim='time')
dco2_cp=dco2_eq.sel(time=cp_dates).mean(dim='time')
dco2_nina=dco2_eq.sel(time=nina_dates).mean(dim='time')
dco2_neutral=dco2_eq.drop_sel(time=ep_dates.append(cp_dates).append(nina_dates)).mean(dim='time')

np_eq=newprod.sel(lat=slice(-latbnd,latbnd)).mean(dim='lat')
np_ep=np_eq.sel(time=ep_dates).mean(dim='time')
np_cp=np_eq.sel(time=cp_dates).mean(dim='time')
np_nina=np_eq.sel(time=nina_dates).mean(dim='time')
np_neutral=np_eq.drop_sel(time=ep_dates.append(cp_dates).append(nina_dates)).mean(dim='time')

sst_eq=sst.sel(lat=slice(latbnd,-latbnd)).mean(dim='lat')
sst_ep=sst_eq.sel(time=ep_dates).mean(dim='time')
sst_cp=sst_eq.sel(time=cp_dates).mean(dim='time')
sst_nina=sst_eq.sel(time=nina_dates).mean(dim='time')
sst_neutral=sst_eq.drop_sel(time=ep_dates.append(cp_dates).append(nina_dates)).mean(dim='time')


#Remove the cp[:-7] and change -latbnd,latbnd to latbnd,-latbnd for NCEP2 rather than CCMP
ws_eq=ws.sel(lat=slice(-latbnd,latbnd)).mean(dim='lat')
ws_eq_NCEP2=ws_ncep2.sel(lat=slice(latbnd,-latbnd)).mean(dim='lat')


ws_ep=ws_eq.sel(time=ep_dates).mean(dim='time')
ws_cp=ws_eq.sel(time=cp_dates[:-7]).mean(dim='time')
ws_nina=ws_eq.sel(time=nina_dates).mean(dim='time')
ws_neutral=ws_eq.drop_sel(time=ep_dates.append(cp_dates[:-7]).append(nina_dates)).mean(dim='time')

chl_eq=chl.sel(lat=slice(-latbnd,latbnd)).mean(dim='lat')
chl_ep=chl_eq.sel(time=ep_dates).mean(dim='time')
chl_cp=chl_eq.sel(time=cp_dates[:-5]).mean(dim='time')
chl_nina=chl_eq.sel(time=nina_dates).mean(dim='time')
chl_neutral=chl_eq.drop_sel(time=ep_dates.append(cp_dates[:-5]).append(nina_dates)).mean(dim='time')


prec_eq=precip.sel(lat=slice(latbnd+0.5,-latbnd-0.5)).mean(dim='lat')
prec_ep=prec_eq.sel(time=ep_dates).mean(dim='time')
prec_cp=prec_eq.sel(time=cp_dates).mean(dim='time')
prec_nina=prec_eq.sel(time=nina_dates).mean(dim='time')
prec_neutral=prec_eq.drop_sel(time=ep_dates.append(cp_dates).append(nina_dates)).mean(dim='time')

# %% ENSO FIGURE

#ENSO BREAKDOWN
positions = (140, 160, 180,200,220,240,260)
labels = ("140$^\circ$E", "160$^\circ$E", "180$^\circ$",'160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W')


plt.figure(figsize=(10,8))

plt.subplot(4,2,1)
    
plt.plot(sst_ep.lon,sst_ep,label='East Pacific',c='darkred')
plt.plot(sst_cp.lon,sst_cp,label='Central Pacific',c='darkred',linestyle='--')
plt.plot(sst_nina.lon,sst_nina, label='La Nina',c='blue')
plt.plot(sst_ep.lon,sst_neutral,label='Neutral',c='k')
#plt.title('SST')
plt.xlim([140,265])
plt.grid()

#plt.xlabel('Mooring')
plt.ylabel('Degree C')
plt.ylim([23,31])
plt.title('a) SST',loc='left')
#print('SST')
plt.xticks(positions, labels)


ax = plt.subplot(4,2,3)
plt.plot(chl_ep.lon,chl_ep,label='East Pacific',c='darkred')
plt.plot(chl_cp.lon,chl_cp,label='Central Pacific',c='darkred',linestyle='--')
plt.plot(chl_nina.lon,chl_nina, label='La Nina',c='blue')
plt.plot(chl_ep.lon,chl_neutral,label='Neutral',c='k')
#plt.title('TPCA Chlorophyll')
plt.xlim([140,265])

plt.ylim([0.05,0.35])
plt.ylabel('mg Chl m$^{-3}$ day$^{-1}$')
plt.title('c)  TPCA Chlorophyll',loc='left')
plt.grid()
plt.xticks(positions, labels)


ax = plt.subplot(4,2,4)
plt.plot(np_ep.lon,np_ep,label='New Production East Pacific',c='darkred')
plt.plot(np_cp.lon,np_cp,label='New Production Central Pacific',c='darkred',linestyle='--')
plt.plot(np_nina.lon,np_nina, label='New Producion Nina',c='blue')
plt.plot(np_ep.lon,np_neutral,label='New Production Neutral',c='k')

plt.xlim([140,265])
plt.xticks(positions, labels)

plt.grid()

plt.ylim([xl0,yl0])
print('NEW PRODUCTION')

plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.title('d)  New production',loc='left')


#ENSO BREAKDOWN
plt.subplot(426)

plt.plot(dco2_ep.lon,dco2_ep,label='East Pacific',c='darkred')
plt.plot(dco2_cp.lon,dco2_cp,label='Central Pacific',c='darkred',linestyle='--')
plt.plot(dco2_nina.lon,dco2_nina, label='La Nina',c='blue')
plt.plot(dco2_ep.lon,dco2_neutral,label='Neutral',c='k')
plt.xlim([140,265])
plt.xticks(positions, labels)

plt.grid()
plt.ylabel('μatm')
plt.ylim([0,110])
plt.title('f) \u0394pCO$_{2}$',loc='left')
print('pCO2')


#ENSO BREAKDOWN

plt.subplot(4,2,2)
plt.plot(ws_ep.lon,ws_ep,label='East Pacific',c='darkred')
plt.plot(ws_cp.lon,ws_cp,label='Central Pacific',c='darkred',linestyle='--')
plt.plot(ws_nina.lon,ws_nina, label='La Nina',c='blue')
plt.plot(ws_ep.lon,ws_neutral,label='Neutral',c='k')
plt.xticks(positions, labels)

plt.grid()
plt.ylim([2,7])
plt.xlim([140,265])
plt.ylabel('m s$^{-1}$')
plt.title('b) Wind speed',loc='left')
print('WINDSPEED')

# Calculate and Put mooring data on just for reference

means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()

ty='month' #Actually month though need to fix this.
fp='processed/combined_dataset/month_data_exports.nc'
moorings=['110W','125W','140W','155W','170W','165E']#[::-1]

lns=[165,190,205,220,235,250][::-1]
moors=[110, 125, 140, 155, 170, 195]
for i, mooring_name in enumerate(lns):
    fp='processed/combined_dataset/month_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(moors[i]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
        
    #plt.subplot(6,1,i+1)
    #pco2_m=pco2_month.sel(lat=0,lon=mooring_name,method='nearest')

    d=dat#.sel(Date=slice(pco2_month.time.min().values,pco2_month.time.max().values))
    
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
    
    #Calculate Windspeed biases vs insitu?
    # # plt.show()    
    # # plt.plot(info.windspeed.index,info.windspeed)
    # # ws_eq.sel(lon=mooring_name,method='nearest').plot()
    # # ws_eq_NCEP2.sel(lon=mooring_name,method='nearest').plot()
    # wws_ccmp=ws_eq.sel(lon=mooring_name,method='nearest')-info.windspeed[:-8]
    # wws_ncep2=ws_eq_NCEP2.sel(lon=mooring_name,method='nearest').sel(time=slice('1997-01-01','2019-12-01'))-info.windspeed[6:]
    # wws_ccmp.plot(),wws_ncep2.plot()
    # plt.show()
   
    #pd.Serie
    wsdf ={#'El Nino':nino1.windspeed,
        'La Nina':nina1.windspeed,
 #       'Cold Modoki':nmodoki.windspeed,
        'Neutral':neutral1.windspeed,
        'CP El Nino':modoki.windspeed,
        'EP El Nino':EP.windspeed
       # 'All time':info.co2.mean()-info.select_model.mean()
        }

    ensoavgs=pd.Series(wsdf,name=moorings[i][0:3]+'\xb0'+moorings[i][3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    print(x[0])
    print(x[0]=='CP El Nino')
    ls='-'
    markr='o'
    markers=5
    if 'Neutral' in x[0]:
        c='black'
    #elif 'Nino' in x[0]:
    #    c='darkred'
    elif 'Nina' in x[0]:
        c='blue'
        markers=7
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
    elif 'CP' in x[0]:
        c='red'
        ls='--'
        markers=8
        markr='x'
            

    plt.plot(lns,x[1],c=c,marker=markr,label=x[0],linewidth=0,alpha=0.9,markersize=markers)



ax = plt.subplot(4,2,5)

plt.plot(prec_ep.lon,prec_ep,label='East Pacific',c='darkred')
plt.plot(prec_cp.lon,prec_cp,label='Central Pacific',c='darkred',linestyle='--')
plt.plot(prec_nina.lon,prec_nina, label='La Nina',c='blue')
plt.plot(prec_neutral.lon,prec_neutral,label='Neutral',c='k')
plt.xlim([140,265])
plt.grid()

plt.ylabel('mm day$^{-1}$')
plt.title('e)  Precipitation',loc='left')
plt.xticks(positions, labels)
#plt.ylim([xl0,yl0])


means=pd.DataFrame()
final_mooring_enso=pd.DataFrame()
final_mooring_enso_avgs=pd.DataFrame()

ty='month' #Actually month though need to fix this.
fp='processed/combined_dataset/month_data_exports.nc'
moorings=['110W','125W','140W','155W','170W','165E']#[::-1]

lns=[165,190,205,220,235,250][::-1]
moors=[110, 125, 140, 155, 170, 195]
for i, mooring_name in enumerate(lns):
    fp='processed/combined_dataset/month_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(moors[i]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
        
    #plt.subplot(6,1,i+1)
    #pco2_m=pco2_month.sel(lat=0,lon=mooring_name,method='nearest')

    d=dat#.sel(Date=slice(pco2_month.time.min().values,pco2_month.time.max().values))
    
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
    wsdf ={#'El Nino':nino1.windspeed,
        'La Nina':nina1.precip,
 #       'Cold Modoki':nmodoki.windspeed,
        'Neutral':neutral1.precip,
        'CP El Nino':modoki.precip,
        'EP El Nino':EP.precip
       # 'All time':info.co2.mean()-info.select_model.mean()
        }

    ensoavgs=pd.Series(wsdf,name=moorings[i][0:3]+'\xb0'+moorings[i][3:4])
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
for x in final_mooring_enso.T.iterrows():
   # if 'All time' in x[0]:
   #     c='k'
    print(x[0])
    print(x[0]=='CP El Nino')
    ls='-'
    markr='o'
    markers=5
    if 'Neutral' in x[0]:
        c='black'
    #elif 'Nino' in x[0]:
    #    c='darkred'
    elif 'Nina' in x[0]:
        c='blue'
        markers=7
    elif 'Cold Modoki' in x[0]:
        c='royalblue'
        ls='--'
    elif 'EP' in x[0]:
        c='darkred'
    elif 'CP' in x[0]:
        c='red'
        ls='--'
        markers=8
        markr='x'
            
            

    plt.plot(lns,x[1]*24,c=c,marker=markr,label=x[0],linewidth=0,alpha=0.9,markersize=markers)


ax = plt.subplot(4,2,7)

plt.plot(co2_ep.lon,co2_ep,label=r'East Pacific El Niño',c='darkred')
plt.plot(co2_cp.lon,co2_cp,label='Central Pacific El Niño',c='darkred',linestyle='--')
plt.plot(co2_nina.lon,co2_nina, label='La Niña',c='blue')
plt.plot(co2_ep.lon,co2_neutral,label='Neutral conditions',c='k')
plt.xlim([140,265])
plt.grid()

plt.title('g)  Air-sea CO$_{2}$ flux',loc='left')
plt.ylabel('gC m$^{-2}$ day$^{-1}$')
plt.xticks(positions, labels)
plt.ylim([xl0,yl0])

plt.tight_layout()
plt.legend(fontsize=12,bbox_to_anchor=(1.35,0.5), loc="center left")#, borderaxespad=0)



# %% 
#plt.savefig('figs/Figure4.png',dpi=200)
#try:
#    plt.savefig('figs/Figure4.jpeg',dpi=300)
#except:
#    pass
#plt.savefig('figs/vector/Figure4.eps')
#plt.savefig('figs/vector/Figure4.pdf')
plt.show()

# %% 
fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)

xx=dat.sel(Mooring=110)

asf=((xx.co2flux4_land_gmyr/365)).groupby('Date.month').mean()
npr=(xx.laws2011a*xx.cafe/1000).groupby('Date.month').mean()
year=(xx.groupby('Date.month').mean()-xx.groupby('Date.month').mean().mean())
 
plt.xcorr(asf,npr,maxlags=6)