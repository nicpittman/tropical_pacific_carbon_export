#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 13:38:01 2020
@author: Nic Pittman


Produces to the total mass of CO2 flux and new production for different boxes.

Requires:
    'datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc'
    'processed/flux/fratios.nc'
    avg_npp_rg_cbpm.nc
    datasets/npp_satellite/avg_npp_rg_cafe.nc
    xr.open_dataarray('processed/earth_m2.nc
    
Produces:
    processed/results/enso_basin_means.csv
    processed/results/carbon_mass.csv
    figs/Figure6_basinavg_pG.png
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
from matplotlib.ticker import ScalarFormatter
import matplotlib
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of 
    magnitude
    https://stackoverflow.com/a/3679918/9965678
    """
    def __init__(self, order_of_mag=0, useOffset=True, useMathText=False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset=useOffset, 
                                 useMathText=useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag

def trends(ax,x,y,c='r'):  
    from scipy.stats import linregress
    mean=np.nanmean(y)
    std=np.nanstd(y)*1

    #x_n=np.arange(0,len(x))
    # x1=np.arange(np.datetime64(x[0],'M'),np.datetime64(x[-1],'M')+np.timedelta64(1,'M'))
    #x1=trd.index.values.astype('datetime64[D]')
    x1=x.astype('datetime64[D]')
    x_n=pd.to_numeric(x1)
    slope, intercept, r_value, p_value, std_err = linregress(x_n,y)
    mn=min(x_n)
    mx=max(x_n)
    x1=np.linspace(mn,mx,len(x))
    y1=slope*x1+intercept
    
    ax.plot(pd.to_datetime(x),y1,':'+c,linewidth=2.5)  
    #ax.text(x1[-1]-(x1[-1]*0.1),y1[-1]-(y1[-1]*0.1),'R2='+str(np.round(r_value**2,3)))
    return slope, intercept, r_value,p_value,std_err
    

seamask=xr.open_dataset('processed/seamask.nc') #Because 2020 version doesn't have it.
seamask= seamask.assign_coords(lon=(seamask.lon % 360)).roll(lon=(seamask.dims['lon']),roll_coords=False).sortby('lon')	

landsch_fp='datasets/co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc'
landschutzer=xr.open_dataset(landsch_fp)
landschutzer= landschutzer.assign_coords(lon=(landschutzer.lon % 360)).roll(lon=(landschutzer.dims['lon']),roll_coords=False).sortby('lon')		#EPIC 1 line fix for the dateline problem.
land_pac=landschutzer.sel(lon=slice(120,290),lat=slice(-20,20))
#land_pac.to_netcdf('processed/fluxmaps/landshutzer.nc')
land_pac=moles_to_carbon(land_pac.fgco2_smoothed)

#JMA=moles_to_carbon(xr.open_mfdataset('datasets/co2/JMA_co2/jma_flux.nc').flux.sel(lon=slice(120,290),lat=slice(-20,20)))
#yasanaka=moles_to_carbon(xr.open_mfdataset('datasets/co2/yasanaka_co2/Yasunaka_pCO2_flux.nc').flux_masked).sel(lon=slice(120,290),lat=slice(-20,20))
 
f_ratios=xr.open_mfdataset('processed/flux/fratios.nc')
ratio=f_ratios.laws2000#laws2000,laws2011a,laws2011b,henson2011

npp=(xr.open_dataset('processed/flux/avg_npp_rg_cbpm.nc').avg_npp/1000*365)
npp=(xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc').avg_npp/1000*365)

#grid=xr.open_dataarray('processed/tropics_size_m2.nc')
grid=xr.open_dataarray('processed/earth_m2.nc')
grid['lon']=grid.lon+180
grid=grid.where(seamask.seamask==1)
#fig=plt.figure(figsize=(8,10))

limits=[['West',165,180],
      ['Central',190,205],
      ['East',230,245],
      #['Basin',180,280]]
      ['Basin',150,280]]
      #['Basin',135,280]] #Le Borgne 2002 Warm Pool
      
lims=10

mass_table=pd.DataFrame({
    'Region':[],
    'Area':[],
    'NP Mean':[],
    'NP Neutral':[],
    'NP El Nino':[],
    'NP La Nina':[],
    'NP El Nino Diff':[],
    'NP La Nina Diff':[],
    'NP Trends':[],
    'NP Pval':[],
    
    'CO2 Mean':[],
    'CO2 Neutral':[],
    'CO2 El Nino':[],
    'CO2 La Nina':[],
    'CO2 El Nino Diff':[],
    'CO2 La Nina Diff':[],
    'CO2 Trends':[],
    'CO2 Pval':[]})

#Calculating overall flux rates.
plt.figure(figsize=(13,10))
for i,ty in enumerate(limits):
    if i <=2:
        ax = plt.subplot(2,3,i+1)
    else:
        ax=plt.subplot(2,1,2)
    
    startl= ty[1]#120#80#160#135#180#135#180 #not 150
    endl=ty[2]#280 #80W #270 #90W
    gs=grid.sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum()
    print('\n'+ty[0])
    print('gridsize= '+str(gs.values/1e13))    
    
    # if i==0:
    #     #ax.text()
    #     plt.text('2000-01-01',1.6e14,'10NS, 165E-180W, '+str(np.round(gs.values/1e13,3))+'x10$^{13}$ m$^2$')
    # elif i==1:
    #     plt.text('2000-01-01',1.6e14,'10NS, 170W-155W, '+str(np.round(gs.values/1e13,3))+'x10$^{13}$ m$^2$')
    # elif i==2:
    #     plt.text('2000-01-01',1.6e14,'10NS, 130W-115W, '+str(np.round(gs.values/1e13,3))+'x10$^{13}$ m$^2$')
    # elif i==3:
    #     plt.text('2007-01-01',trenNP[3]0.5e14,'10NS, 150E-80W, '+str(np.round(gs.values/1e13,3))+'x10$^{13}$ m$^2$')
    

    #Year average of CO2 flux
    CO2=(land_pac*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    #JMC=(JMA*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    #YAS=(yasanaka*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    #plt.show()
    CO2['time']=CO2['time'].astype('datetime64[M]')
    #Year average of CO2 flux
    henson=(npp*f_ratios.henson2011*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    laws2000=(npp*f_ratios.laws2000*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    laws2011a=(npp*f_ratios.laws2011a*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    laws2011b=(npp*f_ratios.laws2011b*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])
    
    dunne=(npp*f_ratios.dunne2005*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])#.plot(label='Dunne 2005')
    trim=(npp*f_ratios.trim*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])#.plot(label='Dunne 2005')
    
    #JMC.plot(label='Iida flux',ax=ax)
    #YAS.plot(label='Yasanaka flux',ax=ax)
    #ax.plot(CO2.time,CO2,c='k')
    CO2.plot(label='CO2 outgassing',ax=ax,c='k')
    
    CO=CO2.sel(time=slice('2000-01-01','2017-12-01'))
    trenCO=trends(ax,CO.time.values,CO.values,c='k')
    annual_rate_of_changeCO=trenCO[0]*365
 
    #henson.plot(label='Henson',ax=ax)
    laws2011a.plot(label='Laws2011a',ax=ax,c='r')
    #laws2011b.plot(label='Laws2011b',ax=ax,c='pink')
            
    mod=laws2011a.sel(time=slice('2000-01-01','2017-12-01'))
    trenNP=trends(ax,mod.time.values,mod.values)
    annual_rate_of_changeNP=trenNP[0]*365
   
   
    #(CO2+laws2011a).plot(label='combined',ax=ax,c='m')
    if i==3:
        trim.plot(label='DeVries and Webber 2017',ax=ax,c='purple')
    #laws2011b.plot(label='Laws2011b',ax=ax,c='slategray',linewidth=2)
        dunne.sel(time=slice('1997-01-01','2019-07-01')).plot(label='Dunne 2005',ax=ax,c='steelblue')
        laws2000.plot(label='Laws2000',ax=ax,c='darkblue')
        #henson.plot(label='Henson',ax=ax,c='deeppink')
    
    
    
    ax.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
    
    ax.set_xlabel('Year')
    
    if i <=2:
        #ax.set_ylim([0,0.27*1e15])131
        
        ax.set_ylim([-0.005*1e15,0.175*1e15])
        ax.set_title(chr(97+i)+') '+ty[0]+' Pacific carbon export',pad=16)
        ax.set_ylabel('Carbon export / outgassing (PgC/yr$^{-1}$)')
        ax.yaxis.set_major_formatter(FixedOrderFormatter(15))

       #y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        #ax.yaxis.set_major_formatter(y_formatter)
    else:
        ax.legend(loc='upper center',ncol=5)
        ax.set_title(chr(97+i)+') Tropical Pacific carbon export')
        ax.set_ylabel('Carbon export / outgassing (PgC/yr$^{-1}$)')
        ax.set_ylim([0,1.25*1e15])
    import matplotlib.patches as patches
      

    dat=CO2.to_dataframe(name='CO2')
    dat['henson']=henson.to_dataframe(name='henson').henson
    dat['laws2011a']=laws2011a.to_dataframe(name='laws2011a').laws2011a
    dat['laws2011b']=laws2011b.to_dataframe(name='laws2011b').laws2011b
    
    dat['laws2000']=laws2000.to_dataframe(name='laws2000').laws2000
    dat['dunne']=dunne.to_dataframe(name='dunne2005').dunne2005
   # dat['JMC']=JMC.to_dataframe(name='jmc').jmc
    
    #gs=grid.sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum()
    #print('gridsize= '+str(gs.values/1e14))
    
    print('CO2: '+str(dat.CO2.mean()/1e15))
    print('henson: '+str(dat.henson.mean()/1e15))
    print('laws2011a: '+str(dat.laws2011a.mean()/1e15))
    print('laws2011b: '+str(dat.laws2011b.mean()/1e15))
    print('laws2000: '+str(dat.laws2000.mean()/1e15))
    print('Dunne: '+str(dat.dunne.mean()/1e15))
    #print('JMC: '+str(dat.JMC.mean()/1e15))
    
    print('NP Trend: '+ str(annual_rate_of_changeNP/1e15)+' ' +u"\u00B1 "+str((trenNP[4]*365)/1e15) +' PgC/yr/yr')
    print('pval= '+str(trenNP[3]))
    print('CO2 Trend: '+str(annual_rate_of_changeCO/1e15)+' ' +u"\u00B1 " +str((trenCO[4]*365)/1e15)+' PgC/yr/yr')
    print('pval= '+str(trenCO[3]))
    
    dat=dat[dat.index>np.datetime64('1997-08')]
    nino=pd.DataFrame()
    nina=pd.DataFrame()
    ensofps=['processed/indexes/el_nino_events.csv','processed/indexes/la_nina_events.csv']
    for whichenso,fp in enumerate(ensofps):
        events=pd.read_csv(fp)
        for ev in events.iterrows():
            endm=np.datetime64(ev[1].end).astype('datetime64[M]')
            endm1=endm-np.timedelta64(1,'M')
            start=np.datetime64(ev[1].start).astype('datetime64[M]')
            if np.datetime64(start)==np.datetime64(endm1):
                pass
                #print('fail',start,endm1,whichenso)
            else:
                #print(start,endm,whichenso)
                if whichenso==0:
                    #if el nino
                    
                    vol=dat[(dat.index>=start) & (dat.index<=endm)]
                    nino=nino.append(vol)
                    patchcol='firebrick'
                elif whichenso==1:
                    patchcol='deepskyblue'
                    #if la nina
                    vol=dat[(dat.index>=start) & (dat.index<=endm)]
                    nina=nina.append(vol)

                rect=patches.Rectangle((start,-0.005*1e15),endm-start,1e16,linewidth=0,alpha=0.2,color=patchcol)
                ax.add_patch(rect)
                rect=patches.Rectangle((start,-0.005*1e15),endm-start,1e16,linewidth=0,alpha=0.2,color=patchcol)
                ax.add_patch(rect)
    neutral=dat.drop(index=nina.index).drop(index=nino.index)
    
    
    
    
    # #Put in the ENSO box regions
    # ensofps=['processed/indexes/el_nino_events.csv','processed/indexes/la_nina_events.csv']
    # for whichenso,fp in enumerate(ensofps):
    #     events=pd.read_csv(fp)
    #     for ev in events.iterrows():
    #         endm=np.datetime64(ev[1].end).astype('datetime64[M]')
    #         endm1=endm-np.timedelta64(1,'M')
    #         start=np.datetime64(ev[1].start).astype('datetime64[M]')
    #         if start==endm1:
    #             pass    
    #         else:
    #             if whichenso==0:
    #                 #if el nino
    #                 patchcol='firebrick'
    #             else:
    #                 #if la nina
    #                 patchcol='deepskyblue'
    #             rect=patches.Rectangle((start,-0.005*1e15),endm-start,1e16,linewidth=0,alpha=0.2,color=patchcol)
    #             ax.add_patch(rect)
    #             rect=patches.Rectangle((start,-0.005*1e15),endm-start,1e16,linewidth=0,alpha=0.2,color=patchcol)
    #             ax.add_patch(rect)
    
    
    na=pd.DataFrame(nino.mean(axis=0)/1e15).T
    nb=pd.DataFrame(nina.mean(axis=0)/1e15).T
    nc=pd.DataFrame(neutral.mean(axis=0)/1e15).T
    nd=pd.DataFrame(dat.mean(axis=0)/1e15).T
    
    na1=pd.DataFrame(nino.std(axis=0)/1e15).T
    nb1=pd.DataFrame(nina.std(axis=0)/1e15).T
    nc1=pd.DataFrame(neutral.std(axis=0)/1e15).T
    nd1=pd.DataFrame(dat.std(axis=0)/1e15).T
    
    na['name']='nino mean'
    na1['name']='nino std'
    nb['name']='nina mean'
    nb1['name']='nina std'
    nc['name']='neutral mean'
    nc1['name']='neutral std'
    nd['name']='all mean'
    nd1['name']='all std'
    
    #print(na,nb,nc,nd)
    means=na.append(na1).append(nb).append(nb1).append(nc).append(nc1).append(nd).append(nd1)
    means=means.round(3)
    print(means)
    
    
    
    
    def perc(a,b):
        #New Number - Original Number / orinal
        c=(a-b)/b
        d=a-b
        return (c*100).round(1),d.round(3)
        
    p0=perc(means.laws2011a.iloc[2],means.laws2011a.iloc[4]) #Nina
    p1=(perc(means.laws2011a.iloc[0],means.laws2011a.iloc[4])) #Nino
    p2=(perc(means.CO2.iloc[2],means.CO2.iloc[4])) #Nina
    p3=(perc(means.CO2.iloc[0],means.CO2.iloc[4])) #Nina
    
    print('Nina change% NP: '+str(p0))
    print('Nino change% NP: '+str(p1))
    print('Nina change% CO: '+str(p2))
    print('Nino change% CO: '+str(p3))
          
    mass_table=mass_table.append(pd.DataFrame({
    'Region':ty[0],
    'Area':str((gs.values/1e12).round(3))+' (10^12m)',
    'NP Mean':str(nd.laws2011a.round(3).values[0])+' PgC/yr',
    'NP Neutral':str(nc.laws2011a.round(3).values[0])+' PgC/yr',
    'NP El Nino':str(na.laws2011a.round(3).values[0])+' PgC/yr',
    'NP La Nina':str(nb.laws2011a.round(3).values[0])+' PgC/yr',
    'NP El Nino Diff':str(p0[0])+'% ,'+str(p0[1])+'PgC/yr',
    'NP La Nina Diff':str(p1[0])+'% ,'+str(p1[1])+'PgC/yr',
    'NP Trends':str((annual_rate_of_changeNP/1e12).round(3))+' ' +u"\u00B1 "+str(((trenNP[4]*365)/1e12).round(3)) +' TgC/yr-2',
    'NP Pval':trenNP[3].round(7),
    'CO2 Mean':str(nd.CO2.round(3).values[0])+' PgC/yr',
    'CO2 Neutral':str(nc.CO2.round(3).values[0])+' PgC/yr',
    'CO2 El Nino':str(na.CO2.round(3).values[0])+' PgC/yr',
    'CO2 La Nina':str(nb.CO2.round(3).values[0])+' PgC/yr',
    'CO2 El Nino Diff':str(p2[0])+'% ,'+str(p2[1])+'PgC/yr',
    'CO2 La Nina Diff':str(p3[0])+'% ,'+str(p3[1])+'PgC/yr',
    'CO2 Trends':str((annual_rate_of_changeCO/1e12).round(3))+' ' +u"\u00B1 "+str(((trenNP[4]*365)/1e12).round(3)) +' TgC/yr-2',
    'CO2 Pval':trenCO[3].round(7)},index=[i]))
        
    
    
    
    
means.to_csv('processed/results/enso_basin_means.csv',index=False)
mass_table=mass_table.T
mass_table.to_csv('processed/results/carbon_mass.csv')
plt.tight_layout()
plt.savefig('figs/Figure6_basinavg_pG.png',dpi=200)
plt.show()

print(gs)
