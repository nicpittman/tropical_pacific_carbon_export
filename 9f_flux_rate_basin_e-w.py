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
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, ScalarFormatter)
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
      ['Central',205,220],
      ['East',235,250],
      #['Basin_check',180,280]]
      ['Basin',150,270]] #As used in the paper

      #['Ishii',135,270]]
      #['Borgne',135,270]] #Le Borgne 2002 Warm Pool
      #['Wyrtiki',180,270]]
#Borgne is lims=1
#Wyrtiki is both lims=5 and 10
lims=15


mass_table=pd.DataFrame({})



check_lag_corr_x=[]
check_lag_corr_y=[]

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
    
    #dunne=(npp*f_ratios.dunne2005_tpca*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])#.plot(label='Dunne 2005')
    dunne=(npp*f_ratios.dunne2005*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])#.plot(label='Dunne 2005')
 
    trim=(npp*f_ratios.trim*grid).sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum(dim=['lat','lon'])#.plot(label='Dunne 2005')
    
    #JMC.plot(label='Iida flux',ax=ax)
    #YAS.plot(label='Yasanaka flux',ax=ax)
    #ax.plot(CO2.time,CO2,c='k')
   
    
    CO=CO2.sel(time=slice('2000-01-01','2019-12-01'))
    trenCO=trends(ax,CO.time.values,CO.values,c='k')
    annual_rate_of_changeCO=trenCO[0]*365
 
    #henson.plot(label='Henson',ax=ax)
    
    #laws2011b.plot(label='Laws2011b',ax=ax,c='pink')
            
    mod=laws2011a.sel(time=slice('2000-01-01','2019-12-01'))
    trenNP=trends(ax,mod.time.values,mod.values)
    annual_rate_of_changeNP=trenNP[0]*365
   
   
    #(CO2+laws2011a).plot(label='combined',ax=ax,c='m')
    if i==3:
        trim.plot(label='DeVries and Webber 2017',ax=ax,c='darkorange',linewidth=2)#,linestyle='--')
    #laws2011b.plot(label='Laws2011b',ax=ax,c='slategray',linewidth=2)
        dunne.sel(time=slice('1997-01-01','2019-07-01')).plot(label='Dunne 2005',ax=ax,c='darkblue')
        laws2000.plot(label='Laws 2000',ax=ax,c='green',linestyle='--')
        laws2011a.plot(label='Laws 2011a',ax=ax,c='r',linewidth=2.5)
        CO2.plot(label='CO$_{2}$ outgassing',ax=ax,c='k',linewidth=2.5)
    else:
        laws2011a.plot(label='Laws 2011a',ax=ax,c='r')
        CO2.plot(label='CO$_{2}$ outgassing',ax=ax,c='k')
        #henson.plot(label='Henson',ax=ax,c='deeppink')
    
    
    
    ax.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
    #ax.xaxis.grid(True, which='both')
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    #plt.grid()
    ax.set_xlabel('Year')
    
    if i <=2:
        #ax.set_ylim([0,0.27*1e15])131
        
        ax.set_ylim([-0.005*1e15,0.25*1e15])
        ax.set_title(chr(97+i)+') '+ty[0]+' Pacific',pad=16)
        ax.set_ylabel('New production and CO$_{2}$ flux (PgC yr$^{-1}$)')
        ax.yaxis.set_major_formatter(FixedOrderFormatter(15))

       #y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        #ax.yaxis.set_major_formatter(y_formatter)
    else:
        ax.legend(loc='lower center',ncol=5)
        ax.set_title(chr(97+i)+') Entire Basin')
        ax.set_ylabel('New Production and CO$_{2}$ flux (PgC yr$^{-1}$)')
        ax.set_ylim([0,1.6*1e15])
    import matplotlib.patches as patches
      

    dat=CO2.to_dataframe(name='CO2')
    dat['henson']=henson.to_dataframe(name='henson').henson
    dat['laws2011a']=laws2011a.to_dataframe(name='laws2011a').laws2011a
    dat['laws2011b']=laws2011b.to_dataframe(name='laws2011b').laws2011b
    
    dat['laws2000']=laws2000.to_dataframe(name='laws2000').laws2000
    dat['dunne']=dunne.to_dataframe(name='dunne2005').dunne2005
    dat['strim']=trim.to_dataframe(name='simpletrim').simpletrim
    
   # dat['JMC']=JMC.to_dataframe(name='jmc').jmc
    
    #gs=grid.sel(lat=slice(-lims,lims)).sel(lon=slice(startl,endl)).sum()
    #print('gridsize= '+str(gs.values/1e14))
    
    print('CO2: '+str(dat.CO2.mean()/1e15))
    print('CO2 STD: '+str(dat.CO2.std()/1e15))
    print('henson: '+str(dat.henson.mean()/1e15))
    print('laws2011a: '+str(dat.laws2011a.mean()/1e15))
    print('laws2011a STD: '+str(dat.laws2011a.std()/1e15))
    
    print('laws2011b: '+str(dat.laws2011b.mean()/1e15))
    print('laws2000: '+str(dat.laws2000.mean()/1e15))
    print('Dunne: '+str(dat.dunne.mean()/1e15))
    #print('JMC: '+str(dat.JMC.mean()/1e15))
    
    print('NP Trend: '+ str(annual_rate_of_changeNP/1e15)+' ' +u"\u00B1 "+str((trenNP[4]*365)/1e15) +' PgC/yr/yr')
    print('pval= '+str(trenNP[3]))
    print('CO2 Trend: '+str(annual_rate_of_changeCO/1e15)+' ' +u"\u00B1 " +str((trenCO[4]*365)/1e15)+' PgC/yr/yr')
    print('pval= '+str(trenCO[3]))
    
    dat=dat[dat.index>np.datetime64('1997-08')]
    
    
    lanina=pd.read_csv('processed/indexes/la_nina_events.csv')
    cp_nino=pd.read_csv('processed/indexes/cp_events.csv')
    #cpc.to_csv('processed/indexes/cold_cp_events.csv')
    ep_nino=pd.read_csv('processed/indexes/ep_events.csv')

    info=dat
    ninaf=pd.DataFrame()
    epf=pd.DataFrame()
    cpf=pd.DataFrame()
    for i in lanina.iterrows(): ninaf=ninaf.append(info[slice(i[1].start,i[1].end)])
    for i in ep_nino.iterrows(): epf=epf.append(info[slice(i[1].start,i[1].end)])
    for i in cp_nino.iterrows(): cpf=cpf.append(info[slice(i[1].start,i[1].end)])
    nina_dates=ninaf.index
    ep_dates=epf.index
    cp_dates=cpf.index
    
    ensofps=['processed/indexes/ep_events.csv','processed/indexes/la_nina_events.csv','processed/indexes/cp_events.csv']
    for whichenso,fp in enumerate(ensofps):
        events=pd.read_csv(fp)
        for ev in events.iterrows():
            endm=np.datetime64(ev[1].end).astype('datetime64[M]')
            endm1=endm-np.timedelta64(1,'M')
            endm2=endm+np.timedelta64(1,'M')
            start=np.datetime64(ev[1].start).astype('datetime64[M]')
          
            if start==endm1: #We don't want to plot events that last for only a month
                pass    
            #elif start==endm2-np.timedelta64(1,'M'): #There was some weirdness with the 2015 event not being continuous, and this fixes it..,
            #    pass
            else:
                if whichenso==0:
                    #if el nino
                    patchcol='darkred'#'firebrick'
                elif whichenso==1:
                    #if la nina
                    patchcol='deepskyblue'
                elif whichenso==2:
                    patchcol='darkorange'
                rect=patches.Rectangle((start,-1*10e15),endm-start,2*10e15,linewidth=0,alpha=0.3,color=patchcol)
                ax.add_patch(rect)
                rect=patches.Rectangle((start,-1*10e15),endm-start,2*10e15,linewidth=0,alpha=0.3,color=patchcol)
                ax.add_patch(rect)
       
    
   
    epnino=info[info.index.isin(ep_dates)]
    cpnino=info[info.index.isin(cp_dates)]
    nina=info[info.index.isin(nina_dates)]
    neutral1=info[~info.index.isin(cp_dates)]
    neutral1=neutral1[~neutral1.index.isin(ep_dates)]
    neutral=neutral1[~neutral1.index.isin(nina_dates)]

    
    
    na=pd.DataFrame(epnino.mean(axis=0)/1e15).T
    nb=pd.DataFrame(nina.mean(axis=0)/1e15).T
    nc=pd.DataFrame(neutral.mean(axis=0)/1e15).T
    nd=pd.DataFrame(dat.mean(axis=0)/1e15).T
    ne=pd.DataFrame(cpnino.mean(axis=0)/1e15).T
    
    
    na1=pd.DataFrame(epnino.std(axis=0)/1e15).T
    nb1=pd.DataFrame(nina.std(axis=0)/1e15).T
    nc1=pd.DataFrame(neutral.std(axis=0)/1e15).T
    nd1=pd.DataFrame(dat.std(axis=0)/1e15).T
    ne1=pd.DataFrame(cpnino.std(axis=0)/1e15).T
    
    na['name']='ep mean'
    na1['name']='ep std'
    nb['name']='nina mean'
    nb1['name']='nina std'
    nc['name']='neutral mean'
    nc1['name']='neutral std'
    nd['name']='all mean'
    nd1['name']='all std'
    ne['name']='cp mean'
    ne1['name']='cp std'
    
    
    #print(na,nb,nc,nd)
    means=na.append(na1).append(nb).append(nb1).append(nc).append(nc1).append(nd).append(nd1).append(ne).append(ne1)
    means=means.round(3)
    print(means)
    
    
    
    
    def perc(a,b):
        #New Number - Original Number / orinal
        c=(a-b)/b
        d=a-b
        return (c*100).round(1),d.round(3)
        
    p0=(perc(means.laws2011a.iloc[2],means.laws2011a.iloc[4])) #Nina
    p1=(perc(means.laws2011a.iloc[0],means.laws2011a.iloc[4])) #EP Nino
    p1a=(perc(means.laws2011a.iloc[8],means.laws2011a.iloc[4])) #CP Nino
    
    p2=(perc(means.CO2.iloc[2],means.CO2.iloc[4])) #Nina
    p3=(perc(means.CO2.iloc[0],means.CO2.iloc[4])) #EP Nino
    p3a=(perc(means.CO2.iloc[8],means.CO2.iloc[4])) #CP Nino
    
     
    #ADD CP here
    
    print('Nina change% NP: '+str(p0))
    print('Nino change% NP: '+str(p1))
    print('Nina change% CO: '+str(p2))
    print('Nino change% CO: '+str(p3))
    print('NP P VAL: '+ str(trenNP[3]))
    mass_table=mass_table.append(pd.DataFrame({
    'Region':ty[0],
    'Area':str((gs.values/1e12).round(3))+' (10^12m)',
    'NP Mean (PgC yr-1)':str(nd.laws2011a.round(3).values[0]),
    'NP Neutral (PgC yr-1)':str(nc.laws2011a.round(3).values[0]),
    'NP EP El Nino (PgC yr-1)':str(na.laws2011a.round(3).values[0]),
    'NP CP El Nino (PgC yr-1)':str(ne.laws2011a.round(3).values[0]),
    'NP La Nina (PgC yr-1)':str(nb.laws2011a.round(3).values[0]),
    'NP EP El Nino Diff':str(p1[0])+'% ,'+str(p1[1]),
    'NP CP El Nino Diff':str(p1a[0])+'% ,'+str(p1a[1]),
    'NP La Nina Diff':str(p0[0])+'% ,'+str(p0[1]),
    'NP Trends (TgC yr-2)':str((annual_rate_of_changeNP/1e12).round(3))+' ' +u"\u00B1 "+str(((trenNP[4]*365)/1e12).round(3)),
    'NP Pval':trenNP[3].round(15),
    
    'CO2 Mean (PgC yr-1)':str(nd.CO2.round(3).values[0]),
    'CO2 Neutral (PgC yr-1)':str(nc.CO2.round(3).values[0]),
    'CO2 EP El Nino (PgC yr-1)':str(na.CO2.round(3).values[0]),
    'CO2 CP El Nino (PgC yr-1)':str(ne.CO2.round(3).values[0]),
    'CO2 La Nina (PgC yr-1)':str(nb.CO2.round(3).values[0]),
    'CO2 EP El Nino Diff':str(p3[0])+'% ,'+str(p3[1]),
    'CO2 CP El Nino Diff':str(p3a[0])+'% ,'+str(p3a[1]),
    'CO2 La Nina Diff':str(p2[0])+'% ,'+str(p2[1]),
    'CO2 Trends (TgC yr-2)':str((annual_rate_of_changeCO/1e12).round(3))+' ' +u"\u00B1 "+str(((trenNP[4]*365)/1e12).round(3)),
    'CO2 Pval':trenCO[3].round(10)},index=[i]))
        
    check_lag_corr_x.append(laws2011a.sel(time=slice('1998-01-01','2019-12-01')).values)
    check_lag_corr_y.append(CO2.sel(time=slice('1998-01-01','2019-12-01')).values)
    
    
means.to_csv('processed/results/enso_basin_means.csv',index=False)
mass_table=mass_table.T
mass_table.to_csv('processed/results/carbon_mass.csv',header=False)
plt.tight_layout()
plt.savefig('figs/Figure6.png',dpi=200)
plt.savefig('figs/vector/Figure6.eps',dpi=300)
plt.savefig('figs/vector/Figure6.pdf',dpi=300)

try:
    plt.savefig('figs/Figure6.jpeg',dpi=300) #Conda install pilliow needed to save to jpeg.
except:
    pass

        
plt.show()
print(means[['name','laws2011a','dunne','strim']].iloc[6])

m=means[['name','laws2011a','dunne','strim']].iloc[6].values[1:].mean()
st=means[['name','laws2011a','dunne','strim']].iloc[6].values[1:].std()
print('basin Mean of L2011a, Dunne and strim: '+str(m)+'+-'+str(st))
#Check Lag coefficient.

for i in range(len(check_lag_corr_x)):
    hh=plt.xcorr(check_lag_corr_x[i],check_lag_corr_y[i],maxlags=3,normed=True)
    coefs=hh[1]
    index_max = max(range(len(coefs)), key=coefs.__getitem__)
    print(hh[0][index_max],hh[1][index_max])
    print(hh[1][3])
    print(hh[1][index_max]-hh[1][3])
print(gs)
