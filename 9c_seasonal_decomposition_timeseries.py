#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:52:00 2020
@author: Nic Pittman

This script reproduces Pittman et al., 2021 Figure 3.
It also produces the same for just new production and just air-sea flux.


Requires:
    processed/flux/landsch_mooring_co2_flux.nc
    processed/flux/npp.nc
    processed/combined_dataset/month_data_exports.nc
    
produces:
    figs/Figure3_decomp_megaplot'+div_factor.name+typ+'.png',dpi=200) (Three different figures)
    All trends are inside the figure.
"""


import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
import matplotlib.patches as patches
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.seasonal import STL
from scipy.stats import linregress
from matplotlib.ticker import AutoMinorLocator


#Not recommended but allows us to retrieve the numbers we want.

import warnings
warnings.filterwarnings("ignore")

moorings=['110W','125W','140W','155W','170W','165E'][::-1]


fs=16
lw=2.5
ms=0
a=0.8


def trends(ax,x,y):  

    #x_n=np.arange(0,len(x))
    # x1=np.arange(np.datetime64(x[0],'M'),np.datetime64(x[-1],'M')+np.timedelta64(1,'M'))
    #x1=trd.index.values.astype('datetime64[D]')
    x1=x.values.astype('datetime64[D]')
    x_n=pd.to_numeric(x1)
    slope, intercept, r_value, p_value, std_err = linregress(x_n,y)
    mn=min(x_n)
    mx=max(x_n)
    x1=np.linspace(mn,mx,len(x))
    y1=slope*x1+intercept
    
    ax.plot(pd.to_datetime(x),y1,':r',linewidth=2.5)  
    #ax.text(x1[-1]-(x1[-1]*0.1),y1[-1]-(y1[-1]*0.1),'R2='+str(np.round(r_value**2,3)))
    return slope, intercept, r_value,p_value,std_err
    

t='magnitude'
#t='percent'
sday='1999-01-01'
sday='1997-09-01'
seasonaltrend=[]

tsd = '2000-01-01' #Trend start at


tyyp=['np','asf','both'] #So we can do all three plot types. Onl one in the paper is 'both'
for typ in tyyp:

#190 mm x 230 mm
    print('\n\n'+typ)
    fig=plt.figure(figsize=((19/2.54)*2,(23/2.54)*2))#,constrained_layout=True)
    s=fig.add_gridspec(3,2,hspace=0.15,wspace=0.25)

    for i, mooring_name in enumerate(moorings):
        xi=0
        yi=0
        if (i+1)%2!=0:
            xi=0
        else:
            xi=1
        row_number=int(np.floor(i/2))    
        yi=row_number
        ax = s[yi,xi].subgridspec(4,1)#.add_gridspec(spec[yi,xi])
        
        ax1=fig.add_subplot(ax[0])
        ax2=fig.add_subplot(ax[1])
        ax3=fig.add_subplot(ax[2])
        ax4=fig.add_subplot(ax[3])
        
        print(mooring_name)
        #JMAco2flux_data=xr.open_mfdataset('processed/flux/JMA_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
        LandSch_co2flux_data=xr.open_mfdataset('processed/flux/landsch_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
        npp=xr.open_mfdataset('processed/flux/npp.nc').sel(Mooring=mooring_name).sel(Date=slice(sday,'2017-12-15'))
        
        ty='month' #Actually month though need to fix this.
        fp='processed/combined_dataset/'+ty+'_data_exports.nc'
        try:
            dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
        except:
            dat=xr.open_mfdataset(fp).sel(Mooring=195)
        dat['Date']=dat.Date.astype('datetime64[M]')
        #dat=dat.sel(Date=slice('1997-09-01','2017-12-15'))
        dat=dat.sel(Date=slice(sday,'2017-12-15'))
    
        div_factor=dat.laws2011a#dat.laws2011b#f_ratio#laws2011b_vgpm
    
        title=str(mooring_name)#+' (CO2 flux - New Production)'
     
    
        epp=pd.DataFrame([npp.mod_eppley.values,npp.sw_eppley.values]).mean()/1000*div_factor
        cbpm=pd.DataFrame([npp.mod_cbpm.values,npp.sw_cbpm.values,npp.viirs_cbpm.values]).mean()/1000*div_factor
        vgpm=pd.DataFrame([npp.sw_vgpm.values,npp.mod_vgpm.values,npp.viirs_vgpm.values]).mean()/1000*div_factor
        cafe=pd.Series([npp.mod_cafe.values,npp.sw_cafe.values]).mean()/1000*div_factor
        
        #avg_npp=npp[['viirs_eppley','viirs_cbpm','viirs_vgpm','sw_eppley','sw_cbpm','sw_vgpm','mod_cafe','mod_eppley','mod_cbpm','mod_vgpm']].to_dataframe().drop(columns='Mooring').mean(axis=1).to_xarray()
        avg_npp=npp[['viirs_cbpm','sw_cbpm','mod_cbpm']].to_dataframe().drop(columns='Mooring').mean(axis=1).to_xarray()
        avg_npp=npp[['sw_cafe','mod_cafe']].to_dataframe().drop(columns='Mooring').mean(axis=1).to_xarray()
        
        
        npa=avg_npp/1000*div_factor
        
        LandSch_co2flux_data['Date']=LandSch_co2flux_data.Date.astype('datetime64[M]')
    
        land_flux=((moles_to_carbon(LandSch_co2flux_data.sel(Date=slice(sday,'2020-01-01'))).fgco2_smoothed.values)/365)
        
        #JMA= moles_to_carbon(JMAco2flux_data.flux/365)
        
       #((x-y)/(y))*100 
        if t=='magnitude':
            
            tyyp=['np','asf','both']
            if typ=='np':
                datset=pd.DataFrame({'nppavg':(npa)})
            elif typ=='asf':
                datset=pd.DataFrame({'nppavg':(land_flux)})
            elif typ=='both':               
                datset=pd.DataFrame({'nppavg':(land_flux-npa)}) #})## JAP
            #datset=pd.DataFrame({'nppavg':(land_flux-npa)}) #})# #LandSchutz
        #elif t=='percent':
        #    datset=pd.DataFrame({'nppavg':((land_flux-npa)/npa)*100}) #})#
        
        datset=datset.set_index(LandSch_co2flux_data.sel(Date=slice(sday,'2020-01-01')).Date.values)
        decomp=seasonal_decompose(datset, model='addative', extrapolate_trend='freq')
        decomp=STL(datset,seasonal=13).fit()
        dates=decomp.resid.index.values.astype('datetime64')
     
        units='\n gC m$^{-2}$ day$^{-1}$'
        col='k'
        ax1.axhline(0,c='gray')#,linestyle=':')
        ax1.plot(decomp.observed,c=col)
        ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax1.set_ylabel('Observations'+units,fontsize=10)
        ax1.grid(axis='y',which='major')
        ax1.grid(axis='x',which='both')
        ax1.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
        #ax1.set_ylim([-0.12,0.12])
        if typ=='both':
            ax1.set_ylim([-0.2,0.12])
     
        trd=decomp.observed[decomp.observed.index>=tsd]
        diff110=trd
       # tren=trends(ax2,trd.index,trd.values)
        tren=trends(ax1,trd.index,trd.values.squeeze())
        
        #annual_rate_of_change1=((tren[0]*pd.to_numeric(trd.index)[-1]+tren[1])-(tren[0]*pd.to_numeric(trd.index)[0]+tren[1]))/len(trd.index)*12*365
        annual_rate_of_change1=((tren[0]*pd.to_numeric(trd.index)[-1]+tren[1])-(tren[0]*pd.to_numeric(trd.index)[0]+tren[1]))/(len(trd.index)/12)*365
        annual_rate_of_change2=((trd.iloc[-1].values[0]-trd.iloc[0].values[0])*365)/((trd.index[-1]-trd.index[0]).days/365)
        #print(annual_rate_of_change1,annual_rate_of_change2)
        annual_rate_of_change=tren[0]*365*1000
        print('tr: '+str(annual_rate_of_change))
        print('stderr '+str(tren[4]*365*1000))
        print('r2 '+str(tren[2]))
        print('pval  '+str(tren[3]))
        print('mean '+str(trd.values.mean()))
        print('std '+str(trd.values.std()))
        
        #print((trd.iloc[-1].values[0],trd.iloc[0].values[0]))
        if t=='magnitude':
            tex=str(np.round(annual_rate_of_change,2))+'mgC m$^{-2}$ day$^{-1}$ yr$^{-1}$'
            #if i>=5:
            if typ=='both':
                ax1.text(np.datetime64('2010-06-01'),-0.16,tex,fontsize=10)
            else:
                ax1.text(np.datetime64('2010-06-01'),0.01,tex,fontsize=10)
            #else:
            #    ax1.text(np.datetime64('2013-01-01'),0.055,tex,fontsize=11)
            
        ax1.set_title(title,fontsize=fs)
        
        if typ=='both':
            ax1.text(np.datetime64('1998-08-01'),0.025,'Larger air-sea flux',fontsize=11)
            #if i>2:
            ax1.text(np.datetime64('1998-08-01'),-0.16,'Larger new production',fontsize=11)
        #else:
        #    ax1.text(np.datetime64('1998-08-01'),-0.12,'Larger Biology',fontsize=11)
        
        lp=20
        ax2.set_ylabel('Trend',labelpad=lp,fontsize=10)
        ax2.axhline(0,c='gray')#,linestyle=':')
        ax2.plot(decomp.trend,c=col)
        #ax2.plot(decomp.seasonal,c='k',linestyle=':')
        trd=decomp.trend[decomp.trend.index>=tsd]
        
        tren=trends(ax2,trd.index,trd.values)
        ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax2.grid(axis='y',which='major')
        ax2.grid(axis='x',which='both')
        #ax2.set_ylim([-0.08,0.07])
        if typ=='both':
            ax2.set_ylim([-0.15,0.07])
     
        annual_rate_of_change3=(((tren[0]*pd.to_numeric(trd.index)[-1]+tren[1])-(tren[0]*pd.to_numeric(trd.index)[0]+tren[1]))/len(trd.index))*12*365
        annual_rate_of_change4=((trd[-1]-trd[0])*365)/((trd.index[-1]-trd.index[0]).days/365)
        annual_rate_of_change=tren[0]*365*1000
        #print(annual_rate_of_change3,annual_rate_of_change4)
        #print(trd[-1],trd[0])
        if typ=='both':
            ax2.text(np.datetime64('1998-08-01'),0.025,'Larger air-sea flux',fontsize=11)
        #if i<4:
            ax2.text(np.datetime64('1998-08-01'),-0.12,'Larger new production',fontsize=11)
        #else:
        #    ax2.text(np.datetime64('1998-08-01'),-0.05,'Larger Biology',fontsize=11)
            
        if t=='magnitude':
            tex=str(np.round(annual_rate_of_change,2))+'mgC m$^{-2}$ day$^{-1}$ yr$^{-1}$'
            #if i>=2:
            if typ=='both':
                ax2.text(np.datetime64('2010-06-01'),-0.14,tex,fontsize=10)
            else:
                ax2.text(np.datetime64('2010-06-01'),0.01,tex,fontsize=10)
           #else:
            #    ax2.text(np.datetime64('2013-01-01'),0.02,tex,fontsize=10)
        else:
            tex=str(np.round(annual_rate_of_change,3))+'%/yr'
            ax2.text(np.datetime64('2013-01-01'),-0.03,tex,fontsize=11)
        
              
        ax2.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
        
       
        ax3.axhline(0,c='gray')#,linestyle=':')
        ax3.plot(decomp.seasonal,c=col)
        ax3.set_ylabel('Seasonality',labelpad=lp,fontsize=10)
        ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax3.grid(axis='y',which='major')
        ax3.grid(axis='x',which='both')
        ax3.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
        if typ=='both':
            ax3.set_ylim([-0.06,0.06])
        
        
        ax4.axhline(0,c='gray')
        ax4.scatter(dates,decomp.resid,c=dat.mei,cmap='bwr')
        ax4.set_ylabel('Residuals',labelpad=lp,fontsize=10)
        
        ax4.grid()
        
        ax4.set_xlim([np.datetime64('1997-06-01'),np.datetime64('2020-01-01')])
        ax4.set_ylim([-0.05,0.05])
        seasonaltrend.append(decomp.seasonal)
        
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])        
        ax3.set_xticklabels([])
        #ax4.set_xtick([])
       #Put in the ENSO box regions
        ensofps=['processed/indexes/el_nino_events.csv','processed/indexes/la_nina_events.csv']
        for whichenso,fp in enumerate(ensofps):
            events=pd.read_csv(fp)
            for ev in events.iterrows():
                endm=np.datetime64(ev[1].end).astype('datetime64[M]')
                endm1=endm-np.timedelta64(1,'M')
                start=np.datetime64(ev[1].start).astype('datetime64[M]')
                if start==endm1:
                    pass    
                else:
                    if whichenso==0:
                        #if el nino
                        patchcol='firebrick'
                    else:
                        #if la nina
                        patchcol='deepskyblue'
                    rect=patches.Rectangle((start,-500),endm-start,1000,linewidth=0,alpha=0.4,color=patchcol)
                    ax1.add_patch(rect)
                    rect=patches.Rectangle((start,-500),endm-start,1000,linewidth=0,alpha=0.4,color=patchcol)
                    ax2.add_patch(rect)
                    rect=patches.Rectangle((start,-500),endm-start,1000,linewidth=0,alpha=0.4,color=patchcol)
                    ax3.add_patch(rect)
                    
                    #ax2.add_patch(rect)
                    #rect=patches.Rectangle((start,-50),endm-start,350,linewidthty=0,alpha=0.2,color=patchcol)
                    #ax2.add_patch(rect)
    
            
    plt.tight_layout()     
    plt.savefig('figs/Figure3_decomp_megaplot'+div_factor.name+typ+'.png',dpi=200)
    plt.show()

