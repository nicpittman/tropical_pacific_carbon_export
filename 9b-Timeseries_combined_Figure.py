#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:52:00 2020
@author: Nic Pittman

This script reproduces Figure 2 from Pittman et al., 2019.

It produces a timeseries of New production and Co2 flux for each of the 6 eqpac moorings.
Includes SST, isotherm depth. 
Calculates averages during el nino, la nina and neutral
Also calculates thermocline slope.

Requires:
    datasets/tao/tao_physics/*
    processed/combined_dataset/month_data_exports.nc
    processed/flux/landsch_mooring_co2_flux.nc
    processed/flux/npp.nc
    
Produces:
    processed/results/means.csv
    figs/Figure2_Co2fluxevents+ratio_name+.png
    processed/results/enso_mooring_avg.csv

"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator,AutoMinorLocator



#Not recommended but helps to get the info that we need out
import warnings
warnings.filterwarnings("ignore")


moorings=['110W','125W','140W','155W','170W','165E'][::-1]

#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)


fs=20
fig=plt.figure(figsize=(28,30))#,constrained_layout=True)
gs=fig.add_gridspec(38,1)
means=pd.DataFrame()
# npp_insitu=pd.read_csv('processed/flux/shipboard_npp.csv',index_col='id') #Redundant 
final_mooring_enso=pd.DataFrame()

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp)
print((dat.sel(Mooring=155).co2flux4_land_gmyr/365).min().values) #Min value at 155W (1998)
for i, mooring_name in enumerate(moorings):
    print(mooring_name)
    #JMAco2flux_data=xr.open_mfdataset('processed/flux/JMA_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
    LandSch_co2flux_data=xr.open_mfdataset('processed/flux/landsch_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
    npp=xr.open_mfdataset('processed/flux/npp.nc').sel(Mooring=mooring_name)
    #mooring_obs_npp=npp_insitu[((npp_insitu.mooring.astype(int)>=int(mooring_name[:-1])-1)&(npp_insitu.mooring.astype(int)<=int(mooring_name[:-1])+1))]
    
    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
   
    #Original all NPPS
    #avg_npp=npp[['viirs_eppley','viirs_cbpm','viirs_vgpm','sw_eppley','sw_cbpm','sw_vgpm','mod_cafe','mod_eppley','mod_cbpm','mod_vgpm']].to_dataframe().drop(columns='Mooring').mean(axis=1)
    #Now only CBPM and CAFE
    #avg_npp=npp[['viirs_cbpm','sw_cbpm','mod_cbpm']].to_dataframe().drop(columns='Mooring').mean(axis=1)
    avg_npp=npp[['sw_cafe','mod_cafe']].to_dataframe().drop(columns='Mooring').mean(axis=1)
    
    temps=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/t0n'+mooring_name.lower()+'_dy.cdf')


    dat['Date'].astype('datetime64[M]')
    co2flux=carbon_flux(dat.sss2,dat.sst2,dat.windspeed,None,None,dat.delta_pCO2)[0]
   
    
    ratio=dat.laws2011a#b#f_ratio#11a#f_ratio#laws2011a#.dunne_zeu1#thE_ratio#f_ratio
    ratio_name='laws2011a'#'Laws 2011b'#'11a'#'Laws2000'#'Laws2011a'#'f_ratio'#'Dunne'#'thE-ratio'
    #ratio_name='0.05'
    
    #Need a .T to transform in ax2.contourf
    #MAX 60 day gap fill. All depth filling is fine though.
    temps_selected=temps.sel(depth=slice(0,250))
    temps_selected=temps_selected.interpolate_na(dim='depth',method='linear')
    temps_selected=temps_selected.interpolate_na(dim='time',method='linear',limit=60) #Stop doing major gap fills but small ones are useful.
    temps_selected['depth']=temps_selected['depth']*-1
    temps_selected=temps_selected.resample(time='M').mean().T_20
    temps_selected['time']=temps_selected.time.astype('datetime64[M]')
    start=i*6
    bignext=start+4
    smallnext=bignext+2

    ax1=fig.add_subplot(gs[start:bignext,0])
    ax1.set_title(str(mooring_name),fontsize=fs+4)
  
    #ax1.plot(JMAco2flux_data.Date,moles_to_carbon(JMAco2flux_data.flux.values)/365,linewidth=4,label='Iida CO$_{2}$ flux',c='mediumblue')#'medium_')
    ax1.plot(LandSch_co2flux_data.Date.astype('datetime64[M]'),moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed.values)/365,linewidth=4,label='Landschutzer CO$_{2}$ flux',c='slategray' )
    ax1.plot(dat.Date.astype('datetime64[M]'),dat.co2flux_gmyr/365,linewidth=5,label='in situ CO$_{2}$ flux',c='k')#'mediumblue')
   
    ax1.plot(avg_npp.index,avg_npp.values/1000*ratio.sel({'Date':avg_npp.index}),linewidth=5,label='New Production',c='orangered')
    ax1.axhline(0,linestyle='--',c='k',linewidth=3,alpha=0.8)


    #ax1.scatter(mooring_obs_npp.time.astype('datetime64').values,abs(mooring_obs_npp.pp.values)/1000*ratio.sel({'Date':mooring_obs_npp.time.astype('datetime64[M]').values}),marker='x',c='k',s=40,alpha=0.8,label='_nolegend_')
    #print('Have :'+str(len(mooring_obs_npp))+' insitu NPP observations for '+mooring_name)
    
    #Put a drawdown indicator
    #co2
    co222=moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed/365)
    #drawdown=dat.co2flux4_land_gmyr.where(dat.co2flux4_land_gmyr<0)
    drawdown=co222.where(co222<0)
    draw_dates=drawdown[~drawdown.isnull()].Date.values
    print('Drawdown Mean: '+str(drawdown.mean().values))
    ax1.scatter(draw_dates.astype('datetime64[M]'),np.zeros(len(draw_dates))+0.015,c='r',s=500,marker=11,label='Drawdown')
    #Drawdown indicator
    
    print(mooring_name)
    print('NPP avg: '+str((avg_npp/1000*ratio.sel({'Date':avg_npp.index})).mean()))
    print('NPP std: '+str((avg_npp/1000*ratio.sel({'Date':avg_npp.index})).std()))
    print('Land avg: '+str((moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed.values)/365).mean()))
    print('Land std: '+str((moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed.values)/365).std()))
    
    #print('Iida avg: '+str((moles_to_carbon(JMAco2flux_data.flux.values)/365).mean()))
    # print('Iida std: '+str((moles_to_carbon(JMAco2flux_data.flux.values)/365).std()))
    print('In Situ avg: '+str((dat.co2flux_gmyr/365).mean().values))
    print('In Situ std: '+str((dat.co2flux_gmyr/365).std().values))
    
    epp=pd.DataFrame([dat.mod_eppley.values,dat.sw_eppley.values]).mean()
    cbpm=pd.DataFrame([dat.mod_cbpm.values,dat.sw_cbpm.values,dat.viirs_cbpm.values]).mean()
    vgpm=pd.DataFrame([dat.sw_vgpm.values,dat.mod_vgpm.values,dat.viirs_vgpm.values]).mean()
    cafe=pd.Series(dat.mod_cafe.values).astype(float)
    
    avgs=pd.Series({'NPP Avg':float((avg_npp/1000*ratio.sel({'Date':avg_npp.index})).mean()),
                  'NPP std':float((avg_npp/1000*ratio.sel({'Date':avg_npp.index})).std()),
                  
                  'EPP Avg':float((epp/1000*ratio).mean()),
                  'EPP std':float((epp/1000*ratio).std()),
                  
                  'CBPM Avg':float((cbpm/1000*ratio).mean()),
                  'CBPM std':float((cbpm/1000*ratio).std()),
                  
                  'VGPM Avg':float((vgpm/1000*ratio).mean()),
                  'VGPM std':float((vgpm/1000*ratio).std()),
                  
                  'CAFE Avg':float((cafe/1000*ratio).mean()),
                  'CAFE std':float((cafe/1000*ratio).std()),
                  
                  'Land avg':float((moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed.values)/365).mean()),
                  'Land std':float((moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed.values)/365).std()),
                  #'Iida avg':float((moles_to_carbon(JMAco2flux_data.flux.values)/365).mean()),
                 # 'Iida std':float((moles_to_carbon(JMAco2flux_data.flux.values)/365).std()),
                  'In Situ avg':float((dat.co2flux_gmyr/365).mean().values),
                  'In Situ std':float((dat.co2flux_gmyr/365).std().values)
                  })
    avgs.name=mooring_name
    means=means.append(avgs)
    ax2=fig.add_subplot(gs[(bignext-1):smallnext,0]) #-1 is to put ontop of each other.

    x=ax2.contourf(temps_selected.time.values.astype(np.datetime64),
                   temps_selected.depth.values,temps_selected.squeeze().values.T,
                   levels=np.arange(10,32,2),cmap='RdYlBu_r',extend='both')
    
    ax2.contour(temps_selected.time.values.astype(np.datetime64),
                temps_selected.depth.values,temps_selected.squeeze().T,
                levels=[20],colors='k',linewidths=2)
    
    
    temps_selected1=temps_selected.sel(time=slice('2000-01-01','2018-01-01'))
    thermocli=ax2.contour(temps_selected1.time.values.astype(np.datetime64),
                temps_selected1.depth.values,temps_selected1.squeeze().T,
                levels=[20],colors='k',linewidths=0)
    
    #Get the thermocline contour
    vvy=[]
    vvx=[]
    for ii, seg in enumerate(thermocli.allsegs[0]):
        #plt.plot(seg[:,0], seg[:,1], '.-', label=ii)
        vvx.extend(list(seg[:,0]))
        vvy.extend(list(seg[:,1]))
    df=pd.DataFrame({'DT':vvx,'Depth':vvy})
    df1=df.sort_values('DT')
    from scipy.stats import linregress
    shoaling=linregress(df.DT,df.Depth).slope*365
    shoaling_std=linregress(df.DT,df.Depth).stderr*365
    pval=linregress(df.DT,df.Depth).pvalue
    print('SHOALING AT: '+str(shoaling) + ' : '+str(shoaling_std))
    print('p=' +str(pval))
    #plt.title(mooring_name+': Upper ocean temperatures')
    
    finyear='2020-01-01'
    ax1.set_ylim([-0.03,0.26])
    ax1.set_xlim([np.datetime64(startyear),np.datetime64(finyear)])
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax1.tick_params(labelbottom=False)
    ax1.grid()
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.tick_params(which='major', length=12,direction='in',color='k',top=True)
    ax1.tick_params(which='minor', length=6,direction='in',top=True)
    #ax1.tick_params(which='minor', length=7,direction='inout')
    
    ax2.set_ylim([-250,-0])
    
    ax2.grid()
    ax2.set_yticks([-50,-150,-250])
    ax2.tick_params(axis='both', which='major', labelsize=fs)
    ax2.set_xlim([np.datetime64(startyear),np.datetime64(finyear)])
    ax2.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax2.tick_params(which='x', width=3)
    ax2.tick_params(which='major', length=12,direction='inout',color='k')
    ax2.tick_params(which='minor', length=6,direction='inout')
    
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
                rect=patches.Rectangle((start,0),endm-start,0.5,linewidth=0,alpha=0.2,color=patchcol)
                ax1.add_patch(rect)
                rect=patches.Rectangle((start,-0.5),endm-start,0.5,linewidth=0,alpha=0.2,color=patchcol)
                ax1.add_patch(rect)
                

    #Put relevant labels on
    if i==5:
        ax2.set_ylabel('Depth (m)',fontsize=fs)
        ax1.set_ylabel('gC m$^{-2}$ Day$^{-1}$',fontsize=fs,labelpad=30)    #Started at 16
        cax=fig.add_subplot(gs[smallnext:smallnext+2])
        # cbar=plt.colorbar(x,cax,orientation='horizontal',aspect=50,pad=0.1)
        # cbar.ax.tick_params(labelsize=fs)
        # cbar.set_label('Temperature (C)',fontsize=fs)
        
        ax=inset_axes(cax,'70%', '35%',loc='upper center')#),borderpad=-6.5)#,bbox_to_anchor=(0,0,1,1))
        cbar=plt.colorbar(x,ax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=fs)
        cbar.set_label('Temperature (C)',fontsize=fs)
        cax.set_visible(False)
        ax1.text(np.datetime64('2009-11-01'),0.22,'CP',fontsize=14)
        ax1.text(np.datetime64('2015-09-01'),0.22,'CP',fontsize=14)
    elif i==0:
        ax1.legend(ncol=4,fontsize=fs,loc='upper center')
        ax2.tick_params(labelbottom=False)
    else:
        ax1.text(np.datetime64('2009-11-01'),0.22,'CP',fontsize=14)
        ax1.text(np.datetime64('2015-09-01'),0.22,'CP',fontsize=14)
        ax2.tick_params(labelbottom=False)
   
    #We want to create a table of NINO, NINA, neutral and all CO2 and NP averages for each mooring

    info=dat.to_dataframe()
    info['select_model']=info.f_ratio*info.cbpmmean/1000
    info['co2']=info.co2flux4_land_gmyr/365
    info=info[~np.isnan(info.select_model)]
    
    nino1=info[info.mei>0.5].mean()
    nina1=info[info.mei<-0.5].mean()
    neutral1=info[(info.mei<0.5)&(info.mei>-0.5)].mean()
    
    x ={'Nino CO2':nino1.co2,
        'Nino NP':nino1.select_model,
        'Nina CO2':nina1.co2,
        'Nina NP':nina1.select_model,
        'Neutral CO2':neutral1.co2,
        'Neutral NP':neutral1.select_model}
    ensoavgs=pd.Series(x,name=mooring_name)
    final_mooring_enso=final_mooring_enso.append(ensoavgs)
    
    
plt.tight_layout()
if startyear==str(1997):
    plt.savefig('figs/Figure2_Co2fluxevents'+ratio_name+'.png',format='png',dpi=100)
else:
    plt.savefig('figs/fig2_Co2fluxevents_alltime.png',format='png',dpi=100)
plt.show()

#Can calculate correlations like:
dat.to_dataframe().corr().dic.abs().sort_values()[::-1].head(25)
mm=means.T
print(mm.to_string())
mm.to_csv('processed/results/means.csv')

final_mooring_enso.to_csv('processed/results/enso_mooring_avg.csv')
print(final_mooring_enso)
#final_mooring_enso.plot()
