#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 11:52:00 2020
@author: Nic Pittman

This script reproduces Figure 2 from Pittman et al., 2021.

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
from scipy.stats import linregress
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error


def check_single_bias(truth,model):
    bias=((model-truth)/truth)*100
    abs_error=abs(bias)
            
    logbias=10**(np.nanmean(np.log10(model)-np.log10(truth)))

    medlogbias=10**(np.nanmedian(np.log10(model)-np.log10(truth))) 
    
    mae=10**(np.nanmean(abs(np.log10(model)-np.log10(truth))))
    
    med_ae=10**(np.nanmedian(abs(np.log10(model)-np.log10(truth))))

    return logbias,medlogbias,mae,med_ae,bias,abs_error
    
    
#Not recommended but helps to get the info that we need out
import warnings
warnings.filterwarnings("ignore")


### OPTIONAL CHANGE
plot_mooring_CO2=1 #CHANGE to 1 to check how the mooring data lines up.



moorings=['110W','125W','140W','155W','170W','165E'][::-1]

#Change startyear to 1980 for full timeseries and will auto save _alltime.
startyear=str(1997)

check_lag_corr_asf=[]
check_lag_corr_np=[]




# # %% Show distribution
# plt.figure(figsize=(10,10))
# plt.subplot(211)
# plt.grid()
# sns.violinplot(data=np.array(d).T,orient='h')#x='Mooring',y='co2flux',data=d.to_dataframe().T)
# plt.title('In situ CO2 flux (gC/m2/day)')
# plt.yticks([0,1,2,3,4,5],['110W','125W','140W','155W','170W','165E'])
# plt.xlim([-0.025,0.2])
# plt.subplot(212)
# plt.grid()
# plt.title('Landschutzer CO2 flux')
# sns.violinplot(data=np.array(d1).T,alpa=2,orient='h')
# plt.yticks([0,1,2,3,4,5],['110W','125W','140W','155W','170W','165E'])
# plt.xlim([-0.025,0.2])
# plt.show()

# %%





#190 mm x 230 mm

fs=12
fig=plt.figure()#figsize=(28,30))#,constrained_layout=True)
fig.set_size_inches(19,24)
gs=fig.add_gridspec(38,8)
means=pd.DataFrame()
# npp_insitu=pd.read_csv('processed/flux/shipboard_npp.csv',index_col='id') #Redundant 
final_mooring_enso=pd.DataFrame()

fp='processed/combined_dataset/month_data_exports.nc'

#dat=xr.open_mfdataset(fp)

#print((dat.sel(Mooring=155).co2flux4_land_gmyr/365).min().values) #Min value at 155W (1998)
for i, mooring_name in enumerate(moorings):
    print(mooring_name)
    #JMAco2flux_data=xr.open_mfdataset('processed/flux/JMA_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
    LandSch_co2flux_data=xr.open_mfdataset('processed/flux/landsch_mooring_co2_flux.nc').rename({'time':'Date'}).sel(Mooring=mooring_name)
    LandSch_co2flux_data['Date']=LandSch_co2flux_data['Date'].astype("datetime64[M]")
    npp=xr.open_mfdataset('processed/flux/npp.nc').sel(Mooring=mooring_name)/12
    #mooring_obs_npp=npp_insitu[((npp_insitu.mooring.astype(int)>=int(mooring_name[:-1])-1)&(npp_insitu.mooring.astype(int)<=int(mooring_name[:-1])+1))]
    
    ty='month' #Actually month though need to fix this.
    fp='processed/combined_dataset/'+ty+'_data_exports.nc'
    try:
        dat=xr.open_mfdataset(fp).sel(Mooring=int(mooring_name[:-1]))
    except:
        dat=xr.open_mfdataset(fp).sel(Mooring=195)
   
    ratio=dat.laws2011a#b#f_ratio#11a#f_ratio#laws2011a#.dunne_zeu1#thE_ratio#f_ratio
    ratio_name='laws2011a'#'Laws 2011b'#'11a'#'Laws2000'#'Laws2011a'#'f_ratio'#'Dunne'#'thE-ratio'

    #Original all NPPS
    #avg_npp=npp[['viirs_eppley','viirs_cbpm','viirs_vgpm','sw_eppley','sw_cbpm','sw_vgpm','mod_cafe','mod_eppley','mod_cbpm','mod_vgpm']].to_dataframe().drop(columns='Mooring').mean(axis=1)
    #Now only CBPM and CAFE
    #avg_npp=npp[['viirs_cbpm','sw_cbpm','mod_cbpm']].to_dataframe().drop(columns='Mooring').mean(axis=1)
    avg_npp=npp[['sw_cafe','mod_cafe']]
    avg_npp=avg_npp.sel(Date=slice(ratio.Date.min(),ratio.Date.max()))
    avg_npp=avg_npp.to_dataframe().drop(columns='Mooring').mean(axis=1)
    
    temps=xr.open_mfdataset('datasets/tao/tao_physics/'+mooring_name+'/t0n'+mooring_name.lower()+'_dy.cdf')


    dat['Date'].astype('datetime64[M]')
    co2flux=carbon_flux(dat.sss2,dat.sst2,dat.windspeed,None,None,dat.delta_pCO2)[0]
   
    
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

    ax1=fig.add_subplot(gs[start:bignext,0:7]) #plt.subplot(111)#
    if i==0:
        ax1.set_title(str(str(mooring_name[0:3]+'\xb0'+mooring_name[3:4])),fontsize=fs+2)
    else:
        ax1.annotate(str(str(mooring_name[0:3]+'\xb0'+mooring_name[3:4])),
                    [np.datetime64('2008-02-01'),17],
                    fontsize=fs+2)
    #ax1.plot(JMAco2flux_data.Date,moles_to_carbon(JMAco2flux_data.flux.values)/365,linewidth=4,label='Iida CO$_{2}$ flux',c='mediumblue')#'medium_')
    
    
    ax1.plot(LandSch_co2flux_data.Date.astype('datetime64[M]'),
             (LandSch_co2flux_data.fgco2_smoothed.values*1000)/365,
             linewidth=3,
             label='Landschutzer CO$_{2}$ flux',
             c='slategray')
    
    #CHECK HERE TO PLOT IN SITU CO2 FLUX
    if plot_mooring_CO2==1:
        ax1.plot(dat.Date.astype('datetime64[M]'),(dat.co2flux_gmyr/365)*1000/12,linewidth=4,label='in situ CO$_{2}$ flux',c='mediumblue')#'mediumblue')
   
    
    #Calculate bias between insitu and product
    landCO2=moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed)/365
    landCO2['Date']=LandSch_co2flux_data.Date.astype('datetime64[M]')
    in_situCO2=dat.co2flux_gmyr/365
    bias=(landCO2-in_situCO2).mean().values
    print('LANDSCHUTZER BIAS: '+str(bias))
    ax1.plot(avg_npp.index,(avg_npp.values)*ratio.sel({'Date':avg_npp.index}),linewidth=3,label='New Production',c='orangered')
    ax1.axhline(0,linestyle='--',c='k',linewidth=3,alpha=0.8)

    asf=(moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed)/365)
    npp=avg_npp.values/12*ratio.sel({'Date':avg_npp.index})
    check_lag_corr_asf.append(asf.sel(Date=slice('1998-01-01','2019-12-31')).values)
    check_lag_corr_np.append(npp.sel(Date=slice('1998-01-01','2019-12-31')).values)
    


    #Put a drawdown indicator
    #co2
    co222=moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed/365)
    #drawdown=dat.co2flux4_land_gmyr.where(dat.co2flux4_land_gmyr<0)
    drawdown=co222.where(co222<0)
    draw_dates=drawdown[~drawdown.isnull()].Date.values
    print('Drawdown Mean: '+str(drawdown.mean().values))
    #if ~np.isnan(drawdown.mean()):
    #    print(drawdown.values)
    ax1.scatter(draw_dates.astype('datetime64[M]'),np.zeros(len(draw_dates))+0.1,c='r',s=500,marker=11,label='Drawdown')
    #Drawdown indicator
    
    
    finyear='2020-01-01'
    ax1.set_ylim([-6,20])
    ax1.set_yticks([0,5,10,15,20])
    ax1.set_xlim([np.datetime64(startyear),np.datetime64(finyear)])
    ax1.tick_params(axis='both', which='major', labelsize=fs)
    ax1.tick_params(labelbottom=False)
    ax1.grid()
    ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax1.tick_params(which='major', length=12,direction='in',color='k',top=True)
    ax1.tick_params(which='minor', length=6,direction='in',top=True)
    #ax1.tick_params(which='minor', length=7,direction='inout')
    
    
    
    
    
    print(mooring_name)
    print('NPP avg: '+str((avg_npp/12*ratio.sel({'Date':avg_npp.index})).mean()))
    print('NPP std: '+str((avg_npp/12*ratio.sel({'Date':avg_npp.index})).std()))
    land_avg_dat=(((LandSch_co2flux_data.fgco2_smoothed.values)*1000)/365).mean()
    print('Land avg: '+str(land_avg_dat))
    print('Land std: '+str((((LandSch_co2flux_data.fgco2_smoothed.values)*1000)/365).std()))
    
    #print('Iida avg: '+str((moles_to_carbon(JMAco2flux_data.flux.values)/365).mean()))
    # print('Iida std: '+str((moles_to_carbon(JMAco2flux_data.flux.values)/365).std()))
    situ_avg=(dat.co2flux_gmyr/365).mean().values*1000/12
    print('In Situ avg: '+str(situ_avg))
    print('In Situ std: '+str((dat.co2flux_gmyr/365).std().values*1000/12))
    
    a=(check_single_bias(dat.co2flux_gmyr/365,moles_to_carbon(LandSch_co2flux_data.fgco2_smoothed)/365))
    print('CO2 abs error ='+str(a[5].mean().values))
    
    def perc_err(a,b):
        return ((a-b)/b)*100
    
    def MAPE(Y_actual,Y_Predicted):
        Ya=Y_actual.dropna(dim='Date')
        Yp=Y_Predicted.sel(Date=Ya.Date)
        mape = np.mean(np.abs((Ya - Yp)/Ya))*100
        mape1 = (np.mean(np.abs(Yp)) - np.mean(np.abs(Ya))/np.mean(np.abs(Yp)))*100
        print(mape1.values)
        return mape.values
    
    def perc_err(a,b):
        
        
        return ((a-b)/b)*100
 

    def MAE(Y_actual,Y_Predicted):
        Ya=Y_actual.dropna(dim='Date')
        Yp=Y_Predicted.sel(Date=Ya.Date)
        mae = np.mean(np.abs((Ya - Yp)))
        return mae.values
    
    print(f"HACK PERC {perc_err(situ_avg,land_avg_dat)}")
    aa=((dat.co2flux_gmyr/12)/365)
    bb=((LandSch_co2flux_data.fgco2_smoothed.sel(Date=dat.Date))/365)
    
    mape=MAPE(aa,bb)#((, (LandSch_co2flux_data.fgco2_smoothed.sel(Date=dat.Date))/365)
    mae=MAE(aa,bb)#(((dat.co2flux_gmyr/12)/365), (LandSch_co2flux_data.fgco2_smoothed.sel(Date=dat.Date))/365)

    print(f"CO2 MAPE: {mape} and mae: {mae}")
    epp=pd.DataFrame([dat.mod_eppley.values,dat.sw_eppley.values]).mean()
    cbpm=pd.DataFrame([dat.mod_cbpm.values,dat.sw_cbpm.values,dat.viirs_cbpm.values]).mean()
    vgpm=pd.DataFrame([dat.sw_vgpm.values,dat.mod_vgpm.values,dat.viirs_vgpm.values]).mean()
    cafe=pd.Series(dat.mod_cafe.values).astype(float)
    
    avgs=pd.Series({'NPP Avg':float((avg_npp/12*ratio.sel({'Date':avg_npp.index})).mean()),
                  'NPP std':float((avg_npp/12*ratio.sel({'Date':avg_npp.index})).std()),
                  
                  'EPP Avg':float((epp/12*ratio).mean()),
                  'EPP std':float((epp/12*ratio).std()),
                  
                  'CBPM Avg':float((cbpm/12*ratio).mean()),
                  'CBPM std':float((cbpm/12*ratio).std()),
                  
                  'VGPM Avg':float((vgpm/12*ratio).mean()),
                  'VGPM std':float((vgpm/12*ratio).std()),
                  
                  'CAFE Avg':float((cafe/12*ratio).mean()),
                  'CAFE std':float((cafe/12*ratio).std()),
                  
                  'Land avg':float(((LandSch_co2flux_data.fgco2_smoothed.values)/365).mean())*1000/12,
                  'Land std':float(((LandSch_co2flux_data.fgco2_smoothed.values)/365).std())*1000/12,
                  #'Iida avg':float((moles_to_carbon(JMAco2flux_data.flux.values)/365).mean()),
                 # 'Iida std':float((moles_to_carbon(JMAco2flux_data.flux.values)/365).std()),
                  'In Situ avg':float(((dat.co2flux_gmyr/365)).mean().values)*1000,
                  'In Situ std':float(((dat.co2flux_gmyr/365)).std().values)*1000
                  })
    
 
    
    
    avgs.name=mooring_name
    means=means.append(avgs)
    ax2=fig.add_subplot(gs[(bignext-1):smallnext,0:7]) #-1 is to put ontop of each other.

    x=ax2.contourf(temps_selected.time.values.astype(np.datetime64),
                   temps_selected.depth.values,temps_selected.squeeze().values.T,
                   levels=np.arange(10,32,2),cmap='RdYlBu_r',extend='both')
    
    ax2.contour(temps_selected.time.values.astype(np.datetime64),
                temps_selected.depth.values,temps_selected.squeeze().T,
                levels=[20],colors='k',linewidths=2)
    
    
    temps_selected1=temps_selected.sel(time=slice('2000-01-01','2020-01-01'))
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

    shoaling=linregress(df.DT,df.Depth).slope*365
    shoaling_std=linregress(df.DT,df.Depth).stderr*365
    pval=linregress(df.DT,df.Depth).pvalue
    print('SHOALING AT: '+str(shoaling) + ' ± '+str(shoaling_std) + ' m')
    print('p=' +str(pval))
    #plt.title(mooring_name+': Upper ocean temperatures')
    
    
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
    #ep_events and cp events
    ensofps=['processed/indexes/el_nino_events.csv','processed/indexes/la_nina_events.csv']
    ensofps=['processed/indexes/el_nino_events.csv','processed/indexes/la_nina_events.csv','processed/indexes/cp_events.csv','processed/indexes/cold_cp_events.csv']
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
                elif whichenso==3:
                    patchcol='navy'
                rect=patches.Rectangle((start,-5),endm-start,25,linewidth=0,alpha=0.3,color=patchcol)
                ax1.add_patch(rect)
                rect=patches.Rectangle((start,-5),endm-start,25,linewidth=0,alpha=0.3,color=patchcol)
                ax1.add_patch(rect)
                

    #Put relevant labels on
    if i==5:
        ax2.set_ylabel('Depth (m)',fontsize=fs)
        ax1.set_ylabel('NP and CO$_{2}$ flux\nmmolC m$^{-2}$ day$^{-1}$',fontsize=fs,labelpad=1)    #Started at 16
        cax=fig.add_subplot(gs[smallnext:smallnext+2,0:7])
        # cbar=plt.colorbar(x,cax,orientation='horizontal',aspect=50,pad=0.1)
        # cbar.ax.tick_params(labelsize=fs)
        # cbar.set_label('Temperature (C)',fontsize=fs)
        
        ax=inset_axes(cax,'70%', '35%',loc='upper center')#),borderpad=-6.5)#,bbox_to_anchor=(0,0,1,1))
        cbar=plt.colorbar(x,ax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=fs)
        cbar.set_label('Temperature (C)',fontsize=fs)
        cax.set_visible(False)

        
    elif i==0:
        ax1.legend(ncol=4,fontsize=fs,loc='upper center')
        ax2.tick_params(labelbottom=False)
    else:
        
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
    
    
    
    
## %%
# Show distribution to somehow plug into the main figure?
# One row for each mooring then three violins for the distribution of each var.

# fp='processed/combined_dataset/month_data_exports.nc'

# dat=xr.open_mfdataset(fp)
    d=(dat.co2flux_gmyr*1000/12)/365
    d1=(dat.co2flux4_land_gmyr*1000/12)/365
    npr=(dat.cafe/12)*dat.laws2011a

# plt.figure(figsize=(3,10))
#for i,x in enumerate(np.arange(5,-1,-1)):
    start=i*6
    bignext=start+4
    smallnext=bignext+2

    ax=fig.add_subplot(gs[start:bignext,7]) #plt.subplot(6,1,i+1)#plt.subplot(111)#
    
    insitu=np.array(d)#[x]
    land=np.array(d1)#[x]
    newpr=np.array(npr)#[x]
    
    
    sns.kdeplot(land,color='grey',fill='grey',alpha=0.5,vertical=True)
    sns.kdeplot(insitu,color='mediumblue',fill='blue',alpha=0.4,vertical=True)
    sns.kdeplot(newpr,color='orangered',fill='orangered',alpha=0.5,vertical=True)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.xticks(ticks=np.arange(0,0.4,0.1),labels=(np.arange(0,0.4,0.1)*100).astype(int))
    #plt.xticks(ax.get_xticks(), ax.get_xticks() * 100)
    ax.set_xlim([0,0.35])
    plt.xlabel('Distribution (%)',fontsize=fs)
    ax.grid(True,which='both')
    ax.set_ylim([-6,20])
    
    # mooring_dist_df=pd.DataFrame({'insitu':insitu,
    #               'landschutzer':land,
    #               'new Production':newpr})
    
    # dummyDf=pd.DataFrame({'insitu':np.nan,'landschutzer':np.nan,'new Production':np.nan},index=[len(mooring_dist_df)+1])
                         
    # df2=mooring_dist_df.append(dummyDf)

    # df2['huecol'] = 0.0
    # df2['huecol'].iloc[-1]= np.nan
    
    # sns.violinplot(data=mooring_dist_df,
    #                #x=mooring_dist_df.keys().values,
    #                orient='v',
    #                ax=ax1,
    #                #hue=['k','gray','orangered'],
    #                split=True,
    #                inner='quartile'
                   
    #                #hue=['k','gray','orangered'],
    #                )
    if i==0:
        plt.title('KDE',fontsize=fs+2)
        #ax.set_yticks([-5,0,5,10,15])
    if i==5:
        plt.legend(['Landschutzer CO$_{2}$','Insitu CO$_{2}$','New Production'],
                   loc='lower middle',
                   fontsize=fs,
                   bbox_to_anchor=(1.03,-0.25),
                   borderaxespad=0)
    # else: plt.xticks([],[])
    # plt.ylim([-0.02,0.26])
    
        
plt.tight_layout()

# %%
if startyear==str(1997):
    #plt.savefig('figs/Figure2_Co2fluxevents'+ratio_name+'.png',format='png',dpi=100)
    try:
        plt.savefig('figs/Figure2.jpeg',dpi=300) #Conda install pilliow needed to save to jpeg.
    except:
        pass
    plt.savefig('figs/vector/Figure2.eps',dpi=300)
    plt.savefig('figs/vector/Figure2.pdf',dpi=300)
    plt.savefig('figs/Figure2.png',dpi=300)
    
else:
    plt.savefig('figs/Figure2_Co2fluxevents_alltime.png',format='png',dpi=100)
plt.show()

#Can calculate correlations like:
dat.to_dataframe().corr().dic.abs().sort_values()[::-1].head(25)
mm=means.T
print(mm.to_string())
mm.to_csv('processed/results/means.csv')

final_mooring_enso.to_csv('processed/results/enso_mooring_avg.csv')
print(final_mooring_enso)
#final_mooring_enso.plot()

check_lag_corr_asf1=check_lag_corr_asf[::-1]
check_lag_corr_np1=check_lag_corr_np[::-1]
moorings1=moorings[::-1]
#Just so it starts with the west


fig=plt.figure(figsize=(12,10))#,constrained_layout=True)
for i in range(len(check_lag_corr_asf)):
    ax=plt.subplot(7,2,i+1)
    hh=ax.xcorr(check_lag_corr_asf1[0],check_lag_corr_asf1[i],maxlags=12,normed=True)
    plt.ylim([0.5,1])
    plt.title('ASF: Mooring '+moorings1[0]+' verse: mooring '+moorings1[i])
    #plt.show()
    coefs=hh[1]
    index_max = max(range(len(coefs)), key=coefs.__getitem__)
    plt.plot([hh[0][index_max],hh[0][index_max]],[0,hh[1][index_max]],c='r')
    print('largest')
    print(hh[0][index_max],hh[1][index_max])
    print('0 lag')
    print(hh[1][3])
    print('Difference')
    print(hh[1][index_max]-hh[1][3])
    print('\n')
print(gs)
#plt.show()

#fig=plt.figure(figsize=(12,10))#,constrained_layout=True)
for i in range(len(check_lag_corr_asf)):
    ax=plt.subplot(7,2,8+i+1)
    hh=ax.xcorr(check_lag_corr_np1[0],check_lag_corr_np1[i],maxlags=12,normed=True)
    plt.ylim([0.75,1])
    plt.title('NP: Mooring '+moorings1[0]+' verse: mooring '+moorings1[i])
    #plt.show()
    coefs=hh[1]
    index_max = max(range(len(coefs)), key=coefs.__getitem__)
    plt.plot([hh[0][index_max],hh[0][index_max]],[0,hh[1][index_max]],c='r')
    print('largest')
    print(hh[0][index_max],hh[1][index_max])
    print('0 lag')
    print(hh[1][3])
    print('Difference')
    print(hh[1][index_max]-hh[1][3])
    print('\n')
print(gs)
plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(12,10))#,constrained_layout=True)
for i in range(len(check_lag_corr_asf)):
    ax=plt.subplot(3,2,i+1)
    #hh=ax.xcorr(check_lag_corr_asf[0],check_lag_corr_asf[i],maxlags=12,normed=True)
    hh=ax.xcorr(check_lag_corr_asf1[i],check_lag_corr_np1[i],maxlags=12,normed=True)
    plt.ylim([0.5,1])
    plt.title('ASF vs new production at: Mooring '+moorings1[i])
    #plt.show()
    coefs=hh[1]
    index_max = max(range(len(coefs)), key=coefs.__getitem__)
    plt.plot([hh[0][index_max],hh[0][index_max]],[0,hh[1][index_max]],c='r')
    print('largest')
    print(hh[0][index_max],hh[1][index_max])
    print('0 lag')
    print(hh[1][3])
    print('Difference')
    print(hh[1][index_max]-hh[1][3])
    print('\n')
print(gs)
plt.show()



