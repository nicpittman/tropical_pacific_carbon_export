#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:09:00 2020
@author: Nic Pittman

Reproduces Figure 1 from Pittman et al., 2021.
Requires the full pipeline run (Specifically, 8 cleanup script npp average regridder.)

Requires:
    processed/flux/avg_npp_rg_vgpm.nc
    processed/flux/avg_npp_rg_cbpm.nc
    processed/flux/avg_npp_rg_eppley.nc
    processed/flux/avg_npp_rg_cafe.nc

    processed/combined_dataset/month_data_exports.nc
    processed/flux/fratios.nc
Saves to: 
   figs/Figure1.png
"""
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from carbon_math import *

# %% From NPP ratio compare site.py

lw=2
#create_npp_avgs
plt.figure(figsize=(7,12))
ax=plt.subplot(411)        
vgpm=xr.open_dataset('processed/flux/avg_npp_rg_vgpm.nc').avg_npp
cbpm=xr.open_dataset('processed/flux/avg_npp_rg_cbpm.nc').avg_npp
eppley=xr.open_dataset('processed/flux/avg_npp_rg_eppley.nc').avg_npp
cafe=xr.open_dataset('processed/flux/avg_npp_rg_cafe.nc').avg_npp

ll=220 #140W

fp='processed/combined_dataset/month_data_exports.nc'
dat=xr.open_mfdataset(fp).sel(Mooring=int(140)) 

title=['vgpm','cbpm','eppley','cafe']
#ll=220
vg=vgpm.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')/1000
cb=cbpm.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')/1000
ep=eppley.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')/1000
cf=cafe.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')/1000


# %% Part One, compare NPP estimates
ax.plot(vg.time,vg.values,label='VGPM')#,linestyle='--')
ax.plot(ep.time,ep.values,label='Eppley')#,linestyle='--'
ax.plot(cb.time,cb.values,label='CbPM')
ax.plot(cf.time,cf.values,linewidth=2,c='k',label='CAFE')
#ax.plot(dat.modis_tpca.Date,np.sqrt(dat.modis_tpca),label='sqrt chl')
print('Panel a Averages')
print('vgpm: '+str(np.round(vg.values.mean(),3)))
print('eppley: '+str(np.round(ep.values.mean(),3)))
print('cbpm: '+str(np.round(cb.values.mean(),3)))
print('cafe: '+str(np.round(cf.values.mean(),3))+'\n')

combine=cb.to_dataframe()
combine['cf']=cf.to_dataframe().avg_npp

avg=np.nanmean(combine.iloc[:,-2:],axis=1)
#ax.plot(combine.index,avg,c='k',label='average: '+str(np.round(avg.mean(),3)))
#ax.axhline(0.9,c='gray',linestyle=':')

ax.legend(loc='upper right',ncol=2)
ax.set_title('a) Primary production estimates',loc='left') #Primary Production estimates at 0$^\circ$N 140$^\circ$W
ax.set_ylabel('gC m$^{-2}$ day$^{-1}$')
ax.set_xlabel('Year')

model=cafe

ratios=xr.open_mfdataset('processed/flux/fratios.nc')

f_ratio=ratios.laws2000#(0.62-(0.02*sst))
th_e_ratio=ratios.henson2011#(0.23*np.exp(-0.08*sst))
laws2011a=ratios.laws2011a#((0.5857-0.0165*sst)*model)/(51.7+model)
laws2011b=ratios.laws2011b#0.04756*(0.78-((0.43*sst)/30))*model**0.307
dunne2005=ratios.dunne2005
trim=ratios.trim
trim_std=ratios.trim_std
modd=model.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')/1000
  
f=f_ratio.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
th=th_e_ratio.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
l11a=laws2011a.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
l11b=laws2011b.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
dunne=dunne2005.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
SIMPLETRIM=trim.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')
trim_std_loc=trim_std.sel(lat=0,method='nearest').sel(lon=ll,method='nearest')

# %% Part Two, compare f-ratio estimates
ax=plt.subplot(412)
ll1=plt.plot(dunne.time,dunne,label='Dunne 2005',linewidth=lw)
ll2=plt.plot(th.time,th,label='Henson 2011',linewidth=lw)#,linestyle='--')
ll3=plt.plot(f.time,f,label='Laws 2000') 
#ll4=plt.scatter(np.datetime64('1997-09-01'),SIMPLETRIM.values,label='DeVries & Webber 2017',c='m',marker='*',s=70)

ll4=plt.errorbar(np.datetime64('1997-09-01'),SIMPLETRIM.values,yerr=trim_std_loc.values,label='DeVries & Webber 2017 (1 std)',c='m',marker='o',ms=6,linewidth=2)

#ll5=plt.plot(l11b.time,l11b,label='Laws2011b',c='r')
ll5=plt.plot(l11a.time,l11a,label='Laws2011a',c='k')
#plt.scatter(np.datetime64('1997-09-01'),SIMPLETRIM.values+(trim_std_loc.values*3),label='1DeVries & Webber 2017',c='m',marker='_',s=70)
#plt.scatter(np.datetime64('1997-09-01'),SIMPLETRIM.values-(trim_std_loc.values*3),label='2DeVries & Webber 2017',c='m',marker='_',s=70)

ax.set_ylabel('f-ratio')
ax.set_xlabel('Year')
plt.title('b) f-ratios',loc='left') # Comparison of f ratios at  0$^\circ$N 140$^\circ$W

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,1,2,4,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],ncol=3,loc='upper right')


print('Panel b: f-ratio averages')
print('l11b: '+str(l11b.mean().values))
print('l11a: '+str(l11a.mean().values))
print('l2000: '+str(f.mean().values))
print('hens: '+str(th.mean().values))
print('dunne: '+str(dunne.mean().values))
print('SimpleTRIM '+str(SIMPLETRIM.mean().values))
print('SimpleTRIM std '+str(trim_std_loc.values))


# %%Part three (c) here, plots a hovmoller of f-ratio across the pacific and time.
ax=plt.subplot(413)
hov=laws2011a.sel(lat=slice(-5,5),lon=slice(150,280.5)).mean(dim='lat').load()
hovmol=hov.where((hov>hov.quantile(0.0003))&(hov<hov.quantile(0.9999)))

#hov=ax.contourf(f_ratio.lon,f_ratio.time.astype('datetime64[M]').values,f_ratio.sel(lat=slice(5,-5)).mean(dim='lat'))
hov=ax.contourf(hov.lon,laws2011a.time.astype('datetime64[M]').values,hovmol,levels=np.arange(0.06,0.22,0.01))

#Quick anti-aliasing fix as per: https://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
for c in hov.collections:
    c.set_edgecolor("face")

print('\nMin f-ratio = '+str(hovmol.min().values))
print('Max f-ratio = '+str(hovmol.max().values))
print('Mean f-ratio = '+str(hovmol.mean().values))
print('Mean (all) f-ratio = '+str(laws2011a.mean().values))


cb=plt.colorbar(hov)
cb.set_label('f-ratio')
plt.axvline(180,c='k',linestyle=':')
plt.title('c) Laws 2011a f-ratio',loc='left') # Laws 2011b f ratio 5N-5S
plt.xlabel('Longitude')
plt.ylabel('Time')
labels=[None,'160$^\circ$E','180$^\circ$E','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W']
ax.set_xticklabels(labels)

# %% Part Four, to plot the new production estimates. ax2=plt.subplot(414)

#Co2 can be plotted here but is omitted
# landsch_fp='processed/fluxmaps/landshutzer.nc'
# land_pac=xr.open_dataset(landsch_fp)
# land_pac=moles_to_carbon(land_pac.fgco2_smoothed)/365
# land_site=land_pac.sel(lat=0,method='nearest').sel(lon=ll,method='nearest').sel(time=slice(modd.time.min().values,modd.time.max().values))
# land_site['time']=land_site.time.astype('datetime64[M]')
#ax2.plot((land_site-f*modd).time.astype('datetime64[M]'),land_site-f*modd,label='Laws2000',linewidth=lw)
#ax2.plot((land_site-th*modd).time,land_site-th*modd,label='Henson2011',linewidth=lw)
#ax2.plot((land_site-l11a*modd).time,land_site=l11a*modd,label='Laws2011a',linewidth=lw)
#ax2.plot((land_site-l11b*modd).time,land_site-l11b*modd,label='Laws2011b',linewidth=lw)
ax2=plt.subplot(414)
ax2.plot(l11b.time,dunne*modd,label='Dunne 2005')
ax2.plot(th.time,th*modd,label='Henson 2011')#,linestyle='--')
ax2.plot(f.time,f*modd,label='Laws 2000')
#ax2.plot(l11a.time,l11a*modd,label='Laws2011a',linewidth=lw)
ax2.plot(modd.time,SIMPLETRIM*modd,label='DeVries & Webber 2017',c='m')
#ax2.plot(l11b.time,l11b*modd,label='Laws 2011b',linewidth=lw,c='r')
ax2.plot(l11a.time,l11a*modd,label='Laws 2011a',linewidth=lw,c='k')

#ax2.plot(land_site.time,land_site,c='k',label='CO2 flux',linewidth=lw-1)
print('\nPanel d: export flux ratios')

print('l11b: '+str((l11b*modd).mean().values))
print('l11a: '+str((l11a*modd).mean().values))
print('l2000: '+str((f*modd).mean().values))
print('hens: '+str((th*modd).mean().values))
print('dunne: '+str((dunne*modd).mean().values))


ax2.set_ylabel('gC m$^{-2}$ day$^{-1}$')
ax2.set_xlabel('Year')
ax2.legend(ncol=3)
ax2.set_title('d) New production estimates',loc='left')  #Comparison of Export Flux at 0$^\circ$N 140$^\circ$W




plt.tight_layout()
plt.savefig('figs/vector/Figure1.eps',dpi=300)
plt.savefig('figs/vector/Figure1.pdf',dpi=300)
plt.savefig('figs/Figure1.png',dpi=100)
try:
    plt.savefig('figs/Figure1.jpeg',dpi=300) #Conda install pillow needed to save to jpeg.
except:
    pass
plt.show()
