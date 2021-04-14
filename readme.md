## Tropical Pacific CO<sub>2</sub> flux and new production rates

###### Pittman et al., 2021; Physical drivers of air-sea CO<sub>2</sub> flux and biological draw-down in the equatorial Pacific (DOI: xxxx).

#### 	1. Overview

This set of python scripts provides all data and figures for this paper. I originally intended for this to be structured and highly reproducible. However, the manuscript changed focus significantly, so there is a lot of redundant and vestigial code. I have cleaned where possible but some things are rough, apologies for the mess. This readme is as detailed as possible. If you have any problems with implementation, please email me.

This repository is organised so that all figures can be produced immediately after download (except you will need to download the Landschutzer CO<sub>2</sub> , windspeed and precipitation products, which you will need to download yourself. A curl script is provided below). The figure scripts are named 9a-f. 

`carbon_math.py` is a useful Python function which can be used to calculate CO<sub>2</sub> fluxes through Scmidt Number and Solubility as described in Wanninkhof, R. (2014).

Due to the manuscript changing focus over time, the mooring compilation step is time consuming and is almost, but not entirely redundant. Thus it is still an essential step but the value of it now is low, and is only used for Figure 2. My recommendation is to just use the supplied processed files.



##### 2. Getting Conda working

Firstly, I assume you have a Unix environment. Clone the repository to where you want.

`git clone https://github.com/nicpittman/tropical_pacific_carbon_export.git`

You will  need to ensure you have all of the pangeo dependencies outlined in *requirements.txt*

You can create an Anaconda environment (with conda-forge) using *requirements.txt* (There may be xESMF compatability issues):

```
conda config --append channels conda-forge
conda create --name pacific_carbon --file requirements.txt
```

If it breaks for some reason, something like below can also work, however the above method is preferred (Spyder not essential):

```
conda create -n pacific_carbon python=3.7 ESMPy=7.1.0 xesmf==0.3 spyder=4.0 xarray=0.15 h5py ipython Markdown numpy pandas scipy matplotlib dask nco netcdf4 statsmodels h5netcdf bs4c joblib pillow lxml basemap cartopy curl cbsyst
```

Otherwise, you can try removing the library versions, but this setup has worked on my system. Things might break if you use different versions than those prescribed. 

You should then activate the environment like so:
```
conda activate pacific_carbon
```


### 	3. Reproduce figures 1-6:

I have noticed the figures turn out differently if run in terminal compared to using Spyder (likely due to graphics software differences, on the supercomputer i get the error: `MESA-LOADER: failed to load driver swrast (search paths /usr/lib64/dri)` - This may be fixable? I don't know. Spyder is included in the dependencies and should reproduce the figures accurately if used.  The numbers reported in the manuscript are printed to console during figure production, and some are saved into *processed/results/\*.csv*.

Figures are embedded at the bottom of this document for reference. 

##### 3.1 Data to download

Most of the processed data is provided in order for the plotting functions (scripts 9[a-f]) to work. However, three sets of data will need to be downloaded manually. SST and the Landschutzer CO<sub>2</sub> flux product (and NCO to convert to a usable format) and CMAP precipitation data. You can perform this by running `sh datasets/to_download.sh`:

**Note** You will need NCO installed for the script to run as it needs to convert the Landschutzer dataset time variable name (t to date) so they can be opened by xarray. NCO is included in the conda environment.

TAO data is included here in `datasets/tao/tao_physics/*`. You can update it yourself here (but not essential, and a little time consuming): https://www.pmel.noaa.gov/tao/drupal/disdel/ using all variables at equatorial Moorings 110W,125W, 140W,155W,170E,165E



**Important** CCMP wind speeds is a little more challenging to download. Main web page is located at: http://www.remss.com/measurements/ccmp/

This download link was used:
http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_3387_f2e3_e359.nc?uwnd[(1997-01-01):1:(2017-12-01T00:00:00Z)][(-20):1:(20)][(0.125):1:(359.875)],vwnd[(1997-01-01):1:(2017-12-01T00:00:00Z)][(-20):1:(20)][(0.125):1:(359.875)],nobs[(1997-01-01):1:(2017-12-01T00:00:00Z)][(-20):1:(20)][(0.125):1:(359.875)]

Which came from a chain of links: http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_3387_f2e3_e359.html From http://apdrc.soest.hawaii.edu/datadoc/ccmp_month_v2.php From: http://www.remss.com/measurements/ccmp/ 

But the first link might not work. So use the first in the chain of links and download 1997-01-01 to 2017-12-01, -20 to 20 and just take the whole globe. Select 'nc' filetype and submit to download. Might take a while. 

You will need to then download 2017 to 2020 monthly data manually from http://data.remss.com/ccmp/v02.0/ (ie http://data.remss.com/ccmp/v02.0/Y2018/M01/CCMP_Wind_Analysis_201801_V02.0_L3.5_RSS.nc) . Its possible to write a curl or even python script to do this. But the above links worked faster.

Once this is done you will need to run `python 8-CCMP_windspeed_process.py`



##### 3.2 Data included

The following data files are included here for the reproduction of figures (All files are created during the 8. Cleanup script):

    #These files here are monthly and 1 degree averaged from oregon state university website.		
    processed/flux/avg_npp_rg_cbpm.nc 	
    processed/flux/avg_npp_rg_cafe.nc					
    processed/flux/avg_npp_rg_vgpm.nc					
    processed/flux/avg_npp_rg_eppley.nc 	
    
    #These files are produced from multiple data sources, see full pipeline for more information.
    processed/flux/fratios.nc 							
    processed/flux/pco2grams.nc
    processed/flux/pco2grams_norm.nc
    processed/flux/tpca.nc						
    processed/flux/landsch_mooring_co2_flux.nc 		
    processed/flux/JMA_mooring_co2_flux.nc
    processed/flux/npp.nc		
    processed/flux/chl.nc
    processed/flux/zeu.nc
    
    processed/combined_dataset/month_data.nc			
    processed/combined_dataset/month_data_exports.nc 	
    #Other preliminary files are also included here
    processed/earth_m2.nc	
    processed/seamask.nc # Processed from Landschutzer 2018 product. It was not included in the 2020 version.
    
    processed/indexes/ep_nino_events.csv
    processed/indexes/cp_nino_events.csv				
    processed/indexes/la_nina_events.csv
    processed/indexes/el_nino_events.csv #all El Nino events, not split by EP/CP
    
    processed/results/carbon_mass.csv #Integrated mass of carbon in each Pacific Box
    processed/results/enso_basin_means.csv #How ENSO changes CO2/new production estimates
    processed/results/enso_mooring_avg.csv #Average ENSO conditions at each mooring
    processed/results/means.csv #Average and std for CO2/primary production at each mooring

##### 3.3 Producing the figures

After you have installed the dependencies, downloaded SST, Landschutzer CO<sub>2</sub> and Precipitation+wind speed data, in terminal, you should be able to run (or preferably in Spyder for correct figure size formatting):

```
python 9a ... (tab)
python 9b ...
python 9c ...
python 9e ...
```



### 	4. Reprocess entire pipeline:

As described above, figures can be reproduced from the provided processed data. However, the entire the pipeline is provided for reproducible research, but is provided as is. I have tried to make the flow as streamlined as possible but it is likely that there will be unique problems on different systems. Since the analysis has changed, the time:reward ratio here is low, but worth including for full reproducability. 

###### Notes: 

- In total, after downloading and processing all required data the folder will be ~80gb (including chlorophyll and new production data), and possibly more during Script 1.  
- Script 1 takes a long time to download 30gb of NPP data. Script 2 then processes this. Script 2 will download TPCA chlorophyll data automatically. You will need to download the NASA chlor_a data manually as described below using the shell script provided.
- Scripts 3-5 can be run (in order) at the same time as script 1 (Separate pipelines [1 & 2, 3-5], and are combined in Script 7). 
- Script 7a and b (parallelised version), combines everything and takes a long time to process. Script 7b is recommended and takes (6 hours per mooring, or ~6 hours x 6 cores total for script 7b) - **so this is a warning here**- I made some pandas inefficiencies. 



Two dataset types are produced.

Firstly, 6 Moorings at 110W, 125W, 140W, 155W, 170W, 165E provide in-situ temperature, windspeed, thermocline depth, pCO<sub>2</sub> in water and atmosphere among other variables. These have been compiled with different NPP products, new production algorithms, SST and different CO<sub>2</sub> flux products. Secondly, data used in the paper uses a 1 degree, 1 month resolution, however products for all, daily and weekly are also produced throughout this process (datasets/combined_data/*nc). 

Secondly, all primary productivity models, f-ratios, co2 flux and sst products are gridded to 1 degree grids and monthly grids since September 1997. You will need the xesmf package for these to work properly (included in the Anaconda environment, and run in the cleanup script). These are used to produce Figure 4 (Script 9d).

##### 4.1 Scripts included and notes for running them

*carbon_math.py* is a series of functions that make it easy to convert between carbon units and calculate carbon outgassing.

Scripts are to be run in numerical order. 

1. Downloads Primary Production data from www.science.oregonstate.edu/ocean.productivity/

   1. Besides Script 7, this is probably the slowest. Downloads ~30gb. 

2. Cuts the NPP data into mooring csv files. 

   1. **NOTE** This script is not actually essential, processed data is provided in processed/npp_mooring_timeseries. The processing here can be pretty memory / CPU intensive. I actually cut out the tropics so it wouldn't get killed by the supercomputer login node for too much memory.

   2. NASA chlor_a needs to be downloaded manually. If you already have it stored locally, make sure you change the path variables in script 2. To download, A script in `datasets/chl/download_NASA_chlora.sh` is provided. You will need authorisation cookie as described at: https://oceancolor.gsfc.nasa.gov/data/download_methods/. Sourced from https://oceandata.sci.gsfc.nasa.gov/ 

   3. Once you create an account with earthdata, set the cookie like this. 

   4. ```
      echo "machine urs.earthdata.nasa.gov login USERNAME password PASSWD" > ~/.netrc ; > ~/.urs_cookies
      chmod  0600 ~/.netrc
      ```
      and then run the download script like ```sh download_NASA_chlora.sh```
      They will be downloaded to `datasets/chl/*nc` and then moved automatically to the relevant folder `/chlor_a/seawifs` or `/modis`. Note that for some reason, some downloaded files are corrupt. You will need to manually download the broken netcdf files. I don't have an easy fix for this.  

   5. TPCA should download automatically during script 2. It will download to `datasets/chl/TPCA/`. Sourced from: https://researchdata.ands.org.au/tropical-pacific-chlorophyll-reprocessing-v10/1438905 

3. Downloads CO<sub>2</sub> flux for several products. (Can run scripts 3-5 simultaneous to script 1 in a separate console)

4. Downloads mooring CO<sub>2</sub> data (Including Japanese JMA product, downloaded automagically).

   1. This data comes from https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Moorings/Pacific.html and https://www.pmel.noaa.gov/co2/timeseries/. 

5. Downloads mooring physics data and combines with the mooring CO<sub>2</sub> data.

   1. Mooring Data was sourced from https://www.pmel.noaa.gov/tao/drupal/disdel/ but can also be sourced from: https://tao.ndbc.noaa.gov/tao/data_download/search_map.shtml

6. Is a shell file with a list of other datasets that need to be manually downloaded (or you can just run the shell commands outlined here).

   - datasets/zeu/download_eu_monthly.sh Needs to be run. Need to set the ~.netrc file as described on https://oceancolor.gsfc.nasa.gov/data/download_methods/
           
     
   - SIMPLE-TRIM export production model from  https://tdevries.eri.ucsb.edu/models-and-data-products/ to datasets/exports

     - To download, you can use this, or run shell script 6-Data_shell_download_SIMPLETRIM_and_indicies.sh:

       ```
       cd datasets | curl https://tdevries.eri.ucsb.edu/wp-content/uploads/2018/03/SIMPLE_TRIM_output.nc --insecure --output SIMPLE_TRIM_output.nc
       ```

   - You should have already downloaded Landschutzer 2020, however 2018 is needed in cleanup script to make the seamask. Not essential as seamask.nc is included in the repo, and can be commented out in script 8 (when the functions are called). Can be downloaded like 

     ```
     cd datasets/cos/landschutzer_co2 | curl https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc --output spco2_MPI_SOM-FFN_v2018.nc
     
     ncrename -v date,t spco2_MPI_SOM-FFN_v2018.nc
     ```

   - Climate indicies MEI, PDO, EMI, SOI are included in the repository. Scripts to download are included in 6, however not essential, and will need some manual processing in excel to be opened by pandas properly.

   - Regarding EMI. The original data file was downloaded from: http://www.jamstec.go.jp/aplinfo/sintexf/DATA/emi.monthly.txt which has now changed to: http://www.jamstec.go.jp/virtualearth/data/SINTEX/SINTEX_EMI.csv.  This is reflected in the scripts.

7. a - Processes all of this data into the mooring timeseries. This is an inefficient script (ie. will take days to run), and I would implement this differently if I was to rewrite this. **Because it is so inefficient, a parallelised version (7b) has been produced to be run** on the Australian supercomputer GADI/NCI but still takes 6 cores and 6 hours (or can be run on 6 or 3 native cores). Oops for pandas inefficiencies. Be aware here. I recommend using 7b (And submitted with `run_combiner.sh` if you need to queue this up. This shell script will need modification for your system. If you have 6 CPU cores accessible, you can run it you usually would.)

8. Is a cleanup script with miscellaneous cleanup functions that i have added on ad-hoc during development. For example, CAFE was released late in the development process, and a function will plug this (and SST) into the Mooring timeseries.  Another part finds ENSO events, and another function converts the mooring csvs into an easy to use netcdf file. These have been debugged but may still contain system specific problems. 

   1. You may need to change some hard-coded file paths if anything breaks.
   2. The carbon conversion( `carbon_uatm_to_grams()`) takes a long time, and I get killed on the supercomputer login nodes, for I assume memory. This can be made into its own script and submitted very easily if you have this problem. Again, not entirely necessary as the relevant files are provided. Two versions are provided (*pco2_grams.nc* and *pco2_grams_tempcorrected.nc*) . We want the normal and not temperature corrected version.

9. Are plotting scripts as discussed above.

   a. Figure 1 - Comparison of NPP, and f-ratios at mooring locations

   b. Integrated plot showing time series for CO2 and the CAFE * Laws2011a product.

   c. Seasonal decomposition. This script also includes Figure 4, ENSO and seasonality.

   d. Spatial maps for new production, CO2 flux and the difference between the two. Provides mean, trends and standard deviation.

   e. Mean seasonality and ENSO for different variables.

   f. Compares the overall volumes of CO2 removal from the equatorial Pacific (East, West, Central and overall).

10. run_combiner.sh will submit 7b to the supercomputer. You may need to change this depending on your system paths, or you may not need it at all. Depending on the number of cores, you can change the bottom script 7b to use say 3 cores, taking approximately 12 hours. Has been left at default of 6 cores, one for each mooring.

##### 4.2 Other scripts:

- *carbon_math.py* is a series of functions that make it easy to convert between carbon units.
- 10-windspeed is a quick correlation assessment of the windspeed vs new production and air-sea flux, with r<sup>2</sup> used in the paper.
- 9z's are old figures that use the mooring data, and there is also a correlation plot. These were not included in final analysis but could be useful for you?

#### 5. Notes

- Things have changed around a lot, with some redundant code still in the repo. I have tried to refactor where possible, but there is still some messiness and possible system dependent issues.  I have tried to fix this the best I can with Anaconda (see Conda, Section 2). I have tried my best to provide at least the minimum standard of code for reproducible science. 
- All data for the figures is provided. However, all scripts to process the data yourself are also provided for complete reproducibility of this paper 



