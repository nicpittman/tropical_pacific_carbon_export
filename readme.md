### Tropical Pacific CO2 flux and NPP rates

##### Nicholas Pittman et al., 2021

This set of python scripts provides all data and figures for Pittman et al., 2021; Physical drivers of air-sea CO2 flux and biological draw-down in the equatorial Pacific.

Some things are rough, apologies for the mess. any problems with implementation email me at nic.pittman [at] utas[dot]edu[dot]au

This repository is set up so that all figures can be produced immediately after download (Except the landschutzer CO2 product you need to download yourself, curl script provided below). These figures are scripts 9a-f. 

Firstly, clone the repository

`git clone https://github.com/nicpittman/tropical_pacific_carbon_export.git`

Secondly, you will need to ensure you have all of the dependencies outlined in *requirements.txt*

You can create a conda environment using requirements.txt):

```
Preferred:
conda create --name pacific_carbon --file requirements.txt
(if this doesnt work, add conda forge like:  conda config --append channels conda-forge)

otherwise something like below can also work:

conda create -n pacific_carbon python=3.7 basemap cartopy curl cbsyst ESMPy=7.1.0 xesmf==0.3 h5py ipython Markdown numpy pandas scipy matplotlib spyder=4.0 xarray dask nco netcdf4 statsmodels h5netcdf bs4c


```



### To reproduce figures:

This is assuming you have a unix environment. I have noticed the figures turn out differently if run in terminal compared to using Spyder (likely due to graphics software differences giving the error: MESA-LOADER: failed to load driver swrast (search paths /usr/lib64/dri) - This may be fixable? I don't know. Spyder is included in the dependencies and should reproduce the figures accurately if used.  Most of the numbers reported in the manuscript are printed to console during figure production, otherwise, some are saved into *processed/results/\*.csv*.

Most of the processed data is provided in order for the plotting functions (9[a-f]) to work. However, two sets of data will need to be downloaded manually. SST and the landschutzer CO2 flux product (and NCO to convert to a usable format). You can perform this as so:

```To download manually:
To be downloaded manually (Also outlined in datasets/what_to_download.txt):
To be run in shell / terminal

ls
For SST:
cd datasets
mkdir sst | curl ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc --output sst/sst.mnmean.nc

Landschutzer 2020 product:
https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html
Download using curl:
mkdir co2 | mkdir co2/landschutzer_co2 | curl https://data.nodc.noaa.gov/ncei/ocads/data/0160558/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc --output co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc

Use NCO to convert the time variable name (t to date) so they can be opened by xarray.  #Need NCO so either download or module load nco
cd co2/landschutzer_co2
ncrename -v date,t spco2_MPI-SOM_FFN_v2020.nc 

TAO data is included here in datasets/tao/tao_physics/*. I was not going to include, but hard to make reproducible without. You can find the data here:
https://www.pmel.noaa.gov/tao/drupal/disdel/
All variables at Moorings 110W,125W, 140W,155W,170E,165E
```



The following data files are included here for the reproduction of figures:

    # All files are created during the cleanup script.
    # These files here are monthly and 1 degree averaged from oregon state.
    processed/flux/avg_npp_rg_cbpm.nc 					
    processed/flux/avg_npp_rg_cafe.nc					
    processed/flux/avg_npp_rg_vgpm.nc					
    processed/flux/avg_npp_rg_eppley.nc 				
    
    #These files are produced from multiple data sources, see full pipeline for more information
    processed/combined_dataset/month_data.nc			
    processed/combined_dataset/month_data_exports.nc 	
    
    processed/flux/fratios.nc 							
    processed/flux/pco2grams_eq1.nc 					
    processed/flux/tpca.nc 								
    processed/earth_m2.nc							
    processed/flux/landsch_mooring_co2_flux.nc 		
    processed/flux/npp.nc								
    processed/indexes/el_nino_events.csv				
    processed/indexes/la_nina_events.csv

After you have installed the dependencies, downloaded SST and Landschutzer CO2, in Terminal, you should be able to run (or preferably in Spyder for correct formatting):

```
python 9a ... (tab)
python 9b ...
python 9c ...
python 9e ...
python 9f ...
```

All numbers used in the manuscript are printed in console.

Note, Figure 4/D uses basemap, a depreciated package. A rough version of this figure is provided using cartopy but is not publication quality, so I recommend using the basemap version if possible. 

### To reprocess entire pipeline:

NOTE: In total, after downloading and processing all required data the folder will be  ~**50gb**., and possibly more during Script 1, during the downloading and consolidating of NPP data.

The entire the pipeline is provided for reproducible research, however is provided as is. I have tried to make the flow as streamlined as possible but it is likely that there will be unique problems on different systems. Note, Script 1 takes a long time to download 30gb of NPP data. Scripts 3-6 can be run (in order) at the same time as script 1. Script 7a and b (parallelised version), combines everything and takes a long time to process (6 hours per mooring, or ~6 hours x 6 cores total for 7b)- so warning here. 

Two dataset types are produced. Firstly, 6 Moorings at 110W, 125W, 140W, 155W, 170W, 165E provide in-situ temperature, windspeed, thermocline depth, pCO2 in water and atmosphere among other variables. These have been compiled with different NPP products and algorithms, SST and different CO2 flux products. This data used in the paper uses a 1 degree, 1 month resolution, however products for all, daily and weekly are also produced throughout this process. 

Secondly, all primary productivity models, f-ratios, co2 flux and sst products are gridded to 1 degree grids and monthly grids since September 1997. You will need xesmf for these to work properly (in cleansup script). These are used to produce Figure 4/D.

#### Files

Scripts are organised to be run in numerical order. 

1. Downloads Primary Production data from www.science.oregonstate.edu/ocean.productivity/

   1. Besides script 7, this is probably the slowest. Downloads ~30gb. 

2. Cuts the NPP data into mooring csv files. 

   1. **NOTE** This script is not actually essential, has been provided in processed/npp_mooring_timeseries.
   2. You may have problems as chl needs to be downloaded manually, see 6. These files are provided in the repository, however they can be reprocessed using the bottom half of Script 2. This is normally commented out, as downloading of all the chlorophyll and TPCA files can be time consuming. This script will need some small changes on a new system to work as intended (ie file paths). See 6 for more information about downloading the relevant data if anything is missing.

3. Downloads CO2 flux for several products. (Can run scripts 3-5 simultaneous to script 1 in a separate console)

4. Downloads mooring CO2 data (Including Japanese JMA product, downloaded automagically)

5. Downloads mooring physics data and combines with the mooring CO2 data.

6. Is a text file with a list of other datasets that need to be manually downloaded.

   Chl needs to be downloaded manually from (or if you have it locally somewhere)

   - https://oceandata.sci.gsfc.nasa.gov/ 
   - TPCA: https://researchdata.ands.org.au/tropical-pacific-chlorophyll-reprocessing-v10/1438905 (Script provided in 6 to download this, but only necessary for Script 2 to run.)
   - SIMPLE-TRIM export production model from  https://tdevries.eri.ucsb.edu/models-and-data-products/ to datasets/exports
   - You can try this:
     cd datasets | curl https://tdevries.eri.ucsb.edu/wp-content/uploads/2018/03/SIMPLE_TRIM_output.nc --output SIMPLE_TRIM_output.nc
     however I sometimes have expired certificate problems. Might need to download manually.

7. a - Processes all of this data into the mooring timeseries. This is inefficient (will take days to run), and I would implement this differently if I was to rewrite this, however to keep data resolution as high as possible, this method is effective. Because it is so inefficient, a parallelised version has been produced to be run on the Australian supercomputer GADI but still 6 cores and 6 hours. Oops for pandas inefficiencies.

   1. run_combiner.sh will submit 7b to the supercomputer. This may need modification to work on your system.

8. Is a cleanup script with miscellaneous cleanup functions that i have added on add hoc during development. For example, CAFE was released late in the development process, and a function will plug this (and SST) into the Mooring timeseries.  Also finds ENSO events and another function converts the mooring csvs into an easy to use netcdf file. These have been debugged but may still contain system specific problems. 

9. Are plotting scripts as discussed above.

   a. Figure 1 - Comparison of NPP, and fratios at mooring locations

   b. Integrated plot showing timeseries for CO2 and the CAFE * Laws2011a product.

   c. Seasonal decomposition. This script also includes Figure 4, ENSO and seasonality.

   d. Spatial maps for new production, CO2 flux and the difference between the two. Provides mean, trends and standard deviation.

   e. Mean seasonality and ENSO for different variables.

   f. Compares the overall volumes of CO2 removal from the equatorial Pacific (East, West, Central and overall).

Other scripts:

- Carbon_math is a series of functions that make it easy to convert between mol and gC
- 10-windspeed is a quick correlation assessment of the windspeed vs new production and air-sea flux, with r2 used in the paper.



#### Notes

- Things have changed around a lot, with some redundant code. I have tried to refactor where possible, but still some messiness and possible system dependent issues. 
- All data for the figures is provided. However, all scripts to process the data yourself are provided. 
- All data should go into datasets/ This is then processed into processed/ as the files we want.
- I have tried my best to provide at least the minimum standard of code for reproducibility. 
