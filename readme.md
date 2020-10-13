## Tropical Pacific CO<sub>2</sub>  flux and NPP rates

###### Pittman et al., 2021; Physical drivers of air-sea CO2 flux and biological draw-down in the equatorial Pacific.



#### 	1. Overview

This set of python scripts provides all data and figures for this paper. Some things are rough, apologies for the mess. This readme is as detailed as possible. If you have any problems with implementation, please email me at nic.pittman [at] utas[dot]edu[dot]au

This repository is organised so that all figures can be produced immediately after download (Except the Landschutzer CO<sub>2</sub> product, which you will need to download yourself. A curl script is provided below). The figures are scripts 9a-f. 

##### Getting Conda working

Firstly, I assume you have a Unix environment. Clone the repository to where you want.

`git clone https://github.com/nicpittman/tropical_pacific_carbon_export.git`

Secondly, you will need to ensure you have all of the dependencies outlined in *requirements.txt*

You can create a Conda environment (with conda-forge) using *requirements.txt*:

```
conda config --append channels conda-forge
conda create --name pacific_carbon --file requirements.txt
```

Or if it breaks for a different reason, something like below can also work (Spyder not essential):

```
conda create -n pacific_carbon python=3.7 ESMPy=7.1.0 xesmf==0.3 spyder=4.0 xarray=0.15 h5py ipython Markdown numpy pandas scipy matplotlib dask nco netcdf4 statsmodels h5netcdf bs4c joblib pillow lxml basemap cartopy curl cbsyst
```

Otherwise, try removing the library versions, but this has worked on my system. Things might break if you use different versions. 

and then activate it like so:
```
conda activate pacific_carbon
```


### 	2. Reproduce figures:

I have noticed the figures turn out differently if run in terminal compared to using Spyder (likely due to graphics software differences, on the supercomputer i get the error: MESA-LOADER: failed to load driver swrast (search paths /usr/lib64/dri) - This may be fixable? I don't know. Spyder is included in the dependencies and should reproduce the figures accurately if used.  The numbers reported in the manuscript are printed to console during figure production, and some are saved into *processed/results/\*.csv*.

##### Data to download

Most of the processed data is provided in order for the plotting functions (9[a-f]) to work. However, two sets of data will need to be downloaded manually. SST and the landschutzer CO<sub>2</sub> flux product (and NCO to convert to a usable format). You can perform this either running `sh datasets/to_download.sh`, or by pasting the following code into console:

SST:

```To download manually:
mkdir datasets/sst | curl ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc --output datasets/sst/sst.mnmean.nc
```

Landschutzer CO<sub>2</sub>  (https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html)

````
mkdir -p co2/landschutzer_co2 | curl https://data.nodc.noaa.gov/ncei/ocads/data/0160558/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc --output co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc
````

You will need to use NCO to convert the time variable name (t to date) so they can be opened by xarray. Included in the conda environment.
```
ncrename -v date,t co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc 
```

TAO data is included here in datasets/tao/tao_physics/*. You can update it yourself here (but not essential): https://www.pmel.noaa.gov/tao/drupal/disdel/ using all variables at Moorings 110W,125W, 140W,155W,170E,165E

##### Data included

The following data files are included here for the reproduction of figures (All files are created during the 8. Cleanup script):

    #These files here are monthly and 1 degree averaged from oregon state university website.		
    processed/flux/avg_npp_rg_cbpm.nc 	
    processed/flux/avg_npp_rg_cafe.nc					
    processed/flux/avg_npp_rg_vgpm.nc					
    processed/flux/avg_npp_rg_eppley.nc 				
    
    #These files are produced from multiple data sources, see full pipeline for more information
    processed/combined_dataset/month_data.nc			
    processed/combined_dataset/month_data_exports.nc 	
    processed/flux/fratios.nc 							
    processed/flux/pco2grams.nc 					
    processed/flux/tpca.nc 								
    processed/earth_m2.nc							
    processed/flux/landsch_mooring_co2_flux.nc 		
    processed/flux/npp.nc								
    processed/indexes/el_nino_events.csv				
    processed/indexes/la_nina_events.csv



##### Producing the figures

After you have installed the dependencies, downloaded SST and Landschutzer CO<sub>2</sub>, in terminal, you should be able to run (or preferably in Spyder for correct figure size formatting):

```
python 9a ... (tab)
python 9b ...
python 9c ...
python 9e ...
python 9f ...
```

Note*, Figure 4 (9d) uses basemap, a depreciated package. A rough version of this figure is provided using cartopy but is not publication quality, so I recommend using the basemap version if possible. It should work in the provided Conda environment. 





### 	3. Reprocess entire pipeline:

The entire the pipeline is provided for reproducible research, however is provided as is. I have tried to make the flow as streamlined as possible but it is likely that there will be unique problems on different systems. 

###### Notes: 

- In total, after downloading and processing all required data the folder will be ~50gb., and possibly more during Script 1, when it downloads and consolidates NPP data.  
- Script 1 takes a long time to download 30gb of NPP data. Script 2 then processes this. Script 2 will download TPCA chlorophyll data automatically, and you will also need to download the NASA chlor_a data as described below (Scripts included, 2)
- Scripts 3-6 can be run (in order) at the same time as script 1. 
- Script 7a and b (parallelised version), combines everything and takes a long time to process (6 hours per mooring, or ~6 hours x 6 cores total for script 7b)- so this is a warning here. 



Two dataset types are produced. Firstly, 6 Moorings at 110W, 125W, 140W, 155W, 170W, 165E provide in-situ temperature, windspeed, thermocline depth, pCO<sub>2</sub> in water and atmosphere among other variables. These have been compiled with different NPP products and algorithms, SST and different CO<sub>2</sub> flux products. Secondly, data used in the paper uses a 1 degree, 1 month resolution, however products for all, daily and weekly are also produced throughout this process (datasets/combined_data/*nc). 

Secondly, all primary productivity models, f-ratios, co2 flux and sst products are gridded to 1 degree grids and monthly grids since September 1997. You will need xesmf for these to work properly (in cleansup script). These are used to produce Figure 4/D.

##### Scripts included

*carbon_math.py* is a series of functions that make it easy to convert between carbon units.

Scripts are organised to be run in numerical order. 

1. Downloads Primary Production data from www.science.oregonstate.edu/ocean.productivity/

   1. Besides Script 7, this is probably the slowest. Downloads ~30gb. 

2. Cuts the NPP data into mooring csv files. 

   1. **NOTE** This script is not actually essential, has been processed data provided in processed/npp_mooring_timeseries.

   2. NASA chlor_a needs to be downloaded manually. A script in 'datasets/chl/download_NASA_chlora.sh' is provided. Will need authorisation cookie as described at: https://oceancolor.gsfc.nasa.gov/data/download_methods/. Sourced from https://oceandata.sci.gsfc.nasa.gov/ 

   3. Once you create an account at earth data, set the cookie like this. 

   4. ```
      echo "machine urs.earthdata.nasa.gov login USERNAME password PASSWD" > ~/.netrc ; > ~/.urs_cookies
      chmod  0600 ~/.netrc
      ```

   5. TPCA should download automatically during script 2. They will be downloaded to datasets/chl/*nc and then moved automatically to the relevant folder /chlor_a/ seawifs or modis. Sourced from: https://researchdata.ands.org.au/tropical-pacific-chlorophyll-reprocessing-v10/1438905 

3. Downloads CO<sub>2</sub> flux for several products. (Can run scripts 3-5 simultaneous to script 1 in a separate console)

4. Downloads mooring CO<sub>2</sub> data (Including Japanese JMA product, downloaded automagically)

5. Downloads mooring physics data and combines with the mooring CO<sub>2</sub> data.

6. Is a shell file with a list of other datasets that need to be manually downloaded (or just the scripts outlined here).

   - SIMPLE-TRIM export production model from  https://tdevries.eri.ucsb.edu/models-and-data-products/ to datasets/exports

     - To download, you can use this, or run shell script 6-Data_shell_download_SIMPLETRIM_and_indicies.sh:

       ```
       cd datasets | curl https://tdevries.eri.ucsb.edu/wp-content/uploads/2018/03/SIMPLE_TRIM_output.nc --insecure --output SIMPLE_TRIM_output.nc
       ```

   - You should have already downloaded Landschutzer 2020, however 2018 is needed in cleanup script to make the seamask. Not essential as seamask.nc is iincluded in the repo, and can be commented out in script 8. Can be downloaded like 

     ```
     cd datasets/cos/landschutzer_co2 | curl https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc --output spco2_MPI_SOM-FFN_v2018.nc
     
     ncrename -v date,t spco2_MPI_SOM-FFN_v2018.nc
     ```

   - Climate indicies MEI, PDO, EMI, SOI are included in the repository. Scripts to download are included in 6, however not essential, and will need some processing to be opened by pandas.

   - Re EMI. The original data file was downloaded from: http://www.jamstec.go.jp/aplinfo/sintexf/DATA/emi.monthly.txt which has now changed to: http://www.jamstec.go.jp/virtualearth/data/SINTEX/SINTEX_EMI.csv.  This is reflected in the scripts.

8. a - Processes all of this data into the mooring timeseries. This is an inefficient script (ie. will take days to run), and I would implement this differently if I was to rewrite this. **Because it is so inefficient, a parallelised version (7b) has been produced to be run** on the Australian supercomputer GADI/NCI but still takes 6 cores and 6 hours (or can be run on 6 or 3 native cores). Oops for pandas inefficiencies. Be aware here. 

   
   
8. Is a cleanup script with miscellaneous cleanup functions that i have added on add hoc during development. For example, CAFE was released late in the development process, and a function will plug this (and SST) into the Mooring timeseries.  Another part finds ENSO events, and another function converts the mooring csvs into an easy to use netcdf file. These have been debugged but may still contain system specific problems. 

   1. You will need to change some filepaths, particularly the `convert_TPCA_to_month()`function. You will also need to create the directory datasets/tpca otherwise it may brexak (unless you have stored TPCA in datasets/tpca, which would make sense.) `mkdir datasets/tpca`
   2. The carbon conversion( `carbon_uatm_to_grams()`) takes a long time, and I get killed on the supercomputer login nodes, for I assume memory. This can be made into its own script and submitted very easily if you have this problem. Again, not entirely necessary as the relevant files are provided.

9. Are plotting scripts as discussed above.

   a. Figure 1 - Comparison of NPP, and fratios at mooring locations

   b. Integrated plot showing timeseries for CO2 and the CAFE * Laws2011a product.

   c. Seasonal decomposition. This script also includes Figure 4, ENSO and seasonality.

   d. Spatial maps for new production, CO2 flux and the difference between the two. Provides mean, trends and standard deviation.

   e. Mean seasonality and ENSO for different variables.

   f. Compares the overall volumes of CO2 removal from the equatorial Pacific (East, West, Central and overall).

1. run_combiner.sh will submit 7b to the supercomputer. You may need to change this depending on your system paths, or you may not need it at all. Depending on the number of cores, you can change the bottom script 7b to use say 3 cores, taking approximately 12 hours. Has been left at default of 6 cores, one for each mooring.

##### Other scripts:

- *carbon_math.py* is a series of functions that make it easy to convert between carbon units.

- 10-windspeed is a quick correlation assessment of the windspeed vs new production and air-sea flux, with r<sup>2</sup> used in the paper.

#### Notes

- Things have changed around a lot, with some redundant code still in the repo. I have tried to refactor where possible, but there is still some messiness and possible system dependent issues.  I have tried to fix this the best I can with Conda (see the conda section). I have tried my best to provide at least the minimum standard of code for reproducibility. 
- All data for the figures is provided. However, all scripts to process the data yourself are also provided for complete reproducability of this paper 
