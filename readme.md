### Tropical Pacific CO2 flux and NPP rates
This set of python scripts provides all data and figures for Pittman 2020; Physical drivers of air-sea CO2 flux and biological draw-down in the equatorial Pacific.

Some things are rough, apologies for the mess. any problems with implementation email me at nic.pittman [at] utas[dot]edu[dot]au

Two dataset types are produced. Firstly, 6 Moorings at 110W, 125W, 140W, 155W, 170W, 165E provide in-situ temperature, windspeed, thermocline depth, pCO2 in water and atmosphere among other variables. These have been compiled with different NPP products and algorithms, SST and different CO2 flux products. This data used in the paper uses a 1 degree, 1 month resolution, however products for all, daily and weekly are also produced throughout this process. 

Secondly, all primary productivity models, f-ratios, co2 flux and sst products are gridded to 1 degree grids and monthly grids since September 1997. These can be used to plot, make hovmollers and calculate longterm trends. 

Some debugging may be required to transition this to your system. Most file paths should work, 

#### Files

Scripts are organised to be run in numerical order. 

1. Downloads Primary Production data from www.science.oregonstate.edu/ocean.productivity/

2. Cuts the NPP data into mooring csv files

3. Downloads CO2 flux for several products.

4. Downloads mooring CO2 data

5. Downloads mooring physics data and combines with the mooring CO2 data.

6. Is a text file with a list of other datasets that need to be manually downloaded.

   Chl needs to be downloaded manually from 

   - https://oceandata.sci.gsfc.nasa.gov/
   - TPCA: https://researchdata.ands.org.au/tropical-pacific-chlorophyll-reprocessing-v10/1438905

7. a - Processes all of this data into the mooring timeseries. This is inefficient (will take days to run), and I would implement this differently if I was to rewrite this. Because it is so inefficient, a parallelised version has been produced to be run on the Australian supercomputer GADI but still 6 cores and 6 hours. Oops for pandas inefficiencies. Easier to plug and play the data provided here (see below). 

   1. Run_combiner.sh will submit 7b to the supercomputer.

8. Is a cleanup script with miscellaneous cleanup functions that i have added on add hoc during development. For example, CAFE was released late in the development process, and a function will plug this (and SST) into the Mooring timeseries.  Also finds ENSO events and another function converts the mooring csvs into an easy to use netcdf file. There are a number of files under 8_ as well. These all need to get run. They are a bit messy

9. Are the plotting scripts

   a. Figure 1 - Comparison of NPP, and fratios at mooring locations

   b. Integrated plot showing timeseries for CO2 and the CAFE * Laws2011a product.

   c. Seasonal decomposition. This script also includes Figure 4, ENSO and seasonality.

   d. Spatial maps for new production, CO2 flux and the difference between the two. Provides mean, trends and standard deviation.

   e. Mean seasonality and ENSO for different variables.

   f. Compares the overall volumes of CO2 removal from the equatorial Pacific (East, West, Central and overall).

Other scripts:

- Carbon_match is a series of functions that make it easy to convert between mol and gC



#### Notes

- Things have changed around a lot, with some redundant code. I have tried to refractor where possible, but still some messiness. 

- All data for the figures is provided. However, all scripts to process the data yourself are provided. 

- All data should go into datasets/ This is then recreated in processed/ as the files we want.

  

### To reproduce figures, some data is provided in this repo.

Some of the processed data is provided in order for the plotting functions (9.x[a-f]) to work. However, three sets of data will need to be downloaded manually.

```To download manually:
To be downloaded manually:
datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc - This will need to be downloaded manually, see source file in folder.
	https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html

datasets/sst/sst.mnmean.nc - Downloaded from, see source file in folder.
	https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html
	
datasets/tao/tao_physics/* - Downloaded from, see source file in folder.
	https://www.pmel.noaa.gov/tao/drupal/disdel/
```

The files in this repo include:

    Included:
    
    #These files here are monthly and 1 degree averaged from oregon state.
    datasets/npp_satellite/avg_npp_rg_cbpm.nc 			- Created in cleanup
    datasets/npp_satellite/avg_npp_rg_cafe.nc			- Created in cleanup
    datasets/npp_satellite/avg_npp_rg_vgpm.nc			- Created in cleanup
    datasets/npp_satellite/avg_npp_rg_eppley.nc 		- Created in cleanup
    
    #These files are produced from multiple sources
    processed/combined_dataset/month_data.nc			- Created in cleanup
    processed/combined_dataset/month_data_exports.nc 	- Created in cleanup
    
    processed/flux/fratios.nc 							- Created in cleanup
    datasets/co2/pco2grams_eq1.nc 						- Created in cleanup
    datasets/tpca/tpca.nc 								- Created in cleanup
    processed/earth_m2.nc								- Created in cleanup
    processed/flux/landsch_mooring_co2_flux.nc 			- Created in cleanup
    processed/flux/npp.nc								- Created in cleanup
    processed/indexes/el_nino_events.csv
    processed/indexes/la_nina_events.csv