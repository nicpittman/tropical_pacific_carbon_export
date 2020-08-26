# As described in readme.md, for figures to run, the following must be downloaded manually:

# datasets/co2/landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc - This will need to be downloaded manually
#	https://www.nodc.noaa.gov/ocads/oceans/SPCO2_1982_present_ETH_SOM_FFN.html

mkdir co2/landschutzer_co2 | curl https://data.nodc.noaa.gov/ncei/ocads/data/0160558/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc --output landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc
cd co2/landschutzer_co2
ncrename -v date,t spco2_MPI-SOM_FFN_v2020.nc 

echo 'You might also want the 2018 version for full pipeline, uncomment if needed.'
# And if you want the 2018 version.
# curl https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc --output  landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc
# ncrename -v date,t spco2_MPI_SOM-FFN_v2018.nc 


#datasets/sst/sst.mnmean.nc 
#	https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html

mkdir sst | curl ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc --output sst/sst.mean.nc


#The following data was originally not included but has been now for reproducability sake.
datasets/tao/tao_physics/* 
#	https://www.pmel.noaa.gov/tao/drupal/disdel/
#	Moorings 110W,125W, 140W,155W,170E,165E and all physical variables.

