# As described in readme.md, for figures to run, the following must be downloaded manually:


echo 'Downloading sst'
mkdir sst | curl ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2/sst.mnmean.nc --output sst/sst.mnmean.nc

echo 'Downloading CMAP precipitation'
curl ftp://ftp.cdc.noaa.gov/Datasets/cmap/enh/precip.mon.mean.nc --output datasets/precip.mon.mean.nc


mkdir co2 | mkdir co2/landschutzer_co2 | curl https://data.nodc.noaa.gov/ncei/ocads/data/0160558/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc --output co2/landschutzer_co2/spco2_MPI-SOM_FFN_v2020.nc
cd co2/landschutzer_co2
ncrename -v t,date spco2_MPI-SOM_FFN_v2020.nc 

echo 'You might also want the 2018 version for full pipeline, uncomment if needed.'
# And if you want the 2018 version.
# curl https://www.nodc.noaa.gov/archive/arc0105/0160558/4.4/data/0-data/MPI_SOM-FFN_v2018/spco2_MPI_SOM-FFN_v2018.nc --output  landschutzer_co2/spco2_MPI_SOM-FFN_v2018.nc
# ncrename -v date,t spco2_MPI_SOM-FFN_v2018.nc 
