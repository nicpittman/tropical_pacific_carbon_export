#REDUNDANT.
#Only SIMPLE-TRIM is required:
echo 'Downloading SIMPLE-TRIM'
cd datasets | curl https://tdevries.eri.ucsb.edu/wp-content/uploads/2018/03/SIMPLE_TRIM_output.nc --insecure --output datasets/SIMPLE_TRIM_output.nc

#Download https://tdevries.eri.ucsb.edu/models-and-data-products/ to datasets/
#Most of the files required have either been downloaded using scripts 1-5 here, or included in the repository.

#Mooring Physics:
#https://www.pmel.noaa.gov/tao/drupal/disdel/
#These were obtained by selecting every variable for single location timeseries. These files are included here by default in datasets/tao/#physics. These can be downloaded yourself if you want.


#Climate Indexes are included here. If you want to download them yourself you can:
# NOTE THAT THE FORMATING NEEDS TO BE FIXED IF YOU DOWNLOAD MANUALLY
#SOI: https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/SOI/
#MEI: https://www.esrl.noaa.gov/psd/enso/mei/
#PDO: http://research.jisao.washington.edu/pdo/PDO.latest

echo 'CLIMATE INDICIES'
echo 'Uncomment these if you need to download them, but note you must fix them in excel before being opened by pandas'

# curl https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/soi.long.data --output soi.csv
# curl https://psl.noaa.gov/enso/mei/data/meiv2.data --output meiv2.csv
# curl http://research.jisao.washington.edu/pdo/PDO.latest --output pdo.csv
# curl http://www.jamstec.go.jp/frsgc/research/d1/iod/DATA/emi.monthly.txt --output emi.csv

#SOI:
#then open in notepad and replace '  ' (two spaces) with a comma, so it is a csv. A slightly modified version is provided here.
#MEIv2:
#Open in notepad and replace '     ' (four spaces) with a comma so it is a csv.
#PDO:
#Reprocess as above. A copy is provided here though.
#EMI
#And process dates. A file of this is included albeit not used.

