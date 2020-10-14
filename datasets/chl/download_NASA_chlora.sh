# Follow NASA instructions to make the authorisation cookie
# https://oceancolor.gsfc.nasa.gov/data/download_methods/

# echo "machine urs.earthdata.nasa.gov login USERNAME password PASSWD" > ~/.netrc ; > ~/.urs_cookies
# chmod  0600 ~/.netrc

#Download Seawifs file urls to seawifs.txt
wget -q --post-data="sensor=seawifs&sdate=1997-09-01&edate=2010-12-01&dtype=L3M&addurl=1&results_as_file=1&search=S*DAY_CHL_chlor_a_9km*nc" -O - https://oceandata.sci.gsfc.nasa.gov/api/file_search > seawifs.txt

#Download MODIS file urls to modis.txt
wget -q --post-data="sensor=modis_aqua&sdate=2002-01-01&edate=2020-12-01&dtype=L3M&addurl=1&results_as_file=1&search=A*DAY_CHL_chlor_a_9km*nc" -O - https://oceandata.sci.gsfc.nasa.gov/api/file_search > modis.txt

#Make our dir 
mkdir chlor_a
mkdir chlor_a/seawifs
mkdir chlor_a/modis

echo 'Downloading Seawifs'
for url in $(cat seawifs.txt); do
    curl -O -b  ~/.urs_cookies -c ~/.urs_cookies -L -n $url
done
mv S*nc chlor_a/seawifs/

echo 'Downloading Modis'
for url in $(cat modis.txt); do
    curl -O -b  ~/.urs_cookies -c ~/.urs_cookies -L -n $url
done
mv A*nc chlor_a/modis/

