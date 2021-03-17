(curl -d "search=A*L3m_MO_ZLEE_Zeu_lee_9km.nc" https://oceandata.sci.gsfc.nasa.gov/api/file_search |grep getfile| cut -d '"' -f 2)> filelist.txt
(curl -d "search=S*L3m_MO_ZLEE_Zeu_lee_9km.nc" https://oceandata.sci.gsfc.nasa.gov/api/file_search |grep getfile| cut -d '"' -f 2 )>> filelist.txt


while IFS='' read -r LINE || [ -n "${LINE}" ]; do
    echo "processing line: ${LINE}"
    curl -b ~/.urs_cookies -c ~/.urs_cookies -L -n --retry 10 -O $LINE;
done < filelist.txt

