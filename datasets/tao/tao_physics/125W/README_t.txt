
To TAO/TRITON, PIRATA, and RAMA data users:

Please acknowledge the GTMBA Project Office of NOAA/PMEL if you 
use these data in publications.  Also, we would appreciate receiving 
a preprint and/or reprint of publications utilizing the data for 
inclusion in the GTMBA Project bibliography.  Relevant publications
should be sent to:

Global Tropical Moored Buoy Array Project Office
NOAA/Pacific Marine Environmental Laboratory
7600 Sand Point Way NE
Seattle, WA 98115

Comments and questions should be directed to atlasrt.pmel@noaa.gov .

Thank you.

Michael J. McPhaden
GTMBA Project Director
-------------------------------------------------------------------

The following topics will be covered in this readme file:

  1. Upper Ocean Temperatures
  2. Time Stamp
  3. 5-day, Monthly and Quarterly Averages
  4. Sampling, Sensors, and Moorings
  5. Quality Codes
  6. Source Codes
  7. Instrument Codes
  8. References

1. Upper Ocean Temperatures:

Upper ocean temperatures between the sea surface (1 m) and
500 m depth begin in 1977 in the Tropical Pacific and continue 
to the present. More recently temperatures are also available
from the Atlantic and Indian oceans. TRITON buoys replaced 
ATLAS buoys in the Pacific west of 160E beginning in 1999. 
For sites recently occupied by TRITON buoys, SST is measured 
at 1.5 m depth, and the deepest temperature is measured at 
750 meters. Temperatures are in units of degrees centigrade.

Three sites in the Indian Ocean, 1.5S,90E, 5S,95E, 
and 8S,95E have data from a TRITON system, called
mTRITON. SST depths for mTRITON data are measured
at 1 meter depth. 

In all ascii files organized by site of daily and hourly SST data 
the depth is specied at the top of each block of data. 

In netcdf files organized by site of daily and hourly SST data a 
time series of the SST depth called ZSST is included to specify 
the SST depth at all times in the file. 

One site in the Indian Ocean at 8S,100E called "BaiLong", 
was first deployed in 2010, and also measures SST at 1 meter.

In August 2011 NDBC Refresh buoys replaced some ATLAS 
buoys in the Pacific. For details on the conversion of 
ATLAS to Refresh, see

  http://www.pmel.noaa.gov/gtmba/refresh-system

In 2015 TFLEX systems were first deployed in the Atlantic 
and Indian Oceans. These buoys measure SST at 1 meter depth.

For standard depths of ATLAS temperature measurements, see

  http://www.pmel.noaa.gov/gtmba/moorings#depths

In certain instances, additional temperature sensors may
have been added to ATLAS moorings for special purposes.

At TAO mooring sites which measured subsurface currents in
addition to temperatures, i.e. on the equator at 156E, 165E,
140W, 125W, 110W, and 95W, and at 7N, 140W, there have been
a wide variety of depths for temperature sensors over time 
in response to changing scientific priorities.

For a qualitative picture of the distribution of depths 
with time, you can view data availability plots on the 
data display and delivery page

     http://www.pmel.noaa.gov/tao/drupal/disdel/

by clicking the green availability button, when you 
have selected T(z) in Time Series mode.

In addition, all depths are clearly indicated in the 
data files themselves, except that in files containing 
multiple sites, and in netcdf files containing both 
ATLAS and TRITON data, all SST's are given the ATLAS 
depth of 1 meter.

If you selected daily data at 8n137e or 6s10w, you may
get more than one file per site. This is because the 
original deployments for these sites were at 7n137e and
5s10w, respectively, while present deployments are at 
8n136e and 6s10w. The file names will clearly indicate
which site the data come from. For more details about 
the mooring locations, you can deliver daily average
position data for most deployments under the data category
of "Buoy Positions" on the delivery page. 

If you selected the site at 0n110w, you may get data in 
separate files from several groups of deployments clustered 
around 0n110w since mooring locations have been in significantly 
different locations at different times. The file names will 
clearly indicate the site locations. See this page for more
information

  http://www.pmel.noaa.gov/tao/drupal/disdel/110w.html

If you selected high resolution data, you may find that 
you have several files, each with a different averaging
interval, for example, hourly, 15 minute, and 10 minute.
The interval in each file is indicated by one of the 
following file name suffixes: "hr" for hourly, "08h" 
for 8 hour, "06h" for 6 hour, "02h" for 2 hour, "04h" 
for 4 hour, "30m" for 30-minute, and similarly for 
"15m", "10m", "7.5m", 3.75m", and "01m".  

Also, some files may have a _tr_ string the name, which 
indicates spot samples rather than averages.

2. Time Stamp:

Time associated with data represent the sample time for single sample
values or the middle of the averaging interval for average values.  For
example, daily averages are computed starting at 0000 GMT and are
assigned an observation "time stamp" of 1200 GMT.

3. 5-day, Monthly and Quarterly Averages:

If you delivered 5-day, monthly, or quarterly averaged data
these definitions are relevant to your files:

5-Day: Average of data collected during consecutive five day 
intervals. A minimum of 2 daily values are required to compute 
a 5-day average.

Monthly: Average of all the data collected during each month.
A minimum of 15 daily values are required to compute a monthly 
average.

Quarterly: Average of 3 monthly values. A minimum of 2 monthly 
values are required to compute a quarterly average. 12 quarterly 
averages are computed for each year, one for each center month, 
which includes the previous month, the center month, and the 
next month in the average.

4. Sampling, Sensors, and Moorings:

For detailed information about sampling, sensors, and moorings,
see these web pages:

  http://www.pmel.noaa.gov/gtmba/sensor-specifications

  http://www.pmel.noaa.gov/gtmba/sampling

  http://www.pmel.noaa.gov/gtmba/moorings


5. Quality Codes:

  0 = datum missing

  1 = highest quality; Pre/post-deployment calibrations agree to within
  sensor specifications.  In most cases only pre-deployment calibrations 
  have been applied

  2 = default quality; Pre-deployment calibrations applied.  Default
  value for sensors presently deployed and for sensors which were either 
  not recovered or not calibratable when recovered.

  3 = adjusted data; Pre/post calibrations differ, or original data do
  not agree with other data sources (e.g., other in situ data or 
  climatology), or original data are noisy.  Data have been adjusted in 
  an attempt to reduce the error.

  4 = lower quality; Pre/post calibrations differ, or data do not agree
  with other data sources (e.g., other in situ data or climatology), or 
  data are noisy.  Data could not be confidently adjusted to correct 
  for error.

  5 = sensor or tube failed

  C (ascii) or -9 (netcdf) = Indicates special adjustments were 
  made to the data. For further information, see the following:

 	Freitag, H.P., M.E. McCarty, C. Nosse, R. Lukas, M.J. McPhaden,
 	and M.F. Cronin, 1999: COARE Seacat data: Calibrations and
 	quality control procedures. NOAA Tech. Memo. ERL PMEL-115,
 	89 pp.

6. Source Codes:

  0 - No Sensor, No Data 
  1 - Real Time (Telemetered Mode)
  2 - Derived from Real Time
  3 - Temporally Interpolated from Real Time
  4 - Source Code Inactive at Present
  5 - Recovered from Instrument RAM (Delayed Mode)
  6 - Derived from RAM
  7 - Temporally Interpolated from RAM

7. Instrument Codes:

   0 - No Sensor
   4 - Conductivity (FSI)
  14 - NextGen Conductivity
  24 - NextGen Conductivity (Firmware version 5.03+)
  70 - Seacat Conductivity
  71 - Microcat Conductivity   
  99 - Unknown

8. References:

For more information about TAO/TRITION, PIRATA, and RAMA, see

McPhaden, M.J., A.J. Busalacchi, R. Cheney, J.R. Donguy,K.S. 
Gage, D. Halpern, M. Ji, P. Julian, G. Meyers, G.T. Mitchum, 
P.P. Niiler, J. Picaut, R.W. Reynolds, N. Smith, K. Takeuchi, 
1998: The Tropical Ocean-Global Atmosphere (TOGA) observing 
system:  A decade of progress. J. Geophys. Res., 103, 14,
169-14,240.

Bourles, B., R. Lumpkin, M.J. McPhaden, F. Hernandez, P. Nobre, 
E.Campos, L. Yu, S. Planton, A. Busalacchi, A.D. Moura, J. 
Servain, and J. Trotte, 2008: The PIRATA Program: History, 
Accomplishments, and Future Directions. Bull. Amer. Meteor. 
Soc., 89, 1111-1125.

McPhaden, M.J., G. Meyers, K. Ando, Y. Masumoto, V.S.N. Murty, M.
Ravichandran, F. Syamsudin, J. Vialard, L. Yu, and W. Yu, 2009: RAMA: The
Research Moored Array for African-Asian-Australian Monsoon Analysis and
Prediction. Bull. Am. Meteorol. Soc., 90, 459-480,
doi:10.1175/2008BAMS2608.1
