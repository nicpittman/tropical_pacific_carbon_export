
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

  1. Surface Meteorological Data
  2. Time Stamp
  3. 5-day, Monthly and Quarterly Averages
  4. Sampling, Sensors, and Moorings
  5. Quality Codes
  6. Source Codes
  7. References

1. Surface Meteorological Data:

Surface met data are delivered in these files as a package 
for the period after 1989 when relative humidity first 
became available. Sea surface temperature, air temperature, 
and winds are available before 1989 at several sites, and since 
1980 at 0,110W. Winds use the oceanographic convention. A wind 
with zonal and meridional components of 1.0 and 1.0 is blowing 
toward the Northeast. Daily averaged wind speeds and directions
are based on daily averaged wind velocity components. Wind speed 
and wind components are in units of meters per second, direction 
is in degrees clockwise from true north, air temperature and sea 
surface temperature are in degrees centigrade, and relative 
humidity is in percent.

ATLAS buoys measure meteorological and oceanographic data at the 
following heights and depths relative to mean sea level

     Winds: 4.0 meters height
     Relative Humidity: 3.0 meters height
     Air Temperature: 3.0 meters height
     Sea Surface Temperature: 1.0 meters depth

TRITON buoys replaced ATLAS buoys in the Pacific west of
160E beginning in 1999.  TRITON buoys measure meteorological
and oceanographic data at the following heights and depths
relative to mean sea level

     Winds: 3.5 meters height
     Relative Humidity: 2.2 meters height
     Air Temperature: 2.2 meters height
     Sea Surface Temperature: 1.5 meters depth

Three sites in the Indian Ocean, 1.5S,90E, 5S,95E, and 8S,95E 
have data from a TRITON system, called mTRITON. mTRITON buoys
measure  meteorological and oceanographic data at the following 
heights and depths relative to mean sea level

     Winds: 3.13 meters height
     Relative Humidity: 2.26 meters height
     Air Temperature: 2.26 meters height
     Sea Surface Temperature: 1.0 meters depth

In ascii files organized by site of daily and hourly data these 
depth and heights are specified at the top of each block of data. 

In netcdf files organized by site of daily and hourly data a
time series of the SST depth called ZSST is included to specify 
the SST depth at all times in the file.

Note that instrument heights are shown as negative depths 
in data files.

One site in the Indian Ocean at 8S,100E called "BaiLong", 
was first deployed in 2010, and also measures SST at 1 meter.

In August 2011 NDBC Refresh buoys replaced some ATLAS 
buoys in the Pacific. For details on the conversion of 
ATLAS to Refresh, see

  http://www.pmel.noaa.gov/gtmba/refresh-system

In 2015 TFLEX systems were first deployed in the Atlantic 
and Indian Oceans. These buoys measure SST at 1 meter depth.

The study by Freitag et al ("Calibration procedures and instrumental 
accuracies for ATLAS wind measurements", NOAA. Tech. Memo. OAR 
PMEL-119, 2001) discovered a systematic error in standard and 
NextGeneration ATLAS wind directions of approximately 6.8 degrees 
in the counterclockwise direction. This error was present possibly 
as far back as 1984. Modifications were made to the NextGeneration
ATLAS system in 2000 to correct this error in subsequent deployments, 
and archived NextGeneration ATLAS wind directions were corrected 
(both daily averages and high resolution datasets) on 28 March 2002. 
See

  http://www.pmel.noaa.gov/gtmba/vane-correction

Standard ATLAS wind directions have not been corrected in the archives 
since the exact time when the error began to affect the measurements 
is unknown. Standard ATLAS were used exclusively between 1984 and 1996 
when NextGeneration ATLAS moorings began to replace them. By November 
2001, the standard ATLAS had been phased out and the array was comprised
entirely of NextGeneration systems. Expected RMS error for standard 
ATLAS wind direction is 7.8 degrees (of which 6.8 degrees is a bias) 
while expected RMS error for NextGeneration ATLAS wind directions is 
about +/- 5 degrees with no appreciable bias.

If you selected daily data at 8n137e or 6s10w, you may
get more than one file per site. This is because the 
original deployments for these sites were at 7n137e and
5s10w, respectively, while present deployments are at 
8n136e and 6s10w. The file names will clearly indicate
which site the data come from. For more details about 
mooring locations, you can deliver daily average position 
data for most deployments under the data category of 
"Buoy Positions" on the delivery page.

If you selected high resolution data, you may find that 
you have several files, each with a different averaging
interval, for example, hourly and 10 minute. The interval 
is indicated by the filename suffixes "hr" or "10m". 

In addition, if you selected the site at 0n110w, you may
get data in separate files from several groups of deployments 
clustered around 0n110w since mooring locations have been 
in significantly different locations at different times.
The file names will clearly indicate the site locations.

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

In ascii format files organized by site, you will find data 
quality and source codes to the right of the data. In NetCDF 
format files organized by site, you will find quality and 
source variables with the same shape as the data.
These codes are defined below.

Using the quality codes you can tune your analysis to 
trade-off between quality and temporal/spatial coverage.
Quality code definitions are listed below

  0 = datum missing

  1 = highest quality; Pre/post-deployment calibrations agree to within
  sensor specifications.  In most cases only pre-deployment calibrations 
  have been applied

  2 = default quality; Pre-deployment calibrations applied.  Default
  value for sensors presently deployed and for sensors which were either 
  not recovered or not calibratable when recovered.

  3 = adjusted data; Pre/post calibrations differ, or original data do
  not agree with other data sources (e.g., other in situ data or 
  climatology), or original data are noisy. Data have been adjusted in 
  an attempt to reduce the error.

  4 = lower quality; Pre/post calibrations differ, or data do not agree
  with other data sources (e.g., other in situ data or climatology), or 
  data are noisy.  Data could not be confidently adjusted to correct 
  for error.

  5 = sensor or tube failed

6. Source Codes:

  0 - No Sensor, No Data 
  1 - Real Time (Telemetered Mode)
  2 - Derived from Real Time
  3 - Temporally Interpolated from Real Time
  4 - Source Code Inactive at Present
  5 - Recovered from Instrument RAM (Delayed Mode)
  6 - Derived from RAM
  7 - Temporally Interpolated from RAM

7. References:

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
