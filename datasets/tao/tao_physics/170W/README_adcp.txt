
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

  1. Acoustic Doppler Current Profiles
  2. Time Stamp
  3. 5-day, Monthly and Quarterly Averages
  4. Sampling, Sensors, and Moorings
  5. Quality Codes
  6. Source Codes
  7. References

1. Acoustic Doppler Current Profiles:

ADCP currents use the oceanographic convention, so that a 
current with zonal and meridional components of 1.0 and 1.0, 
is flowing toward the northeast. Components and speed have
units of centimeters per second while direction has units
of degrees, clockwise from true north. For example, current 
with a direction of 90 degrees is directed to the east.

ADCP data begin in 1988 and continue to the present at 
sites located across the Tropical Pacific, and at a few
sites in the Atlantic and Indian oceans. 

Recent ADCP data from 0,147E in the West Pacific and 0,90E
in the Indian Ocean are provided by the TRITON project of 
JAMSTEC in Japan. ADCP data from 8S,100E in the Indian ocean 
are provided by FIO in China. 

2. Time Stamp:

In general, the time associated with data represent the sample 
time for single sample values or the middle of the averaging 
interval for average values. For example, daily averages are 
computed from hourly data starting at 0000 GMT and are assigned 
an observation "time stamp" of 1200 GMT.

Recent hourly ADCP data provided by GTMBA and NDBC (all sites sites 
except 0,147E, 0,90E, and 8.5S.106.75E) are true hourly averages. 
Regular samples are taken between 30 minutes before the time stamp, 
and 30 minutes after the time stamp and averaged with equal weight. 
Earlier TAO ADCP hourly data consisted of ensemble averages from 
multiple samples starting at the time stamp and lasting for 15 minutes, 
basically similar to how hourly data are computed for the JAMSTEC 
and FIO systems. More details on the TAO ADCP samping are available 
on this page

  http://www.pmel.noaa.gov/gtmba/sampling

JAMSTEC and FIO hourly data are averages of multiple samples 
starting at the time stamp, and lasting for a few minutes. 
The samples are averaged with equal weight. JAMSTEC uses 
27 samples in 3 minutes after the top of the hour. FIO uses 
200 samples between 0000 and 0015. 


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
