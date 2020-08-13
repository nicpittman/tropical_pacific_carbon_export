
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

  1. Fixed Depth Currents
  2. Time Stamp
  3. 5-day, Monthly and Quarterly Averages
  4. Sampling, Sensors, and Moorings
  5. Quality Codes
  6. Source Codes
  7. Instrument Codes
  8. References

1. Fixed Depth Currents:

Currents use the oceanographic convention, so that a current with 
zonal and meridional components of 1 and 1, is flowing toward 
the northeast. Components and speed have units of centimeters
per second, and direction has units of degrees, clockwise from 
true north. For example, current with a direction of 90 degrees
is directed to the east.

Fixed depth current data begin in 1977 and continue 
to the present at sites located across the Tropical 
Pacific, and at a few sites in the Atlantic and Indian
oceans. Fixed depth current measurents prior to 2000 
were made with mechanical current meters, ie., either 
Vector Averaging Current Meters (VACM) which employed 
a Savonious rotor and vane, or Vector Measuring Current 
Meters (VMCM) which employed orthogonal rotors. After 2000, 
measurements were made with Sontek Argonaut single-point, 
acoustic-doppler current meters.

If you selected high resolution currents, you may find that 
you have several files, each with a different recording 
interval. The interval is indicated in each file name by 
one of the following suffixes: "hr" for hourly, "02h" for 
two hours, "20m" for 20-minute, and similarly, "30m", "15m", 
"10m", or "7.5m".

Mechanical current meter data (i.e., prior to 2000) were sampled
thoughout the entire recording interval. Time stamps associated 
with these data are the beginning of the sampling period. Acoustic 
doppler data (i.e., 2000 and later) are means computed over periods 
of 2 to 3 minutes at typical recording intervals of 10 or 20 minutes.  
Time stamps associated with these data are the mid-point of the 
averaging interval.

If you selected the site at 0n110w, you may get data in 
separate files from several groups of deployments clustered 
around 0n110w since mooring locations have been in significantly 
different locations at different times. The file names will 
clearly indicate the site locations. See this page for more
information

  http://www.pmel.noaa.gov/tao/drupal/disdel/110w.html

In ascii format files to the right of the data you will 
find data quality and source codes which use the definitions 
below. Qualities for some inactive sites are associated with
zonal and meridional components, while for all other sites,
qualities are associated with speed and direction. The column 
headers will identify which is which: UVM means zonal and 
meridional qualities and source mode, while SDM means speed 
and direction qualities and source mode.

In NetCDF format files which contain only one site per
file, you will find quality and source variables with 
the same shape as the data, and the variable attributes 
indicate which data the qualities are associated with. 
Netcdf files which contain more than one site ("one big 
file") and which also contain inactive sites, will have a 
mixture of speed/direction qualties, and zonal/meridional
component qualities, and no indication of this in the 
file header. All Active sites have qualities associated 
with speed and direcion only.

2. Time Stamp:

Time associated with data represent the sample time for single sample
values or the middle of the averaging interval for average values.  For
example, daily averages are computed starting at 0000 GMT and are
assigned an observation "time stamp" of 1200 GMT.

Daily mean data are computed from all available high resolution 
data over a 24-hour interval beginning at 0000 GMT. Time stamps 
associated with these data are the middle of the interval 
(1200 GMT).

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
  value for sensors presently deployed and for sensors which were 
  either not recovered or not calibratable when recovered.

  3 = adjusted data; Pre/post calibrations differ, or original data do
  not agree with other data sources (e.g., other in situ data or 
  climatology), or original data are noisy.  Data have been adjusted 
  in an attempt to reduce the error.

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

7. Instrument Codes: 

  40 - Sontek
  48 - VACM/VMCM Current Meter
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
