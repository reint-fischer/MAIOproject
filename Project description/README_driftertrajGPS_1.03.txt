The archive driftertrajGPS_1.03.tar.gz is a gzip-compressed tar archive that contains individual files named dX, one per GPS drifter trajectory with the Global Drifter program ID X.

Each file contains hourly location and velocity of Iridium or Argos GPS-tracked surface drifters from the Global Drifter Program. The method of interpolation is described in Elipot et al. 2016, "A global surface drifter dataset at hourly resolution", JGR-Oceans, doi:10.1002/2016JC011716.  For details on what is new with this version, see http://www.aoml.noaa.gov/phod/gdp/hourly_data.php.

Each file is an ASCII text file with contiguous hourly estimates, one set of hourly estimates per line, unless the interpolating temporal gap between orginal GPS positions is longer than 12 hour. In that case, no estimate is provided. When analyzing these data as time series, be aware of temporal time steps that are larger and multiple of 1 hour. You may want to consider the netcdf file for these data that indicates interruptions in contiguous hourly estimates by "Infinity" values; and separation of individual trajectories by "NaN" values.

1 : identification number from the Global Drifter Program (see deployment log file at http://www.aoml.noaa.gov/phod/dac/deployed.html)
2 : year
3 : month
4 : day
5 : hour of the day (UTC, 0-23)
6 : longitude (degree East)
7 : latitude (degree North)
8 : longitude error (10^-5 degree, missing value: Inf)
9 : latitude error (10^-5 degree, missing value: Inf)
10 : zonal velocity (meter per second)
11 : meridional velocity (meter per second, missing value: Inf)
12 : zonal velocity error (meter per second, missing value: Inf)
13 : meridional velocity error (meter per second)
14 : length of interpolating gap in hours; it is the difference between the time of the first GPS fix forward in time minus the time of the last GPS fix backward in time
15 :  drogue status; 1 = drogue on, 0 = drogue off

Questions about the dataset should be addressed to either Rick.Lumpkin@noaa.gov or Shane Elipot at selipot@rsmas.miami.edu
