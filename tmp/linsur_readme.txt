pcf_all.dat: Nx48 matrix of coefficients of linear fit of each of the N outputs
           : There are total of 47 parameters, hence 48 columns - the first column is the intercept, and the rest are slopes)

xdata_all.txt : Nx5 matrix specifying latitude index (0 to 40), longitude index (0 to 60), qoi (always 0, since we have only considered GPP), month (1 to 12), and 0/1 indicator (0 is mean, 1 is standard deviation - remember the qois are actual month-averages across many years)

outnames_all.txt : is redundant, and contains same information as xdata_all.txt, but more human readable

lons.txt and lats.txt are 61 longitudes and 41 latitudes


N=38400=12 months * 2 (mean and stdev) * 1600 (locations)

Note that 1600<41*61 because not all (lat, lon) pairs are considered, only those that are in land, not in ocean.