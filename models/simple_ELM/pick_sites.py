#!/usr/bin/env python

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

from utils import *
from common import oscm_dir



ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)
dataset = Dataset(ncfile)
qoi='gpp'
twelve = 12
nyears = 24

print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lats = dataset.variables['lat'][:] #.shape
lons = dataset.variables['lon'][:]-360. #.shape

obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
site_names, site_lons, site_lats = read_obsdata(obs_dataset)
site_lon_lat =  pick_sites(lons,lats,site_lons,site_lats)

nsites=len(site_lon_lat)

xdata_full = np.loadtxt('xdata_full.txt')
ind_surr = []
ydata_mean = np.empty((0,))
ydata_std = np.empty((0,))
ydata_all = np.empty((0,nyears))
for site_id,lon_id,lat_id in site_lon_lat:
    print lon_id, lat_id

    inds = []
    inds.append(np.where(xdata_full[:,0]==lon_id)[0])
    inds.append(np.where(xdata_full[:,1]==lat_id)[0])
    inds.append(np.where(xdata_full[:,2]==0.0)[0])
    inds.append(np.where(xdata_full[:,4]==0.0)[0])

    ind_out = np.arange(xdata_full.shape[0])
    for ind in inds:
        ind_out = np.intersect1d(ind_out,ind)
    ind_surr = np.append(ind_surr,ind_out,axis=0)

    var_daily = obs_dataset.variables[qoi.upper()][site_id,:]*24*3600*1000
    ydata = daily_to_monthly(var_daily).reshape(nyears,twelve).T

    var_monthly=np.nanmean(ydata,axis=1)
    var_monthly_std=np.nanstd(ydata,axis=1)
    ydata_mean   = np.append(ydata_mean, var_monthly)
    ydata_std   = np.append(ydata_std, var_monthly_std)
    ydata_all = np.append(ydata_all, ydata,axis=0)


np.savetxt('ind_surr.dat', ind_surr, fmt='%d')

xdata_surr = xdata_full[np.array(ind_surr,dtype=int),:]
np.savetxt('xdata_surr.txt',  xdata_surr,  fmt='%.2f')



np.savetxt('ydata_all_mean.txt', ydata_mean, fmt='%.12f')
np.savetxt('ydata_all_std.txt', ydata_std, fmt='%.12f')
np.savetxt('ydata_all_sam.txt', ydata_all, fmt='%.12f')





# [[ 5 28 40]
#  [22 48 29]
#  [23 13 36]
#  [24 25 27]
#  [25 12 35]
#  [27 14 36]
#  [29 23 35]
#  [31 12 35]]
