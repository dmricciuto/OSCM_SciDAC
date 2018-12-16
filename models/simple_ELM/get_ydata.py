#!/usr/bin/env python

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

from utils import *
from common import oscm_dir

def get_ydata(xdata_full):

    qoi='gpp'
    twelve = 12



    obsdata = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
    sites_info = utils.read_sites(obsdata)
    nsites = len(sites_info)

    ndata_all = xdata_full.shape[0]

    sim_locs=np.unique(xdata_full[:,0:2],axis=0)
    obs_locs = np.array([sites_info[j][1:3] for j in range(nsites)])
    site_xind =  utils.pick_sites2(sim_locs,obs_locs)

    ydata_mean = np.empty((0,))
    ydata_std = np.empty((0,))
    ydata_all=np.empty((0,0))
    ind_data_all=[]

    stat_ind = 0 # mean

    for site_id, lon,lat in site_xind:

        if qoi=='gpp':
            qoi_ind=0
            var_daily = obsdata.variables['GPP'][site_id,:]*24*3600*1000

            ydata = utils.daily_to_monthly(var_daily).reshape(-1,twelve).T

            var_monthly=np.nanmean(ydata,axis=1)
            var_monthly_std=np.nanstd(ydata,axis=1)
            ydata_mean   = np.append(ydata_mean, var_monthly)
            ydata_std   = np.append(ydata_std, var_monthly_std)
            if ydata_all.shape[1]==0:
                ydata_all = ydata.copy()
            else:
                ydata_all = np.append(ydata_all, ydata,axis=0)

            ind_data=utils.pick_ind(xdata_full,[lon,lat,qoi_ind,None,stat_ind])
            ind_data_all.extend(ind_data)

    print ind_data_all
    if stat_ind==0:
        ydata=ydata_mean.copy()

    # ind_calib=np.zeros((ndata_all,))
    # ind_calib[ind_data_all]=1
    # np.savetxt('ind_calib.dat', ind_calib, fmt='%d')

    # xdata_calib = xdata_full[np.array(ind_data,dtype=int),:]
    # np.savetxt('xdata_calib.txt',  xdata_calib,  fmt='%.2f')



    np.savetxt('ydata_all_mean.txt', ydata_mean, fmt='%.12f')
    np.savetxt('ydata_all_std.txt', ydata_std, fmt='%.12f')
    np.savetxt('ydata_all_sam.txt', ydata_all, fmt='%.12f')

    return ydata, ind_data_all, ydata_all, ydata_mean, ydata_std

###################################################################
###################################################################
###################################################################

xdata = np.loadtxt('xdata_full.dat')
ydata, ind_data_all, ydata_all, ydata_mean, ydata_std = get_ydata(xdata)
np.savetxt('ind_calib.dat',ind_data_all)


# [[ 5 28 40]
#  [22 48 29]
#  [23 13 36]
#  [24 25 27]
#  [25 12 35]
#  [27 14 36]
#  [29 23 35]
#  [31 12 35]]
