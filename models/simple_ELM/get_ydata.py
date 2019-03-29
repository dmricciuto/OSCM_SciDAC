#!/usr/bin/env python

import numpy as np

from utils import get_ydata
from common import oscm_dir


###################################################################
###################################################################
###################################################################


#xdata = np.loadtxt('xdata_full.dat')
xdata = np.loadtxt('xdata_all.txt')

# ydata, ind_data_all, ydata_all, ydata_mean, ydata_std = get_ydata(xdata)
# np.savetxt('ind_calib.dat',ind_data_all)


sitenames_select, site_xind, ydata, ind_data_all, ydata_all, ydata_mean, ydata_std = get_ydata(
    xdata)

np.savetxt('ydata_all_mean.txt', ydata_mean, fmt='%.12f')
np.savetxt('ydata_all_std.txt', ydata_std, fmt='%.12f')
np.savetxt('ydata_all_sam.txt', ydata_all, fmt='%.12f')
np.savetxt('sitenames_select.dat', sitenames_select, fmt='%s')
np.savetxt('sitecoords_select.dat', site_xind[:, 1:], fmt='%.2f')
np.savetxt('ind_calib.dat', ind_data_all, fmt='%d')

# [[ 5 28 40]
#  [22 48 29]
#  [23 13 36]
#  [24 25 27]
#  [25 12 35]
#  [27 14 36]
#  [29 23 35]
#  [31 12 35]]
