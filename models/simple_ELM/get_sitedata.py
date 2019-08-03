
import utils
import numpy as np
from netCDF4 import Dataset

from common import oscm_dir

qoi = 'GPP'
twelve = 12


obsdata = Dataset(
    oscm_dir + "/models/site_observations/fluxnet_daily_obs.nc4", 'r', format='NETCDF4')

sites_info = utils.read_sites(obsdata)
nsites = len(sites_info)
sitenames = [sites_info[j][0] for j in range(nsites)]

ydata_mean = np.empty((0,))
ydata_std = np.empty((0,))
ydata_all = np.empty((0, 0))

j = 0
for site in sitenames:

    var_daily = obsdata.variables[qoi][j, :] * 24 * 3600 * 1000

    ydata = utils.daily_to_monthly(var_daily).reshape(-1, twelve).T

    var_monthly = np.nanmean(ydata, axis=1)
    var_monthly_std = np.nanstd(ydata, axis=1)
    ydata_mean = np.append(ydata_mean, var_monthly)
    ydata_std = np.append(ydata_std, var_monthly_std)
    if ydata_all.shape[1] == 0:
        ydata_all = ydata.copy()
    else:
        ydata_all = np.append(ydata_all, ydata, axis=0)

    j += 1

np.savetxt('ydata_all_mean.txt', ydata_mean, fmt='%.12f')
np.savetxt('ydata_all_std.txt', ydata_std, fmt='%.12f')
np.savetxt('ydata_all_sam.txt', ydata_all, fmt='%.12f')
np.savetxt('sitenames_data.dat', sitenames, fmt='%s')
