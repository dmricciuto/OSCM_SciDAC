

import utils
import numpy as np
from common import oscm_dir
from netCDF4 import Dataset


obs_dataset = Dataset(oscm_dir + "/models/site_observations/fluxnet_daily_obs.nc4",
                      'r', format='NETCDF4')
sites_info = utils.read_sites(obs_dataset)
nsites = len(sites_info)
print(sites_info)
np.savetxt('site_coords.txt', [sites_info[j][1:3]
                               for j in range(nsites)], fmt='%.2f')
np.savetxt('site_names.txt', [sites_info[j][0]
                              for j in range(nsites)], fmt='%s')
