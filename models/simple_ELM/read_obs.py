from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os

from utils import *

from common import oscm_dir




obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
sites_info = read_sites(obs_dataset)
nsites = len(sites_info)

#np.savetxt('site_lons.txt',site_lons, fmt='%.2f')
np.savetxt('site_coords.txt',[sites_info[j][1:3] for j in range(nsites)], fmt='%.2f')
np.savetxt('site_names.txt',[sites_info[j][0] for j in range(nsites)],fmt='%s')
