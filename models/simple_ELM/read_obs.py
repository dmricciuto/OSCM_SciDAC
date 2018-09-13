from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os

from utils import *

from common import oscm_dir




obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
site_names, site_lons, site_lats = read_obsdata(obs_dataset)

np.savetxt('site_lons.txt',site_lons)
np.savetxt('site_lats.txt',site_lats)
np.savetxt('site_names.txt',site_names,fmt='%s')
