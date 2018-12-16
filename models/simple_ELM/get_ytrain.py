from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

import utils
from common import oscm_dir



ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)
dataset = Dataset(ncfile)

xdata, outnames, ytrain = utils.read_simdata_ytrain(dataset)

np.savetxt('xdata_all.txt',  xdata,  fmt='%.2f')
np.savetxt('ytrain_all.dat', ytrain, fmt='%.12f')
np.savetxt('outnames_all.txt',outnames,fmt='%s')


