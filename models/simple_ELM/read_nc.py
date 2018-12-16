
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy
import utils

from common import pmin, pmax

ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
#qoi='gpp'

pnames, prange, ptrain, qtrain = utils.read_simdata_input(dataset)

np.savetxt('ptrain_all.dat',ptrain)
np.savetxt('pnames.txt',pnames,fmt='%s')
np.savetxt('prange.dat',prange)
