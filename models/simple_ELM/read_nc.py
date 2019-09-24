
import sys
import utils
import numpy as np
from netCDF4 import Dataset


ncfile = sys.argv[1]
# ncfile='regional_output.nc'
print("Reading " + ncfile)

simdata = Dataset(ncfile)

pnames, prange, ptrain, qtrain = utils.read_simdata_input(simdata)

np.savetxt('ptrain_all.dat', ptrain)
np.savetxt('pnames.txt', pnames, fmt='%s')
np.savetxt('prange.dat', prange)

xdata_all, outnames, ytrain = utils.read_simdata_ytrain(simdata)


np.savetxt('xdata_all.txt', xdata_all, fmt='%.2f')
np.savetxt('ytrain_all.dat', ytrain)
with open('outnames_all.txt', 'w') as f:
    for item in outnames:
        f.write("%s\n" % item)

