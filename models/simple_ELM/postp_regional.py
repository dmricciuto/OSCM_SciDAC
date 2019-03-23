from netCDF4 import Dataset
from pylab import *
import itertools

from utils import *
from common import oscm_dir

ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
qoi='gpp'
ens_id=1
time_id=30
pft_id=0

print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lats=dataset.variables['lat'][:] #.shape
lons=dataset.variables['lon'][:]-360. #.shape


np.savetxt('lons.txt',lons)
np.savetxt('lats.txt',lats)


print(lats)
print(lons)

data2d = dataset.variables[qoi][ens_id,pft_id,time_id,:,:]
print(data2d)

# only works if there is at least 2 cells in each direction
plotMap(data2d, lats, lons, show_map = True, show_dataloc = True)

#xdata_full = np.loadtxt('xdata_all.txt')
#plotMap(data2d, xdata_full[:,1], xdata_full[:,0], show_map = True, show_dataloc = True)


