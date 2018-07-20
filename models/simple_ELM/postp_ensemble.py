from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *

ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
qoi='gpp'
lat_id=30
lon_id=40

print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lons=dataset.variables['lon'][:]-360. #.shape
lats=dataset.variables['lat'][:] #.shape

fig=plt.figure(figsize=(12,6))


for ens_id in range(100):
	thisvar=dataset.variables[qoi][ens_id,0,:,lat_id,lon_id]
	plot(1980+dataset.variables['time'][:]/365,thisvar)
#np.savetxt('gpp.txt',thisvar)
#xlim(2000,2010)
#ylim(0,10)
#savefig('regional.eps')
show()
#sys.exit()

