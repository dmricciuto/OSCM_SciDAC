from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *

dataset = Dataset('regional_output.nc')


print dataset.dimensions.keys() # time(360),lon(5),lat(7)
print dataset.dimensions['lat'] # 7

print dataset.variables.keys()

lons=dataset.variables['lon'][:]-360. #.shape
lats=dataset.variables['lat'][:] #.shape
print lons
print lats
print dataset.variables['time'].units
print dataset.variables['time'][:]