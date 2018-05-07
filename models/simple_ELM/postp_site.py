from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *

dataset = Dataset('model_output.nc')
qoi='gpp'

print dataset.dimensions.keys() # time(360),lon(5),lat(7)

print dataset.variables['time'].units
print dataset.variables['time'][:]

plot(2000+dataset.variables['time'][:]/365,dataset.variables[qoi][:,0,0])
xlim(2000,2010)
ylim(0,9)
savefig('site.eps')
show()