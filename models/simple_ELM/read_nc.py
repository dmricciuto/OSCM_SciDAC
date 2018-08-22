from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *


ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
#qoi='gpp'



print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

print("Variables #######################")
for ivar in dataset.variables.keys():
    print(ivar+str(dataset.variables[ivar].shape))
    for attr in dataset.variables[ivar].ncattrs():
    	print(attr , '=', getattr(dataset.variables[ivar], attr))
    if ivar=='lon':
    	lons=dataset.variables[ivar][:]-360
    elif ivar=='lat':
    	lats=dataset.variables[ivar][:]
    elif ivar=='time':
    	times=dataset.variables[ivar][:]
    elif ivar=='pft_frac':
    	pfts=dataset.variables[ivar][:]
    elif ivar=='gpp':
    	pass
    	#gpps=dataset.variables[ivar][:]
    elif ivar=='lai':
    	pass
    	#lais=dataset.variables[ivar][:]
    else:
    	pass


print(lons)
print(lats)
# print times
# print pfts[:,:,0]