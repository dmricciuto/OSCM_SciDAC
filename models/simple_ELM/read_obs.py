from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os


if os.environ['USER']=='ksargsy' or os.environ['USER']=='root':
  print('Hello Khachik')
  oscm_dir=os.environ['HOME']+'/research/OSCM_SciDAC/'
elif os.environ['USER']=='csafta':
  print('Hello Cosmin')
  oscm_dir=os.environ['HOME']+'/Projects/OSCM_SciDAC.dmr/'
else:
  oscm_dir='../../'




dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')

print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): 
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size))

print("Variables #######################")
for ivar in dataset.variables.keys():
  print(ivar+str(dataset.variables[ivar].shape))
  for attr in dataset.variables[ivar].ncattrs():
    print(attr , '=', getattr(dataset.variables[ivar], attr))
  if ivar=='lon':
  	lons=dataset.variables[ivar][:]
  elif ivar=='lat':
  	lats=dataset.variables[ivar][:]
  elif ivar=='time':
  	times=dataset.variables[ivar][:]
  elif ivar=='site_name':
  	sitenames=dataset.variables[ivar][:]

print(lons)
print(lats)
print(times)
print(sitenames)

sys.exit()

