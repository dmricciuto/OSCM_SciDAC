from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

from utils import *
from common import oscm_dir



ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)
dataset = Dataset(ncfile)

qois = ['gpp']

get_obsdata = False


nqois = len(qois)
twelve = 12
nens   = 2000
monthnames=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']


print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lats = dataset.variables['lat'][:] #.shape
lons = dataset.variables['lon'][:]-360. #.shape

if get_obsdata:
  obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
  site_names, site_lons, site_lats = read_obsdata(obs_dataset)
  site_lon_lat =  pick_sites(lons,lats,site_lons,site_lats)
else:
  site_lon_lat = pick_sites(lons, lats)#[::500]

nsites=len(site_lon_lat)
# print lons, lats
# data2d=np.empty((nsites,3))
# j=0
# for ilon,ilat in site_lon_lat[:,1:]:
#   print ilon, ilat, lons[ilon], lats[ilat]
#   data2d[j,0]=lons[ilon]
#   data2d[j,1]=lats[ilat]
#   data2d[j,2]=lons[ilon]+lats[ilat]
#   j += 1
# plotMap2(data2d)

# sys.exit()


outqois = np.empty((nsites,nqois,2,twelve))
ytrain = np.empty((nens,nsites*nqois*2*twelve))
outnames=[]
xdata = np.empty((0,5))
ydata = np.empty((0,))
jout = 0
for iqoi in range(nqois):
  qoi = qois[iqoi]
  for site_id,lon_id,lat_id in site_lon_lat:
    print lon_id, lat_id
    aa=dataset.variables[qoi][:,0,:,lat_id,lon_id].reshape(nens,-1,twelve)
    ytrain[:,jout:jout+twelve]           = np.average(aa, axis=1)
    ytrain[:,jout+twelve:jout+2*twelve]  = np.std(aa,axis=1)

    if get_obsdata:
      var_daily = obs_dataset.variables[qoi.upper()][site_id,:]*24*3600*1000
      var_monthly=np.nanmean(daily_to_monthly(var_daily).reshape(24,twelve),axis=0)
      ydata   = np.append(ydata, var_monthly)
    for i in range(twelve):
      xdata   = np.append(xdata, [[lat_id,lon_id,iqoi,i,0]], axis=0)
      outnames.append(qoi+'_'+str(lat_id)+'_'+str(lon_id)+'_m '+monthnames[i])
    for i in range(twelve):
      xdata   = np.append(xdata, [[lat_id,lon_id,iqoi,i,1]], axis=0)
      #ydata   = np.append(ydata, -999.)
      outnames.append(qoi+'_'+str(lat_id)+'_'+str(lon_id)+'_s '+monthnames[i])


    jout += 2*twelve


np.savetxt('xdata_full.txt',  xdata,  fmt='%.2f')
np.savetxt('ytrain_full.dat', ytrain, fmt='%.12f')
np.savetxt('outnames_full.txt',outnames,fmt='%s')

if get_obsdata:
  np.savetxt('ydata_all.txt', ydata, fmt='%.12f')





# site_id=29
# sitenames = dataset_obs.variables['site_name'][:]
# nsm_0=sitenames[:,:].shape[0]
# nsm_1=sitenames[:,:].shape[1]
# siteName_str=numpy.array([[sitenames[i,j].decode("utf-8") for j in range(nsm_1)] for i in range(nsm_0)])
# site_name=''.join(siteName_str[site_id])
# print(site_name)
# lat=dataset_obs.variables['lat'][:][site_id]
# lon=dataset_obs.variables['lon'][:][site_id]
# print(lat,lon)


#plot(var_monthly, 'ko', markersize=10)
#print(var_monthly)


# xticks(range(12),monthnames)
# ylabel(qoi)
# title('Lat = '+ str(lats[lat_id])+', Lon = '+str(lons[lon_id]))
# #show()



# [[ 5 28 40]
#  [22 48 29]
#  [23 13 36]
#  [24 25 27]
#  [25 12 35]
#  [27 14 36]
#  [29 23 35]
#  [31 12 35]]
