from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy
import utils


if os.environ['USER']=='ksargsy' or os.environ['USER']=='root':
  print('Hello Khachik')
  oscm_dir=os.environ['HOME']+'/research/OSCM_SciDAC/'
elif os.environ['USER']=='csafta':
  print('Hello Cosmin')
  oscm_dir=os.environ['HOME']+'/Projects/OSCM_SciDAC.dmr/'
else:
  oscm_dir='../../'



ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)

qois = ['gpp']
indlats = [34]
indlons = [23]

nqois = len(qois)
nlats = len(indlats)
nlons = len(indlons)

twelve = 12
nens   = 2000

monthnames=['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lats = dataset.variables['lat'][:] #.shape
lons = dataset.variables['lon'][:]-360. #.shape


outqois = np.empty((nlats,nlons,nqois,2,twelve))
ytrain = np.empty((nens,nlats*nlons*nqois*2*twelve))
outnames=[]
xdata = np.empty((0,5))
jout = 0
for iqoi in range(nqois):
  qoi = qois[iqoi]
  for ilat in range(nlats):
    lat_id = indlats[ilat]
    for ilon in range(nlons):
      lon_id = indlons[ilon]

      fig = plt.figure(figsize=(12,6))
      ytrain_this=np.empty((0,))
      for ens_id in range(nens):
        var_ = dataset.variables[qoi][ens_id,0,:,lat_id,lon_id].reshape(-1,twelve)
        outqois[ilat,ilon,iqoi,0,:] = np.average(var_,axis=0)
        outqois[ilat,ilon,iqoi,1,:] = np.std(var_,axis=0)
        #plot(thisvar_mean,thisvar_std,'o')
        errorbar(range(twelve),outqois[ilat,ilon,iqoi,0],outqois[ilat,ilon,iqoi,1])

        ytrain[ens_id,jout:jout+twelve] = outqois[ilat,ilon,iqoi,0,:]
        ytrain[ens_id,jout+twelve:jout+2*twelve] = outqois[ilat,ilon,iqoi,0,:]

      for i in range(twelve):
        xdata   = np.append(xdata, [[lat_id,lon_id,iqoi,i,0]], axis=0)
        outnames.append(qoi+'_'+str(lat_id)+'_'+str(lon_id)+'_m '+monthnames[i])
      for i in range(twelve):
        xdata   = np.append(xdata, [[lat_id,lon_id,iqoi,i,1]], axis=0)
        outnames.append(qoi+'_'+str(lat_id)+'_'+str(lon_id)+'_s '+monthnames[i])


      jout += 2*twelve

      savefig('ensemble_'+qoi+'_'+str(lat_id)+'_'+str(lon_id)+'.eps')


np.savetxt('xdata.txt', xdata, fmt='%.2f')
np.savetxt('ytrain_all.dat', ytrain, fmt='%.12f')
np.savetxt('outnames.txt',outnames,fmt='%s')


#for ivar in dataset.variables.keys():


sys.exit()


dataset_obs = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
site_id=29

sitenames = dataset_obs.variables['site_name'][:]
nsm_0=sitenames[:,:].shape[0]
nsm_1=sitenames[:,:].shape[1]
siteName_str=numpy.array([[sitenames[i,j].decode("utf-8") for j in range(nsm_1)] for i in range(nsm_0)])
site_name=''.join(siteName_str[site_id])
print(site_name)
lat=dataset_obs.variables['lat'][:][site_id]
lon=dataset_obs.variables['lon'][:][site_id]
print(lat,lon)

var_daily = dataset_obs.variables['GPP'][site_id,:]*24*3600*1000

var_monthly=np.nanmean(utils.daily_to_monthly(var_daily).reshape(24,12),axis=0)
plot(var_monthly, 'ko', markersize=10)
print(var_monthly)

xticks(range(12),monthnames)
ylabel(qoi)
title('Lat = '+ str(lats[lat_id])+', Lon = '+str(lons[lon_id]))
show()




