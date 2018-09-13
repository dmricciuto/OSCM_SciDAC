from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import itertools

from utils import *
from common import oscm_dir

ncfile=sys.argv[1]
#ncfile='regional_output.nc'
print("Reading "+ncfile)

dataset = Dataset(ncfile)
qoi='gpp'
ens_id=111
time_id=30

print("Dimensions #######################")
for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

lons=dataset.variables['lon'][:]-360. #.shape
lats=dataset.variables['lat'][:] #.shape





fig=plt.figure(figsize=(12,6))
#ax=fig.add_axes([0.1,-0.1,0.8,1.3])




m = Basemap(llcrnrlon=-120.,llcrnrlat=20.,urcrnrlon=-60.,urcrnrlat=50.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc',\
            lat_0=40.,lon_0=-20.,lat_ts=0.)

m.drawcoastlines()
#m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(50,30,-10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-120,-60,30),labels=[1,1,0,1])
m.drawlsmask(land_color='Linen', ocean_color='#CCFFFF')
#m.drawstates()
#m.drawcountries()

xx,yy=np.meshgrid(lons,lats)
x,y = m(xx, yy)
#plt.plot(x,y,'o',markersize=2,label='PFT')


obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
site_names, site_lons, site_lats = read_obsdata(obs_dataset)

site_lon_lat =  pick_sites(site_lons,site_lats,lons,lats)




imshow(dataset.variables[qoi][ens_id,0,time_id,::-1,:], \
        extent=(x.min(),x.max(), y.min(),y.max()), \
        interpolation='bilinear', cmap=cm.RdYlGn)
colorbar(fraction=.02,pad=0.1) #location='bottom', pad="10%")

print site_lon_lat

site_lons_mapped,site_lats_mapped= m(site_lons,site_lats)
plt.plot(site_lons_mapped,site_lats_mapped,'ko',markersize=4,zorder=1000)
plt.plot(site_lons_mapped[site_lon_lat[:,0]],site_lats_mapped[site_lon_lat[:,0]],'ro',markersize=4,zorder=2000, label='selected')


savefig('map.eps')
show()