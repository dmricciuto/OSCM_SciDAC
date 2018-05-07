from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *

dataset = Dataset('regional_output.nc')
qoi='gpp'

print dataset.dimensions.keys() # time(360),lon(5),lat(7)
print dataset.dimensions['lat'] # 7

print dataset.variables.keys()

lons=dataset.variables['lon'][:]-360. #.shape
lats=dataset.variables['lat'][:] #.shape
print lons,lats
#print dataset.variables['time'].units
# for attr in dataset.variables['time'].ncattrs():
#     print attr #, '=', getattr(windspeed, attr)
#
for i in range(len(lats)):
    for j in range(len(lons)):
        plot(1980+dataset.variables['time'][:]/365,dataset.variables[qoi][:,i,j])
        xlim(2000,2010)
        ylim(0,9)
        savefig('regional.eps')
        show()
        sys.exit()

fig=plt.figure(figsize=(12,6))
#ax=fig.add_axes([0.1,-0.1,0.8,1.3])




m = Basemap(llcrnrlon=-120.,llcrnrlat=30.,urcrnrlon=-60.,urcrnrlat=50.,\
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
#plt.plot(x,y,'o',markersize=27,label='PFT')



imshow(dataset.variables[qoi][-1,:,:], \
        extent=(x.min(),x.max(), y.min(),y.max()), \
        interpolation='bilinear', cmap=cm.RdYlGn)
colorbar(fraction=.02,pad=0.1) #location='bottom', pad="10%")
#
show()