import numpy
import itertools
from pylab import *
from netCDF4 import Dataset
from common import oscm_dir
from mpl_toolkits.basemap import Basemap
import matplotlib.tri as tri

try:
    import cPickle as pk
except:
    import pickle as pk

#Transform daily data to monthly (or produce a monthly mean seasonal cycle
def daily_to_monthly(var_in, allyearsmean=False):
    ndays_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    doy_end = numpy.cumsum(ndays_month)  #ending day of year for each month
    ndays = len(var_in)
    if (allyearsmean):
        nmonths=12
    else:
        nmonths = int(ndays/365)*12
    #Day of year for each output day
    doy_in = (numpy.cumsum(numpy.ones([ndays], numpy.float))-1) % 365 + 1
    var_out = numpy.zeros([nmonths], numpy.float)
    for d in range(0,ndays):
        thisyear = int(d/365)
        for m in range(0,12):
            if (m == 0 and doy_in[d] <= doy_end[m]):
                thismonth = m
            elif (m > 0 and doy_in[d] <= doy_end[m] and doy_in[d] > doy_end[m-1]):
                thismonth=m
        thismonth_out = thismonth
        if (not allyearsmean):
            thismonth_out = thismonth+thisyear*12
            var_out[thismonth_out] = var_out[thismonth_out] + var_in[d]/(ndays_month[thismonth])
        else:
            var_out[thismonth] = var_out[thismonth] + var_in[d]/(ndays_month[thismonth]*ndays/365)
    return var_out

def read_obsdata(dataset):

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
      snames=[]
      for j in range(sitenames.shape[0]):
        sname=str(sitenames[j][0]+sitenames[j][1]+sitenames[j][2]+sitenames[j][3]+sitenames[j][4]+sitenames[j][5])
        snames.append(sname)

  return snames, lons, lats

def cartes_list(somelists):

    final_list=[]
    for element in itertools.product(*somelists):
        final_list.append(element)

    return final_list


def pick_sites(lons,lats,obs_lons=[],obs_lats=[]):

    if obs_lons!=[]:
        assert(obs_lats!=[])
        nsites=len(obs_lons)
        assert(nsites==len(obs_lats))
        ssind=[]
        for i in range(nsites):
            dists=numpy.linalg.norm(numpy.array(cartes_list([lons-obs_lons[i],lats-obs_lats[i]])),axis=1)
            minind=numpy.argmin(dists)
            indlat=minind%lats.shape[0]
            indlon=int(minind/lats.shape[0])
            if indlon==23 and indlat==35:
                indlat-=1
            if dists[minind]<1.0:
                ssind.append([i, indlon, indlat])
    else:
        nlons = len(lons)
        nlats = len(lats)
        ssind = []
        i = 0
        clist = cartes_list([range(nlons),range(nlats)])
        bm = Basemap()   # default: projection='cyl'
        for indlon, indlat in clist:
            if bm.is_land(lons[indlon],lats[indlat]):
                ssind.append([i, indlon, indlat])
                i = i + 1

    return numpy.array(ssind)


# def create_sites(lons, lats):
#   nlons = len(lons)
#   nlats = len(lats)
#   ssind = []
#   i = 0
#   clist = cartes_list([range(nlons),range(nlats)])
#   bm = Basemap()   # default: projection='cyl'
#   for indlon, indlat in clist:
#     if bm.is_land(lons[indlon],lats[indlat]):
#         ssind.append([i, indlon, indlat])
#         i = i + 1

#   return numpy.array(ssind)





def plotMap(data2d, lats, lons, show_map = False, show_dataloc = False):

    fig=figure(figsize=(12,6))
    #ax=fig.add_axes([0.1,-0.1,0.8,1.3])


    if show_map:
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

    else:
        def m(a,b):
            return a,b




    xx,yy=np.meshgrid(lons,lats)
    x,y = m(xx, yy)
    #plt.plot(x,y,'o',markersize=2,label='PFT')


    #tricontourf(x.flatten(),y.flatten(),data2d.flatten(), 20)
    imshow(data2d, \
            extent=(x.min(),x.max(), y.min(),y.max()), \
            interpolation='bilinear', cmap=cm.RdYlGn, origin='lower')


    colorbar(fraction=.02,pad=0.1) #location='bottom', pad="10%")

    if show_dataloc:
        obs_dataset = Dataset(oscm_dir+"/models/site_observations/fluxnet_daily_obs.nc4",'r',format='NETCDF4')
        site_names, site_lons, site_lats = read_obsdata(obs_dataset)

        site_lon_lat =  pick_sites(lons,lats,site_lons,site_lats)
        print site_lon_lat
        site_lons_mapped,site_lats_mapped= m(site_lons,site_lats)
        plot(site_lons_mapped,site_lats_mapped,'ko',markersize=4,zorder=1000)
        plot(site_lons_mapped[site_lon_lat[:,0]],site_lats_mapped[site_lon_lat[:,0]],'ro',markersize=4,zorder=2000, label='selected')


    savefig('map.eps')
    show()


def plotMap2(data2d):

    fig=figure(figsize=(12,6))
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


    x,y=m(data2d[:,0],data2d[:,1])
    #plot(x,y,'o')
    triang = tri.Triangulation(x, y)
    xmid = x[triang.triangles].mean(axis=1)
    ymid = y[triang.triangles].mean(axis=1)
    aa=[m.is_land(xm,ym) for xm,ym in zip(xmid,ymid)]
    mask = np.where(aa, 0, 1)
    triang.set_mask(mask)

    #tripcolor(x,y,z)

    tricontourf(triang,data2d[:,2], 20)
    colorbar(fraction=.02,pad=0.1) #location='bottom', pad="10%")



    savefig('map2.eps')
    #show()
    return


def savepk(sobj,nameprefix='savestate'):
    pk.dump(sobj,open(nameprefix+'.pk','wb'),-1)

def loadpk(nameprefix='savestate'):
    return pk.load(open(nameprefix+'.pk','rb'))
