import numpy
import itertools
from pylab import *
from netCDF4 import Dataset
from common import oscm_dir
from mpl_toolkits.basemap import Basemap
import matplotlib.tri as tri
from common import pmin, pmax

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

def read_sites(dataset):

  print("Dimensions #######################")
  for ikey in dataset.dimensions.keys():
      print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size))

  print("Variables #######################")
  for ivar in dataset.variables.keys():
    print(ivar+str(dataset.variables[ivar].shape))
    for attr in dataset.variables[ivar].ncattrs():
      print(attr , '=', getattr(dataset.variables[ivar], attr))
  lons=dataset.variables['lon'][:]
  lats=dataset.variables['lat'][:]
  times=dataset.variables['time'][:]
  sitenames=dataset.variables['site_name'][:]

  sites_info=[]
  for j in range(sitenames.shape[0]):
    sname=str(sitenames[j][0]+sitenames[j][1]+sitenames[j][2]+sitenames[j][3]+sitenames[j][4]+sitenames[j][5])
    site_info = [sname, lons[j], lats[j]]
    sites_info.append(site_info)

  return sites_info

def read_simdata_input(dataset):

    print("Dimensions #######################")
    for ikey in dataset.dimensions.keys(): # time(360),lon(5),lat(7)
        print(dataset.dimensions[ikey].name+", size " + str(dataset.dimensions[ikey].size)) # 7

    print("Variables #######################")
    pnames=[]
    nens = dataset.dimensions['ensemble'].size
    ptrain=np.empty((nens,))
    prange_list=[]

    for ivar in dataset.variables.keys():
        print(ivar+str(dataset.variables[ivar].shape))
        for attr in dataset.variables[ivar].ncattrs():
            print(attr , '=', getattr(dataset.variables[ivar], attr))

        # if ivar!='gpp' and ivar!='lai' and \
        #         ivar!='lon' and ivar!='lat' and \
        #         ivar!='time' and ivar!='pft_frac':
        if np.array(dataset.variables[ivar][:]).shape[0]==nens and len(np.array(dataset.variables[ivar][:]).shape)==1:

            pnames.append(ivar)
            print("AAA ", ptrain.shape, np.array(dataset.variables[ivar][:]).shape)
            ptrain=np.vstack((ptrain,np.array(dataset.variables[ivar][:])))
            prange_list.append([pmin[ivar],pmax[ivar]])

    lons=dataset.variables['lon'][:]-360
    lats=dataset.variables['lat'][:]
    times=dataset.variables['time'][:]
    pfts=dataset.variables['pft_frac'][:]

    pnames=np.array(pnames)
    ptrain=ptrain[1:,:].T
    prange=np.array(prange_list)

    qtrain=scaleDomTo01(ptrain, prange)

    return pnames, prange, ptrain, qtrain

def read_simdata_ytrain(dataset):
    qois = ['gpp']
    nqois = len(qois)

    monthnames=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    twelve = len(monthnames)
    print("Dimensions #######################")
    for ikey in dataset.dimensions.keys():
        #ikey==dataset.dimensions[ikey].name
        print(ikey+", size " + str(dataset.dimensions[ikey].size)) # 7

    nens   = dataset.dimensions['ensemble'].size
    npfts  = dataset.dimensions['pft'].size
    nyears = dataset.dimensions['time'].size / twelve

    lons  = dataset.variables['lon'][:]-360. #.shape
    ilons=range(lons.shape[0]) #range(11,14)
    nlons = len(ilons)
    lats  = dataset.variables['lat'][:] #.shape
    ilats=range(lats.shape[0]) #range(36,38)
    nlats = len(ilats)


    bm = Basemap()
    #outqois = np.empty((nsites,nqois,2,twelve))
    ytrain  = np.empty((nens,0))
    outnames=[]
    xdata = np.empty((0,5))
    #ydata = np.empty((0,))
    #jout = 0
    for iqoi in range(nqois):
      qoi = qois[iqoi]
      for ilon in ilons:
        lon = lons[ilon]
        for ilat in ilats:
          lat = lats[ilat]
          if bm.is_land(lon,lat):
            print('lon={:d}, lat={:d} is land'.format(ilon,ilat))
            #print("lon=%d, lat=%d : Land", ilon,ilat)
            aa=dataset.variables[qoi][:,iqoi,:,ilat,ilon].reshape(nens,-1,twelve)
            #print iqoi, ilat, ilon, aa
            ytrain = np.append(ytrain, np.average(aa, axis=1),axis=1)
            ytrain = np.append(ytrain, np.std(aa, axis=1),axis=1)
            for imo in range(1,twelve+1):
              xdata  = np.append(xdata, [[lon,lat,iqoi,imo,0]], axis=0)
              outnames.append(qoi+' '+str(lat)+' '+str(lon)+' mean '+monthnames[imo-1])

            for imo in range(1,twelve+1):
              xdata  = np.append(xdata, [[lon,lat,iqoi,imo,1]], axis=0)
              outnames.append(qoi+' '+str(lat)+' '+str(lon)+' stdev '+monthnames[imo-1])
          else:
            print('lon={:d}, lat={:d} is not land'.format(ilon,ilat))

    return xdata, outnames, ytrain

def scale01ToDom(xx,dom):
    dim=xx.shape[1]
    xxsc=np.zeros((xx.shape[0],xx.shape[1]))
    for i in range(dim):
        xxsc[:,i]=xx[:,i]*(dom[i,1]-dom[i,0])+dom[i,0]

    return xxsc

def scaleDomTo01(xx,dom):
    dim=xx.shape[1]
    xxsc=np.zeros((xx.shape[0],xx.shape[1]))
    for i in range(dim):
        xxsc[:,i]=(xx[:,i]-dom[i,0])/(dom[i,1]-dom[i,0])

    return xxsc

def cartes_list(somelists):

    final_list=[]
    for element in itertools.product(*somelists):
        final_list.append(element)

    return final_list

def pick_ind(data,row):
    ind_data = np.empty((0,),dtype=int)
    assert(data.shape[1]==len(row))
    nr=len(row)
    ndata=data.shape[0]
    inds = []
    for j in range(nr):
        if (row[j]!=None):
            inds.append(np.where(data[:,j]==row[j])[0])
    ind_out = np.arange(ndata)
    for ind in inds:
        ind_out = np.intersect1d(ind_out,ind)
    ind_data = np.append(ind_data,ind_out,axis=0)


    return ind_data

def pick_sites2(lonlats,obs_lonlats):

    nsites=len(obs_lonlats)
    ssind=[]
    for i in range(nsites):
        arr=lonlats-np.array([obs_lonlats[i,0],obs_lonlats[i,1]])
        dists=numpy.linalg.norm(arr,axis=1)
        minind=numpy.argmin(dists)

        if dists[minind]<1.0:
            ssind.append([i, lonlats[minind,0], lonlats[minind,1]])


    return numpy.array(ssind)

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
        print(site_lon_lat)
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
