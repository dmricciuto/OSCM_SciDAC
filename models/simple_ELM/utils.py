import numpy
import itertools

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


def pick_sites(obs_lons,obs_lats,lons,lats):
    nsites=len(obs_lons)
    assert(nsites==len(obs_lats))
    ssind=[]
    for i in range(nsites):
        dists=numpy.linalg.norm(numpy.array(cartes_list([lons-obs_lons[i],lats-obs_lats[i]])),axis=1)
        minind=numpy.argmin(dists)
        indlat=minind%lats.shape[0]
        indlon=int(minind/lats.shape[0])
        if dists[minind]<1.0:
            ssind.append([i, indlon, indlat])
    return ssind


def savepk(sobj,nameprefix='savestate'):
    pk.dump(sobj,open(nameprefix+'.pk','wb'),-1)

def loadpk(nameprefix='savestate'):
    return pk.load(open(nameprefix+'.pk','rb'))
