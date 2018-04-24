import numpy
import cPickle as pk

#Transform daily data to monthly (or produce a monthly mean seasonal cycle
def daily_to_monthly(var_in, allyearsmean=False):
    ndays_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    doy_end = numpy.cumsum(ndays_month)  #ending day of year for each month
    ndays = len(var_in)
    if (allyearsmean):
        nmonths=12
    else:
        nmonths = ndays/365*12
    #Day of year for each output day
    doy_in = (numpy.cumsum(numpy.ones([ndays], numpy.float))-1) % 365 + 1
    var_out = numpy.zeros([nmonths], numpy.float)
    for d in range(0,ndays):
        thisyear = d/365
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


def savepk(sobj,nameprefix='savestate'):
    pk.dump(sobj,open(nameprefix+'.pk','wb'),-1)

def loadpk(nameprefix='savestate'):
    return pk.load(open(nameprefix+'.pk','rb'))
