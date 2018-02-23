import numpy
import model_DALEC as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#Transform daily data to monthly (or produce a monthly mean seasonal cycle
def daily_to_monthly(var_in, allyearsmean=False):
    ndays_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    doy_end = numpy.cumsum(ndays_month)  #ending day of year for each month
    ndays = len(var_in[1:])
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

site = 'US-UMB'
#create model object
model = models.MyModel()

#Load model forcings 
model.load_forcings(site=site)

#Load site observations
model.load_obs(site)

# ------------------- Run and analyze the default model ---------
#Run model with default parameters
model.run(spinup_cycles=6, deciduous=True)
#Ouptut model GPP
numpy.savetxt('gpp_model.txt',model.output['gpp'])
#Output observed GPP
numpy.savetxt('gpp_obs.txt',model.obs['gpp'])
#plot all variables for the default run
model.plot_output(startyear=2000,endyear=2005)

#---------------------- run an ensemble ----------------------------

nens = 50     #ensemble size
p_all = numpy.zeros([nens,model.nparms],numpy.float)
gpp_all = numpy.zeros([nens,12],numpy.float)  #Monthly GPP (mean seasonal cycle over all years)
nee_all = numpy.zeros([nens,12],numpy.float)  #Monthly NEE (mean seasonal cycle over all years)
lai_all = numpy.zeros([nens,12],numpy.float)  #Monthly LAI (mean seasonal cycle over all years)

names_out = open('pnames.txt','w')
for i in range(0,nens):
    print 'Running #'+str(i+1)
    pnum=0
    for p in model.parms:
        #Sample uniformly from the parameter space
        model.parms[p] = numpy.random.uniform(low=model.pmin[p],high=model.pmax[p])
        #Save the parameters
        p_all[i,pnum] = model.parms[p]
        if (i == 0):
            names_out.write(p+'\n')
        pnum=pnum+1
    model.run(spinup_cycles=6,deciduous=True)
    gpp_all[i,:] = daily_to_monthly(model.output['gpp'],allyearsmean=True)
    nee_all[i,:] = daily_to_monthly(model.output['nee'],allyearsmean=True)
    lai_all[i,:] = daily_to_monthly(model.output['lai'],allyearsmean=True)
names_out.close()

#Save "UQ-ready" outputs
numpy.savetxt('p_ensemble.txt',p_all)
numpy.savetxt('gpp_ensemble.txt',gpp_all)
numpy.savetxt('nee_ensemble.txt',gpp_all)
numpy.savetxt('lai_ensemble.txt',gpp_all)
