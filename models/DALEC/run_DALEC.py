import numpy
import model_DALEC as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


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

#run an ensemble
nens = 100     #ensemble size
pall   = numpy.zeros([nens,model.nparms],numpy.float)
gppall = numpy.zeros([nens,365],numpy.float)
for i in range(0,nens):
    print 'Running #'+str(i+1)
    pnum=0
    for p in model.parms:
        #Sample uniformly from the parameter space
        model.parms[p] = numpy.random.uniform(low=model.pmin[p],high=model.pmax[p])
        #Save the parameters
        pall[i,pnum] = model.parms[p]
        pnum=pnum+1
    model.run(spinup_cycles=6,deciduous=True)
    #Save the first year of daily output.  We may want to operate on this (e.g. average to monthly)
    gppall[i,:] = model.output['gpp'][0:365]

#Save "UQ-ready" outputs
numpy.savetxt('ptrain.txt',pall)
numpy.savetxt('ytrain.txt',gppall)
