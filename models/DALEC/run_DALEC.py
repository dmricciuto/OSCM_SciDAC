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
#Run model with default parameters
model.run(spinup_cycles=40, deciduous=True)

#Ouptut model GPP (default parameters)
numpy.savetxt('gpp.txt',model.output['gpp'])

#plot GPP
model.plot_output(startyear=2007,endyear=2008)
