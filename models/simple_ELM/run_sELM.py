import numpy
import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *

site = 'US-UMB'
#create model object
model = models.MyModel()

#Load site observations
model.load_obs(site)
model.run_selm(spinup_cycles=4, lon_bounds=[-110,-50], lat_bounds=[22,62], do_monthly_output=True)
model.load_forcings(site=site)

# ------------------- Run and analyze the default model ---------

#Run model with default parameters
model.load_forcings(site=site)
#model.run_selm(spinup_cycles=6, deciduous=True, prefix=site)
#model.run_selm(spinup_cycles=4, lon_bounds=[-90,-89], lat_bounds=[42,43], do_monthly_output=True)

#Ouptut model GPP
#numpy.savetxt('gpp_model.txt',model.output['gpp'])
#Output observed GPP
#numpy.savetxt('gpp_obs.txt',model.obs['gpp'])
##plot all variables for the default run
#model.plot_output(startyear=2000,endyear=2005)

#---------------------- run an ensemble ----------------------------

pnames_ens = ['slatop', 'leafcn']
n_ensemble = 100

model.generate_ensemble(pnames_ens, n_ensemble)
model.run_selm(spinup_cycles=6, deciduous=True, prefix=site)


#Save "UQ-ready" outputs
numpy.savetxt('p_ensemble.txt',p_all)
numpy.savetxt('gpp_ensemble.txt',gpp_all)
numpy.savetxt('nee_ensemble.txt',nee_all)
numpy.savetxt('lai_ensemble.txt',lai_all)
