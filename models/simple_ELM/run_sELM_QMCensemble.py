import numpy
import model_sELM as models
import os, math, random
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt

from utils import *


#create model object
model = models.MyModel()

#site level data
#site = 'US-UMB'
#load site forcings
#model.load_forcings(site=site)
#Load site observations
#model.load_obs(site)

#---------------------- run an ensemble ----------------------------

#pnames_ens = ['leafcn','slatop_decid','q10_mr','frootcn','froot_leaf','br_mr','crit_dayl','gdd_crit']

#generate an ensemble with all parameters
pnames_ens=[]
for v in model.parms:
    pnames_ens.append(v)

n_ensemble = 2000
myoutvars = ['gpp','lai']
print 'Generating ensemble'

#model.generate_ensemble(n_ensemble, pnames_ens, fname='QMCsample2.dat')
model.generate_ensemble(n_ensemble, pnames_ens)

#model.run_selm(spinup_cycles=6, pft=0, do_monthly_output=True, use_MPI=True)
model.run_selm(spinup_cycles=4, lon_bounds=[-96,-66], lat_bounds=[28,48], pft=0, do_monthly_output=True,  myoutvars=myoutvars, use_MPI=True)

