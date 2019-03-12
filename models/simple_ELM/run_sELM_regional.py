import numpy
import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *


#create model object
model = models.MyModel()

#Create a 10-member, 3-parameter ensemble
pnames_ens = ['nue_grass','leafcn','slatop_decid']
n_ensemble = 10
#Desired outputs
myoutvars = ['gpp','lai','npp']
print('Generating ensemble')
model.generate_ensemble(n_ensemble, pnames_ens)

#Daily output
model.run_selm(spinup_cycles=4, lon_bounds=[-82,-80], lat_bounds=[40,42], prefix='regional', pft=0, use_MPI=True)

#Monthly output
#model.run_selm(spinup_cycles=4, lon_bounds=[-85,-80], lat_bounds=[40,45], prefix='regional', do_monthly_output=True, pft=0, use_MPI=True)
