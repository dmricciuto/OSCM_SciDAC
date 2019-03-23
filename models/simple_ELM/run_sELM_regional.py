import numpy
import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *


# create model object
model = models.MyModel()

# Create a 10-member, 3-parameter ensemble
pnames_ens = ['nue_tree', 'gdd_crit'] #'leafcn', 'slatop_decid']
n_ensemble = 50

# Desired outputs
myoutvars = ['gpp', 'lai', 'npp']
print('Generating ensemble')
model.generate_ensemble(n_ensemble, pnames_ens)

print(model.parm_ensemble)

# Daily output
# model.run_selm(spinup_cycles=4,
#                lon_bounds=[-82, -80], lat_bounds=[40, 42],
#                prefix='regional', pft=0, use_MPI=False)

# Monthly output
# model.run_selm(spinup_cycles=4,
#                lon_bounds=[-85,-80], lat_bounds=[40,45],
#                prefix='regional', do_monthly_output=True,
#                pft=0, use_MPI=True)

site_coords = [-72.17, 42.54]  # -72.17 42.54    b'US-Ha1'
dellon = 0.01
dellat = 0.01
model.run_selm(spinup_cycles=4,
               lon_bounds=[site_coords[0] - dellon, site_coords[0] + dellon],
               lat_bounds=[site_coords[1] - dellat, site_coords[1] + dellat],
               prefix='regional', do_monthly_output=True,
               pft=0, use_MPI=False)
