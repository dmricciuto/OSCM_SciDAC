import numpy
import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *


#create model object
model = models.MyModel()

#model.run_selm(spinup_cycles=4, lon_bounds=[-90,-88], lat_bounds=[40,43], prefix='regional', do_monthly_output=True)

#model.run_selm(spinup_cycles=0, lon_bounds=[-84.5,-84.2], lat_bounds=[45,45.5], prefix='regional', do_monthly_output=True, deciduous=True)
#model.run_selm(spinup_cycles=4, lon_bounds=[-72.1,-71.9], lat_bounds=[42.3,42.6], prefix='regional', do_monthly_output=True, deciduous=True)
model.run_selm(spinup_cycles=4, lon_bounds=[-86.6,-86.1], lat_bounds=[39.1,39.6], prefix='regional', do_monthly_output=True, deciduous=True)
