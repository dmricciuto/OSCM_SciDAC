import numpy
import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *

#create model object
model = models.MyModel()

#Load site observations
model.run_selm(spinup_cycles=4, lon_bounds=[-90,-88], lat_bounds=[40,42], prefix='regional', do_monthly_output=True)

