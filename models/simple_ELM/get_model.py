#!/usr/bin/env python

import numpy
import model_DALEC as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from utils import *


sites=['US-Ha1','US-MMS','US-MOz','US-NR1', 'US-UMB', 'US-WCr']

# model = models.MyModel()
# for site in sites:
# 	print site
# 	model.load_forcings(site=site)
# 	model.run(spinup_cycles=6, deciduous=True)
# 	savepk(model,'model_'+site)

for site in sites:
	print site
	model=loadpk('model_'+site)
	model.plot_output(startyear=model.start_year,endyear=model.end_year,obs=False,figname_postfix='_'+site)
	print model.start_year, model.end_year
	print model.output['gpp'].shape