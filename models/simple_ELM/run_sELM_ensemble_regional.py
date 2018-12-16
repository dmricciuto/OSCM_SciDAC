import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse

from utils import *

nens = 5
model = models.MyModel()

splsNames=[]
for v in model.parms:
    splsNames.append(v)

model.generate_ensemble(nens, splsNames)

# names_out = open('pnames.txt','w')
# for p in splsNames:
#   names_out.write(p+'\n')
# names_out.close()

model.run_selm(spinup_cycles=4, lon_bounds=[-90,-86], lat_bounds=[40,44], prefix='regional', use_MPI=True, do_monthly_output=True)
