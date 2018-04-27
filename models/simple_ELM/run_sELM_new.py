import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse

from utils import *

# default parameters
nens = 5 # ensemble size

parser = argparse.ArgumentParser()

# command-line arguments
parser.add_argument("-s","--site",   dest="site",    default="US-UMB", type=str, help="Site ID")
parser.add_argument("-p","--parspl", dest="parfile", default=None,     type=str, help="Parameter Samples")
options = parser.parse_args()

site = options.site
spls = numpy.array([])
if options.parfile is not None:
    infile = open(options.parfile, 'r')
    firstLine = infile.readline()
    infile.close()
    splsNames=firstLine.split()
    spls = numpy.genfromtxt(options.parfile,skip_header=1)
    spls = numpy.atleast_2d(spls) 
    print(spls.shape)
    nens  = spls.shape[0] # ensemble size
    npars = spls.shape[1] # no. of parameters
    if (len(splsNames)!=npars):
        print('Parameter file %s is not consistent'%(options.parfile))

print('Processing site: %s'%(site))

#create model object
model = models.MyModel()

#Load model forcings 
model.load_forcings(site=site)

#Load site observations
model.load_obs(site)

# ------------------- Run and analyze the default model ---------
#Run model with default parameters
model.run_selm(spinup_cycles=6, deciduous=True)
#Ouptut model GPP
numpy.savetxt('gpp_model.txt',model.output['gpp'])
#Output observed GPP
numpy.savetxt('gpp_obs.txt',model.obs['gpp'])
#plot all variables for the default run
#model.plot_output(startyear=2000,endyear=2005)

#---------------------- run an ensemble ----------------------------

# some more checks
if spls.shape[0]>0:
    if npars!=model.nparms:
        print('Not enough parameters in the samples file')
        quit()
    else:
       model.generate_ensemble(nens, splsNames, file=options.parfile, normalized=True)
else:
  #No file provided, generate ensemble using all parameters
  splsNames=[]
  for v in model.parms:
    splsNames.append(v)
  model.generate_ensemble(nens,splsNames)

names_out = open('pnames.txt','w')
for p in splsNames:
  names_out.write(p+'\n')
names_out.close()

model.run_selm(spinup_cycles=6,deciduous=True, do_monthly_output=True)
