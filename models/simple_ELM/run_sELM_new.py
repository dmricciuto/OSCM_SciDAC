import model_sELM as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import argparse

from utils import *

# default parameters
nens  = 5   # ensemble size
npars = 45 # number of parameters

parser = argparse.ArgumentParser()

# command-line arguments
parser.add_argument("-s","--site",   dest="site",    default="US-UMB", type=str, help="Site ID")
parser.add_argument("-p","--parspl", dest="parfile", default=None,     type=str, help="Parameter Samples")
#parser.add_argument("-m","--parnames", dest="pnamesfile", default=None,     type=str, help="Parameter Names")
options = parser.parse_args()

site = options.site
spls = numpy.array([])
if options.parfile is not None:
    spls = numpy.genfromtxt(options.parfile)
    spls = numpy.atleast_2d(spls)
else:
    spls = numpy.empty((nens,npars))
# if options.pnamesfile is not None:
#     splsNames = numpy.genfromtxt(options.pnamesfile, dtype='str').tolist()


print(spls.shape)
nens  = spls.shape[0] # ensemble size
npars = spls.shape[1] # no. of parameters
# if (len(splsNames)!=npars):
#     print('Parameter values file %s is not consistent with parameter names file %s' \
#           % (options.parfile,options.pnamesfile))

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

splsNames=[]
for v in model.parms:
    splsNames.append(v)

# some more checks
if options.parfile is not None:
    if npars!=model.nparms:
        print('Not enough parameters in the samples file')
        quit()
    else:
       model.generate_ensemble(nens, splsNames, fname=options.parfile, normalized=False)
else:
  #No file provided, generate ensemble using all parameters
  model.generate_ensemble(nens,splsNames)

names_out = open('pnames.txt','w')
for p in splsNames:
  names_out.write(p+'\n')
names_out.close()

model.run_selm(spinup_cycles=6,deciduous=True, do_monthly_output=True)
