import numpy
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
    spls = numpy.genfromtxt(options.parfile,skip_header=0)
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
p_all   = numpy.zeros([nens,model.nparms],numpy.float)
gpp_all = numpy.zeros([nens,12],numpy.float)  #Monthly GPP (mean seasonal cycle over all years)
nee_all = numpy.zeros([nens,12],numpy.float)  #Monthly NEE (mean seasonal cycle over all years)
lai_all = numpy.zeros([nens,12],numpy.float)  #Monthly LAI (mean seasonal cycle over all years)

# some more checks
if spls.shape[0]>0:
    if npars!=model.nparms:
        print('Not enough parameters in the samples file')
        quit()

names_out = open('pnames.txt','w')
if spls.shape[0]>0:
    for p in splsNames:
        names_out.write(p+'\n')
else:
    for p in model.parms:
        names_out.write(p+'\n')

names_out.close()

for i in range(0,nens):
    print 'Running #'+str(i+1)
    if spls.shape[0]>0:
        # Sample in [-1,1] from file
        for k in range(model.nparms):
            p = splsNames[k]
            model.parms[p]=model.pmin[p]+0.5*(spls[i,k]+1)*(model.pmax[p]-model.pmin[p])
            p_all[i,k] = model.parms[p]
    else:
        pnum=0
        for p in model.parms:
            #Sample uniformly from the parameter space
            model.parms[p] = numpy.random.uniform(low=model.pmin[p],high=model.pmax[p])
            #Save the parameters
            p_all[i,pnum] = model.parms[p]
            pnum=pnum+1
    model.run_selm(spinup_cycles=6,deciduous=True)
    gpp_all[i,:] = daily_to_monthly(model.output['gpp'],allyearsmean=True)
    nee_all[i,:] = daily_to_monthly(model.output['nee'],allyearsmean=True)
    lai_all[i,:] = daily_to_monthly(model.output['lai'],allyearsmean=True)

#Save "UQ-ready" outputs
numpy.savetxt('p_ensemble.txt',p_all)
numpy.savetxt('gpp_ensemble.txt',gpp_all)
numpy.savetxt('nee_ensemble.txt',nee_all)
numpy.savetxt('lai_ensemble.txt',lai_all)
