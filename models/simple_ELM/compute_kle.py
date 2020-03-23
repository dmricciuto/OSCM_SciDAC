#!/usr/bin/env python


import sys, os
import numpy as np

#sys.path.append(os.environ['UQTK_SRC'])
#import PyUQTk.plotting.surrogate as ss


sys.path.append(os.environ['WW']+'/run/xuq')
import workflows as wf


import warnings
#from sklearn.gaussian_process.kernels import Matern,RBF,WhiteKernel, ExpSineSquared

from utils import *


from fitpy.utils.common import myrc


rcParams = myrc()

##########################################################
##########################################################
##########################################################
##########################################################

xdata_file    = 'xdata_full.txt'
ytrain_file   = 'ytrain_full.dat'


xdata    = np.atleast_2d(np.loadtxt(xdata_file))
ysam     = np.atleast_2d(np.loadtxt(ytrain_file))

ind_out = np.arange(xdata.shape[0])
inds = []
#inds.append(np.where(xdata[:,0]==40.0)[0])
#inds.append(np.where(xdata[:,1]==28.0)[0])
inds.append(np.where(xdata[:,3]==6.0)[0])
inds.append(np.where(xdata[:,4]==0.0)[0])

for ind in inds:
    ind_out = np.intersect1d(ind_out,ind)
print(ind_out)
meanMode, klModes, eigValues, xiSam, relDiag=wf.kle(xdata[ind_out], ysam[:,ind_out].T, neig = 8)


lats=np.loadtxt('lats.txt')
lons=np.loadtxt('lons.txt')
nsites=ind_out.shape[0]
data2d=np.empty((nsites,3))
for j in range(nsites):
    data2d[j,0]=lons[int(xdata[ind_out[j],1])]
    data2d[j,1]=lats[int(xdata[ind_out[j],0])]
data2d[:,2]=meanMode
print(data2d)
plotMap2(data2d)


