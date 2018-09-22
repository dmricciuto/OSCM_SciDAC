#!/usr/bin/env python


import sys, os
import numpy as np

#sys.path.append(os.environ['UQTK_SRC'])
#import PyUQTk.plotting.surrogate as ss


sys.path.append(os.environ['WW']+'/run/xuq')
import workflows as wf
import myutils as mu

from pylab import *
import cPickle as pk

import warnings
#from sklearn.gaussian_process.kernels import Matern,RBF,WhiteKernel, ExpSineSquared

from utils import *


rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)
#rc('font', family='serif')
rc('figure',max_open_warning=100)

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
print ind_out
meanMode, klModes, eigValues, xiSam, relDiag=wf.kle(xdata[ind_out], ysam[:,ind_out].T, neig = 20)


lats=np.loadtxt('lats.txt')
lons=np.loadtxt('lons.txt')
nsites=ind_out.shape[0]
data2d=np.empty((nsites,3))
for j in range(nsites):
    data2d[j,0]=lons[int(xdata[ind_out[j],1])]
    data2d[j,1]=lats[int(xdata[ind_out[j],0])]
data2d[:,2]=meanMode
print data2d
plotMap2(data2d)


