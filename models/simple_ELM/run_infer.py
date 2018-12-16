#!/usr/bin/env python

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from pylab import *
import os
import numpy

from utils import *
from common import oscm_dir



def infer(xdata_full, ind_calib, ydata, data_err, nmcmc=10000,gamma=0.1,fixindnom=None):
    xdata_all=xdata_full[ind_calib]
    np.savetxt('xdata_all.txt',xdata_all)
    np.savetxt('ydata_calib.txt',ydata)
    np.savetxt('dataerr_calib.dat', data_err)
    np.savetxt('xdata_calib.txt',xdata_all[ind_calib])
    dim = np.loadtxt('mindexp.0.dat').shape[1]

    if fixindnom!=None:
        np.savetxt('fixindnom.dat', fixindnom)
        dim -= fixindnom.shape[0]


    comstr  = 'model_inf -l classical -o 0 -i uniform -t xdata_all.txt -f pcs -c LU -s pci  -u 5 -a -1 -b 1' #uniform_lupci?
    comstr += ' -m '+str(nmcmc)+' -g '+str(gamma)+' -d '+str(dim)
    comstr += ' -z -x xdata_calib.txt -y ydata_calib.txt ' # -e $EE -j chainstart.dat
    if fixindnom!=None:
        comstr += ' -v fixindnom.dat'

    os.command(comstr)

    return

###########################################################################
###########################################################################
###########################################################################

ind_calib = np.loadtxt('ind_calib.dat', dtype=int)

nout_calib = ind_calib.shape[0]
dataerr_calib=np.empty((nout_calib,))
nout_calib = ind_data_all.shape[0]
for iout_calib in range(nout_calib):
    iout=ind_calib[iout_calib]
    np.savetxt('mindexp.'+str(iout)+'.dat', results_all[iout]['mindex'],fmt='%d')
    np.savetxt('pccfp.'+str(iout)+'.dat', results_all[iout]['cfs'])
    dataerr_calib[iout]=results_all[iout]['rmse_val'] #errval[i]*np.linalg.norm(ydata[:,i])/np.sqrt(ndata)


np.savetxt('dataerr_calib.dat',dataerr_calib)


infer(xdata_full, ind_calib, ydata, data_err_calib)
