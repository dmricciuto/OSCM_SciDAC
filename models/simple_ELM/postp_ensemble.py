#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

import utils
import model_sELM as selm


ncfile = sys.argv[1] #'regional_output.nc'



print("Reading " + ncfile)
dataset = Dataset(ncfile)
print("Dimensions #######################")
for ikey in dataset.dimensions.keys():  # time(360),lon(5),lat(7)
    print(dataset.dimensions[ikey].name + ", size " +
          str(dataset.dimensions[ikey].size))  # 7

print("Variables #######################")
for ikey in dataset.variables.keys():
    print(dataset.variables[ikey].name + ", size " +
          str(dataset.variables[ikey].shape))

lons = dataset.variables['lon'][:] - 360.  # .shape
lats = dataset.variables['lat'][:]  # .shape


# create model object
model = selm.MyModel()
qois = model.outvars
qois.remove('ctcpools')


years = 1980 + dataset.variables['time'][:] / 365

for qoi in qois:
    qoi_var = utils.get_qoi_regave(dataset, qoi)

    fig = plt.figure(figsize=(12, 6))
    plt.plot(years, qoi_var.T)

    plt.savefig(qoi+'_ens.png')
    plt.clf()

