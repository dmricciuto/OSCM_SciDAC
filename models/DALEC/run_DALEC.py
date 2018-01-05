import numpy
import model_DALEC as models
import os, math, random
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#create model object
model = models.MyModel()
#Load model forcings for Harvard Forest
model.load_forcings('US-Ha1')
#Run model with default parameters
model.run(model.pdef)

#Ouptut model GPP (default parameters)
numpy.savetxt('gpp.txt',model.gpp)
