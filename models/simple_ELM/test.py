#!/usr/bin/env python

import time
import numpy as np
try:
    import cPickle as pk
except:
    import pickle as pk

modeldata={}
print 1
modeldata['ytrain_full']=np.loadtxt('ytrain_full.dat')
print 2
# pk.dump(modeldata,open('modeldata.pk','wb'),-1)
np.save('ytrain_full.npy', modeldata['ytrain_full'])

# Write random floating point numbers as string on a local CSV file
# with open('ytrain_full.dat', 'w') as fdata:
#     for _ in range(n_samples):
#         fdata.write(str(10*np.random.random())+',')
# Read the CSV in a list, convert to ndarray (reshape just for fun) and time it
# t1=time.time()
# with open('ytrain_full.dat','r') as fdata:
#     datastr=fdata.read()
# lst = datastr.split(' ')
# array_lst=np.array(lst,dtype=float)
# t2=time.time()
# print(t2-t1)