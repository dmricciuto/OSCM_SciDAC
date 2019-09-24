import sys
import numpy as np
from common import pmin, pmax

nens = int(sys.argv[1])
pnames_file = sys.argv[2]
pinput_file = sys.argv[3]

pnames = np.loadtxt(pnames_file, dtype=str)
ndim = pnames.shape[0]

ptrain = np.empty((nens, ndim))
j = 0
for p in pnames:
    ptrain[:, j] = np.random.rand(nens,) * (pmax[p] - pmin[p]) + pmin[p]
    j += 1

np.savetxt(pinput_file, ptrain)
