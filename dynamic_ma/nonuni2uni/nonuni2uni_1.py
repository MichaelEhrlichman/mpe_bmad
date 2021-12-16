#!/usr/bin/env python3

import numpy as np
from sys import argv
from scipy import interpolate as interp

CASdata = np.loadtxt(argv[1])

ndata = CASdata.shape[0]

#Assume negative numbers are artifacts and remove them
CASdata[:,1] = np.maximum(CASdata[:,1],[0]*ndata)

nuni = 1000
fdelta = 0.1
CASuni = np.zeros((nuni,2))
funi = np.zeros(nuni)

for ix in range(nuni):
	funi[ix] = round( (ix+1)*fdelta, 3)

CAS = interp.interp1d(CASdata[:,0],CASdata[:,1])

CASuni[:,0] = funi
CASuni[:,1] = np.maximum(CAS(funi),[0]*nuni)

with open(argv[1]+'.uniform','w') as f:
		for freq,dat in CASuni:
				f.write('{:.4f}   {:13.5e}\n'.format(freq,dat))

