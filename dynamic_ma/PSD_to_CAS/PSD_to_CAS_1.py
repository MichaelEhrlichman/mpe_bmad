#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from scipy import interpolate as interp

PSDdata = np.loadtxt(argv[1])
ndata = PSDdata.shape[0]

CASdata = np.zeros_like(PSDdata)
CASdata[:,0] = PSDdata[:,0]
for ix in range(ndata):
	CASdata[ix,1] = np.sqrt(np.trapz(PSDdata[ix:,1],PSDdata[ix:,0]))

#plt.semilogy(PSDdata[:,0],PSDdata[:,3], 'o')
#plt.semilogy(PSDuni[:,0],PSDuni[:,3], 'o')
#plt.show()

with open(argv[1]+'.CASdata','w') as f:
		for freq,x in CASdata:
				f.write("{} {}\n".format(freq,x,))
