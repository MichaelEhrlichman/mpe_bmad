#!/usr/bin/env python

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from scipy import interpolate as interp

PSDdata = np.loadtxt(argv[1])
ndata = PSDdata.shape[0]

nuni = 1999
fdelta = 0.05
PSDuni = np.zeros((nuni,4))
funi = np.zeros(nuni)

for ix in range(nuni):
	funi[ix] = round( (ix+1)*fdelta, 5)

for ff in PSDdata[:,0]:
	print(ff)

PSDfx = interp.interp1d(PSDdata[:,0],PSDdata[:,1])
PSDfy = interp.interp1d(PSDdata[:,0],PSDdata[:,2])
PSDfz = interp.interp1d(PSDdata[:,0],PSDdata[:,3])

PSDuni[:,0] = funi
PSDuni[:,1] = PSDfx(funi)
PSDuni[:,2] = PSDfy(funi)
PSDuni[:,3] = PSDfz(funi)

CASuni = np.zeros_like(PSDuni)
CASuni[:,0] = PSDuni[:,0]
for ix in range(nuni):
	CASuni[ix,1] = np.trapz(PSDuni[ix:,1],dx=fdelta)
	CASuni[ix,2] = np.trapz(PSDuni[ix:,2],dx=fdelta)
	CASuni[ix,3] = np.trapz(PSDuni[ix:,3],dx=fdelta)
	
CASdata = np.zeros_like(PSDdata)
CASdata[:,0] = PSDdata[:,0]
for ix in range(ndata):
	CASdata[ix,1] = np.trapz(PSDdata[ix:,1],PSDdata[ix:,0])
	CASdata[ix,2] = np.trapz(PSDdata[ix:,2],PSDdata[ix:,0])
	CASdata[ix,3] = np.trapz(PSDdata[ix:,3],PSDdata[ix:,0])

#plt.semilogy(PSDdata[:,0],PSDdata[:,3], 'o')
#plt.semilogy(PSDuni[:,0],PSDuni[:,3], 'o')
#plt.show()

with open(argv[1]+'.CASdata','w') as f:
		for freq,x,y,z in CASdata:
				f.write("{} {} {} {}\n".format(freq,x,y,z))

with open(argv[1]+'.CASuni','w') as f:
		for freq,x,y,z in CASuni:
				f.write("{} {} {} {}\n".format(freq,x,y,z))

