#!/usr/bin/env python

import numpy as np
from sys import argv
from scipy import interpolate as interp

CASdata = np.loadtxt(argv[1])

ndata = CASdata.shape[0]

#Assume negative numbers are artifacts and remove them
CASdata[:,1] = np.maximum(CASdata[:,1],[0]*ndata)
CASdata[:,2] = np.maximum(CASdata[:,2],[0]*ndata)
CASdata[:,3] = np.maximum(CASdata[:,3],[0]*ndata)

nuni = 1982
fdelta = 0.05
CASuni = np.zeros((nuni,4))
funi = np.zeros(nuni)

for ix in range(nuni):
	funi[ix] = round( (ix+1)*fdelta, 5)

CASfx = interp.interp1d(CASdata[:,0],CASdata[:,1])
CASfy = interp.interp1d(CASdata[:,0],CASdata[:,2])
CASfz = interp.interp1d(CASdata[:,0],CASdata[:,3])

CASuni[:,0] = funi
CASuni[:,1] = np.maximum(CASfx(funi),[0]*nuni)
CASuni[:,2] = np.maximum(CASfy(funi),[0]*nuni)
CASuni[:,3] = np.maximum(CASfz(funi),[0]*nuni)

with open(argv[1]+'.uniform','w') as f:
		for freq,x,y,z in CASuni:
				f.write(f"{freq:.4f}   {x:13.5e}   {y:13.5e}   {z:13.5e}\n")

