#!/usr/bin/env python

import numpy as np
from sys import argv
from scipy import interpolate as interp

CASdata = np.loadtxt(argv[1])

ndata = CASdata.shape[0]

#Assume negative numbers are artifacts and remove them
for ix in range(1,7):
	CASdata[:,ix] = np.maximum(CASdata[:,ix],[0]*ndata)

nuni = 24999
fdelta = 0.1
CASuni = np.zeros((nuni,6))
funi = np.zeros(nuni)

for ix in range(nuni):
	funi[ix] = round( (ix+1)*fdelta, 3)

CASfx = interp.interp1d(CASdata[:,0],CASdata[:,1])
CASfxr = interp.interp1d(CASdata[:,0],CASdata[:,2])
CASfy = interp.interp1d(CASdata[:,0],CASdata[:,3])
CASfyr = interp.interp1d(CASdata[:,0],CASdata[:,4])
CASfz = interp.interp1d(CASdata[:,0],CASdata[:,5])
CASfzr = interp.interp1d(CASdata[:,0],CASdata[:,6])

CASuni[:,0] = funi
CASuni[:,1] = np.maximum(CASfx(funi),[0]*nuni)
CASuni[:,2] = np.maximum(CASfxr(funi),[0]*nuni)
CASuni[:,3] = np.maximum(CASfy(funi),[0]*nuni)
CASuni[:,4] = np.maximum(CASfyr(funi),[0]*nuni)
CASuni[:,5] = np.maximum(CASfz(funi),[0]*nuni)
CASuni[:,6] = np.maximum(CASfzr(funi),[0]*nuni)

with open(argv[1]+'.uniform','w') as f:
		for freq,x,xr,y,yr,z,zr in CASuni:
				f.write(f"{freq:.4f}   {x:13.5e}   {xr:13.5e}   {y:13.5e}   {yr:13.5e}   {z:13.5e}   {zr:13.5e}\n")

