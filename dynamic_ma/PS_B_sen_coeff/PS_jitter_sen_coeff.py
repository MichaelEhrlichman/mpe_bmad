#!/usr/bin/env python3

import numpy as np
from sys import argv

CASdata = np.loadtxt(argv[1])
nCASdata = CASdata.shape[0]

senData = np.loadtxt(argv[2])
nSenData = senData.shape[0]

print("nCASdata:  {}\nSenData:   {}\n".format(nCASdata,nSenData))

with open(argv[1]+'.XSPECTRUM','w') as f:
	i=0
	cas = CASdata[0,1]
	print(cas)
	for ix,[loc,senQx,senQy,eleix] in enumerate(senData):
		Qx_sigma = (cas * 1.0e-6 / 10.0 * senQx)
		Qy_sigma = 0
		f.write("{}   {}   {}   {}\n".format(loc,ix,Qx_sigma,Qy_sigma))
