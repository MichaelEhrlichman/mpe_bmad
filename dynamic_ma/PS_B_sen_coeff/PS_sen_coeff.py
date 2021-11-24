#!/usr/bin/env python3

import numpy as np
from sys import argv

CASdata = np.loadtxt(argv[1])
nCASdata = CASdata.shape[0]

senData = np.loadtxt(argv[2])
nSenData = senData.shape[0]

print("nCASdata:  {}\nSenData:   {}\n".format(nCASdata,nSenData))

with open(argv[1]+'.QSPECTRUM','w') as f:
	i=0
	for freq,cas in CASdata:
		i += 1
		if not i%10:
			print("processing {} of {}".format(i,nCASdata))
		Qx_sigma2 = 0.0
		Qy_sigma2 = 0.0
		Qx_direct = 0.0
		Qy_direct = 0.0
		lasteleix = -1
		for loc,length,senQx,senQy,eleix in senData:
			if eleix != -1:
				if eleix != lasteleix:
					lasteleix = eleix
					Qx_sigma2 = Qx_sigma2 + Qx_direct**2
					Qy_sigma2 = Qy_sigma2 + Qy_direct**2
					Qx_direct = 0.0
					Qy_direct = 0.0
				Qx_direct = Qx_direct + (cas * 1.0e-6 / 10.0 * senQx)
				Qy_direct = Qy_direct + (cas * 1.0e-6 / 10.0 * senQy)
		f.write("{}   {}   {}\n".format(freq,np.sqrt(Qx_sigma2),np.sqrt(Qy_sigma2)))
		#print("{}   {}   {}\n".format(freq,np.sqrt(Qx_sigma2),np.sqrt(Qy_sigma2)))
