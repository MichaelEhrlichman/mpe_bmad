#!/usr/bin/env python3

import numpy as np
import sys

ele_ix_str = sys.argv[1].strip()

amp_indexes = []
amp_factors = []
with open('../amp_factors/amp_factors_ele'+ele_ix_str+'.dat','r') as f:
	for l in f:
		if not l.strip().startswith("#"):
			s = [x.strip() for x in l.split()]
			amp_indexes.append(s[0:1][0])
			amp_factors.append(s[1:2][0])  #x
			#amp_factors.append(s[2:3][0])  #y

keys = []
with open('../key.dat','r') as f: 
	for l in f:
		if not l.strip().startswith("#"):
			s = [[x.strip() for x in l.split()][i] for i in (0,2,3)]
			keys.append(s)

first = True
for ix,fname,col in keys:
	if fname != '*':
		CASdata = np.loadtxt(fname,usecols=(0,int(col)))
		if first:
			CAS2 = np.zeros_like(CASdata)  # initialize data structure
			CAS2[:,0] = CASdata[:,0]  # copy over indices
			first = False
		amp_factor = float(amp_factors[amp_indexes.index(ix)])
		CAS2[:,1] = CAS2[:,1] + np.square(amp_factor*CASdata[:,1])

with open('accumulated_CAS_at_extraction_ele'+ele_ix_str+'.dat','w') as f:
	for freq,x in CAS2:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,np.sqrt(x)))

def PSD_attenuation_ratio(f):
	fbw = 200.0
	return (f/fbw)**2 / (1.0+(f/fbw)**2)

# PSD = -d(CAS**2)/df
C0x = CAS2[-1,1]
PSD = np.zeros_like(CAS2)
PSD[:,0] = CAS2[:,0]
filt_PSD = np.zeros_like(CAS2)
filt_PSD[:,0] = CAS2[:,0]
ndata = CAS2.shape[0]
for ix in range(0,ndata-1):
	filt = PSD_attenuation_ratio( CAS2[ix,0] )
	PSD[ix,1] = -( (CAS2[ix+1,1]-CAS2[ix,1]) / (CAS2[ix+1,0]-CAS2[ix,0]) )
	filt_PSD[ix,1] = filt*PSD[ix,1]
PSD[ndata-1,1] = -( (CAS2[ndata-1,1]-CAS2[ndata-2,1]) / (CAS2[ndata-1,0]-CAS2[ndata-2,0]) )
filt = PSD_attenuation_ratio( (CAS2[ndata-1,0]+CAS2[ndata-2,0])/2.0 )
filt_PSD[ndata-1,1] = filt*PSD[ndata-1,1]
filt_C0x = C0x * PSD_attenuation_ratio( PSD[-1,0] )

with open('PSD_at_extraction_ele'+ele_ix_str+'.dat','w') as f:
	for freq,x in PSD:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x))

with open('filtered_PSD_at_extraction_ele'+ele_ix_str+'.dat','w') as f:
	for freq,x in filt_PSD:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x))

filt_CAS = np.zeros_like(filt_PSD)
filt_CAS[:,0] = filt_PSD[:,0]
ndata = filt_CAS.shape[0]
for ix in range(ndata):
	filt_CAS[ix,1] = np.sqrt(np.trapz(filt_PSD[ix:,1],filt_PSD[ix:,0])+filt_C0x) 
#filt_CAS[ndata-1,1] = np.sqrt(0.05*filt_PSD[ndata-1,1])
#for ix in range(ndata-2,-1,-1):
#	filt_CAS[ix,1] = np.sqrt(filt_CAS[ix+1,1]**2 + 0.05*filt_PSD[ix,1])	

with open('filtered_CAS_at_extraction_ele'+ele_ix_str+'.dat','w') as f:
	for freq,x in filt_CAS:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x))









