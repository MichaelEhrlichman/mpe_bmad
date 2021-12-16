#!/usr/bin/env python

import numpy as np

ps_CAS2 = []
with open('PSD_1500W_PS_Unipolar_Slow_no_60Hz.txt.CASdata.uniform','r') as f:
	for line in f:
		afreq,acas = map(float, line.split())
		ps_CAS2.append((afreq,(acas/(1e7))**2))

mech_CAS2 = []
with open('accumulated_CAS_at_extraction.dat.uniform','r') as f:
	for line in f:
		afreq,acas = map(float, line.split())
		mech_CAS2.append((afreq,acas**2))

CAS2_lst = []
for psdat,mechdat in zip(ps_CAS2,mech_CAS2):
	CAS2_lst.append((psdat[0],psdat[1]+mechdat[1]))

CAS2 = np.array(CAS2_lst)

with open('merged_CAS_at_extraction.dat','w') as f:
	for freq,x in CAS2:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x**0.5))

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

with open('merged_PSD_at_extraction.dat','w') as f:
	for freq,x in PSD:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x))

with open('filtered_merged_PSD_at_extraction.dat','w') as f:
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

with open('filtered_merged_CAS_at_extraction.dat','w') as f:
	for freq,x in filt_CAS:
		f.write("{:10.4f}   {:13.5e}\n".format(freq,x))









