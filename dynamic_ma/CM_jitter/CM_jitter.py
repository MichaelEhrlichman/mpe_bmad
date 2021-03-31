#!/usr/bin/env python3

import numpy as np
import sys
import glob

jitterData_x = []
nmeas=0
with open('jitter_files_x.dat','r') as f:
	for jitterfile in f:
		nmeas=nmeas+1
		jitterData_x.append(np.loadtxt(jitterfile.strip())[:,0:2])

jitterData_y = []
nmeas=0
with open('jitter_files_y.dat','r') as f:
	for jitterfile in f:
		nmeas=nmeas+1
		jitterData_y.append(np.loadtxt(jitterfile.strip())[:,0:2])

nspectra = jitterData_x[1].shape[0]

orm_Ap_h = np.loadtxt('orm_Ap_h.dat')
orm_Ap_v = np.loadtxt('orm_Ap_v.dat')

# jitterData[BPMs][spectral component,[f,A]

orbit_vector = np.zeros(nmeas)
cmData_x = nspectra*[np.zeros(nmeas)]
for i in range(nspectra):
	for j in range(nmeas):
		orbit_vector[j] = jitterData_x[j][i,1]
		cmData_x[i] = orm_Ap_h.dot(orbit_vector)

cmData_y = nspectra*[np.zeros(nmeas)]
for i in range(nspectra):
	for j in range(nmeas):
		orbit_vector[j] = jitterData_y[j][i,1]
		cmData_y[i] = orm_Ap_v.dot(orbit_vector)

CAS2_x = np.square(cmData_x)
CAS2_y = np.square(cmData_y)
PSD_x = np.zeros_like(CAS2_x)
PSD_y = np.zeros_like(CAS2_y)
for ix in range(0,nspectra-1):
	PSD_x[ix] = -( (CAS2_x[ix+1]-CAS2_x[ix]) / (jitterData_x[0][ix+1,0]-jitterData_x[0][ix,0]) )
	PSD_y[ix] = -( (CAS2_y[ix+1]-CAS2_y[ix]) / (jitterData_y[0][ix+1,0]-jitterData_y[0][ix,0]) )
PSD_x[nspectra-1] = -( (CAS2_x[nspectra-1]-CAS2_x[nspectra-2]) / (jitterData_x[0][nspectra-1,0]-jitterData_x[0][nspectra-2,0]) )
PSD_y[nspectra-1] = -( (CAS2_y[nspectra-1]-CAS2_y[nspectra-2]) / (jitterData_y[0][nspectra-1,0]-jitterData_y[0][nspectra-2,0]) )

with open('CM_CAS_x.out','w') as f:
	fmt_str = ''.join(['{}   ']*nmeas)
	for i in range(nspectra):
		f.write(fmt_str.format(jitterData_x[0][i,0],*cmData_x[i][:])+'\n')

with open('CM_CAS_y.out','w') as f:
	fmt_str = ''.join(['{}   ']*nmeas)
	for i in range(nspectra):
		f.write(fmt_str.format(jitterData_y[0][i,0],*cmData_y[i][:])+'\n')

with open('CM_PSD_x.out','w') as f:
	fmt_str = ''.join(['{}   ']*nmeas)
	for i in range(nspectra):
		f.write(fmt_str.format(jitterData_x[0][i,0],*PSD_x[i][:])+'\n')

with open('CM_PSD_y.out','w') as f:
	fmt_str = ''.join(['{}   ']*nmeas)
	for i in range(nspectra):
		f.write(fmt_str.format(jitterData_y[0][i,0],*PSD_y[i][:])+'\n')


