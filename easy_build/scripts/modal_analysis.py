#!/usr/bin/env python3

import numpy as np
import scipy.fft as fft
import scipy.optimize as opt
import pylab

#nturns = 84764+1
nturns=5000

# #horizontal
nbunches = 26
#aix = 2
#bix = 3
#vertical
aix = 4
bix = 5

print("!!!!!!!!!!!! nbunches is set to {}".format(nbunches))

data = np.zeros((nturns,nbunches),dtype=np.complex128)

with open('~bunch_vs_number.env','r') as datafile:
	turn=0
	bunch=0
	for line in datafile:
		if turn == nturns:
			break
		if line.strip() and not line.strip().startswith("#"):
			a = float(line.split()[aix])
			b = float(line.split()[bix])
			data[turn,bunch] = complex(a,-b)
			bunch = bunch + 1
			if bunch == nbunches:
				data[turn,:] = data[turn,:] - complex(np.mean(np.real(data[turn,:])), np.mean(np.imag(data[turn,:])))
				bunch = 0
				turn = turn+1

datafft = np.zeros_like(data)

turn =0
for i in range(nturns):
	datafft[i,:] = fft(data[i,:])

#datafft  looks like complex_fft[turn,mode]

with open('~bunch_vs_number.fft_by_bunch','w') as outfile:
	outfile.write('turn mode Re(F) Im(F)')
	for i in range(nturns):
		for j in range(nbunches):
			outfile.write('{}   {}   {:14.4e}   {:14.4e}\n'.format(i,j,np.real(datafft[i,j]),np.imag(datafft[i,j])))
		outfile.write('\n\n')

with open('~bunch_vs_number.fft_by_mode','w') as outfile:
	outfile.write('turn mode Re(F) Im(F)')
	for j in range(nbunches):
		for i in range(nturns):
			outfile.write('{}   {}   {:14.4e}   {:14.4e}\n'.format(i,j,np.real(datafft[i,j]),np.imag(datafft[i,j])))
		outfile.write('\n\n')

x0 = 2999
x1 = 4999
ra=list(range(x0,x1))
mode = 1

def g(x,g0,m):
	return g0 + m*x - m*x0
def g_wrap(args):
	return abs(datafft[ra,mode]) - [g(x,args[0],args[1]) for x in ra]

A0=0.001
def f(x,b,A,t):
	return b-A*(1.0-np.exp((x-x0)/t))
def f_wrap(args):
	return abs(datafft[ra,mode]) - [f(x,args[0],args[1],args[2]) for x in ra]

def fr(x,b,A,r):
	return b-A*(1.0-np.exp((x-x0)*r))
def fr_wrap(args):
	return abs(datafft[ra,mode]) - [fr(x,args[0],args[1],args[2]) for x in ra]

def fralt(x,b,A,er):
	return b-A*(1.0-er**(x-x0))
def fralt_wrap(args):
	return abs(datafft[ra,mode]) - [fralt(x,args[0],args[1],args[2]) for x in ra]

rate=np.zeros(nbunches)
with open('rate_vs_mode.out','w') as ratefile:
	for mode in range(1,nbunches):
		gopt = opt.least_squares(g_wrap,[0.00045,1.0],method='lm')
		#fraltopt = opt.least_squares(fralt_wrap,[gopt.x[0],A0,np.exp(gopt.x[1]/A0)],method='lm',max_nfev=10000,x_scale=[1.0e6,1.0e6,1.0])
		fopt = opt.least_squares(f_wrap,[gopt.x[0],A0,A0/gopt.x[1]],method='lm',max_nfev=10000)
		#rate[mode]=np.log(fraltopt.x[2])
		rate[mode]=1.0/fopt.x[2]
		ratefile.write('{:5d}   {:15.6e}\n'.format(mode, rate[mode]))

#mode=250
#print(tau[mode])
#fdata = np.zeros_like(ra,dtype=np.float)
#gopt, gcov = opt.curve_fit(g,ra,abs(datafft[ra,mode]),[0.00045,1.0])
#fopt, fcov = opt.curve_fit(f,ra,abs(datafft[ra,mode]),[gopt[0],0.0008,A0/gopt[1]])
#print(fopt[0],fopt[1],fopt[2])
#for i in range(len(ra)):
#	fdata[i] = f(ra[i],fopt[0],fopt[1],fopt[2])
#	print(ra[i], fdata[i])
#pylab.plot(ra,fdata)
#pylab.plot(ra,abs(datafft[ra,mode]))
#pylab.show()
