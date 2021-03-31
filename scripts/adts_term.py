#!/usr/bin/env python

N=6 # periodicity

fout=open('tune.dat','w')

Qx0=14.0
step=0.001
for i in range(0,3000):
	Qx=Qx0+step*i
	sm=0
	for q in range(-1000,1000):
		if (not (Qx**2-(q*N)**2)==0) and (not ((3*Qx)**2-(q*N)**2)==0):
			sm=sm+1.0/(Qx**2-(q*N)**2)+1.0/((3*Qx)**2-(q*N)**2)
	fout.write(str(Qx)+' '+str(abs(sm))+'\n') 
