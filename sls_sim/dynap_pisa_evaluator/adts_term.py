#!/opt/python/python-2.7.5/bin/python

N=12 # periodicity

fout=open('tune.dat','w')

Qx0=36.0
step=0.01
for i in range(0,300):
	Qx=Qx0+step*i
	sm=0
	for q in range(-1000,1000):
		if (not (Qx**2-(q*N)**2)==0) and (not ((3*Qx)**2-(q*N)**2)==0):
			sm=sm+1.0/(Qx**2-(q*N)**2)+1.0/((3*Qx)**2-(q*N)**2)
	fout.write(str(Qx)+' '+str(abs(sm))+'\n') 

