#!/usr/bin/env python3

import sys
import matplotlib
matplotlib.use('agg')   #gymnastics for headless ability
import re, pylab, numpy, math
import analytic_tune_plane_mod as analytic_tune_plane

if len(sys.argv) == 2:
  periodicity = int(sys.argv[1])
else:
  periodicity=3

#load simulation results
f = open('00adts/tracker_adts.out','r')
pData = []
for line in f:
	if line.strip():  #stop when blank line reached
		if (line.strip().split()[0] != '#'):
			if (float(line.strip().split()[2]) > 0.0):
				pData.append(line.strip().split())
	else:
		break;
f.close()
npData = numpy.array(pData)
fnpData = npData.astype(numpy.float)

sloppymidQx = fnpData[int(len(fnpData[:,1])/2),2]
sloppymidQy = fnpData[int(len(fnpData[:,1])/2),3]

center = int((fnpData[:,1].size-1)/2+1 - 1)

def additionalItemsToPlot(ax):
	pylab.plot(fnpData[:,2],fnpData[:,3])
	ax.annotate('center', xy = fnpData[center,2:4],xytext=(0,0),textcoords='offset points',color='red')
	pylab.plot(fnpData[center,2],fnpData[center,3],'o',c='blue') 
	skipper=6
	for zxy in fnpData:
		z = zxy[0]*1000
		if skipper == 6:
			#ax.annotate('%.1f'%z, xy = zxy[2:4],xytext=(15,-5),textcoords='offset points',color='red')
			pylab.plot(zxy[2],zxy[3],'o',c='blue') 
			skipper = 1
		else:
			skipper = skipper + 1

#--------------------------------------------
#Plot horizontal data, then plot vertical data
#--------------------------------------------

#fig = pylab.figure(figsize=(10,10),dpi=125)
fig = pylab.figure(figsize=(6,6),dpi=125)
pylab.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.1f'))

pylab.gca().set_autoscale_on(False)

pylab.xlabel('Qx')
pylab.ylabel('Qy')

additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(sloppymidQx,sloppymidQy)
#AR
#Qxmin = 16.0-0.05
#Qxmax = 16.5+0.05
#Qymin = 8.0-0.05
#Qymax = 8.5+0.05
#ALS-U AR
Qxmin = 41.0-0.05
Qxmax = 42.0+0.05
Qymin = 20.0-0.05
Qymax = 21.0+0.05
pylab.xlim(Qxmin,Qxmax)
pylab.ylim(Qymin,Qymax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane

pylab.title('ADTS Footprint along +x with low order and\nallowed higher order resonance lines. '+atp_str)
#pylab.show()
pylab.savefig('adts_fp.eps',format='eps')






