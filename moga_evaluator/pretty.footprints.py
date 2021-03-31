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

#load chromatic footprint data
f = open('00footprint/footprint.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		if (line.strip().split()[1] != '*'):
			if (line.strip().split()[2] != '*'):
				pData.append(line.strip().split())
f.close()
npData = numpy.array(pData)
fnpData = npData.astype(numpy.float)

sloppymidQx = fnpData[int(len(fnpData[:,1])/2),1]
sloppymidQy = fnpData[int(len(fnpData[:,1])/2),2]

#load 1% steps chromatic footprint data
f = open('00footprint/footprint_1pct.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		if (line.strip().split()[1] != '*'):
			if (line.strip().split()[2] != '*'):
				pData.append(line.strip().split())
f.close()
npData_1pct = numpy.array(pData)
fnpData_1pct = npData_1pct.astype(numpy.float)

#load adts footprint results
f = open('00adts/tracker_adts.out','r')
adts_pData = []
cross_adts_pData = []
ixer = 1
inBetween = False
for line in f:
	if line.strip():  # if line is not empty
		inBetween = False
		if (line.strip().split()[0] != '#'):
			if (float(line.strip().split()[2]) > 0.0):
				if ixer == 1:
					adts_pData.append(line.strip().split())
				if ixer == 3:
					cross_adts_pData.append(line.strip().split())
	else:
		if not inBetween:
			ixer = ixer + 1
			inBetween = True
f.close()
adts_npData = numpy.array(adts_pData)
adts_fnpData = adts_npData.astype(numpy.float)
cross_adts_npData = numpy.array(cross_adts_pData)
cross_adts_fnpData = cross_adts_npData.astype(numpy.float)

def additionalItemsToPlot(ax):
	pylab.plot(fnpData[:,1],fnpData[:,2],c='blue',label='Chromatic')
	pylab.plot(fnpData_1pct[:,1],fnpData_1pct[:,2],'o',c='blue')
#	for zxy in fnpData:
#		pylab.plot(zxy[1],zxy[2],'o',c='blue') 
	for zxy in [fnpData[0], fnpData[-1]]:
		z = zxy[0]*100
		zfrac = z - int(z)
		if z < 0:
			ax.annotate('%+d%%'%z, xy = zxy[1:3],xytext=(-30,-5),textcoords='offset points',color='red')
		elif z==0:
			ax.annotate('%+d%%'%z, xy = zxy[1:3],xytext=(-5,-15),textcoords='offset points',color='red')
		else:
			ax.annotate('%+d%%'%z, xy = zxy[1:3],xytext=(3,2),textcoords='offset points',color='red')
	pylab.plot(adts_fnpData[:,2],adts_fnpData[:,3],c='green',label=r'ADTS along $\pm$x')
	#pylab.plot(adts_fnpData[-1,2],adts_fnpData[-1,3],'o',c='green') 
	#pylab.plot(cross_adts_fnpData[:,2],cross_adts_fnpData[:,3])

#--------------------------------------------
#Plot horizontal data, then plot vertical data
#--------------------------------------------

#fig = pylab.figure(figsize=(10,10),dpi=125)
fig = pylab.figure(figsize=(6,6),dpi=125)
pylab.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))

pylab.gca().set_autoscale_on(False)

#pylab.xlabel(r'Q$_\mathrm{x}$')
#pylab.ylabel(r'Q$_\mathrm{y}$')
pylab.xlabel(r'$\mathregular{Q_x}$')
pylab.ylabel(r'$\mathregular{Q_y}$')

sloppymidQx = fnpData[int(len(fnpData[:,1])/2),1]
sloppymidQy = fnpData[int(len(fnpData[:,1])/2),2]
(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(sloppymidQx,sloppymidQy)
# CANDLE
#Qxmin = 24.45
#Qxmax = 24.90
#Qymin = 14.10
# dc01a - full int
#Qxmin = 36.95
#Qxmax = 38.05
#Qymin = 9.95
#Qymax = 11.05
# CESR
#Qxmin = 17.5
#Qxmax = 17.8
#Qymin = 13.5
#Qymax = 13.8
# APS-U
# Qxmin = 94.85
# Qxmax = 96.15
# Qymin = 35.85
# Qymax = 37.15
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
print(Qxmin, Qxmax, Qymin, Qymax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane
additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

#legCro = matplotlib.lines.Line2D([],[],color='blue',marker='o',label='Chromatic Footprint')
#matplotlib.pyplot.legend()

pylab.title('Tune Footprints and Systematic Resonances\nPeriodicity = '+str(periodicity))
#pylab.show()
pylab.savefig('pretty.footprints.eps',format='eps')






