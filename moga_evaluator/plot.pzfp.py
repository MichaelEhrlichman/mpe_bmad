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

#load data
f = open('00footprint/footprint.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		pData.append(line.strip().split())
f.close()
npData = numpy.array(pData)
npData[npData == "*"] = "NaN"
fnpData = npData.astype(numpy.float)

#load 1% steps chromatic footprint data
f = open('00footprint/footprint_1pct.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		pData.append(line.strip().split())
f.close()
npData_1pct = numpy.array(pData)
npData_1pct[npData_1pct == "*"] = "NaN"
fnpData_1pct = npData_1pct.astype(numpy.float)

def additionalItemsToPlot(ax):
	pylab.plot(fnpData[:,1],fnpData[:,2],c='blue')
	pylab.plot(fnpData_1pct[:,1],fnpData_1pct[:,2],'o',c='blue')
	#for zxy in [fnpData_1pct[0], fnpData_1pct[-1]]:
	for zxy in fnpData_1pct:
		z = zxy[0]*100
		zfrac = z - int(z)
		if zfrac < 0:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(-35,-5),textcoords='offset points',color='red')
		elif zfrac==0:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(-5,-15),textcoords='offset points',color='red')
		else:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(5,-5),textcoords='offset points',color='red')

#--------------------------------------------
# Plot full integer
#--------------------------------------------

fig = pylab.figure(figsize=(6,6),dpi=125)
pylab.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))

pylab.gca().set_autoscale_on(False)

pylab.xlabel('Qx')
pylab.ylabel('Qy')

sloppymidQx = fnpData[int(len(fnpData[:,1])/2),1]
sloppymidQy = fnpData[int(len(fnpData[:,1])/2),2]
if numpy.isnan(sloppymidQx) or numpy.isnan(sloppymidQy):
	print("Isnan detected in on-energy fnpData.  Broken.")
	quit()
(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(sloppymidQx,sloppymidQy)
pylab.xlim(Qxmin,Qxmax)
pylab.ylim(Qymin,Qymax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane
additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

#pylab.title('Chromatic Footprint with low order and\nallowed higher order resonance lines. '+atp_str)
pylab.title('Chromatic footprint with systematic resonances\n'+atp_str)
#pylab.show()
pylab.savefig('footprint.eps',format='eps')

pylab.clf()
#--------------------------------------------
# Plot zoom 
#--------------------------------------------

# fig = pylab.figure(figsize=(6,6),dpi=125)
# pylab.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
# ax = fig.add_subplot(111)
# ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
# ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
# 
# pylab.gca().set_autoscale_on(False)
# 
# pylab.xlabel('Qx')
# pylab.ylabel('Qy')
# 
# Qxmin = 17.5 #numpy.amin(fnpData[:,1])-0.1
# Qxmax = 17.8 #numpy.amax(fnpData[:,1])+0.1
# Qymin = 13.5 #numpy.amin(fnpData[:,2])-0.1
# Qymax = 13.8 #numpy.amax(fnpData[:,2])+0.1
# pylab.xlim(Qxmin,Qxmax)
# pylab.ylim(Qymin,Qymax)
# print(Qxmin, Qxmax, Qymin, Qymax)
# atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane
# additionalItemsToPlot(ax)  #Plot additional items specified at top of this file
# 
# pylab.title('Chromatic footprint with systematic resonances\n'+atp_str)
# pylab.savefig('footprint_zoom.eps',format='eps')



