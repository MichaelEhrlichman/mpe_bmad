#!/usr/bin/python

import sys
import matplotlib
matplotlib.use('agg')   #gymnastics for headless ability
import re, pylab, numpy, math
import analytic_tune_plane_mod as analytic_tune_plane

if len(sys.argv) == 2:
  periodicity = int(sys.argv[1])
else:
  periodicity=3
print('Periodicity setting : ', periodicity)

#load chromatic simulation results
f = open('00footprint/footprint.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		pData.append(line.strip().split())
f.close()
npData = numpy.array(pData)
npData[npData == "*"] = "NaN"
cr_fnpData1 = npData.astype(numpy.float)

sloppymidQx = cr_fnpData[len(cr_fnpData[:,1])/2,1]
sloppymidQy = cr_fnpData[len(cr_fnpData[:,1])/2,2]
intQx,d = divmod(sloppymidQx,1)
intQy,d = divmod(sloppymidQy,1)

def additionalItemsToPlot(ax):
	#Plot ADTS footprint data
	for (x,y) in fnpData[:,[0,1]]:
		pylab.plot(intQx+x,intQy+y,'.',c='green',markersize=1.0)
	#Plot chromatic footprint data
	pylab.plot(cr_fnpData[:,1],cr_fnpData[:,2],c='blue')
	pylab.plot(cr_fnpData_1pct[:,1],cr_fnpData_1pct[:,2],'o',c='blue')
	for zxy in cr_fnpData_1pct:
		z = zxy[0]*100
		zfrac = z - int(z)
		if zfrac < 0:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(-35,-5),textcoords='offset points',color='red')
		elif zfrac==0:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(-5,-15),textcoords='offset points',color='red')
		else:
			ax.annotate('%.1f'%z, xy = zxy[1:3],xytext=(5,-5),textcoords='offset points',color='red')

#fig = pylab.figure(figsize=(10,10),dpi=125)
fig = pylab.figure(figsize=(6,6),dpi=125)
pylab.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))

pylab.gca().set_autoscale_on(False)

pylab.xlabel('Qx')
pylab.ylabel('Qy')

additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(intQx,intQy)
Qxmin = 37.0-0.05
Qxmax = 38.0+0.05
Qymin = 10.0-0.05
Qymax = 11.0+0.05
pylab.xlim(Qxmin,Qxmax)
pylab.ylim(Qymin,Qymax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane

pylab.title('Chromatic and ADTS Footprints. '+atp_str)
#pylab.show()
pylab.savefig('adts_footprint.eps',format='eps')






