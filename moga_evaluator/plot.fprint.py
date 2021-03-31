#!/usr/bin/env python3

import sys
import matplotlib
matplotlib.use('agg')   #gymnastics for headless ability
import re, pylab, numpy, math
import analytic_tune_plane_mod as analytic_tune_plane

if len(sys.argv) == 2:
  periodicity = int(sys.argv[1])
else:
  periodicity=1
print('Periodicity setting : ', periodicity)

#load ADTS simulation results
try:
	with open("00da_raster/raster_1fprint.dat",'r') as f:
		pData = []
		for line in f:
			if line.strip():  #stop when blank line reached
				if (line.strip().split()[0] != '#'):
					if (float(line.strip().split()[2]) > -99.0):
						pData.append(line.strip().split())
			else:
				break;
		f.close()
		npData = numpy.array(pData)
		fnpData = npData.astype(numpy.float)
except IOError:
	print("Warning: No ADTS footprint data found.")
	fnpData = None

#load chromatic simulation results
f = open('00footprint/footprint.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		pData.append(line.strip().split())
f.close()
npData = numpy.array(pData)
npData[npData == "*"] = "NaN"
cr_fnpData = npData.astype(numpy.float)

#load 1% steps chromatic footprint data
f = open('00footprint/footprint_1pct.dat','r')
pData = []
for line in f:
	if (line.strip().split()[0] != '#'):
		pData.append(line.strip().split())
f.close()
npData_1pct = numpy.array(pData)
npData_1pct[npData_1pct == "*"] = "NaN"
cr_fnpData_1pct = npData_1pct.astype(numpy.float)

sloppymidQx = cr_fnpData[int(len(cr_fnpData[:,1])/2),1]
sloppymidQy = cr_fnpData[int(len(cr_fnpData[:,1])/2),2]
intQx,d = divmod(sloppymidQx,1)
intQy,d = divmod(sloppymidQy,1)
print(f'Found intQx, intQy: {intQx} {intQy}')

def additionalItemsToPlot(ax):
	#Plot ADTS footprint data
	if not fnpData is None:
		sc=pylab.scatter(intQx+fnpData[:,2],intQy+fnpData[:,3],c=numpy.abs(fnpData[:,4]),cmap=pylab.cm.jet,s=2.0,vmin=5,vmax=20)
		pylab.colorbar(sc)
		#for (x,y,z) in fnpData[:,[2,3,4]]:
			#pylab.plot(intQx+x,intQy+y,'.',c='green',markersize=1.0)
			#if z > -50.0:
				#pylab.scatter(intQx+x,intQy+y,c=numpy.abs(z),cmap='jet')
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
fig = pylab.figure(figsize=(9,6),dpi=125)
pylab.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))

pylab.gca().set_autoscale_on(False)

pylab.xlabel('Qx')
pylab.ylabel('Qy')

additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(intQx,intQy)
pylab.plot(41.310,20.325,'o',c='red')
#AR
#Qxmin = 16.0-0.05
#Qxmax = 16.5+0.05
#Qymin = 8.0-0.05
#Qymax = 8.5+0.05
#ALS-U AR
#Qxmin = 38.0-0.05
#Qxmax = 41.0+0.05
#Qymin = 19.0-0.05
#Qymax = 21.0+0.05
Qxmin = 41.0-0.05
Qxmax = 42.0+0.05
Qymin = 20.0-0.05
Qymax = 21.0+0.05
pylab.xlim(Qxmin,Qxmax)
pylab.ylim(Qymin,Qymax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane

pylab.title('Chromatic and ADTS Footprints. '+atp_str)
#pylab.show()
pylab.savefig('adts_footprint.eps',format='eps')






