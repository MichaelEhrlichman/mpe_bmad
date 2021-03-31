#!/usr/bin/env python

import re, pylab, numpy, math
import analytic_tune_plane_mod as analytic_tune_plane

# SOLEIL
#Qx = 37.92
#Qy = 9.828
#periodicity = 4
#SESAME
# Qx = 24.70 # 24.61
# Qy = 14.358 #14.37
# periodicity = 16
#CHESS
#Qx = 11.285
#Qy = 8.791
#periodicity = 1
#ALS-U AR
Qx = 16.2209
Qy = 8.3285
periodicity = 12

#APS-U
#Qx = 95.1
#Qy = 36.1
#periodicity = 40

# f = open('adts_driver.dat','r')
# adData = []
# for line in f:
# 	adData.append(line.strip().split())
# f.close()
# aadData = numpy.array(adData)
# faadData = aadData.astype(numpy.float)

def additionalItemsToPlot():
	pylab.plot(Qx,Qy,'o') 
	pass;

#--------------------------------------------
#Plot horizontal data, then plot vertical data
#--------------------------------------------

fig = pylab.figure(figsize=(10,10),dpi=125)
pylab.subplots_adjust(left=0.10, right=0.90, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(pylab.FormatStrFormatter('%.2f'))

pylab.gca().set_autoscale_on(False)

(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(Qx,Qy)
Qxmin = 16-0.05
Qxmax = 17.0+0.05
Qymin = 8.0-0.05
Qymax = 9.0+0.05

#Qxmin = 95-0.05
#Qxmax = 96.0+0.05
#Qymin = 36.0-0.05
#Qymax = 37.0+0.05

pylab.xlim(Qxmin,Qxmax)
pylab.ylim(Qymin,Qymax)

pylab.xlabel('Qx')
pylab.ylabel('Qy')

# yarr = numpy.vstack((faadData[:,1],))
# pylab.imshow(yarr,extent=(Qxmin,Qxmax,Qymin,Qymax),cmap='cool',vmin=0,vmax=0.06)

additionalItemsToPlot()  #Plot additional items specified at top of this file
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity)      #Plot an overlay of an analytic tune plane

pylab.title('Tune plane with allowed resonance lines.\n'+atp_str)
pylab.show()






