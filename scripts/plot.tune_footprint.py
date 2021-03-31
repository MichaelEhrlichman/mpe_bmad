#!/usr/bin/env python

import sys
import pylab as plt
import numpy as np
import analytic_tune_plane_mod as analytic_tune_plane

plt.rc('text', usetex=True)

#intQx = 16.
#intQy = 12.

intQx = 16.
intQy = 8.

#if len(sys.argv) == 2:
#  periodicity = int(sys.argv[1])
#else:
periodicity=12
print('Periodicity setting : ', periodicity)

#load ADTS simulation results
try:
  with open("raster_1fprint.dat",'r') as f:
    pData = []
    for line in f:
      if line.strip():  #stop when blank line reached
        if (line.strip().split()[0] != '#'):
          pData.append(line.strip().split())
      else:
        break;
    f.close()
    npData = np.array(pData)
    fnpData = npData.astype(np.float)
except IOError:
  fnpData = None

def additionalItemsToPlot(ax):
  #Plot ADTS footprint data
  if not fnpData is None:
    x = fnpData[:,[2]]
    y = fnpData[:,[3]]
    z = fnpData[:,[4]]
    ix = z > -99.0
    sc = plt.scatter(intQx+x[ix], intQy+y[ix], c=z[ix], s=1, cmap='jet')
    plt.colorbar(sc, label=r'$\log\sqrt{\Delta\nu_x^2+\Delta\nu_y^2}$')
  else:
    print("Warning: No ADTS footprint data found.")

#fig = plt.figure(figsize=(10,10),dpi=125)
fig = plt.figure(figsize=(6,6),dpi=125)
plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

plt.gca().set_autoscale_on(False)

plt.xlabel('Qx')
plt.ylabel('Qy')

additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

(Qxmin,Qxmax,Qymin,Qymax) = analytic_tune_plane.make_min_max(intQx,intQy)
#Qxmin = 16.50 #16.52 #16.0  
#Qxmax = 16.60 #16.56 #16.5  
#Qymin = 12.60 #12.61 #8.0   
#Qymax = 12.65 #12.63 #8.5   
Qxmin = 16.0 - 0.05  
Qxmax = 16.5 + 0.05
Qymin = 8.0 - 0.05
Qymax = 8.5 + 0.05
plt.xlim(Qxmin,Qxmax)
plt.ylim(Qymin,Qymax)
plt.clim(-2.0,-25.0)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=periodicity)      #Plot an overlay of an analytic tune plane

plt.title('Chromatic and ADTS Footprints. '+atp_str)
#plt.show()
plt.savefig('adts_footprint.eps',format='eps')






