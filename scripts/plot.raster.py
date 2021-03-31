#!/usr/bin/env python

import sys
import numpy as np
import pylab as plt

plt.rc('text', usetex=True)

#load ADTS simulation results
try:
  with open("raster_1fprint.dat",'r') as f:
    pData = []
    for line in f:
      if line.strip():  #stop when blank line reached
        if (line.strip().split()[0] != '#'):
          pData.append(line.strip().split())
        else:
          header_line = line.strip().split()
          nx = int(header_line[2])
          ny = int(header_line[4])
      else:
        break;
    f.close()
    npData = np.array(pData)
    fnpData = npData.astype(np.float)
except IOError:
  fnpData = None

#make mesh
datMesh = np.zeros((nx,ny))
k=0
for i in range(nx):
  for j in range(ny):
    datMesh[i,j] = fnpData[k,4]
    k += 1
masked_datMesh = np.ma.masked_where(datMesh<-99.0,datMesh)
xmin = fnpData[0,0]
xmax = fnpData[-1,0]
ymin = fnpData[0,1]
ymax = fnpData[ny-1,1]
datX = np.zeros((nx+1,ny+1))
datY = np.zeros((nx+1,ny+1))
dx = (xmax-xmin)/float(nx-1)
dy = (ymax-ymin)/float(ny-1)
for i in range(nx+1):
  for j in range(ny+1):
    datX[i,j] = xmin-(dx/2.0)+i*dx
    datY[i,j] = ymin-(dy/2.0)+j*dy

def additionalItemsToPlot(ax):
  #Plot ADTS footprint data
  if not fnpData is None:
    sc = ax.pcolormesh(datX,datY,masked_datMesh,cmap='jet')
    plt.colorbar(sc, label=r'$\log\sqrt{\Delta\nu_x^2+\Delta\nu_y^2}$')
    #x = fnpData[:,[0]]
    #y = fnpData[:,[1]]
    #z = fnpData[:,[4]]
    #ix = z>-99.0
    #plt.scatter(x[ix], y[ix], c=z[ix], cmap='jet', verts=(verts), s=47.0, antialiased=None, edgecolors='none')
  else:
    print("Warning: No ADTS footprint data found.")


#fig = plt.figure(figsize=(10,10),dpi=125)
fig = plt.figure(figsize=(6,6),dpi=125)
plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
ax = fig.add_subplot(111)
ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.3f'))

plt.gca().set_autoscale_on(True)

plt.xlabel('x (m)')
plt.ylabel('y (m)')

additionalItemsToPlot(ax)  #Plot additional items specified at top of this file

#plt.xlim(-0.010,0.010)
#plt.ylim(0.0,0.010)

plt.title('Raster with Frequency Map Coloring')
#plt.show()
#plt.savefig('raster.svg')
plt.savefig('raster.eps',format='eps')






