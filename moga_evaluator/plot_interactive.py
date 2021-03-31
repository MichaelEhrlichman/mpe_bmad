#!/usr/bin/env python

import re, pylab, numpy, math
import analytic_tune_plane_mod as analytic_tune_plane

#Properties of simulation data
filename = 'raster_1fprint.dat'
metricMin = -20
metricMax = -1

# AR
intx = 16.
inty = 8.
periodicity = 1
Qxmin = 0.15
Qxmax = 0.5
Qymin = 0.1
Qymax = 0.4
xmax = 0.03 *1e3
ymax = 0.006 *1e3
nxdata = 301
nydata = 301

## SR
#intx = 41.
#inty = 20.
#periodicity = 12
#Qxmin = 0.1
#Qxmax = 0.4
#Qymin = 0.0
#Qymax = 0.4
#xmax = 0.002 *1e3
#ymax = 0.002 *1e3
#nxdata = 301
#nydata = 301

#for x-y plot
dx = (2.0*xmax)/nxdata
dy = ymax/nydata
lx = numpy.linspace(-xmax-dx/2,xmax+dx/2,nxdata+1)
ly = numpy.linspace(ymax/nydata-dy/2,ymax+dy/2,nydata+1)

meshx, meshy = numpy.meshgrid(lx,ly)
gridded = numpy.zeros(nxdata*nydata).reshape(nxdata,nydata)

midx = numpy.linspace(-xmax+dx/2,xmax-dx/2,nxdata)
midy = numpy.linspace(ymax/nydata+dy/2,ymax-dy/2,nydata)
meshmidx, meshmidy = numpy.meshgrid(midx,midy)
flatmidx = meshmidx.flatten('F')
flatmidy = meshmidy.flatten('F')

print()

#for Qx-Qx plot
Qx = numpy.zeros(nxdata*nydata)
Qy = numpy.zeros(nxdata*nydata)
metric = numpy.zeros(nxdata*nydata)

#Read data file
ix = -1 
iy =  0
with open(filename, 'r') as f:
	ii=0
	for line in f:
		fileLine = line.split()
		if len(fileLine) == 0 or fileLine[0][0] == '#': continue
		
		Qx[ii] = float(fileLine[2]) + intx
		Qy[ii] = float(fileLine[3]) + inty
		metric[ii] = float(fileLine[4])

		ix = ix + 1
		if ix > nxdata-1:
			ix = 0
			iy = iy + 1
		gridded[ix,iy] = float(fileLine[4])

		ii = ii + 1

for i in range(len(flatmidx)):
	if Qx[i]<0 or Qy[i]<0:
		flatmidx[i] = float('NaN')
		flatmidy[i] = float('NaN')

#
# Link Plots Class
#
class AnnoteFinder:
	"""
	callback for matplotlib to display an annotation when points are clicked on.	The
	point which is closest to the click and within xtol and ytol is identified.
		
	Register this function like this:
		
	scatter(xdata, ydata)
	af = AnnoteFinder(xdata, ydata, annotes)
	connect('button_press_event', af)
	"""

	def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
		self.data = list(zip(xdata, ydata, annotes))
		if xtol is None:
			xtol = ((max(xdata) - min(xdata))/float(len(xdata)))/2
		if ytol is None:
			ytol = ((max(ydata) - min(ydata))/float(len(ydata)))/2
		self.xtol = xtol
		self.ytol = ytol
		if axis is None:
			self.axis = pylab.gca()
		else:
			self.axis= axis
		self.drawnAnnotations = {}
		self.links = []

	def distance(self, x1, x2, y1, y2):
		"""
		return the distance between two points
		"""
		return math.hypot(x1 - x2, y1 - y2)

	def __call__(self, event):
		if event.inaxes == self.axis:
			if event.button == 1:
				clickX = event.xdata
				clickY = event.ydata
				if self.axis is None or self.axis==event.inaxes:
					annotes = []
					for x,y,a in self.data:
						if not math.isnan(x) and not math.isnan(y):
							if clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
								annotes.append((self.distance(x,clickX,y,clickY),x,y, a) )
					if annotes:
						annotes.sort()
						distance, x, y, annote = annotes[0]
						self.drawAnnote(event.inaxes, x, y, annote)
						for l in self.links:
							l.drawSpecificAnnote(annote)
		if event.button == 3:
			self.clearAnnotes(event.inaxes)

	def clearAnnotes(self, axis):
		for key,val in self.drawnAnnotations.items():
			val[0].remove()
			val[1].remove()
		self.drawnAnnotations = {}
		pylab.gcf().canvas.draw()

	def drawAnnote(self, axis, x, y, annote):
		"""
		Draw the annotation on the plot
		"""
		if (x,y) in self.drawnAnnotations:
			markers = self.drawnAnnotations[(x,y)]
			for m in markers:
				m.set_visible(not m.get_visible())
			self.axis.figure.canvas.draw()
		else:
			t = axis.text(x,y, "({:.2f}, {:.2f}) - {:s}".format(x,y,annote.decode()), color='black', fontsize=15, weight='bold')
			self.axis.set_autoscale_on(False)
			m = axis.scatter([x],[y], marker='d', c='magenta', edgecolors='none', zorder=100)
			self.drawnAnnotations[(x,y)] =(t,m)
			self.axis.figure.canvas.draw()

	def drawSpecificAnnote(self, annote):
		annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
		for x,y,a in annotesToDraw:
			self.drawAnnote(self.axis, x, y, a)

def linkAnnotationFinders(afs):
	for i in range(len(afs)):
		allButSelfAfs = afs[:i]+afs[i+1:]
		afs[i].links.extend(allButSelfAfs)

#--------------------------------------------
#Plot horizontal data, then plot vertical data
#--------------------------------------------

fig = pylab.figure(figsize=(32,14.22),dpi=75)
pylab.subplots_adjust(left=0.05, right=1.00, top=0.9, bottom=0.1)

#make Qx-Qy plot
ax_QxQy = fig.add_subplot(121)
ax_QxQy.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.3f'))

Annotes = numpy.chararray(len(Qx),itemsize=8)
for i in range(len(Qx)):
	Annotes[i] = str(i)

pylab.scatter(Qx,Qy,marker='o',c=metric,s=10, edgecolors='none',cmap='jet')
af1 = AnnoteFinder(Qx,Qy,Annotes,xtol=0.0005,ytol=0.0005,axis=ax_QxQy)
pylab.connect('button_press_event',af1)
pylab.xlabel('Qx')
pylab.ylabel('Qy')
pylab.title('Frequency Map')
pylab.xlim(Qxmin+intx,Qxmax+intx)
pylab.ylim(Qymin+inty,Qymax+inty)
pylab.clim(metricMin, metricMax)
atp_str = analytic_tune_plane.analyticTunePlane(fig,ax_QxQy,Qxmin+intx,Qxmax+intx,Qymin+inty,Qymax+inty,periodicity=periodicity)

#make x-y plot
ax_xy = fig.add_subplot(122)

Zgridded = numpy.ma.masked_where(numpy.isnan(gridded),gridded)
pylab.pcolor(meshx,meshy,Zgridded,cmap='jet')
#print midx
af2 = AnnoteFinder(flatmidx,flatmidy,Annotes,xtol=0.1,ytol=0.1,axis=ax_xy)
pylab.connect('button_press_event',af2)

pylab.colorbar()
pylab.xlabel('x (mm)')
pylab.ylabel('y (mm)')
pylab.title('Frequency Map')
pylab.xlim(-xmax,xmax)
pylab.ylim(0,ymax)
pylab.clim(metricMin, metricMax)

linkAnnotationFinders([af1,af2])

pylab.gca().set_autoscale_on(False)
pylab.show()
