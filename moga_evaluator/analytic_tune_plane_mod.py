#!/usr/bin/env python

import pylab, math
import numpy as np
from fractions import gcd
from operator import add

#Properties of tune plane
nMax=1000
#pMax = 10
#qMax = 10
#pqMax = 15
pMax = 4
qMax = 4
pqMax = 4

def make_min_max(Qx,Qy):
	intx = int(Qx)
	inty = int(Qy)
	Qxmin = intx - 0.05
	Qxmax = intx + 1.0 + 0.05
	Qymin = inty - 0.05
	Qymax = inty + 1.0 + 0.05
	return (Qxmin,Qxmax,Qymin,Qymax)

def get_ax_size(fig,ax):
	bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	width, height = bbox.width, bbox.height
	width *= fig.dpi
	height *= fig.dpi
	return width, height

def vert_visible_test(p,n,xmin,xmax,ymin,ymax):
	if(p==0):
		return False
	x = float(n)/p
	if( x<(xmax+0.001) and x>(xmin-0.001) ):
		return True
	return False

def visible_test(p,q,n,xmin,xmax,ymin,ymax):
	if(q != 0):
		y1 = (float(n)-p*xmin)/q
		if( y1<(ymax+0.00001) and y1>(ymin-0.00001) ):
			return True
		y2 = (float(n)-p*xmax)/q
		if( y2<(ymax+0.00001) and y2>(ymin-0.00001) ):
			return True
	if(p != 0):
		x1 = (float(n)-q*ymin)/p
		if( x1<(xmax+0.00001) and x1>(xmin-0.00001) ):
			return True
		x2 = (float(n)-q*ymax)/p
		if( x2<(xmax+0.00001) and x2>(xmin-0.00001) ):
			return True
	return False

def analyticTunePlane(fig,ax,Qxmin,Qxmax,Qymin,Qymax,periodicity=1):
	resLines = []
	#-------------------------------
	# Plot low-order non-systematic resonances
	#-------------------------------
	# Plot non-systematic integer resonances (systematic loop below will catch systematic integers)
	intxrange = range(int(math.floor(Qxmin)),int(math.ceil(Qxmax))+1)
	for n in intxrange: #plot non-systematic Qx resonances
		if( n%periodicity != 0 ):
			x1 = n; p = 1; q = 0
			pylab.plot([x1,x1],[Qymin,Qymax], linestyle='dashed', color='red', lw=2.00,scalex=False,scaley=False)
			m = float("inf")  #line is vertical
			b = float("nan")  #no y intercept
			resLines.append([m,b,x1,"(%1i,%1i,%1i,%1i)"%(p,0,0,n)]) # [slope,y-int,x-int] x-int is used for vertical lines
	intyrange = range(int(math.floor(Qymin)),int(math.ceil(Qymax))+1)
	for n in intyrange: #plot non-systematic Qy resonances
		if( n%periodicity != 0 ):
			y1 = n; p = 0; q = 1
			pylab.plot([Qxmin,Qxmax],[y1,y1], linestyle='dashed', color='red', lw=2.00,scalex=False,scaley=False)
			m = 0.0
			b = y1
			resLines.append([m,b,0.0,"(%1i,%1i,%1i,%1i)"%(p,q,0,n)])
	# Plot non-systematic lowest order sum (unbounded) coupling resonances.
	mx, my = np.meshgrid(intxrange, intyrange)
	mx = mx.flatten()
	my = my.flatten()
	cplrange = map(add, mx, my) #element-wise list addition
	cplrange = list(set(cplrange)) #removes duplicate entries from cplrange
	p = 1; q = 1; 
	for n in cplrange:
		if( n%periodicity != 0):
			if( visible_test(p,q,n,Qxmin,Qxmax,Qymin,Qymax) ):
				x1 = Qxmin
				x2 = Qxmax
				y1 = (float(n)-x1)
				y2 = (float(n)-x2)
				#pylab.plot([x1,x2],[y1,y2], linestyle='dashed', color='green', lw=0.25, scalex=False, scaley=False)
				pylab.plot([x1,x2],[y1,y2], linestyle=(0,(5,10)), color='green', lw=0.25, scalex=False, scaley=False)
				m = -float(p)/q
				b = float(n)/q
				resLines.append([m,b,0.0,"(%1i,%1i,%1i,%1i)"%(p,q,0,n)])
	# Plot non-systematic lowest order difference (bounded) coupling resonances.
	cplrange = list(map(add, mx, -my)) + list(map(add, -mx, my)) #element-wise list addition
	cplrange = list(set(cplrange)) #removes duplicate entries from cplrange
	for n in cplrange:
		if( n%periodicity != 0):
			for (p,q) in ((1,-1),(-1,1)):
				if( visible_test(p,q,n,Qxmin,Qxmax,Qymin,Qymax) and \
					 (gcd(gcd(p,q),n) == 1) ): 
					x1 = Qxmin
					x2 = Qxmax
					y1 = (float(n)-p*x1)/q
					y2 = (float(n)-p*x2)/q
					pylab.plot([x1,x2],[y1,y2], linestyle='dashed', color='green', lw=0.25, scalex=False, scaley=False)
					m = -float(p)/q
					b = float(n)/q
					resLines.append([m,b,0.0,"(%1i,%1i,%1i,%1i)"%(p,q,0,n)])
	#-------------------------------
	# Plot systematic resonances
	#-------------------------------
	nrange = [x*periodicity for x in range(-nMax,nMax+1)]
	for n in nrange:
		for p in range(-(pMax),(pMax+1)):
			#Process vertical lines.
			q = 0
			if( (p == 0) or \
					(not vert_visible_test(p,n,Qxmin,Qxmax,Qymin,Qymax)) or \
				  (abs(p) > pqMax) or \
					(gcd(p,n) != 1) ):
				pass;
			else:
				x1 = float(n)/p
				if( abs(p)==1 ):
					pylab.plot([x1,x1],[Qymin,Qymax], linestyle='solid', color='red', lw=2.00,scalex=False,scaley=False)
				elif( abs(p)==2 ):
					pylab.plot([x1,x1],[Qymin,Qymax], linestyle='solid', color='orange', lw=2.00,scalex=False,scaley=False)
				else:
					pylab.plot([x1,x1],[Qymin,Qymax], linestyle='solid', color='0.0', lw=0.25,scalex=False,scaley=False)
				#Save line in y=mx+b format, with m = inf, b = nan, and x = x-intercept
				m = float("inf")  #line is vertical
				b = float("nan")  #no y intercept
				resLines.append([m,b,x1,"(%1i,%1i,%1i,%1i)"%(p,0,0,n)]) # [slope,y-int,x-int] x-int is used for vertical lines

			#Process non-vertical lines
			for q in range(-qMax,qMax+1):
				if( (q==0) or \
						(not visible_test(p,q,n,Qxmin,Qxmax,Qymin,Qymax)) or \
						((abs(p)+abs(q)) > pqMax) or \
					  (gcd(gcd(p,q),n) != 1) ): 
					pass;
				else:
					x1 = Qxmin
					x2 = Qxmax
					y1 = (float(n)-p*x1)/q
					y2 = (float(n)-p*x2)/q
					if(p==0 and abs(q)==1):
						pylab.plot([x1,x2],[y1,y2], linestyle='solid', color='red', lw=2.00,scalex=False,scaley=False)
					elif(p==0 and abs(q)==2):
						pylab.plot([x1,x2],[y1,y2], linestyle='solid', color='orange', lw=2.00,scalex=False,scaley=False)
					elif(abs(p)==1 and abs(q)==1):
						pylab.plot([x1,x2],[y1,y2], linestyle='solid', color='green', lw=0.25,scalex=False,scaley=False)
					else:
						pylab.plot([x1,x2],[y1,y2], linestyle='solid', color='0.0', lw=0.25,scalex=False,scaley=False)
					#Save line in y=mx+b format. third element, x, is used for vertical lines only
					m = -float(p)/q
					b = float(n)/q
					resLines.append([m,b,0.0,"(%1i,%1i,%1i,%1i)"%(p,q,0,n)])

	#Log (but don't plot) window borders:  NEEDED FOR LABEL PLACER
	resLines.append([float('inf'),float('nan'),1.00001*(Qxmin),'lb'])  #left border
	resLines.append([float('inf'),float('nan'),0.99999*(Qxmax),'rb'])  #right border
	resLines.append([0.0,1.00001*(Qymin),0.0,'bb']) #bottom border
	resLines.append([0.0,0.99999*(Qymax),0.0,'tp']) #top border

	#Find point on each line that is furthest from any other lines.
	#This point is where the label will be placed
	for i in range(len(resLines)-4): 
		[m1,b1,x1,labelText] = resLines[i]
		locs = []  #Will be populated with intersections

		otherLines = list(range(len(resLines)))
		otherLines.remove(i)

		#Calculate where line i intersects each other line
		if (m1 != float('inf')):
			for j in otherLines:
				[m2,b2,x2,notUsed] = resLines[j]
				if (m2 != float('inf')):
					#neither m1 nor m2 is vertical
					if b1 != b2:
						if m1 != m2:
							x0 = (b2-b1)/(m1-m2)
							y0 = m1*x0+b1
							#Check if intersection is inside window
							if (x0<Qxmax) & (x0>Qxmin):
								if (y0<Qymax) & (y0>Qymin):
									#Use y-intercept of line i as origin
									loc = np.sqrt(x0*x0+(y0-b1)*(y0-b1))
									locs.append(loc)
				else:
					#m1 is not vertical, m2 is vertical
						x0 = x2
						y0 = m1*x0+b1
						if (x0<Qxmax) & (x0>Qxmin):
							if (y0<Qymax) & (y0>Qymin):
								#Use y-intercept of line i as origin
								loc = np.sqrt(x0*x0+(y0-b1)*(y0-b1))
								locs.append(loc)
		else: #m1 is vertical
			for j in otherLines:
				[m2,b2,x2,notUsed] = resLines[j]
				if (m2 != float('inf')):
					#m1 is vertical, m2 is not vertical
					x0 = x1
					y0 = m2*x0+b2
					if (x0<Qxmax) & (x0>Qxmin):
						if (y0<Qymax) & (y0>Qymin):
							#Use x-intercept of line i as origin
							loc = np.sqrt((x0-x1)*(x0-x1)+y0*y0)
							locs.append(loc)
	  
		if (len(locs) > 0):
			#Determine best location for label as location furthest any other line intersection
			locs.sort() #Contains locations of intersections, relative to y-intercept or x-intercept
			dists = []  #Will contain distances between neighboring points
			for i in range(len(locs)-1):
				dists.append( locs[i+1] - locs[i] )

			idx = dists.index(max(dists))
			labelLoc = locs[idx] + dists[idx]/2.0  #Distance along line i from y-intercept or x-intercept to place label

			if (m1 != float('inf')):
				# line i is not vertical
				theta = math.atan(m1)
				xlabel = labelLoc * math.cos(theta)
				ylabel = labelLoc * math.sin(theta) + b1

				#Transform angle from plot to screen coordinate system
				coords = pylab.array((xlabel,ylabel))
				xsize = ax.get_xlim()[1] - ax.get_xlim()[0]
				ysize = ax.get_ylim()[1] - ax.get_ylim()[0]
				ximageExtent = get_ax_size(fig,ax)[0]
				yimageExtent = get_ax_size(fig,ax)[1]
				trans_angle = math.degrees(math.atan(m1*yimageExtent/ximageExtent*xsize/ysize))
			else:
				# line i is vertical
				trans_angle=90
				xlabel = x1
				ylabel = labelLoc

			#Plot the label
			pylab.text(xlabel,ylabel,labelText,size=8.0,ha='center',va='bottom',rotation_mode='anchor',rotation=trans_angle,color='k')
	atp_str = 'Periodicity: '+str(periodicity)
	return atp_str





