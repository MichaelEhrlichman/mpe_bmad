#!/usr/bin/python

import numpy as np
import convex_hull as ch

n = sum(1 for line in open('beats_metric.dat','r'))
print "Number of data points: ", n

data = np.zeros((n,5))

n50 = n//2-1
n90 = n//100*90-1

f = open('beats_metric.dat','r')
i = 0
for line in f:
	data[i,0:4] = line.split()
	i = i + 1
f.close()

mean_beta_x = data[:,0].mean()
min_beta_x = np.amin(data[:,0])
max_beta_x = np.amax(data[:,0])

mean_beta_y = data[:,1].mean()

mean_area = data[:,3].mean()
min_area = np.amin(data[:,3])
max_area = np.amax(data[:,3])

print "Average x beta beat: ", mean_beta_x
print "Average y beta beat: ", mean_beta_y
print "Average DA area: ", mean_area

for i in range(0,n):
	data[i,4] = ((data[i,0]-mean_beta_x)/(max_beta_x-min_beta_x))**2 + ((data[i,3]-mean_area)/(max_area-min_area))**2

sorted = data[data[:,4].argsort()]

sorted_50 = sorted[0:n50]
sorted_90 = sorted[0:n90]

#50
data_alt = np.zeros((2,n50))
data_alt[0,:] = sorted_50[:,0]
data_alt[1,:] = sorted_50[:,3]

hull_pts = ch.convex_hull(data_alt)

f = open('50_pct_hull.dat','w')
for xy in hull_pts:
 	f.write(str(xy[0])+"   "+str(xy[1])+"\n")
xy = hull_pts[0]
f.write(str(xy[0])+"   "+str(xy[1])+"\n")
f.close()

#90
data_alt = np.zeros((2,n90))
data_alt[0,:] = sorted_90[:,0]
data_alt[1,:] = sorted_90[:,3]

hull_pts = ch.convex_hull(data_alt)

f = open('90_pct_hull.dat','w')
for xy in hull_pts:
 	f.write(str(xy[0])+"   "+str(xy[1])+"\n")
xy = hull_pts[0]
f.write(str(xy[0])+"   "+str(xy[1])+"\n")
f.close()

print mean_beta_x, mean_beta_y, mean_area
