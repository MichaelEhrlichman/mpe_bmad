#!/usr/bin/env python

import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as pyplot
import csv

datafile='../moga_results.out'

# Count number of generations
f = open(datafile,'r')
genCounter=1
lastStrip=False
for line in f:
	if (line.strip() == ''):
		if (lastStrip):
			genCounter = genCounter + 1
			lastStrip = False
		else:
			lastStrip = True
f.close()
GenerationNumber = genCounter - 1
print("Last Generation: ", genCounter)

if len(sys.argv) > 1:
	GenerationNumber = int(sys.argv[1])
	print("Using genration ", GenerationNumber, " from command line.")
else:
	print("No generation specified, plotting last genration by default.")

# Count number of columns
with open(datafile) as f:
	reader = csv.reader(f,delimiter=' ',skipinitialspace=True)
	first_row=next(reader)
	num_cols = len(first_row) - 1
f.close()

xcol = num_cols-4
ycol = num_cols-3
zcol = num_cols-2

def on_release(event):	
	if event.button==3:
		fig.canvas.draw()
	return

def on_pick(event):
	global dataseed
	if event.mouseevent.button==3:
		#Delete clicked seed from data array.
		seedindex=event.artist.get_label()
		arrayindex=np.where(dataseed==seedindex)
		dataseed=np.delete(dataseed,arrayindex)
		#Remove clicked seed from plot.
		event.artist.remove()
		#fig.canvas.draw()
	if event.mouseevent.button==1:
		print("Seed ID ", event.artist.get_label())

def press(event):
	sys.stdout.flush()
	if event.key=='r':
		print('rescaling y-axis...')
		ax.relim()
		ax.autoscale_view(False,False,True)
		fig.canvas.draw()
	if event.key=='w':
		print('writing moga_picked.out ...')
		fin = open(datafile,'r')
		fout = open('moga_picked.out','w')
		genCounter=1
		lastStrip=False
		for line in fin:
			if (line.strip() == ''):
				if (lastStrip):
					genCounter = genCounter + 1
					lastStrip = False
				else:
					lastStrip = True
			if (genCounter == GenerationNumber):
				if (line.strip() != ''):
					if (line.strip().split()[0] != '#'):
						oneData = line.strip().split()
						if oneData[0] in dataseed:
							fout.write(line)
		print('done!')
		fin.close()
		fout.close()

print("\n\n\n")
print("Left click to show seed IDs of lines near cursor.")
print("\n")
print("Right click to hide lines near.")
print("\n")
print("Press 'r' to rescale y-axis.")
print("\n")
print("Press 'w' to write moga_picked.out file containing all visible lines.")
print("\n\n\n")

# Read in data
f = open(datafile,'r')
pData = []
genCounter=1
lastStrip=False
for line in f:
	if (line.strip() == ''):
		if (lastStrip):
			genCounter = genCounter + 1
			lastStrip = False
		else:
			lastStrip = True
	if (genCounter == GenerationNumber):
		if (line.strip() != ''):
			if (line.strip().split()[0] != '#'):
				pData.append(line.strip().split())
f.close()
npData = np.array(pData)
fnpData = npData.astype(np.float)
dataseed = npData[:,0]
datax = fnpData[:,xcol]
datay = fnpData[:,ycol]
dataz = fnpData[:,zcol]

# Calculate some particulars of the data
axmax = 1.1*np.amax([datax,datay,dataz])
axmin = 0.9*np.amin([datax,datay,dataz])

# Make the axes
fig=pyplot.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])

# Plot the data
indexes = [1,2,3]
datablock = zip(dataseed,datax,datay,dataz)
for datapoint in datablock:
	ax.plot(indexes,datapoint[1:4],picker=4,label=datapoint[0])

# Add decorations to plot
ax.set_xlabel('Energy Offset')
ax.set_ylabel('Objective Value')
ax.set_xlim([0.5,3.5])
#ax.set_ylim([axmin,axmax])
matplotlib.pyplot.autoscale(axis='y')

fig.canvas.mpl_connect('pick_event', on_pick)
fig.canvas.mpl_connect('key_press_event', press)
fig.canvas.mpl_connect('button_release_event', on_release)

pyplot.show()
