#!/usr/bin/python

import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update({'font.size':26})

tlIdeal = []
tl = []

f=open('../00touschek/tl.dat')
lines=f.readlines()
tlIdeal = float(lines[1])
f.close()

for filename in glob.iglob('seed_*/00touschek/tl.dat'):
	f=open(filename)
	lines=f.readlines()
	tl.append(float(lines[1]))
	f.close()

print len(tl)

tl = sorted(tl)
i=0
for item in tl:
	i = i + 1
	print i, item

plt.hist(tl,bins=9)
plt.title('Touschek Lifetimes of MA & Cor. Lattices')
plt.xlabel('Lifetime (hr)')
plt.ylabel('Frequency')
plt.plot([tlIdeal,tlIdeal],[0,1],linewidth=4)
plt.annotate('ideal',xy=(tlIdeal,1),xytext=(tlIdeal+0.01,1.05))
plt.xlim((3.5,tlIdeal+0.1))
plt.ylim((0,6))
plt.show()


#fig = plt.gcf()
