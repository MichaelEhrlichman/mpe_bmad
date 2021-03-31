#!/usr/bin/env python

import numpy as np
import sys
import matplotlib.pyplot as pyplot
import csv

#datafile='moga_results.out'
#datafile='p12_771403.out' #'moga_picked.out'
#outfile='p12_771403.at' #'moga_picked.at'
datafile='moga_mcc.out'
outfile='moga_mcc.at'

layout='p12'
at_factors =  (0.5,0.5,0.5,0.5, \
               1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, \
               1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
at_ordering = ('sf[k2]','sd[k2]','shh[k2]','shh2[k2]','b31[k1]','b32[k1]','b3[k1]', \
               'qf1[k1]','qf2[k1]','qf3[k1]','qf4[k1]','qf5[k1]','qf6[k1]', 'qd1[k1]', \
							 'b31[t]','b32[t]','b3[t]', \
							 'qf1[t]','qf2[t]','qf3[t]','qf4[t]','qf5[t]','qf6[t]', 'qd1[t]')
defaults =    {'shh[k2]':'0.0','shh2[k2]':'0.0','b31[t]':'.05817755','b32[t]':'.06362920','b3[t]':'.06567145', \
							 'qf1[t]':'0','qf2[t]':'-0.00159621','qf3[t]':'-0.00010845','qf4[t]':'-0.00749382','qf5[t]':'-0.00749382','qf6[t]':'-0.00749382', 'qd1[t]':'0'}

#defaults =    {'shh[k2]':'0.0','shh2[k2]':'0.0','b31[t]':'0.01163551','b32[t]':'0.01163551','b3[t]':'0.01163551', \
#							 'qf1[t]':'0','qf2[t]':'0','qf3[t]':'0','qf4[t]':'0','qf5[t]':'0','qf6[t]':'0', 'qd1[t]':'0'}

#layout='p6i'
#at_factors =  (0.5,0.5,0.5,0.5,0.5,0.5,\
#               1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0, \
#               1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
#at_ordering = ('sfa[k2]','sda[k2]','sfb[k2]','sdb[k2]','shh[k2]','shh2[k2]', \
#               'b31[k1]','b32[k1]','b3[k1]', 'qf1[k1]','qf2[k1]','qf3[k1]','qf4[k1]','qf5[k1]','qf6[k1]', 'qd1[k1]', \
#							 'b31[t]','b32[t]','b3[t]','qf1[t]','qf2[t]','qf3[t]','qf4[t]','qf5[t]','qf6[t]', 'qd1[t]')
#defaults =    {'b31[t]':'.05817755','b32[t]':'.06362920','b3[t]':'.06567145', \
#							 'qf1[t]':'0','qf2[t]':'-0.00159621','qf3[t]':'-0.00010845','qf4[t]':'-0.00749382','qf5[t]':'-0.00749382','qf6[t]':'-0.00749382', 'qd1[t]':'0'}

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
	print("No generation specified, using last genration by default.")
	GenerationNumber=1

# Read in header
with open(datafile,'r') as f:
	header = f.readline().split()
header.remove('#') #remove hash mark

# Make mapping from Bmad to desired AT ordering.
mapping = [None]*len(at_ordering)
for ix,i in enumerate(at_ordering): 
	if i in header:
		print(ix,header.index(i))
		mapping[ix] = header.index(i)
	elif i in defaults:
		print('{} not in header.  Will apply default.'.format(i))

# Read in data
with open(datafile,'r') as f, open(outfile,'w') as o:
	genCounter=1
	lastStrip=False
	o.write("# Layout "+layout+"\n")
	o.write("# ID "+'   '.join(at_ordering)+'\n')
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
					sline = line.split()
					outline = [str(at_factors[ix]*float(sline[m])) if m is not None else defaults[at_ordering[ix]] for ix,m in enumerate(mapping)]
					o.write(sline[0]+'   '+'   '.join(outline)+'\n')




