#!/usr/bin/env python

import sys
import numpy as np

datafile = sys.argv[1]
column = int(sys.argv[2])-1

data = []
with open(datafile,'r') as file:
  for line in file:
    data.append(float(line.split()[column]))

print("Mean: {:9.3f}".format(np.mean(data)))
print("Standard Deviation: {:9.3f}".format(np.std(data)))
