#!/bin/env python3

import sys
import numpy as np
from numpy.linalg import inv

rtape_file = sys.argv[1]

with open(rtape_file,'r') as f:
  lines = f.readlines()

nhead = 2
recsize = 8
nrec = (len(lines)-1.0-nhead)/recsize
assert int(nrec) - nrec == 0, "Math?  What math?"
nrec=int(nrec)

def chunk(string, length=16):
    return [string[0+i:length+i] for i in range(0, len(string), length)]

def chop(a,eps=1e-9):
  a[np.abs(a) < eps] = 0
  return a

def print_mat(mat):
  for row in mat:
    row_ = chop(row)
    print(' '.join(['{:16.9e}']*6).format(*row_))

class Element:
  # name: element name
  # rmat: rmat from beginning to element end
  # rmat1: rmat of just this element
  name = ""
  rmat = np.empty((6,6))
  rmat1 = np.empty((6,6))
  def __init__(self,name,rmat,rmat1):
    self.name = name
    self.rmat = np.array(rmat)
    self.rmat1 = np.array(rmat1)

ix = nhead
record={}
rmat=[None]*6
lat = []
for i in range(nrec):
  rix = i*recsize+nhead
  ele_name = lines[rix].split()[0]
  for j in range(6):
    rmat[j] = [float(x) for x in chunk(lines[rix+j+2])[0:6]]
  if i==0:
    rmat1 = rmat
  else:
    rmat1 = np.matmul(rmat, inv(lat[-1].rmat))
  lat.append(Element(ele_name,rmat,rmat1))

for i in range(nrec):
  print()
  print(lat[i].name)
  #print("full rmat")
  #print_mat(lat[i].rmat)
  print("r55   r56   r65   r66")
  mat = chop(lat[i].rmat1)
  print('{:10.7f}   {:10.7f}   {:10.7f}   {:10.7f}'.format(mat[4,4],mat[4,5],mat[5,4],mat[5,5]))


