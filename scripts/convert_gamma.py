#!/bin/env python3

from numpy import sqrt

me = 0.510998950e6

def from_energy(energy):
  gamma = energy/me
  beta = sqrt(1.0 - 1.0/gamma/gamma)
  print("energy: {}".format(energy))
  print("gamma:  {}".format(gamma))
  print("beta:   {}".format(beta))

from_energy(1.5008701E+07)
from_energy(7.3008701E+07)
